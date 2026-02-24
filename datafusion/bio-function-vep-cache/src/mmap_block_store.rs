//! mmap-backed Arrow window block storage primitives for cache formats v2-v4.

use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{Cursor, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
#[cfg(unix)]
use std::{ffi::c_void, os::fd::AsRawFd, ptr};

use datafusion::arrow::array::{
    ArrayRef, BooleanArray, Float32Array, Float64Array, Int8Array, Int16Array, Int32Array,
    Int64Array, RecordBatch, UInt8Array, UInt16Array, UInt32Array, UInt64Array,
};
use datafusion::arrow::buffer::{
    BooleanBuffer, Buffer, MutableBuffer, NullBuffer, OffsetBuffer, ScalarBuffer,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::common::{DataFusionError, Result};

const BLOCK_FILE_MAGIC: &[u8; 8] = b"VEPBLK2\0";
const FRAME_LEN_BYTES: u64 = 8;
const DEFAULT_BLOCK_FILE_ID: u32 = 0;
pub const WINDOW_BLOCK_REF_ENCODED_LEN: usize = 20;
const COLUMNAR_PAYLOAD_MAGIC: &[u8; 8] = b"V2COLM01";
const COLUMNAR_HEADER_LEN: usize = 8 + 4 + 2; // magic + row_count + col_count
const COLUMNAR_ENTRY_LEN: usize = 2 + 4 + 4; // col_idx + offset + length

#[cfg(unix)]
const PROT_READ: i32 = 0x1;
#[cfg(unix)]
const MAP_PRIVATE: i32 = 0x2;

#[cfg(unix)]
const MADV_SEQUENTIAL: i32 = 2;

#[cfg(unix)]
unsafe extern "C" {
    fn mmap(
        addr: *mut c_void,
        len: usize,
        prot: i32,
        flags: i32,
        fd: i32,
        offset: i64,
    ) -> *mut c_void;
    fn munmap(addr: *mut c_void, len: usize) -> i32;
    fn madvise(addr: *mut c_void, len: usize, advice: i32) -> i32;
}

#[cfg(unix)]
struct MmapRegion {
    ptr: *mut u8,
    len: usize,
}

#[cfg(unix)]
impl MmapRegion {
    fn as_slice(&self) -> &[u8] {
        // SAFETY: mapping is created with `mmap` for `len` bytes and remains
        // valid for the lifetime of this region. We expose read-only bytes.
        unsafe { std::slice::from_raw_parts(self.ptr as *const u8, self.len) }
    }
}

#[cfg(unix)]
impl Drop for MmapRegion {
    fn drop(&mut self) {
        // SAFETY: pointer and length come from successful `mmap`.
        let rc = unsafe { munmap(self.ptr.cast::<c_void>(), self.len) };
        debug_assert_eq!(rc, 0, "munmap failed for mapped block region");
    }
}

#[cfg(unix)]
// SAFETY: mappings are created read-only and exposed only as shared byte slices.
unsafe impl Send for MmapRegion {}

#[cfg(unix)]
// SAFETY: no mutable aliasing is provided; region contents are immutable.
unsafe impl Sync for MmapRegion {}

#[cfg(not(unix))]
struct MmapRegion {
    bytes: Vec<u8>,
}

#[cfg(not(unix))]
impl MmapRegion {
    fn as_slice(&self) -> &[u8] {
        &self.bytes
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowBlockCodec {
    None,
    Lz4,
    Zstd,
}

impl Default for WindowBlockCodec {
    fn default() -> Self {
        Self::None
    }
}

impl WindowBlockCodec {
    pub fn parse(s: &str) -> Result<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "none" => Ok(Self::None),
            "lz4" | "lz4_frame" => Ok(Self::Lz4),
            "zstd" => Ok(Self::Zstd),
            other => Err(DataFusionError::Execution(format!(
                "unsupported v2 block codec '{other}'; supported: none, lz4, zstd"
            ))),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Lz4 => "lz4",
            Self::Zstd => "zstd",
        }
    }

    fn to_ipc_compression(self) -> Option<arrow_ipc::CompressionType> {
        match self {
            Self::None => None,
            Self::Lz4 => Some(arrow_ipc::CompressionType::LZ4_FRAME),
            Self::Zstd => Some(arrow_ipc::CompressionType::ZSTD),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowBlockRef {
    pub file_id: u32,
    pub offset: u64,
    pub length: u32,
    pub row_count: u32,
}

impl WindowBlockRef {
    pub fn to_bytes(self) -> [u8; WINDOW_BLOCK_REF_ENCODED_LEN] {
        let mut out = [0u8; WINDOW_BLOCK_REF_ENCODED_LEN];
        out[0..4].copy_from_slice(&self.file_id.to_le_bytes());
        out[4..12].copy_from_slice(&self.offset.to_le_bytes());
        out[12..16].copy_from_slice(&self.length.to_le_bytes());
        out[16..20].copy_from_slice(&self.row_count.to_le_bytes());
        out
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        if bytes.len() != WINDOW_BLOCK_REF_ENCODED_LEN {
            return Err(DataFusionError::Execution(format!(
                "invalid WindowBlockRef length: expected {}, got {}",
                WINDOW_BLOCK_REF_ENCODED_LEN,
                bytes.len()
            )));
        }
        let file_id = u32::from_le_bytes(bytes[0..4].try_into().unwrap());
        let offset = u64::from_le_bytes(bytes[4..12].try_into().unwrap());
        let length = u32::from_le_bytes(bytes[12..16].try_into().unwrap());
        let row_count = u32::from_le_bytes(bytes[16..20].try_into().unwrap());
        Ok(Self {
            file_id,
            offset,
            length,
            row_count,
        })
    }
}

#[derive(Debug, Clone, Copy)]
struct ColumnSlice {
    column_index: usize,
    offset: usize,
    length: usize,
}

pub struct MmapArrowBlockStore {
    blocks_dir: PathBuf,
    codec: WindowBlockCodec,
    write_lock: Mutex<()>,
    mmap_cache: Mutex<HashMap<u32, Arc<MmapRegion>>>,
}

impl MmapArrowBlockStore {
    pub fn open(blocks_dir: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_codec(blocks_dir, WindowBlockCodec::None)
    }

    pub fn open_with_codec(blocks_dir: impl AsRef<Path>, codec: WindowBlockCodec) -> Result<Self> {
        let blocks_dir = blocks_dir.as_ref().to_path_buf();
        fs::create_dir_all(&blocks_dir).map_err(io_err)?;
        Ok(Self {
            blocks_dir,
            codec,
            write_lock: Mutex::new(()),
            mmap_cache: Mutex::new(HashMap::new()),
        })
    }

    pub fn blocks_dir(&self) -> &Path {
        &self.blocks_dir
    }

    pub fn codec(&self) -> WindowBlockCodec {
        self.codec
    }

    pub fn append_batch(&self, batch: &RecordBatch) -> Result<WindowBlockRef> {
        let payload = serialize_batch_payload(batch, self.codec)?;
        let payload_len_u32 = u32::try_from(payload.len()).map_err(|_| {
            DataFusionError::Execution(format!(
                "window payload too large for u32 length: {} bytes",
                payload.len()
            ))
        })?;
        let row_count = u32::try_from(batch.num_rows()).map_err(|_| {
            DataFusionError::Execution(format!(
                "batch row count too large for u32: {}",
                batch.num_rows()
            ))
        })?;

        let _guard = self.write_lock.lock().unwrap();
        let file_id = DEFAULT_BLOCK_FILE_ID;
        let path = self.block_file_path(file_id);
        let mut file = OpenOptions::new()
            .create(true)
            .read(true)
            .append(true)
            .open(&path)
            .map_err(io_err)?;
        let offset_before = ensure_block_file_header(&mut file)?;

        let frame_len = u64::try_from(payload.len()).map_err(|_| {
            DataFusionError::Execution(format!(
                "window payload length overflow converting usize->u64: {}",
                payload.len()
            ))
        })?;
        file.write_all(&frame_len.to_le_bytes()).map_err(io_err)?;
        let payload_offset = offset_before + FRAME_LEN_BYTES;
        file.write_all(&payload).map_err(io_err)?;
        file.flush().map_err(io_err)?;

        // New append extends file length. Drop stale mmap so next read remaps.
        self.mmap_cache.lock().unwrap().remove(&file_id);

        Ok(WindowBlockRef {
            file_id,
            offset: payload_offset,
            length: payload_len_u32,
            row_count,
        })
    }

    pub fn read_batch(&self, block_ref: WindowBlockRef) -> Result<RecordBatch> {
        let mmap = self.get_or_map_file(block_ref.file_id)?;
        let mmap_bytes = mmap.as_slice();
        let payload = payload_slice(mmap_bytes, block_ref)?;
        let batch = if let Some((row_count, slices)) = parse_columnar_payload(payload)? {
            let mut columns: Vec<Option<ArrayRef>> = vec![None; slices.len()];
            let mut fields: Vec<Option<Field>> = vec![None; slices.len()];
            for slice in slices {
                if slice.column_index >= columns.len() {
                    return Err(DataFusionError::Execution(format!(
                        "columnar payload has non-contiguous column index {} (num_cols={})",
                        slice.column_index,
                        columns.len()
                    )));
                }
                let data = slice_bytes(payload, slice.offset, slice.length)?;
                let (array, field) = deserialize_column_ipc(data)?;
                columns[slice.column_index] = Some(array);
                fields[slice.column_index] = Some(field);
            }

            let arrays = columns
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "columnar payload missing one or more columns".to_string(),
                    )
                })?;
            let schema_fields =
                fields
                    .into_iter()
                    .collect::<Option<Vec<_>>>()
                    .ok_or_else(|| {
                        DataFusionError::Execution(
                            "columnar payload missing one or more column fields".to_string(),
                        )
                    })?;

            if block_ref.row_count != 0 && block_ref.row_count != row_count {
                return Err(DataFusionError::Execution(format!(
                    "window block row_count mismatch: ref={} payload={}",
                    block_ref.row_count, row_count
                )));
            }

            RecordBatch::try_new(Arc::new(Schema::new(schema_fields)), arrays)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
        } else {
            deserialize_batch_ipc(payload)?
        };

        if block_ref.row_count != 0
            && usize::try_from(block_ref.row_count).ok() != Some(batch.num_rows())
        {
            return Err(DataFusionError::Execution(format!(
                "window block row_count mismatch: ref={} actual={}",
                block_ref.row_count,
                batch.num_rows()
            )));
        }
        Ok(batch)
    }

    pub fn read_columns(
        &self,
        block_ref: WindowBlockRef,
        schema: &SchemaRef,
        column_positions: &[usize],
    ) -> Result<Vec<ArrayRef>> {
        let mmap = self.get_or_map_file(block_ref.file_id)?;
        let mmap_bytes = mmap.as_slice();
        let payload = payload_slice(mmap_bytes, block_ref)?;

        if let Some((row_count, slices)) = parse_columnar_payload(payload)? {
            if block_ref.row_count != 0 && block_ref.row_count != row_count {
                return Err(DataFusionError::Execution(format!(
                    "window block row_count mismatch: ref={} payload={}",
                    block_ref.row_count, row_count
                )));
            }

            if is_full_projection(column_positions, schema.fields().len()) {
                let mut ordered = slices;
                ordered.sort_by_key(|s| s.column_index);
                if ordered.iter().enumerate().any(|(i, s)| s.column_index != i) {
                    return Err(DataFusionError::Execution(
                        "columnar payload missing contiguous columns for full projection"
                            .to_string(),
                    ));
                }
                let mut out = Vec::with_capacity(column_positions.len());
                for (col_pos, slice) in ordered.iter().enumerate() {
                    let data = slice_bytes(payload, slice.offset, slice.length)?;
                    let (array, field) = deserialize_column_ipc(data)?;
                    let expected_type = schema.field(col_pos).data_type();
                    if field.data_type() != expected_type {
                        return Err(DataFusionError::Execution(format!(
                            "column type mismatch at index {}: payload={} schema={}",
                            col_pos,
                            field.data_type(),
                            expected_type
                        )));
                    }
                    out.push(array);
                }
                return Ok(out);
            }

            let mut slice_by_col: HashMap<usize, ColumnSlice> =
                HashMap::with_capacity(slices.len());
            for s in slices {
                slice_by_col.insert(s.column_index, s);
            }
            let mut out = Vec::with_capacity(column_positions.len());
            for &col_pos in column_positions {
                if col_pos >= schema.fields().len() {
                    return Err(DataFusionError::Execution(format!(
                        "requested column index {} out of bounds for schema width {}",
                        col_pos,
                        schema.fields().len()
                    )));
                }
                let slice = slice_by_col.get(&col_pos).ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "columnar payload missing requested column index {col_pos}"
                    ))
                })?;
                let data = slice_bytes(payload, slice.offset, slice.length)?;
                let (array, field) = deserialize_column_ipc(data)?;
                let expected_type = schema.field(col_pos).data_type();
                if field.data_type() != expected_type {
                    return Err(DataFusionError::Execution(format!(
                        "column type mismatch at index {}: payload={} schema={}",
                        col_pos,
                        field.data_type(),
                        expected_type
                    )));
                }
                out.push(array);
            }
            return Ok(out);
        }

        let batch = deserialize_batch_ipc(payload)?;
        if block_ref.row_count != 0
            && usize::try_from(block_ref.row_count).ok() != Some(batch.num_rows())
        {
            return Err(DataFusionError::Execution(format!(
                "window block row_count mismatch: ref={} actual={}",
                block_ref.row_count,
                batch.num_rows()
            )));
        }
        let mut out = Vec::with_capacity(column_positions.len());
        for &col_pos in column_positions {
            if col_pos >= batch.num_columns() {
                return Err(DataFusionError::Execution(format!(
                    "requested column index {} out of bounds for batch width {}",
                    col_pos,
                    batch.num_columns()
                )));
            }
            out.push(batch.column(col_pos).clone());
        }
        Ok(out)
    }

    /// Append a whole-batch IPC payload (v3 format — no columnar split).
    ///
    /// Unlike `append_batch` (which writes per-column IPC blobs via the V2COLM01 format),
    /// this writes a single multi-column Arrow IPC stream, eliminating N-1 redundant
    /// StreamReader constructions and compression frames.
    pub fn append_batch_whole(&self, batch: &RecordBatch) -> Result<WindowBlockRef> {
        let payload = serialize_batch_ipc(batch, self.codec)?;
        let payload_len_u32 = u32::try_from(payload.len()).map_err(|_| {
            DataFusionError::Execution(format!(
                "window payload too large for u32 length: {} bytes",
                payload.len()
            ))
        })?;
        let row_count = u32::try_from(batch.num_rows()).map_err(|_| {
            DataFusionError::Execution(format!(
                "batch row count too large for u32: {}",
                batch.num_rows()
            ))
        })?;

        let _guard = self.write_lock.lock().unwrap();
        let file_id = DEFAULT_BLOCK_FILE_ID;
        let path = self.block_file_path(file_id);
        let mut file = OpenOptions::new()
            .create(true)
            .read(true)
            .append(true)
            .open(&path)
            .map_err(io_err)?;
        let offset_before = ensure_block_file_header(&mut file)?;

        let frame_len = u64::try_from(payload.len()).map_err(|_| {
            DataFusionError::Execution(format!(
                "window payload length overflow converting usize->u64: {}",
                payload.len()
            ))
        })?;
        file.write_all(&frame_len.to_le_bytes()).map_err(io_err)?;
        let payload_offset = offset_before + FRAME_LEN_BYTES;
        file.write_all(&payload).map_err(io_err)?;
        file.flush().map_err(io_err)?;

        self.mmap_cache.lock().unwrap().remove(&file_id);

        Ok(WindowBlockRef {
            file_id,
            offset: payload_offset,
            length: payload_len_u32,
            row_count,
        })
    }

    /// Read a batch from a block, only deserializing the requested columns.
    ///
    /// For non-columnar payloads (v3 whole-batch IPC), uses Arrow IPC projection
    /// to skip decompression of unrequested columns. For columnar payloads (v2),
    /// only the requested column blobs are deserialized.
    pub fn read_batch_projected(
        &self,
        block_ref: WindowBlockRef,
        column_positions: &[usize],
    ) -> Result<RecordBatch> {
        let mmap = self.get_or_map_file(block_ref.file_id)?;
        let mmap_bytes = mmap.as_slice();
        let payload = payload_slice(mmap_bytes, block_ref)?;

        if let Some((row_count, slices)) = parse_columnar_payload(payload)? {
            // Columnar payload (v2): selectively deserialize requested columns.
            if block_ref.row_count != 0 && block_ref.row_count != row_count {
                return Err(DataFusionError::Execution(format!(
                    "window block row_count mismatch: ref={} payload={}",
                    block_ref.row_count, row_count
                )));
            }
            let mut slice_by_col: HashMap<usize, ColumnSlice> =
                HashMap::with_capacity(slices.len());
            for s in slices {
                slice_by_col.insert(s.column_index, s);
            }
            let mut arrays = Vec::with_capacity(column_positions.len());
            let mut fields = Vec::with_capacity(column_positions.len());
            for &col_pos in column_positions {
                let slice = slice_by_col.get(&col_pos).ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "columnar payload missing requested column index {col_pos}"
                    ))
                })?;
                let data = slice_bytes(payload, slice.offset, slice.length)?;
                let (array, field) = deserialize_column_ipc(data)?;
                arrays.push(array);
                fields.push(field);
            }
            let schema = Arc::new(Schema::new(fields));
            let batch = RecordBatch::try_new(schema, arrays)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            return Ok(batch);
        }

        // Non-columnar payload (v3): use IPC projection to skip unused columns.
        let batch = deserialize_batch_ipc_projected(payload, column_positions)?;
        if block_ref.row_count != 0
            && usize::try_from(block_ref.row_count).ok() != Some(batch.num_rows())
        {
            return Err(DataFusionError::Execution(format!(
                "window block row_count mismatch: ref={} actual={}",
                block_ref.row_count,
                batch.num_rows()
            )));
        }
        Ok(batch)
    }

    /// Append a V4 (raw Arrow buffer) payload to the block file.
    pub fn append_batch_v4(&self, batch: &RecordBatch) -> Result<WindowBlockRef> {
        let payload = serialize_batch_v4(batch, self.codec)?;
        let payload_len_u32 = u32::try_from(payload.len()).map_err(|_| {
            DataFusionError::Execution(format!(
                "V4 window payload too large: {} bytes",
                payload.len()
            ))
        })?;
        let row_count = u32::try_from(batch.num_rows()).map_err(|_| {
            DataFusionError::Execution(format!(
                "V4 batch row count too large: {}",
                batch.num_rows()
            ))
        })?;

        let _guard = self.write_lock.lock().unwrap();
        let file_id = DEFAULT_BLOCK_FILE_ID;
        let path = self.block_file_path(file_id);
        let mut file = OpenOptions::new()
            .create(true)
            .read(true)
            .append(true)
            .open(&path)
            .map_err(io_err)?;
        let offset_before = ensure_block_file_header(&mut file)?;

        let frame_len = u64::try_from(payload.len()).map_err(|_| {
            DataFusionError::Execution(format!("V4 payload length overflow: {}", payload.len()))
        })?;
        file.write_all(&frame_len.to_le_bytes()).map_err(io_err)?;
        let payload_offset = offset_before + FRAME_LEN_BYTES;
        file.write_all(&payload).map_err(io_err)?;
        file.flush().map_err(io_err)?;

        self.mmap_cache.lock().unwrap().remove(&file_id);

        Ok(WindowBlockRef {
            file_id,
            offset: payload_offset,
            length: payload_len_u32,
            row_count,
        })
    }

    /// Read selected columns from a V4 block, returning Arrow arrays directly.
    pub fn read_columns_v4(
        &self,
        block_ref: WindowBlockRef,
        schema: &SchemaRef,
        column_positions: &[usize],
    ) -> Result<Vec<ArrayRef>> {
        let mmap = self.get_or_map_file(block_ref.file_id)?;
        let mmap_bytes = mmap.as_slice();
        let payload = payload_slice(mmap_bytes, block_ref)?;
        deserialize_batch_v4_projected(payload, schema, column_positions)
    }

    pub fn invalidate_mmap_cache(&self) {
        self.mmap_cache.lock().unwrap().clear();
    }

    fn get_or_map_file(&self, file_id: u32) -> Result<Arc<MmapRegion>> {
        if let Some(m) = self.mmap_cache.lock().unwrap().get(&file_id) {
            return Ok(Arc::clone(m));
        }

        let path = self.block_file_path(file_id);
        let mut file = File::open(&path).map_err(io_err)?;
        verify_block_file_header(&mut file)?;
        let mapped = Arc::new(map_file_readonly(&file)?);

        let mut guard = self.mmap_cache.lock().unwrap();
        let entry = guard.entry(file_id).or_insert_with(|| Arc::clone(&mapped));
        Ok(Arc::clone(entry))
    }

    fn block_file_path(&self, file_id: u32) -> PathBuf {
        self.blocks_dir.join(format!("block_{file_id:06}.arrowblk"))
    }
}

fn io_err(e: std::io::Error) -> DataFusionError {
    DataFusionError::External(Box::new(e))
}

fn arrow_err(e: datafusion::arrow::error::ArrowError) -> DataFusionError {
    DataFusionError::ArrowError(Box::new(e), None)
}

fn ensure_block_file_header(file: &mut File) -> Result<u64> {
    let len = file.metadata().map_err(io_err)?.len();
    if len == 0 {
        file.write_all(BLOCK_FILE_MAGIC).map_err(io_err)?;
        return Ok(BLOCK_FILE_MAGIC.len() as u64);
    }
    if len < BLOCK_FILE_MAGIC.len() as u64 {
        return Err(DataFusionError::Execution(format!(
            "block file too short: expected at least {} bytes, got {len}",
            BLOCK_FILE_MAGIC.len()
        )));
    }
    verify_block_file_header(file)?;
    Ok(len)
}

fn verify_block_file_header(file: &mut File) -> Result<()> {
    file.seek(SeekFrom::Start(0)).map_err(io_err)?;
    let mut hdr = [0u8; 8];
    file.read_exact(&mut hdr).map_err(io_err)?;
    if &hdr != BLOCK_FILE_MAGIC {
        return Err(DataFusionError::Execution(format!(
            "invalid block file header: expected {BLOCK_FILE_MAGIC:?}, got {hdr:?}"
        )));
    }
    Ok(())
}

#[cfg(unix)]
fn map_file_readonly(file: &File) -> Result<MmapRegion> {
    let len_u64 = file.metadata().map_err(io_err)?.len();
    let len = usize::try_from(len_u64).map_err(|_| {
        DataFusionError::Execution(format!("block file length does not fit usize: {len_u64}"))
    })?;
    if len == 0 {
        return Err(DataFusionError::Execution(
            "cannot mmap empty block file".to_string(),
        ));
    }

    let fd = file.as_raw_fd();
    // SAFETY: fd is valid for this file, len > 0, offset is page-aligned 0,
    // mapping is read-only/private and not mutated through this API.
    let ptr = unsafe { mmap(ptr::null_mut(), len, PROT_READ, MAP_PRIVATE, fd, 0) };
    if ptr as isize == -1 {
        return Err(io_err(std::io::Error::last_os_error()));
    }
    // Hint sequential access pattern for read-ahead optimisation.
    unsafe {
        madvise(ptr, len, MADV_SEQUENTIAL);
    }
    Ok(MmapRegion {
        ptr: ptr.cast::<u8>(),
        len,
    })
}

#[cfg(not(unix))]
fn map_file_readonly(file: &File) -> Result<MmapRegion> {
    let mut bytes = Vec::new();
    let mut f = file.try_clone().map_err(io_err)?;
    f.seek(SeekFrom::Start(0)).map_err(io_err)?;
    f.read_to_end(&mut bytes).map_err(io_err)?;
    Ok(MmapRegion { bytes })
}

fn payload_slice(mmap_bytes: &[u8], block_ref: WindowBlockRef) -> Result<&[u8]> {
    let start = usize::try_from(block_ref.offset).map_err(|_| {
        DataFusionError::Execution(format!(
            "window block offset does not fit usize: {}",
            block_ref.offset
        ))
    })?;
    let len = usize::try_from(block_ref.length).map_err(|_| {
        DataFusionError::Execution(format!(
            "window block length does not fit usize: {}",
            block_ref.length
        ))
    })?;
    let end = start.checked_add(len).ok_or_else(|| {
        DataFusionError::Execution(format!(
            "window block slice overflow: start={start} len={len}"
        ))
    })?;
    if end > mmap_bytes.len() {
        return Err(DataFusionError::Execution(format!(
            "window block reference out of bounds: file_id={} offset={} length={} file_len={}",
            block_ref.file_id,
            block_ref.offset,
            block_ref.length,
            mmap_bytes.len()
        )));
    }
    Ok(&mmap_bytes[start..end])
}

fn slice_bytes(data: &[u8], offset: usize, length: usize) -> Result<&[u8]> {
    let end = offset.checked_add(length).ok_or_else(|| {
        DataFusionError::Execution(format!(
            "slice overflow while reading payload: offset={offset} length={length}"
        ))
    })?;
    if end > data.len() {
        return Err(DataFusionError::Execution(format!(
            "slice out of bounds while reading payload: offset={offset} length={length} payload_len={}",
            data.len()
        )));
    }
    Ok(&data[offset..end])
}

fn serialize_batch_payload(batch: &RecordBatch, codec: WindowBlockCodec) -> Result<Vec<u8>> {
    serialize_batch_columnar_payload(batch, codec)
}

fn serialize_batch_columnar_payload(
    batch: &RecordBatch,
    codec: WindowBlockCodec,
) -> Result<Vec<u8>> {
    let num_cols = batch.num_columns();
    let num_cols_u16 = u16::try_from(num_cols).map_err(|_| {
        DataFusionError::Execution(format!(
            "too many columns for columnar payload format: {} (max {})",
            num_cols,
            u16::MAX
        ))
    })?;
    let row_count_u32 = u32::try_from(batch.num_rows()).map_err(|_| {
        DataFusionError::Execution(format!(
            "batch row count too large for columnar payload: {}",
            batch.num_rows()
        ))
    })?;

    let mut blobs: Vec<Vec<u8>> = Vec::with_capacity(num_cols);
    let schema = batch.schema();
    for col_idx in 0..num_cols {
        let field = schema.field(col_idx);
        let blob = serialize_column_ipc(batch.column(col_idx), field, codec)?;
        blobs.push(blob);
    }

    let header_len = COLUMNAR_HEADER_LEN + num_cols * COLUMNAR_ENTRY_LEN;
    let mut total_len = header_len;
    for blob in &blobs {
        total_len = total_len.checked_add(blob.len()).ok_or_else(|| {
            DataFusionError::Execution("columnar payload length overflow".to_string())
        })?;
    }

    let mut out = Vec::with_capacity(total_len);
    out.extend_from_slice(COLUMNAR_PAYLOAD_MAGIC);
    out.extend_from_slice(&row_count_u32.to_le_bytes());
    out.extend_from_slice(&num_cols_u16.to_le_bytes());

    let mut running_offset = header_len;
    for (col_idx, blob) in blobs.iter().enumerate() {
        let col_idx_u16 = u16::try_from(col_idx).map_err(|_| {
            DataFusionError::Execution(format!(
                "column index {col_idx} does not fit u16 in payload"
            ))
        })?;
        let offset_u32 = u32::try_from(running_offset).map_err(|_| {
            DataFusionError::Execution(format!("payload offset {running_offset} does not fit u32"))
        })?;
        let len_u32 = u32::try_from(blob.len()).map_err(|_| {
            DataFusionError::Execution(format!(
                "column payload length {} does not fit u32",
                blob.len()
            ))
        })?;

        out.extend_from_slice(&col_idx_u16.to_le_bytes());
        out.extend_from_slice(&offset_u32.to_le_bytes());
        out.extend_from_slice(&len_u32.to_le_bytes());

        running_offset = running_offset.checked_add(blob.len()).ok_or_else(|| {
            DataFusionError::Execution("columnar payload length overflow".to_string())
        })?;
    }

    for blob in blobs {
        out.extend_from_slice(&blob);
    }

    Ok(out)
}

fn parse_columnar_payload(payload: &[u8]) -> Result<Option<(u32, Vec<ColumnSlice>)>> {
    if payload.len() < COLUMNAR_PAYLOAD_MAGIC.len()
        || &payload[..COLUMNAR_PAYLOAD_MAGIC.len()] != COLUMNAR_PAYLOAD_MAGIC
    {
        return Ok(None);
    }
    if payload.len() < COLUMNAR_HEADER_LEN {
        return Err(DataFusionError::Execution(format!(
            "invalid short columnar payload header: expected at least {} bytes, got {}",
            COLUMNAR_HEADER_LEN,
            payload.len()
        )));
    }

    let row_count = u32::from_le_bytes(payload[8..12].try_into().unwrap());
    let num_cols = u16::from_le_bytes(payload[12..14].try_into().unwrap()) as usize;
    let index_end = COLUMNAR_HEADER_LEN
        .checked_add(num_cols * COLUMNAR_ENTRY_LEN)
        .ok_or_else(|| {
            DataFusionError::Execution("columnar payload index length overflow".to_string())
        })?;
    if index_end > payload.len() {
        return Err(DataFusionError::Execution(format!(
            "columnar payload index out of bounds: index_end={} payload_len={}",
            index_end,
            payload.len()
        )));
    }

    let mut slices = Vec::with_capacity(num_cols);
    let mut pos = COLUMNAR_HEADER_LEN;
    for _ in 0..num_cols {
        let col_idx = u16::from_le_bytes(payload[pos..pos + 2].try_into().unwrap()) as usize;
        let offset = u32::from_le_bytes(payload[pos + 2..pos + 6].try_into().unwrap()) as usize;
        let length = u32::from_le_bytes(payload[pos + 6..pos + 10].try_into().unwrap()) as usize;
        let _ = slice_bytes(payload, offset, length)?;
        slices.push(ColumnSlice {
            column_index: col_idx,
            offset,
            length,
        });
        pos += COLUMNAR_ENTRY_LEN;
    }

    Ok(Some((row_count, slices)))
}

fn is_full_projection(column_positions: &[usize], schema_width: usize) -> bool {
    column_positions.len() == schema_width
        && column_positions
            .iter()
            .enumerate()
            .all(|(i, &pos)| pos == i)
}

fn serialize_batch_ipc(batch: &RecordBatch, codec: WindowBlockCodec) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    {
        let options = arrow_ipc::writer::IpcWriteOptions::default()
            .try_with_compression(codec.to_ipc_compression())
            .map_err(arrow_err)?;
        let mut writer = arrow_ipc::writer::StreamWriter::try_new_with_options(
            &mut out,
            batch.schema_ref(),
            options,
        )
        .map_err(arrow_err)?;
        writer.write(batch).map_err(arrow_err)?;
        writer.finish().map_err(arrow_err)?;
    }
    Ok(out)
}

fn serialize_column_ipc(
    array: &ArrayRef,
    field: &Field,
    codec: WindowBlockCodec,
) -> Result<Vec<u8>> {
    let schema = Arc::new(Schema::new(vec![field.clone()]));
    let batch = RecordBatch::try_new(schema, vec![array.clone()])
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
    serialize_batch_ipc(&batch, codec)
}

fn deserialize_batch_ipc(data: &[u8]) -> Result<RecordBatch> {
    deserialize_batch_ipc_impl(data, None)
}

fn deserialize_batch_ipc_projected(data: &[u8], projection: &[usize]) -> Result<RecordBatch> {
    deserialize_batch_ipc_impl(data, Some(projection.to_vec()))
}

fn deserialize_batch_ipc_impl(data: &[u8], projection: Option<Vec<usize>>) -> Result<RecordBatch> {
    let cursor = Cursor::new(data);
    let reader = arrow_ipc::reader::StreamReader::try_new(cursor, projection).map_err(arrow_err)?;
    let mut batches = Vec::new();
    for batch in reader {
        batches.push(batch.map_err(arrow_err)?);
    }
    match batches.len() {
        0 => Err(DataFusionError::Execution(
            "block payload contained no batches".to_string(),
        )),
        1 => Ok(batches.pop().unwrap()),
        _ => datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches)
            .map_err(arrow_err),
    }
}

fn deserialize_column_ipc(data: &[u8]) -> Result<(ArrayRef, Field)> {
    let batch = deserialize_batch_ipc(data)?;
    if batch.num_columns() != 1 {
        return Err(DataFusionError::Execution(format!(
            "invalid column payload: expected 1 column, got {}",
            batch.num_columns()
        )));
    }
    Ok((batch.column(0).clone(), batch.schema().field(0).clone()))
}

// ---------------------------------------------------------------------------
// FORMAT V4: Zero-copy Arrow buffer cache format
// ---------------------------------------------------------------------------

const V4_MAGIC: &[u8; 8] = b"V4RAW01\0";
const V4_HEADER_LEN: usize = 16; // magic(8) + row_count(4) + num_cols(2) + codec(2)
const V4_TOC_ENTRY_LEN: usize = 32;
const V4_CODEC_NONE: u16 = 0;
const V4_CODEC_LZ4_BLOCK: u16 = 1;
const V4_CODEC_ZSTD_BLOCK: u16 = 2;

/// Max buffers per column in V4 format (null_bitmap + up to 2 data buffers).
const V4_MAX_BUFFERS: usize = 3;

/// Encode an Arrow DataType to a V4 type code.
///
/// `Utf8View` is mapped to `Utf8` (code 2) — the array is cast before buffer extraction.
fn encode_data_type(dt: &DataType) -> Result<u8> {
    match dt {
        DataType::Int32 => Ok(0),
        DataType::Int64 => Ok(1),
        DataType::Utf8 | DataType::Utf8View => Ok(2),
        DataType::Boolean => Ok(3),
        DataType::Float32 => Ok(4),
        DataType::Float64 => Ok(5),
        DataType::UInt32 => Ok(6),
        DataType::UInt64 => Ok(7),
        DataType::LargeUtf8 => Ok(8),
        DataType::Int8 => Ok(9),
        DataType::Int16 => Ok(10),
        DataType::UInt8 => Ok(11),
        DataType::UInt16 => Ok(12),
        other => Err(DataFusionError::Execution(format!(
            "V4 format does not support data type: {other}"
        ))),
    }
}

/// Decode a V4 type code back to an Arrow DataType.
fn decode_data_type(code: u8) -> Result<DataType> {
    match code {
        0 => Ok(DataType::Int32),
        1 => Ok(DataType::Int64),
        2 => Ok(DataType::Utf8),
        3 => Ok(DataType::Boolean),
        4 => Ok(DataType::Float32),
        5 => Ok(DataType::Float64),
        6 => Ok(DataType::UInt32),
        7 => Ok(DataType::UInt64),
        8 => Ok(DataType::LargeUtf8),
        9 => Ok(DataType::Int8),
        10 => Ok(DataType::Int16),
        11 => Ok(DataType::UInt8),
        12 => Ok(DataType::UInt16),
        other => Err(DataFusionError::Execution(format!(
            "V4 unknown data type code: {other}"
        ))),
    }
}

/// Map WindowBlockCodec to V4 codec id.
fn codec_to_v4(codec: WindowBlockCodec) -> u16 {
    match codec {
        WindowBlockCodec::None => V4_CODEC_NONE,
        WindowBlockCodec::Lz4 => V4_CODEC_LZ4_BLOCK,
        WindowBlockCodec::Zstd => V4_CODEC_ZSTD_BLOCK,
    }
}

/// Compress bytes with V4 codec (LZ4 block, zstd block, or none).
fn v4_compress(raw: &[u8], codec: u16) -> Vec<u8> {
    match codec {
        V4_CODEC_LZ4_BLOCK => lz4_flex::block::compress(raw),
        V4_CODEC_ZSTD_BLOCK => zstd::bulk::compress(raw, 3).expect("zstd compress failed"),
        _ => raw.to_vec(),
    }
}

/// Decompress V4 block data into a pre-allocated buffer.
fn v4_decompress_into(compressed: &[u8], output: &mut [u8], codec: u16) -> Result<()> {
    if codec == V4_CODEC_NONE {
        if compressed.len() != output.len() {
            return Err(DataFusionError::Execution(format!(
                "V4 uncompressed buffer size mismatch: expected {}, got {}",
                output.len(),
                compressed.len()
            )));
        }
        output.copy_from_slice(compressed);
        return Ok(());
    }
    if codec == V4_CODEC_ZSTD_BLOCK {
        let decompressed = zstd::bulk::decompress(compressed, output.len()).map_err(|e| {
            DataFusionError::Execution(format!("V4 zstd block decompression failed: {e}"))
        })?;
        if decompressed.len() != output.len() {
            return Err(DataFusionError::Execution(format!(
                "V4 zstd decompressed size mismatch: expected {}, got {}",
                output.len(),
                decompressed.len()
            )));
        }
        output.copy_from_slice(&decompressed);
        return Ok(());
    }
    // LZ4 block decompress
    lz4_flex::block::decompress_into(compressed, output).map_err(|e| {
        DataFusionError::Execution(format!("V4 LZ4 block decompression failed: {e}"))
    })?;
    Ok(())
}

/// Extract raw buffers from an Arrow array column for V4 serialization.
///
/// Returns: (null_bitmap_buffer, data_buffers) as owned Arrow Buffers.
fn extract_column_buffers(
    array: &dyn datafusion::arrow::array::Array,
    type_code: u8,
) -> Result<(Option<Buffer>, Vec<Buffer>)> {
    let null_buf = array.nulls().map(|n| n.inner().inner().clone());

    let data = array.to_data();
    let mut buffers = Vec::new();

    match type_code {
        // Primitive types: single values buffer
        0 | 1 | 4 | 5 | 6 | 7 | 9 | 10 | 11 | 12 => {
            buffers.push(data.buffers()[0].clone());
        }
        // Boolean: bit-packed values buffer
        3 => {
            buffers.push(data.buffers()[0].clone());
        }
        // Utf8: offsets + values
        2 => {
            buffers.push(data.buffers()[0].clone()); // i32 offsets
            buffers.push(data.buffers()[1].clone()); // byte values
        }
        // LargeUtf8: offsets (i64) + values
        8 => {
            buffers.push(data.buffers()[0].clone()); // i64 offsets
            buffers.push(data.buffers()[1].clone()); // byte values
        }
        _ => {
            return Err(DataFusionError::Execution(format!(
                "V4 extract_column_buffers: unsupported type code {type_code}"
            )));
        }
    }

    Ok((null_buf, buffers))
}

/// Serialize a RecordBatch into V4 raw buffer format.
pub(crate) fn serialize_batch_v4(batch: &RecordBatch, codec: WindowBlockCodec) -> Result<Vec<u8>> {
    let num_cols = batch.num_columns();
    let num_cols_u16 = u16::try_from(num_cols)
        .map_err(|_| DataFusionError::Execution(format!("V4: too many columns: {num_cols}")))?;
    let row_count = u32::try_from(batch.num_rows()).map_err(|_| {
        DataFusionError::Execution(format!("V4: row count too large: {}", batch.num_rows()))
    })?;
    let v4_codec = codec_to_v4(codec);

    // Pre-compute all column buffer data
    struct ColData {
        schema_col_idx: u16,
        type_code: u8,
        null_count: u32,
        // (compressed_size, uncompressed_size) per buffer: [null_bitmap, buf0, buf1]
        buf_sizes: Vec<(u32, u32)>,
        compressed_bufs: Vec<Vec<u8>>,
    }

    let schema = batch.schema();
    let mut col_datas = Vec::with_capacity(num_cols);

    for col_idx in 0..num_cols {
        let field = schema.field(col_idx);
        let type_code = encode_data_type(field.data_type())?;
        let array = batch.column(col_idx);

        // Cast Utf8View → Utf8 so we can store using standard Utf8 buffer layout.
        let array = if *field.data_type() == DataType::Utf8View {
            datafusion::arrow::compute::cast(array, &DataType::Utf8)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
        } else {
            array.clone()
        };

        let null_count = u32::try_from(array.null_count())
            .map_err(|_| DataFusionError::Execution("V4: null_count overflow".to_string()))?;

        let (null_bitmap, data_buffers) = extract_column_buffers(array.as_ref(), type_code)?;

        let mut buf_sizes = Vec::new();
        let mut compressed_bufs = Vec::new();

        // Null bitmap buffer (always slot 0)
        if let Some(ref nulls) = null_bitmap {
            let raw = nulls.as_slice();
            let compressed = v4_compress(raw, v4_codec);
            buf_sizes.push((
                u32::try_from(compressed.len()).map_err(|_| {
                    DataFusionError::Execution("V4: compressed size overflow".to_string())
                })?,
                u32::try_from(raw.len()).map_err(|_| {
                    DataFusionError::Execution("V4: uncompressed size overflow".to_string())
                })?,
            ));
            compressed_bufs.push(compressed);
        } else {
            buf_sizes.push((0, 0)); // no null bitmap
        }

        // Data buffers
        for raw_buf in &data_buffers {
            let raw = raw_buf.as_slice();
            let compressed = v4_compress(raw, v4_codec);
            buf_sizes.push((
                u32::try_from(compressed.len()).map_err(|_| {
                    DataFusionError::Execution("V4: compressed size overflow".to_string())
                })?,
                u32::try_from(raw.len()).map_err(|_| {
                    DataFusionError::Execution("V4: uncompressed size overflow".to_string())
                })?,
            ));
            compressed_bufs.push(compressed);
        }

        let num_buffers = 1 + data_buffers.len(); // 1 for null_bitmap slot + data
        debug_assert!(num_buffers <= V4_MAX_BUFFERS);

        col_datas.push(ColData {
            schema_col_idx: u16::try_from(col_idx)
                .map_err(|_| DataFusionError::Execution("V4: column index overflow".to_string()))?,
            type_code,
            null_count,
            buf_sizes,
            compressed_bufs,
        });
    }

    // Calculate total size
    let toc_size = num_cols * V4_TOC_ENTRY_LEN;
    let payload_size: usize = col_datas
        .iter()
        .map(|c| c.compressed_bufs.iter().map(|b| b.len()).sum::<usize>())
        .sum();
    let total_size = V4_HEADER_LEN + toc_size + payload_size;

    let mut out = Vec::with_capacity(total_size);

    // Write header
    out.extend_from_slice(V4_MAGIC);
    out.extend_from_slice(&row_count.to_le_bytes());
    out.extend_from_slice(&num_cols_u16.to_le_bytes());
    out.extend_from_slice(&v4_codec.to_le_bytes());

    // Write TOC entries (32 bytes each)
    for cd in &col_datas {
        let mut entry = [0u8; V4_TOC_ENTRY_LEN];
        entry[0..2].copy_from_slice(&cd.schema_col_idx.to_le_bytes());
        entry[2] = cd.type_code;
        entry[3] = cd.buf_sizes.len() as u8; // num_buffers
        entry[4..8].copy_from_slice(&cd.null_count.to_le_bytes());
        // Per buffer sizes (up to 3 buffers × 8 bytes = 24 bytes, fits in remaining 24)
        let mut pos = 8;
        for &(comp_size, uncomp_size) in &cd.buf_sizes {
            entry[pos..pos + 4].copy_from_slice(&comp_size.to_le_bytes());
            entry[pos + 4..pos + 8].copy_from_slice(&uncomp_size.to_le_bytes());
            pos += 8;
        }
        out.extend_from_slice(&entry);
    }

    // Write compressed buffer payloads
    for cd in &col_datas {
        for buf in &cd.compressed_bufs {
            out.extend_from_slice(buf);
        }
    }

    Ok(out)
}

/// Deserialize selected columns from a V4 payload into Arrow arrays.
///
/// `column_positions` are schema-order indices of the columns to extract.
/// Returns arrays in the same order as `column_positions`.
pub(crate) fn deserialize_batch_v4_projected(
    payload: &[u8],
    schema: &SchemaRef,
    column_positions: &[usize],
) -> Result<Vec<ArrayRef>> {
    if payload.len() < V4_HEADER_LEN {
        return Err(DataFusionError::Execution(format!(
            "V4 payload too short for header: {} < {V4_HEADER_LEN}",
            payload.len()
        )));
    }
    if &payload[..8] != V4_MAGIC {
        return Err(DataFusionError::Execution(
            "V4 payload magic mismatch".to_string(),
        ));
    }

    let row_count = u32::from_le_bytes(payload[8..12].try_into().unwrap()) as usize;
    let num_cols = u16::from_le_bytes(payload[12..14].try_into().unwrap()) as usize;
    let v4_codec = u16::from_le_bytes(payload[14..16].try_into().unwrap());

    let toc_end = V4_HEADER_LEN + num_cols * V4_TOC_ENTRY_LEN;
    if payload.len() < toc_end {
        return Err(DataFusionError::Execution(format!(
            "V4 payload too short for TOC: {} < {toc_end}",
            payload.len()
        )));
    }

    // Parse TOC into column descriptors with buffer offsets
    struct TocEntry {
        schema_col_idx: usize,
        type_code: u8,
        null_count: u32,
        // (compressed_offset_in_payload, compressed_size, uncompressed_size) per buffer
        buffers: Vec<(usize, u32, u32)>,
    }

    let mut entries = Vec::with_capacity(num_cols);
    let mut data_offset = toc_end; // running offset into payload for buffer data

    for i in 0..num_cols {
        let toc_start = V4_HEADER_LEN + i * V4_TOC_ENTRY_LEN;
        let entry_bytes = &payload[toc_start..toc_start + V4_TOC_ENTRY_LEN];

        let schema_col_idx = u16::from_le_bytes(entry_bytes[0..2].try_into().unwrap()) as usize;
        let type_code = entry_bytes[2];
        let num_buffers = entry_bytes[3] as usize;
        let null_count = u32::from_le_bytes(entry_bytes[4..8].try_into().unwrap());

        let mut buffers = Vec::with_capacity(num_buffers);
        let mut pos = 8;
        for _ in 0..num_buffers {
            let comp_size = u32::from_le_bytes(entry_bytes[pos..pos + 4].try_into().unwrap());
            let uncomp_size = u32::from_le_bytes(entry_bytes[pos + 4..pos + 8].try_into().unwrap());
            buffers.push((data_offset, comp_size, uncomp_size));
            data_offset += comp_size as usize;
            pos += 8;
        }

        entries.push(TocEntry {
            schema_col_idx,
            type_code,
            null_count,
            buffers,
        });
    }

    // Build lookup from schema_col_idx -> toc entry index
    let mut col_idx_to_toc: HashMap<usize, usize> = HashMap::with_capacity(num_cols);
    for (toc_idx, entry) in entries.iter().enumerate() {
        col_idx_to_toc.insert(entry.schema_col_idx, toc_idx);
    }

    // Deserialize only requested columns
    let profile = std::env::var_os("VEP_KV_PROFILE").is_some();
    let mut decompress_time = std::time::Duration::ZERO;
    let mut alloc_time = std::time::Duration::ZERO;
    let mut construct_time = std::time::Duration::ZERO;
    let mut total_compressed: u64 = 0;
    let mut total_uncompressed: u64 = 0;
    let mut total_buffers: u64 = 0;

    let mut out = Vec::with_capacity(column_positions.len());
    for &col_pos in column_positions {
        let toc_idx = col_idx_to_toc.get(&col_pos).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "V4: requested column {col_pos} not found in payload"
            ))
        })?;
        let entry = &entries[*toc_idx];

        // Decompress buffers
        let mut decompressed: Vec<Buffer> = Vec::with_capacity(entry.buffers.len());
        for &(offset, comp_size, uncomp_size) in &entry.buffers {
            if comp_size == 0 && uncomp_size == 0 {
                decompressed.push(Buffer::from(MutableBuffer::new(0)));
                continue;
            }
            if profile {
                total_compressed += comp_size as u64;
                total_uncompressed += uncomp_size as u64;
                total_buffers += 1;
            }
            let alloc_start = std::time::Instant::now();
            let compressed = &payload[offset..offset + comp_size as usize];
            let mut buf = MutableBuffer::new(uncomp_size as usize);
            buf.resize(uncomp_size as usize, 0);
            if profile {
                alloc_time += alloc_start.elapsed();
            }
            let dec_start = std::time::Instant::now();
            v4_decompress_into(compressed, buf.as_slice_mut(), v4_codec)?;
            if profile {
                decompress_time += dec_start.elapsed();
            }
            decompressed.push(Buffer::from(buf));
        }

        // Construct typed array from decompressed buffers
        let construct_start = std::time::Instant::now();
        let null_buf = if entry.null_count > 0 && !decompressed[0].is_empty() {
            Some(NullBuffer::new(BooleanBuffer::new(
                decompressed[0].clone(),
                0,
                row_count,
            )))
        } else {
            None
        };

        let array: ArrayRef = match entry.type_code {
            0 => {
                // Int32
                let values = ScalarBuffer::<i32>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(Int32Array::new(values, null_buf))
            }
            1 => {
                // Int64
                let values = ScalarBuffer::<i64>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(Int64Array::new(values, null_buf))
            }
            2 => {
                // Utf8 (StringArray)
                let offsets = ScalarBuffer::<i32>::new(decompressed[1].clone(), 0, row_count + 1);
                let offsets = OffsetBuffer::new(offsets);
                let values = decompressed[2].clone();
                Arc::new(datafusion::arrow::array::StringArray::new(
                    offsets, values, null_buf,
                ))
            }
            3 => {
                // Boolean
                let values = BooleanBuffer::new(decompressed[1].clone(), 0, row_count);
                Arc::new(BooleanArray::new(values, null_buf))
            }
            4 => {
                // Float32
                let values = ScalarBuffer::<f32>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(Float32Array::new(values, null_buf))
            }
            5 => {
                // Float64
                let values = ScalarBuffer::<f64>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(Float64Array::new(values, null_buf))
            }
            6 => {
                // UInt32
                let values = ScalarBuffer::<u32>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(UInt32Array::new(values, null_buf))
            }
            7 => {
                // UInt64
                let values = ScalarBuffer::<u64>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(UInt64Array::new(values, null_buf))
            }
            8 => {
                // LargeUtf8 (LargeStringArray)
                let offsets = ScalarBuffer::<i64>::new(decompressed[1].clone(), 0, row_count + 1);
                let offsets = OffsetBuffer::new(offsets);
                let values = decompressed[2].clone();
                Arc::new(datafusion::arrow::array::LargeStringArray::new(
                    offsets, values, null_buf,
                ))
            }
            9 => {
                // Int8
                let values = ScalarBuffer::<i8>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(Int8Array::new(values, null_buf))
            }
            10 => {
                // Int16
                let values = ScalarBuffer::<i16>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(Int16Array::new(values, null_buf))
            }
            11 => {
                // UInt8
                let values = ScalarBuffer::<u8>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(UInt8Array::new(values, null_buf))
            }
            12 => {
                // UInt16
                let values = ScalarBuffer::<u16>::new(decompressed[1].clone(), 0, row_count);
                Arc::new(UInt16Array::new(values, null_buf))
            }
            other => {
                return Err(DataFusionError::Execution(format!(
                    "V4 deserialize: unsupported type code {other}"
                )));
            }
        };

        // Validate type against schema (Utf8View in schema matches Utf8 in payload).
        let expected_type = schema.field(col_pos).data_type();
        let decoded_type = decode_data_type(entry.type_code)?;
        let types_match = *expected_type == decoded_type
            || (*expected_type == DataType::Utf8View && decoded_type == DataType::Utf8);
        if !types_match {
            return Err(DataFusionError::Execution(format!(
                "V4 column type mismatch at index {col_pos}: payload={decoded_type} schema={expected_type}"
            )));
        }

        if profile {
            construct_time += construct_start.elapsed();
        }
        out.push(array);
    }

    if profile {
        eprintln!(
            "[v4-deser] cols={} bufs={} alloc={:.3}ms decompress={:.3}ms construct={:.3}ms comp={:.1}MB uncomp={:.1}MB ratio={:.2}",
            column_positions.len(),
            total_buffers,
            alloc_time.as_secs_f64() * 1000.0,
            decompress_time.as_secs_f64() * 1000.0,
            construct_time.as_secs_f64() * 1000.0,
            total_compressed as f64 / 1_048_576.0,
            total_uncompressed as f64 / 1_048_576.0,
            if total_compressed > 0 {
                total_uncompressed as f64 / total_compressed as f64
            } else {
                0.0
            },
        );
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    fn test_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_window_block_ref_roundtrip() {
        let r = WindowBlockRef {
            file_id: 7,
            offset: 12345,
            length: 678,
            row_count: 42,
        };
        let bytes = r.to_bytes();
        let got = WindowBlockRef::from_bytes(&bytes).unwrap();
        assert_eq!(r, got);
    }

    #[test]
    fn test_append_and_read_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let store = MmapArrowBlockStore::open(dir.path()).unwrap();
        let batch = test_batch();

        let loc = store.append_batch(&batch).unwrap();
        assert_eq!(loc.row_count, batch.num_rows() as u32);
        let got = store.read_batch(loc).unwrap();
        assert_eq!(got, batch);
    }

    #[test]
    fn test_append_and_read_roundtrip_lz4() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            MmapArrowBlockStore::open_with_codec(dir.path(), WindowBlockCodec::Lz4).unwrap();
        let batch = test_batch();

        let loc = store.append_batch(&batch).unwrap();
        let got = store.read_batch(loc).unwrap();
        assert_eq!(got, batch);
    }

    #[test]
    fn test_append_and_read_roundtrip_zstd() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            MmapArrowBlockStore::open_with_codec(dir.path(), WindowBlockCodec::Zstd).unwrap();
        let batch = test_batch();

        let loc = store.append_batch(&batch).unwrap();
        let got = store.read_batch(loc).unwrap();
        assert_eq!(got, batch);
    }

    #[test]
    fn test_read_columns_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            MmapArrowBlockStore::open_with_codec(dir.path(), WindowBlockCodec::Lz4).unwrap();
        let batch = test_batch();

        let loc = store.append_batch(&batch).unwrap();
        let cols = store
            .read_columns(loc, &batch.schema_ref(), &[0, 3])
            .unwrap();
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].to_data(), batch.column(0).to_data());
        assert_eq!(cols[1].to_data(), batch.column(3).to_data());
    }

    #[test]
    fn test_v4_append_and_read_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let store = MmapArrowBlockStore::open(dir.path()).unwrap();
        let batch = test_batch();

        let loc = store.append_batch_v4(&batch).unwrap();
        assert_eq!(loc.row_count, batch.num_rows() as u32);

        let cols = store
            .read_columns_v4(loc, &batch.schema_ref(), &[0, 1, 2, 3])
            .unwrap();
        assert_eq!(cols.len(), 4);
        for (i, col) in cols.iter().enumerate() {
            assert_eq!(col.to_data(), batch.column(i).to_data());
        }
    }

    #[test]
    fn test_v4_append_and_read_roundtrip_lz4() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            MmapArrowBlockStore::open_with_codec(dir.path(), WindowBlockCodec::Lz4).unwrap();
        let batch = test_batch();

        let loc = store.append_batch_v4(&batch).unwrap();
        let cols = store
            .read_columns_v4(loc, &batch.schema_ref(), &[0, 1, 2, 3])
            .unwrap();
        assert_eq!(cols.len(), 4);
        for (i, col) in cols.iter().enumerate() {
            assert_eq!(col.to_data(), batch.column(i).to_data());
        }
    }

    #[test]
    fn test_v4_projected_read() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            MmapArrowBlockStore::open_with_codec(dir.path(), WindowBlockCodec::Lz4).unwrap();
        let batch = test_batch();

        let loc = store.append_batch_v4(&batch).unwrap();

        // Read only columns 0 and 3 (chrom and allele_string)
        let cols = store
            .read_columns_v4(loc, &batch.schema_ref(), &[0, 3])
            .unwrap();
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].to_data(), batch.column(0).to_data());
        assert_eq!(cols[1].to_data(), batch.column(3).to_data());
    }

    #[test]
    fn test_v4_all_data_types() {
        use datafusion::arrow::array::{
            BooleanArray, Float32Array, Float64Array, LargeStringArray, UInt32Array, UInt64Array,
        };

        let schema = Arc::new(Schema::new(vec![
            Field::new("i32", DataType::Int32, true),
            Field::new("i64", DataType::Int64, false),
            Field::new("utf8", DataType::Utf8, true),
            Field::new("bool", DataType::Boolean, false),
            Field::new("f32", DataType::Float32, true),
            Field::new("f64", DataType::Float64, false),
            Field::new("u32", DataType::UInt32, false),
            Field::new("u64", DataType::UInt64, false),
            Field::new("large_utf8", DataType::LargeUtf8, true),
        ]));

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(Int32Array::from(vec![Some(1), None, Some(3)])),
                Arc::new(Int64Array::from(vec![100, 200, 300])),
                Arc::new(StringArray::from(vec![Some("a"), None, Some("c")])),
                Arc::new(BooleanArray::from(vec![true, false, true])),
                Arc::new(Float32Array::from(vec![Some(1.5), None, Some(3.5)])),
                Arc::new(Float64Array::from(vec![10.0, 20.0, 30.0])),
                Arc::new(UInt32Array::from(vec![1u32, 2, 3])),
                Arc::new(UInt64Array::from(vec![10u64, 20, 30])),
                Arc::new(LargeStringArray::from(vec![Some("xx"), None, Some("zz")])),
            ],
        )
        .unwrap();

        let dir = tempfile::tempdir().unwrap();
        let store =
            MmapArrowBlockStore::open_with_codec(dir.path(), WindowBlockCodec::Lz4).unwrap();

        let loc = store.append_batch_v4(&batch).unwrap();
        let all_cols: Vec<usize> = (0..9).collect();
        let cols = store.read_columns_v4(loc, &schema, &all_cols).unwrap();
        assert_eq!(cols.len(), 9);
        for (i, col) in cols.iter().enumerate() {
            assert_eq!(
                col.to_data(),
                batch.column(i).to_data(),
                "column {i} mismatch"
            );
        }
    }
}
