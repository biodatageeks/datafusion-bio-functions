//! fjall KV store wrapper for reading/writing Arrow IPC batches per window.

use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{ArrayRef, RecordBatch};
use datafusion::arrow::datatypes::{Field, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use fjall::{Config, Keyspace, PartitionCreateOptions, PartitionHandle, PersistMode};

use crate::key_encoding::{EntryType, decode_entry_key, encode_entry_key, encode_window_key};
use crate::mmap_block_store::{MmapArrowBlockStore, WindowBlockCodec, WindowBlockRef};
use crate::position_index::PositionIndex;

const META_PARTITION: &str = "meta";
const DATA_PARTITION: &str = "data";
const SCHEMA_KEY: &[u8] = b"schema";
const WINDOW_SIZE_KEY: &[u8] = b"window_size";
const FORMAT_VERSION_KEY: &[u8] = b"format_version";
const V2_BLOCK_DIR_KEY: &[u8] = b"v2_block_dir";
const V2_BLOCK_CODEC_KEY: &[u8] = b"v2_block_codec";
const V2_WINDOW_REF_PREFIX: &[u8] = b"v2_window_ref:";
const V2_DEFAULT_BLOCK_DIR: &str = "blocks_v2";
const V5_ZSTD_DICT_KEY: &[u8] = b"v5_zstd_dict";

/// Format version 1: columnar (position index + per-column Arrow IPC entries).
pub const FORMAT_V1: u8 = 1;
/// Format version 2: position index in Fjall + mmap-backed Arrow window blocks (per-column IPC).
pub const FORMAT_V2: u8 = 2;
/// Format version 3: position index in Fjall + mmap-backed whole-batch Arrow IPC blocks.
///
/// V3 stores one multi-column Arrow IPC blob per window (instead of N single-column blobs),
/// reducing IPC StreamReader constructions, schema parses, and compression frames from N to 1.
pub const FORMAT_V3: u8 = 3;
/// Format version 4: zero-copy Arrow buffer cache — LZ4-block compressed raw Arrow buffers.
///
/// Eliminates IPC overhead (flatbuffer parsing, StreamReader, double-buffering) by storing
/// raw Arrow buffers that decompress directly into aligned memory.
pub const FORMAT_V4: u8 = 4;
/// Format version 5: position-keyed entries — one fjall entry per genomic position.
///
/// Eliminates window decompression entirely. Each `(chrom, start, end)` maps to
/// a serialized entry containing the allele table and column-major data. Only
/// positions matching VCF queries are fetched.
pub const FORMAT_V5: u8 = 5;

/// Wrapper around fjall `Keyspace` for VEP cache storage.
///
/// Layout:
/// - `meta` partition: schema (Arrow IPC), window size, coordinate metadata
/// - `data` partition: `(chrom, window_id)` -> Arrow IPC serialized `RecordBatch`
pub struct VepKvStore {
    root_path: PathBuf,
    keyspace: Keyspace,
    data: PartitionHandle,
    meta: PartitionHandle,
    schema: SchemaRef,
    window_size: u64,
    format_version: u8,
    v2_blocks: Option<MmapArrowBlockStore>,
    v2_block_codec: Option<WindowBlockCodec>,
    /// Raw zstd dictionary bytes for V5 position entry decompression.
    v5_zstd_dict: Option<Arc<Vec<u8>>>,
}

fn arrow_err(e: datafusion::arrow::error::ArrowError) -> DataFusionError {
    DataFusionError::ArrowError(Box::new(e), None)
}

/// Serialize a RecordBatch to Arrow IPC stream format with LZ4 compression.
///
/// Exposed as `pub(crate)` so that the loader can pre-serialize entries before
/// batch-inserting them into fjall.
pub(crate) fn serialize_ipc_pub(batch: &RecordBatch) -> Result<Vec<u8>> {
    serialize_ipc(batch)
}

pub(crate) fn max_v1_columns() -> usize {
    EntryType::MAX_COLUMN_INDEX as usize + 1
}

pub(crate) fn validate_v1_schema_width(num_columns: usize) -> Result<()> {
    let max_cols = max_v1_columns();
    if num_columns > max_cols {
        return Err(DataFusionError::Execution(format!(
            "v1 cache supports at most {max_cols} columns (max encoded index {}). got {num_columns}",
            EntryType::MAX_COLUMN_INDEX
        )));
    }
    Ok(())
}

pub(crate) fn to_v1_column_index(column_index: usize) -> Result<u8> {
    let max_idx = EntryType::MAX_COLUMN_INDEX as usize;
    if column_index > max_idx {
        return Err(DataFusionError::Execution(format!(
            "v1 cache column index {column_index} exceeds max encoded index {max_idx}"
        )));
    }
    Ok(column_index as u8)
}

pub(crate) fn validate_v1_column_index(column_index: u8) -> Result<()> {
    if column_index > EntryType::MAX_COLUMN_INDEX {
        return Err(DataFusionError::Execution(format!(
            "v1 cache column index {column_index} exceeds max encoded index {}",
            EntryType::MAX_COLUMN_INDEX
        )));
    }
    Ok(())
}

impl VepKvStore {
    /// Open an existing KV store with default settings (256 MB block cache).
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_cache_size(path, 256 * 1024 * 1024)
    }

    /// Open an existing KV store with a custom fjall block cache size (bytes).
    pub fn open_with_cache_size(path: impl AsRef<Path>, cache_size_bytes: u64) -> Result<Self> {
        let root_path = path.as_ref().to_path_buf();
        let keyspace = Config::new(&root_path)
            .cache_size(cache_size_bytes)
            .open()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let data = keyspace
            .open_partition(DATA_PARTITION, PartitionCreateOptions::default())
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        let meta = keyspace
            .open_partition(META_PARTITION, PartitionCreateOptions::default())
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let schema = Self::read_schema(&meta)?;
        let window_size = Self::read_window_size(&meta)?;
        let format_version = Self::read_format_version(&meta)?;
        let (v2_blocks, v2_block_codec) = if (FORMAT_V2..=FORMAT_V4).contains(&format_version) {
            let block_dir = Self::read_v2_block_dir(&meta)?;
            let block_codec = Self::read_v2_block_codec(&meta)?;
            let blocks_path = root_path.join(block_dir);
            (
                Some(MmapArrowBlockStore::open_with_codec(
                    blocks_path,
                    block_codec,
                )?),
                Some(block_codec),
            )
        } else {
            (None, None)
        };

        // Load V5 zstd dictionary if present.
        let v5_zstd_dict = if format_version == FORMAT_V5 {
            meta.get(V5_ZSTD_DICT_KEY)
                .map_err(|e| DataFusionError::External(Box::new(e)))?
                .map(|raw| Arc::new(raw.to_vec()))
        } else {
            None
        };

        Ok(Self {
            root_path,
            keyspace,
            data,
            meta,
            schema,
            window_size,
            format_version,
            v2_blocks,
            v2_block_codec,
            v5_zstd_dict,
        })
    }

    /// Create a new KV store with the given schema and window size.
    pub fn create(path: impl AsRef<Path>, schema: SchemaRef, window_size: u64) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            window_size,
            256 * 1024 * 1024,
            PartitionCreateOptions::default(),
            FORMAT_V1,
            WindowBlockCodec::None,
        )
    }

    /// Create a new KV store with full control over fjall tuning.
    pub fn create_with_options(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        cache_size_bytes: u64,
        partition_opts: PartitionCreateOptions,
    ) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            window_size,
            cache_size_bytes,
            partition_opts,
            FORMAT_V1,
            WindowBlockCodec::None,
        )
    }

    /// Create a new v2 KV store (mmap-backed window payloads + Fjall metadata/index).
    pub fn create_v2(path: impl AsRef<Path>, schema: SchemaRef, window_size: u64) -> Result<Self> {
        Self::create_v2_with_codec(path, schema, window_size, WindowBlockCodec::None)
    }

    /// Create a new v2 KV store with explicit IPC compression codec for window blocks.
    pub fn create_v2_with_codec(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        block_codec: WindowBlockCodec,
    ) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            window_size,
            256 * 1024 * 1024,
            PartitionCreateOptions::default(),
            FORMAT_V2,
            block_codec,
        )
    }

    /// Create a new v2 KV store with full control over fjall tuning.
    pub fn create_v2_with_options(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        cache_size_bytes: u64,
        partition_opts: PartitionCreateOptions,
    ) -> Result<Self> {
        Self::create_v2_with_options_and_codec(
            path,
            schema,
            window_size,
            cache_size_bytes,
            partition_opts,
            WindowBlockCodec::None,
        )
    }

    /// Create a new v3 KV store (whole-batch Arrow IPC per window).
    pub fn create_v3(path: impl AsRef<Path>, schema: SchemaRef, window_size: u64) -> Result<Self> {
        Self::create_v3_with_codec(path, schema, window_size, WindowBlockCodec::None)
    }

    /// Create a new v3 KV store with explicit IPC compression codec.
    pub fn create_v3_with_codec(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        block_codec: WindowBlockCodec,
    ) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            window_size,
            256 * 1024 * 1024,
            PartitionCreateOptions::default(),
            FORMAT_V3,
            block_codec,
        )
    }

    /// Create a new v4 KV store (zero-copy raw Arrow buffers per window).
    pub fn create_v4(path: impl AsRef<Path>, schema: SchemaRef, window_size: u64) -> Result<Self> {
        Self::create_v4_with_codec(path, schema, window_size, WindowBlockCodec::None)
    }

    /// Create a new v4 KV store with explicit LZ4-block compression codec.
    pub fn create_v4_with_codec(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        block_codec: WindowBlockCodec,
    ) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            window_size,
            256 * 1024 * 1024,
            PartitionCreateOptions::default(),
            FORMAT_V4,
            block_codec,
        )
    }

    /// Create a new v5 KV store (position-keyed entries, no mmap block store).
    pub fn create_v5(path: impl AsRef<Path>, schema: SchemaRef) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            0, // window_size not used in V5
            256 * 1024 * 1024,
            PartitionCreateOptions::default(),
            FORMAT_V5,
            WindowBlockCodec::None,
        )
    }

    /// Create a new v2 KV store with full control over fjall tuning and block codec.
    pub fn create_v2_with_options_and_codec(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        cache_size_bytes: u64,
        partition_opts: PartitionCreateOptions,
        block_codec: WindowBlockCodec,
    ) -> Result<Self> {
        Self::create_internal(
            path,
            schema,
            window_size,
            cache_size_bytes,
            partition_opts,
            FORMAT_V2,
            block_codec,
        )
    }

    fn create_internal(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        cache_size_bytes: u64,
        partition_opts: PartitionCreateOptions,
        format_version: u8,
        v2_block_codec: WindowBlockCodec,
    ) -> Result<Self> {
        if format_version == FORMAT_V1 {
            validate_v1_schema_width(schema.fields().len())?;
        } else if format_version != FORMAT_V2
            && format_version != FORMAT_V3
            && format_version != FORMAT_V4
            && format_version != FORMAT_V5
        {
            return Err(DataFusionError::Execution(format!(
                "unsupported cache format version for create: {format_version}"
            )));
        }

        let root_path = path.as_ref().to_path_buf();
        let keyspace = Config::new(&root_path)
            .cache_size(cache_size_bytes)
            .open()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let data = keyspace
            .open_partition(DATA_PARTITION, partition_opts)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        let meta = keyspace
            .open_partition(META_PARTITION, PartitionCreateOptions::default())
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        // Store schema as Arrow IPC (empty batch carries schema).
        let schema_bytes = schema_to_ipc_bytes(&schema)?;
        meta.insert(SCHEMA_KEY, &schema_bytes)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        meta.insert(WINDOW_SIZE_KEY, window_size.to_be_bytes())
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        meta.insert(FORMAT_VERSION_KEY, [format_version])
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let uses_block_store = format_version == FORMAT_V2
            || format_version == FORMAT_V3
            || format_version == FORMAT_V4;
        let v2_blocks = if uses_block_store {
            meta.insert(V2_BLOCK_DIR_KEY, V2_DEFAULT_BLOCK_DIR.as_bytes())
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            meta.insert(V2_BLOCK_CODEC_KEY, v2_block_codec.as_str().as_bytes())
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            Some(MmapArrowBlockStore::open_with_codec(
                root_path.join(V2_DEFAULT_BLOCK_DIR),
                v2_block_codec,
            )?)
        } else {
            None
        };

        keyspace
            .persist(PersistMode::SyncAll)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        Ok(Self {
            root_path,
            keyspace,
            data,
            meta,
            schema,
            window_size,
            format_version,
            v2_blocks,
            v2_block_codec: if uses_block_store {
                Some(v2_block_codec)
            } else {
                None
            },
            v5_zstd_dict: None,
        })
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn window_size(&self) -> u64 {
        self.window_size
    }

    pub fn root_path(&self) -> &Path {
        &self.root_path
    }

    /// Get the format version of this store.
    pub fn format_version(&self) -> u8 {
        self.format_version
    }

    /// Return the v2 block store handle if this cache is opened in format v2.
    pub fn v2_block_store(&self) -> Option<&MmapArrowBlockStore> {
        self.v2_blocks.as_ref()
    }

    /// Return the v2 payload codec if this cache is opened in format v2.
    pub fn v2_block_codec(&self) -> Option<WindowBlockCodec> {
        self.v2_block_codec
    }

    /// Store a v2 window payload reference in metadata.
    pub fn put_window_ref_v2(
        &self,
        chrom: &str,
        window_id: u64,
        block_ref: WindowBlockRef,
    ) -> Result<()> {
        self.ensure_v2()?;
        let key = v2_window_ref_key(chrom, window_id);
        self.meta
            .insert(&key, block_ref.to_bytes())
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Read a v2 window payload reference from metadata.
    pub fn get_window_ref_v2(&self, chrom: &str, window_id: u64) -> Result<Option<WindowBlockRef>> {
        self.ensure_v2()?;
        let key = v2_window_ref_key(chrom, window_id);
        match self.meta.get(&key) {
            Ok(Some(v)) => Ok(Some(WindowBlockRef::from_bytes(&v)?)),
            Ok(None) => Ok(None),
            Err(e) => Err(DataFusionError::External(Box::new(e))),
        }
    }

    /// Append and register a v2 window payload batch.
    pub fn put_window_block_v2(
        &self,
        chrom: &str,
        window_id: u64,
        batch: &RecordBatch,
    ) -> Result<WindowBlockRef> {
        self.ensure_v2()?;
        let blocks = self.v2_blocks.as_ref().unwrap();
        let block_ref = blocks.append_batch(batch)?;
        self.put_window_ref_v2(chrom, window_id, block_ref)?;
        Ok(block_ref)
    }

    /// Read a v2 window payload batch by `(chrom, window_id)`.
    pub fn get_window_block_v2(&self, chrom: &str, window_id: u64) -> Result<Option<RecordBatch>> {
        self.ensure_v2()?;
        let Some(block_ref) = self.get_window_ref_v2(chrom, window_id)? else {
            return Ok(None);
        };
        let blocks = self.v2_blocks.as_ref().unwrap();
        let batch = blocks.read_batch(block_ref)?;
        Ok(Some(batch))
    }

    /// Read selected columns from a v2 window payload by `(chrom, window_id)`.
    ///
    /// `column_positions` are indices in cache schema order.
    pub fn get_window_columns_v2(
        &self,
        chrom: &str,
        window_id: u64,
        column_positions: &[usize],
    ) -> Result<Option<Vec<ArrayRef>>> {
        self.ensure_v2()?;
        let Some(block_ref) = self.get_window_ref_v2(chrom, window_id)? else {
            return Ok(None);
        };
        let blocks = self.v2_blocks.as_ref().unwrap();
        let columns = blocks.read_columns(block_ref, &self.schema, column_positions)?;
        Ok(Some(columns))
    }

    // --- v1 columnar API ---

    /// Write a position index for a window (v1 format).
    pub fn put_position_index(
        &self,
        chrom: &str,
        window_id: u64,
        index: &PositionIndex,
    ) -> Result<()> {
        let key = encode_entry_key(chrom, window_id, EntryType::PositionIndex);
        let value = index.to_bytes();
        self.data
            .insert(&key, &value)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Read a position index for a window (v1 format).
    pub fn get_position_index(&self, chrom: &str, window_id: u64) -> Result<Option<PositionIndex>> {
        let key = encode_entry_key(chrom, window_id, EntryType::PositionIndex);
        match self.data.get(&key) {
            Ok(Some(value)) => {
                let index = PositionIndex::from_bytes(&value)?;
                Ok(Some(index))
            }
            Ok(None) => Ok(None),
            Err(e) => Err(DataFusionError::External(Box::new(e))),
        }
    }

    /// Write a single column's data as Arrow IPC (v1 format).
    pub fn put_column(
        &self,
        chrom: &str,
        window_id: u64,
        col_idx: u8,
        array: &ArrayRef,
        field: &Field,
    ) -> Result<()> {
        self.ensure_v1()?;
        validate_v1_column_index(col_idx)?;
        let key = encode_entry_key(chrom, window_id, EntryType::Column(col_idx));
        let batch = RecordBatch::try_new(
            Arc::new(datafusion::arrow::datatypes::Schema::new(vec![
                field.clone(),
            ])),
            vec![array.clone()],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        let value = serialize_ipc(&batch)?;
        self.data
            .insert(&key, &value)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Read specific columns' data as Arrow arrays (v1 format).
    ///
    /// Returns arrays in the same order as `col_indices`. If the window
    /// has no data, returns `None`.
    pub fn get_columns(
        &self,
        chrom: &str,
        window_id: u64,
        col_indices: &[u8],
    ) -> Result<Option<Vec<(ArrayRef, Field)>>> {
        self.ensure_v1()?;
        let mut results = Vec::with_capacity(col_indices.len());
        for &col_idx in col_indices {
            validate_v1_column_index(col_idx)?;
            let key = encode_entry_key(chrom, window_id, EntryType::Column(col_idx));
            match self.data.get(&key) {
                Ok(Some(value)) => {
                    let batch = deserialize_ipc(&value)?;
                    let field = batch.schema().field(0).clone();
                    results.push((batch.column(0).clone(), field));
                }
                Ok(None) => return Ok(None),
                Err(e) => return Err(DataFusionError::External(Box::new(e))),
            }
        }
        Ok(Some(results))
    }

    /// Persist all data to disk.
    pub fn persist(&self) -> Result<()> {
        self.keyspace
            .persist(PersistMode::SyncAll)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Atomically insert a batch of pre-serialized (key, value) pairs into the data partition.
    ///
    /// Uses fjall's `Batch` API to reduce WAL overhead and L0 segment pressure
    /// compared to individual `insert()` calls.
    pub(crate) fn batch_insert_raw(&self, entries: &[(Vec<u8>, Vec<u8>)]) -> Result<()> {
        if entries.is_empty() {
            return Ok(());
        }
        let mut batch = fjall::Batch::with_capacity(self.keyspace.clone(), entries.len());
        for (key, value) in entries {
            batch.insert(&self.data, key, value);
        }
        batch
            .commit()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Store arbitrary metadata.
    pub fn put_metadata(&self, key: &[u8], value: &[u8]) -> Result<()> {
        self.meta
            .insert(key, value)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Read arbitrary metadata.
    pub fn get_metadata(&self, key: &[u8]) -> Result<Option<Vec<u8>>> {
        match self.meta.get(key) {
            Ok(Some(v)) => Ok(Some(v.to_vec())),
            Ok(None) => Ok(None),
            Err(e) => Err(DataFusionError::External(Box::new(e))),
        }
    }

    /// Iterate all windows for a given chromosome (prefix scan).
    ///
    /// Returns full batches assembled from v1 per-column entries.
    pub fn windows_for_chrom(&self, chrom: &str) -> Result<Vec<(String, u64, RecordBatch)>> {
        self.ensure_v1()?;
        let prefix = crate::key_encoding::chrom_prefix(chrom);
        let mut windows = Vec::new();
        for item in self.data.prefix(&prefix) {
            let (key, _) = item.map_err(|e| DataFusionError::External(Box::new(e)))?;
            if key.len() >= 11 {
                let (c, wid, entry) = decode_entry_key(&key);
                if entry == EntryType::PositionIndex {
                    windows.push((c, wid));
                }
            }
        }
        windows.sort();
        windows.dedup();

        let all_cols: Vec<u8> = (0..self.schema.fields().len())
            .map(to_v1_column_index)
            .collect::<Result<Vec<_>>>()?;
        let mut results = Vec::with_capacity(windows.len());
        for (c, wid) in windows {
            if let Some(cols) = self.get_columns(&c, wid, &all_cols)? {
                let arrays: Vec<ArrayRef> = cols.into_iter().map(|(arr, _)| arr).collect();
                let batch = RecordBatch::try_new(self.schema.clone(), arrays)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                results.push((c, wid, batch));
            }
        }
        Ok(results)
    }

    fn read_schema(meta: &PartitionHandle) -> Result<SchemaRef> {
        let value = meta
            .get(SCHEMA_KEY)
            .map_err(|e| DataFusionError::External(Box::new(e)))?
            .ok_or_else(|| {
                DataFusionError::Execution("KV store missing schema metadata".to_string())
            })?;
        schema_from_ipc_bytes(&value)
    }

    fn read_window_size(meta: &PartitionHandle) -> Result<u64> {
        let value = meta
            .get(WINDOW_SIZE_KEY)
            .map_err(|e| DataFusionError::External(Box::new(e)))?
            .ok_or_else(|| {
                DataFusionError::Execution("KV store missing window_size metadata".to_string())
            })?;
        let bytes: [u8; 8] = value[..8]
            .try_into()
            .map_err(|_| DataFusionError::Execution("Invalid window_size bytes".to_string()))?;
        Ok(u64::from_be_bytes(bytes))
    }

    fn read_v2_block_dir(meta: &PartitionHandle) -> Result<String> {
        let value = meta
            .get(V2_BLOCK_DIR_KEY)
            .map_err(|e| DataFusionError::External(Box::new(e)))?
            .ok_or_else(|| {
                DataFusionError::Execution(
                    "KV store missing v2_block_dir metadata for format v2 cache".to_string(),
                )
            })?;
        let raw = std::str::from_utf8(&value).map_err(|e| {
            DataFusionError::Execution(format!("invalid UTF-8 in v2_block_dir metadata: {e}"))
        })?;
        if raw.is_empty() {
            return Err(DataFusionError::Execution(
                "invalid empty v2_block_dir metadata".to_string(),
            ));
        }
        Ok(raw.to_string())
    }

    fn read_v2_block_codec(meta: &PartitionHandle) -> Result<WindowBlockCodec> {
        let raw = match meta
            .get(V2_BLOCK_CODEC_KEY)
            .map_err(|e| DataFusionError::External(Box::new(e)))?
        {
            Some(v) => v,
            None => {
                // Backward compatibility for older v2 caches.
                return Ok(WindowBlockCodec::None);
            }
        };
        let s = std::str::from_utf8(&raw).map_err(|e| {
            DataFusionError::Execution(format!("invalid UTF-8 in v2_block_codec metadata: {e}"))
        })?;
        WindowBlockCodec::parse(s)
    }

    /// Append and register a v3 window payload batch (whole-batch IPC).
    pub fn put_window_block_v3(
        &self,
        chrom: &str,
        window_id: u64,
        batch: &RecordBatch,
    ) -> Result<WindowBlockRef> {
        self.ensure_v3()?;
        let blocks = self.v2_blocks.as_ref().unwrap();
        let block_ref = blocks.append_batch_whole(batch)?;
        self.put_window_ref_v2(chrom, window_id, block_ref)?;
        Ok(block_ref)
    }

    /// Read a v3 window payload batch by `(chrom, window_id)`.
    pub fn get_window_batch_v3(&self, chrom: &str, window_id: u64) -> Result<Option<RecordBatch>> {
        self.ensure_v3()?;
        let Some(block_ref) = self.get_window_ref_v2(chrom, window_id)? else {
            return Ok(None);
        };
        let blocks = self.v2_blocks.as_ref().unwrap();
        let batch = blocks.read_batch(block_ref)?;
        Ok(Some(batch))
    }

    /// Read selected columns from a v3 window payload by `(chrom, window_id)`.
    ///
    /// Uses Arrow IPC projection to only deserialize requested columns,
    /// skipping decompression of unrequested column buffers.
    pub fn get_window_columns_v3(
        &self,
        chrom: &str,
        window_id: u64,
        column_positions: &[usize],
    ) -> Result<Option<Vec<ArrayRef>>> {
        self.ensure_v3()?;
        let Some(block_ref) = self.get_window_ref_v2(chrom, window_id)? else {
            return Ok(None);
        };
        let blocks = self.v2_blocks.as_ref().unwrap();
        let batch = blocks.read_batch_projected(block_ref, column_positions)?;
        let out: Vec<ArrayRef> = (0..batch.num_columns())
            .map(|i| batch.column(i).clone())
            .collect();
        Ok(Some(out))
    }

    /// Append and register a v4 window payload batch (raw Arrow buffers).
    pub fn put_window_block_v4(
        &self,
        chrom: &str,
        window_id: u64,
        batch: &RecordBatch,
    ) -> Result<WindowBlockRef> {
        self.ensure_v4()?;
        let blocks = self.v2_blocks.as_ref().unwrap();
        let block_ref = blocks.append_batch_v4(batch)?;
        self.put_window_ref_v2(chrom, window_id, block_ref)?;
        Ok(block_ref)
    }

    /// Look up the block reference for a v4 window, returning timing info.
    pub fn get_window_ref_v4_timed(
        &self,
        chrom: &str,
        window_id: u64,
    ) -> Result<(Option<WindowBlockRef>, std::time::Duration)> {
        self.ensure_v4()?;
        let start = std::time::Instant::now();
        let block_ref = self.get_window_ref_v2(chrom, window_id)?;
        let dur = start.elapsed();
        Ok((block_ref, dur))
    }

    /// Read selected columns from a v4 block reference directly (no fjall lookup).
    pub fn read_columns_v4_direct(
        &self,
        block_ref: WindowBlockRef,
        column_positions: &[usize],
    ) -> Result<Vec<ArrayRef>> {
        let blocks = self.v2_blocks.as_ref().unwrap();
        blocks.read_columns_v4(block_ref, &self.schema, column_positions)
    }

    /// Read selected columns from a v4 window payload by `(chrom, window_id)`.
    ///
    /// Returns raw Arrow arrays without IPC deserialization overhead.
    pub fn get_window_columns_v4(
        &self,
        chrom: &str,
        window_id: u64,
        column_positions: &[usize],
    ) -> Result<Option<Vec<ArrayRef>>> {
        self.ensure_v4()?;
        let Some(block_ref) = self.get_window_ref_v2(chrom, window_id)? else {
            return Ok(None);
        };
        let blocks = self.v2_blocks.as_ref().unwrap();
        let columns = blocks.read_columns_v4(block_ref, &self.schema, column_positions)?;
        Ok(Some(columns))
    }

    /// Read a v4 window as a full RecordBatch (all columns).
    pub fn get_window_batch_v4(&self, chrom: &str, window_id: u64) -> Result<Option<RecordBatch>> {
        self.ensure_v4()?;
        let all_cols: Vec<usize> = (0..self.schema.fields().len()).collect();
        let Some(columns) = self.get_window_columns_v4(chrom, window_id, &all_cols)? else {
            return Ok(None);
        };
        let batch = RecordBatch::try_new(self.schema.clone(), columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        Ok(Some(batch))
    }

    // --- v5 position-keyed API ---

    /// Write a position entry (v5 format).
    pub fn put_position_entry(
        &self,
        chrom: &str,
        start: i64,
        end: i64,
        value: &[u8],
    ) -> Result<()> {
        self.ensure_v5()?;
        let key = crate::key_encoding::encode_position_key(chrom, start, end);
        self.data
            .insert(&key, value)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Read a position entry by chrom code + coordinates (v5 format).
    ///
    /// Uses a pre-encoded key buffer for hot-path efficiency.
    pub fn get_position_entry(
        &self,
        chrom_code: u16,
        start: i64,
        end: i64,
    ) -> Result<Option<fjall::UserValue>> {
        self.ensure_v5()?;
        let mut key_buf = Vec::with_capacity(18);
        crate::key_encoding::encode_position_key_buf(chrom_code, start, end, &mut key_buf);
        match self.data.get(&key_buf) {
            Ok(v) => Ok(v),
            Err(e) => Err(DataFusionError::External(Box::new(e))),
        }
    }

    /// Store a trained zstd dictionary for V5 position entries.
    pub fn store_v5_dict(&self, dict_bytes: &[u8]) -> Result<()> {
        self.ensure_v5()?;
        self.meta
            .insert(V5_ZSTD_DICT_KEY, dict_bytes)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Whether this V5 store has a zstd dictionary loaded.
    pub fn has_v5_dict(&self) -> bool {
        self.v5_zstd_dict.is_some()
    }

    /// Read and decompress a V5 position entry using the zstd dictionary.
    ///
    /// If no dictionary is loaded (legacy uncompressed V5 store), returns the
    /// raw bytes without decompression.
    pub fn get_position_entry_decompressed(
        &self,
        chrom_code: u16,
        start: i64,
        end: i64,
    ) -> Result<Option<Vec<u8>>> {
        let raw = self.get_position_entry(chrom_code, start, end)?;
        match raw {
            None => Ok(None),
            Some(compressed) => {
                if let Some(dict) = &self.v5_zstd_dict {
                    let mut decompressor = zstd::bulk::Decompressor::with_dictionary(dict)
                        .map_err(|e| {
                            DataFusionError::Execution(format!(
                                "failed to create zstd decompressor: {e}"
                            ))
                        })?;
                    // Dictionary compression on small entries can achieve high ratios;
                    // use a generous capacity (at least 4KB, or 16x compressed size).
                    let capacity = (compressed.len() * 16).max(4096);
                    let decompressed =
                        decompressor
                            .decompress(&compressed, capacity)
                            .map_err(|e| {
                                DataFusionError::Execution(format!(
                                    "zstd decompression failed (capacity={capacity}): {e}"
                                ))
                            })?;
                    Ok(Some(decompressed))
                } else {
                    // No dictionary — uncompressed V5 store, return raw bytes.
                    Ok(Some(compressed.to_vec()))
                }
            }
        }
    }

    /// Create a reusable zstd decompressor for V5 entries.
    ///
    /// Returns `None` if no dictionary is loaded (uncompressed V5 store).
    /// The caller should create one decompressor and reuse it across many calls
    /// to [`Self::decompress_position_entry`] for best performance.
    pub fn create_v5_decompressor(&self) -> Result<Option<zstd::bulk::Decompressor<'static>>> {
        match &self.v5_zstd_dict {
            Some(dict) => {
                let decompressor =
                    zstd::bulk::Decompressor::with_dictionary(dict).map_err(|e| {
                        DataFusionError::Execution(format!(
                            "failed to create zstd decompressor: {e}"
                        ))
                    })?;
                Ok(Some(decompressor))
            }
            None => Ok(None),
        }
    }

    /// Fetch and decompress a V5 position entry into a reusable buffer.
    ///
    /// Hot-path variant: avoids allocating a new `Vec<u8>` and `Decompressor` per call.
    /// Use [`Self::create_v5_decompressor`] once, then call this in a loop.
    /// Returns `true` if an entry was found (buffer contains decompressed data).
    pub fn get_position_entry_fast(
        &self,
        chrom_code: u16,
        start: i64,
        end: i64,
        decompressor: Option<&mut zstd::bulk::Decompressor<'_>>,
        buf: &mut Vec<u8>,
    ) -> Result<bool> {
        let raw = self.get_position_entry(chrom_code, start, end)?;
        match raw {
            None => Ok(false),
            Some(compressed) => match decompressor {
                Some(dec) => {
                    buf.clear();
                    // Try zero-alloc path first; if buf capacity is too small,
                    // fall back to allocating and swap so future calls succeed.
                    if dec.decompress_to_buffer(&compressed, buf).is_err() {
                        buf.clear();
                        let fresh = dec
                            .decompress(&compressed, compressed.len() * 20 + 4096)
                            .map_err(|e| {
                                DataFusionError::Execution(format!(
                                    "zstd decompression failed: {e}"
                                ))
                            })?;
                        *buf = fresh;
                    }
                    Ok(true)
                }
                None => {
                    buf.clear();
                    buf.extend_from_slice(&compressed);
                    Ok(true)
                }
            },
        }
    }

    /// Batch-insert pre-serialized position entries into the data partition (v5).
    pub(crate) fn batch_insert_position_entries(
        &self,
        entries: &[(Vec<u8>, Vec<u8>)],
    ) -> Result<()> {
        self.ensure_v5()?;
        self.batch_insert_raw(entries)
    }

    fn ensure_v5(&self) -> Result<()> {
        if self.format_version != FORMAT_V5 {
            return Err(DataFusionError::Execution(format!(
                "v5 operation requires format_version={FORMAT_V5}, got {}",
                self.format_version
            )));
        }
        Ok(())
    }

    fn ensure_v1(&self) -> Result<()> {
        if self.format_version != FORMAT_V1 {
            return Err(DataFusionError::Execution(format!(
                "v1 operation requires format_version={FORMAT_V1}, got {}",
                self.format_version
            )));
        }
        Ok(())
    }

    fn ensure_v2(&self) -> Result<()> {
        if self.format_version != FORMAT_V2
            && self.format_version != FORMAT_V3
            && self.format_version != FORMAT_V4
        {
            return Err(DataFusionError::Execution(format!(
                "v2/v3/v4 operation requires format_version={FORMAT_V2}, {FORMAT_V3}, or {FORMAT_V4}, got {}",
                self.format_version
            )));
        }
        if self.v2_blocks.is_none() {
            return Err(DataFusionError::Execution(
                "v2/v3/v4 block store is not initialized".to_string(),
            ));
        }
        Ok(())
    }

    fn ensure_v3(&self) -> Result<()> {
        if self.format_version != FORMAT_V3 {
            return Err(DataFusionError::Execution(format!(
                "v3 operation requires format_version={FORMAT_V3}, got {}",
                self.format_version
            )));
        }
        if self.v2_blocks.is_none() {
            return Err(DataFusionError::Execution(
                "v3 block store is not initialized".to_string(),
            ));
        }
        Ok(())
    }

    fn ensure_v4(&self) -> Result<()> {
        if self.format_version != FORMAT_V4 {
            return Err(DataFusionError::Execution(format!(
                "v4 operation requires format_version={FORMAT_V4}, got {}",
                self.format_version
            )));
        }
        if self.v2_blocks.is_none() {
            return Err(DataFusionError::Execution(
                "v4 block store is not initialized".to_string(),
            ));
        }
        Ok(())
    }

    fn read_format_version(meta: &PartitionHandle) -> Result<u8> {
        let version = match meta
            .get(FORMAT_VERSION_KEY)
            .map_err(|e| DataFusionError::External(Box::new(e)))?
        {
            Some(v) => v[0],
            None => {
                return Err(DataFusionError::Execution(
                    "cache format metadata missing: legacy v0 caches are no longer supported; rebuild cache in v1 format"
                        .to_string(),
                ));
            }
        };
        if version != FORMAT_V1
            && version != FORMAT_V2
            && version != FORMAT_V3
            && version != FORMAT_V4
            && version != FORMAT_V5
        {
            return Err(DataFusionError::Execution(format!(
                "unsupported cache format version {version}: supported versions are {FORMAT_V1}, {FORMAT_V2}, {FORMAT_V3}, {FORMAT_V4}, and {FORMAT_V5}"
            )));
        }
        Ok(version)
    }
}

fn v2_window_ref_key(chrom: &str, window_id: u64) -> Vec<u8> {
    let window_key = encode_window_key(chrom, window_id);
    let mut key = Vec::with_capacity(V2_WINDOW_REF_PREFIX.len() + window_key.len());
    key.extend_from_slice(V2_WINDOW_REF_PREFIX);
    key.extend_from_slice(&window_key);
    key
}

/// Serialize a RecordBatch to Arrow IPC stream format with LZ4 compression.
fn serialize_ipc(batch: &RecordBatch) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    {
        let mut writer = arrow_ipc::writer::StreamWriter::try_new_with_options(
            &mut buf,
            batch.schema_ref(),
            arrow_ipc::writer::IpcWriteOptions::default()
                .try_with_compression(Some(arrow_ipc::CompressionType::LZ4_FRAME))
                .map_err(arrow_err)?,
        )
        .map_err(arrow_err)?;
        writer.write(batch).map_err(arrow_err)?;
        writer.finish().map_err(arrow_err)?;
    }
    Ok(buf)
}

/// Deserialize Arrow IPC stream format into a RecordBatch.
fn deserialize_ipc(data: &[u8]) -> Result<RecordBatch> {
    let cursor = Cursor::new(data);
    let reader = arrow_ipc::reader::StreamReader::try_new(cursor, None).map_err(arrow_err)?;

    let mut batches = Vec::new();
    for batch_result in reader {
        let batch = batch_result.map_err(arrow_err)?;
        batches.push(batch);
    }

    match batches.len() {
        0 => Err(DataFusionError::Execution(
            "IPC data contained no batches".to_string(),
        )),
        1 => Ok(batches.into_iter().next().unwrap()),
        _ => datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches)
            .map_err(arrow_err),
    }
}

/// Serialize schema to IPC bytes by writing an empty batch.
fn schema_to_ipc_bytes(schema: &SchemaRef) -> Result<Vec<u8>> {
    let empty_batch = RecordBatch::new_empty(schema.clone());
    let mut buf = Vec::new();
    {
        let mut writer =
            arrow_ipc::writer::StreamWriter::try_new(&mut buf, schema).map_err(arrow_err)?;
        writer.write(&empty_batch).map_err(arrow_err)?;
        writer.finish().map_err(arrow_err)?;
    }
    Ok(buf)
}

/// Deserialize schema from IPC bytes (empty batch carrying schema).
fn schema_from_ipc_bytes(data: &[u8]) -> Result<SchemaRef> {
    let cursor = Cursor::new(data);
    let reader = arrow_ipc::reader::StreamReader::try_new(cursor, None).map_err(arrow_err)?;
    Ok(reader.schema())
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn test_schema() -> SchemaRef {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ]))
    }

    fn test_batch(schema: &SchemaRef) -> RecordBatch {
        RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(StringArray::from(vec!["rs123", "rs456"])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_ipc_roundtrip() {
        let schema = test_schema();
        let batch = test_batch(&schema);
        let bytes = serialize_ipc(&batch).unwrap();
        let restored = deserialize_ipc(&bytes).unwrap();
        assert_eq!(batch, restored);
    }

    #[test]
    fn test_schema_ipc_roundtrip() {
        let schema = test_schema();
        let bytes = schema_to_ipc_bytes(&schema).unwrap();
        let restored = schema_from_ipc_bytes(&bytes).unwrap();
        assert_eq!(*schema, *restored);
    }

    #[test]
    fn test_open_rejects_unsupported_format() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();

        let store = VepKvStore::create(dir.path(), schema, 1_000_000).unwrap();
        store.put_metadata(FORMAT_VERSION_KEY, &[0]).unwrap();
        store.persist().unwrap();
        drop(store);

        let err = match VepKvStore::open(dir.path()) {
            Ok(_) => panic!("expected open to fail for unsupported format version"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("unsupported cache format version"));
    }

    #[test]
    fn test_kv_store_columnar_roundtrip_v1() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        let store = VepKvStore::create(dir.path(), schema, 1_000_000).unwrap();
        assert_eq!(store.format_version(), FORMAT_V1);

        // Write position index
        let pos_index = crate::position_index::PositionIndex::from_batch(&batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();

        // Write columns
        for (col_idx, field) in batch.schema().fields().iter().enumerate() {
            store
                .put_column(
                    "1",
                    0,
                    to_v1_column_index(col_idx).unwrap(),
                    batch.column(col_idx),
                    field,
                )
                .unwrap();
        }
        store.persist().unwrap();

        // Read back position index
        let loaded_idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(loaded_idx.num_rows(), 2);

        // Read back specific columns
        let cols = store.get_columns("1", 0, &[0, 1, 4]).unwrap().unwrap();
        assert_eq!(cols.len(), 3);
        assert_eq!(cols[0].1.name(), "chrom");
        assert_eq!(cols[1].1.name(), "start");
        assert_eq!(cols[2].1.name(), "allele_string");

        // Missing window returns None
        assert!(store.get_position_index("1", 999).unwrap().is_none());
        assert!(store.get_columns("1", 999, &[0]).unwrap().is_none());
    }

    #[test]
    fn test_put_column_rejects_out_of_range_index() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);
        let store = VepKvStore::create(dir.path(), schema, 1_000_000).unwrap();

        let err = store
            .put_column(
                "1",
                0,
                EntryType::MAX_COLUMN_INDEX + 1,
                batch.column(0),
                batch.schema().field(0),
            )
            .unwrap_err()
            .to_string();
        assert!(err.contains("max encoded index"));
    }

    #[test]
    fn test_get_columns_rejects_out_of_range_index() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let store = VepKvStore::create(dir.path(), schema, 1_000_000).unwrap();

        let err = store
            .get_columns("1", 0, &[EntryType::MAX_COLUMN_INDEX + 1])
            .unwrap_err()
            .to_string();
        assert!(err.contains("max encoded index"));
    }

    #[test]
    fn test_kv_store_window_blocks_roundtrip_v2() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        let store = VepKvStore::create_v2(dir.path(), schema.clone(), 1_000_000).unwrap();
        assert_eq!(store.format_version(), FORMAT_V2);
        assert!(store.v2_block_store().is_some());

        let pos_index = crate::position_index::PositionIndex::from_batch(&batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        let block_ref = store.put_window_block_v2("1", 0, &batch).unwrap();
        assert_eq!(block_ref.row_count, batch.num_rows() as u32);
        store.persist().unwrap();
        drop(store);

        let reopened = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(reopened.format_version(), FORMAT_V2);
        assert!(reopened.v2_block_store().is_some());

        let loaded_idx = reopened.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(loaded_idx.num_rows(), batch.num_rows());

        let loaded_batch = reopened.get_window_block_v2("1", 0).unwrap().unwrap();
        assert_eq!(loaded_batch, batch);

        let selected = reopened
            .get_window_columns_v2("1", 0, &[0, 4])
            .unwrap()
            .unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].to_data(), batch.column(0).to_data());
        assert_eq!(selected[1].to_data(), batch.column(4).to_data());
    }

    #[test]
    fn test_v1_column_api_rejected_for_v2_store() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let store = VepKvStore::create_v2(dir.path(), schema, 1_000_000).unwrap();

        let err = store.get_columns("1", 0, &[0]).unwrap_err().to_string();
        assert!(err.contains("v1 operation requires format_version"));
    }

    #[test]
    fn test_kv_store_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let store = VepKvStore::create(dir.path(), schema, 1_000_000).unwrap();

        store.put_metadata(b"test_key", b"test_value").unwrap();
        let value = store.get_metadata(b"test_key").unwrap().unwrap();
        assert_eq!(value, b"test_value");

        assert!(store.get_metadata(b"missing").unwrap().is_none());
    }

    #[test]
    fn test_kv_store_window_blocks_roundtrip_v3() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        let store = VepKvStore::create_v3(dir.path(), schema.clone(), 1_000_000).unwrap();
        assert_eq!(store.format_version(), FORMAT_V3);
        assert!(store.v2_block_store().is_some());

        let pos_index = crate::position_index::PositionIndex::from_batch(&batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        let block_ref = store.put_window_block_v3("1", 0, &batch).unwrap();
        assert_eq!(block_ref.row_count, batch.num_rows() as u32);
        store.persist().unwrap();
        drop(store);

        let reopened = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(reopened.format_version(), FORMAT_V3);
        assert!(reopened.v2_block_store().is_some());

        let loaded_idx = reopened.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(loaded_idx.num_rows(), batch.num_rows());

        let loaded_batch = reopened.get_window_batch_v3("1", 0).unwrap().unwrap();
        assert_eq!(loaded_batch, batch);

        let selected = reopened
            .get_window_columns_v3("1", 0, &[0, 4])
            .unwrap()
            .unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].to_data(), batch.column(0).to_data());
        assert_eq!(selected[1].to_data(), batch.column(4).to_data());
    }

    #[test]
    fn test_kv_store_window_blocks_roundtrip_v3_lz4() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        let store = VepKvStore::create_v3_with_codec(
            dir.path(),
            schema.clone(),
            1_000_000,
            WindowBlockCodec::Lz4,
        )
        .unwrap();
        assert_eq!(store.format_version(), FORMAT_V3);

        let pos_index = crate::position_index::PositionIndex::from_batch(&batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        store.put_window_block_v3("1", 0, &batch).unwrap();
        store.persist().unwrap();
        drop(store);

        let reopened = VepKvStore::open(dir.path()).unwrap();
        let loaded_batch = reopened.get_window_batch_v3("1", 0).unwrap().unwrap();
        assert_eq!(loaded_batch, batch);
    }

    #[test]
    fn test_kv_store_window_blocks_roundtrip_v4() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        let store = VepKvStore::create_v4(dir.path(), schema.clone(), 1_000_000).unwrap();
        assert_eq!(store.format_version(), FORMAT_V4);
        assert!(store.v2_block_store().is_some());

        let pos_index = crate::position_index::PositionIndex::from_batch(&batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        let block_ref = store.put_window_block_v4("1", 0, &batch).unwrap();
        assert_eq!(block_ref.row_count, batch.num_rows() as u32);
        store.persist().unwrap();
        drop(store);

        let reopened = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(reopened.format_version(), FORMAT_V4);
        assert!(reopened.v2_block_store().is_some());

        let loaded_idx = reopened.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(loaded_idx.num_rows(), batch.num_rows());

        let loaded_batch = reopened.get_window_batch_v4("1", 0).unwrap().unwrap();
        assert_eq!(loaded_batch, batch);

        let selected = reopened
            .get_window_columns_v4("1", 0, &[0, 4])
            .unwrap()
            .unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].to_data(), batch.column(0).to_data());
        assert_eq!(selected[1].to_data(), batch.column(4).to_data());
    }

    #[test]
    fn test_kv_store_window_blocks_roundtrip_v4_lz4() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        let store = VepKvStore::create_v4_with_codec(
            dir.path(),
            schema.clone(),
            1_000_000,
            WindowBlockCodec::Lz4,
        )
        .unwrap();
        assert_eq!(store.format_version(), FORMAT_V4);

        let pos_index = crate::position_index::PositionIndex::from_batch(&batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        store.put_window_block_v4("1", 0, &batch).unwrap();
        store.persist().unwrap();
        drop(store);

        let reopened = VepKvStore::open(dir.path()).unwrap();
        let loaded_batch = reopened.get_window_batch_v4("1", 0).unwrap().unwrap();
        assert_eq!(loaded_batch, batch);

        let cols = reopened
            .get_window_columns_v4("1", 0, &[0, 3])
            .unwrap()
            .unwrap();
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].len(), 2);
    }

    #[test]
    fn test_kv_store_v5_put_get() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();

        let store = VepKvStore::create_v5(dir.path(), schema.clone()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V5);

        let value = b"test_position_data";
        store.put_position_entry("1", 100, 100, value).unwrap();
        store.persist().unwrap();

        let chrom_code = crate::key_encoding::chrom_to_code("1");
        let loaded = store
            .get_position_entry(chrom_code, 100, 100)
            .unwrap()
            .unwrap();
        assert_eq!(&*loaded, value);

        // Missing position returns None.
        let missing = store.get_position_entry(chrom_code, 999, 999).unwrap();
        assert!(missing.is_none());

        // Reopen and verify persistence.
        drop(store);
        let reopened = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(reopened.format_version(), FORMAT_V5);
        let loaded = reopened
            .get_position_entry(chrom_code, 100, 100)
            .unwrap()
            .unwrap();
        assert_eq!(&*loaded, value);
    }

    #[test]
    fn test_v5_zstd_dict_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();

        let store = VepKvStore::create_v5(dir.path(), schema.clone()).unwrap();
        assert!(!store.has_v5_dict());

        // Generate enough samples to train a dictionary.
        let mut samples: Vec<Vec<u8>> = Vec::new();
        for i in 0..200u32 {
            let data = format!("variation_name_{i}\tallele_A/G_{i}\textra_data_{i}");
            samples.push(data.into_bytes());
        }

        // Train dictionary.
        let mut continuous = Vec::new();
        let mut sizes = Vec::with_capacity(samples.len());
        for s in &samples {
            continuous.extend_from_slice(s);
            sizes.push(s.len());
        }
        let dict = zstd::dict::from_continuous(&continuous, &sizes, 16 * 1024).unwrap();
        store.store_v5_dict(&dict).unwrap();

        // Compress and store entries.
        let mut compressor = zstd::bulk::Compressor::with_dictionary(3, &dict).unwrap();
        for (i, sample) in samples.iter().enumerate() {
            let compressed = compressor.compress(sample).unwrap();
            store
                .put_position_entry("1", i as i64, i as i64, &compressed)
                .unwrap();
        }
        store.persist().unwrap();
        drop(store);

        // Reopen — dictionary should be loaded.
        let reopened = VepKvStore::open(dir.path()).unwrap();
        assert!(reopened.has_v5_dict());

        // Verify decompressed roundtrip for each entry.
        let chrom_code = crate::key_encoding::chrom_to_code("1");
        for (i, sample) in samples.iter().enumerate() {
            let decompressed = reopened
                .get_position_entry_decompressed(chrom_code, i as i64, i as i64)
                .unwrap()
                .unwrap();
            assert_eq!(&decompressed, sample, "mismatch at entry {i}");
        }

        // Missing entry still returns None.
        assert!(
            reopened
                .get_position_entry_decompressed(chrom_code, 9999, 9999)
                .unwrap()
                .is_none()
        );
    }
}
