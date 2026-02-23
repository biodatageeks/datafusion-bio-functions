//! fjall KV store wrapper for reading/writing Arrow IPC batches per window.

use std::io::Cursor;
use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{ArrayRef, RecordBatch};
use datafusion::arrow::datatypes::{Field, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use fjall::{Config, Keyspace, PartitionCreateOptions, PartitionHandle, PersistMode};

use crate::key_encoding::{EntryType, decode_window_key, encode_entry_key, encode_window_key};
use crate::position_index::PositionIndex;

const META_PARTITION: &str = "meta";
const DATA_PARTITION: &str = "data";
const SCHEMA_KEY: &[u8] = b"schema";
const WINDOW_SIZE_KEY: &[u8] = b"window_size";
const FORMAT_VERSION_KEY: &[u8] = b"format_version";

/// Format version 0: legacy monolithic (all columns in one Arrow IPC blob per window).
pub const FORMAT_V0: u8 = 0;
/// Format version 1: columnar (position index + per-column Arrow IPC entries).
pub const FORMAT_V1: u8 = 1;

/// Wrapper around fjall `Keyspace` for VEP cache storage.
///
/// Layout:
/// - `meta` partition: schema (Arrow IPC), window size, coordinate metadata
/// - `data` partition: `(chrom, window_id)` -> Arrow IPC serialized `RecordBatch`
pub struct VepKvStore {
    keyspace: Keyspace,
    data: PartitionHandle,
    meta: PartitionHandle,
    schema: SchemaRef,
    window_size: u64,
    format_version: u8,
}

fn arrow_err(e: datafusion::arrow::error::ArrowError) -> DataFusionError {
    DataFusionError::ArrowError(Box::new(e), None)
}

impl VepKvStore {
    /// Open an existing KV store with default settings (256 MB block cache).
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_cache_size(path, 256 * 1024 * 1024)
    }

    /// Open an existing KV store with a custom fjall block cache size (bytes).
    pub fn open_with_cache_size(path: impl AsRef<Path>, cache_size_bytes: u64) -> Result<Self> {
        let keyspace = Config::new(path)
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

        Ok(Self {
            keyspace,
            data,
            meta,
            schema,
            window_size,
            format_version,
        })
    }

    /// Create a new KV store with the given schema, window size, and format version.
    ///
    /// Use `FORMAT_V1` for new caches (columnar), `FORMAT_V0` for legacy.
    pub fn create(path: impl AsRef<Path>, schema: SchemaRef, window_size: u64) -> Result<Self> {
        Self::create_with_version(path, schema, window_size, FORMAT_V1)
    }

    /// Create a new KV store with explicit format version.
    pub fn create_with_version(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        format_version: u8,
    ) -> Result<Self> {
        Self::create_with_options(
            path,
            schema,
            window_size,
            format_version,
            256 * 1024 * 1024,
            PartitionCreateOptions::default(),
        )
    }

    /// Create a new KV store with full control over fjall tuning.
    pub fn create_with_options(
        path: impl AsRef<Path>,
        schema: SchemaRef,
        window_size: u64,
        format_version: u8,
        cache_size_bytes: u64,
        partition_opts: PartitionCreateOptions,
    ) -> Result<Self> {
        let keyspace = Config::new(path)
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

        keyspace
            .persist(PersistMode::SyncAll)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        Ok(Self {
            keyspace,
            data,
            meta,
            schema,
            window_size,
            format_version,
        })
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn window_size(&self) -> u64 {
        self.window_size
    }

    /// Get the format version of this store.
    pub fn format_version(&self) -> u8 {
        self.format_version
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
        let mut results = Vec::with_capacity(col_indices.len());
        for &col_idx in col_indices {
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

    // --- Legacy v0 API ---

    /// Read a window's data as a RecordBatch (v0 legacy format).
    pub fn get_window(&self, chrom: &str, window_id: u64) -> Result<Option<RecordBatch>> {
        let key = encode_window_key(chrom, window_id);
        match self.data.get(&key) {
            Ok(Some(value)) => {
                let batch = deserialize_ipc(&value)?;
                Ok(Some(batch))
            }
            Ok(None) => Ok(None),
            Err(e) => Err(DataFusionError::External(Box::new(e))),
        }
    }

    /// Write a window's data as an Arrow IPC batch.
    pub fn put_window(&self, chrom: &str, window_id: u64, batch: &RecordBatch) -> Result<()> {
        let key = encode_window_key(chrom, window_id);
        let value = serialize_ipc(batch)?;
        self.data
            .insert(&key, &value)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        Ok(())
    }

    /// Persist all data to disk.
    pub fn persist(&self) -> Result<()> {
        self.keyspace
            .persist(PersistMode::SyncAll)
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
    /// For v0 format, returns full monolithic batches.
    /// For v1 format, skips non-legacy entries (position index + columns are separate).
    pub fn windows_for_chrom(&self, chrom: &str) -> Result<Vec<(String, u64, RecordBatch)>> {
        let prefix = crate::key_encoding::chrom_prefix(chrom);
        let mut results = Vec::new();
        for item in self.data.prefix(&prefix) {
            let (key, value) = item.map_err(|e| DataFusionError::External(Box::new(e)))?;
            if key.len() == 10 {
                // v0 legacy 10-byte key
                let (c, wid) = decode_window_key(&key);
                let batch = deserialize_ipc(&value)?;
                results.push((c, wid, batch));
            }
            // v1 keys (11 bytes) are skipped — use get_position_index + get_columns instead
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

    fn read_format_version(meta: &PartitionHandle) -> Result<u8> {
        match meta
            .get(FORMAT_VERSION_KEY)
            .map_err(|e| DataFusionError::External(Box::new(e)))?
        {
            Some(v) => Ok(v[0]),
            None => Ok(FORMAT_V0), // Legacy caches have no version key
        }
    }
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
    fn test_kv_store_roundtrip_v0() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        {
            let store =
                VepKvStore::create_with_version(dir.path(), schema.clone(), 1_000_000, FORMAT_V0)
                    .unwrap();
            store.put_window("1", 0, &batch).unwrap();
            store.persist().unwrap();
        }

        {
            let store = VepKvStore::open(dir.path()).unwrap();
            assert_eq!(store.format_version(), FORMAT_V0);
            assert_eq!(store.schema().as_ref(), schema.as_ref());
            assert_eq!(store.window_size(), 1_000_000);

            let retrieved = store.get_window("1", 0).unwrap().unwrap();
            assert_eq!(batch, retrieved);

            assert!(store.get_window("1", 999).unwrap().is_none());
        }
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
                .put_column("1", 0, col_idx as u8, batch.column(col_idx), field)
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
    fn test_kv_store_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let store = VepKvStore::create(dir.path(), schema, 1_000_000).unwrap();

        store.put_metadata(b"test_key", b"test_value").unwrap();
        let value = store.get_metadata(b"test_key").unwrap().unwrap();
        assert_eq!(value, b"test_value");

        assert!(store.get_metadata(b"missing").unwrap().is_none());
    }
}
