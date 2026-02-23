//! fjall KV store wrapper for reading/writing Arrow IPC batches per window.

use std::io::Cursor;
use std::path::Path;

use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use fjall::{Config, Keyspace, PartitionCreateOptions, PartitionHandle, PersistMode};

use crate::key_encoding::{decode_window_key, encode_window_key};

const META_PARTITION: &str = "meta";
const DATA_PARTITION: &str = "data";
const SCHEMA_KEY: &[u8] = b"schema";
const WINDOW_SIZE_KEY: &[u8] = b"window_size";

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
}

fn arrow_err(e: datafusion::arrow::error::ArrowError) -> DataFusionError {
    DataFusionError::ArrowError(Box::new(e), None)
}

impl VepKvStore {
    /// Open an existing KV store.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let keyspace = Config::new(path)
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

        Ok(Self {
            keyspace,
            data,
            meta,
            schema,
            window_size,
        })
    }

    /// Create a new KV store with the given schema and window size.
    pub fn create(path: impl AsRef<Path>, schema: SchemaRef, window_size: u64) -> Result<Self> {
        let keyspace = Config::new(path)
            .open()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let data = keyspace
            .open_partition(DATA_PARTITION, PartitionCreateOptions::default())
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

        keyspace
            .persist(PersistMode::SyncAll)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        Ok(Self {
            keyspace,
            data,
            meta,
            schema,
            window_size,
        })
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn window_size(&self) -> u64 {
        self.window_size
    }

    /// Read a window's data as a RecordBatch.
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
    pub fn windows_for_chrom(&self, chrom: &str) -> Result<Vec<(String, u64, RecordBatch)>> {
        let prefix = crate::key_encoding::chrom_prefix(chrom);
        let mut results = Vec::new();
        for item in self.data.prefix(&prefix) {
            let (key, value) = item.map_err(|e| DataFusionError::External(Box::new(e)))?;
            let (c, wid) = decode_window_key(&key);
            let batch = deserialize_ipc(&value)?;
            results.push((c, wid, batch));
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
    fn test_kv_store_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let batch = test_batch(&schema);

        {
            let store = VepKvStore::create(dir.path(), schema.clone(), 1_000_000).unwrap();
            store.put_window("1", 0, &batch).unwrap();
            store.persist().unwrap();
        }

        {
            let store = VepKvStore::open(dir.path()).unwrap();
            assert_eq!(store.schema().as_ref(), schema.as_ref());
            assert_eq!(store.window_size(), 1_000_000);

            let retrieved = store.get_window("1", 0).unwrap().unwrap();
            assert_eq!(batch, retrieved);

            assert!(store.get_window("1", 999).unwrap().is_none());
        }
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
