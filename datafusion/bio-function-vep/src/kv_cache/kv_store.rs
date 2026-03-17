//! fjall KV store for position-keyed VEP cache entries.

use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use fjall::{Database, Keyspace, KeyspaceCreateOptions, PersistMode};

const META_KEYSPACE: &str = "meta";
const DATA_KEYSPACE: &str = "data";
const SCHEMA_KEY: &[u8] = b"schema";
const FORMAT_VERSION_KEY: &[u8] = b"format_version";
const ZSTD_DICT_KEY: &[u8] = b"zstd_dict";
const ZSTD_LEVEL_KEY: &[u8] = b"zstd_level";
const MIN_DECOMPRESS_CAPACITY: usize = 4 * 1024;
const MAX_DECOMPRESSED_ENTRY_BYTES: usize = 2 * 1024 * 1024 * 1024;

/// Format version 0: position-keyed entries — one fjall entry per genomic position.
///
/// Each `(chrom, start, end)` maps to a zstd-compressed serialized entry containing
/// the allele table and column-major data. Only positions matching VCF queries are fetched.
pub const FORMAT_V0: u8 = 0;

/// Wrapper around fjall `Database` for VEP cache storage.
///
/// Layout:
/// - `meta` keyspace: schema (Arrow IPC), format version, zstd dictionary
/// - `data` keyspace: `(chrom, start, end)` -> compressed position entry
pub struct VepKvStore {
    root_path: PathBuf,
    db: Database,
    data: Keyspace,
    meta: Keyspace,
    schema: SchemaRef,
    format_version: u8,
    /// Raw zstd dictionary bytes for position entry decompression.
    zstd_dict: Option<Arc<Vec<u8>>>,
}

fn arrow_err(e: datafusion::arrow::error::ArrowError) -> DataFusionError {
    DataFusionError::ArrowError(Box::new(e), None)
}

fn fjall_err(e: fjall::Error) -> DataFusionError {
    DataFusionError::External(Box::new(e))
}

fn is_destination_too_small_msg(msg: &str) -> bool {
    msg.contains("Destination buffer is too small")
}

fn next_capacity(current: usize) -> Option<usize> {
    if current >= MAX_DECOMPRESSED_ENTRY_BYTES {
        None
    } else {
        current
            .checked_mul(2)
            .map(|v| v.min(MAX_DECOMPRESSED_ENTRY_BYTES))
    }
}

/// Decompress a zstd payload into a reusable buffer, growing capacity as needed.
pub(crate) fn decompress_into_buffer_with_retry(
    decompressor: &mut zstd::bulk::Decompressor<'_>,
    compressed: &[u8],
    out: &mut Vec<u8>,
    context: &str,
) -> Result<()> {
    let mut target_capacity = compressed
        .len()
        .saturating_mul(16)
        .max(MIN_DECOMPRESS_CAPACITY)
        .min(MAX_DECOMPRESSED_ENTRY_BYTES);

    if out.capacity() < target_capacity {
        out.reserve(target_capacity - out.capacity());
    }

    loop {
        out.clear();
        match decompressor.decompress_to_buffer(compressed, out) {
            Ok(_) => return Ok(()),
            Err(e) => {
                let err_msg = e.to_string();
                if !is_destination_too_small_msg(&err_msg) {
                    return Err(DataFusionError::Execution(format!("{context}: {e}")));
                }
                let Some(next) = next_capacity(target_capacity) else {
                    return Err(DataFusionError::Execution(format!(
                        "{context}: destination buffer remains too small at {} bytes (compressed={} bytes)",
                        target_capacity,
                        compressed.len()
                    )));
                };
                target_capacity = next;
                if out.capacity() < target_capacity {
                    out.reserve(target_capacity - out.capacity());
                }
            }
        }
    }
}

impl VepKvStore {
    /// Access the underlying fjall Database handle.
    pub fn database(&self) -> &Database {
        &self.db
    }

    /// Open an existing KV store with default settings (256 MB block cache).
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_cache_size(path, 256 * 1024 * 1024)
    }

    /// Open an existing KV store with a custom fjall block cache size (bytes).
    pub fn open_with_cache_size(path: impl AsRef<Path>, cache_size_bytes: u64) -> Result<Self> {
        let root_path = path.as_ref().to_path_buf();
        let db = Database::builder(&root_path)
            .cache_size(cache_size_bytes)
            .open()
            .map_err(fjall_err)?;

        let data = db
            .keyspace(DATA_KEYSPACE, KeyspaceCreateOptions::default)
            .map_err(fjall_err)?;
        let meta = db
            .keyspace(META_KEYSPACE, KeyspaceCreateOptions::default)
            .map_err(fjall_err)?;

        let schema = Self::read_schema(&meta)?;
        let format_version = Self::read_format_version(&meta)?;

        // Load zstd dictionary if present.
        let zstd_dict = meta
            .get(ZSTD_DICT_KEY)
            .map_err(fjall_err)?
            .map(|raw| Arc::new(raw.to_vec()));

        Ok(Self {
            root_path,
            db,
            data,
            meta,
            schema,
            format_version,
            zstd_dict,
        })
    }

    /// Create a new KV store with the given schema.
    ///
    /// Uses write-optimized settings for bulk loading:
    /// - `manual_journal_persist` — skip per-batch fsync (persist once at the end)
    /// - L0 threshold raised to 16 — defer compaction during ingestion
    /// - LZ4 disabled on data blocks — values are already zstd-compressed
    pub fn create(path: impl AsRef<Path>, schema: SchemaRef) -> Result<Self> {
        let root_path = path.as_ref().to_path_buf();
        let db = Database::builder(&root_path)
            .cache_size(256 * 1024 * 1024)
            .manual_journal_persist(true)
            .open()
            .map_err(fjall_err)?;

        // Data keyspace: write-optimized for bulk loading.
        let data_opts = || {
            KeyspaceCreateOptions::default()
                .manual_journal_persist(true)
                .compaction_strategy(Arc::new(
                    fjall::compaction::Leveled::default().with_l0_threshold(16),
                ))
                .data_block_compression_policy(fjall::config::CompressionPolicy::disabled())
        };

        let data = db.keyspace(DATA_KEYSPACE, data_opts).map_err(fjall_err)?;
        let meta = db
            .keyspace(META_KEYSPACE, KeyspaceCreateOptions::default)
            .map_err(fjall_err)?;

        let schema_bytes = schema_to_ipc_bytes(&schema)?;
        meta.insert(SCHEMA_KEY, &schema_bytes).map_err(fjall_err)?;

        meta.insert(FORMAT_VERSION_KEY, [FORMAT_V0])
            .map_err(fjall_err)?;

        db.persist(PersistMode::SyncAll).map_err(fjall_err)?;

        Ok(Self {
            root_path,
            db,
            data,
            meta,
            schema,
            format_version: FORMAT_V0,
            zstd_dict: None,
        })
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn root_path(&self) -> &Path {
        &self.root_path
    }

    /// Get the format version of this store.
    pub fn format_version(&self) -> u8 {
        self.format_version
    }

    /// Return a reference to the data keyspace handle (for diagnostics).
    pub fn data_partition(&self) -> &Keyspace {
        &self.data
    }

    // --- position-keyed API ---

    /// Write a position entry.
    pub fn put_position_entry(
        &self,
        chrom: &str,
        start: i64,
        end: i64,
        value: &[u8],
    ) -> Result<()> {
        let key = super::key_encoding::encode_position_key(chrom, start, end);
        self.data.insert(&key, value).map_err(fjall_err)?;
        Ok(())
    }

    /// Read a position entry by chrom code + coordinates.
    ///
    /// Uses a pre-encoded key buffer for hot-path efficiency.
    pub fn get_position_entry(
        &self,
        chrom_code: u16,
        start: i64,
        end: i64,
    ) -> Result<Option<fjall::UserValue>> {
        let mut key_buf = Vec::with_capacity(18);
        super::key_encoding::encode_position_key_buf(chrom_code, start, end, &mut key_buf);
        match self.data.get(&key_buf) {
            Ok(v) => Ok(v),
            Err(e) => Err(fjall_err(e)),
        }
    }

    /// Store a trained zstd dictionary for position entries.
    pub fn store_dict(&self, dict_bytes: &[u8]) -> Result<()> {
        self.meta
            .insert(ZSTD_DICT_KEY, dict_bytes)
            .map_err(fjall_err)?;
        Ok(())
    }

    /// Store the zstd compression level used when building this cache (informational).
    pub fn store_zstd_level(&self, level: i32) -> Result<()> {
        self.meta
            .insert(ZSTD_LEVEL_KEY, level.to_le_bytes())
            .map_err(fjall_err)?;
        Ok(())
    }

    /// Whether this store has a zstd dictionary loaded.
    pub fn has_dict(&self) -> bool {
        self.zstd_dict.is_some()
    }

    /// Read and decompress a position entry using the zstd dictionary.
    ///
    /// If no dictionary is loaded (legacy uncompressed store), returns the
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
                if let Some(dict) = &self.zstd_dict {
                    let mut decompressor = zstd::bulk::Decompressor::with_dictionary(dict)
                        .map_err(|e| {
                            DataFusionError::Execution(format!(
                                "failed to create zstd decompressor: {e}"
                            ))
                        })?;
                    let mut decompressed = Vec::with_capacity(
                        compressed
                            .len()
                            .saturating_mul(16)
                            .max(MIN_DECOMPRESS_CAPACITY),
                    );
                    decompress_into_buffer_with_retry(
                        &mut decompressor,
                        &compressed,
                        &mut decompressed,
                        "zstd decompression failed",
                    )?;
                    Ok(Some(decompressed))
                } else {
                    Ok(Some(compressed.to_vec()))
                }
            }
        }
    }

    /// Create a reusable zstd decompressor for position entries.
    ///
    /// Returns `None` if no dictionary is loaded (uncompressed store).
    /// The caller should create one decompressor and reuse it across many calls
    /// for best performance.
    pub fn create_decompressor(&self) -> Result<Option<zstd::bulk::Decompressor<'static>>> {
        match &self.zstd_dict {
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

    /// Fetch and decompress a position entry into a reusable buffer.
    ///
    /// Hot-path variant: avoids allocating a new `Vec<u8>` and `Decompressor` per call.
    /// Use [`Self::create_decompressor`] once, then call this in a loop.
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
                    decompress_into_buffer_with_retry(
                        dec,
                        &compressed,
                        buf,
                        "zstd decompression failed",
                    )?;
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

    /// Batch-insert pre-serialized position entries into the data keyspace.
    pub(crate) fn batch_insert_position_entries(
        &self,
        entries: &[(Vec<u8>, Vec<u8>)],
    ) -> Result<()> {
        self.batch_insert_raw(entries)
    }

    /// Persist all data to disk.
    pub fn persist(&self) -> Result<()> {
        self.db.persist(PersistMode::SyncAll).map_err(fjall_err)?;
        Ok(())
    }

    /// Atomically insert a batch of pre-serialized (key, value) pairs into the data keyspace.
    pub(crate) fn batch_insert_raw(&self, entries: &[(Vec<u8>, Vec<u8>)]) -> Result<()> {
        if entries.is_empty() {
            return Ok(());
        }
        let mut batch = fjall::OwnedWriteBatch::with_capacity(self.db.clone(), entries.len());
        for (key, value) in entries {
            batch.insert(&self.data, key, value);
        }
        batch.commit().map_err(fjall_err)?;
        Ok(())
    }

    /// Store arbitrary metadata.
    pub fn put_metadata(&self, key: &[u8], value: &[u8]) -> Result<()> {
        self.meta.insert(key, value).map_err(fjall_err)?;
        Ok(())
    }

    /// Read arbitrary metadata.
    pub fn get_metadata(&self, key: &[u8]) -> Result<Option<Vec<u8>>> {
        match self.meta.get(key) {
            Ok(Some(v)) => Ok(Some(v.to_vec())),
            Ok(None) => Ok(None),
            Err(e) => Err(fjall_err(e)),
        }
    }

    fn read_schema(meta: &Keyspace) -> Result<SchemaRef> {
        let value = meta.get(SCHEMA_KEY).map_err(fjall_err)?.ok_or_else(|| {
            DataFusionError::Execution("KV store missing schema metadata".to_string())
        })?;
        schema_from_ipc_bytes(&value)
    }

    fn read_format_version(meta: &Keyspace) -> Result<u8> {
        let version = match meta.get(FORMAT_VERSION_KEY).map_err(fjall_err)? {
            Some(v) => v[0],
            None => {
                return Err(DataFusionError::Execution(
                    "cache format metadata missing".to_string(),
                ));
            }
        };
        if version != FORMAT_V0 {
            return Err(DataFusionError::Execution(format!(
                "unsupported cache format version {version}: only version {FORMAT_V0} is supported"
            )));
        }
        Ok(version)
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

        let store = VepKvStore::create(dir.path(), schema).unwrap();
        store.put_metadata(FORMAT_VERSION_KEY, &[99]).unwrap();
        store.persist().unwrap();
        drop(store);

        let err = match VepKvStore::open(dir.path()) {
            Ok(_) => panic!("expected open to fail for unsupported format version"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("unsupported cache format version"));
    }

    #[test]
    fn test_kv_store_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();
        let store = VepKvStore::create(dir.path(), schema).unwrap();

        store.put_metadata(b"test_key", b"test_value").unwrap();
        let value = store.get_metadata(b"test_key").unwrap().unwrap();
        assert_eq!(value, b"test_value");

        assert!(store.get_metadata(b"missing").unwrap().is_none());
    }

    #[test]
    fn test_kv_store_put_get() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();

        let store = VepKvStore::create(dir.path(), schema.clone()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V0);

        let value = b"test_position_data";
        store.put_position_entry("1", 100, 100, value).unwrap();
        store.persist().unwrap();

        let chrom_code = crate::kv_cache::key_encoding::chrom_to_code("1");
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
        assert_eq!(reopened.format_version(), FORMAT_V0);
        let loaded = reopened
            .get_position_entry(chrom_code, 100, 100)
            .unwrap()
            .unwrap();
        assert_eq!(&*loaded, value);
    }

    #[test]
    fn test_zstd_dict_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();

        let store = VepKvStore::create(dir.path(), schema.clone()).unwrap();
        assert!(!store.has_dict());

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
        store.store_dict(&dict).unwrap();

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
        assert!(reopened.has_dict());

        // Verify decompressed roundtrip for each entry.
        let chrom_code = crate::kv_cache::key_encoding::chrom_to_code("1");
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

    #[test]
    fn test_zstd_decompression_retries_when_destination_buffer_too_small() {
        let dir = tempfile::tempdir().unwrap();
        let schema = test_schema();

        let store = VepKvStore::create(dir.path(), schema).unwrap();

        // Train a tiny dictionary to exercise dictionary-based decompression.
        let mut samples: Vec<Vec<u8>> = Vec::new();
        for i in 0..256u32 {
            samples.push(format!("sample_{i}_allele_A/G").into_bytes());
        }
        let mut continuous = Vec::new();
        let mut sizes = Vec::with_capacity(samples.len());
        for s in &samples {
            continuous.extend_from_slice(s);
            sizes.push(s.len());
        }
        let dict = zstd::dict::from_continuous(&continuous, &sizes, 8 * 1024).unwrap();
        store.store_dict(&dict).unwrap();

        // Highly compressible payload: very large output relative to compressed bytes.
        let huge = vec![b'X'; 512 * 1024];
        let mut compressor = zstd::bulk::Compressor::with_dictionary(3, &dict).unwrap();
        let compressed = compressor.compress(&huge).unwrap();

        // Guard the regression condition: initial 16x guess is insufficient.
        assert!(
            huge.len() > compressed.len() * 16,
            "test payload is not compressed enough to trigger resize retry"
        );

        store.put_position_entry("1", 42, 42, &compressed).unwrap();
        store.persist().unwrap();
        drop(store);

        let reopened = VepKvStore::open(dir.path()).unwrap();
        let chrom_code = crate::kv_cache::key_encoding::chrom_to_code("1");
        let decompressed = reopened
            .get_position_entry_decompressed(chrom_code, 42, 42)
            .unwrap()
            .unwrap();
        assert_eq!(decompressed, huge);
    }
}
