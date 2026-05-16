//! redb KV cache backend for VEP variant lookup.
//!
//! The on-disk value format matches the fjall-backed store: one serialized
//! position entry per `(chrom, start)` key, with optional zstd dictionary
//! compression recorded in metadata. Annotation opens redb through
//! `ReadOnlyDatabase` so separate annotator processes can read the same cache
//! concurrently.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use async_trait::async_trait;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::Expr;
use redb::{
    Database, ReadOnlyDatabase, ReadTransaction, ReadableDatabase, ReadableTable, TableDefinition,
};

use super::cache_exec::PositionEntryStore;
use super::kv_store::{
    FORMAT_V0, decompress_into_buffer_with_retry, schema_from_ipc_bytes, schema_to_ipc_bytes,
};

const META_TABLE: TableDefinition<&str, &[u8]> = TableDefinition::new("meta");
const DATA_TABLE: TableDefinition<&[u8], &[u8]> = TableDefinition::new("data");

const SCHEMA_KEY: &str = "schema";
const FORMAT_VERSION_KEY: &str = "format_version";
const ZSTD_DICT_KEY: &str = "zstd_dict";
const ZSTD_LEVEL_KEY: &str = "zstd_level";
const CONTIG_CODES_KEY: &str = "contig_codes";

enum RedbHandle {
    Writable(Mutex<Database>),
    ReadOnly(ReadOnlyDatabase),
}

type RedbMetadata = (SchemaRef, u8, Option<Arc<Vec<u8>>>);

/// redb-backed variant cache store.
///
/// Use [`Self::create`] while building caches. Use [`Self::open`] for
/// annotation/lookups; it opens a process-safe read-only database handle.
pub struct VepRedbStore {
    root_path: PathBuf,
    db: RedbHandle,
    schema: SchemaRef,
    format_version: u8,
    zstd_dict: Option<Arc<Vec<u8>>>,
    write_lock: Mutex<()>,
}

fn redb_err(e: impl std::fmt::Display) -> DataFusionError {
    DataFusionError::Execution(format!("redb error: {e}"))
}

impl VepRedbStore {
    /// Open an existing redb cache for annotation using `ReadOnlyDatabase`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let root_path = path.as_ref().to_path_buf();
        let db = ReadOnlyDatabase::open(&root_path).map_err(redb_err)?;
        let (schema, format_version, zstd_dict) = Self::read_metadata(&db)?;
        Ok(Self {
            root_path,
            db: RedbHandle::ReadOnly(db),
            schema,
            format_version,
            zstd_dict,
            write_lock: Mutex::new(()),
        })
    }

    /// Create or open a writable redb cache and initialize metadata.
    pub fn create(path: impl AsRef<Path>, schema: SchemaRef) -> Result<Self> {
        let root_path = path.as_ref().to_path_buf();
        let db = Database::create(&root_path).map_err(redb_err)?;

        {
            let write_txn = db.begin_write().map_err(redb_err)?;
            {
                let mut meta = write_txn.open_table(META_TABLE).map_err(redb_err)?;
                let schema_bytes = schema_to_ipc_bytes(&schema)?;
                meta.insert(SCHEMA_KEY, schema_bytes.as_slice())
                    .map_err(redb_err)?;
                meta.insert(FORMAT_VERSION_KEY, [FORMAT_V0].as_slice())
                    .map_err(redb_err)?;
            }
            {
                let _data = write_txn.open_table(DATA_TABLE).map_err(redb_err)?;
            }
            write_txn.commit().map_err(redb_err)?;
        }

        Ok(Self {
            root_path,
            db: RedbHandle::Writable(Mutex::new(db)),
            schema,
            format_version: FORMAT_V0,
            zstd_dict: None,
            write_lock: Mutex::new(()),
        })
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn root_path(&self) -> &Path {
        &self.root_path
    }

    pub fn format_version(&self) -> u8 {
        self.format_version
    }

    pub fn has_dict(&self) -> bool {
        self.zstd_dict.is_some()
    }

    pub fn put_position_entry(&self, chrom: &str, start: i64, value: &[u8]) -> Result<()> {
        let key = super::key_encoding::encode_position_key(chrom, start);
        self.put_position_entry_raw(&key, value)
    }

    pub(crate) fn put_position_entry_raw(&self, key: &[u8], value: &[u8]) -> Result<()> {
        let _guard = self
            .write_lock
            .lock()
            .map_err(|e| DataFusionError::Execution(format!("redb write mutex poisoned: {e}")))?;
        let RedbHandle::Writable(db) = &self.db else {
            return Err(DataFusionError::Execution(
                "redb cache was opened read-only; writes are not allowed".to_string(),
            ));
        };
        let db = db.lock().map_err(|e| {
            DataFusionError::Execution(format!("redb database mutex poisoned: {e}"))
        })?;
        let write_txn = db.begin_write().map_err(redb_err)?;
        {
            let mut table = write_txn.open_table(DATA_TABLE).map_err(redb_err)?;
            table.insert(key, value).map_err(redb_err)?;
        }
        write_txn.commit().map_err(redb_err)?;
        Ok(())
    }

    pub(crate) fn batch_insert_raw(&self, entries: &[(Vec<u8>, Vec<u8>)]) -> Result<()> {
        if entries.is_empty() {
            return Ok(());
        }
        let _guard = self
            .write_lock
            .lock()
            .map_err(|e| DataFusionError::Execution(format!("redb write mutex poisoned: {e}")))?;
        let RedbHandle::Writable(db) = &self.db else {
            return Err(DataFusionError::Execution(
                "redb cache was opened read-only; writes are not allowed".to_string(),
            ));
        };
        let db = db.lock().map_err(|e| {
            DataFusionError::Execution(format!("redb database mutex poisoned: {e}"))
        })?;
        let write_txn = db.begin_write().map_err(redb_err)?;
        {
            let mut table = write_txn.open_table(DATA_TABLE).map_err(redb_err)?;
            for (key, value) in entries {
                table
                    .insert(key.as_slice(), value.as_slice())
                    .map_err(redb_err)?;
            }
        }
        write_txn.commit().map_err(redb_err)?;
        Ok(())
    }

    pub fn get_position_entry(&self, chrom_code: u16, start: i64) -> Result<Option<Vec<u8>>> {
        let mut key_buf = Vec::with_capacity(10);
        super::key_encoding::encode_position_key_buf(chrom_code, start, &mut key_buf);
        self.get_position_entry_raw(&key_buf)
    }

    pub(crate) fn get_position_entry_raw(&self, key: &[u8]) -> Result<Option<Vec<u8>>> {
        let read_txn = self.begin_read()?;
        let table = read_txn.open_table(DATA_TABLE).map_err(redb_err)?;
        table
            .get(key)
            .map_err(redb_err)
            .map(|value| value.map(|guard| guard.value().to_vec()))
    }

    pub fn get_position_entry_decompressed(
        &self,
        chrom_code: u16,
        start: i64,
    ) -> Result<Option<Vec<u8>>> {
        let raw = self.get_position_entry(chrom_code, start)?;
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
                    let mut decompressed = Vec::with_capacity(4096);
                    decompress_into_buffer_with_retry(
                        &mut decompressor,
                        &compressed,
                        &mut decompressed,
                        "zstd decompression failed",
                    )?;
                    Ok(Some(decompressed))
                } else {
                    Ok(Some(compressed))
                }
            }
        }
    }

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

    pub fn get_position_entry_fast(
        &self,
        chrom_code: u16,
        start: i64,
        decompressor: Option<&mut zstd::bulk::Decompressor<'_>>,
        buf: &mut Vec<u8>,
    ) -> Result<bool> {
        let mut key_buf = Vec::with_capacity(10);
        super::key_encoding::encode_position_key_buf(chrom_code, start, &mut key_buf);
        let read_txn = self.begin_read()?;
        let table = read_txn.open_table(DATA_TABLE).map_err(redb_err)?;
        let Some(raw) = table.get(key_buf.as_slice()).map_err(redb_err)? else {
            return Ok(false);
        };
        let compressed = raw.value();
        match decompressor {
            Some(dec) => {
                decompress_into_buffer_with_retry(
                    dec,
                    compressed,
                    buf,
                    "zstd decompression failed",
                )?;
            }
            None => {
                buf.clear();
                buf.extend_from_slice(compressed);
            }
        }
        Ok(true)
    }

    pub fn store_dict(&self, dict_bytes: &[u8]) -> Result<()> {
        self.put_metadata(ZSTD_DICT_KEY, dict_bytes)
    }

    pub fn store_zstd_level(&self, level: i32) -> Result<()> {
        self.put_metadata(ZSTD_LEVEL_KEY, &level.to_le_bytes())
    }

    pub fn store_contig_codes(&self, mapping: &[(String, u16)]) -> Result<()> {
        let bytes = super::key_encoding::serialize_contig_codes(mapping);
        self.put_metadata(CONTIG_CODES_KEY, &bytes)
    }

    pub fn put_metadata(&self, key: &str, value: &[u8]) -> Result<()> {
        let _guard = self
            .write_lock
            .lock()
            .map_err(|e| DataFusionError::Execution(format!("redb write mutex poisoned: {e}")))?;
        let RedbHandle::Writable(db) = &self.db else {
            return Err(DataFusionError::Execution(
                "redb cache was opened read-only; writes are not allowed".to_string(),
            ));
        };
        let db = db.lock().map_err(|e| {
            DataFusionError::Execution(format!("redb database mutex poisoned: {e}"))
        })?;
        let write_txn = db.begin_write().map_err(redb_err)?;
        {
            let mut meta = write_txn.open_table(META_TABLE).map_err(redb_err)?;
            meta.insert(key, value).map_err(redb_err)?;
        }
        write_txn.commit().map_err(redb_err)?;
        Ok(())
    }

    pub fn get_metadata(&self, key: &str) -> Result<Option<Vec<u8>>> {
        let read_txn = self.begin_read()?;
        let table = read_txn.open_table(META_TABLE).map_err(redb_err)?;
        table
            .get(key)
            .map_err(redb_err)
            .map(|value| value.map(|guard| guard.value().to_vec()))
    }

    pub fn persist(&self) -> Result<()> {
        Ok(())
    }

    pub(crate) fn optimize_after_load(&self) -> Result<()> {
        let RedbHandle::Writable(db) = &self.db else {
            return Ok(());
        };
        let _guard = self
            .write_lock
            .lock()
            .map_err(|e| DataFusionError::Execution(format!("redb write mutex poisoned: {e}")))?;
        let mut db = db.lock().map_err(|e| {
            DataFusionError::Execution(format!("redb database mutex poisoned: {e}"))
        })?;
        db.compact().map_err(redb_err)?;
        Ok(())
    }

    fn begin_read(&self) -> Result<ReadTransaction> {
        match &self.db {
            RedbHandle::Writable(db) => {
                let db = db.lock().map_err(|e| {
                    DataFusionError::Execution(format!("redb database mutex poisoned: {e}"))
                })?;
                db.begin_read().map_err(redb_err)
            }
            RedbHandle::ReadOnly(db) => db.begin_read().map_err(redb_err),
        }
    }

    fn read_metadata(db: &impl ReadableDatabase) -> Result<RedbMetadata> {
        let read_txn = db.begin_read().map_err(redb_err)?;
        let meta = read_txn.open_table(META_TABLE).map_err(redb_err)?;

        let schema_raw = meta.get(SCHEMA_KEY).map_err(redb_err)?.ok_or_else(|| {
            DataFusionError::Execution("redb store missing schema metadata".to_string())
        })?;
        let schema = schema_from_ipc_bytes(schema_raw.value())?;

        let version_raw = meta
            .get(FORMAT_VERSION_KEY)
            .map_err(redb_err)?
            .ok_or_else(|| {
                DataFusionError::Execution("redb cache format metadata missing".to_string())
            })?;
        let version = version_raw.value()[0];
        if version != FORMAT_V0 {
            return Err(DataFusionError::Execution(format!(
                "unsupported redb cache format version {version}: only version {FORMAT_V0} is supported"
            )));
        }

        let zstd_dict = meta
            .get(ZSTD_DICT_KEY)
            .map_err(redb_err)?
            .map(|raw| Arc::new(raw.value().to_vec()));

        if let Some(raw) = meta.get(CONTIG_CODES_KEY).map_err(redb_err)? {
            if let Some(mapping) = super::key_encoding::deserialize_contig_codes(raw.value()) {
                super::key_encoding::load_non_canonical_registry(&mapping);
            }
        }

        Ok((schema, version, zstd_dict))
    }
}

impl PositionEntryStore for VepRedbStore {
    fn schema(&self) -> SchemaRef {
        VepRedbStore::schema(self).clone()
    }

    fn root_path(&self) -> &Path {
        VepRedbStore::root_path(self)
    }

    fn create_decompressor(&self) -> Result<Option<zstd::bulk::Decompressor<'static>>> {
        VepRedbStore::create_decompressor(self)
    }

    fn get_position_entry_fast(
        &self,
        chrom_code: u16,
        start: i64,
        decompressor: Option<&mut zstd::bulk::Decompressor<'_>>,
        buf: &mut Vec<u8>,
    ) -> Result<bool> {
        VepRedbStore::get_position_entry_fast(self, chrom_code, start, decompressor, buf)
    }
}

/// TableProvider backed by a redb read-only VEP variation cache.
pub struct RedbCacheTableProvider {
    store: Arc<VepRedbStore>,
    schema: SchemaRef,
}

impl RedbCacheTableProvider {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let store = Arc::new(VepRedbStore::open(path)?);
        let schema = store.schema().clone();
        Ok(Self { store, schema })
    }

    pub fn from_store(store: Arc<VepRedbStore>) -> Self {
        let schema = store.schema().clone();
        Self { store, schema }
    }

    pub fn store(&self) -> &Arc<VepRedbStore> {
        &self.store
    }
}

impl Debug for RedbCacheTableProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "RedbCacheTableProvider {{ schema: {:?} }}", self.schema)
    }
}

#[async_trait]
impl TableProvider for RedbCacheTableProvider {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn table_type(&self) -> TableType {
        TableType::Base
    }

    async fn scan(
        &self,
        _state: &dyn Session,
        _projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        Err(DataFusionError::NotImplemented(
            "Direct scan of RedbCacheTableProvider is not yet supported. \
             Use lookup_variants() table function instead."
                .to_string(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

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
    fn redb_store_reopens_read_only_and_reads_position_entry() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let store = VepRedbStore::create(file.path(), test_schema()).unwrap();
        store.put_position_entry("1", 100, b"entry").unwrap();
        store.persist().unwrap();
        drop(store);

        let reopened = VepRedbStore::open(file.path()).unwrap();
        assert_eq!(reopened.format_version(), FORMAT_V0);
        let chrom_code = crate::kv_cache::key_encoding::chrom_to_code("1");
        let value = reopened
            .get_position_entry(chrom_code, 100)
            .unwrap()
            .unwrap();
        assert_eq!(value, b"entry");
    }

    #[test]
    fn redb_store_batch_insert_raw_writes_multiple_entries() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let store = VepRedbStore::create(file.path(), test_schema()).unwrap();
        let key_100 = crate::kv_cache::key_encoding::encode_position_key("1", 100);
        let key_200 = crate::kv_cache::key_encoding::encode_position_key("1", 200);

        store
            .batch_insert_raw(&[
                (key_100.clone(), b"entry-100".to_vec()),
                (key_200.clone(), b"entry-200".to_vec()),
            ])
            .unwrap();
        drop(store);

        let reopened = VepRedbStore::open(file.path()).unwrap();
        assert_eq!(
            reopened.get_position_entry_raw(&key_100).unwrap().unwrap(),
            b"entry-100"
        );
        assert_eq!(
            reopened.get_position_entry_raw(&key_200).unwrap().unwrap(),
            b"entry-200"
        );
    }
}
