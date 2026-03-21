//! KvCacheTableProvider: TableProvider for the fjall-backed VEP cache.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::path::Path;
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::catalog::Session;
use datafusion::common::Result;
use datafusion::datasource::{TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::Expr;

use super::kv_store::VepKvStore;

/// TableProvider backed by a fjall KV store containing VEP cache data.
///
/// Register this as a table in a DataFusion session to use it with
/// `lookup_variants()` or directly in SQL queries.
pub struct KvCacheTableProvider {
    store: Arc<VepKvStore>,
    schema: SchemaRef,
}

impl KvCacheTableProvider {
    /// Open an existing KV cache at the given path (default 256MB block cache).
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_cache_size(path, 256 * 1024 * 1024)
    }

    /// Open an existing KV cache with a custom fjall block cache size (bytes).
    pub fn open_with_cache_size(path: impl AsRef<Path>, cache_size_bytes: u64) -> Result<Self> {
        let store = Arc::new(VepKvStore::open_with_cache_size(path, cache_size_bytes)?);
        let schema = store.schema().clone();
        Ok(Self { store, schema })
    }

    /// Wrap an already-opened `VepKvStore` without re-opening the database.
    pub fn from_store(store: Arc<VepKvStore>) -> Self {
        let schema = store.schema().clone();
        Self { store, schema }
    }

    /// Get a reference to the underlying KV store.
    pub fn store(&self) -> &Arc<VepKvStore> {
        &self.store
    }
}

impl Debug for KvCacheTableProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "KvCacheTableProvider {{ schema: {:?} }}", self.schema)
    }
}

#[async_trait]
impl TableProvider for KvCacheTableProvider {
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
        Err(datafusion::common::DataFusionError::NotImplemented(
            "Direct scan of KvCacheTableProvider is not yet supported. \
             Use lookup_variants() table function instead."
                .to_string(),
        ))
    }
}
