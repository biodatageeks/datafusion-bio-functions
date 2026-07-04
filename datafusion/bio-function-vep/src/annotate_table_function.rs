//! Table function registration for consequence annotation.
//!
//! `annotate_vep()` is the high-level consequence annotation entrypoint.

use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::{CatalogProviderList, SchemaProvider, TableFunctionImpl};
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;
use serde_json::Value;

use crate::annotate_provider::AnnotateProvider;
use crate::annotation_store::AnnotationBackend;
use crate::cache_source::{CACHE_SOURCE_METADATA_KEY, CacheSourceType};

/// Table function implementing
/// `annotate_vep(vcf_table, cache_source, backend [, options_json])`.
///
/// Cache source mode is read from Arrow schema metadata on the Parquet cache
/// backend under `{cache_source}/variation.cache`.
pub struct AnnotateFunction {
    session: Arc<SessionContext>,
    /// Catalog list captured at registration time to avoid acquiring
    /// SessionState locks during `call()` (planning time).
    catalog_list: Arc<dyn CatalogProviderList>,
    /// Default catalog name captured at registration time.
    default_catalog: String,
    /// Default schema name captured at registration time.
    default_schema: String,
}

impl AnnotateFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
        // Capture catalog metadata at registration time (before any planner
        // locks exist) so resolve_schema can look up tables without touching
        // the SessionState RwLock.
        let state = session.state();
        let catalog_list = Arc::clone(state.catalog_list());
        let default_catalog = state.config_options().catalog.default_catalog.clone();
        let default_schema = state.config_options().catalog.default_schema.clone();
        Self {
            session,
            catalog_list,
            default_catalog,
            default_schema,
        }
    }
}

impl std::fmt::Debug for AnnotateFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "AnnotateFunction")
    }
}

impl TableFunctionImpl for AnnotateFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.len() < 3 {
            return Err(DataFusionError::Plan(
                "annotate_vep() requires at least 3 arguments: vcf_table, cache_source, backend"
                    .to_string(),
            ));
        }

        let vcf_table = extract_string_arg(&args[0], "vcf_table", "annotate_vep")?;
        let cache_source = extract_string_arg(&args[1], "cache_source", "annotate_vep")?;
        let backend_raw = extract_string_arg(&args[2], "backend", "annotate_vep")?;
        let backend = AnnotationBackend::parse(&backend_raw)?;

        let options_json = if args.len() > 3 {
            Some(extract_string_arg(
                &args[3],
                "options_json",
                "annotate_vep",
            )?)
        } else {
            None
        };
        reject_options_json_source_selectors(options_json.as_deref())?;
        // The cache is always the partitioned Parquet cache; read its source-type
        // metadata off the first variation shard.
        #[cfg(feature = "parquet-cache")]
        let cache_source_type =
            CacheSourceType::from_partitioned_parquet_cache_source(&cache_source)?;
        #[cfg(not(feature = "parquet-cache"))]
        let cache_source_type: CacheSourceType = {
            return Err(DataFusionError::Plan(
                "annotate_vep(): Parquet cache source metadata requires the parquet-cache feature"
                    .to_string(),
            ));
        };

        let vcf_schema = resolve_schema_from_catalog(
            &*self.catalog_list,
            &self.default_catalog,
            &self.default_schema,
            &vcf_table,
        )?;

        Ok(Arc::new(AnnotateProvider::new(
            Arc::clone(&self.session),
            vcf_table,
            cache_source,
            backend,
            cache_source_type,
            options_json,
            vcf_schema,
        )?))
    }
}

/// Extract a string literal from an expression.
fn extract_string_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<String> {
    match arg {
        Expr::Literal(ScalarValue::Utf8(Some(val)), _) => {
            if val.contains('`') {
                return Err(DataFusionError::Plan(format!(
                    "{fn_name}() {name} must not contain backtick characters, got: {val}"
                )));
            }
            Ok(val.clone())
        }
        other => Err(DataFusionError::Plan(format!(
            "{fn_name}() {name} must be a string literal, got: {other}"
        ))),
    }
}

/// Resolve the Arrow schema of a registered table using a pre-captured
/// `CatalogProviderList`, completely bypassing `SessionContext` and its
/// `SessionState` RwLock.
///
/// This is critical for vepyr / polars-bio integration: DataFusion's SQL
/// planner holds a **write lock** on `SessionState` while resolving table
/// functions. Any call to `session.state()`, `session.catalog()`, or
/// `session.table()` from within `TableFunctionImpl::call()` will deadlock
/// because they all acquire a read lock on the same RwLock.
///
/// By capturing the `CatalogProviderList` at registration time (before any
/// planner locks exist), we can look up tables without touching the session.
fn resolve_schema_from_catalog(
    catalog_list: &dyn CatalogProviderList,
    default_catalog: &str,
    default_schema: &str,
    table_name: &str,
) -> Result<Schema> {
    // Support bare names ("vcf"), schema-qualified ("public.vcf"), and
    // fully-qualified ("datafusion.public.vcf") table references.
    let parts: Vec<&str> = table_name.split('.').collect();
    let (cat_name, schema_name, bare_name) = match parts.len() {
        3 => (parts[0], parts[1], parts[2]),
        2 => (default_catalog, parts[0], parts[1]),
        _ => (default_catalog, default_schema, table_name),
    };

    let catalog = catalog_list
        .catalog(cat_name)
        .ok_or_else(|| DataFusionError::Plan(format!("Catalog '{cat_name}' not found")))?;
    let schema_provider = catalog.schema(schema_name).ok_or_else(|| {
        DataFusionError::Plan(format!(
            "Schema '{schema_name}' not found in catalog '{cat_name}'"
        ))
    })?;

    let table_provider = resolve_table_sync(&*schema_provider, bare_name)?;
    Ok(table_provider.schema().as_ref().clone())
}

/// Run `SchemaProvider::table()` synchronously, handling both tokio-context
/// and no-tokio-context cases.
fn resolve_table_sync(
    schema_provider: &dyn SchemaProvider,
    table_name: &str,
) -> Result<Arc<dyn TableProvider>> {
    let result = match tokio::runtime::Handle::try_current() {
        // Multi-thread runtime: `block_in_place` is supported and cheapest.
        Ok(handle) if handle.runtime_flavor() == tokio::runtime::RuntimeFlavor::MultiThread => {
            tokio::task::block_in_place(|| handle.block_on(schema_provider.table(table_name)))
        }
        // Current-thread runtime (e.g. `#[tokio::test]`, embedded callers):
        // `block_in_place` would panic and we cannot nest a runtime on this
        // thread, so resolve on a dedicated OS thread with its own runtime.
        Ok(_handle) => std::thread::scope(|scope| {
            scope
                .spawn(|| {
                    let rt = tokio::runtime::Runtime::new()
                        .map_err(|e| DataFusionError::External(Box::new(e)))?;
                    rt.block_on(schema_provider.table(table_name))
                })
                .join()
                .unwrap_or_else(|_| {
                    Err(DataFusionError::Execution(
                        "table resolution worker thread panicked".to_string(),
                    ))
                })
        }),
        // No runtime in scope: create one here.
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            rt.block_on(schema_provider.table(table_name))
        }
    };
    result
        .map_err(|e| DataFusionError::External(Box::new(e)))?
        .ok_or_else(|| DataFusionError::Plan(format!("Table '{table_name}' not found")))
}

fn options_json_has_key(options_json: Option<&str>, key: &str) -> Result<bool> {
    let Some(raw) = options_json else {
        return Ok(false);
    };
    let value: Value = serde_json::from_str(raw).map_err(|err| {
        DataFusionError::Plan(format!(
            "annotate_vep() options_json must be valid JSON: {err}"
        ))
    })?;
    Ok(value.as_object().is_some_and(|obj| obj.contains_key(key)))
}

fn options_json_string_value(options_json: Option<&str>, key: &str) -> Result<Option<String>> {
    let Some(raw) = options_json else {
        return Ok(None);
    };
    let value: Value = serde_json::from_str(raw).map_err(|err| {
        DataFusionError::Plan(format!(
            "annotate_vep() options_json must be valid JSON: {err}"
        ))
    })?;
    let Some(raw_value) = value.as_object().and_then(|obj| obj.get(key)) else {
        return Ok(None);
    };
    let Some(text) = raw_value.as_str() else {
        return Err(DataFusionError::Plan(format!(
            "annotate_vep(): options_json key '{key}' must be a string"
        )));
    };
    if text.contains('`') {
        return Err(DataFusionError::Plan(format!(
            "annotate_vep(): options_json key '{key}' must not contain backtick characters"
        )));
    }
    Ok(Some(text.to_string()))
}

fn reject_options_json_source_selectors(options_json: Option<&str>) -> Result<()> {
    for key in ["merged", "refseq"] {
        if options_json_has_key(options_json, key)? {
            return Err(DataFusionError::Plan(format!(
                "annotate_vep(): options_json key '{key}' is unsupported; register/export cache tables with schema metadata {CACHE_SOURCE_METADATA_KEY}='{key}'"
            )));
        }
    }
    Ok(())
}
