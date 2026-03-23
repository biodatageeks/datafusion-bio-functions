//! Table function registration for VEP lookup.
//!
//! # Join semantics
//!
//! By default, `lookup_variants()` uses a self-contained left interval join
//! (`VariantLookupExec`) with exact coordinate matching after VEP normalization
//! and `match_allele()` as a post-filter.
//!
//! When `extended_probes = true`, the join uses interval-overlap matching
//! instead of exact coordinate equality, handling insertion-style coordinates
//! and shifted-deletion cache coordinates.

use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::{CatalogProviderList, SchemaProvider, TableFunctionImpl};
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;

use crate::lookup_provider::LookupProvider;
use crate::schema_contract::{COORDINATE_COLUMNS, parse_column_list, validate_requested_columns};

/// Table function implementing
/// `lookup_variants(vcf_table, cache_table [, columns [, match_mode [, extended_probes]]])`.
pub struct LookupFunction {
    session: Arc<SessionContext>,
    /// Catalog list captured at registration time to avoid acquiring
    /// SessionState locks during `call()` (planning time).
    catalog_list: Arc<dyn CatalogProviderList>,
    default_catalog: String,
    default_schema: String,
}

impl LookupFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
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

impl std::fmt::Debug for LookupFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LookupFunction")
    }
}

impl TableFunctionImpl for LookupFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.len() < 2 {
            return Err(DataFusionError::Plan(
                "lookup_variants() requires at least 2 arguments: \
                 vcf_table and variation_cache_table"
                    .to_string(),
            ));
        }

        let vcf_table = extract_string_arg(&args[0], "vcf_table", "lookup_variants")?;
        let cache_table = extract_string_arg(&args[1], "cache_table", "lookup_variants")?;

        // Resolve schemas via pre-captured catalog list (no session state lock).
        let (vcf_schema, cache_schema) = resolve_schemas_from_catalog(
            &*self.catalog_list,
            &self.default_catalog,
            &self.default_schema,
            &vcf_table,
            &cache_table,
        )?;

        // Optional third argument: comma-separated column list.
        // Default: all cache columns except coordinate columns (chrom, start, end)
        // which are already present on the VCF side.
        let columns = if args.len() > 2 {
            let col_str = extract_string_arg(&args[2], "columns", "lookup_variants")?;
            parse_column_list(&col_str)
        } else {
            cache_schema
                .fields()
                .iter()
                .map(|f| f.name().clone())
                .filter(|name| {
                    !COORDINATE_COLUMNS.contains(&name.as_str()) && !name.starts_with("source_")
                })
                .collect()
        };

        // Optional fourth argument: match mode (default: "exact").
        // Only "exact" is supported. The parameter is accepted for backwards
        // compatibility but other values are rejected.
        if args.len() > 3 {
            let mode = extract_string_arg(&args[3], "match_mode", "lookup_variants")?;
            if mode != "exact" {
                return Err(DataFusionError::Plan(format!(
                    "lookup_variants() match_mode must be 'exact'; got: {mode}"
                )));
            }
        }

        // Optional fifth argument: extended_probes (default: false)
        // When true, uses interval-overlap matching and multi-probe KV lookups.
        let extended_probes = if args.len() > 4 {
            extract_bool_arg(&args[4], "extended_probes", "lookup_variants")?
        } else {
            false
        };

        // Validate requested columns exist in cache schema
        let cache_schema_ref = Arc::new(cache_schema.clone());
        validate_requested_columns(&cache_schema_ref, &columns)?;

        Ok(Arc::new(LookupProvider::new(
            Arc::clone(&self.session),
            vcf_table,
            cache_table,
            vcf_schema,
            cache_schema,
            columns,
            extended_probes,
            0,
            None,
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

/// Extract a boolean literal from an expression.
fn extract_bool_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<bool> {
    match arg {
        Expr::Literal(ScalarValue::Boolean(Some(val)), _) => Ok(*val),
        other => Err(DataFusionError::Plan(format!(
            "{fn_name}() {name} must be a boolean literal, got: {other}"
        ))),
    }
}

/// Resolve schemas for both tables using a pre-captured `CatalogProviderList`,
/// completely bypassing `SessionContext` and its `SessionState` RwLock.
///
/// See `annotate_table_function::resolve_schema_from_catalog` for rationale.
fn resolve_schemas_from_catalog(
    catalog_list: &dyn CatalogProviderList,
    default_catalog: &str,
    default_schema: &str,
    vcf_table: &str,
    cache_table: &str,
) -> Result<(Schema, Schema)> {
    let vcf_provider = resolve_table_ref(catalog_list, default_catalog, default_schema, vcf_table)?;
    let cache_provider =
        resolve_table_ref(catalog_list, default_catalog, default_schema, cache_table)?;

    Ok((
        vcf_provider.schema().as_ref().clone(),
        cache_provider.schema().as_ref().clone(),
    ))
}

/// Resolve a table by name, supporting bare, schema-qualified, and fully-qualified references.
fn resolve_table_ref(
    catalog_list: &dyn CatalogProviderList,
    default_catalog: &str,
    default_schema: &str,
    table_name: &str,
) -> Result<Arc<dyn TableProvider>> {
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
    resolve_table_sync(&*schema_provider, bare_name)
}

/// Run `SchemaProvider::table()` synchronously, handling both tokio-context
/// and no-tokio-context cases.
fn resolve_table_sync(
    schema_provider: &dyn SchemaProvider,
    table_name: &str,
) -> Result<Arc<dyn TableProvider>> {
    let result = match tokio::runtime::Handle::try_current() {
        Ok(handle) => {
            tokio::task::block_in_place(|| handle.block_on(schema_provider.table(table_name)))
        }
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
