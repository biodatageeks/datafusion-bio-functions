//! Table function registration for VEP lookup.

use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::TableFunctionImpl;
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;

use crate::lookup_provider::LookupProvider;
use crate::schema_contract::{
    DEFAULT_LOOKUP_COLUMNS, parse_column_list, validate_requested_columns,
};

/// Table function implementing `lookup_variants(vcf_table, cache_table [, columns [, prune]])`.
pub struct LookupFunction {
    session: Arc<SessionContext>,
}

impl LookupFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
        Self { session }
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

        // Optional third argument: comma-separated column list
        let columns = if args.len() > 2 {
            let col_str = extract_string_arg(&args[2], "columns", "lookup_variants")?;
            parse_column_list(&col_str)
        } else {
            DEFAULT_LOOKUP_COLUMNS
                .iter()
                .map(|s| s.to_string())
                .collect()
        };

        // Optional fourth argument: prune all-null columns (default: false)
        let prune_nulls = if args.len() > 3 {
            extract_bool_arg(&args[3], "prune_nulls", "lookup_variants")?
        } else {
            false
        };

        // Resolve schemas from registered tables
        let (vcf_schema, cache_schema) = resolve_schemas(&self.session, &vcf_table, &cache_table)?;

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
            prune_nulls,
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

/// Resolve schemas for both tables, handling tokio context.
fn resolve_schemas(
    session: &SessionContext,
    vcf_table: &str,
    cache_table: &str,
) -> Result<(Schema, Schema)> {
    match tokio::runtime::Handle::try_current() {
        Ok(handle) => tokio::task::block_in_place(|| {
            let vcf = handle.block_on(session.table(vcf_table))?;
            let cache = handle.block_on(session.table(cache_table))?;
            Ok::<_, DataFusionError>((
                vcf.schema().as_arrow().clone(),
                cache.schema().as_arrow().clone(),
            ))
        }),
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            let vcf = rt.block_on(session.table(vcf_table))?;
            let cache = rt.block_on(session.table(cache_table))?;
            Ok((
                vcf.schema().as_arrow().clone(),
                cache.schema().as_arrow().clone(),
            ))
        }
    }
}
