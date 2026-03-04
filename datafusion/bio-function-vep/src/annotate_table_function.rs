//! Table function registration for consequence annotation.
//!
//! `annotate_vep()` is the high-level consequence annotation entrypoint.

use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::TableFunctionImpl;
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;

use crate::annotate_provider::AnnotateProvider;
use crate::annotation_store::AnnotationBackend;

/// Table function implementing
/// `annotate_vep(vcf_table, cache_source, backend [, options_json])`.
pub struct AnnotateFunction {
    session: Arc<SessionContext>,
}

impl AnnotateFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
        Self { session }
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

        let (vcf_schema, _) = resolve_schema(&self.session, &vcf_table)?;

        Ok(Arc::new(AnnotateProvider::new(
            Arc::clone(&self.session),
            vcf_table,
            cache_source,
            backend,
            options_json,
            vcf_schema,
        )))
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

fn resolve_schema(session: &SessionContext, vcf_table: &str) -> Result<(Schema, String)> {
    match tokio::runtime::Handle::try_current() {
        Ok(handle) => tokio::task::block_in_place(|| {
            let vcf = handle.block_on(session.table(vcf_table))?;
            Ok::<_, DataFusionError>((vcf.schema().as_arrow().clone(), vcf_table.to_string()))
        }),
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            let vcf = rt.block_on(session.table(vcf_table))?;
            Ok((vcf.schema().as_arrow().clone(), vcf_table.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::create_vep_session;
    use datafusion::arrow::array::{Int64Array, RecordBatch, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use std::sync::Arc;

    fn vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "2"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![101, 201])),
                Arc::new(StringArray::from(vec!["A", "C"])),
                Arc::new(StringArray::from(vec!["G", "T"])),
            ],
        )
        .expect("valid vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid vcf memtable")
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_appends_annotation_columns() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let df = ctx
            .sql("SELECT * FROM annotate_vep('vcf_data', '/tmp/vep_cache', 'parquet')")
            .await
            .expect("annotate_vep query should parse");

        let batches = df.collect().await.expect("collect annotate_vep");
        let batch = &batches[0];

        assert!(batch.column_by_name("csq").is_some());
        assert!(batch.column_by_name("most_severe_consequence").is_some());
        assert_eq!(batch.num_rows(), 2);

        let csq = batch.column_by_name("csq").expect("csq column exists");
        let most_severe = batch
            .column_by_name("most_severe_consequence")
            .expect("most_severe_consequence column exists");

        assert_eq!(csq.null_count(), batch.num_rows());
        assert_eq!(most_severe.null_count(), batch.num_rows());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_projection_includes_null_placeholder_fields() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let df = ctx
            .sql("SELECT chrom, csq FROM annotate_vep('vcf_data', '/tmp/vep_cache', 'fjall')")
            .await
            .expect("projection query should parse");

        let batches = df.collect().await.expect("collect projected annotate_vep");
        let batch = &batches[0];

        assert_eq!(batch.num_columns(), 2);
        assert_eq!(batch.schema().field(0).name(), "chrom");
        assert_eq!(batch.schema().field(1).name(), "csq");
        assert_eq!(
            batch
                .column_by_name("csq")
                .expect("csq column exists")
                .null_count(),
            batch.num_rows()
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_rejects_unknown_backend() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let err = ctx
            .sql("SELECT * FROM annotate_vep('vcf_data', '/tmp/vep_cache', 'bad_backend')")
            .await
            .expect_err("unknown backend should fail")
            .to_string();

        assert!(err.contains("annotate_vep() backend must be one of"));
    }
}
