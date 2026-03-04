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
    use datafusion::arrow::array::{Array, Float64Array, Int64Array, RecordBatch, StringArray};
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

    fn cache_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(StringArray::from(vec!["rs100"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["benign"])),
                Arc::new(Float64Array::from(vec![0.12_f64])),
            ],
        )
        .expect("valid cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid cache memtable")
    }

    fn transcripts_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["tx1", "tx2"])),
                Arc::new(StringArray::from(vec!["1", "2"])),
                Arc::new(Int64Array::from(vec![50, 150])),
                Arc::new(Int64Array::from(vec![200, 250])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(StringArray::from(vec!["protein_coding", "lincRNA"])),
                Arc::new(Int64Array::from(vec![80, 0])),
                Arc::new(Int64Array::from(vec![180, 0])),
            ],
        )
        .expect("valid transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid transcript memtable")
    }

    fn exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["tx1", "tx2"])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(Int64Array::from(vec![50, 150])),
                Arc::new(Int64Array::from(vec![200, 250])),
            ],
        )
        .expect("valid exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid exon memtable")
    }

    fn string_values(col: &Arc<dyn datafusion::arrow::array::Array>) -> Vec<Option<String>> {
        if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else if let Some(arr) = col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringViewArray>()
        {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else {
            panic!("expected string or string view");
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_appends_annotation_columns() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");

        let df = ctx
            .sql("SELECT * FROM annotate_vep('vcf_data', 'var_cache', 'parquet')")
            .await
            .expect("annotate_vep query should parse");

        let batches = df.collect().await.expect("collect annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);

        let mut csq_values = Vec::new();
        let mut most_values = Vec::new();
        for batch in &batches {
            assert!(batch.column_by_name("csq").is_some());
            assert!(batch.column_by_name("most_severe_consequence").is_some());
            csq_values.extend(string_values(
                batch.column_by_name("csq").expect("csq column exists"),
            ));
            most_values.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }
        assert!(
            csq_values
                .iter()
                .any(|v| v.as_ref().is_some_and(|s| s.contains("sequence_variant")))
        );
        assert!(csq_values.iter().any(|v| v.is_none()));
        assert!(
            most_values
                .iter()
                .any(|v| v.as_ref() == Some(&"sequence_variant".to_string()))
        );
        assert!(most_values.iter().any(|v| v.is_none()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_projection_includes_null_placeholder_fields() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");

        let df = ctx
            .sql("SELECT chrom, csq FROM annotate_vep('vcf_data', 'var_cache', 'parquet')")
            .await
            .expect("projection query should parse");

        let batches = df.collect().await.expect("collect projected annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);
        for batch in &batches {
            assert_eq!(batch.num_columns(), 2);
            assert_eq!(batch.schema().field(0).name(), "chrom");
            assert_eq!(batch.schema().field(1).name(), "csq");
        }
        let mut csq_values = Vec::new();
        for batch in &batches {
            csq_values.extend(string_values(
                batch.column_by_name("csq").expect("csq column exists"),
            ));
        }
        assert!(csq_values.iter().any(|v| v.is_some()));
        assert!(csq_values.iter().any(|v| v.is_none()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_transcript_context_tables_when_available() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");
        ctx.register_table("var_cache_transcripts", Arc::new(transcripts_table()))
            .expect("register transcripts table");
        ctx.register_table("var_cache_exons", Arc::new(exons_table()))
            .expect("register exons table");

        let df = ctx
            .sql(
                "SELECT chrom, csq, most_severe_consequence \
                 FROM annotate_vep('vcf_data', 'var_cache', 'parquet') \
                 ORDER BY chrom",
            )
            .await
            .expect("query should parse");

        let batches = df
            .collect()
            .await
            .expect("collect transcript-aware annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);

        let mut chrom = Vec::new();
        let mut csq = Vec::new();
        let mut most = Vec::new();
        for batch in &batches {
            chrom.extend(string_values(
                batch.column_by_name("chrom").expect("chrom column exists"),
            ));
            csq.extend(string_values(
                batch.column_by_name("csq").expect("csq column exists"),
            ));
            most.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }

        assert_eq!(chrom, vec![Some("1".to_string()), Some("2".to_string())]);
        assert!(csq.iter().all(|v| v.is_some()));
        assert!(most.iter().all(|v| v.is_some()));
        assert!(
            csq[0]
                .as_ref()
                .is_some_and(|s| s.contains("missense_variant"))
        );
        assert!(csq[0].as_ref().is_some_and(|s| s.contains("rs100")));
        assert!(
            csq[1]
                .as_ref()
                .is_some_and(|s| s.contains("non_coding_transcript_exon_variant"))
        );
        assert_eq!(
            most,
            vec![
                Some("missense_variant".to_string()),
                Some("non_coding_transcript_exon_variant".to_string())
            ]
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_options_json_table_overrides_for_transcript_context() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");
        ctx.register_table("tx_ctx", Arc::new(transcripts_table()))
            .expect("register transcript context table");
        ctx.register_table("ex_ctx", Arc::new(exons_table()))
            .expect("register exon context table");

        let df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep( \
                   'vcf_data', \
                   'var_cache', \
                   'parquet', \
                   '{\"transcripts_table\":\"tx_ctx\",\"exons_table\":\"ex_ctx\"}' \
                 )",
            )
            .await
            .expect("query should parse");

        let batches = df
            .collect()
            .await
            .expect("collect transcript-aware annotate_vep");
        let mut csq = Vec::new();
        let mut most = Vec::new();
        for batch in &batches {
            csq.extend(string_values(
                batch.column_by_name("csq").expect("csq column exists"),
            ));
            most.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }
        assert!(csq.iter().all(|v| v.is_some()));
        assert!(most.iter().all(|v| v.is_some()));
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
