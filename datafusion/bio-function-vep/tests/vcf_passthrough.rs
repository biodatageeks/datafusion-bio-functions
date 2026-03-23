//! Test that all VCF input columns are properly passed through annotate_vep().
//!
//! Verifies that the annotation pipeline preserves every column from the input
//! VCF (core columns, INFO fields, FORMAT/sample fields) and adds the expected
//! annotation columns (csq, most_severe_consequence).

use std::sync::Arc;

use datafusion::arrow::array::{Array, StringArray, StringViewArray};
use datafusion::datasource::TableProvider;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::register_vep_functions;

/// Resolve a path relative to the workspace root (two levels up from CARGO_MANIFEST_DIR).
/// Check if a file is a Git LFS pointer (not actual content).
fn is_lfs_pointer(path: &std::path::Path) -> bool {
    std::fs::read_to_string(path)
        .map(|s| s.starts_with("version https://git-lfs.github.com"))
        .unwrap_or(false)
}

fn workspace_path(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

fn count_non_null_strings(array: &dyn Array) -> usize {
    if let Some(arr) = array.as_any().downcast_ref::<StringArray>() {
        arr.len() - arr.null_count()
    } else if let Some(arr) = array.as_any().downcast_ref::<StringViewArray>() {
        arr.len() - arr.null_count()
    } else {
        panic!("expected a string array");
    }
}

async fn build_test_query(
    sql_select: &str,
) -> Option<(datafusion::prelude::DataFrame, Vec<String>)> {
    let input_vcf = workspace_path("vep-benchmark/data/golden/input_1000.vcf");
    if !input_vcf.exists() || is_lfs_pointer(&input_vcf) {
        eprintln!(
            "Skipping: test fixtures not available at {}",
            input_vcf.display()
        );
        return None;
    }
    let input_vcf = input_vcf.to_str().unwrap();
    let cache_path = workspace_path("vep-benchmark/data/golden/cache");
    let cache_path = cache_path.to_str().unwrap();

    let input_str = input_vcf.to_string();
    let vcf_provider = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(input_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    let input_schema = vcf_provider.schema();
    let input_field_names: Vec<String> = input_schema
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();

    let ctx = SessionContext::new();
    register_vep_functions(&ctx);
    ctx.register_table("vcf", Arc::new(vcf_provider)).unwrap();

    let ref_fasta = workspace_path("vep-benchmark/data/golden/reference_chr1.fa");
    let ref_fasta = ref_fasta.to_str().unwrap();
    let sql = format!(
        "{sql_select} FROM annotate_vep('vcf', '{cache_path}', 'parquet', \
         '{{\"partitioned\":true,\"everything\":true,\"extended_probes\":true,\
         \"reference_fasta_path\":\"{ref_fasta}\"}}')"
    );
    let df = ctx.sql(&sql).await.unwrap();
    Some((df, input_field_names))
}

#[tokio::test(flavor = "multi_thread")]
async fn test_vcf_columns_pass_through_annotation() {
    let Some((df, input_field_names)) = build_test_query("SELECT *").await else {
        return;
    };
    let output_schema = df.schema().clone();

    // Verify ALL input columns exist in the output.
    for name in &input_field_names {
        assert!(
            output_schema.field_with_name(None, name).is_ok(),
            "Input column '{name}' missing from annotate_vep output"
        );
    }

    // Verify annotation columns are added.
    assert!(
        output_schema.field_with_name(None, "csq").is_ok(),
        "Expected 'csq' column in annotate_vep output"
    );
    assert!(
        output_schema
            .field_with_name(None, "most_severe_consequence")
            .is_ok(),
        "Expected 'most_severe_consequence' column in annotate_vep output"
    );

    // Verify total column count: input + 89 annotation columns.
    assert_eq!(
        output_schema.fields().len(),
        input_field_names.len() + 89,
        "Output should have input columns + 89 annotation columns"
    );

    // Collect and verify row count.
    let batches = df.collect().await.unwrap();
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 1000, "Expected 1000 annotated rows");

    // Regression check: on an unprojected SELECT *, csq must be populated.
    let non_null_csq: usize = batches
        .iter()
        .map(|batch| {
            let csq_col = batch.column(batch.schema().index_of("csq").unwrap());
            count_non_null_strings(csq_col.as_ref())
        })
        .sum();
    assert!(
        non_null_csq > 0,
        "Expected non-null csq values for SELECT * annotate_vep()"
    );
}

#[tokio::test(flavor = "multi_thread")]
async fn test_projection_including_csq_preserves_values() {
    let Some((df, _)) = build_test_query("SELECT chrom, start, csq").await else {
        return;
    };

    let schema = df.schema().clone();
    assert!(schema.field_with_name(None, "csq").is_ok());

    let batches = df.collect().await.unwrap();
    let non_null_csq: usize = batches
        .iter()
        .map(|batch| {
            let csq_col = batch.column(batch.schema().index_of("csq").unwrap());
            count_non_null_strings(csq_col.as_ref())
        })
        .sum();
    assert!(
        non_null_csq > 0,
        "Expected non-null csq values when csq is explicitly projected"
    );
}

#[tokio::test(flavor = "multi_thread")]
async fn test_projection_omitting_csq_skips_column() {
    let Some((df, _)) = build_test_query("SELECT chrom, start").await else {
        return;
    };

    let schema = df.schema().clone();
    assert!(schema.field_with_name(None, "chrom").is_ok());
    assert!(schema.field_with_name(None, "start").is_ok());
    assert!(schema.field_with_name(None, "csq").is_err());

    let batches = df.collect().await.unwrap();
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 1000, "Expected 1000 projected rows");
}
