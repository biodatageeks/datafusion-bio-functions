//! Test that all VCF input columns are properly passed through annotate_vep().
//!
//! Verifies that the annotation pipeline preserves every column from the input
//! VCF (core columns, INFO fields, FORMAT/sample fields) and adds the expected
//! annotation columns (csq, most_severe_consequence).

use std::sync::Arc;

use datafusion::datasource::TableProvider;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::register_vep_functions;

/// Resolve a path relative to the workspace root (two levels up from CARGO_MANIFEST_DIR).
fn workspace_path(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

#[tokio::test(flavor = "multi_thread")]
async fn test_vcf_columns_pass_through_annotation() {
    let input_vcf = workspace_path("vep-benchmark/data/golden/input_1000.vcf");
    if !input_vcf.exists() {
        eprintln!("Skipping: test fixtures not found at {}", input_vcf.display());
        return;
    }
    let input_vcf = input_vcf.to_str().unwrap();
    let cache_path = workspace_path("vep-benchmark/data/golden/cache");
    let cache_path = cache_path.to_str().unwrap();

    // Read input VCF schema.
    // None = include ALL INFO/FORMAT fields (Some(vec![]) means include NONE).
    // VcfTableProvider::new() uses futures::executor::block_on() internally,
    // which deadlocks inside a tokio runtime — use spawn_blocking.
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

    // Register and annotate.
    let ctx = SessionContext::new();
    register_vep_functions(&ctx);
    ctx.register_table("vcf", Arc::new(vcf_provider)).unwrap();

    let ref_fasta = workspace_path("vep-benchmark/data/golden/reference_chr1.fa");
    let ref_fasta = ref_fasta.to_str().unwrap();
    let sql = format!(
        "SELECT * FROM annotate_vep('vcf', '{cache_path}', 'parquet', \
         '{{\"partitioned\":true,\"everything\":true,\"extended_probes\":true,\
         \"reference_fasta_path\":\"{ref_fasta}\"}}')"
    );
    let df = ctx.sql(&sql).await.unwrap();
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
}
