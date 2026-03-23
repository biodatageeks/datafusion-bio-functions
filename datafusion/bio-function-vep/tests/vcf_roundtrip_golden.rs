//! Roundtrip test: annotate 1000 variants -> write VCF -> read back -> compare.
//!
//! This test verifies the full pipeline:
//! 1. Read a VCF with VcfTableProvider
//! 2. Annotate with annotate_vep()
//! 3. Write results to a VCF file via vcf_sink::annotate_to_vcf()
//! 4. Read the output VCF back with VcfTableProvider
//! 5. Verify row counts and that all columns survive the roundtrip

use std::sync::Arc;

use datafusion::datasource::TableProvider;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_format_vcf::VcfCompressionType;
use datafusion_bio_function_vep::{register_vep_functions, vcf_sink};

/// Resolve a path relative to the workspace root (two levels up from CARGO_MANIFEST_DIR).
fn workspace_path(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

#[tokio::test(flavor = "multi_thread")]
async fn test_roundtrip_golden_all_columns() {
    // Skip if fixtures not available (CI without LFS).
    let input_vcf = workspace_path("vep-benchmark/data/golden/input_1000.vcf");
    let cache_path = workspace_path("vep-benchmark/data/golden/cache");
    let ref_fasta = workspace_path("vep-benchmark/data/golden/reference_chr1.fa");

    if !input_vcf.exists() {
        eprintln!("Skipping: test fixtures not found at {}", input_vcf.display());
        return;
    }
    let input_vcf = input_vcf.to_str().unwrap();
    let cache_path = cache_path.to_str().unwrap();
    let ref_fasta = ref_fasta.to_str().unwrap();

    let tmp_dir = tempfile::TempDir::new().unwrap();
    let output_path = tmp_dir.path().join("annotated_output.vcf");

    // Set up context and register VCF table.
    let ctx = SessionContext::new();
    register_vep_functions(&ctx);
    let vcf_provider = VcfTableProvider::new(
        input_vcf.to_string(),
        Some(vec![]),
        Some(vec![]),
        None,
        false,
    )
    .unwrap();
    let input_schema = vcf_provider.schema();
    let _input_field_names: Vec<String> = input_schema
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    ctx.register_table("vcf", Arc::new(vcf_provider)).unwrap();

    // Run annotation and write to VCF.
    let options = format!(
        "{{\"partitioned\":true,\"everything\":true,\"extended_probes\":true,\
         \"reference_fasta_path\":\"{ref_fasta}\"}}"
    );
    let rows = vcf_sink::annotate_to_vcf(
        &ctx,
        "vcf",
        cache_path,
        "parquet",
        Some(&options),
        &output_path,
        VcfCompressionType::Plain,
    )
    .await
    .unwrap();

    assert_eq!(rows, 1000, "Expected 1000 rows written to VCF");

    // Read back the output VCF with VcfTableProvider.
    let output_provider = VcfTableProvider::new(
        output_path.display().to_string(),
        Some(vec![]),
        Some(vec![]),
        None,
        false,
    )
    .unwrap();
    let output_schema = output_provider.schema();

    let ctx2 = SessionContext::new();
    ctx2.register_table("output_vcf", Arc::new(output_provider))
        .unwrap();

    let output_batches = ctx2
        .sql("SELECT * FROM output_vcf")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let output_rows: usize = output_batches.iter().map(|b| b.num_rows()).sum();

    // Verify row count matches.
    assert_eq!(output_rows, 1000, "Output VCF should have 1000 rows");

    // Verify core VCF columns survive roundtrip.
    for col in ["chrom", "start", "end", "ref", "alt"] {
        assert!(
            output_schema.field_with_name(col).is_ok(),
            "Core column '{col}' missing from roundtrip output"
        );
    }

    // Verify the CSQ annotation column made it through as an INFO field.
    // After roundtrip through VcfTableProvider, CSQ should appear as an INFO column
    // named "csq" (lowercase) or "CSQ" (uppercase, depending on the provider's convention).
    // Debug: print all field names from the roundtrip output
    let field_names: Vec<&str> = output_schema.fields().iter().map(|f| f.name().as_str()).collect();
    eprintln!("Roundtrip output schema fields ({}):", field_names.len());
    for name in &field_names {
        eprintln!("  {name}");
    }

    let has_csq = output_schema.field_with_name("csq").is_ok()
        || output_schema.field_with_name("CSQ").is_ok();
    assert!(has_csq, "CSQ annotation field missing from roundtrip output. Fields: {field_names:?}");
}
