//! Roundtrip correctness test: annotate 1000 variants → write VCF → read back
//! → compare ALL column values against golden Ensembl VEP 115 output.
//!
//! Verifies 100% correctness with 0 mismatches for:
//! - Core VCF columns (chrom, pos, ref, alt, qual, filter, id)
//! - Original INFO fields (platforms, datasets, etc.) preserved from input
//! - CSQ annotation field matches golden VEP 115 (80 pipe-delimited fields)
//! - FORMAT/sample fields (GT, DP, GQ, AD, etc.) preserved from input

use std::sync::Arc;

use datafusion::arrow::array::Array;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::vcf_sink;

/// Check if a file is a Git LFS pointer (not actual content).
fn is_lfs_pointer(path: &std::path::Path) -> bool {
    std::fs::read_to_string(path)
        .map(|s| s.starts_with("version https://git-lfs.github.com"))
        .unwrap_or(false)
}

/// Resolve a path relative to the workspace root.
fn workspace_path(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

#[tokio::test(flavor = "multi_thread")]
async fn test_roundtrip_golden_all_column_values() {
    let input_vcf = workspace_path("vep-benchmark/data/golden/input_1000.vcf");
    let golden_vcf = workspace_path("vep-benchmark/data/golden/golden_1000_vep115.vcf");
    let cache_path = workspace_path("vep-benchmark/data/golden/cache");
    let ref_fasta = workspace_path("vep-benchmark/data/golden/reference_chr1.fa");

    if !input_vcf.exists() || !golden_vcf.exists() || is_lfs_pointer(&input_vcf) {
        eprintln!(
            "Skipping: test fixtures not found at {}",
            input_vcf.display()
        );
        return;
    }
    let input_vcf = input_vcf.to_str().unwrap();
    let golden_vcf = golden_vcf.to_str().unwrap();
    let cache_path = cache_path.to_str().unwrap();
    let ref_fasta = ref_fasta.to_str().unwrap();

    // ── Step 1: Annotate and write to VCF ──────────────────────────
    let tmp_dir = tempfile::TempDir::new().unwrap();
    let output_path = tmp_dir.path().join("annotated.vcf");

    let config = vcf_sink::AnnotateVcfConfig {
        everything: true,
        extended_probes: true,
        reference_fasta_path: Some(ref_fasta.to_string()),
        ..Default::default()
    };
    let rows_written = vcf_sink::annotate_to_vcf(
        input_vcf,
        cache_path,
        "parquet",
        output_path.to_str().unwrap(),
        &config,
    )
    .await
    .unwrap();
    assert_eq!(rows_written, 1000, "Should write 1000 annotated rows");

    // ── Step 1b: Verify CSQ INFO header contains Format: field list ──
    {
        let vcf_text = std::fs::read_to_string(&output_path).unwrap();
        let csq_header = vcf_text
            .lines()
            .find(|l| l.starts_with("##INFO=<ID=CSQ,"))
            .expect("Output VCF must have ##INFO=<ID=CSQ,...> header line");
        assert!(
            csq_header.contains("Format: "),
            "CSQ INFO header must contain 'Format: ' field list, got: {csq_header}"
        );
        // Verify field list matches the --everything constant (80 fields).
        let format_str = csq_header
            .split("Format: ")
            .nth(1)
            .unwrap()
            .trim_end_matches("\">");
        let fields: Vec<&str> = format_str.split('|').collect();
        let expected = datafusion_bio_function_vep::golden_benchmark::CSQ_FIELD_NAMES_EVERYTHING;
        assert_eq!(
            fields.len(),
            expected.len(),
            "CSQ Format field count mismatch: got {}, expected {}",
            fields.len(),
            expected.len()
        );
        assert_eq!(fields, expected, "CSQ Format field names/order mismatch");
    }

    // ── Step 2: Read input, output, and golden with VcfTableProvider ──
    let ctx2 = SessionContext::new();

    let input_str2 = input_vcf.to_string();
    let input_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(input_str2, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    // VcfTableProvider::new() uses futures::executor::block_on() internally,
    // which deadlocks inside a tokio runtime. Use spawn_blocking to avoid.
    let out_path_str = output_path.display().to_string();
    let output_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(out_path_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    let golden_str = golden_vcf.to_string();
    let golden_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(golden_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();

    ctx2.register_table("input_vcf", Arc::new(input_prov))
        .unwrap();
    ctx2.register_table("output_vcf", Arc::new(output_prov))
        .unwrap();
    ctx2.register_table("golden_vcf", Arc::new(golden_prov))
        .unwrap();

    // Collect all three as sorted by position for deterministic comparison.
    let input_batches = ctx2
        .sql("SELECT * FROM input_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let output_batches = ctx2
        .sql("SELECT * FROM output_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let golden_batches = ctx2
        .sql("SELECT * FROM golden_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();

    let input_rows: usize = input_batches.iter().map(|b| b.num_rows()).sum();
    let output_rows: usize = output_batches.iter().map(|b| b.num_rows()).sum();
    let golden_rows: usize = golden_batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(output_rows, 1000, "Output should have 1000 rows");
    assert_eq!(golden_rows, 1000, "Golden should have 1000 rows");
    assert_eq!(input_rows, 1000, "Input should have 1000 rows");

    // Concatenate into single batches for easier row-by-row comparison.
    let input_batch =
        datafusion::arrow::compute::concat_batches(&input_batches[0].schema(), &input_batches)
            .unwrap();
    let output_batch =
        datafusion::arrow::compute::concat_batches(&output_batches[0].schema(), &output_batches)
            .unwrap();
    let golden_batch =
        datafusion::arrow::compute::concat_batches(&golden_batches[0].schema(), &golden_batches)
            .unwrap();

    // ── Step 3: Compare core VCF columns (output vs input) ─────────
    // Use Arrow array equality for fast column-level comparison.
    let core_columns = [
        "chrom", "start", "end", "ref", "alt", "id", "qual", "filter",
    ];
    for col_name in &core_columns {
        let in_idx = input_batch.schema().index_of(col_name);
        let out_idx = output_batch.schema().index_of(col_name);
        if let (Ok(i), Ok(o)) = (in_idx, out_idx) {
            let in_col = input_batch.column(i);
            let out_col = output_batch.column(o);
            assert_eq!(
                in_col.as_ref(),
                out_col.as_ref(),
                "Core column '{col_name}' values differ between input and output"
            );
        }
    }

    // ── Step 4: Compare original INFO fields (output vs input) ─────
    let input_schema = input_batch.schema();
    let info_columns: Vec<String> = input_schema
        .fields()
        .iter()
        .filter(|f| {
            f.metadata()
                .get("bio.vcf.field.field_type")
                .is_some_and(|v| v == "INFO")
        })
        .map(|f| f.name().clone())
        .collect();
    eprintln!(
        "Comparing {} INFO columns: {:?}",
        info_columns.len(),
        info_columns
    );

    let mut info_missing = Vec::new();
    let mut info_mismatched = Vec::new();
    for col_name in &info_columns {
        if output_batch.schema().index_of(col_name).is_err() {
            info_missing.push(col_name.clone());
            continue;
        }
        let in_col = input_batch.column(input_batch.schema().index_of(col_name).unwrap());
        let out_col = output_batch.column(output_batch.schema().index_of(col_name).unwrap());
        if in_col.as_ref() != out_col.as_ref() {
            info_mismatched.push(col_name.clone());
        }
    }
    assert!(
        info_missing.is_empty(),
        "INFO columns missing from output: {info_missing:?}"
    );
    assert!(
        info_mismatched.is_empty(),
        "INFO columns with value mismatches: {info_mismatched:?}"
    );

    // ── Step 5: Compare FORMAT/sample fields (output vs input) ─────
    let format_columns: Vec<String> = input_schema
        .fields()
        .iter()
        .filter(|f| {
            f.metadata()
                .get("bio.vcf.field.field_type")
                .is_some_and(|v| v == "FORMAT")
        })
        .map(|f| f.name().clone())
        .collect();
    eprintln!(
        "Comparing {} FORMAT columns: {:?}",
        format_columns.len(),
        format_columns
    );

    let mut format_missing = Vec::new();
    let mut format_mismatched = Vec::new();
    for col_name in &format_columns {
        if output_batch.schema().index_of(col_name).is_err() {
            format_missing.push(col_name.clone());
            continue;
        }
        let in_col = input_batch.column(input_batch.schema().index_of(col_name).unwrap());
        let out_col = output_batch.column(output_batch.schema().index_of(col_name).unwrap());
        if in_col.as_ref() != out_col.as_ref() {
            format_mismatched.push(col_name.clone());
        }
    }
    assert!(
        format_missing.is_empty(),
        "FORMAT columns missing from output: {format_missing:?}"
    );
    assert!(
        format_mismatched.is_empty(),
        "FORMAT columns with value mismatches: {format_mismatched:?}"
    );
    // ── Step 6: Compare CSQ annotation (output vs golden VEP 115) ──
    let our_csq_col = output_batch.column(
        output_batch
            .schema()
            .index_of("CSQ")
            .expect("CSQ column not found in output"),
    );
    let golden_csq_col = golden_batch.column(
        golden_batch
            .schema()
            .index_of("CSQ")
            .expect("CSQ column not found in golden"),
    );

    // First try fast path: exact array equality.
    if our_csq_col.as_ref() == golden_csq_col.as_ref() {
    } else {
        // Per-row string comparison for detailed diagnostics.
        // CSQ may be StringArray or StringViewArray depending on the provider.
        let get_str = |col: &dyn Array, row: usize| -> String {
            if col.is_null(row) {
                return String::new();
            }
            if let Some(a) = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringArray>()
            {
                return a.value(row).to_string();
            }
            if let Some(a) = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringViewArray>()
            {
                return a.value(row).to_string();
            }
            String::new()
        };

        let mut mismatched_rows = 0;
        for row in 0..1000 {
            let ours = get_str(our_csq_col.as_ref(), row);
            let golden = get_str(golden_csq_col.as_ref(), row);
            if ours != golden {
                mismatched_rows += 1;
                if mismatched_rows <= 3 {
                    eprintln!("CSQ mismatch row {row}:");
                    eprintln!("  golden: {}...", &golden[..golden.len().min(150)]);
                    eprintln!("  ours:   {}...", &ours[..ours.len().min(150)]);
                }
            }
        }
        assert_eq!(
            mismatched_rows, 0,
            "CSQ should have 0 row-level mismatches vs golden VEP 115 ({mismatched_rows} mismatched)"
        );
    }
}
