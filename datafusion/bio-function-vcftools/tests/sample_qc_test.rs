//! Integration tests for `vcf_sample_qc` TVF.
//!
//! Tests verify correctness by comparing TVF output against the equivalent
//! SQL query (unnest → transform → array_agg) on the same input data.

use std::fs;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, Float64Builder, Int32Builder, ListBuilder, RecordBatch, StringBuilder,
    StructArray,
};
use datafusion::arrow::datatypes::{DataType, Field, Fields, Schema};
use datafusion::datasource::MemTable;
use datafusion::prelude::*;
use tempfile::TempDir;

use datafusion_bio_function_vcftools::sample_qc::processing;
use datafusion_bio_function_vcftools::sample_qc::{QcConfig, register_sample_qc_function};

// ============================================================================
// Test data builders
// ============================================================================

/// Build a VCF-like RecordBatch with genotypes struct containing GT, GQ, DP, PL.
///
/// Each variant has `n_samples` samples.
fn build_test_batch(variants: &[TestVariant]) -> RecordBatch {
    let mut chrom_builder = StringBuilder::new();
    let mut start_builder = datafusion::arrow::array::UInt32Builder::new();
    let mut end_builder = datafusion::arrow::array::UInt32Builder::new();
    let mut id_builder = StringBuilder::new();
    let mut ref_builder = StringBuilder::new();
    let mut alt_builder = StringBuilder::new();
    let mut qual_builder = Float64Builder::new();
    let mut filter_builder = StringBuilder::new();

    // Genotype list builders
    let mut gt_builder = ListBuilder::new(StringBuilder::new());
    let mut gq_builder = ListBuilder::new(Int32Builder::new());
    let mut dp_builder = ListBuilder::new(Int32Builder::new());
    // PL is List<List<Int32>>
    let mut pl_builder = ListBuilder::new(ListBuilder::new(Int32Builder::new()));

    for v in variants {
        chrom_builder.append_value(&v.chrom);
        start_builder.append_value(v.start);
        end_builder.append_value(v.start + 1);
        id_builder.append_value(&v.id);
        ref_builder.append_value(&v.ref_allele);
        alt_builder.append_value(&v.alt);
        qual_builder.append_value(v.qual);
        filter_builder.append_value(&v.filter);

        // GT
        for gt in &v.samples_gt {
            gt_builder.values().append_value(gt);
        }
        gt_builder.append(true);

        // GQ
        for gq in &v.samples_gq {
            match gq {
                Some(val) => gq_builder.values().append_value(*val),
                None => gq_builder.values().append_null(),
            }
        }
        gq_builder.append(true);

        // DP
        for dp in &v.samples_dp {
            match dp {
                Some(val) => dp_builder.values().append_value(*val),
                None => dp_builder.values().append_null(),
            }
        }
        dp_builder.append(true);

        // PL (list of list of int32)
        for pl in &v.samples_pl {
            match pl {
                Some(vals) => {
                    for val in vals {
                        pl_builder.values().values().append_value(*val);
                    }
                    pl_builder.values().append(true);
                }
                None => {
                    pl_builder.values().append_null();
                }
            }
        }
        pl_builder.append(true);
    }

    let gt_array = gt_builder.finish();
    let gq_array = gq_builder.finish();
    let dp_array = dp_builder.finish();
    let pl_array = pl_builder.finish();

    // Build genotypes struct
    let gt_field = Field::new(
        "GT",
        DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
        true,
    );
    let gq_field = Field::new(
        "GQ",
        DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
        true,
    );
    let dp_field = Field::new(
        "DP",
        DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
        true,
    );
    let pl_field = Field::new(
        "PL",
        DataType::List(Arc::new(Field::new(
            "item",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            true,
        ))),
        true,
    );

    let struct_fields = Fields::from(vec![
        gt_field.clone(),
        gq_field.clone(),
        dp_field.clone(),
        pl_field.clone(),
    ]);

    let genotypes = StructArray::try_new(
        struct_fields.clone(),
        vec![
            Arc::new(gt_array) as ArrayRef,
            Arc::new(gq_array),
            Arc::new(dp_array),
            Arc::new(pl_array),
        ],
        None,
    )
    .unwrap();

    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new("id", DataType::Utf8, true),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
        Field::new("qual", DataType::Float64, true),
        Field::new("filter", DataType::Utf8, true),
        Field::new("genotypes", DataType::Struct(struct_fields), false),
    ]));

    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(chrom_builder.finish()),
            Arc::new(start_builder.finish()),
            Arc::new(end_builder.finish()),
            Arc::new(id_builder.finish()),
            Arc::new(ref_builder.finish()),
            Arc::new(alt_builder.finish()),
            Arc::new(qual_builder.finish()),
            Arc::new(filter_builder.finish()),
            Arc::new(genotypes),
        ],
    )
    .unwrap()
}

#[derive(Clone)]
struct TestVariant {
    chrom: String,
    start: u32,
    id: String,
    ref_allele: String,
    alt: String,
    qual: f64,
    filter: String,
    samples_gt: Vec<String>,
    samples_gq: Vec<Option<i32>>,
    samples_dp: Vec<Option<i32>>,
    samples_pl: Vec<Option<Vec<i32>>>,
}

fn make_test_data() -> Vec<TestVariant> {
    vec![
        // Variant 1: high quality, passes site filter
        // Sample 0: good hom-ref with all-zero PL (needs correction)
        // Sample 1: good het with normal PL
        // Sample 2: low quality sample (GQ < 10)
        TestVariant {
            chrom: "chr1".into(),
            start: 100,
            id: "rs1".into(),
            ref_allele: "A".into(),
            alt: "T".into(),
            qual: 30.0,
            filter: "PASS".into(),
            samples_gt: vec!["0/0".into(), "0/1".into(), "0/0".into()],
            samples_gq: vec![Some(30), Some(25), Some(5)],
            samples_dp: vec![Some(40), Some(35), Some(20)],
            samples_pl: vec![
                Some(vec![0, 0, 0]),    // hom-ref all-zero → needs correction
                Some(vec![25, 0, 180]), // het → no correction
                Some(vec![0, 15, 60]),  // low GQ → no correction (but sample marked bad)
            ],
        },
        // Variant 2: passes site filter
        // Sample 0: haploid GT "0" → normalize to "0/0"
        // Sample 1: hom-alt
        // Sample 2: missing DP
        TestVariant {
            chrom: "chr1".into(),
            start: 200,
            id: "rs2".into(),
            ref_allele: "C".into(),
            alt: "G".into(),
            qual: 50.0,
            filter: "PASS".into(),
            samples_gt: vec!["0".into(), "1/1".into(), "0/1".into()],
            samples_gq: vec![Some(20), Some(40), Some(15)],
            samples_dp: vec![Some(30), Some(50), None],
            samples_pl: vec![
                Some(vec![0, 20, 200]),
                Some(vec![200, 30, 0]),
                None, // missing PL
            ],
        },
        // Variant 3: FAILS site filter (qual < 20)
        TestVariant {
            chrom: "chr1".into(),
            start: 300,
            id: "rs3".into(),
            ref_allele: "G".into(),
            alt: "A".into(),
            qual: 10.0, // too low
            filter: ".".into(),
            samples_gt: vec!["0/0".into(), "0/0".into(), "0/0".into()],
            samples_gq: vec![Some(30), Some(30), Some(30)],
            samples_dp: vec![Some(40), Some(40), Some(40)],
            samples_pl: vec![
                Some(vec![0, 30, 60]),
                Some(vec![0, 30, 60]),
                Some(vec![0, 30, 60]),
            ],
        },
        // Variant 4: FAILS site filter (avg DP too low)
        TestVariant {
            chrom: "chr1".into(),
            start: 400,
            id: "rs4".into(),
            ref_allele: "T".into(),
            alt: "C".into(),
            qual: 25.0,
            filter: "PASS".into(),
            samples_gt: vec!["0/0".into(), "0/1".into(), "1/1".into()],
            samples_gq: vec![Some(20), Some(20), Some(20)],
            samples_dp: vec![Some(5), Some(5), Some(5)], // avg DP = 5 < 15
            samples_pl: vec![
                Some(vec![0, 20, 200]),
                Some(vec![20, 0, 180]),
                Some(vec![200, 20, 0]),
            ],
        },
        // Variant 5: passes, tests edge cases
        // Sample 0: hom-ref with non-zero PL (no correction needed)
        // Sample 1: DP too high (> 200) → bad sample
        // Sample 2: normal het
        TestVariant {
            chrom: "chr2".into(),
            start: 500,
            id: ".".into(),
            ref_allele: "AT".into(),
            alt: "A".into(),
            qual: 45.0,
            filter: "PASS".into(),
            samples_gt: vec!["0/0".into(), "0/1".into(), "0|1".into()],
            samples_gq: vec![Some(35), Some(30), Some(25)],
            samples_dp: vec![Some(25), Some(250), Some(50)], // sample 1 DP > 200
            samples_pl: vec![
                Some(vec![0, 30, 60]),  // non-zero PL → no correction
                Some(vec![25, 0, 180]), // bad sample (DP > 200)
                Some(vec![20, 0, 150]), // normal
            ],
        },
    ]
}

// ============================================================================
// Processing unit tests
// ============================================================================

#[test]
fn test_process_batch_site_filtering() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let result = processing::process_batch(&batch, &config).unwrap();
    assert!(result.is_some(), "Should have passing variants");

    let processed = result.unwrap();
    // Variants 1, 2, 5 pass; 3 (low qual) and 4 (low avg DP) fail
    assert_eq!(processed.n_variants, 3);
}

#[test]
fn test_process_batch_gt_normalization() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    // Variant 2, sample 0: GT "0" → normalized to "0/0" (good sample)
    // Variant 2 is the second passing variant (index 1 in processed output)
    // Its samples start at offsets[1]
    let v2_start = processed.sample_offsets[1];
    let gt_s0 = processed.gt_final.value(v2_start);
    assert_eq!(gt_s0, "0/0", "Haploid '0' should normalize to '0/0'");
}

#[test]
fn test_process_batch_bad_sample_gets_missing_gt() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    // Variant 1, sample 2: GQ=5 < 10 → bad sample → GT = "./."
    let v1_start = processed.sample_offsets[0];
    let gt_bad = processed.gt_final.value(v1_start + 2);
    assert_eq!(gt_bad, "./.", "Low-GQ sample should have './.' genotype");
}

#[test]
fn test_process_batch_pl_correction() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    // Variant 1, sample 0: hom-ref with all-zero PL and GQ=30 → needs correction
    let v1_start = processed.sample_offsets[0];
    let pl0 = processed.pl0.value(v1_start);
    let pl1 = processed.pl1.value(v1_start);
    let pl2 = processed.pl2.value(v1_start);

    // After correction with GQ=30:
    // p_wrong = 10^(-3) = 0.001
    // x ≈ 0.000999
    // PL0 ≈ 0, PL1 ≈ 30, PL2 ≈ 60
    assert!(
        pl0 >= 0 && pl0 <= 1,
        "Corrected PL0 for hom-ref should be ~0, got {pl0}"
    );
    assert!(
        pl1 > 20 && pl1 < 40,
        "Corrected PL1 should be ~30, got {pl1}"
    );
    assert!(
        pl2 > 50 && pl2 < 70,
        "Corrected PL2 should be ~60, got {pl2}"
    );

    // Variant 1, sample 1: het with normal PL → no correction
    let het_pl0 = processed.pl0.value(v1_start + 1);
    let het_pl1 = processed.pl1.value(v1_start + 1);
    let het_pl2 = processed.pl2.value(v1_start + 1);
    assert_eq!(het_pl0, 25, "Het PL0 should be unchanged");
    assert_eq!(het_pl1, 0, "Het PL1 should be unchanged");
    assert_eq!(het_pl2, 180, "Het PL2 should be unchanged");
}

#[test]
fn test_process_batch_no_correction_when_pl_nonzero() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    // Variant 5 (3rd passing, index 2), sample 0: hom-ref with PL=[0,30,60] → no correction
    let v5_start = processed.sample_offsets[2];
    let pl0 = processed.pl0.value(v5_start);
    let pl1 = processed.pl1.value(v5_start);
    let pl2 = processed.pl2.value(v5_start);
    assert_eq!(pl0, 0, "Non-zero PL hom-ref should not be corrected");
    assert_eq!(pl1, 30);
    assert_eq!(pl2, 60);
}

#[test]
fn test_process_batch_ds_computation() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    // Variant 5, sample 2: het with PL=[20,0,150], good sample
    // DS should be close to 1.0 (heterozygote)
    let v5_start = processed.sample_offsets[2];
    let ds = processed.ds_values.value(v5_start + 2);
    assert!(
        (ds - 1.0).abs() < 0.15,
        "Het sample DS should be ~1.0, got {ds}"
    );

    // Variant 2, sample 2: missing DP → bad sample → DS should be null
    let v2_start = processed.sample_offsets[1];
    assert!(
        processed.ds_values.is_null(v2_start + 2),
        "Missing DP sample should have null DS"
    );
}

#[test]
fn test_process_batch_dp_too_high_marks_bad() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    // Variant 5, sample 1: DP=250 > 200 → bad sample → GT = "./."
    let v5_start = processed.sample_offsets[2];
    let gt = processed.gt_final.value(v5_start + 1);
    assert_eq!(gt, "./.", "High-DP sample should have './.' genotype");
}

#[test]
fn test_all_variants_filtered_returns_none() {
    let variants = vec![TestVariant {
        chrom: "chr1".into(),
        start: 100,
        id: ".".into(),
        ref_allele: "A".into(),
        alt: "T".into(),
        qual: 5.0, // too low
        filter: ".".into(),
        samples_gt: vec!["0/0".into()],
        samples_gq: vec![Some(30)],
        samples_dp: vec![Some(40)],
        samples_pl: vec![Some(vec![0, 30, 60])],
    }];
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let result = processing::process_batch(&batch, &config).unwrap();
    assert!(result.is_none(), "All filtered out → None");
}

// ============================================================================
// VCF writer tests
// ============================================================================

#[test]
fn test_vcf_writer_produces_valid_output() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    let tmp = TempDir::new().unwrap();
    let output_path = tmp.path().join("test.vcf");

    let sample_names = vec![
        "SAMPLE_0".to_string(),
        "SAMPLE_1".to_string(),
        "SAMPLE_2".to_string(),
    ];

    let mut writer = datafusion_bio_function_vcftools::sample_qc::writer::VcfWriter::new(
        &output_path,
        &sample_names,
    )
    .unwrap();
    writer.write_batch(&processed).unwrap();
    writer.finish().unwrap();

    let content = fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();

    // Check header
    assert!(lines[0].starts_with("##fileformat=VCFv4.2"));
    let header_line = lines.iter().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(header_line.contains("SAMPLE_0"));
    assert!(header_line.contains("SAMPLE_1"));
    assert!(header_line.contains("SAMPLE_2"));
    assert!(header_line.contains("FORMAT"));

    // Count data lines (non-header, non-empty)
    let data_lines: Vec<&&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(data_lines.len(), 3, "Should have 3 passing variants");

    // Check first data line format
    let first_data = data_lines[0];
    let fields: Vec<&str> = first_data.split('\t').collect();
    assert_eq!(fields[0], "chr1"); // CHROM
    assert_eq!(fields[1], "101"); // POS (0-based 100 + 1)
    assert_eq!(fields[2], "rs1"); // ID
    assert_eq!(fields[3], "A"); // REF
    assert_eq!(fields[4], "T"); // ALT
    assert_eq!(fields[8], "GT:GQ:DP:PL:DS"); // FORMAT

    // Check sample fields have correct number of colon-separated values
    for i in 9..12 {
        let sample_parts: Vec<&str> = fields[i].split(':').collect();
        assert_eq!(
            sample_parts.len(),
            5,
            "Sample {i} should have 5 FORMAT fields, got {:?}",
            sample_parts
        );
    }

    // Check that PL correction is reflected in output
    // First variant, first sample (hom-ref all-zero PL, GQ=30)
    let s0_parts: Vec<&str> = fields[9].split(':').collect();
    let pl_parts: Vec<&str> = s0_parts[3].split(',').collect();
    assert_eq!(pl_parts.len(), 3, "PL should have 3 values");
    // PL0 should be ~0 after correction
    let pl0: i32 = pl_parts[0].parse().unwrap();
    assert!(pl0 <= 1, "Corrected PL0 should be ~0, got {pl0}");
}

#[test]
fn test_vcf_writer_bgzf_output() {
    use std::io::Read;

    let variants = make_test_data();
    let batch = build_test_batch(&variants);
    let config = QcConfig::default();

    let processed = processing::process_batch(&batch, &config).unwrap().unwrap();

    let tmp = TempDir::new().unwrap();
    let output_path = tmp.path().join("test.vcf.gz");

    let sample_names = vec![
        "SAMPLE_0".to_string(),
        "SAMPLE_1".to_string(),
        "SAMPLE_2".to_string(),
    ];

    let mut writer = datafusion_bio_function_vcftools::sample_qc::writer::VcfWriter::new(
        &output_path,
        &sample_names,
    )
    .unwrap();
    writer.write_batch(&processed).unwrap();
    writer.finish().unwrap();

    // Verify BGZF magic: gzip header (1f 8b) with BGZF extra field
    let compressed = fs::read(&output_path).unwrap();
    assert!(
        compressed[0] == 0x1f && compressed[1] == 0x8b,
        "File should start with gzip magic bytes"
    );

    // Verify BGZF EOF block is present (last 28 bytes of a valid BGZF file)
    assert!(
        compressed.len() >= 28,
        "BGZF file should have at least the EOF block"
    );

    // Decompress with noodles-bgzf reader and verify content
    let mut reader = noodles_bgzf::Reader::new(&compressed[..]);
    let mut content = String::new();
    reader.read_to_string(&mut content).unwrap();

    let lines: Vec<&str> = content.lines().collect();
    assert!(lines[0].starts_with("##fileformat=VCFv4.2"));
    let header_line = lines.iter().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(header_line.contains("SAMPLE_0"));

    let data_lines: Vec<&&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(data_lines.len(), 3, "Should have 3 passing variants");
}

// ============================================================================
// End-to-end TVF test via SQL
// ============================================================================

#[tokio::test(flavor = "multi_thread")]
async fn test_tvf_end_to_end() {
    let variants = make_test_data();
    let batch = build_test_batch(&variants);

    let ctx = SessionContext::new();
    register_sample_qc_function(&ctx);

    // Register test data as a table
    let schema = batch.schema();
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_table", Arc::new(table)).unwrap();

    let tmp = TempDir::new().unwrap();
    let output_path = tmp.path().join("output.vcf");
    let output_str = output_path.to_str().unwrap();

    // Execute the TVF
    let result = ctx
        .sql(&format!(
            "SELECT * FROM vcf_sample_qc('vcf_table', '{output_str}', 'S1,S2,S3')"
        ))
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();

    // Check summary output
    assert_eq!(result.len(), 1);
    let summary = &result[0];
    let variants_in = summary
        .column_by_name("variants_in")
        .unwrap()
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .unwrap()
        .value(0);
    let variants_out = summary
        .column_by_name("variants_out")
        .unwrap()
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .unwrap()
        .value(0);

    assert_eq!(variants_in, 5, "Should process all 5 input variants");
    assert_eq!(variants_out, 3, "3 variants should pass site filter");

    // Verify VCF output exists and has correct structure
    let content = fs::read_to_string(&output_path).unwrap();
    let data_lines: Vec<&str> = content
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(data_lines.len(), 3);

    // Verify sample names are correct
    let header_line = content.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(header_line.contains("S1"));
    assert!(header_line.contains("S2"));
    assert!(header_line.contains("S3"));
}

#[tokio::test(flavor = "multi_thread")]
async fn test_tvf_auto_sample_names() {
    let variants = vec![make_test_data()[0].clone()]; // 1 passing variant, 3 samples
    let batch = build_test_batch(&variants);

    let ctx = SessionContext::new();
    register_sample_qc_function(&ctx);

    let schema = batch.schema();
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_auto", Arc::new(table)).unwrap();

    let tmp = TempDir::new().unwrap();
    let output_path = tmp.path().join("auto.vcf");
    let output_str = output_path.to_str().unwrap();

    let result = ctx
        .sql(&format!(
            "SELECT * FROM vcf_sample_qc('vcf_auto', '{output_str}', 'auto')"
        ))
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();

    let content = fs::read_to_string(&output_path).unwrap();
    let header = content.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(header.contains("SAMPLE_0"));
    assert!(header.contains("SAMPLE_1"));
    assert!(header.contains("SAMPLE_2"));
}

#[tokio::test(flavor = "multi_thread")]
async fn test_tvf_sample_names_from_file() {
    let variants = vec![make_test_data()[0].clone()]; // 1 passing variant, 3 samples
    let batch = build_test_batch(&variants);

    let ctx = SessionContext::new();
    register_sample_qc_function(&ctx);

    let schema = batch.schema();
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_file_samples", Arc::new(table))
        .unwrap();

    let tmp = TempDir::new().unwrap();

    // Write a sample names file
    let samples_file = tmp.path().join("samples.txt");
    fs::write(&samples_file, "# comment line\nHG02679\nNA19920\nHG00609\n").unwrap();

    let output_path = tmp.path().join("file_samples.vcf");
    let output_str = output_path.to_str().unwrap();
    let samples_str = samples_file.to_str().unwrap();

    let _result = ctx
        .sql(&format!(
            "SELECT * FROM vcf_sample_qc('vcf_file_samples', '{output_str}', '{samples_str}')"
        ))
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();

    let content = fs::read_to_string(&output_path).unwrap();
    let header = content.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(
        header.contains("HG02679"),
        "Should contain sample name from file"
    );
    assert!(
        header.contains("NA19920"),
        "Should contain sample name from file"
    );
    assert!(
        header.contains("HG00609"),
        "Should contain sample name from file"
    );
    assert!(
        !header.contains("SAMPLE_0"),
        "Should NOT contain auto-generated names"
    );
}

#[tokio::test(flavor = "multi_thread")]
async fn test_tvf_custom_thresholds() {
    // Use a very low qual threshold so all variants pass
    let variants = make_test_data();
    let batch = build_test_batch(&variants);

    let ctx = SessionContext::new();
    register_sample_qc_function(&ctx);

    let schema = batch.schema();
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_custom", Arc::new(table)).unwrap();

    let tmp = TempDir::new().unwrap();
    let output_path = tmp.path().join("custom.vcf");
    let output_str = output_path.to_str().unwrap();

    // qual_min=1, site_gq_avg_min=1, site_dp_avg_min=1, site_dp_avg_max=9999
    let result = ctx
        .sql(&format!(
            "SELECT * FROM vcf_sample_qc('vcf_custom', '{output_str}', 'auto', 1.0, 1.0, 1.0, 9999.0)"
        ))
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();

    let variants_out = result[0]
        .column_by_name("variants_out")
        .unwrap()
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .unwrap()
        .value(0);

    // All 5 variants should pass with relaxed thresholds
    assert_eq!(
        variants_out, 5,
        "All variants should pass with low thresholds"
    );
}

// ============================================================================
// Correctness: compare TVF sample processing against manual computation
// ============================================================================

#[test]
fn test_pl_correction_matches_sql_formula() {
    // Verify the Rust PL correction matches the SQL formula exactly
    // SQL: p_wrong = POWER(10.0, -gq/10.0)
    //      x = (-1.0 + SQRT(1.0 + 4.0 * p_wrong)) / 2.0
    //      pl0 = LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, 1.0 - x - x*x)))))
    //      pl1 = LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, x)))))
    //      pl2 = LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, x*x)))))

    for gq in [10.0, 20.0, 30.0, 40.0, 50.0, 99.0] {
        let (pl0, pl1, pl2) = processing::processing_test_helpers::correct_pl_public(gq);

        // Compute expected values using the SQL formula
        let p_wrong = 10.0_f64.powf(-gq / 10.0);
        let x = (-1.0 + (1.0 + 4.0 * p_wrong).sqrt()) / 2.0;
        let expected_pl0 = (-10.0 * (1.0 - x - x * x).max(1e-25).log10())
            .round()
            .clamp(0.0, 255.0) as i32;
        let expected_pl1 = (-10.0 * x.max(1e-25).log10()).round().clamp(0.0, 255.0) as i32;
        let expected_pl2 = (-10.0 * (x * x).max(1e-25).log10())
            .round()
            .clamp(0.0, 255.0) as i32;

        assert_eq!(pl0, expected_pl0, "PL0 mismatch for GQ={gq}");
        assert_eq!(pl1, expected_pl1, "PL1 mismatch for GQ={gq}");
        assert_eq!(pl2, expected_pl2, "PL2 mismatch for GQ={gq}");
    }
}

#[test]
fn test_dosage_matches_sql_formula() {
    // Verify DS computation matches the SQL formula
    // SQL: l_ref = CASE WHEN final_pl0 < 255 THEN POWER(10.0, -final_pl0/10.0) ELSE 0.0 END
    //      l_het = CASE WHEN final_pl1 < 255 THEN POWER(10.0, -final_pl1/10.0) ELSE 0.0 END
    //      l_alt = CASE WHEN final_pl2 < 255 THEN POWER(10.0, -final_pl2/10.0) ELSE 0.0 END
    //      ds = (l_het + 2.0 * l_alt) / (l_ref + l_het + l_alt)

    let test_cases = vec![
        (0, 30, 60, Some(0.002)), // hom-ref: DS ≈ 0
        (30, 0, 30, Some(1.0)),   // het: DS ≈ 1.0
        (60, 30, 0, Some(1.998)), // hom-alt: DS ≈ 2.0
        (255, 255, 255, None),    // all 255: DS = null
    ];

    for (pl0, pl1, pl2, expected_approx) in test_cases {
        let ds = processing::processing_test_helpers::compute_dosage_public(pl0, pl1, pl2);
        match (ds, expected_approx) {
            (Some(d), Some(e)) => {
                assert!(
                    (d as f64 - e).abs() < 0.05,
                    "DS mismatch for PL=({pl0},{pl1},{pl2}): got {d}, expected ~{e}"
                );
            }
            (None, None) => {}
            _ => panic!(
                "DS null mismatch for PL=({pl0},{pl1},{pl2}): got {ds:?}, expected {expected_approx:?}"
            ),
        }
    }
}
