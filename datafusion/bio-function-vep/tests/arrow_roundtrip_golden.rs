//! Arrow output correctness test: annotate 1000 variants via `annotate_vep()`
//! table function → collect typed Arrow columns → verify that every transcript-
//! level list column contains the correct values for ALL transcripts, matched
//! against the pipe-delimited CSQ string emitted by the same call.
//!
//! With the list-typed column change, each transcript-level column (SYMBOL, Gene,
//! Feature, IMPACT, etc.) is `List<T>` with one element per transcript. This test
//! verifies every element in every list matches the corresponding CSQ field.
//!
//! Combined with `vcf_roundtrip_golden.rs` (which proves our CSQ == golden
//! Ensembl VEP 115 at 100%), this transitively proves:
//!     typed Arrow columns == golden VEP 115.

use std::collections::HashMap;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, Float32Array, Int8Array, Int64Array, ListArray, StringArray, StringViewArray,
};
use datafusion::arrow::compute::concat_batches;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;

// ── Helpers ──────────────────────────────────────────────────────────

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

fn get_str(col: &dyn Array, row: usize) -> Option<String> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
        return Some(a.value(row).to_string());
    }
    None
}

fn get_list_str(col: &dyn Array, row: usize) -> Option<Vec<String>> {
    if col.is_null(row) {
        return None;
    }
    let list_arr = col.as_any().downcast_ref::<ListArray>()?;
    let offsets = list_arr.offsets();
    let start = offsets[row] as usize;
    let end = offsets[row + 1] as usize;
    let values = list_arr.values();
    let mut out = Vec::with_capacity(end - start);
    for i in start..end {
        if values.is_null(i) {
            out.push(String::new());
        } else if let Some(a) = values.as_any().downcast_ref::<StringArray>() {
            out.push(a.value(i).to_string());
        } else if let Some(a) = values.as_any().downcast_ref::<StringViewArray>() {
            out.push(a.value(i).to_string());
        } else {
            out.push(String::new());
        }
    }
    Some(out)
}

fn get_list_i8(col: &dyn Array, row: usize) -> Option<Vec<Option<i8>>> {
    if col.is_null(row) {
        return None;
    }
    let list_arr = col.as_any().downcast_ref::<ListArray>()?;
    let offsets = list_arr.offsets();
    let start = offsets[row] as usize;
    let end = offsets[row + 1] as usize;
    let values = list_arr.values();
    let i8_arr = values.as_any().downcast_ref::<Int8Array>()?;
    let mut out = Vec::with_capacity(end - start);
    for i in start..end {
        if i8_arr.is_null(i) {
            out.push(None);
        } else {
            out.push(Some(i8_arr.value(i)));
        }
    }
    Some(out)
}

fn get_list_i64(col: &dyn Array, row: usize) -> Option<Vec<Option<i64>>> {
    if col.is_null(row) {
        return None;
    }
    let list_arr = col.as_any().downcast_ref::<ListArray>()?;
    let offsets = list_arr.offsets();
    let start = offsets[row] as usize;
    let end = offsets[row + 1] as usize;
    let values = list_arr.values();
    let i64_arr = values.as_any().downcast_ref::<Int64Array>()?;
    let mut out = Vec::with_capacity(end - start);
    for i in start..end {
        if i64_arr.is_null(i) {
            out.push(None);
        } else {
            out.push(Some(i64_arr.value(i)));
        }
    }
    Some(out)
}

#[allow(dead_code)]
fn get_list_f32(col: &dyn Array, row: usize) -> Option<Vec<Option<f32>>> {
    if col.is_null(row) {
        return None;
    }
    let list_arr = col.as_any().downcast_ref::<ListArray>()?;
    let offsets = list_arr.offsets();
    let start = offsets[row] as usize;
    let end = offsets[row + 1] as usize;
    let values = list_arr.values();
    let f32_arr = values.as_any().downcast_ref::<Float32Array>()?;
    let mut out = Vec::with_capacity(end - start);
    for i in start..end {
        if f32_arr.is_null(i) {
            out.push(None);
        } else {
            out.push(Some(f32_arr.value(i)));
        }
    }
    Some(out)
}

// ── CSQ field index ──────────────────────────────────────────────────

const CSQ_FIELD_NAMES: &[&str] = &[
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "VARIANT_CLASS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "CANONICAL",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    "GENE_PHENO",
    "SIFT",
    "PolyPhen",
    "DOMAINS",
    "miRNA",
    "HGVS_OFFSET",
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_MID_AF",
    "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF",
    "gnomADg_AF",
    "gnomADg_AFR_AF",
    "gnomADg_AMI_AF",
    "gnomADg_AMR_AF",
    "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF",
    "gnomADg_FIN_AF",
    "gnomADg_MID_AF",
    "gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF",
    "gnomADg_SAS_AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
];

fn csq_field_idx(name: &str) -> usize {
    CSQ_FIELD_NAMES
        .iter()
        .position(|&n| n == name)
        .unwrap_or_else(|| panic!("CSQ field '{name}' not found"))
}

fn parse_csq(csq: &str) -> Vec<Vec<&str>> {
    csq.split(',')
        .map(|entry| entry.split('|').collect())
        .collect()
}

// ── Test ─────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_arrow_typed_columns_vs_golden_vep115() {
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

    let cache_str = cache_path.to_str().unwrap();
    let ref_fasta_str = ref_fasta.to_str().unwrap();

    // ── 1. Annotate via annotate_vep() → typed Arrow columns ─────────
    let session_config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(session_config);
    datafusion_bio_function_vep::register_vep_functions(&ctx);

    let vcf_path = input_vcf.to_str().unwrap().to_string();
    let vcf_provider = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(vcf_path, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    ctx.register_table("vcf_input", Arc::new(vcf_provider))
        .unwrap();

    let options_json = format!(
        r#"{{"partitioned":true,"everything":true,"extended_probes":true,"reference_fasta_path":"{}"}}"#,
        ref_fasta_str.replace('\\', "\\\\").replace('"', "\\\"")
    );
    let sql = format!(
        "SELECT * FROM annotate_vep('vcf_input', '{}', 'parquet', '{}')",
        cache_str.replace('\'', "''"),
        options_json.replace('\'', "''"),
    );

    let df = ctx.sql(&sql).await.expect("annotate_vep SQL should parse");
    let batches = df.collect().await.expect("annotate_vep should execute");
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 1000, "Should produce 1000 annotated rows");

    let schema = batches[0].schema();
    let arrow_batch = concat_batches(&schema, &batches).unwrap();
    eprintln!(
        "Arrow output: {} rows, {} columns",
        arrow_batch.num_rows(),
        schema.fields().len()
    );

    // ── 2. Verify expected columns exist ─────────────────────────────
    for col_name in &[
        "CSQ",
        "most_severe_consequence",
        "Allele",
        "Feature",
        "SYMBOL",
        "Gene",
        "IMPACT",
        "Consequence",
        "AF",
        "MAX_AF",
        "CLIN_SIG",
    ] {
        assert!(
            schema.index_of(col_name).is_ok(),
            "Expected column '{col_name}' missing"
        );
    }

    // ── 3. Per-transcript verification ───────────────────────────────
    //
    // For each variant row:
    //   1. Parse the CSQ string into transcript entries
    //   2. Read the Feature list from Arrow (one element per transcript)
    //   3. For each CSQ entry, find its index in the Feature list
    //   4. Verify every list column at that index matches the CSQ field

    /// Column comparison type
    #[derive(Clone, Copy)]
    enum Cmp {
        ListStr,     // List<Utf8> element ↔ CSQ string field
        ListI8,      // List<Int8> element ↔ CSQ string field
        ListI64,     // List<Int64> element ↔ CSQ string field
        Scalar,      // Scalar Utf8 (Allele, VARIANT_CLASS) — same for all entries
        ScalarFloat, // Scalar Float32 (frequency columns)
        ScalarStr,   // Scalar Utf8 (variant-level: SOMATIC, PHENO, etc.)
        CsvList,     // List<Utf8> (CLIN_SIG, PUBMED) — variant-level, comma-split
    }

    let columns_under_test: &[(&str, &str, Cmp)] = &[
        // ── Scalar columns (same for all transcripts) ──
        ("Allele", "Allele", Cmp::Scalar),
        ("VARIANT_CLASS", "VARIANT_CLASS", Cmp::Scalar),
        // ── Transcript-level List columns ──
        ("Consequence", "Consequence", Cmp::ListStr),
        ("IMPACT", "IMPACT", Cmp::ListStr),
        ("SYMBOL", "SYMBOL", Cmp::ListStr),
        ("Gene", "Gene", Cmp::ListStr),
        ("Feature_type", "Feature_type", Cmp::ListStr),
        ("Feature", "Feature", Cmp::ListStr),
        ("BIOTYPE", "BIOTYPE", Cmp::ListStr),
        ("EXON", "EXON", Cmp::ListStr),
        ("INTRON", "INTRON", Cmp::ListStr),
        ("HGVSc", "HGVSc", Cmp::ListStr),
        ("HGVSp", "HGVSp", Cmp::ListStr),
        ("cDNA_position", "cDNA_position", Cmp::ListStr),
        ("CDS_position", "CDS_position", Cmp::ListStr),
        ("Protein_position", "Protein_position", Cmp::ListStr),
        ("Amino_acids", "Amino_acids", Cmp::ListStr),
        ("Codons", "Codons", Cmp::ListStr),
        ("Existing_variation", "Existing_variation", Cmp::CsvList),
        ("DISTANCE", "DISTANCE", Cmp::ListI64),
        ("STRAND", "STRAND", Cmp::ListI8),
        ("FLAGS", "FLAGS", Cmp::ListStr),
        ("SYMBOL_SOURCE", "SYMBOL_SOURCE", Cmp::ListStr),
        ("HGNC_ID", "HGNC_ID", Cmp::ListStr),
        ("CANONICAL", "CANONICAL", Cmp::ListStr),
        ("MANE", "MANE", Cmp::ListStr),
        ("MANE_SELECT", "MANE_SELECT", Cmp::ListStr),
        ("MANE_PLUS_CLINICAL", "MANE_PLUS_CLINICAL", Cmp::ListStr),
        ("TSL", "TSL", Cmp::ListI8),
        ("APPRIS", "APPRIS", Cmp::ListStr),
        ("CCDS", "CCDS", Cmp::ListStr),
        ("ENSP", "ENSP", Cmp::ListStr),
        ("SWISSPROT", "SWISSPROT", Cmp::ListStr),
        ("TREMBL", "TREMBL", Cmp::ListStr),
        ("UNIPARC", "UNIPARC", Cmp::ListStr),
        ("UNIPROT_ISOFORM", "UNIPROT_ISOFORM", Cmp::ListStr),
        ("GENE_PHENO", "GENE_PHENO", Cmp::ListStr),
        ("SIFT", "SIFT", Cmp::ListStr),
        ("PolyPhen", "PolyPhen", Cmp::ListStr),
        ("DOMAINS", "DOMAINS", Cmp::ListStr),
        ("miRNA", "miRNA", Cmp::ListStr),
        ("HGVS_OFFSET", "HGVS_OFFSET", Cmp::ListI64),
        // ── Frequency (scalar, variant-level) ──
        ("AF", "AF", Cmp::ScalarFloat),
        ("AFR_AF", "AFR_AF", Cmp::ScalarFloat),
        ("AMR_AF", "AMR_AF", Cmp::ScalarFloat),
        ("EAS_AF", "EAS_AF", Cmp::ScalarFloat),
        ("EUR_AF", "EUR_AF", Cmp::ScalarFloat),
        ("SAS_AF", "SAS_AF", Cmp::ScalarFloat),
        ("gnomADe_AF", "gnomADe_AF", Cmp::ScalarFloat),
        ("gnomADe_AFR_AF", "gnomADe_AFR_AF", Cmp::ScalarFloat),
        ("gnomADe_AMR_AF", "gnomADe_AMR_AF", Cmp::ScalarFloat),
        ("gnomADe_ASJ_AF", "gnomADe_ASJ_AF", Cmp::ScalarFloat),
        ("gnomADe_EAS_AF", "gnomADe_EAS_AF", Cmp::ScalarFloat),
        ("gnomADe_FIN_AF", "gnomADe_FIN_AF", Cmp::ScalarFloat),
        ("gnomADe_MID_AF", "gnomADe_MID_AF", Cmp::ScalarFloat),
        ("gnomADe_NFE_AF", "gnomADe_NFE_AF", Cmp::ScalarFloat),
        (
            "gnomADe_REMAINING_AF",
            "gnomADe_REMAINING_AF",
            Cmp::ScalarFloat,
        ),
        ("gnomADe_SAS_AF", "gnomADe_SAS_AF", Cmp::ScalarFloat),
        ("gnomADg_AF", "gnomADg_AF", Cmp::ScalarFloat),
        ("gnomADg_AFR_AF", "gnomADg_AFR_AF", Cmp::ScalarFloat),
        ("gnomADg_AMI_AF", "gnomADg_AMI_AF", Cmp::ScalarFloat),
        ("gnomADg_AMR_AF", "gnomADg_AMR_AF", Cmp::ScalarFloat),
        ("gnomADg_ASJ_AF", "gnomADg_ASJ_AF", Cmp::ScalarFloat),
        ("gnomADg_EAS_AF", "gnomADg_EAS_AF", Cmp::ScalarFloat),
        ("gnomADg_FIN_AF", "gnomADg_FIN_AF", Cmp::ScalarFloat),
        ("gnomADg_MID_AF", "gnomADg_MID_AF", Cmp::ScalarFloat),
        ("gnomADg_NFE_AF", "gnomADg_NFE_AF", Cmp::ScalarFloat),
        (
            "gnomADg_REMAINING_AF",
            "gnomADg_REMAINING_AF",
            Cmp::ScalarFloat,
        ),
        ("gnomADg_SAS_AF", "gnomADg_SAS_AF", Cmp::ScalarFloat),
        ("MAX_AF", "MAX_AF", Cmp::ScalarFloat),
        ("MAX_AF_POPS", "MAX_AF_POPS", Cmp::ScalarStr),
        // ── Variant-level ──
        ("CLIN_SIG", "CLIN_SIG", Cmp::CsvList),
        ("SOMATIC", "SOMATIC", Cmp::ScalarStr),
        ("PHENO", "PHENO", Cmp::ScalarStr),
        ("PUBMED", "PUBMED", Cmp::CsvList),
        ("MOTIF_NAME", "MOTIF_NAME", Cmp::ScalarStr),
        ("MOTIF_POS", "MOTIF_POS", Cmp::ScalarStr),
        ("HIGH_INF_POS", "HIGH_INF_POS", Cmp::ScalarStr),
        ("MOTIF_SCORE_CHANGE", "MOTIF_SCORE_CHANGE", Cmp::ScalarFloat),
        (
            "TRANSCRIPTION_FACTORS",
            "TRANSCRIPTION_FACTORS",
            Cmp::CsvList,
        ),
    ];

    let resolved: Vec<(usize, usize, Cmp, &str)> = columns_under_test
        .iter()
        .filter_map(|&(arrow_name, csq_name, cmp)| {
            schema
                .index_of(arrow_name)
                .ok()
                .map(|ai| (ai, csq_field_idx(csq_name), cmp, arrow_name))
        })
        .collect();

    let csq_feature_fi = csq_field_idx("Feature");
    let our_csq_idx = schema.index_of("CSQ").unwrap();
    let feature_col_idx = schema.index_of("Feature").unwrap();

    let mut col_match: HashMap<&str, usize> = HashMap::new();
    let mut col_mismatch: HashMap<&str, usize> = HashMap::new();
    let mut col_samples: HashMap<&str, Vec<String>> = HashMap::new();
    let mut total_entries_checked = 0usize;

    for row in 0..arrow_batch.num_rows() {
        let our_csq = match get_str(arrow_batch.column(our_csq_idx).as_ref(), row) {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let csq_entries = parse_csq(&our_csq);
        if csq_entries.is_empty() {
            continue;
        }

        // Read the Feature list to build transcript index mapping.
        let feature_list =
            get_list_str(arrow_batch.column(feature_col_idx).as_ref(), row).unwrap_or_default();

        // For each CSQ entry, find its position in the Arrow Feature list.
        for csq_entry in &csq_entries {
            let csq_feature = csq_entry.get(csq_feature_fi).unwrap_or(&"");

            // Find the index of this transcript in the Arrow list.
            let tc_idx = feature_list.iter().position(|f| f == csq_feature);
            let Some(tc_idx) = tc_idx else {
                // Feature not found in list — skip (e.g., intergenic with empty feature)
                continue;
            };

            total_entries_checked += 1;

            for &(arrow_idx, csq_fi, cmp, col_name) in &resolved {
                let csq_val = csq_entry.get(csq_fi).unwrap_or(&"");

                let is_match = match cmp {
                    Cmp::ListStr => {
                        let list = get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_default();
                        let ours = list.get(tc_idx).map(|s| s.as_str()).unwrap_or("");
                        ours == *csq_val || (ours.is_empty() && csq_val.is_empty())
                    }
                    Cmp::ListI8 => {
                        let list = get_list_i8(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_default();
                        let ours = list.get(tc_idx).copied().unwrap_or(None);
                        match ours {
                            Some(v) => *csq_val == v.to_string(),
                            None => csq_val.is_empty(),
                        }
                    }
                    Cmp::ListI64 => {
                        let list = get_list_i64(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_default();
                        let ours = list.get(tc_idx).copied().unwrap_or(None);
                        match ours {
                            Some(v) => *csq_val == v.to_string(),
                            None => csq_val.is_empty(),
                        }
                    }
                    Cmp::Scalar => {
                        // Scalar Utf8 — same for all transcripts, just check once
                        let ours = get_str(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_default();
                        ours == *csq_val || (ours.is_empty() && csq_val.is_empty())
                    }
                    Cmp::ScalarFloat => {
                        let ours = arrow_batch
                            .column(arrow_idx)
                            .as_any()
                            .downcast_ref::<Float32Array>()
                            .and_then(|a| {
                                if a.is_null(row) {
                                    None
                                } else {
                                    Some(a.value(row))
                                }
                            });
                        match (ours, csq_val.parse::<f32>().ok()) {
                            (Some(a), Some(b)) => (a - b).abs() < 1e-5,
                            (None, None) => true,
                            (None, _) if csq_val.is_empty() => true,
                            _ => false,
                        }
                    }
                    Cmp::ScalarStr => {
                        let ours = get_str(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_default();
                        ours == *csq_val || (ours.is_empty() && csq_val.is_empty())
                    }
                    Cmp::CsvList => {
                        let ours = get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_default();
                        let joined = ours.join(",");
                        joined == *csq_val || (joined.is_empty() && csq_val.is_empty())
                    }
                };

                if is_match {
                    *col_match.entry(col_name).or_default() += 1;
                } else {
                    let cnt = col_mismatch.entry(col_name).or_default();
                    *cnt += 1;
                    let samples = col_samples.entry(col_name).or_default();
                    if samples.len() < 5 {
                        let ours_dbg = match cmp {
                            Cmp::ListStr => {
                                let list =
                                    get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                                        .unwrap_or_default();
                                list.get(tc_idx)
                                    .cloned()
                                    .unwrap_or_else(|| "OUT_OF_BOUNDS".to_string())
                            }
                            Cmp::ListI8 => {
                                let list = get_list_i8(arrow_batch.column(arrow_idx).as_ref(), row)
                                    .unwrap_or_default();
                                list.get(tc_idx)
                                    .map(|v| {
                                        v.map(|x| x.to_string())
                                            .unwrap_or_else(|| "NULL".to_string())
                                    })
                                    .unwrap_or_else(|| "OUT_OF_BOUNDS".to_string())
                            }
                            Cmp::ListI64 => {
                                let list =
                                    get_list_i64(arrow_batch.column(arrow_idx).as_ref(), row)
                                        .unwrap_or_default();
                                list.get(tc_idx)
                                    .map(|v| {
                                        v.map(|x| x.to_string())
                                            .unwrap_or_else(|| "NULL".to_string())
                                    })
                                    .unwrap_or_else(|| "OUT_OF_BOUNDS".to_string())
                            }
                            Cmp::Scalar | Cmp::ScalarStr => {
                                get_str(arrow_batch.column(arrow_idx).as_ref(), row)
                                    .unwrap_or_else(|| "NULL".to_string())
                            }
                            Cmp::ScalarFloat => arrow_batch
                                .column(arrow_idx)
                                .as_any()
                                .downcast_ref::<Float32Array>()
                                .and_then(|a| {
                                    if a.is_null(row) {
                                        None
                                    } else {
                                        Some(a.value(row).to_string())
                                    }
                                })
                                .unwrap_or_else(|| "NULL".to_string()),
                            Cmp::CsvList => {
                                get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                                    .map(|v| v.join(","))
                                    .unwrap_or_else(|| "NULL".to_string())
                            }
                        };
                        samples.push(format!(
                            "row {row} tc {tc_idx}: typed='{}' csq='{csq_val}'",
                            ours_dbg
                        ));
                    }
                }
            }
        }
    }

    // ── 4. Report ────────────────────────────────────────────────────
    eprintln!("\n=== Typed Arrow columns ↔ own CSQ (all transcripts, 1000 variants) ===");
    eprintln!("Total transcript entries checked: {total_entries_checked}");

    let mut failures: Vec<String> = Vec::new();
    for &(_, _, _, col_name) in &resolved {
        let matches = col_match.get(col_name).copied().unwrap_or(0);
        let mismatches = col_mismatch.get(col_name).copied().unwrap_or(0);
        let total = matches + mismatches;
        if mismatches > 0 {
            let pct = 100.0 * matches as f64 / total as f64;
            eprintln!(
                "  {col_name}: {matches}/{total} ({pct:.1}%)  *** {mismatches} mismatches ***"
            );
            if let Some(samples) = col_samples.get(col_name) {
                for s in samples {
                    eprintln!("    {s}");
                }
            }
            failures.push(format!("{col_name}: {mismatches} mismatches"));
        } else if total > 0 {
            eprintln!("  {col_name}: {matches}/{total} (100.0%)");
        }
    }

    assert!(
        total_entries_checked > 10000,
        "Should check many transcript entries (got {total_entries_checked})"
    );
    assert!(
        failures.is_empty(),
        "Typed columns inconsistent with own CSQ: {:?}",
        failures
    );

    eprintln!(
        "\nArrow typed columns test PASSED: all {} columns at 100% across {total_entries_checked} transcript entries",
        resolved.len()
    );
}
