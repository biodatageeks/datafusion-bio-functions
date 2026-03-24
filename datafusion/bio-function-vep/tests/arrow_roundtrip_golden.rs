//! Arrow output correctness test: annotate 1000 variants via `annotate_vep()`
//! table function → collect typed Arrow columns → verify that every top-level
//! typed column (SYMBOL, Gene, IMPACT, Consequence, Feature, AF, STRAND, etc.)
//! matches the corresponding pipe-delimited field in the variant's own CSQ
//! string for the transcript selected by the most-severe pick.
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

fn get_i64(col: &dyn Array, row: usize) -> Option<i64> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<Int64Array>() {
        return Some(a.value(row));
    }
    if let Some(a) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::UInt32Array>()
    {
        return Some(a.value(row) as i64);
    }
    if let Some(a) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int32Array>()
    {
        return Some(a.value(row) as i64);
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
        if let Some(a) = values.as_any().downcast_ref::<StringArray>() {
            out.push(if a.is_null(i) {
                String::new()
            } else {
                a.value(i).to_string()
            });
        } else if let Some(a) = values.as_any().downcast_ref::<StringViewArray>() {
            out.push(if a.is_null(i) {
                String::new()
            } else {
                a.value(i).to_string()
            });
        }
    }
    Some(out)
}

fn get_f32(col: &dyn Array, row: usize) -> Option<f32> {
    if col.is_null(row) {
        return None;
    }
    col.as_any()
        .downcast_ref::<Float32Array>()
        .map(|a| a.value(row))
}

fn get_i8(col: &dyn Array, row: usize) -> Option<i8> {
    if col.is_null(row) {
        return None;
    }
    col.as_any()
        .downcast_ref::<Int8Array>()
        .map(|a| a.value(row))
}

// ── CSQ field index ──────────────────────────────────────────────────

/// CSQ field ordering for VEP 115 `--everything` mode (80 pipe-delimited fields).
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

struct CsqEntry<'a> {
    fields: Vec<&'a str>,
}

impl<'a> CsqEntry<'a> {
    fn get(&self, idx: usize) -> &str {
        self.fields.get(idx).unwrap_or(&"")
    }
}

fn parse_csq(csq: &str) -> Vec<CsqEntry<'_>> {
    csq.split(',')
        .map(|entry| CsqEntry {
            fields: entry.split('|').collect(),
        })
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

    // ── 2. Verify expected typed columns exist ───────────────────────
    let expected_cols = [
        "CSQ",
        "most_severe_consequence",
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",
        "BIOTYPE",
        "STRAND",
        "VARIANT_CLASS",
        "AF",
        "MAX_AF",
        "CLIN_SIG",
    ];
    for col_name in &expected_cols {
        assert!(
            schema.index_of(col_name).is_ok(),
            "Expected typed column '{col_name}' missing from annotate_vep output"
        );
    }

    // ── 3. Typed columns ↔ own CSQ: 100% consistency ─────────────────
    //
    // For each row, parse the CSQ string emitted by the same annotate_vep()
    // call, locate the entry matching the typed Allele + Feature (the
    // most-severe pick), and verify every typed column matches the
    // corresponding CSQ field exactly. Same code path, same data — must
    // be 100%.

    #[derive(Clone, Copy)]
    enum Cmp {
        Str,
        TermList, // "&"-separated in CSQ ↔ List<Utf8> in Arrow (Consequence)
        AmpList,  // "&"-separated in CSQ ↔ List<Utf8> in Arrow (DOMAINS, TRANSCRIPTION_FACTORS)
        CsvList,  // raw string in CSQ ↔ List<Utf8> built by comma-split in Arrow
        Float,
        StrandI8,
        Int8Str,
        Int64Str,
    }

    // All 80 CSQ fields mapped to their typed Arrow column.
    // Column types: Str (Utf8), TermList (List<Utf8> with & separator in CSQ),
    // CsvList (List<Utf8> with , separator in CSQ), Float (Float32),
    // StrandI8 (Int8 ↔ "1"/"-1"), Int64Str (Int64 ↔ string), Int8Str (Int8 ↔ string).
    let columns_under_test: &[(&str, &str, Cmp)] = &[
        // ── Transcript-level (42) ──
        ("Allele", "Allele", Cmp::Str),
        ("Consequence", "Consequence", Cmp::TermList),
        ("IMPACT", "IMPACT", Cmp::Str),
        ("SYMBOL", "SYMBOL", Cmp::Str),
        ("Gene", "Gene", Cmp::Str),
        ("Feature_type", "Feature_type", Cmp::Str),
        ("Feature", "Feature", Cmp::Str),
        ("BIOTYPE", "BIOTYPE", Cmp::Str),
        ("EXON", "EXON", Cmp::Str),
        ("INTRON", "INTRON", Cmp::Str),
        ("HGVSc", "HGVSc", Cmp::Str),
        ("HGVSp", "HGVSp", Cmp::Str),
        ("cDNA_position", "cDNA_position", Cmp::Str),
        ("CDS_position", "CDS_position", Cmp::Str),
        ("Protein_position", "Protein_position", Cmp::Str),
        ("Amino_acids", "Amino_acids", Cmp::Str),
        ("Codons", "Codons", Cmp::Str),
        ("Existing_variation", "Existing_variation", Cmp::CsvList),
        ("DISTANCE", "DISTANCE", Cmp::Int64Str),
        ("STRAND", "STRAND", Cmp::StrandI8),
        ("FLAGS", "FLAGS", Cmp::Str),
        ("VARIANT_CLASS", "VARIANT_CLASS", Cmp::Str),
        ("SYMBOL_SOURCE", "SYMBOL_SOURCE", Cmp::Str),
        ("HGNC_ID", "HGNC_ID", Cmp::Str),
        ("CANONICAL", "CANONICAL", Cmp::Str),
        ("MANE", "MANE", Cmp::Str),
        ("MANE_SELECT", "MANE_SELECT", Cmp::Str),
        ("MANE_PLUS_CLINICAL", "MANE_PLUS_CLINICAL", Cmp::Str),
        ("TSL", "TSL", Cmp::Int8Str),
        ("APPRIS", "APPRIS", Cmp::Str),
        ("CCDS", "CCDS", Cmp::Str),
        ("ENSP", "ENSP", Cmp::Str),
        ("SWISSPROT", "SWISSPROT", Cmp::Str),
        ("TREMBL", "TREMBL", Cmp::Str),
        ("UNIPARC", "UNIPARC", Cmp::Str),
        ("UNIPROT_ISOFORM", "UNIPROT_ISOFORM", Cmp::Str),
        ("GENE_PHENO", "GENE_PHENO", Cmp::Str),
        ("SIFT", "SIFT", Cmp::Str),
        ("PolyPhen", "PolyPhen", Cmp::Str),
        ("DOMAINS", "DOMAINS", Cmp::AmpList),
        ("miRNA", "miRNA", Cmp::Str),
        ("HGVS_OFFSET", "HGVS_OFFSET", Cmp::Int64Str),
        // ── Frequency (29) ──
        ("AF", "AF", Cmp::Float),
        ("AFR_AF", "AFR_AF", Cmp::Float),
        ("AMR_AF", "AMR_AF", Cmp::Float),
        ("EAS_AF", "EAS_AF", Cmp::Float),
        ("EUR_AF", "EUR_AF", Cmp::Float),
        ("SAS_AF", "SAS_AF", Cmp::Float),
        ("gnomADe_AF", "gnomADe_AF", Cmp::Float),
        ("gnomADe_AFR_AF", "gnomADe_AFR_AF", Cmp::Float),
        ("gnomADe_AMR_AF", "gnomADe_AMR_AF", Cmp::Float),
        ("gnomADe_ASJ_AF", "gnomADe_ASJ_AF", Cmp::Float),
        ("gnomADe_EAS_AF", "gnomADe_EAS_AF", Cmp::Float),
        ("gnomADe_FIN_AF", "gnomADe_FIN_AF", Cmp::Float),
        ("gnomADe_MID_AF", "gnomADe_MID_AF", Cmp::Float),
        ("gnomADe_NFE_AF", "gnomADe_NFE_AF", Cmp::Float),
        ("gnomADe_REMAINING_AF", "gnomADe_REMAINING_AF", Cmp::Float),
        ("gnomADe_SAS_AF", "gnomADe_SAS_AF", Cmp::Float),
        ("gnomADg_AF", "gnomADg_AF", Cmp::Float),
        ("gnomADg_AFR_AF", "gnomADg_AFR_AF", Cmp::Float),
        ("gnomADg_AMI_AF", "gnomADg_AMI_AF", Cmp::Float),
        ("gnomADg_AMR_AF", "gnomADg_AMR_AF", Cmp::Float),
        ("gnomADg_ASJ_AF", "gnomADg_ASJ_AF", Cmp::Float),
        ("gnomADg_EAS_AF", "gnomADg_EAS_AF", Cmp::Float),
        ("gnomADg_FIN_AF", "gnomADg_FIN_AF", Cmp::Float),
        ("gnomADg_MID_AF", "gnomADg_MID_AF", Cmp::Float),
        ("gnomADg_NFE_AF", "gnomADg_NFE_AF", Cmp::Float),
        ("gnomADg_REMAINING_AF", "gnomADg_REMAINING_AF", Cmp::Float),
        ("gnomADg_SAS_AF", "gnomADg_SAS_AF", Cmp::Float),
        ("MAX_AF", "MAX_AF", Cmp::Float),
        ("MAX_AF_POPS", "MAX_AF_POPS", Cmp::Str),
        // ── Variant-level (9) ──
        ("CLIN_SIG", "CLIN_SIG", Cmp::CsvList),
        ("SOMATIC", "SOMATIC", Cmp::Str),
        ("PHENO", "PHENO", Cmp::Str),
        ("PUBMED", "PUBMED", Cmp::CsvList),
        ("MOTIF_NAME", "MOTIF_NAME", Cmp::Str),
        ("MOTIF_POS", "MOTIF_POS", Cmp::Str),
        ("HIGH_INF_POS", "HIGH_INF_POS", Cmp::Str),
        ("MOTIF_SCORE_CHANGE", "MOTIF_SCORE_CHANGE", Cmp::Float),
        (
            "TRANSCRIPTION_FACTORS",
            "TRANSCRIPTION_FACTORS",
            Cmp::AmpList,
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

    let csq_allele_fi = csq_field_idx("Allele");
    let csq_feature_fi = csq_field_idx("Feature");
    let our_csq_idx = schema.index_of("CSQ").unwrap();

    let mut col_match: HashMap<&str, usize> = HashMap::new();
    let mut col_mismatch: HashMap<&str, usize> = HashMap::new();
    let mut col_samples: HashMap<&str, Vec<String>> = HashMap::new();
    let mut matched_rows = 0usize;

    for row in 0..arrow_batch.num_rows() {
        let our_csq = match get_str(arrow_batch.column(our_csq_idx).as_ref(), row) {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let entries = parse_csq(&our_csq);
        if entries.is_empty() {
            continue;
        }

        let our_allele = schema
            .index_of("Allele")
            .ok()
            .and_then(|i| get_str(arrow_batch.column(i).as_ref(), row));
        let our_feature = schema
            .index_of("Feature")
            .ok()
            .and_then(|i| get_str(arrow_batch.column(i).as_ref(), row));

        let matched_entry = entries.iter().find(|e| {
            let allele_ok = our_allele
                .as_deref()
                .is_none_or(|a| e.get(csq_allele_fi) == a);
            let feature_ok = our_feature
                .as_deref()
                .is_none_or(|f| e.get(csq_feature_fi) == f);
            allele_ok && feature_ok
        });

        let Some(entry) = matched_entry else {
            continue;
        };

        matched_rows += 1;

        for &(arrow_idx, csq_fi, cmp, col_name) in &resolved {
            let csq_val = entry.get(csq_fi);
            let is_match = match cmp {
                Cmp::Str => {
                    let ours =
                        get_str(arrow_batch.column(arrow_idx).as_ref(), row).unwrap_or_default();
                    ours == csq_val || (ours.is_empty() && csq_val.is_empty())
                }
                Cmp::TermList | Cmp::AmpList => {
                    let ours = get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                        .unwrap_or_default();
                    let mut ours_sorted = ours.clone();
                    ours_sorted.sort();
                    let mut csq_terms: Vec<String> =
                        csq_val.split('&').map(|s| s.to_string()).collect();
                    csq_terms.sort();
                    ours_sorted == csq_terms
                        || (ours_sorted.is_empty()
                            && csq_terms.len() == 1
                            && csq_terms[0].is_empty())
                }
                Cmp::CsvList => {
                    let ours = get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                        .unwrap_or_default();
                    let joined = ours.join(",");
                    joined == csq_val || (joined.is_empty() && csq_val.is_empty())
                }
                Cmp::Float => {
                    let ours = get_f32(arrow_batch.column(arrow_idx).as_ref(), row);
                    match (ours, csq_val.parse::<f32>().ok()) {
                        (Some(a), Some(b)) => (a - b).abs() < 1e-5,
                        (None, None) => true,
                        (None, _) if csq_val.is_empty() => true,
                        _ => false,
                    }
                }
                Cmp::StrandI8 => {
                    let ours = get_i8(arrow_batch.column(arrow_idx).as_ref(), row);
                    match ours {
                        Some(s) => csq_val.parse::<i8>().is_ok_and(|g| g == s),
                        None => csq_val.is_empty(),
                    }
                }
                Cmp::Int8Str => {
                    let ours = get_i8(arrow_batch.column(arrow_idx).as_ref(), row);
                    match ours {
                        Some(v) => csq_val == v.to_string(),
                        None => csq_val.is_empty(),
                    }
                }
                Cmp::Int64Str => {
                    let ours = get_i64(arrow_batch.column(arrow_idx).as_ref(), row);
                    match ours {
                        Some(v) => csq_val == v.to_string(),
                        None => csq_val.is_empty(),
                    }
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
                        Cmp::Str => get_str(arrow_batch.column(arrow_idx).as_ref(), row)
                            .unwrap_or_else(|| "NULL".to_string()),
                        Cmp::CsvList => get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                            .map(|v| v.join(","))
                            .unwrap_or_else(|| "NULL".to_string()),
                        Cmp::TermList | Cmp::AmpList => {
                            get_list_str(arrow_batch.column(arrow_idx).as_ref(), row)
                                .map(|v| v.join("&"))
                                .unwrap_or_else(|| "NULL".to_string())
                        }
                        Cmp::Float => get_f32(arrow_batch.column(arrow_idx).as_ref(), row)
                            .map(|v| v.to_string())
                            .unwrap_or_else(|| "NULL".to_string()),
                        Cmp::StrandI8 | Cmp::Int8Str => {
                            get_i8(arrow_batch.column(arrow_idx).as_ref(), row)
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "NULL".to_string())
                        }
                        Cmp::Int64Str => get_i64(arrow_batch.column(arrow_idx).as_ref(), row)
                            .map(|v| v.to_string())
                            .unwrap_or_else(|| "NULL".to_string()),
                    };
                    samples.push(format!("row {row}: typed='{}' csq='{csq_val}'", ours_dbg));
                }
            }
        }
    }

    // ── 4. Report ────────────────────────────────────────────────────
    eprintln!("\n=== Typed Arrow columns ↔ own CSQ consistency (1000 variants) ===");
    eprintln!("Rows with matched typed+CSQ entry: {matched_rows}/1000");

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
        } else {
            eprintln!("  {col_name}: {matches}/{total} (100.0%)");
        }
    }

    assert!(
        matched_rows >= 900,
        "Should match typed Allele+Feature in own CSQ for ≥900 rows (got {matched_rows})"
    );
    assert!(
        failures.is_empty(),
        "Typed columns inconsistent with own CSQ: {:?}",
        failures
    );

    eprintln!(
        "\nArrow typed columns test PASSED: all {} columns at 100%",
        resolved.len()
    );
}
