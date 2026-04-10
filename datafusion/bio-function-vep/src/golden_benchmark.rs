//! Utilities for golden benchmark comparison against Ensembl VEP output.
//!
//! This module is used by the benchmark example that:
//! - samples a fixed number of variants from a gzipped VCF,
//! - runs Ensembl VEP in Docker as golden standard,
//! - compares the golden output with `annotate_vep()` output.

use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use datafusion::common::{DataFusionError, Result};
use flate2::read::MultiGzDecoder;

use crate::so_terms::SoTerm;

pub const DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ: &str =
    "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz";
pub const DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ_TBI: &str = "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz.tbi";
pub const DEFAULT_LOCAL_HG002_CHR22_VCF_GZ: &str = "vep-benchmark/data/HG002_chr22.vcf.gz";
pub const DEFAULT_LOCAL_HG002_CHR22_VCF_GZ_TBI: &str = "vep-benchmark/data/HG002_chr22.vcf.gz.tbi";
pub const DEFAULT_EXTERNAL_VEP_CACHE_DIR: &str =
    "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/cache";

#[derive(Debug, Clone, Eq)]
pub struct VariantKey {
    pub chrom: String,
    pub pos: i64,
    pub ref_allele: String,
    pub alt_alleles: String,
}

impl PartialEq for VariantKey {
    fn eq(&self, other: &Self) -> bool {
        self.chrom == other.chrom
            && self.pos == other.pos
            && self.ref_allele == other.ref_allele
            && self.alt_alleles == other.alt_alleles
    }
}

impl Hash for VariantKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chrom.hash(state);
        self.pos.hash(state);
        self.ref_allele.hash(state);
        self.alt_alleles.hash(state);
    }
}

#[derive(Debug, Clone)]
pub struct VariantAnnotation {
    pub key: VariantKey,
    pub csq: Option<String>,
    pub most_severe_consequence: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComparisonReport {
    pub golden_rows: usize,
    pub ours_rows: usize,
    pub intersection_rows: usize,
    pub missing_in_ours: usize,
    pub extra_in_ours: usize,
    pub golden_with_csq: usize,
    pub ours_with_csq: usize,
    pub csq_exact_matches: usize,
    pub most_severe_exact_matches: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TermComparisonReport {
    pub comparable_rows: usize,
    pub golden_with_terms: usize,
    pub ours_with_terms: usize,
    pub term_set_exact_matches: usize,
    pub golden_term_subset_matches: usize,
}

/// Per-variant discrepancy between golden and our annotation.
#[derive(Debug, Clone)]
pub struct VariantDiscrepancy {
    pub key: VariantKey,
    pub golden_most_severe: Option<String>,
    pub ours_most_severe: Option<String>,
    pub golden_terms: BTreeSet<String>,
    pub ours_terms: BTreeSet<String>,
    pub most_severe_match: bool,
    pub term_set_match: bool,
}

/// Ensure benchmark input is available in this repository.
///
/// Copies source file and optional tabix index when destination is missing.
pub fn ensure_local_copy(
    source_vcf_gz: &Path,
    dest_vcf_gz: &Path,
    source_tbi: Option<&Path>,
    dest_tbi: Option<&Path>,
) -> Result<()> {
    if !source_vcf_gz.exists() {
        return Err(DataFusionError::Execution(format!(
            "source VCF gzip does not exist: {}",
            source_vcf_gz.display()
        )));
    }

    if let Some(parent) = dest_vcf_gz.parent() {
        fs::create_dir_all(parent).map_err(io_err)?;
    }
    if !dest_vcf_gz.exists() {
        fs::copy(source_vcf_gz, dest_vcf_gz).map_err(io_err)?;
    }

    if let (Some(src), Some(dst)) = (source_tbi, dest_tbi) {
        if src.exists() {
            if let Some(parent) = dst.parent() {
                fs::create_dir_all(parent).map_err(io_err)?;
            }
            if !dst.exists() {
                fs::copy(src, dst).map_err(io_err)?;
            }
        }
    }

    Ok(())
}

/// Sample first `limit` variant records from a gzipped VCF, preserving headers.
pub fn sample_gz_vcf_first_n(
    input_vcf_gz: &Path,
    output_vcf: &Path,
    limit: usize,
) -> Result<usize> {
    let file = File::open(input_vcf_gz).map_err(io_err)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let out_file = File::create(output_vcf).map_err(io_err)?;
    let mut writer = BufWriter::new(out_file);

    sample_vcf_lines(reader, &mut writer, limit)
}

fn sample_vcf_lines<R: BufRead, W: Write>(
    mut reader: R,
    writer: &mut W,
    limit: usize,
) -> Result<usize> {
    let mut line = String::new();
    let mut count = 0usize;

    loop {
        line.clear();
        let read = reader.read_line(&mut line).map_err(io_err)?;
        if read == 0 {
            break;
        }
        if line.starts_with('#') {
            writer.write_all(line.as_bytes()).map_err(io_err)?;
            continue;
        }
        if count < limit {
            writer.write_all(line.as_bytes()).map_err(io_err)?;
            count += 1;
        } else {
            break;
        }
    }
    writer.flush().map_err(io_err)?;
    Ok(count)
}

pub fn parse_vep_vcf_annotations(vcf_path: &Path) -> Result<Vec<VariantAnnotation>> {
    let file = File::open(vcf_path).map_err(io_err)?;
    let mut out = Vec::new();

    let is_gz = vcf_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        let reader = BufReader::new(MultiGzDecoder::new(file));
        parse_vep_vcf_lines(reader, &mut out)?;
    } else {
        let reader = BufReader::new(file);
        parse_vep_vcf_lines(reader, &mut out)?;
    }
    Ok(out)
}

fn parse_vep_vcf_lines<R: BufRead>(mut reader: R, out: &mut Vec<VariantAnnotation>) -> Result<()> {
    let mut line = String::new();
    loop {
        line.clear();
        let read = reader.read_line(&mut line).map_err(io_err)?;
        if read == 0 {
            break;
        }
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 8 {
            continue;
        }
        let pos = fields[1].parse::<i64>().map_err(|e| {
            DataFusionError::Execution(format!("failed to parse POS '{}': {e}", fields[1]))
        })?;
        let key = VariantKey {
            chrom: normalize_chrom(fields[0]),
            pos,
            ref_allele: fields[3].to_string(),
            alt_alleles: fields[4].to_string(),
        };
        let csq = parse_info_value(fields[7], "CSQ");
        let most = csq.as_deref().and_then(most_severe_from_csq);
        out.push(VariantAnnotation {
            key,
            csq,
            most_severe_consequence: most,
        });
    }
    Ok(())
}

pub fn parse_info_value(info: &str, key: &str) -> Option<String> {
    let prefix = format!("{key}=");
    for kv in info.split(';') {
        if kv.starts_with(&prefix) {
            return Some(kv[prefix.len()..].to_string());
        }
    }
    None
}

pub fn normalize_chrom(chrom: &str) -> String {
    chrom.strip_prefix("chr").unwrap_or(chrom).to_string()
}

pub fn most_severe_from_csq(csq: &str) -> Option<String> {
    let mut best: Option<(SoTerm, String)> = None;

    for ann in csq.split(',') {
        let mut parts = ann.split('|');
        let _allele = parts.next();
        let consequence_field = parts.next().unwrap_or("");
        for term in consequence_field.split('&') {
            let Some(parsed) = SoTerm::from_str(term) else {
                continue;
            };
            match &best {
                None => best = Some((parsed, term.to_string())),
                Some((best_term, _)) if parsed.rank() < best_term.rank() => {
                    best = Some((parsed, term.to_string()))
                }
                _ => {}
            }
        }
    }

    best.map(|(_, term)| term)
}

pub fn compare_annotations(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
) -> ComparisonReport {
    let golden_map: HashMap<VariantKey, &VariantAnnotation> =
        golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> =
        ours.iter().map(|r| (r.key.clone(), r)).collect();

    let mut intersection_rows = 0usize;
    let mut missing_in_ours = 0usize;
    let mut golden_with_csq = 0usize;
    let mut ours_with_csq = 0usize;
    let mut csq_exact_matches = 0usize;
    let mut most_severe_exact_matches = 0usize;

    for (key, golden_row) in &golden_map {
        if golden_row.csq.as_deref().is_some_and(|s| !s.is_empty()) {
            golden_with_csq += 1;
        }

        match ours_map.get(key) {
            Some(ours_row) => {
                intersection_rows += 1;
                if ours_row.csq.as_deref().is_some_and(|s| !s.is_empty()) {
                    ours_with_csq += 1;
                }
                if ours_row.csq == golden_row.csq {
                    csq_exact_matches += 1;
                }
                if ours_row.most_severe_consequence == golden_row.most_severe_consequence {
                    most_severe_exact_matches += 1;
                }
            }
            None => {
                missing_in_ours += 1;
            }
        }
    }

    let extra_in_ours = ours_map
        .keys()
        .filter(|k| !golden_map.contains_key(*k))
        .count();

    ComparisonReport {
        golden_rows: golden_map.len(),
        ours_rows: ours_map.len(),
        intersection_rows,
        missing_in_ours,
        extra_in_ours,
        golden_with_csq,
        ours_with_csq,
        csq_exact_matches,
        most_severe_exact_matches,
    }
}

pub fn extract_csq_term_set(csq: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for ann in csq.split(',') {
        let mut parts = ann.split('|');
        let _allele = parts.next();
        let term_field = parts.next().unwrap_or("");
        for term in term_field.split('&').filter(|t| !t.is_empty()) {
            out.insert(term.to_string());
        }
    }
    out
}

pub fn compare_annotation_terms(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
) -> TermComparisonReport {
    let golden_map: HashMap<VariantKey, &VariantAnnotation> =
        golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> =
        ours.iter().map(|r| (r.key.clone(), r)).collect();

    let mut comparable_rows = 0usize;
    let mut golden_with_terms = 0usize;
    let mut ours_with_terms = 0usize;
    let mut term_set_exact_matches = 0usize;
    let mut golden_term_subset_matches = 0usize;

    for (key, golden_row) in &golden_map {
        let Some(ours_row) = ours_map.get(key) else {
            continue;
        };
        comparable_rows += 1;

        let golden_terms = golden_row
            .csq
            .as_deref()
            .map(extract_csq_term_set)
            .unwrap_or_default();
        let ours_terms = ours_row
            .csq
            .as_deref()
            .map(extract_csq_term_set)
            .unwrap_or_default();

        if !golden_terms.is_empty() {
            golden_with_terms += 1;
        }
        if !ours_terms.is_empty() {
            ours_with_terms += 1;
        }
        if golden_terms == ours_terms {
            term_set_exact_matches += 1;
        }
        if golden_terms.is_subset(&ours_terms) {
            golden_term_subset_matches += 1;
        }
    }

    TermComparisonReport {
        comparable_rows,
        golden_with_terms,
        ours_with_terms,
        term_set_exact_matches,
        golden_term_subset_matches,
    }
}

/// Collect per-variant discrepancies between golden and our annotations.
///
/// Only returns entries where most_severe or term set differ.
pub fn collect_discrepancies(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
) -> Vec<VariantDiscrepancy> {
    let golden_map: HashMap<VariantKey, &VariantAnnotation> =
        golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> =
        ours.iter().map(|r| (r.key.clone(), r)).collect();

    let mut out = Vec::new();

    for (key, golden_row) in &golden_map {
        let Some(ours_row) = ours_map.get(key) else {
            // Missing in ours — always a discrepancy.
            let golden_terms = golden_row
                .csq
                .as_deref()
                .map(extract_csq_term_set)
                .unwrap_or_default();
            out.push(VariantDiscrepancy {
                key: key.clone(),
                golden_most_severe: golden_row.most_severe_consequence.clone(),
                ours_most_severe: None,
                golden_terms,
                ours_terms: BTreeSet::new(),
                most_severe_match: false,
                term_set_match: false,
            });
            continue;
        };

        let golden_terms = golden_row
            .csq
            .as_deref()
            .map(extract_csq_term_set)
            .unwrap_or_default();
        let ours_terms = ours_row
            .csq
            .as_deref()
            .map(extract_csq_term_set)
            .unwrap_or_default();

        let most_severe_match =
            golden_row.most_severe_consequence == ours_row.most_severe_consequence;
        let term_set_match = golden_terms == ours_terms;

        if !most_severe_match || !term_set_match {
            out.push(VariantDiscrepancy {
                key: key.clone(),
                golden_most_severe: golden_row.most_severe_consequence.clone(),
                ours_most_severe: ours_row.most_severe_consequence.clone(),
                golden_terms,
                ours_terms,
                most_severe_match,
                term_set_match,
            });
        }
    }

    // Sort by position for deterministic output.
    out.sort_by(|a, b| {
        a.key
            .chrom
            .cmp(&b.key.chrom)
            .then(a.key.pos.cmp(&b.key.pos))
    });
    out
}

/// CSQ field names matching VEP `--fields` output order (74 fields).
pub const CSQ_FIELD_NAMES: &[&str] = &[
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
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
    "SOURCE",
    // Batch 1 fields.
    "VARIANT_CLASS",
    "CANONICAL",
    "TSL",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "ENSP",
    "GENE_PHENO",
    "CCDS",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    // Batch 3 fields.
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "gnomADe_AF",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
    "gnomADg_AF",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
];

/// CSQ field names for `--everything` mode (80 fields).
///
/// Traceability:
/// - VEP Constants.pm CSQ field order
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Constants.pm#L66-L138>
///
/// Key differences from 74-field layout:
/// - `SOURCE` removed (only present in `--merged` mode)
/// - `VARIANT_CLASS` moved from after SOURCE to after FLAGS
/// - `MANE` generic field added (separate from MANE_SELECT/MANE_PLUS_CLINICAL)
/// - `APPRIS`, `SIFT`, `PolyPhen`, `DOMAINS`, `miRNA`, `HGVS_OFFSET` added
/// - gnomAD sub-population fields have `_AF` suffix
/// - `MOTIF_*` fields moved to end (positions 75-79)
pub const CSQ_FIELD_NAMES_EVERYTHING: &[&str] = &[
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
    // Batch 3 fields (gnomAD sub-pops have _AF suffix in --everything).
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
    // Motif fields moved to end in --everything.
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
];

/// Return the CSQ field names for the selected transcript source mode.
///
/// The non-`everything` profile preserves the existing explicit benchmark field
/// order and swaps the transcript-source block near `SOURCE`:
/// - default Ensembl: `SOURCE`
/// - RefSeq: `REFSEQ_MATCH`, `BAM_EDIT`
/// - merged: `REFSEQ_MATCH`, `BAM_EDIT`, `SOURCE`
///
/// The `everything` profile inserts the RefSeq-specific fields next to the
/// transcript metadata block, after `UNIPROT_ISOFORM`, matching VEP's flag
/// expansion order more closely.
pub fn csq_field_names_for_mode(everything: bool, refseq: bool, merged: bool) -> Vec<&'static str> {
    let mut fields = if everything {
        CSQ_FIELD_NAMES_EVERYTHING.to_vec()
    } else {
        CSQ_FIELD_NAMES.to_vec()
    };

    if everything {
        if let Some(insert_at) = fields.iter().position(|field| *field == "GENE_PHENO") {
            if merged {
                fields.splice(insert_at..insert_at, ["REFSEQ_MATCH", "BAM_EDIT", "SOURCE"]);
            } else if refseq {
                fields.splice(insert_at..insert_at, ["REFSEQ_MATCH", "BAM_EDIT"]);
            }
        }
        return fields;
    }

    if let Some(source_idx) = fields.iter().position(|field| *field == "SOURCE") {
        if merged {
            fields.splice(
                source_idx..=source_idx,
                ["REFSEQ_MATCH", "BAM_EDIT", "SOURCE"],
            );
        } else if refseq {
            fields.splice(source_idx..=source_idx, ["REFSEQ_MATCH", "BAM_EDIT"]);
        }
    }

    fields
}

/// Sample of a field-level mismatch for debugging.
#[derive(Debug, Clone)]
pub struct FieldMismatchSample {
    pub variant_key: VariantKey,
    pub allele: String,
    pub feature: String,
    pub golden_val: String,
    pub ours_val: String,
}

/// Per-field match report for CSQ comparison.
#[derive(Debug, Clone)]
pub struct CsqFieldReport {
    pub field_names: Vec<String>,
    pub field_match_counts: Vec<usize>,
    pub total_entries_compared: usize,
    /// Up to 5 mismatch samples per field index, for debugging.
    pub field_mismatch_samples: HashMap<usize, Vec<FieldMismatchSample>>,
}

/// Diagnostic report for unmatched CSQ entries between golden and ours.
#[derive(Debug, Clone)]
pub struct CsqUnmatchedReport {
    /// Golden entries not found in ours, grouped by Feature_type.
    pub golden_only_by_ft: HashMap<String, usize>,
    /// Our entries not found in golden, grouped by Feature_type.
    pub ours_only_by_ft: HashMap<String, usize>,
    /// Total golden CSQ entries.
    pub golden_total: usize,
    /// Total ours CSQ entries.
    pub ours_total: usize,
    /// Matched by (Allele, Feature) key.
    pub matched: usize,
    /// Sample of unmatched golden Feature IDs (first 30).
    pub golden_only_sample: Vec<(String, String, String)>,
    /// Sample of unmatched ours Feature IDs (first 30).
    pub ours_only_sample: Vec<(String, String, String)>,
}

/// Diagnostic sample for `(Allele, Feature)` multiplicity drift within one variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CsqMultiplicitySample {
    pub variant_key: VariantKey,
    pub allele: String,
    pub feature: String,
    pub feature_type: String,
    pub golden_count: usize,
    pub ours_count: usize,
}

/// Report of CSQ entry multiplicity drift after `(Allele, Feature)` matching.
#[derive(Debug, Clone, Default)]
pub struct CsqMultiplicityReport {
    pub golden_extra_total: usize,
    pub ours_extra_total: usize,
    pub samples: Vec<CsqMultiplicitySample>,
}

/// Diagnose why CSQ entries are unmatched between golden and ours.
pub fn diagnose_unmatched_csq(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
) -> CsqUnmatchedReport {
    let golden_map: HashMap<VariantKey, &VariantAnnotation> =
        golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> =
        ours.iter().map(|r| (r.key.clone(), r)).collect();

    let mut golden_only_by_ft: HashMap<String, usize> = HashMap::new();
    let mut ours_only_by_ft: HashMap<String, usize> = HashMap::new();
    let mut golden_total = 0usize;
    let mut ours_total = 0usize;
    let mut matched = 0usize;
    let mut golden_only_sample: Vec<(String, String, String)> = Vec::new();
    let mut ours_only_sample: Vec<(String, String, String)> = Vec::new();

    for (key, golden_row) in &golden_map {
        let golden_csq = match golden_row.csq.as_deref() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let golden_entries = parse_csq_entries(golden_csq);
        golden_total += golden_entries.len();

        let ours_csq = ours_map
            .get(key)
            .and_then(|r| r.csq.as_deref())
            .unwrap_or("");
        let ours_entries = parse_csq_entries(ours_csq);

        let mut ours_keys: HashMap<(String, String), &CsqEntry> = HashMap::new();
        for e in &ours_entries {
            ours_keys
                .entry((e.allele().to_string(), e.feature().to_string()))
                .or_insert(e);
        }

        for ge in &golden_entries {
            let k = (ge.allele().to_string(), ge.feature().to_string());
            if ours_keys.contains_key(&k) {
                matched += 1;
            } else {
                let ft = ge.fields.get(5).cloned().unwrap_or_default();
                *golden_only_by_ft.entry(ft.clone()).or_default() += 1;
                if golden_only_sample.len() < 30 {
                    golden_only_sample.push((
                        ge.allele().to_string(),
                        ge.feature().to_string(),
                        ft,
                    ));
                }
            }
        }
    }

    // Count ours-only entries
    for (key, ours_row) in &ours_map {
        let ours_csq = match ours_row.csq.as_deref() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let ours_entries = parse_csq_entries(ours_csq);
        ours_total += ours_entries.len();

        let golden_csq = golden_map
            .get(key)
            .and_then(|r| r.csq.as_deref())
            .unwrap_or("");
        let golden_entries = parse_csq_entries(golden_csq);

        let mut golden_keys: HashMap<(String, String), &CsqEntry> = HashMap::new();
        for e in &golden_entries {
            golden_keys
                .entry((e.allele().to_string(), e.feature().to_string()))
                .or_insert(e);
        }

        for oe in &ours_entries {
            let k = (oe.allele().to_string(), oe.feature().to_string());
            if !golden_keys.contains_key(&k) {
                let ft = oe.fields.get(5).cloned().unwrap_or_default();
                *ours_only_by_ft.entry(ft.clone()).or_default() += 1;
                if ours_only_sample.len() < 30 {
                    ours_only_sample.push((oe.allele().to_string(), oe.feature().to_string(), ft));
                }
            }
        }
    }

    CsqUnmatchedReport {
        golden_only_by_ft,
        ours_only_by_ft,
        golden_total,
        ours_total,
        matched,
        golden_only_sample,
        ours_only_sample,
    }
}

/// Diagnose multiplicity drift for entries that share the same
/// `(Allele, Feature)` key within a single variant.
pub fn diagnose_csq_multiplicity(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
) -> CsqMultiplicityReport {
    let golden_map: HashMap<VariantKey, &VariantAnnotation> =
        golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> =
        ours.iter().map(|r| (r.key.clone(), r)).collect();

    let mut report = CsqMultiplicityReport::default();

    for (key, golden_row) in &golden_map {
        let Some(ours_row) = ours_map.get(key) else {
            continue;
        };

        let golden_entries = parse_csq_entries(golden_row.csq.as_deref().unwrap_or(""));
        let ours_entries = parse_csq_entries(ours_row.csq.as_deref().unwrap_or(""));

        let mut golden_counts: HashMap<(String, String), usize> = HashMap::new();
        let mut ours_counts: HashMap<(String, String), usize> = HashMap::new();
        let mut feature_types: HashMap<(String, String), String> = HashMap::new();

        for entry in &golden_entries {
            let key = (entry.allele().to_string(), entry.feature().to_string());
            *golden_counts.entry(key.clone()).or_default() += 1;
            feature_types
                .entry(key)
                .or_insert_with(|| entry.fields.get(5).cloned().unwrap_or_default());
        }

        for entry in &ours_entries {
            let key = (entry.allele().to_string(), entry.feature().to_string());
            *ours_counts.entry(key.clone()).or_default() += 1;
            feature_types
                .entry(key)
                .or_insert_with(|| entry.fields.get(5).cloned().unwrap_or_default());
        }

        let mut all_keys: BTreeSet<(String, String)> = BTreeSet::new();
        all_keys.extend(golden_counts.keys().cloned());
        all_keys.extend(ours_counts.keys().cloned());

        for pair in all_keys {
            let golden_count = *golden_counts.get(&pair).unwrap_or(&0);
            let ours_count = *ours_counts.get(&pair).unwrap_or(&0);
            if golden_count == ours_count {
                continue;
            }

            if ours_count > golden_count {
                report.ours_extra_total += ours_count - golden_count;
            } else {
                report.golden_extra_total += golden_count - ours_count;
            }

            if report.samples.len() < 30 {
                report.samples.push(CsqMultiplicitySample {
                    variant_key: (*key).clone(),
                    allele: pair.0.clone(),
                    feature: pair.1.clone(),
                    feature_type: feature_types.get(&pair).cloned().unwrap_or_default(),
                    golden_count,
                    ours_count,
                });
            }
        }
    }

    report
}

/// Parsed CSQ entry (one pipe-delimited annotation from one transcript).
#[derive(Debug, Clone)]
struct CsqEntry {
    fields: Vec<String>,
}

impl CsqEntry {
    fn allele(&self) -> &str {
        self.fields.first().map(|s| s.as_str()).unwrap_or("")
    }
    fn feature(&self) -> &str {
        self.fields.get(6).map(|s| s.as_str()).unwrap_or("")
    }
}

fn parse_csq_entries(csq: &str) -> Vec<CsqEntry> {
    csq.split(',')
        .map(|ann| CsqEntry {
            fields: ann.split('|').map(|s| s.to_string()).collect(),
        })
        .collect()
}

/// Compare per-field CSQ accuracy between golden VEP output and ours.
///
/// Matches entries by (Allele, Feature) key within each variant, then
/// compares each CSQ field. Uses the provided `field_names` to determine
/// the number of fields (pass `CSQ_FIELD_NAMES` for 74-field mode or
/// `CSQ_FIELD_NAMES_EVERYTHING` for 80-field `--everything` mode).
pub fn compare_csq_fields(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
) -> CsqFieldReport {
    compare_csq_fields_with_names(golden, ours, CSQ_FIELD_NAMES)
}

/// Compare per-field CSQ accuracy using a custom field name list.
pub fn compare_csq_fields_with_names(
    golden: &[VariantAnnotation],
    ours: &[VariantAnnotation],
    field_names: &[&str],
) -> CsqFieldReport {
    let field_count = field_names.len();
    let mut field_match_counts = vec![0usize; field_count];
    let mut field_mismatch_samples: HashMap<usize, Vec<FieldMismatchSample>> = HashMap::new();
    let mut total_entries = 0usize;

    let golden_map: HashMap<VariantKey, &VariantAnnotation> =
        golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> =
        ours.iter().map(|r| (r.key.clone(), r)).collect();

    for (key, golden_row) in &golden_map {
        let Some(ours_row) = ours_map.get(key) else {
            continue;
        };

        let golden_csq = match golden_row.csq.as_deref() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let ours_csq = match ours_row.csq.as_deref() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };

        let golden_entries = parse_csq_entries(golden_csq);
        let ours_entries = parse_csq_entries(ours_csq);

        // Index ours entries by (allele, feature) for lookup.
        let mut ours_by_key: HashMap<(String, String), &CsqEntry> = HashMap::new();
        for entry in &ours_entries {
            let k = (entry.allele().to_string(), entry.feature().to_string());
            ours_by_key.entry(k).or_insert(entry);
        }

        for golden_entry in &golden_entries {
            let k = (
                golden_entry.allele().to_string(),
                golden_entry.feature().to_string(),
            );
            let Some(ours_entry) = ours_by_key.get(&k) else {
                // No matching entry in ours — count as mismatch for all fields.
                total_entries += 1;
                continue;
            };

            total_entries += 1;
            for i in 0..field_count {
                let golden_val = golden_entry.fields.get(i).map(|s| s.as_str()).unwrap_or("");
                let ours_val = ours_entry.fields.get(i).map(|s| s.as_str()).unwrap_or("");
                if golden_val == ours_val {
                    field_match_counts[i] += 1;
                } else {
                    let samples = field_mismatch_samples.entry(i).or_insert_with(Vec::new);
                    if samples.len() < 10 {
                        samples.push(FieldMismatchSample {
                            variant_key: key.clone(),
                            allele: golden_entry.allele().to_string(),
                            feature: golden_entry.feature().to_string(),
                            golden_val: golden_val.to_string(),
                            ours_val: ours_val.to_string(),
                        });
                    }
                }
            }
        }
    }

    CsqFieldReport {
        field_names: field_names.iter().map(|s| s.to_string()).collect(),
        field_match_counts,
        total_entries_compared: total_entries,
        field_mismatch_samples,
    }
}

fn io_err<E: std::fmt::Display>(e: E) -> DataFusionError {
    DataFusionError::Execution(e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn sample_lines_preserves_header_and_limit() {
        let input = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\t.\tA\tG\t.\t.\t.
1\t20\t.\tC\tT\t.\t.\t.
1\t30\t.\tG\tA\t.\t.\t.
";
        let mut out = Vec::new();
        let n = sample_vcf_lines(Cursor::new(input.as_bytes()), &mut out, 2).unwrap();
        assert_eq!(n, 2);
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("##fileformat=VCFv4.2"));
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
        assert!(text.contains("1\t10\t.\tA\tG"));
        assert!(text.contains("1\t20\t.\tC\tT"));
        assert!(!text.contains("1\t30\t.\tG\tA"));
    }

    #[test]
    fn parse_info_value_extracts_key() {
        let info = "AC=1;AF=0.5;CSQ=A|missense_variant|MODERATE;DB=1";
        assert_eq!(
            parse_info_value(info, "CSQ"),
            Some("A|missense_variant|MODERATE".to_string())
        );
        assert_eq!(parse_info_value(info, "NOTHING"), None);
    }

    #[test]
    fn normalize_chrom_strips_chr_prefix() {
        assert_eq!(normalize_chrom("chr22"), "22");
        assert_eq!(normalize_chrom("22"), "22");
        assert_eq!(normalize_chrom("chrX"), "X");
    }

    #[test]
    fn most_severe_from_csq_picks_best_rank() {
        let csq = "A|synonymous_variant|LOW,B|stop_gained|HIGH";
        assert_eq!(most_severe_from_csq(csq), Some("stop_gained".to_string()));
    }

    #[test]
    fn most_severe_from_csq_handles_ampersand_terms() {
        let csq = "A|splice_region_variant&intron_variant|LOW";
        assert_eq!(
            most_severe_from_csq(csq),
            Some("splice_region_variant".to_string())
        );
    }

    #[test]
    fn parse_vep_vcf_lines_extracts_key_and_csq() {
        let input = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|MODERATE
";
        let mut out = Vec::new();
        parse_vep_vcf_lines(Cursor::new(input.as_bytes()), &mut out).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].key.chrom, "22");
        assert_eq!(out[0].key.pos, 100);
        assert_eq!(out[0].key.ref_allele, "A");
        assert_eq!(out[0].key.alt_alleles, "G");
        assert_eq!(
            out[0].most_severe_consequence,
            Some("missense_variant".to_string())
        );
    }

    #[test]
    fn compare_annotations_counts_matches_and_mismatches() {
        let golden = vec![
            VariantAnnotation {
                key: VariantKey {
                    chrom: "22".to_string(),
                    pos: 100,
                    ref_allele: "A".to_string(),
                    alt_alleles: "G".to_string(),
                },
                csq: Some("G|missense_variant|MODERATE".to_string()),
                most_severe_consequence: Some("missense_variant".to_string()),
            },
            VariantAnnotation {
                key: VariantKey {
                    chrom: "22".to_string(),
                    pos: 200,
                    ref_allele: "C".to_string(),
                    alt_alleles: "T".to_string(),
                },
                csq: Some("T|synonymous_variant|LOW".to_string()),
                most_severe_consequence: Some("synonymous_variant".to_string()),
            },
        ];
        let ours = vec![
            VariantAnnotation {
                key: VariantKey {
                    chrom: "22".to_string(),
                    pos: 100,
                    ref_allele: "A".to_string(),
                    alt_alleles: "G".to_string(),
                },
                csq: Some("G|missense_variant|MODERATE".to_string()),
                most_severe_consequence: Some("missense_variant".to_string()),
            },
            VariantAnnotation {
                key: VariantKey {
                    chrom: "22".to_string(),
                    pos: 300,
                    ref_allele: "G".to_string(),
                    alt_alleles: "A".to_string(),
                },
                csq: None,
                most_severe_consequence: None,
            },
        ];

        let report = compare_annotations(&golden, &ours);
        assert_eq!(
            report,
            ComparisonReport {
                golden_rows: 2,
                ours_rows: 2,
                intersection_rows: 1,
                missing_in_ours: 1,
                extra_in_ours: 1,
                golden_with_csq: 2,
                ours_with_csq: 1,
                csq_exact_matches: 1,
                most_severe_exact_matches: 1,
            }
        );
    }

    #[test]
    fn extract_csq_term_set_collects_unique_terms() {
        let csq =
            "A|splice_region_variant&intron_variant|LOW,B|missense_variant&intron_variant|MODERATE";
        let terms = extract_csq_term_set(csq);
        assert_eq!(
            terms,
            BTreeSet::from([
                "intron_variant".to_string(),
                "missense_variant".to_string(),
                "splice_region_variant".to_string()
            ])
        );
    }

    #[test]
    fn compare_annotation_terms_reports_exact_and_subset_matches() {
        let golden = vec![VariantAnnotation {
            key: VariantKey {
                chrom: "1".to_string(),
                pos: 100,
                ref_allele: "A".to_string(),
                alt_alleles: "G".to_string(),
            },
            csq: Some("G|splice_region_variant&intron_variant|LOW".to_string()),
            most_severe_consequence: Some("splice_region_variant".to_string()),
        }];
        let ours = vec![VariantAnnotation {
            key: VariantKey {
                chrom: "1".to_string(),
                pos: 100,
                ref_allele: "A".to_string(),
                alt_alleles: "G".to_string(),
            },
            csq: Some(
                "G|splice_region_variant&intron_variant&coding_transcript_variant|LOW".to_string(),
            ),
            most_severe_consequence: Some("splice_region_variant".to_string()),
        }];

        let report = compare_annotation_terms(&golden, &ours);
        assert_eq!(
            report,
            TermComparisonReport {
                comparable_rows: 1,
                golden_with_terms: 1,
                ours_with_terms: 1,
                term_set_exact_matches: 0,
                golden_term_subset_matches: 1,
            }
        );
    }

    #[test]
    fn csq_field_names_has_74_entries() {
        assert_eq!(CSQ_FIELD_NAMES.len(), 74);
        assert_eq!(CSQ_FIELD_NAMES[0], "Allele");
        assert_eq!(CSQ_FIELD_NAMES[1], "Consequence");
        assert_eq!(CSQ_FIELD_NAMES[6], "Feature");
        assert_eq!(CSQ_FIELD_NAMES[17], "Existing_variation");
        assert_eq!(CSQ_FIELD_NAMES[19], "STRAND");
        assert_eq!(CSQ_FIELD_NAMES[28], "SOURCE");
        assert_eq!(CSQ_FIELD_NAMES[29], "VARIANT_CLASS");
        assert_eq!(CSQ_FIELD_NAMES[40], "UNIPROT_ISOFORM");
        assert_eq!(CSQ_FIELD_NAMES[41], "AF");
        assert_eq!(CSQ_FIELD_NAMES[47], "gnomADe_AF");
        assert_eq!(CSQ_FIELD_NAMES[57], "gnomADg_AF");
        assert_eq!(CSQ_FIELD_NAMES[68], "MAX_AF");
        assert_eq!(CSQ_FIELD_NAMES[73], "PUBMED");
    }

    #[test]
    fn csq_field_names_everything_has_80_entries() {
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING.len(), 80);
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[0], "Allele");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[1], "Consequence");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[6], "Feature");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[17], "Existing_variation");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[19], "STRAND");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[20], "FLAGS");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[21], "VARIANT_CLASS");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[24], "CANONICAL");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[25], "MANE");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[29], "APPRIS");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[35], "UNIPROT_ISOFORM");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[37], "SIFT");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[38], "PolyPhen");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[39], "DOMAINS");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[40], "miRNA");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[41], "HGVS_OFFSET");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[42], "AF");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[49], "gnomADe_AFR_AF");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[58], "gnomADg_AF");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[69], "MAX_AF");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[74], "PUBMED");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[75], "MOTIF_NAME");
        assert_eq!(CSQ_FIELD_NAMES_EVERYTHING[79], "TRANSCRIPTION_FACTORS");
    }

    #[test]
    fn csq_field_names_for_refseq_and_merged_modes_insert_expected_fields() {
        let refseq = csq_field_names_for_mode(false, true, false);
        assert_eq!(refseq.len(), 75);
        assert_eq!(refseq[28], "REFSEQ_MATCH");
        assert_eq!(refseq[29], "BAM_EDIT");
        assert_eq!(refseq[30], "VARIANT_CLASS");

        let merged = csq_field_names_for_mode(false, false, true);
        assert_eq!(merged.len(), 76);
        assert_eq!(merged[28], "REFSEQ_MATCH");
        assert_eq!(merged[29], "BAM_EDIT");
        assert_eq!(merged[30], "SOURCE");
        assert_eq!(merged[31], "VARIANT_CLASS");

        let everything_refseq = csq_field_names_for_mode(true, true, false);
        assert_eq!(everything_refseq.len(), 82);
        assert_eq!(everything_refseq[36], "REFSEQ_MATCH");
        assert_eq!(everything_refseq[37], "BAM_EDIT");
        assert_eq!(everything_refseq[38], "GENE_PHENO");

        let everything_merged = csq_field_names_for_mode(true, false, true);
        assert_eq!(everything_merged.len(), 83);
        assert_eq!(everything_merged[36], "REFSEQ_MATCH");
        assert_eq!(everything_merged[37], "BAM_EDIT");
        assert_eq!(everything_merged[38], "SOURCE");
        assert_eq!(everything_merged[39], "GENE_PHENO");
    }

    #[test]
    fn parse_csq_entries_splits_correctly() {
        let csq = "G|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000269305|protein_coding|rs123|1|HGNC|11998|Ensembl";
        let entries = parse_csq_entries(csq);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].allele(), "G");
        assert_eq!(entries[0].feature(), "ENST00000269305");
        assert_eq!(entries[0].fields.len(), 13);
    }

    #[test]
    fn parse_csq_entries_handles_multi_transcript() {
        let csq = "G|missense_variant|MODERATE|SYM|GENE1|Transcript|TX1|pc|rs1|1|HGNC|100|Ensembl,G|intron_variant|MODIFIER|SYM|GENE1|Transcript|TX2|pc|rs1|-1|HGNC|100|Ensembl";
        let entries = parse_csq_entries(csq);
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].feature(), "TX1");
        assert_eq!(entries[1].feature(), "TX2");
    }

    #[test]
    fn compare_csq_fields_all_match() {
        let csq = "G|missense_variant|MODERATE|TP53|ENSG1|Transcript|ENST1|protein_coding|rs1|1|HGNC|11998|Ensembl";
        let golden = vec![VariantAnnotation {
            key: VariantKey {
                chrom: "22".into(),
                pos: 100,
                ref_allele: "A".into(),
                alt_alleles: "G".into(),
            },
            csq: Some(csq.to_string()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let ours = golden.clone();
        let report = compare_csq_fields(&golden, &ours);
        assert_eq!(report.total_entries_compared, 1);
        assert_eq!(report.field_match_counts.len(), 74);
        for &c in &report.field_match_counts {
            assert_eq!(c, 1, "all fields should match");
        }
    }

    #[test]
    fn compare_csq_fields_partial_mismatch() {
        let golden_csq = "G|missense_variant|MODERATE|TP53|ENSG1|Transcript|ENST1|protein_coding|rs1|1|HGNC|11998|Ensembl";
        let ours_csq = "G|missense_variant|MODERATE|BRCA1|ENSG2|Transcript|ENST1|protein_coding|rs1|-1|HGNC|11998|RefSeq";
        let key = VariantKey {
            chrom: "22".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_alleles: "G".into(),
        };
        let golden = vec![VariantAnnotation {
            key: key.clone(),
            csq: Some(golden_csq.into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let ours = vec![VariantAnnotation {
            key,
            csq: Some(ours_csq.into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let report = compare_csq_fields(&golden, &ours);
        assert_eq!(report.total_entries_compared, 1);
        // Allele(0)=match, Consequence(1)=match, IMPACT(2)=match
        // SYMBOL(3)=mismatch, Gene(4)=mismatch
        // Feature_type(5)=match, Feature(6)=match, BIOTYPE(7)=match
        // Existing_variation(8)=match, STRAND(9)=mismatch
        // SYMBOL_SOURCE(10)=match, HGNC_ID(11)=match, SOURCE(12)=mismatch
        assert_eq!(report.field_match_counts[0], 1); // Allele
        assert_eq!(report.field_match_counts[1], 1); // Consequence
        assert_eq!(report.field_match_counts[2], 1); // IMPACT
        assert_eq!(report.field_match_counts[3], 0); // SYMBOL
        assert_eq!(report.field_match_counts[4], 0); // Gene
        assert_eq!(report.field_match_counts[5], 1); // Feature_type
        assert_eq!(report.field_match_counts[6], 1); // Feature
        assert_eq!(report.field_match_counts[7], 1); // BIOTYPE
        assert_eq!(report.field_match_counts[8], 1); // Existing_variation
        assert_eq!(report.field_match_counts[9], 0); // STRAND
        assert_eq!(report.field_match_counts[10], 1); // SYMBOL_SOURCE
        assert_eq!(report.field_match_counts[11], 1); // HGNC_ID
        assert_eq!(report.field_match_counts[12], 0); // SOURCE
    }

    #[test]
    fn compare_csq_fields_no_matching_variant_skipped() {
        let key_a = VariantKey {
            chrom: "1".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_alleles: "G".into(),
        };
        let key_b = VariantKey {
            chrom: "1".into(),
            pos: 200,
            ref_allele: "C".into(),
            alt_alleles: "T".into(),
        };
        let golden = vec![VariantAnnotation {
            key: key_a,
            csq: Some("G|missense_variant|MODERATE|S|G|Transcript|TX|pc|rs1|1|H|1|E".into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let ours = vec![VariantAnnotation {
            key: key_b,
            csq: Some("T|intron_variant|MODIFIER|S|G|Transcript|TX|pc|rs2|-1|H|1|E".into()),
            most_severe_consequence: Some("intron_variant".into()),
        }];
        let report = compare_csq_fields(&golden, &ours);
        assert_eq!(report.total_entries_compared, 0);
    }

    #[test]
    fn compare_csq_fields_multi_transcript_matching() {
        let csq_golden = "G|missense_variant|MODERATE|S|G1|Transcript|TX1|pc|rs1|1|H|1|E,G|intron_variant|MODIFIER|S|G1|Transcript|TX2|nc|rs1|-1|H|1|E";
        let csq_ours = "G|intron_variant|MODIFIER|S|G1|Transcript|TX2|nc|rs1|-1|H|1|E,G|missense_variant|MODERATE|S|G1|Transcript|TX1|pc|rs1|1|H|1|E";
        let key = VariantKey {
            chrom: "1".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_alleles: "G".into(),
        };
        let golden = vec![VariantAnnotation {
            key: key.clone(),
            csq: Some(csq_golden.into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let ours = vec![VariantAnnotation {
            key,
            csq: Some(csq_ours.into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let report = compare_csq_fields(&golden, &ours);
        // 2 golden entries matched by (Allele, Feature) key, all fields match
        assert_eq!(report.total_entries_compared, 2);
        for &c in &report.field_match_counts {
            assert_eq!(c, 2);
        }
    }

    #[test]
    fn diagnose_unmatched_reports_golden_only_and_ours_only() {
        let key = VariantKey {
            chrom: "22".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_alleles: "G".into(),
        };
        // Golden has TX1 and TX2; ours has TX1 and TX3.
        let golden = vec![VariantAnnotation {
            key: key.clone(),
            csq: Some(
                "G|miss|MOD|S|G1|Transcript|TX1|pc||||||||||||||||||||||||,\
                 G|intr|MOD|S|G1|Transcript|TX2|pc||||||||||||||||||||||||"
                    .into(),
            ),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let ours = vec![VariantAnnotation {
            key,
            csq: Some(
                "G|miss|MOD|S|G1|Transcript|TX1|pc||||||||||||||||||||||||,\
                 G|intr|MOD|S|G1|Transcript|TX3|pc||||||||||||||||||||||||"
                    .into(),
            ),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let report = diagnose_unmatched_csq(&golden, &ours);
        assert_eq!(report.golden_total, 2);
        assert_eq!(report.ours_total, 2);
        assert_eq!(report.matched, 1); // TX1 matched
        // TX2 is golden-only, TX3 is ours-only
        assert_eq!(report.golden_only_by_ft.get("Transcript"), Some(&1));
        assert_eq!(report.ours_only_by_ft.get("Transcript"), Some(&1));
        assert_eq!(report.golden_only_sample.len(), 1);
        assert_eq!(report.golden_only_sample[0].1, "TX2");
        assert_eq!(report.ours_only_sample.len(), 1);
        assert_eq!(report.ours_only_sample[0].1, "TX3");
    }

    #[test]
    fn diagnose_unmatched_all_matched() {
        let key = VariantKey {
            chrom: "22".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_alleles: "G".into(),
        };
        let csq = "G|miss|MOD|S|G1|Transcript|TX1|pc||||||||||||||||||||||||";
        let ann = vec![VariantAnnotation {
            key,
            csq: Some(csq.into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let report = diagnose_unmatched_csq(&ann, &ann);
        assert_eq!(report.matched, 1);
        assert!(report.golden_only_by_ft.is_empty());
        assert!(report.ours_only_by_ft.is_empty());
    }

    #[test]
    fn diagnose_multiplicity_reports_duplicate_ours_entries() {
        let key = VariantKey {
            chrom: "22".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_alleles: "G".into(),
        };
        let golden = vec![VariantAnnotation {
            key: key.clone(),
            csq: Some("G|miss|MOD|S|G1|Transcript|TX1|pc||||||||||||||||||||||||".into()),
            most_severe_consequence: Some("missense_variant".into()),
        }];
        let ours = vec![VariantAnnotation {
            key: key.clone(),
            csq: Some(
                "G|miss|MOD|S|G1|Transcript|TX1|pc||||||||||||||||||||||||,\
                 G|miss|MOD|S|G1|Transcript|TX1|pc||||||||||||||||||||||||"
                    .into(),
            ),
            most_severe_consequence: Some("missense_variant".into()),
        }];

        let report = diagnose_csq_multiplicity(&golden, &ours);
        assert_eq!(report.golden_extra_total, 0);
        assert_eq!(report.ours_extra_total, 1);
        assert_eq!(report.samples.len(), 1);
        assert_eq!(report.samples[0].variant_key, key);
        assert_eq!(report.samples[0].allele, "G");
        assert_eq!(report.samples[0].feature, "TX1");
        assert_eq!(report.samples[0].feature_type, "Transcript");
        assert_eq!(report.samples[0].golden_count, 1);
        assert_eq!(report.samples[0].ours_count, 2);
    }
}
