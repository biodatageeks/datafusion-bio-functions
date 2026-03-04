//! Utilities for golden benchmark comparison against Ensembl VEP output.
//!
//! This module is used by the benchmark example that:
//! - samples a fixed number of variants from a gzipped VCF,
//! - runs Ensembl VEP in Docker as golden standard,
//! - compares the golden output with `annotate_vep()` output.

use std::collections::HashMap;
use std::fs::{self, File};
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use datafusion::common::{DataFusionError, Result};
use flate2::read::MultiGzDecoder;

pub const DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ: &str =
    "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz";
pub const DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ_TBI: &str =
    "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz.tbi";
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
pub fn sample_gz_vcf_first_n(input_vcf_gz: &Path, output_vcf: &Path, limit: usize) -> Result<usize> {
    let file = File::open(input_vcf_gz).map_err(io_err)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let out_file = File::create(output_vcf).map_err(io_err)?;
    let mut writer = BufWriter::new(out_file);

    sample_vcf_lines(reader, &mut writer, limit)
}

fn sample_vcf_lines<R: BufRead, W: Write>(mut reader: R, writer: &mut W, limit: usize) -> Result<usize> {
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

fn parse_vep_vcf_lines<R: BufRead>(
    mut reader: R,
    out: &mut Vec<VariantAnnotation>,
) -> Result<()> {
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
    let mut best: Option<(u8, String)> = None;

    for ann in csq.split(',') {
        let mut parts = ann.split('|');
        let _allele = parts.next();
        let consequence_field = parts.next().unwrap_or("");
        for term in consequence_field.split('&') {
            let rank = so_rank(term);
            if rank == u8::MAX {
                continue;
            }
            match &best {
                None => best = Some((rank, term.to_string())),
                Some((best_rank, _)) if rank < *best_rank => best = Some((rank, term.to_string())),
                _ => {}
            }
        }
    }

    best.map(|(_, term)| term)
}

pub fn compare_annotations(golden: &[VariantAnnotation], ours: &[VariantAnnotation]) -> ComparisonReport {
    let golden_map: HashMap<VariantKey, &VariantAnnotation> = golden.iter().map(|r| (r.key.clone(), r)).collect();
    let ours_map: HashMap<VariantKey, &VariantAnnotation> = ours.iter().map(|r| (r.key.clone(), r)).collect();

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

    let extra_in_ours = ours_map.keys().filter(|k| !golden_map.contains_key(*k)).count();

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

fn so_rank(term: &str) -> u8 {
    match term {
        "transcript_ablation" => 1,
        "splice_acceptor_variant" => 2,
        "splice_donor_variant" => 3,
        "stop_gained" => 4,
        "frameshift_variant" => 5,
        "stop_lost" => 6,
        "start_lost" => 7,
        "transcript_amplification" => 8,
        "feature_elongation" => 9,
        "feature_truncation" => 10,
        "inframe_insertion" => 11,
        "inframe_deletion" => 12,
        "missense_variant" => 13,
        "protein_altering_variant" => 14,
        "splice_donor_5th_base_variant" => 15,
        "splice_region_variant" => 16,
        "splice_donor_region_variant" => 17,
        "splice_polypyrimidine_tract_variant" => 18,
        "incomplete_terminal_codon_variant" => 19,
        "start_retained_variant" => 20,
        "stop_retained_variant" => 21,
        "synonymous_variant" => 22,
        "coding_sequence_variant" => 23,
        "mature_miRNA_variant" => 24,
        "5_prime_UTR_variant" => 25,
        "3_prime_UTR_variant" => 26,
        "non_coding_transcript_exon_variant" => 27,
        "intron_variant" => 28,
        "NMD_transcript_variant" => 29,
        "non_coding_transcript_variant" => 30,
        "coding_transcript_variant" => 31,
        "upstream_gene_variant" => 32,
        "downstream_gene_variant" => 33,
        "TFBS_ablation" => 34,
        "TFBS_amplification" => 35,
        "TF_binding_site_variant" => 36,
        "regulatory_region_ablation" => 37,
        "regulatory_region_amplification" => 38,
        "regulatory_region_variant" => 39,
        "intergenic_variant" => 40,
        "sequence_variant" => 41,
        _ => u8::MAX,
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
}
