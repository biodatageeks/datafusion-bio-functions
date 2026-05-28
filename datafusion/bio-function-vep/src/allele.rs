//! VCF ↔ VEP allele conversion and matching.

use datafusion::arrow::array::builder::StringBuilder;
use datafusion::arrow::array::{Array, ArrayRef, BooleanArray, Int64Array, StringArray};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::Result;
use datafusion::logical_expr::{ColumnarValue, ScalarUDF, Volatility, create_udf};
use std::collections::HashSet;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MatchedVariantAllele {
    pub a_allele: String,
    pub a_index: usize,
    pub b_allele: String,
    pub b_index: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct VariantAlleleInput<'a> {
    pub allele_string: &'a str,
    pub pos: i64,
    pub strand: i8,
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// This is a source-equivalent port of the release/115 allele minimization
/// used by Ensembl VEP's matched-alleles flow.
pub fn trim_sequences_ensembl_with_end(
    ref_allele: &str,
    alt_allele: &str,
    start: i64,
    end: i64,
    end_first: bool,
    strand: i8,
) -> (String, String, i64, i64, bool) {
    let mut ref_allele = ref_allele.as_bytes().to_vec();
    let mut alt_allele = alt_allele.as_bytes().to_vec();
    let mut start = start;
    let mut end = end;
    let mut changed = false;

    if end_first {
        while !ref_allele.is_empty()
            && !alt_allele.is_empty()
            && ref_allele[ref_allele.len() - 1] == alt_allele[alt_allele.len() - 1]
        {
            ref_allele.pop();
            alt_allele.pop();
            if strand == -1 {
                start += 1;
            } else {
                end -= 1;
            }
            changed = true;
        }

        while !ref_allele.is_empty() && !alt_allele.is_empty() && ref_allele[0] == alt_allele[0] {
            ref_allele.remove(0);
            alt_allele.remove(0);
            if strand == -1 {
                end -= 1;
            } else {
                start += 1;
            }
            changed = true;
        }
    } else {
        while !ref_allele.is_empty() && !alt_allele.is_empty() && ref_allele[0] == alt_allele[0] {
            ref_allele.remove(0);
            alt_allele.remove(0);
            if strand == -1 {
                end -= 1;
            } else {
                start += 1;
            }
            changed = true;
        }

        while !ref_allele.is_empty()
            && !alt_allele.is_empty()
            && ref_allele[ref_allele.len() - 1] == alt_allele[alt_allele.len() - 1]
        {
            ref_allele.pop();
            alt_allele.pop();
            if strand == -1 {
                start += 1;
            } else {
                end -= 1;
            }
            changed = true;
        }
    }

    let ref_allele = if ref_allele.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8(ref_allele).unwrap_or_default()
    };
    let alt_allele = if alt_allele.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8(alt_allele).unwrap_or_default()
    };

    (ref_allele, alt_allele, start, end, changed)
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Convenience wrapper for callers that start from a REF allele span where the
/// end coordinate can be derived directly from `start + len(ref) - 1`.
pub fn trim_sequences_ensembl(
    ref_allele: &str,
    alt_allele: &str,
    start: i64,
    end_first: bool,
    strand: i8,
) -> (String, String, i64, i64, bool) {
    let end = start + ref_allele.len() as i64 - 1;
    trim_sequences_ensembl_with_end(ref_allele, alt_allele, start, end, end_first, strand)
}

/// Traceability:
/// - Ensembl Variation `get_matched_variant_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L1098-L1258>
///
/// Upstream matched-allele comparison reverse-complements parser alleles when
/// the compared features are on opposite strands before sequence trimming.
fn reverse_complement_ascii(seq: &str) -> Option<String> {
    let mut out = String::with_capacity(seq.len());
    for b in seq.as_bytes().iter().rev() {
        out.push(match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            b'-' => '-',
            _ => return None,
        });
    }
    Some(out)
}

/// Traceability:
/// - Ensembl Variation `get_matched_variant_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L1098-L1258>
///
/// VEP matched-allele logic operates on `allele_string` in `REF/ALT...`
/// form and ignores malformed/non-variant strings before comparison.
fn parse_variant_allele_string(allele_string: &str) -> Option<(&str, Vec<&str>)> {
    if allele_string.starts_with('/') || !allele_string.contains('/') {
        return None;
    }
    let mut parts = allele_string.split('/');
    let ref_allele = parts.next()?;
    let alts: Vec<&str> = parts.collect();
    if alts.is_empty() {
        return None;
    }
    Some((ref_allele, alts))
}

/// Traceability:
/// - Ensembl Variation `get_matched_variant_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L1098-L1258>
///
/// Upstream allele minimization checks both left-first and right-first trim
/// order for non-SNV alleles and only a single pass for SNVs.
fn trim_directions(ref_allele: &str, alt_allele: &str) -> &'static [bool] {
    if ref_allele.len() > 1 || alt_allele.len() > 1 {
        &[false, true]
    } else {
        &[false]
    }
}

/// Traceability:
/// - Ensembl Variation `get_matched_variant_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L1098-L1258>
///
/// This is a source-equivalent port of the release/115 matched-alleles
/// algorithm used by Ensembl VEP `compare_existing()`.
pub fn get_matched_variant_alleles(
    a: VariantAlleleInput<'_>,
    b: VariantAlleleInput<'_>,
) -> Vec<MatchedVariantAllele> {
    let Some((a_ref_raw, a_alts_raw)) = parse_variant_allele_string(a.allele_string) else {
        return Vec::new();
    };
    let Some((b_ref_raw, b_alts_raw)) = parse_variant_allele_string(b.allele_string) else {
        return Vec::new();
    };
    if a.pos == 0 || b.pos == 0 {
        return Vec::new();
    }

    let mut a_ref = a_ref_raw.to_string();
    if a.strand != b.strand {
        let Some(reverse_ref) = reverse_complement_ascii(&a_ref) else {
            return Vec::new();
        };
        a_ref = reverse_ref;
    }

    let mut minimised_a_alleles: Vec<(String, String, String, usize)> = Vec::new();
    for (a_index, orig_a_alt) in a_alts_raw.iter().enumerate() {
        let mut a_alt = (*orig_a_alt).to_string();
        if a.strand != b.strand {
            let Some(reverse_alt) = reverse_complement_ascii(&a_alt) else {
                return Vec::new();
            };
            a_alt = reverse_alt;
        }

        for &end_first in trim_directions(&a_ref, orig_a_alt) {
            let (trimmed_ref, trimmed_alt, trimmed_pos, _, _) =
                trim_sequences_ensembl(&a_ref, &a_alt, a.pos, end_first, 1);
            minimised_a_alleles.push((
                format!("{trimmed_ref}_{trimmed_alt}_{trimmed_pos}"),
                (*orig_a_alt).to_string(),
                a_alt.clone(),
                a_index,
            ));
        }
    }

    let mut matches = Vec::new();
    let mut seen = HashSet::new();

    for (b_index, orig_b_alt) in b_alts_raw.iter().enumerate() {
        for &end_first in trim_directions(b_ref_raw, orig_b_alt) {
            let (trimmed_ref, trimmed_alt, trimmed_pos, _, _) =
                trim_sequences_ensembl(b_ref_raw, orig_b_alt, b.pos, end_first, 1);
            let key = format!("{trimmed_ref}_{trimmed_alt}_{trimmed_pos}");

            if let Some((_, orig_a_alt, _, a_index)) = minimised_a_alleles
                .iter()
                .find(|(candidate, ..)| candidate == &key)
            {
                let matched = MatchedVariantAllele {
                    a_allele: orig_a_alt.clone(),
                    a_index: *a_index,
                    b_allele: (*orig_b_alt).to_string(),
                    b_index,
                };
                if seen.insert(matched.clone()) {
                    matches.push(matched);
                }
                break;
            }
        }
    }

    matches
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Convert VCF REF/ALT pair to VEP allele format.
///
/// Strips shared prefix and suffix between REF and ALT to produce minimal
/// VEP-style allele representations.
///
/// VCF: REF="ACGT", ALT="A"    → VEP: "CGT/-"   (deletion)
/// VCF: REF="A", ALT="ACGT"    → VEP: "-/CGT"   (insertion)
/// VCF: REF="A", ALT="G"       → VEP: "A/G"     (SNV)
/// VCF: REF="AC", ALT="GT"     → VEP: "AC/GT"   (MNV)
/// VCF: REF="TCAC", ALT="T"    → VEP: "CAC/-"   (deletion, prefix+suffix)
/// VCF: REF="ATCG", ALT="AGCG" → VEP: "TCG/GCG"  (MNV: prefix-only trim, no suffix trim)
///
/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
pub fn vcf_to_vep_allele(ref_allele: &str, alt_allele: &str) -> (String, String) {
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        // SNV
        return (ref_allele.to_string(), alt_allele.to_string());
    }

    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    // Find common prefix
    let prefix_len = ref_bytes
        .iter()
        .zip(alt_bytes.iter())
        .take_while(|(a, b)| a == b)
        .count();

    // Find common suffix (not overlapping with prefix).
    // VEP only suffix-trims indels (different-length alleles), not MNVs (same-length substitutions).
    let mut suffix_len = 0;
    if ref_bytes.len() != alt_bytes.len() {
        let ref_remaining = ref_bytes.len() - prefix_len;
        let alt_remaining = alt_bytes.len() - prefix_len;
        while suffix_len < ref_remaining
            && suffix_len < alt_remaining
            && ref_bytes[ref_bytes.len() - 1 - suffix_len]
                == alt_bytes[alt_bytes.len() - 1 - suffix_len]
        {
            suffix_len += 1;
        }
    }

    let ref_trimmed = &ref_allele[prefix_len..ref_allele.len() - suffix_len];
    let alt_trimmed = &alt_allele[prefix_len..alt_allele.len() - suffix_len];

    let vep_ref = if ref_trimmed.is_empty() {
        "-".to_string()
    } else {
        ref_trimmed.to_string()
    };
    let vep_alt = if alt_trimmed.is_empty() {
        "-".to_string()
    } else {
        alt_trimmed.to_string()
    };

    (vep_ref, vep_alt)
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Multi-ALT analog of `vep_prefix_suffix_len`: returns the longest prefix
/// and suffix shared between REF and EVERY alt in `alts`. For single-ALT
/// inputs the result is identical to the biallelic version.
///
/// Suffix trimming only fires when at least one alt differs in length from
/// REF (`is_indel`), matching VEP's `trim_sequences()` MNV-preserving rule.
fn vep_prefix_suffix_len_multi(ref_allele: &str, alts: &[&str]) -> (usize, usize) {
    if alts.is_empty() {
        return (0, 0);
    }
    let ref_bytes = ref_allele.as_bytes();
    if ref_bytes.is_empty() {
        return (0, 0);
    }

    // Per-VEP: SNV-only multi-ALT (all single-base) gets no trimming.
    if ref_bytes.len() == 1 && alts.iter().all(|a| a.as_bytes().len() == 1) {
        return (0, 0);
    }

    // Longest prefix shared with ref AND every alt.
    let mut prefix_len = 0;
    let max_prefix = std::cmp::min(
        ref_bytes.len(),
        alts.iter().map(|a| a.as_bytes().len()).min().unwrap_or(0),
    );
    while prefix_len < max_prefix {
        let r = ref_bytes[prefix_len];
        if !alts
            .iter()
            .all(|a| a.as_bytes().get(prefix_len) == Some(&r))
        {
            break;
        }
        prefix_len += 1;
    }

    // Longest suffix shared with ref AND every alt, only when at least one
    // alt has a different length from REF (indel context).
    let any_indel = alts.iter().any(|a| a.as_bytes().len() != ref_bytes.len());
    let mut suffix_len = 0;
    if any_indel {
        let ref_remaining = ref_bytes.len() - prefix_len;
        let min_alt_remaining = alts
            .iter()
            .map(|a| a.as_bytes().len().saturating_sub(prefix_len))
            .min()
            .unwrap_or(0);
        let max_suffix = std::cmp::min(ref_remaining, min_alt_remaining);
        while suffix_len < max_suffix {
            let r = ref_bytes[ref_bytes.len() - 1 - suffix_len];
            if !alts.iter().all(|a| {
                let ab = a.as_bytes();
                ab.get(ab.len() - 1 - suffix_len) == Some(&r)
            }) {
                break;
            }
            suffix_len += 1;
        }
    }

    (prefix_len, suffix_len)
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Multi-ALT VEP allele normalization: trims the shared prefix/suffix across
/// REF and ALL alts jointly. This matches VEP's behavior for multi-allelic
/// VCF rows where the joint normalization may leave ALT bases that would
/// have been stripped under per-ALT trimming.
///
/// Example: `C → T,CCGC` has no shared prefix (T differs), so `Allele` for
/// the `CCGC` ALT stays `CCGC` (not the per-ALT-trimmed `CGC`). This is
/// load-bearing for `port_cache_variation` cv_02 and cv_04.
pub fn vcf_to_vep_allele_multi(ref_allele: &str, alts: &[&str]) -> (String, Vec<String>) {
    if alts.is_empty() {
        return (ref_allele.to_string(), Vec::new());
    }
    let (prefix_len, suffix_len) = vep_prefix_suffix_len_multi(ref_allele, alts);

    let trim_dash = |s: &str| -> String {
        if s.is_empty() {
            "-".to_string()
        } else {
            s.to_string()
        }
    };

    let vep_ref = trim_dash(&ref_allele[prefix_len..ref_allele.len() - suffix_len]);
    let vep_alts: Vec<String> = alts
        .iter()
        .map(|a| trim_dash(&a[prefix_len..a.len() - suffix_len]))
        .collect();

    (vep_ref, vep_alts)
}

/// Multi-ALT analog of `vep_norm_start`: shift by the joint common prefix.
pub fn vep_norm_start_multi(vcf_pos: i64, ref_allele: &str, alts: &[&str]) -> i64 {
    let (prefix_len, _) = vep_prefix_suffix_len_multi(ref_allele, alts);
    vcf_pos + prefix_len as i64
}

/// Multi-ALT analog of `vep_norm_end`: shrink by the joint common suffix.
pub fn vep_norm_end_multi(vcf_pos: i64, ref_allele: &str, alts: &[&str]) -> i64 {
    let (_, suffix_len) = vep_prefix_suffix_len_multi(ref_allele, alts);
    vcf_pos + ref_allele.len() as i64 - 1 - suffix_len as i64
}

/// Traceability:
/// - Ensembl VEP `Parser::VCF::create_VariationFeatures()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/Parser/VCF.pm#L321-L345>
///
/// Multi-ALT analog of `vcf_to_vep_input_allele`: strips a shared leading
/// anchor base across REF and ALL alts (parser-level only, no suffix trim).
/// Used to derive the sink-identity allele string for multi-allelic VCF rows
/// so the lookup key matches between the build side (variant_lookup_exec)
/// and the annotation side (annotate_provider).
pub fn vcf_to_vep_input_allele_multi(
    pos: i64,
    ref_allele: &str,
    alts: &[&str],
) -> (String, Vec<String>, i64) {
    if alts.is_empty() {
        return (ref_allele.to_string(), Vec::new(), pos);
    }

    let any_indel = alts.iter().any(|a| a.len() != ref_allele.len());
    let all_anchored = !ref_allele.is_empty()
        && alts
            .iter()
            .all(|a| !a.is_empty() && a.as_bytes()[0] == ref_allele.as_bytes()[0]);

    if any_indel && all_anchored {
        let trim_dash = |s: &str| -> String {
            if s.is_empty() {
                "-".to_string()
            } else {
                s.to_string()
            }
        };
        return (
            trim_dash(&ref_allele[1..]),
            alts.iter().map(|a| trim_dash(&a[1..])).collect(),
            pos + 1,
        );
    }

    (
        ref_allele.to_string(),
        alts.iter().map(|a| a.to_string()).collect(),
        pos,
    )
}

/// Traceability:
/// - Ensembl VEP `Parser::VCF::create_VariationFeatures()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/Parser/VCF.pm#L321-L345>
///
/// Converts a biallelic VCF REF/ALT pair into the parser-level
/// `VariationFeature->allele_string` representation used by upstream
/// `compare_existing()`: indels lose only a shared leading anchor base, and
/// the reported start coordinate is incremented when that happens. Unlike
/// `vcf_to_vep_allele()`, this does not perform suffix trimming.
pub fn vcf_to_vep_input_allele(
    pos: i64,
    ref_allele: &str,
    alt_allele: &str,
) -> (String, String, i64) {
    let is_indel = ref_allele.len() != 1 || alt_allele.len() != 1;
    if is_indel
        && !ref_allele.is_empty()
        && !alt_allele.is_empty()
        && ref_allele.as_bytes()[0] == alt_allele.as_bytes()[0]
    {
        let ref_trimmed = &ref_allele[1..];
        let alt_trimmed = &alt_allele[1..];
        return (
            if ref_trimmed.is_empty() {
                "-".to_string()
            } else {
                ref_trimmed.to_string()
            },
            if alt_trimmed.is_empty() {
                "-".to_string()
            } else {
                alt_trimmed.to_string()
            },
            pos + 1,
        );
    }

    (ref_allele.to_string(), alt_allele.to_string(), pos)
}

/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206>
/// - Ensembl Variation `get_matched_variant_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L1098-L1258>
///
/// Check if a VCF variant matches a VEP-format allele_string.
///
/// The cache `allele_string` is formatted as "REF/ALT1/ALT2" (e.g. "A/G", "A/G/T", "CGT/-").
/// This function converts the VCF REF/ALT to VEP format and checks that:
/// 1. The reference allele in `allele_string` matches the VEP-converted REF
///    (or the original VCF REF, for caches using VCF-format allele strings).
/// 2. At least one ALT allele in `allele_string` matches one of the input ALT alleles.
///
/// Some VCF readers expose multi-allelic rows as a single pipe-joined string
/// (`ALT = "A|T|..."`). We treat `|` as an ALT separator and match any ALT token.
///
/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206>
pub fn allele_matches(vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
    let mut parts = allele_string.split('/');

    let Some(cache_ref) = parts.next() else {
        return false;
    };
    let cache_alts: Vec<&str> = parts.collect();

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        // Left-first trimming (standard VEP normalization).
        let (vep_ref, vep_alt) = vcf_to_vep_allele(vcf_ref, alt);

        if (cache_ref == vep_ref || cache_ref == vcf_ref)
            && cache_alts.iter().any(|&a| a == vep_alt)
        {
            return true;
        }

        // Bidirectional matching: also trim each cache allele right-first
        // to handle cases where the cache stores a different representation
        // of the same variant (e.g., GCC/GCCCAGCC vs -/GCCCA).
        // Traceability: VEP get_matched_variant_alleles uses _get_trim_directions
        // which tries both direction=0 (left-first) and direction=1 (right-first).
        for &cache_alt in &cache_alts {
            let (trimmed_ref, trimmed_alt) = trim_right_first(cache_ref, cache_alt);
            if trimmed_ref == vep_ref && trimmed_alt == vep_alt {
                return true;
            }
        }
    }

    false
}

/// Trim shared suffix first, then shared prefix (right-first trimming).
/// Returns (ref, alt) after trimming, using "-" for empty strings.
fn trim_right_first(ref_allele: &str, alt_allele: &str) -> (String, String) {
    let mut r = ref_allele.as_bytes().to_vec();
    let mut a = alt_allele.as_bytes().to_vec();

    // Right-trim: strip matching chars from the end.
    while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
        r.pop();
        a.pop();
    }

    // Left-trim: strip matching chars from the start.
    while !r.is_empty() && !a.is_empty() && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
    }

    let ref_str = if r.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8(r).unwrap_or_default()
    };
    let alt_str = if a.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8(a).unwrap_or_default()
    };
    (ref_str, alt_str)
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Relaxed allele matching for Existing_variation-style fallback.
///
/// This matcher first applies strict `allele_matches()`. If strict matching fails,
/// it allows indel-only compatibility by event type + event length after
/// canonical trimming (common prefix/suffix), which helps bridge representation
/// differences for co-located indel IDs without consequence computation.
pub fn allele_matches_relaxed(vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
    if allele_matches(vcf_ref, vcf_alt, allele_string) {
        return true;
    }

    let mut parts = allele_string.split('/');
    let Some(cache_ref_allele) = parts.next() else {
        return false;
    };
    let cache_alt_alleles: Vec<&str> = parts.filter(|a| !a.is_empty()).collect();
    if cache_alt_alleles.is_empty() {
        return false;
    }

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (vcf_ref_len, vcf_alt_len) = canonical_event_lengths(vcf_ref, alt);
        let vcf_is_insertion = vcf_ref_len == 0 && vcf_alt_len > 0;
        let vcf_is_deletion = vcf_ref_len > 0 && vcf_alt_len == 0;
        if !vcf_is_insertion && !vcf_is_deletion {
            continue;
        }

        for cache_alt in &cache_alt_alleles {
            let (cache_ref_len, cache_alt_len) =
                canonical_event_lengths(cache_ref_allele, cache_alt);
            let cache_is_insertion = cache_ref_len == 0 && cache_alt_len > 0;
            let cache_is_deletion = cache_ref_len > 0 && cache_alt_len == 0;

            if vcf_is_insertion && cache_is_insertion && vcf_alt_len == cache_alt_len {
                return true;
            }
            if vcf_is_deletion && cache_is_deletion && vcf_ref_len == cache_ref_len {
                return true;
            }
        }
    }

    false
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Canonical event lengths are measured after the same common prefix/suffix
/// trimming VEP uses before comparing indel allele classes.
fn canonical_event_lengths(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_allele = if ref_allele == "-" { "" } else { ref_allele };
    let alt_allele = if alt_allele == "-" { "" } else { alt_allele };

    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    let mut ref_start = 0usize;
    let mut alt_start = 0usize;
    while ref_start < ref_bytes.len()
        && alt_start < alt_bytes.len()
        && ref_bytes[ref_start] == alt_bytes[alt_start]
    {
        ref_start += 1;
        alt_start += 1;
    }

    let mut ref_end = ref_bytes.len();
    let mut alt_end = alt_bytes.len();
    while ref_end > ref_start
        && alt_end > alt_start
        && ref_bytes[ref_end - 1] == alt_bytes[alt_end - 1]
    {
        ref_end -= 1;
        alt_end -= 1;
    }

    (ref_end - ref_start, alt_end - alt_start)
}

/// Create the `match_allele(vcf_ref, vcf_alt, allele_string)` scalar UDF.
///
/// Returns true if the VCF ALT allele matches any allele in the cache allele_string.
///
/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206>
pub fn match_allele_udf() -> ScalarUDF {
    create_udf(
        "match_allele",
        vec![DataType::Utf8, DataType::Utf8, DataType::Utf8],
        DataType::Boolean,
        Volatility::Immutable,
        Arc::new(match_allele_impl),
    )
}

/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206>
fn match_allele_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    match_allele_impl_with(args, allele_matches, "match_allele")
}

/// Create the `match_allele_relaxed(vcf_ref, vcf_alt, allele_string)` scalar UDF.
///
/// Returns true if strict allele matching succeeds, or if indel event type/length
/// are compatible under relaxed canonicalization.
///
/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
pub fn match_allele_relaxed_udf() -> ScalarUDF {
    create_udf(
        "match_allele_relaxed",
        vec![DataType::Utf8, DataType::Utf8, DataType::Utf8],
        DataType::Boolean,
        Volatility::Immutable,
        Arc::new(match_allele_relaxed_impl),
    )
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
fn match_allele_relaxed_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    match_allele_impl_with(args, allele_matches_relaxed, "match_allele_relaxed")
}

/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206>
fn match_allele_impl_with(
    args: &[ColumnarValue],
    matcher: fn(&str, &str, &str) -> bool,
    fn_name: &str,
) -> Result<ColumnarValue> {
    let refs = args[0].to_owned().into_array(1)?;
    let alts = args[1].to_owned().into_array(1)?;
    let allele_strs = args[2].to_owned().into_array(1)?;

    let refs = refs.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(format!("{fn_name}: first arg must be Utf8"))
    })?;
    let alts = alts.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(format!("{fn_name}: second arg must be Utf8"))
    })?;
    let allele_strings = allele_strs
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            datafusion::common::DataFusionError::Internal(format!(
                "{fn_name}: third arg must be Utf8"
            ))
        })?;

    let len = refs.len().max(alts.len()).max(allele_strings.len());
    let mut builder = BooleanArray::builder(len);

    for i in 0..len {
        let ref_idx = if refs.len() == 1 { 0 } else { i };
        let alt_idx = if alts.len() == 1 { 0 } else { i };
        let as_idx = if allele_strings.len() == 1 { 0 } else { i };

        if refs.is_null(ref_idx) || alts.is_null(alt_idx) || allele_strings.is_null(as_idx) {
            builder.append_null();
        } else {
            builder.append_value(matcher(
                refs.value(ref_idx),
                alts.value(alt_idx),
                allele_strings.value(as_idx),
            ));
        }
    }

    Ok(ColumnarValue::Array(Arc::new(builder.finish()) as ArrayRef))
}

/// Create the `vep_allele(vcf_ref, vcf_alt)` scalar UDF.
///
/// Returns the VEP-format allele string "ref/alt" for a VCF REF/ALT pair.
///
/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
pub fn vep_allele_udf() -> ScalarUDF {
    create_udf(
        "vep_allele",
        vec![DataType::Utf8, DataType::Utf8],
        DataType::Utf8,
        Volatility::Immutable,
        Arc::new(vep_allele_impl),
    )
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
fn vep_allele_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    let refs = args[0].to_owned().into_array(1)?;
    let alts = args[1].to_owned().into_array(1)?;

    let refs = refs.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(
            "vep_allele: first arg must be Utf8".to_string(),
        )
    })?;
    let alts = alts.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(
            "vep_allele: second arg must be Utf8".to_string(),
        )
    })?;

    let len = refs.len().max(alts.len());
    let mut builder = StringBuilder::new();

    for i in 0..len {
        let ref_idx = if refs.len() == 1 { 0 } else { i };
        let alt_idx = if alts.len() == 1 { 0 } else { i };

        if refs.is_null(ref_idx) || alts.is_null(alt_idx) {
            builder.append_null();
        } else {
            let (vep_ref, vep_alt) = vcf_to_vep_allele(refs.value(ref_idx), alts.value(alt_idx));
            builder.append_value(format!("{vep_ref}/{vep_alt}"));
        }
    }

    Ok(ColumnarValue::Array(Arc::new(builder.finish()) as ArrayRef))
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Compute the common prefix length between VCF REF and ALT alleles.
///
/// This is the number of leading bases that are identical in both alleles.
/// For SNVs (single-base REF and ALT), the prefix length is 0.
/// For indels, this equals the "anchor base" count that VEP strips when
/// normalizing coordinates.
fn vep_prefix_suffix_len(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    if ref_bytes.len() == 1 && alt_bytes.len() == 1 {
        return (0, 0);
    }

    let prefix_len = ref_bytes
        .iter()
        .zip(alt_bytes.iter())
        .take_while(|(a, b)| a == b)
        .count();

    // VEP only suffix-trims indels (different-length alleles), not MNVs.
    let mut suffix_len = 0;
    if ref_bytes.len() != alt_bytes.len() {
        let ref_remaining = ref_bytes.len() - prefix_len;
        let alt_remaining = alt_bytes.len() - prefix_len;
        while suffix_len < ref_remaining
            && suffix_len < alt_remaining
            && ref_bytes[ref_bytes.len() - 1 - suffix_len]
                == alt_bytes[alt_bytes.len() - 1 - suffix_len]
        {
            suffix_len += 1;
        }
    }

    (prefix_len, suffix_len)
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Compute VEP-normalized start coordinate from VCF POS, REF, ALT.
///
/// VEP strips the common prefix between REF and ALT and shifts the start
/// coordinate accordingly: `vep_start = vcf_pos + prefix_len`.
///
/// Examples:
/// - SNV `A>G` at POS=100: vep_start = 100 (no prefix)
/// - Deletion `CT>C` at POS=100: vep_start = 101 (prefix "C", len=1)
/// - Insertion `C>CT` at POS=100: vep_start = 101 (prefix "C", len=1)
///
/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
pub fn vep_norm_start(vcf_pos: i64, ref_allele: &str, alt_allele: &str) -> i64 {
    let (prefix_len, _) = vep_prefix_suffix_len(ref_allele, alt_allele);
    vcf_pos + prefix_len as i64
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
///
/// Compute VEP-normalized end coordinate from VCF POS, REF, ALT.
///
/// `vep_end = vcf_pos + len(REF) - 1 - suffix_len`
///
/// For insertions this produces `start > end` (VEP convention).
pub fn vep_norm_end(vcf_pos: i64, ref_allele: &str, alt_allele: &str) -> i64 {
    let (_, suffix_len) = vep_prefix_suffix_len(ref_allele, alt_allele);
    vcf_pos + ref_allele.len() as i64 - 1 - suffix_len as i64
}

/// Create the `vep_norm_start(pos, ref, alt)` scalar UDF.
///
/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
pub fn vep_norm_start_udf() -> ScalarUDF {
    create_udf(
        "vep_norm_start",
        vec![DataType::Int64, DataType::Utf8, DataType::Utf8],
        DataType::Int64,
        Volatility::Immutable,
        Arc::new(vep_norm_start_impl),
    )
}

/// Create the `vep_norm_end(pos, ref, alt)` scalar UDF.
///
/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
pub fn vep_norm_end_udf() -> ScalarUDF {
    create_udf(
        "vep_norm_end",
        vec![DataType::Int64, DataType::Utf8, DataType::Utf8],
        DataType::Int64,
        Volatility::Immutable,
        Arc::new(vep_norm_end_impl),
    )
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
fn vep_norm_start_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    vep_norm_coord_impl(args, true, "vep_norm_start")
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
fn vep_norm_end_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    vep_norm_coord_impl(args, false, "vep_norm_end")
}

/// Traceability:
/// - Ensembl Variation `trim_sequences()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L965-L1038>
fn vep_norm_coord_impl(
    args: &[ColumnarValue],
    is_start: bool,
    fn_name: &str,
) -> Result<ColumnarValue> {
    let positions = args[0].to_owned().into_array(1)?;
    let refs = args[1].to_owned().into_array(1)?;
    let alts = args[2].to_owned().into_array(1)?;

    let positions = positions
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| {
            datafusion::common::DataFusionError::Internal(format!(
                "{fn_name}: first arg must be Int64"
            ))
        })?;
    let refs = refs.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(format!("{fn_name}: second arg must be Utf8"))
    })?;
    let alts = alts.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(format!("{fn_name}: third arg must be Utf8"))
    })?;

    let len = positions.len().max(refs.len()).max(alts.len());
    let mut builder = Int64Array::builder(len);

    for i in 0..len {
        let pos_idx = if positions.len() == 1 { 0 } else { i };
        let ref_idx = if refs.len() == 1 { 0 } else { i };
        let alt_idx = if alts.len() == 1 { 0 } else { i };

        if positions.is_null(pos_idx) || refs.is_null(ref_idx) || alts.is_null(alt_idx) {
            builder.append_null();
        } else {
            let pos = positions.value(pos_idx);
            let r = refs.value(ref_idx);
            let a = alts.value(alt_idx);
            if is_start {
                builder.append_value(vep_norm_start(pos, r, a));
            } else {
                builder.append_value(vep_norm_end(pos, r, a));
            }
        }
    }

    Ok(ColumnarValue::Array(Arc::new(builder.finish()) as ArrayRef))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snv_conversion() {
        let (r, a) = vcf_to_vep_allele("A", "G");
        assert_eq!(r, "A");
        assert_eq!(a, "G");
    }

    #[test]
    fn test_trim_sequences_ensembl_left_first_deletion() {
        let (r, a, start, end, changed) = trim_sequences_ensembl("ACGT", "A", 100, false, 1);
        assert_eq!(
            (r.as_str(), a.as_str(), start, end, changed),
            ("CGT", "-", 101, 103, true)
        );
    }

    #[test]
    fn test_trim_sequences_ensembl_right_first_homopolymer() {
        let (r, a, start, end, changed) = trim_sequences_ensembl("AAAA", "AAA", 100, true, 1);
        assert_eq!(
            (r.as_str(), a.as_str(), start, end, changed),
            ("A", "-", 100, 100, true)
        );
    }

    #[test]
    fn test_get_matched_variant_alleles_repeat_shifted_deletion() {
        let matches = get_matched_variant_alleles(
            VariantAlleleInput {
                allele_string: "AAA/A",
                pos: 100,
                strand: 1,
            },
            VariantAlleleInput {
                allele_string: "AA/-",
                pos: 101,
                strand: 1,
            },
        );
        assert_eq!(
            matches,
            vec![MatchedVariantAllele {
                a_allele: "A".to_string(),
                a_index: 0,
                b_allele: "-".to_string(),
                b_index: 0,
            }]
        );
    }

    #[test]
    fn test_get_matched_variant_alleles_multiallelic() {
        let matches = get_matched_variant_alleles(
            VariantAlleleInput {
                allele_string: "A/G/T",
                pos: 100,
                strand: 1,
            },
            VariantAlleleInput {
                allele_string: "A/C/T",
                pos: 100,
                strand: 1,
            },
        );
        assert_eq!(
            matches,
            vec![MatchedVariantAllele {
                a_allele: "T".to_string(),
                a_index: 1,
                b_allele: "T".to_string(),
                b_index: 1,
            }]
        );
    }

    #[test]
    fn test_insertion_conversion() {
        // VCF: REF=A, ALT=ACGT → insertion of CGT
        let (r, a) = vcf_to_vep_allele("A", "ACGT");
        assert_eq!(r, "-");
        assert_eq!(a, "CGT");
    }

    #[test]
    fn test_deletion_conversion() {
        // VCF: REF=ACGT, ALT=A → deletion of CGT
        let (r, a) = vcf_to_vep_allele("ACGT", "A");
        assert_eq!(r, "CGT");
        assert_eq!(a, "-");
    }

    #[test]
    fn test_mnv_conversion() {
        let (r, a) = vcf_to_vep_allele("AC", "GT");
        assert_eq!(r, "AC");
        assert_eq!(a, "GT");
    }

    #[test]
    fn test_complex_indel_common_prefix() {
        // VCF: REF=ATCG, ALT=ATTT → common prefix AT, then CG→TT
        let (r, a) = vcf_to_vep_allele("ATCG", "ATTT");
        assert_eq!(r, "CG");
        assert_eq!(a, "TT");
    }

    #[test]
    fn test_allele_matches_snv() {
        assert!(allele_matches("A", "G", "A/G"));
        assert!(!allele_matches("A", "G", "A/T"));
    }

    #[test]
    fn test_allele_matches_multi_allelic() {
        assert!(allele_matches("A", "G", "A/G/T"));
        assert!(allele_matches("A", "T", "A/G/T"));
        assert!(!allele_matches("A", "C", "A/G/T"));
    }

    #[test]
    fn test_allele_matches_deletion_vep_format() {
        // VCF: REF=ACGT, ALT=A → VEP: ref="CGT", alt="-"
        // VEP cache stores prefix-stripped alleles: "CGT/-"
        assert!(allele_matches("ACGT", "A", "CGT/-"));
    }

    #[test]
    fn test_allele_matches_deletion_vcf_format() {
        // Some caches store the full VCF REF in allele_string: "ACGT/-"
        assert!(allele_matches("ACGT", "A", "ACGT/-"));
    }

    #[test]
    fn test_allele_matches_insertion() {
        // VCF: REF=A, ALT=ACGT → VEP alt = "CGT"
        assert!(allele_matches("A", "ACGT", "-/CGT"));
    }

    #[test]
    fn test_allele_does_not_match_ref() {
        // The first allele in allele_string is the reference — should not match
        assert!(!allele_matches("A", "A", "A/G"));
    }

    #[test]
    fn test_allele_rejects_ref_mismatch() {
        // Same ALT allele but different REF — should not match.
        // This catches false positives at adjacent positions with the same
        // common allele (e.g. "A/G" at position 100 and "C/G" at position 101).
        assert!(!allele_matches("A", "G", "C/G"));
        assert!(!allele_matches("A", "G", "T/G"));
    }

    #[test]
    fn test_allele_matches_mnv() {
        // MNV: REF=AC, ALT=GT → VEP: ref="AC", alt="GT"
        assert!(allele_matches("AC", "GT", "AC/GT"));
        // Wrong ref
        assert!(!allele_matches("AC", "GT", "TC/GT"));
    }

    #[test]
    fn test_allele_matches_pipe_joined_multi_alt_vcf_input() {
        assert!(allele_matches("A", "G|T", "A/G"));
        assert!(allele_matches("A", "G|T", "A/T"));
        assert!(!allele_matches("A", "G|T", "A/C"));
    }

    #[test]
    fn test_allele_matches_comma_joined_multi_alt_vcf_input() {
        assert!(allele_matches("A", "G,T", "A/G"));
        assert!(allele_matches("A", "G,T", "A/T"));
        assert!(!allele_matches("A", "G,T", "A/C"));
    }

    #[test]
    fn test_relaxed_insertion_length_match() {
        assert!(!allele_matches("A", "AT", "-/G"));
        assert!(allele_matches_relaxed("A", "AT", "-/G"));
    }

    #[test]
    fn test_relaxed_deletion_length_match() {
        assert!(!allele_matches("AA", "A", "C/-"));
        assert!(allele_matches_relaxed("AA", "A", "C/-"));
    }

    #[test]
    fn test_relaxed_does_not_relax_snv() {
        assert!(!allele_matches_relaxed("A", "G", "C/T"));
    }

    #[test]
    fn test_vcf_to_vep_allele_suffix_trim_deletion() {
        // REF=TCAC ALT=T → after prefix "T" trimmed: ref="CAC" alt=""
        // No shared suffix needed — pure prefix trim yields "-"
        let (r, a) = vcf_to_vep_allele("TCAC", "T");
        assert_eq!(r, "CAC");
        assert_eq!(a, "-");
    }

    #[test]
    fn test_vcf_to_vep_input_allele_only_trims_anchor_base() {
        let (r, a, pos) =
            vcf_to_vep_input_allele(62689175, "CATACATATATATATATATATATATAT", "CATATATATATATAT");
        assert_eq!(r, "ATACATATATATATATATATATATAT");
        assert_eq!(a, "ATATATATATATAT");
        assert_eq!(pos, 62689176);
    }

    #[test]
    fn test_vcf_to_vep_input_allele_insertion_uses_dash_after_anchor_trim() {
        let (r, a, pos) = vcf_to_vep_input_allele(100, "A", "ATG");
        assert_eq!(r, "-");
        assert_eq!(a, "TG");
        assert_eq!(pos, 101);
    }

    #[test]
    fn test_vcf_to_vep_allele_mnv_no_suffix_trim() {
        // REF=ATCG ALT=AGCG → same length (MNV), only prefix-trim "A"
        // VEP does NOT suffix-trim MNVs
        let (r, a) = vcf_to_vep_allele("ATCG", "AGCG");
        assert_eq!(r, "TCG");
        assert_eq!(a, "GCG");
    }

    #[test]
    fn test_vcf_to_vep_allele_mnv_2bp() {
        // REF=GT ALT=TT → same length, prefix-trim nothing (G != T)
        let (r, a) = vcf_to_vep_allele("GT", "TT");
        assert_eq!(r, "GT");
        assert_eq!(a, "TT");
    }

    #[test]
    fn test_vcf_to_vep_allele_mnv_with_prefix() {
        // REF=CTA ALT=ATA → same length, no common prefix (C != A)
        let (r, a) = vcf_to_vep_allele("CTA", "ATA");
        assert_eq!(r, "CTA");
        assert_eq!(a, "ATA");
    }

    #[test]
    fn test_vcf_to_vep_allele_no_suffix_needed() {
        // SNV: no trimming needed
        let (r, a) = vcf_to_vep_allele("A", "T");
        assert_eq!(r, "A");
        assert_eq!(a, "T");
    }

    #[test]
    fn test_vcf_to_vep_allele_prefix_and_suffix_insertion() {
        // REF=AG ALT=ATCG → shared prefix "A", shared suffix "G"
        // Net: ref="-" alt="TC"
        let (r, a) = vcf_to_vep_allele("AG", "ATCG");
        assert_eq!(r, "-");
        assert_eq!(a, "TC");
    }

    #[test]
    fn test_vcf_to_vep_allele_identity() {
        // REF=A ALT=A → SNV path, no trimming
        let (r, a) = vcf_to_vep_allele("A", "A");
        assert_eq!(r, "A");
        assert_eq!(a, "A");
    }

    #[test]
    fn test_vep_norm_coords_snv() {
        // SNV: A>G at POS=100 → no shift
        assert_eq!(vep_norm_start(100, "A", "G"), 100);
        assert_eq!(vep_norm_end(100, "A", "G"), 100);
    }

    #[test]
    fn test_vep_norm_coords_deletion() {
        // Deletion: CT>C at POS=35295124 → shift start by 1
        assert_eq!(vep_norm_start(35295124, "CT", "C"), 35295125);
        assert_eq!(vep_norm_end(35295124, "CT", "C"), 35295125);
    }

    #[test]
    fn test_vep_norm_coords_insertion() {
        // Insertion: C>CT at POS=32519310 → start > end (VEP convention)
        assert_eq!(vep_norm_start(32519310, "C", "CT"), 32519311);
        assert_eq!(vep_norm_end(32519310, "C", "CT"), 32519310);
    }

    #[test]
    fn test_vep_norm_coords_multi_base_deletion() {
        // Deletion: CCTGGTAGCA>C at POS=17817013 → shift start by 1, end = start + 8
        assert_eq!(vep_norm_start(17817013, "CCTGGTAGCA", "C"), 17817014);
        assert_eq!(vep_norm_end(17817013, "CCTGGTAGCA", "C"), 17817022);
    }

    #[test]
    fn test_vep_norm_coords_deletion_with_suffix() {
        // Deletion with suffix: TCAC>T at POS=100 → prefix "T" (len=1), suffix "C" (len=1 for indel)
        // But wait: vcf_to_vep_allele("TCAC","T") → prefix=1 ("T"), remaining ref="CAC", alt=""
        // No suffix since alt_remaining=0. So: start=101, end=103
        assert_eq!(vep_norm_start(100, "TCAC", "T"), 101);
        assert_eq!(vep_norm_end(100, "TCAC", "T"), 103);
    }

    // ---- Multi-ALT helpers ---------------------------------------------------

    #[test]
    fn test_multi_alt_single_alt_matches_biallelic_snv() {
        // Single ALT path must equal biallelic helper exactly.
        let (r, a) = vcf_to_vep_allele_multi("A", &["G"]);
        assert_eq!(r, "A");
        assert_eq!(a, vec!["G".to_string()]);
    }

    #[test]
    fn test_multi_alt_single_alt_matches_biallelic_insertion() {
        // C → CCGC: shared prefix C → vep_ref="-", alts=["CGC"], start shift = 1.
        let (r, a) = vcf_to_vep_allele_multi("C", &["CCGC"]);
        assert_eq!(r, "-");
        assert_eq!(a, vec!["CGC".to_string()]);
        assert_eq!(vep_norm_start_multi(8987004, "C", &["CCGC"]), 8987005);
    }

    #[test]
    fn test_multi_alt_snv_plus_insertion_no_trim() {
        // cv_02: C → T,CCGC. No shared prefix (T != C in alt0). VEP keeps both
        // alts untrimmed in the multi-allelic Allele string.
        let (r, a) = vcf_to_vep_allele_multi("C", &["T", "CCGC"]);
        assert_eq!(r, "C");
        assert_eq!(a, vec!["T".to_string(), "CCGC".to_string()]);
        // No prefix to strip → start stays put.
        assert_eq!(vep_norm_start_multi(8987004, "C", &["T", "CCGC"]), 8987004);
    }

    #[test]
    fn test_multi_alt_two_insertions_shared_prefix_trim() {
        // cv_04: C → CCGC,CGGC. Shared prefix C (both alts and ref start
        // with C) → trim 1 → ("-", ["CGC", "GGC"]) and start shifts by 1.
        let (r, a) = vcf_to_vep_allele_multi("C", &["CCGC", "CGGC"]);
        assert_eq!(r, "-");
        assert_eq!(a, vec!["CGC".to_string(), "GGC".to_string()]);
        assert_eq!(
            vep_norm_start_multi(8987004, "C", &["CCGC", "CGGC"]),
            8987005
        );
    }

    #[test]
    fn test_multi_alt_insertion_with_shared_suffix() {
        // cv_03: CT → CCGCT. Shared prefix C, shared suffix T (indel).
        // → ("-", ["CGC"]).
        let (r, a) = vcf_to_vep_allele_multi("CT", &["CCGCT"]);
        assert_eq!(r, "-");
        assert_eq!(a, vec!["CGC".to_string()]);
    }

    #[test]
    fn test_multi_alt_input_allele_multi_alt_no_anchor() {
        // C → T,CCGC has no fully-shared anchor base (ALT[0]=T differs).
        // Parser-level input-allele string must NOT strip anything.
        let (r, a, p) = vcf_to_vep_input_allele_multi(8987004, "C", &["T", "CCGC"]);
        assert_eq!(r, "C");
        assert_eq!(a, vec!["T".to_string(), "CCGC".to_string()]);
        assert_eq!(p, 8987004);
    }

    #[test]
    fn test_multi_alt_input_allele_shared_anchor() {
        // C → CCGC,CGGC: both ALTs share anchor C with REF. Parser-level
        // strips one base and shifts position by 1.
        let (r, a, p) = vcf_to_vep_input_allele_multi(8987004, "C", &["CCGC", "CGGC"]);
        assert_eq!(r, "-");
        assert_eq!(a, vec!["CGC".to_string(), "GGC".to_string()]);
        assert_eq!(p, 8987005);
    }

    // ----------------------------------------------------------------
    // Axis B (vepyr-side) tests for `AnnotationSource_Cache_Variation.t`.
    // See `porting-tests/detailed_plans/AnnotationSource_Cache_Variation.md`
    // §Axis-B additions B2 + B3.
    //
    // Perl's `minimise_alleles` trims regardless of REF/ALT shape. Vepyr's
    // `vcf_to_vep_allele_multi` deliberately diverges in two cases:
    //   * SNV-only multi-ALT → NO trimming (preserves the per-ALT base).
    //   * MNV (equal-length REF/ALT) multi-ALT → prefix-only trim, never
    //     suffix-trim (matches the documented comment block at
    //     `allele.rs:28` "VEP keeps the multi-allelic Allele string aligned
    //     to the per-position REF base").
    // ----------------------------------------------------------------

    // Axis B B2: subtest (Variation Phase D Axis B B2) — SNV-only multi-ALT
    // does NO trimming. See `allele.rs:351` for the source comment.
    #[test]
    fn axis_b_b2_vcf_to_vep_allele_multi_snv_only_no_trim() {
        // A → G,T: pure SNVs, no shared prefix beyond ALT[0]'s single base.
        // Perl `minimise_alleles` would trim the shared base; vepyr keeps
        // the per-position REF + per-allele ALT untouched.
        let (r, a) = vcf_to_vep_allele_multi("A", &["G", "T"]);
        assert_eq!(r, "A");
        assert_eq!(a, vec!["G".to_string(), "T".to_string()]);
        assert_eq!(vep_norm_start_multi(100, "A", &["G", "T"]), 100);
    }

    // Axis B B3: subtest (Variation Phase D Axis B B3) — Multi-ALT MNV
    // (equal-length REF/ALT) prefix-only trim, never suffix-trim. See
    // `allele.rs:28` for the design comment.
    #[test]
    fn axis_b_b3_vcf_to_vep_allele_multi_mnv_prefix_only_trim() {
        // AAA → AAT,ACA: equal-length MNV with shared prefix `A` (the
        // first base is common to REF + all ALTs; ALT2 also extends
        // the shared prefix to `A`+anything, but ALT1 diverges at
        // position 2 — so the longest shared prefix is `A`).
        //
        // Suffix is NOT shared (REF ends `A`, ALT1 ends `T`, ALT2 ends
        // `A`) — and even if it were, the MNV rule states prefix-only.
        let (r, a) = vcf_to_vep_allele_multi("AAA", &["AAT", "ACA"]);
        // Expectation: shared prefix `A` trimmed → REF="AA", ALTs=["AT","CA"].
        // If the helper instead strips both prefix `A` AND a shared
        // suffix, the test fails — preserving the Axis B B3 invariant.
        assert_eq!(r, "AA");
        assert_eq!(a, vec!["AT".to_string(), "CA".to_string()]);
    }

    // ─────────────────────────────────────────────────────────────────────
    // Port of `ensembl-vep/t/Parser.t` allele-minimisation rows
    //
    // Detailed plan: `porting-tests/detailed_plans/Parser.md`.
    // v2 paradigm: sztywno 1:1, standalone (no docker, no golden.vcf).
    // Anchor: `~/.claude/skills/port-to-vepyr/references/v2-paradigm.md`.
    //
    // Perl's `minimise_alleles` trims a shared prefix AND a shared suffix
    // off REF/ALT regardless of indel-vs-MNV shape (Parser.pm
    // `minimise_alleles()` walks both ends unconditionally).
    //
    // Vepyr's `vcf_to_vep_allele` deliberately preserves MNV form: when
    // `ref.len() == alt.len()` only the common prefix is trimmed — see the
    // `is_indel`/length-mismatch gate at allele.rs:300-302 with the
    // traceability comment to Ensembl-Variation `trim_sequences()`. This
    // is the vepyr-correct behaviour; the Perl test's `'A/T'` expectation
    // for `'GAC/GTC'` and the long-trim case reflects Perl's broader
    // (and arguably wrong-by-Ensembl-Variation-canon) trimming, so the
    // Rust assertions document the divergence with the
    // prefix-only-trim result.
    // ─────────────────────────────────────────────────────────────────────

    /// Port of `Parser.t` subtest #17 (Perl l.76-92):
    ///   `minimise_alleles('GA/GT')` → ref='A' alt='T', start shifted +1.
    ///
    /// Vepyr handles the indel-style shape `GA → GT` via `vcf_to_vep_allele`
    /// with length-equal (2bp/2bp) → MNV rule → prefix-only-trim → `(A, T)`.
    /// `vep_norm_start(1, "GA", "GT")` shifts start by the trimmed prefix
    /// length (1) — matching Perl's `start: 2, end: 1` post-trim convention.
    #[test]
    fn port_parser_subtest_17_minimise_alleles_ga_gt() {
        // Perl l.76-92: 'GA/GT' minimises to 'A/T' with start+=1.
        let (vep_ref, vep_alt) = vcf_to_vep_allele("GA", "GT");
        assert_eq!(vep_ref, "A");
        assert_eq!(vep_alt, "T");
        // Perl's `start = 2` corresponds to vepyr's `vep_norm_start = 1 + 1`.
        assert_eq!(vep_norm_start(1, "GA", "GT"), 2);
    }

    /// Port of `Parser.t` subtest #18 (Perl l.94-98):
    ///   `minimise_alleles('GAC/GTC')->[0]->{allele_string}` → 'A/T'.
    ///
    /// **Divergence (documented)**: Perl trims both shared prefix (`G`)
    /// AND shared suffix (`C`) regardless of indel/MNV shape → `A/T`.
    /// Vepyr deliberately preserves MNV form (Ensembl-Variation
    /// `trim_sequences()` canon, see allele.rs:300-302) →
    /// prefix-only-trim → `AC/TC`. This is intentional and
    /// load-bearing for vepyr correctness, not a bug.
    #[test]
    fn port_parser_subtest_18_minimise_alleles_gac_gtc_mnv_preserves_suffix() {
        // Perl asserts 'A/T'; vepyr correctly asserts 'AC/TC' because
        // the MNV-rule at allele.rs:300-302 forbids suffix-trim when
        // REF and ALT have equal length.
        let (vep_ref, vep_alt) = vcf_to_vep_allele("GAC", "GTC");
        assert_eq!(vep_ref, "AC", "vepyr preserves MNV suffix (divergent from Perl)");
        assert_eq!(vep_alt, "TC", "vepyr preserves MNV suffix (divergent from Perl)");
    }

    /// Port of `Parser.t` subtest #19 (Perl l.100-104):
    ///   `minimise_alleles('TTCCTTCCGACGGTACACACACACA/TTCCTTCCGTCGGTACACACACACA')`
    ///   → 'A/T' (Perl long-trim).
    ///
    /// **Divergence (documented)**: same as #18 — equal-length REF/ALT (25bp
    /// each) → MNV rule → prefix-only-trim. Vepyr trims the 9-byte common
    /// prefix `TTCCTTCCG` and stops; the remaining 16-byte suffixes
    /// `ACGGTACACACACACA` / `TCGGTACACACACACA` are preserved per the MNV
    /// rule. Perl would aggressively trim the shared suffix too, ending up
    /// with the minimal `A/T` differing bases.
    #[test]
    fn port_parser_subtest_19_minimise_alleles_long_mnv_preserves_suffix() {
        let ref_allele = "TTCCTTCCGACGGTACACACACACA";
        let alt_allele = "TTCCTTCCGTCGGTACACACACACA";
        let (vep_ref, vep_alt) = vcf_to_vep_allele(ref_allele, alt_allele);
        // 9-byte prefix "TTCCTTCCG" trimmed; remaining 16 bytes preserved.
        assert_eq!(vep_ref, "ACGGTACACACACACA");
        assert_eq!(vep_alt, "TCGGTACACACACACA");
    }

    /// Axis-B row B1 (post-Phase-D 2026-05-27): pin vepyr's behaviour on
    /// degenerate `REF == ALT` input.
    ///
    /// Detailed plan: `porting-tests/detailed_plans/Parser.md` row B1.
    /// Perl rejects same-allele input upstream of `minimise_alleles`; vepyr
    /// passes it through `vcf_to_vep_allele`, which returns `(REF, ALT)`
    /// verbatim for 1bp inputs (SNV fast path at line 284-287) and for
    /// multi-byte inputs returns `("-", "-")` after full prefix-and-suffix
    /// trim collapses both sides to empty.
    #[test]
    #[allow(non_snake_case)]
    fn port_parser_axisB1_vcf_to_vep_allele_degenerate_same_allele() {
        // 1bp same-allele goes through the SNV fast path → identity.
        let (r1, a1) = vcf_to_vep_allele("A", "A");
        assert_eq!((r1.as_str(), a1.as_str()), ("A", "A"));

        // 3bp same-allele: equal length → MNV rule → prefix-only-trim, full
        // 3-byte prefix matches → both sides empty → "-" dash.
        let (r3, a3) = vcf_to_vep_allele("AGT", "AGT");
        assert_eq!(
            (r3.as_str(), a3.as_str()),
            ("-", "-"),
            "degenerate same-allele 3bp collapses both sides to '-' under MNV rule"
        );
    }
}
