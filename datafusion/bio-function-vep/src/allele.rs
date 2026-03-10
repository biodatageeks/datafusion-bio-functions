//! VCF ↔ VEP allele conversion and matching.

use datafusion::arrow::array::builder::StringBuilder;
use datafusion::arrow::array::{Array, ArrayRef, BooleanArray, StringArray};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::Result;
use datafusion::logical_expr::{ColumnarValue, ScalarUDF, Volatility, create_udf};
use std::sync::Arc;

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
pub fn allele_matches(vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
    let mut parts = allele_string.split('/');

    let Some(ref_allele) = parts.next() else {
        return false;
    };

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (vep_ref, vep_alt) = vcf_to_vep_allele(vcf_ref, alt);

        // Verify reference allele: accept either VEP-format (prefix-stripped)
        // or VCF-format (full REF) since caches may use either convention.
        if ref_allele != vep_ref && ref_allele != vcf_ref {
            continue;
        }

        // Check if any ALT allele matches this ALT token.
        if parts.clone().any(|a| a == vep_alt) {
            return true;
        }
    }

    false
}

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
pub fn match_allele_udf() -> ScalarUDF {
    create_udf(
        "match_allele",
        vec![DataType::Utf8, DataType::Utf8, DataType::Utf8],
        DataType::Boolean,
        Volatility::Immutable,
        Arc::new(match_allele_impl),
    )
}

fn match_allele_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    match_allele_impl_with(args, allele_matches, "match_allele")
}

/// Create the `match_allele_relaxed(vcf_ref, vcf_alt, allele_string)` scalar UDF.
///
/// Returns true if strict allele matching succeeds, or if indel event type/length
/// are compatible under relaxed canonicalization.
pub fn match_allele_relaxed_udf() -> ScalarUDF {
    create_udf(
        "match_allele_relaxed",
        vec![DataType::Utf8, DataType::Utf8, DataType::Utf8],
        DataType::Boolean,
        Volatility::Immutable,
        Arc::new(match_allele_relaxed_impl),
    )
}

fn match_allele_relaxed_impl(args: &[ColumnarValue]) -> Result<ColumnarValue> {
    match_allele_impl_with(args, allele_matches_relaxed, "match_allele_relaxed")
}

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
pub fn vep_allele_udf() -> ScalarUDF {
    create_udf(
        "vep_allele",
        vec![DataType::Utf8, DataType::Utf8],
        DataType::Utf8,
        Volatility::Immutable,
        Arc::new(vep_allele_impl),
    )
}

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
}
