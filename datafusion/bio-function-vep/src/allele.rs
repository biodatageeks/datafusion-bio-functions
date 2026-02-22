//! VCF ↔ VEP allele conversion and matching.

use datafusion::arrow::array::builder::StringBuilder;
use datafusion::arrow::array::{Array, ArrayRef, BooleanArray, StringArray};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::Result;
use datafusion::logical_expr::{ColumnarValue, ScalarUDF, Volatility, create_udf};
use std::sync::Arc;

/// Convert VCF REF/ALT pair to VEP allele format.
///
/// VCF: REF="ACGT", ALT="A"   → VEP: "CGT/-"   (deletion)
/// VCF: REF="A", ALT="ACGT"   → VEP: "-/CGT"   (insertion)
/// VCF: REF="A", ALT="G"      → VEP: "A/G"     (SNV)
/// VCF: REF="AC", ALT="GT"    → VEP: "AC/GT"   (MNV)
pub fn vcf_to_vep_allele(ref_allele: &str, alt_allele: &str) -> (String, String) {
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        // SNV
        return (ref_allele.to_string(), alt_allele.to_string());
    }

    // Find common prefix
    let prefix_len = ref_allele
        .bytes()
        .zip(alt_allele.bytes())
        .take_while(|(a, b)| a == b)
        .count();

    let ref_trimmed = &ref_allele[prefix_len..];
    let alt_trimmed = &alt_allele[prefix_len..];

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

/// Check if an ALT allele matches any allele in a VEP-format allele_string.
///
/// The cache `allele_string` is formatted as "REF/ALT1/ALT2" (e.g. "A/G", "A/G/T", "ACGT/-").
/// This function converts the VCF ALT to VEP format and checks for a match.
pub fn allele_matches(vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
    let (_, vep_alt) = vcf_to_vep_allele(vcf_ref, vcf_alt);
    allele_string
        .split('/')
        .skip(1) // skip the reference allele
        .any(|a| a == vep_alt)
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
    let refs = args[0].to_owned().into_array(1)?;
    let alts = args[1].to_owned().into_array(1)?;
    let allele_strs = args[2].to_owned().into_array(1)?;

    let refs = refs.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(
            "match_allele: first arg must be Utf8".to_string(),
        )
    })?;
    let alts = alts.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        datafusion::common::DataFusionError::Internal(
            "match_allele: second arg must be Utf8".to_string(),
        )
    })?;
    let allele_strings = allele_strs
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            datafusion::common::DataFusionError::Internal(
                "match_allele: third arg must be Utf8".to_string(),
            )
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
            builder.append_value(allele_matches(
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
    fn test_allele_matches_deletion() {
        // VCF: REF=ACGT, ALT=A → VEP alt = "-"
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
}
