//! Allele matching function type for cache lookup.

/// Function signature for allele matching: `(vcf_ref, vcf_alt, allele_string) -> bool`.
pub type AlleleMatcher = fn(&str, &str, &str) -> bool;
