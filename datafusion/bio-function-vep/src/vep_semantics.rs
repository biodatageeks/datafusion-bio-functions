//! Compiled Ensembl VEP/cache compatibility contract.
//!
//! This is the single canonical support matrix used by the annotation engine
//! and exposed read-only to language bindings. Cache directory names are never
//! consulted when resolving semantics.

use datafusion::common::{DataFusionError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VepSemantics {
    V115,
    V116,
}

impl VepSemantics {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::V115 => "115",
            Self::V116 => "116",
        }
    }

    /// VEP 116 permits HGVS for variants that partially overlap a transcript.
    pub const fn partial_overlap_hgvs(self) -> bool {
        matches!(self, Self::V116)
    }

    /// VEP 116 changed the coordinated stop/inframe/frameshift predicates.
    pub const fn vep116_stop_predicates(self) -> bool {
        matches!(self, Self::V116)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SupportedVepTarget {
    pub cache_version: &'static str,
    pub vep_codebase_version: &'static str,
    pub api_version: &'static str,
    pub ensembl_core_revision: &'static str,
    pub ensembl_variation_revision: &'static str,
    pub semantics: VepSemantics,
}

pub const SUPPORTED_VEP_TARGETS: &[SupportedVepTarget] = &[
    SupportedVepTarget {
        cache_version: "115",
        vep_codebase_version: "115.2",
        api_version: "115",
        ensembl_core_revision: "266b84d",
        ensembl_variation_revision: "b7c2637",
        semantics: VepSemantics::V115,
    },
    SupportedVepTarget {
        cache_version: "116",
        vep_codebase_version: "116.0",
        api_version: "116",
        ensembl_core_revision: "c0cf13d",
        ensembl_variation_revision: "2fb834b",
        semantics: VepSemantics::V116,
    },
];

pub fn supported_vep_targets() -> &'static [SupportedVepTarget] {
    SUPPORTED_VEP_TARGETS
}

pub fn target_for_cache_version(cache_version: &str) -> Result<&'static SupportedVepTarget> {
    if cache_version.is_empty() || !cache_version.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(DataFusionError::Execution(format!(
            "invalid VEP cache version '{cache_version}': expected a decimal Ensembl release"
        )));
    }
    SUPPORTED_VEP_TARGETS
        .iter()
        .find(|target| target.cache_version == cache_version)
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "unsupported VEP cache version '{cache_version}'; supported cache versions: {}",
                SUPPORTED_VEP_TARGETS
                    .iter()
                    .map(|target| target.cache_version)
                    .collect::<Vec<_>>()
                    .join(", ")
            ))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn support_matrix_is_exact_and_unique() {
        assert_eq!(SUPPORTED_VEP_TARGETS.len(), 2);
        assert_eq!(
            SUPPORTED_VEP_TARGETS
                .iter()
                .map(|target| (
                    target.cache_version,
                    target.vep_codebase_version,
                    target.api_version,
                    target.ensembl_core_revision,
                    target.ensembl_variation_revision,
                    target.semantics.as_str(),
                ))
                .collect::<Vec<_>>(),
            vec![
                ("115", "115.2", "115", "266b84d", "b7c2637", "115"),
                ("116", "116.0", "116", "c0cf13d", "2fb834b", "116"),
            ]
        );
    }

    #[test]
    fn cache_version_resolution_rejects_malformed_and_unsupported_values() {
        assert_eq!(
            target_for_cache_version("115").unwrap().semantics,
            VepSemantics::V115
        );
        assert!(target_for_cache_version("115.2").is_err());
        assert!(target_for_cache_version("v116").is_err());
        assert!(target_for_cache_version("117").is_err());
    }

    #[test]
    fn only_two_semantic_policy_gates_exist() {
        assert!(!VepSemantics::V115.partial_overlap_hgvs());
        assert!(!VepSemantics::V115.vep116_stop_predicates());
        assert!(VepSemantics::V116.partial_overlap_hgvs());
        assert!(VepSemantics::V116.vep116_stop_predicates());
    }
}
