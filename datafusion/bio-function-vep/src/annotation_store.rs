//! Backend abstraction for VEP consequence annotation stores.
//!
//! Parquet is the only supported backend.

use datafusion::common::{DataFusionError, Result};

/// Supported annotation backend types. Parquet is the only backend; the variant
/// name is retained from the pre-migration era (a rename is tracked as part of
/// the lance-only dead-code removal).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationBackend {
    /// The Parquet point-lookup cache (the only backend).
    Parquet,
}

impl AnnotationBackend {
    /// Parse the backend from the UDTF argument. `"parquet"` is canonical;
    /// `"lance"` is accepted as a deprecated alias for the same Parquet backend.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "parquet" | "lance" => Ok(Self::Parquet),
            other => Err(DataFusionError::Plan(format!(
                "annotate_vep(): backend must be 'parquet' (deprecated alias 'lance'); got: {other}"
            ))),
        }
    }

    /// Stable display value for logs/debugging.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Parquet => "parquet",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AnnotationBackend;

    #[test]
    fn backend_parse_accepts_parquet_and_lance_alias() {
        assert_eq!(
            AnnotationBackend::parse("parquet").unwrap(),
            AnnotationBackend::Parquet
        );
        assert_eq!(
            AnnotationBackend::parse("lance").unwrap(),
            AnnotationBackend::Parquet
        );
    }

    #[test]
    fn backend_parse_rejects_unknown() {
        let err = AnnotationBackend::parse("local").unwrap_err().to_string();
        assert!(err.contains("backend must be 'parquet'"));
    }
}
