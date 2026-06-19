//! Backend abstraction for VEP consequence annotation stores.
//!
//! Lance is the only supported backend.

use datafusion::common::{DataFusionError, Result};

/// Supported annotation backend types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationBackend {
    Lance,
}

impl AnnotationBackend {
    /// Parse backend from UDTF argument. Only `"lance"` is accepted.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "lance" => Ok(Self::Lance),
            other => Err(DataFusionError::Plan(format!(
                "annotate_vep(): backend must be 'lance'; got: {other}"
            ))),
        }
    }

    /// Stable display value for logs/debugging.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Lance => "lance",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AnnotationBackend;

    #[test]
    fn backend_parse_ok() {
        assert_eq!(
            AnnotationBackend::parse("lance").unwrap(),
            AnnotationBackend::Lance
        );
    }

    #[test]
    fn backend_parse_rejects_unknown() {
        let err = AnnotationBackend::parse("parquet").unwrap_err().to_string();
        assert!(err.contains("backend must be 'lance'"));
    }
}
