//! Backend abstraction for VEP consequence annotation stores.
//!
//! This module defines a backend-neutral contract for reading annotation
//! features from either Parquet datasets or Fjall-backed caches.

use datafusion::common::{DataFusionError, Result};

/// Supported annotation backend types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationBackend {
    Parquet,
    Fjall,
}

impl AnnotationBackend {
    /// Parse backend from UDTF argument.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "parquet" => Ok(Self::Parquet),
            "fjall" => Ok(Self::Fjall),
            other => Err(DataFusionError::Plan(format!(
                "annotate_vep() backend must be one of: parquet, fjall; got: {other}"
            ))),
        }
    }

    /// Stable display value for logs/debugging.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Parquet => "parquet",
            Self::Fjall => "fjall",
        }
    }
}

/// Shared contract for annotation stores.
///
/// This is intentionally minimal in phase 1. The concrete reader APIs
/// for transcript/exon/variation/regulatory/motif/miRNA features will be
/// introduced as consequence calculators are integrated.
pub trait AnnotationStore: std::fmt::Debug + Send + Sync {
    fn backend(&self) -> AnnotationBackend;
    fn source(&self) -> &str;
}

/// Parquet-backed annotation store descriptor.
#[derive(Debug, Clone)]
pub struct ParquetAnnotationStore {
    source: String,
}

impl ParquetAnnotationStore {
    pub fn new(source: String) -> Self {
        Self { source }
    }
}

impl AnnotationStore for ParquetAnnotationStore {
    fn backend(&self) -> AnnotationBackend {
        AnnotationBackend::Parquet
    }

    fn source(&self) -> &str {
        &self.source
    }
}

/// Fjall-backed annotation store descriptor.
#[derive(Debug, Clone)]
pub struct FjallAnnotationStore {
    source: String,
}

impl FjallAnnotationStore {
    pub fn new(source: String) -> Self {
        Self { source }
    }
}

impl AnnotationStore for FjallAnnotationStore {
    fn backend(&self) -> AnnotationBackend {
        AnnotationBackend::Fjall
    }

    fn source(&self) -> &str {
        &self.source
    }
}

/// Build a store descriptor for the selected backend.
pub fn build_store(backend: AnnotationBackend, source: String) -> Box<dyn AnnotationStore> {
    match backend {
        AnnotationBackend::Parquet => Box::new(ParquetAnnotationStore::new(source)),
        AnnotationBackend::Fjall => Box::new(FjallAnnotationStore::new(source)),
    }
}

#[cfg(test)]
mod tests {
    use super::AnnotationBackend;

    #[test]
    fn backend_parse_ok() {
        assert_eq!(
            AnnotationBackend::parse("parquet").unwrap(),
            AnnotationBackend::Parquet
        );
        assert_eq!(
            AnnotationBackend::parse("fjall").unwrap(),
            AnnotationBackend::Fjall
        );
    }

    #[test]
    fn backend_parse_rejects_unknown() {
        let err = AnnotationBackend::parse("unknown").unwrap_err().to_string();
        assert!(err.contains("annotate_vep() backend must be one of"));
    }
}
