//! Coordinate system normalization between 0-based half-open and 1-based closed.

use datafusion::arrow::datatypes::SchemaRef;
use datafusion_bio_function_ranges::FilterOp;

/// Metadata key used to store coordinate system info in Arrow schema.
const COORDINATE_SYSTEM_METADATA_KEY: &str = "bio.coordinate_system_zero_based";

/// Normalizer for converting between coordinate systems.
///
/// Reads `bio.coordinate_system_zero_based` from both VCF and cache table schemas
/// and provides methods to convert coordinates for interval matching.
#[derive(Debug, Clone)]
pub struct CoordinateNormalizer {
    /// Whether the input (VCF) table uses 0-based coordinates.
    pub input_zero_based: bool,
    /// Whether the cache table uses 0-based coordinates.
    pub cache_zero_based: bool,
}

impl CoordinateNormalizer {
    /// Create a normalizer from two table schemas.
    pub fn from_schemas(input_schema: &SchemaRef, cache_schema: &SchemaRef) -> Self {
        Self {
            input_zero_based: is_zero_based(input_schema),
            cache_zero_based: is_zero_based(cache_schema),
        }
    }

    /// Create a normalizer with explicit coordinate system flags.
    pub fn new(input_zero_based: bool, cache_zero_based: bool) -> Self {
        Self {
            input_zero_based,
            cache_zero_based,
        }
    }

    /// Determine the `FilterOp` for interval joins based on coordinate systems.
    ///
    /// When both tables use the same coordinate system, use `Weak` (standard overlap).
    /// When they differ, use `Strict` to account for boundary semantics.
    pub fn filter_op(&self) -> FilterOp {
        if self.input_zero_based == self.cache_zero_based {
            FilterOp::Weak
        } else {
            FilterOp::Strict
        }
    }

    /// Returns true if both tables use the same coordinate system (no conversion needed).
    pub fn same_system(&self) -> bool {
        self.input_zero_based == self.cache_zero_based
    }
}

/// Read the coordinate system from schema metadata.
fn is_zero_based(schema: &SchemaRef) -> bool {
    schema
        .metadata()
        .get(COORDINATE_SYSTEM_METADATA_KEY)
        .is_some_and(|v| v == "true")
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::collections::HashMap;
    use std::sync::Arc;

    fn schema_with_coord(zero_based: bool) -> SchemaRef {
        let mut metadata = HashMap::new();
        metadata.insert(
            COORDINATE_SYSTEM_METADATA_KEY.to_string(),
            zero_based.to_string(),
        );
        Arc::new(Schema::new_with_metadata(
            vec![Field::new("chrom", DataType::Utf8, false)],
            metadata,
        ))
    }

    fn schema_without_coord() -> SchemaRef {
        Arc::new(Schema::new(vec![Field::new(
            "chrom",
            DataType::Utf8,
            false,
        )]))
    }

    #[test]
    fn same_system_both_zero_based() {
        let norm = CoordinateNormalizer::new(true, true);
        assert!(norm.same_system());
        assert_eq!(norm.filter_op(), FilterOp::Weak);
    }

    #[test]
    fn same_system_both_one_based() {
        let norm = CoordinateNormalizer::new(false, false);
        assert!(norm.same_system());
        assert_eq!(norm.filter_op(), FilterOp::Weak);
    }

    #[test]
    fn different_system_zero_to_one() {
        let norm = CoordinateNormalizer::new(true, false);
        assert!(!norm.same_system());
        assert_eq!(norm.filter_op(), FilterOp::Strict);
    }

    #[test]
    fn different_system_one_to_zero() {
        let norm = CoordinateNormalizer::new(false, true);
        assert!(!norm.same_system());
        assert_eq!(norm.filter_op(), FilterOp::Strict);
    }

    #[test]
    fn from_schemas_zero_based() {
        let input = schema_with_coord(true);
        let cache = schema_with_coord(false);
        let norm = CoordinateNormalizer::from_schemas(&input, &cache);
        assert!(norm.input_zero_based);
        assert!(!norm.cache_zero_based);
    }

    #[test]
    fn from_schemas_missing_metadata_defaults_to_one_based() {
        let input = schema_without_coord();
        let cache = schema_with_coord(false);
        let norm = CoordinateNormalizer::from_schemas(&input, &cache);
        assert!(!norm.input_zero_based); // defaults to false (1-based)
        assert!(!norm.cache_zero_based);
    }
}
