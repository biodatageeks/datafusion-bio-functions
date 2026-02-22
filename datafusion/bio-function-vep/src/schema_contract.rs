//! Schema contract validation for VEP cache tables.
//!
//! Defines expected column names and types for variation cache tables,
//! and validates that a `TableProvider` exposes the required columns.

use datafusion::arrow::datatypes::{DataType, SchemaRef};
use datafusion::common::{DataFusionError, Result};

/// Required columns that must be present in the variation cache table.
pub const REQUIRED_VARIATION_COLUMNS: &[(&str, DataType)] = &[
    ("chrom", DataType::Utf8),
    ("start", DataType::Int64),
    ("end", DataType::Int64),
    ("variation_name", DataType::Utf8),
    ("allele_string", DataType::Utf8),
];

/// Default columns returned by `lookup_variants()` when no explicit column list is provided.
pub const DEFAULT_LOOKUP_COLUMNS: &[&str] = &["variation_name", "allele_string", "clin_sig"];

/// Check if two data types are compatible (treating Utf8/Utf8View/LargeUtf8 as
/// interchangeable, since DataFusion 50+ reads parquet strings as Utf8View).
fn types_compatible(actual: &DataType, expected: &DataType) -> bool {
    if actual == expected {
        return true;
    }
    matches!(
        (actual, expected),
        (DataType::Utf8View | DataType::LargeUtf8, DataType::Utf8)
            | (DataType::Utf8, DataType::Utf8View | DataType::LargeUtf8)
    )
}

/// Validate that a variation cache table schema contains all required columns.
pub fn validate_variation_schema(schema: &SchemaRef) -> Result<()> {
    for (name, expected_type) in REQUIRED_VARIATION_COLUMNS {
        match schema.field_with_name(name) {
            Ok(field) => {
                if !types_compatible(field.data_type(), expected_type) {
                    return Err(DataFusionError::Plan(format!(
                        "Variation cache column '{}' has type {:?}, expected {:?}",
                        name,
                        field.data_type(),
                        expected_type,
                    )));
                }
            }
            Err(_) => {
                return Err(DataFusionError::Plan(format!(
                    "Variation cache table is missing required column '{}'. \
                     Required columns: {}",
                    name,
                    REQUIRED_VARIATION_COLUMNS
                        .iter()
                        .map(|(n, _)| *n)
                        .collect::<Vec<_>>()
                        .join(", ")
                )));
            }
        }
    }
    Ok(())
}

/// Parse a comma-separated column list string into a vector of column names.
///
/// Trims whitespace from each column name and filters out empty entries.
pub fn parse_column_list(column_list: &str) -> Vec<String> {
    column_list
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

/// Validate that all requested columns exist in the schema.
pub fn validate_requested_columns(schema: &SchemaRef, columns: &[String]) -> Result<()> {
    for col in columns {
        if schema.field_with_name(col).is_err() {
            let available: Vec<_> = schema.fields().iter().map(|f| f.name().as_str()).collect();
            return Err(DataFusionError::Plan(format!(
                "Requested column '{}' not found in cache table. Available columns: {}",
                col,
                available.join(", ")
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{Field, Schema};
    use std::sync::Arc;

    fn valid_variation_schema() -> SchemaRef {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
        ]))
    }

    #[test]
    fn validate_valid_schema() {
        let schema = valid_variation_schema();
        assert!(validate_variation_schema(&schema).is_ok());
    }

    #[test]
    fn validate_utf8view_accepted() {
        // DataFusion 50+ reads parquet strings as Utf8View
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8View, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8View, false),
            Field::new("allele_string", DataType::Utf8View, false),
        ]));
        assert!(validate_variation_schema(&schema).is_ok());
    }

    #[test]
    fn validate_missing_column() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
        ]));
        let result = validate_variation_schema(&schema);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("end"));
    }

    #[test]
    fn validate_wrong_type() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Utf8, false), // wrong type
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let result = validate_variation_schema(&schema);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("start"));
    }

    #[test]
    fn parse_column_list_basic() {
        let cols = parse_column_list("variation_name,clin_sig,gnomADg_AF");
        assert_eq!(cols, vec!["variation_name", "clin_sig", "gnomADg_AF"]);
    }

    #[test]
    fn parse_column_list_with_whitespace() {
        let cols = parse_column_list(" variation_name , clin_sig ");
        assert_eq!(cols, vec!["variation_name", "clin_sig"]);
    }

    #[test]
    fn parse_column_list_empty() {
        let cols = parse_column_list("");
        assert!(cols.is_empty());
    }

    #[test]
    fn validate_requested_columns_valid() {
        let schema = valid_variation_schema();
        let cols = vec!["chrom".to_string(), "clin_sig".to_string()];
        assert!(validate_requested_columns(&schema, &cols).is_ok());
    }

    #[test]
    fn validate_requested_columns_invalid() {
        let schema = valid_variation_schema();
        let cols = vec!["nonexistent".to_string()];
        let result = validate_requested_columns(&schema, &cols);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("nonexistent"));
    }
}
