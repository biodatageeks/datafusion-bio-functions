use std::collections::HashMap;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};

use crate::cache_source::CACHE_SOURCE_METADATA_KEY;

pub const VARIATION_FORBIDDEN_COLUMNS: &[&str] = &[
    "position_key",
    "variant_keys",
    "region_bin",
    "var_synonyms",
    "strand",
    "source_cache_path",
    "source_file",
];

pub const VARIATION_REQUIRED_COLUMNS: &[&str] = &[
    "chrom",
    "start",
    "end",
    "allele_string",
    "failed",
    "variation_name",
    "clin_sig",
    "clin_sig_allele",
    "clinical_impact",
    "phenotype_or_disease",
    "pubmed",
    "somatic",
    "minor_allele",
    "minor_allele_freq",
    "AF",
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_NFE",
    "gnomADe_SAS",
    "gnomADe_MID",
    "gnomADe_REMAINING",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_SAS",
    "gnomADg_REMAINING",
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
];

pub const LANCE_STRUCTURAL_ENCODING_KEY: &str = "lance-encoding:structural-encoding";
pub const LANCE_STRUCTURAL_ENCODING_MINIBLOCK: &str = "miniblock";
pub const LANCE_STRUCTURAL_ENCODING_FULLZIP: &str = "fullzip";

pub fn lance_field_metadata() -> HashMap<String, String> {
    [
        (
            LANCE_STRUCTURAL_ENCODING_KEY,
            LANCE_STRUCTURAL_ENCODING_MINIBLOCK,
        ),
        ("lance-encoding:compression", "zstd"),
        ("lance-encoding:compression-level", "3"),
        ("lance-encoding:dict-values-compression", "zstd"),
        ("lance-encoding:dict-values-compression-level", "3"),
        ("lance-encoding:rle-threshold", "0.95"),
        ("lance-encoding:dict-size-ratio", "0.99"),
        ("lance-encoding:dict-divisor", "1"),
        ("lance-encoding:minichunk-size", "4096"),
    ]
    .into_iter()
    .map(|(key, value)| (key.to_string(), value.to_string()))
    .collect()
}

pub fn validate_variation_schema(schema: &Schema) -> Result<()> {
    for name in VARIATION_FORBIDDEN_COLUMNS {
        if schema.index_of(name).is_ok() {
            return Err(DataFusionError::Execution(format!(
                "single-path Lance variation schema must not contain {name}"
            )));
        }
    }

    require_type(schema, "chrom", &DataType::Utf8)?;
    require_type(schema, "start", &DataType::UInt32)?;
    require_type(schema, "end", &DataType::UInt32)?;
    require_type(schema, "allele_string", &DataType::Utf8)?;
    require_type(schema, "failed", &DataType::Int8)?;
    require_type(schema, "variation_name", &DataType::Utf8)?;
    // `tier` is a build-time derived column (0 = warm/common, 1 = cold/rare) that
    // clusters frequently-queried rows into a dense prefix. It is not read from the
    // source table, so it is appended by the writer rather than listed in
    // VARIATION_REQUIRED_COLUMNS.
    require_type(schema, "tier", &DataType::Int8)?;
    Ok(())
}

fn require_type(schema: &Schema, name: &str, expected: &DataType) -> Result<()> {
    let field = schema.field_with_name(name).map_err(|_| {
        DataFusionError::Execution(format!("variation schema missing required field {name}"))
    })?;
    if field.data_type() != expected {
        return Err(DataFusionError::Execution(format!(
            "variation field {name} must be {expected:?}, got {:?}",
            field.data_type()
        )));
    }
    Ok(())
}

pub fn with_cache_source_metadata(schema: &Schema, source_type: &str) -> Schema {
    let mut metadata = schema.metadata().clone();
    metadata.insert(
        CACHE_SOURCE_METADATA_KEY.to_string(),
        source_type.to_string(),
    );
    Schema::new_with_metadata(schema.fields().clone(), metadata)
}

pub fn with_lance_field_metadata(schema: &Schema) -> Schema {
    let defaults = lance_field_metadata();
    let fields = schema
        .fields()
        .iter()
        .map(|field| {
            let mut metadata = field.metadata().clone();
            for (key, value) in &defaults {
                metadata.entry(key.clone()).or_insert_with(|| value.clone());
            }
            field.as_ref().clone().with_metadata(metadata)
        })
        .collect::<Vec<Field>>();
    Schema::new_with_metadata(fields, schema.metadata().clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    #[test]
    fn required_variation_columns_exclude_legacy_helpers() {
        assert!(VARIATION_REQUIRED_COLUMNS.contains(&"start"));
        assert!(VARIATION_REQUIRED_COLUMNS.contains(&"end"));
        assert!(VARIATION_REQUIRED_COLUMNS.contains(&"allele_string"));
        assert!(!VARIATION_REQUIRED_COLUMNS.contains(&"position_key"));
        assert!(!VARIATION_REQUIRED_COLUMNS.contains(&"variant_keys"));
        assert!(!VARIATION_REQUIRED_COLUMNS.contains(&"tier"));
    }

    #[test]
    fn variation_projection_requires_start_end_uint32() {
        let schema = Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]);
        validate_variation_schema(&schema).unwrap();
    }

    #[test]
    fn variation_projection_requires_int8_tier() {
        let schema = Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ]);
        let err = validate_variation_schema(&schema).unwrap_err().to_string();
        assert!(
            err.contains("tier"),
            "expected missing tier error, got {err}"
        );

        let wrong_type = Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
            Field::new("tier", DataType::UInt8, false),
        ]);
        let err = validate_variation_schema(&wrong_type)
            .unwrap_err()
            .to_string();
        assert!(err.contains("tier"), "expected tier type error, got {err}");
    }

    #[test]
    fn variation_projection_rejects_position_key() {
        let schema = Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
            Field::new("position_key", DataType::UInt64, false),
        ]);
        let err = validate_variation_schema(&schema).unwrap_err().to_string();
        assert!(err.contains("position_key"));
    }
}
