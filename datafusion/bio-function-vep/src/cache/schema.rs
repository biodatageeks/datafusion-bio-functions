use std::collections::{HashMap, HashSet};

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};

use crate::cache_identity::CACHE_VERSION_METADATA_KEY;
use crate::cache_source::CACHE_SOURCE_METADATA_KEY;

/// Project a source variation Arrow schema to the cache's variation output
/// schema: keep [`VARIATION_REQUIRED_COLUMNS`] plus any present
/// [`VARIATION_OPTIONAL_COLUMNS`] (dropping [`VARIATION_FORBIDDEN_COLUMNS`]),
/// widen `start`/`end` to `UInt32`, and append the derived `tier` column. Used
/// by both the Parquet shard writer and the build driver, so it lives here
/// (under the read-runtime `parquet-cache` feature) rather than in the builder
/// module.
pub(crate) fn variation_projected_schema(
    source_schema: &Schema,
    source_type: &str,
) -> Result<Schema> {
    let mut fields = Vec::new();
    let forbidden = VARIATION_FORBIDDEN_COLUMNS
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    for name in VARIATION_REQUIRED_COLUMNS {
        if forbidden.contains(name) {
            continue;
        }
        let (_, field) = source_schema.column_with_name(name).ok_or_else(|| {
            DataFusionError::Execution(format!("variation source batch missing column {name}"))
        })?;
        if *name == "start" || *name == "end" {
            fields.push(Field::new(*name, DataType::UInt32, field.is_nullable()));
        } else {
            fields.push(field.as_ref().clone());
        }
        if *name == "clin_sig_allele" {
            for optional_name in VARIATION_OPTIONAL_COLUMNS {
                if let Some((_, optional_field)) = source_schema.column_with_name(optional_name) {
                    if optional_field.data_type() != &DataType::Utf8 {
                        return Err(DataFusionError::Execution(format!(
                            "optional variation field {optional_name} must be Utf8, got {:?}",
                            optional_field.data_type()
                        )));
                    }
                    fields.push(Field::new(*optional_name, DataType::Utf8, true));
                }
            }
        }
    }

    // Derived warm/cold tier column (0 = warm/common, 1 = cold/rare). Appended
    // here rather than read from the source table.
    fields.push(Field::new("tier", DataType::Int8, false));

    let projected = Schema::new_with_metadata(fields, source_schema.metadata().clone());
    let target_schema = with_cache_source_metadata(&projected, source_type);
    validate_variation_schema(&target_schema)?;
    Ok(target_schema)
}

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

/// Release-specific variation fields that are preserved when present but are
/// not required. VEP cache 116 added `clin_sig_ref_allele`; cache 115 does not
/// contain it and must retain its existing physical schema and behavior.
pub const VARIATION_OPTIONAL_COLUMNS: &[&str] = &["clin_sig_ref_allele"];

pub fn validate_variation_schema(schema: &Schema) -> Result<()> {
    for name in VARIATION_FORBIDDEN_COLUMNS {
        if schema.index_of(name).is_ok() {
            return Err(DataFusionError::Execution(format!(
                "single-path variation schema must not contain {name}"
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

pub fn with_cache_identity_metadata(
    schema: &Schema,
    source_type: &str,
    cache_version: &str,
) -> Schema {
    let schema = with_cache_source_metadata(schema, source_type);
    let mut metadata = schema.metadata().clone();
    metadata.insert(
        CACHE_VERSION_METADATA_KEY.to_string(),
        cache_version.to_string(),
    );
    Schema::new_with_metadata(schema.fields().clone(), metadata)
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
        assert!(!VARIATION_REQUIRED_COLUMNS.contains(&"clin_sig_ref_allele"));
        assert_eq!(VARIATION_OPTIONAL_COLUMNS, &["clin_sig_ref_allele"]);
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
