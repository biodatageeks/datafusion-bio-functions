use std::collections::HashMap;
use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};

// Input column names (matching datafusion-bio-format-bam)
pub const COL_CHROM: &str = "chrom";
pub const COL_START: &str = "start";
pub const COL_FLAGS: &str = "flags";
pub const COL_CIGAR: &str = "cigar";
pub const COL_MAPPING_QUALITY: &str = "mapping_quality";

// Output column names (matching SeQuiLa coverage schema)
pub const OUT_CONTIG: &str = "contig";
pub const OUT_POS_START: &str = "pos_start";
pub const OUT_POS_END: &str = "pos_end";
pub const OUT_COVERAGE: &str = "coverage";

/// Schema-level metadata key indicating whether coordinates are 0-based.
/// Same value as `datafusion-bio-format-core` uses.
pub const COORDINATE_SYSTEM_METADATA_KEY: &str = "bio.coordinate_system_zero_based";

/// Returns the output schema for coverage results.
///
/// The `zero_based` flag is embedded as schema metadata so downstream consumers
/// (e.g. polars-bio) can interpret the coordinate system.
pub fn coverage_output_schema(zero_based: bool) -> SchemaRef {
    let fields = vec![
        Field::new(OUT_CONTIG, DataType::Utf8, true),
        Field::new(OUT_POS_START, DataType::Int32, false),
        Field::new(OUT_POS_END, DataType::Int32, false),
        Field::new(OUT_COVERAGE, DataType::Int16, false),
    ];
    let mut metadata = HashMap::new();
    metadata.insert(
        COORDINATE_SYSTEM_METADATA_KEY.to_string(),
        zero_based.to_string(),
    );
    Arc::new(Schema::new(fields).with_metadata(metadata))
}
