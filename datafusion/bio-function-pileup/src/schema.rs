use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use std::sync::Arc;

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

/// Returns the output schema for coverage results.
pub fn coverage_output_schema() -> SchemaRef {
    Arc::new(Schema::new(vec![
        Field::new(OUT_CONTIG, DataType::Utf8, true),
        Field::new(OUT_POS_START, DataType::Int32, false),
        Field::new(OUT_POS_END, DataType::Int32, false),
        Field::new(OUT_COVERAGE, DataType::Int16, false),
    ]))
}
