use datafusion::arrow::record_batch::RecordBatch;

use crate::cache::row_index::ResolvedRowIds;
use crate::cache::schema::VARIATION_FORBIDDEN_COLUMNS;

/// The rows taken for a set of resolved positions: the resolution bookkeeping
/// (`resolved`) plus the projected payload `batch`. Produced by the Parquet
/// variation lookup and consumed by the read seam.
#[derive(Debug)]
pub struct TakenVariationRows {
    pub resolved: ResolvedRowIds,
    pub batch: RecordBatch,
}

/// Sanitize a requested variation projection: drop build-only/forbidden columns
/// (notably `tier`, the build clustering column that must never be materialized
/// into annotation output) and always include the coordinate/matcher columns the
/// read seam needs.
pub fn ensure_runtime_projection(projection: Vec<String>) -> Vec<String> {
    let mut sanitized = Vec::with_capacity(projection.len() + 4);
    for column in projection {
        // `tier` is a build-only clustering column: it is stored in the dataset
        // (so it is no longer in VARIATION_FORBIDDEN_COLUMNS) but must never be
        // materialized into annotation output.
        let excluded = column == "tier"
            || VARIATION_FORBIDDEN_COLUMNS
                .iter()
                .any(|forbidden| column == *forbidden);
        if !excluded && !sanitized.iter().any(|existing| existing == &column) {
            sanitized.push(column);
        }
    }
    for required in ["start", "end", "allele_string", "failed"] {
        if !sanitized.iter().any(|column| column == required) {
            sanitized.push(required.to_string());
        }
    }
    sanitized
}
