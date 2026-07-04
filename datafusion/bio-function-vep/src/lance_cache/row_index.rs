/// The outcome of resolving a batch of requested variation positions to dataset
/// row ids: how many positions were requested, how many matched, and the matched
/// row ids (in request order). Produced by the Parquet variation lookup.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedRowIds {
    pub requested_positions: usize,
    pub matched_positions: usize,
    pub row_ids: Vec<u64>,
}
