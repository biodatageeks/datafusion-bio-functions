use datafusion::arrow::array::{
    Array, LargeStringArray, StringArray, StringViewArray, UInt32Array, UInt64Array,
};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use lance::dataset::index::LanceIndexStoreExt;
use lance::index::{DatasetIndexExt, DatasetIndexInternalExt};
use lance::table::format::IndexMetadata;
use lance_index::IndexCriteria;
use lance_index::scalar::{IndexStore, lance_format::LanceIndexStore};

const BTREE_PAGE_DATA_FILE: &str = "page_data.lance";
const BTREE_PAGE_LOOKUP_FILE: &str = "page_lookup.lance";
const BTREE_VALUES_COLUMN: &str = "values";
const BTREE_IDS_COLUMN: &str = "ids";
const PAGE_LOOKUP_MIN_COLUMN: &str = "min";
const PAGE_LOOKUP_PAGE_IDX_COLUMN: &str = "page_idx";
const BATCH_SIZE_META_KEY: &str = "batch_size";
const DEFAULT_BTREE_BATCH_SIZE: u64 = 4096;

/// Fully-materialized u32 position→row_id index. Retained only as the test
/// parity oracle for [`StreamingPositionCursor`]; production resolves positions
/// lazily via the streaming cursor and never builds this.
#[cfg(test)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PositionRowIdIndex {
    positions: Vec<u32>,
    row_ids: Vec<u64>,
    unique_positions: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringRowIdIndex {
    values: Vec<String>,
    row_ids: Vec<u64>,
    unique_values: usize,
}

/// In-memory `u64 key -> row_id` index built from a scalar BTree index's page
/// data. Used by the position-sliced SIFT store, where `key = (transcript_uid
/// << 32) | protein_position`. Keys are unique by construction (one row per
/// `(uid, position)`), so each key maps to a single row id.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct U64RowIdIndex {
    keys: Vec<u64>,
    row_ids: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedRowIds {
    pub requested_positions: usize,
    pub matched_positions: usize,
    pub row_ids: Vec<u64>,
}

#[cfg(test)]
impl PositionRowIdIndex {
    #[cfg(test)]
    pub fn from_pairs_for_test(mut pairs: Vec<(u32, u64)>) -> Self {
        pairs.sort_unstable_by_key(|pair| *pair);
        let mut positions = Vec::with_capacity(pairs.len());
        let mut row_ids = Vec::with_capacity(pairs.len());
        let mut unique_positions = 0;
        let mut last = None;
        for (position, row_id) in pairs {
            if last != Some(position) {
                unique_positions += 1;
                last = Some(position);
            }
            positions.push(position);
            row_ids.push(row_id);
        }
        Self {
            positions,
            row_ids,
            unique_positions,
        }
    }

    pub fn new(positions: Vec<u32>, row_ids: Vec<u64>, unique_positions: usize) -> Result<Self> {
        if positions.len() != row_ids.len() {
            return Err(DataFusionError::Execution(
                "position row-id index arrays have different lengths".into(),
            ));
        }
        if positions.windows(2).any(|pair| pair[0] > pair[1]) {
            return Err(DataFusionError::Execution(
                "position row-id index positions are not sorted".into(),
            ));
        }
        Ok(Self {
            positions,
            row_ids,
            unique_positions,
        })
    }

    pub fn unique_positions(&self) -> usize {
        self.unique_positions
    }

    pub fn row_ids_len(&self) -> usize {
        self.row_ids.len()
    }

    pub fn resolve_sorted_positions_from_cursor(
        &self,
        query_positions: &[u32],
        cursor: &mut usize,
    ) -> ResolvedRowIds {
        let mut matched_positions = 0;
        let mut row_ids = Vec::new();
        *cursor = (*cursor).min(self.positions.len());
        if let Some(first_position) = query_positions.first().copied() {
            let cursor_position = self.positions.get(*cursor).copied();
            if cursor_position.is_none_or(|position| position > first_position) {
                *cursor = self
                    .positions
                    .partition_point(|position| *position < first_position);
            }
        }

        for &position in query_positions {
            while *cursor < self.positions.len() && self.positions[*cursor] < position {
                *cursor += 1;
            }
            if *cursor == self.positions.len() {
                break;
            }
            if self.positions[*cursor] == position {
                matched_positions += 1;
                while *cursor < self.positions.len() && self.positions[*cursor] == position {
                    row_ids.push(self.row_ids[*cursor]);
                    *cursor += 1;
                }
            }
        }

        ResolvedRowIds {
            requested_positions: query_positions.len(),
            matched_positions,
            row_ids,
        }
    }
}

impl StringRowIdIndex {
    pub fn new(values: Vec<String>, row_ids: Vec<u64>, unique_values: usize) -> Result<Self> {
        if values.len() != row_ids.len() {
            return Err(DataFusionError::Execution(
                "string row-id index arrays have different lengths".into(),
            ));
        }
        if values.windows(2).any(|pair| pair[0] > pair[1]) {
            return Err(DataFusionError::Execution(
                "string row-id index values are not sorted".into(),
            ));
        }
        Ok(Self {
            values,
            row_ids,
            unique_values,
        })
    }

    pub fn unique_values(&self) -> usize {
        self.unique_values
    }

    pub fn row_ids_len(&self) -> usize {
        self.row_ids.len()
    }

    pub fn resolve_sorted_values(&self, query_values: &[String]) -> Vec<u64> {
        let mut cursor = 0usize;
        let mut row_ids = Vec::new();

        for value in query_values {
            while cursor < self.values.len() && self.values[cursor].as_str() < value.as_str() {
                cursor += 1;
            }
            if cursor == self.values.len() {
                break;
            }
            if self.values[cursor] == *value {
                while cursor < self.values.len() && self.values[cursor] == *value {
                    row_ids.push(self.row_ids[cursor]);
                    cursor += 1;
                }
            }
        }

        row_ids
    }
}

impl U64RowIdIndex {
    pub fn new(keys: Vec<u64>, row_ids: Vec<u64>) -> Result<Self> {
        if keys.len() != row_ids.len() {
            return Err(DataFusionError::Execution(
                "u64 row-id index arrays have different lengths".into(),
            ));
        }
        // Keys are unique by construction (one row per `(uid, position)`), so
        // require strictly increasing — this also rejects duplicates, which
        // `resolve()`'s binary search would otherwise silently collapse.
        if keys.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(DataFusionError::Execution(
                "u64 row-id index keys are not strictly increasing (unsorted or duplicate)".into(),
            ));
        }
        Ok(Self { keys, row_ids })
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Resolve queried keys to row ids via binary search. Keys absent from the
    /// index (a position with no predictions / not built) are silently skipped.
    /// Returns `(row_ids, present_keys)` aligned by index so the caller can map
    /// a fetched row back to the key it satisfied.
    pub fn resolve(&self, query_keys: &[u64]) -> (Vec<u64>, Vec<u64>) {
        let mut row_ids = Vec::with_capacity(query_keys.len());
        let mut present = Vec::with_capacity(query_keys.len());
        for &key in query_keys {
            if let Ok(idx) = self.keys.binary_search(&key) {
                row_ids.push(self.row_ids[idx]);
                present.push(key);
            }
        }
        (row_ids, present)
    }
}

/// In-memory page directory loaded from `page_lookup.lance` (the BTree
/// sub-index). Maps a u32 position to the page number whose `[min,max]` band may
/// contain it, so a streaming cursor can seek straight to its band without
/// scanning `page_data.lance`. Tiny (one record per 4096-row page).
pub(crate) struct PositionPageDirectory {
    /// `page_mins[i]` = minimum position stored in page `i`; ascending across
    /// pages because `page_data` is globally sorted by position.
    page_mins: Vec<u32>,
    num_pages: u32,
    batch_size: u64,
    num_rows: usize,
}

impl PositionPageDirectory {
    pub(crate) async fn load(dataset: &lance::Dataset) -> Result<Self> {
        let segments = load_btree_segments(dataset, "start", "start_btree_idx").await?;
        let store = LanceIndexStore::from_dataset_for_existing(dataset, &segments[0])
            .await
            .map_err(|e| DataFusionError::Execution(format!("open index store: {e}")))?;

        let lookup = store
            .open_index_file(BTREE_PAGE_LOOKUP_FILE)
            .await
            .map_err(|e| DataFusionError::Execution(format!("open page_lookup: {e}")))?;
        let batch_size = lookup
            .schema()
            .metadata
            .get(BATCH_SIZE_META_KEY)
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(DEFAULT_BTREE_BATCH_SIZE);
        let n = lookup.num_rows();
        let batch = lookup
            .read_range(0..n, None)
            .await
            .map_err(|e| DataFusionError::Execution(format!("read page_lookup: {e}")))?;
        let mins = batch
            .column_by_name(PAGE_LOOKUP_MIN_COLUMN)
            .ok_or_else(|| DataFusionError::Execution("page_lookup missing 'min'".into()))?
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| DataFusionError::Execution("page_lookup 'min' not u32".into()))?;
        let pages = batch
            .column_by_name(PAGE_LOOKUP_PAGE_IDX_COLUMN)
            .ok_or_else(|| DataFusionError::Execution("page_lookup missing 'page_idx'".into()))?
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| DataFusionError::Execution("page_lookup 'page_idx' not u32".into()))?;

        let num_pages = pages
            .values()
            .iter()
            .copied()
            .max()
            .map(|m| m + 1)
            .unwrap_or(0);
        let mut page_mins = vec![u32::MAX; num_pages as usize];
        for i in 0..mins.len() {
            // `min` is nullable in the schema; positions are non-null in practice.
            if mins.is_valid(i) {
                page_mins[pages.value(i) as usize] = mins.value(i);
            }
        }

        let data = store
            .open_index_file(BTREE_PAGE_DATA_FILE)
            .await
            .map_err(|e| DataFusionError::Execution(format!("open page_data: {e}")))?;

        Ok(Self {
            page_mins,
            num_pages,
            batch_size,
            num_rows: data.num_rows(),
        })
    }

    pub(crate) fn num_pages(&self) -> u32 {
        self.num_pages
    }

    pub(crate) fn batch_size(&self) -> u64 {
        self.batch_size
    }

    pub(crate) fn num_rows(&self) -> usize {
        self.num_rows
    }

    /// First page that may contain `position`: the last page whose `min <=
    /// position` (clamped to page 0 when `position` precedes all page minima).
    pub(crate) fn first_page_for(&self, position: u32) -> u32 {
        match self.page_mins.partition_point(|&m| m <= position) {
            0 => 0,
            i => (i - 1) as u32,
        }
    }
}

#[cfg(test)]
pub(crate) async fn load_start_index_from_lance_btree(
    dataset: &lance::Dataset,
) -> Result<PositionRowIdIndex> {
    load_u32_btree_index(dataset, "start", "start_btree_idx").await
}

pub async fn load_sift_key_index_from_lance_btree(
    dataset: &lance::Dataset,
) -> Result<U64RowIdIndex> {
    load_u64_btree_index(dataset, "key", "sift_key_btree_idx").await
}

pub async fn load_transcript_id_index_from_lance_btree(
    dataset: &lance::Dataset,
) -> Result<StringRowIdIndex> {
    load_string_btree_index(dataset, "transcript_id", "transcript_id_btree_idx").await
}

#[cfg(test)]
async fn load_u32_btree_index(
    dataset: &lance::Dataset,
    column: &str,
    index_name: &str,
) -> Result<PositionRowIdIndex> {
    let profile = std::env::var_os("VEP_LANCE_PROFILE").is_some();
    let segments_started = profile.then(std::time::Instant::now);
    let index_segments = load_btree_segments(dataset, column, index_name).await?;
    let segments_s = segments_started.map(|t| t.elapsed());
    let pairs_started = profile.then(std::time::Instant::now);
    let mut pairs = Vec::<(u32, u64)>::new();
    for index_segment in &index_segments {
        append_btree_segment_pairs(dataset, index_segment, &mut pairs).await?;
    }
    let pairs_read_s = pairs_started.map(|t| t.elapsed());
    let sort_started = profile.then(std::time::Instant::now);
    pairs.sort_unstable_by_key(|(position, row_id)| (*position, *row_id));
    let sort_s = sort_started.map(|t| t.elapsed());
    if let (Some(segments_s), Some(pairs_read_s), Some(sort_s)) = (segments_s, pairs_read_s, sort_s)
    {
        eprintln!(
            "[vep-lance-profile] u32_index_breakdown column={column} segments={} pairs={} segments_s={:.3} pairs_read_s={:.3} sort_s={:.3}",
            index_segments.len(),
            pairs.len(),
            segments_s.as_secs_f64(),
            pairs_read_s.as_secs_f64(),
            sort_s.as_secs_f64(),
        );
    }

    let mut positions = Vec::with_capacity(pairs.len());
    let mut row_ids = Vec::with_capacity(pairs.len());
    let mut unique_positions = 0usize;
    let mut previous = None;
    for (position, row_id) in pairs {
        if previous != Some(position) {
            unique_positions += 1;
            previous = Some(position);
        }
        positions.push(position);
        row_ids.push(row_id);
    }

    PositionRowIdIndex::new(positions, row_ids, unique_positions)
}

async fn load_u64_btree_index(
    dataset: &lance::Dataset,
    column: &str,
    index_name: &str,
) -> Result<U64RowIdIndex> {
    let index_segments = load_btree_segments(dataset, column, index_name).await?;
    let mut pairs = Vec::<(u64, u64)>::new();
    for index_segment in &index_segments {
        append_btree_segment_u64_pairs(dataset, index_segment, &mut pairs).await?;
    }
    pairs.sort_unstable_by_key(|(key, row_id)| (*key, *row_id));

    let mut keys = Vec::with_capacity(pairs.len());
    let mut row_ids = Vec::with_capacity(pairs.len());
    for (key, row_id) in pairs {
        keys.push(key);
        row_ids.push(row_id);
    }

    U64RowIdIndex::new(keys, row_ids)
}

async fn load_string_btree_index(
    dataset: &lance::Dataset,
    column: &str,
    index_name: &str,
) -> Result<StringRowIdIndex> {
    let index_segments = load_btree_segments(dataset, column, index_name).await?;
    let mut pairs = Vec::<(String, u64)>::new();
    for index_segment in &index_segments {
        append_btree_segment_string_pairs(dataset, index_segment, &mut pairs).await?;
    }
    pairs.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));

    let mut values = Vec::with_capacity(pairs.len());
    let mut row_ids = Vec::with_capacity(pairs.len());
    let mut unique_values = 0usize;
    let mut previous: Option<String> = None;
    for (value, row_id) in pairs {
        if previous.as_deref() != Some(value.as_str()) {
            unique_values += 1;
        }
        previous = Some(value.clone());
        values.push(value);
        row_ids.push(row_id);
    }

    StringRowIdIndex::new(values, row_ids, unique_values)
}

pub(crate) async fn load_btree_segments(
    dataset: &lance::Dataset,
    column: &str,
    index_name: &str,
) -> Result<Vec<IndexMetadata>> {
    let mut index_segments = dataset
        .load_indices_by_name(index_name)
        .await
        .map_err(|err| DataFusionError::Execution(format!("failed to load {index_name}: {err}")))?;
    if index_segments.is_empty() {
        if let Some(index) = dataset
            .load_scalar_index(
                IndexCriteria::default()
                    .for_column(column)
                    .supports_exact_equality(),
            )
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!(
                    "failed to discover BTree index for {column}: {err}"
                ))
            })?
        {
            index_segments.push(index);
        }
    }
    if index_segments.is_empty() {
        return Err(DataFusionError::Execution(format!(
            "missing scalar BTree index '{index_name}' for column '{column}'"
        )));
    }
    Ok(index_segments)
}

#[cfg(test)]
async fn append_btree_segment_pairs(
    dataset: &lance::Dataset,
    index_segment: &IndexMetadata,
    pairs: &mut Vec<(u32, u64)>,
) -> Result<()> {
    let store = LanceIndexStore::from_dataset_for_existing(dataset, index_segment)
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to open Lance index store: {err}"))
        })?;
    let reader = store
        .open_index_file(BTREE_PAGE_DATA_FILE)
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to open BTree page data: {err}"))
        })?;
    let num_rows = reader.num_rows();
    if num_rows == 0 {
        return Ok(());
    }
    pairs.reserve(num_rows);
    let mut stream = reader
        .read_range_stream(0..num_rows, Some(&[BTREE_VALUES_COLUMN, BTREE_IDS_COLUMN]))
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to stream BTree page data: {err}"))
        })?;
    while let Some(batch) = stream.try_next().await.map_err(|err| {
        DataFusionError::Execution(format!("failed to read BTree page data batch: {err}"))
    })? {
        append_btree_page_row_id_pairs(&batch, pairs)?;
    }
    Ok(())
}

async fn append_btree_segment_u64_pairs(
    dataset: &lance::Dataset,
    index_segment: &IndexMetadata,
    pairs: &mut Vec<(u64, u64)>,
) -> Result<()> {
    let store = LanceIndexStore::from_dataset_for_existing(dataset, index_segment)
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to open Lance index store: {err}"))
        })?;
    let reader = store
        .open_index_file(BTREE_PAGE_DATA_FILE)
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to open BTree page data: {err}"))
        })?;
    let num_rows = reader.num_rows();
    if num_rows == 0 {
        return Ok(());
    }
    pairs.reserve(num_rows);
    let mut stream = reader
        .read_range_stream(0..num_rows, Some(&[BTREE_VALUES_COLUMN, BTREE_IDS_COLUMN]))
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to stream BTree page data: {err}"))
        })?;
    while let Some(batch) = stream.try_next().await.map_err(|err| {
        DataFusionError::Execution(format!("failed to read BTree page data batch: {err}"))
    })? {
        append_btree_page_u64_row_id_pairs(&batch, pairs)?;
    }
    Ok(())
}

async fn append_btree_segment_string_pairs(
    dataset: &lance::Dataset,
    index_segment: &IndexMetadata,
    pairs: &mut Vec<(String, u64)>,
) -> Result<()> {
    let store = LanceIndexStore::from_dataset_for_existing(dataset, index_segment)
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to open Lance index store: {err}"))
        })?;
    let reader = store
        .open_index_file(BTREE_PAGE_DATA_FILE)
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to open BTree page data: {err}"))
        })?;
    let num_rows = reader.num_rows();
    if num_rows == 0 {
        return Ok(());
    }
    pairs.reserve(num_rows);
    let mut stream = reader
        .read_range_stream(0..num_rows, Some(&[BTREE_VALUES_COLUMN, BTREE_IDS_COLUMN]))
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!("failed to stream BTree page data: {err}"))
        })?;
    while let Some(batch) = stream.try_next().await.map_err(|err| {
        DataFusionError::Execution(format!("failed to read BTree page data batch: {err}"))
    })? {
        append_btree_page_string_row_id_pairs(&batch, pairs)?;
    }
    Ok(())
}

#[cfg(test)]
fn append_btree_page_row_id_pairs(batch: &RecordBatch, pairs: &mut Vec<(u32, u64)>) -> Result<()> {
    let value_array = batch
        .column_by_name(BTREE_VALUES_COLUMN)
        .unwrap_or_else(|| batch.column(0));
    let row_ids = batch
        .column_by_name(BTREE_IDS_COLUMN)
        .ok_or_else(|| DataFusionError::Execution("BTree page data missing ids column".into()))?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| DataFusionError::Execution("BTree page ids column must be UInt64".into()))?;
    let positions = value_array
        .as_any()
        .downcast_ref::<UInt32Array>()
        .ok_or_else(|| {
            DataFusionError::Execution("BTree page values column must be UInt32".into())
        })?;

    for row in 0..batch.num_rows() {
        if !positions.is_null(row) {
            pairs.push((positions.value(row), row_ids.value(row)));
        }
    }
    Ok(())
}

fn append_btree_page_u64_row_id_pairs(
    batch: &RecordBatch,
    pairs: &mut Vec<(u64, u64)>,
) -> Result<()> {
    let value_array = batch
        .column_by_name(BTREE_VALUES_COLUMN)
        .unwrap_or_else(|| batch.column(0));
    let row_ids = batch
        .column_by_name(BTREE_IDS_COLUMN)
        .ok_or_else(|| DataFusionError::Execution("BTree page data missing ids column".into()))?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| DataFusionError::Execution("BTree page ids column must be UInt64".into()))?;
    let keys = value_array
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| {
            DataFusionError::Execution("BTree page values column must be UInt64".into())
        })?;

    for row in 0..batch.num_rows() {
        if !keys.is_null(row) {
            pairs.push((keys.value(row), row_ids.value(row)));
        }
    }
    Ok(())
}

fn append_btree_page_string_row_id_pairs(
    batch: &RecordBatch,
    pairs: &mut Vec<(String, u64)>,
) -> Result<()> {
    let value_array = batch
        .column_by_name(BTREE_VALUES_COLUMN)
        .unwrap_or_else(|| batch.column(0));
    let row_ids = batch
        .column_by_name(BTREE_IDS_COLUMN)
        .ok_or_else(|| DataFusionError::Execution("BTree page data missing ids column".into()))?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| DataFusionError::Execution("BTree page ids column must be UInt64".into()))?;

    for row in 0..batch.num_rows() {
        if let Some(value) = string_value_at(value_array.as_ref(), row) {
            pairs.push((value, row_ids.value(row)));
        }
    }
    Ok(())
}

fn string_value_at(array: &dyn Array, row: usize) -> Option<String> {
    if array.is_null(row) {
        return None;
    }
    if let Some(values) = array.as_any().downcast_ref::<StringArray>() {
        return Some(values.value(row).to_string());
    }
    if let Some(values) = array.as_any().downcast_ref::<LargeStringArray>() {
        return Some(values.value(row).to_string());
    }
    if let Some(values) = array.as_any().downcast_ref::<StringViewArray>() {
        return Some(values.value(row).to_string());
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::ArrayRef;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    #[test]
    fn resolves_sorted_positions_from_monotonic_cursor() {
        let index = PositionRowIdIndex::from_pairs_for_test(vec![
            (10, 1),
            (10, 2),
            (20, 3),
            (30, 4),
            (30, 5),
        ]);
        let mut cursor = 0;
        let result = index.resolve_sorted_positions_from_cursor(&[10, 30], &mut cursor);
        assert_eq!(result.matched_positions, 2);
        assert_eq!(result.row_ids, vec![1, 2, 4, 5]);
        assert_eq!(cursor, 5);
    }

    #[test]
    fn skips_missing_positions_without_resetting_cursor() {
        let index = PositionRowIdIndex::from_pairs_for_test(vec![(10, 1), (20, 2), (40, 3)]);
        let mut cursor = 0;
        let result = index.resolve_sorted_positions_from_cursor(&[15, 20, 35], &mut cursor);
        assert_eq!(result.matched_positions, 1);
        assert_eq!(result.row_ids, vec![2]);
        assert_eq!(cursor, 2);
    }

    #[test]
    fn rewinds_cursor_when_extended_probes_overlap_later_window() {
        let index =
            PositionRowIdIndex::from_pairs_for_test(vec![(100, 1), (101, 2), (102, 3), (130, 4)]);
        let mut cursor = 0;

        let first = index.resolve_sorted_positions_from_cursor(&[100, 102, 130], &mut cursor);
        assert_eq!(first.row_ids, vec![1, 3, 4]);
        assert_eq!(cursor, 4);

        let second = index.resolve_sorted_positions_from_cursor(&[101, 102], &mut cursor);
        assert_eq!(second.row_ids, vec![2, 3]);
        assert_eq!(cursor, 3);
    }

    #[test]
    fn extracts_u32_position_row_id_pairs_from_btree_page_data() {
        let schema = Arc::new(Schema::new(vec![
            Field::new(BTREE_VALUES_COLUMN, DataType::UInt32, false),
            Field::new(BTREE_IDS_COLUMN, DataType::UInt64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(vec![10, 10, 30])) as ArrayRef,
                Arc::new(UInt64Array::from(vec![101, 102, 300])) as ArrayRef,
            ],
        )
        .unwrap();
        let mut pairs = Vec::new();

        append_btree_page_row_id_pairs(&batch, &mut pairs).unwrap();

        assert_eq!(pairs, vec![(10, 101), (10, 102), (30, 300)]);
    }

    #[test]
    fn extracts_u64_key_row_id_pairs_from_btree_page_data() {
        let schema = Arc::new(Schema::new(vec![
            Field::new(BTREE_VALUES_COLUMN, DataType::UInt64, false),
            Field::new(BTREE_IDS_COLUMN, DataType::UInt64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![
                    1u64 << 32,
                    (1u64 << 32) | 5,
                    2u64 << 32,
                ])) as ArrayRef,
                Arc::new(UInt64Array::from(vec![0u64, 1, 2])) as ArrayRef,
            ],
        )
        .unwrap();
        let mut pairs = Vec::new();

        append_btree_page_u64_row_id_pairs(&batch, &mut pairs).unwrap();

        assert_eq!(
            pairs,
            vec![(1u64 << 32, 0), ((1u64 << 32) | 5, 1), (2u64 << 32, 2)]
        );
    }

    #[test]
    fn u64_index_resolves_present_keys_and_skips_missing() {
        let index = U64RowIdIndex::new(vec![10, 20, 40], vec![100, 200, 400]).unwrap();
        let (row_ids, present) = index.resolve(&[10, 30, 40, 50]);
        assert_eq!(row_ids, vec![100, 400]);
        assert_eq!(present, vec![10, 40]);
    }

    #[test]
    fn u64_index_rejects_duplicate_keys() {
        // Uniqueness is invariant by construction; duplicates would be silently
        // collapsed by binary search, so the constructor must reject them.
        let err = U64RowIdIndex::new(vec![10, 10, 20], vec![100, 101, 200]).unwrap_err();
        assert!(err.to_string().contains("strictly increasing"));
    }
}
