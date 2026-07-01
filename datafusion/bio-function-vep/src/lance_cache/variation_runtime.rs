use std::collections::VecDeque;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{Array, UInt32Array, UInt64Array};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use lance::dataset::ProjectionRequest;
use lance::dataset::index::LanceIndexStoreExt;
use lance_index::scalar::lance_format::LanceIndexStore;
use lance_index::scalar::{IndexReader, IndexStore};

use crate::lance_cache::row_index::{PositionPageDirectory, ResolvedRowIds, load_btree_segments};
use crate::lance_cache::schema::VARIATION_FORBIDDEN_COLUMNS;

#[derive(Debug)]
pub struct TakenVariationRows {
    pub resolved: ResolvedRowIds,
    pub batch: RecordBatch,
}

/// Per-contig variation lookup: a tiny in-memory page directory over the BTree
/// `start` index plus a streaming `page_data` reader. Position→row_id resolution
/// is performed lazily, per partition, by a [`StreamingPositionCursor`] — the
/// index is never fully materialized.
pub struct SinglePathLanceVariationLookup {
    dataset: lance::Dataset,
    projection: Vec<String>,
    reader: Arc<dyn IndexReader>,
    directory: Arc<PositionPageDirectory>,
}

impl SinglePathLanceVariationLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let profile = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        let dataset_started = profile.then(Instant::now);
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!(
                    "failed to open Lance variation dataset '{}': {err}",
                    path.display()
                ))
            })?;
        let dataset_open = dataset_started.map(|t| t.elapsed());

        let load_started = profile.then(Instant::now);
        let directory = Arc::new(PositionPageDirectory::load(&dataset).await?);
        let reader = open_page_data_reader(&dataset).await?;
        if let (Some(dataset_open), Some(load_started)) = (dataset_open, load_started) {
            eprintln!(
                "[vep-lance-profile] open_breakdown dataset_open_s={:.3} directory_load_s={:.3} rows={} pages={}",
                dataset_open.as_secs_f64(),
                load_started.elapsed().as_secs_f64(),
                directory.num_rows(),
                directory.num_pages(),
            );
        }

        Ok(Self {
            dataset,
            projection: ensure_runtime_projection(projection),
            reader,
            directory,
        })
    }

    /// Create a fresh per-partition streaming cursor.
    pub fn new_cursor(&self) -> StreamingPositionCursor {
        StreamingPositionCursor::new(self.reader.clone(), self.directory.clone())
    }

    pub async fn resolve_and_take(
        &self,
        sorted_unique_starts: &[u32],
        cursor: &mut StreamingPositionCursor,
    ) -> Result<TakenVariationRows> {
        let profile_enabled = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        let resolved = cursor.resolve(sorted_unique_starts).await?;
        if profile_enabled {
            eprintln!(
                "[vep-lance-profile] variation_resolve requested_starts={} matched_positions={} row_ids={}",
                resolved.requested_positions,
                resolved.matched_positions,
                resolved.row_ids.len(),
            );
        }
        // The cache stores the 27 AF columns bundled into 3 fullzip List<Utf8> columns.
        // Project the bundled columns for the take, then expand them back so downstream
        // sees the 27 logical AF columns unchanged.
        let physical_projection =
            crate::lance_cache::af_bundle::bundle_projection(&self.projection);
        let projection_request = ProjectionRequest::from_columns(
            physical_projection.iter().map(String::as_str),
            self.dataset.schema(),
        );
        let take_started = profile_enabled.then(Instant::now);
        let batch = self
            .dataset
            .take_rows(&resolved.row_ids, projection_request)
            .await
            .map_err(|err| DataFusionError::Execution(format!("Lance take_rows failed: {err}")))?;
        let batch = crate::lance_cache::af_bundle::unbundle_af_columns(&batch)?;
        if let Some(started) = take_started {
            eprintln!(
                "[vep-lance-profile] variation_take row_ids={} batch_rows={} seconds={:.3}",
                resolved.row_ids.len(),
                batch.num_rows(),
                started.elapsed().as_secs_f64(),
            );
        }
        Ok(TakenVariationRows { resolved, batch })
    }

    pub fn projection(&self) -> &[String] {
        &self.projection
    }

    pub fn row_ids_len(&self) -> usize {
        self.directory.num_rows()
    }

    pub fn unique_positions(&self) -> usize {
        // Streaming never materializes unique positions; report total rows as an
        // upper bound for the informational profile line.
        self.directory.num_rows()
    }
}

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

/// Number of pages retained in the streaming cursor's sliding window: the page
/// holding the cursor plus ±2 neighbours. Caps resident index RAM per partition
/// at `WINDOW_PAGES × batch_size × 12B` (~240 KB at 4096-row pages).
pub(crate) const WINDOW_PAGES: usize = 5;

/// One decoded BTree leaf page: parallel `(position, row_id)` columns, sorted.
struct DecodedPage {
    positions: Vec<u32>,
    row_ids: Vec<u64>,
}

/// Open the BTree `page_data.lance` leaf reader for the variation dataset's
/// `start` index. Cheap — no data is read until a page is requested.
pub(crate) async fn open_page_data_reader(
    dataset: &lance::Dataset,
) -> Result<Arc<dyn IndexReader>> {
    let segments = load_btree_segments(dataset, "start", "start_btree_idx").await?;
    let store = LanceIndexStore::from_dataset_for_existing(dataset, &segments[0])
        .await
        .map_err(|e| DataFusionError::Execution(format!("open index store: {e}")))?;
    store
        .open_index_file("page_data.lance")
        .await
        .map_err(|e| DataFusionError::Execution(format!("open page_data: {e}")))
}

/// Per-partition forward streaming cursor over `page_data`. Seeks selectively to
/// its band via the page directory, then walks pages with a bounded sliding
/// window (cursor page ± 2), falling back to a `page_lookup` re-read for any
/// backstep beyond the window. Resident RAM is `O(WINDOW_PAGES)`, not the index.
///
/// The per-partition resolution cursor used by `KvLookupStream`; all
/// constructors/methods are crate-internal.
pub struct StreamingPositionCursor {
    reader: Arc<dyn IndexReader>,
    directory: Arc<PositionPageDirectory>,
    /// Sliding window of decoded pages: up to 2 trailing pages (for backsteps),
    /// the cursor page, and any loaded ahead — capped at `WINDOW_PAGES`.
    window: VecDeque<DecodedPage>,
    /// Index within `window` of the page holding the cursor.
    cursor_page: usize,
    /// Offset of the next unconsumed row within `window[cursor_page]`.
    off: usize,
    /// Next page number to pull when the cursor advances past loaded pages.
    next_page_to_load: u32,
    /// First page actually loaded (the band start chosen by the seek); `None`
    /// until the first `resolve`. Exposed for tests/observability.
    first_loaded_page: Option<u32>,
    seeded: bool,
}

impl StreamingPositionCursor {
    pub(crate) fn new(reader: Arc<dyn IndexReader>, directory: Arc<PositionPageDirectory>) -> Self {
        Self {
            reader,
            directory,
            window: VecDeque::new(),
            cursor_page: 0,
            off: 0,
            next_page_to_load: 0,
            first_loaded_page: None,
            seeded: false,
        }
    }

    pub(crate) fn window_len(&self) -> usize {
        self.window.len()
    }

    pub(crate) fn first_loaded_page(&self) -> Option<u32> {
        self.first_loaded_page
    }

    async fn load_page(&self, page_number: u32) -> Result<DecodedPage> {
        let batch = self
            .reader
            .read_record_batch(page_number as u64, self.directory.batch_size())
            .await
            .map_err(|e| DataFusionError::Execution(format!("read page {page_number}: {e}")))?;
        let positions = batch
            .column_by_name("values")
            .ok_or_else(|| DataFusionError::Execution("page_data missing 'values'".into()))?
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| DataFusionError::Execution("page_data 'values' not u32".into()))?
            .values()
            .to_vec();
        let row_ids = batch
            .column_by_name("ids")
            .ok_or_else(|| DataFusionError::Execution("page_data missing 'ids'".into()))?
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| DataFusionError::Execution("page_data 'ids' not u64".into()))?
            .values()
            .to_vec();
        Ok(DecodedPage { positions, row_ids })
    }

    fn cur_pos(&self) -> Option<u32> {
        self.window
            .get(self.cursor_page)
            .and_then(|p| p.positions.get(self.off).copied())
    }

    fn cur_id(&self) -> Option<u64> {
        self.window
            .get(self.cursor_page)
            .and_then(|p| p.row_ids.get(self.off).copied())
    }

    /// Drop pages so at most 2 remain *behind* the cursor page (the backstep
    /// margin) and the window never exceeds `WINDOW_PAGES`.
    fn evict_trailing(&mut self) {
        while self.cursor_page > 2 {
            self.window.pop_front();
            self.cursor_page -= 1;
        }
        while self.window.len() > WINDOW_PAGES && self.cursor_page > 0 {
            self.window.pop_front();
            self.cursor_page -= 1;
        }
    }

    /// Ensure the cursor page has an unconsumed row, moving to / loading the next
    /// page when the current one is exhausted. No-op once all pages are consumed.
    async fn ensure_cursor(&mut self) -> Result<()> {
        loop {
            match self.window.get(self.cursor_page) {
                Some(p) if self.off < p.positions.len() => return Ok(()),
                Some(_) => {
                    self.cursor_page += 1;
                    self.off = 0;
                }
                None => {
                    if self.next_page_to_load >= self.directory.num_pages() {
                        return Ok(());
                    }
                    let page = self.load_page(self.next_page_to_load).await?;
                    self.next_page_to_load += 1;
                    self.window.push_back(page);
                    self.evict_trailing();
                }
            }
        }
    }

    async fn advance(&mut self) -> Result<()> {
        self.off += 1;
        self.ensure_cursor().await
    }

    /// Scan all retained window pages for `q`, appending its row_ids. Used for a
    /// small backward step (extended probe) that lands in a trailing page.
    fn try_match_in_window(&self, q: u32, out: &mut Vec<u64>) -> bool {
        let mut found = false;
        for page in &self.window {
            let mut i = page.positions.partition_point(|&p| p < q);
            while i < page.positions.len() && page.positions[i] == q {
                out.push(page.row_ids[i]);
                i += 1;
                found = true;
            }
        }
        found
    }

    /// True if `q` precedes the smallest position currently retained — a
    /// backstep beyond the window that needs a `page_lookup` re-read.
    fn below_window(&self, q: u32) -> bool {
        self.window
            .front()
            .and_then(|p| p.positions.first().copied())
            .is_some_and(|min| q < min)
    }

    /// Correctness floor: resolve `q` by seeking its page via `page_lookup` and
    /// reading just that page (plus the next if `q`'s run spills over a
    /// boundary), without disturbing the streaming cursor. Rare — logged under
    /// `VEP_LANCE_PROFILE` so the window size can be tuned if it fires often.
    async fn reread_match(&self, q: u32, out: &mut Vec<u64>) -> Result<bool> {
        if std::env::var_os("VEP_LANCE_PROFILE").is_some() {
            eprintln!("[vep-lance-profile] streaming_reread position={q}");
        }
        let mut page_no = self.directory.first_page_for(q);
        let mut found = false;
        loop {
            let page = self.load_page(page_no).await?;
            let mut i = page.positions.partition_point(|&p| p < q);
            let mut hit_end = false;
            while i < page.positions.len() && page.positions[i] == q {
                out.push(page.row_ids[i]);
                i += 1;
                found = true;
                hit_end = i == page.positions.len();
            }
            // `q`'s run reached the page end → it may continue into the next page.
            if hit_end && (page_no as usize + 1) < self.directory.num_pages() as usize {
                page_no += 1;
                continue;
            }
            break;
        }
        Ok(found)
    }

    /// Resolve a batch of position queries (ascending) to row_ids, advancing the
    /// cursor forward. Within a matched position, row_ids are emitted ascending
    /// to match `load_u32_btree_index`'s `(position, row_id)` sort.
    pub(crate) async fn resolve(&mut self, sorted_positions: &[u32]) -> Result<ResolvedRowIds> {
        if !self.seeded {
            // Selective seek: jump straight to the band covering the first query
            // position instead of reading from page 0.
            if let Some(&first) = sorted_positions.first() {
                self.next_page_to_load = self.directory.first_page_for(first);
                self.first_loaded_page = Some(self.next_page_to_load);
            }
            self.seeded = true;
        }
        self.ensure_cursor().await?;
        let mut matched_positions = 0usize;
        let mut row_ids = Vec::new();
        for &q in sorted_positions {
            // Backward step (extended probe): the cursor is already past q. Match
            // against retained trailing pages; if q is below the window, fall
            // back to a selective page_lookup re-read (correctness floor).
            if self.cur_pos().is_some_and(|p| p > q) {
                let mut at_pos: Vec<u64> = Vec::new();
                let mut matched = self.try_match_in_window(q, &mut at_pos);
                if !matched && self.below_window(q) {
                    matched = self.reread_match(q, &mut at_pos).await?;
                }
                if matched {
                    matched_positions += 1;
                    at_pos.sort_unstable();
                    row_ids.extend(at_pos);
                }
                continue;
            }
            // Forward path.
            while self.cur_pos().is_some_and(|p| p < q) {
                self.advance().await?;
            }
            if self.cur_pos().is_none() {
                break;
            }
            if self.cur_pos() == Some(q) {
                matched_positions += 1;
                let mut at_pos: Vec<u64> = Vec::new();
                while self.cur_pos() == Some(q) {
                    if let Some(id) = self.cur_id() {
                        at_pos.push(id);
                    }
                    self.advance().await?;
                }
                at_pos.sort_unstable();
                row_ids.extend(at_pos);
            }
        }
        Ok(ResolvedRowIds {
            requested_positions: sorted_positions.len(),
            matched_positions,
            row_ids,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lance_cache::row_index::load_start_index_from_lance_btree;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    /// Build + write a single-contig indexed dataset with `n` rows at positions
    /// `base + i*step`, plus any extra `(pos, name)` rows prepended (for
    /// multi-allelic positions). Returns the dataset path's parent tempdir.
    async fn write_indexed_dataset(
        path: &std::path::Path,
        extra: &[(u32, &str)],
        n: u32,
        base: u32,
        step: u32,
    ) {
        let owned: Vec<String> = (0..n).map(|i| format!("rs{i}")).collect();
        let mut rows: Vec<(&str, u32, u32, &str, i8, &str)> =
            Vec::with_capacity(extra.len() + n as usize);
        for (pos, name) in extra {
            rows.push(("chr1", *pos, *pos, "A/T", 0, name));
        }
        for (i, name) in owned.iter().enumerate() {
            let p = base + (i as u32) * step;
            rows.push(("chr1", p, p, "A/T", 0, name.as_str()));
        }
        crate::lance_cache::write::write_record_batches_to_lance(
            path,
            vec![test_variation_batch(rows)],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn streaming_open_resolve_and_take_matches_eager() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        let batch = test_variation_batch(vec![
            ("chr1", 10, 10, "A/T", 0, "rs10"),
            ("chr1", 20, 20, "C/G", 0, "rs20"),
            ("chr1", 30, 30, "G/A", 0, "rs30"),
        ]);
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();
        let projection = vec![
            "start".to_string(),
            "end".to_string(),
            "allele_string".to_string(),
            "failed".to_string(),
            "variation_name".to_string(),
        ];

        // Streaming lookup (the only production path).
        let s = SinglePathLanceVariationLookup::open(&path, projection)
            .await
            .unwrap();
        let mut sc = s.new_cursor();
        let st = s.resolve_and_take(&[10, 20, 30], &mut sc).await.unwrap();

        // Cross-check row_ids against the materialized index oracle.
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let materialized = load_start_index_from_lance_btree(&dataset).await.unwrap();
        let mut oc = 0usize;
        let oracle = materialized.resolve_sorted_positions_from_cursor(&[10, 20, 30], &mut oc);

        assert_eq!(st.resolved.matched_positions, 3);
        assert_eq!(st.resolved.matched_positions, oracle.matched_positions);
        assert_eq!(st.resolved.row_ids, oracle.row_ids);
        assert_eq!(st.batch.num_rows(), 3);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn streaming_reread_resolves_out_of_window_backstep() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        write_indexed_dataset(&path, &[], 40_000, 1, 1).await;
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
        let reader = open_page_data_reader(&dataset).await.unwrap();
        let mut cursor = StreamingPositionCursor::new(reader, directory);

        let _ = cursor.resolve(&[39_000]).await.unwrap(); // seek far forward
        let back = cursor.resolve(&[100]).await.unwrap(); // way behind the window
        assert_eq!(
            back.matched_positions, 1,
            "out-of-window backstep must still match via page_lookup re-read"
        );
        assert_eq!(back.row_ids.len(), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn streaming_seek_skips_pages_before_band() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        write_indexed_dataset(&path, &[], 40_000, 1, 1).await;
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
        let reader = open_page_data_reader(&dataset).await.unwrap();
        let mut cursor = StreamingPositionCursor::new(reader, directory.clone());

        let r = cursor.resolve(&[39_000]).await.unwrap();
        assert_eq!(r.matched_positions, 1, "must still match after seeking");
        let expected_first = directory.first_page_for(39_000);
        assert!(
            expected_first > 0,
            "a late position should not start at page 0"
        );
        assert_eq!(
            cursor.first_loaded_page(),
            Some(expected_first),
            "cursor must seek to its band, not read from page 0"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn streaming_window_bounds_pages_and_handles_small_backstep() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        // 40k rows at positions 1..=40_000.
        write_indexed_dataset(&path, &[], 40_000, 1, 1).await;
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let materialized = load_start_index_from_lance_btree(&dataset).await.unwrap();
        let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
        let reader = open_page_data_reader(&dataset).await.unwrap();
        let mut cursor = StreamingPositionCursor::new(reader, directory);

        // Batch 1 deep in the stream; batch 2 dips back one position (extended probe).
        let r1 = cursor.resolve(&[20_000, 20_001]).await.unwrap();
        let r2 = cursor.resolve(&[19_999, 20_002]).await.unwrap();
        assert_eq!(r1.matched_positions, 2);
        assert_eq!(
            r2.matched_positions, 2,
            "small backstep within window must still match"
        );
        assert!(
            cursor.window_len() <= WINDOW_PAGES,
            "window must never exceed {WINDOW_PAGES} pages, got {}",
            cursor.window_len()
        );

        // Oracle parity for batch 1.
        let mut oc = 0usize;
        let o1 = materialized.resolve_sorted_positions_from_cursor(&[20_000, 20_001], &mut oc);
        assert_eq!(r1.row_ids, o1.row_ids);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn streaming_resolve_matches_materialized_with_multiallele_order() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        // Position 20 is multi-allelic (3 rows).
        let batch = test_variation_batch(vec![
            ("chr1", 10, 10, "A/T", 0, "rs10"),
            ("chr1", 20, 20, "C/G", 0, "rs20c"),
            ("chr1", 20, 20, "C/T", 0, "rs20t"),
            ("chr1", 20, 20, "C/A", 0, "rs20a"),
            ("chr1", 30, 30, "G/A", 0, "rs30"),
        ]);
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();

        // Oracle: materialized index.
        let materialized = load_start_index_from_lance_btree(&dataset).await.unwrap();
        let mut oracle_cursor = 0usize;
        let oracle =
            materialized.resolve_sorted_positions_from_cursor(&[10, 20, 30], &mut oracle_cursor);

        // Subject: streaming cursor.
        let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
        let reader = open_page_data_reader(&dataset).await.unwrap();
        let mut cursor = StreamingPositionCursor::new(reader, directory);
        let streamed = cursor.resolve(&[10, 20, 30]).await.unwrap();

        assert_eq!(streamed.requested_positions, oracle.requested_positions);
        assert_eq!(streamed.matched_positions, oracle.matched_positions);
        assert_eq!(
            streamed.row_ids, oracle.row_ids,
            "row_ids must match incl. per-position ascending order"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn page_directory_maps_position_to_first_page() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        // 10k rows at positions 1, 11, 21, ... 99_991.
        write_indexed_dataset(&path, &[], 10_000, 1, 10).await;
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();

        let directory = PositionPageDirectory::load(&dataset).await.unwrap();
        assert!(
            directory.num_pages() >= 2,
            "expected multiple pages for 10k rows"
        );
        assert_eq!(
            directory.batch_size(),
            4096,
            "batch_size read from page_lookup metadata"
        );

        let p_low = directory.first_page_for(1);
        let p_mid = directory.first_page_for(50_000);
        let p_high = directory.first_page_for(99_991);
        assert_eq!(p_low, 0, "first position lands on page 0");
        assert!(p_low <= p_mid, "first_page_for is monotonic non-decreasing");
        assert!(
            p_mid <= p_high,
            "first_page_for is monotonic non-decreasing"
        );
        assert!(p_high > 0, "a late position must seek past page 0");
    }

    #[tokio::test]
    async fn lookup_resolves_start_rows_with_take_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        let batch = test_variation_batch(vec![
            ("chr1", 10, 10, "A/T", 0, "rs10a"),
            ("chr1", 10, 10, "A/G", 0, "rs10b"),
            ("chr1", 20, 20, "C/G", 0, "rs20"),
        ]);
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();
        let lookup = SinglePathLanceVariationLookup::open(
            &path,
            vec![
                "start".into(),
                "end".into(),
                "allele_string".into(),
                "failed".into(),
                "variation_name".into(),
            ],
        )
        .await
        .unwrap();

        let mut cursor = lookup.new_cursor();
        let result = lookup
            .resolve_and_take(&[10, 20], &mut cursor)
            .await
            .unwrap();

        assert_eq!(result.batch.num_rows(), 3);
        assert_eq!(result.resolved.matched_positions, 2);
        assert_eq!(lookup.row_ids_len(), 3);
        // Streaming reports total rows as an upper bound for unique positions
        // (it never materializes the distinct-position count).
        assert_eq!(lookup.unique_positions(), lookup.row_ids_len());
    }

    #[test]
    fn runtime_projection_includes_required_lookup_columns_once() {
        let projection = ensure_runtime_projection(vec![
            "variation_name".into(),
            "start".into(),
            "allele_string".into(),
        ]);

        assert_eq!(
            projection,
            vec!["variation_name", "start", "allele_string", "end", "failed"]
        );
    }

    #[test]
    fn runtime_projection_drops_legacy_warm_cold_columns() {
        let projection = ensure_runtime_projection(vec![
            "position_key".into(),
            "variant_keys".into(),
            "tier".into(),
            "variation_name".into(),
            "position_key".into(),
        ]);

        assert_eq!(
            projection,
            vec!["variation_name", "start", "end", "allele_string", "failed"]
        );
    }

    #[tokio::test]
    async fn lance_lookup_round_trips_bundled_af_columns() {
        use crate::lance_cache::af_bundle::{af_column_order, bundle_af_columns};
        use datafusion::arrow::array::ArrayRef;

        let order = af_column_order();
        // 3 rows: all-present, gnomADg-only, all-absent.
        let af_for = |col: &str, row: usize| -> Option<String> {
            let is_gnomadg = col.starts_with("gnomADg");
            match row {
                0 => Some(format!("v_{col}")),
                1 => {
                    if is_gnomadg {
                        Some(format!("g_{col}"))
                    } else {
                        Some(String::new())
                    }
                }
                _ => Some(String::new()),
            }
        };
        let mut fields = vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ];
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(vec!["chr1", "chr1", "chr1"])),
            Arc::new(UInt32Array::from(vec![10u32, 20, 30])),
            Arc::new(UInt32Array::from(vec![10u32, 20, 30])),
            Arc::new(StringArray::from(vec!["A/T", "C/G", "G/A"])),
            Arc::new(Int8Array::from(vec![Some(0i8), Some(0), Some(0)])),
            Arc::new(StringArray::from(vec!["rs10", "rs20", "rs30"])),
        ];
        for col in &order {
            fields.push(Field::new(*col, DataType::Utf8, true));
            let vals: Vec<Option<String>> = (0..3).map(|r| af_for(col, r)).collect();
            cols.push(Arc::new(StringArray::from(vals)) as ArrayRef);
        }
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), cols).unwrap();
        let bundled = bundle_af_columns(&batch).unwrap();

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("var.lance");
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![bundled],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();

        // Projection includes the 27 AF column names (as the cold-parquet projection would).
        let mut projection: Vec<String> = vec![
            "variation_name".into(),
            "start".into(),
            "allele_string".into(),
        ];
        projection.extend(order.iter().map(|s| s.to_string()));
        let lookup = SinglePathLanceVariationLookup::open(&path, projection)
            .await
            .unwrap();
        let mut cursor = lookup.new_cursor();
        let result = lookup
            .resolve_and_take(&[10, 20, 30], &mut cursor)
            .await
            .unwrap();
        let out = result.batch;
        assert_eq!(out.num_rows(), 3);

        // Every one of the 27 AF columns is present (unbundled) with exact values.
        for col in &order {
            let idx = out
                .schema()
                .index_of(col)
                .expect("AF column present after unbundle");
            let arr = out
                .column(idx)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            assert_eq!(arr.value(0), format!("v_{col}"), "row0 {col}");
            let expected1 = if col.starts_with("gnomADg") {
                format!("g_{col}")
            } else {
                String::new()
            };
            assert_eq!(arr.value(1), expected1, "row1 {col}");
            assert_eq!(arr.value(2), "", "row2 {col}");
        }
    }

    fn test_variation_batch(rows: Vec<(&str, u32, u32, &str, i8, &str)>) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(
                    rows.iter().map(|row| row.0).collect::<Vec<_>>(),
                )),
                Arc::new(UInt32Array::from(
                    rows.iter().map(|row| row.1).collect::<Vec<_>>(),
                )),
                Arc::new(UInt32Array::from(
                    rows.iter().map(|row| row.2).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|row| row.3).collect::<Vec<_>>(),
                )),
                Arc::new(Int8Array::from(
                    rows.iter().map(|row| Some(row.4)).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|row| row.5).collect::<Vec<_>>(),
                )),
            ],
        )
        .unwrap()
    }
}
