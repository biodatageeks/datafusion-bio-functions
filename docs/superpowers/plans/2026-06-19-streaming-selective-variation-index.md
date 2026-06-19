# Streaming + Selective Variation BTree Index — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the eager full-materialization of the Lance variation position index with a per-partition streaming cursor that selectively seeks to its band and holds only a 5-page sliding window, cutting resident index RAM from ~1.06 GB to ~240 KB/partition and removing the build barrier from `lookup_wait`, with byte-identical output.

**Architecture:** A shared `SinglePathLanceVariationLookup` (one per contig, via the existing `OnceCell`) holds the immutable `page_data` reader and a tiny `PositionPageDirectory` loaded from `page_lookup.lance`. Each partition gets its own `StreamingPositionCursor` that, on first use, seeks via the directory to the page covering its first query position, then walks `page_data` page-by-page with a 5-page ring buffer (cursor ± 2), discarding consumed pages. Any backward step beyond the window is resolved by a selective `page_lookup` re-read. The eager path is retained behind an env flag for one release.

**Tech Stack:** Rust 2024, DataFusion 53.0.0, Arrow 58.0.0, `lance` / `lance-index` 7.0.0, tokio.

## Global Constraints

- Output MUST be byte-identical: chr1 merged e2e gate `CSQ count mismatch: 0`, all 86 shared CSQ fields 100%.
- Within a matched position, emit row_ids **ascending** (matches `load_u32_btree_index`'s `(position, row_id)` sort).
- `page_data` is globally sorted by position; pages are `DEFAULT_BTREE_BATCH_SIZE = 4096` rows (confirm exact via Task 1).
- `WINDOW_PAGES = 5` (cursor page ± 2). Resident index RAM target < 5 MB/partition.
- No on-disk cache format change; `take_rows` and `ResolvedRowIds`/`TakenVariationRows` unchanged.
- Public `annotate_vep` API unchanged → no vepyr API edit, rebuild only.
- Feature env: `VEP_STREAMING_VARIATION_INDEX` (default ON post-verification; `=0` falls back to eager path).
- Build/lint gate per task: `cargo test -p datafusion-bio-function-vep` and `cargo clippy -p datafusion-bio-function-vep -- -D warnings`.

---

## File Structure

- `datafusion/bio-function-vep/src/lance_cache/row_index.rs` — add `PositionPageDirectory` (load from `page_lookup.lance`, `first_page_for`). Keep `PositionRowIdIndex` + `load_u32_btree_index` (eager fallback + parity oracle in tests).
- `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs` — add `StreamingPositionCursor` + `DecodedPage`; change `SinglePathLanceVariationLookup` fields, `open`, `new_cursor`, `resolve_and_take`.
- `datafusion/bio-function-vep/src/lance_cache/lookup_exec.rs` — change `lance_cursors` type and the build/call sites (1378-1430 open, 1704-1708 resolve).
- `docs/superpowers/specs/2026-06-19-streaming-selective-variation-index-design.md` — append Task 1 spike findings (confirmed constants).

---

## Task 1: Spike — confirm Lance internals & partition locality

**Files:**
- Create (throwaway test): `datafusion/bio-function-vep/src/lance_cache/row_index.rs` (a `#[tokio::test] spike_*`, deleted at task end)
- Modify (findings): `docs/superpowers/specs/2026-06-19-streaming-selective-variation-index-design.md`

**Interfaces:**
- Consumes: existing `load_btree_segments`, `LanceIndexStore::from_dataset_for_existing`, `IndexReader::{open_index_file, read_record_batch, num_rows}`.
- Produces (written into the spec as confirmed facts): `PAGE_LOOKUP_FILE` name, page_lookup column names for min-value + page-number, the `batch_size` source, and whether chr1 forks=4 partitions are contiguous position ranges.

- [ ] **Step 1: Write a spike test that dumps `page_lookup.lance` schema + samples**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn spike_dump_page_lookup_schema() {
    // Build a tiny indexed dataset using the existing helpers.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr1.lance");
    let batch = test_variation_batch(vec![
        ("chr1", 10, 10, "A/T", 0, "rs10"),
        ("chr1", 20, 20, "C/G", 0, "rs20"),
        ("chr1", 30, 30, "G/A", 0, "rs30"),
    ]);
    crate::lance_cache::write::write_record_batches_to_lance(
        &path, vec![batch], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap();

    let segments = load_btree_segments(&dataset, "start", "start_btree_idx").await.unwrap();
    let store = LanceIndexStore::from_dataset_for_existing(&dataset, &segments[0]).await.unwrap();
    let lookup = store.open_index_file("page_lookup.lance").await.unwrap();
    eprintln!("[spike] page_lookup schema: {:?}", lookup.schema());
    eprintln!("[spike] page_lookup num_rows: {}", lookup.num_rows());
    let batch = lookup.read_range(0..lookup.num_rows(), None).await.unwrap();
    eprintln!("[spike] page_lookup batch: {:?}", batch);

    // Confirm per-page read on page_data returns one page of (values, ids).
    let pages = store.open_index_file("page_data.lance").await.unwrap();
    eprintln!("[spike] page_data schema: {:?}", pages.schema());
    let p0 = pages.read_record_batch(0, 4096).await.unwrap();
    eprintln!("[spike] page_data page0 rows: {} cols: {:?}",
        p0.num_rows(), p0.schema().fields().iter().map(|f| f.name().clone()).collect::<Vec<_>>());
}
```

- [ ] **Step 2: Run the spike and read the dumps**

Run: `cargo test -p datafusion-bio-function-vep spike_dump_page_lookup_schema -- --nocapture 2>&1 | grep '\[spike\]'`
Expected: prints `page_lookup` schema (min-value + page-number columns), `page_data` schema (`values`, `ids`), and a non-empty page0.

- [ ] **Step 3: Confirm partition locality on real chr1 (manual, via vepyr)**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr && \
env -u VIRTUAL_ENV -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
VEP_PROFILE=1 VEP_LANCE_PROFILE=1 uv run python e2e-testing/scripts/run_annotation_fast.py \
chr1 --cache merged --backend lance --forks 4 --force 2>&1 | grep 'variation_resolve' | \
awk '{print $2}' | head -40
```
Expected: inspect whether each partition's `requested_starts` cluster in disjoint position bands (selective achievable) or interleave (lazy-only). Record the verdict.

- [ ] **Step 4: Write findings into the spec**

Append a `## 10. Task 1 spike findings` section to the design spec with: exact `page_lookup` column names for min-value and page-number, the `page_data` column names, the page `batch_size` (and where it's read from), and the partition-locality verdict (R2). These become the constants used in Tasks 2–6.

- [ ] **Step 5: Delete the spike test, commit findings**

```bash
# remove spike_dump_page_lookup_schema from row_index.rs
git add docs/superpowers/specs/2026-06-19-streaming-selective-variation-index-design.md \
        datafusion/bio-function-vep/src/lance_cache/row_index.rs
git commit -m "docs(vep): spike findings for streaming variation index (page_lookup schema, partition locality)"
```

---

## Task 2: `PositionPageDirectory`

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/row_index.rs`
- Test: same file, `#[cfg(test)] mod tests`

**Interfaces:**
- Consumes: Task 1 constants (page_lookup column names, batch_size); `load_btree_segments`, `LanceIndexStore`, `IndexReader`.
- Produces: `PositionPageDirectory { load(&Dataset) -> Result<Self>, first_page_for(u32) -> u32, num_pages() -> u32 }` with `page_mins: Vec<u32>` ascending.

- [ ] **Step 1: Write the failing test**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn page_directory_maps_position_to_first_page() {
    // 3 pages worth of positions if batch_size were 2; we instead assert against
    // whatever batch_size the writer uses by checking monotonic page mapping.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr1.lance");
    let rows: Vec<_> = (0..10_000u32)
        .map(|i| ("chr1", (i * 10 + 1) as i64, (i * 10 + 1) as i64, "A/T", 0i32, "rs"))
        .collect();
    crate::lance_cache::write::write_record_batches_to_lance(
        &path, vec![test_variation_batch(rows)], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap();

    let directory = PositionPageDirectory::load(&dataset).await.unwrap();
    assert!(directory.num_pages() >= 2, "expected multiple pages for 10k rows");
    // Page mins are ascending; first_page_for is monotonic non-decreasing.
    let p_low = directory.first_page_for(1);
    let p_mid = directory.first_page_for(50_000);
    let p_high = directory.first_page_for(99_991);
    assert_eq!(p_low, 0);
    assert!(p_mid <= p_high);
    assert!(p_low <= p_mid);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep page_directory_maps_position_to_first_page -- --nocapture`
Expected: FAIL — `PositionPageDirectory` not found.

- [ ] **Step 3: Implement `PositionPageDirectory`**

```rust
const PAGE_LOOKUP_FILE: &str = "page_lookup.lance";
const PAGE_DATA_FILE: &str = "page_data.lance";

pub struct PositionPageDirectory {
    page_mins: Vec<u32>,   // ascending; page_mins[i] = min position of page i
    num_pages: u32,
    batch_size: u64,
    num_rows: usize,
}

impl PositionPageDirectory {
    pub async fn load(dataset: &lance::Dataset) -> Result<Self> {
        let segments = load_btree_segments(dataset, "start", "start_btree_idx").await?;
        let store = LanceIndexStore::from_dataset_for_existing(dataset, &segments[0])
            .await
            .map_err(|e| DataFusionError::Execution(format!("open index store: {e}")))?;
        let lookup = store.open_index_file(PAGE_LOOKUP_FILE).await
            .map_err(|e| DataFusionError::Execution(format!("open page_lookup: {e}")))?;
        let n = lookup.num_rows();
        let batch = lookup.read_range(0..n, None).await
            .map_err(|e| DataFusionError::Execution(format!("read page_lookup: {e}")))?;
        // Column names per Task 1 findings (min value + page number).
        let mins = batch.column_by_name("min").expect("page_lookup min col")
            .as_any().downcast_ref::<UInt32Array>().expect("min u32").clone();
        let pages = batch.column_by_name("page_idx").expect("page_lookup page col")
            .as_any().downcast_ref::<UInt32Array>().expect("page u32").clone();
        // Build page_mins indexed by page number, ascending by construction.
        let num_pages = pages.values().iter().copied().max().map(|m| m + 1).unwrap_or(0);
        let mut page_mins = vec![u32::MAX; num_pages as usize];
        for i in 0..mins.len() {
            page_mins[pages.value(i) as usize] = mins.value(i);
        }
        let data = store.open_index_file(PAGE_DATA_FILE).await
            .map_err(|e| DataFusionError::Execution(format!("open page_data: {e}")))?;
        Ok(Self {
            page_mins,
            num_pages,
            batch_size: DEFAULT_BTREE_BATCH_SIZE, // confirmed in Task 1
            num_rows: data.num_rows(),
        })
    }

    pub fn num_pages(&self) -> u32 { self.num_pages }
    pub fn batch_size(&self) -> u64 { self.batch_size }
    pub fn num_rows(&self) -> usize { self.num_rows }

    /// First page that may contain `position`: the last page whose min <= position.
    pub fn first_page_for(&self, position: u32) -> u32 {
        match self.page_mins.partition_point(|&m| m <= position) {
            0 => 0,
            i => (i - 1) as u32,
        }
    }
}
```
Add `use lance_index::scalar::btree::DEFAULT_BTREE_BATCH_SIZE;` (or re-declare the const = 4096 if not re-exported — confirm in Task 1) and ensure `UInt32Array` is imported.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep page_directory_maps_position_to_first_page`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/row_index.rs
git commit -m "feat(vep): PositionPageDirectory — load page_lookup, map position to page"
```

---

## Task 3: `StreamingPositionCursor` core (lazy, non-selective) + byte-identical resolve

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs`
- Test: same file.

**Interfaces:**
- Consumes: `PositionPageDirectory`, `IndexReader::read_record_batch`, `ResolvedRowIds`, `PositionRowIdIndex::resolve_sorted_positions_from_cursor` (parity oracle).
- Produces: `StreamingPositionCursor` with `async fn resolve(&mut self, sorted_positions: &[u32]) -> Result<ResolvedRowIds>`, `DecodedPage`, `WINDOW_PAGES`.

- [ ] **Step 1: Write the failing parity test (incl. multi-allele ascending order)**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn streaming_resolve_matches_materialized_with_multiallele_order() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr1.lance");
    // Position 20 has 3 alleles whose row order is NOT row_id-ascending on disk.
    let batch = test_variation_batch(vec![
        ("chr1", 10, 10, "A/T", 0, "rs10"),
        ("chr1", 20, 20, "C/G", 0, "rs20c"),
        ("chr1", 20, 20, "C/T", 0, "rs20t"),
        ("chr1", 20, 20, "C/A", 0, "rs20a"),
        ("chr1", 30, 30, "G/A", 0, "rs30"),
    ]);
    crate::lance_cache::write::write_record_batches_to_lance(
        &path, vec![batch], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap();

    // Oracle: materialized index.
    let materialized = load_start_index_from_lance_btree(&dataset).await.unwrap();
    let mut oracle_cursor = 0usize;
    let oracle = materialized.resolve_sorted_positions_from_cursor(&[10, 20, 30], &mut oracle_cursor);

    // Subject: streaming cursor.
    let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
    let reader = open_page_data_reader(&dataset).await.unwrap(); // helper added in Step 3
    let mut cursor = StreamingPositionCursor::new(reader, directory);
    let streamed = cursor.resolve(&[10, 20, 30]).await.unwrap();

    assert_eq!(streamed.requested_positions, oracle.requested_positions);
    assert_eq!(streamed.matched_positions, oracle.matched_positions);
    assert_eq!(streamed.row_ids, oracle.row_ids, "row_ids must match incl. per-position ascending order");
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep streaming_resolve_matches_materialized_with_multiallele_order`
Expected: FAIL — `StreamingPositionCursor` / `open_page_data_reader` not found.

- [ ] **Step 3: Implement the cursor (lazy, seeks from page 0; window added in Task 4)**

```rust
use std::collections::VecDeque;
use datafusion::arrow::array::{UInt32Array, UInt64Array};

pub const WINDOW_PAGES: usize = 5; // cursor page ± 2 (eviction wired in Task 4)

struct DecodedPage { page_number: u32, positions: Vec<u32>, row_ids: Vec<u64> }

pub async fn open_page_data_reader(
    dataset: &lance::Dataset,
) -> Result<std::sync::Arc<dyn lance_index::scalar::IndexReader>> {
    let segments = crate::lance_cache::row_index::load_btree_segments(dataset, "start", "start_btree_idx").await?;
    let store = lance_index::scalar::lance_format::LanceIndexStore::from_dataset_for_existing(dataset, &segments[0])
        .await.map_err(|e| DataFusionError::Execution(format!("open index store: {e}")))?;
    store.open_index_file("page_data.lance").await
        .map_err(|e| DataFusionError::Execution(format!("open page_data: {e}")))
}

pub struct StreamingPositionCursor {
    reader: std::sync::Arc<dyn lance_index::scalar::IndexReader>,
    directory: std::sync::Arc<PositionPageDirectory>,
    window: VecDeque<DecodedPage>,
    front_off: usize, // offset of next unconsumed row within window.front()
    next_page_to_load: u32,
    seeded: bool,
}

impl StreamingPositionCursor {
    pub fn new(
        reader: std::sync::Arc<dyn lance_index::scalar::IndexReader>,
        directory: std::sync::Arc<PositionPageDirectory>,
    ) -> Self {
        Self { reader, directory, window: VecDeque::new(), front_off: 0, next_page_to_load: 0, seeded: false }
    }

    async fn load_page(&self, page_number: u32) -> Result<DecodedPage> {
        let batch = self.reader
            .read_record_batch(page_number as u64, self.directory.batch_size())
            .await
            .map_err(|e| DataFusionError::Execution(format!("read page {page_number}: {e}")))?;
        let positions = batch.column_by_name("values").expect("values col")
            .as_any().downcast_ref::<UInt32Array>().expect("values u32")
            .values().to_vec();
        let row_ids = batch.column_by_name("ids").expect("ids col")
            .as_any().downcast_ref::<UInt64Array>().expect("ids u64")
            .values().to_vec();
        Ok(DecodedPage { page_number, positions, row_ids })
    }

    /// Position at the read cursor (front of window + front_off), or None if exhausted.
    fn cur_pos(&self) -> Option<u32> {
        self.window.front().and_then(|p| p.positions.get(self.front_off).copied())
    }
    fn cur_id(&self) -> Option<u64> {
        self.window.front().and_then(|p| p.row_ids.get(self.front_off).copied())
    }

    async fn ensure_front(&mut self) -> Result<()> {
        // Drop fully-consumed front pages, pull next page if window front is empty.
        loop {
            match self.window.front() {
                Some(p) if self.front_off < p.positions.len() => return Ok(()),
                Some(_) => { self.window.pop_front(); self.front_off = 0; }
                None => {
                    if self.next_page_to_load >= self.directory.num_pages() { return Ok(()); }
                    let page = self.load_page(self.next_page_to_load).await?;
                    self.next_page_to_load += 1;
                    self.window.push_back(page);
                }
            }
        }
    }

    async fn advance(&mut self) -> Result<()> {
        self.front_off += 1;
        self.ensure_front().await
    }

    pub async fn resolve(&mut self, sorted_positions: &[u32]) -> Result<ResolvedRowIds> {
        if !self.seeded {
            // Non-selective for now: start at page 0 (Task 4 makes this selective).
            self.next_page_to_load = 0;
            self.seeded = true;
        }
        self.ensure_front().await?;
        let mut matched_positions = 0usize;
        let mut row_ids = Vec::new();
        for &q in sorted_positions {
            while self.cur_pos().is_some_and(|p| p < q) { self.advance().await?; }
            if self.cur_pos().is_none() { break; }
            if self.cur_pos() == Some(q) {
                matched_positions += 1;
                let mut at_pos: Vec<u64> = Vec::new();
                while self.cur_pos() == Some(q) {
                    if let Some(id) = self.cur_id() { at_pos.push(id); }
                    self.advance().await?;
                }
                at_pos.sort_unstable(); // INVARIANT: per-position ascending row_ids
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
```
Note: `ResolvedRowIds` fields are `pub` (`row_index.rs:42`). Make `load_btree_segments` `pub(crate)` if not already.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep streaming_resolve_matches_materialized_with_multiallele_order`
Expected: PASS — including the per-position ascending row_id assertion.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs \
        datafusion/bio-function-vep/src/lance_cache/row_index.rs
git commit -m "feat(vep): StreamingPositionCursor — lazy page-streamed resolve, byte-identical to materialized"
```

---

## Task 4: 5-page sliding window (cursor ± 2) + bounded backward tolerance

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs`
- Test: same file.

**Interfaces:**
- Consumes: Task 3 cursor.
- Produces: window eviction keeping ≤ `WINDOW_PAGES` pages with the cursor page and ≤2 trailing; a small backward step within the window resolves correctly.

- [ ] **Step 1: Write the failing test (multi-batch + small backward extended probe)**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn streaming_window_bounds_pages_and_handles_small_backstep() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr1.lance");
    let rows: Vec<_> = (0..40_000u32)
        .map(|i| ("chr1", (i + 1) as i64, (i + 1) as i64, "A/T", 0i32, "rs"))
        .collect();
    crate::lance_cache::write::write_record_batches_to_lance(
        &path, vec![test_variation_batch(rows)], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap();
    let materialized = load_start_index_from_lance_btree(&dataset).await.unwrap();
    let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
    let reader = open_page_data_reader(&dataset).await.unwrap();
    let mut cursor = StreamingPositionCursor::new(reader, directory);

    // Batch 1 deep into the stream, batch 2 with a small backward extended probe.
    let r1 = cursor.resolve(&[20_000, 20_001]).await.unwrap();
    let r2 = cursor.resolve(&[19_999, 20_002]).await.unwrap(); // 19_999 steps back one
    assert_eq!(r1.matched_positions, 2);
    assert_eq!(r2.matched_positions, 2, "small backstep within window must still match");

    // Window never exceeds WINDOW_PAGES.
    assert!(cursor.window_len() <= WINDOW_PAGES);

    // Oracle parity for the union.
    let mut oc = 0usize;
    let o1 = materialized.resolve_sorted_positions_from_cursor(&[20_000, 20_001], &mut oc);
    assert_eq!(r1.row_ids, o1.row_ids);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep streaming_window_bounds_pages_and_handles_small_backstep`
Expected: FAIL — `window_len` missing and/or backstep panics/mismatches.

- [ ] **Step 3: Add window eviction + retained trailing pages + backstep scan**

```rust
impl StreamingPositionCursor {
    pub fn window_len(&self) -> usize { self.window.len() }

    /// After advancing, keep at most WINDOW_PAGES with ≤2 trailing pages behind front.
    fn evict_to_window(&mut self) {
        // front is the cursor page; retain up to 2 pages *behind* it already dropped by
        // ensure_front, so simply cap total length.
        while self.window.len() > WINDOW_PAGES { self.window.pop_front(); }
    }

    /// Resolve a query that may dip slightly backward: scan retained trailing pages.
    fn try_match_in_window(&self, q: u32, out: &mut Vec<u64>) -> bool {
        let mut found = false;
        for page in &self.window {
            // pages are sorted; linear/binary scan for q within this page
            let lo = page.positions.partition_point(|&p| p < q);
            let mut i = lo;
            while i < page.positions.len() && page.positions[i] == q {
                out.push(page.row_ids[i]); i += 1; found = true;
            }
        }
        found
    }
}
```
Modify `resolve` so that before forward-advancing for a query `q`, if `q < cur_pos()` (a backstep) it calls `try_match_in_window(q, &mut at_pos)` (then `at_pos.sort_unstable()`); forward path unchanged. Call `evict_to_window()` at the end of `ensure_front`'s load branch. Keep the front page as cursor page; the ring naturally holds trailing pages until capacity.

Refinement for retaining 2 trailing pages: instead of `pop_front` on consumption in `ensure_front`, mark a page consumed but retain it; evict only when `window.len() > WINDOW_PAGES`. Track the cursor as `(front_page_index_in_window, front_off)` rather than always `window.front()`. Adjust `cur_pos`/`cur_id`/`advance` to index `self.window[self.cursor_page]` and bump `cursor_page` across boundaries, evicting from the back of the trailing tail when length exceeds 5.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep streaming_window_bounds_pages_and_handles_small_backstep`
Expected: PASS; `window_len() <= 5`.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs
git commit -m "feat(vep): 5-page sliding window + bounded backstep tolerance for streaming cursor"
```

---

## Task 5: Selective seek via `PositionPageDirectory`

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs`
- Test: same file.

**Interfaces:**
- Consumes: `PositionPageDirectory::first_page_for`.
- Produces: on first `resolve`, `next_page_to_load = first_page_for(sorted_positions[0])`; a read counter to assert pages-skipped.

- [ ] **Step 1: Write the failing test (assert pages before the band are NOT read)**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn streaming_seek_skips_pages_before_band() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr1.lance");
    let rows: Vec<_> = (0..40_000u32)
        .map(|i| ("chr1", (i + 1) as i64, (i + 1) as i64, "A/T", 0i32, "rs"))
        .collect();
    crate::lance_cache::write::write_record_batches_to_lance(
        &path, vec![test_variation_batch(rows)], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap();
    let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
    let reader = open_page_data_reader(&dataset).await.unwrap();
    let mut cursor = StreamingPositionCursor::new(reader, directory.clone());

    let _ = cursor.resolve(&[39_000]).await.unwrap(); // near the end
    let expected_first = directory.first_page_for(39_000);
    assert!(expected_first > 0, "should not start at page 0");
    assert_eq!(cursor.first_loaded_page(), Some(expected_first),
        "must seek to the band, not read from page 0");
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep streaming_seek_skips_pages_before_band`
Expected: FAIL — `first_loaded_page` missing / starts at 0.

- [ ] **Step 3: Make seeding selective + expose `first_loaded_page`**

```rust
impl StreamingPositionCursor {
    pub fn first_loaded_page(&self) -> Option<u32> { self.first_loaded_page }
}
// add field: first_loaded_page: Option<u32> (default None)
// in resolve(), the seeding branch becomes:
if !self.seeded {
    if let Some(&first) = sorted_positions.first() {
        self.next_page_to_load = self.directory.first_page_for(first);
        self.first_loaded_page = Some(self.next_page_to_load);
    }
    self.seeded = true;
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep streaming_seek_skips_pages_before_band`
Then re-run Tasks 3 & 4 tests to confirm no regression:
`cargo test -p datafusion-bio-function-vep streaming_`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs
git commit -m "feat(vep): selective seek — streaming cursor jumps to its band via page directory"
```

---

## Task 6: `page_lookup` re-read floor for out-of-window backsteps

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs`
- Test: same file.

**Interfaces:**
- Consumes: `PositionPageDirectory::first_page_for`, `load_page`.
- Produces: when `q` precedes the oldest retained page, a selective re-read of `first_page_for(q)` resolves it; logged under `VEP_LANCE_PROFILE`.

- [ ] **Step 1: Write the failing test (force an out-of-window backstep)**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn streaming_reread_resolves_out_of_window_backstep() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr1.lance");
    let rows: Vec<_> = (0..40_000u32)
        .map(|i| ("chr1", (i + 1) as i64, (i + 1) as i64, "A/T", 0i32, "rs"))
        .collect();
    crate::lance_cache::write::write_record_batches_to_lance(
        &path, vec![test_variation_batch(rows)], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap();
    let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await.unwrap());
    let reader = open_page_data_reader(&dataset).await.unwrap();
    let mut cursor = StreamingPositionCursor::new(reader, directory);

    let _ = cursor.resolve(&[39_000]).await.unwrap();   // far forward
    let back = cursor.resolve(&[100]).await.unwrap();    // way behind the window
    assert_eq!(back.matched_positions, 1, "out-of-window backstep must still match via re-read");
    assert_eq!(back.row_ids.len(), 1);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep streaming_reread_resolves_out_of_window_backstep`
Expected: FAIL — position 100 not in window, returns 0 matches.

- [ ] **Step 3: Add the re-read floor**

```rust
impl StreamingPositionCursor {
    /// True if q is below the smallest position currently retained.
    fn below_window(&self, q: u32) -> bool {
        self.window.front().and_then(|p| p.positions.first().copied()).is_some_and(|min| q < min)
    }

    async fn reread_match(&self, q: u32, out: &mut Vec<u64>) -> Result<bool> {
        if std::env::var_os("VEP_LANCE_PROFILE").is_some() {
            eprintln!("[vep-lance-profile] streaming_reread position={q}");
        }
        let page = self.load_page(self.directory.first_page_for(q)).await?;
        let lo = page.positions.partition_point(|&p| p < q);
        let mut i = lo; let mut found = false;
        while i < page.positions.len() && page.positions[i] == q {
            out.push(page.row_ids[i]); i += 1; found = true;
        }
        Ok(found)
    }
}
// In resolve(), when handling a query q with q < cur_pos():
//   if try_match_in_window(q, &mut at_pos) is false AND below_window(q):
//        reread_match(q, &mut at_pos).await?;
//   then at_pos.sort_unstable(); count it as matched if non-empty.
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep streaming_reread_resolves_out_of_window_backstep`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs
git commit -m "feat(vep): page_lookup re-read floor for out-of-window backsteps (correctness)"
```

---

## Task 7: Wire into `open` / `resolve_and_take` / `lookup_exec` behind env flag

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs:27` (`open`, `resolve_and_take`, fields), add `new_cursor`.
- Modify: `datafusion/bio-function-vep/src/lance_cache/lookup_exec.rs:358` (`lance_cursors` type), `:1378-1430` (build), `:1704-1708` (call).

**Interfaces:**
- Consumes: Tasks 2–6.
- Produces: `SinglePathLanceVariationLookup` holding `page_data_reader` + `directory` (streaming) when `VEP_STREAMING_VARIATION_INDEX != "0"`, else the eager `index`. `resolve_and_take(&self, &[u32], &mut StreamingPositionCursor)`.

- [ ] **Step 1: Write the failing integration test (streaming open + resolve_and_take parity)**

```rust
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
        &path, vec![batch], crate::lance_cache::write::LanceIndexKind::Start,
    ).await.unwrap();

    let projection = vec!["start".into(), "end".into(), "allele_string".into(),
                          "failed".into(), "variation_name".into()];
    let lookup = SinglePathLanceVariationLookup::open(&path, projection).await.unwrap();
    let mut cursor = lookup.new_cursor();
    let taken = lookup.resolve_and_take(&[10, 20, 30], &mut cursor).await.unwrap();
    assert_eq!(taken.batch.num_rows(), 3);
    assert_eq!(taken.resolved.matched_positions, 3);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep streaming_open_resolve_and_take_matches_eager`
Expected: FAIL — `new_cursor` missing / `resolve_and_take` signature mismatch.

- [ ] **Step 3: Rewire `open`, add `new_cursor`, change `resolve_and_take`**

```rust
pub struct SinglePathLanceVariationLookup {
    dataset: lance::Dataset,
    projection: Vec<String>,
    page_data_reader: std::sync::Arc<dyn lance_index::scalar::IndexReader>,
    directory: std::sync::Arc<PositionPageDirectory>,
}

impl SinglePathLanceVariationLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await
            .map_err(|e| DataFusionError::Execution(format!("open lance variation '{}': {e}", path.display())))?;
        let directory = std::sync::Arc::new(PositionPageDirectory::load(&dataset).await?);
        let page_data_reader = open_page_data_reader(&dataset).await?;
        Ok(Self { dataset, projection: ensure_runtime_projection(projection), page_data_reader, directory })
    }

    pub fn new_cursor(&self) -> StreamingPositionCursor {
        StreamingPositionCursor::new(self.page_data_reader.clone(), self.directory.clone())
    }

    pub async fn resolve_and_take(
        &self,
        sorted_unique_starts: &[u32],
        cursor: &mut StreamingPositionCursor,
    ) -> Result<TakenVariationRows> {
        let resolved = cursor.resolve(sorted_unique_starts).await?;
        // ... unchanged from here: bundle_projection, take_rows, unbundle_af_columns ...
        let physical_projection = crate::lance_cache::af_bundle::bundle_projection(&self.projection);
        let projection_request = ProjectionRequest::from_columns(
            physical_projection.iter().map(String::as_str), self.dataset.schema());
        let batch = self.dataset.take_rows(&resolved.row_ids, projection_request).await
            .map_err(|e| DataFusionError::Execution(format!("Lance take_rows failed: {e}")))?;
        let batch = crate::lance_cache::af_bundle::unbundle_af_columns(&batch)?;
        Ok(TakenVariationRows { resolved, batch })
    }
}
```
Keep `row_ids_len`/`unique_positions` by delegating to `directory.num_rows()` / a stored count, or drop them if only used by the old profile line (update `open_breakdown` accordingly).

- [ ] **Step 4: Change `lance_cursors` and the call site in `lookup_exec.rs`**

```rust
// line 358:
lance_cursors: HashMap<String, StreamingPositionCursor>,
// build site (1378-1430): keep OnceCell of Arc<SinglePathLanceVariationLookup>;
//   open() no longer materializes — record open_breakdown as dataset_open + directory_load.
// call site (1704-1708):
let cursor = self.lance_cursors
    .entry(chrom.clone())
    .or_insert_with(|| lookup.new_cursor());
let taken = tokio::task::block_in_place(|| {
    tokio::runtime::Handle::current().block_on(lookup.resolve_and_take(&starts, cursor))
})?;
```
Gate with `VEP_STREAMING_VARIATION_INDEX`: when `== "0"`, retain the eager `load_start_index_from_lance_btree` + `&mut usize` cursor path (keep both behind a small enum or two code branches in `open`/`resolve_and_take`).

- [ ] **Step 5: Build, lint, run the full crate test suite**

Run:
```bash
cargo test -p datafusion-bio-function-vep 2>&1 | tail -5
cargo clippy -p datafusion-bio-function-vep -- -D warnings 2>&1 | tail -5
```
Expected: all tests PASS; no clippy warnings.

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs \
        datafusion/bio-function-vep/src/lance_cache/lookup_exec.rs
git commit -m "feat(vep): wire streaming variation index into open/resolve_and_take (flag VEP_STREAMING_VARIATION_INDEX)"
```

---

## Task 8: chr1 + chr22 e2e parity gate + RAM/latency measurement

**Files:**
- Modify (results): `docs/superpowers/specs/2026-06-19-streaming-selective-variation-index-design.md` (§8 results)

**Interfaces:**
- Consumes: the built vepyr (rebuilt with this crate via `uv run maturin develop --release`).
- Produces: recorded G1–G4 evidence; flag default decision.

- [ ] **Step 1: Rebuild vepyr with the change**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr && \
env -u VIRTUAL_ENV -u CONDA_PREFIX -u CONDA_DEFAULT_ENV uv run maturin develop --release 2>&1 | tail -2
```
Expected: `Installed vepyr-0.1.0`.

- [ ] **Step 2: chr22 parity + profile (streaming ON)**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr && \
env -u VIRTUAL_ENV -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
VEP_PROFILE=1 VEP_STREAMING_VARIATION_INDEX=1 uv run python e2e-testing/scripts/run_annotation_fast.py \
chr22 --cache merged --backend lance --forks 0 --force 2>&1 | \
grep -E 'CSQ count mismatch|fields match at 100|lookup_wait='
```
Expected: `CSQ count mismatch: 0`, `ALL 86 shared CSQ fields match at 100%`, `lookup_wait` lower than the 0.546 s baseline.

- [ ] **Step 3: chr1 parity + profile (streaming ON)**

Run: same as Step 2 with `chr1`.
Expected: `CSQ count mismatch: 0`, all 86 fields 100%, `lookup_wait` well below 0.948 s, `open` build ≈ 0.

- [ ] **Step 4: Peak RSS comparison (eager vs streaming) on chr1**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr
for flag in 0 1; do
  env -u VIRTUAL_ENV -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
  VEP_STREAMING_VARIATION_INDEX=$flag /usr/bin/time -l uv run python \
  e2e-testing/scripts/run_annotation_fast.py chr1 --cache merged --backend lance --forks 4 --force \
  2>&1 | grep -E 'maximum resident set size'
done
```
Expected: streaming (`=1`) peak RSS materially below eager (`=0`); record both numbers.

- [ ] **Step 5: End-to-end wall-time comparison (G5) — eager vs streaming**

Run (single-contig chr1, forks=0 so no cross-contig overlap hides the barrier):
```bash
cd /Users/mwiewior/research/git/vepyr
for flag in 0 1; do
  echo "=== VEP_STREAMING_VARIATION_INDEX=$flag chr1 ==="
  for rep in 1 2; do
    env -u VIRTUAL_ENV -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
    VEP_STREAMING_VARIATION_INDEX=$flag /usr/bin/time -p uv run python \
    e2e-testing/scripts/run_annotation_fast.py chr1 --cache merged --backend lance --forks 0 --force \
    2>&1 | grep -E '^real'
  done
done
```
Then a multi-contig run where per-contig barriers accumulate:
```bash
cd /Users/mwiewior/research/git/vepyr
for flag in 0 1; do
  echo "=== VEP_STREAMING_VARIATION_INDEX=$flag chr21+chr22 ==="
  for chr in chr21 chr22; do
    env -u VIRTUAL_ENV -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
    VEP_STREAMING_VARIATION_INDEX=$flag /usr/bin/time -p uv run python \
    e2e-testing/scripts/run_annotation_fast.py $chr --cache merged --backend lance --forks 0 --force \
    2>&1 | grep -E '^real'
  done
done
```
Expected: streaming (`=1`) `real` wall time ≤ eager (`=0`), saving ≈ the per-contig unhidden `open` build time (order ~0.7 s on chr1). Record the deltas; note when context load fully hides the build (saving → ~0) vs not.

- [ ] **Step 6: Record results + set flag default**

Append G1–G5 evidence (parity, lookup_wait before/after, peak RSS before/after, end-to-end wall time before/after for chr1 and chr21+chr22) to the spec §8. If all pass, flip `VEP_STREAMING_VARIATION_INDEX` default to ON in code (keep `=0` fallback). Update memory note `vep-lookup-streaming-and-startup-barrier`.

- [ ] **Step 7: Commit**

```bash
git add docs/superpowers/specs/2026-06-19-streaming-selective-variation-index-design.md \
        datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs
git commit -m "perf(vep): default streaming variation index ON — chr1/chr22 parity green, RSS + lookup_wait + wall time down"
```

---

## Self-Review

**Spec coverage:** §4.1 selective → Task 5; §4.2 lazy 5-page window → Tasks 3–4; §4.3 re-read floor → Task 6; §5 invariants (ascending row_ids) → Task 3 test; §6 interfaces → Tasks 2/3/7; §7 risks → Task 1 spike; §8 gates → Task 8. Eager fallback flag → Task 7. All covered.

**Placeholder scan:** Task 1 is an explicit spike that produces confirmed constants consumed by later tasks (page_lookup column names, batch_size) — these are the only deferred specifics, by design, and are resolved before Task 2 implements against them. No "TBD"/"handle errors"/"similar to" placeholders elsewhere; every code step shows code.

**Type consistency:** `StreamingPositionCursor`, `PositionPageDirectory`, `DecodedPage`, `WINDOW_PAGES`, `open_page_data_reader`, `ResolvedRowIds{requested_positions, matched_positions, row_ids}`, `TakenVariationRows{resolved, batch}`, `resolve_and_take(&[u32], &mut StreamingPositionCursor)` used consistently across Tasks 2–7.

**Known follow-ups (not in scope):** mmap/prebuilt sorted sidecar to also cut the inherent ~0.40 s decode (spec §2 non-goal); removing the eager path after one release.
