//! Footer page-index resolution + a coalescing async reader for Parquet
//! point lookups.
//!
//! Ported from the validated benchmark `the validated `pq_take` encoding benchmark`
//! (402,085 rows byte-identical to Parquet). Production adaptations: errors via
//! `DataFusionError` instead of `anyhow`; bound to the direct `parquet` 58 crate.
//!
//! - [`PageDir`] is the built-in equivalent of Parquet's chunked BTree leaf
//!   directory: a per-page `min`/`max` of the key column + each page's global
//!   row range, read straight from the footer `ColumnIndex`/`OffsetIndex` — no
//!   full-column scan. The key column is written sorted within each tier, so the
//!   pages form 1–2 monotonically-increasing runs; [`PageDir::resolve_ranges`]
//!   binary-searches within each run.
//! - [`CoalescingAsyncReader`] merges nearby page byte-ranges into fewer, larger
//!   reads (like Parquet / `ParquetObjectReader`), which is what makes cold I/O
//!   competitive (measured ~3.3× fewer read-ops than Parquet on chr1).

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use bytes::Bytes;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, RowSelection, RowSelector};

/// Cache-independent read-op / byte counters (ops and bytes are the same warm or
/// cold, so warm measurement reveals the cold I/O profile).
#[derive(Clone, Default)]
pub struct IoCounters {
    bytes: Arc<AtomicU64>,
    ops: Arc<AtomicU64>,
}

impl IoCounters {
    pub fn new() -> Self {
        Self::default()
    }

    /// `(bytes_read, read_ops)`.
    pub fn snapshot(&self) -> (u64, u64) {
        (
            self.bytes.load(Ordering::Relaxed),
            self.ops.load(Ordering::Relaxed),
        )
    }
}

/// A sparse position → page directory built from the footer page index.
///
/// Keys are held as `u64` regardless of the physical key column width so a single
/// implementation serves both the `INT32` `start` key (variation) and the `INT64`
/// `key` (translation_sift). `u32` keys are widened losslessly at build time and
/// callers cast their probes up to `u64`.
pub struct PageDir {
    min: Vec<u64>,
    max: Vec<u64>,
    first_row: Vec<u64>,
    end_row: Vec<u64>,
    /// Contiguous sorted runs `[start, end)` in page order, split where the key
    /// `min` resets downward (the warm→cold tier boundary).
    runs: Vec<(usize, usize)>,
}

impl PageDir {
    /// Build from preloaded reader metadata (must have been loaded with the page
    /// index enabled). `key_leaf` is the leaf column index of the sort key
    /// (`start` for variation, `key` for translation_sift), which must be an
    /// `INT32` or `INT64` physical column. `INT32`/`UInt32` values are widened to
    /// `u64` losslessly (the `as u32` step recovers the unsigned value first).
    pub fn build(meta: &ArrowReaderMetadata, key_leaf: usize) -> Result<Self> {
        use parquet::file::page_index::column_index::ColumnIndexMetaData;
        let md = meta.metadata();
        let ci = md
            .column_index()
            .ok_or_else(|| DataFusionError::Execution("parquet has no column index".into()))?;
        let oi = md
            .offset_index()
            .ok_or_else(|| DataFusionError::Execution("parquet has no offset index".into()))?;
        let (mut min, mut max, mut first_row, mut end_row) = (vec![], vec![], vec![], vec![]);
        let mut rg_start = 0u64;
        for rg in 0..md.num_row_groups() {
            let rg_rows = md.row_group(rg).num_rows() as u64;
            let col = &ci[rg][key_leaf];
            let (mins, maxs): (Vec<u64>, Vec<u64>) = match col {
                ColumnIndexMetaData::INT32(p) => (
                    p.min_values().iter().map(|&v| v as u32 as u64).collect(),
                    p.max_values().iter().map(|&v| v as u32 as u64).collect(),
                ),
                ColumnIndexMetaData::INT64(p) => (
                    p.min_values().iter().map(|&v| v as u64).collect(),
                    p.max_values().iter().map(|&v| v as u64).collect(),
                ),
                _ => {
                    return Err(DataFusionError::Execution(
                        "parquet key column index is not INT32 or INT64".into(),
                    ));
                }
            };
            let locs = oi[rg][key_leaf].page_locations();
            for p in 0..mins.len() {
                let fr = rg_start + locs[p].first_row_index as u64;
                let er = if p + 1 < locs.len() {
                    rg_start + locs[p + 1].first_row_index as u64
                } else {
                    rg_start + rg_rows
                };
                min.push(mins[p]);
                max.push(maxs[p]);
                first_row.push(fr);
                end_row.push(er);
            }
            rg_start += rg_rows;
        }
        Ok(Self::from_parts(min, max, first_row, end_row))
    }

    /// Construct from raw per-page arrays, detecting the sorted runs. Exposed for
    /// unit tests (and reused by [`build`]).
    fn from_parts(min: Vec<u64>, max: Vec<u64>, first_row: Vec<u64>, end_row: Vec<u64>) -> Self {
        let mut runs = vec![];
        let mut s = 0usize;
        for i in 1..min.len() {
            if min[i] < min[i - 1] {
                runs.push((s, i));
                s = i;
            }
        }
        if !min.is_empty() {
            runs.push((s, min.len()));
        }
        Self {
            min,
            max,
            first_row,
            end_row,
            runs,
        }
    }

    /// Merged candidate page row-ranges covering any of `positions`.
    pub fn resolve_ranges(&self, positions: &[u64]) -> Vec<(u64, u64)> {
        let mut pages: Vec<usize> = Vec::new();
        for &p in positions {
            for &(rs, re) in &self.runs {
                // First page in this run with max >= p (lower bound).
                let mut lo = rs;
                let mut hi = re;
                while lo < hi {
                    let m = (lo + hi) / 2;
                    if self.max[m] < p {
                        lo = m + 1;
                    } else {
                        hi = m;
                    }
                }
                // `p` may span several consecutive pages when a run of identical
                // key values straddles page boundaries (a highly multi-allelic
                // position). `max` is non-decreasing within a sorted run, so from
                // `lo` onward `max >= p` holds; include every page while its `min`
                // is still `<= p`. Pushing only `lo` would drop the rows of `p` that
                // live in the following page(s).
                let mut pg = lo;
                while pg < re && self.min[pg] <= p {
                    pages.push(pg);
                    pg += 1;
                }
                // A run's FIRST page can STRADDLE the tier seam: the shard streams
                // the warm tier then the cold tier without forcing a page break at
                // the boundary, so one page holds the warm tail (high keys)
                // followed by the cold head (low keys). Its key `min` resets
                // downward (which is exactly what starts this run) while its `max`
                // stays high (a warm value), so its `[min,max]` spans both tiers.
                // That out-of-order `max` breaks the monotonic-max binary search
                // above (which skips it), dropping the warm rows living in that
                // page. Bracket-test the run's first page directly to recover them.
                // Deduped below, so at most one extra page per batch.
                if self.min[rs] <= p && p <= self.max[rs] {
                    pages.push(rs);
                }
            }
        }
        pages.sort_unstable();
        pages.dedup();
        // Merge adjacent page row-ranges.
        let mut out: Vec<(u64, u64)> = Vec::new();
        for &pg in &pages {
            let (a, b) = (self.first_row[pg], self.end_row[pg]);
            if let Some(last) = out.last_mut()
                && last.1 == a
            {
                last.1 = b;
            } else {
                out.push((a, b));
            }
        }
        out
    }

    #[cfg(test)]
    pub fn num_pages(&self) -> usize {
        self.min.len()
    }

    #[cfg(test)]
    pub fn num_runs(&self) -> usize {
        self.runs.len()
    }
}

/// A `RowSelection` that selects whole page row-ranges (skips the gaps between).
pub fn selection_from_ranges(ranges: &[(u64, u64)]) -> RowSelection {
    let mut sel = Vec::with_capacity(ranges.len() * 2);
    let mut prev = 0u64;
    for &(a, b) in ranges {
        if a > prev {
            sel.push(RowSelector::skip((a - prev) as usize));
        }
        sel.push(RowSelector::select((b - a) as usize));
        prev = b;
    }
    RowSelection::from(sel)
}

/// A `RowSelection` that selects exactly the given (ascending) row offsets.
pub fn selection_from_offsets(offsets: &[u64]) -> RowSelection {
    let mut sel = Vec::with_capacity(offsets.len() * 2);
    let mut prev = 0u64;
    for &off in offsets {
        if off > prev {
            sel.push(RowSelector::skip((off - prev) as usize));
        }
        sel.push(RowSelector::select(1));
        prev = off + 1;
    }
    RowSelection::from(sel)
}

/// An [`AsyncFileReader`](parquet::arrow::async_reader::AsyncFileReader) over a
/// local file that coalesces nearby page byte-ranges into fewer, larger reads.
pub struct CoalescingAsyncReader {
    file: tokio::fs::File,
    counters: IoCounters,
    /// Merge ranges separated by <= `gap` bytes.
    gap: u64,
}

impl CoalescingAsyncReader {
    /// `gap_bytes` is the coalescing window (64 KiB is the measured sweet spot).
    pub fn new(file: tokio::fs::File, counters: IoCounters, gap_bytes: u64) -> Self {
        Self {
            file,
            counters,
            gap: gap_bytes,
        }
    }
}

impl parquet::arrow::async_reader::AsyncFileReader for CoalescingAsyncReader {
    fn get_bytes(
        &mut self,
        range: std::ops::Range<u64>,
    ) -> futures::future::BoxFuture<'_, parquet::errors::Result<Bytes>> {
        use futures::FutureExt;
        use tokio::io::{AsyncReadExt, AsyncSeekExt};
        async move {
            let len = (range.end - range.start) as usize;
            self.file
                .seek(std::io::SeekFrom::Start(range.start))
                .await
                .map_err(|e| parquet::errors::ParquetError::External(Box::new(e)))?;
            let mut buf = vec![0u8; len];
            self.file
                .read_exact(&mut buf)
                .await
                .map_err(|e| parquet::errors::ParquetError::External(Box::new(e)))?;
            self.counters.ops.fetch_add(1, Ordering::Relaxed);
            self.counters.bytes.fetch_add(len as u64, Ordering::Relaxed);
            Ok(Bytes::from(buf))
        }
        .boxed()
    }

    fn get_byte_ranges(
        &mut self,
        ranges: Vec<std::ops::Range<u64>>,
    ) -> futures::future::BoxFuture<'_, parquet::errors::Result<Vec<Bytes>>> {
        use futures::FutureExt;
        use tokio::io::{AsyncReadExt, AsyncSeekExt};
        let gap = self.gap;
        async move {
            let n = ranges.len();
            let mut idx: Vec<usize> = (0..n).collect();
            idx.sort_by_key(|&i| ranges[i].start);
            let mut out: Vec<Option<Bytes>> = vec![None; n];
            let mut i = 0;
            while i < n {
                let mut j = i;
                let span_start = ranges[idx[i]].start;
                let mut span_end = ranges[idx[i]].end;
                while j + 1 < n && ranges[idx[j + 1]].start <= span_end.saturating_add(gap) {
                    j += 1;
                    span_end = span_end.max(ranges[idx[j]].end);
                }
                let len = (span_end - span_start) as usize;
                self.file
                    .seek(std::io::SeekFrom::Start(span_start))
                    .await
                    .map_err(|e| parquet::errors::ParquetError::External(Box::new(e)))?;
                let mut buf = vec![0u8; len];
                self.file
                    .read_exact(&mut buf)
                    .await
                    .map_err(|e| parquet::errors::ParquetError::External(Box::new(e)))?;
                self.counters.ops.fetch_add(1, Ordering::Relaxed);
                self.counters.bytes.fetch_add(len as u64, Ordering::Relaxed);
                let span = Bytes::from(buf);
                for &k in &idx[i..=j] {
                    let off = (ranges[k].start - span_start) as usize;
                    let l = (ranges[k].end - ranges[k].start) as usize;
                    out[k] = Some(span.slice(off..off + l));
                }
                i = j + 1;
            }
            Ok(out
                .into_iter()
                .map(|x| x.expect("every range filled"))
                .collect())
        }
        .boxed()
    }

    fn get_metadata<'a>(
        &'a mut self,
        _options: Option<&'a parquet::arrow::arrow_reader::ArrowReaderOptions>,
    ) -> futures::future::BoxFuture<
        'a,
        parquet::errors::Result<Arc<parquet::file::metadata::ParquetMetaData>>,
    > {
        use futures::FutureExt;
        // Never called: the reader is always built with preloaded metadata
        // (`new_with_metadata`), so this path is unused.
        async move {
            Err(parquet::errors::ParquetError::General(
                "CoalescingAsyncReader::get_metadata is unused (metadata is preloaded)".into(),
            ))
        }
        .boxed()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn page_dir_resolves_within_two_runs() {
        // Two tier runs: pages 0..2 (warm), then min resets at page 2 (cold).
        let pd = PageDir::from_parts(
            vec![10, 30, 5, 40],  // min
            vec![29, 50, 39, 60], // max
            vec![0, 100, 200, 300],
            vec![100, 200, 300, 400],
        );
        assert_eq!(pd.num_pages(), 4);
        assert_eq!(pd.num_runs(), 2); // split at the tier boundary (5 < 30)

        // 35 lives in page 1 (warm run, 30..50) AND page 2 (cold run, 5..39) —
        // both runs must be searched. Pages 1 and 2 are row-adjacent
        // (100..200, 200..300), so they merge into one range.
        let r = pd.resolve_ranges(&[35]);
        assert_eq!(r, vec![(100, 300)]);

        // 45 is in page 1 (30..50) and page 3 (40..60); those are not adjacent.
        let r = pd.resolve_ranges(&[45]);
        assert_eq!(r, vec![(100, 200), (300, 400)]);
    }

    #[test]
    fn page_dir_spans_value_across_page_boundary() {
        // A run of identical key values (50) straddles the page boundary: page 0
        // ends with 50s (min=10,max=50), page 1 begins with 50s (min=50,max=90).
        // Both pages must be returned — else the rows of 50 in page 1 are dropped
        // (the real chr1 multi-allelic under-count bug).
        let pd = PageDir::from_parts(vec![10, 50], vec![50, 90], vec![0, 100], vec![100, 200]);
        assert_eq!(pd.num_runs(), 1);
        // Adjacent pages → merged into one contiguous range covering both.
        assert_eq!(pd.resolve_ranges(&[50]), vec![(0, 200)]);
        // A value only in the first page still returns just that page.
        assert_eq!(pd.resolve_ranges(&[10]), vec![(0, 100)]);
        // A value only in the second page returns just that page.
        assert_eq!(pd.resolve_ranges(&[90]), vec![(100, 200)]);
    }

    #[test]
    fn page_dir_spans_value_across_three_pages() {
        // 50 occupies the tail of page 0, all of page 1, and the head of page 2.
        let pd = PageDir::from_parts(
            vec![10, 50, 50],
            vec![50, 50, 90],
            vec![0, 100, 200],
            vec![100, 200, 300],
        );
        assert_eq!(pd.resolve_ranges(&[50]), vec![(0, 300)]);
    }

    #[test]
    fn page_dir_finds_warm_rows_in_a_tier_straddling_page() {
        // Reproduces the chr16:90202097 miss: the warm tier (pages 0,1 ascending)
        // is followed by the cold tier, but the boundary lands INSIDE a page, so
        // page 2 holds the warm tail (keys up to 300) AND the cold head (min 1).
        // Its min resets to 1 (→ run boundary) while its max stays high (300),
        // which is out of order vs the following pure-cold page 3 (max 50).
        let pd = PageDir::from_parts(
            vec![10, 91, 1, 2],     // min  (reset at page 2 → run split)
            vec![90, 150, 300, 50], // max  (page 2 max=300 is out of order vs page 3=50)
            vec![0, 100, 200, 300],
            vec![100, 200, 300, 400],
        );
        assert_eq!(pd.num_runs(), 2);
        // 200 is a warm key that lives ONLY in the straddling page 2 (its warm
        // tail): > page 1's max (150) and > pure-cold page 3's max (50). The
        // monotonic binary search skips page 2; the run-first-page bracket test
        // must recover it.
        // Without the run-first-page bracket test this returns [] (the bug).
        assert_eq!(pd.resolve_ranges(&[200]), vec![(200, 300)]);
    }

    #[test]
    fn page_dir_merges_adjacent_pages() {
        // Single run, pages 0 and 1 are adjacent in row space (0..100, 100..200).
        let pd = PageDir::from_parts(vec![0, 100], vec![99, 199], vec![0, 100], vec![100, 200]);
        assert_eq!(pd.num_runs(), 1);
        let r = pd.resolve_ranges(&[50, 150]);
        assert_eq!(r, vec![(0, 200)]); // merged into one contiguous range
    }

    #[test]
    fn selection_from_offsets_skips_and_selects_singletons() {
        let sel = selection_from_offsets(&[2, 5]);
        // skip 2, select 1, skip 2, select 1
        let rows: Vec<RowSelector> = sel.into();
        assert_eq!(rows.len(), 4);
        assert!(rows[0].skip);
        assert!(!rows[1].skip && rows[1].row_count == 1);
        assert!(rows[2].skip && rows[2].row_count == 2);
        assert!(!rows[3].skip && rows[3].row_count == 1);
    }
}
