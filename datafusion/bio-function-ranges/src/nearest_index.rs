use std::cmp::Ordering;
use std::sync::OnceLock;

use ahash::AHashSet;
use coitrees::{COITree, Interval, IntervalTree};

pub type Position = usize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntervalRecord {
    pub start: i32,
    pub end: i32,
    pub position: Position,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct IntervalMeta {
    start: i32,
    end: i32,
    position: Position,
}

pub struct NearestIntervalIndex {
    tree: COITree<Position, u32>,
    by_start: Vec<Interval<Position>>,
    /// Eagerly allocated because the k=1 include-overlaps hot path probes it on
    /// every query to find the first deterministic overlap in O(log n).
    prefix_max_end: Vec<i32>,
    /// Lazily initialized: only allocated when a non-overlap search is needed.
    /// For the common k=1 include_overlaps=true hot path where most queries
    /// find overlaps via the COITree, this avoids ~24MB per contig of cache
    /// pressure from a sorted copy that would rarely be accessed.
    by_end: OnceLock<Vec<Interval<Position>>>,
}

impl std::fmt::Debug for NearestIntervalIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NearestIntervalIndex")
            .field("size", &self.by_start.len())
            .finish()
    }
}

impl NearestIntervalIndex {
    pub fn from_records(records: Vec<IntervalRecord>) -> Self {
        let mut by_start = records
            .iter()
            .map(|r| Interval::new(r.start, r.end, r.position))
            .collect::<Vec<_>>();
        by_start.sort_by(|a, b| {
            a.first
                .cmp(&b.first)
                .then_with(|| a.last.cmp(&b.last))
                .then_with(|| a.metadata.cmp(&b.metadata))
        });

        let tree = COITree::new(by_start.iter());
        let mut prefix_max_end = Vec::with_capacity(by_start.len());
        let mut max_end = i32::MIN;
        for interval in &by_start {
            max_end = max_end.max(interval.last);
            prefix_max_end.push(max_end);
        }

        Self {
            tree,
            by_start,
            prefix_max_end,
            by_end: OnceLock::new(),
        }
    }

    /// Returns the by-end sorted intervals, creating them lazily on first use.
    fn by_end(&self) -> &[Interval<Position>] {
        self.by_end.get_or_init(|| {
            let mut by_end = self.by_start.clone();
            by_end.sort_by(|a, b| {
                a.last
                    .cmp(&b.last)
                    .then_with(|| a.first.cmp(&b.first))
                    .then_with(|| a.metadata.cmp(&b.metadata))
            });
            by_end
        })
    }

    pub fn is_empty(&self) -> bool {
        self.by_start.is_empty()
    }

    /// Bounds for comparing a query against the raw stored interval bounds.
    ///
    /// Returns `None` for a genuinely empty query. Under `half_open` the query
    /// is narrowed by one on each side so that inclusive-style comparisons
    /// against raw half-open bounds implement half-open overlap; for a one-base
    /// query that narrowed pair is inverted, which is correct for every
    /// comparison below and must not be mistaken for an empty range.
    fn comparison_bounds(start: i32, end: i32, half_open: bool) -> Option<(i32, i32)> {
        if half_open {
            if end <= start {
                return None;
            }
            Some((start + 1, end - 1))
        } else {
            if end < start {
                return None;
            }
            Some((start, end))
        }
    }

    pub fn nearest_one(
        &self,
        start: i32,
        end: i32,
        include_overlaps: bool,
        half_open: bool,
    ) -> Option<Position> {
        if self.is_empty() {
            return None;
        }
        let (lo, hi) = Self::comparison_bounds(start, end, half_open)?;
        self.nearest_one_by_bounds(lo, hi, include_overlaps)
    }

    /// [`Self::nearest_one`] for callers that already hold comparison bounds
    /// (see [`Self::comparison_bounds`]) rather than raw interval coordinates.
    pub fn nearest_one_by_bounds(
        &self,
        lo: i32,
        hi: i32,
        include_overlaps: bool,
    ) -> Option<Position> {
        if self.is_empty() {
            return None;
        }

        if include_overlaps && let Some(best) = self.first_overlap_by_start(lo, hi) {
            return Some(best.position);
        }

        self.nearest_non_overlap_one(lo, hi).map(|m| m.position)
    }

    pub fn nearest_k(
        &self,
        start: i32,
        end: i32,
        k: usize,
        include_overlaps: bool,
        half_open: bool,
        out: &mut Vec<Position>,
    ) {
        if k == 0 || self.is_empty() {
            return;
        }

        if k == 1 {
            if let Some(p) = self.nearest_one(start, end, include_overlaps, half_open) {
                out.push(p);
            }
            return;
        }

        let Some((start, end)) = Self::comparison_bounds(start, end, half_open) else {
            return;
        };

        let mut seen = AHashSet::<Position>::default();

        if include_overlaps {
            let mut overlaps = Vec::<IntervalMeta>::new();
            self.tree.query(start, end, |node| {
                overlaps.push(extract_coitree_meta(node));
            });
            overlaps.sort_by(cmp_interval_meta);
            for m in overlaps {
                if out.len() == k {
                    return;
                }
                if seen.insert(m.position) {
                    out.push(m.position);
                }
            }
        }

        if out.len() == k {
            return;
        }

        let by_end = self.by_end();
        let mut left_idx = by_end.partition_point(|iv| iv.last < start);
        let mut right_idx = self.by_start.partition_point(|iv| iv.first <= end);

        while out.len() < k {
            let left_meta = if left_idx > 0 {
                Some(interval_meta_from(&by_end[left_idx - 1]))
            } else {
                None
            };
            let right_meta = if right_idx < self.by_start.len() {
                Some(interval_meta_from(&self.by_start[right_idx]))
            } else {
                None
            };

            let next = match (left_meta, right_meta) {
                (None, None) => break,
                (Some(left), None) => {
                    left_idx -= 1;
                    left
                }
                (None, Some(right)) => {
                    right_idx += 1;
                    right
                }
                (Some(left), Some(right)) => {
                    let take_left = cmp_candidate(start, end, &left, &right).is_le();
                    if take_left {
                        left_idx -= 1;
                        left
                    } else {
                        right_idx += 1;
                        right
                    }
                }
            };

            // Comparison bounds again, so this is the plain overlap test: see
            // cmp_candidate on why the convention is immaterial here.
            if !include_overlaps && candidate_distance(start, end, next.start, next.end, true) == 0
            {
                continue;
            }

            if seen.insert(next.position) {
                out.push(next.position);
            }
        }
    }

    fn nearest_non_overlap_one(&self, start: i32, end: i32) -> Option<IntervalMeta> {
        let by_end = self.by_end();
        let left_idx = by_end.partition_point(|iv| iv.last < start);
        let right_idx = self.by_start.partition_point(|iv| iv.first <= end);

        let left = if left_idx > 0 {
            Some(interval_meta_from(&by_end[left_idx - 1]))
        } else {
            None
        };
        let right = if right_idx < self.by_start.len() {
            Some(interval_meta_from(&self.by_start[right_idx]))
        } else {
            None
        };

        match (left, right) {
            (None, None) => None,
            (Some(left), None) => Some(left),
            (None, Some(right)) => Some(right),
            (Some(left), Some(right)) => {
                if cmp_candidate(start, end, &left, &right).is_le() {
                    Some(left)
                } else {
                    Some(right)
                }
            }
        }
    }

    /// `start`/`end` are comparison bounds from [`Self::comparison_bounds`], so
    /// `end < start` here means a one-base half-open query, not an empty one --
    /// every comparison below handles that correctly. Emptiness is already
    /// filtered out by the caller.
    fn first_overlap_by_start(&self, start: i32, end: i32) -> Option<IntervalMeta> {
        let prefix_len = self.by_start.partition_point(|iv| iv.first <= end);
        if prefix_len == 0 || self.prefix_max_end[prefix_len - 1] < start {
            return None;
        }

        let idx = self.prefix_max_end[..prefix_len].partition_point(|&max_end| max_end < start);
        Some(interval_meta_from(&self.by_start[idx]))
    }
}

fn interval_meta_from(iv: &Interval<Position>) -> IntervalMeta {
    IntervalMeta {
        start: iv.first,
        end: iv.last,
        position: iv.metadata,
    }
}

fn cmp_interval_meta(a: &IntervalMeta, b: &IntervalMeta) -> Ordering {
    a.start
        .cmp(&b.start)
        .then_with(|| a.end.cmp(&b.end))
        .then_with(|| a.position.cmp(&b.position))
}

/// Number of bases separating a query from a candidate interval, 0 if they meet.
///
/// `half_open` is `true` for 0-based `[start, end)`, where the exclusive end is
/// already one past the last covered base and the raw coordinate difference is
/// the gap. Under 1-based inclusive coordinates both endpoints are covered, so
/// the gap is one smaller than that difference.
pub fn candidate_distance(
    query_start: i32,
    query_end: i32,
    iv_start: i32,
    iv_end: i32,
    half_open: bool,
) -> i64 {
    let shared_endpoint = i64::from(!half_open);
    if query_end < iv_start {
        i64::from(iv_start) - i64::from(query_end) - shared_endpoint
    } else if iv_end < query_start {
        i64::from(query_start) - i64::from(iv_end) - shared_endpoint
    } else {
        0
    }
}

fn cmp_candidate(start: i32, end: i32, a: &IntervalMeta, b: &IntervalMeta) -> Ordering {
    // Ordering only, on comparison bounds rather than raw ones. The convention
    // shifts every non-zero distance by the same constant, so it cannot change
    // the ranking -- `true` here keeps this a plain coordinate difference.
    let ad = candidate_distance(start, end, a.start, a.end, true);
    let bd = candidate_distance(start, end, b.start, b.end, true);
    ad.cmp(&bd).then_with(|| cmp_interval_meta(a, b))
}

// COITree node accessors — grouped per platform so both functions share the
// same cfg predicate; adding a new target only requires updating one block.

/// x86_64 without AVX (unoptimized builds on Linux/macOS/Windows).
/// COITree uses `IntervalNode` in this configuration.
#[cfg(any(
    all(
        target_os = "linux",
        target_arch = "x86_64",
        not(target_feature = "avx")
    ),
    all(
        target_os = "macos",
        target_arch = "x86_64",
        not(target_feature = "avx")
    ),
    all(
        target_os = "windows",
        target_arch = "x86_64",
        not(target_feature = "avx")
    ),
))]
mod coitree_extract {
    use super::*;

    pub(crate) fn extract_coitree_position(
        node: &coitrees::IntervalNode<Position, u32>,
    ) -> Position {
        node.metadata
    }

    pub(super) fn extract_coitree_meta(
        node: &coitrees::IntervalNode<Position, u32>,
    ) -> IntervalMeta {
        IntervalMeta {
            start: node.first,
            end: node.last,
            position: node.metadata,
        }
    }
}

/// aarch64 (Apple M1+, Linux ARM) and x86_64 with AVX (optimized builds).
/// COITree uses `Interval<&T>` in this configuration.
#[cfg(any(
    all(target_os = "macos", target_arch = "aarch64"),
    all(target_os = "macos", target_arch = "x86_64", target_feature = "avx"),
    all(target_os = "linux", target_arch = "x86_64", target_feature = "avx"),
    all(target_os = "linux", target_arch = "aarch64"),
    all(target_os = "windows", target_arch = "x86_64", target_feature = "avx"),
))]
mod coitree_extract {
    use super::*;

    pub(crate) fn extract_coitree_position(node: &coitrees::Interval<&Position>) -> Position {
        *node.metadata
    }

    pub(super) fn extract_coitree_meta(node: &coitrees::Interval<&Position>) -> IntervalMeta {
        IntervalMeta {
            start: node.first,
            end: node.last,
            position: *node.metadata,
        }
    }
}

use coitree_extract::extract_coitree_meta;
pub(crate) use coitree_extract::extract_coitree_position;

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(records: &[(i32, i32, usize)]) -> NearestIntervalIndex {
        NearestIntervalIndex::from_records(
            records
                .iter()
                .map(|(start, end, position)| IntervalRecord {
                    start: *start,
                    end: *end,
                    position: *position,
                })
                .collect(),
        )
    }

    #[test]
    fn comparison_bounds_reject_only_genuinely_empty_queries() {
        // One-base half-open: inverted bounds, but not empty.
        assert_eq!(
            NearestIntervalIndex::comparison_bounds(17, 18, true),
            Some((18, 17))
        );
        // Zero-length half-open: empty.
        assert_eq!(NearestIntervalIndex::comparison_bounds(17, 17, true), None);
        // 1-based inclusive passes through untouched.
        assert_eq!(
            NearestIntervalIndex::comparison_bounds(17, 17, false),
            Some((17, 17))
        );
        assert_eq!(NearestIntervalIndex::comparison_bounds(18, 17, false), None);
    }

    #[test]
    fn nearest_one_finds_enclosing_interval_for_one_base_half_open_query() {
        let idx = mk(&[(9, 28, 0)]);
        // [17,18) sits inside [9,28); the enclosing interval must be found.
        assert_eq!(idx.nearest_one(17, 18, true, true), Some(0));
        assert_eq!(idx.nearest_one(15, 20, true, true), Some(0));
    }

    #[test]
    fn candidate_distance_counts_bases_between() {
        // 0-based [10,20) to [30,40): bases 20..29, i.e. 10.
        assert_eq!(candidate_distance(10, 20, 30, 40, true), 10);
        // 1-based [10,20] to [31,40]: bases 21..30, i.e. 10.
        assert_eq!(candidate_distance(10, 20, 31, 40, false), 10);
        // 1-based neighbours that touch have nothing between them.
        assert_eq!(candidate_distance(10, 20, 21, 30, false), 0);
        // Overlapping is 0 under either convention.
        assert_eq!(candidate_distance(10, 20, 15, 25, true), 0);
        assert_eq!(candidate_distance(10, 20, 15, 25, false), 0);
    }

    #[test]
    fn nearest_one_finds_true_nearest_not_only_adjacent_start() {
        // Query 1100 is closer to [0,1000] than to [900,905], but
        // [0,1000] is not adjacent by start around binary search pivot.
        let idx = mk(&[(0, 1000, 0), (900, 905, 1), (2000, 2010, 2)]);
        let nearest = idx.nearest_one(1100, 1100, true, false);
        assert_eq!(nearest, Some(0));
    }

    #[test]
    fn nearest_one_prefers_overlap_when_requested() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1)]);
        let nearest = idx.nearest_one(12, 12, true, false);
        assert_eq!(nearest, Some(0));
    }

    #[test]
    fn nearest_one_overlap_uses_deterministic_coordinate_order() {
        let idx = mk(&[(20, 30, 10), (10, 40, 20), (15, 25, 30)]);
        let nearest = idx.nearest_one(21, 21, true, false);
        // all overlap; deterministic tie-break is (start, end, position)
        assert_eq!(nearest, Some(20));
    }

    #[test]
    fn nearest_one_overlap_skips_early_non_overlapping_end() {
        let idx = mk(&[(10, 20, 0), (10, 50, 1), (30, 40, 2)]);
        let nearest = idx.nearest_one(30, 30, true, false);
        assert_eq!(nearest, Some(1));
    }

    #[test]
    fn nearest_k_non_overlap_returns_expected_order() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1), (50, 60, 2), (70, 80, 3)]);
        let mut out = Vec::new();
        idx.nearest_k(22, 22, 3, false, false, &mut out);
        assert_eq!(out, vec![0, 1, 2]);
    }

    #[test]
    fn nearest_k_includes_overlaps_then_fills_nearest() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1), (50, 60, 2)]);
        let mut out = Vec::new();
        idx.nearest_k(35, 35, 2, true, false, &mut out);
        assert_eq!(out, vec![1, 0]);
    }

    #[test]
    fn nearest_k_non_overlap_excludes_overlapping_candidates() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1), (50, 60, 2)]);
        let mut out = Vec::new();
        idx.nearest_k(35, 35, 2, false, false, &mut out);
        assert_eq!(out, vec![0, 2]);
    }

    #[test]
    fn nearest_k_empty_returns_no_candidates() {
        let idx = mk(&[]);
        let mut out = Vec::new();
        idx.nearest_k(10, 20, 3, true, false, &mut out);
        assert!(out.is_empty());
    }
}
