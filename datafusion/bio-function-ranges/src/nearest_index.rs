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

    pub fn nearest_one(&self, start: i32, end: i32, include_overlaps: bool) -> Option<Position> {
        if self.is_empty() {
            return None;
        }

        if include_overlaps && let Some(best) = self.first_overlap_by_start(start, end) {
            return Some(best.position);
        }

        self.nearest_non_overlap_one(start, end).map(|m| m.position)
    }

    pub fn nearest_k(
        &self,
        start: i32,
        end: i32,
        k: usize,
        include_overlaps: bool,
        out: &mut Vec<Position>,
    ) {
        if k == 0 || self.is_empty() {
            return;
        }

        if k == 1 {
            if let Some(p) = self.nearest_one(start, end, include_overlaps) {
                out.push(p);
            }
            return;
        }

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

            if !include_overlaps && candidate_distance(start, end, next.start, next.end) == 0 {
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

    fn first_overlap_by_start(&self, start: i32, end: i32) -> Option<IntervalMeta> {
        if end < start {
            return None;
        }

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

pub fn candidate_distance(query_start: i32, query_end: i32, iv_start: i32, iv_end: i32) -> i64 {
    if query_end < iv_start {
        i64::from(iv_start) - i64::from(query_end)
    } else if iv_end < query_start {
        i64::from(query_start) - i64::from(iv_end)
    } else {
        0
    }
}

fn cmp_candidate(start: i32, end: i32, a: &IntervalMeta, b: &IntervalMeta) -> Ordering {
    let ad = candidate_distance(start, end, a.start, a.end);
    let bd = candidate_distance(start, end, b.start, b.end);
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
    fn nearest_one_finds_true_nearest_not_only_adjacent_start() {
        // Query 1100 is closer to [0,1000] than to [900,905], but
        // [0,1000] is not adjacent by start around binary search pivot.
        let idx = mk(&[(0, 1000, 0), (900, 905, 1), (2000, 2010, 2)]);
        let nearest = idx.nearest_one(1100, 1100, true);
        assert_eq!(nearest, Some(0));
    }

    #[test]
    fn nearest_one_prefers_overlap_when_requested() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1)]);
        let nearest = idx.nearest_one(12, 12, true);
        assert_eq!(nearest, Some(0));
    }

    #[test]
    fn nearest_one_overlap_uses_deterministic_coordinate_order() {
        let idx = mk(&[(20, 30, 10), (10, 40, 20), (15, 25, 30)]);
        let nearest = idx.nearest_one(21, 21, true);
        // all overlap; deterministic tie-break is (start, end, position)
        assert_eq!(nearest, Some(20));
    }

    #[test]
    fn nearest_one_overlap_skips_early_non_overlapping_end() {
        let idx = mk(&[(10, 20, 0), (10, 50, 1), (30, 40, 2)]);
        let nearest = idx.nearest_one(30, 30, true);
        assert_eq!(nearest, Some(1));
    }

    #[test]
    fn nearest_k_non_overlap_returns_expected_order() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1), (50, 60, 2), (70, 80, 3)]);
        let mut out = Vec::new();
        idx.nearest_k(22, 22, 3, false, &mut out);
        assert_eq!(out, vec![0, 1, 2]);
    }

    #[test]
    fn nearest_k_includes_overlaps_then_fills_nearest() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1), (50, 60, 2)]);
        let mut out = Vec::new();
        idx.nearest_k(35, 35, 2, true, &mut out);
        assert_eq!(out, vec![1, 0]);
    }

    #[test]
    fn nearest_k_non_overlap_excludes_overlapping_candidates() {
        let idx = mk(&[(10, 20, 0), (30, 40, 1), (50, 60, 2)]);
        let mut out = Vec::new();
        idx.nearest_k(35, 35, 2, false, &mut out);
        assert_eq!(out, vec![0, 2]);
    }

    #[test]
    fn nearest_k_empty_returns_no_candidates() {
        let idx = mk(&[]);
        let mut out = Vec::new();
        idx.nearest_k(10, 20, 3, true, &mut out);
        assert!(out.is_empty());
    }

    // =========================================================================
    // V2 PORT — TranscriptTree.t
    // =========================================================================
    //
    // Sztywno 1:1 Rust analogues of subtests from ensembl-vep/t/TranscriptTree.t.
    // See porting-tests/detailed_plans/TranscriptTree.md for the full mapping
    // table. Coverage: 20 unit-port / 13 architectural-no-analogue /
    // 3 blocked-future-work / 36 substantive rows = 56%.
    //
    // The COITree direct tests live here because they are pure interval-tree
    // overlap / nearest assertions that do not touch any VEP-specific types.
    // The PreparedContext-surface tests live in transcript_consequence.rs.

    use coitrees::{COITree, GenericInterval, Interval, IntervalTree};

    /// Builds a `COITree<usize, u32>` from a flat list of (first, last, metadata)
    /// triples. Mirrors the construction used by `PreparedFeatureIndex::new` in
    /// `bio-function-vep::transcript_consequence` so the port tests exercise
    /// the same coitrees idiom vepyr uses in production.
    fn mk_coitree(records: &[(i32, i32, usize)]) -> COITree<usize, u32> {
        let intervals: Vec<Interval<usize>> = records
            .iter()
            .map(|(first, last, pos)| Interval::new(*first, *last, *pos))
            .collect();
        COITree::new(&intervals)
    }

    /// Wraps `tree.query` to return a sorted Vec of matched metadata values, so
    /// assertions don't depend on coitrees' internal callback order.
    fn fetch_hits(tree: &COITree<usize, u32>, start: i32, end: i32) -> Vec<usize> {
        let mut out = Vec::new();
        tree.query(start, end, |node| {
            out.push(*GenericInterval::<usize>::metadata(node));
        });
        out.sort_unstable();
        out
    }

    // SUBTEST #13 — Perl: insert(1,5); fetch(2,3) == [[1,5]] (centered overlap).
    // vepyr analogue: COITree::query(2, 3, ...) on a tree containing only
    // Interval::new(1, 5, 10) returns exactly metadata == 10.
    #[test]
    fn tt13_fetch_centered_overlap() {
        let tree = mk_coitree(&[(1, 5, 10)]);
        assert_eq!(fetch_hits(&tree, 2, 3), vec![10]);
    }

    // SUBTEST #14 — Perl: fetch(0,1) == [[1,5]] (left-boundary overlap).
    // coitrees uses inclusive [first, last]; query (0, 1) touches first=1.
    #[test]
    fn tt14_fetch_left_boundary() {
        let tree = mk_coitree(&[(1, 5, 10)]);
        assert_eq!(fetch_hits(&tree, 0, 1), vec![10]);
    }

    // SUBTEST #15 — Perl: fetch(0,6) == [[1,5]] (superset query overlap).
    #[test]
    fn tt15_fetch_superset_query() {
        let tree = mk_coitree(&[(1, 5, 10)]);
        assert_eq!(fetch_hits(&tree, 0, 6), vec![10]);
    }

    // SUBTEST #16 — Perl: fetch(4,6) == [[1,5]] (right-boundary overlap).
    #[test]
    fn tt16_fetch_right_boundary() {
        let tree = mk_coitree(&[(1, 5, 10)]);
        assert_eq!(fetch_hits(&tree, 4, 6), vec![10]);
    }

    // SUBTEST #17 — Perl: fetch(6,7) == [] (no overlap).
    #[test]
    fn tt17_fetch_no_overlap() {
        let tree = mk_coitree(&[(1, 5, 10)]);
        let hits = fetch_hits(&tree, 6, 7);
        assert!(hits.is_empty(), "expected no hits, got {hits:?}");
    }

    // SUBTEST #18 — architectural-no-analogue.
    // Perl: insert(1,5) + insert(4,8); fetch(2,3) == [[1,8]] (merged).
    // vepyr/coitrees does NOT merge overlapping inserts: after two inserts,
    // a query that overlaps both returns BOTH intervals. This is the desired
    // behavior for downstream NEAREST / per-transcript consequence calculation;
    // merging would lose transcript identity. See detailed_plan #18 / #19.

    // SUBTEST #19 — architectural-no-analogue.
    // Same justification as #18. coitrees returns 3 separate intervals after
    // 3 inserts; the Perl test's merged [[1,8],[9,10]] expectation is a
    // wrapper-side artifact of Perl's `insert()` (TranscriptTree.pm:155-185).

    // SUBTEST #24 — Perl: insert('chrobj',5,10,{foo=>'bar'});
    //   fetch(7,8) == [{foo=>'bar', s=>5, e=>10}].
    // vepyr analogue: build a COITree where each interval's metadata is an
    // index into a parallel `Vec<Payload>`; query at (7,8) returns metadata
    // 0; payloads[0] is the recovered "foo=bar" payload. This is exactly the
    // pattern used by PreparedFeatureIndex (transcript_consequence.rs:641):
    // metadata = usize index, payload = &T from a Vec<&T>.
    #[test]
    fn tt24_payload_recovery_by_index() {
        #[derive(Debug, PartialEq)]
        struct Payload {
            tag: &'static str,
        }
        let payloads = vec![Payload { tag: "bar" }];
        let tree = mk_coitree(&[(5, 10, 0)]);
        let hits = fetch_hits(&tree, 7, 8);
        assert_eq!(hits, vec![0]);
        assert_eq!(payloads[hits[0]], Payload { tag: "bar" });
    }

    // SUBTEST #25 — Perl: second insert('chrobj',6,12,{goo=>'car'}) does NOT
    // merge; fetch(7,8) returns BOTH payloads (Perl insert-with-payload
    // skips the merge branch). vepyr/coitrees: 2 inserts → query returns 2
    // metadata values, both deterministically present.
    #[test]
    fn tt25_multi_payload_no_merge() {
        let tree = mk_coitree(&[(5, 10, 0), (6, 12, 1)]);
        let hits = fetch_hits(&tree, 7, 8);
        assert_eq!(hits, vec![0, 1]);
    }

    // SUBTEST #31 — Perl: _get_dist(1, 5, 6, 7) == 1 (non-overlap forward).
    // vepyr analogue: `candidate_distance` from this module returns gap >= 0
    // between two non-overlapping intervals.
    #[test]
    fn tt31_candidate_distance_non_overlap_forward() {
        assert_eq!(candidate_distance(1, 5, 6, 7), 1);
    }

    // SUBTEST #32 — Perl: _get_dist(6, 8, 2, 4) == 2 (non-overlap reverse).
    #[test]
    fn tt32_candidate_distance_non_overlap_reverse() {
        assert_eq!(candidate_distance(6, 8, 2, 4), 2);
    }

    // SUBTEST #33 — Perl: _get_dist(1, 5, 4, 5) == 0 (overlap).
    #[test]
    fn tt33_candidate_distance_overlap_zero() {
        assert_eq!(candidate_distance(1, 5, 4, 5), 0);
    }

    // SUBTEST #34 of TranscriptTree.t (reclassified architectural-no-analogue → unit-port
    // per Phase D 2026-05-27 PORT re-audit). Perl Set::IntervalTree accepts negative
    // coordinates as a side-effect of signed-int internals; the Perl test asserts that
    // `_get_dist(1, 5, -3, -1) == 2` produces the unsigned distance correctly.
    //
    // vepyr's `candidate_distance` is the analogue, signed by design.
    // verified via direct function call (no VEP 115 oracle needed; pure-fn arithmetic)
    #[test]
    fn tt34_candidate_distance_handles_negative_coords() {
        let d = candidate_distance(1, 5, -3, -1);
        assert_eq!(
            d, 2,
            "candidate_distance(1, 5, -3, -1) must equal 2 (negative-coord case)"
        );
    }

    // SUBTEST #35 — Perl: nearest('nearest', 1, 5) == [] (empty tree).
    // vepyr: NearestIntervalIndex built from an empty record list returns
    // None on nearest_one for any query (overlap or non-overlap).
    #[test]
    fn tt35_nearest_on_empty_index_returns_none() {
        let idx = mk(&[]);
        assert_eq!(idx.nearest_one(1, 5, true), None);
        assert_eq!(idx.nearest_one(1, 5, false), None);
    }

    // SUBTEST #36 — Perl: nearest('nearest', 6, 6) returns the LEFT interval
    // (1,5) after inserting (1,5) and (11,16). Query at point 6 is distance 1
    // from (1,5) (right edge) and distance 5 from (11,16) (left edge).
    #[test]
    fn tt36_nearest_one_picks_left() {
        let idx = mk(&[(1, 5, 0), (11, 16, 1)]);
        assert_eq!(idx.nearest_one(6, 6, false), Some(0));
    }

    // SUBTEST #37 — Perl: nearest('nearest', 9, 9) returns the RIGHT interval
    // (11,16). Distance to (1,5) is 4; distance to (11,16) is 2.
    #[test]
    fn tt37_nearest_one_picks_right() {
        let idx = mk(&[(1, 5, 0), (11, 16, 1)]);
        assert_eq!(idx.nearest_one(9, 9, false), Some(1));
    }

    // SUBTEST #38 — Perl: nearest('nearest', 8, 8) returns BOTH intervals
    // (tie at distance 3 from each). vepyr analogue: `nearest_k(8, 8, 2,
    // false, &mut out)` returns both positions, with deterministic tie-break
    // by (start, end, position).
    #[test]
    fn tt38_nearest_k_returns_tied_pair() {
        let idx = mk(&[(1, 5, 0), (11, 16, 1)]);
        let mut out = Vec::new();
        idx.nearest_k(8, 8, 2, false, &mut out);
        assert_eq!(out.len(), 2, "expected both ties; got {out:?}");
        // Deterministic order: left interval (start=1) precedes right (start=11).
        assert_eq!(out, vec![0, 1]);
    }

    // SUBTEST #39 — Perl: nearest('nearest', 7, 9) (non-point query) returns
    // BOTH intervals — both tied at distance 2 from query interval [7,9].
    #[test]
    fn tt39_nearest_k_non_point_query_returns_tied_pair() {
        let idx = mk(&[(1, 5, 0), (11, 16, 1)]);
        let mut out = Vec::new();
        idx.nearest_k(7, 9, 2, false, &mut out);
        assert_eq!(out.len(), 2, "expected both ties; got {out:?}");
        assert_eq!(out, vec![0, 1]);
    }

    // =========================================================================
    // blocked-future-work stubs (TranscriptTree.t)
    // =========================================================================
    //
    // SUBTEST #11 — Perl: valid_chromosomes([21]) setter override.
    // Missing API: PreparedContext::override_valid_chromosomes(set: HashSet<String>).
    // Future-work entry: porting-tests/future-work-vepyr.md
    //   "PreparedContext::override_valid_chromosomes".
    //
    // // #[test]
    // // fn tt11_override_valid_chromosomes() {
    // //     let mut ctx = PreparedContext::new(...);
    // //     ctx.override_valid_chromosomes(HashSet::from(["21".to_string()]));
    // //     assert_eq!(
    // //         ctx.valid_chromosomes(),
    // //         &HashSet::from(["21".to_string()])
    // //     );
    // // }
    //
    // SUBTEST #20 + #21 + #23 — chromosome-synonyms alias table.
    // Missing API: PreparedContext::set_chrom_synonyms(syns: HashMap<String, String>)
    //   + canonical_chrom(alias: &str) -> &str, routed through every tx_trees.get(...)
    //   call site in transcript_consequence.rs.
    // Future-work entry: porting-tests/future-work-vepyr.md
    //   "PreparedContext chromosome-synonyms alias table".
    //
    // // #[test]
    // // fn tt20_load_chrom_synonyms_from_file() {
    // //     let mut ctx = PreparedContext::new(...);
    // //     let syns = HashMap::from([("NC_000021.9".to_string(), "21".to_string())]);
    // //     ctx.set_chrom_synonyms(syns);
    // //     // assert that a query keyed "NC_000021.9" reaches the "21" tree.
    // // }
    //
    // // #[test]
    // // fn tt21_synonym_NC_000021_9_resolves_to_canonical_21() {
    // //     // Build ctx with a tx on "21"; set synonyms; query "NC_000021.9" hits.
    // // }
    //
    // // #[test]
    // // fn tt23_query_21_hits_tree_keyed_chr21() {
    // //     // Build ctx with a tx keyed "chr21"; query "21" should hit. Requires the
    // //     // alias table because normalize_chrom only STRIPS chr, never ADDS it.
    // // }
}
