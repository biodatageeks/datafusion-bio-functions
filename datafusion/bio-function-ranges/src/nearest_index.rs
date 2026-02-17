use std::cmp::Ordering;

use coitrees::{COITree, Interval, IntervalTree};
use fnv::FnvHashMap;

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
    by_end: Vec<Interval<Position>>,
    by_position: FnvHashMap<Position, IntervalMeta>,
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

        let mut by_end = by_start.clone();
        by_end.sort_by(|a, b| {
            a.last
                .cmp(&b.last)
                .then_with(|| a.first.cmp(&b.first))
                .then_with(|| a.metadata.cmp(&b.metadata))
        });

        let by_position = by_start
            .iter()
            .map(|r| {
                (
                    r.metadata,
                    IntervalMeta {
                        start: r.first,
                        end: r.last,
                        position: r.metadata,
                    },
                )
            })
            .collect::<FnvHashMap<_, _>>();

        let tree = COITree::new(by_start.iter());

        Self {
            tree,
            by_start,
            by_end,
            by_position,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.by_start.is_empty()
    }

    pub fn nearest_one(&self, start: i32, end: i32, include_overlaps: bool) -> Option<Position> {
        if self.is_empty() {
            return None;
        }

        if include_overlaps {
            let mut first_overlap = None;
            self.tree.query(start, end, |node| {
                if first_overlap.is_none() {
                    first_overlap = Some(extract_position(node));
                }
            });
            if first_overlap.is_some() {
                return first_overlap;
            }
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

        let mut seen = FnvHashMap::<Position, ()>::default();

        if include_overlaps {
            let mut overlaps = Vec::<IntervalMeta>::new();
            self.tree.query(start, end, |node| {
                let pos = extract_position(node);
                if let Some(meta) = self.by_position.get(&pos) {
                    overlaps.push(*meta);
                }
            });
            overlaps.sort_by(cmp_interval_meta);
            overlaps.dedup_by_key(|m| m.position);
            for m in overlaps {
                if out.len() == k {
                    return;
                }
                if seen.insert(m.position, ()).is_none() {
                    out.push(m.position);
                }
            }
        }

        if out.len() == k {
            return;
        }

        let mut left_idx = self.by_end.partition_point(|iv| iv.last < start);
        let mut right_idx = self.by_start.partition_point(|iv| iv.first <= end);

        while out.len() < k {
            let left_meta = if left_idx > 0 {
                self.interval_meta(&self.by_end[left_idx - 1])
            } else {
                None
            };
            let right_meta = if right_idx < self.by_start.len() {
                self.interval_meta(&self.by_start[right_idx])
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

            if seen.insert(next.position, ()).is_none() {
                out.push(next.position);
            }
        }
    }

    fn nearest_non_overlap_one(&self, start: i32, end: i32) -> Option<IntervalMeta> {
        let left_idx = self.by_end.partition_point(|iv| iv.last < start);
        let right_idx = self.by_start.partition_point(|iv| iv.first <= end);

        let left = if left_idx > 0 {
            self.interval_meta(&self.by_end[left_idx - 1])
        } else {
            None
        };
        let right = if right_idx < self.by_start.len() {
            self.interval_meta(&self.by_start[right_idx])
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

    fn interval_meta(&self, iv: &Interval<Position>) -> Option<IntervalMeta> {
        self.by_position.get(&iv.metadata).copied()
    }
}

fn cmp_interval_meta(a: &IntervalMeta, b: &IntervalMeta) -> Ordering {
    a.start
        .cmp(&b.start)
        .then_with(|| a.end.cmp(&b.end))
        .then_with(|| a.position.cmp(&b.position))
}

fn candidate_distance(start: i32, end: i32, m: &IntervalMeta) -> i64 {
    if end < m.start {
        i64::from(m.start) - i64::from(end)
    } else if m.end < start {
        i64::from(start) - i64::from(m.end)
    } else {
        0
    }
}

fn cmp_candidate(start: i32, end: i32, a: &IntervalMeta, b: &IntervalMeta) -> Ordering {
    let ad = candidate_distance(start, end, a);
    let bd = candidate_distance(start, end, b);
    ad.cmp(&bd).then_with(|| cmp_interval_meta(a, b))
}

/// unoptimized on Linux x64 (without target-cpu=native)
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
fn extract_position(node: &coitrees::IntervalNode<Position, u32>) -> Position {
    node.metadata
}

/// for Apple Intel, Apple M1+(both optimized and not) and optimized (target-cpu=native) on Linux x64 and Linux aarch64
#[cfg(any(
    all(target_os = "macos", target_arch = "aarch64"),
    all(target_os = "macos", target_arch = "x86_64", target_feature = "avx"),
    all(target_os = "linux", target_arch = "x86_64", target_feature = "avx"),
    all(target_os = "linux", target_arch = "aarch64"),
    all(target_os = "windows", target_arch = "x86_64", target_feature = "avx")
))]
fn extract_position(node: &coitrees::Interval<&Position>) -> Position {
    *node.metadata
}

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
    fn nearest_k_empty_returns_no_candidates() {
        let idx = mk(&[]);
        let mut out = Vec::new();
        idx.nearest_k(10, 20, 3, true, &mut out);
        assert!(out.is_empty());
    }
}
