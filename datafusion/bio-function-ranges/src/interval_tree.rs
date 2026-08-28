use std::cmp::{max, min};
use std::sync::Arc;

use ahash::AHashMap;
use coitrees::{COITree, Interval, IntervalTree};
use datafusion::arrow::array::{Int64Array, RecordBatch};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use futures::StreamExt;
use futures::stream::BoxStream;

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;

type IntervalHashMap = AHashMap<String, Vec<Interval<()>>>;

pub struct CountOverlapIndex {
    starts: Vec<i32>,
    ends: Vec<i32>,
}

impl CountOverlapIndex {
    fn new(intervals: Vec<Interval<()>>) -> Self {
        let mut starts = Vec::with_capacity(intervals.len());
        let mut ends = Vec::with_capacity(intervals.len());

        for interval in intervals {
            starts.push(interval.first);
            ends.push(interval.last);
        }

        starts.sort_unstable();
        ends.sort_unstable();

        Self { starts, ends }
    }

    fn query_count(&self, start: i32, end: i32) -> i64 {
        if end < start {
            return 0;
        }

        let started = self.starts.partition_point(|&value| value <= end);
        let ended_before = self.ends.partition_point(|&value| value < start);
        (started - ended_before) as i64
    }
}

pub fn merge_intervals(mut intervals: Vec<Interval<()>>) -> Vec<Interval<()>> {
    if intervals.is_empty() {
        return vec![];
    }

    intervals.sort_by(|a, b| a.first.cmp(&b.first));

    let mut merged = Vec::with_capacity(intervals.len());
    let mut current = intervals[0];

    for interval in intervals.into_iter().skip(1) {
        if interval.first <= current.last {
            current.last = current.last.max(interval.last);
        } else {
            merged.push(current);
            current = interval;
        }
    }
    merged.push(current);

    merged
}

pub fn build_coitree_from_batches(
    batches: Vec<RecordBatch>,
    columns: (&str, &str, &str),
    coverage: bool,
) -> Result<AHashMap<String, COITree<(), u32>>> {
    let mut nodes = IntervalHashMap::default();

    for batch in batches {
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(&batch, (columns.0, columns.1, columns.2))?;
        let start_resolved = start_arr.resolve()?;
        let end_resolved = end_arr.resolve()?;
        let starts = &*start_resolved;
        let ends = &*end_resolved;

        for i in 0..batch.num_rows() {
            let contig = contig_arr.value(i);
            let interval = Interval::new(starts[i], ends[i], ());

            if let Some(seqname_nodes) = nodes.get_mut(contig) {
                seqname_nodes.push(interval);
            } else {
                nodes.insert(contig.to_owned(), vec![interval]);
            }
        }
    }

    let mut trees = AHashMap::<String, COITree<(), u32>>::default();
    for (seqname, seqname_nodes) in nodes {
        if coverage {
            trees.insert(seqname, COITree::new(&merge_intervals(seqname_nodes)));
        } else {
            trees.insert(seqname, COITree::new(&seqname_nodes));
        }
    }
    Ok(trees)
}

pub fn build_count_index_from_batches(
    batches: Vec<RecordBatch>,
    columns: (&str, &str, &str),
) -> Result<AHashMap<String, CountOverlapIndex>> {
    let mut nodes = IntervalHashMap::default();

    for batch in batches {
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(&batch, (columns.0, columns.1, columns.2))?;
        let start_resolved = start_arr.resolve()?;
        let end_resolved = end_arr.resolve()?;
        let starts = &*start_resolved;
        let ends = &*end_resolved;

        for i in 0..batch.num_rows() {
            let contig = contig_arr.value(i);
            let interval = Interval::new(starts[i], ends[i], ());

            if let Some(seqname_nodes) = nodes.get_mut(contig) {
                seqname_nodes.push(interval);
            } else {
                nodes.insert(contig.to_owned(), vec![interval]);
            }
        }
    }

    Ok(nodes
        .into_iter()
        .map(|(seqname, intervals)| (seqname, CountOverlapIndex::new(intervals)))
        .collect())
}

/// Sum the covered bases shared by a query range and a tree of *merged* targets.
///
/// `start`/`end` are the bounds as handed to [`COITree::query`], which means they
/// arrive already narrowed by one base on each side when `half_open` is set --
/// see the `strict_filter` branch in [`get_stream`]. That narrowing exists purely
/// so coitrees' inclusive containment test implements half-open overlap; it must
/// be undone before any length is measured, which is what `query_start` and
/// `query_end` recover below.
///
/// `half_open` names the coordinate convention shared by the query and the tree:
/// `true` for 0-based half-open `[start, end)` (`FilterOp::Strict`), where the
/// length of a range is `end - start`; `false` for 1-based inclusive
/// `[start, end]` (`FilterOp::Weak`), where both endpoints are covered bases and
/// the length is therefore `end - start + 1`.
pub fn get_coverage(tree: &COITree<(), u32>, start: i32, end: i32, half_open: bool) -> i32 {
    let (query_start, query_end) = if half_open {
        (start - 1, end + 1)
    } else {
        (start, end)
    };
    // A 1-based range covers both of its endpoints; a half-open one covers only
    // its start.
    let endpoint = i32::from(!half_open);

    let mut coverage = 0;
    tree.query(start, end, |node| {
        let overlap = min(query_end, node.last) - max(query_start, node.first) + endpoint;
        coverage += max(0, overlap);
    });
    coverage
}

#[allow(clippy::too_many_arguments)]
pub fn get_stream(
    right_plan: Arc<dyn ExecutionPlan>,
    trees: Arc<AHashMap<String, COITree<(), u32>>>,
    new_schema: SchemaRef,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    coverage: bool,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let partition_stream = right_plan.execute(partition, context)?;
    let schema_for_closure = new_schema.clone();
    let strict_filter = filter_op == FilterOp::Strict;

    let iter = partition_stream.map(move |rb| match rb {
        Ok(rb) => {
            let (contig, pos_start, pos_end) =
                get_join_col_arrays(&rb, (&columns_2.0, &columns_2.1, &columns_2.2))?;
            let start_resolved = pos_start.resolve()?;
            let end_resolved = pos_end.resolve()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            let mut count_arr = Vec::with_capacity(rb.num_rows());
            let num_rows = rb.num_rows();
            let mut cached_contig: Option<&str> = None;
            let mut cached_tree: Option<&COITree<(), u32>> = None;
            for i in 0..num_rows {
                let contig = contig.value(i);
                let mut query_start = starts[i];
                let mut query_end = ends[i];
                if strict_filter {
                    query_start += 1;
                    query_end -= 1;
                }

                let tree = if cached_contig == Some(contig) {
                    cached_tree
                } else {
                    cached_contig = Some(contig);
                    cached_tree = trees.get(contig);
                    cached_tree
                };
                let count = match tree {
                    None => 0,
                    Some(tree) => {
                        if coverage {
                            get_coverage(tree, query_start, query_end, strict_filter)
                        } else {
                            tree.query_count(query_start, query_end) as i32
                        }
                    }
                };
                count_arr.push(count as i64);
            }
            let count_arr = Arc::new(Int64Array::from(count_arr));
            let mut columns = Vec::with_capacity(rb.num_columns() + 1);
            columns.extend_from_slice(rb.columns());
            columns.push(count_arr);
            RecordBatch::try_new(schema_for_closure.clone(), columns)
                .map_err(|e| datafusion::common::DataFusionError::ArrowError(Box::new(e), None))
        }
        Err(e) => Err(e),
    });

    let adapted_stream = RecordBatchStreamAdapter::new(new_schema, Box::pin(iter) as BoxStream<_>);
    Ok(Box::pin(adapted_stream))
}

#[allow(clippy::too_many_arguments)]
pub fn get_count_stream(
    right_plan: Arc<dyn ExecutionPlan>,
    indexes: Arc<AHashMap<String, CountOverlapIndex>>,
    new_schema: SchemaRef,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let partition_stream = right_plan.execute(partition, context)?;
    let schema_for_closure = new_schema.clone();
    let strict_filter = filter_op == FilterOp::Strict;

    let iter = partition_stream.map(move |rb| match rb {
        Ok(rb) => {
            let (contig, pos_start, pos_end) =
                get_join_col_arrays(&rb, (&columns_2.0, &columns_2.1, &columns_2.2))?;
            let start_resolved = pos_start.resolve()?;
            let end_resolved = pos_end.resolve()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            let mut count_arr = Vec::with_capacity(rb.num_rows());
            let num_rows = rb.num_rows();
            let mut cached_contig: Option<&str> = None;
            let mut cached_index: Option<&CountOverlapIndex> = None;
            for i in 0..num_rows {
                let contig = contig.value(i);
                let mut query_start = starts[i];
                let mut query_end = ends[i];
                if strict_filter {
                    query_start += 1;
                    query_end -= 1;
                }

                let index = if cached_contig == Some(contig) {
                    cached_index
                } else {
                    cached_contig = Some(contig);
                    cached_index = indexes.get(contig);
                    cached_index
                };
                let count = index.map_or(0, |index| index.query_count(query_start, query_end));
                count_arr.push(count);
            }
            let count_arr = Arc::new(Int64Array::from(count_arr));
            let mut columns = Vec::with_capacity(rb.num_columns() + 1);
            columns.extend_from_slice(rb.columns());
            columns.push(count_arr);
            RecordBatch::try_new(schema_for_closure.clone(), columns)
                .map_err(|e| datafusion::common::DataFusionError::ArrowError(Box::new(e), None))
        }
        Err(e) => Err(e),
    });

    let adapted_stream = RecordBatchStreamAdapter::new(new_schema, Box::pin(iter) as BoxStream<_>);
    Ok(Box::pin(adapted_stream))
}

#[cfg(test)]
mod tests {
    use super::*;
    use coitrees::Interval;

    #[test]
    fn test_merge_intervals_empty() {
        let result = merge_intervals(vec![]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_merge_intervals_non_overlapping() {
        let intervals = vec![
            Interval::new(1, 5, ()),
            Interval::new(10, 15, ()),
            Interval::new(20, 25, ()),
        ];
        let result = merge_intervals(intervals);
        assert_eq!(result.len(), 3);
        assert_eq!((result[0].first, result[0].last), (1, 5));
        assert_eq!((result[1].first, result[1].last), (10, 15));
        assert_eq!((result[2].first, result[2].last), (20, 25));
    }

    #[test]
    fn test_merge_intervals_overlapping() {
        let intervals = vec![
            Interval::new(1, 5, ()),
            Interval::new(3, 8, ()),
            Interval::new(10, 15, ()),
        ];
        let result = merge_intervals(intervals);
        assert_eq!(result.len(), 2);
        assert_eq!((result[0].first, result[0].last), (1, 8));
        assert_eq!((result[1].first, result[1].last), (10, 15));
    }

    #[test]
    fn test_merge_intervals_adjacent() {
        let intervals = vec![Interval::new(1, 5, ()), Interval::new(5, 10, ())];
        let result = merge_intervals(intervals);
        assert_eq!(result.len(), 1);
        assert_eq!((result[0].first, result[0].last), (1, 10));
    }

    #[test]
    fn test_merge_intervals_all_overlapping() {
        let intervals = vec![
            Interval::new(150, 250, ()),
            Interval::new(190, 300, ()),
            Interval::new(300, 501, ()),
            Interval::new(500, 700, ()),
        ];
        let result = merge_intervals(intervals);
        assert_eq!(result.len(), 1);
        assert_eq!((result[0].first, result[0].last), (150, 700));
    }

    #[test]
    fn test_merge_intervals_unsorted() {
        let intervals = vec![
            Interval::new(10, 15, ()),
            Interval::new(1, 5, ()),
            Interval::new(3, 8, ()),
        ];
        let result = merge_intervals(intervals);
        assert_eq!(result.len(), 2);
        assert_eq!((result[0].first, result[0].last), (1, 8));
        assert_eq!((result[1].first, result[1].last), (10, 15));
    }

    #[test]
    fn test_count_overlap_index_counts_inclusive_overlaps() {
        let index = CountOverlapIndex::new(vec![
            Interval::new(10, 20, ()),
            Interval::new(15, 25, ()),
            Interval::new(30, 40, ()),
        ]);

        assert_eq!(index.query_count(5, 9), 0);
        assert_eq!(index.query_count(5, 10), 1);
        assert_eq!(index.query_count(20, 30), 3);
        assert_eq!(index.query_count(26, 29), 0);
        assert_eq!(index.query_count(31, 35), 1);
        assert_eq!(index.query_count(10, 9), 0);
    }

    // Bounds are passed raw, i.e. in the already-narrowed form get_stream
    // hands to get_coverage on the 0-based path; the constants are unchanged
    // from before half_open was threaded through and guard that path.
    #[test]
    fn test_get_coverage_single_interval() {
        let intervals = vec![Interval::new(100, 200, ())];
        let tree = COITree::new(&intervals);

        // Query fully contained in interval
        assert_eq!(get_coverage(&tree, 120, 150, true), 32);
        // Query partially overlapping
        assert_eq!(get_coverage(&tree, 50, 120, true), 21);
        // Query with no overlap
        assert_eq!(get_coverage(&tree, 300, 400, true), 0);
    }

    // Bounds are passed raw, i.e. in the already-narrowed form get_stream
    // hands to get_coverage on the 0-based path; the constants are unchanged
    // from before half_open was threaded through and guard that path.
    #[test]
    fn test_get_coverage_multiple_merged_intervals() {
        // Simulate the chr1 merged reads from the CSV test data
        let merged = merge_intervals(vec![
            Interval::new(150, 250, ()),
            Interval::new(190, 300, ()),
            Interval::new(300, 501, ()),
            Interval::new(500, 700, ()),
        ]);
        assert_eq!(merged.len(), 1);
        assert_eq!((merged[0].first, merged[0].last), (150, 700));

        let tree = COITree::new(&merged);

        // These match the expected coverage values from test_coverage_csv
        assert_eq!(get_coverage(&tree, 100, 190, true), 41);
        assert_eq!(get_coverage(&tree, 200, 290, true), 92);
        assert_eq!(get_coverage(&tree, 400, 600, true), 202);
    }

    /// Mirror what [`get_stream`] does to a query row before calling
    /// [`get_coverage`]: under `FilterOp::Strict` (0-based, half-open) it narrows
    /// the row by one base on each side so that coitrees' inclusive containment
    /// test implements half-open overlap.
    fn coverage_of(tree: &COITree<(), u32>, query: (i32, i32), half_open: bool) -> i32 {
        let (start, end) = if half_open {
            (query.0 + 1, query.1 - 1)
        } else {
            query
        };
        get_coverage(tree, start, end, half_open)
    }

    fn tree_of(target: (i32, i32)) -> COITree<(), u32> {
        COITree::new(&[Interval::new(target.0, target.1, ())])
    }

    #[test]
    fn test_get_coverage_point_interval_one_based() {
        // 1-based inclusive: [15000, 15000] is a single covered base.
        let tree = tree_of((15000, 15000));
        assert_eq!(coverage_of(&tree, (10000, 20000), false), 1);
        assert_eq!(coverage_of(&tree, (15000, 15000), false), 1);
    }

    #[test]
    fn test_get_coverage_point_interval_zero_based() {
        // 0-based half-open: [15000, 15000) is empty and covers nothing.
        let tree = tree_of((15000, 15000));
        assert_eq!(coverage_of(&tree, (10000, 20000), true), 0);
        assert_eq!(coverage_of(&tree, (15000, 15000), true), 0);
    }

    #[test]
    fn test_get_coverage_zero_based_is_half_open_intersection() {
        // Query [10, 20) -- 10 covered bases.
        let query = (10, 20);
        for (target, expected) in [
            ((10, 20), 10), // identical to the query
            ((5, 25), 10),  // strict superset
            ((0, 100), 10), // much larger superset
            ((9, 21), 10),  // overhangs both ends
            ((12, 15), 3),  // contained
            ((5, 15), 5),   // overhangs left
            ((19, 30), 1),  // overhangs right
            ((0, 10), 0),   // bookended, no shared base
            ((20, 30), 0),  // bookended, no shared base
            ((0, 5), 0),    // disjoint
            ((25, 30), 0),  // disjoint
        ] {
            assert_eq!(
                coverage_of(&tree_of(target), query, true),
                expected,
                "0-based query {query:?} against target {target:?}"
            );
        }
    }

    #[test]
    fn test_get_coverage_one_based_is_inclusive_intersection() {
        // Query [10, 20] -- 11 covered bases, both endpoints included.
        let query = (10, 20);
        for (target, expected) in [
            ((10, 20), 11), // identical to the query
            ((5, 25), 11),  // strict superset
            ((0, 100), 11), // much larger superset
            ((9, 21), 11),  // overhangs both ends
            ((12, 15), 4),  // contained
            ((5, 15), 6),   // overhangs left
            ((19, 30), 2),  // overhangs right
            ((0, 10), 1),   // touches at base 10
            ((20, 30), 1),  // touches at base 20
            ((0, 9), 0),    // disjoint
            ((21, 30), 0),  // disjoint
        ] {
            assert_eq!(
                coverage_of(&tree_of(target), query, false),
                expected,
                "1-based query {query:?} against target {target:?}"
            );
        }
    }

    #[test]
    fn test_get_coverage_never_exceeds_query_length() {
        // No target, however large, can cover more bases than the query has.
        let query = (10, 20);
        for half_open in [true, false] {
            let query_len = if half_open {
                query.1 - query.0
            } else {
                query.1 - query.0 + 1
            };
            for start in 0..30 {
                for end in start..30 {
                    let coverage = coverage_of(&tree_of((start, end)), query, half_open);
                    assert!(
                        (0..=query_len).contains(&coverage),
                        "half_open={half_open} target=({start}, {end}) gave {coverage}, \
                         outside 0..={query_len}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_get_coverage_sums_disjoint_targets() {
        let tree = COITree::new(&merge_intervals(vec![
            Interval::new(12, 14, ()),
            Interval::new(16, 18, ()),
        ]));
        // 0-based: [12,14) and [16,18) contribute 2 + 2.
        assert_eq!(coverage_of(&tree, (10, 20), true), 4);
        // 1-based: [12,14] and [16,18] contribute 3 + 3.
        assert_eq!(coverage_of(&tree, (10, 20), false), 6);
    }
}
