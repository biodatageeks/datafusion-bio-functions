use std::cmp::{max, min};
use std::sync::Arc;

use coitrees::{COITree, Interval, IntervalTree};
use datafusion::arrow::array::{Int64Array, RecordBatch};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use fnv::FnvHashMap;
use futures::StreamExt;
use futures::stream::BoxStream;

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;

type IntervalHashMap = FnvHashMap<String, Vec<Interval<()>>>;

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
) -> Result<FnvHashMap<String, COITree<(), u32>>> {
    let mut nodes = IntervalHashMap::default();

    for batch in batches {
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(&batch, (columns.0, columns.1, columns.2))?;

        for i in 0..batch.num_rows() {
            let contig = contig_arr.value(i);
            let pos_start = start_arr.value(i)?;
            let pos_end = end_arr.value(i)?;
            let interval = Interval::new(pos_start, pos_end, ());

            if let Some(seqname_nodes) = nodes.get_mut(contig) {
                seqname_nodes.push(interval);
            } else {
                nodes.insert(contig.to_owned(), vec![interval]);
            }
        }
    }

    let mut trees = FnvHashMap::<String, COITree<(), u32>>::default();
    for (seqname, seqname_nodes) in nodes {
        if coverage {
            trees.insert(seqname, COITree::new(&merge_intervals(seqname_nodes)));
        } else {
            trees.insert(seqname, COITree::new(&seqname_nodes));
        }
    }
    Ok(trees)
}

pub fn get_coverage(tree: &COITree<(), u32>, start: i32, end: i32) -> i32 {
    let mut coverage = 0;
    tree.query(start, end, |node| {
        let overlap = max(1, min(end + 1, node.last) - max(start - 1, node.first));
        coverage += overlap;
    });
    coverage
}

#[allow(clippy::too_many_arguments)]
pub fn get_stream(
    right_plan: Arc<dyn ExecutionPlan>,
    trees: Arc<FnvHashMap<String, COITree<(), u32>>>,
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
            let mut count_arr = Vec::with_capacity(rb.num_rows());
            let num_rows = rb.num_rows();
            let mut cached_contig: Option<&str> = None;
            let mut cached_tree: Option<&COITree<(), u32>> = None;
            for i in 0..num_rows {
                let contig = contig.value(i);
                let mut query_start = pos_start.value(i)?;
                let mut query_end = pos_end.value(i)?;
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
                            get_coverage(tree, query_start, query_end)
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

    let adapted_stream =
        RecordBatchStreamAdapter::new(new_schema, Box::pin(iter) as BoxStream<_>);
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
    fn test_get_coverage_single_interval() {
        let intervals = vec![Interval::new(100, 200, ())];
        let tree = COITree::new(&intervals);

        // Query fully contained in interval
        assert_eq!(get_coverage(&tree, 120, 150), 32);
        // Query partially overlapping
        assert_eq!(get_coverage(&tree, 50, 120), 21);
        // Query with no overlap
        assert_eq!(get_coverage(&tree, 300, 400), 0);
    }

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
        assert_eq!(get_coverage(&tree, 100, 190), 41);
        assert_eq!(get_coverage(&tree, 200, 290), 92);
        assert_eq!(get_coverage(&tree, 400, 600), 202);
    }

    #[test]
    fn test_get_coverage_point_interval() {
        let intervals = vec![Interval::new(15000, 15000, ())];
        let tree = COITree::new(&intervals);

        // Query that contains the point
        assert_eq!(get_coverage(&tree, 10000, 20000), 1);
        // Query at exactly the point
        assert_eq!(get_coverage(&tree, 15000, 15000), 1);
    }
}
