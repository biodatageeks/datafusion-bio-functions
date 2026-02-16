use std::cmp::{max, min};
use std::sync::Arc;

use coitrees::{COITree, Interval, IntervalTree};
use datafusion::arrow::array::{Int64Array, RecordBatch};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::Partitioning;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::prelude::SessionContext;
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

    let mut merged = Vec::new();
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
    columns: (String, String, String),
    coverage: bool,
) -> FnvHashMap<String, COITree<(), u32>> {
    let mut nodes = IntervalHashMap::default();

    for batch in batches {
        let (contig_arr, start_arr, end_arr) = get_join_col_arrays(&batch, columns.clone());

        for i in 0..batch.num_rows() {
            let contig = contig_arr.value(i).to_string();
            let pos_start = start_arr.value(i);
            let pos_end = end_arr.value(i);
            nodes
                .entry(contig)
                .or_default()
                .push(Interval::new(pos_start, pos_end, ()));
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
    trees
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
pub async fn get_stream(
    session: Arc<SessionContext>,
    trees: Arc<FnvHashMap<String, COITree<(), u32>>>,
    right_table: String,
    new_schema: SchemaRef,
    _columns_1: (String, String, String),
    columns_2: (String, String, String),
    filter_op: FilterOp,
    coverage: bool,
    target_partitions: usize,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let right_table = session.table(right_table);
    let table_stream = right_table.await?;
    let plan = table_stream.create_physical_plan().await?;
    let repartition_stream =
        RepartitionExec::try_new(plan, Partitioning::RoundRobinBatch(target_partitions))?;

    let partition_stream = repartition_stream.execute(partition, context)?;
    let new_schema_out = new_schema.clone();

    let iter = partition_stream.map(move |rb| match rb {
        Ok(rb) => {
            let (contig, pos_start, pos_end) = get_join_col_arrays(&rb, columns_2.clone());
            let mut count_arr = Vec::with_capacity(rb.num_rows());
            let num_rows = rb.num_rows();
            for i in 0..num_rows {
                let contig = contig.value(i).to_string();
                let pos_start = pos_start.value(i);
                let pos_end = pos_end.value(i);
                let tree = trees.get(&contig);
                if tree.is_none() {
                    count_arr.push(0);
                    continue;
                }
                let count = if coverage {
                    if filter_op == FilterOp::Strict {
                        get_coverage(tree.unwrap(), pos_start + 1, pos_end - 1)
                    } else {
                        get_coverage(tree.unwrap(), pos_start, pos_end)
                    }
                } else if filter_op == FilterOp::Strict {
                    tree.unwrap().query_count(pos_start + 1, pos_end - 1) as i32
                } else {
                    tree.unwrap().query_count(pos_start, pos_end) as i32
                };
                count_arr.push(count as i64);
            }
            let count_arr = Arc::new(Int64Array::from(count_arr));
            let mut columns = rb.columns().to_vec();
            columns.push(count_arr);
            let new_rb = RecordBatch::try_new(new_schema.clone(), columns).unwrap();
            Ok(new_rb)
        }
        Err(e) => Err(e),
    });

    let adapted_stream =
        RecordBatchStreamAdapter::new(new_schema_out, Box::pin(iter) as BoxStream<_>);
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
