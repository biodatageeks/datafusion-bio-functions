use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use ahash::AHashMap;
use async_trait::async_trait;
use datafusion::arrow::array::{Int64Builder, RecordBatch, StringBuilder};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{MemTable, TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{Expr, SessionContext};

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;

pub struct ClusterProvider {
    session: Arc<SessionContext>,
    table: String,
    columns: (String, String, String),
    min_dist: i64,
    filter_op: FilterOp,
    schema: SchemaRef,
}

impl ClusterProvider {
    pub fn new(
        session: Arc<SessionContext>,
        table: String,
        columns: (String, String, String),
        min_dist: i64,
        filter_op: FilterOp,
    ) -> Self {
        let schema = Arc::new(Schema::new(vec![
            Arc::new(Field::new(&columns.0, DataType::Utf8, false)),
            Arc::new(Field::new(&columns.1, DataType::Int64, false)),
            Arc::new(Field::new(&columns.2, DataType::Int64, false)),
            Arc::new(Field::new("cluster", DataType::Int64, false)),
            Arc::new(Field::new("cluster_start", DataType::Int64, false)),
            Arc::new(Field::new("cluster_end", DataType::Int64, false)),
        ]));
        Self {
            session,
            table,
            columns,
            min_dist,
            filter_op,
            schema,
        }
    }
}

impl Debug for ClusterProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ClusterProvider {{ table: {}, min_dist: {} }}",
            self.table, self.min_dist
        )
    }
}

#[async_trait]
impl TableProvider for ClusterProvider {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn table_type(&self) -> TableType {
        TableType::Temporary
    }

    async fn scan(
        &self,
        state: &dyn Session,
        projection: Option<&Vec<usize>>,
        filters: &[Expr],
        limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let df = self.session.table(&self.table).await?;
        let batches = df.collect().await?;

        // Collect intervals per contig, preserving original indices for cluster assignment
        let mut groups: AHashMap<String, Vec<(i64, i64, usize)>> = AHashMap::default();
        let mut all_rows: Vec<(String, i64, i64)> = Vec::new();

        for batch in &batches {
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(batch, (&self.columns.0, &self.columns.1, &self.columns.2))?;
            let start_resolved = start_arr.resolve_i64()?;
            let end_resolved = end_arr.resolve_i64()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i).to_string();
                let global_idx = all_rows.len();
                groups
                    .entry(contig.clone())
                    .or_default()
                    .push((starts[i], ends[i], global_idx));
                all_rows.push((contig, starts[i], ends[i]));
            }
        }

        let n_rows = all_rows.len();
        // Per-row cluster assignment
        let mut row_cluster_id: Vec<i64> = vec![0; n_rows];
        let mut row_cluster_start: Vec<i64> = vec![0; n_rows];
        let mut row_cluster_end: Vec<i64> = vec![0; n_rows];

        let min_dist = self.min_dist;
        let strict = self.filter_op == FilterOp::Strict;

        // Sort contigs for deterministic output
        let mut contigs: Vec<String> = groups.keys().cloned().collect();
        contigs.sort();

        let mut cluster_id: i64 = 0;

        for contig in &contigs {
            let intervals = groups.get_mut(contig.as_str()).unwrap();
            // Sort by (start, end)
            intervals.sort_unstable_by_key(|&(s, e, _)| (s, e));

            if intervals.is_empty() {
                continue;
            }

            let mut cur_start = intervals[0].0;
            let mut cur_end = intervals[0].1;
            let mut cluster_members: Vec<usize> = vec![intervals[0].2];

            for &(s, e, idx) in &intervals[1..] {
                let merge_condition = if strict {
                    s < cur_end + min_dist
                } else {
                    s <= cur_end + min_dist
                };

                if merge_condition {
                    if e > cur_end {
                        cur_end = e;
                    }
                    cluster_members.push(idx);
                } else {
                    // Finalize current cluster
                    for &member_idx in &cluster_members {
                        row_cluster_id[member_idx] = cluster_id;
                        row_cluster_start[member_idx] = cur_start;
                        row_cluster_end[member_idx] = cur_end;
                    }
                    cluster_id += 1;
                    cur_start = s;
                    cur_end = e;
                    cluster_members.clear();
                    cluster_members.push(idx);
                }
            }

            // Finalize last cluster
            for &member_idx in &cluster_members {
                row_cluster_id[member_idx] = cluster_id;
                row_cluster_start[member_idx] = cur_start;
                row_cluster_end[member_idx] = cur_end;
            }
            cluster_id += 1;
        }

        // Build output in sorted order (contig, start, end) like merge does
        // We already have all_rows indexed, and groups are sorted by contig
        // Re-walk in sorted order to emit rows
        let mut contig_builder = StringBuilder::with_capacity(n_rows, n_rows * 4);
        let mut start_builder = Int64Builder::with_capacity(n_rows);
        let mut end_builder = Int64Builder::with_capacity(n_rows);
        let mut cluster_builder = Int64Builder::with_capacity(n_rows);
        let mut cluster_start_builder = Int64Builder::with_capacity(n_rows);
        let mut cluster_end_builder = Int64Builder::with_capacity(n_rows);

        for contig in &contigs {
            let intervals = groups.get(contig.as_str()).unwrap();
            // intervals are already sorted by (start, end)
            for &(s, e, idx) in intervals {
                contig_builder.append_value(contig);
                start_builder.append_value(s);
                end_builder.append_value(e);
                cluster_builder.append_value(row_cluster_id[idx]);
                cluster_start_builder.append_value(row_cluster_start[idx]);
                cluster_end_builder.append_value(row_cluster_end[idx]);
            }
        }

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(contig_builder.finish()),
                Arc::new(start_builder.finish()),
                Arc::new(end_builder.finish()),
                Arc::new(cluster_builder.finish()),
                Arc::new(cluster_start_builder.finish()),
                Arc::new(cluster_end_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;

        let mem_table = MemTable::try_new(self.schema.clone(), vec![vec![batch]])?;
        mem_table.scan(state, projection, filters, limit).await
    }
}
