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

struct ContigGroup {
    contig: Arc<str>,
    intervals: Vec<(i64, i64)>,
}

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

        // Group intervals by contig while owning each contig string only once.
        let mut groups: Vec<ContigGroup> = Vec::new();
        let mut group_idx_by_contig: AHashMap<Arc<str>, usize> = AHashMap::default();
        let mut n_rows = 0usize;
        let mut contig_bytes = 0usize;

        for batch in &batches {
            n_rows += batch.num_rows();
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(batch, (&self.columns.0, &self.columns.1, &self.columns.2))?;
            let start_resolved = start_arr.resolve_i64()?;
            let end_resolved = end_arr.resolve_i64()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i);
                contig_bytes += contig.len();

                if let Some(&group_idx) = group_idx_by_contig.get(contig) {
                    groups[group_idx].intervals.push((starts[i], ends[i]));
                } else {
                    let contig_key: Arc<str> = Arc::from(contig);
                    let group_idx = groups.len();
                    group_idx_by_contig.insert(Arc::clone(&contig_key), group_idx);
                    groups.push(ContigGroup {
                        contig: contig_key,
                        intervals: vec![(starts[i], ends[i])],
                    });
                }
            }
        }

        let min_dist = self.min_dist;
        let strict = self.filter_op == FilterOp::Strict;

        // Sort contigs for deterministic output
        let mut group_order: Vec<usize> = (0..groups.len()).collect();
        group_order.sort_unstable_by(|&a, &b| groups[a].contig.cmp(&groups[b].contig));

        let mut contig_builder = StringBuilder::with_capacity(n_rows, contig_bytes);
        let mut start_builder = Int64Builder::with_capacity(n_rows);
        let mut end_builder = Int64Builder::with_capacity(n_rows);
        let mut cluster_builder = Int64Builder::with_capacity(n_rows);
        let mut cluster_start_builder = Int64Builder::with_capacity(n_rows);
        let mut cluster_end_builder = Int64Builder::with_capacity(n_rows);

        let mut cluster_id: i64 = 0;

        for group_idx in group_order {
            let group = &mut groups[group_idx];
            let intervals = &mut group.intervals;
            intervals.sort_unstable();

            if intervals.is_empty() {
                continue;
            }

            let mut cluster_start_idx = 0usize;
            let mut cur_start = intervals[0].0;
            let mut cur_end = intervals[0].1;

            for i in 1..intervals.len() {
                let (s, e) = intervals[i];
                let boundary = cur_end.saturating_add(min_dist);
                let merge_condition = if strict { s < boundary } else { s <= boundary };

                if merge_condition {
                    if e > cur_end {
                        cur_end = e;
                    }
                } else {
                    for &(cs, ce) in &intervals[cluster_start_idx..i] {
                        contig_builder.append_value(group.contig.as_ref());
                        start_builder.append_value(cs);
                        end_builder.append_value(ce);
                        cluster_builder.append_value(cluster_id);
                        cluster_start_builder.append_value(cur_start);
                        cluster_end_builder.append_value(cur_end);
                    }
                    cluster_id += 1;
                    cluster_start_idx = i;
                    cur_start = s;
                    cur_end = e;
                }
            }

            for &(cs, ce) in &intervals[cluster_start_idx..] {
                contig_builder.append_value(group.contig.as_ref());
                start_builder.append_value(cs);
                end_builder.append_value(ce);
                cluster_builder.append_value(cluster_id);
                cluster_start_builder.append_value(cur_start);
                cluster_end_builder.append_value(cur_end);
            }
            cluster_id += 1;
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
