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

pub struct MergeProvider {
    session: Arc<SessionContext>,
    table: String,
    columns: (String, String, String),
    min_dist: i64,
    filter_op: FilterOp,
    schema: SchemaRef,
}

impl MergeProvider {
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
            Arc::new(Field::new("n_intervals", DataType::Int64, false)),
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

impl Debug for MergeProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "MergeProvider {{ table: {}, min_dist: {} }}",
            self.table, self.min_dist
        )
    }
}

#[async_trait]
impl TableProvider for MergeProvider {
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

        // Group intervals by contig
        let mut groups: AHashMap<String, Vec<(i32, i32)>> = AHashMap::default();
        for batch in &batches {
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(batch, (&self.columns.0, &self.columns.1, &self.columns.2))?;
            let start_resolved = start_arr.resolve()?;
            let end_resolved = end_arr.resolve()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i).to_string();
                groups.entry(contig).or_default().push((starts[i], ends[i]));
            }
        }

        // Sort contigs for deterministic output
        let mut contigs: Vec<String> = groups.keys().cloned().collect();
        contigs.sort();

        // Estimate output size
        let estimated_rows: usize = groups.values().map(|v| v.len()).sum();

        let mut contig_builder = StringBuilder::with_capacity(estimated_rows, estimated_rows * 4);
        let mut start_builder = Int64Builder::with_capacity(estimated_rows);
        let mut end_builder = Int64Builder::with_capacity(estimated_rows);
        let mut count_builder = Int64Builder::with_capacity(estimated_rows);

        let min_dist = self.min_dist;
        let strict = self.filter_op == FilterOp::Strict;

        for contig in &contigs {
            let intervals = groups.get_mut(contig.as_str()).unwrap();
            intervals.sort_unstable();

            if intervals.is_empty() {
                continue;
            }

            let mut cur_start = intervals[0].0;
            let mut cur_end = intervals[0].1;
            let mut count: i64 = 1;

            for &(s, e) in &intervals[1..] {
                let merge_condition = if strict {
                    i64::from(s) < i64::from(cur_end) + min_dist
                } else {
                    i64::from(s) <= i64::from(cur_end) + min_dist
                };

                if merge_condition {
                    if e > cur_end {
                        cur_end = e;
                    }
                    count += 1;
                } else {
                    contig_builder.append_value(contig);
                    start_builder.append_value(i64::from(cur_start));
                    end_builder.append_value(i64::from(cur_end));
                    count_builder.append_value(count);
                    cur_start = s;
                    cur_end = e;
                    count = 1;
                }
            }

            // Emit final interval
            contig_builder.append_value(contig);
            start_builder.append_value(i64::from(cur_start));
            end_builder.append_value(i64::from(cur_end));
            count_builder.append_value(count);
        }

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(contig_builder.finish()),
                Arc::new(start_builder.finish()),
                Arc::new(end_builder.finish()),
                Arc::new(count_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;

        let mem_table = MemTable::try_new(self.schema.clone(), vec![vec![batch]])?;
        mem_table.scan(state, projection, filters, limit).await
    }
}
