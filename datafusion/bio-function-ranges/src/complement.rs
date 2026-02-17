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

pub struct ComplementProvider {
    session: Arc<SessionContext>,
    table: String,
    view_table: Option<String>,
    columns: (String, String, String),
    view_columns: (String, String, String),
    filter_op: FilterOp,
    schema: SchemaRef,
}

impl ComplementProvider {
    pub fn new(
        session: Arc<SessionContext>,
        table: String,
        view_table: Option<String>,
        columns: (String, String, String),
        view_columns: (String, String, String),
        filter_op: FilterOp,
    ) -> Self {
        let schema = Arc::new(Schema::new(vec![
            Arc::new(Field::new(&columns.0, DataType::Utf8, false)),
            Arc::new(Field::new(&columns.1, DataType::Int64, false)),
            Arc::new(Field::new(&columns.2, DataType::Int64, false)),
        ]));
        Self {
            session,
            table,
            view_table,
            columns,
            view_columns,
            filter_op,
            schema,
        }
    }
}

impl Debug for ComplementProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ComplementProvider {{ table: {}, view: {:?} }}",
            self.table, self.view_table
        )
    }
}

/// Merge overlapping intervals in-place (intervals must be sorted by start).
/// Returns merged (start, end) pairs.
fn merge_sorted_intervals(intervals: &[(i64, i64)], strict: bool) -> Vec<(i64, i64)> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut merged = Vec::new();
    let mut cur_start = intervals[0].0;
    let mut cur_end = intervals[0].1;

    for &(s, e) in &intervals[1..] {
        let merge_condition = if strict { s < cur_end } else { s <= cur_end };
        if merge_condition {
            if e > cur_end {
                cur_end = e;
            }
        } else {
            merged.push((cur_start, cur_end));
            cur_start = s;
            cur_end = e;
        }
    }
    merged.push((cur_start, cur_end));
    merged
}

#[async_trait]
impl TableProvider for ComplementProvider {
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
        let mut groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
        for batch in &batches {
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(batch, (&self.columns.0, &self.columns.1, &self.columns.2))?;
            let start_resolved = start_arr.resolve_i64()?;
            let end_resolved = end_arr.resolve_i64()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i).to_string();
                groups.entry(contig).or_default().push((starts[i], ends[i]));
            }
        }

        // Sort each contig's intervals by start
        for intervals in groups.values_mut() {
            intervals.sort_unstable();
        }

        // Load view boundaries
        let strict = self.filter_op == FilterOp::Strict;
        let view_bounds: AHashMap<String, Vec<(i64, i64)>> =
            if let Some(view_name) = &self.view_table {
                let vdf = self.session.table(view_name).await?;
                let vbatches = vdf.collect().await?;
                let mut vgroups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
                for batch in &vbatches {
                    let (contig_arr, start_arr, end_arr) = get_join_col_arrays(
                        batch,
                        (
                            &self.view_columns.0,
                            &self.view_columns.1,
                            &self.view_columns.2,
                        ),
                    )?;
                    let start_resolved = start_arr.resolve_i64()?;
                    let end_resolved = end_arr.resolve_i64()?;
                    let starts = &*start_resolved;
                    let ends = &*end_resolved;
                    for i in 0..batch.num_rows() {
                        let contig = contig_arr.value(i).to_string();
                        vgroups
                            .entry(contig)
                            .or_default()
                            .push((starts[i], ends[i]));
                    }
                }
                // Sort view intervals by start for each contig
                for intervals in vgroups.values_mut() {
                    intervals.sort_unstable();
                }
                vgroups
            } else {
                // No view: use 0..i64::MAX for each contig present in input
                groups
                    .keys()
                    .map(|c| (c.clone(), vec![(0i64, i64::MAX)]))
                    .collect()
            };

        // Merge overlapping intervals per contig
        let merged: AHashMap<String, Vec<(i64, i64)>> = groups
            .iter()
            .map(|(contig, intervals)| (contig.clone(), merge_sorted_intervals(intervals, strict)))
            .collect();

        // Sort contigs for deterministic output
        let mut contigs: Vec<String> = view_bounds.keys().cloned().collect();
        contigs.sort();

        // Estimate output size
        let estimated_rows: usize = contigs.len() * 2;
        let mut contig_builder = StringBuilder::with_capacity(estimated_rows, estimated_rows * 4);
        let mut start_builder = Int64Builder::with_capacity(estimated_rows);
        let mut end_builder = Int64Builder::with_capacity(estimated_rows);

        for contig in &contigs {
            let view_intervals = &view_bounds[contig];
            let empty_vec = Vec::new();
            let merged_intervals = merged.get(contig.as_str()).unwrap_or(&empty_vec);

            for &(view_start, view_end) in view_intervals {
                // Find merged intervals that fall within this view
                let mut cursor = view_start;

                for &(ms, me) in merged_intervals {
                    // Skip intervals entirely before the view
                    if me <= view_start {
                        continue;
                    }
                    // Stop when intervals are entirely past the view
                    if ms >= view_end {
                        break;
                    }

                    // Clamp to view boundaries
                    let interval_start = ms.max(view_start);
                    let interval_end = me.min(view_end);

                    // Emit gap before this interval
                    if interval_start > cursor {
                        contig_builder.append_value(contig);
                        start_builder.append_value(cursor);
                        end_builder.append_value(interval_start);
                    }
                    cursor = interval_end;
                }

                // Emit trailing gap
                if cursor < view_end {
                    contig_builder.append_value(contig);
                    start_builder.append_value(cursor);
                    end_builder.append_value(view_end);
                }
            }
        }

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(contig_builder.finish()),
                Arc::new(start_builder.finish()),
                Arc::new(end_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;

        let mem_table = MemTable::try_new(self.schema.clone(), vec![vec![batch]])?;
        mem_table.scan(state, projection, filters, limit).await
    }
}
