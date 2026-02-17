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

pub struct SubtractProvider {
    session: Arc<SessionContext>,
    left_table: String,
    right_table: String,
    left_columns: (String, String, String),
    right_columns: (String, String, String),
    filter_op: FilterOp,
    schema: SchemaRef,
}

impl SubtractProvider {
    pub fn new(
        session: Arc<SessionContext>,
        left_table: String,
        right_table: String,
        left_columns: (String, String, String),
        right_columns: (String, String, String),
        filter_op: FilterOp,
    ) -> Self {
        let schema = Arc::new(Schema::new(vec![
            Arc::new(Field::new(&left_columns.0, DataType::Utf8, false)),
            Arc::new(Field::new(&left_columns.1, DataType::Int64, false)),
            Arc::new(Field::new(&left_columns.2, DataType::Int64, false)),
        ]));
        Self {
            session,
            left_table,
            right_table,
            left_columns,
            right_columns,
            filter_op,
            schema,
        }
    }
}

impl Debug for SubtractProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SubtractProvider {{ left: {}, right: {} }}",
            self.left_table, self.right_table
        )
    }
}

#[async_trait]
impl TableProvider for SubtractProvider {
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
        // Load left table
        let left_df = self.session.table(&self.left_table).await?;
        let left_batches = left_df.collect().await?;

        let mut left_groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
        for batch in &left_batches {
            let (contig_arr, start_arr, end_arr) = get_join_col_arrays(
                batch,
                (
                    &self.left_columns.0,
                    &self.left_columns.1,
                    &self.left_columns.2,
                ),
            )?;
            let start_resolved = start_arr.resolve_i64()?;
            let end_resolved = end_arr.resolve_i64()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i).to_string();
                left_groups
                    .entry(contig)
                    .or_default()
                    .push((starts[i], ends[i]));
            }
        }

        // Load right table
        let right_df = self.session.table(&self.right_table).await?;
        let right_batches = right_df.collect().await?;

        let mut right_groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
        for batch in &right_batches {
            let (contig_arr, start_arr, end_arr) = get_join_col_arrays(
                batch,
                (
                    &self.right_columns.0,
                    &self.right_columns.1,
                    &self.right_columns.2,
                ),
            )?;
            let start_resolved = start_arr.resolve_i64()?;
            let end_resolved = end_arr.resolve_i64()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i).to_string();
                right_groups
                    .entry(contig)
                    .or_default()
                    .push((starts[i], ends[i]));
            }
        }

        // Sort intervals by start
        for intervals in left_groups.values_mut() {
            intervals.sort_unstable();
        }
        for intervals in right_groups.values_mut() {
            intervals.sort_unstable();
        }

        let strict = self.filter_op == FilterOp::Strict;

        // Sort contigs for deterministic output
        let mut contigs: Vec<String> = left_groups.keys().cloned().collect();
        contigs.sort();

        let estimated_rows: usize = left_groups.values().map(|v| v.len()).sum();
        let mut contig_builder = StringBuilder::with_capacity(estimated_rows, estimated_rows * 4);
        let mut start_builder = Int64Builder::with_capacity(estimated_rows);
        let mut end_builder = Int64Builder::with_capacity(estimated_rows);

        let empty_vec = Vec::new();

        for contig in &contigs {
            let left_intervals = left_groups.get(contig.as_str()).unwrap();
            let right_intervals = right_groups.get(contig.as_str()).unwrap_or(&empty_vec);

            let mut right_idx = 0;

            for &(ls, le) in left_intervals {
                let mut cursor = ls;

                // Advance right pointer to first interval that could overlap
                // An interval (rs, re) can overlap (ls, le) if re > ls (weak) or re > ls (strict uses same check here)
                while right_idx < right_intervals.len() && right_intervals[right_idx].1 <= ls {
                    // For strict: re must be > ls to overlap
                    // For weak: re must be > ls to overlap (re == ls means touching at boundary,
                    //   which in weak means overlap, but since we subtract we need re > ls)
                    // Actually for weak overlap: intervals overlap if rs <= le && re >= ls
                    // We need to be more careful. Let's just not advance past potential overlaps.
                    if strict {
                        // strict overlap: rs < le && re > ls
                        if right_intervals[right_idx].1 <= ls {
                            right_idx += 1;
                        } else {
                            break;
                        }
                    } else {
                        // weak overlap: rs <= le && re >= ls
                        if right_intervals[right_idx].1 < ls {
                            right_idx += 1;
                        } else {
                            break;
                        }
                    }
                }

                // Walk overlapping right intervals
                let mut j = right_idx;
                while j < right_intervals.len() {
                    let (rs, re) = right_intervals[j];

                    // Check if this right interval starts past the left end
                    let no_overlap = if strict { rs >= le } else { rs > le };
                    if no_overlap {
                        break;
                    }

                    // This right interval overlaps with [cursor, le)
                    if rs > cursor {
                        contig_builder.append_value(contig);
                        start_builder.append_value(cursor);
                        end_builder.append_value(rs);
                    }
                    if re > cursor {
                        cursor = re;
                    }
                    j += 1;
                }

                // Emit remaining portion
                if cursor < le {
                    contig_builder.append_value(contig);
                    start_builder.append_value(cursor);
                    end_builder.append_value(le);
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
