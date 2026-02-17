use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use ahash::AHashMap;
use async_trait::async_trait;
use datafusion::arrow::array::{Int64Builder, RecordBatch, StringBuilder};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::expressions::Column;
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::common::collect;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use futures::stream::{BoxStream, once};

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
        _state: &dyn Session,
        _projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let target_partitions = self
            .session
            .state()
            .config()
            .options()
            .execution
            .target_partitions;

        let left_df = self
            .session
            .table(&self.left_table)
            .await?
            .select_columns(&[
                &self.left_columns.0,
                &self.left_columns.1,
                &self.left_columns.2,
            ])?;
        let left_plan = left_df.create_physical_plan().await?;
        let left_plan: Arc<dyn ExecutionPlan> = if target_partitions > 1 {
            Arc::new(RepartitionExec::try_new(
                left_plan,
                Partitioning::Hash(
                    vec![Arc::new(Column::new(self.left_columns.0.as_str(), 0))],
                    target_partitions,
                ),
            )?)
        } else {
            left_plan
        };

        let right_df = self
            .session
            .table(&self.right_table)
            .await?
            .select_columns(&[
                &self.right_columns.0,
                &self.right_columns.1,
                &self.right_columns.2,
            ])?;
        let right_plan = right_df.create_physical_plan().await?;
        let right_plan: Arc<dyn ExecutionPlan> = if target_partitions > 1 {
            Arc::new(RepartitionExec::try_new(
                right_plan,
                Partitioning::Hash(
                    vec![Arc::new(Column::new(self.right_columns.0.as_str(), 0))],
                    target_partitions,
                ),
            )?)
        } else {
            right_plan
        };

        let output_partitions = left_plan.output_partitioning().partition_count();

        Ok(Arc::new(SubtractExec {
            schema: self.schema.clone(),
            left: left_plan,
            right: right_plan,
            left_columns: Arc::new(self.left_columns.clone()),
            right_columns: Arc::new(self.right_columns.clone()),
            strict: self.filter_op == FilterOp::Strict,
            cache: PlanProperties::new(
                EquivalenceProperties::new(self.schema.clone()),
                Partitioning::UnknownPartitioning(output_partitions),
                EmissionType::Final,
                Boundedness::Bounded,
            ),
        }))
    }
}

#[derive(Debug)]
struct SubtractExec {
    schema: SchemaRef,
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    left_columns: Arc<(String, String, String)>,
    right_columns: Arc<(String, String, String)>,
    strict: bool,
    cache: PlanProperties,
}

impl DisplayAs for SubtractExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "SubtractExec: strict={}", self.strict)
    }
}

impl ExecutionPlan for SubtractExec {
    fn name(&self) -> &str {
        "SubtractExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn properties(&self) -> &PlanProperties {
        &self.cache
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.left, &self.right]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        if children.len() != 2 {
            return Err(DataFusionError::Internal(
                "SubtractExec expects exactly two child plans".to_string(),
            ));
        }

        Ok(Arc::new(Self {
            schema: self.schema.clone(),
            left: Arc::clone(&children[0]),
            right: Arc::clone(&children[1]),
            left_columns: Arc::clone(&self.left_columns),
            right_columns: Arc::clone(&self.right_columns),
            strict: self.strict,
            cache: PlanProperties::new(
                EquivalenceProperties::new(self.schema.clone()),
                Partitioning::UnknownPartitioning(
                    children[0].output_partitioning().partition_count(),
                ),
                EmissionType::Final,
                Boundedness::Bounded,
            ),
        }))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        get_subtract_stream(
            Arc::clone(&self.left),
            Arc::clone(&self.right),
            self.schema.clone(),
            Arc::clone(&self.left_columns),
            Arc::clone(&self.right_columns),
            self.strict,
            partition,
            context,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn get_subtract_stream(
    left_plan: Arc<dyn ExecutionPlan>,
    right_plan: Arc<dyn ExecutionPlan>,
    schema: SchemaRef,
    left_columns: Arc<(String, String, String)>,
    right_columns: Arc<(String, String, String)>,
    strict: bool,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let schema_for_closure = schema.clone();
    let once_stream = once(async move {
        let left_stream = left_plan.execute(partition, Arc::clone(&context))?;
        let left_batches = collect(left_stream).await?;

        let right_stream = right_plan.execute(partition, context)?;
        let right_batches = collect(right_stream).await?;

        build_subtract_batch(
            &left_batches,
            &right_batches,
            &schema_for_closure,
            &left_columns,
            &right_columns,
            strict,
        )
    });

    let adapted_stream =
        RecordBatchStreamAdapter::new(schema, Box::pin(once_stream) as BoxStream<'_, _>);
    Ok(Box::pin(adapted_stream))
}

fn build_subtract_batch(
    left_batches: &[RecordBatch],
    right_batches: &[RecordBatch],
    schema: &SchemaRef,
    left_columns: &(String, String, String),
    right_columns: &(String, String, String),
    strict: bool,
) -> Result<RecordBatch> {
    let mut left_groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
    for batch in left_batches {
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(batch, (&left_columns.0, &left_columns.1, &left_columns.2))?;
        let start_resolved = start_arr.resolve_i64()?;
        let end_resolved = end_arr.resolve_i64()?;
        let starts = &*start_resolved;
        let ends = &*end_resolved;
        for i in 0..batch.num_rows() {
            left_groups
                .entry(contig_arr.value(i).to_string())
                .or_default()
                .push((starts[i], ends[i]));
        }
    }

    let mut right_groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
    for batch in right_batches {
        let (contig_arr, start_arr, end_arr) = get_join_col_arrays(
            batch,
            (&right_columns.0, &right_columns.1, &right_columns.2),
        )?;
        let start_resolved = start_arr.resolve_i64()?;
        let end_resolved = end_arr.resolve_i64()?;
        let starts = &*start_resolved;
        let ends = &*end_resolved;
        for i in 0..batch.num_rows() {
            right_groups
                .entry(contig_arr.value(i).to_string())
                .or_default()
                .push((starts[i], ends[i]));
        }
    }

    for intervals in left_groups.values_mut() {
        intervals.sort_unstable();
    }
    for intervals in right_groups.values_mut() {
        intervals.sort_unstable();
    }

    let mut contigs: Vec<String> = left_groups.keys().cloned().collect();
    contigs.sort_unstable();

    let estimated_rows: usize = left_groups.values().map(Vec::len).sum();
    let mut contig_builder = StringBuilder::with_capacity(estimated_rows, estimated_rows * 8);
    let mut start_builder = Int64Builder::with_capacity(estimated_rows);
    let mut end_builder = Int64Builder::with_capacity(estimated_rows);

    let empty = Vec::new();
    for contig in &contigs {
        let left_intervals = left_groups.get(contig.as_str()).unwrap();
        let right_intervals = right_groups.get(contig.as_str()).unwrap_or(&empty);
        let mut right_idx = 0usize;

        for &(ls, le) in left_intervals {
            let mut cursor = ls;

            while right_idx < right_intervals.len() && right_intervals[right_idx].1 <= ls {
                if strict {
                    if right_intervals[right_idx].1 <= ls {
                        right_idx += 1;
                    } else {
                        break;
                    }
                } else if right_intervals[right_idx].1 < ls {
                    right_idx += 1;
                } else {
                    break;
                }
            }

            let mut j = right_idx;
            while j < right_intervals.len() {
                let (rs, re) = right_intervals[j];
                let no_overlap = if strict { rs >= le } else { rs > le };
                if no_overlap {
                    break;
                }

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

            if cursor < le {
                contig_builder.append_value(contig);
                start_builder.append_value(cursor);
                end_builder.append_value(le);
            }
        }
    }

    RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(contig_builder.finish()),
            Arc::new(start_builder.finish()),
            Arc::new(end_builder.finish()),
        ],
    )
    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}
