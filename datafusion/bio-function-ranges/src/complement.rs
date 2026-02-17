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

        let input_df = self.session.table(&self.table).await?.select_columns(&[
            &self.columns.0,
            &self.columns.1,
            &self.columns.2,
        ])?;
        let input_plan = input_df.create_physical_plan().await?;
        let input_plan: Arc<dyn ExecutionPlan> = if target_partitions > 1 {
            Arc::new(RepartitionExec::try_new(
                input_plan,
                Partitioning::Hash(
                    vec![Arc::new(Column::new(self.columns.0.as_str(), 0))],
                    target_partitions,
                ),
            )?)
        } else {
            input_plan
        };

        let view_plan = if let Some(view_name) = &self.view_table {
            let view_df = self.session.table(view_name).await?.select_columns(&[
                &self.view_columns.0,
                &self.view_columns.1,
                &self.view_columns.2,
            ])?;
            let plan = view_df.create_physical_plan().await?;
            let plan: Arc<dyn ExecutionPlan> = if target_partitions > 1 {
                Arc::new(RepartitionExec::try_new(
                    plan,
                    Partitioning::Hash(
                        vec![Arc::new(Column::new(self.view_columns.0.as_str(), 0))],
                        target_partitions,
                    ),
                )?)
            } else {
                plan
            };
            Some(plan)
        } else {
            None
        };

        let output_partitions = input_plan.output_partitioning().partition_count();

        Ok(Arc::new(ComplementExec {
            schema: self.schema.clone(),
            input: input_plan,
            view: view_plan,
            columns: Arc::new(self.columns.clone()),
            view_columns: Arc::new(self.view_columns.clone()),
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
struct ComplementExec {
    schema: SchemaRef,
    input: Arc<dyn ExecutionPlan>,
    view: Option<Arc<dyn ExecutionPlan>>,
    columns: Arc<(String, String, String)>,
    view_columns: Arc<(String, String, String)>,
    strict: bool,
    cache: PlanProperties,
}

impl DisplayAs for ComplementExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "ComplementExec: strict={}", self.strict)
    }
}

impl ExecutionPlan for ComplementExec {
    fn name(&self) -> &str {
        "ComplementExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn properties(&self) -> &PlanProperties {
        &self.cache
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        let mut children = vec![&self.input];
        if let Some(view) = &self.view {
            children.push(view);
        }
        children
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let expected = if self.view.is_some() { 2 } else { 1 };
        if children.len() != expected {
            return Err(DataFusionError::Internal(format!(
                "ComplementExec expects exactly {expected} child plan(s)"
            )));
        }

        let new_view = if expected == 2 {
            Some(Arc::clone(&children[1]))
        } else {
            None
        };

        Ok(Arc::new(Self {
            schema: self.schema.clone(),
            input: Arc::clone(&children[0]),
            view: new_view,
            columns: Arc::clone(&self.columns),
            view_columns: Arc::clone(&self.view_columns),
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
        get_complement_stream(
            Arc::clone(&self.input),
            self.view.as_ref().map(Arc::clone),
            self.schema.clone(),
            Arc::clone(&self.columns),
            Arc::clone(&self.view_columns),
            self.strict,
            partition,
            context,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn get_complement_stream(
    input_plan: Arc<dyn ExecutionPlan>,
    view_plan: Option<Arc<dyn ExecutionPlan>>,
    schema: SchemaRef,
    columns: Arc<(String, String, String)>,
    view_columns: Arc<(String, String, String)>,
    strict: bool,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let schema_for_closure = schema.clone();
    let once_stream = once(async move {
        let input_stream = input_plan.execute(partition, Arc::clone(&context))?;
        let input_batches = collect(input_stream).await?;

        let view_batches = if let Some(view_plan) = view_plan {
            let view_stream = view_plan.execute(partition, context)?;
            Some(collect(view_stream).await?)
        } else {
            None
        };

        build_complement_batch(
            &input_batches,
            view_batches.as_deref(),
            &schema_for_closure,
            &columns,
            &view_columns,
            strict,
        )
    });

    let adapted_stream =
        RecordBatchStreamAdapter::new(schema, Box::pin(once_stream) as BoxStream<'_, _>);
    Ok(Box::pin(adapted_stream))
}

fn build_complement_batch(
    input_batches: &[RecordBatch],
    view_batches: Option<&[RecordBatch]>,
    schema: &SchemaRef,
    columns: &(String, String, String),
    view_columns: &(String, String, String),
    strict: bool,
) -> Result<RecordBatch> {
    let mut groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
    for batch in input_batches {
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(batch, (&columns.0, &columns.1, &columns.2))?;
        let start_resolved = start_arr.resolve_i64()?;
        let end_resolved = end_arr.resolve_i64()?;
        let starts = &*start_resolved;
        let ends = &*end_resolved;
        for i in 0..batch.num_rows() {
            groups
                .entry(contig_arr.value(i).to_string())
                .or_default()
                .push((starts[i], ends[i]));
        }
    }

    for intervals in groups.values_mut() {
        intervals.sort_unstable();
    }

    let view_bounds: AHashMap<String, Vec<(i64, i64)>> = if let Some(view_batches) = view_batches {
        let mut vgroups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
        for batch in view_batches {
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(batch, (&view_columns.0, &view_columns.1, &view_columns.2))?;
            let start_resolved = start_arr.resolve_i64()?;
            let end_resolved = end_arr.resolve_i64()?;
            let starts = &*start_resolved;
            let ends = &*end_resolved;
            for i in 0..batch.num_rows() {
                vgroups
                    .entry(contig_arr.value(i).to_string())
                    .or_default()
                    .push((starts[i], ends[i]));
            }
        }
        for intervals in vgroups.values_mut() {
            intervals.sort_unstable();
        }
        vgroups
    } else {
        groups
            .keys()
            .map(|contig| (contig.clone(), vec![(0_i64, i64::MAX)]))
            .collect()
    };

    let merged: AHashMap<String, Vec<(i64, i64)>> = groups
        .iter()
        .map(|(contig, intervals)| (contig.clone(), merge_sorted_intervals(intervals, strict)))
        .collect();

    let mut contigs: Vec<String> = view_bounds.keys().cloned().collect();
    contigs.sort_unstable();

    let estimated_rows = contigs.len() * 2;
    let mut contig_builder = StringBuilder::with_capacity(estimated_rows, estimated_rows * 8);
    let mut start_builder = Int64Builder::with_capacity(estimated_rows);
    let mut end_builder = Int64Builder::with_capacity(estimated_rows);

    for contig in &contigs {
        let view_intervals = &view_bounds[contig];
        let empty = Vec::new();
        let merged_intervals = merged.get(contig.as_str()).unwrap_or(&empty);

        for &(view_start, view_end) in view_intervals {
            let mut cursor = view_start;
            for &(ms, me) in merged_intervals {
                if me <= view_start {
                    continue;
                }
                if ms >= view_end {
                    break;
                }
                let interval_start = ms.max(view_start);
                let interval_end = me.min(view_end);
                if interval_start > cursor {
                    contig_builder.append_value(contig);
                    start_builder.append_value(cursor);
                    end_builder.append_value(interval_start);
                }
                cursor = interval_end;
            }
            if cursor < view_end {
                contig_builder.append_value(contig);
                start_builder.append_value(cursor);
                end_builder.append_value(view_end);
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
