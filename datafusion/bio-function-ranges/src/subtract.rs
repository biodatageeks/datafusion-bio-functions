use std::any::Any;
use std::collections::BTreeMap;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use async_trait::async_trait;
use datafusion::arrow::array::{Int64Builder, RecordBatch, StringBuilder};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::expressions::Column;
use datafusion::physical_expr::{EquivalenceProperties, LexOrdering, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::expressions::PhysicalSortExpr;
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::sorts::sort::SortExec;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use futures::{Stream, StreamExt, ready};

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
        let left_partitions = left_plan.output_partitioning().partition_count();
        let left_plan: Arc<dyn ExecutionPlan> = if left_partitions > 1 || target_partitions > 1 {
            Arc::new(RepartitionExec::try_new(
                left_plan,
                Partitioning::Hash(
                    vec![Arc::new(Column::new(self.left_columns.0.as_str(), 0))],
                    target_partitions.max(1),
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
        let right_partitions = right_plan.output_partitioning().partition_count();
        let right_plan: Arc<dyn ExecutionPlan> = if right_partitions > 1 || target_partitions > 1 {
            Arc::new(RepartitionExec::try_new(
                right_plan,
                Partitioning::Hash(
                    vec![Arc::new(Column::new(self.right_columns.0.as_str(), 0))],
                    target_partitions.max(1),
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
                EmissionType::Incremental,
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

impl SubtractExec {
    fn left_sort_ordering(&self) -> LexOrdering {
        let schema = self.left.schema();
        [
            PhysicalSortExpr::new_default(Arc::new(Column::new(
                &self.left_columns.0,
                schema.index_of(&self.left_columns.0).unwrap(),
            ))),
            PhysicalSortExpr::new_default(Arc::new(Column::new(
                &self.left_columns.1,
                schema.index_of(&self.left_columns.1).unwrap(),
            ))),
            PhysicalSortExpr::new_default(Arc::new(Column::new(
                &self.left_columns.2,
                schema.index_of(&self.left_columns.2).unwrap(),
            ))),
        ]
        .into()
    }
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
                EmissionType::Incremental,
                Boundedness::Bounded,
            ),
        }))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let batch_size = context.session_config().batch_size();
        // Sort left at execution time to avoid optimizer stripping SortExec
        let left_sort = SortExec::new(self.left_sort_ordering(), Arc::clone(&self.left))
            .with_preserve_partitioning(true);
        let left = left_sort.execute(partition, Arc::clone(&context))?;
        let right = self.right.execute(partition, context)?;
        Ok(Box::pin(SubtractStream {
            schema: self.schema.clone(),
            left,
            right: Some(right),
            left_columns: Arc::clone(&self.left_columns),
            right_columns: Arc::clone(&self.right_columns),
            strict: self.strict,
            batch_size,
            phase: SubtractPhase::CollectRight,
            right_groups: BTreeMap::new(),
            right_cursors: BTreeMap::new(),
            contig_builder: StringBuilder::new(),
            start_builder: Int64Builder::new(),
            end_builder: Int64Builder::new(),
            pending_rows: 0,
        }))
    }
}

enum SubtractPhase {
    CollectRight,
    StreamLeft,
    Done,
}

struct SubtractStream {
    schema: SchemaRef,
    left: SendableRecordBatchStream,
    right: Option<SendableRecordBatchStream>,
    left_columns: Arc<(String, String, String)>,
    right_columns: Arc<(String, String, String)>,
    strict: bool,
    batch_size: usize,
    phase: SubtractPhase,
    right_groups: BTreeMap<String, Vec<(i64, i64)>>,
    right_cursors: BTreeMap<String, usize>,
    contig_builder: StringBuilder,
    start_builder: Int64Builder,
    end_builder: Int64Builder,
    pending_rows: usize,
}

impl SubtractStream {
    fn flush_builders(&mut self) -> Result<RecordBatch> {
        self.pending_rows = 0;
        RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(self.contig_builder.finish()),
                Arc::new(self.start_builder.finish()),
                Arc::new(self.end_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for SubtractStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                SubtractPhase::CollectRight => {
                    if let Some(ref mut right_stream) = this.right {
                        loop {
                            let batch_opt = ready!(right_stream.poll_next_unpin(cx));
                            match batch_opt {
                                Some(Ok(batch)) => {
                                    if batch.num_rows() == 0 {
                                        continue;
                                    }
                                    let right_cols = &this.right_columns;
                                    let (contig_arr, start_arr, end_arr) = match get_join_col_arrays(
                                        &batch,
                                        (&right_cols.0, &right_cols.1, &right_cols.2),
                                    ) {
                                        Ok(v) => v,
                                        Err(e) => {
                                            this.phase = SubtractPhase::Done;
                                            return Poll::Ready(Some(Err(e)));
                                        }
                                    };
                                    let start_resolved = match start_arr.resolve_i64() {
                                        Ok(v) => v,
                                        Err(e) => {
                                            this.phase = SubtractPhase::Done;
                                            return Poll::Ready(Some(Err(e)));
                                        }
                                    };
                                    let end_resolved = match end_arr.resolve_i64() {
                                        Ok(v) => v,
                                        Err(e) => {
                                            this.phase = SubtractPhase::Done;
                                            return Poll::Ready(Some(Err(e)));
                                        }
                                    };
                                    let starts = &*start_resolved;
                                    let ends = &*end_resolved;
                                    for i in 0..batch.num_rows() {
                                        this.right_groups
                                            .entry(contig_arr.value(i).to_string())
                                            .or_default()
                                            .push((starts[i], ends[i]));
                                    }
                                }
                                Some(Err(e)) => {
                                    this.phase = SubtractPhase::Done;
                                    return Poll::Ready(Some(Err(e)));
                                }
                                None => break,
                            }
                        }
                    }
                    this.right = None;
                    for intervals in this.right_groups.values_mut() {
                        intervals.sort_unstable();
                    }
                    this.phase = SubtractPhase::StreamLeft;
                }
                SubtractPhase::StreamLeft => {
                    loop {
                        let batch_opt = ready!(this.left.poll_next_unpin(cx));

                        match batch_opt {
                            Some(Ok(batch)) => {
                                if batch.num_rows() == 0 {
                                    continue;
                                }
                                let (contig_arr, start_arr, end_arr) = match get_join_col_arrays(
                                    &batch,
                                    (
                                        &this.left_columns.0,
                                        &this.left_columns.1,
                                        &this.left_columns.2,
                                    ),
                                ) {
                                    Ok(v) => v,
                                    Err(e) => {
                                        this.phase = SubtractPhase::Done;
                                        return Poll::Ready(Some(Err(e)));
                                    }
                                };
                                let start_resolved = match start_arr.resolve_i64() {
                                    Ok(v) => v,
                                    Err(e) => {
                                        this.phase = SubtractPhase::Done;
                                        return Poll::Ready(Some(Err(e)));
                                    }
                                };
                                let end_resolved = match end_arr.resolve_i64() {
                                    Ok(v) => v,
                                    Err(e) => {
                                        this.phase = SubtractPhase::Done;
                                        return Poll::Ready(Some(Err(e)));
                                    }
                                };
                                let starts = &*start_resolved;
                                let ends = &*end_resolved;

                                let empty = Vec::new();
                                for i in 0..batch.num_rows() {
                                    let contig = contig_arr.value(i);
                                    let ls = starts[i];
                                    let le = ends[i];

                                    let right_intervals =
                                        this.right_groups.get(contig).unwrap_or(&empty);
                                    let right_idx =
                                        this.right_cursors.entry(contig.to_string()).or_insert(0);

                                    // Advance cursor past right intervals that end before
                                    // left start
                                    while *right_idx < right_intervals.len() {
                                        let skip = if this.strict {
                                            right_intervals[*right_idx].1 <= ls
                                        } else {
                                            right_intervals[*right_idx].1 < ls
                                        };
                                        if skip {
                                            *right_idx += 1;
                                        } else {
                                            break;
                                        }
                                    }

                                    let mut cursor = ls;
                                    let mut j = *right_idx;
                                    while j < right_intervals.len() {
                                        let (rs, re) = right_intervals[j];
                                        let no_overlap =
                                            if this.strict { rs >= le } else { rs > le };
                                        if no_overlap {
                                            break;
                                        }

                                        if rs > cursor {
                                            this.contig_builder.append_value(contig);
                                            this.start_builder.append_value(cursor);
                                            this.end_builder.append_value(rs);
                                            this.pending_rows += 1;
                                        }
                                        if re > cursor {
                                            cursor = re;
                                        }
                                        j += 1;
                                    }

                                    if cursor < le {
                                        this.contig_builder.append_value(contig);
                                        this.start_builder.append_value(cursor);
                                        this.end_builder.append_value(le);
                                        this.pending_rows += 1;
                                    }
                                }

                                // Flush after processing the entire batch
                                if this.pending_rows >= this.batch_size {
                                    return Poll::Ready(Some(this.flush_builders()));
                                }
                            }
                            Some(Err(e)) => {
                                this.phase = SubtractPhase::Done;
                                return Poll::Ready(Some(Err(e)));
                            }
                            None => {
                                this.phase = SubtractPhase::Done;
                                if this.pending_rows > 0 {
                                    return Poll::Ready(Some(this.flush_builders()));
                                }
                                return Poll::Ready(None);
                            }
                        }
                    }
                }
                SubtractPhase::Done => {
                    return Poll::Ready(None);
                }
            }
        }
    }
}

impl RecordBatchStream for SubtractStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}
