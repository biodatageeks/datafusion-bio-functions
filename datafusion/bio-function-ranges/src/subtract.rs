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
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use futures::{Stream, StreamExt, ready};

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;
use crate::grouped_stream::{DEFAULT_BATCH_SIZE, StreamCollector};

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
        let left = self.left.execute(partition, Arc::clone(&context))?;
        let right = self.right.execute(partition, context)?;
        Ok(Box::pin(SubtractStream {
            schema: self.schema.clone(),
            left_collector: StreamCollector::new(left, Arc::clone(&self.left_columns)),
            right: Some(right),
            right_columns: Arc::clone(&self.right_columns),
            strict: self.strict,
            phase: SubtractPhase::CollectRight,
            right_groups: BTreeMap::new(),
            left_groups: Vec::new(),
            group_idx: 0,
            interval_idx: 0,
            contig_builder: StringBuilder::new(),
            start_builder: Int64Builder::new(),
            end_builder: Int64Builder::new(),
            pending_rows: 0,
        }))
    }
}

enum SubtractPhase {
    CollectRight,
    CollectLeft,
    Emit,
    Done,
}

struct SubtractStream {
    schema: SchemaRef,
    left_collector: StreamCollector,
    right: Option<SendableRecordBatchStream>,
    right_columns: Arc<(String, String, String)>,
    strict: bool,
    phase: SubtractPhase,
    right_groups: BTreeMap<String, Vec<(i64, i64)>>,
    left_groups: Vec<(String, Vec<(i64, i64)>)>,
    group_idx: usize,
    interval_idx: usize,
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
                                        let contig = contig_arr.value(i);
                                        if let Some(vec) = this.right_groups.get_mut(contig) {
                                            vec.push((starts[i], ends[i]));
                                        } else {
                                            this.right_groups.insert(
                                                contig.to_string(),
                                                vec![(starts[i], ends[i])],
                                            );
                                        }
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
                    this.phase = SubtractPhase::CollectLeft;
                }
                SubtractPhase::CollectLeft => match ready!(this.left_collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.left_groups = this.left_collector.take_groups().into_iter().collect();
                        this.group_idx = 0;
                        this.interval_idx = 0;
                        this.phase = SubtractPhase::Emit;
                    }
                    Ok(false) => unreachable!(),
                    Err(e) => {
                        this.phase = SubtractPhase::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                },
                SubtractPhase::Emit => {
                    let empty = Vec::new();

                    while this.group_idx < this.left_groups.len() {
                        let (ref contig, ref left_intervals) = this.left_groups[this.group_idx];
                        let right_intervals = this.right_groups.get(contig).unwrap_or(&empty);

                        let mut right_cursor: usize = 0;

                        while this.interval_idx < left_intervals.len() {
                            let (ls, le) = left_intervals[this.interval_idx];
                            this.interval_idx += 1;

                            // Advance cursor past right intervals that end before left start
                            while right_cursor < right_intervals.len() {
                                let skip = if this.strict {
                                    right_intervals[right_cursor].1 <= ls
                                } else {
                                    right_intervals[right_cursor].1 < ls
                                };
                                if skip {
                                    right_cursor += 1;
                                } else {
                                    break;
                                }
                            }

                            let mut cursor = ls;
                            let mut j = right_cursor;
                            while j < right_intervals.len() {
                                let (rs, re) = right_intervals[j];
                                let no_overlap = if this.strict { rs >= le } else { rs > le };
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

                            if this.pending_rows >= DEFAULT_BATCH_SIZE {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;

                        if this.pending_rows >= DEFAULT_BATCH_SIZE {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

                    // All groups done
                    this.phase = SubtractPhase::Done;
                    this.left_groups.clear();
                    if this.pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders()));
                    }
                    return Poll::Ready(None);
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
