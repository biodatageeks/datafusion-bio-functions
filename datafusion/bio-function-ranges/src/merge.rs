use std::any::Any;
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
use futures::{Stream, ready};

use crate::filter_op::FilterOp;
use crate::grouped_stream::{DEFAULT_BATCH_SIZE, StreamCollector};

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
        let input_partitions = input_plan.output_partitioning().partition_count();
        let input_plan: Arc<dyn ExecutionPlan> = if input_partitions > 1 || target_partitions > 1 {
            Arc::new(RepartitionExec::try_new(
                input_plan,
                Partitioning::Hash(
                    vec![Arc::new(Column::new(self.columns.0.as_str(), 0))],
                    target_partitions.max(1),
                ),
            )?)
        } else {
            input_plan
        };

        let output_partitions = input_plan.output_partitioning().partition_count();

        Ok(Arc::new(MergeExec {
            schema: self.schema.clone(),
            input: input_plan,
            columns: Arc::new(self.columns.clone()),
            min_dist: self.min_dist,
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
struct MergeExec {
    schema: SchemaRef,
    input: Arc<dyn ExecutionPlan>,
    columns: Arc<(String, String, String)>,
    min_dist: i64,
    strict: bool,
    cache: PlanProperties,
}

impl DisplayAs for MergeExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "MergeExec: min_dist={}, strict={}",
            self.min_dist, self.strict
        )
    }
}

impl ExecutionPlan for MergeExec {
    fn name(&self) -> &str {
        "MergeExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn properties(&self) -> &PlanProperties {
        &self.cache
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.input]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        if children.len() != 1 {
            return Err(DataFusionError::Internal(
                "MergeExec expects exactly one child plan".to_string(),
            ));
        }

        Ok(Arc::new(Self {
            schema: self.schema.clone(),
            input: Arc::clone(&children[0]),
            columns: Arc::clone(&self.columns),
            min_dist: self.min_dist,
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
        let input = self.input.execute(partition, context)?;
        Ok(Box::pin(MergeStream {
            schema: self.schema.clone(),
            collector: StreamCollector::new(input, Arc::clone(&self.columns)),
            min_dist: self.min_dist,
            strict: self.strict,
            phase: MergePhase::Collecting,
            groups: Vec::new(),
            group_idx: 0,
            interval_idx: 0,
            cur_start: 0,
            cur_end: 0,
            cur_count: 0,
            has_current: false,
            contig_builder: StringBuilder::new(),
            start_builder: Int64Builder::new(),
            end_builder: Int64Builder::new(),
            count_builder: Int64Builder::new(),
            pending_rows: 0,
        }))
    }
}

enum MergePhase {
    Collecting,
    Emitting,
    Done,
}

struct MergeStream {
    schema: SchemaRef,
    collector: StreamCollector,
    min_dist: i64,
    strict: bool,
    phase: MergePhase,
    /// Sorted groups: Vec of (contig, intervals) in BTreeMap order.
    groups: Vec<(String, Vec<(i64, i64)>)>,
    group_idx: usize,
    interval_idx: usize,
    cur_start: i64,
    cur_end: i64,
    cur_count: i64,
    has_current: bool,
    contig_builder: StringBuilder,
    start_builder: Int64Builder,
    end_builder: Int64Builder,
    count_builder: Int64Builder,
    pending_rows: usize,
}

impl MergeStream {
    fn flush_builders(&mut self) -> Result<RecordBatch> {
        self.pending_rows = 0;
        RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(self.contig_builder.finish()),
                Arc::new(self.start_builder.finish()),
                Arc::new(self.end_builder.finish()),
                Arc::new(self.count_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for MergeStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                MergePhase::Collecting => match ready!(this.collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.groups = this.collector.take_groups();
                        this.group_idx = 0;
                        this.interval_idx = 0;
                        this.has_current = false;
                        this.phase = MergePhase::Emitting;
                    }
                    Ok(false) => unreachable!(),
                    Err(e) => {
                        this.phase = MergePhase::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                },
                MergePhase::Emitting => {
                    while this.group_idx < this.groups.len() {
                        let (ref contig, ref intervals) = this.groups[this.group_idx];

                        while this.interval_idx < intervals.len() {
                            let (s, e) = intervals[this.interval_idx];
                            this.interval_idx += 1;

                            if this.has_current {
                                let boundary = this.cur_end.saturating_add(this.min_dist);
                                let merge_condition = if this.strict {
                                    s < boundary
                                } else {
                                    s <= boundary
                                };
                                if merge_condition {
                                    if e > this.cur_end {
                                        this.cur_end = e;
                                    }
                                    this.cur_count += 1;
                                } else {
                                    this.contig_builder.append_value(contig);
                                    this.start_builder.append_value(this.cur_start);
                                    this.end_builder.append_value(this.cur_end);
                                    this.count_builder.append_value(this.cur_count);
                                    this.pending_rows += 1;
                                    this.cur_start = s;
                                    this.cur_end = e;
                                    this.cur_count = 1;
                                }
                            } else {
                                this.cur_start = s;
                                this.cur_end = e;
                                this.cur_count = 1;
                                this.has_current = true;
                            }

                            if this.pending_rows >= DEFAULT_BATCH_SIZE {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        // Contig done — emit current interval
                        if this.has_current {
                            let contig = &this.groups[this.group_idx].0;
                            this.contig_builder.append_value(contig);
                            this.start_builder.append_value(this.cur_start);
                            this.end_builder.append_value(this.cur_end);
                            this.count_builder.append_value(this.cur_count);
                            this.pending_rows += 1;
                            this.has_current = false;
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;

                        if this.pending_rows >= DEFAULT_BATCH_SIZE {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

                    // All groups done
                    this.phase = MergePhase::Done;
                    this.groups.clear();
                    if this.pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders()));
                    }
                    return Poll::Ready(None);
                }
                MergePhase::Done => {
                    return Poll::Ready(None);
                }
            }
        }
    }
}

impl RecordBatchStream for MergeStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}
