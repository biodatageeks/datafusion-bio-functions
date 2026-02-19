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

        Ok(Arc::new(ClusterExec {
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
struct ClusterExec {
    schema: SchemaRef,
    input: Arc<dyn ExecutionPlan>,
    columns: Arc<(String, String, String)>,
    min_dist: i64,
    strict: bool,
    cache: PlanProperties,
}

impl DisplayAs for ClusterExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ClusterExec: min_dist={}, strict={}",
            self.min_dist, self.strict
        )
    }
}

impl ExecutionPlan for ClusterExec {
    fn name(&self) -> &str {
        "ClusterExec"
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
                "ClusterExec expects exactly one child plan".to_string(),
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
        Ok(Box::pin(ClusterStream {
            schema: self.schema.clone(),
            collector: StreamCollector::new(input, Arc::clone(&self.columns)),
            min_dist: self.min_dist,
            strict: self.strict,
            phase: ClusterPhase::Collecting,
            groups: Vec::new(),
            group_idx: 0,
            interval_idx: 0,
            cluster_id: 0,
            cluster_start: 0,
            cluster_end: 0,
            pending_intervals: Vec::new(),
            contig_builder: StringBuilder::new(),
            start_builder: Int64Builder::new(),
            end_builder: Int64Builder::new(),
            cluster_builder: Int64Builder::new(),
            cluster_start_builder: Int64Builder::new(),
            cluster_end_builder: Int64Builder::new(),
            pending_rows: 0,
        }))
    }
}

enum ClusterPhase {
    Collecting,
    Emitting,
    Done,
}

struct ClusterStream {
    schema: SchemaRef,
    collector: StreamCollector,
    min_dist: i64,
    strict: bool,
    phase: ClusterPhase,
    /// Sorted groups: Vec of (contig, intervals) in BTreeMap order.
    groups: Vec<(String, Vec<(i64, i64)>)>,
    group_idx: usize,
    interval_idx: usize,
    cluster_id: i64,
    cluster_start: i64,
    cluster_end: i64,
    pending_intervals: Vec<(i64, i64)>,
    contig_builder: StringBuilder,
    start_builder: Int64Builder,
    end_builder: Int64Builder,
    cluster_builder: Int64Builder,
    cluster_start_builder: Int64Builder,
    cluster_end_builder: Int64Builder,
    pending_rows: usize,
}

impl ClusterStream {
    fn flush_pending_cluster(&mut self, contig: &str) {
        let cluster_id = self.cluster_id;
        let cluster_start = self.cluster_start;
        let cluster_end = self.cluster_end;
        for &(s, e) in &self.pending_intervals {
            self.contig_builder.append_value(contig);
            self.start_builder.append_value(s);
            self.end_builder.append_value(e);
            self.cluster_builder.append_value(cluster_id);
            self.cluster_start_builder.append_value(cluster_start);
            self.cluster_end_builder.append_value(cluster_end);
        }
        self.pending_rows += self.pending_intervals.len();
        self.pending_intervals.clear();
        self.cluster_id += 1;
    }

    fn flush_builders(&mut self) -> Result<RecordBatch> {
        self.pending_rows = 0;
        RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(self.contig_builder.finish()),
                Arc::new(self.start_builder.finish()),
                Arc::new(self.end_builder.finish()),
                Arc::new(self.cluster_builder.finish()),
                Arc::new(self.cluster_start_builder.finish()),
                Arc::new(self.cluster_end_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for ClusterStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                ClusterPhase::Collecting => match ready!(this.collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.groups = this.collector.take_groups();
                        this.group_idx = 0;
                        this.interval_idx = 0;
                        this.phase = ClusterPhase::Emitting;
                    }
                    Ok(false) => unreachable!(),
                    Err(e) => {
                        this.phase = ClusterPhase::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                },
                ClusterPhase::Emitting => {
                    while this.group_idx < this.groups.len() {
                        // Clone contig name up front to avoid borrow conflicts
                        let contig = this.groups[this.group_idx].0.clone();
                        let interval_count = this.groups[this.group_idx].1.len();

                        while this.interval_idx < interval_count {
                            let (s, e) = this.groups[this.group_idx].1[this.interval_idx];
                            this.interval_idx += 1;

                            if this.pending_intervals.is_empty() {
                                this.cluster_start = s;
                                this.cluster_end = e;
                                this.pending_intervals.push((s, e));
                            } else {
                                let boundary = this.cluster_end.saturating_add(this.min_dist);
                                let merge_condition = if this.strict {
                                    s < boundary
                                } else {
                                    s <= boundary
                                };
                                if merge_condition {
                                    if e > this.cluster_end {
                                        this.cluster_end = e;
                                    }
                                    this.pending_intervals.push((s, e));
                                } else {
                                    this.flush_pending_cluster(&contig);
                                    this.cluster_start = s;
                                    this.cluster_end = e;
                                    this.pending_intervals.push((s, e));
                                }
                            }

                            if this.pending_rows >= DEFAULT_BATCH_SIZE {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        // Contig done — emit final cluster
                        if !this.pending_intervals.is_empty() {
                            this.flush_pending_cluster(&contig);
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;

                        if this.pending_rows >= DEFAULT_BATCH_SIZE {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

                    // All groups done
                    this.phase = ClusterPhase::Done;
                    this.groups.clear();
                    if this.pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders()));
                    }
                    return Poll::Ready(None);
                }
                ClusterPhase::Done => {
                    return Poll::Ready(None);
                }
            }
        }
    }
}

impl RecordBatchStream for ClusterStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}
