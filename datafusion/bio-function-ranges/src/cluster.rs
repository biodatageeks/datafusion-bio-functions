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

impl ClusterExec {
    fn sort_ordering(&self) -> LexOrdering {
        let schema = self.input.schema();
        [
            PhysicalSortExpr::new_default(Arc::new(Column::new(
                &self.columns.0,
                schema.index_of(&self.columns.0).unwrap(),
            ))),
            PhysicalSortExpr::new_default(Arc::new(Column::new(
                &self.columns.1,
                schema.index_of(&self.columns.1).unwrap(),
            ))),
            PhysicalSortExpr::new_default(Arc::new(Column::new(
                &self.columns.2,
                schema.index_of(&self.columns.2).unwrap(),
            ))),
        ]
        .into()
    }
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
        let batch_size = context.session_config().batch_size();
        let sort_exec = SortExec::new(self.sort_ordering(), Arc::clone(&self.input))
            .with_preserve_partitioning(true);
        let input = sort_exec.execute(partition, context)?;
        Ok(Box::pin(ClusterStream {
            schema: self.schema.clone(),
            input,
            columns: Arc::clone(&self.columns),
            min_dist: self.min_dist,
            strict: self.strict,
            batch_size,
            current_contig: None,
            cluster_id: 0,
            cluster_start: 0,
            cluster_end: 0,
            pending_intervals: Vec::new(),
            done: false,
        }))
    }
}

struct ClusterStream {
    schema: SchemaRef,
    input: SendableRecordBatchStream,
    columns: Arc<(String, String, String)>,
    min_dist: i64,
    strict: bool,
    batch_size: usize,
    current_contig: Option<String>,
    cluster_id: i64,
    cluster_start: i64,
    cluster_end: i64,
    pending_intervals: Vec<(i64, i64)>,
    done: bool,
}

impl ClusterStream {
    #[allow(clippy::too_many_arguments)]
    fn emit_cluster(
        &self,
        contig_builder: &mut StringBuilder,
        start_builder: &mut Int64Builder,
        end_builder: &mut Int64Builder,
        cluster_builder: &mut Int64Builder,
        cluster_start_builder: &mut Int64Builder,
        cluster_end_builder: &mut Int64Builder,
        contig: &str,
        pending: &[(i64, i64)],
        cluster_id: i64,
        cluster_start: i64,
        cluster_end: i64,
    ) {
        for &(s, e) in pending {
            contig_builder.append_value(contig);
            start_builder.append_value(s);
            end_builder.append_value(e);
            cluster_builder.append_value(cluster_id);
            cluster_start_builder.append_value(cluster_start);
            cluster_end_builder.append_value(cluster_end);
        }
    }

    fn flush_builders(
        &self,
        contig_builder: &mut StringBuilder,
        start_builder: &mut Int64Builder,
        end_builder: &mut Int64Builder,
        cluster_builder: &mut Int64Builder,
        cluster_start_builder: &mut Int64Builder,
        cluster_end_builder: &mut Int64Builder,
    ) -> Result<RecordBatch> {
        RecordBatch::try_new(
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
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for ClusterStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        if this.done {
            return Poll::Ready(None);
        }

        let mut contig_builder = StringBuilder::new();
        let mut start_builder = Int64Builder::new();
        let mut end_builder = Int64Builder::new();
        let mut cluster_builder = Int64Builder::new();
        let mut cluster_start_builder = Int64Builder::new();
        let mut cluster_end_builder = Int64Builder::new();
        let mut pending_rows = 0usize;

        loop {
            let batch_opt = ready!(this.input.poll_next_unpin(cx));

            match batch_opt {
                Some(Ok(batch)) => {
                    if batch.num_rows() == 0 {
                        continue;
                    }
                    let (contig_arr, start_arr, end_arr) = match get_join_col_arrays(
                        &batch,
                        (&this.columns.0, &this.columns.1, &this.columns.2),
                    ) {
                        Ok(v) => v,
                        Err(e) => {
                            this.done = true;
                            return Poll::Ready(Some(Err(e)));
                        }
                    };
                    let start_resolved = match start_arr.resolve_i64() {
                        Ok(v) => v,
                        Err(e) => {
                            this.done = true;
                            return Poll::Ready(Some(Err(e)));
                        }
                    };
                    let end_resolved = match end_arr.resolve_i64() {
                        Ok(v) => v,
                        Err(e) => {
                            this.done = true;
                            return Poll::Ready(Some(Err(e)));
                        }
                    };
                    let starts = &*start_resolved;
                    let ends = &*end_resolved;

                    for i in 0..batch.num_rows() {
                        let contig = contig_arr.value(i);
                        let s = starts[i];
                        let e = ends[i];

                        let same_contig = this
                            .current_contig
                            .as_ref()
                            .is_some_and(|c| c.as_str() == contig);

                        if same_contig {
                            let boundary = this.cluster_end.saturating_add(this.min_dist);
                            let merge_condition = if this.strict {
                                s < boundary
                            } else {
                                s <= boundary
                            };
                            if merge_condition {
                                // Extend cluster
                                if e > this.cluster_end {
                                    this.cluster_end = e;
                                }
                                this.pending_intervals.push((s, e));
                            } else {
                                // Emit current cluster
                                let pending = std::mem::take(&mut this.pending_intervals);
                                this.emit_cluster(
                                    &mut contig_builder,
                                    &mut start_builder,
                                    &mut end_builder,
                                    &mut cluster_builder,
                                    &mut cluster_start_builder,
                                    &mut cluster_end_builder,
                                    this.current_contig.as_ref().unwrap(),
                                    &pending,
                                    this.cluster_id,
                                    this.cluster_start,
                                    this.cluster_end,
                                );
                                pending_rows += pending.len();
                                this.cluster_id += 1;
                                this.cluster_start = s;
                                this.cluster_end = e;
                                this.pending_intervals.push((s, e));
                            }
                        } else {
                            // Contig change: emit current cluster if any
                            if this.current_contig.is_some() {
                                let pending = std::mem::take(&mut this.pending_intervals);
                                this.emit_cluster(
                                    &mut contig_builder,
                                    &mut start_builder,
                                    &mut end_builder,
                                    &mut cluster_builder,
                                    &mut cluster_start_builder,
                                    &mut cluster_end_builder,
                                    this.current_contig.as_ref().unwrap(),
                                    &pending,
                                    this.cluster_id,
                                    this.cluster_start,
                                    this.cluster_end,
                                );
                                pending_rows += pending.len();
                                this.cluster_id += 1;
                            }
                            this.current_contig = Some(contig.to_string());
                            this.cluster_start = s;
                            this.cluster_end = e;
                            this.pending_intervals.push((s, e));
                        }

                        if pending_rows >= this.batch_size {
                            return Poll::Ready(Some(this.flush_builders(
                                &mut contig_builder,
                                &mut start_builder,
                                &mut end_builder,
                                &mut cluster_builder,
                                &mut cluster_start_builder,
                                &mut cluster_end_builder,
                            )));
                        }
                    }
                }
                Some(Err(e)) => {
                    this.done = true;
                    return Poll::Ready(Some(Err(e)));
                }
                None => {
                    // Input exhausted — emit final cluster
                    this.done = true;
                    if this.current_contig.is_some() {
                        let pending = std::mem::take(&mut this.pending_intervals);
                        this.emit_cluster(
                            &mut contig_builder,
                            &mut start_builder,
                            &mut end_builder,
                            &mut cluster_builder,
                            &mut cluster_start_builder,
                            &mut cluster_end_builder,
                            this.current_contig.as_ref().unwrap(),
                            &pending,
                            this.cluster_id,
                            this.cluster_start,
                            this.cluster_end,
                        );
                        pending_rows += pending.len();
                    }
                    if pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders(
                            &mut contig_builder,
                            &mut start_builder,
                            &mut end_builder,
                            &mut cluster_builder,
                            &mut cluster_start_builder,
                            &mut cluster_end_builder,
                        )));
                    }
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
