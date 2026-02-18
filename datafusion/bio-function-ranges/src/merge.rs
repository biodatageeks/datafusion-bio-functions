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

impl MergeExec {
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
        let batch_size = context.session_config().batch_size();
        // Sort at execution time to avoid optimizer stripping the SortExec
        let sort_exec = SortExec::new(self.sort_ordering(), Arc::clone(&self.input))
            .with_preserve_partitioning(true);
        let input = sort_exec.execute(partition, context)?;
        Ok(Box::pin(MergeStream {
            schema: self.schema.clone(),
            input,
            columns: Arc::clone(&self.columns),
            min_dist: self.min_dist,
            strict: self.strict,
            batch_size,
            current_contig: None,
            cur_start: 0,
            cur_end: 0,
            cur_count: 0,
            done: false,
            contig_builder: StringBuilder::new(),
            start_builder: Int64Builder::new(),
            end_builder: Int64Builder::new(),
            count_builder: Int64Builder::new(),
            pending_rows: 0,
        }))
    }
}

struct MergeStream {
    schema: SchemaRef,
    input: SendableRecordBatchStream,
    columns: Arc<(String, String, String)>,
    min_dist: i64,
    strict: bool,
    batch_size: usize,
    current_contig: Option<String>,
    cur_start: i64,
    cur_end: i64,
    cur_count: i64,
    done: bool,
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

        if this.done {
            return Poll::Ready(None);
        }

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
                                // Emit current interval
                                this.contig_builder
                                    .append_value(this.current_contig.as_ref().unwrap());
                                this.start_builder.append_value(this.cur_start);
                                this.end_builder.append_value(this.cur_end);
                                this.count_builder.append_value(this.cur_count);
                                this.pending_rows += 1;

                                this.cur_start = s;
                                this.cur_end = e;
                                this.cur_count = 1;
                            }
                        } else {
                            // Contig change: emit current interval if any
                            if this.current_contig.is_some() {
                                this.contig_builder
                                    .append_value(this.current_contig.as_ref().unwrap());
                                this.start_builder.append_value(this.cur_start);
                                this.end_builder.append_value(this.cur_end);
                                this.count_builder.append_value(this.cur_count);
                                this.pending_rows += 1;
                            }
                            this.current_contig = Some(contig.to_string());
                            this.cur_start = s;
                            this.cur_end = e;
                            this.cur_count = 1;
                        }
                    }

                    // Flush after processing the entire batch
                    if this.pending_rows >= this.batch_size {
                        return Poll::Ready(Some(this.flush_builders()));
                    }
                }
                Some(Err(e)) => {
                    this.done = true;
                    return Poll::Ready(Some(Err(e)));
                }
                None => {
                    // Input exhausted — emit final interval
                    this.done = true;
                    if this.current_contig.is_some() {
                        this.contig_builder
                            .append_value(this.current_contig.as_ref().unwrap());
                        this.start_builder.append_value(this.cur_start);
                        this.end_builder.append_value(this.cur_end);
                        this.count_builder.append_value(this.cur_count);
                        this.pending_rows += 1;
                    }
                    if this.pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders()));
                    }
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
