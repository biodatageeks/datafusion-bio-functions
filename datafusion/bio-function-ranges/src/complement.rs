use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use ahash::AHashMap;
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

        let view_plan = if let Some(view_name) = &self.view_table {
            let view_df = self.session.table(view_name).await?.select_columns(&[
                &self.view_columns.0,
                &self.view_columns.1,
                &self.view_columns.2,
            ])?;
            let plan = view_df.create_physical_plan().await?;
            let view_partitions = plan.output_partitioning().partition_count();
            let plan: Arc<dyn ExecutionPlan> = if view_partitions > 1 || target_partitions > 1 {
                Arc::new(RepartitionExec::try_new(
                    plan,
                    Partitioning::Hash(
                        vec![Arc::new(Column::new(self.view_columns.0.as_str(), 0))],
                        target_partitions.max(1),
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
                EmissionType::Incremental,
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
        let input = self.input.execute(partition, Arc::clone(&context))?;
        let view_collector = self
            .view
            .as_ref()
            .map(|v| -> Result<StreamCollector> {
                let stream = v.execute(partition, Arc::clone(&context))?;
                Ok(StreamCollector::new(stream, Arc::clone(&self.view_columns)))
            })
            .transpose()?;
        Ok(Box::pin(ComplementStream {
            schema: self.schema.clone(),
            collector: StreamCollector::new(input, Arc::clone(&self.columns)),
            view_collector,
            strict: self.strict,
            phase: ComplementPhase::CollectView,
            view_bounds: AHashMap::new(),
            input_groups: Vec::new(),
            group_idx: 0,
            trailing_iter_idx: 0,
            seen_contigs: Vec::new(),
            contig_builder: StringBuilder::new(),
            start_builder: Int64Builder::new(),
            end_builder: Int64Builder::new(),
            pending_rows: 0,
        }))
    }
}

enum ComplementPhase {
    CollectView,
    CollectInput,
    Emit,
    EmitTrailingGaps,
    Done,
}

struct ComplementStream {
    schema: SchemaRef,
    collector: StreamCollector,
    view_collector: Option<StreamCollector>,
    strict: bool,
    phase: ComplementPhase,
    view_bounds: AHashMap<String, Vec<(i64, i64)>>,
    input_groups: Vec<(String, Vec<(i64, i64)>)>,
    group_idx: usize,
    trailing_iter_idx: usize,
    seen_contigs: Vec<String>,
    contig_builder: StringBuilder,
    start_builder: Int64Builder,
    end_builder: Int64Builder,
    pending_rows: usize,
}

impl ComplementStream {
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

    /// Merge overlapping/touching intervals in a sorted slice, respecting `strict` mode.
    fn merge_intervals(intervals: &[(i64, i64)], strict: bool) -> Vec<(i64, i64)> {
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

    /// Compute complement of merged intervals against view intervals for a contig.
    fn emit_contig_complement(
        contig: &str,
        merged_intervals: &[(i64, i64)],
        view_intervals: &[(i64, i64)],
        contig_builder: &mut StringBuilder,
        start_builder: &mut Int64Builder,
        end_builder: &mut Int64Builder,
    ) -> usize {
        let mut rows = 0;
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
                    rows += 1;
                }
                cursor = interval_end;
            }
            if cursor < view_end {
                contig_builder.append_value(contig);
                start_builder.append_value(cursor);
                end_builder.append_value(view_end);
                rows += 1;
            }
        }
        rows
    }
}

impl Stream for ComplementStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                ComplementPhase::CollectView => {
                    if let Some(ref mut vc) = this.view_collector {
                        match ready!(vc.poll_collect(cx)) {
                            Ok(true) => {}
                            Ok(false) => unreachable!(),
                            Err(e) => {
                                this.phase = ComplementPhase::Done;
                                return Poll::Ready(Some(Err(e)));
                            }
                        }
                        this.view_bounds = this.view_collector.take().unwrap().take_groups_as_map();
                    }
                    this.phase = ComplementPhase::CollectInput;
                }
                ComplementPhase::CollectInput => match ready!(this.collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.input_groups = this.collector.take_groups();
                        this.group_idx = 0;
                        this.phase = ComplementPhase::Emit;
                    }
                    Ok(false) => unreachable!(),
                    Err(e) => {
                        this.phase = ComplementPhase::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                },
                ComplementPhase::Emit => {
                    while this.group_idx < this.input_groups.len() {
                        let (ref contig, ref intervals) = this.input_groups[this.group_idx];

                        // For contigs with no explicit view, add implicit (0, i64::MAX)
                        if !this.view_bounds.contains_key(contig) {
                            this.view_bounds.insert(contig.clone(), vec![(0, i64::MAX)]);
                        }

                        let merged = Self::merge_intervals(intervals, this.strict);

                        let view_intervals = this
                            .view_bounds
                            .get(contig)
                            .map(|v| v.as_slice())
                            .unwrap_or(&[]);

                        this.pending_rows += Self::emit_contig_complement(
                            contig,
                            &merged,
                            view_intervals,
                            &mut this.contig_builder,
                            &mut this.start_builder,
                            &mut this.end_builder,
                        );

                        this.seen_contigs.push(contig.clone());
                        this.group_idx += 1;

                        if this.pending_rows >= DEFAULT_BATCH_SIZE {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

                    this.phase = ComplementPhase::EmitTrailingGaps;
                    // Fall through
                }
                ComplementPhase::EmitTrailingGaps => {
                    let mut view_contigs: Vec<String> = this.view_bounds.keys().cloned().collect();
                    view_contigs.sort();

                    while this.trailing_iter_idx < view_contigs.len() {
                        let contig = &view_contigs[this.trailing_iter_idx];
                        this.trailing_iter_idx += 1;

                        if this.seen_contigs.contains(contig) {
                            continue;
                        }

                        let view_intervals = this.view_bounds.get(contig).unwrap();
                        for &(view_start, view_end) in view_intervals {
                            this.contig_builder.append_value(contig);
                            this.start_builder.append_value(view_start);
                            this.end_builder.append_value(view_end);
                            this.pending_rows += 1;
                        }

                        if this.pending_rows >= DEFAULT_BATCH_SIZE {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

                    this.phase = ComplementPhase::Done;
                    this.input_groups.clear();
                    if this.pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders()));
                    }
                    return Poll::Ready(None);
                }
                ComplementPhase::Done => {
                    return Poll::Ready(None);
                }
            }
        }
    }
}

impl RecordBatchStream for ComplementStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}
