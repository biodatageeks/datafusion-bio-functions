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

impl ComplementExec {
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
        let batch_size = context.session_config().batch_size();
        let sort_exec = SortExec::new(self.sort_ordering(), Arc::clone(&self.input))
            .with_preserve_partitioning(true);
        let input = sort_exec.execute(partition, Arc::clone(&context))?;
        let view_stream = self
            .view
            .as_ref()
            .map(|v| v.execute(partition, context))
            .transpose()?;
        Ok(Box::pin(ComplementStream {
            schema: self.schema.clone(),
            input,
            view_stream,
            columns: Arc::clone(&self.columns),
            view_columns: Arc::clone(&self.view_columns),
            strict: self.strict,
            batch_size,
            phase: ComplementPhase::CollectView,
            view_bounds: BTreeMap::new(),
            current_contig: None,
            cur_merged_start: 0,
            cur_merged_end: 0,
            merged_intervals: Vec::new(),
            seen_contigs: Vec::new(),
            trailing_iter_idx: 0,
        }))
    }
}

enum ComplementPhase {
    CollectView,
    StreamInput,
    EmitTrailingGaps,
    Done,
}

struct ComplementStream {
    schema: SchemaRef,
    input: SendableRecordBatchStream,
    view_stream: Option<SendableRecordBatchStream>,
    columns: Arc<(String, String, String)>,
    view_columns: Arc<(String, String, String)>,
    strict: bool,
    batch_size: usize,
    phase: ComplementPhase,
    view_bounds: BTreeMap<String, Vec<(i64, i64)>>,
    current_contig: Option<String>,
    cur_merged_start: i64,
    cur_merged_end: i64,
    /// Collected merged (non-overlapping) intervals for the current contig.
    /// Typically small: one entry per contiguous coverage region.
    merged_intervals: Vec<(i64, i64)>,
    seen_contigs: Vec<String>,
    trailing_iter_idx: usize,
}

impl ComplementStream {
    fn flush_builders(
        &self,
        contig_builder: &mut StringBuilder,
        start_builder: &mut Int64Builder,
        end_builder: &mut Int64Builder,
    ) -> Result<RecordBatch> {
        RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(contig_builder.finish()),
                Arc::new(start_builder.finish()),
                Arc::new(end_builder.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }

    /// Compute complement of merged intervals against view intervals for a contig.
    /// This uses the same algorithm as the original `build_complement_batch`.
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

    /// Finalize the current contig: push last merged interval, compute complement, clear state.
    fn finalize_contig(
        &mut self,
        contig_builder: &mut StringBuilder,
        start_builder: &mut Int64Builder,
        end_builder: &mut Int64Builder,
    ) -> usize {
        if let Some(prev_contig) = self.current_contig.take() {
            // Push the last in-progress merged interval
            self.merged_intervals
                .push((self.cur_merged_start, self.cur_merged_end));

            let view_intervals = self
                .view_bounds
                .get(&prev_contig)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let rows = Self::emit_contig_complement(
                &prev_contig,
                &self.merged_intervals,
                view_intervals,
                contig_builder,
                start_builder,
                end_builder,
            );

            self.merged_intervals.clear();
            self.seen_contigs.push(prev_contig);
            rows
        } else {
            0
        }
    }
}

impl Stream for ComplementStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                ComplementPhase::CollectView => {
                    if let Some(ref mut view_stream) = this.view_stream {
                        loop {
                            let batch_opt = ready!(view_stream.poll_next_unpin(cx));
                            match batch_opt {
                                Some(Ok(batch)) => {
                                    if batch.num_rows() == 0 {
                                        continue;
                                    }
                                    let view_cols = &this.view_columns;
                                    let (contig_arr, start_arr, end_arr) = match get_join_col_arrays(
                                        &batch,
                                        (&view_cols.0, &view_cols.1, &view_cols.2),
                                    ) {
                                        Ok(v) => v,
                                        Err(e) => {
                                            this.phase = ComplementPhase::Done;
                                            return Poll::Ready(Some(Err(e)));
                                        }
                                    };
                                    let start_resolved = match start_arr.resolve_i64() {
                                        Ok(v) => v,
                                        Err(e) => {
                                            this.phase = ComplementPhase::Done;
                                            return Poll::Ready(Some(Err(e)));
                                        }
                                    };
                                    let end_resolved = match end_arr.resolve_i64() {
                                        Ok(v) => v,
                                        Err(e) => {
                                            this.phase = ComplementPhase::Done;
                                            return Poll::Ready(Some(Err(e)));
                                        }
                                    };
                                    let starts = &*start_resolved;
                                    let ends = &*end_resolved;
                                    for i in 0..batch.num_rows() {
                                        this.view_bounds
                                            .entry(contig_arr.value(i).to_string())
                                            .or_default()
                                            .push((starts[i], ends[i]));
                                    }
                                }
                                Some(Err(e)) => {
                                    this.phase = ComplementPhase::Done;
                                    return Poll::Ready(Some(Err(e)));
                                }
                                None => break,
                            }
                        }
                        for intervals in this.view_bounds.values_mut() {
                            intervals.sort_unstable();
                        }
                    }
                    this.view_stream = None;
                    this.phase = ComplementPhase::StreamInput;
                }
                ComplementPhase::StreamInput => {
                    let mut contig_builder = StringBuilder::new();
                    let mut start_builder = Int64Builder::new();
                    let mut end_builder = Int64Builder::new();
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
                                        this.phase = ComplementPhase::Done;
                                        return Poll::Ready(Some(Err(e)));
                                    }
                                };
                                let start_resolved = match start_arr.resolve_i64() {
                                    Ok(v) => v,
                                    Err(e) => {
                                        this.phase = ComplementPhase::Done;
                                        return Poll::Ready(Some(Err(e)));
                                    }
                                };
                                let end_resolved = match end_arr.resolve_i64() {
                                    Ok(v) => v,
                                    Err(e) => {
                                        this.phase = ComplementPhase::Done;
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
                                        let merge_condition = if this.strict {
                                            s < this.cur_merged_end
                                        } else {
                                            s <= this.cur_merged_end
                                        };
                                        if merge_condition {
                                            if e > this.cur_merged_end {
                                                this.cur_merged_end = e;
                                            }
                                        } else {
                                            // Finalize current merged interval, start new
                                            this.merged_intervals
                                                .push((this.cur_merged_start, this.cur_merged_end));
                                            this.cur_merged_start = s;
                                            this.cur_merged_end = e;
                                        }
                                    } else {
                                        // Contig change: compute complement for previous contig
                                        pending_rows += this.finalize_contig(
                                            &mut contig_builder,
                                            &mut start_builder,
                                            &mut end_builder,
                                        );

                                        // If no explicit view, add implicit (0, i64::MAX)
                                        if this.view_bounds.is_empty() {
                                            this.view_bounds
                                                .insert(contig.to_string(), vec![(0, i64::MAX)]);
                                        }

                                        this.current_contig = Some(contig.to_string());
                                        this.cur_merged_start = s;
                                        this.cur_merged_end = e;
                                    }

                                    if pending_rows >= this.batch_size {
                                        return Poll::Ready(Some(this.flush_builders(
                                            &mut contig_builder,
                                            &mut start_builder,
                                            &mut end_builder,
                                        )));
                                    }
                                }
                            }
                            Some(Err(e)) => {
                                this.phase = ComplementPhase::Done;
                                return Poll::Ready(Some(Err(e)));
                            }
                            None => {
                                // Input exhausted: finalize last contig
                                pending_rows += this.finalize_contig(
                                    &mut contig_builder,
                                    &mut start_builder,
                                    &mut end_builder,
                                );

                                this.phase = ComplementPhase::EmitTrailingGaps;

                                if pending_rows > 0 {
                                    return Poll::Ready(Some(this.flush_builders(
                                        &mut contig_builder,
                                        &mut start_builder,
                                        &mut end_builder,
                                    )));
                                }
                                break;
                            }
                        }
                    }
                    // Fall through to EmitTrailingGaps
                }
                ComplementPhase::EmitTrailingGaps => {
                    // Emit full view ranges for contigs never seen in input
                    let mut contig_builder = StringBuilder::new();
                    let mut start_builder = Int64Builder::new();
                    let mut end_builder = Int64Builder::new();
                    let mut pending_rows = 0usize;

                    let view_contigs: Vec<String> = this.view_bounds.keys().cloned().collect();

                    while this.trailing_iter_idx < view_contigs.len() {
                        let contig = &view_contigs[this.trailing_iter_idx];
                        this.trailing_iter_idx += 1;

                        if this.seen_contigs.contains(contig) {
                            continue;
                        }

                        let view_intervals = this.view_bounds.get(contig).unwrap();
                        for &(view_start, view_end) in view_intervals {
                            contig_builder.append_value(contig);
                            start_builder.append_value(view_start);
                            end_builder.append_value(view_end);
                            pending_rows += 1;
                        }

                        if pending_rows >= this.batch_size {
                            return Poll::Ready(Some(this.flush_builders(
                                &mut contig_builder,
                                &mut start_builder,
                                &mut end_builder,
                            )));
                        }
                    }

                    this.phase = ComplementPhase::Done;
                    if pending_rows > 0 {
                        return Poll::Ready(Some(this.flush_builders(
                            &mut contig_builder,
                            &mut start_builder,
                            &mut end_builder,
                        )));
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
