use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use ahash::AHashMap;
use async_trait::async_trait;
use datafusion::arrow::array::{Int64Array, Int64Builder, RecordBatch, StringBuilder, UInt32Array};
use datafusion::arrow::compute::take;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use futures::{Stream, ready};

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;
use crate::grouped_stream::{FullBatchCollector, IndexedGroups, StreamCollector};

pub struct SubtractProvider {
    session: Arc<SessionContext>,
    left_table: String,
    right_table: String,
    left_columns: (String, String, String),
    right_columns: (String, String, String),
    filter_op: FilterOp,
    schema: SchemaRef,
    has_extra_cols: bool,
}

impl SubtractProvider {
    pub fn new(
        session: Arc<SessionContext>,
        left_table: String,
        right_table: String,
        left_columns: (String, String, String),
        right_columns: (String, String, String),
        filter_op: FilterOp,
        left_schema: Schema,
    ) -> Self {
        let has_extra_cols = left_schema.fields().len() > 3;

        let fields: Vec<Arc<Field>> = if has_extra_cols {
            // Preserve all left fields as-is (take kernel returns same types);
            // start/end columns will be replaced with Int64 fragment values.
            let mut flds: Vec<Arc<Field>> = Vec::with_capacity(left_schema.fields().len());
            for f in left_schema.fields() {
                if f.name() == &left_columns.1 || f.name() == &left_columns.2 {
                    flds.push(Arc::new(Field::new(
                        f.name(),
                        DataType::Int64,
                        f.is_nullable(),
                    )));
                } else {
                    flds.push(Arc::clone(f));
                }
            }
            flds
        } else {
            vec![
                Arc::new(Field::new(&left_columns.0, DataType::Utf8, false)),
                Arc::new(Field::new(&left_columns.1, DataType::Int64, false)),
                Arc::new(Field::new(&left_columns.2, DataType::Int64, false)),
            ]
        };

        let schema = Arc::new(Schema::new(fields));
        Self {
            session,
            left_table,
            right_table,
            left_columns,
            right_columns,
            filter_op,
            schema,
            has_extra_cols,
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
        let left_df = if self.has_extra_cols {
            self.session.table(&self.left_table).await?
        } else {
            self.session
                .table(&self.left_table)
                .await?
                .select_columns(&[
                    &self.left_columns.0,
                    &self.left_columns.1,
                    &self.left_columns.2,
                ])?
        };
        let left_plan = left_df.create_physical_plan().await?;

        let right_batches = self
            .session
            .table(&self.right_table)
            .await?
            .select_columns(&[
                &self.right_columns.0,
                &self.right_columns.1,
                &self.right_columns.2,
            ])?
            .collect()
            .await?;
        let right_groups = Arc::new(build_right_groups(right_batches, &self.right_columns)?);

        let output_partitions = left_plan.output_partitioning().partition_count();

        Ok(Arc::new(SubtractExec {
            schema: self.schema.clone(),
            left: left_plan,
            left_columns: Arc::new(self.left_columns.clone()),
            right_groups,
            strict: self.filter_op == FilterOp::Strict,
            has_extra_cols: self.has_extra_cols,
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
    left_columns: Arc<(String, String, String)>,
    right_groups: Arc<AHashMap<String, Vec<(i64, i64)>>>,
    strict: bool,
    has_extra_cols: bool,
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
        vec![&self.left]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        if children.len() != 1 {
            return Err(DataFusionError::Internal(
                "SubtractExec expects exactly one child plan".to_string(),
            ));
        }

        Ok(Arc::new(Self {
            schema: self.schema.clone(),
            left: Arc::clone(&children[0]),
            left_columns: Arc::clone(&self.left_columns),
            right_groups: Arc::clone(&self.right_groups),
            strict: self.strict,
            has_extra_cols: self.has_extra_cols,
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
        let left = self.left.execute(partition, context)?;

        if self.has_extra_cols {
            let left_input_schema = left.schema();
            Ok(Box::pin(SubtractStreamExtra {
                schema: self.schema.clone(),
                left_columns: Arc::clone(&self.left_columns),
                left_collector: FullBatchCollector::new(
                    left,
                    Arc::clone(&self.left_columns),
                    left_input_schema,
                ),
                strict: self.strict,
                phase: SubtractPhase::CollectLeft,
                right_groups: Arc::clone(&self.right_groups),
                left_groups: Vec::new(),
                concatenated: None,
                start_col_idx: 0,
                end_col_idx: 0,
                group_idx: 0,
                interval_idx: 0,
                right_cursor: 0,
                output_row_indices: Vec::new(),
                output_starts: Vec::new(),
                output_ends: Vec::new(),
                pending_rows: 0,
                batch_size,
            }))
        } else {
            Ok(Box::pin(SubtractStream {
                schema: self.schema.clone(),
                left_collector: StreamCollector::new(left, Arc::clone(&self.left_columns)),
                strict: self.strict,
                phase: SubtractPhase::CollectLeft,
                right_groups: Arc::clone(&self.right_groups),
                left_groups: Vec::new(),
                group_idx: 0,
                interval_idx: 0,
                right_cursor: 0,
                contig_builder: StringBuilder::new(),
                start_builder: Int64Builder::new(),
                end_builder: Int64Builder::new(),
                pending_rows: 0,
                batch_size,
            }))
        }
    }
}

enum SubtractPhase {
    CollectLeft,
    Emit,
    Done,
}

fn build_right_groups(
    batches: Vec<RecordBatch>,
    columns: &(String, String, String),
) -> Result<AHashMap<String, Vec<(i64, i64)>>> {
    let mut groups = AHashMap::<String, Vec<(i64, i64)>>::new();

    for batch in batches {
        if batch.num_rows() == 0 {
            continue;
        }

        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(&batch, (&columns.0, &columns.1, &columns.2))?;
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
        merge_mask_intervals(intervals);
    }

    Ok(groups)
}

fn merge_mask_intervals(intervals: &mut Vec<(i64, i64)>) {
    if intervals.len() <= 1 {
        return;
    }

    let mut write_idx = 0usize;
    for read_idx in 1..intervals.len() {
        let (start, end) = intervals[read_idx];
        if start <= intervals[write_idx].1 {
            if end > intervals[write_idx].1 {
                intervals[write_idx].1 = end;
            }
        } else {
            write_idx += 1;
            intervals[write_idx] = (start, end);
        }
    }
    intervals.truncate(write_idx + 1);
}

// ─── Fast path: 3-column left input ─────────────────────────────────────────

struct SubtractStream {
    schema: SchemaRef,
    left_collector: StreamCollector,
    strict: bool,
    phase: SubtractPhase,
    right_groups: Arc<AHashMap<String, Vec<(i64, i64)>>>,
    left_groups: Vec<(String, Vec<(i64, i64)>)>,
    group_idx: usize,
    interval_idx: usize,
    right_cursor: usize,
    contig_builder: StringBuilder,
    start_builder: Int64Builder,
    end_builder: Int64Builder,
    pending_rows: usize,
    batch_size: usize,
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
                SubtractPhase::CollectLeft => match ready!(this.left_collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.left_groups = this.left_collector.take_groups();
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
                    while this.group_idx < this.left_groups.len() {
                        let (ref contig, ref left_intervals) = this.left_groups[this.group_idx];
                        let right_intervals = this
                            .right_groups
                            .get(contig)
                            .map(Vec::as_slice)
                            .unwrap_or(&[]);

                        while this.interval_idx < left_intervals.len() {
                            let (ls, le) = left_intervals[this.interval_idx];
                            this.interval_idx += 1;

                            while this.right_cursor < right_intervals.len() {
                                let skip = if this.strict {
                                    right_intervals[this.right_cursor].1 <= ls
                                } else {
                                    right_intervals[this.right_cursor].1 < ls
                                };
                                if skip {
                                    this.right_cursor += 1;
                                } else {
                                    break;
                                }
                            }

                            let mut cursor = ls;
                            let mut j = this.right_cursor;
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

                            if this.pending_rows >= this.batch_size {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;
                        this.right_cursor = 0;

                        if this.pending_rows >= this.batch_size {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

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

// ─── Extra-columns path: left input has >3 columns ──────────────────────────

struct SubtractStreamExtra {
    schema: SchemaRef,
    left_columns: Arc<(String, String, String)>,
    left_collector: FullBatchCollector,
    strict: bool,
    phase: SubtractPhase,
    right_groups: Arc<AHashMap<String, Vec<(i64, i64)>>>,
    left_groups: IndexedGroups,
    concatenated: Option<RecordBatch>,
    start_col_idx: usize,
    end_col_idx: usize,
    group_idx: usize,
    interval_idx: usize,
    right_cursor: usize,
    output_row_indices: Vec<u32>,
    output_starts: Vec<i64>,
    output_ends: Vec<i64>,
    pending_rows: usize,
    batch_size: usize,
}

impl SubtractStreamExtra {
    fn flush_builders(&mut self) -> Result<RecordBatch> {
        let concatenated = self.concatenated.as_ref().ok_or_else(|| {
            DataFusionError::Internal("FullBatchCollector: no concatenated batch".to_string())
        })?;

        let indices = UInt32Array::from(std::mem::take(&mut self.output_row_indices));
        let starts = Int64Array::from(std::mem::take(&mut self.output_starts));
        let ends = Int64Array::from(std::mem::take(&mut self.output_ends));

        let mut columns: Vec<Arc<dyn datafusion::arrow::array::Array>> =
            Vec::with_capacity(concatenated.num_columns());
        for (col_idx, col) in concatenated.columns().iter().enumerate() {
            if col_idx == self.start_col_idx {
                columns.push(Arc::new(starts.clone()));
            } else if col_idx == self.end_col_idx {
                columns.push(Arc::new(ends.clone()));
            } else {
                columns.push(take(col.as_ref(), &indices, None)?);
            }
        }

        self.pending_rows = 0;
        RecordBatch::try_new(self.schema.clone(), columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for SubtractStreamExtra {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                SubtractPhase::CollectLeft => match ready!(this.left_collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.left_groups = this.left_collector.take_groups();
                        this.concatenated = this.left_collector.take_concatenated();
                        if let Some(ref concat) = this.concatenated {
                            let schema = concat.schema();
                            this.start_col_idx = schema
                                .index_of(&this.left_columns.1)
                                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                            this.end_col_idx = schema
                                .index_of(&this.left_columns.2)
                                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                        }
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
                    while this.group_idx < this.left_groups.len() {
                        let (ref contig, ref left_intervals) = this.left_groups[this.group_idx];
                        let right_intervals = this
                            .right_groups
                            .get(contig)
                            .map(Vec::as_slice)
                            .unwrap_or(&[]);

                        while this.interval_idx < left_intervals.len() {
                            let (ls, le, row_idx) = left_intervals[this.interval_idx];
                            this.interval_idx += 1;

                            while this.right_cursor < right_intervals.len() {
                                let skip = if this.strict {
                                    right_intervals[this.right_cursor].1 <= ls
                                } else {
                                    right_intervals[this.right_cursor].1 < ls
                                };
                                if skip {
                                    this.right_cursor += 1;
                                } else {
                                    break;
                                }
                            }

                            let mut cursor = ls;
                            let mut j = this.right_cursor;
                            while j < right_intervals.len() {
                                let (rs, re) = right_intervals[j];
                                let no_overlap = if this.strict { rs >= le } else { rs > le };
                                if no_overlap {
                                    break;
                                }

                                if rs > cursor {
                                    debug_assert!(
                                        row_idx <= u32::MAX as usize,
                                        "row index {row_idx} exceeds u32::MAX"
                                    );
                                    this.output_row_indices.push(row_idx as u32);
                                    this.output_starts.push(cursor);
                                    this.output_ends.push(rs);
                                    this.pending_rows += 1;
                                }
                                if re > cursor {
                                    cursor = re;
                                }
                                j += 1;
                            }

                            if cursor < le {
                                debug_assert!(
                                    row_idx <= u32::MAX as usize,
                                    "row index {row_idx} exceeds u32::MAX"
                                );
                                this.output_row_indices.push(row_idx as u32);
                                this.output_starts.push(cursor);
                                this.output_ends.push(le);
                                this.pending_rows += 1;
                            }

                            if this.pending_rows >= this.batch_size {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;
                        this.right_cursor = 0;

                        if this.pending_rows >= this.batch_size {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

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

impl RecordBatchStream for SubtractStreamExtra {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}
