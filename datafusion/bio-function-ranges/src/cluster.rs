use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll, Waker};

use ahash::AHashMap;
use async_trait::async_trait;
use datafusion::arrow::array::{Int64Array, Int64Builder, RecordBatch, StringBuilder, UInt32Array};
use datafusion::arrow::compute::take;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::expressions::Column;
use datafusion::physical_expr::{Distribution, EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use futures::{Stream, ready};
use parking_lot::Mutex;

use crate::filter_op::FilterOp;
use crate::grouped_stream::{FullBatchCollector, IndexedGroups, StreamCollector};

pub struct ClusterProvider {
    session: Arc<SessionContext>,
    table: String,
    columns: (String, String, String),
    min_dist: i64,
    filter_op: FilterOp,
    schema: SchemaRef,
    input_schema: Schema,
    has_extra_cols: bool,
}

impl ClusterProvider {
    pub fn new(
        session: Arc<SessionContext>,
        table: String,
        columns: (String, String, String),
        min_dist: i64,
        filter_op: FilterOp,
        input_schema: Schema,
    ) -> Self {
        let has_extra_cols = input_schema.fields().len() > 3;

        let mut fields: Vec<Arc<Field>> = if has_extra_cols {
            // Preserve all input fields as-is (take kernel returns same types)
            input_schema.fields().iter().map(Arc::clone).collect()
        } else {
            vec![
                Arc::new(Field::new(&columns.0, DataType::Utf8, false)),
                Arc::new(Field::new(&columns.1, DataType::Int64, false)),
                Arc::new(Field::new(&columns.2, DataType::Int64, false)),
            ]
        };
        fields.push(Arc::new(Field::new("cluster", DataType::Int64, false)));
        fields.push(Arc::new(Field::new(
            "cluster_start",
            DataType::Int64,
            false,
        )));
        fields.push(Arc::new(Field::new("cluster_end", DataType::Int64, false)));

        let schema = Arc::new(Schema::new(fields));
        Self {
            session,
            table,
            columns,
            min_dist,
            filter_op,
            schema,
            input_schema,
            has_extra_cols,
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
        let input_df = if self.has_extra_cols {
            self.session.table(&self.table).await?
        } else {
            self.session.table(&self.table).await?.select_columns(&[
                &self.columns.0,
                &self.columns.1,
                &self.columns.2,
            ])?
        };
        let input_plan = input_df.create_physical_plan().await?;

        let contig_col_idx = if self.has_extra_cols {
            self.input_schema.index_of(&self.columns.0)?
        } else {
            0
        };

        let output_partitions = input_plan.output_partitioning().partition_count();

        Ok(Arc::new(ClusterExec {
            schema: self.schema.clone(),
            input: input_plan,
            columns: Arc::new(self.columns.clone()),
            cluster_id_coordinator: Arc::new(ClusterIdCoordinator::new(output_partitions)),
            contig_col_idx,
            min_dist: self.min_dist,
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
struct ClusterExec {
    schema: SchemaRef,
    input: Arc<dyn ExecutionPlan>,
    columns: Arc<(String, String, String)>,
    cluster_id_coordinator: Arc<ClusterIdCoordinator>,
    contig_col_idx: usize,
    min_dist: i64,
    strict: bool,
    has_extra_cols: bool,
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

    fn required_input_distribution(&self) -> Vec<Distribution> {
        vec![Distribution::HashPartitioned(vec![Arc::new(Column::new(
            self.columns.0.as_str(),
            self.contig_col_idx,
        ))])]
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
            cluster_id_coordinator: Arc::new(ClusterIdCoordinator::new(
                children[0].output_partitioning().partition_count(),
            )),
            contig_col_idx: self.contig_col_idx,
            min_dist: self.min_dist,
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
        let input = self.input.execute(partition, context)?;
        let cluster_id_coordinator = Arc::clone(&self.cluster_id_coordinator);

        if self.has_extra_cols {
            let input_schema = input.schema();
            Ok(Box::pin(ClusterStreamExtra {
                schema: self.schema.clone(),
                collector: FullBatchCollector::new(input, Arc::clone(&self.columns), input_schema),
                cluster_id_coordinator,
                cluster_counts: None,
                cluster_offsets: None,
                min_dist: self.min_dist,
                strict: self.strict,
                phase: ClusterPhase::Collecting,
                groups: Vec::new(),
                concatenated: None,
                partition,
                group_idx: 0,
                interval_idx: 0,
                cluster_id: 0,
                cluster_start: 0,
                cluster_end: 0,
                pending_intervals: Vec::new(),
                output_row_indices: Vec::new(),
                output_cluster_ids: Vec::new(),
                output_cluster_starts: Vec::new(),
                output_cluster_ends: Vec::new(),
                pending_rows: 0,
                batch_size,
            }))
        } else {
            Ok(Box::pin(ClusterStream {
                schema: self.schema.clone(),
                collector: StreamCollector::new(input, Arc::clone(&self.columns)),
                cluster_id_coordinator,
                cluster_counts: None,
                cluster_offsets: None,
                min_dist: self.min_dist,
                strict: self.strict,
                phase: ClusterPhase::Collecting,
                groups: Vec::new(),
                partition,
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
                batch_size,
            }))
        }
    }
}

enum ClusterPhase {
    Collecting,
    AwaitingOffsets,
    Emitting,
    Done,
}

struct ClusterIdCoordinator {
    state: Mutex<ClusterIdCoordinatorState>,
}

impl ClusterIdCoordinator {
    fn new(partition_count: usize) -> Self {
        Self {
            state: Mutex::new(ClusterIdCoordinatorState {
                reports: vec![None; partition_count],
                registered: 0,
                offsets: None,
                error: None,
                wakers: Vec::new(),
            }),
        }
    }

    fn poll_cluster_offsets(
        &self,
        partition: usize,
        counts: &[(String, usize)],
        cx: &Context<'_>,
    ) -> Poll<Result<Arc<AHashMap<String, i64>>>> {
        let mut state = self.state.lock();

        if let Some(offsets) = &state.offsets {
            return Poll::Ready(Ok(Arc::clone(offsets)));
        }

        if let Some(error) = &state.error {
            return Poll::Ready(Err(DataFusionError::Internal(error.clone())));
        }

        if partition >= state.reports.len() {
            return Poll::Ready(Err(DataFusionError::Internal(format!(
                "cluster partition {partition} out of range for {} output partitions",
                state.reports.len()
            ))));
        }

        if let Some(existing) = &state.reports[partition] {
            if existing.as_slice() != counts {
                return Poll::Ready(Err(DataFusionError::Internal(format!(
                    "cluster partition {partition} reported inconsistent cluster counts"
                ))));
            }
        } else {
            state.reports[partition] = Some(counts.to_vec());
            state.registered += 1;
        }

        if state.registered == state.reports.len() {
            match Self::build_offsets(&state.reports) {
                Ok(offsets) => {
                    let offsets = Arc::new(offsets);
                    state.offsets = Some(Arc::clone(&offsets));
                    for waker in state.wakers.drain(..) {
                        waker.wake();
                    }
                    return Poll::Ready(Ok(offsets));
                }
                Err(error) => {
                    let message = error.to_string();
                    state.error = Some(message.clone());
                    for waker in state.wakers.drain(..) {
                        waker.wake();
                    }
                    return Poll::Ready(Err(DataFusionError::Internal(message)));
                }
            }
        }

        state.wakers.push(cx.waker().clone());
        Poll::Pending
    }

    fn build_offsets(reports: &[Option<Vec<(String, usize)>>]) -> Result<AHashMap<String, i64>> {
        let mut contig_counts = Vec::with_capacity(reports.iter().flatten().map(Vec::len).sum());
        for report in reports {
            let report = report.as_ref().ok_or_else(|| {
                DataFusionError::Internal(
                    "cluster offset coordinator finalized before all partitions reported"
                        .to_string(),
                )
            })?;
            contig_counts.extend(report.iter().cloned());
        }

        contig_counts.sort_unstable_by(|left, right| left.0.cmp(&right.0));

        let mut offsets = AHashMap::with_capacity(contig_counts.len());
        let mut next_cluster_id = 0_i64;

        for (contig, cluster_count) in contig_counts {
            if offsets.contains_key(&contig) {
                return Err(DataFusionError::Internal(format!(
                    "cluster offset coordinator received duplicate contig '{contig}' across partitions"
                )));
            }

            offsets.insert(contig, next_cluster_id);
            let cluster_count = i64::try_from(cluster_count).map_err(|_| {
                DataFusionError::Internal(
                    "cluster count exceeded i64::MAX while computing global offsets".to_string(),
                )
            })?;
            next_cluster_id = next_cluster_id.checked_add(cluster_count).ok_or_else(|| {
                DataFusionError::Internal(
                    "cluster id overflow while computing global offsets".to_string(),
                )
            })?;
        }

        Ok(offsets)
    }
}

impl Debug for ClusterIdCoordinator {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let state = self.state.lock();
        f.debug_struct("ClusterIdCoordinator")
            .field("partition_count", &state.reports.len())
            .field("registered", &state.registered)
            .field("ready", &state.offsets.is_some())
            .finish()
    }
}

struct ClusterIdCoordinatorState {
    reports: Vec<Option<Vec<(String, usize)>>>,
    registered: usize,
    offsets: Option<Arc<AHashMap<String, i64>>>,
    error: Option<String>,
    wakers: Vec<Waker>,
}

fn count_clusters<I>(intervals: I, min_dist: i64, strict: bool) -> usize
where
    I: Iterator<Item = (i64, i64)>,
{
    let mut cluster_count = 0;
    let mut cluster_end = None;

    for (start, end) in intervals {
        match cluster_end {
            None => {
                cluster_count += 1;
                cluster_end = Some(end);
            }
            Some(current_end) => {
                let boundary = current_end.saturating_add(min_dist);
                let merge_condition = if strict {
                    start < boundary
                } else {
                    start <= boundary
                };

                if merge_condition {
                    if end > current_end {
                        cluster_end = Some(end);
                    }
                } else {
                    cluster_count += 1;
                    cluster_end = Some(end);
                }
            }
        }
    }

    cluster_count
}

fn collect_cluster_counts(
    groups: &[(String, Vec<(i64, i64)>)],
    min_dist: i64,
    strict: bool,
) -> Vec<(String, usize)> {
    groups
        .iter()
        .map(|(contig, intervals)| {
            (
                contig.clone(),
                count_clusters(intervals.iter().copied(), min_dist, strict),
            )
        })
        .collect()
}

fn collect_indexed_cluster_counts(
    groups: &IndexedGroups,
    min_dist: i64,
    strict: bool,
) -> Vec<(String, usize)> {
    groups
        .iter()
        .map(|(contig, intervals)| {
            (
                contig.clone(),
                count_clusters(
                    intervals.iter().map(|&(start, end, _)| (start, end)),
                    min_dist,
                    strict,
                ),
            )
        })
        .collect()
}

// ─── Fast path: 3-column input (no extra columns) ───────────────────────────

struct ClusterStream {
    schema: SchemaRef,
    collector: StreamCollector,
    cluster_id_coordinator: Arc<ClusterIdCoordinator>,
    cluster_counts: Option<Vec<(String, usize)>>,
    cluster_offsets: Option<Arc<AHashMap<String, i64>>>,
    min_dist: i64,
    strict: bool,
    phase: ClusterPhase,
    groups: Vec<(String, Vec<(i64, i64)>)>,
    partition: usize,
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
    batch_size: usize,
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

    fn set_group_cluster_base(&mut self, contig: &str) -> Result<()> {
        let cluster_offsets = self.cluster_offsets.as_ref().ok_or_else(|| {
            DataFusionError::Internal("cluster stream missing global cluster offsets".to_string())
        })?;
        self.cluster_id = *cluster_offsets.get(contig).ok_or_else(|| {
            DataFusionError::Internal(format!(
                "cluster stream missing global cluster offset for contig '{contig}'"
            ))
        })?;
        Ok(())
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
                        this.cluster_counts = Some(collect_cluster_counts(
                            &this.groups,
                            this.min_dist,
                            this.strict,
                        ));
                        this.phase = ClusterPhase::AwaitingOffsets;
                    }
                    Ok(false) => unreachable!(),
                    Err(e) => {
                        this.phase = ClusterPhase::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                },
                ClusterPhase::AwaitingOffsets => {
                    let cluster_counts = this.cluster_counts.as_deref().ok_or_else(|| {
                        DataFusionError::Internal(
                            "cluster stream missing local cluster counts".to_string(),
                        )
                    });
                    let cluster_counts = match cluster_counts {
                        Ok(cluster_counts) => cluster_counts,
                        Err(error) => {
                            this.phase = ClusterPhase::Done;
                            return Poll::Ready(Some(Err(error)));
                        }
                    };

                    match this.cluster_id_coordinator.poll_cluster_offsets(
                        this.partition,
                        cluster_counts,
                        cx,
                    ) {
                        Poll::Ready(Ok(cluster_offsets)) => {
                            this.cluster_offsets = Some(cluster_offsets);
                            this.group_idx = 0;
                            this.interval_idx = 0;
                            this.phase = ClusterPhase::Emitting;
                        }
                        Poll::Ready(Err(error)) => {
                            this.phase = ClusterPhase::Done;
                            return Poll::Ready(Some(Err(error)));
                        }
                        Poll::Pending => return Poll::Pending,
                    }
                }
                ClusterPhase::Emitting => {
                    while this.group_idx < this.groups.len() {
                        let contig = this.groups[this.group_idx].0.clone();
                        let interval_count = this.groups[this.group_idx].1.len();

                        if this.interval_idx == 0
                            && this.pending_intervals.is_empty()
                            && let Err(error) = this.set_group_cluster_base(&contig)
                        {
                            this.phase = ClusterPhase::Done;
                            return Poll::Ready(Some(Err(error)));
                        }

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

                            if this.pending_rows >= this.batch_size {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        if !this.pending_intervals.is_empty() {
                            this.flush_pending_cluster(&contig);
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;

                        if this.pending_rows >= this.batch_size {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

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

// ─── Extra-columns path: input has >3 columns ───────────────────────────────

struct ClusterStreamExtra {
    schema: SchemaRef,
    collector: FullBatchCollector,
    cluster_id_coordinator: Arc<ClusterIdCoordinator>,
    cluster_counts: Option<Vec<(String, usize)>>,
    cluster_offsets: Option<Arc<AHashMap<String, i64>>>,
    min_dist: i64,
    strict: bool,
    phase: ClusterPhase,
    groups: IndexedGroups,
    concatenated: Option<RecordBatch>,
    partition: usize,
    group_idx: usize,
    interval_idx: usize,
    cluster_id: i64,
    cluster_start: i64,
    cluster_end: i64,
    pending_intervals: Vec<(i64, i64, usize)>,
    output_row_indices: Vec<u32>,
    output_cluster_ids: Vec<i64>,
    output_cluster_starts: Vec<i64>,
    output_cluster_ends: Vec<i64>,
    pending_rows: usize,
    batch_size: usize,
}

impl ClusterStreamExtra {
    fn flush_pending_cluster(&mut self) {
        let cluster_id = self.cluster_id;
        let cluster_start = self.cluster_start;
        let cluster_end = self.cluster_end;
        for &(_s, _e, row_idx) in &self.pending_intervals {
            debug_assert!(
                row_idx <= u32::MAX as usize,
                "row index {row_idx} exceeds u32::MAX"
            );
            self.output_row_indices.push(row_idx as u32);
            self.output_cluster_ids.push(cluster_id);
            self.output_cluster_starts.push(cluster_start);
            self.output_cluster_ends.push(cluster_end);
        }
        self.pending_rows += self.pending_intervals.len();
        self.pending_intervals.clear();
        self.cluster_id += 1;
    }

    fn set_group_cluster_base(&mut self, contig: &str) -> Result<()> {
        let cluster_offsets = self.cluster_offsets.as_ref().ok_or_else(|| {
            DataFusionError::Internal(
                "cluster stream with extra columns missing global cluster offsets".to_string(),
            )
        })?;
        self.cluster_id = *cluster_offsets.get(contig).ok_or_else(|| {
            DataFusionError::Internal(format!(
                "cluster stream with extra columns missing global cluster offset for contig '{contig}'"
            ))
        })?;
        Ok(())
    }

    fn flush_builders(&mut self) -> Result<RecordBatch> {
        let concatenated = self.concatenated.as_ref().ok_or_else(|| {
            DataFusionError::Internal("FullBatchCollector: no concatenated batch".to_string())
        })?;

        let indices = UInt32Array::from(std::mem::take(&mut self.output_row_indices));
        let mut columns: Vec<Arc<dyn datafusion::arrow::array::Array>> =
            Vec::with_capacity(concatenated.num_columns() + 3);
        for col in concatenated.columns() {
            columns.push(take(col.as_ref(), &indices, None)?);
        }

        columns.push(Arc::new(Int64Array::from(std::mem::take(
            &mut self.output_cluster_ids,
        ))));
        columns.push(Arc::new(Int64Array::from(std::mem::take(
            &mut self.output_cluster_starts,
        ))));
        columns.push(Arc::new(Int64Array::from(std::mem::take(
            &mut self.output_cluster_ends,
        ))));

        self.pending_rows = 0;
        RecordBatch::try_new(self.schema.clone(), columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for ClusterStreamExtra {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match this.phase {
                ClusterPhase::Collecting => match ready!(this.collector.poll_collect(cx)) {
                    Ok(true) => {
                        this.groups = this.collector.take_groups();
                        this.concatenated = this.collector.take_concatenated();
                        this.cluster_counts = Some(collect_indexed_cluster_counts(
                            &this.groups,
                            this.min_dist,
                            this.strict,
                        ));
                        this.phase = ClusterPhase::AwaitingOffsets;
                    }
                    Ok(false) => unreachable!(),
                    Err(e) => {
                        this.phase = ClusterPhase::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                },
                ClusterPhase::AwaitingOffsets => {
                    let cluster_counts = this.cluster_counts.as_deref().ok_or_else(|| {
                        DataFusionError::Internal(
                            "cluster stream with extra columns missing local cluster counts"
                                .to_string(),
                        )
                    });
                    let cluster_counts = match cluster_counts {
                        Ok(cluster_counts) => cluster_counts,
                        Err(error) => {
                            this.phase = ClusterPhase::Done;
                            return Poll::Ready(Some(Err(error)));
                        }
                    };

                    match this.cluster_id_coordinator.poll_cluster_offsets(
                        this.partition,
                        cluster_counts,
                        cx,
                    ) {
                        Poll::Ready(Ok(cluster_offsets)) => {
                            this.cluster_offsets = Some(cluster_offsets);
                            this.group_idx = 0;
                            this.interval_idx = 0;
                            this.phase = ClusterPhase::Emitting;
                        }
                        Poll::Ready(Err(error)) => {
                            this.phase = ClusterPhase::Done;
                            return Poll::Ready(Some(Err(error)));
                        }
                        Poll::Pending => return Poll::Pending,
                    }
                }
                ClusterPhase::Emitting => {
                    while this.group_idx < this.groups.len() {
                        let contig = this.groups[this.group_idx].0.clone();
                        let interval_count = this.groups[this.group_idx].1.len();

                        if this.interval_idx == 0
                            && this.pending_intervals.is_empty()
                            && let Err(error) = this.set_group_cluster_base(&contig)
                        {
                            this.phase = ClusterPhase::Done;
                            return Poll::Ready(Some(Err(error)));
                        }

                        while this.interval_idx < interval_count {
                            let (s, e, row_idx) = this.groups[this.group_idx].1[this.interval_idx];
                            this.interval_idx += 1;

                            if this.pending_intervals.is_empty() {
                                this.cluster_start = s;
                                this.cluster_end = e;
                                this.pending_intervals.push((s, e, row_idx));
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
                                    this.pending_intervals.push((s, e, row_idx));
                                } else {
                                    this.flush_pending_cluster();
                                    this.cluster_start = s;
                                    this.cluster_end = e;
                                    this.pending_intervals.push((s, e, row_idx));
                                }
                            }

                            if this.pending_rows >= this.batch_size {
                                return Poll::Ready(Some(this.flush_builders()));
                            }
                        }

                        if !this.pending_intervals.is_empty() {
                            this.flush_pending_cluster();
                        }

                        this.group_idx += 1;
                        this.interval_idx = 0;

                        if this.pending_rows >= this.batch_size {
                            return Poll::Ready(Some(this.flush_builders()));
                        }
                    }

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

impl RecordBatchStream for ClusterStreamExtra {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}
