use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use ahash::AHashMap;
use async_trait::async_trait;
use datafusion::arrow::array::{Int64Builder, RecordBatch, StringBuilder};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::expressions::Column;
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::common::collect;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use futures::stream::{BoxStream, once};

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
                EmissionType::Final,
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
                EmissionType::Final,
                Boundedness::Bounded,
            ),
        }))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        get_cluster_stream(
            Arc::clone(&self.input),
            self.schema.clone(),
            Arc::clone(&self.columns),
            self.min_dist,
            self.strict,
            partition,
            context,
        )
    }
}

fn get_cluster_stream(
    input_plan: Arc<dyn ExecutionPlan>,
    schema: SchemaRef,
    columns: Arc<(String, String, String)>,
    min_dist: i64,
    strict: bool,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let schema_for_closure = schema.clone();
    let once_stream = once(async move {
        let partition_stream = input_plan.execute(partition, context)?;
        let batches = collect(partition_stream).await?;
        build_cluster_batch(&batches, &schema_for_closure, &columns, min_dist, strict)
    });

    let adapted_stream =
        RecordBatchStreamAdapter::new(schema, Box::pin(once_stream) as BoxStream<'_, _>);
    Ok(Box::pin(adapted_stream))
}

fn build_cluster_batch(
    batches: &[RecordBatch],
    schema: &SchemaRef,
    columns: &(String, String, String),
    min_dist: i64,
    strict: bool,
) -> Result<RecordBatch> {
    let mut groups: AHashMap<String, Vec<(i64, i64)>> = AHashMap::default();
    let mut n_rows = 0usize;
    let mut contig_bytes = 0usize;

    for batch in batches {
        n_rows += batch.num_rows();
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(batch, (&columns.0, &columns.1, &columns.2))?;
        let start_resolved = start_arr.resolve_i64()?;
        let end_resolved = end_arr.resolve_i64()?;
        let starts = &*start_resolved;
        let ends = &*end_resolved;
        for i in 0..batch.num_rows() {
            let contig = contig_arr.value(i);
            contig_bytes += contig.len();
            groups
                .entry(contig.to_string())
                .or_default()
                .push((starts[i], ends[i]));
        }
    }

    let mut contigs: Vec<String> = groups.keys().cloned().collect();
    contigs.sort_unstable();

    let mut contig_builder = StringBuilder::with_capacity(n_rows, contig_bytes);
    let mut start_builder = Int64Builder::with_capacity(n_rows);
    let mut end_builder = Int64Builder::with_capacity(n_rows);
    let mut cluster_builder = Int64Builder::with_capacity(n_rows);
    let mut cluster_start_builder = Int64Builder::with_capacity(n_rows);
    let mut cluster_end_builder = Int64Builder::with_capacity(n_rows);

    let mut cluster_id: i64 = 0;

    for contig in &contigs {
        let intervals = groups.get_mut(contig.as_str()).unwrap();
        intervals.sort_unstable();
        if intervals.is_empty() {
            continue;
        }

        let mut cluster_start_idx = 0usize;
        let mut cur_start = intervals[0].0;
        let mut cur_end = intervals[0].1;

        for i in 1..intervals.len() {
            let (s, e) = intervals[i];
            let boundary = cur_end.saturating_add(min_dist);
            let merge_condition = if strict { s < boundary } else { s <= boundary };
            if merge_condition {
                if e > cur_end {
                    cur_end = e;
                }
            } else {
                for &(cs, ce) in &intervals[cluster_start_idx..i] {
                    contig_builder.append_value(contig);
                    start_builder.append_value(cs);
                    end_builder.append_value(ce);
                    cluster_builder.append_value(cluster_id);
                    cluster_start_builder.append_value(cur_start);
                    cluster_end_builder.append_value(cur_end);
                }
                cluster_id += 1;
                cluster_start_idx = i;
                cur_start = s;
                cur_end = e;
            }
        }

        for &(cs, ce) in &intervals[cluster_start_idx..] {
            contig_builder.append_value(contig);
            start_builder.append_value(cs);
            end_builder.append_value(ce);
            cluster_builder.append_value(cluster_id);
            cluster_start_builder.append_value(cur_start);
            cluster_end_builder.append_value(cur_end);
        }
        cluster_id += 1;
    }

    RecordBatch::try_new(
        schema.clone(),
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
