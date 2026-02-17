use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::array::{Array, PrimitiveArray, RecordBatch};
use datafusion::arrow::buffer::NullBuffer;
use datafusion::arrow::datatypes::{Field, Schema, SchemaRef, UInt32Type};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::logical_expr::Expr;
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::SessionContext;
use fnv::FnvHashMap;
use futures::StreamExt;
use futures::stream::BoxStream;

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;
use crate::nearest_index::{IntervalRecord, NearestIntervalIndex};

pub struct NearestProvider {
    session: Arc<SessionContext>,
    left_table: String,
    right_table: String,
    columns_1: (String, String, String),
    columns_2: (String, String, String),
    filter_op: FilterOp,
    include_overlaps: bool,
    k: usize,
    schema: SchemaRef,
}

impl NearestProvider {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        session: Arc<SessionContext>,
        left_table: String,
        right_table: String,
        left_table_schema: Schema,
        right_table_schema: Schema,
        columns_1: Vec<String>,
        columns_2: Vec<String>,
        filter_op: FilterOp,
        include_overlaps: bool,
        k: usize,
    ) -> Self {
        let mut fields = left_table_schema
            .fields()
            .iter()
            .map(|f| {
                Arc::new(Field::new(
                    format!("left_{}", f.name()),
                    f.data_type().clone(),
                    true,
                ))
            })
            .collect::<Vec<_>>();
        fields.extend(right_table_schema.fields().iter().map(|f| {
            Arc::new(Field::new(
                format!("right_{}", f.name()),
                f.data_type().clone(),
                true,
            ))
        }));
        let schema = Arc::new(Schema::new(fields));

        Self {
            session,
            left_table,
            right_table,
            schema,
            columns_1: (
                columns_1[0].clone(),
                columns_1[1].clone(),
                columns_1[2].clone(),
            ),
            columns_2: (
                columns_2[0].clone(),
                columns_2[1].clone(),
                columns_2[2].clone(),
            ),
            filter_op,
            include_overlaps,
            k,
        }
    }
}

impl Debug for NearestProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "NearestProvider {{ left: {}, right: {}, k: {}, include_overlaps: {} }}",
            self.left_table, self.right_table, self.k, self.include_overlaps
        )
    }
}

#[async_trait]
impl TableProvider for NearestProvider {
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

        let left_df = self.session.table(self.left_table.clone()).await?;
        let left_schema = Arc::new(left_df.schema().as_arrow().clone());
        let left_batches = left_df.collect().await?;
        let left_batch = datafusion::arrow::compute::concat_batches(&left_schema, &left_batches)?;

        let indexes = Arc::new(build_nearest_indexes(
            &left_batch,
            (&self.columns_1.0, &self.columns_1.1, &self.columns_1.2),
        )?);

        let right_df = self.session.table(self.right_table.clone()).await?;
        let right_plan = right_df.create_physical_plan().await?;
        let right_plan: Arc<dyn ExecutionPlan> =
            if right_plan.output_partitioning().partition_count() == target_partitions {
                right_plan
            } else {
                Arc::new(RepartitionExec::try_new(
                    right_plan,
                    Partitioning::RoundRobinBatch(target_partitions),
                )?)
            };
        let output_partitions = right_plan.output_partitioning().partition_count();

        Ok(Arc::new(NearestExec {
            schema: self.schema.clone(),
            left_batch: Arc::new(left_batch),
            indexes,
            right: right_plan,
            columns_2: Arc::new(self.columns_2.clone()),
            filter_op: self.filter_op.clone(),
            include_overlaps: self.include_overlaps,
            k: self.k,
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
struct NearestExec {
    schema: SchemaRef,
    left_batch: Arc<RecordBatch>,
    indexes: Arc<FnvHashMap<String, NearestIntervalIndex>>,
    right: Arc<dyn ExecutionPlan>,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    include_overlaps: bool,
    k: usize,
    cache: PlanProperties,
}

impl DisplayAs for NearestExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "NearestExec: k={}, include_overlaps={}",
            self.k, self.include_overlaps
        )
    }
}

impl ExecutionPlan for NearestExec {
    fn name(&self) -> &str {
        "NearestExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn properties(&self) -> &PlanProperties {
        &self.cache
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.right]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        if children.len() != 1 {
            return Err(DataFusionError::Internal(
                "NearestExec expects exactly one child plan".to_string(),
            ));
        }

        Ok(Arc::new(NearestExec {
            schema: self.schema.clone(),
            left_batch: Arc::clone(&self.left_batch),
            indexes: Arc::clone(&self.indexes),
            right: Arc::clone(&children[0]),
            columns_2: Arc::clone(&self.columns_2),
            filter_op: self.filter_op.clone(),
            include_overlaps: self.include_overlaps,
            k: self.k,
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
        get_nearest_stream(
            Arc::clone(&self.right),
            Arc::clone(&self.left_batch),
            Arc::clone(&self.indexes),
            self.schema.clone(),
            Arc::clone(&self.columns_2),
            self.filter_op.clone(),
            self.include_overlaps,
            self.k,
            partition,
            context,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn get_nearest_stream(
    right_plan: Arc<dyn ExecutionPlan>,
    left_batch: Arc<RecordBatch>,
    indexes: Arc<FnvHashMap<String, NearestIntervalIndex>>,
    new_schema: SchemaRef,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    include_overlaps: bool,
    k: usize,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let partition_stream = right_plan.execute(partition, context)?;
    let schema_for_closure = new_schema.clone();
    let strict_filter = filter_op == FilterOp::Strict;

    let iter = partition_stream.map(move |rb| match rb {
        Ok(rb) => {
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(&rb, (&columns_2.0, &columns_2.1, &columns_2.2))?;

            let mut left_indices = Vec::<u32>::new();
            let mut validity = Vec::<bool>::new();
            let mut right_indices = Vec::<u32>::new();
            let mut nearest_buf = Vec::<usize>::with_capacity(k.max(1));

            for i in 0..rb.num_rows() {
                let contig = contig_arr.value(i);
                let mut query_start = start_arr.value(i)?;
                let mut query_end = end_arr.value(i)?;

                if strict_filter {
                    query_start += 1;
                    query_end -= 1;
                }

                nearest_buf.clear();

                if let Some(index) = indexes.get(contig) {
                    index.nearest_k(
                        query_start,
                        query_end,
                        k,
                        include_overlaps,
                        &mut nearest_buf,
                    );
                }

                if nearest_buf.is_empty() {
                    left_indices.push(0);
                    validity.push(false);
                    right_indices.push(i as u32);
                } else {
                    for &pos in &nearest_buf {
                        let pos_u32 = u32::try_from(pos).map_err(|_| {
                            DataFusionError::Execution(format!(
                                "left row index {pos} exceeds UInt32 index capacity"
                            ))
                        })?;
                        left_indices.push(pos_u32);
                        validity.push(true);
                        right_indices.push(i as u32);
                    }
                }
            }

            let left_index_array = if validity.iter().all(|&v| v) {
                PrimitiveArray::<UInt32Type>::from(left_indices)
            } else {
                let null_buffer = NullBuffer::from(validity);
                PrimitiveArray::<UInt32Type>::new(left_indices.into(), Some(null_buffer))
            };
            let right_index_array = PrimitiveArray::<UInt32Type>::from(right_indices);

            let mut columns: Vec<Arc<dyn Array>> =
                Vec::with_capacity(left_batch.num_columns() + rb.num_columns());

            for array in left_batch.columns() {
                columns.push(datafusion::arrow::compute::take(
                    array,
                    &left_index_array,
                    None,
                )?);
            }
            for array in rb.columns() {
                columns.push(datafusion::arrow::compute::take(
                    array,
                    &right_index_array,
                    None,
                )?);
            }

            RecordBatch::try_new(schema_for_closure.clone(), columns)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
        }
        Err(e) => Err(e),
    });

    let adapted_stream = RecordBatchStreamAdapter::new(new_schema, Box::pin(iter) as BoxStream<_>);
    Ok(Box::pin(adapted_stream))
}

fn build_nearest_indexes(
    batch: &RecordBatch,
    columns: (&str, &str, &str),
) -> Result<FnvHashMap<String, NearestIntervalIndex>> {
    let (contig_arr, start_arr, end_arr) = get_join_col_arrays(batch, columns)?;

    let mut grouped = FnvHashMap::<String, Vec<IntervalRecord>>::default();
    for i in 0..batch.num_rows() {
        let contig = contig_arr.value(i);
        grouped
            .entry(contig.to_string())
            .or_default()
            .push(IntervalRecord {
                start: start_arr.value(i)?,
                end: end_arr.value(i)?,
                position: i,
            });
    }

    Ok(grouped
        .into_iter()
        .map(|(k, records)| (k, NearestIntervalIndex::from_records(records)))
        .collect())
}
