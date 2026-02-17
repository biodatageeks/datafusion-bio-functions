use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
use coitrees::COITree;
use datafusion::arrow::datatypes::{DataType, Field, FieldRef, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use fnv::FnvHashMap;
use futures::TryStreamExt;

use crate::filter_op::FilterOp;
use crate::interval_tree::{build_coitree_from_batches, get_stream};

pub struct CountOverlapsProvider {
    session: Arc<SessionContext>,
    left_table: String,
    right_table: String,
    columns_1: (String, String, String),
    columns_2: (String, String, String),
    filter_op: FilterOp,
    coverage: bool,
    schema: SchemaRef,
}

impl CountOverlapsProvider {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        session: Arc<SessionContext>,
        left_table: String,
        right_table: String,
        right_table_schema: Schema,
        columns_1: Vec<String>,
        columns_2: Vec<String>,
        filter_op: FilterOp,
        coverage: bool,
    ) -> Self {
        Self {
            session,
            left_table,
            right_table,
            schema: {
                let mut fields = right_table_schema.fields().to_vec();
                let name = if coverage { "coverage" } else { "count" };
                let new_field = Field::new(name, DataType::Int64, false);
                fields.push(FieldRef::new(new_field));
                Arc::new(Schema::new(fields))
            },
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
            coverage,
        }
    }
}

impl Debug for CountOverlapsProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CountOverlapsProvider {{ left: {}, right: {}, coverage: {} }}",
            self.left_table, self.right_table, self.coverage
        )
    }
}

#[async_trait]
impl TableProvider for CountOverlapsProvider {
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

        let left_table = self
            .session
            .table(self.left_table.clone())
            .await?
            .select_columns(&[&self.columns_1.0, &self.columns_1.1, &self.columns_1.2])?
            .collect()
            .await?;

        let trees = Arc::new(build_coitree_from_batches(
            left_table,
            (&self.columns_1.0, &self.columns_1.1, &self.columns_1.2),
            self.coverage,
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

        Ok(Arc::new(CountOverlapsExec {
            schema: self.schema().clone(),
            trees,
            right: right_plan,
            columns_2: Arc::new(self.columns_2.clone()),
            filter_op: self.filter_op.clone(),
            coverage: self.coverage,
            cache: PlanProperties::new(
                EquivalenceProperties::new(self.schema().clone()),
                Partitioning::UnknownPartitioning(output_partitions),
                EmissionType::Final,
                Boundedness::Bounded,
            ),
        }))
    }
}

struct CountOverlapsExec {
    schema: SchemaRef,
    trees: Arc<FnvHashMap<String, COITree<(), u32>>>,
    right: Arc<dyn ExecutionPlan>,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    coverage: bool,
    cache: PlanProperties,
}

impl Debug for CountOverlapsExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "CountOverlapsExec")
    }
}

impl DisplayAs for CountOverlapsExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter) -> std::fmt::Result {
        write!(f, "CountOverlapsExec")
    }
}

impl ExecutionPlan for CountOverlapsExec {
    fn name(&self) -> &str {
        "CountOverlapsExec"
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
                "CountOverlapsExec expects exactly one child plan".to_string(),
            ));
        }

        Ok(Arc::new(CountOverlapsExec {
            schema: self.schema.clone(),
            trees: Arc::clone(&self.trees),
            right: Arc::clone(&children[0]),
            columns_2: Arc::clone(&self.columns_2),
            filter_op: self.filter_op.clone(),
            coverage: self.coverage,
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
        let fut = get_stream(
            Arc::clone(&self.right),
            self.trees.clone(),
            self.schema.clone(),
            Arc::clone(&self.columns_2),
            self.filter_op.clone(),
            self.coverage,
            partition,
            context,
        );
        let stream = futures::stream::once(fut).try_flatten();
        let schema = self.schema.clone();
        Ok(Box::pin(RecordBatchStreamAdapter::new(schema, stream)))
    }
}
