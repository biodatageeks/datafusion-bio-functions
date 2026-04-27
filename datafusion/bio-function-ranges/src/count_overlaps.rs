use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use ahash::AHashMap;
use async_trait::async_trait;
use coitrees::COITree;
use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::datatypes::{DataType, Field, FieldRef, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::repartition::RepartitionExec;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, SessionContext};
use log::warn;

#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
use crate::count_overlaps_apple_gpu::{AppleGpuCountOverlapsBackend, get_apple_gpu_stream};
#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
use crate::count_overlaps_rank::CountOverlapsRankIndex;
use crate::filter_op::FilterOp;
use crate::interval_tree::{build_coitree_from_batches, get_stream};
use crate::session_context::{BioConfig, CountOverlapsBackendMode};

// Initial Auto threshold from the local Apple M3 Max 8-7 release benchmark:
// ex-rna (9,944,559 left intervals) vs ex-anno (1,194,285 right intervals).
// Keep this conservative until a broader size/distribution matrix is measured.
const AUTO_APPLE_GPU_MIN_LEFT_INTERVALS: usize = 5_000_000;

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

        let default = BioConfig::default();
        let state = self.session.state();
        let config = state.config();
        let bio_config = config
            .options()
            .extensions
            .get::<BioConfig>()
            .unwrap_or(&default);
        let backend = self.select_backend(left_table, bio_config.count_overlaps_backend)?;

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
            backend,
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

#[derive(Clone)]
enum CountOverlapsBackendKind {
    Cpu {
        trees: Arc<AHashMap<String, COITree<(), u32>>>,
    },
    #[cfg(all(feature = "apple-gpu", target_os = "macos"))]
    AppleGpu {
        backend: Arc<AppleGpuCountOverlapsBackend>,
    },
}

impl CountOverlapsBackendKind {
    fn name(&self) -> &'static str {
        match self {
            CountOverlapsBackendKind::Cpu { .. } => "cpu-coitree",
            #[cfg(all(feature = "apple-gpu", target_os = "macos"))]
            CountOverlapsBackendKind::AppleGpu { .. } => "apple-gpu-rank",
        }
    }
}

impl CountOverlapsProvider {
    fn select_backend(
        &self,
        left_table: Vec<RecordBatch>,
        backend_mode: CountOverlapsBackendMode,
    ) -> Result<CountOverlapsBackendKind> {
        let left_interval_count = left_table.iter().map(RecordBatch::num_rows).sum::<usize>();
        let try_apple_gpu = !self.coverage
            && match backend_mode {
                CountOverlapsBackendMode::AppleGpu => true,
                CountOverlapsBackendMode::Auto => {
                    left_interval_count >= AUTO_APPLE_GPU_MIN_LEFT_INTERVALS
                }
                CountOverlapsBackendMode::Cpu => false,
            };

        if try_apple_gpu {
            #[cfg(all(feature = "apple-gpu", target_os = "macos"))]
            {
                let rank_index = CountOverlapsRankIndex::build_from_batches(
                    &left_table,
                    (&self.columns_1.0, &self.columns_1.1, &self.columns_1.2),
                )?;
                match AppleGpuCountOverlapsBackend::try_new(rank_index) {
                    Ok(backend) => {
                        return Ok(CountOverlapsBackendKind::AppleGpu {
                            backend: Arc::new(backend),
                        });
                    }
                    Err(error) => {
                        warn!(
                            "falling back to CPU count_overlaps backend after Apple GPU initialization failed: {error}"
                        );
                    }
                }
            }

            #[cfg(not(all(feature = "apple-gpu", target_os = "macos")))]
            if backend_mode == CountOverlapsBackendMode::AppleGpu {
                warn!(
                    "falling back to CPU count_overlaps backend because apple-gpu feature is not active on macOS"
                );
            }
        }

        let trees = Arc::new(build_coitree_from_batches(
            left_table,
            (&self.columns_1.0, &self.columns_1.1, &self.columns_1.2),
            self.coverage,
        )?);
        Ok(CountOverlapsBackendKind::Cpu { trees })
    }
}

struct CountOverlapsExec {
    schema: SchemaRef,
    backend: CountOverlapsBackendKind,
    right: Arc<dyn ExecutionPlan>,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    coverage: bool,
    cache: PlanProperties,
}

impl Debug for CountOverlapsExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CountOverlapsExec {{ backend: {} }}",
            self.backend.name()
        )
    }
}

impl DisplayAs for CountOverlapsExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter) -> std::fmt::Result {
        write!(f, "CountOverlapsExec: backend={}", self.backend.name())
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
            backend: self.backend.clone(),
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
        match &self.backend {
            CountOverlapsBackendKind::Cpu { trees } => get_stream(
                Arc::clone(&self.right),
                Arc::clone(trees),
                self.schema.clone(),
                Arc::clone(&self.columns_2),
                self.filter_op.clone(),
                self.coverage,
                partition,
                context,
            ),
            #[cfg(all(feature = "apple-gpu", target_os = "macos"))]
            CountOverlapsBackendKind::AppleGpu { backend } => get_apple_gpu_stream(
                Arc::clone(&self.right),
                Arc::clone(backend),
                self.schema.clone(),
                Arc::clone(&self.columns_2),
                self.filter_op.clone(),
                partition,
                context,
            ),
        }
    }
}
