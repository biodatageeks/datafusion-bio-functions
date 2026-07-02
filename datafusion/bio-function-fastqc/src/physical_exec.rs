use std::any::Any;
use std::fmt;
use std::sync::Arc;

use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::TaskContext;
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, Partitioning, PlanProperties,
    SendableRecordBatchStream,
};
use futures::StreamExt;

use crate::{tidy_schema, ModuleSet};

/// Physical operator: fold FASTQ (sequence, quality_scores) batches through the
/// selected QC modules and emit a single tidy RecordBatch.
pub struct FastqcExec {
    input: Arc<dyn ExecutionPlan>,
    selection: Option<Vec<String>>,
    cache: Arc<PlanProperties>,
}

impl FastqcExec {
    pub fn new(input: Arc<dyn ExecutionPlan>, selection: Option<Vec<String>>) -> Self {
        let cache = Arc::new(PlanProperties::new(
            EquivalenceProperties::new(tidy_schema()),
            Partitioning::UnknownPartitioning(1),
            EmissionType::Final,
            Boundedness::Bounded,
        ));
        Self {
            input,
            selection,
            cache,
        }
    }
}

impl fmt::Debug for FastqcExec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FastqcExec(selection={:?})", self.selection)
    }
}

impl DisplayAs for FastqcExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FastqcExec: modules={:?}", self.selection)
    }
}

impl ExecutionPlan for FastqcExec {
    fn name(&self) -> &str {
        "FastqcExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        tidy_schema()
    }

    fn properties(&self) -> &Arc<PlanProperties> {
        &self.cache
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.input]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        Ok(Arc::new(FastqcExec::new(
            children[0].clone(),
            self.selection.clone(),
        )))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        assert_eq!(partition, 0, "FastqcExec has a single output partition");
        let input = self.input.clone();
        let selection = self.selection.clone();
        let n_parts = input.properties().partitioning.partition_count();
        let schema = tidy_schema();

        let fut = async move {
            let mut tasks = Vec::with_capacity(n_parts);
            for p in 0..n_parts {
                let input = input.clone();
                let selection = selection.clone();
                let ctx = context.clone();
                tasks.push(async move {
                    let mut set = ModuleSet::build(selection.as_deref())?;
                    let mut stream = input.execute(p, ctx)?;
                    while let Some(batch) = stream.next().await {
                        set.update_batch(&batch?)?;
                    }
                    Ok::<ModuleSet, DataFusionError>(set)
                });
            }
            let sets = futures::future::try_join_all(tasks).await?;
            let mut merged = ModuleSet::build(selection.as_deref())?;
            for s in sets {
                merged.merge(s);
            }
            let batch = merged.finalize()?;
            Ok::<RecordBatch, DataFusionError>(batch)
        };

        let stream = futures::stream::once(fut);
        Ok(Box::pin(RecordBatchStreamAdapter::new(schema, stream)))
    }
}

#[cfg(test)]
mod tests {
    use datafusion::arrow::array::{Float64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use datafusion::prelude::*;

    use super::*;

    fn read_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("sequence", DataType::Utf8, true),
            Field::new("quality_scores", DataType::Utf8, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ACGT", "GGCC"])),
                Arc::new(StringArray::from(vec!["IIII", "IIII"])),
            ],
        )
        .unwrap()
    }

    async fn run_n_seq(n_parts: usize) -> f64 {
        let schema = read_batch().schema();
        let partitions: Vec<Vec<RecordBatch>> = (0..n_parts).map(|_| vec![read_batch()]).collect();
        let mem_table = MemTable::try_new(schema, partitions).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();
        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let exec = FastqcExec::new(plan, Some(vec!["basic_stats".into()]));
        assert_eq!(exec.properties().partitioning.partition_count(), 1);

        let mut stream = exec.execute(0, ctx.task_ctx()).unwrap();
        let b = stream.next().await.unwrap().unwrap();
        let metric = b
            .column_by_name("metric")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let value = b
            .column_by_name("value")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let idx = (0..b.num_rows())
            .find(|&i| metric.value(i) == "n_seq")
            .unwrap();
        value.value(idx)
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn same_result_1_vs_n_partitions() {
        // 2 reads per partition; merged n_seq must scale with partition count.
        assert_eq!(run_n_seq(1).await, 2.0);
        assert_eq!(run_n_seq(4).await, 8.0);
    }
}
