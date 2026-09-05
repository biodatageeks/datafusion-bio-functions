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
use tokio::task::JoinSet;

use crate::{ModuleSet, tidy_schema};

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
    fn apply_expressions(
        &self,
        _f: &mut dyn FnMut(
            &std::sync::Arc<dyn datafusion::physical_expr::PhysicalExpr>,
        ) -> datafusion::common::Result<
            datafusion::common::tree_node::TreeNodeRecursion,
        >,
    ) -> datafusion::common::Result<datafusion::common::tree_node::TreeNodeRecursion> {
        Ok(datafusion::common::tree_node::TreeNodeRecursion::Continue)
    }

    fn name(&self) -> &str {
        "FastqcExec"
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

        // Some modules are order-dependent and only match FastQC when reads are
        // seen in file order:
        //   - kmer_content samples every 50th read in file order.
        //   - overrepresented and dup_levels stop tracking NEW distinct keys once
        //     they hit FastQC's OBSERVATION_CUTOFF (100k uniques); which keys
        //     survive the freeze depends on the global read order, so applying the
        //     cutoff per partition can track a different key set than FastQC.
        // Accumulating one ModuleSet per partition and merging would process each
        // partition independently (partition-dependent, non-FastQC). When any
        // order-dependent module is selected we therefore read partitions
        // sequentially in index order (= file order for range-based scans) into a
        // single ModuleSet.
        let order_dependent = match &selection {
            None => true, // None resolves to all modules, which includes them
            Some(sel) => sel.iter().any(|m| {
                matches!(
                    m.as_str(),
                    "kmer_content" | "overrepresented" | "dup_levels"
                )
            }),
        };

        let fut = async move {
            if order_dependent {
                let mut set = ModuleSet::build(selection.as_deref())?;
                for p in 0..n_parts {
                    let mut stream = input.execute(p, context.clone())?;
                    while let Some(batch) = stream.next().await {
                        set.update_batch(&batch?)?;
                    }
                }
                let batch = set.finalize()?;
                return Ok::<RecordBatch, DataFusionError>(batch);
            }

            // No order-dependent module selected: spawn one task per input
            // partition so accumulation runs in parallel across runtime worker threads
            // (CPU-bound work). A plain `try_join_all` would only interleave the
            // futures on a single task, adding overhead without parallelism.
            let mut join_set: JoinSet<Result<ModuleSet>> = JoinSet::new();
            for p in 0..n_parts {
                let input = input.clone();
                let selection = selection.clone();
                let ctx = context.clone();
                join_set.spawn(async move {
                    let mut set = ModuleSet::build(selection.as_deref())?;
                    let mut stream = input.execute(p, ctx)?;
                    while let Some(batch) = stream.next().await {
                        set.update_batch(&batch?)?;
                    }
                    Ok::<ModuleSet, DataFusionError>(set)
                });
            }
            let mut merged = ModuleSet::build(selection.as_deref())?;
            while let Some(joined) = join_set.join_next().await {
                let set = joined
                    .map_err(|e| DataFusionError::Execution(format!("fastqc join error: {e}")))??;
                merged.merge(set);
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
    use datafusion::arrow::array::{Array, Float64Array, StringArray};
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

    /// Collect kmer_content (label -> count) running `total` motif-bearing reads
    /// through FastqcExec split across `n_parts` file-ordered partitions.
    async fn run_kmer(total: usize, n_parts: usize) -> Vec<(String, i64)> {
        let schema = Arc::new(Schema::new(vec![
            Field::new("sequence", DataType::Utf8, true),
            Field::new("quality_scores", DataType::Utf8, true),
        ]));
        // A fixed motif at a fixed position in every read -> strong kmer
        // enrichment; the 2% (every-50th) sample must be identical whether read
        // as one partition or several contiguous ones.
        let seqs: Vec<String> = (0..total)
            .map(|i| format!("ACGTACGTACGATTACGTTTT{:09}", i % 10))
            .collect();
        let per = total.div_ceil(n_parts);
        let partitions: Vec<Vec<RecordBatch>> = (0..n_parts)
            .map(|p| {
                let lo = p * per;
                let hi = ((p + 1) * per).min(total);
                let s: Vec<&str> = seqs[lo..hi].iter().map(|s| s.as_str()).collect();
                let q: Vec<&str> = s.iter().map(|_| "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII").collect();
                vec![
                    RecordBatch::try_new(
                        schema.clone(),
                        vec![
                            Arc::new(StringArray::from(s)),
                            Arc::new(StringArray::from(q)),
                        ],
                    )
                    .unwrap(),
                ]
            })
            .collect();

        let mem = MemTable::try_new(schema, partitions).unwrap();
        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem)).unwrap();
        let plan = ctx
            .table("reads")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();
        let exec = FastqcExec::new(plan, Some(vec!["kmer_content".into()]));
        let mut stream = exec.execute(0, ctx.task_ctx()).unwrap();
        let b = stream.next().await.unwrap().unwrap();
        let label = b
            .column_by_name("label")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
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
        let mut out: Vec<(String, i64)> = (0..b.num_rows())
            .filter(|&i| metric.value(i) == "count" && !label.is_null(i))
            .map(|i| (label.value(i).to_string(), value.value(i) as i64))
            .collect();
        out.sort();
        out
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn kmer_content_is_partition_invariant_when_ordered() {
        // kmer must sample in file order: splitting the same reads across
        // contiguous partitions must not change the reported kmers/counts.
        let one = run_kmer(300, 1).await;
        assert!(!one.is_empty(), "fixture should enrich some kmers");
        assert_eq!(one, run_kmer(300, 3).await);
        assert_eq!(one, run_kmer(300, 7).await);
    }
}
