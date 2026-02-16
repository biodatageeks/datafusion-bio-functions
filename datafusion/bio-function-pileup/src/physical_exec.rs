use std::any::Any;
use std::collections::{HashMap, VecDeque};
use std::fmt;
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::TaskContext;
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, Partitioning, PlanProperties, RecordBatchStream,
    SendableRecordBatchStream,
};
use futures::stream::{Stream, StreamExt};
use tokio::task::JoinSet;

use crate::coverage;
use crate::events;
use crate::events::{ColumnIndices, DenseContigDepth};
use crate::filter::ReadFilter;
use crate::schema::{coverage_output_schema, per_base_output_schema};

/// Controls which depth-accumulation strategy the pileup stream uses.
///
/// - `Auto` — heuristic (currently defaults to sparse; may be enhanced later).
/// - `Force` — always use the dense array when BAM header metadata is available.
/// - `Disable` — always use the sparse BTreeMap, never attempt dense.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum DenseMode {
    /// Heuristic selection (currently equivalent to `Force`).
    #[default]
    Auto,
    /// Force dense array accumulation (best for WGS with high coverage).
    Force,
    /// Always use sparse BTreeMap (best for WES / targeted panels).
    Disable,
}

/// Configuration for the pileup computation.
#[derive(Debug, Clone)]
pub struct PileupConfig {
    pub filter: ReadFilter,
    pub dense_mode: DenseMode,
    /// When `true`, CIGAR data arrives as packed binary (`DataType::Binary`) instead of strings.
    /// The actual format is detected from the schema, but this flag controls what we request
    /// from the upstream `BamTableProvider`.
    pub binary_cigar: bool,
    /// When `true`, output coordinates are 0-based (start inclusive, end inclusive).
    /// When `false` (default), output coordinates are 1-based.
    /// This also controls the coordinate system passed to the BAM reader.
    pub zero_based: bool,
    /// When `true`, emit one row per genomic position (like `samtools depth -a`)
    /// instead of RLE coverage blocks. Requires dense mode (BAM header with
    /// contig lengths). Default: `false`.
    pub per_base: bool,
}

impl Default for PileupConfig {
    fn default() -> Self {
        Self {
            filter: ReadFilter::default(),
            dense_mode: DenseMode::default(),
            binary_cigar: true,
            zero_based: false,
            per_base: false,
        }
    }
}

/// DataFusion ExecutionPlan that computes depth-of-coverage from BAM/alignment data.
///
/// Wraps a child plan (the BAM reader) and computes per-contig coverage blocks.
/// Always produces exactly 1 output partition by merging all input partitions.
/// Phase 1: N input partitions are drained in parallel (one tokio task each).
/// Phase 2: Per-contig delta arrays are merged element-wise, then emitted as
/// coverage blocks or per-base rows.
#[derive(Debug)]
pub struct PileupExec {
    /// The child execution plan (BAM reader).
    input: Arc<dyn ExecutionPlan>,
    /// Pileup configuration.
    config: PileupConfig,
    /// Cached plan properties.
    cache: PlanProperties,
}

impl PileupExec {
    pub fn new(input: Arc<dyn ExecutionPlan>, config: PileupConfig) -> Self {
        let schema = if config.per_base {
            per_base_output_schema(config.zero_based)
        } else {
            coverage_output_schema(config.zero_based)
        };
        let cache = PlanProperties::new(
            EquivalenceProperties::new(schema),
            Partitioning::UnknownPartitioning(1), // always 1 output partition
            EmissionType::Incremental,
            Boundedness::Bounded,
        );
        Self {
            input,
            config,
            cache,
        }
    }
}

impl DisplayAs for PileupExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "PileupExec: filter_flag={}, min_mapq={}, dense={:?}, binary_cigar={}, zero_based={}, per_base={}",
            self.config.filter.filter_flag,
            self.config.filter.min_mapping_quality,
            self.config.dense_mode,
            self.config.binary_cigar,
            self.config.zero_based,
            self.config.per_base
        )
    }
}

impl ExecutionPlan for PileupExec {
    fn name(&self) -> &str {
        "PileupExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.cache.eq_properties.schema().clone()
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
        assert_eq!(children.len(), 1);
        Ok(Arc::new(PileupExec::new(
            children[0].clone(),
            self.config.clone(),
        )))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        assert_eq!(
            partition, 0,
            "PileupExec always produces 1 output partition"
        );

        let num_input_partitions = self.input.properties().partitioning.partition_count();
        let mut input_streams = Vec::with_capacity(num_input_partitions);
        for p in 0..num_input_partitions {
            input_streams.push(self.input.execute(p, context.clone())?);
        }

        let batch_size = context.session_config().batch_size();
        let config = self.config.clone();
        let schema = self.schema();

        Ok(Box::pin(MergingPileupStream::new(
            input_streams,
            config,
            schema,
            batch_size,
        )))
    }
}

// ---------------------------------------------------------------------------
// Merging pileup stream — replaces the old per-partition PileupStream
// ---------------------------------------------------------------------------

/// Result of draining one input partition.
enum PartitionResult {
    /// Dense delta arrays accumulated per contig (BAM header metadata available).
    Dense {
        contig_depths: HashMap<String, DenseContigDepth>,
    },
    /// Sparse event lists per contig (fallback when no contig lengths).
    Sparse {
        contig_events: HashMap<String, events::ContigEvents>,
    },
}

/// Internal state machine for the merging stream.
enum StreamPhase {
    /// Phase 1: parallel accumulation of input partitions.
    Accumulating {
        join_set: JoinSet<Result<PartitionResult>>,
        partial_results: Vec<PartitionResult>,
        total: usize,
    },
    /// Phase 2: merged results being emitted as RecordBatches.
    Emitting {
        pending_batches: VecDeque<RecordBatch>,
        per_base_emitter: Option<coverage::PerBaseEmitter>,
        remaining_contigs: VecDeque<(String, Vec<i32>)>,
    },
    /// Terminal state — no more output.
    Done,
}

/// Stream that merges coverage from all input partitions before emitting.
///
/// Phase 1 spawns one tokio task per input partition; each task drains its
/// stream via [`accumulate_partition`]. Phase 2 element-wise adds delta
/// arrays per contig across all partitions, then converts merged arrays
/// to coverage blocks or per-base rows.
struct MergingPileupStream {
    /// Output schema.
    schema: SchemaRef,
    /// Pileup configuration.
    config: PileupConfig,
    /// Maximum rows per output batch.
    batch_size: usize,
    /// Current phase of the stream.
    phase: StreamPhase,
}

impl MergingPileupStream {
    fn new(
        input_streams: Vec<SendableRecordBatchStream>,
        config: PileupConfig,
        schema: SchemaRef,
        batch_size: usize,
    ) -> Self {
        let total = input_streams.len();
        let mut join_set = JoinSet::new();

        for stream in input_streams {
            let cfg = config.clone();
            join_set.spawn(accumulate_partition(stream, cfg));
        }

        Self {
            schema,
            config,
            batch_size,
            phase: StreamPhase::Accumulating {
                join_set,
                partial_results: Vec::with_capacity(total),
                total,
            },
        }
    }
}

/// Drain one input partition stream, accumulating CIGAR events into
/// dense delta arrays (preferred) or sparse event lists (fallback).
async fn accumulate_partition(
    mut stream: SendableRecordBatchStream,
    config: PileupConfig,
) -> Result<PartitionResult> {
    let mut initialized = false;
    let mut contig_depths: HashMap<String, DenseContigDepth> = HashMap::new();
    let mut contig_events: HashMap<String, events::ContigEvents> = HashMap::new();
    let mut contig_lengths: HashMap<String, usize> = HashMap::new();
    let mut col_idx: Option<ColumnIndices> = None;
    let mut is_dense = false;

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }

        if !initialized {
            col_idx = Some(ColumnIndices::from_schema(&batch.schema()));

            let use_dense = match config.dense_mode {
                DenseMode::Force | DenseMode::Auto => true,
                DenseMode::Disable => false,
            };

            if use_dense {
                if let Some(lengths) = events::extract_contig_lengths(&batch.schema()) {
                    contig_lengths = lengths;
                    is_dense = true;
                }
            }

            if config.per_base && !is_dense {
                return Err(datafusion::common::DataFusionError::Execution(
                    "per_base mode requires dense accumulation (BAM header with contig lengths). \
                     Sparse fallback (e.g. MemTable) is not supported for per_base output."
                        .to_string(),
                ));
            }

            initialized = true;
        }

        let idx = col_idx.as_ref().unwrap();
        if is_dense {
            events::process_batch_dense(
                &batch,
                &config.filter,
                &mut contig_depths,
                &contig_lengths,
                idx,
            );
        } else {
            events::process_batch(&batch, &config.filter, &mut contig_events, idx);
        }
    }

    if is_dense {
        Ok(PartitionResult::Dense { contig_depths })
    } else {
        Ok(PartitionResult::Sparse { contig_events })
    }
}

/// Merge partition results and transition to the emitting phase.
fn merge_results(
    partial_results: Vec<PartitionResult>,
    schema: &SchemaRef,
    batch_size: usize,
    config: &PileupConfig,
) -> StreamPhase {
    if partial_results.is_empty() {
        return StreamPhase::Done;
    }

    // Separate dense and sparse results (empty sparse partitions are discarded)
    let mut dense_results = Vec::new();
    let mut sparse_results = Vec::new();

    for result in partial_results {
        match result {
            PartitionResult::Dense { .. } => dense_results.push(result),
            PartitionResult::Sparse { contig_events } => {
                if !contig_events.is_empty() {
                    sparse_results.push(PartitionResult::Sparse { contig_events });
                }
            }
        }
    }

    if !dense_results.is_empty() {
        merge_dense_results(dense_results, schema, batch_size, config)
    } else if !sparse_results.is_empty() {
        merge_sparse_results(sparse_results, schema, batch_size)
    } else {
        StreamPhase::Done
    }
}

/// Merge dense partition results by element-wise adding delta arrays per contig.
fn merge_dense_results(
    partitions: Vec<PartitionResult>,
    schema: &SchemaRef,
    batch_size: usize,
    config: &PileupConfig,
) -> StreamPhase {
    let mut merged_depths: HashMap<String, DenseContigDepth> = HashMap::new();

    for part in partitions {
        if let PartitionResult::Dense { contig_depths } = part {
            for (contig, depth) in contig_depths {
                let start = depth.touched_start();
                let range = depth.touched_range();
                if range.is_empty() {
                    continue;
                }
                let end = start + range.len() - 1;

                if let Some(existing) = merged_depths.get_mut(&contig) {
                    // Element-wise add the touched region
                    for (i, &delta) in range.iter().enumerate() {
                        if delta != 0 {
                            existing.depth[start + i] += delta;
                        }
                    }
                    existing.update_bounds(start, end);
                } else {
                    merged_depths.insert(contig, depth);
                }
            }
        }
    }

    // Sort contigs for deterministic output
    let mut contigs: Vec<String> = merged_depths.keys().cloned().collect();
    contigs.sort();

    if config.per_base {
        let remaining: VecDeque<(String, Vec<i32>)> = contigs
            .into_iter()
            .filter_map(|c| merged_depths.remove(&c).map(|d| (c, d.depth)))
            .collect();

        StreamPhase::Emitting {
            pending_batches: VecDeque::new(),
            per_base_emitter: None,
            remaining_contigs: remaining,
        }
    } else {
        let mut pending_batches = VecDeque::new();
        for contig in &contigs {
            if let Some(depth) = merged_depths.get(contig) {
                if let Ok(batches) =
                    coverage::dense_depth_to_record_batches(contig, depth, schema, batch_size)
                {
                    pending_batches.extend(batches);
                }
            }
        }
        StreamPhase::Emitting {
            pending_batches,
            per_base_emitter: None,
            remaining_contigs: VecDeque::new(),
        }
    }
}

/// Merge sparse partition results by concatenating event vectors per contig.
fn merge_sparse_results(
    partitions: Vec<PartitionResult>,
    schema: &SchemaRef,
    batch_size: usize,
) -> StreamPhase {
    let mut merged: HashMap<String, events::ContigEvents> = HashMap::new();

    for part in partitions {
        if let PartitionResult::Sparse { contig_events } = part {
            for (contig, ce) in contig_events {
                merged.entry(contig).or_default().events.extend(ce.events);
            }
        }
    }

    if merged.is_empty() {
        return StreamPhase::Done;
    }

    match coverage::all_events_to_record_batches(&mut merged, schema, batch_size) {
        Ok(batches) => StreamPhase::Emitting {
            pending_batches: VecDeque::from(batches),
            per_base_emitter: None,
            remaining_contigs: VecDeque::new(),
        },
        Err(_) => StreamPhase::Done,
    }
}

impl Stream for MergingPileupStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            match &mut this.phase {
                StreamPhase::Accumulating {
                    join_set,
                    partial_results,
                    total,
                } => {
                    match join_set.poll_join_next(cx) {
                        Poll::Ready(Some(Ok(Ok(result)))) => {
                            partial_results.push(result);
                            if partial_results.len() == *total {
                                let results = std::mem::take(partial_results);
                                this.phase = merge_results(
                                    results,
                                    &this.schema,
                                    this.batch_size,
                                    &this.config,
                                );
                            }
                            continue;
                        }
                        Poll::Ready(Some(Ok(Err(e)))) => {
                            // Partition returned a DataFusion error
                            this.phase = StreamPhase::Done;
                            return Poll::Ready(Some(Err(e)));
                        }
                        Poll::Ready(Some(Err(join_err))) => {
                            // Tokio JoinError (panic or cancellation)
                            this.phase = StreamPhase::Done;
                            return Poll::Ready(Some(Err(
                                datafusion::common::DataFusionError::Execution(format!(
                                    "Partition task failed: {join_err}"
                                )),
                            )));
                        }
                        Poll::Ready(None) => {
                            // All tasks completed (JoinSet drained)
                            let results = std::mem::take(partial_results);
                            this.phase =
                                merge_results(results, &this.schema, this.batch_size, &this.config);
                            continue;
                        }
                        Poll::Pending => return Poll::Pending,
                    }
                }

                StreamPhase::Emitting {
                    pending_batches,
                    per_base_emitter,
                    remaining_contigs,
                } => {
                    // 0. Drain per-base emitter lazily (one batch per poll)
                    if let Some(emitter) = per_base_emitter {
                        if let Some(result) = emitter.next_batch(&this.schema, this.batch_size) {
                            return Poll::Ready(Some(result));
                        }
                        *per_base_emitter = None;
                    }

                    // 1. Drain pending batches
                    if let Some(batch) = pending_batches.pop_front() {
                        return Poll::Ready(Some(Ok(batch)));
                    }

                    // 2. Start next contig for per-base mode
                    if let Some((contig, depth_vec)) = remaining_contigs.pop_front() {
                        *per_base_emitter = Some(coverage::PerBaseEmitter::new(
                            contig,
                            depth_vec,
                            this.config.zero_based,
                        ));
                        continue;
                    }

                    // 3. All done
                    this.phase = StreamPhase::Done;
                    return Poll::Ready(None);
                }

                StreamPhase::Done => return Poll::Ready(None),
            }
        }
    }
}

impl RecordBatchStream for MergingPileupStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{BinaryArray, Int16Array, Int32Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use datafusion::prelude::*;

    fn bam_schema() -> SchemaRef {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("start", DataType::UInt32, true),
            Field::new("end", DataType::UInt32, true),
            Field::new("flags", DataType::UInt32, false),
            Field::new("cigar", DataType::Utf8, false),
            Field::new("mapping_quality", DataType::UInt32, false),
        ]))
    }

    fn bam_schema_binary() -> SchemaRef {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("start", DataType::UInt32, true),
            Field::new("end", DataType::UInt32, true),
            Field::new("flags", DataType::UInt32, false),
            Field::new("cigar", DataType::Binary, false),
            Field::new("mapping_quality", DataType::UInt32, false),
        ]))
    }

    fn encode_op(len: u32, code: u32) -> [u8; 4] {
        ((len << 4) | code).to_le_bytes()
    }

    fn make_batch(rows: Vec<(&str, u32, u32, u32, &str, u32)>) -> RecordBatch {
        let schema = bam_schema();
        let chroms: Vec<Option<&str>> = rows.iter().map(|r| Some(r.0)).collect();
        let starts: Vec<Option<u32>> = rows.iter().map(|r| Some(r.1)).collect();
        let ends: Vec<Option<u32>> = rows.iter().map(|r| Some(r.2)).collect();
        let flags: Vec<u32> = rows.iter().map(|r| r.3).collect();
        let cigars: Vec<&str> = rows.iter().map(|r| r.4).collect();
        let mapqs: Vec<u32> = rows.iter().map(|r| r.5).collect();

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(UInt32Array::from(starts)),
                Arc::new(UInt32Array::from(ends)),
                Arc::new(UInt32Array::from(flags)),
                Arc::new(StringArray::from(cigars)),
                Arc::new(UInt32Array::from(mapqs)),
            ],
        )
        .unwrap()
    }

    #[allow(clippy::type_complexity)]
    fn make_batch_binary(rows: Vec<(&str, u32, u32, u32, Vec<u8>, u32)>) -> RecordBatch {
        let schema = bam_schema_binary();
        let chroms: Vec<Option<&str>> = rows.iter().map(|r| Some(r.0)).collect();
        let starts: Vec<Option<u32>> = rows.iter().map(|r| Some(r.1)).collect();
        let ends: Vec<Option<u32>> = rows.iter().map(|r| Some(r.2)).collect();
        let flags: Vec<u32> = rows.iter().map(|r| r.3).collect();
        let cigar_refs: Vec<&[u8]> = rows.iter().map(|r| r.4.as_slice()).collect();
        let mapqs: Vec<u32> = rows.iter().map(|r| r.5).collect();

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(UInt32Array::from(starts)),
                Arc::new(UInt32Array::from(ends)),
                Arc::new(UInt32Array::from(flags)),
                Arc::new(BinaryArray::from(cigar_refs)),
                Arc::new(UInt32Array::from(mapqs)),
            ],
        )
        .unwrap()
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_basic_coverage() {
        let batch = make_batch(vec![("chr1", 0, 10, 0, "10M", 60)]);
        let schema = bam_schema();
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        assert_eq!(batch.num_rows(), 1);

        let contigs = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let starts = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let ends = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let covs = batch
            .column(3)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        assert_eq!(contigs.value(0), "chr1");
        assert_eq!(starts.value(0), 0);
        assert_eq!(ends.value(0), 9);
        assert_eq!(covs.value(0), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_filtering() {
        let batch = make_batch(vec![
            ("chr1", 0, 10, 0, "10M", 60),
            ("chr1", 5, 15, 4, "10M", 60), // unmapped — should be filtered
        ]);
        let schema = bam_schema();
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        assert_eq!(batch.num_rows(), 1);
        let covs = batch
            .column(3)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();
        assert_eq!(covs.value(0), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_empty_input() {
        let schema = bam_schema();
        let batch = RecordBatch::new_empty(schema.clone());
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        assert!(batches.is_empty());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_multi_partition() {
        let batch1 = make_batch(vec![("chr1", 0, 10, 0, "10M", 60)]);
        let batch2 = make_batch(vec![("chr2", 100, 120, 0, "20M", 60)]);
        let schema = bam_schema();
        let mem_table = MemTable::try_new(schema, vec![vec![batch1], vec![batch2]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );

        // Should always produce 1 output partition
        assert_eq!(pileup.properties().partitioning.partition_count(), 1);

        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        // Both contigs should be in the merged output
        let all_contigs: Vec<String> = batches
            .iter()
            .flat_map(|batch| {
                let contigs = batch
                    .column(0)
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .unwrap();
                (0..batch.num_rows()).map(move |i| contigs.value(i).to_string())
            })
            .collect();

        assert!(
            all_contigs.contains(&"chr1".to_string()),
            "merged output should contain chr1"
        );
        assert!(
            all_contigs.contains(&"chr2".to_string()),
            "merged output should contain chr2"
        );
    }

    /// Same contig in 2 partitions — verify merged coverage is correct.
    #[tokio::test(flavor = "multi_thread")]
    async fn test_multi_partition_merge_overlapping() {
        // Partition 0: chr1 read at [0, 10)
        // Partition 1: chr1 read at [5, 15)
        // Merged: chr1 [0,4] cov=1, [5,9] cov=2, [10,14] cov=1
        let batch1 = make_batch(vec![("chr1", 0, 10, 0, "10M", 60)]);
        let batch2 = make_batch(vec![("chr1", 5, 15, 0, "10M", 60)]);
        let schema = bam_schema();
        let mem_table = MemTable::try_new(schema, vec![vec![batch1], vec![batch2]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 3, "Expected 3 coverage blocks");

        let starts: Vec<i32> = batches
            .iter()
            .flat_map(|b| {
                let arr = b.column(1).as_any().downcast_ref::<Int32Array>().unwrap();
                (0..b.num_rows()).map(move |i| arr.value(i))
            })
            .collect();
        let ends: Vec<i32> = batches
            .iter()
            .flat_map(|b| {
                let arr = b.column(2).as_any().downcast_ref::<Int32Array>().unwrap();
                (0..b.num_rows()).map(move |i| arr.value(i))
            })
            .collect();
        let covs: Vec<i16> = batches
            .iter()
            .flat_map(|b| {
                let arr = b.column(3).as_any().downcast_ref::<Int16Array>().unwrap();
                (0..b.num_rows()).map(move |i| arr.value(i))
            })
            .collect();

        assert_eq!(starts, vec![0, 5, 10]);
        assert_eq!(ends, vec![4, 9, 14]);
        assert_eq!(covs, vec![1, 2, 1]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_basic_coverage_binary_cigar() {
        let cigar_10m = encode_op(10, 0); // 10M
        let batch = make_batch_binary(vec![("chr1", 0, 10, 0, cigar_10m.to_vec(), 60)]);
        let schema = bam_schema_binary();
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        assert_eq!(batch.num_rows(), 1);

        let contigs = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let starts = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let ends = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let covs = batch
            .column(3)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        assert_eq!(contigs.value(0), "chr1");
        assert_eq!(starts.value(0), 0);
        assert_eq!(ends.value(0), 9);
        assert_eq!(covs.value(0), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_binary_cigar_filtering() {
        let cigar_10m = encode_op(10, 0); // 10M
        let batch = make_batch_binary(vec![
            ("chr1", 0, 10, 0, cigar_10m.to_vec(), 60),
            ("chr1", 5, 15, 4, cigar_10m.to_vec(), 60), // unmapped — filtered
        ]);
        let schema = bam_schema_binary();
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        assert_eq!(batch.num_rows(), 1);
        let covs = batch
            .column(3)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();
        assert_eq!(covs.value(0), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_batch_size_chunking() {
        // Two overlapping reads produce 3 coverage blocks:
        // [0,4] cov=1, [5,9] cov=2, [10,14] cov=1
        let batch = make_batch(vec![
            ("chr1", 0, 10, 0, "10M", 60),
            ("chr1", 5, 15, 0, "10M", 60),
        ]);
        let schema = bam_schema();
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        // Configure batch_size=2 so 3 rows get split into 2 batches
        let config = SessionConfig::new().with_batch_size(2);
        let ctx = SessionContext::new_with_config(config);
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        // Should produce 2 batches: [2 rows] + [1 row]
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].num_rows(), 2);
        assert_eq!(batches[1].num_rows(), 1);

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 3);
    }

    /// Per-base mode on sparse path (MemTable without contig lengths) must error.
    #[tokio::test(flavor = "multi_thread")]
    async fn test_per_base_sparse_errors() {
        let batch = make_batch(vec![("chr1", 0, 10, 0, "10M", 60)]);
        let schema = bam_schema();
        let mem_table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table("reads", Arc::new(mem_table)).unwrap();

        let df = ctx.table("reads").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let pileup = PileupExec::new(
            plan,
            PileupConfig {
                zero_based: true,
                per_base: true,
                ..PileupConfig::default()
            },
        );
        let task_ctx = ctx.task_ctx();
        let stream = pileup.execute(0, task_ctx).unwrap();
        let results: Vec<Result<RecordBatch>> = stream.collect::<Vec<_>>().await;

        // Should get an error since MemTable lacks contig length metadata
        assert!(results.iter().any(|r| r.is_err()));
        let err = results
            .into_iter()
            .find(|r| r.is_err())
            .unwrap()
            .unwrap_err();
        assert!(
            err.to_string().contains("per_base"),
            "Error should mention per_base: {err}"
        );
    }
}
