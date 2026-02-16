use std::any::Any;
use std::collections::{HashMap, VecDeque};
use std::fmt;
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::arrow::datatypes::{SchemaRef, UInt32Type};
use datafusion::common::Result;
use datafusion::execution::TaskContext;
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, PlanProperties, RecordBatchStream,
    SendableRecordBatchStream,
};
use futures::stream::{Stream, StreamExt};

use crate::cigar;
use crate::coverage;
use crate::events;
use crate::events::{ColumnIndices, DenseContigDepth};
use crate::filter::ReadFilter;
use crate::schema::{coverage_output_schema, per_base_output_schema};

/// Handle a contig transition in the dense accumulator.
///
/// If the new contig differs from the current one, finalizes the previous contig
/// (converting its depth array to a RecordBatch and pushing it to `pending_batches`),
/// then allocates a new `DenseContigDepth` for the incoming contig.
///
/// Returns `true` if the row should be processed (contig is available),
/// `false` if the contig is unknown and should be skipped.
#[allow(clippy::too_many_arguments)]
fn handle_contig_transition(
    chrom: &str,
    current_contig: &mut Option<String>,
    current_depth: &mut Option<DenseContigDepth>,
    contig_lengths: &HashMap<String, usize>,
    pending_batches: &mut VecDeque<RecordBatch>,
    schema: &SchemaRef,
    batch_size: usize,
    per_base_emitter: &mut Option<coverage::PerBaseEmitter>,
    per_base: bool,
    zero_based: bool,
) -> bool {
    let contig_changed = current_contig.as_ref().is_none_or(|c| c.as_str() != chrom);
    if !contig_changed {
        return true;
    }

    // Emit the previous contig if we had one
    if let (Some(prev_contig), Some(prev_depth)) = (current_contig.take(), current_depth.take()) {
        if per_base {
            // Flush any existing per-base emitter (rare multi-contig transition)
            if let Some(emitter) = per_base_emitter.as_mut() {
                if let Ok(batches) = emitter.flush_remaining(schema, batch_size) {
                    pending_batches.extend(batches);
                }
            }
            *per_base_emitter = Some(coverage::PerBaseEmitter::new(
                prev_contig,
                prev_depth.depth,
                zero_based,
            ));
        } else if let Ok(batches) =
            coverage::dense_depth_to_record_batches(&prev_contig, &prev_depth, schema, batch_size)
        {
            pending_batches.extend(batches);
        }
    }

    // Start a new contig
    if let Some(&len) = contig_lengths.get(chrom) {
        *current_contig = Some(chrom.to_string());
        *current_depth = Some(DenseContigDepth::new(len));
        true
    } else {
        false
    }
}

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
/// Each partition runs independently — the child plan handles partition-level parallelism.
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
            input.properties().partitioning.clone(),
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
        let input_stream = self.input.execute(partition, context.clone())?;
        let batch_size = context.session_config().batch_size();
        let config = self.config.clone();
        let schema = self.schema();

        Ok(Box::pin(PileupStream::new(
            input_stream,
            config,
            schema,
            batch_size,
        )))
    }
}

/// Depth accumulator strategy — dense (mosdepth-style) or sparse (BTreeMap fallback).
enum DepthAccumulator {
    /// Dense array accumulation with streaming per-contig emission.
    /// Used when BAM header metadata provides contig lengths.
    Dense {
        /// Currently active contig name.
        current_contig: Option<String>,
        /// Dense depth array for the current contig.
        current_depth: Option<DenseContigDepth>,
        /// Contig lengths from BAM header metadata.
        contig_lengths: HashMap<String, usize>,
    },
    /// Sparse BTreeMap accumulation — fallback when no contig lengths available.
    Sparse {
        contig_events: HashMap<String, events::ContigEvents>,
    },
}

impl DepthAccumulator {
    fn is_dense(&self) -> bool {
        matches!(self, DepthAccumulator::Dense { .. })
    }
}

/// Stream that consumes input batches, accumulates coverage events, and emits
/// coverage blocks. Supports two strategies:
///
/// - **Dense** (preferred): uses a flat `i32[]` array per contig, streaming
///   completed contigs as soon as the BAM's coordinate-sorted order moves
///   to the next contig. Peak memory = one contig at a time.
/// - **Sparse** (fallback): uses `BTreeMap` per contig, emitting all contigs
///   at end. Used when schema lacks BAM header metadata (e.g., MemTable tests).
struct PileupStream {
    /// The input stream from the child plan.
    input: SendableRecordBatchStream,
    /// Pileup configuration.
    config: PileupConfig,
    /// Output schema.
    schema: SchemaRef,
    /// Maximum rows per output batch (from DataFusion session config).
    batch_size: usize,
    /// Whether we've finished processing (no more output).
    finished: bool,
    /// Whether we've consumed all input.
    input_exhausted: bool,
    /// Whether we've checked the input schema for contig lengths.
    initialized: bool,
    /// The depth accumulator (dense or sparse).
    accumulator: DepthAccumulator,
    /// Queue of completed RecordBatches ready to emit.
    pending_batches: VecDeque<RecordBatch>,
    /// Cached column indices (computed once on first batch).
    col_idx: Option<ColumnIndices>,
    /// Lazy per-base emitter for the current/last contig (per_base mode only).
    per_base_emitter: Option<coverage::PerBaseEmitter>,
}

impl PileupStream {
    fn new(
        input: SendableRecordBatchStream,
        config: PileupConfig,
        schema: SchemaRef,
        batch_size: usize,
    ) -> Self {
        Self {
            input,
            config,
            schema,
            batch_size,
            finished: false,
            input_exhausted: false,
            initialized: false,
            accumulator: DepthAccumulator::Sparse {
                contig_events: HashMap::new(),
            },
            pending_batches: VecDeque::new(),
            col_idx: None,
            per_base_emitter: None,
        }
    }

    /// Try to initialize the dense accumulator from the input schema metadata.
    /// Also computes and caches column indices from the first batch's schema.
    /// Called once on the first input batch. Respects `DenseMode` config.
    ///
    /// Returns an error if `per_base` mode is enabled but dense initialization
    /// fails (no contig lengths in schema metadata).
    fn try_init_dense(&mut self, input_schema: &SchemaRef) -> Result<()> {
        // Cache column indices once
        self.col_idx = Some(ColumnIndices::from_schema(input_schema));

        let use_dense = match self.config.dense_mode {
            DenseMode::Force => true,
            DenseMode::Disable => false,
            DenseMode::Auto => true, // default to dense for best performance
        };

        if use_dense {
            if let Some(contig_lengths) = events::extract_contig_lengths(input_schema) {
                self.accumulator = DepthAccumulator::Dense {
                    current_contig: None,
                    current_depth: None,
                    contig_lengths,
                };
            }
        }

        if self.config.per_base && !self.accumulator.is_dense() {
            return Err(datafusion::common::DataFusionError::Execution(
                "per_base mode requires dense accumulation (BAM header with contig lengths). \
                 Sparse fallback (e.g. MemTable) is not supported for per_base output."
                    .to_string(),
            ));
        }

        self.initialized = true;
        Ok(())
    }

    /// Process a single batch through the dense accumulator.
    ///
    /// For each row, applies the CIGAR directly to the current contig's
    /// dense depth array. When a contig transition is detected, the completed
    /// contig's coverage is converted to a RecordBatch and queued for emission.
    ///
    /// Splits into separate binary and string CIGAR loops to avoid per-row
    /// `Option::unwrap()` overhead and branch on CIGAR format.
    fn process_batch_dense(&mut self, batch: &RecordBatch) {
        let col_idx = self.col_idx.as_ref().expect("col_idx not initialized");

        let chrom_arr = batch.column(col_idx.chrom).as_string::<i32>();
        let start_arr = batch.column(col_idx.start).as_primitive::<UInt32Type>();
        let flags_arr = batch.column(col_idx.flags).as_primitive::<UInt32Type>();
        let mapq_arr = batch.column(col_idx.mapq).as_primitive::<UInt32Type>();

        let DepthAccumulator::Dense {
            current_contig,
            current_depth,
            contig_lengths,
        } = &mut self.accumulator
        else {
            return;
        };

        let pending_batches = &mut self.pending_batches;
        let schema = &self.schema;
        let batch_size = self.batch_size;
        let per_base_emitter = &mut self.per_base_emitter;
        let per_base = self.config.per_base;
        let zero_based = self.config.zero_based;

        if col_idx.binary_cigar {
            let cigar_arr = batch.column(col_idx.cigar).as_binary::<i32>();
            for row in 0..batch.num_rows() {
                if chrom_arr.is_null(row) || start_arr.is_null(row) {
                    continue;
                }
                let chrom = chrom_arr.value(row);
                let cigar_bytes = cigar_arr.value(row);
                if cigar_bytes.is_empty() {
                    continue;
                }
                let flags = flags_arr.value(row);
                let mapq = mapq_arr.value(row);
                if !self.config.filter.passes(flags, mapq) {
                    continue;
                }
                if !handle_contig_transition(
                    chrom,
                    current_contig,
                    current_depth,
                    contig_lengths,
                    pending_batches,
                    schema,
                    batch_size,
                    per_base_emitter,
                    per_base,
                    zero_based,
                ) {
                    continue;
                }
                let depth = current_depth.as_mut().unwrap();
                let start = start_arr.value(row);
                if let Some((lo, hi)) =
                    cigar::apply_binary_cigar_to_depth(start, cigar_bytes, &mut depth.depth)
                {
                    depth.update_bounds(lo, hi);
                }
            }
        } else {
            let cigar_arr = batch.column(col_idx.cigar).as_string::<i32>();
            for row in 0..batch.num_rows() {
                if chrom_arr.is_null(row) || start_arr.is_null(row) {
                    continue;
                }
                let chrom = chrom_arr.value(row);
                let cigar_str = cigar_arr.value(row);
                if cigar_str == "*" {
                    continue;
                }
                let flags = flags_arr.value(row);
                let mapq = mapq_arr.value(row);
                if !self.config.filter.passes(flags, mapq) {
                    continue;
                }
                if !handle_contig_transition(
                    chrom,
                    current_contig,
                    current_depth,
                    contig_lengths,
                    pending_batches,
                    schema,
                    batch_size,
                    per_base_emitter,
                    per_base,
                    zero_based,
                ) {
                    continue;
                }
                let depth = current_depth.as_mut().unwrap();
                let start = start_arr.value(row);
                if let Some((lo, hi)) =
                    cigar::apply_cigar_to_depth(start, cigar_str, &mut depth.depth)
                {
                    depth.update_bounds(lo, hi);
                }
            }
        }
    }

    /// Finalize the dense accumulator: emit the last contig's coverage.
    fn finalize_dense(&mut self) {
        let DepthAccumulator::Dense {
            current_contig,
            current_depth,
            ..
        } = &mut self.accumulator
        else {
            return;
        };

        if let (Some(contig), Some(depth)) = (current_contig.take(), current_depth.take()) {
            if self.config.per_base {
                // Flush any existing per-base emitter first (rare multi-contig case)
                if let Some(emitter) = self.per_base_emitter.as_mut() {
                    if let Ok(batches) = emitter.flush_remaining(&self.schema, self.batch_size) {
                        self.pending_batches.extend(batches);
                    }
                }
                self.per_base_emitter = Some(coverage::PerBaseEmitter::new(
                    contig,
                    depth.depth,
                    self.config.zero_based,
                ));
            } else if let Ok(batches) = coverage::dense_depth_to_record_batches(
                &contig,
                &depth,
                &self.schema,
                self.batch_size,
            ) {
                self.pending_batches.extend(batches);
            }
        }
    }
}

impl Stream for PileupStream {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let this = self.get_mut();

        loop {
            // 0. Drain per-base emitter lazily (one batch per poll)
            if let Some(emitter) = &mut this.per_base_emitter {
                if let Some(result) = emitter.next_batch(&this.schema, this.batch_size) {
                    return Poll::Ready(Some(result));
                }
                this.per_base_emitter = None;
            }

            // 1. Drain pending batches (shared queue for both dense and sparse)
            if let Some(batch) = this.pending_batches.pop_front() {
                return Poll::Ready(Some(Ok(batch)));
            }

            // 2. If finished, no more output
            if this.finished {
                return Poll::Ready(None);
            }

            // 3. Consume input
            if !this.input_exhausted {
                match this.input.poll_next_unpin(cx) {
                    Poll::Ready(Some(Ok(batch))) => {
                        if !this.initialized {
                            if let Err(e) = this.try_init_dense(&batch.schema()) {
                                this.finished = true;
                                return Poll::Ready(Some(Err(e)));
                            }
                        }

                        // Use is_dense() to avoid holding a borrow on this.accumulator
                        // when calling methods that need &mut self
                        if this.accumulator.is_dense() {
                            this.process_batch_dense(&batch);
                            // Loop back to drain any pending batches from contig transitions
                            continue;
                        }

                        // Sparse path
                        if let DepthAccumulator::Sparse { contig_events } = &mut this.accumulator {
                            let col_idx = this.col_idx.as_ref().expect("col_idx not initialized");
                            events::process_batch(
                                &batch,
                                &this.config.filter,
                                contig_events,
                                col_idx,
                            );
                        }
                        continue;
                    }
                    Poll::Ready(Some(Err(e))) => {
                        this.finished = true;
                        return Poll::Ready(Some(Err(e)));
                    }
                    Poll::Ready(None) => {
                        this.input_exhausted = true;
                        // Fall through to finalization
                    }
                    Poll::Pending => {
                        return Poll::Pending;
                    }
                }
            }

            // 4. Input exhausted — finalize
            this.finished = true;

            if this.accumulator.is_dense() {
                this.finalize_dense();
                // Loop back to drain pending batches / per-base emitter
                continue;
            }

            // Sparse finalization: emit chunked batches
            if let DepthAccumulator::Sparse { contig_events } = &mut this.accumulator {
                if contig_events.is_empty() {
                    return Poll::Ready(None);
                }
                match coverage::all_events_to_record_batches(
                    contig_events,
                    &this.schema,
                    this.batch_size,
                ) {
                    Ok(batches) => this.pending_batches.extend(batches),
                    Err(e) => return Poll::Ready(Some(Err(e))),
                }
                // Loop back to drain pending batches (step 1), then step 2 returns None
                continue;
            }

            return Poll::Ready(None);
        }
    }
}

impl RecordBatchStream for PileupStream {
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

    fn make_batch_binary(rows: Vec<(&str, u32, u32, u32, Vec<u8>, u32)>) -> RecordBatch {
        let schema = bam_schema_binary();
        let chroms: Vec<Option<&str>> = rows.iter().map(|r| Some(r.0.as_ref())).collect();
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

    #[tokio::test]
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

    #[tokio::test]
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

    #[tokio::test]
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

    #[tokio::test]
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
        let task_ctx = ctx.task_ctx();

        // Partition 0
        let stream0 = pileup.execute(0, task_ctx.clone()).unwrap();
        let batches0: Vec<RecordBatch> = stream0
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();
        assert_eq!(batches0.len(), 1);
        let contigs0 = batches0[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(contigs0.value(0), "chr1");

        // Partition 1
        let stream1 = pileup.execute(1, task_ctx).unwrap();
        let batches1: Vec<RecordBatch> = stream1
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();
        assert_eq!(batches1.len(), 1);
        let contigs1 = batches1[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(contigs1.value(0), "chr2");
    }

    #[tokio::test]
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

    #[tokio::test]
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

    #[tokio::test]
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
    #[tokio::test]
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
