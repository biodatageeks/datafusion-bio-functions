//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! the Parquet columnar variation dataset per-position for annotation.

use std::any::Any;
use std::collections::{HashMap, VecDeque};
use std::fmt::{Debug, Formatter};
use std::ops::Range;
use std::path::PathBuf;
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};
use std::time::{Duration, Instant};

use datafusion::arrow::array::{
    Array, ArrayRef, Int8Array, Int16Array, Int32Array, Int64Array, LargeStringArray, RecordBatch,
    StringArray, StringViewArray, UInt8Array, UInt16Array, UInt32Array, UInt64Array,
};
use datafusion::arrow::datatypes::{DataType, Field, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use futures::{Stream, StreamExt};
use smallvec::SmallVec;

use crate::allele::{
    VariantAlleleInput, allele_matches, get_matched_variant_alleles, vcf_to_vep_allele,
    vcf_to_vep_input_allele, vep_norm_end, vep_norm_start,
};
use crate::cache::key_encoding::chrom_to_code;
use crate::cache::variant_key::{
    position_key_from_code as warm_position_key_from_code,
    variant_key_from_position as warm_variant_key_from_position,
};
use crate::cache_common::AlleleMatcher;
use crate::colocated::{
    AF_COL_NAMES, AfColumns, ColocatedCacheEntry, ColocatedKey, ColocatedSink, ColocatedSinkValue,
    compare_existing_variant_alleles, output_allele_from_allele_string, read_reference_sequence,
};
use crate::parquet_cache::variation_lookup::SinglePathParquetVariationLookup;
use datafusion::arrow::array::BooleanArray;
use tokio::sync::OnceCell;

/// Drive a Parquet future to completion from a synchronous context regardless of
/// the ambient Tokio runtime flavor. `block_in_place` panics on a current-thread
/// runtime (`#[tokio::test]`'s default flavor, embedded single-thread callers),
/// so fall back to a dedicated OS thread there; with no runtime in scope, create
/// one. Keeps the multi-thread fast path (`block_in_place`) for production.
#[cfg(feature = "parquet-cache")]
pub(crate) fn block_on<T, F>(fut: F) -> Result<T>
where
    F: std::future::Future<Output = Result<T>> + Send,
    T: Send,
{
    match tokio::runtime::Handle::try_current() {
        Ok(handle) if handle.runtime_flavor() == tokio::runtime::RuntimeFlavor::MultiThread => {
            tokio::task::block_in_place(|| handle.block_on(fut))
        }
        Ok(_handle) => std::thread::scope(|scope| {
            scope
                .spawn(|| {
                    let rt = tokio::runtime::Runtime::new()
                        .map_err(|e| DataFusionError::External(Box::new(e)))?;
                    rt.block_on(fut)
                })
                .join()
                .unwrap_or_else(|_| {
                    Err(DataFusionError::Execution(
                        "cache lookup worker thread panicked".to_string(),
                    ))
                })
        }),
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            rt.block_on(fut)
        }
    }
}

const DEFAULT_LANCE_LOOKUP_PROCESS_BATCH_ROWS: usize = 5_000;

/// Lookup match mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KvMatchMode {
    /// Exact allele matching only.
    Exact,
}

/// Physical execution plan for KV-backed variant lookup.
///
/// Takes a VCF input plan, probes the Parquet variation dataset per-position,
/// and emits LEFT JOIN output (unmatched VCF rows get NULL cache columns).
pub struct KvLookupExec {
    input: Arc<dyn ExecutionPlan>,
    variation_storage: VariationLookupStorage,
    cache_schema: SchemaRef,
    cache_columns: Vec<String>,
    match_mode: KvMatchMode,
    exact_matcher: AlleleMatcher,
    schema: SchemaRef,
    vcf_has_chr: bool,
    vcf_zero_based: bool,
    cache_zero_based: bool,
    /// When true, probe multiple coordinate encodings (insertion-style,
    /// shifted deletions, tandem repeat window). When false, probe only
    /// the exact normalized interval.
    extended_probes: bool,
    /// Maximum allowed `failed` flag value from the cache.
    allowed_failed: i64,
    properties: Arc<PlanProperties>,
    /// Cache schema column positions for requested cache output columns.
    output_col_positions: Vec<usize>,
    /// Optional sink for co-located data collected during probe phase.
    colocated_sink: Option<ColocatedSink>,
    /// Optional partition-local sinks for co-located data collected during probe phase.
    colocated_partition_sinks: Option<Vec<ColocatedSink>>,
    target_partitions: usize,
    /// Per-contig Parquet variation lookup (shard + footer page directory), built
    /// once and shared across all per-partition streams via the Arc. Single-flight
    /// via OnceCell so simultaneous fan-out probes build it exactly once.
    #[cfg(feature = "parquet-cache")]
    parquet_lookup_cell: Arc<OnceCell<Arc<SinglePathParquetVariationLookup>>>,
}

#[derive(Debug, Clone)]
struct VariationLookupStorage {
    cache_root: PathBuf,
}

impl VariationLookupStorage {
    fn as_str(&self) -> &'static str {
        "parquet"
    }
}

fn build_lookup_output_schema(
    input_schema: SchemaRef,
    cache_schema: SchemaRef,
    cache_columns: &[String],
) -> (SchemaRef, Vec<usize>) {
    let mut output_col_positions = Vec::new();
    let mut fields: Vec<Arc<Field>> = input_schema.fields().iter().cloned().collect();
    for col_name in cache_columns {
        if let Ok(field) = cache_schema.field_with_name(col_name) {
            fields.push(Arc::new(Field::new(
                format!("cache_{}", field.name()),
                normalize_cache_output_type(field.data_type()),
                true,
            )));
            if let Ok(idx) = cache_schema.index_of(col_name) {
                output_col_positions.push(idx);
            }
        }
    }
    (
        Arc::new(datafusion::arrow::datatypes::Schema::new(fields)),
        output_col_positions,
    )
}

impl KvLookupExec {
    /// Construct a variation lookup exec backed by the Parquet cache. The read
    /// seam resolves through [`SinglePathParquetVariationLookup`].
    #[cfg(feature = "parquet-cache")]
    #[allow(clippy::too_many_arguments)]
    pub fn new_parquet(
        input: Arc<dyn ExecutionPlan>,
        cache_root: PathBuf,
        cache_schema: SchemaRef,
        cache_columns: Vec<String>,
        match_mode: KvMatchMode,
        exact_matcher: AlleleMatcher,
        vcf_has_chr: bool,
        vcf_zero_based: bool,
        cache_zero_based: bool,
        extended_probes: bool,
        allowed_failed: i64,
    ) -> Result<Self> {
        let input_schema = input.schema();
        let (schema, output_col_positions) =
            build_lookup_output_schema(input_schema, cache_schema.clone(), &cache_columns);

        let properties = Arc::new(PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            input.output_partitioning().clone(),
            EmissionType::Incremental,
            Boundedness::Bounded,
        ));

        Ok(Self {
            input,
            variation_storage: VariationLookupStorage { cache_root },
            cache_schema,
            cache_columns,
            match_mode,
            exact_matcher,
            schema,
            vcf_has_chr,
            vcf_zero_based,
            cache_zero_based,
            extended_probes,
            allowed_failed,
            properties,
            output_col_positions,
            colocated_sink: None,
            colocated_partition_sinks: None,
            target_partitions: 1,
            #[cfg(feature = "parquet-cache")]
            parquet_lookup_cell: Arc::new(OnceCell::new()),
        })
    }

    pub fn with_target_partitions(mut self, target_partitions: usize) -> Self {
        self.target_partitions = target_partitions.max(1);
        self
    }

    /// Set the co-located data sink for piggybacked collection during probe.
    pub fn with_colocated_sink(mut self, sink: ColocatedSink) -> Self {
        self.colocated_sink = Some(sink);
        self
    }

    /// Set one co-located data sink per input partition.
    pub fn with_colocated_partition_sinks(mut self, sinks: Vec<ColocatedSink>) -> Self {
        self.colocated_partition_sinks = Some(sinks);
        self
    }
}

impl Debug for KvLookupExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "KvLookupExec {{ mode: {:?}, cache_columns: {:?} }}",
            self.match_mode, self.cache_columns
        )
    }
}

impl DisplayAs for KvLookupExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "KvLookupExec: mode={:?}, columns={:?}",
            self.match_mode, self.cache_columns
        )
    }
}

impl ExecutionPlan for KvLookupExec {
    fn name(&self) -> &str {
        "KvLookupExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn properties(&self) -> &Arc<PlanProperties> {
        &self.properties
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.input]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        assert_eq!(children.len(), 1);
        let cache_root = self.variation_storage.cache_root.clone();
        let mut exec = KvLookupExec::new_parquet(
            children[0].clone(),
            cache_root,
            self.cache_schema.clone(),
            self.cache_columns.clone(),
            self.match_mode,
            self.exact_matcher,
            self.vcf_has_chr,
            self.vcf_zero_based,
            self.cache_zero_based,
            self.extended_probes,
            self.allowed_failed,
        )?;
        if let Some(sink) = &self.colocated_sink {
            exec = exec.with_colocated_sink(Arc::clone(sink));
        }
        if let Some(sinks) = &self.colocated_partition_sinks {
            exec = exec.with_colocated_partition_sinks(sinks.clone());
        }
        exec = exec.with_target_partitions(self.target_partitions);
        // Carry the shared per-contig variation lookup forward across re-planning
        // instead of letting new_parquet reset it.
        #[cfg(feature = "parquet-cache")]
        {
            exec.parquet_lookup_cell = Arc::clone(&self.parquet_lookup_cell);
        }
        Ok(Arc::new(exec))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let input_stream = self.input.execute(partition, context)?;
        let colocated_sink = self
            .colocated_partition_sinks
            .as_ref()
            .and_then(|sinks| sinks.get(partition).cloned())
            .or_else(|| self.colocated_sink.clone());

        #[cfg_attr(not(feature = "parquet-cache"), allow(unused_mut))]
        let mut stream = KvLookupStream::new(
            input_stream,
            self.variation_storage.clone(),
            self.cache_schema.clone(),
            self.schema.clone(),
            self.cache_columns.clone(),
            self.match_mode,
            self.exact_matcher,
            self.vcf_has_chr,
            self.vcf_zero_based,
            self.cache_zero_based,
            self.extended_probes,
            self.allowed_failed,
            self.output_col_positions.clone(),
            colocated_sink,
            self.target_partitions,
        );

        // Share the exec's single per-contig variation lookup with this stream;
        // new() seeded a placeholder empty cell.
        #[cfg(feature = "parquet-cache")]
        {
            stream.parquet_lookup_cell = Arc::clone(&self.parquet_lookup_cell);
        }

        Ok(Box::pin(stream))
    }
}

/// Streaming implementation that processes VCF batches and probes the KV store.
///
/// When a colocated sink is present, batches are buffered during the probe
/// phase and only emitted after the input stream is exhausted. This ensures
/// the colocated sink is fully populated before downstream consumers build
/// the colocated map.
struct KvLookupStream {
    input: SendableRecordBatchStream,
    variation_storage: VariationLookupStorage,
    cache_schema: SchemaRef,
    schema: SchemaRef,
    cache_columns: Vec<String>,
    _match_mode: KvMatchMode,
    exact_matcher: AlleleMatcher,
    vcf_has_chr: bool,
    vcf_zero_based: bool,
    cache_zero_based: bool,
    extended_probes: bool,
    allowed_failed: i64,
    output_col_positions: Vec<usize>,
    colocated_sink: Option<ColocatedSink>,
    target_partitions: usize,
    /// Shared (Arc-cloned from the exec) per-contig Parquet variation lookup,
    /// built once.
    #[cfg(feature = "parquet-cache")]
    parquet_lookup_cell: Arc<OnceCell<Arc<SinglePathParquetVariationLookup>>>,
    warm_cold_backend: WarmColdVariationBackend,
    warm_cold_index_mode: WarmColdVariationIndexMode,
    profile_enabled: bool,
    profile_detailed: bool,
    profile_emitted: bool,
    profile: LookupProfile,
    /// Input slices waiting to be processed after a large upstream batch was
    /// split to bound Parquet `take_rows` work.
    pending_input_slices: VecDeque<RecordBatch>,
    /// Buffered matched batches (used when colocated sink is present).
    /// Batches are collected during probe and emitted after input exhaustion.
    matched_batches: VecDeque<RecordBatch>,
    /// True once the input stream is exhausted and we're emitting buffered batches.
    input_exhausted: bool,
}

struct WarmColocIndices {
    variation_name: usize,
    end_col: Option<usize>,
    failed: Option<usize>,
    somatic: Option<usize>,
    pheno: Option<usize>,
    clin_sig: Option<usize>,
    clin_sig_allele: Option<usize>,
    pubmed: Option<usize>,
    af_indices: Vec<Option<usize>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LookupDecision {
    EmitWarmExact,
    SkipFjall,
    UseFjall,
}

#[derive(Debug, Clone)]
struct PendingColdProbe {
    chrom: String,
    probe_start: i64,
    position_key: u64,
    vcf_ref: String,
    vcf_alt: String,
    vcf_iv_start: i64,
    vcf_iv_end: i64,
    vcf_row: u32,
}

/// Outcome of probing a single cold-tier position for a VCF variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ColdProbeResult {
    Match,
    PositionCoveredNoExact,
    NotCovered,
}

/// Build an Arrow `ArrayBuilder` for the given output column type.
fn make_builder(
    dt: &DataType,
    capacity: usize,
) -> Result<Box<dyn datafusion::arrow::array::ArrayBuilder>> {
    use datafusion::arrow::array::{
        BooleanBuilder, Float32Builder, Float64Builder, Int8Builder, Int16Builder, Int32Builder,
        Int64Builder, StringBuilder, UInt8Builder, UInt16Builder, UInt32Builder, UInt64Builder,
    };
    match dt {
        DataType::Int32 => Ok(Box::new(Int32Builder::with_capacity(capacity))),
        DataType::Int64 => Ok(Box::new(Int64Builder::with_capacity(capacity))),
        DataType::Float32 => Ok(Box::new(Float32Builder::with_capacity(capacity))),
        DataType::Float64 => Ok(Box::new(Float64Builder::with_capacity(capacity))),
        DataType::UInt32 => Ok(Box::new(UInt32Builder::with_capacity(capacity))),
        DataType::UInt64 => Ok(Box::new(UInt64Builder::with_capacity(capacity))),
        DataType::Int8 => Ok(Box::new(Int8Builder::with_capacity(capacity))),
        DataType::Int16 => Ok(Box::new(Int16Builder::with_capacity(capacity))),
        DataType::UInt8 => Ok(Box::new(UInt8Builder::with_capacity(capacity))),
        DataType::UInt16 => Ok(Box::new(UInt16Builder::with_capacity(capacity))),
        DataType::Boolean => Ok(Box::new(BooleanBuilder::with_capacity(capacity))),
        // For output, we always use Utf8 (normalized from Utf8View/LargeUtf8).
        DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => Ok(Box::new(
            StringBuilder::with_capacity(capacity, capacity * 16),
        )),
        other => Err(DataFusionError::Execution(format!(
            "make_builder: unsupported type {other}"
        ))),
    }
}

fn push_unique_column(columns: &mut Vec<String>, name: &str) {
    if !columns.iter().any(|existing| existing == name) {
        columns.push(name.to_string());
    }
}

/// Columns to project from the Parquet variation dataset for a cold-tier probe.
fn cold_parquet_projection_columns(
    cache_columns: &[String],
    include_colocated: bool,
) -> Vec<String> {
    let mut columns = Vec::with_capacity(cache_columns.len() + 16);
    for name in ["position_key", "allele_string", "end", "failed"] {
        push_unique_column(&mut columns, name);
    }
    for name in cache_columns {
        push_unique_column(&mut columns, name);
    }
    if include_colocated {
        for name in [
            "variation_name",
            "end",
            "failed",
            "somatic",
            "phenotype_or_disease",
            "clin_sig",
            "clin_sig_allele",
            "pubmed",
        ] {
            push_unique_column(&mut columns, name);
        }
        for name in crate::colocated::AF_COL_NAMES {
            push_unique_column(&mut columns, name);
        }
    }
    columns
}

#[derive(Default)]
struct ColdChunkProbeMetrics {
    append_elapsed: Duration,
    colocated_prepare_elapsed: Duration,
    colocated_match_elapsed: Duration,
    primary_allele_rows: u64,
    exact_match_calls: u64,
    colocated_allele_rows: u64,
    colocated_entries: u64,
    cold_rows_scanned: u64,
    emitted: bool,
}

/// Variation cold-tier backend. The single-path cache is Parquet-only; the enum
/// is retained so the profiling/diagnostic plumbing keeps its `backend` label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WarmColdVariationBackend {
    Parquet,
}

impl WarmColdVariationBackend {
    fn as_str(self) -> &'static str {
        match self {
            Self::Parquet => "parquet",
        }
    }

    #[cfg(feature = "parquet-cache")]
    fn is_parquet(self) -> bool {
        matches!(self, Self::Parquet)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WarmColdVariationIndexMode {
    PositionThenVariantBloom,
}

impl WarmColdVariationIndexMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::PositionThenVariantBloom => "posidx_bloom",
        }
    }
}

fn format_variation_lookup_profile_line(
    storage: &str,
    warm_source: &str,
    warm_cold_backend: WarmColdVariationBackend,
    warm_cold_index_mode: WarmColdVariationIndexMode,
    cache_root: &str,
    warm_dir: Option<&str>,
    cold_dir: Option<&str>,
    position_index_dir: Option<&str>,
    variant_bloom_index_dir: Option<&str>,
) -> String {
    format!(
        "[vep-kv-profile-detail] variation_lookup storage={} warm_source={} warm_cold_backend={} index_mode={} cache_root={} warm_dir={} cold_dir={} position_index_dir={} variant_bloom_index_dir={}",
        storage,
        warm_source,
        warm_cold_backend.as_str(),
        warm_cold_index_mode.as_str(),
        cache_root,
        warm_dir.unwrap_or("-"),
        cold_dir.unwrap_or("-"),
        position_index_dir.unwrap_or("-"),
        variant_bloom_index_dir.unwrap_or("-"),
    )
}

fn warm_source_label(backend: WarmColdVariationBackend) -> &'static str {
    match backend {
        #[cfg(feature = "parquet-cache")]
        WarmColdVariationBackend::Parquet => "parquet",
    }
}

#[cfg(feature = "parquet-cache")]
fn parquet_variation_path_for_chrom(
    cache: &crate::parquet_cache::detect::PartitionedParquetCache,
    chrom: &str,
) -> Option<PathBuf> {
    cache.variation_path(chrom).or_else(|| {
        chrom
            .strip_prefix("chr")
            .and_then(|bare| cache.variation_path(bare))
            .or_else(|| {
                if chrom.starts_with("chr") {
                    None
                } else {
                    cache.variation_path(&format!("chr{chrom}"))
                }
            })
    })
}

#[derive(Default)]
struct LookupProfile {
    batches: u64,
    input_rows: u64,
    output_rows: u64,
    probes: u64,
    point_gets: u64,
    range_prefetch_batches: u64,
    range_prefetch_entries: u64,
    range_prefetch_bytes: u64,
    range_prefetch_skipped: u64,
    range_prefetch_skip_multi_chrom: u64,
    range_prefetch_skip_span: u64,
    range_prefetch_skip_entries: u64,
    range_prefetch_skip_bytes: u64,
    range_prefetch: Duration,
    extract_cols: Duration,
    match_loop: Duration,
    vcf_take: Duration,
    cache_build: Duration,
    probe_build: Duration,
    point_get_raw: Duration,
    prefetch_map_lookup: Duration,
    decompress: Duration,
    reader_init: Duration,
    primary_match: Duration,
    cache_column_append: Duration,
    colocated_prepare: Duration,
    colocated_match: Duration,
    colocated_flush: Duration,
    null_append: Duration,
    warm_probe: Duration,
    warm_chunk_load: Duration,
    position_index_load: Duration,
    variant_bloom_load: Duration,
    raw_get_hits: u64,
    raw_get_misses: u64,
    prefetch_hits: u64,
    prefetch_misses: u64,
    decode_calls: u64,
    compressed_bytes: u64,
    decompressed_bytes: u64,
    primary_allele_rows: u64,
    primary_failed_skips: u64,
    primary_interval_skips: u64,
    exact_match_calls: u64,
    primary_matches: u64,
    colocated_allele_rows: u64,
    colocated_entries: u64,
    null_rows: u64,
    warm_probes: u64,
    warm_matches: u64,
    warm_definitive_misses: u64,
    warm_not_covered: u64,
    warm_boundary_fallbacks: u64,
    warm_candidate_rows: u64,
    warm_rows_scanned: u64,
    warm_chunks_loaded: u64,
    warm_chunk_rows: u64,
    position_index_checks: u64,
    position_index_negative_skips: u64,
    position_index_positive_gets: u64,
    position_index_loaded: u64,
    position_index_persisted_loads: u64,
    position_index_parquet_fallback_loads: u64,
    position_index_rows: u64,
    position_index_bytes: u64,
    variant_bloom_checks: u64,
    variant_bloom_negative_skips: u64,
    variant_bloom_positive_gets: u64,
    variant_bloom_loaded: u64,
    variant_bloom_entries: u64,
    variant_bloom_bits: u64,
    variant_bloom_hashes: u64,
    variant_bloom_bytes: u64,
    cold_parquet_probes: u64,
    cold_parquet_matches: u64,
    cold_parquet_position_misses: u64,
    cold_parquet_not_covered: u64,
    cold_parquet_rows_scanned: u64,
    cold_parquet_row_groups_total: u64,
    cold_parquet_row_groups_unique_touched: u64,
    cold_parquet_row_group_metadata_probes: u64,
    cold_parquet_row_group_current_hits: u64,
    cold_parquet_row_group_previous_hits: u64,
    cold_parquet_row_group_advanced_hits: u64,
    cold_parquet_row_group_binary_search_hits: u64,
    cold_parquet_row_group_metadata_misses: u64,
    cold_parquet_row_group_skipped_ahead: u64,
    cold_parquet_row_groups_loaded: u64,
    cold_parquet_row_group_load_batches: u64,
    cold_parquet_row_group_cache_hits: u64,
    cold_parquet_row_group_cache_misses: u64,
    cold_parquet_rows_loaded: u64,
    cold_parquet_position_page_index_loaded: bool,
    cold_parquet_position_column_index_loaded: bool,
    cold_parquet_position_pages_total: u64,
    cold_parquet_position_bloom_filter_row_groups: u64,
    cold_parquet_position_bloom_filter_bytes: u64,
    cold_parquet_page_index_probes: u64,
    cold_parquet_page_index_available_probes: u64,
    cold_parquet_page_index_unavailable_probes: u64,
    cold_parquet_page_index_pages_in_probed_row_groups: u64,
    cold_parquet_page_index_candidate_pages: u64,
    cold_parquet_page_index_candidate_rows: u64,
    cold_parquet_page_index_candidate_misses: u64,
    cold_parquet_page_index_unique_candidate_pages: u64,
    cold_parquet_page_index_unique_candidate_rows: u64,
    /// Cold-tier read time. Emitted as `cold_tier_load` (and the read-stats line
    /// as `cold_tier`) — the cold tier is backend-agnostic and is Parquet in the
    /// single-path cache; the `cold_parquet*` field names are retained from the
    /// former parquet cold tier and only matter for the parquet-format diagnostics
    /// (row_groups/pages) below, which stay zero under the Parquet backend.
    cold_parquet_load: Duration,
}

impl LookupProfile {
    fn total_known(&self) -> Duration {
        self.extract_cols + self.range_prefetch + self.match_loop + self.vcf_take + self.cache_build
    }

    fn pct(stage: Duration, total: Duration) -> f64 {
        if total.is_zero() {
            0.0
        } else {
            stage.as_secs_f64() * 100.0 / total.as_secs_f64()
        }
    }

    fn detail_known(&self) -> Duration {
        self.probe_build
            + self.point_get_raw
            + self.prefetch_map_lookup
            + self.decompress
            + self.reader_init
            + self.primary_match
            + self.cache_column_append
            + self.colocated_prepare
            + self.colocated_match
            + self.colocated_flush
            + self.null_append
            + self.warm_probe
            + self.warm_chunk_load
            + self.position_index_load
            + self.variant_bloom_load
            + self.cold_parquet_load
    }

    fn detail_lines(&self) -> Vec<String> {
        let detail_total = self.detail_known();
        vec![
            format!(
                "[vep-kv-profile-detail] stages total_s={:.3} probe_build={:.3}s ({:.1}%) point_get_raw={:.3}s ({:.1}%) prefetch_map_lookup={:.3}s ({:.1}%) decompress={:.3}s ({:.1}%) reader_init={:.3}s ({:.1}%) primary_match={:.3}s ({:.1}%) cache_column_append={:.3}s ({:.1}%) colocated_prepare={:.3}s ({:.1}%) colocated_match={:.3}s ({:.1}%) colocated_flush={:.3}s ({:.1}%) null_append={:.3}s ({:.1}%) warm_probe={:.3}s ({:.1}%) warm_chunk_load={:.3}s ({:.1}%) position_index_load={:.3}s ({:.1}%) variant_bloom_load={:.3}s ({:.1}%) cold_tier_load={:.3}s ({:.1}%)",
                detail_total.as_secs_f64(),
                self.probe_build.as_secs_f64(),
                Self::pct(self.probe_build, detail_total),
                self.point_get_raw.as_secs_f64(),
                Self::pct(self.point_get_raw, detail_total),
                self.prefetch_map_lookup.as_secs_f64(),
                Self::pct(self.prefetch_map_lookup, detail_total),
                self.decompress.as_secs_f64(),
                Self::pct(self.decompress, detail_total),
                self.reader_init.as_secs_f64(),
                Self::pct(self.reader_init, detail_total),
                self.primary_match.as_secs_f64(),
                Self::pct(self.primary_match, detail_total),
                self.cache_column_append.as_secs_f64(),
                Self::pct(self.cache_column_append, detail_total),
                self.colocated_prepare.as_secs_f64(),
                Self::pct(self.colocated_prepare, detail_total),
                self.colocated_match.as_secs_f64(),
                Self::pct(self.colocated_match, detail_total),
                self.colocated_flush.as_secs_f64(),
                Self::pct(self.colocated_flush, detail_total),
                self.null_append.as_secs_f64(),
                Self::pct(self.null_append, detail_total),
                self.warm_probe.as_secs_f64(),
                Self::pct(self.warm_probe, detail_total),
                self.warm_chunk_load.as_secs_f64(),
                Self::pct(self.warm_chunk_load, detail_total),
                self.position_index_load.as_secs_f64(),
                Self::pct(self.position_index_load, detail_total),
                self.variant_bloom_load.as_secs_f64(),
                Self::pct(self.variant_bloom_load, detail_total),
                self.cold_parquet_load.as_secs_f64(),
                Self::pct(self.cold_parquet_load, detail_total),
            ),
            format!(
                "[vep-kv-profile-detail] io raw_get_hits={} raw_get_misses={} prefetch_hits={} prefetch_misses={} decode_calls={} compressed_bytes={} decompressed_bytes={}",
                self.raw_get_hits,
                self.raw_get_misses,
                self.prefetch_hits,
                self.prefetch_misses,
                self.decode_calls,
                self.compressed_bytes,
                self.decompressed_bytes,
            ),
            format!(
                "[vep-kv-profile-detail] match primary_allele_rows={} primary_failed_skips={} primary_interval_skips={} exact_match_calls={} primary_matches={} colocated_allele_rows={} colocated_entries={} null_rows={}",
                self.primary_allele_rows,
                self.primary_failed_skips,
                self.primary_interval_skips,
                self.exact_match_calls,
                self.primary_matches,
                self.colocated_allele_rows,
                self.colocated_entries,
                self.null_rows,
            ),
            format!(
                "[vep-kv-profile-detail] warm probes={} matches={} definitive_misses={} not_covered={} boundary_fallbacks={} candidate_rows={} rows_scanned={} chunks_loaded={} chunk_rows={}",
                self.warm_probes,
                self.warm_matches,
                self.warm_definitive_misses,
                self.warm_not_covered,
                self.warm_boundary_fallbacks,
                self.warm_candidate_rows,
                self.warm_rows_scanned,
                self.warm_chunks_loaded,
                self.warm_chunk_rows,
            ),
            format!(
                "[vep-kv-profile-detail] position_index checks={} negative_skips={} positive_gets={} loaded={} persisted_loads={} parquet_fallback_loads={} rows={} bytes={}",
                self.position_index_checks,
                self.position_index_negative_skips,
                self.position_index_positive_gets,
                self.position_index_loaded,
                self.position_index_persisted_loads,
                self.position_index_parquet_fallback_loads,
                self.position_index_rows,
                self.position_index_bytes,
            ),
            format!(
                "[vep-kv-profile-detail] variant_bloom checks={} negative_skips={} positive_gets={} loaded={} entries={} bits={} hashes={} bytes={}",
                self.variant_bloom_checks,
                self.variant_bloom_negative_skips,
                self.variant_bloom_positive_gets,
                self.variant_bloom_loaded,
                self.variant_bloom_entries,
                self.variant_bloom_bits,
                self.variant_bloom_hashes,
                self.variant_bloom_bytes,
            ),
            format!(
                "[vep-kv-profile-detail] cold_tier probes={} matches={} position_misses={} not_covered={} rows_scanned={}",
                self.cold_parquet_probes,
                self.cold_parquet_matches,
                self.cold_parquet_position_misses,
                self.cold_parquet_not_covered,
                self.cold_parquet_rows_scanned,
            ),
            format!(
                "[vep-kv-profile-detail] cold_parquet_row_groups total={} unique_touched={} untouched={} metadata_probes={} current_hits={} previous_hits={} advanced_hits={} binary_search_hits={} metadata_misses={} skipped_ahead={} loaded={} load_batches={} cache_hits={} cache_misses={} rows_loaded={}",
                self.cold_parquet_row_groups_total,
                self.cold_parquet_row_groups_unique_touched,
                self.cold_parquet_row_groups_total
                    .saturating_sub(self.cold_parquet_row_groups_unique_touched),
                self.cold_parquet_row_group_metadata_probes,
                self.cold_parquet_row_group_current_hits,
                self.cold_parquet_row_group_previous_hits,
                self.cold_parquet_row_group_advanced_hits,
                self.cold_parquet_row_group_binary_search_hits,
                self.cold_parquet_row_group_metadata_misses,
                self.cold_parquet_row_group_skipped_ahead,
                self.cold_parquet_row_groups_loaded,
                self.cold_parquet_row_group_load_batches,
                self.cold_parquet_row_group_cache_hits,
                self.cold_parquet_row_group_cache_misses,
                self.cold_parquet_rows_loaded,
            ),
            format!(
                "[vep-kv-profile-detail] cold_parquet_pages position_page_index_loaded={} position_column_index_loaded={} position_pages_total={} position_bloom_filter_row_groups={} position_bloom_filter_bytes={} page_index_probes={} available_probes={} unavailable_probes={} pages_in_probed_row_groups={} candidate_pages={} candidate_rows={} unique_candidate_pages={} unique_candidate_rows={} candidate_misses={}",
                self.cold_parquet_position_page_index_loaded,
                self.cold_parquet_position_column_index_loaded,
                self.cold_parquet_position_pages_total,
                self.cold_parquet_position_bloom_filter_row_groups,
                self.cold_parquet_position_bloom_filter_bytes,
                self.cold_parquet_page_index_probes,
                self.cold_parquet_page_index_available_probes,
                self.cold_parquet_page_index_unavailable_probes,
                self.cold_parquet_page_index_pages_in_probed_row_groups,
                self.cold_parquet_page_index_candidate_pages,
                self.cold_parquet_page_index_candidate_rows,
                self.cold_parquet_page_index_unique_candidate_pages,
                self.cold_parquet_page_index_unique_candidate_rows,
                self.cold_parquet_page_index_candidate_misses,
            ),
        ]
    }

    fn emit(&self, detailed: bool) {
        let total = self.total_known();
        let input_rate = if total.is_zero() {
            0.0
        } else {
            self.input_rows as f64 / total.as_secs_f64()
        };
        let output_rate = if total.is_zero() {
            0.0
        } else {
            self.output_rows as f64 / total.as_secs_f64()
        };
        eprintln!(
            "[vep-kv-profile] batches={} input_rows={} output_rows={} probes={} point_gets={} range_prefetch_batches={} range_prefetch_entries={} range_prefetch_bytes={} range_prefetch_skipped={} total_s={:.3} input_rows_per_s={:.1} output_rows_per_s={:.1}",
            self.batches,
            self.input_rows,
            self.output_rows,
            self.probes,
            self.point_gets,
            self.range_prefetch_batches,
            self.range_prefetch_entries,
            self.range_prefetch_bytes,
            self.range_prefetch_skipped,
            total.as_secs_f64(),
            input_rate,
            output_rate
        );
        eprintln!(
            "[vep-kv-profile] range_prefetch_skip_reasons multi_chrom={} span={} entries={} bytes={}",
            self.range_prefetch_skip_multi_chrom,
            self.range_prefetch_skip_span,
            self.range_prefetch_skip_entries,
            self.range_prefetch_skip_bytes,
        );
        eprintln!(
            "[vep-kv-profile] extract_cols={:.3}s ({:.1}%) range_prefetch={:.3}s ({:.1}%) match_loop={:.3}s ({:.1}%) vcf_take={:.3}s ({:.1}%) cache_build={:.3}s ({:.1}%)",
            self.extract_cols.as_secs_f64(),
            Self::pct(self.extract_cols, total),
            self.range_prefetch.as_secs_f64(),
            Self::pct(self.range_prefetch, total),
            self.match_loop.as_secs_f64(),
            Self::pct(self.match_loop, total),
            self.vcf_take.as_secs_f64(),
            Self::pct(self.vcf_take, total),
            self.cache_build.as_secs_f64(),
            Self::pct(self.cache_build, total),
        );
        if detailed {
            for line in self.detail_lines() {
                eprintln!("{line}");
            }
        }
    }
}

fn kv_profile_enabled() -> bool {
    std::env::var_os("VEP_KV_PROFILE").is_some()
        || std::env::var_os("VEP_KV_PROFILE_DETAILED").is_some()
}

fn kv_profile_detailed_enabled() -> bool {
    std::env::var_os("VEP_KV_PROFILE_DETAILED").is_some()
}

fn lookup_process_batch_rows() -> usize {
    std::env::var("VEP_LANCE_LOOKUP_PROCESS_BATCH_ROWS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|value| *value > 0)
        .unwrap_or(DEFAULT_LANCE_LOOKUP_PROCESS_BATCH_ROWS)
}

fn enqueue_record_batch_slices(
    queue: &mut VecDeque<RecordBatch>,
    batch: RecordBatch,
    max_rows: usize,
) {
    let max_rows = max_rows.max(1);
    if batch.num_rows() <= max_rows {
        queue.push_back(batch);
        return;
    }

    let mut offset = 0usize;
    while offset < batch.num_rows() {
        let len = max_rows.min(batch.num_rows() - offset);
        queue.push_back(batch.slice(offset, len));
        offset += len;
    }
}

enum StringColumnView<'a> {
    Utf8(&'a StringArray),
    Utf8View(&'a StringViewArray),
    LargeUtf8(&'a LargeStringArray),
}

impl<'a> StringColumnView<'a> {
    fn value_or_empty(&self, row: usize) -> &'a str {
        match self {
            Self::Utf8(arr) => {
                if arr.is_null(row) {
                    ""
                } else {
                    arr.value(row)
                }
            }
            Self::Utf8View(arr) => {
                if arr.is_null(row) {
                    ""
                } else {
                    arr.value(row)
                }
            }
            Self::LargeUtf8(arr) => {
                if arr.is_null(row) {
                    ""
                } else {
                    arr.value(row)
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn probe_taken_batch_position(
    batch: &RecordBatch,
    rows: &[u32],
    exact_matcher: AlleleMatcher,
    allowed_failed: i64,
    profile_detailed: bool,
    chrom: &str,
    probe_start: i64,
    vcf_ref: &str,
    vcf_alt: &str,
    vcf_iv_start: i64,
    vcf_iv_end: i64,
    emit_output: bool,
    vcf_row: u32,
    cache_columns: &[String],
    col_map: &[usize],
    builders: &mut [Box<dyn datafusion::arrow::array::ArrayBuilder>],
    vcf_indices: &mut Vec<u32>,
    coloc_buf: Option<&mut HashMap<ColocatedKey, ColocatedSinkValue>>,
) -> Result<(ColdProbeResult, ColdChunkProbeMetrics)> {
    struct PreparedColoc {
        chrom_norm: String,
        input_start: i64,
        input_allele_string: String,
        compare_allele_string: String,
        vep_start: i64,
        vep_end: i64,
        compare_output_allele: Option<String>,
        unshifted_output_allele: Option<String>,
    }

    let mut metrics = ColdChunkProbeMetrics::default();
    if rows.is_empty() {
        return Ok((ColdProbeResult::NotCovered, metrics));
    }
    let indices = BatchProbeIndices::new(batch)?;
    let mut coloc_buf = coloc_buf;
    let prepared_coloc = if coloc_buf.is_some() {
        let prepare_started = Instant::now();
        let (input_ref, input_alt, input_start) =
            vcf_to_vep_input_allele(vcf_iv_start, vcf_ref, vcf_alt);
        let input_allele_string = format!("{input_ref}/{input_alt}");
        let (compare_ref, compare_alt) = vcf_to_vep_allele(vcf_ref, vcf_alt);
        let compare_allele_string = format!("{compare_ref}/{compare_alt}");
        let prepared = PreparedColoc {
            chrom_norm: chrom.to_string(),
            input_start,
            input_allele_string,
            compare_output_allele: output_allele_from_allele_string(&compare_allele_string)
                .map(str::to_string),
            unshifted_output_allele: None,
            compare_allele_string,
            vep_start: vep_norm_start(vcf_iv_start, vcf_ref, vcf_alt),
            vep_end: vep_norm_end(vcf_iv_start, vcf_ref, vcf_alt),
        };
        metrics.colocated_prepare_elapsed += prepare_started.elapsed();
        Some(prepared)
    } else {
        None
    };
    let coloc_visible = prepared_coloc.as_ref().is_some_and(|prepared| {
        probe_start_visible_to_window(probe_start, prepared.vep_start, prepared.vep_end)
    });
    let coloc_indices = if coloc_visible {
        resolve_batch_coloc_indices(batch)
    } else {
        None
    };
    // Build zero-copy AF column refs ONCE per batch (Arc::clone of the 27 AF
    // columns, no string copy); every matched entry shares this by Arc + row index.
    let af_columns = coloc_indices
        .as_ref()
        .map(|ci| AfColumns::from_batch(batch, &ci.af_indices));

    let mut matched = false;
    for &row in rows {
        let row_usize = row as usize;
        let Some(allele_string) =
            batch_string_value(batch, Some(indices.allele_string), row_usize)?
        else {
            continue;
        };

        metrics.cold_rows_scanned += 1;
        if profile_detailed {
            metrics.primary_allele_rows += 1;
            metrics.exact_match_calls += 1;
        }

        if let (Some(buf), Some(prepared), Some(ci), Some(af_cols)) = (
            coloc_buf.as_deref_mut(),
            prepared_coloc.as_ref(),
            coloc_indices.as_ref(),
            af_columns.as_ref(),
        ) {
            let colocated_started = Instant::now();
            if profile_detailed {
                metrics.colocated_allele_rows += 1;
            }

            let failed = batch_i64_value(batch, ci.failed, row_usize).unwrap_or(0);
            if failed <= allowed_failed
                && let Some(var_name) =
                    batch_string_value(batch, Some(ci.variation_name), row_usize)?
                && !var_name.is_empty()
            {
                let existing_end =
                    batch_i64_value(batch, ci.end_col, row_usize).unwrap_or(probe_start);
                if let Some(matched_alleles) = compare_existing_variant_alleles(
                    &prepared.compare_allele_string,
                    prepared.vep_start,
                    prepared.vep_end,
                    None,
                    None,
                    None,
                    &allele_string,
                    probe_start,
                    existing_end,
                ) {
                    let key: ColocatedKey = (
                        prepared.chrom_norm.clone(),
                        prepared.input_start,
                        vcf_iv_end,
                        prepared.input_allele_string.clone(),
                    );
                    let sink_value = buf.entry(key).or_insert_with(|| ColocatedSinkValue {
                        entries: Vec::new(),
                        compare_output_allele: prepared.compare_output_allele.clone(),
                        unshifted_output_allele: prepared.unshifted_output_allele.clone(),
                    });
                    sink_value.entries.push(ColocatedCacheEntry {
                        variation_name: var_name,
                        allele_string: allele_string.clone(),
                        matched_alleles,
                        somatic: batch_i64_value(batch, ci.somatic, row_usize).unwrap_or(0),
                        pheno: batch_i64_value(batch, ci.pheno, row_usize).unwrap_or(0),
                        clin_sig: batch_string_value(batch, ci.clin_sig, row_usize)?,
                        clin_sig_allele: batch_string_value(batch, ci.clin_sig_allele, row_usize)?,
                        pubmed: batch_string_value(batch, ci.pubmed, row_usize)?,
                        af: af_cols.clone(),
                        af_row: row,
                    });
                    if profile_detailed {
                        metrics.colocated_entries += 1;
                    }
                }
            }
            metrics.colocated_match_elapsed += colocated_started.elapsed();
        }

        if matched {
            continue;
        }

        let failed = batch_i64_value(batch, indices.failed, row_usize).unwrap_or(0);
        if failed > allowed_failed {
            continue;
        }

        let existing_end = batch_i64_value(batch, indices.end, row_usize).unwrap_or(probe_start);
        let cache_iv_start = probe_start.min(existing_end);
        let cache_iv_end = probe_start.max(existing_end);
        if cache_iv_start > vcf_iv_end || cache_iv_end < vcf_iv_start {
            continue;
        }

        if exact_matcher(vcf_ref, vcf_alt, &allele_string) {
            if emit_output {
                let append_started = Instant::now();
                let Some(output_indices) = batch_output_indices(batch, cache_columns, col_map)
                else {
                    return Err(DataFusionError::Execution(
                        "Lance batch missing one or more requested cache output columns".into(),
                    ));
                };
                vcf_indices.push(vcf_row);
                append_batch_row_values(batch, row_usize, &output_indices, builders)?;
                metrics.append_elapsed += append_started.elapsed();
                metrics.emitted = true;
            }
            matched = true;
        }
    }

    if matched {
        Ok((ColdProbeResult::Match, metrics))
    } else {
        Ok((ColdProbeResult::PositionCoveredNoExact, metrics))
    }
}

#[cfg(feature = "parquet-cache")]
struct BatchProbeIndices {
    allele_string: usize,
    end: Option<usize>,
    failed: Option<usize>,
}

#[cfg(feature = "parquet-cache")]
impl BatchProbeIndices {
    fn new(batch: &RecordBatch) -> Result<Self> {
        let schema = batch.schema();
        Ok(Self {
            allele_string: schema.index_of("allele_string").map_err(|_| {
                DataFusionError::Execution("Lance batch missing allele_string".into())
            })?,
            end: schema.index_of("end").ok(),
            failed: schema.index_of("failed").ok(),
        })
    }
}

#[cfg(feature = "parquet-cache")]
fn start_row_map(batch: &RecordBatch) -> Result<HashMap<u32, Vec<u32>>> {
    let schema = batch.schema();
    let start_idx = schema
        .index_of("start")
        .map_err(|_| DataFusionError::Execution("Lance batch missing start".into()))?;
    let starts = batch
        .column(start_idx)
        .as_any()
        .downcast_ref::<UInt32Array>()
        .ok_or_else(|| DataFusionError::Execution("Lance batch start must be UInt32".into()))?;
    let mut rows = HashMap::<u32, Vec<u32>>::new();
    for row in 0..batch.num_rows() {
        if !starts.is_null(row) {
            rows.entry(starts.value(row)).or_default().push(row as u32);
        }
    }
    Ok(rows)
}

impl KvLookupStream {
    #[allow(clippy::too_many_arguments)]
    fn new(
        input: SendableRecordBatchStream,
        variation_storage: VariationLookupStorage,
        cache_schema: SchemaRef,
        schema: SchemaRef,
        cache_columns: Vec<String>,
        match_mode: KvMatchMode,
        exact_matcher: AlleleMatcher,
        vcf_has_chr: bool,
        vcf_zero_based: bool,
        cache_zero_based: bool,
        extended_probes: bool,
        allowed_failed: i64,
        output_col_positions: Vec<usize>,
        colocated_sink: Option<ColocatedSink>,
        target_partitions: usize,
    ) -> Self {
        let cache_root = &variation_storage.cache_root;
        let warm_cold_backend = WarmColdVariationBackend::Parquet;
        let warm_cold_index_mode = WarmColdVariationIndexMode::PositionThenVariantBloom;
        let profile_enabled = kv_profile_enabled();
        let profile_detailed = kv_profile_detailed_enabled();
        if profile_detailed || std::env::var_os("VEP_LANCE_PROFILE").is_some() {
            let cache_root_label = cache_root.display().to_string();
            eprintln!(
                "{}",
                format_variation_lookup_profile_line(
                    variation_storage.as_str(),
                    warm_source_label(warm_cold_backend),
                    warm_cold_backend,
                    warm_cold_index_mode,
                    &cache_root_label,
                    None,
                    None,
                    None,
                    None,
                )
            );
        }

        Self {
            input,
            variation_storage: variation_storage.clone(),
            cache_schema,
            schema,
            cache_columns,
            _match_mode: match_mode,
            exact_matcher,
            vcf_has_chr,
            vcf_zero_based,
            cache_zero_based,
            extended_probes,
            allowed_failed,
            output_col_positions,
            colocated_sink,
            target_partitions: target_partitions.max(1),
            // Placeholder; execute() overwrites with the exec's shared cell.
            #[cfg(feature = "parquet-cache")]
            parquet_lookup_cell: Arc::new(OnceCell::new()),
            warm_cold_backend,
            warm_cold_index_mode,
            profile_enabled,
            profile_detailed,
            profile_emitted: false,
            profile: LookupProfile::default(),
            pending_input_slices: VecDeque::new(),
            matched_batches: VecDeque::new(),
            input_exhausted: false,
        }
    }

    /// Single-pass position-keyed lookup.
    ///
    /// For each VCF row, fetch the per-position entry from Parquet, match alleles,
    /// and append matched column values directly into ArrayBuilders.
    #[allow(clippy::too_many_arguments)]
    fn probe_warm_position(
        &mut self,
        chrom: &str,
        _chrom_code: u16,
        _probe_start: i64,
        _vcf_ref: &str,
        _vcf_alt: &str,
        _vcf_iv_start: i64,
        _vcf_iv_end: i64,
        emit_output: bool,
        _vcf_row: u32,
        cache_columns: &[String],
        _col_map: &[usize],
        _builders: &mut [Box<dyn datafusion::arrow::array::ArrayBuilder>],
        _vcf_indices: &mut Vec<u32>,
        coloc_buf: Option<&mut HashMap<ColocatedKey, ColocatedSinkValue>>,
    ) -> Result<Option<LookupDecision>> {
        let collect_colocated = coloc_buf.is_some();
        self.ensure_parquet_lookup(chrom, cache_columns, collect_colocated)?;
        if self.profile_enabled {
            self.profile.warm_probes += 1;
            self.profile.warm_not_covered += 1;
        }

        let _ = emit_output;
        Ok(Some(LookupDecision::UseFjall))
    }

    /// Single-flight build of the shared per-contig Parquet variation lookup.
    #[cfg(feature = "parquet-cache")]
    fn ensure_parquet_lookup(
        &mut self,
        chrom: &str,
        cache_columns: &[String],
        collect_colocated: bool,
    ) -> Result<()> {
        let cache_root = self.variation_storage.cache_root.clone();
        if self.parquet_lookup_cell.get().is_none() {
            let cache = crate::parquet_cache::detect::PartitionedParquetCache::detect(
                cache_root.to_string_lossy().as_ref(),
            )
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "parquet variation lookup but no variation manifest under {}",
                    cache_root.display()
                ))
            })?;
            let path = parquet_variation_path_for_chrom(&cache, chrom).ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "parquet variation lookup but no shard for {chrom} in {}",
                    cache_root.display()
                ))
            })?;
            let projection = cold_parquet_projection_columns(cache_columns, collect_colocated);
            let cell = Arc::clone(&self.parquet_lookup_cell);
            let open_started = Instant::now();
            let open_fut = async {
                cell.get_or_try_init(|| async {
                    SinglePathParquetVariationLookup::open(&path, projection)
                        .await
                        .map(Arc::new)
                })
                .await
                .cloned()
            };
            let _lookup = block_on(open_fut)?;
            if self.profile_enabled {
                self.profile.position_index_load += open_started.elapsed();
                self.profile.position_index_loaded += 1;
            }
        }
        Ok(())
    }

    fn record_cold_chunk_probe_metrics(
        &mut self,
        result: ColdProbeResult,
        metrics: ColdChunkProbeMetrics,
    ) {
        self.profile.cold_parquet_probes += 1;
        self.profile.cache_column_append += metrics.append_elapsed;
        self.profile.colocated_prepare += metrics.colocated_prepare_elapsed;
        self.profile.colocated_match += metrics.colocated_match_elapsed;
        self.profile.cold_parquet_rows_scanned += metrics.cold_rows_scanned;
        match result {
            ColdProbeResult::Match => self.profile.cold_parquet_matches += 1,
            ColdProbeResult::PositionCoveredNoExact => {
                self.profile.cold_parquet_position_misses += 1;
            }
            ColdProbeResult::NotCovered => self.profile.cold_parquet_not_covered += 1,
        }
        if self.profile_detailed {
            self.profile.primary_allele_rows += metrics.primary_allele_rows;
            self.profile.exact_match_calls += metrics.exact_match_calls;
            self.profile.colocated_allele_rows += metrics.colocated_allele_rows;
            self.profile.colocated_entries += metrics.colocated_entries;
            if metrics.emitted {
                self.profile.primary_matches += 1;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn enqueue_input_batch(&mut self, batch: RecordBatch) {
        let max_rows = self.input_slice_rows();
        enqueue_record_batch_slices(&mut self.pending_input_slices, batch, max_rows);
    }

    fn input_slice_rows(&self) -> usize {
        #[cfg(feature = "parquet-cache")]
        if self.warm_cold_backend.is_parquet() {
            return lookup_process_batch_rows();
        }
        usize::MAX
    }

    fn process_next_pending_input_slice(&mut self) -> Option<Result<RecordBatch>> {
        self.pending_input_slices
            .pop_front()
            .map(|batch| self.process_batch(&batch))
    }

    fn process_batch(&mut self, vcf_batch: &RecordBatch) -> Result<RecordBatch> {
        let vcf_schema = vcf_batch.schema();
        let chrom_idx = vcf_schema.index_of("chrom")?;
        let start_idx = vcf_schema.index_of("start")?;
        let end_idx = vcf_schema.index_of("end")?;
        let ref_idx = vcf_schema.index_of("ref")?;
        let alt_idx = vcf_schema.index_of("alt")?;

        let extract_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let chroms = as_string_column(vcf_batch.column(chrom_idx), "chrom")?;
        let starts = get_i32_column(vcf_batch.column(start_idx), "start")?;
        let ends = get_i32_column(vcf_batch.column(end_idx), "end")?;
        let refs = as_string_column(vcf_batch.column(ref_idx), "ref")?;
        let alts = as_string_column(vcf_batch.column(alt_idx), "alt")?;
        if let Some(t0) = extract_started {
            self.profile.extract_cols += t0.elapsed();
        }

        let num_vcf_cols = vcf_schema.fields().len();
        let num_cache_cols = self.cache_columns.len();
        let num_rows = vcf_batch.num_rows();
        if self.profile_enabled {
            self.profile.batches += 1;
            self.profile.input_rows += num_rows as u64;
        }

        // Resolve output column types from output schema.
        let output_col_types: Vec<DataType> = (0..num_cache_cols)
            .map(|i| self.schema.field(num_vcf_cols + i).data_type().clone())
            .collect();

        // Create ArrayBuilders for each cache output column.
        let mut builders: Vec<Box<dyn datafusion::arrow::array::ArrayBuilder>> = output_col_types
            .iter()
            .map(|dt| make_builder(dt, num_rows))
            .collect::<Result<Vec<_>>>()?;

        // VCF row indices for output expansion (one per output row).
        let mut vcf_indices: Vec<u32> = Vec::with_capacity(num_rows);

        // Determine which column indices in the entry correspond to our output columns.
        // Entry stores all columns except chrom/start, in schema order minus those 2.
        // (`end` is stored as a regular column inside the entry.)
        let cache_schema = &self.cache_schema;
        let cache_chrom_idx = cache_schema.index_of("chrom").unwrap_or(usize::MAX);
        let cache_start_idx = cache_schema.index_of("start").unwrap_or(usize::MAX);

        // Build mapping: output_col_positions[i] -> index within the entry's column list.
        let stored_cols: Vec<usize> = (0..cache_schema.fields().len())
            .filter(|&i| i != cache_chrom_idx && i != cache_start_idx)
            .collect();

        // Find end column for interval overlap filtering in the main match loop.
        let end_stored_col_idx: Option<usize> = cache_schema
            .index_of("end")
            .ok()
            .and_then(|schema_idx| stored_cols.iter().position(|&c| c == schema_idx));
        let failed_stored_col_idx: Option<usize> = cache_schema
            .index_of("failed")
            .ok()
            .and_then(|schema_idx| stored_cols.iter().position(|&c| c == schema_idx));

        let col_map: Vec<usize> = self
            .output_col_positions
            .iter()
            .map(|&pos| {
                stored_cols
                    .iter()
                    .position(|&c| c == pos)
                    .unwrap_or(usize::MAX)
            })
            .collect();
        let cache_columns = self.cache_columns.clone();

        let match_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };

        // Local buffer for colocated data (flushed to the shared sink after the loop).
        let mut coloc_buf: Option<HashMap<ColocatedKey, ColocatedSinkValue>> =
            if self.colocated_sink.is_some() {
                Some(HashMap::new())
            } else {
                None
            };
        let mut row_output_emitted = vec![false; num_rows];
        let mut pending_cold_probes: Vec<PendingColdProbe> = Vec::new();
        let mut pending_cold_by_row: Vec<Vec<usize>> = vec![Vec::new(); num_rows];

        for row in 0..num_rows {
            let raw_chrom = chroms.value_or_empty(row);
            let chrom = if self.vcf_has_chr {
                raw_chrom.strip_prefix("chr").unwrap_or(raw_chrom)
            } else {
                raw_chrom
            };

            let vcf_start = starts[row];
            let vcf_end = ends[row];
            let (norm_start, norm_end) = normalize_vcf_coords(
                vcf_start,
                vcf_end,
                self.vcf_zero_based,
                self.cache_zero_based,
            )?;
            let norm_start_i64 = i64::from(norm_start);
            let norm_end_i64 = i64::from(norm_end);

            let vcf_ref = refs.value_or_empty(row);
            let vcf_alt = alts.value_or_empty(row);

            let chrom_code = chrom_to_code(chrom);

            // Probe a small set of start positions used by VEP-style caches.
            // All variants at a given (chrom, start) are in one entry, so we
            // only need to probe distinct start values.
            let probe_build_started = self.profile_detailed.then(Instant::now);
            let probe_starts = build_probe_starts(
                norm_start_i64,
                norm_end_i64,
                vcf_ref,
                vcf_alt,
                self.extended_probes,
            );
            if let Some(t0) = probe_build_started {
                self.profile.probe_build += t0.elapsed();
            }

            let mut emitted_match = false;
            let vcf_iv_start = norm_start_i64.min(norm_end_i64);
            let vcf_iv_end = norm_start_i64.max(norm_end_i64);
            for probe_start in &probe_starts {
                if self.profile_enabled {
                    self.profile.probes += 1;
                }
                if let Some(warm_decision) = self.probe_warm_position(
                    chrom,
                    chrom_code,
                    *probe_start,
                    vcf_ref,
                    vcf_alt,
                    vcf_iv_start,
                    vcf_iv_end,
                    !emitted_match,
                    row as u32,
                    &cache_columns,
                    &col_map,
                    &mut builders,
                    &mut vcf_indices,
                    coloc_buf.as_mut(),
                )? {
                    match warm_decision {
                        LookupDecision::EmitWarmExact => {
                            emitted_match = true;
                            row_output_emitted[row] = true;
                            continue;
                        }
                        LookupDecision::SkipFjall => {
                            continue;
                        }
                        LookupDecision::UseFjall => {
                            let collect_colocated = coloc_buf.is_some();
                            self.ensure_parquet_lookup(chrom, &cache_columns, collect_colocated)?;
                            let pending_idx = pending_cold_probes.len();
                            pending_cold_probes.push(PendingColdProbe {
                                chrom: chrom.to_string(),
                                probe_start: *probe_start,
                                position_key: 0,
                                vcf_ref: vcf_ref.to_string(),
                                vcf_alt: vcf_alt.to_string(),
                                vcf_iv_start,
                                vcf_iv_end,
                                vcf_row: row as u32,
                            });
                            pending_cold_by_row[row].push(pending_idx);
                            continue;
                        }
                    }
                }
            }

            if !emitted_match && pending_cold_by_row[row].is_empty() {
                let null_append_started = self.profile_detailed.then(Instant::now);
                // No coordinate probe matched any allele -> null cache columns.
                vcf_indices.push(row as u32);
                for builder in &mut builders {
                    append_null_to_builder(builder.as_mut())?;
                }
                row_output_emitted[row] = true;
                if let Some(t0) = null_append_started {
                    self.profile.null_append += t0.elapsed();
                    self.profile.null_rows += 1;
                }
            }
        }

        if !pending_cold_probes.is_empty() {
            let mut taken_by_chrom =
                HashMap::<String, (RecordBatch, HashMap<u32, Vec<u32>>)>::new();

            let mut pending_starts_by_chrom: HashMap<String, Vec<u32>> = HashMap::new();
            for pending in &pending_cold_probes {
                if let Ok(start) = u32::try_from(pending.probe_start) {
                    pending_starts_by_chrom
                        .entry(pending.chrom.clone())
                        .or_default()
                        .push(start);
                }
            }

            for (chrom, mut starts) in pending_starts_by_chrom {
                starts.sort_unstable();
                starts.dedup();
                let take_started = self.profile_enabled.then(Instant::now);
                let taken = {
                    let lookup = self.parquet_lookup_cell.get().ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "parquet variation lookup not built for {chrom}"
                        ))
                    })?;
                    // Parquet resolution is stateless (footer PageDir), so a
                    // fresh per-call cursor is sufficient.
                    let mut cursor = lookup.new_cursor();
                    block_on(lookup.resolve_and_take(&starts, &mut cursor))?
                };
                if let Some(t0) = take_started {
                    self.profile.cold_parquet_load += t0.elapsed();
                }
                let row_map = start_row_map(&taken.batch)?;
                taken_by_chrom.insert(chrom, (taken.batch, row_map));
            }

            for pending in &pending_cold_probes {
                let row = pending.vcf_row as usize;
                let emit_output = !row_output_emitted[row];
                let (result, metrics) =
                    if let Some((batch, row_map)) = taken_by_chrom.get(&pending.chrom) {
                        let rows = u32::try_from(pending.probe_start)
                            .ok()
                            .and_then(|start| row_map.get(&start))
                            .map(Vec::as_slice)
                            .unwrap_or(&[]);
                        probe_taken_batch_position(
                            batch,
                            rows,
                            self.exact_matcher,
                            self.allowed_failed,
                            self.profile_detailed,
                            &pending.chrom,
                            pending.probe_start,
                            &pending.vcf_ref,
                            &pending.vcf_alt,
                            pending.vcf_iv_start,
                            pending.vcf_iv_end,
                            emit_output,
                            pending.vcf_row,
                            &cache_columns,
                            &col_map,
                            &mut builders,
                            &mut vcf_indices,
                            coloc_buf.as_mut(),
                        )?
                    } else {
                        (
                            ColdProbeResult::NotCovered,
                            ColdChunkProbeMetrics::default(),
                        )
                    };

                if self.profile_enabled {
                    self.record_cold_chunk_probe_metrics(result, metrics);
                }
                if result == ColdProbeResult::Match && emit_output {
                    row_output_emitted[row] = true;
                }
            }

            for (row, pending_indices) in pending_cold_by_row.iter().enumerate() {
                if pending_indices.is_empty() || row_output_emitted[row] {
                    continue;
                }

                let null_append_started = self.profile_detailed.then(Instant::now);
                vcf_indices.push(row as u32);
                for builder in &mut builders {
                    append_null_to_builder(builder.as_mut())?;
                }
                row_output_emitted[row] = true;
                if let Some(t0) = null_append_started {
                    self.profile.null_append += t0.elapsed();
                    self.profile.null_rows += 1;
                }
            }
        }

        // Flush colocated data to the shared sink.
        let colocated_flush_started = self.profile_detailed.then(Instant::now);
        if let (Some(buf), Some(sink)) = (coloc_buf, &self.colocated_sink) {
            if !buf.is_empty() {
                let mut guard = sink.lock().map_err(|e| {
                    DataFusionError::Execution(format!("colocated sink mutex poisoned: {e}"))
                })?;
                for (key, mut value) in buf {
                    guard
                        .entry(key)
                        .and_modify(|existing| {
                            if existing.compare_output_allele.is_none() {
                                existing.compare_output_allele =
                                    value.compare_output_allele.clone();
                            }
                            if existing.unshifted_output_allele.is_none() {
                                existing.unshifted_output_allele =
                                    value.unshifted_output_allele.clone();
                            }
                            existing.entries.append(&mut value.entries);
                        })
                        .or_insert(value);
                }
            }
        }
        if let Some(t0) = colocated_flush_started {
            self.profile.colocated_flush += t0.elapsed();
        }

        if let Some(t0) = match_started {
            self.profile.match_loop += t0.elapsed();
        }

        if self.profile_enabled {
            self.profile.output_rows += vcf_indices.len() as u64;
        }

        // Take VCF columns using expanded indices.
        let take_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let output_order = output_order_for_vcf_indices(&vcf_indices);
        let ordered_vcf_indices = if let Some(output_order) = output_order.as_ref() {
            output_order
                .iter()
                .map(|&idx| vcf_indices[idx as usize])
                .collect::<Vec<_>>()
        } else {
            vcf_indices
        };
        let cache_reorder_indices = output_order
            .as_ref()
            .map(|output_order| UInt32Array::from(output_order.clone()));
        let take_indices = UInt32Array::from(ordered_vcf_indices);
        let mut output_columns: Vec<ArrayRef> = Vec::with_capacity(num_vcf_cols + num_cache_cols);
        for col_idx in 0..num_vcf_cols {
            let taken =
                datafusion::arrow::compute::take(vcf_batch.column(col_idx), &take_indices, None)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            output_columns.push(taken);
        }
        if let Some(t0) = take_started {
            self.profile.vcf_take += t0.elapsed();
        }

        // Finish builders -> cache output columns.
        let cache_build_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        for builder in &mut builders {
            let array = builder.finish();
            if let Some(cache_reorder_indices) = cache_reorder_indices.as_ref() {
                let reordered =
                    datafusion::arrow::compute::take(array.as_ref(), cache_reorder_indices, None)
                        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                output_columns.push(reordered);
            } else {
                output_columns.push(array);
            }
        }
        if let Some(t0) = cache_build_started {
            self.profile.cache_build += t0.elapsed();
        }

        RecordBatch::try_new(self.schema.clone(), output_columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

fn output_order_for_vcf_indices(vcf_indices: &[u32]) -> Option<Vec<u32>> {
    if vcf_indices.windows(2).all(|window| window[0] <= window[1]) {
        return None;
    }

    let mut output_order: Vec<u32> = (0..vcf_indices.len() as u32).collect();
    output_order.sort_by_key(|&idx| (vcf_indices[idx as usize], idx));
    Some(output_order)
}

#[cfg(feature = "parquet-cache")]
fn batch_output_indices(
    batch: &RecordBatch,
    cache_columns: &[String],
    col_map: &[usize],
) -> Option<Vec<Option<usize>>> {
    cache_columns
        .iter()
        .zip(col_map.iter())
        .map(|(name, entry_idx)| {
            if *entry_idx == usize::MAX {
                Some(None)
            } else {
                batch.schema().index_of(name).ok().map(Some)
            }
        })
        .collect()
}

#[cfg(feature = "parquet-cache")]
fn append_batch_row_values(
    batch: &RecordBatch,
    row: usize,
    output_indices: &[Option<usize>],
    builders: &mut [Box<dyn datafusion::arrow::array::ArrayBuilder>],
) -> Result<()> {
    for (idx, builder) in output_indices.iter().zip(builders.iter_mut()) {
        if let Some(idx) = idx {
            append_array_value_to_builder(batch.column(*idx).as_ref(), row, builder.as_mut())?;
        } else {
            append_null_to_builder(builder.as_mut())?;
        }
    }
    Ok(())
}

#[cfg(feature = "parquet-cache")]
fn batch_string_value(
    batch: &RecordBatch,
    column_idx: Option<usize>,
    row: usize,
) -> Result<Option<String>> {
    let Some(column_idx) = column_idx else {
        return Ok(None);
    };
    let array = batch.column(column_idx);
    if row >= array.len() || array.is_null(row) {
        return Ok(None);
    }

    if let Some(array) = array.as_any().downcast_ref::<StringArray>() {
        Ok(Some(array.value(row).to_string()))
    } else if let Some(array) = array.as_any().downcast_ref::<StringViewArray>() {
        Ok(Some(array.value(row).to_string()))
    } else if let Some(array) = array.as_any().downcast_ref::<LargeStringArray>() {
        Ok(Some(array.value(row).to_string()))
    } else {
        Err(DataFusionError::Execution(format!(
            "Lance batch column expected string array, got {:?}",
            array.data_type()
        )))
    }
}

#[cfg(feature = "parquet-cache")]
fn batch_i64_value(batch: &RecordBatch, column_idx: Option<usize>, row: usize) -> Option<i64> {
    let array = batch.column(column_idx?);
    if row >= array.len() || array.is_null(row) {
        return None;
    }

    if let Some(array) = array.as_any().downcast_ref::<Int8Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<Int16Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<Int32Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        Some(array.value(row))
    } else if let Some(array) = array.as_any().downcast_ref::<UInt8Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<UInt16Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<UInt32Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<BooleanArray>() {
        // Parquet variation cache stores the binary flags (failed/somatic/
        // phenotype_or_disease) as presence Boolean; read true->1, false->0 so
        // they compare identically to the Int8-encoded Parquet cache.
        Some(array.value(row) as i64)
    } else {
        array
            .as_any()
            .downcast_ref::<UInt64Array>()
            .map(|array| array.value(row) as i64)
    }
}

#[cfg(all(test, feature = "parquet-cache"))]
mod batch_i64_boolean_tests {
    use super::*;
    use datafusion::arrow::array::BooleanArray;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    #[test]
    fn batch_i64_value_reads_boolean_presence_flags() {
        let arr = BooleanArray::from(vec![Some(true), Some(false), None]);
        let schema = Arc::new(Schema::new(vec![Field::new(
            "failed",
            DataType::Boolean,
            true,
        )]));
        let batch = RecordBatch::try_new(schema, vec![Arc::new(arr)]).unwrap();
        assert_eq!(batch_i64_value(&batch, Some(0), 0), Some(1)); // true  -> 1
        assert_eq!(batch_i64_value(&batch, Some(0), 1), Some(0)); // false -> 0
        assert_eq!(batch_i64_value(&batch, Some(0), 2), None); //    null  -> None
    }
}

#[cfg(feature = "parquet-cache")]
fn resolve_batch_coloc_indices(batch: &RecordBatch) -> Option<WarmColocIndices> {
    let schema = batch.schema();
    let find = |name: &str| schema.index_of(name).ok();
    Some(WarmColocIndices {
        variation_name: find("variation_name")?,
        end_col: find("end"),
        failed: find("failed"),
        somatic: find("somatic"),
        pheno: find("phenotype_or_disease"),
        clin_sig: find("clin_sig"),
        clin_sig_allele: find("clin_sig_allele"),
        pubmed: find("pubmed"),
        af_indices: AF_COL_NAMES.iter().map(|name| find(name)).collect(),
    })
}

fn append_array_value_to_builder(
    array: &dyn Array,
    row: usize,
    builder: &mut dyn datafusion::arrow::array::ArrayBuilder,
) -> Result<()> {
    use datafusion::arrow::array::*;

    if row >= array.len() || array.is_null(row) {
        return append_null_to_builder(builder);
    }

    match array.data_type() {
        DataType::Utf8 => {
            let array = array
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    DataFusionError::Execution("warm Utf8 column is not StringArray".into())
                })?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<StringBuilder>()
                .ok_or_else(|| DataFusionError::Execution("expected StringBuilder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Utf8View => {
            let array = array
                .as_any()
                .downcast_ref::<StringViewArray>()
                .ok_or_else(|| {
                    DataFusionError::Execution("warm Utf8View column is not StringViewArray".into())
                })?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<StringBuilder>()
                .ok_or_else(|| DataFusionError::Execution("expected StringBuilder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::LargeUtf8 => {
            let array = array
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "warm LargeUtf8 column is not LargeStringArray".into(),
                    )
                })?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<StringBuilder>()
                .ok_or_else(|| DataFusionError::Execution("expected StringBuilder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Int8 => {
            let array = array
                .as_any()
                .downcast_ref::<Int8Array>()
                .ok_or_else(|| DataFusionError::Execution("warm Int8 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<Int8Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected Int8Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Int16 => {
            let array = array
                .as_any()
                .downcast_ref::<Int16Array>()
                .ok_or_else(|| DataFusionError::Execution("warm Int16 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<Int16Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected Int16Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Int32 => {
            let array = array
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or_else(|| DataFusionError::Execution("warm Int32 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<Int32Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected Int32Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Int64 => {
            let array = array
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| DataFusionError::Execution("warm Int64 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<Int64Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected Int64Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::UInt8 => {
            let array = array
                .as_any()
                .downcast_ref::<UInt8Array>()
                .ok_or_else(|| DataFusionError::Execution("warm UInt8 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<UInt8Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected UInt8Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::UInt16 => {
            let array = array
                .as_any()
                .downcast_ref::<UInt16Array>()
                .ok_or_else(|| DataFusionError::Execution("warm UInt16 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<UInt16Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected UInt16Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::UInt32 => {
            let array = array
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("warm UInt32 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<UInt32Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected UInt32Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::UInt64 => {
            let array = array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| DataFusionError::Execution("warm UInt64 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<UInt64Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected UInt64Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Float32 => {
            let array = array
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| DataFusionError::Execution("warm Float32 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<Float32Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected Float32Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Float64 => {
            let array = array
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| DataFusionError::Execution("warm Float64 column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<Float64Builder>()
                .ok_or_else(|| DataFusionError::Execution("expected Float64Builder".into()))?;
            builder.append_value(array.value(row));
        }
        DataType::Boolean => {
            let array = array
                .as_any()
                .downcast_ref::<BooleanArray>()
                .ok_or_else(|| DataFusionError::Execution("warm Boolean column mismatch".into()))?;
            let builder = builder
                .as_any_mut()
                .downcast_mut::<BooleanBuilder>()
                .ok_or_else(|| DataFusionError::Execution("expected BooleanBuilder".into()))?;
            builder.append_value(array.value(row));
        }
        other => {
            return Err(DataFusionError::Execution(format!(
                "warm direct append: unsupported column type {other}"
            )));
        }
    }

    Ok(())
}

/// Append a single null value to any supported ArrayBuilder.
fn append_null_to_builder(builder: &mut dyn datafusion::arrow::array::ArrayBuilder) -> Result<()> {
    use datafusion::arrow::array::*;

    if let Some(b) = builder.as_any_mut().downcast_mut::<StringBuilder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<Int32Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<Int64Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<Float32Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<Float64Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<UInt32Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<UInt64Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<BooleanBuilder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<Int8Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<Int16Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<UInt8Builder>() {
        b.append_null();
    } else if let Some(b) = builder.as_any_mut().downcast_mut::<UInt16Builder>() {
        b.append_null();
    } else {
        return Err(DataFusionError::Execution(
            "unsupported builder type for null append".into(),
        ));
    }
    Ok(())
}

fn normalize_cache_output_type(data_type: &DataType) -> DataType {
    match data_type {
        DataType::Utf8View | DataType::LargeUtf8 => DataType::Utf8,
        other => other.clone(),
    }
}

fn get_i32_column(col: &ArrayRef, column_name: &str) -> Result<Vec<i32>> {
    if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
        return (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    Err(DataFusionError::Execution(format!(
                        "column '{column_name}' contains NULL at row {i}"
                    )))
                } else {
                    Ok(arr.value(i))
                }
            })
            .collect();
    }
    if let Some(arr) = col.as_any().downcast_ref::<UInt32Array>() {
        return (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    Err(DataFusionError::Execution(format!(
                        "column '{column_name}' contains NULL at row {i}"
                    )))
                } else {
                    i32::try_from(arr.value(i)).map_err(|_| {
                        DataFusionError::Execution(format!(
                            "column '{column_name}' value {} at row {i} overflows Int32",
                            arr.value(i)
                        ))
                    })
                }
            })
            .collect();
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    Err(DataFusionError::Execution(format!(
                        "column '{column_name}' contains NULL at row {i}"
                    )))
                } else {
                    i32::try_from(arr.value(i)).map_err(|_| {
                        DataFusionError::Execution(format!(
                            "column '{column_name}' value {} at row {i} overflows Int32",
                            arr.value(i)
                        ))
                    })
                }
            })
            .collect();
    }
    if let Some(arr) = col.as_any().downcast_ref::<UInt64Array>() {
        return (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    Err(DataFusionError::Execution(format!(
                        "column '{column_name}' contains NULL at row {i}"
                    )))
                } else {
                    i32::try_from(arr.value(i)).map_err(|_| {
                        DataFusionError::Execution(format!(
                            "column '{column_name}' value {} at row {i} overflows Int32",
                            arr.value(i)
                        ))
                    })
                }
            })
            .collect();
    }
    Err(DataFusionError::Execution(format!(
        "column '{column_name}' expected Int32/UInt32/Int64/UInt64 array, got {:?}",
        col.data_type()
    )))
}

fn as_string_column<'a>(col: &'a ArrayRef, column_name: &str) -> Result<StringColumnView<'a>> {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        Ok(StringColumnView::Utf8(arr))
    } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        Ok(StringColumnView::Utf8View(arr))
    } else if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
        Ok(StringColumnView::LargeUtf8(arr))
    } else {
        Err(DataFusionError::Execution(format!(
            "column '{column_name}' expected string array, got {:?}",
            col.data_type()
        )))
    }
}

fn normalize_vcf_coords(
    start: i32,
    end: i32,
    vcf_zero_based: bool,
    cache_zero_based: bool,
) -> Result<(i32, i32)> {
    if vcf_zero_based == cache_zero_based {
        Ok((start, end))
    } else if vcf_zero_based {
        let shifted_start = start.checked_add(1).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "failed to normalize VCF coordinates: start {start} overflows Int32 during 0-based -> 1-based conversion"
            ))
        })?;
        Ok((shifted_start, end)) // 0-based half-open -> 1-based closed
    } else {
        let shifted_start = start.checked_sub(1).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "failed to normalize VCF coordinates: start {start} underflows Int32 during 1-based -> 0-based conversion"
            ))
        })?;
        Ok((shifted_start, end)) // 1-based closed -> 0-based half-open
    }
}

#[inline]
fn push_unique_probe_start(probe_starts: &mut Vec<i64>, start: i64) {
    if !probe_starts.contains(&start) {
        probe_starts.push(start);
    }
}

fn build_probe_starts(
    norm_start_i64: i64,
    norm_end_i64: i64,
    vcf_ref: &str,
    vcf_alt: &str,
    extended_probes: bool,
) -> Vec<i64> {
    let mut probe_starts: Vec<i64> = Vec::with_capacity(6);
    push_unique_probe_start(&mut probe_starts, norm_start_i64);

    if !extended_probes {
        return probe_starts;
    }

    if norm_start_i64 == norm_end_i64 {
        push_unique_probe_start(&mut probe_starts, norm_start_i64.saturating_add(1));
    } else {
        push_unique_probe_start(&mut probe_starts, norm_end_i64);
    }

    for alt in vcf_alt.split(['|', ',']).filter(|a| !a.is_empty()) {
        let (_, _, input_start) = vcf_to_vep_input_allele(norm_start_i64, vcf_ref, alt);
        push_unique_probe_start(&mut probe_starts, input_start);

        let shift_usize = common_prefix_len(vcf_ref, alt);
        if shift_usize == 0 {
            continue;
        }
        let shift = shift_usize as i64;
        if let Some(shifted_start) = norm_start_i64.checked_add(shift) {
            push_unique_probe_start(&mut probe_starts, shifted_start);
        }
    }

    // Deletions in tandem repeats may be right/left shifted in cache coordinates.
    // Probe a bounded window of equivalent deletion start positions.
    for alt in vcf_alt.split(['|', ',']).filter(|a| !a.is_empty()) {
        let (ref_event_len, alt_event_len) = canonical_event_lengths(vcf_ref, alt);
        if ref_event_len == 0 || alt_event_len != 0 {
            continue;
        }
        let del_len = ref_event_len as i64;
        let max_shift = del_len.min(32);
        for base_start in [norm_start_i64, norm_start_i64.saturating_sub(1)] {
            for shift in 0..=max_shift {
                let Some(candidate_start) = base_start.checked_add(shift) else {
                    continue;
                };
                let Some(candidate_end) = candidate_start.checked_add(del_len - 1) else {
                    continue;
                };
                // Mirror SQL interval semantics: only consider intervals overlapping
                // the normalized VCF interval.
                if candidate_start > norm_end_i64 || candidate_end < norm_start_i64 {
                    continue;
                }
                push_unique_probe_start(&mut probe_starts, candidate_start);
            }
        }
    }

    probe_starts
}

#[inline]
fn probe_start_visible_to_window(probe_start: i64, compare_start: i64, compare_end: i64) -> bool {
    let vis_start = (compare_start - 1).min(compare_end + 1);
    let vis_end = (compare_start - 1).max(compare_end + 1);
    probe_start >= vis_start && probe_start <= vis_end
}

#[inline]
fn common_prefix_len(left: &str, right: &str) -> usize {
    left.as_bytes()
        .iter()
        .zip(right.as_bytes().iter())
        .take_while(|(a, b)| a == b)
        .count()
}

#[inline]
fn canonical_event_lengths(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    let mut ref_start = 0usize;
    let mut alt_start = 0usize;
    while ref_start < ref_bytes.len()
        && alt_start < alt_bytes.len()
        && ref_bytes[ref_start] == alt_bytes[alt_start]
    {
        ref_start += 1;
        alt_start += 1;
    }

    let mut ref_end = ref_bytes.len();
    let mut alt_end = alt_bytes.len();
    while ref_end > ref_start
        && alt_end > alt_start
        && ref_bytes[ref_end - 1] == alt_bytes[alt_end - 1]
    {
        ref_end -= 1;
        alt_end -= 1;
    }

    (ref_end - ref_start, alt_end - alt_start)
}

/// Resolve column indices within the KV entry for co-located fields.
///
/// The entry stores all cache schema columns except chrom/start, in schema order
/// minus those two. `end` is stored as a regular column inside the entry.
impl Stream for KvLookupStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        if let Some(result) = self.process_next_pending_input_slice() {
            return Poll::Ready(Some(result));
        }

        // KvLookupExec reads ALL alleles at each position in a single point
        // lookup (the position entry contains all alleles), so co-located
        // data for each VCF row is complete immediately — no buffering needed.
        match self.input.poll_next_unpin(cx) {
            Poll::Ready(Some(Ok(batch))) => {
                self.enqueue_input_batch(batch);
                Poll::Ready(self.process_next_pending_input_slice())
            }
            Poll::Ready(Some(Err(e))) => Poll::Ready(Some(Err(e))),
            Poll::Ready(None) => {
                if self.profile_enabled && !self.profile_emitted {
                    self.profile.emit(self.profile_detailed);
                    self.profile_emitted = true;
                }
                Poll::Ready(None)
            }
            Poll::Pending => Poll::Pending,
        }
    }
}

impl RecordBatchStream for KvLookupStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

struct EmptyLookupInput {
    schema: SchemaRef,
}

impl Stream for EmptyLookupInput {
    type Item = Result<RecordBatch>;

    fn poll_next(self: Pin<&mut Self>, _cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        Poll::Ready(None)
    }
}

impl RecordBatchStream for EmptyLookupInput {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_order_restores_vcf_row_order_after_grouped_cold_replay() {
        assert_eq!(output_order_for_vcf_indices(&[0, 1, 2]), None);
        assert_eq!(
            output_order_for_vcf_indices(&[2, 0, 1]),
            Some(vec![1, 2, 0])
        );
    }

    #[test]
    fn build_probe_starts_includes_parser_input_start_for_shifted_insertions() {
        let probe_starts = build_probe_starts(
            215230092,
            215230102,
            "TACACACACAC",
            "TATACACACACACACAC",
            true,
        );

        assert_eq!(probe_starts[0], 215230092);
        assert!(probe_starts.contains(&215230093));
        assert!(probe_starts.contains(&215230094));
    }

    #[test]
    fn build_probe_starts_includes_parser_input_start_for_repeat_insertions() {
        let probe_starts = build_probe_starts(165387539, 165387541, "CTG", "CTCTGTG", true);

        assert_eq!(probe_starts[0], 165387539);
        assert!(probe_starts.contains(&165387540));
        assert!(probe_starts.contains(&165387541));
    }

    #[test]
    fn colocated_visibility_uses_vep_minimized_probe_window() {
        assert!(probe_start_visible_to_window(602114, 602114, 602113));
    }

    #[test]
    fn colocated_visibility_rejects_outside_vep_minimized_probe_window() {
        assert!(!probe_start_visible_to_window(602117, 602114, 602113));
    }

    #[test]
    fn lookup_profile_detail_line_formats_probe_decode_match_breakdown() {
        let mut profile = LookupProfile::default();
        profile.probe_build += Duration::from_millis(1);
        profile.point_get_raw += Duration::from_millis(2);
        profile.prefetch_map_lookup += Duration::from_millis(3);
        profile.decompress += Duration::from_millis(4);
        profile.reader_init += Duration::from_millis(5);
        profile.primary_match += Duration::from_millis(6);
        profile.cache_column_append += Duration::from_millis(7);
        profile.colocated_prepare += Duration::from_millis(8);
        profile.colocated_match += Duration::from_millis(9);
        profile.colocated_flush += Duration::from_millis(10);
        profile.null_append += Duration::from_millis(11);
        profile.warm_probe += Duration::from_millis(12);
        profile.warm_chunk_load += Duration::from_millis(13);
        profile.raw_get_hits = 12;
        profile.raw_get_misses = 13;
        profile.prefetch_hits = 14;
        profile.prefetch_misses = 15;
        profile.decode_calls = 16;
        profile.compressed_bytes = 17;
        profile.decompressed_bytes = 18;
        profile.primary_allele_rows = 19;
        profile.primary_failed_skips = 20;
        profile.primary_interval_skips = 21;
        profile.exact_match_calls = 22;
        profile.primary_matches = 23;
        profile.colocated_allele_rows = 24;
        profile.colocated_entries = 25;
        profile.null_rows = 26;
        profile.warm_probes = 27;
        profile.warm_matches = 28;
        profile.warm_definitive_misses = 29;
        profile.warm_not_covered = 30;
        profile.warm_boundary_fallbacks = 31;
        profile.warm_candidate_rows = 32;
        profile.warm_rows_scanned = 33;
        profile.warm_chunks_loaded = 34;
        profile.warm_chunk_rows = 35;
        profile.position_index_persisted_loads = 1;
        profile.variant_bloom_checks = 2;
        profile.variant_bloom_negative_skips = 3;
        profile.variant_bloom_positive_gets = 4;
        profile.variant_bloom_loaded = 5;
        profile.variant_bloom_entries = 6;
        profile.variant_bloom_bits = 7;
        profile.variant_bloom_hashes = 8;
        profile.variant_bloom_bytes = 9;
        profile.cold_parquet_probes = 36;
        profile.cold_parquet_matches = 37;
        profile.cold_parquet_rows_scanned = 38;
        profile.cold_parquet_row_groups_total = 50;
        profile.cold_parquet_row_groups_unique_touched = 40;
        profile.cold_parquet_row_group_metadata_probes = 41;
        profile.cold_parquet_row_group_current_hits = 42;
        profile.cold_parquet_row_group_previous_hits = 43;
        profile.cold_parquet_row_group_advanced_hits = 44;
        profile.cold_parquet_row_group_binary_search_hits = 45;
        profile.cold_parquet_row_group_metadata_misses = 46;
        profile.cold_parquet_row_group_skipped_ahead = 47;
        profile.cold_parquet_row_groups_loaded = 39;
        profile.cold_parquet_row_group_load_batches = 12;
        profile.cold_parquet_position_page_index_loaded = true;
        profile.cold_parquet_position_column_index_loaded = true;
        profile.cold_parquet_position_pages_total = 51;
        profile.cold_parquet_position_bloom_filter_row_groups = 52;
        profile.cold_parquet_position_bloom_filter_bytes = 53;
        profile.cold_parquet_page_index_probes = 54;
        profile.cold_parquet_page_index_available_probes = 55;
        profile.cold_parquet_page_index_unavailable_probes = 56;
        profile.cold_parquet_page_index_pages_in_probed_row_groups = 57;
        profile.cold_parquet_page_index_candidate_pages = 58;
        profile.cold_parquet_page_index_candidate_rows = 59;
        profile.cold_parquet_page_index_unique_candidate_pages = 60;
        profile.cold_parquet_page_index_unique_candidate_rows = 61;
        profile.cold_parquet_page_index_candidate_misses = 62;

        let lines = profile.detail_lines();

        assert_eq!(lines.len(), 9);
        assert!(lines[0].contains("probe_build=0.001s"));
        assert!(lines[0].contains("point_get_raw=0.002s"));
        assert!(lines[0].contains("decompress=0.004s"));
        assert!(lines[0].contains("warm_probe=0.012s"));
        assert!(lines[0].contains("warm_chunk_load=0.013s"));
        assert!(lines[0].contains("variant_bloom_load=0.000s"));
        assert!(lines[1].contains("raw_get_hits=12"));
        assert!(lines[1].contains("compressed_bytes=17"));
        assert!(lines[2].contains("primary_allele_rows=19"));
        assert!(lines[2].contains("colocated_entries=25"));
        assert!(lines[2].contains("null_rows=26"));
        assert!(lines[3].contains("probes=27"));
        assert!(lines[3].contains("definitive_misses=29"));
        assert!(lines[3].contains("candidate_rows=32"));
        assert!(lines[3].contains("rows_scanned=33"));
        assert!(lines[3].contains("chunk_rows=35"));
        assert!(lines[4].contains("position_index checks=0"));
        assert!(lines[4].contains("persisted_loads=1"));
        assert!(lines[5].contains("variant_bloom checks=2"));
        assert!(lines[5].contains("negative_skips=3"));
        assert!(lines[5].contains("positive_gets=4"));
        assert!(lines[5].contains("entries=6"));
        assert!(lines[6].contains("cold_tier probes=36"));
        assert!(lines[6].contains("matches=37"));
        assert!(lines[6].contains("rows_scanned=38"));
        assert!(lines[7].contains("cold_parquet_row_groups total=50"));
        assert!(lines[7].contains("unique_touched=40"));
        assert!(lines[7].contains("untouched=10"));
        assert!(lines[7].contains("metadata_probes=41"));
        assert!(lines[7].contains("loaded=39"));
        assert!(lines[7].contains("load_batches=12"));
        assert!(lines[8].contains("cold_parquet_pages position_page_index_loaded=true"));
        assert!(lines[8].contains("position_pages_total=51"));
        assert!(lines[8].contains("position_bloom_filter_row_groups=52"));
        assert!(lines[8].contains("page_index_probes=54"));
        assert!(lines[8].contains("unique_candidate_rows=61"));
    }

    #[cfg(feature = "parquet-cache")]
    #[test]
    fn variation_lookup_profile_line_identifies_parquet_backend() {
        let line = format_variation_lookup_profile_line(
            "parquet",
            "parquet",
            WarmColdVariationBackend::Parquet,
            WarmColdVariationIndexMode::PositionThenVariantBloom,
            "/cache/variation",
            None,
            None,
            Some("/cache/variation.position_index"),
            Some("/cache/variation.variant_bloom_index"),
        );

        assert!(line.contains("storage=parquet"));
        assert!(line.contains("warm_source=parquet"));
        assert!(line.contains("warm_cold_backend=parquet"));
        assert!(line.contains("index_mode=posidx_bloom"));
        assert!(line.contains("warm_dir=-"));
        assert!(line.contains("cold_dir=-"));
        assert!(line.contains("variant_bloom_index_dir=/cache/variation.variant_bloom_index"));
    }

    #[cfg(feature = "parquet-cache")]
    #[test]
    fn warm_source_label_is_parquet() {
        assert_eq!(
            warm_source_label(WarmColdVariationBackend::Parquet),
            "parquet"
        );
    }

    #[test]
    fn record_batch_slicer_bounds_large_lookup_batches() {
        let schema = Arc::new(datafusion::arrow::datatypes::Schema::new(vec![
            datafusion::arrow::datatypes::Field::new("start", DataType::Int32, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![Arc::new(Int32Array::from_iter_values(0..12_001)) as ArrayRef],
        )
        .unwrap();
        let mut queue = VecDeque::new();

        enqueue_record_batch_slices(&mut queue, batch, 5_000);

        let lengths = queue.iter().map(RecordBatch::num_rows).collect::<Vec<_>>();
        assert_eq!(lengths, vec![5_000, 5_000, 2_001]);
    }
    #[cfg(feature = "parquet-cache")]
    #[test]
    fn colocated_probe_matches_chr1_homopolymer_deletion() {
        use datafusion::arrow::datatypes::{Field, Schema};
        use std::collections::HashMap as StdHashMap;
        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, false),
            Field::new("gnomADg", DataType::Utf8, true),
            Field::new("gnomADg_EAS", DataType::Utf8, true),
        ]));
        let allele_string =
            "AAAAAAAAAAAAAAAA/AAAAAAAAAAAA/AAAAAAAAAAAAA/AAAAAAAAAAAAAA/AAAAAAAAAAAAAAA";
        let batch = RecordBatch::try_new(
            cache_schema,
            vec![
                Arc::new(UInt32Array::from(vec![244_978_492])) as ArrayRef,
                Arc::new(UInt32Array::from(vec![244_978_507])) as ArrayRef,
                Arc::new(StringArray::from(vec![Some("rs58680543")])) as ArrayRef,
                Arc::new(StringArray::from(vec![allele_string])) as ArrayRef,
                Arc::new(Int8Array::from(vec![0])) as ArrayRef,
                Arc::new(StringArray::from(vec![Some("AAAAAAAAAAAAAAA:0.1017")])) as ArrayRef,
                Arc::new(StringArray::from(vec![Some("AAAAAAAAAAAAAAA:0.2402")])) as ArrayRef,
            ],
        )
        .unwrap();

        let mut coloc = StdHashMap::new();
        let mut builders = Vec::new();
        let mut vcf_indices = Vec::new();
        let row_map = start_row_map(&batch).unwrap();
        let rows = row_map.get(&244_978_492).unwrap();

        let (result, metrics) = probe_taken_batch_position(
            &batch,
            rows,
            allele_matches as fn(&str, &str, &str) -> bool,
            0,
            true,
            "1",
            244_978_492,
            "CA",
            "C",
            244_978_491,
            244_978_492,
            false,
            0,
            &[],
            &[],
            &mut builders,
            &mut vcf_indices,
            Some(&mut coloc),
        )
        .unwrap();

        assert_eq!(result, ColdProbeResult::Match);
        assert_eq!(metrics.colocated_entries, 1);
        let key = ("1".to_string(), 244_978_492, 244_978_492, "A/-".to_string());
        let sink_value = coloc
            .get(&key)
            .expect("colocated entry keyed by parser allele");
        assert_eq!(sink_value.entries.len(), 1);
        assert_eq!(sink_value.entries[0].variation_name, "rs58680543");
        assert_eq!(sink_value.entries[0].af_value(16), "AAAAAAAAAAAAAAA:0.1017");
        assert_eq!(sink_value.entries[0].af_value(21), "AAAAAAAAAAAAAAA:0.2402");
    }
}
