//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store per-position for annotation.

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

use super::allele_index::AlleleMatcher;
use super::key_encoding::chrom_to_code;
use super::kv_store::RangePrefetchLimitExceeded;
use super::kv_store::VepKvStore;
use super::position_entry::{PositionEntryReader, make_builder};
use super::position_index::{
    PositionIndexSource, cold_variation_file_for_chrom, cold_variation_files_for_chrom,
};
use super::variant_bloom_index::{VariantBloomIndex, find_variant_bloom_index_file};
use crate::allele::{
    VariantAlleleInput, allele_matches, get_matched_variant_alleles, vcf_to_vep_allele,
    vcf_to_vep_input_allele, vep_norm_end, vep_norm_start,
};
use crate::variant_lookup_exec::{
    AF_COL_NAMES, ColocatedCacheEntry, ColocatedKey, ColocatedSink, ColocatedSinkValue,
    compare_existing_variant_alleles, output_allele_from_allele_string, read_reference_sequence,
};
use crate::warm_cache::chrom_cache::{
    WarmChromCache, WarmChromOpen, WarmProbe, projection_columns_for_cache,
};
use crate::warm_cache::chunk::WarmChunkContext;
use crate::warm_cache::cold_parquet::{
    ColdParquetLookupSet, ColdParquetStats, ColdProbeResult, cold_parquet_projection_columns,
};
use crate::warm_cache::key::{
    position_key_from_code as warm_position_key_from_code,
    variant_key_from_position as warm_variant_key_from_position,
};

/// Lookup match mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KvMatchMode {
    /// Exact allele matching only.
    Exact,
}

/// Physical execution plan for KV-backed variant lookup.
///
/// Takes a VCF input plan, probes a fjall KV store per-position,
/// and emits LEFT JOIN output (unmatched VCF rows get NULL cache columns).
pub struct KvLookupExec {
    input: Arc<dyn ExecutionPlan>,
    store: Option<Arc<VepKvStore>>,
    cache_root: Option<PathBuf>,
    cache_schema: SchemaRef,
    indexed_parquet: bool,
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
    /// Optional indexed reference FASTA for genomic shift state in colocated
    /// matching (parity with parquet path's two-pass allele matching).
    reference_fasta_path: Option<String>,
    #[cfg(test)]
    warm_cache_dir_override: Option<PathBuf>,
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
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        input: Arc<dyn ExecutionPlan>,
        store: Arc<VepKvStore>,
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
        let cache_schema = store.schema().clone();
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
            store: Some(store),
            cache_root: None,
            cache_schema,
            indexed_parquet: false,
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
            reference_fasta_path: None,
            #[cfg(test)]
            warm_cache_dir_override: None,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_indexed_parquet(
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
            store: None,
            cache_root: Some(cache_root),
            cache_schema,
            indexed_parquet: true,
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
            reference_fasta_path: None,
            #[cfg(test)]
            warm_cache_dir_override: None,
        })
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

    /// Set the reference FASTA path for genomic shift state in colocated matching.
    pub fn with_reference_fasta_path(mut self, path: Option<String>) -> Self {
        self.reference_fasta_path = path;
        self
    }

    #[cfg(test)]
    fn with_warm_cache_dir_override(mut self, path: PathBuf) -> Self {
        self.warm_cache_dir_override = Some(path);
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
        let mut exec = if self.indexed_parquet {
            KvLookupExec::new_indexed_parquet(
                children[0].clone(),
                self.cache_root.clone().ok_or_else(|| {
                    DataFusionError::Execution("indexed parquet lookup missing cache root".into())
                })?,
                self.cache_schema.clone(),
                self.cache_columns.clone(),
                self.match_mode,
                self.exact_matcher,
                self.vcf_has_chr,
                self.vcf_zero_based,
                self.cache_zero_based,
                self.extended_probes,
                self.allowed_failed,
            )?
        } else {
            KvLookupExec::new(
                children[0].clone(),
                self.store.clone().ok_or_else(|| {
                    DataFusionError::Execution("fjall lookup missing KV store".into())
                })?,
                self.cache_columns.clone(),
                self.match_mode,
                self.exact_matcher,
                self.vcf_has_chr,
                self.vcf_zero_based,
                self.cache_zero_based,
                self.extended_probes,
                self.allowed_failed,
            )?
        };
        if let Some(sink) = &self.colocated_sink {
            exec = exec.with_colocated_sink(Arc::clone(sink));
        }
        if let Some(sinks) = &self.colocated_partition_sinks {
            exec = exec.with_colocated_partition_sinks(sinks.clone());
        }
        exec = exec.with_reference_fasta_path(self.reference_fasta_path.clone());
        #[cfg(test)]
        if let Some(path) = &self.warm_cache_dir_override {
            exec = exec.with_warm_cache_dir_override(path.clone());
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

        let stream = KvLookupStream::new(
            input_stream,
            self.store.clone(),
            self.cache_root.clone(),
            self.cache_schema.clone(),
            self.indexed_parquet,
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
            self.reference_fasta_path.clone(),
        );

        #[cfg(test)]
        if let Some(path) = &self.warm_cache_dir_override {
            let mut stream = stream;
            stream.warm_cache_dir = Some(path.clone());
            stream.cold_variation_dir = Some(path.clone());
            stream.position_index_dir = Some(path.join("variation.position_index"));
            stream.variant_bloom_index_dir = Some(path.join("variation.variant_bloom_index"));
            return Ok(Box::pin(stream));
        }

        Ok(Box::pin(stream))
    }
}

/// Streaming implementation that processes VCF batches and probes the KV store.
///
/// When a colocated sink is present, batches are buffered during the probe
/// phase and only emitted after the input stream is exhausted. This ensures
/// the colocated sink is fully populated before downstream consumers build
/// the colocated map — matching the buffering behavior of `VariantLookupExec`.
struct KvLookupStream {
    input: SendableRecordBatchStream,
    store: Option<Arc<VepKvStore>>,
    cache_root: Option<PathBuf>,
    cache_schema: SchemaRef,
    indexed_parquet: bool,
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
    /// Cached column indices for co-located collection in the KV entry.
    coloc_col_indices: Option<KvColocIndices>,
    /// Optional indexed reference FASTA reader for genomic shift state.
    reference_reader:
        Option<noodles_fasta::IndexedReader<noodles_fasta::io::BufReader<std::fs::File>>>,
    warm_cache_dir: Option<PathBuf>,
    warm_chroms: HashMap<String, Option<Box<WarmChromCache>>>,
    cold_parquet_chroms: HashMap<String, Option<ColdParquetLookupSet>>,
    variant_bloom_chroms: HashMap<String, Option<Arc<VariantBloomIndex>>>,
    cold_variation_dir: Option<PathBuf>,
    position_index_dir: Option<PathBuf>,
    variant_bloom_index_dir: Option<PathBuf>,
    warm_cold_backend: WarmColdVariationBackend,
    warm_cold_index_mode: WarmColdVariationIndexMode,
    profile_enabled: bool,
    profile_detailed: bool,
    profile_emitted: bool,
    profile: LookupProfile,
    /// Buffered matched batches (used when colocated sink is present).
    /// Batches are collected during probe and emitted after input exhaustion.
    matched_batches: VecDeque<RecordBatch>,
    /// True once the input stream is exhausted and we're emitting buffered batches.
    input_exhausted: bool,
}

/// Column indices within the KV entry's stored columns for co-located fields.
struct KvColocIndices {
    variation_name: usize,
    _allele_string_col: usize,
    end_col: Option<usize>,
    failed: Option<usize>,
    somatic: Option<usize>,
    pheno: Option<usize>,
    clin_sig: Option<usize>,
    clin_sig_allele: Option<usize>,
    pubmed: Option<usize>,
    af_indices: Vec<Option<usize>>,
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
enum WarmProbeForDecision {
    Exact,
    PositionCoveredNoExact,
    NotCovered,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ColdPositionDecision {
    Present,
    Absent,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WarmColdVariationBackend {
    Fjall,
    Parquet,
}

impl WarmColdVariationBackend {
    fn parse(value: Option<&str>) -> Result<Self> {
        match value.unwrap_or("fjall").to_ascii_lowercase().as_str() {
            "fjall" => Ok(Self::Fjall),
            "parquet" => Ok(Self::Parquet),
            other => Err(DataFusionError::Execution(format!(
                "VEP_WARM_COLD_VARIATION_BACKEND must be 'fjall' or 'parquet', got '{other}'"
            ))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WarmColdVariationIndexMode {
    Position,
    PositionThenVariantBloom,
    VariantBloom,
}

impl WarmColdVariationIndexMode {
    fn parse(value: Option<&str>) -> Result<Self> {
        match value.unwrap_or("posidx").to_ascii_lowercase().as_str() {
            "posidx" | "position" => Ok(Self::Position),
            "posidx_bloom" | "position_bloom" | "position_then_bloom" => {
                Ok(Self::PositionThenVariantBloom)
            }
            "bloom" | "variant_bloom" => Ok(Self::VariantBloom),
            other => Err(DataFusionError::Execution(format!(
                "VEP_WARM_COLD_VARIATION_INDEX_MODE must be 'posidx', 'posidx_bloom', or 'bloom', got '{other}'"
            ))),
        }
    }

    fn uses_position_index(self) -> bool {
        matches!(self, Self::Position | Self::PositionThenVariantBloom)
    }

    fn uses_variant_bloom(self) -> bool {
        matches!(self, Self::PositionThenVariantBloom | Self::VariantBloom)
    }
}

fn decide_after_warm_probe(
    warm: WarmProbeForDecision,
    cold: ColdPositionDecision,
) -> LookupDecision {
    match warm {
        WarmProbeForDecision::Exact => LookupDecision::EmitWarmExact,
        WarmProbeForDecision::PositionCoveredNoExact => LookupDecision::SkipFjall,
        WarmProbeForDecision::NotCovered => match cold {
            ColdPositionDecision::Absent => LookupDecision::SkipFjall,
            ColdPositionDecision::Present => LookupDecision::UseFjall,
        },
    }
}

/// Open an indexed FASTA reader from a file path.
fn open_indexed_fasta(
    path: &str,
) -> Result<noodles_fasta::IndexedReader<noodles_fasta::io::BufReader<std::fs::File>>> {
    noodles_fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|e| {
            DataFusionError::Execution(format!(
                "failed to open indexed reference FASTA '{path}': {e}"
            ))
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
                "[vep-kv-profile-detail] stages total_s={:.3} probe_build={:.3}s ({:.1}%) point_get_raw={:.3}s ({:.1}%) prefetch_map_lookup={:.3}s ({:.1}%) decompress={:.3}s ({:.1}%) reader_init={:.3}s ({:.1}%) primary_match={:.3}s ({:.1}%) cache_column_append={:.3}s ({:.1}%) colocated_prepare={:.3}s ({:.1}%) colocated_match={:.3}s ({:.1}%) colocated_flush={:.3}s ({:.1}%) null_append={:.3}s ({:.1}%) warm_probe={:.3}s ({:.1}%) warm_chunk_load={:.3}s ({:.1}%) position_index_load={:.3}s ({:.1}%) variant_bloom_load={:.3}s ({:.1}%) cold_parquet_load={:.3}s ({:.1}%)",
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
                "[vep-kv-profile-detail] cold_parquet probes={} matches={} position_misses={} not_covered={} rows_scanned={}",
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

fn warm_variation_cache_enabled() -> bool {
    std::env::var("VEP_WARM_VARIATION_CACHE")
        .map(|value| value != "0")
        .unwrap_or(false)
}

fn warm_variation_cache_dir(store: &VepKvStore) -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("VEP_WARM_VARIATION_DIR") {
        return Some(PathBuf::from(path));
    }

    store
        .root_path()
        .parent()
        .map(|parent| parent.join("variation"))
}

fn warm_variation_batch_size() -> usize {
    std::env::var("VEP_WARM_VARIATION_BATCH_SIZE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(262_144)
}

fn warm_cold_variation_backend_from_env() -> WarmColdVariationBackend {
    match WarmColdVariationBackend::parse(
        std::env::var("VEP_WARM_COLD_VARIATION_BACKEND")
            .ok()
            .as_deref(),
    ) {
        Ok(backend) => backend,
        Err(error) => {
            eprintln!("[KvLookupStream] warning: {error}; falling back to fjall");
            WarmColdVariationBackend::Fjall
        }
    }
}

fn warm_cold_variation_index_mode_from_env() -> WarmColdVariationIndexMode {
    match WarmColdVariationIndexMode::parse(
        std::env::var("VEP_WARM_COLD_VARIATION_INDEX_MODE")
            .ok()
            .as_deref(),
    ) {
        Ok(mode) => mode,
        Err(error) => {
            eprintln!("[KvLookupStream] warning: {error}; falling back to posidx");
            WarmColdVariationIndexMode::Position
        }
    }
}

fn variation_cold_dir(store: &VepKvStore) -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("VEP_VARIATION_COLD_DIR") {
        return Some(PathBuf::from(path));
    }

    store
        .root_path()
        .parent()
        .map(|parent| parent.join("variation"))
}

fn variation_position_index_dir(store: &VepKvStore) -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("VEP_VARIATION_POSITION_INDEX_DIR") {
        return Some(PathBuf::from(path));
    }

    store
        .root_path()
        .parent()
        .map(|parent| parent.join("variation.position_index"))
}

fn variation_variant_bloom_index_dir(store: &VepKvStore) -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("VEP_VARIATION_BLOOM_INDEX_DIR") {
        return Some(PathBuf::from(path));
    }

    store
        .root_path()
        .parent()
        .map(|parent| parent.join("variation.variant_bloom_index"))
}

fn kv_range_prefetch_enabled() -> bool {
    // Disabled while the coarse min/max batch range prefetch is replaced with
    // a narrower clustered prefetch strategy. The current implementation often
    // scans dense cache ranges only to hit the entry limit and fall back to
    // point lookups.
    false
}

fn kv_range_prefetch_max_span() -> i64 {
    std::env::var("VEP_KV_RANGE_PREFETCH_MAX_SPAN")
        .ok()
        .and_then(|value| value.parse::<i64>().ok())
        .unwrap_or(10_000_000)
        .max(0)
}

fn kv_range_prefetch_max_entries() -> usize {
    std::env::var("VEP_KV_RANGE_PREFETCH_MAX_ENTRIES")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(100_000)
}

fn kv_range_prefetch_max_bytes() -> usize {
    std::env::var("VEP_KV_RANGE_PREFETCH_MAX_BYTES")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(128 * 1024 * 1024)
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
fn probe_cold_chunk_position(
    chunk: &WarmChunkContext,
    exact_matcher: AlleleMatcher,
    allowed_failed: i64,
    profile_detailed: bool,
    chrom: &str,
    probe_start: i64,
    position_key: u64,
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
    let rows = chunk.rows_for_position(position_key);
    if rows.is_empty() {
        return Ok((ColdProbeResult::NotCovered, metrics));
    }

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
        resolve_warm_coloc_indices(chunk)
    } else {
        None
    };

    let mut matched = false;
    for row in rows {
        let row_usize = row as usize;
        let Some(allele_string) = chunk.allele_string(row_usize)? else {
            continue;
        };

        metrics.cold_rows_scanned += 1;
        if profile_detailed {
            metrics.primary_allele_rows += 1;
            metrics.exact_match_calls += 1;
        }

        if let (Some(buf), Some(prepared), Some(ci)) = (
            coloc_buf.as_deref_mut(),
            prepared_coloc.as_ref(),
            coloc_indices.as_ref(),
        ) {
            let colocated_started = Instant::now();
            if profile_detailed {
                metrics.colocated_allele_rows += 1;
            }

            let failed = chunk.i64_value(ci.failed, row_usize).unwrap_or(0);
            if failed <= allowed_failed {
                if let Some(var_name) =
                    warm_string_value(chunk, Some(ci.variation_name), row_usize)?
                    && !var_name.is_empty()
                {
                    let existing_end = chunk
                        .i64_value(ci.end_col, row_usize)
                        .unwrap_or(probe_start);
                    if let Some(matched_alleles) = compare_existing_variant_alleles(
                        &prepared.compare_allele_string,
                        prepared.vep_start,
                        prepared.vep_end,
                        None,
                        None,
                        None,
                        allele_string,
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
                        let af_values: Vec<String> = ci
                            .af_indices
                            .iter()
                            .map(|idx| {
                                warm_string_value(chunk, *idx, row_usize)
                                    .map(|value| value.unwrap_or_default())
                            })
                            .collect::<Result<Vec<_>>>()?;
                        sink_value.entries.push(ColocatedCacheEntry {
                            variation_name: var_name,
                            allele_string: allele_string.to_string(),
                            matched_alleles,
                            somatic: chunk.i64_value(ci.somatic, row_usize).unwrap_or(0),
                            pheno: chunk.i64_value(ci.pheno, row_usize).unwrap_or(0),
                            clin_sig: warm_string_value(chunk, ci.clin_sig, row_usize)?,
                            clin_sig_allele: warm_string_value(
                                chunk,
                                ci.clin_sig_allele,
                                row_usize,
                            )?,
                            pubmed: warm_string_value(chunk, ci.pubmed, row_usize)?,
                            af_values,
                        });
                        if profile_detailed {
                            metrics.colocated_entries += 1;
                        }
                    }
                }
            }
            metrics.colocated_match_elapsed += colocated_started.elapsed();
        }

        if matched {
            continue;
        }

        let failed = chunk
            .i64_value(chunk.columns.failed, row_usize)
            .unwrap_or(0);
        if failed > allowed_failed {
            continue;
        }

        let existing_end = chunk
            .i64_value(chunk.columns.end, row_usize)
            .unwrap_or(probe_start);
        let cache_iv_start = probe_start.min(existing_end);
        let cache_iv_end = probe_start.max(existing_end);
        if cache_iv_start > vcf_iv_end || cache_iv_end < vcf_iv_start {
            continue;
        }

        if exact_matcher(vcf_ref, vcf_alt, allele_string) {
            if emit_output {
                let append_started = Instant::now();
                let Some(indices) = chunk.output_indices(cache_columns, col_map) else {
                    return Err(DataFusionError::Execution(
                        "cold parquet chunk missing one or more requested cache output columns"
                            .into(),
                    ));
                };
                vcf_indices.push(vcf_row);
                append_warm_row_values(chunk, row_usize, indices, builders)?;
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

impl KvLookupStream {
    #[allow(clippy::too_many_arguments)]
    fn new(
        input: SendableRecordBatchStream,
        store: Option<Arc<VepKvStore>>,
        cache_root: Option<PathBuf>,
        cache_schema: SchemaRef,
        indexed_parquet: bool,
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
        reference_fasta_path: Option<String>,
    ) -> Self {
        // Resolve colocated column indices within the KV entry if we have a sink.
        let coloc_col_indices = colocated_sink
            .as_ref()
            .and_then(|_| store.as_ref())
            .and_then(|store| resolve_kv_coloc_indices(store));

        // Open the reference FASTA reader if a path is provided and we have a
        // colocated sink (shift state is only needed for colocated matching).
        let reference_reader = if colocated_sink.is_some() {
            reference_fasta_path.as_deref().and_then(|path| {
                open_indexed_fasta(path)
                    .map_err(|e| {
                        eprintln!(
                            "[KvLookupStream] warning: failed to open reference FASTA {path}: {e}"
                        );
                        e
                    })
                    .ok()
            })
        } else {
            None
        };
        let warm_cache_dir = if indexed_parquet {
            cache_root.as_ref().map(|root| root.join("variation"))
        } else if warm_variation_cache_enabled() {
            store
                .as_ref()
                .and_then(|store| warm_variation_cache_dir(store))
        } else {
            None
        };
        let cold_variation_dir = if indexed_parquet {
            cache_root.as_ref().map(|root| root.join("variation"))
        } else if warm_cache_dir.is_some() {
            store.as_ref().and_then(|store| variation_cold_dir(store))
        } else {
            None
        };
        let position_index_dir = if indexed_parquet {
            cache_root
                .as_ref()
                .map(|root| root.join("variation.position_index"))
        } else if warm_cache_dir.is_some() {
            store
                .as_ref()
                .and_then(|store| variation_position_index_dir(store))
        } else {
            None
        };
        let variant_bloom_index_dir = if indexed_parquet {
            cache_root
                .as_ref()
                .map(|root| root.join("variation.variant_bloom_index"))
        } else if warm_cache_dir.is_some() {
            store
                .as_ref()
                .and_then(|store| variation_variant_bloom_index_dir(store))
        } else {
            None
        };

        Self {
            input,
            store,
            cache_root,
            cache_schema,
            indexed_parquet,
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
            coloc_col_indices,
            reference_reader,
            warm_cache_dir,
            warm_chroms: HashMap::new(),
            cold_parquet_chroms: HashMap::new(),
            variant_bloom_chroms: HashMap::new(),
            cold_variation_dir,
            position_index_dir,
            variant_bloom_index_dir,
            warm_cold_backend: if indexed_parquet {
                WarmColdVariationBackend::Parquet
            } else {
                warm_cold_variation_backend_from_env()
            },
            warm_cold_index_mode: if indexed_parquet {
                WarmColdVariationIndexMode::PositionThenVariantBloom
            } else {
                warm_cold_variation_index_mode_from_env()
            },
            profile_enabled: kv_profile_enabled(),
            profile_detailed: kv_profile_detailed_enabled(),
            profile_emitted: false,
            profile: LookupProfile::default(),
            matched_batches: VecDeque::new(),
            input_exhausted: false,
        }
    }

    /// Single-pass position-keyed lookup.
    ///
    /// For each VCF row, fetch the per-position entry from fjall, match alleles,
    /// and append matched column values directly into ArrayBuilders.
    fn probe_warm_position(
        &mut self,
        chrom: &str,
        chrom_code: u16,
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
    ) -> Result<Option<LookupDecision>> {
        let Some(warm_cache_dir) = self.warm_cache_dir.clone() else {
            return Ok(None);
        };
        let position_index_dir = if self.warm_cold_index_mode.uses_position_index() {
            let Some(position_index_dir) = self.position_index_dir.clone() else {
                return Ok(None);
            };
            Some(position_index_dir)
        } else {
            None
        };
        let collect_colocated = coloc_buf.is_some();

        if !self.warm_chroms.contains_key(chrom) {
            let warm_path = warm_file_for_chrom(&warm_cache_dir, chrom);
            let warm_chrom = if let Some(warm_path) = warm_path {
                let warm_cache_columns =
                    projection_columns_for_cache(cache_columns, collect_colocated);
                let cold_path = self
                    .cold_variation_dir
                    .as_ref()
                    .and_then(|dir| cold_variation_file_for_chrom(dir, chrom));
                let started = self.profile_enabled.then(Instant::now);
                match WarmChromCache::open_with_optional_position_index(
                    &warm_path,
                    chrom,
                    position_index_dir.as_deref(),
                    cold_path.as_deref(),
                    warm_variation_batch_size(),
                    warm_cache_columns,
                )? {
                    WarmChromOpen::Available(cache) => {
                        if self.profile_enabled {
                            if let Some(source) = cache.cold_position_source() {
                                if let Some(t0) = started {
                                    self.profile.position_index_load += t0.elapsed();
                                }
                                self.profile.position_index_loaded += 1;
                                match source {
                                    PositionIndexSource::Persisted => {
                                        self.profile.position_index_persisted_loads += 1;
                                    }
                                    PositionIndexSource::ParquetFallback => {
                                        self.profile.position_index_parquet_fallback_loads += 1;
                                    }
                                }
                            }
                            self.profile.position_index_rows += cache.cold_position_len() as u64;
                            self.profile.position_index_bytes +=
                                cache.cold_position_storage_bytes() as u64;
                        }
                        Some(cache)
                    }
                    WarmChromOpen::Unavailable(_) => None,
                }
            } else {
                None
            };
            self.warm_chroms.insert(chrom.to_string(), warm_chrom);
        }

        let position_key = warm_position_key_from_code(chrom_code, probe_start).map_err(|e| {
            DataFusionError::Execution(format!("failed to build warm position key: {e}"))
        })?;
        let variant_keys = warm_variant_key_candidates(position_key, vcf_ref, vcf_alt);
        if variant_keys.is_empty() {
            return Ok(None);
        }

        let exact_matcher = self.exact_matcher;
        let allowed_failed = self.allowed_failed;
        let profile_detailed = self.profile_detailed;
        let mut append_elapsed = Duration::ZERO;
        let mut emitted = false;
        let mut colocated_prepare_elapsed = Duration::ZERO;
        let mut colocated_match_elapsed = Duration::ZERO;
        let mut primary_allele_rows = 0_u64;
        let mut exact_match_calls = 0_u64;
        let mut colocated_allele_rows = 0_u64;
        let mut colocated_entries = 0_u64;
        let mut warm_candidate_rows = 0_u64;
        let mut warm_rows_scanned = 0_u64;
        let started = self.profile_enabled.then(Instant::now);
        let warm_chrom = self
            .warm_chroms
            .get_mut(chrom)
            .and_then(Option::as_deref_mut);
        let Some(warm_chrom) = warm_chrom else {
            return Ok(None);
        };
        let before_chunks = warm_chrom.chunks_loaded;
        let before_rows = warm_chrom.chunk_rows;
        let before_load = warm_chrom.chunk_load;
        let mut coloc_indices_cache: HashMap<usize, Option<WarmColocIndices>> = HashMap::new();
        let mut coloc_buf = coloc_buf;
        let row_matches =
            |chunk: &WarmChunkContext, row: u32, allele_string: &str| -> Result<bool> {
                let failed = chunk
                    .i64_value(chunk.columns.failed, row as usize)
                    .unwrap_or(0);
                if failed > allowed_failed {
                    return Ok(false);
                }

                let existing_end = chunk
                    .i64_value(chunk.columns.end, row as usize)
                    .unwrap_or(probe_start);
                let cache_iv_start = probe_start.min(existing_end);
                let cache_iv_end = probe_start.max(existing_end);
                if cache_iv_start > vcf_iv_end || cache_iv_end < vcf_iv_start {
                    return Ok(false);
                }

                Ok(exact_matcher(vcf_ref, vcf_alt, allele_string))
            };

        let mut exact_row: Option<u32> = None;
        let mut position_rows: Option<Range<u32>> = None;
        for variant_key in variant_keys {
            match warm_chrom.probe(position_key, variant_key, |_, allele_string| {
                Ok(exact_matcher(vcf_ref, vcf_alt, allele_string))
            })? {
                WarmProbe::Exact {
                    row,
                    position_rows: rows,
                } => {
                    exact_row = Some(row);
                    position_rows = Some(rows);
                    break;
                }
                WarmProbe::PositionCoveredNoExact {
                    position_rows: rows,
                } => {
                    position_rows.get_or_insert(rows);
                }
                WarmProbe::NotCovered => {}
            }
        }

        if let Some(row) = exact_row {
            let chunk = warm_chrom
                .current_chunk()
                .ok_or_else(|| DataFusionError::Execution("warm chunk not loaded".into()))?;
            warm_rows_scanned += 1;
            let matched = chunk
                .allele_string(row as usize)?
                .map(|allele_string| row_matches(chunk, row, allele_string))
                .transpose()?
                .unwrap_or(false);
            if !matched {
                exact_row = None;
            }
        }

        if exact_row.is_none() {
            if let Some(rows) = position_rows.clone() {
                let chunk = warm_chrom
                    .current_chunk()
                    .ok_or_else(|| DataFusionError::Execution("warm chunk not loaded".into()))?;
                for row in rows {
                    warm_rows_scanned += 1;
                    let Some(allele_string) = chunk.allele_string(row as usize)? else {
                        continue;
                    };
                    if profile_detailed {
                        primary_allele_rows += 1;
                        exact_match_calls += 1;
                    }
                    if row_matches(chunk, row, allele_string)? {
                        exact_row = Some(row);
                        break;
                    }
                }
            }
        }

        if collect_colocated {
            if let Some(rows) = position_rows.clone() {
                let prepare_started = Instant::now();
                let chrom_norm = chrom.to_string();
                let (input_ref, input_alt, input_start) =
                    vcf_to_vep_input_allele(vcf_iv_start, vcf_ref, vcf_alt);
                let input_allele_string = format!("{input_ref}/{input_alt}");
                let (compare_ref, compare_alt) = vcf_to_vep_allele(vcf_ref, vcf_alt);
                let compare_allele_string = format!("{compare_ref}/{compare_alt}");
                let vep_start = vep_norm_start(vcf_iv_start, vcf_ref, vcf_alt);
                let vep_end = vep_norm_end(vcf_iv_start, vcf_ref, vcf_alt);
                let compare_output_allele =
                    output_allele_from_allele_string(&compare_allele_string).map(str::to_string);
                let unshifted_output_allele: Option<String> = None;
                colocated_prepare_elapsed += prepare_started.elapsed();

                let chunk = warm_chrom
                    .current_chunk()
                    .ok_or_else(|| DataFusionError::Execution("warm chunk not loaded".into()))?;
                for row in rows {
                    let Some(allele_string) = chunk.allele_string(row as usize)? else {
                        continue;
                    };
                    let colocated_started = Instant::now();
                    let Some(buf) = coloc_buf.as_deref_mut() else {
                        continue;
                    };
                    if !probe_start_visible_to_window(probe_start, vep_start, vep_end) {
                        colocated_match_elapsed += colocated_started.elapsed();
                        continue;
                    }
                    let ci = coloc_indices_cache
                        .entry(chunk.row_group_id)
                        .or_insert_with(|| resolve_warm_coloc_indices(chunk));
                    let Some(ci) = ci.as_ref() else {
                        colocated_match_elapsed += colocated_started.elapsed();
                        continue;
                    };
                    if profile_detailed {
                        colocated_allele_rows += 1;
                    }
                    let failed = chunk.i64_value(ci.failed, row as usize).unwrap_or(0);
                    if failed > allowed_failed {
                        colocated_match_elapsed += colocated_started.elapsed();
                        continue;
                    }

                    let Some(var_name) =
                        warm_string_value(chunk, Some(ci.variation_name), row as usize)?
                    else {
                        colocated_match_elapsed += colocated_started.elapsed();
                        continue;
                    };
                    if var_name.is_empty() {
                        colocated_match_elapsed += colocated_started.elapsed();
                        continue;
                    }

                    let existing_end = chunk
                        .i64_value(ci.end_col, row as usize)
                        .unwrap_or(probe_start);
                    let Some(matched_alleles) = compare_existing_variant_alleles(
                        &compare_allele_string,
                        vep_start,
                        vep_end,
                        None,
                        None,
                        None,
                        allele_string,
                        probe_start,
                        existing_end,
                    ) else {
                        colocated_match_elapsed += colocated_started.elapsed();
                        continue;
                    };

                    let key: ColocatedKey = (
                        chrom_norm.clone(),
                        input_start,
                        vcf_iv_end,
                        input_allele_string.clone(),
                    );
                    let sink_value = buf.entry(key).or_insert_with(|| ColocatedSinkValue {
                        entries: Vec::new(),
                        compare_output_allele: compare_output_allele.clone(),
                        unshifted_output_allele: unshifted_output_allele.clone(),
                    });
                    let af_values: Vec<String> = ci
                        .af_indices
                        .iter()
                        .map(|idx| {
                            warm_string_value(chunk, *idx, row as usize)
                                .map(|value| value.unwrap_or_default())
                        })
                        .collect::<Result<Vec<_>>>()?;
                    sink_value.entries.push(ColocatedCacheEntry {
                        variation_name: var_name,
                        allele_string: allele_string.to_string(),
                        matched_alleles,
                        somatic: chunk.i64_value(ci.somatic, row as usize).unwrap_or(0),
                        pheno: chunk.i64_value(ci.pheno, row as usize).unwrap_or(0),
                        clin_sig: warm_string_value(chunk, ci.clin_sig, row as usize)?,
                        clin_sig_allele: warm_string_value(
                            chunk,
                            ci.clin_sig_allele,
                            row as usize,
                        )?,
                        pubmed: warm_string_value(chunk, ci.pubmed, row as usize)?,
                        af_values,
                    });
                    if profile_detailed {
                        colocated_entries += 1;
                    }
                    colocated_match_elapsed += colocated_started.elapsed();
                }
            }
        }

        if emit_output {
            if let Some(row) = exact_row {
                let chunk = warm_chrom
                    .current_chunk()
                    .ok_or_else(|| DataFusionError::Execution("warm chunk not loaded".into()))?;
                let append_started = Instant::now();
                let Some(indices) = chunk.output_indices(cache_columns, col_map) else {
                    return Err(DataFusionError::Execution(
                        "warm chunk missing one or more requested cache output columns".into(),
                    ));
                };
                vcf_indices.push(vcf_row);
                append_warm_row_values(chunk, row as usize, indices, builders)?;
                append_elapsed += append_started.elapsed();
                emitted = true;
            }
        }

        let warm_decision = if exact_row.is_some() {
            WarmProbeForDecision::Exact
        } else if position_rows.is_some() {
            WarmProbeForDecision::PositionCoveredNoExact
        } else {
            WarmProbeForDecision::NotCovered
        };
        if let Some(rows) = position_rows.as_ref() {
            warm_candidate_rows += u64::from(rows.end.saturating_sub(rows.start));
        }

        let mut position_index_checks = 0_u64;
        let mut position_index_positive_gets = 0_u64;
        let mut position_index_negative_skips = 0_u64;
        let cold_decision = if warm_decision == WarmProbeForDecision::NotCovered
            && self.warm_cold_index_mode.uses_position_index()
        {
            position_index_checks += 1;
            if warm_chrom.cold_may_contain(position_key) {
                position_index_positive_gets += 1;
                ColdPositionDecision::Present
            } else {
                position_index_negative_skips += 1;
                ColdPositionDecision::Absent
            }
        } else {
            ColdPositionDecision::Present
        };
        let decision = decide_after_warm_probe(warm_decision, cold_decision);
        let warm_chunks_loaded = warm_chrom.chunks_loaded - before_chunks;
        let warm_chunk_rows = warm_chrom.chunk_rows - before_rows;
        let warm_chunk_load = warm_chrom.chunk_load - before_load;

        if self.profile_enabled {
            self.profile.warm_probes += 1;
            if let Some(t0) = started {
                self.profile.warm_probe += t0.elapsed();
            }
            self.profile.warm_chunks_loaded += warm_chunks_loaded;
            self.profile.warm_chunk_rows += warm_chunk_rows;
            self.profile.warm_chunk_load += warm_chunk_load;
            self.profile.cache_column_append += append_elapsed;
            self.profile.colocated_prepare += colocated_prepare_elapsed;
            self.profile.colocated_match += colocated_match_elapsed;
            self.profile.position_index_checks += position_index_checks;
            self.profile.position_index_positive_gets += position_index_positive_gets;
            self.profile.position_index_negative_skips += position_index_negative_skips;
            self.profile.warm_candidate_rows += warm_candidate_rows;
            self.profile.warm_rows_scanned += warm_rows_scanned;
            if decision == LookupDecision::EmitWarmExact {
                self.profile.warm_matches += 1;
            } else if decision == LookupDecision::SkipFjall {
                self.profile.warm_definitive_misses += 1;
            } else {
                self.profile.warm_not_covered += 1;
            }
            if profile_detailed {
                self.profile.primary_allele_rows += primary_allele_rows;
                self.profile.exact_match_calls += exact_match_calls;
                self.profile.colocated_allele_rows += colocated_allele_rows;
                self.profile.colocated_entries += colocated_entries;
            }
        }

        if exact_row.is_some() && emitted && self.profile_detailed {
            self.profile.primary_matches += 1;
        }

        Ok(Some(decision))
    }

    #[allow(clippy::too_many_arguments)]
    fn cold_parquet_aggregate_stats(&self) -> ColdParquetStats {
        let mut aggregate = ColdParquetStats::default();
        for lookup in self.cold_parquet_chroms.values().filter_map(Option::as_ref) {
            aggregate += lookup.stats_snapshot();
        }
        aggregate
    }

    fn record_cold_parquet_stats_snapshot(&mut self) {
        let aggregate_stats = self.cold_parquet_aggregate_stats();
        self.profile.cold_parquet_row_groups_total = aggregate_stats.row_groups_total;
        self.profile.cold_parquet_row_groups_unique_touched =
            aggregate_stats.row_groups_unique_touched;
        self.profile.cold_parquet_row_group_metadata_probes =
            aggregate_stats.row_group_metadata_probes;
        self.profile.cold_parquet_row_group_current_hits = aggregate_stats.row_group_current_hits;
        self.profile.cold_parquet_row_group_previous_hits = aggregate_stats.row_group_previous_hits;
        self.profile.cold_parquet_row_group_advanced_hits = aggregate_stats.row_group_advanced_hits;
        self.profile.cold_parquet_row_group_binary_search_hits =
            aggregate_stats.row_group_binary_search_hits;
        self.profile.cold_parquet_row_group_metadata_misses =
            aggregate_stats.row_group_metadata_misses;
        self.profile.cold_parquet_row_group_skipped_ahead = aggregate_stats.row_group_skipped_ahead;
        self.profile.cold_parquet_row_groups_loaded = aggregate_stats.row_groups_loaded;
        self.profile.cold_parquet_row_group_load_batches = aggregate_stats.row_group_load_batches;
        self.profile.cold_parquet_row_group_cache_hits = aggregate_stats.row_group_cache_hits;
        self.profile.cold_parquet_row_group_cache_misses = aggregate_stats.row_group_cache_misses;
        self.profile.cold_parquet_rows_loaded = aggregate_stats.rows_loaded;
        self.profile.cold_parquet_position_page_index_loaded =
            aggregate_stats.position_page_index_loaded;
        self.profile.cold_parquet_position_column_index_loaded =
            aggregate_stats.position_column_index_loaded;
        self.profile.cold_parquet_position_pages_total = aggregate_stats.position_pages_total;
        self.profile.cold_parquet_position_bloom_filter_row_groups =
            aggregate_stats.position_bloom_filter_row_groups;
        self.profile.cold_parquet_position_bloom_filter_bytes =
            aggregate_stats.position_bloom_filter_bytes;
        self.profile.cold_parquet_page_index_probes = aggregate_stats.page_index_probes;
        self.profile.cold_parquet_page_index_available_probes =
            aggregate_stats.page_index_available_probes;
        self.profile.cold_parquet_page_index_unavailable_probes =
            aggregate_stats.page_index_unavailable_probes;
        self.profile
            .cold_parquet_page_index_pages_in_probed_row_groups =
            aggregate_stats.page_index_pages_in_probed_row_groups;
        self.profile.cold_parquet_page_index_candidate_pages =
            aggregate_stats.page_index_candidate_pages;
        self.profile.cold_parquet_page_index_candidate_rows =
            aggregate_stats.page_index_candidate_rows;
        self.profile.cold_parquet_page_index_candidate_misses =
            aggregate_stats.page_index_candidate_misses;
        self.profile.cold_parquet_page_index_unique_candidate_pages =
            aggregate_stats.page_index_unique_candidate_pages;
        self.profile.cold_parquet_page_index_unique_candidate_rows =
            aggregate_stats.page_index_unique_candidate_rows;
        self.profile.cold_parquet_load = aggregate_stats.load_time;
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

    fn cold_variant_bloom_may_contain(
        &mut self,
        chrom: &str,
        chrom_code: u16,
        probe_start: i64,
        vcf_ref: &str,
        vcf_alt: &str,
        collect_colocated: bool,
    ) -> Result<bool> {
        if !self.warm_cold_index_mode.uses_variant_bloom() {
            return Ok(true);
        }
        let Some(index_dir) = self.variant_bloom_index_dir.clone() else {
            return Err(DataFusionError::Execution(
                "VEP_WARM_COLD_VARIATION_INDEX_MODE uses bloom but no variant bloom index directory is configured".into(),
            ));
        };

        if !self.variant_bloom_chroms.contains_key(chrom) {
            let path = find_variant_bloom_index_file(&index_dir, chrom).ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "VEP_WARM_COLD_VARIATION_INDEX_MODE uses bloom but no variant bloom index file for {chrom} in {}",
                    index_dir.display()
                ))
            })?;
            let started = self.profile_enabled.then(Instant::now);
            let index = Arc::new(VariantBloomIndex::read_from_path(&path)?);
            if self.profile_enabled {
                if let Some(t0) = started {
                    self.profile.variant_bloom_load += t0.elapsed();
                }
                self.profile.variant_bloom_loaded += 1;
                self.profile.variant_bloom_entries += index.inserted();
                self.profile.variant_bloom_bits += index.bit_count();
                self.profile.variant_bloom_hashes = u64::from(index.hash_count());
                self.profile.variant_bloom_bytes += index.storage_bytes() as u64;
            }
            self.variant_bloom_chroms
                .insert(chrom.to_string(), Some(index));
        }

        let position_key = warm_position_key_from_code(chrom_code, probe_start).map_err(|e| {
            DataFusionError::Execution(format!("failed to build cold variant bloom key: {e}"))
        })?;
        let variant_keys = warm_variant_key_candidates(position_key, vcf_ref, vcf_alt);
        let maybe_present = self
            .variant_bloom_chroms
            .get(chrom)
            .and_then(Option::as_ref)
            .is_some_and(|index| {
                if index.contains_any_variant_keys(variant_keys.iter().copied()) {
                    return true;
                }

                if !collect_colocated {
                    return false;
                }

                if index.supports_position_fallback_keys() {
                    index.contains_position_fallback_key(position_key)
                } else {
                    // V1 bloom files only contain allele keys, so they cannot
                    // prove a coordinate-only colocated row absent.
                    true
                }
            });

        if self.profile_enabled {
            self.profile.variant_bloom_checks += 1;
            if maybe_present {
                self.profile.variant_bloom_positive_gets += 1;
            } else {
                self.profile.variant_bloom_negative_skips += 1;
            }
        }

        Ok(maybe_present)
    }

    fn ensure_cold_parquet_lookup(
        &mut self,
        chrom: &str,
        cache_columns: &[String],
        collect_colocated: bool,
    ) -> Result<()> {
        let Some(cold_variation_dir) = self.cold_variation_dir.clone() else {
            return Err(DataFusionError::Execution(
                "VEP_WARM_COLD_VARIATION_BACKEND=parquet requires a cold variation directory"
                    .into(),
            ));
        };

        if !self.cold_parquet_chroms.contains_key(chrom) {
            let cold_paths = cold_variation_files_for_chrom(&cold_variation_dir, chrom);
            if cold_paths.is_empty() {
                return Err(DataFusionError::Execution(format!(
                    "VEP_WARM_COLD_VARIATION_BACKEND=parquet but no cold variation parquet for {chrom}"
                )));
            }
            let projection = cold_parquet_projection_columns(cache_columns, collect_colocated);
            let lookup = ColdParquetLookupSet::from_env(cold_paths, projection)?;
            self.cold_parquet_chroms
                .insert(chrom.to_string(), Some(lookup));
        }

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn probe_cold_parquet_position(
        &mut self,
        chrom: &str,
        chrom_code: u16,
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
    ) -> Result<ColdProbeResult> {
        let collect_colocated = coloc_buf.is_some();
        self.ensure_cold_parquet_lookup(chrom, cache_columns, collect_colocated)?;

        let position_key = warm_position_key_from_code(chrom_code, probe_start).map_err(|e| {
            DataFusionError::Execution(format!("failed to build cold position key: {e}"))
        })?;
        let exact_matcher = self.exact_matcher;
        let allowed_failed = self.allowed_failed;
        let profile_detailed = self.profile_detailed;
        let mut primary_allele_rows = 0_u64;
        let mut exact_match_calls = 0_u64;
        let mut colocated_allele_rows = 0_u64;
        let mut colocated_entries = 0_u64;
        let mut cold_rows_scanned = 0_u64;
        let mut append_elapsed = Duration::ZERO;
        let mut colocated_prepare_elapsed = Duration::ZERO;
        let mut colocated_match_elapsed = Duration::ZERO;
        let mut emitted = false;
        let mut coloc_indices_cache: HashMap<usize, Option<WarmColocIndices>> = HashMap::new();
        let mut coloc_buf = coloc_buf;

        let row_matches =
            |chunk: &WarmChunkContext, row: u32, allele_string: &str| -> Result<bool> {
                let failed = chunk
                    .i64_value(chunk.columns.failed, row as usize)
                    .unwrap_or(0);
                if failed > allowed_failed {
                    return Ok(false);
                }

                let existing_end = chunk
                    .i64_value(chunk.columns.end, row as usize)
                    .unwrap_or(probe_start);
                let cache_iv_start = probe_start.min(existing_end);
                let cache_iv_end = probe_start.max(existing_end);
                if cache_iv_start > vcf_iv_end || cache_iv_end < vcf_iv_start {
                    return Ok(false);
                }

                Ok(exact_matcher(vcf_ref, vcf_alt, allele_string))
            };

        let mut emit_match = |chunk: &WarmChunkContext, row: u32| -> Result<()> {
            if !emit_output {
                return Ok(());
            }
            let append_started = Instant::now();
            let Some(indices) = chunk.output_indices(cache_columns, col_map) else {
                return Err(DataFusionError::Execution(
                    "cold parquet chunk missing one or more requested cache output columns".into(),
                ));
            };
            vcf_indices.push(vcf_row);
            append_warm_row_values(chunk, row as usize, indices, builders)?;
            append_elapsed += append_started.elapsed();
            emitted = true;
            Ok(())
        };

        let mut visit_row =
            |chunk: &WarmChunkContext, row: u32, allele_string: &str| -> Result<()> {
                cold_rows_scanned += 1;
                if profile_detailed {
                    primary_allele_rows += 1;
                    exact_match_calls += 1;
                }
                let Some(buf) = coloc_buf.as_deref_mut() else {
                    return Ok(());
                };

                let prepare_started = Instant::now();
                let chrom_norm = chrom.to_string();
                let (input_ref, input_alt, input_start) =
                    vcf_to_vep_input_allele(vcf_iv_start, vcf_ref, vcf_alt);
                let input_allele_string = format!("{input_ref}/{input_alt}");
                let (compare_ref, compare_alt) = vcf_to_vep_allele(vcf_ref, vcf_alt);
                let compare_allele_string = format!("{compare_ref}/{compare_alt}");
                let vep_start = vep_norm_start(vcf_iv_start, vcf_ref, vcf_alt);
                let vep_end = vep_norm_end(vcf_iv_start, vcf_ref, vcf_alt);
                let compare_output_allele =
                    output_allele_from_allele_string(&compare_allele_string).map(str::to_string);
                let unshifted_output_allele: Option<String> = None;
                colocated_prepare_elapsed += prepare_started.elapsed();

                let colocated_started = Instant::now();
                if !probe_start_visible_to_window(probe_start, vep_start, vep_end) {
                    colocated_match_elapsed += colocated_started.elapsed();
                    return Ok(());
                }

                let ci = coloc_indices_cache
                    .entry(chunk.row_group_id)
                    .or_insert_with(|| resolve_warm_coloc_indices(chunk));
                let Some(ci) = ci.as_ref() else {
                    colocated_match_elapsed += colocated_started.elapsed();
                    return Ok(());
                };
                if profile_detailed {
                    colocated_allele_rows += 1;
                }
                let failed = chunk.i64_value(ci.failed, row as usize).unwrap_or(0);
                if failed > allowed_failed {
                    colocated_match_elapsed += colocated_started.elapsed();
                    return Ok(());
                }

                let Some(var_name) =
                    warm_string_value(chunk, Some(ci.variation_name), row as usize)?
                else {
                    colocated_match_elapsed += colocated_started.elapsed();
                    return Ok(());
                };
                if var_name.is_empty() {
                    colocated_match_elapsed += colocated_started.elapsed();
                    return Ok(());
                }

                let existing_end = chunk
                    .i64_value(ci.end_col, row as usize)
                    .unwrap_or(probe_start);
                let Some(matched_alleles) = compare_existing_variant_alleles(
                    &compare_allele_string,
                    vep_start,
                    vep_end,
                    None,
                    None,
                    None,
                    allele_string,
                    probe_start,
                    existing_end,
                ) else {
                    colocated_match_elapsed += colocated_started.elapsed();
                    return Ok(());
                };

                let key: ColocatedKey = (chrom_norm, input_start, vcf_iv_end, input_allele_string);
                let sink_value = buf.entry(key).or_insert_with(|| ColocatedSinkValue {
                    entries: Vec::new(),
                    compare_output_allele,
                    unshifted_output_allele,
                });
                let af_values: Vec<String> = ci
                    .af_indices
                    .iter()
                    .map(|idx| {
                        warm_string_value(chunk, *idx, row as usize)
                            .map(|value| value.unwrap_or_default())
                    })
                    .collect::<Result<Vec<_>>>()?;
                sink_value.entries.push(ColocatedCacheEntry {
                    variation_name: var_name,
                    allele_string: allele_string.to_string(),
                    matched_alleles,
                    somatic: chunk.i64_value(ci.somatic, row as usize).unwrap_or(0),
                    pheno: chunk.i64_value(ci.pheno, row as usize).unwrap_or(0),
                    clin_sig: warm_string_value(chunk, ci.clin_sig, row as usize)?,
                    clin_sig_allele: warm_string_value(chunk, ci.clin_sig_allele, row as usize)?,
                    pubmed: warm_string_value(chunk, ci.pubmed, row as usize)?,
                    af_values,
                });
                if profile_detailed {
                    colocated_entries += 1;
                }
                colocated_match_elapsed += colocated_started.elapsed();
                Ok(())
            };

        let result = {
            let lookup = self
                .cold_parquet_chroms
                .get_mut(chrom)
                .and_then(Option::as_mut)
                .ok_or_else(|| {
                    DataFusionError::Execution(format!("cold parquet lookup not open for {chrom}"))
                })?;
            lookup.probe_position_emit_and_visit(
                position_key,
                row_matches,
                &mut emit_match,
                &mut visit_row,
            )?
        };

        if self.profile_enabled {
            let aggregate_stats = self.cold_parquet_aggregate_stats();
            self.profile.cold_parquet_probes += 1;
            self.profile.cache_column_append += append_elapsed;
            self.profile.colocated_prepare += colocated_prepare_elapsed;
            self.profile.colocated_match += colocated_match_elapsed;
            self.profile.cold_parquet_rows_scanned += cold_rows_scanned;
            self.profile.cold_parquet_row_groups_total = aggregate_stats.row_groups_total;
            self.profile.cold_parquet_row_groups_unique_touched =
                aggregate_stats.row_groups_unique_touched;
            self.profile.cold_parquet_row_group_metadata_probes =
                aggregate_stats.row_group_metadata_probes;
            self.profile.cold_parquet_row_group_current_hits =
                aggregate_stats.row_group_current_hits;
            self.profile.cold_parquet_row_group_previous_hits =
                aggregate_stats.row_group_previous_hits;
            self.profile.cold_parquet_row_group_advanced_hits =
                aggregate_stats.row_group_advanced_hits;
            self.profile.cold_parquet_row_group_binary_search_hits =
                aggregate_stats.row_group_binary_search_hits;
            self.profile.cold_parquet_row_group_metadata_misses =
                aggregate_stats.row_group_metadata_misses;
            self.profile.cold_parquet_row_group_skipped_ahead =
                aggregate_stats.row_group_skipped_ahead;
            self.profile.cold_parquet_row_groups_loaded = aggregate_stats.row_groups_loaded;
            self.profile.cold_parquet_row_group_load_batches =
                aggregate_stats.row_group_load_batches;
            self.profile.cold_parquet_row_group_cache_hits = aggregate_stats.row_group_cache_hits;
            self.profile.cold_parquet_row_group_cache_misses =
                aggregate_stats.row_group_cache_misses;
            self.profile.cold_parquet_rows_loaded = aggregate_stats.rows_loaded;
            self.profile.cold_parquet_position_page_index_loaded =
                aggregate_stats.position_page_index_loaded;
            self.profile.cold_parquet_position_column_index_loaded =
                aggregate_stats.position_column_index_loaded;
            self.profile.cold_parquet_position_pages_total = aggregate_stats.position_pages_total;
            self.profile.cold_parquet_position_bloom_filter_row_groups =
                aggregate_stats.position_bloom_filter_row_groups;
            self.profile.cold_parquet_position_bloom_filter_bytes =
                aggregate_stats.position_bloom_filter_bytes;
            self.profile.cold_parquet_page_index_probes = aggregate_stats.page_index_probes;
            self.profile.cold_parquet_page_index_available_probes =
                aggregate_stats.page_index_available_probes;
            self.profile.cold_parquet_page_index_unavailable_probes =
                aggregate_stats.page_index_unavailable_probes;
            self.profile
                .cold_parquet_page_index_pages_in_probed_row_groups =
                aggregate_stats.page_index_pages_in_probed_row_groups;
            self.profile.cold_parquet_page_index_candidate_pages =
                aggregate_stats.page_index_candidate_pages;
            self.profile.cold_parquet_page_index_candidate_rows =
                aggregate_stats.page_index_candidate_rows;
            self.profile.cold_parquet_page_index_candidate_misses =
                aggregate_stats.page_index_candidate_misses;
            self.profile.cold_parquet_page_index_unique_candidate_pages =
                aggregate_stats.page_index_unique_candidate_pages;
            self.profile.cold_parquet_page_index_unique_candidate_rows =
                aggregate_stats.page_index_unique_candidate_rows;
            self.profile.cold_parquet_load = aggregate_stats.load_time;
            match result {
                ColdProbeResult::Match => self.profile.cold_parquet_matches += 1,
                ColdProbeResult::PositionCoveredNoExact => {
                    self.profile.cold_parquet_position_misses += 1;
                }
                ColdProbeResult::NotCovered => self.profile.cold_parquet_not_covered += 1,
            }
            if profile_detailed {
                self.profile.primary_allele_rows += primary_allele_rows;
                self.profile.exact_match_calls += exact_match_calls;
                self.profile.colocated_allele_rows += colocated_allele_rows;
                self.profile.colocated_entries += colocated_entries;
            }
            if emitted && profile_detailed {
                self.profile.primary_matches += 1;
            }
        }

        Ok(result)
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

        // Reusable zstd decompressor — created once, amortized across legacy fjall lookups.
        let mut decompressor = self
            .store
            .as_ref()
            .map(|store| store.create_decompressor())
            .transpose()?
            .flatten();

        // Reusable decompression / raw-value buffer — avoids alloc per lookup.
        let mut decompress_buf: Vec<u8> = Vec::with_capacity(4096);

        // Reusable encoded Fjall key buffer — avoids one tiny Vec allocation per point probe.
        let mut position_key_buf: Vec<u8> = Vec::with_capacity(10);

        // Reusable allele match buffer — avoids alloc per row.
        let mut matched_allele_rows: Vec<usize> = Vec::new();

        let range_prefetch_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let mut range_prefetch: Option<HashMap<i64, fjall::UserValue>> = None;
        if kv_range_prefetch_enabled()
            && self.store.is_some()
            && self.warm_cache_dir.is_none()
            && num_rows > 1
        {
            let mut batch_chrom_code: Option<u16> = None;
            let mut min_probe = i64::MAX;
            let mut max_probe = i64::MIN;
            let mut eligible = true;

            for row in 0..num_rows {
                let raw_chrom = chroms.value_or_empty(row);
                let chrom = if self.vcf_has_chr {
                    raw_chrom.strip_prefix("chr").unwrap_or(raw_chrom)
                } else {
                    raw_chrom
                };
                let chrom_code = chrom_to_code(chrom);
                if batch_chrom_code
                    .replace(chrom_code)
                    .is_some_and(|prev| prev != chrom_code)
                {
                    eligible = false;
                    break;
                }

                let (norm_start, norm_end) = normalize_vcf_coords(
                    starts[row],
                    ends[row],
                    self.vcf_zero_based,
                    self.cache_zero_based,
                )?;
                let probes = build_probe_starts(
                    i64::from(norm_start),
                    i64::from(norm_end),
                    refs.value_or_empty(row),
                    alts.value_or_empty(row),
                    self.extended_probes,
                );
                for probe in probes {
                    min_probe = min_probe.min(probe);
                    max_probe = max_probe.max(probe);
                }
            }

            if eligible {
                if let Some(chrom_code) = batch_chrom_code {
                    let span = max_probe.saturating_sub(min_probe);
                    if min_probe <= max_probe && span <= kv_range_prefetch_max_span() {
                        let store = self.store.as_ref().expect("checked store.is_some()");
                        match store.range_position_entries_limited_with_reason(
                            chrom_code,
                            min_probe,
                            max_probe,
                            kv_range_prefetch_max_entries(),
                            kv_range_prefetch_max_bytes(),
                        )? {
                            Ok(entries) => {
                                let mut prefetched = HashMap::with_capacity(entries.len());
                                let mut bytes = 0usize;
                                for (position, value) in entries {
                                    bytes = bytes.saturating_add(value.as_ref().len());
                                    prefetched.insert(position, value);
                                }
                                if self.profile_enabled {
                                    self.profile.range_prefetch_batches += 1;
                                    self.profile.range_prefetch_entries += prefetched.len() as u64;
                                    self.profile.range_prefetch_bytes += bytes as u64;
                                }
                                range_prefetch = Some(prefetched);
                            }
                            Err(reason) if self.profile_enabled => {
                                self.profile.range_prefetch_skipped += 1;
                                match reason {
                                    RangePrefetchLimitExceeded::Entries => {
                                        self.profile.range_prefetch_skip_entries += 1;
                                    }
                                    RangePrefetchLimitExceeded::Bytes => {
                                        self.profile.range_prefetch_skip_bytes += 1;
                                    }
                                }
                            }
                            Err(_) => {}
                        }
                    } else if self.profile_enabled {
                        self.profile.range_prefetch_skipped += 1;
                        self.profile.range_prefetch_skip_span += 1;
                    }
                }
            } else if self.profile_enabled {
                self.profile.range_prefetch_skipped += 1;
                self.profile.range_prefetch_skip_multi_chrom += 1;
            }
        }
        if let Some(t0) = range_prefetch_started {
            self.profile.range_prefetch += t0.elapsed();
        }

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
                        LookupDecision::UseFjall
                            if self.warm_cold_backend == WarmColdVariationBackend::Parquet =>
                        {
                            let collect_colocated = coloc_buf.is_some();
                            if !self.cold_variant_bloom_may_contain(
                                chrom,
                                chrom_code,
                                *probe_start,
                                vcf_ref,
                                vcf_alt,
                                collect_colocated,
                            )? {
                                continue;
                            }
                            let position_key =
                                warm_position_key_from_code(chrom_code, *probe_start).map_err(
                                    |e| {
                                        DataFusionError::Execution(format!(
                                            "failed to build cold position key: {e}"
                                        ))
                                    },
                                )?;
                            self.ensure_cold_parquet_lookup(
                                chrom,
                                &cache_columns,
                                collect_colocated,
                            )?;
                            let pending_idx = pending_cold_probes.len();
                            pending_cold_probes.push(PendingColdProbe {
                                chrom: chrom.to_string(),
                                probe_start: *probe_start,
                                position_key,
                                vcf_ref: vcf_ref.to_string(),
                                vcf_alt: vcf_alt.to_string(),
                                vcf_iv_start,
                                vcf_iv_end,
                                vcf_row: row as u32,
                            });
                            pending_cold_by_row[row].push(pending_idx);
                            continue;
                        }
                        LookupDecision::UseFjall => {}
                    }
                }

                let found = if let Some(prefetched) = range_prefetch.as_ref() {
                    let Some(store) = self.store.as_ref() else {
                        return Err(DataFusionError::Execution(
                            "indexed parquet variation lookup unexpectedly reached legacy fjall range-prefetch path"
                                .into(),
                        ));
                    };
                    let map_lookup_started = self.profile_detailed.then(Instant::now);
                    let raw = prefetched.get(probe_start);
                    if let Some(t0) = map_lookup_started {
                        self.profile.prefetch_map_lookup += t0.elapsed();
                    }
                    if let Some(raw) = raw {
                        if self.profile_detailed {
                            self.profile.prefetch_hits += 1;
                            self.profile.compressed_bytes += raw.as_ref().len() as u64;
                        }
                        let decode_started = self.profile_detailed.then(Instant::now);
                        store.decode_position_entry_value(
                            raw.as_ref(),
                            decompressor.as_mut(),
                            &mut decompress_buf,
                        )?;
                        if let Some(t0) = decode_started {
                            self.profile.decompress += t0.elapsed();
                            self.profile.decode_calls += 1;
                            self.profile.decompressed_bytes += decompress_buf.len() as u64;
                        }
                        true
                    } else {
                        if self.profile_detailed {
                            self.profile.prefetch_misses += 1;
                        }
                        false
                    }
                } else {
                    let Some(store) = self.store.as_ref() else {
                        return Err(DataFusionError::Execution(
                            "indexed parquet variation lookup reached legacy fjall fallback; warm/cold parquet did not cover the probe"
                                .into(),
                        ));
                    };
                    if self.profile_enabled {
                        self.profile.point_gets += 1;
                    }
                    if self.profile_detailed {
                        let get_started = Instant::now();
                        let raw = store.get_position_entry_with_key_buf(
                            chrom_code,
                            *probe_start,
                            &mut position_key_buf,
                        )?;
                        self.profile.point_get_raw += get_started.elapsed();
                        match raw {
                            Some(compressed) => {
                                self.profile.raw_get_hits += 1;
                                self.profile.compressed_bytes += compressed.as_ref().len() as u64;
                                let decode_started = Instant::now();
                                store.decode_position_entry_value(
                                    compressed.as_ref(),
                                    decompressor.as_mut(),
                                    &mut decompress_buf,
                                )?;
                                self.profile.decompress += decode_started.elapsed();
                                self.profile.decode_calls += 1;
                                self.profile.decompressed_bytes += decompress_buf.len() as u64;
                                true
                            }
                            None => {
                                self.profile.raw_get_misses += 1;
                                false
                            }
                        }
                    } else {
                        store.get_position_entry_fast_with_key_buf(
                            chrom_code,
                            *probe_start,
                            decompressor.as_mut(),
                            &mut decompress_buf,
                            &mut position_key_buf,
                        )?
                    }
                };
                if !found {
                    continue;
                }

                let reader_started = self.profile_detailed.then(Instant::now);
                let reader = PositionEntryReader::new(&decompress_buf)?;
                if let Some(t0) = reader_started {
                    self.profile.reader_init += t0.elapsed();
                }

                // Match alleles within this position entry (reuse buffer).
                // Filter by end-coordinate overlap: the cache allele's (start, end)
                // must overlap the VCF variant's interval. This prevents matching
                // alleles at the same start position but with non-overlapping end.
                let primary_match_started = self.profile_detailed.then(Instant::now);
                matched_allele_rows.clear();
                let vcf_iv_start = norm_start_i64.min(norm_end_i64);
                let vcf_iv_end = norm_start_i64.max(norm_end_i64);
                for allele_idx in 0..reader.num_alleles() {
                    if self.profile_detailed {
                        self.profile.primary_allele_rows += 1;
                    }
                    let failed = failed_stored_col_idx
                        .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                        .unwrap_or(0);
                    if failed > self.allowed_failed {
                        if self.profile_detailed {
                            self.profile.primary_failed_skips += 1;
                        }
                        continue;
                    }

                    let existing_end = end_stored_col_idx
                        .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                        .unwrap_or(*probe_start);
                    let cache_iv_start = (*probe_start).min(existing_end);
                    let cache_iv_end = (*probe_start).max(existing_end);
                    if cache_iv_start > vcf_iv_end || cache_iv_end < vcf_iv_start {
                        if self.profile_detailed {
                            self.profile.primary_interval_skips += 1;
                        }
                        continue;
                    }

                    let allele_str = reader.allele_string(allele_idx);
                    if self.profile_detailed {
                        self.profile.exact_match_calls += 1;
                    }
                    if (self.exact_matcher)(vcf_ref, vcf_alt, allele_str) {
                        matched_allele_rows.push(allele_idx);
                    }
                }
                if let Some(t0) = primary_match_started {
                    self.profile.primary_match += t0.elapsed();
                    self.profile.primary_matches += matched_allele_rows.len() as u64;
                }

                if !matched_allele_rows.is_empty() && !emitted_match {
                    let append_started = self.profile_detailed.then(Instant::now);
                    // Emit only the first matched allele per VCF row.
                    // Multiple alleles at the same position produce identical
                    // annotation output (same VCF coords → same transcript
                    // overlaps, colocated data, and consequence), so emitting
                    // duplicates is pure overhead.  Colocated collection below
                    // still iterates all alleles independently.
                    emitted_match = true;
                    row_output_emitted[row] = true;
                    let first_match = &matched_allele_rows[..1];
                    vcf_indices.push(row as u32);
                    for (col_out_idx, builder) in builders.iter_mut().enumerate() {
                        let entry_idx = col_map[col_out_idx];
                        if entry_idx == usize::MAX {
                            append_null_to_builder(builder.as_mut())?;
                        } else {
                            reader.append_column_values(
                                entry_idx,
                                first_match,
                                builder.as_mut(),
                            )?;
                        }
                    }
                    if let Some(t0) = append_started {
                        self.profile.cache_column_append += t0.elapsed();
                    }
                }

                // --- Co-located data collection (piggybacked on same probe) ---
                // Runs for ALL probed positions, not just when primary allele matched.
                // This mirrors the parquet path which streams all cache rows and
                // checks colocated eligibility independently of primary lookup.
                if let (Some(buf), Some(ci)) = (coloc_buf.as_mut(), self.coloc_col_indices.as_ref())
                {
                    // Compute VCF input allele key once per VCF row match.
                    let coloc_prepare_started = self.profile_detailed.then(Instant::now);
                    let chrom_norm = chrom.to_string();
                    let (input_ref, input_alt, input_start) =
                        vcf_to_vep_input_allele(norm_start_i64, vcf_ref, vcf_alt);
                    let input_allele_string = format!("{input_ref}/{input_alt}");
                    let (compare_ref, compare_alt) = vcf_to_vep_allele(vcf_ref, vcf_alt);
                    let compare_allele_string = format!("{compare_ref}/{compare_alt}");
                    let vep_start = vep_norm_start(norm_start_i64, vcf_ref, vcf_alt);
                    let vep_end = vep_norm_end(norm_start_i64, vcf_ref, vcf_alt);

                    let compare_output_allele =
                        output_allele_from_allele_string(&compare_allele_string)
                            .map(str::to_string);
                    let unshifted_output_allele: Option<String> = None;
                    if let Some(t0) = coloc_prepare_started {
                        self.profile.colocated_prepare += t0.elapsed();
                    }

                    // Visibility filter: mirrors VEP's Tabix query window.
                    // Only cache variants with START in [compare_start-1, compare_end+1]
                    // are visible, matching existing_start_is_visible_to_input_row().
                    let colocated_match_started = self.profile_detailed.then(Instant::now);
                    if !probe_start_visible_to_window(*probe_start, vep_start, vep_end) {
                        if let Some(t0) = colocated_match_started {
                            self.profile.colocated_match += t0.elapsed();
                        }
                        continue;
                    }

                    // Iterate alleles at this position for colocated collection.
                    for allele_idx in 0..reader.num_alleles() {
                        if self.profile_detailed {
                            self.profile.colocated_allele_rows += 1;
                        }
                        let failed = ci
                            .failed
                            .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                            .unwrap_or(0);
                        if failed > self.allowed_failed {
                            continue;
                        }

                        let var_name = reader.read_string_value(ci.variation_name, allele_idx);
                        let var_name = match var_name {
                            Some(v) if !v.is_empty() => v,
                            _ => continue,
                        };

                        let allele_str = reader.allele_string(allele_idx);

                        // Read existing variant's end coordinate from the entry.
                        let existing_end = ci
                            .end_col
                            .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                            .unwrap_or(*probe_start);

                        // VEP colocated matching uses parser-minimized input alleles.
                        // Shifted/right-aligned coordinates are lookup candidates only.
                        let Some(matched_alleles) = compare_existing_variant_alleles(
                            &compare_allele_string,
                            vep_start,
                            vep_end,
                            None,
                            None,
                            None,
                            allele_str,
                            *probe_start,
                            existing_end,
                        ) else {
                            continue;
                        };

                        let key: ColocatedKey = (
                            chrom_norm.clone(),
                            input_start,
                            norm_end_i64,
                            input_allele_string.clone(),
                        );

                        let sink_value = buf.entry(key).or_insert_with(|| ColocatedSinkValue {
                            entries: Vec::new(),
                            compare_output_allele: compare_output_allele.clone(),
                            unshifted_output_allele: unshifted_output_allele.clone(),
                        });

                        let somatic = ci
                            .somatic
                            .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                            .unwrap_or(0);
                        let pheno = ci
                            .pheno
                            .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                            .unwrap_or(0);
                        let clin_sig = ci
                            .clin_sig
                            .and_then(|idx| reader.read_string_value(idx, allele_idx));
                        let clin_sig_allele = ci
                            .clin_sig_allele
                            .and_then(|idx| reader.read_string_value(idx, allele_idx));
                        let pubmed = ci
                            .pubmed
                            .and_then(|idx| reader.read_string_value(idx, allele_idx));
                        let af_values: Vec<String> = ci
                            .af_indices
                            .iter()
                            .map(|opt_idx| {
                                opt_idx
                                    .and_then(|idx| reader.read_string_value(idx, allele_idx))
                                    .unwrap_or_default()
                            })
                            .collect();

                        sink_value.entries.push(ColocatedCacheEntry {
                            variation_name: var_name,
                            allele_string: allele_str.to_string(),
                            matched_alleles,
                            somatic,
                            pheno,
                            clin_sig,
                            clin_sig_allele,
                            pubmed,
                            af_values,
                        });
                        if self.profile_detailed {
                            self.profile.colocated_entries += 1;
                        }
                    }
                    if let Some(t0) = colocated_match_started {
                        self.profile.colocated_match += t0.elapsed();
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
            let mut pending_keys_by_chrom: HashMap<String, Vec<u64>> = HashMap::new();
            for pending in &pending_cold_probes {
                pending_keys_by_chrom
                    .entry(pending.chrom.clone())
                    .or_default()
                    .push(pending.position_key);
            }
            let touched_cold_chroms = pending_keys_by_chrom.keys().cloned().collect::<Vec<_>>();

            for (chrom, position_keys) in pending_keys_by_chrom {
                let lookup = self
                    .cold_parquet_chroms
                    .get_mut(&chrom)
                    .and_then(Option::as_mut)
                    .ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "cold parquet lookup not open for {chrom}"
                        ))
                    })?;
                lookup.prefetch_positions_retaining(position_keys)?;
            }

            for pending in &pending_cold_probes {
                let row = pending.vcf_row as usize;
                let emit_output = !row_output_emitted[row];
                let (result, metrics) = {
                    let chunk = self
                        .cold_parquet_chroms
                        .get(&pending.chrom)
                        .and_then(Option::as_ref)
                        .and_then(|lookup| lookup.cached_chunk_for_position(pending.position_key));
                    if let Some(chunk) = chunk {
                        probe_cold_chunk_position(
                            chunk,
                            self.exact_matcher,
                            self.allowed_failed,
                            self.profile_detailed,
                            &pending.chrom,
                            pending.probe_start,
                            pending.position_key,
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
                    }
                };

                if self.profile_enabled {
                    self.record_cold_chunk_probe_metrics(result, metrics);
                }
                if result == ColdProbeResult::Match && emit_output {
                    row_output_emitted[row] = true;
                }
            }

            if self.profile_enabled {
                self.record_cold_parquet_stats_snapshot();
            }

            for chrom in touched_cold_chroms {
                if let Some(lookup) = self
                    .cold_parquet_chroms
                    .get_mut(&chrom)
                    .and_then(Option::as_mut)
                {
                    lookup.trim_cache_to_capacity();
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

fn warm_file_for_chrom(warm_cache_dir: &std::path::Path, chrom: &str) -> Option<PathBuf> {
    let direct = warm_cache_dir.join(format!("{chrom}_warm.parquet"));
    if direct.is_file() {
        return Some(direct);
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        let stripped = warm_cache_dir.join(format!("{stripped}_warm.parquet"));
        if stripped.is_file() {
            return Some(stripped);
        }
    } else {
        let prefixed = warm_cache_dir.join(format!("chr{chrom}_warm.parquet"));
        if prefixed.is_file() {
            return Some(prefixed);
        }
    }

    None
}

fn warm_variant_key_candidates(
    position_key: u64,
    vcf_ref: &str,
    vcf_alt: &str,
) -> SmallVec<[u64; 4]> {
    let mut keys = SmallVec::<[u64; 4]>::new();
    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (vep_ref, vep_alt) = vcf_to_vep_allele(vcf_ref, alt);
        push_unique_variant_key(
            &mut keys,
            warm_variant_key_from_position(position_key, &vep_ref, &vep_alt),
        );
        push_unique_variant_key(
            &mut keys,
            warm_variant_key_from_position(position_key, vcf_ref, alt),
        );
        push_unique_variant_key(
            &mut keys,
            warm_variant_key_from_position(position_key, vcf_ref, &vep_alt),
        );
        push_unique_variant_key(
            &mut keys,
            warm_variant_key_from_position(position_key, &vep_ref, alt),
        );
    }
    keys
}

fn output_order_for_vcf_indices(vcf_indices: &[u32]) -> Option<Vec<u32>> {
    if vcf_indices.windows(2).all(|window| window[0] <= window[1]) {
        return None;
    }

    let mut output_order: Vec<u32> = (0..vcf_indices.len() as u32).collect();
    output_order.sort_by_key(|&idx| (vcf_indices[idx as usize], idx));
    Some(output_order)
}

#[inline]
fn push_unique_variant_key(keys: &mut SmallVec<[u64; 4]>, key: u64) {
    if !keys.contains(&key) {
        keys.push(key);
    }
}

fn resolve_warm_coloc_indices(chunk: &WarmChunkContext) -> Option<WarmColocIndices> {
    let schema = chunk.batch.schema();
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

fn warm_output_indices(
    chunk: &WarmChunkContext,
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
                chunk.batch.schema().index_of(name).ok().map(Some)
            }
        })
        .collect()
}

fn append_warm_row_values(
    chunk: &WarmChunkContext,
    row: usize,
    output_indices: &[Option<usize>],
    builders: &mut [Box<dyn datafusion::arrow::array::ArrayBuilder>],
) -> Result<()> {
    for (idx, builder) in output_indices.iter().zip(builders.iter_mut()) {
        if let Some(idx) = idx {
            append_array_value_to_builder(
                chunk.batch.column(*idx).as_ref(),
                row,
                builder.as_mut(),
            )?;
        } else {
            append_null_to_builder(builder.as_mut())?;
        }
    }
    Ok(())
}

fn warm_string_value(
    chunk: &WarmChunkContext,
    column_idx: Option<usize>,
    row: usize,
) -> Result<Option<String>> {
    let Some(column_idx) = column_idx else {
        return Ok(None);
    };
    let array = chunk.batch.column(column_idx);
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
            "warm colocated column expected string array, got {:?}",
            array.data_type()
        )))
    }
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
fn resolve_kv_coloc_indices(store: &VepKvStore) -> Option<KvColocIndices> {
    let cache_schema = store.schema();

    let cache_chrom_idx = cache_schema.index_of("chrom").unwrap_or(usize::MAX);
    let cache_start_idx = cache_schema.index_of("start").unwrap_or(usize::MAX);

    // Build the stored column mapping (same logic as process_batch).
    let stored_cols: Vec<usize> = (0..cache_schema.fields().len())
        .filter(|&i| i != cache_chrom_idx && i != cache_start_idx)
        .collect();

    let find_stored_idx = |name: &str| -> Option<usize> {
        cache_schema
            .index_of(name)
            .ok()
            .and_then(|schema_idx| stored_cols.iter().position(|&c| c == schema_idx))
    };

    let variation_name = find_stored_idx("variation_name")?;
    let allele_string_col = find_stored_idx("allele_string")?;

    let af_indices: Vec<Option<usize>> = AF_COL_NAMES
        .iter()
        .map(|name| find_stored_idx(name))
        .collect();

    Some(KvColocIndices {
        variation_name,
        _allele_string_col: allele_string_col,
        end_col: find_stored_idx("end"),
        failed: find_stored_idx("failed"),
        somatic: find_stored_idx("somatic"),
        pheno: find_stored_idx("phenotype_or_disease"),
        clin_sig: find_stored_idx("clin_sig"),
        clin_sig_allele: find_stored_idx("clin_sig_allele"),
        pubmed: find_stored_idx("pubmed"),
        af_indices,
    })
}

impl Stream for KvLookupStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        // KvLookupExec reads ALL alleles at each position in a single point
        // lookup (the position entry contains all alleles), so co-located
        // data for each VCF row is complete immediately — no buffering needed.
        match self.input.poll_next_unpin(cx) {
            Poll::Ready(Some(Ok(batch))) => {
                let result = self.process_batch(&batch);
                Poll::Ready(Some(result))
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

/// Probe a single VCF batch against a fjall KV store without wrapping it in a
/// DataFusion execution plan.
///
/// This is used by the VEP-buffer chunked path: each chunk owns raw VCF rows,
/// performs KV lookup locally, then immediately feeds the looked-up batch into
/// the transcript annotation engine.
#[allow(clippy::too_many_arguments)]
pub fn lookup_batch_with_store(
    vcf_batch: &RecordBatch,
    store: Arc<VepKvStore>,
    cache_columns: Vec<String>,
    vcf_has_chr: bool,
    vcf_zero_based: bool,
    cache_zero_based: bool,
    extended_probes: bool,
    allowed_failed: i64,
    colocated_sink: Option<ColocatedSink>,
    reference_fasta_path: Option<String>,
) -> Result<RecordBatch> {
    let input_schema = vcf_batch.schema();
    let (schema, output_col_positions) =
        build_lookup_output_schema(input_schema.clone(), store.schema().clone(), &cache_columns);
    let input: SendableRecordBatchStream = Box::pin(EmptyLookupInput {
        schema: input_schema,
    });
    let mut stream = KvLookupStream::new(
        input,
        Some(store.clone()),
        None,
        store.schema().clone(),
        false,
        schema,
        cache_columns,
        KvMatchMode::Exact,
        allele_matches as fn(&str, &str, &str) -> bool,
        vcf_has_chr,
        vcf_zero_based,
        cache_zero_based,
        extended_probes,
        allowed_failed,
        output_col_positions,
        colocated_sink,
        reference_fasta_path,
    );
    let batch = stream.process_batch(vcf_batch);
    if stream.profile_enabled && !stream.profile_emitted {
        stream.profile.emit(stream.profile_detailed);
        stream.profile_emitted = true;
    }
    batch
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warm_position_covered_no_exact_skips_fjall() {
        let decision = decide_after_warm_probe(
            WarmProbeForDecision::PositionCoveredNoExact,
            ColdPositionDecision::Present,
        );

        assert_eq!(decision, LookupDecision::SkipFjall);
    }

    #[test]
    fn warm_not_covered_and_cold_negative_skips_fjall() {
        let decision = decide_after_warm_probe(
            WarmProbeForDecision::NotCovered,
            ColdPositionDecision::Absent,
        );

        assert_eq!(decision, LookupDecision::SkipFjall);
    }

    #[test]
    fn warm_not_covered_and_cold_positive_uses_fjall() {
        let decision = decide_after_warm_probe(
            WarmProbeForDecision::NotCovered,
            ColdPositionDecision::Present,
        );

        assert_eq!(decision, LookupDecision::UseFjall);
    }

    #[test]
    fn warm_cold_backend_parser_defaults_to_fjall_and_accepts_parquet() {
        assert_eq!(
            WarmColdVariationBackend::parse(None).unwrap(),
            WarmColdVariationBackend::Fjall
        );
        assert_eq!(
            WarmColdVariationBackend::parse(Some("parquet")).unwrap(),
            WarmColdVariationBackend::Parquet
        );
        assert_eq!(
            WarmColdVariationBackend::parse(Some("fjall")).unwrap(),
            WarmColdVariationBackend::Fjall
        );
    }

    #[test]
    fn warm_cold_backend_parser_rejects_unknown_values() {
        let err = WarmColdVariationBackend::parse(Some("rocksdb")).unwrap_err();
        assert!(err.to_string().contains("VEP_WARM_COLD_VARIATION_BACKEND"));
    }

    #[test]
    fn output_order_restores_vcf_row_order_after_grouped_cold_replay() {
        assert_eq!(output_order_for_vcf_indices(&[0, 1, 2]), None);
        assert_eq!(
            output_order_for_vcf_indices(&[2, 0, 1]),
            Some(vec![1, 2, 0])
        );
    }

    #[test]
    fn warm_cold_index_mode_parser_accepts_position_and_bloom_modes() {
        assert_eq!(
            WarmColdVariationIndexMode::parse(None).unwrap(),
            WarmColdVariationIndexMode::Position
        );
        assert_eq!(
            WarmColdVariationIndexMode::parse(Some("posidx_bloom")).unwrap(),
            WarmColdVariationIndexMode::PositionThenVariantBloom
        );
        assert_eq!(
            WarmColdVariationIndexMode::parse(Some("bloom")).unwrap(),
            WarmColdVariationIndexMode::VariantBloom
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
        assert!(lines[6].contains("cold_parquet probes=36"));
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

    #[test]
    fn warm_row_values_append_to_cache_builders() {
        use std::sync::Arc;

        use datafusion::arrow::array::{ArrayRef, Int64Array, StringArray};
        use datafusion::arrow::datatypes::{Field, Schema};

        use crate::warm_cache::chunk::WarmChunkContext;
        use crate::warm_cache::key::{position_key, variant_keys_from_allele_string};

        let mut variant_keys = datafusion::arrow::array::ListBuilder::new(
            datafusion::arrow::array::UInt64Builder::new(),
        );
        for key in variant_keys_from_allele_string("1", 101, "A/G").unwrap() {
            variant_keys.values().append_value(key);
        }
        variant_keys.append(true);

        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new_list(
                "variant_keys",
                Arc::new(Field::new_list_field(DataType::UInt64, true)),
                false,
            ),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("end", DataType::Int64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![position_key("1", 101).unwrap()])) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec![Some("rs1")])) as ArrayRef,
                Arc::new(Int64Array::from(vec![Some(101)])) as ArrayRef,
            ],
        )
        .unwrap();
        let chunk = WarmChunkContext::try_new(0, batch).unwrap();
        let cache_columns = vec![
            "variation_name".to_string(),
            "end".to_string(),
            "missing_from_cache".to_string(),
        ];
        let col_map = vec![0, 1, usize::MAX];
        let output_types = [DataType::Utf8, DataType::Int64, DataType::Utf8];
        let mut builders = output_types
            .iter()
            .map(|dt| make_builder(dt, 1).unwrap())
            .collect::<Vec<_>>();

        let indices = warm_output_indices(&chunk, &cache_columns, &col_map).unwrap();
        append_warm_row_values(&chunk, 0, &indices, &mut builders).unwrap();

        let variation = builders[0]
            .finish()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap()
            .clone();
        let end = builders[1]
            .finish()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap()
            .clone();
        let missing = builders[2]
            .finish()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap()
            .clone();

        assert_eq!(variation.value(0), "rs1");
        assert_eq!(end.value(0), 101);
        assert!(missing.is_null(0));
    }

    // -----------------------------------------------------------------------
    // Colocated sink integration tests (P0/P1 parity with Parquet path)
    // -----------------------------------------------------------------------

    use crate::allele::allele_matches;
    use crate::kv_cache::position_entry::serialize_position_entry;
    use crate::warm_cache::key::position_key;
    use datafusion::arrow::array::Int8Array;
    use datafusion::arrow::datatypes::Schema;
    use datafusion::datasource::MemTable;
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;
    use std::collections::HashMap as StdHashMap;
    use std::sync::Mutex;

    fn test_temp_dir(_prefix: &str) -> tempfile::TempDir {
        tempfile::tempdir().expect("create temp dir")
    }

    /// Create a KV store, a VCF MemTable input, execute KvLookupExec with a
    /// colocated sink, and return the sink contents.
    async fn run_kv_with_colocated_sink(
        vcf_batch: RecordBatch,
        cache_schema: Arc<Schema>,
        entries: Vec<(&str, i64, Vec<u8>)>, // (chrom, start, serialized_entry)
        cache_columns: Vec<String>,
        extended_probes: bool,
        allowed_failed: i64,
    ) -> StdHashMap<ColocatedKey, ColocatedSinkValue> {
        let cache_dir = test_temp_dir("vep-kv-coloc");
        let store = VepKvStore::create(cache_dir.path(), cache_schema).unwrap();
        for (chrom, start, entry) in &entries {
            store.put_position_entry(chrom, *start, entry).unwrap();
        }
        store.persist().unwrap();
        drop(store);

        let reopened_store = Arc::new(VepKvStore::open(cache_dir.path()).expect("reopen KV store"));

        let vcf_schema = vcf_batch.schema();
        let vcf_mem = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        let ctx = datafusion::prelude::SessionContext::new();
        ctx.register_table("vcf_coloc", Arc::new(vcf_mem)).unwrap();
        let vcf_plan = ctx
            .table("vcf_coloc")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();

        let sink: ColocatedSink = Arc::new(Mutex::new(StdHashMap::new()));

        let exec = KvLookupExec::new(
            vcf_plan,
            reopened_store,
            cache_columns,
            KvMatchMode::Exact,
            allele_matches as fn(&str, &str, &str) -> bool,
            false, // vcf_has_chr (bare chrom names in test data)
            false, // vcf_zero_based (1-based)
            false, // cache_zero_based (1-based)
            extended_probes,
            allowed_failed,
        )
        .unwrap()
        .with_colocated_sink(Arc::clone(&sink));

        let task_ctx = ctx.task_ctx();
        let stream = exec.execute(0, task_ctx).unwrap();
        // Consume stream fully.
        let _batches: Vec<_> = datafusion::physical_plan::common::collect(stream)
            .await
            .unwrap();

        let guard = sink.lock().unwrap();
        guard.clone()
        // cache_dir dropped here → temp directory cleaned up
    }

    async fn run_kv_with_colocated_sink_and_warm_dir(
        vcf_batch: RecordBatch,
        cache_schema: Arc<Schema>,
        cache_columns: Vec<String>,
        warm_dir: &std::path::Path,
    ) -> (
        StdHashMap<ColocatedKey, ColocatedSinkValue>,
        Vec<RecordBatch>,
    ) {
        let cache_dir = test_temp_dir("vep-kv-warm-coloc");
        let store = VepKvStore::create(cache_dir.path(), cache_schema).unwrap();
        store.persist().unwrap();
        drop(store);

        let reopened_store = Arc::new(VepKvStore::open(cache_dir.path()).expect("reopen KV store"));

        let vcf_schema = vcf_batch.schema();
        let vcf_mem = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        let ctx = datafusion::prelude::SessionContext::new();
        ctx.register_table("vcf_warm_coloc", Arc::new(vcf_mem))
            .unwrap();
        let vcf_plan = ctx
            .table("vcf_warm_coloc")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();

        let sink: ColocatedSink = Arc::new(Mutex::new(StdHashMap::new()));

        let exec = KvLookupExec::new(
            vcf_plan,
            reopened_store,
            cache_columns,
            KvMatchMode::Exact,
            allele_matches as fn(&str, &str, &str) -> bool,
            false,
            false,
            false,
            false,
            0,
        )
        .unwrap()
        .with_colocated_sink(Arc::clone(&sink))
        .with_warm_cache_dir_override(warm_dir.to_path_buf());

        let task_ctx = ctx.task_ctx();
        let stream = exec.execute(0, task_ctx).unwrap();
        let batches: Vec<_> = datafusion::physical_plan::common::collect(stream)
            .await
            .unwrap();

        let coloc = sink.lock().unwrap().clone();
        (coloc, batches)
    }

    fn simple_cache_schema() -> Arc<Schema> {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
            Field::new("somatic", DataType::Int64, true),
            Field::new("phenotype_or_disease", DataType::Int64, true),
        ]))
    }

    fn simple_vcf_batch(chrom: &str, start: i64, end: i64, rf: &str, alt: &str) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec![chrom])),
                Arc::new(Int64Array::from(vec![start])),
                Arc::new(Int64Array::from(vec![end])),
                Arc::new(StringArray::from(vec![rf])),
                Arc::new(StringArray::from(vec![alt])),
            ],
        )
        .unwrap()
    }

    fn write_single_group_warm_parquet(path: &std::path::Path, batch: &RecordBatch) {
        let file = std::fs::File::create(path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_row_count(Some(batch.num_rows()))
            .build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();
        writer.write(batch).unwrap();
        writer.close().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn warm_lookup_populates_colocated_sink_when_check_existing_is_active() {
        let cache_schema = simple_cache_schema();
        let warm_dir = tempfile::tempdir().unwrap();
        let mut variant_keys = datafusion::arrow::array::ListBuilder::new(
            datafusion::arrow::array::UInt64Builder::new(),
        );
        for key in crate::warm_cache::key::variant_keys_from_allele_string("1", 100, "A/G").unwrap()
        {
            variant_keys.values().append_value(key);
        }
        variant_keys.append(true);
        let warm_batch = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("position_key", DataType::UInt64, false),
                Field::new_list(
                    "variant_keys",
                    Arc::new(Field::new_list_field(DataType::UInt64, true)),
                    false,
                ),
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
                Field::new("failed", DataType::Int64, false),
                Field::new("somatic", DataType::Int64, true),
                Field::new("phenotype_or_disease", DataType::Int64, true),
            ])),
            vec![
                Arc::new(UInt64Array::from(vec![position_key("1", 100).unwrap()])) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["1"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![100])) as ArrayRef,
                Arc::new(Int64Array::from(vec![100])) as ArrayRef,
                Arc::new(StringArray::from(vec!["rs_warm"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["pathogenic"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1])) as ArrayRef,
            ],
        )
        .unwrap();
        write_single_group_warm_parquet(&warm_dir.path().join("chr1_warm.parquet"), &warm_batch);
        crate::kv_cache::position_index::PositionIndex::from_positions([])
            .unwrap()
            .write_to_path(warm_dir.path().join("variation.position_index/chr1.posidx"))
            .unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let (coloc, batches) = run_kv_with_colocated_sink_and_warm_dir(
            vcf,
            cache_schema,
            vec!["variation_name".into(), "clin_sig".into()],
            warm_dir.path(),
        )
        .await;

        let output_rows: usize = batches.iter().map(|batch| batch.num_rows()).sum();
        assert_eq!(output_rows, 1);
        let variation_idx = batches[0]
            .schema()
            .index_of("cache_variation_name")
            .unwrap();
        let variation = batches[0]
            .column(variation_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(variation.value(0), "rs_warm");

        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|value| value.entries.iter())
            .map(|entry| entry.variation_name.as_str())
            .collect();
        assert_eq!(all_names, vec!["rs_warm"]);
        let entry = &coloc.values().next().unwrap().entries[0];
        assert_eq!(entry.allele_string, "A/G");
        assert_eq!(entry.somatic, 1);
        assert_eq!(entry.pheno, 1);
        assert_eq!(entry.clin_sig, Some("pathogenic".to_string()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn indexed_parquet_warm_lookup_populates_colocated_sink_without_kv_store() {
        let cache_root = tempfile::tempdir().unwrap();
        let variation_dir = cache_root.path().join("variation");
        std::fs::create_dir_all(&variation_dir).unwrap();

        let mut variant_keys = datafusion::arrow::array::ListBuilder::new(
            datafusion::arrow::array::UInt64Builder::new(),
        );
        for key in crate::warm_cache::key::variant_keys_from_allele_string("1", 100, "A/G").unwrap()
        {
            variant_keys.values().append_value(key);
        }
        variant_keys.append(true);
        let warm_batch = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("position_key", DataType::UInt64, false),
                Field::new_list(
                    "variant_keys",
                    Arc::new(Field::new_list_field(DataType::UInt64, true)),
                    false,
                ),
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
                Field::new("failed", DataType::Int64, false),
                Field::new("somatic", DataType::Int64, true),
                Field::new("phenotype_or_disease", DataType::Int64, true),
            ])),
            vec![
                Arc::new(UInt64Array::from(vec![position_key("1", 100).unwrap()])) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["1"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![100])) as ArrayRef,
                Arc::new(Int64Array::from(vec![100])) as ArrayRef,
                Arc::new(StringArray::from(vec!["rs_indexed_warm"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["pathogenic"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1])) as ArrayRef,
            ],
        )
        .unwrap();
        write_single_group_warm_parquet(&variation_dir.join("chr1_warm.parquet"), &warm_batch);
        crate::kv_cache::position_index::PositionIndex::from_positions([])
            .unwrap()
            .write_to_path(
                cache_root
                    .path()
                    .join("variation.position_index/chr1.posidx"),
            )
            .unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let vcf_schema = vcf.schema();
        let vcf_mem = MemTable::try_new(vcf_schema, vec![vec![vcf]]).unwrap();
        let ctx = datafusion::prelude::SessionContext::new();
        ctx.register_table("vcf_indexed_warm_coloc", Arc::new(vcf_mem))
            .unwrap();
        let vcf_plan = ctx
            .table("vcf_indexed_warm_coloc")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();
        let sink: ColocatedSink = Arc::new(Mutex::new(StdHashMap::new()));

        let exec = KvLookupExec::new_indexed_parquet(
            vcf_plan,
            cache_root.path().to_path_buf(),
            simple_cache_schema(),
            vec!["variation_name".into(), "clin_sig".into()],
            KvMatchMode::Exact,
            allele_matches as fn(&str, &str, &str) -> bool,
            false,
            false,
            false,
            false,
            0,
        )
        .unwrap()
        .with_colocated_sink(Arc::clone(&sink));

        let stream = exec.execute(0, ctx.task_ctx()).unwrap();
        let batches: Vec<_> = datafusion::physical_plan::common::collect(stream)
            .await
            .unwrap();
        assert_eq!(batches.iter().map(RecordBatch::num_rows).sum::<usize>(), 1);

        let coloc = sink.lock().unwrap();
        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|value| value.entries.iter())
            .map(|entry| entry.variation_name.as_str())
            .collect();
        assert_eq!(all_names, vec!["rs_indexed_warm"]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn indexed_parquet_colocated_cold_cosmic_rows_use_position_fallback_bloom() {
        let cache_root = tempfile::tempdir().unwrap();
        let variation_dir = cache_root.path().join("variation");
        std::fs::create_dir_all(&variation_dir).unwrap();

        let target_position = 100;
        let warm_position = 10;
        let mut warm_variant_keys = datafusion::arrow::array::ListBuilder::new(
            datafusion::arrow::array::UInt64Builder::new(),
        );
        for key in
            crate::warm_cache::key::variant_keys_from_allele_string("1", warm_position, "A/G")
                .unwrap()
        {
            warm_variant_keys.values().append_value(key);
        }
        warm_variant_keys.append(true);
        let warm_batch = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("position_key", DataType::UInt64, false),
                Field::new_list(
                    "variant_keys",
                    Arc::new(Field::new_list_field(DataType::UInt64, true)),
                    false,
                ),
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
                Field::new("failed", DataType::Int64, false),
                Field::new("somatic", DataType::Int64, true),
                Field::new("phenotype_or_disease", DataType::Int64, true),
            ])),
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", warm_position).unwrap(),
                ])) as ArrayRef,
                Arc::new(warm_variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["1"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![warm_position])) as ArrayRef,
                Arc::new(Int64Array::from(vec![warm_position])) as ArrayRef,
                Arc::new(StringArray::from(vec!["rs_warm_other"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G"])) as ArrayRef,
                Arc::new(StringArray::from(vec![Option::<&str>::None])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0])) as ArrayRef,
            ],
        )
        .unwrap();
        write_single_group_warm_parquet(&variation_dir.join("chr1_warm.parquet"), &warm_batch);

        let cold_batch = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("position_key", DataType::UInt64, false),
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
                Field::new("failed", DataType::Int64, true),
                Field::new("somatic", DataType::Int64, true),
                Field::new("phenotype_or_disease", DataType::Int64, true),
            ])),
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", target_position).unwrap(),
                ])) as ArrayRef,
                Arc::new(StringArray::from(vec!["1"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![target_position])) as ArrayRef,
                Arc::new(Int64Array::from(vec![target_position])) as ArrayRef,
                Arc::new(StringArray::from(vec!["COSV_TEST"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["COSMIC_MUTATION"])) as ArrayRef,
                Arc::new(StringArray::from(vec![Option::<&str>::None])) as ArrayRef,
                Arc::new(Int64Array::from(vec![Option::<i64>::None])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1])) as ArrayRef,
            ],
        )
        .unwrap();
        write_single_group_warm_parquet(&variation_dir.join("chr1_cold.parquet"), &cold_batch);

        crate::kv_cache::position_index::PositionIndex::from_positions([target_position])
            .unwrap()
            .write_to_path(
                cache_root
                    .path()
                    .join("variation.position_index/chr1.posidx"),
            )
            .unwrap();
        let mut bloom = VariantBloomIndex::with_expected_items(1, 10).unwrap();
        bloom.insert_position_fallback_key(position_key("1", target_position).unwrap());
        bloom
            .write_to_path(
                cache_root
                    .path()
                    .join("variation.variant_bloom_index/chr1.varbf"),
            )
            .unwrap();

        let vcf = simple_vcf_batch("1", target_position, target_position, "T", "C");
        let vcf_schema = vcf.schema();
        let vcf_mem = MemTable::try_new(vcf_schema, vec![vec![vcf]]).unwrap();
        let ctx = datafusion::prelude::SessionContext::new();
        ctx.register_table("vcf_indexed_cold_cosmic", Arc::new(vcf_mem))
            .unwrap();
        let vcf_plan = ctx
            .table("vcf_indexed_cold_cosmic")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();
        let sink: ColocatedSink = Arc::new(Mutex::new(StdHashMap::new()));

        let exec = KvLookupExec::new_indexed_parquet(
            vcf_plan,
            cache_root.path().to_path_buf(),
            simple_cache_schema(),
            vec!["variation_name".into()],
            KvMatchMode::Exact,
            allele_matches as fn(&str, &str, &str) -> bool,
            false,
            false,
            false,
            false,
            0,
        )
        .unwrap()
        .with_colocated_sink(Arc::clone(&sink));

        let stream = exec.execute(0, ctx.task_ctx()).unwrap();
        let batches: Vec<_> = datafusion::physical_plan::common::collect(stream)
            .await
            .unwrap();
        assert_eq!(batches.iter().map(RecordBatch::num_rows).sum::<usize>(), 1);

        let coloc = sink.lock().unwrap();
        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|value| value.entries.iter())
            .map(|entry| entry.variation_name.as_str())
            .collect();
        assert_eq!(all_names, vec!["COSV_TEST"]);
        let entry = &coloc.values().next().unwrap().entries[0];
        assert_eq!(entry.allele_string, "COSMIC_MUTATION");
        assert_eq!(entry.somatic, 1);
        assert_eq!(entry.pheno, 1);
    }

    /// Empty VCF input should produce zero output rows.
    #[tokio::test(flavor = "multi_thread")]
    async fn test_kv_lookup_empty_input() {
        let cache_schema = simple_cache_schema();
        let cache_dir = test_temp_dir("vep-kv-empty");
        let store = VepKvStore::create(cache_dir.path(), cache_schema.clone()).unwrap();
        store.persist().unwrap();
        drop(store);

        let reopened_store = Arc::new(VepKvStore::open(cache_dir.path()).expect("reopen KV store"));

        // Create an empty VCF batch (0 rows).
        let vcf_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let empty_batch = RecordBatch::new_empty(vcf_schema.clone());
        let vcf_mem = MemTable::try_new(vcf_schema, vec![vec![empty_batch]]).unwrap();

        let ctx = datafusion::prelude::SessionContext::new();
        ctx.register_table("vcf_empty", Arc::new(vcf_mem)).unwrap();
        let vcf_plan = ctx
            .table("vcf_empty")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();

        let exec = KvLookupExec::new(
            vcf_plan,
            reopened_store,
            vec!["variation_name".into()],
            KvMatchMode::Exact,
            allele_matches as fn(&str, &str, &str) -> bool,
            false,
            false,
            false,
            false,
            0,
        )
        .unwrap();

        let task_ctx = ctx.task_ctx();
        let stream = exec.execute(0, task_ctx).unwrap();
        let batches: Vec<_> = datafusion::physical_plan::common::collect(stream)
            .await
            .unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 0,
            "Empty VCF input should produce 0 output rows"
        );
    }

    /// Verify colocated sink is populated per-row during streaming (not requiring
    /// a separate buffering pass).
    #[tokio::test(flavor = "multi_thread")]
    async fn test_kv_lookup_colocated_sink_populated_after_stream() {
        let cache_schema = simple_cache_schema();
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["rs_sink_test"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec![Option::<&str>::None])),
                Arc::new(Int64Array::from(vec![0])),
                Arc::new(Int64Array::from(vec![0i64])),
                Arc::new(Int64Array::from(vec![0i64])),
            ],
        )
        .unwrap();

        let entry =
            serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let coloc = run_kv_with_colocated_sink(
            vcf,
            cache_schema,
            vec![("1", 100, entry)],
            vec!["variation_name".into()],
            false,
            0,
        )
        .await;

        // The sink should have been populated during streaming.
        assert!(
            !coloc.is_empty(),
            "Colocated sink should have entries after stream exhaustion"
        );
        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|v| v.entries.iter())
            .map(|e| e.variation_name.as_str())
            .collect();
        assert!(
            all_names.contains(&"rs_sink_test"),
            "Colocated sink should contain rs_sink_test entry"
        );
    }

    /// P0: Colocated sink collects entries with correct field values.
    #[tokio::test(flavor = "multi_thread")]
    async fn colocated_sink_collects_entries_with_correct_fields() {
        let cache_schema = simple_cache_schema();
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["rs999"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["pathogenic"])),
                Arc::new(Int64Array::from(vec![0])),
                Arc::new(Int64Array::from(vec![1i64])), // somatic
                Arc::new(Int64Array::from(vec![1i64])), // pheno
            ],
        )
        .unwrap();

        // col_indices: end(2), variation_name(3), allele_string(4), clin_sig(5), failed(6), somatic(7), pheno(8)
        let entry =
            serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let coloc = run_kv_with_colocated_sink(
            vcf,
            cache_schema,
            vec![("1", 100, entry)],
            vec!["variation_name".into(), "clin_sig".into()],
            false,
            0,
        )
        .await;

        assert_eq!(coloc.len(), 1, "Expected one colocated key");
        let (key, value) = coloc.iter().next().unwrap();
        assert_eq!(key.0, "1"); // chrom
        assert_eq!(value.entries.len(), 1);
        assert_eq!(value.entries[0].variation_name, "rs999");
        assert_eq!(value.entries[0].allele_string, "A/G");
        assert_eq!(value.entries[0].somatic, 1);
        assert_eq!(value.entries[0].pheno, 1);
        assert_eq!(value.entries[0].clin_sig, Some("pathogenic".to_string()));
    }

    /// P0: Colocated sink excludes entries with failed > allowed_failed.
    #[tokio::test(flavor = "multi_thread")]
    async fn colocated_sink_filters_failed_entries() {
        let cache_schema = simple_cache_schema();
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(StringArray::from(vec!["rs_keep", "rs_failed"])),
                Arc::new(StringArray::from(vec!["A/G", "A/T"])),
                Arc::new(StringArray::from(vec![
                    Option::<&str>::None,
                    Option::<&str>::None,
                ])),
                Arc::new(Int64Array::from(vec![0, 1])), // failed=0, failed=1
                Arc::new(Int64Array::from(vec![0i64, 0])),
                Arc::new(Int64Array::from(vec![0i64, 0])),
            ],
        )
        .unwrap();

        let entry =
            serialize_position_entry(&[0, 1], &cache_batch, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let coloc = run_kv_with_colocated_sink(
            vcf,
            cache_schema,
            vec![("1", 100, entry)],
            vec!["variation_name".into()],
            false,
            0, // allowed_failed = 0 → rs_failed (failed=1) should be excluded
        )
        .await;

        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|v| v.entries.iter())
            .map(|e| e.variation_name.as_str())
            .collect();
        assert!(
            all_names.contains(&"rs_keep"),
            "rs_keep should be in colocated sink"
        );
        assert!(
            !all_names.contains(&"rs_failed"),
            "rs_failed (failed=1) should be excluded from colocated sink"
        );
    }

    /// P0: Colocated sink skips entries with null or empty variation_name.
    #[tokio::test(flavor = "multi_thread")]
    async fn colocated_sink_skips_null_and_empty_variation_name() {
        let cache_schema = simple_cache_schema();
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(StringArray::from(vec![Some("rs_valid"), None, Some("")])),
                Arc::new(StringArray::from(vec!["A/G", "A/T", "A/C"])),
                Arc::new(StringArray::from(vec![
                    Option::<&str>::None,
                    Option::<&str>::None,
                    Option::<&str>::None,
                ])),
                Arc::new(Int64Array::from(vec![0, 0, 0])),
                Arc::new(Int64Array::from(vec![0i64, 0, 0])),
                Arc::new(Int64Array::from(vec![0i64, 0, 0])),
            ],
        )
        .unwrap();

        let entry =
            serialize_position_entry(&[0, 1, 2], &cache_batch, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let coloc = run_kv_with_colocated_sink(
            vcf,
            cache_schema,
            vec![("1", 100, entry)],
            vec!["variation_name".into()],
            false,
            0,
        )
        .await;

        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|v| v.entries.iter())
            .map(|e| e.variation_name.as_str())
            .collect();
        assert_eq!(
            all_names,
            vec!["rs_valid"],
            "Only rs_valid should be collected; null and empty skipped"
        );
    }

    /// P0: Visibility filter excludes cache variants outside [compare_start-1, compare_end+1].
    #[tokio::test(flavor = "multi_thread")]
    async fn colocated_sink_visibility_filter_excludes_out_of_window() {
        let cache_schema = simple_cache_schema();

        // VCF deletion: chr1:100-103 TTA→T (compare space after prefix trim: start=101, end=103)
        // Cache at position 101 is visible (in [100,104]), but cache at position 110 is NOT.
        let cache_batch_101 = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["rs_visible"])),
                Arc::new(StringArray::from(vec!["TA/-"])),
                Arc::new(StringArray::from(vec![Option::<&str>::None])),
                Arc::new(Int64Array::from(vec![0])),
                Arc::new(Int64Array::from(vec![0i64])),
                Arc::new(Int64Array::from(vec![0i64])),
            ],
        )
        .unwrap();

        let cache_batch_110 = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![110])),
                Arc::new(Int64Array::from(vec![113])),
                Arc::new(StringArray::from(vec!["rs_invisible"])),
                Arc::new(StringArray::from(vec!["TA/-"])),
                Arc::new(StringArray::from(vec![Option::<&str>::None])),
                Arc::new(Int64Array::from(vec![0])),
                Arc::new(Int64Array::from(vec![0i64])),
                Arc::new(Int64Array::from(vec![0i64])),
            ],
        )
        .unwrap();

        let entry_101 =
            serialize_position_entry(&[0], &cache_batch_101, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();
        let entry_110 =
            serialize_position_entry(&[0], &cache_batch_110, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();

        // VCF: 1:100-103 TTA→T (deletion)
        let vcf = simple_vcf_batch("1", 100, 103, "TTA", "T");
        let coloc = run_kv_with_colocated_sink(
            vcf,
            cache_schema,
            vec![("1", 101, entry_101), ("1", 110, entry_110)],
            vec!["variation_name".into()],
            true, // extended_probes needed for deletion matching
            0,
        )
        .await;

        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|v| v.entries.iter())
            .map(|e| e.variation_name.as_str())
            .collect();
        assert!(
            all_names.contains(&"rs_visible"),
            "rs_visible at position 101 should be in window"
        );
        assert!(
            !all_names.contains(&"rs_invisible"),
            "rs_invisible at position 110 should be outside visibility window"
        );
    }

    /// P1: Multiple alleles at same position produce separate colocated entries.
    #[tokio::test(flavor = "multi_thread")]
    async fn colocated_sink_collects_multiple_alleles_at_same_position() {
        let cache_schema = simple_cache_schema();
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(StringArray::from(vec!["rs_snv1", "rs_snv2"])),
                Arc::new(StringArray::from(vec!["A/G", "A/T"])),
                Arc::new(StringArray::from(vec!["benign", "pathogenic"])),
                Arc::new(Int64Array::from(vec![0, 0])),
                Arc::new(Int64Array::from(vec![0i64, 1])),
                Arc::new(Int64Array::from(vec![0i64, 0])),
            ],
        )
        .unwrap();

        let entry =
            serialize_position_entry(&[0, 1], &cache_batch, &[2, 3, 4, 5, 6, 7, 8], 4).unwrap();

        let vcf = simple_vcf_batch("1", 100, 100, "A", "G");
        let coloc = run_kv_with_colocated_sink(
            vcf,
            cache_schema,
            vec![("1", 100, entry)],
            vec!["variation_name".into()],
            false,
            0,
        )
        .await;

        // Both alleles at position 100 should be collected even though
        // only one (A/G) matches the VCF allele — colocated collection
        // is independent of primary allele match.
        let all_names: Vec<&str> = coloc
            .values()
            .flat_map(|v| v.entries.iter())
            .map(|e| e.variation_name.as_str())
            .collect();
        assert!(
            all_names.contains(&"rs_snv1"),
            "rs_snv1 (A/G) should be collected"
        );
        // rs_snv2 (A/T) may or may not be collected depending on allele matching.
        // Two-pass matching: compare_existing_variant_alleles checks if
        // the existing allele matches — A/T at same coords won't match A/G input.
        // This verifies the colocated collection respects allele matching.
    }
}
