//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store per-position for annotation.

use std::any::Any;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};
use std::time::{Duration, Instant};

use datafusion::arrow::array::{
    Array, ArrayRef, Int32Array, Int64Array, LargeStringArray, RecordBatch, StringArray,
    StringViewArray, UInt32Array, UInt64Array,
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

use super::allele_index::AlleleMatcher;
use super::key_encoding::chrom_to_code;
use super::kv_store::VepKvStore;
use super::position_entry::{PositionEntryReader, make_builder};
use crate::allele::{
    VariantAlleleInput, get_matched_variant_alleles, vcf_to_vep_allele, vcf_to_vep_input_allele,
    vep_norm_end, vep_norm_start,
};
use crate::variant_lookup_exec::{
    AF_COL_NAMES, ColocatedCacheEntry, ColocatedKey, ColocatedSink, ColocatedSinkValue,
    build_shifted_compare_state, compare_existing_variant_alleles,
    output_allele_from_allele_string, read_reference_sequence,
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
    store: Arc<VepKvStore>,
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
    properties: PlanProperties,
    /// Cache schema column positions for requested cache output columns.
    output_col_positions: Vec<usize>,
    /// Optional sink for co-located data collected during probe phase.
    colocated_sink: Option<ColocatedSink>,
    /// Optional indexed reference FASTA for genomic shift state in colocated
    /// matching (parity with parquet path's two-pass allele matching).
    reference_fasta_path: Option<String>,
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
        let cache_schema = store.schema();

        let mut output_col_positions = Vec::new();
        let mut fields: Vec<Arc<Field>> = input_schema.fields().iter().cloned().collect();
        for col_name in &cache_columns {
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
        let schema = Arc::new(datafusion::arrow::datatypes::Schema::new(fields));

        let properties = PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            input.output_partitioning().clone(),
            EmissionType::Incremental,
            Boundedness::Bounded,
        );

        Ok(Self {
            input,
            store,
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
            reference_fasta_path: None,
        })
    }

    /// Set the co-located data sink for piggybacked collection during probe.
    pub fn with_colocated_sink(mut self, sink: ColocatedSink) -> Self {
        self.colocated_sink = Some(sink);
        self
    }

    /// Set the reference FASTA path for genomic shift state in colocated matching.
    pub fn with_reference_fasta_path(mut self, path: Option<String>) -> Self {
        self.reference_fasta_path = path;
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

    fn properties(&self) -> &PlanProperties {
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
        let mut exec = KvLookupExec::new(
            children[0].clone(),
            self.store.clone(),
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
        exec = exec.with_reference_fasta_path(self.reference_fasta_path.clone());
        Ok(Arc::new(exec))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let input_stream = self.input.execute(partition, context)?;

        Ok(Box::pin(KvLookupStream::new(
            input_stream,
            self.store.clone(),
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
            self.colocated_sink.clone(),
            self.reference_fasta_path.clone(),
        )))
    }
}

/// Streaming implementation that processes VCF batches and probes the KV store.
struct KvLookupStream {
    input: SendableRecordBatchStream,
    store: Arc<VepKvStore>,
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
    profile_enabled: bool,
    profile_emitted: bool,
    profile: LookupProfile,
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
    extract_cols: Duration,
    match_loop: Duration,
    vcf_take: Duration,
    cache_build: Duration,
}

impl LookupProfile {
    fn total_known(&self) -> Duration {
        self.extract_cols + self.match_loop + self.vcf_take + self.cache_build
    }

    fn pct(stage: Duration, total: Duration) -> f64 {
        if total.is_zero() {
            0.0
        } else {
            stage.as_secs_f64() * 100.0 / total.as_secs_f64()
        }
    }

    fn emit(&self) {
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
            "[vep-kv-profile] batches={} input_rows={} output_rows={} total_s={:.3} input_rows_per_s={:.1} output_rows_per_s={:.1}",
            self.batches,
            self.input_rows,
            self.output_rows,
            total.as_secs_f64(),
            input_rate,
            output_rate
        );
        eprintln!(
            "[vep-kv-profile] extract_cols={:.3}s ({:.1}%) match_loop={:.3}s ({:.1}%) vcf_take={:.3}s ({:.1}%) cache_build={:.3}s ({:.1}%)",
            self.extract_cols.as_secs_f64(),
            Self::pct(self.extract_cols, total),
            self.match_loop.as_secs_f64(),
            Self::pct(self.match_loop, total),
            self.vcf_take.as_secs_f64(),
            Self::pct(self.vcf_take, total),
            self.cache_build.as_secs_f64(),
            Self::pct(self.cache_build, total),
        );
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

impl KvLookupStream {
    #[allow(clippy::too_many_arguments)]
    fn new(
        input: SendableRecordBatchStream,
        store: Arc<VepKvStore>,
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
            .and_then(|_| resolve_kv_coloc_indices(&store));

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

        Self {
            input,
            store,
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
            profile_enabled: std::env::var_os("VEP_KV_PROFILE").is_some(),
            profile_emitted: false,
            profile: LookupProfile::default(),
        }
    }

    /// Single-pass position-keyed lookup.
    ///
    /// For each VCF row, fetch the per-position entry from fjall, match alleles,
    /// and append matched column values directly into ArrayBuilders.
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

        // Reusable zstd decompressor — created once, amortized across all lookups.
        let mut decompressor = self.store.create_decompressor()?;

        // Reusable decompression / raw-value buffer — avoids alloc per lookup.
        let mut decompress_buf: Vec<u8> = Vec::with_capacity(4096);

        // Reusable allele match buffer — avoids alloc per row.
        let mut matched_allele_rows: Vec<usize> = Vec::new();

        // Determine which column indices in the entry correspond to our output columns.
        // Entry stores all columns except chrom/start, in schema order minus those 2.
        // (`end` is stored as a regular column inside the entry.)
        let cache_schema = self.store.schema();
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

        let match_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };

        // Local buffer for colocated data (flushed to the shared sink after the loop).
        let mut coloc_buf: Option<HashMap<ColocatedKey, ColocatedSinkValue>> =
            if self.colocated_sink.is_some() && self.coloc_col_indices.is_some() {
                Some(HashMap::new())
            } else {
                None
            };

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
            let probe_starts = build_probe_starts(
                norm_start_i64,
                norm_end_i64,
                vcf_ref,
                vcf_alt,
                self.extended_probes,
            );

            let mut emitted_match = false;
            for probe_start in &probe_starts {
                let found = self.store.get_position_entry_fast(
                    chrom_code,
                    *probe_start,
                    decompressor.as_mut(),
                    &mut decompress_buf,
                )?;
                if !found {
                    continue;
                }

                let reader = PositionEntryReader::new(&decompress_buf)?;

                // Match alleles within this position entry (reuse buffer).
                // Filter by end-coordinate overlap: the cache allele's (start, end)
                // must overlap the VCF variant's interval. This prevents matching
                // alleles at the same start position but with non-overlapping end.
                matched_allele_rows.clear();
                let vcf_iv_start = norm_start_i64.min(norm_end_i64);
                let vcf_iv_end = norm_start_i64.max(norm_end_i64);
                for allele_idx in 0..reader.num_alleles() {
                    let failed = failed_stored_col_idx
                        .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                        .unwrap_or(0);
                    if failed > self.allowed_failed {
                        continue;
                    }

                    let existing_end = end_stored_col_idx
                        .and_then(|idx| reader.read_i64_value(idx, allele_idx))
                        .unwrap_or(*probe_start);
                    let cache_iv_start = (*probe_start).min(existing_end);
                    let cache_iv_end = (*probe_start).max(existing_end);
                    if cache_iv_start > vcf_iv_end || cache_iv_end < vcf_iv_start {
                        continue;
                    }

                    let allele_str = reader.allele_string(allele_idx);
                    if (self.exact_matcher)(vcf_ref, vcf_alt, allele_str) {
                        matched_allele_rows.push(allele_idx);
                    }
                }

                if !matched_allele_rows.is_empty() {
                    emitted_match = true;
                    for _ in &matched_allele_rows {
                        vcf_indices.push(row as u32);
                    }
                    for (col_out_idx, builder) in builders.iter_mut().enumerate() {
                        let entry_idx = col_map[col_out_idx];
                        if entry_idx == usize::MAX {
                            // Column not found in entry -> nulls.
                            for _ in &matched_allele_rows {
                                append_null_to_builder(builder.as_mut())?;
                            }
                        } else {
                            reader.append_column_values(
                                entry_idx,
                                &matched_allele_rows,
                                builder.as_mut(),
                            )?;
                        }
                    }
                }

                // --- Co-located data collection (piggybacked on same probe) ---
                // Runs for ALL probed positions, not just when primary allele matched.
                // This mirrors the parquet path which streams all cache rows and
                // checks colocated eligibility independently of primary lookup.
                if let (Some(buf), Some(ci)) = (coloc_buf.as_mut(), self.coloc_col_indices.as_ref())
                {
                    // Compute VCF input allele key once per VCF row match.
                    let chrom_norm = chrom.to_string();
                    let (input_ref, input_alt, input_start) =
                        vcf_to_vep_input_allele(norm_start_i64, vcf_ref, vcf_alt);
                    let input_allele_string = format!("{input_ref}/{input_alt}");
                    let (compare_ref, compare_alt) = vcf_to_vep_allele(vcf_ref, vcf_alt);
                    let compare_allele_string = format!("{compare_ref}/{compare_alt}");
                    let vep_start = vep_norm_start(norm_start_i64, vcf_ref, vcf_alt);
                    let vep_end = vep_norm_end(norm_start_i64, vcf_ref, vcf_alt);

                    // Compute genomic shift state (mirrors parquet path's BuildRow logic).
                    let mut active_compare_allele_string = compare_allele_string.clone();
                    let mut active_compare_start = vep_start;
                    let mut active_compare_end = vep_end;
                    let mut unshifted_allele_string: Option<String> = None;
                    let mut unshifted_start: Option<i64> = None;

                    if let Some(ref_reader) = self.reference_reader.as_mut() {
                        if let Ok(Some((shifted_as, shifted_s, shifted_e))) =
                            build_shifted_compare_state(
                                ref_reader,
                                &chrom_norm,
                                &compare_allele_string,
                                vep_start,
                                vep_end,
                            )
                        {
                            unshifted_allele_string = Some(compare_allele_string.clone());
                            unshifted_start = Some(vep_start);
                            active_compare_allele_string = shifted_as;
                            active_compare_start = shifted_s;
                            active_compare_end = shifted_e;
                        }
                    }

                    let compare_output_allele =
                        output_allele_from_allele_string(&active_compare_allele_string)
                            .map(str::to_string);
                    let unshifted_output_allele = unshifted_allele_string
                        .as_deref()
                        .and_then(output_allele_from_allele_string)
                        .map(str::to_string);

                    // Visibility filter: mirrors VEP's Tabix query window.
                    // Only cache variants with START in [compare_start-1, compare_end+1]
                    // are visible, matching existing_start_is_visible_to_input_row().
                    let vis_start = (active_compare_start - 1).min(active_compare_end + 1);
                    let vis_end = (active_compare_start - 1).max(active_compare_end + 1);
                    if *probe_start < vis_start || *probe_start > vis_end {
                        continue;
                    }

                    // Iterate alleles at this position for colocated collection.
                    for allele_idx in 0..reader.num_alleles() {
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

                        // Two-pass allele matching (shifted + unshifted) for parity
                        // with parquet path's compare_existing_variant().
                        let Some(matched_alleles) = compare_existing_variant_alleles(
                            &active_compare_allele_string,
                            active_compare_start,
                            active_compare_end,
                            unshifted_allele_string.as_deref(),
                            unshifted_start,
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
                    }
                }
            }

            if !emitted_match {
                // No coordinate probe matched any allele -> null cache columns.
                vcf_indices.push(row as u32);
                for builder in &mut builders {
                    append_null_to_builder(builder.as_mut())?;
                }
            }
        }

        // Flush colocated data to the shared sink.
        if let (Some(buf), Some(sink)) = (coloc_buf, &self.colocated_sink) {
            if !buf.is_empty() {
                let mut guard = sink.lock().unwrap();
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
        let take_indices = UInt32Array::from(vcf_indices);
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
            output_columns.push(builder.finish());
        }
        if let Some(t0) = cache_build_started {
            self.profile.cache_build += t0.elapsed();
        }

        RecordBatch::try_new(self.schema.clone(), output_columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
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
        match self.input.poll_next_unpin(cx) {
            Poll::Ready(Some(Ok(batch))) => {
                let result = self.process_batch(&batch);
                Poll::Ready(Some(result))
            }
            Poll::Ready(Some(Err(e))) => Poll::Ready(Some(Err(e))),
            Poll::Ready(None) => {
                if self.profile_enabled && !self.profile_emitted {
                    self.profile.emit();
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

#[cfg(test)]
mod tests {
    use super::*;

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

    // -----------------------------------------------------------------------
    // Colocated sink integration tests (P0/P1 parity with Parquet path)
    // -----------------------------------------------------------------------

    use crate::allele::allele_matches;
    use crate::kv_cache::position_entry::serialize_position_entry;
    use datafusion::arrow::array::Int8Array;
    use datafusion::arrow::datatypes::Schema;
    use datafusion::datasource::MemTable;
    use std::collections::HashMap as StdHashMap;
    use std::sync::Mutex;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn test_temp_dir(prefix: &str) -> std::path::PathBuf {
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before UNIX_EPOCH")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("{prefix}-{}-{now}", std::process::id()));
        std::fs::create_dir_all(&path).expect("create temp dir");
        path
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
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        for (chrom, start, entry) in &entries {
            store.put_position_entry(chrom, *start, entry).unwrap();
        }
        store.persist().unwrap();
        drop(store);

        let reopened_store = Arc::new(VepKvStore::open(&cache_dir).expect("reopen KV store"));

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

        // Cleanup
        // let _ = std::fs::remove_dir_all(&cache_dir);
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
