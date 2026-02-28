//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store per-position for annotation.

use std::any::Any;
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

use crate::allele_index::AlleleMatcher;
use crate::key_encoding::chrom_to_code;
use crate::kv_store::VepKvStore;
use crate::position_entry::{PositionEntryReader, make_builder};

/// Lookup match mode (mirrors MatchMode from bio-function-vep).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KvMatchMode {
    /// Exact allele matching only.
    Exact,
    /// Exact first, then positional co-location fallback.
    ExactOrColocated,
    /// Exact first, then relaxed (indel-aware) fallback.
    ExactOrRelaxed,
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
    relaxed_matcher: Option<AlleleMatcher>,
    schema: SchemaRef,
    vcf_has_chr: bool,
    vcf_zero_based: bool,
    cache_zero_based: bool,
    /// When true, probe multiple coordinate encodings (insertion-style,
    /// shifted deletions, tandem repeat window). When false, probe only
    /// the exact normalized interval.
    extended_probes: bool,
    properties: PlanProperties,
    /// Cache schema column positions for requested cache output columns.
    output_col_positions: Vec<usize>,
}

impl KvLookupExec {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        input: Arc<dyn ExecutionPlan>,
        store: Arc<VepKvStore>,
        cache_columns: Vec<String>,
        match_mode: KvMatchMode,
        exact_matcher: AlleleMatcher,
        relaxed_matcher: Option<AlleleMatcher>,
        vcf_has_chr: bool,
        vcf_zero_based: bool,
        cache_zero_based: bool,
        extended_probes: bool,
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
            relaxed_matcher,
            schema,
            vcf_has_chr,
            vcf_zero_based,
            cache_zero_based,
            extended_probes,
            properties,
            output_col_positions,
        })
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
        Ok(Arc::new(KvLookupExec::new(
            children[0].clone(),
            self.store.clone(),
            self.cache_columns.clone(),
            self.match_mode,
            self.exact_matcher,
            self.relaxed_matcher,
            self.vcf_has_chr,
            self.vcf_zero_based,
            self.cache_zero_based,
            self.extended_probes,
        )?))
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
            self.relaxed_matcher,
            self.vcf_has_chr,
            self.vcf_zero_based,
            self.cache_zero_based,
            self.extended_probes,
            self.output_col_positions.clone(),
        )))
    }
}

/// Streaming implementation that processes VCF batches and probes the KV store.
struct KvLookupStream {
    input: SendableRecordBatchStream,
    store: Arc<VepKvStore>,
    schema: SchemaRef,
    cache_columns: Vec<String>,
    match_mode: KvMatchMode,
    exact_matcher: AlleleMatcher,
    relaxed_matcher: Option<AlleleMatcher>,
    vcf_has_chr: bool,
    vcf_zero_based: bool,
    cache_zero_based: bool,
    extended_probes: bool,
    output_col_positions: Vec<usize>,
    profile_enabled: bool,
    profile_emitted: bool,
    profile: LookupProfile,
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
        relaxed_matcher: Option<AlleleMatcher>,
        vcf_has_chr: bool,
        vcf_zero_based: bool,
        cache_zero_based: bool,
        extended_probes: bool,
        output_col_positions: Vec<usize>,
    ) -> Self {
        Self {
            input,
            store,
            schema,
            cache_columns,
            match_mode,
            exact_matcher,
            relaxed_matcher,
            vcf_has_chr,
            vcf_zero_based,
            cache_zero_based,
            extended_probes,
            output_col_positions,
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
        // Entry stores all columns except chrom/start/end, in schema order minus those 3.
        let cache_schema = self.store.schema();
        let cache_chrom_idx = cache_schema.index_of("chrom").unwrap_or(usize::MAX);
        let cache_start_idx = cache_schema.index_of("start").unwrap_or(usize::MAX);
        let cache_end_idx = cache_schema.index_of("end").unwrap_or(usize::MAX);

        // Build mapping: output_col_positions[i] -> index within the entry's column list.
        let stored_cols: Vec<usize> = (0..cache_schema.fields().len())
            .filter(|&i| i != cache_chrom_idx && i != cache_start_idx && i != cache_end_idx)
            .collect();
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

            // Probe a small set of coordinate encodings used by VEP-style caches:
            // - exact normalized interval
            // - indel point encodings at interval boundaries
            // - insertion-style start>end form for point variants
            let mut probe_keys: Vec<(i64, i64)> = Vec::with_capacity(4);
            probe_keys.push((norm_start_i64, norm_end_i64));
            if self.extended_probes {
                if norm_start == norm_end {
                    let shifted = i64::from(norm_start.saturating_add(1));
                    if !probe_keys.contains(&(shifted, norm_start_i64)) {
                        probe_keys.push((shifted, norm_start_i64));
                    }
                } else {
                    if !probe_keys.contains(&(norm_end_i64, norm_end_i64)) {
                        probe_keys.push((norm_end_i64, norm_end_i64));
                    }
                    if !probe_keys.contains(&(norm_start_i64, norm_start_i64)) {
                        probe_keys.push((norm_start_i64, norm_start_i64));
                    }
                    if !probe_keys.contains(&(norm_end_i64, norm_start_i64)) {
                        probe_keys.push((norm_end_i64, norm_start_i64));
                    }
                }
                // Probe prefix-trimmed coordinates used by VEP allele normalization
                // (e.g. REF=TTA ALT=T -> cache allele TA/- at shifted start).
                for alt in vcf_alt.split(['|', ',']).filter(|a| !a.is_empty()) {
                    let shift_usize = common_prefix_len(vcf_ref, alt);
                    if shift_usize == 0 {
                        continue;
                    }
                    // Apply shifted probes only to deletion-like events. For insertions,
                    // probing shifted point keys produces false positives.
                    let ref_remaining = vcf_ref.len().saturating_sub(shift_usize);
                    let alt_remaining = alt.len().saturating_sub(shift_usize);
                    if ref_remaining <= alt_remaining {
                        continue;
                    }
                    let shift = shift_usize as i64;
                    if let Some(shifted_start) = norm_start_i64.checked_add(shift) {
                        if !probe_keys.contains(&(shifted_start, norm_end_i64)) {
                            probe_keys.push((shifted_start, norm_end_i64));
                        }
                        if !probe_keys.contains(&(shifted_start, shifted_start)) {
                            probe_keys.push((shifted_start, shifted_start));
                        }
                        if !probe_keys.contains(&(norm_end_i64, shifted_start)) {
                            probe_keys.push((norm_end_i64, shifted_start));
                        }
                    }
                }
                // Deletions in tandem repeats may be right/left shifted in cache coordinates.
                // Probe a bounded window of equivalent deletion intervals that still overlap
                // the normalized VCF interval.
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
                            let Some(candidate_end) =
                                candidate_start.checked_add(del_len - 1)
                            else {
                                continue;
                            };
                            // Mirror SQL interval semantics: only consider intervals overlapping
                            // the normalized VCF interval.
                            if candidate_start > norm_end_i64 || candidate_end < norm_start_i64 {
                                continue;
                            }
                            if !probe_keys.contains(&(candidate_start, candidate_end)) {
                                probe_keys.push((candidate_start, candidate_end));
                            }
                        }
                    }
                }
            }

            let mut emitted_match = false;
            for (probe_start, probe_end) in probe_keys {
                let found = self.store.get_position_entry_fast(
                    chrom_code,
                    probe_start,
                    probe_end,
                    decompressor.as_mut(),
                    &mut decompress_buf,
                )?;
                if !found {
                    continue;
                }

                let reader = PositionEntryReader::new(&decompress_buf)?;

                // Match alleles within this position entry (reuse buffer).
                matched_allele_rows.clear();
                match self.match_mode {
                    KvMatchMode::Exact => {
                        for allele_idx in 0..reader.num_alleles() {
                            let allele_str = reader.allele_string(allele_idx);
                            if (self.exact_matcher)(vcf_ref, vcf_alt, allele_str) {
                                matched_allele_rows.push(allele_idx);
                            }
                        }
                    }
                    KvMatchMode::ExactOrColocated => {
                        for allele_idx in 0..reader.num_alleles() {
                            let allele_str = reader.allele_string(allele_idx);
                            if (self.exact_matcher)(vcf_ref, vcf_alt, allele_str) {
                                matched_allele_rows.push(allele_idx);
                            }
                        }
                        if matched_allele_rows.is_empty() {
                            // Colocated fallback: all alleles at this position.
                            matched_allele_rows.extend(0..reader.num_alleles());
                        }
                    }
                    KvMatchMode::ExactOrRelaxed => {
                        for allele_idx in 0..reader.num_alleles() {
                            let allele_str = reader.allele_string(allele_idx);
                            if (self.exact_matcher)(vcf_ref, vcf_alt, allele_str) {
                                matched_allele_rows.push(allele_idx);
                            }
                        }
                        if matched_allele_rows.is_empty() {
                            if let Some(relaxed) = self.relaxed_matcher {
                                for allele_idx in 0..reader.num_alleles() {
                                    let allele_str = reader.allele_string(allele_idx);
                                    if relaxed(vcf_ref, vcf_alt, allele_str) {
                                        matched_allele_rows.push(allele_idx);
                                    }
                                }
                            }
                        }
                    }
                }

                if matched_allele_rows.is_empty() {
                    continue;
                }

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

            if !emitted_match {
                // No coordinate probe matched any allele -> null cache columns.
                vcf_indices.push(row as u32);
                for builder in &mut builders {
                    append_null_to_builder(builder.as_mut())?;
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
