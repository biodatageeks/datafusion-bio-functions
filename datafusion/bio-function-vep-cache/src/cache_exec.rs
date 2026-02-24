//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store window-by-window for annotation.

use std::any::Any;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};
use std::time::{Duration, Instant};

use datafusion::arrow::array::builder::{
    BooleanBuilder, Float32Builder, Float64Builder, Int8Builder, Int16Builder, Int32Builder,
    Int64Builder, StringBuilder, UInt8Builder, UInt16Builder, UInt32Builder, UInt64Builder,
};
use datafusion::arrow::array::{
    Array, ArrayRef, BooleanArray, Float32Array, Float64Array, Int8Array, Int16Array, Int32Array,
    Int64Array, LargeStringArray, RecordBatch, StringArray, StringViewArray, UInt8Array,
    UInt16Array, UInt32Array, UInt64Array,
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

use crate::allele_index::{AlleleMatcher, WindowAlleleIndex};
use crate::key_encoding::window_id_for_position;
use crate::kv_store::{FORMAT_V1, VepKvStore, to_v1_column_index, validate_v1_schema_width};

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
/// Takes a sorted VCF input plan, probes a fjall KV store per-window,
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
    properties: PlanProperties,
    /// For v1 format: which cache column indices to fetch for output.
    output_col_indices: Vec<u8>,
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
    ) -> Result<Self> {
        let input_schema = input.schema();
        let cache_schema = store.schema();
        if store.format_version() != FORMAT_V1 {
            return Err(DataFusionError::Execution(format!(
                "unsupported cache format version {}: only v1 is supported",
                store.format_version()
            )));
        }
        validate_v1_schema_width(cache_schema.fields().len())?;

        // Compute output column indices for v1 format.
        let mut output_col_indices = Vec::new();
        let mut fields: Vec<Arc<Field>> = input_schema.fields().iter().cloned().collect();
        for col_name in &cache_columns {
            if let Ok(field) = cache_schema.field_with_name(col_name) {
                fields.push(Arc::new(Field::new(
                    format!("cache_{}", field.name()),
                    normalize_cache_output_type(field.data_type()),
                    true,
                )));
                if let Ok(idx) = cache_schema.index_of(col_name) {
                    output_col_indices.push(to_v1_column_index(idx)?);
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
            properties,
            output_col_indices,
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
            self.store.window_size(),
            self.output_col_indices.clone(),
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
    window_size: u64,
    output_col_indices: Vec<u8>,
    current_window: Option<(String, u64)>,
    current_index: Option<WindowAlleleIndex>,
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
    fetch_columns: Duration,
    vcf_take: Duration,
    cache_build: Duration,
    record_batch_build: Duration,
    window_loads: u64,
    window_load_time: Duration,
}

impl LookupProfile {
    fn total_known(&self) -> Duration {
        self.extract_cols
            + self.match_loop
            + self.fetch_columns
            + self.vcf_take
            + self.cache_build
            + self.record_batch_build
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
            "[vep-kv-profile] extract_cols={:.3}s ({:.1}%) match_loop={:.3}s ({:.1}%) fetch_columns={:.3}s ({:.1}%) vcf_take={:.3}s ({:.1}%) cache_build={:.3}s ({:.1}%) record_batch={:.3}s ({:.1}%)",
            self.extract_cols.as_secs_f64(),
            Self::pct(self.extract_cols, total),
            self.match_loop.as_secs_f64(),
            Self::pct(self.match_loop, total),
            self.fetch_columns.as_secs_f64(),
            Self::pct(self.fetch_columns, total),
            self.vcf_take.as_secs_f64(),
            Self::pct(self.vcf_take, total),
            self.cache_build.as_secs_f64(),
            Self::pct(self.cache_build, total),
            self.record_batch_build.as_secs_f64(),
            Self::pct(self.record_batch_build, total),
        );
        eprintln!(
            "[vep-kv-profile] window_loads={} window_load_time_s={:.3}",
            self.window_loads,
            self.window_load_time.as_secs_f64()
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
        window_size: u64,
        output_col_indices: Vec<u8>,
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
            window_size,
            output_col_indices,
            current_window: None,
            current_index: None,
            profile_enabled: std::env::var_os("VEP_KV_PROFILE").is_some(),
            profile_emitted: false,
            profile: LookupProfile::default(),
        }
    }

    /// Phase A: Load only the position index (v1).
    fn ensure_window(&mut self, chrom: &str, pos: i32) -> Result<()> {
        let wid = window_id_for_position(i64::from(pos), self.window_size);
        let need_load = match &self.current_window {
            Some((c, w)) => c != chrom || *w != wid,
            None => true,
        };

        if need_load {
            let started = if self.profile_enabled {
                Some(Instant::now())
            } else {
                None
            };
            let pos_index = self.store.get_position_index(chrom, wid)?;
            if let Some(t0) = started {
                self.profile.window_loads += 1;
                self.profile.window_load_time += t0.elapsed();
            }
            self.current_window = Some((chrom.to_string(), wid));
            self.current_index = pos_index.map(WindowAlleleIndex::from_position_index);
        }
        Ok(())
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
        let num_rows = vcf_batch.num_rows();
        if self.profile_enabled {
            self.profile.batches += 1;
            self.profile.input_rows += num_rows as u64;
        }

        let mut matched_entries: Vec<MatchInfoStatic> = Vec::new();
        let mut matched_cache_rows: Vec<usize> = Vec::new();
        let mut unmatched_vcf_rows: Vec<usize> = Vec::new();

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

            self.ensure_window(chrom, norm_start)?;

            let row_start = matched_cache_rows.len();
            if let Some(index) = &self.current_index {
                match self.match_mode {
                    KvMatchMode::Exact => {
                        index.append_matches(
                            norm_start_i64,
                            norm_end_i64,
                            vcf_ref,
                            vcf_alt,
                            self.exact_matcher,
                            &mut matched_cache_rows,
                        );
                        if matched_cache_rows.len() == row_start && norm_start == norm_end {
                            index.append_matches(
                                i64::from(norm_start.saturating_add(1)),
                                norm_start_i64,
                                vcf_ref,
                                vcf_alt,
                                self.exact_matcher,
                                &mut matched_cache_rows,
                            );
                        }
                    }
                    KvMatchMode::ExactOrColocated => {
                        let mut matches = index.find_matches(
                            norm_start_i64,
                            norm_end_i64,
                            vcf_ref,
                            vcf_alt,
                            self.exact_matcher,
                        );
                        if matches.is_empty() && norm_start == norm_end {
                            matches = index.find_matches(
                                i64::from(norm_start.saturating_add(1)),
                                norm_start_i64,
                                vcf_ref,
                                vcf_alt,
                                self.exact_matcher,
                            );
                        }
                        if matches.is_empty() {
                            let mut coloc = index.find_colocated(norm_start_i64, norm_end_i64);
                            if coloc.is_empty() && norm_start == norm_end {
                                coloc = index.find_colocated(
                                    i64::from(norm_start.saturating_add(1)),
                                    norm_start_i64,
                                );
                            }
                            matches = coloc;
                        }
                        matched_cache_rows.extend(matches);
                    }
                    KvMatchMode::ExactOrRelaxed => {
                        let mut matches = index.find_matches(
                            norm_start_i64,
                            norm_end_i64,
                            vcf_ref,
                            vcf_alt,
                            self.exact_matcher,
                        );
                        if matches.is_empty() && norm_start == norm_end {
                            matches = index.find_matches(
                                i64::from(norm_start.saturating_add(1)),
                                norm_start_i64,
                                vcf_ref,
                                vcf_alt,
                                self.exact_matcher,
                            );
                        }
                        if matches.is_empty() {
                            if let Some(relaxed) = self.relaxed_matcher {
                                let mut rel = index.find_matches(
                                    norm_start_i64,
                                    norm_end_i64,
                                    vcf_ref,
                                    vcf_alt,
                                    relaxed,
                                );
                                if rel.is_empty() && norm_start == norm_end {
                                    rel = index.find_matches(
                                        i64::from(norm_start.saturating_add(1)),
                                        norm_start_i64,
                                        vcf_ref,
                                        vcf_alt,
                                        relaxed,
                                    );
                                }
                                matches = rel;
                            }
                        }
                        matched_cache_rows.extend(matches);
                    }
                }
            }

            let row_end = matched_cache_rows.len();
            if row_end > row_start {
                let (c, wid) = self.current_window.as_ref().unwrap();
                matched_entries.push(MatchInfoStatic {
                    vcf_row: row,
                    chrom: c.clone(),
                    window_id: *wid,
                    cache_rows_start: row_start,
                    cache_rows_end: row_end,
                });
            } else {
                unmatched_vcf_rows.push(row);
            }
        }
        if let Some(t0) = match_started {
            self.profile.match_loop += t0.elapsed();
        }

        // Phase B: fetch needed columns for matched windows.
        let total_output_rows: usize = matched_entries
            .iter()
            .map(|m| m.cache_rows_end - m.cache_rows_start)
            .sum::<usize>()
            + unmatched_vcf_rows.len();
        if self.profile_enabled {
            self.profile.output_rows += total_output_rows as u64;
        }

        // Build VCF output indices (expanded for matches).
        let mut vcf_output_indices: Vec<u32> = Vec::with_capacity(total_output_rows);

        // We'll interleave matched and unmatched in original VCF order.
        // Build a combined list sorted by original vcf_row.
        let mut combined: Vec<OutputEntry> =
            Vec::with_capacity(matched_entries.len() + unmatched_vcf_rows.len());
        let mut mi = 0;
        let mut ui = 0;
        while mi < matched_entries.len() || ui < unmatched_vcf_rows.len() {
            let m_row = if mi < matched_entries.len() {
                matched_entries[mi].vcf_row
            } else {
                usize::MAX
            };
            let u_row = if ui < unmatched_vcf_rows.len() {
                unmatched_vcf_rows[ui]
            } else {
                usize::MAX
            };
            if m_row <= u_row {
                combined.push(OutputEntry::Matched(&matched_entries[mi]));
                mi += 1;
            } else {
                combined.push(OutputEntry::Unmatched(unmatched_vcf_rows[ui]));
                ui += 1;
            }
        }

        // Fetch columns for each matched window (batch-fetch per window).
        // Cache: (chrom, window_id) -> fetched column arrays
        let mut window_columns: HashMap<(String, u64), Vec<ArrayRef>> = HashMap::new();

        let fetch_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        if !self.output_col_indices.is_empty() {
            // Collect unique windows that have matches.
            let mut needed_windows: std::collections::HashSet<(String, u64)> =
                std::collections::HashSet::new();
            for entry in &matched_entries {
                needed_windows.insert((entry.chrom.clone(), entry.window_id));
            }
            for (chrom, wid) in needed_windows {
                if let Some(cols) = self
                    .store
                    .get_columns(&chrom, wid, &self.output_col_indices)?
                {
                    window_columns.insert((chrom, wid), cols.into_iter().map(|(a, _)| a).collect());
                }
            }
        }
        if let Some(t0) = fetch_started {
            self.profile.fetch_columns += t0.elapsed();
        }

        // Build expanded VCF indices.
        for entry in &combined {
            match entry {
                OutputEntry::Matched(m) => {
                    for _ in m.cache_rows_start..m.cache_rows_end {
                        vcf_output_indices.push(m.vcf_row as u32);
                    }
                }
                OutputEntry::Unmatched(row) => {
                    vcf_output_indices.push(*row as u32);
                }
            }
        }

        // Take VCF columns.
        let mut output_columns: Vec<ArrayRef> =
            Vec::with_capacity(num_vcf_cols + self.cache_columns.len());
        let take_indices = datafusion::arrow::array::UInt32Array::from(vcf_output_indices);
        let take_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        for col_idx in 0..num_vcf_cols {
            let taken =
                datafusion::arrow::compute::take(vcf_batch.column(col_idx), &take_indices, None)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            output_columns.push(taken);
        }
        if let Some(t0) = take_started {
            self.profile.vcf_take += t0.elapsed();
        }

        // Build cache columns.
        let cache_build_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        for (out_idx, cache_col_name) in self.cache_columns.iter().enumerate() {
            let output_schema_idx = num_vcf_cols + out_idx;
            let data_type = self.schema.field(output_schema_idx).data_type().clone();
            let output_col = build_cache_column_v1(
                &combined,
                &matched_cache_rows,
                &window_columns,
                out_idx,
                &data_type,
                total_output_rows,
                cache_col_name,
            )?;
            output_columns.push(output_col);
        }
        if let Some(t0) = cache_build_started {
            self.profile.cache_build += t0.elapsed();
        }

        let batch_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let out = RecordBatch::try_new(self.schema.clone(), output_columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None));
        if let Some(t0) = batch_started {
            self.profile.record_batch_build += t0.elapsed();
        }
        out
    }
}

enum OutputEntry<'a> {
    Matched(&'a MatchInfoStatic),
    Unmatched(usize),
}

struct MatchInfoStatic {
    vcf_row: usize,
    chrom: String,
    window_id: u64,
    cache_rows_start: usize,
    cache_rows_end: usize,
}

/// Build a cache output column using v1 fetched column arrays.
fn build_cache_column_v1(
    entries: &[OutputEntry<'_>],
    matched_cache_rows: &[usize],
    window_columns: &HashMap<(String, u64), Vec<ArrayRef>>,
    out_col_idx: usize,
    data_type: &DataType,
    total_rows: usize,
    cache_col_name: &str,
) -> Result<ArrayRef> {
    match data_type {
        DataType::Utf8 => {
            let mut builder = StringBuilder::with_capacity(total_rows, total_rows * 16);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_string_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int64 => {
            let mut builder = Int64Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_int64_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int32 => {
            let mut builder = Int32Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_int32_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int16 => {
            let mut builder = Int16Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_int16_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int8 => {
            let mut builder = Int8Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_int8_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::UInt64 => {
            let mut builder = UInt64Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_u64_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::UInt32 => {
            let mut builder = UInt32Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_u32_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::UInt16 => {
            let mut builder = UInt16Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_u16_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::UInt8 => {
            let mut builder = UInt8Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_u8_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Float64 => {
            let mut builder = Float64Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_f64_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Float32 => {
            let mut builder = Float32Builder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_f32_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Boolean => {
            let mut builder = BooleanBuilder::with_capacity(total_rows);
            visit_output_rows(
                entries,
                matched_cache_rows,
                window_columns,
                out_col_idx,
                |value| {
                    if let Some((col, i)) = value {
                        if col.is_null(i) {
                            builder.append_null();
                        } else {
                            builder.append_value(get_bool_value(col, i, cache_col_name)?);
                        }
                    } else {
                        builder.append_null();
                    }
                    Ok(())
                },
            )?;
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        _ => Err(DataFusionError::Execution(format!(
            "cache column '{cache_col_name}' has unsupported output type {data_type:?}"
        ))),
    }
}

fn normalize_cache_output_type(data_type: &DataType) -> DataType {
    match data_type {
        DataType::Utf8View | DataType::LargeUtf8 => DataType::Utf8,
        other => other.clone(),
    }
}

fn visit_output_rows(
    entries: &[OutputEntry<'_>],
    matched_cache_rows: &[usize],
    window_columns: &HashMap<(String, u64), Vec<ArrayRef>>,
    out_col_idx: usize,
    mut f: impl FnMut(Option<(&ArrayRef, usize)>) -> Result<()>,
) -> Result<()> {
    for entry in entries {
        match entry {
            OutputEntry::Matched(m) => {
                let key = (m.chrom.clone(), m.window_id);
                let col = window_columns
                    .get(&key)
                    .and_then(|cols| cols.get(out_col_idx));
                for pos in m.cache_rows_start..m.cache_rows_end {
                    let i = matched_cache_rows[pos];
                    f(col.map(|c| (c, i)))?;
                }
            }
            OutputEntry::Unmatched(_) => f(None)?,
        }
    }
    Ok(())
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

fn get_int64_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<i64> {
    col.as_any()
        .downcast_ref::<Int64Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Int64 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_int32_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<i32> {
    col.as_any()
        .downcast_ref::<Int32Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Int32 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_int16_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<i16> {
    col.as_any()
        .downcast_ref::<Int16Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Int16 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_int8_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<i8> {
    col.as_any()
        .downcast_ref::<Int8Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Int8 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_u64_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<u64> {
    col.as_any()
        .downcast_ref::<UInt64Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected UInt64 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_u32_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<u32> {
    col.as_any()
        .downcast_ref::<UInt32Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected UInt32 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_u16_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<u16> {
    col.as_any()
        .downcast_ref::<UInt16Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected UInt16 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_u8_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<u8> {
    col.as_any()
        .downcast_ref::<UInt8Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected UInt8 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_f64_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<f64> {
    col.as_any()
        .downcast_ref::<Float64Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Float64 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_f32_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<f32> {
    col.as_any()
        .downcast_ref::<Float32Array>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Float32 array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_bool_value(col: &ArrayRef, i: usize, column_name: &str) -> Result<bool> {
    col.as_any()
        .downcast_ref::<BooleanArray>()
        .map(|a| a.value(i))
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "column '{column_name}' expected Boolean array, got {:?}",
                col.data_type()
            ))
        })
}

fn get_string_value<'a>(col: &'a ArrayRef, i: usize, column_name: &str) -> Result<&'a str> {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        Ok(arr.value(i))
    } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        Ok(arr.value(i))
    } else if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
        Ok(arr.value(i))
    } else {
        Err(DataFusionError::Execution(format!(
            "column '{column_name}' expected string array, got {:?}",
            col.data_type()
        )))
    }
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
    use datafusion::arrow::array::{Int32Array, StringArray};
    use datafusion::arrow::datatypes::Schema;
    use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
    use futures::stream;

    fn passthrough_matcher(_: &str, _: &str, _: &str) -> bool {
        true
    }

    fn empty_stream(schema: SchemaRef) -> SendableRecordBatchStream {
        let stream = stream::iter(Vec::<Result<RecordBatch>>::new());
        Box::pin(RecordBatchStreamAdapter::new(schema, Box::pin(stream)))
    }

    fn minimal_cache_schema() -> SchemaRef {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ]))
    }

    #[test]
    fn test_build_cache_column_v1_type_mismatch_returns_error() {
        let matched = MatchInfoStatic {
            vcf_row: 0,
            chrom: "1".to_string(),
            window_id: 0,
            cache_rows_start: 0,
            cache_rows_end: 1,
        };
        let entries = vec![OutputEntry::Matched(&matched)];
        let matched_cache_rows = vec![0usize];
        let mut window_columns: HashMap<(String, u64), Vec<ArrayRef>> = HashMap::new();
        window_columns.insert(
            ("1".to_string(), 0),
            vec![Arc::new(StringArray::from(vec!["not_an_int64"])) as ArrayRef],
        );

        let err = build_cache_column_v1(
            &entries,
            &matched_cache_rows,
            &window_columns,
            0,
            &DataType::Int64,
            1,
            "cache_test_col",
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("cache_test_col"));
        assert!(err.contains("Int64"));
    }

    #[test]
    fn test_process_batch_accepts_int32_start_column() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            Arc::new(VepKvStore::create(dir.path(), minimal_cache_schema(), 1_000_000).unwrap());

        let vcf_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int32, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int32Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();

        let mut stream = KvLookupStream::new(
            empty_stream(vcf_schema.clone()),
            store,
            vcf_schema,
            vec![],
            KvMatchMode::Exact,
            passthrough_matcher,
            None,
            false,
            false,
            false,
            1_000_000,
            vec![],
        );

        let output = stream.process_batch(&vcf_batch).unwrap();
        assert_eq!(output.num_rows(), 1);
        assert_eq!(output.num_columns(), 5);
    }
}
