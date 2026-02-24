//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store window-by-window for annotation.

use std::any::Any;
use std::collections::{HashMap, VecDeque};
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

use crate::allele_index::{AlleleMatcher, WindowAlleleIndex};
use crate::key_encoding::{chrom_to_code, encode_position_key_buf, window_id_for_position};
use crate::kv_store::{
    FORMAT_V1, FORMAT_V2, FORMAT_V3, FORMAT_V4, FORMAT_V5, VepKvStore, to_v1_column_index,
    validate_v1_schema_width,
};
use crate::position_entry::{PositionEntryReader, make_builder};

const DEFAULT_WINDOW_COLUMNS_CACHE_WINDOWS: usize = 8;

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
    cache_format_version: u8,
    /// Cache schema column positions for requested cache output columns.
    output_col_positions: Vec<usize>,
    /// For v1 format only: encoded cache column indices used by `get_columns`.
    output_col_indices_v1: Vec<u8>,
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
        let cache_format_version = store.format_version();
        if cache_format_version != FORMAT_V1
            && cache_format_version != FORMAT_V2
            && cache_format_version != FORMAT_V3
            && cache_format_version != FORMAT_V4
            && cache_format_version != FORMAT_V5
        {
            return Err(DataFusionError::Execution(format!(
                "unsupported cache format version {cache_format_version}: supported versions are {FORMAT_V1}, {FORMAT_V2}, {FORMAT_V3}, {FORMAT_V4}, and {FORMAT_V5}"
            )));
        }
        if cache_format_version == FORMAT_V1 {
            validate_v1_schema_width(cache_schema.fields().len())?;
        }

        let mut output_col_positions = Vec::new();
        let mut output_col_indices_v1 = Vec::new();
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
                    if cache_format_version == FORMAT_V1 {
                        output_col_indices_v1.push(to_v1_column_index(idx)?);
                    }
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
            cache_format_version,
            output_col_positions,
            output_col_indices_v1,
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
            self.cache_format_version,
            self.output_col_positions.clone(),
            self.output_col_indices_v1.clone(),
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
    cache_format_version: u8,
    output_col_positions: Vec<usize>,
    output_col_indices_v1: Vec<u8>,
    current_window: Option<(u16, u64)>,
    current_window_chrom_str: Option<String>,
    current_index: Option<WindowAlleleIndex>,
    window_columns_cache: HashMap<(u16, u64), Arc<Vec<ArrayRef>>>,
    window_columns_cache_lru: VecDeque<(u16, u64)>,
    window_columns_cache_capacity: usize,
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
    fetch_ref_time: Duration,
    fetch_deser_time: Duration,
    fetch_total_compressed: u64,
    fetch_total_uncompressed: u64,
    vcf_take: Duration,
    cache_build: Duration,
    record_batch_build: Duration,
    window_loads: u64,
    window_load_time: Duration,
    window_col_cache_hits: u64,
    window_col_cache_misses: u64,
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
        eprintln!(
            "[vep-kv-profile] window_col_cache hits={} misses={}",
            self.window_col_cache_hits, self.window_col_cache_misses
        );
        if self.fetch_deser_time.as_nanos() > 0 || self.fetch_ref_time.as_nanos() > 0 {
            eprintln!(
                "[vep-kv-profile] fetch_ref={:.3}s fetch_deser={:.3}s (comp={:.1}MB uncomp={:.1}MB)",
                self.fetch_ref_time.as_secs_f64(),
                self.fetch_deser_time.as_secs_f64(),
                self.fetch_total_compressed as f64 / 1_048_576.0,
                self.fetch_total_uncompressed as f64 / 1_048_576.0,
            );
        }
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
        cache_format_version: u8,
        output_col_positions: Vec<usize>,
        output_col_indices_v1: Vec<u8>,
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
            cache_format_version,
            output_col_positions,
            output_col_indices_v1,
            current_window: None,
            current_window_chrom_str: None,
            current_index: None,
            window_columns_cache: HashMap::new(),
            window_columns_cache_lru: VecDeque::new(),
            window_columns_cache_capacity: window_columns_cache_capacity(),
            profile_enabled: std::env::var_os("VEP_KV_PROFILE").is_some(),
            profile_emitted: false,
            profile: LookupProfile::default(),
        }
    }

    fn get_cached_window_columns(&mut self, key: &(u16, u64)) -> Option<Arc<Vec<ArrayRef>>> {
        if self.window_columns_cache_capacity == 0 {
            return None;
        }
        let cols = self.window_columns_cache.get(key).cloned();
        if cols.is_some() {
            self.touch_window_cache_key(key);
        }
        cols
    }

    fn put_cached_window_columns(&mut self, key: (u16, u64), columns: Arc<Vec<ArrayRef>>) {
        if self.window_columns_cache_capacity == 0 {
            return;
        }

        let exists = self.window_columns_cache.contains_key(&key);
        if !exists {
            while self.window_columns_cache.len() >= self.window_columns_cache_capacity {
                let Some(oldest) = self.window_columns_cache_lru.pop_front() else {
                    break;
                };
                self.window_columns_cache.remove(&oldest);
            }
        }

        self.window_columns_cache.insert(key, columns);
        self.touch_window_cache_key(&key);
    }

    fn touch_window_cache_key(&mut self, key: &(u16, u64)) {
        if self.window_columns_cache_capacity == 0 {
            return;
        }
        if let Some(pos) = self.window_columns_cache_lru.iter().position(|k| k == key) {
            self.window_columns_cache_lru.remove(pos);
        }
        self.window_columns_cache_lru.push_back(*key);
    }

    /// Phase A: Load only the position index.
    fn ensure_window(&mut self, chrom: &str, chrom_code: u16, pos: i32) -> Result<()> {
        let wid = window_id_for_position(i64::from(pos), self.window_size);
        let need_load = match &self.current_window {
            Some((c, w)) => *c != chrom_code || *w != wid,
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
            self.current_window = Some((chrom_code, wid));
            self.current_window_chrom_str = Some(chrom.to_string());
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

            let chrom_code = chrom_to_code(chrom);
            self.ensure_window(chrom, chrom_code, norm_start)?;

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
                let &(cc, wid) = self.current_window.as_ref().unwrap();
                matched_entries.push(MatchInfoStatic {
                    vcf_row: row,
                    chrom_code: cc,
                    window_id: wid,
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
        // Cache: (chrom_code, window_id) -> fetched column arrays
        let mut window_columns: HashMap<(u16, u64), Arc<Vec<ArrayRef>>> = HashMap::new();
        // Map chrom_code -> chrom string for KV lookups (only populated on cache misses).
        let mut code_to_chrom_str: HashMap<u16, String> = HashMap::new();

        let fetch_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        if !self.output_col_positions.is_empty() {
            // Collect unique windows that have matches.
            let mut needed_windows: std::collections::HashSet<(u16, u64)> =
                std::collections::HashSet::new();
            for entry in &matched_entries {
                needed_windows.insert((entry.chrom_code, entry.window_id));
            }
            // Build code-to-string map. Use key_encoding reverse lookup first,
            // fall back to current_window_chrom_str for the last window.
            for &(cc, _wid) in &needed_windows {
                if let std::collections::hash_map::Entry::Vacant(e) = code_to_chrom_str.entry(cc) {
                    if let Some(name) = crate::key_encoding::code_to_chrom(cc) {
                        e.insert(name.to_string());
                    } else if let Some(ref s) = self.current_window_chrom_str {
                        if let Some(&(cur_cc, _)) = self.current_window.as_ref() {
                            if cur_cc == cc {
                                e.insert(s.clone());
                            }
                        }
                    }
                }
            }
            for (cc, wid) in needed_windows {
                let key = (cc, wid);
                if let Some(cols) = self.get_cached_window_columns(&key) {
                    self.profile.window_col_cache_hits += 1;
                    window_columns.insert(key, cols);
                    continue;
                }
                self.profile.window_col_cache_misses += 1;

                let chrom = code_to_chrom_str.get(&cc).ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "no chrom string for code {cc:#06x} in fetch path"
                    ))
                })?;

                let ref_started = if self.profile_enabled {
                    Some(Instant::now())
                } else {
                    None
                };
                let fetched = match self.cache_format_version {
                    FORMAT_V1 => {
                        if let Some(t0) = ref_started {
                            self.profile.fetch_ref_time += t0.elapsed();
                        }
                        self.store
                            .get_columns(chrom, wid, &self.output_col_indices_v1)?
                            .map(|cols| cols.into_iter().map(|(a, _)| a).collect::<Vec<_>>())
                    }
                    FORMAT_V2 => {
                        if let Some(t0) = ref_started {
                            self.profile.fetch_ref_time += t0.elapsed();
                        }
                        self.store
                            .get_window_columns_v2(chrom, wid, &self.output_col_positions)?
                    }
                    FORMAT_V3 => {
                        if let Some(t0) = ref_started {
                            self.profile.fetch_ref_time += t0.elapsed();
                        }
                        self.store
                            .get_window_columns_v3(chrom, wid, &self.output_col_positions)?
                    }
                    FORMAT_V4 => {
                        let (block_ref, ref_dur) =
                            self.store.get_window_ref_v4_timed(chrom, wid)?;
                        if self.profile_enabled {
                            self.profile.fetch_ref_time += ref_dur;
                        }
                        match block_ref {
                            Some(br) => {
                                let deser_start = Instant::now();
                                let cols = self
                                    .store
                                    .read_columns_v4_direct(br, &self.output_col_positions)?;
                                if self.profile_enabled {
                                    self.profile.fetch_deser_time += deser_start.elapsed();
                                }
                                Some(cols)
                            }
                            None => None,
                        }
                    }
                    other => {
                        return Err(DataFusionError::Execution(format!(
                            "unsupported cache format version in fetch path: {other}"
                        )));
                    }
                };

                if let Some(cols) = fetched {
                    let cols = Arc::new(cols);
                    self.put_cached_window_columns(key, cols.clone());
                    window_columns.insert(key, cols);
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

        // Build cache columns using vectorized take-based assembly.
        let cache_build_started = if self.profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let cache_cols = build_cache_columns_take(
            &combined,
            &matched_cache_rows,
            &window_columns,
            self.cache_columns.len(),
            total_output_rows,
            &self.schema,
            num_vcf_cols,
        )?;
        output_columns.extend(cache_cols);
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

    /// V5 read path: single-pass position-keyed lookup.
    ///
    /// For each VCF row, fetch the per-position entry from fjall, match alleles,
    /// and append matched column values directly into ArrayBuilders.
    fn process_batch_v5(&mut self, vcf_batch: &RecordBatch) -> Result<RecordBatch> {
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

        // Reusable key buffer.
        let mut key_buf: Vec<u8> = Vec::with_capacity(18);

        // Reusable zstd decompressor — created once, amortized across all lookups.
        let mut v5_decompressor = self.store.create_v5_decompressor()?;

        // Reusable decompression / raw-value buffer — avoids alloc per lookup.
        let mut decompress_buf: Vec<u8> = Vec::with_capacity(4096);

        // Reusable allele match buffer — avoids alloc per row.
        let mut matched_allele_rows: Vec<usize> = Vec::new();

        // Determine which column indices in the V5 entry correspond to our output columns.
        // V5 entry stores all columns except chrom/start/end, in schema order minus those 3.
        let cache_schema = self.store.schema();
        let cache_chrom_idx = cache_schema.index_of("chrom").unwrap_or(usize::MAX);
        let cache_start_idx = cache_schema.index_of("start").unwrap_or(usize::MAX);
        let cache_end_idx = cache_schema.index_of("end").unwrap_or(usize::MAX);

        // Build mapping: output_col_positions[i] -> index within the V5 entry's column list.
        let v5_stored_cols: Vec<usize> = (0..cache_schema.fields().len())
            .filter(|&i| i != cache_chrom_idx && i != cache_start_idx && i != cache_end_idx)
            .collect();
        let v5_col_map: Vec<usize> = self
            .output_col_positions
            .iter()
            .map(|&pos| {
                v5_stored_cols
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

            // Try fetching the position entry (with zstd decompression if dict exists).
            encode_position_key_buf(chrom_code, norm_start_i64, norm_end_i64, &mut key_buf);
            let mut found = self.store.get_position_entry_fast(
                chrom_code,
                norm_start_i64,
                norm_end_i64,
                v5_decompressor.as_mut(),
                &mut decompress_buf,
            )?;

            // SNP retry: if no match and start==end, try (start+1, start).
            if !found && norm_start == norm_end {
                let shifted = i64::from(norm_start.saturating_add(1));
                found = self.store.get_position_entry_fast(
                    chrom_code,
                    shifted,
                    norm_start_i64,
                    v5_decompressor.as_mut(),
                    &mut decompress_buf,
                )?;
            }

            if found {
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
                    // Position exists but no allele match -> null cache columns.
                    vcf_indices.push(row as u32);
                    for builder in &mut builders {
                        append_null_to_builder(builder.as_mut())?;
                    }
                } else {
                    // Append matched rows to builders.
                    for _ in &matched_allele_rows {
                        vcf_indices.push(row as u32);
                    }
                    for (col_out_idx, builder) in builders.iter_mut().enumerate() {
                        let v5_idx = v5_col_map[col_out_idx];
                        if v5_idx == usize::MAX {
                            // Column not found in V5 entry -> nulls.
                            for _ in &matched_allele_rows {
                                append_null_to_builder(builder.as_mut())?;
                            }
                        } else {
                            reader.append_column_values(
                                v5_idx,
                                &matched_allele_rows,
                                builder.as_mut(),
                            )?;
                        }
                    }
                }
            } else {
                // Position not found -> null cache columns.
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

enum OutputEntry<'a> {
    Matched(&'a MatchInfoStatic),
    Unmatched(usize),
}

struct MatchInfoStatic {
    vcf_row: usize,
    chrom_code: u16,
    window_id: u64,
    cache_rows_start: usize,
    cache_rows_end: usize,
}

/// Build all cache output columns using vectorized `arrow::compute::take` + `concat`.
///
/// Algorithm:
/// 1. Walk output entries, grouping matched rows by (chrom_code, window_id).
/// 2. Per window: build a `UInt32Array` of cache row indices and a parallel
///    `UInt32Array` of output position indices.
/// 3. Per output column: `take(window_col, &cache_indices)` then scatter
///    into the correct output positions; unmatched positions get null.
/// 4. Uses `arrow::compute::interleave` to combine per-window slices
///    and null slices into the final output column.
fn build_cache_columns_take(
    entries: &[OutputEntry<'_>],
    matched_cache_rows: &[usize],
    window_columns: &HashMap<(u16, u64), Arc<Vec<ArrayRef>>>,
    num_cache_cols: usize,
    total_output_rows: usize,
    output_schema: &SchemaRef,
    num_vcf_cols: usize,
) -> Result<Vec<ArrayRef>> {
    if num_cache_cols == 0 {
        return Ok(vec![]);
    }

    // Build a list of (window_key, cache_row_indices, output_positions) per window,
    // plus track null (unmatched) output row ranges.
    // We'll build interleave instructions: a vec of (array_index, row_within_array)
    // per cache output column.

    // Step 1: collect per-window take indices and output position mapping.
    // We process entries in output order and build "segments" — contiguous runs
    // from the same window or null runs — that we'll later interleave.
    struct Segment {
        /// Index into `window_arrays` vec, or usize::MAX for null segment.
        source_idx: usize,
        /// Row indices within the source array for this segment.
        rows: Vec<u32>,
    }

    // Map from window key to index in window_arrays.
    let mut window_key_to_idx: HashMap<(u16, u64), usize> = HashMap::new();
    // Collected per-window arrays (one per unique window).
    let mut window_array_keys: Vec<(u16, u64)> = Vec::new();

    let mut segments: Vec<Segment> = Vec::new();

    for entry in entries {
        match entry {
            OutputEntry::Matched(m) => {
                let key = (m.chrom_code, m.window_id);
                let source_idx = match window_key_to_idx.get(&key) {
                    Some(&idx) => idx,
                    None => {
                        let idx = window_array_keys.len();
                        window_key_to_idx.insert(key, idx);
                        window_array_keys.push(key);
                        idx
                    }
                };
                let rows: Vec<u32> = (m.cache_rows_start..m.cache_rows_end)
                    .map(|pos| matched_cache_rows[pos] as u32)
                    .collect();
                segments.push(Segment { source_idx, rows });
            }
            OutputEntry::Unmatched(_) => {
                segments.push(Segment {
                    source_idx: usize::MAX,
                    rows: vec![0], // single null row placeholder
                });
            }
        }
    }

    // Step 2: per cache column, build taken arrays per window, then interleave.
    let mut result = Vec::with_capacity(num_cache_cols);

    for out_col_idx in 0..num_cache_cols {
        let output_schema_idx = num_vcf_cols + out_col_idx;
        let output_data_type = output_schema.field(output_schema_idx).data_type();

        // Pre-take: for each window, take the needed rows from this column.
        let mut taken_per_window: Vec<ArrayRef> = Vec::with_capacity(window_array_keys.len());
        for key in &window_array_keys {
            let cols = window_columns.get(key);
            let col = cols.and_then(|c| c.get(out_col_idx));
            if let Some(col) = col {
                // Collect ALL cache row indices needed from this window across all segments.
                let mut all_indices: Vec<u32> = Vec::new();
                for seg in &segments {
                    if seg.source_idx == window_key_to_idx[key] {
                        all_indices.extend_from_slice(&seg.rows);
                    }
                }
                let indices = UInt32Array::from(all_indices);
                let taken = datafusion::arrow::compute::take(col.as_ref(), &indices, None)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                // Cast to output type if needed (e.g. StringView -> Utf8).
                let taken = cast_to_output_type(&taken, output_data_type)?;
                taken_per_window.push(taken);
            } else {
                // Window had no data for this column — will produce nulls.
                let mut null_count = 0;
                for seg in &segments {
                    if seg.source_idx == window_key_to_idx[key] {
                        null_count += seg.rows.len();
                    }
                }
                taken_per_window.push(datafusion::arrow::array::new_null_array(
                    output_data_type,
                    null_count,
                ));
            }
        }

        // Build a single null array for unmatched rows.
        let null_count: usize = segments
            .iter()
            .filter(|s| s.source_idx == usize::MAX)
            .map(|s| s.rows.len())
            .sum();
        let null_array = datafusion::arrow::array::new_null_array(output_data_type, null_count);

        // Now build interleave spec: (array_index, row_index) pairs.
        // array layout: [taken_per_window[0], taken_per_window[1], ..., null_array]
        let null_array_idx = taken_per_window.len();
        let all_arrays: Vec<&dyn Array> = taken_per_window
            .iter()
            .map(|a| a.as_ref())
            .chain(std::iter::once(null_array.as_ref()))
            .collect();

        // Track how many rows we've consumed from each taken array and from the null array.
        let mut consumed_per_window: Vec<usize> = vec![0; taken_per_window.len()];
        let mut null_consumed: usize = 0;

        let mut interleave_indices: Vec<(usize, usize)> = Vec::with_capacity(total_output_rows);
        for seg in &segments {
            if seg.source_idx == usize::MAX {
                // Null segment.
                for _ in &seg.rows {
                    interleave_indices.push((null_array_idx, null_consumed));
                    null_consumed += 1;
                }
            } else {
                let arr_idx = seg.source_idx;
                for _ in &seg.rows {
                    interleave_indices.push((arr_idx, consumed_per_window[arr_idx]));
                    consumed_per_window[arr_idx] += 1;
                }
            }
        }

        let output_col = datafusion::arrow::compute::interleave(&all_arrays, &interleave_indices)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;

        result.push(output_col);
    }

    Ok(result)
}

/// Cast an array to the output type if needed (e.g., StringView/LargeUtf8 -> Utf8).
fn cast_to_output_type(array: &ArrayRef, target: &DataType) -> Result<ArrayRef> {
    if array.data_type() == target {
        return Ok(array.clone());
    }
    datafusion::arrow::compute::cast(array, target)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
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
            "V5: unsupported builder type for null append".into(),
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

fn window_columns_cache_capacity() -> usize {
    std::env::var("VEP_KV_WINDOW_CACHE_WINDOWS")
        .ok()
        .and_then(|raw| raw.parse::<usize>().ok())
        .unwrap_or(DEFAULT_WINDOW_COLUMNS_CACHE_WINDOWS)
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

impl Stream for KvLookupStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        match self.input.poll_next_unpin(cx) {
            Poll::Ready(Some(Ok(batch))) => {
                let result = if self.cache_format_version == FORMAT_V5 {
                    self.process_batch_v5(&batch)
                } else {
                    self.process_batch(&batch)
                };
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
    use crate::PositionIndex;
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
    fn test_build_cache_columns_take_produces_correct_output() {
        let chrom_code = chrom_to_code("1");
        let matched = MatchInfoStatic {
            vcf_row: 0,
            chrom_code,
            window_id: 0,
            cache_rows_start: 0,
            cache_rows_end: 1,
        };
        let entries = vec![OutputEntry::Matched(&matched)];
        let matched_cache_rows = vec![0usize];
        let mut window_columns: HashMap<(u16, u64), Arc<Vec<ArrayRef>>> = HashMap::new();
        window_columns.insert(
            (chrom_code, 0),
            Arc::new(vec![Arc::new(StringArray::from(vec!["hello"])) as ArrayRef]),
        );

        // Output schema: 0 VCF cols + 1 cache col (Utf8)
        let output_schema = Arc::new(Schema::new(vec![Field::new(
            "cache_test_col",
            DataType::Utf8,
            true,
        )]));

        let cols = build_cache_columns_take(
            &entries,
            &matched_cache_rows,
            &window_columns,
            1,
            1,
            &output_schema,
            0,
        )
        .unwrap();
        assert_eq!(cols.len(), 1);
        let arr = cols[0].as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(arr.value(0), "hello");
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
            FORMAT_V1,
            vec![],
            vec![],
        );

        let output = stream.process_batch(&vcf_batch).unwrap();
        assert_eq!(output.num_rows(), 1);
        assert_eq!(output.num_columns(), 5);
    }

    #[test]
    fn test_process_batch_reads_cache_columns_from_v2_window_blocks() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            Arc::new(VepKvStore::create_v2(dir.path(), minimal_cache_schema(), 1_000_000).unwrap());

        let cache_batch = RecordBatch::try_new(
            minimal_cache_schema(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec![Some("rs123")])),
                Arc::new(StringArray::from(vec!["A/G"])),
            ],
        )
        .unwrap();
        let pos_index = PositionIndex::from_batch(&cache_batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        store.put_window_block_v2("1", 0, &cache_batch).unwrap();

        let vcf_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int32, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let output_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int32, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
            Field::new("cache_variation_name", DataType::Utf8, true),
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
            output_schema,
            vec!["variation_name".to_string()],
            KvMatchMode::Exact,
            passthrough_matcher,
            None,
            false,
            false,
            false,
            1_000_000,
            FORMAT_V2,
            vec![3],
            vec![],
        );

        let output = stream.process_batch(&vcf_batch).unwrap();
        assert_eq!(output.num_rows(), 1);
        assert_eq!(output.num_columns(), 6);

        let col = output
            .column(5)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(col.value(0), "rs123");
    }

    #[test]
    fn test_process_batch_reads_cache_columns_from_v3_window_blocks() {
        let dir = tempfile::tempdir().unwrap();
        let store =
            Arc::new(VepKvStore::create_v3(dir.path(), minimal_cache_schema(), 1_000_000).unwrap());

        let cache_batch = RecordBatch::try_new(
            minimal_cache_schema(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec![Some("rs789")])),
                Arc::new(StringArray::from(vec!["C/T"])),
            ],
        )
        .unwrap();
        let pos_index = PositionIndex::from_batch(&cache_batch).unwrap();
        store.put_position_index("1", 0, &pos_index).unwrap();
        store.put_window_block_v3("1", 0, &cache_batch).unwrap();

        let vcf_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int32, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let output_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int32, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
            Field::new("cache_variation_name", DataType::Utf8, true),
        ]));

        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int32Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["C"])),
                Arc::new(StringArray::from(vec!["T"])),
            ],
        )
        .unwrap();

        let mut stream = KvLookupStream::new(
            empty_stream(vcf_schema.clone()),
            store,
            output_schema,
            vec!["variation_name".to_string()],
            KvMatchMode::Exact,
            passthrough_matcher,
            None,
            false,
            false,
            false,
            1_000_000,
            FORMAT_V3,
            vec![3],
            vec![],
        );

        let output = stream.process_batch(&vcf_batch).unwrap();
        assert_eq!(output.num_rows(), 1);
        assert_eq!(output.num_columns(), 6);

        let col = output
            .column(5)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(col.value(0), "rs789");
    }
}
