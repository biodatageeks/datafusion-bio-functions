//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store window-by-window for annotation.

use std::any::Any;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::arrow::array::builder::{Int8Builder, Int64Builder, StringBuilder};
use datafusion::arrow::array::{Array, ArrayRef, Int8Array, Int64Array, NullArray, RecordBatch};
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
        if store.format_version() >= FORMAT_V1 {
            validate_v1_schema_width(cache_schema.fields().len())?;
        }

        // Compute output column indices for v1 format.
        let mut output_col_indices = Vec::new();
        let mut fields: Vec<Arc<Field>> = input_schema.fields().iter().cloned().collect();
        for col_name in &cache_columns {
            if let Ok(field) = cache_schema.field_with_name(col_name) {
                fields.push(Arc::new(Field::new(
                    format!("cache_{}", field.name()),
                    field.data_type().clone(),
                    true,
                )));
                if store.format_version() >= FORMAT_V1 {
                    if let Ok(idx) = cache_schema.index_of(col_name) {
                        output_col_indices.push(to_v1_column_index(idx)?);
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
            self.store.format_version(),
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
    format_version: u8,
    output_col_indices: Vec<u8>,
    current_window: Option<(String, u64)>,
    current_index: Option<WindowAlleleIndex>,
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
        format_version: u8,
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
            format_version,
            output_col_indices,
            current_window: None,
            current_index: None,
        }
    }

    /// Phase A: Load only the position index (v1) or full batch (v0).
    fn ensure_window(&mut self, chrom: &str, pos: i64) -> Result<()> {
        let wid = window_id_for_position(pos, self.window_size);
        let need_load = match &self.current_window {
            Some((c, w)) => c != chrom || *w != wid,
            None => true,
        };

        if need_load {
            if self.format_version >= FORMAT_V1 {
                // v1 fast path: load only position index (~2-8KB)
                let pos_index = self.store.get_position_index(chrom, wid)?;
                self.current_window = Some((chrom.to_string(), wid));
                self.current_index = pos_index.map(WindowAlleleIndex::from_position_index);
            } else {
                // v0 legacy path: load full monolithic batch
                let batch = self.store.get_window(chrom, wid)?;
                self.current_window = Some((chrom.to_string(), wid));
                self.current_index = match batch {
                    Some(b) => Some(WindowAlleleIndex::from_batch(b)?),
                    None => None,
                };
            }
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

        let chroms = get_string_column(vcf_batch.column(chrom_idx), "chrom")?;
        let starts = get_int64_column(vcf_batch.column(start_idx), "start")?;
        let ends = get_int64_column(vcf_batch.column(end_idx), "end")?;
        let refs = get_string_column(vcf_batch.column(ref_idx), "ref")?;
        let alts = get_string_column(vcf_batch.column(alt_idx), "alt")?;

        let cache_schema = self.store.schema().clone();
        let num_vcf_cols = vcf_schema.fields().len();
        let num_rows = vcf_batch.num_rows();

        let mut matched_entries: Vec<MatchInfoStatic> = Vec::new();
        let mut unmatched_vcf_rows: Vec<usize> = Vec::new();

        for row in 0..num_rows {
            let raw_chrom = chroms[row].as_deref().unwrap_or("");
            let chrom = if self.vcf_has_chr {
                raw_chrom.strip_prefix("chr").unwrap_or(raw_chrom)
            } else {
                raw_chrom
            };

            let vcf_start = starts.value(row);
            let vcf_end = ends.value(row);
            let (norm_start, norm_end) = normalize_vcf_coords(
                vcf_start,
                vcf_end,
                self.vcf_zero_based,
                self.cache_zero_based,
            );

            let vcf_ref = refs[row].as_deref().unwrap_or("");
            let vcf_alt = alts[row].as_deref().unwrap_or("");

            self.ensure_window(chrom, norm_start)?;

            let matched = if let Some(index) = &self.current_index {
                let mut matches =
                    index.find_matches(norm_start, norm_end, vcf_ref, vcf_alt, self.exact_matcher);

                // For insertions, also check cache coords where start=end+1.
                if matches.is_empty() && norm_start == norm_end {
                    matches = index.find_matches(
                        norm_start + 1,
                        norm_start,
                        vcf_ref,
                        vcf_alt,
                        self.exact_matcher,
                    );
                }

                if matches.is_empty() {
                    match self.match_mode {
                        KvMatchMode::Exact => {}
                        KvMatchMode::ExactOrColocated => {
                            let mut coloc = index.find_colocated(norm_start, norm_end);
                            if coloc.is_empty() && norm_start == norm_end {
                                coloc = index.find_colocated(norm_start + 1, norm_start);
                            }
                            matches = coloc;
                        }
                        KvMatchMode::ExactOrRelaxed => {
                            if let Some(relaxed) = self.relaxed_matcher {
                                let mut rel = index
                                    .find_matches(norm_start, norm_end, vcf_ref, vcf_alt, relaxed);
                                if rel.is_empty() && norm_start == norm_end {
                                    rel = index.find_matches(
                                        norm_start + 1,
                                        norm_start,
                                        vcf_ref,
                                        vcf_alt,
                                        relaxed,
                                    );
                                }
                                matches = rel;
                            }
                        }
                    }
                }

                if matches.is_empty() {
                    None
                } else {
                    Some(matches)
                }
            } else {
                None
            };

            match matched {
                Some(cache_rows) => {
                    let (c, wid) = self.current_window.as_ref().unwrap();
                    matched_entries.push(MatchInfoStatic {
                        vcf_row: row,
                        chrom: c.clone(),
                        window_id: *wid,
                        cache_rows,
                    });
                }
                None => {
                    unmatched_vcf_rows.push(row);
                }
            }
        }

        // Phase B: fetch needed columns for matched windows.
        let total_output_rows: usize = matched_entries
            .iter()
            .map(|m| m.cache_rows.len())
            .sum::<usize>()
            + unmatched_vcf_rows.len();

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

        if self.format_version >= FORMAT_V1 && !self.output_col_indices.is_empty() {
            // Collect unique windows that have matches.
            let mut needed_windows: Vec<(String, u64)> = Vec::new();
            for entry in &matched_entries {
                let key = (entry.chrom.clone(), entry.window_id);
                if !needed_windows.contains(&key) {
                    needed_windows.push(key);
                }
            }
            for (chrom, wid) in &needed_windows {
                if let Some(cols) = self
                    .store
                    .get_columns(chrom, *wid, &self.output_col_indices)?
                {
                    window_columns.insert(
                        (chrom.clone(), *wid),
                        cols.into_iter().map(|(a, _)| a).collect(),
                    );
                }
            }
        }

        // Build expanded VCF indices.
        for entry in &combined {
            match entry {
                OutputEntry::Matched(m) => {
                    for _ in &m.cache_rows {
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
        for col_idx in 0..num_vcf_cols {
            let taken =
                datafusion::arrow::compute::take(vcf_batch.column(col_idx), &take_indices, None)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            output_columns.push(taken);
        }

        // Build cache columns.
        if self.format_version >= FORMAT_V1 {
            // v1 path: use fetched column arrays directly.
            for (out_idx, cache_col_name) in self.cache_columns.iter().enumerate() {
                let data_type = cache_schema
                    .field_with_name(cache_col_name)
                    .map(|f| f.data_type().clone())
                    .unwrap_or(DataType::Utf8);

                let output_col = build_cache_column_v1(
                    &combined,
                    &window_columns,
                    out_idx,
                    &data_type,
                    total_output_rows,
                    cache_col_name,
                )?;
                output_columns.push(output_col);
            }
        } else {
            // v0 legacy path: we don't have column data in v0 through this path.
            // For v0, we need to fall back to the old approach of reading from the batch.
            // But since WindowAlleleIndex no longer holds a batch, we need to read the
            // full batch from the store for v0 format.
            // However, this code path only triggers for v0 stores, where we need to
            // reconstruct the batch-based approach.
            for cache_col_name in &self.cache_columns {
                let cache_col_idx = cache_schema.index_of(cache_col_name).ok();
                let output_col = build_cache_column_v0_fallback(
                    &combined,
                    &self.store,
                    cache_col_idx,
                    cache_col_name,
                    cache_schema.as_ref(),
                    total_output_rows,
                )?;
                output_columns.push(output_col);
            }
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
    chrom: String,
    window_id: u64,
    cache_rows: Vec<usize>,
}

/// Build a cache output column using v1 fetched column arrays.
fn build_cache_column_v1(
    entries: &[OutputEntry<'_>],
    window_columns: &HashMap<(String, u64), Vec<ArrayRef>>,
    out_col_idx: usize,
    data_type: &DataType,
    total_rows: usize,
    cache_col_name: &str,
) -> Result<ArrayRef> {
    match data_type {
        DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => {
            let mut builder = StringBuilder::with_capacity(total_rows, total_rows * 16);
            for entry in entries {
                match entry {
                    OutputEntry::Matched(m) => {
                        let key = (m.chrom.clone(), m.window_id);
                        let col = window_columns
                            .get(&key)
                            .and_then(|cols| cols.get(out_col_idx));
                        for &i in &m.cache_rows {
                            if let Some(col) = col {
                                if col.is_null(i) {
                                    builder.append_null();
                                } else {
                                    builder.append_value(get_string_value(col, i, cache_col_name)?);
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    OutputEntry::Unmatched(_) => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int64 => {
            let mut builder = Int64Builder::with_capacity(total_rows);
            for entry in entries {
                match entry {
                    OutputEntry::Matched(m) => {
                        let key = (m.chrom.clone(), m.window_id);
                        let col = window_columns
                            .get(&key)
                            .and_then(|cols| cols.get(out_col_idx));
                        for &i in &m.cache_rows {
                            if let Some(col) = col {
                                if col.is_null(i) {
                                    builder.append_null();
                                } else {
                                    builder.append_value(
                                        get_int64_column(col, cache_col_name)?.value(i),
                                    );
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    OutputEntry::Unmatched(_) => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int8 => {
            let mut builder = Int8Builder::with_capacity(total_rows);
            for entry in entries {
                match entry {
                    OutputEntry::Matched(m) => {
                        let key = (m.chrom.clone(), m.window_id);
                        let col = window_columns
                            .get(&key)
                            .and_then(|cols| cols.get(out_col_idx));
                        for &i in &m.cache_rows {
                            if let Some(col) = col {
                                if col.is_null(i) {
                                    builder.append_null();
                                } else {
                                    builder.append_value(
                                        get_int8_column(col, cache_col_name)?.value(i),
                                    );
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    OutputEntry::Unmatched(_) => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        _ => Ok(Arc::new(NullArray::new(total_rows)) as ArrayRef),
    }
}

/// Build a cache output column for v0 format (reads full window batch from store).
fn build_cache_column_v0_fallback(
    entries: &[OutputEntry<'_>],
    store: &VepKvStore,
    cache_col_idx: Option<usize>,
    col_name: &str,
    cache_schema: &datafusion::arrow::datatypes::Schema,
    total_rows: usize,
) -> Result<ArrayRef> {
    let data_type = cache_schema
        .field_with_name(col_name)
        .map(|f| f.data_type().clone())
        .unwrap_or(DataType::Utf8);

    // Cache for loaded window batches.
    let mut batch_cache: HashMap<(String, u64), RecordBatch> = HashMap::new();

    match &data_type {
        DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => {
            let mut builder = StringBuilder::with_capacity(total_rows, total_rows * 16);
            for entry in entries {
                match entry {
                    OutputEntry::Matched(m) => {
                        let key = (m.chrom.clone(), m.window_id);
                        if !batch_cache.contains_key(&key) {
                            if let Some(batch) = store.get_window(&m.chrom, m.window_id)? {
                                batch_cache.insert(key.clone(), batch);
                            }
                        }
                        let batch = batch_cache.get(&key);
                        for &i in &m.cache_rows {
                            if let Some(batch) = batch {
                                let col = cache_col_idx.map(|idx| batch.column(idx));
                                if let Some(col) = col {
                                    if col.is_null(i) {
                                        builder.append_null();
                                    } else {
                                        builder.append_value(get_string_value(col, i, col_name)?);
                                    }
                                } else {
                                    builder.append_null();
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    OutputEntry::Unmatched(_) => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int64 => {
            let mut builder = Int64Builder::with_capacity(total_rows);
            for entry in entries {
                match entry {
                    OutputEntry::Matched(m) => {
                        let key = (m.chrom.clone(), m.window_id);
                        if !batch_cache.contains_key(&key) {
                            if let Some(batch) = store.get_window(&m.chrom, m.window_id)? {
                                batch_cache.insert(key.clone(), batch);
                            }
                        }
                        let batch = batch_cache.get(&key);
                        for &i in &m.cache_rows {
                            if let Some(batch) = batch {
                                let col = cache_col_idx.map(|idx| batch.column(idx));
                                if let Some(col) = col {
                                    if col.is_null(i) {
                                        builder.append_null();
                                    } else {
                                        builder.append_value(
                                            get_int64_column(col, col_name)?.value(i),
                                        );
                                    }
                                } else {
                                    builder.append_null();
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    OutputEntry::Unmatched(_) => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int8 => {
            let mut builder = Int8Builder::with_capacity(total_rows);
            for entry in entries {
                match entry {
                    OutputEntry::Matched(m) => {
                        let key = (m.chrom.clone(), m.window_id);
                        if !batch_cache.contains_key(&key) {
                            if let Some(batch) = store.get_window(&m.chrom, m.window_id)? {
                                batch_cache.insert(key.clone(), batch);
                            }
                        }
                        let batch = batch_cache.get(&key);
                        for &i in &m.cache_rows {
                            if let Some(batch) = batch {
                                let col = cache_col_idx.map(|idx| batch.column(idx));
                                if let Some(col) = col {
                                    if col.is_null(i) {
                                        builder.append_null();
                                    } else {
                                        builder
                                            .append_value(get_int8_column(col, col_name)?.value(i));
                                    }
                                } else {
                                    builder.append_null();
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    OutputEntry::Unmatched(_) => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        _ => Ok(Arc::new(NullArray::new(total_rows)) as ArrayRef),
    }
}

fn get_int64_column<'a>(col: &'a ArrayRef, column_name: &str) -> Result<&'a Int64Array> {
    col.as_any().downcast_ref::<Int64Array>().ok_or_else(|| {
        DataFusionError::Execution(format!(
            "column '{column_name}' expected Int64 array, got {:?}",
            col.data_type()
        ))
    })
}

fn get_int8_column<'a>(col: &'a ArrayRef, column_name: &str) -> Result<&'a Int8Array> {
    col.as_any().downcast_ref::<Int8Array>().ok_or_else(|| {
        DataFusionError::Execution(format!(
            "column '{column_name}' expected Int8 array, got {:?}",
            col.data_type()
        ))
    })
}

fn get_string_value<'a>(col: &'a ArrayRef, i: usize, column_name: &str) -> Result<&'a str> {
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        Ok(arr.value(i))
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        Ok(arr.value(i))
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::LargeStringArray>()
    {
        Ok(arr.value(i))
    } else {
        Err(DataFusionError::Execution(format!(
            "column '{column_name}' expected string array, got {:?}",
            col.data_type()
        )))
    }
}

fn get_string_column(col: &ArrayRef, column_name: &str) -> Result<Vec<Option<String>>> {
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect())
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect())
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::LargeStringArray>()
    {
        Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect())
    } else {
        Err(DataFusionError::Execution(format!(
            "column '{column_name}' expected string array, got {:?}",
            col.data_type()
        )))
    }
}

fn normalize_vcf_coords(
    start: i64,
    end: i64,
    vcf_zero_based: bool,
    cache_zero_based: bool,
) -> (i64, i64) {
    if vcf_zero_based == cache_zero_based {
        (start, end)
    } else if vcf_zero_based {
        (start + 1, end) // 0-based half-open -> 1-based closed
    } else {
        (start - 1, end) // 1-based closed -> 0-based half-open
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
            Poll::Ready(None) => Poll::Ready(None),
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
            cache_rows: vec![0],
        };
        let entries = vec![OutputEntry::Matched(&matched)];
        let mut window_columns: HashMap<(String, u64), Vec<ArrayRef>> = HashMap::new();
        window_columns.insert(
            ("1".to_string(), 0),
            vec![Arc::new(StringArray::from(vec!["not_an_int64"])) as ArrayRef],
        );

        let err = build_cache_column_v1(
            &entries,
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
    fn test_process_batch_rejects_non_int64_start_column() {
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
        );

        let err = stream.process_batch(&vcf_batch).unwrap_err().to_string();
        assert!(err.contains("start"));
        assert!(err.contains("Int64"));
    }
}
