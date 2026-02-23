//! KvLookupExec: ExecutionPlan that streams VCF batches and probes
//! a fjall KV store window-by-window for annotation.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::arrow::array::builder::{Int8Builder, Int64Builder, StringBuilder};
use datafusion::arrow::array::{Array, ArrayRef, AsArray, NullArray, RecordBatch};
use datafusion::arrow::datatypes::{DataType, Field, Int64Type, SchemaRef};
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
use crate::kv_store::VepKvStore;

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

        let mut fields: Vec<Arc<Field>> = input_schema.fields().iter().cloned().collect();
        for col_name in &cache_columns {
            if let Ok(field) = cache_schema.field_with_name(col_name) {
                fields.push(Arc::new(Field::new(
                    format!("cache_{}", field.name()),
                    field.data_type().clone(),
                    true,
                )));
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
            current_window: None,
            current_index: None,
        }
    }

    fn ensure_window(&mut self, chrom: &str, pos: i64) -> Result<()> {
        let wid = window_id_for_position(pos, self.window_size);
        let need_load = match &self.current_window {
            Some((c, w)) => c != chrom || *w != wid,
            None => true,
        };

        if need_load {
            let batch = self.store.get_window(chrom, wid)?;
            self.current_window = Some((chrom.to_string(), wid));
            self.current_index = match batch {
                Some(b) => Some(WindowAlleleIndex::from_batch(b)?),
                None => None,
            };
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

        let chroms = get_string_column(vcf_batch.column(chrom_idx));
        let starts = vcf_batch.column(start_idx).as_primitive::<Int64Type>();
        let ends = vcf_batch.column(end_idx).as_primitive::<Int64Type>();
        let refs = get_string_column(vcf_batch.column(ref_idx));
        let alts = get_string_column(vcf_batch.column(alt_idx));

        let cache_schema = self.store.schema().clone();
        let num_vcf_cols = vcf_schema.fields().len();
        let num_rows = vcf_batch.num_rows();

        // (vcf_row_idx, Option<(cache_batch, matched_row_indices)>).
        type MatchEntry = (usize, Option<(RecordBatch, Vec<usize>)>);
        let mut output_entries: Vec<MatchEntry> = Vec::new();

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
                    Some((index.batch().clone(), matches))
                }
            } else {
                None
            };

            output_entries.push((row, matched));
        }

        // Build output arrays.
        let total_output_rows: usize = output_entries
            .iter()
            .map(|(_, m)| match m {
                Some((_, indices)) => indices.len(),
                None => 1,
            })
            .sum();

        // Expand VCF row indices.
        let mut vcf_output_indices: Vec<u32> = Vec::with_capacity(total_output_rows);
        for (vcf_row, matched) in &output_entries {
            match matched {
                Some((_, indices)) => {
                    for _ in indices {
                        vcf_output_indices.push(*vcf_row as u32);
                    }
                }
                None => {
                    vcf_output_indices.push(*vcf_row as u32);
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
        for cache_col_name in &self.cache_columns {
            let cache_col_idx = cache_schema.index_of(cache_col_name).ok();
            let output_col = build_cache_column(
                &output_entries,
                cache_col_idx,
                cache_col_name,
                cache_schema.as_ref(),
                total_output_rows,
            )?;
            output_columns.push(output_col);
        }

        RecordBatch::try_new(self.schema.clone(), output_columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

type MatchEntry = (usize, Option<(RecordBatch, Vec<usize>)>);

fn build_cache_column(
    entries: &[MatchEntry],
    cache_col_idx: Option<usize>,
    col_name: &str,
    cache_schema: &datafusion::arrow::datatypes::Schema,
    total_rows: usize,
) -> Result<ArrayRef> {
    let data_type = cache_schema
        .field_with_name(col_name)
        .map(|f| f.data_type().clone())
        .unwrap_or(DataType::Utf8);

    match &data_type {
        DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => {
            let mut builder = StringBuilder::with_capacity(total_rows, total_rows * 16);
            for (_, matched) in entries {
                match matched {
                    Some((batch, indices)) => {
                        let col = cache_col_idx.map(|idx| batch.column(idx));
                        for &i in indices {
                            if let Some(col) = col {
                                if col.is_null(i) {
                                    builder.append_null();
                                } else {
                                    builder.append_value(get_string_value(col, i));
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    None => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int64 => {
            let mut builder = Int64Builder::with_capacity(total_rows);
            for (_, matched) in entries {
                match matched {
                    Some((batch, indices)) => {
                        let col = cache_col_idx.map(|idx| batch.column(idx));
                        for &i in indices {
                            if let Some(col) = col {
                                if col.is_null(i) {
                                    builder.append_null();
                                } else {
                                    builder.append_value(col.as_primitive::<Int64Type>().value(i));
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    None => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        DataType::Int8 => {
            let mut builder = Int8Builder::with_capacity(total_rows);
            for (_, matched) in entries {
                match matched {
                    Some((batch, indices)) => {
                        let col = cache_col_idx.map(|idx| batch.column(idx));
                        for &i in indices {
                            if let Some(col) = col {
                                if col.is_null(i) {
                                    builder.append_null();
                                } else {
                                    builder.append_value(
                                        col.as_primitive::<datafusion::arrow::datatypes::Int8Type>(
                                        )
                                        .value(i),
                                    );
                                }
                            } else {
                                builder.append_null();
                            }
                        }
                    }
                    None => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()) as ArrayRef)
        }
        _ => Ok(Arc::new(NullArray::new(total_rows)) as ArrayRef),
    }
}

fn get_string_value(col: &ArrayRef, i: usize) -> &str {
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        arr.value(i)
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        arr.value(i)
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::LargeStringArray>()
    {
        arr.value(i)
    } else {
        ""
    }
}

fn get_string_column(col: &ArrayRef) -> Vec<Option<String>> {
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect()
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect()
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::LargeStringArray>()
    {
        (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect()
    } else {
        vec![None; col.len()]
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
