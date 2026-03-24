//! FusedArrayTransformExec - Physical execution plan for fused array transformation.
//!
//! This execution plan processes arrays element-wise without unnesting to intermediate rows,
//! achieving bounded memory usage regardless of input size.

use std::any::Any;
use std::fmt::{self, Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};
use std::time::Instant;

use datafusion::arrow::array::{
    Array, ArrayRef, BooleanArray, ListArray, MutableArrayData, RecordBatch,
};
use datafusion::arrow::buffer::OffsetBuffer;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::{EquivalenceProperties, PhysicalExpr};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::metrics::{BaselineMetrics, ExecutionPlanMetricsSet, MetricsSet};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use futures::{Stream, StreamExt, ready};

use crate::common::{build_element_schema_from_arrow, extract_list_element_type};

/// Physical execution plan for fused array transformation.
///
/// This operator takes an input with array columns and applies transformations
/// element-wise, producing output arrays without materializing intermediate unnested rows.
///
/// ## Physical-Level CSE
///
/// Before evaluating transform expressions, the operator detects common
/// sub-expression trees across all expressions (by comparing their Display
/// representation). Common sub-expressions are pre-evaluated once and added
/// as temporary columns to the batch, so subsequent expression evaluations
/// reference the cached result instead of recomputing.
pub struct FusedArrayTransformExec {
    /// Input execution plan
    input: Arc<dyn ExecutionPlan>,
    /// Names of array columns to transform
    array_columns: Vec<String>,
    /// Names of columns to pass through unchanged
    passthrough_columns: Vec<String>,
    /// Names for output columns
    output_columns: Vec<String>,
    /// Physical transformation expressions
    transform_exprs: Vec<Arc<dyn PhysicalExpr>>,
    /// Output schema
    schema: SchemaRef,
    /// Cached plan properties
    cache: PlanProperties,
    /// Execution metrics
    metrics: ExecutionPlanMetricsSet,
}

impl Debug for FusedArrayTransformExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("FusedArrayTransformExec")
            .field("array_columns", &self.array_columns)
            .field("passthrough_columns", &self.passthrough_columns)
            .field("output_columns", &self.output_columns)
            .field("schema", &self.schema)
            .finish()
    }
}

impl FusedArrayTransformExec {
    /// Create a new FusedArrayTransformExec.
    pub fn try_new(
        input: Arc<dyn ExecutionPlan>,
        array_columns: Vec<String>,
        passthrough_columns: Vec<String>,
        output_columns: Vec<String>,
        transform_exprs: Vec<Arc<dyn PhysicalExpr>>,
    ) -> Result<Self> {
        let input_schema = input.schema();
        let mut fields: Vec<Field> = Vec::new();

        // Add passthrough columns
        for col_name in &passthrough_columns {
            let idx = input_schema.index_of(col_name)?;
            let field = input_schema.field(idx).clone();
            fields.push(field);
        }

        // Build element-level schema for type inference
        let element_schema = build_element_schema_from_arrow(&input_schema, &array_columns)?;

        // Add output array columns - infer type from corresponding transform expression
        for (i, output_col) in output_columns.iter().enumerate() {
            let element_type = infer_output_element_type(
                i,
                output_col,
                &transform_exprs,
                &element_schema,
                &array_columns,
                &input_schema,
            )?;

            fields.push(Field::new(
                output_col,
                DataType::List(Arc::new(Field::new("item", element_type, true))),
                true,
            ));
        }

        let schema = Arc::new(Schema::new(fields));
        let cache = Self::compute_properties(&input, schema.clone());

        Ok(Self {
            input,
            array_columns,
            passthrough_columns,
            output_columns,
            transform_exprs,
            schema,
            cache,
            metrics: ExecutionPlanMetricsSet::new(),
        })
    }

    fn compute_properties(input: &Arc<dyn ExecutionPlan>, schema: SchemaRef) -> PlanProperties {
        let eq_properties = EquivalenceProperties::new(schema);
        PlanProperties::new(
            eq_properties,
            input.output_partitioning().clone(),
            EmissionType::Incremental,
            Boundedness::Bounded,
        )
    }
}

impl DisplayAs for FusedArrayTransformExec {
    fn fmt_as(&self, t: DisplayFormatType, f: &mut Formatter) -> fmt::Result {
        match t {
            DisplayFormatType::Default
            | DisplayFormatType::Verbose
            | DisplayFormatType::TreeRender => {
                write!(
                    f,
                    "FusedArrayTransformExec: arrays=[{}], passthrough=[{}], outputs=[{}]",
                    self.array_columns.join(", "),
                    self.passthrough_columns.join(", "),
                    self.output_columns.join(", ")
                )
            }
        }
    }
}

impl ExecutionPlan for FusedArrayTransformExec {
    fn name(&self) -> &str {
        "FusedArrayTransformExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
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
        if children.len() != 1 {
            return Err(DataFusionError::Plan(
                "FusedArrayTransformExec requires exactly 1 child".to_string(),
            ));
        }
        Ok(Arc::new(Self::try_new(
            children[0].clone(),
            self.array_columns.clone(),
            self.passthrough_columns.clone(),
            self.output_columns.clone(),
            self.transform_exprs.clone(),
        )?))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let baseline_metrics = BaselineMetrics::new(&self.metrics, partition);
        let input_stream = self.input.execute(partition, context)?;
        Ok(Box::pin(FusedArrayTransformStream::new(
            input_stream,
            self.schema.clone(),
            self.array_columns.clone(),
            self.passthrough_columns.clone(),
            self.transform_exprs.clone(),
            baseline_metrics,
        )))
    }

    fn metrics(&self) -> Option<MetricsSet> {
        Some(self.metrics.clone_inner())
    }
}

/// Stream that processes batches with fused array transformation.
struct FusedArrayTransformStream {
    input: SendableRecordBatchStream,
    schema: SchemaRef,
    array_columns: Vec<String>,
    passthrough_columns: Vec<String>,
    transform_exprs: Vec<Arc<dyn PhysicalExpr>>,
    baseline_metrics: BaselineMetrics,
}

impl FusedArrayTransformStream {
    fn new(
        input: SendableRecordBatchStream,
        schema: SchemaRef,
        array_columns: Vec<String>,
        passthrough_columns: Vec<String>,
        transform_exprs: Vec<Arc<dyn PhysicalExpr>>,
        baseline_metrics: BaselineMetrics,
    ) -> Self {
        Self {
            input,
            schema,
            array_columns,
            passthrough_columns,
            transform_exprs,
            baseline_metrics,
        }
    }

    /// Process a single input batch.
    fn process_batch(&self, batch: &RecordBatch) -> Result<RecordBatch> {
        let _timer = self.baseline_metrics.elapsed_compute().timer();
        let t_start = Instant::now();

        let num_rows = batch.num_rows();
        if num_rows == 0 {
            return Ok(RecordBatch::new_empty(self.schema.clone()));
        }
        let input_schema = batch.schema();

        // Collect passthrough columns (we may later filter some rows out)
        let mut passthrough_arrays: Vec<ArrayRef> =
            Vec::with_capacity(self.passthrough_columns.len());

        // Collect passthrough columns
        for col_name in &self.passthrough_columns {
            let idx = input_schema.index_of(col_name)?;
            passthrough_arrays.push(batch.column(idx).clone());
        }

        // If there are no array columns, just passthrough
        if self.array_columns.is_empty() {
            return RecordBatch::try_new(self.schema.clone(), passthrough_arrays)
                .map_err(DataFusionError::from);
        }

        // Extract all array columns (supports List, LargeList, FixedSizeList)
        let list_arrays: Vec<(&String, ArrayRef)> = self
            .array_columns
            .iter()
            .map(|col_name| {
                let idx = input_schema.index_of(col_name)?;
                let col = batch.column(idx).clone();
                // Verify it's a list type
                match col.data_type() {
                    DataType::List(_) | DataType::LargeList(_) | DataType::FixedSizeList(_, _) => {
                        Ok((col_name, col))
                    }
                    other => Err(DataFusionError::Plan(format!(
                        "Column '{col_name}' has type {other}, expected List, LargeList, or FixedSizeList"
                    ))),
                }
            })
            .collect::<Result<Vec<_>>>()?;

        // Compute, for each input row, the maximum array length across all unnest columns.
        // Rows where this maximum length is zero correspond to rows that would produce
        // no output rows after UNNEST + GROUP BY and thus must be dropped to match
        // baseline semantics.
        //
        // VECTORIZED: use ListArray::offsets() directly instead of per-row extract_list_values()
        // which creates a sliced ArrayRef per row. Offset arithmetic gives us lengths in O(1).
        let mut row_max_lengths: Vec<usize> = vec![0; num_rows];
        let mut row_mask: Vec<bool> = vec![false; num_rows];

        for (_col_name, array) in &list_arrays {
            if let Some(list_arr) = array.as_any().downcast_ref::<ListArray>() {
                let offsets = list_arr.offsets();
                for row_idx in 0..num_rows {
                    if !list_arr.is_null(row_idx) {
                        let len = (offsets[row_idx + 1] - offsets[row_idx]) as usize;
                        if len > row_max_lengths[row_idx] {
                            row_max_lengths[row_idx] = len;
                        }
                    }
                }
            } else {
                // Fallback for LargeList/FixedSizeList: per-row extraction
                for row_idx in 0..num_rows {
                    let values_opt =
                        crate::common::extract_list_values(array.as_ref(), row_idx, _col_name)?;
                    if let Some(values) = values_opt {
                        let len = values.len();
                        if len > row_max_lengths[row_idx] {
                            row_max_lengths[row_idx] = len;
                        }
                    }
                }
            }
        }
        for row_idx in 0..num_rows {
            row_mask[row_idx] = row_max_lengths[row_idx] > 0;
        }

        let num_output_rows = row_mask.iter().filter(|&&b| b).count();

        // If no rows would survive UNNEST, return an empty batch with the correct schema.
        if num_output_rows == 0 {
            let empty_columns: Vec<ArrayRef> = self
                .schema
                .fields()
                .iter()
                .map(|f| datafusion::arrow::array::new_empty_array(f.data_type()))
                .collect();
            return RecordBatch::try_new(self.schema.clone(), empty_columns)
                .map_err(DataFusionError::from);
        }

        // Filter passthrough columns to only include rows that have at least one element
        // after UNNEST (max_len > 0).
        let bool_mask = BooleanArray::from(row_mask.clone());
        let mut output_columns: Vec<ArrayRef> =
            Vec::with_capacity(self.passthrough_columns.len() + self.transform_exprs.len());

        for array in passthrough_arrays {
            let filtered =
                datafusion::arrow::compute::filter(array.as_ref(), &bool_mask).map_err(|e| {
                    DataFusionError::Execution(format!("Error filtering passthrough column: {e}"))
                })?;
            output_columns.push(filtered);
        }

        let t_after_passthrough = Instant::now();

        // Check if we have expressions to evaluate or do identity transform
        if self.transform_exprs.is_empty() {
            // No transform expressions: identity transform (pass through arrays unchanged)
            for col_name in &self.array_columns {
                let idx = input_schema.index_of(col_name)?;
                let col = batch.column(idx);

                let filtered_col = datafusion::arrow::compute::filter(col.as_ref(), &bool_mask)?;
                let output_array = self.apply_identity_transform(&filtered_col, col_name)?;
                output_columns.push(output_array);
            }
        } else {
            // Transform expressions: evaluate element-wise
            let transformed =
                self.apply_transforms(batch, num_rows, &list_arrays, &row_max_lengths, &row_mask)?;
            output_columns.extend(transformed);
        }

        let t_after_transform = Instant::now();
        let result =
            RecordBatch::try_new(self.schema.clone(), output_columns).map_err(DataFusionError::from);

        let t_end = Instant::now();
        log::debug!(
            "FusedArrayTransformExec::process_batch: num_rows={}, \
             passthrough={:.3}ms, transform={:.3}ms, total={:.3}ms",
            num_rows,
            (t_after_passthrough - t_start).as_secs_f64() * 1000.0,
            (t_after_transform - t_after_passthrough).as_secs_f64() * 1000.0,
            (t_end - t_start).as_secs_f64() * 1000.0,
        );

        self.baseline_metrics.record_output(result.as_ref().map(|b| b.num_rows()).unwrap_or(0));
        result
    }

    /// Apply transform expressions using batched evaluation.
    ///
    /// Instead of building a mini-batch per row and evaluating expressions N times,
    /// this builds ONE large "unnested" RecordBatch for the entire input, evaluates
    /// each expression once, then splits results back into ListArrays using offsets.
    fn apply_transforms(
        &self,
        batch: &RecordBatch,
        num_rows: usize,
        list_arrays: &[(&String, ArrayRef)],
        row_max_lengths: &[usize],
        row_mask: &[bool],
    ) -> Result<Vec<ArrayRef>> {
        // Build the unnested schema once
        let unnested_schema =
            build_mini_batch_schema(list_arrays, batch, &self.passthrough_columns)?;

        // Compute total unnested rows and per-row offsets
        let mut offsets: Vec<i32> = Vec::with_capacity(num_rows + 1);
        let mut null_bitmap: Vec<bool> = Vec::with_capacity(num_rows);
        offsets.push(0);
        let mut total_unnested_rows: usize = 0;

        for row_idx in 0..num_rows {
            if !row_mask[row_idx] {
                continue;
            }
            let max_len = row_max_lengths[row_idx];
            total_unnested_rows += max_len;
            offsets.push(total_unnested_rows as i32);
            null_bitmap.push(true);
        }

        if total_unnested_rows == 0 {
            // All rows filtered out, return empty list arrays
            return self
                .transform_exprs
                .iter()
                .enumerate()
                .map(|(expr_idx, _)| {
                    self.build_list_array_from_values(expr_idx, vec![], vec![0], vec![])
                })
                .collect();
        }

        // Build the single large unnested RecordBatch
        let t_unnest_start = Instant::now();
        let mut unnested_batch = build_batched_unnest(
            list_arrays,
            row_max_lengths,
            row_mask,
            num_rows,
            total_unnested_rows,
            &unnested_schema,
            batch,
            &self.passthrough_columns,
        )?;
        let t_unnest_end = Instant::now();

        // Physical-level CSE: find common sub-expression trees across all
        // transform expressions, pre-evaluate them, and add as cached columns.
        // Returns the enriched batch and a cache mapping Display strings to column indices.
        let t_cse_start = Instant::now();
        let (enriched_batch, cse_cache) = if !self.transform_exprs.is_empty() {
            apply_physical_cse(&self.transform_exprs, unnested_batch)?
        } else {
            (unnested_batch, std::collections::HashMap::new())
        };
        unnested_batch = enriched_batch;
        let t_cse_end = Instant::now();

        log::debug!(
            "FusedArrayTransformExec::apply_transforms: total_unnested_rows={}, \
             unnest_build={:.3}ms, cse={:.3}ms (cached {} sub-exprs), batch_cols={}, num_exprs={}",
            total_unnested_rows,
            (t_unnest_end - t_unnest_start).as_secs_f64() * 1000.0,
            (t_cse_end - t_cse_start).as_secs_f64() * 1000.0,
            cse_cache.len(),
            unnested_batch.num_columns(),
            self.transform_exprs.len(),
        );

        // Rewrite transform expressions to reference cached columns
        let final_exprs: Vec<Arc<dyn PhysicalExpr>> = if cse_cache.is_empty() {
            self.transform_exprs.clone()
        } else {
            self.transform_exprs
                .iter()
                .map(|expr| rewrite_expr_with_cache(expr, &cse_cache))
                .collect()
        };

        // Evaluate each final transform expression on the CSE-enriched batch
        final_exprs
            .iter()
            .enumerate()
            .map(|(expr_idx, expr)| {
                let t_eval_start = Instant::now();
                let result = expr.evaluate(&unnested_batch)?;
                let result_array = result.into_array(unnested_batch.num_rows())?;

                // Split the flat result array back into a ListArray using offsets
                let offset_buffer = OffsetBuffer::new(offsets.clone().into());
                let list_field =
                    Arc::new(Field::new("item", result_array.data_type().clone(), true));
                let null_buffer =
                    datafusion::arrow::buffer::NullBuffer::from(null_bitmap.clone());
                let list_array = ListArray::try_new(
                    list_field,
                    offset_buffer,
                    result_array,
                    Some(null_buffer),
                )?;
                log::debug!(
                    "  expr[{}] eval={:.3}ms",
                    expr_idx,
                    (Instant::now() - t_eval_start).as_secs_f64() * 1000.0,
                );
                Ok(Arc::new(list_array) as ArrayRef)
            })
            .collect()
    }

    /// Build a ListArray from collected value arrays.
    fn build_list_array_from_values(
        &self,
        expr_idx: usize,
        all_values: Vec<ArrayRef>,
        offsets: Vec<i32>,
        null_bitmap: Vec<bool>,
    ) -> Result<ArrayRef> {
        let values_array = if all_values.is_empty() {
            let output_field = self.schema.field(self.passthrough_columns.len() + expr_idx);
            let element_type = match output_field.data_type() {
                DataType::List(inner) => inner.data_type().clone(),
                other => {
                    return Err(DataFusionError::Plan(format!(
                        "Expected List type, got {other}"
                    )));
                }
            };
            datafusion::arrow::array::new_empty_array(&element_type)
        } else {
            let refs: Vec<&dyn Array> = all_values.iter().map(|a| a.as_ref()).collect();
            datafusion::arrow::compute::concat(&refs)?
        };

        let offset_buffer = OffsetBuffer::new(offsets.into());
        let list_field = Arc::new(Field::new("item", values_array.data_type().clone(), true));
        let null_buffer = datafusion::arrow::buffer::NullBuffer::from(null_bitmap);

        let list_array =
            ListArray::try_new(list_field, offset_buffer, values_array, Some(null_buffer))?;

        Ok(Arc::new(list_array))
    }

    /// Apply identity transformation (pass through arrays unchanged).
    /// Supports List, LargeList, and FixedSizeList types.
    fn apply_identity_transform(&self, array: &ArrayRef, col_name: &str) -> Result<ArrayRef> {
        // Verify the array is a list type
        match array.data_type() {
            DataType::List(_) | DataType::LargeList(_) | DataType::FixedSizeList(_, _) => {
                Ok(array.clone())
            }
            other => Err(DataFusionError::Plan(format!(
                "Column '{col_name}' has type {other}, expected List, LargeList, or FixedSizeList"
            ))),
        }
    }
}

impl Stream for FusedArrayTransformStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        loop {
            let batch_opt = ready!(self.input.poll_next_unpin(cx));
            match batch_opt {
                Some(Ok(batch)) => {
                    if batch.num_rows() == 0 {
                        continue;
                    }
                    let result = self.process_batch(&batch);
                    return Poll::Ready(Some(result));
                }
                Some(Err(e)) => return Poll::Ready(Some(Err(e))),
                None => return Poll::Ready(None),
            }
        }
    }
}

impl RecordBatchStream for FusedArrayTransformStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

/// Build the schema for mini-batches. This is computed once and reused for all rows.
///
/// Since the schema structure (column names and element types) is the same for all rows,
/// we avoid creating a new Schema object for every row.
///
/// The schema includes both unnested array columns and passthrough columns,
/// allowing transform expressions to reference passthrough columns.
fn build_mini_batch_schema(
    list_arrays: &[(&String, ArrayRef)],
    batch: &RecordBatch,
    passthrough_columns: &[String],
) -> Result<SchemaRef> {
    let mut fields: Vec<Field> = Vec::with_capacity(list_arrays.len() + passthrough_columns.len());

    // Add unnested array column element types
    for (col_name, array) in list_arrays {
        let element_type = extract_list_element_type(array.data_type(), col_name)?;
        fields.push(Field::new(col_name.as_str(), element_type, true));
    }

    // Add passthrough columns with their original types
    let input_schema = batch.schema();
    for col_name in passthrough_columns {
        let idx = input_schema.index_of(col_name)?;
        let field = input_schema.field(idx);
        fields.push(Field::new(
            col_name.as_str(),
            field.data_type().clone(),
            true,
        ));
    }

    Ok(Arc::new(Schema::new(fields)))
}

/// Build a single large unnested RecordBatch for all rows at once.
///
/// Instead of creating one mini-batch per row, this concatenates all rows' array
/// elements into a single RecordBatch and broadcasts passthrough columns in bulk.
/// This enables single-pass expression evaluation instead of per-row evaluation.
fn build_batched_unnest(
    list_arrays: &[(&String, ArrayRef)],
    row_max_lengths: &[usize],
    row_mask: &[bool],
    num_rows: usize,
    total_unnested_rows: usize,
    schema: &SchemaRef,
    batch: &RecordBatch,
    passthrough_columns: &[String],
) -> Result<RecordBatch> {
    let mut columns: Vec<ArrayRef> =
        Vec::with_capacity(list_arrays.len() + passthrough_columns.len());

    // Build unnested array columns by concatenating all rows' elements
    for (col_idx, (col_name, array)) in list_arrays.iter().enumerate() {
        let element_type = schema.field(col_idx).data_type();
        let unnested_col = build_batched_unnest_column(
            array.as_ref(),
            row_max_lengths,
            row_mask,
            num_rows,
            total_unnested_rows,
            element_type,
            col_name,
        )?;
        columns.push(unnested_col);
    }

    // Build broadcasted passthrough columns in bulk using a single take() per column
    let input_schema = batch.schema();

    // Pre-compute the indices array once: for each surviving row, repeat its index max_len times
    let mut indices: Vec<u64> = Vec::with_capacity(total_unnested_rows);
    for row_idx in 0..num_rows {
        if !row_mask[row_idx] {
            continue;
        }
        let max_len = row_max_lengths[row_idx];
        indices.extend(std::iter::repeat(row_idx as u64).take(max_len));
    }
    let indices_array = datafusion::arrow::array::UInt64Array::from(indices);

    for col_name in passthrough_columns {
        let idx = input_schema.index_of(col_name)?;
        let col = batch.column(idx);
        let broadcasted =
            datafusion::arrow::compute::kernels::take::take(col.as_ref(), &indices_array, None)?;
        columns.push(broadcasted);
    }

    RecordBatch::try_new(schema.clone(), columns).map_err(DataFusionError::from)
}

/// Build a single unnested column by concatenating all rows' array elements.
///
/// For each surviving row, extracts the list values and pads shorter arrays with NULLs.
/// Uses MutableArrayData for efficient incremental construction.
///
/// VECTORIZED: For ListArray (the common case), uses offsets and the flat values buffer
/// directly with a single MutableArrayData source, avoiding per-row ArrayRef creation.
fn build_batched_unnest_column(
    array: &dyn Array,
    row_max_lengths: &[usize],
    row_mask: &[bool],
    num_rows: usize,
    total_unnested_rows: usize,
    element_type: &DataType,
    col_name: &str,
) -> Result<ArrayRef> {
    // Fast path for ListArray: use offsets + flat values directly
    if let Some(list_arr) = array.as_any().downcast_ref::<ListArray>() {
        let offsets = list_arr.offsets();
        let flat_values = list_arr.values();
        let flat_data = flat_values.to_data();

        // Single MutableArrayData over the flat values buffer (index 0) + null source
        let mut mutable = MutableArrayData::new(vec![&flat_data], true, total_unnested_rows);

        for row_idx in 0..num_rows {
            if !row_mask[row_idx] {
                continue;
            }
            let max_len = row_max_lengths[row_idx];
            if list_arr.is_null(row_idx) {
                mutable.extend_nulls(max_len);
            } else {
                let start = offsets[row_idx] as usize;
                let end = offsets[row_idx + 1] as usize;
                let values_len = end - start;
                mutable.extend(0, start, end);
                if max_len > values_len {
                    mutable.extend_nulls(max_len - values_len);
                }
            }
        }

        return Ok(datafusion::arrow::array::make_array(mutable.freeze()));
    }

    // Fallback for LargeList/FixedSizeList: per-row extraction
    let mut row_values: Vec<Option<ArrayRef>> = Vec::new();
    let mut row_lengths: Vec<usize> = Vec::new();

    for row_idx in 0..num_rows {
        if !row_mask[row_idx] {
            continue;
        }
        let max_len = row_max_lengths[row_idx];
        let values_opt = crate::common::extract_list_values(array, row_idx, col_name)?;
        row_values.push(values_opt);
        row_lengths.push(max_len);
    }

    if row_values.is_empty() {
        return Ok(datafusion::arrow::array::new_empty_array(element_type));
    }

    let template = datafusion::arrow::array::new_empty_array(element_type);
    let template_data = template.to_data();

    let value_data: Vec<_> = row_values
        .iter()
        .filter_map(|v| v.as_ref().map(|a| a.to_data()))
        .collect();
    let mut all_data_refs: Vec<&datafusion::arrow::array::ArrayData> = vec![&template_data];
    for d in &value_data {
        all_data_refs.push(d);
    }

    let mut mutable = MutableArrayData::new(all_data_refs, true, total_unnested_rows);

    let mut data_idx = 1_usize;
    for (i, values_opt) in row_values.iter().enumerate() {
        let max_len = row_lengths[i];
        match values_opt {
            None => {
                mutable.extend_nulls(max_len);
            }
            Some(values) => {
                let values_len = values.len();
                mutable.extend(data_idx, 0, values_len);
                if max_len > values_len {
                    mutable.extend_nulls(max_len - values_len);
                }
                data_idx += 1;
            }
        }
    }

    Ok(datafusion::arrow::array::make_array(mutable.freeze()))
}

/// Infer the output element type for a transformed array column.
///
/// Attempts type inference in the following order:
/// 1. From the transform expression's data type (if expression exists)
/// 2. From the corresponding array column's element type (as fallback)
///
/// # Errors
///
/// Returns an error if type inference fails and no fallback is available.
/// This is intentional to avoid silent schema mismatches caused by Float64 fallback.
fn infer_output_element_type(
    index: usize,
    output_col: &str,
    transform_exprs: &[Arc<dyn PhysicalExpr>],
    element_schema: &Schema,
    array_columns: &[String],
    input_schema: &SchemaRef,
) -> Result<DataType> {
    if index < transform_exprs.len() {
        // Try to infer from transform expression
        match transform_exprs[index].data_type(element_schema) {
            Ok(dt) => return Ok(dt),
            Err(e) => {
                // Try falling back to array column's element type
                if let Some(array_col) = array_columns.get(index) {
                    return get_array_element_type(array_col, input_schema).map_err(|_| {
                        DataFusionError::Plan(format!(
                            "Type inference failed for output '{output_col}': {e}. \
                             Could not determine element type from array column '{array_col}'."
                        ))
                    });
                }
                return Err(DataFusionError::Plan(format!(
                    "Type inference failed for output '{output_col}': {e}. \
                     No corresponding array column to fall back to."
                )));
            }
        }
    }

    if index < array_columns.len() {
        // No transform expression, use array column's element type
        return get_array_element_type(&array_columns[index], input_schema);
    }

    // No fallback available
    Err(DataFusionError::Plan(format!(
        "No type information for output '{output_col}' at index {index}. \
         Number of outputs exceeds number of array columns ({}).",
        array_columns.len()
    )))
}

/// Get the element type from an array column in the schema.
fn get_array_element_type(array_col: &str, schema: &SchemaRef) -> Result<DataType> {
    let idx = schema.index_of(array_col)?;
    let field = schema.field(idx);
    extract_list_element_type(field.data_type(), array_col)
}

/// Collect all sub-expressions from a physical expression tree, with their
/// Display representation as key. Returns (display_string, expr_arc) pairs.
fn collect_sub_exprs(
    expr: &Arc<dyn PhysicalExpr>,
    out: &mut Vec<(String, Arc<dyn PhysicalExpr>)>,
) {
    let display = format!("{expr}");
    // Only collect non-trivial sub-expressions (not leaf Column refs)
    if !expr.children().is_empty() {
        out.push((display, expr.clone()));
        for child in expr.children() {
            collect_sub_exprs(&child, out);
        }
    }
}

/// Rewrite a physical expression tree, replacing sub-expressions that match
/// a cached result with a Column reference to the pre-computed column.
fn rewrite_expr_with_cache(
    expr: &Arc<dyn PhysicalExpr>,
    cache: &std::collections::HashMap<String, usize>, // display_string -> column_index
) -> Arc<dyn PhysicalExpr> {
    let display = format!("{expr}");
    if let Some(&col_idx) = cache.get(&display) {
        // Replace with a Column reference to the pre-computed column
        return Arc::new(datafusion::physical_expr::expressions::Column::new(
            &format!("__cse_{col_idx}"),
            col_idx,
        ));
    }

    // Recursively rewrite children
    let children = expr.children();
    if children.is_empty() {
        return expr.clone();
    }

    let new_children: Vec<Arc<dyn PhysicalExpr>> = children
        .iter()
        .map(|child| rewrite_expr_with_cache(child, cache))
        .collect();

    expr.clone().with_new_children(new_children).unwrap_or_else(|_| expr.clone())
}

/// Apply physical-level Common Sub-Expression Elimination.
///
/// Walks all transform expression trees, finds sub-expressions that appear
/// 2+ times (either across different expressions OR within the same expression
/// tree), evaluates them once, adds results as new columns, and rewrites
/// expressions to reference the cached columns.
///
/// This handles both:
/// - Cross-expression CSE: e.g., `nc/pw/x` shared between PL and DS outputs
/// - Intra-expression CSE: e.g., `nc/pw/x` appearing 3× inside `make_array(fp0,fp1,fp2)`
fn apply_physical_cse(
    transform_exprs: &[Arc<dyn PhysicalExpr>],
    mut batch: RecordBatch,
) -> Result<(RecordBatch, std::collections::HashMap<String, usize>)> {
    use std::collections::HashMap;

    // Step 1: Collect all sub-expressions from all expression trees,
    // counting TOTAL occurrences (not unique sources). This catches both
    // cross-expression and intra-expression duplication.
    let mut occurrence_count: HashMap<String, (Arc<dyn PhysicalExpr>, usize)> = HashMap::new();

    for expr in transform_exprs.iter() {
        let mut sub_exprs = Vec::new();
        collect_sub_exprs(expr, &mut sub_exprs);

        for (display, sub_expr) in sub_exprs {
            occurrence_count
                .entry(display)
                .and_modify(|(_, count)| *count += 1)
                .or_insert_with(|| (sub_expr, 1));
        }
    }

    // Step 2: Find sub-expressions that appear 2+ times total.
    let mut shared: Vec<(String, Arc<dyn PhysicalExpr>)> = occurrence_count
        .into_iter()
        .filter(|(_, (_, count))| *count > 1)
        .map(|(display, (expr, _))| (display, expr))
        .collect();

    if shared.is_empty() {
        return Ok((batch, HashMap::new()));
    }

    // Sort by display string length ascending (evaluate simpler/inner exprs first)
    shared.sort_by_key(|(display, _)| display.len());

    // Step 3: Evaluate shared sub-expressions and add as columns.
    // Use a cache to avoid evaluating sub-expressions of already-cached expressions.
    let mut cache: HashMap<String, usize> = HashMap::new();
    let mut new_columns: Vec<ArrayRef> = batch.columns().to_vec();
    let mut new_fields: Vec<Arc<Field>> = batch.schema().fields().iter().cloned().collect();

    for (display, expr) in &shared {
        // Skip if this exact expression was already cached
        if cache.contains_key(display) {
            continue;
        }

        // Rewrite the expression to use any already-cached sub-expressions
        let rewritten = rewrite_expr_with_cache(expr, &cache);

        // Evaluate on current batch (which may have cached columns from prior iterations)
        let temp_schema = Arc::new(Schema::new(new_fields.clone()));
        let temp_batch = RecordBatch::try_new(temp_schema, new_columns.clone())?;

        let result = rewritten.evaluate(&temp_batch)?;
        let result_array = result.into_array(temp_batch.num_rows())?;

        let col_idx = new_fields.len();
        new_fields.push(Arc::new(Field::new(
            &format!("__cse_{col_idx}"),
            result_array.data_type().clone(),
            true,
        )));
        new_columns.push(result_array);
        cache.insert(display.clone(), col_idx);
    }

    if cache.is_empty() {
        return Ok((batch, cache));
    }

    log::debug!(
        "  CSE: cached {} shared sub-expressions as columns",
        cache.len()
    );

    // Build the enriched batch with cached columns
    let enriched_schema = Arc::new(Schema::new(new_fields));
    batch = RecordBatch::try_new(enriched_schema, new_columns)?;

    Ok((batch, cache))
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float64Builder, ListBuilder, StringArray};
    use datafusion::arrow::datatypes::{Field, Schema};
    use datafusion::physical_plan::test::TestMemoryExec;

    fn create_test_batch(
        row0_lista: Vec<f64>,
        row0_listb: Vec<f64>,
        row1_lista: Vec<f64>,
        row1_listb: Vec<f64>,
    ) -> RecordBatch {
        let mut list_builder_a = ListBuilder::new(Float64Builder::new());
        let mut list_builder_b = ListBuilder::new(Float64Builder::new());

        // Row 0
        for val in row0_lista {
            list_builder_a.values().append_value(val);
        }
        list_builder_a.append(true);

        for val in row0_listb {
            list_builder_b.values().append_value(val);
        }
        list_builder_b.append(true);

        // Row 1
        for val in row1_lista {
            list_builder_a.values().append_value(val);
        }
        list_builder_a.append(true);

        for val in row1_listb {
            list_builder_b.values().append_value(val);
        }
        list_builder_b.append(true);

        let arr_a = list_builder_a.finish();
        let arr_b = list_builder_b.finish();

        let metadata = StringArray::from(vec!["meta1", "meta2"]);

        let schema = Arc::new(Schema::new(vec![
            Field::new("metadata", DataType::Utf8, false),
            Field::new(
                "values_a",
                DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
                true,
            ),
            Field::new(
                "values_b",
                DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
                true,
            ),
        ]));

        RecordBatch::try_new(
            schema,
            vec![Arc::new(metadata), Arc::new(arr_a), Arc::new(arr_b)],
        )
        .unwrap()
    }

    macro_rules! create_test_batch {
        ($row0_lista: expr, $row0_listb: expr, $row1_lista: expr, $row1_listb: expr) => {
            create_test_batch($row0_lista, $row0_listb, $row1_lista, $row1_listb)
        };
        () => {
            create_test_batch(
                vec![1.0, 2.0, 3.0],
                vec![10.0, 20.0, 30.0],
                vec![4.0, 5.0],
                vec![40.0, 50.0],
            )
        };
    }

    #[tokio::test]
    async fn test_identity_transform() {
        let batch = create_test_batch!();
        let schema = batch.schema();

        let mem_exec = TestMemoryExec::try_new(&[vec![batch.clone()]], schema, None).unwrap();

        let fused = FusedArrayTransformExec::try_new(
            Arc::new(mem_exec),
            vec!["values_a".to_string()],
            vec!["metadata".to_string()],
            vec!["values_a_out".to_string()],
            vec![],
        )
        .unwrap();

        // Schema should have 2 fields: metadata + values_a_out
        assert_eq!(fused.schema().fields().len(), 2);

        let ctx = Arc::new(TaskContext::default());
        let mut stream = fused.execute(0, ctx).unwrap();
        let result_batch = stream.next().await.unwrap().unwrap();
        assert_eq!(result_batch.num_rows(), 2);
    }

    #[tokio::test]
    async fn test_identity_transform_with_empty_array() {
        let batch = create_test_batch!(vec![], vec![], vec![4.0, 5.0], vec![40.0, 50.0]);
        let schema = batch.schema();

        let mem_exec = TestMemoryExec::try_new(&[vec![batch.clone()]], schema, None).unwrap();

        let fused = FusedArrayTransformExec::try_new(
            Arc::new(mem_exec),
            vec!["values_a".to_string(), "values_b".to_string()],
            vec!["metadata".to_string()],
            vec!["values_a_out".to_string(), "values_b_out".to_string()],
            vec![],
        )
        .unwrap();

        // Schema should have 3 fields: metadata + values_a_out + values_b_out
        assert_eq!(fused.schema().fields().len(), 3);
        // row0 should be filtered out due to empty array, so only row1 remains
        let ctx = Arc::new(TaskContext::default());
        let mut stream = fused.execute(0, ctx).unwrap();
        let result_batch = stream.next().await.unwrap().unwrap();
        assert_eq!(result_batch.num_rows(), 1);
    }

    #[tokio::test]
    async fn test_execution() {
        let batch = create_test_batch!();
        let schema = batch.schema();

        let mem_exec = TestMemoryExec::try_new(&[vec![batch.clone()]], schema, None).unwrap();

        let fused = FusedArrayTransformExec::try_new(
            Arc::new(mem_exec),
            vec!["values_a".to_string(), "values_b".to_string()],
            vec!["metadata".to_string()],
            vec!["values_a_out".to_string(), "values_b_out".to_string()],
            vec![],
        )
        .unwrap();

        let ctx = Arc::new(TaskContext::default());
        let mut stream = fused.execute(0, ctx).unwrap();

        let result_batch = stream.next().await.unwrap().unwrap();

        assert_eq!(result_batch.num_rows(), 2);
        assert_eq!(result_batch.num_columns(), 3); // metadata + 2 outputs
    }

    /// Test that the vectorized offset-based row_max_lengths computation
    /// correctly handles null lists and variable-length arrays.
    #[tokio::test]
    async fn test_vectorized_row_max_lengths_with_nulls() {
        // Row 0: values_a=[1.0, 2.0], values_b=[10.0, 20.0, 30.0] → max_len=3
        // Row 1: values_a=NULL, values_b=[40.0] → max_len=1
        let mut list_a = ListBuilder::new(Float64Builder::new());
        list_a.values().append_value(1.0);
        list_a.values().append_value(2.0);
        list_a.append(true);
        list_a.append(false); // NULL

        let mut list_b = ListBuilder::new(Float64Builder::new());
        list_b.values().append_value(10.0);
        list_b.values().append_value(20.0);
        list_b.values().append_value(30.0);
        list_b.append(true);
        list_b.values().append_value(40.0);
        list_b.append(true);

        let metadata = StringArray::from(vec!["m1", "m2"]);

        let schema = Arc::new(Schema::new(vec![
            Field::new("metadata", DataType::Utf8, false),
            Field::new("va", DataType::List(Arc::new(Field::new("item", DataType::Float64, true))), true),
            Field::new("vb", DataType::List(Arc::new(Field::new("item", DataType::Float64, true))), true),
        ]));

        let batch = RecordBatch::try_new(
            schema,
            vec![Arc::new(metadata), Arc::new(list_a.finish()), Arc::new(list_b.finish())],
        ).unwrap();

        let batch_schema = batch.schema();
        let mem_exec = TestMemoryExec::try_new(&[vec![batch]], batch_schema, None).unwrap();

        let fused = FusedArrayTransformExec::try_new(
            Arc::new(mem_exec),
            vec!["va".to_string(), "vb".to_string()],
            vec!["metadata".to_string()],
            vec!["va_out".to_string(), "vb_out".to_string()],
            vec![],
        ).unwrap();

        let ctx = Arc::new(TaskContext::default());
        let mut stream = fused.execute(0, ctx).unwrap();
        let result = stream.next().await.unwrap().unwrap();

        // Both rows should survive (both have at least one non-null array)
        assert_eq!(result.num_rows(), 2);

        // Identity transform: arrays passed through with original lengths (no padding)
        let va_out = result.column(1).as_any().downcast_ref::<ListArray>().unwrap();
        assert_eq!(va_out.value(0).len(), 2); // original [1,2]
        assert!(va_out.is_null(1)); // NULL list stays NULL

        let vb_out = result.column(2).as_any().downcast_ref::<ListArray>().unwrap();
        assert_eq!(vb_out.value(0).len(), 3); // original [10,20,30]
        assert_eq!(vb_out.value(1).len(), 1); // original [40]
    }

    /// Test that the fast-path ListArray unnest column construction
    /// produces correct values using offsets directly.
    #[tokio::test]
    async fn test_vectorized_unnest_column_values() {
        // Row 0: [10.0, 20.0, 30.0]
        // Row 1: [40.0, 50.0]
        let batch = create_test_batch!(
            vec![1.0, 2.0, 3.0], vec![10.0, 20.0, 30.0],
            vec![4.0, 5.0], vec![40.0, 50.0]
        );
        let schema = batch.schema();

        let mem_exec = TestMemoryExec::try_new(&[vec![batch.clone()]], schema, None).unwrap();

        let fused = FusedArrayTransformExec::try_new(
            Arc::new(mem_exec),
            vec!["values_a".to_string()],
            vec!["metadata".to_string()],
            vec!["values_a_out".to_string()],
            vec![],
        ).unwrap();

        let ctx = Arc::new(TaskContext::default());
        let mut stream = fused.execute(0, ctx).unwrap();
        let result = stream.next().await.unwrap().unwrap();

        assert_eq!(result.num_rows(), 2);

        // values_a_out should preserve the original values
        let out = result.column(1).as_any().downcast_ref::<ListArray>().unwrap();

        // Row 0: [1.0, 2.0, 3.0]
        let row0 = out.value(0);
        let row0_f64 = row0.as_any().downcast_ref::<datafusion::arrow::array::Float64Array>().unwrap();
        assert_eq!(row0_f64.len(), 3);
        assert!((row0_f64.value(0) - 1.0).abs() < 1e-10);
        assert!((row0_f64.value(1) - 2.0).abs() < 1e-10);
        assert!((row0_f64.value(2) - 3.0).abs() < 1e-10);

        // Row 1: [4.0, 5.0]
        let row1 = out.value(1);
        let row1_f64 = row1.as_any().downcast_ref::<datafusion::arrow::array::Float64Array>().unwrap();
        assert_eq!(row1_f64.len(), 2);
        assert!((row1_f64.value(0) - 4.0).abs() < 1e-10);
        assert!((row1_f64.value(1) - 5.0).abs() < 1e-10);
    }
}
