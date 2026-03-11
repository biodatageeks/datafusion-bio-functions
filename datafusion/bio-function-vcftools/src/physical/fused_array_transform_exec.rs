//! FusedArrayTransformExec - Physical execution plan for fused array transformation.
//!
//! This execution plan processes arrays element-wise without unnesting to intermediate rows,
//! achieving bounded memory usage regardless of input size.

use std::any::Any;
use std::fmt::{self, Debug, Formatter};
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::arrow::array::{
    Array, ArrayRef, BooleanArray, ListArray, MutableArrayData, RecordBatch,
};
use datafusion::arrow::buffer::OffsetBuffer;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::{EquivalenceProperties, PhysicalExpr};
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use futures::{Stream, StreamExt, ready};
use log::warn;

/// Physical execution plan for fused array transformation.
///
/// This operator takes an input with array columns and applies transformations
/// element-wise, producing output arrays without materializing intermediate unnested rows.
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
        let element_schema = build_element_schema(&input_schema, &array_columns)?;

        // Helper to get element type from array column
        let get_array_element_type = |array_col: &str| -> Result<DataType> {
            let idx = input_schema.index_of(array_col)?;
            let input_field = input_schema.field(idx);
            match input_field.data_type() {
                DataType::List(inner) | DataType::LargeList(inner) => Ok(inner.data_type().clone()),
                DataType::FixedSizeList(inner, _) => Ok(inner.data_type().clone()),
                other => Err(DataFusionError::Plan(format!(
                    "Column '{array_col}' has type {other}, expected List type"
                ))),
            }
        };

        // Add output array columns - infer type from corresponding transform expression
        for (i, output_col) in output_columns.iter().enumerate() {
            // Infer element type from the transform expression
            // If type inference fails (e.g., column not in element schema), fall back to array column's type
            let element_type = if i < transform_exprs.len() {
                let inferred = transform_exprs[i].data_type(&element_schema);
                match inferred {
                    Ok(dt) => dt,
                    Err(e) => {
                        // Try falling back to array column's element type
                        if let Some(array_col) = array_columns.get(i) {
                            if let Ok(dt) = get_array_element_type(array_col) {
                                warn!(
                                    "Type inference failed for output '{output_col}': {e}. Falling back to array column type: {dt:?}"
                                );
                                dt
                            } else {
                                warn!(
                                    "Type inference failed for output '{output_col}': {e}. Using Float64 fallback."
                                );
                                DataType::Float64
                            }
                        } else {
                            warn!(
                                "Type inference failed for output '{output_col}': {e}. Using Float64 fallback."
                            );
                            DataType::Float64
                        }
                    }
                }
            } else if i < array_columns.len() {
                // No transform expression, use array column's element type
                get_array_element_type(&array_columns[i])?
            } else {
                // Default fallback for extra outputs without corresponding array columns
                warn!(
                    "No type information for output '{output_col}' at index {i}. Using Float64 fallback."
                );
                DataType::Float64
            };

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
        let input_stream = self.input.execute(partition, context)?;
        Ok(Box::pin(FusedArrayTransformStream::new(
            input_stream,
            self.schema.clone(),
            self.array_columns.clone(),
            self.passthrough_columns.clone(),
            self.transform_exprs.clone(),
        )))
    }
}

/// Stream that processes batches with fused array transformation.
struct FusedArrayTransformStream {
    input: SendableRecordBatchStream,
    schema: SchemaRef,
    array_columns: Vec<String>,
    passthrough_columns: Vec<String>,
    transform_exprs: Vec<Arc<dyn PhysicalExpr>>,
}

impl FusedArrayTransformStream {
    fn new(
        input: SendableRecordBatchStream,
        schema: SchemaRef,
        array_columns: Vec<String>,
        passthrough_columns: Vec<String>,
        transform_exprs: Vec<Arc<dyn PhysicalExpr>>,
    ) -> Self {
        Self {
            input,
            schema,
            array_columns,
            passthrough_columns,
            transform_exprs,
        }
    }

    /// Process a single input batch.
    fn process_batch(&self, batch: &RecordBatch) -> Result<RecordBatch> {
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

        // Extract all array columns as ListArrays
        let list_arrays: Vec<(&String, &ListArray)> = self
            .array_columns
            .iter()
            .map(|col_name| {
                let idx = input_schema.index_of(col_name)?;
                let col = batch.column(idx);
                let list_arr = col.as_any().downcast_ref::<ListArray>().ok_or_else(|| {
                    DataFusionError::Plan(format!("Column '{col_name}' is not a List type"))
                })?;
                Ok((col_name, list_arr))
            })
            .collect::<Result<Vec<_>>>()?;

        // Compute, for each input row, the maximum array length across all unnest columns.
        // Rows where this maximum length is zero correspond to rows that would produce
        // no output rows after UNNEST + GROUP BY and thus must be dropped to match
        // baseline semantics.
        let mut row_max_lengths: Vec<usize> = Vec::with_capacity(num_rows);
        let mut row_mask: Vec<bool> = Vec::with_capacity(num_rows);

        for row_idx in 0..num_rows {
            let mut max_len = 0_usize;

            for (_, list_arr) in &list_arrays {
                if list_arr.is_null(row_idx) {
                    continue;
                }
                let values = list_arr.value(row_idx);
                let len = values.len();
                if len > max_len {
                    max_len = len;
                }
            }

            row_max_lengths.push(max_len);
            row_mask.push(max_len > 0);
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
        let mut output_columns: Vec<ArrayRef> = Vec::with_capacity(
            self.passthrough_columns.len() + self.transform_exprs.len(),
        );

        for array in passthrough_arrays {
            let filtered =
                datafusion::arrow::compute::filter(array.as_ref(), &bool_mask).map_err(
                    |e| DataFusionError::Execution(format!("Error filtering passthrough column: {e}")),
                )?;
            output_columns.push(filtered);
        }

        // Check if we have expressions to evaluate or do identity transform
        if self.transform_exprs.is_empty() {
            // No transform expressions: identity transform (pass through arrays unchanged)
            for col_name in &self.array_columns {
                let idx = input_schema.index_of(col_name)?;
                let col = batch.column(idx);
                let list_arr = col.as_any().downcast_ref::<ListArray>().ok_or_else(|| {
                    DataFusionError::Plan(format!("Column '{col_name}' is not a List type"))
                })?;

                let output_array = self.apply_identity_transform(list_arr, num_rows)?;
                output_columns.push(output_array);
            }
        } else {
            // Transform expressions: evaluate element-wise
            let transformed =
                self.apply_transforms(batch, num_rows, &list_arrays, &row_max_lengths, &row_mask)?;
            output_columns.extend(transformed);
        }

        RecordBatch::try_new(self.schema.clone(), output_columns).map_err(DataFusionError::from)
    }

    /// Apply transform expressions element-wise to array columns.
    fn apply_transforms(
        &self,
        batch: &RecordBatch,
        num_rows: usize,
        list_arrays: &[(&String, &ListArray)],
        row_max_lengths: &[usize],
        row_mask: &[bool],
    ) -> Result<Vec<ArrayRef>> {
        // For each output expression, build a ListArray by evaluating element-wise
        let mut output_arrays: Vec<ArrayRef> = Vec::with_capacity(self.transform_exprs.len());

        for (expr_idx, expr) in self.transform_exprs.iter().enumerate() {
            // Collect evaluated arrays and offsets for each *included* row
            let mut all_values: Vec<ArrayRef> = Vec::with_capacity(num_rows);
            let mut offsets: Vec<i32> = vec![0];
            let mut null_bitmap: Vec<bool> = Vec::with_capacity(num_rows);
            let mut current_offset: i32 = 0;

            for row_idx in 0..num_rows {
                if !row_mask[row_idx] {
                    // This row would produce no UNNEST output, so it is not present
                    // in the fused result at all.
                    continue;
                }

                let max_len = row_max_lengths[row_idx];
                if max_len == 0 {
                    // Should not happen given the row_mask, but guard just in case.
                    continue;
                }

                // Build a mini-batch for this row's array elements, following
                // UNNEST semantics:
                // - Use the longest array length across all columns.
                // - Shorter arrays are padded with NULLs.
                // - NULL arrays become all-null arrays of length max_len.
                let mut columns: Vec<ArrayRef> = Vec::with_capacity(list_arrays.len());
                let mut fields: Vec<Field> = Vec::with_capacity(list_arrays.len());

                for (col_name, list_arr) in list_arrays {
                    let element_type = match list_arr.data_type() {
                        DataType::List(inner) => inner.data_type().clone(),
                        DataType::LargeList(inner) => inner.data_type().clone(),
                        DataType::FixedSizeList(inner, _) => inner.data_type().clone(),
                        other => {
                            return Err(DataFusionError::Plan(format!(
                                "Expected list type for '{col_name}', got {other}"
                            )));
                        }
                    };

                    let padded_values: ArrayRef = if list_arr.is_null(row_idx) {
                        // Entire array is NULL -> all-null values after UNNEST padding.
                        datafusion::arrow::array::new_null_array(&element_type, max_len)
                    } else {
                        let values = list_arr.value(row_idx);
                        let values_len = values.len();

                        if values_len == max_len {
                            values
                        } else {
                            // Pad shorter arrays with NULLs to reach max_len.
                            let array_data = values.to_data();
                            let mut mutable =
                                MutableArrayData::new(vec![&array_data], true, max_len);
                            // Copy existing values
                            mutable.extend(0, 0, values_len);
                            // Append NULLs for padding
                            if max_len > values_len {
                                mutable.extend_nulls(max_len - values_len);
                            }
                            datafusion::arrow::array::make_array(mutable.freeze())
                        }
                    };

                    fields.push(Field::new(col_name.as_str(), element_type, true));
                    columns.push(padded_values);
                }

                let schema = Arc::new(Schema::new(fields));
                let mini_batch =
                    RecordBatch::try_new(schema, columns).map_err(DataFusionError::from)?;

                // Evaluate the expression on the mini-batch
                let result = expr.evaluate(&mini_batch)?;
                let result_array = result.into_array(mini_batch.num_rows())?;

                current_offset += result_array.len() as i32;
                all_values.push(result_array);
                // The list itself is present (non-null) for this row; element nulls
                // are represented inside result_array.
                null_bitmap.push(true);
                offsets.push(current_offset);
            }

            // Concatenate all value arrays
            let values_array = if all_values.is_empty() {
                // Determine element type from schema
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

            // Build the ListArray
            let offset_buffer = OffsetBuffer::new(offsets.into());
            let list_field = Arc::new(Field::new("item", values_array.data_type().clone(), true));

            // Build null buffer
            let null_buffer = datafusion::arrow::buffer::NullBuffer::from(null_bitmap);

            let list_array =
                ListArray::try_new(list_field, offset_buffer, values_array, Some(null_buffer))?;

            output_arrays.push(Arc::new(list_array));
        }

        Ok(output_arrays)
    }

    /// Apply identity transformation (pass through arrays unchanged).
    fn apply_identity_transform(&self, list_arr: &ListArray, _num_rows: usize) -> Result<ArrayRef> {
        // For identity transform, just clone the array
        Ok(Arc::new(list_arr.clone()))
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

/// Build a schema representing the element types of array columns.
///
/// For each array column, we extract the element type and create a
/// column with the same name but the scalar (element) type. This schema
/// is used for inferring the output type of transform expressions.
fn build_element_schema(input_schema: &SchemaRef, array_columns: &[String]) -> Result<Schema> {
    let mut fields: Vec<Field> = Vec::with_capacity(array_columns.len());

    for col_name in array_columns {
        let idx = input_schema.index_of(col_name)?;
        let input_field = input_schema.field(idx);

        // Extract element type from List/LargeList/FixedSizeList
        let element_type = match input_field.data_type() {
            DataType::List(inner) | DataType::LargeList(inner) => inner.data_type().clone(),
            DataType::FixedSizeList(inner, _) => inner.data_type().clone(),
            other => {
                return Err(DataFusionError::Plan(format!(
                    "Column '{col_name}' has type {other}, expected List type"
                )));
            }
        };

        fields.push(Field::new(col_name, element_type, true));
    }

    Ok(Schema::new(fields))
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float64Builder, ListBuilder, StringArray};
    use datafusion::arrow::datatypes::{Field, Schema};
    use datafusion::physical_plan::test::TestMemoryExec;

    fn create_test_batch() -> RecordBatch {
        let mut list_builder_a = ListBuilder::new(Float64Builder::new());
        let mut list_builder_b = ListBuilder::new(Float64Builder::new());

        // Row 0: [1.0, 2.0, 3.0], [10.0, 20.0, 30.0]
        list_builder_a.values().append_value(1.0);
        list_builder_a.values().append_value(2.0);
        list_builder_a.values().append_value(3.0);
        list_builder_a.append(true);

        list_builder_b.values().append_value(10.0);
        list_builder_b.values().append_value(20.0);
        list_builder_b.values().append_value(30.0);
        list_builder_b.append(true);

        // Row 1: [4.0, 5.0], [40.0, 50.0]
        list_builder_a.values().append_value(4.0);
        list_builder_a.values().append_value(5.0);
        list_builder_a.append(true);

        list_builder_b.values().append_value(40.0);
        list_builder_b.values().append_value(50.0);
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

    #[tokio::test]
    async fn test_identity_transform() {
        let batch = create_test_batch();
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
    }

    #[tokio::test]
    async fn test_execution() {
        let batch = create_test_batch();
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
}
