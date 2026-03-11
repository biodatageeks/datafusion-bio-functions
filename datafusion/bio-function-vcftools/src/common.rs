//! Common utilities shared between logical and physical modules.
//!
//! This module contains shared types and functions used by both the logical
//! and physical components of the fused array transform optimization.

use std::sync::Arc;

use datafusion::arrow::array::{Array, ArrayRef, ListArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DFSchema, DataFusionError, Result};

/// Build a schema representing the element types of array columns.
///
/// For each array column, we extract the element type and create a
/// column with the same name but the scalar (element) type. This schema
/// is used for inferring the output type of transform expressions and
/// for converting logical expressions to physical expressions.
///
/// # Arguments
///
/// * `input_schema` - The schema containing array columns
/// * `array_columns` - Names of columns to extract element types from
///
/// # Errors
///
/// Returns an error if:
/// * A column name is not found in the schema
/// * A column has a non-list type (not List, LargeList, or FixedSizeList)
pub fn build_element_schema_from_arrow(
    input_schema: &Arc<Schema>,
    array_columns: &[String],
) -> Result<Schema> {
    let mut fields: Vec<Field> = Vec::with_capacity(array_columns.len());

    for col_name in array_columns {
        let idx = input_schema.index_of(col_name)?;
        let input_field = input_schema.field(idx);

        let element_type = extract_list_element_type(input_field.data_type(), col_name)?;
        fields.push(Field::new(col_name, element_type, true));
    }

    Ok(Schema::new(fields))
}

/// Build a schema representing the element types of array columns from a DFSchema.
///
/// Same as `build_element_schema_from_arrow` but works with DataFusion's DFSchema.
pub fn build_element_schema_from_df(
    input_schema: &DFSchema,
    array_columns: &[String],
) -> Result<DFSchema> {
    let mut fields: Vec<Field> = Vec::with_capacity(array_columns.len());

    for col_name in array_columns {
        let input_field = input_schema.field_with_unqualified_name(col_name)?;

        let element_type = extract_list_element_type(input_field.data_type(), col_name)?;
        fields.push(Field::new(col_name, element_type, true));
    }

    DFSchema::try_from(Schema::new(fields))
}

/// Extract the element type from a List/LargeList/FixedSizeList data type.
///
/// # Arguments
///
/// * `data_type` - The data type to extract from
/// * `col_name` - Column name for error messages
///
/// # Errors
///
/// Returns an error if the data type is not a list type.
pub fn extract_list_element_type(data_type: &DataType, col_name: &str) -> Result<DataType> {
    match data_type {
        DataType::List(inner) | DataType::LargeList(inner) => Ok(inner.data_type().clone()),
        DataType::FixedSizeList(inner, _) => Ok(inner.data_type().clone()),
        other => Err(DataFusionError::Plan(format!(
            "Column '{col_name}' has type {other}, expected List, LargeList, or FixedSizeList type"
        ))),
    }
}

/// Extract values from an array column at a specific row index, handling all list types.
///
/// Supports List, LargeList, and FixedSizeList array types.
///
/// # Arguments
///
/// * `column` - The array column containing list data
/// * `row_idx` - The row index to extract values from
/// * `col_name` - Column name for error messages
///
/// # Returns
///
/// * `None` if the value at row_idx is null
/// * `Some(array)` with the values at that row index
///
/// # Errors
///
/// Returns an error if the column is not a supported list type.
pub fn extract_list_values(
    column: &dyn Array,
    row_idx: usize,
    col_name: &str,
) -> Result<Option<ArrayRef>> {
    use datafusion::arrow::array::{FixedSizeListArray, LargeListArray};

    // Try List first (most common)
    if let Some(list_arr) = column.as_any().downcast_ref::<ListArray>() {
        if list_arr.is_null(row_idx) {
            return Ok(None);
        }
        return Ok(Some(list_arr.value(row_idx)));
    }

    // Try LargeList
    if let Some(list_arr) = column.as_any().downcast_ref::<LargeListArray>() {
        if list_arr.is_null(row_idx) {
            return Ok(None);
        }
        return Ok(Some(list_arr.value(row_idx)));
    }

    // Try FixedSizeList
    if let Some(list_arr) = column.as_any().downcast_ref::<FixedSizeListArray>() {
        if list_arr.is_null(row_idx) {
            return Ok(None);
        }
        return Ok(Some(list_arr.value(row_idx)));
    }

    Err(DataFusionError::Plan(format!(
        "Column '{col_name}' is not a List, LargeList, or FixedSizeList type"
    )))
}

/// Check if a row in a list array is null, supporting all list types.
pub fn is_list_null(column: &dyn Array, row_idx: usize, col_name: &str) -> Result<bool> {
    use datafusion::arrow::array::{FixedSizeListArray, LargeListArray};

    if let Some(list_arr) = column.as_any().downcast_ref::<ListArray>() {
        return Ok(list_arr.is_null(row_idx));
    }

    if let Some(list_arr) = column.as_any().downcast_ref::<LargeListArray>() {
        return Ok(list_arr.is_null(row_idx));
    }

    if let Some(list_arr) = column.as_any().downcast_ref::<FixedSizeListArray>() {
        return Ok(list_arr.is_null(row_idx));
    }

    Err(DataFusionError::Plan(format!(
        "Column '{col_name}' is not a List, LargeList, or FixedSizeList type"
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float64Array, Float64Builder, ListBuilder};

    #[test]
    fn test_extract_list_element_type_list() {
        let dt = DataType::List(Arc::new(Field::new("item", DataType::Float64, true)));
        let result = extract_list_element_type(&dt, "test").unwrap();
        assert_eq!(result, DataType::Float64);
    }

    #[test]
    fn test_extract_list_element_type_large_list() {
        let dt = DataType::LargeList(Arc::new(Field::new("item", DataType::Int32, true)));
        let result = extract_list_element_type(&dt, "test").unwrap();
        assert_eq!(result, DataType::Int32);
    }

    #[test]
    fn test_extract_list_element_type_fixed_size_list() {
        let dt = DataType::FixedSizeList(Arc::new(Field::new("item", DataType::Utf8, true)), 5);
        let result = extract_list_element_type(&dt, "test").unwrap();
        assert_eq!(result, DataType::Utf8);
    }

    #[test]
    fn test_extract_list_element_type_non_list() {
        let dt = DataType::Float64;
        let result = extract_list_element_type(&dt, "test");
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_list_values() {
        let mut builder = ListBuilder::new(Float64Builder::new());
        builder.values().append_value(1.0);
        builder.values().append_value(2.0);
        builder.append(true);
        builder.append(false); // null
        let list_arr = builder.finish();

        // Non-null row
        let values = extract_list_values(&list_arr, 0, "test").unwrap();
        assert!(values.is_some());
        let arr = values.unwrap();
        let f64_arr = arr.as_any().downcast_ref::<Float64Array>().unwrap();
        assert_eq!(f64_arr.len(), 2);

        // Null row
        let values = extract_list_values(&list_arr, 1, "test").unwrap();
        assert!(values.is_none());
    }
}
