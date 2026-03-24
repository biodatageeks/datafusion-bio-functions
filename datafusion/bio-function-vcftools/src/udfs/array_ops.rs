//! Array manipulation UDFs: list_index, list_replace_null, list_make_array.

use std::any::Any;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, AsArray, Float64Array, Float64Builder, Int32Array, Int32Builder, ListArray,
    StringArray,
};
use datafusion::arrow::buffer::OffsetBuffer;
use datafusion::arrow::datatypes::{DataType, Field, Float64Type, Int32Type};
use datafusion::common::{DataFusionError, Result};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, Volatility,
};

use super::{f64_field, i32_field, scalar_to_f64, scalar_to_i64, scalar_to_utf8, utf8_field};

// ============================================================================
// list_index(List<List<T>>, Int) → List<T>
//   Extract element at 1-based index from each inner list.
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListIndexUdf {
    signature: Signature,
}

impl ListIndexUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(2, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListIndexUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_index"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        // Input is List<List<T>>, output is List<T>
        match &arg_types[0] {
            DataType::List(inner) => Ok(inner.data_type().clone()),
            dt => Err(DataFusionError::Internal(format!(
                "list_index: expected List type, got {dt}"
            ))),
        }
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arr = args.args[0].clone().into_array(1)?;
        let outer_list = arr.as_list::<i32>();
        let idx = scalar_to_i64(&args.args[1])?;
        let k = (idx - 1) as usize; // 1-based → 0-based

        // outer_list.values() is the inner ListArray
        let inner_list = outer_list
            .values()
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| {
                DataFusionError::Internal("list_index: expected List<List<T>>".into())
            })?;

        let inner_offsets = inner_list.offsets();
        let flat_values = inner_list.values();
        let total_inner = inner_list.len();

        let result: ArrayRef = match flat_values.data_type() {
            DataType::Int32 => {
                let flat = flat_values.as_primitive::<Int32Type>();
                let mut builder = Int32Builder::with_capacity(total_inner);
                for i in 0..total_inner {
                    if inner_list.is_null(i) {
                        builder.append_null();
                        continue;
                    }
                    let start = inner_offsets[i] as usize;
                    let end = inner_offsets[i + 1] as usize;
                    if k < (end - start) && !flat.is_null(start + k) {
                        builder.append_value(flat.value(start + k));
                    } else {
                        builder.append_null();
                    }
                }
                Arc::new(ListArray::new(
                    i32_field(),
                    outer_list.offsets().clone(),
                    Arc::new(builder.finish()),
                    outer_list.nulls().cloned(),
                ))
            }
            DataType::Float64 => {
                let flat = flat_values.as_primitive::<Float64Type>();
                let mut builder = Float64Builder::with_capacity(total_inner);
                for i in 0..total_inner {
                    if inner_list.is_null(i) {
                        builder.append_null();
                        continue;
                    }
                    let start = inner_offsets[i] as usize;
                    let end = inner_offsets[i + 1] as usize;
                    if k < (end - start) && !flat.is_null(start + k) {
                        builder.append_value(flat.value(start + k));
                    } else {
                        builder.append_null();
                    }
                }
                Arc::new(ListArray::new(
                    f64_field(),
                    outer_list.offsets().clone(),
                    Arc::new(builder.finish()),
                    outer_list.nulls().cloned(),
                ))
            }
            dt => {
                return Err(DataFusionError::Internal(format!(
                    "list_index: unsupported inner element type {dt}"
                )));
            }
        };

        Ok(ColumnarValue::Array(result))
    }
}

pub fn list_index_udf() -> ScalarUDF {
    ScalarUDF::from(ListIndexUdf::new())
}

// ============================================================================
// list_replace_null(List<T>, T) → List<T>
//   Replace null elements with a default value.
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListReplaceNullUdf {
    signature: Signature,
}

impl ListReplaceNullUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(2, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListReplaceNullUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_replace_null"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        Ok(arg_types[0].clone())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arr = args.args[0].clone().into_array(1)?;
        let list = arr.as_list::<i32>();
        let flat = list.values();

        let result: ArrayRef = match flat.data_type() {
            DataType::Int32 => {
                let default = scalar_to_f64(&args.args[1])? as i32;
                let values = flat.as_primitive::<Int32Type>();
                let out: Int32Array = (0..values.len())
                    .map(|i| Some(if values.is_null(i) { default } else { values.value(i) }))
                    .collect();
                Arc::new(ListArray::new(
                    i32_field(),
                    list.offsets().clone(),
                    Arc::new(out),
                    list.nulls().cloned(),
                ))
            }
            DataType::Float64 => {
                let default = scalar_to_f64(&args.args[1])?;
                let values = flat.as_primitive::<Float64Type>();
                let out: Float64Array = (0..values.len())
                    .map(|i| Some(if values.is_null(i) { default } else { values.value(i) }))
                    .collect();
                Arc::new(ListArray::new(
                    f64_field(),
                    list.offsets().clone(),
                    Arc::new(out),
                    list.nulls().cloned(),
                ))
            }
            DataType::Utf8 => {
                let default = scalar_to_utf8(&args.args[1])?;
                let values = flat
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .ok_or_else(|| {
                        DataFusionError::Internal("list_replace_null: expected StringArray".into())
                    })?;
                let out: StringArray = (0..values.len())
                    .map(|i| {
                        Some(if values.is_null(i) {
                            default.clone()
                        } else {
                            values.value(i).to_string()
                        })
                    })
                    .collect();
                Arc::new(ListArray::new(
                    utf8_field(),
                    list.offsets().clone(),
                    Arc::new(out),
                    list.nulls().cloned(),
                ))
            }
            dt => {
                return Err(DataFusionError::Internal(format!(
                    "list_replace_null: unsupported element type {dt}"
                )));
            }
        };

        Ok(ColumnarValue::Array(result))
    }
}

pub fn list_replace_null_udf() -> ScalarUDF {
    ScalarUDF::from(ListReplaceNullUdf::new())
}

// ============================================================================
// list_make_array(List<T>, List<T>, ...) → List<List<T>>
//   Zip N lists element-wise: for each element i, produce [a[i], b[i], ...].
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListMakeArrayUdf {
    signature: Signature,
}

impl ListMakeArrayUdf {
    fn new() -> Self {
        Self {
            signature: Signature::variadic_any(Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListMakeArrayUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_make_array"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        // Input: List<T>, ..., Output: List<List<T>>
        match &arg_types[0] {
            DataType::List(inner) => {
                let inner_list_type = DataType::List(Arc::new(Field::new(
                    "item",
                    inner.data_type().clone(),
                    true,
                )));
                Ok(DataType::List(Arc::new(Field::new(
                    "item",
                    inner_list_type,
                    true,
                ))))
            }
            dt => Err(DataFusionError::Internal(format!(
                "list_make_array: expected List type, got {dt}"
            ))),
        }
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        if args.args.is_empty() {
            return Err(DataFusionError::Internal(
                "list_make_array requires at least 1 argument".into(),
            ));
        }

        let n = args.args.len();

        // Get all input ListArrays
        let arrays: Vec<ArrayRef> = args
            .args
            .iter()
            .map(|a| a.clone().into_array(1))
            .collect::<Result<Vec<_>>>()?;
        let lists: Vec<&ListArray> = arrays
            .iter()
            .map(|a| {
                a.as_any()
                    .downcast_ref::<ListArray>()
                    .ok_or_else(|| {
                        DataFusionError::Internal(
                            "list_make_array: all args must be ListArrays".into(),
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        let total_elements = lists[0].values().len();
        let element_type = lists[0].value_type();

        // Verify all lists have same structure
        for (i, list) in lists.iter().enumerate().skip(1) {
            if list.values().len() != total_elements {
                return Err(DataFusionError::Internal(format!(
                    "list_make_array: arg {i} has {} elements, expected {total_elements}",
                    list.values().len()
                )));
            }
        }

        let result: ArrayRef = match element_type {
            DataType::Int32 => {
                let flat_arrays: Vec<&Int32Array> = lists
                    .iter()
                    .map(|l| l.values().as_primitive::<Int32Type>())
                    .collect();

                // Interleave flat values
                let mut interleaved = Int32Builder::with_capacity(total_elements * n);
                for elem_idx in 0..total_elements {
                    for arr in &flat_arrays {
                        if arr.is_null(elem_idx) {
                            interleaved.append_null();
                        } else {
                            interleaved.append_value(arr.value(elem_idx));
                        }
                    }
                }
                let interleaved_arr = interleaved.finish();

                // Inner offsets: [0, n, 2n, 3n, ...]
                let inner_offsets: Vec<i32> =
                    (0..=total_elements).map(|i| (i * n) as i32).collect();
                let inner_offsets = OffsetBuffer::new(inner_offsets.into());

                let inner_list = ListArray::new(
                    i32_field(),
                    inner_offsets,
                    Arc::new(interleaved_arr),
                    None, // inner null bitmap inherited from element nulls
                );

                let outer_field = Arc::new(Field::new(
                    "item",
                    DataType::List(i32_field()),
                    true,
                ));
                Arc::new(ListArray::new(
                    outer_field,
                    lists[0].offsets().clone(),
                    Arc::new(inner_list),
                    lists[0].nulls().cloned(),
                ))
            }
            DataType::Float64 => {
                let flat_arrays: Vec<&Float64Array> = lists
                    .iter()
                    .map(|l| l.values().as_primitive::<Float64Type>())
                    .collect();

                let mut interleaved = Float64Builder::with_capacity(total_elements * n);
                for elem_idx in 0..total_elements {
                    for arr in &flat_arrays {
                        if arr.is_null(elem_idx) {
                            interleaved.append_null();
                        } else {
                            interleaved.append_value(arr.value(elem_idx));
                        }
                    }
                }
                let interleaved_arr = interleaved.finish();

                let inner_offsets: Vec<i32> =
                    (0..=total_elements).map(|i| (i * n) as i32).collect();
                let inner_offsets = OffsetBuffer::new(inner_offsets.into());

                let inner_list = ListArray::new(
                    f64_field(),
                    inner_offsets,
                    Arc::new(interleaved_arr),
                    None,
                );

                let outer_field = Arc::new(Field::new(
                    "item",
                    DataType::List(f64_field()),
                    true,
                ));
                Arc::new(ListArray::new(
                    outer_field,
                    lists[0].offsets().clone(),
                    Arc::new(inner_list),
                    lists[0].nulls().cloned(),
                ))
            }
            dt => {
                return Err(DataFusionError::Internal(format!(
                    "list_make_array: unsupported element type {dt}"
                )));
            }
        };

        Ok(ColumnarValue::Array(result))
    }
}

pub fn list_make_array_udf() -> ScalarUDF {
    ScalarUDF::from(ListMakeArrayUdf::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_args;
    use datafusion::arrow::array::{Array, AsArray, ListBuilder};
    use datafusion::common::ScalarValue;

    fn make_nested_i32_list(rows: &[Vec<Vec<Option<i32>>>]) -> ListArray {
        let mut outer = ListBuilder::new(ListBuilder::new(Int32Builder::new()));
        for row in rows {
            let inner_builder = outer.values();
            for inner_list in row {
                let val_builder = inner_builder.values();
                for v in inner_list {
                    match v {
                        Some(x) => val_builder.append_value(*x),
                        None => val_builder.append_null(),
                    }
                }
                inner_builder.append(true);
            }
            outer.append(true);
        }
        outer.finish()
    }

    fn make_i32_list(values: &[Option<i32>]) -> ListArray {
        let mut builder = ListBuilder::new(Int32Builder::new());
        for v in values {
            match v {
                Some(x) => builder.values().append_value(*x),
                None => builder.values().append_null(),
            }
        }
        builder.append(true);
        builder.finish()
    }

    #[test]
    fn test_list_index() {
        let nested = make_nested_i32_list(&[vec![
            vec![Some(10), Some(20), Some(30)],
            vec![Some(40), Some(50), Some(60)],
            vec![Some(70), Some(80), Some(90)],
        ]]);

        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(nested)),
                ColumnarValue::Scalar(ScalarValue::Int64(Some(2))),
            ],
            DataType::List(i32_field()),
        );

        let result = ListIndexUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let vals = rl.values().as_primitive::<Int32Type>();
            assert_eq!(vals.value(0), 20);
            assert_eq!(vals.value(1), 50);
            assert_eq!(vals.value(2), 80);
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_replace_null_i32() {
        let list = make_i32_list(&[Some(1), None, Some(3), None]);
        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(list)),
                ColumnarValue::Scalar(ScalarValue::Int32(Some(0))),
            ],
            DataType::List(i32_field()),
        );

        let result = ListReplaceNullUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let vals = rl.values().as_primitive::<Int32Type>();
            assert_eq!(vals.value(0), 1);
            assert_eq!(vals.value(1), 0);
            assert_eq!(vals.value(2), 3);
            assert_eq!(vals.value(3), 0);
            assert!(!vals.is_null(1));
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_make_array() {
        let a = make_i32_list(&[Some(1), Some(2), Some(3)]);
        let b = make_i32_list(&[Some(10), Some(20), Some(30)]);
        let c = make_i32_list(&[Some(100), Some(200), Some(300)]);

        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(a)),
                ColumnarValue::Array(Arc::new(b)),
                ColumnarValue::Array(Arc::new(c)),
            ],
            DataType::List(Arc::new(Field::new(
                "item",
                DataType::List(i32_field()),
                true,
            ))),
        );

        let result = ListMakeArrayUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let outer = arr.as_list::<i32>();
            assert_eq!(outer.len(), 1);

            let inner = outer.values().as_list::<i32>();
            assert_eq!(inner.len(), 3);

            let e0 = inner.value(0);
            let e0_ints = e0.as_primitive::<Int32Type>();
            assert_eq!(e0_ints.value(0), 1);
            assert_eq!(e0_ints.value(1), 10);
            assert_eq!(e0_ints.value(2), 100);

            let e1 = inner.value(1);
            let e1_ints = e1.as_primitive::<Int32Type>();
            assert_eq!(e1_ints.value(0), 2);
            assert_eq!(e1_ints.value(1), 20);
            assert_eq!(e1_ints.value(2), 200);
        } else {
            panic!("Expected Array");
        }
    }
}
