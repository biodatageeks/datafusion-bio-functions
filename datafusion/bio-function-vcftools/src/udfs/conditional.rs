//! Element-wise conditional UDFs: list_where, list_case.

use std::any::Any;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, AsArray, BooleanArray, Float64Array, Int32Array, ListArray, StringArray,
};
use datafusion::arrow::datatypes::{DataType, Field, Float64Type, Int32Type};
use datafusion::common::{DataFusionError, Result};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, Volatility,
};

use super::{f64_field, i32_field, scalar_to_utf8, utf8_field};

// ============================================================================
// list_where(List<Bool>, List<T>, List<T>) → List<T>
//   Element-wise: mask[i] ? true_val[i] : false_val[i]
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListWhereUdf {
    signature: Signature,
}

impl ListWhereUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(3, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListWhereUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_where"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        // Output type matches the true branch (arg[1])
        Ok(arg_types[1].clone())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let mask_arr = args.args[0].clone().into_array(1)?;
        let true_arr = args.args[1].clone().into_array(1)?;
        let false_arr = args.args[2].clone().into_array(1)?;

        let mask_list = mask_arr.as_list::<i32>();
        let true_list = true_arr.as_list::<i32>();
        let false_list = false_arr.as_list::<i32>();

        let mask_flat = mask_list
            .values()
            .as_any()
            .downcast_ref::<BooleanArray>()
            .ok_or_else(|| {
                DataFusionError::Internal("list_where: first arg must be List<Bool>".into())
            })?;

        let element_type = true_list.value_type();

        let result: ListArray = match element_type {
            DataType::Int32 => {
                let t = true_list.values().as_primitive::<Int32Type>();
                let f = false_list.values().as_primitive::<Int32Type>();
                let out: Int32Array = (0..mask_flat.len())
                    .map(|i| {
                        if mask_flat.is_null(i) {
                            return None;
                        }
                        let src = if mask_flat.value(i) { t } else { f };
                        if src.is_null(i) {
                            None
                        } else {
                            Some(src.value(i))
                        }
                    })
                    .collect();
                ListArray::new(
                    i32_field(),
                    true_list.offsets().clone(),
                    Arc::new(out),
                    true_list.nulls().cloned(),
                )
            }
            DataType::Float64 => {
                let t = true_list.values().as_primitive::<Float64Type>();
                let f = false_list.values().as_primitive::<Float64Type>();
                let out: Float64Array = (0..mask_flat.len())
                    .map(|i| {
                        if mask_flat.is_null(i) {
                            return None;
                        }
                        let src = if mask_flat.value(i) { t } else { f };
                        if src.is_null(i) {
                            None
                        } else {
                            Some(src.value(i))
                        }
                    })
                    .collect();
                ListArray::new(
                    f64_field(),
                    true_list.offsets().clone(),
                    Arc::new(out),
                    true_list.nulls().cloned(),
                )
            }
            DataType::Utf8 => {
                let t = true_list
                    .values()
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .unwrap();
                let f = false_list
                    .values()
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .unwrap();
                let out: StringArray = (0..mask_flat.len())
                    .map(|i| {
                        if mask_flat.is_null(i) {
                            return None;
                        }
                        if mask_flat.value(i) {
                            if t.is_null(i) {
                                None
                            } else {
                                Some(t.value(i).to_string())
                            }
                        } else if f.is_null(i) {
                            None
                        } else {
                            Some(f.value(i).to_string())
                        }
                    })
                    .collect();
                ListArray::new(
                    utf8_field(),
                    true_list.offsets().clone(),
                    Arc::new(out),
                    true_list.nulls().cloned(),
                )
            }
            dt => {
                return Err(DataFusionError::Internal(format!(
                    "list_where: unsupported element type {dt}"
                )));
            }
        };

        Ok(ColumnarValue::Array(Arc::new(result)))
    }
}

pub fn list_where_udf() -> ScalarUDF {
    ScalarUDF::from(ListWhereUdf::new())
}

// ============================================================================
// list_case(List<Utf8>, when1, then1, when2, then2, ..., [else])
//   Element-wise CASE: first match wins, else default or original.
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListCaseUdf {
    signature: Signature,
}

impl ListCaseUdf {
    fn new() -> Self {
        Self {
            signature: Signature::variadic_any(Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListCaseUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_case"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(DataType::List(Arc::new(Field::new(
            "item",
            DataType::Utf8,
            true,
        ))))
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        if args.args.len() < 3 {
            return Err(DataFusionError::Internal(
                "list_case requires at least 3 arguments: list, when, then".into(),
            ));
        }

        let list_arr = args.args[0].clone().into_array(1)?;
        let list = list_arr.as_list::<i32>();
        let flat = list
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| {
                DataFusionError::Internal("list_case: expected List<Utf8>".into())
            })?;

        // Parse when/then pairs and optional else from remaining args
        let rest = &args.args[1..];
        let num_pairs = rest.len() / 2;
        let has_else = rest.len() % 2 == 1;

        let mut pairs: Vec<(String, String)> = Vec::with_capacity(num_pairs);
        for i in 0..num_pairs {
            let when_val = scalar_to_utf8(&rest[i * 2])?;
            let then_val = scalar_to_utf8(&rest[i * 2 + 1])?;
            pairs.push((when_val, then_val));
        }
        let else_val = if has_else {
            Some(scalar_to_utf8(rest.last().unwrap())?)
        } else {
            None
        };

        let result: StringArray = (0..flat.len())
            .map(|i| {
                if flat.is_null(i) {
                    return None;
                }
                let val = flat.value(i);
                for (when, then) in &pairs {
                    if val == when.as_str() {
                        return Some(then.clone());
                    }
                }
                // No match: use else or keep original
                Some(
                    else_val
                        .as_ref()
                        .cloned()
                        .unwrap_or_else(|| val.to_string()),
                )
            })
            .collect();

        Ok(ColumnarValue::Array(Arc::new(ListArray::new(
            utf8_field(),
            list.offsets().clone(),
            Arc::new(result),
            list.nulls().cloned(),
        ))))
    }
}

pub fn list_case_udf() -> ScalarUDF {
    ScalarUDF::from(ListCaseUdf::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_args;
    use datafusion::arrow::array::{BooleanBuilder, Int32Builder, ListBuilder, StringBuilder};
    use datafusion::common::ScalarValue;

    #[test]
    fn test_list_where_int32() {
        let mut mask_b = ListBuilder::new(BooleanBuilder::new());
        mask_b.values().append_value(true);
        mask_b.values().append_value(false);
        mask_b.values().append_value(true);
        mask_b.append(true);

        let mut tb = ListBuilder::new(Int32Builder::new());
        tb.values().append_value(10);
        tb.values().append_value(20);
        tb.values().append_value(30);
        tb.append(true);

        let mut fb = ListBuilder::new(Int32Builder::new());
        fb.values().append_value(100);
        fb.values().append_value(200);
        fb.values().append_value(300);
        fb.append(true);

        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(mask_b.finish())),
                ColumnarValue::Array(Arc::new(tb.finish())),
                ColumnarValue::Array(Arc::new(fb.finish())),
            ],
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
        );

        let result = ListWhereUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let ints = rl.values().as_primitive::<Int32Type>();
            assert_eq!(ints.value(0), 10);
            assert_eq!(ints.value(1), 200);
            assert_eq!(ints.value(2), 30);
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_case() {
        let mut builder = ListBuilder::new(StringBuilder::new());
        builder.values().append_value("0");
        builder.values().append_value("1");
        builder.values().append_value("0/1");
        builder.values().append_value("2");
        builder.append(true);

        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(builder.finish())),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("0".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("0/0".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("1".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("1/1".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("2".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("2/2".into()))),
            ],
            DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
        );

        let result = ListCaseUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let strs = rl.values().as_any().downcast_ref::<StringArray>().unwrap();
            assert_eq!(strs.value(0), "0/0");
            assert_eq!(strs.value(1), "1/1");
            assert_eq!(strs.value(2), "0/1");
            assert_eq!(strs.value(3), "2/2");
        } else {
            panic!("Expected Array");
        }
    }
}
