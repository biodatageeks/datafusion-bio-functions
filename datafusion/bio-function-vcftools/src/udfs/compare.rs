//! Element-wise comparison UDFs: list_eq, list_in.

use std::any::Any;
use std::collections::HashSet;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, AsArray, BooleanArray, ListArray, StringArray,
};
use datafusion::arrow::datatypes::{DataType, Float64Type, Int32Type};
use datafusion::common::{DataFusionError, Result};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, Volatility,
};

use super::{bool_field, bool_list_type, scalar_to_f64, scalar_to_utf8};

// ============================================================================
// list_eq(List<T>, T) → List<Bool>  — element-wise equality
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListEqUdf {
    signature: Signature,
}

impl ListEqUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(2, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListEqUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_eq"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(bool_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let list_arr = args.args[0].clone().into_array(1)?;
        let list = list_arr.as_list::<i32>();
        let flat = list.values();

        let result: BooleanArray = match flat.data_type() {
            DataType::Int32 => {
                let scalar = scalar_to_f64(&args.args[1])? as i32;
                let arr = flat.as_primitive::<Int32Type>();
                arr.iter()
                    .map(|v| v.map(|x| x == scalar))
                    .collect()
            }
            DataType::Float64 => {
                let scalar = scalar_to_f64(&args.args[1])?;
                let arr = flat.as_primitive::<Float64Type>();
                arr.iter()
                    .map(|v| v.map(|x| (x - scalar).abs() < f64::EPSILON))
                    .collect()
            }
            DataType::Utf8 => {
                let scalar = scalar_to_utf8(&args.args[1])?;
                let arr = flat
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .ok_or_else(|| DataFusionError::Internal("Expected StringArray".into()))?;
                (0..arr.len())
                    .map(|i| {
                        if arr.is_null(i) {
                            None
                        } else {
                            Some(arr.value(i) == scalar.as_str())
                        }
                    })
                    .collect()
            }
            dt => {
                return Err(DataFusionError::Internal(format!(
                    "list_eq: unsupported element type {dt}"
                )));
            }
        };

        Ok(ColumnarValue::Array(Arc::new(ListArray::new(
            bool_field(),
            list.offsets().clone(),
            Arc::new(result),
            list.nulls().cloned(),
        ))))
    }
}

pub fn list_eq_udf() -> ScalarUDF {
    ScalarUDF::from(ListEqUdf::new())
}

// ============================================================================
// list_in(List<Utf8>, Utf8, Utf8, ...) → List<Bool>  — element-wise IN
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListInUdf {
    signature: Signature,
}

impl ListInUdf {
    fn new() -> Self {
        Self {
            signature: Signature::variadic_any(Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListInUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_in"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(bool_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        if args.args.len() < 2 {
            return Err(DataFusionError::Internal(
                "list_in requires at least 2 arguments".into(),
            ));
        }

        let list_arr = args.args[0].clone().into_array(1)?;
        let list = list_arr.as_list::<i32>();

        // Collect match values into a set
        let mut match_set = HashSet::new();
        for arg in &args.args[1..] {
            match_set.insert(scalar_to_utf8(arg)?);
        }

        let flat = list
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| DataFusionError::Internal("list_in: expected List<Utf8>".into()))?;

        let result: BooleanArray = (0..flat.len())
            .map(|i| {
                if flat.is_null(i) {
                    None
                } else {
                    Some(match_set.contains(flat.value(i)))
                }
            })
            .collect();

        Ok(ColumnarValue::Array(Arc::new(ListArray::new(
            bool_field(),
            list.offsets().clone(),
            Arc::new(result),
            list.nulls().cloned(),
        ))))
    }
}

pub fn list_in_udf() -> ScalarUDF {
    ScalarUDF::from(ListInUdf::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_args;
    use datafusion::arrow::array::{Int32Builder, ListBuilder, StringBuilder};
    use datafusion::common::ScalarValue;

    #[test]
    fn test_list_eq_int32() {
        let mut builder = ListBuilder::new(Int32Builder::new());
        builder.values().append_value(1);
        builder.values().append_value(2);
        builder.values().append_value(1);
        builder.append(true);
        let list = builder.finish();

        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(list)),
                ColumnarValue::Scalar(ScalarValue::Int32(Some(1))),
            ],
            bool_list_type(),
        );

        let result = ListEqUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let result_list = arr.as_list::<i32>();
            let bools = result_list.values().as_any().downcast_ref::<BooleanArray>().unwrap();
            assert!(bools.value(0));
            assert!(!bools.value(1));
            assert!(bools.value(2));
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_in_utf8() {
        let mut builder = ListBuilder::new(StringBuilder::new());
        builder.values().append_value("0/0");
        builder.values().append_value("0/1");
        builder.values().append_value("0|0");
        builder.values().append_value("1/1");
        builder.append(true);
        let list = builder.finish();

        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(list)),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("0/0".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("0|0".into()))),
                ColumnarValue::Scalar(ScalarValue::Utf8(Some("0".into()))),
            ],
            bool_list_type(),
        );

        let result = ListInUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let result_list = arr.as_list::<i32>();
            let bools = result_list.values().as_any().downcast_ref::<BooleanArray>().unwrap();
            assert!(bools.value(0));   // 0/0 in set
            assert!(!bools.value(1));  // 0/1 not in set
            assert!(bools.value(2));   // 0|0 in set
            assert!(!bools.value(3));  // 1/1 not in set
        } else {
            panic!("Expected Array");
        }
    }
}
