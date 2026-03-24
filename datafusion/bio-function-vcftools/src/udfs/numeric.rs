//! Element-wise arithmetic, math, and bounds UDFs for List columns.

use std::any::Any;
use std::sync::Arc;

use datafusion::arrow::array::{AsArray, Float64Array, ListArray};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::Result;
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, Volatility,
};

use super::{f64_list_type, list_values_as_f64, rebuild_f64_list, scalar_to_f64};

// ============================================================================
// Binary arithmetic: list_add, list_sub, list_mul, list_div
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum ArithOp {
    Add,
    Sub,
    Mul,
    Div,
}

impl ArithOp {
    fn apply(self, a: f64, b: f64) -> f64 {
        match self {
            Self::Add => a + b,
            Self::Sub => a - b,
            Self::Mul => a * b,
            Self::Div => a / b,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Add => "list_add",
            Self::Sub => "list_sub",
            Self::Mul => "list_mul",
            Self::Div => "list_div",
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListBinaryArithUdf {
    op: ArithOp,
    signature: Signature,
}

impl ListBinaryArithUdf {
    fn new(op: ArithOp) -> Self {
        Self {
            op,
            signature: Signature::any(2, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListBinaryArithUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        self.op.name()
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(f64_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let op = self.op;

        // Detect which arg is the list and which might be scalar
        let lhs_is_list = matches!(&args.args[0], ColumnarValue::Array(a) if a.as_any().downcast_ref::<ListArray>().is_some());
        let rhs_is_list = matches!(&args.args[1], ColumnarValue::Array(a) if a.as_any().downcast_ref::<ListArray>().is_some());

        if lhs_is_list {
            let lhs_arr = args.args[0].clone().into_array(1)?;
            let lhs_list = lhs_arr.as_list::<i32>();
            let lhs_vals = list_values_as_f64(lhs_list)?;

            let result: Float64Array = if rhs_is_list {
                let rhs_arr = args.args[1].clone().into_array(1)?;
                let rhs_list = rhs_arr.as_list::<i32>();
                let rhs_vals = list_values_as_f64(rhs_list)?;
                lhs_vals.iter().zip(rhs_vals.iter())
                    .map(|(a, b)| match (a, b) {
                        (Some(a), Some(b)) => Some(op.apply(a, b)),
                        _ => None,
                    })
                    .collect()
            } else {
                let scalar = scalar_to_f64(&args.args[1])?;
                lhs_vals.iter().map(|v| v.map(|x| op.apply(x, scalar))).collect()
            };
            Ok(ColumnarValue::Array(Arc::new(rebuild_f64_list(lhs_list, result))))
        } else if rhs_is_list {
            // Scalar op List → broadcast scalar on left: scalar op list[i]
            let scalar = scalar_to_f64(&args.args[0])?;
            let rhs_arr = args.args[1].clone().into_array(1)?;
            let rhs_list = rhs_arr.as_list::<i32>();
            let rhs_vals = list_values_as_f64(rhs_list)?;
            let result: Float64Array = rhs_vals.iter()
                .map(|v| v.map(|x| op.apply(scalar, x)))
                .collect();
            Ok(ColumnarValue::Array(Arc::new(rebuild_f64_list(rhs_list, result))))
        } else {
            Err(datafusion::common::DataFusionError::Internal(
                format!("{}: at least one argument must be a List", op.name()),
            ))
        }
    }
}

pub fn list_add_udf() -> ScalarUDF {
    ScalarUDF::from(ListBinaryArithUdf::new(ArithOp::Add))
}
pub fn list_sub_udf() -> ScalarUDF {
    ScalarUDF::from(ListBinaryArithUdf::new(ArithOp::Sub))
}
pub fn list_mul_udf() -> ScalarUDF {
    ScalarUDF::from(ListBinaryArithUdf::new(ArithOp::Mul))
}
pub fn list_div_udf() -> ScalarUDF {
    ScalarUDF::from(ListBinaryArithUdf::new(ArithOp::Div))
}

// ============================================================================
// Unary: list_neg
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListNegUdf {
    signature: Signature,
}

impl ListNegUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(1, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListNegUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_neg"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(f64_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arr = args.args[0].clone().into_array(1)?;
        let list = arr.as_list::<i32>();
        let vals = list_values_as_f64(list)?;
        let result: Float64Array = vals.iter().map(|v| v.map(|x| -x)).collect();
        Ok(ColumnarValue::Array(Arc::new(rebuild_f64_list(
            list, result,
        ))))
    }
}

pub fn list_neg_udf() -> ScalarUDF {
    ScalarUDF::from(ListNegUdf::new())
}

// ============================================================================
// Unary math: list_log10, list_sqrt, list_pow10, list_round
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum MathOp {
    Log10,
    Sqrt,
    Pow10,
    Round,
}

impl MathOp {
    fn apply(self, x: f64) -> Option<f64> {
        match self {
            Self::Log10 => Some(x.log10()),
            Self::Sqrt => Some(x.sqrt()),
            Self::Pow10 => Some(10.0_f64.powf(x)),
            Self::Round => Some(x.round()),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Log10 => "list_log10",
            Self::Sqrt => "list_sqrt",
            Self::Pow10 => "list_pow10",
            Self::Round => "list_round",
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListUnaryMathUdf {
    op: MathOp,
    signature: Signature,
}

impl ListUnaryMathUdf {
    fn new(op: MathOp) -> Self {
        Self {
            op,
            signature: Signature::any(1, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListUnaryMathUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        self.op.name()
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(f64_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arr = args.args[0].clone().into_array(1)?;
        let list = arr.as_list::<i32>();
        let vals = list_values_as_f64(list)?;
        let op = self.op;
        let result: Float64Array = vals.iter().map(|v| v.and_then(|x| op.apply(x))).collect();
        Ok(ColumnarValue::Array(Arc::new(rebuild_f64_list(
            list, result,
        ))))
    }
}

pub fn list_log10_udf() -> ScalarUDF {
    ScalarUDF::from(ListUnaryMathUdf::new(MathOp::Log10))
}
pub fn list_sqrt_udf() -> ScalarUDF {
    ScalarUDF::from(ListUnaryMathUdf::new(MathOp::Sqrt))
}
pub fn list_pow10_udf() -> ScalarUDF {
    ScalarUDF::from(ListUnaryMathUdf::new(MathOp::Pow10))
}
pub fn list_round_udf() -> ScalarUDF {
    ScalarUDF::from(ListUnaryMathUdf::new(MathOp::Round))
}

// ============================================================================
// Bounds: list_greatest, list_least, list_clamp
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum BoundOp {
    Greatest,
    Least,
}

impl BoundOp {
    fn name(self) -> &'static str {
        match self {
            Self::Greatest => "list_greatest",
            Self::Least => "list_least",
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListBoundUdf {
    op: BoundOp,
    signature: Signature,
}

impl ListBoundUdf {
    fn new(op: BoundOp) -> Self {
        Self {
            op,
            signature: Signature::any(2, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListBoundUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        self.op.name()
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(f64_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arr = args.args[0].clone().into_array(1)?;
        let list = arr.as_list::<i32>();
        let vals = list_values_as_f64(list)?;
        let bound = scalar_to_f64(&args.args[1])?;

        let result: Float64Array = match self.op {
            BoundOp::Greatest => vals.iter().map(|v| v.map(|x| x.max(bound))).collect(),
            BoundOp::Least => vals.iter().map(|v| v.map(|x| x.min(bound))).collect(),
        };

        Ok(ColumnarValue::Array(Arc::new(rebuild_f64_list(
            list, result,
        ))))
    }
}

pub fn list_greatest_udf() -> ScalarUDF {
    ScalarUDF::from(ListBoundUdf::new(BoundOp::Greatest))
}
pub fn list_least_udf() -> ScalarUDF {
    ScalarUDF::from(ListBoundUdf::new(BoundOp::Least))
}

// ============================================================================
// list_clamp(List<T>, min, max) → List<Float64>
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct ListClampUdf {
    signature: Signature,
}

impl ListClampUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(3, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for ListClampUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "list_clamp"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(f64_list_type())
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arr = args.args[0].clone().into_array(1)?;
        let list = arr.as_list::<i32>();
        let vals = list_values_as_f64(list)?;
        let min_val = scalar_to_f64(&args.args[1])?;
        let max_val = scalar_to_f64(&args.args[2])?;

        let result: Float64Array = vals
            .iter()
            .map(|v| v.map(|x| x.max(min_val).min(max_val)))
            .collect();

        Ok(ColumnarValue::Array(Arc::new(rebuild_f64_list(
            list, result,
        ))))
    }
}

pub fn list_clamp_udf() -> ScalarUDF {
    ScalarUDF::from(ListClampUdf::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_args;
    use datafusion::arrow::array::{Array, Float64Builder, Int32Builder, ListBuilder};
    use datafusion::arrow::datatypes::Float64Type;
    use datafusion::common::ScalarValue;

    fn make_f64_list(values: &[Option<f64>]) -> ListArray {
        let mut builder = ListBuilder::new(Float64Builder::new());
        for v in values {
            match v {
                Some(x) => builder.values().append_value(*x),
                None => builder.values().append_null(),
            }
        }
        builder.append(true);
        builder.finish()
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
    fn test_list_add_list_scalar() {
        let list = make_i32_list(&[Some(10), Some(20), None, Some(40)]);
        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(list)),
                ColumnarValue::Scalar(ScalarValue::Float64(Some(5.0))),
            ],
            f64_list_type(),
        );
        let result = ListBinaryArithUdf::new(ArithOp::Add)
            .invoke_with_args(args)
            .unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let vals = rl.values().as_primitive::<Float64Type>();
            assert_eq!(vals.value(0), 15.0);
            assert_eq!(vals.value(1), 25.0);
            assert!(vals.is_null(2));
            assert_eq!(vals.value(3), 45.0);
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_mul_two_lists() {
        let a = make_f64_list(&[Some(2.0), Some(3.0), Some(4.0)]);
        let b = make_f64_list(&[Some(10.0), Some(20.0), Some(30.0)]);
        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(a)),
                ColumnarValue::Array(Arc::new(b)),
            ],
            f64_list_type(),
        );
        let result = ListBinaryArithUdf::new(ArithOp::Mul)
            .invoke_with_args(args)
            .unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let vals = rl.values().as_primitive::<Float64Type>();
            assert_eq!(vals.value(0), 20.0);
            assert_eq!(vals.value(1), 60.0);
            assert_eq!(vals.value(2), 120.0);
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_log10() {
        let list = make_f64_list(&[Some(1.0), Some(10.0), Some(100.0), Some(0.001)]);
        let args = test_args(
            vec![ColumnarValue::Array(Arc::new(list))],
            f64_list_type(),
        );
        let result = ListUnaryMathUdf::new(MathOp::Log10)
            .invoke_with_args(args)
            .unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let vals = rl.values().as_primitive::<Float64Type>();
            assert!((vals.value(0) - 0.0).abs() < 1e-10);
            assert!((vals.value(1) - 1.0).abs() < 1e-10);
            assert!((vals.value(2) - 2.0).abs() < 1e-10);
            assert!((vals.value(3) - (-3.0)).abs() < 1e-10);
        } else {
            panic!("Expected Array");
        }
    }

    #[test]
    fn test_list_clamp() {
        let list = make_f64_list(&[Some(-5.0), Some(50.0), Some(300.0), None]);
        let args = test_args(
            vec![
                ColumnarValue::Array(Arc::new(list)),
                ColumnarValue::Scalar(ScalarValue::Float64(Some(0.0))),
                ColumnarValue::Scalar(ScalarValue::Float64(Some(255.0))),
            ],
            f64_list_type(),
        );
        let result = ListClampUdf::new().invoke_with_args(args).unwrap();
        if let ColumnarValue::Array(arr) = result {
            let rl = arr.as_list::<i32>();
            let vals = rl.values().as_primitive::<Float64Type>();
            assert_eq!(vals.value(0), 0.0);
            assert_eq!(vals.value(1), 50.0);
            assert_eq!(vals.value(2), 255.0);
            assert!(vals.is_null(3));
        } else {
            panic!("Expected Array");
        }
    }
}
