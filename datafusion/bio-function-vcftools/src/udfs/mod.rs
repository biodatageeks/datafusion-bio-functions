//! Generic element-wise list UDFs for Apache DataFusion.
//!
//! Operate directly on `List<T>` columns without unnesting,
//! avoiding the unnest → transform → array_agg pattern entirely.

pub mod array_ops;
pub mod compare;
pub mod conditional;
pub mod numeric;

use std::sync::Arc;

use datafusion::arrow::array::{Array, AsArray, Float64Array, ListArray};
use datafusion::arrow::datatypes::{DataType, Field, Float32Type, Float64Type, Int32Type};
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::logical_expr::ColumnarValue;
use datafusion::prelude::SessionContext;

// Re-exports
pub use array_ops::{list_index_udf, list_make_array_udf, list_replace_null_udf};
pub use compare::{list_eq_udf, list_in_udf};
pub use conditional::{list_case_udf, list_where_udf};
pub use numeric::{
    list_add_udf, list_clamp_udf, list_div_udf, list_greatest_udf, list_least_udf,
    list_log10_udf, list_mul_udf, list_neg_udf, list_pow10_udf, list_round_udf, list_sqrt_udf,
    list_sub_udf,
};

/// Register all generic list UDFs with the given session context.
pub fn register_list_udfs(ctx: &SessionContext) {
    // Comparison
    ctx.register_udf(list_eq_udf());
    ctx.register_udf(list_in_udf());
    // Conditional
    ctx.register_udf(list_where_udf());
    ctx.register_udf(list_case_udf());
    // Arithmetic
    ctx.register_udf(list_add_udf());
    ctx.register_udf(list_sub_udf());
    ctx.register_udf(list_mul_udf());
    ctx.register_udf(list_div_udf());
    ctx.register_udf(list_neg_udf());
    // Math
    ctx.register_udf(list_log10_udf());
    ctx.register_udf(list_sqrt_udf());
    ctx.register_udf(list_pow10_udf());
    ctx.register_udf(list_round_udf());
    // Bounds
    ctx.register_udf(list_greatest_udf());
    ctx.register_udf(list_least_udf());
    ctx.register_udf(list_clamp_udf());
    // Array ops
    ctx.register_udf(list_index_udf());
    ctx.register_udf(list_replace_null_udf());
    ctx.register_udf(list_make_array_udf());
}

// ============================================================================
// Shared helpers
// ============================================================================

/// List<Float64> DataType.
pub(crate) fn f64_list_type() -> DataType {
    DataType::List(f64_field())
}

/// List<Boolean> DataType.
pub(crate) fn bool_list_type() -> DataType {
    DataType::List(bool_field())
}

pub(crate) fn f64_field() -> Arc<Field> {
    Arc::new(Field::new("item", DataType::Float64, true))
}

pub(crate) fn bool_field() -> Arc<Field> {
    Arc::new(Field::new("item", DataType::Boolean, true))
}

pub(crate) fn i32_field() -> Arc<Field> {
    Arc::new(Field::new("item", DataType::Int32, true))
}

pub(crate) fn utf8_field() -> Arc<Field> {
    Arc::new(Field::new("item", DataType::Utf8, true))
}

/// Convert a ListArray's flat values to Float64Array, casting Int32/Float32 if needed.
pub(crate) fn list_values_as_f64(list: &ListArray) -> Result<Float64Array> {
    let values = list.values();
    match values.data_type() {
        DataType::Float64 => Ok(values.as_primitive::<Float64Type>().clone()),
        DataType::Float32 => {
            let arr = values.as_primitive::<Float32Type>();
            Ok(arr.iter().map(|v| v.map(|x| x as f64)).collect())
        }
        DataType::Int32 => {
            let arr = values.as_primitive::<Int32Type>();
            Ok(arr.iter().map(|v| v.map(|x| x as f64)).collect())
        }
        dt => Err(DataFusionError::Internal(format!(
            "Cannot convert list element type {dt} to Float64"
        ))),
    }
}

/// Rebuild a ListArray<Float64> reusing offsets and null bitmap from the source.
pub(crate) fn rebuild_f64_list(src: &ListArray, values: Float64Array) -> ListArray {
    ListArray::new(
        f64_field(),
        src.offsets().clone(),
        Arc::new(values),
        src.nulls().cloned(),
    )
}

/// Extract a scalar f64 from a ColumnarValue.
pub(crate) fn scalar_to_f64(val: &ColumnarValue) -> Result<f64> {
    match val {
        ColumnarValue::Scalar(sv) => match sv {
            ScalarValue::Float64(Some(v)) => Ok(*v),
            ScalarValue::Float32(Some(v)) => Ok(*v as f64),
            ScalarValue::Int64(Some(v)) => Ok(*v as f64),
            ScalarValue::Int32(Some(v)) => Ok(*v as f64),
            ScalarValue::Int16(Some(v)) => Ok(*v as f64),
            ScalarValue::UInt32(Some(v)) => Ok(*v as f64),
            _ => Err(DataFusionError::Internal(format!(
                "Cannot convert scalar {sv} to f64"
            ))),
        },
        ColumnarValue::Array(a) if a.len() == 1 => {
            if let Some(arr) = a.as_any().downcast_ref::<Float64Array>() {
                if !arr.is_null(0) {
                    return Ok(arr.value(0));
                }
            }
            Err(DataFusionError::Internal(
                "Cannot extract scalar f64".into(),
            ))
        }
        _ => Err(DataFusionError::Internal("Expected scalar f64".into())),
    }
}

/// Extract a scalar string from a ColumnarValue.
pub(crate) fn scalar_to_utf8(val: &ColumnarValue) -> Result<String> {
    match val {
        ColumnarValue::Scalar(ScalarValue::Utf8(Some(s))) => Ok(s.clone()),
        ColumnarValue::Scalar(ScalarValue::LargeUtf8(Some(s))) => Ok(s.clone()),
        _ => Err(DataFusionError::Internal("Expected scalar Utf8".into())),
    }
}

/// Build a ScalarFunctionArgs for tests.
#[cfg(test)]
pub(crate) fn test_args(args: Vec<ColumnarValue>, return_type: DataType) -> datafusion::logical_expr::ScalarFunctionArgs {
    use datafusion::config::ConfigOptions;
    datafusion::logical_expr::ScalarFunctionArgs {
        args,
        arg_fields: vec![], // not inspected by our UDFs
        number_rows: 1,
        return_field: Arc::new(Field::new("output", return_type, true)),
        config_options: Arc::new(ConfigOptions::new()),
    }
}

/// Extract a scalar i64 from a ColumnarValue.
pub(crate) fn scalar_to_i64(val: &ColumnarValue) -> Result<i64> {
    match val {
        ColumnarValue::Scalar(sv) => match sv {
            ScalarValue::Int64(Some(v)) => Ok(*v),
            ScalarValue::Int32(Some(v)) => Ok(*v as i64),
            ScalarValue::Int16(Some(v)) => Ok(*v as i64),
            ScalarValue::Int8(Some(v)) => Ok(*v as i64),
            ScalarValue::UInt64(Some(v)) => Ok(*v as i64),
            ScalarValue::UInt32(Some(v)) => Ok(*v as i64),
            _ => Err(DataFusionError::Internal(format!(
                "Cannot convert scalar {sv} to i64"
            ))),
        },
        _ => Err(DataFusionError::Internal("Expected scalar i64".into())),
    }
}
