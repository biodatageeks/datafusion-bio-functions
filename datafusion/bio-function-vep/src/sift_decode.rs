//! Public DataFusion UDF that decodes position-sliced SIFT/PolyPhen `Binary`
//! blobs (v2 de-interleaved fixed-divisor layout) into structured rows, so
//! tools other than the VEP engine can read `translation_sift.cache`.

use std::any::Any;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, BinaryArray, Float32Builder, ListBuilder, StructBuilder, UInt8Builder,
};
use datafusion::arrow::datatypes::{DataType, Field, Fields};
use datafusion::common::cast::as_string_array;
use datafusion::common::{DataFusionError, Result, exec_err};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, TypeSignature,
    Volatility,
};

use crate::cache_common::{POLY_CODEC, SIFT_CODEC, deserialize_position_entries_v2};

fn item_fields() -> Fields {
    Fields::from(vec![
        Field::new("amino_acid", DataType::UInt8, false),
        Field::new("prediction", DataType::UInt8, false),
        Field::new("score", DataType::Float32, false),
    ])
}

fn return_item_field() -> Field {
    Field::new("item", DataType::Struct(item_fields()), true)
}

/// `vep_decode_sift_predictions(blob: Binary, predictor: Utf8)` decoding one
/// position-sliced SIFT/PolyPhen blob into a list of `{amino_acid, prediction,
/// score}` structs.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct SiftDecodeUdf {
    signature: Signature,
}

impl Default for SiftDecodeUdf {
    fn default() -> Self {
        Self {
            signature: Signature::new(
                TypeSignature::Exact(vec![DataType::Binary, DataType::Utf8]),
                Volatility::Immutable,
            ),
        }
    }
}

impl ScalarUDFImpl for SiftDecodeUdf {
    fn name(&self) -> &str {
        "vep_decode_sift_predictions"
    }

    fn signature(&self) -> &Signature {
        &self.signature
    }

    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(DataType::List(Arc::new(return_item_field())))
    }

    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arrays = ColumnarValue::values_to_arrays(&args.args)?;
        let blobs = arrays[0]
            .as_any()
            .downcast_ref::<BinaryArray>()
            .ok_or_else(|| {
                DataFusionError::Execution(
                    "vep_decode_sift_predictions: arg 0 must be Binary".into(),
                )
            })?;
        let predictor = as_string_array(&arrays[1])?;

        let values_builder = StructBuilder::new(
            item_fields(),
            vec![
                Box::new(UInt8Builder::new()),
                Box::new(UInt8Builder::new()),
                Box::new(Float32Builder::new()),
            ],
        );
        let mut list = ListBuilder::new(values_builder);

        for row in 0..blobs.len() {
            // A null blob or a null predictor (nullable column, `CAST(NULL AS
            // Utf8)`) yields a null list rather than panicking on `value(row)`.
            if blobs.is_null(row) || predictor.is_null(row) {
                list.append_null();
                continue;
            }
            let codec = match predictor.value(row) {
                "sift" => SIFT_CODEC,
                "polyphen" | "poly" => POLY_CODEC,
                other => return exec_err!("unknown predictor '{other}' (expected sift|polyphen)"),
            };
            // The position is not stored in the blob (it lives in the row key),
            // so we pass 0 as a placeholder. The decoded structs' `position`
            // field is therefore always 0 here, but it is intentionally not
            // surfaced: the output schema (`item_fields`) exposes only
            // amino_acid/prediction/score.
            let entries = deserialize_position_entries_v2(0, blobs.value(row), codec)?;
            let sb = list.values();
            for e in &entries {
                sb.field_builder::<UInt8Builder>(0)
                    .unwrap()
                    .append_value(e.amino_acid);
                sb.field_builder::<UInt8Builder>(1)
                    .unwrap()
                    .append_value(e.prediction);
                sb.field_builder::<Float32Builder>(2)
                    .unwrap()
                    .append_value(e.score);
                sb.append(true);
            }
            list.append(true);
        }
        Ok(ColumnarValue::Array(Arc::new(list.finish())))
    }
}

/// Build the `vep_decode_sift_predictions` scalar UDF.
pub fn vep_decode_sift_predictions_udf() -> ScalarUDF {
    ScalarUDF::from(SiftDecodeUdf::default())
}

/// Register the VEP SIFT decode UDF on a session context.
pub fn register_vep_sift_functions(ctx: &datafusion::prelude::SessionContext) {
    ctx.register_udf(vep_decode_sift_predictions_udf());
}

#[cfg(test)]
fn invoke_for_test(udf: &ScalarUDF, args: Vec<ColumnarValue>, number_rows: usize) -> ArrayRef {
    use datafusion::config::ConfigOptions;
    let arg_fields = args
        .iter()
        .map(|a| Arc::new(Field::new("a", a.data_type(), true)))
        .collect::<Vec<_>>();
    let return_field = Arc::new(Field::new(
        "out",
        udf.return_type(&[DataType::Binary, DataType::Utf8])
            .unwrap(),
        true,
    ));
    let res = udf
        .invoke_with_args(ScalarFunctionArgs {
            args,
            arg_fields,
            number_rows,
            return_field,
            config_options: Arc::new(ConfigOptions::default()),
        })
        .unwrap();
    match res {
        ColumnarValue::Array(a) => a,
        ColumnarValue::Scalar(s) => s.to_array().unwrap(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{
        Array, BinaryArray, Float32Array, ListArray, StringArray, StructArray, UInt8Array,
    };
    use datafusion::logical_expr::ColumnarValue;

    #[test]
    fn decodes_sift_blob_to_struct_list() {
        use crate::cache_common::{SIFT_CODEC, serialize_position_entries_v2};
        use crate::transcript_consequence::CompactPrediction;
        let blob = serialize_position_entries_v2(
            &[
                CompactPrediction {
                    position: 0,
                    amino_acid: 4,
                    prediction: 1,
                    score: 0.02,
                },
                CompactPrediction {
                    position: 0,
                    amino_acid: 9,
                    prediction: 0,
                    score: 1.0,
                },
            ],
            SIFT_CODEC,
        )
        .unwrap();
        let udf = vep_decode_sift_predictions_udf();
        let args = vec![
            ColumnarValue::Array(Arc::new(BinaryArray::from(vec![blob.as_slice()]))),
            ColumnarValue::Array(Arc::new(StringArray::from(vec!["sift"]))),
        ];
        let out = invoke_for_test(&udf, args, 1);
        let list = out.as_any().downcast_ref::<ListArray>().unwrap();
        let inner = list.value(0);
        let s = inner.as_any().downcast_ref::<StructArray>().unwrap();
        assert_eq!(s.len(), 2);
        let aa = s.column(0).as_any().downcast_ref::<UInt8Array>().unwrap();
        let score = s.column(2).as_any().downcast_ref::<Float32Array>().unwrap();
        assert_eq!(aa.value(0), 4);
        assert_eq!(score.value(1).to_bits(), 1.0f32.to_bits());
    }

    #[test]
    fn null_predictor_yields_null_without_panicking() {
        use crate::cache_common::{SIFT_CODEC, serialize_position_entries_v2};
        use crate::transcript_consequence::CompactPrediction;
        let blob = serialize_position_entries_v2(
            &[CompactPrediction {
                position: 0,
                amino_acid: 4,
                prediction: 1,
                score: 0.02,
            }],
            SIFT_CODEC,
        )
        .unwrap();
        let udf = vep_decode_sift_predictions_udf();
        // Row 0: valid blob + valid predictor. Row 1: valid blob + NULL predictor.
        let out = invoke_for_test(
            &udf,
            vec![
                ColumnarValue::Array(Arc::new(BinaryArray::from(vec![
                    Some(blob.as_slice()),
                    Some(blob.as_slice()),
                ]))),
                ColumnarValue::Array(Arc::new(StringArray::from(vec![Some("sift"), None]))),
            ],
            2,
        );
        let list = out.as_any().downcast_ref::<ListArray>().unwrap();
        assert!(!list.is_null(0), "valid predictor decodes");
        assert!(list.is_null(1), "null predictor yields null, not a panic");
    }

    #[test]
    fn udf_matches_core_decode_polyphen() {
        use crate::cache_common::{
            POLY_CODEC, deserialize_position_entries_v2, serialize_position_entries_v2,
        };
        use crate::transcript_consequence::CompactPrediction;

        let entries: Vec<CompactPrediction> = (0..19u8)
            .map(|i| CompactPrediction {
                position: 0,
                amino_acid: i,
                prediction: i % 3,
                // On the k/1000 grid (PolyPhen resolution), so the lossless codec accepts it.
                score: (i as f32) / 1000.0,
            })
            .collect();
        let blob = serialize_position_entries_v2(&entries, POLY_CODEC).unwrap();
        let expected = deserialize_position_entries_v2(0, &blob, POLY_CODEC).unwrap();

        let udf = vep_decode_sift_predictions_udf();
        let out = invoke_for_test(
            &udf,
            vec![
                ColumnarValue::Array(Arc::new(BinaryArray::from(vec![blob.as_slice()]))),
                ColumnarValue::Array(Arc::new(StringArray::from(vec!["polyphen"]))),
            ],
            1,
        );
        let list = out.as_any().downcast_ref::<ListArray>().unwrap();
        let s = list.value(0);
        let s = s.as_any().downcast_ref::<StructArray>().unwrap();
        let score = s.column(2).as_any().downcast_ref::<Float32Array>().unwrap();
        for (i, e) in expected.iter().enumerate() {
            assert_eq!(score.value(i).to_bits(), e.score.to_bits());
        }
    }
}
