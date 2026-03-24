//! Fused per-sample genotype QC UDF: GT masking, PL correction, DS computation
//! in a single pass with zero intermediate List columns.

use std::any::Any;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, AsArray, Float32Builder, Int32Builder, ListArray, StringBuilder, StructArray,
};
use datafusion::arrow::buffer::OffsetBuffer;
use datafusion::arrow::datatypes::{DataType, Field, Float32Type, Int32Type};
use datafusion::common::{DataFusionError, Result};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, Volatility,
};

use super::{scalar_to_f64, scalar_to_utf8};

// ============================================================================
// vcf_genotype_qc(genotypes, gq_min, dp_min, dp_max, mask_gt)
//
// Single-pass per-sample genotype QC:
//   1. GT normalization (haploid → diploid)
//   2. Quality masking (GQ/DP thresholds)
//   3. Hom-ref PL correction (zero-PL with sufficient GQ)
//   4. DS (dosage) computation from final PL
//
// Returns: Struct { GT, GQ, DP, PL, DS }
// ============================================================================

#[derive(Debug, PartialEq, Eq, Hash)]
struct VcfGenotypeQcUdf {
    signature: Signature,
}

impl VcfGenotypeQcUdf {
    fn new() -> Self {
        Self {
            signature: Signature::any(5, Volatility::Immutable),
        }
    }
}

impl ScalarUDFImpl for VcfGenotypeQcUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "vcf_genotype_qc"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        if let DataType::Struct(fields) = &arg_types[0] {
            let mut out_fields: Vec<Arc<Field>> = fields.iter().cloned().collect();
            out_fields.push(Arc::new(Field::new(
                "DS",
                DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
                true,
            )));
            Ok(DataType::Struct(out_fields.into()))
        } else {
            Err(DataFusionError::Plan(
                "vcf_genotype_qc: first argument must be a Struct (genotypes)".into(),
            ))
        }
    }

    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        // Extract parameters
        let gq_min = scalar_to_f64(&args.args[1])? as i32;
        let dp_min = scalar_to_f64(&args.args[2])? as i32;
        let dp_max = scalar_to_f64(&args.args[3])? as i32;
        let mask_gt = scalar_to_utf8(&args.args[4])?;

        // Extract struct
        let struct_arr = args.args[0].clone().into_array(1)?;
        let sa = struct_arr
            .as_any()
            .downcast_ref::<StructArray>()
            .ok_or_else(|| DataFusionError::Internal("Expected StructArray".into()))?;

        let gt_col = sa
            .column_by_name("GT")
            .ok_or_else(|| DataFusionError::Internal("Missing GT".into()))?;
        let gq_col = sa
            .column_by_name("GQ")
            .ok_or_else(|| DataFusionError::Internal("Missing GQ".into()))?;
        let dp_col = sa
            .column_by_name("DP")
            .ok_or_else(|| DataFusionError::Internal("Missing DP".into()))?;
        let pl_col = sa
            .column_by_name("PL")
            .ok_or_else(|| DataFusionError::Internal("Missing PL".into()))?;

        let gt_list = gt_col.as_list::<i32>();
        let gq_list = gq_col.as_list::<i32>();
        let dp_list = dp_col.as_list::<i32>();
        let pl_list = pl_col.as_list::<i32>();

        // Flat value arrays
        let gt_flat = gt_list
            .values()
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringArray>()
            .ok_or_else(|| DataFusionError::Internal("GT: expected StringArray".into()))?;
        let gq_flat = gq_list.values().as_primitive::<Int32Type>();
        let dp_flat = dp_list.values().as_primitive::<Int32Type>();

        // PL nested: List<List<Int32>>
        let pl_inner = pl_list
            .values()
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| DataFusionError::Internal("PL: expected List<List<Int32>>".into()))?;
        let pl_inner_off = pl_inner.offsets();
        let pl_flat = pl_inner.values().as_primitive::<Int32Type>();

        let total = gt_flat.len();

        // Output builders
        let mut gt_out = StringBuilder::with_capacity(total, total * 4);
        let mut pl_vals = Int32Builder::with_capacity(total * 3);
        let mut pl_inner_off_out: Vec<i32> = Vec::with_capacity(total + 1);
        pl_inner_off_out.push(0);
        let mut ds_out = Float32Builder::with_capacity(total);

        for i in 0..total {
            // ---- GQ / DP ----
            let gq = if gq_flat.is_null(i) {
                None
            } else {
                Some(gq_flat.value(i))
            };
            let dp = if dp_flat.is_null(i) {
                None
            } else {
                Some(dp_flat.value(i))
            };

            let is_good = gq.is_some_and(|g| g >= gq_min)
                && dp.is_some_and(|d| d >= dp_min && d <= dp_max);

            // ---- GT normalization + masking ----
            let gt_raw = if gt_flat.is_null(i) {
                None
            } else {
                Some(gt_flat.value(i))
            };

            let is_hom_ref = matches!(gt_raw, Some("0/0" | "0|0" | "0"));

            if is_good {
                match gt_raw {
                    Some("0") => gt_out.append_value("0/0"),
                    Some("1") => gt_out.append_value("1/1"),
                    Some("2") => gt_out.append_value("2/2"),
                    Some(v) => gt_out.append_value(v),
                    None => gt_out.append_null(),
                }
            } else {
                gt_out.append_value(&mask_gt);
            }

            // ---- PL extraction ----
            let (pl0, pl1, pl2, pl_valid) = if pl_inner.is_null(i) {
                (0, 0, 0, false)
            } else {
                let s = pl_inner_off[i] as usize;
                let e = pl_inner_off[i + 1] as usize;
                let len = e - s;
                if len >= 3 {
                    (
                        if pl_flat.is_null(s) { 0 } else { pl_flat.value(s) },
                        if pl_flat.is_null(s + 1) { 0 } else { pl_flat.value(s + 1) },
                        if pl_flat.is_null(s + 2) { 0 } else { pl_flat.value(s + 2) },
                        true,
                    )
                } else {
                    (0, 0, 0, false)
                }
            };

            let gq_sufficient = gq.is_some_and(|g| g >= gq_min);

            // ---- PL correction ----
            let needs_correction =
                is_hom_ref && pl_valid && pl0 == 0 && pl1 == 0 && pl2 == 0 && gq_sufficient;

            let (fp0, fp1, fp2) = if needs_correction {
                let gq_f = gq.unwrap_or(0) as f64;
                let p_wrong = 10.0_f64.powf(-gq_f / 10.0);
                let x = (-1.0 + (1.0 + 4.0 * p_wrong).sqrt()) / 2.0;
                let l_ref = 1.0 - x - x * x;
                (
                    pl_phred(l_ref),
                    pl_phred(x),
                    pl_phred(x * x),
                )
            } else {
                (pl0, pl1, pl2)
            };

            pl_vals.append_value(fp0);
            pl_vals.append_value(fp1);
            pl_vals.append_value(fp2);
            pl_inner_off_out.push(*pl_inner_off_out.last().unwrap() + 3);

            // ---- DS ----
            if pl_valid && gq_sufficient {
                let lr = if fp0 < 255 { 10.0_f64.powf(-fp0 as f64 / 10.0) } else { 0.0 };
                let lh = if fp1 < 255 { 10.0_f64.powf(-fp1 as f64 / 10.0) } else { 0.0 };
                let la = if fp2 < 255 { 10.0_f64.powf(-fp2 as f64 / 10.0) } else { 0.0 };
                let sum = lr + lh + la;
                if sum > 0.0 {
                    ds_out.append_value(((lh + 2.0 * la) / sum) as f32);
                } else {
                    ds_out.append_null();
                }
            } else {
                ds_out.append_null();
            }
        }

        // ---- Assemble output ListArrays ----
        let utf8_field = Arc::new(Field::new("item", DataType::Utf8, true));
        let i32_field = Arc::new(Field::new("item", DataType::Int32, true));
        let f32_field = Arc::new(Field::new("item", DataType::Float32, true));

        let gt_out_list = ListArray::new(
            utf8_field,
            gt_list.offsets().clone(),
            Arc::new(gt_out.finish()),
            gt_list.nulls().cloned(),
        );

        let pl_inner_out = ListArray::new(
            i32_field.clone(),
            OffsetBuffer::new(pl_inner_off_out.into()),
            Arc::new(pl_vals.finish()),
            None,
        );
        let pl_out_list = ListArray::new(
            Arc::new(Field::new("item", DataType::List(i32_field), true)),
            pl_list.offsets().clone(),
            Arc::new(pl_inner_out),
            pl_list.nulls().cloned(),
        );

        let ds_out_list = ListArray::new(
            f32_field,
            gt_list.offsets().clone(),
            Arc::new(ds_out.finish()),
            gt_list.nulls().cloned(),
        );

        // ---- Output struct: preserve input field order + metadata, add DS ----
        let input_fields = sa.fields();
        let mut out_fields: Vec<Arc<Field>> = Vec::with_capacity(input_fields.len() + 1);
        let mut out_columns: Vec<ArrayRef> = Vec::with_capacity(input_fields.len() + 1);

        for (fi, field) in input_fields.iter().enumerate() {
            match field.name().as_str() {
                "GT" => {
                    // Preserve metadata from input field
                    out_fields.push(Arc::new(
                        field.as_ref().clone().with_data_type(gt_out_list.data_type().clone()),
                    ));
                    out_columns.push(Arc::new(gt_out_list.clone()));
                }
                "PL" => {
                    out_fields.push(Arc::new(
                        field.as_ref().clone().with_data_type(pl_out_list.data_type().clone()),
                    ));
                    out_columns.push(Arc::new(pl_out_list.clone()));
                }
                _ => {
                    // Passthrough: GQ, DP, or any other field
                    out_fields.push(field.clone());
                    out_columns.push(sa.column(fi).clone());
                }
            }
        }
        // Add DS
        out_fields.push(Arc::new(Field::new(
            "DS",
            ds_out_list.data_type().clone(),
            true,
        )));
        out_columns.push(Arc::new(ds_out_list));

        let out_struct = StructArray::try_new(out_fields.into(), out_columns, sa.nulls().cloned())?;

        Ok(ColumnarValue::Array(Arc::new(out_struct)))
    }
}

/// Convert a likelihood to phred scale, clamped to [0, 255].
#[inline(always)]
fn pl_phred(likelihood: f64) -> i32 {
    (-10.0 * likelihood.max(1e-25).log10())
        .round()
        .max(0.0)
        .min(255.0) as i32
}

pub fn vcf_genotype_qc_udf() -> ScalarUDF {
    ScalarUDF::from(VcfGenotypeQcUdf::new())
}
