//! Vectorized Arrow processing for VCF sample QC.
//!
//! Operates directly on flat value arrays from List columns, avoiding
//! UNNEST/GROUP BY materialization entirely.

use datafusion::arrow::array::{
    Array, ArrayRef, AsArray, BooleanArray, Float32Array, Float32Builder, Float64Array, Int32Array,
    ListArray, RecordBatch, StringArray, StringBuilder, StructArray,
};
use datafusion::arrow::buffer::{BooleanBuffer, NullBuffer};
use datafusion::arrow::compute;
use datafusion::arrow::datatypes::{DataType, Float64Type};
use datafusion::common::{DataFusionError, Result};

use super::QcConfig;

/// Result of processing a single batch through site filtering + sample QC.
pub struct ProcessedBatch {
    /// Number of variants that passed site filtering.
    pub n_variants: usize,

    // ---- Site-level columns (one value per variant) ----
    pub chrom: ArrayRef,
    pub start: ArrayRef,
    pub id: ArrayRef,
    pub ref_allele: ArrayRef,
    pub alt: ArrayRef,
    pub qual: ArrayRef,
    pub filter: ArrayRef,

    // ---- Sample-level flat arrays ----
    /// Offsets into the flat arrays (length = n_variants + 1).
    /// Samples for variant i are at indices offsets[i]..offsets[i+1].
    pub sample_offsets: Vec<usize>,

    pub gt_final: StringArray,
    pub gq_values: Int32Array,
    pub dp_values: Int32Array,
    pub pl0: Int32Array,
    pub pl1: Int32Array,
    pub pl2: Int32Array,
    pub ds_values: Float32Array,
}

/// Process a RecordBatch through site filtering and sample QC.
///
/// Returns `None` if no variants pass site filtering.
pub fn process_batch(batch: &RecordBatch, config: &QcConfig) -> Result<Option<ProcessedBatch>> {
    let n_rows = batch.num_rows();
    if n_rows == 0 {
        return Ok(None);
    }

    // ---- Extract site-level columns ----
    let qual_col = get_column_f64(batch, "qual")?;

    // ---- Extract genotypes struct ----
    let genotypes = get_genotypes_struct(batch)?;
    let gt_list = get_struct_list_field(genotypes, "GT")?;
    let gq_list = get_struct_list_field(genotypes, "GQ")?;
    let dp_list = get_struct_list_field(genotypes, "DP")?;
    let pl_list = get_struct_list_field(genotypes, "PL")?;

    // ---- Site filtering ----
    let avg_gq = list_avg_f64(&gq_list)?;
    let avg_dp = list_avg_f64(&dp_list)?;
    let site_mask = compute_site_mask(&qual_col, &avg_gq, &avg_dp, config)?;

    let n_passing = site_mask.iter().filter(|v| *v == Some(true)).count();
    if n_passing == 0 {
        return Ok(None);
    }

    // ---- Filter site-level columns ----
    let chrom = filter_col(batch, "chrom", &site_mask)?;
    let start = filter_col(batch, "start", &site_mask)?;
    let id = filter_col(batch, "id", &site_mask)?;
    let ref_allele = filter_col(batch, "ref", &site_mask)?;
    let alt = filter_col(batch, "alt", &site_mask)?;
    let qual = compute::filter(
        batch
            .column_by_name("qual")
            .ok_or_else(|| DataFusionError::Execution("Missing 'qual' column".into()))?,
        &site_mask,
    )?;
    let filter_col_arr = filter_col(batch, "filter", &site_mask)?;

    // ---- Filter genotype list arrays ----
    let gt_filtered = filter_list(&gt_list, &site_mask)?;
    let gq_filtered = filter_list(&gq_list, &site_mask)?;
    let dp_filtered = filter_list(&dp_list, &site_mask)?;
    let pl_filtered = filter_list(&pl_list, &site_mask)?;

    // ---- Compute sample offsets from GT list (all lists share same offsets) ----
    let sample_offsets = compute_offsets(&gt_filtered);
    let total_samples = *sample_offsets.last().unwrap_or(&0);

    // ---- Extract flat value arrays ----
    let flat_gt = get_flat_string_values(&gt_filtered)?;
    let flat_gq = get_flat_f64_values(&gq_filtered)?;
    let flat_dp = get_flat_f64_values(&dp_filtered)?;

    // ---- Process all samples in one vectorized pass ----
    let mut gt_builder = StringBuilder::with_capacity(total_samples, total_samples * 4);
    let mut gq_out = Vec::with_capacity(total_samples);
    let mut gq_null = Vec::with_capacity(total_samples);
    let mut dp_out = Vec::with_capacity(total_samples);
    let mut dp_null = Vec::with_capacity(total_samples);
    let mut pl0_out = Vec::with_capacity(total_samples);
    let mut pl1_out = Vec::with_capacity(total_samples);
    let mut pl2_out = Vec::with_capacity(total_samples);
    let mut ds_builder = Float32Builder::with_capacity(total_samples);

    for i in 0..total_samples {
        // ---- GT normalization ----
        let gt_raw = if flat_gt.is_null(i) {
            None
        } else {
            Some(flat_gt.value(i))
        };
        let gt_norm = normalize_gt(gt_raw);

        // ---- GQ / DP values ----
        let gq_val = if flat_gq.is_null(i) {
            None
        } else {
            Some(flat_gq.value(i))
        };
        let dp_val = if flat_dp.is_null(i) {
            None
        } else {
            Some(flat_dp.value(i))
        };

        // ---- Sample quality ----
        let is_good = gq_val.is_some_and(|g| g >= config.sample_gq_min)
            && dp_val.is_some_and(|d| d >= config.sample_dp_min && d <= config.sample_dp_max);

        let gt_final = if is_good { gt_norm } else { "./." };
        gt_builder.append_value(gt_final);

        // ---- GQ / DP output (preserve originals) ----
        if let Some(g) = gq_val {
            gq_out.push(g as i32);
            gq_null.push(true);
        } else {
            gq_out.push(0);
            gq_null.push(false);
        }
        if let Some(d) = dp_val {
            dp_out.push(d as i32);
            dp_null.push(true);
        } else {
            dp_out.push(0);
            dp_null.push(false);
        }

        // ---- PL extraction and correction ----
        let (raw_pl0, raw_pl1, raw_pl2, pl_valid) =
            extract_pl_triple(&pl_filtered, &sample_offsets, i);

        let is_hom_ref = matches!(gt_raw, Some("0/0" | "0|0" | "0"));
        let gq_sufficient = gq_val.is_some_and(|g| g >= config.sample_gq_min);

        let needs_correction =
            is_hom_ref && pl_valid && raw_pl0 == 0 && raw_pl1 == 0 && raw_pl2 == 0 && gq_sufficient;

        let (final_pl0, final_pl1, final_pl2) = if needs_correction {
            correct_pl(gq_val.unwrap_or(0.0))
        } else {
            (raw_pl0, raw_pl1, raw_pl2)
        };

        pl0_out.push(final_pl0);
        pl1_out.push(final_pl1);
        pl2_out.push(final_pl2);

        // ---- DS (dosage) ----
        if pl_valid && gq_sufficient {
            let ds = compute_dosage(final_pl0, final_pl1, final_pl2);
            match ds {
                Some(d) => ds_builder.append_value(d),
                None => ds_builder.append_null(),
            }
        } else {
            ds_builder.append_null();
        }
    }

    // ---- Build output arrays ----
    let gq_array = Int32Array::from(gq_out);
    let gq_nulls = NullBuffer::from(BooleanBuffer::from(gq_null));
    let gq_array = Int32Array::new(gq_array.into_parts().1, Some(gq_nulls));

    let dp_array = Int32Array::from(dp_out);
    let dp_nulls = NullBuffer::from(BooleanBuffer::from(dp_null));
    let dp_array = Int32Array::new(dp_array.into_parts().1, Some(dp_nulls));

    Ok(Some(ProcessedBatch {
        n_variants: n_passing,
        chrom,
        start,
        id,
        ref_allele,
        alt,
        qual,
        filter: filter_col_arr,
        sample_offsets,
        gt_final: gt_builder.finish(),
        gq_values: gq_array,
        dp_values: dp_array,
        pl0: Int32Array::from(pl0_out),
        pl1: Int32Array::from(pl1_out),
        pl2: Int32Array::from(pl2_out),
        ds_values: ds_builder.finish(),
    }))
}

// ============================================================================
// Helper functions
// ============================================================================

/// Get the `genotypes` struct column from a batch.
fn get_genotypes_struct(batch: &RecordBatch) -> Result<&StructArray> {
    let col = batch
        .column_by_name("genotypes")
        .ok_or_else(|| DataFusionError::Execution("Missing 'genotypes' column".into()))?;
    col.as_any()
        .downcast_ref::<StructArray>()
        .ok_or_else(|| DataFusionError::Execution("'genotypes' is not a Struct".into()))
}

/// Get a List field from a StructArray.
fn get_struct_list_field(s: &StructArray, name: &str) -> Result<ListArray> {
    let col = s
        .column_by_name(name)
        .ok_or_else(|| DataFusionError::Execution(format!("Missing genotypes.'{name}'")))?;
    col.as_any()
        .downcast_ref::<ListArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution(format!("genotypes.'{name}' is not a List")))
}

/// Get a column as Float64Array (with casting if needed).
fn get_column_f64(batch: &RecordBatch, name: &str) -> Result<Float64Array> {
    let col = batch
        .column_by_name(name)
        .ok_or_else(|| DataFusionError::Execution(format!("Missing '{name}' column")))?;
    if let Some(f) = col.as_any().downcast_ref::<Float64Array>() {
        return Ok(f.clone());
    }
    let casted = compute::cast(col, &DataType::Float64)?;
    Ok(casted.as_primitive::<Float64Type>().clone())
}

/// Compute per-variant average of a List<numeric> column, returning Float64Array.
fn list_avg_f64(list: &ListArray) -> Result<Float64Array> {
    let offsets = list.offsets();
    let values = list.values();
    let f64_values = if values.data_type() == &DataType::Float64 {
        values.as_primitive::<Float64Type>().clone()
    } else {
        let casted = compute::cast(values, &DataType::Float64)?;
        casted.as_primitive::<Float64Type>().clone()
    };

    let n = list.len();
    let mut avgs = Vec::with_capacity(n);
    let mut nulls = Vec::with_capacity(n);

    for i in 0..n {
        if list.is_null(i) {
            avgs.push(0.0);
            nulls.push(false);
            continue;
        }
        let start = offsets[i] as usize;
        let end = offsets[i + 1] as usize;
        let len = end - start;
        if len == 0 {
            avgs.push(0.0);
            nulls.push(false);
            continue;
        }
        let mut sum = 0.0;
        let mut count = 0usize;
        for j in start..end {
            if !f64_values.is_null(j) {
                sum += f64_values.value(j);
                count += 1;
            }
        }
        if count > 0 {
            avgs.push(sum / count as f64);
            nulls.push(true);
        } else {
            avgs.push(0.0);
            nulls.push(false);
        }
    }

    let values_buf = Float64Array::from(avgs);
    let null_buf = NullBuffer::from(BooleanBuffer::from(nulls));
    Ok(Float64Array::new(values_buf.into_parts().1, Some(null_buf)))
}

/// Compute the site-level filter mask.
fn compute_site_mask(
    qual: &Float64Array,
    avg_gq: &Float64Array,
    avg_dp: &Float64Array,
    config: &QcConfig,
) -> Result<BooleanArray> {
    let n = qual.len();
    let mut mask = Vec::with_capacity(n);
    for i in 0..n {
        let passes = !qual.is_null(i)
            && qual.value(i) >= config.qual_min
            && !avg_gq.is_null(i)
            && avg_gq.value(i) >= config.site_gq_avg_min
            && !avg_dp.is_null(i)
            && avg_dp.value(i) >= config.site_dp_avg_min
            && avg_dp.value(i) <= config.site_dp_avg_max;
        mask.push(passes);
    }
    Ok(BooleanArray::from(mask))
}

/// Filter a column from a batch by a boolean mask.
fn filter_col(batch: &RecordBatch, name: &str, mask: &BooleanArray) -> Result<ArrayRef> {
    let col = batch
        .column_by_name(name)
        .ok_or_else(|| DataFusionError::Execution(format!("Missing '{name}' column")))?;
    compute::filter(col, mask)
        .map_err(|e| DataFusionError::Execution(format!("filter {name}: {e}")))
}

/// Filter a ListArray by a boolean mask.
fn filter_list(list: &ListArray, mask: &BooleanArray) -> Result<ListArray> {
    let filtered = compute::filter(list, mask)?;
    filtered
        .as_any()
        .downcast_ref::<ListArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution("filter_list: downcast failed".into()))
}

/// Compute sample offsets from a filtered ListArray.
fn compute_offsets(list: &ListArray) -> Vec<usize> {
    let offsets = list.offsets();
    let base = offsets[0] as usize;
    offsets.iter().map(|o| *o as usize - base).collect()
}

/// Get flat string values from a List<Utf8> or List<LargeUtf8>.
fn get_flat_string_values(list: &ListArray) -> Result<StringArray> {
    let values = list.values();
    if let Some(s) = values.as_any().downcast_ref::<StringArray>() {
        return Ok(s.clone());
    }
    // Try casting from other string types
    let casted = compute::cast(values, &DataType::Utf8)?;
    casted
        .as_any()
        .downcast_ref::<StringArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution("Cannot get string values from GT list".into()))
}

/// Get flat f64 values from a List<numeric>.
fn get_flat_f64_values(list: &ListArray) -> Result<Float64Array> {
    let values = list.values();
    if let Some(f) = values.as_any().downcast_ref::<Float64Array>() {
        return Ok(f.clone());
    }
    let casted = compute::cast(values, &DataType::Float64)?;
    Ok(casted.as_primitive::<Float64Type>().clone())
}

/// Normalize a genotype string: haploid shorthand to diploid.
#[inline]
fn normalize_gt(gt: Option<&str>) -> &'static str {
    match gt {
        Some("0") => "0/0",
        Some("1") => "1/1",
        Some("2") => "2/2",
        // Return a static string for common genotypes; for uncommon ones we need allocation.
        // Since we write to VCF immediately, we use static strings for the common cases.
        Some("0/0") => "0/0",
        Some("0|0") => "0|0",
        Some("0/1") => "0/1",
        Some("0|1") => "0|1",
        Some("1/0") => "1/0",
        Some("1|0") => "1|0",
        Some("1/1") => "1/1",
        Some("1|1") => "1|1",
        Some("./." | ".|.") => "./.",
        Some(".") => "./.",
        None => "./.",
        // For rare genotypes, return as-is (caller handles non-static case)
        Some(_) => "OTHER",
    }
}

/// Extract PL triple (pl0, pl1, pl2) for a sample at flat index `sample_idx`.
///
/// The `pl_filtered` is the filtered PL ListArray (List<List<Int32>>).
/// We need to index into the inner list for this sample.
fn extract_pl_triple(
    pl_filtered: &ListArray,
    sample_offsets: &[usize],
    sample_idx: usize,
) -> (i32, i32, i32, bool) {
    // Find which variant this sample belongs to
    let variant_idx = sample_offsets
        .windows(2)
        .position(|w| sample_idx >= w[0] && sample_idx < w[1])
        .unwrap_or(0);

    // Get the PL list for this variant
    let pl_offsets = pl_filtered.offsets();
    let base = pl_offsets[0] as usize;
    let variant_pl_start = pl_offsets[variant_idx] as usize - base;
    let local_sample = sample_idx - sample_offsets[variant_idx];
    let inner_idx = variant_pl_start + local_sample;

    // The values of pl_filtered is itself a ListArray (the per-sample PL arrays)
    let inner_lists = pl_filtered.values();
    if let Some(inner_list) = inner_lists.as_any().downcast_ref::<ListArray>() {
        if inner_idx >= inner_list.len() || inner_list.is_null(inner_idx) {
            return (0, 0, 0, false);
        }
        let pl_values = inner_list.value(inner_idx);
        let pl_len = pl_values.len();
        if pl_len < 3 {
            return (0, 0, 0, false);
        }
        if let Some(int_arr) = pl_values.as_any().downcast_ref::<Int32Array>() {
            let pl0 = if int_arr.is_null(0) {
                0
            } else {
                int_arr.value(0)
            };
            let pl1 = if int_arr.is_null(1) {
                0
            } else {
                int_arr.value(1)
            };
            let pl2 = if int_arr.is_null(2) {
                0
            } else {
                int_arr.value(2)
            };
            return (pl0, pl1, pl2, true);
        }
        // Try f64 and cast
        if let Some(f_arr) = pl_values.as_any().downcast_ref::<Float64Array>() {
            let pl0 = if f_arr.is_null(0) {
                0
            } else {
                f_arr.value(0) as i32
            };
            let pl1 = if f_arr.is_null(1) {
                0
            } else {
                f_arr.value(1) as i32
            };
            let pl2 = if f_arr.is_null(2) {
                0
            } else {
                f_arr.value(2) as i32
            };
            return (pl0, pl1, pl2, true);
        }
    }
    (0, 0, 0, false)
}

/// Correct PL values for homozygous-reference sites with all-zero PLs.
///
/// Implements the quadratic correction from GQ.
#[inline]
fn correct_pl(gq: f64) -> (i32, i32, i32) {
    let p_wrong = 10.0_f64.powf(-gq / 10.0);
    let x = (-1.0 + (1.0 + 4.0 * p_wrong).sqrt()) / 2.0;

    let val0 = (1.0 - x - x * x).max(1e-25);
    let val1 = x.max(1e-25);
    let val2 = (x * x).max(1e-25);

    let pl0 = (-10.0 * val0.log10()).round().clamp(0.0, 255.0) as i32;
    let pl1 = (-10.0 * val1.log10()).round().clamp(0.0, 255.0) as i32;
    let pl2 = (-10.0 * val2.log10()).round().clamp(0.0, 255.0) as i32;

    (pl0, pl1, pl2)
}

/// Compute dosage from final PL values.
///
/// Returns `None` if the total likelihood is zero.
#[inline]
fn compute_dosage(pl0: i32, pl1: i32, pl2: i32) -> Option<f32> {
    let l_ref = if pl0 < 255 {
        10.0_f64.powf(-pl0 as f64 / 10.0)
    } else {
        0.0
    };
    let l_het = if pl1 < 255 {
        10.0_f64.powf(-pl1 as f64 / 10.0)
    } else {
        0.0
    };
    let l_alt = if pl2 < 255 {
        10.0_f64.powf(-pl2 as f64 / 10.0)
    } else {
        0.0
    };
    let total = l_ref + l_het + l_alt;
    if total > 0.0 {
        Some(((l_het + 2.0 * l_alt) / total) as f32)
    } else {
        None
    }
}

/// Public test helpers for integration tests.
pub mod processing_test_helpers {
    use super::*;

    pub fn correct_pl_public(gq: f64) -> (i32, i32, i32) {
        correct_pl(gq)
    }

    pub fn compute_dosage_public(pl0: i32, pl1: i32, pl2: i32) -> Option<f32> {
        compute_dosage(pl0, pl1, pl2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_gt() {
        assert_eq!(normalize_gt(Some("0")), "0/0");
        assert_eq!(normalize_gt(Some("1")), "1/1");
        assert_eq!(normalize_gt(Some("2")), "2/2");
        assert_eq!(normalize_gt(Some("0/0")), "0/0");
        assert_eq!(normalize_gt(Some("0/1")), "0/1");
        assert_eq!(normalize_gt(Some("./.")), "./.");
        assert_eq!(normalize_gt(None), "./.");
    }

    #[test]
    fn test_correct_pl() {
        let (pl0, pl1, pl2) = correct_pl(30.0);
        // With GQ=30: p_wrong=0.001, x≈0.000999
        // PL0 should be ~0 (high confidence ref)
        // PL1 should be ~30 (het penalty from GQ)
        // PL2 should be ~60 (alt penalty = 2*GQ approx)
        assert!(pl0 >= 0 && pl0 <= 1);
        assert!(pl1 > 20 && pl1 < 40);
        assert!(pl2 > 50 && pl2 < 70);
    }

    #[test]
    fn test_compute_dosage() {
        // All-zero PLs: equal likelihoods → DS = 1.0
        assert!((compute_dosage(0, 0, 0).unwrap() - 1.0).abs() < 0.01);

        // Hom-ref (PL0=0, PL1=30, PL2=60): DS ≈ 0
        let ds = compute_dosage(0, 30, 60).unwrap();
        assert!(ds < 0.01);

        // Hom-alt (PL0=60, PL1=30, PL2=0): DS ≈ 2
        let ds = compute_dosage(60, 30, 0).unwrap();
        assert!(ds > 1.99);

        // All 255: total = 0 → None
        assert!(compute_dosage(255, 255, 255).is_none());
    }
}
