//! Vectorized Arrow processing for VCF sample QC.
//!
//! Operates directly on flat value arrays from List columns, avoiding
//! UNNEST/GROUP BY materialization entirely.

use std::sync::OnceLock;

use datafusion::arrow::array::{
    Array, ArrayRef, AsArray, BooleanArray, Float32Array, Float32Builder, Float64Array, Int32Array,
    ListArray, RecordBatch, StringArray, StringBuilder, StructArray,
};
use datafusion::arrow::buffer::{BooleanBuffer, NullBuffer};
use datafusion::arrow::compute;
use datafusion::arrow::datatypes::{DataType, Float64Type, Int32Type};
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

// ============================================================================
// Lookup table for 10^(-PL/10) — avoids 3× powf() per sample (fix #6)
// ============================================================================

/// Pre-computed `10^(-i/10)` for PL values 0..=255.
/// Index 255 maps to 0.0 (sentinel for "no information").
fn pl_likelihood_table() -> &'static [f64; 256] {
    static TABLE: OnceLock<[f64; 256]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = [0.0f64; 256];
        for (i, entry) in t.iter_mut().enumerate().take(255) {
            *entry = 10.0_f64.powf(-(i as f64) / 10.0);
        }
        // t[255] = 0.0 already from initialization
        t
    })
}

// ============================================================================
// Pre-extracted PL data
// ============================================================================

/// Pre-extracted PL data for O(1) access per sample.
struct PlFlatData {
    pl0: Vec<i32>,
    pl1: Vec<i32>,
    pl2: Vec<i32>,
    valid: Vec<bool>,
}

impl PlFlatData {
    /// Pre-extract all PL triples from a filtered `List<List<Int32>>` in one pass.
    fn extract(pl_filtered: &ListArray, total_samples: usize) -> Result<Self> {
        let mut pl0 = Vec::with_capacity(total_samples);
        let mut pl1 = Vec::with_capacity(total_samples);
        let mut pl2 = Vec::with_capacity(total_samples);
        let mut valid = Vec::with_capacity(total_samples);

        let inner_lists = pl_filtered.values();
        let Some(inner_list) = inner_lists.as_any().downcast_ref::<ListArray>() else {
            pl0.resize(total_samples, 0);
            pl1.resize(total_samples, 0);
            pl2.resize(total_samples, 0);
            valid.resize(total_samples, false);
            return Ok(Self {
                pl0,
                pl1,
                pl2,
                valid,
            });
        };

        let flat_values = inner_list.values();
        let inner_offsets = inner_list.offsets();
        let as_i32 = flat_values.as_any().downcast_ref::<Int32Array>();
        let as_f64 = flat_values.as_any().downcast_ref::<Float64Array>();

        for i in 0..total_samples {
            if i >= inner_list.len() || inner_list.is_null(i) {
                pl0.push(0);
                pl1.push(0);
                pl2.push(0);
                valid.push(false);
                continue;
            }

            let start = inner_offsets[i] as usize;
            let end = inner_offsets[i + 1] as usize;
            if end - start < 3 {
                pl0.push(0);
                pl1.push(0);
                pl2.push(0);
                valid.push(false);
                continue;
            }

            if let Some(arr) = as_i32 {
                pl0.push(if arr.is_null(start) {
                    0
                } else {
                    arr.value(start)
                });
                pl1.push(if arr.is_null(start + 1) {
                    0
                } else {
                    arr.value(start + 1)
                });
                pl2.push(if arr.is_null(start + 2) {
                    0
                } else {
                    arr.value(start + 2)
                });
                valid.push(true);
            } else if let Some(arr) = as_f64 {
                pl0.push(if arr.is_null(start) {
                    0
                } else {
                    arr.value(start) as i32
                });
                pl1.push(if arr.is_null(start + 1) {
                    0
                } else {
                    arr.value(start + 1) as i32
                });
                pl2.push(if arr.is_null(start + 2) {
                    0
                } else {
                    arr.value(start + 2) as i32
                });
                valid.push(true);
            } else {
                pl0.push(0);
                pl1.push(0);
                pl2.push(0);
                valid.push(false);
            }
        }

        Ok(Self {
            pl0,
            pl1,
            pl2,
            valid,
        })
    }

    #[inline]
    fn get(&self, idx: usize) -> (i32, i32, i32, bool) {
        (self.pl0[idx], self.pl1[idx], self.pl2[idx], self.valid[idx])
    }
}

// ============================================================================
// GQ/DP flat accessor — avoids f64 cast (fix #8)
// ============================================================================

/// Flat accessor for GQ or DP values, keeping native Int32 when possible.
enum FlatNumeric {
    I32(Int32Array),
    F64(Float64Array),
}

impl FlatNumeric {
    fn from_list(list: &ListArray) -> Result<Self> {
        let values = list.values();
        if let Some(arr) = values.as_any().downcast_ref::<Int32Array>() {
            return Ok(Self::I32(arr.clone()));
        }
        if let Some(arr) = values.as_any().downcast_ref::<Float64Array>() {
            return Ok(Self::F64(arr.clone()));
        }
        // Try casting to Int32 first (cheaper), fall back to Float64
        if let Ok(casted) = compute::cast(values, &DataType::Int32) {
            return Ok(Self::I32(casted.as_primitive::<Int32Type>().clone()));
        }
        let casted = compute::cast(values, &DataType::Float64)?;
        Ok(Self::F64(casted.as_primitive::<Float64Type>().clone()))
    }

    #[inline]
    fn is_null(&self, i: usize) -> bool {
        match self {
            Self::I32(a) => a.is_null(i),
            Self::F64(a) => a.is_null(i),
        }
    }

    /// Get value as i32 (for output). No-op for Int32, truncates Float64.
    #[inline]
    fn value_i32(&self, i: usize) -> i32 {
        match self {
            Self::I32(a) => a.value(i),
            Self::F64(a) => a.value(i) as i32,
        }
    }

    /// Compare value >= threshold (threshold is f64 from config).
    #[inline]
    fn gte(&self, i: usize, threshold: f64) -> bool {
        match self {
            Self::I32(a) => a.value(i) as f64 >= threshold,
            Self::F64(a) => a.value(i) >= threshold,
        }
    }

    /// Compare value <= threshold.
    #[inline]
    fn lte(&self, i: usize, threshold: f64) -> bool {
        match self {
            Self::I32(a) => (a.value(i) as f64) <= threshold,
            Self::F64(a) => a.value(i) <= threshold,
        }
    }

    /// Get value as f64 (for averaging).
    #[inline]
    fn value_f64(&self, i: usize) -> f64 {
        match self {
            Self::I32(a) => a.value(i) as f64,
            Self::F64(a) => a.value(i),
        }
    }
}

// ============================================================================
// Main processing function
// ============================================================================

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

    // ---- Site filtering (fix #3: native Int32 averaging) ----
    let avg_gq = list_avg(&gq_list)?;
    let avg_dp = list_avg(&dp_list)?;
    let site_mask = compute_site_mask(&qual_col, &avg_gq, &avg_dp, config)?;

    let n_passing = site_mask.iter().filter(|v| *v == Some(true)).count();
    if n_passing == 0 {
        return Ok(None);
    }

    // ---- Fix #4: skip filtering when ALL variants pass ----
    let all_pass = n_passing == n_rows;

    let (chrom, start, id, ref_allele, alt, qual, filter_col_arr) = if all_pass {
        (
            get_col(batch, "chrom")?,
            get_col(batch, "start")?,
            get_col(batch, "id")?,
            get_col(batch, "ref")?,
            get_col(batch, "alt")?,
            get_col(batch, "qual")?,
            get_col(batch, "filter")?,
        )
    } else {
        (
            filter_col(batch, "chrom", &site_mask)?,
            filter_col(batch, "start", &site_mask)?,
            filter_col(batch, "id", &site_mask)?,
            filter_col(batch, "ref", &site_mask)?,
            filter_col(batch, "alt", &site_mask)?,
            filter_col(batch, "qual", &site_mask)?,
            filter_col(batch, "filter", &site_mask)?,
        )
    };

    let gt_filtered = if all_pass {
        gt_list
    } else {
        filter_list(&gt_list, &site_mask)?
    };
    let gq_filtered = if all_pass {
        gq_list
    } else {
        filter_list(&gq_list, &site_mask)?
    };
    let dp_filtered = if all_pass {
        dp_list
    } else {
        filter_list(&dp_list, &site_mask)?
    };
    let pl_filtered = if all_pass {
        pl_list
    } else {
        filter_list(&pl_list, &site_mask)?
    };

    // ---- Compute sample offsets ----
    let sample_offsets = compute_offsets(&gt_filtered);
    let total_samples = *sample_offsets.last().unwrap_or(&0);

    // ---- Extract flat value arrays (fix #8: keep native Int32 for GQ/DP) ----
    let flat_gt = get_flat_string_values(&gt_filtered)?;
    let flat_gq = FlatNumeric::from_list(&gq_filtered)?;
    let flat_dp = FlatNumeric::from_list(&dp_filtered)?;

    // ---- Pre-extract PL triples ----
    let pl_data = PlFlatData::extract(&pl_filtered, total_samples)?;

    // ---- Allocate output buffers ----
    let mut gt_builder = StringBuilder::with_capacity(total_samples, total_samples * 4);
    let mut gq_out = Vec::with_capacity(total_samples);
    let mut gq_null = Vec::with_capacity(total_samples);
    let mut dp_out = Vec::with_capacity(total_samples);
    let mut dp_null = Vec::with_capacity(total_samples);
    let mut pl0_out = Vec::with_capacity(total_samples);
    let mut pl1_out = Vec::with_capacity(total_samples);
    let mut pl2_out = Vec::with_capacity(total_samples);
    let mut ds_builder = Float32Builder::with_capacity(total_samples);

    // ---- Process samples: variant-major loop ----
    for var_idx in 0..n_passing {
        let s_start = sample_offsets[var_idx];
        let s_end = sample_offsets[var_idx + 1];

        for i in s_start..s_end {
            // ---- GT normalization ----
            let gt_raw = if flat_gt.is_null(i) {
                None
            } else {
                Some(flat_gt.value(i))
            };
            let gt_norm = normalize_gt(gt_raw);

            // ---- GQ / DP: native Int32 comparison (fix #8) ----
            let gq_present = !flat_gq.is_null(i);
            let dp_present = !flat_dp.is_null(i);

            let is_good = gq_present
                && flat_gq.gte(i, config.sample_gq_min)
                && dp_present
                && flat_dp.gte(i, config.sample_dp_min)
                && flat_dp.lte(i, config.sample_dp_max);

            let gt_final = if is_good { gt_norm } else { "./." };
            gt_builder.append_value(gt_final);

            // ---- GQ / DP output: direct i32, no f64 round-trip ----
            if gq_present {
                gq_out.push(flat_gq.value_i32(i));
                gq_null.push(true);
            } else {
                gq_out.push(0);
                gq_null.push(false);
            }
            if dp_present {
                dp_out.push(flat_dp.value_i32(i));
                dp_null.push(true);
            } else {
                dp_out.push(0);
                dp_null.push(false);
            }

            // ---- PL extraction: O(1) ----
            let (raw_pl0, raw_pl1, raw_pl2, pl_valid) = pl_data.get(i);

            let is_hom_ref = matches!(gt_raw, Some("0/0" | "0|0" | "0"));
            let gq_sufficient = gq_present && flat_gq.gte(i, config.sample_gq_min);

            let needs_correction = is_hom_ref
                && pl_valid
                && raw_pl0 == 0
                && raw_pl1 == 0
                && raw_pl2 == 0
                && gq_sufficient;

            let (final_pl0, final_pl1, final_pl2) = if needs_correction {
                correct_pl(flat_gq.value_f64(i))
            } else {
                (raw_pl0, raw_pl1, raw_pl2)
            };

            pl0_out.push(final_pl0);
            pl1_out.push(final_pl1);
            pl2_out.push(final_pl2);

            // ---- DS: lookup table instead of powf (fix #6) ----
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

fn get_genotypes_struct(batch: &RecordBatch) -> Result<&StructArray> {
    let col = batch
        .column_by_name("genotypes")
        .ok_or_else(|| DataFusionError::Execution("Missing 'genotypes' column".into()))?;
    col.as_any()
        .downcast_ref::<StructArray>()
        .ok_or_else(|| DataFusionError::Execution("'genotypes' is not a Struct".into()))
}

fn get_struct_list_field(s: &StructArray, name: &str) -> Result<ListArray> {
    let col = s
        .column_by_name(name)
        .ok_or_else(|| DataFusionError::Execution(format!("Missing genotypes.'{name}'")))?;
    col.as_any()
        .downcast_ref::<ListArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution(format!("genotypes.'{name}' is not a List")))
}

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

/// Get a column by name (no filtering).
fn get_col(batch: &RecordBatch, name: &str) -> Result<ArrayRef> {
    batch
        .column_by_name(name)
        .cloned()
        .ok_or_else(|| DataFusionError::Execution(format!("Missing '{name}' column")))
}

/// Compute per-variant average of a List<numeric> column.
///
/// Handles Int32 natively (sum as i64) to avoid f64 cast overhead (fix #3).
fn list_avg(list: &ListArray) -> Result<Float64Array> {
    let offsets = list.offsets();
    let values = list.values();
    let n = list.len();
    let mut avgs = Vec::with_capacity(n);
    let mut nulls = Vec::with_capacity(n);

    // Int32 fast path: sum as i64, divide once
    if let Some(i32_values) = values.as_any().downcast_ref::<Int32Array>() {
        for i in 0..n {
            if list.is_null(i) {
                avgs.push(0.0);
                nulls.push(false);
                continue;
            }
            let start = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            if end == start {
                avgs.push(0.0);
                nulls.push(false);
                continue;
            }
            let mut sum = 0i64;
            let mut count = 0usize;
            for j in start..end {
                if !i32_values.is_null(j) {
                    sum += i32_values.value(j) as i64;
                    count += 1;
                }
            }
            if count > 0 {
                avgs.push(sum as f64 / count as f64);
                nulls.push(true);
            } else {
                avgs.push(0.0);
                nulls.push(false);
            }
        }
    } else {
        // Float64 path (fallback)
        let f64_values = if values.data_type() == &DataType::Float64 {
            values.as_primitive::<Float64Type>().clone()
        } else {
            let casted = compute::cast(values, &DataType::Float64)?;
            casted.as_primitive::<Float64Type>().clone()
        };
        for i in 0..n {
            if list.is_null(i) {
                avgs.push(0.0);
                nulls.push(false);
                continue;
            }
            let start = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            if end == start {
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
    }

    let values_buf = Float64Array::from(avgs);
    let null_buf = NullBuffer::from(BooleanBuffer::from(nulls));
    Ok(Float64Array::new(values_buf.into_parts().1, Some(null_buf)))
}

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

fn filter_col(batch: &RecordBatch, name: &str, mask: &BooleanArray) -> Result<ArrayRef> {
    let col = batch
        .column_by_name(name)
        .ok_or_else(|| DataFusionError::Execution(format!("Missing '{name}' column")))?;
    compute::filter(col, mask)
        .map_err(|e| DataFusionError::Execution(format!("filter {name}: {e}")))
}

fn filter_list(list: &ListArray, mask: &BooleanArray) -> Result<ListArray> {
    let filtered = compute::filter(list, mask)?;
    filtered
        .as_any()
        .downcast_ref::<ListArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution("filter_list: downcast failed".into()))
}

fn compute_offsets(list: &ListArray) -> Vec<usize> {
    let offsets = list.offsets();
    let base = offsets[0] as usize;
    offsets.iter().map(|o| *o as usize - base).collect()
}

fn get_flat_string_values(list: &ListArray) -> Result<StringArray> {
    let values = list.values();
    if let Some(s) = values.as_any().downcast_ref::<StringArray>() {
        return Ok(s.clone());
    }
    let casted = compute::cast(values, &DataType::Utf8)?;
    casted
        .as_any()
        .downcast_ref::<StringArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution("Cannot get string values from GT list".into()))
}

#[inline]
fn normalize_gt(gt: Option<&str>) -> &'static str {
    match gt {
        Some("0") => "0/0",
        Some("1") => "1/1",
        Some("2") => "2/2",
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
        Some(_) => "OTHER",
    }
}

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

/// Compute dosage from final PL values using a lookup table (fix #6).
///
/// Avoids 3× `powf()` per sample by indexing into a pre-computed table.
#[inline]
fn compute_dosage(pl0: i32, pl1: i32, pl2: i32) -> Option<f32> {
    let table = pl_likelihood_table();
    let l_ref = table[pl0.clamp(0, 255) as usize];
    let l_het = table[pl1.clamp(0, 255) as usize];
    let l_alt = table[pl2.clamp(0, 255) as usize];
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
        assert!(pl0 >= 0 && pl0 <= 1);
        assert!(pl1 > 20 && pl1 < 40);
        assert!(pl2 > 50 && pl2 < 70);
    }

    #[test]
    fn test_compute_dosage() {
        assert!((compute_dosage(0, 0, 0).unwrap() - 1.0).abs() < 0.01);

        let ds = compute_dosage(0, 30, 60).unwrap();
        assert!(ds < 0.01);

        let ds = compute_dosage(60, 30, 0).unwrap();
        assert!(ds > 1.99);

        assert!(compute_dosage(255, 255, 255).is_none());
    }

    #[test]
    fn test_pl_likelihood_table() {
        let table = pl_likelihood_table();
        assert!((table[0] - 1.0).abs() < 1e-10);
        assert!((table[10] - 0.1).abs() < 1e-10);
        assert!((table[20] - 0.01).abs() < 1e-10);
        assert_eq!(table[255], 0.0);
    }
}
