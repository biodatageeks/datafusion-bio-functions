//! Variation-cache column encoders for the Parquet backend.
//!
//! Rust equivalent of the validated chr1 encoder (measured smallest + fastest +
//! lossless). Three transforms + one formatter:
//! - [`presence_boolean`] — binary flags (`Int8` presence) → non-nullable `Boolean`;
//! - [`encode_af_2array`] — `allele:freq|...` AF string → struct-of-arrays
//!   (`List<Utf8>` alleles + `List<Float32>` freqs, positional per population);
//! - [`dedup_variation_name`] — null `variation_name` where it equals `dbsnp_ids`
//!   (reconstructed on read via `coalesce`);
//! - [`format_g4`] — C `printf("%.4g", …)` equivalent; AF is stored to 4
//!   significant figures and the reader reproduces the exact CSQ text with this.

use datafusion::arrow::array::{
    Array, BooleanArray, Float32Builder, Int8Array, ListArray, ListBuilder, StringArray,
    StringBuilder,
};
use datafusion::common::{DataFusionError, Result};

/// Format `f` like C `printf("%.4g", f)`: 4 significant figures, trailing zeros
/// (and a trailing `.`) stripped, scientific notation when the decimal exponent
/// is `< -4` or `>= 4` (exponent as `e-05` / `e+10`, signed, min 2 digits).
/// `0.0` → `"0"`.
pub fn format_g4(f: f32) -> String {
    format_g(f, 4)
}

fn format_g(f: f32, p: usize) -> String {
    if f == 0.0 {
        return "0".to_string();
    }
    // Scientific with `p-1` mantissa decimals rounds to `p` significant figures
    // and exposes the decimal exponent of the leading digit.
    let sci = format!("{:.*e}", p - 1, f); // e.g. 0.08806 -> "8.806e-2"
    let (mant, exp) = sci.split_once('e').expect("scientific format has 'e'");
    let e: i32 = exp.parse().expect("exponent parses");

    if e < -4 || e >= p as i32 {
        // Scientific: strip mantissa trailing zeros, format exponent.
        format!("{}e{}", strip_trailing_zeros(mant), format_exp(e))
    } else {
        // Fixed with `p-1-e` decimals (formatting the original value rounds
        // consistently with the exponent picked above).
        let decimals = (p as i32 - 1 - e).max(0) as usize;
        strip_trailing_zeros(&format!("{:.*}", decimals, f))
    }
}

/// Strip trailing zeros and a dangling decimal point from a decimal string.
fn strip_trailing_zeros(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

/// Format a decimal exponent as `e`-suffix style: sign + at least 2 digits.
fn format_exp(e: i32) -> String {
    let sign = if e < 0 { '-' } else { '+' };
    format!("{}{:02}", sign, e.abs())
}

/// Convert a binary variation flag column (`Int8`, values `null` or `1`) to a
/// non-nullable presence `Boolean` (non-null → `true`, null → `false`).
///
/// Fails if the column is not `Int8` or any non-null value is not exactly `1`
/// (the data is verified to contain only `null`/`1` across chr1–22; the guard
/// prevents a future source silently corrupting the presence mapping).
pub fn presence_boolean(col: &dyn Array, name: &str) -> Result<BooleanArray> {
    let arr = col.as_any().downcast_ref::<Int8Array>().ok_or_else(|| {
        DataFusionError::Execution(format!(
            "variation flag {name} must be Int8, got {:?}",
            col.data_type()
        ))
    })?;
    let mut values = Vec::with_capacity(arr.len());
    for i in 0..arr.len() {
        if arr.is_null(i) {
            values.push(false);
        } else {
            let v = arr.value(i);
            if v != 1 {
                return Err(DataFusionError::Execution(format!(
                    "variation flag {name} has non-1 value {v}: presence invariant broken"
                )));
            }
            values.push(true);
        }
    }
    Ok(BooleanArray::from(values))
}

/// The struct-of-arrays AF encoding for one source column.
pub struct AfArrays {
    /// `List<Utf8>` — the alt allele(s) present in this row (usually one).
    pub alleles: ListArray,
    /// `List<Float32>` — allele-major, `n_pops` positional slots per allele
    /// (`null` where a population is missing or lacks that allele).
    pub freqs: ListArray,
}

/// Encode an AF source string column (`allele:freq|...`, pipe-separated per
/// population, comma-separated per allele, empty entry = missing population)
/// into the [`AfArrays`] struct-of-arrays.
///
/// Losslessness is asserted per freq token: `format_g4(parsed) == token`
/// (chr1 is 100% clean; the `af_overflow_raw` fallback is a later task).
pub fn encode_af_2array(col: &StringArray, n_pops: usize) -> Result<AfArrays> {
    let mut alleles = ListBuilder::new(StringBuilder::new());
    let mut freqs = ListBuilder::new(Float32Builder::new());

    for i in 0..col.len() {
        if col.is_null(i) || col.value(i).is_empty() {
            alleles.append(false);
            freqs.append(false);
            continue;
        }
        let s = col.value(i);

        // Per-population allele→freq lists (`None` = missing population).
        let mut pops: Vec<Option<Vec<(String, f32)>>> = Vec::new();
        for entry in s.split('|') {
            if entry.is_empty() {
                pops.push(None);
                continue;
            }
            let mut pairs = Vec::new();
            for pair in entry.split(',') {
                let (allele, freq) = pair.rsplit_once(':').ok_or_else(|| {
                    DataFusionError::Execution(format!("malformed AF pair '{pair}'"))
                })?;
                let value: f32 = freq.parse().map_err(|_| {
                    DataFusionError::Execution(format!("unparseable AF frequency '{freq}'"))
                })?;
                if format_g4(value) != freq {
                    return Err(DataFusionError::Execution(format!(
                        "AF frequency '{freq}' is not %.4g-round-trippable (got '{}')",
                        format_g4(value)
                    )));
                }
                pairs.push((allele.to_string(), value));
            }
            pops.push(Some(pairs));
        }

        // Ordered allele set (first-seen order across populations).
        let mut order: Vec<String> = Vec::new();
        for pop in pops.iter().flatten() {
            for (allele, _) in pop {
                if !order.iter().any(|a| a == allele) {
                    order.push(allele.clone());
                }
            }
        }

        for allele in &order {
            alleles.values().append_value(allele);
        }
        alleles.append(true);

        // Allele-major: for each allele, `n_pops` positional freq slots.
        for allele in &order {
            for pop_idx in 0..n_pops {
                match pops.get(pop_idx) {
                    Some(Some(pairs)) => match pairs.iter().find(|(a, _)| a == allele) {
                        Some((_, f)) => freqs.values().append_value(*f),
                        None => freqs.values().append_null(),
                    },
                    _ => freqs.values().append_null(),
                }
            }
        }
        freqs.append(true);
    }

    Ok(AfArrays {
        alleles: alleles.finish(),
        freqs: freqs.finish(),
    })
}

/// Null out `variation_name` entries that exactly equal `dbsnp_ids` (they are
/// losslessly reconstructed on read via `coalesce(variation_name, dbsnp_ids)`).
pub fn dedup_variation_name(vn: &StringArray, dbsnp: &StringArray) -> StringArray {
    let out: Vec<Option<&str>> = (0..vn.len())
        .map(|i| {
            if vn.is_null(i) {
                None
            } else {
                let v = vn.value(i);
                if !dbsnp.is_null(i) && dbsnp.value(i) == v {
                    None
                } else {
                    Some(v)
                }
            }
        })
        .collect();
    StringArray::from(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_g4_matches_measured_tokens() {
        let cases: &[(f32, &str)] = &[
            (0.08806, "0.08806"),
            (0.4253, "0.4253"),
            (0.0367, "0.0367"),
            (0.9969, "0.9969"),
            (0.006912, "0.006912"),
            (0.0, "0"),
            (2.682e-05, "2.682e-05"),
            (9.911e-05, "9.911e-05"),
            (4.248e-06, "4.248e-06"),
            (5.896e-05, "5.896e-05"),
            (0.5, "0.5"),
            (1.0, "1"),
        ];
        for (f, expected) in cases {
            assert_eq!(&format_g4(*f), expected, "format_g4({f})");
        }
    }

    #[test]
    fn presence_boolean_maps_presence_and_guards_invariant() {
        let arr = Int8Array::from(vec![None, Some(1), Some(1), None]);
        let b = presence_boolean(&arr, "failed").unwrap();
        assert_eq!(b.len(), 4);
        assert_eq!(b.null_count(), 0); // non-nullable
        assert!(!b.value(0) && b.value(1) && b.value(2) && !b.value(3));

        assert!(presence_boolean(&Int8Array::from(vec![Some(2)]), "somatic").is_err());
        assert!(presence_boolean(&Int8Array::from(vec![Some(0)]), "somatic").is_err());
        // wrong type
        assert!(presence_boolean(&StringArray::from(vec!["x"]), "failed").is_err());
    }

    fn list_row_f32(arr: &ListArray, row: usize) -> Vec<Option<f32>> {
        use datafusion::arrow::array::Float32Array;
        let v = arr.value(row);
        let a = v.as_any().downcast_ref::<Float32Array>().unwrap();
        (0..a.len())
            .map(|i| if a.is_null(i) { None } else { Some(a.value(i)) })
            .collect()
    }
    fn list_row_str(arr: &ListArray, row: usize) -> Vec<String> {
        let v = arr.value(row);
        let a = v.as_any().downcast_ref::<StringArray>().unwrap();
        (0..a.len()).map(|i| a.value(i).to_string()).collect()
    }

    #[test]
    fn encode_af_single_allele_round_trips() {
        let col = StringArray::from(vec![Some("A:0.1|A:0.2|A:0")]);
        let af = encode_af_2array(&col, 3).unwrap();
        assert_eq!(list_row_str(&af.alleles, 0), vec!["A"]);
        assert_eq!(
            list_row_f32(&af.freqs, 0),
            vec![Some(0.1), Some(0.2), Some(0.0)]
        );
    }

    #[test]
    fn encode_af_multiallelic_is_allele_major() {
        let col = StringArray::from(vec![Some("A:0.1,G:0.2|A:0.3,G:0.4")]);
        let af = encode_af_2array(&col, 2).unwrap();
        assert_eq!(list_row_str(&af.alleles, 0), vec!["A", "G"]);
        // A over 2 pops, then G over 2 pops.
        assert_eq!(
            list_row_f32(&af.freqs, 0),
            vec![Some(0.1), Some(0.3), Some(0.2), Some(0.4)]
        );
    }

    #[test]
    fn encode_af_missing_population_is_null() {
        let col = StringArray::from(vec![Some("A:0.1||A:0.3")]);
        let af = encode_af_2array(&col, 3).unwrap();
        assert_eq!(list_row_str(&af.alleles, 0), vec!["A"]);
        assert_eq!(list_row_f32(&af.freqs, 0), vec![Some(0.1), None, Some(0.3)]);
    }

    #[test]
    fn encode_af_null_and_empty_rows_are_null_elements() {
        let col = StringArray::from(vec![None, Some("")]);
        let af = encode_af_2array(&col, 2).unwrap();
        assert!(af.alleles.is_null(0) && af.alleles.is_null(1));
        assert!(af.freqs.is_null(0) && af.freqs.is_null(1));
    }

    #[test]
    fn dedup_variation_name_nulls_exact_matches() {
        let vn = StringArray::from(vec![Some("rs1"), Some("rs2"), Some("x")]);
        let db = StringArray::from(vec![Some("rs1"), Some("rsZ"), Some("x")]);
        let out = dedup_variation_name(&vn, &db);
        assert!(out.is_null(0));
        assert_eq!(out.value(1), "rs2");
        assert!(out.is_null(2));
    }
}
