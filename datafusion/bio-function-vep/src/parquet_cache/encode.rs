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

use std::collections::HashMap;
use std::fmt::{self, Write as _};

use datafusion::arrow::array::{
    Array, BooleanArray, Float32Array, Float32Builder, Int8Array, ListArray, ListBuilder,
    StringArray, StringBuilder,
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

        // Master allele order: the population with the most (allele:freq) pairs,
        // PRESERVING duplicate allele strings and order. The read side
        // (`reconstruct_af_group_string`) is positional (allele-major, `idx =
        // a*n_pops + p`), so a single master order carries every population; each
        // population's own alleles are an in-order subsequence of it (VEP lists a
        // canonical allele order per variant and populations omit alleles they
        // lack data for). A first-seen DEDUP here (the old behavior) silently
        // dropped the 2nd+ freq of a repeated allele string (e.g. microsatellites
        // where distinct true alleles trim to an identical repeat string).
        let master: Vec<String> = pops
            .iter()
            .flatten()
            .max_by_key(|pairs| pairs.len())
            .map(|pairs| pairs.iter().map(|(a, _)| a.clone()).collect())
            .unwrap_or_default();

        for allele in &master {
            alleles.values().append_value(allele);
        }
        alleles.append(true);

        // Align each population's pairs to `master` as an in-order subsequence,
        // recording the freq at each master position (`None` where the population
        // has no entry there).
        let mut pop_rows: Vec<Vec<Option<f32>>> = Vec::with_capacity(n_pops);
        for pop_idx in 0..n_pops {
            let mut row = vec![None; master.len()];
            if let Some(Some(pairs)) = pops.get(pop_idx) {
                let mut j = 0usize;
                for (i, mname) in master.iter().enumerate() {
                    if j < pairs.len() && &pairs[j].0 == mname {
                        row[i] = Some(pairs[j].1);
                        j += 1;
                    }
                }
                if j != pairs.len() {
                    return Err(DataFusionError::Execution(format!(
                        "AF group population alleles are not an in-order subsequence \
                         of the master allele order (positional 2-array encoding \
                         cannot represent this row): '{s}'"
                    )));
                }
            }
            pop_rows.push(row);
        }

        // Allele-major emission: for each master allele position, one freq slot per
        // population.
        for pos in 0..master.len() {
            for row in &pop_rows {
                match row[pos] {
                    Some(f) => freqs.values().append_value(f),
                    None => freqs.values().append_null(),
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

/// Inverse of [`encode_af_2array`]: rebuild the pipe-joined AF group string
/// (`pop0|pop1|...`, each population `allele:freq[,allele:freq]`) from the
/// struct-of-arrays. `n_pops` is the group's population count (6/10/11). The
/// result is the scalar `Utf8` group column that `af_bundle::unbundle_af_columns`
/// splits into the per-population logical AF columns.
///
/// Frequencies are formatted with [`format_g4`] (the source is 4 significant
/// figures), which reproduces the original CSQ text byte-for-byte.
///
/// The read path no longer goes through this: [`reconstruct_af_members`] builds
/// the per-population columns directly. This stays as the readable definition of
/// the round trip, and as the oracle that function is tested against.
pub fn reconstruct_af_group_string(
    alleles: &ListArray,
    freqs: &ListArray,
    n_pops: usize,
) -> Result<StringArray> {
    let n = alleles.len();
    let mut b = StringBuilder::new();
    for r in 0..n {
        // A null list element means the whole group is absent for this row
        // (matches `concat_group`'s null-when-all-absent).
        if alleles.is_null(r) || freqs.is_null(r) {
            b.append_null();
            continue;
        }
        let al = alleles.value(r);
        let al = al.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
            DataFusionError::Execution("AF alleles list element must be Utf8".to_string())
        })?;
        let fr = freqs.value(r);
        let fr = fr.as_any().downcast_ref::<Float32Array>().ok_or_else(|| {
            DataFusionError::Execution("AF freqs list element must be Float32".to_string())
        })?;
        let n_alleles = al.len();
        let mut segments: Vec<String> = Vec::with_capacity(n_pops);
        for p in 0..n_pops {
            let mut parts: Vec<String> = Vec::new();
            for a in 0..n_alleles {
                let idx = a * n_pops + p;
                if idx < fr.len() && !fr.is_null(idx) {
                    parts.push(format!("{}:{}", al.value(a), format_g4(fr.value(idx))));
                }
            }
            segments.push(parts.join(","));
        }
        b.append_value(segments.join("|"));
    }
    Ok(b.finish())
}

/// Longest string the fast path of [`write_g4`] produces: `0.000dddd` and
/// `d.ddde-NN` are both nine bytes.
const G4_FAST_MAX: usize = 12;

/// A fixed-capacity `fmt::Write` sink on the stack, so formatting a frequency
/// allocates nothing.
struct StackBuf<const N: usize> {
    bytes: [u8; N],
    len: usize,
}

impl<const N: usize> StackBuf<N> {
    fn new() -> Self {
        Self {
            bytes: [0; N],
            len: 0,
        }
    }

    fn push(&mut self, byte: u8) {
        self.bytes[self.len] = byte;
        self.len += 1;
    }

    fn as_bytes(&self) -> &[u8] {
        &self.bytes[..self.len]
    }

    fn as_str(&self) -> &str {
        // Only ASCII digits and punctuation, or whole `&str`s, are ever written.
        std::str::from_utf8(self.as_bytes()).expect("StackBuf holds UTF-8")
    }
}

impl<const N: usize> fmt::Write for StackBuf<N> {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        let end = self.len + s.len();
        if end > N {
            return Err(fmt::Error);
        }
        self.bytes[self.len..end].copy_from_slice(s.as_bytes());
        self.len = end;
        Ok(())
    }
}

/// [`format_g4`] without the heap: byte-identical output, written straight into
/// `out`.
///
/// `format_g4` formats the value twice and allocates three or four `String`s on
/// the way. This makes the same single `{:.3e}` call into a stack buffer -- so
/// the rounding is the standard library's in both -- reads the decimal exponent
/// off the bytes, and places the four significant digits by hand. Anything
/// outside the domain an allele frequency lives in (not finite, negative, or a
/// decimal exponent of 4 and up) is handed to `format_g4` itself, which keeps
/// its behaviour there, including the panic on a non-finite value.
///
/// `format_g4` stays the definition: the cache builder rejects any frequency it
/// cannot reproduce (`encode_af_2array`), so this must equal it, never replace
/// it. The tests compare the two over every four-digit decimal and its
/// neighbouring floats.
fn write_g4<W: fmt::Write>(out: &mut W, f: f32) -> fmt::Result {
    if f == 0.0 {
        return out.write_char('0');
    }
    let mut sci = StackBuf::<16>::new();
    if !(f.is_finite() && f > 0.0) || write!(sci, "{f:.3e}").is_err() {
        return out.write_str(&format_g4(f));
    }
    // `d.ddde[-]N`
    let s = sci.as_bytes();
    if s.len() < 7 || s[1] != b'.' || s[5] != b'e' {
        return out.write_str(&format_g4(f));
    }
    let (negative, exp_digits) = match s[6] {
        b'-' => (true, &s[7..]),
        _ => (false, &s[6..]),
    };
    if exp_digits.is_empty() || exp_digits.len() > 2 || !exp_digits.iter().all(u8::is_ascii_digit) {
        return out.write_str(&format_g4(f));
    }
    let magnitude = exp_digits
        .iter()
        .fold(0i32, |acc, d| acc * 10 + i32::from(d - b'0'));
    let e = if negative { -magnitude } else { magnitude };
    if e >= 4 {
        return out.write_str(&format_g4(f));
    }

    let digits = [s[0], s[2], s[3], s[4]];
    // `%g` drops trailing zeros; at least one significant digit always stays.
    let mut kept = 4;
    while kept > 1 && digits[kept - 1] == b'0' {
        kept -= 1;
    }

    let mut text = StackBuf::<G4_FAST_MAX>::new();
    if e < -4 {
        // Scientific: `d[.ddd]e-NN`. The exponent is negative here by
        // construction, and an f32 never needs a third exponent digit.
        text.push(digits[0]);
        if kept > 1 {
            text.push(b'.');
            for &d in &digits[1..kept] {
                text.push(d);
            }
        }
        text.push(b'e');
        text.push(b'-');
        text.push(b'0' + (magnitude / 10) as u8);
        text.push(b'0' + (magnitude % 10) as u8);
    } else if e < 0 {
        // `0.` + leading zeros + the significant digits.
        text.push(b'0');
        text.push(b'.');
        for _ in 0..(-e - 1) {
            text.push(b'0');
        }
        for &d in &digits[..kept] {
            text.push(d);
        }
    } else {
        // `e + 1` integer digits (never trimmed: `1200` keeps its zeros), then
        // whatever significant digits are left after the point.
        let int_digits = (e + 1) as usize;
        for &d in &digits[..int_digits] {
            text.push(d);
        }
        if kept > int_digits {
            text.push(b'.');
            for &d in &digits[int_digits..kept] {
                text.push(d);
            }
        }
    }
    out.write_str(text.as_str())
}

/// Formatted frequencies seen so far, keyed on the float's bits.
///
/// A shard stores frequencies to four significant figures, so the values repeat
/// heavily: on HG002 chr1 the rebuild formats 19.7 M frequencies but only ~53 k
/// distinct ones. Scoped to one rebuild call on purpose -- the lookup that owns
/// the rebuild is shared across workers, and a per-call map needs neither a
/// lock nor a size bound.
struct G4Memo {
    seen: HashMap<u32, (u8, [u8; G4_FAST_MAX])>,
}

impl G4Memo {
    fn new() -> Self {
        Self {
            seen: HashMap::new(),
        }
    }

    fn write<W: fmt::Write>(&mut self, out: &mut W, f: f32) -> fmt::Result {
        let bits = f.to_bits();
        if let Some((len, bytes)) = self.seen.get(&bits) {
            let text = std::str::from_utf8(&bytes[..usize::from(*len)]).expect("memo holds UTF-8");
            return out.write_str(text);
        }
        let mut text = StackBuf::<48>::new();
        write_g4(&mut text, f)?;
        if text.len <= G4_FAST_MAX {
            let mut bytes = [0u8; G4_FAST_MAX];
            bytes[..text.len].copy_from_slice(text.as_bytes());
            self.seen.insert(bits, (text.len as u8, bytes));
        }
        out.write_str(text.as_str())
    }
}

/// Rebuild one AF group straight into its per-population `Utf8` columns.
///
/// Equal, cell for cell, to `split_group(&reconstruct_af_group_string(..)?, n_pops)`
/// -- the pipe-joined group string is never built, because nothing downstream
/// reads it: `unbundle_af_columns` only ever split it back apart. Every cell is
/// non-null; an absent group, population or allele contributes `""`, exactly as
/// `split_group` yields.
///
/// Reads the two child arrays through the list offsets instead of slicing a
/// fresh array per row. The `idx < row_end` bound is what keeps a freqs list
/// shorter than `n_alleles * n_pops` from reading the next row's values.
pub fn reconstruct_af_members(
    alleles: &ListArray,
    freqs: &ListArray,
    n_pops: usize,
) -> Result<Vec<StringArray>> {
    let n = alleles.len();
    let allele_values = alleles.values().as_any().downcast_ref::<StringArray>();
    let freq_values = freqs.values().as_any().downcast_ref::<Float32Array>();
    let (allele_offsets, freq_offsets) = (alleles.value_offsets(), freqs.value_offsets());

    let mut members: Vec<StringBuilder> = (0..n_pops)
        .map(|_| StringBuilder::with_capacity(n, n * 8))
        .collect();
    let mut memo = G4Memo::new();
    let write_failed =
        |_: fmt::Error| DataFusionError::Execution("failed to write an AF value".to_string());

    for r in 0..n {
        // A null list slot may still span values, so test the row first.
        if alleles.is_null(r) || freqs.is_null(r) {
            for member in &mut members {
                member.append_value("");
            }
            continue;
        }
        let al = allele_values.ok_or_else(|| {
            DataFusionError::Execution("AF alleles list element must be Utf8".to_string())
        })?;
        let fr = freq_values.ok_or_else(|| {
            DataFusionError::Execution("AF freqs list element must be Float32".to_string())
        })?;
        let (al_start, al_end) = (allele_offsets[r] as usize, allele_offsets[r + 1] as usize);
        let (fr_start, fr_end) = (freq_offsets[r] as usize, freq_offsets[r + 1] as usize);

        for (p, member) in members.iter_mut().enumerate() {
            let mut first = true;
            for a in 0..(al_end - al_start) {
                let idx = fr_start + a * n_pops + p;
                if idx < fr_end && !fr.is_null(idx) {
                    if !first {
                        member.write_char(',').map_err(write_failed)?;
                    }
                    first = false;
                    member
                        .write_str(al.value(al_start + a))
                        .map_err(write_failed)?;
                    member.write_char(':').map_err(write_failed)?;
                    memo.write(member, fr.value(idx)).map_err(write_failed)?;
                }
            }
            member.append_value("");
        }
    }
    Ok(members.into_iter().map(|mut b| b.finish()).collect())
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
    fn encode_af_preserves_duplicate_allele_strings() {
        // Microsatellite case: the same trimmed allele string ("X") appears twice
        // with different freqs. The old first-seen dedup dropped the 2nd, losing
        // the 0.2712 freq (see the parquet-af-duplicate-allele-encoder-bug). Both
        // populations share the same allele sequence [X, Y, X].
        let col = StringArray::from(vec![Some(
            "X:0.004792,Y:0.02556,X:0.2712|X:0,Y:0.0106,X:0.4319",
        )]);
        let af = encode_af_2array(&col, 2).unwrap();
        assert_eq!(list_row_str(&af.alleles, 0), vec!["X", "Y", "X"]);
        // allele-major over 2 pops: X, X, Y, Y, X', X'.
        assert_eq!(
            list_row_f32(&af.freqs, 0),
            vec![
                Some(0.004792),
                Some(0.0),
                Some(0.02556),
                Some(0.0106),
                Some(0.2712),
                Some(0.4319),
            ]
        );
        // Round-trips back to the exact source string (the losslessness the
        // downstream CSQ text relies on).
        let s = reconstruct_af_group_string(&af.alleles, &af.freqs, 2).unwrap();
        assert_eq!(
            s.value(0),
            "X:0.004792,Y:0.02556,X:0.2712|X:0,Y:0.0106,X:0.4319"
        );
    }

    #[test]
    fn encode_af_rejects_divergent_population_allele_order() {
        // The positional 2-array encoding requires a shared allele sequence across
        // populations; a divergent order must fail loud rather than corrupt.
        let col = StringArray::from(vec![Some("A:0.1,G:0.2|G:0.3,A:0.4")]);
        assert!(encode_af_2array(&col, 2).is_err());
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

    #[test]
    fn af_2array_reconstruct_round_trips() {
        // Single-allele, multi-allelic, a missing population, and a null row —
        // encode then reconstruct must reproduce the original group string.
        let inputs = vec![
            Some("A:0.08806|A:0.0367|A:2.682e-05"), // 3 pops, single allele
            Some("A:0.1,G:0.2|A:0.3,G:0.4|A:0.5,G:0.6"), // multi-allelic
            Some("A:0.1||A:0.3"),                   // middle pop missing
            None,                                   // null row
        ];
        let col = StringArray::from(inputs.clone());
        let af = encode_af_2array(&col, 3).unwrap();
        let back = reconstruct_af_group_string(&af.alleles, &af.freqs, 3).unwrap();
        assert_eq!(back.len(), inputs.len());
        for (i, want) in inputs.iter().enumerate() {
            match want {
                Some(s) => assert_eq!(back.value(i), *s, "row {i}"),
                None => assert!(back.is_null(i), "row {i} should be null"),
            }
        }
    }

    fn g4(f: f32) -> String {
        let mut s = String::new();
        write_g4(&mut s, f).unwrap();
        s
    }

    fn xorshift(state: &mut u32) -> u32 {
        *state ^= *state << 13;
        *state ^= *state >> 17;
        *state ^= *state << 5;
        *state
    }

    #[test]
    fn write_g4_matches_the_measured_tokens() {
        for (f, expected) in [
            (0.08806f32, "0.08806"),
            (0.0, "0"),
            (-0.0, "0"),
            (1.0, "1"),
            (0.5, "0.5"),
            (2.682e-05, "2.682e-05"),
            (4.248e-06, "4.248e-06"),
            (1e-05, "1e-05"),
            (0.0001, "0.0001"),
            (9.9996e-3, "0.01"),
            (0.99996, "1"),
            (9.9996e-5, "0.0001"),
            (12.5, "12.5"),
            (1200.0, "1200"),
            (9999.4, "9999"),
        ] {
            assert_eq!(g4(f), expected, "write_g4({f:e})");
            assert_eq!(format_g4(f), expected, "format_g4({f:e})");
        }
    }

    /// Every four-significant-digit decimal an AF can be, plus the float on
    /// either side of it, plus a fixed sample of arbitrary finite floats. The
    /// neighbours are what reach the rounding carries (`9.9996e-3` -> `0.01`)
    /// and the `e == -4` / `e == -5` and `e == 3` / `e == 4` switches.
    #[test]
    fn write_g4_equals_format_g4_over_the_af_domain() {
        let mut checked = 0u32;
        for e in -12..=4 {
            for m in 1000..=9999u32 {
                let base: f32 = format!("{}.{:03}e{e}", m / 1000, m % 1000).parse().unwrap();
                for step in [-1i64, 0, 1] {
                    let f = f32::from_bits((i64::from(base.to_bits()) + step) as u32);
                    assert_eq!(g4(f), format_g4(f), "{f:e} (bits {:#x})", f.to_bits());
                    checked += 1;
                }
            }
        }
        let mut state = 0x9E37_79B9u32;
        for _ in 0..2_000_000 {
            // Any finite float, either sign, subnormals included.
            let f = f32::from_bits(xorshift(&mut state));
            if f.is_finite() {
                assert_eq!(g4(f), format_g4(f), "{f:e} (bits {:#x})", f.to_bits());
                checked += 1;
            }
        }
        assert!(checked > 2_400_000);
    }

    /// All 4.28 billion finite floats. Minutes in release, so not part of the
    /// default run: `cargo test --release -- --ignored write_g4_equals_format_g4_exhaustively`.
    #[test]
    #[ignore = "exhaustive f32 sweep; run in release"]
    fn write_g4_equals_format_g4_exhaustively() {
        for bits in 0..=u32::MAX {
            let f = f32::from_bits(bits);
            if f.is_finite() {
                assert_eq!(g4(f), format_g4(f), "bits {bits:#x}");
            }
        }
    }

    /// `format_g4` panics on a non-finite value and the cache builder calls it
    /// on every frequency, so none can be stored. The fast path must not turn
    /// that into a quietly different string.
    #[test]
    fn write_g4_keeps_format_g4s_panic_on_non_finite_values() {
        for f in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            assert!(std::panic::catch_unwind(|| format_g4(f)).is_err());
            assert!(std::panic::catch_unwind(|| g4(f)).is_err());
        }
    }

    const ALLELE_POOL: [&str; 7] = ["A", "C", "G", "T", "-", "AAAAAAAAAA", "TG"];
    const FREQ_POOL: [f32; 14] = [
        0.0, 1.0, 0.5, 0.08806, 2.682e-05, 9.911e-05, 4.248e-06, 0.0003994, 0.9969, 0.006912,
        9.9996e-3, 0.99996, 1234.0, 12.5,
    ];

    /// Rows of every shape the reader has to survive, including the ones the
    /// writer never produces: a null slot that still spans values, one list null
    /// and the other not, an empty list, and a freqs list shorter or longer than
    /// `n_alleles * n_pops`.
    fn generated_af_lists(seed: u32, n_rows: usize, n_pops: usize) -> (ListArray, ListArray) {
        let mut state = seed;
        let mut alleles = ListBuilder::new(StringBuilder::new());
        let mut freqs = ListBuilder::new(Float32Builder::new());
        for _ in 0..n_rows {
            let shape = xorshift(&mut state) % 10;
            let n_alleles = if shape == 3 {
                0
            } else {
                1 + (xorshift(&mut state) % 3) as usize
            };
            let n_freqs = match shape {
                4 => (n_alleles * n_pops).saturating_sub(1 + (xorshift(&mut state) % 4) as usize),
                5 => n_alleles * n_pops + 1 + (xorshift(&mut state) % 3) as usize,
                _ => n_alleles * n_pops,
            };
            for _ in 0..n_alleles {
                let allele = ALLELE_POOL[(xorshift(&mut state) as usize) % ALLELE_POOL.len()];
                alleles.values().append_value(allele);
            }
            for _ in 0..n_freqs {
                if xorshift(&mut state).is_multiple_of(4) {
                    freqs.values().append_null();
                } else {
                    let f = FREQ_POOL[(xorshift(&mut state) as usize) % FREQ_POOL.len()];
                    freqs.values().append_value(f);
                }
            }
            alleles.append(!matches!(shape, 0 | 1));
            freqs.append(!matches!(shape, 0 | 2));
        }
        (alleles.finish(), freqs.finish())
    }

    fn assert_members_match_the_group_string_path(a: &ListArray, f: &ListArray, n_pops: usize) {
        let joined = reconstruct_af_group_string(a, f, n_pops).unwrap();
        let expected = crate::cache::af_bundle::split_group(&joined, n_pops);
        let actual = reconstruct_af_members(a, f, n_pops).unwrap();
        assert_eq!(actual.len(), n_pops);
        for (p, (want, got)) in expected.iter().zip(&actual).enumerate() {
            assert_eq!(got.null_count(), 0, "population {p}");
            assert_eq!(got.len(), a.len(), "population {p}");
            for r in 0..a.len() {
                assert_eq!(got.value(r), want.value(r), "population {p} row {r}");
            }
        }
    }

    #[test]
    fn af_members_equal_the_group_string_path_on_generated_rows() {
        for (seed, n_pops) in [(1u32, 6usize), (2, 10), (3, 11), (4, 1)] {
            let (a, f) = generated_af_lists(seed, 2_000, n_pops);
            assert_members_match_the_group_string_path(&a, &f, n_pops);
        }
        // Zero rows.
        let (a, f) = generated_af_lists(9, 0, 6);
        assert_members_match_the_group_string_path(&a, &f, 6);
    }

    /// The lookup hands over a `concat_batches` result, and a slice keeps the
    /// parent's child array with non-zero offsets -- the two ways direct child
    /// indexing goes wrong if it assumes a row starts at 0.
    #[test]
    fn af_members_equal_the_group_string_path_on_sliced_and_concatenated_lists() {
        use datafusion::arrow::compute::concat;
        let n_pops = 11;
        let (a, f) = generated_af_lists(5, 500, n_pops);
        let (a_mid, f_mid) = (a.slice(37, 401), f.slice(37, 401));
        assert_members_match_the_group_string_path(&a_mid, &f_mid, n_pops);

        let (a2, f2) = generated_af_lists(6, 300, n_pops);
        let a_cat = concat(&[&a_mid, &a2]).unwrap();
        let f_cat = concat(&[&f_mid, &f2]).unwrap();
        assert_members_match_the_group_string_path(
            a_cat.as_any().downcast_ref::<ListArray>().unwrap(),
            f_cat.as_any().downcast_ref::<ListArray>().unwrap(),
            n_pops,
        );
    }

    #[test]
    fn af_members_reject_the_same_wrong_child_types() {
        use datafusion::arrow::array::Int32Builder;
        let mut wrong = ListBuilder::new(Int32Builder::new());
        wrong.values().append_value(1);
        wrong.append(true);
        let wrong = wrong.finish();
        let (a, f) = generated_af_lists(7, 1, 1);
        // Row 0 of this seed is a normal row, so both paths reach the downcast.
        assert!(!a.is_null(0) && !f.is_null(0));
        assert!(reconstruct_af_group_string(&wrong, &f, 1).is_err());
        assert!(reconstruct_af_members(&wrong, &f, 1).is_err());
        assert!(reconstruct_af_group_string(&a, &wrong, 1).is_err());
        assert!(reconstruct_af_members(&a, &wrong, 1).is_err());
    }
}
