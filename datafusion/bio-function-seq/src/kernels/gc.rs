//! GC content.

use super::KernelError;

/// Base counts needed for GC content.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GcCounts {
    /// Sequence length in bytes, gaps and unknown bytes included.
    pub len: usize,
    /// Count of `G`, `C`, `g` and `c`. `S` is not counted, like every oracle tool.
    pub gc: usize,
    /// Count of `A`, `C`, `G` and `T` in either case.
    pub acgt: usize,
}

/// What GC content is divided by.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum GcDenominator {
    /// The full sequence length: seqkit `fx2tab -B GC` and `stats`, seqfu,
    /// seqtk `comp`.
    #[default]
    Length,
    /// Only A, C, G and T: seqkit `fx2tab -g`.
    Acgt,
}

impl GcDenominator {
    /// Parses `"length"` or `"acgt"`, ignoring case.
    pub fn parse(value: &str) -> Result<Self, KernelError> {
        match value.to_ascii_lowercase().as_str() {
            "length" => Ok(GcDenominator::Length),
            "acgt" => Ok(GcDenominator::Acgt),
            _ => Err(KernelError::InvalidOption {
                name: "denominator",
                value: value.to_string(),
                expected: "'length' or 'acgt'",
            }),
        }
    }
}

const OTHER: u8 = 0;
const AT: u8 = 1;
const GC: u8 = 2;

const fn class_table() -> [u8; 256] {
    let mut table = [OTHER; 256];
    table[b'A' as usize] = AT;
    table[b'a' as usize] = AT;
    table[b'T' as usize] = AT;
    table[b't' as usize] = AT;
    table[b'G' as usize] = GC;
    table[b'g' as usize] = GC;
    table[b'C' as usize] = GC;
    table[b'c' as usize] = GC;
    table
}

static CLASS: [u8; 256] = class_table();

/// Counts the bases that GC content is computed from.
pub fn gc_counts(seq: &[u8]) -> GcCounts {
    let mut at = 0usize;
    let mut gc = 0usize;
    for &b in seq {
        let class = CLASS[b as usize];
        at += (class == AT) as usize;
        gc += (class == GC) as usize;
    }
    GcCounts {
        len: seq.len(),
        gc,
        acgt: at + gc,
    }
}

/// GC fraction between 0 and 1, or `None` when the denominator is zero.
pub fn gc_content(seq: &[u8], denominator: GcDenominator) -> Option<f64> {
    let counts = gc_counts(seq);
    let total = match denominator {
        GcDenominator::Length => counts.len,
        GcDenominator::Acgt => counts.acgt,
    };
    (total > 0).then(|| counts.gc as f64 / total as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pct(seq: &str, d: GcDenominator) -> Option<f64> {
        gc_content(seq.as_bytes(), d).map(|v| v * 100.0)
    }

    #[test]
    fn seqkit_probes() {
        // (input, seqkit fx2tab -g, seqkit fx2tab -B GC)
        for (seq, acgt, length) in [
            ("GGNN", 100.0, 50.0),
            ("GG--", 100.0, 50.0),
            ("GGRR", 100.0, 50.0),
            ("SSAA", 0.0, 0.0),
            ("ggaa", 50.0, 50.0),
        ] {
            assert_eq!(pct(seq, GcDenominator::Acgt), Some(acgt), "{seq} acgt");
            assert_eq!(
                pct(seq, GcDenominator::Length),
                Some(length),
                "{seq} length"
            );
        }
    }

    #[test]
    fn zero_denominator_is_none() {
        assert_eq!(gc_content(b"", GcDenominator::Length), None);
        assert_eq!(gc_content(b"", GcDenominator::Acgt), None);
        assert_eq!(gc_content(b"NNN", GcDenominator::Acgt), None);
        assert_eq!(gc_content(b"NNN", GcDenominator::Length), Some(0.0));
    }

    #[test]
    fn uracil_is_not_acgt() {
        assert_eq!(
            gc_counts(b"ACGU"),
            GcCounts {
                len: 4,
                gc: 2,
                acgt: 3
            }
        );
    }

    #[test]
    fn denominator_parse() {
        assert_eq!(GcDenominator::parse("ACGT").unwrap(), GcDenominator::Acgt);
        assert!(GcDenominator::parse("ignore_n").is_err());
    }
}
