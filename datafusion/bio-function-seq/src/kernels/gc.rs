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
    /// The full sequence length in bytes: seqkit `fx2tab -B GC` and `stats`,
    /// seqfu, seqtk `comp`.
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

/// Windows over a sequence. A window starts at every multiple of `step`
/// below the sequence length and is cut at the sequence end, so the last
/// windows can be shorter. With `step == window` this is the windowing of
/// Biopython's `GC_skew`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Windows {
    window: usize,
    step: usize,
}

impl Windows {
    /// `window` and `step` must be at least 1. `step` defaults to `window`,
    /// which gives windows that do not overlap.
    pub fn new(window: i64, step: Option<i64>) -> Result<Self, KernelError> {
        let positive = |name: &'static str, value: i64| {
            usize::try_from(value)
                .ok()
                .filter(|v| *v > 0)
                .ok_or_else(|| KernelError::InvalidOption {
                    name,
                    value: value.to_string(),
                    expected: "a positive integer",
                })
        };
        let window = positive("window", window)?;
        let step = match step {
            Some(step) => positive("step", step)?,
            None => window,
        };
        Ok(Self { window, step })
    }

    fn ranges(self, len: usize) -> impl Iterator<Item = std::ops::Range<usize>> {
        (0..len)
            .step_by(self.step)
            .map(move |start| start..start.saturating_add(self.window).min(len))
    }
}

/// Appends the GC content of every window of `seq`, each computed as
/// [`gc_content`] does: `None` for a window whose denominator is zero.
pub fn gc_content_windows_into(
    seq: &[u8],
    windows: Windows,
    denominator: GcDenominator,
    out: &mut Vec<Option<f64>>,
) {
    out.extend(
        windows
            .ranges(seq.len())
            .map(|range| gc_content(&seq[range], denominator)),
    );
}

/// Appends the GC skew `(G - C) / (G + C)` of every window of `seq`,
/// ignoring case. A window with no G or C gives 0, as in Biopython's
/// `GC_skew`.
pub fn gc_skew_into(seq: &[u8], windows: Windows, out: &mut Vec<f64>) {
    out.extend(windows.ranges(seq.len()).map(|range| {
        let (mut g, mut c) = (0i64, 0i64);
        for &b in &seq[range] {
            match b {
                b'G' | b'g' => g += 1,
                b'C' | b'c' => c += 1,
                _ => {}
            }
        }
        if g + c == 0 {
            0.0
        } else {
            (g - c) as f64 / (g + c) as f64
        }
    }));
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

#[cfg(test)]
mod window_tests {
    use super::*;

    fn skew(seq: &str, window: i64, step: Option<i64>) -> Vec<f64> {
        let mut out = Vec::new();
        gc_skew_into(
            seq.as_bytes(),
            Windows::new(window, step).unwrap(),
            &mut out,
        );
        out
    }

    fn content(seq: &str, window: i64, step: Option<i64>, d: GcDenominator) -> Vec<Option<f64>> {
        let mut out = Vec::new();
        gc_content_windows_into(
            seq.as_bytes(),
            Windows::new(window, step).unwrap(),
            d,
            &mut out,
        );
        out
    }

    #[test]
    fn biopython_gc_skew_probes() {
        // Bio.SeqUtils.GC_skew(seq, window), Biopython 1.88.
        assert_eq!(skew("GGGGCCCCGGGCAAAA", 4, None), [1.0, -1.0, 0.5, 0.0]);
        assert_eq!(skew("GGCC", 100, None), [0.0]);
        assert_eq!(skew("ggcc", 2, None), [1.0, -1.0]);
        // The last window is partial.
        assert_eq!(skew("AAAAGG", 4, None), [0.0, 1.0]);
        // S is not counted; a window with no G or C gives 0.
        assert_eq!(skew("GGNNSS", 3, None), [1.0, 0.0]);
        assert_eq!(skew("", 4, None), [] as [f64; 0]);
    }

    #[test]
    fn windows_start_at_every_multiple_of_the_step() {
        // Starts 0, 2, 4: GGCC, CCGG, GG (cut at the end).
        assert_eq!(skew("GGCCGG", 4, Some(2)), [0.0, 0.0, 1.0]);
        // A step larger than the window skips bases: GG, CC.
        assert_eq!(skew("GGAACC", 2, Some(4)), [1.0, -1.0]);
    }

    #[test]
    fn gc_content_per_window() {
        let length = GcDenominator::Length;
        assert_eq!(content("GGCCAAAT", 4, None, length), [Some(1.0), Some(0.0)]);
        // The last window is partial: AA.
        assert_eq!(content("GGCCAA", 4, None, length), [Some(1.0), Some(0.0)]);
        assert_eq!(content("GGNN", 4, None, length), [Some(0.5)]);
        // A zero denominator gives a null element.
        assert_eq!(
            content("GGNNNNNN", 4, None, GcDenominator::Acgt),
            [Some(1.0), None]
        );
        assert_eq!(content("", 4, None, length), []);
    }

    #[test]
    fn one_window_equals_gc_content() {
        for seq in ["GATTACA", "ggcc", "NNNN", "G"] {
            for d in [GcDenominator::Length, GcDenominator::Acgt] {
                assert_eq!(
                    content(seq, seq.len() as i64, None, d),
                    [gc_content(seq.as_bytes(), d)]
                );
            }
        }
    }

    #[test]
    fn invalid_windows() {
        for (window, step, name) in [
            (0, None, "window"),
            (-1, None, "window"),
            (4, Some(0), "step"),
        ] {
            let err = Windows::new(window, step).unwrap_err();
            assert!(
                matches!(err, KernelError::InvalidOption { name: n, .. } if n == name),
                "{err:?}"
            );
        }
    }
}
