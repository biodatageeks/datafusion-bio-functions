//! Quality trimming, after Trimmomatic's `LEADING`, `TRAILING` and
//! `SLIDINGWINDOW` steps.

use super::KernelError;

/// Trims a read by quality with Trimmomatic's `LEADING`, `TRAILING` and
/// `SLIDINGWINDOW` rules, applied in that order. Verified against
/// Trimmomatic 0.39.
///
/// As in Trimmomatic, an `N` base (uppercase) counts as quality 0 whatever
/// its quality character.
#[derive(Debug, Clone, Copy)]
pub struct QualityTrimmer {
    leading: Option<i64>,
    trailing: Option<i64>,
    window: Option<(usize, f32)>,
    offset: u8,
}

/// The span of a read that Trimmomatic would drop.
const DROPPED: (usize, usize) = (0, 0);

impl QualityTrimmer {
    /// - `leading`: remove bases from the start until one has at least this
    ///   quality.
    /// - `trailing`: remove bases from the end until one has at least this
    ///   quality.
    /// - `window`: `(size, quality)`; cut the read where the mean quality of
    ///   a window of `size` bases falls below `quality`.
    ///
    /// At least one step is required. `offset` is the Phred offset of the
    /// quality string.
    pub fn new(
        leading: Option<i64>,
        trailing: Option<i64>,
        window: Option<(i64, f64)>,
        offset: u8,
    ) -> Result<Self, KernelError> {
        if leading.is_none() && trailing.is_none() && window.is_none() {
            return Err(KernelError::InvalidOption {
                name: "steps",
                value: "none".to_string(),
                expected: "at least one of leading, trailing and sliding_window",
            });
        }
        let window = match window {
            Some((size, quality)) => {
                let size = usize::try_from(size)
                    .ok()
                    .filter(|size| *size > 0)
                    .ok_or_else(|| KernelError::InvalidOption {
                        name: "sliding_window",
                        value: size.to_string(),
                        expected: "a window size of at least 1",
                    })?;
                // Trimmomatic parses the quality as a 32-bit float.
                Some((size, quality as f32))
            }
            None => None,
        };
        Ok(Self {
            leading,
            trailing,
            window,
            offset,
        })
    }

    /// The byte range `start..end` of the read that is kept. `(0, 0)` when
    /// nothing is kept, which is when Trimmomatic drops the read.
    ///
    /// Two Trimmomatic details: `TRAILING` never examines the first base, so
    /// a read whose only passing base is its first is dropped; and
    /// `SLIDINGWINDOW` drops a read shorter than the window, and after
    /// cutting at a failing window also removes the bases below the window
    /// quality from the end of what is left, keeping at least one base.
    pub fn span(&self, seq: &[u8], qual: &[u8]) -> Result<(usize, usize), KernelError> {
        if seq.len() != qual.len() {
            return Err(KernelError::LengthMismatch {
                sequence: seq.len(),
                quality: qual.len(),
            });
        }
        if let Some(position) = qual.iter().position(|&b| b < self.offset) {
            return Err(KernelError::QualityBelowOffset {
                byte: qual[position],
                offset: self.offset,
                position,
            });
        }
        let score = |i: usize| -> i64 {
            if seq[i] == b'N' {
                0
            } else {
                i64::from(qual[i] - self.offset)
            }
        };
        let (mut start, mut end) = (0, qual.len());
        if let Some(quality) = self.leading {
            match (start..end).find(|&i| score(i) >= quality) {
                Some(i) => start = i,
                None => return Ok(DROPPED),
            }
        }
        if let Some(quality) = self.trailing {
            match (start + 1..end).rev().find(|&i| score(i) >= quality) {
                Some(i) => end = i + 1,
                None => return Ok(DROPPED),
            }
        }
        if let Some((size, quality)) = self.window {
            let len = end - start;
            if len < size {
                return Ok(DROPPED);
            }
            let required = quality * size as f32;
            let mut total: i64 = (start..start + size).map(score).sum();
            if (total as f32) < required {
                return Ok(DROPPED);
            }
            let mut keep = len;
            for i in 0..len - size {
                total += score(start + i + size) - score(start + i);
                if (total as f32) < required {
                    keep = i + size;
                    break;
                }
            }
            while keep > 1 && (score(start + keep - 1) as f32) < quality {
                keep -= 1;
            }
            end = start + keep;
        }
        Ok((start, end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Trimmed = Option<(&'static str, &'static str)>;

    /// Reads probed with Trimmomatic 0.39 (`SE -phred33`).
    const READS: [(&str, &str); 11] = [
        ("ACGTA", "I!!!!"),               // r1: only the first base is good
        ("ACGTAC", "!!II!!"),             // r2
        ("ACGTACGTACGT", "IIII!!!!IIII"), // r3
        ("ACG", "III"),                   // r4: shorter than the window
        ("ANNNA", "IIIII"),               // r5: N counts as quality 0
        ("NNACG", "IIIII"),               // r6
        ("ACGTA", "+++++"),               // r7: equal to a threshold of 10
        ("A", "I"),                       // r8: a single base
        ("ACGNN", "IIIII"),               // r9
        ("ACGTACGT", "!!!!IIII"),         // r10
        ("nnACG", "IIIII"),               // r11: lowercase n keeps its quality
    ];

    fn trim(
        trimmer: &QualityTrimmer,
        seq: &'static str,
        qual: &'static str,
    ) -> Option<(&'static str, &'static str)> {
        let (start, end) = trimmer.span(seq.as_bytes(), qual.as_bytes()).unwrap();
        (start != end).then(|| (&seq[start..end], &qual[start..end]))
    }

    #[track_caller]
    fn check(
        leading: Option<i64>,
        trailing: Option<i64>,
        window: Option<(i64, f64)>,
        expected: [Trimmed; 11],
    ) {
        let trimmer = QualityTrimmer::new(leading, trailing, window, 33).unwrap();
        for (i, ((seq, qual), expected)) in READS.into_iter().zip(expected).enumerate() {
            assert_eq!(trim(&trimmer, seq, qual), expected, "r{}", i + 1);
        }
    }

    const ALL: [Trimmed; 11] = [
        Some(("ACGTA", "I!!!!")),
        Some(("ACGTAC", "!!II!!")),
        Some(("ACGTACGTACGT", "IIII!!!!IIII")),
        Some(("ACG", "III")),
        Some(("ANNNA", "IIIII")),
        Some(("NNACG", "IIIII")),
        Some(("ACGTA", "+++++")),
        Some(("A", "I")),
        Some(("ACGNN", "IIIII")),
        Some(("ACGTACGT", "!!!!IIII")),
        Some(("nnACG", "IIIII")),
    ];

    #[test]
    fn leading() {
        // LEADING:10
        let mut expected = ALL;
        expected[1] = Some(("GTAC", "II!!"));
        expected[5] = Some(("ACG", "III"));
        expected[9] = Some(("ACGT", "IIII"));
        check(Some(10), None, None, expected);
    }

    #[test]
    fn trailing() {
        // TRAILING:10. Trimmomatic never checks the first base, so a read
        // whose only good base is the first one is dropped (r1, r8).
        let mut expected = ALL;
        expected[0] = None;
        expected[1] = Some(("ACGT", "!!II"));
        expected[7] = None;
        expected[8] = Some(("ACG", "III"));
        check(None, Some(10), None, expected);
    }

    #[test]
    fn sliding_window() {
        // SLIDINGWINDOW:4:15
        let mut expected = [None; 11];
        expected[1] = Some(("ACGT", "!!II"));
        expected[2] = Some(("ACGT", "IIII"));
        expected[5] = ALL[5];
        expected[8] = Some(("ACG", "III"));
        expected[10] = ALL[10];
        check(None, None, Some((4, 15.0)), expected);

        // SLIDINGWINDOW:4:10: a window equal to the threshold passes.
        let mut expected = expected;
        expected[0] = Some(("A", "I"));
        expected[4] = ALL[4];
        expected[6] = ALL[6];
        check(None, None, Some((4, 10.0)), expected);
    }

    #[test]
    fn steps_run_in_order() {
        // LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15
        let mut expected = [None; 11];
        expected[2] = Some(("ACGT", "IIII"));
        expected[9] = Some(("ACGT", "IIII"));
        expected[10] = ALL[10];
        check(Some(10), Some(10), Some((4, 15.0)), expected);
    }

    #[test]
    fn span_is_a_range_of_the_read() {
        let trimmer = QualityTrimmer::new(Some(10), Some(10), Some((4, 15.0)), 33).unwrap();
        assert_eq!(trimmer.span(b"ACGTACGT", b"!!!!IIII").unwrap(), (4, 8));
        assert_eq!(trimmer.span(b"", b"").unwrap(), (0, 0));
    }

    #[test]
    fn phred64() {
        // h is Q40 and B is Q2 with offset 64.
        let trimmer = QualityTrimmer::new(Some(10), None, None, 64).unwrap();
        assert_eq!(trimmer.span(b"ACGT", b"BBhh").unwrap(), (2, 4));
    }

    #[test]
    fn errors() {
        let trimmer = QualityTrimmer::new(Some(10), None, None, 33).unwrap();
        assert_eq!(
            trimmer.span(b"ACG", b"II").unwrap_err(),
            KernelError::LengthMismatch {
                sequence: 3,
                quality: 2
            }
        );
        assert!(matches!(
            trimmer.span(b"ACG", b"I I").unwrap_err(),
            KernelError::QualityBelowOffset { position: 1, .. }
        ));
        assert!(matches!(
            QualityTrimmer::new(None, None, None, 33).unwrap_err(),
            KernelError::InvalidOption { .. }
        ));
        assert!(matches!(
            QualityTrimmer::new(None, None, Some((0, 15.0)), 33).unwrap_err(),
            KernelError::InvalidOption {
                name: "sliding_window",
                ..
            }
        ));
    }
}
