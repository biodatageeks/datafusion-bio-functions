//! Mean Phred quality.

use std::sync::LazyLock;

use super::KernelError;

/// How per-base Phred scores are averaged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum QualityMethod {
    /// `-10·log10(mean(10^(-Q/10)))`, the mean of the error probabilities
    /// expressed as a Phred score. seqkit uses this on purpose (seqkit #448).
    #[default]
    ErrorProb,
    /// The plain mean of the Phred scores, as in FastQC.
    Arithmetic,
}

impl QualityMethod {
    /// Parses `"error_prob"` or `"arithmetic"`, ignoring case.
    pub fn parse(value: &str) -> Result<Self, KernelError> {
        match value.to_ascii_lowercase().as_str() {
            "error_prob" => Ok(QualityMethod::ErrorProb),
            "arithmetic" => Ok(QualityMethod::Arithmetic),
            _ => Err(KernelError::InvalidOption {
                name: "method",
                value: value.to_string(),
                expected: "'error_prob' or 'arithmetic'",
            }),
        }
    }
}

/// `10^(-Q/10)` for every possible Phred score.
pub(crate) static ERROR_PROB: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0; 256];
    for (q, p) in table.iter_mut().enumerate() {
        *p = 10f64.powf(q as f64 / -10.0);
    }
    table
});

/// Converts a quality byte to a Phred score.
#[inline]
pub(crate) fn phred(byte: u8, offset: u8, position: usize) -> Result<u8, KernelError> {
    byte.checked_sub(offset)
        .ok_or(KernelError::QualityBelowOffset {
            byte,
            offset,
            position,
        })
}

/// Mean quality of a read, or `None` for an empty quality string.
///
/// A byte below `offset` is an error rather than a silent negative score.
pub fn phred_mean(
    qual: &[u8],
    offset: u8,
    method: QualityMethod,
) -> Result<Option<f64>, KernelError> {
    if qual.is_empty() {
        return Ok(None);
    }
    let n = qual.len() as f64;
    match method {
        QualityMethod::Arithmetic => {
            let mut sum = 0u64;
            for (i, &b) in qual.iter().enumerate() {
                sum += phred(b, offset, i)? as u64;
            }
            Ok(Some(sum as f64 / n))
        }
        QualityMethod::ErrorProb => {
            let table = &*ERROR_PROB;
            let mut sum = 0.0;
            for (i, &b) in qual.iter().enumerate() {
                sum += table[phred(b, offset, i)? as usize];
            }
            Ok(Some(-10.0 * (sum / n).log10()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round2(v: f64) -> f64 {
        (v * 100.0).round() / 100.0
    }

    #[test]
    fn seqkit_probes() {
        // seqkit fx2tab -q
        for (qual, expected) in [("I#", 5.01), ("5I", 22.97), ("????", 30.0)] {
            let got = phred_mean(qual.as_bytes(), 33, QualityMethod::ErrorProb)
                .unwrap()
                .unwrap();
            assert_eq!(round2(got), expected, "{qual}");
        }
    }

    #[test]
    fn arithmetic_mean() {
        assert_eq!(
            phred_mean(b"I#", 33, QualityMethod::Arithmetic).unwrap(),
            Some(21.0)
        );
        // Phred+64
        assert_eq!(
            phred_mean(b"hh", 64, QualityMethod::Arithmetic).unwrap(),
            Some(40.0)
        );
    }

    #[test]
    fn empty_is_none() {
        assert_eq!(phred_mean(b"", 33, QualityMethod::ErrorProb).unwrap(), None);
    }

    #[test]
    fn below_offset_is_an_error() {
        let err = phred_mean(b"II ", 33, QualityMethod::ErrorProb).unwrap_err();
        assert_eq!(
            err,
            KernelError::QualityBelowOffset {
                byte: b' ',
                offset: 33,
                position: 2
            }
        );
        assert!(phred_mean(b"5", 64, QualityMethod::Arithmetic).is_err());
    }

    #[test]
    fn method_parse() {
        assert_eq!(
            QualityMethod::parse("Arithmetic").unwrap(),
            QualityMethod::Arithmetic
        );
        assert!(QualityMethod::parse("mean").is_err());
    }
}
