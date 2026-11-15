//! Per-read metrics that whole-file statistics (`seqkit stats -a`) are built from.

use super::KernelError;
use super::gc::gc_counts;
use super::quality::{ERROR_PROB, phred};

/// One read's contribution to whole-file statistics.
///
/// Summing these over a file and adding a length histogram gives every
/// `seqkit stats -a` column.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReadMetrics {
    /// Sequence length in bytes, gaps included.
    pub len: u64,
    /// Count of G and C in either case.
    pub gc: u64,
    /// Count of `N`/`n`, seqkit's `sum_n` for nucleotide input.
    pub n: u64,
    /// Count of `X`/`x`, seqkit's `sum_n` for protein input.
    pub x: u64,
    /// Count of the gap characters `-`, `.` and space (seqkit's default `-G`).
    pub gap: u64,
    /// Bases with a Phred score of at least 20.
    pub q20: u64,
    /// Bases with a Phred score of at least 30.
    pub q30: u64,
    /// Sum of Phred scores.
    pub qual_sum: u64,
    /// Sum of per-base error probabilities `10^(-Q/10)`.
    pub err_sum: f64,
}

/// Computes [`ReadMetrics`] for one read. `qual` is `None` for FASTA input, in
/// which case the quality fields are zero.
pub fn read_metrics(
    seq: &[u8],
    qual: Option<&[u8]>,
    offset: u8,
) -> Result<ReadMetrics, KernelError> {
    let gc = gc_counts(seq);
    let mut m = ReadMetrics {
        len: seq.len() as u64,
        gc: gc.gc as u64,
        ..Default::default()
    };
    for &b in seq {
        m.n += matches!(b, b'N' | b'n') as u64;
        m.x += matches!(b, b'X' | b'x') as u64;
        m.gap += matches!(b, b'-' | b'.' | b' ') as u64;
    }
    if let Some(qual) = qual {
        let table = &*ERROR_PROB;
        for (i, &b) in qual.iter().enumerate() {
            let q = phred(b, offset, i)?;
            m.q20 += (q >= 20) as u64;
            m.q30 += (q >= 30) as u64;
            m.qual_sum += q as u64;
            m.err_sum += table[q as usize];
        }
    }
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts() {
        let m = read_metrics(b"ACGN-. nX", Some(b"I5+#I5+#I"), 33).unwrap();
        assert_eq!(m.len, 9);
        assert_eq!(m.gc, 2);
        assert_eq!(m.n, 2);
        assert_eq!(m.x, 1);
        assert_eq!(m.gap, 3);
        // Q40, Q20, Q10, Q2, ...: >= 20 counts Q40 and Q20.
        assert_eq!(m.q20, 5);
        assert_eq!(m.q30, 3);
        assert_eq!(m.qual_sum, 40 * 3 + 20 * 2 + 10 * 2 + 2 * 2);
    }

    #[test]
    fn seqkit_q20_probe() {
        // 7 of 8 bases at Q20 or above: seqkit prints 88 (87.5 rounded).
        let m = read_metrics(b"ACGTACGT", Some(b"5555555#"), 33).unwrap();
        assert_eq!(m.q20, 7);
    }

    #[test]
    fn fasta_has_no_quality() {
        let m = read_metrics(b"ACGT", None, 33).unwrap();
        assert_eq!((m.q20, m.q30, m.qual_sum, m.err_sum), (0, 0, 0, 0.0));
    }

    #[test]
    fn quality_below_offset() {
        assert!(read_metrics(b"A", Some(b" "), 33).is_err());
    }
}
