//! IUPAC motif search, after Biopython's `Bio.SeqUtils.nt_search`.

use super::KernelError;
use super::translate::{A, C, G, T, base_set};

/// The base bit of an unambiguous sequence byte (`ACGTU` in either case,
/// with `U` read as `T`); 0 for anything else, which matches no pattern.
const fn sequence_bits() -> [u8; 256] {
    let mut m = [0u8; 256];
    let bases: [(u8, u8); 5] = [(b'A', A), (b'C', C), (b'G', G), (b'T', T), (b'U', T)];
    let mut i = 0;
    while i < bases.len() {
        let (base, bit) = bases[i];
        m[base as usize] = bit;
        m[base.to_ascii_lowercase() as usize] = bit;
        i += 1;
    }
    m
}

static SEQUENCE_BITS: [u8; 256] = sequence_bits();

/// A compiled nucleotide motif whose positions may be IUPAC ambiguity codes.
///
/// Ambiguity applies to the pattern only: a pattern code matches a sequence
/// base that is an unambiguous `A`, `C`, `G`, `T` or `U` in the code's set.
/// An ambiguous or unknown base in the sequence matches nothing, as in
/// Biopython's `nt_search`. Unlike `nt_search`, matching ignores the case of
/// the sequence and reads `U` as `T`.
#[derive(Debug, Clone)]
pub struct Motif {
    forward: Vec<u8>,
    /// The reverse complement, when both strands are searched and it differs
    /// from `forward`.
    reverse: Option<Vec<u8>>,
    overlapping: bool,
}

impl Motif {
    /// `pattern` is a non-empty string of IUPAC nucleotide codes in either
    /// case. With `both_strands`, the reverse complement of the pattern is
    /// searched too. Without `overlapping`, the search resumes after the end
    /// of each hit; the two orientations are then scanned independently.
    pub fn new(pattern: &str, both_strands: bool, overlapping: bool) -> Result<Self, KernelError> {
        let forward: Vec<u8> = pattern.bytes().map(base_set).collect();
        if forward.is_empty() || forward.contains(&0) {
            return Err(KernelError::InvalidOption {
                name: "pattern",
                value: pattern.to_string(),
                expected: "a non-empty string of IUPAC nucleotide codes (ACGTURYSWKMBDHVN)",
            });
        }
        let reverse = both_strands
            .then(|| forward.iter().rev().map(|&set| complement(set)).collect())
            .filter(|reverse: &Vec<u8>| *reverse != forward);
        Ok(Self {
            forward,
            reverse,
            overlapping,
        })
    }

    /// Replaces `out` with the 0-based start of every hit in `seq`, sorted
    /// and without duplicates. A hit of the reverse complement is reported
    /// at its leftmost position on `seq`, like a forward hit.
    pub fn find_into(&self, seq: &[u8], out: &mut Vec<u64>) {
        out.clear();
        scan(&self.forward, seq, self.overlapping, out);
        if let Some(reverse) = &self.reverse {
            let forward_hits = out.len();
            scan(reverse, seq, self.overlapping, out);
            if forward_hits > 0 && out.len() > forward_hits {
                out.sort_unstable();
                out.dedup();
            }
        }
    }
}

/// The complement of a set of bases.
fn complement(set: u8) -> u8 {
    let mut out = 0;
    for (base, pair) in [(A, T), (T, A), (C, G), (G, C)] {
        if set & base != 0 {
            out |= pair;
        }
    }
    out
}

fn scan(pattern: &[u8], seq: &[u8], overlapping: bool, out: &mut Vec<u64>) {
    let n = pattern.len();
    if seq.len() < n {
        return;
    }
    let mut i = 0;
    while i + n <= seq.len() {
        let hit = pattern
            .iter()
            .zip(&seq[i..])
            .all(|(&set, &base)| SEQUENCE_BITS[base as usize] & set != 0);
        if hit {
            out.push(i as u64);
        }
        i += if hit && !overlapping { n } else { 1 };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::{Alphabet, revcomp_into};

    fn find(pattern: &str, seq: &str, both_strands: bool, overlapping: bool) -> Vec<u64> {
        let motif = Motif::new(pattern, both_strands, overlapping).unwrap();
        let mut out = vec![99]; // find_into clears it
        motif.find_into(seq.as_bytes(), &mut out);
        out
    }

    #[test]
    fn nt_search_probes() {
        // Bio.SeqUtils.nt_search(seq, pattern)[1:], Biopython 1.88.
        assert_eq!(find("ACG", "ACGTACGTAAAA", false, true), [0, 4]);
        assert_eq!(find("AA", "AAAAA", false, true), [0, 1, 2, 3]);
        assert_eq!(
            find("N", "ACGTNACGT", false, true),
            [0, 1, 2, 3, 5, 6, 7, 8]
        );
        assert_eq!(find("ACGTN", "ACGTNACGT", false, true), [] as [u64; 0]);
        assert_eq!(find("R", "ACGTRYACGT", false, true), [0, 2, 6, 8]);
        assert_eq!(find("GAATTC", "GAATTCGAATTC", false, true), [0, 6]);
    }

    #[test]
    fn case_and_uracil_are_folded() {
        // Deliberate differences from nt_search, which finds nothing here.
        assert_eq!(find("ACG", "acgtacgt", false, true), [0, 4]);
        assert_eq!(find("ACGT", "ACGUACGU", false, true), [0, 4]);
        assert_eq!(find("acg", "ACGTACGT", false, true), [0, 4]);
        assert_eq!(find("ACGU", "ACGTACGT", false, true), [0, 4]);
    }

    #[test]
    fn non_overlapping() {
        assert_eq!(find("AA", "AAAAA", false, false), [0, 2]);
    }

    #[test]
    fn both_strands() {
        // A palindromic site is reported once.
        assert_eq!(find("GAATTC", "TTGAATTCTT", true, true), [2]);
        // AAC forward at 0; its reverse complement GTT at 3.
        assert_eq!(find("AAC", "AACGTT", true, true), [0, 3]);
        assert_eq!(find("AAC", "AACGTT", false, true), [0]);
        // An ambiguous pattern: RAC is AAC or GAC; reverse complement GTY.
        assert_eq!(find("RAC", "GACTTGTC", true, true), [0, 5]);
    }

    #[test]
    fn both_strands_non_overlapping_are_independent_scans() {
        // ACG at 0 and its reverse complement CGT at 1 overlap each other.
        assert_eq!(find("ACG", "ACGT", true, false), [0, 1]);
        // Within one orientation hits do not overlap: AA forward at 0, 2;
        // TT (reverse complement) nowhere.
        assert_eq!(find("AA", "AAAAA", true, false), [0, 2]);
    }

    #[test]
    fn short_sequences() {
        assert_eq!(find("ACGT", "ACG", true, true), [] as [u64; 0]);
        assert_eq!(find("A", "", true, true), [] as [u64; 0]);
    }

    #[test]
    fn invalid_patterns() {
        for pattern in ["", "AC!", "ACé", "A-C", "AXC"] {
            let err = Motif::new(pattern, false, true).unwrap_err();
            assert!(
                matches!(
                    err,
                    KernelError::InvalidOption {
                        name: "pattern",
                        ..
                    }
                ),
                "{pattern:?}: {err:?}"
            );
        }
    }

    /// The bases an IUPAC code stands for.
    fn expansion(code: u8) -> &'static [u8] {
        match code.to_ascii_uppercase() {
            b'A' => b"A",
            b'C' => b"C",
            b'G' => b"G",
            b'T' | b'U' => b"T",
            b'R' => b"AG",
            b'Y' => b"CT",
            b'S' => b"CG",
            b'W' => b"AT",
            b'K' => b"GT",
            b'M' => b"AC",
            b'B' => b"CGT",
            b'D' => b"AGT",
            b'H' => b"ACT",
            b'V' => b"ACG",
            b'N' => b"ACGT",
            _ => b"",
        }
    }

    /// A direct transcription of the matching rule, one position at a time.
    fn scan(pattern: &[u8], seq: &[u8], overlapping: bool) -> Vec<u64> {
        let mut hits = Vec::new();
        let mut i = 0;
        while i + pattern.len() <= seq.len() {
            let hit = pattern.iter().zip(&seq[i..]).all(|(&p, &b)| {
                let b = match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    b => b,
                };
                expansion(p).contains(&b)
            });
            if hit {
                hits.push(i as u64);
            }
            i += if hit && !overlapping {
                pattern.len()
            } else {
                1
            };
        }
        hits
    }

    fn reference(pattern: &str, seq: &str, both_strands: bool, overlapping: bool) -> Vec<u64> {
        let mut hits = scan(pattern.as_bytes(), seq.as_bytes(), overlapping);
        if both_strands {
            let mut rc = Vec::new();
            revcomp_into(pattern.as_bytes(), Alphabet::Dna, false, &mut rc).unwrap();
            hits.extend(scan(&rc, seq.as_bytes(), overlapping));
            hits.sort_unstable();
            hits.dedup();
        }
        hits
    }

    mod props {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn agrees_with_the_reference_scan(
                seq in "[ACGTacgtNRYU-]{0,60}",
                pattern in "[ACGTRYSWKMBDHVN]{1,5}",
                both_strands: bool,
                overlapping: bool,
            ) {
                prop_assert_eq!(
                    find(&pattern, &seq, both_strands, overlapping),
                    reference(&pattern, &seq, both_strands, overlapping)
                );
            }

            #[test]
            fn hits_are_sorted_and_unique(seq in "[ACGT]{0,80}", pattern in "[ACGTN]{1,3}") {
                let hits = find(&pattern, &seq, true, true);
                prop_assert!(hits.windows(2).all(|w| w[0] < w[1]));
            }
        }
    }
}
