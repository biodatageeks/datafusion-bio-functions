//! Translation with NCBI genetic codes.

use super::KernelError;
use super::codon_tables::NCBI_TABLES;
use super::complement::{Alphabet, revcomp_into};

/// The NCBI translation table ids that [`Translator::new`] accepts.
pub fn codon_table_ids() -> impl Iterator<Item = u8> {
    NCBI_TABLES.iter().map(|(id, _)| *id)
}

/// Bit for each base, in NCBI codon order (T, C, A, G).
const T: u8 = 1;
const C: u8 = 2;
const A: u8 = 4;
const G: u8 = 8;
/// Marks a gap character.
const GAP: u8 = 16;

/// Bases each IUPAC code can stand for; 0 for an unknown byte.
const fn base_masks() -> [u8; 256] {
    let mut m = [0u8; 256];
    let codes: [(u8, u8); 16] = [
        (b'T', T),
        (b'U', T),
        (b'C', C),
        (b'A', A),
        (b'G', G),
        (b'R', A | G),
        (b'Y', C | T),
        (b'S', G | C),
        (b'W', A | T),
        (b'K', G | T),
        (b'M', A | C),
        (b'B', C | G | T),
        (b'D', A | G | T),
        (b'H', A | C | T),
        (b'V', A | C | G),
        (b'N', A | C | G | T),
    ];
    let mut i = 0;
    while i < codes.len() {
        let (code, mask) = codes[i];
        m[code as usize] = mask;
        m[code.to_ascii_lowercase() as usize] = mask;
        i += 1;
    }
    m[b'-' as usize] = GAP;
    m[b'.' as usize] = GAP;
    m
}

static MASKS: [u8; 256] = base_masks();

/// Translates nucleotide sequences in one reading frame with one genetic code.
///
/// Holds a scratch buffer for negative frames, so reuse one instance across a
/// batch.
#[derive(Debug, Clone)]
pub struct Translator {
    table: &'static [u8; 64],
    frame: i8,
    to_stop: bool,
    scratch: Vec<u8>,
}

impl Translator {
    /// `frame` is ±1, ±2 or ±3; negative frames read the reverse complement.
    /// `table_id` is an NCBI translation table id. With `to_stop`, output ends
    /// before the first stop codon.
    pub fn new(frame: i64, table_id: i64, to_stop: bool) -> Result<Self, KernelError> {
        if !matches!(frame, -3..=-1 | 1..=3) {
            return Err(KernelError::InvalidFrame(frame));
        }
        let table = NCBI_TABLES
            .iter()
            .find(|(id, _)| i64::from(*id) == table_id)
            .map(|(_, t)| *t)
            .ok_or(KernelError::UnknownCodonTable(table_id))?;
        Ok(Self {
            table,
            frame: frame as i8,
            to_stop,
            scratch: Vec::new(),
        })
    }

    /// Appends the protein translation of `seq` to `out`.
    ///
    /// - A trailing partial codon is dropped.
    /// - Stops are `*`.
    /// - An ambiguous codon translates to an amino acid only if every base it
    ///   can stand for gives the same one, otherwise `X`. (Biopython would
    ///   return `B`, `Z` or `J` for some; seqkit returns `X`.)
    /// - Three gap characters (`-` or `.`) give `-`.
    /// - Any other codon, such as one with an unknown byte or a partial gap,
    ///   gives `X`.
    pub fn translate_into(&mut self, seq: &[u8], out: &mut Vec<u8>) {
        let offset = (self.frame.unsigned_abs() - 1) as usize;
        let input: &[u8] = if self.frame < 0 {
            self.scratch.clear();
            // Non-strict reverse complement cannot fail.
            let _ = revcomp_into(seq, Alphabet::Dna, false, &mut self.scratch);
            &self.scratch
        } else {
            seq
        };
        let Some(coding) = input.get(offset..) else {
            return;
        };
        out.reserve(coding.len() / 3);
        for codon in coding.chunks_exact(3) {
            let aa = translate_codon(self.table, codon[0], codon[1], codon[2]);
            if self.to_stop && aa == b'*' {
                break;
            }
            out.push(aa);
        }
    }
}

#[inline]
fn index(b1: u8, b2: u8, b3: u8) -> usize {
    // Single-bit masks: trailing_zeros gives the NCBI base index 0..4.
    (b1.trailing_zeros() * 16 + b2.trailing_zeros() * 4 + b3.trailing_zeros()) as usize
}

fn translate_codon(table: &[u8; 64], c1: u8, c2: u8, c3: u8) -> u8 {
    let (m1, m2, m3) = (MASKS[c1 as usize], MASKS[c2 as usize], MASKS[c3 as usize]);
    if m1.count_ones() == 1 && m2.count_ones() == 1 && m3.count_ones() == 1 {
        if (m1 | m2 | m3) & GAP != 0 {
            return if m1 & m2 & m3 == GAP { b'-' } else { b'X' };
        }
        return table[index(m1, m2, m3)];
    }
    if m1 == 0 || m2 == 0 || m3 == 0 || (m1 | m2 | m3) & GAP != 0 {
        return b'X';
    }
    let mut result = None;
    for b1 in bits(m1) {
        for b2 in bits(m2) {
            for b3 in bits(m3) {
                let aa = table[index(b1, b2, b3)];
                match result {
                    None => result = Some(aa),
                    Some(r) if r != aa => return b'X',
                    _ => {}
                }
            }
        }
    }
    result.unwrap_or(b'X')
}

/// The single-bit masks set in `mask`.
fn bits(mask: u8) -> impl Iterator<Item = u8> {
    [T, C, A, G].into_iter().filter(move |b| mask & b != 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tr(seq: &str, frame: i64, table: i64, to_stop: bool) -> String {
        let mut t = Translator::new(frame, table, to_stop).unwrap();
        let mut out = Vec::new();
        t.translate_into(seq.as_bytes(), &mut out);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn oracle_probes() {
        // seqkit translate and Biopython Seq.translate agree on all of these.
        assert_eq!(tr("ATGGCCTAAGC", 1, 1, false), "MA*");
        assert_eq!(tr("GCNAARTARNNNAAN", 1, 1, false), "AK*XX");
        assert_eq!(tr("atggcc", 1, 1, false), "MA");
        assert_eq!(tr("AGATGA", 1, 2, false), "*W");
        assert_eq!(tr("ATG---GCC", 1, 1, false), "M-A");
        assert_eq!(tr("ATGGCCTAAGC", -1, 1, false), "A*A");
        assert_eq!(tr("ATGGCCTAAGC", 1, 1, true), "MA");
        assert_eq!(tr("TRA", 1, 1, false), "*");
        assert_eq!(tr("YTR", 1, 1, false), "L");
        assert_eq!(tr("TAN", 1, 1, false), "X");
        assert_eq!(tr("ATGUUU", 1, 1, false), "MF");
    }

    #[test]
    fn seqkit_not_biopython_for_ambiguity_letters() {
        // Biopython: B, Z, J. seqkit: X.
        assert_eq!(tr("RATSARMTT", 1, 1, false), "XXX");
    }

    #[test]
    fn partial_gaps_and_unknown_bytes_give_x() {
        assert_eq!(tr("ATGA-GAT-", 1, 1, false), "MXX");
        assert_eq!(tr("ATGZZZ", 1, 1, false), "MX");
        assert_eq!(tr("ATG...GCC", 1, 1, false), "M-A");
    }

    #[test]
    fn frames() {
        let seq = "ATGGCCTAAGC";
        // Biopython on the sliced (and reverse-complemented) sequence.
        assert_eq!(tr(seq, 2, 1, false), "WPK");
        assert_eq!(tr(seq, 3, 1, false), "GLS");
        assert_eq!(tr(seq, -2, 1, false), "LRP");
        assert_eq!(tr(seq, -3, 1, false), "LGH");
        // An ambiguous codon that is a stop either way ends to_stop output.
        assert_eq!(tr("ATGTARGCC", 1, 1, true), "M");
        assert_eq!(tr("AT", 3, 1, false), "");
        assert_eq!(tr("", -1, 1, false), "");
    }

    #[test]
    fn dual_meaning_codons_are_amino_acids() {
        // Tables 27, 28 and 31: TAA/TAG/TGA can be a stop or an amino acid.
        assert_eq!(tr("TAATAGTGA", 1, 27, false), "QQW");
        assert_eq!(tr("TAATAGTGA", 1, 31, true), "EEW");
    }

    #[test]
    fn invalid_options() {
        assert_eq!(
            Translator::new(0, 1, false).unwrap_err(),
            KernelError::InvalidFrame(0)
        );
        assert_eq!(
            Translator::new(4, 1, false).unwrap_err(),
            KernelError::InvalidFrame(4)
        );
        assert_eq!(
            Translator::new(1, 7, false).unwrap_err(),
            KernelError::UnknownCodonTable(7)
        );
        assert!(
            KernelError::UnknownCodonTable(7)
                .to_string()
                .contains("[1, 2, 3")
        );
    }

    mod props {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn frames_are_offsets(seq in "[ACGTN]{0,120}") {
                // Frame k is frame 1 of the sequence without its first k-1 bases.
                for k in 1..=3usize {
                    let shifted = seq.get(k - 1..).unwrap_or("");
                    prop_assert_eq!(tr(&seq, k as i64, 1, false), tr(shifted, 1, 1, false));
                }
            }

            #[test]
            fn negative_frames_read_the_reverse_complement(seq in "[ACGTRYN]{0,120}") {
                let mut rc = Vec::new();
                revcomp_into(seq.as_bytes(), Alphabet::Dna, false, &mut rc).unwrap();
                let rc = String::from_utf8(rc).unwrap();
                for k in 1..=3i64 {
                    prop_assert_eq!(tr(&seq, -k, 1, false), tr(&rc, k, 1, false));
                }
            }

            #[test]
            fn output_length(seq in "[ACGTacgt]{0,120}", frame in 1..=3i64) {
                let expected = seq.len().saturating_sub(frame as usize - 1) / 3;
                prop_assert_eq!(tr(&seq, frame, 1, false).len(), expected);
            }

            #[test]
            fn to_stop_is_a_prefix(seq in "[ACGT]{0,120}") {
                let full = tr(&seq, 1, 1, false);
                let stopped = tr(&seq, 1, 1, true);
                prop_assert!(!stopped.contains('*'));
                prop_assert_eq!(full.split('*').next().unwrap(), stopped);
            }
        }
    }
}
