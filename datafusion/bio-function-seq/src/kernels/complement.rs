//! Complement and reverse complement over the IUPAC nucleotide alphabet.

use super::KernelError;

/// Which base `A` pairs with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Alphabet {
    /// `A` pairs with `T`, and `U` with `A` (seqfu, seqtk, Biopython).
    #[default]
    Dna,
    /// `A` pairs with `U`, and `U` with `A` (Biopython `reverse_complement_rna`).
    Rna,
}

impl Alphabet {
    /// Parses `"dna"` or `"rna"`, ignoring case.
    pub fn parse(value: &str) -> Result<Self, KernelError> {
        match value.to_ascii_lowercase().as_str() {
            "dna" => Ok(Alphabet::Dna),
            "rna" => Ok(Alphabet::Rna),
            _ => Err(KernelError::InvalidOption {
                name: "alphabet",
                value: value.to_string(),
                expected: "'dna' or 'rna'",
            }),
        }
    }

    fn table(self) -> &'static [u8; 256] {
        match self {
            Alphabet::Dna => &DNA_COMPLEMENT,
            Alphabet::Rna => &RNA_COMPLEMENT,
        }
    }
}

/// Upper-case complement pairs. S, W and N complement to themselves.
const PAIRS: [(u8, u8); 16] = [
    (b'A', b'T'),
    (b'C', b'G'),
    (b'G', b'C'),
    (b'T', b'A'),
    (b'U', b'A'),
    (b'R', b'Y'),
    (b'Y', b'R'),
    (b'S', b'S'),
    (b'W', b'W'),
    (b'K', b'M'),
    (b'M', b'K'),
    (b'B', b'V'),
    (b'V', b'B'),
    (b'D', b'H'),
    (b'H', b'D'),
    (b'N', b'N'),
];

const fn complement_table(rna: bool) -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }
    let mut p = 0;
    while p < PAIRS.len() {
        let (from, mut to) = PAIRS[p];
        if rna && from == b'A' {
            to = b'U';
        }
        table[from as usize] = to;
        table[from.to_ascii_lowercase() as usize] = to.to_ascii_lowercase();
        p += 1;
    }
    table
}

const fn iupac_table() -> [bool; 256] {
    let mut table = [false; 256];
    let mut p = 0;
    while p < PAIRS.len() {
        let from = PAIRS[p].0;
        table[from as usize] = true;
        table[from.to_ascii_lowercase() as usize] = true;
        p += 1;
    }
    table
}

static DNA_COMPLEMENT: [u8; 256] = complement_table(false);
static RNA_COMPLEMENT: [u8; 256] = complement_table(true);
static IUPAC: [bool; 256] = iupac_table();

/// The character `bytes` starts with, or U+FFFD if it is not valid UTF-8.
fn first_char(bytes: &[u8]) -> char {
    bytes
        .utf8_chunks()
        .next()
        .and_then(|chunk| chunk.valid().chars().next())
        .unwrap_or(char::REPLACEMENT_CHARACTER)
}

fn check_alphabet(seq: &[u8]) -> Result<(), KernelError> {
    match seq.iter().position(|&b| !IUPAC[b as usize]) {
        Some(position) => Err(KernelError::InvalidBase {
            symbol: first_char(&seq[position..]),
            position,
        }),
        None => Ok(()),
    }
}

/// Appends the complement of `seq` to `out`.
///
/// Case is preserved. Bytes outside the IUPAC alphabet, including the gap
/// characters `-` and `.`, are copied unchanged unless `strict` is set, in
/// which case the first one is reported and nothing is written.
pub fn complement_into(
    seq: &[u8],
    alphabet: Alphabet,
    strict: bool,
    out: &mut Vec<u8>,
) -> Result<(), KernelError> {
    if strict {
        check_alphabet(seq)?;
    }
    let table = alphabet.table();
    out.extend(seq.iter().map(|&b| table[b as usize]));
    Ok(())
}

/// Appends the reverse complement of `seq` to `out`.
///
/// Same rules as [`complement_into`]. A multi-byte UTF-8 character, which can
/// only be an unknown symbol, is moved as a whole so the output stays valid
/// UTF-8 whenever the input is.
pub fn revcomp_into(
    seq: &[u8],
    alphabet: Alphabet,
    strict: bool,
    out: &mut Vec<u8>,
) -> Result<(), KernelError> {
    if strict {
        check_alphabet(seq)?;
    }
    let table = alphabet.table();
    if seq.is_ascii() {
        out.extend(seq.iter().rev().map(|&b| table[b as usize]));
        return Ok(());
    }
    out.reserve(seq.len());
    let mut end = seq.len();
    while end > 0 {
        let mut start = end - 1;
        // Step back over UTF-8 continuation bytes to the start of the character.
        while start > 0 && seq[start] & 0xC0 == 0x80 {
            start -= 1;
        }
        if end - start == 1 {
            out.push(table[seq[start] as usize]);
        } else {
            out.extend_from_slice(&seq[start..end]);
        }
        end = start;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rc(seq: &str, alphabet: Alphabet) -> String {
        let mut out = Vec::new();
        revcomp_into(seq.as_bytes(), alphabet, false, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    fn comp(seq: &str, alphabet: Alphabet) -> String {
        let mut out = Vec::new();
        complement_into(seq.as_bytes(), alphabet, false, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn iupac_table_matches_every_oracle() {
        // seqkit seq -r -p, seqfu rc, seqtk seq -r and Biopython all agree.
        assert_eq!(rc("ACGTRYSWKMBDHVN", Alphabet::Dna), "NBDHVKMWSRYACGT");
    }

    #[test]
    fn case_is_preserved() {
        assert_eq!(rc("acgtrysw", Alphabet::Dna), "wsryacgt");
        assert_eq!(rc("ACGTN", Alphabet::Dna), "NACGT");
        assert_eq!(rc("aacc", Alphabet::Dna), "ggtt");
    }

    #[test]
    fn gaps_and_unknown_bytes_pass_through() {
        assert_eq!(rc("AC-G.T", Alphabet::Dna), "A.C-GT");
        assert_eq!(rc("ACGZ", Alphabet::Dna), "ZCGT");
        assert_eq!(comp("ACGZ", Alphabet::Dna), "TGCZ");
    }

    #[test]
    fn uracil_handling() {
        assert_eq!(rc("ACGU", Alphabet::Dna), "ACGT");
        assert_eq!(rc("ACGU", Alphabet::Rna), "ACGU");
        assert_eq!(comp("acgu", Alphabet::Rna), "ugca");
    }

    #[test]
    fn multibyte_characters_stay_intact() {
        assert_eq!(rc("AéC", Alphabet::Dna), "GéT");
        assert_eq!(rc("€A", Alphabet::Dna), "T€");
    }

    #[test]
    fn strict_mode_reports_the_first_invalid_byte() {
        let mut out = Vec::new();
        let err = revcomp_into(b"ACGZ-", Alphabet::Dna, true, &mut out).unwrap_err();
        assert_eq!(
            err,
            KernelError::InvalidBase {
                symbol: 'Z',
                position: 3
            }
        );
        assert!(out.is_empty());
        assert!(err.to_string().contains("'Z'"));
        let err = complement_into(b"AC-G", Alphabet::Dna, true, &mut out).unwrap_err();
        assert_eq!(
            err,
            KernelError::InvalidBase {
                symbol: '-',
                position: 2
            }
        );
        revcomp_into(b"ACGTRYSWKMBDHVNacgtu", Alphabet::Dna, true, &mut out).unwrap();
        // A multi-byte character is reported whole, at its byte offset.
        let err = revcomp_into("ACéT".as_bytes(), Alphabet::Dna, true, &mut out).unwrap_err();
        assert_eq!(
            err,
            KernelError::InvalidBase {
                symbol: 'é',
                position: 2
            }
        );
        assert!(err.to_string().contains("'é' at byte offset 2"), "{err}");
        let err = revcomp_into(b"A\tC", Alphabet::Dna, true, &mut out).unwrap_err();
        assert!(err.to_string().contains("U+0009"), "{err}");
    }

    #[test]
    fn empty_input() {
        assert_eq!(rc("", Alphabet::Dna), "");
    }

    #[test]
    fn alphabet_parse() {
        assert_eq!(Alphabet::parse("DNA").unwrap(), Alphabet::Dna);
        assert_eq!(Alphabet::parse("rna").unwrap(), Alphabet::Rna);
        assert!(Alphabet::parse("protein").is_err());
    }

    mod props {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn revcomp_is_an_involution_without_uracil(seq in "[ACGTRYSWKMBDHVNacgtryswkmbdhvn.\\-XZ*]{0,200}") {
                let mut once = Vec::new();
                revcomp_into(seq.as_bytes(), Alphabet::Dna, false, &mut once).unwrap();
                let mut twice = Vec::new();
                revcomp_into(&once, Alphabet::Dna, false, &mut twice).unwrap();
                prop_assert_eq!(twice, seq.as_bytes());
            }

            #[test]
            fn revcomp_is_reversed_complement(seq in "\\PC{0,100}") {
                let mut rc = Vec::new();
                revcomp_into(seq.as_bytes(), Alphabet::Dna, false, &mut rc).unwrap();
                let mut c = Vec::new();
                complement_into(seq.as_bytes(), Alphabet::Dna, false, &mut c).unwrap();
                let c = String::from_utf8(c).unwrap();
                let expected: String = c.chars().rev().collect();
                prop_assert_eq!(String::from_utf8(rc).unwrap(), expected);
            }
        }
    }
}
