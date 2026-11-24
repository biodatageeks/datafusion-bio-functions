//! Translation with NCBI genetic codes.

use super::KernelError;
use super::codon_tables::{CodonTable, NCBI_TABLES};
use super::complement::{Alphabet, revcomp_into};

/// The NCBI translation table ids that [`Translator::new`] accepts.
pub fn codon_table_ids() -> impl Iterator<Item = u8> {
    NCBI_TABLES.iter().map(|t| t.id)
}

/// Looks up an NCBI translation table by id.
pub(crate) fn find_table(table_id: i64) -> Result<&'static CodonTable, KernelError> {
    NCBI_TABLES
        .iter()
        .find(|t| i64::from(t.id) == table_id)
        .ok_or(KernelError::UnknownCodonTable(table_id))
}

/// Bit for each base, in NCBI codon order (T, C, A, G).
pub(crate) const T: u8 = 1;
pub(crate) const C: u8 = 2;
pub(crate) const A: u8 = 4;
pub(crate) const G: u8 = 8;
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

/// The bases an IUPAC nucleotide code stands for, as a set of the `T`, `C`,
/// `A` and `G` bits. 0 for a gap or any byte that is not an IUPAC code.
#[inline]
pub(crate) fn base_set(byte: u8) -> u8 {
    MASKS[byte as usize] & (T | C | A | G)
}

/// What an ambiguous codon translates to when its expansions disagree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Ambiguity {
    /// `X`, as seqkit does.
    #[default]
    X,
    /// `B` when the expansions encode exactly D and N, `Z` for exactly E and
    /// Q, `J` for exactly I and L, otherwise `X`, as Biopython's
    /// `Seq.translate` does.
    Biopython,
}

impl Ambiguity {
    /// Parses `"x"` or `"biopython"`, ignoring case.
    pub fn parse(value: &str) -> Result<Self, KernelError> {
        match value.to_ascii_lowercase().as_str() {
            "x" => Ok(Ambiguity::X),
            "biopython" => Ok(Ambiguity::Biopython),
            _ => Err(KernelError::InvalidOption {
                name: "ambiguous",
                value: value.to_string(),
                expected: "'x' or 'biopython'",
            }),
        }
    }
}

/// The first check a sequence fails as a complete coding sequence (CDS), in
/// the order Biopython's `Seq.translate(cds=True)` applies them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CdsError {
    /// The first codon is not a start codon of the table.
    NoStartCodon,
    /// The length is not a multiple of three.
    LengthNotMultipleOf3,
    /// The final codon is not a stop codon of the table.
    NoStopCodon,
    /// A codon before the final one is a stop.
    InternalStopCodon,
}

impl CdsError {
    /// A stable name: `no_start_codon`, `length_not_multiple_of_3`,
    /// `no_stop_codon` or `internal_stop_codon`.
    pub fn as_str(self) -> &'static str {
        match self {
            CdsError::NoStartCodon => "no_start_codon",
            CdsError::LengthNotMultipleOf3 => "length_not_multiple_of_3",
            CdsError::NoStopCodon => "no_stop_codon",
            CdsError::InternalStopCodon => "internal_stop_codon",
        }
    }
}

/// Translates nucleotide sequences in one reading frame with one genetic code.
///
/// Holds a scratch buffer for negative frames, so reuse one instance across a
/// batch.
#[derive(Debug, Clone)]
pub struct Translator {
    table: &'static CodonTable,
    frame: i8,
    to_stop: bool,
    stop_symbol: Vec<u8>,
    cds: bool,
    ambiguity: Ambiguity,
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
        Ok(Self {
            table: find_table(table_id)?,
            frame: frame as i8,
            to_stop,
            stop_symbol: b"*".to_vec(),
            cds: false,
            ambiguity: Ambiguity::X,
            scratch: Vec::new(),
        })
    }

    /// The text emitted for a stop codon (default `*`). It may be empty or
    /// longer than one character, as in Biopython.
    pub fn with_stop_symbol(mut self, stop_symbol: &str) -> Self {
        self.stop_symbol = stop_symbol.as_bytes().to_vec();
        self
    }

    /// Validates every sequence as a complete CDS, as Biopython's
    /// `Seq.translate(cds=True)` does. See [`Translator::translate_into`].
    pub fn with_cds(mut self, cds: bool) -> Self {
        self.cds = cds;
        self
    }

    /// How an ambiguous codon is translated when its expansions disagree.
    pub fn with_ambiguity(mut self, ambiguity: Ambiguity) -> Self {
        self.ambiguity = ambiguity;
        self
    }

    /// Appends the protein translation of `seq` to `out`.
    ///
    /// - A trailing partial codon is dropped.
    /// - Stops are the stop symbol (`*` by default).
    /// - An ambiguous codon translates to an amino acid only if every base it
    ///   can stand for gives the same one, otherwise `X` (or `B`, `Z`, `J`
    ///   with [`Ambiguity::Biopython`]).
    /// - Three gap characters (`-` or `.`) give `-`.
    /// - Any other codon, such as one with an unknown byte or a partial gap,
    ///   gives `X`.
    ///
    /// With `cds`, the sequence in this frame must be a complete CDS: it
    /// starts with a start codon of the table, its length is a multiple of
    /// three, it ends with a stop codon, and no other codon is a stop. The
    /// first codon is then emitted as `M` and the final stop is not emitted.
    /// Otherwise this returns [`KernelError::InvalidCds`] and writes nothing.
    /// An ambiguous first or final codon passes only if every expansion is a
    /// start or stop codon. Without `cds` this never fails.
    pub fn translate_into(&mut self, seq: &[u8], out: &mut Vec<u8>) -> Result<(), KernelError> {
        let table = self.table;
        let coding = frame_input(self.frame, seq, &mut self.scratch);
        if self.cds {
            if let Some(error) = cds_ends_error(table, coding) {
                return Err(invalid_cds(error, table, coding, &[]));
            }
            let start_len = out.len();
            out.reserve(coding.len() / 3);
            out.push(b'M');
            for codon in interior(coding).chunks_exact(3) {
                let aa = translate_codon(table.aa, codon[0], codon[1], codon[2], self.ambiguity);
                if aa == b'*' {
                    out.truncate(start_len);
                    let error = CdsError::InternalStopCodon;
                    return Err(invalid_cds(error, table, coding, codon));
                }
                out.push(aa);
            }
            return Ok(());
        }
        out.reserve(coding.len() / 3);
        for codon in coding.chunks_exact(3) {
            let aa = translate_codon(table.aa, codon[0], codon[1], codon[2], self.ambiguity);
            if aa == b'*' {
                if self.to_stop {
                    break;
                }
                out.extend_from_slice(&self.stop_symbol);
            } else {
                out.push(aa);
            }
        }
        Ok(())
    }

    /// The first CDS check `seq` fails in this frame and table, or `None`
    /// when [`Translator::translate_into`] with `cds` would succeed.
    pub fn cds_status(&mut self, seq: &[u8]) -> Option<CdsError> {
        let table = self.table;
        let coding = frame_input(self.frame, seq, &mut self.scratch);
        if let Some(error) = cds_ends_error(table, coding) {
            return Some(error);
        }
        interior(coding)
            .chunks_exact(3)
            .any(|c| translate_codon(table.aa, c[0], c[1], c[2], Ambiguity::X) == b'*')
            .then_some(CdsError::InternalStopCodon)
    }
}

/// The part of `seq` a frame reads: the reverse complement for a negative
/// frame, without the first one or two bases for frames ±2 and ±3.
fn frame_input<'a>(frame: i8, seq: &'a [u8], scratch: &'a mut Vec<u8>) -> &'a [u8] {
    let offset = (frame.unsigned_abs() - 1) as usize;
    let input: &[u8] = if frame < 0 {
        scratch.clear();
        // Non-strict reverse complement cannot fail.
        let _ = revcomp_into(seq, Alphabet::Dna, false, scratch);
        scratch
    } else {
        seq
    };
    input.get(offset..).unwrap_or(&[])
}

/// The codons of a CDS between its first and its final codon.
fn interior(coding: &[u8]) -> &[u8] {
    coding.get(3..coding.len().saturating_sub(3)).unwrap_or(&[])
}

/// The first of the start-codon, length and final-stop checks that fails.
fn cds_ends_error(table: &CodonTable, coding: &[u8]) -> Option<CdsError> {
    match coding.get(..3) {
        Some(c) if every_expansion_in(table.starts, c[0], c[1], c[2]) => {}
        _ => return Some(CdsError::NoStartCodon),
    }
    if !coding.len().is_multiple_of(3) {
        return Some(CdsError::LengthNotMultipleOf3);
    }
    let last = &coding[coding.len() - 3..];
    if !every_expansion_in(table.stops, last[0], last[1], last[2]) {
        return Some(CdsError::NoStopCodon);
    }
    None
}

fn invalid_cds(error: CdsError, table: &CodonTable, coding: &[u8], codon: &[u8]) -> KernelError {
    let text = |codon: &[u8]| String::from_utf8_lossy(codon).into_owned();
    let detail = match error {
        CdsError::NoStartCodon => format!(
            "first codon '{}' is not a start codon in table {}",
            text(&coding[..coding.len().min(3)]),
            table.id
        ),
        CdsError::LengthNotMultipleOf3 => format!(
            "sequence length {} is not a multiple of three",
            coding.len()
        ),
        CdsError::NoStopCodon => format!(
            "final codon '{}' is not a stop codon in table {}",
            text(&coding[coding.len().saturating_sub(3)..]),
            table.id
        ),
        CdsError::InternalStopCodon => {
            format!("extra in-frame stop codon '{}' found", text(codon))
        }
    };
    KernelError::InvalidCds { error, detail }
}

#[inline]
fn index(b1: u8, b2: u8, b3: u8) -> usize {
    // Single-bit masks: trailing_zeros gives the NCBI base index 0..4.
    (b1.trailing_zeros() * 16 + b2.trailing_zeros() * 4 + b3.trailing_zeros()) as usize
}

/// Bit of an amino-acid letter, or of `*`, in a set of translations.
#[inline]
fn aa_bit(aa: u8) -> u32 {
    if aa == b'*' {
        1 << 26
    } else {
        1 << (aa - b'A')
    }
}

const fn aa_pair(a: u8, b: u8) -> u32 {
    (1 << (a - b'A')) | (1 << (b - b'A'))
}

const ASP_ASN: u32 = aa_pair(b'D', b'N');
const GLU_GLN: u32 = aa_pair(b'E', b'Q');
const ILE_LEU: u32 = aa_pair(b'I', b'L');

pub(crate) fn translate_codon(
    table: &[u8; 64],
    c1: u8,
    c2: u8,
    c3: u8,
    ambiguity: Ambiguity,
) -> u8 {
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
    let mut first = b'X';
    let mut seen = 0u32;
    for b1 in bits(m1) {
        for b2 in bits(m2) {
            for b3 in bits(m3) {
                let aa = table[index(b1, b2, b3)];
                if seen == 0 {
                    first = aa;
                }
                seen |= aa_bit(aa);
                if ambiguity == Ambiguity::X && seen.count_ones() > 1 {
                    return b'X';
                }
            }
        }
    }
    match (seen.count_ones(), seen) {
        (1, _) => first,
        (_, ASP_ASN) => b'B',
        (_, GLU_GLN) => b'Z',
        (_, ILE_LEU) => b'J',
        _ => b'X',
    }
}

/// Whether every codon an IUPAC codon can stand for is in `set`, a bit set
/// over the 64 codon indices such as a table's start or stop codons. False
/// for a codon with a gap or a byte that is not an IUPAC code.
pub(crate) fn every_expansion_in(set: u64, c1: u8, c2: u8, c3: u8) -> bool {
    let (m1, m2, m3) = (MASKS[c1 as usize], MASKS[c2 as usize], MASKS[c3 as usize]);
    if m1 == 0 || m2 == 0 || m3 == 0 || (m1 | m2 | m3) & GAP != 0 {
        return false;
    }
    if m1.count_ones() == 1 && m2.count_ones() == 1 && m3.count_ones() == 1 {
        return set >> index(m1, m2, m3) & 1 == 1;
    }
    bits(m1).all(|b1| bits(m2).all(|b2| bits(m3).all(|b3| set >> index(b1, b2, b3) & 1 == 1)))
}

/// The single-bit masks set in `mask`.
fn bits(mask: u8) -> impl Iterator<Item = u8> {
    [T, C, A, G].into_iter().filter(move |b| mask & b != 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tr(seq: &str, frame: i64, table: i64, to_stop: bool) -> String {
        tr_opts(seq, Translator::new(frame, table, to_stop).unwrap()).unwrap()
    }

    fn tr_opts(seq: &str, mut t: Translator) -> Result<String, KernelError> {
        let mut out = Vec::new();
        t.translate_into(seq.as_bytes(), &mut out)?;
        Ok(String::from_utf8(out).unwrap())
    }

    fn plain(table: i64) -> Translator {
        Translator::new(1, table, false).unwrap()
    }

    fn cds(table: i64) -> Translator {
        plain(table).with_cds(true)
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

    #[test]
    fn stop_symbol() {
        // Biopython Seq.translate(stop_symbol=...)
        let t = |s| plain(1).with_stop_symbol(s);
        assert_eq!(tr_opts("ATGTAAGCCTGA", t("#")).unwrap(), "M#A#");
        assert_eq!(tr_opts("ATGTAAGCCTGA", t("")).unwrap(), "MA");
        assert_eq!(tr_opts("ATGTAAGCCTGA", t("**")).unwrap(), "M**A**");
        assert_eq!(tr_opts("ATGTARGCC", t("#")).unwrap(), "M#A");
        // to_stop still ends at the stop, whatever its symbol.
        let stopped = Translator::new(1, 1, true).unwrap().with_stop_symbol("#");
        assert_eq!(tr_opts("ATGTAAGCC", stopped).unwrap(), "M");
    }

    #[test]
    fn biopython_ambiguity_letters() {
        let bio = |table| plain(table).with_ambiguity(Ambiguity::Biopython);
        assert_eq!(
            tr_opts("RATSARMTTYTAATHTANNNNTRA", bio(1)).unwrap(),
            "BZJLIXX*"
        );
        assert_eq!(tr_opts("RATSARMTT", bio(2)).unwrap(), "BZJ");
        assert_eq!(tr_opts("ratsarmtt", bio(1)).unwrap(), "BZJ");
        assert_eq!(tr_opts("RAYSAGMTY", bio(1)).unwrap(), "BZJ");
        assert_eq!(tr_opts("WTTHTT", bio(1)).unwrap(), "XX");
        // The default mode is unchanged.
        assert_eq!(tr_opts("RATSARMTTWTT", plain(1)).unwrap(), "XXXX");
        assert_eq!(Ambiguity::parse("Biopython").unwrap(), Ambiguity::Biopython);
        assert_eq!(Ambiguity::parse("x").unwrap(), Ambiguity::X);
        assert!(Ambiguity::parse("bzj").is_err());
    }

    #[test]
    fn cds_valid() {
        // Biopython Seq.translate(cds=True) with the given table.
        for (seq, table, expected) in [
            ("ATGGCCTAA", 1, "MA"),
            ("GTGGCCTAA", 11, "MA"),
            ("TTGGCCTAA", 11, "MA"),
            ("CTGGCCTAA", 1, "MA"),
            ("ATTGCCTAA", 2, "MA"),
            ("ATGTAA", 1, "M"),
            ("atggcctaa", 1, "MA"),
            ("AUGGCCUAA", 1, "MA"),
            // Ambiguous first and final codons: every expansion qualifies.
            ("ATGGCCTAR", 1, "MA"),
            ("RTGGCCTAA", 11, "MA"),
            ("NTGGCCTAA", 11, "MA"),
            // Interior codons translate as usual.
            ("ATGNNNTAA", 1, "MX"),
            ("ATGTANTAA", 1, "MX"),
            ("ATG---TAA", 1, "M-"),
            // Dual-meaning codons: a stop when final, an amino acid inside.
            ("ATGGCCTGA", 27, "MA"),
            ("ATGTGAGCCTGA", 27, "MWA"),
        ] {
            assert_eq!(
                tr_opts(seq, cds(table)).unwrap(),
                expected,
                "{seq} table {table}"
            );
            assert_eq!(
                cds(table).cds_status(seq.as_bytes()),
                None,
                "{seq} table {table}"
            );
        }
        // B, Z and J apply inside a CDS too.
        let bio = cds(1).with_ambiguity(Ambiguity::Biopython);
        assert_eq!(tr_opts("ATGRATTAA", bio).unwrap(), "MB");
        // to_stop and stop_symbol change nothing for a valid CDS.
        let t = Translator::new(1, 1, true)
            .unwrap()
            .with_cds(true)
            .with_stop_symbol("#");
        assert_eq!(tr_opts("ATGGCCTAA", t).unwrap(), "MA");
    }

    #[test]
    fn cds_invalid() {
        // The first failing check, in Biopython's order.
        for (seq, error, text) in [
            (
                "GTGGCCTAA",
                CdsError::NoStartCodon,
                "first codon 'GTG' is not a start codon",
            ),
            (
                "",
                CdsError::NoStartCodon,
                "first codon '' is not a start codon",
            ),
            (
                "AT",
                CdsError::NoStartCodon,
                "first codon 'AT' is not a start codon",
            ),
            (
                "ATGGCCTA",
                CdsError::LengthNotMultipleOf3,
                "length 8 is not a multiple of three",
            ),
            (
                "ATGGCCGCC",
                CdsError::NoStopCodon,
                "final codon 'GCC' is not a stop codon",
            ),
            (
                "ATG",
                CdsError::NoStopCodon,
                "final codon 'ATG' is not a stop codon",
            ),
            ("ATGGCCTAN", CdsError::NoStopCodon, "final codon 'TAN'"),
            (
                "ATGTAAGCCTAA",
                CdsError::InternalStopCodon,
                "in-frame stop codon 'TAA'",
            ),
            (
                "ATGGCCTAATAA",
                CdsError::InternalStopCodon,
                "in-frame stop codon 'TAA'",
            ),
            (
                "ATGTARGCCTAA",
                CdsError::InternalStopCodon,
                "in-frame stop codon 'TAR'",
            ),
            // The start codon check comes before the length check.
            ("GTGGCCTA", CdsError::NoStartCodon, "first codon 'GTG'"),
        ] {
            assert_eq!(cds(1).cds_status(seq.as_bytes()), Some(error), "{seq}");
            let err = tr_opts(seq, cds(1)).unwrap_err();
            let message = err.to_string();
            assert!(
                matches!(&err, KernelError::InvalidCds { error: e, .. } if *e == error),
                "{seq}: {err:?}"
            );
            assert!(message.contains(text), "{seq}: {message}");
            assert!(message.contains("cds_status"), "{seq}: {message}");
        }
    }

    #[test]
    fn cds_error_writes_nothing() {
        let mut out = b"keep".to_vec();
        assert!(cds(1).translate_into(b"ATGGCCTAATAA", &mut out).is_err());
        assert_eq!(out, b"keep");
    }

    #[test]
    fn cds_applies_to_the_requested_frame() {
        // The reverse complement of ATGGCCTAA, read in frame -1.
        let t = Translator::new(-1, 1, false).unwrap().with_cds(true);
        assert_eq!(tr_opts("TTAGGCCAT", t.clone()).unwrap(), "MA");
        assert_eq!(t.clone().cds_status(b"TTAGGCCAT"), None);
        // Frame 2 skips one base first.
        let t = Translator::new(2, 1, false).unwrap().with_cds(true);
        assert_eq!(tr_opts("CATGGCCTAA", t).unwrap(), "MA");
    }

    #[test]
    fn cds_error_names() {
        assert_eq!(CdsError::NoStartCodon.as_str(), "no_start_codon");
        assert_eq!(
            CdsError::LengthNotMultipleOf3.as_str(),
            "length_not_multiple_of_3"
        );
        assert_eq!(CdsError::NoStopCodon.as_str(), "no_stop_codon");
        assert_eq!(CdsError::InternalStopCodon.as_str(), "internal_stop_codon");
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
            fn cds_status_agrees_with_cds_translation(
                seq in "(ATG|GTG|TAA|TGA|GCC|[ACGTN]){0,12}",
                table in prop::sample::select(vec![1i64, 2, 11, 27, 28]),
            ) {
                let mut t = Translator::new(1, table, false).unwrap().with_cds(true);
                let mut out = Vec::new();
                let translated = t.translate_into(seq.as_bytes(), &mut out);
                prop_assert_eq!(t.cds_status(seq.as_bytes()).is_none(), translated.is_ok());
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
