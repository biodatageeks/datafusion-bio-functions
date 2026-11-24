//! Open reading frame (ORF) finding.

use super::KernelError;
use super::codon_tables::CodonTable;
use super::complement::{Alphabet, revcomp_into};
use super::translate::{Ambiguity, every_expansion_in, find_table, translate_codon};

/// The strands an [`OrfFinder`] searches.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Strands {
    #[default]
    Both,
    Forward,
    Reverse,
}

impl Strands {
    /// Parses `"both"`, `"+"` or `"-"`.
    pub fn parse(value: &str) -> Result<Self, KernelError> {
        match value.to_ascii_lowercase().as_str() {
            "both" => Ok(Strands::Both),
            "+" => Ok(Strands::Forward),
            "-" => Ok(Strands::Reverse),
            _ => Err(KernelError::InvalidOption {
                name: "strands",
                value: value.to_string(),
                expected: "'both', '+' or '-'",
            }),
        }
    }
}

/// The codons an ORF may start at.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum StartCodons {
    /// The start codons of the translation table.
    #[default]
    Table,
    /// `ATG` only.
    Atg,
}

impl StartCodons {
    /// Parses `"table"` or `"ATG"`, ignoring case.
    pub fn parse(value: &str) -> Result<Self, KernelError> {
        match value.to_ascii_lowercase().as_str() {
            "table" => Ok(StartCodons::Table),
            "atg" => Ok(StartCodons::Atg),
            _ => Err(KernelError::InvalidOption {
                name: "start_codons",
                value: value.to_string(),
                expected: "'table' or 'ATG'",
            }),
        }
    }
}

/// The codon set holding only `ATG` (codon index 35 in NCBI order).
const ATG_ONLY: u64 = 1 << 35;

/// One open reading frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Orf {
    /// 0-based start on the sequence as given, for either strand.
    pub start: u64,
    /// Exclusive end on the sequence as given. The interval includes the
    /// stop codon when the ORF has one.
    pub end: u64,
    /// Whether the ORF is on the reverse strand.
    pub reverse: bool,
    /// 1, 2, 3 or -1, -2, -3, with the meaning of the translation frame.
    pub frame: i8,
    /// Amino acids, without the stop.
    pub aa_len: u64,
    /// False for an ORF that begins at the sequence start instead of a start
    /// codon. Only possible with `allow_partial`.
    pub has_start: bool,
    /// False for an ORF that runs to the sequence end without a stop codon.
    /// Only possible with `allow_partial`.
    pub has_stop: bool,
    protein: (usize, usize),
}

/// Finds open reading frames.
///
/// In each reading frame, stop codons split the codons into segments. A
/// complete ORF runs from the first start codon of a segment to the stop that
/// ends it, so a stop yields at most one ORF, the longest. Every codon in the
/// table's stop-codon set ends an ORF, including the codons that tables 27,
/// 28 and 31 also read as an amino acid. A codon counts as a start or a stop
/// only if every expansion of its IUPAC codes is one, so an ORF runs through
/// codons containing `N`.
///
/// With `allow_partial`, ORFs bounded by a sequence end are added; every
/// complete ORF is still reported unchanged:
/// - from the first complete codon of a frame to its first stop, when no
///   start codon precedes that stop (`has_start` false);
/// - from the first start codon after the last stop of a frame to its last
///   complete codon (`has_stop` false);
/// - the whole frame, when it has neither a start nor a stop codon.
///
/// Holds its buffers, so reuse one instance across a batch.
#[derive(Debug, Clone)]
pub struct OrfFinder {
    table: &'static CodonTable,
    starts: u64,
    min_len: u64,
    strands: Strands,
    allow_partial: bool,
    include_protein: bool,
    scratch: Vec<u8>,
    orfs: Vec<Orf>,
    proteins: Vec<u8>,
}

impl OrfFinder {
    /// `table_id` is an NCBI translation table id. ORFs with fewer than
    /// `min_len` amino acids are skipped; an ORF always has at least one.
    /// With `include_protein`, [`OrfFinder::protein`] returns each ORF's
    /// translation.
    pub fn new(
        table_id: i64,
        min_len: i64,
        strands: Strands,
        start_codons: StartCodons,
        allow_partial: bool,
        include_protein: bool,
    ) -> Result<Self, KernelError> {
        let table = find_table(table_id)?;
        let min_len = u64::try_from(min_len).map_err(|_| KernelError::InvalidOption {
            name: "min_len",
            value: min_len.to_string(),
            expected: "a non-negative integer",
        })?;
        Ok(Self {
            table,
            starts: match start_codons {
                StartCodons::Table => table.starts,
                StartCodons::Atg => ATG_ONLY,
            },
            min_len: min_len.max(1),
            strands,
            allow_partial,
            include_protein,
            scratch: Vec::new(),
            orfs: Vec::new(),
            proteins: Vec::new(),
        })
    }

    /// The ORFs of `seq`, sorted by start, then end, then strand (forward
    /// first). The slice is valid until the next call.
    pub fn find(&mut self, seq: &[u8]) -> &[Orf] {
        self.orfs.clear();
        self.proteins.clear();
        if self.strands != Strands::Reverse {
            self.scan_strand(seq, false);
        }
        if self.strands != Strands::Forward {
            let mut reverse = std::mem::take(&mut self.scratch);
            reverse.clear();
            // Non-strict reverse complement cannot fail.
            let _ = revcomp_into(seq, Alphabet::Dna, false, &mut reverse);
            self.scan_strand(&reverse, true);
            self.scratch = reverse;
        }
        self.orfs
            .sort_unstable_by_key(|orf| (orf.start, orf.end, orf.reverse));
        &self.orfs
    }

    /// The translation of an ORF returned by the last [`OrfFinder::find`]:
    /// its first codon as `M` when the ORF has a start codon, an ambiguous
    /// codon as `X` unless every expansion agrees, and no stop. Empty unless
    /// the finder was created with `include_protein`.
    pub fn protein(&self, orf: &Orf) -> &[u8] {
        &self.proteins[orf.protein.0..orf.protein.1]
    }

    fn scan_strand(&mut self, strand: &[u8], reverse: bool) {
        for offset in 0..3 {
            let mut first_start = None;
            let mut first_segment = true;
            let mut pos = offset;
            while pos + 3 <= strand.len() {
                let (c1, c2, c3) = (strand[pos], strand[pos + 1], strand[pos + 2]);
                if every_expansion_in(self.table.stops, c1, c2, c3) {
                    match first_start {
                        Some(start) => self.push(strand, reverse, start, pos + 3, true, true),
                        None if self.allow_partial && first_segment => {
                            self.push(strand, reverse, offset, pos + 3, false, true)
                        }
                        None => {}
                    }
                    first_start = None;
                    first_segment = false;
                } else if first_start.is_none() && every_expansion_in(self.starts, c1, c2, c3) {
                    first_start = Some(pos);
                }
                pos += 3;
            }
            // `pos` is now the end of the frame's last complete codon.
            if self.allow_partial {
                match first_start {
                    Some(start) => self.push(strand, reverse, start, pos, true, false),
                    None if first_segment => self.push(strand, reverse, offset, pos, false, false),
                    None => {}
                }
            }
        }
    }

    /// Records the ORF at `start..end` of `strand` unless it is too short.
    fn push(
        &mut self,
        strand: &[u8],
        reverse: bool,
        start: usize,
        end: usize,
        has_start: bool,
        has_stop: bool,
    ) {
        let coding_end = if has_stop { end - 3 } else { end };
        let aa_len = (coding_end.saturating_sub(start) / 3) as u64;
        if aa_len < self.min_len {
            return;
        }
        let protein_start = self.proteins.len();
        if self.include_protein {
            for (i, codon) in strand[start..coding_end].chunks_exact(3).enumerate() {
                self.proteins.push(if i == 0 && has_start {
                    b'M'
                } else {
                    translate_codon(self.table.aa, codon[0], codon[1], codon[2], Ambiguity::X)
                });
            }
        }
        let frame = (start % 3) as i8 + 1;
        let (start, end, frame) = if reverse {
            (strand.len() - end, strand.len() - start, -frame)
        } else {
            (start, end, frame)
        };
        self.orfs.push(Orf {
            start: start as u64,
            end: end as u64,
            reverse,
            frame,
            aa_len,
            has_start,
            has_stop,
            protein: (protein_start, self.proteins.len()),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::{Alphabet, Translator, revcomp_into};

    /// (start, end, strand, frame, aa_len, has_start, has_stop, protein)
    type Row = (u64, u64, char, i8, u64, bool, bool, String);

    fn rows(finder: &mut OrfFinder, seq: &str) -> Vec<Row> {
        let orfs = finder.find(seq.as_bytes()).to_vec();
        orfs.iter()
            .map(|o| {
                (
                    o.start,
                    o.end,
                    if o.reverse { '-' } else { '+' },
                    o.frame,
                    o.aa_len,
                    o.has_start,
                    o.has_stop,
                    String::from_utf8(finder.protein(o).to_vec()).unwrap(),
                )
            })
            .collect()
    }

    fn finder(table: i64, min_len: i64, starts: StartCodons, partial: bool) -> OrfFinder {
        OrfFinder::new(table, min_len, Strands::Both, starts, partial, true).unwrap()
    }

    /// ORFs with ATG-only starts, both strands, complete only.
    fn atg(seq: &str, min_len: i64) -> Vec<Row> {
        rows(&mut finder(1, min_len, StartCodons::Atg, false), seq)
    }

    fn row(start: u64, end: u64, strand: char, frame: i8, aa: &str, s: bool, e: bool) -> Row {
        (
            start,
            end,
            strand,
            frame,
            aa.len() as u64,
            s,
            e,
            aa.to_string(),
        )
    }

    #[test]
    fn single_forward_orf() {
        assert_eq!(atg("ATGGCCTAA", 1), [row(0, 9, '+', 1, "MA", true, true)]);
    }

    #[test]
    fn one_orf_per_stop_from_the_first_start() {
        // Two ATGs before the same stop: one ORF, from the first.
        assert_eq!(
            atg("ATGAAAATGCCCTAA", 1),
            [row(0, 15, '+', 1, "MKMP", true, true)]
        );
    }

    #[test]
    fn two_orfs_in_one_frame() {
        // ATG GCC TAA | ATG AAA TGA
        assert_eq!(
            atg("ATGGCCTAAATGAAATGA", 1),
            [
                row(0, 9, '+', 1, "MA", true, true),
                row(9, 18, '+', 1, "MK", true, true)
            ]
        );
    }

    #[test]
    fn reverse_strand_coordinates_are_on_the_given_sequence() {
        // The reverse complement of ATGGCCTAA.
        assert_eq!(atg("TTAGGCCAT", 1), [row(0, 9, '-', -1, "MA", true, true)]);
        // With flanks: C + TTAGGCCAT + GG. On the reverse complement the ORF
        // starts at offset 2, so the frame is -3 and the interval is [1, 10).
        assert_eq!(
            atg("CTTAGGCCATGG", 1),
            [row(1, 10, '-', -3, "MA", true, true)]
        );
    }

    #[test]
    fn frames_two_and_three() {
        assert_eq!(atg("CATGGCCTAA", 1), [row(1, 10, '+', 2, "MA", true, true)]);
        assert_eq!(
            atg("CCATGGCCTAAC", 1),
            [row(2, 11, '+', 3, "MA", true, true)]
        );
    }

    #[test]
    fn min_len_counts_amino_acids_without_the_stop() {
        assert_eq!(atg("ATGGCCTAA", 2).len(), 1);
        assert!(atg("ATGGCCTAA", 3).is_empty());
        // min_len 0 still needs one amino acid: a lone stop is not an ORF.
        let mut forward =
            OrfFinder::new(1, 0, Strands::Forward, StartCodons::Atg, true, true).unwrap();
        assert!(rows(&mut forward, "TAA").is_empty());
    }

    #[test]
    fn strands_option() {
        let seq = "ATGGCCTAAACCTTAGGCCAT"; // forward ORF, spacer, reverse ORF
        let run = |strands| {
            let mut f = OrfFinder::new(1, 1, strands, StartCodons::Atg, false, true).unwrap();
            rows(&mut f, seq)
        };
        let both = run(Strands::Both);
        assert_eq!(
            both,
            [
                row(0, 9, '+', 1, "MA", true, true),
                row(12, 21, '-', -1, "MA", true, true)
            ]
        );
        assert_eq!(run(Strands::Forward), both[..1]);
        assert_eq!(run(Strands::Reverse), both[1..]);
    }

    #[test]
    fn table_start_codons() {
        // CTG is a start codon of table 1 and is emitted as M, as cds=True does.
        let seq = "CTGGCCTAA";
        assert_eq!(
            rows(&mut finder(1, 1, StartCodons::Table, false), seq),
            [row(0, 9, '+', 1, "MA", true, true)]
        );
        assert!(atg(seq, 1).is_empty());
    }

    #[test]
    fn ambiguous_codons() {
        // NNN is neither a start nor a stop; TAN is only possibly a stop.
        assert_eq!(atg("ATGNNNTAA", 1), [row(0, 9, '+', 1, "MX", true, true)]);
        assert_eq!(
            atg("ATGTANGCCTAA", 1),
            [row(0, 12, '+', 1, "MXA", true, true)]
        );
        // TAR is a stop either way.
        assert_eq!(atg("ATGGCCTAR", 1), [row(0, 9, '+', 1, "MA", true, true)]);
        // Lowercase and RNA.
        assert_eq!(atg("auggccuaa", 1), [row(0, 9, '+', 1, "MA", true, true)]);
    }

    #[test]
    fn dual_meaning_codons_end_orfs() {
        // Table 27: TGA is a stop or W. An ORF finder treats it as a stop.
        assert_eq!(
            rows(
                &mut finder(27, 1, StartCodons::Table, false),
                "ATGGCCTGAGCC"
            ),
            [row(0, 9, '+', 1, "MA", true, true)]
        );
    }

    fn partial(seq: &str) -> Vec<Row> {
        rows(&mut finder(1, 1, StartCodons::Atg, true), seq)
    }

    fn frame1(rows: Vec<Row>) -> Vec<Row> {
        rows.into_iter().filter(|r| r.3 == 1).collect()
    }

    #[test]
    fn partial_orfs() {
        // No start before the first stop.
        assert_eq!(
            frame1(partial("GCCGCCTAA")),
            [row(0, 9, '+', 1, "AA", false, true)]
        );
        assert!(frame1(atg("GCCGCCTAA", 1)).is_empty());
        // A start with no stop after it. The trailing partial codon is not part of it.
        assert_eq!(
            frame1(partial("ATGGCCGCCGC")),
            [row(0, 9, '+', 1, "MAA", true, false)]
        );
        assert!(frame1(atg("ATGGCCGCC", 1)).is_empty());
        // Neither a start nor a stop in the frame.
        assert_eq!(
            frame1(partial("GCCGCCGCC")),
            [row(0, 9, '+', 1, "AAA", false, false)]
        );
    }

    #[test]
    fn partial_only_adds() {
        // The first segment has a start: the complete ORF is reported, not
        // a longer one from the sequence start.
        assert_eq!(
            frame1(partial("GCCATGGCCTAA")),
            [row(3, 12, '+', 1, "MA", true, true)]
        );
        // After the last stop: start-less tails are not reported.
        assert_eq!(
            frame1(partial("ATGGCCTAAGCCGCC")),
            [row(0, 9, '+', 1, "MA", true, true)]
        );
    }

    #[test]
    fn short_and_empty_sequences() {
        for seq in ["", "A", "AT"] {
            assert!(partial(seq).is_empty(), "{seq:?}");
        }
    }

    #[test]
    fn without_proteins() {
        let mut f = OrfFinder::new(1, 1, Strands::Both, StartCodons::Atg, false, false).unwrap();
        let orfs = f.find(b"ATGGCCTAA").to_vec();
        assert_eq!(orfs.len(), 1);
        assert_eq!(f.protein(&orfs[0]), b"");
    }

    #[test]
    fn finder_is_reusable() {
        let mut f = finder(1, 1, StartCodons::Atg, false);
        assert_eq!(rows(&mut f, "ATGGCCTAA").len(), 1);
        assert!(rows(&mut f, "GGGGGG").is_empty());
        assert_eq!(rows(&mut f, "TTAGGCCAT")[0].2, '-');
    }

    #[test]
    fn invalid_options() {
        let new = |table, min_len| {
            OrfFinder::new(
                table,
                min_len,
                Strands::Both,
                StartCodons::Table,
                false,
                false,
            )
        };
        assert_eq!(new(7, 1).unwrap_err(), KernelError::UnknownCodonTable(7));
        assert!(matches!(
            new(1, -1).unwrap_err(),
            KernelError::InvalidOption {
                name: "min_len",
                ..
            }
        ));
        assert_eq!(Strands::parse("+").unwrap(), Strands::Forward);
        assert_eq!(Strands::parse("-").unwrap(), Strands::Reverse);
        assert_eq!(Strands::parse("Both").unwrap(), Strands::Both);
        assert!(Strands::parse("plus").is_err());
        assert_eq!(StartCodons::parse("atg").unwrap(), StartCodons::Atg);
        assert_eq!(StartCodons::parse("table").unwrap(), StartCodons::Table);
        assert!(StartCodons::parse("ATA").is_err());
    }

    fn revcomp(seq: &str) -> String {
        let mut out = Vec::new();
        revcomp_into(seq.as_bytes(), Alphabet::Dna, false, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// ORFs of one strand found from the definition, stop by stop, using
    /// the translation kernel to classify codons instead of the finder's
    /// own scan. `starts` lists the start codons as text.
    fn reference_strand(
        strand_seq: &str,
        starts: &[&str],
        partial: bool,
    ) -> Vec<(usize, usize, bool, bool)> {
        let mut found = Vec::new();
        for offset in 0..3 {
            let codons: Vec<&str> = (offset..strand_seq.len().saturating_sub(2))
                .step_by(3)
                .map(|i| &strand_seq[i..i + 3])
                .collect();
            let is_stop = |c: &str| {
                let mut out = Vec::new();
                Translator::new(1, 1, false)
                    .unwrap()
                    .translate_into(c.as_bytes(), &mut out)
                    .unwrap();
                out == b"*"
            };
            let is_start = |c: &str| starts.contains(&c.to_ascii_uppercase().as_str());
            let pos = |codon_index: usize| offset + 3 * codon_index;
            let stops: Vec<usize> = (0..codons.len()).filter(|&i| is_stop(codons[i])).collect();
            for (n, &stop) in stops.iter().enumerate() {
                let segment_start = if n == 0 { 0 } else { stops[n - 1] + 1 };
                match (segment_start..stop).find(|&i| is_start(codons[i])) {
                    Some(first) => found.push((pos(first), pos(stop) + 3, true, true)),
                    None if partial && n == 0 && stop > 0 => {
                        found.push((pos(0), pos(stop) + 3, false, true))
                    }
                    None => {}
                }
            }
            if partial {
                let tail_start = stops.last().map_or(0, |s| s + 1);
                match (tail_start..codons.len()).find(|&i| is_start(codons[i])) {
                    Some(first) => found.push((pos(first), pos(codons.len()), true, false)),
                    None if stops.is_empty() && !codons.is_empty() => {
                        found.push((pos(0), pos(codons.len()), false, false))
                    }
                    None => {}
                }
            }
        }
        found
    }

    fn reference(seq: &str, starts: &[&str], partial: bool) -> Vec<(u64, u64, char, bool, bool)> {
        let n = seq.len();
        let mut all: Vec<(u64, u64, char, bool, bool)> = reference_strand(seq, starts, partial)
            .into_iter()
            .map(|(s, e, hs, he)| (s as u64, e as u64, '+', hs, he))
            .collect();
        all.extend(
            reference_strand(&revcomp(seq), starts, partial)
                .into_iter()
                .map(|(s, e, hs, he)| ((n - e) as u64, (n - s) as u64, '-', hs, he)),
        );
        all.sort();
        all
    }

    mod props {
        use super::*;
        use proptest::prelude::*;

        // A small alphabet rich in starts and stops, so short sequences
        // contain ORFs on both strands.
        const DNA: &str = "(ATG|TAA|TAG|TGA|CAT|TTA|CTA|TCA|[ACGT]){0,24}[ACGT]{0,2}";

        proptest! {
            #[test]
            fn agrees_with_the_definition(seq in DNA, partial: bool, table_starts: bool) {
                let (mode, starts): (_, &[&str]) = if table_starts {
                    (StartCodons::Table, &["ATG", "CTG", "TTG"])
                } else {
                    (StartCodons::Atg, &["ATG"])
                };
                let got: Vec<_> = rows(&mut finder(1, 1, mode, partial), &seq)
                    .into_iter()
                    .map(|r| (r.0, r.1, r.2, r.5, r.6))
                    .collect();
                prop_assert_eq!(got, reference(&seq, starts, partial));
            }

            #[test]
            fn partial_output_contains_the_default_output(seq in DNA) {
                let complete = rows(&mut finder(1, 1, StartCodons::Table, false), &seq);
                let with_partial = rows(&mut finder(1, 1, StartCodons::Table, true), &seq);
                for orf in &complete {
                    prop_assert!(with_partial.contains(orf), "{:?} missing", orf);
                }
            }

            #[test]
            fn a_complete_orf_is_a_valid_cds(seq in DNA, table in prop::sample::select(vec![1i64, 2, 11, 27, 28, 31])) {
                let mut f = finder(table, 1, StartCodons::Table, false);
                for (start, end, strand, frame, aa_len, _, _, protein) in rows(&mut f, &seq) {
                    let slice = &seq[start as usize..end as usize];
                    let cds = if strand == '-' { revcomp(slice) } else { slice.to_string() };
                    let mut t = Translator::new(1, table, false).unwrap().with_cds(true);
                    prop_assert_eq!(t.cds_status(cds.as_bytes()), None, "{}", cds);
                    let mut out = Vec::new();
                    t.translate_into(cds.as_bytes(), &mut out).unwrap();
                    prop_assert_eq!(String::from_utf8(out).unwrap(), protein.clone());
                    prop_assert_eq!(aa_len as usize, protein.len());
                    prop_assert_eq!(cds.len(), 3 * (protein.len() + 1));
                    // `frame` means what it means for translate(): the full
                    // sequence read in that frame contains the protein.
                    let mut whole = Vec::new();
                    Translator::new(frame as i64, table, false)
                        .unwrap()
                        .translate_into(seq.as_bytes(), &mut whole)
                        .unwrap();
                    let whole = String::from_utf8(whole).unwrap();
                    prop_assert!(whole.contains(&protein[1..]), "{} not in {}", protein, whole);
                }
            }

            #[test]
            fn output_is_sorted(seq in DNA) {
                let got = rows(&mut finder(1, 1, StartCodons::Table, true), &seq);
                let keys: Vec<_> = got.iter().map(|r| (r.0, r.1, r.2 == '-')).collect();
                let mut sorted = keys.clone();
                sorted.sort();
                prop_assert_eq!(keys, sorted);
            }
        }
    }
}
