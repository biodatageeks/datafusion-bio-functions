//! Position frequency matrix scoring behind VEP's `HIGH_INF_POS` and
//! `MOTIF_SCORE_CHANGE`.
//!
//! Mirrors, in order:
//! - `Bio::EnsEMBL::Funcgen::BindingMatrix::Converter` — the
//!   frequencies → probabilities → {weights, bits} conversions,
//! - `Bio::EnsEMBL::Funcgen::BindingMatrix` — `is_position_informative` and
//!   `relative_sequence_similarity_score`,
//! - `Bio::EnsEMBL::Variation::MotifFeatureVariationAllele::motif_score_delta`.
//!
//! Matrices arrive from the cache as the flat `binding_matrix_elements`
//! encoding: `A,C,G,T` frequencies per position, positions joined by `;`.

/// Pseudocount added to every frequency before normalising to probabilities.
/// `Converter::from_frequencies_to_probabilities` defaults to 0.1.
const PSEUDOCOUNT: f64 = 0.1;

/// Background probability per nucleotide, the Converter's default expected
/// frequency. Weights are `log2(observed / expected)`.
const EXPECTED_FREQUENCY: f64 = 0.25;

/// Information content, in bits, above which a position counts as informative.
/// `BindingMatrix::is_position_informative` defaults to 1.5 and VEP never
/// overrides it.
const INFORMATIVE_THRESHOLD_BITS: f64 = 1.5;

/// The unit `binding_matrix_elements` must be in for scoring to be meaningful.
/// Ensembl caches only ship frequency matrices, and the Converter refuses any
/// other unit outright.
pub(crate) const FREQUENCIES_UNIT: &str = "Frequencies";

/// A binding matrix as frequencies, one `[A, C, G, T]` column per position.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MotifMatrix {
    columns: Vec<[f64; 4]>,
}

/// Index of a nucleotide in a matrix column, or `None` for anything outside
/// `ACGT` (`N` and IUPAC codes included — Perl throws on those, we decline to
/// score).
fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl MotifMatrix {
    /// Parses the flat cache encoding. Yields `None` if any column is
    /// malformed, so a half-read matrix never produces a score.
    pub(crate) fn parse(elements: &str) -> Option<Self> {
        let mut columns = Vec::new();
        for column in elements.split(';') {
            let mut frequencies = [0.0f64; 4];
            let mut count = 0usize;
            for (index, value) in column.split(',').enumerate() {
                if index >= 4 {
                    return None;
                }
                frequencies[index] = value.trim().parse::<f64>().ok()?;
                count += 1;
            }
            if count != 4 {
                return None;
            }
            columns.push(frequencies);
        }
        (!columns.is_empty()).then_some(Self { columns })
    }

    pub(crate) fn len(&self) -> usize {
        self.columns.len()
    }

    /// `Converter::from_frequencies_to_probabilities`:
    /// `(frequency + pseudocount) / (column_sum + 4 * pseudocount)`.
    fn probabilities(&self, index: usize) -> [f64; 4] {
        let column = self.columns[index];
        let sum: f64 = column.iter().sum();
        let denominator = sum + 4.0 * PSEUDOCOUNT;
        let mut probabilities = [0.0f64; 4];
        for (slot, frequency) in probabilities.iter_mut().zip(column.iter()) {
            *slot = (frequency + PSEUDOCOUNT) / denominator;
        }
        probabilities
    }

    /// `Converter::from_probabilities_to_weights`: `log2(p / 0.25)`.
    fn weights(&self, index: usize) -> [f64; 4] {
        let mut weights = self.probabilities(index);
        for weight in weights.iter_mut() {
            *weight = (*weight / EXPECTED_FREQUENCY).log2();
        }
        weights
    }

    /// `BindingMatrix::is_position_informative`, on a 1-based position.
    ///
    /// The bits matrix is `p * (2 - H)` per nucleotide with
    /// `H = -sum(p * log2(p))`; the position is informative when those bits sum
    /// above 1.5.
    pub(crate) fn is_position_informative(&self, position: usize) -> bool {
        if position < 1 || position > self.len() {
            return false;
        }
        let probabilities = self.probabilities(position - 1);
        let entropy: f64 = probabilities
            .iter()
            .map(|probability| -probability * probability.log2())
            .sum();
        let information_content = 2.0 - entropy;
        let bits: f64 = probabilities
            .iter()
            .map(|probability| probability * information_content)
            .sum();
        bits > INFORMATIVE_THRESHOLD_BITS
    }

    /// `BindingMatrix::relative_sequence_similarity_score` on the log scale
    /// (VEP never passes `$linear`).
    ///
    /// Perl throws when the sequence has the wrong length or non-`ACGT`
    /// characters; both yield `None` here, which surfaces as an empty field
    /// rather than a fatal error.
    pub(crate) fn relative_sequence_similarity_score(&self, sequence: &str) -> Option<f64> {
        let bases = sequence.as_bytes();
        if bases.len() != self.len() {
            return None;
        }

        let mut score = 0.0f64;
        let mut minimum = 0.0f64;
        let mut maximum = 0.0f64;
        for (index, base) in bases.iter().enumerate() {
            let weights = self.weights(index);
            score += weights[base_index(*base)?];
            minimum += weights.iter().copied().fold(f64::INFINITY, f64::min);
            maximum += weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        }

        let span = maximum - minimum;
        (span != 0.0).then(|| (score - minimum) / span)
    }
}

/// Where a variant sits inside a motif, in motif orientation, as
/// `MotifFeatureVariationAllele::motif_start` / `motif_end` produce it.
///
/// On the reverse strand both are mirrored, which puts `end` *before* `start`
/// for multi-base variants. Ensembl release 116 deliberately still splices at
/// `motif_start - 1`, not at the lower of the two mirrored coordinates:
/// <https://github.com/Ensembl/ensembl-variation/blob/release/116/modules/Bio/EnsEMBL/Variation/MotifFeatureVariationAllele.pm#L211-L251>.
/// Preserve that behavior even though it looks counter-intuitive; this module
/// implements VEP concordance rather than an independent PWM correction.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MotifPlacement {
    pub start: i64,
    pub end: i64,
}

/// `MotifFeatureVariationAllele::motif_score_delta`: the change in relative
/// binding affinity between the reference and variant motif sequences.
///
/// `reference_allele` and `variant_allele` must already be in motif
/// orientation (reverse-complemented for reverse-strand motifs) and in VEP's
/// allele convention, where `-` marks an empty allele.
///
/// Yields `None` — an empty `MOTIF_SCORE_CHANGE` — whenever Ensembl bails out:
/// an indel allele, a length change, a placement outside the motif, or a
/// sequence the matrix cannot score.
pub(crate) fn motif_score_delta(
    matrix: &MotifMatrix,
    motif_seq: &str,
    motif_span_length: i64,
    placement: MotifPlacement,
    reference_allele: &str,
    variant_allele: &str,
) -> Option<f64> {
    // "we can't call a score because the sequence will change length"
    if reference_allele == "-"
        || variant_allele == "-"
        || reference_allele.len() != variant_allele.len()
    {
        return None;
    }

    let sequence_length = motif_seq.len() as i64;
    let mut start = placement.start;
    let mut allele = variant_allele;

    // Ensembl trims the allele when the variant runs off either end of the
    // motif. `motif_start` already rejects anything before position 1, so only
    // the trailing overhang is reachable.
    if placement.end > sequence_length {
        let keep = sequence_length - start + 1;
        if keep <= 0 {
            return None;
        }
        allele = allele.get(..keep as usize)?;
    }
    if start < 1 {
        allele = allele.get((1 - start) as usize..)?;
        start = 1;
    }

    let allele_length = allele.len() as i64;
    if allele_length > motif_span_length {
        return None;
    }

    let reference_affinity = matrix.relative_sequence_similarity_score(motif_seq)?;

    // Splice the variant allele in (0-based). Perl would grow the string if the
    // replacement ran past the end and then discard the result on its
    // length check; refuse up front instead.
    let offset = (start - 1) as usize;
    if offset + allele.len() > motif_seq.len() {
        return None;
    }
    let mut variant_seq = String::with_capacity(motif_seq.len());
    variant_seq.push_str(motif_seq.get(..offset)?);
    variant_seq.push_str(allele);
    variant_seq.push_str(motif_seq.get(offset + allele.len()..)?);

    let variant_affinity = matrix.relative_sequence_similarity_score(&variant_seq)?;
    Some(variant_affinity - reference_affinity)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ENSPFM0085 (EBF1), the release-116 matrix for the chr22 motifs used in
    /// the parity runs.
    fn ebf1_matrix() -> MotifMatrix {
        MotifMatrix::parse(
            "980,99,249,216;234,210,463,636;37,146,3,1544;5,1544,0,1;3,1544,9,61;\
             0,1544,0,6;503,314,71,656;649,69,285,541;21,7,1544,0;63,2,1544,0;\
             7,2,1544,15;1544,13,170,24;646,456,210,232;187,281,97,979",
        )
        .unwrap()
    }

    #[test]
    fn parses_the_flat_cache_encoding() {
        let matrix = MotifMatrix::parse("980,99,249,216;234,210,463,636").unwrap();
        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix.columns[0], [980.0, 99.0, 249.0, 216.0]);
        assert_eq!(matrix.columns[1], [234.0, 210.0, 463.0, 636.0]);
    }

    #[test]
    fn refuses_malformed_matrices() {
        assert_eq!(MotifMatrix::parse(""), None);
        assert_eq!(MotifMatrix::parse("980,99,249"), None);
        assert_eq!(MotifMatrix::parse("980,99,249,216,3"), None);
        assert_eq!(MotifMatrix::parse("980,99,249,x"), None);
    }

    /// Position 4 of EBF1 is almost pure C (5,1544,0,1) and position 1 is
    /// spread across all four bases, so one is informative and the other is
    /// not.
    #[test]
    fn informative_positions_track_information_content() {
        let matrix = ebf1_matrix();
        assert!(matrix.is_position_informative(4));
        assert!(!matrix.is_position_informative(1));
    }

    #[test]
    fn positions_outside_the_matrix_are_not_informative() {
        let matrix = ebf1_matrix();
        assert!(!matrix.is_position_informative(0));
        assert!(!matrix.is_position_informative(matrix.len() + 1));
    }

    /// The consensus sequence must score at or near 1, and a sequence built
    /// from each position's rarest base at or near 0.
    #[test]
    fn relative_score_spans_zero_to_one() {
        let matrix = ebf1_matrix();
        let best = matrix.relative_sequence_similarity_score("ATCCCCTAGGGACT");
        let worst = matrix.relative_sequence_similarity_score("GGGGGGGGTTTTGT");
        let best = best.unwrap();
        let worst = worst.unwrap();
        assert!(best > worst);
        assert!((0.0..=1.0).contains(&best));
        assert!((0.0..=1.0).contains(&worst));
    }

    #[test]
    fn relative_score_declines_wrong_length_or_ambiguous_sequence() {
        let matrix = ebf1_matrix();
        assert_eq!(matrix.relative_sequence_similarity_score("ACGT"), None);
        assert_eq!(
            matrix.relative_sequence_similarity_score("NCTCTCCAGGGATT"),
            None
        );
    }

    #[test]
    fn score_delta_is_the_change_in_relative_affinity() {
        let matrix = ebf1_matrix();
        let motif_seq = "ACTCTCCAGGGATT";
        let placement = MotifPlacement { start: 4, end: 4 };
        let delta = motif_score_delta(&matrix, motif_seq, 14, placement, "C", "A").unwrap();

        let reference = matrix
            .relative_sequence_similarity_score(motif_seq)
            .unwrap();
        let variant = matrix
            .relative_sequence_similarity_score("ACTATCCAGGGATT")
            .unwrap();
        assert!((delta - (variant - reference)).abs() < 1e-12);
        // Replacing the near-invariant C at position 4 must lose affinity.
        assert!(delta < 0.0);
    }

    #[test]
    fn score_delta_declines_indels_and_length_changes() {
        let matrix = ebf1_matrix();
        let motif_seq = "ACTCTCCAGGGATT";
        let placement = MotifPlacement { start: 4, end: 4 };
        assert_eq!(
            motif_score_delta(&matrix, motif_seq, 14, placement, "C", "-"),
            None
        );
        assert_eq!(
            motif_score_delta(&matrix, motif_seq, 14, placement, "-", "A"),
            None
        );
        assert_eq!(
            motif_score_delta(&matrix, motif_seq, 14, placement, "C", "AT"),
            None
        );
    }

    /// A variant whose allele runs past the motif end is trimmed to fit, the
    /// way Ensembl trims it.
    #[test]
    fn score_delta_trims_an_allele_overhanging_the_motif() {
        let matrix = ebf1_matrix();
        let motif_seq = "ACTCTCCAGGGATT";
        let placement = MotifPlacement { start: 13, end: 15 };
        let delta = motif_score_delta(&matrix, motif_seq, 14, placement, "TTG", "AAA").unwrap();

        let reference = matrix
            .relative_sequence_similarity_score(motif_seq)
            .unwrap();
        let variant = matrix
            .relative_sequence_similarity_score("ACTCTCCAGGGAAA")
            .unwrap();
        assert!((delta - (variant - reference)).abs() < 1e-12);
    }

    /// Release 116 reverse-complements the allele but still uses
    /// `motif_start - 1` when `motif_end < motif_start` for a reverse-strand
    /// MNV. Pin that exact VEP behavior rather than "correcting" the offset to
    /// `min(start, end)`.
    #[test]
    fn score_delta_reverse_strand_mnv_uses_vep_116_motif_start() {
        let matrix = ebf1_matrix();
        let motif_seq = "ACTCTCCAGGGATT";
        let placement = MotifPlacement { start: 8, end: 7 };
        let delta = motif_score_delta(&matrix, motif_seq, 14, placement, "AC", "AA").unwrap();

        let reference = matrix
            .relative_sequence_similarity_score(motif_seq)
            .unwrap();
        let vep_variant = matrix
            .relative_sequence_similarity_score("ACTCTCCAAGGATT")
            .unwrap();
        let lower_coordinate_variant = matrix
            .relative_sequence_similarity_score("ACTCTCAAGGGATT")
            .unwrap();
        assert!((delta - (vep_variant - reference)).abs() < 1e-12);
        assert!((delta - (lower_coordinate_variant - reference)).abs() > 1e-12);
    }

    #[test]
    fn score_delta_declines_an_allele_longer_than_the_motif() {
        let matrix = ebf1_matrix();
        let placement = MotifPlacement { start: 1, end: 3 };
        assert_eq!(
            motif_score_delta(&matrix, "ACTCTCCAGGGATT", 2, placement, "ACT", "TGA"),
            None
        );
    }
}
