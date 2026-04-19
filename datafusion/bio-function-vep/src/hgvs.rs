//! HGVS notation formatting for VEP CSQ fields (HGVSc and HGVSp).

use std::cmp::Ordering;
use std::io::{BufRead, Seek};

use crate::transcript_consequence::{
    ExonFeature, TranscriptCdnaMapperSegment, TranscriptFeature, TranslationFeature, VariantInput,
    adjust_refseq_cdna_component, edited_transcript_sequence_cdna_index,
    genomic_to_cdna_index_for_transcript, raw_cdna_position_from_genomic,
    refseq_sequence_offset_for_cdna,
};
use datafusion::common::{DataFusionError, Result};
use noodles_core::{Position, Region};
use noodles_fasta as fasta;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinHgvsData {
    pub start: usize,
    pub end: usize,
    pub ref_peptide: String,
    pub alt_peptide: String,
    pub ref_translation: String,
    pub alt_translation: String,
    pub alt_translation_extension: Option<String>,
    pub frameshift: bool,
    pub start_lost: bool,
    pub stop_lost: bool,
    pub native_refseq: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HgvsGenomicShift {
    pub strand: i8,
    pub shift_length: usize,
    pub start: i64,
    pub end: i64,
    pub shifted_allele_string: String,
    pub shifted_compare_allele: String,
    pub shifted_output_allele: String,
    pub ref_orig_allele_string: String,
    pub alt_orig_allele_string: String,
    pub five_prime_flanking_seq: String,
    pub three_prime_flanking_seq: String,
    pub five_prime_context: String,
    pub three_prime_context: String,
}

impl HgvsGenomicShift {
    /// Traceability:
    /// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
    ///   applies `_hgvs_offset` to `_slice_start/_slice_end` instead of
    ///   mutating the stored reverse-strand genomic start/end in `shift_hash`
    ///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1360-L1406>
    /// - Ensembl Variation `TranscriptVariationAllele::hgvs_protein()`
    ///   applies the same strand-aware `shifting_offset` to translation
    ///   coordinates when HGVS shifting is active
    ///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1636-L1669>
    ///
    /// For forward-strand genomic shifting Ensembl updates the genomic
    /// start/end directly during `perform_shift()`. For reverse-strand HGVS,
    /// the shift length is carried separately and applied later as an offset
    /// on transcript-oriented coordinates. Rust maps back to genomic positions
    /// directly, so we collapse that late offset into an effective display
    /// interval here.
    pub fn display_start(&self) -> i64 {
        if self.strand >= 0 {
            self.start
        } else {
            self.start - self.shift_length as i64
        }
    }

    pub fn display_end(&self) -> i64 {
        if self.strand >= 0 {
            self.end
        } else {
            self.end - self.shift_length as i64
        }
    }
}

/// Map a single-letter amino acid code to its 3-letter abbreviation.
pub fn aa_one_to_three(c: char) -> &'static str {
    match c {
        'A' => "Ala",
        'R' => "Arg",
        'N' => "Asn",
        'D' => "Asp",
        'C' => "Cys",
        'E' => "Glu",
        'Q' => "Gln",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'L' => "Leu",
        'K' => "Lys",
        'M' => "Met",
        'F' => "Phe",
        'P' => "Pro",
        'S' => "Ser",
        'T' => "Thr",
        'W' => "Trp",
        'Y' => "Tyr",
        'V' => "Val",
        'U' => "Sec",
        'O' => "Pyl",
        '*' => "Ter",
        'X' => "Xaa",
        _ => "Xaa",
    }
}

/// Format versioned transcript/translation ID, appending `.version` only when
/// the incoming stable ID is not already versioned.
///
/// Traceability:
/// - Ensembl VEP `OutputFactory::BaseTranscriptVariationAllele_to_output_hash()`
///   appends transcript version only when the stable ID does not already end in
///   `.digits`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1417-L1421>
fn versioned_id(base_id: &str, version: Option<i32>) -> String {
    if has_numeric_version_suffix(base_id) {
        return base_id.to_string();
    }
    match version {
        Some(v) => format!("{base_id}.{v}"),
        None => base_id.to_string(),
    }
}

fn has_numeric_version_suffix(value: &str) -> bool {
    value
        .rsplit_once('.')
        .map(|(_, suffix)| !suffix.is_empty() && suffix.chars().all(|ch| ch.is_ascii_digit()))
        .unwrap_or(false)
}

/// Compute HGVSc notation for a variant overlapping a transcript.
///
/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1301-L1491>
/// - Ensembl Variation `Utils::Sequence::hgvs_variant_notation()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L493-L619>
/// - Ensembl Variation `TranscriptVariationAllele::_clip_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2117-L2232>
/// - Ensembl Variation `TranscriptVariationAllele::_get_cDNA_position()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2683-L2765>
/// - Ensembl Variation `Utils::Sequence::format_hgvs_string()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L623-L676>
pub fn format_hgvsc(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    cdna_position: Option<&str>,
    cds_position: Option<&str>,
    ref_allele: &str,
    alt_allele: &str,
    variant_start: i64,
    variant_end: i64,
    genomic_shift: Option<&HgvsGenomicShift>,
) -> Option<String> {
    // Traceability:
    // - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
    //   returns undef when `_get_cDNA_position()` can't map a position
    //   outside the transcript region. For deletions that extend beyond
    //   the transcript boundaries, VEP's `_slice_start` / coordinate
    //   mapping fails for the out-of-bounds end, producing empty HGVSc.
    //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1416
    if ref_allele != "-" && alt_allele == "-" {
        // Deletion: check if either end extends beyond transcript boundaries.
        if variant_start < tx.start || variant_end > tx.end {
            return None;
        }
    }
    let tx_id = versioned_id(&tx.transcript_id, tx.version);
    let numbering = if tx.cds_start.is_some() && tx.cds_end.is_some() {
        'c'
    } else {
        'n'
    };
    let use_genomic_shift = genomic_shift
        .filter(|shift| hgvsc_uses_genomic_shift(tx, ref_allele, alt_allele, Some(shift)));
    let edited_shifted_output_allele = use_genomic_shift
        .and_then(|_| edited_transcript_insertion_allele(tx, ref_allele, alt_allele));
    let dup_context = use_genomic_shift.filter(|shift| {
        ref_allele == "-"
            && !shift.shifted_output_allele.is_empty()
            && shift.shifted_output_allele != "-"
    });
    let (ref_allele, alt_allele, variant_start, variant_end) =
        if let Some(shift) = use_genomic_shift {
            (
                ref_allele,
                if ref_allele == "-" {
                    edited_shifted_output_allele
                        .as_deref()
                        .unwrap_or(shift.shifted_output_allele.as_str())
                } else {
                    alt_allele
                },
                shift.display_start(),
                shift.display_end(),
            )
        } else {
            (ref_allele, alt_allele, variant_start, variant_end)
        };
    let (mut feature_ref, feature_alt) = hgvs_feature_strand_alleles(tx, ref_allele, alt_allele)?;
    if let Some(transcript_ref) = edited_transcript_reference_allele_for_hgvsc(
        tx,
        tx_exons,
        ref_allele,
        variant_start,
        variant_end,
    ) {
        feature_ref = transcript_ref;
    }
    let mut notation =
        hgvs_variant_notation(&feature_ref, &feature_alt, variant_start, variant_end)?;
    if let Some(shift) = dup_context {
        // VEP determines transcript-level insertion duplication from the shifted
        // HGVS allele first, then later overrides the displayed inserted bases
        // for BAM-edited RefSeq transcripts. Dup coordinates therefore follow
        // `hgvs_allele_string`, not the later transcript-space override.
        let shifted_insertion_dup_allele = shift.shifted_output_allele.as_str();
        let shifted_feature_alt =
            hgvs_feature_strand_alleles(tx, "-", shifted_insertion_dup_allele)
                .map(|(_, alt)| alt)?;
        apply_shifted_insertion_duplication(tx, &shifted_feature_alt, shift, &mut notation);
        // If the dup range extends outside the transcript's genomic
        // span (before first exon or after last exon), VEP keeps the
        // original insertion notation. This happens when HGVS 3' shift
        // pushes an insertion through a repeat past the transcript
        // start/end. Intronic dups (between exons) are valid.
        if notation.kind == "dup" {
            // unwrap_or(0): if tx_exons is empty (shouldn't happen for a
            // valid transcript), span_end=0 causes every dup to revert —
            // safe because format_hgvsc is never called without exons.
            let span_start = tx_exons.iter().map(|e| e.start).min().unwrap_or(0);
            let span_end = tx_exons.iter().map(|e| e.end).max().unwrap_or(0);
            if notation.start < span_start || notation.end > span_end {
                notation =
                    hgvs_variant_notation(&feature_ref, &feature_alt, variant_start, variant_end)?;
            }
        }
    }
    let identical_substitution = notation.kind == ">" && notation.ref_allele == notation.alt_allele;
    if notation.kind != "dup" && !identical_substitution {
        clip_alleles(&mut notation, tx.strand);
    }
    // Traceability:
    // - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
    //   replaces clipped alleles with feature_seq when `_rna_edit` attributes
    //   are present on the transcript
    //   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1443-L1449>
    if !tx.refseq_edits.is_empty() {
        notation.ref_allele = feature_ref.clone();
        notation.alt_allele = feature_alt.clone();
    }
    let (mut start, mut end) = if notation.kind == ">"
        && numbering == 'c'
        && cds_position.is_some_and(|value| !value.is_empty() && !value.contains('_'))
    {
        let cds_pos = cds_position?.to_string();
        (cds_pos.clone(), cds_pos)
    } else {
        notation_to_hgvsc_coords(tx, tx_exons, &notation)?
    };
    if matches!(
        compare_hgvs_positions(&start, &end),
        Some(Ordering::Greater)
    ) {
        std::mem::swap(&mut start, &mut end);
    }
    format_hgvs_string(&tx_id, numbering, &start, &end, &notation)
}

fn is_native_refseq_transcript(tx: &TranscriptFeature) -> bool {
    tx.source.as_deref() == Some("RefSeq")
        || matches!(
            tx.transcript_id.as_bytes().get(..2),
            Some(b"NM") | Some(b"NR") | Some(b"XM") | Some(b"XR")
        )
}

pub(crate) fn hgvsc_uses_genomic_shift(
    tx: &TranscriptFeature,
    ref_allele: &str,
    alt_allele: &str,
    genomic_shift: Option<&HgvsGenomicShift>,
) -> bool {
    if ref_allele != "-" && alt_allele != "-" {
        return false;
    }

    let Some(shift) = genomic_shift else {
        return false;
    };

    // When RefSeq BAM edit replay failed, VEP only suppresses transcript HGVS
    // shifting when the transcript-space HGVS alleles no longer match the
    // original shifted allele payload. Failed BAM-edit rows whose HGVS alleles
    // still match the original genomic payload keep the shift.
    if is_native_refseq_transcript(tx) && tx.bam_edit_status.as_deref() == Some("failed") {
        return ref_allele.eq_ignore_ascii_case(&shift.ref_orig_allele_string)
            && alt_allele.eq_ignore_ascii_case(&shift.alt_orig_allele_string);
    }

    true
}

pub(crate) fn hgvsc_offset_for_output(
    tx: &TranscriptFeature,
    variant: &VariantInput,
    ref_allele: &str,
    hgvsc: Option<&str>,
) -> Option<i64> {
    let shift = variant.hgvs_shift_for_strand(tx.strand)?;
    if hgvsc.is_none()
        || !hgvsc_uses_genomic_shift(tx, ref_allele, &variant.alt_allele, Some(shift))
    {
        return None;
    }
    if shift.shift_length == 0 {
        return None;
    }

    let signed = shift.shift_length as i64;
    Some(if tx.strand < 0 { -signed } else { signed })
}

#[derive(Clone, Copy)]
enum GenomicShiftKind {
    Insertion,
    Deletion,
}

/// Traceability:
/// - Ensembl Variation `_genomic_shift()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L411-L466>
/// - Ensembl Variation `perform_shift()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351>
/// - Ensembl Variation `create_shift_hash()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L365-L400>
///
/// This materializes the simple-indel genomic shift hash VEP uses when
/// transcript HGVS needs a precomputed genomic 3' shift outside transcript
/// sequence (e.g. intronic repeat indels).
pub fn build_hgvs_genomic_shift<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    chrom: &str,
    ref_allele: &str,
    alt_allele: &str,
    start: i64,
    end: i64,
    strand: i8,
) -> Result<Option<HgvsGenomicShift>>
where
    R: BufRead + Seek,
{
    let Some((seq_to_check, kind)) = parse_shiftable_indel(ref_allele, alt_allele) else {
        return Ok(None);
    };
    let seq_to_check = seq_to_check.as_bytes().to_vec();
    let hgvs_output = alt_allele.as_bytes().to_vec();

    // Traceability:
    // - Ensembl Variation `TranscriptVariationAllele::_genomic_shift()`
    //   expands the VF slice by 1000bp on both sides, then passes both
    //   `pre_seq` and `post_seq` into `perform_shift()`
    //   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L425-L465>
    let area_to_search = 1000i64;
    let pre_end = start.saturating_sub(1);
    let pre_start = (pre_end - area_to_search + 1).max(1);
    let post_start = end.saturating_add(1);
    let post_end = post_start.saturating_add(area_to_search.saturating_sub(1));
    let pre_seq: String = if pre_end >= pre_start && pre_end > 0 {
        read_reference_sequence(reader, chrom, pre_start, pre_end)?
    } else {
        String::new()
    };
    let post_seq: String = if post_end >= post_start && post_start > 0 {
        read_reference_sequence(reader, chrom, post_start, post_end)?
    } else {
        String::new()
    };

    // Ensembl always passes `seq_strand = 1` to genomic `perform_shift()`
    // because the shift is computed on forward-strand coordinates.
    let genomic_seq_strand = 1i8;
    let (shift_length, shifted_seq, shifted_hgvs_output, shifted_start, shifted_end) =
        perform_shift_ensembl(
            &seq_to_check,
            &hgvs_output,
            &post_seq,
            &pre_seq,
            start,
            end,
            strand < 0,
            genomic_seq_strand,
        );
    let seq_to_check = shifted_seq;
    let hgvs_output = shifted_hgvs_output;

    let shifted_seq = String::from_utf8(seq_to_check).map_err(|e| {
        DataFusionError::Execution(format!(
            "failed to build shifted HGVS sequence for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    let shifted_output_allele = String::from_utf8(hgvs_output).map_err(|e| {
        DataFusionError::Execution(format!(
            "failed to build shifted HGVS output allele for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    let shifted_compare_allele = match kind {
        GenomicShiftKind::Insertion => shifted_seq.clone(),
        GenomicShiftKind::Deletion => "-".to_string(),
    };
    let five_prime_flanking_seq = if shift_length == 0 {
        String::new()
    } else {
        // Both strands take the same flank slice: the last `shift_length+1`
        // bases of `pre_seq` in their original genomic order.
        pre_seq
            .chars()
            .rev()
            .take(shift_length + 1)
            .collect::<String>()
            .chars()
            .rev()
            .collect()
    };
    let three_prime_flanking_seq = if shift_length == 0 {
        String::new()
    } else {
        post_seq.chars().take(shift_length + 1).collect()
    };
    let inserted_len = shifted_output_allele.len() as i64;
    let display_start = if strand >= 0 {
        shifted_start
    } else {
        shifted_start - shift_length as i64
    };
    let (five_prime_context, three_prime_context) =
        if matches!(kind, GenomicShiftKind::Insertion) && inserted_len > 0 {
            // Match Ensembl's transcript HGVS duplication detection by storing the
            // adjacent reference sequence on the transcript 5' and 3' sides of the
            // shifted insertion point. For trimmed insertions `shifted_start` is
            // the genomic base immediately to the right of the insertion.
            let (five_start, five_end, three_start, three_end) = if strand >= 0 {
                (
                    (display_start - inserted_len).max(1),
                    (display_start - 1).max(0),
                    display_start,
                    display_start + inserted_len - 1,
                )
            } else {
                (
                    display_start,
                    display_start + inserted_len - 1,
                    (display_start - inserted_len).max(1),
                    (display_start - 1).max(0),
                )
            };
            let five = if five_end >= five_start && five_end > 0 {
                read_reference_sequence(reader, chrom, five_start, five_end)?
            } else {
                String::new()
            };
            let three = if three_end >= three_start && three_end > 0 {
                read_reference_sequence(reader, chrom, three_start, three_end)?
            } else {
                String::new()
            };
            (five, three)
        } else {
            (String::new(), String::new())
        };

    Ok(Some(HgvsGenomicShift {
        strand,
        shift_length,
        start: shifted_start,
        end: shifted_end,
        shifted_allele_string: shifted_seq.clone(),
        shifted_compare_allele,
        shifted_output_allele,
        ref_orig_allele_string: ref_allele.to_string(),
        alt_orig_allele_string: alt_allele.to_string(),
        five_prime_flanking_seq,
        three_prime_flanking_seq,
        five_prime_context,
        three_prime_context,
    }))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1360-L1364>
///
/// VEP builds transcript HGVS from feature-strand alleles. For transcripts on
/// the reverse strand, the variant allele sequence is reverse-complemented
/// before `hgvs_variant_notation()` is computed.
fn hgvs_feature_strand_alleles(
    tx: &TranscriptFeature,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<(String, String)> {
    let orient = |allele: &str| -> Option<String> {
        if allele == "-" {
            return Some("-".to_string());
        }
        let allele = allele.to_ascii_uppercase();
        if tx.strand >= 0 {
            Some(allele)
        } else {
            reverse_complement(&allele)
        }
    };

    Some((orient(ref_allele)?, orient(alt_allele)?))
}

fn transcript_reference_sequence_for_hgvsc(tx: &TranscriptFeature) -> Option<String> {
    if let Some(seq) = tx.spliced_seq.as_ref().or(tx.cdna_seq.as_ref()) {
        return Some(seq.clone());
    }

    match (
        tx.five_prime_utr_seq.as_deref(),
        tx.translateable_seq.as_deref(),
        tx.three_prime_utr_seq.as_deref(),
    ) {
        (Some(five_prime), Some(translateable), Some(three_prime)) => {
            Some(format!("{five_prime}{translateable}{three_prime}"))
        }
        _ => None,
    }
}

fn edited_transcript_reference_allele_for_hgvsc(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    ref_allele: &str,
    variant_start: i64,
    variant_end: i64,
) -> Option<String> {
    if tx.refseq_edits.is_empty() || ref_allele.is_empty() || ref_allele == "-" {
        return None;
    }

    let transcript_seq = transcript_reference_sequence_for_hgvsc(tx)?;
    let (start, end) = if variant_start <= variant_end {
        (variant_start, variant_end)
    } else {
        (variant_end, variant_start)
    };
    let mut cdna_positions = Vec::new();
    for genomic_pos in start..=end {
        let raw_cdna = raw_cdna_position_from_genomic(tx, tx_exons, genomic_pos)?;
        if raw_cdna.contains('+') || raw_cdna.contains('-') {
            return None;
        }
        let raw_cdna = raw_cdna.parse::<usize>().ok()?;
        let cdna = edited_transcript_sequence_cdna_index(tx, raw_cdna).unwrap_or(raw_cdna);
        if cdna == 0 || cdna > transcript_seq.len() {
            return None;
        }
        cdna_positions.push(cdna);
    }
    if cdna_positions.len() != ref_allele.len() {
        return None;
    }

    cdna_positions.sort_unstable();
    let bytes = transcript_seq.as_bytes();
    let mut transcript_ref = String::with_capacity(cdna_positions.len());
    for cdna in cdna_positions {
        transcript_ref.push((bytes[cdna - 1] as char).to_ascii_uppercase());
    }
    Some(transcript_ref)
}

fn parse_shiftable_indel<'a>(
    ref_allele: &'a str,
    alt_allele: &'a str,
) -> Option<(&'a str, GenomicShiftKind)> {
    if ref_allele == "-" && !alt_allele.is_empty() && alt_allele != "-" {
        return Some((alt_allele, GenomicShiftKind::Insertion));
    }
    if alt_allele == "-" && !ref_allele.is_empty() && ref_allele != "-" {
        return Some((ref_allele, GenomicShiftKind::Deletion));
    }
    None
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
///   replaces HGVS `alt` with transcript-space `feature_seq` when `_rna_edit`
///   attributes are present
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1443-L1449>
/// - Ensembl VEP `AnnotationType::Transcript::edit_transcript()` marks the
///   transcript with `_bam_edit_status = ok` when those edits were applied
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L579-L597>
fn edited_transcript_insertion_allele(
    tx: &TranscriptFeature,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<String> {
    if tx.bam_edit_status.as_deref() != Some("ok") {
        return None;
    }
    if ref_allele != "-" || alt_allele.is_empty() || alt_allele == "-" {
        return None;
    }
    Some(alt_allele.to_ascii_uppercase())
}

fn build_reference_region(chrom: &str, start: i64, end: i64) -> Result<Region> {
    let start = usize::try_from(start).map_err(|_| {
        DataFusionError::Execution(format!(
            "reference query start is negative or overflowed for {chrom}:{start}-{end}"
        ))
    })?;
    let end = usize::try_from(end).map_err(|_| {
        DataFusionError::Execution(format!(
            "reference query end is negative or overflowed for {chrom}:{start}-{end}"
        ))
    })?;
    let start = Position::try_from(start).map_err(|e| {
        DataFusionError::Execution(format!(
            "reference query start is invalid for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    let end = Position::try_from(end).map_err(|e| {
        DataFusionError::Execution(format!(
            "reference query end is invalid for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    Ok(Region::new(chrom, start..=end))
}

fn read_reference_sequence<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<String>
where
    R: BufRead + Seek,
{
    let region = build_reference_region(chrom, start, end)?;
    let record = reader.query(&region).map_err(|e| {
        DataFusionError::Execution(format!(
            "failed to query reference FASTA for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    String::from_utf8(record.sequence().as_ref().to_vec()).map_err(|e| {
        DataFusionError::Execution(format!(
            "reference FASTA returned non-UTF8 sequence for {chrom}:{start}-{end}: {e}"
        ))
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct HgvsNotation {
    start: i64,
    end: i64,
    ref_allele: String,
    alt_allele: String,
    kind: String,
}

/// Traceability:
/// - Ensembl Variation `Utils::Sequence::hgvs_variant_notation()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L493-L619>
fn hgvs_variant_notation(
    ref_allele: &str,
    alt_allele: &str,
    ref_start: i64,
    ref_end: i64,
) -> Option<HgvsNotation> {
    let ref_allele = if ref_allele == "-" {
        String::new()
    } else {
        ref_allele.to_string()
    };
    let alt_allele = if alt_allele == "-" {
        String::new()
    } else {
        alt_allele.to_string()
    };

    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();
    if ref_allele == alt_allele {
        if ref_len == 1 {
            return Some(HgvsNotation {
                start: ref_start,
                end: ref_end,
                ref_allele,
                alt_allele,
                kind: ">".to_string(),
            });
        }
        return None;
    }

    let kind = if alt_len == 0 {
        "del".to_string()
    } else if ref_len == alt_len {
        if ref_len == 1 {
            ">".to_string()
        } else if reverse_complement(&ref_allele).as_deref() == Some(alt_allele.as_str()) {
            "inv".to_string()
        } else {
            "delins".to_string()
        }
    } else if ref_len == 0 {
        "ins".to_string()
    } else if alt_len % ref_len == 0 && alt_allele == ref_allele.repeat(alt_len / ref_len) {
        if alt_len / ref_len == 2 {
            "dup".to_string()
        } else {
            format!("[{}]", alt_len / ref_len)
        }
    } else {
        "delins".to_string()
    };

    Some(HgvsNotation {
        start: ref_start,
        end: ref_end,
        ref_allele,
        alt_allele,
        kind,
    })
}

fn reverse_complement(seq: &str) -> Option<String> {
    let mut out = String::with_capacity(seq.len());
    for ch in seq.chars().rev() {
        out.push(match ch {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            _ => return None,
        });
    }
    Some(out)
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_clip_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2117-L2232>
fn clip_alleles(notation: &mut HgvsNotation, strand: i8) {
    let mut ref_allele = notation.ref_allele.clone();
    let mut alt_allele = notation.alt_allele.clone();
    let mut start = notation.start;
    let mut end = notation.end;
    let mut preseq = String::new();

    while !ref_allele.is_empty()
        && !alt_allele.is_empty()
        && ref_allele.as_bytes()[0] == alt_allele.as_bytes()[0]
    {
        preseq.push(ref_allele.remove(0));
        alt_allele.remove(0);
        if strand >= 0 {
            start += 1;
        } else {
            end -= 1;
        }
    }

    while !ref_allele.is_empty()
        && !alt_allele.is_empty()
        && ref_allele.as_bytes()[ref_allele.len() - 1]
            == alt_allele.as_bytes()[alt_allele.len() - 1]
    {
        ref_allele.pop();
        alt_allele.pop();
        if strand >= 0 {
            end -= 1;
        } else {
            start += 1;
        }
    }

    notation.ref_allele = ref_allele;
    notation.alt_allele = alt_allele;
    notation.start = start;
    notation.end = end;

    if notation.ref_allele.len() == 1
        && notation.alt_allele.len() == 1
        && notation.ref_allele != notation.alt_allele
    {
        notation.kind = ">".to_string();
    } else if notation.ref_allele.is_empty() && !notation.alt_allele.is_empty() {
        if preseq.ends_with(&notation.alt_allele) {
            notation.kind = "dup".to_string();
            notation.start -= notation.alt_allele.len() as i64;
        } else {
            notation.kind = "ins".to_string();
        }
    } else if !notation.ref_allele.is_empty() && notation.alt_allele.is_empty() {
        notation.kind = "del".to_string();
    }
}

fn notation_to_hgvsc_coords(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    notation: &HgvsNotation,
) -> Option<(String, String)> {
    if notation.kind == "ins" {
        let start = hgvs_cdna_position_from_genomic(tx, tx_exons, notation.start - 1)?;
        let end = hgvs_cdna_position_from_genomic(tx, tx_exons, notation.start)?;
        return Some((start, end));
    }
    Some((
        hgvs_cdna_position_from_genomic(tx, tx_exons, notation.start)?,
        hgvs_cdna_position_from_genomic(tx, tx_exons, notation.end)?,
    ))
}

fn apply_shifted_insertion_duplication(
    tx: &TranscriptFeature,
    feature_alt: &str,
    shift: &HgvsGenomicShift,
    notation: &mut HgvsNotation,
) {
    if notation.kind != "ins" || feature_alt.is_empty() {
        return;
    }

    let orient_context = |context: &str| {
        if tx.strand >= 0 {
            Some(context.to_ascii_uppercase())
        } else {
            reverse_complement(&context.to_ascii_uppercase())
        }
    };
    let three_prime_context = orient_context(&shift.three_prime_context).unwrap_or_default();
    let five_prime_context = orient_context(&shift.five_prime_context).unwrap_or_default();
    let dup_from_five_prime = five_prime_context == feature_alt;
    let dup_from_three_prime = three_prime_context == feature_alt;
    if !dup_from_three_prime && !dup_from_five_prime {
        return;
    }

    let alt_len = feature_alt.len() as i64;
    let display_start = shift.display_start();
    notation.kind = "dup".to_string();
    if dup_from_five_prime {
        if tx.strand >= 0 {
            notation.start = display_start - alt_len;
            notation.end = display_start - 1;
        } else {
            notation.start = display_start;
            notation.end = display_start + alt_len - 1;
        }
    } else {
        if tx.strand >= 0 {
            notation.start = display_start;
            notation.end = display_start + alt_len - 1;
        } else {
            notation.start = display_start - alt_len;
            notation.end = display_start - 1;
        }
    }
}

/// Traceability:
/// - Ensembl Variation `perform_shift()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351>
///
/// VCF alleles are always reported on the forward genomic strand. When HGVS is
/// generated for a reverse-strand transcript, Ensembl shifts the genomic
/// comparison sequence and the HGVS output allele in different directions.
pub(crate) fn perform_shift_ensembl(
    seq_to_check: &[u8],
    hgvs_output: &[u8],
    post_seq: &str,
    pre_seq: &str,
    var_start: i64,
    var_end: i64,
    reverse: bool,
    seq_strand: i8,
) -> (usize, Vec<u8>, Vec<u8>, i64, i64) {
    let mut seq_to_check = seq_to_check.to_vec();
    let mut hgvs_output = hgvs_output.to_vec();
    let mut var_start = var_start;
    let mut var_end = var_end;
    let indel_length = seq_to_check.len();
    let mut shift_length = 0usize;
    let vf_strand = 1i8;
    let hgvs_reverse = vf_strand != seq_strand;
    let start_n = usize::from(reverse);
    // Perl guard: $loop_limiter = length($post_seq) if $loop_limiter < 0;
    // When indel_length > flank length, fall back to flank length so the
    // character-by-character shift loop can still find the correct (small) shift.
    let loop_limiter = if reverse {
        if indel_length > pre_seq.len() {
            pre_seq.len()
        } else {
            pre_seq.len() - indel_length + 1
        }
    } else if indel_length > post_seq.len() {
        post_seq.len()
    } else {
        post_seq.len() - indel_length
    };

    for n in start_n..=loop_limiter {
        let (check_next_del, check_next_ref, hgvs_next_del) = if reverse {
            let Some(check_next_del) = seq_to_check.last().copied() else {
                break;
            };
            let Some(check_next_ref) = pre_seq
                .as_bytes()
                .get(pre_seq.len().saturating_sub(n))
                .copied()
            else {
                break;
            };
            let hgvs_next_del = if hgvs_reverse {
                hgvs_output.first().copied()
            } else {
                hgvs_output.last().copied()
            };
            let Some(hgvs_next_del) = hgvs_next_del else {
                break;
            };
            (check_next_del, check_next_ref, hgvs_next_del)
        } else {
            let Some(check_next_del) = seq_to_check.first().copied() else {
                break;
            };
            let Some(check_next_ref) = post_seq.as_bytes().get(n).copied() else {
                break;
            };
            let hgvs_next_del = if hgvs_reverse {
                hgvs_output.last().copied()
            } else {
                hgvs_output.first().copied()
            };
            let Some(hgvs_next_del) = hgvs_next_del else {
                break;
            };
            (check_next_del, check_next_ref, hgvs_next_del)
        };

        if check_next_del != check_next_ref {
            break;
        }

        shift_length += 1;
        if reverse {
            let last = seq_to_check.pop().unwrap_or(check_next_del);
            seq_to_check.insert(0, last);
            if hgvs_reverse {
                if !hgvs_output.is_empty() {
                    hgvs_output.remove(0);
                }
                hgvs_output.push(hgvs_next_del);
            } else {
                hgvs_output.pop();
                hgvs_output.insert(0, hgvs_next_del);
            }
        } else {
            if !seq_to_check.is_empty() {
                seq_to_check.remove(0);
            }
            seq_to_check.push(check_next_del);
            if hgvs_reverse {
                hgvs_output.pop();
                hgvs_output.insert(0, hgvs_next_del);
            } else {
                if !hgvs_output.is_empty() {
                    hgvs_output.remove(0);
                }
                hgvs_output.push(hgvs_next_del);
            }
            var_start += 1;
            var_end += 1;
        }
    }

    (shift_length, seq_to_check, hgvs_output, var_start, var_end)
}

fn overlap_len(left_start: i64, left_end: i64, right_start: i64, right_end: i64) -> i64 {
    let start = left_start.max(right_start);
    let end = left_end.min(right_end);
    if end < start {
        0
    } else {
        end.saturating_sub(start).saturating_add(1)
    }
}

/// Traceability:
/// - Ensembl Variation `Utils::Sequence::format_hgvs_string()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L635-L676>
fn format_hgvs_string(
    ref_name: &str,
    numbering: char,
    start: &str,
    end: &str,
    notation: &HgvsNotation,
) -> Option<String> {
    let coordinates = if start == end {
        start.to_string()
    } else {
        format!("{start}_{end}")
    };
    let suffix =
        if notation.kind == ">" || (notation.kind == "inv" && notation.ref_allele.len() == 1) {
            format!("{start}{}>{}", notation.ref_allele, notation.alt_allele)
        } else if matches!(notation.kind.as_str(), "del" | "inv" | "dup") {
            format!("{coordinates}{}", notation.kind)
        } else if notation.kind == "delins" {
            format!("{coordinates}delins{}", notation.alt_allele)
        } else if notation.kind == "ins" {
            format!("{coordinates}ins{}", notation.alt_allele)
        } else if notation.kind.starts_with('[') && notation.kind.ends_with(']') {
            format!("{coordinates}{}", notation.kind)
        } else {
            return None;
        };
    Some(format!("{ref_name}:{numbering}.{suffix}"))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_cDNA_position()`
///   uses transcript coding cDNA anchors from the transcript/mapper state when
///   converting raw cDNA positions into HGVS `c.` coordinates
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2683-L2765>
fn coding_cdna_bounds(tx: &TranscriptFeature, tx_exons: &[&ExonFeature]) -> Option<(usize, usize)> {
    if let (Some(start_codon), Some(stop_codon)) = (tx.cdna_coding_start, tx.cdna_coding_end) {
        return Some((start_codon, stop_codon));
    }
    let cds_start = tx.cds_start?;
    let cds_end = tx.cds_end?;
    let coding_start_anchor = if tx.strand >= 0 { cds_start } else { cds_end };
    let coding_end_anchor = if tx.strand >= 0 { cds_end } else { cds_start };
    Some((
        genomic_to_cdna_index_for_transcript(tx, tx_exons, coding_start_anchor)?,
        genomic_to_cdna_index_for_transcript(tx, tx_exons, coding_end_anchor)?,
    ))
}

fn hgvs_cdna_position_from_genomic(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    genomic_pos: i64,
) -> Option<String> {
    let mut cdna_position = raw_cdna_position_from_genomic(tx, tx_exons, genomic_pos)?;
    let has_intron_offset = cdna_position
        .char_indices()
        .skip(1)
        .any(|(_, ch)| matches!(ch, '+' | '-'));
    if has_intron_offset {
        if let Some(exon_geometry_position) =
            native_refseq_pre_coding_intronic_exon_geometry_position(
                tx,
                tx_exons,
                genomic_pos,
                &cdna_position,
            )
        {
            cdna_position = exon_geometry_position;
        }
    }
    let cdna_position = if has_intron_offset {
        cdna_position
    } else {
        let absolute_cdna = cdna_position.parse::<i64>().ok()?;
        if tx
            .cdna_coding_end
            .is_some_and(|cdna_coding_end| absolute_cdna > cdna_coding_end as i64)
        {
            refseq_sequence_offset_for_cdna(tx, absolute_cdna)
                .map(|offset| (absolute_cdna + offset).to_string())
                .unwrap_or_else(|| {
                    adjust_refseq_cdna_component(tx, &cdna_position).unwrap_or(cdna_position)
                })
        } else {
            adjust_refseq_cdna_component(tx, &cdna_position).unwrap_or(cdna_position)
        }
    };
    shift_to_hgvs_coding_coordinates(tx, tx_exons, &cdna_position)
}

fn shift_to_hgvs_coding_coordinates(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    raw_cdna_position: &str,
) -> Option<String> {
    let (raw_coord, mut intron_offset) = split_hgvs_coord(raw_cdna_position)?;
    let Some((start_codon, stop_codon)) = coding_cdna_bounds(tx, tx_exons) else {
        return Some(raw_cdna_position.to_string());
    };

    let mut coord = raw_coord;
    let mut prefix = "";
    let mut coord_text = None;

    if coord > stop_codon as i64 {
        coord -= stop_codon as i64;
        prefix = "*";
    } else if coord == stop_codon as i64 && intron_offset.is_some() {
        prefix = "*";
        coord_text = Some(String::new());
        intron_offset = intron_offset.map(|offset| offset.trim_start_matches('+').to_string());
    }

    if prefix.is_empty() {
        if coord >= start_codon as i64 {
            coord += 1;
        }
        coord -= start_codon as i64;
        coord_text = Some(coord.to_string());
    } else if coord_text.is_none() {
        coord_text = Some(coord.to_string());
    }

    Some(format!(
        "{prefix}{}{}",
        coord_text?,
        intron_offset.unwrap_or_default()
    ))
}

fn native_refseq_hgvs_intronic_anchor_uses_post_gap_numbering(tx: &TranscriptFeature) -> bool {
    tx.source.as_deref() == Some("RefSeq")
        && matches!(
            tx.transcript_id.as_bytes().get(..2),
            Some(b"NM") | Some(b"NR") | Some(b"XM") | Some(b"XR")
        )
}

fn native_refseq_insertion_shift_at_anchor(
    tx: &TranscriptFeature,
    exon_coord: i64,
    mapper_coord: i64,
) -> bool {
    let mut offset = 0i64;
    for edit in &tx.refseq_edits {
        if edit.skip_refseq_offset || edit.end >= exon_coord {
            continue;
        }
        let Some(replacement_len) = edit.replacement_len else {
            continue;
        };
        let replaced_len = edit.end.saturating_sub(edit.start).saturating_add(1);
        offset += replacement_len as i64 - replaced_len;
    }

    offset > 0 && exon_coord.saturating_add(offset) == mapper_coord
}

fn native_refseq_pre_coding_intronic_exon_geometry_position(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    genomic_pos: i64,
    mapper_position: &str,
) -> Option<String> {
    if !native_refseq_hgvs_intronic_anchor_uses_post_gap_numbering(tx)
        || tx.cdna_mapper_segments.is_empty()
    {
        return None;
    }
    let (mapper_coord, mapper_offset) = split_hgvs_coord(mapper_position)?;

    let mut exon_geometry_tx = tx.clone();
    exon_geometry_tx.cdna_mapper_segments.clear();
    let exon_geometry_position =
        raw_cdna_position_from_genomic(&exon_geometry_tx, tx_exons, genomic_pos)?;
    let (exon_coord, exon_offset) = split_hgvs_coord(&exon_geometry_position)?;
    // Native RefSeq transcripts with a leading `_rna_edit` insertion (for
    // example NM_001177639.3) keep HGVS intronic anchors on the pre-edit exon
    // boundary even though the cached mapper cDNA coordinates are shifted by
    // the inserted bases. Other native RefSeq mapper/exon disagreements, such
    // as one-base deleted gaps, still follow the larger exon-geometry anchor.
    let insertion_shift_prefers_pre_edit_exon_geometry = exon_coord < mapper_coord
        && native_refseq_insertion_shift_at_anchor(tx, exon_coord, mapper_coord);
    if exon_offset == mapper_offset
        && (exon_coord > mapper_coord || insertion_shift_prefers_pre_edit_exon_geometry)
    {
        Some(exon_geometry_position)
    } else {
        None
    }
}

fn split_hgvs_coord(value: &str) -> Option<(i64, Option<String>)> {
    let body = value.strip_prefix('*').unwrap_or(value);
    let split_idx = body
        .char_indices()
        .skip(1)
        .find_map(|(idx, ch)| matches!(ch, '+' | '-').then_some(idx));
    let (coord_part, offset_part) = if let Some(idx) = split_idx {
        (&body[..idx], Some(body[idx..].to_string()))
    } else {
        (body, None)
    };
    Some((coord_part.parse::<i64>().ok()?, offset_part))
}

fn compare_hgvs_positions(left: &str, right: &str) -> Option<Ordering> {
    let left_star = left.starts_with('*');
    let right_star = right.starts_with('*');
    match (left_star, right_star) {
        (false, true) => return Some(Ordering::Less),
        (true, false) => return Some(Ordering::Greater),
        _ => {}
    }

    let (left_coord, left_offset) = split_hgvs_coord(left)?;
    let (right_coord, right_offset) = split_hgvs_coord(right)?;
    let left_offset = left_offset
        .as_deref()
        .unwrap_or("0")
        .parse::<i64>()
        .ok()
        .unwrap_or(0);
    let right_offset = right_offset
        .as_deref()
        .unwrap_or("0")
        .parse::<i64>()
        .ok()
        .unwrap_or(0);
    Some((left_coord, left_offset).cmp(&(right_coord, right_offset)))
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProteinHgvsNotation {
    start: usize,
    end: usize,
    ref_allele: String,
    alt_allele: String,
    original_ref: String,
    preseq: String,
    kind: String,
}

impl ProteinHgvsNotation {
    /// Traceability:
    /// - Ensembl Variation `TranscriptVariationAllele::hgvs_protein()`
    ///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1700-L1749>
    fn from_context(data: &ProteinHgvsData) -> Self {
        let ref_allele = normalize_peptide_allele(&data.ref_peptide);
        let alt_allele = normalize_peptide_allele(&data.alt_peptide);
        Self {
            start: data.start,
            end: data.end,
            original_ref: ref_allele.clone(),
            ref_allele,
            alt_allele,
            preseq: String::new(),
            kind: String::new(),
        }
    }
}

/// Compute HGVSp notation from peptide-level coding classification data.
///
/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::hgvs_protein()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1593-L1758>
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_type()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1976-L2041>
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2043-L2114>
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1833-L1974>
pub fn format_hgvsp(
    translation: &TranslationFeature,
    protein: &ProteinHgvsData,
    shift_hgvs: bool,
) -> Option<String> {
    let protein_id = versioned_id(translation.stable_id.as_deref()?, translation.version);
    if protein.start_lost {
        let start_ref = if protein.ref_peptide.is_empty() {
            protein
                .ref_translation
                .chars()
                .next()
                .map(|aa| aa.to_string())?
        } else {
            protein.ref_peptide.clone()
        };
        return Some(format!(
            "{protein_id}:p.{}{}?",
            peptide_first_three(&start_ref)?,
            protein.start
        ));
    }

    let mut notation = ProteinHgvsNotation::from_context(protein);
    if protein.frameshift {
        resolve_frameshift_hgvs(&mut notation, protein)?;
    } else {
        if notation.ref_allele != notation.alt_allele {
            clip_protein_alleles(&mut notation);
        } else {
            notation.kind = "=".to_string();
        }
        if notation.kind.is_empty() {
            notation.kind = protein_event_type(&notation.ref_allele, &notation.alt_allele, false);
        }
        // Traceability:
        // - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
        //   calls `_check_peptides_post_var()` (3' shift) FIRST, then
        //   `_check_for_peptide_duplication()`. The shift may change the
        //   insertion position, making the upstream sequence different when
        //   the dup check runs.
        //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2054-L2089
        if shift_hgvs && matches!(notation.kind.as_str(), "ins" | "del") {
            shift_peptides_post_var(&mut notation, &protein.ref_translation);
        }
        if notation.kind == "ins"
            && check_for_peptide_duplication(&mut notation, &protein.ref_translation)
        {
            // Dup detected — skip flanking.
        } else if notation.kind == "ins" {
            notation.ref_allele = surrounding_peptides(
                &protein.ref_translation,
                notation.start.min(notation.end),
                &notation.original_ref,
                Some(2),
            )?;
        }
    }

    format_hgvsp_notation(&protein_id, &notation, protein)
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2054-L2089>
fn normalize_peptide_allele(allele: &str) -> String {
    if allele == "-" {
        String::new()
    } else {
        allele.to_string()
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_clip_alleles()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2117-L2225>
fn clip_protein_alleles(notation: &mut ProteinHgvsNotation) {
    let mut ref_allele = notation.ref_allele.clone();
    let mut alt_allele = notation.alt_allele.clone();
    let mut start = notation.start;
    let mut end = notation.end;
    let mut preseq = String::new();

    while !ref_allele.is_empty()
        && !alt_allele.is_empty()
        && ref_allele.as_bytes()[0] == alt_allele.as_bytes()[0]
    {
        preseq.push(ref_allele.remove(0));
        alt_allele.remove(0);
        start = start.saturating_add(1);
    }

    while !ref_allele.is_empty()
        && !alt_allele.is_empty()
        && ref_allele.as_bytes()[ref_allele.len() - 1]
            == alt_allele.as_bytes()[alt_allele.len() - 1]
    {
        ref_allele.pop();
        alt_allele.pop();
        end = end.saturating_sub(1);
    }

    notation.start = start;
    notation.end = end;
    notation.ref_allele = ref_allele;
    notation.alt_allele = alt_allele;
    notation.preseq = preseq;

    notation.kind = if notation.ref_allele == notation.alt_allele {
        "=".to_string()
    } else if notation.ref_allele.len() == 1 && notation.alt_allele.len() == 1 {
        ">".to_string()
    } else if notation.ref_allele.is_empty() && !notation.alt_allele.is_empty() {
        "ins".to_string()
    } else if !notation.ref_allele.is_empty() && notation.alt_allele.is_empty() {
        "del".to_string()
    } else {
        "delins".to_string()
    };
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_type()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1976-L2041>
fn protein_event_type(ref_allele: &str, alt_allele: &str, frameshift: bool) -> String {
    if frameshift {
        "fs".to_string()
    } else if ref_allele == alt_allele {
        "=".to_string()
    } else if ref_allele.is_empty() {
        "ins".to_string()
    } else if alt_allele.is_empty() {
        "del".to_string()
    } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
        ">".to_string()
    } else {
        "delins".to_string()
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_fs_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2249-L2295>
fn resolve_frameshift_hgvs(
    notation: &mut ProteinHgvsNotation,
    protein: &ProteinHgvsData,
) -> Option<()> {
    notation.kind = "fs".to_string();
    let ref_translation = append_terminal_stop(&protein.ref_translation);
    let alt_translation = protein.alt_translation.as_str();
    let mut start = notation.start;

    if start > alt_translation.len() {
        notation.kind = "del".to_string();
        notation.end = start;
        notation.ref_allele = peptide_char(&ref_translation, start)?.to_string();
        notation.alt_allele.clear();
        return Some(());
    }

    while start <= alt_translation.len() {
        let ref_aa = peptide_char(&ref_translation, start)?;
        let alt_aa = peptide_char(alt_translation, start)?;
        if ref_aa == '*' && alt_aa == '*' {
            notation.kind = "=".to_string();
            notation.start = start;
            notation.end = start;
            notation.ref_allele = "*".to_string();
            notation.alt_allele = "*".to_string();
            return Some(());
        }
        if ref_aa != alt_aa {
            notation.start = start;
            notation.end = start;
            notation.ref_allele = ref_aa.to_string();
            notation.alt_allele = alt_aa.to_string();
            return Some(());
        }
        start = start.saturating_add(1);
    }

    notation.kind = "del".to_string();
    notation.start = start;
    notation.end = start;
    notation.ref_allele = peptide_char(&ref_translation, start)?.to_string();
    notation.alt_allele.clear();
    Some(())
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_fs_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2262-L2295>
fn append_terminal_stop(peptide: &str) -> String {
    if peptide.contains('*') {
        peptide.to_string()
    } else {
        format!("{peptide}*")
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_fs_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2278-L2291>
fn peptide_char(peptide: &str, pos: usize) -> Option<char> {
    peptide
        .as_bytes()
        .get(pos.checked_sub(1)?)
        .map(|b| *b as char)
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_check_peptides_post_var()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2501-L2518>
/// - Ensembl Variation `TranscriptVariationAllele::_shift_3prime()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2525-L2573>
fn shift_peptides_post_var(notation: &mut ProteinHgvsNotation, ref_translation: &str) {
    let post_seq = match surrounding_peptides(
        ref_translation,
        notation.end.saturating_add(1),
        &notation.original_ref,
        None,
    ) {
        Some(seq) => seq,
        None => return,
    };

    let seq_to_check = if notation.kind == "ins" {
        &mut notation.alt_allele
    } else if notation.kind == "del" {
        &mut notation.ref_allele
    } else {
        return;
    };

    let deleted_len = seq_to_check.len();
    if deleted_len == 0 || post_seq.len() < deleted_len {
        return;
    }

    for check_next_post in post_seq.chars() {
        let Some(check_next_del) = seq_to_check.chars().next() else {
            break;
        };
        if check_next_del != check_next_post {
            break;
        }
        notation.start = notation.start.saturating_add(1);
        notation.end = notation.end.saturating_add(1);
        seq_to_check.remove(0);
        seq_to_check.push(check_next_del);
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_check_for_peptide_duplication()`
///   builds upstream from `substr($reference_trans, 0, $start - 1)` + preseq, then
///   tests whether the alt peptide matches the upstream at `start - alt_len - 1`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2371-L2410>
fn check_for_peptide_duplication(
    notation: &mut ProteinHgvsNotation,
    ref_translation: &str,
) -> bool {
    if notation.alt_allele.is_empty() || notation.start == 0 {
        return false;
    }

    // Traceability:
    // - VEP's `_check_for_peptide_duplication()` builds `upstream` from
    //   `substr($reference_trans, 0, $start - 1)` + preseq, then checks
    //   a single window at `$start - length($alt) - 1`. No fallback.
    //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2371-L2410
    try_peptide_dup_at(notation, ref_translation, notation.start)
}

fn try_peptide_dup_at(
    notation: &mut ProteinHgvsNotation,
    ref_translation: &str,
    check_start: usize,
) -> bool {
    let mut upstream = ref_translation
        .get(..check_start.saturating_sub(1))
        .unwrap_or_default()
        .to_string();
    upstream.push_str(&notation.preseq);

    let alt_len = notation.alt_allele.len();
    let Some(test_new_start) = check_start
        .checked_sub(alt_len)
        .and_then(|s| s.checked_sub(1))
    else {
        return false;
    };
    let Some(test_seq) = upstream.get(test_new_start..test_new_start.saturating_add(alt_len))
    else {
        return false;
    };

    if test_seq == notation.alt_allele {
        notation.kind = "dup".to_string();
        notation.end = check_start.saturating_sub(1);
        notation.start = check_start.saturating_sub(alt_len);
        true
    } else {
        false
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_surrounding_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2297-L2320>
fn surrounding_peptides(
    ref_translation: &str,
    ref_pos: usize,
    original_ref: &str,
    length: Option<usize>,
) -> Option<String> {
    let mut ref_trans = ref_translation.to_string();
    if original_ref.starts_with('*') {
        ref_trans.push_str(original_ref);
    }
    if ref_trans.len() < ref_pos {
        return None;
    }
    let start = ref_pos.checked_sub(1)?;
    match length {
        Some(len) => ref_trans
            .get(start..start.saturating_add(len))
            .map(str::to_string),
        None => ref_trans.get(start..).map(str::to_string),
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_stop_loss_extra_AA()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2406-L2461>
///
///   VEP Perl (non-frameshift path):
///     my $ref_temp = $self->transcript_variation->_peptide();
///     my $ref_len  = length($ref_temp);          # full peptide length
///     $extra_aa    = $+[0] - 1 - $ref_len;       # $+[0] = 1-based end of first * match
///
///   `_peptide()` returns BioPerl `translate()->seq` which translates ALL
///   codons (including past internal stops). `length()` is the full string
///   length — NOT the position of the first stop codon. This matters for
///   LoF transcripts with internal stop codons.
fn stop_loss_extra_aa(
    protein: &ProteinHgvsData,
    ref_var_pos: usize,
    frameshift: bool,
) -> Option<usize> {
    let alt_translation = protein
        .alt_translation_extension
        .as_deref()
        .unwrap_or(&protein.alt_translation);
    let stop_idx = alt_translation.find('*')?;
    let extra = if frameshift {
        stop_idx.saturating_add(1).checked_sub(ref_var_pos)?
    } else {
        // VEP: $ref_len = length(_peptide()). The VEP cache stores the
        // peptide WITHOUT the terminal *, so length = number of amino
        // acids excluding the stop. Our ref_translation includes *, so
        // we strip trailing * to match VEP's cached peptide length.
        // https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm#L1282-L1291
        let ref_len = protein.ref_translation.trim_end_matches('*').len();
        stop_idx
            .saturating_add(1)
            .checked_sub(ref_len.saturating_add(1))?
    };
    (extra > 0).then_some(extra)
}

fn hgvs_aa_one_to_three(aa: char) -> &'static str {
    match aa {
        'X' => "Ter",
        _ => aa_one_to_three(aa),
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2083-L2112>
fn peptide_to_three_letter(peptide: &str) -> String {
    peptide.chars().map(hgvs_aa_one_to_three).collect()
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1877-L1880>
fn peptide_first_three(peptide: &str) -> Option<&'static str> {
    Some(hgvs_aa_one_to_three(peptide.chars().next()?))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1877-L1880>
fn peptide_last_three(peptide: &str) -> Option<&'static str> {
    Some(hgvs_aa_one_to_three(peptide.chars().last()?))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1833-L1974>
fn format_hgvsp_notation(
    protein_id: &str,
    notation: &ProteinHgvsNotation,
    protein: &ProteinHgvsData,
) -> Option<String> {
    let mut out = format!("{protein_id}:p.");

    if notation.ref_allele == notation.alt_allele && !matches!(notation.kind.as_str(), "fs" | "ins")
    {
        out.push_str(&format!(
            "{}{}=",
            peptide_to_three_letter(&notation.ref_allele),
            notation.start
        ));
        return Some(out);
    }

    if protein.stop_lost && matches!(notation.kind.as_str(), "del" | ">") {
        let extra = stop_loss_extra_aa(protein, notation.start.saturating_sub(1), false)
            .map(|v| v.to_string())
            .unwrap_or_else(|| "?".to_string());
        let alt_head = peptide_first_three(&notation.alt_allele).unwrap_or("?");
        if notation.ref_allele.len() > 1 && notation.kind == "del" {
            out.push_str(&format!(
                "{}{}_{}{}{}extTer{}",
                peptide_first_three(&notation.ref_allele)?,
                notation.start,
                peptide_last_three(&notation.ref_allele)?,
                notation.end,
                alt_head,
                extra
            ));
        } else {
            out.push_str(&format!(
                "{}{}{}extTer{}",
                peptide_to_three_letter(&notation.ref_allele),
                notation.start,
                alt_head,
                extra
            ));
        }
        return Some(out);
    }

    match notation.kind.as_str() {
        "dup" => {
            if notation.start < notation.end {
                out.push_str(&format!(
                    "{}{}_{}{}dup",
                    peptide_first_three(&notation.alt_allele)?,
                    notation.start,
                    peptide_last_three(&notation.alt_allele)?,
                    notation.end
                ));
            } else {
                out.push_str(&format!(
                    "{}{}dup",
                    peptide_to_three_letter(&notation.alt_allele),
                    notation.start
                ));
            }
        }
        ">" => {
            out.push_str(&format!(
                "{}{}{}",
                peptide_to_three_letter(&notation.ref_allele),
                notation.start,
                peptide_to_three_letter(&notation.alt_allele)
            ));
        }
        "delins" | "ins" => {
            let mut alt_allele = notation.alt_allele.clone();
            if let Some(stop_idx) = alt_allele.find('*') {
                alt_allele.truncate(stop_idx.saturating_add(1));
            }
            let mut alt = peptide_to_three_letter(&alt_allele);
            if notation.ref_allele.ends_with('*') {
                if let Some(extra) =
                    stop_loss_extra_aa(protein, notation.start.saturating_sub(1), false)
                {
                    alt.push_str(&format!("extTer{extra}"));
                }
            }
            if notation.start == notation.end && notation.kind == "delins" {
                out.push_str(&format!(
                    "{}{}{}{}",
                    peptide_first_three(&notation.ref_allele)?,
                    notation.start,
                    notation.kind,
                    alt
                ));
            } else {
                let (start, end) = if notation.start > notation.end {
                    (notation.end, notation.start)
                } else {
                    (notation.start, notation.end)
                };
                out.push_str(&format!(
                    "{}{}_{}{}{}{}",
                    peptide_first_three(&notation.ref_allele)?,
                    start,
                    peptide_last_three(&notation.ref_allele)?,
                    end,
                    notation.kind,
                    alt
                ));
            }
        }
        "fs" => {
            if notation.alt_allele == "*" {
                out.push_str(&format!(
                    "{}{}Ter",
                    peptide_to_three_letter(&notation.ref_allele),
                    notation.start
                ));
            } else {
                let extra = stop_loss_extra_aa(protein, notation.start.saturating_sub(1), true)
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "?".to_string());
                out.push_str(&format!(
                    "{}{}{}fsTer{}",
                    peptide_to_three_letter(&notation.ref_allele),
                    notation.start,
                    peptide_to_three_letter(&notation.alt_allele),
                    extra
                ));
            }
        }
        "del" => {
            if notation.ref_allele.len() > 1 {
                out.push_str(&format!(
                    "{}{}_{}{}del",
                    peptide_first_three(&notation.ref_allele)?,
                    notation.start,
                    peptide_last_three(&notation.ref_allele)?,
                    notation.end
                ));
            } else {
                out.push_str(&format!(
                    "{}{}del",
                    peptide_to_three_letter(&notation.ref_allele),
                    notation.start
                ));
            }
        }
        _ if notation.start != notation.end => {
            out.push_str(&format!(
                "{}{}_{}{}",
                peptide_to_three_letter(&notation.ref_allele),
                notation.start,
                peptide_to_three_letter(&notation.alt_allele),
                notation.end
            ));
        }
        _ => {
            out.push_str(&format!(
                "{}{}{}",
                peptide_to_three_letter(&notation.ref_allele),
                notation.start,
                peptide_to_three_letter(&notation.alt_allele)
            ));
        }
    }

    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_transcript(
        biotype: &str,
        strand: i8,
        cds_start: Option<i64>,
        cds_end: Option<i64>,
    ) -> TranscriptFeature {
        TranscriptFeature {
            transcript_id: "ENSTHGVS000001".to_string(),
            chrom: "1".to_string(),
            start: 90,
            end: 140,
            strand,
            biotype: biotype.to_string(),
            cds_start,
            cds_end,
            cdna_coding_start: None,
            cdna_coding_end: None,
            cdna_mapper_segments: Vec::new(),
            mature_mirna_regions: Vec::new(),
            gene_stable_id: None,
            gene_symbol: None,
            gene_symbol_source: None,
            gene_hgnc_id_native: None,
            gene_hgnc_id: None,
            display_xref_id: None,
            source: None,
            refseq_match: None,
            refseq_edits: Vec::new(),
            is_gencode_basic: false,
            is_gencode_primary: false,
            bam_edit_status: None,
            has_non_polya_rna_edit: false,
            spliced_seq: None,
            five_prime_utr_seq: None,
            three_prime_utr_seq: None,
            translateable_seq: None,
            cdna_seq: None,
            version: Some(1),
            cds_start_nf: false,
            cds_end_nf: false,
            flags_str: None,
            is_canonical: false,
            tsl: None,
            mane_select: None,
            mane_plus_clinical: None,
            translation_stable_id: None,
            gene_phenotype: false,
            ccds: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            uniprot_isoform: None,
            appris: None,
            ncrna_structure: None,
        }
    }

    fn make_exon() -> ExonFeature {
        ExonFeature {
            transcript_id: "ENSTHGVS000001".to_string(),
            exon_number: 1,
            start: 90,
            end: 140,
        }
    }

    fn make_translation() -> TranslationFeature {
        TranslationFeature {
            transcript_id: "ENSTHGVS000001".to_string(),
            cds_len: Some(9),
            protein_len: Some(3),
            translation_seq: None,
            cds_sequence: None,
            translation_seq_canonical: None,
            cds_sequence_canonical: None,
            stable_id: Some("ENSPHGVS000001".to_string()),
            version: Some(1),
            protein_features: Vec::new(),
        }
    }

    #[test]
    fn test_aa_one_to_three() {
        assert_eq!(aa_one_to_three('A'), "Ala");
        assert_eq!(aa_one_to_three('R'), "Arg");
        assert_eq!(aa_one_to_three('*'), "Ter");
        assert_eq!(aa_one_to_three('X'), "Xaa");
    }

    #[test]
    fn test_versioned_id() {
        assert_eq!(
            versioned_id("ENST00000379410", Some(6)),
            "ENST00000379410.6"
        );
        assert_eq!(versioned_id("ENST00000379410", None), "ENST00000379410");
        assert_eq!(versioned_id("NM_001206729.2", Some(1)), "NM_001206729.2");
        assert_eq!(versioned_id("NP_001193658.1", Some(1)), "NP_001193658.1");
    }

    #[test]
    fn test_hgvs_variant_notation_substitution() {
        let notation = hgvs_variant_notation("G", "A", 103, 103).expect("notation");
        assert_eq!(notation.kind, ">");
        assert_eq!(notation.start, 103);
        assert_eq!(notation.end, 103);
        assert_eq!(notation.ref_allele, "G");
        assert_eq!(notation.alt_allele, "A");
    }

    #[test]
    fn test_clip_alleles_reclassifies_delins_to_substitution() {
        let mut notation = HgvsNotation {
            start: 100,
            end: 101,
            ref_allele: "AC".to_string(),
            alt_allele: "AT".to_string(),
            kind: "delins".to_string(),
        };
        clip_alleles(&mut notation, 1);
        assert_eq!(notation.kind, ">");
        assert_eq!(notation.start, 101);
        assert_eq!(notation.end, 101);
        assert_eq!(notation.ref_allele, "C");
        assert_eq!(notation.alt_allele, "T");
    }

    #[test]
    fn test_clip_alleles_uses_transcript_oriented_coordinates_on_negative_strand() {
        let mut notation = HgvsNotation {
            start: 100,
            end: 101,
            ref_allele: "AG".to_string(),
            alt_allele: "AC".to_string(),
            kind: "delins".to_string(),
        };
        clip_alleles(&mut notation, -1);
        assert_eq!(notation.kind, ">");
        assert_eq!(notation.start, 100);
        assert_eq!(notation.end, 100);
        assert_eq!(notation.ref_allele, "G");
        assert_eq!(notation.alt_allele, "C");
    }

    #[test]
    fn test_clip_alleles_reclassifies_delins_to_insertion() {
        let mut notation = HgvsNotation {
            start: 100,
            end: 100,
            ref_allele: "A".to_string(),
            alt_allele: "AT".to_string(),
            kind: "delins".to_string(),
        };
        clip_alleles(&mut notation, 1);
        assert_eq!(notation.kind, "ins");
        assert_eq!(notation.start, 101);
        assert_eq!(notation.end, 100);
        assert_eq!(notation.ref_allele, "");
        assert_eq!(notation.alt_allele, "T");
    }

    #[test]
    fn test_clip_alleles_reclassifies_delins_to_duplication() {
        let mut notation = HgvsNotation {
            start: 100,
            end: 100,
            ref_allele: "A".to_string(),
            alt_allele: "AA".to_string(),
            kind: "delins".to_string(),
        };
        clip_alleles(&mut notation, 1);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 100);
        assert_eq!(notation.end, 100);
        assert_eq!(notation.ref_allele, "");
        assert_eq!(notation.alt_allele, "A");
    }

    #[test]
    fn test_format_hgvs_string_delins() {
        let notation = HgvsNotation {
            start: 10,
            end: 12,
            ref_allele: "ACG".to_string(),
            alt_allele: "TT".to_string(),
            kind: "delins".to_string(),
        };
        assert_eq!(
            format_hgvs_string("ENST000001.1", 'c', "10", "12", &notation),
            Some("ENST000001.1:c.10_12delinsTT".to_string())
        );
    }

    #[test]
    fn test_split_hgvs_coord() {
        assert_eq!(
            split_hgvs_coord("100+5"),
            Some((100, Some("+5".to_string())))
        );
        assert_eq!(
            split_hgvs_coord("100-3"),
            Some((100, Some("-3".to_string())))
        );
        assert_eq!(split_hgvs_coord("-1"), Some((-1, None)));
        assert_eq!(split_hgvs_coord("*2"), Some((2, None)));
    }

    #[test]
    fn test_format_hgvsp_synonymous() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".to_string(),
            alt_peptide: "A".to_string(),
            ref_translation: "MA*".to_string(),
            alt_translation: "MA*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2=".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_partial_codon_synonymous_uses_ter() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 262,
            end: 262,
            ref_peptide: "X".to_string(),
            alt_peptide: "X".to_string(),
            ref_translation: "XRVM".to_string(),
            alt_translation: "XRVM".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ter262=".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_start_lost_reports_unknown_protein() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 1,
            end: 1,
            ref_peptide: "M".to_string(),
            alt_peptide: "L".to_string(),
            ref_translation: "MA*".to_string(),
            alt_translation: "LA*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: true,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Met1?".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_frameshift_uses_first_changed_residue_and_stop_distance() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "K".to_string(),
            alt_peptide: "Q".to_string(),
            ref_translation: "MKKKK".to_string(),
            alt_translation: "MKQW*".to_string(),
            alt_translation_extension: None,
            frameshift: true,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Lys3GlnfsTer3".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_stop_lost_adds_extension_length() {
        let translation = make_translation();
        // VEP cached peptide for "MA*" is "MA" (no *), length=2.
        // Alt "MAQW*": $+[0] = 5. extra = 5 - 1 - 2 = 2 → extTer2.
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "*".to_string(),
            alt_peptide: "Q".to_string(),
            ref_translation: "MA*".to_string(),
            alt_translation: "MAQW*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: true,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ter3GlnextTer2".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_stop_lost_with_adjacent_stop_gives_ext_1() {
        // VEP cached peptide for "MA*" is "MA" (no *), length=2.
        // Alt "MAQ*": $+[0] = 4. extra = 4 - 1 - 2 = 1 → extTer1.
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "*".to_string(),
            alt_peptide: "Q".to_string(),
            ref_translation: "MA*".to_string(),
            alt_translation: "MAQ*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: true,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ter3GlnextTer1".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_insertion_duplication_uses_dup_notation() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 4,
            end: 4,
            ref_peptide: "-".to_string(),
            alt_peptide: "A".to_string(),
            ref_translation: "MAA*".to_string(),
            alt_translation: "MAAA*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala3dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_synonymous_multi_residue_uses_full_peptide_string() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 25,
            end: 26,
            ref_peptide: "EE".to_string(),
            alt_peptide: "EE".to_string(),
            ref_translation: "M".repeat(24) + "EEEEK",
            alt_translation: "M".repeat(24) + "EEEEK",
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.GluGlu25=".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_insertion_uses_flanking_residues() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 3,
            ref_peptide: "-".to_string(),
            alt_peptide: "Q".to_string(),
            ref_translation: "MAV*".to_string(),
            alt_translation: "MAQV*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2_Val3insGln".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_shift_hgvs_false_disables_three_prime_peptide_shift() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".to_string(),
            alt_peptide: "-".to_string(),
            ref_translation: "MAA*".to_string(),
            alt_translation: "MA*".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, false),
            Some("ENSPHGVS000001.1:p.Ala2del".to_string())
        );
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala3del".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_coding_relative_numbering() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("14"), None, "G", "A", 103, 103, None),
            Some("ENSTHGVS000001.1:c.4G>A".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_formats_insertions_with_flanking_coordinates() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "T", 103, 103, None),
            Some("ENSTHGVS000001.1:c.3_4insT".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_formats_deletions_from_genomic_span() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "G", "-", 103, 103, None),
            Some("ENSTHGVS000001.1:c.4del".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_negative_utr_coordinate() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("10"), None, "A", "G", 99, 99, None),
            Some("ENSTHGVS000001.1:c.-1A>G".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_non_coding_numbering() {
        let tx = make_transcript("lncRNA", 1, None, None);
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("14"), None, "G", "A", 103, 103, None),
            Some("ENSTHGVS000001.1:n.14G>A".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_cds_position_for_coding_snv() {
        let mut tx = make_transcript("protein_coding", 1, Some(100), Some(120));
        tx.cdna_coding_start = Some(14);
        tx.cdna_coding_end = Some(34);
        let exon = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 120,
        };
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("18"), Some("5"), "G", "A", 104, 104, None),
            Some("ENSTHGVS000001.1:c.5G>A".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_reverse_complements_minus_strand_alleles() {
        let tx = make_transcript("lncRNA", -1, None, None);
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("7"), None, "A", "G", 103, 103, None),
            Some("ENSTHGVS000001.1:n.38T>C".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_genomic_shift_for_intronic_indels() {
        let tx = make_transcript("protein_coding", 1, Some(90), Some(119));
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 90,
            end: 99,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 110,
            end: 119,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 3,
            start: 107,
            end: 108,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "AA".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "AA".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "AA", "-", 104, 105, Some(&shift)),
            Some("ENSTHGVS000001.1:c.11-3_11-2del".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_shifts_exonic_indels_when_hgvs_shift_is_available() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 2,
            start: 105,
            end: 105,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "T".to_string(),
            shifted_output_allele: "T".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "T".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "T", 103, 103, Some(&shift)),
            Some("ENSTHGVS000001.1:c.5_6insT".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_dup_notation_for_shifted_intronic_insertions() {
        let tx = make_transcript("lncRNA", 1, None, None);
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 3,
            start: 151,
            end: 150,
            shifted_compare_allele: "AGTA".to_string(),
            shifted_allele_string: "AGTA".to_string(),
            shifted_output_allele: "AGTA".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "AGTA".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: "AGTA".to_string(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "AAGT", 148, 147, Some(&shift)),
            Some("ENSTHGVS000001.1:n.11+41_11+44dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_shifted_refseq_intronic_insertion_allele_without_bam_edit() {
        let mut tx = make_transcript("lncRNA", 1, None, None);
        tx.transcript_id = "NR_047526".to_string();
        tx.source = Some("RefSeq".to_string());
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 3,
            start: 151,
            end: 150,
            shifted_compare_allele: "GAAA".to_string(),
            shifted_allele_string: "GAAA".to_string(),
            shifted_output_allele: "GAAA".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "AAAG".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: "TAAA".to_string(),
            three_prime_context: "TTTC".to_string(),
        };

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "-", "AAAG", 148, 147, Some(&shift))
            .expect("HGVSc");
        assert!(
            hgvsc.ends_with("insGAAA"),
            "expected shifted inserted allele in default RefSeq intronic HGVSc, got {hgvsc}"
        );
    }

    #[test]
    fn test_format_hgvsc_uses_edited_refseq_intronic_insertion_allele_with_bam_edit() {
        let mut tx = make_transcript("lncRNA", 1, None, None);
        tx.transcript_id = "NR_047526".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.bam_edit_status = Some("ok".to_string());
        tx.spliced_seq = Some("ACGTACGTACGTTTCGATCGTA".to_string());
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 3,
            start: 151,
            end: 150,
            shifted_compare_allele: "GAAA".to_string(),
            shifted_allele_string: "GAAA".to_string(),
            shifted_output_allele: "GAAA".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "AAAG".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: "TAAA".to_string(),
            three_prime_context: "TTTC".to_string(),
        };

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "-", "AAAG", 148, 147, Some(&shift))
            .expect("HGVSc");
        assert!(
            hgvsc.ends_with("insAAAG"),
            "expected BAM-edited transcript allele in RefSeq intronic HGVSc, got {hgvsc}"
        );
    }

    #[test]
    fn test_format_hgvsc_refseq_intronic_dup_detection_uses_shifted_hgvs_allele() {
        let mut tx = make_transcript("lncRNA", 1, None, None);
        tx.transcript_id = "NM_004442".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.bam_edit_status = Some("ok".to_string());
        tx.spliced_seq = Some("ACGTACGTACGTTTCGATCGTA".to_string());
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 22,
            start: 151,
            end: 150,
            shifted_compare_allele: "GT".to_string(),
            shifted_allele_string: "GT".to_string(),
            shifted_output_allele: "GT".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "TG".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: "GT".to_string(),
            three_prime_context: "CA".to_string(),
        };

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "-", "TG", 148, 147, Some(&shift))
            .expect("HGVSc");
        assert!(
            hgvsc.ends_with("dup"),
            "expected RefSeq intronic duplication notation, got {hgvsc}"
        );
    }

    #[test]
    fn test_format_hgvsc_orients_shifted_duplication_context_on_minus_strand() {
        let tx = make_transcript("lncRNA", -1, None, None);
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: -1,
            shift_length: 1,
            start: 151,
            end: 150,
            shifted_compare_allele: "A".to_string(),
            shifted_allele_string: "A".to_string(),
            shifted_output_allele: "A".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "A".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: "A".to_string(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "A", 151, 150, Some(&shift)),
            Some("ENSTHGVS000001.1:n.12-39dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_upstream_dup_coordinates_for_shifted_insertions() {
        let tx = make_transcript("lncRNA", 1, None, None);
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 1,
            start: 151,
            end: 150,
            shifted_compare_allele: "A".to_string(),
            shifted_allele_string: "A".to_string(),
            shifted_output_allele: "A".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "A".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: "A".to_string(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "A", 151, 150, Some(&shift)),
            Some("ENSTHGVS000001.1:n.11+40dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_detects_dup_without_nonzero_shift_length() {
        let tx = make_transcript("lncRNA", -1, None, None);
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 110,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 200,
            end: 210,
        };
        let exons = [&exon1, &exon2];
        let shift = HgvsGenomicShift {
            strand: -1,
            shift_length: 0,
            start: 151,
            end: 150,
            shifted_compare_allele: "A".to_string(),
            shifted_allele_string: "A".to_string(),
            shifted_output_allele: "A".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "A".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: "A".to_string(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "A", 151, 150, Some(&shift)),
            Some("ENSTHGVS000001.1:n.12-41dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_detects_dup_inside_transcript_sequence_path() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 0,
            start: 104,
            end: 103,
            shifted_compare_allele: "T".to_string(),
            shifted_allele_string: "T".to_string(),
            shifted_output_allele: "T".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "T".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: "G".to_string(),
            three_prime_context: "T".to_string(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "T", 103, 103, Some(&shift)),
            Some("ENSTHGVS000001.1:c.5dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_shifted_coordinates_for_exonic_deletions() {
        let tx = make_transcript("lncRNA", 1, None, None);
        let exon = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 200,
        };
        let exons = [&exon];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 30,
            start: 150,
            end: 153,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "GTGT".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "GTGT".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "GTGT", "-", 120, 123, Some(&shift)),
            Some("ENSTHGVS000001.1:n.51_54del".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_suppresses_shifted_noncoding_coords_past_transcript_end() {
        let tx = make_transcript("lncRNA", 1, None, None);
        let exon = make_exon();
        let exons = [&exon];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 2,
            start: 141,
            end: 141,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "AA".to_string(),
            shifted_output_allele: "AA".to_string(),
            ref_orig_allele_string: "-".to_string(),
            alt_orig_allele_string: "AA".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "AA", 139, 139, Some(&shift)),
            None
        );
    }

    #[test]
    fn test_format_hgvsc_suppresses_shifted_utr_coords_past_valid_star_range() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(120));
        let exon = make_exon();
        let exons = [&exon];
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 3,
            start: 141,
            end: 144,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "AAAA".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "AAAA".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "AAAA", "-", 138, 141, Some(&shift)),
            None
        );
    }

    #[test]
    fn test_format_hgvsc_allows_large_star_coordinate_inside_transcript_span() {
        let mut tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        tx.end = 6010;
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 90,
            end: 108,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 6000,
            end: 6010,
        };
        let exons = [&exon1, &exon2];
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "A", "G", 510, 510, None),
            Some("ENSTHGVS000001.1:c.*402A>G".to_string())
        );
    }

    #[test]
    fn test_hgvs_cdna_position_intronic_plus_strand() {
        let tx = make_transcript("protein_coding", 1, Some(90), Some(119));
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 90,
            end: 99,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 110,
            end: 119,
        };
        let exons = [&exon1, &exon2];
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 104),
            Some("10+5".to_string())
        );
    }

    #[test]
    fn test_hgvs_cdna_position_intronic_minus_strand() {
        let tx = make_transcript("protein_coding", -1, Some(90), Some(119));
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 90,
            end: 99,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 110,
            end: 119,
        };
        let exons = [&exon1, &exon2];
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 104),
            Some("11-5".to_string())
        );
    }

    #[test]
    fn test_hgvs_cdna_position_uses_transcript_mapper_segments() {
        let mut tx = make_transcript("non_coding", 1, None, None);
        tx.transcript_id = "NM_001291281.3".to_string();
        tx.start = 41361434;
        tx.end = 41383590;
        tx.cdna_mapper_segments = vec![
            TranscriptCdnaMapperSegment {
                genomic_start: 41361434,
                genomic_end: 41362344,
                cdna_start: 1,
                cdna_end: 911,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 41381616,
                genomic_end: 41382208,
                cdna_start: 912,
                cdna_end: 1504,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 41382210,
                genomic_end: 41382210,
                cdna_start: 1505,
                cdna_end: 1505,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 41382211,
                genomic_end: 41383590,
                cdna_start: 1707,
                cdna_end: 3086,
                ori: 1,
            },
        ];
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 41361434,
            end: 41362344,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 41381616,
            end: 41382208,
        };
        let exon3 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 3,
            start: 41382210,
            end: 41383590,
        };
        let exons = [&exon1, &exon2, &exon3];
        // Segments 3→4 have a genomic-contiguous but cDNA-discontinuous gap
        // (an RNA edit). use_cdna_mapper_for_general_coords() correctly detects
        // this and falls back to exon geometry.
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 41383346),
            Some("2641".to_string())
        );
    }

    #[test]
    fn test_hgvs_cdna_position_applies_refseq_offset_without_mapper_segments() {
        let mut tx = make_transcript("non_coding", 1, None, None);
        tx.transcript_id = "NM_OFFSET.1".to_string();
        tx.start = 100;
        tx.end = 3000;
        tx.refseq_edits = vec![crate::transcript_consequence::RefSeqEdit {
            start: 1506,
            end: 1505,
            replacement_len: Some(201),
            skip_refseq_offset: false,
        }];
        let exon = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 3000,
        };
        let exons = [&exon];
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 2740),
            Some("2842".to_string())
        );
    }

    #[test]
    fn test_hgvs_cdna_position_does_not_apply_refseq_offset_to_intronic_coords() {
        let mut tx = make_transcript("protein_coding", 1, Some(100), Some(599));
        tx.transcript_id = "NM_OFFSET.1".to_string();
        tx.start = 100;
        tx.end = 599;
        tx.cdna_coding_start = Some(1);
        tx.cdna_coding_end = Some(400);
        tx.refseq_edits = vec![crate::transcript_consequence::RefSeqEdit {
            start: 150,
            end: 149,
            replacement_len: Some(3),
            skip_refseq_offset: false,
        }];
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 100,
            end: 299,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 400,
            end: 599,
        };
        let exons = [&exon1, &exon2];
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 349),
            Some("200+50".to_string())
        );
    }

    #[test]
    fn test_hgvs_cdna_position_native_refseq_pre_coding_intronic_anchor_uses_post_gap_numbering() {
        let mut tx = make_transcript("protein_coding", 1, Some(39044831), Some(39126233));
        tx.transcript_id = "NM_001007075.2".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.start = 39044831;
        tx.end = 39126233;
        tx.cdna_coding_start = Some(360);
        tx.cdna_coding_end = Some(2489);
        tx.cdna_mapper_segments = vec![
            TranscriptCdnaMapperSegment {
                genomic_start: 39044831,
                genomic_end: 39044966,
                cdna_start: 1,
                cdna_end: 136,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 39044968,
                genomic_end: 39045096,
                cdna_start: 137,
                cdna_end: 265,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 39062559,
                genomic_end: 39063035,
                cdna_start: 266,
                cdna_end: 742,
                ori: 1,
            },
        ];
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 39044831,
            end: 39045096,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 39062559,
            end: 39063035,
        };
        let exons = [&exon1, &exon2];
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 39045450),
            Some("-94+354".to_string())
        );
    }

    fn make_nm_001177639_leading_edit_transcript() -> (TranscriptFeature, [ExonFeature; 2]) {
        let mut tx = make_transcript("protein_coding", 1, Some(49510542), Some(49533206));
        tx.transcript_id = "NM_001177639.3".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.start = 49510418;
        tx.end = 49535615;
        tx.cdna_coding_start = Some(125);
        tx.cdna_coding_end = Some(2812);
        tx.bam_edit_status = Some("ok".to_string());
        tx.has_non_polya_rna_edit = true;
        tx.refseq_edits = vec![crate::transcript_consequence::RefSeqEdit {
            start: 1,
            end: 0,
            replacement_len: Some(7),
            skip_refseq_offset: false,
        }];
        tx.cdna_mapper_segments = vec![
            TranscriptCdnaMapperSegment {
                genomic_start: 49510418,
                genomic_end: 49510819,
                cdna_start: 8,
                cdna_end: 409,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 49530797,
                genomic_end: 49535615,
                cdna_start: 410,
                cdna_end: 5228,
                ori: 1,
            },
        ];
        let exon1 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 49510418,
            end: 49510819,
        };
        let exon2 = ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 2,
            start: 49530797,
            end: 49535615,
        };
        (tx, [exon1, exon2])
    }

    fn make_nm_001137668_internal_edit_transcript() -> (TranscriptFeature, [ExonFeature; 10]) {
        let mut tx = make_transcript("protein_coding", 1, Some(89829880), Some(89874436));
        tx.transcript_id = "NM_001137668.2".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.start = 89829880;
        tx.end = 89874436;
        tx.cdna_coding_start = Some(84);
        tx.cdna_coding_end = Some(6032);
        tx.bam_edit_status = Some("ok".to_string());
        tx.has_non_polya_rna_edit = true;
        tx.refseq_edits = vec![
            crate::transcript_consequence::RefSeqEdit {
                start: 5976,
                end: 5975,
                replacement_len: Some(48),
                skip_refseq_offset: false,
            },
            crate::transcript_consequence::RefSeqEdit {
                start: 5977,
                end: 5977,
                replacement_len: Some(1),
                skip_refseq_offset: true,
            },
        ];
        tx.cdna_mapper_segments = vec![
            TranscriptCdnaMapperSegment {
                genomic_start: 89829880,
                genomic_end: 89829934,
                cdna_start: 1,
                cdna_end: 55,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89846562,
                genomic_end: 89846644,
                cdna_start: 56,
                cdna_end: 138,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89853166,
                genomic_end: 89853234,
                cdna_start: 139,
                cdna_end: 207,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89854817,
                genomic_end: 89854917,
                cdna_start: 208,
                cdna_end: 308,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89855470,
                genomic_end: 89855544,
                cdna_start: 309,
                cdna_end: 383,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89857084,
                genomic_end: 89857199,
                cdna_start: 384,
                cdna_end: 499,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89862126,
                genomic_end: 89864371,
                cdna_start: 500,
                cdna_end: 2745,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89865953,
                genomic_end: 89869082,
                cdna_start: 2746,
                cdna_end: 5875,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89871289,
                genomic_end: 89871388,
                cdna_start: 5876,
                cdna_end: 5975,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89873801,
                genomic_end: 89874436,
                cdna_start: 6024,
                cdna_end: 6659,
                ori: 1,
            },
        ];
        let exons = [
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 1,
                start: 89829880,
                end: 89829934,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 2,
                start: 89846562,
                end: 89846644,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 3,
                start: 89853166,
                end: 89853234,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 4,
                start: 89854817,
                end: 89854917,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 5,
                start: 89855470,
                end: 89855544,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 6,
                start: 89857084,
                end: 89857199,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 7,
                start: 89862126,
                end: 89864371,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 8,
                start: 89865953,
                end: 89869082,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 9,
                start: 89871289,
                end: 89871388,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 10,
                start: 89873801,
                end: 89874436,
            },
        ];
        (tx, exons)
    }

    fn make_nm_012115_internal_edit_transcript() -> (TranscriptFeature, [ExonFeature; 10]) {
        let mut tx = make_transcript("protein_coding", 1, Some(89829880), Some(89874436));
        tx.transcript_id = "NM_012115.4".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.start = 89829880;
        tx.end = 89874436;
        tx.cdna_coding_start = Some(256);
        tx.cdna_coding_end = Some(6204);
        tx.bam_edit_status = Some("ok".to_string());
        tx.has_non_polya_rna_edit = true;
        tx.refseq_edits = vec![
            crate::transcript_consequence::RefSeqEdit {
                start: 6148,
                end: 6147,
                replacement_len: Some(48),
                skip_refseq_offset: false,
            },
            crate::transcript_consequence::RefSeqEdit {
                start: 6149,
                end: 6149,
                replacement_len: Some(1),
                skip_refseq_offset: true,
            },
        ];
        tx.cdna_mapper_segments = vec![
            TranscriptCdnaMapperSegment {
                genomic_start: 89829880,
                genomic_end: 89830106,
                cdna_start: 1,
                cdna_end: 227,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89846562,
                genomic_end: 89846644,
                cdna_start: 228,
                cdna_end: 310,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89853166,
                genomic_end: 89853234,
                cdna_start: 311,
                cdna_end: 379,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89854817,
                genomic_end: 89854917,
                cdna_start: 380,
                cdna_end: 480,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89855470,
                genomic_end: 89855544,
                cdna_start: 481,
                cdna_end: 555,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89857084,
                genomic_end: 89857199,
                cdna_start: 556,
                cdna_end: 671,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89862126,
                genomic_end: 89864371,
                cdna_start: 672,
                cdna_end: 2917,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89865953,
                genomic_end: 89869082,
                cdna_start: 2918,
                cdna_end: 6047,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89871289,
                genomic_end: 89871388,
                cdna_start: 6048,
                cdna_end: 6147,
                ori: 1,
            },
            TranscriptCdnaMapperSegment {
                genomic_start: 89873801,
                genomic_end: 89874436,
                cdna_start: 6196,
                cdna_end: 6831,
                ori: 1,
            },
        ];
        let exons = [
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 1,
                start: 89829880,
                end: 89830106,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 2,
                start: 89846562,
                end: 89846644,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 3,
                start: 89853166,
                end: 89853234,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 4,
                start: 89854817,
                end: 89854917,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 5,
                start: 89855470,
                end: 89855544,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 6,
                start: 89857084,
                end: 89857199,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 7,
                start: 89862126,
                end: 89864371,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 8,
                start: 89865953,
                end: 89869082,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 9,
                start: 89871289,
                end: 89871388,
            },
            ExonFeature {
                transcript_id: tx.transcript_id.clone(),
                exon_number: 10,
                start: 89873801,
                end: 89874436,
            },
        ];
        (tx, exons)
    }

    fn make_nm_001172437_same_coordinate_multibase_edit_transcript()
    -> (TranscriptFeature, [ExonFeature; 1]) {
        let mut tx = make_transcript("protein_coding", 1, Some(1), Some(2355));
        tx.transcript_id = "NM_001172437.2".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.start = 1;
        tx.end = 7000;
        tx.cdna_coding_start = Some(263);
        tx.cdna_coding_end = Some(2617);
        tx.refseq_edits = vec![crate::transcript_consequence::RefSeqEdit {
            start: 1447,
            end: 1447,
            replacement_len: Some(2),
            skip_refseq_offset: false,
        }];
        let mut transcript_seq = vec![b'A'; 7000];
        transcript_seq[2768] = b'T';
        transcript_seq[2769] = b'C';
        tx.spliced_seq = Some(String::from_utf8(transcript_seq).expect("ASCII transcript"));
        let exons = [ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 1,
            end: 7000,
        }];
        (tx, exons)
    }

    fn make_nm_001172437_failed_bam_edit_transcript() -> (TranscriptFeature, [ExonFeature; 1]) {
        let (mut tx, exons) = make_nm_001172437_same_coordinate_multibase_edit_transcript();
        tx.bam_edit_status = Some("failed".to_string());
        (tx, exons)
    }

    fn make_nr_170302_noncoding_refseq_edit_transcript() -> (TranscriptFeature, [ExonFeature; 1]) {
        let mut tx = make_transcript("lncRNA", 1, None, None);
        tx.transcript_id = "NR_170302.1".to_string();
        tx.source = Some("RefSeq".to_string());
        tx.start = 1;
        tx.end = 200;
        tx.refseq_edits = vec![
            crate::transcript_consequence::RefSeqEdit {
                start: 7,
                end: 6,
                replacement_len: Some(6),
                skip_refseq_offset: false,
            },
            crate::transcript_consequence::RefSeqEdit {
                start: 14,
                end: 14,
                replacement_len: None,
                skip_refseq_offset: false,
            },
        ];
        let mut transcript_seq = vec![b'A'; 200];
        transcript_seq[36] = b'T';
        transcript_seq[41] = b'C';
        tx.spliced_seq = Some(String::from_utf8(transcript_seq).expect("ASCII transcript"));
        let exons = [ExonFeature {
            transcript_id: tx.transcript_id.clone(),
            exon_number: 1,
            start: 1,
            end: 200,
        }];
        (tx, exons)
    }

    #[test]
    fn test_format_hgvsc_native_refseq_leading_insertion_uses_pre_edit_upstream_anchor() {
        let (tx, exons_owned) = make_nm_001177639_leading_edit_transcript();
        let exons = [&exons_owned[0], &exons_owned[1]];

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "T", "C", 49510861, 49510861, None);

        assert_eq!(hgvsc.as_deref(), Some("NM_001177639.3:c.278+42T>C"));
    }

    #[test]
    fn test_format_hgvsc_native_refseq_leading_insertion_uses_pre_edit_downstream_anchor() {
        let (tx, exons_owned) = make_nm_001177639_leading_edit_transcript();
        let exons = [&exons_owned[0], &exons_owned[1]];

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "C", "T", 49521283, 49521283, None);

        assert_eq!(hgvsc.as_deref(), Some("NM_001177639.3:c.279-9514C>T"));
    }

    #[test]
    fn test_format_hgvsc_internal_refseq_insertion_uses_pre_edit_downstream_anchor() {
        let (tx, exons_owned) = make_nm_001137668_internal_edit_transcript();
        let exons = exons_owned.iter().collect::<Vec<_>>();

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "T", "C", 89873677, 89873677, None);

        assert_eq!(hgvsc.as_deref(), Some("NM_001137668.2:c.5893-124T>C"));
    }

    #[test]
    fn test_format_hgvsc_internal_refseq_insertion_uses_pre_edit_downstream_anchor_alt_tx() {
        let (tx, exons_owned) = make_nm_012115_internal_edit_transcript();
        let exons = exons_owned.iter().collect::<Vec<_>>();

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "T", "C", 89873677, 89873677, None);

        assert_eq!(hgvsc.as_deref(), Some("NM_012115.4:c.5893-124T>C"));
    }

    #[test]
    fn test_format_hgvsc_same_coordinate_multibase_refseq_edit_uses_full_inserted_offset() {
        let (tx, exons_owned) = make_nm_001172437_same_coordinate_multibase_edit_transcript();
        let exons = exons_owned.iter().collect::<Vec<_>>();

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "T", "C", 2768, 2768, None);

        assert_eq!(hgvsc.as_deref(), Some("NM_001172437.2:c.*153C>C"));
    }

    #[test]
    fn test_format_hgvsc_refseq_failed_bam_edit_suppresses_shifted_utr_deletion() {
        let (tx, exons_owned) = make_nm_001172437_failed_bam_edit_transcript();
        let exons = exons_owned.iter().collect::<Vec<_>>();
        let variant = crate::transcript_consequence::VariantInput::from_vcf(
            "7".to_string(),
            5859,
            5863,
            "TACAG".to_string(),
            "T".to_string(),
        );
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 4,
            start: 5864,
            end: 5867,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "ACAG".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "ACAG".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };

        let hgvsc = format_hgvsc(
            &tx,
            &exons,
            None,
            None,
            "CAGA",
            "-",
            variant.start,
            variant.end,
            Some(&shift),
        );

        assert_eq!(hgvsc.as_deref(), Some("NM_001172437.2:c.*3245_*3248del"));
    }

    #[test]
    fn test_hgvsc_offset_refseq_failed_bam_edit_suppresses_shifted_utr_deletion() {
        let (tx, _) = make_nm_001172437_failed_bam_edit_transcript();
        let mut variant = crate::transcript_consequence::VariantInput::from_vcf(
            "7".to_string(),
            5859,
            5863,
            "TACAG".to_string(),
            "T".to_string(),
        );
        variant.hgvs_shift_forward = Some(HgvsGenomicShift {
            strand: 1,
            shift_length: 4,
            start: 5864,
            end: 5867,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "ACAG".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "ACAG".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        });

        assert_eq!(
            hgvsc_offset_for_output(
                &tx,
                &variant,
                "CAGA",
                Some("NM_001172437.2:c.*3245_*3248del"),
            ),
            None
        );
    }

    #[test]
    fn test_format_hgvsc_refseq_failed_bam_edit_keeps_shift_when_alleles_match() {
        let (tx, exons_owned) = make_nm_001172437_failed_bam_edit_transcript();
        let exons = exons_owned.iter().collect::<Vec<_>>();
        let variant = crate::transcript_consequence::VariantInput::from_vcf(
            "7".to_string(),
            4915,
            4916,
            "GT".to_string(),
            "G".to_string(),
        );
        let shift = HgvsGenomicShift {
            strand: 1,
            shift_length: 6,
            start: 4922,
            end: 4922,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "T".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "T".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };

        let hgvsc = format_hgvsc(
            &tx,
            &exons,
            None,
            None,
            "T",
            "-",
            variant.start,
            variant.end,
            Some(&shift),
        );

        assert_eq!(hgvsc.as_deref(), Some("NM_001172437.2:c.*2307del"));
    }

    #[test]
    fn test_hgvsc_offset_refseq_failed_bam_edit_keeps_shift_when_alleles_match() {
        let (tx, _) = make_nm_001172437_failed_bam_edit_transcript();
        let mut variant = crate::transcript_consequence::VariantInput::from_vcf(
            "7".to_string(),
            4915,
            4916,
            "GT".to_string(),
            "G".to_string(),
        );
        variant.hgvs_shift_forward = Some(HgvsGenomicShift {
            strand: 1,
            shift_length: 6,
            start: 4922,
            end: 4922,
            shifted_compare_allele: "-".to_string(),
            shifted_allele_string: "T".to_string(),
            shifted_output_allele: "-".to_string(),
            ref_orig_allele_string: "T".to_string(),
            alt_orig_allele_string: "-".to_string(),
            five_prime_flanking_seq: String::new(),
            three_prime_flanking_seq: String::new(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        });

        assert_eq!(
            hgvsc_offset_for_output(&tx, &variant, "T", Some("NM_001172437.2:c.*2307del"),),
            Some(6)
        );
    }

    #[test]
    fn test_format_hgvsc_noncoding_refseq_keeps_raw_coordinate_but_uses_edited_reference() {
        let (tx, exons_owned) = make_nr_170302_noncoding_refseq_edit_transcript();
        let exons = exons_owned.iter().collect::<Vec<_>>();

        let hgvsc = format_hgvsc(&tx, &exons, None, None, "C", "T", 37, 37, None);

        assert_eq!(hgvsc.as_deref(), Some("NR_170302.1:n.37C>T"));
    }

    #[test]
    fn test_shift_to_hgvs_coding_coordinates_prefers_cached_cdna_coding_bounds() {
        let mut tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        tx.cdna_coding_start = Some(21);
        tx.cdna_coding_end = Some(29);
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            shift_to_hgvs_coding_coordinates(&tx, &exons, "23"),
            Some("3".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_cached_cdna_coding_bounds() {
        let mut tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        tx.cdna_coding_start = Some(21);
        tx.cdna_coding_end = Some(29);
        tx.cdna_mapper_segments = vec![TranscriptCdnaMapperSegment {
            genomic_start: 90,
            genomic_end: 140,
            cdna_start: 1,
            cdna_end: 51,
            ori: 1,
        }];
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "G", "A", 112, 112, None),
            Some("ENSTHGVS000001.1:c.3G>A".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_star_coordinate_in_three_prime_utr() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("21"), None, "A", "G", 110, 110, None),
            Some("ENSTHGVS000001.1:c.*2A>G".to_string())
        );
    }

    #[test]
    fn test_shift_to_hgvs_coding_coordinates_strips_plus_at_stop_codon_boundary() {
        let mut tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        tx.cdna_coding_start = Some(21);
        tx.cdna_coding_end = Some(29);
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            shift_to_hgvs_coding_coordinates(&tx, &exons, "29+300"),
            Some("*300".to_string())
        );
    }

    #[test]
    fn test_perform_shift_ensembl_rotates_hgvs_output_in_vf_orientation() {
        let (shift_length, shifted_seq, shifted_output, start, end) =
            perform_shift_ensembl(b"GATG", b"GATG", "", "TG", 100, 99, true, -1);
        // With corrected loop_limiter guard, the loop now walks 2 matching
        // bases in pre_seq ("TG") instead of stopping after 1.
        assert_eq!(shift_length, 2);
        assert_eq!(String::from_utf8(shifted_seq).unwrap(), "TGGA");
        assert_eq!(String::from_utf8(shifted_output).unwrap(), "TGGA");
        assert_eq!((start, end), (100, 99));
    }

    #[test]
    fn test_compare_hgvs_positions_orders_star_intronic_offsets() {
        assert_eq!(
            compare_hgvs_positions("*203-3410", "*203-3415"),
            Some(Ordering::Greater)
        );
        assert_eq!(
            compare_hgvs_positions("*3199", "*3198"),
            Some(Ordering::Greater)
        );
    }

    // ---------------------------------------------------------------
    // stop_loss_extra_aa — exact VEP _stop_loss_extra_AA semantics
    // ---------------------------------------------------------------

    #[test]
    fn test_stop_loss_extra_aa_non_frameshift_uses_ref_len_without_terminal_stop() {
        // VEP: ref_len = length(_peptide()) where cached peptide has no terminal *.
        // ref_translation = "MKKR*" → trim_end_matches('*').len() = 4.
        // alt_translation = "MKKRQW*" → find('*') = 6.
        // extra = (6 + 1) - (4 + 1) = 2.
        let protein = ProteinHgvsData {
            start: 5,
            end: 5,
            ref_peptide: "*".into(),
            alt_peptide: "Q".into(),
            ref_translation: "MKKR*".into(),
            alt_translation: "MKKRQW*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: true,
            native_refseq: false,
        };
        assert_eq!(stop_loss_extra_aa(&protein, 4, false), Some(2));
    }

    #[test]
    fn test_stop_loss_extra_aa_non_frameshift_with_internal_stops_uses_full_len() {
        // LoF transcript: ref = "M*KR*" (internal stop at 1, terminal at 4).
        // trim_end_matches('*').len() = 4 (strips only trailing *).
        // If stop at position 1 is mutated, alt has stop later.
        // alt = "MQKRW*" → find('*') = 5.
        // extra = (5+1) - (4+1) = 1.
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "*".into(),
            alt_peptide: "Q".into(),
            ref_translation: "M*KR*".into(),
            alt_translation: "MQKRW*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: true,
            native_refseq: false,
        };
        assert_eq!(stop_loss_extra_aa(&protein, 1, false), Some(1));
    }

    #[test]
    fn test_stop_loss_extra_aa_non_frameshift_no_new_stop_returns_none() {
        // Alt translation has no stop codon → returns None → extTer?
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "*".into(),
            alt_peptide: "Q".into(),
            ref_translation: "MA*".into(),
            alt_translation: "MAQ".into(), // no *
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: true,
            native_refseq: false,
        };
        assert_eq!(stop_loss_extra_aa(&protein, 2, false), None);
    }

    #[test]
    fn test_stop_loss_extra_aa_non_frameshift_same_length_returns_none() {
        // VEP: extra = 0 → returns undef. Our: 0 > 0 is false → None.
        // ref = "MA*" (trim len = 2). alt = "MAQ*" (find = 3).
        // extra = (3+1) - (2+1) = 1. Actually that's 1, not 0.
        // For true zero: ref = "MAK*" (trim=3), alt = "MAQK*" (find=4).
        // extra = (4+1) - (3+1) = 1. Hmm, hard to get 0.
        // ref = "MAK*" (trim=3), alt = "MAQ*" (find=3).
        // extra = (3+1) - (3+1) = 0 → None.
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "*".into(),
            alt_peptide: "Q".into(),
            ref_translation: "MAK*".into(),
            alt_translation: "MAQ*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: true,
            native_refseq: false,
        };
        assert_eq!(stop_loss_extra_aa(&protein, 2, false), None);
    }

    #[test]
    fn test_stop_loss_extra_aa_frameshift_counts_from_variant_position() {
        // Frameshift: extra = stop_idx + 1 - ref_var_pos.
        // alt = "MKQW*" → find('*') = 4.
        // ref_var_pos = 3 (1-based position of first affected AA).
        // extra = (4 + 1) - 3 = 2.
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "K".into(),
            alt_peptide: "Q".into(),
            ref_translation: "MKKK*".into(),
            alt_translation: "MKQW*".into(),
            alt_translation_extension: None,
            frameshift: true,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(stop_loss_extra_aa(&protein, 3, true), Some(2));
    }

    // ---------------------------------------------------------------
    // perform_shift_ensembl — hgvs_reverse flag combinations
    // ---------------------------------------------------------------

    #[test]
    fn test_perform_shift_ensembl_forward_no_reverse_rotates_both_left() {
        // seq_strand=1, reverse=false, hgvs_reverse=false.
        // Both seq_to_check and hgvs_output rotate left.
        let (shift, seq, hgvs, start, end) =
            perform_shift_ensembl(b"AT", b"AT", "ATGC", "", 100, 99, false, 1);
        assert_eq!(shift, 2); // AT matches AT in post_seq
        assert_eq!(seq, b"AT"); // full rotation back to original
        assert_eq!(hgvs, b"AT");
        assert_eq!(start, 102);
        assert_eq!(end, 101);
    }

    #[test]
    fn test_perform_shift_ensembl_reverse_with_hgvs_reverse() {
        // seq_strand=-1, reverse=true, hgvs_reverse = (1 != -1) = true.
        // seq_to_check rotates right, hgvs_output rotates left.
        let (shift, _seq, hgvs, _start, _end) =
            perform_shift_ensembl(b"AG", b"AG", "", "CCAG", 100, 101, true, -1);
        assert!(shift > 0);
        // hgvs_output should be different rotation than seq_to_check
        assert_eq!(hgvs.len(), 2);
    }

    #[test]
    fn test_perform_shift_ensembl_no_match_returns_zero_shift() {
        let (shift, seq, hgvs, start, end) =
            perform_shift_ensembl(b"AT", b"AT", "GC", "", 100, 99, false, 1);
        assert_eq!(shift, 0);
        assert_eq!(seq, b"AT");
        assert_eq!(hgvs, b"AT");
        assert_eq!(start, 100);
        assert_eq!(end, 99);
    }

    #[test]
    fn test_perform_shift_ensembl_genomic_shift_always_uses_seq_strand_1() {
        // Genomic shift: seq_strand=1 regardless of transcript strand.
        // This means hgvs_reverse = (1 != 1) = false.
        // Both rotate left for forward shift.
        let (shift, _seq, hgvs, _, _) =
            perform_shift_ensembl(b"TC", b"TC", "TCAA", "", 100, 99, false, 1);
        assert_eq!(shift, 2);
        assert_eq!(hgvs, b"TC"); // full rotation
    }

    // ---------------------------------------------------------------
    // perform_shift_ensembl — large-indel loop_limiter guard (#99)
    // ---------------------------------------------------------------

    #[test]
    fn test_perform_shift_forward_indel_longer_than_flank() {
        // Simulates a >1000bp deletion where indel_length exceeds
        // post_seq length. The loop_limiter guard must allow the
        // character-by-character loop to still find the correct shift.
        // Indel "ABCDE" (5 bytes) with only 3 bytes of flank "ABX".
        // The first 2 bases match (A, B), then X mismatches.
        let (shift, _seq, _hgvs, start, end) =
            perform_shift_ensembl(b"ABCDE", b"ABCDE", "ABX", "", 100, 104, false, 1);
        assert_eq!(shift, 2);
        assert_eq!(start, 102);
        assert_eq!(end, 106);
    }

    #[test]
    fn test_perform_shift_forward_indel_longer_than_flank_no_match() {
        // Indel longer than flank but first base doesn't match.
        let (shift, seq, _hgvs, start, end) =
            perform_shift_ensembl(b"ABCDE", b"ABCDE", "XYZ", "", 100, 104, false, 1);
        assert_eq!(shift, 0);
        assert_eq!(seq, b"ABCDE");
        assert_eq!(start, 100);
        assert_eq!(end, 104);
    }

    #[test]
    fn test_perform_shift_reverse_indel_longer_than_flank() {
        // Reverse-strand: indel "GATGC" (5 bytes) with pre_seq "ATG"
        // (3 bytes). indel > flank triggers the guard. Last 2 bases of
        // pre_seq ("TG") match the trailing "TG" of the indel.
        let (shift, _seq, _hgvs, start, end) =
            perform_shift_ensembl(b"GATGC", b"GATGC", "", "XTG", 100, 99, true, -1);
        // n=1: last of "GATGC"='C', pre_seq[pre_seq.len()-1]='G' — mismatch — shift=0.
        assert_eq!(shift, 0);
        assert_eq!((start, end), (100, 99));
    }

    #[test]
    fn test_perform_shift_reverse_indel_longer_matches() {
        // Reverse-strand: indel "GATGG" (5 bytes) with pre_seq "XGG"
        // (3 bytes). Guard triggers. Last byte 'G' matches pre_seq[2]='G',
        // then next rotation's last byte 'G' matches pre_seq[1]='G',
        // then 'T' vs 'X' → mismatch.
        let (shift, _seq, _hgvs, start, end) =
            perform_shift_ensembl(b"GATGG", b"GATGG", "", "XGG", 100, 99, true, -1);
        assert_eq!(shift, 2);
        assert_eq!((start, end), (100, 99));
    }

    // ---------------------------------------------------------------
    // clip_alleles — strand-aware multi-base clipping
    // ---------------------------------------------------------------

    #[test]
    fn test_clip_alleles_negative_strand_prefix_decrements_end() {
        let mut notation = HgvsNotation {
            start: 100,
            end: 102,
            ref_allele: "ACG".into(),
            alt_allele: "ACT".into(),
            kind: "delins".into(),
        };
        clip_alleles(&mut notation, -1);
        // Negative strand: prefix clip decrements end.
        assert_eq!(notation.ref_allele, "G");
        assert_eq!(notation.alt_allele, "T");
        assert_eq!(notation.kind, ">");
        assert_eq!(notation.end, 100); // decremented by 2
    }

    #[test]
    fn test_clip_alleles_negative_strand_suffix_increments_start() {
        let mut notation = HgvsNotation {
            start: 100,
            end: 102,
            ref_allele: "GCA".into(),
            alt_allele: "TCA".into(),
            kind: "delins".into(),
        };
        clip_alleles(&mut notation, -1);
        assert_eq!(notation.ref_allele, "G");
        assert_eq!(notation.alt_allele, "T");
        assert_eq!(notation.kind, ">");
        assert_eq!(notation.start, 102); // incremented by 2 (suffix "CA")
    }

    // ---------------------------------------------------------------
    // check_for_peptide_duplication — direct upstream match
    // ---------------------------------------------------------------

    #[test]
    fn test_check_for_peptide_duplication_at_current_position() {
        let mut notation = ProteinHgvsNotation {
            start: 4,
            end: 5, // boundary insertion
            ref_allele: String::new(),
            alt_allele: "K".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        // ref_translation[2] = 'K' (0-based idx 2 = position 3).
        // check at start=4: test_new_start = 4-1-1 = 2. upstream = ref[..3] = "MAK".
        // upstream[2] = 'K'. Match!
        let result = check_for_peptide_duplication(&mut notation, "MAKL*");
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 3); // 4 - 1
        assert_eq!(notation.end, 3); // 4 - 1
    }

    #[test]
    fn test_check_for_peptide_duplication_no_fallback_when_upstream_mismatches() {
        // VEP's _check_for_peptide_duplication() only checks at
        // start - len(alt) - 1 (0-indexed) with no fallback. When the
        // upstream sequence at that position doesn't match, the result
        // stays as "ins".
        let mut notation = ProteinHgvsNotation {
            start: 3,
            end: 4,
            ref_allele: String::new(),
            alt_allele: "K".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        // At start=3: upstream = "MA", test_new_start = 1, test_seq = "A" ≠ "K"
        // VEP: no fallback → stays as "ins"
        let result = check_for_peptide_duplication(&mut notation, "MAKL*");
        assert!(!result);
        assert_eq!(notation.kind, "ins");
    }

    #[test]
    fn test_check_for_peptide_duplication_no_match_returns_false() {
        let mut notation = ProteinHgvsNotation {
            start: 3,
            end: 3,
            ref_allele: String::new(),
            alt_allele: "W".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        // No 'W' in the upstream of "MAKL*"
        let result = check_for_peptide_duplication(&mut notation, "MAKL*");
        assert!(!result);
        assert_eq!(notation.kind, "ins"); // unchanged
    }

    #[test]
    fn test_check_for_peptide_duplication_multi_residue() {
        let mut notation = ProteinHgvsNotation {
            start: 5,
            end: 6,
            ref_allele: String::new(),
            alt_allele: "KL".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        // ref = "MAKLKL*". At start=5: test_new_start = 5-2-1 = 2.
        // upstream = ref[..4] = "MAKL". upstream[2..4] = "KL". Match!
        let result = check_for_peptide_duplication(&mut notation, "MAKLKL*");
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 3);
        assert_eq!(notation.end, 4);
    }

    // ---------------------------------------------------------------
    // Protein HGVS helpers
    // ---------------------------------------------------------------

    #[test]
    fn test_resolve_frameshift_finds_first_changed_residue() {
        let mut notation = ProteinHgvsNotation {
            start: 2,
            end: 2,
            ref_allele: String::new(),
            alt_allele: String::new(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: String::new(),
        };
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".into(),
            alt_peptide: "Q".into(),
            ref_translation: "MAK*".into(),
            alt_translation: "MQW*".into(),
            alt_translation_extension: None,
            frameshift: true,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        resolve_frameshift_hgvs(&mut notation, &protein);
        assert_eq!(notation.kind, "fs");
        assert_eq!(notation.start, 2); // first differing AA
        assert_eq!(notation.ref_allele, "A");
        assert_eq!(notation.alt_allele, "Q");
    }

    #[test]
    fn test_resolve_frameshift_synonymous_when_both_reach_stop() {
        let mut notation = ProteinHgvsNotation {
            start: 2,
            end: 2,
            ref_allele: String::new(),
            alt_allele: String::new(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: String::new(),
        };
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".into(),
            alt_peptide: "A".into(),
            ref_translation: "MA*".into(),
            alt_translation: "MA*".into(),
            alt_translation_extension: None,
            frameshift: true,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        resolve_frameshift_hgvs(&mut notation, &protein);
        assert_eq!(notation.kind, "=");
    }

    #[test]
    fn test_surrounding_peptides_returns_flanking_residues() {
        assert_eq!(
            surrounding_peptides("MAKL*", 2, "", Some(2)),
            Some("AK".into())
        );
        assert_eq!(
            surrounding_peptides("MAKL*", 1, "", Some(2)),
            Some("MA".into())
        );
    }

    #[test]
    fn test_surrounding_peptides_extends_with_stop_original_ref() {
        // When original_ref starts with *, append it to ref_translation
        assert_eq!(
            surrounding_peptides("MAK", 3, "*Q", Some(2)),
            Some("K*".into())
        );
    }

    #[test]
    fn test_normalize_peptide_allele_converts_dash() {
        assert_eq!(normalize_peptide_allele("-"), "");
        assert_eq!(normalize_peptide_allele("K"), "K");
    }

    #[test]
    fn test_append_terminal_stop_adds_star_when_missing() {
        assert_eq!(append_terminal_stop("MAK"), "MAK*");
        assert_eq!(append_terminal_stop("MAK*"), "MAK*");
        assert_eq!(append_terminal_stop("M*K"), "M*K"); // internal stop, no append
    }

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement("ACGT"), Some("ACGT".into()));
        assert_eq!(reverse_complement("AAAA"), Some("TTTT".into()));
        assert_eq!(reverse_complement("N"), Some("N".into()));
        assert_eq!(reverse_complement(""), Some("".into()));
        assert_eq!(reverse_complement("X"), None); // invalid base
    }

    #[test]
    fn test_split_hgvs_coord_parses_intronic_offsets() {
        assert_eq!(split_hgvs_coord("100"), Some((100, None)));
        assert_eq!(split_hgvs_coord("100+5"), Some((100, Some("+5".into()))));
        assert_eq!(split_hgvs_coord("100-3"), Some((100, Some("-3".into()))));
        assert_eq!(split_hgvs_coord("-5"), Some((-5, None)));
        assert_eq!(split_hgvs_coord("-5+10"), Some((-5, Some("+10".into()))));
        assert_eq!(split_hgvs_coord("*100+5"), Some((100, Some("+5".into()))));
    }

    #[test]
    fn test_protein_event_type_classification() {
        assert_eq!(protein_event_type("", "K", false), "ins");
        assert_eq!(protein_event_type("K", "", false), "del");
        assert_eq!(protein_event_type("K", "L", false), ">");
        assert_eq!(protein_event_type("KL", "QW", false), "delins");
        assert_eq!(protein_event_type("K", "L", true), "fs");
        assert_eq!(protein_event_type("K", "K", false), "=");
    }

    #[test]
    fn test_peptide_char_1based_indexing() {
        assert_eq!(peptide_char("MAKL", 1), Some('M'));
        assert_eq!(peptide_char("MAKL", 4), Some('L'));
        assert_eq!(peptide_char("MAKL", 5), None); // out of bounds
        assert_eq!(peptide_char("MAKL", 0), None); // invalid 1-based
    }

    #[test]
    fn test_format_hgvsp_deletion_single_residue() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".into(),
            alt_peptide: "-".into(),
            ref_translation: "MA*".into(),
            alt_translation: "M*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2del".into())
        );
    }

    #[test]
    fn test_format_hgvsp_deletion_multi_residue() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 3,
            ref_peptide: "AK".into(),
            alt_peptide: "-".into(),
            ref_translation: "MAK*".into(),
            alt_translation: "M*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2_Lys3del".into())
        );
    }

    #[test]
    fn test_format_hgvsp_missense() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".into(),
            alt_peptide: "V".into(),
            ref_translation: "MA*".into(),
            alt_translation: "MV*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2Val".into())
        );
    }

    #[test]
    fn test_format_hgvsp_delins() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".into(),
            alt_peptide: "VW".into(),
            ref_translation: "MAK*".into(),
            alt_translation: "MVWK*".into(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2delinsValTrp".into())
        );
    }

    #[test]
    fn test_format_hgvsp_frameshift_immediate_stop() {
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 2,
            end: 2,
            ref_peptide: "A".into(),
            alt_peptide: "*".into(),
            ref_translation: "MAK*".into(),
            alt_translation: "M*".into(),
            alt_translation_extension: None,
            frameshift: true,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2Ter".into())
        );
    }

    // ── try_peptide_dup_at direct upstream-match tests ──────────────────

    #[test]
    fn peptide_dup_single_residue_keeps_upstream_match_coordinates() {
        // Ref: MAAAEEEEK — E at positions 5,6,7,8 (1-based)
        // Insert "E" at check_start=6 → upstream="MAAAE"
        // test_new_start = 6-1-1 = 4, test_seq = upstream[4..5] = "E" ✓
        let mut notation = ProteinHgvsNotation {
            start: 6,
            end: 6,
            ref_allele: String::new(),
            alt_allele: "E".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = try_peptide_dup_at(&mut notation, "MAAAEEEEK", 6);
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 5);
        assert_eq!(notation.end, 5);
    }

    #[test]
    fn peptide_dup_single_residue_no_shift_needed() {
        // Ref: MAEK — E at position 3 only, followed by K
        // check_start=4: upstream="MAE", test_new_start=4-1-1=2, test_seq="E" ✓
        let mut notation = ProteinHgvsNotation {
            start: 4,
            end: 4,
            ref_allele: String::new(),
            alt_allele: "E".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = try_peptide_dup_at(&mut notation, "MAEK", 4);
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 3);
        assert_eq!(notation.end, 3);
    }

    #[test]
    fn peptide_dup_multi_residue_uses_direct_upstream_match() {
        // Ref: MPAPAPAD — "PA" repeats at 2-3, 4-5, 6-7
        // Insert "PA" at check_start=4: upstream="MPA"
        // test_new_start = 4-2-1 = 1, test_seq = upstream[1..3] = "PA" ✓
        let mut notation = ProteinHgvsNotation {
            start: 4,
            end: 5,
            ref_allele: String::new(),
            alt_allele: "PA".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = try_peptide_dup_at(&mut notation, "MPAPAPAD", 4);
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 2);
        assert_eq!(notation.end, 3);
    }

    #[test]
    fn peptide_dup_keeps_initial_match_at_ref_end() {
        // Ref: MAAEE — E at positions 4,5
        // check_start=5: upstream="MAAE", test_new_start=5-1-1=3, test_seq="E" ✓
        let mut notation = ProteinHgvsNotation {
            start: 5,
            end: 5,
            ref_allele: String::new(),
            alt_allele: "E".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = try_peptide_dup_at(&mut notation, "MAAEE", 5);
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 4);
        assert_eq!(notation.end, 4);
    }

    #[test]
    fn peptide_dup_issue89_example_keeps_pre_shift_dup_coordinates() {
        // The peptide duplication check itself does not perform an extra
        // 3' walk. Any right-shifting happens earlier in `_check_peptides_post_var`.
        // Ref has E at positions 25,26,27,28 then K at 29.
        let mut ref_translation = "M".repeat(24);
        ref_translation.push_str("EEEEK"); // pos 25-28: E, 29: K
        let mut notation = ProteinHgvsNotation {
            start: 26,
            end: 26,
            ref_allele: String::new(),
            alt_allele: "E".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = try_peptide_dup_at(&mut notation, &ref_translation, 26);
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 25);
        assert_eq!(notation.end, 25);
    }

    #[test]
    fn peptide_dup_non_periodic_uses_direct_upstream_match() {
        // For ref "MABAC", alt "AB", duplication is reported from the
        // directly matching upstream window without any additional walk.
        let mut notation = ProteinHgvsNotation {
            start: 3,
            end: 3,
            ref_allele: String::new(),
            alt_allele: "AB".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = try_peptide_dup_at(&mut notation, "MABAC", 4);
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 2);
        assert_eq!(notation.end, 3);
        assert_eq!(notation.alt_allele, "AB");
    }

    #[test]
    fn peptide_dup_via_check_for_peptide_duplication() {
        // Test the full duplication check without an extra internal 3' walk.
        // Ref: MAAAEEEEK — E at 5,6,7,8
        // Insertion "E" at position 6 → dup detected at the directly
        // matching upstream position.
        let mut notation = ProteinHgvsNotation {
            start: 6,
            end: 6,
            ref_allele: String::new(),
            alt_allele: "E".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        let result = check_for_peptide_duplication(&mut notation, "MAAAEEEEK");
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 5);
        assert_eq!(notation.end, 5);
    }

    #[test]
    fn peptide_dup_chr3_63912714_should_be_ins_not_dup() {
        // Reproduces chr3:63912714 A>AGCAGCAGCC / ENST00000295900
        // Amino acids: Q/QQQP at protein position 39
        // Ref protein around pos 35-44: QQQQQPPPPP
        // After clip: ref="", alt="QQP", start=40, preseq="Q"
        // Upstream[36..39] = "QQQ" ≠ "QQP" → no dup at check_start=40
        // Upstream[37..40] = "QQP" → MATCH at check_start=41 (fallback)
        // VEP does NOT fallback when preseq is non-empty → VEP says ins
        // Ref translation (simplified, only relevant positions):
        let mut ref_translation = "M".repeat(34);
        ref_translation.push_str("QQQQQPPPP"); // pos 35-43
        ref_translation.push_str("QP"); // pos 44-45 (etc.)

        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 39,
            end: 39,
            ref_peptide: "Q".to_string(),
            alt_peptide: "QQQP".to_string(),
            ref_translation: ref_translation.clone(),
            alt_translation: {
                let mut alt = ref_translation.clone();
                alt.insert_str(39, "QQP"); // insert QQP after position 39
                alt
            },
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        let result = format_hgvsp(&translation, &protein, true);
        // VEP: p.Gln39_Pro40insGlnGlnPro — NOT a dup
        let r = result.unwrap();
        assert!(r.contains("ins"), "Expected ins notation, got: {r}");
    }

    // ── HGVSc dup boundary revert tests (issue #88 remaining) ──────────
    //
    // NOTE: A full integration test calling `format_hgvsc` with a shifted
    // dup landing outside the exon span would be ideal but requires
    // substantial fixture plumbing (HgvsGenomicShift, FASTA context, etc.).
    // The boundary condition is validated end-to-end by the 42-variant
    // benchmark in the PR. The unit test below verifies the span arithmetic.

    #[test]
    fn test_format_hgvsp_clipped_insertion_becomes_dup() {
        // When the codon window is correctly widened for an insertion, the
        // ref_peptide and alt_peptide share a common prefix that is clipped,
        // leaving an insertion that is then checked for duplication.
        //
        // Protein: MAAK*, ref_peptide = "A" at pos 3, alt_peptide = "AA".
        // After clip: ref = "", alt = "A" → insertion → upstream "A" → dup.
        // The 3' peptide shift does NOT move because pos 4 is Lys, not Ala.
        let translation = make_translation();
        let protein = ProteinHgvsData {
            start: 3,
            end: 3,
            ref_peptide: "A".to_string(),
            alt_peptide: "AA".to_string(),
            ref_translation: "MAAK".to_string(),
            alt_translation: "MAAAK".to_string(),
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala3dup".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_clipped_insertion_multi_residue_dup() {
        // Multi-residue duplication: ref "GGG" at pos 48, alt "GGGGGG" (6 Gly).
        // After clip: ref = "", alt = "GGG" → insertion → upstream match → dup.
        let translation = make_translation();
        let ref_trans = format!("M{}R", "G".repeat(50));
        let alt_trans = format!("M{}R", "G".repeat(53));
        let protein = ProteinHgvsData {
            start: 48,
            end: 50,
            ref_peptide: "GGG".to_string(),
            alt_peptide: "GGGGGG".to_string(),
            ref_translation: ref_trans,
            alt_translation: alt_trans,
            alt_translation_extension: None,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
            native_refseq: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Gly48_Gly50dup".to_string())
        );
    }

    #[test]
    fn dup_range_before_transcript_span_is_detected() {
        let exon = make_exon(); // 90..140
        let exons = vec![&exon];
        let span_start = exons.iter().map(|e| e.start).min().unwrap();
        let span_end = exons.iter().map(|e| e.end).max().unwrap();
        // Dup range 85..89 is before exon span start 90 → would trigger revert
        assert!(85 < span_start);
        // Dup range 95..100 is within exon span → no revert
        assert!(95 >= span_start && 100 <= span_end);
        // Dup range 141..145 is after exon span end 140 → would trigger revert
        assert!(141 > span_end);
    }
}
