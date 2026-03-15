//! HGVS notation formatting for VEP CSQ fields (HGVSc and HGVSp).

use std::cmp::Ordering;
use std::io::{BufRead, Seek};

use crate::transcript_consequence::{
    ExonFeature, TranscriptCdnaMapperSegment, TranscriptFeature, TranslationFeature,
    genomic_to_cdna_index_for_transcript, raw_cdna_position_from_genomic,
    unshifted_cdna_bounds_for_hgvs_shift,
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
    pub frameshift: bool,
    pub start_lost: bool,
    pub stop_lost: bool,
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
    pub alt_orig_allele_string: String,
    pub five_prime_context: String,
    pub three_prime_context: String,
}

impl HgvsGenomicShift {
    /// Traceability:
    /// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
    ///   applies `_hgvs_offset` to `_slice_start/_slice_end` instead of
    ///   mutating the stored reverse-strand genomic start/end in `shift_hash`
    ///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1360-L1406
    /// - Ensembl Variation `TranscriptVariationAllele::hgvs_protein()`
    ///   applies the same strand-aware `shifting_offset` to translation
    ///   coordinates when HGVS shifting is active
    ///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1636-L1669
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
///   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1417-L1421
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1301-L1491
/// - Ensembl Variation `Utils::Sequence::hgvs_variant_notation()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L493-L619
/// - Ensembl Variation `TranscriptVariationAllele::_clip_alleles()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2117-L2232
/// - Ensembl Variation `TranscriptVariationAllele::_get_cDNA_position()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2683-L2765
/// - Ensembl Variation `Utils::Sequence::format_hgvs_string()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L623-L676
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
    let use_genomic_shift = genomic_shift.filter(|_| ref_allele == "-" || alt_allele == "-");
    let edited_shifted_output_allele = use_genomic_shift.and_then(|_| {
        edited_transcript_shifted_output_allele(
            tx,
            tx_exons,
            ref_allele,
            alt_allele,
            variant_start,
            variant_end,
        )
    });
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
    let (feature_ref, feature_alt) = hgvs_feature_strand_alleles(tx, ref_allele, alt_allele)?;
    let mut notation =
        hgvs_variant_notation(&feature_ref, &feature_alt, variant_start, variant_end)?;
    if let Some(shift) = dup_context {
        let shifted_feature_alt = hgvs_feature_strand_alleles(
            tx,
            "-",
            edited_shifted_output_allele
                .as_deref()
                .unwrap_or(&shift.shifted_output_allele),
        )
        .map(|(_, alt)| alt)?;
        apply_shifted_insertion_duplication(tx, &shifted_feature_alt, shift, &mut notation);
    }
    if notation.kind != "dup" {
        clip_alleles(&mut notation, tx.strand);
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

#[derive(Clone, Copy)]
enum GenomicShiftKind {
    Insertion,
    Deletion,
}

/// Traceability:
/// - Ensembl Variation `_genomic_shift()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L411-L466
/// - Ensembl Variation `perform_shift()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351
/// - Ensembl Variation `create_shift_hash()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L365-L400
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
    let mut seq_to_check = seq_to_check.as_bytes().to_vec();
    let mut hgvs_output = alt_allele.as_bytes().to_vec();
    let mut shifted_start = start;
    let mut shifted_end = end;
    let shift_length;

    // Traceability:
    // - Ensembl Variation `TranscriptVariationAllele::_genomic_shift()`
    //   always passes `seq_strand = 1` to `perform_shift()` because the
    //   genomic shift operates on forward-strand coordinates; the HGVS
    //   output allele stays in VF orientation
    //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L431
    let genomic_seq_strand = 1i8;
    if strand >= 0 {
        let flank_start = end + 1;
        if flank_start <= 0 {
            return Ok(None);
        }
        let post_seq = read_reference_sequence(reader, chrom, flank_start, flank_start + 999)?;
        let pre_seq = String::new();
        let (computed_shift_length, shifted_seq, shifted_hgvs_output, new_start, new_end) =
            perform_shift_ensembl(
                &seq_to_check,
                &hgvs_output,
                &post_seq,
                &pre_seq,
                shifted_start,
                shifted_end,
                false,
                genomic_seq_strand,
            );
        shift_length = computed_shift_length;
        seq_to_check = shifted_seq;
        hgvs_output = shifted_hgvs_output;
        shifted_start = new_start;
        shifted_end = new_end;
    } else {
        let flank_end = start - 1;
        if flank_end <= 0 {
            return Ok(None);
        }
        let flank_start = (flank_end - 999).max(1);
        let pre_seq = read_reference_sequence(reader, chrom, flank_start, flank_end)?;
        let post_seq = String::new();
        let (computed_shift_length, shifted_seq, shifted_hgvs_output, new_start, new_end) =
            perform_shift_ensembl(
                &seq_to_check,
                &hgvs_output,
                &post_seq,
                &pre_seq,
                shifted_start,
                shifted_end,
                true,
                genomic_seq_strand,
            );
        shift_length = computed_shift_length;
        seq_to_check = shifted_seq;
        hgvs_output = shifted_hgvs_output;
        shifted_start = new_start;
        shifted_end = new_end;
    }

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
        alt_orig_allele_string: alt_allele.to_string(),
        five_prime_context,
        three_prime_context,
    }))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1360-L1364
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
/// - Ensembl Variation `TranscriptVariationAllele::_return_3prime()`
///   uses edited transcript `spliced_seq` for RefSeq HGVS shifting when RNA
///   edit attributes are present
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L149-L235
/// - Ensembl VEP `AnnotationType::Transcript::edit_transcript()` stores the
///   edited sequence on `_variation_effect_feature_cache.spliced_seq`
///   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L596-L597
fn edited_transcript_shifted_output_allele(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    ref_allele: &str,
    alt_allele: &str,
    variant_start: i64,
    variant_end: i64,
) -> Option<String> {
    if tx.bam_edit_status.as_deref() != Some("ok") || !tx.has_non_polya_rna_edit {
        return None;
    }
    if ref_allele != "-" || alt_allele.is_empty() || alt_allele == "-" {
        return None;
    }

    let edited_transcript_seq = tx.spliced_seq.as_deref()?;
    let (start, end) = unshifted_cdna_bounds_for_hgvs_shift(
        tx,
        tx_exons,
        variant_start,
        variant_end,
        ref_allele,
        alt_allele,
    )?;
    if start == 0 || end == 0 || start > end || end > edited_transcript_seq.len() {
        return None;
    }

    let search_start = start.saturating_sub(1001);
    let search_end = edited_transcript_seq.len().min(end.saturating_add(1000));
    let pre_seq = edited_transcript_seq[search_start..start.saturating_sub(1)].to_ascii_uppercase();
    let post_seq = edited_transcript_seq[end..search_end].to_ascii_uppercase();
    shift_output_allele_across_transcript(alt_allele, &pre_seq, &post_seq, tx.strand)
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_return_3prime()`
///   calls `perform_shift()` for edited RefSeq transcripts with the
///   transcript strand as `seq_strand` and the reverse flag derived
///   from the strand direction
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L237-L243
/// - Ensembl Variation `TranscriptVariationAllele::perform_shift()`
///   rotates the transcript-side allele and HGVS output allele separately
///   depending on transcript strand and variation-feature strand
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351
fn shift_output_allele_across_transcript(
    alt_allele: &str,
    pre_seq: &str,
    post_seq: &str,
    transcript_strand: i8,
) -> Option<String> {
    let vf_allele = alt_allele.to_ascii_uppercase();
    let feature_allele = if transcript_strand >= 0 {
        vf_allele.clone()
    } else {
        reverse_complement(&vf_allele)?
    };
    let seq_to_check = feature_allele.into_bytes();
    let hgvs_output = vf_allele.into_bytes();

    // VEP: $reverse = (-1 * ($strand - 1)) / 2
    let reverse = transcript_strand < 0;
    // VEP: perform_shift(..., $strand) — transcript strand as seq_strand
    let (_, _, shifted_output, _, _) = perform_shift_ensembl(
        &seq_to_check,
        &hgvs_output,
        post_seq,
        pre_seq,
        0, // var_start/end not used for output allele
        0,
        reverse,
        transcript_strand,
    );

    String::from_utf8(shifted_output).ok()
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L493-L619
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

    if ref_allele == alt_allele {
        return None;
    }

    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2117-L2232
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

fn variant_lies_within_transcript_sequence(
    tx_exons: &[&ExonFeature],
    start: i64,
    end: i64,
    is_insertion: bool,
) -> bool {
    if is_insertion {
        return tx_exons
            .iter()
            .any(|exon| start > exon.start && start <= exon.end);
    }

    let target_len = end.saturating_sub(start).saturating_add(1);
    let covered = tx_exons
        .iter()
        .map(|exon| overlap_len(start, end, exon.start, exon.end))
        .sum::<i64>();
    covered == target_len
}

/// Traceability:
/// - Ensembl Variation `perform_shift()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351
///
/// VCF alleles are always reported on the forward genomic strand. When HGVS is
/// generated for a reverse-strand transcript, Ensembl shifts the genomic
/// comparison sequence and the HGVS output allele in different directions.
fn perform_shift_ensembl(
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
    let loop_limiter = if reverse {
        pre_seq.len().saturating_sub(indel_length).saturating_add(1)
    } else {
        post_seq.len().saturating_sub(indel_length)
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L635-L676
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2683-L2765
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
    let cdna_position = raw_cdna_position_from_genomic(tx, tx_exons, genomic_pos)?;
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
    ///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1700-L1749
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1593-L1758
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_type()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1976-L2041
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2043-L2114
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1833-L1974
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
        // - Ensembl Variation `TranscriptVariationAllele::hgvs_protein()`
        //   checks for peptide duplication BEFORE 3' shifting
        //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1700-L1758
        if notation.kind == "ins"
            && check_for_peptide_duplication(&mut notation, &protein.ref_translation)
        {
            // Dup detected — skip shift and flanking.
        } else {
            if shift_hgvs && matches!(notation.kind.as_str(), "ins" | "del") {
                shift_peptides_post_var(&mut notation, &protein.ref_translation);
            }
            if notation.kind == "ins" {
                notation.ref_allele = surrounding_peptides(
                    &protein.ref_translation,
                    notation.start.min(notation.end),
                    &notation.original_ref,
                    Some(2),
                )?;
            }
        }
    }

    format_hgvsp_notation(&protein_id, &notation, protein)
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2054-L2089
fn normalize_peptide_allele(allele: &str) -> String {
    if allele == "-" {
        String::new()
    } else {
        allele.to_string()
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_clip_alleles()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2117-L2225
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1976-L2041
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2249-L2295
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2262-L2295
fn append_terminal_stop(peptide: &str) -> String {
    if peptide.contains('*') {
        peptide.to_string()
    } else {
        format!("{peptide}*")
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_fs_peptides()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2278-L2291
fn peptide_char(peptide: &str, pos: usize) -> Option<char> {
    peptide
        .as_bytes()
        .get(pos.checked_sub(1)?)
        .map(|b| *b as char)
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_check_peptides_post_var()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2501-L2518
/// - Ensembl Variation `TranscriptVariationAllele::_shift_3prime()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2525-L2573
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2371-L2410
/// - VEP's `$hgvs_notation->{start}` comes from `translation_start()` via the
///   `genomic2pep()` mapper, which for codon-boundary insertions can return a
///   position 1 higher than our `protein_position_start`. We compensate by also
///   trying with `start + 1` when the insertion has no clipped prefix (boundary).
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm#L467-L499
fn check_for_peptide_duplication(
    notation: &mut ProteinHgvsNotation,
    ref_translation: &str,
) -> bool {
    if notation.alt_allele.is_empty() || notation.start == 0 {
        return false;
    }

    // Try dup check at the current position first.
    if try_peptide_dup_at(notation, ref_translation, notation.start) {
        return true;
    }
    // For codon-boundary insertions (empty preseq = no clipped prefix),
    // also try one position forward to match VEP's genomic2pep mapper.
    if notation.preseq.is_empty() {
        if try_peptide_dup_at(notation, ref_translation, notation.start.saturating_add(1)) {
            return true;
        }
    }
    false
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2297-L2320
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2406-L2461
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
    let stop_idx = protein.alt_translation.find('*')?;
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

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_peptides()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2083-L2112
fn peptide_to_three_letter(peptide: &str) -> String {
    peptide.chars().map(aa_one_to_three).collect()
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1877-L1880
fn peptide_first_three(peptide: &str) -> Option<&'static str> {
    Some(aa_one_to_three(peptide.chars().next()?))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1877-L1880
fn peptide_last_three(peptide: &str) -> Option<&'static str> {
    Some(aa_one_to_three(peptide.chars().last()?))
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_get_hgvs_protein_format()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1833-L1974
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
            peptide_first_three(&notation.ref_allele)?,
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
            gene_hgnc_id: None,
            source: None,
            bam_edit_status: None,
            has_non_polya_rna_edit: false,
            spliced_seq: None,
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2=".to_string())
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
            frameshift: false,
            start_lost: true,
            stop_lost: false,
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
            frameshift: true,
            start_lost: false,
            stop_lost: false,
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
            frameshift: false,
            start_lost: false,
            stop_lost: true,
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
            frameshift: false,
            start_lost: false,
            stop_lost: true,
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
            start: 3,
            end: 3,
            ref_peptide: "-".to_string(),
            alt_peptide: "A".to_string(),
            ref_translation: "MAA*".to_string(),
            alt_translation: "MAAA*".to_string(),
            frameshift: false,
            start_lost: false,
            stop_lost: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2dup".to_string())
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
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
            alt_orig_allele_string: "-".to_string(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "AA", "-", 104, 105, Some(&shift)),
            Some("ENSTHGVS000001.1:c.11-3_11-2del".to_string())
        );
    }

    #[test]
    fn test_edited_transcript_shifted_output_allele_uses_spliced_seq_for_intronic_insertions() {
        let mut tx = make_transcript("protein_coding", 1, Some(95), Some(118));
        tx.bam_edit_status = Some("ok".to_string());
        tx.has_non_polya_rna_edit = true;
        tx.spliced_seq = Some("TTTTTTTTTTTGCTTTTTTT".to_string());
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
            edited_transcript_shifted_output_allele(&tx, &exons, "-", "GAA", 104, 104),
            Some("AAG".to_string())
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
            alt_orig_allele_string: "T".to_string(),
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
            alt_orig_allele_string: "AGTA".to_string(),
            five_prime_context: String::new(),
            three_prime_context: "AGTA".to_string(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "-", "AAGT", 148, 147, Some(&shift)),
            Some("ENSTHGVS000001.1:n.11+41_11+44dup".to_string())
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
            alt_orig_allele_string: "A".to_string(),
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
            alt_orig_allele_string: "A".to_string(),
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
            alt_orig_allele_string: "A".to_string(),
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
            alt_orig_allele_string: "T".to_string(),
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
            alt_orig_allele_string: "-".to_string(),
            five_prime_context: String::new(),
            three_prime_context: String::new(),
        };
        assert_eq!(
            format_hgvsc(&tx, &exons, None, None, "GTGT", "-", 120, 123, Some(&shift)),
            Some("ENSTHGVS000001.1:n.51_54del".to_string())
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
        assert_eq!(
            hgvs_cdna_position_from_genomic(&tx, &exons, 41383346),
            Some("2842".to_string())
        );
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
        assert_eq!(shift_length, 1);
        assert_eq!(String::from_utf8(shifted_seq).unwrap(), "GGAT");
        assert_eq!(String::from_utf8(shifted_output).unwrap(), "ATGG");
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
            frameshift: false,
            start_lost: false,
            stop_lost: true,
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
            frameshift: false,
            start_lost: false,
            stop_lost: true,
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
            frameshift: false,
            start_lost: false,
            stop_lost: true,
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
            frameshift: false,
            start_lost: false,
            stop_lost: true,
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
            frameshift: true,
            start_lost: false,
            stop_lost: false,
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
    // check_for_peptide_duplication — fallback positions
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
    fn test_check_for_peptide_duplication_fallback_at_offset_2() {
        let mut notation = ProteinHgvsNotation {
            start: 3,
            end: 4,
            ref_allele: String::new(),
            alt_allele: "K".into(),
            original_ref: String::new(),
            preseq: String::new(),
            kind: "ins".into(),
        };
        // At start=3: test_new_start = 3-1-1 = 1. upstream = ref[..2] = "MA".
        // upstream[1] = 'A'. Not 'K'. Fail.
        // At start+1=4: test_new_start = 4-1-1 = 2. upstream = ref[..3] = "MAK".
        // upstream[2] = 'K'. Match!
        // At start+2=5: test_new_start = 5-1-1 = 3. upstream = ref[..4] = "MAKL".
        // upstream[3] = 'L'. Not 'K'. Fail.
        // Fallback tries offset 2 first (reversed), then 1. Offset 2 (start=5) fails.
        // Offset 1 (start=4) succeeds.
        let result = check_for_peptide_duplication(&mut notation, "MAKL*");
        assert!(result);
        assert_eq!(notation.kind, "dup");
        assert_eq!(notation.start, 3); // 4 - 1
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
        assert_eq!(notation.start, 3); // 5 - 2
        assert_eq!(notation.end, 4); // 5 - 1
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
            frameshift: true,
            start_lost: false,
            stop_lost: false,
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
            frameshift: true,
            start_lost: false,
            stop_lost: false,
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
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
            frameshift: false,
            start_lost: false,
            stop_lost: false,
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
            frameshift: true,
            start_lost: false,
            stop_lost: false,
        };
        assert_eq!(
            format_hgvsp(&translation, &protein, true),
            Some("ENSPHGVS000001.1:p.Ala2Ter".into())
        );
    }
}
