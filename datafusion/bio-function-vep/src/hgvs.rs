//! HGVS notation formatting for VEP CSQ fields (HGVSc and HGVSp).

use std::cmp::Ordering;

use crate::transcript_consequence::{ExonFeature, TranscriptFeature, TranslationFeature};

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

/// Format versioned transcript ID: "ENST00000379410.6".
fn versioned_id(base_id: &str, version: Option<i32>) -> String {
    match version {
        Some(v) => format!("{base_id}.{v}"),
        None => base_id.to_string(),
    }
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
    _cdna_position: Option<&str>,
    ref_allele: &str,
    alt_allele: &str,
    variant_start: i64,
    variant_end: i64,
) -> Option<String> {
    let tx_id = versioned_id(&tx.transcript_id, tx.version);
    let numbering = if tx.cds_start.is_some() && tx.cds_end.is_some() {
        'c'
    } else {
        'n'
    };
    let mut notation =
        hgvs_variant_notation(ref_allele, alt_allele, variant_start, variant_end)?;
    if notation.kind != "dup" {
        clip_alleles(&mut notation);
    }
    let (mut start, mut end) = notation_to_hgvsc_coords(tx, tx_exons, &notation)?;
    if !end.starts_with('*')
        && matches!(
            compare_hgvs_positions(&start, &end),
            Some(Ordering::Greater)
        )
    {
        std::mem::swap(&mut start, &mut end);
    }
    format_hgvs_string(&tx_id, numbering, &start, &end, &notation)
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
fn clip_alleles(notation: &mut HgvsNotation) {
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
        start += 1;
    }

    while !ref_allele.is_empty()
        && !alt_allele.is_empty()
        && ref_allele.as_bytes()[ref_allele.len() - 1] == alt_allele.as_bytes()[alt_allele.len() - 1]
    {
        ref_allele.pop();
        alt_allele.pop();
        end -= 1;
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
    let suffix = if notation.kind == ">" || (notation.kind == "inv" && notation.ref_allele.len() == 1)
    {
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

fn coding_cdna_bounds(tx: &TranscriptFeature, tx_exons: &[&ExonFeature]) -> Option<(usize, usize)> {
    let cds_start = tx.cds_start?;
    let cds_end = tx.cds_end?;
    let coding_start_anchor = if tx.strand >= 0 { cds_start } else { cds_end };
    let coding_end_anchor = if tx.strand >= 0 { cds_end } else { cds_start };
    Some((
        genomic_to_cdna_index_pub(tx_exons, tx.strand, coding_start_anchor)?,
        genomic_to_cdna_index_pub(tx_exons, tx.strand, coding_end_anchor)?,
    ))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExonCdnaCoord {
    start: i64,
    end: i64,
    cdna_start: usize,
    cdna_end: usize,
}

fn exon_cdna_coords(tx_exons: &[&ExonFeature], strand: i8) -> Option<Vec<ExonCdnaCoord>> {
    let mut exons: Vec<&ExonFeature> = tx_exons.to_vec();
    exons.sort_by_key(|exon| exon.start);
    if exons.is_empty() {
        return None;
    }

    let exon_lens: Vec<usize> = exons
        .iter()
        .map(|exon| usize::try_from(exon.end.saturating_sub(exon.start).saturating_add(1)).ok())
        .collect::<Option<Vec<_>>>()?;
    let total_len = exon_lens.iter().copied().sum::<usize>();

    let mut coords = Vec::with_capacity(exons.len());
    if strand >= 0 {
        let mut offset = 0usize;
        for (exon, exon_len) in exons.into_iter().zip(exon_lens) {
            let cdna_start = offset + 1;
            let cdna_end = offset + exon_len;
            coords.push(ExonCdnaCoord {
                start: exon.start,
                end: exon.end,
                cdna_start,
                cdna_end,
            });
            offset = cdna_end;
        }
    } else {
        let mut consumed = 0usize;
        for (exon, exon_len) in exons.into_iter().zip(exon_lens) {
            let cdna_end = total_len.saturating_sub(consumed);
            let cdna_start = cdna_end.saturating_sub(exon_len).saturating_add(1);
            coords.push(ExonCdnaCoord {
                start: exon.start,
                end: exon.end,
                cdna_start,
                cdna_end,
            });
            consumed = consumed.saturating_add(exon_len);
        }
    }

    Some(coords)
}

fn hgvs_cdna_position_from_genomic(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    genomic_pos: i64,
) -> Option<String> {
    let exons = exon_cdna_coords(tx_exons, tx.strand)?;
    let mut cdna_position = None;

    for (i, exon) in exons.iter().enumerate() {
        if genomic_pos > exon.end {
            continue;
        }

        if genomic_pos >= exon.start {
            let coord = if tx.strand >= 0 {
                exon.cdna_start as i64 + (genomic_pos - exon.start)
            } else {
                exon.cdna_start as i64 + (exon.end - genomic_pos)
            };
            cdna_position = Some(coord.to_string());
            break;
        }

        let prev_exon = exons.get(i.checked_sub(1)?)?;
        let updist = (genomic_pos - prev_exon.end).abs();
        let downdist = (exon.start - genomic_pos).abs();
        cdna_position = Some(
            if updist < downdist || (updist == downdist && tx.strand >= 0) {
                if tx.strand >= 0 {
                    format!("{}+{}", prev_exon.cdna_end, updist)
                } else {
                    format!("{}-{}", prev_exon.cdna_start, updist)
                }
            } else if tx.strand >= 0 {
                format!("{}-{}", exon.cdna_start, downdist)
            } else {
                format!("{}+{}", exon.cdna_end, downdist)
            },
        );
        break;
    }

    shift_to_hgvs_coding_coordinates(tx, tx_exons, &cdna_position?)
}

fn shift_to_hgvs_coding_coordinates(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    raw_cdna_position: &str,
) -> Option<String> {
    let (raw_coord, intron_offset) = split_hgvs_coord(raw_cdna_position)?;
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

/// Public wrapper for genomic_to_cdna_index (which is private in transcript_consequence).
fn genomic_to_cdna_index_pub(tx_exons: &[&ExonFeature], strand: i8, pos: i64) -> Option<usize> {
    let mut segments: Vec<(i64, i64)> = tx_exons.iter().map(|e| (e.start, e.end)).collect();
    segments.sort_by_key(|(s, _)| *s);
    if strand < 0 {
        segments.reverse();
    }
    let mut offset = 0usize;
    for (seg_start, seg_end) in segments {
        let seg_len = usize::try_from(seg_end.saturating_sub(seg_start).saturating_add(1)).ok()?;
        if pos >= seg_start && pos <= seg_end {
            let local = if strand >= 0 {
                usize::try_from(pos.saturating_sub(seg_start)).ok()?
            } else {
                usize::try_from(seg_end.saturating_sub(pos)).ok()?
            };
            return Some(offset + local + 1); // 1-based
        }
        offset = offset.saturating_add(seg_len);
    }
    None
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
        return Some(format!("{protein_id}:p.?"));
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
        if shift_hgvs && matches!(notation.kind.as_str(), "ins" | "del") {
            shift_peptides_post_var(&mut notation, &protein.ref_translation);
        }
        if notation.kind == "ins" && !check_for_peptide_duplication(&mut notation, &protein.ref_translation) {
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
        && ref_allele.as_bytes()[ref_allele.len() - 1] == alt_allele.as_bytes()[alt_allele.len() - 1]
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
    peptide.as_bytes().get(pos.checked_sub(1)?).map(|b| *b as char)
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
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2371-L2404
fn check_for_peptide_duplication(notation: &mut ProteinHgvsNotation, ref_translation: &str) -> bool {
    if notation.alt_allele.is_empty() || notation.start == 0 {
        return false;
    }

    let mut upstream = ref_translation
        .get(..notation.start.saturating_sub(1))
        .unwrap_or_default()
        .to_string();
    upstream.push_str(&notation.preseq);

    let alt_len = notation.alt_allele.len();
    let Some(test_new_start) = notation
        .start
        .checked_sub(alt_len)
        .and_then(|s| s.checked_sub(1))
    else {
        return false;
    };
    let Some(test_seq) = upstream.get(test_new_start..test_new_start.saturating_add(alt_len)) else {
        return false;
    };

    if test_seq == notation.alt_allele {
        notation.kind = "dup".to_string();
        notation.end = notation.start.saturating_sub(1);
        notation.start = notation.start.saturating_sub(alt_len);
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
        Some(len) => ref_trans.get(start..start.saturating_add(len)).map(str::to_string),
        None => ref_trans.get(start..).map(str::to_string),
    }
}

/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::_stop_loss_extra_AA()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2406-L2455
fn stop_loss_extra_aa(protein: &ProteinHgvsData, ref_var_pos: usize, frameshift: bool) -> Option<usize> {
    let stop_idx = protein.alt_translation.find('*')?;
    let extra = if frameshift {
        stop_idx.saturating_add(1).checked_sub(ref_var_pos)?
    } else {
        let ref_len = protein
            .ref_translation
            .find('*')
            .unwrap_or(protein.ref_translation.len());
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

    if notation.ref_allele == notation.alt_allele && !matches!(notation.kind.as_str(), "fs" | "ins") {
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
            mature_mirna_regions: Vec::new(),
            gene_stable_id: None,
            gene_symbol: None,
            gene_symbol_source: None,
            gene_hgnc_id: None,
            source: None,
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
        clip_alleles(&mut notation);
        assert_eq!(notation.kind, ">");
        assert_eq!(notation.start, 101);
        assert_eq!(notation.end, 101);
        assert_eq!(notation.ref_allele, "C");
        assert_eq!(notation.alt_allele, "T");
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
        clip_alleles(&mut notation);
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
        clip_alleles(&mut notation);
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
            Some("ENSPHGVS000001.1:p.?".to_string())
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
            format_hgvsc(&tx, &exons, Some("14"), "G", "A", 103, 103),
            Some("ENSTHGVS000001.1:c.4G>A".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_formats_insertions_with_flanking_coordinates() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, None, "-", "T", 103, 103),
            Some("ENSTHGVS000001.1:c.3_4insT".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_formats_deletions_from_genomic_span() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, None, "G", "-", 103, 103),
            Some("ENSTHGVS000001.1:c.4del".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_negative_utr_coordinate() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("10"), "A", "G", 99, 99),
            Some("ENSTHGVS000001.1:c.-1A>G".to_string())
        );
    }

    #[test]
    fn test_format_hgvsc_uses_non_coding_numbering() {
        let tx = make_transcript("lncRNA", 1, None, None);
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("14"), "G", "A", 103, 103),
            Some("ENSTHGVS000001.1:n.14G>A".to_string())
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
    fn test_format_hgvsc_uses_star_coordinate_in_three_prime_utr() {
        let tx = make_transcript("protein_coding", 1, Some(100), Some(108));
        let exon = make_exon();
        let exons = [&exon];
        assert_eq!(
            format_hgvsc(&tx, &exons, Some("21"), "A", "G", 110, 110),
            Some("ENSTHGVS000001.1:c.*2A>G".to_string())
        );
    }
}
