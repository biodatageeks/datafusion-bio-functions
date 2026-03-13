//! HGVS notation formatting for VEP CSQ fields (HGVSc and HGVSp).

use std::cmp::Ordering;

use crate::transcript_consequence::{ExonFeature, TranscriptFeature, TranslationFeature};

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
/// For exonic variants, uses cDNA position. For intronic variants, computes
/// offset from nearest exon boundary.
///
/// Traceability:
/// - Ensembl Variation `TranscriptVariationAllele::hgvs_transcript()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1301-L1491
/// - Ensembl Variation `TranscriptVariationAllele::_get_cDNA_position()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2683-L2765
/// - Ensembl Variation `Utils::Sequence::format_hgvs_string()`
///   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L623-L676
pub fn format_hgvsc(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    cdna_position: Option<&str>,
    ref_allele: &str,
    alt_allele: &str,
    variant_start: i64,
    variant_end: i64,
) -> Option<String> {
    let tx_id = versioned_id(&tx.transcript_id, tx.version);
    let is_ins = ref_allele == "-";
    let is_del = alt_allele == "-";
    let numbering = if tx.cds_start.is_some() && tx.cds_end.is_some() {
        'c'
    } else {
        'n'
    };

    // Try exonic first (we have a cDNA position).
    if let Some(cdna_pos) = cdna_position {
        let hgvs_pos = exonic_hgvsc_position(tx, tx_exons, cdna_pos)?;
        let change = format_cdna_change(&hgvs_pos, ref_allele, alt_allele, is_ins, is_del);
        return Some(format!("{tx_id}:{numbering}.{change}"));
    }

    // Intronic: compute offset from nearest exon boundary.
    if let Some(intronic) = compute_intronic_hgvsc(
        tx,
        tx_exons,
        variant_start,
        variant_end,
        ref_allele,
        alt_allele,
    ) {
        return Some(format!("{tx_id}:{numbering}.{intronic}"));
    }

    None
}

fn exonic_hgvsc_position(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    cdna_position: &str,
) -> Option<String> {
    if tx.cds_start.is_none() || tx.cds_end.is_none() {
        return Some(cdna_position.replace('-', "_"));
    }

    let (coding_start, coding_end) = coding_cdna_bounds(tx, tx_exons)?;
    let mut converted = Vec::new();
    for raw_part in cdna_position.split('-') {
        let raw = raw_part.parse::<usize>().ok()?;
        converted.push(convert_cdna_index_to_hgvsc(raw, coding_start, coding_end));
    }
    match converted.as_slice() {
        [single] => Some(single.clone()),
        [start, end] => Some(format!("{start}_{end}")),
        _ => None,
    }
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

fn convert_cdna_index_to_hgvsc(raw: usize, coding_start: usize, coding_end: usize) -> String {
    if raw > coding_end {
        format!("*{}", raw - coding_end)
    } else if raw >= coding_start {
        (raw - coding_start + 1).to_string()
    } else {
        (raw as i64 - coding_start as i64).to_string()
    }
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

/// Format the cDNA change part for exonic variants.
fn format_cdna_change(
    cdna_pos: &str,
    ref_allele: &str,
    alt_allele: &str,
    is_ins: bool,
    is_del: bool,
) -> String {
    if is_ins {
        // Insertion: c.10_11insACG
        format!("{cdna_pos}ins{alt_allele}")
    } else if is_del && ref_allele.len() == 1 {
        // Single-base deletion: c.10del
        format!("{cdna_pos}del")
    } else if is_del {
        // Multi-base deletion: c.10_12del
        format!("{cdna_pos}del")
    } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
        // SNV: c.1043G>A
        format!("{cdna_pos}{ref_allele}>{alt_allele}")
    } else if ref_allele.len() > 1 && alt_allele.len() > 1 {
        // Delins (MNV or complex): c.10_12delinsACG
        format!("{cdna_pos}delins{alt_allele}")
    } else if alt_allele.len() > ref_allele.len() {
        // Insertion-like (after trimming this shouldn't happen, but handle it)
        format!("{cdna_pos}delins{alt_allele}")
    } else {
        // Deletion-like with alt shorter than ref
        format!("{cdna_pos}delins{alt_allele}")
    }
}

/// Compute intronic HGVS notation: c.{cdna_boundary}{+/-}{offset}{ref}>{alt}
///
/// Finds the nearest exon boundary and computes the signed offset from it.
fn compute_intronic_hgvsc(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    variant_start: i64,
    variant_end: i64,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<String> {
    if tx_exons.is_empty() {
        return None;
    }

    let is_ins = ref_allele == "-";
    let is_del = alt_allele == "-";

    if is_ins {
        let mut pos1 = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start - 1)?;
        let mut pos2 = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start)?;
        if !pos1.starts_with('*')
            && !pos2.starts_with('*')
            && matches!(
                compare_hgvs_positions(&pos1, &pos2),
                Some(Ordering::Greater)
            )
        {
            std::mem::swap(&mut pos1, &mut pos2);
        }
        Some(format!("{pos1}_{pos2}ins{alt_allele}"))
    } else if is_del {
        if variant_start == variant_end {
            let pos = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start)?;
            Some(format!("{pos}del"))
        } else {
            let mut pos1 = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start)?;
            let mut pos2 = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_end)?;
            if !pos1.starts_with('*')
                && !pos2.starts_with('*')
                && matches!(
                    compare_hgvs_positions(&pos1, &pos2),
                    Some(Ordering::Greater)
                )
            {
                std::mem::swap(&mut pos1, &mut pos2);
            }
            Some(format!("{pos1}_{pos2}del"))
        }
    } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
        let pos = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start)?;
        Some(format!("{pos}{ref_allele}>{alt_allele}"))
    } else {
        if variant_start == variant_end {
            let pos = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start)?;
            Some(format!("{pos}delins{alt_allele}"))
        } else {
            let mut pos1 = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_start)?;
            let mut pos2 = hgvs_cdna_position_from_genomic(tx, tx_exons, variant_end)?;
            if !pos1.starts_with('*')
                && !pos2.starts_with('*')
                && matches!(
                    compare_hgvs_positions(&pos1, &pos2),
                    Some(Ordering::Greater)
                )
            {
                std::mem::swap(&mut pos1, &mut pos2);
            }
            Some(format!("{pos1}_{pos2}delins{alt_allele}"))
        }
    }
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

/// Compute HGVSp notation from coding classification data.
///
/// Format: `ENSP00000368698.2:p.Arg348His`
pub fn format_hgvsp(
    translation: &TranslationFeature,
    protein_position: Option<&str>,
    amino_acids: Option<&str>,
    terms: &[crate::so_terms::SoTerm],
) -> Option<String> {
    use crate::so_terms::SoTerm;

    let protein_id = versioned_id(translation.stable_id.as_deref()?, translation.version);

    let aa = amino_acids?;
    let pos_str = protein_position?;

    // Parse amino acids: "V/I" → (ref_aa, alt_aa)
    let parts: Vec<&str> = aa.split('/').collect();
    if parts.is_empty() {
        return None;
    }
    let ref_aa = parts[0];
    let alt_aa = if parts.len() > 1 { parts[1] } else { ref_aa };

    // Check for frameshift
    let is_frameshift = terms.contains(&SoTerm::FrameshiftVariant);
    if is_frameshift {
        return format_hgvsp_frameshift(&protein_id, pos_str, ref_aa, alt_aa);
    }

    // Check for stop gained
    if terms.contains(&SoTerm::StopGained) {
        return format_hgvsp_stop_gained(&protein_id, pos_str, ref_aa, alt_aa);
    }

    // Check for synonymous / stop retained / start retained
    let is_synonymous = terms.contains(&SoTerm::SynonymousVariant)
        || terms.contains(&SoTerm::StopRetainedVariant)
        || terms.contains(&SoTerm::StartRetainedVariant);
    if is_synonymous {
        return format_hgvsp_synonymous(&protein_id, pos_str, ref_aa);
    }

    // Check for inframe deletion
    if terms.contains(&SoTerm::InframeDeletion) {
        return format_hgvsp_inframe_del(&protein_id, pos_str, ref_aa);
    }

    // Check for inframe insertion
    if terms.contains(&SoTerm::InframeInsertion) {
        return format_hgvsp_inframe_ins(&protein_id, pos_str, ref_aa, alt_aa);
    }

    // Missense (default for protein-altering)
    if ref_aa.len() == 1 && alt_aa.len() == 1 {
        let ref3 = aa_one_to_three(ref_aa.chars().next()?);
        let alt3 = aa_one_to_three(alt_aa.chars().next()?);
        // Parse position (may be "348" or "348-350")
        let pos_num = pos_str.split('-').next()?;
        return Some(format!("{protein_id}:p.{ref3}{pos_num}{alt3}"));
    }

    // Multi-AA missense / complex
    if let Some(pos_num) = pos_str.split('-').next() {
        let ref3: String = ref_aa.chars().map(aa_one_to_three).collect();
        let alt3: String = alt_aa.chars().map(aa_one_to_three).collect();
        return Some(format!("{protein_id}:p.{ref3}{pos_num}{alt3}"));
    }

    None
}

fn format_hgvsp_frameshift(
    protein_id: &str,
    pos_str: &str,
    ref_aa: &str,
    alt_aa: &str,
) -> Option<String> {
    let pos_num = pos_str.split('-').next()?;
    let ref3 = aa_one_to_three(ref_aa.chars().next()?);
    // VEP shows the alt AA for frameshifts, defaulting to Ter if it's a stop
    let alt3 = if alt_aa == "*" || alt_aa.is_empty() {
        "Ter"
    } else {
        aa_one_to_three(alt_aa.chars().next()?)
    };
    // VEP format: p.Arg348AlafsTer? (we use Ter? since we don't compute exact stop)
    Some(format!("{protein_id}:p.{ref3}{pos_num}{alt3}fsTer?"))
}

fn format_hgvsp_stop_gained(
    protein_id: &str,
    pos_str: &str,
    ref_aa: &str,
    _alt_aa: &str,
) -> Option<String> {
    let pos_num = pos_str.split('-').next()?;
    let ref3 = aa_one_to_three(ref_aa.chars().next()?);
    Some(format!("{protein_id}:p.{ref3}{pos_num}Ter"))
}

fn format_hgvsp_synonymous(protein_id: &str, pos_str: &str, ref_aa: &str) -> Option<String> {
    let pos_num = pos_str.split('-').next()?;
    let ref3 = aa_one_to_three(ref_aa.chars().next()?);
    // VEP URL-encodes the "=" sign as %3D
    Some(format!("{protein_id}:p.{ref3}{pos_num}%3D"))
}

fn format_hgvsp_inframe_del(protein_id: &str, pos_str: &str, ref_aa: &str) -> Option<String> {
    let parts: Vec<&str> = pos_str.split('-').collect();
    if parts.len() == 2 {
        let start_pos = parts[0];
        let end_pos = parts[1];
        let first_ref3 = aa_one_to_three(ref_aa.chars().next()?);
        let last_ref3 = aa_one_to_three(ref_aa.chars().last()?);
        Some(format!(
            "{protein_id}:p.{first_ref3}{start_pos}_{last_ref3}{end_pos}del"
        ))
    } else {
        let ref3 = aa_one_to_three(ref_aa.chars().next()?);
        Some(format!("{protein_id}:p.{ref3}{pos_str}del"))
    }
}

fn format_hgvsp_inframe_ins(
    protein_id: &str,
    pos_str: &str,
    ref_aa: &str,
    alt_aa: &str,
) -> Option<String> {
    let parts: Vec<&str> = pos_str.split('-').collect();
    // The inserted sequence is the extra AAs in alt beyond ref
    let ins_seq: String = if alt_aa.len() > ref_aa.len() {
        alt_aa[ref_aa.len()..]
            .chars()
            .map(aa_one_to_three)
            .collect()
    } else {
        alt_aa.chars().map(aa_one_to_three).collect()
    };
    if parts.len() == 2 {
        let start_pos = parts[0];
        let end_pos = parts[1];
        let first_ref3 = aa_one_to_three(ref_aa.chars().next()?);
        let last_ref3 = if ref_aa.len() > 1 {
            aa_one_to_three(ref_aa.chars().last()?)
        } else {
            // For single-AA position, use the next position
            first_ref3
        };
        Some(format!(
            "{protein_id}:p.{first_ref3}{start_pos}_{last_ref3}{end_pos}ins{ins_seq}"
        ))
    } else {
        let ref3 = aa_one_to_three(ref_aa.chars().next()?);
        let pos_num: usize = pos_str.parse().ok()?;
        let next_pos = pos_num + 1;
        Some(format!(
            "{protein_id}:p.{ref3}{pos_num}_{ref3}{next_pos}ins{ins_seq}"
        ))
    }
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
    fn test_format_cdna_change_snv() {
        assert_eq!(
            format_cdna_change("1043", "G", "A", false, false),
            "1043G>A"
        );
    }

    #[test]
    fn test_format_cdna_change_insertion() {
        assert_eq!(
            format_cdna_change("10-11", "-", "ACG", true, false),
            "10-11insACG"
        );
    }

    #[test]
    fn test_format_cdna_change_deletion() {
        assert_eq!(format_cdna_change("10", "A", "-", false, true), "10del");
        assert_eq!(
            format_cdna_change("10-12", "ACG", "-", false, true),
            "10-12del"
        );
    }

    #[test]
    fn test_format_cdna_change_delins() {
        assert_eq!(
            format_cdna_change("10-12", "ACG", "TT", false, false),
            "10-12delinsTT"
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
        assert_eq!(
            format_hgvsp_synonymous("ENSP00000368698.2", "100", "V"),
            Some("ENSP00000368698.2:p.Val100%3D".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_stop_gained() {
        assert_eq!(
            format_hgvsp_stop_gained("ENSP00000368698.2", "100", "R", "*"),
            Some("ENSP00000368698.2:p.Arg100Ter".to_string())
        );
    }

    #[test]
    fn test_format_hgvsp_frameshift() {
        assert_eq!(
            format_hgvsp_frameshift("ENSP00000368698.2", "100", "R", "A"),
            Some("ENSP00000368698.2:p.Arg100AlafsTer?".to_string())
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
