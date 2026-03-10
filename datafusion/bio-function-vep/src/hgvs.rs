//! HGVS notation formatting for VEP CSQ fields (HGVSc and HGVSp).

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

    // Try exonic first (we have a cDNA position).
    if let Some(cdna_pos) = cdna_position {
        let change = format_cdna_change(cdna_pos, ref_allele, alt_allele, is_ins, is_del);
        return Some(format!("{tx_id}:c.{change}"));
    }

    // Intronic: compute offset from nearest exon boundary.
    if let Some(intronic) = compute_intronic_hgvsc(tx, tx_exons, variant_start, variant_end, ref_allele, alt_allele) {
        return Some(format!("{tx_id}:c.{intronic}"));
    }

    None
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

    // Get sorted exon boundaries.
    let mut exon_boundaries: Vec<(i64, i64)> = tx_exons.iter().map(|e| (e.start, e.end)).collect();
    exon_boundaries.sort_by_key(|(s, _)| *s);

    // Find the nearest exon boundary and compute the offset.
    let (cdna_boundary, offset) = nearest_exon_boundary(
        tx, tx_exons, &exon_boundaries, variant_start, variant_end,
    )?;

    if is_ins {
        // Intronic insertion
        let (cdna_boundary2, offset2) = nearest_exon_boundary(
            tx, tx_exons, &exon_boundaries, variant_start - 1, variant_start - 1,
        )?;
        let pos1 = format_intronic_pos(cdna_boundary2, offset2);
        let pos2 = format_intronic_pos(cdna_boundary, offset);
        // Order positions correctly
        Some(format!("{pos1}_{pos2}ins{alt_allele}"))
    } else if is_del {
        if variant_start == variant_end {
            let pos = format_intronic_pos(cdna_boundary, offset);
            Some(format!("{pos}del"))
        } else {
            let (cdna_boundary2, offset2) = nearest_exon_boundary(
                tx, tx_exons, &exon_boundaries, variant_end, variant_end,
            )?;
            let pos1 = format_intronic_pos(cdna_boundary, offset);
            let pos2 = format_intronic_pos(cdna_boundary2, offset2);
            Some(format!("{pos1}_{pos2}del"))
        }
    } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
        // SNV in intron
        let pos = format_intronic_pos(cdna_boundary, offset);
        Some(format!("{pos}{ref_allele}>{alt_allele}"))
    } else {
        // Complex intronic change
        let pos = format_intronic_pos(cdna_boundary, offset);
        if variant_start == variant_end {
            Some(format!("{pos}delins{alt_allele}"))
        } else {
            let (cdna_boundary2, offset2) = nearest_exon_boundary(
                tx, tx_exons, &exon_boundaries, variant_end, variant_end,
            )?;
            let pos2 = format_intronic_pos(cdna_boundary2, offset2);
            Some(format!("{pos}_{pos2}delins{alt_allele}"))
        }
    }
}

/// Format an intronic position as "123+45" or "123-45".
fn format_intronic_pos(cdna_boundary: usize, offset: i64) -> String {
    if offset > 0 {
        format!("{cdna_boundary}+{offset}")
    } else if offset < 0 {
        format!("{cdna_boundary}{offset}")
    } else {
        // Offset of 0 means right at the boundary (shouldn't happen for intronic).
        format!("{cdna_boundary}")
    }
}

/// Find the nearest exon boundary to a genomic position and return
/// (cdna_position_at_boundary, signed_offset_from_boundary).
///
/// Positive offset = downstream of the exon end (donor side, +).
/// Negative offset = upstream of the exon start (acceptor side, -).
fn nearest_exon_boundary(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    sorted_exon_bounds: &[(i64, i64)],
    variant_start: i64,
    _variant_end: i64,
) -> Option<(usize, i64)> {
    let pos = variant_start;
    let strand = tx.strand;

    // Build ordered exon list in transcript direction.
    let mut segments: Vec<(i64, i64)> = sorted_exon_bounds.to_vec();
    if strand < 0 {
        segments.reverse();
    }

    // For each intron gap between consecutive exons, check if position falls in it.
    for i in 0..segments.len().saturating_sub(1) {
        let (prev_start, prev_end) = segments[i];
        let (next_start, next_end) = segments[i + 1];

        let (intron_start, intron_end) = if strand >= 0 {
            (prev_end + 1, next_start - 1)
        } else {
            (next_end + 1, prev_start - 1)
        };

        if pos >= intron_start && pos <= intron_end {
            // Position is in this intron gap.
            if strand >= 0 {
                // Distance from donor (prev exon end) and acceptor (next exon start).
                let dist_to_donor = pos - prev_end;
                let dist_to_acceptor = next_start - pos;

                if dist_to_donor <= dist_to_acceptor {
                    // Closer to donor (upstream exon end).
                    let cdna = genomic_to_cdna_index_pub(tx_exons, strand, prev_end)?;
                    return Some((cdna, dist_to_donor));
                } else {
                    // Closer to acceptor (downstream exon start).
                    let cdna = genomic_to_cdna_index_pub(tx_exons, strand, next_start)?;
                    return Some((cdna, -dist_to_acceptor));
                }
            } else {
                // Negative strand: donor is at next_end, acceptor is at prev_start.
                let dist_to_donor = next_end - pos;
                let dist_to_acceptor = pos - prev_start;

                if dist_to_donor.abs() <= dist_to_acceptor.abs() {
                    // Closer to donor.
                    let cdna = genomic_to_cdna_index_pub(tx_exons, strand, next_end)?;
                    return Some((cdna, dist_to_donor.abs()));
                } else {
                    // Closer to acceptor.
                    let cdna = genomic_to_cdna_index_pub(tx_exons, strand, prev_start)?;
                    return Some((cdna, -(dist_to_acceptor.abs())));
                }
            }
        }
    }

    // Position may be beyond all exons (UTR-intronic or upstream/downstream).
    None
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

    let protein_id = versioned_id(
        translation.stable_id.as_deref()?,
        translation.version,
    );

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

fn format_hgvsp_synonymous(
    protein_id: &str,
    pos_str: &str,
    ref_aa: &str,
) -> Option<String> {
    let pos_num = pos_str.split('-').next()?;
    let ref3 = aa_one_to_three(ref_aa.chars().next()?);
    // VEP URL-encodes the "=" sign as %3D
    Some(format!("{protein_id}:p.{ref3}{pos_num}%3D"))
}

fn format_hgvsp_inframe_del(
    protein_id: &str,
    pos_str: &str,
    ref_aa: &str,
) -> Option<String> {
    let parts: Vec<&str> = pos_str.split('-').collect();
    if parts.len() == 2 {
        let start_pos = parts[0];
        let end_pos = parts[1];
        let first_ref3 = aa_one_to_three(ref_aa.chars().next()?);
        let last_ref3 = aa_one_to_three(ref_aa.chars().last()?);
        Some(format!("{protein_id}:p.{first_ref3}{start_pos}_{last_ref3}{end_pos}del"))
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
        alt_aa[ref_aa.len()..].chars().map(aa_one_to_three).collect()
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
        Some(format!("{protein_id}:p.{ref3}{pos_num}_{ref3}{next_pos}ins{ins_seq}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aa_one_to_three() {
        assert_eq!(aa_one_to_three('A'), "Ala");
        assert_eq!(aa_one_to_three('R'), "Arg");
        assert_eq!(aa_one_to_three('*'), "Ter");
        assert_eq!(aa_one_to_three('X'), "Xaa");
    }

    #[test]
    fn test_versioned_id() {
        assert_eq!(versioned_id("ENST00000379410", Some(6)), "ENST00000379410.6");
        assert_eq!(versioned_id("ENST00000379410", None), "ENST00000379410");
    }

    #[test]
    fn test_format_cdna_change_snv() {
        assert_eq!(format_cdna_change("1043", "G", "A", false, false), "1043G>A");
    }

    #[test]
    fn test_format_cdna_change_insertion() {
        assert_eq!(format_cdna_change("10-11", "-", "ACG", true, false), "10-11insACG");
    }

    #[test]
    fn test_format_cdna_change_deletion() {
        assert_eq!(format_cdna_change("10", "A", "-", false, true), "10del");
        assert_eq!(format_cdna_change("10-12", "ACG", "-", false, true), "10-12del");
    }

    #[test]
    fn test_format_cdna_change_delins() {
        assert_eq!(
            format_cdna_change("10-12", "ACG", "TT", false, false),
            "10-12delinsTT"
        );
    }

    #[test]
    fn test_format_intronic_pos() {
        assert_eq!(format_intronic_pos(100, 5), "100+5");
        assert_eq!(format_intronic_pos(100, -3), "100-3");
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
}
