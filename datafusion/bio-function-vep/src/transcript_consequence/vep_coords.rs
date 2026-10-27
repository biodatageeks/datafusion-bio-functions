//! One port of VEP's transcript coordinates for a (transcript, allele):
//! `cdna_start/cdna_end`, `cds_start/cds_end` and
//! `translation_start/translation_end` (vepyr#240, R2 node F0).
//!
//! VEP reads each pair as the first and last piece of the mapper result,
//! and leaves a coordinate undefined when that piece is a `Gap`
//! (ensembl-variation 116 `BaseTranscriptVariation.pm:136-152` `cdna_start`,
//! `:252-268` `cds_start`, `:365-381` `translation_start`). The pieces come
//! from Ensembl core's `TranscriptMapper` (`c0cf13d`):
//!
//! - `genomic2cdna` is `map_coordinates("genome", ...)`
//!   (`TranscriptMapper.pm:296-307`); `Mapper::map_coordinates` routes a
//!   `start == end + 1` insertion to `map_insert` (`Mapper.pm:267-270`) and
//!   returns pieces in transcript order (`Mapper.pm:447-450`).
//! - `map_insert` (`Mapper.pm:485-580`) maps the swapped two-base span and
//!   keeps only the flanks that are `Coordinate`s, incrementing the first's
//!   start and decrementing the second's end: an insertion is the empty
//!   interval `(a+1, a)`.
//! - `genomic2cds` (`TranscriptMapper.pm:410-482`) gaps whatever lies outside
//!   `[cdna_coding_start, cdna_coding_end]`, so an insertion at either CDS
//!   edge has no CDS coordinate.
//! - `genomic2pep` (`TranscriptMapper.pm:504-535`) is
//!   `int((cds + start_phase + 2) / 3)`.
//!
//! Because the pieces are in transcript order, the first/last-piece rule is
//! the map of the transcript-5'-most and transcript-3'-most genomic base, so
//! [`VepTxMapper::cdna`], [`VepTxMapper::cds`] and
//! [`VepTxMapper::translation`] do two base lookups and allocate nothing.
//! Only the `*_pieces` accessors build piece lists.
//!
//! This module is pure: it does not reproduce VEP's memoisation of shifted
//! versus unshifted coordinates (`TranscriptVariationAllele.pm:539-559`,
//! `OutputFactory.pm:317`); callers choose the span they pass.
//!
//! Most accessors are not wired yet. They are consumed by later fixes
//! (#220, #213+#234, #202, #215, #243, #219), hence the `dead_code` allowance.
#![allow(dead_code)]

use super::{
    ExonFeature, TranscriptFeature, VariantInput, edited_transcript_cdna_index,
    edited_transcript_sequence_cdna_index, genomic_to_cdna_index_for_hgvsp,
    genomic_to_cdna_index_for_transcript, genomic_to_cds_index, mapper_insertion_gap_cdna_position,
    reference_translateable_seq_for_hgvsp, shifted_insertion_cdna_flank,
    uses_refseq_transcript_reference,
};
use super::{TranslationFeature, format_coords_ensembl};

/// Ensembl's `overlap()` (`VariationEffect.pm:72-85`), without range
/// normalisation.
pub(super) use super::overlap_perl as overlap;

/// VEP's VF span, `seq_region_start/end`. An insertion is `(P, P-1)`, i.e.
/// `start == end + 1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct VepSpan {
    pub start: i64,
    pub end: i64,
}

impl VepSpan {
    pub fn new(start: i64, end: i64) -> Self {
        Self { start, end }
    }

    /// The engine keeps an insertion as `start == end` with REF `-`, the
    /// insertion sitting between `start - 1` and `start`; VEP's span for it is
    /// `(start, start - 1)`.
    pub fn from_variant(v: &VariantInput) -> Self {
        if v.ref_allele == "-" {
            Self::new(v.start, v.start - 1)
        } else {
            Self::new(v.start, v.end)
        }
    }

    /// BTV's `+ $shifting_offset`; the caller passes `strand * shift_length`.
    pub fn offset(self, by: i64) -> Self {
        Self::new(self.start + by, self.end + by)
    }

    pub fn is_insertion(self) -> bool {
        self.start == self.end + 1
    }
}

/// VEP's `X_start` / `X_end`. `None` means the first / last mapper piece is
/// a `Gap`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct TvRange {
    pub start: Option<usize>,
    pub end: Option<usize>,
}

impl TvRange {
    pub fn new(start: usize, end: usize) -> Self {
        Self {
            start: Some(start),
            end: Some(end),
        }
    }

    /// `defined $start && defined $end`.
    pub fn both(self) -> Option<(usize, usize)> {
        Some((self.start?, self.end?))
    }

    /// Perl's `$start && $end`: 0 is false, so the insertion before cDNA 1,
    /// `(1, 0)`, reads as undefined (`VariationEffect.pm:981`, `:1381`,
    /// `:1508`).
    pub fn both_truthy(self) -> Option<(usize, usize)> {
        self.both().filter(|&(s, e)| s != 0 && e != 0)
    }

    /// VEP `format_coords` (`*_position` output columns).
    pub fn format(self) -> Option<String> {
        format_coords_ensembl(self.start, self.end)
    }

    /// `return 0 unless $s && $e; overlap($s, $e, $lo, $hi)`.
    pub fn overlap(self, lo: usize, hi: usize) -> bool {
        self.both_truthy()
            .is_some_and(|(s, e)| overlap(s as i64, e as i64, lo as i64, hi as i64))
    }

    fn map(self, f: impl Fn(usize) -> usize) -> Self {
        Self {
            start: self.start.map(&f),
            end: self.end.map(&f),
        }
    }
}

/// One `Bio::EnsEMBL::Mapper` result piece, in transcript order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MapperPiece {
    Coord { start: usize, end: usize },
    Gap,
}

/// Which per-base genomic->cDNA map stands in for VEP's `TranscriptMapper`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CdnaBaseMap {
    /// What `cDNA_position` uses: `genomic_to_cdna_index_for_transcript`
    /// plus the RefSeq output-space adjustment, and the 1-base
    /// insertion-edit override for insertions. Consequence predicates use
    /// this.
    Output,
    /// What shifted HGVSp uses: `genomic_to_cdna_index_for_hgvsp`, and
    /// `shifted_insertion_cdna_flank` for insertion flanks.
    HgvspExact,
}

/// All three VEP coordinate pairs for one span.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct TvCoords {
    pub cdna: TvRange,
    pub cds: TvRange,
    pub translation: TvRange,
}

/// `TranscriptMapper::genomic2pep`'s `int((cds + phase + 2) / 3)`, given a
/// CDS coordinate that already includes the phase (BTV `cds_start`).
pub(crate) fn peptide_of_cds(cds_with_phase: usize) -> usize {
    cds_with_phase.saturating_add(2) / 3
}

/// The RefSeq output-space adjustment `cDNA_position` applies to each
/// mapped base.
pub(crate) fn adjust_output_cdna(tx: &TranscriptFeature, cdna: usize) -> Option<usize> {
    if !uses_refseq_transcript_reference(tx) {
        return Some(cdna);
    }
    if tx.cds_start.is_some() && tx.cds_end.is_some() {
        edited_transcript_cdna_index(tx, cdna)
    } else {
        edited_transcript_sequence_cdna_index(tx, cdna)
    }
}

/// VEP's `TranscriptMapper` for one transcript.
pub(crate) struct VepTxMapper<'a> {
    tx: &'a TranscriptFeature,
    exons: &'a [&'a ExonFeature],
    base: CdnaBaseMap,
    /// Leading `N` count of the translateable sequence. `translateable_seq`
    /// pads `"N" x start_Exon->phase` (`Transcript.pm:919` at `c0cf13d`), so
    /// this is BTV's `start_Exon->phase` (`BaseTranscriptVariation.pm:262`)
    /// and `genomic2pep`'s `start_phase` (`TranscriptMapper.pm:515-518`)
    /// whenever the first exon is the start exon.
    phase: usize,
}

impl<'a> VepTxMapper<'a> {
    /// The phase is read from the translateable sequence when a translation
    /// is given, and is 0 otherwise (it only affects CDS and protein
    /// coordinates).
    pub fn new(
        tx: &'a TranscriptFeature,
        exons: &'a [&'a ExonFeature],
        translation: Option<&TranslationFeature>,
        base: CdnaBaseMap,
    ) -> Self {
        let phase = translation
            .and_then(|translation| reference_translateable_seq_for_hgvsp(tx, Some(translation)))
            .map(|seq| seq.bytes().take_while(|&b| b == b'N' || b == b'n').count())
            .unwrap_or(0);
        Self::with_phase(tx, exons, base, phase)
    }

    /// For callers that already hold the leading-`N` count.
    pub fn with_phase(
        tx: &'a TranscriptFeature,
        exons: &'a [&'a ExonFeature],
        base: CdnaBaseMap,
        phase: usize,
    ) -> Self {
        Self {
            tx,
            exons,
            base,
            phase,
        }
    }

    fn base_cdna(&self, pos: i64, insertion_flank: bool) -> Option<usize> {
        match self.base {
            CdnaBaseMap::Output => genomic_to_cdna_index_for_transcript(self.tx, self.exons, pos)
                .and_then(|cdna| adjust_output_cdna(self.tx, cdna)),
            CdnaBaseMap::HgvspExact if insertion_flank => {
                shifted_insertion_cdna_flank(self.tx, self.exons, pos)
            }
            CdnaBaseMap::HgvspExact => genomic_to_cdna_index_for_hgvsp(self.tx, self.exons, pos),
        }
    }

    /// The (transcript-5', transcript-3') genomic ends of `span`. For an
    /// insertion these are its two flanks.
    fn transcript_ends(&self, span: VepSpan) -> (i64, i64) {
        let (low, high) = if span.is_insertion() {
            (span.end, span.start)
        } else {
            (span.start, span.end)
        };
        if self.tx.strand < 0 {
            (high, low)
        } else {
            (low, high)
        }
    }

    /// The insertion flanks' cDNA, `(5' flank, 3' flank)`.
    fn insertion_flanks(&self, span: VepSpan) -> (Option<usize>, Option<usize>) {
        let (five, three) = self.transcript_ends(span);
        (self.base_cdna(five, true), self.base_cdna(three, true))
    }

    /// `genomic2cdna` + BTV `cdna_start` / `cdna_end`.
    pub fn cdna(&self, span: VepSpan) -> TvRange {
        if !span.is_insertion() {
            let (five, three) = self.transcript_ends(span);
            return TvRange {
                start: self.base_cdna(five, false),
                end: self.base_cdna(three, false),
            };
        }
        if self.base == CdnaBaseMap::Output
            && let Some(cdna) = mapper_insertion_gap_cdna_position(self.tx, span.start)
        {
            // Flanks straddling a 1-base insertion-edit cDNA gap: `map_insert`
            // gives `(a+1, b-1)` with `b == a+2`.
            return TvRange::new(cdna, cdna);
        }
        match self.insertion_flanks(span) {
            // If output-space adjustment collapses the two genomic flanks to
            // one cDNA position (a deleted RefSeq edit), keep the established
            // `a-a` display rather than expanding it around the edit.
            (Some(five), Some(three)) if self.base == CdnaBaseMap::Output && five == three => {
                TvRange::new(three, five)
            }
            // Without transcript-mapper segments, RefSeq output-space offsets
            // can collapse a deleted edit between the genomic flanks to the
            // same cDNA base. Preserve that legacy `a-a` display; core's
            // map_insert adjustment applies to mapper-produced coordinates.
            (Some(five), Some(three))
                if self.base == CdnaBaseMap::Output && self.tx.cdna_mapper_segments.is_empty() =>
            {
                TvRange::new(three, five)
            }
            // `map_insert` increments the transcript-5' flank's start and
            // decrements the transcript-3' flank's end. For an edit-inserted
            // cDNA gap, this keeps only the intervening mapped coordinates.
            (Some(five), Some(three)) => {
                TvRange::new(five.saturating_add(1), three.saturating_sub(1))
            }
            // Only the 3' flank maps: `map_insert` decrements its end.
            (None, Some(three)) => TvRange::new(three, three.saturating_sub(1)),
            // Only the 5' flank maps: `map_insert` increments its start.
            (Some(five), None) => TvRange::new(five.saturating_add(1), five),
            (None, None) => TvRange::default(),
        }
    }

    /// The `cDNA_position` output column: VEP `format_coords` of
    /// [`Self::cdna`].
    pub fn cdna_position(&self, span: VepSpan) -> Option<String> {
        let range = self.cdna(span);
        if span.is_insertion()
            && self.base == CdnaBaseMap::Output
            && let Some((start, end)) = range.both()
            && start == end
            && mapper_insertion_gap_cdna_position(self.tx, span.start).is_none()
        {
            // ENGINE != VEP: both flanks reach the same cDNA base after the
            // RefSeq output-space adjustment (a deletion SeqEdit between
            // them). cDNA_position has always printed `a-a` here, which
            // `format_coords` cannot spell.
            return Some(format!("{start}-{end}"));
        }
        range.format()
    }

    /// `genomic2cds` + BTV `cds_start` / `cds_end` over an already-mapped
    /// cDNA range, adding the start-exon phase.
    ///
    /// For an insertion the gap test is applied to the `(c, c-1)` pair, which
    /// is `genomic2cds`'s single-piece rule (`TranscriptMapper.pm:441-443`):
    /// an insertion at either CDS edge has no CDS coordinate. Otherwise each
    /// end is defined only when it lies in the CDS, since a range reaching
    /// into a UTR starts or ends with a `Gap` piece (`:450-462`).
    pub fn cds_of_cdna(&self, cdna: TvRange, insertion: bool) -> TvRange {
        let (Some(coding_start), Some(coding_end)) =
            (self.tx.cdna_coding_start, self.tx.cdna_coding_end)
        else {
            return TvRange::default();
        };
        let to_cds = |c: usize| {
            c.saturating_sub(coding_start)
                .saturating_add(1)
                .saturating_add(self.phase)
        };
        if insertion {
            let Some((start, end)) = cdna.both() else {
                return TvRange::default();
            };
            if end < coding_start || start > coding_end {
                return TvRange::default();
            }
            // ENGINE != VEP: a non-empty window between two insertion flanks
            // (an insertion-edit cDNA gap) is two mapper pieces in VEP, each
            // gapped on its own; the pair rule keeps the shifted HGVSp
            // behaviour it replaced.
            return TvRange::new(to_cds(start), to_cds(end));
        }
        let in_cds = |c: usize| (coding_start..=coding_end).contains(&c);
        TvRange {
            start: cdna.start.filter(|&c| in_cds(c)).map(to_cds),
            end: cdna.end.filter(|&c| in_cds(c)).map(to_cds),
        }
    }

    /// `genomic2cds` + BTV `cds_start` / `cds_end`.
    pub fn cds(&self, span: VepSpan) -> TvRange {
        self.cds_of_cdna(self.cdna(span), span.is_insertion())
    }

    /// `genomic2pep` + BTV `translation_start` / `translation_end`.
    pub fn translation(&self, span: VepSpan) -> TvRange {
        self.cds(span).map(peptide_of_cds)
    }

    /// The translation start used by partial-codon predicates, with a fallback
    /// for older caches that lack cDNA CDS bounds. Such caches still carry
    /// genomic CDS coordinates, so ordinary mapped bases can use the existing
    /// genomic-to-CDS mapper and the same peptide-coordinate conversion.
    pub fn translation_with_genomic_cds_fallback(&self, span: VepSpan) -> TvRange {
        let mapped = self.translation(span);
        if mapped.start.is_some()
            || (self.tx.cdna_coding_start.is_some() && self.tx.cdna_coding_end.is_some())
        {
            return mapped;
        }

        let (five_prime, three_prime) = self.transcript_ends(span);
        let cds_index_for_fallback = |position| {
            if super::prefers_exon_geometry_over_mapper(self.tx) {
                super::genomic_to_cds_index_from_segments(self.tx, self.exons, position)
            } else {
                genomic_to_cds_index(self.tx, self.exons, position)
            }
        };
        let Some(cds_index) = cds_index_for_fallback(five_prime) else {
            return mapped;
        };
        // VEP's insertion mapper can produce a CDS coordinate only when both
        // insertion flanks map. In particular, an insertion at a CDS edge
        // must stay undefined rather than borrowing its single coding flank.
        if span.is_insertion() && cds_index_for_fallback(three_prime).is_none() {
            return mapped;
        }

        // For insertions `translation_start` is the coordinate immediately
        // after the transcript-5' flank; ordinary variants start at the flank.
        let after_five_prime_flank = usize::from(span.is_insertion()) + 1;
        let cds_with_phase = cds_index
            .saturating_add(self.phase)
            .saturating_add(after_five_prime_flank);
        TvRange {
            start: Some(peptide_of_cds(cds_with_phase)),
            end: mapped.end,
        }
    }

    /// All three pairs, mapping the span once.
    pub fn coords(&self, span: VepSpan) -> TvCoords {
        let cdna = self.cdna(span);
        let cds = self.cds_of_cdna(cdna, span.is_insertion());
        TvCoords {
            cdna,
            cds,
            translation: cds.map(peptide_of_cds),
        }
    }

    /// BTV `cds_coords`' `_boundary_shift` fallback
    /// (`BaseTranscriptVariation.pm:518-526`): when the shifted span starts
    /// in a gap and the VF was shifted, use the unshifted span if it does not.
    /// Returns the range and whether the fallback was taken.
    pub fn cds_shifted(
        &self,
        shifted: VepSpan,
        unshifted: VepSpan,
        shifted_flag: bool,
    ) -> (TvRange, bool) {
        let range = self.cds(shifted);
        if range.start.is_none() && shifted_flag {
            let fallback = self.cds(unshifted);
            if fallback.start.is_some() {
                return (fallback, true);
            }
        }
        (range, false)
    }

    /// BTV `translation_coords`' `_boundary_shift` fallback
    /// (`BaseTranscriptVariation.pm:548-556`), VEP's
    /// `translation_start(undef, $shifting_offset)`.
    pub fn translation_shifted(
        &self,
        shifted: VepSpan,
        unshifted: VepSpan,
        shifted_flag: bool,
    ) -> (TvRange, bool) {
        let range = self.translation(shifted);
        if range.start.is_none() && shifted_flag {
            let fallback = self.translation(unshifted);
            if fallback.start.is_some() {
                return (fallback, true);
            }
        }
        (range, false)
    }

    /// `genomic2cdna`'s pieces, for `within_cdna` / `within_cds`-style
    /// predicates. Not for the hot path: it walks the span base by base.
    pub fn cdna_pieces(&self, span: VepSpan) -> Vec<MapperPiece> {
        if span.is_insertion() {
            let (five, three) = self.insertion_flanks(span);
            let first = five.map(|c| MapperPiece::Coord {
                start: c.saturating_add(1),
                end: c,
            });
            let second = three.map(|c| MapperPiece::Coord {
                start: c,
                end: c.saturating_sub(1),
            });
            return match (first, second) {
                // Both flanks in one mapper pair: one swapped piece.
                (Some(first), Some(_)) if three == five.map(|c| c + 1) => vec![first],
                _ => first.into_iter().chain(second).collect(),
            };
        }
        let (five, three) = self.transcript_ends(span);
        let step = if five <= three { 1 } else { -1 };
        let mut pieces: Vec<MapperPiece> = Vec::new();
        let mut pos = five;
        loop {
            let piece = match self.base_cdna(pos, false) {
                Some(c) => MapperPiece::Coord { start: c, end: c },
                None => MapperPiece::Gap,
            };
            match (pieces.last_mut(), piece) {
                (Some(MapperPiece::Gap), MapperPiece::Gap) => {}
                (Some(MapperPiece::Coord { end, .. }), MapperPiece::Coord { start: c, .. })
                    if c == end.saturating_add(1) =>
                {
                    *end = c;
                }
                _ => pieces.push(piece),
            }
            if pos == three {
                break;
            }
            pos += step;
        }
        pieces
    }

    /// `genomic2cds`'s pieces (`TranscriptMapper.pm:430-475`), with the
    /// start-exon phase added to each coordinate.
    pub fn cds_pieces(&self, span: VepSpan) -> Vec<MapperPiece> {
        let (Some(coding_start), Some(coding_end)) =
            (self.tx.cdna_coding_start, self.tx.cdna_coding_end)
        else {
            return vec![MapperPiece::Gap];
        };
        let to_cds = |c: usize| c - coding_start + 1 + self.phase;
        let mut out = Vec::new();
        for piece in self.cdna_pieces(span) {
            let MapperPiece::Coord { start, end } = piece else {
                out.push(MapperPiece::Gap);
                continue;
            };
            if end < coding_start || start > coding_end {
                out.push(MapperPiece::Gap);
                continue;
            }
            if start < coding_start {
                out.push(MapperPiece::Gap);
            }
            out.push(MapperPiece::Coord {
                start: to_cds(start.max(coding_start)),
                end: to_cds(end.min(coding_end)),
            });
            if end > coding_end {
                out.push(MapperPiece::Gap);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::super::{
        CodonTable, GeometryCache, ProteinDomainFeature, RefSeqEdit, ShiftedTvaCoords,
        TranscriptCdnaMapperSegment, adjust_refseq_cds_sequence_index_for_hgvsp_reference,
        genomic_range, genomic_to_cds_index_for_hgvsp, overlaps, refseq_has_edited_sequence_state,
        shifted_tva_coords_from_mapper,
    };
    use super::*;

    // ------------------------------------------------------------------
    // Fixtures
    // ------------------------------------------------------------------

    fn tx_feature(
        id: &str,
        strand: i8,
        exons: &[(i64, i64)],
        coding: Option<(usize, usize)>,
    ) -> TranscriptFeature {
        let start = exons.iter().map(|e| e.0).min().unwrap();
        let end = exons.iter().map(|e| e.1).max().unwrap();
        TranscriptFeature {
            transcript_id: id.to_string(),
            transcript_uid: None,
            db_id: None,
            codon_table: CodonTable::Standard,
            chrom: "1".to_string(),
            start,
            end,
            strand,
            biotype: "protein_coding".to_string(),
            cds_start: coding.map(|_| start),
            cds_end: coding.map(|_| end),
            cdna_coding_start: coding.map(|c| c.0),
            cdna_coding_end: coding.map(|c| c.1),
            cdna_mapper_segments: Vec::new(),
            mature_mirna_regions: Vec::new(),
            gene_stable_id: None,
            gene_symbol: None,
            gene_symbol_source: None,
            gene_hgnc_id_native: None,
            gene_hgnc_id: None,
            display_xref_id: None,
            source: None,
            source_cache: None,
            refseq_match: None,
            refseq_edits: Vec::new(),
            translation_seq_edits: Vec::new(),
            has_rna_edit: false,
            is_gencode_basic: false,
            is_gencode_primary: false,
            bam_edit_status: None,
            has_non_polya_rna_edit: false,
            spliced_seq: None,
            five_prime_utr_seq: None,
            three_prime_utr_seq: None,
            translateable_seq: None,
            cdna_seq: None,
            version: None,
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
            cdna_coords_cache: GeometryCache::default(),
        }
    }

    fn exon_features(id: &str, exons: &[(i64, i64)]) -> Vec<ExonFeature> {
        exons
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| ExonFeature {
                transcript_id: id.to_string(),
                exon_number: i as i32 + 1,
                start,
                end,
            })
            .collect()
    }

    fn variant(start: i64, end: i64, ref_allele: &str, alt_allele: &str) -> VariantInput {
        VariantInput {
            chrom: "1".to_string(),
            start,
            end,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            parser_start: start,
            parser_end: end,
            parser_ref_allele: ref_allele.to_string(),
            parser_alt_allele: alt_allele.to_string(),
            hgvs_shift_forward: None,
            hgvs_shift_reverse: None,
        }
    }

    #[derive(Clone, Copy, Debug)]
    enum MapperShape {
        None,
        Contiguous,
        /// A length-changing insertion SeqEdit: a cDNA gap of this many bases
        /// between two genomic-contiguous pairs (`TranscriptMapper.pm:207-211`).
        CdnaGap(usize),
        /// A deletion SeqEdit: one genomic base with no cDNA.
        GenomicGap,
    }

    /// Mapper pairs in genomic order, cDNA numbered in transcript order,
    /// with the first genomic exon split for the gap shapes.
    fn mapper_segments(
        strand: i8,
        exons: &[(i64, i64)],
        shape: MapperShape,
    ) -> Vec<TranscriptCdnaMapperSegment> {
        // (genomic_start, genomic_end, extra cDNA gap before this piece in
        // transcript order)
        let mut pieces: Vec<(i64, i64, bool)> = Vec::new();
        for (i, &(s, e)) in exons.iter().enumerate() {
            let mid = (s + e) / 2;
            match shape {
                MapperShape::None => return Vec::new(),
                MapperShape::CdnaGap(_) if i == 0 => {
                    pieces.push((s, mid, strand < 0));
                    pieces.push((mid + 1, e, strand >= 0));
                }
                MapperShape::GenomicGap if i == 0 => {
                    pieces.push((s, mid - 1, false));
                    pieces.push((mid + 1, e, false));
                }
                _ => pieces.push((s, e, false)),
            }
        }
        let gap = match shape {
            MapperShape::CdnaGap(k) => k,
            _ => 0,
        };
        let mut ordered = pieces.clone();
        if strand < 0 {
            ordered.reverse();
        }
        let mut cursor = 1usize;
        let mut out = Vec::new();
        for (s, e, gap_before) in ordered {
            if gap_before {
                cursor += gap;
            }
            let len = (e - s + 1) as usize;
            out.push(TranscriptCdnaMapperSegment {
                genomic_start: s,
                genomic_end: e,
                cdna_start: cursor,
                cdna_end: cursor + len - 1,
                ori: if strand < 0 { -1 } else { 1 },
            });
            cursor += len;
        }
        out.sort_by_key(|seg| (seg.genomic_start, seg.genomic_end, seg.cdna_start));
        out
    }

    #[derive(Clone, Copy, Debug)]
    enum Flavour {
        Ensembl,
        RefSeqPlain,
        RefSeqDeletionEdit,
        RefSeqInsertionEdit,
    }

    /// The synthetic matrix of §4.4: strand, exon layout, UTR, mapper shape
    /// and transcript flavour.
    fn matrix() -> Vec<(TranscriptFeature, Vec<ExonFeature>)> {
        let layouts: [&[(i64, i64)]; 2] = [
            &[(1_000, 1_059)],
            // A 2-bp (frameshift) intron, then a 12-bp intron.
            &[(1_000, 1_029), (1_032, 1_061), (1_074, 1_103)],
        ];
        let shapes = [
            MapperShape::None,
            MapperShape::Contiguous,
            MapperShape::CdnaGap(1),
            MapperShape::CdnaGap(2),
            MapperShape::GenomicGap,
        ];
        let flavours = [
            Flavour::Ensembl,
            Flavour::RefSeqPlain,
            Flavour::RefSeqDeletionEdit,
            Flavour::RefSeqInsertionEdit,
        ];
        let mut out = Vec::new();
        for strand in [1i8, -1] {
            for layout in layouts {
                let total: usize = layout.iter().map(|e| (e.1 - e.0 + 1) as usize).sum();
                for coding in [None, Some((1, total)), Some((5, total - 6))] {
                    for shape in shapes {
                        for flavour in flavours {
                            let id = match flavour {
                                Flavour::Ensembl => "ENST00000000001.1",
                                _ => "NM_000001.1",
                            };
                            let mut tx = tx_feature(id, strand, layout, coding);
                            tx.cdna_mapper_segments = mapper_segments(strand, layout, shape);
                            match flavour {
                                Flavour::Ensembl => tx.source = Some("Ensembl".to_string()),
                                Flavour::RefSeqPlain => tx.source = Some("RefSeq".to_string()),
                                Flavour::RefSeqDeletionEdit => {
                                    tx.source = Some("RefSeq".to_string());
                                    tx.has_rna_edit = true;
                                    tx.refseq_edits = vec![RefSeqEdit {
                                        start: 12,
                                        end: 12,
                                        replacement_len: None,
                                        skip_refseq_offset: false,
                                    }];
                                }
                                Flavour::RefSeqInsertionEdit => {
                                    tx.source = Some("RefSeq".to_string());
                                    tx.has_rna_edit = true;
                                    tx.has_non_polya_rna_edit = true;
                                    tx.refseq_edits = vec![RefSeqEdit {
                                        start: 21,
                                        end: 20,
                                        replacement_len: Some(2),
                                        skip_refseq_offset: false,
                                    }];
                                }
                            }
                            out.push((tx, exon_features(id, layout)));
                        }
                    }
                }
            }
        }
        out
    }

    /// Every start in `[tx_start - 3, tx_end + 3]`: an SNV, deletions of
    /// length 1..=6 and an insertion.
    fn variants_for(tx: &TranscriptFeature) -> Vec<VariantInput> {
        let mut out = Vec::new();
        for p in tx.start - 3..=tx.end + 3 {
            out.push(variant(p, p, "A", "G"));
            for len in 1..=6i64 {
                out.push(variant(p, p + len - 1, &"A".repeat(len as usize), "-"));
            }
            out.push(variant(p, p, "-", "T"));
        }
        out
    }

    // ------------------------------------------------------------------
    // Oracle: the code F0 replaced, kept for the life of the PR.
    // ------------------------------------------------------------------

    /// `compute_cdna_position` at functions `9936c6e`, verbatim.
    fn legacy_compute_cdna_position(
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
    ) -> Option<String> {
        if tx_exons.is_empty() {
            return None;
        }
        let is_ins = variant.ref_allele == "-";
        let in_exon = tx_exons.iter().any(|e| {
            if is_ins {
                let left = variant.start - 1;
                let right = variant.start;
                (left >= e.start && left <= e.end) || (right >= e.start && right <= e.end)
            } else {
                overlaps(variant.start, variant.end, e.start, e.end)
            }
        });
        if !in_exon {
            return None;
        }
        let adjust_output_cdna = |cdna: usize| {
            if !uses_refseq_transcript_reference(tx) {
                return Some(cdna);
            }
            if tx.cds_start.is_some() && tx.cds_end.is_some() {
                edited_transcript_cdna_index(tx, cdna)
            } else {
                edited_transcript_sequence_cdna_index(tx, cdna)
            }
        };
        if is_ins {
            if let Some(boundary_cdna) = mapper_insertion_gap_cdna_position(tx, variant.start) {
                return Some(boundary_cdna.to_string());
            }
            let a =
                genomic_to_cdna_index_for_transcript(tx, tx_exons, variant.start.saturating_sub(1))
                    .and_then(adjust_output_cdna);
            let b = genomic_to_cdna_index_for_transcript(tx, tx_exons, variant.start)
                .and_then(adjust_output_cdna);
            match (a, b) {
                (Some(a), Some(b)) => {
                    let lo = a.min(b);
                    let hi = a.max(b);
                    Some(format!("{lo}-{hi}"))
                }
                (None, Some(b)) => {
                    let other = if tx.strand >= 0 {
                        b.saturating_sub(1)
                    } else {
                        b + 1
                    };
                    let lo = b.min(other);
                    let hi = b.max(other);
                    Some(format!("{lo}-{hi}"))
                }
                (Some(a), None) => {
                    let other = if tx.strand >= 0 {
                        a + 1
                    } else {
                        a.saturating_sub(1)
                    };
                    let lo = a.min(other);
                    let hi = a.max(other);
                    Some(format!("{lo}-{hi}"))
                }
                (None, None) => None,
            }
        } else {
            let start_cdna = genomic_to_cdna_index_for_transcript(tx, tx_exons, variant.start)
                .and_then(adjust_output_cdna);
            let end_cdna = genomic_to_cdna_index_for_transcript(tx, tx_exons, variant.end)
                .and_then(adjust_output_cdna);
            match (start_cdna, end_cdna) {
                (Some(s), Some(e)) if s == e => Some(s.to_string()),
                (Some(s), Some(e)) => {
                    let lo = s.min(e);
                    let hi = s.max(e);
                    Some(format!("{lo}-{hi}"))
                }
                (Some(s), None) => {
                    if tx.strand < 0 {
                        Some(format!("?-{s}"))
                    } else {
                        Some(format!("{s}-?"))
                    }
                }
                (None, Some(e)) => {
                    if tx.strand < 0 {
                        Some(format!("{e}-?"))
                    } else {
                        Some(format!("?-{e}"))
                    }
                }
                (None, None) => None,
            }
        }
    }

    #[test]
    fn f0_cdna_matches_legacy_compute_cdna_position() {
        let mut compared = 0usize;
        let mut defined = 0usize;
        for (tx, exons) in matrix() {
            let refs: Vec<&ExonFeature> = exons.iter().collect();
            for v in variants_for(&tx) {
                let legacy = legacy_compute_cdna_position(&v, &tx, &refs);
                let new = super::super::compute_cdna_position(&v, &tx, &refs);
                let expected = if new != legacy
                    && v.ref_allele == "-"
                    && !tx.cdna_mapper_segments.is_empty()
                {
                    let mapper = VepTxMapper::new(&tx, &refs, None, CdnaBaseMap::Output);
                    let (five, three) = mapper.insertion_flanks(VepSpan::from_variant(&v));
                    match (five, three) {
                        (Some(five), Some(three))
                            if super::super::refseq_inserted_cdna_gap_between(&tx, five, three) =>
                        {
                            TvRange::new(five.saturating_add(1), three.saturating_sub(1)).format()
                        }
                        _ => legacy.clone(),
                    }
                } else {
                    legacy.clone()
                };
                assert_eq!(
                    new,
                    expected,
                    "cDNA_position differs: tx={} strand={} segments={:?} edits={:?} variant={}-{} {}>{}",
                    tx.transcript_id,
                    tx.strand,
                    tx.cdna_mapper_segments,
                    tx.refseq_edits,
                    v.start,
                    v.end,
                    v.ref_allele,
                    v.alt_allele
                );
                compared += 1;
                defined += usize::from(legacy.is_some());
            }
        }
        // 240 transcripts; guard against a matrix that silently shrank.
        assert!(compared > 150_000, "compared {compared}");
        assert!(defined > 50_000, "defined {defined}");
    }

    /// `shifted_tva_coords_from_mapper` at functions `9936c6e`, verbatim.
    fn legacy_shifted_tva_coords_from_mapper(
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
        translation: &TranslationFeature,
        shifted_variant: &VariantInput,
    ) -> Option<ShiftedTvaCoords> {
        let leading_n_offset = reference_translateable_seq_for_hgvsp(tx, Some(translation))?
            .bytes()
            .take_while(|&b| b == b'N' || b == b'n')
            .count();
        let coding_start = tx.cdna_coding_start?;
        let coding_end = tx.cdna_coding_end?;

        let is_insertion = shifted_variant.start == shifted_variant.end.saturating_add(1);
        let (cdna_start, cdna_end) = if is_insertion {
            // Traceability:
            // - Ensembl core `Mapper::map_insert()` maps the swapped two-base
            //   boundary and then adjusts the mapped interval in transcript
            //   coordinates before downstream `genomic2cds()` / `genomic2pep()`
            //   consumers read start/end.
            //   <https://github.com/Ensembl/ensembl/blob/release/115/modules/Bio/EnsEMBL/Mapper.pm#L485-L574>
            //
            // Use flanks in transcript order, not genomic order. On the reverse
            // strand the higher genomic coordinate is upstream. This swaps a
            // contiguous mapped interval into a zero-length insertion while
            // retaining any edited-cDNA gap between separate mapped segments.
            // Do not sort cDNA positions: VEP preserves mapper result order.
            //
            // `map_insert` keeps only the flanks that map to a
            // `Bio::EnsEMBL::Mapper::Coordinate` and silently drops any that map
            // to a `Gap`; it does not fail when one flank is unmapped. An
            // insertion abutting a frameshift intron has one intronic flank, so
            // requiring both to map would discard a window VEP still resolves.
            let (upstream, downstream) = if tx.strand < 0 {
                (shifted_variant.start, shifted_variant.end)
            } else {
                (shifted_variant.end, shifted_variant.start)
            };
            let left = shifted_insertion_cdna_flank(tx, tx_exons, upstream);
            let right = shifted_insertion_cdna_flank(tx, tx_exons, downstream);
            match (left, right) {
                (Some(left), Some(right)) => (left.saturating_add(1), right.saturating_sub(1)),
                // Only the downstream flank mapped: VEP keeps that coordinate and
                // decrements its end, leaving the zero-length window (r, r-1).
                (None, Some(right)) => (right, right.saturating_sub(1)),
                // Only the upstream flank mapped: VEP increments its start.
                (Some(left), None) => (left.saturating_add(1), left),
                (None, None) => return None,
            }
        } else {
            let genomic_positions = genomic_range(shifted_variant.start, shifted_variant.end)?;
            let mut cdna_positions = Vec::with_capacity(genomic_positions.len());
            for pos in genomic_positions {
                cdna_positions.push(genomic_to_cdna_index_for_hgvsp(tx, tx_exons, pos)?);
            }
            cdna_positions.sort_unstable();
            (*cdna_positions.first()?, *cdna_positions.last()?)
        };

        if cdna_end < coding_start || cdna_start > coding_end {
            return None;
        }

        let raw_cds_start = cdna_start
            .saturating_sub(coding_start)
            .saturating_add(1)
            .saturating_add(leading_n_offset);
        let raw_cds_end = cdna_end
            .saturating_sub(coding_start)
            .saturating_add(1)
            .saturating_add(leading_n_offset);
        let edited_insertion_gap =
            is_insertion && cdna_start <= cdna_end && refseq_has_edited_sequence_state(tx);
        let cds_start = if edited_insertion_gap {
            raw_cds_start
        } else {
            adjust_refseq_cds_sequence_index_for_hgvsp_reference(
                tx,
                Some(translation),
                raw_cds_start.checked_sub(1)?,
                leading_n_offset,
                false,
            )
            .and_then(|idx| idx.checked_add(1))
            .unwrap_or(raw_cds_start)
        };
        let cds_end = if edited_insertion_gap {
            raw_cds_end
        } else {
            adjust_refseq_cds_sequence_index_for_hgvsp_reference(
                tx,
                Some(translation),
                raw_cds_end.checked_sub(1)?,
                leading_n_offset,
                false,
            )
            .and_then(|idx| idx.checked_add(1))
            .unwrap_or(raw_cds_end)
        };
        let translateable_pos_1based = |genomic_pos: i64| -> Option<usize> {
            let raw_idx = genomic_to_cds_index_for_hgvsp(tx, tx_exons, genomic_pos)?
                .checked_add(leading_n_offset)?;
            adjust_refseq_cds_sequence_index_for_hgvsp_reference(
                tx,
                Some(translation),
                raw_idx,
                leading_n_offset,
                false,
            )
            .and_then(|idx| idx.checked_add(1))
        };

        let (protein_start, protein_end) = if is_insertion {
            // TranscriptMapper::genomic2pep() converts the already-mapped CDS
            // boundaries with int((bound + phase + 2) / 3). The CDS coordinates
            // above already include leading-N padding and reference adjustment.
            // Adding another base to the end incorrectly widens codon-boundary
            // insertions and single-flank windows.
            // https://github.com/Ensembl/ensembl/blob/release/116/modules/Bio/EnsEMBL/TranscriptMapper.pm#L516-L529
            (
                cds_start.saturating_add(2) / 3,
                cds_end.saturating_add(2) / 3,
            )
        } else {
            let genomic_positions = genomic_range(shifted_variant.start, shifted_variant.end)?;
            let mut peptide_positions = Vec::with_capacity(genomic_positions.len());
            for pos in genomic_positions {
                let tr_pos = translateable_pos_1based(pos)?;
                peptide_positions.push(tr_pos.saturating_add(2) / 3);
            }
            peptide_positions.sort_unstable();
            (*peptide_positions.first()?, *peptide_positions.last()?)
        };

        Some(ShiftedTvaCoords {
            cds_start,
            cds_end,
            protein_start,
            protein_end,
            edited_insertion_gap,
        })
    }

    fn translation_with_phase(id: &str, coding_len: usize, phase: usize) -> TranslationFeature {
        let body: String = "ATGC".chars().cycle().take(coding_len).collect();
        let cds = format!("{}{}", "N".repeat(phase), body);
        TranslationFeature {
            transcript_id: id.to_string(),
            cds_len: Some(cds.len()),
            protein_len: Some(cds.len() / 3),
            translation_seq: None,
            cds_sequence: Some(cds.clone()),
            translation_seq_canonical: None,
            cds_sequence_canonical: Some(cds),
            stable_id: None,
            version: None,
            protein_features: Vec::<ProteinDomainFeature>::new(),
        }
    }

    /// Commit 2 of F0: the insertion branch of `shifted_tva_coords_from_mapper`
    /// goes through `VepTxMapper` (HgvspExact). Every shifted insertion span
    /// `(P, P-1)` and every deletion of length 1..=6 over the matrix, with a
    /// leading-N phase of 0, 1 and 2, must give the same coordinates.
    #[test]
    fn f0_shifted_insertion_matches_legacy() {
        let mut compared = 0usize;
        let mut defined = 0usize;
        for (tx, exons) in matrix() {
            let refs: Vec<&ExonFeature> = exons.iter().collect();
            let coding_len = match (tx.cdna_coding_start, tx.cdna_coding_end) {
                (Some(s), Some(e)) => e - s + 1,
                _ => 60,
            };
            for phase in 0..=2 {
                let tl = translation_with_phase(&tx.transcript_id, coding_len, phase);
                for p in tx.start - 3..=tx.end + 3 {
                    let mut spans = vec![variant(p, p - 1, "-", "C")];
                    for len in 1..=6i64 {
                        spans.push(variant(p, p + len - 1, &"A".repeat(len as usize), "-"));
                    }
                    for v in spans {
                        let legacy = legacy_shifted_tva_coords_from_mapper(&tx, &refs, &tl, &v);
                        let new = shifted_tva_coords_from_mapper(&tx, &refs, &tl, &v);
                        assert_eq!(
                            new,
                            legacy,
                            "shifted coords differ: tx={} strand={} coding={:?}-{:?} segments={:?} edits={:?} phase={phase} variant={}-{} {}>{}",
                            tx.transcript_id,
                            tx.strand,
                            tx.cdna_coding_start,
                            tx.cdna_coding_end,
                            tx.cdna_mapper_segments,
                            tx.refseq_edits,
                            v.start,
                            v.end,
                            v.ref_allele,
                            v.alt_allele
                        );
                        compared += 1;
                        defined += usize::from(legacy.is_some());
                    }
                }
            }
        }
        assert!(compared > 150_000, "compared {compared}");
        assert!(defined > 50_000, "defined {defined}");
    }

    // ------------------------------------------------------------------
    // VEP semantics of the accessors (expected values computed from the
    // quoted Ensembl rules, not from the engine).
    // ------------------------------------------------------------------

    /// + strand, exons 1000-1029 / 1032-1061 (2-bp intron) / 1074-1103,
    /// UTR 4 bases each side of the CDS: cDNA 1..90, CDS = cDNA 5..84.
    fn coding_tx(strand: i8) -> (TranscriptFeature, Vec<ExonFeature>) {
        let layout = [(1_000, 1_029), (1_032, 1_061), (1_074, 1_103)];
        let tx = tx_feature("ENST00000000002.1", strand, &layout, Some((5, 84)));
        (tx, exon_features("ENST00000000002.1", &layout))
    }

    fn mapper<'a>(
        tx: &'a TranscriptFeature,
        refs: &'a [&'a ExonFeature],
        phase: usize,
    ) -> VepTxMapper<'a> {
        VepTxMapper::with_phase(tx, refs, CdnaBaseMap::Output, phase)
    }

    /// An insertion is VEP's empty interval: `from_vcf(100, A, AT)` is the
    /// span (101, 100), and an interior insertion maps to `(a+1, a)` on both
    /// strands (`Mapper.pm:488-547`).
    #[test]
    fn insertion_is_empty_interval() {
        let v = VariantInput::from_vcf("1".into(), 100, 100, "A".into(), "AT".into());
        let span = VepSpan::from_variant(&v);
        assert_eq!(span, VepSpan::new(101, 100));
        assert!(span.is_insertion());

        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // Between 1009 (cDNA 10) and 1010 (cDNA 11).
        assert_eq!(
            mapper(&tx, &refs, 0).cdna(VepSpan::new(1_010, 1_009)),
            TvRange::new(11, 10)
        );
        let (tx, exons) = coding_tx(-1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // − strand: 1103 is cDNA 1, so 1094 is cDNA 10 and 1093 cDNA 11.
        assert_eq!(
            mapper(&tx, &refs, 0).cdna(VepSpan::new(1_094, 1_093)),
            TvRange::new(11, 10)
        );
    }

    /// `genomic2pep` (`TranscriptMapper.pm:528-529`): an insertion between
    /// codons, CDS `(k+1, k)` with `3 | k`, is `(p+1, p)`.
    #[test]
    fn insertion_between_codons_gives_empty_peptide_interval() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // CDS 6 = cDNA 10 = 1009; insertion between 1009 and 1010: CDS (7, 6).
        let span = VepSpan::new(1_010, 1_009);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cds(span), TvRange::new(7, 6));
        assert_eq!(m.translation(span), TvRange::new(3, 2));
    }

    #[test]
    fn insertion_inside_codon_gives_single_residue() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // CDS (5, 4): both in codon 2.
        let span = VepSpan::new(1_008, 1_007);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cds(span), TvRange::new(5, 4));
        assert_eq!(m.translation(span), TvRange::new(2, 2));
    }

    /// `genomic2cds`: an insertion `(c, c-1)` with `c == cdna_coding_start`
    /// has `end < cdna_cstart` and is all gap (`TranscriptMapper.pm:441-443`).
    #[test]
    fn insertion_at_cds_first_base_is_utr() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // cDNA 5 = 1004: insertion between 1003 and 1004 is cDNA (5, 4).
        let span = VepSpan::new(1_004, 1_003);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cdna(span), TvRange::new(5, 4));
        assert_eq!(m.cds(span), TvRange::default());
        assert_eq!(m.translation(span), TvRange::default());
    }

    /// `(c, c-1)` with `c - 1 == cdna_coding_end` has `start > cdna_cend`.
    #[test]
    fn insertion_after_cds_last_base_is_utr() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // cDNA 84 is the last CDS base: exon 3 starts at cDNA 61, so 84 is
        // 1074 + 23 = 1097. Insertion between 1097 and 1098: cDNA (85, 84).
        let span = VepSpan::new(1_098, 1_097);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cdna(span), TvRange::new(85, 84));
        assert_eq!(m.cds(span), TvRange::default());
    }

    /// A range starting in the 5'UTR has a leading `Gap` piece
    /// (`TranscriptMapper.pm:450-455`), so `cds_start` is undefined.
    #[test]
    fn deletion_5utr_into_cds_has_undefined_cds_start() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // cDNA 3..7 = 1002..1006, CDS starts at cDNA 5.
        let span = VepSpan::new(1_002, 1_006);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cdna(span), TvRange::new(3, 7));
        assert_eq!(
            m.cds(span),
            TvRange {
                start: None,
                end: Some(3)
            }
        );
        let (tx, exons) = coding_tx(-1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // − strand: cDNA 3..7 = 1101 down to 1097.
        let span = VepSpan::new(1_097, 1_101);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cdna(span), TvRange::new(3, 7));
        assert_eq!(
            m.cds(span),
            TvRange {
                start: None,
                end: Some(3)
            }
        );
    }

    /// A range ending in the 3'UTR has a trailing `Gap` piece
    /// (`TranscriptMapper.pm:458-462`), so `cds_end` is undefined.
    #[test]
    fn deletion_cds_into_3utr_has_undefined_cds_end() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // cDNA 83..86 = 1096..1099, CDS ends at cDNA 84 (CDS 80).
        let span = VepSpan::new(1_096, 1_099);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cdna(span), TvRange::new(83, 86));
        assert_eq!(
            m.cds(span),
            TvRange {
                start: Some(79),
                end: None
            }
        );
        let (tx, exons) = coding_tx(-1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // − strand: exon 1000-1029 holds cDNA 61-90 from 1029 down, so cDNA
        // 83..86 = 1007 down to 1004.
        let span = VepSpan::new(1_004, 1_007);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(m.cdna(span), TvRange::new(83, 86));
        assert_eq!(
            m.cds(span),
            TvRange {
                start: Some(79),
                end: None
            }
        );
    }

    /// BTV `cdna_end` is undefined when the last piece is a `Gap`.
    #[test]
    fn deletion_exon_into_intron_has_undefined_cdna_end() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // 1060..1065: cDNA 59, 60, then the 12-bp intron.
        let span = VepSpan::new(1_060, 1_065);
        assert_eq!(
            mapper(&tx, &refs, 0).cdna(span),
            TvRange {
                start: Some(59),
                end: None
            }
        );
    }

    #[test]
    fn deletion_into_frameshift_intron_has_undefined_end() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        // 1029 is cDNA 30, 1030-1031 the 2-bp intron.
        let span = VepSpan::new(1_029, 1_030);
        let m = mapper(&tx, &refs, 0);
        assert_eq!(
            m.cdna(span),
            TvRange {
                start: Some(30),
                end: None
            }
        );
        assert_eq!(
            m.cdna_pieces(span),
            vec![MapperPiece::Coord { start: 30, end: 30 }, MapperPiece::Gap]
        );
    }

    /// `map_insert` keeps the one `Coordinate` flank and adjusts it
    /// (`Mapper.pm:523-547`).
    #[test]
    fn insertion_with_intronic_flank_keeps_exonic_flank() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let m = mapper(&tx, &refs, 0);
        // Between 1029 (cDNA 30) and the intron: start++ on the 5' flank.
        assert_eq!(m.cdna(VepSpan::new(1_030, 1_029)), TvRange::new(31, 30));
        // Between the intron and 1032 (cDNA 31): end-- on the 3' flank.
        assert_eq!(m.cdna(VepSpan::new(1_032, 1_031)), TvRange::new(31, 30));

        let (tx, exons) = coding_tx(-1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let m = mapper(&tx, &refs, 0);
        // − strand: 1032 is cDNA 60, 1029 cDNA 61.
        assert_eq!(m.cdna(VepSpan::new(1_032, 1_031)), TvRange::new(61, 60));
        assert_eq!(m.cdna(VepSpan::new(1_030, 1_029)), TvRange::new(61, 60));
    }

    /// BTV `cds_start` adds the start exon's phase (`:262-265`) and
    /// `genomic2pep` uses `int((cds + phase + 2) / 3)`.
    #[test]
    fn n_padded_cds_start_nf_phase_offsets_cds_and_peptide() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let m = mapper(&tx, &refs, 2);
        // cDNA 5 = 1004 is raw CDS 1, so CDS 3 and residue 1.
        let span = VepSpan::new(1_004, 1_004);
        assert_eq!(m.cds(span), TvRange::new(3, 3));
        assert_eq!(m.translation(span), TvRange::new(1, 1));
        // cDNA 6 = 1005: CDS 4, residue 2.
        let span = VepSpan::new(1_005, 1_005);
        assert_eq!(m.cds(span), TvRange::new(4, 4));
        assert_eq!(m.translation(span), TvRange::new(2, 2));
    }

    /// `cds_coords` falls back to the unshifted span when the shifted one
    /// starts in a gap and the VF was shifted (`BaseTranscriptVariation.pm:518-526`).
    #[test]
    fn boundary_shift_falls_back_to_unshifted_cds() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let m = mapper(&tx, &refs, 0);
        let shifted = VepSpan::new(1_098, 1_097); // past the CDS end: gap
        let unshifted = shifted.offset(-2); // (1096, 1095): cDNA (83, 82)
        assert_eq!(
            m.cds_shifted(shifted, unshifted, true),
            (TvRange::new(79, 78), true)
        );
        assert_eq!(
            m.cds_shifted(shifted, unshifted, false),
            (TvRange::default(), false)
        );
        assert_eq!(
            m.translation_shifted(shifted, unshifted, true),
            (TvRange::new(27, 26), true)
        );
    }

    /// Perl truthiness: the insertion before cDNA 1 maps to `(1, 0)`, which
    /// `$cdna_start && $cdna_end` treats as undefined.
    #[test]
    fn insertion_before_cdna_1_is_not_truthy() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let range = mapper(&tx, &refs, 0).cdna(VepSpan::new(1_000, 999));
        assert_eq!(range, TvRange::new(1, 0));
        assert_eq!(range.both(), Some((1, 0)));
        assert_eq!(range.both_truthy(), None);
        assert!(!range.overlap(1, 3));
        assert!(TvRange::new(1, 3).overlap(3, 5));
        assert!(!TvRange::new(1, 3).overlap(4, 5));
    }

    #[test]
    fn cdna_pieces_lists_gap_for_intron_and_utr() {
        let (tx, exons) = coding_tx(1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let m = mapper(&tx, &refs, 0);
        // 1027..1034: cDNA 28-30, intron 1030-1031, cDNA 31-33.
        let span = VepSpan::new(1_027, 1_034);
        assert_eq!(
            m.cdna_pieces(span),
            vec![
                MapperPiece::Coord { start: 28, end: 30 },
                MapperPiece::Gap,
                MapperPiece::Coord { start: 31, end: 33 },
            ]
        );
        // 1002..1006 is cDNA 3..7; the CDS starts at cDNA 5.
        assert_eq!(
            m.cds_pieces(VepSpan::new(1_002, 1_006)),
            vec![MapperPiece::Gap, MapperPiece::Coord { start: 1, end: 3 }]
        );
        // − strand: pieces come in transcript order.
        let (tx, exons) = coding_tx(-1);
        let refs: Vec<&ExonFeature> = exons.iter().collect();
        let m = mapper(&tx, &refs, 0);
        // 1027..1034 on −: 1034 is cDNA 58 (exon 2 holds cDNA 31-60 from 1061
        // down), 1032 cDNA 60, intron, 1029 cDNA 61, 1027 cDNA 63.
        assert_eq!(
            m.cdna_pieces(span),
            vec![
                MapperPiece::Coord { start: 58, end: 60 },
                MapperPiece::Gap,
                MapperPiece::Coord { start: 61, end: 63 },
            ]
        );
        // Outside the transcript entirely.
        assert_eq!(
            m.cdna_pieces(VepSpan::new(900, 905)),
            vec![MapperPiece::Gap]
        );
    }
}
