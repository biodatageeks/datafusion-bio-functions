//! Transcript/exon-driven consequence evaluation (phase 2).

use std::collections::{BTreeSet, HashMap};

use crate::so_terms::{ALL_SO_TERMS, SoTerm, unique_sorted_terms};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantInput {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub ref_allele: String,
    pub alt_allele: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptFeature {
    pub transcript_id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: i8,
    pub biotype: String,
    pub cds_start: Option<i64>,
    pub cds_end: Option<i64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExonFeature {
    pub transcript_id: String,
    pub exon_number: i32,
    pub start: i64,
    pub end: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranslationFeature {
    pub transcript_id: String,
    pub cds_len: Option<usize>,
    pub protein_len: Option<usize>,
    pub translation_seq: Option<String>,
    pub cds_sequence: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegulatoryFeature {
    pub feature_id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MotifFeature {
    pub motif_id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MirnaFeature {
    pub mirna_id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvFeatureKind {
    Transcript,
    Regulatory,
    Tfbs,
    Generic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvEventKind {
    Ablation,
    Amplification,
    Elongation,
    Truncation,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuralFeature {
    pub feature_id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub feature_kind: SvFeatureKind,
    pub event_kind: SvEventKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptConsequence {
    pub transcript_id: Option<String>,
    pub terms: Vec<SoTerm>,
}

#[derive(Debug, Clone)]
pub struct TranscriptConsequenceEngine {
    upstream_distance: i64,
    downstream_distance: i64,
}

impl Default for TranscriptConsequenceEngine {
    fn default() -> Self {
        Self::new(5000, 5000)
    }
}

impl TranscriptConsequenceEngine {
    pub fn new(upstream_distance: i64, downstream_distance: i64) -> Self {
        Self {
            upstream_distance,
            downstream_distance,
        }
    }

    pub fn evaluate_variant(
        &self,
        variant: &VariantInput,
        transcripts: &[TranscriptFeature],
        exons: &[ExonFeature],
    ) -> Vec<TranscriptConsequence> {
        self.evaluate_variant_with_context(variant, transcripts, exons, &[], &[], &[], &[], &[])
    }

    pub fn evaluate_variant_with_context(
        &self,
        variant: &VariantInput,
        transcripts: &[TranscriptFeature],
        exons: &[ExonFeature],
        translations: &[TranslationFeature],
        regulatory: &[RegulatoryFeature],
        motifs: &[MotifFeature],
        mirnas: &[MirnaFeature],
        structural: &[StructuralFeature],
    ) -> Vec<TranscriptConsequence> {
        let mut exons_by_tx: HashMap<&str, Vec<&ExonFeature>> = HashMap::new();
        for exon in exons {
            exons_by_tx
                .entry(exon.transcript_id.as_str())
                .or_default()
                .push(exon);
        }
        for tx_exons in exons_by_tx.values_mut() {
            tx_exons.sort_by_key(|e| e.exon_number);
        }
        let translation_by_tx: HashMap<&str, &TranslationFeature> = translations
            .iter()
            .map(|t| (t.transcript_id.as_str(), t))
            .collect();

        let mut out = Vec::new();
        let variant_chrom = normalize_chrom(&variant.chrom);
        for tx in transcripts {
            if normalize_chrom(&tx.chrom) != variant_chrom {
                continue;
            }

            let tx_exons = exons_by_tx
                .get(tx.transcript_id.as_str())
                .cloned()
                .unwrap_or_default();
            let tx_translation = translation_by_tx.get(tx.transcript_id.as_str()).copied();

            if overlaps(variant.start, variant.end, tx.start, tx.end) {
                let terms =
                    self.evaluate_transcript_overlap(variant, tx, &tx_exons, tx_translation);
                if !terms.is_empty() {
                    out.push(TranscriptConsequence {
                        transcript_id: Some(tx.transcript_id.clone()),
                        terms,
                    });
                }
            } else if let Some(term) = self.upstream_downstream_term(variant, tx) {
                out.push(TranscriptConsequence {
                    transcript_id: Some(tx.transcript_id.clone()),
                    terms: vec![term],
                });
            }
        }

        self.append_regulatory_terms(&mut out, variant, regulatory, structural);
        self.append_tfbs_terms(&mut out, variant, motifs, structural);
        self.append_mirna_terms(&mut out, variant, mirnas);
        self.append_structural_transcript_terms(&mut out, variant, transcripts, structural);

        if out.is_empty() {
            out.push(TranscriptConsequence {
                transcript_id: None,
                terms: vec![SoTerm::IntergenicVariant],
            });
        }
        out
    }

    pub fn collapse_variant_terms(assignments: &[TranscriptConsequence]) -> Vec<SoTerm> {
        let mut terms = Vec::new();
        for a in assignments {
            terms.extend_from_slice(&a.terms);
        }
        unique_sorted_terms(terms)
    }

    fn evaluate_transcript_overlap(
        &self,
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
        tx_translation: Option<&TranslationFeature>,
    ) -> Vec<SoTerm> {
        let mut terms = BTreeSet::new();

        let overlaps_exon = tx_exons
            .iter()
            .any(|e| overlaps(variant.start, variant.end, e.start, e.end));

        // Exonic boundary +/- 3bp contributes splice_region.
        if tx_exons.iter().any(|e| {
            overlaps(variant.start, variant.end, e.start, e.start + 2)
                || overlaps(variant.start, variant.end, e.end - 2, e.end)
        }) {
            terms.insert(SoTerm::SpliceRegionVariant);
        }

        if !overlaps_exon {
            terms.insert(SoTerm::IntronVariant);
            self.add_intronic_splice_terms(&mut terms, variant, tx, tx_exons);
        } else if is_non_coding_biotype(&tx.biotype) {
            terms.insert(SoTerm::NonCodingTranscriptExonVariant);
        } else if self.overlaps_cds(variant, tx) {
            self.add_coding_terms(&mut terms, variant, tx, tx_exons, tx_translation);
        } else if let Some(utr_term) = self.utr_term(variant, tx) {
            terms.insert(utr_term);
        } else {
            terms.insert(SoTerm::CodingTranscriptVariant);
        }

        if tx.biotype == "nonsense_mediated_decay" {
            terms.insert(SoTerm::NmdTranscriptVariant);
        }
        if is_non_coding_biotype(&tx.biotype) {
            terms.insert(SoTerm::NonCodingTranscriptVariant);
        } else {
            terms.insert(SoTerm::CodingTranscriptVariant);
        }

        let mut terms_vec: Vec<SoTerm> = terms.into_iter().collect();
        terms_vec.sort_by_key(|t| t.rank());
        terms_vec
    }

    fn append_regulatory_terms(
        &self,
        out: &mut Vec<TranscriptConsequence>,
        variant: &VariantInput,
        regulatory: &[RegulatoryFeature],
        structural: &[StructuralFeature],
    ) {
        let mut terms = BTreeSet::new();
        let chrom = normalize_chrom(&variant.chrom);
        let overlaps_reg = regulatory.iter().any(|r| {
            normalize_chrom(&r.chrom) == chrom
                && overlaps(variant.start, variant.end, r.start, r.end)
        });
        if overlaps_reg {
            terms.insert(SoTerm::RegulatoryRegionVariant);
        }
        for sv in structural {
            if normalize_chrom(&sv.chrom) != chrom
                || !overlaps(variant.start, variant.end, sv.start, sv.end)
            {
                continue;
            }
            if sv.feature_kind != SvFeatureKind::Regulatory {
                continue;
            }
            match sv.event_kind {
                SvEventKind::Ablation => {
                    terms.insert(SoTerm::RegulatoryRegionAblation);
                }
                SvEventKind::Amplification => {
                    terms.insert(SoTerm::RegulatoryRegionAmplification);
                }
                SvEventKind::Elongation | SvEventKind::Truncation => {}
            }
        }
        if !terms.is_empty() {
            let mut ordered: Vec<SoTerm> = terms.into_iter().collect();
            ordered.sort_by_key(|t| t.rank());
            out.push(TranscriptConsequence {
                transcript_id: None,
                terms: ordered,
            });
        }
    }

    fn append_tfbs_terms(
        &self,
        out: &mut Vec<TranscriptConsequence>,
        variant: &VariantInput,
        motifs: &[MotifFeature],
        structural: &[StructuralFeature],
    ) {
        let mut terms = BTreeSet::new();
        let chrom = normalize_chrom(&variant.chrom);
        let overlaps_tfbs = motifs.iter().any(|m| {
            normalize_chrom(&m.chrom) == chrom
                && overlaps(variant.start, variant.end, m.start, m.end)
        });
        if overlaps_tfbs {
            terms.insert(SoTerm::TfBindingSiteVariant);
        }
        for sv in structural {
            if normalize_chrom(&sv.chrom) != chrom
                || !overlaps(variant.start, variant.end, sv.start, sv.end)
            {
                continue;
            }
            if sv.feature_kind != SvFeatureKind::Tfbs {
                continue;
            }
            match sv.event_kind {
                SvEventKind::Ablation => {
                    terms.insert(SoTerm::TfbsAblation);
                }
                SvEventKind::Amplification => {
                    terms.insert(SoTerm::TfbsAmplification);
                }
                SvEventKind::Elongation | SvEventKind::Truncation => {}
            }
        }
        if !terms.is_empty() {
            let mut ordered: Vec<SoTerm> = terms.into_iter().collect();
            ordered.sort_by_key(|t| t.rank());
            out.push(TranscriptConsequence {
                transcript_id: None,
                terms: ordered,
            });
        }
    }

    fn append_mirna_terms(
        &self,
        out: &mut Vec<TranscriptConsequence>,
        variant: &VariantInput,
        mirnas: &[MirnaFeature],
    ) {
        let chrom = normalize_chrom(&variant.chrom);
        let overlaps_mirna = mirnas.iter().any(|m| {
            normalize_chrom(&m.chrom) == chrom
                && overlaps(variant.start, variant.end, m.start, m.end)
        });
        if overlaps_mirna {
            out.push(TranscriptConsequence {
                transcript_id: None,
                terms: vec![SoTerm::MatureMirnaVariant],
            });
        }
    }

    fn append_structural_transcript_terms(
        &self,
        out: &mut Vec<TranscriptConsequence>,
        variant: &VariantInput,
        transcripts: &[TranscriptFeature],
        structural: &[StructuralFeature],
    ) {
        let chrom = normalize_chrom(&variant.chrom);
        let mut terms = BTreeSet::new();
        for sv in structural {
            if normalize_chrom(&sv.chrom) != chrom
                || !overlaps(variant.start, variant.end, sv.start, sv.end)
            {
                continue;
            }
            match sv.feature_kind {
                SvFeatureKind::Transcript => match sv.event_kind {
                    SvEventKind::Ablation => {
                        terms.insert(SoTerm::TranscriptAblation);
                    }
                    SvEventKind::Amplification => {
                        terms.insert(SoTerm::TranscriptAmplification);
                    }
                    SvEventKind::Elongation => {
                        terms.insert(SoTerm::FeatureElongation);
                    }
                    SvEventKind::Truncation => {
                        terms.insert(SoTerm::FeatureTruncation);
                    }
                },
                SvFeatureKind::Generic => match sv.event_kind {
                    SvEventKind::Elongation => {
                        terms.insert(SoTerm::FeatureElongation);
                    }
                    SvEventKind::Truncation => {
                        terms.insert(SoTerm::FeatureTruncation);
                    }
                    SvEventKind::Ablation | SvEventKind::Amplification => {}
                },
                SvFeatureKind::Regulatory | SvFeatureKind::Tfbs => {}
            }
        }
        if !terms.is_empty()
            && transcripts
                .iter()
                .any(|tx| normalize_chrom(&tx.chrom) == chrom)
        {
            let mut ordered: Vec<SoTerm> = terms.into_iter().collect();
            ordered.sort_by_key(|t| t.rank());
            out.push(TranscriptConsequence {
                transcript_id: None,
                terms: ordered,
            });
        }
    }

    fn overlaps_cds(&self, variant: &VariantInput, tx: &TranscriptFeature) -> bool {
        let Some(cds_start) = tx.cds_start else {
            return false;
        };
        let Some(cds_end) = tx.cds_end else {
            return false;
        };
        if cds_start <= 0 || cds_end <= 0 {
            return false;
        }
        overlaps(variant.start, variant.end, cds_start, cds_end)
    }

    fn add_coding_terms(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
        tx_translation: Option<&TranslationFeature>,
    ) {
        let (ref_len, alt_len) = allele_lengths(&variant.ref_allele, &variant.alt_allele);

        if ref_len != alt_len {
            let diff = ref_len.abs_diff(alt_len);
            if diff % 3 == 0 {
                if alt_len > ref_len {
                    terms.insert(SoTerm::InframeInsertion);
                } else {
                    terms.insert(SoTerm::InframeDeletion);
                }
            } else {
                terms.insert(SoTerm::FrameshiftVariant);
            }
            terms.insert(SoTerm::ProteinAlteringVariant);
            terms.insert(SoTerm::CodingSequenceVariant);
            return;
        }

        terms.insert(SoTerm::CodingSequenceVariant);

        if cds_is_incomplete(tx, tx_translation) && self.overlaps_stop_codon(variant, tx) {
            terms.insert(SoTerm::IncompleteTerminalCodonVariant);
        }

        if ref_len == 0 {
            return;
        }

        if ref_len == alt_len {
            if let Some(classification) =
                classify_coding_substitution(tx, tx_exons, tx_translation, variant)
            {
                apply_codon_classification(terms, classification);
                return;
            }
        }

        if self.overlaps_start_codon(variant, tx) {
            if is_start_codon(&variant.ref_allele) && is_start_codon(&variant.alt_allele) {
                terms.insert(SoTerm::StartRetainedVariant);
            } else {
                terms.insert(SoTerm::StartLost);
            }
        }
        if self.overlaps_stop_codon(variant, tx) {
            if is_stop_codon(&variant.ref_allele) && is_stop_codon(&variant.alt_allele) {
                terms.insert(SoTerm::StopRetainedVariant);
            } else if !is_stop_codon(&variant.ref_allele) && is_stop_codon(&variant.alt_allele) {
                terms.insert(SoTerm::StopGained);
            } else {
                terms.insert(SoTerm::StopLost);
            }
        }

        if variant.ref_allele.eq_ignore_ascii_case(&variant.alt_allele) {
            terms.insert(SoTerm::SynonymousVariant);
        } else if ref_len == 1 && alt_len == 1 {
            terms.insert(SoTerm::MissenseVariant);
            terms.insert(SoTerm::ProteinAlteringVariant);
        } else if ref_len % 3 == 0 && alt_len % 3 == 0 {
            if !is_stop_codon(&variant.ref_allele) && is_stop_codon(&variant.alt_allele) {
                terms.insert(SoTerm::StopGained);
            } else {
                terms.insert(SoTerm::ProteinAlteringVariant);
            }
        } else {
            terms.insert(SoTerm::ProteinAlteringVariant);
        }
    }

    fn overlaps_start_codon(&self, variant: &VariantInput, tx: &TranscriptFeature) -> bool {
        let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
            return false;
        };
        let (s, e) = if tx.strand >= 0 {
            (cds_start, cds_start + 2)
        } else {
            (cds_end - 2, cds_end)
        };
        overlaps(variant.start, variant.end, s, e)
    }

    fn overlaps_stop_codon(&self, variant: &VariantInput, tx: &TranscriptFeature) -> bool {
        let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
            return false;
        };
        let (s, e) = if tx.strand >= 0 {
            (cds_end - 2, cds_end)
        } else {
            (cds_start, cds_start + 2)
        };
        overlaps(variant.start, variant.end, s, e)
    }

    fn utr_term(&self, variant: &VariantInput, tx: &TranscriptFeature) -> Option<SoTerm> {
        let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
            return None;
        };
        if tx.strand >= 0 {
            if variant.end < cds_start {
                return Some(SoTerm::FivePrimeUtrVariant);
            }
            if variant.start > cds_end {
                return Some(SoTerm::ThreePrimeUtrVariant);
            }
        } else {
            if variant.end < cds_start {
                return Some(SoTerm::ThreePrimeUtrVariant);
            }
            if variant.start > cds_end {
                return Some(SoTerm::FivePrimeUtrVariant);
            }
        }
        None
    }

    fn upstream_downstream_term(
        &self,
        variant: &VariantInput,
        tx: &TranscriptFeature,
    ) -> Option<SoTerm> {
        if tx.strand >= 0 {
            let up_start = tx.start.saturating_sub(self.upstream_distance);
            let up_end = tx.start.saturating_sub(1);
            if overlaps(variant.start, variant.end, up_start, up_end) {
                return Some(SoTerm::UpstreamGeneVariant);
            }
            let down_start = tx.end.saturating_add(1);
            let down_end = tx.end.saturating_add(self.downstream_distance);
            if overlaps(variant.start, variant.end, down_start, down_end) {
                return Some(SoTerm::DownstreamGeneVariant);
            }
        } else {
            let up_start = tx.end.saturating_add(1);
            let up_end = tx.end.saturating_add(self.upstream_distance);
            if overlaps(variant.start, variant.end, up_start, up_end) {
                return Some(SoTerm::UpstreamGeneVariant);
            }
            let down_start = tx.start.saturating_sub(self.downstream_distance);
            let down_end = tx.start.saturating_sub(1);
            if overlaps(variant.start, variant.end, down_start, down_end) {
                return Some(SoTerm::DownstreamGeneVariant);
            }
        }
        None
    }

    fn add_intronic_splice_terms(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
    ) {
        for exon in tx_exons {
            if tx.strand >= 0 {
                // Donor (intron after exon end): +1/+2, +5, +3..+6, +7/+8
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 1,
                    exon.end + 2,
                    SoTerm::SpliceDonorVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 5,
                    exon.end + 5,
                    SoTerm::SpliceDonor5thBaseVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 3,
                    exon.end + 6,
                    SoTerm::SpliceDonorRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 7,
                    exon.end + 8,
                    SoTerm::SpliceRegionVariant,
                );

                // Acceptor (before exon start): -1/-2, -3..-8, -3..-17
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 2,
                    exon.start - 1,
                    SoTerm::SpliceAcceptorVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 8,
                    exon.start - 3,
                    SoTerm::SpliceRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 17,
                    exon.start - 3,
                    SoTerm::SplicePolypyrimidineTractVariant,
                );
            } else {
                // Donor for negative strand lives before exon start.
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 2,
                    exon.start - 1,
                    SoTerm::SpliceDonorVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 5,
                    exon.start - 5,
                    SoTerm::SpliceDonor5thBaseVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 6,
                    exon.start - 3,
                    SoTerm::SpliceDonorRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.start - 8,
                    exon.start - 7,
                    SoTerm::SpliceRegionVariant,
                );

                // Acceptor for negative strand lives after exon end.
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 1,
                    exon.end + 2,
                    SoTerm::SpliceAcceptorVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 3,
                    exon.end + 8,
                    SoTerm::SpliceRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    variant,
                    exon.end + 3,
                    exon.end + 17,
                    SoTerm::SplicePolypyrimidineTractVariant,
                );
            }
        }
    }
}

pub fn implemented_phase2_terms() -> Vec<SoTerm> {
    ALL_SO_TERMS.to_vec()
}

pub fn missing_phase2_terms() -> Vec<SoTerm> {
    let implemented: BTreeSet<SoTerm> = implemented_phase2_terms().into_iter().collect();
    ALL_SO_TERMS
        .iter()
        .copied()
        .filter(|t| !implemented.contains(t))
        .collect()
}

fn normalize_chrom(chrom: &str) -> &str {
    chrom.strip_prefix("chr").unwrap_or(chrom)
}

fn overlaps(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> bool {
    let (a_start, a_end) = if a_start <= a_end {
        (a_start, a_end)
    } else {
        (a_end, a_start)
    };
    let (b_start, b_end) = if b_start <= b_end {
        (b_start, b_end)
    } else {
        (b_end, b_start)
    };
    a_start <= b_end && a_end >= b_start
}

fn add_if_overlaps(
    terms: &mut BTreeSet<SoTerm>,
    variant: &VariantInput,
    start: i64,
    end: i64,
    term: SoTerm,
) {
    if overlaps(variant.start, variant.end, start, end) {
        terms.insert(term);
    }
}

fn is_non_coding_biotype(biotype: &str) -> bool {
    biotype != "protein_coding"
}

fn allele_lengths(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_len = if ref_allele == "-" {
        0
    } else {
        ref_allele.len()
    };
    let alt_len = if alt_allele == "-" {
        0
    } else {
        alt_allele.len()
    };
    (ref_len, alt_len)
}

fn cds_is_incomplete(tx: &TranscriptFeature, tx_translation: Option<&TranslationFeature>) -> bool {
    if let Some(translation) = tx_translation {
        if let Some(cds_len) = translation.cds_len {
            return cds_len % 3 != 0;
        }
    }
    let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
        return false;
    };
    if cds_start <= 0 || cds_end <= 0 || cds_end < cds_start {
        return false;
    }
    let cds_len = cds_end - cds_start + 1;
    cds_len % 3 != 0
}

fn is_start_codon(allele: &str) -> bool {
    allele.eq_ignore_ascii_case("ATG")
}

fn is_stop_codon(allele: &str) -> bool {
    matches!(
        allele.to_ascii_uppercase().as_str(),
        "TAA" | "TAG" | "TGA" | "*"
    )
}

#[derive(Debug, Clone, Copy, Default)]
struct CodingClassification {
    synonymous: bool,
    missense: bool,
    stop_gained: bool,
    stop_lost: bool,
    stop_retained: bool,
    start_lost: bool,
    start_retained: bool,
}

fn apply_codon_classification(terms: &mut BTreeSet<SoTerm>, c: CodingClassification) {
    if c.start_lost {
        terms.insert(SoTerm::StartLost);
    }
    if c.start_retained {
        terms.insert(SoTerm::StartRetainedVariant);
    }
    if c.stop_gained {
        terms.insert(SoTerm::StopGained);
    }
    if c.stop_lost {
        terms.insert(SoTerm::StopLost);
    }
    if c.stop_retained {
        terms.insert(SoTerm::StopRetainedVariant);
    }
    if c.synonymous {
        terms.insert(SoTerm::SynonymousVariant);
    }
    if c.missense {
        terms.insert(SoTerm::MissenseVariant);
        terms.insert(SoTerm::ProteinAlteringVariant);
    }
}

fn classify_coding_substitution(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    tx_translation: Option<&TranslationFeature>,
    variant: &VariantInput,
) -> Option<CodingClassification> {
    let translation = tx_translation?;
    let cds_seq = translation.cds_sequence.as_deref()?.to_ascii_uppercase();
    let ref_len = variant.ref_allele.len();
    let alt_len = variant.alt_allele.len();
    if ref_len == 0 || ref_len != alt_len {
        return None;
    }

    let genomic_positions = genomic_range(variant.start, variant.end)?;
    if genomic_positions.len() != ref_len {
        return None;
    }

    let mut cds_indices = Vec::with_capacity(genomic_positions.len());
    for pos in &genomic_positions {
        cds_indices.push(genomic_to_cds_index(tx, tx_exons, *pos)?);
    }
    cds_indices.sort_unstable();
    if cds_indices
        .windows(2)
        .any(|w| w[1] != w[0].saturating_add(1))
    {
        return None;
    }
    let start_idx = *cds_indices.first()?;
    let end_idx = *cds_indices.last()?;
    if end_idx >= cds_seq.len() {
        return None;
    }

    let ref_tx = if tx.strand >= 0 {
        variant.ref_allele.to_ascii_uppercase()
    } else {
        reverse_complement(&variant.ref_allele)?.to_ascii_uppercase()
    };
    let alt_tx = if tx.strand >= 0 {
        variant.alt_allele.to_ascii_uppercase()
    } else {
        reverse_complement(&variant.alt_allele)?.to_ascii_uppercase()
    };

    let ref_seq_slice = &cds_seq[start_idx..=end_idx];
    if ref_seq_slice != ref_tx {
        return None;
    }

    let mut mutated = cds_seq.as_bytes().to_vec();
    let alt_bytes = alt_tx.as_bytes();
    if alt_bytes.len() != ref_len {
        return None;
    }
    for (i, b) in alt_bytes.iter().enumerate() {
        mutated[start_idx + i] = *b;
    }

    let codon_start = start_idx / 3;
    let codon_end = end_idx / 3;
    let mut class = CodingClassification::default();
    let mut aa_changed = false;
    for codon_idx in codon_start..=codon_end {
        let codon_base = codon_idx * 3;
        if codon_base + 2 >= mutated.len() {
            continue;
        }
        let old_codon = &cds_seq[codon_base..codon_base + 3];
        let new_codon = std::str::from_utf8(&mutated[codon_base..codon_base + 3]).ok()?;
        let old_aa = translate_codon(old_codon)?;
        let new_aa = translate_codon(new_codon)?;

        if codon_idx == 0 && old_aa == 'M' {
            if new_aa == 'M' {
                class.start_retained = true;
            } else {
                class.start_lost = true;
            }
        }
        if old_aa == '*' && new_aa == '*' {
            class.stop_retained = true;
        } else if old_aa == '*' && new_aa != '*' {
            class.stop_lost = true;
        } else if old_aa != '*' && new_aa == '*' {
            class.stop_gained = true;
        }
        if old_aa != new_aa {
            aa_changed = true;
        }
    }

    if aa_changed
        && !class.stop_gained
        && !class.stop_lost
        && !class.start_lost
        && !class.stop_retained
    {
        class.missense = true;
    } else if !aa_changed && !class.stop_retained && !class.start_retained {
        class.synonymous = true;
    }

    Some(class)
}

fn genomic_range(start: i64, end: i64) -> Option<Vec<i64>> {
    let (s, e) = if start <= end {
        (start, end)
    } else {
        (end, start)
    };
    let len = e.saturating_sub(s).saturating_add(1);
    let len_usize = usize::try_from(len).ok()?;
    let mut out = Vec::with_capacity(len_usize);
    for p in s..=e {
        out.push(p);
    }
    Some(out)
}

fn genomic_to_cds_index(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    pos: i64,
) -> Option<usize> {
    let segments = coding_segments(tx, tx_exons)?;
    let mut offset = 0usize;
    for (seg_start, seg_end) in segments {
        let seg_len = usize::try_from(seg_end.saturating_sub(seg_start).saturating_add(1)).ok()?;
        if pos >= seg_start && pos <= seg_end {
            let local = if tx.strand >= 0 {
                usize::try_from(pos.saturating_sub(seg_start)).ok()?
            } else {
                usize::try_from(seg_end.saturating_sub(pos)).ok()?
            };
            return Some(offset + local);
        }
        offset = offset.saturating_add(seg_len);
    }
    None
}

fn coding_segments(tx: &TranscriptFeature, tx_exons: &[&ExonFeature]) -> Option<Vec<(i64, i64)>> {
    let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
        return None;
    };
    if cds_start <= 0 || cds_end <= 0 || cds_end < cds_start {
        return None;
    }
    let mut segments = Vec::new();
    for exon in tx_exons {
        let s = exon.start.max(cds_start);
        let e = exon.end.min(cds_end);
        if s <= e {
            segments.push((s, e));
        }
    }
    if segments.is_empty() {
        return None;
    }
    segments.sort_by_key(|(s, _)| *s);
    if tx.strand < 0 {
        segments.reverse();
    }
    Some(segments)
}

fn reverse_complement(seq: &str) -> Option<String> {
    let mut out = String::with_capacity(seq.len());
    for b in seq.as_bytes().iter().rev() {
        let c = match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            _ => return None,
        };
        out.push(c);
    }
    Some(out)
}

fn translate_codon(codon: &str) -> Option<char> {
    match codon {
        "TTT" | "TTC" => Some('F'),
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => Some('L'),
        "ATT" | "ATC" | "ATA" => Some('I'),
        "ATG" => Some('M'),
        "GTT" | "GTC" | "GTA" | "GTG" => Some('V'),
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => Some('S'),
        "CCT" | "CCC" | "CCA" | "CCG" => Some('P'),
        "ACT" | "ACC" | "ACA" | "ACG" => Some('T'),
        "GCT" | "GCC" | "GCA" | "GCG" => Some('A'),
        "TAT" | "TAC" => Some('Y'),
        "CAT" | "CAC" => Some('H'),
        "CAA" | "CAG" => Some('Q'),
        "AAT" | "AAC" => Some('N'),
        "AAA" | "AAG" => Some('K'),
        "GAT" | "GAC" => Some('D'),
        "GAA" | "GAG" => Some('E'),
        "TGT" | "TGC" => Some('C'),
        "TGG" => Some('W'),
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Some('R'),
        "GGT" | "GGC" | "GGA" | "GGG" => Some('G'),
        "TAA" | "TAG" | "TGA" => Some('*'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(
        id: &str,
        chrom: &str,
        start: i64,
        end: i64,
        strand: i8,
        biotype: &str,
        cds_start: Option<i64>,
        cds_end: Option<i64>,
    ) -> TranscriptFeature {
        TranscriptFeature {
            transcript_id: id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            biotype: biotype.to_string(),
            cds_start,
            cds_end,
        }
    }

    fn exon(tx_id: &str, num: i32, start: i64, end: i64) -> ExonFeature {
        ExonFeature {
            transcript_id: tx_id.to_string(),
            exon_number: num,
            start,
            end,
        }
    }

    fn var(chrom: &str, start: i64, end: i64, r: &str, a: &str) -> VariantInput {
        VariantInput {
            chrom: chrom.to_string(),
            start,
            end,
            ref_allele: r.to_string(),
            alt_allele: a.to_string(),
        }
    }

    fn regulatory(id: &str, chrom: &str, start: i64, end: i64) -> RegulatoryFeature {
        RegulatoryFeature {
            feature_id: id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    fn motif(id: &str, chrom: &str, start: i64, end: i64) -> MotifFeature {
        MotifFeature {
            motif_id: id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    fn mirna(id: &str, chrom: &str, start: i64, end: i64) -> MirnaFeature {
        MirnaFeature {
            mirna_id: id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    fn sv(
        id: &str,
        chrom: &str,
        start: i64,
        end: i64,
        feature_kind: SvFeatureKind,
        event_kind: SvEventKind,
    ) -> StructuralFeature {
        StructuralFeature {
            feature_id: id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            feature_kind,
            event_kind,
        }
    }

    #[test]
    fn intergenic_when_no_transcripts_overlap_or_nearby() {
        let engine = TranscriptConsequenceEngine::default();
        let assignments = engine.evaluate_variant(
            &var("22", 1_000_000, 1_000_000, "A", "G"),
            &[tx(
                "tx1",
                "22",
                100,
                200,
                1,
                "protein_coding",
                Some(120),
                Some(180),
            )],
            &[exon("tx1", 1, 100, 200)],
        );
        assert_eq!(assignments.len(), 1);
        assert_eq!(assignments[0].terms, vec![SoTerm::IntergenicVariant]);
    }

    #[test]
    fn upstream_downstream_strand_aware() {
        let engine = TranscriptConsequenceEngine::new(5000, 5000);
        let positive = tx(
            "txp",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1100),
            Some(1900),
        );
        let negative = tx(
            "txn",
            "22",
            3000,
            4000,
            -1,
            "protein_coding",
            Some(3100),
            Some(3900),
        );

        let up_p = engine.evaluate_variant(
            &var("22", 900, 900, "A", "G"),
            std::slice::from_ref(&positive),
            &[],
        );
        assert_eq!(up_p[0].terms, vec![SoTerm::UpstreamGeneVariant]);

        let down_p = engine.evaluate_variant(
            &var("22", 2100, 2100, "A", "G"),
            std::slice::from_ref(&positive),
            &[],
        );
        assert_eq!(down_p[0].terms, vec![SoTerm::DownstreamGeneVariant]);

        let up_n = engine.evaluate_variant(
            &var("22", 4100, 4100, "A", "G"),
            std::slice::from_ref(&negative),
            &[],
        );
        assert_eq!(up_n[0].terms, vec![SoTerm::UpstreamGeneVariant]);

        let down_n = engine.evaluate_variant(
            &var("22", 2900, 2900, "A", "G"),
            std::slice::from_ref(&negative),
            &[],
        );
        assert_eq!(down_n[0].terms, vec![SoTerm::DownstreamGeneVariant]);
    }

    #[test]
    fn non_coding_exon_and_intron_terms() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx("lnc", "22", 100, 300, 1, "lincRNA", None, None);
        let exons = vec![exon("lnc", 1, 100, 150), exon("lnc", 2, 250, 300)];

        let exonic = engine.evaluate_variant(
            &var("22", 120, 120, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        let terms = &exonic[0].terms;
        assert!(terms.contains(&SoTerm::NonCodingTranscriptExonVariant));
        assert!(terms.contains(&SoTerm::NonCodingTranscriptVariant));

        let intronic = engine.evaluate_variant(
            &var("22", 200, 200, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        let terms = &intronic[0].terms;
        assert!(terms.contains(&SoTerm::IntronVariant));
        assert!(terms.contains(&SoTerm::NonCodingTranscriptVariant));
    }

    #[test]
    fn coding_indels_emit_frameshift_or_inframe() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            100,
            300,
            1,
            "protein_coding",
            Some(120),
            Some(280),
        );
        let exons = vec![exon("pc", 1, 100, 300)];

        let frameshift = engine.evaluate_variant(
            &var("22", 150, 150, "A", "AT"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(frameshift[0].terms.contains(&SoTerm::FrameshiftVariant));

        let inframe_ins = engine.evaluate_variant(
            &var("22", 150, 150, "A", "ATGC"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(inframe_ins[0].terms.contains(&SoTerm::InframeInsertion));

        let inframe_del = engine.evaluate_variant(
            &var("22", 150, 152, "ATGC", "A"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(inframe_del[0].terms.contains(&SoTerm::InframeDeletion));
    }

    #[test]
    fn coding_substitution_emits_missense_variant() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            100,
            300,
            1,
            "protein_coding",
            Some(120),
            Some(280),
        );
        let exons = vec![exon("pc", 1, 100, 300)];
        let assignments = engine.evaluate_variant(
            &var("22", 150, 150, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(assignments[0].terms.contains(&SoTerm::MissenseVariant));
    }

    #[test]
    fn utr_terms_are_strand_aware() {
        let engine = TranscriptConsequenceEngine::default();
        let tx_pos = tx(
            "pcp",
            "22",
            100,
            300,
            1,
            "protein_coding",
            Some(150),
            Some(250),
        );
        let tx_neg = tx(
            "pcn",
            "22",
            100,
            300,
            -1,
            "protein_coding",
            Some(150),
            Some(250),
        );
        let exons = vec![exon("pcp", 1, 100, 300), exon("pcn", 1, 100, 300)];

        let pos_5p = engine.evaluate_variant(
            &var("22", 120, 120, "A", "G"),
            std::slice::from_ref(&tx_pos),
            &exons,
        );
        assert!(pos_5p[0].terms.contains(&SoTerm::FivePrimeUtrVariant));
        let pos_3p = engine.evaluate_variant(
            &var("22", 280, 280, "A", "G"),
            std::slice::from_ref(&tx_pos),
            &exons,
        );
        assert!(pos_3p[0].terms.contains(&SoTerm::ThreePrimeUtrVariant));

        let neg_5p = engine.evaluate_variant(
            &var("22", 280, 280, "A", "G"),
            std::slice::from_ref(&tx_neg),
            &exons,
        );
        assert!(neg_5p[0].terms.contains(&SoTerm::FivePrimeUtrVariant));
        let neg_3p = engine.evaluate_variant(
            &var("22", 120, 120, "A", "G"),
            std::slice::from_ref(&tx_neg),
            &exons,
        );
        assert!(neg_3p[0].terms.contains(&SoTerm::ThreePrimeUtrVariant));
    }

    #[test]
    fn splice_terms_from_intronic_offsets() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            100,
            300,
            1,
            "protein_coding",
            Some(120),
            Some(280),
        );
        let exons = vec![exon("pc", 1, 100, 150), exon("pc", 2, 250, 300)];

        let donor = engine.evaluate_variant(
            &var("22", 151, 151, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(donor[0].terms.contains(&SoTerm::SpliceDonorVariant));
        let acceptor = engine.evaluate_variant(
            &var("22", 248, 248, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(acceptor[0].terms.contains(&SoTerm::SpliceAcceptorVariant));
    }

    #[test]
    fn start_and_stop_regions_flagged_in_cds() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            100,
            300,
            1,
            "protein_coding",
            Some(150),
            Some(240),
        );
        let exons = vec![exon("pc", 1, 100, 300)];

        let start_hit = engine.evaluate_variant(
            &var("22", 151, 151, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(start_hit[0].terms.contains(&SoTerm::StartLost));
        let stop_hit = engine.evaluate_variant(
            &var("22", 239, 239, "A", "G"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(stop_hit[0].terms.contains(&SoTerm::StopLost));
    }

    #[test]
    fn retained_and_gained_stop_terms_are_emitted() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            100,
            400,
            1,
            "protein_coding",
            Some(150),
            Some(360),
        );
        let exons = vec![exon("pc", 1, 100, 400)];

        let stop_retained = engine.evaluate_variant(
            &var("22", 359, 361, "TAA", "TGA"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(
            stop_retained[0]
                .terms
                .contains(&SoTerm::StopRetainedVariant)
        );

        let stop_gained = engine.evaluate_variant(
            &var("22", 220, 222, "CAA", "TAA"),
            std::slice::from_ref(&tx),
            &exons,
        );
        assert!(stop_gained[0].terms.contains(&SoTerm::StopGained));
    }

    #[test]
    fn start_retained_and_incomplete_terminal_codon_terms() {
        let engine = TranscriptConsequenceEngine::default();
        let tx_complete = tx(
            "pc",
            "22",
            100,
            350,
            1,
            "protein_coding",
            Some(151),
            Some(240), // 90bp (complete)
        );
        let tx_incomplete = tx(
            "pc2",
            "22",
            100,
            350,
            1,
            "protein_coding",
            Some(151),
            Some(241), // 91bp (incomplete terminal codon)
        );
        let exons = vec![exon("pc", 1, 100, 350), exon("pc2", 1, 100, 350)];

        let start_retained = engine.evaluate_variant(
            &var("22", 151, 153, "ATG", "ATG"),
            std::slice::from_ref(&tx_complete),
            &exons,
        );
        assert!(
            start_retained[0]
                .terms
                .contains(&SoTerm::StartRetainedVariant)
        );

        let incomplete = engine.evaluate_variant(
            &var("22", 239, 241, "TAA", "TAG"),
            std::slice::from_ref(&tx_incomplete),
            &exons,
        );
        assert!(
            incomplete[0]
                .terms
                .contains(&SoTerm::IncompleteTerminalCodonVariant)
        );
    }

    #[test]
    fn context_features_emit_regulatory_tfbs_mirna_and_sv_terms() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "tx1",
            "22",
            100,
            250,
            1,
            "protein_coding",
            Some(120),
            Some(240),
        );
        let exons = vec![exon("tx1", 1, 100, 250)];
        let regulatory_features = vec![regulatory("reg1", "22", 140, 180)];
        let motifs = vec![motif("motif1", "22", 150, 160)];
        let mirnas = vec![mirna("mir1", "22", 155, 165)];
        let structural = vec![
            sv(
                "sv_tx_del",
                "22",
                150,
                160,
                SvFeatureKind::Transcript,
                SvEventKind::Ablation,
            ),
            sv(
                "sv_tx_amp",
                "22",
                150,
                160,
                SvFeatureKind::Transcript,
                SvEventKind::Amplification,
            ),
            sv(
                "sv_tx_elong",
                "22",
                150,
                160,
                SvFeatureKind::Transcript,
                SvEventKind::Elongation,
            ),
            sv(
                "sv_tx_trunc",
                "22",
                150,
                160,
                SvFeatureKind::Transcript,
                SvEventKind::Truncation,
            ),
            sv(
                "sv_reg_del",
                "22",
                150,
                160,
                SvFeatureKind::Regulatory,
                SvEventKind::Ablation,
            ),
            sv(
                "sv_reg_amp",
                "22",
                150,
                160,
                SvFeatureKind::Regulatory,
                SvEventKind::Amplification,
            ),
            sv(
                "sv_tfbs_del",
                "22",
                150,
                160,
                SvFeatureKind::Tfbs,
                SvEventKind::Ablation,
            ),
            sv(
                "sv_tfbs_amp",
                "22",
                150,
                160,
                SvFeatureKind::Tfbs,
                SvEventKind::Amplification,
            ),
        ];

        let assignments = engine.evaluate_variant_with_context(
            &var("22", 155, 155, "A", "G"),
            &[tx],
            &exons,
            &[],
            &regulatory_features,
            &motifs,
            &mirnas,
            &structural,
        );
        let collapsed = TranscriptConsequenceEngine::collapse_variant_terms(&assignments);

        assert!(collapsed.contains(&SoTerm::RegulatoryRegionVariant));
        assert!(collapsed.contains(&SoTerm::RegulatoryRegionAblation));
        assert!(collapsed.contains(&SoTerm::RegulatoryRegionAmplification));
        assert!(collapsed.contains(&SoTerm::TfBindingSiteVariant));
        assert!(collapsed.contains(&SoTerm::TfbsAblation));
        assert!(collapsed.contains(&SoTerm::TfbsAmplification));
        assert!(collapsed.contains(&SoTerm::MatureMirnaVariant));
        assert!(collapsed.contains(&SoTerm::TranscriptAblation));
        assert!(collapsed.contains(&SoTerm::TranscriptAmplification));
        assert!(collapsed.contains(&SoTerm::FeatureElongation));
        assert!(collapsed.contains(&SoTerm::FeatureTruncation));
    }

    #[test]
    fn phase2_coverage_lists_sum_to_41() {
        let implemented = implemented_phase2_terms();
        let missing = missing_phase2_terms();
        assert_eq!(implemented.len() + missing.len(), 41);
        assert!(implemented.contains(&SoTerm::FrameshiftVariant));
        assert!(implemented.contains(&SoTerm::MissenseVariant));
        assert!(implemented.contains(&SoTerm::RegulatoryRegionVariant));
        assert!(missing.is_empty());
    }
}
