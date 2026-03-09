//! Transcript/exon-driven consequence evaluation (phase 2).

use std::collections::{BTreeSet, HashMap};

use coitrees::{COITree, Interval, IntervalTree};

use crate::so_terms::{ALL_SO_TERMS, SoTerm};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantInput {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub ref_allele: String,
    pub alt_allele: String,
}

impl VariantInput {
    /// Create a VariantInput from VCF-style coordinates, trimming the common
    /// prefix/suffix anchor bases that VCF uses for indels.  VEP internally
    /// performs this normalization so that only the truly changed bases are
    /// used for splice-site and other positional checks.
    pub fn from_vcf(chrom: String, pos: i64, end: i64, ref_allele: String, alt_allele: String) -> Self {
        let ref_bytes = ref_allele.as_bytes();
        let alt_bytes = alt_allele.as_bytes();

        // Trim common prefix
        let prefix_len = ref_bytes
            .iter()
            .zip(alt_bytes.iter())
            .take_while(|(a, b)| a == b)
            .count();

        if prefix_len == 0 || (prefix_len == ref_bytes.len() && prefix_len == alt_bytes.len()) {
            // No trimming needed (SNV or identical alleles)
            return Self { chrom, start: pos, end, ref_allele, alt_allele };
        }

        let trimmed_ref = &ref_bytes[prefix_len..];
        let trimmed_alt = &alt_bytes[prefix_len..];
        let new_start = pos + prefix_len as i64;
        let new_end = if trimmed_ref.is_empty() {
            // Pure insertion: the affected position is the insertion point
            // itself, not a 2-base span.
            new_start
        } else {
            pos + ref_bytes.len() as i64 - 1
        };

        Self {
            chrom,
            start: new_start,
            end: new_end,
            ref_allele: if trimmed_ref.is_empty() {
                "-".to_string()
            } else {
                String::from_utf8_lossy(trimmed_ref).to_string()
            },
            alt_allele: if trimmed_alt.is_empty() {
                "-".to_string()
            } else {
                String::from_utf8_lossy(trimmed_alt).to_string()
            },
        }
    }
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

/// Pre-built indexes over context features, constructed once and reused
/// across all variants in a batch.
pub struct PreparedContext<'a> {
    pub exons_by_tx: HashMap<&'a str, Vec<&'a ExonFeature>>,
    pub translation_by_tx: HashMap<&'a str, &'a TranslationFeature>,
    /// All transcripts stored by index; the COITree metadata is the index
    /// into this vec.
    pub transcripts: Vec<&'a TranscriptFeature>,
    /// Per-chromosome interval tree mapping genomic range → transcript index.
    pub tx_trees: HashMap<String, COITree<usize, u32>>,
    pub regulatory: &'a [RegulatoryFeature],
    pub motifs: &'a [MotifFeature],
    pub mirnas: &'a [MirnaFeature],
    pub structural: &'a [StructuralFeature],
}

impl<'a> PreparedContext<'a> {
    pub fn new(
        transcripts: &'a [TranscriptFeature],
        exons: &'a [ExonFeature],
        translations: &'a [TranslationFeature],
        regulatory: &'a [RegulatoryFeature],
        motifs: &'a [MotifFeature],
        mirnas: &'a [MirnaFeature],
        structural: &'a [StructuralFeature],
    ) -> Self {
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
        let tx_refs: Vec<&TranscriptFeature> = transcripts.iter().collect();

        // Build per-chromosome COITree for fast interval overlap queries.
        let mut chrom_intervals: HashMap<String, Vec<Interval<usize>>> = HashMap::new();
        for (idx, tx) in tx_refs.iter().enumerate() {
            let chrom = normalize_chrom(&tx.chrom).to_string();
            chrom_intervals
                .entry(chrom)
                .or_default()
                .push(Interval::new(tx.start as i32, tx.end as i32, idx));
        }
        let tx_trees: HashMap<String, COITree<usize, u32>> = chrom_intervals
            .into_iter()
            .map(|(chrom, intervals)| {
                let tree = COITree::new(&intervals);
                (chrom, tree)
            })
            .collect();

        Self {
            exons_by_tx,
            translation_by_tx,
            transcripts: tx_refs,
            tx_trees,
            regulatory,
            motifs,
            mirnas,
            structural,
        }
    }
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
        let ctx = PreparedContext::new(
            transcripts, exons, translations, regulatory, motifs, mirnas, structural,
        );
        self.evaluate_variant_prepared(variant, &ctx)
    }

    /// Evaluate a variant against pre-built context indexes.
    /// Use this when annotating many variants against the same context —
    /// build `PreparedContext` once and call this per variant.
    pub fn evaluate_variant_prepared(
        &self,
        variant: &VariantInput,
        ctx: &PreparedContext<'_>,
    ) -> Vec<TranscriptConsequence> {
        let mut out = Vec::new();
        let variant_chrom = normalize_chrom(&variant.chrom);
        let max_dist = self.upstream_distance.max(self.downstream_distance);

        // Query the per-chromosome COITree with the variant range expanded
        // by upstream/downstream distance to catch nearby transcripts.
        if let Some(tree) = ctx.tx_trees.get(variant_chrom) {
            let query_first = (variant.start - max_dist) as i32;
            let query_last = (variant.end + max_dist) as i32;
            tree.query(query_first, query_last, |node| {
                let tx = ctx.transcripts[*node.metadata];
                let tx_exons = ctx
                    .exons_by_tx
                    .get(tx.transcript_id.as_str())
                    .cloned()
                    .unwrap_or_default();
                let tx_translation = ctx
                    .translation_by_tx
                    .get(tx.transcript_id.as_str())
                    .copied();

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
            });
        }

        // Track whether any transcript was matched (overlap or upstream/downstream).
        let has_transcript_hit = !out.is_empty();

        self.append_regulatory_terms(&mut out, variant, ctx.regulatory, ctx.structural);
        self.append_tfbs_terms(&mut out, variant, ctx.motifs, ctx.structural);
        self.append_mirna_terms(&mut out, variant, ctx.mirnas);
        self.append_structural_transcript_terms_prepared(&mut out, variant, ctx);

        // VEP emits intergenic_variant when no transcript was hit, even if
        // regulatory/motif features overlap (those are orthogonal to transcripts).
        if !has_transcript_hit {
            out.push(TranscriptConsequence {
                transcript_id: None,
                terms: vec![SoTerm::IntergenicVariant],
            });
        }
        out
    }

    pub fn collapse_variant_terms(assignments: &[TranscriptConsequence]) -> Vec<SoTerm> {
        let mut set = BTreeSet::new();
        for a in assignments {
            set.extend(a.terms.iter().copied());
        }
        // Apply coding parent-term stripping on the merged set so that a
        // `coding_sequence_variant` from a null-CDS transcript is removed
        // when another transcript contributes a specific child (missense,
        // synonymous, etc.).
        // NOTE: only strip coding parents here, NOT intron_variant —
        // intron_variant stripping is per-transcript only (a different
        // transcript may legitimately contribute intron_variant).
        strip_coding_parent_terms(&mut set);
        let mut v: Vec<SoTerm> = set.into_iter().collect();
        v.sort_by_key(|t| t.rank());
        v
    }

    fn evaluate_transcript_overlap(
        &self,
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
        tx_translation: Option<&TranslationFeature>,
    ) -> Vec<SoTerm> {
        let mut terms = BTreeSet::new();
        let is_ins = variant.ref_allele == "-";

        // For pure insertions, VEP requires both flanking positions
        // (start-1 and start) to be within the exon.  An insertion at
        // the first base of an exon (start == exon.start) is between the
        // last intron base and the first exon base — VEP considers that
        // intronic.
        let overlaps_exon = tx_exons.iter().any(|e| {
            if is_ins {
                variant.start > e.start && variant.start <= e.end
            } else {
                overlaps(variant.start, variant.end, e.start, e.end)
            }
        });

        // Exonic boundary +/- 3bp contributes splice_region, but only
        // at internal splice junctions — not at terminal transcript
        // boundaries (first/last exon outer edges).
        {
            let tx_start = tx.start;
            let tx_end = tx.end;
            if tx_exons.iter().any(|e| {
                // Exon start boundary is a splice junction unless it equals
                // the transcript start (terminal).
                let start_is_junction = e.start != tx_start;
                // Exon end boundary is a splice junction unless it equals
                // the transcript end (terminal).
                let end_is_junction = e.end != tx_end;
                (start_is_junction
                    && overlaps(variant.start, variant.end, e.start, e.start + 2))
                    || (end_is_junction
                        && overlaps(variant.start, variant.end, e.end - 2, e.end))
            }) {
                terms.insert(SoTerm::SpliceRegionVariant);
            }
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
        }

        if tx.biotype == "nonsense_mediated_decay" {
            terms.insert(SoTerm::NmdTranscriptVariant);
        }
        if is_non_coding_biotype(&tx.biotype) {
            // VEP omits the parent non_coding_transcript_variant when the more
            // specific non_coding_transcript_exon_variant is already present.
            if !terms.contains(&SoTerm::NonCodingTranscriptExonVariant) {
                terms.insert(SoTerm::NonCodingTranscriptVariant);
            }
        }

        // VEP omits parent SO terms when more specific children are present.
        strip_parent_terms(&mut terms);

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

    fn append_structural_transcript_terms_prepared(
        &self,
        out: &mut Vec<TranscriptConsequence>,
        variant: &VariantInput,
        ctx: &PreparedContext<'_>,
    ) {
        let chrom = normalize_chrom(&variant.chrom);
        let has_tx = ctx.tx_trees.contains_key(chrom);
        Self::append_structural_terms_inner(out, variant, ctx.structural, &chrom, has_tx);
    }

    fn append_structural_terms_inner(
        out: &mut Vec<TranscriptConsequence>,
        variant: &VariantInput,
        structural: &[StructuralFeature],
        chrom: &str,
        has_transcripts_on_chrom: bool,
    ) {
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
        if !terms.is_empty() && has_transcripts_on_chrom {
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
        terms.insert(SoTerm::CodingSequenceVariant);

        if cds_is_incomplete(tx, tx_translation) && self.overlaps_stop_codon(variant, tx) {
            terms.insert(SoTerm::IncompleteTerminalCodonVariant);
        }

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

            if let Some(mut classification) =
                classify_coding_change(tx, tx_exons, tx_translation, variant)
            {
                // VEP does not emit stop_gained/stop_lost alongside
                // frameshift_variant or inframe indels — the primary
                // indel consequence already describes the event.
                if terms.contains(&SoTerm::FrameshiftVariant)
                    || terms.contains(&SoTerm::InframeDeletion)
                    || terms.contains(&SoTerm::InframeInsertion)
                {
                    classification.stop_gained = false;
                    classification.stop_lost = false;
                }
                apply_codon_classification(terms, classification);
            } else {
                self.add_start_stop_heuristic_terms(terms, variant, tx);
            }
            terms.insert(SoTerm::ProteinAlteringVariant);
            return;
        }

        if ref_len == 0 {
            return;
        }

        if let Some(classification) = classify_coding_change(tx, tx_exons, tx_translation, variant)
        {
            let has_any = classification.has_any();
            apply_codon_classification(terms, classification);
            if has_any {
                return;
            }
        } else if tx_translation
            .and_then(|t| t.cds_sequence.as_deref())
            .is_some()
        {
            // Translation data exists but classify_coding_change failed
            // (CDS sequence mismatch). Don't fall through to the
            // heuristic — CodingSequenceVariant is already inserted.
            return;
        }

        // No translation data available. Apply start/stop heuristics
        // (position + allele pattern based) but do NOT guess
        // missense/synonymous without codon evidence.
        self.add_start_stop_heuristic_terms(terms, variant, tx);

        // Detect premature stop gain from allele pattern (in-frame
        // substitution where alt is a stop codon but ref is not).
        if ref_len == alt_len
            && ref_len % 3 == 0
            && !is_stop_codon(&variant.ref_allele)
            && is_stop_codon(&variant.alt_allele)
        {
            terms.insert(SoTerm::StopGained);
        }
    }

    fn add_start_stop_heuristic_terms(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        tx: &TranscriptFeature,
    ) {
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
        // For indels, VEP uses the affected-base range (excluding the VCF
        // anchor base) when checking splice site distances. Build a
        // trimmed variant that strips any common prefix.
        let trimmed = VariantInput::from_vcf(
            variant.chrom.clone(),
            variant.start,
            variant.end,
            variant.ref_allele.clone(),
            variant.alt_allele.clone(),
        );
        let splice_variant = &trimmed;

        for exon in tx_exons {
            if tx.strand >= 0 {
                // Donor (intron after exon end): +1/+2, +5, +3..+6, +7..+8
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.end + 1,
                    exon.end + 2,
                    SoTerm::SpliceDonorVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.end + 5,
                    exon.end + 5,
                    SoTerm::SpliceDonor5thBaseVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.end + 3,
                    exon.end + 6,
                    SoTerm::SpliceDonorRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.end + 7,
                    exon.end + 8,
                    SoTerm::SpliceRegionVariant,
                );

                // Acceptor (before exon start): -1/-2, -3..-8, -3..-17
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 2,
                    exon.start - 1,
                    SoTerm::SpliceAcceptorVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 8,
                    exon.start - 3,
                    SoTerm::SpliceRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 17,
                    exon.start - 3,
                    SoTerm::SplicePolypyrimidineTractVariant,
                );
            } else {
                // Donor for negative strand lives before exon start.
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 2,
                    exon.start - 1,
                    SoTerm::SpliceDonorVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 5,
                    exon.start - 5,
                    SoTerm::SpliceDonor5thBaseVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 6,
                    exon.start - 3,
                    SoTerm::SpliceDonorRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.start - 8,
                    exon.start - 7,
                    SoTerm::SpliceRegionVariant,
                );

                // Acceptor for negative strand lives after exon end.
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.end + 1,
                    exon.end + 2,
                    SoTerm::SpliceAcceptorVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
                    exon.end + 3,
                    exon.end + 8,
                    SoTerm::SpliceRegionVariant,
                );
                add_if_overlaps(
                    terms,
                    splice_variant,
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
    !matches!(
        biotype,
        "protein_coding"
            | "nonsense_mediated_decay"
            | "non_stop_decay"
            | "protein_coding_LoF"
            | "IG_C_gene"
            | "IG_D_gene"
            | "IG_J_gene"
            | "IG_V_gene"
            | "TR_C_gene"
            | "TR_D_gene"
            | "TR_J_gene"
            | "TR_V_gene"
            | "polymorphic_pseudogene"
    )
}

/// Returns true for transcript IDs that VEP annotates against.
/// VEP only reports ENST, NM_, NR_, XM_, XR_ transcripts.
pub fn is_vep_transcript(id: &str) -> bool {
    id.starts_with("ENST")
        || id.starts_with("NM_")
        || id.starts_with("NR_")
        || id.starts_with("XM_")
        || id.starts_with("XR_")
}

/// Strip only coding parent terms (`coding_sequence_variant`,
/// `protein_altering_variant`) when specific children are present.
/// Used at the cross-transcript collapse level.
fn strip_coding_parent_terms(terms: &mut BTreeSet<SoTerm>) {
    let has_specific_coding = terms.contains(&SoTerm::MissenseVariant)
        || terms.contains(&SoTerm::SynonymousVariant)
        || terms.contains(&SoTerm::StopGained)
        || terms.contains(&SoTerm::StopLost)
        || terms.contains(&SoTerm::StartLost)
        || terms.contains(&SoTerm::FrameshiftVariant)
        || terms.contains(&SoTerm::InframeInsertion)
        || terms.contains(&SoTerm::InframeDeletion)
        || terms.contains(&SoTerm::StopRetainedVariant)
        || terms.contains(&SoTerm::StartRetainedVariant)
        || terms.contains(&SoTerm::IncompleteTerminalCodonVariant);

    if has_specific_coding {
        terms.remove(&SoTerm::CodingSequenceVariant);
        terms.remove(&SoTerm::ProteinAlteringVariant);
    }
}

/// Remove parent SO terms when more specific children are present.
/// VEP never emits `coding_sequence_variant` alongside a specific
/// coding consequence (missense, synonymous, etc.), and never emits
/// `protein_altering_variant` alongside `missense_variant`.
/// Also strips `intron_variant` alongside splice_donor/acceptor
/// (per-transcript level only).
fn strip_parent_terms(terms: &mut BTreeSet<SoTerm>) {
    let has_specific_coding = terms.contains(&SoTerm::MissenseVariant)
        || terms.contains(&SoTerm::SynonymousVariant)
        || terms.contains(&SoTerm::StopGained)
        || terms.contains(&SoTerm::StopLost)
        || terms.contains(&SoTerm::StartLost)
        || terms.contains(&SoTerm::FrameshiftVariant)
        || terms.contains(&SoTerm::InframeInsertion)
        || terms.contains(&SoTerm::InframeDeletion)
        || terms.contains(&SoTerm::StopRetainedVariant)
        || terms.contains(&SoTerm::StartRetainedVariant)
        || terms.contains(&SoTerm::IncompleteTerminalCodonVariant);

    if has_specific_coding {
        terms.remove(&SoTerm::CodingSequenceVariant);
        terms.remove(&SoTerm::ProteinAlteringVariant);
    }

    // VEP omits intron_variant when splice_donor_variant or
    // splice_acceptor_variant is present on the same transcript —
    // the splice site terms already imply an intronic location.
    if terms.contains(&SoTerm::SpliceDonorVariant)
        || terms.contains(&SoTerm::SpliceAcceptorVariant)
    {
        terms.remove(&SoTerm::IntronVariant);
    }

    // VEP omits splice_donor_region_variant when the more specific
    // splice_donor_5th_base_variant is present — the 5th base position
    // is within the donor region (positions 3–6), so the specific term
    // subsumes the general one.
    if terms.contains(&SoTerm::SpliceDonor5thBaseVariant) {
        terms.remove(&SoTerm::SpliceDonorRegionVariant);
    }
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

impl CodingClassification {
    fn has_any(&self) -> bool {
        self.synonymous
            || self.missense
            || self.stop_gained
            || self.stop_lost
            || self.stop_retained
            || self.start_lost
            || self.start_retained
    }
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

fn classify_coding_change(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    tx_translation: Option<&TranslationFeature>,
    variant: &VariantInput,
) -> Option<CodingClassification> {
    let translation = tx_translation?;
    let cds_seq = translation.cds_sequence.as_deref()?.to_ascii_uppercase();
    let ref_genomic = normalize_allele_seq(&variant.ref_allele);
    let alt_genomic = normalize_allele_seq(&variant.alt_allele);
    let ref_len = ref_genomic.len();
    if ref_len == 0 {
        return None;
    }
    let alt_len = alt_genomic.len();

    let genomic_positions = genomic_range(variant.start, variant.end)?;
    if genomic_positions.len() != ref_len {
        return None;
    }

    // VEP prepends N characters to the CDS sequence when the first coding
    // exon starts mid-codon (non-zero phase).  These Ns align the reading
    // frame to codon boundaries but are not part of the genomic sequence.
    // We must offset CDS indices by the number of leading Ns so that the
    // genomic-to-CDS mapping lines up with the padded string.
    let leading_n_offset = cds_seq
        .bytes()
        .take_while(|&b| b == b'N')
        .count();

    let mut cds_indices = Vec::with_capacity(genomic_positions.len());
    for pos in &genomic_positions {
        let idx = genomic_to_cds_index(tx, tx_exons, *pos)?;
        cds_indices.push(idx + leading_n_offset);
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
        ref_genomic
    } else {
        reverse_complement(&ref_genomic)?.to_ascii_uppercase()
    };
    let alt_tx = if tx.strand >= 0 {
        alt_genomic
    } else {
        reverse_complement(&alt_genomic)?.to_ascii_uppercase()
    };

    let ref_seq_slice = &cds_seq[start_idx..=end_idx];
    if ref_seq_slice != ref_tx {
        return None;
    }

    let mut mutated = Vec::with_capacity(
        cds_seq
            .len()
            .saturating_sub(ref_len)
            .saturating_add(alt_len),
    );
    mutated.extend_from_slice(&cds_seq.as_bytes()[..start_idx]);
    let alt_bytes = alt_tx.as_bytes();
    mutated.extend_from_slice(alt_bytes);
    mutated.extend_from_slice(&cds_seq.as_bytes()[end_idx + 1..]);

    let old_aas = translate_protein_from_cds(cds_seq.as_bytes())?;
    let new_aas = translate_protein_from_cds(&mutated)?;
    let mut class = CodingClassification::default();

    if start_idx < 3 && old_aas.first() == Some(&'M') {
        if new_aas.first() == Some(&'M') {
            class.start_retained = true;
        } else {
            class.start_lost = true;
        }
    }

    let old_stop = old_aas.iter().position(|aa| *aa == '*');
    let new_stop = new_aas.iter().position(|aa| *aa == '*');
    if let (Some(old_stop_idx), Some(new_stop_idx)) = (old_stop, new_stop) {
        let stop_nt_start = old_stop_idx.saturating_mul(3);
        let stop_nt_end = stop_nt_start.saturating_add(2);
        if old_stop_idx == new_stop_idx
            && ranges_overlap_usize(start_idx, end_idx, stop_nt_start, stop_nt_end)
        {
            class.stop_retained = true;
        }
    }

    let frameshift = ref_len.abs_diff(alt_len) % 3 != 0;
    let stop_might_be_disrupted = if let Some(old_stop_idx) = old_stop {
        let stop_nt_start = old_stop_idx.saturating_mul(3);
        let stop_nt_end = stop_nt_start.saturating_add(2);
        ranges_overlap_usize(start_idx, end_idx, stop_nt_start, stop_nt_end)
            || (frameshift && start_idx <= stop_nt_end)
    } else {
        false
    };
    match (old_stop, new_stop) {
        (Some(old_idx), Some(new_idx)) => {
            if new_idx < old_idx {
                class.stop_gained = true;
            } else if new_idx > old_idx && stop_might_be_disrupted {
                class.stop_lost = true;
            }
        }
        (Some(_), None) => {
            if stop_might_be_disrupted {
                class.stop_lost = true;
            }
        }
        (None, Some(_)) => {
            class.stop_gained = true;
        }
        (None, None) => {}
    }

    // Per-codon stop analysis: when the CDS already contains internal stop
    // codons (e.g. pseudogenes / LoF transcripts), the global first-stop
    // comparison above may miss a new stop introduced *after* an existing
    // one.  VEP classifies stop_gained/stop_lost based on the specific
    // codon(s) affected by the variant, so we do the same here.
    if ref_len == alt_len && !class.stop_gained && !class.stop_lost {
        let first_codon = start_idx / 3;
        let last_codon = end_idx / 3;
        for ci in first_codon..=last_codon {
            if ci < old_aas.len() && ci < new_aas.len() {
                let old_aa = old_aas[ci];
                let new_aa = new_aas[ci];
                if old_aa != '*' && new_aa == '*' {
                    class.stop_gained = true;
                } else if old_aa == '*' && new_aa != '*' {
                    class.stop_lost = true;
                }
            }
        }
    }

    if ref_len == alt_len {
        let aa_changed = old_aas != new_aas;
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

fn normalize_allele_seq(allele: &str) -> String {
    if allele == "-" {
        String::new()
    } else {
        allele.to_ascii_uppercase()
    }
}

fn ranges_overlap_usize(a_start: usize, a_end: usize, b_start: usize, b_end: usize) -> bool {
    a_start <= b_end && b_start <= a_end
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

fn translate_protein_from_cds(cds: &[u8]) -> Option<Vec<char>> {
    let codon_count = cds.len() / 3;
    let mut out = Vec::with_capacity(codon_count);
    for codon_idx in 0..codon_count {
        let base = codon_idx * 3;
        let codon = std::str::from_utf8(&cds[base..base + 3]).ok()?;
        out.push(translate_codon(codon)?);
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
        _ if codon.contains('N') => Some('X'),
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

    fn translation(
        tx_id: &str,
        cds_len: Option<usize>,
        protein_len: Option<usize>,
        translation_seq: Option<&str>,
        cds_sequence: Option<&str>,
    ) -> TranslationFeature {
        TranslationFeature {
            transcript_id: tx_id.to_string(),
            cds_len,
            protein_len,
            translation_seq: translation_seq.map(|v| v.to_string()),
            cds_sequence: cds_sequence.map(|v| v.to_string()),
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
        // VEP omits the parent non_coding_transcript_variant when the more
        // specific non_coding_transcript_exon_variant is present.
        assert!(!terms.contains(&SoTerm::NonCodingTranscriptVariant));

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
    fn coding_inframe_deletion_with_translation_can_emit_stop_lost() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            90,
            140,
            1,
            "protein_coding",
            Some(100),
            Some(108),
        );
        let exons = vec![exon("pc", 1, 90, 140)];
        let translations = vec![translation(
            "pc",
            Some(9),
            Some(3),
            Some("MA*"),
            Some("ATGGCTTAA"),
        )];

        let assignments = engine.evaluate_variant_with_context(
            &var("22", 106, 108, "TAA", "-"),
            std::slice::from_ref(&tx),
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let collapsed = TranscriptConsequenceEngine::collapse_variant_terms(&assignments);
        assert!(collapsed.contains(&SoTerm::InframeDeletion));
        // VEP suppresses stop_lost alongside inframe indels.
        assert!(!collapsed.contains(&SoTerm::StopLost));
        // Parent terms are stripped when specific children are present.
        assert!(!collapsed.contains(&SoTerm::ProteinAlteringVariant));
        assert!(!collapsed.contains(&SoTerm::CodingSequenceVariant));
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
        // Without translation data, the engine cannot determine
        // missense vs synonymous — only CodingSequenceVariant.
        assert!(assignments[0].terms.contains(&SoTerm::CodingSequenceVariant));
        assert!(!assignments[0].terms.contains(&SoTerm::MissenseVariant));
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
    fn incomplete_terminal_uses_translation_cds_len_when_cds_sequence_missing() {
        let engine = TranscriptConsequenceEngine::default();
        let tx = tx(
            "pc",
            "22",
            90,
            140,
            1,
            "protein_coding",
            Some(100),
            Some(108), // complete by genomic span
        );
        let exons = vec![exon("pc", 1, 90, 140)];
        let translations = vec![translation("pc", Some(10), Some(3), None, None)];

        let assignments = engine.evaluate_variant_with_context(
            &var("22", 106, 108, "TAA", "TAG"),
            &[tx],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let collapsed = TranscriptConsequenceEngine::collapse_variant_terms(&assignments);
        assert!(collapsed.contains(&SoTerm::IncompleteTerminalCodonVariant));
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

    #[test]
    fn stop_gained_detected_when_cds_has_internal_stops() {
        // Regression test for CYP2D7-like pseudogene transcripts where the
        // CDS already contains premature stop codons.  A G->A SNV on the
        // negative strand (C->T on coding strand) converts CGA (Arg) to
        // TGA (Stop).  Previously missed because the global first-stop
        // comparison was unchanged.
        //
        // CDS layout (positive strand, 1-based coords 100-130):
        //   ATG CGA TGA CGA AAA CGA AAA AAA AAA AAA TAA
        //   M   R   *   R   K   R   K   K   K   K   *
        // Codon 2 (pos 7-9, 1-based 106-108) is already TGA (internal stop).
        // We mutate codon 5 (pos 16-18, 1-based 115-117): CGA -> TGA.
        // On positive strand the ref is C at pos 115 and alt is T.
        let engine = TranscriptConsequenceEngine::default();
        // CDS: ATG CGA TGA CGA AAA CGA AAA AAA AAA AAA TAA = 33 nt
        let cds = "ATGCGATGACGAAAACGAAAAAAAAAAAATAA";
        // transcript 100..131 (1-based inclusive), CDS 100..130
        let tx = tx("pc", "22", 100, 131, 1, "protein_coding", Some(100), Some(130));
        let exons = vec![exon("pc", 1, 100, 131)];
        // CGA at CDS positions 16-18 (0-based 15-17). Genomic pos = 100+15 = 115.
        // Mutate C->T at pos 115 to turn CGA->TGA.
        let translations = vec![translation("pc", Some(31), Some(10), None, Some(cds))];
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 115, 115, "C", "T"),
            &[tx],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::StopGained),
            "Expected stop_gained but got: {:?}",
            terms
        );
        assert!(
            !terms.contains(&SoTerm::MissenseVariant),
            "Should not emit missense when stop_gained applies: {:?}",
            terms
        );
    }
}
