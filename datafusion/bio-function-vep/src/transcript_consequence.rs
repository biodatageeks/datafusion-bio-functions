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
    pub fn from_vcf(
        chrom: String,
        pos: i64,
        end: i64,
        ref_allele: String,
        alt_allele: String,
    ) -> Self {
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
            return Self {
                chrom,
                start: pos,
                end,
                ref_allele,
                alt_allele,
            };
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
    /// Mature miRNA genomic regions (mapped from cDNA attributes in raw_object_json).
    pub mature_mirna_regions: Vec<(i64, i64)>,
    // Gene metadata for per-transcript CSQ serialization.
    pub gene_stable_id: Option<String>,
    pub gene_symbol: Option<String>,
    pub gene_symbol_source: Option<String>,
    pub gene_hgnc_id: Option<String>,
    pub source: Option<String>,
    /// CDS start not found (incomplete 5' end).
    pub cds_start_nf: bool,
    /// CDS end not found (incomplete 3' end).
    pub cds_end_nf: bool,
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

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TranscriptConsequence {
    pub transcript_id: Option<String>,
    pub terms: Vec<SoTerm>,
    /// Index into `PreparedContext.transcripts` for efficient metadata lookup
    /// without cloning. `None` for non-transcript features (regulatory, motif,
    /// intergenic, etc.).
    pub transcript_idx: Option<usize>,
    /// Feature type label for CSQ serialization.
    pub feature_type: FeatureType,
    /// EXON field: "N/total" when variant overlaps an exon.
    pub exon_str: Option<String>,
    /// INTRON field: "N/total" when variant overlaps an intron.
    pub intron_str: Option<String>,
    /// cDNA_position: 1-based position within the cDNA.
    pub cdna_position: Option<String>,
    /// CDS_position: 1-based position within the CDS.
    pub cds_position: Option<String>,
    /// Protein_position: 1-based amino acid position.
    pub protein_position: Option<String>,
    /// Amino_acids: e.g. "V/I" for missense.
    pub amino_acids: Option<String>,
    /// Codons: e.g. "gTt/gAt" with VEP case convention.
    pub codons: Option<String>,
    /// Distance to transcript for upstream/downstream variants.
    pub distance: Option<i64>,
    /// FLAGS: e.g. "cds_start_NF", "cds_end_NF".
    pub flags: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FeatureType {
    Transcript,
    RegulatoryFeature,
    MotifFeature,
    #[default]
    None,
}

impl FeatureType {
    pub fn as_str(&self) -> &'static str {
        match self {
            FeatureType::Transcript => "Transcript",
            FeatureType::RegulatoryFeature => "RegulatoryFeature",
            FeatureType::MotifFeature => "MotifFeature",
            FeatureType::None => "",
        }
    }
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
            transcripts,
            exons,
            translations,
            regulatory,
            motifs,
            mirnas,
            structural,
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
        let is_ins = variant.ref_allele == "-";
        let max_dist = self.upstream_distance.max(self.downstream_distance);

        // Query the per-chromosome COITree with the variant range expanded
        // by upstream/downstream distance to catch nearby transcripts.
        if let Some(tree) = ctx.tx_trees.get(variant_chrom) {
            let query_first = (variant.start - max_dist) as i32;
            let query_last = (variant.end + max_dist) as i32;
            tree.query(query_first, query_last, |node| {
                let tx_idx = *node.metadata;
                let tx = ctx.transcripts[tx_idx];
                let tx_exons = ctx
                    .exons_by_tx
                    .get(tx.transcript_id.as_str())
                    .cloned()
                    .unwrap_or_default();
                let tx_translation = ctx
                    .translation_by_tx
                    .get(tx.transcript_id.as_str())
                    .copied();

                let variant_overlaps_tx = if is_ins {
                    // For insertions, require both flanking positions
                    // to be within the transcript (same logic as exon check).
                    variant.start > tx.start && variant.start <= tx.end
                } else {
                    overlaps(variant.start, variant.end, tx.start, tx.end)
                };
                if variant_overlaps_tx {
                    let (terms, coding_class) =
                        self.evaluate_transcript_overlap(variant, tx, &tx_exons, tx_translation);
                    if !terms.is_empty() {
                        let exon_str = which_exon_str(variant, &tx_exons);
                        let intron_str = which_intron_str(variant, &tx_exons, tx.strand);
                        let cdna_position = compute_cdna_position(variant, &tx_exons, tx.strand);
                        let (cds_position, protein_position, amino_acids, codons) =
                            if let Some(ref cc) = coding_class {
                                let cds_pos = match (cc.cds_position_start, cc.cds_position_end) {
                                    (Some(s), Some(e)) if s == e => Some(s.to_string()),
                                    (Some(s), Some(e)) => Some(format!("{s}-{e}")),
                                    _ => None,
                                };
                                let prot_pos = match (cc.protein_position_start, cc.protein_position_end) {
                                    (Some(s), Some(e)) if s == e => Some(s.to_string()),
                                    (Some(s), Some(e)) => Some(format!("{s}-{e}")),
                                    _ => None,
                                };
                                (cds_pos, prot_pos, cc.amino_acids.clone(), cc.codons.clone())
                            } else {
                                (None, None, None, None)
                            };
                        let flags = compute_flags(tx);
                        out.push(TranscriptConsequence {
                            transcript_id: Some(tx.transcript_id.clone()),
                            transcript_idx: Some(tx_idx),
                            feature_type: FeatureType::Transcript,
                            terms,
                            exon_str,
                            intron_str,
                            cdna_position,
                            cds_position,
                            protein_position,
                            amino_acids,
                            codons,
                            distance: None,
                            flags,
                        });
                    }
                } else if let Some((term, dist)) = self.upstream_downstream_term(variant, tx) {
                    out.push(TranscriptConsequence {
                        transcript_id: Some(tx.transcript_id.clone()),
                        transcript_idx: Some(tx_idx),
                        feature_type: FeatureType::Transcript,
                        terms: vec![term],
                        distance: Some(dist),
                        flags: compute_flags(tx),
                        ..Default::default()
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
                terms: vec![SoTerm::IntergenicVariant],
                ..Default::default()
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
    ) -> (Vec<SoTerm>, Option<CodingClassification>) {
        let mut terms = BTreeSet::new();
        let is_ins = variant.ref_allele == "-";
        let mut coding_class = None;

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

        // VEP adds intron_variant when the variant overlaps the intron body
        // (excluding splice site positions at the first/last 2bp).  The
        // narrower range in variant_overlaps_intron handles this correctly
        // for both purely-intronic variants and exon-spanning deletions.
        let overlaps_intron = self.variant_overlaps_intron(variant, tx_exons);
        if overlaps_intron {
            terms.insert(SoTerm::IntronVariant);
        }

        // VEP treats variants falling in frameshift introns (≤13bp) as part
        // of the surrounding coding context.  If the variant is within CDS
        // bounds and sits in a frameshift intron, emit coding_sequence_variant.
        let in_frameshift_intron = !overlaps_exon && self.in_frameshift_intron(variant, tx_exons);

        if is_non_coding_biotype(&tx.biotype) && overlaps_exon {
            // VEP: mature_miRNA_variant for miRNA transcripts where variant
            // overlaps a mature miRNA region (cDNA-mapped from attributes).
            // When within a mature miRNA region, VEP suppresses
            // non_coding_transcript_exon_variant (predicate returns 0).
            let mut in_mature_mirna = false;
            if tx.biotype == "miRNA" {
                for &(mstart, mend) in &tx.mature_mirna_regions {
                    if overlaps(variant.start, variant.end, mstart, mend) {
                        terms.insert(SoTerm::MatureMirnaVariant);
                        in_mature_mirna = true;
                        break;
                    }
                }
            }
            if !in_mature_mirna {
                terms.insert(SoTerm::NonCodingTranscriptExonVariant);
            }
        } else if in_frameshift_intron && self.overlaps_cds(variant, tx) {
            // VEP treats frameshift intron variants as coding but cannot
            // determine the precise codon change (the CDS includes the
            // intron bases).  Emit only coding_sequence_variant — no
            // specific child (frameshift/inframe) — so it survives stripping.
            terms.insert(SoTerm::CodingSequenceVariant);
        } else if overlaps_exon && self.overlaps_cds(variant, tx) {
            coding_class = self.add_coding_terms(&mut terms, variant, tx, tx_exons, tx_translation);
        } else if overlaps_exon {
            if let Some(utr_term) = self.utr_term(variant, tx) {
                terms.insert(utr_term);
            }
        }

        // VEP checks splice sites using intron boundary regions that extend
        // into exons (intron_start-4..intron_start+7 and intron_end-8..intron_end+3).
        // This naturally handles both fully-intronic variants AND exon-spanning
        // deletions that reach into splice sites.
        self.add_intron_splice_terms(&mut terms, variant, tx, tx_exons);

        if tx.biotype == "nonsense_mediated_decay" {
            terms.insert(SoTerm::NmdTranscriptVariant);
        }
        if is_non_coding_biotype(&tx.biotype) {
            // VEP omits the parent non_coding_transcript_variant when either:
            // - the more specific non_coding_transcript_exon_variant is present, or
            // - the variant falls within a mature miRNA region (VEP's
            //   within_non_coding_gene predicate checks within_mature_miRNA).
            if !terms.contains(&SoTerm::NonCodingTranscriptExonVariant)
                && !terms.contains(&SoTerm::MatureMirnaVariant)
            {
                terms.insert(SoTerm::NonCodingTranscriptVariant);
            }
        }

        // VEP omits parent SO terms when more specific children are present.
        strip_parent_terms(&mut terms);

        let mut terms_vec: Vec<SoTerm> = terms.into_iter().collect();
        terms_vec.sort_by_key(|t| t.rank());
        (terms_vec, coding_class)
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
                feature_type: FeatureType::RegulatoryFeature,
                terms: ordered,
                ..Default::default()
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
                feature_type: FeatureType::MotifFeature,
                terms: ordered,
                ..Default::default()
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
                terms: vec![SoTerm::MatureMirnaVariant],
                ..Default::default()
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
                terms: ordered,
                ..Default::default()
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

    /// Returns true when a deletion extends beyond the exons it overlaps
    /// into a non-frameshift intron.  VEP treats these as "complex indels"
    /// where CDS consequences can't be determined — only
    /// `coding_sequence_variant`.  Frameshift introns (≤12bp) are treated
    /// as part of the coding context, so they don't trigger complex indel.
    fn is_complex_indel(&self, variant: &VariantInput, tx_exons: &[&ExonFeature]) -> bool {
        // Insertions aren't complex
        if variant.ref_allele == "-" {
            return false;
        }
        let mut sorted: Vec<&ExonFeature> = tx_exons.to_vec();
        sorted.sort_by_key(|e| e.start);
        for e in &sorted {
            if !overlaps(variant.start, variant.end, e.start, e.end) {
                continue;
            }
            // Variant overlaps this exon — does it extend beyond?
            if variant.start < e.start || variant.end > e.end {
                // Check if the intron it extends into is a frameshift intron
                if self.extends_into_frameshift_intron_only(variant, &sorted, e) {
                    continue;
                }
                return true;
            }
        }
        false
    }

    /// Returns true if the variant extends only into frameshift introns
    /// (≤12bp) from the given exon.
    fn extends_into_frameshift_intron_only(
        &self,
        variant: &VariantInput,
        sorted_exons: &[&ExonFeature],
        exon: &ExonFeature,
    ) -> bool {
        for pair in sorted_exons.windows(2) {
            let intron_start = pair[0].end + 1;
            let intron_end = pair[1].start - 1;
            if intron_start > intron_end {
                continue;
            }
            // Check if variant extends into this intron from the given exon
            let adjacent = (pair[0].start == exon.start && pair[0].end == exon.end)
                || (pair[1].start == exon.start && pair[1].end == exon.end);
            if !adjacent {
                continue;
            }
            if overlaps(variant.start, variant.end, intron_start, intron_end) {
                if (intron_end - intron_start).abs() > 12 {
                    return false; // Non-frameshift intron
                }
            }
        }
        true
    }

    fn add_coding_terms(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
        tx_translation: Option<&TranslationFeature>,
    ) -> Option<CodingClassification> {
        let (ref_len, alt_len) = allele_lengths(&variant.ref_allele, &variant.alt_allele);
        terms.insert(SoTerm::CodingSequenceVariant);

        // VEP: complex indel — deletion spans exon→intron boundary.
        // VEP can't compute protein change; only coding_sequence_variant.
        if self.is_complex_indel(variant, tx_exons) {
            return None;
        }

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
                apply_codon_classification(terms, &classification);
                terms.insert(SoTerm::ProteinAlteringVariant);
                return Some(classification);
            } else {
                self.add_start_stop_heuristic_terms(terms, variant, tx);
            }
            terms.insert(SoTerm::ProteinAlteringVariant);
            return None;
        }

        if ref_len == 0 {
            return None;
        }

        if let Some(classification) = classify_coding_change(tx, tx_exons, tx_translation, variant)
        {
            let has_any = classification.has_any();
            apply_codon_classification(terms, &classification);
            if has_any {
                return Some(classification);
            }
            return Some(classification);
        } else if tx_translation
            .and_then(|t| t.cds_sequence.as_deref())
            .is_some()
        {
            // Translation data exists but classify_coding_change failed
            // (CDS sequence mismatch). Don't fall through to the
            // heuristic — CodingSequenceVariant is already inserted.
            return None;
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
        None
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
    ) -> Option<(SoTerm, i64)> {
        // For insertions, use the left flanking position (start - 1) so
        // that an insertion at the transcript boundary is detected as
        // upstream/downstream.  The insertion is between start-1 and start;
        // if start-1 is in the upstream window, the insertion is upstream.
        let check_start = if variant.ref_allele == "-" {
            variant.start.saturating_sub(1)
        } else {
            variant.start
        };
        let check_end = variant.end;
        if tx.strand >= 0 {
            let up_start = tx.start.saturating_sub(self.upstream_distance);
            let up_end = tx.start.saturating_sub(1);
            if overlaps(check_start, check_end, up_start, up_end) {
                let dist = tx.start.saturating_sub(check_end).max(0);
                return Some((SoTerm::UpstreamGeneVariant, dist));
            }
            let down_start = tx.end.saturating_add(1);
            let down_end = tx.end.saturating_add(self.downstream_distance);
            if overlaps(check_start, check_end, down_start, down_end) {
                let dist = check_start.saturating_sub(tx.end).max(0);
                return Some((SoTerm::DownstreamGeneVariant, dist));
            }
        } else {
            let up_start = tx.end.saturating_add(1);
            let up_end = tx.end.saturating_add(self.upstream_distance);
            if overlaps(check_start, check_end, up_start, up_end) {
                let dist = check_start.saturating_sub(tx.end).max(0);
                return Some((SoTerm::UpstreamGeneVariant, dist));
            }
            let down_start = tx.start.saturating_sub(self.downstream_distance);
            let down_end = tx.start.saturating_sub(1);
            if overlaps(check_start, check_end, down_start, down_end) {
                let dist = tx.start.saturating_sub(check_end).max(0);
                return Some((SoTerm::DownstreamGeneVariant, dist));
            }
        }
        None
    }

    /// Check splice site terms using intron-based iteration, matching VEP's
    /// algorithm.  VEP builds an intron boundary interval tree with extended
    /// regions (intron_start-4..intron_start+7 and intron_end-8..intron_end+3)
    /// that reach into exons.  This naturally handles both fully-intronic
    /// variants AND exon-spanning deletions.
    /// Returns true if the variant overlaps any intron body (gap between
    /// adjacent exons).  VEP's `intronic` flag effectively excludes splice
    /// site positions (the first and last 2 bases of each intron), so we
    /// use `intron_start + 2 .. intron_end - 2`.  This prevents false
    /// `intron_variant` at splice donor/acceptor sites without needing
    /// post-hoc stripping.  Frameshift introns (≤12bp span) are also
    /// skipped, matching VEP behaviour.
    fn variant_overlaps_intron(&self, variant: &VariantInput, tx_exons: &[&ExonFeature]) -> bool {
        if tx_exons.len() < 2 {
            return false;
        }
        let is_ins = variant.ref_allele == "-";
        let mut sorted: Vec<&ExonFeature> = tx_exons.to_vec();
        sorted.sort_by_key(|e| e.start);
        for pair in sorted.windows(2) {
            let intron_start = pair[0].end + 1;
            let intron_end = pair[1].start - 1;
            if intron_start > intron_end {
                continue;
            }
            // Skip frameshift introns (≤12bp)
            if (intron_end - intron_start).abs() <= 12 {
                continue;
            }
            // Narrower range excluding splice site positions
            let inner_start = intron_start + 2;
            let inner_end = intron_end - 2;
            if inner_start > inner_end {
                continue;
            }
            let hit = if is_ins {
                // For insertions, VEP requires both flanking positions
                // (start-1 and start) within the intron body — but using
                // the full intron range, not the narrower splice-excluded
                // range.  This means an insertion at intron_start (first
                // intron base) is NOT intronic (left flank is the exon).
                variant.start > intron_start && variant.start <= intron_end
            } else {
                overlaps(variant.start, variant.end, inner_start, inner_end)
            };
            if hit {
                return true;
            }
        }
        false
    }

    /// Returns true if the variant falls within a frameshift intron (≤13bp).
    /// VEP treats such variants as part of the surrounding coding context.
    fn in_frameshift_intron(&self, variant: &VariantInput, tx_exons: &[&ExonFeature]) -> bool {
        if tx_exons.len() < 2 {
            return false;
        }
        let mut sorted: Vec<&ExonFeature> = tx_exons.to_vec();
        sorted.sort_by_key(|e| e.start);
        for pair in sorted.windows(2) {
            let intron_start = pair[0].end + 1;
            let intron_end = pair[1].start - 1;
            if intron_start > intron_end {
                continue;
            }
            if (intron_end - intron_start).abs() <= 12
                && overlaps(variant.start, variant.end, intron_start, intron_end)
            {
                return true;
            }
        }
        false
    }

    fn add_intron_splice_terms(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        tx: &TranscriptFeature,
        tx_exons: &[&ExonFeature],
    ) {
        if tx_exons.len() < 2 {
            return; // Single-exon transcript has no introns.
        }

        // For indels, VEP uses the affected-base range (excluding the VCF
        // anchor base) when checking splice site distances.
        let trimmed = VariantInput::from_vcf(
            variant.chrom.clone(),
            variant.start,
            variant.end,
            variant.ref_allele.clone(),
            variant.alt_allele.clone(),
        );
        let sv = &trimmed;
        let is_ins = sv.ref_allele == "-";
        let (sv_min, sv_max) = if sv.start <= sv.end {
            (sv.start, sv.end)
        } else {
            (sv.end, sv.start)
        };

        // Derive introns from sorted exon pairs.
        let mut sorted_exons: Vec<&ExonFeature> = tx_exons.to_vec();
        sorted_exons.sort_by_key(|e| e.start);

        for pair in sorted_exons.windows(2) {
            let intron_start = pair[0].end + 1;
            let intron_end = pair[1].start - 1;
            if intron_start > intron_end {
                continue; // Degenerate / overlapping exons.
            }

            // VEP skips splice checks for "frameshift introns" — very short
            // introns (≤13bp) where abs(intron_end - intron_start) <= 12 —
            // but only when the variant actually overlaps the intron body.
            // A variant in the exonic boundary region of a frameshift intron
            // can still receive splice_region_variant.
            let is_frameshift_intron = (intron_end - intron_start).abs() <= 12;
            if is_frameshift_intron && overlaps(sv_min, sv_max, intron_start, intron_end) {
                continue;
            }

            // Quick check: does variant overlap the intron or its extended
            // boundary margins?  VEP uses two separate interval trees:
            //   tree 1: [intron_start-4, intron_end+3]  (full span)
            //   tree 2: [intron_start-4, intron_start+7] (5' boundary)
            //           [intron_end-8,   intron_end+3]   (3' boundary)
            // For small introns, the 3' boundary can extend past intron_start-4
            // and the 5' boundary can extend past intron_end+3.
            let boundary_min = (intron_start - 4).min(intron_end - 8);
            let boundary_max = (intron_end + 3).max(intron_start + 7);
            if !overlaps(sv_min, sv_max, boundary_min, boundary_max) {
                continue;
            }

            // For frameshift introns where the variant is in the exonic
            // boundary region (not overlapping the intron body), VEP only
            // checks splice_region via _intron_overlap — the first loop
            // (which checks donor/acceptor/PPT/etc.) requires overlap with
            // the full intron span tree, which small introns don't trigger.
            if is_frameshift_intron {
                // Only splice_region checks (exonic boundary of intron)
                self.add_splice_region_only(terms, sv, is_ins, intron_start, intron_end);
            } else if tx.strand >= 0 {
                self.add_splice_for_intron_positive(terms, sv, is_ins, intron_start, intron_end);
            } else {
                self.add_splice_for_intron_negative(terms, sv, is_ins, intron_start, intron_end);
            }
        }
    }

    /// Splice checks for a single intron on the positive strand.
    /// Intron coordinates are 1-based inclusive [intron_start, intron_end].
    ///
    /// VEP Perl offsets (from `_intron_effects`):
    ///   donor (start):  +0/+1 = splice_donor, +4 = 5th_base, +2..+5 = donor_region
    ///   acceptor (end): -0/-1 = splice_acceptor
    ///   splice_region:  _intron_overlap checks +2..+7, -7..-2, start-3..start-1, end+1..end+3
    ///   polypyrimidine: -16..-2 (from intron_end)
    /// For frameshift introns where the variant is in the exonic boundary
    /// region, VEP only checks splice_region (via `_intron_overlap` in the
    /// boundary loop).  Strand is irrelevant for _intron_overlap — it checks
    /// the same genomic coordinate ranges on both strands.
    fn add_splice_region_only(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        is_ins: bool,
        intron_start: i64,
        intron_end: i64,
    ) {
        if is_ins {
            let p = variant.start;
            // _intron_overlap insertion checks for splice_region
            if (intron_start + 3..=intron_start + 7).contains(&p)
                || (intron_end - 6..=intron_end - 2).contains(&p)
                || (intron_start - 2..=intron_start - 1).contains(&p)
                || (intron_end + 2..=intron_end + 3).contains(&p)
                || p == intron_start
                || p == intron_end + 1
                || p == intron_start + 2
                || p == intron_end - 1
            {
                terms.insert(SoTerm::SpliceRegionVariant);
            }
        } else {
            // splice_region (intronic): intron_start+2..intron_start+7, intron_end-7..intron_end-2
            add_if_overlaps(
                terms,
                variant,
                intron_start + 2,
                intron_start + 7,
                SoTerm::SpliceRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end - 7,
                intron_end - 2,
                SoTerm::SpliceRegionVariant,
            );
            // splice_region (exonic): intron_start-3..intron_start-1, intron_end+1..intron_end+3
            add_if_overlaps(
                terms,
                variant,
                intron_start - 3,
                intron_start - 1,
                SoTerm::SpliceRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end + 1,
                intron_end + 3,
                SoTerm::SpliceRegionVariant,
            );
        }
    }

    fn add_splice_for_intron_positive(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        is_ins: bool,
        intron_start: i64,
        intron_end: i64,
    ) {
        if is_ins {
            // VEP insertion coordinate model: overlap(P, P-1, X, Y) where P
            // is the insertion point.  This effectively checks P in [X+1, Y],
            // which is narrower than a simple point-in-range check.
            // Single-position checks overlap(P, P-1, X, X) are impossible.
            let p = variant.start;

            // Donor: overlap(P, P-1, intron_start, intron_start+1)
            //   → P in [intron_start+1, intron_start+1] = P == intron_start+1
            if p == intron_start + 1 {
                terms.insert(SoTerm::SpliceDonorVariant);
            }
            // 5th base: overlap(P, P-1, X, X) is impossible for insertions.
            // (single-position range never matches inverted coords)

            // Donor region: overlap(P, P-1, intron_start+2, intron_start+5)
            //   → P in [intron_start+3, intron_start+5]
            if (intron_start + 3..=intron_start + 5).contains(&p) {
                terms.insert(SoTerm::SpliceDonorRegionVariant);
            }

            // Acceptor: overlap(P, P-1, intron_end-1, intron_end)
            //   → P in [intron_end, intron_end] = P == intron_end
            if p == intron_end {
                terms.insert(SoTerm::SpliceAcceptorVariant);
            }

            // Polypyrimidine: VEP checks in two loops with different coord models.
            //   _overlapped_introns (normalized P-1,P): P in [intron_end-16, intron_end-1]
            //   _overlapped_introns_boundary (raw P,P-1): P in [intron_end-15, intron_end-2]
            //   Union: P in [intron_end-16, intron_end-1]
            if (intron_end - 16..=intron_end - 1).contains(&p) {
                terms.insert(SoTerm::SplicePolypyrimidineTractVariant);
            }

            // VEP _intron_overlap splice_region for insertions uses explicit
            // checks on vf_start/vf_end plus overlap on boundary ranges.
            // vf_start = P, vf_end = P-1.
            //
            // Intronic splice_region overlap ranges:
            //   overlap(P, P-1, intron_start+2, intron_start+7) → P in [intron_start+3, intron_start+7]
            //   overlap(P, P-1, intron_end-7, intron_end-2)     → P in [intron_end-6, intron_end-2]
            // Exonic splice_region overlap ranges:
            //   overlap(P, P-1, intron_start-3, intron_start-1) → P in [intron_start-2, intron_start-1]
            //   overlap(P, P-1, intron_end+1, intron_end+3)     → P in [intron_end+2, intron_end+3]
            // Explicit insertion checks from _intron_overlap:
            //   vf_start == intron_start     → P == intron_start
            //   vf_end   == intron_end       → P-1 == intron_end → P == intron_end+1
            //   vf_start == intron_start+2   → P == intron_start+2
            //   vf_end   == intron_end-2     → P-1 == intron_end-2 → P == intron_end-1
            if (intron_start + 3..=intron_start + 7).contains(&p)
                || (intron_end - 6..=intron_end - 2).contains(&p)
                || (intron_start - 2..=intron_start - 1).contains(&p)
                || (intron_end + 2..=intron_end + 3).contains(&p)
                || p == intron_start
                || p == intron_end + 1
                || p == intron_start + 2
                || p == intron_end - 1
            {
                terms.insert(SoTerm::SpliceRegionVariant);
            }
        } else {
            // Deletions/SNVs: range overlap
            // Donor: intron_start..intron_start+1
            add_if_overlaps(
                terms,
                variant,
                intron_start,
                intron_start + 1,
                SoTerm::SpliceDonorVariant,
            );
            // 5th base: intron_start+4
            add_if_overlaps(
                terms,
                variant,
                intron_start + 4,
                intron_start + 4,
                SoTerm::SpliceDonor5thBaseVariant,
            );
            // Donor region: intron_start+2..intron_start+5
            add_if_overlaps(
                terms,
                variant,
                intron_start + 2,
                intron_start + 5,
                SoTerm::SpliceDonorRegionVariant,
            );
            // splice_region (intronic): intron_start+2..intron_start+7, intron_end-7..intron_end-2
            add_if_overlaps(
                terms,
                variant,
                intron_start + 2,
                intron_start + 7,
                SoTerm::SpliceRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end - 7,
                intron_end - 2,
                SoTerm::SpliceRegionVariant,
            );
            // splice_region (exonic): intron_start-3..intron_start-1, intron_end+1..intron_end+3
            add_if_overlaps(
                terms,
                variant,
                intron_start - 3,
                intron_start - 1,
                SoTerm::SpliceRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end + 1,
                intron_end + 3,
                SoTerm::SpliceRegionVariant,
            );

            // Acceptor: intron_end-1..intron_end
            add_if_overlaps(
                terms,
                variant,
                intron_end - 1,
                intron_end,
                SoTerm::SpliceAcceptorVariant,
            );
            // Polypyrimidine: intron_end-16..intron_end-2
            add_if_overlaps(
                terms,
                variant,
                intron_end - 16,
                intron_end - 2,
                SoTerm::SplicePolypyrimidineTractVariant,
            );
        }
    }

    /// Splice checks for a single intron on the negative strand.
    /// On negative strand, the donor is at the intron END and acceptor at
    /// the intron START (genomic coordinates).
    fn add_splice_for_intron_negative(
        &self,
        terms: &mut BTreeSet<SoTerm>,
        variant: &VariantInput,
        is_ins: bool,
        intron_start: i64,
        intron_end: i64,
    ) {
        if is_ins {
            // VEP insertion coordinate model: overlap(P, P-1, X, Y) → P in [X+1, Y].
            // Single-position overlap(P, P-1, X, X) is impossible for insertions.
            let p = variant.start;

            // Donor (negative strand = intron end):
            //   overlap(P, P-1, intron_end-1, intron_end) → P == intron_end
            if p == intron_end {
                terms.insert(SoTerm::SpliceDonorVariant);
            }
            // 5th base reverse: single-position → impossible for insertions

            // Donor region reverse: overlap(P, P-1, intron_end-5, intron_end-2)
            //   → P in [intron_end-4, intron_end-2]
            if (intron_end - 4..=intron_end - 2).contains(&p) {
                terms.insert(SoTerm::SpliceDonorRegionVariant);
            }

            // Acceptor (negative strand = intron start):
            //   overlap(P, P-1, intron_start, intron_start+1) → P == intron_start+1
            if p == intron_start + 1 {
                terms.insert(SoTerm::SpliceAcceptorVariant);
            }

            // Polypyrimidine reverse: VEP checks in two loops.
            //   _overlapped_introns (normalized): P in [intron_start+2, intron_start+17]
            //   _overlapped_introns_boundary (raw): P in [intron_start+3, intron_start+16]
            //   Union: P in [intron_start+2, intron_start+17]
            if (intron_start + 2..=intron_start + 17).contains(&p) {
                terms.insert(SoTerm::SplicePolypyrimidineTractVariant);
            }

            // splice_region: same _intron_overlap logic as positive strand
            // (strand-independent in VEP)
            if (intron_start + 3..=intron_start + 7).contains(&p)
                || (intron_end - 6..=intron_end - 2).contains(&p)
                || (intron_start - 2..=intron_start - 1).contains(&p)
                || (intron_end + 2..=intron_end + 3).contains(&p)
                || p == intron_start
                || p == intron_end + 1
                || p == intron_start + 2
                || p == intron_end - 1
            {
                terms.insert(SoTerm::SpliceRegionVariant);
            }
        } else {
            // Donor (intron end, negative strand): intron_end-1..intron_end
            add_if_overlaps(
                terms,
                variant,
                intron_end - 1,
                intron_end,
                SoTerm::SpliceDonorVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end - 4,
                intron_end - 4,
                SoTerm::SpliceDonor5thBaseVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end - 5,
                intron_end - 2,
                SoTerm::SpliceDonorRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end - 7,
                intron_end - 2,
                SoTerm::SpliceRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_end + 1,
                intron_end + 3,
                SoTerm::SpliceRegionVariant,
            );

            // Acceptor (intron start, negative strand): intron_start..intron_start+1
            add_if_overlaps(
                terms,
                variant,
                intron_start,
                intron_start + 1,
                SoTerm::SpliceAcceptorVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_start + 2,
                intron_start + 7,
                SoTerm::SpliceRegionVariant,
            );
            add_if_overlaps(
                terms,
                variant,
                intron_start - 3,
                intron_start - 1,
                SoTerm::SpliceRegionVariant,
            );
            // Polypyrimidine (negative strand): intron_start+2..intron_start+16
            add_if_overlaps(
                terms,
                variant,
                intron_start + 2,
                intron_start + 16,
                SoTerm::SplicePolypyrimidineTractVariant,
            );
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
/// When `merged` is true (VEP `--merged` mode), both Ensembl and RefSeq
/// transcripts are accepted.  Default (non-merged) accepts only ENST.
pub fn is_vep_transcript(id: &str, merged: bool) -> bool {
    if merged {
        id.starts_with("ENST")
            || id.starts_with("NM_")
            || id.starts_with("NR_")
            || id.starts_with("XM_")
            || id.starts_with("XR_")
    } else {
        id.starts_with("ENST")
    }
}

/// Strip only `protein_altering_variant` when specific children are present.
/// Used at the cross-transcript collapse level.
///
/// NOTE: `coding_sequence_variant` is NOT stripped here because it may
/// legitimately come from a different transcript than the specific coding
/// term (e.g., transcript A gives frameshift, transcript B gives
/// coding_sequence_variant because the codon change can't be computed).
/// Per-transcript stripping in `strip_parent_terms` handles the case
/// where both appear on the SAME transcript.
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

    // VEP omits splice_donor_region_variant when the more specific
    // splice_donor_5th_base_variant is present — the 5th base position
    // is within the donor region (positions 3–6), so the specific term
    // subsumes the general one.
    if terms.contains(&SoTerm::SpliceDonor5thBaseVariant) {
        terms.remove(&SoTerm::SpliceDonorRegionVariant);
    }

    // VEP's splice_region function suppresses splice_region_variant when
    // any of the more specific splice site terms are present.
    if terms.contains(&SoTerm::SpliceDonorVariant)
        || terms.contains(&SoTerm::SpliceAcceptorVariant)
        || terms.contains(&SoTerm::SpliceDonorRegionVariant)
        || terms.contains(&SoTerm::SpliceDonor5thBaseVariant)
    {
        terms.remove(&SoTerm::SpliceRegionVariant);
    }

    // VEP's splice_polypyrimidine_tract_variant predicate returns 0 when
    // splice_donor or splice_acceptor is present for the same transcript.
    if terms.contains(&SoTerm::SpliceDonorVariant) || terms.contains(&SoTerm::SpliceAcceptorVariant)
    {
        terms.remove(&SoTerm::SplicePolypyrimidineTractVariant);
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

#[derive(Debug, Clone, Default)]
struct CodingClassification {
    synonymous: bool,
    missense: bool,
    stop_gained: bool,
    stop_lost: bool,
    stop_retained: bool,
    start_lost: bool,
    start_retained: bool,
    /// Amino acids: "REF/ALT" (e.g. "V/I").
    amino_acids: Option<String>,
    /// Codons with VEP case convention (e.g. "gTt/gAt").
    codons: Option<String>,
    /// 1-based CDS start position.
    cds_position_start: Option<usize>,
    /// 1-based CDS end position.
    cds_position_end: Option<usize>,
    /// 1-based protein start position.
    protein_position_start: Option<usize>,
    /// 1-based protein end position.
    protein_position_end: Option<usize>,
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

fn apply_codon_classification(terms: &mut BTreeSet<SoTerm>, c: &CodingClassification) {
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
    let leading_n_offset = cds_seq.bytes().take_while(|&b| b == b'N').count();

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

    // Populate CSQ metadata fields from the computed data.
    let first_codon = start_idx / 3;
    let last_codon = end_idx / 3;

    // CDS position (1-based)
    class.cds_position_start = Some(start_idx + 1);
    class.cds_position_end = Some(end_idx + 1);

    // Protein position (1-based)
    class.protein_position_start = Some(first_codon + 1);
    class.protein_position_end = Some(last_codon + 1);

    // Amino acids
    if first_codon < old_aas.len() {
        let ref_end = (last_codon + 1).min(old_aas.len());
        let ref_aa: String = old_aas[first_codon..ref_end].iter().collect();
        if ref_len == alt_len {
            let alt_end = (last_codon + 1).min(new_aas.len());
            if first_codon < new_aas.len() {
                let alt_aa: String = new_aas[first_codon..alt_end].iter().collect();
                if ref_aa == alt_aa {
                    class.amino_acids = Some(ref_aa);
                } else {
                    class.amino_acids = Some(format!("{ref_aa}/{alt_aa}"));
                }
            }
        } else {
            // Indels: show affected ref AAs and alt AAs
            let alt_first_codon = first_codon;
            let alt_last_codon = if mutated.len() >= 3 {
                (start_idx + alt_len).saturating_sub(1) / 3
            } else {
                0
            };
            let alt_end = (alt_last_codon + 1).min(new_aas.len());
            if alt_first_codon < new_aas.len() {
                let alt_aa: String = new_aas[alt_first_codon..alt_end].iter().collect();
                if alt_aa.is_empty() {
                    class.amino_acids = Some(format!("{ref_aa}/-"));
                } else {
                    class.amino_acids = Some(format!("{ref_aa}/{alt_aa}"));
                }
            } else {
                class.amino_acids = Some(format!("{ref_aa}/-"));
            }
        }
    }

    // Codons with VEP case convention (changed bases uppercase, context lowercase)
    if ref_len == alt_len {
        let codon_nt_start = first_codon * 3;
        let codon_nt_end = ((last_codon + 1) * 3).min(cds_seq.len());
        let codon_nt_end_alt = ((last_codon + 1) * 3).min(mutated.len());
        if codon_nt_end <= cds_seq.len() && codon_nt_end_alt <= mutated.len() {
            let ref_codon = format_codon_display(
                &cds_seq.as_bytes()[codon_nt_start..codon_nt_end],
                start_idx,
                end_idx,
                codon_nt_start,
            );
            let alt_codon = format_codon_display(
                &mutated[codon_nt_start..codon_nt_end_alt],
                start_idx,
                end_idx,
                codon_nt_start,
            );
            class.codons = Some(format!("{ref_codon}/{alt_codon}"));
        }
    }

    Some(class)
}

/// Format codon bases with VEP case convention: changed positions uppercase,
/// unchanged positions lowercase.
fn format_codon_display(bases: &[u8], changed_start: usize, changed_end: usize, codon_nt_offset: usize) -> String {
    bases.iter().enumerate().map(|(i, &b)| {
        let abs_pos = codon_nt_offset + i;
        if abs_pos >= changed_start && abs_pos <= changed_end {
            (b as char).to_ascii_uppercase()
        } else {
            (b as char).to_ascii_lowercase()
        }
    }).collect()
}

/// Determine which exon (if any) the variant overlaps.
/// Returns "N/total" string for the EXON CSQ field.
fn which_exon_str(variant: &VariantInput, tx_exons: &[&ExonFeature]) -> Option<String> {
    if tx_exons.is_empty() {
        return None;
    }
    let is_ins = variant.ref_allele == "-";
    let total = tx_exons.len();
    for exon in tx_exons {
        let hit = if is_ins {
            variant.start > exon.start && variant.start <= exon.end
        } else {
            overlaps(variant.start, variant.end, exon.start, exon.end)
        };
        if hit {
            return Some(format!("{}/{}", exon.exon_number, total));
        }
    }
    None
}

/// Determine which intron (if any) the variant overlaps.
/// Returns "N/total" string for the INTRON CSQ field.
fn which_intron_str(variant: &VariantInput, tx_exons: &[&ExonFeature], strand: i8) -> Option<String> {
    if tx_exons.len() < 2 {
        return None;
    }
    let mut sorted: Vec<&ExonFeature> = tx_exons.to_vec();
    sorted.sort_by_key(|e| e.start);
    let total_introns = sorted.len() - 1;

    for (i, pair) in sorted.windows(2).enumerate() {
        let intron_start = pair[0].end + 1;
        let intron_end = pair[1].start - 1;
        if intron_start > intron_end {
            continue;
        }
        let hit = if variant.ref_allele == "-" {
            variant.start >= intron_start && variant.start <= intron_end
        } else {
            overlaps(variant.start, variant.end, intron_start, intron_end)
        };
        if hit {
            let intron_num = if strand >= 0 {
                i + 1
            } else {
                total_introns - i
            };
            return Some(format!("{intron_num}/{total_introns}"));
        }
    }
    None
}

/// Build ordered exon segments for cDNA coordinate mapping.
/// Returns (start, end) pairs ordered by transcript direction.
fn exon_segments(tx_exons: &[&ExonFeature], strand: i8) -> Option<Vec<(i64, i64)>> {
    if tx_exons.is_empty() {
        return None;
    }
    let mut segments: Vec<(i64, i64)> = tx_exons.iter().map(|e| (e.start, e.end)).collect();
    segments.sort_by_key(|(s, _)| *s);
    if strand < 0 {
        segments.reverse();
    }
    Some(segments)
}

/// Map a genomic position to 1-based cDNA index (counting all exonic bases).
fn genomic_to_cdna_index(tx_exons: &[&ExonFeature], strand: i8, pos: i64) -> Option<usize> {
    let segments = exon_segments(tx_exons, strand)?;
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

/// Compute the cDNA_position string for a variant overlapping an exon.
fn compute_cdna_position(variant: &VariantInput, tx_exons: &[&ExonFeature], strand: i8) -> Option<String> {
    if tx_exons.is_empty() {
        return None;
    }
    let is_ins = variant.ref_allele == "-";
    // Check if variant overlaps any exon
    let in_exon = tx_exons.iter().any(|e| {
        if is_ins {
            variant.start > e.start && variant.start <= e.end
        } else {
            overlaps(variant.start, variant.end, e.start, e.end)
        }
    });
    if !in_exon {
        return None;
    }
    if is_ins {
        let a = genomic_to_cdna_index(tx_exons, strand, variant.start.saturating_sub(1));
        let b = genomic_to_cdna_index(tx_exons, strand, variant.start);
        match (a, b) {
            (Some(a), Some(b)) => {
                let lo = a.min(b);
                let hi = a.max(b);
                Some(format!("{lo}-{hi}"))
            }
            _ => None,
        }
    } else {
        let start_cdna = genomic_to_cdna_index(tx_exons, strand, variant.start)?;
        let end_cdna = genomic_to_cdna_index(tx_exons, strand, variant.end)?;
        if start_cdna == end_cdna {
            Some(start_cdna.to_string())
        } else {
            let lo = start_cdna.min(end_cdna);
            let hi = start_cdna.max(end_cdna);
            Some(format!("{lo}-{hi}"))
        }
    }
}

/// Compute FLAGS field from transcript attributes.
fn compute_flags(tx: &TranscriptFeature) -> Option<String> {
    match (tx.cds_start_nf, tx.cds_end_nf) {
        (true, true) => Some("cds_start_NF&cds_end_NF".to_string()),
        (true, false) => Some("cds_start_NF".to_string()),
        (false, true) => Some("cds_end_NF".to_string()),
        (false, false) => None,
    }
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
            mature_mirna_regions: Vec::new(),
            gene_stable_id: None,
            gene_symbol: None,
            gene_symbol_source: None,
            gene_hgnc_id: None,
            source: None,
            cds_start_nf: false,
            cds_end_nf: false,
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
        assert!(
            assignments[0]
                .terms
                .contains(&SoTerm::CodingSequenceVariant)
        );
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
        let tx = tx(
            "pc",
            "22",
            100,
            131,
            1,
            "protein_coding",
            Some(100),
            Some(130),
        );
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

    // ── Approach A tests: insertion exact-matching for splice sites ──

    #[test]
    fn insertion_at_splice_acceptor_exact_match() {
        // Transcript: 1000..2000, two exons with intron between them.
        // Exon 1: 1000..1200, Exon 2: 1400..2000
        // Intron: 1201..1399
        // Acceptor site for exon 2: positions 1398 (-2) and 1399 (-1)
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];

        // Insertion exactly at intron position 1399 (exon2.start - 1) → splice_acceptor
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1399, 1399, "-", "AAAA"),
            &[t.clone()],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::SpliceAcceptorVariant),
            "Insertion at intron pos -1 should get splice_acceptor: {:?}",
            terms
        );

        // Insertion at intron position 1397 (exon2.start - 3) → NOT splice_acceptor
        let assignments2 = engine.evaluate_variant_with_context(
            &var("22", 1397, 1397, "-", "AAAA"),
            &[t],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms2 = &assignments2[0].terms;
        assert!(
            !terms2.contains(&SoTerm::SpliceAcceptorVariant),
            "Insertion at intron pos -3 should NOT get splice_acceptor: {:?}",
            terms2
        );
    }

    #[test]
    fn insertion_splice_donor_region_uses_exact_position() {
        // Exon 1: 1000..1200, Exon 2: 1400..2000
        // Donor region for exon 1 (positive strand): +3..+6 = 1203..1206
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];

        // Insertion at position 1204 (+4): should get splice_region (from exonic
        // boundary check) but NOT splice_donor_region for insertions
        // Actually +4 is in range +3..+6 so with exact matching it DOES match.
        // The key insight: insertion at +4 should get splice_donor_region only
        // if it's an exact match to the donor region positions.
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1204, 1204, "-", "ACGCACCGCGCACCG"),
            &[t.clone()],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        // Position 1204 = exon.end + 4 = in donor region range → should match
        assert!(
            terms.contains(&SoTerm::SpliceDonorRegionVariant),
            "Insertion exactly at +4 should get splice_donor_region: {:?}",
            terms
        );

        // Insertion at +7 should get splice_region but NOT donor_region
        let assignments2 = engine.evaluate_variant_with_context(
            &var("22", 1207, 1207, "-", "ACGC"),
            &[t],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms2 = &assignments2[0].terms;
        assert!(
            terms2.contains(&SoTerm::SpliceRegionVariant),
            "Insertion at +7 should get splice_region: {:?}",
            terms2
        );
        assert!(
            !terms2.contains(&SoTerm::SpliceDonorRegionVariant),
            "Insertion at +7 should NOT get splice_donor_region: {:?}",
            terms2
        );
    }

    #[test]
    fn insertion_splice_donor_5th_base_exact_match() {
        // Exon 1: 1000..1200, intron starts at 1201
        // +5 position = 1205
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];

        // Insertion at +5 (1205): VEP's overlap(P, P-1, X, X) is impossible
        // for insertions — single-position checks never fire. So 5th_base
        // should NOT be emitted for insertions. Instead, donor_region fires
        // because overlap(P, P-1, intron_start+2, intron_start+5) →
        // P in [intron_start+3, intron_start+5] and 1205 = 1201+4 = intron_start+4.
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1205, 1205, "-", "ACGC"),
            &[t.clone()],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            !terms.contains(&SoTerm::SpliceDonor5thBaseVariant),
            "Insertion at +5 should NOT get splice_donor_5th_base (VEP single-pos impossible): {:?}",
            terms
        );
        assert!(
            terms.contains(&SoTerm::SpliceDonorRegionVariant),
            "Insertion at +5 should get splice_donor_region: {:?}",
            terms
        );

        // Insertion at +6 (1206): should NOT get splice_donor_5th_base
        let assignments2 = engine.evaluate_variant_with_context(
            &var("22", 1206, 1206, "-", "ACGC"),
            &[t],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms2 = &assignments2[0].terms;
        assert!(
            !terms2.contains(&SoTerm::SpliceDonor5thBaseVariant),
            "Insertion at +6 should NOT get splice_donor_5th_base: {:?}",
            terms2
        );
    }

    // ── Approach B tests: exon-spanning indel splice detection ──

    #[test]
    fn deletion_spanning_exon_intron_boundary_gets_splice_donor() {
        // Exon 1: 1000..1200, Exon 2: 1400..2000
        // Deletion from 1198..1202 spans exon 1 end into intron
        // Donor site: +1/+2 = 1201/1202
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];
        let cds = &"ATG".repeat(67); // enough CDS
        let translations = vec![translation("T1", Some(201), Some(66), None, Some(cds))];

        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1198, 1202, "NNNNN", "-"),
            &[t],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::SpliceDonorVariant),
            "Deletion spanning exon→intron should get splice_donor: {:?}",
            terms
        );
    }

    #[test]
    fn large_deletion_spanning_into_intron_gets_splice_acceptor() {
        // Exon 1: 1000..1200, Exon 2: 1400..2000
        // Deletion from 1380..1420 spans into intron before exon 2
        // Acceptor site for exon 2: -1/-2 = 1399/1398
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];
        let cds = &"ATG".repeat(267); // enough CDS
        let translations = vec![translation("T1", Some(801), Some(266), None, Some(cds))];

        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1380, 1420, &"N".repeat(41), "-"),
            &[t],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::SpliceAcceptorVariant),
            "Large deletion spanning into intron should get splice_acceptor: {:?}",
            terms
        );
    }

    #[test]
    fn deletion_near_tiny_intron_skips_splice_frameshift_intron() {
        // Exon 1: 1000..1200, Exon 2: 1210..2000 (tiny 9bp intron: 1201..1209)
        // VEP skips splice checks for "frameshift introns" (≤13bp).
        // The 9bp intron should be skipped — no splice terms emitted.
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1210, 2000)];
        let cds = &"ATG".repeat(67);
        let translations = vec![translation("T1", Some(201), Some(66), None, Some(cds))];

        // Deletion from 1199..1203 spans into the tiny intron
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1199, 1203, "NNNNN", "-"),
            &[t],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            !terms.contains(&SoTerm::SpliceDonorVariant),
            "Frameshift intron (≤13bp) should NOT get splice_donor: {:?}",
            terms
        );
        assert!(
            terms.contains(&SoTerm::FrameshiftVariant),
            "Should still get frameshift from CDS overlap: {:?}",
            terms
        );
    }

    #[test]
    fn variant_in_mirna_transcript_gets_mature_mirna_variant() {
        // miRNA transcript on minus strand, single exon 100..200.
        // Mature miRNA region mapped from cDNA "42-59" on minus strand:
        //   genomic_start = 200 - 59 + 1 = 142
        //   genomic_end   = 200 - 42 + 1 = 159
        let engine = TranscriptConsequenceEngine::default();
        let mut t = tx("ENST_MIRNA", "22", 100, 200, -1, "miRNA", None, None);
        t.mature_mirna_regions = vec![(142, 159)];
        let exons = vec![exon("ENST_MIRNA", 1, 100, 200)];

        // SNV at 150 — within the mature miRNA region
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 150, 150, "A", "G"),
            &[t.clone()],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::MatureMirnaVariant),
            "SNV in mature miRNA region should get mature_miRNA_variant: {:?}",
            terms
        );
        assert!(
            !terms.contains(&SoTerm::NonCodingTranscriptExonVariant),
            "VEP suppresses non_coding_transcript_exon_variant inside mature miRNA region: {:?}",
            terms
        );

        // SNV at 120 — within transcript exon but outside mature miRNA region
        let assignments2 = engine.evaluate_variant_with_context(
            &var("22", 120, 120, "A", "G"),
            &[t],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms2 = &assignments2[0].terms;
        assert!(
            !terms2.contains(&SoTerm::MatureMirnaVariant),
            "SNV outside mature miRNA region should NOT get mature_miRNA_variant: {:?}",
            terms2
        );
        assert!(
            terms2.contains(&SoTerm::NonCodingTranscriptExonVariant),
            "miRNA variant outside mature region should still get non_coding_transcript_exon_variant: {:?}",
            terms2
        );
    }

    #[test]
    fn complex_indel_spanning_exon_intron_gets_coding_sequence_variant_only() {
        // Exon 1: 1000..1050, Exon 2: 1200..1400
        // Deletion from 1045..1060 spans exon 1 into intron (non-frameshift)
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            1400,
            1,
            "protein_coding",
            Some(1000),
            Some(1400),
        );
        let exons = vec![exon("T1", 1, 1000, 1050), exon("T1", 2, 1200, 1400)];
        let cds = &"ATG".repeat(84);
        let translations = vec![translation("T1", Some(252), Some(83), None, Some(cds))];

        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1045, 1060, &"N".repeat(16), "-"),
            &[t],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::CodingSequenceVariant),
            "Complex indel should get coding_sequence_variant: {:?}",
            terms
        );
        assert!(
            !terms.contains(&SoTerm::InframeDeletion),
            "Complex indel should NOT get inframe_deletion: {:?}",
            terms
        );
        assert!(
            !terms.contains(&SoTerm::FrameshiftVariant),
            "Complex indel should NOT get frameshift_variant: {:?}",
            terms
        );
    }

    #[test]
    fn intron_variant_not_emitted_at_splice_donor_position() {
        // Variant at intron_start (splice donor +1 position) should NOT get
        // intron_variant — only splice_donor_variant.
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];

        // SNV at intron_start (1201) — splice donor +1
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1201, 1201, "A", "G"),
            &[t],
            &exons,
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::SpliceDonorVariant),
            "Should get splice_donor_variant: {:?}",
            terms
        );
        assert!(
            !terms.contains(&SoTerm::IntronVariant),
            "Splice donor position should NOT get intron_variant: {:?}",
            terms
        );
    }

    #[test]
    fn large_deletion_spanning_exon_intron_gets_intron_variant() {
        // A large deletion that starts in an exon and extends well into
        // an intron should get intron_variant (overlaps narrower range).
        let engine = TranscriptConsequenceEngine::default();
        let t = tx(
            "T1",
            "22",
            1000,
            2000,
            1,
            "protein_coding",
            Some(1000),
            Some(2000),
        );
        let exons = vec![exon("T1", 1, 1000, 1200), exon("T1", 2, 1400, 2000)];
        let cds = &"ATG".repeat(267);
        let translations = vec![translation("T1", Some(801), Some(266), None, Some(cds))];

        // Deletion from 1195..1250: overlaps exon 1, extends deep into intron
        let assignments = engine.evaluate_variant_with_context(
            &var("22", 1195, 1250, &"N".repeat(56), "-"),
            &[t],
            &exons,
            &translations,
            &[],
            &[],
            &[],
            &[],
        );
        let terms = &assignments[0].terms;
        assert!(
            terms.contains(&SoTerm::IntronVariant),
            "Large exon-spanning deletion should get intron_variant: {:?}",
            terms
        );
    }
}
