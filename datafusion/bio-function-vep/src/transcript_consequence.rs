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

            if overlaps(variant.start, variant.end, tx.start, tx.end) {
                let terms = self.evaluate_transcript_overlap(variant, tx, &tx_exons);
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
            self.add_coding_terms(&mut terms, variant, tx);
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

        if self.overlaps_start_codon(variant, tx) {
            terms.insert(SoTerm::StartLost);
        }
        if self.overlaps_stop_codon(variant, tx) {
            terms.insert(SoTerm::StopLost);
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
    vec![
        SoTerm::IntergenicVariant,
        SoTerm::UpstreamGeneVariant,
        SoTerm::DownstreamGeneVariant,
        SoTerm::IntronVariant,
        SoTerm::SpliceAcceptorVariant,
        SoTerm::SpliceDonorVariant,
        SoTerm::SpliceDonor5thBaseVariant,
        SoTerm::SpliceDonorRegionVariant,
        SoTerm::SpliceRegionVariant,
        SoTerm::SplicePolypyrimidineTractVariant,
        SoTerm::FivePrimeUtrVariant,
        SoTerm::ThreePrimeUtrVariant,
        SoTerm::NonCodingTranscriptExonVariant,
        SoTerm::NonCodingTranscriptVariant,
        SoTerm::NmdTranscriptVariant,
        SoTerm::CodingTranscriptVariant,
        SoTerm::CodingSequenceVariant,
        SoTerm::FrameshiftVariant,
        SoTerm::InframeInsertion,
        SoTerm::InframeDeletion,
        SoTerm::ProteinAlteringVariant,
        SoTerm::StartLost,
        SoTerm::StopLost,
        SoTerm::SequenceVariant,
    ]
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
    fn phase2_coverage_lists_sum_to_41() {
        let implemented = implemented_phase2_terms();
        let missing = missing_phase2_terms();
        assert_eq!(implemented.len() + missing.len(), 41);
        assert!(implemented.contains(&SoTerm::FrameshiftVariant));
        assert!(missing.contains(&SoTerm::MissenseVariant));
    }
}
