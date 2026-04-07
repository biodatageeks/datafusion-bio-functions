//! Sequence Ontology (SO) consequence terms used by VEP.

use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SoImpact {
    High,
    Moderate,
    Low,
    Modifier,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SoTerm {
    TranscriptAblation,
    SpliceAcceptorVariant,
    SpliceDonorVariant,
    StopGained,
    FrameshiftVariant,
    StopLost,
    StartLost,
    TranscriptAmplification,
    FeatureElongation,
    FeatureTruncation,
    InframeInsertion,
    InframeDeletion,
    MissenseVariant,
    ProteinAlteringVariant,
    SpliceDonor5thBaseVariant,
    SpliceRegionVariant,
    SpliceDonorRegionVariant,
    SplicePolypyrimidineTractVariant,
    IncompleteTerminalCodonVariant,
    StartRetainedVariant,
    StopRetainedVariant,
    SynonymousVariant,
    CodingSequenceVariant,
    MatureMirnaVariant,
    FivePrimeUtrVariant,
    ThreePrimeUtrVariant,
    NonCodingTranscriptExonVariant,
    IntronVariant,
    NmdTranscriptVariant,
    NonCodingTranscriptVariant,
    CodingTranscriptVariant,
    UpstreamGeneVariant,
    DownstreamGeneVariant,
    TfbsAblation,
    TfbsAmplification,
    TfBindingSiteVariant,
    RegulatoryRegionAblation,
    RegulatoryRegionAmplification,
    RegulatoryRegionVariant,
    IntergenicVariant,
    SequenceVariant,
}

pub const ALL_SO_TERMS: [SoTerm; 41] = [
    SoTerm::TranscriptAblation,
    SoTerm::SpliceAcceptorVariant,
    SoTerm::SpliceDonorVariant,
    SoTerm::StopGained,
    SoTerm::FrameshiftVariant,
    SoTerm::StopLost,
    SoTerm::StartLost,
    SoTerm::TranscriptAmplification,
    SoTerm::FeatureElongation,
    SoTerm::FeatureTruncation,
    SoTerm::InframeInsertion,
    SoTerm::InframeDeletion,
    SoTerm::MissenseVariant,
    SoTerm::ProteinAlteringVariant,
    SoTerm::SpliceDonor5thBaseVariant,
    SoTerm::SpliceRegionVariant,
    SoTerm::SpliceDonorRegionVariant,
    SoTerm::SplicePolypyrimidineTractVariant,
    SoTerm::IncompleteTerminalCodonVariant,
    SoTerm::StartRetainedVariant,
    SoTerm::StopRetainedVariant,
    SoTerm::SynonymousVariant,
    SoTerm::CodingSequenceVariant,
    SoTerm::MatureMirnaVariant,
    SoTerm::FivePrimeUtrVariant,
    SoTerm::ThreePrimeUtrVariant,
    SoTerm::NonCodingTranscriptExonVariant,
    SoTerm::IntronVariant,
    SoTerm::NmdTranscriptVariant,
    SoTerm::NonCodingTranscriptVariant,
    SoTerm::CodingTranscriptVariant,
    SoTerm::UpstreamGeneVariant,
    SoTerm::DownstreamGeneVariant,
    SoTerm::TfbsAblation,
    SoTerm::TfbsAmplification,
    SoTerm::TfBindingSiteVariant,
    SoTerm::RegulatoryRegionAblation,
    SoTerm::RegulatoryRegionAmplification,
    SoTerm::RegulatoryRegionVariant,
    SoTerm::IntergenicVariant,
    SoTerm::SequenceVariant,
];

impl SoTerm {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::TranscriptAblation => "transcript_ablation",
            Self::SpliceAcceptorVariant => "splice_acceptor_variant",
            Self::SpliceDonorVariant => "splice_donor_variant",
            Self::StopGained => "stop_gained",
            Self::FrameshiftVariant => "frameshift_variant",
            Self::StopLost => "stop_lost",
            Self::StartLost => "start_lost",
            Self::TranscriptAmplification => "transcript_amplification",
            Self::FeatureElongation => "feature_elongation",
            Self::FeatureTruncation => "feature_truncation",
            Self::InframeInsertion => "inframe_insertion",
            Self::InframeDeletion => "inframe_deletion",
            Self::MissenseVariant => "missense_variant",
            Self::ProteinAlteringVariant => "protein_altering_variant",
            Self::SpliceDonor5thBaseVariant => "splice_donor_5th_base_variant",
            Self::SpliceRegionVariant => "splice_region_variant",
            Self::SpliceDonorRegionVariant => "splice_donor_region_variant",
            Self::SplicePolypyrimidineTractVariant => "splice_polypyrimidine_tract_variant",
            Self::IncompleteTerminalCodonVariant => "incomplete_terminal_codon_variant",
            Self::StartRetainedVariant => "start_retained_variant",
            Self::StopRetainedVariant => "stop_retained_variant",
            Self::SynonymousVariant => "synonymous_variant",
            Self::CodingSequenceVariant => "coding_sequence_variant",
            Self::MatureMirnaVariant => "mature_miRNA_variant",
            Self::FivePrimeUtrVariant => "5_prime_UTR_variant",
            Self::ThreePrimeUtrVariant => "3_prime_UTR_variant",
            Self::NonCodingTranscriptExonVariant => "non_coding_transcript_exon_variant",
            Self::IntronVariant => "intron_variant",
            Self::NmdTranscriptVariant => "NMD_transcript_variant",
            Self::NonCodingTranscriptVariant => "non_coding_transcript_variant",
            Self::CodingTranscriptVariant => "coding_transcript_variant",
            Self::UpstreamGeneVariant => "upstream_gene_variant",
            Self::DownstreamGeneVariant => "downstream_gene_variant",
            Self::TfbsAblation => "TFBS_ablation",
            Self::TfbsAmplification => "TFBS_amplification",
            Self::TfBindingSiteVariant => "TF_binding_site_variant",
            Self::RegulatoryRegionAblation => "regulatory_region_ablation",
            Self::RegulatoryRegionAmplification => "regulatory_region_amplification",
            Self::RegulatoryRegionVariant => "regulatory_region_variant",
            Self::IntergenicVariant => "intergenic_variant",
            Self::SequenceVariant => "sequence_variant",
        }
    }

    pub fn from_str(value: &str) -> Option<Self> {
        let term = match value {
            "transcript_ablation" => Self::TranscriptAblation,
            "splice_acceptor_variant" => Self::SpliceAcceptorVariant,
            "splice_donor_variant" => Self::SpliceDonorVariant,
            "stop_gained" => Self::StopGained,
            "frameshift_variant" => Self::FrameshiftVariant,
            "stop_lost" => Self::StopLost,
            "start_lost" => Self::StartLost,
            "transcript_amplification" => Self::TranscriptAmplification,
            "feature_elongation" => Self::FeatureElongation,
            "feature_truncation" => Self::FeatureTruncation,
            "inframe_insertion" => Self::InframeInsertion,
            "inframe_deletion" => Self::InframeDeletion,
            "missense_variant" => Self::MissenseVariant,
            "protein_altering_variant" => Self::ProteinAlteringVariant,
            "splice_donor_5th_base_variant" => Self::SpliceDonor5thBaseVariant,
            "splice_region_variant" => Self::SpliceRegionVariant,
            "splice_donor_region_variant" => Self::SpliceDonorRegionVariant,
            "splice_polypyrimidine_tract_variant" => Self::SplicePolypyrimidineTractVariant,
            "incomplete_terminal_codon_variant" => Self::IncompleteTerminalCodonVariant,
            "start_retained_variant" => Self::StartRetainedVariant,
            "stop_retained_variant" => Self::StopRetainedVariant,
            "synonymous_variant" => Self::SynonymousVariant,
            "coding_sequence_variant" => Self::CodingSequenceVariant,
            "mature_miRNA_variant" => Self::MatureMirnaVariant,
            "5_prime_UTR_variant" => Self::FivePrimeUtrVariant,
            "3_prime_UTR_variant" => Self::ThreePrimeUtrVariant,
            "non_coding_transcript_exon_variant" => Self::NonCodingTranscriptExonVariant,
            "intron_variant" => Self::IntronVariant,
            "NMD_transcript_variant" => Self::NmdTranscriptVariant,
            "non_coding_transcript_variant" => Self::NonCodingTranscriptVariant,
            "coding_transcript_variant" => Self::CodingTranscriptVariant,
            "upstream_gene_variant" => Self::UpstreamGeneVariant,
            "downstream_gene_variant" => Self::DownstreamGeneVariant,
            "TFBS_ablation" => Self::TfbsAblation,
            "TFBS_amplification" => Self::TfbsAmplification,
            "TF_binding_site_variant" => Self::TfBindingSiteVariant,
            "regulatory_region_ablation" => Self::RegulatoryRegionAblation,
            "regulatory_region_amplification" => Self::RegulatoryRegionAmplification,
            "regulatory_region_variant" => Self::RegulatoryRegionVariant,
            "intergenic_variant" => Self::IntergenicVariant,
            "sequence_variant" => Self::SequenceVariant,
            _ => return None,
        };
        Some(term)
    }

    pub fn rank(self) -> u8 {
        match self {
            Self::TranscriptAblation => 1,
            Self::SpliceAcceptorVariant => 2,
            Self::SpliceDonorVariant => 3,
            Self::StopGained => 4,
            Self::FrameshiftVariant => 5,
            Self::StopLost => 6,
            Self::StartLost => 7,
            Self::TranscriptAmplification => 8,
            Self::FeatureElongation => 9,
            Self::FeatureTruncation => 10,
            Self::InframeInsertion => 11,
            Self::InframeDeletion => 12,
            Self::MissenseVariant => 13,
            Self::ProteinAlteringVariant => 14,
            Self::SpliceDonor5thBaseVariant => 15,
            Self::SpliceRegionVariant => 16,
            Self::SpliceDonorRegionVariant => 17,
            Self::SplicePolypyrimidineTractVariant => 18,
            Self::IncompleteTerminalCodonVariant => 19,
            Self::StartRetainedVariant => 20,
            Self::StopRetainedVariant => 21,
            Self::SynonymousVariant => 22,
            Self::CodingSequenceVariant => 23,
            Self::MatureMirnaVariant => 24,
            Self::FivePrimeUtrVariant => 25,
            Self::ThreePrimeUtrVariant => 26,
            Self::NonCodingTranscriptExonVariant => 27,
            Self::IntronVariant => 28,
            Self::NmdTranscriptVariant => 29,
            Self::NonCodingTranscriptVariant => 30,
            Self::CodingTranscriptVariant => 31,
            Self::UpstreamGeneVariant => 32,
            Self::DownstreamGeneVariant => 33,
            Self::TfbsAblation => 34,
            Self::TfbsAmplification => 35,
            Self::TfBindingSiteVariant => 36,
            Self::RegulatoryRegionAblation => 37,
            Self::RegulatoryRegionAmplification => 38,
            Self::RegulatoryRegionVariant => 39,
            Self::IntergenicVariant => 40,
            Self::SequenceVariant => 41,
        }
    }

    pub fn impact(self) -> SoImpact {
        match self {
            Self::TranscriptAblation
            | Self::SpliceAcceptorVariant
            | Self::SpliceDonorVariant
            | Self::StopGained
            | Self::FrameshiftVariant
            | Self::StopLost
            | Self::StartLost
            | Self::TranscriptAmplification
            | Self::FeatureElongation
            | Self::FeatureTruncation => SoImpact::High,

            Self::InframeInsertion
            | Self::InframeDeletion
            | Self::MissenseVariant
            | Self::ProteinAlteringVariant
            | Self::TfbsAblation => SoImpact::Moderate,

            Self::SpliceDonor5thBaseVariant
            | Self::SpliceRegionVariant
            | Self::SpliceDonorRegionVariant
            | Self::SplicePolypyrimidineTractVariant
            | Self::IncompleteTerminalCodonVariant
            | Self::StartRetainedVariant
            | Self::StopRetainedVariant
            | Self::SynonymousVariant => SoImpact::Low,

            Self::CodingSequenceVariant => SoImpact::Modifier,

            Self::MatureMirnaVariant
            | Self::FivePrimeUtrVariant
            | Self::ThreePrimeUtrVariant
            | Self::NonCodingTranscriptExonVariant
            | Self::IntronVariant
            | Self::NmdTranscriptVariant
            | Self::NonCodingTranscriptVariant
            | Self::CodingTranscriptVariant
            | Self::UpstreamGeneVariant
            | Self::DownstreamGeneVariant
            | Self::TfbsAmplification
            | Self::TfBindingSiteVariant
            | Self::RegulatoryRegionAblation
            | Self::RegulatoryRegionAmplification
            | Self::RegulatoryRegionVariant
            | Self::IntergenicVariant
            | Self::SequenceVariant => SoImpact::Modifier,
        }
    }
}

pub fn most_severe_term<'a, I>(terms: I) -> Option<SoTerm>
where
    I: IntoIterator<Item = &'a SoTerm>,
{
    terms.into_iter().min_by_key(|t| t.rank()).copied()
}

pub fn unique_sorted_terms(mut terms: Vec<SoTerm>) -> Vec<SoTerm> {
    let mut uniq = BTreeSet::new();
    for term in terms.drain(..) {
        uniq.insert(term);
    }
    let mut out: Vec<SoTerm> = uniq.into_iter().collect();
    out.sort_by_key(|t| t.rank());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_terms_count_is_41() {
        assert_eq!(ALL_SO_TERMS.len(), 41);
    }

    #[test]
    fn all_terms_roundtrip() {
        for term in ALL_SO_TERMS {
            assert_eq!(SoTerm::from_str(term.as_str()), Some(term));
        }
    }

    #[test]
    fn rank_order_selects_stop_gained_over_synonymous() {
        let terms = vec![SoTerm::SynonymousVariant, SoTerm::StopGained];
        assert_eq!(most_severe_term(&terms), Some(SoTerm::StopGained));
    }

    #[test]
    fn unique_sorted_orders_by_rank() {
        let terms = vec![
            SoTerm::SynonymousVariant,
            SoTerm::StopGained,
            SoTerm::SynonymousVariant,
        ];
        assert_eq!(
            unique_sorted_terms(terms),
            vec![SoTerm::StopGained, SoTerm::SynonymousVariant]
        );
    }
}
