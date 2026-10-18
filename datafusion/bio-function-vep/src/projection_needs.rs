//! Projection pushdown: which annotation outputs a scan's projection reads,
//! and therefore which flag groups and pipeline stages it needs.
//!
//! `annotate_vep()`'s output columns are the input VCF columns, `CSQ`,
//! `most_severe_consequence` and the typed annotation columns. A query that
//! projects a subset of them should pay only for the work that subset reads:
//! a flag group no projected column depends on is switched off, and a stage
//! whose outputs are all unprojected is skipped. Every returned column keeps
//! the value a full run gives it.
//!
//! With no projection, or when the raw `CSQ` string is read without a
//! `fields` selection, everything is needed and nothing changes.

use std::collections::HashSet;

/// Typed columns only `everything` fills; reading one keeps `everything`
/// whole. The motif columns are among them because the engine leaves them
/// null outside the `everything` CSQ layout. Verified column by column on
/// HG002 chr22 against the release-116 Ensembl cache (vepyr
/// `TestProjectionPruning`, where this table lived before it moved here).
pub(crate) const EVERYTHING_ONLY_COLUMNS: &[&str] = &[
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
    "MANE",
    "APPRIS",
    "SIFT",
    "PolyPhen",
    "DOMAINS",
    "miRNA",
    "HGVS_OFFSET",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_MID_AF",
    "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF",
    "gnomADg_AFR_AF",
    "gnomADg_AMI_AF",
    "gnomADg_AMR_AF",
    "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF",
    "gnomADg_FIN_AF",
    "gnomADg_MID_AF",
    "gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF",
    "gnomADg_SAS_AF",
];

/// Columns whose value depends on the HGVS computation, `everything` or not:
/// the two notations, the shift they report, the RefSeq misalignment offset
/// (filled only with HGVSc), and the reference-allele fields the shift also
/// feeds on RefSeq/BAM-edited transcripts.
/// `HGVS_OFFSET` is also `everything`-only, so narrowing keeps `everything`
/// for it before this set is consulted; it is listed here for the HGVS
/// pruning and the shift-prelude gate, which apply with or without it.
const HGVS_DEPENDENT_COLUMNS: &[&str] = &[
    "HGVSc",
    "HGVSp",
    "HGVS_OFFSET",
    "REFSEQ_OFFSET",
    "USED_REF",
    "GIVEN_REF",
    "BAM_EDIT",
];

const HGVS_OPTIONS: &[&str] = &[
    "hgvs",
    "hgvsc",
    "hgvsp",
    "shift_hgvs",
    "no_escape",
    "remove_hgvsp_version",
    "hgvsp_use_prediction",
];

/// Columns the co-located (existing variant) group fills. The cache-only
/// columns come from the same lookup and are kept here conservatively.
pub(crate) const COLOCATED_COLUMNS: &[&str] = &[
    "Existing_variation",
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "clin_sig_allele",
    "clinical_impact",
    "minor_allele",
    "minor_allele_freq",
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
];

const COLOCATED_OPTIONS: &[&str] = &[
    "check_existing",
    "af",
    "af_1kg",
    "af_gnomade",
    "af_gnomadg",
    "max_af",
    "pubmed",
];

/// Typed columns read only from the variation lookup output; nothing in the
/// CSQ string comes from them.
pub(crate) const CACHE_ONLY_COLUMNS: &[&str] = &[
    "clin_sig_allele",
    "clinical_impact",
    "minor_allele",
    "minor_allele_freq",
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
];

/// Typed columns computed per variant before the transcript engine runs:
/// the normalised allele, the variant class, and everything the variation
/// lookup supplies. A projection of only these (and input columns) needs no
/// transcript context.
const VARIANT_LEVEL_COLUMNS: &[&str] = &[
    "Allele",
    "VARIANT_CLASS",
    "Existing_variation",
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_MID_AF",
    "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF",
    "gnomADg_AF",
    "gnomADg_AFR_AF",
    "gnomADg_AMI_AF",
    "gnomADg_AMR_AF",
    "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF",
    "gnomADg_FIN_AF",
    "gnomADg_MID_AF",
    "gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF",
    "gnomADg_SAS_AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "clin_sig_allele",
    "clinical_impact",
    "minor_allele",
    "minor_allele_freq",
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
];

/// How much of the CSQ string a projection reads.
#[derive(Debug, Clone, PartialEq, Eq)]
enum CsqNeed {
    /// `CSQ` is not projected.
    Absent,
    /// The whole string, in the active layout.
    Full,
    /// Only these `fields` (plugin fields always trail them).
    Fields(HashSet<String>),
}

/// The annotation outputs a scan reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct OutputNeeds {
    /// Projected output column names; `None` reads every column.
    columns: Option<HashSet<String>>,
    csq: CsqNeed,
}

impl OutputNeeds {
    /// Every output is read: no projection.
    pub(crate) fn all() -> Self {
        Self {
            columns: None,
            csq: CsqNeed::Full,
        }
    }

    /// From the names of the projected output columns (`None`: no
    /// projection) and the `fields` option, which selects the CSQ fields.
    pub(crate) fn from_projection(
        projected: Option<Vec<String>>,
        csq_fields: Option<&[String]>,
    ) -> Self {
        let Some(projected) = projected else {
            return Self::all();
        };
        let columns: HashSet<String> = projected.into_iter().collect();
        let csq = if !columns.contains("CSQ") {
            CsqNeed::Absent
        } else {
            match csq_fields {
                Some(fields) => CsqNeed::Fields(fields.iter().cloned().collect()),
                None => CsqNeed::Full,
            }
        };
        Self {
            columns: Some(columns),
            csq,
        }
    }

    /// These needs plus CSQ `fields` the CSQ string consumes implicitly: the
    /// fields custom-plugin match templates read (`{HGVSc}`), which must be
    /// computed when the string carries those plugins' values.
    pub(crate) fn with_csq_fields<I: IntoIterator<Item = String>>(mut self, extra: I) -> Self {
        if let CsqNeed::Fields(fields) = &mut self.csq {
            fields.extend(extra);
        }
        self
    }

    /// Whether every output is read, so no group or stage can be dropped.
    pub(crate) fn is_all(&self) -> bool {
        self.columns.is_none() || self.csq == CsqNeed::Full
    }

    /// Whether the output column `name` is projected.
    pub(crate) fn column(&self, name: &str) -> bool {
        self.columns
            .as_ref()
            .is_none_or(|columns| columns.contains(name))
    }

    /// Whether the CSQ string is built at all.
    pub(crate) fn csq(&self) -> bool {
        self.csq != CsqNeed::Absent
    }

    fn csq_field(&self, name: &str) -> bool {
        match &self.csq {
            CsqNeed::Absent => false,
            CsqNeed::Full => true,
            CsqNeed::Fields(fields) => fields.contains(name),
        }
    }

    /// Whether the value `name` is read, as a typed column or a CSQ field.
    pub(crate) fn value(&self, name: &str) -> bool {
        self.column(name) || self.csq_field(name)
    }

    fn any_value(&self, names: &[&str]) -> bool {
        names.iter().any(|name| self.value(name))
    }

    /// Whether any output needs the transcript engine and its context:
    /// the CSQ string, the most severe consequence, or a typed column that is
    /// not computed per variant. `typed_columns` are the provider's typed
    /// annotation column names.
    pub(crate) fn transcript_engine(&self, typed_columns: &[&str]) -> bool {
        self.is_all()
            || self.csq()
            || self.column("most_severe_consequence")
            || typed_columns
                .iter()
                .any(|name| !VARIANT_LEVEL_COLUMNS.contains(name) && self.column(name))
    }

    /// Whether any typed column filled per consequence entry is projected:
    /// one that is not computed per variant. The CSQ string is built
    /// separately and does not count.
    pub(crate) fn transcript_columns(&self, typed_columns: &[&str]) -> bool {
        self.columns.is_none()
            || typed_columns
                .iter()
                .any(|name| !VARIANT_LEVEL_COLUMNS.contains(name) && self.column(name))
    }

    /// Whether the variation lookup has a consumer: the co-located sink
    /// (`check_existing`, after pruning) or a cache-only column.
    pub(crate) fn variation_lookup(&self, check_existing: bool) -> bool {
        check_existing || CACHE_ONLY_COLUMNS.iter().any(|name| self.column(name))
    }

    /// Whether SIFT/PolyPhen predictions are read.
    pub(crate) fn sift_polyphen(&self) -> bool {
        self.value("SIFT") || self.value("PolyPhen")
    }

    /// Whether protein domains are read (and with them `protein_features`).
    pub(crate) fn domains(&self) -> bool {
        self.value("DOMAINS")
    }

    /// Whether miRNA secondary-structure annotation is read.
    pub(crate) fn mirna(&self) -> bool {
        self.value("miRNA")
    }

    /// Whether any read value depends on the HGVS computation.
    pub(crate) fn hgvs_dependent(&self) -> bool {
        self.any_value(HGVS_DEPENDENT_COLUMNS)
    }

    /// Narrow an `everything` option set to the groups this projection reads.
    ///
    /// Returns `None` when nothing changes: every output is read, the options
    /// do not set `everything`, or a column only `everything` fills is read.
    /// Otherwise `everything` is replaced by the groups that change column
    /// values: HGVS when any column the HGVS computation feeds is read (the
    /// shift also gives `USED_REF` of a shifted indel), the co-located flags
    /// when a co-located column is read. Each group alone gives the same values as
    /// `everything` does. `reference_fasta_path` is kept: transcript
    /// hydration reads it for consequence calling, not only for HGVS.
    pub(crate) fn narrow_options(&self, options_json: Option<&str>) -> Option<String> {
        if self.is_all() {
            return None;
        }
        let mut options: serde_json::Map<String, serde_json::Value> =
            match serde_json::from_str(options_json?) {
                Ok(options) => options,
                Err(error) => {
                    log::debug!(
                        "projection pushdown: options not narrowed, unparsable JSON: {error}"
                    );
                    return None;
                }
            };
        if options.get("everything") != Some(&serde_json::Value::Bool(true)) {
            return None;
        }
        if self.any_value(EVERYTHING_ONLY_COLUMNS) {
            return None;
        }
        options.remove("everything");
        if self.hgvs_dependent() {
            options.insert("hgvs".into(), serde_json::Value::Bool(true));
        } else {
            for key in HGVS_OPTIONS {
                options.remove(*key);
            }
        }
        if self.any_value(COLOCATED_COLUMNS) {
            for key in COLOCATED_OPTIONS {
                options.insert((*key).into(), serde_json::Value::Bool(true));
            }
        }
        serde_json::to_string(&options).ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const KEY: [&str; 4] = ["chrom", "start", "ref", "alt"];
    const TYPED: &[&str] = &[
        "Allele",
        "Consequence",
        "SYMBOL",
        "Existing_variation",
        "AF",
        "SIFT",
        "DOMAINS",
        "clinvar_ids",
        "VARIANT_CLASS",
    ];
    const EVERYTHING: &str = r#"{"everything":true,"reference_fasta_path":"/ref.fa"}"#;

    fn needs(columns: &[&str]) -> OutputNeeds {
        let projected = KEY.iter().chain(columns).map(|s| s.to_string()).collect();
        OutputNeeds::from_projection(Some(projected), None)
    }

    fn narrowed(needs: &OutputNeeds) -> Option<serde_json::Value> {
        needs
            .narrow_options(Some(EVERYTHING))
            .map(|json| serde_json::from_str(&json).unwrap())
    }

    #[test]
    fn no_projection_reads_everything() {
        let all = OutputNeeds::from_projection(None, None);
        assert!(all.is_all());
        assert_eq!(all.narrow_options(Some(EVERYTHING)), None);
        assert!(all.transcript_engine(TYPED));
        assert!(all.variation_lookup(false));
        assert!(all.sift_polyphen() && all.domains() && all.hgvs_dependent());
    }

    #[test]
    fn raw_csq_without_fields_reads_everything() {
        let n = needs(&["CSQ"]);
        assert!(n.is_all());
        assert_eq!(n.narrow_options(Some(EVERYTHING)), None);
    }

    #[test]
    fn transcript_columns_drop_everything_and_the_lookup() {
        let n = needs(&["Consequence", "SYMBOL"]);
        assert_eq!(
            narrowed(&n).unwrap(),
            serde_json::json!({"reference_fasta_path": "/ref.fa"})
        );
        assert!(n.transcript_engine(TYPED));
        assert!(!n.variation_lookup(false));
        assert!(!n.sift_polyphen() && !n.domains() && !n.hgvs_dependent());
        assert!(!n.csq());
    }

    #[test]
    fn colocated_columns_keep_only_the_colocated_group() {
        let n = needs(&["Consequence", "Existing_variation", "AF"]);
        let json = narrowed(&n).unwrap();
        assert_eq!(json.get("everything"), None);
        assert_eq!(json.get("hgvs"), None);
        for key in COLOCATED_OPTIONS {
            assert_eq!(
                json.get(*key),
                Some(&serde_json::Value::Bool(true)),
                "{key}"
            );
        }
    }

    #[test]
    fn hgvs_columns_keep_only_the_hgvs_group() {
        let n = needs(&["HGVSc"]);
        let json = narrowed(&n).unwrap();
        assert_eq!(json.get("hgvs"), Some(&serde_json::Value::Bool(true)));
        assert_eq!(json.get("check_existing"), None);
        assert!(n.hgvs_dependent());
    }

    #[test]
    fn reference_allele_columns_keep_the_hgvs_group() {
        // USED_REF of a shifted indel is the 3'-shifted reference, which the
        // HGVS shift computes: narrowing must keep HGVS on for it.
        for column in ["USED_REF", "GIVEN_REF", "BAM_EDIT", "REFSEQ_OFFSET"] {
            let n = needs(&["Consequence", column]);
            let json = narrowed(&n).unwrap();
            assert_eq!(
                json.get("hgvs"),
                Some(&serde_json::Value::Bool(true)),
                "{column}"
            );
            assert_eq!(json.get("everything"), None, "{column}");
        }
    }

    #[test]
    fn everything_only_columns_keep_everything() {
        for column in ["SIFT", "DOMAINS", "MANE", "MOTIF_NAME", "gnomADg_NFE_AF"] {
            let n = needs(&["Consequence", column]);
            assert_eq!(n.narrow_options(Some(EVERYTHING)), None, "{column}");
        }
        let domains = needs(&["Consequence", "DOMAINS"]);
        assert!(domains.domains());
        assert!(
            !domains.sift_polyphen(),
            "DOMAINS alone reads no SIFT/PolyPhen"
        );
        assert!(!domains.hgvs_dependent());
        assert!(needs(&["USED_REF"]).hgvs_dependent());
        assert!(needs(&["SIFT"]).sift_polyphen());
        assert!(needs(&["PolyPhen"]).sift_polyphen());
    }

    #[test]
    fn options_without_everything_are_left_as_given() {
        let n = needs(&["Consequence"]);
        assert_eq!(n.narrow_options(Some(r#"{"af":true}"#)), None);
        assert_eq!(n.narrow_options(None), None);
    }

    #[test]
    fn variant_level_projection_needs_no_transcript_engine() {
        let n = needs(&["Allele", "VARIANT_CLASS", "Existing_variation", "AF"]);
        assert!(!n.transcript_engine(TYPED));
        assert!(n.variation_lookup(true));
        // A count projects one input column and needs neither.
        let count = needs(&[]);
        assert!(!count.transcript_engine(TYPED));
        assert!(!count.variation_lookup(false));
        // The most severe consequence is a transcript-engine output.
        assert!(needs(&["most_severe_consequence"]).transcript_engine(TYPED));
    }

    #[test]
    fn transcript_columns_are_the_per_entry_typed_columns() {
        assert!(OutputNeeds::all().transcript_columns(TYPED));
        assert!(needs(&["Consequence"]).transcript_columns(TYPED));
        assert!(needs(&["SIFT"]).transcript_columns(TYPED));
        assert!(!needs(&["Allele", "AF", "clinvar_ids"]).transcript_columns(TYPED));
        // The CSQ string and the most severe consequence are not typed lists.
        assert!(!needs(&["CSQ"]).transcript_columns(TYPED));
        assert!(!needs(&["most_severe_consequence"]).transcript_columns(TYPED));
    }

    #[test]
    fn cache_only_columns_keep_the_lookup() {
        let n = needs(&["clinvar_ids"]);
        assert!(n.variation_lookup(false));
        assert!(!n.transcript_engine(TYPED));
    }

    #[test]
    fn plugin_template_fields_join_the_csq_fields() {
        let fields = ["Allele".to_string()];
        let projected = KEY.iter().chain(&["CSQ"]).map(|s| s.to_string()).collect();
        let n = OutputNeeds::from_projection(Some(projected), Some(&fields));
        assert!(!n.hgvs_dependent());
        let n = n.with_csq_fields(["HGVSc".to_string()]);
        assert!(n.hgvs_dependent(), "a {{HGVSc}} template keeps HGVS");
        assert_eq!(
            narrowed(&n).unwrap().get("hgvs"),
            Some(&serde_json::Value::Bool(true))
        );
        // Without CSQ there is nothing to add to; with the whole CSQ, nothing to add.
        assert!(
            !needs(&["Consequence"])
                .with_csq_fields(["HGVSc".into()])
                .hgvs_dependent()
        );
        assert!(needs(&["CSQ"]).with_csq_fields(["HGVSc".into()]).is_all());
    }

    #[test]
    fn csq_fields_drive_the_groups() {
        let fields = ["Allele".to_string(), "HGVSc".to_string()];
        let projected = KEY.iter().chain(&["CSQ"]).map(|s| s.to_string()).collect();
        let n = OutputNeeds::from_projection(Some(projected), Some(&fields));
        assert!(!n.is_all());
        assert!(n.csq());
        assert!(n.transcript_engine(TYPED));
        assert!(n.hgvs_dependent());
        assert!(!n.sift_polyphen());
        let json = narrowed(&n).unwrap();
        assert_eq!(json.get("hgvs"), Some(&serde_json::Value::Bool(true)));
        assert_eq!(json.get("everything"), None);
    }
}
