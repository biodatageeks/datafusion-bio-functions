//! VEP plugin registry — defines supported external annotation plugins,
//! their output schemas, and runtime configuration.

use std::path::{Path, PathBuf};

use datafusion::arrow::datatypes::{DataType, Field};

/// Supported VEP plugin types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PluginKind {
    ClinVar,
    Cadd,
    SpliceAI,
    AlphaMissense,
    DbNSFP,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PluginSourceKind {
    ClinVar,
    Cadd,
    SpliceAI,
    AlphaMissense,
    DbNSFP,
}

impl PluginKind {
    /// Parse plugin name from a string (case-insensitive).
    pub fn from_name(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "clinvar" => Some(Self::ClinVar),
            "cadd" => Some(Self::Cadd),
            "spliceai" => Some(Self::SpliceAI),
            "alphamissense" => Some(Self::AlphaMissense),
            "dbnsfp" => Some(Self::DbNSFP),
            _ => None,
        }
    }

    /// Cache subdirectory name for this logical plugin.
    pub fn dir_name(&self) -> &'static str {
        match self {
            Self::ClinVar => "clinvar",
            Self::Cadd => "cadd",
            Self::SpliceAI => "spliceai",
            Self::AlphaMissense => "alphamissense",
            Self::DbNSFP => "dbnsfp",
        }
    }

    /// Output columns specific to this plugin (excluding join key columns).
    pub fn output_fields(&self) -> Vec<Field> {
        match self {
            Self::ClinVar => vec![
                Field::new("ClinVar", DataType::Utf8, true),
                Field::new("ClinVar_CLNSIG", DataType::Utf8, true),
                Field::new("ClinVar_CLNREVSTAT", DataType::Utf8, true),
                Field::new("ClinVar_CLNDN", DataType::Utf8, true),
                Field::new("ClinVar_CLNVC", DataType::Utf8, true),
                Field::new("ClinVar_CLNVI", DataType::Utf8, true),
            ],
            Self::Cadd => vec![
                Field::new("raw_score", DataType::Float32, true),
                Field::new("phred_score", DataType::Float32, true),
            ],
            Self::SpliceAI => vec![
                Field::new("symbol", DataType::Utf8, true),
                Field::new("ds_ag", DataType::Float32, true),
                Field::new("ds_al", DataType::Float32, true),
                Field::new("ds_dg", DataType::Float32, true),
                Field::new("ds_dl", DataType::Float32, true),
                Field::new("dp_ag", DataType::Int32, true),
                Field::new("dp_al", DataType::Int32, true),
                Field::new("dp_dg", DataType::Int32, true),
                Field::new("dp_dl", DataType::Int32, true),
            ],
            Self::AlphaMissense => vec![
                Field::new("genome", DataType::Utf8, true),
                Field::new("uniprot_id", DataType::Utf8, true),
                Field::new("transcript_id", DataType::Utf8, true),
                Field::new("protein_variant", DataType::Utf8, true),
                Field::new("am_pathogenicity", DataType::Float32, true),
                Field::new("am_class", DataType::Utf8, true),
            ],
            Self::DbNSFP => vec![
                Field::new("sift4g_score", DataType::Utf8, true),
                Field::new("sift4g_pred", DataType::Utf8, true),
                Field::new("polyphen2_hdiv_score", DataType::Utf8, true),
                Field::new("polyphen2_hvar_score", DataType::Utf8, true),
                Field::new("lrt_score", DataType::Utf8, true),
                Field::new("lrt_pred", DataType::Utf8, true),
                Field::new("mutationtaster_score", DataType::Utf8, true),
                Field::new("mutationtaster_pred", DataType::Utf8, true),
                Field::new("fathmm_score", DataType::Utf8, true),
                Field::new("fathmm_pred", DataType::Utf8, true),
                Field::new("provean_score", DataType::Utf8, true),
                Field::new("provean_pred", DataType::Utf8, true),
                Field::new("vest4_score", DataType::Utf8, true),
                Field::new("metasvm_score", DataType::Utf8, true),
                Field::new("metasvm_pred", DataType::Utf8, true),
                Field::new("metalr_score", DataType::Utf8, true),
                Field::new("metalr_pred", DataType::Utf8, true),
                Field::new("revel_score", DataType::Utf8, true),
                Field::new("gerp_rs", DataType::Utf8, true),
                Field::new("phylop100way", DataType::Utf8, true),
                Field::new("phylop30way", DataType::Utf8, true),
                Field::new("phastcons100way", DataType::Utf8, true),
                Field::new("phastcons30way", DataType::Utf8, true),
                Field::new("siphy_29way", DataType::Utf8, true),
                Field::new("cadd_raw", DataType::Utf8, true),
                Field::new("cadd_phred", DataType::Utf8, true),
            ],
        }
    }

    /// Column names that form the join key (shared across all plugins).
    pub fn join_key_columns() -> &'static [&'static str] {
        &["chrom", "pos", "ref", "alt"]
    }

    /// CSQ field names appended when this plugin is enabled.
    pub fn csq_field_names(&self) -> Vec<&'static str> {
        match self {
            Self::ClinVar => vec![
                "ClinVar",
                "ClinVar_CLNSIG",
                "ClinVar_CLNREVSTAT",
                "ClinVar_CLNDN",
                "ClinVar_CLNVC",
                "ClinVar_CLNVI",
            ],
            Self::Cadd => vec!["raw_score", "phred_score"],
            Self::SpliceAI => vec![
                "symbol", "ds_ag", "ds_al", "ds_dg", "ds_dl", "dp_ag", "dp_al", "dp_dg", "dp_dl",
            ],
            Self::AlphaMissense => vec![
                "genome",
                "uniprot_id",
                "transcript_id",
                "protein_variant",
                "am_pathogenicity",
                "am_class",
            ],
            Self::DbNSFP => vec![
                "sift4g_score",
                "sift4g_pred",
                "polyphen2_hdiv_score",
                "polyphen2_hvar_score",
                "lrt_score",
                "lrt_pred",
                "mutationtaster_score",
                "mutationtaster_pred",
                "fathmm_score",
                "fathmm_pred",
                "provean_score",
                "provean_pred",
                "vest4_score",
                "metasvm_score",
                "metasvm_pred",
                "metalr_score",
                "metalr_pred",
                "revel_score",
                "gerp_rs",
                "phylop100way",
                "phylop30way",
                "phastcons100way",
                "phastcons30way",
                "siphy_29way",
                "cadd_raw",
                "cadd_phred",
            ],
        }
    }
}

impl PluginSourceKind {
    pub fn dir_name(&self) -> &'static str {
        match self {
            Self::ClinVar => "clinvar",
            Self::Cadd => "cadd",
            Self::SpliceAI => "spliceai",
            Self::AlphaMissense => "alphamissense",
            Self::DbNSFP => "dbnsfp",
        }
    }

    pub fn plugin_kind(&self) -> PluginKind {
        match self {
            Self::ClinVar => PluginKind::ClinVar,
            Self::Cadd => PluginKind::Cadd,
            Self::SpliceAI => PluginKind::SpliceAI,
            Self::AlphaMissense => PluginKind::AlphaMissense,
            Self::DbNSFP => PluginKind::DbNSFP,
        }
    }
}

/// Per-plugin configuration resolved at annotation time.
#[derive(Debug, Clone)]
pub struct PluginConfig {
    pub kind: PluginKind,
    pub source_dirs: Vec<PluginSourceConfig>,
}

#[derive(Debug, Clone)]
pub struct PluginSourceConfig {
    pub kind: PluginSourceKind,
    pub source_dir: PathBuf,
}

/// Resolved set of active plugins for an annotation run.
#[derive(Debug, Clone, Default)]
pub struct ActivePlugins {
    pub configs: Vec<PluginConfig>,
}

impl ActivePlugins {
    /// Discover active plugins from a cache root directory.
    /// Only includes plugins whose subdirectory exists.
    pub fn discover(plugins_dir: &Path) -> Self {
        let all_kinds = [
            PluginKind::ClinVar,
            PluginKind::Cadd,
            PluginKind::SpliceAI,
            PluginKind::AlphaMissense,
            PluginKind::DbNSFP,
        ];
        let configs: Vec<PluginConfig> = all_kinds
            .iter()
            .filter_map(|kind| Self::config_for_kind(*kind, plugins_dir))
            .collect();
        Self { configs }
    }

    /// Create from an explicit list of plugin names and a cache-root directory.
    pub fn from_names(names: &[String], plugins_dir: &Path) -> Self {
        let configs: Vec<PluginConfig> = names
            .iter()
            .filter_map(|name| {
                let kind = PluginKind::from_name(name)?;
                Self::config_for_kind(kind, plugins_dir)
            })
            .collect();
        Self { configs }
    }

    fn config_for_kind(kind: PluginKind, plugins_dir: &Path) -> Option<PluginConfig> {
        let source_kinds: &[PluginSourceKind] = match kind {
            PluginKind::ClinVar => &[PluginSourceKind::ClinVar],
            PluginKind::Cadd => &[PluginSourceKind::Cadd],
            PluginKind::SpliceAI => &[PluginSourceKind::SpliceAI],
            PluginKind::AlphaMissense => &[PluginSourceKind::AlphaMissense],
            PluginKind::DbNSFP => &[PluginSourceKind::DbNSFP],
        };

        let source_dirs: Option<Vec<PluginSourceConfig>> = source_kinds
            .iter()
            .map(|source_kind| {
                let dir = plugins_dir.join(source_kind.dir_name());
                dir.is_dir().then(|| PluginSourceConfig {
                    kind: *source_kind,
                    source_dir: dir,
                })
            })
            .collect();

        source_dirs.map(|source_dirs| PluginConfig { kind, source_dirs })
    }

    pub fn is_empty(&self) -> bool {
        self.configs.is_empty()
    }

    pub fn has_kind(&self, kind: PluginKind) -> bool {
        self.configs.iter().any(|cfg| cfg.kind == kind)
    }

    /// Collect all output fields from active plugins.
    pub fn output_fields(&self) -> Vec<Field> {
        self.configs
            .iter()
            .flat_map(|cfg| cfg.kind.output_fields())
            .collect()
    }

    /// Collect CSQ field names from active plugins in output order.
    pub fn csq_field_names(&self) -> Vec<&'static str> {
        self.configs
            .iter()
            .flat_map(|cfg| cfg.kind.csq_field_names())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{ActivePlugins, PluginKind};

    #[test]
    fn plugin_kind_parses_case_insensitively() {
        assert_eq!(PluginKind::from_name("ClinVar"), Some(PluginKind::ClinVar));
        assert_eq!(
            PluginKind::from_name("SPLICEAI"),
            Some(PluginKind::SpliceAI)
        );
        assert_eq!(PluginKind::from_name("missing"), None);
    }

    #[test]
    fn active_plugins_from_names_preserves_declared_order() {
        let temp = tempfile::tempdir().expect("tempdir");
        std::fs::create_dir_all(temp.path().join("spliceai")).expect("spliceai dir");
        std::fs::create_dir_all(temp.path().join("clinvar")).expect("clinvar dir");

        let plugins = ActivePlugins::from_names(
            &["spliceai".to_string(), "clinvar".to_string()],
            temp.path(),
        );

        assert_eq!(plugins.configs.len(), 2);
        assert_eq!(plugins.configs[0].kind, PluginKind::SpliceAI);
        assert_eq!(plugins.configs[1].kind, PluginKind::ClinVar);
        assert_eq!(plugins.csq_field_names()[..3], ["symbol", "ds_ag", "ds_al"]);
    }

    #[test]
    fn cadd_is_discovered_from_single_directory() {
        let temp = tempfile::tempdir().expect("tempdir");
        std::fs::create_dir_all(temp.path().join("cadd")).expect("cadd dir");

        let plugins = ActivePlugins::discover(temp.path());
        assert_eq!(plugins.configs.len(), 1);
        assert_eq!(plugins.configs[0].kind, PluginKind::Cadd);
        assert_eq!(plugins.configs[0].source_dirs.len(), 1);
    }

    #[test]
    fn active_plugins_has_kind_matches_enabled_plugins() {
        let temp = tempfile::tempdir().expect("tempdir");
        std::fs::create_dir_all(temp.path().join("spliceai")).expect("spliceai dir");
        std::fs::create_dir_all(temp.path().join("clinvar")).expect("clinvar dir");

        let plugins = ActivePlugins::from_names(
            &["spliceai".to_string(), "clinvar".to_string()],
            temp.path(),
        );

        assert!(plugins.has_kind(PluginKind::ClinVar));
        assert!(plugins.has_kind(PluginKind::SpliceAI));
        assert!(!plugins.has_kind(PluginKind::Cadd));
    }
}
