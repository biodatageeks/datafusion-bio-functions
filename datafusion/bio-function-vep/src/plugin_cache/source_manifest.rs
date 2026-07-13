//! Declarative source manifest (build input) — TOML from the `vepyr-plugins`
//! catalog. Declares the table provider(s), input schema, ingest `SELECT`,
//! coordinate system, and value→CSQ mapping for one plugin. Tiering is inherited
//! from the variation cache at build time (see `plugin_cache::join`).

use datafusion::common::{DataFusionError, Result};
use serde::Deserialize;

/// Coordinate basis of the raw source. Drives the build-time `start` shift so
/// ingested coordinates match the variation cache's 1-based `start`/`end`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub enum CoordinateSystem {
    #[serde(rename = "1-based")]
    OneBased,
    #[serde(rename = "0-based-half-open")]
    ZeroBasedHalfOpen,
}

/// Table provider backing a raw source. `vcf`/`bed` delegate to bio-formats
/// (sibling crate); `csv`/`tsv`/`parquet` use builtin DataFusion providers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ProviderKind {
    Vcf,
    Csv,
    Tsv,
    Parquet,
    Bed,
}

/// Arrow value type for a declared column.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
pub enum ValueType {
    Utf8,
    Float32,
    Int32,
}

/// One field of an explicit input schema (headerless/typed TSV).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SchemaField {
    pub name: String,
    #[serde(rename = "type")]
    pub ty: ValueType,
}

/// CSV/TSV provider parameters.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CsvParams {
    #[serde(default = "default_delim")]
    pub delimiter: String,
    #[serde(default)]
    pub has_header: bool,
    #[serde(default)]
    pub comment: Option<String>,
    #[serde(default)]
    pub compression: Option<String>,
    #[serde(default)]
    pub schema: Vec<SchemaField>,
}

fn default_delim() -> String {
    "\t".into()
}

/// VCF provider parameters.
///
/// `info_fields` selects which INFO keys are materialized as columns; omit it to
/// take every INFO key declared in the VCF header. NOTE: the reader exposes INFO
/// keys as **bare, case-sensitive column names** (`AF`, `ALLELE_ID`) — not
/// `info_af` as `bio-format-vcf`'s crate docs claim — so `ingest_sql` must
/// backtick them: ``SELECT `AF` FROM plugin_x_src``.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VcfParams {
    #[serde(default)]
    pub info_fields: Option<Vec<String>>,
}

/// One raw source file registered as `plugin_<name>_src[_<part>]`.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceSpec {
    #[serde(default)]
    pub part: Option<String>,
    pub provider: ProviderKind,
    pub path: String,
    #[serde(default)]
    pub csv: Option<CsvParams>,
    #[serde(default)]
    pub vcf: Option<VcfParams>,
}

impl SourceSpec {
    /// The registered table name: `plugin_<plugin>_src` or, for multi-file
    /// sources, `plugin_<plugin>_src_<part>`.
    pub fn table_name(&self, plugin_name: &str) -> String {
        match &self.part {
            Some(p) => format!("plugin_{plugin_name}_src_{p}"),
            None => format!("plugin_{plugin_name}_src"),
        }
    }
}

/// A value column produced by `ingest_sql`, mapped to a CSQ output field.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValueColumn {
    pub column: String,
    pub csq_field: String,
    #[serde(rename = "type")]
    pub ty: ValueType,
}

/// A per-transcript match discriminator: an extra key column (produced by
/// `ingest_sql`, stored in the shard key) matched at runtime against a
/// discriminator built from a `template` over the engine-attribute namespace
/// (see `plugin_cache::template`). Empty for pure per-variant plugins.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatchColumn {
    /// Cache column name (also the column `ingest_sql` must produce).
    pub column: String,
    /// Runtime discriminator template, e.g. `"{ref_aa}{Protein_position}{alt_aa}"`.
    pub template: String,
}

/// The full source manifest for one plugin.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceManifest {
    pub plugin_name: String,
    pub coordinate_system: CoordinateSystem,
    #[serde(rename = "source")]
    pub sources: Vec<SourceSpec>,
    pub ingest_sql: String,
    #[serde(rename = "value_columns")]
    pub value_columns: Vec<ValueColumn>,
    /// Optional per-transcript match discriminators (§3.4). Empty = per-variant.
    #[serde(default, rename = "match_column")]
    pub match_columns: Vec<MatchColumn>,
}

impl SourceManifest {
    /// Match-column names, in order (the discriminator part of the key).
    pub fn match_column_names(&self) -> Vec<String> {
        self.match_columns
            .iter()
            .map(|m| m.column.clone())
            .collect()
    }
}

impl SourceManifest {
    /// Parse a source manifest from a TOML file.
    pub fn load(path: &std::path::Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| {
            DataFusionError::Execution(format!("read source manifest '{}': {e}", path.display()))
        })?;
        toml::from_str(&text).map_err(|e| {
            DataFusionError::Execution(format!("parse source manifest '{}': {e}", path.display()))
        })
    }

    /// The ingest view name: `plugin_<name>_ingest`.
    pub fn ingest_view_name(&self) -> String {
        format!("plugin_{}_ingest", self.plugin_name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NOTE: top-level scalar keys (plugin_name, coordinate_system, ingest_sql)
    // must precede any `[[table]]`/`[table]` header, else TOML absorbs them into
    // the preceding table. Keep this ordering in real manifests too.
    const CADD_LIKE: &str = r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INTEGER) AS start, CAST(pos AS INTEGER) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src_snv
"""

[[source]]
part = "snv"
provider = "csv"
path = "/tmp/snv.tsv.gz"
  [source.csv]
  delimiter = "\t"
  has_header = false
  comment = "#"
  compression = "gzip"
  schema = [
    { name = "chrom", type = "Utf8" },
    { name = "pos",   type = "Utf8" },
    { name = "ref",   type = "Utf8" },
    { name = "alt",   type = "Utf8" },
    { name = "score", type = "Utf8" },
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO_SCORE"
type = "Float32"
"##;

    #[test]
    fn parses_source_manifest() {
        let m: SourceManifest = toml::from_str(CADD_LIKE).unwrap();
        assert_eq!(m.plugin_name, "demo");
        assert_eq!(m.coordinate_system, CoordinateSystem::OneBased);
        assert_eq!(m.sources.len(), 1);
        assert_eq!(
            m.sources[0].table_name(&m.plugin_name),
            "plugin_demo_src_snv"
        );
        assert_eq!(m.sources[0].csv.as_ref().unwrap().schema.len(), 5);
        assert_eq!(m.value_columns[0].csq_field, "DEMO_SCORE");
        assert_eq!(m.value_columns[0].ty, ValueType::Float32);
        assert_eq!(m.ingest_view_name(), "plugin_demo_ingest");
    }

    #[test]
    fn single_source_table_name_has_no_part_suffix() {
        let src = SourceSpec {
            part: None,
            provider: ProviderKind::Csv,
            path: "x".into(),
            csv: None,
            vcf: None,
        };
        assert_eq!(src.table_name("cadd"), "plugin_cadd_src");
    }

    #[test]
    fn parses_vcf_source_with_selected_info_fields() {
        let src = r##"
plugin_name = "mastermind"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "/tmp/mastermind.vcf.gz"
  [source.vcf]
  info_fields = ["MMID3", "MMCNT1"]

[[value_columns]]
column = "mmid3"
csq_field = "MM_MMID3"
type = "Utf8"
"##;
        let m: SourceManifest = toml::from_str(src).unwrap();
        assert_eq!(m.sources[0].provider, ProviderKind::Vcf);
        let vcf = m.sources[0].vcf.as_ref().expect("[source.vcf] must parse");
        assert_eq!(
            vcf.info_fields.as_deref(),
            Some(["MMID3".to_string(), "MMCNT1".to_string()].as_slice())
        );
    }

    #[test]
    fn vcf_source_without_vcf_table_takes_every_info_key() {
        // No `[source.vcf]` at all: `vcf` stays `None`, which the provider reads as
        // "materialize every INFO key declared in the VCF header".
        let src = r##"
plugin_name = "mastermind"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "/tmp/mastermind.vcf.gz"

[[value_columns]]
column = "mmid3"
csq_field = "MM_MMID3"
type = "Utf8"
"##;
        let m: SourceManifest = toml::from_str(src).unwrap();
        assert_eq!(m.sources[0].provider, ProviderKind::Vcf);
        assert!(
            m.sources[0].vcf.is_none(),
            "an omitted [source.vcf] must leave `vcf` as None (= take all INFO keys)"
        );
    }

    #[test]
    fn rejects_unknown_key_inside_vcf_table() {
        // `info_field` (singular) is the plausible typo for `info_fields`. Without
        // deny_unknown_fields on VcfParams this parses happily into `info_fields: None`,
        // i.e. silently indistinguishable from omitting [source.vcf] entirely — the
        // whole INFO selection would be dropped with no diagnostic.
        let src = r##"
plugin_name = "mastermind"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "/tmp/mastermind.vcf.gz"
  [source.vcf]
  info_field = ["MMID3", "MMCNT1"]

[[value_columns]]
column = "mmid3"
csq_field = "MM_MMID3"
type = "Utf8"
"##;
        let err = toml::from_str::<SourceManifest>(src)
            .expect_err("unknown key `info_field` inside [source.vcf] must be rejected");
        assert!(
            err.to_string().contains("info_field"),
            "the error must name the offending key, got: {err}"
        );
    }

    #[test]
    fn rejects_unknown_key_instead_of_ignoring_it() {
        // `[tier]` is documented in old handoffs but does not exist in SourceManifest.
        // Before deny_unknown_fields this parsed happily and did nothing.
        let src = r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "csv"
path = "/tmp/x.tsv"

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"

[tier]
threshold = 0.01
"##;
        let err = toml::from_str::<SourceManifest>(src)
            .expect_err("unknown key [tier] must be rejected");
        assert!(
            err.to_string().contains("tier"),
            "the error must name the offending key, got: {err}"
        );
    }
}
