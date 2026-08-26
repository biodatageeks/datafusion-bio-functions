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
pub struct SchemaField {
    pub name: String,
    #[serde(rename = "type")]
    pub ty: ValueType,
}

/// CSV/TSV provider parameters.
#[derive(Debug, Clone, Deserialize)]
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

/// One raw source file registered as `plugin_<name>_src[_<part>]`.
#[derive(Debug, Clone, Deserialize)]
pub struct SourceSpec {
    #[serde(default)]
    pub part: Option<String>,
    pub provider: ProviderKind,
    pub path: String,
    #[serde(default)]
    pub csv: Option<CsvParams>,
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
pub struct ValueColumn {
    pub column: String,
    pub csq_field: String,
    #[serde(rename = "type")]
    pub ty: ValueType,
    /// Text for this field's `##<CSQ_FIELD>=<description>` VCF header line,
    /// mirroring what an Ensembl plugin returns from `get_header_info()`.
    #[serde(default)]
    pub description: Option<String>,
}

/// A per-transcript match discriminator: an extra key column (produced by
/// `ingest_sql`, stored in the shard key) matched at runtime against a
/// discriminator built from a `template` over the engine-attribute namespace
/// (see `plugin_cache::template`). Empty for pure per-variant plugins.
#[derive(Debug, Clone, Deserialize)]
pub struct MatchColumn {
    /// Cache column name (also the column `ingest_sql` must produce).
    pub column: String,
    /// Runtime discriminator template, e.g. `"{ref_aa}{Protein_position}{alt_aa}"`.
    pub template: String,
}

/// The full source manifest for one plugin.
#[derive(Debug, Clone, Deserialize)]
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
    /// How this plugin's Ensembl implementation compares a variant to a data
    /// row. Set `minimised` only for plugins that call
    /// `get_matched_variant_alleles()` (CADD, AlphaMissense); the default
    /// `exact` matches SpliceAI's and dbNSFP's verbatim comparison.
    #[serde(default)]
    pub allele_match: crate::plugin_cache::cache_manifest::AlleleMatch,
    /// Position of this plugin's field block in the CSQ string (lower first).
    #[serde(default = "crate::plugin_cache::cache_manifest::default_csq_rank_pub")]
    pub csq_rank: u32,
    /// Field order within this plugin's block.
    #[serde(default)]
    pub field_order: crate::plugin_cache::cache_manifest::FieldOrder,
    /// Skip the build-time keep-first dedup pass (`dedup::dedup_keep_first`)
    /// when the source is structurally guaranteed to never emit two rows
    /// sharing a runtime probe key `(start, allele_string, <match cols>)` —
    /// e.g. SpliceAI's masked release (one prediction per variant) or CADD's
    /// SNV+indel files (disjoint allele-string shapes, so they can't collide
    /// with each other or within themselves). Dedup's `HashSet<String>` costs
    /// one heap-allocated key per row and is the dominant memory cost on the
    /// largest chromosomes; skipping it for sources that provably have no
    /// duplicates avoids that cost entirely. Do NOT set this for sources that
    /// can legitimately repeat a key (e.g. AlphaMissense's overlapping
    /// UniProt entries) — the manifest author must justify this per-plugin,
    /// it is not a generic performance knob.
    #[serde(default)]
    pub assume_unique: bool,
}

/// Standard VCF meta-information keys and provenance keys owned by this sink.
/// A plugin field using one would either produce an invalid declaration (for
/// structured keys such as `INFO`) or replace file-owned scalar metadata (for
/// keys such as the mandatory `fileformat`). Arbitrary custom keys are valid;
/// their structured declarations are disambiguated in `vcf_sink` instead.
const RESERVED_VCF_META_KEYS: &[&str] = &[
    "fileformat",
    "fileDate",
    "source",
    "reference",
    "phasing",
    "assembly",
    "pedigreeDB",
    "INFO",
    "FILTER",
    "FORMAT",
    "ALT",
    "contig",
    "SAMPLE",
    "PEDIGREE",
    "META",
    "VEP",
    "VEP-command-line",
];

pub(crate) fn validate_csq_field_name(plugin_name: &str, field: &str) -> Result<()> {
    let mut chars = field.chars();
    let valid_first = chars
        .next()
        .is_some_and(|c| c.is_ascii_alphanumeric() || c == '_');
    let valid_rest = chars.all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '+' | '-'));
    if !valid_first || !valid_rest {
        return Err(DataFusionError::Execution(format!(
            "plugin '{plugin_name}' has invalid CSQ field name {field:?}; field names must start \
             with an ASCII alphanumeric character or '_' and contain only ASCII alphanumeric \
             characters, '_', '.', '+', or '-'"
        )));
    }

    let engine_key = env!("CARGO_PKG_NAME");
    if RESERVED_VCF_META_KEYS.contains(&field)
        || field == engine_key
        || field
            .strip_suffix("-command-line")
            .is_some_and(|prefix| prefix == engine_key)
    {
        return Err(DataFusionError::Execution(format!(
            "plugin '{plugin_name}' CSQ field '{field}' collides with a reserved VCF \
             meta-information key"
        )));
    }
    Ok(())
}

impl SourceManifest {
    /// Validate constraints shared by build-time source manifests and runtime
    /// built-cache manifests.
    pub fn validate(&self) -> Result<()> {
        for value in &self.value_columns {
            validate_csq_field_name(&self.plugin_name, &value.csq_field)?;
        }
        Ok(())
    }

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
        let manifest: Self = toml::from_str(&text).map_err(|e| {
            DataFusionError::Execution(format!("parse source manifest '{}': {e}", path.display()))
        })?;
        manifest.validate()?;
        Ok(manifest)
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
    fn rejects_reserved_vcf_keys_but_accepts_arbitrary_custom_fields() {
        let base: SourceManifest = toml::from_str(CADD_LIKE).unwrap();
        for reserved in RESERVED_VCF_META_KEYS.iter().copied().chain([
            env!("CARGO_PKG_NAME"),
            concat!(env!("CARGO_PKG_NAME"), "-command-line"),
        ]) {
            let mut manifest = base.clone();
            manifest.value_columns[0].csq_field = reserved.to_string();
            let error = manifest.validate().unwrap_err().to_string();
            assert!(
                error.contains("reserved VCF meta-information key"),
                "{error}"
            );
        }

        let mut manifest = base;
        manifest.value_columns[0].csq_field = "SCORE".to_string();
        manifest.value_columns[0].description = Some("<threshold>".to_string());
        manifest.validate().unwrap();
    }

    #[test]
    fn enforces_vcf_safe_csq_field_identifiers() {
        let base: SourceManifest = toml::from_str(CADD_LIKE).unwrap();
        for unsafe_name in [
            "",
            " ",
            "CADD RAW",
            "A|B",
            "A=B",
            "A\nB",
            "A\rB",
            "A\tB",
            "\u{1}SCORE",
            ".SCORE",
            "SCORE/RAW",
            "café",
        ] {
            let mut manifest = base.clone();
            manifest.value_columns[0].csq_field = unsafe_name.to_string();
            let error = manifest.validate().unwrap_err().to_string();
            assert!(error.contains("invalid CSQ field name"), "{error}");
            assert!(
                !error.contains('\n'),
                "unsafe name must be escaped: {error:?}"
            );
        }

        for safe_name in ["CADD_RAW", "GERP++_RS", "SCORE.1", "CUSTOM-FIELD", "1000G"] {
            let mut manifest = base.clone();
            manifest.value_columns[0].csq_field = safe_name.to_string();
            manifest.validate().unwrap();
        }
    }

    #[test]
    fn single_source_table_name_has_no_part_suffix() {
        let src = SourceSpec {
            part: None,
            provider: ProviderKind::Csv,
            path: "x".into(),
            csv: None,
        };
        assert_eq!(src.table_name("cadd"), "plugin_cadd_src");
    }
}
