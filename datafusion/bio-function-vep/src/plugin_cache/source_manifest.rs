//! Declarative source manifest (build input) — TOML from the `vepyr-plugins`
//! catalog. Declares the table provider(s), input schema, ingest `SELECT`,
//! coordinate system, and value→CSQ mapping for one plugin. Tiering is inherited
//! from the variation cache at build time (see `plugin_cache::join`).

use datafusion::common::{DataFusionError, Result};
use serde::Deserialize;

/// Coordinate basis of the raw source. Drives the build-time `start` shift so
/// ingested coordinates match the variation cache's 1-based `start`/`end`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
pub enum CoordinateSystem {
    #[serde(rename = "1-based")]
    OneBased,
    #[serde(rename = "0-based-half-open")]
    ZeroBasedHalfOpen,
}

impl CoordinateSystem {
    /// The manifest spelling, for diagnostics.
    pub const fn as_manifest_str(self) -> &'static str {
        match self {
            Self::OneBased => "1-based",
            Self::ZeroBasedHalfOpen => "0-based-half-open",
        }
    }
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

impl ProviderKind {
    /// The manifest spelling, for diagnostics.
    pub const fn as_manifest_str(self) -> &'static str {
        match self {
            Self::Vcf => "vcf",
            Self::Csv => "csv",
            Self::Tsv => "tsv",
            Self::Parquet => "parquet",
            Self::Bed => "bed",
        }
    }
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
///
/// The reader's core columns are `chrom`, `start`, `end`, `id`, `ref`, `alt`, `qual`,
/// `filter` — there is NO `pos`. VCF POS is exposed as `start` (1-based, matching the
/// cache's coordinate system). So `ingest_sql` must select `start`/`end`, not `pos`.
///
/// # `end` is not always `start`
///
/// `end` equals `start` only for a plain single-ALT ACGT SNV. Otherwise the reader
/// reports the variant end — `POS + len(REF) - 1` for indels/MNVs (and the INFO `END`
/// for symbolic alleles). Do not assume `end = start` in `ingest_sql`.
///
/// # `alt` is pipe-joined for multi-allelic records — SPLIT IT
///
/// The reader joins a record's ALT alleles into ONE `Utf8` column separated by `|`
/// (`physical_exec.rs`: `join_into(…, '|')`). A multi-allelic record
/// `chr1 100 . A G,T` therefore arrives as a single row with `alt = "G|T"`.
///
/// This is a trap, because the obvious mapping
///
/// ```sql
/// concat(ref, '/', alt) AS allele_string   -- WRONG for multi-allelic records
/// ```
///
/// yields `A/G|T`, which can never equal the runtime probe key (`{ref}/{alt}`, one
/// ALT at a time) nor the variation cache's per-allele `allele_string`. Such rows are
/// written to the shard and are dead forever: they match nothing, and nothing warns.
///
/// `ingest_sql` must therefore split `alt` into ONE ROW PER ALT ALLELE, e.g.
///
/// ```sql
/// SELECT chrom, `start`, `end`,
///        concat(`ref`, '/', alt_one) AS allele_string, …
/// FROM (
///   SELECT *, unnest(string_to_array(alt, '|')) AS alt_one FROM plugin_x_src
/// )
/// ```
///
/// A source whose every record is bi-allelic can skip this, but nothing enforces
/// that — so prefer the split unless the input is known bi-allelic by construction.
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

    /// How to name this source in a diagnostic: its `part` if it has one, else its
    /// 1-based position (so an error always points at exactly one `[[source]]`).
    fn label(&self, index: usize) -> String {
        match &self.part {
            Some(p) => format!("[[source]] part = \"{p}\""),
            None => format!("[[source]] #{}", index + 1),
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
    /// Parse a source manifest from a TOML file, then [`validate`](Self::validate) it.
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

    /// Reject manifests whose declarations contradict each other. `deny_unknown_fields`
    /// catches keys that do not exist; this catches keys that exist but do not belong
    /// together — the ones the build would otherwise accept and silently ignore.
    ///
    /// Rules:
    /// - a `vcf` source requires `coordinate_system = "1-based"`;
    /// - `[source.csv]` is only valid for `provider = "csv" | "tsv"`;
    /// - `[source.vcf]` is only valid for `provider = "vcf"`;
    /// - hence `parquet` (and `bed`) accept neither params block.
    ///
    /// Called from [`load`](Self::load) — the file path every real manifest takes —
    /// and deliberately NOT from `Deserialize`: unit tests across this crate build
    /// manifests straight from TOML strings via `toml::from_str`, and hooking
    /// `Deserialize` would force every one of them through validation.
    pub fn validate(&self) -> Result<()> {
        let plugin = &self.plugin_name;
        for (i, spec) in self.sources.iter().enumerate() {
            let src = spec.label(i);
            let provider = spec.provider.as_manifest_str();

            // provider.rs registers the VCF reader with coordinate_system_zero_based =
            // false, because VCF POS is 1-based. If the manifest nonetheless claims
            // 0-based, the builder believes the manifest and shifts every start by +1
            // (normalize.rs). The build then succeeds, rows > 0, the cache manifest looks
            // healthy — and every runtime key misses, so the plugin annotates nothing.
            if spec.provider == ProviderKind::Vcf
                && self.coordinate_system != CoordinateSystem::OneBased
            {
                return Err(DataFusionError::Execution(format!(
                    "plugin '{plugin}': {src} declares provider = \"vcf\", which requires \
                     coordinate_system = \"1-based\" (VCF POS is 1-based), but the manifest \
                     declares coordinate_system = \"{}\". The build would shift every start \
                     by +1, and the plugin would silently annotate nothing.",
                    self.coordinate_system.as_manifest_str()
                )));
            }

            // A params block that does not match its provider is read by nobody: the
            // provider factory reaches for `spec.csv` only for csv/tsv and `spec.vcf`
            // only for vcf. Declaring the other one is always a mistake, and silently
            // dropping it costs (e.g.) the whole `info_fields` INFO selection.
            if spec.csv.is_some() && !matches!(spec.provider, ProviderKind::Csv | ProviderKind::Tsv)
            {
                return Err(DataFusionError::Execution(format!(
                    "plugin '{plugin}': {src} declares provider = \"{provider}\" with a \
                     [source.csv] block. [source.csv] is only valid for provider = \"csv\" \
                     or \"tsv\"; here it would be silently ignored."
                )));
            }
            if spec.vcf.is_some() && spec.provider != ProviderKind::Vcf {
                return Err(DataFusionError::Execution(format!(
                    "plugin '{plugin}': {src} declares provider = \"{provider}\" with a \
                     [source.vcf] block. [source.vcf] (e.g. info_fields) is only valid for \
                     provider = \"vcf\"; here it would be silently ignored."
                )));
            }
        }
        Ok(())
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

    /// `validate()` runs on the file path, so these go through a real temp file.
    fn load_str(body: &str) -> Result<SourceManifest> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.source.toml");
        std::fs::write(&path, body).unwrap();
        SourceManifest::load(&path)
    }

    /// A manifest with one `[[source]]` of the given provider, coordinate system and
    /// optional params block.
    fn manifest_with(provider: &str, coord: &str, params: &str) -> String {
        format!(
            r##"
plugin_name = "demo"
coordinate_system = "{coord}"
ingest_sql = "SELECT 1"

[[source]]
provider = "{provider}"
path = "/tmp/x"
{params}

[[value_columns]]
column = "s"
csq_field = "S"
type = "Float32"
"##
        )
    }

    // provider.rs hard-wires the VCF reader to `coordinate_system_zero_based = false`
    // because VCF POS is 1-based. Nothing stopped a manifest from ALSO declaring
    // `0-based-half-open`, in which case the builder trusts the manifest and shifts
    // every start by +1 (normalize.rs). The build then succeeds with rows > 0 and a
    // healthy-looking manifest, while every key misses and the plugin annotates
    // nothing. Total, silent failure — so it must be rejected at load.
    #[test]
    fn vcf_source_with_zero_based_coordinates_is_rejected() {
        let err = load_str(&manifest_with("vcf", "0-based-half-open", ""))
            .expect_err("vcf + 0-based-half-open must not load");
        let msg = err.to_string();
        assert!(msg.contains("demo"), "must name the plugin: {msg}");
        assert!(msg.contains("vcf"), "must name the provider: {msg}");
        assert!(
            msg.contains("0-based-half-open"),
            "must name the offending value: {msg}"
        );
        assert!(msg.contains("1-based"), "must say what is required: {msg}");
    }

    #[test]
    fn vcf_source_with_one_based_coordinates_loads() {
        let m = load_str(&manifest_with("vcf", "1-based", "")).expect("vcf + 1-based is valid");
        assert_eq!(m.sources[0].provider, ProviderKind::Vcf);
    }

    // A [source.csv] block under provider = "vcf" parses fine and is then silently
    // discarded by the provider factory. Say so instead.
    #[test]
    fn csv_params_on_a_vcf_source_are_rejected() {
        let err = load_str(&manifest_with(
            "vcf",
            "1-based",
            "  [source.csv]\n  delimiter = \"\\t\"",
        ))
        .expect_err("[source.csv] under provider = \"vcf\" must not load");
        let msg = err.to_string();
        assert!(msg.contains("[source.csv]"), "must name the block: {msg}");
        assert!(msg.contains("vcf"), "must name the provider: {msg}");
    }

    // The mirror: `[source.vcf] info_fields = [...]` under provider = "csv" silently
    // drops the INFO selection.
    #[test]
    fn vcf_params_on_a_csv_source_are_rejected() {
        let err = load_str(&manifest_with(
            "csv",
            "1-based",
            "  [source.vcf]\n  info_fields = [\"AF\"]",
        ))
        .expect_err("[source.vcf] under provider = \"csv\" must not load");
        let msg = err.to_string();
        assert!(msg.contains("[source.vcf]"), "must name the block: {msg}");
        assert!(msg.contains("csv"), "must name the provider: {msg}");
    }

    #[test]
    fn parquet_source_accepts_neither_params_block() {
        let csv_err = load_str(&manifest_with(
            "parquet",
            "1-based",
            "  [source.csv]\n  delimiter = \"\\t\"",
        ))
        .expect_err("parquet + [source.csv] must not load");
        assert!(csv_err.to_string().contains("[source.csv]"), "{csv_err}");

        let vcf_err = load_str(&manifest_with(
            "parquet",
            "1-based",
            "  [source.vcf]\n  info_fields = [\"AF\"]",
        ))
        .expect_err("parquet + [source.vcf] must not load");
        assert!(vcf_err.to_string().contains("[source.vcf]"), "{vcf_err}");

        load_str(&manifest_with("parquet", "1-based", "")).expect("bare parquet source is valid");
    }

    // With several sources, a diagnostic that does not say WHICH one is nearly
    // useless — name the offending `part`.
    #[test]
    fn validation_error_names_the_offending_part() {
        let src = r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
part = "snv"
provider = "csv"
path = "/tmp/snv.tsv"
  [source.csv]
  delimiter = "\t"

[[source]]
part = "sites"
provider = "csv"
path = "/tmp/sites.vcf"
  [source.vcf]
  info_fields = ["AF"]

[[value_columns]]
column = "s"
csq_field = "S"
type = "Float32"
"##;
        let err = load_str(src).expect_err("the second source is misconfigured");
        let msg = err.to_string();
        assert!(
            msg.contains("sites"),
            "must name the offending part, not just the plugin: {msg}"
        );
    }

    // The shape of the real AlphaMissense manifest: tsv + [source.csv] + 1-based.
    // It must keep loading.
    #[test]
    fn tsv_source_with_csv_params_stays_valid() {
        let m = load_str(&manifest_with(
            "tsv",
            "1-based",
            "  [source.csv]\n  delimiter = \"\\t\"\n  compression = \"gzip\"",
        ))
        .expect("tsv + [source.csv] + 1-based is the AlphaMissense shape and must stay valid");
        assert_eq!(m.sources[0].provider, ProviderKind::Tsv);
        assert!(m.sources[0].csv.is_some());
    }

    // A csv/tsv source may legitimately be 0-based; only `vcf` is pinned to 1-based.
    #[test]
    fn zero_based_is_still_allowed_for_non_vcf_providers() {
        let m = load_str(&manifest_with("csv", "0-based-half-open", ""))
            .expect("a 0-based csv source is legitimate (the builder shifts it)");
        assert_eq!(m.coordinate_system, CoordinateSystem::ZeroBasedHalfOpen);
    }

    // validate() hangs off load(), NOT Deserialize: the tests above and elsewhere
    // build manifests straight from TOML strings, and must keep working.
    #[test]
    fn toml_from_str_does_not_validate() {
        let m: SourceManifest =
            toml::from_str(&manifest_with("vcf", "0-based-half-open", "")).unwrap();
        assert!(
            m.validate().is_err(),
            "the same manifest must fail an explicit validate()"
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
        let err =
            toml::from_str::<SourceManifest>(src).expect_err("unknown key [tier] must be rejected");
        assert!(
            err.to_string().contains("tier"),
            "the error must name the offending key, got: {err}"
        );
    }
}
