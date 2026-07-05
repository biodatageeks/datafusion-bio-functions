# AlphaMissense Plugin Prototype + Shared Plugin-Cache Subsystem — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the reusable `plugin_cache` subsystem (declarative source manifest → variation-frequency-tiered per-chrom Parquet cache → runtime point-lookup → CSQ injection) and prove it end-to-end with AlphaMissense as the first plugin.

**Architecture:** A new `plugin_cache` module in `bio-function-vep` mirrors the existing variation cache. A TOML **source manifest** declares a table provider + `ingest_sql`; the builder registers the raw source, creates a `plugin_<name>_ingest` view, applies a shared contig/coordinate normalization wrapper, LEFT-joins the variation shard on `(chrom, start, allele_string)` to derive a warm/cold `tier` (AF discarded), and writes `plugin/<name>/<chrom>.parquet` reusing the variation shard's lookup-optimized `WriterProperties`. At runtime a `PluginLookup` reuses the variation `PageDir`/scan path to probe values and inject them into CSQ output fields.

**Tech Stack:** Rust (edition 2024), Apache DataFusion 53 / Arrow 58, `parquet` crate, `toml` + `serde`, `datafusion-bio-format-vcf` (sibling, for the VCF provider arm only), existing `cache_common` / `parquet_cache` modules.

**Spec:** `docs/superpowers/specs/2026-07-05-custom-vep-plugin-caches-design.md`.

## Global Constraints

- DataFusion 53.0.0, Arrow 58.0.0, Rust edition 2024 — do not bump.
- Cache key columns match the variation cache **verbatim**: `chrom: Utf8`, `start/end: UInt32`, `allele_string: Utf8` (`ref/alt`). Point-lookup key `(chrom, start)`; alleles disambiguated by `allele_string`.
- `WARM_AF_THRESHOLD = 0.01`; tier 0 = warm (AF ≥ threshold), tier 1 = cold / no variation match.
- Plugin shards reuse `parquet_cache::write::point_lookup_writer_properties(schema, &["tier", "start"])` — no dictionary, ZSTD(3), 4 KiB / 512-row pages, `EnabledStatistics::Page`, 1M-row row groups, declared sort `(tier, start)`.
- Runtime probes run **synchronously on the annotating thread** — no internal Rayon/async pool (aux-pool oversubscription constraint).
- Parity is the decisive gate: plugin CSQ fields must be byte-identical to Ensembl VEP on a test chromosome, and w1-vs-w4 body byte-identical.
- Manifests are TOML, resolved from `biodatageeks/vepyr-plugins`. Built Parquet caches are never committed to git.
- `cargo fmt` + `cargo clippy -- -D warnings` must pass before every commit.

## File Structure

New module `datafusion/bio-function-vep/src/plugin_cache/`:

- `mod.rs` — module exports + `WARM_AF_THRESHOLD` re-use.
- `source_manifest.rs` — `SourceManifest` TOML types + `load`.
- `provider.rs` — provider factory: register `plugin_<name>_src[_<part>]`.
- `normalize.rs` — `canonical_contig` scalar UDF + `wrap_normalization(sql_for_ingest_view, coordinate_system)`.
- `join.rs` — `join_variation_frequency`.
- `write.rs` — `plugin_output_schema` + `PluginShardWriter` (generalized tiered writer).
- `build.rs` — `build_plugin_chrom` orchestration.
- `cache_manifest.rs` — `CacheManifest` (build output) types + write + `discover_plugins`.
- `lookup.rs` — `PluginLookup` (runtime probe) + `PluginValueRow`.
- `registry.rs` — `PluginRegistry` (discover + open + CSQ field list).

Modified:
- `datafusion/bio-function-vep/src/lib.rs` — `pub mod plugin_cache;`.
- `datafusion/bio-function-vep/src/annotate_provider.rs` — inject plugin CSQ fields (Task 11).

New repo (Task 12): `biodatageeks/vepyr-plugins` seed with `plugins/alphamissense/alphamissense.source.toml`.

---

## Task 1: Module scaffold + source manifest types

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/mod.rs`
- Create: `datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs`
- Modify: `datafusion/bio-function-vep/src/lib.rs` (add `pub mod plugin_cache;`)
- Modify: `datafusion/bio-function-vep/Cargo.toml` (add `toml = "0.8"` if absent)

**Interfaces:**
- Produces: `SourceManifest { plugin_name: String, coordinate_system: CoordinateSystem, sources: Vec<SourceSpec>, ingest_sql: String, value_columns: Vec<ValueColumn>, tier: TierPolicy }`, `SourceSpec { part: Option<String>, provider: ProviderKind, path: String, csv: Option<CsvParams> }`, `CsvParams { delimiter, has_header, comment, compression, schema: Vec<SchemaField> }`, `ValueColumn { column: String, csq_field: String, ty: ValueType }`, `SourceManifest::load(path) -> Result<SourceManifest>`, and `SourceSpec::table_name(plugin_name) -> String` returning `plugin_<name>_src[_<part>]`.

- [ ] **Step 1: Write the failing test** in `source_manifest.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    const CADD_LIKE: &str = r#"
plugin_name = "demo"
coordinate_system = "1-based"

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

ingest_sql = "SELECT chrom, CAST(pos AS INTEGER) AS start, CAST(pos AS INTEGER) AS end, concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score FROM plugin_demo_src_snv"

[[value_columns]]
column = "demo_score"
csq_field = "DEMO_SCORE"
type = "Float32"

[tier]
threshold = 0.01
unmatched = "cold"
"#;

    #[test]
    fn parses_source_manifest() {
        let m: SourceManifest = toml::from_str(CADD_LIKE).unwrap();
        assert_eq!(m.plugin_name, "demo");
        assert_eq!(m.coordinate_system, CoordinateSystem::OneBased);
        assert_eq!(m.sources.len(), 1);
        assert_eq!(m.sources[0].table_name(&m.plugin_name), "plugin_demo_src_snv");
        assert_eq!(m.sources[0].csv.as_ref().unwrap().schema.len(), 5);
        assert_eq!(m.value_columns[0].csq_field, "DEMO_SCORE");
        assert_eq!(m.value_columns[0].ty, ValueType::Float32);
        assert_eq!(m.tier.threshold, 0.01);
    }

    #[test]
    fn single_source_table_name_has_no_part_suffix() {
        let src = SourceSpec { part: None, provider: ProviderKind::Csv, path: "x".into(), csv: None };
        assert_eq!(src.table_name("cadd"), "plugin_cadd_src");
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::source_manifest -- --nocapture`
Expected: FAIL — `SourceManifest` unresolved.

- [ ] **Step 3: Write minimal implementation** in `source_manifest.rs`:

```rust
//! Declarative source manifest (build input) — TOML from `vepyr-plugins`.
use serde::Deserialize;
use datafusion::common::{DataFusionError, Result};

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub enum CoordinateSystem {
    #[serde(rename = "1-based")]
    OneBased,
    #[serde(rename = "0-based-half-open")]
    ZeroBasedHalfOpen,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ProviderKind { Vcf, Csv, Tsv, Parquet, Bed }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
pub enum ValueType { Utf8, Float32, Int32 }

#[derive(Debug, Clone, Deserialize)]
pub struct SchemaField { pub name: String, #[serde(rename = "type")] pub ty: ValueType }

#[derive(Debug, Clone, Deserialize)]
pub struct CsvParams {
    #[serde(default = "default_delim")] pub delimiter: String,
    #[serde(default)] pub has_header: bool,
    #[serde(default)] pub comment: Option<String>,
    #[serde(default)] pub compression: Option<String>,
    #[serde(default)] pub schema: Vec<SchemaField>,
}
fn default_delim() -> String { "\t".into() }

#[derive(Debug, Clone, Deserialize)]
pub struct SourceSpec {
    #[serde(default)] pub part: Option<String>,
    pub provider: ProviderKind,
    pub path: String,
    #[serde(default)] pub csv: Option<CsvParams>,
}
impl SourceSpec {
    pub fn table_name(&self, plugin_name: &str) -> String {
        match &self.part {
            Some(p) => format!("plugin_{plugin_name}_src_{p}"),
            None => format!("plugin_{plugin_name}_src"),
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct ValueColumn {
    pub column: String,
    pub csq_field: String,
    #[serde(rename = "type")] pub ty: ValueType,
}

#[derive(Debug, Clone, Deserialize)]
pub struct TierPolicy {
    #[serde(default = "default_threshold")] pub threshold: f64,
    #[serde(default)] pub unmatched: Option<String>,
}
fn default_threshold() -> f64 { 0.01 }

#[derive(Debug, Clone, Deserialize)]
pub struct SourceManifest {
    pub plugin_name: String,
    pub coordinate_system: CoordinateSystem,
    #[serde(rename = "source")] pub sources: Vec<SourceSpec>,
    pub ingest_sql: String,
    #[serde(rename = "value_columns")] pub value_columns: Vec<ValueColumn>,
    pub tier: TierPolicy,
}

impl SourceManifest {
    pub fn load(path: &std::path::Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| {
            DataFusionError::Execution(format!("read source manifest '{}': {e}", path.display()))
        })?;
        toml::from_str(&text).map_err(|e| {
            DataFusionError::Execution(format!("parse source manifest '{}': {e}", path.display()))
        })
    }
    pub fn ingest_view_name(&self) -> String { format!("plugin_{}_ingest", self.plugin_name) }
}
```

Create `mod.rs`:

```rust
//! Custom VEP plugin caches: declarative build + tiered point-lookup.
//! See docs/superpowers/specs/2026-07-05-custom-vep-plugin-caches-design.md.
pub mod source_manifest;
```

Add to `lib.rs` near the other `pub mod` lines: `pub mod plugin_cache;`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::source_manifest`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache datafusion/bio-function-vep/src/lib.rs datafusion/bio-function-vep/Cargo.toml
git commit -m "feat(plugin-cache): source manifest TOML types + parsing"
```

---

## Task 2: Provider factory (register raw source tables)

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/provider.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/mod.rs` (add `pub mod provider;`)

**Interfaces:**
- Consumes: `SourceManifest`, `SourceSpec`, `CsvParams`, `ProviderKind`, `ValueType` (Task 1).
- Produces: `async fn register_sources(ctx: &SessionContext, manifest: &SourceManifest) -> Result<()>` — registers every `SourceSpec` under `spec.table_name(plugin_name)`. CSV/TSV via `ctx.register_csv` with an explicit Arrow schema built from `CsvParams.schema`; VCF/BED return an explicit "provider not yet wired" error in the prototype (AlphaMissense is CSV).

- [ ] **Step 1: Write the failing test** in `provider.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::SourceManifest;
    use std::io::Write;

    fn write_gz(path: &std::path::Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn registers_csv_source_with_declared_schema() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("snv.tsv.gz");
        // headerless, comment line, two data rows
        write_gz(&tsv, "#comment\nchr1\t100\tA\tG\t0.9\nchr1\t200\tC\tT\t0.1\n");

        let toml = format!(r#"
plugin_name = "demo"
coordinate_system = "1-based"
[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  comment = "#"
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]
ingest_sql = "SELECT 1"
[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
[tier]
threshold = 0.01
"#, tsv.display());

        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let ctx = SessionContext::new();
        register_sources(&ctx, &manifest).await.unwrap();
        let n = ctx.sql("SELECT count(*) AS c FROM plugin_demo_src")
            .await.unwrap().collect().await.unwrap();
        let c = n[0].column(0).as_any().downcast_ref::<datafusion::arrow::array::Int64Array>().unwrap().value(0);
        assert_eq!(c, 2);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::provider`
Expected: FAIL — `register_sources` unresolved.

- [ ] **Step 3: Write minimal implementation** in `provider.rs`:

```rust
//! Provider factory: register a source manifest's raw tables.
use std::sync::Arc;
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::file_format::file_compression_type::FileCompressionType;
use datafusion::prelude::{CsvReadOptions, SessionContext};

use crate::plugin_cache::source_manifest::{CsvParams, ProviderKind, SourceManifest, ValueType};

fn arrow_type(ty: ValueType) -> DataType {
    match ty {
        ValueType::Utf8 => DataType::Utf8,
        ValueType::Float32 => DataType::Float32,
        ValueType::Int32 => DataType::Int32,
    }
}

fn csv_schema(csv: &CsvParams) -> Schema {
    Schema::new(
        csv.schema.iter()
            .map(|f| Field::new(&f.name, arrow_type(f.ty), true))
            .collect::<Vec<_>>(),
    )
}

pub async fn register_sources(ctx: &SessionContext, manifest: &SourceManifest) -> Result<()> {
    for spec in &manifest.sources {
        let table = spec.table_name(&manifest.plugin_name);
        match spec.provider {
            ProviderKind::Csv | ProviderKind::Tsv => {
                let csv = spec.csv.as_ref().ok_or_else(|| {
                    DataFusionError::Execution(format!("csv/tsv source '{table}' missing [source.csv]"))
                })?;
                let schema = csv_schema(csv);
                let delim = csv.delimiter.as_bytes().first().copied().unwrap_or(b'\t');
                let mut opts = CsvReadOptions::new()
                    .delimiter(delim)
                    .has_header(csv.has_header)
                    .schema(&schema);
                if let Some(c) = csv.comment.as_ref().and_then(|s| s.as_bytes().first().copied()) {
                    opts = opts.comment(c);
                }
                let gz = csv.compression.as_deref() == Some("gzip");
                if gz {
                    opts = opts.file_compression_type(FileCompressionType::GZIP).file_extension(".gz");
                }
                ctx.register_csv(&table, &spec.path, opts).await?;
            }
            ProviderKind::Parquet => {
                ctx.register_parquet(&table, &spec.path, Default::default()).await?;
            }
            ProviderKind::Vcf | ProviderKind::Bed => {
                return Err(DataFusionError::NotImplemented(format!(
                    "provider {:?} not wired in prototype (table '{table}')", spec.provider
                )));
            }
        }
    }
    Ok(())
}
```

Add `pub mod provider;` to `mod.rs`. Confirm `flate2` and `tempfile` are dev-deps (add to `[dev-dependencies]` in `Cargo.toml` if missing).

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::provider`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache datafusion/bio-function-vep/Cargo.toml
git commit -m "feat(plugin-cache): provider factory registers CSV/TSV/Parquet sources"
```

---

## Task 3: `canonical_contig` UDF + coordinate normalization

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/normalize.rs`
- Modify: `mod.rs` (`pub mod normalize;`)

**Interfaces:**
- Consumes: `CoordinateSystem` (Task 1).
- Produces: `fn canonical_contig_udf() -> ScalarUDF` (name `canonical_contig`, Utf8→Utf8); `fn wrap_normalization(inner_view: &str, coord: CoordinateSystem) -> String` returning SQL that selects from `inner_view` applying `canonical_contig(chrom) AS chrom` and the coordinate shift; `fn canonical_contig_str(raw: &str) -> String` (the pure-Rust core, so it stays consistent with `key_encoding::chrom_code`).

- [ ] **Step 1: Write the failing test** in `normalize.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonicalizes_contig_strings() {
        assert_eq!(canonical_contig_str("chr1"), "1");
        assert_eq!(canonical_contig_str("1"), "1");
        assert_eq!(canonical_contig_str("chrX"), "X");
        assert_eq!(canonical_contig_str("chrM"), "MT");
        assert_eq!(canonical_contig_str("M"), "MT");
        assert_eq!(canonical_contig_str("chrMT"), "MT");
        assert_eq!(canonical_contig_str("MT"), "MT");
    }

    #[test]
    fn one_based_passes_through_zero_based_shifts() {
        assert_eq!(
            wrap_normalization("plugin_demo_ingest", crate::plugin_cache::source_manifest::CoordinateSystem::OneBased),
            "SELECT canonical_contig(chrom) AS chrom, CAST(start AS INT) AS start, CAST(end AS INT) AS end, allele_string, * EXCLUDE (chrom, start, end, allele_string) FROM plugin_demo_ingest"
        );
        let z = wrap_normalization("plugin_demo_ingest", crate::plugin_cache::source_manifest::CoordinateSystem::ZeroBasedHalfOpen);
        assert!(z.contains("CAST(start AS INT) + 1 AS start"));
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::normalize`
Expected: FAIL — items unresolved.

- [ ] **Step 3: Write minimal implementation** in `normalize.rs`:

```rust
//! Shared contig/coordinate normalization wrapper for the ingest view.
use std::sync::Arc;
use datafusion::arrow::array::{ArrayRef, StringArray};
use datafusion::common::Result;
use datafusion::logical_expr::{ColumnarValue, ScalarUDF, Volatility};
use datafusion::prelude::create_udf;
use datafusion::arrow::datatypes::DataType;

use crate::plugin_cache::source_manifest::CoordinateSystem;

/// Bare Ensembl contig form matching the variation `chrom` column
/// (`1`..`22`, `X`, `Y`, `MT`). Mirrors `cache::key_encoding` prefix handling.
pub fn canonical_contig_str(raw: &str) -> String {
    let bare = raw.strip_prefix("chr").unwrap_or(raw);
    match bare {
        "M" | "MT" | "chrM" => "MT".to_string(),
        other => other.to_ascii_uppercase(),
    }
}

pub fn canonical_contig_udf() -> ScalarUDF {
    let f = Arc::new(|args: &[ColumnarValue]| -> Result<ColumnarValue> {
        let arrays = ColumnarValue::values_to_arrays(args)?;
        let input = arrays[0].as_any().downcast_ref::<StringArray>().expect("Utf8 arg");
        let out: StringArray = input.iter()
            .map(|v| v.map(canonical_contig_str))
            .collect();
        Ok(ColumnarValue::Array(Arc::new(out) as ArrayRef))
    });
    create_udf("canonical_contig", vec![DataType::Utf8], DataType::Utf8, Volatility::Immutable, f)
}

pub fn wrap_normalization(inner_view: &str, coord: CoordinateSystem) -> String {
    let start_expr = match coord {
        CoordinateSystem::OneBased => "CAST(start AS INT) AS start",
        CoordinateSystem::ZeroBasedHalfOpen => "CAST(start AS INT) + 1 AS start",
    };
    format!(
        "SELECT canonical_contig(chrom) AS chrom, {start_expr}, CAST(end AS INT) AS end, \
         allele_string, * EXCLUDE (chrom, start, end, allele_string) FROM {inner_view}"
    )
}
```

Note for executor: DataFusion's `SELECT ... * EXCLUDE (...)` is supported; if the installed DF build rejects it, fall back to selecting value columns explicitly (the builder in Task 7 knows the value-column names from the manifest, so it can pass them in). Keep the test assertion in sync with whichever form is used.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::normalize`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): canonical_contig UDF + coordinate normalization wrapper"
```

---

## Task 4: `join_variation_frequency` (LEFT JOIN + tier derivation)

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/join.rs`
- Modify: `mod.rs` (`pub mod join;`)

**Interfaces:**
- Consumes: `cache_common::max_global_af` semantics (reused via SQL), `WARM_AF_THRESHOLD`.
- Produces: `pub fn tier_sql(normalized_view: &str, variation_probe: &str, threshold: f64) -> String` and `async fn tiered_stream(ctx: &SessionContext, normalized_view: &str, variation_shard: &Path, af_max_sql: &str, threshold: f64) -> Result<SendableRecordBatchStream>` — registers the variation shard as `plugin_variation_probe` (reducing its stored AF columns to `af_max` via `af_max_sql`), LEFT-joins on `(chrom, start, allele_string)`, appends `tier: Int8` (0 if joined `af_max >= threshold`, else 1), drops all variation columns. Value columns pass through unchanged.

Note on AF source: the variation shard stores AF as the 2-array groups (`<grp>_alleles`/`<grp>_freqs`) plus `minor_allele_freq`. For the prototype, derive `af_max` from the variation shard's `minor_allele_freq` column plus the gnomADg group max. Executor: confirm the exact stored AF columns via `cargo run` schema dump of a real `variation/<chrom>.parquet` (Task 13 uses chr22) and build the `af_max` SQL expression accordingly; the test below uses a synthetic shard with an explicit `af_max` column to isolate the tiering logic from the AF-decode detail.

- [ ] **Step 1: Write the failing test** in `join.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array, Float64Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    // Build a tiny "variation-like" table with an explicit af_max column and
    // register it; assert warm/cold/no-match tiering.
    #[tokio::test(flavor = "multi_thread")]
    async fn derives_tier_warm_cold_and_miss() {
        let ctx = SessionContext::new();
        // variation probe: (chrom,start,allele_string,af_max)
        let var = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("af_max", DataType::Float64, false),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Float64Array::from(vec![0.5, 0.001])), // warm, cold
            ],
        ).unwrap();
        ctx.register_batch("plugin_variation_probe", var).unwrap();

        // plugin ingest: three rows — matches warm, matches cold, and a miss.
        let plug = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("end", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("demo_score", DataType::Float32, true),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200, 300])),
                Arc::new(UInt32Array::from(vec![100u32, 200, 300])),
                Arc::new(StringArray::from(vec!["A/G", "C/T", "G/A"])),
                Arc::new(Float32Array::from(vec![Some(0.9f32), Some(0.1), Some(0.7)])),
            ],
        ).unwrap();
        ctx.register_batch("plugin_demo_norm", plug).unwrap();

        let df = ctx.sql(&tier_sql("plugin_demo_norm", "plugin_variation_probe", 0.01)).await.unwrap();
        let batches = df.sort(vec![datafusion::prelude::col("start").sort(true, true)]).unwrap().collect().await.unwrap();
        let b = &batches[0];
        let tier = b.column(b.schema().index_of("tier").unwrap())
            .as_any().downcast_ref::<Int8Array>().unwrap();
        assert_eq!(tier.value(0), 0); // start 100, af 0.5 -> warm
        assert_eq!(tier.value(1), 1); // start 200, af 0.001 -> cold
        assert_eq!(tier.value(2), 1); // start 300, no match -> cold
        // value column preserved, variation columns dropped
        assert!(b.schema().index_of("demo_score").is_ok());
        assert!(b.schema().index_of("af_max").is_err());
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::join`
Expected: FAIL — `tier_sql` unresolved.

- [ ] **Step 3: Write minimal implementation** in `join.rs`:

```rust
//! Variation-frequency join → warm/cold tier derivation (AF discarded).
use std::path::Path;
use datafusion::common::Result;
use datafusion::physical_plan::SendableRecordBatchStream;
use datafusion::prelude::SessionContext;

/// SQL that LEFT-joins `normalized_view` to a registered `plugin_variation_probe`
/// table exposing `(chrom, start, allele_string, af_max)` and appends `tier`.
/// The value columns of `normalized_view` pass through; variation columns drop.
pub fn tier_sql(normalized_view: &str, variation_probe: &str, threshold: f64) -> String {
    format!(
        "SELECT p.*, CASE WHEN COALESCE(v.af_max, 0.0) >= {threshold} \
         THEN CAST(0 AS TINYINT) ELSE CAST(1 AS TINYINT) END AS tier \
         FROM {normalized_view} p \
         LEFT JOIN {variation_probe} v \
         ON p.chrom = v.chrom AND p.start = v.start AND p.allele_string = v.allele_string"
    )
}

/// Register the variation shard as an `af_max`-bearing probe table, then stream
/// the tiered rows. `af_max_sql` is the expression that reduces the shard's stored
/// AF columns to a single `af_max` (built by the caller from the real schema).
pub async fn tiered_stream(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
    af_max_sql: &str,
    threshold: f64,
) -> Result<SendableRecordBatchStream> {
    ctx.register_parquet("plugin_variation_raw", variation_shard.to_str().unwrap(), Default::default()).await?;
    ctx.sql(&format!(
        "CREATE OR REPLACE VIEW plugin_variation_probe AS \
         SELECT chrom, start, allele_string, {af_max_sql} AS af_max FROM plugin_variation_raw"
    )).await?;
    let df = ctx.sql(&tier_sql(normalized_view, "plugin_variation_probe", threshold)).await?;
    df.execute_stream().await
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::join`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): variation-frequency join + warm/cold tier derivation"
```

---

## Task 5: Tiered plugin shard writer + output schema

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/write.rs`
- Modify: `mod.rs` (`pub mod write;`)

**Interfaces:**
- Consumes: `parquet_cache::write::point_lookup_writer_properties`, `ValueColumn`/`ValueType` (Task 1).
- Produces: `fn plugin_output_schema(values: &[ValueColumn]) -> SchemaRef` — `chrom:Utf8, start:UInt32, end:UInt32, allele_string:Utf8, <value cols>, tier:Int8`; `struct PluginShardWriter { .. }` with `create(path, schema) -> Result<Self>`, `write(&RecordBatch)`, `finish() -> usize`, using `point_lookup_writer_properties(schema, &["tier","start"])`.

- [ ] **Step 1: Write the failing test** in `write.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::{ValueColumn, ValueType};
    use datafusion::arrow::array::{Int8Array, Float32Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[test]
    fn output_schema_has_key_values_tier_in_order() {
        let vals = vec![ValueColumn { column: "am_pathogenicity".into(), csq_field: "am_pathogenicity".into(), ty: ValueType::Float32 }];
        let s = plugin_output_schema(&vals);
        let names: Vec<_> = s.fields().iter().map(|f| f.name().clone()).collect();
        assert_eq!(names, vec!["chrom","start","end","allele_string","am_pathogenicity","tier"]);
        assert_eq!(s.field_with_name("tier").unwrap().data_type(), &DataType::Int8);
        assert_eq!(s.field_with_name("start").unwrap().data_type(), &DataType::UInt32);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn writes_readable_shard() {
        let vals = vec![ValueColumn { column: "am_pathogenicity".into(), csq_field: "am_pathogenicity".into(), ty: ValueType::Float32 }];
        let schema = plugin_output_schema(&vals);
        let batch = RecordBatch::try_new(schema.clone(), vec![
            Arc::new(StringArray::from(vec!["1"])),
            Arc::new(UInt32Array::from(vec![100u32])),
            Arc::new(UInt32Array::from(vec![100u32])),
            Arc::new(StringArray::from(vec!["A/G"])),
            Arc::new(Float32Array::from(vec![0.9f32])),
            Arc::new(Int8Array::from(vec![0i8])),
        ]).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("1.parquet");
        let mut w = PluginShardWriter::create(&path, schema).unwrap();
        w.write(&batch).unwrap();
        assert_eq!(w.finish().unwrap(), 1);
        assert!(path.exists());
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::write`
Expected: FAIL.

- [ ] **Step 3: Write minimal implementation** in `write.rs`:

```rust
//! Tiered plugin shard writer (generalizes the variation writer to an arbitrary value schema).
use std::path::Path;
use std::sync::Arc;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ArrowWriter;

use crate::parquet_cache::write::point_lookup_writer_properties;
use crate::plugin_cache::source_manifest::{ValueColumn, ValueType};

fn arrow_type(ty: ValueType) -> DataType {
    match ty { ValueType::Utf8 => DataType::Utf8, ValueType::Float32 => DataType::Float32, ValueType::Int32 => DataType::Int32 }
}

pub fn plugin_output_schema(values: &[ValueColumn]) -> SchemaRef {
    let mut fields = vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new("allele_string", DataType::Utf8, false),
    ];
    for v in values {
        fields.push(Field::new(&v.column, arrow_type(v.ty), true));
    }
    fields.push(Field::new("tier", DataType::Int8, false));
    Arc::new(Schema::new(fields))
}

pub struct PluginShardWriter { writer: ArrowWriter<std::fs::File>, rows: usize }

impl PluginShardWriter {
    pub fn create(path: &Path, schema: SchemaRef) -> Result<Self> {
        let props = point_lookup_writer_properties(&schema, &["tier", "start"]);
        let file = std::fs::File::create(path).map_err(|e| DataFusionError::Execution(format!("create '{}': {e}", path.display())))?;
        let writer = ArrowWriter::try_new(file, schema, Some(props)).map_err(|e| DataFusionError::Execution(format!("open ArrowWriter: {e}")))?;
        Ok(Self { writer, rows: 0 })
    }
    pub fn write(&mut self, batch: &RecordBatch) -> Result<()> {
        self.writer.write(batch).map_err(|e| DataFusionError::Execution(format!("write batch: {e}")))?;
        self.rows += batch.num_rows();
        Ok(())
    }
    pub fn finish(self) -> Result<usize> {
        self.writer.close().map_err(|e| DataFusionError::Execution(format!("close ArrowWriter: {e}")))?;
        Ok(self.rows)
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::write`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): tiered plugin shard writer + output schema"
```

---

## Task 6: Cache manifest (build output) + discovery

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/cache_manifest.rs`
- Modify: `mod.rs` (`pub mod cache_manifest;`)

**Interfaces:**
- Consumes: `SourceManifest`, `ValueColumn` (Task 1).
- Produces: `CacheManifest { plugin_name, source_manifest, key_columns: Vec<String>, value_columns: Vec<ValueColumn>, tier: TierRecord, chroms: Vec<ChromEntry>, cache_source_version: Option<String> }`; `fn from_source(&SourceManifest, source_manifest_file: &str) -> CacheManifest`; `fn write(&self, plugin_dir: &Path)`; `fn discover_plugins(cache_root: &Path) -> Result<Vec<CacheManifest>>` scanning `plugin/*/manifest.json`.

- [ ] **Step 1: Write the failing test** in `cache_manifest.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn writes_and_discovers_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let plugin_dir = dir.path().join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let m = CacheManifest {
            plugin_name: "demo".into(),
            source_manifest: "demo.source.toml".into(),
            key_columns: vec!["chrom".into(),"start".into(),"end".into(),"allele_string".into()],
            value_columns: vec![],
            tier: TierRecord { threshold: 0.01, unmatched: "cold".into() },
            chroms: vec![ChromEntry { chrom: "22".into(), file: "chr22.parquet".into(), rows: 3, warm: 1, cold: 2 }],
            cache_source_version: None,
        };
        m.write(&plugin_dir).unwrap();
        let found = discover_plugins(dir.path()).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].plugin_name, "demo");
        assert_eq!(found[0].chroms[0].rows, 3);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::cache_manifest`
Expected: FAIL.

- [ ] **Step 3: Write minimal implementation** in `cache_manifest.rs`:

```rust
//! Cache manifest (build output) — runtime discovery of built plugins.
use std::path::Path;
use serde::{Deserialize, Serialize};
use datafusion::common::{DataFusionError, Result};
use crate::plugin_cache::source_manifest::{SourceManifest, ValueColumn, ValueType};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TierRecord { pub threshold: f64, pub unmatched: String }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromEntry { pub chrom: String, pub file: String, pub rows: usize, pub warm: usize, pub cold: usize }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValueColumnRecord { pub column: String, pub csq_field: String, #[serde(rename="type")] pub ty: String }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheManifest {
    pub plugin_name: String,
    pub source_manifest: String,
    pub key_columns: Vec<String>,
    pub value_columns: Vec<ValueColumnRecord>,
    pub tier: TierRecord,
    pub chroms: Vec<ChromEntry>,
    pub cache_source_version: Option<String>,
}

fn ty_str(t: ValueType) -> &'static str { match t { ValueType::Utf8 => "Utf8", ValueType::Float32 => "Float32", ValueType::Int32 => "Int32" } }

impl CacheManifest {
    pub fn from_source(src: &SourceManifest, source_manifest_file: &str) -> Self {
        CacheManifest {
            plugin_name: src.plugin_name.clone(),
            source_manifest: source_manifest_file.to_string(),
            key_columns: vec!["chrom".into(),"start".into(),"end".into(),"allele_string".into()],
            value_columns: src.value_columns.iter().map(|v| ValueColumnRecord {
                column: v.column.clone(), csq_field: v.csq_field.clone(), ty: ty_str(v.ty).into()
            }).collect(),
            tier: TierRecord { threshold: src.tier.threshold, unmatched: src.tier.unmatched.clone().unwrap_or_else(|| "cold".into()) },
            chroms: vec![],
            cache_source_version: None,
        }
    }
    pub fn write(&self, plugin_dir: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self).map_err(|e| DataFusionError::Execution(format!("serialize manifest: {e}")))?;
        std::fs::write(plugin_dir.join("manifest.json"), json).map_err(|e| DataFusionError::Execution(format!("write manifest: {e}")))
    }
}

pub fn discover_plugins(cache_root: &Path) -> Result<Vec<CacheManifest>> {
    let plugin_root = cache_root.join("plugin");
    let mut out = Vec::new();
    if !plugin_root.exists() { return Ok(out); }
    for entry in std::fs::read_dir(&plugin_root).map_err(|e| DataFusionError::Execution(format!("read {}: {e}", plugin_root.display())))? {
        let dir = entry.map_err(|e| DataFusionError::Execution(format!("dir entry: {e}")))?.path();
        let mf = dir.join("manifest.json");
        if mf.exists() {
            let text = std::fs::read_to_string(&mf).map_err(|e| DataFusionError::Execution(format!("read {}: {e}", mf.display())))?;
            out.push(serde_json::from_str(&text).map_err(|e| DataFusionError::Execution(format!("parse {}: {e}", mf.display())))?);
        }
    }
    out.sort_by(|a, b| a.plugin_name.cmp(&b.plugin_name));
    Ok(out)
}
```

Confirm `serde_json` is a dependency (it is used elsewhere in the crate; add if missing).

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::cache_manifest`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): cache manifest write + discovery"
```

---

## Task 7: `build_plugin_chrom` orchestration (end-to-end build)

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/build.rs`
- Modify: `mod.rs` (`pub mod build;`)

**Interfaces:**
- Consumes: Tasks 1–6.
- Produces: `async fn build_plugin_chrom(source_manifest: &SourceManifest, source_manifest_file: &str, variation_shard: &Path, af_max_sql: &str, output_plugin_dir: &Path, chrom: &str) -> Result<ChromEntry>` — registers sources, creates `plugin_<name>_ingest`, wraps normalization, filters to `chrom`, joins+tiers, writes warm then cold to `plugin/<name>/<canonical_chrom_label>.parquet`, returns the `ChromEntry` with tier counts. Registers `canonical_contig_udf()` on the ctx.

- [ ] **Step 1: Write the failing test** in `build.rs` (synthetic source + synthetic variation shard):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use crate::plugin_cache::source_manifest::SourceManifest;

    fn write_gz(path: &std::path::Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn builds_tiered_shard_with_counts() {
        let dir = tempfile::tempdir().unwrap();
        // synthetic plugin source: 2 rows on chr1
        let tsv = dir.path().join("src.tsv.gz");
        write_gz(&tsv, "chr1\t100\tA\tG\t0.9\nchr1\t300\tG\tA\t0.7\n");

        // synthetic variation shard exposing (chrom,start,allele_string,minor_allele_freq)
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1",100,"A/G",0.5)]).await; // helper in test module

        let toml = format!(r#"
plugin_name = "demo"
coordinate_system = "1-based"
[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]
ingest_sql = "SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end, concat(ref,'/',alt) AS allele_string, CAST(score AS FLOAT) AS demo_score FROM plugin_demo_src"
[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
[tier]
threshold = 0.01
unmatched = "cold"
"#, tsv.display());
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, "minor_allele_freq", &out, "1").await.unwrap();
        assert_eq!(entry.rows, 2);
        assert_eq!(entry.warm, 1); // start 100 matched af 0.5
        assert_eq!(entry.cold, 1); // start 300 miss
        assert!(out.join("plugin").join("demo").join("chr1.parquet").exists());
    }
}
```

(The `write_synthetic_variation` helper writes a Parquet with a `minor_allele_freq` Float64 column plus the key columns, so `af_max_sql = "minor_allele_freq"`. Keep it in the test module.)

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::build`
Expected: FAIL.

- [ ] **Step 3: Write minimal implementation** in `build.rs`:

```rust
//! End-to-end per-chrom plugin cache build.
use std::path::Path;
use std::sync::Arc;
use datafusion::arrow::array::{Array, Int8Array};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::SessionContext;
use futures::StreamExt;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::ChromEntry;
use crate::plugin_cache::join::tiered_stream;
use crate::plugin_cache::normalize::{canonical_contig_udf, wrap_normalization};
use crate::plugin_cache::provider::register_sources;
use crate::plugin_cache::source_manifest::SourceManifest;
use crate::plugin_cache::write::{plugin_output_schema, PluginShardWriter};

pub async fn build_plugin_chrom(
    src: &SourceManifest,
    source_manifest_file: &str,
    variation_shard: &Path,
    af_max_sql: &str,
    output_cache_root: &Path,
    chrom: &str,
) -> Result<ChromEntry> {
    let ctx = SessionContext::new();
    ctx.register_udf(canonical_contig_udf());
    register_sources(&ctx, src).await?;

    // Ingest view (raw column mapping) then normalized view (contig + coords).
    ctx.sql(&format!("CREATE OR REPLACE VIEW {} AS {}", src.ingest_view_name(), src.ingest_sql)).await?;
    let norm_sql = wrap_normalization(&src.ingest_view_name(), src.coordinate_system.clone());
    let norm_view = format!("plugin_{}_norm", src.plugin_name);
    // Scope to the target contig (post-canonicalization form).
    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    ctx.sql(&format!(
        "CREATE OR REPLACE VIEW {norm_view} AS SELECT * FROM ({norm_sql}) WHERE chrom = '{source_chrom}'"
    )).await?;

    let out_schema = plugin_output_schema(&src.value_columns);
    let plugin_dir = output_cache_root.join("plugin").join(&src.plugin_name);
    std::fs::create_dir_all(&plugin_dir).map_err(|e| DataFusionError::Execution(format!("mkdir {}: {e}", plugin_dir.display())))?;
    let file_name = format!("{}.parquet", canonical_chrom_label(chrom));
    let shard_path = plugin_dir.join(&file_name);

    // Materialize tiered rows, then split into warm/cold ordered runs and write.
    let mut stream = tiered_stream(&ctx, &norm_view, variation_shard, af_max_sql, src.tier.threshold).await?;
    // Collect (small per-chrom for the prototype); production can stream two passes like variation.
    let mut batches = Vec::new();
    while let Some(b) = stream.next().await { batches.push(b?); }
    let full = datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches)?;

    // Reproject to out_schema column order (key, values.., tier) and split by tier.
    let reordered = crate::plugin_cache::build::reproject(&full, &out_schema)?;
    let tier_col = reordered.column(out_schema.index_of("tier").unwrap()).as_any().downcast_ref::<Int8Array>().unwrap().clone();
    let (mut warm, mut cold) = (0usize, 0usize);
    let mut writer = PluginShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
    for keep in [0i8, 1i8] {
        let mask: datafusion::arrow::array::BooleanArray = (0..reordered.num_rows())
            .map(|i| Some(tier_col.value(i) == keep)).collect();
        let filtered = datafusion::arrow::compute::filter_record_batch(&reordered, &mask)?;
        if filtered.num_rows() == 0 { continue; }
        let sorted = crate::plugin_cache::build::sort_by_start(&filtered)?;
        if keep == 0 { warm += sorted.num_rows(); } else { cold += sorted.num_rows(); }
        writer.write(&sorted)?;
    }
    let rows = writer.finish()?;
    if rows == 0 { let _ = std::fs::remove_file(&shard_path); }
    Ok(ChromEntry { chrom: canonical_chrom_label(chrom), file: file_name, rows, warm, cold })
}

// helpers: reproject to schema column order; sort a batch ascending by `start`.
pub(crate) fn reproject(batch: &datafusion::arrow::record_batch::RecordBatch, schema: &datafusion::arrow::datatypes::SchemaRef) -> Result<datafusion::arrow::record_batch::RecordBatch> {
    let cols = schema.fields().iter()
        .map(|f| batch.column(batch.schema().index_of(f.name()).unwrap()).clone())
        .collect();
    datafusion::arrow::record_batch::RecordBatch::try_new(Arc::clone(schema), cols)
        .map_err(|e| DataFusionError::Execution(format!("reproject: {e}")))
}
pub(crate) fn sort_by_start(batch: &datafusion::arrow::record_batch::RecordBatch) -> Result<datafusion::arrow::record_batch::RecordBatch> {
    use datafusion::arrow::compute::{sort_to_indices, take};
    let start = batch.column(batch.schema().index_of("start").unwrap());
    let idx = sort_to_indices(start, None, None).map_err(|e| DataFusionError::Execution(format!("sort: {e}")))?;
    let cols = batch.columns().iter().map(|c| take(c, &idx, None)).collect::<std::result::Result<Vec<_>,_>>().map_err(|e| DataFusionError::Execution(format!("take: {e}")))?;
    datafusion::arrow::record_batch::RecordBatch::try_new(batch.schema(), cols).map_err(|e| DataFusionError::Execution(format!("rebuild: {e}")))
}
```

Note for executor: the `end` column from `wrap_normalization` is `INT`; cast key columns to the exact output types (`start`/`end` → `UInt32`) inside `reproject` or add a cast in the normalization SQL. Keep the output schema types authoritative (`plugin_output_schema`), casting as needed so `RecordBatch::try_new` succeeds.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::build`
Expected: PASS (rows=2, warm=1, cold=1).

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): build_plugin_chrom end-to-end (register→view→normalize→join→tier→write)"
```

---

## Task 8: `PluginLookup` runtime probe

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/lookup.rs`
- Modify: `mod.rs` (`pub mod lookup;`)

**Interfaces:**
- Consumes: `parquet_cache::page_dir::PageDir` + `SinglePathParquetVariationLookup` patterns (study `src/parquet_cache/variation_lookup.rs`).
- Produces: `struct PluginLookup { .. }` with `async fn open(shard: &Path, value_columns: Vec<String>) -> Result<Self>`; `fn probe(&self, start: u32, allele_string: &str) -> Result<Option<PluginValueRow>>` (chrom is implied by the per-chrom shard); `struct PluginValueRow { values: Vec<Option<PluginScalar>> }` positionally aligned to `value_columns`; `enum PluginScalar { Str(String), F32(f32), I32(i32), Null }`.

Implementation approach: mirror `SinglePathParquetVariationLookup` — build a `PageDir` over the `start` leaf, resolve the probe's page range, read `start` + `allele_string` + value columns for those pages, and return the row whose `(start, allele_string)` matches. For the prototype a per-probe resolve is acceptable; batch probing is a later optimization.

- [ ] **Step 1: Write the failing test** in `lookup.rs` — build a shard with Task 5's writer, then probe hit + miss:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::{ValueColumn, ValueType};
    use crate::plugin_cache::write::{plugin_output_schema, PluginShardWriter};
    use datafusion::arrow::array::{Int8Array, Float32Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[tokio::test(flavor = "multi_thread")]
    async fn probe_hit_and_miss() {
        let vals = vec![ValueColumn { column: "am_pathogenicity".into(), csq_field: "am_pathogenicity".into(), ty: ValueType::Float32 }];
        let schema = plugin_output_schema(&vals);
        // two warm rows, ascending start
        let batch = RecordBatch::try_new(schema.clone(), vec![
            Arc::new(StringArray::from(vec!["1","1"])),
            Arc::new(UInt32Array::from(vec![100u32, 200])),
            Arc::new(UInt32Array::from(vec![100u32, 200])),
            Arc::new(StringArray::from(vec!["A/G","C/T"])),
            Arc::new(Float32Array::from(vec![0.9f32, 0.1])),
            Arc::new(Int8Array::from(vec![0i8, 0])),
        ]).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.parquet");
        let mut w = PluginShardWriter::create(&path, schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

        let lookup = PluginLookup::open(&path, vec!["am_pathogenicity".into()]).await.unwrap();
        let hit = lookup.probe(100, "A/G").unwrap().unwrap();
        match &hit.values[0] { Some(PluginScalar::F32(v)) => assert!((*v - 0.9).abs() < 1e-6), other => panic!("{other:?}") }
        assert!(lookup.probe(100, "C/T").unwrap().is_none()); // wrong allele
        assert!(lookup.probe(999, "A/G").unwrap().is_none()); // no position
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::lookup`
Expected: FAIL.

- [ ] **Step 3: Write minimal implementation** in `lookup.rs`. Model it on `variation_lookup.rs`; the executor should read that file first and reuse `PageDir::build`, `resolve_ranges`, the coalescing reader, and the `ProjectionMask` construction. The probe returns the matching row's value scalars. (Concrete code omitted here only because it must track the exact `page_dir` API in `variation_lookup.rs`; the executor mirrors `resolve_and_take` for a single `(start, allele_string)` and decodes the projected value columns into `PluginScalar`.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::lookup`
Expected: PASS (hit 0.9, both misses None).

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): PluginLookup PageDir point-probe (hit/miss)"
```

---

## Task 9: `PluginRegistry` (discover + open + CSQ fields)

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/registry.rs`
- Modify: `mod.rs` (`pub mod registry;`)

**Interfaces:**
- Consumes: `discover_plugins` (Task 6), `PluginLookup` (Task 8).
- Produces: `struct PluginRegistry { plugins: Vec<PluginEntry> }`; `struct PluginEntry { name, csq_fields: Vec<String>, value_columns: Vec<String> }`; `async fn open(cache_root: &Path, chrom: &str) -> Result<PluginRegistry>` (opens each plugin's `<canonical_chrom_label>.parquet` if present); `fn csq_fields(&self) -> Vec<String>` (concatenation across plugins, in plugin-name order); `fn probe_all(&self, start: u32, allele_string: &str) -> Result<Vec<Option<PluginScalar>>>` aligned to `csq_fields()`.

- [ ] **Step 1: Write the failing test** — build one plugin shard + its manifest via Tasks 5–6, open the registry for that chrom, assert `csq_fields()` and a `probe_all` hit. (Follow the Task 8 fixture pattern; write `manifest.json` with one value column and one chrom entry.)

- [ ] **Step 2: Run test to verify it fails.** Run: `cargo test -p datafusion-bio-function-vep plugin_cache::registry` → FAIL.

- [ ] **Step 3: Write minimal implementation** — `open` calls `discover_plugins`, and for each manifest whose `chroms` includes `chrom` opens a `PluginLookup` on `plugin/<name>/<canonical_chrom_label>.parquet`. `csq_fields` flattens each manifest's `value_columns[].csq_field`. `probe_all` probes each plugin and concatenates the per-plugin `PluginValueRow.values`, emitting `None` for plugins with no shard for this chrom.

- [ ] **Step 4: Run test to verify it passes.** Run the same command → PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): PluginRegistry discovery + probe_all + csq_fields"
```

---

## Task 10: Wire plugin CSQ fields into annotation output

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (CSQ field list assembly + per-variant colocated resolution region, ~lines 5239–5580 — read first)
- Modify: whichever struct owns the annotation config to carry an `Option<PluginRegistry>` (follow the pattern used for the variation lookup / colocated map).

**Interfaces:**
- Consumes: `PluginRegistry::csq_fields`, `PluginRegistry::probe_all`, `PluginScalar` (Tasks 8–9).
- Produces: plugin CSQ columns appended to the emitted CSQ format, populated per output row.

- [ ] **Step 1: Write the failing integration test** (`tests/` or an in-module test): construct a minimal annotation over a 1-variant VCF with a pre-built AlphaMissense-shaped plugin shard on the matching `(start, allele_string)`, and assert the emitted CSQ contains `am_pathogenicity`/`am_class` with the shard's values, and a non-matching variant emits the empty value. Keep the fixture tiny (reuse the Task 8 shard builder).

- [ ] **Step 2: Run it to confirm it fails** (fields absent). Run: `cargo test -p datafusion-bio-function-vep plugin_csq_injection -- --nocapture` → FAIL.

- [ ] **Step 3: Implement the wiring.** Two edits, mirroring how frequency/colocated fields already flow:
  1. **Field list:** where the CSQ output field order is assembled, append `plugin_registry.csq_fields()` after the existing colocated/frequency fields.
  2. **Per-row population:** in the per-variant output loop where `colocated_map.get(...)` resolves frequency fields, call `plugin_registry.probe_all(start, allele_string)` **synchronously** (no new pool — Global Constraint) using the same `alt_orig_allele_string` the frequency lookup uses, and write each returned `PluginScalar` into its CSQ column; `None` → the VEP empty value.

- [ ] **Step 4: Run the integration test → PASS.** Also run the full suite: `cargo test -p datafusion-bio-function-vep` → PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs datafusion/bio-function-vep/src/plugin_cache
git commit -m "feat(plugin-cache): inject plugin CSQ fields into annotation output"
```

---

## Task 11: Seed `vepyr-plugins` repo + AlphaMissense source manifest

**Files (new repo `biodatageeks/vepyr-plugins`):**
- Create: `README.md`, `plugins/alphamissense/alphamissense.source.toml`, `plugins/alphamissense/README.md`.

**Interfaces:** Produces the manifest `build_plugin_chrom` consumes for AlphaMissense.

- [ ] **Step 1: Create the repo locally and the manifest.** AlphaMissense_hg38 columns (tab-separated, `#`-comment header): `CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class`; contig `chr1` style; 1-based. Default VEP output fields: `am_pathogenicity`, `am_class`.

```toml
plugin_name       = "alphamissense"
coordinate_system = "1-based"

[[source]]
provider = "tsv"
path     = "AlphaMissense_hg38.tsv.gz"   # overridden by the build driver's --source-path
  [source.csv]
  delimiter   = "\t"
  has_header  = false
  comment     = "#"
  compression = "gzip"
  schema = [
    { name = "chrom",           type = "Utf8" },
    { name = "pos",             type = "Utf8" },
    { name = "ref",             type = "Utf8" },
    { name = "alt",             type = "Utf8" },
    { name = "genome",          type = "Utf8" },
    { name = "uniprot_id",      type = "Utf8" },
    { name = "transcript_id",   type = "Utf8" },
    { name = "protein_variant", type = "Utf8" },
    { name = "am_pathogenicity",type = "Utf8" },
    { name = "am_class",        type = "Utf8" },
  ]

ingest_sql = """
SELECT chrom,
       CAST(pos AS INT) AS start,
       CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string,
       CAST(am_pathogenicity AS FLOAT) AS am_pathogenicity,
       am_class AS am_class
FROM plugin_alphamissense_src
"""

[[value_columns]]
column = "am_pathogenicity"
csq_field = "am_pathogenicity"
type = "Float32"
[[value_columns]]
column = "am_class"
csq_field = "am_class"
type = "Utf8"

[tier]
threshold = 0.01
unmatched = "cold"
```

- [ ] **Step 2: Validate it parses** by pointing a throwaway `SourceManifest::load` at the file (a `cargo test` in this repo or `cargo run` snippet). Expected: parses, 1 source, 2 value columns.

- [ ] **Step 3: Push the repo.**

```bash
gh repo create biodatageeks/vepyr-plugins --public --source . --push   # or add remote + push
```

- [ ] **Step 4: Record the pinned rev** in this repo's plan notes / a `PLUGINS_REV` constant for reproducibility.

- [ ] **Step 5: Commit** (in `vepyr-plugins`):

```bash
git add . && git commit -m "feat: AlphaMissense source manifest (first plugin)"
```

---

## Task 12: Build driver + build AlphaMissense for one chromosome

**Files:**
- Create: `datafusion/bio-function-vep/examples/build_plugin.rs` (small CLI: `--manifest <toml> --source-path <tsv.gz> --variation-shard <parquet> --out <cache_root> --chrom <c> --af-max-sql <expr>`).

**Interfaces:** Consumes Task 7 `build_plugin_chrom`; writes `plugin/alphamissense/chr22.parquet` + `manifest.json`.

- [ ] **Step 1: Download a chr-scoped slice** (chr22 keeps the cycle fast). AlphaMissense_hg38.tsv.gz is one file; extract chr22 rows to a small local `.tsv.gz`:

```bash
curl -sL https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz \
 | zcat | awk -F'\t' '$1=="chr22" || $1 ~ /^#/' | gzip > /tmp/am_chr22.tsv.gz
```

- [ ] **Step 2: Confirm a variation shard for chr22 exists** (from the existing variation cache). Dump its schema to build `--af-max-sql`:

```bash
# identify the AF columns actually stored, then set e.g.
# --af-max-sql "greatest(coalesce(minor_allele_freq,0), <gnomADg group max expr>)"
```

- [ ] **Step 3: Run the build driver** for chr22:

```bash
cargo run -p datafusion-bio-function-vep --example build_plugin -- \
  --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
  --source-path /tmp/am_chr22.tsv.gz \
  --variation-shard <cache>/variation/chr22.parquet \
  --out /tmp/plugin_cache --chrom 22 \
  --af-max-sql "coalesce(minor_allele_freq,0)"
```

Expected: `plugin/alphamissense/chr22.parquet` + `manifest.json` written; log prints rows/warm/cold.

- [ ] **Step 4: Sanity-check the shard** with a quick `cargo run` or DataFusion query: row count > 0, `tier` present, a spot value matches the TSV for one `(pos, ref/alt)`.

- [ ] **Step 5: Commit** the example + a short `datafusion/bio-function-vep/scripts/build_alphamissense_chr22.sh` documenting the exact commands.

```bash
git add datafusion/bio-function-vep/examples/build_plugin.rs datafusion/bio-function-vep/scripts/build_alphamissense_chr22.sh
git commit -m "feat(plugin-cache): build driver + AlphaMissense chr22 build script"
```

---

## Task 13: Golden parity gate (AlphaMissense chr22)

**Files:**
- Create: `datafusion/bio-function-vep/scripts/parity_alphamissense_chr22.sh` (+ a small comparison harness or reuse the existing golden-VEP comparison used for `fast_chr1_22_merged`).

**Interfaces:** Consumes the built chr22 plugin cache + the annotation wiring (Task 10).

- [ ] **Step 1: Produce our annotation** for a chr22 test VCF with the AlphaMissense plugin enabled (the runtime discovers it via `PluginRegistry`), emitting `am_pathogenicity`/`am_class`.

- [ ] **Step 2: Produce the Ensembl VEP golden** for the same VCF with `--plugin AlphaMissense,AlphaMissense_hg38.tsv.gz` (or reuse an existing golden run if available).

- [ ] **Step 3: Compare only the `am_pathogenicity` and `am_class` fields** per `(chrom,pos,ref,alt,transcript)`. Expected: 100% match on populated rows; where VEP has no AlphaMissense value we also emit empty. Record mismatches like the existing summary reports.

- [ ] **Step 4: Run the w1-vs-w4 body byte-identity check** on the plugin-enabled annotation (existing invariant). Expected: identical bodies.

- [ ] **Step 5: Commit** the parity script + a short results note under `datafusion/bio-function-vep/`.

```bash
git add datafusion/bio-function-vep/scripts/parity_alphamissense_chr22.sh
git commit -m "test(plugin-cache): AlphaMissense chr22 golden parity + worker byte-identity gate"
```

---

## Task 14: `vep-add-plugin` skill (thin workflow wrapper)

**Files:**
- Create: `.claude/skills/vep-add-plugin/SKILL.md` (+ optional `references/` with the manifest template).

**Interfaces:** Documents the Tasks 11→12→13 flow.

- [ ] **Step 1: Write `SKILL.md`** with the decision tree from spec §7: (1) author a TOML source manifest in `vepyr-plugins` (provider + schema + `ingest_sql` → shared key cols + values, `coordinate_system`, value→CSQ, tier); (2) new *format* only → add a provider arm (§6.1 factory); (3) build per-chrom via the driver; (4) parity-gate before production. Include the AlphaMissense manifest as the worked example and exact commands from Tasks 12–13.

- [ ] **Step 2: Verify the skill triggers** (description mentions "add a VEP plugin", "custom plugin cache").

- [ ] **Step 3: Commit.**

```bash
git add .claude/skills/vep-add-plugin
git commit -m "docs(plugin-cache): vep-add-plugin skill (author manifest → build → parity)"
```

---

## Self-Review Notes

- **Spec coverage:** §3.1 naming (Tasks 1–2, 7), §3.2 key columns (Tasks 5, 7), §3.3 contig/coord normalization (Task 3, wired in 7), §4.1 join (Task 4), §4.2 build (Task 7), §4.3 Parquet params (Task 5 via `point_lookup_writer_properties`), §5 runtime lookup + injection (Tasks 8–10), §6.1 source manifest (Task 1), §6.2 cache manifest (Task 6), §6.3 vepyr-plugins repo (Task 11), §7 skill (Task 14), §8 testing/parity (Tasks 4,7,8,13). Scope non-goals (per-gene/interval) intentionally not implemented.
- **Prototype pragmatism:** Task 7 collects per-chrom in memory (fine at chr-scale) instead of the two-pass streaming write the production variation builder uses; noted as a later optimization. Task 8 probes per-variant rather than batched; noted.
- **Known executor unknowns flagged inline:** exact variation AF columns → `--af-max-sql` (Task 4/12), `* EXCLUDE` support (Task 3), the CSQ injection site in `annotate_provider.rs` (Task 10) — each names the file/region to read rather than guessing.
