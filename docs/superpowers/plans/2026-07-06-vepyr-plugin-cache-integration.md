# vepyr Plugin-Cache Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Surface the upstream custom-VEP-plugin-cache subsystem through vepyr's three workflows (cache builder, annotation, full chr1–22 e2e validation), while making it fully declarative (tier sourced from the variation cache; runtime discriminator declared as a manifest template — no per-plugin engine code).

**Architecture:** Upstream (`datafusion-bio-functions`) gains a `PluginCacheBuilder` (owns the per-chrom loop), inherits `tier` from the variation cache instead of computing an AF threshold, and replaces the closed `engine_attr` enum with a compiled template over a generic engine-attribute namespace. vepyr gains a thin `_core.build_plugin_cache` PyO3 wrapper + `vepyr.build_plugin_cache()` Python API (with git-tag manifest resolution from the public `vepyr-plugins` repo), a `plugin_cache_root` kwarg on `annotate()` (riding the existing `options_json` channel), and a full chr1–22 plugin-parity profile.

**Tech Stack:** Rust (DataFusion 53 / Arrow 58, edition 2024, PyO3/maturin), Python 3 (polars, uv), TOML manifests, Ensembl VEP r115 docker (golden), Parquet point-lookup cache.

**Spec:** `docs/superpowers/specs/2026-07-06-vepyr-plugin-cache-integration-design.md`

## Global Constraints

- DataFusion **53.0.0**, Arrow **58.0.0**, Rust edition **2024** — keep in sync across `datafusion-bio-functions`, `datafusion-bio-formats`, `vepyr`.
- `datafusion-bio-function-vep` is built `default-features=false` **without** the `compression` feature (xz2/liblzma link collision with noodles-cram) — plugin sources decompress gzip to temp; do not enable DF `compression`.
- `cargo clippy -- -D warnings` must stay clean; run `cargo fmt` before every commit.
- The plugin-cache module ships under the default `parquet-cache` feature (vepyr's `cache-builder` already pulls it in) — no new Cargo feature.
- **No-plugin runs must stay byte-identical** — `plugin_cache_root: None` disables everything.
- The variation cache layout is `<cache_dir>/variation/<canonical_chrom>.parquet` (dir `variation/`, files `chr1.parquet` … `chr22.parquet`), 1-based `start`/`end`, with a required `tier` Int8 column (0=warm, 1=cold).
- vepyr uses a **local path dep** to `/Users/mwiewior/research/git/datafusion-bio-functions` during development; it MUST be converted to a git tag/rev before any vepyr commit (Phase 10).
- Plugin CSQ output field order = manifest `value_columns` declaration order.

---

## File Structure

**datafusion-bio-functions (upstream) — `datafusion/bio-function-vep/src/plugin_cache/`:**
- `join.rs` — MODIFY: tier inherited from variation (`COALESCE(v.tier,1)`), drop `af_max_sql`/threshold.
- `build.rs` — MODIFY: `build_plugin_chrom` drops `af_max_sql`.
- `source_manifest.rs` — MODIFY: remove `TierPolicy`/`tier`; `MatchColumn.engine_attr` → `template`.
- `cache_manifest.rs` — MODIFY: remove `TierRecord`/`tier`; `MatchColumnRecord.engine_attr` → `template`.
- `template.rs` — CREATE: `ATTR_NAMES`, `CompiledTemplate`, `build_attr_namespace`.
- `registry.rs` — MODIFY: `EngineAttrs` → attribute-namespace array; `PluginEntry` holds compiled templates; `probe_all` evaluates templates.
- `csq.rs` — MODIFY: remove `amino_acid_change`.
- `builder.rs` — CREATE: `PluginCacheBuilder`.
- `mod.rs` — MODIFY: `pub mod template; pub mod builder;`.
- `../annotate_provider.rs` — MODIFY (~line 6211): build the attribute namespace instead of `EngineAttrs`.
- `../../examples/build_plugin.rs` — MODIFY: drop `--af-max-sql`.

**vepyr:**
- `Cargo.toml` — MODIFY: dep pin (path dep during dev).
- `src/lib.rs` — MODIFY: `build_plugin_cache` pyfunction + registration.
- `src/vepyr/__init__.py` — MODIFY: `build_plugin_cache()` + manifest resolution helper; `plugin_cache_root` kwarg on `annotate()`.
- `src/vepyr/_core.pyi` — MODIFY: `build_plugin_cache` stub.
- `docs/plugins.md` — MODIFY: manifest reference + engine-attribute table.
- `e2e-testing/vep-docker.md` — MODIFY: `--plugin AlphaMissense` golden commands.
- `e2e-testing/scripts/run_annotation_fast.py` — MODIFY: plugin profile.

**vepyr-plugins:**
- `plugins/alphamissense/alphamissense.source.toml` — MODIFY: drop `[tier]`; `engine_attr` → `template`; new tag.

---

# Phase 1 — Upstream: tier inherited from the variation cache

### Task 1.1: `tiered_stream`/`tier_sql` inherit `v.tier`; drop `af_max_sql`

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/join.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/build.rs:58-107` (call site) and its test (`:186-283`)

**Interfaces:**
- Produces: `tier_sql(normalized_view: &str, variation_probe: &str) -> String` (threshold arg removed); `tiered_stream(ctx: &SessionContext, normalized_view: &str, variation_shard: &Path) -> Result<SendableRecordBatchStream>` (af_max_sql + threshold removed).

- [ ] **Step 1: Update the `join.rs` unit test to the new behaviour (failing test)**

Replace the test body's variation probe + call so it uses a `tier` column instead of `af_max`, and calls `tier_sql` with two args:

```rust
    #[tokio::test(flavor = "multi_thread")]
    async fn inherits_tier_from_variation_and_miss_is_cold() {
        let ctx = SessionContext::new();
        // variation probe now carries `tier` directly (0=warm, 1=cold).
        let var = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("tier", DataType::Int8, false),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Int8Array::from(vec![0i8, 1])), // warm, cold
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_variation_probe", var).unwrap();

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
        )
        .unwrap();
        ctx.register_batch("plugin_demo_norm", plug).unwrap();

        let df = ctx
            .sql(&tier_sql("plugin_demo_norm", "plugin_variation_probe"))
            .await
            .unwrap()
            .sort(vec![col("start").sort(true, true)])
            .unwrap();
        let batches = df.collect().await.unwrap();
        let b = &batches[0];
        let tier = b
            .column(b.schema().index_of("tier").unwrap())
            .as_any()
            .downcast_ref::<Int8Array>()
            .unwrap();
        assert_eq!(tier.value(0), 0); // 100 -> warm (inherited)
        assert_eq!(tier.value(1), 1); // 200 -> cold (inherited)
        assert_eq!(tier.value(2), 1); // 300 -> no match -> cold
        assert!(b.schema().index_of("demo_score").is_ok());
        assert!(b.schema().index_of("tier").is_ok());
    }
```

Remove the now-unused `Float64Array` import from the test module if clippy flags it.

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::join`
Expected: FAIL (compile error — `tier_sql` still takes 3 args).

- [ ] **Step 3: Implement the new `tier_sql` and `tiered_stream`**

Replace the two functions and the module doc comment in `join.rs`:

```rust
//! Variation-tier inheritance: LEFT-joins the normalized plugin rows against the
//! variation shard on `(chrom, start, allele_string)` and inherits the variation
//! record's `tier` (`COALESCE(v.tier, 1)` — no match → cold). Variation columns
//! drop; only the plugin's value columns plus `tier` survive.

use std::path::Path;

use datafusion::common::Result;
use datafusion::physical_plan::SendableRecordBatchStream;
use datafusion::prelude::SessionContext;

/// SQL that LEFT-joins `normalized_view` to a registered `variation_probe`
/// exposing `(chrom, start, allele_string, tier)` and inherits `tier`. The value
/// columns of `normalized_view` pass through (`p.*`); variation columns drop.
pub fn tier_sql(normalized_view: &str, variation_probe: &str) -> String {
    format!(
        "SELECT p.*, CAST(COALESCE(v.tier, 1) AS TINYINT) AS tier \
         FROM {normalized_view} p \
         LEFT JOIN {variation_probe} v \
         ON p.chrom = v.chrom AND p.start = v.start AND p.allele_string = v.allele_string"
    )
}

/// Register the variation shard as a `tier`-bearing probe view, then stream the
/// tiered rows (tier inherited from the variation record).
pub async fn tiered_stream(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
) -> Result<SendableRecordBatchStream> {
    let shard = variation_shard.to_string_lossy();
    ctx.register_parquet(
        "plugin_variation_raw",
        shard.as_ref(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await?;
    ctx.sql(
        "CREATE OR REPLACE VIEW plugin_variation_probe AS \
         SELECT chrom, start, allele_string, tier FROM plugin_variation_raw",
    )
    .await?;
    let df = ctx
        .sql(&tier_sql(normalized_view, "plugin_variation_probe"))
        .await?;
    df.execute_stream().await
}
```

- [ ] **Step 4: Fix the `build.rs` call site**

In `build.rs`, change the `build_plugin_chrom` signature (remove `af_max_sql: &str,` at line 62) and the `tiered_stream` call (lines 100-107):

```rust
    // Materialize tiered rows (tier inherited from the variation cache).
    let mut stream = tiered_stream(&ctx, &norm_view, variation_shard).await?;
```

- [ ] **Step 5: Update the `build.rs` test to match**

In `build.rs` tests: (a) `write_synthetic_variation` must emit a `tier` column instead of `minor_allele_freq`; (b) the call drops the `af_max_sql` arg. Replace the helper and call:

```rust
    /// Minimal variation-like shard: (chrom, start, allele_string, tier).
    fn write_synthetic_variation(path: &std::path::Path, rows: &[(&str, u32, &str, i8)]) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(rows.iter().map(|r| r.0).collect::<Vec<_>>())),
                Arc::new(UInt32Array::from(rows.iter().map(|r| r.1).collect::<Vec<_>>())),
                Arc::new(StringArray::from(rows.iter().map(|r| r.2).collect::<Vec<_>>())),
                Arc::new(Int8Array::from(rows.iter().map(|r| r.3).collect::<Vec<_>>())),
            ],
        )
        .unwrap();
        let file = std::fs::File::create(path).unwrap();
        let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();
    }
```

Update imports in the test module (`Float64Array` → `Int8Array`). In `builds_tiered_shard_with_counts`, change the variation row to a warm tier and drop the `af_max_sql` argument:

```rust
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]); // 100 warm (tier 0)
        // ...
        let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .unwrap();
        assert_eq!(entry.rows, 2);
        assert_eq!(entry.warm, 1); // start 100 inherited tier 0
        assert_eq!(entry.cold, 1); // start 300 miss -> cold
```

- [ ] **Step 6: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::join plugin_cache::build && cargo clippy -p datafusion-bio-function-vep -- -D warnings`
Expected: PASS, clippy clean.

- [ ] **Step 7: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/src/plugin_cache/join.rs datafusion/bio-function-vep/src/plugin_cache/build.rs
git commit -m "refactor(plugin-cache): inherit tier from variation cache; drop af_max_sql

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

# Phase 2 — Upstream: remove tier/AF from the manifest schema

### Task 2.1: Remove `TierPolicy` (source manifest) and `TierRecord` (cache manifest)

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/cache_manifest.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/build.rs` (test manifest `[tier]` block)
- Modify: `datafusion/bio-function-vep/src/plugin_cache/registry.rs` (test manifest `TierRecord`)

**Interfaces:**
- Produces: `SourceManifest` with no `tier` field; `CacheManifest` with no `tier` field (both compile without any tier types).

- [ ] **Step 1: Edit `source_manifest.rs` test to drop `[tier]` and its assertion (failing test)**

In the `CADD_LIKE` const remove the trailing block:
```toml
[tier]
threshold = 0.01
unmatched = "cold"
```
and delete the assertion line `assert_eq!(m.tier.threshold, 0.01);` from `parses_source_manifest`.

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::source_manifest`
Expected: FAIL (compile error — `SourceManifest` still requires `tier`, and `CADD_LIKE` no longer supplies it → serde parse error at runtime, or `m.tier` still referenced).

- [ ] **Step 3: Remove `TierPolicy` from `source_manifest.rs`**

Delete the `TierPolicy` struct (`:108-115`), `default_threshold` (`:117-119`), and the `pub tier: TierPolicy,` field (`:134`) from `SourceManifest`. Update the module doc comment (line 3) to drop "and tier policy".

- [ ] **Step 4: Remove `TierRecord` from `cache_manifest.rs`**

Delete the `TierRecord` struct (`:13-18`), the `pub tier: TierRecord,` field (`:56`) in `CacheManifest`, and the `tier: TierRecord { … }` block in `from_source` (`:98-101`). In the `cache_manifest.rs` test (`writes_and_discovers_manifest`), delete the `tier: TierRecord { … }` literal. Keep `ChromEntry.warm`/`cold` (still meaningful build counts).

- [ ] **Step 5: Fix the `build.rs` and `registry.rs` test manifests**

- `build.rs` test toml: remove the trailing `[tier]` block.
- `registry.rs` test: remove `TierRecord` from the `use` import (`:194`) and delete the `tier: TierRecord { … }` literal from the `CacheManifest { … }` in `discovers_takes_and_probes`.

- [ ] **Step 6: Run all plugin_cache tests + clippy**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache && cargo clippy -p datafusion-bio-function-vep -- -D warnings`
Expected: PASS, clippy clean.

- [ ] **Step 7: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/src/plugin_cache/
git commit -m "refactor(plugin-cache): remove tier policy from source & cache manifest schema

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

# Phase 3 — Upstream: declarative template discriminator

### Task 3.1: Template engine (`template.rs`)

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/template.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/mod.rs` (add `pub mod template;`)

**Interfaces:**
- Produces:
  - `pub const ATTR_NAMES: &[&str]` (the namespace, fixed order).
  - `pub fn attr_index(name: &str) -> Option<usize>`.
  - `pub struct CompiledTemplate` with `pub fn compile(template: &str) -> Result<CompiledTemplate>` and `pub fn eval(&self, attrs: &[Option<&str>]) -> Option<String>`.
  - `pub fn build_attr_namespace<'a>(consequence: &'a str, gene: &'a str, feature_type: &'a str, feature: &'a str, biotype: &'a str, hgvsc: &'a str, hgvsp: &'a str, cdna_pos: &'a str, cds_pos: &'a str, protein_pos: &'a str, amino_acids: &'a str, codons: &'a str, ref_allele: &'a str, alt_allele: &'a str) -> [Option<&'a str>; ATTR_NAMES.len()]` applying the aa-change validity rules.

- [ ] **Step 1: Write the failing tests**

Create `template.rs` with tests first:

```rust
//! Declarative match-discriminator templates over a fixed engine-attribute
//! namespace. A manifest `[[match_column]].template` (e.g.
//! `"{ref_aa}{Protein_position}{alt_aa}"`) is compiled once per plugin into
//! index-resolved segments and evaluated per transcript consequence. If any
//! referenced attribute is `None`, the discriminator is `None` (probe miss →
//! empty output — the per-transcript gate).

use datafusion::common::{DataFusionError, Result};

/// The engine attributes a template may reference, in fixed order. Extend only
/// when a plugin needs a value not already exposed by the transcript engine.
pub const ATTR_NAMES: &[&str] = &[
    "Consequence", "Gene", "Feature_type", "Feature", "BIOTYPE", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
    "ref_aa", "alt_aa", "ref", "alt",
];

/// Resolve an attribute name to its index in [`ATTR_NAMES`].
pub fn attr_index(name: &str) -> Option<usize> {
    ATTR_NAMES.iter().position(|&n| n == name)
}

#[derive(Debug, Clone, PartialEq)]
enum Segment {
    Lit(String),
    Attr(usize),
}

/// A template compiled to index-resolved segments (no per-record parsing).
#[derive(Debug, Clone)]
pub struct CompiledTemplate {
    segments: Vec<Segment>,
}

impl CompiledTemplate {
    /// Parse a `{name}` template; unknown attribute names are a hard error.
    pub fn compile(template: &str) -> Result<CompiledTemplate> {
        let mut segments = Vec::new();
        let mut rest = template;
        while let Some(open) = rest.find('{') {
            if open > 0 {
                segments.push(Segment::Lit(rest[..open].to_string()));
            }
            let close = rest[open..].find('}').ok_or_else(|| {
                DataFusionError::Execution(format!("unterminated '{{' in template '{template}'"))
            })? + open;
            let name = &rest[open + 1..close];
            let idx = attr_index(name).ok_or_else(|| {
                DataFusionError::Execution(format!("unknown template attribute '{name}'"))
            })?;
            segments.push(Segment::Attr(idx));
            rest = &rest[close + 1..];
        }
        if !rest.is_empty() {
            segments.push(Segment::Lit(rest.to_string()));
        }
        Ok(CompiledTemplate { segments })
    }

    /// Evaluate against a namespace array (same order as [`ATTR_NAMES`]). Any
    /// referenced attribute `None` → whole discriminator `None`.
    pub fn eval(&self, attrs: &[Option<&str>]) -> Option<String> {
        let mut out = String::new();
        for seg in &self.segments {
            match seg {
                Segment::Lit(s) => out.push_str(s),
                Segment::Attr(i) => out.push_str(attrs[*i]?),
            }
        }
        Some(out)
    }
}

/// Build the per-consequence namespace array from the engine's local values,
/// applying the amino-acid-change validity rules that preserve missense gating:
/// `ref_aa`/`alt_aa` are set only for a clean single-residue `X/Y` `amino_acids`;
/// `Protein_position` is passed through only when non-empty and not a range.
#[allow(clippy::too_many_arguments)]
pub fn build_attr_namespace<'a>(
    consequence: &'a str,
    gene: &'a str,
    feature_type: &'a str,
    feature: &'a str,
    biotype: &'a str,
    hgvsc: &'a str,
    hgvsp: &'a str,
    cdna_pos: &'a str,
    cds_pos: &'a str,
    protein_pos: &'a str,
    amino_acids: &'a str,
    codons: &'a str,
    ref_allele: &'a str,
    alt_allele: &'a str,
) -> [Option<&'a str>; ATTR_NAMES.len()] {
    let non_empty = |s: &'a str| if s.is_empty() { None } else { Some(s) };
    // ref_aa/alt_aa only for a clean single-residue X/Y substitution.
    let (ref_aa, alt_aa) = match amino_acids.split_once('/') {
        Some((r, a)) if r.len() == 1 && a.len() == 1 => (Some(r), Some(a)),
        _ => (None, None),
    };
    let protein = match non_empty(protein_pos) {
        Some(p) if !p.contains('-') => Some(p),
        _ => None,
    };
    [
        non_empty(consequence),
        non_empty(gene),
        non_empty(feature_type),
        non_empty(feature),
        non_empty(biotype),
        non_empty(hgvsc),
        non_empty(hgvsp),
        non_empty(cdna_pos),
        non_empty(cds_pos),
        protein,
        non_empty(amino_acids),
        non_empty(codons),
        ref_aa,
        alt_aa,
        non_empty(ref_allele),
        non_empty(alt_allele),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ns_missense() -> [Option<&'static str>; ATTR_NAMES.len()] {
        build_attr_namespace(
            "missense_variant", "ENSG1", "Transcript", "ENST1", "protein_coding",
            "", "", "", "", "320", "W/R", "TGG/CGG", "C", "T",
        )
    }

    #[test]
    fn compiles_and_evaluates_amino_acid_change() {
        let t = CompiledTemplate::compile("{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&ns_missense()).as_deref(), Some("W320R"));
    }

    #[test]
    fn feature_only_template() {
        let t = CompiledTemplate::compile("{Feature}").unwrap();
        assert_eq!(t.eval(&ns_missense()).as_deref(), Some("ENST1"));
    }

    #[test]
    fn non_missense_gates_to_none() {
        // synonymous: no aa change, no protein position
        let ns = build_attr_namespace(
            "synonymous_variant", "ENSG1", "Transcript", "ENST1", "protein_coding",
            "", "", "", "", "", "", "", "C", "T",
        );
        let t = CompiledTemplate::compile("{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&ns), None);
    }

    #[test]
    fn range_position_and_multi_residue_gate() {
        let range = build_attr_namespace(
            "inframe", "", "", "", "", "", "", "", "", "550-551", "VV/A", "", "C", "T",
        );
        let t = CompiledTemplate::compile("{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&range), None);
    }

    #[test]
    fn unknown_attribute_is_error() {
        assert!(CompiledTemplate::compile("{NoSuchAttr}").is_err());
    }
}
```

- [ ] **Step 2: Register the module**

In `datafusion/bio-function-vep/src/plugin_cache/mod.rs` add (keeping alphabetical grouping): `pub mod template;`.

- [ ] **Step 3: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::template`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/src/plugin_cache/template.rs datafusion/bio-function-vep/src/plugin_cache/mod.rs
git commit -m "feat(plugin-cache): compiled template + engine-attribute namespace

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 3.2: Manifest `engine_attr` → `template`

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs` (`MatchColumn`)
- Modify: `datafusion/bio-function-vep/src/plugin_cache/cache_manifest.rs` (`MatchColumnRecord`, `from_source`)

**Interfaces:**
- Produces: `MatchColumn { column: String, template: String }`; `MatchColumnRecord { column: String, template: String }`.

- [ ] **Step 1: Rename the field in `MatchColumn` (`source_manifest.rs:99-106`)**

```rust
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
```

- [ ] **Step 2: Rename in `MatchColumnRecord` + `from_source` (`cache_manifest.rs:41-45,81-88`)**

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchColumnRecord {
    pub column: String,
    pub template: String,
}
```
and in `from_source`:
```rust
            match_columns: src
                .match_columns
                .iter()
                .map(|m| MatchColumnRecord {
                    column: m.column.clone(),
                    template: m.template.clone(),
                })
                .collect(),
```

- [ ] **Step 3: Fix compile errors in dependent tests**

`registry.rs` test (`discovers_takes_and_probes`) constructs `MatchColumn`/`MatchColumnRecord` with `engine_attr` — this is fixed in Task 3.3 (registry rewrite). For now this task's build will fail at `registry.rs`; that is expected and resolved in 3.3. To keep this task independently green, also apply the field rename in the `registry.rs` test literals here (`engine_attr: "amino_acid_change".into()` → `template: "{ref_aa}{Protein_position}{alt_aa}".into()`), and in the `EngineAttrs`-based assertions leave them (3.3 rewrites `probe_all`). If the crate does not compile standalone after this rename, fold Steps of 3.2 and 3.3 into a single commit.

- [ ] **Step 4: Run**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::source_manifest plugin_cache::cache_manifest`
Expected: PASS.

- [ ] **Step 5: Commit** (or defer to end of 3.3 if the crate doesn't compile standalone)

```bash
cargo fmt
git add datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs datafusion/bio-function-vep/src/plugin_cache/cache_manifest.rs
git commit -m "refactor(plugin-cache): manifest match_column engine_attr -> template

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 3.3: `registry.rs` — compile templates, evaluate at probe time

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/registry.rs`

**Interfaces:**
- Consumes: `CompiledTemplate`, `ATTR_NAMES` (Task 3.1); `MatchColumnRecord.template` (Task 3.2).
- Produces: `BufferSlices::probe_all(&self, start: u32, allele_string: &str, attrs: &[Option<&str>]) -> Vec<PluginScalar>` (was `attrs: &EngineAttrs`). `EngineAttrs` deleted.

- [ ] **Step 1: Rewrite the top of `registry.rs`**

Delete the `EngineAttrs` struct + impl (`:14-32`). Add import:
```rust
use crate::plugin_cache::template::CompiledTemplate;
```
Change `PluginEntry.match_engine_attrs: Vec<String>` → `match_templates: Vec<CompiledTemplate>` and `SliceEntry.match_engine_attrs: Vec<String>` → `match_templates: Vec<CompiledTemplate>`.

- [ ] **Step 2: Compile templates in `open()`**

Replace the `match_engine_attrs` construction (`:66-70`) with:
```rust
            let match_templates: Vec<CompiledTemplate> = m
                .match_columns
                .iter()
                .map(|mc| CompiledTemplate::compile(&mc.template))
                .collect::<Result<Vec<_>>>()?;
```
and set `match_templates` on the pushed `PluginEntry` (rename the field there and in `take_buffer_all`'s `SliceEntry` construction: `match_templates: p.match_templates.clone()`).

- [ ] **Step 3: Evaluate templates in `probe_all`**

Replace the `probe_all` signature and body (`:167-187`):
```rust
    pub fn probe_all(
        &self,
        start: u32,
        allele_string: &str,
        attrs: &[Option<&str>],
    ) -> Vec<PluginScalar> {
        let mut out = Vec::new();
        for e in &self.entries {
            let match_values: Vec<Option<String>> =
                e.match_templates.iter().map(|t| t.eval(attrs)).collect();
            let hit = e
                .slice
                .as_ref()
                .and_then(|s| s.probe(start, allele_string, &match_values));
            match hit {
                Some(values) => out.extend(values),
                None => out.extend(std::iter::repeat_n(PluginScalar::Null, e.csq_fields_len)),
            }
        }
        out
    }
```
Update the doc comment above `probe_all` ("built by resolving its match columns' `engine_attr`s" → "built by evaluating its match columns' templates against the namespace").

- [ ] **Step 4: Update the `registry.rs` test to drive templates + namespace**

In `discovers_takes_and_probes`:
- `use` line: replace `MatchColumnRecord, TierRecord` handling per Task 3.2/2.1; import `crate::plugin_cache::template::build_attr_namespace`.
- `MatchColumn`/`MatchColumnRecord` literals use `template: "{ref_aa}{Protein_position}{alt_aa}".into()`.
- Replace the two `EngineAttrs`-based probes with namespace arrays. The shard stores `protein_variant` values `"R12G"`/`"R78G"`, so a hit needs `amino_acids="R/G"`, `protein_pos="78"`:
```rust
        let ns_hit = build_attr_namespace(
            "missense_variant", "", "", "", "", "", "", "", "", "78", "R/G", "", "A", "G",
        );
        let hit = slices.probe_all(100, "A/G", &ns_hit);
        match hit[0] {
            PluginScalar::F32(v) => assert!((v - 0.0427).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        // non-missense (no aa change) → Null (gate)
        let ns_miss = build_attr_namespace(
            "synonymous_variant", "", "", "", "", "", "", "", "", "", "", "", "A", "G",
        );
        assert_eq!(slices.probe_all(100, "A/G", &ns_miss), vec![PluginScalar::Null]);
```

- [ ] **Step 5: Run**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::registry`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/src/plugin_cache/registry.rs
git commit -m "feat(plugin-cache): evaluate compiled templates at probe time; drop EngineAttrs

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 3.4: Probe site — build the namespace; remove `amino_acid_change`

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:6210-6234`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/csq.rs` (remove `amino_acid_change` + its test)

**Interfaces:**
- Consumes: `build_attr_namespace` (Task 3.1), `probe_all(&[Option<&str>])` (Task 3.3).

- [ ] **Step 1: Replace the `EngineAttrs` construction at the probe site**

Replace lines 6211-6217 (the `if plugin_n_fields > 0 { let attrs = EngineAttrs {...}; ... }` head) with a namespace built from the CSQ locals already in scope. Use the same local variable names present in the two `write!(csq_buf, …)` blocks above (`terms_str`, `gene`, `feature_type`, `feature`, `biotype`, `hgvsc`, `hgvsp`, `cdna_pos`, `cds_pos`, `protein_pos`, `amino_acids`, `codons_str`, `input_ref`, `input_alt`):

```rust
                        #[cfg(feature = "parquet-cache")]
                        if plugin_n_fields > 0 {
                            let ns = crate::plugin_cache::template::build_attr_namespace(
                                terms_str, gene, feature_type, feature, biotype, hgvsc, hgvsp,
                                cdna_pos, cds_pos, protein_pos, amino_acids, codons_str,
                                input_ref, input_alt,
                            );
                            let plugin_allele = format!("{input_ref}/{input_alt}");
                            let scalars = plugin_slices
                                .as_ref()
                                .map(|s| s.probe_all(
                                    u32::try_from(start_val).unwrap_or(0),
                                    &plugin_allele,
                                    &ns,
                                ))
                                .unwrap_or_else(|| {
                                    vec![
                                        crate::plugin_cache::lookup::PluginScalar::Null;
                                        plugin_n_fields
                                    ]
                                });
                            csq_buf.push_str(&crate::plugin_cache::csq::field_suffix(&scalars));
```

**Verification note:** the exact local identifiers must match those in scope at `annotate_provider.rs:6138-6205`. Confirm each of `terms_str, gene, feature_type, feature, biotype, hgvsc, hgvsp, cdna_pos, cds_pos, protein_pos, amino_acids, codons_str, input_ref, input_alt` is a `&str`-coercible local there (they are the interpolated CSQ values); adjust any that are `String` with `.as_str()`.

- [ ] **Step 2: Remove `amino_acid_change` from `csq.rs`**

Delete the `amino_acid_change` fn (`csq.rs:15-24`) and the `amino_acid_change_forms_discriminator` test (`:59-69`). Its behaviour now lives in `build_attr_namespace` + `CompiledTemplate` (Task 3.1 tests cover it).

- [ ] **Step 3: Build the whole crate + run plugin tests**

Run: `cargo build -p datafusion-bio-function-vep && cargo test -p datafusion-bio-function-vep --lib plugin_cache && cargo clippy -p datafusion-bio-function-vep -- -D warnings`
Expected: PASS, clippy clean.

- [ ] **Step 4: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/src/annotate_provider.rs datafusion/bio-function-vep/src/plugin_cache/csq.rs
git commit -m "feat(plugin-cache): build engine-attribute namespace at probe site; remove amino_acid_change

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

# Phase 4 — Upstream: `PluginCacheBuilder` + example

### Task 4.1: `PluginCacheBuilder` (owns the chrom loop)

**Files:**
- Create: `datafusion/bio-function-vep/src/plugin_cache/builder.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/mod.rs` (`pub mod builder;`)

**Interfaces:**
- Consumes: `build_plugin_chrom` (Phase 1), `CacheManifest::from_source`/`write` (Phase 2), `ChromEntry`.
- Produces:
```rust
PluginCacheBuilder::new(manifest: &SourceManifest, manifest_file: impl Into<String>,
                        variation_cache_dir: impl Into<PathBuf>, out: impl Into<PathBuf>)
    .with_chrom_filter<I: IntoIterator<Item = S>, S: Into<String>>(chroms: I) -> Self
    .with_overwrite(overwrite: bool) -> Self
    .build_all() -> Result<CacheManifest>   // loops chroms, writes plugin/<name>/manifest.json
```

- [ ] **Step 1: Write the failing test (2-chrom build → 2 manifest entries)**

Create `builder.rs`. The test reuses `build.rs`'s synthetic-source pattern (gzip TSV + a `tier`-bearing variation shard per chrom under `variation/`):

```rust
//! Full plugin-cache build across chromosomes. Owns the per-chrom loop
//! (mirrors `cache_builder::CacheBuilder`): calls `build_plugin_chrom` per
//! chrom against `<variation_cache_dir>/variation/<chrom>.parquet`, accumulates
//! the per-chrom entries into one `CacheManifest`, and writes
//! `plugin/<name>/manifest.json`.

use std::path::{Path, PathBuf};

use datafusion::common::{DataFusionError, Result};

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::build::build_plugin_chrom;
use crate::plugin_cache::cache_manifest::CacheManifest;
use crate::plugin_cache::source_manifest::SourceManifest;

pub struct PluginCacheBuilder<'a> {
    manifest: &'a SourceManifest,
    manifest_file: String,
    variation_cache_dir: PathBuf,
    out: PathBuf,
    chrom_filter: Option<Vec<String>>,
    overwrite: bool,
}

impl<'a> PluginCacheBuilder<'a> {
    pub fn new(
        manifest: &'a SourceManifest,
        manifest_file: impl Into<String>,
        variation_cache_dir: impl Into<PathBuf>,
        out: impl Into<PathBuf>,
    ) -> Self {
        Self {
            manifest,
            manifest_file: manifest_file.into(),
            variation_cache_dir: variation_cache_dir.into(),
            out: out.into(),
            chrom_filter: None,
            overwrite: false,
        }
    }

    pub fn with_chrom_filter<I, S>(mut self, chroms: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let chroms: Vec<String> = chroms.into_iter().map(Into::into).collect();
        self.chrom_filter = (!chroms.is_empty()).then_some(chroms);
        self
    }

    pub fn with_overwrite(mut self, overwrite: bool) -> Self {
        self.overwrite = overwrite;
        self
    }

    /// Per-chrom shard path: `<variation_cache_dir>/variation/<canonical>.parquet`.
    fn variation_shard(&self, chrom: &str) -> PathBuf {
        self.variation_cache_dir
            .join("variation")
            .join(format!("{}.parquet", canonical_chrom_label(chrom)))
    }

    /// Chroms to build: the explicit filter, else every `variation/chr*.parquet`.
    fn resolve_chroms(&self) -> Result<Vec<String>> {
        if let Some(f) = &self.chrom_filter {
            return Ok(f.clone());
        }
        let var_dir = self.variation_cache_dir.join("variation");
        let mut chroms = Vec::new();
        for e in std::fs::read_dir(&var_dir).map_err(|e| {
            DataFusionError::Execution(format!("read {}: {e}", var_dir.display()))
        })? {
            let path = e
                .map_err(|e| DataFusionError::Execution(format!("dir entry: {e}")))?
                .path();
            if path.extension().and_then(|s| s.to_str()) == Some("parquet")
                && let Some(stem) = path.file_stem().and_then(|s| s.to_str())
            {
                chroms.push(stem.to_string());
            }
        }
        chroms.sort();
        Ok(chroms)
    }

    pub async fn build_all(&self) -> Result<CacheManifest> {
        let plugin_dir = self.out.join("plugin").join(&self.manifest.plugin_name);
        let mut cache = CacheManifest::from_source(self.manifest, &self.manifest_file);
        cache.chroms.clear();
        for chrom in self.resolve_chroms()? {
            let shard = self.variation_shard(&chrom);
            if !shard.exists() {
                return Err(DataFusionError::Execution(format!(
                    "variation shard not found: {}",
                    shard.display()
                )));
            }
            let entry =
                build_plugin_chrom(self.manifest, &self.manifest_file, &shard, &self.out, &chrom)
                    .await?;
            cache.chroms.push(entry);
        }
        std::fs::create_dir_all(&plugin_dir).map_err(|e| {
            DataFusionError::Execution(format!("mkdir {}: {e}", plugin_dir.display()))
        })?;
        cache.write(&plugin_dir)?;
        Ok(cache)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int8Array, RecordBatch, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use std::io::Write;
    use std::sync::Arc;

    fn write_gz(path: &Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    fn write_variation(path: &Path, rows: &[(&str, u32, &str, i8)]) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(rows.iter().map(|r| r.0).collect::<Vec<_>>())),
                Arc::new(UInt32Array::from(rows.iter().map(|r| r.1).collect::<Vec<_>>())),
                Arc::new(StringArray::from(rows.iter().map(|r| r.2).collect::<Vec<_>>())),
                Arc::new(Int8Array::from(rows.iter().map(|r| r.3).collect::<Vec<_>>())),
            ],
        )
        .unwrap();
        let file = std::fs::File::create(path).unwrap();
        let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn builds_all_chroms_and_writes_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("src.tsv.gz");
        // chr1: pos 100 ; chr2: pos 200
        write_gz(&src, "chr1\t100\tA\tG\t0.9\nchr2\t200\tC\tT\t0.5\n");

        let var_dir = dir.path().join("cache").join("variation");
        std::fs::create_dir_all(&var_dir).unwrap();
        write_variation(&var_dir.join("chr1.parquet"), &[("1", 100, "A/G", 0)]);
        write_variation(&var_dir.join("chr2.parquet"), &[("2", 200, "C/T", 1)]);

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            src.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", dir.path().join("cache"), &out)
            .with_chrom_filter(["1", "2"])
            .build_all()
            .await
            .unwrap();
        assert_eq!(cache.chroms.len(), 2);
        assert!(out.join("plugin").join("demo").join("manifest.json").exists());
        assert!(out.join("plugin").join("demo").join("chr1.parquet").exists());
        assert!(out.join("plugin").join("demo").join("chr2.parquet").exists());
    }
}
```

- [ ] **Step 2: Register the module**

In `plugin_cache/mod.rs` add `pub mod builder;`.

- [ ] **Step 3: Run**

Run: `cargo test -p datafusion-bio-function-vep --lib plugin_cache::builder`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/src/plugin_cache/builder.rs datafusion/bio-function-vep/src/plugin_cache/mod.rs
git commit -m "feat(plugin-cache): PluginCacheBuilder owns the per-chrom build loop

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 4.2: Update `build_plugin` example

**Files:**
- Modify: `datafusion/bio-function-vep/examples/build_plugin.rs`

- [ ] **Step 1: Rewrite the example to use `PluginCacheBuilder`**

Replace the argument parsing + `build_plugin_chrom` call so it drops `--af-max-sql`, takes `--variation-cache-dir` (the cache root, not a single shard) and optional `--chrom` (repeatable → filter; absent → all), and calls the builder. Minimal version:

```rust
//! Build a plugin cache from a source manifest (all chroms, or a filtered set).
//!
//! ```text
//! cargo run -p datafusion-bio-function-vep --example build_plugin -- \
//!   --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
//!   --source-path /tmp/AlphaMissense_hg38.tsv.gz \
//!   --variation-cache-dir <cache root containing variation/> \
//!   --out /tmp/plugin_cache [--chrom 22]
//! ```

use std::path::PathBuf;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::plugin_cache::builder::PluginCacheBuilder;
use datafusion_bio_function_vep::plugin_cache::source_manifest::SourceManifest;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1)).cloned()
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let manifest_path = arg(&args, "--manifest")
        .ok_or_else(|| DataFusionError::Execution("--manifest required".into()))?;
    let variation_cache_dir = PathBuf::from(
        arg(&args, "--variation-cache-dir")
            .ok_or_else(|| DataFusionError::Execution("--variation-cache-dir required".into()))?,
    );
    let out = PathBuf::from(
        arg(&args, "--out").ok_or_else(|| DataFusionError::Execution("--out required".into()))?,
    );

    let mut manifest = SourceManifest::load(&PathBuf::from(&manifest_path))?;
    if let Some(source_path) = arg(&args, "--source-path")
        && let Some(first) = manifest.sources.first_mut()
    {
        first.path = source_path;
    }
    let manifest_file = PathBuf::from(&manifest_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| manifest_path.clone());

    let mut builder =
        PluginCacheBuilder::new(&manifest, &manifest_file, &variation_cache_dir, &out);
    if let Some(chrom) = arg(&args, "--chrom") {
        builder = builder.with_chrom_filter([chrom]);
    }
    let cache = builder.build_all().await?;
    for c in &cache.chroms {
        println!("  {} rows={} warm={} cold={}", c.chrom, c.rows, c.warm, c.cold);
    }
    println!("  wrote plugin/{}/manifest.json", manifest.plugin_name);
    Ok(())
}
```

- [ ] **Step 2: Build the example**

Run: `cargo build -p datafusion-bio-function-vep --example build_plugin`
Expected: compiles.

- [ ] **Step 3: Full upstream test sweep**

Run: `cargo test -p datafusion-bio-function-vep && cargo clippy -p datafusion-bio-function-vep --all-targets -- -D warnings`
Expected: all green (the 803+ crate/integration tests).

- [ ] **Step 4: Commit**

```bash
cargo fmt
git add datafusion/bio-function-vep/examples/build_plugin.rs
git commit -m "chore(plugin-cache): build_plugin example uses PluginCacheBuilder (all chroms)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

# Phase 5 — vepyr-plugins: AlphaMissense manifest update

### Task 5.1: Update the AM manifest to the new schema (local edit)

**Files:**
- Modify: `/Users/mwiewior/workspace/vepyr-plugins/plugins/alphamissense/alphamissense.source.toml`

- [ ] **Step 1: Edit the `[[match_column]]` block (lines 54-56)**

```toml
[[match_column]]
column   = "protein_variant"
template = "{ref_aa}{Protein_position}{alt_aa}"
```

- [ ] **Step 2: Delete the `[tier]` block (lines 69-71)** entirely.

- [ ] **Step 3: Update the header comment** (lines 6-10) to note tiering is inherited from the variation cache and matching uses a template (drop the `engine_attr`/`transcript_match` phrasing as needed).

- [ ] **Step 4: Verify it parses against the new upstream schema**

From the datafusion-bio-functions repo, add a throwaway assertion or reuse the example's `SourceManifest::load`:

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo run -q -p datafusion-bio-function-vep --example build_plugin -- --manifest /Users/mwiewior/workspace/vepyr-plugins/plugins/alphamissense/alphamissense.source.toml --variation-cache-dir /nonexistent --out /tmp/x 2>&1 | head`
Expected: it should get *past* manifest parsing (fail later on the missing variation dir), proving the TOML parses under the new schema. A parse error here means the manifest still has stale fields.

- [ ] **Step 5: Stage but DO NOT push yet**

```bash
cd /Users/mwiewior/workspace/vepyr-plugins
git add plugins/alphamissense/alphamissense.source.toml
git commit -m "feat(alphamissense): template discriminator; drop tier block (variation-inherited)"
```
**Do not `git push` or tag here.** The tagged release + push happen in Phase 10 (Step: vepyr-plugins release), confirmed with the user, so no manifest tag lands before the upstream code that reads it.

---

# Phase 6 — vepyr: cache builder

### Task 6.1: Point vepyr at the local upstream (dev path dep)

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/Cargo.toml:46`

- [ ] **Step 1: Replace the tag pin with a path dep**

```toml
datafusion-bio-function-vep = { path = "/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep", features = ["cache-builder"] }
```

- [ ] **Step 2: Rebuild the extension**

Run: `cd /Users/mwiewior/research/git/vepyr && uv sync --reinstall-package vepyr`
Expected: builds `_core` against local upstream. (If the two `datafusion-bio-format-*` tags conflict with the local upstream's pinned formats revs, align them — see Global Constraints.)

- [ ] **Step 3: Do NOT commit `Cargo.toml`** yet — this path dep is reverted to a tag in Phase 10. Note it in the working tree only.

### Task 6.2: `_core.build_plugin_cache` PyO3 wrapper

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/src/lib.rs`

**Interfaces:**
- Consumes: `PluginCacheBuilder`, `SourceManifest` (upstream).
- Produces: `_core.build_plugin_cache(manifest_path, source_path, variation_cache_dir, plugin_cache_root, chroms=None, overwrite=False) -> list[(chrom, rows, warm, cold)]`.

- [ ] **Step 1: Add imports (lib.rs top)**

```rust
use datafusion_bio_function_vep::plugin_cache::builder::PluginCacheBuilder;
use datafusion_bio_function_vep::plugin_cache::source_manifest::SourceManifest;
```

- [ ] **Step 2: Add the pyfunction (after `build_cache`)**

```rust
/// Build a plugin cache (all chroms, or a filtered set) from a source manifest.
/// Returns per-chrom `(chrom, rows, warm, cold)` tuples.
#[pyfunction]
#[pyo3(signature = (manifest_path, source_path, variation_cache_dir, plugin_cache_root, chroms=None, overwrite=false))]
fn build_plugin_cache(
    py: Python<'_>,
    manifest_path: &str,
    source_path: &str,
    variation_cache_dir: &str,
    plugin_cache_root: &str,
    chroms: Option<Vec<String>>,
    overwrite: bool,
) -> PyResult<Vec<(String, usize, usize, usize)>> {
    let mut manifest = SourceManifest::load(std::path::Path::new(manifest_path))
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("load manifest: {e}")))?;
    if let Some(first) = manifest.sources.first_mut() {
        first.path = source_path.to_string();
    }
    let manifest_file = std::path::Path::new(manifest_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| manifest_path.to_string());

    let rt = tokio::runtime::Builder::new_multi_thread()
        .enable_all()
        .build()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let cache = py.detach(|| {
        rt.block_on(async {
            let mut b = PluginCacheBuilder::new(
                &manifest,
                &manifest_file,
                variation_cache_dir,
                plugin_cache_root,
            )
            .with_overwrite(overwrite);
            if let Some(cs) = chroms {
                b = b.with_chrom_filter(cs);
            }
            b.build_all().await
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("plugin build failed: {e}")))
    })?;

    Ok(cache
        .chroms
        .into_iter()
        .map(|c| (c.chrom, c.rows, c.warm, c.cold))
        .collect())
}
```

- [ ] **Step 3: Register it in `_core` (line 150 area)**

```rust
    m.add_function(wrap_pyfunction!(build_plugin_cache, m)?)?;
```

- [ ] **Step 4: Add the `.pyi` stub (`_core.pyi`)**

```python
def build_plugin_cache(
    manifest_path: str,
    source_path: str,
    variation_cache_dir: str,
    plugin_cache_root: str,
    chroms: list[str] | None = None,
    overwrite: bool = False,
) -> list[tuple[str, int, int, int]]:
    """Build a plugin cache from a source manifest. Returns (chrom, rows, warm, cold)."""
    ...
```

- [ ] **Step 5: Build**

Run: `cd /Users/mwiewior/research/git/vepyr && uv sync --reinstall-package vepyr`
Expected: compiles, `_core.build_plugin_cache` importable.

- [ ] **Step 6: Commit (vepyr repo; Cargo.toml stays unstaged)**

```bash
cd /Users/mwiewior/research/git/vepyr
git add src/lib.rs src/vepyr/_core.pyi
git commit -m "feat: _core.build_plugin_cache PyO3 wrapper over PluginCacheBuilder"
```

### Task 6.3: `vepyr.build_plugin_cache()` Python API + manifest resolution

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/src/vepyr/__init__.py`
- Test: `/Users/mwiewior/research/git/vepyr/tests/test_build_plugin_cache.py`

**Interfaces:**
- Consumes: `_core.build_plugin_cache` (Task 6.2).
- Produces: `vepyr.build_plugin_cache(plugin, version, source_path, cache_dir, plugin_cache_root, chroms=None, plugins_repo=None, overwrite=False) -> list[tuple[str,int,int,int]]` and `_resolve_plugin_manifest(plugin, version, plugins_repo=None, repo_url=DEFAULT) -> str`.

- [ ] **Step 1: Write the failing test for manifest resolution (offline path)**

Create `tests/test_build_plugin_cache.py`. The offline branch (a provided `plugins_repo` clone) is unit-testable without network — seed a tiny git repo with a tag:

```python
import subprocess
from pathlib import Path

import vepyr


def _init_plugins_repo(root: Path) -> Path:
    repo = root / "vepyr-plugins"
    (repo / "plugins" / "demo").mkdir(parents=True)
    (repo / "plugins" / "demo" / "demo.source.toml").write_text('plugin_name = "demo"\n')
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t",
                    "commit", "-qm", "init"], cwd=repo, check=True)
    subprocess.run(["git", "tag", "v0.1.0"], cwd=repo, check=True)
    (repo / "plugins" / "demo" / "demo.source.toml").write_text('plugin_name = "demo2"\n')
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t",
                    "commit", "-qm", "v2"], cwd=repo, check=True)
    return repo


def test_resolve_manifest_offline_checks_out_tag(tmp_path):
    repo = _init_plugins_repo(tmp_path)
    path = vepyr._resolve_plugin_manifest("demo", "v0.1.0", plugins_repo=str(repo))
    assert Path(path).read_text().strip() == 'plugin_name = "demo"'  # the tagged version
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /Users/mwiewior/research/git/vepyr && uv run pytest tests/test_build_plugin_cache.py -q`
Expected: FAIL (`_resolve_plugin_manifest` undefined).

- [ ] **Step 3: Implement resolution + the Python API**

Add near the top of `__init__.py` (after imports):

```python
import subprocess
import tempfile

DEFAULT_PLUGINS_REPO_URL = "https://github.com/biodatageeks/vepyr-plugins.git"


def _resolve_plugin_manifest(
    plugin: str,
    version: str,
    *,
    plugins_repo: str | None = None,
    repo_url: str = DEFAULT_PLUGINS_REPO_URL,
) -> str:
    """Resolve `plugins/<plugin>/<plugin>.source.toml` at git tag `version`.

    Offline: reuse a provided local clone (`plugins_repo`). Online: clone the
    public repo into a temp dir. Either way, materialize the file at `version`
    via `git worktree` (never disturbs the caller's checkout).
    """
    if plugins_repo is not None:
        repo = plugins_repo
    else:
        repo = tempfile.mkdtemp(prefix="vepyr-plugins-")
        subprocess.run(["git", "clone", "--quiet", repo_url, repo], check=True)
    worktree = tempfile.mkdtemp(prefix="vepyr-plugins-wt-")
    subprocess.run(
        ["git", "-C", repo, "worktree", "add", "--quiet", "--detach", worktree, version],
        check=True,
    )
    rel = os.path.join("plugins", plugin, f"{plugin}.source.toml")
    manifest = os.path.join(worktree, rel)
    if not os.path.exists(manifest):
        raise FileNotFoundError(f"{rel} not found at {version} in {repo}")
    return manifest


def build_plugin_cache(
    plugin: str,
    version: str,
    *,
    source_path: str,
    cache_dir: str,
    plugin_cache_root: str,
    chroms: list[str] | None = None,
    plugins_repo: str | None = None,
    overwrite: bool = False,
) -> list[tuple[str, int, int, int]]:
    """Build a per-chromosome plugin cache.

    `plugin`/`version` select `plugins/<plugin>/<plugin>.source.toml` from the
    public vepyr-plugins repo at that git tag (or `plugins_repo` for offline).
    Tiering is inherited from the variation cache at `cache_dir`.
    """
    manifest_path = _resolve_plugin_manifest(plugin, version, plugins_repo=plugins_repo)
    return _core.build_plugin_cache(
        manifest_path,
        source_path,
        cache_dir,
        plugin_cache_root,
        chroms,
        overwrite,
    )
```

Ensure `from . import _core` (or the existing import alias used for `_build_cache`/`_annotate_vcf`) exposes `_core.build_plugin_cache`; match the existing import style in the file.

- [ ] **Step 4: Run the test**

Run: `cd /Users/mwiewior/research/git/vepyr && uv run pytest tests/test_build_plugin_cache.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /Users/mwiewior/research/git/vepyr
git add src/vepyr/__init__.py tests/test_build_plugin_cache.py
git commit -m "feat: vepyr.build_plugin_cache() with git-tag manifest resolution"
```

---

# Phase 7 — vepyr: annotation `plugin_cache_root`

### Task 7.1: Thread `plugin_cache_root` through annotate

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/src/annotate.rs` (config build ~lines 188-286)
- Modify: `/Users/mwiewior/research/git/vepyr/src/vepyr/__init__.py` (`annotate` signature + opts dict)
- Test: `/Users/mwiewior/research/git/vepyr/tests/test_annotate.py` (add a plugin-root option test)

**Interfaces:**
- Consumes: `AnnotateVcfConfig.plugin_cache_root: Option<PathBuf>` (upstream, present).
- Produces: `annotate(..., plugin_cache_root: str | None = None)`; `opts["plugin_cache_root"]` when set.

- [ ] **Step 1: Read the option in `annotate.rs`**

In the `AnnotateVcfConfig { … }` literal add (after `on_batch_written: callback,`):

```rust
        plugin_cache_root: opts
            .get("plugin_cache_root")
            .and_then(|v| v.as_str())
            .map(std::path::PathBuf::from),
```

(Confirm the field name matches upstream `vcf_sink::AnnotateVcfConfig`; it is `plugin_cache_root: Option<PathBuf>`.)

- [ ] **Step 2: Add the kwarg + opts entry in `__init__.py`**

In `annotate(...)` signature add near the output-mode block: `plugin_cache_root: str | None = None,`. In the opts-dict construction (after the `workers` block, ~line 726) add:

```python
    if plugin_cache_root is not None:
        opts["plugin_cache_root"] = plugin_cache_root
```

- [ ] **Step 3: Write a round-trip test (options carry the key)**

Add to `tests/test_annotate.py` a test that a plugin root reaches `options_json`. If the suite has no seam to inspect `options_json` directly, assert behaviorally against a tiny fixture cache: annotate a 1-variant VCF with and without `plugin_cache_root` and assert the CSQ header/body gains the plugin field only when set. Minimal seam-based version if `opts` building is refactorable — otherwise use the fixture approach used elsewhere in `tests/`.

```python
def test_plugin_cache_root_reaches_options(monkeypatch):
    captured = {}
    import vepyr._core as core
    monkeypatch.setattr(core, "annotate_vcf",
                        lambda *a, **k: captured.setdefault("options_json", a[3]) or 0)
    vepyr.annotate("in.vcf", "cache", output_vcf="out.vcf",
                   plugin_cache_root="/tmp/pc")
    assert '"plugin_cache_root":"/tmp/pc"' in captured["options_json"]
```

(Adjust the positional index `a[3]` to match `_core.annotate_vcf(vcf, cache, out, options_json, …)`.)

- [ ] **Step 4: Build + run**

Run: `cd /Users/mwiewior/research/git/vepyr && uv sync --reinstall-package vepyr && uv run pytest tests/test_annotate.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /Users/mwiewior/research/git/vepyr
git add src/annotate.rs src/vepyr/__init__.py tests/test_annotate.py
git commit -m "feat: plugin_cache_root kwarg on vepyr.annotate() (rides options_json)"
```

---

# Phase 8 — vepyr: documentation

### Task 8.1: Manifest reference + engine-attribute table in `docs/plugins.md`

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/docs/plugins.md`

- [ ] **Step 1: Add a "Building a plugin cache" section**

Document `vepyr.build_plugin_cache(plugin, version, source_path=…, cache_dir=…, plugin_cache_root=…, chroms=None, plugins_repo=None)` and `annotate(..., plugin_cache_root=…)`, with the AlphaMissense end-to-end example (build → annotate).

- [ ] **Step 2: Add the "Manifest structure" reference (spec §11.1)**

Document each section of `<plugin>.source.toml`: top-level scalars (`plugin_name`, `coordinate_system`, `ingest_sql`) + the TOML ordering rule; `ingest_sql` required key columns (`chrom`, `start`, `end`, `allele_string`) + discriminator + value columns; `[[source]]`/`[source.csv]`; `[[match_column]]` (`column` + `template`); `[[value_columns]]` (declaration order = CSQ order); and the note that there is **no `[tier]`** block (tier is inherited from the variation cache). Include the full AlphaMissense manifest and a minimal per-variant example.

- [ ] **Step 3: Add the "Engine-attribute namespace" table (spec §11.2)**

Copy the table from spec §11.2 verbatim (the 16 attributes + descriptions), the None-propagation/gating rule, and the two worked template examples (`{ref_aa}{Protein_position}{alt_aa}`, `{Feature}`). Keep the list in sync with `template::ATTR_NAMES`.

- [ ] **Step 4: Commit**

```bash
cd /Users/mwiewior/research/git/vepyr
git add docs/plugins.md
git commit -m "docs: plugin manifest structure reference + engine-attribute table"
```

---

# Phase 9 — E2E: full chr1–22 AlphaMissense parity

### Task 9.1: Golden VEP `--plugin AlphaMissense` generation

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/e2e-testing/vep-docker.md`

- [ ] **Step 1: Add the golden command block**

Document the `docker run … ensemblorg/ensembl-vep:release_115.2 … --everything --hgvs --plugin AlphaMissense,<data>` invocation (mounting the AlphaMissense tabix + `.tbi`) producing a chr1–22 golden VCF whose CSQ carries `am_class`/`am_pathogenicity`. Reference the existing merged-cache golden runbook.

- [ ] **Step 2: Generate the golden VCF** (operational — requires docker + AM data)

Run the documented docker command to produce `HG002_..._merged_am.vcf` under `DATA_DIR`. Record the output path.

- [ ] **Step 3: Commit the runbook change**

```bash
cd /Users/mwiewior/research/git/vepyr
git add e2e-testing/vep-docker.md
git commit -m "docs(e2e): AlphaMissense --plugin golden generation for chr1-22"
```

### Task 9.2: Full-genome AlphaMissense plugin cache build

- [ ] **Step 1: Build all chroms** (operational)

```python
import vepyr
vepyr.build_plugin_cache(
    "alphamissense", "<vepyr-plugins tag>",
    source_path=".../AlphaMissense_hg38.tsv.gz",
    cache_dir="/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged",
    plugin_cache_root="/Users/mwiewior/workspace/data_vepyr/plugin_cache",
    chroms=[str(n) for n in range(1, 23)],
    plugins_repo="/Users/mwiewior/workspace/vepyr-plugins",  # offline, pre-release
)
```
Expected: `plugin/alphamissense/chr1.parquet … chr22.parquet` + `manifest.json`.

### Task 9.3: Plugin profile in `run_annotation_fast.py`

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/e2e-testing/scripts/run_annotation_fast.py`

- [ ] **Step 1: Add a profile entry** to `_CACHE_PROFILES` (after `merged`), threading the plugin root through `annotate_kwargs` (flows via `**args.annotate_kwargs` at the `vepyr.annotate(...)` call, line 731):

```python
    "merged_am": {
        "cache_dir": os.path.join(DATA_DIR, "115_GRCh38_merged"),
        "vep_vcf": os.path.join(
            DATA_DIR, "HG002_annotated_wgs_everything_hgvs_merged_am.vcf"
        ),
        "annotate_kwargs": {
            "plugin_cache_root": os.path.join(DATA_DIR, "plugin_cache"),
        },
        "suffix": "_merged_am",
    },
```

- [ ] **Step 2: Sanity-run one chromosome**

Run: `cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts && uv run python run_annotation_fast.py chr22 --profile merged_am`
Expected: a `reports/fast_chr22_merged_am_report.json` whose `am_class`/`am_pathogenicity` field mismatch counts are **0** (and no over-emit on non-missense lines — the gate).

- [ ] **Step 3: Commit**

```bash
cd /Users/mwiewior/research/git/vepyr
git add e2e-testing/scripts/run_annotation_fast.py
git commit -m "test(e2e): merged_am plugin-parity profile (AlphaMissense)"
```

### Task 9.4: Full chr1–22 parity gate

- [ ] **Step 1: Run all chromosomes** (operational)

Run: `cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts && uv run python run_annotation_fast_all.py --profile merged_am`
Expected: the aggregated `reports/fast_chr1_22_merged_am_summary_*.md` shows **0 mismatches and 0 over-emit** on populated `am_class`/`am_pathogenicity` CSQ lines across chr1–22 (extending the proven chr22 result genome-wide). Investigate any non-zero field before proceeding.

---

# Phase 10 — Release & dependency cutover

### Task 10.1: Upstream release

**Files:**
- `datafusion-bio-functions` (branch `feat/plugin-cache-alphamissense` / PR #190)

- [ ] **Step 1: Full test + clippy sweep** on the upstream branch.

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo test -p datafusion-bio-function-vep && cargo clippy --all-targets -- -D warnings`
Expected: green.

- [ ] **Step 2: Merge PR #190 to `master`** (via the normal review/merge flow), then cut and push a new tag:

```bash
git checkout master && git pull
git tag v0.14.0
git push origin v0.14.0
```
(Confirm the version bump convention with the maintainer before tagging.)

### Task 10.2: vepyr dependency cutover

**Files:**
- Modify: `/Users/mwiewior/research/git/vepyr/Cargo.toml:46`

- [ ] **Step 1: Replace the path dep with the new tag**

```toml
datafusion-bio-function-vep = { git = "https://github.com/biodatageeks/datafusion-bio-functions.git", tag = "v0.14.0", features = ["cache-builder"] }
```

- [ ] **Step 2: Rebuild + full vepyr test suite**

Run: `cd /Users/mwiewior/research/git/vepyr && uv sync --reinstall-package vepyr && uv run pytest -q`
Expected: green (golden suite + build/annotate tests).

- [ ] **Step 3: Commit the pin**

```bash
cd /Users/mwiewior/research/git/vepyr
git add Cargo.toml
git commit -m "chore: pin datafusion-bio-function-vep v0.14.0 (plugin cache)"
```

### Task 10.3: vepyr-plugins tagged release (confirm before push)

**Files:**
- `/Users/mwiewior/workspace/vepyr-plugins`

- [ ] **Step 1: Confirm with the user** that the AM manifest change (Phase 5) is ready to publish, then tag + push:

```bash
cd /Users/mwiewior/workspace/vepyr-plugins
git tag <plugin release tag>   # e.g. alphamissense-v0.2.0 or v0.2.0
git push origin master --tags
```
(This is the outward, versioned release — do not push until the upstream code that reads the new schema is tagged/consumed, i.e. after Task 10.1.)

- [ ] **Step 2: Update the e2e build to use the released tag**

Change the `version=` / drop `plugins_repo=` in the Task 9.2 build command to consume the public tag online, and re-verify one chromosome (Task 9.3 Step 2).

---

## Self-Review Notes

- **Spec coverage:** §3 dep→Phase 6.1/10; §5 builder→Phase 4 + 6; §6 tier→Phase 1; §7.1 schema→Phase 2; §7.2 template→Phase 3; §7.4 repush→Phase 5 + 10.3; §8 annotate→Phase 7; §9 e2e→Phase 9; §11 docs→Phase 8. All covered.
- **Type consistency:** `tier_sql`/`tiered_stream` new arities used consistently in build.rs; `MatchColumn.template`/`MatchColumnRecord.template` consistent across source_manifest/cache_manifest/registry; `probe_all(&[Option<&str>])` matches `build_attr_namespace` return; `PluginCacheBuilder` API identical in builder.rs, example, and `_core.build_plugin_cache`.
- **Ordering caveat:** if the upstream crate does not compile standalone between Task 3.2 and 3.3 (the `engine_attr`→`template` rename transiently breaks `registry.rs`), squash 3.2+3.3 into one commit — noted inline.
