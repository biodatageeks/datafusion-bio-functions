# Plugin Factory — Plan A: Engine (VCF provider wiring + manifest hardening)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `plugin_cache` able to ingest VCF-sourced plugins and make a malformed manifest fail loudly instead of silently, so that Wave-1 plugin ports really are "manifest only".

**Architecture:** `plugin_cache` already registers CSV/TSV/Parquet sources through DataFusion's builtin providers and returns `NotImplemented` for `vcf`/`bed`. The `bio-function-vep` crate **already depends on `datafusion-bio-format-vcf`** (`Cargo.toml:33`, tag `v1.8.8`), which ships a complete `VcfTableProvider` with INFO projection and BGZF support. So this is **wiring, not implementation**. Alongside it, four correctness fixes to the manifest/driver surface that currently fail silently or lie.

**Tech Stack:** Rust, DataFusion, serde/toml, `datafusion-bio-format-vcf`.

**Source spec:** `docs/superpowers/specs/2026-07-13-vep-plugin-port-factory-design.md` (§7).

---

## Scope note (deviation from the spec, deliberate)

Spec §7.1 says wire **both** `vcf` and `bed`. This plan wires **`vcf` only**.

Reason: no Wave-0 client needs BED, and `bio-format-bed` is **not** currently a dependency of
`bio-function-vep` (only `bio-format-vcf` is), so wiring it means adding a dependency for zero
callers. YAGNI. `Bed` keeps returning `NotImplemented`, but Task 4 makes that error honest.
When a BED-sourced plugin appears (Wave 2), wiring it is a copy of Task 3.

---

## Prerequisite 0 — sync `dev-test` with `master` (blocks everything)

**RESOLVED — this prerequisite is obsolete. Read it anyway: it explains why the base branch is what
it is.**

The original plan said "sync `dev-test` with `master`; expect no conflicts, `dev-test` is strictly
behind." **That was wrong.** The branches diverged at `v0.10.0` and both moved:

- `master` took `a6e19ad` — *Lance-only cache, grid-aligned parallel annotation, SIFT v2 + AF
  zero-copy, engine-redundancy trims (#181)* — which rewrote `annotate_provider.rs`
  (**+6907/−2446**) and **deleted** `variant_lookup_exec.rs`.
- `dev-test` has **45 commits** of its own, including a real engine feature that exists **nowhere
  else**: multi-ALT CSQ per-allele expansion (`PerAltCtx`, `vcf_to_vep_allele_multi`). `master` still
  treats `alt` as a single allele (`annotate_provider.rs:5325`).

So merging them is not conflict resolution — it is **porting a feature into a rewritten hot path**.
149 files merge cleanly; the 5 conflicting hunks are all CSQ per-allele assembly. Get it wrong and
nothing crashes: every multi-allelic variant is just silently mis-annotated. That deserves its own PR
and its own parity test on multi-allelic sites, and it has **nothing to do with the plugin factory**.

**Decision (2026-07-13): decouple.**

- Base branch for all factory work is **`master-sitekwb`** — cut from `master`, so it already
  contains `plugin_cache` (14 files). It is treated as our `main`: every PR targets it.
- **Never commit to `master` or `main`.**
- The multi-ALT port from `dev-test` is tracked separately and does not block anything here.

---

## Working branch for Tasks 1–7

```bash
cd /Users/wojtek/Documents/vepyr/datafusion-bio-functions
git fetch origin
git worktree add ../datafusion-bio-functions-worktrees/plugin-engine \
  -b feat/plugin-cache-vcf-provider origin/master-sitekwb
cd ../datafusion-bio-functions-worktrees/plugin-engine
```

Every PR in this plan targets `master-sitekwb` (`gh pr create --base master-sitekwb ...`).

All paths below are relative to that worktree.

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs` | The TOML manifest schema — the contract with `vepyr-plugins` | Modify: `deny_unknown_fields`, new `VcfParams`, new `SourceSpec.vcf` field |
| `datafusion/bio-function-vep/src/plugin_cache/provider.rs` | Registers each source as a DataFusion table | Modify: wire `ProviderKind::Vcf`; sharpen the `Bed` error |
| `datafusion/bio-function-vep/src/plugin_cache/builder.rs` | Per-chrom build loop, manifest UPSERT | Modify: make `overwrite` mean something |
| `datafusion/bio-function-vep/examples/build_plugin.rs` | The build driver | Modify: `--source-path` for all sources, repeatable `--chrom` |
| `datafusion/bio-function-vep/Cargo.toml` | Deps | Modify: add `datafusion-bio-format-core` (for `ObjectStorageOptions`) |

---

### Task 1: Make an unknown manifest key a hard error

Today `SourceManifest` has no `deny_unknown_fields`. A typo, or the `[tier]` block that the design
docs still show and the struct does not have, is **silently ignored** — the manifest "works" and
quietly does something other than what it says. Proof that this already bit us: the existing test in
`provider.rs:173-174` carries a live `[tier]` block that does nothing.

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs`
- Modify: `datafusion/bio-function-vep/src/plugin_cache/provider.rs:146-177` (test fixture)
- Test: same files (`#[cfg(test)] mod tests`)

- [ ] **Step 1: Write the failing test**

Add to the `mod tests` block at the bottom of `source_manifest.rs`:

```rust
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
```

- [ ] **Step 2: Run it and watch it fail**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache \
  plugin_cache::source_manifest::tests::rejects_unknown_key_instead_of_ignoring_it
```
Expected: FAIL — `unknown key [tier] must be rejected` (the parse currently succeeds).

- [ ] **Step 3: Add `deny_unknown_fields` to every manifest struct**

In `source_manifest.rs`, add the attribute to each of these six structs (keep all existing
attributes; just add one line under the `#[derive(...)]`):

```rust
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SchemaField { /* unchanged */ }

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CsvParams { /* unchanged */ }

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceSpec { /* unchanged */ }

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValueColumn { /* unchanged */ }

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatchColumn { /* unchanged */ }

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceManifest { /* unchanged */ }
```

- [ ] **Step 4: Run the new test — it passes, and one old test now fails**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache plugin_cache
```
Expected: `rejects_unknown_key_instead_of_ignoring_it` PASSES, and
`provider::tests::registers_csv_source_with_declared_schema` now **FAILS** on `unknown field 'tier'`.
That failure is the point: it proves the dead block was being silently swallowed.

- [ ] **Step 5: Delete the dead `[tier]` block from the provider test fixture**

In `provider.rs`, inside `registers_csv_source_with_declared_schema`, remove these two lines from the
TOML literal (they sit after the `[[value_columns]]` block):

```toml
[tier]
threshold = 0.01
```

- [ ] **Step 6: Run the whole plugin_cache suite green**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache plugin_cache
```
Expected: PASS, all tests.

- [ ] **Step 7: Commit**

```bash
git add datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs \
        datafusion/bio-function-vep/src/plugin_cache/provider.rs
git commit -m "fix(plugin_cache): reject unknown manifest keys instead of silently ignoring them

A typo — or the [tier] block that the design docs still show and SourceManifest
never had — parsed happily and did nothing. The provider test fixture carried
exactly such a dead [tier] block, which is why this went unnoticed."
```

---

### Task 2: Parse `[source.vcf]` (INFO field selection)

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs`
- Test: same file

- [ ] **Step 1: Write the failing test**

Add to `mod tests` in `source_manifest.rs`:

```rust
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
```

- [ ] **Step 2: Run it and watch it fail**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache \
  plugin_cache::source_manifest::tests::parses_vcf_source_with_selected_info_fields
```
Expected: FAIL — compile error, `no field 'vcf' on SourceSpec`.

- [ ] **Step 3: Add `VcfParams` and the `SourceSpec.vcf` field**

In `source_manifest.rs`, add after `CsvParams`/`default_delim`:

```rust
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
```

and extend `SourceSpec` with one field (everything else unchanged):

```rust
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
```

Fix the one constructor that builds `SourceSpec` literally — in the same file's tests,
`single_source_table_name_has_no_part_suffix` must gain `vcf: None,`:

```rust
        let src = SourceSpec {
            part: None,
            provider: ProviderKind::Csv,
            path: "x".into(),
            csv: None,
            vcf: None,
        };
```

- [ ] **Step 4: Run it green**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache plugin_cache::source_manifest
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs
git commit -m "feat(plugin_cache): parse [source.vcf] with optional info_fields selection"
```

---

### Task 3: Wire `VcfTableProvider` into `register_sources`

`ProviderKind::Vcf` already passes manifest validation and then explodes at build time with
`NotImplemented` (`provider.rs:111-116`) — the worst possible ordering. The reader already exists and
is already a dependency; this is wiring.

**Files:**
- Modify: `datafusion/bio-function-vep/Cargo.toml` (add `datafusion-bio-format-core`)
- Modify: `datafusion/bio-function-vep/src/plugin_cache/provider.rs`
- Test: same file

- [ ] **Step 1: Write the failing test**

Add to `mod tests` in `provider.rs`:

```rust
    #[tokio::test(flavor = "multi_thread")]
    async fn registers_vcf_source_and_projects_info_fields() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("demo.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr22\t22893742\t.\tC\tG\t.\t.\tSCORE=0.9\n\
             chr22\t22893800\t.\tA\tT\t.\t.\tSCORE=0.1\n",
        )
        .unwrap();

        let toml_src = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "{}"
  [source.vcf]
  info_fields = ["SCORE"]

[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
"##,
            vcf.display()
        );

        let manifest: SourceManifest = toml::from_str(&toml_src).unwrap();
        let ctx = SessionContext::new();
        let _temps = register_sources(&ctx, &manifest).await.unwrap();

        // INFO columns are bare, case-sensitive keys -> must be backticked.
        let batches = ctx
            .sql("SELECT chrom, pos, `SCORE` FROM plugin_demo_src ORDER BY pos")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(rows, 2, "both VCF records must be visible");

        // The reader must report 1-based POS (22893742), matching the cache's 1-based start.
        let pos = batches[0]
            .column(1)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::UInt32Array>()
            .expect("pos is UInt32");
        assert_eq!(pos.value(0), 22_893_742);
    }
```

> If `pos` turns out to be a different integer width, the downcast panics with a message naming the
> actual type — fix the `downcast_ref` to that type and keep the value assertion. Do not weaken the
> assertion to "some number".

- [ ] **Step 2: Run it and watch it fail**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache \
  plugin_cache::provider::tests::registers_vcf_source_and_projects_info_fields
```
Expected: FAIL — `NotImplemented: provider Vcf not wired in prototype`.

- [ ] **Step 3: Add the `bio-format-core` dependency**

In `datafusion/bio-function-vep/Cargo.toml`, next to the existing `datafusion-bio-format-vcf` line
(line 33), add — **same git tag**, or the two crates will disagree on `ObjectStorageOptions`:

```toml
datafusion-bio-format-core = { git = "https://github.com/biodatageeks/datafusion-bio-formats.git", tag = "v1.8.8" }
```

- [ ] **Step 4: Wire the provider**

In `provider.rs`, add to the imports at the top:

```rust
use std::sync::Arc;

use datafusion_bio_format_core::object_storage::ObjectStorageOptions;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
```

and replace the `ProviderKind::Vcf | ProviderKind::Bed => { ... }` arm (currently `provider.rs:111-116`)
with these two arms. The surrounding loop already binds `let table = spec.table_name(&manifest.plugin_name);`
(a `String`), so `table.as_str()` is the table name to register under:

```rust
            ProviderKind::Vcf => {
                let info_fields = spec.vcf.as_ref().and_then(|v| v.info_fields.clone());
                // ObjectStorageOptions::default() sets compression_type = AUTO, which is what
                // lets one code path read both plain `.vcf` and BGZF `.vcf.gz`. Passing `None`
                // is NOT equivalent — the reader's own tests always pass explicit options.
                //
                // coordinate_system_zero_based = false: VCF POS is 1-based and the plugin cache
                // stores 1-based start/end, so the reader must not shift. The manifest's
                // `coordinate_system` remains the single source of truth for any shift the
                // builder applies (see plugin_cache::build::wrap_normalization).
                let vcf_table = VcfTableProvider::new(
                    spec.path.clone(),
                    info_fields,
                    None,
                    Some(ObjectStorageOptions::default()),
                    false,
                )?;
                ctx.register_table(table.as_str(), Arc::new(vcf_table))?;
            }
            ProviderKind::Bed => {
                return Err(DataFusionError::NotImplemented(format!(
                    "provider 'bed' is not wired yet (table '{table}'); \
                     no plugin needs it — wire datafusion-bio-format-bed the way \
                     ProviderKind::Vcf is wired when one does"
                )));
            }
```

Also update the module doc comment on line 1-3, which currently claims VCF/BED "are not wired in the
prototype":

```rust
//! Provider factory: register a source manifest's raw tables under their
//! `plugin_<name>_src[_<part>]` names. CSV/TSV/Parquet use builtin DataFusion
//! providers; VCF uses bio-formats' `VcfTableProvider`. BED is not wired yet
//! (no plugin needs it).
```

- [ ] **Step 5: Run it green**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache plugin_cache::provider
```
Expected: PASS — both the CSV test and the new VCF test.

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/Cargo.toml Cargo.lock \
        datafusion/bio-function-vep/src/plugin_cache/provider.rs
git commit -m "feat(plugin_cache): wire VcfTableProvider for provider = \"vcf\"

Unblocks the VCF-sourced Bucket A plugins (Mastermind, gnomADMt, EVE, Geno2MP,
SubsetVCF), which the feasibility analysis had classified as 'manifest only'
while the provider actually returned NotImplemented."
```

---

### Task 4: Make `--overwrite` mean something

`builder.rs:106` reads, literally, `let _ = self.overwrite;` — the flag is plumbed through the
builder and then dropped. A user asking for a clean rebuild silently gets an incremental one.

The one thing it should control is the **preservation of prior chrom entries**: `build_all` reuses
the previous `manifest.json`'s chroms so a filtered rebuild UPSERTs instead of dropping. With
`overwrite = true`, start from an empty chrom list.

**Files:**
- Modify: `datafusion/bio-function-vep/src/plugin_cache/builder.rs:90-106`
- Test: same file

- [ ] **Step 1: Write the failing test**

Add to `mod tests` in `builder.rs`. This is the exact mirror of the existing
`filtered_rebuild_preserves_other_chroms` test (`builder.rs:280-344`) — same `write_gz` /
`write_variation` helpers, same fixture — but asserting the opposite outcome under `--overwrite`.
Note the builder takes chrom filters as raw labels (`"1"`, `"2"`) and canonicalises them to
`chr1`/`chr2` in the manifest.

```rust
    // The mirror image of `filtered_rebuild_preserves_other_chroms`: with --overwrite,
    // a filtered rebuild must NOT preserve the previously built chroms.
    #[tokio::test(flavor = "multi_thread")]
    async fn overwrite_drops_chroms_from_a_previous_build() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("src.tsv.gz");
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
        let cache_dir = dir.path().join("cache");
        let out = dir.path().join("out");

        // Baseline: build both chroms.
        let first = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1", "2"])
            .build_all()
            .await
            .unwrap();
        assert_eq!(first.chroms.len(), 2, "baseline must contain both chroms");

        // Rebuild chr2 alone, with overwrite: chr1 must be gone from the manifest.
        let second = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .with_overwrite(true)
            .build_all()
            .await
            .unwrap();

        let chroms: Vec<&str> = second.chroms.iter().map(|c| c.chrom.as_str()).collect();
        assert_eq!(chroms, ["chr2"], "overwrite must not preserve chr1: {chroms:?}");
    }
```

- [ ] **Step 2: Run it and watch it fail**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache \
  plugin_cache::builder::tests::overwrite_drops_chroms_from_a_previous_build
```
Expected: FAIL — `overwrite must not preserve chr1: ["chr1", "chr2"]`.

- [ ] **Step 3: Honour the flag**

In `builder.rs`, in `build_all`, replace the chrom-preservation block (currently lines ~101-106,
ending in `let _ = self.overwrite;`) with:

```rust
        // Preserve chromosomes from a prior build (their shards remain on disk), so a
        // filtered/incremental build UPSERTs the rebuilt chroms rather than dropping the
        // others from the manifest that runtime discovery relies on. Two things suppress
        // that: an explicit --overwrite (the user asked for a clean build), and a schema
        // change (old shards would be misprojected under the new schema).
        let mut chroms: Vec<ChromEntry> = if self.overwrite {
            Vec::new()
        } else {
            std::fs::read_to_string(plugin_dir.join("manifest.json"))
                .ok()
                .and_then(|t| serde_json::from_str::<CacheManifest>(&t).ok())
                .filter(|old| schema_matches(old, &cache))
                .map(|m| m.chroms)
                .unwrap_or_default()
        };
```

(Delete the `let _ = self.overwrite;` line.)

- [ ] **Step 4: Run it green**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache plugin_cache::builder
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/plugin_cache/builder.rs
git commit -m "fix(plugin_cache): --overwrite was a no-op (let _ = self.overwrite)

It now drops prior chrom entries, so a clean rebuild is actually clean instead of
silently incremental."
```

---

### Task 5: `--source-path` must reach every source, not just the first

`examples/build_plugin.rs:38-42` applies the override to `manifest.sources.first_mut()`. A multi-part
manifest (`[[source]] part = "snv"` + `part = "indel"`) silently keeps the stale path for every
source but the first — and the build then reads whatever the manifest's baked-in path points at.

New semantics, explicit rather than clever:
- `--source-path <path>` — allowed **only** when the manifest has exactly one source. With several,
  it is an error telling you to use the part form.
- `--source-path <part>=<path>` — repeatable; overrides the source with that `part`.

**Files:**
- Modify: `datafusion/bio-function-vep/examples/build_plugin.rs`

- [ ] **Step 1: Add a repeatable-arg helper and the override logic**

Replace the `arg()` helper (lines 17-22) with `arg()` plus an `args_all()`:

```rust
fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

/// Every occurrence of `--key <value>` (the flag may repeat).
fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .enumerate()
        .filter(|(_, a)| a.as_str() == key)
        .filter_map(|(i, _)| args.get(i + 1).cloned())
        .collect()
}
```

and replace the `--source-path` block (lines 38-42) with:

```rust
    for spec in args_all(&args, "--source-path") {
        match spec.split_once('=') {
            // --source-path <part>=<path>
            Some((part, path)) => {
                let target = manifest
                    .sources
                    .iter_mut()
                    .find(|s| s.part.as_deref() == Some(part))
                    .ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "--source-path '{part}=...': no [[source]] with part = \"{part}\""
                        ))
                    })?;
                target.path = path.to_string();
            }
            // --source-path <path> — unambiguous only for a single-source manifest
            None => {
                if manifest.sources.len() != 1 {
                    return Err(DataFusionError::Execution(format!(
                        "--source-path <path> is ambiguous: the manifest has {} sources; \
                         use --source-path <part>=<path> (parts: {})",
                        manifest.sources.len(),
                        manifest
                            .sources
                            .iter()
                            .map(|s| s.part.as_deref().unwrap_or("<none>"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    )));
                }
                manifest.sources[0].path = spec;
            }
        }
    }
```

- [ ] **Step 2: Verify it compiles**

```bash
cargo build -p datafusion-bio-function-vep --features parquet-cache --example build_plugin
```
Expected: builds clean.

- [ ] **Step 3: Verify the ambiguity guard fires**

Point it at the AlphaMissense manifest (single source → the bare form must be accepted) and confirm
the run gets past argument parsing to the "variation shard not found" error, which proves the
override was applied:

```bash
cargo run -p datafusion-bio-function-vep --features parquet-cache --example build_plugin -- \
  --manifest /Users/wojtek/Documents/vepyr/vepyr-plugins/plugins/alphamissense/alphamissense.source.toml \
  --source-path /tmp/does-not-exist.tsv.gz \
  --variation-cache-dir /tmp/nonexistent-cache \
  --out /tmp/plugin_out
```
Expected: it prints `Building plugin 'alphamissense' from '/tmp/does-not-exist.tsv.gz'` — the path
override reached the source — and then fails with `variation shard not found`.

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/examples/build_plugin.rs
git commit -m "fix(build_plugin): --source-path reached only the first source

Multi-part manifests silently kept their baked-in paths. Bare form is now an error
when the manifest has several sources; use --source-path <part>=<path>."
```

---

### Task 6: `--chrom` must be repeatable

`arg()` returns only the first occurrence, so `--chrom 21 --chrom 22` silently builds chr21 alone —
a filtered build that quietly does less than asked.

**Files:**
- Modify: `datafusion/bio-function-vep/examples/build_plugin.rs:54-56`

- [ ] **Step 1: Use the `args_all` helper from Task 5**

Replace:

```rust
    if let Some(chrom) = arg(&args, "--chrom") {
        builder = builder.with_chrom_filter([chrom]);
    }
```

with:

```rust
    let chroms = args_all(&args, "--chrom");
    if !chroms.is_empty() {
        builder = builder.with_chrom_filter(chroms);
    }
```

(`with_chrom_filter` already takes an `IntoIterator` and treats an empty vec as "no filter", so this
is a straight swap.)

- [ ] **Step 2: Verify it compiles and both chroms are requested**

```bash
cargo run -p datafusion-bio-function-vep --features parquet-cache --example build_plugin -- \
  --manifest /Users/wojtek/Documents/vepyr/vepyr-plugins/plugins/alphamissense/alphamissense.source.toml \
  --variation-cache-dir /tmp/nonexistent-cache \
  --out /tmp/plugin_out --chrom 21 --chrom 22
```
Expected: fails with `variation shard not found: /tmp/nonexistent-cache/variation/chr21.parquet` —
i.e. it resolved the *filter*, not just the first value. (Before the fix the message is identical, so
confirm the change by also running with `--chrom 22 --chrom 21` and seeing chr22 named first.)

- [ ] **Step 3: Commit**

```bash
git add datafusion/bio-function-vep/examples/build_plugin.rs
git commit -m "fix(build_plugin): --chrom now accepts repetition (only the first was read)"
```

---

### Task 7: Update the driver's usage docs, then open the PR

The example's header comment (lines 3-9) documents the old single-`--chrom`, first-source-only
behaviour.

**Files:**
- Modify: `datafusion/bio-function-vep/examples/build_plugin.rs:1-9`

- [ ] **Step 1: Rewrite the doc header**

```rust
//! Build a plugin cache from a source manifest (all chroms, or a filtered set).
//!
//! ```text
//! cargo run -p datafusion-bio-function-vep --features parquet-cache --example build_plugin -- \
//!   --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
//!   --source-path /tmp/AlphaMissense_hg38.tsv.gz \
//!   --variation-cache-dir <cache root containing variation/> \
//!   --out /tmp/plugin_cache \
//!   [--chrom 21 --chrom 22]        # repeatable; omit to build every chrom in the cache
//! ```
//!
//! `--source-path` takes a bare path only when the manifest declares a single
//! `[[source]]`; for multi-part manifests use `--source-path <part>=<path>`, repeated.
```

- [ ] **Step 2: Run the full suite**

```bash
cargo test -p datafusion-bio-function-vep --features parquet-cache
cargo clippy -p datafusion-bio-function-vep --features parquet-cache -- -D warnings
cargo fmt --check
```
Expected: all green. Fix anything that is not before proceeding — do not open the PR on a red suite.

- [ ] **Step 3: Commit and open the PR into `dev-test`**

```bash
git add datafusion/bio-function-vep/examples/build_plugin.rs
git commit -m "docs(build_plugin): document repeatable --chrom and part-scoped --source-path"
git push -u origin feat/plugin-cache-vcf-provider
gh pr create --base master-sitekwb \
  --title "feat(plugin_cache): wire VCF provider + harden the manifest surface" \
  --body "Plan A of the plugin port factory (spec: docs/superpowers/specs/2026-07-13-vep-plugin-port-factory-design.md).

Wires \`provider = \"vcf\"\` to bio-formats' VcfTableProvider — already a dependency, so this is wiring, not implementation. Unblocks the VCF-sourced Bucket A plugins (Mastermind, gnomADMt, EVE, Geno2MP, SubsetVCF), which the feasibility analysis had filed under 'manifest only' while the provider in fact returned NotImplemented.

Plus four fixes to things that failed silently:
- unknown manifest keys were ignored (a live, dead \`[tier]\` block in our own test fixture proved it)
- \`--overwrite\` was literally \`let _ = self.overwrite;\`
- \`--source-path\` reached only the first source, so multi-part manifests kept stale paths
- \`--chrom\` read only its first occurrence

BED stays NotImplemented on purpose: no plugin needs it and it is not a dependency. The error now says so."
```

---

## What this plan does NOT do (and which plan does)

- `vepyr.build_plugin_cache` binding, the `vepyr.parity` comparator, the mini-cache slicer →
  **Plan B (toolkit)**.
- The parity harness, `parity.toml`, CI, and the three clients (AlphaMissense / REVEL /
  Mastermind-or-gnomADMt) → **Plan C (catalogue)**, which depends on A and B.
- Gene-keyed lookup, interval lookup, `match_column` fallback, wider `ValueType`, bigwig →
  out of scope by §3 and §9 of the spec.

**A VCF provider with no client is an unproven provider.** Task 3's unit test shows the reader
registers and projects; only Plan C's parity gate on a real VCF-sourced plugin proves it end-to-end.
That is exactly why the spec makes one VCF plugin part of the Wave-0 definition of done.
