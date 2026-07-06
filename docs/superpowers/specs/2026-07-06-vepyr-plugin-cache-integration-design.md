# vepyr Plugin-Cache Integration — Design

**Date:** 2026-07-06
**Repos:** `vepyr` (integration target), `datafusion-bio-functions` (upstream, this repo), `vepyr-plugins` (public manifest repo)
**Predecessor:** `docs/superpowers/plans/2026-07-06-plugin-cache-handoff.md` (the upstream plugin-cache subsystem, PR #190) and its design spec `2026-07-05-custom-vep-plugin-caches-design.md`.

## 1. Goal

Surface the upstream custom-VEP-plugin-cache subsystem (proven with AlphaMissense, exact chr22 golden parity) through **vepyr** across all three of its workflows:

1. **Cache builder** — a `vepyr.build_plugin_cache(...)` Python API that builds a per-chromosome plugin cache.
2. **Annotation** — a `plugin_cache_root` knob on `vepyr.annotate(...)` that enables plugin CSQ fields.
3. **E2E validation** — a full chr1–22 plugin-parity profile in the existing `run_annotation_fast*.py` harness, gated against Ensembl VEP `--plugin AlphaMissense` golden output.

Along the way, three upstream refinements (all justified below) make the subsystem cleaner and fully declarative:
- **Tiering is sourced from the ensembl-vep variation cache**, not a per-plugin AF expression (removes `af_max_sql`).
- **The manifest tier/AF config is removed** (`TierPolicy` gone).
- **The runtime match discriminator is declared as a template** over a generic engine-attribute namespace, replacing the closed `engine_attr` enum — so a new plugin needs **no engine code change**.

## 2. Non-goals

- Full DataFusion-SQL discriminator expressions (the template covers every real discriminator; SQL-on-the-hot-path is explicitly deferred).
- Per-gene / transcript-ID-only plugins and allele-agnostic per-position tracks (out of scope for this cache shape, per the predecessor spec §9).
- Automatic download of plugin *source data* (e.g. `AlphaMissense_hg38.tsv.gz`) — the caller supplies `source_path`.

## 3. Dependency & sequencing

vepyr pins `datafusion-bio-function-vep` by git **tag** (currently `v0.13.1`, which predates all plugin-cache work). The plugin-cache code plus the refinements in this spec must be in a buildable pin before vepyr can consume them.

**Decision:** during development, vepyr's `Cargo.toml` uses a **local path dep** to `/Users/mwiewior/research/git/datafusion-bio-functions` for fast iteration; `uv sync --reinstall-package vepyr` rebuilds `_core`. The plugin-cache module ships under the default `parquet-cache` feature, so vepyr's existing `cache-builder` feature already pulls it in — no feature change needed.

**Hard rule:** the local path dep is converted to a real git **tag/rev before any vepyr commit**. Upstream release task: land the changes in §5–§7 on `feat/plugin-cache-alphamissense` (PR #190), merge to `master`, cut a new tag (e.g. `v0.14.0`); vepyr then pins that tag.

## 4. Architecture overview

```
vepyr (Python)                         vepyr (Rust _core)          datafusion-bio-functions (upstream)
──────────────────────────────────    ────────────────────────    ─────────────────────────────────────
build_plugin_cache(plugin, version,    _core.build_plugin_cache    PluginCacheBuilder::new(..)
  source_path, cache_dir, ...)   ─────▶  (thin wrapper)      ─────▶   .with_chrom_filter(..).build_all()
  └ resolve manifest from                                              └ per-chrom loop + tier-from-variation
    vepyr-plugins @ version                                              join + manifest.json accumulation

annotate(.., plugin_cache_root=..) ──▶ opts dict → options_json ──▶ AnnotateVcfConfig.plugin_cache_root
                                                                     (workers path)  / annotate_vep UDTF
                                                                     (streaming path, already parses it)
```

## 5. Cache builder

### 5.1 Upstream: `PluginCacheBuilder` (owns the chrom loop — option A)

Mirrors `CacheBuilder` (the ensembl variation-cache builder), whose whole per-chrom loop lives upstream and whose vepyr binding is a one-call wrapper. New type in `datafusion/bio-function-vep/src/plugin_cache/`:

```rust
PluginCacheBuilder::new(manifest: &SourceManifest, manifest_file: &str,
                        variation_cache_dir: &Path, out: &Path)
    .with_chrom_filter([..])          // default: all chroms present under <cache_dir>/variation/
    .with_overwrite(bool)
    .build_all().await -> Result<CacheManifest>
```

`build_all()`:
- iterates the chrom filter (or discovers chroms from the variation cache layout),
- calls the existing `build_plugin_chrom` per chrom (deriving `<cache_dir>/variation/chr{N}.parquet` as the tier source),
- accumulates each `ChromEntry` into a single `CacheManifest`,
- writes `plugin/<name>/manifest.json` once.

The existing single-chrom `build_plugin` example is retained for debugging; the `build_plugin_chrom` **signature loses `af_max_sql`** (see §6).

### 5.2 vepyr Rust: `_core.build_plugin_cache` (thin)

New `#[pyfunction]` in `vepyr/src/lib.rs`, structured like `build_cache`:

```
_core.build_plugin_cache(manifest_path, source_path, variation_cache_dir,
                         plugin_cache_root, chroms=None, overwrite=False, on_progress=None)
  → PluginCacheBuilder::new(SourceManifest::load(manifest_path), .., variation_cache_dir, plugin_cache_root)
      .with_chrom_filter(chroms).with_overwrite(overwrite).build_all()
```

It contributes **no loop logic** — symmetric with `build_cache`/`CacheBuilder`. `source_path` overrides the manifest's `[[source]].path` (as the `build_plugin` example does today).

### 5.3 vepyr Python: `vepyr.build_plugin_cache(...)`

```python
vepyr.build_plugin_cache(
    plugin: str,                 # dir name in vepyr-plugins (e.g. "alphamissense")
    version: str,                # git tag of vepyr-plugins to read THIS plugin's manifest from
    source_path: str,            # actual source DATA file (e.g. AlphaMissense_hg38.tsv.gz)
    cache_dir: str,              # ensembl-vep variation cache root (supplies tier)
    plugin_cache_root: str,      # output root -> plugin/<plugin>/chr{N}.parquet + manifest.json
    chroms: list[str] | None = None,
    plugins_repo: str | None = None,   # optional local clone of vepyr-plugins for OFFLINE builds
    on_progress=None,
    overwrite: bool = False,
) -> dict                        # per-chrom counts (rows/warm/cold)
```

**Per-plugin versioning.** One call builds one plugin at one version. Mixing versions (AlphaMissense `0.1.1` + ClinVar `0.2.0`) = call `build_plugin_cache` once per plugin into the **same `plugin_cache_root`**; each writes an independent `plugin/<name>/manifest.json`.

**Manifest resolution (Python-side, uniform online/offline).** Python owns all git/network logic (the Rust layer only ever sees a resolved local `manifest_path`):
- if `plugins_repo` is given → reuse that local clone (**full offline**, no network);
- else → clone the **public** `vepyr-plugins` repo into a local cache dir (no auth);
- either way, check out `version` for that plugin using a **detached checkout / `git worktree`** (never disturbs the caller's working branch; lets different plugins sit at different tags in one clone), then read `plugins/<plugin>/<plugin>.source.toml`.

## 6. Tiering sourced from the variation cache

The plugin cache remains **frequency-tiered** (warm/cold `tier` column → disk locality for point-lookups), but the tier is **inherited from the ensembl-vep variation cache** rather than recomputed from a per-plugin AF expression.

- The variation shard already carries a required, validated `tier` Int8 column (0=warm/common, 1=cold/rare) — `cache/schema.rs:39,122`.
- The plugin build LEFT JOINs plugin rows onto the variation shard on `(chrom, start, allele_string)` and selects **`COALESCE(v.tier, 1) AS tier`** (no variation match → cold).
- **Removed:** `af_max_sql` and the AF threshold — from the vepyr API, `PluginCacheBuilder`, `build_plugin_chrom`, and `join.rs`'s `tier_sql`/`tiered_stream`.
- **Retained:** `cache_dir` (the variation cache) — it is now the single source of truth for tier, shared with variation tiering.

Rationale: the AF threshold was a copy of the variation-cache tiering policy; sourcing tier directly from the variation record makes plugin and variation tiering consistent by construction and removes a redundant, per-plugin knob. (The runtime lookup never reads `tier` — it is purely a build-time page-ordering input — so this change is lookup-invariant.)

## 7. Manifest schema change + declarative discriminator

### 7.1 Remove tier/AF from the manifest schema

- `source_manifest.rs`: delete `TierPolicy` (its `threshold` field + `default_threshold`) and the `SourceManifest.tier` field.
- AM manifest (`vepyr-plugins/plugins/alphamissense/alphamissense.source.toml`): delete the `[tier]` block.

### 7.2 Declarative match discriminator (template over an engine-attribute namespace)

**Problem being fixed:** today the manifest's `[[match_column]].engine_attr` is a value from a **closed enum** (`EngineAttrs`, currently only `amino_acid_change`), so every new plugin discriminator requires an engine code change — defeating the "declarative, no per-plugin Rust" goal.

**Key fact:** the discriminator value is *engine-derived*, not a cache column. `amino_acid_change` (`"W320R"`) does not exist in the ensembl-vep cache; the engine computes it per transcript consequence from `Amino_acids` (`"W/R"`) + `Protein_position` (`"320"`) at the probe site (`annotate_provider.rs:6212`). So the declarative expression operates over the **engine's per-consequence attributes** (the CSQ values the engine already assembles), not raw cache columns.

**Design:**
1. Replace the closed `EngineAttrs { amino_acid_change }` struct with a **generic attribute namespace** — a name→`Option<String>` lookup populated at the probe site from values the engine already has local: `ref_aa`, `alt_aa` (from splitting `Amino_acids` `"W/R"`), `Amino_acids`, `Protein_position`, `Feature` (transcript id), `Gene`, `Consequence`, `ref`, `alt`, … (extend as cheaply-available fields are needed). Cheap: the values are already in scope in the output loop.
2. `[[match_column]]` declares a **template** string with `{name}` placeholders instead of `engine_attr`:
   ```toml
   [[match_column]]
   column   = "protein_variant"
   template = "{ref_aa}{Protein_position}{alt_aa}"
   ```
3. The template is parsed **once per plugin** (from the manifest) and evaluated per consequence via plain Rust string interpolation — **no DataFusion on the hot path**.
4. **Gating rule (preserves missense-only behaviour):** if **any** referenced placeholder resolves to `None`/empty, the whole discriminator is `None` → probe miss → empty output. A non-missense consequence has no `Amino_acids`/`Protein_position` → discriminator `None` → empty — exactly the existing gate.
5. `csq::amino_acid_change` special-case is removed; the `"W/R"` → `ref_aa`/`alt_aa` split happens when populating the namespace.

A new plugin needs engine changes **only** if it references an attribute the namespace doesn't yet expose; the common discriminators (aa-change, transcript id, gene, position) are already computed for CSQ, so in practice adding a plugin is manifest-only.

### 7.3 Fixed positional key (unchanged)

The full runtime compound key is `(start, allele_string [, discriminator…])`. `(start, allele_string)` are **fixed key columns of the cache format** (produced by `ingest_sql`, not declared as lookup config); only the discriminator part is declared, via `[[match_column]]`. No `[[match_column]]` → per-variant lookup on `(start, allele_string)`; with it → per-(variant×transcript) lookup.

### 7.4 vepyr-plugins repush (sequenced, confirmed)

The AM-manifest edits in §7.1 + §7.2 (`engine_attr = "amino_acid_change"` → `template = "{ref_aa}{Protein_position}{alt_aa}"`; drop `[tier]`) are **coupled** to the upstream code change: a manifest without `[tier]`/with `template` won't parse against current released code, and the code won't compile until `join.rs`/`build.rs`/`source_manifest.rs` change together. Therefore:
- the AM-manifest edit + a **new tagged** `vepyr-plugins` release are an explicit implementation task, sequenced **with** the upstream schema/build change;
- the actual `git push` to the `vepyr-plugins` remote is confirmed with the user immediately before pushing (outward, versioned release). No mismatched manifest is pushed mid-implementation.

## 8. Annotation

`plugin_cache_root` rides the existing `options_json` channel — **no PyO3 signature change**:
- `src/vepyr/__init__.py`: add `plugin_cache_root: str | None = None` kwarg to `annotate()`; insert into the `opts` dict (~lines 654–726) so it lands in `options_json`.
- `src/annotate.rs::annotate_to_vcf_file`: read `opts.get("plugin_cache_root")` → set `AnnotateVcfConfig.plugin_cache_root`.
- Streaming/LazyFrame path needs nothing new — the upstream `annotate_vep` UDTF already parses `plugin_cache_root` from `options_json` (`annotate_table_function.rs:105`).
- Update `_core.pyi` stubs and `docs/plugins.md`. `None` = disabled = byte-identical to a no-plugin run.

The plugin CSQ fields (header + per-transcript body) are emitted by the upstream engine purely from the manifests discovered under `plugin_cache_root`; Python only passes the path.

## 9. E2E validation (full chr1–22)

- **Golden generation:** extend `e2e-testing/vep-docker.md` with the `--plugin AlphaMissense,<data>` invocation (+ mount the AlphaMissense tabix) for chr1–22, producing golden VCFs whose CSQ carries `am_class`/`am_pathogenicity`.
- **Plugin cache:** full-genome AlphaMissense build via `vepyr.build_plugin_cache` (§5) over chr1–22.
- **Profile:** add an entry to `_CACHE_PROFILES` in `e2e-testing/scripts/run_annotation_fast.py` pairing the merged `cache_dir`, the AM-enabled golden, and `annotate_kwargs` carrying `plugin_cache_root`.
- **Comparison:** `compare_vcfs()` already diffs **all** CSQ fields generically and tracks populated-vs-empty, so the plugin fields — including the missense gating (empty on non-missense lines) — are validated automatically once golden + header carry them. Expected: no compare-logic change; verify gating counts surface.
- **All-chrom:** `run_annotation_fast_all.py --profile <plugin profile>` aggregates chr1–22 into the timestamped Markdown summary, covering the plugin fields per chrom.
- **Success gate:** 0 mismatches and 0 over-emit on populated plugin CSQ lines across chr1–22 (the chr22 result — 1,912/1,912 — extended genome-wide).

## 10. Prerequisites (operational)

- Full `AlphaMissense_hg38.tsv.gz`.
- An existing chr1–22 ensembl-vep **merged** variation cache (tier source + annotation cache).
- Docker Ensembl VEP r115 + AlphaMissense plugin data for golden generation.

## 11. Documentation deliverable

Extend vepyr's `docs/plugins.md` (the user-facing "add a plugin" reference) with two sections. This ships as part of the integration, not a follow-up.

### 11.1 Manifest structure reference

Document the `<plugin>.source.toml` layout as an authoritative reference, section by section:

- **Top-level scalars** — `plugin_name`, `coordinate_system` (`"1-based"` / `"0-based"`), `ingest_sql`. Note the TOML ordering rule: these MUST precede any `[[table]]` header or TOML absorbs them.
- **`ingest_sql`** — must project the fixed key columns `chrom`, `start`, `end`, `allele_string` (`ref/alt`), plus any discriminator column(s) and the value column(s), from the raw source view `plugin_<name>_src`.
- **`[[source]]`** — `provider` (`tsv`/`csv`/`parquet`), `path` (overridden at build time by `source_path`), and the `[source.csv]` block (`delimiter`, `has_header`, `comment`, `compression`, `schema` = ordered `{name,type}` list).
- **`[[match_column]]`** (optional, 0+) — `column` (the stored discriminator column, build-time) + `template` (the runtime expression over the engine-attribute namespace, §11.2). Omit entirely for per-variant plugins.
- **`[[value_columns]]`** (1+) — `column`, `csq_field` (output field name), `type` (`Utf8`/`Float32`); **declaration order = CSQ output order**.
- Note the removed sections: there is **no** `[tier]` block — tiering is inherited from the variation cache (§6).

Include the full AlphaMissense manifest as a worked example, and a minimal per-variant example (no `[[match_column]]`).

### 11.2 Engine-attribute namespace table

Document the names a `[[match_column]].template` may reference (the per-consequence values the engine exposes at the probe site). Each is `Option`: **if any placeholder a template references is absent, the discriminator is `None` → probe miss → empty output** (this is how missense-only gating works).

| Attribute | Description |
|---|---|
| `Consequence` | Consequence type(s) for the transcript (e.g. `missense_variant`). |
| `Gene` | Ensembl gene stable ID. |
| `Feature_type` | Feature type (e.g. `Transcript`). |
| `Feature` | Transcript stable ID — the transcript-id discriminator (e.g. dbNSFP). |
| `BIOTYPE` | Transcript biotype (e.g. `protein_coding`). |
| `HGVSc` | HGVS coding-sequence notation. |
| `HGVSp` | HGVS protein notation. |
| `cDNA_position` | Position in cDNA. |
| `CDS_position` | Position in the CDS. |
| `Protein_position` | 1-based amino-acid position. |
| `Amino_acids` | Reference/alternate amino acids as `ref/alt` (e.g. `W/R`); single value when unchanged. |
| `Codons` | Reference/alternate codons. |
| `ref_aa` | Reference amino acid (left of `/` in `Amino_acids`). |
| `alt_aa` | Alternate amino acid (right of `/` in `Amino_acids`). |
| `ref` | VCF reference allele. |
| `alt` | VCF alternate allele. |

Worked examples: AlphaMissense (aa-change) `template = "{ref_aa}{Protein_position}{alt_aa}"` → `W320R`; a transcript-keyed plugin `template = "{Feature}"`. State the rule for extending the namespace: a new attribute is added upstream only when a plugin needs a value not already listed here (the common discriminators are all present, so most plugins are manifest-only).

## 12. Affected files (map)

**datafusion-bio-functions (upstream):**
- `plugin_cache/source_manifest.rs` — remove `TierPolicy`/`tier`; `MatchColumn`: `engine_attr` → `template`.
- `plugin_cache/join.rs` — `tier_sql`/`tiered_stream` inherit `COALESCE(v.tier,1)`; drop `af_max_sql`/threshold.
- `plugin_cache/build.rs` — `build_plugin_chrom` drops `af_max_sql`.
- `plugin_cache/mod.rs` (new `builder.rs`) — `PluginCacheBuilder`.
- `plugin_cache/registry.rs` — `EngineAttrs` → generic attribute namespace; template eval.
- `plugin_cache/csq.rs` — remove `amino_acid_change` special-case; template helper.
- `annotate_provider.rs:~6212` — populate the attribute namespace from per-consequence locals.
- `examples/build_plugin.rs` — drop `--af-max-sql`.

**vepyr:**
- `Cargo.toml` — dep pin (path → tag).
- `src/lib.rs` — `_core.build_plugin_cache` pyfunction + module registration.
- `src/vepyr/__init__.py` — `build_plugin_cache()` + manifest resolution; `plugin_cache_root` kwarg on `annotate()`.
- `src/vepyr/_core.pyi` — stubs.
- `e2e-testing/vep-docker.md`, `e2e-testing/scripts/run_annotation_fast.py` (`_CACHE_PROFILES`).
- `docs/plugins.md` — `plugin_cache_root` on `annotate()`, `build_plugin_cache()` usage, **manifest structure reference + engine-attribute namespace table (§11)**.

**vepyr-plugins:**
- `plugins/alphamissense/alphamissense.source.toml` — drop `[tier]`; `engine_attr` → `template`; new tag.
