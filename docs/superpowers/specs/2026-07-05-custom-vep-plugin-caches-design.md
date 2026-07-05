# Custom VEP Plugin Caches — Design

**Date:** 2026-07-05
**Status:** Approved design (pre-implementation)
**Repos touched:** `datafusion-bio-functions` (this repo, primary), `datafusion-bio-formats` (source layer, PR #177)

## 1. Purpose

Let a developer add a **custom VEP plugin** (CADD, ClinVar, SpliceAI, AlphaMissense,
dbNSFP, or a brand-new source) whose annotations are emitted as VEP CSQ output
fields. Plugin data is stored in a **per-chromosome, frequency-tiered Parquet
cache** that mirrors the existing `variation/<chrom>.parquet` layout, and is
point-looked-up per variant at annotation time.

The tier of each plugin variant is derived by **joining against the variation
cache** for its chromosome to obtain that variant's allele frequency. The build
phase exposes a **reusable variation-frequency join component** so every plugin —
and any future tiered consumer — shares one frequency/tiering definition.

The deliverable is co-designed: (a) the Rust subsystem (build + join + runtime
lookup + manifest) and (b) a Claude Code skill (`vep-add-plugin`) packaging the
repeatable "add a plugin" workflow on top of it.

## 2. Locked decisions

| Decision | Choice |
|---|---|
| Deliverable | Rust subsystem **and** Claude Code skill, co-designed |
| Plugin data model | Generic / self-describing: plugin declares key columns + value columns |
| Ingestion | DataFusion table via sibling `datafusion-bio-format-vep-plugin` (`PluginSourceKind` / `register_plugin_source` / `select_query`), or an arbitrary pre-registered table/SQL |
| Tiering on variation miss | **Cold** (tier 1). Warm (tier 0) = joins variation with max global AF ≥ `WARM_AF_THRESHOLD` (0.01). Consistent with variation. |
| Frequency join granularity | **Allele-level** LEFT JOIN on `(chrom, pos, ref, alt)`; AF used only to derive tier, **not stored** in the plugin cache |
| Runtime layout | **One tiered dataset per plugin** — `plugin/<name>/<chrom>.parquet`, independent PageDir point-lookup (clone of variation lookup path) |
| Declaration | **Manifest per plugin cache** (`plugin/<name>/manifest.json`), discovered at runtime by scanning `plugin/*/manifest.json` (mirrors `CHROM_MANIFEST_FILE`) |

## 3. Component boundaries & data flow

```
                 sibling repo                          this repo
┌─────────────────────────────────┐   ┌──────────────────────────────────────────┐
│ datafusion-bio-format-vep-plugin │   │  bio-function-vep                        │
│  PluginSourceKind / register_    │   │                                          │
│  plugin_source / select_query    │   │  ┌── plugin_cache/ (new module) ──────┐  │
│  → normalized DataFusion table   │──▶│  │ build.rs   : tiered per-chrom build │  │
│    (chrom,pos,ref,alt,values…)   │   │  │ join.rs    : variation-freq join    │  │
└─────────────────────────────────┘   │  │ lookup.rs  : per-variant point-look │  │
                                        │  │ manifest.rs: plugin manifest r/w    │  │
                                        │  └─────────────────────────────────────┘  │
                                        │        │ reuses cache_common::             │
                                        │        │  {max_global_af, select_warm…}   │
                                        │        │ reuses parquet_cache::{write,     │
                                        │        │  page_dir, scan} point-lookup     │
                                        └──────────────────────────────────────────┘
```

**Build data flow:** sibling crate registers the raw source →
`select_query` normalizes to `(chrom, pos, ref, alt, values…)` →
`plugin_cache::join::join_variation_frequency` LEFT-joins each row's
`(chrom, pos, ref, alt)` to the variation shard for that chrom, pulls
`max_global_af`, derives `tier` (AF ≥ 0.01 → warm/0, else / no-match → cold/1),
**discards AF** → streaming writer emits `plugin/<name>/<chrom>.parquet` physically
ordered `(tier, pos, ref, alt)` and PageDir-indexed (§4.3) → per-plugin
`manifest.json` records source kind, key columns, value→CSQ mapping, and per-chrom
row/tier counts.

**Runtime data flow:** at annotation, for each enabled plugin (discovered by
manifest scan), an independent PageDir point-lookup — a clone of the variation
lookup path — fetches the variant's value row and injects its values into the
plugin's declared CSQ output fields; a miss emits VEP's empty/`.` value.

## 4. Build pipeline & reusable join component

### 4.1 `plugin_cache::join::join_variation_frequency`

```
join_variation_frequency(
    plugin_stream: SendableRecordBatchStream,  // (chrom,pos,ref,alt,values…)
    variation_shard: &Path,                     // variation/<chrom>.parquet
    threshold: f64,                             // WARM_AF_THRESHOLD (0.01)
) -> SendableRecordBatchStream                  // …values… + tier:Int8
```

- Registers the variation shard and the plugin stream in a `SessionContext`,
  runs a `LEFT JOIN` on `(chrom, pos, ref, alt)`.
- `tier = CASE WHEN COALESCE(max_global_af(...), 0) >= threshold THEN 0 ELSE 1 END`,
  where the AF reduction reuses **`cache_common::max_global_af`** so plugin
  tiering and variation tiering share one frequency definition.
- Output carries the plugin's value columns plus derived `tier:Int8`, and nothing
  else from variation (AF discarded).

### 4.2 `plugin_cache::build`

Mirrors `build_parquet_variation_chrom` (`cache/build.rs`) almost exactly:

- Two ordered passes: warm (0) first, cold (1) second; a chrom with no warm rows
  writes a single cold pass — same `passes` logic as the variation builder.
- A `VariationParquetShardWriter`-style streaming writer **generalized to an
  arbitrary value schema** (value columns come from the manifest, not the fixed
  variation column set).
- Rows physically written warm-run then cold-run, each `pos`-monotonic, so the
  read-side PageDir binary search holds (§4.3); overwrite/skip semantics and
  empty-shard cleanup copied verbatim.
- Per-chrom build; `chrom_filter` honored for scoped/parity builds.

### 4.3 Parquet shard layout parameters

Plugin shards reuse the **exact** lookup-optimized `WriterProperties` the
variation and `translation_sift` shards use —
`parquet_cache::write::point_lookup_writer_properties(schema, &["tier", "pos"])`
— so plugin lookups inherit the variation shard's tuned point-lookup behavior
verbatim. Concrete parameters:

| Parameter | Value | Why |
|---|---|---|
| Compression | `ZSTD` level 3 | Recovers size lost by disabling the dictionary; no-dict file is actually smaller |
| Dictionary | **disabled** | A dictionary forces loading the whole row-group dictionary to read any single row — the per-take "dictionary tax"; disabled for point lookups |
| Data page size limit | **4 KiB** (`4 * 1024`) | Small pages → fine-grained page index → cheap point-lookup resolution |
| Data page row count limit | **512 rows** | Same: bounds rows decoded per resolved page |
| Page statistics | `EnabledStatistics::Page` | Emits the `ColumnIndex` + `OffsetIndex` in the footer that the read-side `PageDir` resolves position→page against |
| Max row-group rows | **1,000,000** | Amortizes footer metadata while keeping row-group pruning useful |
| Declared sort (`SortingColumn`) | `(tier ASC, pos ASC)`, nulls last | The two-run (warm/cold) structure with `pos` monotonic *within each run* that the `PageDir` binary search relies on |

**Sort-key contract.** The PageDir search key is `pos` within a tier run —
identical to variation's `(tier, start)`. `ref`/`alt` are **not** PageDir search
keys: rows sharing a `pos` are additionally ordered by `(ref, alt)` for
deterministic, byte-stable output, but allele disambiguation happens in the
`scan` decode step (probe filters the candidate rows by `ref`/`alt`), exactly as
the variation lookup filters alleles after resolving position. Only `tier` and
`pos` are declared as `SortingColumn`s. All sort keys are top-level primitives, so
leaf index == field index (the `point_lookup_writer_properties` assumption holds).

The physical shape mirrors variation's two-pass write: warm run (tier 0) first,
cold run (tier 1) second, each written in `pos`-ascending order by the source
`ORDER BY` — the writer never re-sorts across the run boundary, which is what
gives the PageDir its two monotonic segments.

## 5. Runtime lookup & CSQ injection

- **`plugin_cache::lookup::PluginLookup`** wraps one plugin's per-chrom
  PageDir-indexed shard; exposes `probe(chrom, pos, ref, alt) -> Option<PluginValueRow>`
  reusing `parquet_cache::page_dir` binary search and `parquet_cache::scan` decode
  verbatim. Same `(tier, pos)` run structure → same warm-first fast path; a
  cold-miss is a bounded page reject, exactly like variation.
- **Wiring point:** plugins resolve in the same per-variant colocated step where
  variation frequency is resolved today (`colocated.rs` / frequency resolution in
  `annotate_provider.rs`). Each plugin's `PluginValueRow` maps positionally onto
  its declared CSQ output fields; a miss emits VEP's empty value.
- **Allele matching:** the probe key reuses `alt_orig_allele_string`
  (`colocated.rs:308`), so plugin allele matching is identical to frequency allele
  matching — no new allele-normalization surface.
- **Output schema:** the plugin's CSQ field names are appended to the CSQ format
  list. One-dataset-per-plugin means enabling/disabling a plugin adds/removes
  exactly its columns, with no effect on variation or other plugins.

### 5.1 Concurrency constraint (from prior findings)

Per the aux-pool-oversubscription lesson, `PluginLookup` probes **synchronously on
the annotating thread** — no internal Rayon/async pool — sharing the engine's
thread budget exactly as the variation point-lookup does. N enabled plugins add N
independent, cheap-on-cold-miss probes per variant; they must not spawn their own
pools.

## 6. Plugin manifest

One `manifest.json` per plugin under `plugin/<name>/`, discovered at runtime by
scanning `plugin/*/manifest.json`:

```json
{
  "plugin_name": "cadd",
  "source_kind": "cadd",
  "key_columns": ["chrom", "pos", "ref", "alt"],
  "value_columns": [
    { "column": "raw_score",   "csq_field": "CADD_RAW",   "type": "Float32" },
    { "column": "phred_score", "csq_field": "CADD_PHRED", "type": "Float32" }
  ],
  "tier": { "threshold": 0.01, "unmatched": "cold" },
  "chroms": [ { "chrom": "1", "file": "1.parquet", "rows": 1234, "warm": 12, "cold": 1222 } ],
  "cache_source_version": "…"
}
```

- `source_kind` is a `PluginSourceKind` in the sibling crate, or `"table"` for an
  arbitrary pre-registered/SQL source.
- `value_columns` is the **single source of truth** for both the shard's value
  schema and the emitted CSQ fields — build and runtime read the same list, so
  they cannot drift.
- `tier` records the policy actually used, so a rebuild with a different threshold
  is self-describing.

## 7. Claude Code skill: `vep-add-plugin`

A markdown workflow skill (sibling to `vep-perf-profiling`) packaging the
end-to-end process as a short decision tree plus exact file:line touch-points:

1. **Decide the source.** One of the 5 known `PluginSourceKind`s → skip to step 3;
   a new format → step 2.
2. **(New source only) add a source kind** in `datafusion-bio-format-vep-plugin`:
   `PluginSourceKind` variant + `register_plugin_source` arm + `select_query`
   normalizing to `(chrom, pos, ref, alt, values…)`. Skill provides the checklist
   and SQL-projection pattern.
3. **Author the manifest** (§6) — value→CSQ mapping and tier policy.
4. **Build** per-chrom for the target chroms; skill documents the driver
   invocation.
5. **Parity-gate** (§8): golden VEP comparison for the plugin's CSQ fields on a
   test chrom; confirm 100% before wiring into production.

The skill stays a thin wrapper over the real subsystem — it does not embed logic
that belongs in Rust.

## 8. Testing & parity strategy

- **Unit:** `join_variation_frequency` tiering table-test (warm / cold / no-match),
  a clone of the existing `transform_variation_tier_batch` tests.
- **Round-trip:** build a tiny synthetic plugin shard; probe hits and misses via
  `PluginLookup`; assert value + tier.
- **Golden parity (decisive gate):** chr-scoped run comparing our plugin CSQ
  fields against real Ensembl VEP output (same harness as the
  `fast_chr1_22_merged` reports). CADD and ClinVar are good first targets — they
  have unambiguous golden fields.
- **Byte-identity across workers:** w1 vs w4 body-identical, per the existing
  parallel-annotation invariants.

## 9. Non-goals (YAGNI)

- No merged/bundled wide plugin shard (Section-3 option 2/3): revisit only if
  per-variant N-lookup cost is measured to matter for the hot production set.
- No storing of joined AF in the plugin cache: runtime already gets frequency from
  the variation lookup.
- No plugin-defined custom transform code in this repo: source-format quirks live
  in the sibling crate's `select_query`; this repo treats the normalized table
  generically.
```
