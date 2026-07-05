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
| Key columns | **Shared with the variation cache verbatim**: `chrom` (Utf8), `start`/`end` (UInt32), `allele_string` (Utf8, `ref/alt`). Point-lookup key `(chrom, start)` via `encode_position_key`; alleles disambiguated by `allele_string`. |
| Frequency join granularity | **Allele-level** LEFT JOIN on the shared key `(chrom, start, allele_string)`; AF used only to derive tier, **not stored** in the plugin cache |
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
`select_query` normalizes to the **shared key columns** plus values
(`chrom, start, end, allele_string, values…` — §3.2) →
`plugin_cache::join::join_variation_frequency` LEFT-joins each row's
`(chrom, start, allele_string)` to the variation shard for that chrom, pulls
`max_global_af`, derives `tier` (AF ≥ 0.01 → warm/0, else / no-match → cold/1),
**discards AF** → streaming writer emits `plugin/<name>/<chrom>.parquet` physically
ordered `(tier, start)` and PageDir-indexed (§4.3) → per-plugin
`manifest.json` records source kind, key columns, value→CSQ mapping, and per-chrom
row/tier counts.

**Runtime data flow:** at annotation, for each enabled plugin (discovered by
manifest scan), an independent PageDir point-lookup — a clone of the variation
lookup path — fetches the variant's value row and injects its values into the
plugin's declared CSQ output fields; a miss emits VEP's empty/`.` value.

### 3.1 Table & view naming convention

Each plugin source registers **two** objects in the build `SessionContext`: a raw
table over the source file(s), and a view exposing only the ingest-ready columns.
The view is the single contract the cache builder reads from — all source-format
quirks stay quarantined behind it.

**1. Raw source table** — registered by the sibling crate directly over the source
file(s), columns verbatim (lowercased header fields, source quirks intact — e.g.
`pos(1-based)`, `gerp++_rs`, VCF INFO arrays):

```
plugin_<name>_src                    # single-file source
plugin_<name>_src_<part>             # multi-file source (part ∈ snv, indel, …)
```

e.g. `plugin_cadd_src_snv`, `plugin_cadd_src_indel`, `plugin_clinvar_src`,
`plugin_dbnsfp_src`.

**2. Ingest view** — a `CREATE VIEW` applying the source's normalization projection
(`select_query` / `cadd_union_query`) over the raw table(s), exposing **only** the
ingest-ready columns and nothing else. The projection maps the raw source columns
onto the shared key convention (§3.2) — `pos → start`/`end`, `ref`/`alt` →
`allele_string` — **standardizes contig names and coordinates to match the
variation cache** (§3.3), then appends the plugin's value columns:

```
plugin_<name>_ingest                 # (chrom, start, end, allele_string, <value columns…>)
```

e.g. `plugin_cadd_ingest`, `plugin_dbnsfp_ingest`.

**Rationale.**

- **`plugin_` namespace** avoids collision with `variation` / context /
  `translation_sift` tables when several sources share one `SessionContext`.
- **`_src` / `_ingest` suffix pair** makes the raw→ready relationship read at a
  glance; the `_<part>` slot handles multi-file sources (CADD SNV+indel) that
  union into the single `_ingest` view.
- **The cache builder touches only the view.** `join_variation_frequency` reads
  `plugin_<name>_ingest`; it never sees a raw column name. The view's column set is
  exactly the manifest's `key_columns` + `value_columns` — a mismatch is a
  build-time assertion.

**Follow-up on sibling PR #177.** `select_query` / `cadd_union_query` currently
hardcode `FROM source` / `FROM source_snv`. To honor this convention they must
reference `plugin_<name>_src[...]` — either templated on the registered raw-table
name, or the builder registers the raw table(s) under those names before creating
the view.

### 3.2 Key column convention (shared with the variation cache)

The plugin cache adopts the variation cache's key columns **verbatim**, so the
join, the physical sort, the PageDir position key, and allele matching are all
literally the same code path as variation — no parallel key vocabulary.

| Column | Type | Meaning |
|---|---|---|
| `chrom` | `Utf8` | Contig label (shard is per-chrom; column retained for the join key) |
| `start` | `UInt32` | 1-based genomic start. Point plugins set `start == end == pos`; interval plugins use the true span |
| `end` | `UInt32` | 1-based genomic end |
| `allele_string` | `Utf8` | `ref/alt` (e.g. `A/G`), the variation cache's allele encoding |

- **Point-lookup key** is `(chrom, start)` via `cache::key_encoding::encode_position_key`
  — the exact encoder the variation lookup uses. Alleles are disambiguated by
  `allele_string` in the scan step, exactly as variation does.
- **`allele_string` is optional per plugin.** Allele-specific plugins (CADD,
  ClinVar, dbNSFP) populate it; a purely position-level plugin (e.g. a
  conservation track) may declare only `(chrom, start, end)` in its manifest and
  match on position alone. The manifest's `key_columns` records which subset
  applies.
- The `_ingest` view (§3.1) is where the raw source columns are renamed/derived
  into this convention; downstream (`join_variation_frequency`, the writer, the
  runtime probe) only ever sees these names.

**Why `allele_string`, not separate `ref`/`alt` columns in the cache.** The
variation cache's own key is `allele_string`, so the frequency join *requires* an
`allele_string` to match against — `ref`/`alt` would have to be merged for the join
regardless. The merge is therefore done **once**, in the `_ingest` view
(`concat(ref, '/', alt) AS allele_string`, as the sibling crate's projection
already does), rather than stored split and re-merged at both join time and
runtime. `ref`/`alt` carry no downstream value either: a plugin emits its score
columns as CSQ fields, and the annotated VCF already carries the alleles. A plugin
that genuinely needs `ref`/`alt` as *output* can still declare them as
`value_columns`; as the *key*, `allele_string` is the single shared form.

### 3.3 Contig & coordinate normalization (in the `_ingest` view)

The frequency join (§4.1) matches plugin rows against the variation shard **on the
`chrom` column literally**, and the runtime probe shares the variation position-key
space — so the `_ingest` view must emit `chrom`, `start`, `end` in exactly the form
the variation cache uses, not the source's raw form. This normalization lives in
the view (SQL), so no downstream stage ever sees a source-specific contig style or
coordinate basis.

**Contig standardization.** Plugin sources arrive in mixed styles — `chr1` (UCSC /
many VCFs), `1` (Ensembl), `chrM` / `M` / `chrMT` (mitochondria). The view maps
each to the **bare Ensembl form the variation `chrom` column stores** (`1`…`22`,
`X`, `Y`, `MT`): strip any `chr` prefix, fold `M`/`chrM`/`chrMT` → `MT`, uppercase
`X`/`Y`. (Note: the *shard filename* still uses `canonical_chrom_label`
→ `chr1.parquet`, exactly as variation splits stored-column form from
filename form.)

- **Reusable component.** A single `canonical_contig(chrom) -> Utf8` scalar UDF,
  backed by the same canonicalization as `cache::key_encoding::chrom_code` /
  `manifest::canonical_chrom_label`, is registered once and used by **every**
  `_ingest` view. This guarantees plugin and variation agree on contig spelling by
  construction rather than by each source's SQL getting it right independently. It
  is the contig-side sibling of `join_variation_frequency`'s shared tiering.
- **Correctness note.** `encode_position_key` already strips `chr` at runtime, so a
  probe would tolerate a mismatched prefix — but the *build-time* join is a plain
  column equality and would silently drop every row on a `chr1` vs `1` mismatch.
  Contig standardization is therefore a **join-correctness requirement**, not
  cosmetics; the parity gate (§8) is what catches a regression here (all-null
  plugin fields).

**Coordinate standardization.** Variation `start`/`end` are **1-based** `UInt32`.
The view converts the source's basis to match: 1-based sources (VCF POS, CADD,
dbNSFP `pos(1-based)`, AlphaMissense) pass through with a `UInt32` cast; a 0-based
half-open source (BED-like) shifts `start + 1`. Point plugins set
`start == end == pos`; interval plugins carry the true 1-based span. The source's
coordinate basis is declared per source kind so the view applies the right shift.

## 4. Build pipeline & reusable join component

### 4.1 `plugin_cache::join::join_variation_frequency`

```
join_variation_frequency(
    plugin_stream: SendableRecordBatchStream,  // from plugin_<name>_ingest: (chrom,start,end,allele_string,values…)
    variation_shard: &Path,                     // variation/<chrom>.parquet
    threshold: f64,                             // WARM_AF_THRESHOLD (0.01)
) -> SendableRecordBatchStream                  // …values… + tier:Int8
```

- Registers the variation shard and the plugin stream in a `SessionContext`,
  runs a `LEFT JOIN` on the shared key `(chrom, start, allele_string)`.
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
- Rows physically written warm-run then cold-run, each `start`-monotonic, so the
  read-side PageDir binary search holds (§4.3); overwrite/skip semantics and
  empty-shard cleanup copied verbatim.
- Per-chrom build; `chrom_filter` honored for scoped/parity builds.

### 4.3 Parquet shard layout parameters

Plugin shards reuse the **exact** lookup-optimized `WriterProperties` the
variation and `translation_sift` shards use —
`parquet_cache::write::point_lookup_writer_properties(schema, &["tier", "start"])`
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
| Declared sort (`SortingColumn`) | `(tier ASC, start ASC)`, nulls last | The two-run (warm/cold) structure with `start` monotonic *within each run* that the `PageDir` binary search relies on |

**Sort-key contract.** The PageDir search key is `start` within a tier run —
**identical to variation's `(tier, start)`**. `allele_string` is **not** a PageDir
search key: rows sharing a `start` are additionally ordered by `allele_string` for
deterministic, byte-stable output, but allele disambiguation happens in the
`scan` decode step (probe filters the candidate rows by `allele_string`), exactly
as the variation lookup filters alleles after resolving position. Only `tier` and
`start` are declared as `SortingColumn`s. All sort keys are top-level primitives,
so leaf index == field index (the `point_lookup_writer_properties` assumption
holds).

The physical shape mirrors variation's two-pass write: warm run (tier 0) first,
cold run (tier 1) second, each written in `start`-ascending order by the source
`ORDER BY` — the writer never re-sorts across the run boundary, which is what
gives the PageDir its two monotonic segments.

## 5. Runtime lookup & CSQ injection

- **`plugin_cache::lookup::PluginLookup`** wraps one plugin's per-chrom
  PageDir-indexed shard; exposes
  `probe(chrom, start, allele_string) -> Option<PluginValueRow>` reusing
  `cache::key_encoding::encode_position_key`, `parquet_cache::page_dir` binary
  search, and `parquet_cache::scan` decode verbatim. Same `(tier, start)` run
  structure → same warm-first fast path; a cold-miss is a bounded page reject,
  exactly like variation.
- **Wiring point:** plugins resolve in the same per-variant colocated step where
  variation frequency is resolved today (`colocated.rs` / frequency resolution in
  `annotate_provider.rs`). Each plugin's `PluginValueRow` maps positionally onto
  its declared CSQ output fields; a miss emits VEP's empty value.
- **Allele matching:** the probe's `allele_string` reuses the same
  `alt_orig_allele_string` (`colocated.rs:308`) the frequency lookup already
  computes, so plugin allele matching is identical to frequency allele matching —
  no new allele-normalization surface.
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
  "key_columns": ["chrom", "start", "end", "allele_string"],
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
   `PluginSourceKind` variant + `register_plugin_source` arm registering
   `plugin_<name>_src[...]` (§3.1) + `select_query` whose `_ingest` view normalizes
   to the shared key columns `(chrom, start, end, allele_string, values…)` (§3.2),
   merging `ref`/`alt` into `allele_string`, applying `canonical_contig(...)` and
   the source's coordinate shift (§3.3). The checklist requires declaring the
   source's **contig style** and **coordinate basis** so the view normalizes both
   to the variation convention.
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
- **Unit:** `canonical_contig` UDF table-test (`chr1`→`1`, `chrM`/`M`/`chrMT`→`MT`,
  `X`/`Y` passthrough) plus a coordinate-shift case (0-based source → 1-based
  `start`), asserting the `_ingest` view matches variation's `chrom`/`start` form —
  the join-correctness guard from §3.3.
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
