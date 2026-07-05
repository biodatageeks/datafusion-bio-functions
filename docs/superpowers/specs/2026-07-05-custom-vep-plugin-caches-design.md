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
| Ingestion | **Declarative source manifest** (§6.1): names a table provider (bio-formats `VcfTableProvider`, or builtin DataFusion CSV/TSV/Parquet) + provider params + input schema + ingest `SELECT` + coordinate system. The cache builder standardizes any declared source into tiered Parquet — no per-plugin Rust for the common case. The sibling `PluginSourceKind`s become shipped **preset** source manifests. |
| Tiering on variation miss | **Cold** (tier 1). Warm (tier 0) = joins variation with max global AF ≥ `WARM_AF_THRESHOLD` (0.01). Consistent with variation. |
| Key columns | **Shared with the variation cache verbatim**: `chrom` (Utf8), `start`/`end` (UInt32), `allele_string` (Utf8, `ref/alt`). Point-lookup key `(chrom, start)` via `encode_position_key`; alleles disambiguated by `allele_string`. |
| Frequency join granularity | **Allele-level** LEFT JOIN on the shared key `(chrom, start, allele_string)`; AF used only to derive tier, **not stored** in the plugin cache |
| Runtime layout | **One tiered dataset per plugin** — `plugin/<name>/<chrom>.parquet`, independent PageDir point-lookup (clone of variation lookup path) |
| Declaration | **Manifest per plugin cache** (`plugin/<name>/manifest.json`), discovered at runtime by scanning `plugin/*/manifest.json` (mirrors `CHROM_MANIFEST_FILE`) |

## 3. Component boundaries & data flow

```
           source manifest (§6.1)                     this repo
┌─────────────────────────────────┐   ┌──────────────────────────────────────────┐
│ provider + params + input schema │   │  bio-function-vep                        │
│ + ingest_sql + coordinate_system │   │  ┌── plugin_cache/ (new module) ──────┐  │
│                                  │   │  │ source.rs  : manifest + provider    │  │
│ providers resolve via factory:   │──▶│  │              factory + ingest view  │  │
│  vcf/bed → bio-formats (sibling) │   │  │ build.rs   : tiered per-chrom build │  │
│  csv/tsv/parquet → builtin DF    │   │  │ join.rs    : variation-freq join    │  │
└─────────────────────────────────┘   │  │ lookup.rs  : per-variant point-look │  │
                                        │  │ manifest.rs: cache manifest r/w     │  │
                                        │  └─────────────────────────────────────┘  │
                                        │        │ reuses cache_common::             │
                                        │        │  {max_global_af, select_warm…}   │
                                        │        │ reuses parquet_cache::{write,     │
                                        │        │  page_dir, scan} point-lookup     │
                                        └──────────────────────────────────────────┘
```

**Build data flow:** the source manifest's provider factory registers the raw
table(s) → `ingest_sql` normalizes to the **shared key columns** plus values
(`chrom, start, end, allele_string, values…` — §3.2) → the builder's shared wrapper
applies `canonical_contig` + coordinate shift (§3.3) →
`plugin_cache::join::join_variation_frequency` LEFT-joins each row's
`(chrom, start, allele_string)` to the variation shard for that chrom, pulls
`max_global_af`, derives `tier` (AF ≥ 0.01 → warm/0, else / no-match → cold/1),
**discards AF** → streaming writer emits `plugin/<name>/<chrom>.parquet` physically
ordered `(tier, start)` and PageDir-indexed (§4.3) → the **cache manifest** (§6.2)
records value→CSQ mapping and per-chrom row/tier counts.

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

**2. Ingest view** — a `CREATE VIEW` whose definition is the **source manifest's
`ingest_sql`** (§6.1) over the raw table(s), exposing **only** the ingest-ready
columns and nothing else. The `SELECT` maps the raw source columns onto the shared
key convention (§3.2) — `pos → start`/`end`, `ref`/`alt` → `allele_string` — and
appends the plugin's value columns. The cache builder then wraps this view in a
**standard normalization projection** that applies `canonical_contig(...)` and the
coordinate shift (§3.3), so contig/coordinate standardization is one shared step
rather than re-authored per source:

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

**Follow-up on sibling PR #177.** Its `select_query` / `cadd_union_query` SQL is the
raw material for the 5 **preset** source manifests' `ingest_sql` (§6.1), but it
hardcodes `FROM source` / `FROM source_snv`. When lifted into preset manifests it
must reference `plugin_<name>_src[...]` to match the registered raw-table names.

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

### 3.3 Contig & coordinate normalization (builder wrapper over the `_ingest` view)

The frequency join (§4.1) matches plugin rows against the variation shard **on the
`chrom` column literally**, and the runtime probe shares the variation position-key
space — so the ingested rows must carry `chrom`, `start`, `end` in exactly the form
the variation cache uses, not the source's raw form. The cache builder applies this
as a **single shared normalization projection wrapped over the `_ingest` view**,
driven by the source manifest's `coordinate_system` (§6.1) and the
`canonical_contig` UDF — so it is written *once* and applied identically to every
plugin, and no downstream stage ever sees a source-specific contig style or
coordinate basis. `ingest_sql` therefore stays focused on column mapping; it does
not re-implement contig/coordinate rules.

**Contig standardization.** Plugin sources arrive in mixed styles — `chr1` (UCSC /
many VCFs), `1` (Ensembl), `chrM` / `M` / `chrMT` (mitochondria). The wrapper maps
each to the **bare Ensembl form the variation `chrom` column stores** (`1`…`22`,
`X`, `Y`, `MT`): strip any `chr` prefix, fold `M`/`chrM`/`chrMT` → `MT`, uppercase
`X`/`Y`. (Note: the *shard filename* still uses `canonical_chrom_label`
→ `chr1.parquet`, exactly as variation splits stored-column form from
filename form.)

- **Reusable component.** A single `canonical_contig(chrom) -> Utf8` scalar UDF,
  backed by the same canonicalization as `cache::key_encoding::chrom_code` /
  `manifest::canonical_chrom_label`, is registered once and applied by the shared
  normalization wrapper for **every** plugin. This guarantees plugin and variation
  agree on contig spelling by construction rather than by each source's SQL getting
  it right independently. It is the contig-side sibling of
  `join_variation_frequency`'s shared tiering.
- **Correctness note.** `encode_position_key` already strips `chr` at runtime, so a
  probe would tolerate a mismatched prefix — but the *build-time* join is a plain
  column equality and would silently drop every row on a `chr1` vs `1` mismatch.
  Contig standardization is therefore a **join-correctness requirement**, not
  cosmetics; the parity gate (§8) is what catches a regression here (all-null
  plugin fields).

**Coordinate standardization.** Variation `start`/`end` are **1-based** `UInt32`.
The wrapper converts the source's basis to match, driven by the manifest's
`coordinate_system`: `1-based` sources (VCF POS, CADD, dbNSFP `pos(1-based)`,
AlphaMissense) pass through with a `UInt32` cast; a `0-based-half-open` source
(BED-like) shifts `start + 1`. Point plugins set `start == end == pos`; interval
plugins carry the true 1-based span.

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

Driven by the source manifest (§6.1), per chrom:

1. **Register sources.** For each `[[source]]`, the provider factory (§6.1)
   registers `plugin_<name>_src[_<part>]` from the declared provider + params +
   input schema.
2. **Create the ingest view.** `CREATE VIEW plugin_<name>_ingest AS <ingest_sql>`.
3. **Normalize.** Wrap the view in the shared normalization projection
   (`canonical_contig`, coordinate shift from `coordinate_system` — §3.3), filtered
   to the target chrom.
4. **Join + tier.** Feed the normalized stream through `join_variation_frequency`
   (§4.1) to derive `tier` (AF discarded).
5. **Write.** Stream to `plugin/<name>/<chrom>.parquet` via the tiered writer.
6. **Emit the cache manifest** (§6.2) with per-chrom row/tier counts.

Steps 4–6 mirror `build_parquet_variation_chrom` (`cache/build.rs`) almost exactly:

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

## 6. Manifests

Two manifests bracket the build: the **source manifest** is the declarative build
*input* (how to read and normalize a source), and the **cache manifest** is the
build *output* (what was built, for runtime discovery). Keeping them separate means
the runtime never parses provider/SQL details, and a rebuild is fully described by
the source manifest.

### 6.1 Source manifest (build input)

A per-plugin TOML declaring the provider(s), input schema, ingest `SELECT`,
coordinate system, value→CSQ mapping, and tier policy. The cache builder consumes
it end-to-end — **no per-plugin Rust for the common case**.

```toml
plugin_name       = "cadd"
coordinate_system = "1-based"          # or "0-based-half-open"; drives the start shift (§3.3)

# One or more raw sources → registered as plugin_<name>_src[_<part>] (§3.1).
[[source]]
part     = "snv"                       # optional; omit for a single-file source
provider = "csv"                       # "vcf" (bio-formats) | "csv"/"tsv" | "parquet" | "bed"
path     = "…/whole_genome_SNVs.tsv.gz"
  [source.csv]                         # provider params (CSV/TSV family)
  delimiter   = "\t"
  has_header  = false
  comment     = "#"
  compression = "gzip"
  # Explicit input schema for headerless/typed TSV (dbNSFP, CADD).
  schema = [
    { name = "chrom", type = "Utf8" },
    { name = "pos",   type = "Utf8" },
    { name = "ref",   type = "Utf8" },
    { name = "alt",   type = "Utf8" },
    { name = "rawscore", type = "Utf8" },
    { name = "phred",    type = "Utf8" },
  ]

[[source]]
part = "indel"
provider = "csv"
path = "…/InDels.tsv.gz"
  # …same csv params + schema…

# The ingest view SELECT (§3.1). Maps raw cols → shared key cols + values, in the
# source's own contig-style/coordinate-basis; the builder applies canonical_contig
# + the coordinate shift as a shared wrapper (§3.3), so this SQL stays focused on
# column mapping (INFO parsing, split_part, unions, unnest for multi-allelic, …).
ingest_sql = """
SELECT chrom, CAST(pos AS INTEGER) AS start, CAST(pos AS INTEGER) AS end,
       concat(ref, '/', alt) AS allele_string,
       CAST(rawscore AS FLOAT) AS raw_score, CAST(phred AS FLOAT) AS phred_score
FROM plugin_cadd_src_snv
UNION ALL
SELECT chrom, CAST(pos AS INTEGER) AS start, CAST(pos AS INTEGER) AS end,
       concat(ref, '/', alt) AS allele_string,
       CAST(rawscore AS FLOAT) AS raw_score, CAST(phred AS FLOAT) AS phred_score
FROM plugin_cadd_src_indel
"""

[[value_columns]]
column = "raw_score";   csq_field = "CADD_RAW";   type = "Float32"
[[value_columns]]
column = "phred_score"; csq_field = "CADD_PHRED"; type = "Float32"

[tier]
threshold = 0.01
unmatched = "cold"
```

**Provider factory.** `provider` names map to a constructor via a small factory:

| `provider` | Backed by | Params (`[source.<provider>]`) |
|---|---|---|
| `vcf` | bio-formats `VcfTableProvider` (sibling crate) | INFO fields, etc. |
| `csv` / `tsv` | builtin `ctx.register_csv` + `CsvReadOptions` | `delimiter`, `has_header`, `comment`, `compression`, `schema` |
| `parquet` | builtin `ctx.register_parquet` | — |
| `bed` | bio-formats BED provider (sibling crate) | — |

Bio-formats providers delegate to the sibling crate (which owns those deps);
builtin providers are registered directly by this repo's build code. The factory
is the single extension point — a brand-new *format* adds a provider arm here (and,
if it's a bio-formats reader, in the sibling crate), not a new plugin-specific code
path. The 5 sibling `PluginSourceKind`s ship as **preset source manifests** built
on this factory.

**Escape hatch.** Row transforms not expressible in `ingest_sql` (SQL `unnest` for
multi-allelic split covers most) fall back to a sibling-crate Rust normalizer
referenced by name — the manifest-first-with-trait-escape spirit, used only when
SQL cannot.

### 6.2 Cache manifest (build output)

One `manifest.json` per plugin under `plugin/<name>/`, **emitted by the build** and
discovered at runtime by scanning `plugin/*/manifest.json`:

```json
{
  "plugin_name": "cadd",
  "source_manifest": "cadd.source.toml",
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

- `value_columns` (copied from the source manifest) is the **single source of
  truth** for both the shard's value schema and the emitted CSQ fields — build and
  runtime read the same list, so they cannot drift.
- `tier` records the policy actually used, so a rebuild is self-describing.
- `source_manifest` back-references the input that produced the cache.

## 7. Claude Code skill: `vep-add-plugin`

A markdown workflow skill (sibling to `vep-perf-profiling`) packaging the
end-to-end process as a short decision tree plus exact file:line touch-points:

1. **Author the source manifest** (§6.1) — the primary and usually only step:
   declare provider(s) + params, input schema (for CSV/TSV), `ingest_sql` mapping
   raw columns → `(chrom, start, end, allele_string, values…)`, `coordinate_system`,
   value→CSQ mapping, and tier policy. Ship it as a preset or keep it alongside the
   plugin's data.
2. **(New *format* only) add a provider arm** to the factory (§6.1) — and, if it's
   a bio-formats reader, in `datafusion-bio-format-vep-plugin`. Not needed when the
   source is VCF / CSV / TSV / Parquet / BED (already covered).
3. **Build** per-chrom for the target chroms; skill documents the driver invocation
   (which reads the source manifest and emits the cache manifest).
4. **Parity-gate** (§8): golden VEP comparison for the plugin's CSQ fields on a
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

## 9. Scope & non-goals (YAGNI)

**Target class.** The Ensembl VEP_plugins survey shows the dominant, highest-value
class is **per-variant scoring** (CADD, REVEL, dbNSFP, SpliceAI, AlphaMissense,
PrimateAI, gnomAD-style) — tabular VCF/TSV/BED keyed by `chrom/pos/ref/alt`. This
design targets that class directly via the shared `(chrom, start, allele_string)`
key and the frequency join.

**Non-goals:**

- **Per-gene / per-transcript plugins** (LOEUF, pLI, G2P, GO, Phenotypes) are **out
  of scope**: they key by gene/transcript ID, not genomic position, so neither the
  variation frequency join nor the position point-lookup applies. They would need a
  separate gene-keyed mechanism.
- **Per-interval / overlap plugins** (Conservation, regulatory tracks) are
  structurally supported by `start`/`end`, but frequency-tiering degrades to the
  interval's start position; treat as a secondary case, not the primary target.
- No merged/bundled wide plugin shard (Section-3 option 2/3): revisit only if
  per-variant N-lookup cost is measured to matter for the hot production set.
- No storing of joined AF in the plugin cache: runtime already gets frequency from
  the variation lookup.
- No per-plugin Rust transform in this repo for the common case: normalization is
  declarative (`ingest_sql` + the shared wrapper); a sibling-crate Rust normalizer
  is the escape hatch only when SQL cannot express the transform (§6.1).
```
