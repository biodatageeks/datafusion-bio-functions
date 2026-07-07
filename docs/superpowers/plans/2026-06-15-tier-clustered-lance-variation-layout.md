# Tier-Clustered Single-Path Lance Variation Layout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cluster the single-path Lance variation dataset by query-hotness (warm common variants first, cold rare variants second) so point lookups dominated by common variants touch ~3x fewer miniblock chunks, recovering the ~2.5-3x `take_rows` regression introduced by the flat position-sorted layout — as a **build-only** change with **no runtime lookup changes**.

**Architecture:** During `build_lance_variation_chrom`, do a cheap pre-pass to classify each genomic `start` as warm/cold using the existing allele-frequency selector (`max_global_af >= 0.01`, ±1 bp radius), then write the dataset in **two tier-ordered passes** — warm rows first (`tier = 0`), cold rows appended after (`tier = 1`) — each pass internally `start`-sorted. A stored `tier: Int8` column records the split and a `tier` bitmap index is created for forward compatibility, but the runtime `start` BTree + `take_rows()` path is unchanged: it benefits automatically because the hot common-variant rows are now physically contiguous in a dense block instead of diluted ~24:1 among rare variants. Warm classification is computed inline from the source table (no upstream warm/cold Parquet dependency).

**Tech Stack:** Rust 2024, DataFusion 53, Arrow 58, Lance 7.0.0, existing `warm_cache::split` selector (`select_warm_positions`, `max_global_af`, `FrequencyFields`), existing `lance_cache` streaming writer and BTree/bitmap index creation.

---

## Background / Evidence

Measured on chr1 (`115_GRCh38_merged/variation.lance` flat layout vs `pylance_variation_chr1_v21_position_u32` tiered layout), 319,349 HG002 positions, all-tiers, identical projection:

| phase (402,085 matched rows) | flat (current build) | tiered (target) |
|---|---|---|
| `sidecar_take_everything` | 17.95 s / 3015 MB | 7.86 s / 1182 MB |

Cause: the 402,085 matched rows hit **20,026** distinct 4096-value miniblock chunks in the flat layout (20 rows/chunk) vs **6,192** in the tiered layout (65 rows/chunk). 97.8% of matches are warm common variants; in the tiered layout they pack into 871 dense chunks (451 rows/chunk). The flat `SELECT * FROM var WHERE chrom='1' ORDER BY chrom, start` interleaves warm and cold rows, destroying that locality. A cold-only workload (`chr1_e2e_bloom`, 0% warm) shows parity, confirming the win is specific to warm-heavy inputs.

Warm classification already exists for the Parquet tier: `warm_cache/split.rs::select_warm_positions(rows, af_threshold=0.01, position_radius=1) -> BTreeSet<i64>`, where `max_global_af` takes the max across `minor_allele_freq`, `AF`, `gnomADg`, `gnomADe`. It keys on genomic `start`, matching the Lance build's key.

Scope decisions (confirmed): **build-time clustering only** (no runtime in-memory warm serving in this plan), **warm split computed inline** from the source table.

---

## Implementation Status (2026-06-15)

Tasks 1–4 implemented; Task 5 (perf validation against a real rebuild) pending source-cache access.

- **schema.rs** — `tier` removed from `VARIATION_FORBIDDEN_COLUMNS`; `validate_variation_schema` requires `tier: Int8`; tests updated.
- **build.rs** — `collect_warm_starts` pre-pass (reuses `warm_cache::split::{select_warm_positions, max_global_af, FrequencyFields, PositionFrequency}`, constants `WARM_AF_THRESHOLD=0.01`, `WARM_POSITION_RADIUS=1`); `variation_projected_schema` appends `tier`; `transform_variation_tier_batch` derives tier + filters to one tier; `write_variation_tiered_to_lance`/`write_variation_tiered_with_ctx`/`write_variation_tier_pass` do the warm-first (`Overwrite`) then cold (`Append`) two-pass write and create `start` BTree + `tier` bitmap once at the end. Integration test `tiered_variation_write_clusters_warm_rows_first_with_indexes` proves warm-first storage order (`tier=[0,0,1,1]`), row preservation, and both indexes.
- **write.rs** — `write_record_batch_stream_to_lance_with_mode_and_version` (mode-parameterized), `create_tier_bitmap_index`, `TIER_BITMAP_INDEX_NAME`.
- **variation_runtime.rs** — `ensure_runtime_projection` explicitly excludes `tier` (build-only column never materialized into annotation output); existing test passes.

**Deviation from the original task list:** `tier` is *not* added to `VARIATION_REQUIRED_COLUMNS` (that constant enumerates columns read from the source table, which has no `tier`). Instead `tier` is derived and appended in `variation_projected_schema` + the transform, and required only on the validated output schema. Correspondingly, the runtime projection now excludes `tier` explicitly rather than relying on the forbidden list.

**Verified:** `cargo check`, `cargo clippy --lib` (clean), `cargo fmt`, all 34 `lance_cache` lib tests pass. The 30 failing `annotate_table_function` tests are pre-existing and environment-dependent (missing local `indexed_parquet` warm/cold fixture data) — confirmed identical failure on the clean tree.

**Pending (Task 5):** rebuild chr1 via `examples/build_lance_variation_chrom` (`--cache-root <source ensembl cache> --output-dir <fresh dir> --chrom chr1 --cache-source-type merged --overwrite`) and run `lance_sandbox bench` on the warm-heavy (`chr1_hg002_all_319349`) and cold-only (`chr1_e2e_bloom`) position sets; expect ~2.5–3× win on the former and parity on the latter.

---

## File Structure

- Modify: `datafusion/bio-function-vep/src/lance_cache/schema.rs`
  - Remove `tier` from `VARIATION_FORBIDDEN_COLUMNS`.
  - Add `tier` to `VARIATION_REQUIRED_COLUMNS`.
  - Require `tier: Int8` (non-null) in `validate_variation_schema`.
- Modify: `datafusion/bio-function-vep/src/lance_cache/build.rs`
  - Add a warm-position pre-pass per chromosome reusing `warm_cache::split::select_warm_positions`.
  - Add `tier: Int8` to `variation_projected_schema` and emit it in `transform_variation_batch`.
  - Write the variation dataset in two tier-ordered passes (warm `Overwrite`, cold `Append`) instead of one StreamFullDataset pass.
  - Create the `tier` bitmap index after writing (in addition to the `start` BTree).
- Modify: `datafusion/bio-function-vep/src/lance_cache/write.rs`
  - Add a `LanceIndexKind`-style or helper to create a bitmap index on `tier` (reuse `create_required_index` pattern with `BuiltinIndexType::Bitmap`).
- Verify only (no functional change expected): `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs`
  - Confirm runtime open/lookup tolerates the extra stored `tier` column and never projects it into annotation output.
- Benchmark (outside repo build artifacts): rerun `research/lance_encoding_sandbox/crates/lance_sandbox` `bench` against the rebuilt dataset.

---

## Task 1: Schema — admit a stored `tier` column

**Step 1: Move `tier` from forbidden to required**
- [ ] In `schema.rs`, delete `"tier"` from `VARIATION_FORBIDDEN_COLUMNS`.
- [ ] Add `"tier"` to `VARIATION_REQUIRED_COLUMNS` (place after `failed`/`variation_name`; order only affects column layout, not correctness).
- [ ] In `validate_variation_schema`, add `require_type(schema, "tier", &DataType::Int8)?;`.
- [ ] Update unit tests in `schema.rs`: `required_variation_columns_exclude_legacy_helpers` and `variation_projection_rejects_position_key` must still pass; add a case asserting a schema **with** `tier: Int8` validates and one with `tier` of the wrong type fails.

**Verification:** `cargo test -p datafusion-bio-function-vep lance_cache::schema` passes; `rg -n '"tier"' src/lance_cache/schema.rs` shows it in required, not forbidden.

---

## Task 2: Build — compute warm split inline

**Step 1: Add a warm-position pre-pass**
- [ ] In `build.rs`, add `async fn collect_warm_starts(ctx, table_name, source_chrom) -> Result<Arc<BTreeSet<i64>>>` that runs `SELECT start, minor_allele_freq, AF, gnomADg, gnomADe FROM <table> WHERE chrom = '<chrom>'` and feeds rows into `warm_cache::split::select_warm_positions(rows, 0.01, 1)`.
  - Build `FrequencyFields` per row from the projected columns; reuse `max_global_af` (do not reimplement the AF parsing).
  - Hardcode `af_threshold = 0.01`, `position_radius = 1` to match the Parquet tier; leave a `// TODO` to thread these through `LanceCacheBuildOptions` if tuning is needed later.
- [ ] Call this once at the top of `build_lance_variation_chrom` before the write; share the resulting `Arc<BTreeSet<i64>>` with the transform closure.

**Step 2: Emit the `tier` column in the variation transform**
- [ ] Extend `variation_projected_schema` to append `Field::new("tier", DataType::Int8, false)` (with `with_lance_field_metadata` applied like other fields).
- [ ] In `transform_variation_batch`, build an `Int8Array` where `tier = if warm_starts.contains(&(start as i64)) { 0 } else { 1 }`, reading `start` from the (already uint32-cast) start column. Append it as the last column.
- [ ] Update `transform_variation_batches` test helper and `variation_transform_drops_legacy_helpers_and_casts_positions` to expect the `tier` column.

**Verification:** A unit test feeds a batch with one high-AF and one low-AF row through `transform_variation_batch` with a warm set and asserts `tier = [0, 1]` and `tier` dtype `Int8`.

---

## Task 3: Build — write in two tier-ordered passes

**Step 1: Split the variation write into warm-first + cold-append**
- [ ] Replace the single `write_query_stream_to_lance(... StreamFullDataset ...)` call in `build_lance_variation_chrom` with two streamed passes over the same `SELECT * ... ORDER BY chrom, start` plan:
  - Pass A (warm): transform tags tier, **drops rows where `tier != 0`**, writes with `WriteMode::Overwrite`.
  - Pass B (cold): same plan re-executed, transform **drops rows where `tier != 1`**, writes with `WriteMode::Append`.
  - Each pass stays `start`-sorted because the source plan is `ORDER BY start` and filtering is stable.
- [ ] This mirrors the validated notebook approach (`WRITE_TIERS_SEPARATELY`) and bounds memory (no full-dataset buffering). Accept the cost of two source scans per chromosome; log it.
- [ ] Factor the existing `write_stream_to_lance_streaming` so it can run in append mode for pass B without recreating the dataset, or generalize `write_record_batch_stream_to_lance_with_version` to accept a `WriteMode`.
- [ ] Create the `start` BTree index only **after** both passes complete (so it covers warm+cold row ids).

**Step 2: Add the `tier` bitmap index**
- [ ] After index creation for `start`, create a `tier` bitmap index using `BuiltinIndexType::Bitmap` (extend `create_required_index`/`LanceIndexKind` or add a dedicated helper in `write.rs`).
- [ ] Mark bitmap creation non-fatal-optional behind the same build path; it is not used by the current runtime but enables a future in-memory warm tier (Layer 2) without a rebuild.

**Verification:** Build a tiny synthetic 2-chrom dataset in a `#[tokio::test]`; assert (a) row count preserved vs flat build, (b) all `tier = 0` rows have lower row ids than all `tier = 1` rows (warm-first), (c) `start_btree_idx` and `tier_bitmap_idx` both present in `load_indices()`.

---

## Task 4: Runtime — confirm no regression from the extra column

- [ ] Open a tiered-built dataset through `variation_runtime` and run an existing point-lookup unit/integration test; confirm the stored `tier` column is ignored (not projected, not required by the BTree/take path).
- [ ] Run the single-path plan's existing guard greps; ensure the **runtime** still issues no `tier = ...` filter and no `tier` projection (the column is build-only): `rg -n "tier =|project.*tier" src/lance_cache/variation_runtime.rs src/kv_cache/cache_exec.rs` shows no new production filter/projection.
- [ ] chr1 + chr4 e2e parity vs a flat-built dataset: zero mismatches (existing acceptance gate from the single-path backend plan).

---

## Task 5: Validate the performance win

- [ ] Rebuild chr1 variation via the build path into a fresh dataset dir.
- [ ] Run `lance_sandbox bench --dataset-path <rebuilt> --positions-file inputs/chr1_hg002_all_319349_positions_u32.txt --lookup-tier all` and compare `sidecar_take_everything` / `btree_direct_take` seconds and `bytes_read` against the current flat dataset.
  - Expected: ~2.5-3x fewer bytes and ~2x faster on the warm-heavy set, matching the notebook tiered numbers (~7.9 s / ~1.18 GB).
- [ ] Run the same bench with `inputs/chr1_e2e_bloom_positions_u32.txt` (cold-only): expect parity with the flat layout (no regression on rare-variant workloads).
- [ ] Record results in a short report under `research/reports/` and link from this plan.

---

## Risks & Notes

- **Build cost:** two extra source scans per chromosome (one freq pre-pass + one extra full pass for the second tier). Acceptable for an offline cache build; log durations.
- **Warm fraction drift:** if a chromosome has very few common variants, the warm block is small and the win shrinks proportionally — this is expected and self-correcting (cold-only inputs already at parity).
- **Threshold coupling:** `0.01 / ±1` must stay consistent with `warm_cache/split.rs` so the Lance tier matches any other tiered consumer; centralize the constants if both paths diverge.
- **Forward path (out of scope here):** the stored `tier` column + bitmap index let a later change revive `warm_cache/lance_variation.rs::load_warm_tier` to serve the ~95% warm hits from memory with zero `take_rows` I/O (Layer 2). This plan intentionally stops at physical clustering.
- **Do not** reintroduce `position_key` — it remains forbidden; tiering keys on `start`, which the build already casts to UInt32 and indexes.
