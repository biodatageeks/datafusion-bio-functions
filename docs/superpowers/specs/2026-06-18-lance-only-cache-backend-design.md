# Lance-only Cache Backend — Design

**Date:** 2026-06-18
**Repo:** `datafusion-bio-functions` (`datafusion/bio-function-vep`). Downstream: `polars-bio`, `vepyr`.
**Status:** Design approved; ready for implementation plan.
**Branch:** `parallel-cache-redesign`.

---

## 1. Problem & Goal

The crate carries **three** variation/SIFT cache backends — fjall (`kv-cache`), indexed/warm-cold parquet (`kv-cache` + `cache-builder`), and Lance (`lance-cache`). They are entangled: `lance-cache` *depends on* `kv-cache`, and the default feature is `kv-cache`. The parquet/fjall surface is ~20k+ LOC (`kv_cache/` ~12k, `warm_cache/` ~5k, `cache_builder.rs` ~5k) and carries ~36 lib tests that currently fail on this branch (the cache redesign made the tiered `indexed_parquet` layout the default, which the legacy tests don't build).

**Goal:** Lance is the sole production cache backend. Permanently delete the fjall and indexed/warm-cold parquet backends — code, tests, examples, and the `kv-cache` feature — leaving a clean Lance-only crate. The public API (`annotate_vep`, `annotate_to_vcf`) keeps its `backend`/`cache_format` arguments but accepts only `lance`.

**Success criteria:**
- `cargo test --all` (default = `lance-cache`) green — no parquet/fjall failures because that code/tests are gone.
- `cargo clippy --all-targets --all-features -- -D warnings` green (CI installs `protoc` for `lance-index`/`prost`).
- The chr1 200k Lance `--everything` e2e output stays **byte-identical** to before this change (the sharded-output parity gate) — proves the relocated Lance lookup engine is behavior-preserving.
- `lance-cache` no longer depends on `kv-cache`; `kv-cache` feature removed.

## 2. Key finding (why this isn't a flat delete)

The **Lance variation-lookup runtime does not live in `lance_cache/`** — it is inside `kv_cache::cache_exec::KvLookupExec` (the same struct that runs parquet/warm-cold lookup), reached via `KvLookupExec::new_lance(...)` (`kv_cache/cache_exec.rs:271`), storing a `SinglePathLanceVariationLookup` in `lance_lookup_cell` and driven from the lance branch (`resolve_and_take`, ~`cache_exec.rs:3681`). The Lance path also borrows, from doomed modules:
- the `SiftPredictionStore` trait + `SiftPredictionStoreRef` (`kv_cache/sift_store.rs:32`; `Arc<dyn …>` aliased at `annotate_provider.rs:206`), which the Lance SIFT stores implement;
- `serialize_position_entries` / `deserialize_position_predictions` (`kv_cache/sift_store.rs`), used by `lance_cache/build.rs:1501,1504,2169`;
- the frequency helpers in `warm_cache/split.rs` (`FrequencyFields`, `PositionFrequency`, `max_global_af`, `select_warm_positions`), used by `lance_cache/build.rs:42`;
- scaffolding types `WarmChunkContext` / `ColdProbeResult` that the lance branch of `KvLookupExec` references (`warm_cache/{chunk,cold_parquet}.rs`).

`variant_lookup_exec.rs` (`ColocatedSink`, …) is already gate-free shared infra and stays. Everything in `lance_cache/*` is otherwise clean and self-contained.

So: **relocate** the borrowed pieces out of the doomed modules, **free** the Lance lookup engine from `cache_exec.rs`, then delete the rest.

## 3. Target module layout

**Keep:** `lance_cache/*` (all); `partitioned_cache.rs` (Lance half — delete `PartitionedParquetCache`); `variant_lookup_exec.rs`; `cache_source.rs`; `annotation_store.rs` (simplified).

**New `lance_cache/lookup_exec.rs`:** the Lance variation-lookup engine relocated out of `kv_cache::cache_exec` — `new_lance` + `lance_lookup_cell` + the lance `resolve_and_take` branch — with all parquet/fjall branches removed. Named for Lance ownership (replaces the misleading `KvLookupExec`).

**New gate-free `cache_common.rs`** (name finalizable in plan): the small set the Lance path borrows — `SiftPredictionStore` trait + `SiftPredictionStoreRef`, the two position-entry blob serde fns, the `warm_cache::split` frequency helpers, and the `WarmChunkContext`/`ColdProbeResult` scaffolding the lance branch needs. (These are no longer "warm"; relocating de-couples them from the deleted warm/kv modules.)

**Delete:** `kv_cache/*` entirely; `warm_cache/{build,reader,cold_parquet,chunk,chrom_cache,key,lance_variation}.rs`; the parquet/fjall code paths in `cache_builder.rs` (keep only the Lance build entry, `lance_cache::build::build_lance_entity`); `kv_cache/sift_parquet_store.rs`; `PartitionedParquetCache`; the fjall/parquet examples (`build_warm_variation_cache`, `build_variation_position_index`, `build_variation_variant_bloom_index`, `build_sift_lookup_parquet`, `bench_warm_tier_chr1`, `rewrite_cold_variation_layout`, `bench_arrow_mmap_cold_lookup`, `bench_sift_parquet_lookup_ids`, and similar — exact list pinned in the plan).

## 4. Cargo / features

- `default = ["lance-cache"]` (was `kv-cache`).
- `lance-cache` drops the `"kv-cache"` dependency and re-adds the optional deps the Lance engine actually uses (today behind `kv-cache`: `ahash`, `hashbrown`, `nohash-hasher`, `rapidhash`, `smallvec`, plus any IPC/compression crates the relocated exec needs — pinned by compile errors during implementation).
- Remove the `kv-cache` feature. `cache-builder` becomes `["lance-cache", "dep:datafusion-bio-format-ensembl-cache"]` (still builds Lance caches).
- Remove `fjall` and parquet-only deps no longer referenced.

## 5. API collapse (signatures preserved)

- `AnnotationBackend` collapses to a single `Lance` variant; `parse()` returns `Lance` for `"lance"` and errors otherwise. Delete the dead `build_store()` (`scan()` already discards it at `annotate_provider.rs:12793`).
- `scan()` (`annotate_provider.rs:12786`): remove `cache_format`/`use_fjall`/`use_indexed_parquet` derivation and the parquet detection branch; Lance is the only path. A `cache_format` option is accepted only when `"lance"`.
- `to_options_json` / `AnnotateVcfConfig` (`vcf_sink.rs`): drop `use_fjall` and the parquet/fjall option emission; keep the Lance/threads options. `backend`/`cache_format` args remain on `annotate_vep`/`annotate_to_vcf`.

## 6. Shared engine structs

- `ContigAnnotationConfig` / `AnnotationWorkerState`: drop the kv-only fields `use_fjall`, `indexed_parquet_cache_root`, `indexed_variation_schema`, `kv_store`. Keep `lance_cache_root` and `sift_prediction_store` (re-gated to `lance-cache`, holding the Lance SIFT store) and `sift_direct`.
- Flip every `#[cfg(all(feature = "kv-cache", feature = "lance-cache"))]` block (e.g. the Lance SIFT store impls at `annotate_provider.rs:3153–3505`) to `#[cfg(feature = "lance-cache")]`.

## 7. Testing & gates

- Delete the ~36 parquet/fjall lib tests (the `test_annotate_vep_*` partitioned-parquet set and `test_stateful_buffer_local_transcripts_*` fixture tests that target the removed path), plus `tests/kv_lookup_partitions.rs` and `tests/vcf_roundtrip_golden_fjall.rs`. Keep the Lance golden roundtrip and Lance/engine tests.
- Gates:
  1. `cargo build -p datafusion-bio-function-vep --features lance-cache` and `--all-features` clean.
  2. `cargo clippy --all-targets --all-features -- -D warnings` clean.
  3. `cargo test --all` green (default now `lance-cache`).
  4. chr1 200k Lance `--everything` e2e: data rows byte-identical to the pre-change output at threads {1,8} (reuses the sharded-output parity harness). **Hard gate** for the engine relocation.
- CI (`.github/workflows/ci.yml`) keeps the `setup-protoc` step (already added) — needed because `--all-features`/default now build `lance-index` → `prost`.

## 8. Risks

- **Freeing the Lance engine from `cache_exec.rs` (6,223 lines).** The lance branch is interleaved with parquet/warm-cold branches and borrows shared scaffolding. Mitigation: relocate-then-strip (move the file/lance subset to `lance_cache/lookup_exec.rs`, delete parquet branches) rather than surgical extraction; the byte-identical chr1 e2e parity is the backstop that the lance lookup still behaves identically.
- **Hidden kv references.** Some `#[cfg(feature="kv-cache")]` items the Lance path uses at runtime (esp. `sift_prediction_store`, the `SiftPredictionStore` trait). Covered by §6 relocation + flipping cfgs; the compiler enumerates the rest once `kv-cache` is gone.
- **Downstream callers.** `polars-bio`/`vepyr` must pass `backend="lance"`. Preserved arg + clear error on non-lance avoids silent breakage; a coordinated downstream check is out of scope for this crate but noted.
- **Dep set drift.** Optional deps move from `kv-cache` to `lance-cache`; missing one shows as a compile error and is added — no silent behavior change.

## 9. Sequencing

1. Relocate shared bits (`SiftPredictionStore` trait + serde fns, `split` helpers, `WarmChunkContext`/`ColdProbeResult`) → `cache_common`.
2. Relocate the Lance lookup engine → `lance_cache/lookup_exec.rs`; strip parquet/fjall branches.
3. Flip `cfg(all(kv-cache, lance-cache))` → `cfg(lance-cache)`; rewire imports.
4. Cargo: `default = ["lance-cache"]`; `lance-cache` no longer needs `kv-cache`; re-add needed optional deps; `cache-builder` → lance-based; remove `kv-cache` feature + fjall/parquet deps.
5. Delete `kv_cache/*`, the parquet `warm_cache/*`, parquet/fjall paths in `cache_builder.rs`, `PartitionedParquetCache`, fjall/parquet examples and tests.
6. Collapse `AnnotationBackend` + `scan()` + `to_options_json`.
7. Verify: build, clippy `--all-features`, `cargo test --all`, chr1 Lance e2e byte-identical parity. Commit per step.
