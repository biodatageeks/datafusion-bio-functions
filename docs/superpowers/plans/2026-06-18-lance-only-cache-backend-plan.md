# Lance-only Cache Backend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Lance the sole variation/SIFT cache backend — permanently delete the fjall and indexed/warm-cold parquet backends (code, tests, examples, the `kv-cache` feature) — leaving a clean Lance-only crate whose public API still accepts only `backend="lance"`.

**Architecture:** This is a behavior-preserving deletion. The Lance lookup *runtime* is currently trapped inside `kv_cache::cache_exec::KvLookupExec` (shared with parquet) and borrows a few types from doomed modules. We **relocate** those borrowed pieces and the Lance engine into Lance-owned modules *while both features still build*, then **decouple** the `lance-cache` feature from `kv-cache`, then **delete** the parquet/fjall modules, then **collapse** the backend API. Correctness is proven at each risky step by a byte-identical chr1 Lance e2e diff.

**Tech Stack:** Rust 2024, DataFusion 53.0, Arrow/parquet 58, tokio, `lance`/`lance-file`/`lance-index` 7.0. Crate `datafusion-bio-function-vep`. Spec: `docs/superpowers/specs/2026-06-18-lance-only-cache-backend-design.md`.

**TDD adaptation:** No new behavior is added, so the primary gates are (a) compiles, (b) `clippy --all-targets --all-features -- -D warnings` clean, (c) `cargo test --all` green, and (d) the chr1 200k Lance `--everything` e2e output stays **byte-identical** (data rows) to the Task 0 baseline at threads {1,8}. New pure relocations are verified by the build; the engine relocation (Task 2) is verified by the e2e parity gate.

## Global Constraints

- DataFusion `=53.0.0`, arrow-ipc/`parquet` `=58.0.0`, lance `7.0.0` (unchanged; do not bump).
- Public signatures of `annotate_vep` (UDTF) and `vcf_sink::annotate_to_vcf` MUST NOT change. The `backend`/`cache_format` arguments remain; only `lance` is accepted (clear error otherwise).
- Every commit MUST pass the pre-commit hook (`cargo fmt -- --check`, `cargo clippy ... -D warnings`). The hook runs clippy workspace-wide; keep it clean.
- Stage ONLY files this work touches (`git add <paths>`), never `git add -A` — the working tree has unrelated untracked research/bench files.
- Branch: `parallel-cache-redesign`. CI runs `cargo clippy --all-targets --all-features` + `cargo test --all`; the `setup-protoc` step is already present (commit `d832cfb`).

---

## Reference: keep / delete / relocate map (from spec §2–§3)

| Module / item | Action |
|---|---|
| `src/lance_cache/*` | KEEP |
| `src/variant_lookup_exec.rs` (`ColocatedSink`, `ColocatedCacheEntry`, …) | KEEP (gate-free shared infra) |
| `src/cache_source.rs` | KEEP (lance-gated blocks only) |
| `src/annotation_store.rs` | KEEP, simplify (drop `Fjall`/`Parquet`, delete `build_store`) |
| `src/partitioned_cache.rs` | KEEP `PartitionedLanceCache`; DELETE `PartitionedParquetCache` |
| `SiftPredictionStore` trait + `SiftPredictionStoreRef` (`kv_cache/sift_store.rs`, alias `annotate_provider.rs`) | RELOCATE → `cache_common` |
| `serialize_position_entries` / `deserialize_position_predictions` (`kv_cache/sift_store.rs`) | RELOCATE → `cache_common` |
| `warm_cache/split.rs` frequency helpers (`FrequencyFields`, `PositionFrequency`, `max_global_af`, `select_warm_positions`) | RELOCATE → `cache_common` |
| `WarmChunkContext` / `ColdProbeResult` (the subset the lance branch borrows) | RELOCATE → `cache_common` |
| Lance lookup runtime in `kv_cache::cache_exec::KvLookupExec` (`new_lance`, `lance_lookup_cell`, lance `resolve_and_take`) | RELOCATE → `lance_cache/lookup_exec.rs` |
| `src/kv_cache/*` (all other) | DELETE |
| `src/warm_cache/{build,reader,cold_parquet,chunk,chrom_cache,key,lance_variation}.rs` | DELETE |
| `src/cache_builder.rs` parquet/fjall paths | DELETE (keep only `lance_cache::build` entry) |
| `kv_cache/sift_parquet_store.rs` | DELETE |
| fjall/parquet examples + `tests/kv_lookup_partitions.rs`, `tests/vcf_roundtrip_golden_fjall.rs`, ~36 parquet lib tests | DELETE |

Line numbers in the spec/exploration drift — locate symbols with `grep -n "<symbol>" <file>` before editing.

---

## Task 0: Pin the parity baseline + confirm CI protoc

**Files:** none modified (baseline capture only).

The chr1 200k Lance e2e is the correctness backstop for the whole refactor. Capture it from the **current** (pre-refactor) binary now.

- [ ] **Step 1: Build the current bench binary.**

Run: `cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1`
Expected: `Finished`.

- [ ] **Step 2: Capture baseline outputs at threads {1,8}.**

```bash
BIN=target/release/examples/bench_annotate_vcf
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
REF=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa
for t in 1 8; do
  $BIN --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/lanceonly_base_t$t.vcf \
    --backend lance --everything --reference-fasta $REF --threads $t --no-progress 2>&1 | grep -E "rows:"
done
grep -v '^#' /tmp/lanceonly_base_t1.vcf > /tmp/lanceonly_base_t1.body
grep -v '^#' /tmp/lanceonly_base_t8.vcf > /tmp/lanceonly_base_t8.body
wc -l /tmp/lanceonly_base_t1.body /tmp/lanceonly_base_t8.body
```
Expected: both runs `rows: 200000`; the two `.body` files have equal line counts. These are the frozen reference for Task 2 and Task 6. If the cache/ref/input paths differ on the machine, record the working ones here before proceeding.

- [ ] **Step 3: Confirm the protoc CI step exists.**

Run: `grep -n "setup-protoc" .github/workflows/ci.yml`
Expected: one match (commit `d832cfb`). If absent, add it (see spec §7) before continuing — `--all-features` builds `lance-index`/`prost`.

- [ ] **Step 4: No commit** (baseline capture only).

---

## Task 1: Relocate shared bits into `cache_common` (both features still build)

**Files:**
- Create: `src/cache_common.rs`
- Modify: `src/lib.rs` (add `mod cache_common;`), `src/kv_cache/sift_store.rs`, `src/warm_cache/split.rs`, `src/warm_cache/chunk.rs`, `src/warm_cache/cold_parquet.rs`, `src/annotate_provider.rs`, `src/lance_cache/build.rs` (import paths).

**Interfaces:**
- Produces: `crate::cache_common::{SiftPredictionStore, SiftPredictionStoreRef, serialize_position_entries, deserialize_position_predictions, FrequencyFields, PositionFrequency, max_global_af, select_warm_positions, WarmChunkContext, ColdProbeResult}` — same names/signatures as today, just a new home.

This is additive: move definitions, leave `pub use` re-exports at the old paths so nothing else breaks yet. Keep `cache_common` gate-free (or `#[cfg(any(feature="kv-cache", feature="lance-cache"))]` if it pulls feature-only deps).

- [ ] **Step 1: Locate the items.**

```bash
grep -n "trait SiftPredictionStore\|type SiftPredictionStoreRef\|fn serialize_position_entries\|fn deserialize_position_predictions" src/kv_cache/sift_store.rs src/annotate_provider.rs
grep -n "struct FrequencyFields\|struct PositionFrequency\|fn max_global_af\|fn select_warm_positions" src/warm_cache/split.rs
grep -n "struct WarmChunkContext" src/warm_cache/chunk.rs
grep -n "enum ColdProbeResult\|struct ColdProbeResult" src/warm_cache/cold_parquet.rs
```
Expected: each symbol found. Note exact signatures/derives to copy verbatim.

- [ ] **Step 2: Create `src/cache_common.rs`** containing the moved definitions (copy the exact bodies, derives, and doc comments from the locations in Step 1). Add `mod cache_common;` (or `pub(crate) mod cache_common;`) to `src/lib.rs` near the other module decls.

- [ ] **Step 3: Replace the originals with re-exports.** In each source file, delete the moved definition and add `pub(crate) use crate::cache_common::<Item>;` so existing references (`crate::kv_cache::sift_store::serialize_position_entries`, `crate::warm_cache::split::FrequencyFields`, etc.) still resolve.

- [ ] **Step 4: Build both feature sets.**

Run:
```bash
cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error|Finished" | head
cargo build -p datafusion-bio-function-vep --features cache-builder 2>&1 | grep -E "error|Finished" | head
```
Expected: both `Finished`.

- [ ] **Step 5: Clippy + commit.**

```bash
cargo clippy -p datafusion-bio-function-vep --features lance-cache --lib 2>&1 | grep -E "error|warning|Finished" | head
git add src/cache_common.rs src/lib.rs src/kv_cache/sift_store.rs src/warm_cache/split.rs src/warm_cache/chunk.rs src/warm_cache/cold_parquet.rs src/annotate_provider.rs src/lance_cache/build.rs
git commit -m "refactor(vep): relocate Lance-shared cache types into cache_common"
```
Expected: clippy clean; commit succeeds (pre-commit hook passes).

---

## Task 2: Relocate the Lance lookup engine → `lance_cache/lookup_exec.rs`

**Files:**
- Create: `src/lance_cache/lookup_exec.rs`
- Modify: `src/lance_cache/mod.rs` (add `pub mod lookup_exec;`), `src/lookup_provider.rs` (lance dispatch), `src/annotate_provider.rs` (if it references the lance exec type).

**Interfaces:**
- Consumes: `crate::cache_common::{WarmChunkContext, ColdProbeResult}`, `crate::variant_lookup_exec::{ColocatedSink, ColocatedCacheEntry}`, `crate::lance_cache::variation_runtime::SinglePathLanceVariationLookup`.
- Produces: `crate::lance_cache::lookup_exec::LanceLookupExec` with the same constructor inputs that `KvLookupExec::new_lance(...)` takes today and an identical `ExecutionPlan` impl (schema, `properties() -> &Arc<PlanProperties>`, `execute` producing the lance lookup stream).

This is the risky step. Create a Lance-only exec that contains ONLY the lance code path; leave `KvLookupExec` (parquet/fjall) intact for now so both still build. Verify by byte-identical e2e parity.

- [ ] **Step 1: Map the lance path inside `cache_exec.rs`.**

```bash
grep -n "fn new_lance\|lance_lookup_cell\|SinglePathLanceVariationLookup\|fn resolve_and_take\|KvMatchMode" src/kv_cache/cache_exec.rs | head -40
```
Identify: the `new_lance` constructor, the `lance_lookup_cell` field, the lance branch of `resolve_and_take` (and any helpers it calls that are lance-specific vs parquet-specific). Note the `ExecutionPlan` trait methods and the struct fields the lance path actually reads.

- [ ] **Step 2: Create `LanceLookupExec`** in `src/lance_cache/lookup_exec.rs` — a struct holding only the fields the lance path uses (input plan, schema, `output_col_positions`, `Arc<PlanProperties>`, `lance_lookup_cell: OnceCell<SinglePathLanceVariationLookup>`, colocated sink, match mode, the lance cache root/handles). Port `new_lance`'s body verbatim into `LanceLookupExec::new`, and port the lance branch of `resolve_and_take` into this struct's execute path. Reuse `cache_common::{WarmChunkContext, ColdProbeResult}` and `variant_lookup_exec::ColocatedSink` for scaffolding. Keep the SAME logic — no behavior changes.

- [ ] **Step 3: Point the lance dispatch at the new exec.**

```bash
grep -n "KvLookupExec\|new_lance\|set_lance_cache_root\|new_indexed_parquet" src/lookup_provider.rs | head
```
In `lookup_provider.rs`, change the lance dispatch (`new_lance` path) to construct `crate::lance_cache::lookup_exec::LanceLookupExec` instead of `KvLookupExec::new_lance`. Leave the indexed-parquet dispatch untouched.

- [ ] **Step 4: Build (lance feature).**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error|Finished" | head`
Expected: `Finished` (KvLookupExec's lance bits may now be dead — that's fine; they get deleted in Task 4).

- [ ] **Step 5: MANDATORY e2e parity gate.**

```bash
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1
BIN=target/release/examples/bench_annotate_vcf
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
REF=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa
for t in 1 8; do
  $BIN --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/lanceonly_t2_t$t.vcf \
    --backend lance --everything --reference-fasta $REF --threads $t --no-progress 2>&1 | grep rows:
done
diff <(grep -v '^#' /tmp/lanceonly_t2_t1.vcf) /tmp/lanceonly_base_t1.body | head
diff <(grep -v '^#' /tmp/lanceonly_t2_t8.vcf) /tmp/lanceonly_base_t8.body | head
```
Expected: both `diff`s **empty** (byte-identical to Task 0 baseline). Any mismatch → STOP, debug with superpowers:systematic-debugging; do not proceed.

- [ ] **Step 6: Commit.**

```bash
git add src/lance_cache/lookup_exec.rs src/lance_cache/mod.rs src/lookup_provider.rs src/annotate_provider.rs
git commit -m "refactor(vep): relocate Lance variation lookup engine to lance_cache::lookup_exec"
```

---

## Task 3: Decouple the `lance-cache` feature from `kv-cache`

**Files:** Modify `datafusion/bio-function-vep/Cargo.toml`, `src/annotate_provider.rs` (cfg flips), and any `src/lance_cache/*` / relocated files still importing under `kv-cache` gates.

**Interfaces:** after this task, `cargo build --no-default-features --features lance-cache` must succeed (Lance with NO kv-cache).

- [ ] **Step 1: Flip the cfg gates on Lance code.**

```bash
grep -rn 'cfg(all(feature = "kv-cache", feature = "lance-cache"))' src/ | head -40
```
Change each `#[cfg(all(feature = "kv-cache", feature = "lance-cache"))]` that guards Lance code (e.g. the Lance SIFT store impls `PositionSlicedLanceSiftStore`, `LanceBinarySiftPredictionStore`, `InMemorySiftPredictionStore`, and `load_lance_sift_prediction_store_for_chrom` in `annotate_provider.rs`) to `#[cfg(feature = "lance-cache")]`. Re-gate `cache_common` to `#[cfg(any(feature = "kv-cache", feature = "lance-cache"))]` if needed.

- [ ] **Step 2: Make `lance-cache` self-sufficient in Cargo.toml.**

Change `lance-cache = ["kv-cache", "dep:lance", "dep:lance-file", "dep:lance-index"]` to drop `"kv-cache"` and add the optional deps the Lance engine uses (formerly kv-gated): `lance-cache = ["dep:lance", "dep:lance-file", "dep:lance-index", "dep:ahash", "dep:hashbrown", "dep:nohash-hasher", "dep:rapidhash", "dep:smallvec", "dep:arrow-ipc", "dep:lz4_flex", "dep:zstd"]`. (Exact set is whatever the next step's compiler errors demand — add/remove to make it build.)

- [ ] **Step 3: Build Lance WITHOUT kv-cache.**

Run: `cargo build -p datafusion-bio-function-vep --no-default-features --features lance-cache 2>&1 | grep -E "error\[|error:|Finished" | head -40`
Expected: iterate — each `error[E0433]/unresolved`/`cfg` points at a Lance item still under a `kv-cache` gate or a missing dep. Fix by flipping the gate (Step 1) or adding the dep (Step 2). Loop until `Finished`.

- [ ] **Step 4: Commit.**

```bash
git add datafusion/bio-function-vep/Cargo.toml src/annotate_provider.rs src/cache_common.rs <other touched>
git commit -m "refactor(vep): make lance-cache feature independent of kv-cache"
```

---

## Task 4: Flip default to lance-cache; delete kv/warm/parquet modules, tests, examples

**Files:** `Cargo.toml` (features, deps, `[[example]]`/`[[test]]` blocks), `src/lib.rs` (module decls), delete `src/kv_cache/`, parquet `src/warm_cache/*`, `src/cache_builder.rs` parquet paths, `PartitionedParquetCache` in `src/partitioned_cache.rs`, parquet/fjall examples and tests.

- [ ] **Step 1: Flip default + features.** In `Cargo.toml`: `default = ["lance-cache"]`; remove the `kv-cache` feature; change `cache-builder = ["lance-cache", "dep:datafusion-bio-format-ensembl-cache"]`; remove `dep:fjall` and any parquet/fjall-only optional deps no longer referenced.

- [ ] **Step 2: Delete the modules.**

```bash
git rm -r src/kv_cache
git rm src/warm_cache/build.rs src/warm_cache/reader.rs src/warm_cache/cold_parquet.rs src/warm_cache/chunk.rs src/warm_cache/chrom_cache.rs src/warm_cache/key.rs src/warm_cache/lance_variation.rs
```
(If `warm_cache/mod.rs` now has no members, `git rm src/warm_cache/mod.rs` and remove its `mod warm_cache;`. Any type the Lance path still needs was already moved to `cache_common` in Task 1.) Update `src/lib.rs`: remove `mod kv_cache;` and `mod warm_cache;` (and their cfg gates).

- [ ] **Step 3: Strip parquet from `cache_builder.rs` and `partitioned_cache.rs`.** In `cache_builder.rs` keep only the `lance_cache::build::build_lance_entity` path; delete the warm/indexed-parquet build functions and their `kv_cache`/`warm_cache` imports. In `partitioned_cache.rs` delete `PartitionedParquetCache` and its `detect`; keep `PartitionedLanceCache`. Update `PartitionedAnnotationCache` (annotate_provider.rs) to drop the `Parquet` variant.

- [ ] **Step 4: Delete parquet/fjall examples + tests.**

```bash
grep -ln "kv_cache\|warm_cache\|indexed_parquet\|use_fjall\|backend.*parquet\|backend.*fjall\|PartitionedParquetCache" tests/*.rs examples/*.rs
```
`git rm` the fjall/parquet integration tests (`tests/kv_lookup_partitions.rs`, `tests/vcf_roundtrip_golden_fjall.rs`) and the fjall/parquet examples; remove their `[[example]]`/`[[test]]` blocks from `Cargo.toml`. For in-`src` `#[cfg(test)]` modules, delete the ~36 parquet/fjall tests (the `test_annotate_vep_*` partitioned-parquet set + `test_stateful_buffer_local_transcripts_*` that target removed paths). Keep Lance + golden-Lance tests.

- [ ] **Step 5: Build default + all-features.**

Run:
```bash
cargo build -p datafusion-bio-function-vep 2>&1 | grep -E "error|Finished" | head -40
cargo build -p datafusion-bio-function-vep --all-features 2>&1 | grep -E "error|Finished" | head -40
```
Expected: both `Finished`. Resolve dangling references to deleted items as they surface.

- [ ] **Step 6: Commit.**

```bash
git add -u && git add Cargo.toml src/lib.rs
git commit -m "refactor(vep): delete fjall + parquet cache backends; default to lance-cache"
```

---

## Task 5: Collapse the backend API

**Files:** `src/annotation_store.rs`, `src/annotate_provider.rs` (`scan`), `src/annotate_table_function.rs`, `src/vcf_sink.rs`.

- [ ] **Step 1: Collapse `AnnotationBackend`.** In `annotation_store.rs`: reduce the enum to a single `Lance` variant; `parse(s)` returns `Ok(Lance)` for `"lance"` else `Err(DataFusionError::Plan("backend must be 'lance'"))`; delete `build_store` + `AnnotationStore` trait + `FjallAnnotationStore`/`ParquetAnnotationStore` impls. Remove `let _store = build_store(...)` in `annotate_provider.rs::scan`.

- [ ] **Step 2: Strip `scan()`.** Remove `legacy_use_fjall`/`use_fjall`/`cache_format`/`use_indexed_parquet` derivation and the parquet detection branch; keep only the Lance path (`PartitionedLanceCache::detect` + the lance `indexed_variation_schema`/SIFT setup). Accept a `cache_format` option only if it equals `"lance"`, else error.

- [ ] **Step 3: Strip `ContigAnnotationConfig`.** Remove the now-unused fields `use_fjall`, `indexed_parquet_cache_root`, `indexed_variation_schema`, `kv_store` and every reference/initializer (incl. the test `minimal_contig_annotation_config`). Keep `lance_cache_root`, `sift_prediction_store`, `sift_direct`.

- [ ] **Step 4: Strip `to_options_json` / `AnnotateVcfConfig`** (`vcf_sink.rs`): remove `use_fjall` field + emission and any parquet/fjall option keys; keep lance + threads. Leave `annotate_to_vcf`/UDTF signatures unchanged.

- [ ] **Step 5: Build + clippy.**

```bash
cargo build -p datafusion-bio-function-vep --all-features 2>&1 | grep -E "error|Finished" | head
cargo clippy -p datafusion-bio-function-vep --all-features --lib 2>&1 | grep -E "error|warning|Finished" | head
```
Expected: `Finished`, clippy clean.

- [ ] **Step 6: Commit.**

```bash
git add -u
git commit -m "refactor(vep): collapse AnnotationBackend to Lance-only; drop fjall/parquet options"
```

---

## Task 6: Final gates

- [ ] **Step 1: Full clippy (exact CI command).**

Run: `cargo clippy --all-targets --all-features -- -D warnings 2>&1 | grep -E "error|warning|Finished" | head`
Expected: `Finished`, no warnings. (`protoc` must be installed locally; CI installs it.)

- [ ] **Step 2: Full test suite (exact CI command).**

Run: `cargo test --all 2>&1 | grep -E "test result|error" | tail -20`
Expected: all `test result: ok` — **0 failures**. The ~36 parquet/fjall failures are gone because that code/tests were deleted.

- [ ] **Step 3: MANDATORY e2e parity vs Task 0 baseline.**

```bash
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1
BIN=target/release/examples/bench_annotate_vcf
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
REF=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa
for t in 1 8; do
  $BIN --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/lanceonly_final_t$t.vcf \
    --backend lance --everything --reference-fasta $REF --threads $t --no-progress 2>&1 | grep rows:
done
diff <(grep -v '^#' /tmp/lanceonly_final_t1.vcf) /tmp/lanceonly_base_t1.body | head
diff <(grep -v '^#' /tmp/lanceonly_final_t8.vcf) /tmp/lanceonly_base_t8.body | head
```
Expected: both `diff`s **empty**.

- [ ] **Step 4: Confirm non-lance backend errors cleanly.**

```bash
$BIN --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/x.vcf --backend parquet --no-progress 2>&1 | grep -iE "backend must be 'lance'|error" | head
```
Expected: a clear "backend must be 'lance'" error, not a panic.

- [ ] **Step 5: Push + watch CI.**

```bash
git push
gh run list --branch parallel-cache-redesign -L 1
```
Expected: CI run starts; after it completes, `gh run view <id>` shows clippy + tests green.

---

## Self-Review

- **Spec coverage:** §1 goal → Tasks 4–6; §2 finding (Lance engine in kv_cache) → Task 2; §3 layout (cache_common, lookup_exec, deletes) → Tasks 1,2,4; §4 Cargo/features → Tasks 3,4; §5 API collapse → Task 5; §6 shared structs → Tasks 3,5; §7 testing/gates → Task 6; §9 sequencing → task order. All covered.
- **Placeholder scan:** the only deferred specifics are (a) the exact optional-dep set for `lance-cache` (Task 3 Step 2 — resolved by compiler), (b) the exact example/test filenames to delete (Task 4 Step 4 — resolved by the grep). Both are pinned by a concrete command, not left vague.
- **Type consistency:** `LanceLookupExec` (Task 2) is the name used in Task 4's deletion of `KvLookupExec`'s lance bits; `cache_common::{...}` names (Task 1) match the consumers in Tasks 2–3; `AnnotationBackend::Lance` (Task 5) matches the `parse` contract.
- **Risk backstop:** the engine relocation (Task 2) and the final state (Task 6) both gate on byte-identical chr1 Lance e2e — captured once in Task 0.
