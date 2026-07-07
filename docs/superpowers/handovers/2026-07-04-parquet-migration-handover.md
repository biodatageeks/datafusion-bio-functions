# Parquet Cache Backend — Session Handover (2026-07-04)

Handover for a fresh session to continue the Phase 2 Parquet-cache migration. **Tasks 1–11 are DONE, committed, and validated end-to-end. Only Task 12 (delete Lance / drop `lance*` deps) remains.**

---

## 0. TL;DR / current state

- The **entire Ensembl cache is migrated to Parquet** (variation + translation_sift point-lookups; transcript/exon/regulatory/motif/translation_core scans) and **validated byte-identical** to the Lance backend end-to-end on chr1/chr2/chr4, and 100% vs Ensembl VEP 115.
- Tree is **GREEN at commit `7ea6d86`** on branch `worktree-parquet-cache-backend` in the worktree. `cargo build -p datafusion-bio-function-vep --features lance-cache,cache-builder` is clean.
- Runtime is a **hybrid**: variation + sift + context all read from Parquet when the shard exists, else fall back to Lance. `cache_format="parquet"` activates it.
- **Task 12** (make Parquet the ONLY backend, delete all Lance code, drop `lance`/`lance-*` deps) is the last step — a ~2000-line deletion across ~12 files + ~40 tests. A 2026-07-04 attempt was **reverted** (ran out of context after 1 edit). Full deletion map is in memory `parquet-task12-lance-deletion-map.md` and summarized in §7 below.

---

## 1. Locations

| Thing | Path |
|---|---|
| **Worktree** (do all work here) | `/Users/mwiewior/research/git/datafusion-bio-functions/.claude/worktrees/parquet-cache-backend` |
| Branch | `worktree-parquet-cache-backend` (off master `a6e19ad`) |
| Crate | `datafusion/bio-function-vep` |
| Parquet code | `datafusion/bio-function-vep/src/parquet_cache/` |
| Phase 2 plan | `docs/superpowers/plans/2026-07-02-parquet-cache-backend-phase2.md` **(in MAIN checkout, not the worktree)** |
| **vepyr** (PyO3 consumer) | `/Users/mwiewior/research/git/vepyr` (branch `feat/annotation-target-partitions`) |
| Cache dir (built) | `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged` — holds BOTH `*.lance/` and `parquet.*/` dirs for chr1/2/4 |
| Source cache | `/Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38` |
| Golden VEP refs | `~/workspace/data_vepyr/HG002_annotated_wgs_everything_hgvs_merged.vcf` (+ per-chrom `vep_chr*_merged.vcf` in vepyr `e2e-testing/results/`) |

**IMPORTANT:** A different agent works `fastqc` in the MAIN checkout (`.../datafusion-bio-functions`, branch `feat/fastqc-phase1`). **Stay in the worktree; do not touch the main checkout.** vepyr is a separate repo — editing it is fine.

---

## 2. Commits on the branch (green baseline `7ea6d86`)

Tasks 1–7 (earlier): `6b697da` async feat → `cce36e1` writer props → `be6a276` PageDir/CoalescingAsyncReader port → `fd18efe` encode primitives → `43bb741` variation writer+roundtrip → `16b1005` variation reader → `3e31dec` detect → `41fddc2` read-seam → `ede392a`/`d869df5`/`33e0139` backend-select + hybrid.

This session (Tasks 2 driver + 8–11):
- `954f105` Task 2 build driver (variation shards)
- `82e2e98` fix: `variation_name` must be nullable (dedup nulls it)
- `dd89322` fix: PageDir must return ALL pages a key spans + Task 8 parity gate
- `3b90765` PageDir generalized to u64 keys (INT32 + INT64)
- `72ed78e` translation_sift no-dict Parquet point-lookup + driver + gate
- `0318077` wire Parquet translation_sift store into prepare_contig_context
- `7ea6d86` scan entities (transcript/exon/regulatory/motif/translation_core) dict-Parquet + wiring

---

## 3. Architecture (how the Parquet backend works)

**Module `src/parquet_cache/`:**
- `page_dir.rs` — `PageDir` (footer ColumnIndex/OffsetIndex → position→page directory; **u64 keys**, handles INT32 `start` and INT64 `key`) + `CoalescingAsyncReader` (merges page byte-ranges) + `selection_from_ranges`/`selection_from_offsets`.
- `encode.rs` — variation column encoders: `format_g4` (C `%.4g`, the losslessness linchpin), `presence_boolean` (Int8 flags → non-null Boolean), `encode_af_2array`/`reconstruct_af_group_string` (AF ↔ `List<Utf8>` alleles + `List<Float32>` freqs), `dedup_variation_name`.
- `write.rs` — `point_lookup_writer_properties(schema, sort_cols)` (no-dict, 4KiB/512-row pages, page index, zstd3) shared by variation `(tier,start)` and sift `(key)`; `variation_output_schema`; `encode_variation_batch`; `VariationParquetShardWriter`; `write_variation_parquet`.
- `variation_lookup.rs` — `SinglePathParquetVariationLookup::{open, resolve_and_take}` → `TakenVariationRows` (3-phase: PageDir resolve → exact-offset key match → payload take; reconstructs AF + coalesces `variation_name`).
- `sift.rs` — `write_sift_parquet` + `SinglePathParquetSiftLookup::take_keys(&[u64]) → (RecordBatch, present)` (same contract as Lance `KeyU64LanceLookup`).
- `scan.rs` — `context_writer_properties` (**dict ENABLED**, zstd3, 50k-row groups), `ContextParquetShardWriter`, `read_context_parquet(path, cols)` (projected full scan, tolerates missing cols, quote-stripped — mirrors Lance `scan_projected_existing_columns`).
- `detect.rs` — `PartitionedParquetCache::{detect, variation_path, context_path, base_dir, available_chroms}` (keys off `parquet.variation/chrom_manifest.json`).

**On-disk layout** (per entity, under the cache base dir): `parquet.<entity>/chrom_manifest.json` + `parquet.<entity>/<chrom>.parquet`. Entities: `variation`, `translation_sift`, `transcript`, `exon`, `regulatory`, `motif`, `translation_core`.

**Backend selection** (`annotate_provider.rs`): `PartitionedAnnotationCache::Parquet { variation: PartitionedParquetCache, lance_context: PartitionedLanceCache }` (hybrid). Detection builds both handles. Variation → Parquet `KvLookupExec`. Context (`scan_lance_context_entity`) + SIFT (`prepare_contig_context` ~12590) prefer the Parquet shard when present, else Lance. `cache_format="parquet"` sets `use_parquet_backend`.

**Build drivers** (all in `src/lance_cache/build.rs`, reuse the Lance query/transform helpers): `build_parquet_variation_chrom`, `build_parquet_translation_sift_chrom`, `build_parquet_context_entity_chrom` (transcript/exon/reg/motif), `build_parquet_translation_core_chrom`. NOTE: these are only invoked from `examples/`, **not** wired into `cache_builder.rs` yet (Task 12 should do that).

---

## 4. Commands

```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions/.claude/worktrees/parquet-cache-backend
# Build + test + lint (feature set is lance-cache,cache-builder today)
cargo build   -p datafusion-bio-function-vep --features lance-cache,cache-builder
cargo test    -p datafusion-bio-function-vep --features lance-cache,cache-builder parquet_cache::
cargo clippy  -p datafusion-bio-function-vep --all-targets --features lance-cache,cache-builder -- -D warnings

# Build a Parquet shard for one chromosome (examples; all need --features lance-cache,cache-builder)
CR=/Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38
OUT=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
./target/release/examples/build_parquet_variation_chrom        --cache-root $CR --output-dir $OUT --chrom chr1 --cache-source-type merged --overwrite
./target/release/examples/build_parquet_translation_sift_chrom --cache-root $CR --output-dir $OUT --chrom chr1 --cache-source-type merged --overwrite
./target/release/examples/build_parquet_context_chrom --entity transcript --cache-root $CR --output-dir $OUT --chrom chr1 --cache-source-type merged --overwrite
#   --entity ∈ {transcript, exon, regulatory, motif, translation_core}

# Read-parity gates (Parquet vs Lance)
./target/release/examples/variation_backend_diff --parquet $OUT/parquet.variation/chr1.parquet --lance $OUT/variation.lance/chr1.lance --positions <probes.txt>
./target/release/examples/sift_backend_diff      --parquet $OUT/parquet.translation_sift/chr4.parquet --lance $OUT/translation_sift.lance/chr4.lance --keys <keys.txt>
./target/release/examples/context_backend_diff   --parquet $OUT/parquet.transcript/chr4.parquet --lance $OUT/transcript.lance/chr4.lance --sort-key stable_id
```

**vepyr e2e** (rebuild required after any crate change; ~3–4 min):
```bash
cd /Users/mwiewior/research/git/vepyr
RUSTFLAGS="-C target-cpu=native" uv sync --reinstall-package vepyr
uv run python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend parquet --skip-compare --force
#   --backend lance|parquet ; drop --skip-compare to diff against golden VEP.
#   output: e2e-testing/results/fast_<chrom>/vepyr_<backend>_<chrom>_merged.vcf
# Data-only parity: cmp <(grep -v '^#' lance.vcf) <(grep -v '^#' parquet.vcf)   (timeout ≥180s; 1.6GB files)
```

Probe files used this session live in the scratchpad: `chr1_probe_positions.txt` (8000 warm/cold/multiallelic), `chr4_sift_keys.txt` (6003 keys). Regenerate with pyarrow if gone.

---

## 5. Validation evidence (what "done" means)

- Unit tests: `parquet_cache` 24/24; lib suite 827+/0.
- Storage (chr4): variation 4.2G→2.0G (−52%), sift 363M→237M (−35%), transcript 79M→23M, exon 22M→6.4M, regulatory 15M→3.4M. **translation_core 8.3M→13M (+57%, see §8 follow-up).**
- Read parity gates: variation 42/42 cols byte-identical; sift blobs byte-identical; context live columns identical (transcript 79/79 sorted, translation_core 10/10).
- **Full vepyr e2e (variation + sift + all context on Parquet vs all-Lance): chr1/chr2/chr4 data byte-identical** (323,430 / 331,324 / 307,295 lines).
- Against Ensembl VEP 115 (Parquet backend): chr1/chr2/chr4 all 86 CSQ fields 100%, 0 mismatches.

---

## 6. Gotchas / lessons (read before touching anything)

1. **`PageIndexPolicy` is private in parquet 58** → must load metadata via `ArrowReaderMetadata::load(.., ArrowReaderOptions::new().with_page_index(true))` under `#[allow(deprecated)]`.
2. **`variation_name` is non-nullable in the source but dedup nulls it** → the Parquet output schema declares it nullable; the reader coalesces `coalesce(variation_name, dbsnp_ids)`.
3. **PageDir must return ALL pages a key spans.** Multi-allelic positions whose identical-key run straddles a 512-row page boundary need every page `while min[pg] <= p` (not just the first). This was a real chr1 bug (dropped 7 rows).
4. **Scan-entity row ORDER differs between backends but row SETS are identical** (source `ORDER BY` isn't total; Lance datasets built at a different time). The annotation is **order-insensitive** (engine sorts internally) → e2e byte-identical. The `context_backend_diff` gate row-by-row FAILS on order; use `--sort-key <unique col>` for a set comparison.
5. **Dead columns** `raw_object_json`/`object_hash`/`source_file` differ between Parquet (fresh) and some `.lance` (older serialization) — but they're **unused by annotation** (dead-weight). Compare only live columns.
6. **vepyr conflates two "backend" concepts.** `annotate_vep`'s `backend` arg = the annotation STORE (`AnnotationBackend`, Lance-only by design). The VARIATION backend = `cache_format` in `options_json`. vepyr was passing `cache_format` as the store `backend` arg → `annotate_vep(): backend must be 'lance'; got: parquet`. Fix (already applied in vepyr `src/annotate.rs:153,317`): pass `backend="lance"`, keep `cache_format` in options.
7. **`translation_core` +57% larger on Parquet** — big high-cardinality sequences; dict doesn't help. Acceptable under zero-Lance but worth a follow-up (§8).

---

## 7. TASK 12 — the remaining work (delete Lance, Parquet-only, drop deps)

Full file:line surface is in memory `parquet-task12-lance-deletion-map.md`. Summary:

**Strategy:** delete-in-place. Keep the `lance_cache/` directory, delete lance-coupled code within each file, keep format-agnostic survivors. Optional final cosmetic `git mv lance_cache → cache_common`. **Phase order is forced by dependency direction.**

**Format-agnostic survivors `parquet_cache` needs (keep):** `manifest.rs` (ChromManifest), `af_bundle.rs`, `schema.rs` (VARIATION_* consts), `key_encoding.rs`, `variant_key.rs`; split-keep `ResolvedRowIds`/`PositionPageDirectory`-struct (row_index.rs), `TakenVariationRows`/`ensure_runtime_projection` (variation_runtime.rs), all `build_parquet_*` drivers + query/transform helpers + `variation_projected_schema` (build.rs), `block_on_lance`→rename `block_on` (lookup_exec.rs, tokio-only), `position_predictions_from_batch` + `parse_lance_*_batches` (annotate_provider.rs, agnostic decoders).

**Delete entirely (pure lance):** `lance_cache/write.rs`, `lance_cache/context_runtime.rs`, `PartitionedLanceCache` (partitioned_cache.rs), Lance fns in cache_source.rs, all `build_lance_*` builders + `write_*_to_lance*` + `LanceIndexKind` (build.rs), Lance lookups (`SinglePathLanceVariationLookup`/`StreamingPositionCursor` in variation_runtime.rs; `load_*_from_lance_btree` in row_index.rs), `PositionSlicedLanceSiftStore` + `load_lance_sift_prediction_store_for_chrom` + `read_lance_dataset_schema` (annotate_provider.rs).

- **Phase 1 — collapse selection to Parquet-only (annotate_provider.rs + lookup_exec.rs):** enum `PartitionedAnnotationCache` → struct `{variation: PartitionedParquetCache}` (drop Lance arm, `lance_context`, `as_lance()`); `cache_format` parse accept parquet only; `detect` → `PartitionedParquetCache::detect` only; `scan_lance_context_entity`/`load_lance_contig_context` drop the Lance param+fallback; SIFT store selection → parquet only; `VariationLookupStorage` drop Lance arm, delete `new_lance`/`ensure_lance_lookup`/`lance_lookup_cell`/`lance_cursors`, rewrite `new_parquet`. Then `cargo build` and fix every error the compiler surfaces (~15+ sites). *This is where the 2026-07-04 attempt stopped — the enum edit cascaded; needs to be finished in one green pass.*
- **Phase 2 — delete the now-unreachable Lance code** (the "delete entirely" + "split" lists above).
- **Phase 3 — Cargo.toml + lib.rs:** delete `lance`/`lance-file`/`lance-index` deps + `lance-cache` feature (rename → `parquet-cache` or fold into `cache-builder`); verify the optional deps `arrow-ipc`/`lz4_flex`/`zstd`/`ahash`/`hashbrown`/`nohash-hasher`/`rapidhash`/`smallvec` have no other consumer before dropping; fix `#[cfg(feature="lance-cache")]` gates in lib.rs; delete `build_lance_*` examples, repoint `build_parquet_*`/`*_diff` `required-features`; repoint `cache_builder.rs::build_entity` to the `build_parquet_*` drivers.
- **Phase 4 — tests + gates:** rework/delete ~40 tests (full list in the memory map). **Acceptance:** `cargo tree -p datafusion-bio-function-vep --features cache-builder | grep -i lance` → EMPTY; `cargo test` + `cargo clippy --all-targets -- -D warnings` clean; rebuild vepyr (update its `Cargo.toml` feature name from `lance-cache` → the new one) + chr1/2/4 e2e byte-identical to the last all-Lance baseline.

**Caution:** after Phase 1, annotation only works with a FULL Parquet cache. Only chr1/2/4 are built to Parquet — other contigs live only in Lance. Option (b) (chosen by the user) accepts rebuilding the cache separately; do not expect other contigs to annotate until rebuilt.

---

## 8. Open follow-ups (not Task 12)

- **`translation_core` Parquet +57% vs Lance** (chr4 8.3M→13M). Big high-cardinality sequence columns where dict is useless. Try: no-dict + higher zstd on the four sequence columns (`translation_seq`, `cds_sequence`, `translation_seq_canonical`, `cds_sequence_canonical`), or column-specific encoding. Measure vs Lance.
- **vepyr wiring is uncommitted** on branch `feat/annotation-target-partitions`: `Cargo.toml:34` path dep points at the worktree (revert with `git checkout Cargo.toml` + rebuild to restore main), `src/vepyr/__init__.py:644` + `src/annotate.rs:43` accept `"parquet"`, `run_annotation_fast.py:34` `BACKENDS += "parquet"`, `src/annotate.rs:153/317` decouple backend/cache_format. Decide whether to commit these on vepyr's branch.
- **Filed vepyr issue #19** — `compare_vcfs` CSQ-count accounting drops empty-CSQ variants (`match + mismatch ≠ variants_compared`); backend-independent, not a Parquet bug.

---

## 9. Related memory files

- `parquet-backend-execution-state.md` — running execution log (Tasks 1–11 detail).
- `parquet-task12-lance-deletion-map.md` — the full Task 12 deletion surface with file:line anchors.
