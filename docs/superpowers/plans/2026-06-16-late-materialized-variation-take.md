# Late-Materialized (Two-Phase) Variation Take — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or superpowers:executing-plans. Steps use `- [ ]` checkboxes.
> **Phase 0 is a GO/NO-GO measurement gate.** Do not implement Phases 2–3 until Phase 0 proves the win.

**Goal:** Decode the heavy variation columns (`variation_name` + the 27 frequency columns + `clin_sig*`/`somatic`/`pheno`/`pubmed`) only for rows that actually match, instead of for all ~735k taken rows — cutting cold-tier decode by skipping the heavy-column miniblocks that are touched *only* by speculative extended-probe positions.

**Architecture:** Today the Lance variation lookup does one `take_rows` with the full projection, then matches. Two-phase: **Phase A** takes only the 4 light matcher columns (`start,end,allele_string,failed`) for all resolved row_ids and runs `compare_existing_variant_alleles` to find the matched row_ids; **Phase B** does a second `take_rows` of *only the matched row_ids* with the heavy projection, and builds the colocated entries from that. Output is byte-identical (same rows, same values; only the decode set shrinks).

**Tech Stack:** Rust 2024, Lance 2.1 `take_rows`, Arrow 58, DataFusion 53.

**Key prior facts:** rows taken 735,017; primary_matches 311,435; colocated_allele_rows 741,592; cold decode (post-CacheOnly) `cold_tier_load` 10.12s of `variation_lookup` 16.85s @ true 1-CPU. The matcher needs only `start,end,allele_string,failed`; the matched allele's row position is `MatchedVariantAllele.b_index`. **Caveat that Phase 0 resolves:** matched rows share miniblocks with their wrong-allele siblings *and* with ±32 bp deletion-shift neighbors, so the saving comes only from speculative positions with *zero* matches — possibly few distinct miniblocks.

---

## File structure

- `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs` — add a heavy-column `take_rows(row_ids, projection)` entry and split the light projection used by `resolve_and_take`.
- `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs` — restructure the Lance take+match block (~L3621–3826): Phase A match on the light batch, collect matched row_ids + per-match context, Phase B heavy take, build colocated entries from the heavy batch.
- Measurement-only (Phase 0): a throwaway dump env var + a pylance analysis script under `/tmp` (not committed).

---

## Phase 0 — Measurement gate (GO/NO-GO, ~1 hour, no production code)

### Task 0.1: Dump the *matched* row_ids from a real e2e run

**Files:**
- Modify (temporary, revert after): `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs` — at the point a row is confirmed a primary match (search `matched_allele_rows.push(allele_idx)` and the colocated `sink_value.entries.push`), append the row's resolved `row_id` to a file behind `VEP_DUMP_MATCHED_ROWIDS`.

- [ ] **Step 1:** Add an env-gated append (mirror the earlier `VEP_DUMP_ROWIDS` dump in `variation_runtime.rs`, git history commit on this branch). Write the **matched** lance `row_id` (fragment-addressed u64) for each row that produces a colocated entry, one per line.
- [ ] **Step 2:** Rebuild vepyr: `cd ~/research/git/vepyr && RUSTFLAGS="-C target-cpu=native" uv sync --reinstall-package vepyr`
- [ ] **Step 3:** Run chr1 e2e with both dumps:
```
env -u VIRTUAL_ENV -u CONDA_PREFIX LANCE_CPU_THREADS=1 LANCE_IO_THREADS=1 RAYON_NUM_THREADS=1 \
  VEP_DUMP_ROWIDS=/tmp/taken_rowids.txt VEP_DUMP_MATCHED_ROWIDS=/tmp/matched_rowids.txt \
  uv run python run_annotation_fast.py chr1 --cache merged --forks 0 --force --backend lance
```
(from `~/research/git/vepyr/e2e-testing/scripts`)
- [ ] **Step 4:** Revert the temporary dump edit (`git checkout -- datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`), rebuild vepyr.

### Task 0.2: Compute the heavy-column miniblock saving

**Files:** `/tmp/lever_b_gate.py` (pylance; use `research/lance_encoding_sandbox/.venv-pylance/bin/python`)

- [ ] **Step 1:** Map both row_id files to global indices (frag = `id>>32`, off = `id & 0xffffffff`, global = cumulative-fragment-offset + off — same as the earlier analysis). Then, with `CHUNK=4096`, compute distinct miniblocks `(fragment, off//CHUNK)` for **(a)** all taken rows and **(b)** matched rows only, split by tier (hot < 3,628,123 cold ≥).
- [ ] **Step 2:** Compute the **gate metric**: `cold_miniblocks(matched) / cold_miniblocks(taken)`. Also time, at 1-CPU, `ds.take(taken_gidx, HEAVY_COLS)` vs `ds.take(matched_gidx, HEAVY_COLS)` (HEAVY = `variation_name` + 27 AF + `clin_sig,clin_sig_allele,somatic,phenotype_or_disease,pubmed`).
- [ ] **Step 3: GO/NO-GO.**
  - **GO** if matched-only decodes ≲ ~70% of taken (i.e. ≥ ~30% heavy-decode reduction): proceed to Phase 1.
  - **NO-GO** if matched-only ≈ taken (miniblocks shared, < ~15% saving): **stop** — record the result in this file and close Lever B. Do not build the two-phase take.

---

## Phase 1 — Two-phase Lance take (only if Phase 0 = GO)

### Task 1.1: Add light + heavy projections to the Lance lookup

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs`

- [ ] **Step 1: Write the failing test** (extend the existing `lookup_resolves_start_rows_with_take_rows` test) asserting a new method `resolve_and_take_light(&starts, cursor)` returns a batch containing exactly `["start","end","allele_string","failed"]`, and `take_heavy(&row_ids)` returns the heavy projection for given row_ids.
- [ ] **Step 2: Run → FAIL** (methods missing).
- [ ] **Step 3: Implement.** Split the stored projection: `light_projection = ["start","end","allele_string","failed"]`; `heavy_projection = ensure_runtime_projection(projection) minus the light set`. `resolve_and_take_light` = current `resolve_and_take` but with `ProjectionRequest::from_columns(light)`; returns `ResolvedRowIds` + light batch. `take_heavy(&[u64]) -> RecordBatch` = `dataset.take_rows(row_ids, heavy_projection)`. Keep `resolve_and_take` for the non-two-phase callers (or migrate all).
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit** `feat(vep): light/heavy split for two-phase variation take`.

### Task 1.2: Restructure the Lance match block into two phases

**Files:**
- Modify: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs` (~L3621 lance take; ~L3500–3603 the per-probe match + colocated build).

- [ ] **Step 1:** Phase A — replace `resolve_and_take` with `resolve_and_take_light`; build `lance_start_row_map` from the light batch; run the existing per-probe coordinate/allele matching (`compare_existing_variant_alleles`) using only `allele_string`/coords from the light batch. For each row that matches, record `(row_id, probe context, matched_alleles, output-allele info)` into a `Vec` (do **not** build the colocated entry yet). Collect the unique matched `row_id`s in order.
- [ ] **Step 2:** Phase B — `let heavy = lookup.take_heavy(&matched_row_ids)?;` build a `row_id -> heavy-batch-row` map. For each recorded match, read `variation_name/clin_sig/clin_sig_allele/somatic/pheno/pubmed/af_indices` from the heavy batch row and push the `ColocatedCacheEntry` (unchanged struct). The `matched_alleles`/output context come from the Phase-A record.
- [ ] **Step 3:** Ensure the `null_append` / no-match path is unchanged (rows with no match still null their cache columns — they're simply absent from `matched_row_ids`).
- [ ] **Step 4: Build + run the variation unit tests** (`cargo test -p datafusion-bio-function-vep --features lance-cache variation`).
- [ ] **Step 5: Commit** `feat(vep): two-phase late-materialized Lance variation take`.

---

## Phase 2 — Parity + decode verification (acceptance gate)

### Task 2.1: e2e CSQ parity + decode

- [ ] **Step 1:** Rebuild vepyr; run chr1 e2e (true 1-CPU, full profiling) per the standard command.
- [ ] **Step 2: GATE:** `ALL 86 shared CSQ fields match at 100%`, 0 count/order mismatches. (Two-phase must not change output — same matched rows, same values.) If any mismatch: STOP — likely a Phase-B row-mapping bug (matched row_id → heavy batch row).
- [ ] **Step 3:** Record `cold_tier_load` / `variation_lookup` deltas vs baseline 10.12s / 16.85s; expect the Phase-0-predicted reduction.
- [ ] **Step 4:** Repeat chr4 (parity 100%).
- [ ] **Step 5:** Append results to this file; commit.

---

## Notes / risks
- **Phase 0 is the whole bet.** If matched-only touches ~the same cold miniblocks as all-taken (because no-match speculative positions share blocks with matched ones), Lever B yields little — stop at Phase 0.
- **Two extra round-trips:** Phase B is a second `take_rows`; the resolve/IO of Phase A is cheap (4 narrow columns). Net win only if Phase-B heavy decode shrinks more than Phase-A light decode adds.
- **Composes with float-encoding** (if revisited): independent — float makes each heavy decode cheaper; this decodes fewer of them.
- **Scope:** Lance path only. fjall/parquet readers keep the single-pass path (gate the two-phase code on `warm_cold_backend.is_lance()`).
- **Row-id stability:** Phase B re-takes by the same `row_id`s Phase A resolved within the same batch — no cross-batch concern.
