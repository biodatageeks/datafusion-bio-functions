# Variation Row-ID Dedup Cache Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop the cold-tier Lance variation lookup from `take_rows`-ing the same variation row more than once. Today the per-batch streaming path re-takes overlapping rows across VCF batches, fetching **~735k rows for chr1 when the unique set is ~402k (1.83×)** — directly inflating the variation read (17.4s at true 1‑CPU). Target: take each unique row at most once (≈402k), roughly halving the variation take, with **byte-identical** annotation output.

## Background / Evidence (measured, chr1, HG002, true 1‑CPU: `LANCE_CPU_THREADS=1 LANCE_IO_THREADS=1 RAYON_NUM_THREADS=1`)

- Standalone lance_sandbox benchmark (single global pass, same `variation.lance/chr1.lance`, 44 projected cols): **402,085 unique rows** for 310,668 selected positions; `take_rows` ≈ 9.3s cold → 6.1s warm.
- vepyr e2e variation lookup: **735,017 rows** across **79** `take_rows` calls (one per VCF batch flush); `variation_take` total **17.4s**. Matched positions ≈ 315,586 — essentially the same as the benchmark, so the extra rows are **re-takes, not new data**.
- Per-row take cost is at parity (vepyr 23.7 vs benchmark cold 23.1 µs/1k rows); the gap is entirely the 1.83× row count.

Root cause (code): `kv_cache/cache_exec.rs` cold-tier lance branch (~L3638-3661) dedups **positions per batch only** (`starts.sort_unstable(); starts.dedup()`), then `lookup.resolve_and_take(&starts, cursor)`. There is **no cross-batch row_id dedup**, and `PositionRowIdIndex::resolve_sorted_positions_from_cursor` **rewinds the cursor** when extended probes overlap a later window (see test `rewinds_cursor_when_extended_probes_overlap_later_window`) — so positions/rows taken in an earlier batch are re-resolved and re-taken. Adjacent batches that share a position (incl. extended/deletion probes) also re-take.

## Architecture

Insert a **contig-scoped `row_id → row` cache** between resolve and take, so `take_rows` only fetches row_ids not already cached. Consumers are served from the cache (cached ∪ newly-taken). The cache is bounded by **eviction on the streaming frontier**: variants arrive in genomic order, so a cached row whose variation `end < current_window_start` will never be referenced again and is dropped — the same pattern the SIFT/translation caches already use. If a needed row_id was evicted (frontier moved past then a rewind re-references it), fall back to re-taking it (correctness preserved, rare).

Key invariant: **the set of rows handed to the colocated-matching/consumer code for a given batch is exactly the same as today** — only the source (cache vs fresh take) changes. This is what keeps output byte-identical.

## Tech Stack

Rust 2024, Lance 7.0.0 `Dataset::take_rows`, existing `kv_cache::cache_exec` cold-tier lance path, `PositionRowIdIndex` (`lance_cache::row_index`), `lance_start_row_map`. Profiling via `VEP_LANCE_PROFILE` (`variation_take`) and `VEP_PROFILE` pipeline `cold_tier_load`.

## File Structure

- Modify: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`
  - The lance cold-tier flush (~L3626-3661): resolve → split row_ids into cached/missing → `take_rows(missing)` → insert → build the batch handed to `lance_start_row_map` from the full (cached ∪ taken) set.
  - Add a `VariationRowCache` (row_id → row payload) field on the cold-tier lance executor state, plus frontier eviction hooked to the existing window/cursor advance.
- Possibly modify: the lance variation lookup type that owns `resolve_and_take` — split into `resolve_row_ids(&starts, cursor) -> Vec<u64>` and `take_row_ids(&[u64]) -> RecordBatch` so the cache can sit between them (keep `resolve_and_take` as a thin wrapper for any other caller).
- Add: profiling counters — rows requested vs rows taken vs cache hits — emitted on the existing `variation_take` / kv-detail lines so the win is measurable.

## Task 1: Split resolve from take

- [ ] Factor the lance variation lookup so position→row_id **resolution** is separable from **`take_rows`**: `resolve_row_ids(&sorted_starts, cursor) -> ResolvedRowIds` and `take_row_ids(&[u64], projection) -> RecordBatch`. Keep `resolve_and_take` as `take_row_ids(resolve_row_ids(...))` for unchanged callers.
- [ ] Unit test: resolve+take composition returns the same batch as the current `resolve_and_take` for a synthetic dataset.

## Task 2: VariationRowCache + dedup take

- [ ] Add `VariationRowCache`: `HashMap<u64 /*row_id*/, RowSlot>` where `RowSlot` holds the decoded row (or a reference into a retained batch) plus the row's genomic `end` for eviction. Provide `get_missing(&[u64]) -> Vec<u64>`, `insert_batch(row_ids, batch)`, `assemble(&[u64]) -> RecordBatch` (rebuild a batch in the requested row_id order from the cache).
- [ ] In the cold-tier flush: `resolve_row_ids` → `missing = cache.get_missing(row_ids)` → `take_row_ids(missing)` → `cache.insert_batch(...)` → `batch = cache.assemble(row_ids)` → `lance_start_row_map(&batch)`. The assembled batch must match (schema, row order, values) what `resolve_and_take` produced before.
- [ ] Assert/justify row order: today `resolve_and_take` returns rows in the order lance yields them for `take_rows(row_ids)`; `assemble` must reproduce that order (request order) so `lance_start_row_map` keys identically.
- [ ] Unit test: two overlapping batches (shared positions / a rewind) take each unique row_id exactly once; assembled batches equal the non-cached path.

## Task 3: Frontier eviction (bounded memory)

- [ ] Track the current streaming window start (already available where SIFT/translation windows are evicted). After each batch, drop cache entries whose `end < window_start`.
- [ ] If `get_missing` is asked for a row_id that was evicted, it is simply treated as missing and re-taken (correctness preserved). Add a counter for evicted-then-refetched.
- [ ] Test: eviction frees entries behind the frontier; a contrived rewind past the frontier re-takes correctly.

## Task 4: Profiling + validation

- [ ] Emit `rows_requested / rows_taken / cache_hits / refetched` on the kv-detail and `variation_take` profile lines.
- [ ] `cargo clippy --lib` clean; `cargo fmt`.
- [ ] **Parity gate (acceptance):** chr1 (and chr4) e2e CSQ diff vs current — must stay **100% on all shared fields** (this is the gate; the change must not alter output).
- [ ] Re-profile true‑1‑CPU chr1: `variation_take` rows should drop from ~735k toward ~402k and time from ~17.4s toward ~9–10s (cold); confirm `cache_hits` ≈ 333k.

## Risks & Notes

- **Output parity is the whole point** — the assembled batch must be order- and value-identical to the current take. The parity gate is non-negotiable; if any field diverges, the assemble/order logic is wrong.
- **Memory**: bounded by the frontier eviction; worst case a very dense region holds more rows transiently. Cap the cache and treat a cache miss as a re-take (degrades to today's behaviour, never incorrect).
- **Cross-chrom**: cache is per-contig (matches `lance_chroms` scoping); cleared/!reused across contigs.
- **Orthogonal to** the probe-expansion (727k probes from 319k positions): those mostly *miss*, so they cost lookup CPU but not rows — out of scope here; a separate probe-dedup could cut `pair_compare`/resolve time later.
- **Not** the SIFT path — SIFT position-sliced takes are already deduped/memoized by key.
