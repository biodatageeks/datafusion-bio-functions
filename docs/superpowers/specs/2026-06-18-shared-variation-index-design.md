# Shared Variation-Lookup Index — Design

**Date:** 2026-06-18
**Repo:** `datafusion-bio-functions` (`datafusion/bio-function-vep`). Downstream: `polars-bio`, `vepyr`.
**Status:** Design approved; ready for implementation plan.
**Related:** RSS root-cause analysis in memory `vep-t8-rss-root-cause.md`; builds on the within-contig parallel annotation path (`docs/superpowers/specs/2026-06-17-within-contig-parallel-annotation-design.md`).

---

## 1. Problem

In the threads>1 (within-contig parallel) VEP annotation path, peak RSS scales ~4.4× from t1 to t8 (chr1 200k `--everything`, lance backend: **~5 GB → ~23 GB**). Empirical decomposition (knob sweeps + input-size sweeps + code trace) proved:

- The growth is **~2.6 GB of FIXED overhead per lookup partition**, invariant to input row count (25k vs 200k rows → same RSS), region width, and `--everything`. Ruled out: output channel (`out_cap`), mimalloc retention, colocated/AF payload, data-proportional accumulation, shared context (Arc-shared, not multiplied).
- The fixed cost is the **`PositionRowIdIndex`** — the variation cache's `start` BTree index, materialized **once per lookup partition**. Each partition's `KvLookupStream` (`cache_exec.rs:497`) starts with an empty `lance_chroms` map and, on first probe (`ensure_cold_parquet_lookup`, `cache_exec.rs:2589`), independently calls `SinglePathLanceVariationLookup::open` (`lance_cache/variation_runtime.rs:27`) → `load_start_index_from_lance_btree` → `load_u32_btree_index` (`lance_cache/row_index.rs:246`).
- `PositionRowIdIndex { positions: Vec<u32>, row_ids: Vec<u64> }` (`row_index.rs:18`) holds **one entry per cache row** (chr1 = 88.15M rows): ~1.0 GB resident (336 MB + 673 MB) plus a ~2.4–2.6 GB transient at build (intermediate `pairs: Vec<(u32,u64)>` alive with the final vectors). All N partitions first-probe at startup, so their transients overlap → `N × ~2.6 GB ≈ 18 GB` of inflation.

The index is **read-only reference data**; only the per-partition walk `cursor` is genuinely per-partition. So building it N times is pure waste.

## 2. Goal

Build each contig's lance variation lookup (dataset + `PositionRowIdIndex`, and any other lance BTree indexes) **once**, and share it across all workers via `Arc`, uniformly for threads=1 and threads>1.

**Success criteria:**
- Byte-identical annotation output at threads ∈ {1,2,4,8} and vs the pre-change serial output.
- Peak RSS at t8 drops from ~23 GB toward **~6–7 GB** (`/usr/bin/time -l`); t1 unchanged.
- Lib unit-suite failure set unchanged vs baseline (35 pre-existing failures, unrelated).

## 3. Non-Goals (explicit, out of scope)

- **Backend removal.** The `Fjall` and `IndexedParquet` backends and their `warm_chroms` / `cold_parquet_chroms` / `variant_bloom_chroms` machinery are untouched (separate follow-up pending cross-repo confirmation they're unused). Those maps stay empty in the lance path and contribute nothing to lance RSS.
- **Region-scoped / sliding index loading.** Considered and rejected for now: range-scoping the index per partition only helps with eviction, saving ~0.8 GB more for materially more code (BTree range pushdown + per-window load + eviction). Against the ~18 GB the shared build already recovers, it's YAGNI; revisit only if the shared ~1.0 GB index ever bites.
- **Variation row data.** Already on-demand via `take_rows` per buffer; unchanged.
- **SIFT lookup (checked — already correct, no change).** The lance SIFT store (`load_lance_sift_prediction_store_for_chrom`, `annotate_provider.rs:3349`, used for `--everything` + lance) *does* carry a sizable lance BTree index analogous to `PositionRowIdIndex` — one of `PositionSlicedLanceSiftStore`/`KeyU64LanceLookup` (`U64RowIdIndex`), `LanceBinarySiftPredictionStore`/`TranscriptIdLanceLookup` (`StringRowIdIndex`), or a fully in-memory `InMemorySiftPredictionStore`. **But it is built once per contig and lives in the `Arc`-shared context** (`SharedContigAnnotationContext.sift_prediction_store: Arc<dyn SiftPredictionStore>`), referenced by every worker via `&shared.sift_prediction_store` (`:11326`) — *not* cloned into per-partition `KvLookupStream` state; even its lazy decode cache is `Arc<Mutex<HashMap>>` (shared). That is exactly the shared-per-contig model this spec applies to variation, which is why SIFT contributed nothing to the per-partition RSS growth. The **variation lookup is the one lance lookup that lives in per-partition stream state** (`KvLookupStream.lance_chroms`) instead of the shared context — the asymmetry this design corrects. (The parquet `SiftDirectReader` at `:1154` is the non-lance path and is likewise small + `Arc`-shared.) No SIFT work needed.

## 4. Architecture

Hoist the lance variation lookup from **per-partition-lazy** to **per-contig-shared**, mirroring how `base_transcripts: Arc<Vec<TranscriptFeature>>` is already loaded once per contig and shared.

```
prepare_contig_context(chrom)            [async, once per contig]
  ├─ register variation table
  ├─ load base_transcripts/exons/...     (existing Arc-shared context)
  ├─ build Arc<SinglePathLanceVariationLookup>   <-- NEW: open dataset + load BTree index ONCE
  └─ build KvLookupExec carrying that Arc

KvLookupExec  (Arc::clone'd to each of N workers)
  └─ execute(partition_id)               [N times]
       └─ KvLookupStream {
            lance_chroms: <shared Arc, cloned>,   <-- was: empty, self-built on first probe
            lance_cursors: <own, per-partition>,  <-- unchanged
            ...
          }
```

- The lookup is **immutable after build**. `resolve_and_take(&self, starts, cursor: &mut usize)` (`variation_runtime.rs:46`) already takes `&self`; the call site already uses `lance_chroms.get(&chrom)` (immutable, `cache_exec.rs:3642`). The only mutable state is the external per-stream `cursor`.
- `lance::Dataset` reads (`take_rows`, `&self`, async) are `Send+Sync` and concurrency-safe — lance's intended usage. Sharing one `Dataset` also collapses lance's internal per-dataset read caches from N copies to 1 (secondary RSS bonus).

## 4a. Worker flow after the change (illustration)

```
                          prepare_contig_context(chrom)        [async, ONCE per contig]
                          ─────────────────────────────
   reads each lance BTree index once, materializes the shared per-contig context:

        ┌───────────────────────  Arc<SharedContigAnnotationContext>  ───────────────────────┐
        │   base_transcripts: Arc<Vec<TranscriptFeature>>          (whole chrom, shared)       │
        │   exons / translations / regulatory: Arc<Vec<…>>         (shared)                     │
        │   sift_prediction_store: Arc<dyn SiftPredictionStore>    (BTree index, shared)  ←already
        │   ┌─────────────────────────────────────────────────────────────────────────────┐   │
        │   │  lance variation lookup  Arc<SinglePathLanceVariationLookup>      ← NEW(shared)│  │
        │   │     dataset: lance::Dataset            (read-only, concurrent take_rows)       │  │
        │   │     index:   PositionRowIdIndex        (~1.0 GB, 88M (pos→row_id), READ-ONLY)  │  │
        │   └─────────────────────────────────────────────────────────────────────────────┘   │
        └───────────────────────────────────────────────────────────────────────────────────┘
                          │ Arc::clone (refcount bump — NO data copy)
        ┌─────────────────┼─────────────────┬─────────────────┬──────────────── … ┐
        ▼                 ▼                 ▼                 ▼                     ▼
   ┌─────────┐       ┌─────────┐       ┌─────────┐       ┌─────────┐         ┌─────────┐
   │worker 0 │       │worker 1 │       │worker 2 │       │worker 3 │   …     │worker 7 │
   │partition│       │partition│       │partition│       │partition│         │partition│
   │ [a0,a1) │       │ [a1,a2) │       │ [a2,a3) │       │ [a3,a4) │         │ [a7,end)│
   ├─────────┤       ├─────────┤       ├─────────┤       ├─────────┤         ├─────────┤
   │PER-WORKER (small, NOT shared):                                                     │
   │ lance_cursor=c0 │ cursor=c1 │     │ cursor=c2 │     │ cursor=c3 │       │ cursor=c7│
   │ window_buffer   │ window_buf│     │ window_buf│     │ window_buf│       │ window_buf
   │ colocated_sink  │ …         │     │ …         │     │ …         │       │ …        │
   └────┬────┘       └────┬────┘       └────┬────┘       └────┬────┘         └────┬────┘
        │  every 5000-row buffer, sliding forward along its sub-range:             │
        │   resolve_and_take(starts, &mut cursor)                                   │
        │     ├─ index.resolve(...)         → reads SHARED index at cursor offset    │
        │     │                               (read-only; only c_i mutates)          │
        │     └─ dataset.take_rows(row_ids) → fetches ONLY matched variation rows    │
        │                                     on-demand (concurrent, shared dataset) │
        │   + base_transcripts: scan shared Vec for buffer range                     │
        │   + sift_prediction_store: query shared store                              │
        ▼                                                                            ▼
   annotated RecordBatches ──────────────►  drained in partition order (0,1,2,…,7)  ──► output
                                            (ParallelContigState.current)

   RESIDENT lance index = 1 × ~1.0 GB (shared)   — was N × ~1.0 GB + N overlapping ~2.6 GB build transients
```

**Key invariants shown:** one index built once (no per-worker rebuild, no build race); each worker holds only a cheap `Arc` handle plus its own `cursor`; all N cursors walk the *same immutable* index at different offsets as their buffers slide; variation rows and SIFT/transcript queries are all reads against shared, read-only structures. The index drops when the contig's context drops (next contig builds its own).

## 5. Components & interface changes

**`SinglePathLanceVariationLookup`** (`lance_cache/variation_runtime.rs:20`) — unchanged internally. Already `{ dataset, projection, index }`, all read-only post-`open`; `resolve_and_take` is `&self`. It becomes shared as `Arc<SinglePathLanceVariationLookup>`.

**`KvLookupExec`** (`cache_exec.rs:76`) — gains a field carrying the pre-built shared lance lookup for the contig, e.g.:
```rust
#[cfg(feature = "lance-cache")]
lance_lookup: Option<Arc<SinglePathLanceVariationLookup>>,  // built once, per contig chrom
```
Populated when the exec is constructed for the contig (lance backend only). Because `KvLookupExec::new_lance` (`cache_exec.rs:260`) is currently **sync** while opening the dataset + reading the BTree is **async**, the build happens in the async plan-construction path (`TableProvider::scan`, reached from `prepare_contig_context` at `annotate_provider.rs:12397`, which already has the chrom via `set_vcf_filter` and the cache root via `VariationLookupStorage::Lance { cache_root }`). The pre-built `Arc` is threaded into the exec constructor. (Exact threading is an implementation detail for the plan.)

**`KvLookupStream`** (`cache_exec.rs:497`):
- `lance_chroms: HashMap<String, Option<SinglePathLanceVariationLookup>>` → holds the shared `Arc<SinglePathLanceVariationLookup>` (cloned from the exec at `execute()` time, keyed by the contig chrom), instead of being built per-stream.
- `KvLookupExec::execute` (init at `cache_exec.rs:1903`) seeds `lance_chroms` with the shared `Arc` instead of `HashMap::new()`.
- `ensure_cold_parquet_lookup` (`cache_exec.rs:2589`): the lance branch's lazy `SinglePathLanceVariationLookup::open(...)` + `lance_chroms.insert(...)` is **removed** — the entry is already present (shared). Non-lance branches unchanged.
- `lance_cursors: HashMap<String, usize>` — **unchanged**, stays per-stream (the only per-partition lance state).
- `resolve_and_take` call site (`cache_exec.rs:3642`) — unchanged (already `lance_chroms.get` + per-stream cursor).

**threads=1 path.** Uses the same shared build (built once in `prepare_contig_context`, one worker clones the one `Arc`). This unifies the code path with threads>1 rather than branching, and removes the lazy first-probe build there too — no behavioral change (output identical; build moves slightly earlier, into context prep).

## 6. Data flow / lifecycle

- **Per contig**, `prepare_contig_context` builds one `Arc<SinglePathLanceVariationLookup>` for the contig's chrom and embeds it in that contig's `KvLookupExec`. The streaming `ContigAnnotationExec` processes contigs sequentially, so **only one contig's index is resident at a time** — no cross-contig accumulation. The `Arc` (and its ~1.0 GB index) drops when the contig's exec/streams drop.
- **Per partition**, `execute(partition_id)` clones the `Arc` (refcount bump, no data copy) into the stream and keeps its own `lance_cursors`. All N cursors walk the same immutable index at their own offsets as their 5000-row buffers slide forward.
- **Per buffer**, `resolve_and_take(starts, &mut cursor)` resolves the buffer's variant positions against the shared index (read-only) and `take_rows` fetches only the matched variation rows on-demand (unchanged).

## 7. Concurrency safety

- Index: read-only after build; `resolve_sorted_positions_from_cursor(&self, …, cursor)` (`row_index.rs:97`) mutates only the caller-owned `cursor`. N concurrent readers over one `Arc<index>` with N independent cursors → no contention, no lock.
- Dataset: `lance::Dataset::take_rows(&self, …)` async, `Send+Sync`, concurrency-safe; shared via the `Arc<lookup>`.
- No interior mutability, no `Mutex`/`RwLock` introduced.
- The eager build sidesteps the lazy-build thundering-herd race (N simultaneous first-probes each building a ~2.6 GB index) entirely — the index is built once before fan-out and only an immutable `Arc` fans out.

## 8. Error handling

- Build failure (dataset open / BTree read) surfaces in `prepare_contig_context` as a `DataFusionError` (same error types as today's `SinglePathLanceVariationLookup::open`), failing the contig before workers spawn — strictly earlier and cleaner than today's per-stream first-probe failure.
- If the contig has no variation table, the existing "skip contig" path is unchanged (no lookup built).

## 9. Testing & verification gates

1. **Byte-identical parity (mandatory).** threads ∈ {1,2,4,8} all produce identical data rows (`diff <(grep -v '^#' …)`), and t1 equals the pre-change serial output. Repro: the `bench_annotate_vcf` chr1 200k `--everything` lance run (same harness as the transcript-clone work).
2. **Peak-RSS gate.** `/usr/bin/time -l` t8 peak RSS drops from ~23 GB to **~6–7 GB**; t1 unchanged (~5 GB). A/B vs the pre-change binary on the same machine.
3. **Regression set.** Lib unit-suite failure set identical to baseline (923 pass / 35 pre-existing fail). Use the e2e parity gate as the authority (the branch's 35 unit failures are pre-existing and unrelated — see `transcript-clone-elimination.md`).
4. **Startup wall-time check (informational).** The redundant 7× BTree decode is eliminated; expect t8 startup/`context_load` no worse, likely slightly better. Not a gate.

## Results (implemented 2026-06-18; plan `docs/superpowers/plans/2026-06-18-shared-variation-index-plan.md`)

Implemented via a shared `Arc<tokio::sync::OnceCell<Arc<SinglePathLanceVariationLookup>>>` on `KvLookupExec`, cloned into each `KvLookupStream`; `get_or_try_init` single-flights the build (reuses the existing path/projection logic verbatim — lowest parity risk). Per-partition `lance_cursors` unchanged.

**Parity (mandatory gate) — PASS:** 0 mismatches in data rows at threads ∈ {1,2,4,8}, and 0 mismatches vs the pre-change serial output. Lib suite 923 pass / 35 pre-existing failures (identical set; unrelated).

**RSS + timing (chr1 200k `--everything`, lance, mimalloc; A/B same machine, min of 2 reps):**

| threads | BEFORE wall / RSS | AFTER wall / RSS |
|---|---|---|
| 1 | 16.8s / 5.30 GB | 16.8s / 5.33 GB |
| 2 | 7.8s / 8.00 GB | 7.7s / 7.45 GB |
| 4 | 6.0s / 12.88 GB | 5.5s / 8.03 GB |
| 8 | 5.9s / **20.47 GB** | 4.3s / **7.90 GB** |

- **t8 peak RSS 20.47 → 7.90 GB** (−12.6 GB, ~2.6×). RSS is now near-flat across thread counts (5.3→7.5→8.0→7.9) vs the prior steep linear growth — the per-partition 88M-entry index duplication is eliminated; the residual ~2.6 GB above t1 is per-partition working buffers/output, not the index.
- **Bonus: t8 wall 5.9 → 4.3s (~27% faster)** — the 7 redundant BTree index decodes at startup are gone.
- t1 unchanged (one worker builds once either way).
- t8 landed slightly above the ~6–7 GB projection (7.9 GB) — the remaining footprint is shared context (~4.5 GB) + one shared index (~1 GB) + per-partition buffers + the single build transient.

## 10. Risks

- **IO concurrency shift.** N partitions now `take_rows` against one shared `Dataset` instead of N independent readers. Not a correctness risk (lance pools/serializes internally); worst case is mild read contention. Caught by the parity gate (output) and observed via wall-time.
- **Build-timing shift.** The index builds during contig prep instead of first probe — a few hundred ms moved earlier in the contig. No output impact; covered by parity + wall-time.
- **Multi-contig keying.** The shared lookup is per-contig (one chrom). `lance_chroms` keyed by chrom must match the chrom the contig's stream probes. Low risk (one chrom per contig exec), covered by parity across a multi-contig input if exercised.
