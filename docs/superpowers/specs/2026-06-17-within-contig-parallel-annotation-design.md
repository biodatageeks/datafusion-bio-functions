# Within-Contig Parallel Annotation via Position Windows

**Date:** 2026-06-17
**Status:** Design — REVISED after implementation findings (see §0)
**Goal:** Make annotation of a *single contig* (e.g. chr1) as fast as possible by parallelizing within it, beating the structural weaknesses of Ensembl VEP's `--fork` model.

---

## 0. Implementation findings & REVISED design (supersedes §2–§7 where conflicting)

The original design (per-window work-stealing) was partly implemented and benchmarked against the real chr1 lance cache. Findings:

**Implemented & committed (verified byte-identical):**
- Parallel **annotation**: each ~`input_buffer_size`-row window annotated by a fresh worker on `spawn_blocking`, drained FIFO (`AwaitingWindow` state, `inflight` queue). `annotate_window_owned`.
- Parallel **output formatting**: `stream_vcf_partition_to_writer` pipelines `threads` `format_vcf_body_chunk` tasks.
- Single `threads` knob (`AnnotateVcfConfig.threads` → options-json `"threads"` → `ContigAnnotationConfig.annotation_threads`).

**Measured (chr1 HG002, --everything, lance): ~2.3× plateau at 4–8 threads.** Diagnosed via `VEP_PROFILE`:
- The **lookup ran single-partition (serial, ~4s)** because the *session* `target_partitions` came from the old `VepConcurrencyPlan` (driven by `forks`, =1), not `threads` — so the VCF scan was never split.
- Forcing `target_partitions=threads` on an **unindexed** (sampled) VCF made DataFusion insert a **RoundRobin** repartition → output **scrambled** (parity DIFFERS). Round-robin breaks the ordered fan-in.
- On a **bgzipped + `.tbi`-indexed** VCF, the provider partitions by genomic region via `balance_partitions` (in `bio-format-core/partition_balancer.rs`), which is a **linear scan placing boundaries at cumulative byte thresholds → contiguous, position-ordered partitions**. Result: parallel lookup **byte-identical** (verified), but wall still ~2.2× because `send_wait≈15.7s` — the **single `poll_next` task is the serial spine** (drains lookup + merges one cumulative `colocated_map` + dispatches + **blocks in `AwaitingWindow`** + drains results), so lookup backs up behind it and stages don't overlap.

**REVISED target design (decided):** replace the central poll-loop with a **per-contig fan-out of `N = threads` independent range pipelines**:
- Keep the **per-contig loop**; load each contig's `SharedContigAnnotationContext` **once** (shared `Arc` across the N pipelines). *Not* DataFusion-static genome-slice partitions (those split contigs and duplicate context).
- `provider.set_vcf_filter(chrom==X)` already makes the provider yield **N contiguous, position-ordered lookup partitions for that contig** (linear-scan, reusing existing `partition_colocated_sinks`).
- Each range pipeline (independent task) consumes **its own** lookup partition stream → a **LOCAL** `colocated_map` (no cumulative map, no `Arc::make_mut` clone, no central merge) → annotates windows with the shared `Arc` context → formats. Within a pipeline, drain-its-lookup-then-annotate (or per-partition lookup barrier) so each window reads an immutable local map.
- **Ordered output** by partition id (= position order, contiguous) — concatenate pipe 0, then 1, …; no sort-merge.
- **Delete:** central fan-in, cumulative `colocated_map` + merge + clone, `AwaitingWindow` blocking. **Reuse:** context load, lookup machinery + per-partition sinks, annotation engine, VCF formatter.

**Decisions:**
- **Indexed input required for parallelism** — parallel path requires a bgzipped+`.tbi` VCF; auto-ensure/create an index, else error clearly. (Unindexed → provider gives 1 partition → would be serial.)
- **N = threads**, byte-balanced ranges (`balance_partitions` balances bytes ≈ density). Revisit oversubscription later.
- **Rely on partition-id ordering** (proven by the `balance_partitions` linear scan) + **guard the degenerate scrambling paths** (`total_bytes==0` round-robin; zero-byte regions) by falling back to serial. No defensive sort-preserving merge.
- Fixed per-contig **context load (~0.8s)** stays serial; amortized at scale, hidden across contigs by context-prefetch for whole-genome.

---

## 1. Motivation & background

### 1.1 What Ensembl VEP 115 `--fork` actually does

Analysis of the `release/115` source (`Runner.pm`, `InputBuffer.pm`, `Config.pm`, `AnnotationSource*.pm`):

- **Raw `fork()`** + `socketpair(AF_UNIX, SOCK_STREAM)` + `IO::Select`, with `Storable::freeze/thaw` IPC. *Not* `Parallel::ForkManager`. (`Runner.pm:449,483,622`)
- **Parallel unit is a sub-buffer chunk, not a contig.** The parent pulls one `InputBuffer::next()` (default `buffer_size=5000`, `Config.pm:295`), and the buffer is **hard-cut at contig boundaries** (`InputBuffer.pm:173`). It then splices the ≤5000 variant features into ~625-variant chunks (`maxForkSize = buffer/(2·forks)`) fed to ≤N forks (`Runner.pm:461-492`).
- **Hard parent-side barrier per buffer**: the next 5000-variant fill does not start until every fork of the current fill reports back. No work-stealing — a gene-dense chunk stalls the whole fill while sibling forks idle.
- **Parent is a single-threaded merge bottleneck**: every result row crosses a socket as a `Storable` blob, is `thaw`ed, and merged (output + `plugin_data` + stats + `$main::_VEP_CACHE` via `merge_hashes`). This is why VEP fork scaling plateaus; docs recommend **4 forks**.
- **Redundant cache deserialization**: each fork keeps its own `$self->{_cache}` and independently inflates the same per-region cache blocks when chunks overlap a hot region (`Cache/BaseSerialized.pm:243`). Nothing is shared laterally between siblings.
- **Order preservation is collect-then-concat**: the parent buffers each child's *entire* output in `%by_pid`, then emits in spawn order, `map {@{$by_pid{$_}}} @pids` (`Runner.pm:591`).

**Two structural weaknesses:** (1) the per-buffer barrier with no work-stealing → load imbalance; (2) the single-threaded `Storable` serialize/merge → scaling plateau.

### 1.2 Where our engine already differs

| Dimension | VEP 115 | Us (today) |
|---|---|---|
| Parallel granularity | sub-buffer chunk (~625 vars) | contig (1 contig/partition when parallelism>1) |
| Within-contig parallelism | yes (chunks), but barriered | none — serial annotation per contig |
| Shared annotation context | each fork re-deserializes regions | `SharedContigAnnotationContext` Arc-shared |
| Lookup parallelism | none | fjall/lance/indexed lookup partitions → ordered fan-in |
| Merge | single-threaded Storable | shared memory, no serialization |
| Order | collect-then-concat per buffer | streaming ordered drain by `partition_id` |

We already avoid VEP's merge bottleneck (shared memory) and we parallelize lookup. Our gap is **within-contig annotation is serial**, and the current hash-partitioned lookup shatters position locality.

---

## 2. Core model

Replace "serial annotation per contig + hash-partitioned lookup" with a **work-stealing pool of `N` window workers**. A **window** is a contiguous position range (one input buffer ≈ 5000 sorted variants, configurable, never spanning contigs). Each worker runs the **full fused per-window pipeline — lookup → annotate → format — for its range**, then emits to an ordered drain. Shared-nothing except the immutable `Arc` context.

```
                 ┌─ contig stream (position-sorted) ─┐
 WindowPlanner ──┤  cut into windows W0,W1,W2,...     │   (lazy, bounded lookahead)
                 └───────────────┬───────────────────┘
                                 │ shared work-stealing queue
        ┌────────────┬───────────┼───────────┬────────────┐
     Worker0      Worker1     Worker2  ...  Worker(N-1)     (each: lookup→annotate→format, fused)
        └────────────┴───────────┼───────────┴────────────┘
                                 │ completed windows (out of order)
                          OrderedDrain (emit by window index)
                                 │
                          VCF byte output  +  partial-stats merge at contig end
```

This directly fixes VEP's two weaknesses:
- **No per-buffer barrier** — work-stealing means idle workers grab the next window instead of waiting on a straggler.
- **No serialize/merge bottleneck** — shared-memory `Arc` context, no `Storable`.

And it preserves the locality advantage VEP lacks: each worker builds transcript geometry once over a *contiguous* range, and the lookup backend is queried over a *contiguous position range* (better Lance/fjall chunk locality).

---

## 3. Why windows are independent (correctness precondition)

For **base VEP annotation**, annotating a variant is a pure function of:
1. the variant itself (position, alleles),
2. transcripts/regulatory features overlapping its position (+ up/downstream flank distance),
3. the reference sequence for that region (FASTA),
4. colocated known variants at that position (lookup DB, keyed by position).

None depend on *other input variants* in another window — there is no left-to-right data dependency. VEP's per-buffer barrier exists because of its serial-read + collect-then-concat design, **not** a real data dependency; we do not inherit that constraint.

Apparent dependencies and why they don't block:
- **Stats/counters** — a commutative/associative reduction; per-worker partial counters merged at contig end (same semantics as VEP's `merge_hashes`, but lock-free).
- **Boundary transcripts** — handled by overlap-fetch from the shared read-only context; a worker reads any transcript overlapping its variants even if it extends past the window edge. Read-only sharing, not a compute dependency.
- **Output order** — a presentation constraint, solved by ordered drain on window index; computation order is free.

**Genuine exception (out of scope, see §9):** haplotype/phase-aware features (`--individual`, TranscriptHaplotypes) look across multiple variants on the same transcript and need per-transcript grouping, not position windows.

---

## 4. Components

- **`WindowPlanner`** — consumes the position-sorted contig stream and cuts it into windows at ~5000-variant boundaries (configurable; never spanning contigs). Produces windows lazily with bounded lookahead so memory stays `O(N + slack)` windows, not the whole contig.
- **Work-stealing queue** — a shared MPMC queue of window descriptors; workers pull as they finish. chr1 yields hundreds of windows, so stealing stays well-balanced.
- **`WindowWorker` (×N)** — owns:
  - a **worker-local geometry/cDNA cache** `HashMap<TranscriptId, CdnaCoordTable>` (see §5),
  - **partial stats counters**.
  For each window: query the lookup backend (fjall/Lance/indexed-parquet) over the window's contiguous position range; annotate each variant against overlapping transcripts from the shared context; format to VCF bytes. **All three steps are one fused task** (§6).
- **`OrderedDrain`** — emits window `k` once it is complete and all `< k` are emitted; holds out-of-order completions in a small reorder map bounded by max-in-flight. Under work-stealing, completion order ≈ window order, so the reorder buffer stays small.
- **`StatsMerge`** — merges per-worker partial counters at contig end.
- **`SharedContigAnnotationContext`** (exists, `annotate_provider.rs:9257`) — stays `Arc`-shared and **immutable**: transcript structure, sequence refs, regulatory features. The derived geometry cache moves *out* of it (§5).

---

## 5. Geometry/cDNA cache — no change needed (already concurrency-safe)

Commit `8a7d4ca` already implemented the cDNA-coord cache as a **`OnceLock`-backed** memo on `TranscriptFeature`:

```rust
// transcript_consequence.rs:143
#[derive(Debug, Default)]
pub struct GeometryCache(std::sync::OnceLock<Option<Vec<TranscriptCdnaCoord>>>);
```

`OnceLock` makes `TranscriptFeature` **`Sync`**, and the memoization site (`transcript_consequence.rs:7545`, `get_or_init`) already builds each transcript's coord table **exactly once and shares it lock-free** across concurrent workers. The transcript already lives in `base_transcripts: Arc<Vec<TranscriptFeature>>` in the shared context.

**Decision: keep it as-is.** No refactor. The shared `OnceLock` is strictly better than a worker-local table on this CPU-bound engine:
- The build (`transcript_cdna_coords_uncached`) is the geometry recomputation flagged as a top CPU cost in the perf handoff (issue C). Shared = **one build per transcript per contig, ever**; worker-local would rebuild it per worker (and 2–3× for boundary-spanning transcripts).
- The read path after init is a **lock-free atomic load** — not a per-access lock. The only synchronization is a one-time init; a second worker hitting the same uninitialized transcript blocks only for that single build, then reads lock-free.
- With position-window locality a transcript overlaps only ~1–3 windows, so init contention is rare and brief regardless.

(Orthogonal future optimization, out of scope: the reader does `.get_or_init(...).clone()` — a `Vec` clone per read, `transcript_consequence.rs:7550` — which could become a borrow.)

---

## 6. Fused single code path (all thread counts)

Lookup → annotate → format is **one fused task per window at every thread count, including `threads=1`**. There is no separate serial path vs. parallel path — `threads=1` is simply one worker pulling windows sequentially and running the same fused task. This gives **one code path**, so the parity gate (§8) exercises identical logic at N=1 and N>1, dramatically lowering parity risk.

**Accepted trade-off:** fusing removes the lookup/annotation IO-overlap that today's separate spawned lookup task provides.
- At `threads ≥ 2`, overlap is recovered for free: with more in-flight windows than workers, while one worker blocks on lookup IO another does CPU annotation.
- At `threads = 1`, the fused worker stalls on lookup IO then does CPU — no overlap. This is accepted: the engine is CPU-bound with a warm cache (per the perf handoff), and `threads=1` is the degenerate/debug/fallback case, not a production speed mode.

---

## 7. Thread budget (unifies plan `2026-05-22`)

Because a window owns its full pipeline, the budget is a **single knob `threads = N`**:
- `threads=1` → one worker, fused task, no queue contention / no spawn fan-out (still the fused code path).
- `threads=N` → N window workers; lookup and format live *inside* each worker, so there is **no separate lookup/format pool to size**.

`threads` is the **only** parallelism knob the engine accepts in its options JSON. It **replaces** (clean break — see §7.1) the cluster the parser reads today at `annotate_provider.rs:5487`: `worker_forks`, `target_partitions`, `chromosome_lanes`, `inline_lookup`.

**Across contigs (inter-chromosome model):** process **one contig at a time with all N workers on its windows**, with **context prefetch** (plan `2026-05-23`) hiding the next contig's immutable-context load behind the current contig's tail. The worker pool persists across contigs; the queue is refilled per contig.

This **replaces today's inter-chromosome parallelism** (`partition_contigs_for_execution`, `annotate_provider.rs:8881`, which today runs N contigs concurrently, each serial inside). The choice is deliberate — it maximizes per-contig throughput (the stated goal) and keeps only ~1–2 `Arc` contexts resident. **Accepted trade-offs:**
- **Tiny contigs underutilize cores.** A contig with fewer windows than `N` (e.g. chrM) cannot fill the worker pool; the spare workers idle until the next contig is admitted.
- **No whole-genome long-tail smoothing.** Contigs are processed serially, so total wall-clock is the sum of per-contig (parallelized) times plus per-contig setup gaps, rather than a single globally load-balanced sweep.

These are acceptable because the optimization target is per-contig throughput, not whole-genome wall-clock. See §10 for the global-queue alternative if that target changes.

### 7.1 Downstream integration (vepyr / polars-bio) — clean break

The new model collapses the engine's option surface to a **single `threads: usize`** knob. This is a **clean break**: the engine no longer accepts the old cluster, and the downstream consumer is updated in lockstep.

**Engine side** (`annotate_provider.rs:5487`):
- Parse only `threads` from the options JSON. Default `threads=1` (fused single-worker path).
- **Remove** parsing of `worker_forks`, `target_partitions`, `chromosome_lanes`, `inline_lookup`.
- **Fail loud, not silent**: if the options JSON still contains any removed key, return a clear error ("`worker_forks` is removed; use `threads`") rather than ignoring it and silently running serial. This catches stale callers immediately.

**Downstream side** (vepyr / polars-bio Python wrapper):
- Expose a single user-facing parameter `threads=N` (Pythonic, distinct from VEP's process-based `--fork`).
- Emit only `threads` in the options JSON handed to the engine; stop emitting the four removed keys.
- Window size, lookahead/`slack`, and reorder-buffer bounds are **internal engine defaults**, not user-facing (clean break, no advanced overrides in v1).

**Cross-repo rollout** (this repo is pinned by commit `rev` in polars-bio, per CLAUDE.md):
1. Land the engine change here (new `threads` parser + removed keys + fail-loud).
2. Update **all in-repo callers** in the same change: example benches (`examples/bench_annotate_vcf.rs`, `bench_annotate_debug.rs`, the `bench_*` partition/fjall examples) and any scripts (`scripts/run_annotation_fast.py` and siblings) that pass the old options.
3. Bump the pinned `rev` in vepyr / polars-bio and switch the Python knob to `threads=N` atomically with that bump.

Because it is a clean break, the engine version bump and the downstream `rev` bump must ship together — there is no transition window where mixed versions interoperate.

---

## 8. Ordering, backpressure, memory

- **Output order**: by window index in `OrderedDrain`; computation order is free.
- **Backpressure**: `WindowPlanner` blocks when in-flight windows reach `N + slack` (precedent: existing `LOOKUP_PARTITION_QUEUE_BATCHES = 2`).
- **Memory**: `O(N + slack)` windows of variants + N worker-local geometry caches + a small reorder buffer. No whole-contig materialization.

---

## 9. Testing

- **Parity gate**: consequence-identical (ideally byte-identical) output vs. the current serial path on a fixed chr1 slice, across `threads ∈ {1, 2, 4, 8}` — proves window independence and ordering. Single fused code path means N=1 is the reference.
- **Determinism**: identical output regardless of worker count or completion order.
- **Stats correctness**: merged counters equal serial counters.
- **Boundary correctness**: variants near window edges whose transcripts straddle the boundary annotate identically to serial.
- **Scaling bench**: per-contig wall-clock (chr1) at `threads ∈ {1,2,4,8,16}` → scaling curve, compared to the engine-bound baseline in the perf handoff (`2026-06-17-vep-engine-perf-handoff.md`).
- **Option-surface (clean break, §7.1)**: engine accepts `threads` and returns a clear error on any removed key (`worker_forks`/`target_partitions`/`chromosome_lanes`/`inline_lookup`); downstream emits only `threads`. End-to-end smoke through vepyr/polars-bio with `threads=4`.

---

## 10. Out of scope (v1)

- **Haplotype/phase-aware features** (`--individual`, TranscriptHaplotypes): need per-transcript grouping, not position windows. Unsupported under parallel mode in v1 — fall back to `threads=1` when enabled, or revisit later.
- **Global cross-contig window queue** (distributing one contig's spare cores onto another contig): deferred. This is the natural generalization — a single work-stealing queue of windows spanning *all* contigs, with `N` total workers pulling from any admitted contig. It would eliminate both accepted trade-offs in §7 (tiny-contig underutilization and the whole-genome long tail: a worker that finishes chrM steals chr1's remaining windows), at the cost of ~2–3 resident `Arc` contexts instead of 1–2. **Adopt this if the optimization target shifts from per-contig throughput to whole-genome wall-clock.** Bounded memory is preserved by admitting contigs in order and evicting a contig's context once all its windows drain (prefetch ~1 ahead).
- **Density-balanced static partitioning**: not needed — dynamic work-stealing handles uneven gene density without a cost estimate.

---

## 11. Relationship to existing plans

This design **unifies and refines** three not-yet-started plans:
- `2026-05-21-parallel-annotation-partitions.md` — parallel annotation workers + ordered drain (now: position-window workers).
- `2026-05-22-vep-threads-concurrency-budget.md` — the `threads` knob (now: single knob, fused pipeline, clean-break downstream surface — §7.1).
- `2026-05-23-contig-context-prefetch.md` — hide next-contig context load behind the current contig's tail.

It supersedes the current branch's hash-partitioned lookup fan-in with range-windowed, worker-owned lookup.

---

## 12. How this beats VEP 115 (summary)

| VEP 115 weakness | This design |
|---|---|
| Per-buffer barrier, no work-stealing → idle forks | Dynamic work-stealing queue, no barrier |
| Single-threaded `Storable` serialize/merge | Shared-memory `Arc` context, lock-free stats reduce |
| Redundant per-fork cache deserialization | Shared immutable context; geometry built once per transcript per contig (shared `OnceLock`), not per worker |
| Sub-buffer granularity, not contig-parallel | True within-contig parallelism over position windows |
| Practical 4-fork plateau | Scales with cores until per-contig window supply / memory bound |
