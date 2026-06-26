# Bounded-Overlap Stateful Parallel Annotation Design

**Date:** 2026-06-25
**Status:** Design proposal (pre-implementation)
**Repo:** `datafusion-bio-functions` (`datafusion/bio-function-vep`)
**Supersedes / corrects:** the "conservative workers=1 fallback" framing in `2026-06-24-stateful-hgnc-mismatch-root-cause-and-performance.md`
**Related:** `2026-06-17-within-contig-parallel-annotation-design.md`, `2026-06-18-sharded-vcf-output-design.md`

## 1. Problem

The fused parallel annotation model splits a contig into N contiguous position
partitions, one annotation worker each (`StreamState::AnnotatingParallel`, per-worker
`AnnotationWorkerState`). For Merged/RefSeq caches this produced an `HGNC_ID`
regression: the `fast_chr1_22_merged_summary_20260623_0712.md` e2e report shows
**9,792 `HGNC_ID` mismatches** (everything else at 100%), clustered on gene-dense
chromosomes (**chr4 = 8,281**, chr10 = 600, chr16 = 307, chr7 = 270, ...), all of the
false-positive shape `vepyr = HGNC:…, VEP = (empty)`.

`master` is clean because it never runs Merged/RefSeq annotation from partition-local
state. The regression lives entirely in the parallel boundary model, not in HGNC
parsing, the consequence engine, Lance lookup, or VCF comparison.

## 2. Root cause

`HGNC_ID` is **not** a pure function of a row and its local transcripts. Ensembl VEP's
`merge_features()` copies HGNC metadata between transcripts that co-reside in the same
**input buffer** (default `buffer_size = 5000` parsed VariationFeatures), mutates the
transcript object in place, and keeps it in a per-region feature cache that persists
across adjacent buffers while the ~1 Mb region stays resident. So a row's `HGNC_ID`
depends on (a) which transcripts shared its buffer and (b) the accumulated left-to-right
mutation history while the region was active.

Two boundary hazards arise when this is parallelised:

- **Hazard 1 — buffer membership.** Which transcripts are "in the buffer" is defined by
  the buffer's variant `[min..max]` range. A seam that splits a buffer gives both
  adjacent workers a *partial* buffer with a different donor set than VEP's whole buffer.
- **Hazard 2 — carried state.** A donated HGNC persists into later buffers via the region
  cache. A worker that starts fresh at a seam loses mutations VEP carried across it.

The current code diverges on the **donor scope** axis. Recipients use the buffer range
(`buffer_variant_bounds`, `annotate_provider.rs:10217`) but donors are drawn **region-wide**
(`region_donors`, `annotate_provider.rs:10623-10631` — the full active ~1 Mb cache-region
transcript set). That region-wide scope was introduced to fix the *opposite* failure (a
false negative where a partition-local donor outside a worker's range was missed), but it
overshoots: a recipient receives HGNC from a same-gene donor that is active in the 1 Mb
region but **not** co-resident in VEP's 5,000-variant buffer. VEP leaves such a recipient
empty; the branch fills it. That is the chr4 TENM3 `HGNC:29944` over-donation — 8,281 of
the 9,792 mismatches.

The two failure shapes are a single defect seen from both sides: donor scope too narrow →
false negative; donor scope too wide → false positive. The correct scope is **neither** —
it is VEP's global 5,000-variant input-buffer range.

### 2.1 Validated assumptions

- **Per-contig buffering.** `AnnotationWorkerState::new` (`annotate_provider.rs:9793`)
  resets `input_buffer_accumulator` and `next_input_buffer_id = 0`, and a worker is created
  fresh at every `StartContig`. So buffers always restart at each contig's first variant;
  VEP's cross-contig buffer straddle (caveat 3 in the 06-24 note) cannot occur in our
  model. The report corroborates this: 11 chromosomes are at exactly 0 `HGNC_ID`
  mismatches — a global buffer-offset bug would smear noise across all of them.
- **Donation does not chain.** `collect_hgnc_donors` reads `gene_hgnc_id_native`;
  `apply_hgnc_donors` writes the effective `gene_hgnc_id`. A received value is never
  re-donated, so a donor→recipient hop is single, not transitive. This bounds reach.

## 3. Corrected invariant

> For Merged/RefSeq, parallel work may not change VEP's global input-buffer boundaries or
> the carried transcript-mutation state visible at each buffer. Correctness comes from
> reconstructing VEP's global buffer boundary and per-buffer state, never from collapsing
> to one worker.

Concretely:

1. Seams must fall on VEP global input-buffer boundaries (multiples of 5,000 variants per
   contig), so no buffer straddles a seam (Hazard 1).
2. Donor scope per buffer must be the buffer's `[min..max]` range, not the 1 Mb region.
3. The carried state at a seam must equal serial VEP's state at that exact boundary
   (Hazard 2).

## 4. Options considered

All three keep the fused N-worker model and satisfy invariant points 1–2; they differ only
in how a worker obtains the carried state at its seam (point 3).

- **Option Q — snap seams to quiescent gaps.** Place seams only where no region is active
  across them, so a fresh worker's empty state already equals VEP's. No extra machinery,
  but seam placement is constrained → load skew on gene-dense chromosomes, and the usable
  worker count is capped by the number of quiescent gaps. Rejected: balance cost falls
  exactly on the busiest chromosomes.
- **Option S — boundary-state snapshots.** A worker tells the next worker its end-of-range
  state at the seam. Allows arbitrary equal-count seams, zero recompute, but requires
  serialising the full carried state (mutated transcript overlays + active region
  lifetimes) and **proving completeness** — a missed field is a rare silent mismatch.
  Held in reserve as a runtime-cost optimisation if ever needed.
- **Option O — bounded overlap (chosen).** Each worker warm-starts a bounded distance
  before its seam, in state-only mode, rebuilding the carried state by re-running the same
  serial donation code over real variants, then discards the warm-up output and emits its
  range. Equal-count seams, no pre-pass, no synchronisation, no state serialisation;
  correctness is by construction (it reuses the serial path).

## 5. Chosen design: bounded overlap

### 5.1 Partitioning

Per contig: `B = ceil(variants / 5000)` whole VEP buffers. Give each of N workers ~`B/N`
consecutive buffers; every seam lands on a buffer boundary (a multiple of 5,000 variants
from the contig's first variant). Because buffer boundaries are dense, equal-whole-buffer
slices are within one buffer of perfectly balanced — preserving today's balance.

The equal-count split already counts variants, so the cumulative variant rank at each seam
(a multiple of 5,000) is known as a byproduct; each worker is told its start rank.

### 5.2 Per-worker warm-up

Worker `k` (covering buffers `[s_k, s_{k+1})`) processes from `overlap_start_k` instead of
its seam:

```
overlap_start_k = seam_position_k − overlap_width
```

It runs **state-only** prepare over the warm-up region (build buffer-local transcripts +
HGNC donation + persist; i.e. `prepare_buffer_annotation_context` /
`build_stateful_buffer_local_transcripts_cow`), producing no output, then continues into
its real range with full annotation. The warm-up:

- reads variant **positions** only — it does **not** issue the Lance variation/colocated
  lookup (carried HGNC depends only on transcripts);
- skips the consequence engine, HGVS, SIFT/PolyPhen, and VCF formatting;
- discards all rows before `seam_position_k` (worker `k-1` owns and emits them).

Worker 0 has no predecessor and starts at the contig's first variant with empty state (the
true VEP start), so it needs no warm-up.

### 5.3 Donor scope correction

Independently of overlap, revert the donor source from region-wide (`region_donors`,
`annotate_provider.rs:10623-10631`) to the buffer `[min..max]` range, matching VEP's
"filter features to the buffer's min/max, then `merge_features()`". With seams on buffer
boundaries and the warm-up supplying pre-seam co-residence, buffer-range donors no longer
miss legitimate donors (the original false negative) and no longer over-donate (the chr4
false positive).

### 5.4 Overlap width

`overlap_width = max(transcript span on the contig) + one cache-region width`.

- It is a true geometric **upper bound**: a donor cannot influence a seam from farther than
  the recipient's region stays active, and a region cannot outlast its transcript's span
  (plus the multi-region `start-region` reset enforced by
  `reset_persisted_hgnc_effective_values_outside_start_region`). Donation does not chain
  (§2.1), so there is no transitive blow-up.
- It is computed by a `max()` reduction over the already-loaded `base_transcripts` at
  context-load time — **zero extra IO**, and identical in every worker (same shared data),
  so it needs **no synchronisation** and cannot under-read.
- The conservative per-contig `max(span)` is sufficient; per-seam tightening is explicitly
  **out of scope** (YAGNI) given the measured cost below.

### 5.5 Measured cost

`max(transcript span)` per contig from `transcript.lance` (115 GRCh38 merged), translated
to redundant work via the report's per-chromosome variant density:

| metric | worst contig | median contig | smallest |
|---|---|---|---|
| max span | 2.47 Mb (chr16) | ~1.5 Mb | 0.49 Mb (chr19) |
| redundant warm-up | **0.68 buffer** (chr7/16) | **0.49 buffer** | 0.15 buffer (chr19) |

The conservative warm-up is **under one 5,000-variant buffer on every chromosome**, and it
is the *cheap* state-only stage. For N=8 on chr1: 7 seams × ~0.4 buffer ≈ 2.8 buffers of
state-prepare against ~65 total buffers — under ~1% of wall time. Redundancy is negligible.

### 5.6 Implementation sequencing

The donor-scope revert (§5.3) and the overlap (§5.2) are jointly necessary for workers>1
but fix different halves, so they are implemented and validated as two gated steps within
one change. They are coupled — neither alone makes workers>1 correct (overlap alone still
over-donates with region-wide scope; donor revert alone reintroduces the carried-state
false negative at seams) — so the steps land together, or workers>1 stays gated between
them.

1. **Step 1 — donor-scope revert** (region-wide → buffer-range). Independently confirms the
   false-positive root cause. **Gate: chr4 + chr2 merged Lance at `workers=1` → 0 `HGNC_ID`
   mismatches.** At `workers=1` there are no seams, so buffer-range donors equal VEP's
   buffer exactly; the carried-state hazard does not arise. This is the cheapest,
   highest-confidence validation and is run first against a captured `workers=1` baseline
   on the same build.
2. **Step 2 — bounded overlap** (§5.2/§5.4). Adds the warm-up so carried state crosses
   seams. **Gate: chr4 + chr2 merged Lance at `workers={1,4}` → 0 mismatches, byte-identical
   bodies**, plus the warm-up-state debug assertion.

## 6. Correctness argument

- **Hazard 1** is removed by seams on buffer boundaries: no buffer straddles a seam, so
  every worker's locally-cut buffers equal VEP's global buffers within its range.
- **Hazard 2** is removed by the warm-up: starting empty at `seam − overlap_width` and
  re-running the real donation code, any wrong state from the empty start is pruned before
  the seam (regions active at the start but not reaching the seam are evicted within
  `overlap_width`), and every region active *at* the seam was seeded within the warm-up, so
  its carried state is reconstructed exactly as serial VEP had it.
- **Donor scope** matches VEP exactly (buffer `[min..max]`), removing both the false
  positive and the false negative.
- Output ordering is unchanged: workers emit contiguous buffer ranges in order; warm-up
  rows are dropped, never emitted.

## 7. Affected code (feasibility)

- `annotate_provider.rs`
  - `AnnotationWorkerState` / `ParallelContigState` / `StreamState::AnnotatingParallel`:
    add `warm_up_until` (seam position) and a state-only prepare loop that consumes the
    warm-up region before the first emitted buffer.
  - `prepare_buffer_annotation_context` / `build_stateful_buffer_local_transcripts_cow`:
    reused as-is for warm-up (already the state-transition unit).
  - `collect_hgnc_donors` call site (`region_donors`, ~`:10623`): change donor source to
    buffer-range transcripts.
  - context load: compute `overlap_width` from `base_transcripts` (one `max()`).
- partition/seam selection (lookup partitioning + `vcf_sink.rs`): choose seams on 5,000-
  multiple buffer boundaries and pass each worker its start rank and `warm_up_until`.
- `window_planner.rs`: warm-up windows are state-only (no lookup, no output).

No new dependency (`rayon` etc.); no change to the sharded VCF output contract.

## 8. Acceptance gates

A change does not land unless all pass:

- `cargo test -p datafusion-bio-function-vep --lib`.
- New unit tests asserting `workers ∈ {1,4}` produce **identical** `HGNC_ID` for:
  - chr4 TENM3 (false-positive class: donor outside the buffer must not donate);
  - chr2 PDK1/NEMP2 (false-negative class: same-buffer donor across a seam must donate);
  - state carry across two adjacent buffers + region eviction.
- A debug assertion (test builds) that a warm-up worker's reconstructed state at its seam
  equals the serial state at that boundary.
- E2E chr4 + chr2 merged Lance at `workers=1` **and** `workers=4`:
  `field_mismatch_counts={}`, byte-identical bodies.
- One larger run (chr1-22 merged Lance) at `workers=4` reproducing the previous report's
  mismatch class at 0.
- Performance report: annotation wall time and variants/s for `workers {1,2,4,8}`,
  confirming recovery vs the ordered single-worker fallback.
- Fail-closed: if `overlap_width` cannot be computed (no transcript cache) or seams cannot
  be aligned to buffer boundaries (non-indexed input under `workers>1`), keep the existing
  `workers>1` guard error rather than silently using partition-local bounds.

## 9. Out of scope

- Per-seam overlap tightening (static per-contig bound is already <1 buffer).
- Boundary-state snapshots (Option S) — revisit only if a future workload makes warm-up
  redundancy material.
- Ensembl-only (stateless) cache sources, which may keep the simpler independent-shard path.

## 10. Implementation status (2026-06-26)

Diffing `e12e647` (fast parallel, HGNC-wrong) → HEAD `2ca37e6` shows the three
correctness commits (`eb308f1`, `adfd2ba`, `2ca37e6`) **already landed two of this
design's three invariant fixes**, both parallel-safe and committed:

| Invariant (§3) | `e12e647` (before) | HEAD (now) | Parallel-safe? |
|---|---|---|---|
| 1. Seams on 5,000-unit buffer grid | whole-batch (8192) row cuts (`drain(..window_end)`) | `drain_window_input_units` cuts exact input-buffer units, slicing the boundary batch | ✅ committed |
| 2. Donor scope = buffer `[min..max]` | `region_donors` (full 1 Mb active region) | `apply_buffer_local_hgnc_propagation` (buffer min/max) | ✅ committed |
| 3. Carried state at seam | none (each window starts empty) | sequential `persisted_seed` threaded window N→N+1 | ❌ **serial-only** |

So the design's **Step 1 (donor-scope revert, §5.3) and the grid-alignment are DONE.**
The *only* reason workers>1 is gated to serial today is invariant 3: the carry is
threaded sequentially (`annotate_window_owned`'s `persisted_seed` in/out, dispatch loop
`:11551`/`:11617`), which is correct only with one window in flight — hence
`annotation_workers = 1` for Merged/RefSeq (`:5276`, the perf regression).

### Remaining work = Step 2 (bounded overlap) + drop the serial gate

The before-mechanism we are restoring is the `spawn_blocking` **whole-window** parallel
dispatch (`ContigAnnotationStream` poll loop, `:11540`–`:11562`). It already ran N windows
concurrently on the tokio blocking pool, so the existing demand-fetch SIFT path
(`get_position_predictions` → `block_in_place(Handle::current().block_on(take_keys))`,
`:3166`) **works unchanged on those threads — no pre-warm, no `std::thread` lanes, no SIFT
machinery is needed.** (This is why Option A's "SIFT blocker" does not arise here: it was an
artifact of Option A switching to raw `std::thread::scope` lanes that lose the runtime
context.)

Concrete integration points:

1. **Remove the serial gate** (`:5266`–`:5276`): stop forcing `annotation_workers = 1` for
   Merged/RefSeq; let it equal `requested_workers`, fail-closed per §8 if `overlap_width`
   is uncomputable.
2. **Replace the sequential carry with self-reconstructed carry.** Drop the
   `persisted_seed` in/out threading between windows. Instead each window task warm-starts:
   given its real range `[seam_k, seam_{k+1})`, it first runs **state-only** prepare
   (`prepare_buffer_annotation_context` / `build_stateful_buffer_local_transcripts_cow`,
   donation + persist, output discarded) over the whole buffers covering
   `[seam_k − overlap_width, seam_k)`, then annotates its real range. Worker 0 (contig
   start) skips warm-up.
3. **Feed each window its warm-up prefix.** The dispatch loop drains windows sequentially
   from `window_buffer`; retain a sliding tail of the last `overlap_width` input units (real
   variant positions only — no Lance lookup/colocated needed for warm-up) and pass it,
   read-only, ahead of window k's batches. Warm-up rows are dropped before emission.
4. **Compute `overlap_width`** once at context-load via a `max()` over `base_transcripts`
   spans (§5.4) and store it on the shared/config struct.
5. **Debug assertion (§8):** in test builds, assert a warm-up window's reconstructed
   persisted state at its seam equals the serial state at that boundary.

Steps 2+3 land together with step 1; gates per §8 (chr4 + chr2 merged Lance workers∈{1,4}
→ 0 `HGNC_ID` mismatches, byte-identical, on 3 consecutive runs to catch intermittent
truncation), then the chr1-22 run and the workers {1,2,4,8} perf report.
