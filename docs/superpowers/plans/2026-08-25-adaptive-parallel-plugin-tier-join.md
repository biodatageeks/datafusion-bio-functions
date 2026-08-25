# Adaptive Parallel Plugin Tier Join — Implementation Plan

**Date:** 2026-08-25
**Status:** Full design retained for reference; the user selected the minimal
custom implementation on 2026-08-25.
**Root-cause evidence:** `docs/superpowers/handovers/2026-08-24-cadd-tier-join-reservation-leak-RESOLVED.md` and `.planning/debug/cadd-tier-join-algorithm-selection.md`

## Minimal implementation selected

The implemented first slice intentionally avoids the separate statistics
catalog and context-factory layers proposed later in this document. It:

- asks DataFusion to build an optimized CollectLeft HashJoin candidate;
- reads row and byte statistics from the physical build child DataFusion
  actually chose;
- compares DataFusion 53's estimated input/hash-table/bitmap reservation with
  the active pool's currently available bytes;
- executes that exact HashJoin plan when it fits, otherwise flips only
  `prefer_hash_join` and verifies DataFusion emitted SortMergeJoin;
- applies the decision independently to probe construction and the final tier
  join, using runtime-default parallelism with a minimum of two partitions;
- retries once with SMJ only for a `ResourcesExhausted` root error that both
  names and is traced to `HashJoinInput`; and
- always requests explicit final ordering for parallel output.

Unknown statistics select SMJ. Inexact statistics remain advisory and are
protected by the attributed runtime fallback. The richer projected-Parquet
statistics model and real CADD benchmark matrix below remain possible follow-up
work, not claims about this minimal slice.

## Goal

Make both joins in the per-chromosome plugin-cache tiering pipeline run with
DataFusion parallelism and automatically select the fastest strategy that is
safe for the active memory pool:

- use HashJoin when its non-spillable build side is predicted to fit;
- use SortMergeJoin when the hash build cannot be shown to fit;
- retry once with SortMergeJoin if an estimated-safe HashJoin is rejected by
  the memory pool at runtime;
- base every decision on input statistics, the runtime memory limit, and the
  runtime's default parallelism — never on chromosome names or chromosome-row
  thresholds.

"Best" has a precise meaning in this plan: prefer the linear-time HashJoin when
it satisfies the memory proof; otherwise use the bounded/spillable SMJ. The real
chr21 comparison supports that policy: HashJoin took 118.2 s, while SMJ produced
identical output under the same 8 GiB pool in 123.7 s (+4.7%).

## Non-goals

- Do not bump DataFusion. DataFusion 55 fixes the shared-buffer-accounting bug,
  but does not add the adaptive HashJoin-vs-SMJ decision required here.
- Do not declare either input sorted. Neither current input is globally ordered
  on `(chrom, start, allele_string)`.
- Do not select strategies by plugin name, chromosome, file name, or fixed row
  count.
- Do not run multiple chromosomes concurrently in the first implementation.
  Parallelism is within each chromosome. Cross-chromosome concurrency needs a
  separate global memory-permit scheduler or it can multiply the memory budget.
- Do not use Partitioned HashJoin as the memory fallback. It parallelizes work
  but remains non-spillable and does not reduce total live build bytes.

## Required architectural change

The current `tier_join_config()` applies one session configuration to probe
construction, aggregation, the final join, and the output sort. With
`target_partitions=1`, DataFusion 53 always creates CollectLeft HashJoin;
`prefer_hash_join=false` and the collection thresholds cannot change that.

Replace that single context with two independently planned stages:

```text
deduped plugin RecordBatches
        |
        | exact projected batch statistics
        v
Adaptive decision A: build variation probe
        |-- HashJoin context (parallel, bounded pool)
        `-- SortMergeJoin context (parallel, bounded/spillable)
        |
        v
owned, execution-sized probe RecordBatches
        |
        | exact plugin + probe statistics
        v
Adaptive decision B: final LEFT tier join
        |-- HashJoin context (parallel, bounded pool)
        `-- SortMergeJoin context (parallel, bounded/spillable)
        |
        v
explicit ORDER BY start -> warm/cold writers
```

The two decisions are mandatory. Adapting only the final join leaves the probe
construction HashJoin able to collect the raw variation shard: 0.95 GiB on
chr21, projected to about 9.44 GiB on chr1 and 10.04 GiB on chr2.

## Strategy model

Create `datafusion/bio-function-vep/src/plugin_cache/join_strategy.rs` with a
pure, deterministic selector. It must not read environment variables, inspect
chromosome labels, or create a DataFusion context.

```rust
enum JoinAlgorithm {
    Hash,
    SortMerge,
}

enum StatsConfidence {
    Exact,
    Estimated,
    Unknown,
}

struct JoinInputStats {
    rows: Option<usize>,
    arrow_bytes: Option<usize>,
    largest_batch_bytes: Option<usize>,
    confidence: StatsConfidence,
}

struct JoinResources {
    memory_limit_bytes: Option<usize>,
    memory_reserved_bytes: usize,
    target_partitions: usize,
    execution_batch_rows: usize,
}

struct JoinDecision {
    algorithm: JoinAlgorithm,
    target_partitions: usize,
    estimated_hash_build_bytes: Option<usize>,
    available_bytes: Option<usize>,
    reason: JoinDecisionReason,
}
```

### Statistics inputs

1. **Deduped plugin MemTable:** compute exact rows and Arrow bytes from the
   projected join-key batches with DataFusion/Arrow's own memory-size functions.
   Do not use the plugin's output row width for the probe-construction join.
   Materialize `plugin_variation_keys` before probe construction so the
   physical join sees those exact statistics too. When there are no match
   columns, the existing keep-first dedup already makes the projected variation
   key unique. When match columns exist, run `SELECT DISTINCT chrom, start,
   allele_string` under the finite spill pool, collect it, and rechunk it into
   owned execution-sized batches.
2. **Raw variation Parquet:** obtain exact row count from Parquet metadata and
   projected byte statistics for `chrom`, `start`, `allele_string`, and `tier`.
   If projected Arrow bytes are not exact, scan an execution batch to derive a
   per-row Arrow-width estimate and mark it `Estimated`.
3. **Materialized variation probe:** compute exact rows and bytes after
   `rechunk_probe_owned`; these batches already own their buffers.
4. **Runtime resources:** read `MemoryPool::memory_limit()` and `reserved()`.
   Use `SessionConfig::new().target_partitions()` for machine-default
   parallelism. Hash may cap the count by non-empty probe batches/scan
   partitions. SMJ must use `max(2, runtime_default)` logical partitions because
   DataFusion 53 otherwise silently emits CollectLeft HashJoin; a one-core
   runtime still gets a correct SMJ plan even though the scheduler cannot run
   its two tasks simultaneously. Do not introduce a row-count threshold.

### Hash-build estimate

Hash mode deliberately uses one shared CollectLeft build plus parallel probing,
not Partitioned HashJoin. The materialized inputs have exact physical
statistics; raw Parquet may still have only a sampled Arrow-width estimate. Set
both Auto-to-CollectLeft thresholds to `usize::MAX`, let
JoinSelection swap the smaller side, create the optimized physical plan, and
inspect it before execution. The plan is executable as Hash only if it contains
`HashJoinExec(mode=CollectLeft)` and its physical left/build child can be
identified as one of the measured inputs. Recompute the hash estimate for that
actual child (using exact MemTable bytes or the Parquet Arrow-width estimate)
and execute only if that actual build fits. Unknown build identity/statistics or
an actual build that does not fit is a pre-execution reason to replan as SMJ.

Physical validation is confidence-aware:

- exact MemTable input: require relation/schema identity plus exact row and byte
  equality with the measured batches;
- sampled Parquet input: require relation/schema identity and exact metadata row
  count, retain the external Arrow-byte estimate and `Estimated` confidence,
  and do not compare that estimate with Parquet file-byte statistics.

Calculate the estimate from that confirmed physical build child using the same
components reserved by DataFusion 53's `collect_left_input`:

```text
input RecordBatch memory
+ datafusion_common::utils::memory::estimate_memory_size<(u32, u64)>(
    rows,
    size_of::<JoinHashMapU32>(),
  )
  (or the u64-index variant above u32::MAX rows)
+ visited-index bitmap for outer joins
+ one data-derived execution-batch workspace per active partition
```

All arithmetic must be checked/saturating. The selector compares the smaller
valid candidate estimate with:

```text
available = finite_pool_limit - currently_reserved
```

Decision rules, in order:

1. Empty input: HashJoin.
2. Infinite pool with usable statistics: HashJoin.
3. Unknown pool limit, missing row count, or missing byte estimate:
   SortMergeJoin — an unspillable build cannot be proven safe.
4. Estimated hash build plus streaming workspace fits available memory:
   HashJoin.
5. Otherwise: SortMergeJoin.

No percentage or per-chromosome safety threshold is needed. Runtime fallback is
the guard against sampling error, Parquet-vs-Arrow size differences, and future
DataFusion hash-table accounting changes.

`JoinHashMapU32` and `JoinHashMapU64` are available from DataFusion 53's public
`physical_plan::joins::join_hash_map` module. Using their actual `size_of`
matches DataFusion's fixed map-structure overhead without a copied numeric
constant.

### DataFusion configuration

Build a fresh context for each selected stage over the same finite
`FairSpillPool`:

| Decision | Required configuration | Expected physical join |
|---|---|---|
| Hash | runtime `target_partitions > 1`, `prefer_hash_join=true`, Auto collection thresholds set to `usize::MAX` | verified `HashJoinExec(mode=CollectLeft)` with one shared smaller build and parallel probing |
| SortMerge | `target_partitions=max(2, runtime_default)`, `prefer_hash_join=false` | verified `SortMergeJoinExec` with required repartition/sorts |

Keep `schema_force_view_types=false` in both configurations; it prevents the
already-proven shared Utf8View accounting problem. Remove the ineffective zero
collection-threshold settings.

## Parallel data layout and output ordering

- Register plugin and materialized-probe MemTables with multiple partitions.
  Assign complete owned batches to partitions by accumulated batch bytes, using
  the selected target partition count. Do not slice batches or share newly
  created parent buffers.
- Let DataFusion repartition Parquet scans and join inputs using its runtime
  target partition count.
- Parallel execution invalidates the current implicit-order fast path. Always
  use the explicit `ORDER BY p.start` stream for parallel contexts. The sort is
  spillable and preserves the writer's `(tier, start)` invariant after the
  stream is split into warm and cold runs.
- Remove the fast-attempt/order-violation/retry double execution once the
  parallel sorted path is verified. Keep `assert_start_monotonic` as a final
  corruption guard, not as strategy control flow.

## Runtime fallback

Preflight statistics are advisory; the memory pool remains authoritative.

For each stage selected as HashJoin:

1. Execute it against the finite pool.
2. Have the stage-local pool wrapper record the stable consumer label whenever
   its inner pool rejects `try_grow`. If `DataFusionError::find_root()` is
   `ResourcesExhausted`, retry only when that recorded label is the stage's
   `HashJoinInput` reservation. Aggregate, sorter, repartition, disk, and other
   resource failures propagate unchanged.
3. Re-register the same immutable owned input batches in the new context.
4. Propagate all non-memory errors unchanged.
5. Never fall back from SMJ to HashJoin.

For the final tier write, the existing `ScratchGuard` and truncating temporary
writers make replay safe. Probe construction must similarly discard a partially
collected probe before retrying.

Add an optional diagnostic override `auto | hash | sort_merge`, defaulting to
`auto`. It is for tests and benchmarking only; it must not accept chromosome
rules. Forced HashJoin still runs inside the finite pool and fails rather than
silently escaping the memory bound.

`TracingPool` must forward `memory_limit()` to the wrapped pool. Enabling
`VEPYR_TRACE_MEMORY` must not change the strategy decision.

## Worktree preservation gate

This worktree already contains user-owned experimental edits in `build.rs` and
`mem_trace.rs`. Before implementation:

1. Record `git status --short` and capture the exact diff of both files in the
   execution notes.
2. Do not use `git reset`, `git restore`, or `git checkout` on either file.
3. Supersede only the known-invalid configuration experiment: the comments
   claiming both inputs are sorted, `prefer_hash_join=false` under
   `target_partitions=1`, and the two zero collection thresholds.
4. Preserve and integrate `mem_trace::describe_batches` and its call site. If a
   new statistics helper replaces it, retain equivalent charged-vs-distinct
   diagnostics and tests rather than deleting the capability.
5. Compare the final overlapping-file diff with the captured baseline before
   committing.

## Observability

Emit one `info` record per decision and one per fallback:

```text
plugin_cache[cadd/chr21]: join_strategy stage=probe_build algorithm=hash
  left_rows=14574753 left_bytes=... right_rows=121809062 right_bytes=...
  hash_build_estimate=... pool_available=... partitions=... confidence=estimated
  reason=hash_build_fits
```

The log fields must include stage, both input statistics and their confidence,
pool limit/reserved/available, target partitions, chosen algorithm, estimated
hash bytes, decision reason, and whether the choice is an automatic fallback.
Keep `VEPYR_EXPLAIN_TIER_JOIN` and `VEPYR_TRACE_MEMORY` as opt-in deeper evidence.

## Implementation tasks

### Task 1: Pure statistics and strategy selector

**Files:**
- Create `datafusion/bio-function-vep/src/plugin_cache/join_strategy.rs`
- Modify `datafusion/bio-function-vep/src/plugin_cache/mod.rs`

**Work:**
- Add the types and pure decision function above.
- Mirror DataFusion's hash-build reservation components using public
  `datafusion_common::utils::memory` helpers.
- Represent exact, estimated, and unknown statistics explicitly.
- Add structured decision reasons; do not store free-form policy strings.

**Tests:**
- Same inputs under two pool sizes choose Hash then SMJ.
- Swapping input sizes picks the smaller hash-build estimate.
- Unknown/inexact-without-sample statistics choose SMJ.
- Existing reservations reduce available memory.
- Partition count derives from runtime default and available input partitions.
- Tests never pass a chromosome or plugin name to the selector.

### Task 2: Collect projected input statistics

**Files:**
- Modify `datafusion/bio-function-vep/src/plugin_cache/join_strategy.rs`
- Modify `datafusion/bio-function-vep/src/plugin_cache/join.rs`

**Work:**
- Add exact `RecordBatch`/MemTable statistics collection by projection.
- Materialize the plugin's variation-key projection as owned MemTable batches;
  apply spillable `DISTINCT` only when match columns make duplicate variation
  keys possible.
- Add Parquet row-count/projected-width collection and bounded sampling when
  Arrow byte size is not exact.
- Calculate statistics before planning probe construction, then recalculate
  exact statistics after probe materialization.
- Verify the sampled scan does not materialize the full shard and uses only the
  join-key projection.

**Tests:**
- Synthetic Parquet metadata reports the exact row count.
- Variable-width strings produce a non-zero projected Arrow estimate.
- Projection excludes plugin value columns.
- Owned probe statistics equal the sum of batch memory sizes.
- The physical `plugin_variation_keys` table reports exact rows and bytes.

### Task 3: Strategy-specific parallel contexts

**Files:**
- Modify `datafusion/bio-function-vep/src/plugin_cache/build.rs`
- Modify `datafusion/bio-function-vep/src/plugin_cache/join.rs`

**Work:**
- Replace `tier_join_config()` with a constructor accepting `JoinAlgorithm` and
  runtime target partitions.
- Keep `schema_force_view_types=false` for every strategy.
- Use one finite `FairSpillPool` for the tiering stages; remove the unbounded
  fast join context.
- Partition MemTable batches by accumulated bytes without copying or slicing.
- Configure Hash and SMJ exactly as specified in the configuration table.
- Inspect the optimized Hash physical plan and reject it before execution unless
  it is CollectLeft, its build child maps to a measured input, and the
  recomputed estimate for that actual child fits.

**Plan tests:**
- Hash decision yields verified `HashJoinExec(mode=CollectLeft)` with an
  identified/measured build schema and parallel probe partitions.
- SMJ decision yields `SortMergeJoinExec`, hash repartitioning, and input sorts.
- One-batch inputs and a configured runtime default of one still yield genuine
  `SortMergeJoinExec` with two logical partitions.
- Neither plan claims pre-existing full-key ordering.

### Task 4: Apply selection to both joins

**Files:**
- Modify `datafusion/bio-function-vep/src/plugin_cache/join.rs`
- Modify `datafusion/bio-function-vep/src/plugin_cache/build.rs`

**Work:**
- Split probe creation from final tier streaming so each accepts its own
  `JoinDecision` and context.
- Select the probe-build strategy from raw variation plus projected plugin-key
  statistics.
- Materialize/rechunk the probe, then select the final tier strategy from exact
  probe plus exact plugin MemTable statistics.
- Log both decisions independently.
- Inspect physical plans before execution and in tests to ensure the requested
  algorithm, Hash partition mode, build-side identity, and recomputed build
  estimate match the decision.

**Tests:**
- A large variation/small pool fixture selects SMJ while a smaller fixture on
  the same code path selects Hash.
- Probe-build and final-tier decisions can differ.
- Duplicate/conflicting variation tiers still collapse with `MIN(tier)` and
  warm wins.
- Output rows and tier counts match the current implementation.

### Task 5: One-pass sorted writer and memory fallback

**Files:**
- Modify `datafusion/bio-function-vep/src/plugin_cache/build.rs`
- Modify `datafusion/bio-function-vep/src/plugin_cache/join.rs`

**Work:**
- Use explicit sorted output for all parallel final-join plans.
- Retain `assert_start_monotonic` as an invariant check.
- Remove order-violation-driven replay and replace it with an attributed
  `HashJoinInput`-rejection-driven Hash-to-SMJ retry around both adaptive stages.
- Ensure failed streams, reservations, temporary probe batches, and writer
  handles are dropped before replay.
- Allow at most one fallback per stage.

**Tests:**
- Inject deliberately optimistic statistics in a test so auto selects Hash,
  then use a stage-local test pool to reject the `HashJoinInput` reservation;
  assert SMJ replay succeeds with identical rows.
- Aggregate and sorter `ResourcesExhausted` errors do not trigger join fallback.
- Non-memory errors do not trigger fallback.
- Warm and cold temporary outputs are safely overwritten after replay.
- Final shard remains sorted by `(tier, start)` and page lookup returns every
  written row.

### Task 6: Diagnostics and regression coverage

**Files:**
- Modify `datafusion/bio-function-vep/src/plugin_cache/mem_trace.rs`
- Modify tests in `build.rs` and `join.rs`

**Work:**
- Add stable structured decision logging and forced-strategy test override.
- Forward `memory_limit()` through `TracingPool` and record the consumer label
  for rejected reservations.
- Preserve the owned execution-sized probe regression from `92f55f9`.
- Supersede only the current invalid configuration experiment claiming that
  `prefer_hash_join=false` works with `target_partitions=1` and that both inputs
  are already sorted; preserve/integrate `describe_batches` per the worktree
  preservation gate.

**Tests:**
- Selector log formatting contains all required fields.
- Tracing enabled and disabled produce identical `JoinDecision` values.
- A rejected aggregate/sorter reservation records a non-Hash consumer label and
  is propagated unchanged.
- `cargo test -p datafusion-bio-function-vep plugin_cache --lib` passes.
- `cargo clippy -p datafusion-bio-function-vep --all-targets --all-features -- -D warnings`
  passes.

### Task 7: Real-data gate and tuning evidence

Run two separate gates:

1. CADD chr21 under the same 8 GiB pool: forced Hash, forced SMJ, and auto must
   all succeed and produce equivalent output.
2. One large chromosome under a constrained pool: auto and forced SMJ must
   succeed; forced Hash must fail specifically with an attributed
   `HashJoinInput` exhaustion. Optionally rerun forced Hash under a separately
   documented larger pool for output parity.

Acceptance criteria:

- Auto selects Hash for chr21 if the computed estimate fits.
- Auto selects SMJ before allocation when a large chromosome's estimated hash
  build exceeds the pool.
- If an estimate is optimistic, automatic fallback completes without manual
  rerun.
- All successful modes produce identical total/warm/cold counts and
  byte-equivalent plugin values.
- `HashJoinInput` never exceeds the finite pool; the large-input SMJ path shows
  spillable sort reservations instead.
- Parallel plans have more than one output partition before the final global
  ordering merge.
- Record wall time, peak RSS, spill bytes, chosen strategies, input estimates,
  and actual peak reservations in the handover.
- Auto chr21 wall time is no more than 10% above the recorded 118.2-second
  HashJoin baseline on the same host/input, or the regression is investigated
  before rollout.

## Executable plan split

Execute this design as three reviewable plans; overlapping Rust files are owned
serially, while the implemented DataFusion pipeline is parallel at runtime.

| Plan | Tasks | Result |
|---|---|---|
| A — statistics contract | 1, 2 | Pure selector, exact projected keys, Parquet estimates |
| B — execution contract | 3, 4, 5 | Verified physical operator/build side, both joins adaptive, attributed replay |
| C — integration gate | 6, 7 | Worktree-safe diagnostics, regression suite, real-data validation |

## Definition of done

- No production strategy branch contains a chromosome or plugin-name match.
- No fixed row-count threshold controls HashJoin vs SMJ.
- Both joins use finite-pool, parallel contexts.
- The decision is reproducible from logged statistics and pool state.
- Unknown statistics choose the bounded path.
- An attributed HashJoin memory exhaustion automatically retries once with SMJ;
  other resource errors do not.
- Output order, tier semantics, dedup behavior, cache schema, and point lookup
  remain unchanged.
- CADD chr21 auto mode stays within 10% of the 118.2-second baseline on the
  same host/input; a large chromosome completes under the configured pool
  without manual strategy changes.
