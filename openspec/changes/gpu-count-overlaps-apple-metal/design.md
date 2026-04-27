# Apple GPU Count Overlaps Backend -- Design

## Context

The current `count_overlaps` table function builds CPU COITrees from the full left table and probes each right-side batch on CPU. The relevant flow is:

1. `CountOverlapsProvider::scan()` collects the left table into memory.
2. `build_coitree_from_batches()` groups intervals by contig and builds `COITree<(), u32>` values.
3. `CountOverlapsExec::execute()` streams right partitions.
4. `get_stream()` probes each right row and appends an `Int64` `count` column.

That implementation is efficient on CPU and must remain the correctness baseline. The GPU backend should not replace the public operator. It should provide an alternative physical implementation for large workloads where fixed-output rank computation can amortize device setup and kernel launch overhead.

## Goals / Non-Goals

**Goals:**

- Preserve exact `count_overlaps` output semantics for weak and strict modes.
- Keep output shape and row order identical to the existing CPU implementation.
- Add a feature-gated Apple GPU backend for large count-overlaps workloads.
- Use a GPU-friendly endpoint rank formulation:
  - count left starts at or before each query end
  - count left ends before each query start
  - subtract ranks per query
- Keep CPU COITrees as fallback and correctness oracle.
- Build reusable buffer, device detection, and benchmark infrastructure for later GPU interval operations.

**Non-Goals:**

- Implementing all-pairs `overlap()` on GPU.
- Implementing `coverage`, `nearest`, `merge`, or `cluster` in this change.
- Replacing COITrees as the default small-workload path.
- Requiring GPU-side sorting in the first implementation.
- Changing SQL syntax, output schema, or null behavior.

## Algorithm

For inclusive intervals, a left interval `[l_start, l_end]` overlaps a right interval `[r_start, r_end]` when:

```text
l_start <= r_end AND l_end >= r_start
```

For a fixed contig, if left starts and ends are sorted independently:

```text
started     = upper_bound(left_starts, r_end)
ended_prior = lower_bound(left_ends, r_start)
count       = started - ended_prior
```

This is containment-correct for counts because it counts every interval satisfying the overlap predicate without needing interval IDs. AIList-style sublists are not required for `count_overlaps`; they become relevant for future pair enumeration where interval identity matters.

Strict mode should be normalized to match the existing CPU code:

```text
strict query start = right.start + 1
strict query end   = right.end - 1
```

Then the same rank formula is applied. Empty strict intervals must produce zero.

## Detailed Design

### 1. Backend boundary

Add an internal trait for count-overlaps execution:

```rust
trait CountOverlapsBackend: Send + Sync {
    fn name(&self) -> &'static str;
    fn execute_partition(
        &self,
        right_plan: Arc<dyn ExecutionPlan>,
        schema: SchemaRef,
        columns_2: Arc<(String, String, String)>,
        filter_op: FilterOp,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream>;
}
```

The existing `get_stream()` path can become the CPU backend implementation. The provider chooses the backend during `scan()` after building or preparing the left-side index.

The initial backend enum should be internal:

```rust
enum CountOverlapsBackendMode {
    Cpu,
    AppleGpu,
    Auto,
}
```

Expose it through a session option such as:

```sql
SET bio.count_overlaps_backend = 'auto';      -- default
SET bio.count_overlaps_backend = 'cpu';
SET bio.count_overlaps_backend = 'apple_gpu';
```

`auto` only selects GPU when all eligibility checks pass.

### 2. Feature gate and platform gating

Add a cargo feature to `datafusion/bio-function-ranges`:

```toml
[features]
apple-gpu = []
```

The implementation should compile cleanly without the feature on every platform. Metal-specific code must be behind:

```rust
#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
```

Use `objc2-metal` as the initial Metal binding, with companion `objc2` and `objc2-foundation` dependencies where needed. Keep it behind a small local wrapper so the rest of the crate depends on project-owned types, not directly on the binding API.

### 3. GPU index layout

Create a left-side GPU index:

```rust
struct GpuCountOverlapsIndex {
    chrom_to_id: AHashMap<String, u32>,
    contigs: Vec<ContigRange>,
    starts_sorted: Vec<i64>,
    ends_sorted: Vec<i64>,
}

struct ContigRange {
    chrom_id: u32,
    start_offset: u32,
    end_offset: u32,
    len: u32,
}
```

The first implementation builds this index on CPU from collected left batches:

1. Read left `(contig, start, end)` columns using existing array helpers.
2. Normalize coordinates to the same type used by the CPU path.
3. Group by contig.
4. Sort starts and ends independently for each contig.
5. Concatenate all contig groups into global SoA buffers.
6. Store offsets/lengths in `contigs`.

Use `i64` in the GPU-facing buffers unless benchmarks prove `i32` is required. The current CPU COITree path narrows to `i32` in some places, but the SQL/table function layer works with `Int64` data. Planning the GPU path around `i64` avoids silent coordinate truncation.

### 4. Right batch encoding

For each right batch, prepare GPU input arrays:

```rust
struct GpuQueryBatch {
    chrom_ids: Vec<u32>,
    starts: Vec<i64>,
    ends: Vec<i64>,
}
```

Encoding rules:

- map unknown contigs to `u32::MAX`
- apply strict-mode start/end adjustment before dispatch
- mark invalid strict intervals where adjusted start > adjusted end
- keep physical row order unchanged
- produce count `0` for unknown contigs or invalid intervals

The GPU output buffer is one `i64` count per input row.

### 5. Kernel strategy

#### Phase 1 kernel: binary-rank per query

Start with one thread per right row:

1. Look up the contig offset and length.
2. Binary-search starts for `upper_bound(end)`.
3. Binary-search ends for `lower_bound(start)`.
4. Write `started - ended_prior`.

This has predictable fixed output and no atomics. It is the simplest production-worthy kernel and avoids sorting every right batch.

#### Phase 2 kernel: Merge Path rank-sweep for large sorted batches

Add a second path only after the binary-rank path is correct and benchmarked.

The Merge Path path should be selected when:

- a right-side batch or partition is already sorted by `(contig, endpoint)`, or
- benchmarking shows that sorting endpoint views pays for the target workload size

It computes the same two ranks with two sorted merges:

1. Merge left starts with right query ends to compute `started`.
2. Merge left ends with right query starts to compute `ended_prior`.
3. Subtract ranks by original row index.

Merge Path partitioning assigns independent cross-diagonal merge segments to threadgroups. For `count_overlaps`, this is a rank-sweep, not an active-list pair-emission sweep. That distinction keeps memory bounded and avoids the 32 KB threadgroup active-list problem that makes all-pairs overlap hard on Apple GPU.

Do not sort every small or unsorted right batch just to use Merge Path. That would likely erase the GPU win.

### 6. Execution integration

`CountOverlapsProvider::scan()` should become:

1. collect/select left columns as today
2. inspect session config for backend mode
3. build CPU COITree backend if mode is `cpu` or GPU is ineligible
4. build GPU SoA index if mode is `apple_gpu` or eligible `auto`
5. create a `CountOverlapsExec` that owns the selected backend

The GPU backend should execute per right partition, same as the CPU path. It should not coalesce all right partitions unless a later Merge Path implementation needs whole-partition sorting and benchmarks justify it.

### 7. Eligibility and fallback

`auto` mode should require:

- `target_os = macos`
- `apple-gpu` feature enabled
- Metal device initialization succeeds
- left row count above a benchmark-derived threshold
- right batch or partition size above a benchmark-derived threshold
- supported coordinate types
- no unsupported null semantics

Any runtime failure before emitting a batch should fall back to CPU. Once a GPU stream has emitted data for a partition, failures should surface as execution errors rather than silently mixing partial outputs from different engines.

### 8. Observability

Add lightweight metrics:

- selected backend
- GPU eligibility decision and fallback reason
- left index build time
- right rows processed on GPU
- GPU kernel elapsed time
- CPU fallback rows

Metrics should be exposed through DataFusion execution metrics where practical and optionally logged at debug level.

## Validation Plan

### Correctness

Add tests comparing CPU and GPU-count-reference behavior for:

- simple overlaps
- no overlaps
- point intervals
- nested/contained intervals
- duplicate intervals
- missing contigs
- strict versus weak mode
- multi-partition right input
- output row order preservation

The tests should run without a real GPU by validating the shared rank implementation against the CPU COITree path. Metal-specific integration tests can be gated and ignored by default when no Apple GPU is available.

### Performance

Add ignored benchmarks for:

- CPU COITree baseline
- CPU rank-reference baseline over sorted endpoint arrays
- GPU binary-rank path
- GPU Merge Path path when implemented

Benchmark matrix:

- left rows: 1e5, 1e6, 1e7, 5e7
- right rows: 1e5, 1e6, 1e7
- contig distributions: single contig, human-chromosome-like, skewed contig
- interval distributions: short uniform, nested/contained, mixed long-tail
- right ordering: sorted, mostly sorted, random

The auto-dispatch threshold should be set from these benchmarks, not guessed.
