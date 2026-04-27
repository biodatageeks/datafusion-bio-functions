# Apple GPU Count Overlaps Backend

## Why

`count_overlaps(left, right)` is the safest first interval operation to move onto Apple GPU:

- output cardinality is fixed: exactly one count per right-side row
- no overlap-pair materialization is required
- no dynamic output buffers or skew-sensitive pair emission are required
- the core computation can be expressed as sorted endpoint rank queries:
  - `started = count(left.start <= right.end)`
  - `ended_before = count(left.end < right.start)`
  - `overlaps = started - ended_before`

That makes `count_overlaps` a better first Metal backend than `overlap()`. It exercises the same data movement, per-contig indexing, and kernel dispatch infrastructure that a future GPU overlap join will need, while avoiding the hard all-pairs output problem.

The current implementation collects the left side into per-contig CPU COITrees, then probes right batches row-by-row in Rust. That is a strong CPU path, but it is still lookup-oriented and does not expose the GPU-friendly sort/rank shape of count-overlaps.

This change plans an optional Apple GPU backend that preserves existing CPU semantics and falls back automatically when GPU execution is not beneficial or unavailable.

## What Changes

### 1. Add a count-overlaps backend abstraction

Introduce an internal backend boundary for `CountOverlapsExec`:

- `CpuCoitreeCountBackend`: existing implementation
- `AppleGpuCountBackend`: optional Metal implementation, compiled only on supported Apple targets and behind a cargo feature
- `AutoCountBackend`: dispatches to GPU only when the device, input shape, and workload size make sense

The public SQL/table function contract does not change.

### 2. Build a GPU-ready left index

For each contig in the left table, build structure-of-arrays buffers:

- sorted starts
- sorted ends
- per-contig offsets and lengths
- dictionary-encoded contig IDs

The first implementation should build these arrays on CPU and pass them to Metal using unified/shared memory. GPU-side radix sort is useful follow-up infrastructure, but it is not required to prove the count-overlaps backend.

### 3. Execute right batches through Metal rank kernels

For each right-side batch:

- dictionary-encode contigs to the left-index contig IDs
- normalize strict/weak interval semantics exactly as the CPU path does
- dispatch GPU rank kernels
- append the resulting `Int64` count column to the original right batch
- preserve right input row order exactly

The implementation should start with a simple one-thread-per-query binary-rank kernel because it avoids sorting every right batch. For large sorted batches, add a Merge Path rank-sweep path that computes the same endpoint ranks with less per-query divergence.

### 4. Keep CPU as the default safety net

GPU execution is opt-in or auto-selected only above a measured threshold. The CPU COITree path remains the default fallback for:

- non-Apple platforms
- unsupported column types
- very small inputs or batches
- GPU device/runtime initialization failure
- memory pressure or allocation failure
- correctness/debug mode

## Impact

- Affected code:
  - `datafusion/bio-function-ranges/src/count_overlaps.rs`
  - `datafusion/bio-function-ranges/src/interval_tree.rs`
  - `datafusion/bio-function-ranges/src/session_context.rs`
  - new Apple GPU backend module under `datafusion/bio-function-ranges/src/`
  - `datafusion/bio-function-ranges/Cargo.toml`
  - integration tests under `datafusion/bio-function-ranges/tests/`
- Public behavior:
  - no SQL syntax change
  - existing CPU output remains the correctness oracle
  - optional session config controls backend selection
- Expected outcome:
  - large `count_overlaps` workloads on M3 Max/Ultra can use GPU bandwidth while small workloads continue using CPU COITrees
  - the repo gains reusable Apple GPU buffer/kernels infrastructure for future range operations
