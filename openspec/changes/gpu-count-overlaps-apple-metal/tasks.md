# Implementation Tasks

## Status Key
- [ ] Not started
- [x] Completed

## 0. Phase 0: Baseline and Semantics

- [x] 0.1 Add focused CPU tests documenting current `count_overlaps` weak-mode semantics
- [x] 0.2 Add focused CPU tests documenting current `count_overlaps` strict-mode semantics
- [x] 0.3 Add tests for nested/contained intervals to prove count semantics are containment-correct
- [x] 0.4 Add tests for missing contigs returning zero counts
- [x] 0.5 Add tests proving output row order matches right input order
- [x] 0.6 Add a small benchmark harness for current CPU COITree `count_overlaps`

## 1. Phase 1: Backend Selection Surface

- [x] 1.1 Add `CountOverlapsBackendMode` enum: `Auto`, `Cpu`, `AppleGpu`
- [x] 1.2 Add `bio.count_overlaps_backend` session config with default `Auto`
- [x] 1.3 Parse accepted values: `auto`, `cpu`, `apple_gpu`
- [x] 1.4 Add display/debug output showing the selected count-overlaps backend in physical plans
- [x] 1.5 Keep current behavior exactly when mode is `Cpu`
- [x] 1.6 Add tests for session config parsing and invalid values

## 2. Phase 2: CPU Rank Reference and GPU-Ready Index

- [x] 2.1 Add `GpuCountOverlapsIndex`-style SoA index type behind a backend-neutral module name
- [x] 2.2 Build per-contig `starts_sorted` and `ends_sorted` arrays from collected left batches
- [x] 2.3 Dictionary-encode contigs to stable `u32` IDs
- [x] 2.4 Use `i64` coordinate buffers and explicitly document any conversion from existing CPU helpers
- [x] 2.5 Implement CPU rank-reference count: `upper_bound(starts, end) - lower_bound(ends, start)`
- [ ] 2.6 Compare CPU rank-reference output against existing COITree output across all Phase 0 fixtures
- [ ] 2.7 Add property-style tests for random interval/query sets comparing rank-reference and COITree counts
- [ ] 2.8 Measure CPU index build time and memory footprint versus COITree build

## 3. Phase 3: Feature Gate and Metal Runtime Wrapper

- [x] 3.1 Add `apple-gpu` cargo feature to `datafusion/bio-function-ranges`
- [x] 3.2 Add `#[cfg(all(feature = "apple-gpu", target_os = "macos"))]` module boundary
- [x] 3.3 Add `objc2`, `objc2-foundation`, and `objc2-metal` dependencies behind the `apple-gpu` feature
- [x] 3.4 Wrap Metal device, queue, pipeline, and buffer allocation behind project-owned types
- [ ] 3.5 Implement device detection returning structured ineligibility reasons
- [x] 3.6 Verify non-macOS and no-feature builds compile without Metal dependencies
- [x] 3.7 Add a tiny ignored Metal smoke test that initializes the device and runs a no-op kernel

## 4. Phase 4: GPU Binary-Rank Kernel

- [x] 4.1 Define Metal buffer layouts for contig ranges, sorted starts, sorted ends, query chrom IDs, query starts, query ends, and output counts
- [x] 4.2 Implement one-thread-per-query `upper_bound(starts, query_end)` kernel logic
- [x] 4.3 Implement one-thread-per-query `lower_bound(ends, query_start)` kernel logic
- [x] 4.4 Combine both ranks in a single count kernel when profiling shows it is faster than two kernels
- [x] 4.5 Return zero for missing contigs and invalid strict intervals
- [x] 4.6 Preserve output order by writing one count per original right row
- [x] 4.7 Add GPU-vs-CPU correctness tests gated on Apple GPU availability
- [ ] 4.8 Add debug assertions comparing a sample of GPU counts to CPU rank-reference counts in test builds

## 5. Phase 5: `CountOverlapsExec` Integration

- [ ] 5.1 Refactor existing `get_stream()` path into `CpuCoitreeCountBackend`
- [x] 5.2 Add `AppleGpuCountBackend` that owns the SoA index and Metal runtime wrapper
- [x] 5.3 Update `CountOverlapsProvider::scan()` to select CPU/GPU/auto backend
- [x] 5.4 Encode right batches into GPU query buffers per partition
- [x] 5.5 Append GPU counts as an `Int64Array` to each original right batch
- [x] 5.6 Preserve existing repartition behavior for right plans
- [x] 5.7 Add fallback before first batch emission when GPU initialization or allocation fails
- [ ] 5.8 Emit execution metrics for backend choice, fallback reason, rows processed, and kernel time

## 6. Phase 6: Auto Dispatch Thresholds

- [x] 6.1 Add ignored benchmark for CPU COITree count-overlaps baseline
- [ ] 6.2 Add ignored benchmark for CPU rank-reference baseline
- [x] 6.3 Add ignored benchmark for GPU binary-rank path
- [ ] 6.4 Benchmark left sizes: 1e5, 1e6, 1e7, 5e7 where locally feasible
- [ ] 6.5 Benchmark right sizes: 1e5, 1e6, 1e7 where locally feasible
- [ ] 6.6 Benchmark sorted, mostly sorted, and random right input order
- [ ] 6.7 Benchmark nested/contained interval distributions
- [x] 6.8 Set initial `Auto` thresholds from measured crossover points
- [x] 6.9 Document threshold rationale in code comments and README/API docs

## 7. Phase 7: Merge Path Rank-Sweep Follow-Up

- [ ] 7.1 Detect whether a right partition is already sorted by `(contig, endpoint)` using physical ordering or a cheap runtime check
- [ ] 7.2 Build endpoint views for query starts and query ends with original row IDs
- [ ] 7.3 Benchmark CPU endpoint-view sorting cost versus binary-rank kernel cost
- [ ] 7.4 Implement Merge Path partition calculation for one contig-local merge
- [ ] 7.5 Compute `started` by merging left starts with sorted query ends
- [ ] 7.6 Compute `ended_prior` by merging left ends with sorted query starts
- [ ] 7.7 Scatter ranks back to original row IDs and subtract counts
- [ ] 7.8 Dispatch Merge Path only when benchmarks show it beats binary-rank for the current input shape

## 8. Phase 8: Documentation and Release Readiness

- [ ] 8.1 Document `bio.count_overlaps_backend`
- [ ] 8.2 Document platform and feature requirements for Apple GPU support
- [ ] 8.3 Document that CPU remains the correctness baseline and fallback path
- [ ] 8.4 Add troubleshooting notes for GPU ineligibility/fallback reasons
- [x] 8.5 Run `cargo test -p datafusion-bio-function-ranges`
- [x] 8.6 Run tests with `--features apple-gpu` on Apple Silicon
- [x] 8.7 Run ignored benchmarks on at least one M3-class machine before enabling `Auto` GPU dispatch by default
