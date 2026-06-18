# Within-Contig Parallel Annotation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make annotation of a single contig fast by replacing the one-worker-per-contig serial annotation with a work-stealing pool of `N` window workers, each running a fused lookup→annotate→format pipeline over a contiguous position window, reassembled into output order.

**Architecture:** A `WindowPlanner` cuts a contig's position-sorted input batches into windows (~one input buffer each). A work-stealing pool of `N` `WindowWorker`s each take a window, run their own lookup over that window's variants (self-contained `colocated_map`), annotate via the existing `annotate_worker_window` engine path, and format. An `OrderedWindowDrain` releases completed windows in window-index order. A single `threads = N` knob replaces the `forks`/`target_partitions`/`contig_parallelism`/`inline_lookup` option cluster (clean break).

**Tech Stack:** Rust 2024, Apache DataFusion 53 / Arrow 58, Tokio (async streams + `tokio::spawn` + `mpsc`), the existing `TranscriptConsequenceEngine` and `ColocatedSink` lookup machinery.

**Design spec:** `docs/superpowers/specs/2026-06-17-within-contig-parallel-annotation-design.md`

---

## File Structure

| File | Responsibility | New/Modify |
|---|---|---|
| `datafusion/bio-function-vep/src/window_planner.rs` | Pure: group position-sorted input batches into indexed `Window`s | **New** |
| `datafusion/bio-function-vep/src/ordered_drain.rs` | Pure: reassemble out-of-order window outputs into window-index order | **New** |
| `datafusion/bio-function-vep/src/lib.rs` (or crate root with `mod` decls) | Register the two new modules | Modify |
| `datafusion/bio-function-vep/src/annotate_provider.rs` | Per-window fused lookup helper; work-stealing pool driver; `threads` parser; replace `AnnotatingContig` serial logic | Modify |
| `datafusion/bio-function-vep/src/vcf_sink.rs` | `AnnotateVcfConfig`: drop `forks`/`workers`/`target_partitions`, add `threads` | Modify |
| `datafusion/bio-function-vep/examples/*.rs`, `scripts/*.py` | Update callers to the `threads` knob | Modify |

**Module placement rationale:** `WindowPlanner` and `OrderedWindowDrain` are pure and standalone → their own files with full unit tests. The fused worker + pool driver must touch private items in `annotate_provider.rs` (`AnnotationWorkerState`, `annotate_worker_window`, `ColocatedSink`, `drain_colocated_sink`) — adding them in-module avoids a large visibility refactor of an already-12k-line file.

---

## Phase 1 — `WindowPlanner` (pure, new module)

### Task 1.1: Create the `WindowPlanner` module with a failing test

**Files:**
- Create: `datafusion/bio-function-vep/src/window_planner.rs`
- Modify: crate root module list (`datafusion/bio-function-vep/src/lib.rs`)

- [ ] **Step 1: Register the module**

In the crate root (`lib.rs`), add alongside the other `mod` declarations:

```rust
mod ordered_drain;
mod window_planner;
```

- [ ] **Step 2: Write the module with the type and a failing test**

Create `datafusion/bio-function-vep/src/window_planner.rs`:

```rust
//! Groups a contig's position-sorted input batches into indexed windows.
//!
//! A window is a contiguous run of whole input batches totalling at least
//! `target_variants` rows (no mid-batch slicing). The final window may be
//! smaller. Windows are the unit of work handed to the parallel window pool.

use datafusion::arrow::record_batch::RecordBatch;

/// A contiguous group of position-sorted input batches for one contig.
pub(crate) struct Window {
    pub index: usize,
    pub batches: Vec<RecordBatch>,
}

/// Accumulates input batches and emits a [`Window`] each time the pending
/// group reaches `target_variants` rows.
pub(crate) struct WindowPlanner {
    target_variants: usize,
    pending: Vec<RecordBatch>,
    pending_rows: usize,
    next_index: usize,
}

impl WindowPlanner {
    pub(crate) fn new(target_variants: usize) -> Self {
        Self {
            target_variants: target_variants.max(1),
            pending: Vec::new(),
            pending_rows: 0,
            next_index: 0,
        }
    }

    /// Push one input batch. Returns a [`Window`] when the pending group
    /// reaches `target_variants`. Empty batches are ignored.
    pub(crate) fn push(&mut self, batch: RecordBatch) -> Option<Window> {
        let rows = batch.num_rows();
        if rows == 0 {
            return None;
        }
        self.pending.push(batch);
        self.pending_rows += rows;
        if self.pending_rows >= self.target_variants {
            Some(self.flush())
        } else {
            None
        }
    }

    /// Emit any remaining pending batches as a final (possibly smaller) window.
    pub(crate) fn finish(&mut self) -> Option<Window> {
        if self.pending.is_empty() {
            None
        } else {
            Some(self.flush())
        }
    }

    fn flush(&mut self) -> Window {
        let index = self.next_index;
        self.next_index += 1;
        let batches = std::mem::take(&mut self.pending);
        self.pending_rows = 0;
        Window { index, batches }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::Int64Array;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn batch(rows: usize) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![Field::new("pos", DataType::Int64, false)]));
        let arr = Int64Array::from((0..rows as i64).collect::<Vec<_>>());
        RecordBatch::try_new(schema, vec![Arc::new(arr)]).unwrap()
    }

    #[test]
    fn groups_whole_batches_to_target_and_flushes_remainder() {
        let mut planner = WindowPlanner::new(5000);
        // 3000 pending -> below target, no window yet.
        assert!(planner.push(batch(3000)).is_none());
        // +3000 = 6000 >= 5000 -> emit window 0 with both batches.
        let w0 = planner.push(batch(3000)).expect("window 0");
        assert_eq!(w0.index, 0);
        assert_eq!(w0.batches.iter().map(|b| b.num_rows()).sum::<usize>(), 6000);
        // Remaining 3000 emitted by finish() as window 1.
        assert!(planner.push(batch(3000)).is_none());
        let w1 = planner.finish().expect("window 1");
        assert_eq!(w1.index, 1);
        assert_eq!(w1.batches.iter().map(|b| b.num_rows()).sum::<usize>(), 3000);
        // Nothing left.
        assert!(planner.finish().is_none());
    }

    #[test]
    fn ignores_empty_batches() {
        let mut planner = WindowPlanner::new(10);
        assert!(planner.push(batch(0)).is_none());
        assert!(planner.finish().is_none());
    }
}
```

- [ ] **Step 3: Run the tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep window_planner -- --nocapture`
Expected: PASS (2 tests). If the package name differs, find it with `cargo metadata --no-deps --format-version 1 | grep -o '"name":"[^"]*vep[^"]*"'`.

- [ ] **Step 4: Lint and commit**

```bash
cargo fmt
cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/window_planner.rs datafusion/bio-function-vep/src/lib.rs
git commit -m "feat(vep): add WindowPlanner to group contig batches into windows"
```

---

## Phase 2 — `OrderedWindowDrain` (pure, new module)

### Task 2.1: Create the ordered-drain reorder buffer with tests

**Files:**
- Create: `datafusion/bio-function-vep/src/ordered_drain.rs`

- [ ] **Step 1: Write the module with type and failing tests**

Create `datafusion/bio-function-vep/src/ordered_drain.rs`:

```rust
//! Reassembles window-worker outputs into window-index order.
//!
//! Workers complete out of order (work-stealing). Each completed window is
//! recorded with its index; the drain releases batches only once every
//! lower-indexed window has been released, so output order is deterministic
//! regardless of completion order.

use std::collections::HashMap;

use datafusion::arrow::record_batch::RecordBatch;

pub(crate) struct OrderedWindowDrain {
    next_emit: usize,
    pending: HashMap<usize, Vec<RecordBatch>>,
}

impl OrderedWindowDrain {
    pub(crate) fn new() -> Self {
        Self {
            next_emit: 0,
            pending: HashMap::new(),
        }
    }

    /// Record a completed window. Returns the batches now releasable in order
    /// (this window if it is `next_emit`, plus any contiguous buffered
    /// successors). Returns empty if this window must wait for an earlier one.
    pub(crate) fn complete(&mut self, index: usize, batches: Vec<RecordBatch>) -> Vec<RecordBatch> {
        self.pending.insert(index, batches);
        let mut released = Vec::new();
        while let Some(batches) = self.pending.remove(&self.next_emit) {
            released.extend(batches);
            self.next_emit += 1;
        }
        released
    }

    /// Number of completed-but-not-yet-released windows held in the buffer.
    pub(crate) fn buffered_windows(&self) -> usize {
        self.pending.len()
    }
}

impl Default for OrderedWindowDrain {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::Int64Array;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn batch(tag: i64) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![Field::new("t", DataType::Int64, false)]));
        RecordBatch::try_new(schema, vec![Arc::new(Int64Array::from(vec![tag]))]).unwrap()
    }

    fn tags(batches: &[RecordBatch]) -> Vec<i64> {
        batches
            .iter()
            .map(|b| {
                b.column(0)
                    .as_any()
                    .downcast_ref::<Int64Array>()
                    .unwrap()
                    .value(0)
            })
            .collect()
    }

    #[test]
    fn releases_in_order_when_completed_in_order() {
        let mut drain = OrderedWindowDrain::new();
        assert_eq!(tags(&drain.complete(0, vec![batch(0)])), vec![0]);
        assert_eq!(tags(&drain.complete(1, vec![batch(1)])), vec![1]);
        assert_eq!(drain.buffered_windows(), 0);
    }

    #[test]
    fn buffers_out_of_order_then_releases_contiguously() {
        let mut drain = OrderedWindowDrain::new();
        // Window 2 finishes first -> nothing released, buffered.
        assert!(drain.complete(2, vec![batch(2)]).is_empty());
        // Window 1 finishes next -> still waiting on 0.
        assert!(drain.complete(1, vec![batch(1)]).is_empty());
        assert_eq!(drain.buffered_windows(), 2);
        // Window 0 finishes -> releases 0, 1, 2 in order.
        assert_eq!(tags(&drain.complete(0, vec![batch(0)])), vec![0, 1, 2]);
        assert_eq!(drain.buffered_windows(), 0);
    }
}
```

- [ ] **Step 2: Run the tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep ordered_drain -- --nocapture`
Expected: PASS (2 tests).

- [ ] **Step 3: Lint and commit**

```bash
cargo fmt
cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/ordered_drain.rs
git commit -m "feat(vep): add OrderedWindowDrain reorder buffer for parallel windows"
```

---

## Phase 3 — Per-window fused lookup helper

**Context:** Today `prepare_contig_context` (`annotate_provider.rs:11532`) builds **one** lookup `ExecutionPlan` for the whole contig and spawns partition workers that feed a shared `colocated_map`. For the fused model, each window must run its *own* lookup over *its own* input variants, producing a self-contained `colocated_map`. The lookup is input-driven (the input VCF positions drive the `ColocatedSink`), so "lookup for a window" = run the same lookup plan with the window's batches as the input table.

### Task 3.1: Extract the lookup-plan construction into a reusable helper

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (the lookup-plan construction inside `prepare_contig_context`, the block that builds `plan` before the spawn loop at `:11798`)

- [ ] **Step 1: Identify the exact construction block**

Run: `grep -n "fn prepare_contig_context" datafusion/bio-function-vep/src/annotate_provider.rs` then read `:11532`–`:11798`. The block to extract is everything that constructs the lookup `plan: Arc<dyn ExecutionPlan>` and its `partition_coloc_sinks` from `chrom`, `config`, `cache`, and the registered ephemeral input table.

- [ ] **Step 2: Introduce the helper signature (no behavior change yet)**

Extract that block verbatim into a new function in the same module:

```rust
/// Build the variation/SIFT lookup plan and its per-partition colocated sinks
/// for a given registered input table. Pulled out of `prepare_contig_context`
/// so a per-window worker can reuse it with a window-scoped input table.
async fn build_lookup_plan(
    session: &SessionContext,
    cache: &PartitionedAnnotationCache,
    chrom: &str,
    config: &ContigAnnotationConfig,
    input_table_name: &str,
) -> Result<(Arc<dyn ExecutionPlan>, Vec<ColocatedSink>)> {
    // (moved verbatim from prepare_contig_context: build the lookup plan over
    //  `input_table_name` and the per-partition ColocatedSinks)
    todo!("paste the extracted construction block, returning (plan, partition_coloc_sinks)")
}
```

Then call it from `prepare_contig_context` where the block used to be. **This is a pure refactor — the `todo!` is a paste target, not a placeholder to leave in.**

- [ ] **Step 3: Verify no behavior change with the existing suite**

Run: `cargo test -p datafusion-bio-function-vep`
Expected: PASS — same tests, refactor is behavior-preserving.

- [ ] **Step 4: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "refactor(vep): extract build_lookup_plan from prepare_contig_context"
```

### Task 3.2: Add `lookup_window` to run lookup for one window's batches

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs`

- [ ] **Step 1: Write a failing integration test**

Add to the `mod tests` block in `annotate_provider.rs` (`:12468`). Use the smallest existing cache fixture used by other annotate tests (find one via `grep -n "PartitionedAnnotationCache" datafusion/bio-function-vep/src/annotate_provider.rs` in the test module; reuse that setup helper). The test runs the whole-contig lookup and a per-window lookup over the same variants and asserts the colocated maps are equal:

```rust
#[tokio::test(flavor = "multi_thread")]
async fn lookup_window_matches_whole_contig_lookup() {
    // Arrange: build a SessionContext + PartitionedAnnotationCache + a small
    // set of input variant batches for one chrom, using the existing test
    // fixture helper used elsewhere in this module.
    let fixture = build_small_annotation_fixture().await; // existing/added helper
    let chrom = fixture.chrom.clone();
    let input_batches = fixture.input_batches.clone();

    // Act: per-window lookup over exactly these batches.
    let window_map = lookup_window(
        &fixture.session,
        &fixture.cache,
        &chrom,
        &fixture.config,
        &input_batches,
    )
    .await
    .unwrap();

    // Assert: every input variant key present, values non-empty where the
    // whole-contig path produced colocated data. (Compare against the map the
    // existing per-contig path yields for the same variants.)
    let whole = fixture.whole_contig_colocated_map().await;
    assert_eq!(window_map, whole);
}
```

If `build_small_annotation_fixture` / `whole_contig_colocated_map` do not already exist, add them in the test module mirroring the setup in the nearest existing async annotate test (locate with `grep -n "#\[tokio::test" datafusion/bio-function-vep/src/annotate_provider.rs`).

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep lookup_window_matches_whole_contig_lookup`
Expected: FAIL — `lookup_window` not defined.

- [ ] **Step 3: Implement `lookup_window`**

```rust
/// Run the variation/SIFT lookup for one window's input batches and return the
/// self-contained colocated map for those variants. Registers the window's
/// batches as a transient MemTable, builds the lookup plan over it, executes
/// every partition to completion, and drains the colocated sinks.
async fn lookup_window(
    session: &SessionContext,
    cache: &PartitionedAnnotationCache,
    chrom: &str,
    config: &ContigAnnotationConfig,
    window_batches: &[RecordBatch],
) -> Result<HashMap<ColocatedKey, ColocatedData>> {
    use datafusion::datasource::MemTable;

    let schema = window_batches
        .first()
        .map(|b| b.schema())
        .ok_or_else(|| DataFusionError::Internal("empty window".into()))?;
    let table_name = format!("__vep_window_input_{chrom}");
    let mem = MemTable::try_new(schema, vec![window_batches.to_vec()])?;
    session.register_table(&table_name, Arc::new(mem))?;

    let (plan, sinks) = build_lookup_plan(session, cache, chrom, config, &table_name).await?;
    let task_ctx = session.task_ctx();

    let mut map: HashMap<ColocatedKey, ColocatedData> = HashMap::new();
    let mut target = Arc::new(std::mem::take(&mut map));
    for partition_id in 0..plan.output_partitioning().partition_count() {
        let mut stream = plan.execute(partition_id, Arc::clone(&task_ctx))?;
        while let Some(batch) = stream.next().await {
            let _ = batch?; // batch contents are captured via the sink
        }
        if let Some(sink) = sinks.get(partition_id) {
            let delta = drain_colocated_sink(sink)?;
            merge_colocated_delta(&mut target, delta);
        }
    }
    session.deregister_table(&table_name)?;
    Ok(Arc::try_unwrap(target).unwrap_or_else(|arc| (*arc).clone()))
}
```

(`merge_colocated_delta` is at `:2261`, `drain_colocated_sink` is the existing sink drainer used by `spawn_lookup_partition_worker` at `:9719`, `ColocatedKey`/`ColocatedData` are the existing types.)

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep lookup_window_matches_whole_contig_lookup`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): add lookup_window for self-contained per-window colocated lookup"
```

---

## Phase 4 — Fused `WindowWorker` (lookup → annotate → format for one window)

### Task 4.1: Add `annotate_window_fused` producing a window's output batches

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs`

- [ ] **Step 1: Write a failing test asserting parity with the serial path**

Add to the `annotate_provider.rs` test module:

```rust
#[tokio::test(flavor = "multi_thread")]
async fn annotate_window_fused_matches_serial_window() {
    // Arrange: shared context for one contig + a window of input batches.
    let fixture = build_small_annotation_fixture().await;
    let shared = fixture.shared_context().await; // Arc<SharedContigAnnotationContext>
    let window_batches = fixture.input_batches.clone();

    // Act: fused path (own lookup + own worker state).
    let fused = annotate_window_fused(
        Arc::clone(&shared),
        &fixture.session,
        &fixture.cache,
        &shared.config.clone(),
        &fixture.chrom,
        window_batches.clone(),
        None, // no projection
    )
    .await
    .unwrap();

    // Reference: existing serial path over the same window.
    let mut worker = AnnotationWorkerState::new(Arc::clone(&shared)).unwrap();
    worker.colocated_map = Arc::new(
        lookup_window(&fixture.session, &fixture.cache, &fixture.chrom, &shared.config, &window_batches)
            .await
            .unwrap(),
    );
    worker.lookup_done = true;
    let serial = annotate_worker_window(&mut worker, &window_batches, None).unwrap();

    // Assert: identical row counts and identical CSQ column bytes.
    let fused_rows: usize = fused.iter().map(|b| b.num_rows()).sum();
    let serial_rows: usize = serial.iter().map(|b| b.num_rows()).sum();
    assert_eq!(fused_rows, serial_rows);
    assert_eq!(concat_csq_strings(&fused), concat_csq_strings(&serial));
}
```

`concat_csq_strings` is a small test helper that concatenates the CSQ string column across batches (add it to the test module; the CSQ column index is `tmp_provider.vcf_field_count()`, see `:5639`).

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep annotate_window_fused_matches_serial_window`
Expected: FAIL — `annotate_window_fused` not defined.

- [ ] **Step 3: Implement `annotate_window_fused`**

```rust
/// Fused per-window pipeline: run the window's own lookup, populate a fresh
/// worker's colocated map, then annotate the window via the existing engine
/// path. Returns the window's annotated (and optionally projected) batches.
async fn annotate_window_fused(
    shared: Arc<SharedContigAnnotationContext>,
    session: &SessionContext,
    cache: &PartitionedAnnotationCache,
    config: &ContigAnnotationConfig,
    chrom: &str,
    window_batches: Vec<RecordBatch>,
    projection: Option<Vec<usize>>,
) -> Result<VecDeque<RecordBatch>> {
    // 1. Window-local lookup -> self-contained colocated map.
    let colocated = lookup_window(session, cache, chrom, config, &window_batches).await?;

    // 2. Fresh worker state seeded with the window's lookup result.
    let mut worker = AnnotationWorkerState::new(shared)?;
    worker.colocated_map = Arc::new(colocated);
    worker.lookup_done = true;

    // 3. Annotate via the existing window engine path (CPU-bound).
    annotate_worker_window(&mut worker, &window_batches, projection.as_deref())
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep annotate_window_fused_matches_serial_window`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): add fused per-window annotate (lookup+annotate) helper"
```

---

## Phase 5 — Work-stealing pool wired into the contig stream

**Context:** `StreamState::AnnotatingContig` (`annotate_provider.rs:11301`) currently drives a single worker. For `threads > 1` we instead: read the contig's input batches, feed `WindowPlanner`, dispatch windows to a bounded pool of `N` `tokio::spawn`ed `annotate_window_fused` tasks, and release results via `OrderedWindowDrain`. For `threads == 1` the existing single-worker path is kept (it *is* the fused path with N=1, lowest risk).

### Task 5.1: Add a parallel window driver behind `threads`

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs`

- [ ] **Step 1: Add a `threads` field to `ContigAnnotationConfig`**

In `ContigAnnotationConfig` (`:8671`) add:

```rust
    /// Number of parallel window workers within a contig (1 = serial fused path).
    threads: usize,
```

Set it in the construction site (`:5575`) to a temporary `1` (the real wiring lands in Phase 6):

```rust
    threads: 1,
```

- [ ] **Step 2: Write a failing parity test across worker counts**

```rust
#[tokio::test(flavor = "multi_thread")]
async fn parallel_windows_match_serial_for_contig() {
    let fixture = build_small_annotation_fixture().await; // multi-window input
    // Serial reference (threads = 1).
    let serial = run_contig_annotation(&fixture, 1).await;
    for threads in [2usize, 4] {
        let parallel = run_contig_annotation(&fixture, threads).await;
        assert_eq!(
            concat_csq_strings(&serial),
            concat_csq_strings(&parallel),
            "threads={threads} output diverged from serial",
        );
    }
}
```

`run_contig_annotation(&fixture, threads)` is a test helper that builds the contig stream end-to-end with `config.threads = threads` and collects all output batches (mirror the nearest existing end-to-end annotate test; locate with `grep -n "ContigAnnotationExec\|scan_with_transcript_engine" datafusion/bio-function-vep/src/annotate_provider.rs`).

- [ ] **Step 3: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep parallel_windows_match_serial_for_contig`
Expected: FAIL — parallel path not yet implemented (or helper missing).

- [ ] **Step 4: Implement the parallel driver**

Add a parallel branch used when `config.threads > 1`. Drive it as an owned async task set rather than threading new state through `StreamState` (simplest correct structure): collect the contig's input batches from the lookup-free input scan, plan windows, and run a bounded `tokio` pool.

```rust
/// Annotate one contig with `threads` parallel window workers, returning all
/// output batches in window order. Used when config.threads > 1.
async fn annotate_contig_parallel(
    shared: Arc<SharedContigAnnotationContext>,
    session: Arc<SessionContext>,
    cache: Arc<PartitionedAnnotationCache>,
    config: ContigAnnotationConfig,
    chrom: String,
    mut input_stream: SendableRecordBatchStream,
    projection: Option<Vec<usize>>,
) -> Result<VecDeque<RecordBatch>> {
    use tokio::sync::Semaphore;

    let mut planner = WindowPlanner::new(config.input_buffer_size.max(1));
    let mut drain = OrderedWindowDrain::new();
    let mut output: VecDeque<RecordBatch> = VecDeque::new();
    let permits = Arc::new(Semaphore::new(config.threads));
    let mut tasks: futures::stream::FuturesUnordered<_> = Default::default();

    let mut dispatch = |window: Window,
                        tasks: &mut futures::stream::FuturesUnordered<_>| {
        let permits = Arc::clone(&permits);
        let shared = Arc::clone(&shared);
        let session = Arc::clone(&session);
        let cache = Arc::clone(&cache);
        let config = config.clone();
        let chrom = chrom.clone();
        let projection = projection.clone();
        tasks.push(tokio::spawn(async move {
            let _permit = permits.acquire_owned().await.unwrap();
            let out = annotate_window_fused(
                shared, &session, &cache, &config, &chrom, window.batches, projection,
            )
            .await?;
            Ok::<_, DataFusionError>((window.index, out))
        }));
    };

    while let Some(batch) = input_stream.next().await {
        if let Some(window) = planner.push(batch?) {
            dispatch(window, &mut tasks);
        }
        // Drain finished tasks opportunistically to bound memory.
        while let Some(done) = tasks.try_next().now_or_never().flatten() {
            let (index, out) = done.map_err(|e| DataFusionError::External(Box::new(e)))??;
            output.extend(drain.complete(index, out.into_iter().collect()));
        }
    }
    if let Some(window) = planner.finish() {
        dispatch(window, &mut tasks);
    }
    while let Some(done) = tasks.next().await {
        let (index, out) = done.map_err(|e| DataFusionError::External(Box::new(e)))??;
        output.extend(drain.complete(index, out.into_iter().collect()));
    }
    Ok(output)
}
```

Wire it into the stream: in the `PreparingContig`→`AnnotatingContig` transition (`:11282`), when `config.threads > 1`, run `annotate_contig_parallel` (block the stream on it via a future stored in a new `StreamState` variant `AnnotatingContigParallel(Pin<Box<dyn Future<...>>>)`, mirroring `PrepareFuture`), then move its `VecDeque` into a `DrainingWindow` state for emission. When `config.threads == 1`, keep the existing serial `AnnotatingContig` path unchanged.

Note: `input_stream` here is the contig's raw input scan (the table the lookup is built over), obtained the same way `prepare_contig_context` registers the input — reuse that registration; the windows drive their own lookups, so the parallel path does **not** spawn the per-contig lookup partition workers.

- [ ] **Step 5: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep parallel_windows_match_serial_for_contig`
Expected: PASS for threads ∈ {1, 2, 4}.

- [ ] **Step 6: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): parallel window pool with ordered drain for threads>1"
```

---

## Phase 6 — `threads` knob, clean break, caller updates

### Task 6.1: Replace the option parser with a single `threads` key

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`:5487`–`:5516`)

- [ ] **Step 1: Write failing tests for the new parser + fail-loud**

Add to the `annotate_provider.rs` test module:

```rust
#[test]
fn parses_threads_option() {
    assert_eq!(parse_threads_option(r#"{"threads": 8}"#), 8);
    assert_eq!(parse_threads_option(r#"{}"#), 1); // default
    assert_eq!(parse_threads_option(r#"{"threads": 0}"#), 1); // clamp to >=1
}

#[test]
fn rejects_removed_options() {
    for key in ["forks", "target_partitions", "contig_parallelism", "inline_lookup"] {
        let json = format!(r#"{{"{key}": 4}}"#);
        let err = check_no_removed_options(&json).unwrap_err();
        assert!(err.to_string().contains(key), "error should name {key}");
        assert!(err.to_string().contains("threads"), "error should point to threads");
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p datafusion-bio-function-vep parses_threads_option rejects_removed_options`
Expected: FAIL — functions not defined.

- [ ] **Step 3: Implement parser + fail-loud, delete old parsing**

Add:

```rust
fn parse_threads_option(options_json: &str) -> usize {
    Self::parse_json_i64_option(options_json, "threads")
        .and_then(|v| usize::try_from(v).ok())
        .filter(|v| *v > 0)
        .unwrap_or(1)
}

fn check_no_removed_options(options_json: &str) -> Result<()> {
    for key in ["forks", "target_partitions", "contig_parallelism", "inline_lookup"] {
        if options_json.contains(&format!("\"{key}\"")) {
            return Err(DataFusionError::Plan(format!(
                "VEP option `{key}` was removed; use `threads` instead"
            )));
        }
    }
    Ok(())
}
```

Replace the parsing block at `:5487`–`:5516` (the `worker_forks`/`target_partitions`/`chromosome_lanes`/`inline_lookup` lets) with:

```rust
    if let Some(opts) = self.options_json.as_deref() {
        Self::check_no_removed_options(opts)?;
    }
    let threads = self
        .options_json
        .as_deref()
        .map(Self::parse_threads_option)
        .unwrap_or(1);
```

Update the `ContigAnnotationConfig` construction (`:5575`) to set `threads,` and remove the now-dead `target_partitions`/`chromosome_lanes`/`inline_lookup` fields **from the struct and all its uses**. For the lookup-partition fallback (`inline_lookup`-gated branches at `:11798`/`:11831`): the parallel path no longer uses them; keep the serial `threads==1` path running the existing single inline lookup (equivalent to the old `inline_lookup=true`).

- [ ] **Step 4: Run the full suite**

Run: `cargo test -p datafusion-bio-function-vep`
Expected: PASS (after fixing all references to removed fields the compiler flags).

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): single threads knob, fail loud on removed options (clean break)"
```

### Task 6.2: Update `AnnotateVcfConfig` and all in-repo callers

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (`:670`–`:743`)
- Modify: examples that set old knobs — `examples/bench_annotate_debug.rs`, `examples/lookup_parquet_bench.rs`, `examples/load_cache.rs`, and any that build the options JSON
- Modify: `scripts/run_annotation_fast.py` and siblings under `scripts/` that pass these options

- [ ] **Step 1: Update `AnnotateVcfConfig`**

In `vcf_sink.rs` remove the `forks`, `workers`, and `target_partitions` fields and add:

```rust
    /// Number of parallel window workers per contig (1 = serial).
    pub threads: usize,
```

Update wherever the options JSON is assembled from `AnnotateVcfConfig` to emit `"threads": <n>` and stop emitting `forks`/`target_partitions`/`contig_parallelism`/`inline_lookup`. Find with `grep -rn "forks\|target_partitions\|contig_parallelism\|inline_lookup" datafusion/bio-function-vep/src`.

- [ ] **Step 2: Update example/bench callers**

For each example flagged by `grep -rn "with_target_partitions\|\"forks\"\|target_partitions\|inline_lookup" datafusion/bio-function-vep/examples`, replace the old option/CLI arg with a `threads` arg. (`with_target_partitions` on `SessionConfig` is a DataFusion call unrelated to our knob — leave those that configure DataFusion's own runtime, only change places that fed our VEP options JSON or `AnnotateVcfConfig`.)

- [ ] **Step 3: Build all examples**

Run: `cargo build -p datafusion-bio-function-vep --examples`
Expected: PASS — no example references removed fields.

- [ ] **Step 4: Update Python scripts**

Edit `scripts/run_annotation_fast.py` (and siblings) to pass `threads=N` and drop the old keys. Run a smoke invocation if a fixture is available (otherwise grep-verify no old keys remain): `grep -rn "forks\|target_partitions\|inline_lookup" datafusion/bio-function-vep/scripts` → expect no option-passing hits.

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add -A datafusion/bio-function-vep
git commit -m "feat(vep): migrate AnnotateVcfConfig + callers to threads knob"
```

---

## Phase 7 — Parity gate + scaling bench

### Task 7.1: chr1 end-to-end parity across thread counts

**Files:**
- Create/Modify: `datafusion/bio-function-vep/examples/bench_parallel_windows.rs` (or extend an existing chr1 bench)

- [ ] **Step 1: Add a parity harness**

Write an example that annotates a fixed chr1 input (use the cache/input the existing chr1 benches use; find via `grep -rln "chr1" datafusion/bio-function-vep/examples`) at `threads ∈ {1,2,4,8}`, writes each output VCF, and diffs the CSQ columns against the `threads=1` run, asserting zero differences. Print PASS/FAIL.

- [ ] **Step 2: Run the parity harness**

Run: `cargo run --release -p datafusion-bio-function-vep --example bench_parallel_windows -- <chr1_input> <cache_root>`
Expected: `PASS` — all thread counts byte-identical to serial on the CSQ column.

- [ ] **Step 3: Record the scaling curve**

Extend the harness to print wall-clock per thread count. Run it and capture the chr1 scaling curve (`threads → seconds`), comparing the `threads=1` number to the engine-bound baseline in `docs/superpowers/plans/2026-06-17-vep-engine-perf-handoff.md`.

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/examples/bench_parallel_windows.rs
git commit -m "test(vep): chr1 parity + scaling bench for parallel window annotation"
```

### Task 7.2: Update memory + spec status

- [ ] **Step 1: Record results in project memory**

Add a memory file noting the chr1 scaling numbers and that parallel within-contig annotation shipped via the `threads` knob; link `[[vep-annotation-bottlenecks]]`. Update `MEMORY.md` index.

- [ ] **Step 2: Commit**

```bash
git add -A
git commit -m "docs(vep): record parallel-window scaling results in memory"
```

---

## Downstream coordination (separate repo — vepyr / polars-bio)

Not part of this repo's tasks, but required to ship (per spec §7.1, clean break):

1. After this plan lands, bump the pinned `rev` of `datafusion-bio-functions` in `polars-bio` / vepyr.
2. Replace the Python parallelism parameter(s) with a single `threads=N`, emitting only `"threads"` in the options JSON.
3. Because it is a clean break, the `rev` bump and the Python knob change must ship together — the engine now errors on the removed keys.

---

## Self-Review notes

- **Spec coverage:** §2 core model → Phases 4–5; §3 independence → exploited by per-window self-contained `colocated_map` (Phase 3–4); §5 geometry cache → no task (already `OnceLock`, confirmed); §6 fused single path → Phases 4–5 (threads=1 is the same fused helper); §7 threads budget + §7.1 clean break → Phase 6; §8 ordering/backpressure → `OrderedWindowDrain` (Phase 2) + `Semaphore` bound (Phase 5); §9 tests → Phases 4,5,7; §10 out-of-scope (haplotype, global queue) → not implemented by design.
- **Known approximations (require reading the named code at execution time, not placeholders to ship):** Task 3.1 `build_lookup_plan` is a verbatim extraction of an existing block; Task 5.4 reuses the existing contig input-scan registration. Both name exact functions/lines to mirror; the executing worker fills the extracted body from the cited source. Test helpers (`build_small_annotation_fixture`, `run_contig_annotation`, `concat_csq_strings`) are to be modelled on the nearest existing annotate test (grep anchors given).
- **Type consistency:** `Window`/`WindowPlanner` (Phase 1), `OrderedWindowDrain::complete` (Phase 2), `lookup_window` (Phase 3), `annotate_window_fused` (Phase 4), `annotate_contig_parallel` (Phase 5), `parse_threads_option`/`check_no_removed_options` (Phase 6) names are used consistently across phases.
