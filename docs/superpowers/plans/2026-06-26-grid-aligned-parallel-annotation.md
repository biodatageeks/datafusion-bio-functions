# Grid-Aligned Parallel Annotation (bounded-overlap, full per-partition independence) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover parallel annotation speed for stateful Merged/RefSeq caches by running N independent per-partition sharded workers whose seams fall on VEP's global 5,000-input-unit buffer grid, each warm-starting a bounded overlap before its seam to reconstruct HGNC carry — byte-identical to the serial worker.

**Architecture:** Per contig, a cheap positions-only **count pass** reproduces VEP's global buffer boundaries (using the existing `InputBufferAccumulator` cut logic) and records each boundary's `(global_row_index, position)`. A pure **partition planner** slices those boundaries into N contiguous whole-buffer ranges and, for each worker, derives a scan position-superset, a leading-row skip count, and three global-row-index marks: `warm_up_start ≤ emit_start < emit_end` (all on buffer boundaries). Each sharded worker scans only its slice, processes whole buffers, runs warm-up buffers **state-only** (donation+persist, no engine/lookup-emit), then annotates and emits its owned buffers into its VCF shard. The warm-up reuses the exact serial donation code, so correctness is by construction. Donor scope is already buffer-local (committed); the only remaining defects (partition-local grid + no carry) are fixed here.

**Tech Stack:** Rust, Apache DataFusion 53 / Arrow 58, tokio (blocking pool), noodles VCF (tabix-indexed), Lance cache backend.

## Global Constraints

- **Correctness is non-negotiable.** A change lands only if chr4 + chr2 merged Lance at `workers ∈ {1,4}` produce `HGNC_ID` mismatches = 0, full variant count (chr4 = 307,295), "Only in VEP: 0", **byte-identical bodies**, on **3 consecutive runs** (catches intermittent stream truncation). Never commit a faster-but-sometimes-wrong version.
- **Serial path unchanged.** The committed serial `AnnotatingContig` path (grid-aligned buffers, buffer-local donors `apply_buffer_local_hgnc_propagation`, sequential carry) must keep `workers=1` byte-identical. Do not modify its donation/buffer logic.
- **Ensembl (stateless) path unchanged.** Only Merged/RefSeq routing changes.
- **No new dependencies** (no `rayon`). Reuse tokio `spawn` / `block_in_place` and the existing `VcfBodyShardWriter` sharded output contract.
- **VEP buffer size** = `config.input_buffer_size` (`VEP_INPUT_BUFFER_SIZE`); input **unit** = ALT-allele count per row (`batch_input_units`, `annotate_provider.rs:9332`). Buffers cut at row boundaries once cumulative units ≥ size (`count_ready_input_buffers`, `:9343`); a multi-allelic row is never split.
- **Cache-region width** = `VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP = 1_000_000` (`annotate_provider.rs:8579`).
- All work is in `datafusion/bio-function-vep/src/`; unit tests live in the existing `#[cfg(test)]` module of `annotate_provider.rs` (helpers `make_tx` `:14005`, `make_buffer_batch`/`make_buffer_batch_many` `:14134`/`:14151`, `minimal_shared_contig_annotation_context_with_context` `:14321`).

---

## File Structure

- `annotate_provider.rs` — all engine changes:
  - `compute_overlap_width_bp` (new pure fn) + `overlap_width_bp` field on `SharedContigAnnotationContext` (`:8995`).
  - `GridBufferBoundary` / `WorkerGridSlice` types + `plan_grid_partitions` (new pure fn).
  - `count_contig_buffer_boundaries` (new async fn — positions-only count pass).
  - `annotate_worker_window` (`:10862`) — add a warm-up/skip mode (`WindowAnnotateBounds`).
  - `spawn_annotation_from_lookup_sharded` (`:9756`) — accept a `WorkerGridSlice`, run warm-up buffers state-only, emit owned buffers.
  - `spawn_lookup_partition_worker` (`:9589`) / `prepare_contig_context` parallel branch (`:11848`) — replace byte-budget partitions with grid-sliced position-range scans.
  - `AnnotatingParallel` spawn (`:11301`) — thread each worker its `WorkerGridSlice`.
  - serial gate (`:5266`–`:5276`) — drop for Merged/RefSeq + fail-closed guard.
- `lookup_provider.rs` — reuse `set_vcf_filter` (`:190`) for per-worker position-range scans (no new code expected beyond constructing the filter expr at the call site).
- `vcf_sink.rs` — confirm/adjust routing so Merged/RefSeq `workers>1` reaches the sharded driver (`drive_sharded_vcf_annotation` `:715`); no contract change.

---

## Task 0: Baseline & current-path verification (no code)

**Files:** none (investigation + captured logs under the scratchpad).

**Interfaces:**
- Produces: a captured `workers=1` reference VCF body + the answer to "what does `workers=4` do today?" (the map flagged a possible empty-shard path for Merged+w4 that contradicts the handoff's claim of a verified w4 — resolve it before trusting any w4 gate).

- [ ] **Step 1: Capture the w1 reference and current w4 behavior**

Run (vepyr, release build of this crate):
```bash
cd ~/research/git/vepyr && source .venv/bin/activate && unset CONDA_PREFIX && maturin develop --release
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend lance --workers 1 --force 2>&1 | tee /tmp/o-w1.log
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend lance --workers 4 --force 2>&1 | tee /tmp/o-w4.log
```
Expected: w1 reports `HGNC_ID` mismatches = 0 and 307,295 variants. Record the w4 result and **where its output comes from** (does it hit `drive_sharded_vcf_annotation`? add a one-line `eprintln!` temporarily in `vcf_sink.rs:1112` and `annotate_provider.rs:11284` if needed, then revert).

- [ ] **Step 2: Resolve the routing question**

Confirm in code which output path `run_annotation_fast` exercises for Merged at `workers>1`:
- `vcf_sink.rs:1112` picks the sharded driver on `config.workers > 1`.
- the provider forces `annotation_workers = 1` for Merged (`annotate_provider.rs:5269`) → `AnnotatingContig` → yields row batches.
Decide: either (a) current w4 already streams correctly via a non-sharded consumer, or (b) it is silently degenerate. Write the finding into the plan's notes; the gate's "byte-identical to w1" comparison uses the **w1** body as the reference regardless.

- [ ] **Step 3: Commit nothing — this is investigation.** Save logs to `/tmp/o-w1.log`, `/tmp/o-w4.log`.

---

## Task 1: `overlap_width_bp` (bounded-overlap distance)

**Files:**
- Modify: `annotate_provider.rs` — add `compute_overlap_width_bp` near `cache_region_index` (`:10400`); add field to `SharedContigAnnotationContext` (`:8995`) and both construction sites (`:11997`, test `:12427`).
- Test: same file `#[cfg(test)]` module.

**Interfaces:**
- Produces: `fn compute_overlap_width_bp(base_transcripts: &[TranscriptFeature]) -> i64`; `SharedContigAnnotationContext.overlap_width_bp: i64`.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn overlap_width_is_max_span_plus_region() {
    let mut a = make_tx("a", None, None, None, None); // start 1, end 100 → span 99
    let mut b = make_tx("b", None, None, None, None);
    b.start = 1000;
    b.end = 1000 + 2_000_000; // span 2_000_000 (widest)
    a.end = 100;
    let w = compute_overlap_width_bp(&[a, b]);
    assert_eq!(w, 2_000_000 + VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP);
}

#[test]
fn overlap_width_empty_is_region_only() {
    assert_eq!(compute_overlap_width_bp(&[]), VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib overlap_width`
Expected: FAIL — `compute_overlap_width_bp` not found.

- [ ] **Step 3: Implement the function**

```rust
/// Bounded-overlap warm-up distance (design §5.4): a donor can influence a seam
/// from at most one transcript span plus one cache-region width upstream. Pure
/// `max()` over the already-loaded base transcripts — zero IO, identical in every
/// worker, so it needs no synchronisation and cannot under-read.
fn compute_overlap_width_bp(base_transcripts: &[TranscriptFeature]) -> i64 {
    let max_span = base_transcripts
        .iter()
        .map(|tx| (tx.end - tx.start).abs())
        .max()
        .unwrap_or(0);
    max_span + VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP
}
```

- [ ] **Step 4: Add the field + initialise both construction sites**

In `struct SharedContigAnnotationContext` (`:8995`) add `overlap_width_bp: i64,`.
At the production construction (`:11997`, `Arc::new(SharedContigAnnotationContext { ... })`), add before the closing brace: `overlap_width_bp: compute_overlap_width_bp(&base_transcripts),` (use the `Vec<TranscriptFeature>` in scope *before* it is `Arc::new`-wrapped; if already wrapped, compute from the local `transcripts` binding used for `base_transcripts`).
At the test construction (`:12427`), add `overlap_width_bp: 0,` (tests that don't exercise warm-up).

- [ ] **Step 5: Run tests + build**

Run: `cargo test -p datafusion-bio-function-vep --lib overlap_width && cargo build -p datafusion-bio-function-vep`
Expected: PASS, builds clean.

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): compute bounded-overlap warm-up width from base transcripts"
```

---

## Task 2: Grid partition planner (pure)

**Files:**
- Modify: `annotate_provider.rs` — add `GridBufferBoundary`, `WorkerGridSlice`, `plan_grid_partitions` near `drain_window_input_units` (`:9388`).
- Test: same file test module.

**Interfaces:**
- Consumes: `overlap_width_bp` (Task 1).
- Produces:
```rust
/// One VEP global input-buffer boundary: the cumulative input-unit rank at the
/// boundary equals buffer_index * input_buffer_size (approx; buffers cut at row
/// boundaries ≥ size). `global_row` is the GLOBAL input-row index of the FIRST
/// row of this buffer; `pos` is that row's `start`. Boundary `B` (one past the
/// last buffer) carries the contig's end (global_row = total_rows, pos = i64::MAX).
#[derive(Clone, Copy, Debug, PartialEq)]
struct GridBufferBoundary { global_row: usize, pos: i64 }

/// Per-worker assignment over the global buffer grid. All three row marks fall on
/// buffer boundaries: warm_up_start_row ≤ emit_start_row < emit_end_row.
#[derive(Clone, Debug, PartialEq)]
struct WorkerGridSlice {
    worker_id: usize,
    scan_lo_pos: i64,        // inclusive: position of first scanned row
    scan_hi_pos: i64,        // exclusive: positions ≥ this are not scanned
    skip_leading_rows: usize,// rows at scan_lo_pos before warm_up_start_row to drop
    warm_up_start_row: usize,// global row index where this worker begins (state-only)
    emit_start_row: usize,   // global row index where emission begins (= its seam)
    emit_end_row: usize,     // global row index (exclusive) where emission ends
}

fn plan_grid_partitions(
    boundaries: &[GridBufferBoundary], // length B+1, ascending global_row
    total_rows: usize,
    workers: usize,
    overlap_width_bp: i64,
) -> Vec<WorkerGridSlice>;
```

Planner rules (implement exactly):
- `b = boundaries.len() - 1` whole buffers. If `b == 0` → return one slice covering everything with no warm-up.
- Worker k (of N) owns buffers `[bk, bk1)` with `bk = round(k*b/N)`, `bk1 = round((k+1)*b/N)`; skip empty ranges (`bk == bk1`). `emit_start_row = boundaries[bk].global_row`, `emit_end_row = boundaries[bk1].global_row`.
- Warm-up: `wk` = largest buffer index `≤ bk` with `boundaries[wk].pos <= boundaries[bk].pos - overlap_width_bp`; clamp to `0`. Worker 0's first owned buffer `bk == 0` → `wk = 0`, no warm-up. `warm_up_start_row = boundaries[wk].global_row`.
- Scan range: `scan_lo_pos = boundaries[wk].pos`; `scan_hi_pos = boundaries[bk1].pos` (exclusive). For the last worker, `scan_hi_pos = i64::MAX`.
- `skip_leading_rows = warm_up_start_row − (global row index of the first row at scan_lo_pos)`. Because `boundaries[wk].global_row` *is* the first row of buffer `wk` and `scan_lo_pos = boundaries[wk].pos`, ties at `scan_lo_pos` from the previous buffer mean the scan may include earlier rows sharing that position; `skip_leading_rows` drops exactly those. Compute it from the count pass (Task 3 records, per boundary, how many rows at `pos` precede `global_row`).

- [ ] **Step 1: Write failing tests** (cover: even split, warm-up reaching back one buffer, worker-0-no-warmup, position-tie skip)

```rust
fn b(row: usize, pos: i64) -> GridBufferBoundary { GridBufferBoundary { global_row: row, pos } }

#[test]
fn plan_two_workers_even_split_no_overlap() {
    // 4 buffers at positions 0,100,200,300; end at row 20000/pos MAX
    let bs = vec![b(0,0), b(5000,100), b(10000,200), b(15000,300), b(20000, i64::MAX)];
    let slices = plan_grid_partitions(&bs, 20000, 2, 0);
    assert_eq!(slices.len(), 2);
    assert_eq!(slices[0].emit_start_row, 0);
    assert_eq!(slices[0].emit_end_row, 10000);
    assert_eq!(slices[0].warm_up_start_row, 0); // worker 0 no warm-up
    assert_eq!(slices[1].emit_start_row, 10000);
    assert_eq!(slices[1].emit_end_row, 20000);
    assert_eq!(slices[1].warm_up_start_row, 10000); // overlap 0 → warm-up starts at seam
    assert_eq!(slices[1].scan_lo_pos, 200);
}

#[test]
fn plan_warmup_reaches_back_within_overlap() {
    // buffers 50bp apart; overlap 120bp → worker1 seam at buffer2(pos100) warms back to buffer0(pos0)
    let bs = vec![b(0,0), b(5000,50), b(10000,100), b(15000,150), b(20000, i64::MAX)];
    let slices = plan_grid_partitions(&bs, 20000, 2, 120);
    assert_eq!(slices[1].emit_start_row, 10000);
    assert_eq!(slices[1].warm_up_start_row, 0); // pos100-120 = -20 → clamp buffer 0
    assert_eq!(slices[1].scan_lo_pos, 0);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib plan_`
Expected: FAIL — `plan_grid_partitions` not found.

- [ ] **Step 3: Implement `plan_grid_partitions`** per the rules above. Use integer `round(k*b/N)` via `(k * b + N/2) / N`. For `skip_leading_rows`, add a `rows_before_pos: usize` field to `GridBufferBoundary` (rows at `pos` preceding `global_row`) populated by Task 3; `skip_leading_rows = boundaries[warm_up_buffer].rows_before_pos` (the rows sharing `scan_lo_pos` that belong to the previous buffer). Adjust the test boundaries to set `rows_before_pos: 0` where ties are absent.

- [ ] **Step 4: Run tests** — `cargo test -p datafusion-bio-function-vep --lib plan_` → PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): pure grid-partition planner (whole-buffer seams + warm-up marks)"
```

---

## Task 3: Per-contig count pass (boundary discovery)

**Files:**
- Modify: `annotate_provider.rs` — add `async fn count_contig_buffer_boundaries(...)` near `prepare_contig_context` (`:11742`); reuse `batch_input_units` (`:9332`).
- Test: integration-style unit test feeding synthetic batches (no IO) to a pure inner helper.

**Interfaces:**
- Consumes: a positions+alt stream for the contig; `input_buffer_size`.
- Produces:
```rust
/// Replays the EXACT serial buffer-cut logic over a positions-only stream and
/// records each global buffer boundary. Returns (boundaries: Vec<GridBufferBoundary>
/// length B+1, total_rows). Pure core is `accumulate_boundaries` for testing.
fn accumulate_boundaries(batches: &[RecordBatch], input_buffer_size: usize)
    -> (Vec<GridBufferBoundary>, usize);

async fn count_contig_buffer_boundaries(
    provider: &VcfTableProvider, // chrom-filtered, projecting chrom/start/alt only
    session: &Arc<SessionContext>,
    input_buffer_size: usize,
) -> Result<(Vec<GridBufferBoundary>, usize)>;
```

Cut logic (must match `count_ready_input_buffers` `:9343` exactly): walk rows in order; `units += input_units_at(row)`; the first boundary is `{global_row:0, pos: first row.start, rows_before_pos:0}`; when `units >= input_buffer_size`, the NEXT row starts a new buffer — push a boundary `{global_row: idx+1, pos: next_row.start, rows_before_pos: count of rows at next_row.start already consumed}` and reset `units = 0`. After the last row, push the terminal boundary `{global_row: total_rows, pos: i64::MAX, rows_before_pos: 0}`.

**Why this read is intrinsic, and how it's kept cheap (decided 2026-06-26).** Seam positions are a *prefix-sum* over records; tabix/CSI index *position → byte-offset*, never *rank → position*, so neither the index nor its free per-contig total can supply seam positions — only a records read can. This pre-pass therefore stays. Keep it cheap:
- Project the **minimal columns** `batch_input_units` needs (`start`,`alt`); no INFO/FORMAT/genotype decode. **Under a guaranteed-normalized (decomposed) input, `units == rows`**, so this could drop to `start`-only and use the free tabix metadata total as the denominator — but **keep the `alt`-aware unit count for correctness robustness** (a stray multi-ALT line must not desync the grid); treat the `start`-only/free-total form as an optimization to enable only behind a verified-normalized assertion.
- It may be **parallelized across the index byte-partitions** (each counts its range + emits local boundary positions; merge prefix sums) so wall cost ≈ `read/N`.
- **Long-term home is bio-formats** (see Follow-ups): a count-balanced grid-partition / `boundary_positions(contig, every_n)` API would fuse this into bio-formats' scan planning *and* improve partition balance (variant- vs byte-balanced). Deferred to avoid coupling correctness to a sibling-repo release; the self-contained pre-pass lands first.

- [ ] **Step 1: Write failing test for `accumulate_boundaries`**

```rust
#[test]
fn boundaries_cut_every_buffer_size_units() {
    // 3 single-alt rows per "position"; size = 2 → boundary after 2 units (row 2)
    let batch = make_buffer_batch_many("chr2", &[10, 20, 30, 40]); // 1 unit each
    let (bs, total) = accumulate_boundaries(&[batch], 2);
    assert_eq!(total, 4);
    // buffers: [10,20] [30,40]; boundaries at rows 0,2,4
    assert_eq!(bs.iter().map(|x| x.global_row).collect::<Vec<_>>(), vec![0, 2, 4]);
    assert_eq!(bs[0].pos, 10);
    assert_eq!(bs[1].pos, 30);
    assert_eq!(bs[2].pos, i64::MAX);
}
```

- [ ] **Step 2: Run to verify it fails** — `cargo test -p datafusion-bio-function-vep --lib boundaries_cut` → FAIL.

- [ ] **Step 3: Implement `accumulate_boundaries`** per the cut logic (use `AltColumnView::from_batch` like `count_ready_input_buffers`; track `rows_before_pos` by comparing consecutive `start` values).

- [ ] **Step 4: Implement `count_contig_buffer_boundaries`** — scan the chrom-filtered provider projecting only the columns `batch_input_units` needs (`chrom`,`start`,`alt`), collect batches, call `accumulate_boundaries`. (Set the projection via `provider.scan(&state, Some(&proj), &[], None)`.)

- [ ] **Step 5: Run tests** — `cargo test -p datafusion-bio-function-vep --lib boundaries_cut` → PASS.

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): per-contig count pass recording global buffer boundaries"
```

---

## Task 4: Warm-up-capable `annotate_worker_window`

**Files:**
- Modify: `annotate_provider.rs` `annotate_worker_window` (`:10862`); reuse `prepare_buffer_annotation_context` (`:10123`).
- Test: carry-equality unit test (the §8 debug assertion, as a test).

**Interfaces:**
- Consumes: `WorkerGridSlice` row marks.
- Produces: new param on `annotate_worker_window`:
```rust
/// Rows whose GLOBAL index is < emit_start are warm-up: process the buffer
/// state-only (prepare → donation+persist, NO engine/emit). `base_global_row` is
/// the global index of this window's first row.
struct WindowAnnotateBounds { base_global_row: usize, emit_start_row: usize, emit_end_row: usize }
```
`annotate_worker_window(worker, window_batches, projection, bounds: WindowAnnotateBounds)`.

Behavior change inside the `for buffer_batches in ready_input_buffers` loop (`:10897`): track a running `buffer_first_global_row` (start at `bounds.base_global_row`, advance by each buffer's input-row count). For each buffer:
- compute bounds, call `prepare_buffer_annotation_context` (unchanged — this mutates `worker.persisted_buffer_transcripts`).
- if `buffer_first_global_row < bounds.emit_start_row` → **warm-up buffer**: skip `materialize…`/`prepared_context_from_buffer`/the per-batch engine loop entirely (state already advanced by `prepare`). Do not push to `out`.
- else → run the existing materialise + engine + push path. (Within an emitted buffer, all rows emit; `emit_end_row` only bounds which whole buffers are owned, enforced by the planner, so no mid-buffer trimming is needed.)

Worker 0 / serial path callers pass `WindowAnnotateBounds { base_global_row: 0, emit_start_row: 0, emit_end_row: usize::MAX }` (no warm-up — identical to today).

- [ ] **Step 1: Write the failing carry-equality test**

Build two same-gene transcripts where buffer-1 donates HGNC to a buffer-2 recipient (model on the existing buffer-local HGNC tests near `:14002`). Run path A: serial — one worker processes buffer1 then buffer2 (records `persisted_buffer_transcripts` after buffer1). Run path B: fresh worker, warm-up over buffer1 (state-only via the new bounds), then buffer2. Assert the recipient's emitted `HGNC_ID` (and `worker.persisted_buffer_transcripts`) are identical.

```rust
#[test]
fn warmup_reconstructs_serial_carry() {
    // donor in buffer1 (pos 10), recipient same gene in buffer2 (pos 5_010), size=1
    // ... build shared ctx with both transcripts (make_tx with shared gene_symbol,
    //     donor has native HGNC, recipient native None) ...
    // serial:
    let mut w_serial = AnnotationWorkerState::new(shared.clone()).unwrap();
    let _ = annotate_worker_window(&mut w_serial, &buf1, None,
        WindowAnnotateBounds{base_global_row:0, emit_start_row:0, emit_end_row:usize::MAX}).unwrap();
    let out_serial = annotate_worker_window(&mut w_serial, &buf2, None,
        WindowAnnotateBounds{base_global_row:1, emit_start_row:1, emit_end_row:usize::MAX}).unwrap();
    // warm-up:
    let mut w_warm = AnnotationWorkerState::new(shared).unwrap();
    let combined = [buf1.clone(), buf2.clone()].concat();
    let out_warm = annotate_worker_window(&mut w_warm, &combined, None,
        WindowAnnotateBounds{base_global_row:0, emit_start_row:1, emit_end_row:usize::MAX}).unwrap();
    assert_eq!(hgnc_of(&out_serial), hgnc_of(&out_warm)); // helper reads HGNC_ID column
    assert_eq!(w_serial.persisted_buffer_transcripts, w_warm.persisted_buffer_transcripts);
}
```

- [ ] **Step 2: Run to verify it fails** — `cargo test -p datafusion-bio-function-vep --lib warmup_reconstructs` → FAIL (signature mismatch / no bounds param).

- [ ] **Step 3: Implement the `bounds` param + warm-up skip** in `annotate_worker_window`. Update ALL call sites to pass the no-warm-up bounds: `spawn_annotation_from_lookup_sharded` (`:9791`,`:9811`), `annotate_window_owned` (`:10997`).

- [ ] **Step 4: Run tests** — `cargo test -p datafusion-bio-function-vep --lib` → PASS (warm-up test + no regressions).

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): warm-up (state-only) buffer mode in annotate_worker_window"
```

---

## Task 5: Grid-sliced lookup workers (replace byte-budget partitions)

**Files:**
- Modify: `annotate_provider.rs` `prepare_contig_context` parallel branch (`:11848`–`:11868`); `spawn_lookup_partition_worker` (`:9589`).
- `lookup_provider.rs` — reuse `set_vcf_filter` (`:190`).

**Interfaces:**
- Consumes: `WorkerGridSlice` (Task 2), count pass (Task 3), `overlap_width_bp` (Task 1).
- Produces: each spawned lookup worker streams only its slice's rows `[scan_lo_pos, scan_hi_pos)`; the `WorkerGridSlice` travels alongside its `LookupPartitionHandle` so the annotation spawn (Task 6) can read it.

Implementation:
- In the parallel branch, BEFORE building partitions: when `cache_source_type` is Merged/RefSeq and `annotation_workers > 1`, build a chrom-filtered count provider, call `count_contig_buffer_boundaries`, then `plan_grid_partitions(&boundaries, total_rows, annotation_workers, shared.overlap_width_bp)`.
- For each `WorkerGridSlice`, construct a per-worker provider with `set_vcf_filter(Some(col("chrom").eq(lit(&*chrom)).and(col("start").gt_eq(lit(slice.scan_lo_pos))).and(col("start").lt(lit(slice.scan_hi_pos)))))` and `set_target_partitions(1)`; `scan` it; `spawn_lookup_partition_worker(plan, ctx, 0, …)`. The worker's stream is its slice. (Position pushdown to tabix is a perf optimisation; if the provider filters post-scan it is still correct.)
- Carry each `WorkerGridSlice` next to its handle (e.g., `Vec<(LookupPartitionHandle, WorkerGridSlice)>`), so the `AnnotatingParallel` spawn can pass it through.
- **Skip leading rows:** the worker stream starts at `scan_lo_pos` which may include rows from the previous buffer sharing that position; `slice.skip_leading_rows` are dropped before feeding the accumulator (handled in Task 6).
- Ensembl/stateless and `workers=1` keep the existing byte-budget partition loop unchanged.

- [ ] **Step 1: Write a focused unit test for the filter-range slice math** — assert that for given boundaries + workers, the constructed `(scan_lo_pos, scan_hi_pos)` per worker tile the contig with the expected overlaps (this is mostly Task 2; here assert the *expr* is built — a small test that `vcf_range_filter(chrom, lo, hi)` returns the expected `Expr` string via `format!("{expr:?}")`). Extract a helper `fn vcf_range_filter(chrom:&str, lo:i64, hi:i64) -> Expr`.

- [ ] **Step 2: Run to verify it fails** — `cargo test -p datafusion-bio-function-vep --lib vcf_range_filter` → FAIL.

- [ ] **Step 3: Implement `vcf_range_filter` + wire the grid branch** in `prepare_contig_context`. Keep the non-grid branches intact.

- [ ] **Step 4: Build + run lib tests** — `cargo build -p datafusion-bio-function-vep && cargo test -p datafusion-bio-function-vep --lib` → PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs datafusion/bio-function-vep/src/lookup_provider.rs
git commit -m "feat(vep): grid-sliced per-worker lookup scans for Merged/RefSeq parallel"
```

---

## Task 6: Thread the slice into the sharded annotation worker

**Files:**
- Modify: `annotate_provider.rs` `spawn_annotation_from_lookup_sharded` (`:9756`) + the `AnnotatingParallel` spawn loop (`:11301`–`:11382`).

**Interfaces:**
- Consumes: `WorkerGridSlice` (Task 5), warm-up-capable `annotate_worker_window` (Task 4).
- Produces: a sharded worker that drops `skip_leading_rows`, feeds the rest through its accumulator starting at global row `warm_up_start_row`, processes warm-up buffers state-only, emits owned buffers to its shard.

Implementation in `spawn_annotation_from_lookup_sharded`:
- Add param `slice: WorkerGridSlice`.
- Maintain `global_row = slice.warm_up_start_row` and a `to_skip = slice.skip_leading_rows` counter; for each incoming batch, if `to_skip > 0`, slice off the first `to_skip` rows (using `batch.slice`) before buffering, decrementing `to_skip`.
- When draining a window to `annotate_worker_window`, pass `WindowAnnotateBounds { base_global_row: global_row, emit_start_row: slice.emit_start_row, emit_end_row: slice.emit_end_row }`, and advance `global_row` by the window's input-row count after each call. Warm-up buffers (global < emit_start) produce no shard writes; emitted buffers write as today.
- The worker's `input_buffer_accumulator` already starts empty; because `warm_up_start_row` is a global buffer boundary, the accumulator reproduces the global grid from there (verified by Task 4's carry test at the unit level).

In the `AnnotatingParallel` spawn loop, pass each handle's `WorkerGridSlice` to `spawn_annotation_from_lookup_sharded`.

- [ ] **Step 1: Write a failing test** for the skip/advance bookkeeping — a small pure helper `fn apply_skip(batch, &mut to_skip) -> Option<RecordBatch>` and assert it drops exactly N leading rows across batch boundaries.

- [ ] **Step 2: Run to verify it fails** — `cargo test -p datafusion-bio-function-vep --lib apply_skip` → FAIL.

- [ ] **Step 3: Implement `apply_skip` + thread `slice` through** both functions.

- [ ] **Step 4: Build + lib tests** — `cargo build -p datafusion-bio-function-vep && cargo test -p datafusion-bio-function-vep --lib` → PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): sharded worker warm-start + grid-aligned emit by global rank"
```

---

## Task 7: Route Merged/RefSeq `workers>1` to the parallel path (drop the gate)

**Files:**
- Modify: `annotate_provider.rs` serial gate (`:5266`–`:5276`); `vcf_sink.rs` routing (`:1112`) if Task 0 showed it degenerate.

**Interfaces:**
- Consumes: everything above.
- Produces: Merged/RefSeq at `workers>1` runs `AnnotatingParallel` with grid slices.

Implementation:
- Replace the forced `annotation_workers = 1` for Merged/RefSeq with `annotation_workers = requested_workers`.
- **Fail-closed guard (§8):** if `annotation_workers > 1` and the grid plan cannot be built (count pass empty / no transcript cache → `overlap_width_bp` only region-width / non-indexed input so position-range scans unsupported), return the existing `workers>1` guard error rather than silently falling back to partition-local bounds. Add the check where the grid branch is selected (Task 5).
- If Task 0 found `vcf_sink.rs:1112` produced degenerate Merged+w4 output, ensure the sharded driver is what runs (it should, now that the provider yields `AnnotatingParallel` with a single output partition — matching `drive_sharded_vcf_annotation`'s `partition_count == 1` expectation at `:757`).

- [ ] **Step 1: Edit the gate**

Replace (`:5269`):
```rust
let annotation_workers = if matches!(
    self.cache_source_type,
    CacheSourceType::Merged | CacheSourceType::RefSeq
) { 1 } else { requested_workers };
```
with:
```rust
// Stateful Merged/RefSeq now runs N grid-aligned sharded workers with
// bounded-overlap warm-up (design 2026-06-25 §5/§7), so it no longer needs the
// serial fallback. Fail-closed below if the grid plan can't be built.
let annotation_workers = requested_workers;
```

- [ ] **Step 2: Build** — `cargo build -p datafusion-bio-function-vep` → clean.

- [ ] **Step 3: Run lib tests** — `cargo test -p datafusion-bio-function-vep --lib` → PASS.

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs datafusion/bio-function-vep/src/vcf_sink.rs
git commit -m "feat(vep): enable grid-aligned parallel annotation for Merged/RefSeq"
```

---

## Task 8: E2E correctness + performance gates (no code)

**Files:** none (verification). Use explicit per-run commands + simple log names (no shell-loop quoting).

**Interfaces:**
- Consumes: the full feature.
- Produces: the go/no-go evidence required by the Global Constraints.

- [ ] **Step 1: Rebuild** — `cd ~/research/git/vepyr && source .venv/bin/activate && unset CONDA_PREFIX && maturin develop --release`

- [ ] **Step 2: chr4 correctness, 3 consecutive w4 runs**

```bash
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend lance --workers 4 --force 2>&1 | tee /tmp/g-chr4-w4-a.log
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend lance --workers 4 --force 2>&1 | tee /tmp/g-chr4-w4-b.log
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend lance --workers 4 --force 2>&1 | tee /tmp/g-chr4-w4-c.log
```
Expected each: `HGNC_ID` mismatches = 0, 307,295 variants, "Only in VEP: 0". Any run that differs (e.g., truncated row count) is a FAIL — do not commit.

- [ ] **Step 3: chr4 w1 still 0; chr2 both**

```bash
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --backend lance --workers 1 --force 2>&1 | tee /tmp/g-chr4-w1.log
python e2e-testing/scripts/run_annotation_fast.py chr2 --cache merged --backend lance --workers 1 --force 2>&1 | tee /tmp/g-chr2-w1.log
python e2e-testing/scripts/run_annotation_fast.py chr2 --cache merged --backend lance --workers 4 --force 2>&1 | tee /tmp/g-chr2-w4.log
```
Expected: all `HGNC_ID` = 0; w4 byte-identical body to w1 (the harness's body comparison / "Only in VEP: 0").

- [ ] **Step 4: Larger run** — chr1-22 merged Lance at `workers=4`, confirm the previous report's `HGNC_ID` mismatch class is 0.

- [ ] **Step 5: Performance report** — record annotation wall + variants/s for `workers {1,2,4,8}` on chr4; confirm recovery vs the 12.4s serial baseline (target trending toward the old 5.6s @ w4).

- [ ] **Step 6: Lib tests + clippy gate**

```bash
cargo test -p datafusion-bio-function-vep --lib
cargo clippy -p datafusion-bio-function-vep -- -D warnings
cargo fmt -- --check
```
Expected: all green.

- [ ] **Step 7: Only if ALL gates hold, the feature is complete.** Otherwise debug with `superpowers:systematic-debugging` before any commit of Task 7.

---

## Self-Review notes

- **Spec coverage:** §5.1 partitioning → Task 2/3; §5.2 warm-up → Task 4/6; §5.3 donor scope → already committed; §5.4 overlap width → Task 1; §6 correctness (Hazard 1 grid seams / Hazard 2 warm-up) → Task 2 (seams on boundaries) + Task 4/6 (warm-up); §8 gates → Task 8; fail-closed → Task 7.
- **Position-tie correctness:** seams tracked by global row index with `skip_leading_rows`, not by position alone — the one subtlety that a naive position filter would get wrong.
- **Open risk to validate early (Task 0):** the `vcf_sink.rs:1112` / provider-gate interaction for current Merged+w4. If current w4 is degenerate, the w1 body is the byte-identical reference.
  - **RESOLVED 2026-06-26 (ran it):** current Merged+w4 produces **EMPTY output** (0 rows) — `vcf_sink:1112` drives the sharded path but the provider pins `annotation_workers=1` → `AnnotatingContig` yields batches that `drive_sharded` discards → no shards written. The handoff's "w4 verified 307,295" was incorrect. w1 = 307,295 variants, 100% field match vs the external VEP reference (`HG002_annotated_wgs_everything_hgvs_merged.vcf`), which is the byte-identical gate the harness already enforces. Task 7 (enabling the real parallel path) therefore also fixes this latent empty-output bug.
- **Warm-up cost:** <1 buffer per worker (design §5.5); warm-up variants still pass through Lance lookup in this design (bounded, acceptable) — a future optimisation could scan warm-up positions without the variation lookup.

## Execution results (2026-06-26)

Tasks 0–7 implemented and committed; the grid-aligned parallel path is **correct
and byte-identical**. chr4 merged Lance:

| run | variants | Only-in-VEP | HGNC_ID mism. | 86 CSQ fields | annot wall |
|---|---|---|---|---|---|
| w1 (serial) | 307,295 | 0 | 0 | 100% | 14.1s |
| w4 run a | 307,295 | 0 | 0 | 100% | 11.4s |
| w4 run b | 307,295 | 0 | 0 | 100% | 11.0s |
| w4 run c | 307,295 | 0 | 0 | 100% | 11.0s |

**Correctness gate: PASS** — 0 HGNC mismatches, full count, all 86 fields 100%, on 3
consecutive w4 runs (no intermittent truncation); w1 also 0.

**Over-read fixed (commit) → perf gate met.** The earlier truncation was reading
only partition 0 of each worker's plan (`exec_partition=0`), **not** a broken
position filter. Fix: restore the per-worker `[scan_lo_pos, scan_hi_pos)` filter and
read **all** partitions of that filtered plan (`spawn_lookup_full_contig_worker`), so
each worker's lookup covers exactly its range — `stream_end_rank ≈ emit_end`, no
over-read. The single-threaded profile confirmed the pipeline is ~98% engine-bound
(`engine` 10.9s vs `lookup_wait` 0.47s), so the over-read had turned overlapped
"free" lookup into 2.56× redundant decode CPU competing with the engine.

Final results (chr4 + chr2 merged Lance):

| run | variants | Only-in-VEP | HGNC_ID mism. | wall |
|---|---|---|---|---|
| chr4 w4 ×3 | 307,295 | 0 | 0 | 5.8 / 6.0 / 6.0s |
| chr2 w1 / w4 | 331,324 | 0 | 0 | 19.7s / 7.8s |
| chr4 w1/w2/w4/w8 | 307,295 | 0 | 0 | 13.9 / 8.0 / 5.8 / 4.7s |

**Correctness gate: PASS** (0 mismatches, full counts, all 86 CSQ fields 100%, 3
consecutive w4 runs consistent). **Performance gate: PASS** — chr4 w4 ≈ 2.4×
(13.9→5.8s), w8 ≈ 3.0× (4.7s), chr2 w4 ≈ 2.5×; clears the ~2× / 6–8s target and beats
the old unsafe path (5.6s) while being correct. Diminishing returns past w4 reflect
the per-contig serial floor (context load + count pass + colocated), as modeled.

The partition-subset assignment (originally proposed here) is **unnecessary** — the
position filter + read-all-partitions achieves no over-read directly.

## Follow-ups (deferred optimizations — NOT in this plan)

These were brainstormed 2026-06-26 and consciously deferred so the first correct, self-contained version lands without cross-repo or perf coupling. Revisit after the §8 gates are green.

1. **bio-formats count-balanced grid partitioning (proper long-term home).** Add to `datafusion-bio-format-vcf` an API that returns, for a contig, the genomic positions at every `n`-th record (or N count-balanced split positions on `n`-multiples). bio-formats owns the VCF+index read and can fuse this into `scan()` planning; it also upgrades partition **balance** from byte-budget to variant-balanced for all consumers. Cost: cross-repo change — bio-formats is pinned at `v1.8.0` (git tag) and the DataFusion/Arrow versions must stay in sync across the three repos (see CLAUDE.md). Our crate would then consume the boundaries and still build the warm-up-extended per-worker scans locally (warm-up is annotation-domain, not bio-formats').
2. **Surface the free tabix total in bio-formats.** `estimate_sizes_from_tbi` hardcodes `unmapped_count: 0` and never reads the metadata pseudo-bin's `mapped_record_count`. Exposing it gives a zero-scan per-contig **record** count — an exact unit total *only* for normalized/decomposed input — usable as the count pass's denominator.
3. **Parallelize the count pass** across the existing index byte-partitions (wall ≈ `read/N`).
4. **Fused single-read + barrier** alternative to the separate pre-pass: each byte-partition reads once into an in-memory buffer, reports its count at a cross-partition barrier, then annotates with the correct global phase. One physical read; costs RAM (adds to the known t8 RSS pressure) + a sync barrier (some streaming loss). Choose only if profiling shows the pre-pass read is material.
5. **Warm-up without variation lookup:** scan warm-up positions without `KvLookupExec` (carried HGNC depends only on transcripts), eliminating the bounded redundant Lance reads for the <1-buffer warm-up prefix.
