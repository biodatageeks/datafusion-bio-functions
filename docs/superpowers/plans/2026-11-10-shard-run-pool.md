# Shard Run Pool Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** VCF-file output at workers > 1 schedules many small grid-ordered runs dynamically instead of one static slice per worker, and hands each finished run's shard to the assembler immediately.

**Architecture:** The streaming path's `RunPoolState` already cuts a contig into runs and starts them as worker slots free up. This plan lets a run write a VCF body shard (`RunOutput::Shard`) instead of sending batches, numbers shards in run order, and publishes a one-shard `ContigShardRange` when the head run is released. Runs are cut by the existing cost weights where they exist, so a heavy buffer becomes its own run. The static grid path stays reachable with `VEP_SHARD_RUN_POOL=0`.

**Tech Stack:** Rust 2024, tokio, DataFusion 53; all changes in `datafusion/bio-function-vep/src/annotate_provider.rs` and `src/vcf_sink.rs`.

**Spec:** `docs/superpowers/specs/2026-11-10-parallel-vcf-output-bottlenecks-design.md` (sections 1.1, 1.2, 3.1, 4)

## Global Constraints

- Branch `perf/shard-run-pool` off `origin/master`; plan 0 (`2026-11-10-local-gate.md`) is merged first.
- Output bytes are identical to the static path for every input, cache type and worker count.
- Runs are whole VEP input buffers; buffer boundaries, warm-up replay (`overlap_width_bp`) and output order do not change.
- The streaming `RecordBatch` path keeps its run length, look-ahead and VCF scan planning.
- Gate targets (spec section 4): 8-worker sum of `time_s` at least 15% lower on ClinVar and 5% lower on HG002.
- Line numbers below refer to `68d0f7b`.

## Review Focus

- A contig with fewer buffers than workers (ClinVar Y and MT have one): exactly one run, one shard, correct output (gate contigs Y and MT; unit test in Task 2).
- A run that emits no rows still creates its shard file, because the assembler opens every published id (unit test in Task 4).
- Several contigs in one invocation: shard ids stay globally ascending across contigs (Task 6, step 3).
- An Ensembl cache, where vepyr#282 came from a cut off the buffer grid: identical output with the pool on and off (Task 6, step 4).
- A run failing after earlier shards were published: the assembler must discard `<output>.part` (existing test `assembler_discards_partial_output_when_not_finalized`; the pool's error path goes through `fail_run_pool`, which drops the sender without setting `finalize`; confirm in review).

---

### Task 1: Mode switch

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`contig_modes` at 16079, its call at 16225, `ContigAnnotationConfig` at 10243 and its construction at 6234, `activate_run_lookup` at 16782, tests at 21598)

**Interfaces:**
- Produces:
  - `fn shard_run_pool_enabled(env: Option<&str>) -> bool` — `true` when unset; otherwise `pipeline_trace::enabled_from_env_value`. Variable: `VEP_SHARD_RUN_POOL`.
  - `ContigAnnotationConfig::shard_run_pool: bool`, set from the environment where the config is built.
  - `fn contig_modes(cache_enabled: bool, annotation_workers: usize, has_shard_ctx: bool, cache_source_type: CacheSourceType, has_regions: bool, shard_run_pool: bool) -> ContigModes`
  - `ContigAnnotationConfig::shard_pool(&self) -> bool` — `annotation_workers > 1 && vcf_shard_ctx.is_some() && shard_run_pool`.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn shard_run_pool_defaults_on_and_zero_disables() {
    assert!(shard_run_pool_enabled(None));
    assert!(shard_run_pool_enabled(Some("1")));
    assert!(!shard_run_pool_enabled(Some("0")));
    assert!(!shard_run_pool_enabled(Some("false")));
}

#[test]
fn sink_uses_the_run_pool_unless_disabled() {
    for source in [CacheSourceType::Merged, CacheSourceType::RefSeq, CacheSourceType::Ensembl] {
        let on = contig_modes(true, 4, true, source, false, true);
        assert!(on.stream_parallel && !on.stateful_parallel && !on.stateful_runs);
        let off = contig_modes(true, 4, true, source, false, false);
        assert!(!off.stream_parallel && off.stateful_parallel);
    }
    // One worker and the streaming path do not depend on the switch.
    assert_eq!(contig_modes(true, 1, true, CacheSourceType::Merged, false, true),
               contig_modes(true, 1, true, CacheSourceType::Merged, false, false));
    assert_eq!(contig_modes(true, 4, false, CacheSourceType::Merged, false, true),
               contig_modes(true, 4, false, CacheSourceType::Merged, false, false));
}
```
Update `contig_modes_are_exclusive` to pass `false` as the new last argument.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib shard_run_pool sink_uses_the_run_pool`
Expected: compile error, `shard_run_pool_enabled` not found.

- [ ] **Step 3: Implement**

`stream_parallel = parallel && (!has_shard_ctx || shard_run_pool)`; `stateful_parallel = parallel && has_shard_ctx && !shard_run_pool`; `stateful_runs` keeps its formula. In `activate_run_lookup` replace the `pool_run` condition by `cache_enabled && config.annotation_workers > 1 && (config.vcf_shard_ctx.is_none() || config.shard_run_pool)`, and in that arm call `provider.set_single_vcf_scan(true)` only when `config.shard_pool()`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep --lib contig_modes shard_run_pool sink_uses_the_run_pool`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): mode switch for the shard run pool"
```

### Task 2: Weighted run cuts

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (next to `plan_stream_runs` at 11808; tests next to `plan_stream_runs_cuts_open_range_into_fixed_pieces` at 21507)

**Interfaces:**
- Consumes: `build_grid_slices`, `plan_stream_runs`, `GridBufferBoundary`, `RunBounds`, `ContigRun`.
- Produces:
  - `const SHARD_RUNS_PER_WORKER: usize = 4;`
  - `fn shard_runs_per_worker(env: Option<&str>) -> usize` — a positive integer from `VEP_SHARD_RUNS_PER_WORKER`, else the constant.
  - `fn weighted_run_cuts(weights: &[u64], bk: usize, bk1: usize, cap: u64) -> Vec<usize>` — ascending buffer indices from `bk` to `bk1` inclusive.
  - `fn plan_shard_runs(boundaries: &[GridBufferBoundary], ranges: &[(Vec<RunBounds>, (usize, usize))], weights: Option<&[u64]>, workers: usize, runs_per_worker: usize, overlap_width_bp: i64) -> Vec<ContigRun>`

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn weighted_run_cuts_isolate_a_heavy_buffer() {
    assert_eq!(weighted_run_cuts(&[1, 1, 10, 1, 1, 1], 0, 6, 3), vec![0, 2, 3, 6]);
}

#[test]
fn weighted_run_cuts_fill_even_weights_to_the_cap() {
    assert_eq!(weighted_run_cuts(&[1; 8], 0, 8, 2), vec![0, 2, 4, 6, 8]);
}

#[test]
fn weighted_run_cuts_stay_inside_the_range() {
    assert_eq!(weighted_run_cuts(&[5, 5, 1, 1, 1, 5], 2, 5, 2), vec![2, 4, 5]);
}

#[test]
fn weighted_run_cuts_one_buffer_is_one_run() {
    assert_eq!(weighted_run_cuts(&[7], 0, 1, 1), vec![0, 1]);
}

#[test]
fn plan_shard_runs_without_weights_cuts_fixed_pieces() {
    // 10 buffers, 2 workers x 4 runs: ceil(10 / 8) = 2 buffers per run.
    let boundaries = grid(10);            // helper used by the plan_stream_runs tests
    let runs = plan_shard_runs(&boundaries, &[(vec![RunBounds::OPEN], (0, 10))], None, 2, 4, 0);
    assert_eq!(runs.len(), 5);
}

#[test]
fn plan_shard_runs_cover_every_row_once_in_order() {
    let boundaries = grid(6);
    let weights = [1, 1, 10, 1, 1, 1];
    let runs = plan_shard_runs(&boundaries, &[(vec![RunBounds::OPEN], (0, 6))], Some(&weights), 2, 2, 0);
    let slices: Vec<_> = runs.iter().map(|r| r.slice.clone().unwrap()).collect();
    assert_eq!(slices.first().unwrap().emit_start_row, boundaries[0].global_row);
    assert_eq!(slices.last().unwrap().emit_end_row, boundaries[6].global_row);
    for pair in slices.windows(2) {
        assert_eq!(pair[0].emit_end_row, pair[1].emit_start_row);
    }
    // total 15, cap ceil(15 / 4) = 4: the weight-10 buffer is alone.
    assert!(slices.iter().any(|s| s.emit_start_row == boundaries[2].global_row
        && s.emit_end_row == boundaries[3].global_row));
}
```
If the existing tests build their grid inline, extract that into `fn grid(buffers: usize) -> Vec<GridBufferBoundary>` (5,000 rows per buffer, one position per row) and use it in both.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib weighted_run_cuts plan_shard_runs`
Expected: compile error, `weighted_run_cuts` not found.

- [ ] **Step 3: Implement**

`weighted_run_cuts`: walk `bk..bk1` accumulating weight; before adding buffer `i`, if the current run is non-empty and `acc + weights[i] > cap`, push `i` as a cut and reset `acc`. `plan_shard_runs` with weights: `cap = max(1, ceil(sum of weights over the covered buffers / (workers * runs_per_worker)))`, cuts per range from `weighted_run_cuts`, slices from `build_grid_slices`, one `ContigRun` per slice with the range's bounds. Without weights, or when `weights.len() != boundaries.len() - 1`: `plan_stream_runs(boundaries, ranges, max(1, ceil(covered / (workers * runs_per_worker))), overlap_width_bp)`. No four-buffer floor on either branch.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep --lib weighted_run_cuts plan_shard_runs plan_stream_runs`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): weight-capped run cuts for sharded output"
```

### Task 3: Plan shard-pool runs in `prepare_contig_data`

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`prepare_contig_data`: `grid_balance_enabled` at 16302, the count pass call at 16299, the weight computation at 16443-16476, the `stream_parallel` run planning at 16519-16562)

**Interfaces:**
- Consumes: `plan_shard_runs`, `shard_runs_per_worker` (Task 2); `ContigAnnotationConfig::shard_pool` (Task 1).
- Produces: `fn grid_buffer_weights(boundaries: &[GridBufferBoundary], positions: &[(i64, u32)], translations: &[TranslationFeature], transcripts: &[TranscriptFeature]) -> Vec<u64>` — the weight computation now inline in the `stateful_parallel` branch, moved unchanged. Trace event `run_pool plan` gains fields `weighted` (0 or 1) and `shard` (0 or 1).

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn grid_buffer_weights_match_the_inline_computation() {
    // Two buffers; one transcript with a 7,680 nt CDS covers only buffer 1.
    // coding_length_cost(7680) = 1 + 7680 / 768 = 11; the floor is total / B / 8.
    let boundaries = grid(2);
    let positions: Vec<(i64, u32)> = (0..10_000).map(|i| (i as i64 + 1, 1)).collect();
    let (tx, tl) = one_coding_transcript(5_001, 10_000, 7_680);   // start, end, cds_len
    let weights = grid_buffer_weights(&boundaries, &positions, &[tl], &[tx]);
    let total = 5_000u64 * 11;
    let fixed = total / 2 / 8;
    assert_eq!(weights, vec![fixed + 1, total + fixed + 1]);
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib grid_buffer_weights_match`
Expected: compile error, `grid_buffer_weights` not found.

- [ ] **Step 3: Implement**

Move the `coding_costs` map and the `buffer_transcript_cost_weights` call into `grid_buffer_weights` and call it from the `stateful_parallel` branch. Set `grid_balance_enabled = (stateful_parallel || config.shard_pool()) && use_weighted_grid_split(..)`; the count pass collects positions under the same condition as before, now through this flag. In the `stream_parallel` branch, when `config.shard_pool()`: compute `weights` if `grid_balance_enabled` and call `plan_shard_runs(&boundaries, &ranges, weights.as_deref(), config.annotation_workers, shard_runs_per_worker(env), overlap)`; otherwise keep `plan_stream_runs`.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib grid_ plan_ buffer_transcript`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): plan weight-capped runs for sharded output"
```

### Task 4: Runs write shards and publish on release

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`RunStats` at 13070, `RunPoolState` at 12343, `spawn_run_task` at 12416, pool construction at 15272-15316, the release arm at 15905-15925)
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (assembler test at 2690)

**Interfaces:**
- Consumes: `VcfBodyShardWriter::create`, `RunOutput::Shard`, `ContigShardRange`, `VcfShardContext` (existing); `ContigAnnotationConfig::shard_pool` (Task 1).
- Produces:
  - `struct ShardPoolOutput { ctx: Arc<crate::vcf_sink::VcfShardContext>, first_shard_id: usize }`
  - `RunPoolState::shard: Option<ShardPoolOutput>`
  - `RunStats::output_lines: usize` — shard body lines; 0 for batch output.
  - `fn run_shard_path(tempdir: &Path, first_shard_id: usize, index: usize) -> PathBuf` — `partition_{:04}.vcf.body` for `first_shard_id + index`.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn run_shard_path_numbers_shards_in_run_order() {
    let dir = std::path::Path::new("/t");
    assert_eq!(run_shard_path(dir, 0, 0), dir.join("partition_0000.vcf.body"));
    assert_eq!(run_shard_path(dir, 24, 7), dir.join("partition_0031.vcf.body"));
    assert_eq!(run_shard_path(dir, 9_999, 1), dir.join("partition_10000.vcf.body"));
}

#[test]
fn shard_pool_admission_is_not_bounded_by_the_release_window() {
    // Head run 0 still running, 8 workers, 40 runs, lookahead = total runs.
    assert!(admission_allows(7, 0, 30, 40, 8, 40));
    assert!(!admission_allows(8, 0, 30, 40, 8, 40)); // no free worker
}
```
In `vcf_sink.rs`, extend `assembler_plain_appends_ranges_in_order_and_counts_lines` with a zero-byte shard published as its own one-shard range between two non-empty ones; the assembled body must equal the concatenation of the non-empty shards.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib run_shard_path shard_pool_admission assembler_plain`
Expected: compile error, `run_shard_path` not found.

- [ ] **Step 3: Implement**

At pool construction: build the pool when `config.vcf_shard_ctx.is_none() || config.shard_pool()`. For a shard pool set `shard = Some(ShardPoolOutput { ctx, first_shard_id: self.next_global_shard_id })`, advance `self.next_global_shard_id` by `total_runs`, and use `lookahead = total_runs`. In `spawn_run_task`, with `pool.shard`: create the `VcfBodyShardWriter` at `run_shard_path(..)` before activating the lookup, annotate into `RunOutput::Shard { shard, rows_done: Arc::clone(&ctx.rows_done) }`, then take the writer back, copy `lines` into `stats.output_lines`, and call `finish()`. Keep the task's `tx` alive until the task returns, so the head's receiver closes exactly when its run ends. In the release arm, for a shard pool: add `stats.input_rows` to `pool.contig_rows` and send `ContigShardRange { first_shard_id: first + head, end_shard_id: first + head + 1, output_lines: stats.output_lines }`, ignoring a closed channel.

- [ ] **Step 4: Run the whole crate's tests**

Run: `cargo test -p datafusion-bio-function-vep`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): run pool writes VCF shards and publishes them in order"
```

### Task 5: Assembler trace events

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (`run_assembler_thread` at 1478)

**Interfaces:**
- Produces: trace events `vcf_assemble range` (`first`, `end`, `bytes`, `elapsed_ms`) after each range and `vcf_assemble finalize` (`bytes`, `elapsed_ms`: time from the first range received to the rename). Plan 2's entry criterion reads them.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn assembler_range_trace_line_has_stable_fields() {
    let line = assembler_range_trace(3, 4, 1_048_576, std::time::Duration::from_millis(12));
    assert_eq!(line, vec![
        ("first", TraceValue::Usize(3)), ("end", TraceValue::Usize(4)),
        ("bytes", TraceValue::Usize(1_048_576)),
        ("elapsed", TraceValue::Duration(std::time::Duration::from_millis(12))),
    ]);
}
```
(`PipelineTraceValue` needs `PartialEq` derived for the comparison.)

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib assembler_range_trace`
Expected: compile error.

- [ ] **Step 3: Implement** `fn assembler_range_trace(first: usize, end: usize, bytes: usize, elapsed: Duration) -> Vec<(&'static str, TraceValue<'static>)>` and emit both events from `run_assembler_thread`.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib assembler`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): trace assembler ranges"
```

### Task 6: Smoke, focus checks, tuning and the gate

**Files:**
- Modify: `datafusion/bio-function-vep/PERFORMANCE.md` (new section "Shard run pool")
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`SHARD_RUNS_PER_WORKER`, only if step 5 selects another value)

**Interfaces:**
- Consumes: `build_venv.sh`, `run_arm.sh`, `report.py`, `gate.sh` (plan 0).

- [ ] **Step 1: Lint and unit tests**

Run: `cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all`
Expected: all pass.

- [ ] **Step 2: Smoke on chr22 and ClinVar Y, MT**

```bash
G="$HOME/research/benchmarks/local-gate"; L=datafusion/bio-function-vep/scripts/local_gate
$L/build_venv.sh "$PWD" "$G/p1/venv-candidate"
$L/run_arm.sh "$G/p1/venv-candidate" "$G/p1/smoke/pool" hg002 8 parity --chroms chr22
$L/run_arm.sh "$G/p1/venv-candidate" "$G/p1/smoke/static" hg002 8 time --chroms chr22 VEP_SHARD_RUN_POOL=0
$L/run_arm.sh "$G/p1/venv-candidate" "$G/p1/smoke/pool" clinvar 8 parity --chroms "Y MT 21"
$L/run_arm.sh "$G/p1/venv-candidate" "$G/p1/smoke/static" clinvar 8 time --chroms "Y MT 21" VEP_SHARD_RUN_POOL=0
```
Expected: `verdict` `PASS` on every `pool` row and the same `body_md5` in `pool` and `static` for each contig.

- [ ] **Step 3: Two contigs in one invocation**

```bash
D="$HOME/workspace/data_vepyr"; T="$G/p1/multi"; mkdir -p "$T"
tabix -h "$D/input/HG002_normalized.vcf.gz" chr21 chr22 | bgzip > "$T/in.vcf.gz" && tabix -p vcf "$T/in.vcf.gz"
for v in 1 0; do VEP_SHARD_RUN_POOL=$v "$G/p1/venv-candidate/bin/python" -P -c 'import sys,vepyr
vepyr.annotate(sys.argv[1], sys.argv[2], reference_fasta=sys.argv[3], output_vcf=sys.argv[4], workers=8, compression="plain", preserve_record_layout=True, show_progress=False)' \
  "$T/in.vcf.gz" "$D/cache/116_GRCh38_merged" "$D/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa" "$T/out.$v.vcf"; grep -v '^#' "$T/out.$v.vcf" | md5; done
```
Expected: two identical digests; `grep -vc '^#'` gives 106673 for both.

- [ ] **Step 4: Ensembl identity**

Repeat step 3's loop with `"$D/cache/116_GRCh38_ensembl"` on two inputs: `tabix -h "$D/input/clinvar/clinvar_20261011.norm.vcf.gz" 9` and `tabix -h "$D/input/HG002_normalized.vcf.gz" chr1`.
Expected: identical digests per input with the pool on and off.

- [ ] **Step 5: Choose runs per worker**

For `VEP_SHARD_RUNS_PER_WORKER` in 2, 4, 8: `run_arm.sh … clinvar 8 time --chroms "2 17 6" VEP_SHARD_RUNS_PER_WORKER=<v>` and `run_arm.sh … hg002 8 time --chroms chr1 VEP_SHARD_RUNS_PER_WORKER=<v>`, two arms each, order 2, 4, 8, 8, 4, 2. Sum the mean `time_s` over the four contigs. Keep 4 unless another value is more than 2% lower; then set `SHARD_RUNS_PER_WORKER` to it and commit. Record the table in `PERFORMANCE.md`.

- [ ] **Step 6: Full gate**

```bash
$L/gate.sh "$G/functions-baseline" "$PWD" "$G/p1/gate" --target-pct HG002=5,ClinVar=15 --bgzf-identity
cat "$G/p1/gate/verdict.txt"
```
Expected: every rule `PASS`. If `target` fails while parity and identity pass, stop and report the per-contig table; do not relax the target.

- [ ] **Step 7: Residual copy tail for plan 2**

Run ClinVar contigs 1, 2 and 17 once each with `VEP_PIPELINE_TRACE=1` through `run_arm.sh … clinvar 8 time`. For each log compute `tail = t(vcf_assemble finalize) - t(last run_pool done)` and `wall = t(vcf_assemble finalize)`. Record `tail`, `wall` and `tail / wall` in `PERFORMANCE.md`.

- [ ] **Step 8: Document and open the PR**

Write the "Shard run pool" section: what changed, `VEP_SHARD_RUN_POOL` and `VEP_SHARD_RUNS_PER_WORKER`, build identities from both `GATE_IDENTITY.json` files, the gate's per-dataset totals and verdict, the step 5 and step 7 tables.

```bash
git commit -am "docs(vep): shard run pool measurements"
```
After the PR merges, a vepyr pin PR and its `/parity-wgs` run follow the usual release flow.
