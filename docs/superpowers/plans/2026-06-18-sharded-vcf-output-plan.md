# Sharded VCF Output (per-worker streaming shards) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `threads>1` VEP path fully parallel from lookup to VCF write — each fused worker streams its own VCF body shard (no ordered drain, no output channel) — then a final merge concatenates shards in position order with one header.

**Architecture:** The worker's loop is extended: `lookup → hydrate → annotate → format_vcf_body_chunk → write_all` into its own `BufWriter<File>` shard. `ParallelContigState`'s ordered drain and the `spawn_annotation_from_lookup` output channel are removed. `vcf_sink` derives the VCF formatter context (info fields / format tags / sample names) and the tempdir, drives the shard-writing annotation, then `write_header` once + concat shards in ascending (position-ordered) index — reusing the existing `format_vcf_body_chunk` + `copy_body_file_to_writer` code.

**Tech Stack:** Rust 2024, DataFusion 52.1, Arrow 57, tokio, `lance`. Crate `datafusion-bio-function-vep`, feature `lance-cache`. Spec: `docs/superpowers/specs/2026-06-18-sharded-vcf-output-design.md`.

**TDD adaptation:** Output is byte-identical, so the primary gates are the e2e parity diff + `send_wait≈0` + the RSS/timing benchmark, plus "lib failure set unchanged vs baseline" (35 pre-existing). New pure helpers (the reusable body-serializer) get a focused unit test.

**Verify each task** with `cargo build -p datafusion-bio-function-vep --features lance-cache` and (for behavior) the parity gate in Task 8. Commit after each task.

---

## File Structure

| File | Responsibility / change |
|---|---|
| `src/vcf_sink.rs` | Factor the reusable per-batch body serializer (`format_vcf_body_chunk` already exists; extract a `VcfBodyShardWriter` that owns a `BufWriter<File>` + formatter context and exposes `write_batch(&RecordBatch)` + `finish()`). Add the `threads>1` orchestration entry that derives the formatter context, drives shard-writing annotation, and does header + concat (reusing `write_header` + the `next_write_partition_id` concat loop). |
| `src/annotate_provider.rs` | `spawn_annotation_from_lookup`: write shard via `VcfBodyShardWriter` instead of `tx.send`. `ParallelContigState` / `AnnotatingParallel` drain: removed; collect shard descriptors. Thread the shard output context (formatter ctx + tempdir) into the worker spawn. |

The reusable serializer lives in `vcf_sink.rs` (its formatting helpers do) and is called from the engine worker — `annotate_provider.rs` depends on `vcf_sink` (or a shared `vcf_format` submodule) for it.

---

## Task 0: De-risk the engine↔sink handoff (design spike, no behavior change)

**The crux:** `vcf_sink::annotate_to_vcf` builds the annotation via SQL over the `annotate_vep` UDTF (`vcf_sink.rs:1205-1215`), so the post-schema formatter context (sample names / info fields / format tags, derived at `:1214-1255`) cannot be injected through SQL into the workers. Resolve this first.

**Decision (implement this):** For the `threads>1` sharded path, `vcf_sink` bypasses the SQL/DataFrame execution and drives the annotation through a **direct entry** that accepts the formatter context + tempdir and returns shard descriptors. Concretely add:

- `struct VcfShardContext { vcf_info_fields: Arc<Vec<String>>, unique_format_tags: Arc<Vec<String>>, sample_names: Arc<Vec<String>>, coordinate_zero_based: bool, tempdir: PathBuf }` in `vcf_sink.rs` (or a shared module), `Clone`.
- A field `vcf_shard_ctx: Option<Arc<VcfShardContext>>` on `ContigAnnotationConfig` (default `None`; set only by the sharded entry). Plumb it from the `annotate_vep` UDTF path is NOT possible via SQL — therefore the sharded path constructs the `ContigAnnotationExec` plan directly (the UDTF's plan-construction function is reused) with `vcf_shard_ctx = Some(...)`, instead of `ctx.sql(...)`.

- [ ] **Step 1: Locate the plan-construction reuse point.**

Run: `grep -n "scan_with_transcript_engine_partitioned\|ContigAnnotationExec::new\|fn annotate_vep\|register_vep_functions\|TableFunctionImpl\|fn call" src/*.rs | head` (in `datafusion/bio-function-vep`).
Expected: find the function that builds `ContigAnnotationExec` (the UDTF `call`/scan). Confirm it can be invoked directly by `vcf_sink` with a `ContigAnnotationConfig` (it is `scan_with_transcript_engine_partitioned`, `annotate_provider.rs:5449`).

- [ ] **Step 2: Add the `VcfShardContext` type + config field (no use yet).**

In `vcf_sink.rs` add the `VcfShardContext` struct above `annotate_to_vcf`. In `annotate_provider.rs` add `vcf_shard_ctx: Option<Arc<crate::vcf_sink::VcfShardContext>>` to `ContigAnnotationConfig` (find the struct; add `#[derive]`-compatible field, default `None` in every constructor — the compiler will name them).

- [ ] **Step 3: Build.**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error|Finished" | head`
Expected: `Finished` (field added, unused).

- [ ] **Step 4: Commit.**

```bash
git add -A && git commit -m "feat(vep): add VcfShardContext + ContigAnnotationConfig.vcf_shard_ctx (unused)"
```

---

## Task 1: Reusable `VcfBodyShardWriter` (extracted from `write_vcf_partition_body`)

**Files:** `src/vcf_sink.rs` (refactor `write_vcf_partition_body` `:320` to use the new writer; add the writer struct).

- [ ] **Step 1: Add the struct + methods.**

```rust
// vcf_sink.rs — reusable streaming VCF body shard writer (used by the engine
// workers AND by write_vcf_partition_body).
pub(crate) struct VcfBodyShardWriter {
    writer: std::io::BufWriter<std::fs::File>,
    ctx: std::sync::Arc<VcfShardContext>,
    batch_id: usize,
    pub(crate) lines: usize,
    pub(crate) input_rows: usize,
    pub(crate) bytes: usize,
}

impl VcfBodyShardWriter {
    pub(crate) fn create(path: &std::path::Path, ctx: std::sync::Arc<VcfShardContext>) -> Result<Self> {
        let file = std::fs::File::create(path).map_err(|e| {
            DataFusionError::Execution(format!("Failed to create VCF body shard {}: {e}", path.display()))
        })?;
        Ok(Self { writer: std::io::BufWriter::new(file), ctx, batch_id: 0, lines: 0, input_rows: 0, bytes: 0 })
    }

    /// Format one batch and stream it to the shard file. No buffering of the batch.
    pub(crate) fn write_batch(&mut self, batch: datafusion::arrow::record_batch::RecordBatch) -> Result<()> {
        use std::io::Write;
        let rows = batch.num_rows();
        let formatted = format_vcf_body_chunk(
            self.batch_id,
            batch,
            std::sync::Arc::clone(&self.ctx.vcf_info_fields),
            std::sync::Arc::clone(&self.ctx.unique_format_tags),
            std::sync::Arc::clone(&self.ctx.sample_names),
            self.ctx.coordinate_zero_based,
        )?;
        self.batch_id += 1;
        self.writer.write_all(&formatted.bytes).map_err(|e| {
            DataFusionError::Execution(format!("Failed to write VCF body shard: {e}"))
        })?;
        self.lines += formatted.lines;
        self.input_rows += rows;
        self.bytes += formatted.bytes.len();
        Ok(())
    }

    pub(crate) fn finish(mut self) -> Result<()> {
        use std::io::Write;
        self.writer.flush().map_err(|e| {
            DataFusionError::Execution(format!("Failed to flush VCF body shard: {e}"))
        })
    }
}
```

(Note: `format_vcf_body_chunk` currently takes `vcf_info_fields/unique_format_tags/sample_names` as `Arc<Vec<String>>` — confirm the exact signature at its definition and match it. If it lives in this module already, no import needed.)

- [ ] **Step 2: Refactor `write_vcf_partition_body` to use it (DRY check).**

Replace the inline `BufWriter` + per-batch `format_vcf_body_chunk` + `write_all` loop in `write_vcf_partition_body` (`:320-389`) with `VcfBodyShardWriter::create(...)` + a loop of `write_batch` + `finish`, preserving the `PartitionedVcfBody` descriptor it returns. This proves the extraction is behavior-identical.

- [ ] **Step 3: Build + run the existing vcf_sink tests.**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache --lib vcf 2>&1 | grep -E "test result|error" | head`
Expected: the vcf_sink tests pass (or match the pre-existing baseline set — none of the 35 known failures are vcf-body-writer tests; verify no NEW failure).

- [ ] **Step 4: Commit.**

```bash
git add -A && git commit -m "refactor(vep): extract VcfBodyShardWriter (streaming) reused by write_vcf_partition_body"
```

---

## Task 2: Worker writes its shard instead of sending to a channel

**Files:** `src/annotate_provider.rs` — `spawn_annotation_from_lookup` (`:10019`).

- [ ] **Step 1: Change the worker output side.**

Add a sharded variant (or branch) of `spawn_annotation_from_lookup` that takes the `Arc<VcfShardContext>` + this worker's shard path + `partition_id`, and returns `JoinHandle<Result<PartitionedVcfBody>>` (the shard descriptor) instead of `(Receiver, JoinHandle<()>)`. Replace the output side: where the current code does `for b in out { tx.send(Ok(b)).await }` (`:10065`) and the final flush (`:10078-10097`), instead:

```rust
let mut shard = crate::vcf_sink::VcfBodyShardWriter::create(&shard_path, Arc::clone(&shard_ctx))?;
// in the annotate loop, for each produced batch `b`:
shard.write_batch(b)?;            // streaming format + write, no channel
// after the loop + final partial-window flush:
shard.finish()?;
// return the descriptor (partition_id, path, lines, bytes)
```

Keep the lookup→hydrate→annotate body (`block_in_place(hydrate_worker_window + annotate_worker_window)`) exactly as-is — only the consumption of `out` changes from channel-send to `shard.write_batch`.

- [ ] **Step 2: Build (callers will break — expected; fixed in Task 3).**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error\[" | head`
Expected: errors only at the `ParallelContigState` wiring (Task 3). Compile-unit with Task 3.

- [ ] **Step 3: Commit with Task 3.**

---

## Task 3: Remove the ordered drain; collect shard descriptors

**Files:** `src/annotate_provider.rs` — `ParallelContigState` (`:9725`), the threads>1 spawn block (`:11718-11778`), the `AnnotatingParallel` drain arm (`:11811-11867`).

- [ ] **Step 1: Rework the threads>1 spawn block.**

In the `config.annotation_threads > 1` branch (`:11718`), require `config.vcf_shard_ctx` to be `Some` (sharded output is the only parallel mode now). For each lookup partition handle, call the sharded `spawn_annotation_from_lookup` with `partition_id` and shard path `tempdir.join(format!("partition_{global_id:04}.vcf.body"))` where `global_id` is globally position-ordered (contig-major, worker-minor). Collect the `JoinHandle<Result<PartitionedVcfBody>>`s. Store them in a reworked `ParallelContigState { shard_jobs: Vec<JoinHandle<Result<PartitionedVcfBody>>>, lookup_join, ephemeral_tables, chrom, config, session, t_contig, shared }` (drop `receivers`/`current`/`annotate_join`).

- [ ] **Step 2: Replace the `AnnotatingParallel` drain arm.**

The arm no longer yields row batches. It `await`s all `shard_jobs` (join), collects the shard descriptors into the contig's result, emits the pipeline profile, then transitions to cleanup. Since `poll_next` is sync, drive the joins via a future stored in the state (mirror `AwaitingWindow`): a `JoinHandle`/`JoinSet` future that resolves to `Vec<PartitionedVcfBody>`; on ready, stash the descriptors (see Task 4 for how they reach the sink) and go to `CleaningUp`. Delete the `state.receivers[current]` polling, `current` advance, and `send_wait`/ordered-drain logic for this path.

- [ ] **Step 3: Build.**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error|Finished" | head`
Expected: resolve remaining type errors; `Finished`.

- [ ] **Step 4: Commit (Tasks 2+3).**

```bash
git add -A && git commit -m "feat(vep): fused workers stream VCF shards; remove ordered drain"
```

---

## Task 4: Surface shard descriptors to the sink + final merge

**Files:** `src/annotate_provider.rs` (how the stream reports shards), `src/vcf_sink.rs` (the threads>1 orchestration entry + header/concat).

- [ ] **Step 1: Decide the descriptor channel.**

The sharded annotation produces shard files in `vcf_shard_ctx.tempdir`; the descriptors are `partition_{id}.vcf.body` with deterministic ids. So the sink does **not** need a data stream — after driving the annotation to completion it can enumerate `tempdir` for `partition_*.vcf.body` and concat by id. Implement: the engine output stream for the sharded path yields **zero data batches** (it only drives work to completion); the sink knows the tempdir and the expected id range.

- [ ] **Step 2: Add the sink orchestration for `threads>1`.**

In `annotate_to_vcf`, where it currently gates `if concurrency_plan.chromosome_lanes > 1 { stream_partitioned_vcf_body } else { serial }` (`:1270`), add a branch for `config.threads > 1`:
- Derive the formatter context (it is already computed for the header at `:1214-1255`): build `Arc<VcfShardContext> { vcf_info_fields, unique_format_tags, sample_names, coordinate_zero_based, tempdir }`.
- Construct the annotation plan **directly** (not via SQL) using the reused plan builder (`scan_with_transcript_engine_partitioned`, found in Task 0) with `ContigAnnotationConfig.vcf_shard_ctx = Some(ctx)`.
- Drive it to completion: `let mut s = plan.execute(0, task_ctx)?; while s.next().await.transpose()?.is_some() {}` (no data; just completes the shard writes).
- `writer.write_header(...)` once.
- Concat shards in ascending id order: reuse the `next_write_partition_id` + `copy_body_file_to_writer` loop (`:648-684`) over `tempdir.join(format!("partition_{id:04}.vcf.body"))` for `id in 0..total_shards`. `total_shards` = sum of per-contig worker counts (known after the run, or derivable from the descriptor files present).

- [ ] **Step 3: Build + targeted run.**

Run: `cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1`
Expected: `Finished`.

- [ ] **Step 4: Commit.**

```bash
git add -A && git commit -m "feat(vep): vcf_sink drives sharded annotation + final shard merge (header+concat)"
```

---

## Task 5: fmt + clippy

- [ ] `cargo fmt`; `cargo clippy -p datafusion-bio-function-vep --features lance-cache --lib --example bench_annotate_vcf -- -D warnings 2>&1 | grep -E "error|warning|Finished" | head` → clean. Commit.

---

## Task 6: Build bench + ensure input

- [ ] `cargo build --release --features lance-cache --example bench_annotate_vcf` → Finished.
- [ ] Ensure `/tmp/chr1_200k.vcf.gz` present (bgzip+tabix the chr1 first-200k slice if missing — same repro as prior plans).

---

## Task 7: lib suite — failure set unchanged

- [ ] `cargo test -p datafusion-bio-function-vep --features lance-cache --lib` → 923 pass / 35 pre-existing fail, **0 new** (compare to `/tmp/b.txt` baseline; `comm -13` empty).

---

## Task 8: MANDATORY parity gate — 0 mismatches threads {1,2,4,8} + vs pre-change serial

- [ ] Build the pre-change baseline binary (the commit before Task 0) → `/tmp/sh_c1_baseline.vcf` at threads=1.
- [ ] Produce `/tmp/sh_c{1,2,4,8}.vcf` with the current binary (`--everything` lance chr1 200k).
- [ ] GATE cross-thread: `diff <(grep -v '^#' /tmp/sh_c1.vcf) <(grep -v '^#' /tmp/sh_c{2,4,8}.vcf)` each EMPTY.
- [ ] GATE serial: `diff <(grep -v '^#' /tmp/sh_c1_baseline.vcf) <(grep -v '^#' /tmp/sh_c1.vcf)` EMPTY.
- [ ] `VEP_PROFILE=1` t8 shows `send_wait ≈ 0`.

Any mismatch → STOP, debug with superpowers:systematic-debugging.

---

## Task 9: Benchmark — RSS + wall, threads {1,2,4,8}, A/B vs pre-change

- [ ] Stage before/after binaries (checkout pre-Task-0 `annotate_provider.rs`+`vcf_sink.rs`, build, cp; restore, build, cp).
- [ ] `/usr/bin/time -l` both, threads {1,2,4,8}, min of 2: record wall + peak RSS.
- [ ] Expect: t8 wall improved (send_wait→0 + parallel serialization); RSS no worse. Record results in the spec's Results section; commit.

---

## Self-Review

- **Spec coverage:** §4 worker-writes-shard → Tasks 1-3; §5 reusable serializer + handoff + concat → Tasks 1,4; final merge → Task 4; remove drain → Task 3; gates §9 → Tasks 7-9.
- **Placeholder scan:** Task 0 and Task 4 carry the one genuinely-open piece — the direct-plan-construction handoff (bypassing SQL to inject `vcf_shard_ctx`). Task 0 is the de-risk spike that pins the exact reuse point (`scan_with_transcript_engine_partitioned`); if that function is not cleanly reusable from `vcf_sink`, STOP and revisit the handoff design with the user before Task 4.
- **Type consistency:** `VcfShardContext`, `VcfBodyShardWriter`, `vcf_shard_ctx`, `PartitionedVcfBody` used consistently across tasks. `format_vcf_body_chunk` signature must be confirmed at Task 1 Step 1 and matched.
- **Risk:** largest change of the sequence; the byte-identical gate (Task 8) is the backstop; `threads==1` untouched. If Task 0/4 reveal the SQL-bypass handoff is too invasive, fall back to the A1 interim (workers write Arrow-IPC shards, single serialize) — but that loses parallel serialization; raise with the user.
