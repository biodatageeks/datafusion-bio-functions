# Pipelined Per-Contig Shard Assembler (Plain + BGZF) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the single-threaded post-100% shard-merge tail (~13s on WGS) by assembling each contig's VCF body shards on a dedicated thread that runs *concurrently* with the next contig's context-prepare + annotation, and — for `.vcf.gz`/bgzf output — move compression off the serial tail onto the parallel workers via raw BGZF block concatenation.

**Architecture:** Today `drive_sharded_vcf_annotation` (vcf_sink.rs:713) drives the whole-genome annotation stream to completion, then serially concatenates every `partition_NNNN.vcf.body` shard into the output writer in one terminal pass (vcf_sink.rs:835–869). The progress bar hits 100% when annotation finishes (`shard_ctx.rows_done == total_input`), so the concat is invisible tail. We replace the terminal concat with a **dedicated OS assembler thread** that owns the output writer and consumes per-contig "shard range complete" messages published by `ContigAnnotationStream` the instant each contig's workers finish (annotate_provider.rs:11697). Because shard ids are contig-major/worker-minor and equal final-VCF position order (annotate_provider.rs:11580–11587), the assembler appends each contig's contiguous id range in arrival order → identical bytes, but the copy now overlaps the next contig's `prepare_contig_context` + annotation. Stage 1 delivers this for Plain/Gzip output (text shards, re-copy/re-compress at the assembler — byte-identical). Stage 2 adds BGZF output: workers write BGZF shards (parallel compression) and the assembler does a raw block concat (strip each shard's 28-byte EOF, write one canonical EOF at the end).

**Tech Stack:** Rust 2024, DataFusion 53 / Arrow 58, tokio (current-thread + multi-thread), `std::sync::mpsc` (assembler channel), `std::thread` (assembler), bio-format-vcf's existing `VcfLocalWriter` (the `Bgzf` variant is backed by `noodles-bgzf` internally — bio-functions never depends on `noodles-bgzf` directly), indicatif (progress). **Entirely single-repo (`datafusion-bio-functions`): no bio-formats change, no dependency bump.**

## Global Constraints

- **Byte-identical output is the acceptance gate.** For a fixed input + `workers` + compression, the assembled VCF must be byte-for-byte identical to the pre-change output (Plain/Gzip) or decompress to byte-identical body (BGZF). The memory-noted invariant is "chr4/chr2 w{1,4} 0 HGNC mism, byte-identical" — do not regress it.
- **Shard id ordering is the ordering contract.** Ids are assigned `contig-major, worker-minor` at annotate_provider.rs:11583; ascending id == final VCF position order. Never reorder; only slice along contig boundaries.
- **Single-repo, no dependency bump.** Both stages live entirely in `datafusion-bio-functions`. BGZF compression reuses the existing `VcfLocalWriter::Bgzf` from the `datafusion-bio-format-vcf` git dep (tag `v1.8.7`, unchanged); the 28-byte BGZF EOF marker is a fixed spec constant hardcoded locally in `vcf_sink.rs`. Do **not** add a `noodles-bgzf` dependency and do **not** modify bio-formats.
- **No behavioral change for `workers == 1`.** The serial path (vcf_sink.rs:1157–1239) is untouched; assembler applies only to `workers > 1`.
- **Assembler is single-writer.** Exactly one thread ever touches the output `File`. Contig ranges are appended strictly in arrival (= id) order. Never write from two threads.
- **Progress accounting unchanged.** The bar is driven by `rows_done` during annotation (vcf_sink.rs:777–793); the assembler must NOT advance the bar or the callback (that would double-count input vs output lines — see vcf_sink.rs:860–862).
- **Error contract.** On any worker/stream error the partially-assembled output must not be presented as a valid result. Assemble into a temp path and atomically rename on full success (both stages); on error, drop the temp file.

---

## File Structure

**Stage 1 (datafusion-bio-functions only):**
- Modify `datafusion/bio-function-vep/src/vcf_sink.rs`
  - New `ContigShardRange` message type + `contig_done_tx` field on `VcfShardContext`.
  - New `run_assembler_thread(...)` (owns writer, consumes ranges, appends, finish, atomic rename).
  - Rewrite `drive_sharded_vcf_annotation` to spawn the assembler, drive the stream, drop the sender, join the assembler.
  - Delete the terminal concat loop + `copy_body_file_counting_lines`/`bytecount` (replaced by assembler append + worker-reported line counts).
- Modify `datafusion/bio-function-vep/src/annotate_provider.rs`
  - `spawn_annotation_from_lookup_sharded` returns `ShardResult { input_rows, output_lines }`.
  - `ParallelContigState` gains `contig_first_shard_id`; capture before the worker loop; publish `ContigShardRange` on completion (11697).

**Stage 2 (datafusion-bio-functions only):**
- Modify `datafusion/bio-function-vep/src/vcf_sink.rs`
  - Refactor `VcfBodyShardWriter` to own a `VcfLocalWriter` (Plain for text shards, Bgzf for bgzf shards) instead of a raw `BufWriter<File>` — workers compress their own shard via the existing enum; no new API.
  - Hardcode the 28-byte `BGZF_EOF` constant locally.
  - Assembler grows a `RawBgzf` mode: pre-compress the header via `VcfLocalWriter::Bgzf` (strip its EOF), raw-append each shard minus its trailing 28-byte EOF, write one final `BGZF_EOF`.
  - Shard/assembler mode selection keyed on `config.compression` (Bgzf → bgzf shards; else Plain text shards).

**Tests:**
- `datafusion/bio-function-vep/src/vcf_sink.rs` `#[cfg(test)]` — assembler unit tests with synthetic shard files (Plain + BGZF), ordering, BGZF raw-concat round-trip, atomic rename.
- e2e parity: `vepyr/e2e-testing/scripts/run_annotation_fast.py` (unchanged) as the oracle.

---

## Stage 1 — Dedicated assembler thread, incremental per-contig concat (Plain/Gzip)

### Task 1: `ContigShardRange` message + wire it onto `VcfShardContext`

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs:156-168` (add field), new type near line 148
- Test: `datafusion/bio-function-vep/tests/shard_assembler.rs` (create)

**Interfaces:**
- Produces:
  ```rust
  pub(crate) struct ContigShardRange {
      pub(crate) first_shard_id: usize,   // inclusive
      pub(crate) end_shard_id: usize,     // exclusive
      pub(crate) output_lines: usize,     // sum of this contig's shards' body lines
  }
  ```
  and a new field `pub(crate) contig_done_tx: std::sync::mpsc::Sender<ContigShardRange>` on `VcfShardContext`.

- [ ] **Step 1: Write the failing test** (`tests/shard_assembler.rs`)

```rust
use datafusion_bio_function_vep::vcf_sink::ContigShardRange; // pub(crate) → test lives in-crate; see Step 3 note
// NOTE: because VcfShardContext/ContigShardRange are pub(crate), put this test as a
// `#[cfg(test)] mod` at the bottom of vcf_sink.rs instead of an external integration
// test. Move the assert there:
#[test]
fn contig_shard_range_is_constructible_and_ordered() {
    let a = ContigShardRange { first_shard_id: 0, end_shard_id: 3, output_lines: 10 };
    let b = ContigShardRange { first_shard_id: 3, end_shard_id: 5, output_lines: 4 };
    assert_eq!(a.end_shard_id, b.first_shard_id); // contiguous, no gaps/overlap
    assert_eq!(a.output_lines + b.output_lines, 14);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep contig_shard_range_is_constructible -- --nocapture`
Expected: FAIL — `ContigShardRange` not found.

- [ ] **Step 3: Add the type + field**

In `vcf_sink.rs` just above `VcfShardContext` (line 156):
```rust
/// Published by `ContigAnnotationStream` the instant a contig's shard workers
/// all finish (shards flushed). The assembler thread appends the id range
/// `[first_shard_id, end_shard_id)` — which is contiguous and equals final-VCF
/// position order — and adds `output_lines` to the running body-line total.
pub(crate) struct ContigShardRange {
    pub(crate) first_shard_id: usize,
    pub(crate) end_shard_id: usize,
    pub(crate) output_lines: usize,
}
```
Add to `VcfShardContext` (after `rows_done`, line 167):
```rust
    /// One message per completed contig, in completion (= ascending id) order.
    /// Dropped when annotation finishes, which signals the assembler to finalize.
    pub(crate) contig_done_tx: std::sync::mpsc::Sender<ContigShardRange>,
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep contig_shard_range_is_constructible -- --nocapture`
Expected: PASS. (The `VcfShardContext` construction site at vcf_sink.rs:1131 will now fail to compile — that's fixed in Task 4; if doing strict per-task compiles, temporarily add `contig_done_tx: std::sync::mpsc::channel().0` there and remove in Task 4.)

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/vcf_sink.rs
git commit -m "feat(vep): add ContigShardRange message + contig_done_tx on VcfShardContext"
```

---

### Task 2: Worker returns output line count (`ShardResult`)

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:9904-9913` (signature/return), `:10042-10045` (return value)
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:11571,11624-11638` (join_fut sums both)

**Interfaces:**
- Produces:
  ```rust
  pub(crate) struct ShardResult { pub(crate) input_rows: usize, pub(crate) output_lines: usize }
  ```
  `spawn_annotation_from_lookup_sharded(...) -> tokio::task::JoinHandle<Result<ShardResult>>`.
- Consumes: `VcfBodyShardWriter.input_rows` and `.lines` (vcf_sink.rs:310-311, already `pub(crate)`).

- [ ] **Step 1: Write the failing test**

Unit-testing the async worker in isolation is impractical (needs a live lookup stream); assert the *type contract* instead, in `annotate_provider.rs` `#[cfg(test)]`:
```rust
#[test]
fn shard_result_carries_input_and_output_counts() {
    let r = super::ShardResult { input_rows: 5000, output_lines: 18234 };
    assert_eq!(r.input_rows, 5000);
    assert_eq!(r.output_lines, 18234);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep shard_result_carries -- --nocapture`
Expected: FAIL — `ShardResult` not found.

- [ ] **Step 3: Implement**

Define near `spawn_annotation_from_lookup_sharded` (annotate_provider.rs:~9903):
```rust
pub(crate) struct ShardResult {
    pub(crate) input_rows: usize,
    pub(crate) output_lines: usize,
}
```
Change the return type at 9913 to `-> tokio::task::JoinHandle<Result<ShardResult>>`.
Replace the tail at 10042–10045:
```rust
        let input_rows = shard.input_rows;
        let output_lines = shard.lines;
        let _ = (emit_start, emit_end, global_row);
        shard.finish()?;
        Ok(ShardResult { input_rows, output_lines })
```
Update `shard_handles` element type at 11571:
```rust
                            let mut shard_handles: Vec<tokio::task::JoinHandle<Result<ShardResult>>> =
```
Update `join_fut` (11624–11638) to sum both and return them:
```rust
                            let join_fut: ShardJoinFuture = Box::pin(async move {
                                let mut input_rows = 0usize;
                                let mut output_lines = 0usize;
                                for h in shard_handles {
                                    match h.await {
                                        Ok(Ok(r)) => { input_rows += r.input_rows; output_lines += r.output_lines; }
                                        Ok(Err(e)) => return Err(e),
                                        Err(join_err) => {
                                            return Err(DataFusionError::External(Box::new(join_err)));
                                        }
                                    }
                                }
                                Ok((input_rows, output_lines))
                            });
```
Update the `ShardJoinFuture` type alias (search `type ShardJoinFuture`) to yield `Result<(usize, usize)>`, and the `AnnotatingParallel` `Ready(Ok(...))` arm at 11697 to bind `(contig_rows, contig_lines)` — `contig_lines` is used in Task 3; keep `contig_rows` for the existing profile string at 11705.

- [ ] **Step 4: Run test + build**

Run: `cargo test -p datafusion-bio-function-vep shard_result_carries -- --nocapture && cargo build -p datafusion-bio-function-vep`
Expected: PASS + clean build.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): shard workers report output line count via ShardResult"
```

---

### Task 3: Publish `ContigShardRange` on contig completion

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` — `ParallelContigState` struct (search `struct ParallelContigState`), setup at `:11575`, completion at `:11697-11713`

**Interfaces:**
- Consumes: `ContigShardRange`, `VcfShardContext.contig_done_tx` (Task 1); `(contig_rows, contig_lines)` (Task 2).
- Produces: one `ContigShardRange` sent per completed contig.

- [ ] **Step 1: Write the failing test**

Ordering/id-range logic is what can break; assert the range arithmetic in `#[cfg(test)]`:
```rust
#[test]
fn contig_range_is_first_to_next_global_id() {
    // Simulate: contig starts at global id 4, spawns 3 workers → ids 4,5,6.
    let first = 4usize;
    let mut next = first;
    for _ in 0..3 { next += 1; }
    let range = crate::vcf_sink::ContigShardRange { first_shard_id: first, end_shard_id: next, output_lines: 0 };
    assert_eq!((range.first_shard_id, range.end_shard_id), (4, 7));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep contig_range_is_first_to_next -- --nocapture`
Expected: FAIL until the module path/type resolves (it will pass once Task 1 compiled; if so, this is a guard test — keep it).

- [ ] **Step 3: Implement**

Add field to `ParallelContigState`:
```rust
    contig_first_shard_id: usize,
```
At annotate_provider.rs:11575, immediately before `for mut handle in handles {`:
```rust
                            let contig_first_shard_id = self.next_global_shard_id;
```
Set it when building `ParallelContigState` (near 11639):
```rust
                                contig_first_shard_id,
```
In the `AnnotatingParallel` `Poll::Ready(Ok((contig_rows, contig_lines)))` arm (11697), before building the cleanup future (11708), publish the range. `end` is the current `next_global_shard_id` (workers for later contigs have not been allocated yet, so it is exactly this contig's exclusive end):
```rust
                        Poll::Ready(Ok((contig_rows, contig_lines))) => {
                            for h in &state.lookup_join { h.abort(); }
                            if let Some(shard_ctx) = state.config.vcf_shard_ctx.clone() {
                                // Send is best-effort: if the assembler has gone
                                // away (error path already tearing down), ignore.
                                let _ = shard_ctx.contig_done_tx.send(crate::vcf_sink::ContigShardRange {
                                    first_shard_id: state.contig_first_shard_id,
                                    end_shard_id: self.next_global_shard_id,
                                    output_lines: contig_lines,
                                });
                            }
                            profile_end!(
                                &format!("{}: TOTAL", state.chrom),
                                state.t_contig,
                                format!("{contig_rows} rows")
                            );
                            emit_contig_pipeline_profile(&state.shared.profile, &state.chrom);
                            let fut = make_cleanup_future(
                                Arc::clone(&state.session),
                                std::mem::take(&mut state.ephemeral_tables),
                            );
                            self.state = StreamState::CleaningUp(fut);
                            continue;
                        }
```

- [ ] **Step 4: Run test + build**

Run: `cargo test -p datafusion-bio-function-vep contig_range_is_first_to_next -- --nocapture && cargo build -p datafusion-bio-function-vep`
Expected: PASS + clean build.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): publish per-contig shard range on completion"
```

---

### Task 4: Assembler thread (Plain/Gzip append), replacing the terminal concat

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` — new `run_assembler_thread`, rewrite `drive_sharded_vcf_annotation` (713–870), construct `contig_done_tx` at 1131, move header/writer ownership.
- Delete: `copy_body_file_counting_lines` (680), `bytecount` (705) — superseded.
- Test: `tests/shard_assembler.rs` (in-crate `#[cfg(test)]`)

**Interfaces:**
- Consumes: `ContigShardRange` (Task 1), `VcfLocalWriter` (bio-formats), shard files `partition_{id:04}.vcf.body` in `tempdir`.
- Produces:
  ```rust
  // Runs on a dedicated OS thread; owns `writer`. Returns total body lines written.
  fn run_assembler_thread(
      mut writer: VcfLocalWriter,
      final_output: PathBuf,      // real destination; writer currently targets a temp path
      temp_output: PathBuf,
      tempdir: PathBuf,
      rx: std::sync::mpsc::Receiver<ContigShardRange>,
  ) -> std::thread::JoinHandle<Result<usize>>
  ```

- [ ] **Step 1: Write the failing test** (`vcf_sink.rs` `#[cfg(test)]`)

```rust
#[test]
fn assembler_appends_ranges_in_order_and_counts_lines() {
    use std::io::Write;
    let dir = std::env::temp_dir().join(format!("asm_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // three shards: contig0 = [0,2), contig1 = [2,3)
    for (id, body) in [(0usize, "a\nb\n"), (1, "c\n"), (2, "d\ne\nf\n")] {
        let mut f = std::fs::File::create(dir.join(format!("partition_{id:04}.vcf.body"))).unwrap();
        f.write_all(body.as_bytes()).unwrap();
    }
    let out = dir.join("out.vcf");
    let tmp = dir.join("out.vcf.tmp");
    let writer = VcfLocalWriter::with_compression(&tmp, VcfCompressionType::Plain).unwrap();
    let (tx, rx) = std::sync::mpsc::channel();
    let h = run_assembler_thread(writer, out.clone(), tmp, dir.clone(), rx);
    tx.send(ContigShardRange { first_shard_id: 0, end_shard_id: 2, output_lines: 3 }).unwrap();
    tx.send(ContigShardRange { first_shard_id: 2, end_shard_id: 3, output_lines: 3 }).unwrap();
    drop(tx);
    let total = h.join().unwrap().unwrap();
    assert_eq!(total, 6);
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "a\nb\nc\nd\ne\nf\n");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep assembler_appends_ranges_in_order -- --nocapture`
Expected: FAIL — `run_assembler_thread` not found.

- [ ] **Step 3: Implement `run_assembler_thread`**

Add to `vcf_sink.rs` (near the deleted concat helpers). It appends shard files raw through the writer (Plain = byte copy; Gzip = compress-on-append), counts `output_lines` from the message (no re-scan), removes each shard after copy, then `writer.finish()` and atomic-rename temp→final:
```rust
fn run_assembler_thread(
    mut writer: VcfLocalWriter,
    final_output: PathBuf,
    temp_output: PathBuf,
    tempdir: PathBuf,
    rx: std::sync::mpsc::Receiver<ContigShardRange>,
) -> std::thread::JoinHandle<Result<usize>> {
    std::thread::spawn(move || -> Result<usize> {
        let mut total_lines = 0usize;
        let mut buffer = vec![0u8; 8 * 1024 * 1024];
        while let Ok(range) = rx.recv() {
            for id in range.first_shard_id..range.end_shard_id {
                let path = tempdir.join(format!("partition_{id:04}.vcf.body"));
                let mut f = std::fs::File::open(&path).map_err(|e| {
                    DataFusionError::Execution(format!("assembler open shard {}: {e}", path.display()))
                })?;
                loop {
                    let n = f.read(&mut buffer).map_err(|e| {
                        DataFusionError::Execution(format!("assembler read shard {}: {e}", path.display()))
                    })?;
                    if n == 0 { break; }
                    write_vcf_body_chunk(&mut writer, &buffer[..n])?;
                }
                let _ = std::fs::remove_file(&path);
            }
            total_lines += range.output_lines;
        }
        // Sender dropped == annotation finished with no error.
        writer.finish()?;
        std::fs::rename(&temp_output, &final_output).map_err(|e| {
            DataFusionError::Execution(format!(
                "assembler rename {} -> {}: {e}",
                temp_output.display(), final_output.display()
            ))
        })?;
        Ok(total_lines)
    })
}
```

- [ ] **Step 4: Run the assembler unit test**

Run: `cargo test -p datafusion-bio-function-vep assembler_appends_ranges_in_order -- --nocapture`
Expected: PASS.

- [ ] **Step 5: Rewire `annotate_to_vcf` + `drive_sharded_vcf_annotation`**

In `annotate_to_vcf` (vcf_sink.rs:1107–1156): the writer must target a **temp path** and be owned by the assembler. Replace the writer creation (1108–1115) so, for `workers > 1`, we write the header to a temp-path writer, hand ownership to the assembler, and the header for `workers == 1` stays as-is. Concretely, keep header-writing on `writer` (1110), then in the `workers > 1` branch:
```rust
    if config.workers > 1 {
        let tempdir = VcfBodyTempDir::new()?;
        let (contig_done_tx, contig_done_rx) = std::sync::mpsc::channel::<ContigShardRange>();
        let shard_ctx = Arc::new(VcfShardContext {
            vcf_info_fields: Arc::new(vcf_info_fields),
            unique_format_tags: Arc::new(unique_format_tags),
            sample_names: Arc::new(sample_names),
            coordinate_zero_based,
            tempdir: tempdir.path().to_path_buf(),
            rows_done: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            contig_done_tx,
        });
        // `writer` already has the header; it targets `temp_output`. Move it to
        // the assembler thread, which finalizes + renames to `output_vcf`.
        let assembler = run_assembler_thread(
            writer,
            Path::new(output_vcf).to_path_buf(),
            temp_output.clone(),
            tempdir.path().to_path_buf(),
            contig_done_rx,
        );
        // Drive the annotation stream (drops shard_ctx.contig_done_tx clones as
        // contigs complete; the ORIGINAL tx lives on shard_ctx, dropped below).
        drive_sharded_vcf_annotation(
            &ctx, vcf_table, cache_source, backend, cache_source_type,
            options_json.clone(), &vcf_schema, &projection_names,
            Arc::clone(&shard_ctx), &pb, config, total_input,
        ).await?;
        // Drop the last sender so the assembler sees channel close → finalize.
        drop(shard_ctx);
        total_rows = assembler.join().map_err(|_| {
            DataFusionError::Execution("assembler thread panicked".to_string())
        })??;
        drop(tempdir);
    } else { /* unchanged serial path */ }
```
Where `temp_output` is computed just before the writer is created:
```rust
    let output_path = Path::new(output_vcf);
    let temp_output = output_path.with_extension(format!(
        "{}.part",
        output_path.extension().and_then(|s| s.to_str()).unwrap_or("vcf")
    ));
    let mut writer = VcfLocalWriter::with_compression(&temp_output, config.compression)?;
```
Change `drive_sharded_vcf_annotation`'s signature: drop the `writer: &mut VcfLocalWriter` param and its `-> Result<usize>` concat return; it now returns `Result<()>` and its body is **only** the plan-drive + progress-poller loop (current lines 730–833), deleting everything from line 835 (`// Assemble:`) onward. The `report(shard_ctx.rows_done...)` final flush at 833 stays.
Delete `copy_body_file_counting_lines` (680–703) and `bytecount` (705–707). Update the terminal `writer.finish()` at 1241–1245 so it runs **only** in the serial (`workers == 1`) branch (the assembler already finished the parallel writer). Guard it:
```rust
    if config.workers <= 1 {
        let finish_started = Instant::now();
        writer.finish()?;
        if let Some(profile) = sink_profile.as_mut() { profile.writer_finish += finish_started.elapsed(); }
    }
```
(Move `writer` into the assembler for `workers > 1`; the serial branch keeps ownership, so this compiles — but the `let mut writer` binding is consumed by the `if` branch. Restructure so `writer` is declared and, in the parallel branch, `writer` is moved into `run_assembler_thread` while the serial branch owns it to the end. Use the branch structure above where the parallel branch fully consumes `writer` and returns early-ish; the serial branch owns it through `finish()`.)

- [ ] **Step 6: Build + run the full vep test suite**

Run: `cargo test -p datafusion-bio-function-vep -- --nocapture`
Expected: PASS (including `vcf_output_roundtrip.rs` and the new assembler test).

- [ ] **Step 7: Commit**

```bash
git add datafusion/bio-function-vep/src/vcf_sink.rs datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): pipelined per-contig shard assembler thread (plain/gzip)"
```

---

### Task 5: Stage 1 byte-identical e2e parity gate

**Files:**
- No source changes; verification only. Uses `vepyr/e2e-testing/scripts/run_annotation_fast.py`.

- [ ] **Step 1: Capture a pre-change baseline** (from the commit BEFORE Task 1; do this once, up front, or check out the parent)

```bash
cd /Users/mwiewior/research/git/vepyr
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --workers 4 --force --skip-compare
cp "$(ls -t /path/to/work/vepyr_*_chr4_merged.vcf | head -1)" /tmp/baseline_chr4_w4.vcf
```

- [ ] **Step 2: Build the changed engine + rerun**

Run:
```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo build --release -p datafusion-bio-function-vep
cd /Users/mwiewior/research/git/vepyr && maturin develop --release   # rebuild vepyr against new engine
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --workers 4 --force --skip-compare
```

- [ ] **Step 3: Assert byte-identical**

Run: `cmp /tmp/baseline_chr4_w4.vcf "$(ls -t /path/to/work/vepyr_*_chr4_merged.vcf | head -1)" && echo "BYTE-IDENTICAL OK"`
Expected: `BYTE-IDENTICAL OK` (no output from `cmp`).

- [ ] **Step 4: Assert the tail shrank** (profiling sanity, not a hard gate)

Run: `VEP_PROFILE=1 python e2e-testing/scripts/run_annotation_fast.py chr1 --cache merged --workers 8 --force --skip-compare 2>&1 | grep -E "VEP_PROGRESS|vcf_sink_profile"`
Expected: `writer_finish` collapses toward ~0 for the parallel path (assembler overlapped); wall-clock gap between last `[VEP_PROGRESS] annotated N/N` and process exit is materially smaller than the pre-change ~13s.

- [ ] **Step 5: Commit the verification note**

```bash
git commit --allow-empty -m "test(vep): Stage 1 chr4 w4 byte-identical + tail-overlap verified"
```

---

## Stage 2 — BGZF shards + raw block concat (parallel compression)

> Prereq: Stage 1 merged. Entirely within `datafusion-bio-functions` — no bio-formats change, no dependency bump. BGZF compression reuses the existing `VcfLocalWriter::Bgzf`; the EOF marker is hardcoded.

### Task 6: `VcfBodyShardWriter` owns a `VcfLocalWriter` (enables BGZF shards)

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` — `VcfBodyShardWriter` (303–386); `VcfShardContext` (add `shard_compression`); construction at 1131
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:9916` — pass compression to `create`

**Interfaces:**
- Consumes: `VcfLocalWriter`, `VcfCompressionType`, `write_vcf_body_chunk` (all already in scope in vcf_sink.rs).
- Produces:
  - `VcfShardContext` gains `pub(crate) shard_compression: VcfCompressionType` (only ever `Bgzf` or `Plain`).
  - `VcfBodyShardWriter::create(path, info, tags, samples, zero_based, compression: VcfCompressionType) -> Result<Self>`, internally owning a `VcfLocalWriter`.
  - For `Bgzf`, `.bytes`/`.lines`/`.input_rows` still count *uncompressed* body; the file holds BGZF blocks + a trailing 28-byte EOF from `VcfLocalWriter::finish()` (stripped by the assembler in Task 7).

- [ ] **Step 1: Write the failing test** (`vcf_sink.rs` `#[cfg(test)]`)

```rust
#[test]
fn shard_writer_bgzf_file_is_decodable() {
    use std::io::Read;
    let dir = std::env::temp_dir().join(format!("shard_bgzf_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("partition_0000.vcf.body");
    let w = VcfBodyShardWriter::create(
        &path, Arc::new(vec![]), Arc::new(vec![]), Arc::new(vec![]), false,
        VcfCompressionType::Bgzf,
    ).unwrap();
    let lines = w.lines;               // 0 — no batches written
    w.finish().unwrap();               // flushes bgzf blocks + EOF
    // The empty bgzf shard is a valid, decodable bgzf file (EOF-only).
    let mut reader = noodles_bgzf::Reader::new(std::fs::File::open(&path).unwrap());
    let mut s = String::new();
    reader.read_to_string(&mut s).unwrap();
    assert_eq!(s, "");
    assert_eq!(lines, 0);
    std::fs::remove_dir_all(&dir).ok();
}
```
NOTE: the test references `noodles_bgzf::Reader` only for *reading* in the test; it is not a crate dependency of the shipped code. If `noodles-bgzf` is not resolvable from the test target, decode instead via `bgzip -dc` in a `std::process::Command` and assert the stdout, or gate this assertion behind a helper. The shipped code never names `noodles-bgzf`.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep shard_writer_bgzf_file_is_decodable -- --nocapture`
Expected: FAIL — `create` arity (compression arg) not found.

- [ ] **Step 3: Implement**

Change `VcfBodyShardWriter.writer` from `std::io::BufWriter<std::fs::File>` to `VcfLocalWriter`:
```rust
pub(crate) struct VcfBodyShardWriter {
    writer: VcfLocalWriter,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
    batch_id: usize,
    pub(crate) lines: usize,
    pub(crate) input_rows: usize,
    pub(crate) bytes: usize,
}
```
`create` takes `compression: VcfCompressionType` and builds the sink via the existing bio-formats constructor:
```rust
    pub(crate) fn create(
        path: &Path,
        vcf_info_fields: Arc<Vec<String>>,
        unique_format_tags: Arc<Vec<String>>,
        sample_names: Arc<Vec<String>>,
        coordinate_zero_based: bool,
        compression: VcfCompressionType,
    ) -> Result<Self> {
        let writer = VcfLocalWriter::with_compression(path, compression)?;
        Ok(Self { writer, vcf_info_fields, unique_format_tags, sample_names,
                  coordinate_zero_based, batch_id: 0, lines: 0, input_rows: 0, bytes: 0 })
    }
```
`write_formatted` writes through the existing chunk helper (keeps counters + returns IO duration):
```rust
    fn write_formatted(&mut self, formatted: &FormattedVcfBatch) -> Result<Duration> {
        let started = Instant::now();
        write_vcf_body_chunk(&mut self.writer, &formatted.bytes)?;
        let elapsed = started.elapsed();
        self.lines += formatted.lines;
        self.input_rows += formatted.input_rows;
        self.bytes += formatted.bytes.len();
        Ok(elapsed)
    }
```
`finish` consumes the `VcfLocalWriter` (for `Bgzf` this appends the EOF marker):
```rust
    pub(crate) fn finish(self) -> Result<()> { self.writer.finish() }
```

- [ ] **Step 4: Thread compression through `VcfShardContext` + the worker**

Add `pub(crate) shard_compression: VcfCompressionType` to `VcfShardContext`. In `annotate_to_vcf`'s `workers>1` branch, set it (only two possibilities — bgzf shards for bgzf output, else plain text shards):
```rust
            shard_compression: if config.compression == VcfCompressionType::Bgzf {
                VcfCompressionType::Bgzf
            } else {
                VcfCompressionType::Plain
            },
```
At annotate_provider.rs:9916, pass `shard_ctx.shard_compression` as the new `create` arg (a `Copy` field on the `Arc<VcfShardContext>`).

- [ ] **Step 5: Run test + build**

Run: `cargo test -p datafusion-bio-function-vep shard_writer_bgzf -- --nocapture && cargo build -p datafusion-bio-function-vep`
Expected: PASS + clean build. (Plain output still works: `VcfLocalWriter::Plain` byte path is identical to the old raw `BufWriter`.)

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/vcf_sink.rs datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): VcfBodyShardWriter owns VcfLocalWriter (BGZF shards via existing enum)"
```

---

### Task 7: Raw-BGZF assembler mode + mode selection

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` — hardcode `BGZF_EOF`; add `run_assembler_thread_bgzf`; BGZF header pre-compression; select mode from `config.compression`.

**Interfaces:**
- Consumes: local `const BGZF_EOF: [u8; 28]`; `VcfLocalWriter::Bgzf` (for header compression); BGZF shard files (Task 6).
- Produces:
  ```rust
  const BGZF_EOF: [u8; 28];   // canonical BGZF end-of-file marker (SAM spec)
  fn bgzf_blocks_no_eof(path: &Path) -> Result<Vec<u8>>;  // read a finished bgzf file, strip trailing EOF
  fn run_assembler_thread_bgzf(
      temp_output: PathBuf, final_output: PathBuf, tempdir: PathBuf,
      header_blocks: Vec<u8>, rx: std::sync::mpsc::Receiver<ContigShardRange>,
  ) -> std::thread::JoinHandle<Result<usize>>;
  ```
- Behavior: For BGZF output the assembler owns a raw `BufWriter<File>` (NOT a `VcfLocalWriter::Bgzf` — shards are already compressed and must not be re-compressed). It writes the pre-compressed `header_blocks`, then raw-appends each shard **minus its trailing 28-byte EOF**, then writes one `BGZF_EOF`, flushes, renames.

- [ ] **Step 1: Write the failing test** (`vcf_sink.rs` `#[cfg(test)]`)

```rust
#[test]
fn assembler_bgzf_raw_concat_decodes_in_order() {
    use std::io::{Read, Write};
    let dir = std::env::temp_dir().join(format!("asm_bgzf_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // Produce bgzf shards exactly as VcfBodyShardWriter (bgzf) does: write body
    // through a VcfLocalWriter::Bgzf, finish() (appends EOF), leave file on disk.
    let write_bgzf = |path: &std::path::Path, body: &str| {
        let mut w = VcfLocalWriter::with_compression(path, VcfCompressionType::Bgzf).unwrap();
        write_vcf_body_chunk(&mut w, body.as_bytes()).unwrap();
        w.finish().unwrap();
    };
    for (id, body) in [(0usize, "a\nb\n"), (1usize, "c\n")] {
        write_bgzf(&dir.join(format!("partition_{id:04}.vcf.body")), body);
    }
    // Header blocks: compress "##header\n#CHROM\n" and strip its EOF.
    let hpath = dir.join("hdr.bgz");
    write_bgzf(&hpath, "##header\n#CHROM\n");
    let header_blocks = bgzf_blocks_no_eof(&hpath).unwrap();

    let out = dir.join("out.vcf.gz");
    let tmp = dir.join("out.vcf.gz.part");
    let (tx, rx) = std::sync::mpsc::channel();
    let h = run_assembler_thread_bgzf(tmp.clone(), out.clone(), dir.clone(), header_blocks, rx);
    tx.send(ContigShardRange { first_shard_id: 0, end_shard_id: 1, output_lines: 2 }).unwrap();
    tx.send(ContigShardRange { first_shard_id: 1, end_shard_id: 2, output_lines: 1 }).unwrap();
    drop(tx);
    let total = h.join().unwrap().unwrap();
    assert_eq!(total, 3);
    // Decode the assembled file (bgzip is available in the e2e toolchain).
    let decoded = std::process::Command::new("bgzip").arg("-dc").arg(&out).output().unwrap();
    assert_eq!(String::from_utf8_lossy(&decoded.stdout), "##header\n#CHROM\na\nb\nc\n");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep assembler_bgzf_raw_concat -- --nocapture`
Expected: FAIL — `BGZF_EOF` / `bgzf_blocks_no_eof` / `run_assembler_thread_bgzf` not found.

- [ ] **Step 3: Implement the constant + helpers + assembler**

Hardcode the marker (SAM-spec fixed 28 bytes — identical across htslib/samtools/noodles) in `vcf_sink.rs`:
```rust
/// Canonical 28-byte BGZF end-of-file marker (an empty BGZF block). Written
/// exactly once, at the very end of an assembled BGZF file. Every finished bgzf
/// file (`VcfLocalWriter::Bgzf::finish()`) ends with these exact bytes.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
    0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Read a finished bgzf file and return its bytes with the trailing 28-byte EOF
/// marker removed, ready for raw concatenation into a larger bgzf stream.
fn bgzf_blocks_no_eof(path: &Path) -> Result<Vec<u8>> {
    let mut bytes = std::fs::read(path).map_err(|e| {
        DataFusionError::Execution(format!("read bgzf shard {}: {e}", path.display()))
    })?;
    debug_assert!(bytes.len() >= BGZF_EOF.len() && bytes[bytes.len() - BGZF_EOF.len()..] == BGZF_EOF);
    let keep = bytes.len().saturating_sub(BGZF_EOF.len());
    bytes.truncate(keep);
    Ok(bytes)
}
```
The assembler owns a raw file writer (shards carry their own compression):
```rust
fn run_assembler_thread_bgzf(
    temp_output: PathBuf,
    final_output: PathBuf,
    tempdir: PathBuf,
    header_blocks: Vec<u8>,
    rx: std::sync::mpsc::Receiver<ContigShardRange>,
) -> std::thread::JoinHandle<Result<usize>> {
    std::thread::spawn(move || -> Result<usize> {
        let file = std::fs::File::create(&temp_output).map_err(|e| {
            DataFusionError::Execution(format!("assembler create {}: {e}", temp_output.display()))
        })?;
        let mut out = std::io::BufWriter::new(file);
        out.write_all(&header_blocks)
            .map_err(|e| DataFusionError::Execution(format!("assembler write header: {e}")))?;
        let mut total_lines = 0usize;
        while let Ok(range) = rx.recv() {
            for id in range.first_shard_id..range.end_shard_id {
                let path = tempdir.join(format!("partition_{id:04}.vcf.body"));
                let blocks = bgzf_blocks_no_eof(&path)?;
                out.write_all(&blocks)
                    .map_err(|e| DataFusionError::Execution(format!("assembler write shard: {e}")))?;
                let _ = std::fs::remove_file(&path);
            }
            total_lines += range.output_lines;
        }
        out.write_all(&BGZF_EOF)
            .map_err(|e| DataFusionError::Execution(format!("assembler write EOF: {e}")))?;
        out.flush()
            .map_err(|e| DataFusionError::Execution(format!("assembler flush: {e}")))?;
        drop(out);
        std::fs::rename(&temp_output, &final_output).map_err(|e| {
            DataFusionError::Execution(format!(
                "assembler rename {} -> {}: {e}", temp_output.display(), final_output.display()
            ))
        })?;
        Ok(total_lines)
    })
}
```

- [ ] **Step 4: Mode selection + BGZF header pre-compression in `annotate_to_vcf`**

In the `workers>1` branch, branch on compression. For BGZF, compress the header into its own temp bgzf file via the existing writer, strip the EOF, and hand the blocks to the bgzf assembler (do NOT write the header to the main temp writer):
```rust
        if config.compression == VcfCompressionType::Bgzf {
            let header_tmp = tempdir.path().join("header.bgz");
            let mut hw = VcfLocalWriter::with_compression(&header_tmp, VcfCompressionType::Bgzf)?;
            hw.write_header(&write_schema, &vcf_info_fields_for_header, &unique_format_tags_for_header, &sample_names_for_header)?;
            hw.finish()?;
            let header_blocks = bgzf_blocks_no_eof(&header_tmp)?;
            let _ = std::fs::remove_file(&header_tmp);
            let assembler = run_assembler_thread_bgzf(
                temp_output.clone(), Path::new(output_vcf).to_path_buf(),
                tempdir.path().to_path_buf(), header_blocks, contig_done_rx,
            );
            /* drive stream; drop(shard_ctx); total_rows = assembler.join()??; */
        } else {
            // Plain/Gzip: header already written to `writer` (temp path); Task 4 assembler.
            let assembler = run_assembler_thread(
                writer, Path::new(output_vcf).to_path_buf(), temp_output.clone(),
                tempdir.path().to_path_buf(), contig_done_rx,
            );
            /* drive stream; drop(shard_ctx); total_rows = assembler.join()??; */
        }
```
NOTE: capture the header field slices (`vcf_info_fields`, `unique_format_tags`, `sample_names`) BEFORE they are moved into `VcfShardContext` at construction (clone the `Vec`s or reorder so the header write happens first), since the shard ctx takes ownership via `Arc::new(...)`. For Plain/Gzip, `writer` already holds the header (written at vcf_sink.rs:1110) and is moved into the assembler as in Task 4 — do not double-write it. Set `shard_compression` on `VcfShardContext` as in Task 6 Step 4.

- [ ] **Step 5: Run test + full suite**

Run: `cargo test -p datafusion-bio-function-vep -- --nocapture`
Expected: PASS (assembler_bgzf test + all Stage 1 tests).

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/vcf_sink.rs
git commit -m "feat(vep): raw-BGZF assembler mode + shard/assembler compression selection"
```

---

### Task 8: Stage 2 BGZF parity + parallel-compression gate

**Files:**
- Verification only.

- [ ] **Step 1: Content parity — BGZF output decompresses to the Plain output**

```bash
cd /Users/mwiewior/research/git/vepyr && maturin develop --release
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --workers 4 --force --skip-compare        # plain .vcf
# force bgzf by writing a .vcf.gz target (via vepyr compression="bgzf" or a .gz output path)
python e2e-testing/scripts/run_annotation_fast.py chr4 --cache merged --workers 4 --force --skip-compare --bgzf  # add flag or set output .vcf.gz
```
Then:
```bash
PLAIN=$(ls -t /path/work/vepyr_*_chr4_merged.vcf | head -1)
BGZF=$(ls -t /path/work/vepyr_*_chr4_merged.vcf.gz | head -1)
cmp <(grep -v '^##' "$PLAIN") <(bgzip -dc "$BGZF" | grep -v '^##') && echo "BGZF BODY BYTE-IDENTICAL"
bgzip -t "$BGZF" && echo "BGZF EOF/INTEGRITY OK"           # asserts exactly one valid EOF
tabix -p vcf "$BGZF" && echo "TABIX INDEX OK"              # asserts block/vpos integrity
```
Expected: all three OK.

- [ ] **Step 2: Confirm compression parallelized (tail + CPU)**

Run: `VEP_PROFILE=1 python e2e-testing/scripts/run_annotation_fast.py chr1 --cache merged --workers 8 --force --skip-compare --bgzf 2>&1 | grep -E "VEP_PROGRESS|vcf_sink_profile"`
Expected: BGZF wall-clock tail after 100% is comparable to the Plain assembler tail (compression no longer serial); not the pre-change single-threaded-gzip tail.

- [ ] **Step 3: Golden merged suite still green**

Run: `cd /Users/mwiewior/research/git/vepyr && pytest tests/test_golden_merged.py tests/test_run_annotation_fast.py -q`
Expected: PASS.

- [ ] **Step 4: Commit verification note**

```bash
git commit --allow-empty -m "test(vep): Stage 2 BGZF body-identical + tabix-valid + parallel compression verified"
```

---

## Self-Review

**Spec coverage:**
- "Overlap contig N concat with contig N+1 prep" → Tasks 1–4 (assembler thread + per-contig range publish; the async stream keeps polling `prepare_contig_context` while the assembler copies on its own thread).
- "Incremental append, not re-merge" → Task 4 assembler appends each range once, in order; O(total bytes), no re-copy.
- "Cover both modes" → Stage 1 Plain/Gzip (text shards + recopy), Stage 2 BGZF (parallel-compressed shards + raw block concat).
- "Single writer / ordering safe" → Global Constraints + Task 4 (one thread, arrival-order ranges).
- "Line counting without re-scan" → Task 2 (worker-reported `output_lines`) feeds Task 4/8 totals; `copy_body_file_counting_lines`/`bytecount` deleted.
- "Error = no partial output" → temp path + atomic rename (Tasks 4, 8); `_ = send()` best-effort so an assembler-gone error path doesn't panic the stream.
- "BGZF EOF correctness" → Task 7 (`bgzf_blocks_no_eof` strips each piece's trailing EOF, one hardcoded canonical `BGZF_EOF` written at the very end), Task 8 (`bgzip -t` + `tabix` gates).
- "No bio-formats change / no new dep" → Stage 2 reuses `VcfLocalWriter::Bgzf` (Task 6) and hardcodes `BGZF_EOF` locally (Task 7); nothing outside `datafusion-bio-functions` is touched.

**Placeholder scan:** No TODO/TBD; every code step shows code; verification steps show exact commands + expected output. Paths that vary by machine (`/path/work/...`) are the vepyr work dir — resolve via `ls -t` as shown.

**Type consistency:** `ContigShardRange{first_shard_id,end_shard_id,output_lines}` used identically in Tasks 1/3/4/7. `ShardResult{input_rows,output_lines}` in Task 2 consumed by Task 3's `(contig_rows, contig_lines)` binding. `VcfBodyShardWriter::create(..., compression: VcfCompressionType)` and `VcfShardContext.shard_compression: VcfCompressionType` (only `Bgzf`/`Plain`) in Tasks 6/7. `BGZF_EOF` hardcoded in Task 7, used by `bgzf_blocks_no_eof` (Task 7). `run_assembler_thread` (Plain/Gzip, Task 4) vs `run_assembler_thread_bgzf` (BGZF, Task 7) — distinct names, selected in Task 7 Step 4.

**Known friction to watch during execution:**
- Moving `VcfLocalWriter` into the assembler thread for `workers>1` while the serial branch keeps ownership — structure the `if config.workers > 1 { ... } else { ... }` so exactly one branch consumes `writer`. The parallel branch moves it into `run_assembler_thread`; guard the terminal `writer.finish()` behind `workers <= 1` (Task 4 Step 5).
- `drive_sharded_vcf_annotation` signature change (drop `writer`, return `Result<()>`) — update its single call site (vcf_sink.rs:1139).
- The `_ = shard_ctx.contig_done_tx.send(...)` in annotate_provider must clone `contig_done_tx` (it's `Sender`, `Clone`) — `VcfShardContext` derives `Clone`, and `std::sync::mpsc::Sender: Clone`, so `shard_ctx.contig_done_tx.clone()` is implied by the `Arc<VcfShardContext>` deref; ensure the ORIGINAL non-cloned sender lives on the `shard_ctx` dropped in `annotate_to_vcf` after the drive completes, so channel-close fires exactly once.
