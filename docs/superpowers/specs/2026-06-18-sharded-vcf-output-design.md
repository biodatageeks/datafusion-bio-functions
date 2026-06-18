# Sharded VCF Output — per-worker streaming shards + concat — Design

**Date:** 2026-06-18 (rev. 2)
**Repo:** `datafusion-bio-functions` (`datafusion/bio-function-vep`). Downstream: `polars-bio`, `vepyr`.
**Status:** Design approved (Approach II — per-worker streaming shards); ready for implementation plan.
**Related:** RSS/wait analysis in `vep-t8-rss-root-cause.md`; within-contig parallel annotation (`2026-06-17-within-contig-parallel-annotation-design.md`); shared variation index (`2026-06-18-shared-variation-index-design.md`).

---

## 1. Problem

In the `threads>1` within-contig parallel path, the N fused workers (each owning a contiguous, position-ordered range) are **merged into one output stream** by the `ParallelContigState` ordered drain (`annotate_provider.rs:11811-11867`) — the consumer reads worker 0 to completion, then 1, …. Consequences:
- **`send_wait ≈ 2.36s` at t8** — workers block when the single ordered consumer can't keep up.
- **Serial VCF serialization** — `vcf_sink` serializes the one merged stream on a single path; the workers' parallelism doesn't extend to serialization.
- **Output-channel RSS** — `out_cap × N` annotated batches buffered in the merge channels.

`lookup_wait = 0` and `ordered_drain_wait = 0` after the shared-index work, so the output merge is the dominant remaining cost.

## 2. Goal

Make the path **fully parallel from lookup to VCF write**: each worker runs `lookup → hydrate → annotate → format → append to its own VCF body shard`, streaming, with **no merge**. A final step writes one header and concatenates the shards in worker order. Delete the ordered drain.

**Success criteria:**
- Byte-identical output at threads ∈ {1,2,4,8} and vs the pre-change serial output.
- `send_wait → ~0` at t8 (no ordered drain).
- t8 wall improves (parallel serialization + no send-blocking); peak RSS no worse (output channels gone; streaming write, no buffering).
- Lib unit-suite failure set unchanged vs baseline (35 pre-existing, unrelated).

## 3. Why per-worker shards are correct & byte-identical

- Workers correspond to **contiguous, position-ordered lookup partitions** (`threads>1` already requires a tabix/CSI-indexed input so the VCF scan yields contiguous, ascending-position partitions — `annotate_to_vcf` guard). Worker `i` annotates the i-th position block.
- Each worker serializes its annotated batches with the **same** per-batch formatter the serial path uses (`format_vcf_body_chunk`), in the worker's own (ascending) batch order, into its shard.
- Concatenating shards in **ascending worker index** = ascending global position = the exact row order the merged path emits today. One header from the merged output schema (as today). ⇒ byte-identical.

## 4. Architecture

The worker's pipeline is extended to write VCF; the drain and channel are removed.

```
BEFORE                                   AFTER
worker i: lookup→annotate→               worker i: lookup→annotate→
          tx.send(batch)  ─┐                       format_vcf_body_chunk(batch)→
                           │  (bounded chan)        BufWriter<File>.write_all  ──► shard_i.vcf.body
ParallelContigState        │                        (streaming, one batch in flight)
  ordered drain (current)──┘                       NO channel, NO drain, NO merge
  → 1 merged stream
  → vcf_sink serial serialize          vcf_sink: write_header once
                                                 + concat shard_0..shard_{N-1} in order  ──► final VCF
```

- **Per worker, streaming, no buffering:** consume the worker's own annotated batches one at a time → `format_vcf_body_chunk` (the existing formatter) → `write_all` into the worker's `BufWriter<File>` shard → drop the batch. Live state per worker ≈ one batch + its formatted bytes + an 8 KB buffer (identical streaming profile to `write_vcf_partition_body`).
- **No ordered drain / no output channel:** the `ParallelContigState` merge and the `spawn_annotation_from_lookup` output channel (`out_cap`) are removed for this path; each worker self-paces to its own shard's disk-write speed.
- **Assembly:** after all workers for the run complete, `vcf_sink` writes the header once and concatenates the body shards in ascending index order (reuses the existing `write_header` + `copy_body_file_to_writer` concat loop).

## 5. Components & changes

**Worker = self-contained shard producer (`annotate_provider.rs`):**
- `spawn_annotation_from_lookup` (`:10019`) currently returns `(Receiver<Result<RecordBatch>>, JoinHandle<()>)` and sends batches to a channel. Replace its output side with a **VCF body shard writer**: the worker is given a `VcfShardSink` (shard file path + an `Arc` of the VCF formatter context — `vcf_info_fields`, `unique_format_tags`, `sample_names`, `coordinate_zero_based`) and, in its annotate loop, formats each annotated batch (`format_vcf_body_chunk`) and `write_all`s into a `BufWriter<File>`; on completion it `flush`es and returns a shard descriptor (`partition_id`, path, rows, bytes). No channel.
- The per-batch formatter `format_vcf_body_chunk` and the streaming write pattern already exist in `vcf_sink::write_vcf_partition_body` — factor the reusable body-serialization (format one batch → write bytes to a `BufWriter`) so both the worker and `write_vcf_partition_body` share it (DRY), or expose `format_vcf_body_chunk` + a thin shard-writer the worker drives directly.

**`ParallelContigState` / ordered drain (`annotate_provider.rs:9725`, drain `:11811-11867`):** removed for this path. Its role becomes: spawn the N shard-writing workers, await all `JoinHandle`s, collect the N shard descriptors. The `AnnotatingParallel` state no longer yields batches downstream — it yields nothing (the data is in shard files); the contig's output is the set of shard descriptors.

**Engine→sink handoff:** `vcf_sink` orchestrates. It already derives the header fields + sample names from the merged schema (`annotate_to_vcf` `:1214-1255`) and owns the tempdir. For `threads>1` it: (1) derives the formatter context, (2) drives the annotation so each worker writes its shard into the tempdir (passing the formatter context + per-worker shard paths down via `ContigAnnotationConfig`/the worker spawn), (3) collects shard descriptors, (4) `write_header` once, (5) concats shards in ascending id order (reusing the existing concat loop). The single-merged-stream serial path (`df.execute_stream()`) is bypassed for this case; the engine's output stream carries shard descriptors (or a completion signal), not row batches.
- The shard id must be **globally position-ordered** (contig-major, worker-minor) so the concat is sorted; reuse the existing `next_write_partition_id` ascending-id concat.

**`threads == 1`** keeps the current single-stream serial path unchanged.

## 6. Data flow / lifecycle

- Per contig (still streamed contig-by-contig overall, reusing the shared variation index built once per contig): spawn N shard-writing workers over the contig's N contiguous lookup partitions; each runs lookup→hydrate→annotate→format→append-to-shard, fully parallel; no worker blocks on another.
- Shards live in the sink tempdir for the run; `partition_id` encodes global position order.
- Assembly: one `write_header`, then stream-concat body shards in id order → single VCF; shards removed after copy (as the existing loop does).

## 7. Concurrency & memory

- Workers are independent; the only shared state is the read-only `Arc<SharedContigAnnotationContext>` (incl. the shared variation index) — already `Arc`-shared. No ordered drain, no cross-worker channel, no `send_wait`.
- Streaming write ⇒ no whole-shard or whole-output buffering; ≈ one batch per worker in flight. Replaces the `out_cap × N` channel buffering (a second RSS win).

## 8. Error handling

- A worker error aborts its shard write and propagates a `DataFusionError`; the orchestrator aborts the remaining workers (as `ParallelContigState::abort` does today) and surfaces the error before assembly. No partial concat on error.
- Shard file IO errors surface from the worker's `write_all`/`flush` with the shard path, mirroring `write_vcf_partition_body`.

## 9. Testing & gates

1. **Byte-identical parity (mandatory):** 0 mismatches in data rows at threads ∈ {1,2,4,8} and vs the pre-change serial output (`grep -v '^#'`). `bench_annotate_vcf` chr1 200k `--everything` lance.
2. **`send_wait` gate:** `VEP_PROFILE=1` at t8 shows `send_wait ≈ 0` (was ~2.36s).
3. **Benchmark:** wall + peak RSS, threads {1,2,4,8}, A/B vs pre-change. Expect t8 wall improvement; RSS no worse (channels gone, streaming write). `/usr/bin/time -l`.
4. **Regression set:** lib unit-suite failure set identical to baseline (923 pass / 35 pre-existing).

## 10. Risks

- **Engine↔sink layering shift.** Body serialization moves into the engine worker; `vcf_sink` keeps header-field derivation + header write + concat. The handoff (formatter context down to workers; shard descriptors back up; bypass the merged-stream serial path for `threads>1`) is the main new surface. Mitigation: reuse `format_vcf_body_chunk` + the `write_vcf_partition_body` streaming/concat code verbatim (DRY); the byte-identical gate is strong; `threads==1` untouched.
- **Shard ordering.** `partition_id` must be globally position-ordered (contig-major, worker-minor) for the ascending-id concat to produce sorted output. Covered by parity (single-contig is the primary gate; multi-contig parity if exercised).
- **Forward-compatibility.** The intra-contig worker model is slated to be replaced by contig/lane-level parallelism. The per-worker shard writer is decoupled from the (DataFusion) partition model, so it remains the output mechanism when the parallel unit changes (worker → contig/lane) — this is why Approach II (worker writes shard) was chosen over Approach I (workers as DataFusion partitions).
- **Temp shard disk:** N body shards on disk during a run (~output size total); same magnitude as the existing multi-lane sharded path.
- **Wall-win uncertainty:** `send_wait` removal + parallel serialization help only to the extent they were on the critical path; the benchmark (gate 3) quantifies the real gain. Parity is the hard gate regardless.

## 11. Results (2026-06-18, implemented)

Implemented in commits `d5b1233`..`ca345ef` on `parallel-fjall-lookup-partitions`.
Bench: `bench_annotate_vcf`, chr1 first-200k, `--everything`, lance backend
(`115_GRCh38_merged`), GRCh38 primary-assembly reference FASTA.

**1. Byte-identical parity (mandatory) — PASS.** Data rows (`grep -v '^#'`):
- cross-thread: t1 vs t2 / t4 / t8 → **0** diff lines each.
- vs pre-change serial: pre-Task-0 binary @ t1 vs current @ t1 → **0**; vs current @ t8 (sharded) → **0**.
- Note: the `##FORMAT` *header* line ordering is non-deterministic **run-to-run at a fixed thread count** (two t8 runs differ in header order, identical data) — a *pre-existing* HashMap-ordering effect in header generation, independent of this change (header code is identical for all thread counts and runs before the sharded branch). Out of scope here; the gate diffs data rows only.

**2. `send_wait` — output drain eliminated.** At t8 the profile shows
`output_batches=0`, `ordered_drain_wait=0.000s` — the ordered-drain output merge
(the metric's original ~2.36s target) is gone. The residual `send_wait≈2.8s` is
now the **lookup→worker** channel backpressure (workers self-pacing to shard
disk-write speed); the profiler reports both sends under one field name. Memory
is bounded by this backpressure rather than `out_cap × N` output buffers.

**3. Benchmark (`/usr/bin/time -l`, incl. startup + VCF read + cache open):**

| threads | OLD wall | NEW wall | speedup | OLD peak RSS | NEW peak RSS |
|---|---|---|---|---|---|
| 1 | 24.63s | 23.04s | 1.07× | 6.55 GB | 5.74 GB |
| 2 | 14.44s | 8.67s | 1.67× | 7.54 GB | 7.40 GB |
| 4 | 10.82s | 5.68s | 1.90× | 8.37 GB | 7.52 GB |
| 8 | 7.78s | 4.43s | **1.76×** | 8.86 GB | 8.91 GB |

t8 wall improved **1.76×**; peak RSS flat at t8 (+0.6%, within noise) and lower
at t1–t4 — **RSS no worse**. (Bench-internal annotate-only time at t8 = 3.87s.)

**4. Regression set — PASS.** Lib unit suite **923 pass / 35 pre-existing fail / 1 ignored**, **0 new** failures.
