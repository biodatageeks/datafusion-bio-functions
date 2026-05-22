# VEP-Style Annotation Parallelism Implementation Plan

> **For implementation:** run this in one agent/session. Do not resurrect the
> removed partition-local annotation worker path from
> `2026-05-21-parallel-annotation-partitions.md`.

## Goal

Add real parallelism to `annotate_vep()` annotation while preserving VEP-like
input-buffer semantics and strict VCF output order.

The pipeline should be:

```text
indexed VCF/DataFusion partitions
  -> parallel fjall lookup partitions
  -> strict ordered drain by source partition id
  -> global input-unit chunker, VEP buffer_size units
  -> bounded annotation jobs
  -> ordered drain by buffer id and subchunk id
  -> VCF writer
```

## Why This Plan

The removed path ran one `AnnotationWorkerState` per lookup partition. That
looked parallel, but it made each DataFusion partition define its own VEP input
buffers. VEP does not do that. VEP builds one ordered parent `InputBuffer`, then
forks chunks from that buffer.

The current surviving path already has the right semantic boundary:
`InputBufferAccumulator` creates global ordered buffers by ALT/input-unit count.
The missing performance is that the stream tends to stop after the first ready
global buffer, so `target_partitions=2` often gives only one annotation job.

## Current Cleanup

- [x] Removed the disabled partition-local annotation branch from
  `datafusion/bio-function-vep/src/annotate_provider.rs`.
- [x] Removed `should_parallelize_annotation()`,
  `AnnotationPartitionHandle`, `ParallelAnnotationDrainState`,
  `run_annotation_partition_worker()`, `spawn_annotation_partition_worker()`,
  and `poll_next_annotated_partition_batch()`.
- [x] Simplified contig setup to always enter the global input-buffer
  annotation state.
- [x] Removed the unit test that only asserted the dead path stayed disabled.

Verified after cleanup:

```bash
cargo fmt --check
cargo check -p datafusion-bio-function-vep --features kv-cache
cargo test -p datafusion-bio-function-vep --features kv-cache test_fjall_annotation_parallelizes_global_input_buffers -- --nocapture
cargo test -p datafusion-bio-function-vep --features kv-cache input_buffer_accumulator -- --nocapture
```

## Design Constraints

- Keep `ContigAnnotationExec` as one output partition. The VCF writer remains a
  single ordered sink.
- Source partitions may be read/probed in parallel, but the coordinator drains
  lookup output strictly by source partition id.
- VEP `buffer_size` is counted in parsed input units, not Arrow rows. A
  multi-ALT row must stay intact even if it crosses an exact chunk boundary.
- The global buffer boundary is semantic. Parallel workers may split work
  inside a global buffer only after buffer-local transcripts/exons/translations
  are materialized.
- Bounded memory: no full chromosome buffering. At most roughly
  `target_partitions * input_buffer_size` input units should be staged for
  annotation, plus bounded worker outputs.

## Shared vs Local State

Safe shared immutable state:

- `SharedContigAnnotationContext`
- `AnnotateProvider`
- `TranscriptConsequenceEngine`
- base transcripts and translations
- exons, regulatory features, motifs, miRNAs, structural features
- transcript cache-region map
- translateable sequence map
- fjall stores and SIFT kv store handles
- per-buffer `colocated_map` snapshot

Keep local per annotation job or subchunk:

- `PreparedContext` construction from buffer-local features
- `SiftPolyphenCache`
- loaded SIFT window set
- FASTA indexed reader for HGVS
- output `RecordBatch` queue

Keep local to the ordered parent coordinator:

- `InputBufferAccumulator`
- `persisted_buffer_transcripts`
- transcript and translation overrides
- colocated map merge from ordered lookup deltas
- `next_input_buffer_id`

## Task 1: Make Read-Ahead Produce Enough Global Buffers

**Files:**

- `datafusion/bio-function-vep/src/annotate_provider.rs`

**Status:** Implemented.

**Change:**

Today `AnnotatingContig` stops polling lookup as soon as
`InputBufferAccumulator::has_ready_input_buffer()` is true. That commonly means
one 5000-unit buffer is annotated at a time, so `target_partitions=2` has no
work to split across whole-buffer jobs.

Add an accumulator helper:

```rust
fn ready_input_buffer_count(&self, input_unit_limit: usize) -> usize
```

Then change the lookup polling condition to continue draining ordered lookup
partitions until one of these is true:

- lookup is done
- LIMIT pushdown has enough buffered rows
- `ready_input_buffer_count(input_buffer_size) >= target_partitions`
- `window_buffer.len() >= HYDRATION_WINDOW_SIZE`

This creates a bounded batch of global VEP buffers, not a chromosome buffer.
If lookup is pending before `target_partitions` complete buffers are staged,
the stream now starts annotation as soon as at least one complete global buffer
is available, so the subchunk scheduler can still use parallel workers.

**Tests:**

- Unit test `InputBufferAccumulator::ready_input_buffer_count()` for exact,
  partial, and multi-ALT cases.
- Existing global-id test should still prove buffer ids are assigned in order.

## Task 2: Fix Whole-Buffer Job Scheduling

**Files:**

- `datafusion/bio-function-vep/src/annotate_provider.rs`

**Status:** Implemented with `tokio::task::JoinSet`.

**Change:**

Keep `prepare_input_buffer_annotation_jobs()` as the semantic producer, but
make `annotate_window_parallel_input_buffers()` collect completed jobs without
blocking on the oldest handle first.

Use `FuturesUnordered` or `JoinSet` plus an ordered ready map:

```text
spawn buffer jobs up to target_partitions
as each job completes:
  store by buffer_id
  spawn the next pending job if any
  emit only when next_output_buffer_id is ready
```

This removes head-of-line wait when a later buffer finishes early, while still
emitting strictly by `buffer_id`.

**Tests:**

- Unit-level scheduler test with artificial delayed jobs:
  - buffer 1 completes before buffer 0
  - output remains 0 then 1
  - no more than `target_partitions` jobs are in flight

## Task 3: Add VEP-Fork-Style Subchunks Inside One Global Buffer

**Files:**

- `datafusion/bio-function-vep/src/annotate_provider.rs`

**Status:** Implemented.

**New structs:**

```rust
struct AnnotationSubchunkJob {
    buffer_id: usize,
    chunk_id: usize,
    batches: Vec<RecordBatch>,
    transcripts: Arc<Vec<TranscriptFeature>>,
    exons: Arc<Vec<ExonFeature>>,
    translations: Arc<Vec<TranslationFeature>>,
    colocated_map: Arc<HashMap<ColocatedKey, ColocatedData>>,
}

struct AnnotatedSubchunk {
    buffer_id: usize,
    chunk_id: usize,
    batches: VecDeque<RecordBatch>,
}
```

**Change:**

After `prepare_input_buffer_annotation_jobs()` creates a global 5000-unit job,
split `job.batches` into bounded subchunks by input-unit count.
The implementation distributes at most `target_partitions` subchunks across
the prepared global buffers for the current window. A single ready global
buffer can therefore use all workers, while multiple ready global buffers do
not multiply into `target_partitions * target_partitions` subjobs.

Rules:

- Preserve row order.
- Do not split a VCF row across subchunks.
- Keep empty subchunks out.
- Each subchunk owns only Arrow slices, not copied rows.
- The buffer-local transcripts/exons/translations are shared via `Arc<Vec<_>>`.

Each subchunk builds its own `PreparedContext` from the shared per-buffer
feature vectors, then annotates its own batch slices. This duplicates
`PreparedContext::new()` per subchunk initially, which is acceptable for a
first implementation because it avoids unsafe borrowed context sharing across
spawned tasks. Profile it before optimizing.

**Tests:**

- Splitter test with rows whose ALT unit counts are `[1, 2, 1, 4, 1]`.
- Assert no row is split and concatenating subchunks equals original row order.
- Assert `target_partitions=1` returns one subchunk identical to the input.

## Task 4: Ordered Drain Across Buffers And Subchunks

**Files:**

- `datafusion/bio-function-vep/src/annotate_provider.rs`

**Status:** Implemented.

**Change:**

Unify whole-buffer and subchunk outputs under ordered keys:

```text
(buffer_id, chunk_id)
```

Drain policy:

- next expected buffer starts at the first prepared buffer id
- next expected chunk starts at `0`
- only append batches to output when the exact next key is ready
- advance to the next buffer after all chunks for the current buffer emitted

Use a bounded ready map:

```rust
BTreeMap<(usize, usize), AnnotatedSubchunk>
```

The map can never grow beyond the number of in-flight jobs if the scheduler
keeps `target_partitions` bounded.

**Tests:**

- Completion order `(0,1)`, `(1,0)`, `(0,0)` emits `(0,0)`, `(0,1)`, `(1,0)`.
- Error from any subchunk aborts remaining jobs and enters cleanup.

## Task 5: Keep Lookup, Chunking, And Annotation Overlapped

**Files:**

- `datafusion/bio-function-vep/src/annotate_provider.rs`

**Status:** Follow-up. Tasks 1-4 are complete; this remains the larger
coordinator rewrite.

**Change:**

The first implementation can keep the existing window future boundary. After
Tasks 1-4 pass, remove the remaining large stall by replacing the
`AnnotatingWindow` wave with a coordinator state that owns:

- ordered lookup partition handles
- `InputBufferAccumulator`
- pending annotation jobs
- in-flight annotation tasks
- ordered ready output map

The coordinator should:

1. Poll lookup until it has enough ready buffers or workers are saturated.
2. Spawn annotation subchunk jobs while capacity exists.
3. Poll completed annotation jobs.
4. Emit ordered batches when available.
5. Continue lookup while previous annotation jobs are running.

This is the largest state-machine change, so do it after the smaller scheduling
and subchunk tests are green.

## Task 6: Profiling And Metrics

**Files:**

- `datafusion/bio-function-vep/src/annotate_provider.rs`
- `datafusion/bio-function-vep/examples/bench_annotate_debug.rs`

**Status:** Follow-up.

Add profile counters:

- `ready_input_buffers`
- `annotation_subchunks`
- `annotation_inflight_max`
- `annotation_spawn_wait`
- `annotation_join_wait`
- `ordered_subchunk_drain_wait`

Benchmark commands should avoid `bench_annotate_vcf --limit` if that creates an
unindexed temporary VCF. Use the original indexed `.vcf.gz` and SQL LIMIT, or
run full input:

```bash
VEP_PROFILE=1 cargo run --release --features kv-cache \
  -p datafusion-bio-function-vep --example bench_annotate_debug -- \
  --cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --vcf /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --target-partitions 1

VEP_PROFILE=1 cargo run --release --features kv-cache \
  -p datafusion-bio-function-vep --example bench_annotate_debug -- \
  --cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --vcf /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --target-partitions 2

VEP_PROFILE=1 cargo run --release --features kv-cache \
  -p datafusion-bio-function-vep --example bench_annotate_debug -- \
  --cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --vcf /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --target-partitions 4
```

## Task 7: Correctness Verification

**Status:** Implemented verification for the completed code changes.

Run:

```bash
cargo fmt --check
cargo check -p datafusion-bio-function-vep --features kv-cache
cargo test -p datafusion-bio-function-vep --features kv-cache input_buffer -- --nocapture
cargo test -p datafusion-bio-function-vep --features kv-cache vcf_roundtrip_golden_fjall -- --nocapture
```

Actual verification run:

```bash
cargo fmt --check
cargo check -p datafusion-bio-function-vep --features kv-cache
cargo test -p datafusion-bio-function-vep --features kv-cache input_buffer -- --nocapture
cargo test -p datafusion-bio-function-vep --features kv-cache split_batches_by_input_units -- --nocapture
cargo test -p datafusion-bio-function-vep --features kv-cache ordered_subchunk_drain -- --nocapture
cargo test -p datafusion-bio-function-vep --features kv-cache test_fjall_parallel_lookup_and_annotation_output_is_invariant_across_target_partitions -- --nocapture
cargo test -p datafusion-bio-function-vep --features kv-cache test_roundtrip_golden_fjall_all_column_values -- --nocapture
```

Then run an output equivalence check:

```text
target_partitions=1 output == target_partitions=2 output
target_partitions=1 output == target_partitions=4 output
```

For any mismatch, first inspect:

- HGNC fields near transcript cache-region boundaries
- multi-ALT variants crossing chunk boundaries
- `check_existing` colocated fields
- HGVS/SIFT fields when `--everything` is enabled

## Expected Impact

Task 1 should fix the immediate "only one CPU busy" symptom when there is
enough ordered lookup output to stage multiple global VEP buffers.

Tasks 3-4 should improve cases where only one global buffer is ready but it has
enough rows to split across annotation workers, matching VEP's fork model more
closely.

Task 5 should reduce idle time between lookup, annotation, and writing. This is
where the full pipeline becomes streaming parallel rather than wave-based
parallel.
