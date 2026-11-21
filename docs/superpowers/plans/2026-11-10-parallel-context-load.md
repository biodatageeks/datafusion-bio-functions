# Parallel Context Load Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** At workers > 1, the per-chromosome context (transcripts, exons, translations, regulatory features, motifs) loads on several threads instead of one, so workers start sooner.

**Architecture:** Today `prepare_contig_data` scans and parses transcripts on the driving task, then joins the four other scans with `tokio::try_join!` on that same task, so every Parquet decode is serial. This plan decodes row groups on spawned tasks, starts the four other scans together with the transcript scan (they need the transcript ids only to filter after parsing), and parses transcript batches on blocking threads. At one worker, and with `VEP_CTX_PARALLEL=0`, the load stays serial.

**Tech Stack:** Rust 2024, tokio, parquet `ParquetRecordBatchStreamBuilder`; changes in `datafusion/bio-function-vep/src/parquet_cache/scan.rs` and `src/annotate_provider.rs`.

**Spec:** `docs/superpowers/specs/2026-11-10-parallel-vcf-output-bottlenecks-design.md` (sections 1.3, 3.3, 4)

## Global Constraints

- Branch `perf/parallel-context-load` off `origin/master`.
- The loaded context is identical: same features in the same order, same `translateable_seq` map.
- `context_parallel_enabled(annotation_workers, env)` (17248) remains the only switch: serial at one worker, `VEP_CTX_PARALLEL` overrides.
- Dropping the prepare future (an aborted prefetch, a cancelled query) stops every spawned scan and parse.
- Gate targets (spec section 4): 8-worker sum of `time_s` at least 2% lower on HG002; ClinVar target `-2` (not slower beyond noise). Additionally the summed `context done` trace time over HG002 chr1-chr22 at 8 workers is at least 35% lower.
- Measured on HG002 chr1 (merged) before the change: `context_tx_scan=0.146s context_tx=0.143s context_exons_scan=0.188s context_tl_scan=0.170s`, 0.59 s in total.
- Line numbers refer to `68d0f7b`.

## Review Focus

- A shard with one row group, and a shard that does not exist (`motif` on most contigs): same result as today (tests in Tasks 1 and 2).
- A transcript id that occurs in two batches: the later batch wins in the `translateable_seq` map, as in the serial parser (test in Task 3).
- The prefetch of the next contig runs these scans while the current contig annotates: peak RSS must stay within the gate's +10% (gate rule; ClinVar chr1 is the largest context).
- A scan error on one entity must surface as that error, not as a join error or a hang of the others (test in Task 2).
- One worker: no spawned scan, no blocking-pool parse (test in Task 2 asserts the serial arm is taken).

---

### Task 1: Row-group-parallel context scan

**Files:**
- Modify: `datafusion/bio-function-vep/src/parquet_cache/scan.rs` (`read_context_parquet` at 82; tests at 126)

**Interfaces:**
- Consumes: `existing_projection_indices` (existing).
- Produces: `pub async fn read_context_parquet_parallel(path: &Path, requested_columns: &[&str]) -> Result<Vec<RecordBatch>>` — the same batches, in file order, as `read_context_parquet`, with each row group decoded by its own `tokio::spawn`ed task. A file with at most one row group delegates to `read_context_parquet`.

- [ ] **Step 1: Write the failing tests**

```rust
#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn parallel_context_scan_matches_the_serial_scan() {
    // 10,000 rows (id: Int64 0..10000, name: Utf8 "t{id}"), max_row_group_size 2,500 -> 4 row groups.
    let path = write_fixture(10_000, 2_500);
    let serial = concat(&read_context_parquet(&path, &["id", "name"]).await.unwrap());
    let parallel = concat(&read_context_parquet_parallel(&path, &["id", "name"]).await.unwrap());
    assert_eq!(serial, parallel);
    assert_eq!(parallel.num_rows(), 10_000);
}

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn parallel_context_scan_keeps_projection_rules() {
    // A missing column is dropped and a quoted name matches, as in the serial reader.
    let path = write_fixture(100, 25);
    let out = concat(&read_context_parquet_parallel(&path, &["\"name\"", "absent"]).await.unwrap());
    assert_eq!(out.num_columns(), 1);
    assert_eq!(out.schema().field(0).name(), "name");
}

#[tokio::test]
async fn parallel_context_scan_single_row_group() {
    let path = write_fixture(10, 100);
    assert_eq!(read_context_parquet_parallel(&path, &["id"]).await.unwrap().iter().map(|b| b.num_rows()).sum::<usize>(), 10);
}
```
`write_fixture(rows, max_row_group_size) -> PathBuf` and `concat(&[RecordBatch]) -> RecordBatch` are test helpers in the same module.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib parallel_context_scan`
Expected: compile error, `read_context_parquet_parallel` not found.

- [ ] **Step 3: Implement**

Open the file once for metadata (`ArrowReaderMetadata::load_async`), compute the projection mask once, then per row group `i` spawn a task that opens its own `tokio::fs::File`, builds `ParquetRecordBatchStreamBuilder::new_with_metadata(file, metadata.clone()).with_projection(mask.clone()).with_row_groups(vec![i])` and collects its batches. Await the handles in row-group order; a failed task returns its own error; dropping the returned future aborts the tasks (hold the handles in the existing `AbortOnDrop` wrapper from `annotate_provider.rs`, moved to `crate::` scope if needed).

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib context_scan read_context`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): decode context shard row groups on separate tasks"
```

### Task 2: Start the other four scans with the transcript scan

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`load_contig_context_rest` at 10894, `scan_context_entity` at 11038, `prepare_contig_data` at 16252-16410)

**Interfaces:**
- Consumes: `read_context_parquet_parallel` (Task 1), `context_parallel_enabled`, `parse_context_batches`, `ContextParseCancel` (existing).
- Produces:
  - `enum ContextRestScans { Started(Box<[tokio::task::JoinHandle<Result<(Vec<RecordBatch>, Duration)>>; 4]>), Deferred, Skipped }` — order: exon, translation_core, regulatory, motif.
  - `fn start_context_rest_scans(cache: &PartitionedParquetCache, chrom: &str, concurrent: bool, needs_transcript_context: bool, read_protein_features: bool) -> ContextRestScans` — `Started` when `concurrent`, `Deferred` when serial, `Skipped` when no transcript context is needed.
  - `async fn finish_context_rest(scans: ContextRestScans, cache: &PartitionedParquetCache, chrom: &str, profile: &Option<SharedContigPipelineProfile>, tx_ids: &HashSet<String>, concurrent: bool, read_protein_features: bool) -> Result<ContigContextRest>` — the scans' results, parsed and filtered exactly as `load_contig_context_rest` does today. `load_contig_context_rest` is removed.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn context_rest_scans_are_deferred_when_serial_and_skipped_without_context() {
    let cache = empty_parquet_cache();          // helper: a cache root with no shards
    assert!(matches!(start_context_rest_scans(&cache, "1", false, true, true), ContextRestScans::Deferred));
    assert!(matches!(start_context_rest_scans(&cache, "1", true, false, true), ContextRestScans::Skipped));
}

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn started_and_deferred_context_rest_agree() {
    // Fixture cache with exon and translation_core shards for contig "1" (3 row groups each),
    // no regulatory or motif shard. tx_ids keeps two of three transcripts.
    let (cache, tx_ids) = context_fixture_cache();
    let started = finish_context_rest(start_context_rest_scans(&cache, "1", true, true, true),
        &cache, "1", &None, &tx_ids, true, true).await.unwrap();
    let deferred = finish_context_rest(ContextRestScans::Deferred,
        &cache, "1", &None, &tx_ids, false, true).await.unwrap();
    assert_eq!(started.0, deferred.0);   // exons
    assert_eq!(started.1, deferred.1);   // translations
    assert!(started.2.is_empty() && started.3.is_empty());
}

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn a_failing_context_scan_returns_its_own_error() {
    // The exon shard is a text file, not Parquet.
    let (cache, tx_ids) = context_fixture_cache_with_corrupt_exons();
    let err = finish_context_rest(start_context_rest_scans(&cache, "1", true, true, true),
        &cache, "1", &None, &tx_ids, true, true).await.unwrap_err();
    assert!(err.to_string().contains("context parquet"), "{err}");
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib context_rest`
Expected: compile error, `ContextRestScans` not found.

- [ ] **Step 3: Implement**

`start_context_rest_scans` resolves the four paths with `cache.context_path(entity, chrom)` and spawns one task per entity that returns `(read_context_parquet_parallel(..), elapsed)`, or an empty result for a missing shard. In `prepare_contig_data`, call it immediately before `context_fut` is first polled and keep the value in an abort-on-drop guard; after the transcripts are parsed and `tx_ids` is built, call `finish_context_rest`. The `context done` trace event and the `context_*_scan` profile fields keep their meaning (a scan's elapsed time is measured inside its task).

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib context`
Expected: PASS, including `context_parse_by_batch_matches_whole_parse_and_stops_when_cancelled`.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): overlap the non-transcript context scans with the transcript load"
```

### Task 3: Parse transcript batches on blocking threads

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (`load_contig_transcripts` at 10864, `AnnotateProvider::parse_transcript_batches` at 5201)

**Interfaces:**
- Consumes: `AnnotateProvider::parse_transcript_batches`, `ContextParseCancel`, `read_context_parquet_parallel` (Task 1).
- Produces: `async fn parse_transcript_batches_parallel(batches: Vec<RecordBatch>, tasks: usize, cancel: &ContextParseCancel) -> Result<(Vec<TranscriptFeature>, HashMap<String, String>)>` — splits `batches` into at most `tasks` consecutive chunks, parses each with `parse_transcript_batches` under `spawn_blocking`, appends the vectors in chunk order and extends the map in chunk order. `load_contig_transcripts` uses it, and `read_context_parquet_parallel` for the scan, when `context_parallel_enabled`; `tasks` is `min(4, config.annotation_workers)`.

- [ ] **Step 1: Write the failing test**

```rust
#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn parallel_transcript_parse_matches_the_serial_parse() {
    // Six batches of 50 transcripts; "ENST_DUP" appears in batches 1 and 4 with different
    // translateable sequences ("AAA" then "CCC").
    let batches = transcript_fixture_batches();
    let serial = AnnotateProvider::parse_transcript_batches("transcript", &batches).unwrap();
    for tasks in [1, 2, 4, 8] {
        let parallel = parse_transcript_batches_parallel(batches.clone(), tasks, &ContextParseCancel::default())
            .await.unwrap();
        assert_eq!(parallel.0, serial.0, "tasks={tasks}");
        assert_eq!(parallel.1, serial.1, "tasks={tasks}");
    }
    assert_eq!(serial.1.get("ENST_DUP").map(String::as_str), Some("CCC"));
}
```
`transcript_fixture_batches()` builds batches with every column `parse_transcript_batches` requires (it names each missing column in its error) plus the translateable-sequence column.

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib parallel_transcript_parse`
Expected: compile error, `parse_transcript_batches_parallel` not found.

- [ ] **Step 3: Implement** as described under Interfaces; each blocking task checks the cancel flag between batches, like `parse_context_batches`.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): parse transcript batches on blocking threads"
```

### Task 4: Smoke, trace check and the gate

**Files:**
- Modify: `datafusion/bio-function-vep/PERFORMANCE.md` (new section "Parallel context load")

**Interfaces:**
- Consumes: `build_venv.sh`, `run_arm.sh`, `gate.sh` (plan 0).

- [ ] **Step 1: Lint and unit tests**

Run: `cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all`
Expected: all pass.

- [ ] **Step 2: Smoke**

```bash
G="$HOME/research/benchmarks/local-gate"; L=datafusion/bio-function-vep/scripts/local_gate
git worktree add --detach "$G/p3/functions-baseline" origin/master
$L/build_venv.sh "$PWD" "$G/p3/venv-candidate"
$L/run_arm.sh "$G/p3/venv-candidate" "$G/p3/smoke/par" hg002 8 parity --chroms "chr1 chr22"
$L/run_arm.sh "$G/p3/venv-candidate" "$G/p3/smoke/ser" hg002 8 time --chroms "chr1 chr22" VEP_CTX_PARALLEL=0
$L/run_arm.sh "$G/p3/venv-candidate" "$G/p3/smoke/par" clinvar 8 parity --chroms "1 MT"
```
Expected: `PASS` verdicts; equal `body_md5` in `par` and `ser`.

- [ ] **Step 3: Context time from the trace**

Build `venv-baseline` (`$L/build_venv.sh "$G/p3/functions-baseline" "$G/p3/venv-baseline"`). For each venv run `run_arm.sh <venv> "$G/p3/trace/<name>" hg002 8 time VEP_PIPELINE_TRACE=1`. From each contig log take `t_ms` of the `stage=context event=done` line and sum over the 22 contigs.
Expected: candidate sum at most 65% of the baseline sum. Record both sums.

- [ ] **Step 4: Full gate**

```bash
$L/gate.sh "$G/p3/functions-baseline" "$PWD" "$G/p3/gate" --target-pct HG002=2,ClinVar=-2
cat "$G/p3/gate/verdict.txt"
```
Expected: every rule `PASS`, including `peak rss`.

- [ ] **Step 5: Document, commit and open the PR**

Section "Parallel context load": what changed, the step 3 sums, the gate totals and verdict, both build identities, and the note that a multi-contig run prefetches the next contig's context, so the gain there is on the first contig only.

```bash
git commit -am "docs(vep): parallel context load measurements"
```
