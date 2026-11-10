# Parallel VCF output: four measured bottlenecks and their fixes

Date: 2026-11-10. Measured on `datafusion-bio-functions` `4c6cc45` (src identical to
`68d0f7b`, v0.24.1) through vepyr `0b02a64`, macOS, Apple M3 Max (12 P + 4 E cores),
release 116 caches, 8 workers, one process per chromosome, plain VCF output.

This document is the spec for five plans in `docs/superpowers/plans/`:

| Order | Plan | Bottleneck |
|---|---|---|
| 0 | `2026-11-10-local-gate.md` | none; the validation harness the others use |
| 1 | `2026-11-10-shard-run-pool.md` | idle workers, and most of the output copy |
| 2 | `2026-11-10-assembler-follow-head-shard.md` | what is left of the output copy (conditional) |
| 3 | `2026-11-10-parallel-context-load.md` | serial context load |
| 4 | `2026-11-10-long-cds-coding-change.md` | whole-sequence work per variant on long transcripts |

Each plan ships on its own branch and PR and must pass the gate in section 4 before the
next one starts, so every measurement has one cause.

## 1. Evidence

Sums over per-chromosome runs, seconds. "Idle workers" is annotation-phase wall minus
annotation-phase CPU / 8. "Copy tail" is the time from the last formatted batch to
process exit. The harness (thread CPU polling, `VEP_PIPELINE_TRACE=1`, shard-file size
polling) runs 5-10% slower than an untraced benchmark; use these for attribution.

| Panel | Wall | Idle workers | Copy tail | Context load |
|---|---:|---:|---:|---:|
| ClinVar, merged (25 contigs) | 272.4 | 81.2 (30%) | 39.6 (15%) | 5.9 (2%) |
| ClinVar, Ensembl (22 autosomes) | 153.0 | 32.0 (21%) | 20.5 (13%) | 3.9 (3%) |
| HG002, merged (22 contigs) | 73.8 | 9.3 (13%) | 12.2 (17%) | 5.5 (7%) |
| HG002, Ensembl (22 contigs) | 53.0 | 3.8 (7%) | 7.7 (15%) | 3.9 (7%) |

### 1.1 Idle workers

`plan_grid_partitions` / `plan_grid_partitions_balanced` give each worker one static,
contiguous slice of whole 5,000-unit buffers. A worker that finishes does nothing.

- ClinVar chr2, merged: worker 0 planned at 14.8% of the weight, ran 40.6 s (30% of all
  worker time). The TTN single-buffer slices, planned at 8% each, ran 8-9 s.
- ClinVar chr6, merged: one buffer planned at 21.7% took 3.1 s (7.5%). Buffers are
  atomic, so it set the min-max limit and three workers were packed up to it
  (7.5-13.0 s) while five finished within 3.1 s.
- ClinVar chr2, Ensembl (uniform split): seven workers 7.7-9.3 s, the TTN worker 20.1 s.
- Paired runs, 23 merged ClinVar contigs: weighted 266.1 s, uniform 292.5 s, best of the
  two per contig 252.6 s. No static policy closes the gap. Sum of longest-worker time
  212.1 s against sum of mean-worker time 143.2 s (HG002 merged: 52.3 against 41.1).
- CPU per unit of work rises when all workers are busy: ClinVar chr6 produced identical
  output with 40.6 CPU-s (lopsided split) and 51.3 CPU-s (even split). The gap figures
  are therefore upper bounds on the gain.

The streaming (non-file) path already schedules dynamically: `RunPoolState` cuts a contig
into more runs than workers (`plan_stream_runs`) and starts them as slots free up.

### 1.2 Output copy

Each grid worker writes `partition_NNNN.vcf.body` in `std::env::temp_dir()`.
`run_assembler_thread` copies them into `<output>.part` in id order, single-threaded,
and only after `ContigShardRange` arrives, which is sent when the whole contig is done.

- In every measured run the copy started within 0.1 s of the last worker finishing.
- ClinVar merged: 36.5 s for 78.0 GB (2.1 GB/s); HG002 merged: 7.0 s for 26.1 GB.
- With `compression="bgzf"` the copy is 0.14-0.21 s; ClinVar chr1 23.0 s plain against
  19.1 s BGZF, chr17 26.6 against 22.9.
- Replaying the recorded shard growth: publishing each shard as soon as it and its
  predecessors are complete removes only 23% of the tail on the static split (36.5 to
  28.0 s), because eight equal slices finish together. Following the growing head shard
  as well removes 37%. Many small in-order runs are what makes early publishing pay.

### 1.3 Context load

`prepare_contig_data` runs on the task that drives the stream. On HG002 chr1 (merged) one
thread is at 1.0 core for about 0.55 s before any worker starts:
`context_tx_scan=0.146s context_tx=0.143s context_exons_scan=0.188s
context_tl_scan=0.170s` (`VEP_PROFILE=1`). The four non-transcript scans are joined with
`tokio::try_join!` on one task, so their Parquet decode is serial, and they start only
after the transcripts are parsed. chr1 shards: transcript 2 row groups, exon 15,
translation_core 2. In a multi-contig run the next contig's load is prefetched, so this
cost shows on the first contig and in per-chromosome invocations.

### 1.4 Long transcripts

CPU per CSQ entry is 5.4-6.4 us on all four panels and 7.9 us on ClinVar chr2 (merged).
In the chr2 stack profile, coding-change classification with protein translation is 33%
of busy samples (chr17: 17%). In `classify_coding_change_with_spliced_span_cached`:
`mutated` (a full copy of the CDS) is built for every variant-transcript pair;
`substitute_peptide_window` returns a full `Vec<char>` peptide; any indel takes
`translate_protein_for_coding_classification` over the whole CDS; frameshift and
stop-lost variants then call `alternate_translation_for_vep_hgvs`, which rebuilds and
retranslates CDS plus 3' UTR. TTN (`2:178525989-178830802`) holds 41,279 ClinVar records.

## 2. Goals and non-goals

Goals: remove the idle-worker loss and the copy tail for VCF-file output at workers > 1;
shorten the serial start; cut whole-sequence work on long transcripts. Output bytes do
not change for any input, cache type or worker count.

Non-goals: cache format or Python API changes; the streaming `RecordBatch` path's
scheduling policy; the VCF scan count of the streaming run pool; the
`balance_partitions` behaviour in `datafusion-bio-formats`; BGZF assembly.

## 3. Designs

1. **Shard run pool.** VCF-file output at workers > 1 uses `RunPoolState`. Each run
   writes its own shard; shard ids follow run order; a run's shard is published to the
   assembler when the run is released at the head. Runs are cut by cost weight where
   weights exist (Merged, RefSeq), by buffer count otherwise. `VEP_SHARD_RUN_POOL`
   selects the path (`0` keeps the static grid path) so one binary gives both arms.
2. **Follow the head shard** (only if the copy tail is still at least 3% of wall after
   plan 1). The Plain/Gzip assembler also appends the bytes already written to the
   lowest unpublished shard, so the copy overlaps annotation.
3. **Parallel context load.** Start the four non-transcript scans with the transcript
   scan, decode row groups on separate tasks, parse transcript batches on blocking
   threads. At one worker everything stays serial.
4. **Long-CDS coding change.** Profile the TTN block first; then remove whole-CDS and
   whole-peptide copies from the substitution path and translate indels from the first
   affected codon.

## 4. Local gate (every plan)

Driven by `datafusion/bio-function-vep/scripts/local_gate/` (plan 0), which builds vepyr
against a functions checkout and calls vepyr's `e2e-testing/scripts/run_comparison.py`
once per chromosome. Merged cache only; data under `$DATA_VEPYR_DIR`
(`~/workspace/data_vepyr`).

| Dataset | Profile | Input | Contigs |
|---|---|---|---|
| HG002 | `merged` | `input/HG002_normalized.vcf.gz`, `--no-normalize` | chr1-chr22 |
| ClinVar | `clinvar_merged` | `input/clinvar/clinvar_20261011.norm.vcf.gz`, `--no-normalize` | 1-22, X, Y, MT |

**Parity** (candidate build, 8 workers, plain output):
`--comparison-mode md5 --md5-mode both` for HG002 and `--md5-mode strict-clin-sig` for
ClinVar. Every contig must print `PASS`; the runner must exit 0.

**Identity** (raw body MD5, `grep -v '^#' | md5`, per contig): the candidate at 1, 4 and
8 workers and the baseline at 8 workers must all agree. For plans 1 and 2 also the
candidate with BGZF output, decompressed.

**Performance** (`--skip-compare`, sum of `annotation.time_s` over contigs, order
baseline, candidate, candidate, baseline; mean of two per arm; each process under
`/usr/bin/time -l`):

| Check | Rule |
|---|---|
| 8 workers, each dataset | meets the plan's own target |
| any single contig, 8 workers | not slower by more than 10% and 0.5 s |
| peak RSS, max over contigs | at most +10% |
| 1 worker, each dataset (one run per arm) | not slower by more than 3% |

Runs need at least 20 GiB free disk and background CPU under 200%; the scripts delete
each contig's VCFs before the next. Docker is not needed (the VEP goldens are on disk).
The CI `/parity-wgs` run on the vepyr pin PR stays the release gate; the local gate uses
a path dependency and never qualifies a release.

## 5. Global constraints

- Rust 1.91.1, `RUSTFLAGS="-C target-cpu=native"` for every measured build.
- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings`
  and `cargo test --all` pass before any gate run.
- New environment switches are diagnostic, documented in `PERFORMANCE.md`, and parsed
  with `pipeline_trace::enabled_from_env_value`.
- No timed run while another benchmark or build is running.
- Results, with build identity and commands, go into `PERFORMANCE.md` in the same PR.
