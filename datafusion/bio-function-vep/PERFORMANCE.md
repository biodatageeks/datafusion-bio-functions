# Annotation reuse and memory limits

Dense coding inputs repeatedly query the same transcript sequences and prediction
pages. The annotation engine reuses immutable reference sequences and decoded
Parquet prediction pages within an annotation execution. Neither changes the
cache format or the annotation API.

| Internal diagnostic setting | Default | Meaning |
|---|---:|---|
| `VEP_REFERENCE_CACHE_MB` | 0 at one worker; 128 MiB otherwise | Admission limit for newly owned reference CDS and peptide payloads |
| `VEP_SIFT_CACHE_MB` | 128 MiB at one worker; 256 MiB otherwise | Admission limit for decoded prediction pages, their slots, and residency tables |

The limits are shared across workers and active/prefetched contigs. They are
ceilings, not eager allocations. Zero disables the corresponding reuse. Invalid
or overflowing values return an error. A disabled-cache run is an ablation of
the new implementation, not a reconstruction of the old implementation.

Reference slots distinguish the actual transcript and translation objects.
Hydrated transcripts and translation-only overrides receive local slots instead
of reusing a stable-ID match. Eligible, unambiguous substitutions translate only
their affected codons. Edited references, ambiguous/partial codons and indels
retain the existing full-translation fallback. Peptide assembly still copies a
whole peptide; this optimization does not make all coding work constant-time.

Prediction reuse retains original versioned binary blobs in Arrow batches;
scores, categories and noncanonical amino acids are not requantized. The primary
residency unit is a key-column page. Requested predictions are decoded once per
variant row and shared by CSQ and typed output. No permanent per-position memo
grows with the number of annotated variants.

Each prediction store holds at most 4,096 page slots. Its residency table has a
512 KiB conservative charge, and every slot has a 512-byte charge. Payload
leases remain charged while readers retain an evicted page. Per-page loading is
single-flight; I/O and decompression occur outside the residency lock. Missing
keys in resident pages require no further load. Read errors propagate rather
than becoming cached missing predictions.

Bulk fills use two shared 16 MiB scratch permits. Footer byte counts bound a
fill before reading payloads; oversized pages and files without suitable bounds
use the original point reader. Busy or exhausted budgets also fall back without
waiting for memory. The scratch limit applies to additional bulk cache fills,
not to every allocation made by the existing point reader.

Admission limits are not RSS limits. Fixed reference-slot metadata, transient
row results, the annotation pipeline, and allocator retention also consume
memory. Measure peak RSS in a fresh process at each worker count.

`VEP_ENGINE_PROFILE=1` emits reference reuse and prediction I/O counters. The
existing `collapse_pick_sort` timer retains its inclusive meaning for comparison
with older binaries. `sift_warm` and `collapse_pick_sort_excluding_sift_warm`
separate prediction loading from the remaining work. Global execution peak
counters repeat in each store's report and must be combined with `max`, not
summed.

Use the same native release build, input, output fields/compression, reference
and cache for each comparison. Discard a warm-up, verify record identity/order
and exact output-body digests, and compare both wall time and peak RSS. Competitor
timings require the same input and worker count; vep-rs has known differences in
regulatory and edited RefSeq support, so it is a timing reference rather than the
correctness oracle.

## Whole-buffer worker balancing

Merged and RefSeq VCF output now balances whole input buffers by transcript
coverage rather than assigning an equal number of buffers to each worker.
This preserves the existing buffer boundaries, warm-up replay and output order.
Set `VEP_GRID_BALANCE=0` to reproduce uniform scheduling for diagnostic ablation.
At one worker the scheduling policy is unchanged. Planning temporarily retains
one `(position, input_units)` pair per input row on the contig.

ClinVar chr17 (277,114 records, eight workers) measured 53.79 seconds process
wall with balancing versus 132.65 seconds for uniform partitions using the
same A+B native binary. Both matched the original output body exactly. The
vep-rs v0.5.0 reference on the same input took 32.72 seconds; its annotation
feature coverage differs, so it is a timing reference, not a parity oracle.
These are one measured run after warm-up, not whole-ClinVar results.

Serial output uses at most 1,024 input rows per Arrow batch; parallel VCF shards
use at most 640 rows. Both prepare the full VEP input-buffer context first.
Each shard formats and writes a completed batch before annotating the next one,
so a whole window's annotated output is no longer retained. The generic parallel
RecordBatch stream preserves its existing batch and asynchronous budget policy.
This does not change the 5,000-input-unit state transitions, transcript hydration,
HGNC persistence, warm-up or partition seams. Serial annotation sends each
completed output batch through its existing bounded one-batch queue. The next buffer still waits for the current buffer's completed carry
state. Dropping the consumer closes the queue and wakes a blocked producer.
One very large record remains unbounded in bytes; batch count is not a byte cap.

The prediction cache defaults are 128 MiB for serial work and 256 MiB shared
across parallel workers. A 64 MiB experiment on ClinVar chr17
caused repeated eviction and decoding (about 121 GB of decoded pages) and was
rejected despite exact output parity. Smaller payload limits do not guarantee
lower process RSS when repeated decoding increases allocator churn.

Serial annotation uses the indexed full translator without retaining reference
CDS/peptide payloads by default. That policy was selected using an earlier chr17
configuration: disabling reuse took 166.11 seconds and 4.23 GiB peak RSS, versus
146.92 seconds and 5.10 GiB with reuse. Subsequent measurements on the current
buffer implementation, below, no longer support that memory rationale. The
zero budget also disables the local-codon optimization and penalizes long
transcripts such as TTN. This remains a documented follow-up; the diagnostic
cache override has not been promoted to a new default in this change.

## Parallel output lifetime diagnostics

Projection shares Arrow arrays. The original batch must be dropped before its
projected output reaches the sink; otherwise the formatter cannot release CSQ
storage before writing owned lines or allocating a contiguous output body. Regression tests inspect
the backing buffer's ownership after the consumer drops its projected batch,
and exercise bounded chunks, error short-circuiting, full/partial/empty shard
windows, multiallelic ties, progress and both Tokio runtime flavors.

Set `VEP_OUTPUT_MEMORY=1` in a fresh process for formatter capacity diagnostics.
The log reports simultaneous CSQ, line-vector and body capacity estimates across
active workers, including a concurrent combined peak. CSQ is an upper bound if
another caller retains an Arrow alias; line capacity is conservatively charged
through conversion. Engine builders, context/cache metadata and allocator
retention are outside this scope. Do not add separate worker or category peaks
to infer a simultaneous total, or use these counters as a process RSS bound.
Diagnostics are disabled for performance acceptance runs.

The pipeline profile's `annotate` span remains inclusive of its output callback.
`output_emission` reports the callback portion separately; inline shard writing
therefore appears inside `annotate` as well as its existing formatter trace.
The established VCF trace boundaries and performance gate are unchanged.

The acceptance rule for this work is lower peak RSS with no repeatable slowdown
versus the preserved native Stage G binary, with exact output bytes and order.
The standard gate is a coarse guard; its tolerances do not authorize exchanging
annotation speed for lower memory. Whole-ClinVar validation requires separate
user approval after the chromosome-17 implementation loop.

### Chromosome-17 memory measurements (2026-11-06)

The 640-row bulk-write implementation passed the original phase/time/RSS gate
on dense ClinVar and sparse HG002 chr17, at one and eight workers, against fresh
Stage G controls. It also passed the no-slowdown application comparison. No
thresholds were relaxed.

| Engine | Median application time | Median peak RSS |
|---|---:|---:|
| Stage G baseline | 29.34 s | 6.561 GiB |
| 640-row inline bulk writes | 29.16 s | 4.159 GiB |
| Official vep-rs v0.5.0 | 33.06 s | 4.198 GiB |

Median RSS fell 36.6%, with a 0.6% lower median runtime. All three paired
candidate-minus-baseline time deltas were negative (-0.21, -0.26, -0.02 s).
Candidate times ranged from 27.92 to 29.32 s; Stage G from 28.13 to 29.42 s.
The small timing gain should be read as preserved throughput, not a claim of
statistical significance. Relative to vep-rs, median time was 11.8% lower and
RSS about 0.9% lower. These are local chr17 measurements, not a guarantee for
other inputs or tail latency.

The dense eight-worker control reduced annotation time by 2.2% and peak RSS by
38.0%. Its formatting phase increased 0.9%, within the existing 5% limit.
The one-worker control was 5.2% faster, with 4.4% higher peak RSS (within the
10% guard); no serial memory reduction is claimed. Sparse controls were 0.4%
and 2.1% faster at one/eight workers, with RSS 0.5% and 5.2% lower.

Earlier experiments are retained as rejected or incomplete evidence. Each row
has its own paired baseline and official vep-rs v0.5.0 reference. Times are
application medians including startup; output validation is outside timers.
The 2,048-row row has six measurements per arm; other rows have three.

| Experiment | Stage G | Candidate | vep-rs | Candidate peak RSS | Result |
|---|---:|---:|---:|---:|---|
| 1,024-row bulk writes | 28.36 s | 28.47 s | 33.71 s | 4.679 GiB | Slower in every pair; rejected |
| 2,048-row bulk writes | 28.055 s | 27.825 s | 32.815 s | 5.275 GiB | 19.1% RSS saving missed 20% target |
| 512-row bulk writes | 28.11 s | 26.77 s | 32.51 s | 4.140 GiB | Application passed; formatting phase +5.3% failed |
| 512-row vectored writes | 27.52 s | 28.66 s | 32.86 s | 4.143 GiB | 4.1% slower; rejected |
| 512-row fixed write buffer | 29.76 s | 30.62 s | 33.95 s | 4.142 GiB | 2.9% slower; rejected |

The writer-copy experiments are absent from the current implementation. Bulk
writes preserve throughput while bounded synchronous emission limits live
output allocations. Neither cache shrinking nor reduced annotation coverage
was used. An individual very large consequence row is not bounded in bytes.

Every measured ClinVar candidate retained body MD5
`d2f003d60f8c85a91dbc2247d21d855a` and original input identity/order across
277,114 records. The sparse control has 108,376 records and body MD5
`927aa175d903867883ca3cf5a5a80510`. Warm-ups are excluded by a fixed protocol;
all measured runs, including slow results and observed host contention, remain
recorded. vep-rs coverage differs for regulatory and edited RefSeq annotations,
so it is a timing reference rather than a correctness oracle. No whole-ClinVar
result is claimed.


### Additional chromosome comparisons and remaining limits (2026-11-06)

Official vep-rs v0.5.0 and the current implementation, with identical normalized
chromosome inputs, reference, worker counts, full annotation flags and plain VCF
output. Each row is a median of three fresh processes after discarded warm-ups;
validation is outside timing. The chr17/eight-worker row reuses the preceding
paired session. The other rows were collected subsequently on the shared host.

| Chromosome | Workers | Current time | vep-rs time | Current peak RSS | vep-rs peak RSS |
|---|---:|---:|---:|---:|---:|
| 17 | 1 | 144.69 s | 210.73 s | 3.154 GiB | 3.782 GiB |
| 17 | 8 | 29.16 s | 33.06 s | 4.159 GiB | 4.198 GiB |
| 2 | 1 | 383.39 s | 313.19 s | 3.448 GiB | 3.948 GiB |
| 2 | 8 | 138.74 s | 49.84 s | 4.883 GiB | 4.072 GiB |

All current outputs retain original body digests and input row identity/order.
Chr2 contains 401,230 records; its body MD5 is
`4df6ae8cdc714236b95943abd465752d`. Competitor outputs reproduce their own
body digests across repetitions and worker counts. Biological coverage differs,
so these timings do not establish cross-engine annotation equivalence.

Two limitations remain:

- The serial reference-cache default is too broad. On the same 1,000-record TTN
  block, three-run medians were 6.80 s before the output-memory change, 6.83 s
  after it, 4.93 s with the existing `VEP_REFERENCE_CACHE_MB=128` override, and
  4.56 s for vep-rs. The vepyr arms had effectively equal peak RSS and exact
  output. Separate profiles reduced coding time from 4.86 to 3.00 s, restored
  23,630 local-translation operations, and retained identical prediction-page
  activity. The policy therefore predates the latest memory patch and disables
  a useful optimization for TTN.
- The parallel VCF path gives each worker a fixed whole-buffer slice. Transcript
  count weighting omits coding length, and idle workers cannot take a peer's
  remaining work. Chr2's eight-worker warm-up used fewer user CPU-seconds than
  vep-rs (295.75 versus 339.14) but took much longer (133.60 versus 47.92 s),
  with roughly one core busy through much of the tail. vep-rs divides batches
  into smaller parallel tasks. Existing streaming-run controls do not change
  this VCF-shard scheduling path.

One exploratory whole-chromosome run each with the reference-cache override
produced 130.72 s / 3.006 GiB on chr17 and 290.10 s / 3.216 GiB on chr2, with
exact original body and row-order matches. These are single diagnostic runs,
not acceptance medians or a changed production default. Repeat paired controls
before choosing the follow-up cache policy. Neither these experiments nor the
parallel scheduling diagnosis changes the implementation in this PR.

The host was shared, including intermittent external CPU activity. All measured
runs were retained. A requested diagnostic paused orchestration only after a
timed child exited, and a disk guard stopped the queue before chr2/eight-worker
measurements; the pending case resumed after cache cleanup. No timed process
was suspended and no performance or disk thresholds were relaxed.

Local verification for the implementation: 1,631 Rust tests passed (3 ignored),
2,635 Python tests passed (10 skipped, 1 expected failure), and pre-commit checks
passed, including Clippy. Full-input validation after these changes remains
pending; no whole-ClinVar performance or parity claim is made.
