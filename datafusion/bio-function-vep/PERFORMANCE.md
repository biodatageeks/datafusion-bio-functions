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

At one worker, output is subdivided into at most 1,024 input rows per Arrow
batch after preparing the full VEP input-buffer context. Parallel annotation
retains the original output batches. This does not change the 5,000-input-unit
state transitions, transcript hydration, HGNC persistence, warm-up or partition
seams. Serial annotation sends each completed output batch through a bounded
one-batch queue instead of retaining every output batch for the whole input
buffer. The next buffer still waits for the current buffer's completed carry
state. Dropping the consumer closes the queue and wakes a blocked producer.
One very large record remains unbounded in bytes; batch count is not a byte cap.

The prediction cache defaults are 128 MiB for serial work and 256 MiB shared
across parallel workers. A 64 MiB experiment on ClinVar chr17
caused repeated eviction and decoding (about 121 GB of decoded pages) and was
rejected despite exact output parity. Smaller payload limits do not guarantee
lower process RSS when repeated decoding increases allocator churn.

Serial annotation uses the indexed full translator without retaining reference
CDS/peptide payloads by default. A chr17 diagnostic with reference reuse disabled
took 166.11 seconds and 4.23 GiB peak RSS, versus 146.92 seconds and 5.10 GiB with
reuse; both matched the original output exactly. The fresh original-binary
serial baseline took 793.62 seconds and 4.37 GiB. These measurements motivate
the default; they do not establish memory bounds or whole-ClinVar performance.
