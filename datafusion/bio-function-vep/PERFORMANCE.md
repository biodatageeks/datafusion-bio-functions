# Annotation reuse and memory limits

Dense coding inputs repeatedly query the same transcript sequences and prediction
pages. The annotation engine reuses immutable reference sequences and decoded
Parquet prediction pages within an annotation execution. Neither changes the
cache format or the annotation API.

| Internal diagnostic setting | Default | Meaning |
|---|---:|---|
| `VEP_REFERENCE_CACHE_MB` | 128 MiB | Admission limit for newly owned reference CDS and peptide payloads |
| `VEP_SIFT_CACHE_MB` | 256 MiB | Admission limit for decoded prediction pages, their slots, and residency tables |

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
