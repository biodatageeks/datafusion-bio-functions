# Problem: eliminate per-buffer transcript clone to improve within-contig parallel scaling

**Date:** 2026-06-18
**Repo:** `datafusion-bio-functions` (engine). Downstream: `/Users/mwiewior/research/git/vepyr`.
**Status:** Diagnosis complete & validated; fix scoped but NOT started (deferred for a clean budget).

---

## 1. Background / what already exists (all committed, byte-identical at every thread count)

Within-contig parallel VEP annotation via a single `threads` knob is implemented and correct:

- `threads>1` → N **shared-nothing per-partition pipelines**: each lookup partition (contiguous, position-ordered, from a tabix-indexed input) feeds an independent worker with a **local** colocated map; annotation runs inline via `block_in_place`; output drained in partition order. Central poll-loop / cumulative colocated map / `AwaitingWindow` were deleted.
- Key code: `spawn_annotation_from_lookup`, `StreamState::AnnotatingParallel`, `ParallelContigState` in `datafusion/bio-function-vep/src/annotate_provider.rs`.
- Index requirement guard in `vcf_sink::annotate_to_vcf` (threads>1 needs `.tbi`/`.csi`).
- HGNC-donor partition-boundary bug fixed: donor map built from full active cache-region transcripts (`collect_hgnc_donors`/`apply_hgnc_donors`).
- **mimalloc** is the global allocator in the bench (`examples/bench_annotate_vcf.rs`) and in vepyr's cdylib — it was the dominant cap (1.67× even single-threaded). REQUIRED to reproduce numbers.
- Relevant commits: `ee2d395` (per-partition + HGNC fix), `5ba47c0` (mimalloc bench). Design + findings: `docs/superpowers/specs/2026-06-17-within-contig-parallel-annotation-design.md` §0. Memory: `within-contig-parallel-annotation.md`.

**Headline result (chr1 200k, --everything, indexed lance cache, mimalloc):**
`t1 16.0s → t8 4.94s` (within-mimalloc 3.25×); **5.4× vs the original system-malloc serial baseline (26.8s)**.

## 2. The problem to solve

Scaling plateaus at ~3.25× by t8. Root, validated by profiling (`VEP_PROFILE`, `VEP_ENGINE_PROFILE`):

- **No serial-wait bottleneck**: `lookup_wait=0`, `ordered_drain_wait=0`, `send_wait≈0.8s`; only serial section is `context_load≈0.7s`.
- The cap is **engine CPU SUM rising +37% at t8 for identical work** (7.57s→10.40s; same variant/transcript op count). Decomposed:
  - **~1/3 = boost-vs-all-core CPU clock scaling** (t1 single-core ~4GHz, t8 all-core ~3.5GHz on this 12P+4E Apple Silicon). Wall-time metric surfaces it as +CPU. **Unfixable in software.**
  - **~2/3 = shared L2/SLC cache contention** from the shared-nothing design's **~8× duplicated per-worker working set** (cloned `buffer_transcripts` + COITree indexes + local colocated map + output buffers).
- Inflation is **uniform 1.20–1.28× across ALL engine sub-stages** (`tx_query_total`, `transcript_output_materialize`, `transcript_overlap_eval`, `transcript_hgvsc`) → cache-contention signature.
- **Ruled out by experiment**: NOT E-cores (user-interactive QoS pinning `pthread_set_qos_class_self_np(0x21)` → zero change), NOT locks (uniform, not one stage spiking; geometry `OnceLock` is per-worker because clones reset it; SIFT reader `sift_*≈0`), NOT the output funnel (out_cap=16384 → no change), NOT Lance decode (`lookup_wait=0`), NOT tokio worker count (`TOKIO_WORKER_THREADS` → ~no change).

### 2a. VERIFIED evidence: per-worker mutable-state duplication + RSS growth (2026-06-18)

Code claim CONFIRMED: each per-partition `AnnotationWorkerState` (one per worker, `AnnotationWorkerState::new`) owns and duplicates: `transcript_overrides: HashMap<String, TranscriptPartitionState>` — `TranscriptPartitionState` (`annotate_provider.rs:10270`) holds owned sequence Strings `spliced_seq`/`five_prime_utr_seq`/`three_prime_utr_seq`/`translateable_seq`/`cdna_seq`/`refseq_edits`; plus `translation_overrides`, `persisted_buffer_transcripts`, `colocated_map`, `sift_cache`, `hgvs_reader` (FASTA), `window_buffer`.

MEASURED peak RSS (chr1 200k, --everything, mimalloc, `/usr/bin/time -l`):
- **t1 = 5.92 GB → t8 = 25.7 GB = ~4.3×** (sub-linear, NOT ~8×, but large). Corroborates the cache/memory-bandwidth mechanism behind the +37% engine-CPU-sum.

CORRECTIONS / caveats to the original hypothesis:
- "RSS scales **almost linearly**" overstates it — measured **~4.3× at 8 workers**, not 8×.
- The **sequence-String** fields (`spliced_seq` etc.) appear **inactive in this workload** (`hydrate=0`, profiler "0 CDS transcripts"), so `transcript_overrides` is largely empty here. The RSS growth is dominated instead by **cloned `buffer_transcripts` + COITree indexes + local colocated maps + the annotated-output channel buffer** (`PARALLEL_ANNOTATE_OUT_CAP=1024` × N partitions of large `--everything` batches). The hydration-string duplication would matter more on a hydration-heavy workload (e.g. HGVS-heavy / RefSeq).
- Therefore the fix below (immutable shared transcripts) addresses the **transcript-clone** portion; the **output-buffer** portion is a separate axis (tune `PARALLEL_ANNOTATE_OUT_CAP` / temp-file output) if RSS becomes the constraint.

## 3. The fix (the recoverable 2/3)

**Eliminate the per-buffer transcript clone → immutable shared `Arc` transcripts + HGNC side-table.**

Why it clones today: `select_buffer_local_transcripts` does `.cloned()` because `build_stateful_buffer_local_transcripts` **mutates** `gene_hgnc_id` (donation) on the buffer copies, and the geometry `OnceLock` is per-instance (clone resets it).

Change:
- `select_buffer_local_transcripts` / `build_*_buffer_local_transcripts`: return **`&TranscriptFeature` refs into `shared.base_transcripts`** instead of owned clones.
- Move the donated HGNC result into a **side-table keyed by `transcript_id`** (e.g. `HashMap<&str, String>`), produced by `apply_hgnc_donors` over region-complete donors.
- The engine's per-transcript HGNC read (`tx.gene_hgnc_id` inside `evaluate_variant_prepared`, `transcript_consequence.rs`) consults the side-table instead of the field.
- Keep the geometry `OnceLock` on the **shared `Arc` transcripts** → built once globally, lock-free reads (we established this is NOT contended).
- Per-window `transcript_overrides`/`translation_overrides` similarly via side-tables.

Effect: 8 workers share one cache-resident transcript set (not 8× copies) → less cache pressure (scaling) AND less allocation (serial speed). Attacks both halves of the absolute engine cost.

**Touch points:** `annotate_provider.rs` — `select_buffer_local_transcripts` (~10719), `build_buffer_local_transcripts` (~10749), `build_stateful_buffer_local_transcripts` (~10770), `apply_buffer_local_hgnc_propagation`/`collect_hgnc_donors`/`apply_hgnc_donors`, `prepared_context_from_buffer` (~10435), `PreparedContext`. `transcript_consequence.rs` — `evaluate_variant_prepared` transcript/HGNC access.

## 4. Verification gate (MANDATORY — high regression risk)

This is the same code area as the partition-boundary HGNC bug we just fixed. Output MUST stay byte-identical.

Reproduce the indexed input (if `/tmp/chr1_200k.vcf.gz` is gone):
```
bcftools view /Users/mwiewior/workspace/data_vepyr/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr1 \
 | awk 'BEGIN{c=0} /^#/{print;next} {if(c<200000){print;c++}else{exit}}' \
 | bgzip > /tmp/chr1_200k.vcf.gz && tabix -p vcf -f /tmp/chr1_200k.vcf.gz
```
Build + parity + scaling:
```
cargo build --release --features lance-cache --example bench_annotate_vcf
# run threads=1,2,4,8 to /tmp/c{T}.vcf:
./target/release/examples/bench_annotate_vcf --input /tmp/chr1_200k.vcf.gz \
  --cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --output /tmp/cT.vcf --backend lance --everything \
  --reference-fasta /Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --threads T --no-progress
# GATE: diff <(grep -v '^#' /tmp/c1.vcf) <(grep -v '^#' /tmp/c8.vcf) must be EMPTY at every T,
#       AND f1 must equal the pre-change serial output (no serial regression).
```
Profiling: `VEP_PROFILE=1` (pipeline stages incl. engine sum), `VEP_ENGINE_PROFILE=1` (sub-stage breakdown). Confirm the engine-CPU-sum +37% shrinks toward the ~1/3 clock floor.

## 5. Success criteria

- Byte-identical output at threads ∈ {1,2,4,8} vs the current serial path.
- `t1` no slower (ideally faster from less allocation).
- `t8` engine-CPU-sum inflation drops from +37% toward ~+12% (the clock-only floor); wall-clock t8 improves below 4.94s.

## 6. Out of scope / notes

- The ~1/3 clock-scaling component is physics, not fixable.
- Whole-genome (multi-contig) is the easier scaling case (more independent work, fixed context-load amortized) — the single-contig 200k bench is the hardest case.
- Production allocator win (mimalloc `#[global_allocator]`) already wired in vepyr cdylib; a library crate cannot set it, the final binary must.
