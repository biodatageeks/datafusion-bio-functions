---
name: vep-perf-profiling
description: >-
  Profile the VEP annotation pipeline end-to-end and present a structured timing
  breakdown (pipeline stages, Lance variation/SIFT lookup, engine + transcript-engine
  stages) to find bottlenecks. Use this whenever the user wants to profile, benchmark,
  or get a performance/timing breakdown of vepyr annotation, the variation cache lookup,
  the VEP engine, or asks "where is the time going / what's the bottleneck" for chr-scale
  annotation — even if they don't say the word "profile". Also use after a cache-layout or
  encoding change to confirm lookup cost and that the pipeline is still engine-bound.
---

# VEP performance profiling & breakdown

Real end-to-end profiling runs through the sibling **vepyr** repo against a chr-scale cache,
not this repo's unit tests. The pipeline emits nested per-batch timers to stderr; this skill
runs it single-threaded (so timings are attributable, not smeared across worker threads) and
turns the raw stderr into a clean breakdown.

Background: see the `vep-profiling-setup` memory and [[vep-annotation-bottlenecks]].

## 1. Run it (single-thread, all profile levels)

Pin threads to 1 so each stage's time is real, not overlapped. Enable every profile level.
`--skip-compare` keeps it focused on timing (drop it only when you also want a parity check).

```bash
cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts
LANCE_CPU_THREADS=1 LANCE_IO_THREADS=1 RAYON_NUM_THREADS=1 \
VEP_PROFILE=1 VEP_ENGINE_PROFILE=1 VEP_TX_ENGINE_PROFILE=1 VEP_ENGINE_PROFILE_HGVSC=1 VEP_LANCE_PROFILE=1 \
  env -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py \
  chr1 --cache merged --forks 0 --force --backend lance --skip-compare \
  >/tmp/prof.stdout 2>/tmp/prof.profile
```

Profile knobs (each prints to stderr, aggregated per batch):
- `VEP_PROFILE` → `pipeline_profile` (annotate/engine/lookup_wait/hydrate/context_*) + vcf_sink + concurrency_plan.
- `VEP_LANCE_PROFILE` → `[vep-lance-profile]` per-batch: `variation_take`, `sift_key_take`, dataset `open`, `variation_resolve`.
- `VEP_ENGINE_PROFILE` → per-batch engine stages (`evaluate_prepared`, `colocated_fields`, `csq_format`, `collapse_pick_sort`, …).
- `VEP_TX_ENGINE_PROFILE` → transcript-engine stages (`tx_query_total`, `transcript_output_materialize`, `transcript_hgvsc`, `transcript_overlap_eval`, …).
- `VEP_ENGINE_PROFILE_HGVSC` → extra HGVSc sub-timers (coord_map/simple_fast/fallback); 0 if off.

Note: the tqdm `Done: N variants` line is a display artifact under heavy stderr profiling — trust
`output_rows=` in `pipeline_profile` for the real count.

## 2. Extract the breakdown

Run the script — it parses `/tmp/prof.profile` into the four sections.

```bash
bash <skill-dir>/scripts/breakdown.sh /tmp/prof.profile
```

It emits: pipeline stage line, per-event Lance lookup totals + `variation_take` per-batch
distribution, the full engine stage list, and the full transcript-engine stage list (each summed
across all batches, sorted by time).

## 3. Present it

Use this structure (scale the tables to what's interesting). The goal is to show **where the wall
time goes and whether the pipeline is engine-bound or lookup-bound** — that's the decision the
breakdown exists to support.

```
## Pipeline (VEP_PROFILE)
annotate / engine / lookup_wait / hydrate / context_load

## Lance lookup (VEP_LANCE_PROFILE)
table: event | n | total | (variation_take per-batch min/median/mean/max)

## Engine — top level (VEP_ENGINE_PROFILE), additive
evaluate_prepared (the transcript-engine call), colocated_fields, csq_format, collapse_pick_sort, ...

## Transcript engine — inside evaluate_prepared (VEP_TX_ENGINE_PROFILE), nested
tx_query_total -> transcript_output_materialize / transcript_hgvsc (coord_map/simple_fast/fallback)
                  / transcript_overlap_eval (splice_terms) / position_fields / ...

## Reading it
- engine-bound? (engine / wall). Lookup the bottleneck? (lookup_wait vs engine).
- name the top 3-4 stages to target if optimizing.
```

### Key interpretation rules (so the numbers aren't just dumped)
- **`evaluate_prepared` ≈ `tx_query_total`** — the top-level engine umbrella *is* the transcript
  engine; don't double-count them as separate costs. Engine top-level stages (`evaluate_prepared`,
  `colocated_fields`, `csq_format`, `collapse_pick_sort`, `hgvs_shift`) are additive; the tx-engine
  stages are *nested inside* `evaluate_prepared`.
- **Lookup is usually not the bottleneck.** `variation_take` may be several seconds of CPU, but it
  overlaps engine work on tokio workers — the real cost the engine pays is `lookup_wait` (typically
  ~1 s). If `lookup_wait` ≪ `engine`, the pipeline is engine-bound and lookup/encoding changes won't
  move the end-to-end number. Say so explicitly.
- This is the canonical profile: **engine-bound (~60% of wall)**, dominated by
  `transcript_output_materialize`, `transcript_hgvsc`, `colocated_fields`, `csq_format`.

## When comparing two cache layouts / encodings
Run the profile on each and compare `variation_take` (total + per-batch), `lookup_wait`, `engine`,
and AF/dataset on-disk size. For raw read-amplification (bytes read, IOPS) — which `VEP_LANCE_PROFILE`
does NOT give — use the Lance sandbox `take-batches-bench` against captured row_ids instead; that's a
separate read-amp measurement, not part of this timing skill.
