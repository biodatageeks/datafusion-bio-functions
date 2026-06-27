# Handoff: attack VEP engine redundancy (per-transcript geometry rebuild)

**Date:** 2026-06-26 · **Updated:** 2026-06-27 (measured — see banner)
**Branch:** `hgnc-grid-aligned-buffers-prototype` (or a fresh branch off it)
**Status:** grid-aligned parallel annotation DONE (correct + 2.4× w4); levers below
RE-MEASURED 2026-06-27 — the original estimates were ~3× too optimistic.

---

## ⚠️ UPDATE 2026-06-27 — levers measured against code; estimates corrected

The attack order and impact estimates in the original sections below were
ground-truthed against the actual code and **measured** on complete chr4 runs.
Memory: `vep-engine-redundancy-levers-measured`. Net: realistic eliminable work
is **~1.0–1.5s off ~11.4s engine (9–13%)**, NOT the 3–5s implied below.

**Attribution traps the original profile table hides:**
- `transcript_output_materialize` (3.44s) is a *superset* of `transcript_hgvsc`
  (2.19s) — **nested, not additive**. Non-HGVSc residual ≈ 1.2s. Don't sum them.
- `csq_format` (1.91s) *includes* the nested `sift_polyphen`/`domains`/`mirna`
  lookups; true formatting cost is far less. `collapse_pick_sort` (0.94s) is
  mostly a mislabeled SIFT prefetch — the integer-key sort itself is cheap.

**Per-lever verdicts (replacing the "Attack order" section):**
- **C (geometry rebuild): DISPROVEN as a big lever.** `GeometryCache` already
  memoizes the cDNA coord table. (1) Routing the sorted-exon/mapper fast paths
  through the cache = **0 measurable** — RefSeq returns early at
  `mapper_deleted_gap_cdna_index` / `refseq_deleted_edit_cdna_index` *before* the
  dispatch I changed; the offset arithmetic is trivial + short-circuits.
  Implemented, measured zero, REVERTED. (2) The only real C-shaped waste = the
  per-call `collect+sort` of `cdna_mapper_segments`; deduped via a `Cow` helper
  `sorted_mapper_segments` (**commit `30dfa09`**) but the gain (~0.05–0.1s) is
  **below the ~0.2–0.4s noise floor** — kept as alloc/dup cleanup only. The
  "~2–4s for C" below is WRONG; the big stages are dominated by the coord walk +
  necessary VEP-spec string formatting, not cheap redundancy.
- **E (AF clones): STALE, not 1–2s.** The `to_string()` storm was already fixed
  (borrowed `&str` keys in `freq_data`). ~1 dup `to_string()` per *populated*
  column remains → **~0.15–0.35s**, all in `colocated_fields`, none in `csq_format`.
- **B (HGVSc gating): FALSE PREMISE, ~0s.** Intergenic/up/down already skip
  `format_hgvsc`; intronic HGVSc is VEP-correct output. No output pays full
  `coord_map` then discards a non-empty result. Drop it (verify via
  `hgvsc_calls`≈`hgvsc_outputs` counter if you want empirical closure).
- **D (fallback round-trip): ~0.1–0.15s**, overlaps C. Small.
- **A (materialize-then-pick): 0 under `--everything`** (profiled mode). Unchanged.
- **hgvs_shift (NOT in the original order — the cleanest real win): DONE.**
  `build_hgvs_genomic_shift` was called twice per indel (fwd+rev) re-reading the
  same strand-independent ~1000bp flanks; split prelude (parse+reads, once) from
  per-strand build (**commit `e83ce59`**). Measured **0.580s → 0.464s** (~0.12s;
  low end — OS page cache softens the dup reads). Byte-identical w1+w4.

**Method lesson:** profile stage names mislead (nested timers; early-return paths
bypass the code you "optimize"). MEASURE each lever on a COMPLETE run — incomplete
profile captures show round `output_rows` (95000/275000/290000) vs the true
**307295** for chr4; trust only the exact-count runs.

**Where remaining engine time actually is** (if pushing further): the coord-walk
+ string-formatting in `transcript_output_materialize` / `transcript_hgvsc`, most
of which is necessary VEP-spec output generation. E (~0.15–0.35s) is the only
clean unclaimed lever left; everything else is sub-noise or correctness-risky.

*Everything below is the ORIGINAL 2026-06-26 analysis, retained for context but
SUPERSEDED by this banner where they conflict.*

---

## Why now

Parallel annotation is complete and correct (see `2026-06-26-grid-aligned-parallel-annotation.md`:
chr4/chr2 w{1,4} 0 HGNC mismatches, byte-identical, w4 2.4× / w8 3.0×). The lookup over-read
that capped scaling is fixed. **The pipeline is now ~98% engine-bound, and the engine is
parallelized across grid workers — so the per-worker engine time is the floor.** Every redundant
engine op is paid by *all* workers; cutting it multiplies across the fleet AND lowers the serial
w1 number. This is the highest-leverage remaining work.

## Root cause (from the prior analysis — verify against current code first)

Spec: `docs/superpowers/specs/2026-06-17-vep-everything-redundancy-analysis.md`
(memory: `vep-everything-redundancy-analysis`, `vep-output-gen-scale-and-coord-walks`,
`vep-annotation-bottlenecks`). Ensembl VEP `release/115` Perl is cloned at
`~/research/git/{ensembl-vep,ensembl-variation,ensembl}` (consequence+HGVS engine lives in
ensembl-variation/ensembl core, NOT ensembl-vep).

**Inverted loop nesting.** VEP loops transcript-outer / variant-inner and caches per-transcript
geometry ONCE on `Transcript->{_variation_effect_feature_cache}` (mapper, sorted exons, coord
table, CDS bounds), reused by every overlapping variant + allele. Ours loops
row(=position×single-ALT)-outer / transcript-inner and **rebuilds geometry per (variant×transcript)**.
The `vep-output-gen-scale-and-coord-walks` note measured **~2.09 genomic→cDNA exon-walks per
output**, different output flavors (HGVSc / cDNA / CDS / protein position) re-walking the same
geometry — partially memoizable.

## Fresh profile (chr4 serial, 2026-06-26, post-parallel — single-thread, all VEP_*_PROFILE)

```
annotate 11.18s   engine 10.93s (98%)   lookup_wait 0.47s   context_load 0.38s
ENGINE top-level (additive): evaluate_prepared 5.79 · csq_format 1.86 · colocated_fields 1.43 · collapse_pick_sort 0.87 · hgvs_shift 0.54
TX ENGINE (nested in evaluate_prepared):
  transcript_output_materialize 3.37   ← TOP
  transcript_hgvsc 2.16  (coord_map 1.55 · simple_fast 1.40 · fallback 0.63)   ← coord walks
  transcript_overlap_eval 1.66
```

Reproduce: `vep-perf-profiling` skill (chr4, `--workers 1`).

## Attack order (unchanged from prior analysis; behavior-preserving first)

1. **C — per-transcript geometry rebuild (TOP, ~2–4s, behavior-preserving).** ⚠️ DISPROVEN/REVERTED — see banner (measured ~0). Cache the per-transcript
   geometry (sorted exons, cdna mapper/coord table, CDS bounds) once per transcript per buffer and
   reuse across all overlapping variants/alleles, instead of rebuilding in `transcript_output_materialize`
   + the `transcript_hgvsc` coord walks. Subsumes per-ALT recompute (#2) and half of D. There is
   already a `GeometryCache` (OnceLock) on `TranscriptFeature.cdna_coords_cache` and it is `Sync`
   (safe across grid workers) — **check what it currently memoizes vs what's still rebuilt per output**;
   the win is widening that cache to the full per-output coord walk.
2. **E — frequency string clones (~1–2s, Rust-only waste).** ⚠️ STALE — measured ~0.15–0.35s (storm already fixed); see banner. VEP also re-parses AF per row, so only our
   3–4× AF string clones are genuine waste. Cut the clones (the `b_af` 27-column build path).
3. **B — HGVSc gating.** VEP hard-gates HGVSc behind `within_feature` before `coord_map`; we pay
   `coord_map` on non-cDNA outputs. **Confirm the counter first** (how many coord_maps are on
   non-cDNA outputs) before optimizing.
4. **D — fallback string round-trip (~0.3–0.8s).** `transcript_hgvsc_fallback` 0.63s; half subsumed by C.
5. **A — materialize-then-pick.** **0 under `--everything`** (VEP materializes all transcripts too;
   the 4.74M-row/0-mismatch parity gate proves the 12× fan-out is VEP-correct). Only ~3–5s IF the
   product switches to `--pick`. Do NOT touch unless pick becomes the target.

## Verification gate (non-negotiable — C/B/D touch HGVSc correctness)

- After EACH step: re-profile (`vep-perf-profiling`) + **e2e parity WITHOUT `--skip-compare`**:
  `run_annotation_fast.py chr4 --cache merged --backend lance --workers 1 --force` →
  0 mismatches on all 86 CSQ fields (esp. HGVSC/HGVSP/cDNA_position/CDS_position/Protein_position).
- Then chr4 **w4** (the parallel path) must stay byte-identical too — the geometry cache is shared
  read-only across grid workers via the `OnceLock`, so confirm no cross-worker races (run w4 ×3).
- `cargo test -p datafusion-bio-function-vep --lib` + `cargo clippy -- -D warnings` green.
- Commit only when parity holds; never commit a faster-but-divergent engine.

## Parallel-path interactions to keep in mind

- The grid workers share `Arc<SharedContigAnnotationContext>` incl. `base_transcripts` whose
  `cdna_coords_cache: GeometryCache(OnceLock)` is concurrently filled — widening this cache is
  safe (OnceLock) but verify the value computed is identical regardless of which worker fills it
  first (it is geometry, position-independent — should be).
- Engine wins lower BOTH the serial w1 floor and every grid worker, so the speedup compounds with
  the existing 2.4× rather than competing with it.

## Quick-start for the next session

**Superseded by the 2026-06-27 banner.** C and hgvs_shift are done/closed; do NOT
re-chase C (measured zero — see banner). If continuing:

1. `vep-perf-profiling` chr4 w1 → confirm on a COMPLETE run (`output_rows=307295`).
2. **E** is the only clean unclaimed lever (~0.15–0.35s): in `frequency_fields()`
   (`annotate_provider.rs` ~:2047–2059) drop the duplicate `to_string()` — move
   the owned `chosen` into `per_column` instead of re-stringifying. Parity gate; commit.
3. Beyond E, remaining engine time is mostly necessary VEP-spec coord-walk +
   string formatting; expect diminishing, sub-noise returns. Re-scope before investing.

*(Original step 3 "Implement C" was attempted and reverted — zero measurable gain.)*
