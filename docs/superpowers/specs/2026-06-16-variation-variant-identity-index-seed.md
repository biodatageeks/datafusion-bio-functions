# Variation Variant-Identity Index — Brainstorm Seed

> Status: **pre-brainstorm seed** (not a finished design). Start a `superpowers:brainstorming` session from here; brainstorm BOTH options A and B below.

## Problem & headline metric (measured, chr1 HG002, single-path lance, true 1-CPU)

The cold-tier variation lookup resolves by **position**, takes **all alleles at a matched position**, then `colocated_match` filters to the input variant's allele — so most taken rows are decoded only to be discarded.

```
rows taken/decoded (exact_match_calls) = 775,056
useful (primary_matches / colocated_entries) ≈ 311k–341k   (~44%)
WASTED (taken+decoded → discarded by allele filter) ≈ 434k  (~56%)
cold_tier: probes=727363 matches=315586 position_misses=233304 not_covered=178473
variation_take = 79 calls, 17.4s @ 1-CPU (8.0s uncapped); the take/decode is the #1 lookup cost.
```

**Goal:** stop taking the ~56% wrong-allele rows by filtering on **variant identity** before the take.

## Corrected facts (don't re-derive — these were verified this session)

- **Single take path now.** The `warm`/`cold` profiler labels are **vestigial** from the old two-tier (warm-fjall + cold-parquet) design. `warm matches=0, not_covered=727363, warm_*_load=0s` → warm tier does nothing. Everything flows through one lance `take_rows`.
- **The variant bloom is INACTIVE**: `variant_bloom checks=0 loaded=0 entries=0` (despite `index_mode=posidx_bloom`). So there is **zero variant-identity filtering today** — an `variation.variant_bloom_index` dir exists but isn't built/wired for this cache.
- Filtering today is **position-level only**, via the loaded start BTree (`position_index loaded=1 rows=88,153,966`): "is there ANY variant at this position?" (231k miss, 178k not-covered → no take). A wrong-allele-at-existing-position variant is **not** filtered → it takes all alleles.
- The row gap vs the lance_sandbox benchmark (vepyr 735k vs 402k) is **real distinct rows**, NOT cross-batch re-takes — the row-id dedup-cache hypothesis was implemented and **falsified** (`cache_hits=5`/79 batches; see `docs/superpowers/plans/2026-06-16-variation-rowid-dedup-cache.md`). Take-side dedup is a dead end.
- Per-row take cost is at parity with the sandbox (~23µs/1k cold); the lever is **row count (breadth)**, not take speed. (Clustering = a separate, speed-only lever; measure in `lance_encoding_sandbox` if pursued — orthogonal.)

## The two options to brainstorm

**A — variant-identity NEGATIVE filter (safe, YAGNI).** Build/activate a variant-level bloom keyed on the normalized variant. Probe: bloom "definitely absent" → **skip take**; "maybe" → **fall back to today's position-take + allele-filter**. A false positive just does the old thing ⇒ **cannot break parity**. Captures most of the 56% for non-existing / no-allele-match variants. Smaller change.

**B — variant-keyed PRIMARY resolution (aggressive).** Resolve `row_id`s by variant identity → return only allele-matching rows. Bigger win (also cuts `exact_match_calls`/`colocated_match` work), but **must exactly replicate VEP's match set** or it silently drops `Existing_variation`/frequencies.

## The parity landmine (applies to both, critical for B)

VEP matches an input variant against existing variants by **output allele OR unshifted/shifted allele representations** (this is *why* `extended_probes` exists — see `cache_exec.rs` `extended_probes`, "needed for deletion matching"). Any variant key/bloom must cover **all** allele representations VEP would match (output + unshifted, normalized), per chrom/strand. Get it wrong → missed colocated annotations.

## Key code & validation

- Variation lookup / cold-tier take: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs` (lance flush ~L3626-3661; `match`/`colocated_match`; the `variant_bloom`/`position_index` profile lines).
- Existing builder for the variant bloom: example `build_variation_variant_bloom_index` + `VEP_VARIATION_BLOOM_INDEX_DIR`.
- Allele matching semantics: `ColocatedData::variant_fields` / `matching_allele` / `matches_output_allele` in `annotate_provider.rs`.
- **Acceptance gate (non-negotiable):** chr1 (+ chr4) e2e CSQ diff = **100% on all shared fields**. Profile target: `variation_take` rows 775k → ~341k, decode/`exact_match_calls` proportional.

## Environment notes for the next session

- New-format (position-sliced sift) chr1 cache is swapped into `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged` (sift = key/sift/poly); restore via `/tmp/restore_chr1_legacy.sh` if needed.
- e2e (true 1-CPU): from `~/research/git/vepyr/e2e-testing/scripts`, `env -u VIRTUAL_ENV -u CONDA_PREFIX LANCE_CPU_THREADS=1 LANCE_IO_THREADS=1 RAYON_NUM_THREADS=1 VEP_PROFILE=1 VEP_LANCE_PROFILE=1 VEP_KV_PROFILE_DETAILED=1 uv run python run_annotation_fast.py chr1 --cache merged --forks 0 --force --backend lance`. Rebuild vepyr (path dep): `cd ~/research/git/vepyr && RUSTFLAGS="-C target-cpu=native" uv sync --reinstall-package vepyr`.
- e2e run can print a spurious `0 shared CSQ fields` (comparison-step artifact) — verify the output VCF is non-empty before trusting it.
