# Kickoff — golden-VEP parity for the plugin caches

## Task

Run the `merged_plugins` comparison profile against the golden Ensembl VEP
output and report per-field CSQ match rates for the five plugin fields sets
(ClinVar, SpliceAI, CADD, AlphaMissense, dbNSFP).

**This has never been run.** Everything validated so far checks that a shard is
*internally well-formed* — sorted within tier, tiers grouped, no duplicate probe
keys. Nothing yet checks that the annotations are **correct**. That is this task.

## Why it matters

The plugin caches exist so vepyr reproduces Ensembl VEP's plugin output. A cache
can be perfectly well-formed and still carry wrong values — wrong key
normalization, a bad `ingest_sql` transform, a mis-set `assume_unique`, an
anchor-shift off by one. Only a field-level comparison against real VEP output
catches that.

## Start here

| thing | where |
|---|---|
| Comparison harness + `merged_plugins` profile | `vepyr` PR #48, branch `plugins-fixes` |
| Engine + built caches | `datafusion-bio-functions` PR #217, branch `plugins-fixes` |
| Multi-part sources, build notebook | `vepyr` PR #49, branch `plugin-cache-build-validation` |
| Manifests (CADD is two `[[source]]` parts) | `vepyr-plugins` PR #4 |
| How to build/slice a chromosome | `docs/superpowers/handovers/2026-08-25-plugin-cache-testing-handover.md` |
| Engine memory investigation (context only) | `docs/superpowers/handovers/2026-08-24-cadd-tier-join-reservation-leak-RESOLVED.md` |

Roughly:

```bash
cd <vepyr>/e2e-testing/scripts
uv run python run_comparison.py --release 116 --profile merged_plugins --chroms 21
```

`merged_plugins` attaches the plugin cache and compares against a VEP reference
built with the same five plugins; `merged_plugins_base` attaches nothing and
isolates whether a core-field difference predates the plugin machinery. Confirm
the profile resolves before running anything long — it needs the plugin cache at
`$DATA/cache/plugin_cache_116` (or `--plugin-cache`) and the reference at
`output/116/HG002_annotated_wgs_everything_hgvs_merged_clinvar_spliceai_cadd_am_dbnsfp.vcf[.gz]`.

## Watch out for

1. **chr21 first, then chr1.** Both are already built and validated, so a
   mismatch is a comparison or annotation problem, not a cache-build problem.
2. **`verify_parity_gate.py` refuses plugin profiles by design** — it pins the
   Ensembl core CSQ contract and a plugin run emits fields outside it. Use the
   comparison report, not the gate. Extending the gate is separate work.
3. **Representation-only differences are expected and already handled.** PR #48
   counts decimal padding (`0` vs `0.00`) and VEP's `.` marker separately from
   real mismatches. If plugin float fields show ~0% match, check you are on that
   branch before investigating the cache.
4. **Do not pipe per-plugin output through `tail -N`** — it truncates the
   diagnostics and has already produced one wrong conclusion here.
5. **A prior 5-plugin chr22 run reported 37/38 fields at 100%**, but with a
   different comparator (the old index-pairing one). Treat those numbers as
   unverified until reproduced on the current pairing logic.

## Done when

- `merged_plugins` runs clean on chr21 and chr1, with per-field match rates
  reported for every plugin CSQ field.
- Every field below 100% is explained: a real cache defect, a known VEP
  inconsistency, or a comparator artifact — named, not hand-waved.
- Findings posted to PR #217 (the caches) and/or #48 (the comparator),
  whichever the cause implicates.

## Ground rules

Measure before concluding. The engine-side investigation in this repo burned
five plausible-but-wrong hypotheses that each died to the next measurement; the
ones that survived came from instrumenting and reading numbers. If a parity gap
appears, find the specific variant and inspect it in the raw source and the
shard before theorising about the transform.
