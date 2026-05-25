# v115 Tier 1 cache region design

**Spec source**: porting-tests/REVIEW_NOTES.md §B1, B2, B3.
**Pinned spec SHA**: 7ee8669.

## Regions included

| Region | Why | Source rows |
|---|---|---|
| chr21:25,000,000–26,000,000 | Canonical test region from upstream Perl `t/AnnotationSource_Cache_Transcript.t` — diverse transcripts (MRPL39, JAM2, GABPA) + regulatory features. | CORE Cache_Transcript / Cache_RegFeat / Parser_VCF rows; OBSERVABLE-SUBSET Runner / OutputFactory rows |
| chr21:8,960,000–9,010,000 | Variation `nastiness 1-4` matched-allele cases (REVIEW_NOTES §B1). Required for AnnotationSource_Cache_Variation rows 36-39. | Cache_Variation nastiness rows |

Total span: ~1.05 Mb across two slices, both on chr21.

For Tier 1 simplicity this fixture stores **whole chr21 parquet** (not just
the two slices) because the slice-and-merge is cheaper to do at test time
(via `start IN (RANGE …)` filters) than to maintain region-by-region
parquet files. The on-disk footprint stays bounded because chr21 is the
smallest autosome (~33 Mb of cache rows). MT is also included to support
OutputFactory_VCF's MT row (a CORE OBSERVABLE-SUBSET requirement).

## Curated content checklist (B3 acceptance)

- [ ] ≥1 MANE-select transcript
- [ ] ≥1 canonical transcript
- [ ] ≥1 miRNA stem-loop transcript (Cache_RegFeat scope)
- [ ] ≥1 MT-RefSeq transcript (OutputFactory_VCF MT case) — covered by `MT/` inclusion
- [ ] ≥1 merged-cache variant (OutputFactory cluster H REFSEQ_MATCH) — `gap-justified` for now; vepyr uses ensembl-only cache (refseq/merged out of scope per A1-style decisions)
- [ ] ≥1 regulatory feature with TRANSCRIPTION_FACTORS populated (Cache_RegFeat MOTIF rows)
- [ ] rs3989369 OR equivalent v115 variant with regulatory overlap (B2)

## Source

Built by `datafusion-bio-functions/scripts/port/build_v115_parquet.sh`,
which slices chr21 from a pre-built whole-genome vepyr parquet cache
(`$VEPYR_WHOLE_GENOME_PARQUET`, default
`/Users/wojtek/Documents/vepyr/_cache_v115/parquet/115_GRCh38_vep`) plus
the native Ensembl cache subset for the oracle (chr21 + MT + info.txt).

Deviation from the original Batch 0 plan note: the plan presupposed
`vepyr.build_cache(..., chromosomes=["21"])`, but that kwarg is not in
the installed vepyr version (verified 2026-05-26). Whole-genome build +
chr21 slice via file copy is the working alternative.

## Build script

`datafusion-bio-functions/scripts/port/build_v115_parquet.sh` (this branch).

## B1 / B2 / B3 outcomes

| Question | Outcome | Recorded |
|---|---|---|
| B1: extend region for variation nastiness 1-4? | EXTEND to include chr21:8.96-9.01 Mb (covered by whole-chr21 inclusion) | 2026-05-26 |
| B2: rs3989369 v115 regulatory overlap | pending verification (requires v115 cache; recorded after Task 11) | TBD post-cache-build |
| B3: merged/miRNA/MT richness | partial — miRNA + MT covered by whole-chr21 + MT inclusion; merged-cache variant is gap-justified (vepyr cache is ensembl-only) | 2026-05-26 |
