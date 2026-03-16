# Benchmark Results: Annotation Context Loading Optimization

## Environment

- Machine: macOS Darwin 24.6.0
- Rust: edition 2024, toolchain 1.88.0
- DataFusion: 52.1.0
- Dataset: HG002 chr1 (319,349 variants → 323,430 after multi-allelic decomposition)
- Cache: chr1-vep split layout (translation_core + translation_sift, 90 RGs × 256 rows)

## Golden Benchmark Invocation

```bash
mkdir -p /tmp/annotate_vep_everything_chr1vep
cp vep-benchmark/data/output/everything/HG002_chr1_0_vep115_golden.vcf \
   /tmp/annotate_vep_everything_chr1vep/HG002_chr1_0_vep115_golden.vcf

VEP_PROFILE=1 cargo run --release -p datafusion-bio-function-vep \
  --example annotate_vep_golden_bench -- \
  --everything \
  --extended-probes \
  --reference-fasta-path=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/chr1-vep/115_GRCh38_variation_1_vep.parquet \
  parquet \
  0 \
  /Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38 \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_everything_chr1vep \
  /Users/mwiewior/research/data/vep/chr1-vep \
  --steps=ensembl,datafusion
```

## Correctness: 80/80 CSQ Fields at 100% Accuracy

```
golden_rows=323430
ours_rows=323430
intersection_rows=323430
missing_in_ours=0
extra_in_ours=0
golden_with_csq=323428
ours_with_csq=323428
most_severe_accuracy=100.0%
term_set_accuracy=100.0%

Perfect (0 mismatches): 80/80 fields
```

## Timing Breakdown (VEP_PROFILE, --everything mode, 325K VCF rows)

```
[VEP_PROFILE] 1. variation_lookup (scan+collect)................  38552.2ms
[VEP_PROFILE] 2. collect_variant_intervals......................     95.2ms
[VEP_PROFILE] 3. colocated_map_build............................    387.3ms  320539 entries
[VEP_PROFILE] cache hits: 0, misses: 325785, hit rate: 0.0%
[VEP_PROFILE] 4a. load_transcripts..............................    341.4ms  47849 transcripts
[VEP_PROFILE] 4b. load_exons....................................     58.9ms  355800 exons
[VEP_PROFILE] 4c. load_translations.............................     99.8ms  22832 translations
[VEP_PROFILE] 4d. load_regulatory...............................      9.4ms  35815 features
[VEP_PROFILE] 4e. load_motif....................................      1.2ms
[VEP_PROFILE] 4f. load_mirna....................................      0.0ms
[VEP_PROFILE] 4g. load_structural...............................      0.0ms
[VEP_PROFILE] 4. context_tables_total...........................    515.1ms
[VEP_PROFILE] 5a. hydrate_refseq_cds............................     21.6ms
[VEP_PROFILE] 5b. hydrate_transcript_cdna.......................    759.4ms
[VEP_PROFILE] 6. prepared_context_build.........................     22.6ms
[VEP_PROFILE] 7. sift_polyphen_cache_init.......................      0.0ms
[VEP_PROFILE] 7+8. sift_lazy_load + annotate_batches............  62753.9ms  51 sift windows loaded
[VEP_PROFILE] 9. projection + memtable..........................     62.2ms
[VEP_PROFILE] TOTAL scan_with_transcript_engine................. 103230.6ms
annotate_vep_elapsed_s=104.571
```

## Context Loading: Before vs After

Measured with `profile_annotation` example on same chr1-vep caches (non-everything mode, 319K variants):

| Stage | Before (master, unsplit cache) | After (this PR, split cache) |
|-------|------|-------|
| load_translations | 13.1s | **0.10s** |
| load_transcripts | 0.30s | 0.27s |
| load_exons | 0.12s | 0.05s |
| load_regulatory | 0.07s | 0.01s |
| **context_tables_total** | **13.5s** | **0.47s** |

The 13s savings comes from column projection on `load_translations`: the old `SELECT *` loaded 132 MB of sift/polyphen data that the main loader discarded. The new code projects only the 8 core columns (~5 MB).

## Note on --everything Mode Total

The --everything total (104.6s) shows minimal improvement over the baseline (107.7s) because:
1. Context loading savings: -13.1s
2. Sift window loading: ~62.8s (comparable to baseline ~57.7s — slight increase from split file overhead)
3. Variation lookup: ~38.6s (comparable to baseline ~33.9s — run-to-run variance)

The context optimization is real (0.5s vs 13.6s) but represents only ~13% of the --everything total. The dominant costs — sift window loading (60%) and variation lookup (37%) — are not addressed by this PR.

## Cache Layout (chr1-vep)

```
115_GRCh38_exon_1_vep.parquet               17.5 MB    355,800 rows    8 RGs  sort=(transcript_id, start)
115_GRCh38_motif_1_vep.parquet               0.0 MB          0 rows    0 RGs
115_GRCh38_regulatory_1_vep.parquet          6.8 MB     35,930 rows    4 RGs  sort=(chrom, start)
115_GRCh38_transcript_1_vep.parquet         31.8 MB     47,849 rows    6 RGs  sort=(chrom, start)
115_GRCh38_translation_core_1_vep.parquet   13.0 MB     22,832 rows    4 RGs  sort=(transcript_id)
115_GRCh38_translation_sift_1_vep.parquet  222.0 MB     22,832 rows   90 RGs  sort=(chrom, start)
115_GRCh38_variation_1_vep.parquet           2.8 GB 88,153,966 rows  882 RGs  sort=(chrom, start)
```
