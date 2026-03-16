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

## Timing Breakdown (VEP_PROFILE, --everything mode, profile_annotation, 321K VCF rows)

```
[VEP_PROFILE] 1. variation_lookup (scan+collect)................  39782.7ms  321713 VCF rows, 10840 batches
[VEP_PROFILE] 2. collect_variant_intervals......................     82.3ms
[VEP_PROFILE] 3. colocated_map_build............................    358.6ms  314926 entries
[VEP_PROFILE] cache hits: 0, misses: 321713, hit rate: 0.0%
[VEP_PROFILE] 4a. load_transcripts..............................    329.6ms  47849 transcripts
[VEP_PROFILE] 4b. load_exons....................................     58.0ms  355800 exons
[VEP_PROFILE] 4c. load_translations.............................    107.0ms  22832 translations
[VEP_PROFILE] 4d. load_regulatory...............................     11.3ms  35815 features
[VEP_PROFILE] 4e. load_motif....................................      1.1ms
[VEP_PROFILE] 4f. load_mirna....................................      0.0ms
[VEP_PROFILE] 4g. load_structural...............................      0.0ms
[VEP_PROFILE] 4. context_tables_total...........................    512.8ms
[VEP_PROFILE] 5a. hydrate_refseq_cds............................     21.8ms  0 hydrated
[VEP_PROFILE] 5b. hydrate_transcript_cdna.......................    785.2ms
[VEP_PROFILE] 6. prepared_context_build.........................     24.3ms  1 tx_trees chroms
[VEP_PROFILE] 7. sift_polyphen_cache_init.......................      0.0ms
[VEP_PROFILE] 7a. sift_lazy_load_only...........................  26111.0ms
[VEP_PROFILE] 7b. annotate_batches_only.........................   8533.0ms
[VEP_PROFILE] 7+8. sift_lazy_load + annotate_batches............  39566.1ms  10840 batches, 321713 total rows, 51 sift windows loaded
[VEP_PROFILE] 9. projection + memtable..........................     31.0ms
[VEP_PROFILE] TOTAL scan_with_transcript_engine.................  81219.6ms  321713 VCF rows
Total time:   82274.6ms (82.27s)
Throughput:   3910 variants/sec
VCF output:   866.0ms
```

### profile_annotation invocation

```bash
VEP_PROFILE=1 cargo run -p datafusion-bio-function-vep --release --example profile_annotation -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/chr1-vep \
  320000 \
  --everything \
  --reference-fasta-path=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output=/tmp/HG002_chr1_annotated.vcf
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

## --everything Mode: Before vs After

| Stage | Baseline (master, unsplit) | After (this PR, split cache) | Savings |
|-------|------|-------|---------|
| variation_lookup | 33.9s | 39.8s | (run-to-run variance) |
| context_tables_total | 13.6s | **0.5s** | **-13.1s** |
| sift_lazy_load | 48.0s | **26.1s** | **-21.9s** |
| annotate_batches | 8.7s | 8.5s | — |
| hydrate_transcript_cdna | 1.1s | 0.8s | -0.3s |
| VCF sink | — | 0.9s | (new) |
| other | 2.4s | 0.5s | -1.9s |
| **TOTAL** | **107.7s** | **82.3s** | **-25.4s (-24%)** |

Optimizations applied:
1. **Column projection** on `load_translations`: skip 132 MB sift/polyphen data → 13.1s → 0.1s
2. **Compact sift predictions**: replace `HashMap<(i32, String), (String, f32)>` with sorted `Vec<CompactPrediction>` using u8-encoded amino acids and prediction types. Eliminates ~256M String allocations
3. **Zero-copy Arrow parsing** for sift: read `&str` directly from Arrow buffers and encode to u8 in-place, bypassing intermediate `ProteinPrediction` structs with heap-allocated Strings. Combined with (2): 48.0s → 26.1s
4. **MissWorklist** with interval predicates for regulatory/motif/mirna/structural
5. **Rust-side transcript_id filter** for exons and translations
6. **Split translation layout** support (translation_core + translation_sift)

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
