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
[VEP_PROFILE] 1. variation_lookup (scan+collect)................  38642.3ms  321713 VCF rows, 10840 batches
[VEP_PROFILE] 2. collect_variant_intervals......................     94.0ms
[VEP_PROFILE] 3. colocated_map_build............................    365.7ms  314926 entries
[VEP_PROFILE] cache hits: 0, misses: 321713, hit rate: 0.0%
[VEP_PROFILE] 4a. load_transcripts..............................    301.2ms  47849 transcripts
[VEP_PROFILE] 4b. load_exons....................................     49.2ms  355800 exons
[VEP_PROFILE] 4c. load_translations.............................     93.6ms  22832 translations
[VEP_PROFILE] 4d. load_regulatory...............................     10.9ms  35815 features
[VEP_PROFILE] 4e. load_motif....................................      1.1ms
[VEP_PROFILE] 4f. load_mirna....................................      0.0ms
[VEP_PROFILE] 4g. load_structural...............................      0.0ms
[VEP_PROFILE] 4. context_tables_total...........................    461.1ms
[VEP_PROFILE] 5a. hydrate_refseq_cds............................     20.6ms  0 hydrated
[VEP_PROFILE] 5b. hydrate_transcript_cdna.......................    895.9ms
[VEP_PROFILE] 6. prepared_context_build.........................     25.5ms  1 tx_trees chroms
[VEP_PROFILE] sift_direct_reader: enabled (bypassing DataFusion SQL)
[VEP_PROFILE] 7. sift_polyphen_cache_init.......................      1.7ms
[VEP_PROFILE] 7a. sift_lazy_load_only...........................  23725.0ms
[VEP_PROFILE] 7b. annotate_batches_only.........................   8395.0ms
[VEP_PROFILE] 7+8. sift_lazy_load + annotate_batches............  36995.0ms  10840 batches, 321713 total rows, 51 sift windows loaded
[VEP_PROFILE] 9. projection + memtable..........................     34.4ms
[VEP_PROFILE] TOTAL scan_with_transcript_engine.................  77589.1ms  321713 VCF rows
Total time:   78525.4ms (78.53s)
Throughput:   4097 variants/sec
VCF output:   884.1ms
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
| load_translations | 13.1s | **0.09s** |
| load_transcripts | 0.30s | 0.30s |
| load_exons | 0.12s | 0.05s |
| load_regulatory | 0.07s | 0.01s |
| **context_tables_total** | **13.5s** | **0.46s** |

The 13s savings comes from column projection on `load_translations`: the old `SELECT *` loaded 132 MB of sift/polyphen data that the main loader discarded. The new code projects only the 8 core columns (~5 MB).

## --everything Mode: Before vs After

| Stage | Baseline (master, unsplit) | After (this PR, split cache) | Savings |
|-------|------|-------|---------|
| variation_lookup | 33.9s | 38.6s | (run-to-run variance) |
| context_tables_total | 13.6s | **0.5s** | **-13.1s** |
| sift_lazy_load | 48.0s | **23.7s** | **-24.3s** |
| annotate_batches | 8.7s | 8.4s | — |
| hydrate_transcript_cdna | 1.1s | 0.9s | -0.2s |
| VCF sink | — | 0.9s | (new) |
| other | 2.4s | 0.5s | -1.9s |
| **TOTAL** | **107.7s** | **78.5s** | **-29.2s (-27%)** |

### Cumulative optimization impact

| Optimization | Sift load | Context load | Total pipeline |
|---|---|---|---|
| **Baseline (master)** | 48.0s | 13.6s | **107.7s** |
| + Column projection | 48.0s | **0.5s** | ~94s |
| + Compact predictions (u8 enum + sorted Vec) | 35.8s | 0.5s | ~92s |
| + Zero-copy Arrow parsing | 26.1s | 0.5s | ~82s |
| + Direct parquet-rs reader (cached metadata) | **23.7s** | 0.5s | **78.5s** |

### Optimizations applied

1. **Column projection** on `load_translations`: skip 132 MB sift/polyphen data → 13.1s → 0.1s
2. **Compact sift predictions**: replace `HashMap<(i32, String), (String, f32)>` with sorted `Vec<CompactPrediction>` using u8-encoded amino acids and prediction types. Eliminates ~256M String allocations
3. **Zero-copy Arrow parsing** for sift: read `&str` directly from Arrow buffers and encode to u8 in-place, bypassing intermediate `ProteinPrediction` structs with heap-allocated Strings. Combined with (2): 48.0s → 26.1s
4. **Direct parquet-rs sift reader** with cached `ArrowReaderMetadata`: bypass DataFusion SQL planning for all 51 window queries. Read footer once, reuse metadata. 26.1s → 23.7s
5. **MissWorklist** with interval predicates for regulatory/motif/mirna/structural
6. **Rust-side transcript_id filter** for exons and translations
7. **Split translation layout** support (translation_core + translation_sift)

### Remaining bottlenecks

| Stage | Time | % of total | Addressable by |
|---|---|---|---|
| variation_lookup | 38.6s | 49% | fjall KV cache (separate PR) |
| sift I/O (parquet decompression) | 23.7s | 30% | ~21.5s is irreducible parquet I/O; single-pass bulk load could save ~2-3s more |
| annotate_batches | 8.4s | 11% | already efficient (~26μs/row) |
| VCF sink | 0.9s | 1% | — |
| context + hydrate + other | 1.9s | 2% | — |

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
