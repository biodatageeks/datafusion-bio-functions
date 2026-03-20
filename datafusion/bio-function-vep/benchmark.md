# Benchmark Results: VEP Annotation (Fjall vs Parquet)

## Environment

- Machine: macOS Darwin 24.6.0 (Apple Silicon)
- Rust: edition 2024, toolchain 1.88.0
- DataFusion: 52.1.0
- fjall: 3.1.1
- Dataset: HG002 chr1 (319,349 variants -> 323,430 after multi-allelic decomposition)
- Mode: `--everything` (80-field CSQ schema, HGVSc/HGVSp, SIFT, PolyPhen, DOMAINS, all AFs)

## Correctness

Both backends verified 100% accurate against Ensembl VEP 115 golden truth:

```
Variants: 323430 golden, 323430 ours, 323430 matched, 0 missing, 0 extra
Accuracy: most_severe=100.00%, term_set=100.00%
CSQ per-field: 2,997,504 entries across 80 fields — Perfect (0 mismatches): 80/80 fields
```

## End-to-End Performance

| Backend | Cold (after `sudo purge`) | Warm | Speedup vs Parquet |
|---------|--------------------------|------|--------------------|
| **Fjall** (all keyspaces) | **28s** | **19s** | **3.5x** |
| **Parquet** | **70s** | **70s** | baseline |

## Fjall Timing Breakdown (Cold vs Warm)

Profiled with `VEP_PROFILE=1 VEP_KV_PROFILE=1`, `--everything` mode, 323,430 variants.

| Stage | Warm | Cold | Delta | % of cold gap |
|-------|------|------|-------|---------------|
| **KV lookup** (variation_lookup) | 3.9s | **11.2s** | **+7.3s** | **87%** |
| Colocated map build | 0.4s | 0.4s | 0 | |
| Context tables total | 0.5s | **0.9s** | +0.4s | 5% |
| - load_transcripts (parquet) | 0.3s | 0.4s | +0.1s | |
| - load_exons (kv) | 0.06s | 0.1s | +0.04s | |
| - load_translations (kv) | 0.09s | **0.5s** | +0.4s | |
| Hydrate transcript cDNA (FASTA) | 0.6s | **0.8s** | +0.2s | 3% |
| **Annotation engine** | 13.3s | 12.9s | -0.4s | (noise) |
| Other | 1.2s | 1.4s | +0.2s | |
| **TOTAL (engine)** | **18.8s** | **26.4s** | **+7.6s** | |
| **TOTAL (wall)** | **19.9s** | **27.6s** | **+7.7s** | |

### Key observations

- **Annotation engine (13s, 67%)** is the dominant cost — pure CPU, no I/O optimization possible
- **KV lookup cold penalty (+7.3s)** is the main cold-start bottleneck — fjall data block I/O
- Context tables via fjall KV (exons + translations) load in <0.5s warm, <1s cold
- SIFT/PolyPhen via fjall keyspace: 0ms (vs 23.7s parquet windowed loading in baseline)

## Fjall Configuration (Design Spec Decisions 1-7)

```rust
// Cold-start-optimized keyspace options:
.filter_block_partitioning_policy(PartitioningPolicy::all(false))   // Decision 1
.index_block_partitioning_policy(PartitioningPolicy::all(false))    // Decision 1
.filter_block_pinning_policy(PinningPolicy::all(true))              // Decision 1
.index_block_pinning_policy(PinningPolicy::all(true))               // Decision 1
.expect_point_read_hits(true)                                       // Decision 2
.data_block_hash_ratio_policy(HashRatioPolicy::all(0.75))           // Decision 3
.data_block_size_policy(BlockSizePolicy::all(8 * 1024))             // Decision 4
.filter_policy(FilterPolicy::all(FilterPolicyEntry::Bloom(          // Decision 6
    BloomConstructionPolicy::FalsePositiveRate(0.0001))))
.cache_size(1024 * 1024 * 1024)                                     // Decision 7: 1 GB
```

Post-load `major_compact()` ensures bloom filters and optimized LSM levels from first read.

## Fjall Keyspaces (chr1)

| Keyspace | Size | Entries | Key | Purpose |
|----------|------|---------|-----|---------|
| data (variation) | 7.3 GB | 83M positions | `[chrom_code 2B][start 8B]` | Position-keyed variant lookup |
| sift | 2.4 GB | 15K transcripts | `transcript_id` | SIFT/PolyPhen predictions |
| translations | 128 MB | 22K transcripts | `transcript_id` | CDS sequence, protein features |
| exons | 7.8 MB | 22K transcripts | `transcript_id` | Exon boundaries per transcript |
| meta | 132 KB | — | — | Schema, zstd dict, format version |

## Cache Build Time

```
Variation:    206s (4 threads, with major compaction)
SIFT:          24s
Exons:         <1s
Translations:  <1s
Total:        ~4 min
```

## Remaining Bottlenecks

| Stage | Time (warm) | % | Addressable by |
|-------|-------------|---|----------------|
| Annotation engine | 13.3s | 67% | Parallelization (per-chromosome partitions) |
| KV lookup | 3.9s | 20% | Already fast; cold-start prefetch for +7s gap |
| Hydrate cDNA | 0.6s | 3% | Background FASTA prefetch |
| Context tables | 0.5s | 3% | Already fast via KV |
| Other | 1.5s | 7% | Session setup, VCF loading |

## Golden Benchmark Invocation

```bash
# Fjall backend
./target/release/examples/annotate_vep_golden_bench \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /path/to/chr1-vep/fjall_variation_cache \
  fjall 0 \
  /path/to/homo_sapiens \
  /tmp/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_chr1 \
  /path/to/chr1-vep \
  --steps=ensembl,datafusion --extended-probes --everything \
  --reference-fasta-path=/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Parquet backend
./target/release/examples/annotate_vep_golden_bench \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /path/to/chr1-vep/115_GRCh38_variation_1_vep.parquet \
  parquet 0 \
  /path/to/homo_sapiens \
  /tmp/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_chr1_parquet \
  /path/to/chr1-vep \
  --steps=ensembl,datafusion --extended-probes --everything \
  --reference-fasta-path=/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

## Cache Regeneration

```bash
# 1/4: Variation cache
./target/release/examples/load_cache_full \
  /path/to/115_GRCh38_variation_1_vep.parquet \
  /path/to/fjall_variation_cache 4

# 2/4: SIFT/PolyPhen
./target/release/examples/load_sift_cache \
  /path/to/115_GRCh38_translation_sift_1_vep.parquet \
  /path/to/fjall_variation_cache

# 3/4 + 4/4: Exons + Translations
./target/release/examples/load_context_kv \
  /path/to/fjall_variation_cache \
  /path/to/115_GRCh38_exon_1_vep.parquet \
  /path/to/115_GRCh38_translation_core_1_vep.parquet
```
