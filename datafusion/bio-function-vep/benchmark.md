# Benchmark Results: VEP Annotation (Fjall vs Parquet)

## Environment

- Machine: macOS Darwin 24.6.0 (Apple Silicon)
- Rust: edition 2024, toolchain 1.88.0
- DataFusion: 52.1.0
- fjall: 3.1.1

## Datasets

| Dataset | VCF | Variants | Contigs |
|---------|-----|----------|---------|
| HG002 chr1 | `vep-benchmark/data/HG002_chr1.vcf.gz` | 323,430 | 1 |
| HG002 WGS chr1-22 | `vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` | 4,096,123 | 22 |

## Cache Paths

```
# Parquet cache (partitioned per-chromosome)
PARQUET_CACHE=/Users/mwiewior/research/data/vep/wgs/parquet/115_GRCh38_vep

# Fjall cache (variation.fjall + translation_sift.fjall + symlinked context parquets)
FJALL_CACHE=/Users/mwiewior/research/data/vep/wgs/fjall/115_GRCh38_vep

# Reference FASTA (required for --everything / --hgvs)
REF_FASTA=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Ensembl VEP offline cache (for golden truth generation)
VEP_CACHE_DIR=/Users/mwiewior/research/data/vep/homo_sapiens
```

## Cache Directory Layout

```
# Parquet cache
115_GRCh38_vep/
├── variation/chr{1..22,X,Y,other}.parquet     # 37 GB
├── transcript/chr*.parquet                     # 352 MB
├── exon/chr*.parquet                           # 196 MB
├── translation_core/chr*.parquet               # 145 MB
├── translation_sift/chr*.parquet               # 2.4 GB
├── regulatory/chr*.parquet                     # 72 MB
└── motif/chr*.parquet                          # 900 KB

# Fjall cache (variation + sift in fjall, context from parquet via symlinks)
115_GRCh38_vep/
├── variation.fjall/                            # 109 GB (1.17B variants, 1.11B positions)
├── translation_sift.fjall/                     # 26 GB (SIFT/PolyPhen predictions)
├── variation -> ../parquet/115_GRCh38_vep/variation
├── transcript -> ../parquet/115_GRCh38_vep/transcript
├── exon -> ../parquet/115_GRCh38_vep/exon
├── translation_core -> ../parquet/115_GRCh38_vep/translation_core
├── translation_sift -> ../parquet/115_GRCh38_vep/translation_sift
├── regulatory -> ../parquet/115_GRCh38_vep/regulatory
└── motif -> ../parquet/115_GRCh38_vep/motif
```

## Correctness

Both backends verified 100% accurate against Ensembl VEP 115 golden truth (chr1, 323,430 variants):

```
Variants: 323430 golden, 323430 ours, 323430 matched, 0 missing, 0 extra
Accuracy: most_severe=100.00%, term_set=100.00%
CSQ per-field: 2,997,504 entries across 80 fields — Perfect (0 mismatches): 80/80 fields
```

## End-to-End Performance

### HG002 chr1 (323,430 variants, `--everything`)

| Backend | Time | Throughput |
|---------|------|------------|
| **Fjall** | **20s** | **16K variants/s** |
| **Parquet** | **70s** | **4.6K variants/s** |

### HG002 WGS chr1-22 (4,096,123 variants, `--everything`)

| Backend | Time | Throughput |
|---------|------|------------|
| **Fjall** | **365s (6.1 min)** | **11.2K variants/s** |

### Fjall Per-Contig Breakdown (WGS chr1-22)

| Contig | Variants | Time | Throughput |
|--------|----------|------|------------|
| chr1 | 323,430 | 24.5s | 13.2K/s |
| chr2 | 331,324 | 29.1s | 11.4K/s |
| chr3 | 288,531 | 26.3s | 11.0K/s |
| chr4 | 307,295 | 23.8s | 12.9K/s |
| chr5 | 264,411 | 20.8s | 12.7K/s |
| chr6 | 271,966 | 21.6s | 12.6K/s |
| chr7 | 234,522 | 20.1s | 11.7K/s |
| chr8 | 225,240 | 18.9s | 11.9K/s |
| chr9 | 176,111 | 15.8s | 11.1K/s |
| chr10 | 213,466 | 18.6s | 11.5K/s |
| chr11 | 206,822 | 17.6s | 11.8K/s |
| chr12 | 197,815 | 17.0s | 11.6K/s |
| chr13 | 161,419 | 12.2s | 13.2K/s |
| chr14 | 133,199 | 12.5s | 10.7K/s |
| chr15 | 125,179 | 12.2s | 10.3K/s |
| chr16 | 123,358 | 11.9s | 10.4K/s |
| chr17 | 108,376 | 12.0s | 9.0K/s |
| chr18 | 119,383 | 10.8s | 11.1K/s |
| chr19 | 90,699 | 10.6s | 8.6K/s |
| chr20 | 86,904 | 8.8s | 9.9K/s |
| chr21 | 55,812 | 6.2s | 9.0K/s |
| chr22 | 50,861 | 7.0s | 7.3K/s |
| **Total** | **4,096,123** | **365s** | **11.2K/s** |

## Golden Benchmark Invocation

```bash
# Build (release, with fjall support)
cargo build --release --features kv-cache --example annotate_vep_golden_bench

# ── Fjall backend ────────────────────────────────────────────────

# Chr1 — correctness check (reuses existing golden VEP output)
VEP_PROFILE=1 cargo run --release --features kv-cache --example annotate_vep_golden_bench -- \
  --everything --use-fjall --extended-probes \
  --steps=ensembl,datafusion \
  --reference-fasta-path=$REF_FASTA \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  $FJALL_CACHE \
  parquet 0 /dev/null \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  vep-benchmark/data/output/everything

# Chr1 — annotation only (no comparison, fastest)
cargo run --release --features kv-cache --example annotate_vep_golden_bench -- \
  --everything --use-fjall --extended-probes \
  --steps=datafusion \
  --reference-fasta-path=$REF_FASTA \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  $FJALL_CACHE \
  parquet 0 /dev/null \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  vep-benchmark/data/output/everything

# WGS chr1-22 — annotation only
cargo run --release --features kv-cache --example annotate_vep_golden_bench -- \
  --everything --use-fjall --extended-probes \
  --steps=datafusion \
  --reference-fasta-path=$REF_FASTA \
  vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  $FJALL_CACHE \
  parquet 0 /dev/null \
  vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  vep-benchmark/data/output/wgs_fjall

# ── Parquet backend ──────────────────────────────────────────────

# Chr1 — correctness check
VEP_PROFILE=1 cargo run --release --features kv-cache --example annotate_vep_golden_bench -- \
  --everything --extended-probes \
  --steps=ensembl,datafusion \
  --reference-fasta-path=$REF_FASTA \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  $PARQUET_CACHE \
  parquet 0 /dev/null \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  vep-benchmark/data/output/everything

# WGS chr1-22 — annotation only
cargo run --release --features kv-cache --example annotate_vep_golden_bench -- \
  --everything --extended-probes \
  --steps=datafusion \
  --reference-fasta-path=$REF_FASTA \
  vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  $PARQUET_CACHE \
  parquet 0 /dev/null \
  vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  vep-benchmark/data/output/wgs_parquet
```

## Cache Build (Fjall)

```bash
# 1/2: Variation cache (37 GB parquet -> 109 GB fjall, ~76 min with 8 threads)
cargo run --release --features kv-cache --example load_cache_full -- \
  $PARQUET_CACHE/variation/ \
  $FJALL_CACHE/variation.fjall \
  8

# 2/2: SIFT/PolyPhen (per-chromosome parquets -> standalone fjall DB, ~5 min)
for f in $PARQUET_CACHE/translation_sift/*.parquet; do
  cargo run --release --features kv-cache --example load_sift_cache -- \
    "$f" $FJALL_CACHE/translation_sift.fjall
done

# Context tables: symlink from parquet cache (no conversion needed)
for dir in variation transcript exon translation_core translation_sift regulatory motif; do
  ln -sfn $PARQUET_CACHE/$dir $FJALL_CACHE/$dir
done
```

## Architecture

Both backends use the same streaming contig-by-contig pipeline (`ContigAnnotationExec`):

```
VCF (tabix-indexed) → contig discovery → per-contig:
  ├── variation lookup [PARQUET: COITree interval join | FJALL: point lookups]
  ├── context load (transcripts, exons, translations from parquet)
  ├── window-based annotation (1000 batches per window)
  │   ├── HGVS hydration (reference FASTA)
  │   ├── SIFT/PolyPhen [PARQUET: SiftDirectReader | FJALL: SiftKvStore]
  │   └── consequence calling + CSQ assembly
  └── memory reclamation (context + lookup data freed after contig)
```

The `use_fjall` option (`options_json: {"use_fjall": true}`) swaps variation lookup and SIFT
to fjall point lookups while keeping context tables on partitioned parquet.
