# datafusion-bio-function-vep

VEP-oriented DataFusion functions and benchmark tooling.

## Golden Benchmark: `annotate_vep` vs Ensembl VEP 115

### Overview

The `annotate_vep_golden_bench` example compares our `annotate_vep()` table function against Ensembl VEP 115 on real whole-genome sequencing data (HG002/NA24385, GIAB).

**Pipeline:**
1. Samples variants from a bgzipped VCF (or uses all rows with `sample_limit=0`).
2. Decomposes multi-allelic sites into single-alt records.
3. *(Optional)* Runs Ensembl VEP 115 in Docker to produce golden output.
4. Runs `annotate_vep(...)` on the same normalized input.
5. Compares per-variant CSQ strings across 74 fields and writes reports.

### Latest chr1 results (March 2026)

| Metric | Value |
|--------|-------|
| Input variants | 319,349 (323,430 after multi-allelic decomposition) |
| CSQ entries compared | 2,997,504 |
| **most_severe accuracy** | **100.00%** (323,430/323,430) |
| **term_set accuracy** | **100.00%** (323,430/323,430) |
| Discrepancies (variant-level) | 0 |
| Perfect CSQ fields (0 mismatches) | 53/74 |
| Runtime | ~84s |

**Per-field accuracy (74 CSQ fields):**

53 fields at 100.00% (0 mismatches): Allele, IMPACT, SYMBOL, Gene, Feature_type, Feature, BIOTYPE, EXON, HGVSc, HGVSp, cDNA_position, Codons, DISTANCE, STRAND, SYMBOL_SOURCE, HGNC_ID, MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE, TRANSCRIPTION_FACTORS, SOURCE, VARIANT_CLASS, CANONICAL, TSL, MANE_SELECT, MANE_PLUS_CLINICAL, ENSP, GENE_PHENO, CCDS, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, gnomADe_AFR, gnomADe_AMR, gnomADe_ASJ, gnomADe_EAS, gnomADe_FIN, gnomADe_MID, gnomADe_NFE, gnomADe_REMAINING, gnomADe_SAS, gnomADg_AFR, gnomADg_AMI, gnomADg_AMR, gnomADg_ASJ, gnomADg_EAS, gnomADg_FIN, gnomADg_MID, gnomADg_NFE, gnomADg_REMAINING, gnomADg_SAS

21 fields with remaining mismatches:

| Field | Matched | Mismatches | Accuracy | Root cause |
|-------|---------|------------|----------|------------|
| Consequence | 2,997,489 | 15 | 100.00% | Deletion at splice_acceptor: VEP has intron_variant, we have PPT |
| INTRON | 2,997,485 | 19 | 100.00% | Splice-site deletion: VEP suppresses intron number |
| CDS_position | 2,997,502 | 2 | 100.00% | `?-N` prefix for cds_start_nf edge case |
| Protein_position | 2,997,502 | 2 | 100.00% | Same `?-N` edge case |
| Amino_acids | 2,997,461 | 43 | 100.00% | Frameshift boundary amino acid formatting |
| FLAGS | 2,997,076 | 428 | 99.99% | Sort order: VEP Perl uses non-deterministic hash |
| Existing_variation | 2,685,888 | 311,616 | 89.60% | Shifted indel co-located matching |
| AF..SAS_AF (6 fields) | ~2,706,840 | ~290,660 | 90.30% | AF from co-located shifted indels |
| gnomADe_AF | 2,953,405 | 44,099 | 98.53% | Same co-located AF issue |
| gnomADg_AF | 2,554,927 | 442,577 | 85.24% | Same co-located AF issue |
| MAX_AF | 2,553,486 | 444,018 | 85.19% | Derived from above |
| MAX_AF_POPS | 2,555,419 | 442,085 | 85.25% | Derived from above |
| CLIN_SIG | 2,991,851 | 5,653 | 99.81% | Co-located metadata for shifted indels |
| SOMATIC | 2,994,303 | 3,201 | 99.89% | Co-located metadata for shifted indels |
| PHENO | 2,982,158 | 15,346 | 99.49% | Co-located metadata for shifted indels |
| PUBMED | 2,989,573 | 7,931 | 99.74% | Co-located metadata for shifted indels |

### Prerequisites

- VEP parquet cache files for the target chromosome (e.g., `115_GRCh38_variation_1.parquet`)
- Context parquet files: transcripts, exons, translations, regulatory, motif (in the same directory)
- Input VCF (bgzipped + tabix-indexed)
- *(Optional, for generating golden standard)* Docker with Ensembl VEP 115 image and cache

### Positional arguments

```
cargo run --release --example annotate_vep_golden_bench -- \
  <source_vcf_gz>       # 1: Input VCF (bgzipped)
  <cache_source>        # 2: Variation cache path (parquet file or fjall directory)
  <backend>             # 3: "parquet" or "fjall"
  <sample_limit>        # 4: Number of variants to sample (0 = all)
  <vep_cache_dir>       # 5: Directory with golden VEP output / Docker VEP cache
  <local_copy_vcf_gz>   # 6: Local copy path for the input VCF
  <work_dir>            # 7: Output directory for reports and intermediate files
  <context_dir>         # 8: Directory with context parquet files (transcripts, exons, etc.)
  [--steps=ensembl,datafusion]  # Which steps to run (default: both)
  [--extended-probes]           # Enable interval-overlap fallback for shifted indels
  [--merged]                    # Use VEP --merged flag for merged Ensembl+RefSeq cache
```

### Running with a precomputed golden standard

When a golden VEP output file already exists in the work directory, the benchmark reuses it instead of re-running Docker. This is the recommended workflow to avoid the slow Docker VEP step.

**Step 1: Set up the golden file (one-time)**

The golden file must be named `HG002_chr22_{sample_limit}_vep115_golden.vcf` in the work directory. For `sample_limit=0` (all variants), the file is `HG002_chr22_0_vep115_golden.vcf`.

Create a symlink to a precomputed golden file:

```bash
mkdir -p /tmp/annotate_vep_golden_bench
ln -sf /path/to/HG002_chr1_full_vep115_golden.vcf \
  /tmp/annotate_vep_golden_bench/HG002_chr22_0_vep115_golden.vcf
```

**Step 2: Run the benchmark (both ensembl parse + datafusion annotation + comparison)**

```bash
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

The `--steps=ensembl,datafusion` flag is required for comparison. The `ensembl` step parses the golden VEP output (reusing the existing file), and `datafusion` runs our engine. Comparison only runs when both steps are enabled.

**Step 3: Run only our engine (skip golden parsing)**

To iterate on the annotation engine without re-parsing the golden file:

```bash
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /Users/mwiewior/research/data/vep \
  --steps=datafusion \
  --extended-probes
```

This skips golden parsing and comparison, only running `annotate_vep()` and printing row counts and timing.

### Output files

All output goes to `<work_dir>` (default: `/tmp/annotate_vep_golden_bench`):

| File | Description |
|------|-------------|
| `HG002_chr22_{N}.vcf` | Sampled VCF (first N variants, or all if N=0) |
| `HG002_chr22_{N}_norm.vcf` | Normalized VCF (multi-allelic decomposed) |
| `HG002_chr22_{N}_vep115_golden.vcf` | Golden VEP 115 output (or symlink to precomputed) |
| `HG002_chr22_{N}_comparison_report.txt` | Full comparison report with per-field accuracy |
| `HG002_chr22_{N}_discrepancies.txt` | Per-variant discrepancies (most_severe/term_set mismatches) |

### Benchmark input data

| File | Description |
|------|-------------|
| `vep-benchmark/data/HG002_chr22.vcf.gz` | HG002 chr22 variants (default small benchmark) |
| `vep-benchmark/data/HG002_chr22.vcf.gz.tbi` | Tabix index |
| `vep-benchmark/data/HG002_chr1.vcf.gz` | HG002 chr1 variants (full benchmark, 319k variants) |
| `vep-benchmark/data/HG002_chr1.vcf.gz.tbi` | Tabix index |

### Quick chr22 run (1000 variants, ~5s)

```bash
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_22.parquet \
  parquet \
  1000 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

### Generating a new golden standard (requires Docker)

If no golden file exists, the `ensembl` step runs VEP 115 in Docker:

```bash
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /path/to/vep/cache \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

This takes significantly longer (~30-60 min for chr1) as Docker VEP processes all 323k variants.

## Other benchmark examples

| Example | Description |
|---------|-------------|
| `bench_annotate` | End-to-end `annotate_vep()` performance benchmark |
| `lookup_parquet_bench` | `lookup_variants()` with parquet backend |
| `lookup_kv_bench` | `lookup_variants()` with fjall KV backend |
| `bench_fjall_gets` | Raw fjall KV store read throughput |

Run with:

```bash
cargo run --release --example <example_name> -- <args>
```
