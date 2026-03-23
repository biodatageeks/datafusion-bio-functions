# datafusion-bio-function-vep

VEP-oriented DataFusion functions and benchmark tooling.

## `annotate_vep()` Output Schema

The output of `annotate_vep()` consists of **VCF input columns** (pass-through) followed by **annotation columns** (89 columns total).

All 80 CSQ fields are exposed as **individual top-level Arrow columns** so Arrow-native consumers (vepyr, polars-bio) can access structured data without parsing the legacy pipe-delimited CSQ string. Per-transcript fields (Consequence, SYMBOL, Gene, SIFT, etc.) contain values from the **most-severe transcript** for each variant.

### VCF input columns (pass-through)

The first columns in the output are pass-through from the input VCF table. The schema depends on the VCF file and the table provider used. The following 5 columns are **required** by `annotate_vep()`:

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chrom` | `Utf8` | yes | Chromosome identifier |
| `start` | `Int64` | yes | Start position (0-based) |
| `end` | `Int64` | yes | End position (0-based, exclusive) |
| `ref` | `Utf8` | yes | Reference allele |
| `alt` | `Utf8` | yes | Alternate allele |

The `datafusion-bio-format-vcf` provider (used by polars-bio) decomposes VCF fields into typed columns. Common additional columns:

| Column | Type | Description |
|--------|------|-------------|
| `id` | `Utf8` | Variant ID (rsID) |
| `qual` | `Float64` | Variant quality score |
| `filter` | `Utf8` | Filter status (PASS, etc.) |
| INFO fields | varies | Decomposed into individual typed columns (e.g., `DP: Int32`) |
| Sample fields | varies | Decomposed genotype fields (e.g., `GT: Utf8`, `AD: List<Int32>`, `GQ: Int32`) |

All VCF columns are passed through unmodified. The exact schema depends on the VCF file's header.

### Annotation columns (89 total)

Appended after VCF columns. Column indices start at `vcf_field_count`.

**Meta columns (2):**

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `csq` | `Utf8` | Legacy pipe-delimited CSQ string. **Skipped (NULL) by default** — only computed when explicitly projected |
| 2 | `most_severe_consequence` | `Utf8` | SO term for the most severe consequence across all transcripts |

**Transcript-level columns (42) — from the most-severe transcript (`--everything` mode):**

| # | Column | Type | Source | Description |
|---|--------|------|--------|-------------|
| 3 | `Allele` | `Utf8` | engine | Variant allele used to calculate consequence |
| 4 | `Consequence` | `List<Utf8>` | engine | Consequence types (SO terms, can be multiple per transcript) |
| 5 | `IMPACT` | `Utf8` | engine | Impact rating: HIGH, MODERATE, LOW, MODIFIER |
| 6 | `SYMBOL` | `Utf8` | transcript | Gene symbol |
| 7 | `Gene` | `Utf8` | transcript | Ensembl gene ID |
| 8 | `Feature_type` | `Utf8` | engine | Transcript, RegulatoryFeature, or MotifFeature |
| 9 | `Feature` | `Utf8` | engine | Ensembl feature ID (transcript/regulatory/motif) |
| 10 | `BIOTYPE` | `Utf8` | transcript | Biotype of transcript (protein_coding, lncRNA, etc.) |
| 11 | `EXON` | `Utf8` | engine | Exon number / total exons |
| 12 | `INTRON` | `Utf8` | engine | Intron number / total introns |
| 13 | `HGVSc` | `Utf8` | engine | HGVS coding sequence notation |
| 14 | `HGVSp` | `Utf8` | engine | HGVS protein sequence notation |
| 15 | `cDNA_position` | `Utf8` | engine | Position in cDNA |
| 16 | `CDS_position` | `Utf8` | engine | Position in CDS |
| 17 | `Protein_position` | `Utf8` | engine | Position in protein |
| 18 | `Amino_acids` | `Utf8` | engine | Reference/variant amino acids |
| 19 | `Codons` | `Utf8` | engine | Reference/variant codons |
| 20 | `Existing_variation` | `List<Utf8>` | cache | Known variant identifiers (rsIDs, COSMIC, etc.) |
| 21 | `DISTANCE` | `Int64` | engine | Distance to transcript (upstream/downstream) |
| 22 | `STRAND` | `Int8` | transcript | Strand of feature (1 or -1) |
| 23 | `FLAGS` | `Utf8` | transcript | Transcript quality flags |
| 24 | `VARIANT_CLASS` | `Utf8` | engine | SO term variant class (SNV, deletion, insertion, etc.) |
| 25 | `SYMBOL_SOURCE` | `Utf8` | transcript | Source of gene symbol (HGNC, EntrezGene, etc.) |
| 26 | `HGNC_ID` | `Utf8` | transcript | HGNC gene identifier |
| 27 | `CANONICAL` | `Utf8` | transcript | Canonical transcript flag (YES or empty) |
| 28 | `MANE` | `Utf8` | transcript | MANE transcript type (MANE_Select or MANE_Plus_Clinical) |
| 29 | `MANE_SELECT` | `Utf8` | transcript | MANE Select transcript ID |
| 30 | `MANE_PLUS_CLINICAL` | `Utf8` | transcript | MANE Plus Clinical transcript ID |
| 31 | `TSL` | `Int8` | transcript | Transcript support level (1-5) |
| 32 | `APPRIS` | `Utf8` | transcript | APPRIS annotation (P1-P5, A1-A2) |
| 33 | `CCDS` | `Utf8` | transcript | CCDS identifier |
| 34 | `ENSP` | `Utf8` | transcript | Ensembl protein ID |
| 35 | `SWISSPROT` | `Utf8` | transcript | UniProtKB/Swiss-Prot ID |
| 36 | `TREMBL` | `Utf8` | transcript | UniProtKB/TrEMBL ID |
| 37 | `UNIPARC` | `Utf8` | transcript | UniParc ID |
| 38 | `UNIPROT_ISOFORM` | `Utf8` | transcript | UniProt isoform |
| 39 | `GENE_PHENO` | `Utf8` | transcript | Gene has known phenotype association |
| 40 | `SIFT` | `Utf8` | sift cache | SIFT prediction and score, e.g. `tolerated(0.23)` |
| 41 | `PolyPhen` | `Utf8` | sift cache | PolyPhen prediction and score, e.g. `benign(0.01)` |
| 42 | `DOMAINS` | `List<Utf8>` | engine | Overlapping protein domains (`source:id` pairs) |
| 43 | `miRNA` | `Utf8` | engine | miRNA secondary structure position |
| 44 | `HGVS_OFFSET` | `Int64` | engine | HGVS offset for indels (signed, strand-aware) |

**Population frequency columns (29) — resolved for the matching allele:**

| # | Column | Type | Source | Description |
|---|--------|------|--------|-------------|
| 45 | `AF` | `Float32` | cache | 1000 Genomes global allele frequency |
| 46 | `AFR_AF` | `Float32` | cache | 1000 Genomes African AF |
| 47 | `AMR_AF` | `Float32` | cache | 1000 Genomes Ad Mixed American AF |
| 48 | `EAS_AF` | `Float32` | cache | 1000 Genomes East Asian AF |
| 49 | `EUR_AF` | `Float32` | cache | 1000 Genomes European AF |
| 50 | `SAS_AF` | `Float32` | cache | 1000 Genomes South Asian AF |
| 51 | `gnomADe_AF` | `Float32` | cache | gnomAD exome global AF |
| 52 | `gnomADe_AFR_AF` | `Float32` | cache | gnomAD exome African/African American AF |
| 53 | `gnomADe_AMR_AF` | `Float32` | cache | gnomAD exome Latino/Admixed American AF |
| 54 | `gnomADe_ASJ_AF` | `Float32` | cache | gnomAD exome Ashkenazi Jewish AF |
| 55 | `gnomADe_EAS_AF` | `Float32` | cache | gnomAD exome East Asian AF |
| 56 | `gnomADe_FIN_AF` | `Float32` | cache | gnomAD exome Finnish AF |
| 57 | `gnomADe_MID_AF` | `Float32` | cache | gnomAD exome Middle Eastern AF |
| 58 | `gnomADe_NFE_AF` | `Float32` | cache | gnomAD exome Non-Finnish European AF |
| 59 | `gnomADe_REMAINING_AF` | `Float32` | cache | gnomAD exome remaining populations AF |
| 60 | `gnomADe_SAS_AF` | `Float32` | cache | gnomAD exome South Asian AF |
| 61 | `gnomADg_AF` | `Float32` | cache | gnomAD genome global AF |
| 62 | `gnomADg_AFR_AF` | `Float32` | cache | gnomAD genome African/African American AF |
| 63 | `gnomADg_AMI_AF` | `Float32` | cache | gnomAD genome Amish AF |
| 64 | `gnomADg_AMR_AF` | `Float32` | cache | gnomAD genome Latino/Admixed American AF |
| 65 | `gnomADg_ASJ_AF` | `Float32` | cache | gnomAD genome Ashkenazi Jewish AF |
| 66 | `gnomADg_EAS_AF` | `Float32` | cache | gnomAD genome East Asian AF |
| 67 | `gnomADg_FIN_AF` | `Float32` | cache | gnomAD genome Finnish AF |
| 68 | `gnomADg_MID_AF` | `Float32` | cache | gnomAD genome Middle Eastern AF |
| 69 | `gnomADg_NFE_AF` | `Float32` | cache | gnomAD genome Non-Finnish European AF |
| 70 | `gnomADg_REMAINING_AF` | `Float32` | cache | gnomAD genome remaining populations AF |
| 71 | `gnomADg_SAS_AF` | `Float32` | cache | gnomAD genome South Asian AF |
| 72 | `MAX_AF` | `Float32` | computed | Maximum AF across all populations |
| 73 | `MAX_AF_POPS` | `Utf8` | computed | Population(s) with maximum AF (`&`-separated if tied) |

**Variant-level annotation columns (9) — from variation cache lookup:**

| # | Column | Type | Source | Description |
|---|--------|------|--------|-------------|
| 74 | `CLIN_SIG` | `List<Utf8>` | cache | ClinVar clinical significance (multi-valued) |
| 75 | `SOMATIC` | `Utf8` | cache | Somatic variant flag |
| 76 | `PHENO` | `Utf8` | cache | Phenotype/disease association flag |
| 77 | `PUBMED` | `List<Utf8>` | cache | PubMed IDs |
| 78 | `MOTIF_NAME` | `Utf8` | engine | Motif feature name |
| 79 | `MOTIF_POS` | `Utf8` | engine | Position in motif |
| 80 | `HIGH_INF_POS` | `Utf8` | engine | High information position in motif |
| 81 | `MOTIF_SCORE_CHANGE` | `Float32` | engine | Change in motif score |
| 82 | `TRANSCRIPTION_FACTORS` | `List<Utf8>` | engine | Transcription factors binding the motif |

**Cache-only columns (7) — not in CSQ, from variation cache:**

| # | Column | Type | Source | Description |
|---|--------|------|--------|-------------|
| 83 | `clin_sig_allele` | `List<Utf8>` | cache | Clinical significance per allele (`;`-separated) |
| 84 | `clinical_impact` | `Utf8` | cache | Clinical impact assessment |
| 85 | `minor_allele` | `Utf8` | cache | Minor allele |
| 86 | `minor_allele_freq` | `Float32` | cache | Minor allele frequency |
| 87 | `clinvar_ids` | `List<Utf8>` | cache | ClinVar accession IDs |
| 88 | `cosmic_ids` | `List<Utf8>` | cache | COSMIC IDs |
| 89 | `dbsnp_ids` | `List<Utf8>` | cache | dbSNP IDs |

### Notes

- **Per-transcript fields** (columns 3-44): One variant may overlap multiple transcripts. Top-level columns contain values from the transcript with the **most severe consequence**. The legacy `csq` column (if projected) contains all transcript entries.
- **AF resolution**: Cache stores `allele:frequency` pairs (e.g., `"T:0.9301"`, multi-allelic: `"A:0.006,G:0.994"`). Top-level AF columns contain the **resolved Float32 frequency** for the matching allele via `extract_af_for_allele()`.
- **Float32 precision**: ~7.2 significant digits covers all observed AF values (max 7 decimal places in Ensembl 115).
- **`csq` column**: Skipped by default. Only assembled when explicitly projected. Legacy VEP pipe-delimited format.
- **`--everything` mode**: Enables all 80 CSQ fields. Without it, 6 fields are absent: MANE, APPRIS, SIFT, PolyPhen, DOMAINS, miRNA, HGVS_OFFSET.

### VEP flags (`options_json`)

| Flag | Effect | Implied by `everything` |
|------|--------|------------------------|
| `everything` | Enables all flags, all 80 CSQ fields populated | - |
| `check_existing` | Variation lookup in cache | Yes (also implied by any AF flag) |
| `af` | Global AF | Yes |
| `af_1kg` | 1000 Genomes AF sub-populations | Yes |
| `af_gnomade` | gnomAD exome sub-populations | Yes |
| `af_gnomadg` | gnomAD genome sub-populations | Yes |
| `max_af` | MAX_AF / MAX_AF_POPS | Yes |
| `pubmed` | PubMed IDs | Yes |
| `hgvs` | HGVSc / HGVSp notation | Yes |

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

Output filenames are derived from the source VCF stem. For `HG002_chr1.vcf.gz` with `sample_limit=0`, the golden file must be named `HG002_chr1_0_vep115_golden.vcf` in the work directory.

Create a symlink to a precomputed golden file:

```bash
mkdir -p /tmp/annotate_vep_golden_bench
ln -sf /path/to/HG002_chr1_full_vep115_golden.vcf \
  /tmp/annotate_vep_golden_bench/HG002_chr1_0_vep115_golden.vcf
```

**Step 2: Run the benchmark (both ensembl parse + datafusion annotation + comparison)**

```bash
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr1.vcf.gz \
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
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /Users/mwiewior/research/data/vep \
  --steps=datafusion \
  --extended-probes
```

This skips golden parsing and comparison, only running `annotate_vep()` and printing row counts and timing.

### Output files

All output goes to `<work_dir>` (default: `/tmp/annotate_vep_golden_bench`). Filenames are derived from the source VCF stem (e.g., `HG002_chr1` from `HG002_chr1.vcf.gz`):

| File | Description |
|------|-------------|
| `{stem}_{N}.vcf` | Sampled VCF (first N variants, or all if N=0) |
| `{stem}_{N}_norm.vcf` | Normalized VCF (multi-allelic decomposed) |
| `{stem}_{N}_vep115_golden.vcf` | Golden VEP 115 output (or symlink to precomputed) |
| `{stem}_{N}_comparison_report.txt` | Full comparison report with per-field accuracy |
| `{stem}_{N}_discrepancies.txt` | Per-variant discrepancies (most_severe/term_set mismatches) |

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
  --steps=ensembl,datafusion
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
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

This takes significantly longer (~30-60 min for chr1) as Docker VEP processes all 323k variants.

## E2E Annotation Benchmark (`bench_annotate_vcf`)

Full pipeline benchmark: read VCF → annotate with transcript engine → write annotated VCF. Measures end-to-end throughput including I/O.

```bash
cargo run --release --features kv-cache --example bench_annotate_vcf -- \
  --input <vcf_path> \
  --cache <cache_dir> \
  --output <output.vcf> \
  [--backend parquet|fjall] \
  [--everything] \
  [--extended-probes] \
  [--reference-fasta <path>] \
  [--compression none|gzip|bgzf] \
  [--limit <n>]
```

### Examples

```bash
# Parquet backend, --everything, 1000 variants from golden fixtures
cargo run --release --example bench_annotate_vcf -- \
  --input vep-benchmark/data/golden/input_1000.vcf \
  --cache vep-benchmark/data/golden/cache \
  --output /tmp/annotated.vcf \
  --everything --extended-probes \
  --reference-fasta vep-benchmark/data/golden/reference_chr1.fa

# Fjall backend, full chr1 (323K variants), gzip output
cargo run --release --features kv-cache --example bench_annotate_vcf -- \
  --input vep-benchmark/data/HG002_chr1.vcf.gz \
  --cache /data/vep/wgs/fjall/115_GRCh38_vep \
  --output /tmp/annotated_chr1.vcf.gz \
  --backend fjall --everything --extended-probes \
  --reference-fasta /data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --compression gzip

# Quick smoke test with LIMIT
cargo run --release --example bench_annotate_vcf -- \
  --input vep-benchmark/data/golden/input_1000.vcf \
  --cache vep-benchmark/data/golden/cache \
  --output /tmp/quick.vcf \
  --limit 100
```

### Output

```
=== bench_annotate_vcf ===
  input:      vep-benchmark/data/golden/input_1000.vcf
  cache:      vep-benchmark/data/golden/cache
  backend:    parquet
  everything: true
  ...

[0.7s] Annotation + VCF write complete

=== Results ===
  rows:       1000
  time:       0.70s
  throughput: 1419 variants/s
  output:     13.1 MB
  total:      0.71s (including VCF read)
```

## Other benchmark examples

| Example | Description |
|---------|-------------|
| `bench_annotate_vcf` | E2E pipeline: VCF in → annotate → VCF out (parquet or fjall) |
| `bench_annotate` | `annotate_vep()` Arrow-level performance (fjall only) |
| `annotate_vep_golden_bench` | Golden comparison: per-field CSQ accuracy vs Ensembl VEP 115 |
| `lookup_parquet_bench` | `lookup_variants()` with parquet backend |
| `lookup_kv_bench` | `lookup_variants()` with fjall KV backend |
| `bench_fjall_gets` | Raw fjall KV store read throughput |

Run with:

```bash
cargo run --release --example <example_name> -- <args>
```
