# datafusion-bio-function-vep

VEP-oriented DataFusion functions and benchmark tooling.

## `annotate_vep()` Output Schema

The output of `annotate_vep()` consists of **VCF input columns** (pass-through) followed by **annotation columns**.

### Annotation columns

Appended after VCF columns. Column indices start at `vcf_field_count`.

| # | Column | Type | Nullable | Description |
|---|--------|------|----------|-------------|
| 1 | `csq` | `Utf8` | yes | Pipe-delimited CSQ string (VEP legacy format). **Skipped (NULL) by default** — only computed when explicitly included in projection |
| 2 | `most_severe_consequence` | `Utf8` | yes | SO term for the most severe consequence across all transcripts |
| 3 | `variation_name` | `Utf8` | yes | dbSNP rsID or known variant identifier |
| 4 | `clin_sig` | `List<Utf8>` | yes | ClinVar clinical significance (multi-valued, `,`-separated) |
| 5 | `clin_sig_allele` | `List<Utf8>` | yes | Clinical significance per allele (`;`-separated) |
| 6 | `clinical_impact` | `Utf8` | yes | Clinical impact assessment |
| 7 | `phenotype_or_disease` | `Utf8` | yes | Associated phenotype/disease flag |
| 8 | `pubmed` | `List<Utf8>` | yes | PubMed IDs (multi-valued, `,`-separated) |
| 9 | `somatic` | `Utf8` | yes | Somatic variant flag |
| 10 | `minor_allele` | `Utf8` | yes | Minor allele |
| 11 | `minor_allele_freq` | `Utf8` | yes | Minor allele frequency |
| 12 | `AF` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | 1000 Genomes global allele frequency |
| 13 | `AFR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | 1000 Genomes African population AF |
| 14 | `AMR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | 1000 Genomes Ad Mixed American AF |
| 15 | `EAS` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | 1000 Genomes East Asian AF |
| 16 | `EUR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | 1000 Genomes European AF |
| 17 | `SAS` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | 1000 Genomes South Asian AF |
| 18 | `gnomADe` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome global AF |
| 19 | `gnomADe_AFR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome African/African American AF |
| 20 | `gnomADe_AMR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome Latino/Admixed American AF |
| 21 | `gnomADe_ASJ` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome Ashkenazi Jewish AF |
| 22 | `gnomADe_EAS` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome East Asian AF |
| 23 | `gnomADe_FIN` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome Finnish AF |
| 24 | `gnomADe_NFE` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome Non-Finnish European AF |
| 25 | `gnomADe_SAS` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome South Asian AF |
| 26 | `gnomADe_MID` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome Middle Eastern AF |
| 27 | `gnomADe_REMAINING` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD exome remaining populations AF |
| 28 | `gnomADg` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome global AF |
| 29 | `gnomADg_AFR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome African/African American AF |
| 30 | `gnomADg_AMI` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome Amish AF |
| 31 | `gnomADg_AMR` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome Latino/Admixed American AF |
| 32 | `gnomADg_ASJ` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome Ashkenazi Jewish AF |
| 33 | `gnomADg_EAS` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome East Asian AF |
| 34 | `gnomADg_FIN` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome Finnish AF |
| 35 | `gnomADg_MID` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome Middle Eastern AF |
| 36 | `gnomADg_NFE` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome Non-Finnish European AF |
| 37 | `gnomADg_SAS` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome South Asian AF |
| 38 | `gnomADg_REMAINING` | `List<Struct<allele: Utf8, freq: Float32>>` | yes | gnomAD genome remaining populations AF |
| 39 | `clinvar_ids` | `List<Utf8>` | yes | ClinVar accession IDs (`,`-separated) |
| 40 | `cosmic_ids` | `List<Utf8>` | yes | COSMIC IDs (`,`-separated) |
| 41 | `dbsnp_ids` | `List<Utf8>` | yes | dbSNP IDs (`,`-separated) |

**Notes:**
- AF columns in the source parquet cache store `allele:frequency` pairs (e.g., `"T:0.9301"`, multi-allelic: `"A:0.006,G:0.994"`). The output parses these into structured `List<Struct<allele, freq>>` arrays.
- Float32 precision (~7.2 significant digits) covers all observed AF values (max 7 decimal places in Ensembl 115).
- The `csq` column is skipped by default for Arrow-native consumers. It is only assembled when explicitly included in the query projection.

### CSQ sub-fields: default mode (74 fields)

When `csq` is projected, it contains pipe-delimited fields per transcript. One `csq` value may contain multiple comma-separated transcript entries.

| # | Field | Description |
|---|-------|-------------|
| 1 | Allele | Variant allele used to calculate consequence |
| 2 | Consequence | Consequence type (SO term) |
| 3 | IMPACT | Impact rating (HIGH, MODERATE, LOW, MODIFIER) |
| 4 | SYMBOL | Gene symbol |
| 5 | Gene | Ensembl gene ID |
| 6 | Feature_type | Feature type (Transcript, RegulatoryFeature, MotifFeature) |
| 7 | Feature | Ensembl feature ID |
| 8 | BIOTYPE | Biotype of transcript |
| 9 | EXON | Exon number / total exons |
| 10 | INTRON | Intron number / total introns |
| 11 | HGVSc | HGVS coding sequence notation |
| 12 | HGVSp | HGVS protein sequence notation |
| 13 | cDNA_position | Position in cDNA |
| 14 | CDS_position | Position in CDS |
| 15 | Protein_position | Position in protein |
| 16 | Amino_acids | Reference/variant amino acids |
| 17 | Codons | Reference/variant codons |
| 18 | Existing_variation | Known variant identifiers |
| 19 | DISTANCE | Distance to transcript |
| 20 | STRAND | Strand of feature (+1/-1) |
| 21 | FLAGS | Transcript quality flags |
| 22 | SYMBOL_SOURCE | Source of gene symbol (HGNC, EntrezGene, etc.) |
| 23 | HGNC_ID | HGNC gene identifier |
| 24 | MOTIF_NAME | Motif feature name |
| 25 | MOTIF_POS | Position in motif |
| 26 | HIGH_INF_POS | High information position in motif |
| 27 | MOTIF_SCORE_CHANGE | Change in motif score |
| 28 | TRANSCRIPTION_FACTORS | Transcription factors binding the motif |
| 29 | SOURCE | Source of transcript (Ensembl, RefSeq — only in --merged) |
| 30 | VARIANT_CLASS | SO term variant class |
| 31 | CANONICAL | Canonical transcript flag |
| 32 | TSL | Transcript support level |
| 33 | MANE_SELECT | MANE Select transcript ID |
| 34 | MANE_PLUS_CLINICAL | MANE Plus Clinical transcript ID |
| 35 | ENSP | Ensembl protein ID |
| 36 | GENE_PHENO | Gene has known phenotype association |
| 37 | CCDS | CCDS identifier |
| 38 | SWISSPROT | UniProtKB/Swiss-Prot ID |
| 39 | TREMBL | UniProtKB/TrEMBL ID |
| 40 | UNIPARC | UniParc ID |
| 41 | UNIPROT_ISOFORM | UniProt isoform |
| 42 | AF | 1000 Genomes global allele frequency |
| 43 | AFR_AF | 1000 Genomes African AF |
| 44 | AMR_AF | 1000 Genomes American AF |
| 45 | EAS_AF | 1000 Genomes East Asian AF |
| 46 | EUR_AF | 1000 Genomes European AF |
| 47 | SAS_AF | 1000 Genomes South Asian AF |
| 48 | gnomADe_AF | gnomAD exome global AF |
| 49 | gnomADe_AFR | gnomAD exome African AF |
| 50 | gnomADe_AMR | gnomAD exome American AF |
| 51 | gnomADe_ASJ | gnomAD exome Ashkenazi Jewish AF |
| 52 | gnomADe_EAS | gnomAD exome East Asian AF |
| 53 | gnomADe_FIN | gnomAD exome Finnish AF |
| 54 | gnomADe_MID | gnomAD exome Middle Eastern AF |
| 55 | gnomADe_NFE | gnomAD exome Non-Finnish European AF |
| 56 | gnomADe_REMAINING | gnomAD exome remaining pops AF |
| 57 | gnomADe_SAS | gnomAD exome South Asian AF |
| 58 | gnomADg_AF | gnomAD genome global AF |
| 59 | gnomADg_AFR | gnomAD genome African AF |
| 60 | gnomADg_AMI | gnomAD genome Amish AF |
| 61 | gnomADg_AMR | gnomAD genome American AF |
| 62 | gnomADg_ASJ | gnomAD genome Ashkenazi Jewish AF |
| 63 | gnomADg_EAS | gnomAD genome East Asian AF |
| 64 | gnomADg_FIN | gnomAD genome Finnish AF |
| 65 | gnomADg_MID | gnomAD genome Middle Eastern AF |
| 66 | gnomADg_NFE | gnomAD genome Non-Finnish European AF |
| 67 | gnomADg_REMAINING | gnomAD genome remaining pops AF |
| 68 | gnomADg_SAS | gnomAD genome South Asian AF |
| 69 | MAX_AF | Maximum AF across all populations |
| 70 | MAX_AF_POPS | Population(s) with maximum AF |
| 71 | CLIN_SIG | ClinVar clinical significance |
| 72 | SOMATIC | Somatic variant flag |
| 73 | PHENO | Phenotype/disease association flag |
| 74 | PUBMED | PubMed IDs |

### CSQ sub-fields: `--everything` mode (80 fields)

Uses 80-field layout when `options_json` contains `"everything": true`. Key differences from default: VARIANT_CLASS moved after FLAGS, SOURCE removed, MANE/APPRIS/SIFT/PolyPhen/DOMAINS/miRNA/HGVS_OFFSET added, gnomAD sub-pops gain `_AF` suffix, MOTIF fields moved to end.

| # | Field | Description |
|---|-------|-------------|
| 1 | Allele | Variant allele used to calculate consequence |
| 2 | Consequence | Consequence type (SO term) |
| 3 | IMPACT | Impact rating (HIGH, MODERATE, LOW, MODIFIER) |
| 4 | SYMBOL | Gene symbol |
| 5 | Gene | Ensembl gene ID |
| 6 | Feature_type | Feature type |
| 7 | Feature | Ensembl feature ID |
| 8 | BIOTYPE | Biotype of transcript |
| 9 | EXON | Exon number / total exons |
| 10 | INTRON | Intron number / total introns |
| 11 | HGVSc | HGVS coding sequence notation |
| 12 | HGVSp | HGVS protein sequence notation |
| 13 | cDNA_position | Position in cDNA |
| 14 | CDS_position | Position in CDS |
| 15 | Protein_position | Position in protein |
| 16 | Amino_acids | Reference/variant amino acids |
| 17 | Codons | Reference/variant codons |
| 18 | Existing_variation | Known variant identifiers |
| 19 | DISTANCE | Distance to transcript |
| 20 | STRAND | Strand of feature |
| 21 | FLAGS | Transcript quality flags |
| 22 | VARIANT_CLASS | SO term variant class |
| 23 | SYMBOL_SOURCE | Source of gene symbol |
| 24 | HGNC_ID | HGNC gene identifier |
| 25 | CANONICAL | Canonical transcript flag |
| 26 | MANE | MANE transcript ID (generic) |
| 27 | MANE_SELECT | MANE Select transcript ID |
| 28 | MANE_PLUS_CLINICAL | MANE Plus Clinical transcript ID |
| 29 | TSL | Transcript support level |
| 30 | APPRIS | APPRIS annotation |
| 31 | CCDS | CCDS identifier |
| 32 | ENSP | Ensembl protein ID |
| 33 | SWISSPROT | UniProtKB/Swiss-Prot ID |
| 34 | TREMBL | UniProtKB/TrEMBL ID |
| 35 | UNIPARC | UniParc ID |
| 36 | UNIPROT_ISOFORM | UniProt isoform |
| 37 | GENE_PHENO | Gene has known phenotype association |
| 38 | SIFT | SIFT prediction and score |
| 39 | PolyPhen | PolyPhen prediction and score |
| 40 | DOMAINS | Overlapping protein domains |
| 41 | miRNA | miRNA secondary structure position |
| 42 | HGVS_OFFSET | HGVS offset for indels |
| 43 | AF | 1000 Genomes global AF |
| 44 | AFR_AF | 1000 Genomes African AF |
| 45 | AMR_AF | 1000 Genomes American AF |
| 46 | EAS_AF | 1000 Genomes East Asian AF |
| 47 | EUR_AF | 1000 Genomes European AF |
| 48 | SAS_AF | 1000 Genomes South Asian AF |
| 49 | gnomADe_AF | gnomAD exome global AF |
| 50 | gnomADe_AFR_AF | gnomAD exome African AF |
| 51 | gnomADe_AMR_AF | gnomAD exome American AF |
| 52 | gnomADe_ASJ_AF | gnomAD exome Ashkenazi Jewish AF |
| 53 | gnomADe_EAS_AF | gnomAD exome East Asian AF |
| 54 | gnomADe_FIN_AF | gnomAD exome Finnish AF |
| 55 | gnomADe_MID_AF | gnomAD exome Middle Eastern AF |
| 56 | gnomADe_NFE_AF | gnomAD exome Non-Finnish European AF |
| 57 | gnomADe_REMAINING_AF | gnomAD exome remaining pops AF |
| 58 | gnomADe_SAS_AF | gnomAD exome South Asian AF |
| 59 | gnomADg_AF | gnomAD genome global AF |
| 60 | gnomADg_AFR_AF | gnomAD genome African AF |
| 61 | gnomADg_AMI_AF | gnomAD genome Amish AF |
| 62 | gnomADg_AMR_AF | gnomAD genome American AF |
| 63 | gnomADg_ASJ_AF | gnomAD genome Ashkenazi Jewish AF |
| 64 | gnomADg_EAS_AF | gnomAD genome East Asian AF |
| 65 | gnomADg_FIN_AF | gnomAD genome Finnish AF |
| 66 | gnomADg_MID_AF | gnomAD genome Middle Eastern AF |
| 67 | gnomADg_NFE_AF | gnomAD genome Non-Finnish European AF |
| 68 | gnomADg_REMAINING_AF | gnomAD genome remaining pops AF |
| 69 | gnomADg_SAS_AF | gnomAD genome South Asian AF |
| 70 | MAX_AF | Maximum AF across all populations |
| 71 | MAX_AF_POPS | Population(s) with maximum AF |
| 72 | CLIN_SIG | ClinVar clinical significance |
| 73 | SOMATIC | Somatic variant flag |
| 74 | PHENO | Phenotype/disease association flag |
| 75 | PUBMED | PubMed IDs |
| 76 | MOTIF_NAME | Motif feature name |
| 77 | MOTIF_POS | Position in motif |
| 78 | HIGH_INF_POS | High information position in motif |
| 79 | MOTIF_SCORE_CHANGE | Change in motif score |
| 80 | TRANSCRIPTION_FACTORS | Transcription factors binding the motif |

### VEP flags (`options_json`)

| Flag | Effect | Implied by `everything` |
|------|--------|------------------------|
| `everything` | Enables all flags, uses 80-field CSQ layout | - |
| `check_existing` | Variation lookup in cache | Yes (also implied by any AF flag) |
| `af` | Global AF in CSQ | Yes |
| `af_1kg` | 1000 Genomes AF sub-populations | Yes |
| `af_gnomade` | gnomAD exome sub-populations | Yes |
| `af_gnomadg` | gnomAD genome sub-populations | Yes |
| `max_af` | MAX_AF / MAX_AF_POPS fields | Yes |
| `pubmed` | PubMed IDs field | Yes |
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
