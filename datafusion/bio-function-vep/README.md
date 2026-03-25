# datafusion-bio-function-vep

Apache DataFusion UDFs for variant effect prediction (VEP) annotation. Provides scalar functions for allele manipulation and table functions for variant lookup and full VEP-style annotation against Ensembl variation caches.

## Functions

### Scalar UDFs

| Function | Signature | Description |
|----------|-----------|-------------|
| `match_allele` | `(ref, alt, allele_string) → Boolean` | Check if a VCF variant matches a cache `allele_string` |
| `match_allele_relaxed` | `(ref, alt, allele_string) → Boolean` | Relaxed allele matching (fallback for shifted indels) |
| `vep_allele` | `(ref, alt) → Utf8` | Convert VCF REF/ALT to VEP allele format |
| `vep_norm_start` | `(pos, ref, alt) → Int64` | Normalize variant start position (1-based, Ensembl-style) |
| `vep_norm_end` | `(pos, ref, alt) → Int64` | Normalize variant end position |

### Table Functions

#### `lookup_variants`

Low-level known-variant lookup via left interval join against a variation cache.

```sql
SELECT * FROM lookup_variants(
  'my_vcf_table',
  'my_variation_cache',
  'variation_name,clin_sig,AF',   -- columns to return (optional)
  'exact',                         -- match mode (optional)
  'true'                           -- extended probes (optional)
)
```

#### `annotate_vep`

Full VEP-style annotation: consequence prediction, population frequencies, clinical significance, HGVS notation, SIFT/PolyPhen predictions.

```sql
SELECT * FROM annotate_vep(
  'my_vcf_table',
  '/path/to/cache',
  'parquet',                                          -- backend: parquet | fjall
  '{"everything":true,"extended_probes":true,
    "reference_fasta_path":"/path/to/reference.fa"}'  -- options (optional)
)
```

**Options JSON:**

| Flag | Effect | Implied by `everything` |
|------|--------|------------------------|
| `everything` | Enable all annotation fields (80 CSQ columns) | — |
| `check_existing` | Variation lookup in cache | Yes (also implied by any AF flag) |
| `af` | Global allele frequency | Yes |
| `af_1kg` | 1000 Genomes sub-population AFs | Yes |
| `af_gnomade` | gnomAD exome sub-population AFs | Yes |
| `af_gnomadg` | gnomAD genome sub-population AFs | Yes |
| `max_af` | MAX_AF / MAX_AF_POPS | Yes |
| `pubmed` | PubMed IDs | Yes |
| `hgvs` | HGVSc / HGVSp notation | Yes |
| `extended_probes` | Interval-overlap fallback for shifted indels | No |
| `reference_fasta_path` | Indexed FASTA for HGVS genomic shifting | No |
| `partitioned` | Use per-chromosome partitioned cache | No |
| `merged` | Use VEP `--merged` Ensembl+RefSeq cache | No |

## Output Schema

The output of `annotate_vep()` consists of **VCF pass-through columns** followed by **annotation columns** (up to 89 total).

All CSQ fields are exposed as individual top-level Arrow columns. Per-transcript fields contain values from the **most-severe transcript** for each variant.

### VCF input columns (pass-through)

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chrom` | `Utf8` | yes | Chromosome identifier |
| `start` | `Int64` | yes | Start position (0-based) |
| `end` | `Int64` | yes | End position (0-based, exclusive) |
| `ref` | `Utf8` | yes | Reference allele |
| `alt` | `Utf8` | yes | Alternate allele |

Additional VCF columns (`id`, `qual`, `filter`, decomposed INFO/sample fields) are passed through unmodified.

### Meta columns (2)

| Column | Type | Description |
|--------|------|-------------|
| `csq` | `Utf8` | Legacy pipe-delimited CSQ string (NULL unless explicitly projected) |
| `most_severe_consequence` | `Utf8` | SO term for the most severe consequence across all transcripts |

### Transcript-level columns (42)

Values from the transcript with the most severe consequence. Enabled with `everything: true`.

| Column | Type | Description |
|--------|------|-------------|
| `Allele` | `Utf8` | Variant allele used to calculate consequence |
| `Consequence` | `List<Utf8>` | Consequence types (SO terms) |
| `IMPACT` | `Utf8` | Impact rating: HIGH, MODERATE, LOW, MODIFIER |
| `SYMBOL` | `Utf8` | Gene symbol |
| `Gene` | `Utf8` | Ensembl gene ID |
| `Feature_type` | `Utf8` | Transcript, RegulatoryFeature, or MotifFeature |
| `Feature` | `Utf8` | Ensembl feature ID |
| `BIOTYPE` | `Utf8` | Biotype (protein_coding, lncRNA, etc.) |
| `EXON` | `Utf8` | Exon number / total exons |
| `INTRON` | `Utf8` | Intron number / total introns |
| `HGVSc` | `Utf8` | HGVS coding sequence notation |
| `HGVSp` | `Utf8` | HGVS protein sequence notation |
| `cDNA_position` | `Utf8` | Position in cDNA |
| `CDS_position` | `Utf8` | Position in CDS |
| `Protein_position` | `Utf8` | Position in protein |
| `Amino_acids` | `Utf8` | Reference/variant amino acids |
| `Codons` | `Utf8` | Reference/variant codons |
| `Existing_variation` | `List<Utf8>` | Known variant identifiers (rsIDs, COSMIC, etc.) |
| `DISTANCE` | `Int64` | Distance to transcript |
| `STRAND` | `Int8` | Strand of feature (1 or -1) |
| `FLAGS` | `Utf8` | Transcript quality flags |
| `VARIANT_CLASS` | `Utf8` | SO variant class (SNV, deletion, insertion, etc.) |
| `SYMBOL_SOURCE` | `Utf8` | Source of gene symbol (HGNC, EntrezGene, etc.) |
| `HGNC_ID` | `Utf8` | HGNC gene identifier |
| `CANONICAL` | `Utf8` | Canonical transcript flag |
| `MANE` | `Utf8` | MANE transcript type |
| `MANE_SELECT` | `Utf8` | MANE Select transcript ID |
| `MANE_PLUS_CLINICAL` | `Utf8` | MANE Plus Clinical transcript ID |
| `TSL` | `Int8` | Transcript support level (1-5) |
| `APPRIS` | `Utf8` | APPRIS annotation |
| `CCDS` | `Utf8` | CCDS identifier |
| `ENSP` | `Utf8` | Ensembl protein ID |
| `SWISSPROT` | `Utf8` | UniProtKB/Swiss-Prot ID |
| `TREMBL` | `Utf8` | UniProtKB/TrEMBL ID |
| `UNIPARC` | `Utf8` | UniParc ID |
| `UNIPROT_ISOFORM` | `Utf8` | UniProt isoform |
| `GENE_PHENO` | `Utf8` | Gene has known phenotype association |
| `SIFT` | `Utf8` | SIFT prediction and score |
| `PolyPhen` | `Utf8` | PolyPhen prediction and score |
| `DOMAINS` | `List<Utf8>` | Overlapping protein domains |
| `miRNA` | `Utf8` | miRNA secondary structure position |
| `HGVS_OFFSET` | `Int64` | HGVS offset for indels |

### Population frequency columns (29)

Resolved for the matching allele. Stored as `Float32` (~7.2 significant digits).

| Column | Source | Description |
|--------|--------|-------------|
| `AF` | 1000 Genomes | Global allele frequency |
| `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF` | 1000 Genomes | Sub-population AFs |
| `gnomADe_AF` | gnomAD exome | Global exome AF |
| `gnomADe_{AFR,AMR,ASJ,EAS,FIN,MID,NFE,REMAINING,SAS}_AF` | gnomAD exome | Sub-population exome AFs |
| `gnomADg_AF` | gnomAD genome | Global genome AF |
| `gnomADg_{AFR,AMI,AMR,ASJ,EAS,FIN,MID,NFE,REMAINING,SAS}_AF` | gnomAD genome | Sub-population genome AFs |
| `MAX_AF` | computed | Maximum AF across all populations |
| `MAX_AF_POPS` | computed | Population(s) with maximum AF |

### Variant-level annotation columns (9)

| Column | Type | Description |
|--------|------|-------------|
| `CLIN_SIG` | `List<Utf8>` | ClinVar clinical significance |
| `SOMATIC` | `Utf8` | Somatic variant flag |
| `PHENO` | `Utf8` | Phenotype/disease association flag |
| `PUBMED` | `List<Utf8>` | PubMed IDs |
| `MOTIF_NAME` | `Utf8` | Motif feature name |
| `MOTIF_POS` | `Utf8` | Position in motif |
| `HIGH_INF_POS` | `Utf8` | High information position in motif |
| `MOTIF_SCORE_CHANGE` | `Float32` | Change in motif score |
| `TRANSCRIPTION_FACTORS` | `List<Utf8>` | Transcription factors binding the motif |

### Cache-only columns (7)

Additional fields from the variation cache, not present in the CSQ string.

| Column | Type | Description |
|--------|------|-------------|
| `clin_sig_allele` | `List<Utf8>` | Clinical significance per allele |
| `clinical_impact` | `Utf8` | Clinical impact assessment |
| `minor_allele` | `Utf8` | Minor allele |
| `minor_allele_freq` | `Float32` | Minor allele frequency |
| `clinvar_ids` | `List<Utf8>` | ClinVar accession IDs |
| `cosmic_ids` | `List<Utf8>` | COSMIC IDs |
| `dbsnp_ids` | `List<Utf8>` | dbSNP IDs |

## Cache Backends

### Parquet (default)

Ensembl variation cache converted to Parquet format. Supports per-chromosome partitioning.

### Fjall (KV store)

LSM-tree backend using [fjall](https://crates.io/crates/fjall) with zstd compression. Requires the `kv-cache` feature (enabled by default).

Configuration via `AnnotationConfig` (session-level):

| Option | Default | Description |
|--------|---------|-------------|
| `cache_size_mb` | 1024 | Fjall block cache size |
| `zstd_level` | 3 | Compression level (1-19) |
| `dict_size_kb` | 112 | Dictionary size for compression |

## Prerequisites

- Ensembl VEP parquet cache files for the target chromosome(s)
- Context parquet files: transcripts, exons, translations, regulatory, motif features
- Input VCF (bgzipped + tabix-indexed for large files)
- *(Optional)* Indexed reference FASTA for HGVS genomic shifting
- *(Optional)* Docker with Ensembl VEP 115 image for golden benchmark generation

## Benchmarks

### Golden benchmark (`annotate_vep_golden_bench`)

Compares `annotate_vep()` output against Ensembl VEP 115 on real whole-genome data (HG002/NA24385, GIAB). Validates per-field CSQ accuracy across 74 fields.

```bash
# Full chr1 benchmark (323K variants, ~84s)
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /path/to/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /path/to/vep_cache \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /path/to/context_parquets \
  --steps=ensembl,datafusion \
  --extended-probes

# Quick chr22 run (1000 variants, ~5s)
cargo run --release --example annotate_vep_golden_bench -- \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /path/to/115_GRCh38_variation_22.parquet \
  parquet \
  1000 \
  /path/to/vep_cache \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /tmp/annotate_vep_golden_bench \
  /path/to/context_parquets \
  --steps=ensembl,datafusion
```

**Positional arguments:**

| # | Argument | Description |
|---|----------|-------------|
| 1 | `source_vcf_gz` | Input VCF (bgzipped) |
| 2 | `cache_source` | Variation cache path (parquet file or fjall directory) |
| 3 | `backend` | `parquet` or `fjall` |
| 4 | `sample_limit` | Number of variants to sample (0 = all) |
| 5 | `vep_cache_dir` | Directory with golden VEP output / Docker VEP cache |
| 6 | `local_copy_vcf_gz` | Local copy path for the input VCF |
| 7 | `work_dir` | Output directory for reports and intermediate files |
| 8 | `context_dir` | Directory with context parquet files |

**Optional flags:** `--steps=ensembl,datafusion`, `--extended-probes`, `--merged`

**Output files** (in `work_dir`, named by VCF stem):

| File | Description |
|------|-------------|
| `{stem}_{N}.vcf` | Sampled VCF |
| `{stem}_{N}_norm.vcf` | Normalized VCF (multi-allelic decomposed) |
| `{stem}_{N}_vep115_golden.vcf` | Golden VEP 115 output |
| `{stem}_{N}_comparison_report.txt` | Per-field accuracy report |
| `{stem}_{N}_discrepancies.txt` | Per-variant discrepancies |

### E2E annotation benchmark (`bench_annotate_vcf`)

Full pipeline: read VCF, annotate, write annotated VCF. Measures end-to-end throughput.

```bash
cargo run --release --example bench_annotate_vcf -- \
  --input vep-benchmark/data/golden/input_1000.vcf \
  --cache vep-benchmark/data/golden/cache \
  --output /tmp/annotated.vcf \
  --everything --extended-probes \
  --reference-fasta vep-benchmark/data/golden/reference_chr1.fa

# With fjall backend and gzip output
cargo run --release --features kv-cache --example bench_annotate_vcf -- \
  --input vep-benchmark/data/HG002_chr1.vcf.gz \
  --cache /data/vep/wgs/fjall/115_GRCh38_vep \
  --output /tmp/annotated_chr1.vcf.gz \
  --backend fjall --everything --extended-probes \
  --reference-fasta /data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --compression gzip
```

**Options:** `--input`, `--cache`, `--output`, `--backend parquet|fjall`, `--everything`, `--extended-probes`, `--reference-fasta`, `--compression none|gzip|bgzf`, `--limit`

### Other benchmarks

| Example | Feature | Description |
|---------|---------|-------------|
| `bench_annotate` | `kv-cache` | Arrow-level annotation performance (fjall) |
| `lookup_parquet_bench` | — | `lookup_variants()` with parquet backend |
| `lookup_kv_bench` | `kv-cache` | `lookup_variants()` with fjall backend |
| `bench_fjall_gets` | `kv-cache` | Raw fjall KV store read throughput |
| `bench_sift_queries` | `kv-cache` | SIFT prediction query performance |

Run any benchmark with:

```bash
cargo run --release [--features kv-cache] --example <name> -- <args>
```

## Cargo features

| Feature | Default | Description |
|---------|---------|-------------|
| `kv-cache` | yes | Fjall LSM-tree backend with zstd compression |

## License

Apache 2.0
