# datafusion-bio-function-vep

VEP-oriented DataFusion functions and benchmark tooling.

## Implementation Status (as of March 4, 2026)

| Area | Status | Notes |
|---|---|---|
| `lookup_variants(...)` table function | Implemented | Parquet + Fjall backends, match modes (`exact`, `exact_or_colocated_ids`, `exact_or_vep_existing`) |
| `annotate_vep(...)` table function | Phase 2 (transcript/exon baseline) | Uses `lookup_variants` metadata + transcript/exon context (when available) for ranked SO output; falls back to phase-1.5 placeholder CSQ when context tables are absent |
| Ensembl VEP 115 golden benchmark pipeline | Implemented | Samples 1000 variants, runs VEP 115 in Docker, compares with `annotate_vep` output |
| Full consequence engine parity (41 SO terms) | In progress | Design and cache contract documented in `PORTING_PLAN.md` |

## Golden Benchmark: `annotate_vep` vs Ensembl VEP 115

### What it does

The `annotate_vep_golden_bench` example:
1. Ensures local benchmark input VCF is available.
2. Samples the first 1000 variants from `HG002_chr22.vcf.gz`.
3. Runs Ensembl VEP release `115.2` in Docker as golden output.
4. Runs `annotate_vep(...)` on the same sampled input.
5. Writes a comparison report with row/key overlap and CSQ agreement metrics.

### Prerequisites

- Docker installed and running.
- VEP cache directory available (default used below):
  - `/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/cache`
- Annotation cache source for `annotate_vep` (example parquet path):
  - `/Users/mwiewior/research/data/vep/115_GRCh38_variants.parquet`

### Benchmark input in this repo

- `vep-benchmark/data/HG002_chr22.vcf.gz`
- `vep-benchmark/data/HG002_chr22.vcf.gz.tbi`

### Run command

```bash
cargo run -p datafusion-bio-function-vep --example annotate_vep_golden_bench --release -- \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variants.parquet \
  parquet \
  1000 \
  /Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/cache \
  vep-benchmark/data/HG002_chr22.vcf.gz \
  /tmp/annotate_vep_golden_bench
```

### Output

Report file:
- `/tmp/annotate_vep_golden_bench/HG002_chr22_1000_comparison_report.txt`

Current expected behavior:
- key overlap should be high (variant identity matches),
- `ours_with_csq` should be non-zero for cache-matched variants (and typically higher when transcript/exon context tables are provided),
- `csq_exact_matches` should still be low because full transcript/regulatory consequence computation is not implemented yet.

## Other benchmark examples

- `lookup_parquet_bench`
- `lookup_kv_bench`
- `bench_annotate`
- `bench_fjall_gets`

Run with:

```bash
cargo run -p datafusion-bio-function-vep --example <example_name> --release -- <args>
```
