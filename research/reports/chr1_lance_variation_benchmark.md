# Lance Variation chr1 Benchmark

## Summary

- Lance Python package: `7.0.0`
- PyArrow: `22.0.0`
- Cache dir: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`
- Chromosome: `chr1`
- Variants: `2.1-unpacked, 2.2-packed`
- Cold sample size: `2000` position keys
- VEPyr workload VCF: `/Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr1/input_chr1.vcf.gz`
- VEPyr workload batch size: `2000` probe keys
- Cold fragment row range: `512, 1024, 2048, 4096, 8192, 16384`
- Row limit: `full chr1`

## Results

| variant | tier | operation | rows | seconds | rows/s | artifact GiB |
|---|---|---:|---:|---:|---:|---:|
| parquet-current | warm | full_scan | 3628123 | 0.201 | 18037205 | 0.330 |
| parquet-current | cold | point_lookup_2000 | 2212 | 5.294 | 418 | 3.303 |
| 2.1-unpacked | warm | full_scan | 3628123 | 0.266 | 13623110 | 0.836 |
| 2.1-unpacked | warm | vepyr_warm_probe_batches_2000 | 668681 | 6.391 | 104633 | 0.836 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.575 | 619 | 4.650 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 66496 | 13.447 | 4945 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.569 | 620 | 4.650 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 66496 | 13.526 | 4916 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.323 | 666 | 4.650 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 66496 | 13.512 | 4921 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.413 | 648 | 4.650 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 66496 | 13.599 | 4890 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.599 | 615 | 4.650 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 66496 | 14.189 | 4687 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.647 | 607 | 4.650 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 66496 | 14.599 | 4555 | 4.650 |
| 2.2-packed | warm | full_scan | 3628123 | 0.251 | 14478249 | 0.927 |
| 2.2-packed | warm | vepyr_warm_probe_batches_2000 | 668681 | 7.158 | 93421 | 0.927 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.691 | 472 | 13.493 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 66496 | 22.126 | 3005 | 13.493 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.825 | 458 | 13.493 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 66496 | 22.850 | 2910 | 13.493 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 5.008 | 442 | 13.494 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 66496 | 22.820 | 2914 | 13.494 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 5.054 | 438 | 13.494 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 66496 | 23.208 | 2865 | 13.494 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.540 | 487 | 13.493 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 66496 | 23.838 | 2789 | 13.493 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.424 | 500 | 13.494 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 66496 | 25.120 | 2647 | 13.494 |

## VEPyr Workload Detail

| variant | tier | operation | row-group rows | probe keys | unique probes | selected positions | selected rows | batches | seconds |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| 2.1-unpacked | warm | vepyr_warm_probe_batches_2000 |  | 726462 | 709301 | 485407 | 668681 | 364 | 6.391 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 512 | 229285 | 223894 | 46785 | 66496 | 115 | 13.447 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 1024 | 229285 | 223894 | 46785 | 66496 | 115 | 13.526 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 2048 | 229285 | 223894 | 46785 | 66496 | 115 | 13.512 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 4096 | 229285 | 223894 | 46785 | 66496 | 115 | 13.599 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 8192 | 229285 | 223894 | 46785 | 66496 | 115 | 14.189 |
| 2.1-unpacked | cold | vepyr_cold_probe_batches_2000 | 16384 | 229285 | 223894 | 46785 | 66496 | 115 | 14.599 |
| 2.2-packed | warm | vepyr_warm_probe_batches_2000 |  | 726462 | 709301 | 485407 | 668681 | 364 | 7.158 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 512 | 229285 | 223894 | 46785 | 66496 | 115 | 22.126 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 1024 | 229285 | 223894 | 46785 | 66496 | 115 | 22.850 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 2048 | 229285 | 223894 | 46785 | 66496 | 115 | 22.820 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 4096 | 229285 | 223894 | 46785 | 66496 | 115 | 23.208 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 8192 | 229285 | 223894 | 46785 | 66496 | 115 | 23.838 |
| 2.2-packed | cold | vepyr_cold_probe_batches_2000 | 16384 | 229285 | 223894 | 46785 | 66496 | 115 | 25.120 |

The VEPyr workload uses `input_chr1.vcf.gz` to build extended probe positions, probes warm first, then sends only warm-miss positions to cold in 2,000-key batches. It is a position-level I/O proxy and does not apply the Rust allele matcher or cold variant Bloom sidecar.
Parquet VEPyr workload measurement: `skipped`.

## Cold Fragment Detail

| variant | row-group rows | lookup rows | seconds | rows/s | artifact GiB |
|---|---:|---:|---:|---:|---:|
| 2.1-unpacked | 512 | 2212 | 3.575 | 619 | 4.650 |
| 2.1-unpacked | 512 | 66496 | 13.447 | 4945 | 4.650 |
| 2.1-unpacked | 1024 | 2212 | 3.569 | 620 | 4.650 |
| 2.1-unpacked | 1024 | 66496 | 13.526 | 4916 | 4.650 |
| 2.1-unpacked | 2048 | 2212 | 3.323 | 666 | 4.650 |
| 2.1-unpacked | 2048 | 66496 | 13.512 | 4921 | 4.650 |
| 2.1-unpacked | 4096 | 2212 | 3.413 | 648 | 4.650 |
| 2.1-unpacked | 4096 | 66496 | 13.599 | 4890 | 4.650 |
| 2.1-unpacked | 8192 | 2212 | 3.599 | 615 | 4.650 |
| 2.1-unpacked | 8192 | 66496 | 14.189 | 4687 | 4.650 |
| 2.1-unpacked | 16384 | 2212 | 3.647 | 607 | 4.650 |
| 2.1-unpacked | 16384 | 66496 | 14.599 | 4555 | 4.650 |
| 2.2-packed | 512 | 2212 | 4.691 | 472 | 13.493 |
| 2.2-packed | 512 | 66496 | 22.126 | 3005 | 13.493 |
| 2.2-packed | 1024 | 2212 | 4.825 | 458 | 13.493 |
| 2.2-packed | 1024 | 66496 | 22.850 | 2910 | 13.493 |
| 2.2-packed | 2048 | 2212 | 5.008 | 442 | 13.494 |
| 2.2-packed | 2048 | 66496 | 22.820 | 2914 | 13.494 |
| 2.2-packed | 4096 | 2212 | 5.054 | 438 | 13.494 |
| 2.2-packed | 4096 | 66496 | 23.208 | 2865 | 13.494 |
| 2.2-packed | 8192 | 2212 | 4.540 | 487 | 13.493 |
| 2.2-packed | 8192 | 66496 | 23.838 | 2789 | 13.493 |
| 2.2-packed | 16384 | 2212 | 4.424 | 500 | 13.494 |
| 2.2-packed | 16384 | 66496 | 25.120 | 2647 | 13.494 |

Parquet cold baseline: `5.294s` for `2212` rows from `2000` sampled keys.

## Build Artifacts

| variant | tier | fragment rows | rows | build s | index s | artifact GiB | built |
|---|---|---:|---:|---:|---:|---:|---|
| 2.1-unpacked | warm | 500000 | 3628123 | 0.000 | 0.000 | 0.836 | False |
| 2.1-unpacked | cold | 512 | 84525843 | 0.000 | 0.000 | 4.650 | False |
| 2.1-unpacked | cold | 1024 | 84525843 | 0.000 | 0.000 | 4.650 | False |
| 2.1-unpacked | cold | 2048 | 84525843 | 0.000 | 0.000 | 4.650 | False |
| 2.1-unpacked | cold | 4096 | 84525843 | 0.000 | 0.000 | 4.650 | False |
| 2.1-unpacked | cold | 8192 | 84525843 | 0.000 | 0.000 | 4.650 | False |
| 2.1-unpacked | cold | 16384 | 84525843 | 0.000 | 0.000 | 4.650 | False |
| 2.2-packed | warm | 500000 | 3628123 | 0.000 | 0.000 | 0.927 | False |
| 2.2-packed | cold | 512 | 84525843 | 0.000 | 0.000 | 13.493 | False |
| 2.2-packed | cold | 1024 | 84525843 | 0.000 | 0.000 | 13.493 | False |
| 2.2-packed | cold | 2048 | 84525843 | 0.000 | 0.000 | 13.494 | False |
| 2.2-packed | cold | 4096 | 84525843 | 0.000 | 0.000 | 13.494 | False |
| 2.2-packed | cold | 8192 | 84525843 | 0.000 | 0.000 | 13.493 | False |
| 2.2-packed | cold | 16384 | 84525843 | 0.000 | 0.000 | 13.494 | False |

## Layout Notes

- `2.1-unpacked` keeps the hot projection as top-level Lance columns using storage version 2.1.
- `2.2-packed` keeps `position_key`, `variant_keys`, `chrom`, and `start` top-level and packs match, identity, clinical, flag, and AF groups into structs.
- With pylance/lance 7.0.0, the writer-safe default strips only the `lance-encoding:packed=true` metadata key because it panics on packed structs with all-null children; pass `--unsafe-packed-metadata` to reproduce the intended metadata path.
- `minor_allele` and `minor_allele_freq` remain logical fields and are marked as constant-null candidates.
- Cold Lance lookups use a `BTREE` scalar index on `position_key` with `use_scalar_index=True`.

## Sources

- Lance encoding strategy: <https://lance.org/format/file/encoding/>
- Lance versioning: <https://lance.org/format/file/versioning/>
- Lance table format: <https://lance.org/format/table/>
- LanceDB scalar indexes: <https://docs.lancedb.com/indexing/scalar-index>
