# Lance Variation chr1 Benchmark

## Summary

- Lance Python package: `7.0.0`
- PyArrow: `22.0.0`
- Cache dir: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`
- Chromosome: `chr1`
- Variants: `2.1-unpacked, 2.2-packed`
- Cold sample size: `2000` position keys
- Cold fragment row range: `512, 1024, 2048, 4096, 8192, 16384`
- Row limit: `full chr1`

## Results

| variant | tier | operation | rows | seconds | rows/s | artifact GiB |
|---|---|---:|---:|---:|---:|---:|
| parquet-current | warm | full_scan | 3628123 | 0.196 | 18490431 | 0.330 |
| parquet-current | cold | point_lookup_2000 | 2212 | 5.307 | 417 | 3.303 |
| 2.1-unpacked | warm | full_scan | 3628123 | 0.239 | 15165968 | 0.836 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 1.563 | 1416 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.505 | 631 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.329 | 664 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.205 | 690 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.428 | 645 | 4.650 |
| 2.1-unpacked | cold | point_lookup_2000 | 2212 | 3.271 | 676 | 4.650 |
| 2.2-packed | warm | full_scan | 3628123 | 0.215 | 16864310 | 0.927 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.267 | 518 | 13.493 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.134 | 535 | 13.493 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.226 | 523 | 13.494 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.297 | 515 | 13.494 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 4.008 | 552 | 13.493 |
| 2.2-packed | cold | point_lookup_2000 | 2212 | 3.757 | 589 | 13.494 |

## Cold Fragment Detail

| variant | row-group rows | lookup rows | seconds | rows/s | artifact GiB |
|---|---:|---:|---:|---:|---:|
| 2.1-unpacked | 512 | 2212 | 1.563 | 1416 | 4.650 |
| 2.1-unpacked | 1024 | 2212 | 3.505 | 631 | 4.650 |
| 2.1-unpacked | 2048 | 2212 | 3.329 | 664 | 4.650 |
| 2.1-unpacked | 4096 | 2212 | 3.205 | 690 | 4.650 |
| 2.1-unpacked | 8192 | 2212 | 3.428 | 645 | 4.650 |
| 2.1-unpacked | 16384 | 2212 | 3.271 | 676 | 4.650 |
| 2.2-packed | 512 | 2212 | 4.267 | 518 | 13.493 |
| 2.2-packed | 1024 | 2212 | 4.134 | 535 | 13.493 |
| 2.2-packed | 2048 | 2212 | 4.226 | 523 | 13.494 |
| 2.2-packed | 4096 | 2212 | 4.297 | 515 | 13.494 |
| 2.2-packed | 8192 | 2212 | 4.008 | 552 | 13.493 |
| 2.2-packed | 16384 | 2212 | 3.757 | 589 | 13.494 |

Parquet cold baseline: `5.307s` for `2212` rows from `2000` sampled keys.

## Build Artifacts

| variant | tier | fragment rows | rows | build s | index s | artifact GiB | built |
|---|---|---:|---:|---:|---:|---:|---|
| 2.1-unpacked | warm | 500000 | 3628123 | 11.640 | 0.231 | 0.836 | True |
| 2.1-unpacked | cold | 512 | 84525843 | 28.850 | 5.274 | 4.650 | True |
| 2.1-unpacked | cold | 1024 | 84525843 | 29.316 | 5.416 | 4.650 | True |
| 2.1-unpacked | cold | 2048 | 84525843 | 29.209 | 5.183 | 4.650 | True |
| 2.1-unpacked | cold | 4096 | 84525843 | 29.246 | 5.202 | 4.650 | True |
| 2.1-unpacked | cold | 8192 | 84525843 | 29.503 | 5.536 | 4.650 | True |
| 2.1-unpacked | cold | 16384 | 84525843 | 29.003 | 5.259 | 4.650 | True |
| 2.2-packed | warm | 500000 | 3628123 | 2.231 | 0.209 | 0.927 | True |
| 2.2-packed | cold | 512 | 84525843 | 32.091 | 5.157 | 13.493 | True |
| 2.2-packed | cold | 1024 | 84525843 | 32.211 | 5.768 | 13.493 | True |
| 2.2-packed | cold | 2048 | 84525843 | 31.611 | 5.314 | 13.494 | True |
| 2.2-packed | cold | 4096 | 84525843 | 31.342 | 7.097 | 13.494 | True |
| 2.2-packed | cold | 8192 | 84525843 | 37.562 | 5.659 | 13.493 | True |
| 2.2-packed | cold | 16384 | 84525843 | 40.794 | 6.271 | 13.494 | True |

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
