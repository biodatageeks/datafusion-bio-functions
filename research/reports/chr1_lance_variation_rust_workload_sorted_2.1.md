# Rust Lance chr1 Variation Workload Benchmark

- Lance crate: `7.0.0`
- Cache root: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`
- Lance root: `/Users/mwiewior/workspace/data_vepyr/lance_variation_chr1_bench`
- VCF: `/Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr1/input_chr1.vcf.gz`
- Position index: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.position_index/chr1.posidx` (0.009s load)
- Variant bloom index: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.variant_bloom_index/chr1.varbf` (0.117s load)
- Workload query batch size: `2000`
- Warm scan batch size: `65536`
- Cold scan batch size: `8192`
- Cold positions sorted before Lance batches: `true`

## Summary

- Position index reduces the cold lookup set to `46785` unique positions, `50813` probe attempts, and `24` 2k batches.
- Adding variant bloom admission reduces that further to `9497` unique positions, `9638` probe attempts, and `5` 2k batches.
- Best warm full scan: `2.1-unpacked` at `0.397s` for `3628123` selected rows.
- Best default cold `posidx` scan: `2.1-unpacked` with `1024` fragment rows at `8.046s`, selecting `66461` rows across `46785` positions.
- Best optional `posidx_bloom` scan: `2.1-unpacked` with `8192` fragment rows at `3.544s`, selecting `11463` rows across `9497` positions.
- Current result favors `2.1-unpacked` for cold point lookups; `2.2-packed` is about 2x slower in this all-column Rust workload.

## Workload

| metric | value |
|---|---:|
| chrom | chr1 |
| chrom code | 1 |
| parsed VCF variants | 323430 |
| raw extended probe attempts | 727363 |
| unique probe positions | 709301 |
| warm-miss unique positions | 223894 |
| warm-miss probe attempts | 229286 |
| position-index admitted unique positions | 46785 |
| position-index admitted probe attempts | 50813 |
| variant-bloom admitted unique cold positions | 9497 |
| variant-bloom admitted cold probe attempts | 9638 |
| 2k cold Lance filter batches | 24 |
| 2k cold Lance filter batches after bloom | 5 |

## Rust Lance Results

| phase | layout | fragment rows | query batch | probe positions | selected rows | selected positions | result batches | seconds | rows/s |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| warm_full_scan_all_columns | 2.1-unpacked |  | full_scan |  | 3628123 | 2655552 | 58 | 0.397 | 9135868 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 512 | 2000 | 46785 | 66461 | 46785 | 106 | 8.123 | 8181 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 512 | 2000 | 9497 | 11463 | 9497 | 87 | 3.774 | 3038 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 1024 | 2000 | 46785 | 66461 | 46785 | 106 | 8.046 | 8261 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 1024 | 2000 | 9497 | 11463 | 9497 | 87 | 3.681 | 3114 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 2048 | 2000 | 46785 | 66461 | 46785 | 106 | 8.163 | 8142 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 2048 | 2000 | 9497 | 11463 | 9497 | 87 | 3.817 | 3003 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 4096 | 2000 | 46785 | 66461 | 46785 | 106 | 8.169 | 8136 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 4096 | 2000 | 9497 | 11463 | 9497 | 87 | 3.714 | 3087 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 8192 | 2000 | 46785 | 66461 | 46785 | 106 | 8.199 | 8106 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 8192 | 2000 | 9497 | 11463 | 9497 | 87 | 3.544 | 3235 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 16384 | 2000 | 46785 | 66461 | 46785 | 106 | 8.214 | 8091 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 16384 | 2000 | 9497 | 11463 | 9497 | 87 | 5.328 | 2151 |
