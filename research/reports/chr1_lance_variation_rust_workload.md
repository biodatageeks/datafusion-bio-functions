# Rust Lance chr1 Variation Workload Benchmark

- Lance crate: `7.0.0`
- Cache root: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`
- Lance root: `/Users/mwiewior/workspace/data_vepyr/lance_variation_chr1_bench`
- VCF: `/Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr1/input_chr1.vcf.gz`
- Position index: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.position_index/chr1.posidx` (0.007s load)
- Variant bloom index: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.variant_bloom_index/chr1.varbf` (0.098s load)
- Workload query batch size: `2000`
- Warm scan batch size: `65536`
- Cold scan batch size: `8192`

## Summary

- Position index reduces the cold lookup set to `46785` unique positions, `50813` probe attempts, and `24` 2k batches.
- Adding variant bloom admission reduces that further to `9497` unique positions, `9638` probe attempts, and `5` 2k batches.
- Best warm full scan: `2.2-packed` at `0.405s` for `3628123` selected rows.
- Best default cold `posidx` scan: `2.1-unpacked` with `512` fragment rows at `7.934s`, selecting `66461` rows across `46785` positions.
- Best optional `posidx_bloom` scan: `2.1-unpacked` with `1024` fragment rows at `3.575s`, selecting `11463` rows across `9497` positions.
- Current result favors `2.1-unpacked` for cold point lookups; `2.2-packed` is about 2x slower in this all-column Rust workload.
- Implemented merged 2.1-unpacked Lance materializer output for chr1 is `4.6G` total with `93` data files after projecting the 47 runtime `everything` columns and applying miniblock/zstd field encoding metadata.

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

## Implementation Validation

| metric | value |
|---|---:|
| materialized dataset | `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/chr1.lance` |
| layout | `2.1-unpacked` |
| projected runtime columns | 47 |
| warm rows | 3628123 |
| cold rows | 84525843 |
| warm rows/file | 500000 |
| warm row-group rows | 262144 |
| cold rows/file | 1000000 |
| cold row-group rows | 1024 |
| build elapsed | 64.72s |
| artifact size | 4.6G |
| data file count | 93 |
| total file count | 98 |
| e2e command | `uv run python run_annotation_fast.py chr1 --cache merged --backend lance --forks 0 --force` |
| e2e annotated variants | 323430 |
| e2e annotation elapsed | 51.3s |
| e2e annotation throughput | 6309 variants/s |
| e2e peak memory | 1937 MB |
| e2e only in vepyr | 0 |
| e2e only in VEP | 0 |
| e2e CSQ count mismatches | 0 |
| e2e CSQ order mismatches | 0 |
| e2e shared CSQ fields at 100% | 86 |

## Rust Lance Results

| phase | layout | fragment rows | query batch | probe positions | selected rows | selected positions | result batches | seconds | rows/s |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| warm_full_scan_all_columns | 2.1-unpacked |  | full_scan |  | 3628123 | 2655552 | 58 | 0.405 | 8966202 |
| warm_full_scan_all_columns | 2.2-packed |  | full_scan |  | 3628123 | 2655552 | 58 | 0.405 | 8968230 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 512 | 2000 | 46785 | 66461 | 46785 | 106 | 7.934 | 8377 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 512 | 2000 | 9497 | 11463 | 9497 | 87 | 3.617 | 3169 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 1024 | 2000 | 46785 | 66461 | 46785 | 106 | 7.976 | 8333 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 1024 | 2000 | 9497 | 11463 | 9497 | 87 | 3.575 | 3207 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 2048 | 2000 | 46785 | 66461 | 46785 | 106 | 7.989 | 8319 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 2048 | 2000 | 9497 | 11463 | 9497 | 87 | 3.604 | 3180 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 4096 | 2000 | 46785 | 66461 | 46785 | 106 | 8.096 | 8209 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 4096 | 2000 | 9497 | 11463 | 9497 | 87 | 3.732 | 3071 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 8192 | 2000 | 46785 | 66461 | 46785 | 106 | 8.295 | 8012 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 8192 | 2000 | 9497 | 11463 | 9497 | 87 | 5.004 | 2291 |
| cold_posidx_probe_batches_all_columns | 2.1-unpacked | 16384 | 2000 | 46785 | 66461 | 46785 | 106 | 8.117 | 8188 |
| cold_posidx_bloom_probe_batches_all_columns | 2.1-unpacked | 16384 | 2000 | 9497 | 11463 | 9497 | 87 | 5.043 | 2273 |
| cold_posidx_probe_batches_all_columns | 2.2-packed | 512 | 2000 | 46785 | 66461 | 46785 | 106 | 16.462 | 4037 |
| cold_posidx_bloom_probe_batches_all_columns | 2.2-packed | 512 | 2000 | 9497 | 11463 | 9497 | 87 | 10.419 | 1100 |
| cold_posidx_probe_batches_all_columns | 2.2-packed | 1024 | 2000 | 46785 | 66461 | 46785 | 106 | 16.764 | 3964 |
| cold_posidx_bloom_probe_batches_all_columns | 2.2-packed | 1024 | 2000 | 9497 | 11463 | 9497 | 87 | 10.808 | 1061 |
| cold_posidx_probe_batches_all_columns | 2.2-packed | 2048 | 2000 | 46785 | 66461 | 46785 | 106 | 16.944 | 3922 |
| cold_posidx_bloom_probe_batches_all_columns | 2.2-packed | 2048 | 2000 | 9497 | 11463 | 9497 | 87 | 10.790 | 1062 |
| cold_posidx_probe_batches_all_columns | 2.2-packed | 4096 | 2000 | 46785 | 66461 | 46785 | 106 | 16.926 | 3926 |
| cold_posidx_bloom_probe_batches_all_columns | 2.2-packed | 4096 | 2000 | 9497 | 11463 | 9497 | 87 | 10.849 | 1057 |
| cold_posidx_probe_batches_all_columns | 2.2-packed | 8192 | 2000 | 46785 | 66461 | 46785 | 106 | 17.649 | 3766 |
| cold_posidx_bloom_probe_batches_all_columns | 2.2-packed | 8192 | 2000 | 9497 | 11463 | 9497 | 87 | 10.975 | 1044 |
| cold_posidx_probe_batches_all_columns | 2.2-packed | 16384 | 2000 | 46785 | 66461 | 46785 | 106 | 19.602 | 3390 |
| cold_posidx_bloom_probe_batches_all_columns | 2.2-packed | 16384 | 2000 | 9497 | 11463 | 9497 | 87 | 11.922 | 961 |
