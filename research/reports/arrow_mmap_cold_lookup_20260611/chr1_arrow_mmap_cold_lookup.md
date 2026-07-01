# Arrow Mmap Cold Lookup Benchmark

## Workload

- Cache root: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`
- VCF: `/Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr1/input_chr1.vcf.gz`
- Lance dataset: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/chr1.lance`
- Arrow store root: `/Users/mwiewior/workspace/data_vepyr/arrow_mmap_chr1_bench`
- Chrom: `chr1`
- Parsed variants: `323430`
- Raw probe attempts: `390457`
- Warm-miss unique positions: `25557`
- Position-index admitted unique positions: `16487`
- Bloom-admitted unique positions: `8749`
- Bloom-admitted attempts: `8780`
- Projected columns: `43`
- Row-location index kind: `packed_zstd`
- Row-location index path: `/Users/mwiewior/workspace/data_vepyr/arrow_mmap_chr1_bench/chr1.range_index.pzstd`
- Row-location index block entries: `4096`
- Omitted requested columns absent from cold parquet: `variant_keys`

## Store Sizes

Total size is data-store directory size plus the shared packed range index.

| compression | build seconds | rows | unique positions | batches | data size | index size | total size |
|---|---:|---:|---:|---:|---:|---:|---:|
| uncompressed | 467.641 | 84525843 | 80564179 | 5160 | 20862376872 | 377771759 | 21240148631 |
| lz4 | 478.483 | 84525843 | 80564179 | 5160 | 7351782576 | 377771759 | 7729554335 |
| zstd | 493.170 | 84525843 | 80564179 | 5160 | 4072286664 | 377771759 | 4450058423 |

## Lance Cold Lookup

| plan | query batch | scans | filter keys | rows | result batches | selected positions | seconds | bytes read | requests | ranges | fragments | parts loaded | index comparisons |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| small_chunks | 238 | 37 | 8749 | 10036 | 119 | 8749 | 5.745 | 1239707055 | 11341 | 8469 | 119 | 5667 | 23269413 |
| large_chunks | 2000 | 5 | 8749 | 10036 | 87 | 8749 | 3.491 | 442509704 | 5025 | 8468 | 87 | 0 | 23220229 |

## Arrow Mmap Lookup

| compression | mmap direct | total size | data size | index size | open seconds | lookup seconds | positions requested | positions found | selected rows | projected columns | estimated bytes touched | batch files opened | index gets | missing positions | checksum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uncompressed | true | 21240148631 | 20862376872 | 377771759 | 0.002 | 27.397 | 8749 | 8749 | 10036 | 43 | 1846555 | 3515 | 8749 | 0 | 736780182109964908 |
| lz4 | false | 7729554335 | 7351782576 | 377771759 | 0.001 | 21.287 | 8749 | 8749 | 10036 | 43 | 1846555 | 3515 | 8749 | 0 | 736780182109964908 |
| zstd | false | 4450058423 | 4072286664 | 377771759 | 0.001 | 18.715 | 8749 | 8749 | 10036 | 43 | 1846555 | 3515 | 8749 | 0 | 736780182109964908 |

## Verification

- Arrow row counts match: `true`
- Arrow found-position counts match: `true`
- Arrow checksums match: `true`
- Lance large-chunks rows vs Arrow rows: `Some(10036)` vs `Some(10036)`
- Lance large-chunks positions vs Arrow found positions: `Some(8749)` vs `Some(8749)`
