# Warm Cache Cold Parquet Variation Lookup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Benchmark and integrate the indexed-parquet warm/cold variation path where warm misses that are present in the cold position index are resolved from `chrN_cold.parquet` in buffer-sized batches instead of by per-position `variation.fjall` gets, with cold `variant_keys` removed from the cache format.

**Architecture:** Keep the current warm tier and `.posidx` routing. The cold variation parquet schema stores scalar `position_key`, allele fields, and requested annotation columns, but does not store `variant_keys: List<UInt64>`. Runtime cold lookup collects per-buffer cold-positive warm misses, groups them by sorted `position_key` and Parquet row group/page, reads projected cold Parquet ranges, performs existing allele verification, and emits the same cache columns as the Fjall baseline.

**Tech Stack:** Rust, DataFusion 50/53 APIs in this workspace, Arrow 56, parquet-rs 56, existing `warm_cache` modules, existing `kv_cache/cache_exec.rs`, indexed-parquet cache generation.

---

## Source Context

DataFusion's Parquet pruning pipeline is directly relevant to this benchmark. The DataFusion blog "Parquet Pruning in DataFusion: Read Only What Matters" describes these reader-side mechanisms:

- metadata read and metadata caching
- projection pruning
- row-group pruning from min/max statistics
- row-group pruning from Bloom filters for equality predicates
- page-level pruning when page indexes exist
- final byte-range access planning
- row-level filter pushdown as an experimental/coming feature rather than a default dependency

Reference: <https://datafusion.apache.org/blog/2025/03/20/parquet-pruning/>

Echtvar is also directly relevant because its fast path is a genomic-bin cursor over sorted query variants. Its README states that it chunks the genome into `1<<20` bases and encodes compact variants into 32-bit integers. The NAR paper describes the runtime behavior: when the query variant enters a different `2^20` base bin, echtvar reads that bin's `var32`, supplemental long-variant table, and fields into memory; if the next query stays in the same bin, the in-memory data is reused; exact lookup is a binary search over the in-memory encoded variant table.

References:

- <https://github.com/brentp/echtvar>
- <https://academic.oup.com/nar/article/51/1/e3/6775383>

## Current State Observed

Input and cache:

```text
input VCF:
  /Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr4/input_chr4.vcf.gz
  307,295 data rows

cache:
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
```

Existing chr4 variation files:

```text
variation/chr4_warm.parquet        284 MiB, 3,128,467 rows, 7 row groups
variation/chr4_cold.parquet        2.5 GiB, 70,335,895 rows, 141 row groups
variation.position_index/chr4.posidx 23 MiB
variation.fjall                    110 GiB for all chromosomes
translation_sift.fjall              51 GiB for all chromosomes
```

The current warm decision is in `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`:

```text
warm exact hit                         -> emit warm row, skip variation.fjall
warm position-covered no exact match   -> definitive miss, skip variation.fjall
warm not-covered + cold index absent   -> definitive miss, skip variation.fjall
warm not-covered + cold index present  -> use variation.fjall
```

The proposed benchmark replaces only the last branch with cold Parquet lookup. Translation/SIFT still uses the current Fjall path when `--everything` and `use_fjall` are active.

One local baseline attempt with an existing release binary returned zero rows because the binary appears stale relative to the current code that canonicalizes `_warm/_cold` variation filenames. Rebuild before timing:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release \
  -p datafusion-bio-function-vep --features kv-cache \
  --example bench_annotate_vcf
```

## Recommended Approach

Use a two-phase implementation.

Phase 1 is a standalone benchmark binary. It can be built and validated without perturbing the annotation pipeline:

```text
examples/bench_warm_cold_parquet_chr4.rs
```

It parses chr4 VCF probes, opens `chr4_warm.parquet`, `chr4_cold.parquet`, `.posidx`, and `variation.fjall`, then runs two lookup modes over the same probes:

```text
mode=warm-fjall:
  current behavior; cold-positive probes call variation.fjall

mode=warm-cold-parquet:
  warm exact/miss/negative behavior unchanged; cold-positive probes collected per buffer and resolved from chr4_cold.parquet
```

Phase 2 wires the successful backend into `cache_exec.rs` behind the indexed-parquet runtime path:

```text
VEP_WARM_COLD_VARIATION_BACKEND=parquet
```

The indexed-parquet path is the target behavior. Do not add compatibility logic that accepts old cold parquet files containing `variant_keys`; rebuild those caches instead.

## Files

- Create: `datafusion/bio-function-vep/src/warm_cache/cold_parquet.rs`
  - Cold Parquet reader, row-group metadata, per-buffer probe grouping, projected row-group loading, allele verification.
- Modify: `datafusion/bio-function-vep/src/warm_cache/mod.rs`
  - Export `cold_parquet`.
- Create: `datafusion/bio-function-vep/examples/bench_warm_cold_parquet_chr4.rs`
  - Microbenchmark for chr4: current warm+Fjall vs warm+cold-Parquet.
- Modify: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`
  - Add opt-in backend switch and buffer-level cold probe collection after benchmark proves correctness/performance.
- Modify: `datafusion/bio-function-vep/src/warm_cache/build.rs`
  - Add layout knobs for row-group size, page index/statistics verification, and optional Bloom filters on scalar lookup columns.
  - Project `variant_keys` out of the cold variation parquet output.
  - Keep warm `variant_keys` because warm chunk lookup still builds an in-memory exact variant index from it.
- Optional create: `datafusion/bio-function-vep/examples/rewrite_cold_variation_layout.rs`
  - Rewrite selected chromosomes with alternative row-group sizes and writer properties for layout experiments.

## Implementation Tasks

### Task 1: Standalone Benchmark Harness

- [ ] Create `bench_warm_cold_parquet_chr4.rs`.
- [ ] Reuse probe parsing from `bench_warm_tier_chr1.rs`, but make chromosome configurable and default to `chr4`.
- [ ] Add arguments:

```text
--vcf /Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr4/input_chr4.vcf.gz
--cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
--chrom chr4
--buffer-size 5000
--max-variants optional
--mode warm-fjall|warm-cold-parquet|both
```

- [ ] Report:

```text
variants
probes
warm_hits
warm_definitive_misses
warm_not_covered
cold_index_negative_skips
cold_index_positive_lookups
fjall_gets
fjall_bytes
cold_parquet_row_groups_read
cold_parquet_rows_decoded
cold_parquet_projected_bytes_estimate
matched_alleles
elapsed_s
probes_per_s
```

### Task 2: Cold Parquet Reader

- [ ] Implement `ColdParquetLookup::open(path, projection, batch_size)`.
- [x] Treat `variant_keys` in a cold variation parquet schema as invalid for the indexed-parquet format.

Expected cold-schema validation:

```text
if schema.contains("variant_keys"):
  error: cold variation parquet contains deprecated variant_keys; rebuild indexed_parquet cache
```

- [ ] Load `ArrowReaderMetadata` once and keep it for the chromosome lifetime.
- [ ] Extract row-group ranges from `position_key` min/max statistics.
- [ ] Extract an echtvar-style genomic bin for each row group:

```text
echtvar_bin = (position_key & 0x0000_FFFF_FFFF_FFFF) >> 20
```

- [ ] Build a compact in-memory `Vec<Range<usize>>` from `echtvar_bin` to candidate row-group indices for the chromosome. This is the Parquet analogue of echtvar's `chrom/bin` directory lookup.
- [ ] Keep a monotonic cursor over row groups because VCF probes are position sorted.
- [ ] Project only:

```text
position_key
allele_string
end
failed
variation_name
requested cache columns
colocated columns only when needed
```

- [x] Do not read or tolerate `variant_keys` in cold Parquet. For cold positives, `position_key` narrows candidates; exact correctness comes from `allele_string`, `end`, and `failed` verification.
- [ ] Add a bounded current/next row-group cache; a 5000-variant buffer should never reload the same cold row group more than once.

### Task 3: Per-Buffer Probe Collection

- [ ] In the benchmark, after warm probing:

```text
if warm exact:
  count hit
elif warm position covered:
  count definitive miss
elif cold posidx absent:
  count cold negative skip
else:
  push ColdProbe
```

- [ ] Sort and deduplicate `ColdProbe` by:

```text
position_key, normalized ref, normalized alt, vcf_row, probe_start
```

- [ ] Group probes by row group using row-group min/max and the monotonic cursor.
- [ ] For each row group, load projected rows once and binary-search `position_key` to candidate row ranges.
- [ ] Verify with existing `allele_matches` and failed/end interval checks from `cache_exec.rs`.
- [ ] Assert `warm-fjall` and `warm-cold-parquet` produce equal `matched_alleles` and equal emitted cache rows in the microbenchmark.

### Task 4: Runtime Integration

- [ ] Make indexed parquet the normal variation path for `cache_format=indexed_parquet`.
- [ ] Keep any Fjall comparison path under `cache_format=legacy_fjall` only. Do not silently fall back from indexed parquet to `variation.fjall`.
- [ ] Keep only runtime tuning config for the parquet path:

```text
VEP_COLD_PARQUET_BATCH_SIZE=262144
VEP_COLD_PARQUET_ROW_GROUP_CACHE=2
```

- [ ] In `KvLookupStream::poll_next`, collect cold-positive probes for the current input batch instead of immediately falling through to Fjall under indexed parquet.
- [ ] Resolve the collected probes before finalizing cache output builders for the batch.
- [ ] Do not open `variation.fjall` under indexed parquet.
- [ ] Use compact parquet `translation_sift/` under indexed parquet; `translation_sift.fjall` belongs only to `cache_format=legacy_fjall`.

### Task 4b: Remove Cold `variant_keys` From Cache Generation

- [x] Keep key augmentation internally while partitioning/building the warm/cold tier:

```text
input partition rows
  -> append position_key
  -> append variant_keys temporarily
  -> split warm/cold by position
```

- [x] Write warm parquet with `variant_keys` because warm lookup uses it for the in-memory exact variant index.
- [x] Write cold parquet without `variant_keys`.
- [x] Build `variation.variant_bloom_index/chrN.varbf` from cold `position_key` + `allele_string` after dropping `variant_keys` from the cold output schema.
- [x] Do not add a reader fallback for old cold files that still contain `variant_keys`.
- [x] Do not add a migration reader. The migration path is to rebuild the cache as `cache_format=indexed_parquet`.
- [x] Add a schema test that fails if `chrN_cold.parquet` includes `variant_keys`.
- [x] Add a schema test that confirms `chrN_warm.parquet` still includes `variant_keys`.
- [ ] Add an end-to-end hash test showing indexed parquet without cold `variant_keys` matches the Fjall baseline output.

### Task 5: Layout Experiments

Current `chr4_cold.parquet` has 141 row groups of roughly 500k rows. For sparse warm misses, this may decode too much per buffer. Benchmark these cold layouts:

```text
4k rows/row-group
8k rows/row-group
16k rows/row-group
64k rows/row-group
128k rows/row-group
250k rows/row-group
500k rows/row-group (current)
```

All layouts must remain position-aligned:

```text
no position_key may span two row groups
sort order: position_key, start, end, allele_string
```

Also test an echtvar-aligned layout variant for each row-group size:

```text
do not let a row group cross a 2^20-base genomic bin
within each bin, split by row-group target rows: 4k, 8k, 16k, 64k, 128k, 250k, 500k
write metadata: vepyr.echtvar_bin_bits = 20
```

This separates two effects:

```text
row-count granularity controls decode size
2^20 genomic binning controls cursor locality and avoids row groups spanning natural query windows
```

For every rewritten cold layout, rebuild the persisted cold position index from that exact cold Parquet file:

```text
variation.rowgroup_4k/chr4_cold.parquet
variation.position_index.rowgroup_4k/chr4.posidx

variation.rowgroup_8k/chr4_cold.parquet
variation.position_index.rowgroup_8k/chr4.posidx

...
```

The bitset contents should be logically identical across row-group sizes when the cold row set is unchanged, but rebuilding from the candidate file is required for benchmark isolation. It proves the layout is self-contained, catches rewrite bugs, and keeps the benchmark from accidentally mixing a 500k-layout parquet with a stale or differently generated index. Record `positions`, `storage_bytes`, and a checksum of the serialized `.posidx` in the benchmark output.

Writer properties to test:

```rust
WriterProperties::builder()
    .set_compression(Compression::ZSTD(Default::default()))
    .set_statistics_enabled(EnabledStatistics::Page)
    .set_column_index_truncate_length(Some(64))
    .set_offset_index_disabled(false)
    .set_column_bloom_filter_enabled(ColumnPath::from("position_key"), true)
    .set_column_bloom_filter_fpp(ColumnPath::from("position_key"), 0.01)
```

Bloom filters should be tested on scalar `position_key`. Do not write nested `variant_keys` to cold parquet, and do not add an exploded scalar `variant_key` sidecar in this iteration. The existing `variation.variant_bloom_index` sidecar is built during cache generation from temporary augmented keys, then the cold parquet payload is written without those keys.

### Task 6: DataFusion Reader Feature Experiments

Benchmark two cold lookup engines:

1. Direct parquet-rs row-group reader:
   - minimal overhead
   - explicit monotonic cursor
   - no SQL planning
   - uses row-group stats manually

2. DataFusion scan with predicates:
   - query shape: `position_key IN (...)`
   - enables DataFusion row-group stat pruning and Bloom filter pruning when the file has Bloom filters
   - can use page index pruning when page stats exist
   - exposes metrics such as row groups pruned by statistics/Bloom/page index

Expected decision rule:

```text
direct reader wins when row groups are already narrowly selected by monotonic cursor
DataFusion wins if Bloom/page pruning skips enough row groups/pages to offset planning overhead
```

### Task 6b: Echtvar-Style Cursor Mapping

- [ ] Add `ColdBinIndex` for a chromosome:

```rust
struct ColdBinIndex {
    bin_bits: u8,                 // 20 for echtvar-compatible bins
    row_groups_by_bin: Vec<Range<usize>>,
}
```

- [ ] Populate it from row-group `position_key` min/max:

```text
min_bin = position_from_key(min_position_key) >> 20
max_bin = position_from_key(max_position_key) >> 20
```

- [ ] For the strict echtvar-aligned layout, assert `min_bin == max_bin` for every row group.
- [ ] For the current 500k layout, allow row groups to span bins and insert the row group into every overlapped bin range.
- [ ] During lookup, compute `probe_bin = position >> 20`, fetch candidate row groups from `row_groups_by_bin[probe_bin]`, then use the existing monotonic row-group cursor inside that candidate range.
- [ ] Keep a current-bin decoded cache:

```text
current_bin
decoded row groups for current_bin, bounded by VEP_COLD_PARQUET_ROW_GROUP_CACHE
```

- [ ] Measure:

```text
bin_changes
same_bin_probes
row_groups_per_bin_p50/p95/max
row_group_cache_hits
row_group_cache_misses
```

- [ ] Compare three layouts:

```text
row-count only: current position-aligned row groups, may cross 2^20 bins
echtvar-aligned row groups: never cross 2^20 bins
indexed-parquet default: 8k row groups with scalar position_key pruning and external variant bloom index
```

### Task 7: Correctness Verification

- [ ] Unit tests for row-group selection:

```text
probe before first range -> no row group
probe inside range -> expected row group
probe between ranges -> no row group
monotonic cursor never moves backward unless previous-row-group cache handles boundary probe
```

- [ ] Unit tests for cold lookup matching:

```text
exact SNV match emits same variation_name as Fjall fixture
position present but allele mismatch is definitive miss
failed > allowed_failed is skipped
end/start interval visibility matches cache_exec.rs
```

- [ ] Microbenchmark correctness:

```bash
cargo run --release -p datafusion-bio-function-vep --features kv-cache \
  --example bench_warm_cold_parquet_chr4 -- \
  --mode both \
  --vcf /Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr4/input_chr4.vcf.gz \
  --cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --chrom chr4 \
  --buffer-size 5000
```

Expected:

```text
matched_alleles and emitted cache rows equal between modes
warm/cold routing counters equal except fjall_gets vs cold_parquet_reads
```

- [ ] Full annotation hash/diff against the legacy baseline:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release \
  -p datafusion-bio-function-vep --features kv-cache --example bench_annotate_vcf

target/release/examples/bench_annotate_vcf \
  --input /Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr4/input_chr4.vcf.gz \
  --cache /Users/mwiewior/workspace/data_vepyr/legacy_fjall/115_GRCh38_merged \
  --output /tmp/chr4_legacy_fjall.vcf \
  --cache-format legacy_fjall \
  --everything \
  --extended-probes \
  --reference-fasta /Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --forks 1 \
  --buffer-size 5000 \
  --compression none \
  --no-progress

target/release/examples/bench_annotate_vcf \
  --input /Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr4/input_chr4.vcf.gz \
  --cache /Users/mwiewior/workspace/data_vepyr/parquet/115_GRCh38_merged \
  --output /tmp/chr4_indexed_parquet.vcf \
  --cache-format indexed_parquet \
  --everything \
  --extended-probes \
  --reference-fasta /Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --forks 1 \
  --buffer-size 5000 \
  --compression none \
  --no-progress

shasum -a 256 /tmp/chr4_legacy_fjall.vcf /tmp/chr4_indexed_parquet.vcf
diff -u /tmp/chr4_legacy_fjall.vcf /tmp/chr4_indexed_parquet.vcf
```

Expected:

```text
matching SHA-256 hashes and no diff
```

## Benchmark Matrix

Run all cases after rebuilding.

| Case | Variation path | Translation/SIFT | Purpose |
|---|---|---|---|
| A | `cache_format=legacy_fjall`, current warm + `variation.fjall` | `translation_sift.fjall` | Correctness/performance baseline only |
| B | `cache_format=indexed_parquet`, warm + cold Parquet direct reader, no cold `variant_keys` | compact `translation_sift/` | Target default |
| C | `cache_format=indexed_parquet`, warm + cold Parquet through DataFusion predicates, no cold `variant_keys` | compact `translation_sift/` | Test Bloom/page pruning benefit |
| D | Parquet-only existing `lookup_variants` | compact `translation_sift/` | Broader baseline, not expected to win |

Repeat with:

```text
buffer_size = 1000, 5000, 20000
cold row group rows = 4k, 8k, 16k, 64k, 128k, 250k, 500k
with and without position_key Bloom filters
with and without page indexes/statistics
```

## Result Table Template

| Case | Layout | Buffer | Total s | variation_lookup s | warm hit % | cold positive probes | fjall gets | cold RG read | cold rows decoded | output diff |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| A warm+fjall | 500k | 5000 | benchmark output `time` | `VEP_PROFILE` | benchmark counter | benchmark counter | benchmark counter | 0 | 0 | baseline |
| B direct cold parquet | 500k | 5000 | benchmark output `time` | `VEP_PROFILE` | benchmark counter | benchmark counter | 0 for variation | benchmark counter | benchmark counter | none |
| C DataFusion cold parquet | 500k+bloom | 5000 | benchmark output `time` | `VEP_PROFILE` | benchmark counter | benchmark counter | 0 for variation | benchmark counter | benchmark counter | none |

## Optimization Notes

### Parquet Layout

- Smaller cold row groups are likely the main lever. Current 500k-row groups are good for scan throughput but expensive for sparse batched point lookup.
- Keep row groups position-aligned. Splitting the same `position_key` across row groups breaks definitive miss semantics and forces boundary fallbacks.
- Add an echtvar-style `2^20` genomic bin boundary constraint to the layout experiments. Echtvar's chunking is by genomic span, not row count; applying that to Parquet means row groups should not cross bin boundaries when testing cursor locality.
- Do not assume 4k/8k/16k row groups are always better. They reduce decoded candidate rows but increase row-group count and footer/page-index metadata. The benchmark must report metadata load/decode time separately.
- Store `position_key` as a first-class scalar column with page statistics and Bloom filters.
- Remove `variant_keys: List<UInt64>` from cold parquet entirely. It remains only in warm parquet, where warm chunk lookup needs it for the in-memory exact variant index.
- Build the existing `variation.variant_bloom_index` sidecar from temporary generated variant keys during cache generation, then drop the `variant_keys` column before writing cold parquet.
- Do not keep backward compatibility for old cold parquet schemas with `variant_keys`; fail fast and require rebuilding the indexed-parquet cache.
- Keep cold lookup projection narrow. The DataFusion blog calls projection pruning one of the simplest high-impact optimizations; the cold lookup should never read all 80 variation columns.
- Cache `ArrowReaderMetadata` per chromosome; repeated metadata decode is explicitly called out by the DataFusion pruning pipeline as latency-sensitive.

### Monotonic Cursor / Echtvar-Style Lookup

- VCF probes and cold Parquet are both position-sorted, so the reader should advance through row groups monotonically instead of binary searching from the start for every probe.
- Mirror echtvar's "load bin once, reuse while queries stay in the same bin" behavior with `probe_bin = position >> 20` and a bounded decoded cache for row groups in the active bin.
- Use a `ColdProbeBatch` sorted by `position_key`, then process row groups in order.
- Maintain current and previous row-group caches to handle deletion/shift probes that step backward by a small amount.
- Within a loaded row group, binary-search the first `position_key`, then advance a row cursor for all probes in that group.
- Use packed/stable variant-key generation only while building the `variation.variant_bloom_index` sidecar. Cold parquet lookup itself stays position-key driven and must verify alleles from `allele_string`.
- For large buffers, deduplicate probes before IO. Multiple VEP probe starts can map to the same `(position_key, ref, alt)`.

## Expected Outcome

The parquet approach should win if the number of cold-positive probes per buffer is large enough to amortize row-group reads and if the cold row-group layout is narrow enough. It will lose against Fjall when cold-positive probes are very sparse and each buffer touches many wide 500k-row groups. The first benchmark should therefore answer two questions separately:

1. Is cold Parquet faster than Fjall on the current 500k-row chr4 layout?
2. If not, does a 4k/8k/16k/64k/128k position-aligned cold layout with scalar `position_key` Bloom/page pruning become faster?
