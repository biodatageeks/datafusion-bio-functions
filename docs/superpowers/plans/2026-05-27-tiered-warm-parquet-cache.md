# Tiered Warm Parquet Cache Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in warm Parquet tier that uses echtvar-style keys and row-group chunks to avoid most Fjall point lookups while preserving Fjall fallback correctness.

**Architecture:** Parquet remains the master cache. Cache build writes `chrN_warm.parquet` and `chrN_cold.parquet`; runtime loads one projected warm row group per chromosome lane, builds a small key index, and falls back to Fjall for positions not covered by the warm chunk.

**Tech Stack:** Rust, DataFusion, Arrow, Parquet, Fjall, `HashMap`, `SmallVec`, existing VEPyR chromosome-lane VCF sink.

---

## File Structure

- Create: `datafusion/bio-function-vep/src/warm_cache/key.rs`
  - Deterministic `position_key` and `variant_keys` generation.
- Create: `datafusion/bio-function-vep/src/warm_cache/chunk.rs`
  - `WarmChunkContext`, Arrow projection handling, and chunk-local indexes.
- Create: `datafusion/bio-function-vep/src/warm_cache/split.rs`
  - Warm/cold split logic from chromosome Parquet files.
- Create: `datafusion/bio-function-vep/src/warm_cache/mod.rs`
  - Public module boundary.
- Modify: `datafusion/bio-function-vep/src/lib.rs`
  - Export `warm_cache`.
- Modify: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`
  - Probe warm chunk before Fjall fallback.
- Modify: `datafusion/bio-function-vep/src/cache_builder.rs`
  - Add opt-in warm/cold Parquet writer or route to the splitter module.
- Create: `datafusion/bio-function-vep/examples/build_warm_variation_cache.rs`
  - CLI for generating `_warm` and `_cold` files from existing per-chromosome Parquet.
- Create: `datafusion/bio-function-vep/examples/bench_warm_tier_chr1.rs`
  - Isolated benchmark: current Fjall lookup vs warm chunk + Fjall fallback.

## Task 1: Key Encoding

**Files:**
- Create: `datafusion/bio-function-vep/src/warm_cache/key.rs`
- Create: `datafusion/bio-function-vep/src/warm_cache/mod.rs`
- Modify: `datafusion/bio-function-vep/src/lib.rs`

- [ ] **Step 1: Add failing tests for position key ordering**

Add tests that prove keys sort by chromosome and start:

```rust
#[test]
fn position_key_orders_by_chrom_and_start() {
    let chr1_100 = position_key("1", 100);
    let chr1_200 = position_key("1", 200);
    let chr2_1 = position_key("2", 1);
    assert!(chr1_100 < chr1_200);
    assert!(chr1_200 < chr2_1);
}
```

- [ ] **Step 2: Add failing tests for multi-allelic variant keys**

```rust
#[test]
fn variant_keys_include_each_alt_from_allele_string() {
    let keys = variant_keys_from_allele_string("1", 101, "A/G/T").unwrap();
    assert_eq!(keys.len(), 2);
    assert_ne!(keys[0], keys[1]);
}
```

- [ ] **Step 3: Implement `position_key`**

Use the same chromosome code semantics as `kv_cache::key_encoding::chrom_to_code`.

```rust
pub fn position_key(chrom: &str, start: i64) -> u64 {
    let chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let chrom_code = crate::kv_cache::key_encoding::chrom_to_code(chrom) as u64;
    let start = u64::try_from(start).expect("start must be non-negative");
    (chrom_code << 48) | (start & 0x0000_FFFF_FFFF_FFFF)
}
```

- [ ] **Step 4: Implement `variant_keys_from_allele_string`**

Generate one key per alternate allele. For MVP, use a stable deterministic hash over `(position_key, ref, alt)` and always verify `allele_string` after lookup.

```rust
pub fn variant_keys_from_allele_string(
    chrom: &str,
    start: i64,
    allele_string: &str,
) -> Result<Vec<u64>> {
    let Some((reference, alts)) = allele_string.split_once('/') else {
        return Ok(Vec::new());
    };
    let pos = position_key(chrom, start);
    Ok(alts
        .split('/')
        .filter(|alt| !alt.is_empty())
        .map(|alt| stable_variant_hash(pos, reference, alt))
        .collect())
}
```

- [ ] **Step 5: Run tests**

Run:

```bash
cargo test -p datafusion-bio-function-vep --features kv-cache warm_cache::key
```

Expected: tests pass.

## Task 2: Warm/Cold Splitter

**Files:**
- Create: `datafusion/bio-function-vep/src/warm_cache/split.rs`
- Create: `datafusion/bio-function-vep/examples/build_warm_variation_cache.rs`

- [ ] **Step 1: Add tests for AF parsing**

Cover scalar and allele-pair string forms:

```rust
#[test]
fn max_af_parses_pair_lists() {
    assert_eq!(max_af_from_pairs(Some("A:0.1,C:2.6e-05")), 0.1);
    assert_eq!(max_af_from_pairs(Some("A:0,C:0")), 0.0);
    assert_eq!(max_af_from_pairs(None), 0.0);
}
```

- [ ] **Step 2: Implement warm position selection**

Given rows with `start` and `max_global_af`, select complete positions:

```text
hot = pos_max_af >= 0.01
warm = hot - 1, hot, hot + 1 if those positions exist in the source chromosome Parquet
```

- [ ] **Step 3: Write `_warm` and `_cold` Parquet files**

Output:

```text
variation/chrN_warm.parquet
variation/chrN_cold.parquet
```

Both files must include generated columns:

```text
position_key: UInt64
variant_keys: List<UInt64>
```

Use row-group size from CLI:

```text
--row-group-rows 250000
```

- [ ] **Step 4: Add metadata**

Write Parquet file metadata:

```text
vepyr.cache_tier = warm|cold
vepyr.warm_selector = max_global_af>=0.01,+/-1
vepyr.key_version = 1
vepyr.row_group_rows = 250000
```

- [ ] **Step 5: Run chr1 splitter smoke test**

Run:

```bash
RUSTFLAGS="-C target-cpu=native" cargo run --release \
  -p datafusion-bio-function-vep \
  --features kv-cache \
  --example build_warm_variation_cache -- \
  --input /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation/chr1.parquet \
  --output-dir /tmp/vepyr-warm-prototype/variation \
  --af-threshold 0.01 \
  --position-radius 1 \
  --row-group-rows 250000
```

Expected:

```text
chr1_warm.parquet rows ~= 3.6M
chr1_warm.parquet size ~= 275-300 MiB
chr1_cold.parquet rows + warm rows = source rows
```

## Task 3: Warm Chunk Context

**Files:**
- Create: `datafusion/bio-function-vep/src/warm_cache/chunk.rs`

- [ ] **Step 1: Add tests for chunk indexing**

Create a small `RecordBatch` with:

```text
position_key: [101, 101, 205]
variant_keys: [[k1], [k2], [k3]]
allele_string: ["A/G", "A/T", "C/A"]
```

Assert:

```rust
assert!(chunk.contains_position(101));
assert_eq!(chunk.lookup_variant(k2).len(), 1);
assert!(chunk.lookup_variant(missing).is_empty());
```

- [ ] **Step 2: Implement `WarmChunkContext`**

```rust
pub struct WarmChunkContext {
    pub row_group_id: usize,
    pub min_position_key: u64,
    pub max_position_key: u64,
    pub position_keys: Vec<u64>,
    pub variant_index: HashMap<u64, SmallVec<[u32; 1]>>,
    pub batch: RecordBatch,
    pub columns: WarmColumnIndices,
}
```

- [ ] **Step 3: Implement index build**

Read `position_key` and `variant_keys` from the loaded row group. Deduplicate and sort `position_keys`; add one `variant_index` entry per key in the list column.

- [ ] **Step 4: Implement chunk lookup primitives**

```rust
pub enum WarmProbeResult {
    Hit(SmallVec<[u32; 1]>),
    DefinitiveMiss,
    NotCovered,
}
```

Logic:

```text
variant key hit -> Hit(rows)
variant key miss + position present -> DefinitiveMiss
position absent/out of chunk range -> NotCovered
```

## Task 4: Warm Chunk Reader

**Files:**
- Create: `datafusion/bio-function-vep/src/warm_cache/reader.rs`
- Modify: `datafusion/bio-function-vep/src/warm_cache/mod.rs`

- [ ] **Step 1: Add row-group range tests**

Given fake row-group stats:

```text
rg0 min=100 max=200
rg1 min=201 max=300
```

Assert that a probe range `150..250` selects `rg0` and `rg1`.

- [ ] **Step 2: Implement metadata open**

Read `chrN_warm.parquet` metadata and collect each row group's min/max `position_key`.

- [ ] **Step 3: Implement projected row-group load**

Projection must include only runtime columns:

```text
position_key, variant_keys, chrom, start, end, allele_string, variation_name,
failed, somatic, strand, minor_allele, minor_allele_freq, clin_sig,
phenotype_or_disease, clinical_impact, clin_sig_allele, AF, gnomADg, gnomADe
```

- [ ] **Step 4: Keep at most two active chunks**

The reader should keep the current row group and the next row group only when a VCF buffer crosses a boundary.

## Task 5: Tiered Lookup Integration

**Files:**
- Modify: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`
- Modify: `datafusion/bio-function-vep/src/lookup_provider.rs`

- [ ] **Step 1: Add opt-in configuration**

Use environment flag first:

```text
VEP_WARM_VARIATION_CACHE=1
```

Auto-detect warm files:

```text
{cache_source}/variation/chrN_warm.parquet
```

- [ ] **Step 2: Route probes through warm chunk before Fjall**

For every probe currently sent to Fjall:

```text
try warm chunk
  Hit -> verify allele_string and append output values
  DefinitiveMiss -> skip Fjall for this probe
  NotCovered -> existing Fjall path
```

- [ ] **Step 3: Preserve current Fjall fallback**

Do not remove existing Fjall code. The warm path must be bypassable at runtime.

- [ ] **Step 4: Add tests for fallback semantics**

Cases:

```text
variant_key hit -> no Fjall call
variant_key miss + position_key present -> no Fjall call
variant_key miss + position_key absent -> Fjall call
warm disabled -> Fjall call
```

Use a fake Fjall lookup counter in unit tests.

## Task 6: Profiling Counters

**Files:**
- Modify: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`

- [ ] **Step 1: Extend profile struct**

Add:

```text
warm_chunks_loaded
warm_chunk_rows
warm_chunk_load_s
warm_chunk_index_s
warm_variant_key_hits
warm_variant_key_misses
warm_verified_hits
warm_definitive_misses
warm_position_misses
fjall_fallbacks
fjall_fallback_saved
```

- [ ] **Step 2: Emit profile line**

Under `VEP_KV_PROFILE=1`, print:

```text
[vep-warm-profile] chunks=... rows=... load_s=... index_s=... verified_hits=... definitive_misses=... fjall_fallbacks=...
```

- [ ] **Step 3: Test profile formatting**

Add unit test matching the existing `lookup_profile_detail_line_formats_probe_decode_match_breakdown` style.

## Task 7: Benchmarks

**Files:**
- Create: `datafusion/bio-function-vep/examples/bench_warm_tier_chr1.rs`
- Add logs under `/tmp`, not repository.

- [ ] **Step 1: Build native release**

Run:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release \
  -p datafusion-bio-function-vep \
  --features kv-cache \
  --example bench_annotate_vcf \
  --example bench_warm_tier_chr1
```

- [ ] **Step 2: Baseline current Fjall full run**

Run:

```bash
VEP_PROFILE=1 VEP_KV_PROFILE=1 VEP_KV_RANGE_PREFETCH=0 \
target/release/examples/bench_annotate_vcf \
  --input /Users/mwiewior/workspace/data_vepyr/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --cache /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --output /tmp/vepyr_baseline_fjall.vcf \
  --backend fjall \
  --everything \
  --extended-probes \
  --reference-fasta /Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --forks 6 \
  --contig-parallelism 1 \
  --buffer-size 5000 \
  --no-progress \
  2>&1 | tee /tmp/vepyr_baseline_fjall.log
```

Record:

```text
wall time
rows/s
point_gets
max RSS
```

- [ ] **Step 3: Warm tier full run**

Run the same command with:

```text
VEP_WARM_VARIATION_CACHE=1
```

Expected counters:

```text
fjall_fallbacks ~= 2.8M
fjall_fallback_saved ~= 6.2M
warm verified hits + definitive misses covers most saved probes
```

- [ ] **Step 4: Chunk-size sweep**

Generate warm files with:

```text
--row-group-rows 100000
--row-group-rows 250000
--row-group-rows 500000
```

Benchmark chr1 and full VCF. Report:

```text
chunk rows
wall time
warm chunk load/index time
max RSS
fjall fallback count
```

- [ ] **Step 5: Correctness validation**

Run VEPyR e2e comparison for:

```text
merged
ensembl
refseq
```

Acceptance:

```text
0 mismatches for all caches
```

## Acceptance Criteria

- Warm/cold split writes valid Parquet files for all chromosomes.
- Full benchmark keeps warm rows near 4.2% and warm size near 9.7% of variation Parquet.
- Warm path avoids about 65-70% of Fjall point lookups on the full benchmark VCF.
- Full VCF annotation is at least 15% faster than Fjall-only baseline.
- Stretch target: full VCF annotation reaches 35-45 seconds on the local native build.
- Max RSS stays under 6 GB for the chosen production-like settings.
- Merged, Ensembl, and RefSeq comparisons report 0 mismatches.

## Self-Review

- Spec coverage: key design, warm/cold split, row-group chunking, runtime probing, fallback behavior, metrics, and benchmarking are covered.
- Placeholder scan: no placeholder implementation steps remain.
- Type consistency: plan consistently uses `position_key: u64`, `variant_keys: List<UInt64>`, `WarmChunkContext`, `WarmProbeResult`, and the `VEP_WARM_VARIATION_CACHE` opt-in flag.
