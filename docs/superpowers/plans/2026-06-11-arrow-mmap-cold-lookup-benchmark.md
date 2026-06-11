# Arrow Mmap Cold Lookup Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Rust benchmark that compares Lance cold lookup with a compact LMDB position-range index plus mmap-opened Arrow IPC batch files for the chr1 VEPyr cold workload.

**Architecture:** Add a benchmark-only example that builds three Arrow stores for chr1 cold rows: uncompressed IPC, LZ4 IPC, and ZSTD IPC. Each store has a `heed` LMDB index mapping `position_key` to one contiguous `(batch_id, start_row, len)` range and Arrow IPC batch files containing all projected `everything` columns.

**Tech Stack:** Rust, Arrow 58 IPC, parquet Arrow reader, Lance 7.0.0, `heed` LMDB, `memmap2`, `serde_json`, existing VEP position and bloom indexes.

---

### Task 1: Cargo Wiring

**Files:**
- Modify: `datafusion/bio-function-vep/Cargo.toml`

- [ ] **Step 1: Add dev dependencies**

Add these dev dependencies to `datafusion/bio-function-vep/Cargo.toml`:

```toml
heed = "0.22"
memmap2 = "0.9"
```

If Cargo resolves a newer compatible `heed` version and the API differs, pin to the resolved version that compiles with `EnvOpenOptions`, `Database`, `types::Bytes`, and `types::U64`.

- [ ] **Step 2: Register the benchmark example**

Add the example entry:

```toml
[[example]]
name = "bench_arrow_mmap_cold_lookup"
required-features = ["lance-cache"]
```

- [ ] **Step 3: Verify Cargo metadata**

Run:

```bash
cargo metadata -p datafusion-bio-function-vep --format-version 1 >/tmp/arrow_mmap_cargo_metadata.json
```

Expected: command exits successfully.

### Task 2: Example Skeleton And Shared Workload

**Files:**
- Create: `datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs`

- [ ] **Step 1: Create the example with CLI parsing and shared workload helpers**

The example accepts:

```text
bench_arrow_mmap_cold_lookup <cache_root> <input.vcf.gz> <merged_lance_chr_dataset> <arrow_store_root> [report.md]
```

It reads these environment variables:

```text
CHROM=chr1
ARROW_MMAP_BATCH_ROWS=16384
COLD_SCAN_BATCH_SIZE=2000
SMALL_QUERY_BATCH_SIZE=238
LARGE_QUERY_BATCH_SIZE=2000
SORT_COLD_POSITIONS=1
REBUILD_ARROW_MMAP=0
```

Copy the workload construction from `bench_lance_cold_batching.rs` into local helpers:

```rust
fn load_vcf_workload(path: &Path, chrom: &str) -> BenchResult<Workload>
fn build_gate(...)
fn projected_everything_columns(dataset: &Dataset) -> Vec<String>
```

Keep these copied helpers local to the example to avoid changing production library APIs.

- [ ] **Step 2: Add empty smoke test**

Add a tiny test that only validates `PositionRange` encode/decode once Task 3 exists:

```rust
#[test]
fn position_range_round_trips_bytes() {
    let range = PositionRange {
        batch_id: 7,
        start_row: 11,
        len: 13,
    };
    assert_eq!(PositionRange::decode(&range.encode()).unwrap(), range);
}
```

- [ ] **Step 3: Run the example test and expect compile failure before Task 3**

Run:

```bash
cargo test -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup position_range_round_trips_bytes
```

Expected before Task 3: compile fails because `PositionRange` is not implemented.

### Task 3: Position Range Encoding And Store Metadata

**Files:**
- Modify: `datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs`

- [ ] **Step 1: Implement `PositionRange`**

Add:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PositionRange {
    batch_id: u32,
    start_row: u32,
    len: u32,
}

impl PositionRange {
    const BYTE_LEN: usize = 12;

    fn encode(self) -> [u8; Self::BYTE_LEN] {
        let mut out = [0_u8; Self::BYTE_LEN];
        out[0..4].copy_from_slice(&self.batch_id.to_le_bytes());
        out[4..8].copy_from_slice(&self.start_row.to_le_bytes());
        out[8..12].copy_from_slice(&self.len.to_le_bytes());
        out
    }

    fn decode(bytes: &[u8]) -> BenchResult<Self> {
        if bytes.len() != Self::BYTE_LEN {
            return Err(format!("invalid PositionRange byte length: {}", bytes.len()).into());
        }
        Ok(Self {
            batch_id: u32::from_le_bytes(bytes[0..4].try_into()?),
            start_row: u32::from_le_bytes(bytes[4..8].try_into()?),
            len: u32::from_le_bytes(bytes[8..12].try_into()?),
        })
    }
}
```

- [ ] **Step 2: Implement compression and manifest structs**

Add:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArrowStoreCompression {
    Uncompressed,
    Lz4,
    Zstd,
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct ArrowMmapManifest {
    format: String,
    version: u32,
    chrom: String,
    tier: String,
    compression: String,
    projected_columns: Vec<String>,
    target_batch_rows: usize,
    batch_count: usize,
    total_rows: usize,
    unique_positions: usize,
    max_position_group_len: usize,
    data_file_size: u64,
    index_size: u64,
    total_size: u64,
}
```

- [ ] **Step 3: Run the range test**

Run:

```bash
cargo test -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup position_range_round_trips_bytes
```

Expected: test passes.

### Task 4: Arrow Store Builder

**Files:**
- Modify: `datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs`

- [ ] **Step 1: Implement cold parquet file discovery**

Use existing cache layout conventions to find chr cold files under:

```text
<cache_root>/variation/
```

Support both:

```text
chr1_cold.parquet
1_cold.parquet
chr1_cold_*.parquet
1_cold_*.parquet
```

- [ ] **Step 2: Implement projected parquet reads**

Use `ParquetRecordBatchReaderBuilder` and `projection_for_existing_roots`-style logic locally:

```rust
fn projection_mask_for_columns(
    arrow_schema: &Schema,
    parquet_schema: &parquet::schema::types::SchemaDescriptor,
    columns: &[String],
) -> BenchResult<ProjectionMask>
```

Fail if any projected column is missing.

- [ ] **Step 3: Implement position-group buffering**

Read projected batches in order. Accumulate full row groups with the same `position_key` and flush groups without splitting a position across batch files.

If a later `position_key` is lower than the previous key, return:

```text
rows out of position_key order
```

- [ ] **Step 4: Implement Arrow IPC writer**

Write one IPC file per output batch:

```rust
fn write_arrow_batch_file(
    path: &Path,
    batch: &RecordBatch,
    compression: ArrowStoreCompression,
) -> BenchResult<u64>
```

Use Arrow IPC file format. Use no compression for `Uncompressed`, LZ4 for `Lz4`, and ZSTD for `Zstd`.

- [ ] **Step 5: Implement `heed` index writes**

Open an LMDB environment in `<store>/index.mdb`, create database `position_ranges`, and insert one key/value per unique `position_key`.

Use big-endian or little-endian consistently for `u64` keys. Values use `PositionRange::encode()`.

- [ ] **Step 6: Add builder tests**

Create in-example tests that build a tiny store from synthetic batches and assert:

```rust
#[test]
fn builder_keeps_position_group_inside_one_batch_file()

#[test]
fn builder_rejects_out_of_order_position_key()
```

- [ ] **Step 7: Run builder tests**

Run:

```bash
cargo test -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup builder_
```

Expected: both builder tests pass.

### Task 5: Arrow Store Reader And Column Touch

**Files:**
- Modify: `datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs`

- [ ] **Step 1: Implement store open**

Load `manifest.json`, open LMDB read-only, and lazily mmap Arrow batch files on demand.

Cache decoded `RecordBatch` values by `batch_id` inside the benchmark reader. This keeps the first lookup honest while avoiding repeated IPC decoding for the same file.

- [ ] **Step 2: Implement lookup**

For each position:

```rust
fn lookup_position(&mut self, position_key: u64) -> BenchResult<Option<PositionRange>>
```

Count one index get per requested position and one missing position per absent key.

- [ ] **Step 3: Implement deterministic column touching**

Implement:

```rust
fn touch_rows(batch: &RecordBatch, range: PositionRange, acc: &mut TouchAccumulator) -> BenchResult<()>
```

Handle numeric primitive arrays, boolean, string, large string, binary, large binary, list, large list, struct, dictionary, and null arrays. For unknown types, include validity and array length in the checksum and return an explicit error only if the type cannot be safely traversed.

- [ ] **Step 4: Add reader tests**

Add:

```rust
#[test]
fn arrow_variants_return_same_checksum()

#[test]
fn missing_position_is_counted_as_miss()
```

- [ ] **Step 5: Run reader tests**

Run:

```bash
cargo test -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup arrow_variants_return_same_checksum missing_position_is_counted_as_miss
```

Expected: tests pass.

### Task 6: Lance Comparison And Report

**Files:**
- Modify: `datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs`

- [ ] **Step 1: Reuse Lance scan comparison**

Port `scan_plan` and `ScanIoStats` from `bench_lance_cold_batching.rs` into this example. Keep `small_chunks` and `large_chunks` Lance plans.

- [ ] **Step 2: Implement Arrow benchmark result rows**

For each store variant, report:

```text
format
compression
total_size
data_size
index_size
build_s
open_s
lookup_s
positions_requested
positions_found
selected_rows
projected_columns
estimated_bytes_touched
batch_files_opened
index_gets
missing_positions
checksum
```

- [ ] **Step 3: Implement markdown report**

Report sections:

```text
# Arrow Mmap Cold Lookup Benchmark
## Workload
## Store Sizes
## Lance Cold Lookup
## Arrow Mmap Lookup
## Verification
```

The verification section must state whether all Arrow variants agree and whether Lance selected rows/positions match Arrow.

- [ ] **Step 4: Run compile check**

Run:

```bash
RUSTFLAGS="-C target-cpu=native" cargo check -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup
```

Expected: command passes.

### Task 7: Chr1 Benchmark Run

**Files:**
- Output: `research/reports/arrow_mmap_cold_lookup_20260611/`

- [ ] **Step 1: Run the benchmark on chr1**

Run:

```bash
mkdir -p research/reports/arrow_mmap_cold_lookup_20260611
RUSTFLAGS="-C target-cpu=native" \
REBUILD_ARROW_MMAP=1 \
SORT_COLD_POSITIONS=1 \
cargo run --release -p datafusion-bio-function-vep --features lance-cache \
  --example bench_arrow_mmap_cold_lookup -- \
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  /Users/mwiewior/research/git/vepyr/e2e-testing/results/fast_chr1/input_chr1.vcf.gz \
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/chr1.lance \
  /Users/mwiewior/workspace/data_vepyr/arrow_mmap_chr1_bench \
  research/reports/arrow_mmap_cold_lookup_20260611/chr1_arrow_mmap_cold_lookup.md \
  2>&1 | tee research/reports/arrow_mmap_cold_lookup_20260611/chr1_arrow_mmap_cold_lookup.log
```

Expected: report is written and includes all four formats.

- [ ] **Step 2: Compare rows and checksums**

Open the report and verify:

```text
Arrow uncompressed, LZ4, and ZSTD selected rows match.
Arrow uncompressed, LZ4, and ZSTD checksums match.
Any Lance mismatch is explicitly reported.
```

### Task 8: Final Verification And Commit

**Files:**
- Modify: `datafusion/bio-function-vep/Cargo.toml`
- Create: `datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs`
- Output: `research/reports/arrow_mmap_cold_lookup_20260611/`

- [ ] **Step 1: Run focused tests**

Run:

```bash
cargo test -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup
```

Expected: all example tests pass.

- [ ] **Step 2: Run compile check**

Run:

```bash
RUSTFLAGS="-C target-cpu=native" cargo check -p datafusion-bio-function-vep --features lance-cache --example bench_arrow_mmap_cold_lookup
```

Expected: check passes.

- [ ] **Step 3: Commit implementation**

Stage only the benchmark implementation, Cargo changes, and generated report:

```bash
git add datafusion/bio-function-vep/Cargo.toml datafusion/bio-function-vep/examples/bench_arrow_mmap_cold_lookup.rs research/reports/arrow_mmap_cold_lookup_20260611
git commit -m "bench: add arrow mmap cold lookup comparison"
```

Expected: commit succeeds without staging unrelated existing dirty files.
