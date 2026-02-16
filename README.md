# datafusion-bio-functions

[![CI](https://github.com/biodatageeks/datafusion-bio-functions/actions/workflows/ci.yml/badge.svg)](https://github.com/biodatageeks/datafusion-bio-functions/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

Apache DataFusion user-defined functions for bioinformatics, enabling SQL-based genomic analysis.

## Overview

This workspace provides a collection of Rust crates that implement DataFusion UDFs (User-Defined Functions) for bioinformatics. It is the companion to [datafusion-bio-formats](https://github.com/biodatageeks/datafusion-bio-formats) (which handles bioinformatic file format readers). Both are consumed by [polars-bio](https://github.com/biodatageeks/polars-bio).

## Crates

| Crate | Description | Status |
|-------|-------------|--------|
| **[datafusion-bio-function-pileup](datafusion/bio-function-pileup)** | Depth-of-coverage (pileup) computation from BAM alignments | ✅ |

## Features

- **Depth-of-Coverage**: Compute per-base depth from BAM alignment data using an efficient event-based algorithm
- **SQL Interface**: Query coverage via SQL table function: `SELECT * FROM depth('file.bam')`
- **Per-Base Output**: Optional one-row-per-position mode (like `samtools depth -a`) with lazy streaming to keep memory at O(batch_size)
- **DataFusion Integration**: Native `ExecutionPlan` implementation with partition-parallel execution
- **Mosdepth-Compatible**: Fast-mode behavior matching [mosdepth](https://github.com/brentp/mosdepth) (no mate-pair overlap deduplication)
- **Dense Accumulation**: Mosdepth-style flat `i32[]` depth arrays with O(1) writes and streaming per-contig emission

## Installation

Add the crates you need to your `Cargo.toml`:

```toml
[dependencies]
datafusion = "50.3.0"
datafusion-bio-function-pileup = { git = "https://github.com/biodatageeks/datafusion-bio-functions" }
```

## Quick Start

### Coverage via SQL Table Function

The simplest way to compute coverage is via the SQL `depth()` table function:

```rust
use datafusion::prelude::*;
use datafusion_bio_function_pileup::register_pileup_functions;

#[tokio::main]
async fn main() -> datafusion::error::Result<()> {
    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    // 1-based coordinates (default, matches polars-bio)
    let df = ctx.sql("SELECT * FROM depth('path/to/alignments.bam')").await?;
    df.show().await?;

    // 0-based coordinates (explicit)
    let df = ctx.sql("SELECT * FROM depth('path/to/alignments.bam', true)").await?;
    df.show().await?;

    // Per-base output: one row per genomic position (like samtools depth -a)
    let df = ctx.sql("SELECT * FROM depth('path/to/alignments.bam', true, true)").await?;
    df.show().await?;

    Ok(())
}
```

The optional arguments control behavior:

```sql
SELECT * FROM depth('file.bam')                -- 1-based, block output (default)
SELECT * FROM depth('file.bam', true)          -- 0-based, block output
SELECT * FROM depth('file.bam', false)         -- 1-based, block output (explicit)
SELECT * FROM depth('file.bam', true, true)    -- 0-based, per-base output
SELECT * FROM depth('file.bam', false, true)   -- 1-based, per-base output
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `path` | string | (required) | Path to BAM file |
| `zero_based` | bool | `false` | Use 0-based coordinates |
| `per_base` | bool | `false` | Emit one row per position instead of RLE blocks |

#### Block output schema (default)

| Column | Type | Description |
|--------|------|-------------|
| `contig` | Utf8 | Chromosome/contig name |
| `pos_start` | Int32 | Block start position (inclusive) |
| `pos_end` | Int32 | Block end position (inclusive) |
| `coverage` | Int16 | Read depth in this block |

#### Per-base output schema (`per_base=true`)

| Column | Type | Description |
|--------|------|-------------|
| `contig` | Utf8 | Chromosome/contig name |
| `pos` | Int32 | Genomic position |
| `coverage` | Int16 | Read depth at this position |

Per-base mode emits one row for every genomic position in each contig (including 0-coverage positions), similar to `samtools depth -a`. It requires dense mode (BAM header with contig lengths) and uses a lazy `PerBaseEmitter` that streams positions in batch-sized chunks, keeping memory at O(batch_size).

Both schemas include metadata key `bio.coordinate_system_zero_based` (`"true"` or `"false"`) so downstream consumers can interpret the coordinate system.

### Coverage via Programmatic API

For more control, use `PileupExec` directly:

```rust
use std::sync::Arc;
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use datafusion_bio_function_pileup::physical_exec::{PileupConfig, PileupExec};
use futures::StreamExt;

#[tokio::main]
async fn main() -> datafusion::error::Result<()> {
    // Register BAM file as a table (zero_based=false for 1-based coordinates)
    let table = BamTableProvider::new("path/to/alignments.bam", None, false, None).await?;
    let ctx = SessionContext::new();
    ctx.register_table("reads", Arc::new(table))?;

    // Create physical plan and wrap with PileupExec
    let df = ctx.table("reads").await?;
    let plan = df.create_physical_plan().await?;
    let pileup = PileupExec::new(plan, PileupConfig::default()); // zero_based=false by default

    // Execute and collect results
    let task_ctx = ctx.task_ctx();
    let stream = pileup.execute(0, task_ctx)?;
    let batches: Vec<_> = stream.collect::<Vec<_>>().await
        .into_iter()
        .filter_map(|r| r.ok())
        .collect();

    for batch in &batches {
        println!("{:?}", batch);
    }
    Ok(())
}
```

### Filtering

Coverage computation automatically filters reads using SAM flags (default: exclude unmapped, secondary, failed QC, duplicate reads) and minimum mapping quality:

```rust
use datafusion_bio_function_pileup::ReadFilter;
use datafusion_bio_function_pileup::physical_exec::{PileupConfig, PileupExec};

let config = PileupConfig {
    filter: ReadFilter {
        filter_flag: 1796,        // unmapped | secondary | failed_qc | duplicate
        min_mapping_quality: 20,  // minimum MAPQ
    },
    ..PileupConfig::default()     // zero_based=false, dense_mode=Auto
};
let pileup = PileupExec::new(plan, config);
```

## Downstream Library Integration

This crate is designed to be consumed by downstream libraries such as [polars-bio](https://github.com/biodatageeks/polars-bio). This section documents the key APIs and configuration options.

### Dependency Setup

Pin to a specific commit to ensure reproducible builds:

```toml
[dependencies]
datafusion-bio-function-pileup = { git = "https://github.com/biodatageeks/datafusion-bio-functions", rev = "<commit-hash>" }
datafusion-bio-format-bam = { git = "https://github.com/biodatageeks/datafusion-bio-formats", rev = "<commit-hash>" }
```

Keep `datafusion-bio-formats` and `datafusion-bio-functions` revisions in sync — they must use the same DataFusion/Arrow versions (currently DataFusion 50.3.0, Arrow 56.x).

### Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `bam` | yes | Enables the `table_function` module (`DepthFunction`, `register_pileup_functions`) which depends on `datafusion-bio-format-bam` |

Downstream libraries that provide their own `TableProvider` (e.g., polars-bio) can disable the default feature to avoid pulling in `datafusion-bio-format-bam` as a runtime dependency:

```toml
[dependencies]
datafusion-bio-function-pileup = { git = "...", rev = "...", default-features = false }
```

The core pileup API (`PileupExec`, `PileupConfig`, `cigar::*`, `events::*`, `coverage::*`) is always available regardless of features.

### BamTableProvider API

```rust
BamTableProvider::new(
    file_path,           // String: path to BAM file (local or object-store URL)
    storage_opts,        // Option<ObjectStorageOptions>: S3/GCS/Azure credentials
    zero_based,          // bool: false for 1-based (default), true for 0-based
    tag_fields,          // Option<Vec<String>>: optional BAM tag fields to include
).await?
```

The pileup crate auto-detects the CIGAR format from the schema (`DataType::Binary` for binary CIGAR, `DataType::Utf8` for string).

### PileupConfig Defaults

`PileupConfig::default()` uses the fastest settings:

| Field | Default | Description |
|-------|---------|-------------|
| `zero_based` | `false` | 1-based coordinates (matches polars-bio default) |
| `per_base` | `false` | RLE block output; set `true` for one-row-per-position |
| `binary_cigar` | `false` | Auto-detected from schema; `true` when upstream provides binary CIGAR |
| `dense_mode` | `Auto` (dense) | Dense `i32[]` array accumulation with streaming emission |
| `filter.filter_flag` | `1796` | Exclude unmapped, secondary, failed QC, duplicate |
| `filter.min_mapping_quality` | `0` | No MAPQ threshold |

To override specific fields:

```rust
let config = PileupConfig {
    zero_based: true,                       // 0-based coordinates
    per_base: true,                         // one row per position (requires dense mode)
    binary_cigar: false,                    // force string CIGAR (e.g., for debugging)
    dense_mode: DenseMode::Disable,         // force sparse BTreeMap accumulation
    filter: ReadFilter {
        min_mapping_quality: 20,
        ..ReadFilter::default()
    },
    ..PileupConfig::default()
};
```

### Accumulation Strategies

| Strategy | Best For | Memory | How |
|----------|----------|--------|-----|
| **Dense** (default) | WGS, high-coverage | One contig at a time | Flat `i32[]` indexed by position |
| **Sparse** | WES, targeted panels | Proportional to read count | `BTreeMap<pos, delta>` per contig |

Dense mode requires BAM header metadata (contig lengths) in the schema. If unavailable (e.g., `MemTable` in tests), it falls back to sparse automatically.

### End-to-End Example (polars-bio style)

```rust
use std::sync::Arc;
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use datafusion_bio_function_pileup::physical_exec::{PileupConfig, PileupExec};
use futures::StreamExt;

async fn compute_coverage(bam_path: &str, partitions: usize) -> Vec<RecordBatch> {
    let table = BamTableProvider::new(
        bam_path.to_string(), None, false, None,  // zero_based=false
    ).await.unwrap();

    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_table("reads", Arc::new(table)).unwrap();

    // Select only columns needed for coverage (avoids decoding sequence, quality, tags)
    let df = ctx.sql("SELECT chrom, start, flags, cigar, mapping_quality FROM reads")
        .await.unwrap();
    let plan = df.create_physical_plan().await.unwrap();

    let pileup = Arc::new(PileupExec::new(plan, PileupConfig::default()));
    let task_ctx = ctx.task_ctx();

    // Run all partitions concurrently
    let mut handles = Vec::new();
    for partition in 0..pileup.properties().partitioning.partition_count() {
        let pileup = pileup.clone();
        let task_ctx = task_ctx.clone();
        handles.push(tokio::spawn(async move {
            pileup.execute(partition, task_ctx).unwrap()
                .collect::<Vec<_>>().await
                .into_iter()
                .filter_map(|r| r.ok())
                .collect::<Vec<_>>()
        }));
    }

    let mut all_batches = Vec::new();
    for handle in handles {
        all_batches.extend(handle.await.unwrap());
    }
    all_batches
}
```

## Algorithm

The coverage algorithm is ported from [SeQuiLa](https://github.com/biodatageeks/sequila) and uses an event-based approach:

1. **CIGAR Processing**: Each aligned read's CIGAR generates +1 (start) and -1 (end) events at reference positions. Both string and binary CIGAR formats are supported (auto-detected from the Arrow schema).
2. **Dense Accumulation** (default): Events are written to a flat `i32[]` array per contig (O(1) writes). Completed contigs are streamed out as soon as the BAM's coordinate-sorted order moves to the next contig — peak memory is one contig at a time.
3. **Output Generation**:
   - **Block mode** (default): A sweep-line pass computes running coverage and emits run-length encoded blocks where coverage changes. Zero-coverage gaps are skipped.
   - **Per-base mode** (`per_base=true`): A lazy `PerBaseEmitter` walks the delta array with a prefix sum, emitting one row per genomic position (including 0-coverage). Batches are streamed one at a time to keep memory at O(batch_size).

The dense path naturally handles complex CIGAR operations (insertions, deletions, skipped regions) and provides excellent cache locality. A sparse `BTreeMap` fallback is available for targeted sequencing panels where contig-level arrays would be wasteful.

## Development

### Build

```bash
cargo build
```

### Test

```bash
cargo test --all
```

### Format

```bash
cargo fmt --all
```

### Lint

```bash
cargo clippy --all-targets --all-features -- -D warnings
```

### Documentation

```bash
cargo doc --no-deps --all-features --open
```

## Requirements

- Rust 1.88.0 or later (specified in `rust-toolchain.toml`)
- Git (for building crates with git dependencies)

## Related Projects

- [datafusion-bio-formats](https://github.com/biodatageeks/datafusion-bio-formats) - DataFusion table providers for bioinformatics file formats
- [polars-bio](https://github.com/biodatageeks/polars-bio) - Polars extensions for bioinformatics
- [SeQuiLa](https://github.com/biodatageeks/sequila) - Original Scala implementation of the pileup algorithm

## License

Licensed under Apache-2.0. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Acknowledgments

This project builds upon:
- [Apache DataFusion](https://datafusion.apache.org/) - Fast, extensible query engine
- [mosdepth](https://github.com/brentp/mosdepth) - Test data and reference implementation for coverage validation
- [SeQuiLa](https://github.com/biodatageeks/sequila) - Original pileup algorithm design
