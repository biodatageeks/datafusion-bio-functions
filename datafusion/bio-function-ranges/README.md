# datafusion-bio-function-ranges

Interval join, coverage, count-overlaps, and nearest-neighbor operations for Apache DataFusion.

This crate provides optimized genomic interval operations as DataFusion extensions:

- **Interval join** — SQL joins with range overlap conditions, optimized via multiple interval tree algorithms
- **Coverage** — base-pair overlap depth between two interval sets
- **Count overlaps** — number of overlapping intervals per region
- **Nearest** — nearest-neighbor interval matching

## Quick Start

```rust
use datafusion_bio_function_ranges::{create_bio_session, register_ranges_functions};

// Option 1: Create a fully configured session (recommended)
let ctx = create_bio_session();

// Option 2: Register functions on an existing bio-configured session
use datafusion::config::ConfigOptions;
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_ranges::{BioConfig, BioSessionExt};

let config = SessionConfig::from(ConfigOptions::new())
    .with_option_extension(BioConfig::default())
    .with_information_schema(true)
    .with_repartition_joins(false);
let ctx = SessionContext::new_with_bio(config);
register_ranges_functions(&ctx);
```

## Registering Functions

### `register_ranges_functions(ctx)`

Registers the `coverage`, `count_overlaps`, and `nearest` SQL table functions on an existing `SessionContext`. This is analogous to `register_pileup_functions` in the pileup crate.

```rust
use datafusion_bio_function_ranges::register_ranges_functions;

register_ranges_functions(&ctx);
```

### `create_bio_session()`

Convenience function that creates a `SessionContext` with:
- Custom query planner for automatic interval join detection
- Physical optimizer rule that converts hash/nested-loop joins to interval joins
- `BioConfig` extension for algorithm selection via `SET bio.*` statements
- `coverage()`, `count_overlaps()`, and `nearest()` SQL table functions

```rust
use datafusion_bio_function_ranges::create_bio_session;

let ctx = create_bio_session();
```

## SQL Table Functions

### `coverage(left_table, right_table [, columns...] [, filter_op])`

Computes base-pair coverage depth. Builds an interval tree from `left_table`, then for each row in `right_table`, computes the total overlap in base pairs with the merged intervals.

```sql
-- Default column names: contig, pos_start, pos_end
SELECT * FROM coverage('reads', 'targets')

-- Custom shared column names
SELECT * FROM coverage('reads', 'targets', 'chrom', 'start', 'end')

-- Separate column names for left and right tables
SELECT * FROM coverage('reads', 'targets', 'chrom', 'start', 'end', 'contig', 'pos_start', 'pos_end')

-- For 0-based half-open coordinates (adjusts boundaries with +1/-1)
SELECT * FROM coverage('reads', 'targets', 'contig', 'pos_start', 'pos_end', 'strict')
```

### `count_overlaps(left_table, right_table [, columns...] [, filter_op])`

Counts overlapping intervals. Same interface as `coverage`, but returns the count of overlapping (non-merged) intervals instead of base-pair overlap.

```sql
SELECT * FROM count_overlaps('reads', 'targets')
```

### `nearest(left_table, right_table [, k] [, overlap] [, columns...] [, filter_op])`

Returns up to `k` nearest left intervals for each right interval.

- `k` default: `1` (must be `>= 1`)
- `overlap` default: `true`
  - `true`: overlapping intervals are returned first, then nearest non-overlaps if needed
  - `false`: overlaps are ignored, only nearest non-overlaps are returned

Output columns are prefixed with `left_` and `right_` to avoid ambiguity.

```sql
-- Default k=1, include overlaps
SELECT * FROM nearest('targets', 'reads')

-- Top-3 nearest per right interval, ignoring overlaps
SELECT * FROM nearest('targets', 'reads', 3, false)
```

### Filter Operations

| Value | Description | Use when |
|-------|-------------|----------|
| `'weak'` (default) | Standard overlap: `start <= end AND end >= start` | 1-based inclusive coordinates |
| `'strict'` | Adjusted boundaries: queries with `start+1, end-1` | 0-based half-open coordinates |

## Interval Join (SQL)

When using a bio-configured session (`create_bio_session()` or `BioSessionExt::new_with_bio()`), SQL joins with range overlap conditions are automatically optimized:

```sql
-- Automatically detected and optimized as interval join
SELECT *
FROM reads
JOIN targets
  ON reads.contig = targets.contig
  AND reads.pos_start <= targets.pos_end
  AND reads.pos_end >= targets.pos_start
```

### Algorithm Selection

```sql
-- Select interval join algorithm
SET bio.interval_join_algorithm = Coitrees;  -- default, best general performance

-- Available algorithms:
-- Coitrees, IntervalTree, ArrayIntervalTree, Lapper, SuperIntervals
-- CoitreesNearest (1 nearest match per right-side row)
```

### Nearest Join

```sql
SET bio.interval_join_algorithm = CoitreesNearest;

SELECT *
FROM targets
JOIN reads
  ON targets.contig = reads.contig
  AND targets.pos_start <= reads.pos_end
  AND targets.pos_end >= reads.pos_start
```

Returns exactly one match per right-side row: the overlapping interval if one exists, otherwise the nearest interval by distance.

## Programmatic API

For direct Rust usage without SQL:

```rust
use std::sync::Arc;
use datafusion_bio_function_ranges::{CountOverlapsProvider, FilterOp};

let provider = CountOverlapsProvider::new(
    Arc::new(ctx.clone()),
    "reads".to_string(),          // left table (built into interval tree)
    "targets".to_string(),        // right table (gets count/coverage column)
    targets_schema,               // Schema of the right table
    vec!["contig".into(), "pos_start".into(), "pos_end".into()],  // left columns
    vec!["contig".into(), "pos_start".into(), "pos_end".into()],  // right columns
    FilterOp::Weak,               // or FilterOp::Strict for 0-based half-open
    true,                         // true = coverage, false = count_overlaps
);
ctx.register_table("result", Arc::new(provider))?;
let df = ctx.sql("SELECT * FROM result").await?;
```

## Migration from sequila-native

This crate replaces the `sequila-core` crate from the [sequila-native](https://github.com/biodatageeks/sequila-native) repository. The functionality is identical; only names and the module structure have changed.

### Type Renames

| sequila-native | datafusion-bio-function-ranges |
|----------------|-------------------------------|
| `sequila_core::session_context::SeQuiLaSessionExt` | `BioSessionExt` |
| `sequila_core::session_context::SequilaConfig` | `BioConfig` |
| `sequila_core::session_context::Algorithm` | `Algorithm` |
| `SessionContext::new_with_sequila(config)` | `SessionContext::new_with_bio(config)` |

### Configuration Namespace

| sequila-native | datafusion-bio-function-ranges |
|----------------|-------------------------------|
| `SET sequila.prefer_interval_join = true` | `SET bio.prefer_interval_join = true` |
| `SET sequila.interval_join_algorithm = Coitrees` | `SET bio.interval_join_algorithm = Coitrees` |
| `SET sequila.interval_join_low_memory = true` | `SET bio.interval_join_low_memory = true` |

### Registration Pattern

**Before (sequila-native):**
```rust
use sequila_core::session_context::{SeQuiLaSessionExt, SequilaConfig};

let mut sequila_config = SequilaConfig::default();
sequila_config.prefer_interval_join = true;

let config = SessionConfig::from(options)
    .with_option_extension(sequila_config);

let ctx = SessionContext::new_with_sequila(config);
```

**After (datafusion-bio-function-ranges):**
```rust
use datafusion_bio_function_ranges::{create_bio_session, register_ranges_functions};

// Simple: creates context with everything configured
let ctx = create_bio_session();

// Or manually:
use datafusion_bio_function_ranges::{BioConfig, BioSessionExt};

let config = SessionConfig::from(options)
    .with_option_extension(BioConfig::default());
let ctx = SessionContext::new_with_bio(config);
register_ranges_functions(&ctx);  // registers coverage() and count_overlaps() UDTFs
```

### New SQL Table Functions

The `coverage`, `count_overlaps`, and `nearest` operations are now available as SQL table functions:

```sql
SELECT * FROM coverage('reads', 'targets')
SELECT * FROM count_overlaps('reads', 'targets')
SELECT * FROM nearest('targets', 'reads')
```

### Dependency Update

**Before:**
```toml
sequila-core = { git = "https://github.com/biodatageeks/sequila-native.git", rev = "..." }
```

**After:**
```toml
datafusion-bio-function-ranges = { git = "https://github.com/biodatageeks/datafusion-bio-functions.git", rev = "..." }
```

## Version Compatibility

| Dependency | Version |
|-----------|---------|
| DataFusion | 50.3.0 |
| Arrow | 56.1.0 |
| Rust edition | 2024 |

These versions must stay in sync with `datafusion-bio-formats` and `polars-bio`.

## License

This crate is licensed under the **Apache License 2.0**, consistent with the rest of the `datafusion-bio-functions` workspace.

The vendored `superintervals` sub-crate (in `superintervals/`) is licensed under the **MIT License** by Kez Cleal. MIT is a permissive license fully compatible with Apache 2.0 — MIT-licensed code can be included in Apache 2.0 projects without restriction.
