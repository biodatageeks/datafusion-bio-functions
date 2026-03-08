# datafusion-bio-function-vcftools

Fused Array Transform optimization for Apache DataFusion — a memory-efficient solution for queries that unnest arrays, transform elements, and re-aggregate.

## Problem

Queries following the pattern **unnest → transform → array_agg** cause memory explosion in DataFusion because the GROUP BY operator materializes all unnested rows before aggregating.

Example problematic query:

```sql
WITH indexed AS (
    SELECT ROW_NUMBER() OVER () as row_idx, values_a, values_b FROM input
),
unnested AS (
    SELECT row_idx, unnest(values_a) as val_a, unnest(values_b) as val_b FROM indexed
),
transformed AS (
    SELECT row_idx, val_a * val_b AS product FROM unnested
)
SELECT row_idx, array_agg(product) FROM transformed GROUP BY row_idx
```

**Memory impact**: 20,000 rows with 2,000-element arrays becomes 40 million intermediate rows, consuming 55+ GB of RAM for a conceptually simple operation.

## Solution

This crate provides a **FusedArrayTransform** operator that replaces the entire `Unnest → Project → Aggregate` chain with a single streaming operator that:

1. Processes rows one at a time (or in configurable batches)
2. Iterates array elements in lockstep (zip semantics)
3. Applies transformations using DataFusion's expression evaluator
4. Builds output arrays directly with bounded memory

**Memory complexity**: O(batch_size × max_array_size) instead of O(total_rows × avg_array_size)

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
datafusion-bio-function-vcftools = { git = "https://github.com/biodatageeks/datafusion-bio-functions" }
```

## Usage

### Basic Setup

```rust
use std::sync::Arc;
use datafusion::prelude::*;
use datafusion::execution::session_state::SessionStateBuilder;
use datafusion_bio_function_vcftools::{
    FusedArrayTransformOptimizerRule,
    VcfQueryPlanner,
};

// Build SessionState with the custom query planner
let state = SessionStateBuilder::new()
    .with_default_features()
    .with_query_planner(Arc::new(VcfQueryPlanner::new()))
    .build();

let ctx = SessionContext::new_with_state(state);

// Register the optimizer rule
ctx.add_optimizer_rule(Arc::new(FusedArrayTransformOptimizerRule::new()));

// Enable the optimization via environment variable
std::env::set_var("BIO_FUSED_ARRAY_TRANSFORM", "1");

// Execute queries - the optimizer will automatically detect and optimize matching patterns
let df = ctx.sql("...").await?;
```

### Configuration

The optimization is controlled via the `BIO_FUSED_ARRAY_TRANSFORM` environment variable:

| Value | Behavior |
|-------|----------|
| `1` or `true` | Optimization enabled |
| Not set or any other value | Optimization disabled |

## Supported Query Patterns

The optimizer detects queries matching this structure:

```
Aggregate (all array_agg functions)
  → Projection* (optional transformation layers)
    → Unnest (one or more array columns)
      → Input
```

### Transforms

```sql
WITH unnested AS (
    SELECT row_idx, unnest(values_a) as val_a, unnest(values_b) as val_b FROM indexed
),
transformed AS (
    SELECT row_idx,
           CASE WHEN val_a > 5 THEN val_a * val_b ELSE val_a + val_b END AS result
    FROM unnested
)
SELECT row_idx, array_agg(result) FROM transformed GROUP BY row_idx
```

### Struct Field Access (VCF-style data)

Works with struct fields containing arrays:

```sql
WITH indexed AS (
    SELECT ROW_NUMBER() OVER () as row_idx, sample_id, genotypes FROM vcf_data
),
unnested AS (
    SELECT row_idx, sample_id,
           unnest(genotypes."GT") as gt,
           unnest(genotypes."DP") as dp
    FROM indexed
)
SELECT row_idx, sample_id, array_agg(gt), array_agg(dp)
FROM unnested
GROUP BY row_idx, sample_id
```

## API Reference

### Public Types

| Type | Description |
|------|-------------|
| `FusedArrayTransform` | Custom logical node representing the fused operation |
| `FusedArrayTransformExec` | Streaming physical operator |
| `FusedArrayTransformOptimizerRule` | Logical optimizer rule that detects the pattern |
| `FusedArrayTransformPlanner` | Extension planner for converting logical to physical |
| `VcfQueryPlanner` | Query planner with FusedArrayTransform support |

### Inspecting the Optimization

Use `EXPLAIN` to verify the optimization is applied:

```sql
EXPLAIN SELECT ...
```

Look for `FusedArrayTransformExec` in the physical plan output.

## Limitations

- All aggregate functions must be `array_agg` (mixed aggregates not supported)
- Requires a row identifier column (typically from `ROW_NUMBER()`)

## File Structure

```
datafusion/bio-function-vcftools/
├── Cargo.toml
├── README.md                              # This file
├── FUSED_ARRAY_TRANSFORM.md               # Detailed design document
├── src/
│   ├── lib.rs                             # Crate root, re-exports
│   ├── logical/
│   │   ├── fused_array_transform.rs       # Logical node
│   │   └── optimizer_rule.rs              # Pattern detection
│   └── physical/
│       ├── fused_array_transform_exec.rs  # Streaming executor
│       ├── extension_planner.rs           # Expr → PhysicalExpr conversion
│       └── query_planner.rs               # Query planner integration
└── tests/
    └── integration_test.rs                # End-to-end tests
```

## License

Apache 2.0
