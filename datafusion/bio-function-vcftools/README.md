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
    enable_fused_array_transform,
};

// Build SessionState with the custom query planner
let state = SessionStateBuilder::new()
    .with_default_features()
    .with_query_planner(Arc::new(VcfQueryPlanner::new()))
    .build();

let ctx = SessionContext::new_with_state(state);

// Register the optimizer rule
ctx.add_optimizer_rule(Arc::new(FusedArrayTransformOptimizerRule::new()));

// Enable the optimization
enable_fused_array_transform();

// Execute queries - the optimizer will automatically detect and optimize matching patterns
let df = ctx.sql("...").await?;
```

### Configuration

The optimization is disabled by default. Use the API functions to control it:

```rust
use datafusion_bio_function_vcftools::{
    enable_fused_array_transform,
    disable_fused_array_transform,
    set_fused_array_transform_enabled,
};

// Enable globally
enable_fused_array_transform();

// Disable globally
disable_fused_array_transform();

// Conditional enabling
set_fused_array_transform_enabled(config.use_fused_transform);
```

For backward compatibility, the environment variable `BIO_FUSED_ARRAY_TRANSFORM` is still supported:

| Value | Behavior |
|-------|----------|
| `1` or `true` | Optimization enabled |
| Not set or any other value | Optimization disabled |

**Note:** The API takes precedence over the environment variable.

## Supported Query Patterns

The optimizer detects queries matching this structure:

```
Aggregate (all array_agg functions)
  → SubqueryAlias* (optional, transparent wrappers)
    → Projection* (optional transformation layers)
      → SubqueryAlias* (optional)
        → Unnest (one or more array columns)
          → Input
```

### Requirements

The optimization is applied **only when all** of the following conditions are met:

1. **Aggregate node with only `array_agg`/`list_agg` functions** — all aggregate expressions must be `array_agg` or `list_agg`. Mixed aggregates (e.g., `array_agg` + `sum`) are not supported.

2. **No unsupported modifiers** — the `array_agg` functions must not use:
   - `DISTINCT` (`array_agg(DISTINCT col)`)
   - `FILTER` (`array_agg(col) FILTER (WHERE ...)`) 
   - `ORDER BY` inside the aggregate (`array_agg(col ORDER BY other_col)`)

3. **Unnest node reachable** — below the Aggregate, through zero or more Projection and SubqueryAlias nodes, there must be an `Unnest` node.

4. **At least one array column** — the Unnest must operate on at least one array column.

5. **Row-identity column in GROUP BY** — the GROUP BY clause must include a column that uniquely identifies each original input row created using `ROW_NUMBER() OVER ()` with an **empty PARTITION BY** clause.

   **Why this matters:** The `FusedArrayTransformExec` operator emits at most one output row per input row and does not perform cross-row aggregation. If the GROUP BY doesn't uniquely identify rows, the baseline `array_agg` would merge values across rows, but the optimized plan would return separate arrays per row — producing incorrect results.

   **Valid example (optimization applied):**
   ```sql
   WITH indexed AS (
       SELECT ROW_NUMBER() OVER () as row_idx, metadata, values FROM input
   ),
   unnested AS (
       SELECT row_idx, metadata, unnest(values) as val FROM indexed
   )
   SELECT row_idx, metadata, array_agg(val)
   FROM unnested
   GROUP BY row_idx, metadata  -- row_idx uniquely identifies rows
   ```

   **Invalid example (optimization NOT applied, falls back to standard DataFusion):**
   ```sql
   WITH unnested AS (
       SELECT metadata, unnest(values) as val FROM input
   )
   SELECT metadata, array_agg(val)
   FROM unnested
   GROUP BY metadata  -- metadata may not be unique across rows!
   ```

   The optimizer detects `ROW_NUMBER() OVER ()` (without PARTITION BY) in the plan and tracks its column name through projections and aliases to verify it appears in the GROUP BY.

### Order Preservation

The fused transformation **preserves element order** within arrays. Unlike the SQL-based `unnest → transform → array_agg` pattern (where `array_agg` order depends on `GROUP BY` execution and is not guaranteed without `ORDER BY` inside the aggregate), this optimization processes elements sequentially and outputs them in the same order as the input arrays.

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

- All aggregate functions must be `array_agg` or `list_agg` (mixed aggregates not supported)
- `DISTINCT`, `FILTER`, and `ORDER BY` modifiers on `array_agg` are not supported
- **GROUP BY must include a row-identity column** (from `ROW_NUMBER() OVER ()` without PARTITION BY) to ensure each input row maps to exactly one output row
- The query must follow the `Aggregate → Projection* → Unnest` structure

## File Structure

```
datafusion/bio-function-vcftools/
├── Cargo.toml
├── README.md                              # This file
├── src/
│   ├── lib.rs                             # Crate root, re-exports
│   ├── common.rs                          # Common functions
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
