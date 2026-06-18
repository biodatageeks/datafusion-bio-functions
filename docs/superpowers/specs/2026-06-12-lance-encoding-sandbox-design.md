# Lance Encoding Sandbox Design

## Goal

Build a research-only sandbox for comparing Lance 2.1 and 2.2 variation-cache layouts on chr1 with configurable per-field encoding metadata, optional 2.2 struct packing, real post-build Lance inspection, and Rust-only warm/cold benchmarks that always read the full VEPyr `--everything` column set.

## Scope

The sandbox lives under `research/lance_encoding_sandbox/` and does not change the production VEPyr Lance path. It may depend on `datafusion-bio-function-vep` as a local crate dependency to reuse key encoding, VEP projection constants, and external position/bloom index readers.

Only Lance file versions `2.1` and `2.2` are supported. Lance `2.0` is intentionally excluded from the sandbox to keep the matrix focused on the current encoding path and 2.2 struct-packing experiments.

## Dataset Model

The default sandbox layout is chromosome-local. Since each dataset is already split by chromosome, lookup uses:

```text
position = start cast to UInt32
```

The original `start Int64` column remains in the logical `everything` projection for output compatibility. The encoded `position_key UInt64` column is not written in the optimized sandbox layout unless a config explicitly selects `position_key` compatibility mode.

The generated dataset is tiered:

- `tier = 0`: warm rows
- `tier = 1`: cold rows

Indexes are configurable but default to:

- BTree on `position`
- bitmap on `tier`

## Configuration

Each run uses a TOML config. The default config mirrors the current Lance metadata where possible:

- Lance `2.1`
- `miniblock`
- `zstd` level `3`
- minichunk size `16384`
- RLE threshold `0.95`
- dictionary values compressed with `zstd` level `3`

Per-field overrides can set compression such as `fsst` for selected string columns. Struct packing is configured with named groups and requires Lance `2.2`.

Benchmark lookup batch size is explicit:

- `lookup_batch_sizes`: number of keys per `position IN (...)` filter

The sandbox intentionally does not expose row-file sizing, scanner batch sizing, `key.source`, or benchmark projection:

- warm and cold rows use Lance's default `WriteParams.max_rows_per_file`, currently `1_048_576` rows in Lance 7.0.0
- warm, cold, and key-extraction scans use Lance's default scanner batch size
- `position` mode always derives `position UInt32` from `start`
- benchmarks always read the full VEPyr `--everything` logical projection

Unknown TOML fields are rejected so removed knobs fail fast instead of being silently ignored.

## Build Flow

1. Read warm/cold chr parquet files from `<cache_root>/variation`.
2. Project the configured logical `everything` fields.
3. Add `tier UInt8`.
4. Add `position UInt32` from `start`.
5. Optionally pack configured fields into struct columns for Lance 2.2.
6. Apply default and per-field Lance metadata to the actual Arrow schema.
7. Write warm and cold rows with Lance's default `WriteParams.max_rows_per_file`.
8. Create configured BTree/bitmap indexes.

## Inspection Flow

After every build, the CLI reopens the generated Lance dataset and reports observed state, not requested config:

- schema field path
- Arrow type
- nullability
- actual metadata on each field
- physical encoding debug strings observed from Lance file metadata
- compressed bytes per physical/logical column where Lance file metadata can map them
- page count
- data/index/metadata/total filesystem sizes

The inspector also emits a requested-vs-observed table so FSST, RLE, dictionary, packed struct, and general compression decisions are visible after writing.

## Benchmark Flow

The benchmark always reads all columns required by the sandbox `everything` projection. For packed 2.2 layouts, logical child fields are mapped to their struct containers, and the physical projection is recorded in the report.

Warm benchmark:

```sql
tier = 0
```

Cold benchmark:

```sql
tier = 1 AND position IN (...)
```

Cold benchmarks sweep every configured `lookup_batch_size`.

Each run records:

- selected rows
- selected positions
- result batches
- elapsed seconds
- Lance bytes read
- Lance requests
- IOPS
- ranges scanned
- fragments scanned
- rows scanned
- indexes/parts loaded
- index comparisons

## Outputs

Each run writes:

```text
<output_root>/<dataset.name>/
  <chrom>.lance/
  reports/
    config.effective.toml
    build.json
    inspect.json
    benchmark.json
    summary.md
```

A Jupyter notebook under `notebooks/` runs the Rust CLI and renders comparison tables from JSON. Python is only orchestration and visualization; benchmark timing is produced by Rust.

## Validation

The CLI fails before benchmarking if:

- Lance version is not `2.1` or `2.2`
- struct packing is enabled for a non-2.2 config
- `position UInt32` cannot be derived from `start`
- any required logical `everything` field cannot be resolved to a top-level field or packed struct
- the configured key/index column is absent from the generated dataset
