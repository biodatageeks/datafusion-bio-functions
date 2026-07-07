# Arrow Mmap Cold Lookup Benchmark Design

## Context

The current Lance cold lookup benchmark for chr1 shows large read amplification for a sparse VEPyr workload: about 9,497 cold positions and 11,463 returned rows out of tens of millions of cache rows. The first Lance scan can also load a large row-level BTree index. We need a Rust benchmark that compares Lance against a purpose-built cold lookup layout:

```text
position_key -> (batch_id, start_row, len)
```

with Arrow IPC batch files opened through mmap-style access.

This is a benchmark-only design. It is not a production VEP backend until the benchmark shows that the layout is worth integrating.

## Goals

- Compare Lance cold lookup with a compact external range index plus mmap Arrow files.
- Use the same chr1 E2E-derived cold workload as `bench_lance_cold_batching`.
- Read and touch all projected columns required by the `everything` path.
- Measure both performance and disk footprint.
- Build three Arrow variants: uncompressed, LZ4-compressed IPC, and ZSTD-compressed IPC.
- Keep the layout cold-first while leaving room to add warm-tier full-scan support later.

## Non-Goals

- Do not add a production VEP annotation backend in this phase.
- Do not replace the existing position or bloom sidecars.
- Do not implement a custom raw Arrow-buffer file format in this phase.
- Do not claim mmap estimated touched bytes are directly equivalent to Lance `bytes_read`.

## Artifacts

The benchmark writes one store per compression mode:

```text
chr1.arrow_mmap_uncompressed/
  manifest.json
  index.mdb/
  batches/batch_000000.arrow
  batches/batch_000001.arrow

chr1.arrow_mmap_lz4/
  manifest.json
  index.mdb/
  batches/batch_000000.arrow

chr1.arrow_mmap_zstd/
  manifest.json
  index.mdb/
  batches/batch_000000.arrow
```

The benchmark should add `heed` as a dev/example dependency and reuse the existing Arrow 58 IPC dependency already present in the crate.

## Index Format

The LMDB key is the `position_key` encoded as an 8-byte integer in a stable byte order.

The value is a fixed-width record:

```rust
struct PositionRange {
    batch_id: u32,
    start_row: u32,
    len: u32,
}
```

The writer must preserve this invariant:

```text
All rows for one position_key are contiguous and live in exactly one Arrow batch file.
```

If a position group would cross a batch boundary, the writer starts a new Arrow batch file before writing that group. If one position group is larger than the configured batch row target, the group is written as a single oversized batch file and the manifest records the maximum observed group length.

## Arrow Data Format

Each batch file stores a standard Arrow IPC file with the projected cold cache schema. The first implementation writes the same projected `everything` columns used by the Lance benchmark. The store variants differ only by Arrow IPC compression:

- `uncompressed`: no IPC compression, expected to be the closest to direct mmap buffer access.
- `lz4`: Arrow IPC LZ4 compression, expected to trade CPU for lower size.
- `zstd`: Arrow IPC ZSTD compression, expected to minimize size but require more CPU.

Compressed variants are not treated as pure direct mmap reads because selected data requires IPC decompression.

## Manifest

`manifest.json` records:

- format name and version
- chromosome
- tier, initially `cold`
- compression mode
- Arrow schema metadata
- projected column names
- target batch rows
- batch count
- total rows
- unique indexed positions
- maximum position group length
- data file size
- LMDB index size
- total store size
- build parameters

The reader validates the manifest before lookup. Compression, schema, projected columns, and batch file count must match the files on disk.

## Build Flow

1. Read chr1 cold variation parquet rows from the merged cache.
2. Project the same columns used for the Lance `everything` workload.
3. Require rows to be sorted by `position_key`.
4. Buffer one complete `position_key` group at a time.
5. Write groups into Arrow IPC batch files without splitting a group.
6. Insert one LMDB range entry per unique `position_key`.
7. Write `manifest.json` after all files and index entries are complete.

The builder should emit build timing and final size metrics for each compression mode.

## Benchmark Flow

1. Parse the same chr1 VCF workload used by `bench_lance_cold_batching`.
2. Load warm positions for gating.
3. Apply the existing external position index.
4. Apply the existing variant bloom index.
5. Use the resulting bloom-admitted cold positions as the common lookup set.
6. Run the existing Lance cold scan benchmark.
7. Run each Arrow mmap store against the same lookup positions.
8. For every returned row, touch all projected columns.
9. Emit one markdown report comparing Lance, uncompressed Arrow, LZ4 Arrow, and ZSTD Arrow.

## Column Touch Semantics

The benchmark must do enough work to prevent the compiler from optimizing the lookup away:

- numeric columns: read scalar values and mix them into a checksum
- string/binary columns: read offsets and value bytes
- list columns: read offsets and child values or child byte spans
- null arrays: read validity state

Each Arrow variant must produce the same selected row count, selected position count, and deterministic checksum.

## Metrics

Common metrics:

```text
format
compression
total size
data size
index size
manifest/schema size
build seconds
open seconds
lookup seconds
positions requested
positions found
selected rows
projected columns
estimated bytes touched
batch files opened
index gets
missing positions
checksum
```

Lance-specific metrics retained from the current benchmark:

```text
bytes_read
requests
ranges_scanned
fragments_scanned
parts_loaded
index_comparisons
```

Arrow mmap reports `estimated_bytes_touched`, not `bytes_read`, because benchmark code does not observe OS page-cache reads directly.

## Error Handling

Build errors are fatal for:

- missing `position_key`
- missing projected column
- rows out of `position_key` order
- Arrow IPC write failure
- LMDB write failure
- duplicate LMDB key with a different range
- manifest write failure

Lookup errors are fatal for:

- missing manifest
- unsupported compression
- missing batch file
- range outside the target batch
- schema or projected-column mismatch

An LMDB key miss is not fatal. It increments the benchmark miss counter.

## Tests

Small Rust tests should cover:

- position groups are not split across batch files
- LMDB key lookup returns the expected `(batch_id, start_row, len)`
- uncompressed, LZ4, and ZSTD stores return identical rows and checksums for the same positions
- missing position returns a miss
- manifest validation rejects wrong compression, schema, or batch count

Benchmark verification should compare Arrow selected rows and selected positions against Lance for the same gated workload and report any divergence.

## Acceptance Criteria

- A Rust benchmark can build the three Arrow mmap stores for chr1 cold rows.
- The benchmark runs against the same chr1 gated cold workload as the Lance side benchmark.
- The report includes Lance and all three Arrow variants.
- Arrow variants agree with each other on row count, position count, and checksum.
- Divergence from Lance row or position counts is reported explicitly.
- Disk size is broken down into data, index, manifest/schema, and total size.
- The benchmark clearly labels compressed Arrow variants as decompression-backed rather than pure direct mmap.
