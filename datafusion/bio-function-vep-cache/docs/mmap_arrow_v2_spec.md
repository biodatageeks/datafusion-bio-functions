# VEP Cache Storage v2 Spec (mmap-backed Arrow Blocks)

## Status
- Draft
- Target crate: `datafusion/bio-function-vep-cache`
- Scope: on-disk storage and read/write APIs for v2 cache format
- Current implementation uses true OS mmap on unix in `MmapArrowBlockStore`.

## Motivation
Current v1 layout stores one Arrow IPC value per `(chrom, window, column)` in Fjall.
For wide projections, read path cost is dominated by:
- `windows * selected_columns` KV lookups
- `windows * selected_columns` IPC decode calls

v2 shifts window payloads to mmap-backed Arrow blocks and keeps Fjall as metadata/index.
Goal is to make window fetch cost scale mainly with windows, not columns.

## Goals
- Make mmap-backed Arrow blocks the primary data source for cache payloads.
- Keep Fjall as authoritative metadata + lookup index store.
- Preserve fast position index lookup for matching.
- Support incremental migration: v1 and v2 caches can coexist.
- Keep API surface backward compatible for existing v1 users.

## Non-goals (first phase)
- Full v2 integration in `KvLookupExec` and `CacheLoader` hot paths.
- Query-time dynamic re-partitioning.
- Cross-process lock-free writers.
- Perfect zero-copy for every Arrow logical type in every path.

## Format Versioning
- `format_version=1`: existing per-column Fjall entries.
- `format_version=2`: mmap-backed Arrow window blocks + Fjall metadata/index.

Opening behavior:
- v1 and v2 are accepted.
- any other version is rejected.

## v2 Storage Layout
Cache root directory:

```text
<cache_root>/
  fjall/                    # existing keyspace files (managed by Fjall)
  blocks_v2/                # mmap payload files
    block_000000.arrowblk
    block_000001.arrowblk
    ...
```

Notes:
- `blocks_v2` is configurable via metadata key (`v2_block_dir`), default `blocks_v2`.
- phase 1 may write only `block_000000.arrowblk`; rotation is optional.

## Fjall Metadata Keys
- `schema`: Arrow IPC schema bytes (existing)
- `window_size`: u64 big endian (existing)
- `format_version`: u8 (existing)
- `v2_block_dir`: utf8 path relative to cache root, required for v2
- `v2_block_codec`: `none | lz4 | zstd` for payload IPC blobs
- `v2_window_ref:<window_key>`: encoded block location for a window

Where:
- `<window_key>` = existing 10-byte `(chrom, window_id)` encoding.

## Window Payload Block Format (`.arrowblk`)

File header:
- magic bytes: `VEPBLK2\0`
- repeated frames:
  - `frame_len`: u64 little endian
  - `frame_payload`: columnar payload bytes (`V2COLM01`)

Columnar payload (`V2COLM01`):
- header:
  - `magic`: `V2COLM01`
  - `row_count`: u32 LE
  - `num_cols`: u16 LE
  - repeated entries (`num_cols`):
    - `col_idx`: u16 LE (cache schema position)
    - `offset`: u32 LE (payload-relative)
    - `length`: u32 LE
- body:
  - concatenated one-column Arrow IPC blobs in schema order
  - each blob may use IPC compression from `v2_block_codec`

## Window Reference Encoding
Each `v2_window_ref` value is fixed-width:

```text
file_id   : u32 LE
offset    : u64 LE   # payload start (not frame_len start)
length    : u32 LE   # payload byte length
row_count : u32 LE
```

Total: 20 bytes.

## Read Path (v2)
1. Lookup position index in Fjall (`PositionIndex` entry).
2. Resolve window payload reference from `v2_window_ref:<window_key>`.
3. Open/mmap payload file (`file_id`).
4. Slice `[offset..offset+length]`.
5. Parse columnar payload header and decode only requested column blobs.

Performance expectation:
- One metadata lookup + one mmap slice per window.
- Decode cost scales with selected columns (not full window width).

## Write Path (v2)
For each merged `(chrom, window)` RecordBatch:
1. Build and write `PositionIndex` into Fjall.
2. Serialize each cache column as one-column IPC blob and append as one columnar frame.
3. Store `v2_window_ref` entry in Fjall metadata.

Atomicity model in phase 1:
- Position index and window ref are committed through Fjall.
- Payload append is file-system append; recovery strategy is based on ref validity.

## Recovery and Consistency
- A `v2_window_ref` pointing outside file bounds is treated as corruption.
- Missing payload file for a valid ref is treated as corruption.
- Loader should optionally support verification scan:
  - check each ref decodes and row_count matches batch rows.

## Migration Strategy
- No in-place conversion required.
- Create v2 cache from source parquet using loader mode switch.
- Existing v1 caches remain readable without changes.

## API Additions (phase 1)
- `FORMAT_V2` constant.
- v2 block store abstraction for append/read via mmap.
- `VepKvStore::create_v2(...)`.
- `VepKvStore::put_window_block_v2(...)`.
- `VepKvStore::get_window_block_v2(...)`.
- `VepKvStore::get_window_columns_v2(...)` for projection-aware fetch.

## Performance Impact Estimates
Using measured baseline (`all columns`, 1 thread):
- baseline: ~4k rows/s
- fetch stage: ~85.7% of total

Current measured direction:
- all-columns: near v1 parity (~3.8k-3.9k rows/s vs ~4.1k rows/s baseline)
- three-column projection: large improvement vs row-batch v2 path (order-of-magnitude)

## Implementation Plan
1. Add v2 storage primitives (block files, refs, mmap reads).
2. Add `format_version=2` support in KV store open/create.
3. Add loader mode to write v2 payloads + refs.
4. Add lookup execution mode to fetch cache columns from v2 window batches.
5. Add benchmark gates for v1/v2 parity and v2 throughput targets.
