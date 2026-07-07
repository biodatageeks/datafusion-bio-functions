# translation_sift.lance storage optimization — de-interleaved fixed-divisor blob + decode UDF

Date: 2026-06-22
Status: Design (approved, pending implementation plan)
Branch context: `parallel-cache-redesign`

## Problem

`translation_sift.lance/<chrom>.lance` (position-sliced SIFT/PolyPhen cache) is fast to
query but far larger than it should be. On chr1: **2.9 GB total** (data 2.3 GB,
`_versions` 464 MB, `_indices` 147 MB, 3321 fragments). The payload barely compresses:
~226 raw bytes/row × 13.6M rows = 3.08 GB raw → 2.46 GB on disk = only **~1.25×**.

We want to cut storage substantially **while preserving (or improving) point-lookup
latency**, **without losing precision**, and ideally letting other tools read the cache.

## Root cause (measured)

Current schema: `key: UInt64, sift: Binary, poly: Binary`. Each blob holds one protein
position's predictions as **6-byte interleaved entries**: `[aa u8][pred u8][score f32]`,
~18.7 sift + ~17.7 poly entries/row.

The f32 score is 4 of every 6 bytes (67%). It is **not inherently incompressible** — there
are only **101 distinct SIFT scores** (the values `k/100`) and **≤1001 PolyPhen scores**
(`k/1000`), and 58% of SIFT scores are exactly `0.0`. The isolated score column compresses
~6.5× with plain zlib. **The interleaving is what defeats zstd**: the repeating score byte
patterns are chopped up by the `aa`/`pred` bytes wedged in every 6 bytes, so the
compressor's match-finder can't see the runs.

## Options evaluated (all measured on real chr1 data, Lance 7.0.0, format 2.1)

| Approach | Total | Lossless | Point lookup (batch 64) | Verdict |
|---|---|---|---|---|
| Today (interleaved f32 blob) | 2.9 GB | — | 9.6 ms | baseline |
| Compaction + cleanup only | 2.3 GB | ✅ | unchanged | free −20%, but payload untouched |
| Columnar 6× `List` + item-zstd | ~1.7 GB | ✅ | 28.6 ms (**3× slower**) | column count regresses take |
| Two-table flat-dict | ~1.3 GB | ✅ | 23.5 ms (**2.5× slower**) | range-take + redesign |
| u16 quantization (lossy) | ~1.14 GB | output-only | ~7 ms | precision concern; superseded |
| **De-interleaved fixed-divisor blob** | **~0.95 GB** | ✅ **bit-exact** | **7.1 ms (faster)** | **chosen** |

### Why the alternatives lose

Lance's good encodings (dictionary 4.5×, bit-packing) **only fire on flat scalar
columns**. Every per-position structure that keeps the fast single-`take` either can't be
dict-compressed or crashes:

- `List<Dictionary>` → **deterministic panic** `lance-encoding/src/compression.rs:650:
  unreachable: Per-value compression not yet supported for block type: Dictionary`.
  Reproduced with per-batch dict, unified dict, zstd, forced miniblock, **and
  `compression=none`** — the "per-value compression" is Lance's internal *structural*
  choice for variable-length lists, not our byte-compressor, so no setting avoids it.
- `List<f32>` + zstd on the **item** field → works but caps at 2.5×, and multi-column
  layouts regress the point lookup 2.5–3× (one page fetch per column per `take`).
- `FixedSizeList<f32>` + zstd → zstd **silently ignored** (no compression); raw padded
  floats, worse than today.
- `FixedSizeList<Dictionary>` → **panic** `encoder.rs:572: not yet implemented`.
- Format 2.2 gives no benefit here (slightly worse) — stay on **2.1** (= production).

The only point in the space that is simultaneously lossless, smaller, and fast is to keep
the data in **2 opaque `Binary` columns** (so the take stays fast) and do the
de-interleave + index transform Lance won't, then let Lance zstd the result.

## Chosen design

### 1. Storage format (build)

Schema is **unchanged**: `key: UInt64, sift: Binary, poly: Binary`. The key BTree index
(`sift_key_btree_idx`) and the runtime lookup/query path are untouched.

Blob internal layout changes to **v2 (de-interleaved, fixed-divisor index)**:

```
per row (one protein position, n entries):
  [aa  u8        × n]   amino-acid slot index (0..24), as today
  [pred u8       × n]   prediction enum, as today
  [score_idx     × n]   fixed-divisor index

  sift: score_idx = u8  = round(score × 100)    →  decode score = idx / 100.0
  poly: score_idx = u16 = round(score × 1000)   →  decode score = idx / 1000.0
```

`n` is recovered from the blob length: `n = len / (1 + 1 + idx_width)` where
`idx_width = 1` for sift, `2` for poly.

- **Bit-exactness is verified at build time**, per entry: `round(score × div) / div`
  reconstructs the original f32 exactly (confirmed 0/254.7M sift, 0/245M poly on chr1).
  If any score ever violates the `k/div` grid, the build **fails loudly** rather than
  silently losing precision (preserves the lossless guarantee for future cache sources).
- Index width bounds verified: sift idx ≤ 100 (u8), poly idx ≤ 1000 (u16) — universal for
  SIFT (0.01 resolution) / PolyPhen (0.001 resolution).
- A schema-metadata flag (e.g. `bio.vep.sift_blob_version = "2"`) marks the format so the
  runtime, the UDF, and legacy readers branch correctly. Absence ⇒ v1 (legacy
  interleaved f32).
- **Build finalize**: after writing, run `compact_files` + `cleanup_old_versions`
  (preserving field encoding metadata — verified that pylance compaction preserves it;
  the Rust builder should write few large fragments directly and avoid the
  append-3321-then-compact pattern). This reclaims the ~464 MB of stale `_versions`
  churn (the free −20%) and is applied across all contigs.

Result (chr1): payload **2.2 GB → 788 MB**, total **2.9 GB → ~0.95 GB (−67%)**, lossless,
point lookup **7.1 ms** (vs 9.6 ms today — smaller blobs, less I/O).

### 2. VEP runtime decode

Single source of truth in `cache_common.rs`. Evolve `deserialize_position_entries` /
`deserialize_position_predictions` to read the v2 de-interleaved layout (parse the three
contiguous sub-arrays, map `score_idx → score` via the divisor). Branch on the blob-version
flag; **keep the v1 interleaved-f32 path** for back-compat with un-migrated caches.

`position_predictions_from_batch` (annotate_provider.rs) is unchanged structurally — it
still reads `key`/`sift`/`poly` and delegates to the shared decode core.

### 3. Decode UDF (new, this repo)

A reusable DataFusion `ScalarUDF` so polars-bio / SQL / any DataFusion consumer can read
the cache without reimplementing the byte format:

```
vep_decode_sift_predictions(blob: Binary, predictor: Utf8)  ->  List<Struct<
    amino_acid: UInt8,
    prediction: UInt8,
    score:      Float32,
>>
```

- `predictor` ∈ {`"sift"`, `"polyphen"`} selects the index width and divisor; the UDF is
  **self-contained** (no schema-metadata dependency — divisors are fixed constants).
- Implemented on top of the **same** decode core as the runtime (`cache_common.rs`), so
  there is one canonical decoder. The UDF wraps it to emit Arrow structured output.
- Registered following the repo's `ScalarUDFImpl` conventions.

(Out of scope for v1, candidate follow-up: a `prediction`-code → label string helper.)

### 4. Parity / verification gate

- Unit test: bit-exact roundtrip (encode → decode reproduces original `aa`/`pred`/f32
  `score`) over a fixture and over real chr1 sample rows.
- Unit test: UDF output matches the runtime decode for the same blobs.
- chr1 e2e SIFT/PolyPhen parity vs the VEP reference (all shared CSQ fields), using the
  existing parity infra (`position_sliced_sift_matches_legacy_blob_parity` style).
- Lookup-latency check: v2 take ≤ v1 take (no point-lookup regression).
- Size check: per-contig payload regression guard (v2 < v1).

## Non-goals

- No change to the `translation_sift` schema, key scheme, or BTree index.
- No change to the query/lookup execution path (still a single `take_rows` by key).
- No Lance fork or format-2.2 migration.
- No columnar / two-table / FixedSizeList redesign (all measured-worse on the
  size-vs-lookup trade).

## Risks & mitigations

- **Grid assumption** (scores are exactly `k/100` / `k/1000`): mitigated by the build-time
  bit-exactness assertion — a non-conforming score fails the build instead of losing data.
- **Hand-rolled decode coupling**: mitigated by the shared decode core + the public UDF
  (external tools no longer reimplement the format).
- **Compaction inflating data if encoding metadata is dropped**: builder writes few large
  fragments with field encoding metadata up front; finalize-compaction guarded by the
  size-regression check.
- **Back-compat**: v1 path retained and selected by the version flag; un-migrated caches
  keep working until rebuilt.
