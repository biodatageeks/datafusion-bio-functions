# Position-Sliced SIFT/PolyPhen Predictions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-transcript SIFT/PolyPhen prediction blob in `translation_sift.lance` with a position-sliced layout so a missense lookup reads only the needed amino-acid position (~tens–hundreds of bytes) instead of the transcript's entire ~127 KB substitution matrix — driving the whole-chromosome `sift_load` cost (~11 s single-threaded on chr1) to well under 1 s, with byte-identical annotation output.

**Architecture:** Each transcript gets a dense `transcript_uid: u32` written into the **`transcript`** table (a transcript-wide identity, reusable by other per-transcript datasets later). `translation_sift` is rebuilt as **one row per `(transcript, protein_position)`**, keyed by a single bit-packed `key: u64 = (transcript_uid << 32) | position` with a BTree scalar index. The uid is assigned by a shared deterministic helper (dense rank over the chrom's sorted-unique transcript ids from the **transcript source**), so the transcript build and the sift build derive identical uids without depending on each other's build order or artifacts. At query time the engine reads `transcript_uid` from the in-memory transcript record (no sidecar dict), builds `key` for each missense `(transcript, position)`, and fetches just those rows via batched `take_rows`. The runtime detects the dataset format (legacy per-transcript blob vs new position-sliced) from the schema, so existing caches keep working.

**Tech Stack:** Rust 2024, DataFusion 53, Arrow 58, Lance 7.0.0, existing `lance_cache` build/runtime, `kv_cache::sift_store` serialization, `transcript_consequence::{CachedPredictions, CompactPrediction}`.

---

## Background / Evidence (measured chr1, single-threaded)

- `translation_sift.lance/chr1.lance`: 39,340 transcripts (12,200 non-coding/zero predictions), `predictions` blob **mean 127 KB / median 83 KB / max 2.1 MB**, **1.33 GB compressed**, FullZip + zstd. Per transcript: dense `[protein_length × 19 AA]` matrices for SIFT and PolyPhen (mean 12,703 entries; protein length mean 341, max 7,969 on chr1; genome-wide max ≈ 34,350 = Titin/TTN on chr2).
- Whole-chromosome `--everything` annotation touches most coding transcripts → reads most of the 1.33 GB even though each missense variant needs **one `(transcript, position, alt_aa)`** value. After the demand-driven SIFT change ([[sift-load-bottleneck]]) the cost merely relocated into the engine (~11 s); the dominant cost is **blob volume + whole-matrix granularity**.
- Lance 7 supports single-column scalar BTree only, so the composite key is bit-packed into one `u64` (`Dataset::take_rows` resolves it directly). `u32 + u32` makes position-overflow impossible by construction (vs sizing position bits from a single chromosome).

Brainstorm decisions (settled):
- **Key:** `u64 = (transcript_uid as u64) << 32 | (position as u64)`. uid in high 32 bits (orders by uid then position), position in low 32 bits. No bit-budget guard needed beyond "both are u32".
- **Mapping source:** `transcript_uid` rides on the **`transcript`** table, loaded into the in-memory transcript record — no separate sidecar file or query-time map. Both builds assign it identically via a shared rank-over-transcript-source helper (no build-order/artifact coupling).
- **Payload:** keep exact `f32` scores for output parity (score quantization is a later, orthogonal optimization). Position is implicit from `key`; AA stays explicit per entry.
- **Granularity:** one row per `(transcript, position)` (all ~19 alt-AAs packed) — not per-AA (too many rows) and not per-transcript (today's problem).

---

## Row Format (new `translation_sift` position-sliced schema)

Per `(transcript, position)` row:

| column | type | notes |
|---|---|---|
| `key` | `UInt64` | `(transcript_uid << 32) | position`; BTree-indexed; rows sorted by `key` |
| `sift` | `Binary` | packed SIFT entries for this position |
| `poly` | `Binary` | packed PolyPhen entries for this position |

`sift`/`poly` payload = sequence of 6-byte entries `(amino_acid: u8, prediction: u8, score: f32 LE)`; entry count = `len / 6`. Position is taken from `key`, so it is **not** stored per entry (drops 4 bytes/entry vs the legacy format). Empty column = no predictions of that kind at this position.

The `transcript` table gains one column:

| column | type | notes |
|---|---|---|
| `transcript_uid` | `UInt32` | dense rank of the transcript stable id within the chrom's sorted-unique transcript set (from the transcript source) |

> Note: the `transcript` table's identity column is `stable_id`; `translation_sift`/`translation_core` use `transcript_id`. The plan assumes these are the same Ensembl transcript stable id namespace — **verify during Task 2** (assert the sift transcript ids are a subset of the transcript-table ids) before relying on it.

---

## File Structure

- Modify: `datafusion/bio-function-vep/src/kv_cache/sift_store.rs`
  - Add `serialize_position_predictions` / `deserialize_position_predictions` (position-implicit, 6-byte entries) and a small `assign_transcript_uids(sorted_unique_ids) -> HashMap<String,u32>` shared helper.
- Modify: `datafusion/bio-function-vep/src/lance_cache/build.rs`
  - Shared `load_transcript_uid_map(options, chrom) -> HashMap<String,u32>`: query the **transcript source** for the chrom's sorted-unique transcript ids and `assign_transcript_uids`. Used by both the transcript build and the sift build → identical uids, no build-order coupling.
  - `build_lance_context_entity` (Transcript path): emit `transcript_uid` column from the shared map.
  - `build_lance_translation_sift`: load the uid map, explode predictions into position rows keyed by `u64`; new schema; BTree on `key`; validation.
  - New `compact_translation_sift_position_schema` + `transform_translation_sift_position_batch`.
- Modify: `datafusion/bio-function-vep/src/lance_cache/write.rs`
  - Add a `LanceIndexKind::SiftKey` (BTree on `key`) or a generic BTree-on-column helper.
- Modify: `datafusion/bio-function-vep/src/lance_cache/context_runtime.rs`
  - Add a `KeyU64LanceLookup` mirroring `TranscriptIdLanceLookup`: load the `key` BTree page data into an in-memory `key → row_id` index, resolve queried keys → batched `take_rows`. No SQL `IN`/filter-scan (direct `take_rows` is ~3–5× faster per the variation bench).
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs`
  - `load_lance_sift_prediction_store_for_chrom`: detect `key`/`sift`/`poly` schema → build a `PositionSlicedLanceSiftStore`; keep the legacy `predictions` and in-memory branches.
  - Carry `transcript_uid` into the loaded **transcript** record (transcript context load), available at the `lookup_sift_polyphen` call site.
  - `lookup_sift_polyphen`: build `key` from `(transcript_uid, position)`, fetch position payload (memoized), look up `alt_aa`.
  - Batched demand fetch: collect a batch's missense `key`s → one `get_many`.
- Tests: `build.rs`, `sift_store.rs`, `annotate_provider.rs` unit tests; chr1 parity + bench validation.

---

## Task 1: Serialization + uid assignment helpers

**Step 1: Position-implicit payload codec**
- [ ] In `sift_store.rs`, add `serialize_position_predictions(sift: &[CompactPrediction], poly: &[CompactPrediction]) -> (Vec<u8>, Vec<u8>)` emitting 6-byte entries `(aa u8, pred u8, score f32 LE)` (no position). All entries passed in must share one position.
- [ ] Add `deserialize_position_predictions(position: i32, sift: &[u8], poly: &[u8]) -> CachedPredictions` that reconstructs `CompactPrediction { position, amino_acid, prediction, score }`, sorted by `(position, aa)` (it already is, within a position) for `lookup_in`'s binary search.
- [ ] Unit test: round-trip a position's entries; assert lengths `% 6 == 0` and decoded values match.

**Step 2: Deterministic uid assignment**
- [ ] Add `assign_transcript_uids(sorted_unique_ids: &[String]) -> HashMap<String, u32>` = dense rank 0..N over the **sorted-unique** id list. Pure function of the id list, so the transcript build and the sift build agree as long as both feed it the same source id set.
- [ ] Unit test: stable, dense, unique; same input → same map.

---

## Task 2: Build `transcript_uid` into the `transcript` table

- [ ] Add `load_transcript_uid_map(options, chrom)`: query the transcript source (`tx`) for the chrom's sorted-unique stable ids → `assign_transcript_uids`. This is the single shared rule both builds use.
- [ ] **Verify id namespace:** assert/confirm the transcript table's identity column (`stable_id`) is the same namespace as `translation_sift.transcript_id` (Ensembl transcript stable id). If they differ, resolve the join column before proceeding.
- [ ] In the Transcript path of `build_lance_context_entity`, extend the projected schema with `transcript_uid: UInt32` (append like the variation `tier` column was) and populate it from the map (keyed by `stable_id`).
- [ ] Keep existing indexes; `transcript_uid` is a plain column (read with the row, no index needed).
- [ ] Test: a synthetic transcript build emits dense, unique `transcript_uid` matching the sorted stable-id order.

---

## Task 3: Build position-sliced `translation_sift`

**Step 1: New schema + explode transform**
- [ ] Add `compact_translation_sift_position_schema(source_type)` → fields `key: UInt64` (BTree-friendly miniblock + zstd), `sift: Binary`, `poly: Binary` (miniblock + zstd; small payloads, no FullZip).
- [ ] In `build_lance_translation_sift`, call `load_transcript_uid_map(options, chrom)` (same shared helper as Task 2 → identical uids; no dependency on the transcript dataset being built).
- [ ] Add `transform_translation_sift_position_batch(batch, source_type, &uid_map)`:
  - For each transcript: read `transcript_id`, look up `transcript_uid` (error if absent — see validation), decode the source SIFT/PolyPhen entries.
  - Group entries by `position` (entries are already sorted by `(position, aa)`); for each position emit one output row: `key = ((uid as u64) << 32) | (position as u64)`, with `serialize_position_predictions` for that position's SIFT/PolyPhen slices.
  - Skip transcripts/positions with no predictions (non-coding → no rows).
  - Assert `position >= 0` and fits (it's ≤ ~34k; `u32` is automatic). Assert `uid` came from the map.
- [ ] Rows are emitted per source transcript in `(transcript_id)` order; since `uid` follows sorted `transcript_id` and positions are ascending, `key` is naturally ascending → emit sorted, or sort the stream by `key` before write (use the existing coalesce/sort path) so the BTree + locality hold.

**Step 2: Index + wiring**
- [ ] Build a BTree index on `key` (add `LanceIndexKind::SiftKey` → column `"key"`, name e.g. `sift_key_btree_idx`, or a generic helper).
- [ ] Route `build_lance_translation_sift` through the new transform + schema + index. Build order already runs core before sift; pass the `uid_map` (recomputed identically via `assign_transcript_uids`, or read from core — prefer recompute from the same dedup to avoid a core read).

**Step 3: Build-time validation**
- [ ] Assert every sift `transcript_id` is present in the uid map (i.e. ⊆ the transcript source id set); hard-error on any missing id.
- [ ] Assert `key`s are unique and ascending (no duplicate `(uid, position)`).
- [ ] Integration test (synthetic): 2 transcripts × few positions → assert row count = total positions, `key` decode round-trips `(uid, position)`, BTree present, and a `take_rows` by `key` returns the right position payload.

---

## Task 4: Runtime — position-keyed SIFT store + engine wiring

**Step 1: `KeyU64LanceLookup` — load BTree index, batched `take_rows` (same mechanism as `TranscriptIdLanceLookup`)**
- [ ] In `context_runtime.rs`, add a lookup that at `open` **loads the `key` BTree page data into an in-memory `key → row_id` index** — exactly how `TranscriptIdLanceLookup` builds a `StringRowIdIndex` via `load_transcript_id_index_from_lance_btree`, but for `u64` keys. Add a `row_index`-style `U64RowIdIndex` (sorted `(key, row_id)` pairs from the index pages; the real `row_id`s come straight from the BTree, so this does not depend on storage order).
- [ ] `resolve(keys)`: look each queried `key` up in the in-memory index → `row_id`; collect and **`take_rows(row_ids, [key, sift, poly])` in batches**. No `scan().filter("key IN …")` — direct `take_rows` is ~3–5× faster (variation bench).
- [ ] Resolution tolerates keys not present (a position with no predictions / not built) → skip.
- [ ] Memory note: the in-memory index is bounded (~13.4 M `(key,row_id)` entries for chr1, comparable to / smaller than variation's existing start→row_id preload) and is loaded once per active chrom (prefetchable).

**Step 2: `PositionSlicedLanceSiftStore`**
- [ ] In `annotate_provider.rs`, add a store holding the `KeyU64LanceLookup` + a memoizing `HashMap<u64, Option<CachedPredictions>>` (per `(uid,position)` key).
- [ ] `get_position_predictions(keys: &[u64]) -> HashMap<u64, CachedPredictions>`: dedup uncached keys, one batched `take_rows`, `deserialize_position_predictions` (position from `key & 0xFFFF_FFFF`), memoize.
- [ ] In `load_lance_sift_prediction_store_for_chrom`, add a schema branch: `key` + `sift` + `poly` present → `PositionSlicedLanceSiftStore`; keep the legacy `predictions` branch and the in-memory scan branch.

**Step 3: Carry `transcript_uid` to the engine (no runtime map)**
- [ ] Add `transcript_uid: u32` to the in-memory `Transcript` struct and populate it when loading the `transcript` context (extend the loaded projection).
- [ ] At the SIFT call sites (annotate_provider.rs ~6178/6636), `tx_opt: Option<&Transcript>` is already in scope (used for GENE_PHENO immediately above) — pass `tx_opt.map(|tx| tx.transcript_uid)` into `lookup_sift_polyphen`. No `transcript_id → uid` map is needed; the uid rides on the record.
- [ ] Confirm `Transcript.transcript_id` (loaded from the transcript table's `stable_id`) equals the translation `transcript_id` namespace — same check as Task 2 (the build assert is the hard gate).

**Step 4: `lookup_sift_polyphen` position path**
- [ ] When the store is position-sliced and `transcript_uid` is `Some`: parse `(protein_position, alt_aa)`; build `key = (uid << 32) | position`; fetch the position's `CachedPredictions` (memoized by `key`) via the store; `lookup_sift`/`lookup_polyphen` for `alt_aa`. `transcript_id` is **not** used here.
- [ ] If `transcript_uid` is `None` (no transcript record in scope), skip SIFT rather than guess.
- [ ] Preserve the legacy `transcript_id`-keyed path when the store is the old per-transcript kind (format detection drives which path).

**Step 5: Batch the demand fetch (folds in the earlier "batch demand fetch" idea)**
- [ ] Per engine batch, collect the missense `key`s, call `get_position_predictions` once, then format from cache — avoiding per-variant `block_in_place`/single-row takes.

---

## Task 5: Validation

- [ ] Unit + integration tests above pass; `cargo clippy --lib` clean; `cargo fmt`.
- [ ] Rebuild chr1 `translation_core` + `translation_sift` (fresh dir first, then in place) via the build path; assert row counts (positions) and indexes.
- [ ] Re-run the single-threaded profiled annotation (`LANCE_*_THREADS=1 RAYON_NUM_THREADS=1`, `VEP_PROFILE=1`, `--backend lance --forks 0`); compare `sift_load`/engine vs today (target: ~11 s → < 1 s; bytes read for SIFT from ~1.33 GB → a few MB).
- [ ] **Parity gate:** diff annotation output (SIFT/PolyPhen `prediction(score)` fields) against the legacy per-transcript build for chr1 (and chr4) — must be byte-identical. This is the acceptance criterion.

---

## Risks & Notes

- **uid consistency** is the main correctness risk: the transcript build and the sift build must assign identical uids. Mitigation: both feed the **same shared `load_transcript_uid_map`** (rank over the transcript source's sorted-unique ids) — a pure function of the same source, so no build-order/artifact coupling — plus the build-time subset + key-uniqueness asserts.
- **Id namespace:** transcript table keys on `stable_id`, sift on `transcript_id`. The plan assumes one Ensembl stable-id namespace; Task 2 verifies it. If they diverge (e.g. version suffixes), normalize the join key.
- **Row count** grows from ~39 k (transcripts) to ~13.4 M (positions) per chrom (~170 M genome-wide), but each row is tiny and only hit position **payloads** are read. The `key → row_id` BTree index is **loaded into memory** per active chrom (~13.4 M entries; comparable to / smaller than variation's existing start→row_id preload) — a bounded, prefetchable startup cost, in exchange for the fastest resolve + batched `take_rows` (no SQL). Storage stays ~comparable (legacy blob volume moves into many small rows; the `key` column delta-compresses well since it increments by 1 within a transcript).
- **Backward compatibility:** schema-driven format detection keeps legacy `predictions` datasets working; rollout is per-entity rebuild.
- **Genome-wide TTN** (~34,350 aa) and any future longer isoform fit trivially in the 32-bit position field — no per-chromosome sizing assumptions.
- **Follow-ups (orthogonal, not in this plan):** quantize `score` `f32 → u8/u16` (shrinks payload + dataset further; must re-verify parity to VEP's printed precision); apply the same position-keying to any other per-transcript blob stores.
