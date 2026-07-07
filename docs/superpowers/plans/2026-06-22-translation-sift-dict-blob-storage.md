# translation_sift de-interleaved fixed-divisor blob + decode UDF — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Shrink `translation_sift.lance` ~67% losslessly by changing the SIFT/PolyPhen blob's internal byte layout (de-interleaved + fixed-divisor index), keeping the 2-column schema and fast single-`take` lookup, and add a reusable DataFusion UDF so other tools can decode the cache.

**Architecture:** The `key:UInt64, sift:Binary, poly:Binary` schema, the key BTree index, and the runtime lookup path are unchanged. Only the bytes inside each blob change from interleaved `[aa][pred][score f32]` (6B) to de-interleaved `[aa…][pred…][score_idx…]` where `score_idx` is a fixed-divisor integer (sift `u8`=round(score×100), poly `u16`=round(score×1000)). A schema-metadata flag (`bio.vep.sift_blob_version`) selects v1/v2 at decode time. One shared decode core in `cache_common.rs` is used by both the annotation runtime and the new public `vep_decode_sift_predictions` ScalarUDF.

**Tech Stack:** Rust (edition 2024), DataFusion 53.0.0, Arrow 58.0.0, Lance 7.0.0 (`lance`, `lance-index`), cargo workspace crate `datafusion/bio-function-vep`.

## Global Constraints

- DataFusion `=53.0.0`, Arrow `58.0.0`, Rust edition 2024 (do not bump).
- Lance file format **2.1** (production); do **not** migrate to 2.2.
- SIFT score index width `u8`, divisor `100.0`; PolyPhen index width `u16`, divisor `1000.0` — exact verbatim values.
- Schema metadata key: `bio.vep.sift_blob_version`, value `"2"` for the new format; absent ⇒ v1.
- Losslessness is **bit-exact f32** and is enforced at build time — a score not representable as `k/divisor` must fail the build, never silently round.
- `cargo clippy --all-targets --all-features -- -D warnings` and `cargo fmt -- --check` must pass.
- The v1 interleaved-f32 decode path must remain for back-compat with un-migrated caches.

---

## File Structure

- `datafusion/bio-function-vep/src/cache_common.rs` — **shared decode/encode core** (v1 + new v2 codec, version enum, metadata helper). Single source of truth.
- `datafusion/bio-function-vep/src/lance_cache/build.rs` — build emits v2 blobs; schema carries the version flag.
- `datafusion/bio-function-vep/src/lance_cache/write.rs` — build finalize: compact + cleanup old versions.
- `datafusion/bio-function-vep/src/annotate_provider.rs` — runtime decode threads the blob version into `position_predictions_from_batch`.
- `datafusion/bio-function-vep/src/sift_decode.rs` — **new**, public `ScalarUDFImpl` wrapping the shared decode core.
- `datafusion/bio-function-vep/src/lib.rs` — register the new module + re-export.

---

## Task 1: Shared v2 codec in `cache_common.rs`

**Files:**
- Modify: `datafusion/bio-function-vep/src/cache_common.rs` (add after `deserialize_position_predictions`, ~line 178)
- Test: same file, `#[cfg(test)] mod tests` (add if absent, else append)

**Interfaces:**
- Consumes: `CompactPrediction { position: i32, amino_acid: u8, prediction: u8, score: f32 }`, `CachedPredictions { sift: Vec<CompactPrediction>, polyphen: Vec<CompactPrediction> }` from `crate::transcript_consequence`.
- Produces:
  - `pub(crate) const SIFT_BLOB_VERSION_KEY: &str = "bio.vep.sift_blob_version";`
  - `pub(crate) struct PredictorCodec { pub idx_width: usize, pub divisor: f32 }`
  - `pub(crate) const SIFT_CODEC: PredictorCodec` (idx_width 1, divisor 100.0)
  - `pub(crate) const POLY_CODEC: PredictorCodec` (idx_width 2, divisor 1000.0)
  - `pub(crate) enum SiftBlobVersion { V1Interleaved, V2DivIndex }` (derives `Clone, Copy, PartialEq, Eq, Debug`)
  - `pub(crate) fn serialize_position_entries_v2(entries: &[CompactPrediction], codec: PredictorCodec) -> Result<Vec<u8>>`
  - `pub(crate) fn deserialize_position_entries_v2(position: i32, data: &[u8], codec: PredictorCodec) -> Result<Vec<CompactPrediction>>`
  - `pub(crate) fn deserialize_position_predictions_versioned(position: i32, sift: &[u8], poly: &[u8], version: SiftBlobVersion) -> Result<CachedPredictions>`
  - `pub(crate) fn sift_blob_version_from_metadata(meta: &std::collections::HashMap<String, String>) -> SiftBlobVersion`

- [ ] **Step 1: Write the failing tests**

Append to `cache_common.rs`:

```rust
#[cfg(test)]
mod v2_codec_tests {
    use super::*;
    use crate::transcript_consequence::CompactPrediction;

    fn cp(aa: u8, pred: u8, score: f32) -> CompactPrediction {
        CompactPrediction { position: 7, amino_acid: aa, prediction: pred, score }
    }

    #[test]
    fn v2_sift_roundtrip_bit_exact() {
        let entries = vec![cp(0, 1, 0.0), cp(3, 0, 0.02), cp(19, 1, 1.0), cp(7, 0, 0.55)];
        let bytes = serialize_position_entries_v2(&entries, SIFT_CODEC).unwrap();
        // de-interleaved: n*(1+1+1) bytes
        assert_eq!(bytes.len(), entries.len() * 3);
        let back = deserialize_position_entries_v2(7, &bytes, SIFT_CODEC).unwrap();
        assert_eq!(back.len(), entries.len());
        for (a, b) in entries.iter().zip(back.iter()) {
            assert_eq!(a.amino_acid, b.amino_acid);
            assert_eq!(a.prediction, b.prediction);
            assert_eq!(a.score.to_bits(), b.score.to_bits(), "score not bit-exact");
            assert_eq!(b.position, 7);
        }
    }

    #[test]
    fn v2_poly_roundtrip_bit_exact_u16() {
        let entries = vec![cp(0, 2, 0.998), cp(5, 0, 0.001), cp(12, 1, 0.5)];
        let bytes = serialize_position_entries_v2(&entries, POLY_CODEC).unwrap();
        assert_eq!(bytes.len(), entries.len() * 4); // 1+1+2
        let back = deserialize_position_entries_v2(7, &bytes, POLY_CODEC).unwrap();
        for (a, b) in entries.iter().zip(back.iter()) {
            assert_eq!(a.score.to_bits(), b.score.to_bits());
        }
    }

    #[test]
    fn v2_empty_blob() {
        let bytes = serialize_position_entries_v2(&[], SIFT_CODEC).unwrap();
        assert!(bytes.is_empty());
        assert!(deserialize_position_entries_v2(1, &bytes, SIFT_CODEC).unwrap().is_empty());
    }

    #[test]
    fn v2_rejects_off_grid_score() {
        // 0.005 is not k/100 → must error rather than silently round.
        let entries = vec![cp(0, 1, 0.005)];
        assert!(serialize_position_entries_v2(&entries, SIFT_CODEC).is_err());
    }

    #[test]
    fn version_from_metadata() {
        let mut m = std::collections::HashMap::new();
        assert_eq!(sift_blob_version_from_metadata(&m), SiftBlobVersion::V1Interleaved);
        m.insert(SIFT_BLOB_VERSION_KEY.to_string(), "2".to_string());
        assert_eq!(sift_blob_version_from_metadata(&m), SiftBlobVersion::V2DivIndex);
    }

    #[test]
    fn versioned_dispatch_v2() {
        let sift = serialize_position_entries_v2(&[cp(1, 0, 0.04)], SIFT_CODEC).unwrap();
        let poly = serialize_position_entries_v2(&[cp(1, 2, 0.7)], POLY_CODEC).unwrap();
        let out = deserialize_position_predictions_versioned(9, &sift, &poly, SiftBlobVersion::V2DivIndex).unwrap();
        assert_eq!(out.sift[0].score.to_bits(), 0.04f32.to_bits());
        assert_eq!(out.polyphen[0].score.to_bits(), 0.7f32.to_bits());
        assert_eq!(out.sift[0].position, 9);
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep v2_codec_tests 2>&1 | tail -20`
Expected: FAIL — `serialize_position_entries_v2` etc. not found (compile error).

- [ ] **Step 3: Implement the codec**

Insert after `deserialize_position_predictions` (~line 178) in `cache_common.rs`:

```rust
/// Schema-metadata key marking the position-sliced SIFT blob layout version.
pub(crate) const SIFT_BLOB_VERSION_KEY: &str = "bio.vep.sift_blob_version";

/// Per-predictor parameters for the v2 de-interleaved fixed-divisor layout.
#[derive(Clone, Copy)]
pub(crate) struct PredictorCodec {
    /// Bytes per score index: 1 (`u8`, SIFT) or 2 (`u16` LE, PolyPhen).
    pub idx_width: usize,
    /// Score = index / divisor. 100.0 (SIFT, 0.01 grid) or 1000.0 (PolyPhen, 0.001 grid).
    pub divisor: f32,
}

pub(crate) const SIFT_CODEC: PredictorCodec = PredictorCodec { idx_width: 1, divisor: 100.0 };
pub(crate) const POLY_CODEC: PredictorCodec = PredictorCodec { idx_width: 2, divisor: 1000.0 };

/// On-disk layout of the position-sliced SIFT/PolyPhen `Binary` blobs.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum SiftBlobVersion {
    /// Interleaved 6-byte entries `[aa u8][pred u8][score f32 LE]`.
    V1Interleaved,
    /// De-interleaved `[aa u8 × n][pred u8 × n][score_idx × n]` with fixed-divisor index.
    V2DivIndex,
}

pub(crate) fn sift_blob_version_from_metadata(
    meta: &std::collections::HashMap<String, String>,
) -> SiftBlobVersion {
    match meta.get(SIFT_BLOB_VERSION_KEY).map(String::as_str) {
        Some("2") => SiftBlobVersion::V2DivIndex,
        _ => SiftBlobVersion::V1Interleaved,
    }
}

/// Serialize one protein position's predictions for a single predictor in the v2
/// de-interleaved fixed-divisor layout: `[aa × n][pred × n][score_idx × n]`.
/// Returns an error if any score is not bit-exactly representable as `k / divisor`
/// or the index does not fit `idx_width` (preserves the lossless guarantee).
pub(crate) fn serialize_position_entries_v2(
    entries: &[CompactPrediction],
    codec: PredictorCodec,
) -> Result<Vec<u8>> {
    let n = entries.len();
    let mut buf = Vec::with_capacity(n * (2 + codec.idx_width));
    for p in entries {
        buf.push(p.amino_acid);
    }
    for p in entries {
        buf.push(p.prediction);
    }
    let max = if codec.idx_width == 1 { u8::MAX as i64 } else { u16::MAX as i64 };
    for p in entries {
        let idx = (p.score * codec.divisor).round() as i64;
        if idx < 0 || idx > max {
            return Err(DataFusionError::Execution(format!(
                "SIFT/PolyPhen score {} out of range for divisor {}",
                p.score, codec.divisor
            )));
        }
        if (idx as f32 / codec.divisor).to_bits() != p.score.to_bits() {
            return Err(DataFusionError::Execution(format!(
                "SIFT/PolyPhen score {} is not bit-exactly representable as k/{}",
                p.score, codec.divisor
            )));
        }
        if codec.idx_width == 1 {
            buf.push(idx as u8);
        } else {
            buf.extend_from_slice(&(idx as u16).to_le_bytes());
        }
    }
    Ok(buf)
}

/// Inverse of [`serialize_position_entries_v2`].
pub(crate) fn deserialize_position_entries_v2(
    position: i32,
    data: &[u8],
    codec: PredictorCodec,
) -> Result<Vec<CompactPrediction>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    let stride = 2 + codec.idx_width;
    if data.len() % stride != 0 {
        return Err(DataFusionError::Execution(format!(
            "v2 SIFT payload length {} is not a multiple of stride {stride}",
            data.len()
        )));
    }
    let n = data.len() / stride;
    let aa = &data[0..n];
    let pred = &data[n..2 * n];
    let idx = &data[2 * n..];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let raw = if codec.idx_width == 1 {
            idx[i] as u16
        } else {
            u16::from_le_bytes([idx[2 * i], idx[2 * i + 1]])
        };
        out.push(CompactPrediction {
            position,
            amino_acid: aa[i],
            prediction: pred[i],
            score: raw as f32 / codec.divisor,
        });
    }
    Ok(out)
}

/// Decode both predictor blobs for one position, selecting the layout by version.
pub(crate) fn deserialize_position_predictions_versioned(
    position: i32,
    sift: &[u8],
    poly: &[u8],
    version: SiftBlobVersion,
) -> Result<CachedPredictions> {
    match version {
        SiftBlobVersion::V1Interleaved => deserialize_position_predictions(position, sift, poly),
        SiftBlobVersion::V2DivIndex => Ok(CachedPredictions {
            sift: deserialize_position_entries_v2(position, sift, SIFT_CODEC)?,
            polyphen: deserialize_position_entries_v2(position, poly, POLY_CODEC)?,
        }),
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep v2_codec_tests 2>&1 | tail -20`
Expected: PASS (6 tests).

- [ ] **Step 5: Lint + commit**

```bash
cargo fmt -p datafusion-bio-function-vep
cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -5
git add datafusion/bio-function-vep/src/cache_common.rs
git commit -m "feat(vep): v2 de-interleaved fixed-divisor SIFT blob codec"
```

---

## Task 2: Build emits v2 blobs + schema version flag

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/build.rs` — `compact_translation_sift_position_schema` (~1391) and `append_translation_sift_position_rows` (~1461)
- Test: `build.rs` existing `#[cfg(test)]` module (there is a test near line 2161 using `deserialize_position_predictions`)

**Interfaces:**
- Consumes: `serialize_position_entries_v2`, `SIFT_CODEC`, `POLY_CODEC`, `SIFT_BLOB_VERSION_KEY` (Task 1).
- Produces: the on-disk `translation_sift` dataset now stores v2 blobs and the schema metadata `bio.vep.sift_blob_version = "2"`.

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)]` module in `build.rs`:

```rust
#[test]
fn append_translation_sift_position_rows_emits_v2() {
    use crate::cache_common::{deserialize_position_entries_v2, SIFT_CODEC, POLY_CODEC};
    use crate::transcript_consequence::{CachedPredictions, CompactPrediction};
    let preds = CachedPredictions {
        sift: vec![
            CompactPrediction { position: 1, amino_acid: 0, prediction: 1, score: 0.03 },
            CompactPrediction { position: 1, amino_acid: 5, prediction: 0, score: 1.0 },
        ],
        polyphen: vec![
            CompactPrediction { position: 1, amino_acid: 0, prediction: 2, score: 0.998 },
        ],
    };
    let mut keys = Vec::new();
    let mut sift_blobs = Vec::new();
    let mut poly_blobs = Vec::new();
    append_translation_sift_position_rows(42, &preds, &mut keys, &mut sift_blobs, &mut poly_blobs).unwrap();
    assert_eq!(keys, vec![(42u64 << 32) | 1]);
    // v2 sift: 2 entries * 3 bytes
    assert_eq!(sift_blobs[0].len(), 6);
    let s = deserialize_position_entries_v2(1, &sift_blobs[0], SIFT_CODEC).unwrap();
    assert_eq!(s[1].score.to_bits(), 1.0f32.to_bits());
    let p = deserialize_position_entries_v2(1, &poly_blobs[0], POLY_CODEC).unwrap();
    assert_eq!(p[0].score.to_bits(), 0.998f32.to_bits());
}

#[test]
fn sift_position_schema_carries_v2_flag() {
    let schema = compact_translation_sift_position_schema("ensembl");
    assert_eq!(
        schema.metadata().get(crate::cache_common::SIFT_BLOB_VERSION_KEY).map(String::as_str),
        Some("2")
    );
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --all-features append_translation_sift_position_rows_emits_v2 sift_position_schema_carries_v2_flag 2>&1 | tail -20`
Expected: FAIL — blob lengths are 6-byte interleaved (sift_blobs[0].len()==12), schema flag absent.

- [ ] **Step 3: Implement the build change**

In `append_translation_sift_position_rows` (build.rs ~1493), replace the two `serialize_position_entries` calls:

```rust
        keys.push(((uid as u64) << 32) | (next_pos as u64));
        sift_blobs.push(crate::cache_common::serialize_position_entries_v2(
            &sift[s_start..si],
            crate::cache_common::SIFT_CODEC,
        )?);
        poly_blobs.push(crate::cache_common::serialize_position_entries_v2(
            &poly[p_start..pi],
            crate::cache_common::POLY_CODEC,
        )?);
```

In `compact_translation_sift_position_schema` (build.rs ~1391), add the version flag to the schema metadata after `with_cache_source_metadata`:

```rust
fn compact_translation_sift_position_schema(source_type: &str) -> Schema {
    let meta = lance_field_metadata();
    let key = Field::new("key", DataType::UInt64, false).with_metadata(meta.clone());
    let sift = Field::new("sift", DataType::Binary, false).with_metadata(meta.clone());
    let poly = Field::new("poly", DataType::Binary, false).with_metadata(meta);
    let schema = with_cache_source_metadata(&Schema::new(vec![key, sift, poly]), source_type);
    let mut md = schema.metadata().clone();
    md.insert(
        crate::cache_common::SIFT_BLOB_VERSION_KEY.to_string(),
        "2".to_string(),
    );
    Schema::new_with_metadata(schema.fields().clone(), md)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep --all-features append_translation_sift_position_rows_emits_v2 sift_position_schema_carries_v2_flag 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Fix the pre-existing build round-trip test (if any)**

The existing test near build.rs:2193 builds blobs with `append_translation_sift_position_rows` then decodes with `deserialize_position_predictions` (v1). Update it to decode via `deserialize_position_predictions_versioned(pos, &sift, &poly, SiftBlobVersion::V2DivIndex)`.

Run: `cargo test -p datafusion-bio-function-vep --all-features --lib lance_cache::build 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 6: Lint + commit**

```bash
cargo fmt -p datafusion-bio-function-vep
cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -5
git add datafusion/bio-function-vep/src/lance_cache/build.rs
git commit -m "feat(vep): build translation_sift as v2 blobs + schema version flag"
```

---

## Task 3: Runtime decode branches on blob version

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` — `position_predictions_from_batch` (~3188), `PositionSlicedLanceSiftStore` struct + its `get_position_predictions`, and `load_lance_sift_prediction_store_for_chrom` (~3223)

**Interfaces:**
- Consumes: `SiftBlobVersion`, `sift_blob_version_from_metadata`, `deserialize_position_predictions_versioned` (Task 1).
- Produces: runtime correctly decodes both v1 and v2 caches.

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)]` module in `annotate_provider.rs` (near the existing `position_predictions_from_batch_decodes_key_and_payloads` test ~13568):

```rust
#[test]
fn position_predictions_from_batch_decodes_v2() {
    use crate::cache_common::{serialize_position_entries_v2, SiftBlobVersion, SIFT_CODEC, POLY_CODEC};
    use crate::transcript_consequence::CompactPrediction;
    use datafusion::arrow::array::{BinaryArray, UInt64Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    let key = (3u64 << 32) | 11;
    let sift = serialize_position_entries_v2(
        &[CompactPrediction { position: 11, amino_acid: 2, prediction: 1, score: 0.07 }],
        SIFT_CODEC,
    ).unwrap();
    let poly = serialize_position_entries_v2(
        &[CompactPrediction { position: 11, amino_acid: 2, prediction: 2, score: 0.951 }],
        POLY_CODEC,
    ).unwrap();
    let schema = Arc::new(Schema::new(vec![
        Field::new("key", DataType::UInt64, false),
        Field::new("sift", DataType::Binary, false),
        Field::new("poly", DataType::Binary, false),
    ]));
    let batch = datafusion::arrow::record_batch::RecordBatch::try_new(
        schema,
        vec![
            Arc::new(UInt64Array::from(vec![key])),
            Arc::new(BinaryArray::from(vec![sift.as_slice()])),
            Arc::new(BinaryArray::from(vec![poly.as_slice()])),
        ],
    ).unwrap();
    let out = position_predictions_from_batch(&batch, SiftBlobVersion::V2DivIndex).unwrap();
    let cached = out.get(&key).unwrap();
    assert_eq!(cached.sift[0].score.to_bits(), 0.07f32.to_bits());
    assert_eq!(cached.polyphen[0].score.to_bits(), 0.951f32.to_bits());
    assert_eq!(cached.sift[0].position, 11);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --all-features position_predictions_from_batch_decodes_v2 2>&1 | tail -20`
Expected: FAIL — `position_predictions_from_batch` takes 1 arg, not 2.

- [ ] **Step 3: Thread the version through**

In `annotate_provider.rs`, change `position_predictions_from_batch` signature (~3188) and the decode call:

```rust
fn position_predictions_from_batch(
    batch: &RecordBatch,
    version: crate::cache_common::SiftBlobVersion,
) -> Result<HashMap<u64, CachedPredictions>> {
```

and replace the `deserialize_position_predictions(position, sift, poly)` call (~3216) with:

```rust
        let cached = crate::cache_common::deserialize_position_predictions_versioned(
            position, sift, poly, version,
        )?;
```

Add a `blob_version` field to `PositionSlicedLanceSiftStore` (find its `struct` definition; add `blob_version: crate::cache_common::SiftBlobVersion,`) and pass it in its `get_position_predictions` where it calls `position_predictions_from_batch(&batch)` → `position_predictions_from_batch(&batch, self.blob_version)`.

In `load_lance_sift_prediction_store_for_chrom` (~3245), capture the version from the schema and set the field:

```rust
        let blob_version = crate::cache_common::sift_blob_version_from_metadata(
            &schema.metadata().clone().into_iter().collect(),
        );
        return Ok(Some(Arc::new(PositionSlicedLanceSiftStore {
            lookup: Arc::new(lookup),
            cache: Arc::new(Mutex::new(HashMap::new())),
            blob_version,
        }) as SiftPredictionStoreRef));
```

(`schema` here is the `arrow::datatypes::Schema` returned by `read_lance_dataset_schema`; its `.metadata()` is a `HashMap<String, String>` — pass it directly if the type already matches, dropping the `.clone().into_iter().collect()` adapter.)

Update the other call site in the existing test `position_predictions_from_batch_decodes_key_and_payloads` to pass `SiftBlobVersion::V1Interleaved`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p datafusion-bio-function-vep --all-features position_predictions_from_batch 2>&1 | tail -20`
Expected: PASS (both v1 and v2 tests).

- [ ] **Step 5: Lint + commit**

```bash
cargo fmt -p datafusion-bio-function-vep
cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -5
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): runtime decodes v2 SIFT blobs (version from schema metadata)"
```

---

## Task 4: Decode UDF (`vep_decode_sift_predictions`)

**Files:**
- Create: `datafusion/bio-function-vep/src/sift_decode.rs`
- Modify: `datafusion/bio-function-vep/src/lib.rs` (add `pub mod sift_decode;` near line 62 and re-export)
- Test: in `sift_decode.rs` `#[cfg(test)]`

**Interfaces:**
- Consumes: `deserialize_position_entries_v2`, `SIFT_CODEC`, `POLY_CODEC` (Task 1).
- Produces:
  - `pub fn vep_decode_sift_predictions_udf() -> datafusion::logical_expr::ScalarUDF`
  - `pub fn register_vep_sift_functions(ctx: &datafusion::prelude::SessionContext)`
  - UDF signature: `vep_decode_sift_predictions(blob: Binary, predictor: Utf8) -> List<Struct<amino_acid: UInt8, prediction: UInt8, score: Float32>>`. `predictor` ∈ {`"sift"`, `"polyphen"`/`"poly"`}.

- [ ] **Step 1: Write the failing test**

Create `sift_decode.rs` with only the test first:

```rust
//! Public DataFusion UDF that decodes position-sliced SIFT/PolyPhen `Binary`
//! blobs (v2 de-interleaved fixed-divisor layout) into structured rows, so
//! tools other than the VEP engine can read `translation_sift.lance`.

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Array, BinaryArray, Float32Array, ListArray, StringArray, StructArray, UInt8Array};
    use datafusion::logical_expr::ColumnarValue;
    use std::sync::Arc;

    #[test]
    fn decodes_sift_blob_to_struct_list() {
        use crate::cache_common::{serialize_position_entries_v2, SIFT_CODEC};
        use crate::transcript_consequence::CompactPrediction;
        let blob = serialize_position_entries_v2(
            &[
                CompactPrediction { position: 0, amino_acid: 4, prediction: 1, score: 0.02 },
                CompactPrediction { position: 0, amino_acid: 9, prediction: 0, score: 1.0 },
            ],
            SIFT_CODEC,
        ).unwrap();
        let udf = vep_decode_sift_predictions_udf();
        let args = vec![
            ColumnarValue::Array(Arc::new(BinaryArray::from(vec![blob.as_slice()]))),
            ColumnarValue::Array(Arc::new(StringArray::from(vec!["sift"]))),
        ];
        let out = invoke_for_test(&udf, args, 1);
        let list = out.as_any().downcast_ref::<ListArray>().unwrap();
        let inner = list.value(0);
        let s = inner.as_any().downcast_ref::<StructArray>().unwrap();
        assert_eq!(s.len(), 2);
        let aa = s.column(0).as_any().downcast_ref::<UInt8Array>().unwrap();
        let score = s.column(2).as_any().downcast_ref::<Float32Array>().unwrap();
        assert_eq!(aa.value(0), 4);
        assert_eq!(score.value(1).to_bits(), 1.0f32.to_bits());
    }
}
```

> NOTE: `invoke_for_test` is a thin helper you add in Step 3 that builds a
> `ScalarFunctionArgs` and calls `ScalarUDFImpl::invoke_with_args`, returning the
> result `ArrayRef`. Defining it next to the impl keeps the test honest about the
> DF 53 calling convention.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --all-features sift_decode 2>&1 | tail -20`
Expected: FAIL — module/types not found (also requires `pub mod sift_decode;` in lib.rs to compile; add it now).

- [ ] **Step 3: Implement the UDF**

Prepend to `sift_decode.rs` (above the test module):

```rust
use std::any::Any;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, BinaryArray, Float32Builder, ListBuilder, StructBuilder, UInt8Builder,
};
use datafusion::arrow::datatypes::{DataType, Field, Fields};
use datafusion::common::{exec_err, Result};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, TypeSignature,
    Volatility,
};

use crate::cache_common::{deserialize_position_entries_v2, PredictorCodec, POLY_CODEC, SIFT_CODEC};

fn item_fields() -> Fields {
    Fields::from(vec![
        Field::new("amino_acid", DataType::UInt8, false),
        Field::new("prediction", DataType::UInt8, false),
        Field::new("score", DataType::Float32, false),
    ])
}

fn return_field() -> Field {
    Field::new(
        "item",
        DataType::Struct(item_fields()),
        true,
    )
}

#[derive(Debug)]
pub struct SiftDecodeUdf {
    signature: Signature,
}

impl Default for SiftDecodeUdf {
    fn default() -> Self {
        Self {
            signature: Signature::new(
                TypeSignature::Exact(vec![DataType::Binary, DataType::Utf8]),
                Volatility::Immutable,
            ),
        }
    }
}

impl ScalarUDFImpl for SiftDecodeUdf {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn name(&self) -> &str {
        "vep_decode_sift_predictions"
    }
    fn signature(&self) -> &Signature {
        &self.signature
    }
    fn return_type(&self, _arg_types: &[DataType]) -> Result<DataType> {
        Ok(DataType::List(Arc::new(return_field())))
    }
    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let arrays = ColumnarValue::values_to_arrays(&args.args)?;
        let blobs = arrays[0]
            .as_any()
            .downcast_ref::<BinaryArray>()
            .ok_or_else(|| datafusion::common::DataFusionError::Execution(
                "vep_decode_sift_predictions: arg 0 must be Binary".into(),
            ))?;
        let predictor = datafusion::common::cast::as_string_array(&arrays[1])?;

        let values_builder = StructBuilder::new(
            item_fields(),
            vec![
                Box::new(UInt8Builder::new()),
                Box::new(UInt8Builder::new()),
                Box::new(Float32Builder::new()),
            ],
        );
        let mut list = ListBuilder::new(values_builder);

        for row in 0..blobs.len() {
            if blobs.is_null(row) {
                list.append_null();
                continue;
            }
            let codec = match predictor.value(row) {
                "sift" => SIFT_CODEC,
                "polyphen" | "poly" => POLY_CODEC,
                other => return exec_err!("unknown predictor '{other}' (expected sift|polyphen)"),
            };
            let entries = deserialize_position_entries_v2(0, blobs.value(row), codec)?;
            let sb = list.values();
            for e in &entries {
                sb.field_builder::<UInt8Builder>(0).unwrap().append_value(e.amino_acid);
                sb.field_builder::<UInt8Builder>(1).unwrap().append_value(e.prediction);
                sb.field_builder::<Float32Builder>(2).unwrap().append_value(e.score);
                sb.append(true);
            }
            list.append(true);
        }
        Ok(ColumnarValue::Array(Arc::new(list.finish())))
    }
}

/// Build the `vep_decode_sift_predictions` scalar UDF.
pub fn vep_decode_sift_predictions_udf() -> ScalarUDF {
    ScalarUDF::from(SiftDecodeUdf::default())
}

/// Register the VEP SIFT decode UDF on a session context.
pub fn register_vep_sift_functions(ctx: &datafusion::prelude::SessionContext) {
    ctx.register_udf(vep_decode_sift_predictions_udf());
}

#[cfg(test)]
fn invoke_for_test(udf: &ScalarUDF, args: Vec<ColumnarValue>, number_rows: usize) -> ArrayRef {
    let arg_fields = args
        .iter()
        .map(|a| Arc::new(Field::new("a", a.data_type(), true)))
        .collect::<Vec<_>>();
    let return_field = Arc::new(Field::new("out", udf.return_type(&[DataType::Binary, DataType::Utf8]).unwrap(), true));
    let res = udf
        .inner()
        .invoke_with_args(ScalarFunctionArgs { args, arg_fields, number_rows, return_field })
        .unwrap();
    match res {
        ColumnarValue::Array(a) => a,
        ColumnarValue::Scalar(s) => s.to_array().unwrap(),
    }
}
```

Add to `lib.rs` (after `pub mod schema_contract;` / near line 62):

```rust
pub mod sift_decode;
```

and a re-export near the other `pub use` lines:

```rust
pub use sift_decode::{register_vep_sift_functions, vep_decode_sift_predictions_udf};
```

> NOTE on DF 53 API: `ScalarFunctionArgs` fields are `{ args, arg_fields, number_rows, return_field }`. If a field name differs in the pinned version, run `cargo doc -p datafusion --open` or grep the datafusion source under `~/.cargo` for `struct ScalarFunctionArgs` and adjust the struct literal in `invoke_for_test` and the impl accordingly. Do not change the UDF behavior.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep --all-features sift_decode 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Lint + commit**

```bash
cargo fmt -p datafusion-bio-function-vep
cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -5
git add datafusion/bio-function-vep/src/sift_decode.rs datafusion/bio-function-vep/src/lib.rs
git commit -m "feat(vep): public vep_decode_sift_predictions ScalarUDF"
```

---

## Task 5: Build finalize — compaction + version cleanup

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/write.rs` (the function that writes the SIFT dataset and creates the `sift_key_btree_idx`, ~lines 160-210)

**Interfaces:**
- Consumes: an open `lance::Dataset` for the translation_sift path.
- Produces: `async fn compact_and_cleanup_sift_dataset(path: &Path) -> Result<()>` called at the end of the SIFT write/index step.

- [ ] **Step 1: Confirm the Lance Rust compaction API**

Run: `grep -rn "compact_files\|cleanup_old_versions\|CompactionOptions" ~/.cargo/registry/src/*/lance-*/src 2>/dev/null | head` (or `cargo doc -p lance`).
Record the exact paths (expected: `lance::dataset::optimize::compact_files(&mut dataset, options, None).await` and `dataset.cleanup_old_versions(duration, None, None).await` or similar). If the API differs, adapt the call in Step 3 — the behavior (merge fragments to few, drop old versions) is the requirement.

- [ ] **Step 2: Write the failing test**

Add to `write.rs` `#[cfg(test)]` (there are existing tests opening datasets ~line 258):

```rust
#[tokio::test(flavor = "multi_thread")]
async fn compaction_reduces_fragment_count() {
    use datafusion::arrow::array::{BinaryArray, UInt64Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.lance");
    let schema = Arc::new(Schema::new(vec![
        Field::new("key", DataType::UInt64, false),
        Field::new("sift", DataType::Binary, false),
    ]));
    // write several small fragments
    for i in 0..5u64 {
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![i])),
                Arc::new(BinaryArray::from(vec![b"x".as_ref()])),
            ],
        ).unwrap();
        let reader = Box::new(datafusion::arrow::record_batch::RecordBatchIterator::new(
            vec![Ok(batch)], schema.clone(),
        ));
        let mode = if i == 0 { lance::dataset::WriteMode::Create } else { lance::dataset::WriteMode::Append };
        lance::Dataset::write(reader, path.to_string_lossy().as_ref(), Some(lance::dataset::WriteParams { mode, ..Default::default() })).await.unwrap();
    }
    let before = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap().get_fragments().len();
    assert!(before > 1);
    compact_and_cleanup_sift_dataset(&path).await.unwrap();
    let after = lance::Dataset::open(path.to_string_lossy().as_ref()).await.unwrap().get_fragments().len();
    assert!(after < before, "expected compaction to reduce fragments {before} -> {after}");
}
```

- [ ] **Step 3: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --all-features compaction_reduces_fragment_count 2>&1 | tail -20`
Expected: FAIL — `compact_and_cleanup_sift_dataset` not found.

- [ ] **Step 4: Implement (adapt API names from Step 1)**

Add to `write.rs`:

```rust
/// Merge fragments and drop superseded versions for a freshly built SIFT dataset.
/// Field encoding metadata is carried by Lance through compaction.
pub(crate) async fn compact_and_cleanup_sift_dataset(path: &std::path::Path) -> Result<()> {
    use lance::dataset::optimize::{compact_files, CompactionOptions};
    let mut dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
        .await
        .map_err(|e| DataFusionError::Execution(format!("open for compaction: {e}")))?;
    compact_files(&mut dataset, CompactionOptions::default(), None)
        .await
        .map_err(|e| DataFusionError::Execution(format!("compact_files: {e}")))?;
    dataset
        .cleanup_old_versions(std::time::Duration::from_secs(0), Some(true), None)
        .await
        .map_err(|e| DataFusionError::Execution(format!("cleanup_old_versions: {e}")))?;
    Ok(())
}
```

Call it at the end of the SIFT dataset write+index function (after `create_index` for `sift_key_btree_idx`, ~line 210), guarded so it only runs for the translation_sift dataset.

- [ ] **Step 5: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep --all-features compaction_reduces_fragment_count 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 6: Lint + commit**

```bash
cargo fmt -p datafusion-bio-function-vep
cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -5
git add datafusion/bio-function-vep/src/lance_cache/write.rs
git commit -m "feat(vep): compact + cleanup translation_sift after build"
```

---

## Task 6: Parity gate (UDF ↔ runtime) + full regression

**Files:**
- Test: `datafusion/bio-function-vep/src/sift_decode.rs` `#[cfg(test)]`

**Interfaces:**
- Consumes: `vep_decode_sift_predictions_udf` (Task 4), `serialize_position_entries_v2`/`deserialize_position_entries_v2` (Task 1).

- [ ] **Step 1: Write the cross-check test**

Add to `sift_decode.rs` tests:

```rust
#[test]
fn udf_matches_core_decode_polyphen() {
    use crate::cache_common::{serialize_position_entries_v2, deserialize_position_entries_v2, POLY_CODEC};
    use crate::transcript_consequence::CompactPrediction;
    use datafusion::arrow::array::{BinaryArray, Float32Array, ListArray, StringArray, StructArray};
    use datafusion::logical_expr::ColumnarValue;

    let entries: Vec<CompactPrediction> = (0..19u8)
        .map(|i| CompactPrediction { position: 0, amino_acid: i, prediction: i % 3, score: (i as f32) * 0.05 })
        .collect();
    let blob = serialize_position_entries_v2(&entries, POLY_CODEC).unwrap();
    let expected = deserialize_position_entries_v2(0, &blob, POLY_CODEC).unwrap();

    let udf = vep_decode_sift_predictions_udf();
    let out = invoke_for_test(&udf, vec![
        ColumnarValue::Array(std::sync::Arc::new(BinaryArray::from(vec![blob.as_slice()]))),
        ColumnarValue::Array(std::sync::Arc::new(StringArray::from(vec!["polyphen"]))),
    ], 1);
    let list = out.as_any().downcast_ref::<ListArray>().unwrap();
    let s = list.value(0);
    let s = s.as_any().downcast_ref::<StructArray>().unwrap();
    let score = s.column(2).as_any().downcast_ref::<Float32Array>().unwrap();
    for (i, e) in expected.iter().enumerate() {
        assert_eq!(score.value(i).to_bits(), e.score.to_bits());
    }
}
```

- [ ] **Step 2: Run + verify pass**

Run: `cargo test -p datafusion-bio-function-vep --all-features sift_decode 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 3: Full crate test + lint (regression sweep)**

Run: `cargo test -p datafusion-bio-function-vep --all-features 2>&1 | tail -30`
Expected: no NEW failures vs baseline. (Baseline has known pre-existing `annotate_table_function::tests` fixture failures requiring real parquet/cache data — confirm the failure set is unchanged.)

Run: `cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -5`
Expected: clean.

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/src/sift_decode.rs
git commit -m "test(vep): UDF/runtime SIFT decode parity"
```

---

## Task 7: Manual end-to-end validation (documented, not CI)

**Not a code task** — record the procedure and outcome.

- [ ] **Step 1: Rebuild chr1 translation_sift in v2** via the existing `build_lance_chrom_context` example (see memory `sift-position-slicing-resume`):

```bash
cargo run --release --features lance-cache,cache-builder --example build_lance_chrom_context -- \
  --cache-root /Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38 \
  --output-dir /Users/mwiewior/workspace/data_vepyr/scratch_chr1_v2 \
  --chrom chr1 --cache-source-type merged --partitions 8 --overwrite
```

- [ ] **Step 2: Verify size + lossless** — `du -sh translation_sift.lance/chr1.lance` should be ~0.9–1.0 GB (vs 2.9 GB). Spot-check `vep_decode_sift_predictions` output against a v1 reference build for bit-exact scores.

- [ ] **Step 3: chr1 e2e parity** vs the VEP reference (vepyr `run_annotation_fast.py chr1 --cache merged --backend lance`), confirming all shared SIFT/PolyPhen CSQ fields match 100% and lookup timing is not regressed. Record results in the design spec / memory.

---

## Self-Review

**Spec coverage:** storage format change → Tasks 1–2; build finalize/compaction → Task 5; runtime decode (v1+v2) → Task 3; decode UDF → Task 4; parity gate → Task 6; e2e → Task 7. All spec sections covered.

**Placeholder scan:** no TBD/TODO; every code step shows full code; API-uncertainty notes (DF 53 `ScalarFunctionArgs`, Lance compaction) include exact verification commands rather than hand-waving.

**Type consistency:** `serialize_position_entries_v2`/`deserialize_position_entries_v2`/`PredictorCodec`/`SIFT_CODEC`/`POLY_CODEC`/`SiftBlobVersion`/`deserialize_position_predictions_versioned`/`sift_blob_version_from_metadata`/`SIFT_BLOB_VERSION_KEY` are defined in Task 1 and used with identical names/signatures in Tasks 2–6. UDF `vep_decode_sift_predictions_udf`/`register_vep_sift_functions` defined in Task 4, reused in Task 6.
