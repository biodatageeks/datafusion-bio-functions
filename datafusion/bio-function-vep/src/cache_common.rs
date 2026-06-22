//! Backend-neutral cache helpers shared by the annotation engine and the Lance
//! cache backend.
//!
//! These types were extracted from the (now Lance-only) `kv_cache`/`warm_cache`
//! modules so the Lance path does not depend on the deleted fjall/parquet code:
//! the `SiftPredictionStore` trait the Lance SIFT stores implement, the
//! prediction (de)serialization family used by the Lance SIFT layout and
//! builder, and the allele-frequency helpers used to select warm positions
//! during Lance cache construction.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

use datafusion::common::{DataFusionError, Result};

use crate::transcript_consequence::{CachedPredictions, CompactPrediction};

/// Allele matching predicate: `(vcf_ref, vcf_alt, allele_string) -> bool`.
pub type AlleleMatcher = fn(&str, &str, &str) -> bool;

// ---------------------------------------------------------------------------
// SIFT/PolyPhen prediction store interface
// ---------------------------------------------------------------------------

/// A store of SIFT/PolyPhen predictions. Implemented by the Lance-backed SIFT
/// stores (transcript-id-keyed and position-sliced).
pub trait SiftPredictionStore: Send + Sync {
    fn get_many(&self, transcript_ids: &[String]) -> Result<HashMap<String, CachedPredictions>>;

    /// True when this store resolves predictions by the position-sliced packed
    /// key `(transcript_uid << 32) | protein_position` rather than by
    /// transcript id. The annotation engine routes to [`Self::get_position_predictions`]
    /// instead of [`Self::get_many`] for these stores.
    fn is_position_sliced(&self) -> bool {
        false
    }

    /// Resolve position-sliced predictions for the given packed keys. Each
    /// returned [`CachedPredictions`] holds only the entries for that key's
    /// single protein position. Keys with no predictions are simply absent from
    /// the map. Transcript-id-keyed stores do not implement this.
    fn get_position_predictions(&self, _keys: &[u64]) -> Result<HashMap<u64, CachedPredictions>> {
        Err(DataFusionError::Execution(
            "get_position_predictions called on a transcript-id-keyed SIFT store".into(),
        ))
    }
}

/// Shared handle to a [`SiftPredictionStore`].
pub(crate) type SiftPredictionStoreRef = Arc<dyn SiftPredictionStore>;

// ---------------------------------------------------------------------------
// Prediction (de)serialization
// ---------------------------------------------------------------------------

/// Transcript-id-keyed blob: `[4B sift_count][4B poly_count]` then 10-byte
/// `CompactPrediction` records (`[position i32 LE][amino_acid u8][prediction u8][score f32 LE]`).
pub(crate) fn serialize_predictions(preds: &CachedPredictions) -> Vec<u8> {
    let sift_count = preds.sift.len() as u32;
    let polyphen_count = preds.polyphen.len() as u32;
    let mut buf = Vec::with_capacity(8 + (sift_count + polyphen_count) as usize * 10);

    buf.extend_from_slice(&sift_count.to_le_bytes());
    buf.extend_from_slice(&polyphen_count.to_le_bytes());

    for p in &preds.sift {
        buf.extend_from_slice(&p.position.to_le_bytes());
        buf.push(p.amino_acid);
        buf.push(p.prediction);
        buf.extend_from_slice(&p.score.to_le_bytes());
    }
    for p in &preds.polyphen {
        buf.extend_from_slice(&p.position.to_le_bytes());
        buf.push(p.amino_acid);
        buf.push(p.prediction);
        buf.extend_from_slice(&p.score.to_le_bytes());
    }
    buf
}

pub(crate) fn deserialize_predictions(data: &[u8]) -> Result<CachedPredictions> {
    if data.len() < 8 {
        return Err(DataFusionError::Execution(
            "sift entry too short".to_string(),
        ));
    }
    let sift_count = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    let polyphen_count = u32::from_le_bytes(data[4..8].try_into().unwrap()) as usize;

    let expected_len = 8 + (sift_count + polyphen_count) * 10;
    if data.len() < expected_len {
        return Err(DataFusionError::Execution(format!(
            "sift entry too short: expected {expected_len}, got {}",
            data.len()
        )));
    }

    let mut offset = 8;
    let mut sift = Vec::with_capacity(sift_count);
    for _ in 0..sift_count {
        sift.push(CompactPrediction {
            position: i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()),
            amino_acid: data[offset + 4],
            prediction: data[offset + 5],
            score: f32::from_le_bytes(data[offset + 6..offset + 10].try_into().unwrap()),
        });
        offset += 10;
    }

    let mut polyphen = Vec::with_capacity(polyphen_count);
    for _ in 0..polyphen_count {
        polyphen.push(CompactPrediction {
            position: i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()),
            amino_acid: data[offset + 4],
            prediction: data[offset + 5],
            score: f32::from_le_bytes(data[offset + 6..offset + 10].try_into().unwrap()),
        });
        offset += 10;
    }

    // Already sorted (stored sorted during ingestion)
    Ok(CachedPredictions { sift, polyphen })
}

/// Serialize one protein position's predictions for a single predictor, *without*
/// the position field — in the position-sliced Lance layout the position is
/// implicit from the row key, so each entry is 6 bytes:
///   [amino_acid u8][prediction u8][score f32 LE]
/// Entry count is recovered from the byte length (`len / 6`). The caller groups
/// entries by position; entries are expected pre-sorted by `amino_acid`.
pub(crate) fn serialize_position_entries(entries: &[CompactPrediction]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(entries.len() * 6);
    for p in entries {
        buf.push(p.amino_acid);
        buf.push(p.prediction);
        buf.extend_from_slice(&p.score.to_le_bytes());
    }
    buf
}

/// Inverse of [`serialize_position_entries`]: reconstruct entries for one
/// position from a 6-bytes-per-entry payload.
pub(crate) fn deserialize_position_entries(
    position: i32,
    data: &[u8],
) -> Result<Vec<CompactPrediction>> {
    if !data.len().is_multiple_of(6) {
        return Err(DataFusionError::Execution(format!(
            "position prediction payload length {} is not a multiple of 6",
            data.len()
        )));
    }
    let count = data.len() / 6;
    let mut out = Vec::with_capacity(count);
    let mut offset = 0;
    for _ in 0..count {
        out.push(CompactPrediction {
            position,
            amino_acid: data[offset],
            prediction: data[offset + 1],
            score: f32::from_le_bytes(data[offset + 2..offset + 6].try_into().unwrap()),
        });
        offset += 6;
    }
    Ok(out)
}

/// Build a single-position [`CachedPredictions`] from its `sift`/`poly` payloads.
pub(crate) fn deserialize_position_predictions(
    position: i32,
    sift: &[u8],
    poly: &[u8],
) -> Result<CachedPredictions> {
    Ok(CachedPredictions {
        sift: deserialize_position_entries(position, sift)?,
        polyphen: deserialize_position_entries(position, poly)?,
    })
}

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

pub(crate) const SIFT_CODEC: PredictorCodec = PredictorCodec {
    idx_width: 1,
    divisor: 100.0,
};
pub(crate) const POLY_CODEC: PredictorCodec = PredictorCodec {
    idx_width: 2,
    divisor: 1000.0,
};

/// On-disk layout of the position-sliced SIFT/PolyPhen `Binary` blobs.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum SiftBlobVersion {
    /// Interleaved 6-byte entries `[aa u8][pred u8][score f32 LE]`.
    V1Interleaved,
    /// De-interleaved `[aa u8 × n][pred u8 × n][score_idx × n]` with fixed-divisor index.
    V2DivIndex,
}

pub(crate) fn sift_blob_version_from_metadata(meta: &HashMap<String, String>) -> SiftBlobVersion {
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
    let max = if codec.idx_width == 1 {
        u8::MAX as i64
    } else {
        u16::MAX as i64
    };
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
    if !data.len().is_multiple_of(stride) {
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

// ---------------------------------------------------------------------------
// Allele-frequency helpers (warm-position selection for Lance cache build)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct FrequencyFields<'a> {
    pub minor_allele_freq: Option<f64>,
    pub af: Option<&'a str>,
    pub gnomadg: Option<&'a str>,
    pub gnomade: Option<&'a str>,
}

#[derive(Debug, Clone, Copy)]
pub struct PositionFrequency {
    pub start: i64,
    pub max_af: f64,
}

pub fn max_af_from_pairs(value: Option<&str>) -> f64 {
    let Some(value) = value else {
        return 0.0;
    };

    value
        .split(',')
        .filter_map(|part| {
            let raw = part
                .rsplit_once(':')
                .map(|(_, freq)| freq)
                .unwrap_or(part)
                .trim();
            raw.parse::<f64>().ok().filter(|v| v.is_finite())
        })
        .fold(0.0, f64::max)
}

pub fn max_global_af(fields: &FrequencyFields<'_>) -> f64 {
    [
        fields
            .minor_allele_freq
            .filter(|value| value.is_finite())
            .unwrap_or(0.0),
        max_af_from_pairs(fields.af),
        max_af_from_pairs(fields.gnomadg),
        max_af_from_pairs(fields.gnomade),
    ]
    .into_iter()
    .fold(0.0, f64::max)
}

pub fn select_warm_positions<I>(rows: I, af_threshold: f64, position_radius: i64) -> BTreeSet<i64>
where
    I: IntoIterator<Item = PositionFrequency>,
{
    let mut max_af_by_position: BTreeMap<i64, f64> = BTreeMap::new();
    for row in rows {
        let max_af = if row.max_af.is_finite() {
            row.max_af
        } else {
            0.0
        };
        max_af_by_position
            .entry(row.start)
            .and_modify(|existing| *existing = existing.max(max_af))
            .or_insert(max_af);
    }

    let source_positions: BTreeSet<i64> = max_af_by_position.keys().copied().collect();
    let radius = position_radius.max(0);
    let mut warm_positions = BTreeSet::new();

    for (start, max_af) in max_af_by_position {
        if max_af < af_threshold {
            continue;
        }

        for offset in -radius..=radius {
            let candidate = start.saturating_add(offset);
            if source_positions.contains(&candidate) {
                warm_positions.insert(candidate);
            }
        }
    }

    warm_positions
}

#[cfg(test)]
mod v2_codec_tests {
    use super::*;
    use crate::transcript_consequence::CompactPrediction;

    fn cp(aa: u8, pred: u8, score: f32) -> CompactPrediction {
        CompactPrediction {
            position: 7,
            amino_acid: aa,
            prediction: pred,
            score,
        }
    }

    #[test]
    fn v2_sift_roundtrip_bit_exact() {
        let entries = vec![
            cp(0, 1, 0.0),
            cp(3, 0, 0.02),
            cp(19, 1, 1.0),
            cp(7, 0, 0.55),
        ];
        let bytes = serialize_position_entries_v2(&entries, SIFT_CODEC).unwrap();
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
        assert_eq!(bytes.len(), entries.len() * 4);
        let back = deserialize_position_entries_v2(7, &bytes, POLY_CODEC).unwrap();
        for (a, b) in entries.iter().zip(back.iter()) {
            assert_eq!(a.score.to_bits(), b.score.to_bits());
        }
    }

    #[test]
    fn v2_empty_blob() {
        let bytes = serialize_position_entries_v2(&[], SIFT_CODEC).unwrap();
        assert!(bytes.is_empty());
        assert!(
            deserialize_position_entries_v2(1, &bytes, SIFT_CODEC)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn v2_rejects_off_grid_score() {
        let entries = vec![cp(0, 1, 0.005)];
        assert!(serialize_position_entries_v2(&entries, SIFT_CODEC).is_err());
    }

    #[test]
    fn version_from_metadata() {
        let mut m = HashMap::new();
        assert_eq!(
            sift_blob_version_from_metadata(&m),
            SiftBlobVersion::V1Interleaved
        );
        m.insert(SIFT_BLOB_VERSION_KEY.to_string(), "2".to_string());
        assert_eq!(
            sift_blob_version_from_metadata(&m),
            SiftBlobVersion::V2DivIndex
        );
    }

    #[test]
    fn versioned_dispatch_v2() {
        let sift = serialize_position_entries_v2(&[cp(1, 0, 0.04)], SIFT_CODEC).unwrap();
        let poly = serialize_position_entries_v2(&[cp(1, 2, 0.7)], POLY_CODEC).unwrap();
        let out = deserialize_position_predictions_versioned(
            9,
            &sift,
            &poly,
            SiftBlobVersion::V2DivIndex,
        )
        .unwrap();
        assert_eq!(out.sift[0].score.to_bits(), 0.04f32.to_bits());
        assert_eq!(out.polyphen[0].score.to_bits(), 0.7f32.to_bits());
        assert_eq!(out.sift[0].position, 9);
    }
}
