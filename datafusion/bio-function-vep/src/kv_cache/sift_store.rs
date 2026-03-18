//! Fjall-backed SIFT/PolyPhen prediction store.
//!
//! Key: transcript_id (UTF-8 bytes)
//! Value: serialized CompactPrediction entries (sift + polyphen)
//!
//! Binary format per value:
//!   [4B sift_count LE] [4B polyphen_count LE]
//!   [sift_count × 10B CompactPrediction]
//!   [polyphen_count × 10B CompactPrediction]

use std::path::Path;
use std::sync::Arc;

use datafusion::common::{DataFusionError, Result};
use fjall::{Database, Keyspace, KeyspaceCreateOptions};

use crate::transcript_consequence::{CachedPredictions, CompactPrediction};

const SIFT_KEYSPACE: &str = "sift";

fn fjall_err(e: fjall::Error) -> DataFusionError {
    DataFusionError::External(Box::new(e))
}

/// Store for SIFT/PolyPhen predictions keyed by transcript_id.
pub struct SiftKvStore {
    sift_ks: Keyspace,
}

impl SiftKvStore {
    /// Open sift keyspace from an existing fjall database.
    pub fn open(db: &Database) -> Result<Option<Self>> {
        // Only return Some if the keyspace exists (don't create on read).
        match db.keyspace(SIFT_KEYSPACE, KeyspaceCreateOptions::default) {
            Ok(ks) => {
                // Check if it has any data
                if ks.is_empty().unwrap_or(true) {
                    Ok(None)
                } else {
                    Ok(Some(Self { sift_ks: ks }))
                }
            }
            Err(_) => Ok(None),
        }
    }

    /// Create sift keyspace in a fjall database for bulk loading.
    pub fn create(db: &Database) -> Result<Self> {
        let sift_ks = db
            .keyspace(SIFT_KEYSPACE, || {
                KeyspaceCreateOptions::default()
                    .manual_journal_persist(true)
                    .compaction_strategy(Arc::new(
                        fjall::compaction::Leveled::default().with_l0_threshold(16),
                    ))
                    .data_block_compression_policy(fjall::config::CompressionPolicy::disabled())
            })
            .map_err(fjall_err)?;
        Ok(Self { sift_ks })
    }

    /// Store predictions for a transcript.
    pub fn put(&self, transcript_id: &str, preds: &CachedPredictions) -> Result<()> {
        let value = serialize_predictions(preds);
        self.sift_ks
            .insert(transcript_id.as_bytes(), value)
            .map_err(fjall_err)?;
        Ok(())
    }

    /// Retrieve predictions for a transcript. Returns None on miss.
    pub fn get(&self, transcript_id: &str) -> Result<Option<CachedPredictions>> {
        let Some(raw) = self
            .sift_ks
            .get(transcript_id.as_bytes())
            .map_err(fjall_err)?
        else {
            return Ok(None);
        };
        deserialize_predictions(&raw).map(Some)
    }
}

fn serialize_predictions(preds: &CachedPredictions) -> Vec<u8> {
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

fn deserialize_predictions(data: &[u8]) -> Result<CachedPredictions> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript_consequence::{CachedPredictions, CompactPrediction};

    fn sample_predictions() -> CachedPredictions {
        CachedPredictions {
            sift: vec![
                CompactPrediction {
                    position: 42,
                    amino_acid: b'A',
                    prediction: 1,
                    score: 0.05,
                },
                CompactPrediction {
                    position: 100,
                    amino_acid: b'G',
                    prediction: 0,
                    score: 0.95,
                },
            ],
            polyphen: vec![CompactPrediction {
                position: 42,
                amino_acid: b'A',
                prediction: 2,
                score: 0.99,
            }],
        }
    }

    #[test]
    fn test_sift_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let store = SiftKvStore::create(&db).unwrap();

        let preds = sample_predictions();
        store.put("ENST00000123456", &preds).unwrap();

        let loaded = store.get("ENST00000123456").unwrap().unwrap();
        assert_eq!(loaded.sift.len(), 2);
        assert_eq!(loaded.polyphen.len(), 1);
        assert_eq!(loaded.sift[0].position, 42);
        assert_eq!(loaded.sift[0].amino_acid, b'A');
        assert!((loaded.sift[0].score - 0.05).abs() < f32::EPSILON);
        assert_eq!(loaded.polyphen[0].prediction, 2);
    }

    #[test]
    fn test_sift_missing_key_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let store = SiftKvStore::create(&db).unwrap();
        assert!(store.get("NONEXISTENT").unwrap().is_none());
    }

    #[test]
    fn test_deserialize_truncated_buffer_returns_error() {
        // Buffer too short for header (< 8 bytes).
        let short = vec![0u8; 4];
        let result = deserialize_predictions(&short);
        assert!(result.is_err());

        // Header claims data but buffer is truncated.
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u32.to_le_bytes()); // sift_count = 1
        buf.extend_from_slice(&0u32.to_le_bytes()); // polyphen_count = 0
        // Should need 8 + 10 = 18 bytes, but we only have 8.
        let result = deserialize_predictions(&buf);
        assert!(result.is_err());
    }
}
