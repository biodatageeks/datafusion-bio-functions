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
#[derive(Clone)]
pub struct SiftKvStore {
    sift_ks: Keyspace,
}

impl SiftKvStore {
    /// Open a standalone SIFT fjall database at the given path.
    ///
    /// Used when SIFT predictions are stored in a separate `translation_sift.fjall/`
    /// directory (matching the parquet naming convention).
    pub fn open_path(path: impl AsRef<std::path::Path>) -> Result<Option<Self>> {
        let path = path.as_ref();
        if !path.exists() {
            return Ok(None);
        }
        let db = Database::builder(path)
            .cache_size(64 * 1024 * 1024)
            .open()
            .map_err(fjall_err)?;
        Self::open(&db)
    }

    /// Open sift keyspace from an existing fjall database.
    /// Returns `None` if the keyspace doesn't exist or is empty.
    pub fn open(db: &Database) -> Result<Option<Self>> {
        if !db.keyspace_exists(SIFT_KEYSPACE) {
            return Ok(None);
        }
        let ks = db
            .keyspace(SIFT_KEYSPACE, KeyspaceCreateOptions::default)
            .map_err(fjall_err)?;
        if ks.is_empty().unwrap_or(true) {
            Ok(None)
        } else {
            Ok(Some(Self { sift_ks: ks }))
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

    /// Bulk-load from a sorted-by-transcript_id iterator using fjall ingestion.
    ///
    /// Input **must** be sorted in ascending `transcript_id` order (returns `Err` otherwise).
    /// Uses `start_ingestion()` for maximum bulk load speed.
    pub fn ingest_sorted(
        db: &Database,
        sorted_iter: impl Iterator<Item = (String, CachedPredictions)>,
    ) -> Result<Self> {
        let store = Self::create(db)?;
        let mut ingestion = store.sift_ks.start_ingestion().map_err(fjall_err)?;
        for (transcript_id, preds) in sorted_iter {
            let value = serialize_predictions(&preds);
            ingestion
                .write(transcript_id.as_bytes(), value)
                .map_err(fjall_err)?;
        }
        ingestion.finish().map_err(fjall_err)?;
        Ok(store)
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

    fn make_predictions() -> CachedPredictions {
        CachedPredictions {
            sift: vec![
                CompactPrediction {
                    position: 10,
                    amino_acid: 1,
                    prediction: 0,
                    score: 0.05,
                },
                CompactPrediction {
                    position: 20,
                    amino_acid: 2,
                    prediction: 1,
                    score: 0.95,
                },
            ],
            polyphen: vec![CompactPrediction {
                position: 10,
                amino_acid: 1,
                prediction: 2,
                score: 0.88,
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
        let preds = make_predictions();
        store.put("ENST00000123456", &preds).unwrap();

        let loaded = store.get("ENST00000123456").unwrap().unwrap();
        assert_eq!(loaded.sift.len(), 2);
        assert_eq!(loaded.polyphen.len(), 1);
        assert_eq!(loaded.sift[0].position, 10);
        assert_eq!(loaded.sift[1].position, 20);
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
        store.put("ENST00000123456", &make_predictions()).unwrap();

        assert!(store.get("MISSING").unwrap().is_none());
    }

    #[test]
    fn test_deserialize_truncated_header() {
        // Only 4 bytes — too short for the 8-byte header.
        let result = deserialize_predictions(&[0, 0, 0, 0]);
        assert!(result.is_err());
    }

    #[test]
    fn test_deserialize_truncated_body() {
        // Header claims 1 sift + 1 polyphen (needs 8 + 20 = 28 bytes)
        // but only 8 bytes provided.
        let mut data = vec![0u8; 8];
        data[0..4].copy_from_slice(&1u32.to_le_bytes()); // sift_count = 1
        data[4..8].copy_from_slice(&1u32.to_le_bytes()); // polyphen_count = 1
        let result = deserialize_predictions(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_open_empty_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        assert!(SiftKvStore::open(&db).unwrap().is_none());
    }

    #[test]
    fn test_open_path_nonexistent() {
        let result = SiftKvStore::open_path("/nonexistent/path/to/sift.fjall").unwrap();
        assert!(result.is_none(), "Non-existent path should return Ok(None)");
    }

    #[test]
    fn test_open_path_roundtrip() {
        let dir = tempfile::tempdir().unwrap();

        // Write phase: create a standalone SIFT DB and insert predictions.
        {
            let db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();
            let store = SiftKvStore::create(&db).unwrap();
            store.put("ENST00000111111", &make_predictions()).unwrap();
            db.persist(fjall::PersistMode::SyncAll).unwrap();
        }

        // Read phase: reopen via open_path and verify data.
        let store = SiftKvStore::open_path(dir.path())
            .unwrap()
            .expect("open_path should return Some for a valid sift DB");
        let preds = store
            .get("ENST00000111111")
            .unwrap()
            .expect("predictions should be present");
        assert_eq!(preds.sift.len(), 2);
        assert_eq!(preds.polyphen.len(), 1);
        assert_eq!(preds.sift[0].position, 10);
    }

    #[test]
    fn test_clone_shares_data() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let store = SiftKvStore::create(&db).unwrap();
        store.put("ENST00000222222", &make_predictions()).unwrap();

        let cloned = store.clone();

        // Both original and clone should see the same data.
        let from_original = store.get("ENST00000222222").unwrap().unwrap();
        let from_clone = cloned.get("ENST00000222222").unwrap().unwrap();
        assert_eq!(from_original.sift.len(), from_clone.sift.len());
        assert_eq!(from_original.polyphen.len(), from_clone.polyphen.len());

        // Write through clone should be visible from original.
        let extra = CachedPredictions {
            sift: vec![CompactPrediction {
                position: 99,
                amino_acid: 5,
                prediction: 1,
                score: 0.42,
            }],
            polyphen: vec![],
        };
        cloned.put("ENST00000333333", &extra).unwrap();
        assert!(store.get("ENST00000333333").unwrap().is_some());
    }

    #[test]
    fn test_open_keyspace_not_exists() {
        let dir = tempfile::tempdir().unwrap();

        // Create a DB without the sift keyspace (just an empty fjall DB).
        {
            let _db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();
            // Don't create any keyspaces — just open and close.
        }

        // Reopen and try to open sift keyspace — should return None.
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();
        let result = SiftKvStore::open(&db).unwrap();
        assert!(
            result.is_none(),
            "Opening sift keyspace on DB without it should return None"
        );
    }
}
