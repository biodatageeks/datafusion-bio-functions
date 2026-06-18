//! Fjall-backed SIFT/PolyPhen prediction store.
//!
//! Key: transcript_id (UTF-8 bytes)
//! Value: serialized CompactPrediction entries (sift + polyphen)
//!
//! Binary format per value:
//!   [4B sift_count LE] [4B polyphen_count LE]
//!   [sift_count × 10B CompactPrediction]
//!   [polyphen_count × 10B CompactPrediction]

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, LazyLock, Mutex, Weak};

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
    inner: Arc<SiftKvStoreInner>,
}

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

struct SiftKvStoreInner {
    sift_ks: Keyspace,
}

static SHARED_SIFT_STORES: LazyLock<Mutex<HashMap<PathBuf, Weak<SiftKvStoreInner>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

fn canonical_store_path(path: &Path) -> PathBuf {
    std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf())
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
        let path = canonical_store_path(path);
        let mut stores = SHARED_SIFT_STORES.lock().map_err(|e| {
            DataFusionError::Execution(format!("shared fjall sift registry lock poisoned: {e}"))
        })?;
        stores.retain(|_, weak| weak.strong_count() > 0);
        if let Some(inner) = stores.get(&path).and_then(Weak::upgrade) {
            return Ok(Some(Self { inner }));
        }

        let db = Database::builder(&path)
            .cache_size(64 * 1024 * 1024)
            .worker_threads(1)
            .open()
            .map_err(fjall_err)?;
        let store = Self::open(&db)?;
        if let Some(store) = &store {
            stores.insert(path, Arc::downgrade(&store.inner));
        }
        Ok(store)
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
            Ok(Some(Self {
                inner: Arc::new(SiftKvStoreInner { sift_ks: ks }),
            }))
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
        Ok(Self {
            inner: Arc::new(SiftKvStoreInner { sift_ks }),
        })
    }

    /// Access the underlying keyspace (e.g. for compaction).
    pub fn keyspace(&self) -> &Keyspace {
        &self.inner.sift_ks
    }

    /// Store predictions for a transcript.
    pub fn put(&self, transcript_id: &str, preds: &CachedPredictions) -> Result<()> {
        let value = serialize_predictions(preds);
        self.inner
            .sift_ks
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
        let mut ingestion = store.inner.sift_ks.start_ingestion().map_err(fjall_err)?;
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
            .inner
            .sift_ks
            .get(transcript_id.as_bytes())
            .map_err(fjall_err)?
        else {
            return Ok(None);
        };
        deserialize_predictions(&raw).map(Some)
    }
}

impl SiftPredictionStore for SiftKvStore {
    fn get_many(&self, transcript_ids: &[String]) -> Result<HashMap<String, CachedPredictions>> {
        let mut out = HashMap::with_capacity(transcript_ids.len());
        for transcript_id in transcript_ids {
            if out.contains_key(transcript_id) {
                continue;
            }
            if let Some(predictions) = self.get(transcript_id)? {
                out.insert(transcript_id.clone(), predictions);
            }
        }
        Ok(out)
    }
}

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
    fn position_entries_roundtrip_without_position_field() {
        let sift = vec![
            CompactPrediction {
                position: 42,
                amino_acid: 0,
                prediction: 0,
                score: 0.30,
            },
            CompactPrediction {
                position: 42,
                amino_acid: 2,
                prediction: 1,
                score: 0.01,
            },
        ];
        let poly = vec![CompactPrediction {
            position: 42,
            amino_acid: 0,
            prediction: 4,
            score: 0.05,
        }];

        let sift_bytes = serialize_position_entries(&sift);
        let poly_bytes = serialize_position_entries(&poly);
        // 6 bytes per entry, position not stored.
        assert_eq!(sift_bytes.len(), sift.len() * 6);
        assert_eq!(poly_bytes.len(), poly.len() * 6);

        let decoded = deserialize_position_predictions(42, &sift_bytes, &poly_bytes).unwrap();
        assert_eq!(decoded.sift.len(), 2);
        assert_eq!(decoded.polyphen.len(), 1);
        for p in decoded.sift.iter().chain(decoded.polyphen.iter()) {
            assert_eq!(p.position, 42);
        }
        assert_eq!(decoded.sift[1].amino_acid, 2);
        assert_eq!(decoded.sift[1].prediction, 1);
        assert!((decoded.sift[1].score - 0.01).abs() < f32::EPSILON);
        assert_eq!(decoded.polyphen[0].prediction, 4);

        // Empty payload (a position with no predictions of that kind) is valid.
        let empty = deserialize_position_predictions(7, &[], &[]).unwrap();
        assert!(empty.sift.is_empty() && empty.polyphen.is_empty());

        // Malformed length is rejected.
        assert!(deserialize_position_entries(1, &[0u8; 5]).is_err());
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
    fn test_open_path_reuses_live_store() {
        let dir = tempfile::tempdir().unwrap();

        {
            let db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();
            let store = SiftKvStore::create(&db).unwrap();
            store.put("ENST00000111111", &make_predictions()).unwrap();
            db.persist(fjall::PersistMode::SyncAll).unwrap();
        }

        let first = SiftKvStore::open_path(dir.path()).unwrap().unwrap();
        let second = SiftKvStore::open_path(dir.path()).unwrap().unwrap();

        assert!(Arc::ptr_eq(&first.inner, &second.inner));
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
    fn test_sift_prediction_store_get_many_returns_found_and_skips_missing() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let store = SiftKvStore::create(&db).unwrap();
        store.put("ENST00000111111", &make_predictions()).unwrap();
        store.put("ENST00000222222", &make_predictions()).unwrap();

        let ids = vec![
            "ENST00000111111".to_string(),
            "missing".to_string(),
            "ENST00000111111".to_string(),
            "ENST00000222222".to_string(),
        ];
        let found = SiftPredictionStore::get_many(&store, &ids).unwrap();

        assert_eq!(found.len(), 2);
        assert!(found.contains_key("ENST00000111111"));
        assert!(found.contains_key("ENST00000222222"));
        assert!(!found.contains_key("missing"));
    }

    /// Verify that `SiftKvStore::create()` only creates the "sift" keyspace
    /// (plus fjall's internal default). Extra keyspaces degrade performance.
    #[test]
    fn test_create_produces_minimal_keyspaces() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let _store = SiftKvStore::create(&db).unwrap();
        db.persist(fjall::PersistMode::SyncAll).unwrap();

        // Only the "sift" keyspace should exist (no meta/data/translations/exons).
        assert!(db.keyspace_exists(SIFT_KEYSPACE));
        for name in ["translations", "exons"] {
            assert!(
                !db.keyspace_exists(name),
                "SiftKvStore::create() created unexpected keyspace '{name}'"
            );
        }

        // fjall default + sift = 2 dirs on disk.
        let ks_dir = dir.path().join("keyspaces");
        let count = std::fs::read_dir(&ks_dir)
            .unwrap()
            .filter(|e| {
                e.as_ref()
                    .ok()
                    .and_then(|e| e.file_name().to_string_lossy().parse::<u32>().ok())
                    .is_some()
            })
            .count();
        assert_eq!(
            count, 2,
            "SiftKvStore::create() must produce exactly 2 keyspace dirs (fjall default + sift), found {count}"
        );
    }

    /// Verify that `ingest_sorted()` only creates the "sift" keyspace.
    #[test]
    fn test_ingest_sorted_produces_minimal_keyspaces() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let entries = vec![
            ("ENST00000000001".to_string(), make_predictions()),
            ("ENST00000000002".to_string(), make_predictions()),
        ];
        let _store = SiftKvStore::ingest_sorted(&db, entries.into_iter()).unwrap();
        db.persist(fjall::PersistMode::SyncAll).unwrap();

        assert!(db.keyspace_exists(SIFT_KEYSPACE));
        for name in ["translations", "exons"] {
            assert!(
                !db.keyspace_exists(name),
                "ingest_sorted() created unexpected keyspace '{name}'"
            );
        }

        let ks_dir = dir.path().join("keyspaces");
        let count = std::fs::read_dir(&ks_dir)
            .unwrap()
            .filter(|e| {
                e.as_ref()
                    .ok()
                    .and_then(|e| e.file_name().to_string_lossy().parse::<u32>().ok())
                    .is_some()
            })
            .count();
        assert_eq!(
            count, 2,
            "ingest_sorted() must produce exactly 2 keyspace dirs (fjall default + sift), found {count}"
        );
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
