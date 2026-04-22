//! Fjall-backed stores for translation and exon context tables.
//!
//! These stores live as additional keyspaces within the same fjall `Database`
//! used by `VepKvStore`, keyed by `transcript_id` (UTF-8 bytes).
//!
//! ## Translation entry binary format
//!
//! ```text
//! [1B flags]
//!   bit 0: has cds_len
//!   bit 1: has protein_len
//!   bit 2: has translation_seq                   (BAM-edited / vef_cache.peptide)
//!   bit 3: has cds_sequence                      (BAM-edited / vef_cache.translateable_seq)
//!   bit 4: has stable_id
//!   bit 5: has version
//!   bit 6: has translation_seq_canonical         (upstream d26e370, translation.primary_seq)
//!   bit 7: has cds_sequence_canonical            (upstream d26e370, transcript.translateable_seq)
//! [optional 8B cds_len u64 LE]
//! [optional 8B protein_len u64 LE]
//! [optional 4B translation_seq_len u32 LE + bytes]
//! [optional 4B cds_sequence_len u32 LE + bytes]
//! [optional 4B stable_id_len u32 LE + bytes]
//! [optional 4B version i32 LE]
//! [optional 4B translation_seq_canonical_len u32 LE + bytes]
//! [optional 4B cds_sequence_canonical_len u32 LE + bytes]
//! [4B protein_features_count u32 LE]
//! For each protein feature:
//!   [1B flags: bit 0 = has analysis, bit 1 = has hseqname]
//!   [optional 4B analysis_len u32 LE + bytes]
//!   [optional 4B hseqname_len u32 LE + bytes]
//!   [8B start i64 LE]
//!   [8B end i64 LE]
//! ```
//!
//! The canonical bits are driven purely off `Option::is_some`: writers emit
//! the payload whenever the field is populated, readers deliver `Some(_)`
//! iff the flag is set and `None` otherwise. There is no legacy fallback
//! path — fjall databases written by earlier versions must be regenerated.
//!
//! ## Exon entry binary format (per transcript_id)
//!
//! ```text
//! [4B exon_count u32 LE]
//! For each exon (sorted by exon_number):
//!   [4B exon_number i32 LE]
//!   [8B start i64 LE]
//!   [8B end i64 LE]
//! ```

use std::sync::Arc;

use datafusion::common::{DataFusionError, Result};
use fjall::{Database, Keyspace, KeyspaceCreateOptions};

use crate::transcript_consequence::{ExonFeature, ProteinDomainFeature, TranslationFeature};

const TRANSLATIONS_KEYSPACE: &str = "translations";
const EXONS_KEYSPACE: &str = "exons";

fn fjall_err(e: fjall::Error) -> DataFusionError {
    DataFusionError::External(Box::new(e))
}

// ---------------------------------------------------------------------------
// TranslationKvStore
// ---------------------------------------------------------------------------

/// Fjall-backed store for translation features keyed by transcript_id.
pub struct TranslationKvStore {
    ks: Keyspace,
}

impl TranslationKvStore {
    /// Open the translations keyspace from an existing fjall database.
    /// Returns `None` if the keyspace doesn't exist or is empty.
    pub fn open(db: &Database) -> Result<Option<Self>> {
        if !db.keyspace_exists(TRANSLATIONS_KEYSPACE) {
            return Ok(None);
        }
        let ks = db
            .keyspace(TRANSLATIONS_KEYSPACE, KeyspaceCreateOptions::default)
            .map_err(fjall_err)?;
        if ks.is_empty().unwrap_or(true) {
            Ok(None)
        } else {
            Ok(Some(Self { ks }))
        }
    }

    /// Create translations keyspace for bulk loading.
    pub fn create(db: &Database) -> Result<Self> {
        let ks = db
            .keyspace(TRANSLATIONS_KEYSPACE, || {
                KeyspaceCreateOptions::default()
                    .manual_journal_persist(true)
                    .compaction_strategy(Arc::new(
                        fjall::compaction::Leveled::default().with_l0_threshold(16),
                    ))
                    .data_block_compression_policy(fjall::config::CompressionPolicy::disabled())
            })
            .map_err(fjall_err)?;
        Ok(Self { ks })
    }

    /// Store a translation feature.
    pub fn put(&self, feature: &TranslationFeature) -> Result<()> {
        let value = serialize_translation(feature);
        self.ks
            .insert(feature.transcript_id.as_bytes(), value)
            .map_err(fjall_err)?;
        Ok(())
    }

    /// Retrieve a translation feature by transcript_id. Returns `None` on miss.
    pub fn get(&self, transcript_id: &str) -> Result<Option<TranslationFeature>> {
        let Some(raw) = self.ks.get(transcript_id.as_bytes()).map_err(fjall_err)? else {
            return Ok(None);
        };
        deserialize_translation(transcript_id, &raw).map(Some)
    }

    /// Retrieve translations for a set of transcript_ids.
    /// Returns only those that exist in the store.
    pub fn get_many(&self, transcript_ids: &[&str]) -> Result<Vec<TranslationFeature>> {
        let mut out = Vec::with_capacity(transcript_ids.len());
        for &tid in transcript_ids {
            if let Some(t) = self.get(tid)? {
                out.push(t);
            }
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// ExonKvStore
// ---------------------------------------------------------------------------

/// Fjall-backed store for exon features keyed by transcript_id.
/// Each key stores all exons for that transcript as a single entry.
pub struct ExonKvStore {
    ks: Keyspace,
}

impl ExonKvStore {
    /// Open the exons keyspace from an existing fjall database.
    /// Returns `None` if the keyspace doesn't exist or is empty.
    pub fn open(db: &Database) -> Result<Option<Self>> {
        if !db.keyspace_exists(EXONS_KEYSPACE) {
            return Ok(None);
        }
        let ks = db
            .keyspace(EXONS_KEYSPACE, KeyspaceCreateOptions::default)
            .map_err(fjall_err)?;
        if ks.is_empty().unwrap_or(true) {
            Ok(None)
        } else {
            Ok(Some(Self { ks }))
        }
    }

    /// Create exons keyspace for bulk loading.
    pub fn create(db: &Database) -> Result<Self> {
        let ks = db
            .keyspace(EXONS_KEYSPACE, || {
                KeyspaceCreateOptions::default()
                    .manual_journal_persist(true)
                    .compaction_strategy(Arc::new(
                        fjall::compaction::Leveled::default().with_l0_threshold(16),
                    ))
                    .data_block_compression_policy(fjall::config::CompressionPolicy::disabled())
            })
            .map_err(fjall_err)?;
        Ok(Self { ks })
    }

    /// Store all exons for a transcript. Exons are sorted by exon_number before storing.
    pub fn put(&self, transcript_id: &str, exons: &mut [ExonFeature]) -> Result<()> {
        exons.sort_by_key(|e| e.exon_number);
        let value = serialize_exons(exons);
        self.ks
            .insert(transcript_id.as_bytes(), value)
            .map_err(fjall_err)?;
        Ok(())
    }

    /// Retrieve exons for a transcript_id. Returns empty vec on miss.
    pub fn get(&self, transcript_id: &str) -> Result<Vec<ExonFeature>> {
        let Some(raw) = self.ks.get(transcript_id.as_bytes()).map_err(fjall_err)? else {
            return Ok(Vec::new());
        };
        deserialize_exons(transcript_id, &raw)
    }

    /// Retrieve exons for a set of transcript_ids.
    pub fn get_many(&self, transcript_ids: &[&str]) -> Result<Vec<ExonFeature>> {
        let mut out = Vec::new();
        for &tid in transcript_ids {
            out.extend(self.get(tid)?);
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// Translation serialization
// ---------------------------------------------------------------------------

const FLAG_HAS_CDS_LEN: u8 = 1 << 0;
const FLAG_HAS_PROTEIN_LEN: u8 = 1 << 1;
const FLAG_HAS_TRANSLATION_SEQ: u8 = 1 << 2;
const FLAG_HAS_CDS_SEQUENCE: u8 = 1 << 3;
const FLAG_HAS_STABLE_ID: u8 = 1 << 4;
const FLAG_HAS_VERSION: u8 = 1 << 5;
const FLAG_HAS_TRANSLATION_SEQ_CANONICAL: u8 = 1 << 6;
const FLAG_HAS_CDS_SEQUENCE_CANONICAL: u8 = 1 << 7;

fn serialize_translation(t: &TranslationFeature) -> Vec<u8> {
    let mut buf = Vec::with_capacity(128);

    let mut flags: u8 = 0;
    if t.cds_len.is_some() {
        flags |= FLAG_HAS_CDS_LEN;
    }
    if t.protein_len.is_some() {
        flags |= FLAG_HAS_PROTEIN_LEN;
    }
    if t.translation_seq.is_some() {
        flags |= FLAG_HAS_TRANSLATION_SEQ;
    }
    if t.cds_sequence.is_some() {
        flags |= FLAG_HAS_CDS_SEQUENCE;
    }
    if t.stable_id.is_some() {
        flags |= FLAG_HAS_STABLE_ID;
    }
    if t.version.is_some() {
        flags |= FLAG_HAS_VERSION;
    }
    if t.translation_seq_canonical.is_some() {
        flags |= FLAG_HAS_TRANSLATION_SEQ_CANONICAL;
    }
    if t.cds_sequence_canonical.is_some() {
        flags |= FLAG_HAS_CDS_SEQUENCE_CANONICAL;
    }
    buf.push(flags);

    if let Some(v) = t.cds_len {
        buf.extend_from_slice(&(v as u64).to_le_bytes());
    }
    if let Some(v) = t.protein_len {
        buf.extend_from_slice(&(v as u64).to_le_bytes());
    }
    if let Some(ref s) = t.translation_seq {
        write_str(&mut buf, s);
    }
    if let Some(ref s) = t.cds_sequence {
        write_str(&mut buf, s);
    }
    if let Some(ref s) = t.stable_id {
        write_str(&mut buf, s);
    }
    if let Some(v) = t.version {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    if let Some(ref s) = t.translation_seq_canonical {
        write_str(&mut buf, s);
    }
    if let Some(ref s) = t.cds_sequence_canonical {
        write_str(&mut buf, s);
    }

    // Protein features.
    buf.extend_from_slice(&(t.protein_features.len() as u32).to_le_bytes());
    for pf in &t.protein_features {
        let mut pf_flags: u8 = 0;
        if pf.analysis.is_some() {
            pf_flags |= 1;
        }
        if pf.hseqname.is_some() {
            pf_flags |= 2;
        }
        buf.push(pf_flags);
        if let Some(ref s) = pf.analysis {
            write_str(&mut buf, s);
        }
        if let Some(ref s) = pf.hseqname {
            write_str(&mut buf, s);
        }
        buf.extend_from_slice(&pf.start.to_le_bytes());
        buf.extend_from_slice(&pf.end.to_le_bytes());
    }

    buf
}

fn deserialize_translation(transcript_id: &str, data: &[u8]) -> Result<TranslationFeature> {
    if data.is_empty() {
        return Err(DataFusionError::Execution(
            "translation entry is empty".into(),
        ));
    }
    let mut off = 0;
    let flags = data[off];
    off += 1;

    let cds_len = if flags & FLAG_HAS_CDS_LEN != 0 {
        let v = read_u64(data, &mut off)?;
        Some(v as usize)
    } else {
        None
    };

    let protein_len = if flags & FLAG_HAS_PROTEIN_LEN != 0 {
        let v = read_u64(data, &mut off)?;
        Some(v as usize)
    } else {
        None
    };

    let translation_seq = if flags & FLAG_HAS_TRANSLATION_SEQ != 0 {
        Some(read_string(data, &mut off)?)
    } else {
        None
    };

    let cds_sequence = if flags & FLAG_HAS_CDS_SEQUENCE != 0 {
        Some(read_string(data, &mut off)?)
    } else {
        None
    };

    let stable_id = if flags & FLAG_HAS_STABLE_ID != 0 {
        Some(read_string(data, &mut off)?)
    } else {
        None
    };

    let version = if flags & FLAG_HAS_VERSION != 0 {
        let v = read_i32(data, &mut off)?;
        Some(v)
    } else {
        None
    };

    let translation_seq_canonical = if flags & FLAG_HAS_TRANSLATION_SEQ_CANONICAL != 0 {
        Some(read_string(data, &mut off)?)
    } else {
        None
    };

    let cds_sequence_canonical = if flags & FLAG_HAS_CDS_SEQUENCE_CANONICAL != 0 {
        Some(read_string(data, &mut off)?)
    } else {
        None
    };

    let pf_count = read_u32(data, &mut off)? as usize;
    let mut protein_features = Vec::with_capacity(pf_count);
    for _ in 0..pf_count {
        check_len(data, off, 1)?;
        let pf_flags = data[off];
        off += 1;
        let analysis = if pf_flags & 1 != 0 {
            Some(read_string(data, &mut off)?)
        } else {
            None
        };
        let hseqname = if pf_flags & 2 != 0 {
            Some(read_string(data, &mut off)?)
        } else {
            None
        };
        let start = read_i64(data, &mut off)?;
        let end = read_i64(data, &mut off)?;
        protein_features.push(ProteinDomainFeature {
            analysis,
            hseqname,
            start,
            end,
        });
    }

    Ok(TranslationFeature {
        transcript_id: transcript_id.to_string(),
        cds_len,
        protein_len,
        translation_seq,
        cds_sequence,
        translation_seq_canonical,
        cds_sequence_canonical,
        stable_id,
        version,
        protein_features,
    })
}

// ---------------------------------------------------------------------------
// Exon serialization
// ---------------------------------------------------------------------------

fn serialize_exons(exons: &[ExonFeature]) -> Vec<u8> {
    // 4 + count * 20
    let mut buf = Vec::with_capacity(4 + exons.len() * 20);
    buf.extend_from_slice(&(exons.len() as u32).to_le_bytes());
    for e in exons {
        buf.extend_from_slice(&e.exon_number.to_le_bytes());
        buf.extend_from_slice(&e.start.to_le_bytes());
        buf.extend_from_slice(&e.end.to_le_bytes());
    }
    buf
}

fn deserialize_exons(transcript_id: &str, data: &[u8]) -> Result<Vec<ExonFeature>> {
    if data.len() < 4 {
        return Err(DataFusionError::Execution("exon entry too short".into()));
    }
    let mut off = 0;
    let count = read_u32(data, &mut off)? as usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let exon_number = read_i32(data, &mut off)?;
        let start = read_i64(data, &mut off)?;
        let end = read_i64(data, &mut off)?;
        out.push(ExonFeature {
            transcript_id: transcript_id.to_string(),
            exon_number,
            start,
            end,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Primitive read/write helpers
// ---------------------------------------------------------------------------

fn write_str(buf: &mut Vec<u8>, s: &str) {
    buf.extend_from_slice(&(s.len() as u32).to_le_bytes());
    buf.extend_from_slice(s.as_bytes());
}

fn check_len(data: &[u8], off: usize, need: usize) -> Result<()> {
    if off + need > data.len() {
        Err(DataFusionError::Execution(format!(
            "context entry truncated at offset {off}, need {need} more bytes, have {}",
            data.len() - off
        )))
    } else {
        Ok(())
    }
}

fn read_u32(data: &[u8], off: &mut usize) -> Result<u32> {
    check_len(data, *off, 4)?;
    let v = u32::from_le_bytes(data[*off..*off + 4].try_into().unwrap());
    *off += 4;
    Ok(v)
}

fn read_i32(data: &[u8], off: &mut usize) -> Result<i32> {
    check_len(data, *off, 4)?;
    let v = i32::from_le_bytes(data[*off..*off + 4].try_into().unwrap());
    *off += 4;
    Ok(v)
}

fn read_u64(data: &[u8], off: &mut usize) -> Result<u64> {
    check_len(data, *off, 8)?;
    let v = u64::from_le_bytes(data[*off..*off + 8].try_into().unwrap());
    *off += 8;
    Ok(v)
}

fn read_i64(data: &[u8], off: &mut usize) -> Result<i64> {
    check_len(data, *off, 8)?;
    let v = i64::from_le_bytes(data[*off..*off + 8].try_into().unwrap());
    *off += 8;
    Ok(v)
}

fn read_string(data: &[u8], off: &mut usize) -> Result<String> {
    let len = read_u32(data, off)? as usize;
    check_len(data, *off, len)?;
    let s = std::str::from_utf8(&data[*off..*off + len])
        .map_err(|e| DataFusionError::Execution(format!("invalid UTF-8 in context entry: {e}")))?;
    *off += len;
    Ok(s.to_string())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_translation() -> TranslationFeature {
        TranslationFeature {
            transcript_id: "ENST00000123456".to_string(),
            cds_len: Some(1200),
            protein_len: Some(400),
            translation_seq: Some("MSEQ...".to_string()),
            cds_sequence: Some("ATGCDS...".to_string()),
            translation_seq_canonical: Some("MSEQ...".to_string()),
            cds_sequence_canonical: Some("ATGCDS...".to_string()),
            stable_id: Some("ENSP00000123456".to_string()),
            version: Some(3),
            protein_features: vec![
                ProteinDomainFeature {
                    analysis: Some("Pfam".to_string()),
                    hseqname: Some("PF00001".to_string()),
                    start: 10,
                    end: 100,
                },
                ProteinDomainFeature {
                    analysis: None,
                    hseqname: Some("SSF12345".to_string()),
                    start: 50,
                    end: 200,
                },
            ],
        }
    }

    fn make_translation_minimal() -> TranslationFeature {
        TranslationFeature {
            transcript_id: "ENST00000999999".to_string(),
            cds_len: None,
            protein_len: None,
            translation_seq: None,
            cds_sequence: None,
            translation_seq_canonical: None,
            cds_sequence_canonical: None,
            stable_id: None,
            version: None,
            protein_features: Vec::new(),
        }
    }

    fn make_exons() -> Vec<ExonFeature> {
        vec![
            ExonFeature {
                transcript_id: "ENST00000123456".to_string(),
                exon_number: 3,
                start: 300,
                end: 400,
            },
            ExonFeature {
                transcript_id: "ENST00000123456".to_string(),
                exon_number: 1,
                start: 100,
                end: 200,
            },
            ExonFeature {
                transcript_id: "ENST00000123456".to_string(),
                exon_number: 2,
                start: 220,
                end: 280,
            },
        ]
    }

    #[test]
    fn test_translation_roundtrip_full() {
        let t = make_translation();
        let bytes = serialize_translation(&t);
        let restored = deserialize_translation(&t.transcript_id, &bytes).unwrap();
        assert_eq!(t, restored);
    }

    /// BAM-edited RefSeq transcripts: canonical ≠ edited. Both must survive
    /// the roundtrip with their original values intact (no mirror collapse).
    #[test]
    fn test_translation_roundtrip_distinct_canonical() {
        let t = TranslationFeature {
            transcript_id: "NM_002111.8".to_string(),
            cds_len: Some(9435),
            protein_len: Some(3144),
            translation_seq: Some("MATLEKLMKAFESLKSFQQQQQ".to_string()), // BAM-edited (5 Qs)
            cds_sequence: Some(
                "ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAG".to_string(),
            ),
            translation_seq_canonical: Some("MATLEKLMKAFESLKSFQQQ".to_string()), // canonical (3 Qs)
            cds_sequence_canonical: Some(
                "ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAG".to_string(),
            ),
            stable_id: Some("NP_002102.4".to_string()),
            version: Some(4),
            protein_features: Vec::new(),
        };
        let bytes = serialize_translation(&t);
        let restored = deserialize_translation(&t.transcript_id, &bytes).unwrap();
        assert_eq!(t, restored);
        // Specifically guarantee canonical and edited stayed distinct.
        assert_ne!(restored.translation_seq, restored.translation_seq_canonical);
        assert_ne!(restored.cds_sequence, restored.cds_sequence_canonical);
    }

    /// Absent canonical sequences must serialize with the flag bits cleared
    /// and deserialize back to `None`. Guards against a regression that would
    /// synthesize canonical values from the BAM-edited ones.
    #[test]
    fn test_translation_roundtrip_absent_canonical_stays_none() {
        let t = TranslationFeature {
            transcript_id: "ENSTNOCANON0001".to_string(),
            cds_len: Some(900),
            protein_len: Some(300),
            translation_seq: Some("MSEQ...".to_string()),
            cds_sequence: Some("ATGCDS...".to_string()),
            translation_seq_canonical: None,
            cds_sequence_canonical: None,
            stable_id: Some("ENSPNOCANON0001".to_string()),
            version: Some(1),
            protein_features: Vec::new(),
        };
        let bytes = serialize_translation(&t);

        let flags = bytes[0];
        assert_eq!(flags & FLAG_HAS_TRANSLATION_SEQ_CANONICAL, 0);
        assert_eq!(flags & FLAG_HAS_CDS_SEQUENCE_CANONICAL, 0);

        let restored = deserialize_translation(&t.transcript_id, &bytes).unwrap();
        assert_eq!(restored.translation_seq_canonical, None);
        assert_eq!(restored.cds_sequence_canonical, None);
        assert_eq!(t, restored);
    }

    #[test]
    fn test_translation_roundtrip_minimal() {
        let t = make_translation_minimal();
        let bytes = serialize_translation(&t);
        let restored = deserialize_translation(&t.transcript_id, &bytes).unwrap();
        assert_eq!(t, restored);
    }

    #[test]
    fn test_exon_roundtrip() {
        let mut exons = make_exons();
        // serialize_exons expects pre-sorted input; put() sorts for us.
        exons.sort_by_key(|e| e.exon_number);
        let bytes = serialize_exons(&exons);
        let restored = deserialize_exons("ENST00000123456", &bytes).unwrap();
        assert_eq!(exons, restored);
    }

    #[test]
    fn test_exon_empty() {
        let exons: Vec<ExonFeature> = Vec::new();
        let bytes = serialize_exons(&exons);
        let restored = deserialize_exons("ENST00000000000", &bytes).unwrap();
        assert!(restored.is_empty());
    }

    #[test]
    fn test_translation_kv_store() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let store = TranslationKvStore::create(&db).unwrap();

        let t1 = make_translation();
        let t2 = make_translation_minimal();
        store.put(&t1).unwrap();
        store.put(&t2).unwrap();

        let loaded = store.get("ENST00000123456").unwrap().unwrap();
        assert_eq!(loaded, t1);

        let loaded = store.get("ENST00000999999").unwrap().unwrap();
        assert_eq!(loaded, t2);

        assert!(store.get("MISSING").unwrap().is_none());

        let many = store
            .get_many(&["ENST00000123456", "MISSING", "ENST00000999999"])
            .unwrap();
        assert_eq!(many.len(), 2);
        assert_eq!(many[0].transcript_id, "ENST00000123456");
        assert_eq!(many[1].transcript_id, "ENST00000999999");
    }

    #[test]
    fn test_exon_kv_store() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let store = ExonKvStore::create(&db).unwrap();

        let mut exons = make_exons();
        store.put("ENST00000123456", &mut exons).unwrap();

        let loaded = store.get("ENST00000123456").unwrap();
        assert_eq!(loaded.len(), 3);
        // Should be sorted by exon_number.
        assert_eq!(loaded[0].exon_number, 1);
        assert_eq!(loaded[1].exon_number, 2);
        assert_eq!(loaded[2].exon_number, 3);

        assert!(store.get("MISSING").unwrap().is_empty());
    }

    #[test]
    fn test_open_empty_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        assert!(TranslationKvStore::open(&db).unwrap().is_none());
        assert!(ExonKvStore::open(&db).unwrap().is_none());
    }

    #[test]
    fn test_open_without_keyspace_returns_none() {
        let dir = tempfile::tempdir().unwrap();

        // Create a DB, add an unrelated keyspace, then close.
        {
            let db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();
            // Create a dummy keyspace that is NOT "translations" or "exons".
            let _ks = db
                .keyspace("dummy", fjall::KeyspaceCreateOptions::default)
                .unwrap();
            db.persist(fjall::PersistMode::SyncAll).unwrap();
        }

        // Reopen and verify translations/exons return None (without creating phantom keyspaces).
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();
        assert!(
            TranslationKvStore::open(&db).unwrap().is_none(),
            "TranslationKvStore::open should return None when keyspace doesn't exist"
        );
        assert!(
            ExonKvStore::open(&db).unwrap().is_none(),
            "ExonKvStore::open should return None when keyspace doesn't exist"
        );
        // Verify the keyspaces were NOT created as a side effect.
        assert!(!db.keyspace_exists(TRANSLATIONS_KEYSPACE));
        assert!(!db.keyspace_exists(EXONS_KEYSPACE));
    }

    #[test]
    fn test_translation_get_many_empty_input() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();
        let store = TranslationKvStore::create(&db).unwrap();
        store.put(&make_translation()).unwrap();

        let result = store.get_many(&[]).unwrap();
        assert!(
            result.is_empty(),
            "get_many with empty input should return empty vec"
        );
    }

    #[test]
    fn test_translation_get_many_all_missing() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();
        let store = TranslationKvStore::create(&db).unwrap();
        store.put(&make_translation()).unwrap();

        let result = store.get_many(&["MISSING1", "MISSING2"]).unwrap();
        assert!(
            result.is_empty(),
            "get_many with all missing IDs should return empty vec"
        );
    }

    #[test]
    fn test_exon_get_many_empty_input() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();
        let store = ExonKvStore::create(&db).unwrap();

        let result = store.get_many(&[]).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_exon_get_many_all_missing() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();
        let store = ExonKvStore::create(&db).unwrap();

        let result = store.get_many(&["MISSING1", "MISSING2"]).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_reopen_persistence() {
        let dir = tempfile::tempdir().unwrap();

        // Write phase.
        {
            let db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();

            let tl_store = TranslationKvStore::create(&db).unwrap();
            tl_store.put(&make_translation()).unwrap();

            let ex_store = ExonKvStore::create(&db).unwrap();
            let mut exons = make_exons();
            ex_store.put("ENST00000123456", &mut exons).unwrap();

            db.persist(fjall::PersistMode::SyncAll).unwrap();
        }

        // Read phase (new Database handle).
        {
            let db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();

            let tl_store = TranslationKvStore::open(&db).unwrap().unwrap();
            let loaded = tl_store.get("ENST00000123456").unwrap().unwrap();
            assert_eq!(loaded.cds_len, Some(1200));

            let ex_store = ExonKvStore::open(&db).unwrap().unwrap();
            let loaded = ex_store.get("ENST00000123456").unwrap();
            assert_eq!(loaded.len(), 3);
        }
    }
}
