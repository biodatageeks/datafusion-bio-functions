//! Key encoding for (chrom, window_id) pairs in the fjall KV store.
//!
//! Keys are encoded as `[chrom_canonical_2bytes][window_id_be64]` so that
//! lexicographic byte ordering matches genomic ordering: chromosomes in
//! canonical order (1-22, X, Y, MT), windows in ascending position order.

use std::collections::HashMap;
use std::sync::LazyLock;

/// Default window size in base pairs (10 Kb).
///
/// Benchmarked on chr22 (15.1M VEP variants, 78 columns):
///   10kb → ~4K variants/window, 13ms read+index
///  100kb → ~37K variants/window, 84ms read+index
///    1Mb → ~369K variants/window, 453ms read+index
pub const DEFAULT_WINDOW_SIZE: u64 = 10_000;

/// Canonical chromosome order: 1-22, X, Y, MT.
/// Maps chromosome name (without "chr" prefix) to a 2-byte big-endian code.
static CHROM_TO_CODE: LazyLock<HashMap<&'static str, u16>> = LazyLock::new(|| {
    let mut m = HashMap::new();
    for i in 1..=22u16 {
        m.insert(CHROM_NAMES[i as usize - 1], i);
    }
    m.insert("X", 0x0017);
    m.insert("Y", 0x0018);
    m.insert("MT", 0x0019);
    m
});

static CODE_TO_CHROM: LazyLock<HashMap<u16, &'static str>> =
    LazyLock::new(|| CHROM_TO_CODE.iter().map(|(&k, &v)| (v, k)).collect());

const CHROM_NAMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22",
];

/// Non-canonical chromosome code: used for contigs like "GL000220.1".
/// These are sorted lexicographically after canonical chromosomes.
const NON_CANONICAL_PREFIX: u16 = 0x8000;

/// Entry type stored in the last byte of a v1 key.
///
/// v1 keys are 11 bytes: `[2B chrom][8B window_id][1B entry_type]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntryType {
    /// Compact binary position index (start, end, allele offsets).
    PositionIndex,
    /// Single column stored as Arrow IPC. The `u8` is the column index (0-based).
    Column(u8),
    /// Legacy monolithic Arrow IPC batch (all columns). Used for v0 backward compat reads.
    LegacyMonolithic,
}

impl EntryType {
    /// Maximum column index supported by v1 key encoding.
    ///
    /// `entry_type` byte reserves:
    /// - `0x00` for `PositionIndex`
    /// - `0xFF` for `LegacyMonolithic`
    ///
    /// So column entries can use `0x01..=0xFE`, which maps to
    /// column indices `0..=253` (254 total columns).
    pub const MAX_COLUMN_INDEX: u8 = 0xFD;

    fn to_byte(self) -> u8 {
        match self {
            Self::PositionIndex => 0x00,
            Self::Column(idx) => idx + 1, // 0x01..=0xFE
            Self::LegacyMonolithic => 0xFF,
        }
    }

    fn from_byte(b: u8) -> Self {
        match b {
            0x00 => Self::PositionIndex,
            0xFF => Self::LegacyMonolithic,
            col => Self::Column(col - 1), // 0x01..=0xFE → column 0..=253
        }
    }
}

/// Encode a (chrom, window_id, entry_type) triple into an 11-byte key.
///
/// Format: `[2-byte chrom code][8-byte big-endian window_id][1-byte entry_type]`
pub fn encode_entry_key(chrom: &str, window_id: u64, entry: EntryType) -> Vec<u8> {
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let code = CHROM_TO_CODE
        .get(chrom_bare)
        .copied()
        .unwrap_or_else(|| non_canonical_code(chrom_bare));

    let mut key = Vec::with_capacity(11);
    key.extend_from_slice(&code.to_be_bytes());
    key.extend_from_slice(&window_id.to_be_bytes());
    key.push(entry.to_byte());
    key
}

/// Decode an 11-byte key back into (chrom, window_id, entry_type).
pub fn decode_entry_key(key: &[u8]) -> (String, u64, EntryType) {
    assert!(key.len() >= 11, "entry key must be at least 11 bytes");
    let code = u16::from_be_bytes([key[0], key[1]]);
    let window_id = u64::from_be_bytes(key[2..10].try_into().unwrap());
    let entry = EntryType::from_byte(key[10]);

    let chrom = CODE_TO_CHROM
        .get(&code)
        .map(|s| s.to_string())
        .unwrap_or_else(|| format!("unk_{code:04x}"));

    (chrom, window_id, entry)
}

/// Compute the window ID for a given genomic position.
pub fn window_id_for_position(pos: i64, window_size: u64) -> u64 {
    if pos <= 0 {
        return 0;
    }
    pos as u64 / window_size
}

/// Encode a (chrom, window_id) pair into a byte key.
///
/// Format: `[2-byte chrom code][8-byte big-endian window_id]`
pub fn encode_window_key(chrom: &str, window_id: u64) -> Vec<u8> {
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let code = CHROM_TO_CODE
        .get(chrom_bare)
        .copied()
        .unwrap_or_else(|| non_canonical_code(chrom_bare));

    let mut key = Vec::with_capacity(10);
    key.extend_from_slice(&code.to_be_bytes());
    key.extend_from_slice(&window_id.to_be_bytes());
    key
}

/// Decode a byte key back into (chrom, window_id).
pub fn decode_window_key(key: &[u8]) -> (String, u64) {
    assert!(key.len() >= 10, "key must be at least 10 bytes");
    let code = u16::from_be_bytes([key[0], key[1]]);
    let window_id = u64::from_be_bytes(key[2..10].try_into().unwrap());

    let chrom = CODE_TO_CHROM
        .get(&code)
        .map(|s| s.to_string())
        .unwrap_or_else(|| format!("unk_{code:04x}"));

    (chrom, window_id)
}

/// Encode a 2-byte prefix for a chromosome (for prefix scans).
pub fn chrom_prefix(chrom: &str) -> [u8; 2] {
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let code = CHROM_TO_CODE
        .get(chrom_bare)
        .copied()
        .unwrap_or_else(|| non_canonical_code(chrom_bare));
    code.to_be_bytes()
}

/// Deterministic code for non-canonical chromosomes.
/// Uses a hash-based approach to assign stable codes >= 0x8000.
fn non_canonical_code(chrom: &str) -> u16 {
    let mut hash: u16 = 0;
    for b in chrom.bytes() {
        hash = hash.wrapping_mul(31).wrapping_add(b as u16);
    }
    NON_CANONICAL_PREFIX | (hash & 0x7FFF)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_window_id_for_position() {
        const W: u64 = 1_000_000;
        assert_eq!(window_id_for_position(0, W), 0);
        assert_eq!(window_id_for_position(1, W), 0);
        assert_eq!(window_id_for_position(999_999, W), 0);
        assert_eq!(window_id_for_position(1_000_000, W), 1);
        assert_eq!(window_id_for_position(1_500_000, W), 1);
        assert_eq!(window_id_for_position(2_000_000, W), 2);

        // Also test with the actual default
        assert_eq!(window_id_for_position(0, DEFAULT_WINDOW_SIZE), 0);
        assert_eq!(window_id_for_position(9_999, DEFAULT_WINDOW_SIZE), 0);
        assert_eq!(window_id_for_position(10_000, DEFAULT_WINDOW_SIZE), 1);
    }

    #[test]
    fn test_roundtrip_canonical_chromosomes() {
        for chrom in CHROM_NAMES.iter().chain(&["X", "Y", "MT"]) {
            let key = encode_window_key(chrom, 42);
            let (decoded_chrom, decoded_wid) = decode_window_key(&key);
            assert_eq!(&decoded_chrom, chrom, "roundtrip failed for {chrom}");
            assert_eq!(decoded_wid, 42);
        }
    }

    #[test]
    fn test_chr_prefix_stripped() {
        let key1 = encode_window_key("chr1", 5);
        let key2 = encode_window_key("1", 5);
        assert_eq!(key1, key2);
    }

    #[test]
    fn test_lexicographic_ordering() {
        // Chromosome 1 window 0 < Chromosome 2 window 0
        let k1 = encode_window_key("1", 0);
        let k2 = encode_window_key("2", 0);
        assert!(k1 < k2);

        // Same chromosome, window 0 < window 1
        let k3 = encode_window_key("1", 0);
        let k4 = encode_window_key("1", 1);
        assert!(k3 < k4);

        // Chromosome 22 < X < Y < MT
        let k22 = encode_window_key("22", 0);
        let kx = encode_window_key("X", 0);
        let ky = encode_window_key("Y", 0);
        let kmt = encode_window_key("MT", 0);
        assert!(k22 < kx);
        assert!(kx < ky);
        assert!(ky < kmt);
    }

    #[test]
    fn test_chrom_prefix_scan() {
        let prefix = chrom_prefix("1");
        let key_in = encode_window_key("1", 100);
        let key_out = encode_window_key("2", 0);
        assert!(key_in.starts_with(&prefix));
        assert!(!key_out.starts_with(&prefix));
    }

    #[test]
    fn test_negative_position() {
        assert_eq!(window_id_for_position(-1, DEFAULT_WINDOW_SIZE), 0);
    }

    #[test]
    fn test_entry_type_roundtrip() {
        let types = [
            EntryType::PositionIndex,
            EntryType::Column(0),
            EntryType::Column(42),
            EntryType::Column(EntryType::MAX_COLUMN_INDEX),
            EntryType::LegacyMonolithic,
        ];
        for et in types {
            assert_eq!(EntryType::from_byte(et.to_byte()), et);
        }
    }

    #[test]
    fn test_max_column_index_maps_to_last_non_reserved_byte() {
        assert_eq!(
            EntryType::Column(EntryType::MAX_COLUMN_INDEX).to_byte(),
            0xFE
        );
        assert_eq!(
            EntryType::from_byte(0xFE),
            EntryType::Column(EntryType::MAX_COLUMN_INDEX)
        );
    }

    #[test]
    fn test_entry_key_roundtrip() {
        let key = encode_entry_key("1", 42, EntryType::PositionIndex);
        assert_eq!(key.len(), 11);
        let (chrom, wid, entry) = decode_entry_key(&key);
        assert_eq!(chrom, "1");
        assert_eq!(wid, 42);
        assert_eq!(entry, EntryType::PositionIndex);

        let key = encode_entry_key("chr22", 100, EntryType::Column(5));
        let (chrom, wid, entry) = decode_entry_key(&key);
        assert_eq!(chrom, "22");
        assert_eq!(wid, 100);
        assert_eq!(entry, EntryType::Column(5));
    }

    #[test]
    fn test_entry_key_ordering() {
        // Same chrom + window, different entry types: position index < column < legacy
        let k_pos = encode_entry_key("1", 0, EntryType::PositionIndex);
        let k_col0 = encode_entry_key("1", 0, EntryType::Column(0));
        let k_col5 = encode_entry_key("1", 0, EntryType::Column(5));
        let k_legacy = encode_entry_key("1", 0, EntryType::LegacyMonolithic);
        assert!(k_pos < k_col0);
        assert!(k_col0 < k_col5);
        assert!(k_col5 < k_legacy);
    }
}
