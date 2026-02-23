//! Key encoding for (chrom, window_id) pairs in the fjall KV store.
//!
//! Keys are encoded as `[chrom_canonical_2bytes][window_id_be64]` so that
//! lexicographic byte ordering matches genomic ordering: chromosomes in
//! canonical order (1-22, X, Y, MT), windows in ascending position order.

use std::collections::HashMap;
use std::sync::LazyLock;

/// Default window size in base pairs (1 Mb).
pub const DEFAULT_WINDOW_SIZE: u64 = 1_000_000;

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
        assert_eq!(window_id_for_position(0, DEFAULT_WINDOW_SIZE), 0);
        assert_eq!(window_id_for_position(1, DEFAULT_WINDOW_SIZE), 0);
        assert_eq!(window_id_for_position(999_999, DEFAULT_WINDOW_SIZE), 0);
        assert_eq!(window_id_for_position(1_000_000, DEFAULT_WINDOW_SIZE), 1);
        assert_eq!(window_id_for_position(1_500_000, DEFAULT_WINDOW_SIZE), 1);
        assert_eq!(window_id_for_position(2_000_000, DEFAULT_WINDOW_SIZE), 2);
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
}
