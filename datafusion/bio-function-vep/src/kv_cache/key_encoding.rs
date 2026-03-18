//! Key encoding for position-keyed entries in the fjall KV store.
//!
//! Keys are encoded as `[chrom_code_2B][start_BE_8B]` (10 bytes)
//! so that lexicographic byte ordering matches genomic ordering: chromosomes
//! in canonical order (1-22, X, Y, MT), positions in ascending order.
//! All variants at the same `(chrom, start)` are merged into a single entry.

use ahash::AHashMap;
use std::collections::HashMap;
use std::sync::{LazyLock, Mutex};

/// Canonical chromosome order: 1-22, X, Y, MT.
/// Maps chromosome name (without "chr" prefix) to a 2-byte big-endian code.
static CHROM_TO_CODE: LazyLock<AHashMap<&'static str, u16>> = LazyLock::new(|| {
    let mut m = AHashMap::default();
    for i in 1..=22u16 {
        m.insert(CHROM_NAMES[i as usize - 1], i);
    }
    m.insert("X", 0x0017);
    m.insert("Y", 0x0018);
    m.insert("MT", 0x0019);
    m
});

static CODE_TO_CHROM: LazyLock<AHashMap<u16, &'static str>> =
    LazyLock::new(|| CHROM_TO_CODE.iter().map(|(&k, &v)| (v, k)).collect());

const CHROM_NAMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22",
];

/// First code assigned to non-canonical contigs (canonical use 1..=25).
const NON_CANONICAL_START: u16 = 26;

/// Map a chromosome name (with or without "chr" prefix) to its 2-byte code.
///
/// Canonical chromosomes (1-22, X, Y, MT) get codes 0x0001..0x0019.
/// Non-canonical contigs get deterministic codes >= 0x8000.
pub fn chrom_to_code(chrom: &str) -> u16 {
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    CHROM_TO_CODE
        .get(chrom_bare)
        .copied()
        .unwrap_or_else(|| non_canonical_code(chrom_bare))
}

/// Map a 2-byte chromosome code back to its name.
///
/// Returns `None` for non-canonical codes.
pub fn code_to_chrom(code: u16) -> Option<&'static str> {
    CODE_TO_CHROM.get(&code).copied()
}

/// Encode a position key: `[2B chrom_code BE][8B start BE]` = 10 bytes.
///
/// Lexicographic byte ordering matches genomic ordering.
pub fn encode_position_key(chrom: &str, start: i64) -> Vec<u8> {
    let code = chrom_to_code(chrom);
    let mut key = Vec::with_capacity(10);
    key.extend_from_slice(&code.to_be_bytes());
    key.extend_from_slice(&start.to_be_bytes());
    key
}

/// Decode a position key back into `(chrom, start)`.
pub fn decode_position_key(key: &[u8]) -> (String, i64) {
    assert!(key.len() >= 10, "position key must be at least 10 bytes");
    let code = u16::from_be_bytes([key[0], key[1]]);
    let start = i64::from_be_bytes(key[2..10].try_into().unwrap());

    let chrom = CODE_TO_CHROM
        .get(&code)
        .map(|s| s.to_string())
        .unwrap_or_else(|| format!("unk_{code:04x}"));

    (chrom, start)
}

/// Encode a position key into a reusable buffer (hot-path variant).
///
/// The buffer is cleared and filled with 10 bytes.
pub fn encode_position_key_buf(chrom_code: u16, start: i64, buf: &mut Vec<u8>) {
    buf.clear();
    buf.extend_from_slice(&chrom_code.to_be_bytes());
    buf.extend_from_slice(&start.to_be_bytes());
}

/// Global registry tracking which contig name was assigned which non-canonical code.
/// Used to detect hash collisions during cache creation.
static NON_CANONICAL_REGISTRY: LazyLock<Mutex<HashMap<u16, String>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

/// Deterministic code for non-canonical chromosomes.
/// Uses a hash-based approach to assign stable codes in the range
/// `NON_CANONICAL_START..=u16::MAX` (65510 possible values), which is
/// much wider than the previous 15-bit (32768) space and reduces
/// collision probability for non-canonical contig names.
///
/// Panics if two different contig names hash to the same code (collision).
fn non_canonical_code(chrom: &str) -> u16 {
    // FNV-1a 32-bit hash for better distribution.
    let mut hash: u32 = 0x811c_9dc5;
    for b in chrom.bytes() {
        hash ^= b as u32;
        hash = hash.wrapping_mul(0x0100_0193);
    }
    let range = (u16::MAX as u32) - (NON_CANONICAL_START as u32) + 1;
    let code = NON_CANONICAL_START + (hash % range) as u16;

    // Check for collisions: different contig name mapping to the same code.
    let mut registry = NON_CANONICAL_REGISTRY.lock().unwrap();
    match registry.entry(code) {
        std::collections::hash_map::Entry::Occupied(entry) => {
            if entry.get() != chrom {
                panic!(
                    "non-canonical contig hash collision: '{}' and '{}' both map to code {code}",
                    entry.get(),
                    chrom
                );
            }
        }
        std::collections::hash_map::Entry::Vacant(entry) => {
            entry.insert(chrom.to_string());
        }
    }

    code
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_position_key_roundtrip() {
        let key = encode_position_key("1", 12345);
        assert_eq!(key.len(), 10);
        let (chrom, start) = decode_position_key(&key);
        assert_eq!(chrom, "1");
        assert_eq!(start, 12345);

        let key = encode_position_key("chr22", 100);
        let (chrom, start) = decode_position_key(&key);
        assert_eq!(chrom, "22");
        assert_eq!(start, 100);
    }

    #[test]
    fn test_position_key_chr_prefix_stripped() {
        let k1 = encode_position_key("chr1", 100);
        let k2 = encode_position_key("1", 100);
        assert_eq!(k1, k2);
    }

    #[test]
    fn test_position_key_lexicographic_ordering() {
        // Same chrom, start ordering
        let k1 = encode_position_key("1", 100);
        let k2 = encode_position_key("1", 200);
        assert!(k1 < k2);

        // Different chroms
        let k1 = encode_position_key("1", 0);
        let k2 = encode_position_key("2", 0);
        assert!(k1 < k2);

        // Chromosome 22 < X < Y < MT
        let k22 = encode_position_key("22", 0);
        let kx = encode_position_key("X", 0);
        let ky = encode_position_key("Y", 0);
        let kmt = encode_position_key("MT", 0);
        assert!(k22 < kx);
        assert!(kx < ky);
        assert!(ky < kmt);
    }

    #[test]
    fn test_position_key_buf() {
        let mut buf = Vec::new();
        let code = chrom_to_code("22");
        encode_position_key_buf(code, 12345, &mut buf);
        assert_eq!(buf.len(), 10);
        let (chrom, start) = decode_position_key(&buf);
        assert_eq!(chrom, "22");
        assert_eq!(start, 12345);

        // Reuse buffer
        encode_position_key_buf(code, 99999, &mut buf);
        let (_, start) = decode_position_key(&buf);
        assert_eq!(start, 99999);
    }

    #[test]
    fn test_chrom_code_roundtrip() {
        for chrom in CHROM_NAMES.iter().chain(&["X", "Y", "MT"]) {
            let code = chrom_to_code(chrom);
            let decoded = code_to_chrom(code).unwrap();
            assert_eq!(&decoded, chrom, "roundtrip failed for {chrom}");
        }
    }

    #[test]
    fn test_chrom_code_chr_prefix_stripped() {
        assert_eq!(chrom_to_code("chr1"), chrom_to_code("1"));
        assert_eq!(chrom_to_code("chrX"), chrom_to_code("X"));
    }

    #[test]
    fn test_no_canonical_code_in_non_canonical_range() {
        for chrom in CHROM_NAMES.iter().chain(&["X", "Y", "MT"]) {
            let code = chrom_to_code(chrom);
            assert!(
                code < NON_CANONICAL_START,
                "canonical chrom '{chrom}' has code {code} which falls in non-canonical range (>= {NON_CANONICAL_START})"
            );
        }
    }
}
