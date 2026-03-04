//! Key encoding for position-keyed entries in the fjall KV store.
//!
//! Keys are encoded as `[chrom_code_2B][start_BE_8B][end_BE_8B]` (18 bytes)
//! so that lexicographic byte ordering matches genomic ordering: chromosomes
//! in canonical order (1-22, X, Y, MT), positions in ascending order.

use ahash::AHashMap;
use std::sync::LazyLock;

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

/// Non-canonical chromosome code: used for contigs like "GL000220.1".
/// These are sorted lexicographically after canonical chromosomes.
const NON_CANONICAL_PREFIX: u16 = 0x8000;

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

/// Encode a position key: `[2B chrom_code BE][8B start BE][8B end BE]` = 18 bytes.
///
/// Lexicographic byte ordering matches genomic ordering.
pub fn encode_position_key(chrom: &str, start: i64, end: i64) -> Vec<u8> {
    let code = chrom_to_code(chrom);
    let mut key = Vec::with_capacity(18);
    key.extend_from_slice(&code.to_be_bytes());
    key.extend_from_slice(&start.to_be_bytes());
    key.extend_from_slice(&end.to_be_bytes());
    key
}

/// Decode a position key back into `(chrom, start, end)`.
pub fn decode_position_key(key: &[u8]) -> (String, i64, i64) {
    assert!(key.len() >= 18, "position key must be at least 18 bytes");
    let code = u16::from_be_bytes([key[0], key[1]]);
    let start = i64::from_be_bytes(key[2..10].try_into().unwrap());
    let end = i64::from_be_bytes(key[10..18].try_into().unwrap());

    let chrom = CODE_TO_CHROM
        .get(&code)
        .map(|s| s.to_string())
        .unwrap_or_else(|| format!("unk_{code:04x}"));

    (chrom, start, end)
}

/// Encode a position key into a reusable buffer (hot-path variant).
///
/// The buffer is cleared and filled with 18 bytes.
pub fn encode_position_key_buf(chrom_code: u16, start: i64, end: i64, buf: &mut Vec<u8>) {
    buf.clear();
    buf.extend_from_slice(&chrom_code.to_be_bytes());
    buf.extend_from_slice(&start.to_be_bytes());
    buf.extend_from_slice(&end.to_be_bytes());
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
    fn test_position_key_roundtrip() {
        let key = encode_position_key("1", 12345, 12345);
        assert_eq!(key.len(), 18);
        let (chrom, start, end) = decode_position_key(&key);
        assert_eq!(chrom, "1");
        assert_eq!(start, 12345);
        assert_eq!(end, 12345);

        let key = encode_position_key("chr22", 100, 200);
        let (chrom, start, end) = decode_position_key(&key);
        assert_eq!(chrom, "22");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
    }

    #[test]
    fn test_position_key_chr_prefix_stripped() {
        let k1 = encode_position_key("chr1", 100, 100);
        let k2 = encode_position_key("1", 100, 100);
        assert_eq!(k1, k2);
    }

    #[test]
    fn test_position_key_lexicographic_ordering() {
        // Same chrom, start ordering
        let k1 = encode_position_key("1", 100, 100);
        let k2 = encode_position_key("1", 200, 200);
        assert!(k1 < k2);

        // Different chroms
        let k1 = encode_position_key("1", 0, 0);
        let k2 = encode_position_key("2", 0, 0);
        assert!(k1 < k2);

        // Same start, different end
        let k1 = encode_position_key("1", 100, 100);
        let k2 = encode_position_key("1", 100, 200);
        assert!(k1 < k2);

        // Chromosome 22 < X < Y < MT
        let k22 = encode_position_key("22", 0, 0);
        let kx = encode_position_key("X", 0, 0);
        let ky = encode_position_key("Y", 0, 0);
        let kmt = encode_position_key("MT", 0, 0);
        assert!(k22 < kx);
        assert!(kx < ky);
        assert!(ky < kmt);
    }

    #[test]
    fn test_position_key_buf() {
        let mut buf = Vec::new();
        let code = chrom_to_code("22");
        encode_position_key_buf(code, 12345, 12345, &mut buf);
        assert_eq!(buf.len(), 18);
        let (chrom, start, end) = decode_position_key(&buf);
        assert_eq!(chrom, "22");
        assert_eq!(start, 12345);
        assert_eq!(end, 12345);

        // Reuse buffer
        encode_position_key_buf(code, 99999, 99999, &mut buf);
        let (_, start, end) = decode_position_key(&buf);
        assert_eq!(start, 99999);
        assert_eq!(end, 99999);
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
}
