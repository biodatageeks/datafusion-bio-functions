//! Key encoding for position-keyed entries in the fjall KV store.
//!
//! Keys are encoded as `[chrom_code_2B][start_BE_8B]` (10 bytes)
//! so that lexicographic byte ordering matches genomic ordering: chromosomes
//! in canonical order (1-22, X, Y, MT), positions in ascending order.
//! All variants at the same `(chrom, start)` are merged into a single entry.

use ahash::AHashMap;
use std::sync::{LazyLock, RwLock};

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

/// First code assigned to non-canonical contigs.
/// Canonical codes: 1..=22 (autosomes), 0x0017=X, 0x0018=Y, 0x0019=MT (=25).
/// Non-canonical range starts at 26 to avoid colliding with MT.
const NON_CANONICAL_START: u16 = 26;

/// Process-global registry for non-canonical contig → code mappings.
///
/// Populated at cache-build time via [`register_non_canonical_contigs`] and
/// at cache-open time via [`load_non_canonical_registry`]. Lookup is O(1)
/// via the RwLock read path.
///
/// **Limitation:** This is process-global, not per-store. Opening multiple
/// caches with different contig sets in the same process will overwrite each
/// other's mappings. In practice this is fine because `build_cache` and VEP
/// annotation operate on a single cache at a time. If multi-cache support is
/// needed, this should be refactored to per-store scoping.
static NON_CANONICAL_REGISTRY: LazyLock<RwLock<AHashMap<String, u16>>> =
    LazyLock::new(|| RwLock::new(AHashMap::new()));

/// Reverse mapping: code → contig name for non-canonical contigs.
static NON_CANONICAL_CODE_TO_NAME: LazyLock<RwLock<AHashMap<u16, String>>> =
    LazyLock::new(|| RwLock::new(AHashMap::new()));

/// Map a chromosome name (with or without "chr" prefix) to its 2-byte code.
///
/// Canonical chromosomes (1-22, X, Y, MT) get codes 0x0001..0x0019.
/// Non-canonical contigs are looked up in the registry (populated at build/open time).
/// Falls back to FNV-1a hash if the registry has not been populated.
pub fn chrom_to_code(chrom: &str) -> u16 {
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Some(&code) = CHROM_TO_CODE.get(chrom_bare) {
        return code;
    }
    // Check the non-canonical registry
    if let Ok(reg) = NON_CANONICAL_REGISTRY.read() {
        if let Some(&code) = reg.get(chrom_bare) {
            return code;
        }
    }
    // Fallback: FNV-1a hash (used before registry is populated, e.g. during
    // initial schema discovery or if the cache was built without a registry)
    fnv1a_code(chrom_bare)
}

/// Map a 2-byte chromosome code back to its name.
///
/// Returns `None` for codes not found in either the canonical or
/// non-canonical registry.
pub fn code_to_chrom(code: u16) -> Option<String> {
    if let Some(&name) = CODE_TO_CHROM.get(&code) {
        return Some(name.to_string());
    }
    if let Ok(reg) = NON_CANONICAL_CODE_TO_NAME.read() {
        if let Some(name) = reg.get(&code) {
            return Some(name.clone());
        }
    }
    None
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

    let chrom = code_to_chrom(code).unwrap_or_else(|| format!("unk_{code:04x}"));

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

/// Register non-canonical contigs with collision-free sequential codes.
///
/// Contigs are sorted lexicographically and assigned codes starting from
/// [`NON_CANONICAL_START`] (26). This guarantees no collisions and produces
/// a deterministic mapping for the same input set.
///
/// Returns the mapping as a `Vec<(String, u16)>` for persistence.
///
/// # Panics
///
/// Panics if there are more non-canonical contigs than fit in the u16 range
/// (65,510 slots), which is practically impossible.
pub fn register_non_canonical_contigs(contigs: &[&str]) -> Vec<(String, u16)> {
    let mut non_canonical: Vec<&str> = contigs
        .iter()
        .copied()
        .filter(|c| !CHROM_TO_CODE.contains_key(c))
        .collect();
    non_canonical.sort();
    non_canonical.dedup();

    assert!(
        non_canonical.len() <= (u16::MAX - NON_CANONICAL_START + 1) as usize,
        "too many non-canonical contigs: {} (max {})",
        non_canonical.len(),
        u16::MAX - NON_CANONICAL_START + 1
    );

    let mapping: Vec<(String, u16)> = non_canonical
        .iter()
        .enumerate()
        .map(|(i, &name)| (name.to_string(), NON_CANONICAL_START + i as u16))
        .collect();

    // Populate the global registry
    let mut reg = NON_CANONICAL_REGISTRY.write().unwrap();
    reg.clear();
    let mut rev = NON_CANONICAL_CODE_TO_NAME.write().unwrap();
    rev.clear();
    for (name, code) in &mapping {
        reg.insert(name.clone(), *code);
        rev.insert(*code, name.clone());
    }

    mapping
}

/// Load a previously persisted contig→code mapping into the global registry.
///
/// Called at cache-open time from `VepKvStore::open()`.
pub fn load_non_canonical_registry(mapping: &[(String, u16)]) {
    let mut reg = NON_CANONICAL_REGISTRY.write().unwrap();
    reg.clear();
    let mut rev = NON_CANONICAL_CODE_TO_NAME.write().unwrap();
    rev.clear();
    for (name, code) in mapping {
        reg.insert(name.clone(), *code);
        rev.insert(*code, name.clone());
    }
}

/// Serialize a contig→code mapping to JSON bytes for storage in fjall meta.
pub fn serialize_contig_codes(mapping: &[(String, u16)]) -> Vec<u8> {
    serde_json::to_vec(mapping).expect("contig codes should be JSON-serializable")
}

/// Deserialize a contig→code mapping from JSON bytes.
pub fn deserialize_contig_codes(bytes: &[u8]) -> Option<Vec<(String, u16)>> {
    serde_json::from_slice(bytes).ok()
}

/// FNV-1a hash for non-canonical chromosomes (legacy fallback).
///
/// Used when the registry has not been populated (e.g. caches built before
/// the registry was introduced). May produce collisions with large contig sets.
fn fnv1a_code(chrom: &str) -> u16 {
    let mut hash: u32 = 0x811c_9dc5;
    for b in chrom.bytes() {
        hash ^= b as u32;
        hash = hash.wrapping_mul(0x0100_0193);
    }
    let range = (u16::MAX as u32) - (NON_CANONICAL_START as u32) + 1;
    NON_CANONICAL_START + (hash % range) as u16
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
            assert_eq!(decoded, *chrom, "roundtrip failed for {chrom}");
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
                "canonical chrom '{chrom}' has code {code} >= NON_CANONICAL_START {NON_CANONICAL_START}"
            );
        }
    }

    #[test]
    fn test_register_non_canonical_contigs_no_collisions() {
        let contigs = vec!["LRG_1278", "LRG_8", "GL000220.1", "KI270733.1"];
        let mapping = register_non_canonical_contigs(&contigs);

        // All codes must be unique
        let codes: Vec<u16> = mapping.iter().map(|(_, c)| *c).collect();
        let unique: std::collections::HashSet<u16> = codes.iter().copied().collect();
        assert_eq!(codes.len(), unique.len(), "codes must be unique");

        // All codes >= NON_CANONICAL_START
        assert!(codes.iter().all(|&c| c >= NON_CANONICAL_START));

        // Codes are sequential (sorted input)
        let mut sorted_names: Vec<&str> = contigs.to_vec();
        sorted_names.sort();
        sorted_names.dedup();
        for (i, (name, code)) in mapping.iter().enumerate() {
            assert_eq!(*code, NON_CANONICAL_START + i as u16);
            assert_eq!(name, sorted_names[i]);
        }
    }

    #[test]
    fn test_register_solves_real_collision() {
        // LRG_1278 and LRG_8 collide under FNV-1a
        assert_eq!(fnv1a_code("LRG_1278"), fnv1a_code("LRG_8"));

        // After registration, they get distinct codes
        let mapping = register_non_canonical_contigs(&["LRG_1278", "LRG_8"]);
        let code_1278 = mapping.iter().find(|(n, _)| n == "LRG_1278").unwrap().1;
        let code_8 = mapping.iter().find(|(n, _)| n == "LRG_8").unwrap().1;
        assert_ne!(code_1278, code_8);

        // chrom_to_code() now returns registry codes
        assert_eq!(chrom_to_code("LRG_1278"), code_1278);
        assert_eq!(chrom_to_code("LRG_8"), code_8);
    }

    #[test]
    fn test_register_skips_canonical_contigs() {
        let mapping = register_non_canonical_contigs(&["1", "GL000220.1", "X"]);
        // Only GL000220.1 should be in the mapping
        assert_eq!(mapping.len(), 1);
        assert_eq!(mapping[0].0, "GL000220.1");
    }

    #[test]
    fn test_register_deduplicates() {
        let mapping = register_non_canonical_contigs(&["GL000220.1", "GL000220.1", "KI270733.1"]);
        assert_eq!(mapping.len(), 2);
    }

    #[test]
    fn test_load_registry_roundtrip() {
        let contigs = vec!["LRG_1278", "LRG_8", "GL000220.1"];
        let mapping = register_non_canonical_contigs(&contigs);

        // Serialize and deserialize
        let bytes = serialize_contig_codes(&mapping);
        let loaded = deserialize_contig_codes(&bytes).unwrap();
        assert_eq!(mapping, loaded);

        // Load into registry and verify
        load_non_canonical_registry(&loaded);
        for (name, code) in &loaded {
            assert_eq!(chrom_to_code(name), *code);
            assert_eq!(code_to_chrom(*code).unwrap(), *name);
        }
    }

    #[test]
    fn test_sequential_codes_are_ascending() {
        let contigs = vec!["ZZZ", "AAA", "MMM"];
        let mapping = register_non_canonical_contigs(&contigs);
        // Sorted lexicographically: AAA, MMM, ZZZ
        let codes: Vec<u16> = mapping.iter().map(|(_, c)| *c).collect();
        for w in codes.windows(2) {
            assert!(w[0] < w[1], "codes must be strictly ascending");
        }
    }

    #[test]
    fn test_code_to_chrom_non_canonical_after_register() {
        register_non_canonical_contigs(&["HG1012_PATCH"]);
        let code = chrom_to_code("HG1012_PATCH");
        let name = code_to_chrom(code).unwrap();
        assert_eq!(name, "HG1012_PATCH");
    }
}
