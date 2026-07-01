use std::fmt;

use hashbrown::HashMap;
use nohash_hasher::BuildNoHashHasher;
use smallvec::SmallVec;

use crate::lance_cache::key_encoding::chrom_to_code;

const POSITION_BITS: u64 = 48;
const MAX_POSITION: u64 = (1_u64 << POSITION_BITS) - 1;

pub type VariantIndex = HashMap<u64, SmallVec<[u32; 1]>, BuildNoHashHasher<u64>>;
pub type PositionIndex = HashMap<u64, SmallVec<[u32; 4]>, BuildNoHashHasher<u64>>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WarmKeyError {
    NegativeStart(i64),
    StartOutOfRange(i64),
}

impl fmt::Display for WarmKeyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NegativeStart(start) => write!(f, "variant start must be non-negative: {start}"),
            Self::StartOutOfRange(start) => {
                write!(
                    f,
                    "variant start exceeds 48-bit warm-cache key range: {start}"
                )
            }
        }
    }
}

impl std::error::Error for WarmKeyError {}

pub type Result<T> = std::result::Result<T, WarmKeyError>;

pub fn position_key(chrom: &str, start: i64) -> Result<u64> {
    let chrom_code = chrom_to_code(chrom.strip_prefix("chr").unwrap_or(chrom));
    position_key_from_code(chrom_code, start)
}

pub fn position_key_from_code(chrom_code: u16, start: i64) -> Result<u64> {
    if start < 0 {
        return Err(WarmKeyError::NegativeStart(start));
    }

    let start = start as u64;
    if start > MAX_POSITION {
        return Err(WarmKeyError::StartOutOfRange(start as i64));
    }

    Ok(((chrom_code as u64) << POSITION_BITS) | start)
}

pub fn variant_key(chrom: &str, start: i64, reference: &str, alternate: &str) -> Result<u64> {
    let pos = position_key(chrom, start)?;
    Ok(variant_key_from_position(pos, reference, alternate))
}

pub fn variant_key_from_position(position_key: u64, reference: &str, alternate: &str) -> u64 {
    const STACK_CAP: usize = 256;
    let needed = 8 + 4 + reference.len() + 4 + alternate.len();

    if needed <= STACK_CAP {
        let mut bytes = [0_u8; STACK_CAP];
        let mut offset = 0usize;
        bytes[offset..offset + 8].copy_from_slice(&position_key.to_le_bytes());
        offset += 8;
        bytes[offset..offset + 4].copy_from_slice(&(reference.len() as u32).to_le_bytes());
        offset += 4;
        bytes[offset..offset + reference.len()].copy_from_slice(reference.as_bytes());
        offset += reference.len();
        bytes[offset..offset + 4].copy_from_slice(&(alternate.len() as u32).to_le_bytes());
        offset += 4;
        bytes[offset..offset + alternate.len()].copy_from_slice(alternate.as_bytes());
        offset += alternate.len();
        rapidhash::v3::rapidhash_v3(&bytes[..offset])
    } else {
        stable_variant_hash(position_key, reference, alternate)
    }
}

pub fn variant_keys_from_allele_string(
    chrom: &str,
    start: i64,
    allele_string: &str,
) -> Result<Vec<u64>> {
    let Some((reference, alternates)) = allele_string.split_once('/') else {
        return Ok(Vec::new());
    };

    let pos = position_key(chrom, start)?;
    Ok(alternates
        .split('/')
        .filter(|alternate| !alternate.is_empty())
        .map(|alternate| variant_key_from_position(pos, reference, alternate))
        .collect())
}

pub fn new_variant_index(capacity: usize) -> VariantIndex {
    VariantIndex::with_capacity_and_hasher(capacity, BuildNoHashHasher::default())
}

pub fn new_position_index(capacity: usize) -> PositionIndex {
    PositionIndex::with_capacity_and_hasher(capacity, BuildNoHashHasher::default())
}

fn stable_variant_hash(position_key: u64, reference: &str, alternate: &str) -> u64 {
    let mut bytes = Vec::with_capacity(16 + reference.len() + alternate.len());
    bytes.extend_from_slice(&position_key.to_le_bytes());
    bytes.extend_from_slice(&(reference.len() as u32).to_le_bytes());
    bytes.extend_from_slice(reference.as_bytes());
    bytes.extend_from_slice(&(alternate.len() as u32).to_le_bytes());
    bytes.extend_from_slice(alternate.as_bytes());
    rapidhash::v3::rapidhash_v3(&bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn position_key_orders_by_chrom_and_start() {
        let chr1_100 = position_key("1", 100).unwrap();
        let chr1_200 = position_key("1", 200).unwrap();
        let chr2_1 = position_key("2", 1).unwrap();

        assert!(chr1_100 < chr1_200);
        assert!(chr1_200 < chr2_1);
    }

    #[test]
    fn position_key_accepts_chr_prefix() {
        assert_eq!(position_key("1", 100), position_key("chr1", 100));
    }

    #[test]
    fn position_key_rejects_negative_start() {
        assert!(position_key("1", -1).is_err());
    }

    #[test]
    fn variant_keys_include_each_alt_from_allele_string() {
        let keys = variant_keys_from_allele_string("1", 101, "A/G/T").unwrap();

        assert_eq!(keys.len(), 2);
        assert_ne!(keys[0], keys[1]);
    }

    #[test]
    fn variant_key_is_stable_for_same_tuple() {
        let key1 = variant_key("1", 101, "A", "G").unwrap();
        let key2 = variant_key("chr1", 101, "A", "G").unwrap();

        assert_eq!(key1, key2);
    }

    #[test]
    fn variant_key_from_position_matches_stable_hash() {
        let pos = position_key("1", 101).unwrap();
        assert_eq!(
            variant_key_from_position(pos, "A", "G"),
            stable_variant_hash(pos, "A", "G")
        );
    }

    #[test]
    fn warm_variant_index_uses_existing_variant_keys() {
        let mut index = new_variant_index(2);
        index.entry(42).or_default().push(7);
        index.entry(42).or_default().push(11);

        assert_eq!(index.get(&42).unwrap().as_slice(), &[7, 11]);
        assert!(index.get(&43).is_none());
    }

    #[test]
    fn warm_position_index_uses_existing_position_keys() {
        let mut index = new_position_index(2);
        index.entry(42).or_default().push(7);
        index.entry(42).or_default().push(11);

        assert_eq!(index.get(&42).unwrap().as_slice(), &[7, 11]);
        assert!(index.get(&43).is_none());
    }
}
