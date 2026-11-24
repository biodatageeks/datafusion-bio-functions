//! Sequence hashing for deduplication.

use twox_hash::XxHash3_64;

use super::complement::{Alphabet, revcomp_into};

/// The XXH3 64-bit hash (seed 0) of the sequence bytes.
///
/// With `canonical`, the smaller of the hash of `seq` and the hash of its
/// DNA reverse complement, so a sequence and its reverse complement hash
/// equally. `scratch` holds the reverse complement between calls.
///
/// Case is not folded. The values are stable across releases and platforms.
pub fn seq_hash(seq: &[u8], canonical: bool, scratch: &mut Vec<u8>) -> u64 {
    let forward = XxHash3_64::oneshot(seq);
    if !canonical {
        return forward;
    }
    scratch.clear();
    // Non-strict reverse complement cannot fail.
    let _ = revcomp_into(seq, Alphabet::Dna, false, scratch);
    forward.min(XxHash3_64::oneshot(scratch))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::{Alphabet, revcomp_into};

    fn hash(seq: &str, canonical: bool) -> u64 {
        seq_hash(seq.as_bytes(), canonical, &mut Vec::new())
    }

    #[test]
    fn xxh3_64_reference_values() {
        // xxhash.xxh3_64_intdigest (python-xxhash 4.1.0, libxxhash 0.8.4).
        // These values are part of the API and must never change.
        for (seq, expected) in [
            ("", 3244421341483603138u64),
            ("A", 15047818145317598341),
            ("ACGT", 17060466519357990319),
            ("acgt", 251730389294227419),
            ("AAAACCCCGGGGTTTT", 71406446880636177),
            ("TTAGGCCAT", 17056848131934351943),
            ("ATGGCCTAA", 2721565146873327973),
        ] {
            assert_eq!(hash(seq, false), expected, "{seq:?}");
        }
        assert_eq!(hash(&"ACGTN".repeat(40), false), 2378561995346465887);
    }

    #[test]
    fn canonical_is_the_smaller_of_both_strands() {
        // TTAGGCCAT is the reverse complement of ATGGCCTAA.
        assert_eq!(hash("ATGGCCTAA", true), 2721565146873327973);
        assert_eq!(hash("TTAGGCCAT", true), 2721565146873327973);
        // A reverse-complement palindrome.
        assert_eq!(hash("ACGT", true), hash("ACGT", false));
        assert_eq!(hash("", true), hash("", false));
    }

    #[test]
    fn scratch_is_reusable() {
        let mut scratch = Vec::new();
        let a = seq_hash(b"TTAGGCCAT", true, &mut scratch);
        let b = seq_hash(b"TT", true, &mut scratch);
        assert_eq!(a, 2721565146873327973);
        assert_eq!(b, hash("AA", false).min(hash("TT", false)));
    }

    mod props {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn canonical_ignores_strand(seq in "[ACGTNacgt]{0,80}") {
                let mut rc = Vec::new();
                revcomp_into(seq.as_bytes(), Alphabet::Dna, false, &mut rc).unwrap();
                let rc = String::from_utf8(rc).unwrap();
                prop_assert_eq!(hash(&seq, true), hash(&rc, true));
                prop_assert_eq!(hash(&seq, true), hash(&seq, false).min(hash(&rc, false)));
            }
        }
    }
}
