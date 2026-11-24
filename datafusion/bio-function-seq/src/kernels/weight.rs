//! Molecular weight, after Biopython's `Bio.SeqUtils.molecular_weight`.

use super::KernelError;
use super::weights_data::{
    DNA_WEIGHTS, DNA_WEIGHTS_MONOISOTOPIC, PROTEIN_WEIGHTS, PROTEIN_WEIGHTS_MONOISOTOPIC,
    RNA_WEIGHTS, RNA_WEIGHTS_MONOISOTOPIC, WATER, WATER_MONOISOTOPIC,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SeqType {
    Dna,
    Rna,
    Protein,
}

/// Computes molecular weights the way Biopython's
/// `Bio.SeqUtils.molecular_weight` does: the sum of the residue weights minus
/// one water per bond, for a 5'-phosphorylated nucleic acid or a protein.
#[derive(Debug, Clone)]
pub struct MolecularWeight {
    /// Weight per byte, in either case; NaN for a byte Biopython rejects.
    weights: Box<[f64; 256]>,
    /// Weight of each byte's complement, when the molecule is double-stranded.
    complement_weights: Option<Box<[f64; 256]>>,
    water: f64,
    circular: bool,
}

impl MolecularWeight {
    /// `seq_type` is `"DNA"`, `"RNA"` or `"protein"`, ignoring case.
    /// `double_stranded` adds the complementary strand and is not allowed
    /// for proteins. `circular` removes the water of the closing bond.
    /// `monoisotopic` uses monoisotopic instead of average masses.
    pub fn new(
        seq_type: &str,
        double_stranded: bool,
        circular: bool,
        monoisotopic: bool,
    ) -> Result<Self, KernelError> {
        let seq_type = match seq_type.to_ascii_lowercase().as_str() {
            "dna" => SeqType::Dna,
            "rna" => SeqType::Rna,
            "protein" => SeqType::Protein,
            _ => {
                return Err(KernelError::InvalidOption {
                    name: "seq_type",
                    value: seq_type.to_string(),
                    expected: "'DNA', 'RNA' or 'protein'",
                });
            }
        };
        if double_stranded && seq_type == SeqType::Protein {
            return Err(KernelError::InvalidOption {
                name: "double_stranded",
                value: "true".to_string(),
                expected: "false for a protein",
            });
        }
        let data = match (seq_type, monoisotopic) {
            (SeqType::Dna, false) => DNA_WEIGHTS,
            (SeqType::Dna, true) => DNA_WEIGHTS_MONOISOTOPIC,
            (SeqType::Rna, false) => RNA_WEIGHTS,
            (SeqType::Rna, true) => RNA_WEIGHTS_MONOISOTOPIC,
            (SeqType::Protein, false) => PROTEIN_WEIGHTS,
            (SeqType::Protein, true) => PROTEIN_WEIGHTS_MONOISOTOPIC,
        };
        let mut weights = Box::new([f64::NAN; 256]);
        for &(letter, weight) in data {
            weights[letter as usize] = weight;
            weights[letter.to_ascii_lowercase() as usize] = weight;
        }
        let complement_weights = double_stranded.then(|| {
            let a_pair = if seq_type == SeqType::Rna { b'U' } else { b'T' };
            let mut table = Box::new([f64::NAN; 256]);
            for (base, pair) in [(b'A', a_pair), (a_pair, b'A'), (b'C', b'G'), (b'G', b'C')] {
                table[base as usize] = weights[pair as usize];
                table[base.to_ascii_lowercase() as usize] = weights[pair as usize];
            }
            table
        });
        Ok(Self {
            weights,
            complement_weights,
            water: if monoisotopic {
                WATER_MONOISOTOPIC
            } else {
                WATER
            },
            circular,
        })
    }

    /// The molecular weight of `seq`, or `None` if it has a letter Biopython
    /// rejects (an ambiguous or unknown one). Whitespace is skipped. An
    /// empty sequence gives the formula's value, one water for a linear
    /// single strand, as Biopython does.
    pub fn weight(&self, seq: &[u8]) -> Option<f64> {
        let (sum, len) = strand_sum(&self.weights, seq)?;
        let bonds = (len - 1) as f64;
        let mut weight = sum - bonds * self.water;
        if self.circular {
            weight -= self.water;
        }
        if let Some(complement) = &self.complement_weights {
            let (sum, _) = strand_sum(complement, seq)?;
            weight += sum - bonds * self.water;
            if self.circular {
                weight -= self.water;
            }
        }
        Some(weight)
    }
}

/// The summed weight and the number of letters, or `None` for a letter that
/// has no weight. Uses Neumaier's compensated summation, as Python's `sum`
/// does for floats (3.12 and later), so results agree with Biopython.
fn strand_sum(table: &[f64; 256], seq: &[u8]) -> Option<(f64, i64)> {
    let (mut sum, mut compensation, mut len) = (0.0f64, 0.0f64, 0i64);
    for &b in seq {
        // The ASCII characters Python's str.split() treats as whitespace.
        if matches!(b, b'\t'..=b'\r' | 0x1c..=0x1f | b' ') {
            continue;
        }
        let x = table[b as usize];
        if x.is_nan() {
            return None;
        }
        let t = sum + x;
        compensation += if sum.abs() >= x.abs() {
            (sum - t) + x
        } else {
            (x - t) + sum
        };
        sum = t;
        len += 1;
    }
    Some((sum + compensation, len))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mw(seq: &str, seq_type: &str, ds: bool, circular: bool, mono: bool) -> Option<f64> {
        MolecularWeight::new(seq_type, ds, circular, mono)
            .unwrap()
            .weight(seq.as_bytes())
    }

    #[track_caller]
    fn assert_close(got: Option<f64>, expected: f64) {
        let got = got.expect("a weight");
        assert!(
            (got - expected).abs() <= 1e-12 * expected.abs().max(1.0),
            "{got} != {expected}"
        );
    }

    #[test]
    fn biopython_probes() {
        // Bio.SeqUtils.molecular_weight, Biopython 1.88.
        assert_close(mw("ACGT", "DNA", false, false, false), 1253.8027);
        assert_close(mw("ACGT", "DNA", true, false, false), 2507.6054);
        assert_close(mw("ACGT", "DNA", false, true, false), 1235.7874);
        assert_close(mw("ACGT", "DNA", false, false, true), 1253.213097);
        assert_close(mw("ACGU", "RNA", false, false, false), 1303.7737);
        assert_close(
            mw("MKV", "protein", false, false, false),
            376.51460000000003,
        );
        assert_close(mw("acgt", "DNA", false, false, false), 1253.8027);
        // Biopython reads these four letters as amino acids.
        assert_close(
            mw("ACGT", "protein", false, false, false),
            350.39129999999994,
        );
        // The formula gives one water for an empty sequence.
        assert_close(mw("", "DNA", false, false, false), 18.0153);
    }

    #[test]
    fn letters_biopython_rejects_give_none() {
        assert_eq!(mw("ACGN", "DNA", false, false, false), None);
        assert_eq!(mw("ACGU", "DNA", false, false, false), None);
        assert_eq!(mw("ACGT", "RNA", false, false, false), None);
        assert_eq!(mw("MKV*", "protein", false, false, false), None);
        assert_eq!(mw("MKX", "protein", false, false, false), None);
        assert_eq!(mw("AC-GT", "DNA", false, false, false), None);
        assert_eq!(mw("ACé", "DNA", false, false, false), None);
    }

    #[test]
    fn whitespace_is_ignored() {
        // Biopython removes whitespace before weighing.
        assert_eq!(
            mw("AC GT\n", "DNA", false, false, false),
            mw("ACGT", "DNA", false, false, false)
        );
    }

    #[test]
    fn options() {
        assert!(MolecularWeight::new("dna", false, false, false).is_ok());
        assert!(MolecularWeight::new("Protein", false, true, true).is_ok());
        assert!(matches!(
            MolecularWeight::new("xna", false, false, false).unwrap_err(),
            KernelError::InvalidOption {
                name: "seq_type",
                ..
            }
        ));
        assert!(matches!(
            MolecularWeight::new("protein", true, false, false).unwrap_err(),
            KernelError::InvalidOption {
                name: "double_stranded",
                ..
            }
        ));
    }
}
