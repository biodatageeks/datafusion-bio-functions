//! The allocation-free helpers must return exactly what the copying code they
//! replaced returned. Each test keeps that code as a reference implementation
//! and compares both over generated inputs, including empty and incomplete
//! sequences, out-of-range windows, and non-ACGT symbols.
use super::*;

/// Small deterministic generator (no extra test dependency).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn bytes(&mut self, alphabet: &[u8], len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| alphabet[self.below(alphabet.len())])
            .collect()
    }
}

const SEMANTICS: [crate::vep_semantics::VepSemantics; 2] = [
    crate::vep_semantics::VepSemantics::V115,
    crate::vep_semantics::VepSemantics::V116,
];

#[test]
fn string_translation_matches_collected_vec_translation() {
    let mut rng = Lcg(1);
    let alphabets: [&[u8]; 3] = [b"ACGT", b"ACGTN", b"ACGTNacgtnRYKM*-X"];
    for table in [CodonTable::Standard, CodonTable::VertebrateMitochondrial] {
        for alphabet in alphabets {
            for len in 0..40 {
                for _ in 0..40 {
                    let cds = rng.bytes(alphabet, len);
                    let expected = translate_protein_from_cds_in(table, &cds)
                        .map(|aas| aas.into_iter().collect::<String>());
                    assert_eq!(translate_protein_string_in(table, &cds), expected);
                    if table == CodonTable::Standard {
                        assert_eq!(translate_protein_string(&cds), expected);
                    }
                }
            }
        }
    }
    // Invalid UTF-8 inside a non-ACGT codon fails both translators.
    let invalid = [b'A', b'C', b'G', 0xff, b'A', b'C'];
    assert_eq!(translate_protein_from_cds(&invalid), None);
    assert_eq!(translate_protein_string(&invalid), None);
}

#[test]
fn trimmed_cds_len_matches_trim_incomplete_codon() {
    for len in 0..16 {
        let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[i % 4]).collect();
        let expected = vep116_trim_incomplete_codon(&seq);
        let actual = vep116_trimmed_cds_len(&seq).map(|keep| seq[..keep].to_vec());
        assert_eq!(actual, expected, "len {len}");
    }
}

/// The pre-change body of `vep_alternate_codon_window`.
fn alternate_codon_window_reference(
    tx: &TranscriptFeature,
    mutated: &[u8],
    codon_start: usize,
    window_len: usize,
    semantics: crate::vep_semantics::VepSemantics,
) -> Vec<u8> {
    let mut alternate_cds = mutated.to_vec();
    if semantics.vep116_stop_predicates() {
        let Some(trimmed) = vep116_trim_incomplete_codon(&alternate_cds) else {
            return Vec::new();
        };
        alternate_cds = trimmed;
    }
    let utr = three_prime_utr_seq(tx).unwrap_or_default();
    alternate_cds.extend_from_slice(utr.as_bytes());
    alternate_cds
        .iter()
        .skip(codon_start)
        .take(window_len)
        .copied()
        .collect()
}

#[test]
fn alternate_codon_window_matches_concatenated_cds_and_utr() {
    let mut rng = Lcg(2);
    for biotype in ["protein_coding", "protein_coding_LoF"] {
        for utr in [None, Some(""), Some("acgTAG"), Some("GGGCCCAAATTT")] {
            let mut transcript = tx("T1", "1", 1, 100, 1, biotype, Some(1), Some(100));
            transcript.three_prime_utr_seq = utr.map(str::to_string);
            for semantics in SEMANTICS {
                for len in 0..14 {
                    let mutated = rng.bytes(b"ACGT", len);
                    for codon_start in 0..len + 16 {
                        for window_len in [0, 1, 2, 3, 4, 5, 6, 9, 17, usize::MAX] {
                            assert_eq!(
                                vep_alternate_codon_window(
                                    &transcript,
                                    &mutated,
                                    codon_start,
                                    window_len,
                                    semantics,
                                ),
                                alternate_codon_window_reference(
                                    &transcript,
                                    &mutated,
                                    codon_start,
                                    window_len,
                                    semantics,
                                ),
                                "{biotype} utr={utr:?} {semantics:?} len={len} \
                                 start={codon_start} window={window_len}"
                            );
                        }
                    }
                }
            }
        }
    }
}

/// The pre-change body of `vep116_ref_eq_alt_sequence_final_stop`.
fn final_stop_reference(
    ref_translation: &[char],
    translation_start: usize,
    translation_end: usize,
    ref_peptide: &str,
    alt_peptide: &str,
) -> bool {
    if ref_peptide == "X" && alt_peptide == "X" {
        return false;
    }
    if translation_start >= ref_translation.len() && alt_peptide.starts_with('*') {
        return true;
    }
    let same_stop_index = ref_peptide
        .find('*')
        .is_some_and(|stop_index| alt_peptide.find('*') == Some(stop_index));
    let mut mutated = ref_translation.to_vec();
    let splice_start = translation_start.min(mutated.len());
    let splice_len = translation_end
        .checked_sub(translation_start)
        .and_then(|span| span.checked_add(1))
        .unwrap_or(0);
    let splice_end = translation_start
        .saturating_add(splice_len)
        .min(mutated.len())
        .max(splice_start);
    let replacement = if alt_peptide == "-" {
        Vec::new()
    } else {
        alt_peptide.chars().collect()
    };
    mutated.splice(splice_start..splice_end, replacement);
    same_stop_index
        || (mutated.starts_with(ref_translation)
            && mutated
                .get(ref_translation.len())
                .is_some_and(|residue| *residue == '*'))
}

#[test]
fn final_stop_check_matches_built_spliced_peptide() {
    let mut rng = Lcg(3);
    // A two-letter alphabet plus stop makes repeated (duplication-like)
    // insertions, and therefore `true` answers, common.
    let residues = ['A', 'B', '*'];
    // A random peptide of up to `max_len - 1` residues.
    let peptide = |rng: &mut Lcg, max_len: usize| -> String {
        let len = rng.below(max_len);
        (0..len)
            .map(|_| residues[rng.below(residues.len())])
            .collect()
    };
    let mut trues = 0;
    for _ in 0..200_000 {
        let reference: Vec<char> = peptide(&mut rng, 9).chars().collect();
        let start = rng.below(reference.len() + 3);
        let end = rng.below(reference.len() + 3);
        let ref_peptide = match rng.below(4) {
            0 => "-".to_string(),
            1 => "X".to_string(),
            _ => peptide(&mut rng, 4),
        };
        let alt_peptide = match rng.below(5) {
            0 => "-".to_string(),
            1 => "X".to_string(),
            _ => peptide(&mut rng, 7),
        };
        let expected = final_stop_reference(&reference, start, end, &ref_peptide, &alt_peptide);
        trues += usize::from(expected);
        assert_eq!(
            vep116_ref_eq_alt_sequence_final_stop(
                &reference,
                start,
                end,
                &ref_peptide,
                &alt_peptide
            ),
            expected,
            "ref={reference:?} start={start} end={end} ref_pep={ref_peptide} alt_pep={alt_peptide}"
        );
    }
    assert!(trues > 1_000, "too few positive cases exercised: {trues}");
}

#[test]
fn chunked_lowercase_scan_matches_byte_any() {
    for len in [0, 1, 63, 64, 65, 127, 128, 129, 300] {
        let upper = vec![b'A'; len];
        assert!(!has_ascii_lowercase(&upper));
        for pos in 0..len {
            for symbol in [b'a', b'z', b'n', b'Z', b'*', 0xc3] {
                let mut seq = upper.clone();
                seq[pos] = symbol;
                assert_eq!(
                    has_ascii_lowercase(&seq),
                    seq.iter().any(u8::is_ascii_lowercase),
                    "len {len} pos {pos} symbol {symbol}"
                );
            }
        }
    }
}

#[test]
fn hgvsp_reference_matches_uppercased_copy_and_borrows_uppercase_text() {
    let transcript = tx("T1", "1", 1, 100, 1, "protein_coding", Some(1), Some(100));
    for cds in ["ATGAAATAA", "atgAAAtaa", "ATGNNNTAA", "ATGaaNTAA"] {
        let tl = translation("T1", Some(9), Some(2), Some("MK"), Some(cds));
        let seq = reference_translateable_seq_for_hgvsp(&transcript, Some(&tl)).unwrap();
        assert_eq!(seq.as_ref(), cds.to_ascii_uppercase());
        assert_eq!(
            matches!(seq, Cow::Borrowed(_)),
            !cds.bytes().any(|b| b.is_ascii_lowercase()),
            "{cds}"
        );
    }
}

#[test]
fn lazy_canonical_reference_skips_fallback_when_canonical_text_exists() {
    let transcript = tx("T1", "1", 1, 100, 1, "protein_coding", Some(1), Some(100));
    let tl = translation("T1", Some(9), Some(2), Some("MK"), Some("ATGAAATAA"));
    let mut built = false;
    let canonical = canonical_ref_translation_for_hgvsp_with(&transcript, Some(&tl), || {
        built = true;
        "fallback".to_string()
    });
    assert_eq!(canonical, "MK");
    assert!(!built, "fallback built although canonical text exists");
    assert_eq!(
        canonical,
        canonical_ref_translation_for_hgvsp(&transcript, Some(&tl), "fallback".to_string())
    );
    let without = translation("T1", Some(9), None, None, Some("ATGAAATAA"));
    assert_eq!(
        canonical_ref_translation_for_hgvsp_with(&transcript, Some(&without), || {
            "fallback".to_string()
        }),
        "fallback"
    );
}
