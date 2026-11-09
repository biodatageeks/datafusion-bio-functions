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

fn segment(
    genomic_start: i64,
    genomic_end: i64,
    cdna_start: usize,
    ori: i8,
) -> TranscriptCdnaMapperSegment {
    TranscriptCdnaMapperSegment {
        genomic_start,
        genomic_end,
        cdna_start,
        cdna_end: cdna_start + usize::try_from((genomic_end - genomic_start).max(0)).unwrap(),
        ori,
    }
}

/// Generated segment lists of the shapes the memo and bisection must handle:
/// disjoint and sorted (the shape of every real transcript), plus contiguous
/// runs with cDNA gaps, overlaps, unsorted order, malformed spans and empty.
fn generated_segments(rng: &mut Lcg, count: usize) -> Vec<TranscriptCdnaMapperSegment> {
    let shape = rng.below(5);
    let ori = if rng.below(2) == 0 { 1 } else { -1 };
    let mut segments = Vec::with_capacity(count);
    let mut next_start = 100i64;
    let mut cdna = 1usize;
    for _ in 0..count {
        let len = 1 + rng.below(6) as i64;
        let start = match shape {
            // Overlaps the previous segment.
            1 => next_start - rng.below(4) as i64,
            // Contiguous with the previous segment.
            2 => next_start,
            _ => next_start + rng.below(5) as i64,
        };
        // Malformed spans (end before start) occur only in the last shape.
        let end = if shape == 4 && rng.below(3) == 0 {
            start - 1
        } else {
            start + len - 1
        };
        // A contiguous pair with a cDNA jump rejects the general mapper.
        let jump = if shape == 2 && rng.below(3) == 0 {
            2
        } else {
            0
        };
        segments.push(segment(start, end, cdna + jump, ori));
        next_start = end + 1;
        cdna += len as usize + jump;
    }
    if shape == 3 {
        for i in (1..segments.len()).rev() {
            segments.swap(i, rng.below(i + 1));
        }
    }
    segments
}

/// The pre-change body of `sorted_mapper_segments`.
fn sorted_segments_reference(
    segments: &[TranscriptCdnaMapperSegment],
) -> Vec<TranscriptCdnaMapperSegment> {
    let mut owned = segments.to_vec();
    owned.sort_by_key(|s| (s.genomic_start, s.genomic_end, s.cdna_start));
    owned
}

/// The pre-change body of `use_cdna_mapper_for_general_coords`.
fn general_coords_reference(segments: &[TranscriptCdnaMapperSegment]) -> bool {
    if segments.is_empty() {
        return false;
    }
    for window in sorted_segments_reference(segments).windows(2) {
        let (prev, next) = (&window[0], &window[1]);
        if next.genomic_start == prev.genomic_end.saturating_add(1)
            && next.cdna_start != prev.cdna_end.saturating_add(1)
        {
            return false;
        }
    }
    true
}

#[test]
fn memoized_mapper_facts_match_direct_computation() {
    let mut rng = Lcg(4);
    let (mut general, mut rejected, mut unsorted) = (0, 0, 0);
    for _ in 0..30_000 {
        let mut transcript = tx("T1", "1", 1, 1000, 1, "protein_coding", Some(1), Some(1000));
        let count = rng.below(9);
        transcript.cdna_mapper_segments = generated_segments(&mut rng, count);
        let expected = general_coords_reference(&transcript.cdna_mapper_segments);
        let sorted = sorted_segments_reference(&transcript.cdna_mapper_segments);
        for _ in 0..2 {
            assert_eq!(use_cdna_mapper_for_general_coords(&transcript), expected);
            assert_eq!(sorted_mapper_segments(&transcript).as_ref(), sorted);
        }
        general += usize::from(expected);
        rejected += usize::from(!expected && count > 0);
        unsorted += usize::from(sorted != transcript.cdna_mapper_segments);
    }
    assert!(general > 1_000 && rejected > 1_000 && unsorted > 1_000);
}

#[test]
fn mapper_facts_do_not_outlive_their_segments() {
    let mut transcript = tx("T1", "1", 1, 1000, 1, "protein_coding", Some(1), Some(1000));
    // A contiguous pair with a cDNA jump rejects the general mapper.
    transcript.cdna_mapper_segments = vec![segment(100, 109, 1, 1), segment(110, 119, 20, 1)];
    assert!(!use_cdna_mapper_for_general_coords(&transcript));

    let mut clone = transcript.clone();
    clone.cdna_mapper_segments = vec![segment(100, 109, 1, 1), segment(110, 119, 11, 1)];
    assert!(use_cdna_mapper_for_general_coords(&clone));

    // Replacing or clearing the segments on a transcript whose memo is filled.
    transcript.cdna_mapper_segments = vec![segment(100, 109, 1, 1), segment(110, 119, 11, 1)];
    assert!(use_cdna_mapper_for_general_coords(&transcript));
    assert_eq!(
        exact_mapper_cdna_index_for_transcript(&transcript, 115),
        Some(16)
    );
    transcript.cdna_mapper_segments.clear();
    assert!(!use_cdna_mapper_for_general_coords(&transcript));
    assert_eq!(
        exact_mapper_cdna_index_for_transcript(&transcript, 115),
        None
    );
}

#[test]
fn bisected_mapper_lookup_matches_first_match_scan() {
    let mut rng = Lcg(5);
    let mut bisected = 0;
    for _ in 0..20_000 {
        let mut transcript = tx("T1", "1", 1, 1000, 1, "protein_coding", Some(1), Some(1000));
        let count = rng.below(12);
        transcript.cdna_mapper_segments = generated_segments(&mut rng, count);
        bisected += usize::from(count > 1 && mapper_segment_facts(&transcript).disjoint);
        for pos in 90..=220 {
            assert_eq!(
                exact_mapper_cdna_index_for_transcript(&transcript, pos),
                transcript
                    .cdna_mapper_segments
                    .iter()
                    .find_map(|segment| mapper_segment_cdna_index(segment, pos)),
                "pos {pos} segments {:?}",
                transcript.cdna_mapper_segments
            );
        }
    }
    assert!(bisected > 1_000, "too few disjoint lists: {bisected}");
}

#[test]
fn bisected_coord_table_matches_linear_scan() {
    let mut rng = Lcg(6);
    let mut bisected = 0;
    for _ in 0..20_000 {
        let count = rng.below(12);
        let coords: Vec<TranscriptCdnaCoord> = generated_segments(&mut rng, count)
            .into_iter()
            .map(|s| TranscriptCdnaCoord {
                start: s.genomic_start,
                end: s.genomic_end,
                cdna_start: s.cdna_start,
                cdna_end: s.cdna_end,
            })
            .collect();
        let table = CdnaCoordTable::new(coords.clone());
        bisected += usize::from(count > 1 && table.disjoint);
        for strand in [1i8, -1] {
            for pos in 90..=220 {
                assert_eq!(
                    table.cdna_index(strand, pos),
                    cdna_index_from_sorted_coords(strand, pos, &coords),
                    "strand {strand} pos {pos} coords {coords:?}"
                );
            }
        }
    }
    assert!(bisected > 1_000, "too few disjoint tables: {bisected}");
}

#[test]
fn alt_translation_text_matches_collected_residues() {
    let transcript = tx("T1", "1", 1, 100, 1, "protein_coding", Some(1), Some(100));
    let class = CodingClassification {
        protein_position_start: Some(1),
        protein_position_end: Some(1),
        ..Default::default()
    };
    let old_aas = ['M', 'A', '*'];
    let new_aas = ['K', 'L', 'V', '*'];
    let text: String = new_aas.iter().collect();
    for semantics in SEMANTICS {
        let build = |new_aas: &[char], alt_text: Option<&str>| {
            build_protein_hgvs_data(
                &transcript,
                &class,
                Some("M/K"),
                &old_aas,
                new_aas,
                alt_text,
                true,
                None,
                semantics,
                None,
            )
        };
        assert_eq!(build(&new_aas, None), build(&[], Some(text.as_str())));
    }
}
