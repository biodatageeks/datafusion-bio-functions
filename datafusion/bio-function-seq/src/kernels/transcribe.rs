//! Transcription and back-transcription.

/// Appends `seq` with `T` replaced by `U` and `t` by `u`, as Biopython's
/// `Seq.transcribe` does. Every other byte is copied unchanged.
pub fn transcribe_into(seq: &[u8], out: &mut Vec<u8>) {
    out.extend(seq.iter().map(|&b| match b {
        b'T' => b'U',
        b't' => b'u',
        b => b,
    }));
}

/// Appends `seq` with `U` replaced by `T` and `u` by `t`, as Biopython's
/// `Seq.back_transcribe` does. Every other byte is copied unchanged.
pub fn back_transcribe_into(seq: &[u8], out: &mut Vec<u8>) {
    out.extend(seq.iter().map(|&b| match b {
        b'U' => b'T',
        b'u' => b't',
        b => b,
    }));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(f: fn(&[u8], &mut Vec<u8>), seq: &str) -> String {
        let mut out = Vec::new();
        f(seq.as_bytes(), &mut out);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn biopython_probes() {
        // Seq.transcribe() / Seq.back_transcribe(), Biopython 1.88.
        assert_eq!(run(transcribe_into, "ACGTacgtN-U"), "ACGUacguN-U");
        assert_eq!(run(back_transcribe_into, "ACGTacgtN-U"), "ACGTacgtN-T");
        assert_eq!(run(transcribe_into, "ACGUacgu"), "ACGUacgu");
        assert_eq!(run(back_transcribe_into, "ACGUacgu"), "ACGTacgt");
    }

    #[test]
    fn empty_and_appending() {
        assert_eq!(run(transcribe_into, ""), "");
        let mut out = b"x".to_vec();
        transcribe_into(b"T", &mut out);
        assert_eq!(out, b"xU");
    }

    #[test]
    fn multibyte_characters_are_kept_whole() {
        // 'é' is two bytes, neither of which is T/t/U/u.
        assert_eq!(run(transcribe_into, "TéT"), "UéU");
        assert_eq!(run(back_transcribe_into, "UéU"), "TéT");
    }

    mod props {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn round_trip_without_uracil(seq in "[ACGTNacgtn-]{0,80}") {
                let rna = run(transcribe_into, &seq);
                prop_assert!(!rna.contains('T') && !rna.contains('t'));
                prop_assert_eq!(run(back_transcribe_into, &rna), seq);
            }
        }
    }
}
