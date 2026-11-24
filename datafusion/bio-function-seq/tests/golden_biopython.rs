//! Golden tests: the kernels against Biopython 1.88 and Trimmomatic 0.39
//! outputs in `tests/data` (see `make_golden_biopython.py` there).

use std::collections::HashMap;
use std::path::PathBuf;

use datafusion_bio_function_seq::kernels::{
    Ambiguity, KernelError, MolecularWeight, Motif, QualityTrimmer, Translator, Windows,
    back_transcribe_into, gc_skew_into, transcribe_into,
};

fn data(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Rows without the header, padded to the header's width: the
/// trailing-whitespace hook strips trailing empty fields.
fn tsv(name: &str) -> Vec<Vec<String>> {
    let text = data(name);
    let mut lines = text.lines();
    let width = lines.next().unwrap().split('\t').count();
    lines
        .map(|l| {
            let mut row: Vec<String> = l.split('\t').map(str::to_string).collect();
            assert!(row.len() <= width, "{name}: too many fields in {l:?}");
            row.resize(width, String::new());
            row
        })
        .collect()
}

/// `id -> (seq, qual)` from FASTQ, or `id -> (seq, "")` from FASTA.
fn records(name: &str) -> Vec<(String, String, String)> {
    let text = data(name);
    let lines: Vec<&str> = text.lines().collect();
    let step = if name.ends_with(".fastq") { 4 } else { 2 };
    lines
        .chunks(step)
        .map(|r| {
            let qual = r.get(3).copied().unwrap_or("");
            (r[0][1..].to_string(), r[1].to_string(), qual.to_string())
        })
        .collect()
}

fn list<T: std::str::FromStr>(field: &str) -> Vec<T>
where
    T::Err: std::fmt::Debug,
{
    field
        .split(',')
        .filter(|v| !v.is_empty())
        .map(|v| v.parse().unwrap())
        .collect()
}

#[test]
fn cds_translation_and_status() {
    // Seq.translate(table=t, cds=True): the protein, or `!<failed check>`.
    let rows = tsv("bio_cds.tsv");
    assert!(rows.len() > 2000);
    let (mut proteins, mut failures) = (0, 0);
    for row in &rows {
        let (table, seq, expected) = (row[0].parse::<i64>().unwrap(), &row[1], &row[2]);
        let mut translator = Translator::new(1, table, false)
            .unwrap()
            .with_cds(true)
            .with_ambiguity(Ambiguity::Biopython);
        let mut out = Vec::new();
        let got = match translator.translate_into(seq.as_bytes(), &mut out) {
            Ok(()) => String::from_utf8(out).unwrap(),
            Err(KernelError::InvalidCds { error, .. }) => format!("!{}", error.as_str()),
            Err(other) => panic!("table {table} {seq}: {other}"),
        };
        assert_eq!(&got, expected, "table {table} {seq}");
        let status = translator.cds_status(seq.as_bytes());
        match expected.strip_prefix('!') {
            Some(check) => {
                assert_eq!(
                    status.map(|e| e.as_str()),
                    Some(check),
                    "table {table} {seq}"
                );
                failures += 1;
            }
            None => {
                assert_eq!(status, None, "table {table} {seq}");
                proteins += 1;
            }
        }
    }
    assert!(proteins > 500 && failures > 1000, "{proteins} {failures}");
}

#[test]
fn translation_with_biopython_ambiguity_letters() {
    // Seq.translate per frame and table, B, Z and J left in place.
    let seqs: HashMap<String, String> = records("coding.fasta")
        .into_iter()
        .map(|(id, seq, _)| (id, seq))
        .collect();
    let rows = tsv("bio_translate_raw.tsv");
    assert_eq!(rows.len(), 40 * 27 * 6);
    let mut with_letters = 0;
    for row in &rows {
        let (frame, table) = (row[1].parse().unwrap(), row[2].parse().unwrap());
        let mut translator = Translator::new(frame, table, false)
            .unwrap()
            .with_ambiguity(Ambiguity::Biopython);
        let mut out = Vec::new();
        translator
            .translate_into(seqs[&row[0]].as_bytes(), &mut out)
            .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            row[3],
            "{} frame {frame} table {table}",
            row[0]
        );
        with_letters += row[3].contains(['B', 'Z', 'J']) as usize;
    }
    assert!(with_letters > 1000, "{with_letters}");
}

#[test]
fn transcription() {
    let seqs: HashMap<String, String> = records("corpus.fastq")
        .into_iter()
        .map(|(id, seq, _)| (id, seq))
        .collect();
    let rows = tsv("bio_transcribe.tsv");
    assert_eq!(rows.len(), 300);
    for row in &rows {
        let seq = seqs[&row[0]].as_bytes();
        let (mut rna, mut dna) = (Vec::new(), Vec::new());
        transcribe_into(seq, &mut rna);
        back_transcribe_into(seq, &mut dna);
        assert_eq!(String::from_utf8(rna).unwrap(), row[1], "{}", row[0]);
        assert_eq!(String::from_utf8(dna).unwrap(), row[2], "{}", row[0]);
    }
}

#[test]
fn motif_search() {
    // Bio.SeqUtils.nt_search on uppercase DNA: 0-based, overlapping hits.
    let rows = tsv("bio_motif.tsv");
    assert_eq!(rows.len(), 300);
    let (mut hits, mut out) = (0, Vec::new());
    for row in &rows {
        let motif = Motif::new(&row[1], false, true).unwrap();
        motif.find_into(row[0].as_bytes(), &mut out);
        assert_eq!(out, list::<u64>(&row[2]), "{} in {}", row[1], row[0]);
        hits += out.len();
    }
    assert!(hits > 500, "{hits}");
}

#[test]
fn molecular_weight() {
    // Bio.SeqUtils.molecular_weight; an empty weight where it raises.
    let rows = tsv("bio_weight.tsv");
    assert!(rows.len() > 800);
    let (mut rejected, mut exact) = (0, 0);
    for row in &rows {
        let flag = |i: usize| row[i] == "true";
        let weigher = MolecularWeight::new(&row[0], flag(1), flag(2), flag(3)).unwrap();
        let got = weigher.weight(row[5].as_bytes());
        if row[4].is_empty() {
            assert_eq!(got, None, "{row:?}");
            rejected += 1;
            continue;
        }
        let expected: f64 = row[4].parse().unwrap();
        let got = got.unwrap_or_else(|| panic!("no weight for {row:?}"));
        assert!(
            (got - expected).abs() <= 1e-12 * expected.abs().max(1.0),
            "{got} != {expected} for {row:?}"
        );
        exact += (got == expected) as usize;
    }
    assert!(rejected > 20, "{rejected}");
    // The summation follows Python's, so nearly every value is bit-identical.
    assert!(exact * 10 >= (rows.len() - rejected) * 9, "{exact}");
}

#[test]
fn gc_skew() {
    let seqs: HashMap<String, String> = records("corpus.fastq")
        .into_iter()
        .map(|(id, seq, _)| (id, seq))
        .collect();
    let rows = tsv("bio_gc_skew.tsv");
    assert_eq!(rows.len(), 300 * 4);
    let mut out = Vec::new();
    for row in &rows {
        out.clear();
        let windows = Windows::new(row[1].parse().unwrap(), None).unwrap();
        gc_skew_into(seqs[&row[0]].as_bytes(), windows, &mut out);
        assert_eq!(out, list::<f64>(&row[2]), "{} window {}", row[0], row[1]);
    }
}

#[test]
fn quality_trimming() {
    // Trimmomatic 0.39 `SE -phred33` with LEADING, TRAILING and SLIDINGWINDOW
    // in that order. A read missing from the file was dropped.
    let reads = records("corpus.fastq");
    type Config = (String, String, String, String);
    let mut kept: HashMap<Config, HashMap<String, (String, String)>> = HashMap::new();
    for row in tsv("trimmomatic.tsv") {
        let config = (
            row[0].clone(),
            row[1].clone(),
            row[2].clone(),
            row[3].clone(),
        );
        kept.entry(config)
            .or_default()
            .insert(row[4].clone(), (row[5].clone(), row[6].clone()));
    }
    assert_eq!(kept.len(), 9);
    let (mut trimmed, mut dropped) = (0, 0);
    for ((leading, trailing, size, quality), expected) in &kept {
        let int = |s: &String| (!s.is_empty()).then(|| s.parse::<i64>().unwrap());
        let window = int(size).map(|size| (size, quality.parse::<f64>().unwrap()));
        let trimmer = QualityTrimmer::new(int(leading), int(trailing), window, 33).unwrap();
        for (id, seq, qual) in &reads {
            let (start, end) = trimmer.span(seq.as_bytes(), qual.as_bytes()).unwrap();
            let got = (start != end).then(|| (&seq[start..end], &qual[start..end]));
            let want = expected.get(id).map(|(s, q)| (s.as_str(), q.as_str()));
            assert_eq!(
                got, want,
                "{id} with leading={leading} trailing={trailing} window={size}:{quality}"
            );
            match got {
                Some((s, _)) if s.len() < seq.len() => trimmed += 1,
                None => dropped += 1,
                _ => {}
            }
        }
    }
    assert!(trimmed > 500 && dropped > 500, "{trimmed} {dropped}");
}
