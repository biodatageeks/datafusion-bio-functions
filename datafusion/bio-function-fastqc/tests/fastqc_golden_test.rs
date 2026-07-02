//! End-to-end QC parity against committed FastQC 0.12.1 goldens, driving the
//! accumulators directly (no polars-bio / Python needed).
//!
//! Goldens are real FastQC `--nogroup` output. We assert bit-for-bit agreement
//! on per_base_quality, per_seq_gc and dup_levels, and n_seq / floor(%GC) on
//! basic_stats.

use std::collections::HashMap;
use std::sync::Arc;

use datafusion::arrow::array::{Array, Float64Array, Int32Array, RecordBatch, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion_bio_function_fastqc::ModuleSet;

fn data(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

/// Parse a 4-line-per-record FASTQ into (sequence, quality) vectors.
fn read_fastq(path: &str) -> (Vec<String>, Vec<String>) {
    let text = std::fs::read_to_string(path).expect("read fastq");
    let lines: Vec<&str> = text.lines().collect();
    let mut seqs = Vec::new();
    let mut quals = Vec::new();
    let mut i = 0;
    while i + 3 < lines.len() {
        seqs.push(lines[i + 1].to_string());
        quals.push(lines[i + 3].to_string());
        i += 4;
    }
    (seqs, quals)
}

/// Run all modules over a fixture and return the tidy result as a lookup keyed
/// by (module, label, position, metric) -> value. Null label -> "", null
/// position -> i32::MIN.
fn run_tidy(path: &str) -> HashMap<(String, String, i32, String), f64> {
    let (seqs, quals) = read_fastq(path);
    let schema = Arc::new(Schema::new(vec![
        Field::new("sequence", DataType::Utf8, true),
        Field::new("quality_scores", DataType::Utf8, true),
    ]));
    let batch = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(StringArray::from_iter_values(seqs.iter())),
            Arc::new(StringArray::from_iter_values(quals.iter())),
        ],
    )
    .unwrap();

    let mut set = ModuleSet::build(None).unwrap();
    set.update_batch(&batch).unwrap();
    let out = set.finalize().unwrap();

    let module = out
        .column_by_name("module")
        .unwrap()
        .as_any()
        .downcast_ref::<StringArray>()
        .unwrap();
    let label = out
        .column_by_name("label")
        .unwrap()
        .as_any()
        .downcast_ref::<StringArray>()
        .unwrap();
    let position = out
        .column_by_name("position")
        .unwrap()
        .as_any()
        .downcast_ref::<Int32Array>()
        .unwrap();
    let metric = out
        .column_by_name("metric")
        .unwrap()
        .as_any()
        .downcast_ref::<StringArray>()
        .unwrap();
    let value = out
        .column_by_name("value")
        .unwrap()
        .as_any()
        .downcast_ref::<Float64Array>()
        .unwrap();

    let mut map = HashMap::new();
    for i in 0..out.num_rows() {
        if value.is_null(i) {
            continue; // status rows carry only value_str
        }
        let key = (
            module.value(i).to_string(),
            if label.is_null(i) {
                String::new()
            } else {
                label.value(i).to_string()
            },
            if position.is_null(i) {
                i32::MIN
            } else {
                position.value(i)
            },
            metric.value(i).to_string(),
        );
        map.insert(key, value.value(i));
    }
    map
}

/// A tolerant lookup of a golden fastqc_data.txt for one module. Returns rows
/// as (label, position, metric) -> value using the same conventions as `run_tidy`.
fn golden(path: &str) -> GoldenData {
    let text = std::fs::read_to_string(path).expect("read golden");
    let mut per_base: Vec<[f64; 6]> = Vec::new(); // position-1 -> [mean,median,q1,q3,p10,p90]
    let mut per_seq_gc: HashMap<i32, f64> = HashMap::new();
    let mut dup_pct: HashMap<String, f64> = HashMap::new();
    let mut total_dedup = f64::NAN;
    let mut n_seq = f64::NAN;
    let mut gc_pct = f64::NAN;

    let mut per_seq_q: HashMap<i32, f64> = HashMap::new();
    let mut seq_len: HashMap<i32, f64> = HashMap::new();
    let mut per_base_n: HashMap<i32, f64> = HashMap::new();
    let mut per_base_content: HashMap<(i32, String), f64> = HashMap::new();
    let mut overrep: HashMap<String, (f64, f64)> = HashMap::new(); // seq -> (count, pct)
    let mut adapter: HashMap<(i32, String), f64> = HashMap::new(); // (pos, adapter) -> pct

    // Column order of the Adapter Content table (filled from its header line).
    let mut adapter_names: Vec<String> = Vec::new();

    let mut module = "";
    for line in text.lines() {
        if line.starts_with(">>END_MODULE") {
            module = "";
            continue;
        }
        if let Some(rest) = line.strip_prefix(">>") {
            module = match rest.rsplit_once('\t').map(|(m, _)| m).unwrap_or(rest).trim() {
                "Basic Statistics" => "basic_stats",
                "Per base sequence quality" => "per_base_quality",
                "Per sequence quality scores" => "per_seq_quality",
                "Per base sequence content" => "per_base_content",
                "Per sequence GC content" => "per_seq_gc",
                "Per base N content" => "per_base_n",
                "Sequence Length Distribution" => "seq_length",
                "Overrepresented sequences" => "overrepresented",
                "Adapter Content" => "adapter_content",
                "Sequence Duplication Levels" => "dup_levels",
                _ => "",
            };
            continue;
        }
        if module.is_empty() {
            continue;
        }
        // The Adapter Content header names the adapter columns.
        if module == "adapter_content" && line.starts_with("#Position") {
            adapter_names = line.split('\t').skip(1).map(|s| s.to_string()).collect();
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        match module {
            "basic_stats" if line.starts_with("Total Sequences") => {
                n_seq = cols[1].parse().unwrap()
            },
            "basic_stats" if line.starts_with("%GC") => gc_pct = cols[1].parse().unwrap(),
            "per_base_quality" if !line.starts_with('#') => {
                let pos: usize = cols[0].split('-').next().unwrap().parse().unwrap();
                while per_base.len() < pos {
                    per_base.push([0.0; 6]);
                }
                for j in 0..6 {
                    per_base[pos - 1][j] = cols[j + 1].parse().unwrap();
                }
            },
            "per_seq_gc" if !line.starts_with('#') => {
                per_seq_gc.insert(cols[0].parse::<f64>().unwrap() as i32, cols[1].parse().unwrap());
            },
            "per_seq_quality" if !line.starts_with('#') => {
                per_seq_q.insert(cols[0].parse().unwrap(), cols[1].parse().unwrap());
            },
            "seq_length" if !line.starts_with('#') => {
                // FastQC may print a length range "min-max"; our fixtures are uniform.
                let pos: i32 = cols[0].split('-').next().unwrap().parse().unwrap();
                seq_len.insert(pos, cols[1].parse().unwrap());
            },
            "per_base_n" if !line.starts_with('#') => {
                let pos: i32 = cols[0].split('-').next().unwrap().parse().unwrap();
                per_base_n.insert(pos, cols[1].parse().unwrap());
            },
            "per_base_content" if !line.starts_with('#') => {
                // columns: Base, %G, %A, %T, %C
                let pos: i32 = cols[0].split('-').next().unwrap().parse().unwrap();
                for (k, base) in ["G", "A", "T", "C"].iter().enumerate() {
                    per_base_content.insert((pos, base.to_string()), cols[k + 1].parse().unwrap());
                }
            },
            "overrepresented" if !line.starts_with('#') => {
                // columns: Sequence, Count, Percentage, Possible Source
                overrep.insert(
                    cols[0].to_string(),
                    (cols[1].parse().unwrap(), cols[2].parse().unwrap()),
                );
            },
            "adapter_content" if !line.starts_with('#') => {
                let pos: i32 = cols[0].split('-').next().unwrap().parse().unwrap();
                for (k, name) in adapter_names.iter().enumerate() {
                    adapter.insert((pos, name.clone()), cols[k + 1].parse().unwrap());
                }
            },
            "dup_levels" if line.starts_with("#Total Deduplicated") => {
                total_dedup = cols[1].parse().unwrap()
            },
            "dup_levels" if !line.starts_with('#') => {
                dup_pct.insert(cols[0].to_string(), cols[1].parse().unwrap());
            },
            _ => {},
        }
    }
    GoldenData {
        per_base,
        per_seq_gc,
        dup_pct,
        total_dedup,
        n_seq,
        gc_pct,
        per_seq_q,
        seq_len,
        per_base_n,
        per_base_content,
        overrep,
        adapter,
    }
}

struct GoldenData {
    per_base: Vec<[f64; 6]>,
    per_seq_gc: HashMap<i32, f64>,
    dup_pct: HashMap<String, f64>,
    total_dedup: f64,
    n_seq: f64,
    gc_pct: f64,
    per_seq_q: HashMap<i32, f64>,
    seq_len: HashMap<i32, f64>,
    per_base_n: HashMap<i32, f64>,
    per_base_content: HashMap<(i32, String), f64>,
    overrep: HashMap<String, (f64, f64)>,
    adapter: HashMap<(i32, String), f64>,
}

#[test]
fn per_base_quality_matches_fastqc_exactly() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    let metrics = ["mean", "median", "q1", "q3", "p10", "p90"];
    assert!(!g.per_base.is_empty());
    for (i, expected) in g.per_base.iter().enumerate() {
        let pos = i as i32 + 1;
        for (j, metric) in metrics.iter().enumerate() {
            let got = ours[&(
                "per_base_quality".into(),
                String::new(),
                pos,
                metric.to_string(),
            )];
            assert!(
                (got - expected[j]).abs() <= 1e-6,
                "per_base_quality pos {pos} {metric}: {got} vs {}",
                expected[j]
            );
        }
    }
}

#[test]
fn per_seq_gc_matches_fastqc_exactly() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    for (gc, expected) in &g.per_seq_gc {
        let got = ours[&("per_seq_gc".into(), String::new(), *gc, "count".into())];
        assert!(
            (got - expected).abs() <= 1e-9,
            "per_seq_gc bin {gc}: {got} vs {expected}"
        );
    }
}

#[test]
fn dup_levels_match_fastqc_exactly() {
    let ours = run_tidy(&data("dup_mix.fastq"));
    let g = golden(&data("dup_mix.nogroup.fastqc_data.txt"));
    assert_eq!(g.dup_pct.len(), 16);
    for (label, expected) in &g.dup_pct {
        let got = ours[&("dup_levels".into(), label.clone(), i32::MIN, "pct".into())];
        assert!(
            (got - expected).abs() <= 1e-9,
            "dup_levels {label}: {got} vs {expected}"
        );
    }
    let got_dedup = ours[&(
        "dup_levels".into(),
        String::new(),
        i32::MIN,
        "total_dedup_pct".into(),
    )];
    assert!((got_dedup - g.total_dedup).abs() <= 1e-9);
}

#[test]
fn per_seq_quality_matches_fastqc_exactly() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    assert!(!g.per_seq_q.is_empty());
    for (phred, expected) in &g.per_seq_q {
        let got = ours[&("per_seq_quality".into(), String::new(), *phred, "count".into())];
        assert!(
            (got - expected).abs() <= 1e-9,
            "per_seq_quality phred {phred}: {got} vs {expected}"
        );
    }
}

#[test]
fn seq_length_matches_fastqc_exactly() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    assert!(!g.seq_len.is_empty());
    for (len, expected) in &g.seq_len {
        let got = ours[&("seq_length".into(), String::new(), *len, "count".into())];
        assert!(
            (got - expected).abs() <= 1e-9,
            "seq_length {len}: {got} vs {expected}"
        );
    }
}

#[test]
fn per_base_n_matches_fastqc_exactly() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    assert!(!g.per_base_n.is_empty());
    for (pos, expected) in &g.per_base_n {
        let got = ours[&("per_base_n".into(), String::new(), *pos, "pct".into())];
        assert!(
            (got - expected).abs() <= 1e-6,
            "per_base_n pos {pos}: {got} vs {expected}"
        );
    }
}

#[test]
fn per_base_content_matches_fastqc_exactly() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    assert!(!g.per_base_content.is_empty());
    for ((pos, base), expected) in &g.per_base_content {
        let got = ours[&("per_base_content".into(), String::new(), *pos, base.clone())];
        assert!(
            (got - expected).abs() <= 1e-6,
            "per_base_content pos {pos} {base}: {got} vs {expected}"
        );
    }
}

#[test]
fn overrepresented_matches_fastqc_exactly() {
    // dup_mix has genuine overrepresented sequences (up to 4.76%).
    let ours = run_tidy(&data("dup_mix.fastq"));
    let g = golden(&data("dup_mix.nogroup.fastqc_data.txt"));
    assert!(!g.overrep.is_empty());
    for (seq, (exp_count, exp_pct)) in &g.overrep {
        let got_count = ours[&(
            "overrepresented".into(),
            seq.clone(),
            i32::MIN,
            "count".into(),
        )];
        let got_pct = ours[&("overrepresented".into(), seq.clone(), i32::MIN, "pct".into())];
        assert_eq!(got_count, *exp_count, "overrep {seq} count");
        assert!(
            (got_pct - exp_pct).abs() <= 1e-9,
            "overrep {seq} pct: {got_pct} vs {exp_pct}"
        );
    }
    // and no extra sequences beyond what FastQC reported
    let our_seqs: std::collections::HashSet<String> = ours
        .keys()
        .filter(|(m, _, _, metric)| m == "overrepresented" && metric == "count")
        .map(|(_, label, _, _)| label.clone())
        .collect();
    assert_eq!(our_seqs.len(), g.overrep.len(), "overrepresented set size");
}

#[test]
fn adapter_content_matches_fastqc_exactly() {
    // adapter_mix has 20/50 reads carrying the Illumina Universal Adapter.
    let ours = run_tidy(&data("adapter_mix.fastq"));
    let g = golden(&data("adapter_mix.nogroup.fastqc_data.txt"));
    assert!(!g.adapter.is_empty());
    let mut checked_nonzero = false;
    for ((pos, name), expected) in &g.adapter {
        let got = ours[&("adapter_content".into(), name.clone(), *pos, "pct".into())];
        assert!(
            (got - expected).abs() <= 1e-9,
            "adapter {name} pos {pos}: {got} vs {expected}"
        );
        if *expected > 0.0 {
            checked_nonzero = true;
        }
    }
    assert!(checked_nonzero, "fixture should exercise a non-zero adapter curve");
}

#[test]
fn basic_stats_match_fastqc() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    let n = ours[&("basic_stats".into(), String::new(), i32::MIN, "n_seq".into())];
    assert_eq!(n, g.n_seq);
    let gc = ours[&("basic_stats".into(), String::new(), i32::MIN, "gc_pct".into())];
    // FastQC prints floor((G+C)*100/(A+T+G+C)); our precise value must floor to it.
    assert_eq!(gc.floor(), g.gc_pct);
}
