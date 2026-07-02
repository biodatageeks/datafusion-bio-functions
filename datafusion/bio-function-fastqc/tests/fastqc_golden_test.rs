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
                "Per sequence GC content" => "per_seq_gc",
                "Sequence Duplication Levels" => "dup_levels",
                _ => "",
            };
            continue;
        }
        if module.is_empty() {
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
    }
}

struct GoldenData {
    per_base: Vec<[f64; 6]>,
    per_seq_gc: HashMap<i32, f64>,
    dup_pct: HashMap<String, f64>,
    total_dedup: f64,
    n_seq: f64,
    gc_pct: f64,
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
fn basic_stats_match_fastqc() {
    let ours = run_tidy(&data("example.fastq"));
    let g = golden(&data("example.nogroup.fastqc_data.txt"));
    let n = ours[&("basic_stats".into(), String::new(), i32::MIN, "n_seq".into())];
    assert_eq!(n, g.n_seq);
    let gc = ours[&("basic_stats".into(), String::new(), i32::MIN, "gc_pct".into())];
    // FastQC prints floor((G+C)*100/(A+T+G+C)); our precise value must floor to it.
    assert_eq!(gc.floor(), g.gc_pct);
}
