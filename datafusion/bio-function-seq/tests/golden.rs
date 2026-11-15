//! Golden tests: the `seq_*` SQL functions against seqkit 2.14.0 and
//! Biopython 1.88 outputs in `tests/data` (see `make_golden.py` there).

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use datafusion::arrow::array::{ArrayRef, AsArray, RecordBatch, StringArray};
use datafusion::arrow::compute::cast;
use datafusion::arrow::datatypes::{DataType, Field, Float64Type, Schema};
use datafusion::datasource::MemTable;
use datafusion::prelude::SessionContext;
use datafusion_bio_function_seq::register_seq_functions;

const STRING_TYPES: [DataType; 3] = [DataType::Utf8, DataType::LargeUtf8, DataType::Utf8View];

fn data(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Rows without the header. Trailing empty fields may have been stripped (the
/// trailing-whitespace hook removes trailing tabs), so rows are padded back to
/// the header's width.
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

/// `(id, seq, qual)` from FASTQ, or `(id, seq, "")` from FASTA.
fn records(name: &str) -> Vec<(String, String, String)> {
    let text = data(name);
    let lines: Vec<&str> = text.lines().collect();
    if name.ends_with(".fastq") {
        lines
            .chunks(4)
            .map(|r| (r[0][1..].to_string(), r[1].to_string(), r[3].to_string()))
            .collect()
    } else {
        lines
            .chunks(2)
            .map(|r| (r[0][1..].to_string(), r[1].to_string(), String::new()))
            .collect()
    }
}

/// A session with table `reads(id, seq, qual)`, where `seq` and `qual` have type `dt`.
fn session(name: &str, dt: &DataType) -> SessionContext {
    let recs = records(name);
    let col = |f: fn(&(String, String, String)) -> &str, dt: &DataType| -> ArrayRef {
        let utf8: ArrayRef = Arc::new(StringArray::from_iter_values(recs.iter().map(f)));
        cast(&utf8, dt).unwrap()
    };
    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("seq", dt.clone(), false),
        Field::new("qual", dt.clone(), false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            col(|r| &r.0, &DataType::Utf8),
            col(|r| &r.1, dt),
            col(|r| &r.2, dt),
        ],
    )
    .unwrap();
    let ctx = SessionContext::new();
    register_seq_functions(&ctx);
    ctx.register_table(
        "reads",
        Arc::new(MemTable::try_new(schema, vec![vec![batch]]).unwrap()),
    )
    .unwrap();
    ctx
}

async fn query(ctx: &SessionContext, sql: &str) -> Vec<RecordBatch> {
    ctx.sql(sql)
        .await
        .unwrap_or_else(|e| panic!("{sql}: {e}"))
        .collect()
        .await
        .unwrap_or_else(|e| panic!("{sql}: {e}"))
}

/// `id -> value` for a two-column string result.
async fn strings(ctx: &SessionContext, sql: &str, dt: &DataType) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for batch in query(ctx, sql).await {
        assert_eq!(batch.column(1).data_type(), dt, "{sql}: output type");
        let values = cast(batch.column(1), &DataType::Utf8).unwrap();
        let ids = batch.column(0).as_string::<i32>();
        for (id, v) in ids.iter().zip(values.as_string::<i32>().iter()) {
            out.insert(id.unwrap().to_string(), v.unwrap().to_string());
        }
    }
    out
}

/// `id -> value` for a two-column Float64 result.
async fn floats(ctx: &SessionContext, sql: &str) -> HashMap<String, Option<f64>> {
    let mut out = HashMap::new();
    for batch in query(ctx, sql).await {
        let ids = batch.column(0).as_string::<i32>();
        let values = batch.column(1).as_primitive::<Float64Type>();
        for (id, v) in ids.iter().zip(values.iter()) {
            out.insert(id.unwrap().to_string(), v);
        }
    }
    out
}

/// seqkit prints two decimals; accept anything that rounds to the printed value.
fn assert_printed(expected: &str, got: Option<f64>, what: &str) {
    if expected == "NaN" {
        // seqkit prints NaN for a zero denominator; the kernel returns null.
        assert_eq!(got, None, "{what}: seqkit printed NaN");
        return;
    }
    let expected: f64 = expected.parse().unwrap();
    match got {
        Some(v) => assert!(
            (v - expected).abs() <= 0.005 + 1e-9,
            "{what}: got {v}, seqkit printed {expected}"
        ),
        // seqkit prints 0.00 where the denominator is zero; the kernel returns null.
        None => assert_eq!(expected, 0.0, "{what}: got null"),
    }
}

#[tokio::test]
async fn revcomp_matches_seqkit_and_biopython() {
    let golden = tsv("revcomp.tsv");
    for dt in &STRING_TYPES {
        let ctx = session("corpus.fastq", dt);
        let got = strings(&ctx, "SELECT id, seq_revcomp(seq) FROM reads", dt).await;
        for row in &golden {
            assert_eq!(got[&row[0]], row[1], "{dt} {} vs seqkit", row[0]);
            assert_eq!(got[&row[0]], row[2], "{dt} {} vs Biopython", row[0]);
        }
    }
}

#[tokio::test]
async fn gc_content_matches_seqkit() {
    let golden = tsv("gc.tsv");
    for dt in &STRING_TYPES {
        let ctx = session("corpus.fastq", dt);
        let acgt = floats(
            &ctx,
            "SELECT id, seq_gc_content(seq, 'acgt') * 100 FROM reads",
        )
        .await;
        let length = floats(&ctx, "SELECT id, seq_gc_content(seq) * 100 FROM reads").await;
        for row in &golden {
            assert_printed(&row[1], acgt[&row[0]], &format!("{dt} {} acgt", row[0]));
            assert_printed(&row[2], length[&row[0]], &format!("{dt} {} length", row[0]));
        }
    }
}

#[tokio::test]
async fn mean_quality_matches_seqkit() {
    let golden = tsv("qual.tsv");
    for dt in &STRING_TYPES {
        let ctx = session("corpus.fastq", dt);
        let got = floats(&ctx, "SELECT id, seq_mean_quality(qual) FROM reads").await;
        for row in &golden {
            let v = got[&row[0]];
            assert!(v.is_some(), "{dt} {}: null", row[0]);
            assert_printed(&row[1], v, &format!("{dt} {}", row[0]));
        }
    }
}

#[tokio::test]
async fn translate_matches_biopython_and_seqkit() {
    let golden = tsv("translate.tsv");
    let mut combos: Vec<(String, String, String)> = golden
        .iter()
        .map(|r| (r[1].clone(), r[2].clone(), r[3].clone()))
        .collect();
    combos.sort();
    combos.dedup();
    let mut checked_seqkit = 0;
    for dt in &STRING_TYPES {
        let ctx = session("coding.fasta", dt);
        for (frame, table, to_stop) in &combos {
            let sql =
                format!("SELECT id, seq_translate(seq, {frame}, {table}, {to_stop}) FROM reads");
            let got = strings(&ctx, &sql, dt).await;
            for row in golden
                .iter()
                .filter(|r| (&r[1], &r[2], &r[3]) == (frame, table, to_stop))
            {
                let what = format!(
                    "{dt} {} frame {frame} table {table} to_stop {to_stop}",
                    row[0]
                );
                assert_eq!(got[&row[0]], row[5], "{what} vs Biopython");
                // seqkit misses some tables, has no to_stop, and does not
                // complement ambiguity codes in negative frames.
                if !row[4].is_empty() && !frame.starts_with('-') {
                    assert_eq!(got[&row[0]], row[4], "{what} vs seqkit");
                    checked_seqkit += 1;
                }
            }
        }
    }
    assert!(
        checked_seqkit > 1000,
        "only {checked_seqkit} seqkit comparisons"
    );
}
