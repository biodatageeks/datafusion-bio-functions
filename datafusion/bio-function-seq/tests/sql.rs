//! SQL behaviour of the `seq_*` functions: nulls, options, errors, types.

use std::sync::Arc;

use datafusion::arrow::array::{Array, AsArray, LargeStringArray, RecordBatch, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Float64Type, Schema};
use datafusion::datasource::MemTable;
use datafusion::prelude::SessionContext;
use datafusion_bio_function_seq::register_seq_functions;

fn ctx() -> SessionContext {
    let ctx = SessionContext::new();
    register_seq_functions(&ctx);
    ctx
}

async fn one_string(sql: &str) -> Option<String> {
    let batches = ctx().sql(sql).await.unwrap().collect().await.unwrap();
    let col = batches[0].column(0);
    let col = datafusion::arrow::compute::cast(col, &DataType::Utf8).unwrap();
    let a = col.as_string::<i32>();
    (!a.is_null(0)).then(|| a.value(0).to_string())
}

async fn one_float(sql: &str) -> Option<f64> {
    let batches = ctx().sql(sql).await.unwrap().collect().await.unwrap();
    let a = batches[0].column(0).as_primitive::<Float64Type>();
    (!a.is_null(0)).then(|| a.value(0))
}

async fn error(sql: &str) -> String {
    match ctx().sql(sql).await {
        Err(e) => e.to_string(),
        Ok(df) => df.collect().await.unwrap_err().to_string(),
    }
}

#[tokio::test]
async fn scalar_calls_and_defaults() {
    assert_eq!(
        one_string("SELECT seq_revcomp('ACGTN')").await.unwrap(),
        "NACGT"
    );
    assert_eq!(
        one_string("SELECT seq_revcomp('aacc')").await.unwrap(),
        "ggtt"
    );
    assert_eq!(
        one_string("SELECT seq_complement('ACGU', 'rna')")
            .await
            .unwrap(),
        "UGCA"
    );
    assert_eq!(
        one_string("SELECT seq_revcomp('ACGZ')").await.unwrap(),
        "ZCGT"
    );
    assert_eq!(
        one_string("SELECT seq_reverse('II5#')").await.unwrap(),
        "#5II"
    );
    assert_eq!(one_string("SELECT seq_reverse(NULL)").await, None);
    assert_eq!(
        one_string("SELECT seq_translate('ATGGCCTAAGC')")
            .await
            .unwrap(),
        "MA*"
    );
    assert_eq!(
        one_string("SELECT seq_translate('ATGGCCTAAGC', -1, 1, false)")
            .await
            .unwrap(),
        "A*A"
    );
    assert_eq!(
        one_string("SELECT seq_translate('ATGGCCTAAGC', 1, 1, true)")
            .await
            .unwrap(),
        "MA"
    );
    assert_eq!(one_float("SELECT seq_gc_content('GGNN')").await, Some(0.5));
    assert_eq!(
        one_float("SELECT seq_gc_content('GGNN', 'acgt')").await,
        Some(1.0)
    );
    let q = one_float("SELECT seq_mean_quality('I#')").await.unwrap();
    assert!((q - 5.01).abs() < 0.005);
    assert_eq!(
        one_float("SELECT seq_mean_quality('I#', 33, 'arithmetic')").await,
        Some(21.0)
    );
    assert_eq!(
        one_float("SELECT seq_mean_quality('hh', 64, 'arithmetic')").await,
        Some(40.0)
    );
}

#[tokio::test]
async fn nulls_and_zero_denominators() {
    assert_eq!(one_string("SELECT seq_revcomp(NULL)").await, None);
    assert_eq!(
        one_string("SELECT seq_translate(CAST(NULL AS VARCHAR))").await,
        None
    );
    assert_eq!(one_float("SELECT seq_gc_content(NULL)").await, None);
    assert_eq!(one_float("SELECT seq_gc_content('')").await, None);
    assert_eq!(
        one_float("SELECT seq_gc_content('NNN', 'acgt')").await,
        None
    );
    assert_eq!(one_float("SELECT seq_mean_quality('')").await, None);
    // A null gc never passes a filter, unlike NaN.
    let batches = ctx()
        .sql("SELECT s FROM (VALUES ('GGG'), (''), ('AAA')) t(s) WHERE seq_gc_content(s) > 0.3")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    assert_eq!(batches.iter().map(|b| b.num_rows()).sum::<usize>(), 1);
}

#[tokio::test]
async fn output_type_follows_input_type() {
    let batches = ctx()
        .sql(
            "SELECT seq_revcomp(arrow_cast('AC', 'Utf8View')) a, \
                    seq_translate(arrow_cast('ATG', 'LargeUtf8')) b, \
                    seq_complement('AC') c, \
                    seq_gc_content(arrow_cast('GC', 'Utf8View')) d",
        )
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let s = batches[0].schema();
    assert_eq!(s.field(0).data_type(), &DataType::Utf8View);
    assert_eq!(s.field(1).data_type(), &DataType::LargeUtf8);
    assert_eq!(s.field(2).data_type(), &DataType::Utf8);
    assert_eq!(s.field(3).data_type(), &DataType::Float64);
}

#[tokio::test]
async fn kernel_errors() {
    let e = error("SELECT seq_revcomp('ACGZ', 'dna', true)").await;
    assert!(e.contains("seq_revcomp") && e.contains("'Z'"), "{e}");
    let e = error("SELECT seq_mean_quality('II ')").await;
    assert!(e.contains("below the Phred offset 33"), "{e}");
    let e = error("SELECT seq_translate('ATG', 4)").await;
    assert!(e.contains("invalid frame 4"), "{e}");
    let e = error("SELECT seq_translate('ATG', 1, 7)").await;
    assert!(e.contains("unknown NCBI translation table 7"), "{e}");
    let e = error("SELECT seq_gc_content('ACGT', 'ignore_n')").await;
    assert!(e.contains("invalid denominator 'ignore_n'"), "{e}");
    let e = error("SELECT seq_revcomp('ACGT', 'protein')").await;
    assert!(e.contains("invalid alphabet"), "{e}");
    let e = error("SELECT seq_mean_quality('II', 300)").await;
    assert!(e.contains("offset must be between 0 and 255"), "{e}");
}

#[tokio::test]
async fn argument_errors() {
    let e = error("SELECT seq_revcomp(42)").await;
    assert!(e.contains("first argument must be a string"), "{e}");
    let e = error("SELECT seq_revcomp('A', 'dna', true, 1)").await;
    assert!(e.contains("takes 1 to 3 arguments"), "{e}");
    let e = error("SELECT seq_translate('ATG', 'one')").await;
    assert!(e.contains("frame must be an integer"), "{e}");
    let e = error("SELECT seq_revcomp('ACGT', NULL)").await;
    assert!(e.contains("alphabet must not be null"), "{e}");
    let e = error("SELECT seq_translate(s, f) FROM (VALUES ('ATG', 1)) t(s, f)").await;
    assert!(e.contains("frame must be a constant"), "{e}");
}

#[tokio::test]
async fn builtin_translate_is_untouched() {
    assert_eq!(
        one_string("SELECT translate('ACGT', 'AT', 'TA')")
            .await
            .unwrap(),
        "TCGA"
    );
}

#[tokio::test]
async fn sliced_input_arrays() {
    // A sliced array shares its value buffer with rows outside the slice.
    let utf8 = StringArray::from(vec!["AAAA", "ACGT", "ATGGCC", "TTTTTTTT"]);
    let large = LargeStringArray::from(vec!["AAAA", "ACGT", "ATGGCC", "TTTTTTTT"]);
    let schema = Arc::new(Schema::new(vec![
        Field::new("s", DataType::Utf8, false),
        Field::new("l", DataType::LargeUtf8, false),
    ]));
    let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(utf8), Arc::new(large)])
        .unwrap()
        .slice(1, 2);
    let ctx = ctx();
    ctx.register_table(
        "t",
        Arc::new(MemTable::try_new(schema, vec![vec![batch]]).unwrap()),
    )
    .unwrap();
    let batches = ctx
        .sql("SELECT seq_revcomp(s), seq_translate(l) FROM t")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let rc: Vec<_> = batches[0].column(0).as_string::<i32>().iter().collect();
    let aa: Vec<_> = batches[0].column(1).as_string::<i64>().iter().collect();
    assert_eq!(rc, vec![Some("ACGT"), Some("GGCCAT")]);
    assert_eq!(aa, vec![Some("T"), Some("MA")]);
}

#[tokio::test]
async fn all_q0_mean_quality_is_positive_zero() {
    let q = one_float("SELECT seq_mean_quality('!!!!')").await.unwrap();
    assert!(q == 0.0 && q.is_sign_positive(), "{q}");
}

// ---------------------------------------------------------------------------
// Biopython-parity functions (biodatageeks/polars-bio#512)
// ---------------------------------------------------------------------------

/// The first value of the first column, formatted by Arrow; `None` for null.
async fn one_display(sql: &str) -> Option<String> {
    let batches = ctx().sql(sql).await.unwrap().collect().await.unwrap();
    let col = batches[0].column(0);
    (!col.is_null(0))
        .then(|| datafusion::arrow::util::display::array_value_to_string(col, 0).unwrap())
}

async fn one_type(sql: &str) -> DataType {
    let batches = ctx().sql(sql).await.unwrap().collect().await.unwrap();
    batches[0].column(0).data_type().clone()
}

#[track_caller]
fn assert_close(got: Option<f64>, expected: f64) {
    let got = got.expect("a value");
    assert!(
        (got - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{got} != {expected}"
    );
}

#[tokio::test]
async fn translate_options() {
    let tr = |args: &'static str| async move {
        one_string(&format!("SELECT seq_translate({args})"))
            .await
            .unwrap()
    };
    assert_eq!(tr("'ATGTAAGCCTGA', 1, 1, false, '#'").await, "M#A#");
    assert_eq!(tr("'ATGTAAGCCTGA', 1, 1, false, ''").await, "MA");
    assert_eq!(tr("'GTGGCCTAA', 1, 11, false, '*', true").await, "MA");
    assert_eq!(
        tr("'RATSARMTTWTT', 1, 1, false, '*', false, 'biopython'").await,
        "BZJX"
    );
    assert_eq!(
        tr("'ATGRATTAA', 1, 1, false, '*', true, 'biopython'").await,
        "MB"
    );
    // Existing calls are unchanged.
    assert_eq!(tr("'RATSARMTTWTT'").await, "XXXX");
    let e = error("SELECT seq_translate('ATGTAAGCCTAA', 1, 1, false, '*', true)").await;
    assert!(e.contains("in-frame stop codon 'TAA'"), "{e}");
    assert!(e.contains("cds_status"), "{e}");
    let e = error("SELECT seq_translate('ATG', 1, 1, false, '*', false, 'bzj')").await;
    assert!(e.contains("invalid ambiguous 'bzj'"), "{e}");
    let e = error("SELECT seq_translate('ATG', 1, 1, false, '*', false, 'x', 1)").await;
    assert!(e.contains("takes 1 to 7 arguments"), "{e}");
}

#[tokio::test]
async fn cds_status() {
    for (seq, expected) in [
        ("ATGGCCTAA", None),
        ("GTGGCCTAA", Some("no_start_codon")),
        ("ATGGCCTA", Some("length_not_multiple_of_3")),
        ("ATGGCCGCC", Some("no_stop_codon")),
        ("ATGTAAGCCTAA", Some("internal_stop_codon")),
    ] {
        let got = one_string(&format!("SELECT seq_cds_status('{seq}')")).await;
        assert_eq!(got.as_deref(), expected, "{seq}");
    }
    assert_eq!(one_string("SELECT seq_cds_status(NULL)").await, None);
    // Arguments are frame, then table, as in seq_translate. TTAGGCCAC is
    // the reverse complement of GTGGCCTAA, a CDS in table 11 only.
    assert_eq!(
        one_string("SELECT seq_cds_status('TTAGGCCAC', -1, 11)").await,
        None
    );
    assert_eq!(
        one_string("SELECT seq_cds_status('TTAGGCCAC', -1, 1)")
            .await
            .unwrap(),
        "no_start_codon"
    );
    assert_eq!(
        one_type("SELECT seq_cds_status(arrow_cast('ATG', 'Utf8View'))").await,
        DataType::Utf8
    );
}

#[tokio::test]
async fn transcription() {
    assert_eq!(
        one_string("SELECT seq_transcribe('ACGTacgtN-U')")
            .await
            .unwrap(),
        "ACGUacguN-U"
    );
    assert_eq!(
        one_string("SELECT seq_back_transcribe('ACGUacgu')")
            .await
            .unwrap(),
        "ACGTacgt"
    );
    assert_eq!(one_string("SELECT seq_transcribe(NULL)").await, None);
    assert_eq!(
        one_type("SELECT seq_transcribe(arrow_cast('ACGT', 'LargeUtf8'))").await,
        DataType::LargeUtf8
    );
}

#[tokio::test]
async fn quality_summaries() {
    assert_close(
        one_float("SELECT seq_expected_errors('I#')").await,
        1e-4 + 10f64.powf(-0.2),
    );
    assert_eq!(one_float("SELECT seq_expected_errors('')").await, Some(0.0));
    assert_eq!(one_float("SELECT seq_expected_errors(NULL)").await, None);
    assert_close(
        one_float("SELECT seq_expected_errors('hh', 64)").await,
        2e-4,
    );
    assert_eq!(
        one_display("SELECT seq_min_quality('I5#')").await.unwrap(),
        "2"
    );
    assert_eq!(
        one_display("SELECT seq_max_quality('I5#')").await.unwrap(),
        "40"
    );
    assert_eq!(
        one_display("SELECT seq_max_quality('hB', 64)")
            .await
            .unwrap(),
        "40"
    );
    assert_eq!(one_display("SELECT seq_min_quality('')").await, None);
    assert_eq!(one_display("SELECT seq_max_quality(NULL)").await, None);
    assert_eq!(
        one_type("SELECT seq_min_quality('I')").await,
        DataType::Int32
    );
    let e = error("SELECT seq_min_quality('I ')").await;
    assert!(e.contains("below the Phred offset"), "{e}");
    let e = error("SELECT seq_expected_errors('I ')").await;
    assert!(e.contains("below the Phred offset"), "{e}");
}

#[tokio::test]
async fn orfs() {
    // 1-based closed coordinates by default: the stack's default.
    assert_eq!(
        one_display("SELECT seq_orfs('ATGGCCTAA', 1, 1)")
            .await
            .unwrap(),
        "[{start: 1, end: 9, strand: +, frame: 1, aa_len: 2, has_start: true, has_stop: true}]"
    );
    // table, min_len, strands, start_codons, allow_partial, include_protein, zero_based
    assert_eq!(
        one_display("SELECT seq_orfs('TTAGGCCAT', 1, 1, 'both', 'ATG', false, true, true)")
            .await
            .unwrap(),
        "[{start: 0, end: 9, strand: -, frame: -1, aa_len: 2, has_start: true, has_stop: true, protein: MA}]"
    );
    assert_eq!(
        one_display("SELECT seq_orfs('ATGGCCTAA', 1, 1, '-')")
            .await
            .unwrap(),
        "[]"
    );
    // Partial ORFs of the forward strand, one per frame, sorted by start.
    assert_eq!(
        one_display("SELECT seq_orfs('GCCGCCTAA', 1, 1, '+', 'ATG', true)")
            .await
            .unwrap(),
        "[{start: 1, end: 9, strand: +, frame: 1, aa_len: 2, has_start: false, has_stop: true}, \
          {start: 2, end: 7, strand: +, frame: 2, aa_len: 2, has_start: false, has_stop: false}, \
          {start: 3, end: 8, strand: +, frame: 3, aa_len: 2, has_start: false, has_stop: false}]"
    );
    // The default min_len is 30 amino acids.
    assert_eq!(
        one_display("SELECT seq_orfs('ATGGCCTAA')").await.unwrap(),
        "[]"
    );
    assert_eq!(one_display("SELECT seq_orfs(NULL)").await, None);
    let e = error("SELECT seq_orfs('ATG', 7)").await;
    assert!(e.contains("unknown NCBI translation table 7"), "{e}");
    let e = error("SELECT seq_orfs('ATG', 1, 1, 'plus')").await;
    assert!(e.contains("invalid strands 'plus'"), "{e}");

    let fields = |t: DataType| match t {
        DataType::List(item) => match item.data_type() {
            DataType::Struct(fields) => fields
                .iter()
                .map(|f| format!("{}:{}", f.name(), f.data_type()))
                .collect::<Vec<_>>()
                .join(","),
            other => panic!("{other}"),
        },
        other => panic!("{other}"),
    };
    let base = "start:Int64,end:Int64,strand:Utf8,frame:Int32,aa_len:Int64,has_start:Boolean,has_stop:Boolean";
    assert_eq!(fields(one_type("SELECT seq_orfs('ATG')").await), base);
    assert_eq!(
        fields(one_type("SELECT seq_orfs('ATG', 1, 30, 'both', 'table', false, true)").await),
        format!("{base},protein:Utf8")
    );
}

#[tokio::test]
async fn find_motif() {
    let find = |args: &'static str| async move {
        one_display(&format!("SELECT seq_find_motif({args})")).await
    };
    // 1-based by default; zero_based is the last argument.
    assert_eq!(find("'ACGTACGT', 'ACG'").await.unwrap(), "[1, 5]");
    assert_eq!(
        find("'ACGTACGT', 'ACG', false, true, true").await.unwrap(),
        "[0, 4]"
    );
    assert_eq!(
        find("'AAAAA', 'AA', false, false, true").await.unwrap(),
        "[0, 2]"
    );
    assert_eq!(
        find("'AACGTT', 'AAC', true, true, true").await.unwrap(),
        "[0, 3]"
    );
    assert_eq!(
        find("'ACGTNACGT', 'N', false, true, true").await.unwrap(),
        "[0, 1, 2, 3, 5, 6, 7, 8]"
    );
    assert_eq!(find("'ACGT', 'TTT'").await.unwrap(), "[]");
    assert_eq!(find("NULL, 'ACG'").await, None);
    assert_eq!(
        one_type("SELECT seq_find_motif('ACGT', 'A')")
            .await
            .to_string(),
        "List(Int64)"
    );
    let e = error("SELECT seq_find_motif('ACGT')").await;
    assert!(e.contains("takes 2 to 5 arguments"), "{e}");
    for pattern in ["", "AC!"] {
        let e = error(&format!("SELECT seq_find_motif('ACGT', '{pattern}')")).await;
        assert!(e.contains("invalid pattern"), "{e}");
    }
}

#[tokio::test]
async fn windows() {
    assert_eq!(
        one_display("SELECT seq_gc_skew('GGGGCCCCGGGCAAAA', 4)")
            .await
            .unwrap(),
        "[1.0, -1.0, 0.5, 0.0]"
    );
    // The default window is 100, as in Biopython.
    assert_eq!(
        one_display("SELECT seq_gc_skew('GGCC')").await.unwrap(),
        "[0.0]"
    );
    assert_eq!(
        one_display("SELECT seq_gc_skew('GGCCGG', 4, 2)")
            .await
            .unwrap(),
        "[0.0, 0.0, 1.0]"
    );
    assert_eq!(
        one_display("SELECT seq_gc_content_windows('GGCCAAAT', 4)")
            .await
            .unwrap(),
        "[1.0, 0.0]"
    );
    assert_eq!(one_display("SELECT seq_gc_skew(NULL, 4)").await, None);
    assert_eq!(
        one_display("SELECT seq_gc_skew('', 4)").await.unwrap(),
        "[]"
    );

    // A NULL step means step = window; a zero denominator is a null element.
    let batches = ctx()
        .sql("SELECT seq_gc_content_windows('GGNNNNNN', 4, NULL, 'acgt')")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let list = batches[0].column(0).as_list::<i32>().value(0);
    let values: Vec<_> = list.as_primitive::<Float64Type>().iter().collect();
    assert_eq!(values, [Some(1.0), None]);

    let e = error("SELECT seq_gc_content_windows('ACGT')").await;
    assert!(e.contains("takes 2 to 4 arguments"), "{e}");
    let e = error("SELECT seq_gc_content_windows('ACGT', 0)").await;
    assert!(e.contains("invalid window '0'"), "{e}");
    let e = error("SELECT seq_gc_skew('ACGT', 4, 0)").await;
    assert!(e.contains("invalid step '0'"), "{e}");
}

#[tokio::test]
async fn molecular_weight() {
    let mw = |args: &'static str| async move {
        one_float(&format!("SELECT seq_molecular_weight({args})")).await
    };
    // Bio.SeqUtils.molecular_weight, Biopython 1.88.
    assert_close(mw("'ACGT'").await, 1253.8027);
    assert_close(mw("'ACGT', 'DNA', true").await, 2507.6054);
    assert_close(mw("'ACGT', 'DNA', false, true").await, 1235.7874);
    assert_close(mw("'ACGT', 'DNA', false, false, true").await, 1253.213097);
    assert_close(mw("'ACGU', 'RNA'").await, 1303.7737);
    assert_close(mw("'MKV', 'protein'").await, 376.51460000000003);
    assert_eq!(mw("'ACGN'").await, None);
    assert_eq!(mw("NULL").await, None);
    let e = error("SELECT seq_molecular_weight('MKV', 'protein', true)").await;
    assert!(e.contains("invalid double_stranded"), "{e}");
}

#[tokio::test]
async fn hash() {
    // XXH3-64, seed 0. These values must never change.
    assert_eq!(
        one_display("SELECT seq_hash('ACGT')").await.unwrap(),
        "17060466519357990319"
    );
    assert_eq!(
        one_display("SELECT seq_hash('TTAGGCCAT')").await.unwrap(),
        "17056848131934351943"
    );
    assert_eq!(
        one_display("SELECT seq_hash('TTAGGCCAT', true)")
            .await
            .unwrap(),
        "2721565146873327973"
    );
    assert_eq!(one_display("SELECT seq_hash(NULL)").await, None);
    assert_eq!(one_type("SELECT seq_hash('A')").await, DataType::UInt64);
}

#[tokio::test]
async fn trim_quality() {
    let trim = |args: &'static str| async move {
        one_display(&format!("SELECT seq_trim_quality({args})")).await
    };
    // seq, qual, leading, trailing, window_size, window_quality, offset
    assert_eq!(
        trim("'ACGTAC', '!!II!!', 10").await.unwrap(),
        "{sequence: GTAC, quality_scores: II!!}"
    );
    assert_eq!(
        trim("'ACGTAC', '!!II!!', NULL, 10").await.unwrap(),
        "{sequence: ACGT, quality_scores: !!II}"
    );
    assert_eq!(
        trim("'ACGTACGTACGT', 'IIII!!!!IIII', NULL, NULL, 4, 15")
            .await
            .unwrap(),
        "{sequence: ACGT, quality_scores: IIII}"
    );
    assert_eq!(
        trim("'ACGT', 'BBhh', 10, NULL, NULL, NULL, 64")
            .await
            .unwrap(),
        "{sequence: GT, quality_scores: hh}"
    );
    // A read Trimmomatic would drop gives two empty strings.
    assert_eq!(
        trim("'ACGTA', 'I!!!!', NULL, 10").await.unwrap(),
        "{sequence: , quality_scores: }"
    );
    // A null in either input gives null.
    assert_eq!(trim("NULL, 'II', 10").await, None);
    assert_eq!(trim("'AC', NULL, 10").await, None);

    let e = error("SELECT seq_trim_quality('AC', 'II')").await;
    assert!(e.contains("at least one of leading"), "{e}");
    let e = error("SELECT seq_trim_quality('AC', 'II', NULL, NULL, 4)").await;
    assert!(e.contains("window_size and window_quality"), "{e}");
    let e = error("SELECT seq_trim_quality('ACG', 'II', 10)").await;
    assert!(e.contains("differs from quality length"), "{e}");
    let e = error("SELECT seq_trim_quality('ACG')").await;
    assert!(e.contains("takes 2 to 7 arguments"), "{e}");
    let e = error("SELECT seq_trim_quality('ACG', 42, 10)").await;
    assert!(e.contains("must be a string"), "{e}");
}

#[tokio::test]
async fn new_functions_on_sliced_and_typed_arrays() {
    // A sliced batch with a null row; sequence as Utf8View, quality as LargeUtf8.
    use datafusion::arrow::array::StringViewArray;
    let seq = StringViewArray::from(vec![Some("TTTT"), Some("ATGGCCTAA"), None, Some("ACGTAC")]);
    let qual = LargeStringArray::from(vec![
        Some("IIII"),
        Some("IIIIIIII!"),
        Some("II"),
        Some("!!II!!"),
    ]);
    let schema = Arc::new(Schema::new(vec![
        Field::new("s", DataType::Utf8View, true),
        Field::new("q", DataType::LargeUtf8, true),
    ]));
    let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(seq), Arc::new(qual)])
        .unwrap()
        .slice(1, 3);
    let ctx = ctx();
    ctx.register_table(
        "t",
        Arc::new(MemTable::try_new(schema, vec![vec![batch]]).unwrap()),
    )
    .unwrap();
    let batches = ctx
        .sql(
            "SELECT seq_orfs(s, 1, 1) o, seq_find_motif(s, 'AC') m, seq_trim_quality(s, q, 10) t, \
                    seq_hash(s) h, seq_cds_status(s) c, seq_transcribe(s) r FROM t",
        )
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let b = &batches[0];
    let show = |col: usize| -> Vec<Option<String>> {
        (0..b.num_rows())
            .map(|row| {
                let a = b.column(col);
                (!a.is_null(row)).then(|| {
                    datafusion::arrow::util::display::array_value_to_string(a, row).unwrap()
                })
            })
            .collect()
    };
    let some = |s: &str| Some(s.to_string());
    assert_eq!(
        show(0),
        [
            some(
                "[{start: 1, end: 9, strand: +, frame: 1, aa_len: 2, has_start: true, has_stop: true}]"
            ),
            None,
            some("[]")
        ]
    );
    assert_eq!(show(1), [some("[]"), None, some("[1, 5]")]);
    assert_eq!(
        show(2),
        [
            some("{sequence: ATGGCCTAA, quality_scores: IIIIIIII!}"),
            None,
            some("{sequence: GTAC, quality_scores: II!!}")
        ]
    );
    assert_eq!(show(3)[1], None);
    assert_eq!(show(3)[0], some("2721565146873327973"));
    assert_eq!(show(4), [None, None, some("no_start_codon")]);
    assert_eq!(show(5), [some("AUGGCCUAA"), None, some("ACGUAC")]);
    // The struct fields keep the input string types.
    let fields = match b.column(2).data_type() {
        DataType::Struct(fields) => fields.clone(),
        other => panic!("{other}"),
    };
    assert_eq!(fields[0].name(), "sequence");
    assert_eq!(fields[0].data_type(), &DataType::Utf8View);
    assert_eq!(fields[1].name(), "quality_scores");
    assert_eq!(fields[1].data_type(), &DataType::LargeUtf8);
}

#[tokio::test]
async fn new_function_argument_errors() {
    let e = error("SELECT seq_hash('ACGT', NULL)").await;
    assert!(e.contains("canonical must not be null"), "{e}");
    let e = error("SELECT seq_find_motif(s, p) FROM (VALUES ('ACGT', 'A')) t(s, p)").await;
    assert!(e.contains("pattern must be a constant"), "{e}");
    let e = error("SELECT seq_gc_skew('ACGT', 'four')").await;
    assert!(e.contains("window must be an integer"), "{e}");
    let e = error("SELECT seq_min_quality('II', 300)").await;
    assert!(e.contains("offset must be between 0 and 255"), "{e}");
}
