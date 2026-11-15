//! SQL behaviour of the `seq_*` functions: nulls, options, errors, types.

use datafusion::arrow::array::{Array, AsArray};
use datafusion::arrow::datatypes::{DataType, Float64Type};
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
