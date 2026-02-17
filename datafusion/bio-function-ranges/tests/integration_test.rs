use std::sync::Arc;

use datafusion::arrow::array::{Int64Array, RecordBatch};
use datafusion::arrow::util::pretty::pretty_format_batches;
use datafusion::assert_batches_sorted_eq;
use datafusion::common::assert_contains;
use datafusion::config::ConfigOptions;
use datafusion::error::Result;
use datafusion::prelude::{SessionConfig, SessionContext};

use datafusion_bio_function_ranges::session_context::{Algorithm, BioConfig, BioSessionExt};
use datafusion_bio_function_ranges::{
    CountOverlapsProvider, FilterOp, create_bio_session, register_ranges_functions,
};

const READS_PATH: &str = "../../testing/data/interval/reads.csv";
const TARGETS_PATH: &str = "../../testing/data/interval/targets.csv";
const RANGES_READS_PATH: &str = "../../testing/data/ranges/reads.csv";
const RANGES_TARGETS_PATH: &str = "../../testing/data/ranges/targets.csv";
const EXONS_PATH: &str = "../../testing/data/ranges/exons/";
const FBRAIN_PATH: &str = "../../testing/data/ranges/fBrain-DS14718/";
const EXPECTED_COVERAGE_PATH: &str = "../../testing/data/ranges/expected_coverage.parquet";

#[rstest::fixture]
fn ctx() -> SessionContext {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false);

    SessionContext::new_with_bio(config)
}

async fn init_tables(ctx: &SessionContext) -> Result<()> {
    let reads = format!(
        "CREATE EXTERNAL TABLE reads STORED AS CSV LOCATION '{READS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(reads.as_str()).await?;

    let targets = format!(
        "CREATE EXTERNAL TABLE targets STORED AS CSV LOCATION '{TARGETS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(targets.as_str()).await?;

    Ok(())
}

#[rstest::fixture]
#[once]
fn expected_equi() -> [&'static str; 20] {
    [
        "+--------+-----------+---------+--------+-----------+---------+",
        "| contig | pos_start | pos_end | contig | pos_start | pos_end |",
        "+--------+-----------+---------+--------+-----------+---------+",
        "| chr1   | 150       | 250     | chr1   | 100       | 190     |",
        "| chr1   | 150       | 250     | chr1   | 200       | 290     |",
        "| chr1   | 190       | 300     | chr1   | 100       | 190     |",
        "| chr1   | 190       | 300     | chr1   | 200       | 290     |",
        "| chr1   | 300       | 501     | chr1   | 400       | 600     |",
        "| chr1   | 500       | 700     | chr1   | 400       | 600     |",
        "| chr1   | 15000     | 15000   | chr1   | 10000     | 20000   |",
        "| chr1   | 22000     | 22300   | chr1   | 22100     | 22100   |",
        "| chr2   | 150       | 250     | chr2   | 100       | 190     |",
        "| chr2   | 150       | 250     | chr2   | 200       | 290     |",
        "| chr2   | 190       | 300     | chr2   | 100       | 190     |",
        "| chr2   | 190       | 300     | chr2   | 200       | 290     |",
        "| chr2   | 300       | 500     | chr2   | 400       | 600     |",
        "| chr2   | 500       | 700     | chr2   | 400       | 600     |",
        "| chr2   | 15000     | 15000   | chr2   | 10000     | 20000   |",
        "| chr2   | 22000     | 22300   | chr2   | 22100     | 22100   |",
        "+--------+-----------+---------+--------+-----------+---------+",
    ]
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
#[case::hash_join(None)]
#[case::interval_join_coitrees(Some(Algorithm::Coitrees))]
#[case::interval_join_interval_tree(Some(Algorithm::IntervalTree))]
#[case::interval_join_array_interval_tree(Some(Algorithm::ArrayIntervalTree))]
#[case::interval_join_lapper(Some(Algorithm::Lapper))]
#[case::interval_join_superintervals(Some(Algorithm::SuperIntervals))]
async fn test_equi_and_range_condition(
    #[case] algorithm: Option<Algorithm>,
    ctx: SessionContext,
    expected_equi: &[&str; 20],
) -> Result<()> {
    init_tables(&ctx).await?;

    ctx.sql(format!("SET bio.prefer_interval_join = {}", algorithm.is_some()).as_str())
        .await?;
    ctx.sql(
        format!(
            "SET bio.interval_join_algorithm = {}",
            algorithm.unwrap_or_default()
        )
        .as_str(),
    )
    .await?;

    let query = r#"SELECT *
               FROM reads
               JOIN targets
               ON reads.contig = targets.contig
                  AND reads.pos_start <= targets.pos_end
                  AND reads.pos_end >= targets.pos_start
               ORDER BY reads.contig, reads.pos_start, reads.pos_end,
                        targets.contig, targets.pos_start, targets.pos_end"#;

    let plan: Vec<RecordBatch> = ctx
        .sql(format!("EXPLAIN {query}").as_str())
        .await?
        .collect()
        .await?;
    let formatted: String = pretty_format_batches(&plan)?.to_string();
    let expected_plan = match algorithm {
        None => "HashJoinExec: mode=CollectLeft, join_type=Inner, on=[(contig@0, contig@0)], filter=pos_start@0 <= pos_end@3 AND pos_end@1 >= pos_start@2".to_string(),
        Some(alg) => format!("IntervalJoinExec: mode=CollectLeft, join_type=Inner, on=[(contig@0, contig@0)], filter=pos_start@0 <= pos_end@3 AND pos_end@1 >= pos_start@2, alg={alg}"),
    };
    assert_contains!(formatted, expected_plan);

    let result: Vec<RecordBatch> = ctx.sql(query).await?.collect().await?;
    assert_batches_sorted_eq!(expected_equi, &result);

    Ok(())
}

#[rstest::fixture]
#[once]
fn expected_range() -> [&'static str; 36] {
    [
        "+--------+-----------+---------+--------+-----------+---------+",
        "| contig | pos_start | pos_end | contig | pos_start | pos_end |",
        "+--------+-----------+---------+--------+-----------+---------+",
        "| chr1   | 150       | 250     | chr1   | 100       | 190     |",
        "| chr1   | 150       | 250     | chr1   | 200       | 290     |",
        "| chr1   | 150       | 250     | chr2   | 100       | 190     |",
        "| chr1   | 150       | 250     | chr2   | 200       | 290     |",
        "| chr1   | 190       | 300     | chr1   | 100       | 190     |",
        "| chr1   | 190       | 300     | chr1   | 200       | 290     |",
        "| chr1   | 190       | 300     | chr2   | 100       | 190     |",
        "| chr1   | 190       | 300     | chr2   | 200       | 290     |",
        "| chr1   | 300       | 501     | chr1   | 400       | 600     |",
        "| chr1   | 300       | 501     | chr2   | 400       | 600     |",
        "| chr1   | 500       | 700     | chr1   | 400       | 600     |",
        "| chr1   | 500       | 700     | chr2   | 400       | 600     |",
        "| chr1   | 15000     | 15000   | chr1   | 10000     | 20000   |",
        "| chr1   | 15000     | 15000   | chr2   | 10000     | 20000   |",
        "| chr1   | 22000     | 22300   | chr1   | 22100     | 22100   |",
        "| chr1   | 22000     | 22300   | chr2   | 22100     | 22100   |",
        "| chr2   | 150       | 250     | chr1   | 100       | 190     |",
        "| chr2   | 150       | 250     | chr1   | 200       | 290     |",
        "| chr2   | 150       | 250     | chr2   | 100       | 190     |",
        "| chr2   | 150       | 250     | chr2   | 200       | 290     |",
        "| chr2   | 190       | 300     | chr1   | 100       | 190     |",
        "| chr2   | 190       | 300     | chr1   | 200       | 290     |",
        "| chr2   | 190       | 300     | chr2   | 100       | 190     |",
        "| chr2   | 190       | 300     | chr2   | 200       | 290     |",
        "| chr2   | 300       | 500     | chr1   | 400       | 600     |",
        "| chr2   | 300       | 500     | chr2   | 400       | 600     |",
        "| chr2   | 500       | 700     | chr1   | 400       | 600     |",
        "| chr2   | 500       | 700     | chr2   | 400       | 600     |",
        "| chr2   | 15000     | 15000   | chr1   | 10000     | 20000   |",
        "| chr2   | 15000     | 15000   | chr2   | 10000     | 20000   |",
        "| chr2   | 22000     | 22300   | chr1   | 22100     | 22100   |",
        "| chr2   | 22000     | 22300   | chr2   | 22100     | 22100   |",
        "+--------+-----------+---------+--------+-----------+---------+",
    ]
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
#[case::nested_loop_join(None)]
#[case::interval_join_coitrees(Some(Algorithm::Coitrees))]
#[case::interval_join_interval_tree(Some(Algorithm::IntervalTree))]
#[case::interval_join_array_interval_tree(Some(Algorithm::ArrayIntervalTree))]
#[case::interval_join_lapper(Some(Algorithm::Lapper))]
#[case::interval_join_superintervals(Some(Algorithm::SuperIntervals))]
async fn test_range_condition(
    #[case] algorithm: Option<Algorithm>,
    ctx: SessionContext,
    expected_range: &[&str; 36],
) -> Result<()> {
    init_tables(&ctx).await?;

    ctx.sql(format!("SET bio.prefer_interval_join = {}", algorithm.is_some()).as_str())
        .await?;
    ctx.sql(
        format!(
            "SET bio.interval_join_algorithm = {}",
            algorithm.unwrap_or_default()
        )
        .as_str(),
    )
    .await?;

    let query = r#"SELECT *
               FROM reads
               JOIN targets
               ON reads.pos_start <= targets.pos_end AND reads.pos_end >= targets.pos_start
               ORDER BY reads.contig, reads.pos_start, reads.pos_end,
                        targets.contig, targets.pos_start, targets.pos_end"#;

    let plan: Vec<RecordBatch> = ctx
        .sql(format!("EXPLAIN {query}").as_str())
        .await?
        .collect()
        .await?;
    let formatted: String = pretty_format_batches(&plan)?.to_string();
    let expected_plan = match algorithm {
        None => "NestedLoopJoinExec: join_type=Inner, filter=pos_start@0 <= pos_end@3 AND pos_end@1 >= pos_start@2".to_string(),
        Some(alg) => format!("IntervalJoinExec: mode=CollectLeft, join_type=Inner, on=[(1, 1)], filter=pos_start@0 <= pos_end@3 AND pos_end@1 >= pos_start@2, alg={alg}"),
    };
    assert_contains!(formatted, expected_plan);

    let result: Vec<RecordBatch> = ctx.sql(query).await?.collect().await?;
    assert_batches_sorted_eq!(expected_range, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_all_gteq_lteq_conditions(ctx: SessionContext) -> Result<()> {
    let a = r#"
        CREATE TABLE a (contig TEXT, start INTEGER, end INTEGER) AS VALUES
        ('a', 5, 10)
    "#;

    let b = r#"
        CREATE TABLE b (contig TEXT, start INTEGER, end INTEGER) AS VALUES
        ('a', 11, 15),
        ('a', 10, 15),
        ('a', 10, 10),
        ('a',  9, 15),
        ('a',  5, 15),
        ('a',  4, 15),
        ('a',  4, 10),
        ('a',  6, 8),
        ('a',  4, 8),
        ('a',  4, 5),
        ('a',  5, 5),
        ('a',  4, 4)
    "#;

    let q0 = r#"
        SELECT * FROM a JOIN b
        ON a.contig = b.contig AND a.start <= b.end AND a.end >= b.start
    "#;

    let q1 = r#"
        SELECT a.*, b.* FROM b JOIN a
        ON a.contig = b.contig AND a.start <= b.end AND a.end >= b.start
    "#;

    let q2 = r#"
        SELECT a.*, b.* FROM a, b
        WHERE a.contig = b.contig AND a.start <= b.end AND a.end >= b.start
    "#;

    let q3 = r#"
        SELECT a.*, b.* FROM b, a
        WHERE a.contig = b.contig AND b.start <= a.end AND b.end >= a.start
    "#;

    ctx.sql(a).await?;
    ctx.sql(b).await?;

    let expected = vec![
        "+--------+-------+-----+--------+-------+-----+",
        "| contig | start | end | contig | start | end |",
        "+--------+-------+-----+--------+-------+-----+",
        "| a      | 5     | 10  | a      | 10    | 15  |",
        "| a      | 5     | 10  | a      | 10    | 10  |",
        "| a      | 5     | 10  | a      | 9     | 15  |",
        "| a      | 5     | 10  | a      | 5     | 15  |",
        "| a      | 5     | 10  | a      | 4     | 15  |",
        "| a      | 5     | 10  | a      | 4     | 10  |",
        "| a      | 5     | 10  | a      | 6     | 8   |",
        "| a      | 5     | 10  | a      | 4     | 8   |",
        "| a      | 5     | 10  | a      | 5     | 5   |",
        "| a      | 5     | 10  | a      | 4     | 5   |",
        "+--------+-------+-----+--------+-------+-----+",
    ];

    let results = ctx.sql(q0).await?.collect().await?;
    assert_batches_sorted_eq!(expected, &results);

    let results = ctx.sql(q1).await?.collect().await?;
    assert_batches_sorted_eq!(expected, &results);

    let results = ctx.sql(q2).await?.collect().await?;
    assert_batches_sorted_eq!(expected, &results);

    let results = ctx.sql(q3).await?.collect().await?;
    assert_batches_sorted_eq!(expected, &results);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_all_gt_lt_conditions(ctx: SessionContext) -> Result<()> {
    let a = r#"
        CREATE TABLE a (contig TEXT, start INTEGER, end INTEGER) AS VALUES
        ('a', 5, 10)
    "#;

    let b = r#"
        CREATE TABLE b (contig TEXT, start INTEGER, end INTEGER) AS VALUES
        ('a', 11, 15),
        ('a', 10, 15),
        ('a', 10, 10),
        ('a',  9, 15),
        ('a',  5, 15),
        ('a',  4, 15),
        ('a',  4, 10),
        ('a',  6, 8),
        ('a',  4, 8),
        ('a',  4, 5),
        ('a',  5, 5),
        ('a',  4, 4)
    "#;

    let q0 = r#"
        SELECT * FROM a JOIN b
        ON a.contig = b.contig AND a.start < b.end AND a.end > b.start
    "#;

    let q1 = r#"
        SELECT a.*, b.* FROM b JOIN a
        ON a.contig = b.contig AND a.end > b.start AND a.start < b.end
    "#;

    ctx.sql(a).await?;
    ctx.sql(b).await?;

    let expected = [
        "+--------+-------+-----+--------+-------+-----+",
        "| contig | start | end | contig | start | end |",
        "+--------+-------+-----+--------+-------+-----+",
        "| a      | 5     | 10  | a      | 9     | 15  |",
        "| a      | 5     | 10  | a      | 5     | 15  |",
        "| a      | 5     | 10  | a      | 4     | 15  |",
        "| a      | 5     | 10  | a      | 4     | 10  |",
        "| a      | 5     | 10  | a      | 6     | 8   |",
        "| a      | 5     | 10  | a      | 4     | 8   |",
        "+--------+-------+-----+--------+-------+-----+",
    ];

    let results = ctx.sql(q0).await?.collect().await?;
    assert_batches_sorted_eq!(expected, &results);

    let results = ctx.sql(q1).await?.collect().await?;
    assert_batches_sorted_eq!(expected, &results);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_nearest(ctx: SessionContext) -> Result<()> {
    let a = r#"
        CREATE TABLE a (contig TEXT, strand TEXT, start INTEGER, end INTEGER) AS VALUES
        ('a', 's', 5, 10)
    "#;

    let b = r#"
        CREATE TABLE b (contig TEXT, strand TEXT, start INTEGER, end INTEGER) AS VALUES
        ('a', 's', 11, 13),
        ('a', 's', 20, 21),
        ('a', 'x', 0, 1),
        ('b', 's', 1, 2)
    "#;

    ctx.sql("SET bio.interval_join_algorithm TO CoitreesNearest")
        .await?;

    ctx.sql(a).await?;
    ctx.sql(b).await?;

    let q = r#"
        SELECT * FROM a JOIN b
        ON a.contig = b.contig AND a.strand = b.strand
            AND a.start < b.end AND a.end > b.start
    "#;

    let result = ctx.sql(q).await?;
    result.clone().show().await?;

    let results = result.collect().await?;

    let expected = [
        "+--------+--------+-------+-----+--------+--------+-------+-----+",
        "| contig | strand | start | end | contig | strand | start | end |",
        "+--------+--------+-------+-----+--------+--------+-------+-----+",
        "|        |        |       |     | a      | x      | 0     | 1   |",
        "|        |        |       |     | b      | s      | 1     | 2   |",
        "| a      | s      | 5     | 10  | a      | s      | 11    | 13  |",
        "| a      | s      | 5     | 10  | a      | s      | 20    | 21  |",
        "+--------+--------+-------+-----+--------+--------+-------+-----+",
    ];

    assert_batches_sorted_eq!(expected, &results);

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Count overlaps, coverage, and nearest tests using ranges test data
// ─────────────────────────────────────────────────────────────────────────────

async fn init_ranges_tables(ctx: &SessionContext) -> Result<()> {
    let reads = format!(
        "CREATE EXTERNAL TABLE reads STORED AS CSV LOCATION '{RANGES_READS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(reads.as_str()).await?;

    let targets = format!(
        "CREATE EXTERNAL TABLE targets STORED AS CSV LOCATION '{RANGES_TARGETS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(targets.as_str()).await?;

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_count_overlaps(ctx: SessionContext) -> Result<()> {
    init_ranges_tables(&ctx).await?;

    let targets_schema = ctx.table("targets").await?.schema().as_arrow().clone();

    let provider = CountOverlapsProvider::new(
        Arc::new(ctx.clone()),
        "reads".to_string(),
        "targets".to_string(),
        targets_schema,
        vec![
            "contig".to_string(),
            "pos_start".to_string(),
            "pos_end".to_string(),
        ],
        vec![
            "contig".to_string(),
            "pos_start".to_string(),
            "pos_end".to_string(),
        ],
        FilterOp::Weak,
        false,
    );

    ctx.register_table("count_overlaps_result", Arc::new(provider))?;
    let result: Vec<RecordBatch> = ctx
        .sql("SELECT * FROM count_overlaps_result ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------+",
        "| contig | pos_start | pos_end | count |",
        "+--------+-----------+---------+-------+",
        "| chr1   | 100       | 190     | 2     |",
        "| chr1   | 200       | 290     | 2     |",
        "| chr1   | 400       | 600     | 2     |",
        "| chr1   | 10000     | 20000   | 1     |",
        "| chr1   | 22100     | 22100   | 1     |",
        "| chr2   | 100       | 190     | 2     |",
        "| chr2   | 200       | 290     | 2     |",
        "| chr2   | 400       | 600     | 2     |",
        "| chr2   | 10000     | 20000   | 1     |",
        "| chr2   | 22100     | 22100   | 1     |",
        "| chr3   | 100       | 200     | 0     |",
        "+--------+-----------+---------+-------+",
    ];

    assert_batches_sorted_eq!(expected, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_coverage_csv(ctx: SessionContext) -> Result<()> {
    init_ranges_tables(&ctx).await?;

    let targets_schema = ctx.table("targets").await?.schema().as_arrow().clone();

    let provider = CountOverlapsProvider::new(
        Arc::new(ctx.clone()),
        "reads".to_string(),
        "targets".to_string(),
        targets_schema,
        vec![
            "contig".to_string(),
            "pos_start".to_string(),
            "pos_end".to_string(),
        ],
        vec![
            "contig".to_string(),
            "pos_start".to_string(),
            "pos_end".to_string(),
        ],
        FilterOp::Weak,
        true,
    );

    ctx.register_table("coverage_result", Arc::new(provider))?;
    let result: Vec<RecordBatch> = ctx
        .sql("SELECT * FROM coverage_result ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    // Coverage values computed by get_coverage formula on merged reads:
    // chr1 merged reads: [150,700], [15000,15000], [22000,22300]
    // chr2 merged reads: [150,700], [15000,15000], [22000,22300]
    // chr3 merged reads: [234,300]
    let expected = [
        "+--------+-----------+---------+----------+",
        "| contig | pos_start | pos_end | coverage |",
        "+--------+-----------+---------+----------+",
        "| chr1   | 100       | 190     | 41       |",
        "| chr1   | 200       | 290     | 92       |",
        "| chr1   | 400       | 600     | 202      |",
        "| chr1   | 10000     | 20000   | 1        |",
        "| chr1   | 22100     | 22100   | 2        |",
        "| chr2   | 100       | 190     | 41       |",
        "| chr2   | 200       | 290     | 92       |",
        "| chr2   | 400       | 600     | 202      |",
        "| chr2   | 10000     | 20000   | 1        |",
        "| chr2   | 22100     | 22100   | 2        |",
        "| chr3   | 100       | 200     | 0        |",
        "+--------+-----------+---------+----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_coverage_parquet(ctx: SessionContext) -> Result<()> {
    // Register parquet tables
    let fbrain = format!("CREATE EXTERNAL TABLE fbrain STORED AS PARQUET LOCATION '{FBRAIN_PATH}'");
    ctx.sql(fbrain.as_str()).await?;

    let exons = format!("CREATE EXTERNAL TABLE exons STORED AS PARQUET LOCATION '{EXONS_PATH}'");
    ctx.sql(exons.as_str()).await?;

    let exons_schema = ctx.table("exons").await?.schema().as_arrow().clone();

    let provider = CountOverlapsProvider::new(
        Arc::new(ctx.clone()),
        "fbrain".to_string(),
        "exons".to_string(),
        exons_schema,
        vec![
            "contig".to_string(),
            "pos_start".to_string(),
            "pos_end".to_string(),
        ],
        vec![
            "contig".to_string(),
            "pos_start".to_string(),
            "pos_end".to_string(),
        ],
        FilterOp::Strict,
        true,
    );

    ctx.register_table("coverage_result", Arc::new(provider))?;
    let actual: Vec<RecordBatch> = ctx
        .sql("SELECT contig, pos_start, pos_end, coverage FROM coverage_result ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    // Load expected results generated by polars-bio
    let expected_table = format!(
        "CREATE EXTERNAL TABLE expected STORED AS PARQUET LOCATION '{EXPECTED_COVERAGE_PATH}'"
    );
    ctx.sql(expected_table.as_str()).await?;
    let expected: Vec<RecordBatch> = ctx
        .sql("SELECT contig, pos_start, pos_end, coverage FROM expected ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    // Compare row counts
    let actual_rows: usize = actual.iter().map(|b| b.num_rows()).sum();
    let expected_rows: usize = expected.iter().map(|b| b.num_rows()).sum();
    assert_eq!(actual_rows, expected_rows, "Row count mismatch");
    assert_eq!(actual_rows, 438_694, "Expected 438,694 rows");

    // Collect all coverage values for comparison
    let actual_coverage: Vec<i64> = actual
        .iter()
        .flat_map(|b| {
            b.column_by_name("coverage")
                .unwrap()
                .as_any()
                .downcast_ref::<Int64Array>()
                .unwrap()
                .values()
                .to_vec()
        })
        .collect();

    let expected_coverage: Vec<i64> = expected
        .iter()
        .flat_map(|b| {
            b.column_by_name("coverage")
                .unwrap()
                .as_any()
                .downcast_ref::<Int64Array>()
                .unwrap()
                .values()
                .to_vec()
        })
        .collect();

    assert_eq!(
        actual_coverage.len(),
        expected_coverage.len(),
        "Coverage vector length mismatch"
    );
    assert_eq!(
        actual_coverage, expected_coverage,
        "Coverage values mismatch"
    );

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
#[rstest::rstest]
async fn test_nearest_csv(ctx: SessionContext) -> Result<()> {
    init_ranges_tables(&ctx).await?;

    ctx.sql("SET bio.interval_join_algorithm TO CoitreesNearest")
        .await?;

    let query = r#"SELECT *
               FROM targets
               JOIN reads
               ON targets.contig = reads.contig
                  AND targets.pos_start <= reads.pos_end
                  AND targets.pos_end >= reads.pos_start
               ORDER BY targets.contig, targets.pos_start, targets.pos_end,
                        reads.contig, reads.pos_start, reads.pos_end"#;

    let result: Vec<RecordBatch> = ctx.sql(query).await?.collect().await?;

    // CoitreesNearest returns exactly 1 match per right-side (reads) row:
    // the nearest left-side (targets) interval. For overlapping reads,
    // the first overlapping target is returned. For chr3 read 234-300,
    // no target overlaps so the nearest target 100-200 is returned.
    let expected = [
        "+--------+-----------+---------+--------+-----------+---------+",
        "| contig | pos_start | pos_end | contig | pos_start | pos_end |",
        "+--------+-----------+---------+--------+-----------+---------+",
        "| chr1   | 100       | 190     | chr1   | 150       | 250     |",
        "| chr1   | 100       | 190     | chr1   | 190       | 300     |",
        "| chr1   | 10000     | 20000   | chr1   | 15000     | 15000   |",
        "| chr1   | 22100     | 22100   | chr1   | 22000     | 22300   |",
        "| chr1   | 400       | 600     | chr1   | 300       | 501     |",
        "| chr1   | 400       | 600     | chr1   | 500       | 700     |",
        "| chr2   | 100       | 190     | chr2   | 150       | 250     |",
        "| chr2   | 100       | 190     | chr2   | 190       | 300     |",
        "| chr2   | 10000     | 20000   | chr2   | 15000     | 15000   |",
        "| chr2   | 22100     | 22100   | chr2   | 22000     | 22300   |",
        "| chr2   | 400       | 600     | chr2   | 300       | 500     |",
        "| chr2   | 400       | 600     | chr2   | 500       | 700     |",
        "| chr3   | 100       | 200     | chr3   | 234       | 300     |",
        "+--------+-----------+---------+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_nearest_udtf_k1_matches_legacy() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    ctx.sql("SET bio.interval_join_algorithm TO CoitreesNearest")
        .await?;

    let legacy_query = r#"SELECT
               targets.contig AS left_contig,
               targets.pos_start AS left_pos_start,
               targets.pos_end AS left_pos_end,
               reads.contig AS right_contig,
               reads.pos_start AS right_pos_start,
               reads.pos_end AS right_pos_end
               FROM targets
               JOIN reads
               ON targets.contig = reads.contig
                  AND targets.pos_start <= reads.pos_end
                  AND targets.pos_end >= reads.pos_start
               ORDER BY right_contig, right_pos_start, right_pos_end, left_contig, left_pos_start, left_pos_end"#;

    let udtf_query = r#"SELECT
              left_contig, left_pos_start, left_pos_end,
              right_contig, right_pos_start, right_pos_end
              FROM nearest('targets', 'reads', 1, true, false)
              ORDER BY right_contig, right_pos_start, right_pos_end, left_contig, left_pos_start, left_pos_end"#;

    let legacy = ctx.sql(legacy_query).await?.collect().await?;
    let nearest = ctx.sql(udtf_query).await?.collect().await?;
    assert_eq!(
        pretty_format_batches(&legacy)?.to_string(),
        pretty_format_batches(&nearest)?.to_string()
    );

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_nearest_udtf_k2_overlap_false_and_null_match() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE l (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 10, 20),
        ('a', 30, 40),
        ('a', 50, 60)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE r (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 22, 22),
        ('a', 37, 37),
        ('b', 1, 1)
    "#,
    )
    .await?;

    let query = r#"SELECT *
            FROM nearest('l', 'r', 2, false)
            ORDER BY right_contig, right_pos_start, right_pos_end, left_pos_start NULLS LAST, left_contig"#;

    let result = ctx.sql(query).await?.collect().await?;

    let expected = [
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
        "| left_contig | left_pos_start | left_pos_end | right_contig | right_pos_start | right_pos_end | distance |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
        "| a           | 10             | 20           | a            | 22              | 22            | 2        |",
        "| a           | 30             | 40           | a            | 22              | 22            | 8        |",
        "| a           | 10             | 20           | a            | 37              | 37            | 17       |",
        "| a           | 50             | 60           | a            | 37              | 37            | 13       |",
        "|             |                |              | b            | 1               | 1             |          |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Bioframe parity ports (from polars-bio/tests/test_bioframe.py)
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_overlap_count() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql(
            r#"SELECT
                  reads.contig AS contig_1,
                  reads.pos_start AS pos_start_1,
                  reads.pos_end AS pos_end_1,
                  targets.contig AS contig_2,
                  targets.pos_start AS pos_start_2,
                  targets.pos_end AS pos_end_2
               FROM reads
               JOIN targets
                 ON reads.contig = targets.contig
                AND reads.pos_start <= targets.pos_end
                AND reads.pos_end >= targets.pos_start"#,
        )
        .await?
        .collect()
        .await?;

    // Ported from bioframe parity tests: 16 overlapping rows.
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 16);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_overlap_schema_rows() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql(
            r#"SELECT
                  reads.contig AS contig_1,
                  reads.pos_start AS pos_start_1,
                  reads.pos_end AS pos_end_1,
                  targets.contig AS contig_2,
                  targets.pos_start AS pos_start_2,
                  targets.pos_end AS pos_end_2
               FROM reads
               JOIN targets
                 ON reads.contig = targets.contig
                AND reads.pos_start <= targets.pos_end
                AND reads.pos_end >= targets.pos_start
               ORDER BY contig_1, pos_start_1, pos_end_1, contig_2, pos_start_2, pos_end_2"#,
        )
        .await?
        .collect()
        .await?;

    let expected = [
        "+----------+-------------+-----------+----------+-------------+-----------+",
        "| contig_1 | pos_start_1 | pos_end_1 | contig_2 | pos_start_2 | pos_end_2 |",
        "+----------+-------------+-----------+----------+-------------+-----------+",
        "| chr1     | 150         | 250       | chr1     | 100         | 190       |",
        "| chr1     | 150         | 250       | chr1     | 200         | 290       |",
        "| chr1     | 190         | 300       | chr1     | 100         | 190       |",
        "| chr1     | 190         | 300       | chr1     | 200         | 290       |",
        "| chr1     | 300         | 501       | chr1     | 400         | 600       |",
        "| chr1     | 500         | 700       | chr1     | 400         | 600       |",
        "| chr1     | 15000       | 15000     | chr1     | 10000       | 20000     |",
        "| chr1     | 22000       | 22300     | chr1     | 22100       | 22100     |",
        "| chr2     | 150         | 250       | chr2     | 100         | 190       |",
        "| chr2     | 150         | 250       | chr2     | 200         | 290       |",
        "| chr2     | 190         | 300       | chr2     | 100         | 190       |",
        "| chr2     | 190         | 300       | chr2     | 200         | 290       |",
        "| chr2     | 300         | 500       | chr2     | 400         | 600       |",
        "| chr2     | 500         | 700       | chr2     | 400         | 600       |",
        "| chr2     | 15000       | 15000     | chr2     | 10000       | 20000     |",
        "| chr2     | 22000       | 22300     | chr2     | 22100       | 22100     |",
        "+----------+-------------+-----------+----------+-------------+-----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_nearest_k1_count() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql(
            r#"SELECT
                   right_contig AS contig_1,
                   right_pos_start AS pos_start_1,
                   right_pos_end AS pos_end_1,
                   left_contig AS contig_2,
                   left_pos_start AS pos_start_2,
                   left_pos_end AS pos_end_2,
                   distance
               FROM nearest('reads', 'targets', 1, true)
               ORDER BY contig_1, pos_start_1, pos_end_1, distance, contig_2, pos_start_2, pos_end_2"#,
        )
        .await?
        .collect()
        .await?;

    // Ported from bioframe parity tests: one nearest row per target.
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 11);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_nearest_k1_schema_rows() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql(
            r#"SELECT
                   right_contig AS contig_1,
                   right_pos_start AS pos_start_1,
                   right_pos_end AS pos_end_1,
                   left_contig AS contig_2,
                   left_pos_start AS pos_start_2,
                   left_pos_end AS pos_end_2,
                   distance
               FROM nearest('reads', 'targets', 1, true)
               ORDER BY contig_1, pos_start_1, pos_end_1, distance, contig_2, pos_start_2, pos_end_2"#,
        )
        .await?
        .collect()
        .await?;

    // Expected values are ported from polars-bio _expected.PD_DF_NEAREST
    // (which is generated from bioframe.closest for k=1).
    let expected = [
        "+----------+-------------+-----------+----------+-------------+-----------+----------+",
        "| contig_1 | pos_start_1 | pos_end_1 | contig_2 | pos_start_2 | pos_end_2 | distance |",
        "+----------+-------------+-----------+----------+-------------+-----------+----------+",
        "| chr1     | 100         | 190       | chr1     | 150         | 250       | 0        |",
        "| chr1     | 200         | 290       | chr1     | 150         | 250       | 0        |",
        "| chr1     | 400         | 600       | chr1     | 300         | 501       | 0        |",
        "| chr1     | 10000       | 20000     | chr1     | 15000       | 15000     | 0        |",
        "| chr1     | 22100       | 22100     | chr1     | 22000       | 22300     | 0        |",
        "| chr2     | 100         | 190       | chr2     | 150         | 250       | 0        |",
        "| chr2     | 200         | 290       | chr2     | 150         | 250       | 0        |",
        "| chr2     | 400         | 600       | chr2     | 300         | 500       | 0        |",
        "| chr2     | 10000       | 20000     | chr2     | 15000       | 15000     | 0        |",
        "| chr2     | 22100       | 22100     | chr2     | 22000       | 22300     | 0        |",
        "| chr3     | 100         | 200       | chr3     | 234         | 300       | 34       |",
        "+----------+-------------+-----------+----------+-------------+-----------+----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_count_overlaps_count() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM count_overlaps('reads', 'targets')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 11);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_count_overlaps_schema_rows() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM count_overlaps('reads', 'targets') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    // Expected values are ported from polars-bio _expected.PD_DF_COUNT_OVERLAPS
    let expected = [
        "+--------+-----------+---------+-------+",
        "| contig | pos_start | pos_end | count |",
        "+--------+-----------+---------+-------+",
        "| chr1   | 100       | 190     | 2     |",
        "| chr1   | 200       | 290     | 2     |",
        "| chr1   | 400       | 600     | 2     |",
        "| chr1   | 10000     | 20000   | 1     |",
        "| chr1   | 22100     | 22100   | 1     |",
        "| chr2   | 100       | 190     | 2     |",
        "| chr2   | 200       | 290     | 2     |",
        "| chr2   | 400       | 600     | 2     |",
        "| chr2   | 10000     | 20000   | 1     |",
        "| chr2   | 22100     | 22100   | 1     |",
        "| chr3   | 100       | 200     | 0     |",
        "+--------+-----------+---------+-------+",
    ];
    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_coverage_count() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM coverage('reads', 'targets')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 11);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_coverage_schema_rows() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM coverage('reads', 'targets') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+----------+",
        "| contig | pos_start | pos_end | coverage |",
        "+--------+-----------+---------+----------+",
        "| chr1   | 100       | 190     | 41       |",
        "| chr1   | 200       | 290     | 92       |",
        "| chr1   | 400       | 600     | 202      |",
        "| chr1   | 10000     | 20000   | 1        |",
        "| chr1   | 22100     | 22100   | 2        |",
        "| chr2   | 100       | 190     | 41       |",
        "| chr2   | 200       | 290     | 92       |",
        "| chr2   | 400       | 600     | 202      |",
        "| chr2   | 10000     | 20000   | 1        |",
        "| chr2   | 22100     | 22100   | 2        |",
        "| chr3   | 100       | 200     | 0        |",
        "+--------+-----------+---------+----------+",
    ];
    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_count_overlaps_udtf_strict_zero_based_boundary() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE reads (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 190, 300)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE targets (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 190)
    "#,
    )
    .await?;

    let weak = ctx
        .sql("SELECT * FROM count_overlaps('reads', 'targets')")
        .await?
        .collect()
        .await?;
    let strict = ctx
        .sql("SELECT * FROM count_overlaps('reads', 'targets', 'contig', 'pos_start', 'pos_end', 'strict')")
        .await?
        .collect()
        .await?;

    let weak_expected = [
        "+--------+-----------+---------+-------+",
        "| contig | pos_start | pos_end | count |",
        "+--------+-----------+---------+-------+",
        "| a      | 100       | 190     | 1     |",
        "+--------+-----------+---------+-------+",
    ];
    let strict_expected = [
        "+--------+-----------+---------+-------+",
        "| contig | pos_start | pos_end | count |",
        "+--------+-----------+---------+-------+",
        "| a      | 100       | 190     | 0     |",
        "+--------+-----------+---------+-------+",
    ];

    assert_batches_sorted_eq!(weak_expected, &weak);
    assert_batches_sorted_eq!(strict_expected, &strict);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_nearest_udtf_strict_zero_based_boundary_distance() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE reads (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 190, 190)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE targets (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 190)
    "#,
    )
    .await?;

    let result = ctx
        .sql(
            r#"SELECT
                   left_contig,
                   left_pos_start,
                   left_pos_end,
                   right_contig,
                   right_pos_start,
                   right_pos_end,
                   distance
               FROM nearest('reads', 'targets', 1, true, true, 'contig', 'pos_start', 'pos_end', 'strict')"#,
        )
        .await?
        .collect()
        .await?;

    let expected = [
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
        "| left_contig | left_pos_start | left_pos_end | right_contig | right_pos_start | right_pos_end | distance |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
        "| a           | 190            | 190          | a            | 100             | 190           | 1        |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_nearest_udtf_empty_left_emits_null_rows() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE l AS
        SELECT
            CAST('a' AS TEXT) AS contig,
            CAST(1 AS INTEGER) AS pos_start,
            CAST(2 AS INTEGER) AS pos_end
        WHERE FALSE
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE r (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 110),
        ('b', 200, 210)
    "#,
    )
    .await?;

    let result = ctx
        .sql(
            r#"SELECT *
               FROM nearest('l', 'r')
               ORDER BY right_contig, right_pos_start"#,
        )
        .await?
        .collect()
        .await?;

    let expected = [
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
        "| left_contig | left_pos_start | left_pos_end | right_contig | right_pos_start | right_pos_end | distance |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
        "|             |                |              | a            | 100             | 110           |          |",
        "|             |                |              | b            | 200             | 210           |          |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_nearest_udtf_compute_distance_false() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE l (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 10, 20),
        ('a', 30, 40)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE r (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 22, 22)
    "#,
    )
    .await?;

    let result = ctx
        .sql(
            r#"SELECT *
               FROM nearest('l', 'r', 1, true, false)
               ORDER BY left_pos_start"#,
        )
        .await?
        .collect()
        .await?;

    // No distance column when compute_distance=false
    let expected = [
        "+-------------+----------------+--------------+--------------+-----------------+---------------+",
        "| left_contig | left_pos_start | left_pos_end | right_contig | right_pos_start | right_pos_end |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+",
        "| a           | 10             | 20           | a            | 22              | 22            |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Table function (UDTF) tests for coverage() and count_overlaps()
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_count_overlaps_udtf() -> Result<()> {
    let ctx = create_bio_session();

    let reads = format!(
        "CREATE EXTERNAL TABLE reads STORED AS CSV LOCATION '{RANGES_READS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(reads.as_str()).await?;

    let targets = format!(
        "CREATE EXTERNAL TABLE targets STORED AS CSV LOCATION '{RANGES_TARGETS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(targets.as_str()).await?;

    let result: Vec<RecordBatch> = ctx
        .sql("SELECT * FROM count_overlaps('reads', 'targets') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------+",
        "| contig | pos_start | pos_end | count |",
        "+--------+-----------+---------+-------+",
        "| chr1   | 100       | 190     | 2     |",
        "| chr1   | 200       | 290     | 2     |",
        "| chr1   | 400       | 600     | 2     |",
        "| chr1   | 10000     | 20000   | 1     |",
        "| chr1   | 22100     | 22100   | 1     |",
        "| chr2   | 100       | 190     | 2     |",
        "| chr2   | 200       | 290     | 2     |",
        "| chr2   | 400       | 600     | 2     |",
        "| chr2   | 10000     | 20000   | 1     |",
        "| chr2   | 22100     | 22100   | 1     |",
        "| chr3   | 100       | 200     | 0     |",
        "+--------+-----------+---------+-------+",
    ];

    assert_batches_sorted_eq!(expected, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_coverage_udtf() -> Result<()> {
    let ctx = create_bio_session();

    let reads = format!(
        "CREATE EXTERNAL TABLE reads STORED AS CSV LOCATION '{RANGES_READS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(reads.as_str()).await?;

    let targets = format!(
        "CREATE EXTERNAL TABLE targets STORED AS CSV LOCATION '{RANGES_TARGETS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(targets.as_str()).await?;

    let result: Vec<RecordBatch> = ctx
        .sql("SELECT * FROM coverage('reads', 'targets') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+----------+",
        "| contig | pos_start | pos_end | coverage |",
        "+--------+-----------+---------+----------+",
        "| chr1   | 100       | 190     | 41       |",
        "| chr1   | 200       | 290     | 92       |",
        "| chr1   | 400       | 600     | 202      |",
        "| chr1   | 10000     | 20000   | 1        |",
        "| chr1   | 22100     | 22100   | 2        |",
        "| chr2   | 100       | 190     | 41       |",
        "| chr2   | 200       | 290     | 92       |",
        "| chr2   | 400       | 600     | 202      |",
        "| chr2   | 10000     | 20000   | 1        |",
        "| chr2   | 22100     | 22100   | 2        |",
        "| chr3   | 100       | 200     | 0        |",
        "+--------+-----------+---------+----------+",
    ];

    assert_batches_sorted_eq!(expected, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_register_count_overlaps_on_existing_context() -> Result<()> {
    // Test that register_ranges_functions works on an existing bio context
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false);
    let ctx = SessionContext::new_with_bio(config);
    register_ranges_functions(&ctx);

    let reads = format!(
        "CREATE EXTERNAL TABLE reads STORED AS CSV LOCATION '{RANGES_READS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(reads.as_str()).await?;

    let targets = format!(
        "CREATE EXTERNAL TABLE targets STORED AS CSV LOCATION '{RANGES_TARGETS_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(targets.as_str()).await?;

    let result: Vec<RecordBatch> = ctx
        .sql("SELECT * FROM count_overlaps('reads', 'targets') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let actual_rows: usize = result.iter().map(|b| b.num_rows()).sum();
    assert_eq!(actual_rows, 11);

    Ok(())
}
