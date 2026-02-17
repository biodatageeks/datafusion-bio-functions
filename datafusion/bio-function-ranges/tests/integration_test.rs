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

const MERGE_INPUT_PATH: &str = "../../testing/data/merge/input.csv";

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

fn create_bio_session_with_target_partitions(target_partitions: usize) -> SessionContext {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false)
        .with_target_partitions(target_partitions);
    let ctx = SessionContext::new_with_bio(config);
    register_ranges_functions(&ctx);
    ctx
}

async fn collect_udtf_query_with_partitions(
    target_partitions: usize,
    query: &str,
) -> Result<Vec<RecordBatch>> {
    let ctx = create_bio_session_with_target_partitions(target_partitions);
    init_ranges_tables(&ctx).await?;
    ctx.sql(query).await?.collect().await
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
        "| a           | 190            | 190          | a            | 100             | 190           | 0        |",
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

// ─────────────────────────────────────────────────────────────────────────────
// Overlap UDTF tests
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_overlap_udtf_count() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM overlap('reads', 'targets')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 16);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_overlap_udtf_schema_rows() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql(
            r#"SELECT * FROM overlap('reads', 'targets')
               ORDER BY left_contig, left_pos_start, left_pos_end,
                        right_contig, right_pos_start, right_pos_end"#,
        )
        .await?
        .collect()
        .await?;

    let expected = [
        "+-------------+----------------+--------------+--------------+-----------------+---------------+",
        "| left_contig | left_pos_start | left_pos_end | right_contig | right_pos_start | right_pos_end |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+",
        "| chr1        | 150            | 250          | chr1         | 100             | 190           |",
        "| chr1        | 150            | 250          | chr1         | 200             | 290           |",
        "| chr1        | 190            | 300          | chr1         | 100             | 190           |",
        "| chr1        | 190            | 300          | chr1         | 200             | 290           |",
        "| chr1        | 300            | 501          | chr1         | 400             | 600           |",
        "| chr1        | 500            | 700          | chr1         | 400             | 600           |",
        "| chr1        | 15000          | 15000        | chr1         | 10000           | 20000         |",
        "| chr1        | 22000          | 22300        | chr1         | 22100           | 22100         |",
        "| chr2        | 150            | 250          | chr2         | 100             | 190           |",
        "| chr2        | 150            | 250          | chr2         | 200             | 290           |",
        "| chr2        | 190            | 300          | chr2         | 100             | 190           |",
        "| chr2        | 190            | 300          | chr2         | 200             | 290           |",
        "| chr2        | 300            | 500          | chr2         | 400             | 600           |",
        "| chr2        | 500            | 700          | chr2         | 400             | 600           |",
        "| chr2        | 15000          | 15000        | chr2         | 10000           | 20000         |",
        "| chr2        | 22000          | 22300        | chr2         | 22100           | 22100         |",
        "+-------------+----------------+--------------+--------------+-----------------+---------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_parity_with_join() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let join_query = r#"SELECT
          reads.contig AS left_contig,
          reads.pos_start AS left_pos_start,
          reads.pos_end AS left_pos_end,
          targets.contig AS right_contig,
          targets.pos_start AS right_pos_start,
          targets.pos_end AS right_pos_end
       FROM reads
       JOIN targets
         ON reads.contig = targets.contig
        AND reads.pos_start <= targets.pos_end
        AND reads.pos_end >= targets.pos_start
       ORDER BY left_contig, left_pos_start, left_pos_end,
                right_contig, right_pos_start, right_pos_end"#;

    let udtf_query = r#"SELECT
          left_contig, left_pos_start, left_pos_end,
          right_contig, right_pos_start, right_pos_end
       FROM overlap('reads', 'targets')
       ORDER BY left_contig, left_pos_start, left_pos_end,
                right_contig, right_pos_start, right_pos_end"#;

    let join_result = ctx.sql(join_query).await?.collect().await?;
    let udtf_result = ctx.sql(udtf_query).await?.collect().await?;

    assert_eq!(
        pretty_format_batches(&join_result)?.to_string(),
        pretty_format_batches(&udtf_result)?.to_string()
    );
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_strict_boundary() -> Result<()> {
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

    // Weak: touching intervals overlap (190 <= 190)
    let weak = ctx
        .sql("SELECT * FROM overlap('reads', 'targets')")
        .await?
        .collect()
        .await?;
    assert_eq!(weak.iter().map(|b| b.num_rows()).sum::<usize>(), 1);

    // Strict: touching intervals do NOT overlap (190 < 190 is false)
    let strict = ctx
        .sql(
            "SELECT * FROM overlap('reads', 'targets', 'contig', 'pos_start', 'pos_end', 'strict')",
        )
        .await?
        .collect()
        .await?;
    assert_eq!(strict.iter().map(|b| b.num_rows()).sum::<usize>(), 0);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_adjacent_zero_based_no_overlap() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE a (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE b (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 200, 300)
    "#,
    )
    .await?;

    // Strict (0-based half-open): [100,200) and [200,300) don't overlap
    let result = ctx
        .sql("SELECT * FROM overlap('a', 'b', 'contig', 'pos_start', 'pos_end', 'strict')")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_adjacent_one_based_overlap() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE a (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE b (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 200, 300)
    "#,
    )
    .await?;

    // Weak (1-based closed): [100,200] and [200,300] share position 200
    let result = ctx
        .sql("SELECT * FROM overlap('a', 'b')")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 1);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_same_interval() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE a (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE b (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM overlap('a', 'b')")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 1);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_contained() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE a (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE b (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 150, 180)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM overlap('a', 'b')")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 1);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_no_matches() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE a (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE b (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 300, 400),
        ('b', 100, 200)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM overlap('a', 'b')")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_overlap_udtf_custom_columns() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE a (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;
    ctx.sql(
        r#"
        CREATE TABLE b (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 150, 250)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM overlap('a', 'b', 'chr', 's', 'e')")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 1);

    // Verify left_/right_ prefixed column names
    let schema = result[0].schema();
    assert_eq!(schema.field(0).name(), "left_chr");
    assert_eq!(schema.field(3).name(), "right_chr");

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Merge UDTF tests
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_merge_udtf_count() -> Result<()> {
    let ctx = create_bio_session();

    let create = format!(
        "CREATE EXTERNAL TABLE intervals STORED AS CSV LOCATION '{MERGE_INPUT_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(create.as_str()).await?;

    // bioframe uses 0-based half-open coordinates by default → 'strict'
    let result = ctx
        .sql("SELECT * FROM merge('intervals', 0, 'contig', 'pos_start', 'pos_end', 'strict')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 8);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_bioframe_merge_udtf_schema_rows() -> Result<()> {
    let ctx = create_bio_session();

    let create = format!(
        "CREATE EXTERNAL TABLE intervals STORED AS CSV LOCATION '{MERGE_INPUT_PATH}' OPTIONS ('has_header' 'true')"
    );
    ctx.sql(create.as_str()).await?;

    // bioframe uses 0-based half-open coordinates by default → 'strict'
    let result = ctx
        .sql("SELECT * FROM merge('intervals', 0, 'contig', 'pos_start', 'pos_end', 'strict') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| chr1   | 100       | 300     | 4           |",
        "| chr1   | 300       | 700     | 3           |",
        "| chr1   | 10000     | 20000   | 2           |",
        "| chr1   | 22000     | 22300   | 2           |",
        "| chr2   | 100       | 300     | 4           |",
        "| chr2   | 300       | 700     | 3           |",
        "| chr2   | 10000     | 20000   | 2           |",
        "| chr2   | 22000     | 22300   | 2           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_basic() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 150, 250),
        ('a', 300, 400)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM merge('intervals') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| a      | 100       | 250     | 2           |",
        "| a      | 300       | 400     | 1           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_min_dist() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 201, 300)
    "#,
    )
    .await?;

    // min_dist=0: 201 <= 200 + 0 is false → separate
    let result = ctx
        .sql("SELECT * FROM merge('intervals', 0) ORDER BY pos_start")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 2);

    // min_dist=1: 201 <= 200 + 1 is true → merge
    let result = ctx
        .sql("SELECT * FROM merge('intervals', 1) ORDER BY pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| a      | 100       | 300     | 2           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_all_overlap() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 300),
        ('a', 150, 250),
        ('a', 200, 400)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM merge('intervals')")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| a      | 100       | 400     | 3           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_single_row() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM merge('intervals')")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| a      | 100       | 200     | 1           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_empty() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals AS
        SELECT
            CAST('a' AS TEXT) AS contig,
            CAST(1 AS INTEGER) AS pos_start,
            CAST(2 AS INTEGER) AS pos_end
        WHERE FALSE
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM merge('intervals')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_custom_columns() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 150, 250)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM merge('intervals', 0, 'chr', 's', 'e')")
        .await?
        .collect()
        .await?;

    let expected = [
        "+-----+-----+-----+-------------+",
        "| chr | s   | e   | n_intervals |",
        "+-----+-----+-----+-------------+",
        "| a   | 100 | 250 | 2           |",
        "+-----+-----+-----+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);

    // Verify column names come from args
    let schema = result[0].schema();
    assert_eq!(schema.field(0).name(), "chr");
    assert_eq!(schema.field(1).name(), "s");
    assert_eq!(schema.field(2).name(), "e");

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_reads_csv() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM merge('reads') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| chr1   | 150       | 700     | 4           |",
        "| chr1   | 15000     | 15000   | 1           |",
        "| chr1   | 22000     | 22300   | 1           |",
        "| chr2   | 150       | 700     | 4           |",
        "| chr2   | 15000     | 15000   | 1           |",
        "| chr2   | 22000     | 22300   | 1           |",
        "| chr3   | 234       | 300     | 1           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_adjacent_strict() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 150),
        ('a', 150, 200)
    "#,
    )
    .await?;

    // Strict: 150 < 150 is false → should NOT merge
    let result = ctx
        .sql("SELECT * FROM merge('intervals', 0, 'contig', 'pos_start', 'pos_end', 'strict') ORDER BY pos_start")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 2);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_merge_udtf_adjacent_weak() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 150),
        ('a', 150, 200)
    "#,
    )
    .await?;

    // Weak: 150 <= 150 is true → SHOULD merge
    let result = ctx
        .sql("SELECT * FROM merge('intervals')")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+-------------+",
        "| contig | pos_start | pos_end | n_intervals |",
        "+--------+-----------+---------+-------------+",
        "| a      | 100       | 200     | 2           |",
        "+--------+-----------+---------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Cluster UDTF tests
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_basic() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 150, 250),
        ('a', 400, 500)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM cluster('intervals') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 200     | 0       | 100           | 250         |",
        "| a      | 150       | 250     | 0       | 100           | 250         |",
        "| a      | 400       | 500     | 1       | 400           | 500         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_min_dist() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 210, 300)
    "#,
    )
    .await?;

    // min_dist=0: gap of 10 → separate clusters
    let result = ctx
        .sql("SELECT * FROM cluster('intervals', 0) ORDER BY pos_start")
        .await?
        .collect()
        .await?;
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 2);

    // Check that they are in different clusters
    let expected_separate = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 200     | 0       | 100           | 200         |",
        "| a      | 210       | 300     | 1       | 210           | 300         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];
    assert_batches_sorted_eq!(expected_separate, &result);

    // min_dist=10: 210 <= 200 + 10 → same cluster
    let result = ctx
        .sql("SELECT * FROM cluster('intervals', 10) ORDER BY pos_start")
        .await?
        .collect()
        .await?;

    let expected_merged = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 200     | 0       | 100           | 300         |",
        "| a      | 210       | 300     | 0       | 100           | 300         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];
    assert_batches_sorted_eq!(expected_merged, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_multi_contig() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 150, 250),
        ('b', 100, 200),
        ('b', 300, 400)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM cluster('intervals') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 200     | 0       | 100           | 250         |",
        "| a      | 150       | 250     | 0       | 100           | 250         |",
        "| b      | 100       | 200     | 1       | 100           | 200         |",
        "| b      | 300       | 400     | 2       | 300           | 400         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_single_row() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM cluster('intervals')")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 200     | 0       | 100           | 200         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_empty() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals AS
        SELECT
            CAST('a' AS TEXT) AS contig,
            CAST(1 AS INTEGER) AS pos_start,
            CAST(2 AS INTEGER) AS pos_end
        WHERE FALSE
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM cluster('intervals')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_strict() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 150),
        ('a', 150, 200)
    "#,
    )
    .await?;

    // Weak: 150 <= 150 → same cluster
    let result = ctx
        .sql("SELECT * FROM cluster('intervals')")
        .await?
        .collect()
        .await?;

    let expected_weak = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 150     | 0       | 100           | 200         |",
        "| a      | 150       | 200     | 0       | 100           | 200         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];
    assert_batches_sorted_eq!(expected_weak, &result);

    // Strict: 150 < 150 is false → different clusters
    let result = ctx
        .sql("SELECT * FROM cluster('intervals', 0, 'contig', 'pos_start', 'pos_end', 'strict')")
        .await?
        .collect()
        .await?;

    let expected_strict = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| a      | 100       | 150     | 0       | 100           | 150         |",
        "| a      | 150       | 200     | 1       | 150           | 200         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];
    assert_batches_sorted_eq!(expected_strict, &result);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_large_min_dist_no_overflow() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start BIGINT, pos_end BIGINT) AS VALUES
        ('a', CAST(9223372036854775800 AS BIGINT), CAST(9223372036854775806 AS BIGINT)),
        ('a', CAST(9223372036854775807 AS BIGINT), CAST(9223372036854775807 AS BIGINT))
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM cluster('intervals', 100) ORDER BY pos_start")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 2);

    let batch = &result[0];

    let cluster = batch
        .column_by_name("cluster")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();
    let cluster_start = batch
        .column_by_name("cluster_start")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();
    let cluster_end = batch
        .column_by_name("cluster_end")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();

    assert_eq!(cluster.value(0), 0);
    assert_eq!(cluster.value(1), 0);
    assert_eq!(cluster_start.value(0), 9_223_372_036_854_775_800);
    assert_eq!(cluster_start.value(1), 9_223_372_036_854_775_800);
    assert_eq!(cluster_end.value(0), i64::MAX);
    assert_eq!(cluster_end.value(1), i64::MAX);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_custom_columns() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 150, 250)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM cluster('intervals', 0, 'chr', 's', 'e')")
        .await?
        .collect()
        .await?;

    // Verify column names come from args
    let schema = result[0].schema();
    assert_eq!(schema.field(0).name(), "chr");
    assert_eq!(schema.field(1).name(), "s");
    assert_eq!(schema.field(2).name(), "e");
    assert_eq!(schema.field(3).name(), "cluster");

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 2);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_cluster_udtf_reads_csv() -> Result<()> {
    let ctx = create_bio_session();
    init_ranges_tables(&ctx).await?;

    let result = ctx
        .sql("SELECT * FROM cluster('reads') ORDER BY contig, pos_start, pos_end")
        .await?
        .collect()
        .await?;

    // reads.csv has:
    // chr1: [150,250],[190,300],[300,501],[500,700] → cluster 0 (all overlap/touch)
    //       [15000,15000] → cluster 1
    //       [22000,22300] → cluster 2
    // chr2: same pattern → clusters 3,4,5
    // chr3: [234,300] → cluster 6
    let expected = [
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| contig | pos_start | pos_end | cluster | cluster_start | cluster_end |",
        "+--------+-----------+---------+---------+---------------+-------------+",
        "| chr1   | 150       | 250     | 0       | 150           | 700         |",
        "| chr1   | 190       | 300     | 0       | 150           | 700         |",
        "| chr1   | 300       | 501     | 0       | 150           | 700         |",
        "| chr1   | 500       | 700     | 0       | 150           | 700         |",
        "| chr1   | 15000     | 15000   | 1       | 15000         | 15000       |",
        "| chr1   | 22000     | 22300   | 2       | 22000         | 22300       |",
        "| chr2   | 150       | 250     | 3       | 150           | 700         |",
        "| chr2   | 190       | 300     | 3       | 150           | 700         |",
        "| chr2   | 300       | 500     | 3       | 150           | 700         |",
        "| chr2   | 500       | 700     | 3       | 150           | 700         |",
        "| chr2   | 15000     | 15000   | 4       | 15000         | 15000       |",
        "| chr2   | 22000     | 22300   | 5       | 22000         | 22300       |",
        "| chr3   | 234       | 300     | 6       | 234           | 300         |",
        "+--------+-----------+---------+---------+---------------+-------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Complement UDTF tests
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_basic_no_view() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 300, 400)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    // Without view: view is [0, i64::MAX) per contig
    // Gaps: [0, 100), [200, 300), [400, i64::MAX)
    let expected = [
        "+--------+-----------+---------------------+",
        "| contig | pos_start | pos_end             |",
        "+--------+-----------+---------------------+",
        "| a      | 0         | 100                 |",
        "| a      | 200       | 300                 |",
        "| a      | 400       | 9223372036854775807 |",
        "+--------+-----------+---------------------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_with_view() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('a', 300, 400)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE chromsizes (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 0         | 100     |",
        "| a      | 200       | 300     |",
        "| a      | 400       | 500     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_no_gaps() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE chromsizes (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_empty_input() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals AS
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
        CREATE TABLE chromsizes (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    // Empty input → entire view is a gap
    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 0         | 500     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_multi_contig() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200),
        ('b', 300, 400)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE chromsizes (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500),
        ('b', 0, 600)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 0         | 100     |",
        "| a      | 200       | 500     |",
        "| b      | 0         | 300     |",
        "| b      | 400       | 600     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_overlapping_input() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 250),
        ('a', 200, 400)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE chromsizes (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    // Merged input: [100, 400], gaps: [0, 100), [400, 500)
    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 0         | 100     |",
        "| a      | 400       | 500     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_custom_columns() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE chromsizes (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 0, 500)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes', 'chr', 's', 'e') ORDER BY chr, s")
        .await?
        .collect()
        .await?;

    // Verify column names come from input cols
    let schema = result[0].schema();
    assert_eq!(schema.field(0).name(), "chr");
    assert_eq!(schema.field(1).name(), "s");
    assert_eq!(schema.field(2).name(), "e");

    let expected = [
        "+-----+-----+-----+",
        "| chr | s   | e   |",
        "+-----+-----+-----+",
        "| a   | 0   | 100 |",
        "| a   | 200 | 500 |",
        "+-----+-----+-----+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_complement_udtf_view_contig_no_input() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE intervals (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE chromsizes (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 0, 500),
        ('b', 0, 300)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM complement('intervals', 'chromsizes') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    // 'a' has interval [100,200] → gaps [0,100), [200,500)
    // 'b' has no intervals → entire view [0,300) is gap
    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 0         | 100     |",
        "| a      | 200       | 500     |",
        "| b      | 0         | 300     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Subtract UDTF tests
// ─────────────────────────────────────────────────────────────────────────────

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_basic() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 400)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 200, 300)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 200     |",
        "| a      | 300       | 400     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_no_overlap() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 300, 400)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 200     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_complete_removal() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 50, 300)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t')")
        .await?
        .collect()
        .await?;

    // Left [100,200) entirely covered by right [50,300)
    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_multiple_right() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 500)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 150, 200),
        ('a', 300, 350)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 150     |",
        "| a      | 200       | 300     |",
        "| a      | 350       | 500     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_multi_contig() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 300),
        ('b', 100, 300)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 150, 250)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    // 'a': [100,150), [250,300)
    // 'b': [100,300) (no right intervals for 'b')
    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 150     |",
        "| a      | 250       | 300     |",
        "| b      | 100       | 300     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_empty_left() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t AS
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
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t')")
        .await?
        .collect()
        .await?;

    assert_eq!(result.iter().map(|b| b.num_rows()).sum::<usize>(), 0);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_empty_right() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 200)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t AS
        SELECT
            CAST('a' AS TEXT) AS contig,
            CAST(1 AS INTEGER) AS pos_start,
            CAST(2 AS INTEGER) AS pos_end
        WHERE FALSE
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 200     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_custom_columns() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 100, 400)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (chr TEXT, s INTEGER, e INTEGER) AS VALUES
        ('a', 200, 300)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t', 'chr', 's', 'e') ORDER BY chr, s")
        .await?
        .collect()
        .await?;

    let schema = result[0].schema();
    assert_eq!(schema.field(0).name(), "chr");
    assert_eq!(schema.field(1).name(), "s");
    assert_eq!(schema.field(2).name(), "e");

    let expected = [
        "+-----+-----+-----+",
        "| chr | s   | e   |",
        "+-----+-----+-----+",
        "| a   | 100 | 200 |",
        "| a   | 300 | 400 |",
        "+-----+-----+-----+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_strict_boundary() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 300)
    "#,
    )
    .await?;

    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 300, 400)
    "#,
    )
    .await?;

    // Weak: right [300,400) overlaps left [100,300) at boundary → subtract at 300
    let weak = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY pos_start")
        .await?
        .collect()
        .await?;

    let expected_weak = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 300     |",
        "+--------+-----------+---------+",
    ];
    assert_batches_sorted_eq!(expected_weak, &weak);

    // Strict: right [300,400) does NOT overlap left [100,300) (300 >= 300 → no overlap)
    let strict = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t', 'contig', 'pos_start', 'pos_end', 'strict') ORDER BY pos_start")
        .await?
        .collect()
        .await?;

    let expected_strict = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 300     |",
        "+--------+-----------+---------+",
    ];
    assert_batches_sorted_eq!(expected_strict, &strict);

    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_subtract_udtf_overlapping_right() -> Result<()> {
    let ctx = create_bio_session();

    ctx.sql(
        r#"
        CREATE TABLE left_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 100, 500)
    "#,
    )
    .await?;

    // Overlapping right intervals
    ctx.sql(
        r#"
        CREATE TABLE right_t (contig TEXT, pos_start INTEGER, pos_end INTEGER) AS VALUES
        ('a', 150, 250),
        ('a', 200, 350)
    "#,
    )
    .await?;

    let result = ctx
        .sql("SELECT * FROM subtract('left_t', 'right_t') ORDER BY contig, pos_start")
        .await?
        .collect()
        .await?;

    // Right covers [150,350) after considering overlaps
    let expected = [
        "+--------+-----------+---------+",
        "| contig | pos_start | pos_end |",
        "+--------+-----------+---------+",
        "| a      | 100       | 150     |",
        "| a      | 350       | 500     |",
        "+--------+-----------+---------+",
    ];

    assert_batches_sorted_eq!(expected, &result);
    Ok(())
}

#[tokio::test(flavor = "multi_thread")]
async fn test_range_udtfs_target_partitions_invariant() -> Result<()> {
    let cases = [
        (
            "merge",
            "SELECT * FROM merge('reads') ORDER BY contig, pos_start, pos_end, n_intervals",
        ),
        (
            "cluster",
            "SELECT * FROM cluster('reads') ORDER BY contig, pos_start, pos_end, cluster",
        ),
        (
            "complement",
            "SELECT * FROM complement('reads', 'targets') ORDER BY contig, pos_start, pos_end",
        ),
        (
            "subtract",
            "SELECT * FROM subtract('reads', 'targets') ORDER BY contig, pos_start, pos_end",
        ),
    ];

    for (name, query) in cases {
        let result_1 = collect_udtf_query_with_partitions(1, query).await?;
        let result_4 = collect_udtf_query_with_partitions(4, query).await?;

        let formatted_1 = pretty_format_batches(&result_1)?.to_string();
        let formatted_4 = pretty_format_batches(&result_4)?.to_string();
        assert_eq!(
            formatted_1, formatted_4,
            "partition-invariance failed for {name}"
        );
    }

    Ok(())
}
