use std::sync::Arc;

use datafusion::arrow::array::{Int16Array, Int32Array, StringArray};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use futures::StreamExt;

use datafusion_bio_function_pileup::physical_exec::{PileupConfig, PileupExec};
use datafusion_bio_function_pileup::register_pileup_functions;

/// Helper to collect coverage results from a PileupExec over all partitions.
/// Uses `zero_based: true` so positions match 0-based test expectations.
async fn collect_coverage(ctx: &SessionContext, table_name: &str) -> Vec<(String, i32, i32, i16)> {
    let df = ctx.table(table_name).await.unwrap();
    let plan = df.create_physical_plan().await.unwrap();
    let pileup = PileupExec::new(
        plan,
        PileupConfig {
            zero_based: true,
            ..PileupConfig::default()
        },
    );

    let task_ctx = ctx.task_ctx();
    let num_partitions = pileup.properties().partitioning.partition_count();

    let mut all_rows = Vec::new();
    for partition in 0..num_partitions {
        let stream = pileup.execute(partition, task_ctx.clone()).unwrap();
        let batches: Vec<RecordBatch> = stream
            .collect::<Vec<_>>()
            .await
            .into_iter()
            .filter_map(|r| r.ok())
            .collect();

        for batch in batches {
            let contigs = batch
                .column(0)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let starts = batch
                .column(1)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap();
            let ends = batch
                .column(2)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap();
            let covs = batch
                .column(3)
                .as_any()
                .downcast_ref::<Int16Array>()
                .unwrap();

            for i in 0..batch.num_rows() {
                all_rows.push((
                    contigs.value(i).to_string(),
                    starts.value(i),
                    ends.value(i),
                    covs.value(i),
                ));
            }
        }
    }
    all_rows
}

/// Helper to collect coverage results from SQL UDTF.
async fn collect_coverage_sql(ctx: &SessionContext, sql: &str) -> Vec<(String, i32, i32, i16)> {
    let df = ctx.sql(sql).await.unwrap();
    let batches = df.collect().await.unwrap();

    let mut all_rows = Vec::new();
    for batch in batches {
        let contigs = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let starts = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let ends = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let covs = batch
            .column(3)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        for i in 0..batch.num_rows() {
            all_rows.push((
                contigs.value(i).to_string(),
                starts.value(i),
                ends.value(i),
                covs.value(i),
            ));
        }
    }
    all_rows
}

/// Matches mosdepth `overlapFastMode` test.
///
/// Our algorithm does NOT do mate-pair overlap deduplication, so it matches
/// mosdepth's `--fast-mode` behavior.
///
/// Mosdepth fast mode output (BED 0-based half-open): MT 0-6→1, 6-42→2, 42-80→1
/// Our output (0-based closed): MT 0-5→1, 6-41→2, 42-79→1
#[tokio::test]
async fn test_ovl_fast_mode() {
    let bam_path = format!("{}/tests/data/ovl.bam", env!("CARGO_MANIFEST_DIR"));
    let table = BamTableProvider::new(bam_path, None, true, None)
        .await
        .unwrap();

    let ctx = SessionContext::new();
    ctx.register_table("reads", Arc::new(table)).unwrap();

    let rows = collect_coverage(&ctx, "reads").await;

    // Filter to MT contig
    let mt_rows: Vec<_> = rows.iter().filter(|r| r.0 == "MT").collect();

    assert_eq!(
        mt_rows.len(),
        3,
        "Expected 3 coverage blocks on MT, got {}: {:?}",
        mt_rows.len(),
        mt_rows
    );

    assert_eq!(mt_rows[0], &("MT".to_string(), 0, 5, 1));
    assert_eq!(mt_rows[1], &("MT".to_string(), 6, 41, 2));
    assert_eq!(mt_rows[2], &("MT".to_string(), 42, 79, 1));
}

/// Tests overlapping mate pairs (no overlap deduplication).
///
/// The BAM file has two reads (a mate pair) both aligned to chr "1" at position
/// 565173 with 80M CIGAR. Since our algorithm is fast-mode (no mate-pair overlap
/// dedup), both reads contribute coverage, giving coverage=2.
///
/// Mosdepth default mode deduplicates overlapping mates, yielding coverage=1.
/// Our output (0-based closed): 1:565173-565252 → coverage=2
#[tokio::test]
async fn test_overlapping_pairs() {
    let bam_path = format!(
        "{}/tests/data/overlapping-pairs.bam",
        env!("CARGO_MANIFEST_DIR")
    );
    let table = BamTableProvider::new(bam_path, None, true, None)
        .await
        .unwrap();

    let ctx = SessionContext::new();
    ctx.register_table("reads", Arc::new(table)).unwrap();

    let rows = collect_coverage(&ctx, "reads").await;

    // BAM uses "1" not "chr1" as the reference name
    let chr1_rows: Vec<_> = rows.iter().filter(|r| r.0 == "1").collect();

    assert_eq!(
        chr1_rows.len(),
        1,
        "Expected 1 coverage block on 1, got {}: {:?}",
        chr1_rows.len(),
        chr1_rows
    );

    // Both mates align to same region — coverage=2 (no overlap dedup)
    assert_eq!(chr1_rows[0], &("1".to_string(), 565173, 565252, 2));
}

/// Test coverage via SQL UDTF: SELECT * FROM depth('path/to/file.bam', true)
/// Uses `true` (0-based) to match mosdepth 0-based expectations.
#[tokio::test(flavor = "multi_thread")]
async fn test_ovl_fast_mode_sql() {
    let bam_path = format!("{}/tests/data/ovl.bam", env!("CARGO_MANIFEST_DIR"));

    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    let sql = format!("SELECT * FROM depth('{bam_path}', true)");
    let rows = collect_coverage_sql(&ctx, &sql).await;

    let mt_rows: Vec<_> = rows.iter().filter(|r| r.0 == "MT").collect();

    assert_eq!(mt_rows.len(), 3);
    assert_eq!(mt_rows[0], &("MT".to_string(), 0, 5, 1));
    assert_eq!(mt_rows[1], &("MT".to_string(), 6, 41, 2));
    assert_eq!(mt_rows[2], &("MT".to_string(), 42, 79, 1));
}

/// Test coverage via SQL UDTF for overlapping pairs.
/// Uses `true` (0-based) to match mosdepth 0-based expectations.
#[tokio::test(flavor = "multi_thread")]
async fn test_overlapping_pairs_sql() {
    let bam_path = format!(
        "{}/tests/data/overlapping-pairs.bam",
        env!("CARGO_MANIFEST_DIR")
    );

    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    let sql = format!("SELECT * FROM depth('{bam_path}', true)");
    let rows = collect_coverage_sql(&ctx, &sql).await;

    let chr1_rows: Vec<_> = rows.iter().filter(|r| r.0 == "1").collect();

    assert_eq!(chr1_rows.len(), 1);
    assert_eq!(chr1_rows[0], &("1".to_string(), 565173, 565252, 2));
}

/// Test that the default coordinate system (1-based) works via SQL.
#[tokio::test(flavor = "multi_thread")]
async fn test_ovl_fast_mode_sql_default_one_based() {
    let bam_path = format!("{}/tests/data/ovl.bam", env!("CARGO_MANIFEST_DIR"));

    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    // No second argument → default 1-based coordinates
    let sql = format!("SELECT * FROM depth('{bam_path}')");
    let rows = collect_coverage_sql(&ctx, &sql).await;

    let mt_rows: Vec<_> = rows.iter().filter(|r| r.0 == "MT").collect();

    assert_eq!(mt_rows.len(), 3);
    // 1-based: positions are shifted by +1 compared to 0-based
    assert_eq!(mt_rows[0], &("MT".to_string(), 1, 6, 1));
    assert_eq!(mt_rows[1], &("MT".to_string(), 7, 42, 2));
    assert_eq!(mt_rows[2], &("MT".to_string(), 43, 80, 1));
}

/// Helper to collect per-base coverage results from SQL UDTF.
async fn collect_per_base_sql(ctx: &SessionContext, sql: &str) -> Vec<(String, i32, i16)> {
    let df = ctx.sql(sql).await.unwrap();
    let batches = df.collect().await.unwrap();

    let mut all_rows = Vec::new();
    for batch in batches {
        let contigs = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let positions = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let covs = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        for i in 0..batch.num_rows() {
            all_rows.push((
                contigs.value(i).to_string(),
                positions.value(i),
                covs.value(i),
            ));
        }
    }
    all_rows
}

/// Test per-base output via SQL UDTF: `SELECT * FROM depth('path', true, true)`
///
/// The ovl.bam file has reads on MT contig. With per_base=true, we should get
/// one row per genomic position for each contig that has reads, covering the
/// full contig length from the BAM header.
#[tokio::test(flavor = "multi_thread")]
async fn test_per_base_sql() {
    let bam_path = format!("{}/tests/data/ovl.bam", env!("CARGO_MANIFEST_DIR"));

    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    let sql = format!("SELECT * FROM depth('{bam_path}', true, true)");
    let rows = collect_per_base_sql(&ctx, &sql).await;

    // Should have 3 columns per row, not 4
    assert!(!rows.is_empty(), "Expected per-base rows");

    // Filter to MT contig
    let mt_rows: Vec<_> = rows.iter().filter(|r| r.0 == "MT").collect();

    // MT contig should have rows for every position (full contig length from header)
    assert!(
        mt_rows.len() > 80,
        "Per-base MT should cover full contig length, got {} rows",
        mt_rows.len()
    );

    // Verify 0-based positions start at 0
    assert_eq!(mt_rows[0].1, 0, "First MT position should be 0 (0-based)");

    // Verify positions are sequential
    for (i, row) in mt_rows.iter().enumerate() {
        assert_eq!(
            row.1, i as i32,
            "Position should be sequential at index {i}"
        );
    }

    // Verify coverage at known positions from block-based test:
    // Block: MT 0-5 cov=1, 6-41 cov=2, 42-79 cov=1
    assert_eq!(mt_rows[0].2, 1, "pos 0 should have coverage 1");
    assert_eq!(mt_rows[5].2, 1, "pos 5 should have coverage 1");
    assert_eq!(mt_rows[6].2, 2, "pos 6 should have coverage 2");
    assert_eq!(mt_rows[41].2, 2, "pos 41 should have coverage 2");
    assert_eq!(mt_rows[42].2, 1, "pos 42 should have coverage 1");
    assert_eq!(mt_rows[79].2, 1, "pos 79 should have coverage 1");

    // Positions beyond the reads should have coverage 0
    assert_eq!(mt_rows[80].2, 0, "pos 80 should have coverage 0");
}

/// Test per-base output schema has correct columns.
#[tokio::test(flavor = "multi_thread")]
async fn test_per_base_schema() {
    use datafusion_bio_function_pileup::schema::COORDINATE_SYSTEM_METADATA_KEY;

    let bam_path = format!("{}/tests/data/ovl.bam", env!("CARGO_MANIFEST_DIR"));

    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    let sql = format!("SELECT * FROM depth('{bam_path}', true, true)");
    let df = ctx.sql(&sql).await.unwrap();
    let schema = df.schema().inner().clone();

    assert_eq!(
        schema.fields().len(),
        3,
        "Per-base schema should have 3 columns"
    );
    assert_eq!(schema.field(0).name(), "contig");
    assert_eq!(schema.field(1).name(), "pos");
    assert_eq!(schema.field(2).name(), "coverage");
    assert_eq!(
        schema.metadata().get(COORDINATE_SYSTEM_METADATA_KEY),
        Some(&"true".to_string()),
    );
}

/// Test that coordinate metadata is present in the output schema.
#[tokio::test(flavor = "multi_thread")]
async fn test_coordinate_metadata_in_schema() {
    use datafusion_bio_function_pileup::schema::COORDINATE_SYSTEM_METADATA_KEY;

    let bam_path = format!("{}/tests/data/ovl.bam", env!("CARGO_MANIFEST_DIR"));

    let ctx = SessionContext::new();
    register_pileup_functions(&ctx);

    // 0-based
    let sql = format!("SELECT * FROM depth('{bam_path}', true)");
    let df = ctx.sql(&sql).await.unwrap();
    let schema = df.schema().inner().clone();
    assert_eq!(
        schema.metadata().get(COORDINATE_SYSTEM_METADATA_KEY),
        Some(&"true".to_string()),
    );

    // 1-based (default)
    let sql = format!("SELECT * FROM depth('{bam_path}')");
    let df = ctx.sql(&sql).await.unwrap();
    let schema = df.schema().inner().clone();
    assert_eq!(
        schema.metadata().get(COORDINATE_SYSTEM_METADATA_KEY),
        Some(&"false".to_string()),
    );
}
