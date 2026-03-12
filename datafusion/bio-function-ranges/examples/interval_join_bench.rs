use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{ArrayRef, Int32Array, Int64Array, RecordBatch, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::util::pretty::pretty_format_batches;
use datafusion::config::ConfigOptions;
use datafusion::datasource::MemTable;
use datafusion::error::Result;
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_ranges::session_context::{Algorithm, BioConfig, BioSessionExt};

const CONTIGS: usize = 24;
const BUCKETS: usize = 64;
const ROWS_PER_GROUP: usize = 32;
const BATCH_ROWS: usize = 4096;
const PAYLOAD_COLS: usize = 4;
const SPARSE_SHIFT: i32 = 100;
const NO_OVERLAP_SHIFT: i32 = 20_000;

#[derive(Clone)]
struct TableData {
    schema: SchemaRef,
    batches: Vec<RecordBatch>,
}

#[derive(Clone, Copy)]
struct Case {
    name: &'static str,
    prefer_interval_join: bool,
    algorithm: Algorithm,
    low_memory: bool,
}

#[derive(Clone, Copy)]
struct QueryCase {
    name: &'static str,
    sql: &'static str,
    expected_plan_fragment: &'static str,
    parse_count: bool,
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let left = make_table(0, 0)?;
    let right_sparse = make_table(SPARSE_SHIFT, 10_000)?;
    let right_no_overlap = make_table(NO_OVERLAP_SHIFT, 20_000)?;

    let cases = [
        Case {
            name: "hash_join",
            prefer_interval_join: false,
            algorithm: Algorithm::Coitrees,
            low_memory: false,
        },
        Case {
            name: "interval_coitrees",
            prefer_interval_join: true,
            algorithm: Algorithm::Coitrees,
            low_memory: false,
        },
        Case {
            name: "interval_coitrees_low_memory",
            prefer_interval_join: true,
            algorithm: Algorithm::Coitrees,
            low_memory: true,
        },
        Case {
            name: "interval_lapper",
            prefer_interval_join: true,
            algorithm: Algorithm::Lapper,
            low_memory: false,
        },
        Case {
            name: "interval_superintervals",
            prefer_interval_join: true,
            algorithm: Algorithm::SuperIntervals,
            low_memory: false,
        },
    ];

    let queries = [
        QueryCase {
            name: "bucketed_no_overlap_count",
            sql: concat!(
                "SELECT COUNT(*) AS n FROM reads JOIN targets ON ",
                "reads.contig = targets.contig ",
                "AND reads.bucket = targets.bucket ",
                "AND reads.pos_start <= targets.pos_end ",
                "AND reads.pos_end >= targets.pos_start"
            ),
            expected_plan_fragment: "bucket@1, bucket@1",
            parse_count: true,
        },
        QueryCase {
            name: "bucketed_sparse_count",
            sql: concat!(
                "SELECT COUNT(*) AS n FROM reads JOIN targets ON ",
                "reads.contig = targets.contig ",
                "AND reads.bucket = targets.bucket ",
                "AND reads.pos_start <= targets.pos_end ",
                "AND reads.pos_end >= targets.pos_start"
            ),
            expected_plan_fragment: "bucket@1, bucket@1",
            parse_count: true,
        },
        QueryCase {
            name: "bucketed_sparse_wide",
            sql: concat!(
                "SELECT reads.*, targets.* FROM reads JOIN targets ON ",
                "reads.contig = targets.contig ",
                "AND reads.bucket = targets.bucket ",
                "AND reads.pos_start <= targets.pos_end ",
                "AND reads.pos_end >= targets.pos_start"
            ),
            expected_plan_fragment: "bucket@1, bucket@1",
            parse_count: false,
        },
        QueryCase {
            name: "contig_only_sparse_count",
            sql: concat!(
                "SELECT COUNT(*) AS n FROM reads JOIN targets ON ",
                "reads.contig = targets.contig ",
                "AND reads.pos_start <= targets.pos_end ",
                "AND reads.pos_end >= targets.pos_start"
            ),
            expected_plan_fragment: "contig@0, contig@0",
            parse_count: true,
        },
    ];

    println!(
        "dataset rows: left={}, right={}, batches_per_table={}",
        total_rows(),
        total_rows(),
        left.batches.len()
    );

    for query in &queries {
        let target = if query.name.contains("no_overlap") {
            &right_no_overlap
        } else {
            &right_sparse
        };

        println!("\n== {} ==", query.name);
        for case in &cases {
            let ctx = build_context(*case, &left, target).await?;
            let explain = explain_plan(&ctx, query.sql).await?;
            let expected_plan_kind = if case.prefer_interval_join {
                "IntervalJoinExec"
            } else {
                "HashJoinExec"
            };

            assert!(
                explain.contains(expected_plan_kind),
                "expected {expected_plan_kind} in plan:\n{explain}"
            );
            assert!(
                explain.contains(query.expected_plan_fragment),
                "expected plan fragment {:?} in plan:\n{}",
                query.expected_plan_fragment,
                explain
            );

            let (elapsed_ms, value, rows) = run_query(&ctx, query.sql, query.parse_count).await?;
            println!(
                "{:<28} elapsed_ms={:<8} value={:<10} output_rows={:<8} low_memory={}",
                case.name, elapsed_ms, value, rows, case.low_memory
            );
        }
    }

    println!("\n== explain_analyze bucketed_count ==");
    for case in &cases[..2] {
        let ctx = build_context(*case, &left, &right_sparse).await?;
        let analyzed = explain_analyze(
            &ctx,
            concat!(
                "SELECT COUNT(*) AS n FROM reads JOIN targets ON ",
                "reads.contig = targets.contig ",
                "AND reads.bucket = targets.bucket ",
                "AND reads.pos_start <= targets.pos_end ",
                "AND reads.pos_end >= targets.pos_start"
            ),
        )
        .await?;
        println!("\n-- {} --\n{}", case.name, analyzed);
    }

    Ok(())
}

fn total_rows() -> usize {
    CONTIGS * BUCKETS * ROWS_PER_GROUP
}

fn make_table(interval_shift: i32, payload_offset: i32) -> Result<TableData> {
    let schema = make_schema();
    let total_rows = total_rows();

    let mut contigs = Vec::with_capacity(total_rows);
    let mut buckets = Vec::with_capacity(total_rows);
    let mut starts = Vec::with_capacity(total_rows);
    let mut ends = Vec::with_capacity(total_rows);
    let mut payloads = (0..PAYLOAD_COLS)
        .map(|_| Vec::with_capacity(total_rows))
        .collect::<Vec<_>>();

    let mut row_id = 0i32;
    for contig_idx in 0..CONTIGS {
        let contig = format!("chr{}", contig_idx + 1);
        for bucket_idx in 0..BUCKETS {
            let bucket_base = (bucket_idx as i32) * 1_000_000;
            for slot in 0..ROWS_PER_GROUP {
                let start = bucket_base + (slot as i32) * 200 + interval_shift;
                let end = start + 150;

                contigs.push(contig.clone());
                buckets.push(bucket_idx as i32);
                starts.push(start);
                ends.push(end);

                for (payload_idx, col) in payloads.iter_mut().enumerate() {
                    col.push(payload_offset + row_id * (payload_idx as i32 + 1));
                }
                row_id += 1;
            }
        }
    }

    let arrays = build_arrays(contigs, buckets, starts, ends, payloads);
    let mut batches = Vec::new();
    for offset in (0..total_rows).step_by(BATCH_ROWS) {
        let len = usize::min(BATCH_ROWS, total_rows - offset);
        let sliced = arrays
            .iter()
            .map(|array| array.slice(offset, len))
            .collect::<Vec<_>>();
        batches.push(RecordBatch::try_new(schema.clone(), sliced)?);
    }

    Ok(TableData { schema, batches })
}

fn make_schema() -> SchemaRef {
    let mut fields = vec![
        Field::new("contig", DataType::Utf8, false),
        Field::new("bucket", DataType::Int32, false),
        Field::new("pos_start", DataType::Int32, false),
        Field::new("pos_end", DataType::Int32, false),
    ];
    for i in 0..PAYLOAD_COLS {
        fields.push(Field::new(format!("payload{i}"), DataType::Int32, false));
    }
    Arc::new(Schema::new(fields))
}

fn build_arrays(
    contigs: Vec<String>,
    buckets: Vec<i32>,
    starts: Vec<i32>,
    ends: Vec<i32>,
    payloads: Vec<Vec<i32>>,
) -> Vec<ArrayRef> {
    let mut arrays: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(contigs)),
        Arc::new(Int32Array::from(buckets)),
        Arc::new(Int32Array::from(starts)),
        Arc::new(Int32Array::from(ends)),
    ];
    for payload in payloads {
        arrays.push(Arc::new(Int32Array::from(payload)));
    }
    arrays
}

async fn build_context(case: Case, left: &TableData, right: &TableData) -> Result<SessionContext> {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false)
        .with_target_partitions(1);
    let ctx = SessionContext::new_with_bio(config);

    ctx.sql(
        format!(
            "SET bio.prefer_interval_join = {}",
            case.prefer_interval_join
        )
        .as_str(),
    )
    .await?;
    ctx.sql(format!("SET bio.interval_join_algorithm = {}", case.algorithm).as_str())
        .await?;
    ctx.sql(format!("SET bio.interval_join_low_memory = {}", case.low_memory).as_str())
        .await?;

    let left_provider = MemTable::try_new(left.schema.clone(), vec![left.batches.clone()])?;
    let right_provider = MemTable::try_new(right.schema.clone(), vec![right.batches.clone()])?;

    ctx.register_table("reads", Arc::new(left_provider))?;
    ctx.register_table("targets", Arc::new(right_provider))?;
    Ok(ctx)
}

async fn explain_plan(ctx: &SessionContext, sql: &str) -> Result<String> {
    let batches = ctx.sql(&format!("EXPLAIN {sql}")).await?.collect().await?;
    Ok(pretty_format_batches(&batches)?.to_string())
}

async fn explain_analyze(ctx: &SessionContext, sql: &str) -> Result<String> {
    let batches = ctx
        .sql(&format!("EXPLAIN ANALYZE {sql}"))
        .await?
        .collect()
        .await?;
    Ok(pretty_format_batches(&batches)?.to_string())
}

async fn run_query(
    ctx: &SessionContext,
    sql: &str,
    parse_count: bool,
) -> Result<(u128, i64, usize)> {
    let start = Instant::now();
    let batches = ctx.sql(sql).await?.collect().await?;
    let elapsed_ms = start.elapsed().as_millis();
    let rows = batches.iter().map(RecordBatch::num_rows).sum::<usize>();

    let value = if parse_count {
        let batch = &batches[0];
        let array = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int64Array>()
            .expect("count(*) should produce Int64Array");
        array.value(0)
    } else {
        rows as i64
    };

    Ok((elapsed_ms, value, rows))
}
