use std::env;
use std::time::Instant;

use datafusion::arrow::array::{Int64Array, UInt32Array, UInt64Array};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::arrow::util::pretty::pretty_format_batches;
use datafusion::common::{DataFusionError, Result};
use datafusion::config::ConfigOptions;
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_ranges::{
    Algorithm, BioConfig, BioSessionExt, register_ranges_functions,
};
use futures::StreamExt;

fn count_from_batches(batches: &[RecordBatch]) -> Result<u64> {
    let batch = batches.first().ok_or_else(|| {
        DataFusionError::Execution("COUNT query returned no record batches".to_string())
    })?;
    let column = batch.column(0);

    if let Some(values) = column.as_any().downcast_ref::<UInt64Array>() {
        return Ok(values.value(0));
    }
    if let Some(values) = column.as_any().downcast_ref::<Int64Array>() {
        return Ok(values.value(0) as u64);
    }
    if let Some(values) = column.as_any().downcast_ref::<UInt32Array>() {
        return Ok(values.value(0) as u64);
    }

    Err(DataFusionError::Execution(format!(
        "unsupported COUNT type: {:?}",
        column.data_type()
    )))
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let args = env::args().collect::<Vec<_>>();
    let left = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("/tmp/polars-bio-bench/databio/ex-rna");
    let right = args
        .get(2)
        .map(String::as_str)
        .unwrap_or("/tmp/polars-bio-bench/databio/ex-anno");
    let runs = args
        .get(3)
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(3);
    let target_partitions = args
        .get(4)
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(1);
    let batch_size = args
        .get(5)
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(8192);
    let mode = args.get(6).map(String::as_str).unwrap_or("aggregate");

    for run in 0..runs {
        let mut options = ConfigOptions::new();
        options.set("datafusion.execution.coalesce_batches", "false")?;
        options.set("datafusion.optimizer.repartition_joins", "false")?;

        let mut bio_config = BioConfig::default();
        bio_config.prefer_interval_join = true;
        bio_config.interval_join_algorithm = Algorithm::Coitrees;
        bio_config.interval_join_low_memory = false;
        let config = SessionConfig::from(options)
            .with_option_extension(bio_config)
            .with_information_schema(true)
            .with_repartition_joins(false)
            .with_target_partitions(target_partitions)
            .with_batch_size(batch_size);
        let ctx = SessionContext::new_with_bio(config);
        register_ranges_functions(&ctx);

        ctx.sql(&format!(
            "CREATE EXTERNAL TABLE left_ranges STORED AS PARQUET LOCATION '{left}'"
        ))
        .await?;
        ctx.sql(&format!(
            "CREATE EXTERNAL TABLE right_ranges STORED AS PARQUET LOCATION '{right}'"
        ))
        .await?;

        if mode == "explain" {
            let batches = ctx
                .sql(
                    "EXPLAIN ANALYZE
                     SELECT COUNT(*) AS count
                     FROM overlap(
                       'left_ranges',
                       'right_ranges',
                       'contig',
                       'pos_start',
                       'pos_end',
                       'strict'
                     )",
                )
                .await?
                .collect()
                .await?;
            println!("{}", pretty_format_batches(&batches)?);
            continue;
        }

        let started = Instant::now();
        let count = if mode == "stream" {
            let mut stream = ctx
                .sql(
                    "SELECT *
                     FROM overlap(
                       'left_ranges',
                       'right_ranges',
                       'contig',
                       'pos_start',
                       'pos_end',
                       'strict'
                     )",
                )
                .await?
                .execute_stream()
                .await?;
            let mut count = 0_u64;
            while let Some(batch) = stream.next().await {
                count += batch?.num_rows() as u64;
            }
            count
        } else {
            let batches = ctx
                .sql(
                    "SELECT COUNT(*) AS count
                     FROM overlap(
                       'left_ranges',
                       'right_ranges',
                       'contig',
                       'pos_start',
                       'pos_end',
                       'strict'
                     )",
                )
                .await?
                .collect()
                .await?;
            count_from_batches(&batches)?
        };
        let elapsed = started.elapsed();

        println!(
            "datafusion {mode} run {}: seconds={:.3} count={count}",
            run + 1,
            elapsed.as_secs_f64()
        );
    }

    Ok(())
}
