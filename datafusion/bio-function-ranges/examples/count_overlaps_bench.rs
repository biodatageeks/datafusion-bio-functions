use std::error::Error;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;

use datafusion::arrow::array::{Array, Int64Array, UInt64Array};
use datafusion::config::ConfigOptions;
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_ranges::{
    BioConfig, BioSessionExt, CountOverlapsBackendMode, register_ranges_functions,
};

#[derive(Debug)]
struct Args {
    left: PathBuf,
    right: PathBuf,
    backend: CountOverlapsBackendMode,
    threads: usize,
    repeats: usize,
    warmups: usize,
    batch_size: Option<usize>,
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse()?;
    let mut elapsed = Vec::with_capacity(args.repeats);

    for warmup in 0..args.warmups {
        let run = run_once(&args).await?;
        eprintln!(
            "warmup={},backend={},threads={},elapsed_ms={:.3}",
            warmup, args.backend, args.threads, run.elapsed_ms
        );
    }

    for repeat in 0..args.repeats {
        let run = run_once(&args).await?;
        elapsed.push(run.elapsed_ms);
        println!(
            "repeat={},backend={},threads={},elapsed_ms={:.3}",
            repeat, args.backend, args.threads, run.elapsed_ms
        );
        println!("{}", run.result);
    }

    if !elapsed.is_empty() {
        let mean = elapsed.iter().sum::<f64>() / elapsed.len() as f64;
        println!(
            "summary,backend={},threads={},repeats={},mean_elapsed_ms={:.3}",
            args.backend, args.threads, args.repeats, mean
        );
    }

    Ok(())
}

#[derive(Debug)]
struct RunResult {
    elapsed_ms: f64,
    result: String,
}

async fn run_once(args: &Args) -> Result<RunResult, Box<dyn Error>> {
    let mut bio_config = BioConfig::default();
    bio_config.count_overlaps_backend = args.backend;

    let mut config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(bio_config)
        .with_information_schema(true)
        .with_repartition_joins(false)
        .with_target_partitions(args.threads);
    if let Some(batch_size) = args.batch_size {
        config = config.with_batch_size(batch_size);
    }

    let ctx = SessionContext::new_with_bio(config);
    register_ranges_functions(&ctx);
    register_parquet(&ctx, "reads", &args.left).await?;
    register_parquet(&ctx, "targets", &args.right).await?;

    let start = Instant::now();
    let batches = ctx
        .sql(
            r#"SELECT COUNT(*) AS n_rows, SUM("count") AS total_overlaps
               FROM count_overlaps('reads', 'targets')"#,
        )
        .await?
        .collect()
        .await?;
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
    let batch = batches
        .first()
        .ok_or_else(|| "benchmark query returned no batches".to_string())?;
    let n_rows = batch
        .column_by_name("n_rows")
        .ok_or_else(|| "n_rows column not found in benchmark output".to_string())
        .and_then(|array| scalar_i64(array.as_ref(), "n_rows"))?;
    let total_overlaps = batch
        .column_by_name("total_overlaps")
        .ok_or_else(|| "total_overlaps column not found in benchmark output".to_string())
        .and_then(|array| scalar_i64(array.as_ref(), "total_overlaps"))?;
    let result = format!("n_rows={n_rows},total_overlaps={total_overlaps}");

    Ok(RunResult { elapsed_ms, result })
}

fn scalar_i64(array: &dyn Array, name: &str) -> Result<i64, String> {
    if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        return Ok(array.value(0));
    }
    if let Some(array) = array.as_any().downcast_ref::<UInt64Array>() {
        return i64::try_from(array.value(0)).map_err(|_| format!("{name} value overflows i64"));
    }
    Err(format!(
        "{name} column has unsupported type {}",
        array.data_type()
    ))
}

async fn register_parquet(
    ctx: &SessionContext,
    table: &str,
    path: &Path,
) -> Result<(), Box<dyn Error>> {
    let path = path
        .to_str()
        .ok_or_else(|| format!("path for table {table} is not valid UTF-8"))?
        .replace('\'', "''");
    ctx.sql(&format!(
        "CREATE EXTERNAL TABLE {table} STORED AS PARQUET LOCATION '{path}'"
    ))
    .await?;
    Ok(())
}

impl Args {
    fn parse() -> Result<Self, Box<dyn Error>> {
        let mut left = PathBuf::from("/tmp/polars-bio-bench/databio/ex-rna");
        let mut right = PathBuf::from("/tmp/polars-bio-bench/databio/ex-anno");
        let mut backend = CountOverlapsBackendMode::Auto;
        let mut threads = 1usize;
        let mut repeats = 3usize;
        let mut warmups = 1usize;
        let mut batch_size = None;

        let mut args = std::env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--left" => left = PathBuf::from(next_value(&mut args, "--left")?),
                "--right" => right = PathBuf::from(next_value(&mut args, "--right")?),
                "--backend" => {
                    backend =
                        CountOverlapsBackendMode::from_str(&next_value(&mut args, "--backend")?)?
                }
                "--threads" => threads = next_value(&mut args, "--threads")?.parse()?,
                "--repeats" => repeats = next_value(&mut args, "--repeats")?.parse()?,
                "--warmups" => warmups = next_value(&mut args, "--warmups")?.parse()?,
                "--batch-size" => {
                    batch_size = Some(next_value(&mut args, "--batch-size")?.parse()?)
                }
                "--help" | "-h" => {
                    print_usage();
                    std::process::exit(0);
                }
                other => return Err(format!("unknown argument: {other}").into()),
            }
        }

        if threads == 0 {
            return Err("--threads must be greater than zero".into());
        }

        Ok(Self {
            left,
            right,
            backend,
            threads,
            repeats,
            warmups,
            batch_size,
        })
    }
}

fn next_value(
    args: &mut impl Iterator<Item = String>,
    flag: &'static str,
) -> Result<String, Box<dyn Error>> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}").into())
}

fn print_usage() {
    println!(
        "Usage: count_overlaps_bench [--left PATH] [--right PATH] [--backend auto|cpu|apple_gpu] [--threads N] [--repeats N] [--warmups N] [--batch-size N]"
    );
}
