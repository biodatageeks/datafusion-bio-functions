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
    timings: bool,
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
        if args.timings {
            eprintln!("{}", run.timing_line("warmup", warmup, &args));
            if let Some(line) = run.count_overlaps_timing_line("warmup", warmup, &args) {
                eprintln!("{line}");
            }
            if let Some(line) = run.gpu_timing_line("warmup", warmup, &args) {
                eprintln!("{line}");
            }
        }
    }

    for repeat in 0..args.repeats {
        let run = run_once(&args).await?;
        elapsed.push(run.elapsed_ms);
        println!(
            "repeat={},backend={},threads={},elapsed_ms={:.3}",
            repeat, args.backend, args.threads, run.elapsed_ms
        );
        if args.timings {
            println!("{}", run.timing_line("repeat", repeat, &args));
            if let Some(line) = run.count_overlaps_timing_line("repeat", repeat, &args) {
                println!("{line}");
            }
            if let Some(line) = run.gpu_timing_line("repeat", repeat, &args) {
                println!("{line}");
            }
        }
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
    setup_ms: f64,
    sql_ms: f64,
    collect_ms: f64,
    result_ms: f64,
    count_overlaps_timings: Option<CountOverlapsTimingResult>,
    gpu_timings: Option<GpuTimingResult>,
    result: String,
}

#[derive(Debug, Clone, Copy)]
struct CountOverlapsTimingResult {
    scans: u64,
    left_rows: u64,
    left_collect_ms: f64,
    backend_build_ms: f64,
    right_plan_ms: f64,
    scan_total_ms: f64,
}

#[derive(Debug, Clone, Copy)]
struct GpuTimingResult {
    batches: u64,
    rows: u64,
    encode_ms: f64,
    buffer_ms: f64,
    command_ms: f64,
    wait_ms: f64,
    readback_ms: f64,
    output_ms: f64,
    total_ms: f64,
}

async fn run_once(args: &Args) -> Result<RunResult, Box<dyn Error>> {
    let setup_start = Instant::now();
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
    let setup_ms = setup_start.elapsed().as_secs_f64() * 1000.0;

    let start = Instant::now();
    let dataframe = ctx
        .sql(
            r#"SELECT COUNT(*) AS n_rows, SUM("count") AS total_overlaps
               FROM count_overlaps('reads', 'targets')"#,
        )
        .await?;
    let sql_ms = start.elapsed().as_secs_f64() * 1000.0;

    let collect_start = Instant::now();
    datafusion_bio_function_ranges::reset_count_overlaps_timings();
    reset_gpu_timing_stats();
    let batches = dataframe.collect().await?;
    let count_overlaps_timings = snapshot_count_overlaps_timing_stats();
    let gpu_timings = snapshot_gpu_timing_stats();
    let collect_ms = collect_start.elapsed().as_secs_f64() * 1000.0;
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    let result_start = Instant::now();
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
    let result_ms = result_start.elapsed().as_secs_f64() * 1000.0;

    Ok(RunResult {
        elapsed_ms,
        setup_ms,
        sql_ms,
        collect_ms,
        result_ms,
        count_overlaps_timings,
        gpu_timings,
        result,
    })
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
        let mut timings = false;

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
                "--timings" => timings = true,
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
            timings,
        })
    }
}

impl RunResult {
    fn timing_line(&self, phase: &str, iteration: usize, args: &Args) -> String {
        format!(
            "timings,{phase}={iteration},backend={},threads={},setup_ms={:.3},sql_ms={:.3},collect_ms={:.3},result_ms={:.3},elapsed_ms={:.3}",
            args.backend,
            args.threads,
            self.setup_ms,
            self.sql_ms,
            self.collect_ms,
            self.result_ms,
            self.elapsed_ms
        )
    }

    fn gpu_timing_line(&self, phase: &str, iteration: usize, args: &Args) -> Option<String> {
        let timings = self.gpu_timings?;
        Some(format!(
            "gpu_timings,{phase}={iteration},backend={},threads={},batches={},rows={},encode_ms={:.3},buffer_ms={:.3},command_ms={:.3},wait_ms={:.3},readback_ms={:.3},output_ms={:.3},total_ms={:.3}",
            args.backend,
            args.threads,
            timings.batches,
            timings.rows,
            timings.encode_ms,
            timings.buffer_ms,
            timings.command_ms,
            timings.wait_ms,
            timings.readback_ms,
            timings.output_ms,
            timings.total_ms
        ))
    }

    fn count_overlaps_timing_line(
        &self,
        phase: &str,
        iteration: usize,
        args: &Args,
    ) -> Option<String> {
        let timings = self.count_overlaps_timings?;
        Some(format!(
            "count_overlaps_timings,{phase}={iteration},backend={},threads={},scans={},left_rows={},left_collect_ms={:.3},backend_build_ms={:.3},right_plan_ms={:.3},scan_total_ms={:.3}",
            args.backend,
            args.threads,
            timings.scans,
            timings.left_rows,
            timings.left_collect_ms,
            timings.backend_build_ms,
            timings.right_plan_ms,
            timings.scan_total_ms
        ))
    }
}

fn snapshot_count_overlaps_timing_stats() -> Option<CountOverlapsTimingResult> {
    let snapshot = datafusion_bio_function_ranges::snapshot_count_overlaps_timings();
    (snapshot.scans > 0).then(|| CountOverlapsTimingResult {
        scans: snapshot.scans,
        left_rows: snapshot.left_rows,
        left_collect_ms: ns_to_ms(snapshot.left_collect_ns),
        backend_build_ms: ns_to_ms(snapshot.backend_build_ns),
        right_plan_ms: ns_to_ms(snapshot.right_plan_ns),
        scan_total_ms: ns_to_ms(snapshot.scan_total_ns),
    })
}

#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
fn reset_gpu_timing_stats() {
    datafusion_bio_function_ranges::reset_apple_gpu_timings();
}

#[cfg(not(all(feature = "apple-gpu", target_os = "macos")))]
fn reset_gpu_timing_stats() {}

#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
fn snapshot_gpu_timing_stats() -> Option<GpuTimingResult> {
    let snapshot = datafusion_bio_function_ranges::snapshot_apple_gpu_timings();
    (snapshot.batches > 0).then(|| GpuTimingResult {
        batches: snapshot.batches,
        rows: snapshot.rows,
        encode_ms: ns_to_ms(snapshot.encode_ns),
        buffer_ms: ns_to_ms(snapshot.buffer_ns),
        command_ms: ns_to_ms(snapshot.command_ns),
        wait_ms: ns_to_ms(snapshot.wait_ns),
        readback_ms: ns_to_ms(snapshot.readback_ns),
        output_ms: ns_to_ms(snapshot.output_ns),
        total_ms: ns_to_ms(snapshot.total_ns),
    })
}

#[cfg(not(all(feature = "apple-gpu", target_os = "macos")))]
fn snapshot_gpu_timing_stats() -> Option<GpuTimingResult> {
    None
}

fn ns_to_ms(ns: u64) -> f64 {
    ns as f64 / 1_000_000.0
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
        "Usage: count_overlaps_bench [--left PATH] [--right PATH] [--backend auto|cpu|apple_gpu] [--threads N] [--repeats N] [--warmups N] [--batch-size N] [--timings]"
    );
}
