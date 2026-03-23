//! E2E annotation benchmark: read VCF → annotate → write VCF output.
//!
//! Exercises the full pipeline including VCF reading, annotation with
//! transcript engine, and VCF writing via streaming sink with progress bar.
//!
//! Usage:
//!   cargo run --release --features kv-cache --example bench_annotate_vcf -- \
//!     --input <vcf_path> \
//!     --cache <cache_dir> \
//!     --output <output.vcf> \
//!     [--backend parquet|fjall] \
//!     [--everything] \
//!     [--extended-probes] \
//!     [--reference-fasta <path>] \
//!     [--compression none|gzip|bgzf] \
//!     [--limit <n>]

use std::sync::Arc;
use std::time::Instant;

use datafusion::common::Result;
use datafusion::datasource::TableProvider;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::VcfCompressionType;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::{register_vep_functions, vcf_sink};
use indicatif::{ProgressBar, ProgressStyle};

struct Args {
    input: String,
    cache: String,
    output: String,
    backend: String,
    everything: bool,
    extended_probes: bool,
    reference_fasta: Option<String>,
    compression: VcfCompressionType,
    limit: Option<usize>,
}

fn parse_args() -> Args {
    let args: Vec<String> = std::env::args().collect();
    let mut input = None;
    let mut cache = None;
    let mut output = None;
    let mut backend = "parquet".to_string();
    let mut everything = false;
    let mut extended_probes = false;
    let mut reference_fasta = None;
    let mut compression = VcfCompressionType::Plain;
    let mut limit = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--input" => {
                i += 1;
                input = Some(args[i].clone());
            }
            "--cache" => {
                i += 1;
                cache = Some(args[i].clone());
            }
            "--output" => {
                i += 1;
                output = Some(args[i].clone());
            }
            "--backend" => {
                i += 1;
                backend = args[i].clone();
            }
            "--everything" => everything = true,
            "--extended-probes" => extended_probes = true,
            "--reference-fasta" => {
                i += 1;
                reference_fasta = Some(args[i].clone());
            }
            "--compression" => {
                i += 1;
                compression = match args[i].as_str() {
                    "gzip" => VcfCompressionType::Gzip,
                    "bgzf" => VcfCompressionType::Bgzf,
                    _ => VcfCompressionType::Plain,
                };
            }
            "--limit" => {
                i += 1;
                limit = args[i].parse().ok();
            }
            other => {
                eprintln!("Unknown argument: {other}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let input = input.unwrap_or_else(|| {
        eprintln!(
            "Usage: {} --input <vcf> --cache <dir> --output <vcf> [options]",
            args[0]
        );
        std::process::exit(1);
    });
    let cache = cache.unwrap_or_else(|| {
        eprintln!("--cache is required");
        std::process::exit(1);
    });
    let output = output.unwrap_or_else(|| {
        eprintln!("--output is required");
        std::process::exit(1);
    });

    Args {
        input,
        cache,
        output,
        backend,
        everything,
        extended_probes,
        reference_fasta,
        compression,
        limit,
    }
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = parse_args();

    eprintln!("=== bench_annotate_vcf ===");
    eprintln!("  input:      {}", args.input);
    eprintln!("  cache:      {}", args.cache);
    eprintln!("  output:     {}", args.output);
    eprintln!("  backend:    {}", args.backend);
    eprintln!("  everything: {}", args.everything);
    eprintln!("  ext_probes: {}", args.extended_probes);
    eprintln!(
        "  ref_fasta:  {}",
        args.reference_fasta.as_deref().unwrap_or("(none)")
    );
    eprintln!(
        "  compress:   {}",
        match args.compression {
            VcfCompressionType::Gzip => "gzip",
            VcfCompressionType::Bgzf => "bgzf",
            _ => "none",
        }
    );
    if let Some(n) = args.limit {
        eprintln!("  limit:      {n}");
    }

    // ── Step 1: Read VCF ──
    let t0 = Instant::now();
    let input_path = args.input.clone();
    let vcf_provider = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(input_path, None, None, None, false)
    })
    .await
    .map_err(|e| datafusion::common::DataFusionError::External(Box::new(e)))??;
    let vcf_fields = vcf_provider.schema().fields().len();

    // Single partition — required for correct annotation pipeline.
    let config = datafusion::prelude::SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    register_vep_functions(&ctx);
    ctx.register_table("vcf", Arc::new(vcf_provider))?;

    // Count input rows for progress bar.
    let total_input = ctx
        .sql("SELECT COUNT(*) AS n FROM vcf")
        .await?
        .collect()
        .await?[0]
        .column(0)
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .map(|a| a.value(0) as u64)
        .unwrap_or(0);

    eprintln!(
        "\n[{:.1}s] VCF: {} columns, {} variants",
        t0.elapsed().as_secs_f64(),
        vcf_fields,
        total_input
    );

    // ── Step 2: Annotate + stream to VCF with progress ──
    // Show the progress bar immediately (at 0%) so the user sees activity
    // during the variation lookup phase before any batches arrive.
    let pb = ProgressBar::new(total_input);
    pb.set_style(
        ProgressStyle::with_template(
            "  {spinner:.green} {bar:40.cyan/blue} {pos}/{len} [{elapsed_precise}] (eta {eta})",
        )
        .unwrap()
        .progress_chars("##-"),
    );
    pb.set_message("annotating...");
    pb.enable_steady_tick(std::time::Duration::from_millis(200));
    let pb_cb = pb.clone();

    let annotate_config = vcf_sink::AnnotateVcfConfig {
        everything: args.everything,
        extended_probes: args.extended_probes,
        reference_fasta_path: args.reference_fasta.clone(),
        use_fjall: args.backend == "fjall",
        compression: args.compression,
        on_batch_written: Some(Box::new(move |n, _| {
            pb_cb.inc(n as u64);
        })),
        ..Default::default()
    };

    let t_annotate = Instant::now();
    let output_path = std::path::Path::new(&args.output);
    let rows = vcf_sink::annotate_to_vcf(
        &ctx,
        "vcf",
        &args.cache,
        "parquet",
        output_path,
        &annotate_config,
    )
    .await?;
    pb.finish_and_clear();
    let annotate_secs = t_annotate.elapsed().as_secs_f64();

    let output_size = std::fs::metadata(output_path).map(|m| m.len()).unwrap_or(0);

    eprintln!("\n=== Results ===");
    eprintln!("  rows:       {rows}");
    eprintln!("  time:       {annotate_secs:.2}s");
    eprintln!(
        "  throughput: {:.0} variants/s",
        rows as f64 / annotate_secs
    );
    eprintln!(
        "  output:     {:.1} MB",
        output_size as f64 / 1024.0 / 1024.0
    );
    eprintln!(
        "  total:      {:.2}s (including VCF read)",
        t0.elapsed().as_secs_f64()
    );

    Ok(())
}
