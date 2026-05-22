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
//!     [--target-partitions <n>] \
//!     [--forks <n>] \
//!     [--buffer-size <n>] \
//!     [--compression none|gzip|bgzf] \
//!     [--limit <n>] \
//!     [--no-progress]

use std::time::Instant;

use datafusion::common::Result;
use datafusion_bio_format_vcf::VcfCompressionType;
use datafusion_bio_function_vep::golden_benchmark::sample_gz_vcf_first_n;
use datafusion_bio_function_vep::vcf_sink;

struct Args {
    input: String,
    cache: String,
    output: String,
    backend: String,
    everything: bool,
    extended_probes: bool,
    reference_fasta: Option<String>,
    target_partitions: usize,
    forks: Option<usize>,
    buffer_size: usize,
    compression: VcfCompressionType,
    limit: Option<usize>,
    show_progress: bool,
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
    let mut target_partitions = 1;
    let mut forks = None;
    let mut buffer_size = vcf_sink::VEP_DEFAULT_BUFFER_SIZE;
    let mut compression = VcfCompressionType::Plain;
    let mut limit = None;
    let mut show_progress = true;

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
            "--target-partitions" => {
                i += 1;
                target_partitions = args[i].parse().unwrap_or(1);
            }
            "--forks" => {
                i += 1;
                forks = args[i].parse().ok();
            }
            "--buffer-size" => {
                i += 1;
                buffer_size = args[i].parse().unwrap_or(vcf_sink::VEP_DEFAULT_BUFFER_SIZE);
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
            "--no-progress" => show_progress = false,
            "--progress" => show_progress = true,
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
        target_partitions,
        forks,
        buffer_size,
        compression,
        limit,
        show_progress,
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
    eprintln!("  target_partitions: {}", args.target_partitions);
    eprintln!(
        "  forks:      {}",
        args.forks
            .map(|value| value.to_string())
            .unwrap_or_else(|| "(legacy target_partitions)".to_string())
    );
    eprintln!("  buffer_size: {}", args.buffer_size);
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
    eprintln!("  progress:   {}", args.show_progress);

    let mut input_for_run = args.input.clone();
    let _sample_dir = if let Some(limit) = args.limit {
        let sample_dir = tempfile::tempdir().map_err(|e| {
            datafusion::common::DataFusionError::Execution(format!(
                "failed to create sample tempdir: {e}"
            ))
        })?;
        let sampled_vcf = sample_dir.path().join("sampled.vcf");
        let sampled =
            sample_gz_vcf_first_n(std::path::Path::new(&args.input), &sampled_vcf, limit)?;
        eprintln!(
            "  sampled:    {sampled} variants -> {}",
            sampled_vcf.display()
        );
        input_for_run = sampled_vcf.display().to_string();
        Some(sample_dir)
    } else {
        None
    };

    // ── Annotate + stream to VCF ──
    let t0 = Instant::now();
    let annotate_config = vcf_sink::AnnotateVcfConfig {
        everything: args.everything,
        extended_probes: args.extended_probes,
        reference_fasta_path: args.reference_fasta.clone(),
        use_fjall: args.backend == "fjall",
        target_partitions: args.target_partitions,
        forks: args.forks,
        buffer_size: args.buffer_size,
        compression: args.compression,
        show_progress: args.show_progress,
        ..Default::default()
    };

    let t_annotate = Instant::now();
    let rows = vcf_sink::annotate_to_vcf(
        &input_for_run,
        &args.cache,
        "parquet",
        &args.output,
        &annotate_config,
    )
    .await?;
    let annotate_secs = t_annotate.elapsed().as_secs_f64();

    let output_size = std::fs::metadata(&args.output)
        .map(|m| m.len())
        .unwrap_or(0);

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
