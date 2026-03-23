//! E2E annotation benchmark: read VCF → annotate → write VCF output.
//!
//! Exercises the full pipeline including VCF reading, annotation with
//! transcript engine, and VCF writing via vcf_sink.
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
//!
//! Examples:
//!   # Parquet backend, default options
//!   cargo run --release --example bench_annotate_vcf -- \
//!     --input input.vcf --cache /data/cache --output out.vcf
//!
//!   # Fjall backend, --everything, gzip output
//!   cargo run --release --features kv-cache --example bench_annotate_vcf -- \
//!     --input input.vcf --cache /data/fjall_cache --output out.vcf.gz \
//!     --backend fjall --everything --extended-probes \
//!     --reference-fasta /ref/GRCh38.fa --compression gzip

use std::sync::Arc;
use std::time::Instant;

use datafusion::common::Result;
use datafusion::datasource::TableProvider;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::VcfCompressionType;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::{register_vep_functions, vcf_sink};

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

    println!("=== bench_annotate_vcf ===");
    println!("  input:      {}", args.input);
    println!("  cache:      {}", args.cache);
    println!("  output:     {}", args.output);
    println!("  backend:    {}", args.backend);
    println!("  everything: {}", args.everything);
    println!("  ext_probes: {}", args.extended_probes);
    println!(
        "  ref_fasta:  {}",
        args.reference_fasta.as_deref().unwrap_or("(none)")
    );
    println!(
        "  compress:   {}",
        match args.compression {
            VcfCompressionType::Gzip => "gzip",
            VcfCompressionType::Bgzf => "bgzf",
            _ => "none",
        }
    );
    if let Some(n) = args.limit {
        println!("  limit:      {n}");
    }

    // ── Read VCF ──
    // VcfTableProvider::new() uses futures::executor::block_on() internally,
    // which can deadlock inside a tokio runtime — use spawn_blocking.
    let t0 = Instant::now();
    let input_path = args.input.clone();
    let vcf_provider = tokio::task::spawn_blocking(move || {
        // None = include ALL INFO/FORMAT fields for VCF pass-through.
        VcfTableProvider::new(input_path, None, None, None, false)
    })
    .await
    .map_err(|e| datafusion::common::DataFusionError::External(Box::new(e)))??;
    let vcf_schema = vcf_provider.schema();
    let vcf_fields = vcf_schema.fields().len();
    println!(
        "\n[{:.1}s] VCF loaded: {} columns",
        t0.elapsed().as_secs_f64(),
        vcf_fields
    );

    // ── Set up context ──
    // Single partition: VcfTableProvider doesn't support multi-partition reads.
    let config = datafusion::prelude::SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    register_vep_functions(&ctx);
    ctx.register_table("vcf", Arc::new(vcf_provider))?;

    // ── Build annotation config ──
    let config = vcf_sink::AnnotateVcfConfig {
        everything: args.everything,
        extended_probes: args.extended_probes,
        reference_fasta_path: args.reference_fasta.clone(),
        use_fjall: args.backend == "fjall",
        compression: args.compression,
        ..Default::default()
    };

    // ── Annotate and write VCF ──
    let t_annotate = Instant::now();
    let output_path = std::path::Path::new(&args.output);
    let rows = vcf_sink::annotate_to_vcf(&ctx, "vcf", &args.cache, "parquet", output_path, &config)
        .await?;
    let annotate_secs = t_annotate.elapsed().as_secs_f64();

    let output_size = std::fs::metadata(output_path).map(|m| m.len()).unwrap_or(0);

    println!("[{:.1}s] Annotation + VCF write complete", annotate_secs);

    println!("\n=== Results ===");
    println!("  rows:       {rows}");
    println!("  time:       {annotate_secs:.2}s");
    println!(
        "  throughput: {:.0} variants/s",
        rows as f64 / annotate_secs
    );
    println!(
        "  output:     {:.1} MB",
        output_size as f64 / 1024.0 / 1024.0
    );
    println!(
        "  total:      {:.2}s (including VCF read)",
        t0.elapsed().as_secs_f64()
    );

    Ok(())
}
