//! Annotate a VCF end-to-end (optionally with a custom plugin cache) — the V1
//! golden-parity vehicle.
//!
//! ```text
//! cargo run --release -p datafusion-bio-function-vep --example annotate_vcf -- \
//!   --input <in.vcf> --cache <parquet cache root> --out <out.vcf> \
//!   --fasta <ref.fa> [--plugin-cache <plugin cache root>] [--everything] [--hgvs] \\
//!   [--preserve-record-layout]
//! ```

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::vcf_sink::{AnnotateVcfConfig, annotate_to_vcf};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
fn flag(args: &[String], key: &str) -> bool {
    args.iter().any(|a| a == key)
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "--input")
        .ok_or_else(|| DataFusionError::Execution("--input required".into()))?;
    let cache = arg(&args, "--cache")
        .ok_or_else(|| DataFusionError::Execution("--cache required".into()))?;
    let out =
        arg(&args, "--out").ok_or_else(|| DataFusionError::Execution("--out required".into()))?;
    let fasta = arg(&args, "--fasta");
    let plugin_cache = arg(&args, "--plugin-cache");

    // `AnnotateVcfConfig` is `#[non_exhaustive]`, so it is built by assignment
    // rather than a struct literal: a field this crate adds then defaults here
    // instead of breaking every caller's build.
    let mut config = AnnotateVcfConfig::default();
    config.everything = flag(&args, "--everything");
    config.hgvs = flag(&args, "--hgvs");
    config.reference_fasta_path = fasta;
    config.workers = arg(&args, "--workers")
        .and_then(|w| w.parse().ok())
        .unwrap_or(1);
    config.target_partitions = 1;
    config.plugin_cache_root = plugin_cache.map(std::path::PathBuf::from);
    config.preserve_record_layout = flag(&args, "--preserve-record-layout");

    eprintln!(
        "annotate: input={input} cache={cache} plugin_cache={:?} everything={} hgvs={}",
        config.plugin_cache_root, config.everything, config.hgvs
    );
    let n = annotate_to_vcf(&input, &cache, "parquet", &out, &config).await?;
    eprintln!("wrote {n} records -> {out}");
    Ok(())
}
