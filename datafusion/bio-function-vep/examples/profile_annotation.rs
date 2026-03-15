//! Profile the VEP annotation pipeline with fine-grained timing.
//!
//! Usage:
//!   VEP_PROFILE=1 cargo run --release --example profile_annotation -- \
//!     /path/to/input.vcf.gz \
//!     /path/to/cache_dir \
//!     [sample_limit]
//!
//! The cache_dir should contain:
//!   115_GRCh38_variation_1_vep.parquet
//!   115_GRCh38_transcript_1_vep.parquet
//!   115_GRCh38_exon_1_vep.parquet
//!   115_GRCh38_translation_1_vep.parquet
//!   115_GRCh38_regulatory_1_vep.parquet
//!   115_GRCh38_motif_1_vep.parquet
//!
//! Example (1000 variants from chr1, --everything mode):
//!   VEP_PROFILE=1 cargo run --release --example profile_annotation -- \
//!     vep-benchmark/data/HG002_chr1.vcf.gz \
//!     /Users/mwiewior/research/data/vep/chr1-vep \
//!     1000 \
//!     --everything \
//!     --reference-fasta-path=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{ParquetReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::golden_benchmark::sample_gz_vcf_first_n;
use datafusion_bio_function_vep::register_vep_functions;

fn sql_literal(s: &str) -> String {
    s.replace('\'', "''")
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    // Ensure profiling is on.
    if std::env::var("VEP_PROFILE").is_err() {
        // SAFETY: called before any threads are spawned.
        unsafe { std::env::set_var("VEP_PROFILE", "1") };
    }

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: VEP_PROFILE=1 cargo run --release --example profile_annotation -- \
             <vcf.gz> <cache_dir> [sample_limit] [--everything] \
             [--reference-fasta-path=<path>] [--hgvs]"
        );
        std::process::exit(1);
    }

    let vcf_gz = PathBuf::from(&args[1]);
    let cache_dir = PathBuf::from(&args[2]);
    let sample_limit: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1000);

    let everything = args.iter().any(|a| a == "--everything");
    let hgvs = args.iter().any(|a| a == "--hgvs") || everything;
    let reference_fasta_path = args
        .iter()
        .find_map(|a| a.strip_prefix("--reference-fasta-path="))
        .map(|s| s.to_string());

    eprintln!("=== VEP Annotation Profiler ===");
    eprintln!("VCF:          {}", vcf_gz.display());
    eprintln!("Cache dir:    {}", cache_dir.display());
    eprintln!("Sample limit: {}", sample_limit);
    eprintln!("Everything:   {}", everything);
    eprintln!("HGVS:         {}", hgvs);
    eprintln!(
        "Ref FASTA:    {}",
        reference_fasta_path.as_deref().unwrap_or("none")
    );
    eprintln!();

    // Sample VCF.
    let work_dir = tempfile::tempdir()
        .map_err(|e| DataFusionError::Execution(format!("failed to create temp dir: {e}")))?;
    let sampled_vcf = work_dir.path().join("sampled.vcf");
    let t_sample = Instant::now();
    let n_sampled = sample_gz_vcf_first_n(&vcf_gz, &sampled_vcf, sample_limit)?;
    eprintln!(
        "[PREP] Sampled {} variants in {:.1}ms",
        n_sampled,
        t_sample.elapsed().as_secs_f64() * 1000.0
    );

    // Detect cache file naming (e.g. "115_GRCh38_variation_1_vep").
    let variation_file = find_parquet(&cache_dir, "variation")?;
    let base = extract_base(&variation_file);
    eprintln!("[PREP] Base name: {}", base);

    // Set up session.
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    register_vep_functions(&ctx);

    // Register VCF.
    let vcf = VcfTableProvider::new(
        sampled_vcf.display().to_string(),
        Some(vec![]),
        Some(vec![]),
        None,
        false,
    )?;
    ctx.register_table("sampled_vcf", Arc::new(vcf))?;

    // Register cache parquet tables.
    let variation_path = cache_dir.join(&variation_file);
    ctx.register_parquet(
        "variation_cache",
        variation_path.display().to_string().as_str(),
        ParquetReadOptions::default(),
    )
    .await?;

    let mut options = Vec::new();

    // Register context tables.
    let context_tables = [
        ("transcript", "transcripts_table"),
        ("exon", "exons_table"),
        ("translation", "translations_table"),
        ("regulatory", "regulatory_table"),
        ("motif", "motif_table"),
    ];
    for (file_stem, json_key) in &context_tables {
        if let Ok(filename) = find_parquet_stem(&cache_dir, file_stem, &base) {
            let table_name = filename.trim_end_matches(".parquet");
            let path = cache_dir.join(&filename);
            ctx.register_parquet(
                table_name,
                path.display().to_string().as_str(),
                ParquetReadOptions::default(),
            )
            .await?;
            options.push(format!("\"{}\":\"{}\"", json_key, sql_literal(table_name)));
            eprintln!("[PREP] Registered {} -> {}", json_key, filename);
        }
    }

    options.push("\"extended_probes\":true".to_string());
    options.push("\"check_existing\":true".to_string());
    options.push("\"af\":true".to_string());
    options.push("\"af_1kg\":true".to_string());
    options.push("\"af_gnomade\":true".to_string());
    options.push("\"af_gnomadg\":true".to_string());
    options.push("\"max_af\":true".to_string());
    options.push("\"pubmed\":true".to_string());

    if everything {
        options.push("\"everything\":true".to_string());
    }
    if hgvs || everything {
        let fasta = reference_fasta_path.as_deref().ok_or_else(|| {
            DataFusionError::Execution(
                "--hgvs/--everything requires --reference-fasta-path=<path>".to_string(),
            )
        })?;
        if !everything {
            options.push("\"hgvs\":true".to_string());
        }
        options.push(format!(
            "\"reference_fasta_path\":\"{}\"",
            sql_literal(fasta)
        ));
    }

    let options_json = format!("{{{}}}", options.join(","));
    eprintln!("[PREP] options_json: {}", options_json);
    eprintln!();

    // Run annotation with profiling.
    let sql = format!(
        "SELECT * FROM annotate_vep('sampled_vcf', '{}', '{}', '{}')",
        sql_literal(variation_path.display().to_string().as_str()),
        sql_literal("parquet"),
        sql_literal(&options_json)
    );

    eprintln!("=== Running annotation ===");
    let t_total = Instant::now();
    let batches = ctx.sql(&sql).await?.collect().await?;
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    let total_ms = t_total.elapsed().as_secs_f64() * 1000.0;
    eprintln!();
    eprintln!("=== Results ===");
    eprintln!("Output rows:  {}", total_rows);
    eprintln!(
        "Total time:   {:.1}ms ({:.2}s)",
        total_ms,
        total_ms / 1000.0
    );
    eprintln!(
        "Throughput:   {:.0} variants/sec",
        total_rows as f64 / (total_ms / 1000.0)
    );

    // Print schema.
    if let Some(batch) = batches.first() {
        eprintln!("Output cols:  {}", batch.schema().fields().len());
    }

    Ok(())
}

fn find_parquet(dir: &Path, stem_contains: &str) -> Result<String> {
    for entry in std::fs::read_dir(dir).map_err(|e| {
        DataFusionError::Execution(format!("cannot read dir {}: {e}", dir.display()))
    })? {
        let entry =
            entry.map_err(|e| DataFusionError::Execution(format!("dir entry error: {e}")))?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.contains(stem_contains) && name.ends_with(".parquet") {
            return Ok(name);
        }
    }
    Err(DataFusionError::Execution(format!(
        "no parquet file containing '{}' found in {}",
        stem_contains,
        dir.display()
    )))
}

fn find_parquet_stem(dir: &Path, stem: &str, base: &str) -> Result<String> {
    // Try "{base}_{stem}_*_vep.parquet" pattern first, then "{base}_{stem}.parquet".
    for entry in std::fs::read_dir(dir).map_err(|e| {
        DataFusionError::Execution(format!("cannot read dir {}: {e}", dir.display()))
    })? {
        let entry =
            entry.map_err(|e| DataFusionError::Execution(format!("dir entry error: {e}")))?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with(base) && name.contains(stem) && name.ends_with(".parquet") {
            return Ok(name);
        }
    }
    Err(DataFusionError::Execution(format!(
        "no parquet file for '{}' with base '{}' in {}",
        stem,
        base,
        dir.display()
    )))
}

fn extract_base(variation_filename: &str) -> String {
    // "115_GRCh38_variation_1_vep.parquet" -> "115_GRCh38"
    if let Some(idx) = variation_filename.find("_variation") {
        variation_filename[..idx].to_string()
    } else {
        variation_filename.trim_end_matches(".parquet").to_string()
    }
}
