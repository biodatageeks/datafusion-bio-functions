use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{
    Array, ArrayRef, Int32Array, Int64Array, LargeStringArray, StringArray, StringViewArray,
    UInt32Array, UInt64Array,
};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{ParquetReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::golden_benchmark::{
    ComparisonReport, CsqFieldReport, CsqUnmatchedReport, DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ,
    DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ_TBI, DEFAULT_EXTERNAL_VEP_CACHE_DIR,
    DEFAULT_LOCAL_HG002_CHR22_VCF_GZ, DEFAULT_LOCAL_HG002_CHR22_VCF_GZ_TBI, TermComparisonReport,
    VariantAnnotation, VariantDiscrepancy, VariantKey, collect_discrepancies,
    compare_annotation_terms, compare_annotations, compare_csq_fields, diagnose_unmatched_csq,
    ensure_local_copy, normalize_chrom, parse_vep_vcf_annotations, sample_gz_vcf_first_n,
};
use datafusion_bio_function_vep::register_vep_functions;

const DEFAULT_CACHE_SOURCE: &str = "/Users/mwiewior/research/data/vep/115_GRCh38_variants.parquet";
const DEFAULT_BACKEND: &str = "parquet";
const DEFAULT_SAMPLE_LIMIT: usize = 1000;
const DEFAULT_WORK_DIR: &str = "/tmp/annotate_vep_golden_bench";
#[derive(Debug, Clone)]
struct Args {
    source_vcf_gz: PathBuf,
    cache_source: String,
    backend: String,
    sample_limit: usize,
    vep_cache_dir: PathBuf,
    local_copy_vcf_gz: PathBuf,
    work_dir: PathBuf,
    /// Directory containing context parquet files (transcripts, exons, etc.).
    context_dir: Option<PathBuf>,
    /// Whether to use VEP --merged flag (for merged Ensembl+RefSeq cache).
    merged: bool,
}

impl Args {
    fn parse() -> Self {
        let args: Vec<String> = std::env::args().collect();

        // Check for --merged flag anywhere in args.
        let merged = args.iter().any(|a| a == "--merged");
        // Filter out flags for positional parsing.
        let positional: Vec<&String> = args.iter().filter(|a| !a.starts_with("--")).collect();

        Self {
            source_vcf_gz: PathBuf::from(
                positional
                    .get(1)
                    .map(|s| s.as_str())
                    .unwrap_or(DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ),
            ),
            cache_source: positional
                .get(2)
                .map(|s| s.to_string())
                .unwrap_or_else(|| DEFAULT_CACHE_SOURCE.to_string()),
            backend: positional
                .get(3)
                .map(|s| s.to_string())
                .unwrap_or_else(|| DEFAULT_BACKEND.to_string()),
            sample_limit: positional
                .get(4)
                .and_then(|s| s.parse::<usize>().ok())
                .unwrap_or(DEFAULT_SAMPLE_LIMIT),
            vep_cache_dir: PathBuf::from(
                positional
                    .get(5)
                    .map(|s| s.as_str())
                    .unwrap_or(DEFAULT_EXTERNAL_VEP_CACHE_DIR),
            ),
            local_copy_vcf_gz: PathBuf::from(
                positional
                    .get(6)
                    .map(|s| s.as_str())
                    .unwrap_or(DEFAULT_LOCAL_HG002_CHR22_VCF_GZ),
            ),
            work_dir: PathBuf::from(
                positional
                    .get(7)
                    .map(|s| s.as_str())
                    .unwrap_or(DEFAULT_WORK_DIR),
            ),
            context_dir: positional.get(8).map(|s| PathBuf::from(s.as_str())),
            merged,
        }
    }
}

fn sql_literal(value: &str) -> String {
    value.replace('\'', "''")
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse();

    println!("annotate_vep golden benchmark");
    println!("  source_vcf_gz: {}", args.source_vcf_gz.display());
    println!("  local_copy_vcf_gz: {}", args.local_copy_vcf_gz.display());
    println!("  sample_limit: {}", args.sample_limit);
    println!("  cache_source: {}", args.cache_source);
    println!("  backend: {}", args.backend);
    println!("  vep_cache_dir: {}", args.vep_cache_dir.display());
    println!("  work_dir: {}", args.work_dir.display());
    println!(
        "  context_dir: {}",
        args.context_dir
            .as_deref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "(auto-discover)".to_string())
    );
    println!("  merged: {}", args.merged);

    fs::create_dir_all(&args.work_dir).map_err(io_err)?;

    let source_tbi = PathBuf::from(DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ_TBI);
    let local_tbi = PathBuf::from(DEFAULT_LOCAL_HG002_CHR22_VCF_GZ_TBI);
    ensure_local_copy(
        &args.source_vcf_gz,
        &args.local_copy_vcf_gz,
        Some(&source_tbi),
        Some(&local_tbi),
    )?;

    let sampled_vcf = args
        .work_dir
        .join(format!("HG002_chr22_{}.vcf", args.sample_limit));
    let golden_vcf = args.work_dir.join(format!(
        "HG002_chr22_{}_vep115_golden.vcf",
        args.sample_limit
    ));
    let report_path = args.work_dir.join(format!(
        "HG002_chr22_{}_comparison_report.txt",
        args.sample_limit
    ));
    let diff_path = args.work_dir.join(format!(
        "HG002_chr22_{}_discrepancies.txt",
        args.sample_limit
    ));

    let sampling_start = Instant::now();
    let sampled_rows =
        sample_gz_vcf_first_n(&args.local_copy_vcf_gz, &sampled_vcf, args.sample_limit)?;
    let sampling_elapsed = sampling_start.elapsed().as_secs_f64();
    println!(
        "sampled rows: {} (elapsed {:.3}s) -> {}",
        sampled_rows,
        sampling_elapsed,
        sampled_vcf.display()
    );

    // Decompose multi-allelic sites so both VEP and our engine see single-alt records.
    let normalized_vcf = args
        .work_dir
        .join(format!("HG002_chr22_{}_norm.vcf", args.sample_limit));
    let norm_rows = normalize_vcf(&sampled_vcf, &normalized_vcf)?;
    println!(
        "normalized: {} -> {} rows (multi-allelic decomposed) -> {}",
        sampled_rows,
        norm_rows,
        normalized_vcf.display()
    );

    let docker_elapsed = if golden_vcf.exists() {
        println!(
            "reusing existing golden VEP output -> {}",
            golden_vcf.display()
        );
        0.0
    } else {
        let docker_start = Instant::now();
        run_vep_docker(
            &normalized_vcf,
            &golden_vcf,
            &args.work_dir,
            &args.vep_cache_dir,
            args.merged,
        )?;
        let elapsed = docker_start.elapsed().as_secs_f64();
        println!(
            "golden VEP output generated (elapsed {:.3}s) -> {}",
            elapsed,
            golden_vcf.display()
        );
        elapsed
    };

    let parse_start = Instant::now();
    let golden_annotations = parse_vep_vcf_annotations(&golden_vcf)?;
    let parse_elapsed = parse_start.elapsed().as_secs_f64();
    println!(
        "golden parsed rows: {} (elapsed {:.3}s)",
        golden_annotations.len(),
        parse_elapsed
    );

    // Build options_json if context_dir is provided.
    let options_json = build_options_json(&args);

    let effective_context_dir = args
        .context_dir
        .as_deref()
        .or_else(|| Path::new(&args.cache_source).parent());
    let ours_start = Instant::now();
    let ours_annotations = run_annotate_vep(
        &normalized_vcf,
        &args.cache_source,
        &args.backend,
        options_json.as_deref(),
        effective_context_dir,
    )
    .await?;
    let ours_elapsed = ours_start.elapsed().as_secs_f64();
    println!(
        "annotate_vep rows: {} (elapsed {:.3}s)",
        ours_annotations.len(),
        ours_elapsed
    );

    let report = compare_annotations(&golden_annotations, &ours_annotations);
    let term_report = compare_annotation_terms(&golden_annotations, &ours_annotations);
    let csq_field_report = compare_csq_fields(&golden_annotations, &ours_annotations);
    let discrepancies = collect_discrepancies(&golden_annotations, &ours_annotations);
    let unmatched_report = diagnose_unmatched_csq(&golden_annotations, &ours_annotations);

    print_report(&report, &term_report);
    print_csq_field_report(&csq_field_report);
    print_unmatched_report(&unmatched_report);
    print_discrepancy_summary(&discrepancies, &report);

    write_report(
        &report_path,
        &report,
        &term_report,
        &csq_field_report,
        sampled_rows,
        sampling_elapsed,
        docker_elapsed,
        parse_elapsed,
        ours_elapsed,
    )?;
    write_discrepancies(&diff_path, &discrepancies)?;
    println!("report written to {}", report_path.display());
    println!(
        "discrepancies written to {} ({} variants)",
        diff_path.display(),
        discrepancies.len()
    );

    Ok(())
}

/// Build options_json from context_dir by discovering parquet files.
///
/// Scans context_dir for files matching known patterns and builds the
/// JSON string to pass as the 4th argument to `annotate_vep()`.
fn build_options_json(args: &Args) -> Option<String> {
    // Use explicit context_dir, or derive from cache_source parent directory.
    let context_dir = args
        .context_dir
        .as_deref()
        .or_else(|| Path::new(&args.cache_source).parent())
        .filter(|p| p.is_dir())?;

    // Derive base name from cache_source (e.g. "115_GRCh38_variation_22" -> "115_GRCh38")
    // by stripping "_variation*" or "_variants*" suffix.
    let cache_stem = Path::new(&args.cache_source)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("");

    let base = if let Some(idx) = cache_stem.find("_variation") {
        &cache_stem[..idx]
    } else if let Some(idx) = cache_stem.find("_variants") {
        &cache_stem[..idx]
    } else {
        cache_stem
    };

    // Derive suffix (e.g. "_22" from "115_GRCh38_variation_22").
    let suffix = cache_stem
        .strip_prefix(base)
        .and_then(|rest| rest.strip_prefix("_variation"))
        .or_else(|| {
            cache_stem
                .strip_prefix(base)
                .and_then(|rest| rest.strip_prefix("_variants"))
        })
        .unwrap_or("");

    // Build table name -> parquet path mapping for known context tables.
    let table_specs: Vec<(&str, &str)> = vec![
        ("transcripts_table", "transcript"),
        ("exons_table", "exon"),
        ("translations_table", "translation"),
        ("regulatory_table", "regulatory"),
        ("motif_table", "motif"),
    ];

    let mut entries = Vec::new();
    for (json_key, file_stem) in &table_specs {
        let parquet_name = format!("{base}_{file_stem}{suffix}.parquet");
        let parquet_path = context_dir.join(&parquet_name);
        if parquet_path.exists() {
            // Table name used for DataFusion registration (no .parquet extension).
            let table_name = format!("{base}_{file_stem}{suffix}");
            entries.push(format!("\"{json_key}\":\"{table_name}\""));
        }
    }

    // Disable extended probes for faster equi-join.
    entries.push("\"extended_probes\":false".to_string());

    if args.merged {
        entries.push("\"merged\":true".to_string());
    }

    if entries.is_empty() {
        return None;
    }

    Some(format!("{{{}}}", entries.join(",")))
}

/// Register context parquet files from context_dir in the DataFusion session.
async fn register_context_tables(
    ctx: &SessionContext,
    context_dir: &Path,
    options_json: &str,
) -> Result<()> {
    // Parse table names from options_json and register corresponding parquet files.
    let keys = [
        "transcripts_table",
        "exons_table",
        "translations_table",
        "regulatory_table",
        "motif_table",
    ];
    for key in &keys {
        if let Some(table_name) = parse_json_string_value(options_json, key) {
            let parquet_path = context_dir.join(format!("{table_name}.parquet"));
            if parquet_path.exists() {
                ctx.register_parquet(
                    &table_name,
                    parquet_path.display().to_string().as_str(),
                    ParquetReadOptions::default(),
                )
                .await?;
            }
        }
    }
    Ok(())
}

/// Decompose multi-allelic sites using `bcftools norm -m -both`.
/// Returns the number of data rows in the normalized VCF.
fn normalize_vcf(input: &Path, output: &Path) -> Result<usize> {
    let result = Command::new("bcftools")
        .arg("norm")
        .arg("-m")
        .arg("-both")
        .arg("-o")
        .arg(output.as_os_str())
        .arg(input.as_os_str())
        .output()
        .map_err(io_err)?;

    if !result.status.success() {
        return Err(DataFusionError::Execution(format!(
            "bcftools norm failed: {}",
            String::from_utf8_lossy(&result.stderr)
        )));
    }

    // Count data lines in the output VCF.
    let contents = fs::read_to_string(output).map_err(io_err)?;
    let count = contents.lines().filter(|l| !l.starts_with('#')).count();
    Ok(count)
}

/// Minimal JSON string value parser (matches the one in AnnotateProvider).
fn parse_json_string_value(json: &str, key: &str) -> Option<String> {
    let needle = format!("\"{key}\"");
    let start = json.find(&needle)?;
    let rest = &json[start + needle.len()..];
    let colon = rest.find(':')?;
    let after_colon = rest[colon + 1..].trim_start();
    let after_quote = after_colon.strip_prefix('"')?;
    let end_quote = after_quote.find('"')?;
    let value = &after_quote[..end_quote];
    if value.is_empty() {
        return None;
    }
    Some(value.to_string())
}

fn run_vep_docker(
    sampled_vcf: &Path,
    golden_vcf: &Path,
    work_dir: &Path,
    vep_cache_dir: &Path,
    merged: bool,
) -> Result<()> {
    let sampled_name = sampled_vcf
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| DataFusionError::Execution("invalid sampled VCF file name".to_string()))?;
    let golden_name = golden_vcf
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| DataFusionError::Execution("invalid golden VCF file name".to_string()))?;

    let volume_cache = format!("{}:/opt/vep/.vep", vep_cache_dir.display());
    let volume_work = format!("{}:/work", work_dir.display());

    let mut cmd = Command::new("docker");
    cmd.arg("run")
        .arg("--rm")
        .arg("-v")
        .arg(volume_cache)
        .arg("-v")
        .arg(volume_work)
        .arg("ensemblorg/ensembl-vep:release_115.2")
        .arg("vep")
        .arg("--dir")
        .arg("/opt/vep/.vep")
        .arg("--cache")
        .arg("--offline")
        .arg("--assembly")
        .arg("GRCh38")
        .arg("--input_file")
        .arg(format!("/work/{sampled_name}"))
        .arg("--output_file")
        .arg(format!("/work/{golden_name}"))
        .arg("--vcf")
        .arg("--force_overwrite")
        .arg("--no_stats");

    // Enable regulatory feature annotations (VEP skips them by default).
    cmd.arg("--regulatory");

    // Use all default VEP CSQ fields for full per-field comparison.
    cmd.arg("--fields")
        .arg("Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS,SOURCE");

    if merged {
        cmd.arg("--merged");
        // Merged caches auto-enable --use_transcript_ref which requires FASTA.
        // Use --use_given_ref to skip the FASTA requirement.
        cmd.arg("--use_given_ref");
    }

    let output = cmd.output().map_err(io_err)?;

    if !output.status.success() {
        return Err(DataFusionError::Execution(format!(
            "docker VEP command failed with status {}.\nstdout:\n{}\nstderr:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        )));
    }

    if !golden_vcf.exists() {
        return Err(DataFusionError::Execution(format!(
            "golden output file was not created: {}",
            golden_vcf.display()
        )));
    }

    Ok(())
}

async fn run_annotate_vep(
    sampled_vcf: &Path,
    cache_source: &str,
    backend: &str,
    options_json: Option<&str>,
    context_dir: Option<&Path>,
) -> Result<Vec<VariantAnnotation>> {
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    register_vep_functions(&ctx);

    let vcf = VcfTableProvider::new(
        sampled_vcf.display().to_string(),
        Some(vec![]),
        Some(vec![]),
        None,
        false,
    )?;
    ctx.register_table("sampled_vcf", Arc::new(vcf))?;

    // Register context parquet tables if options_json and context_dir are provided.
    if let (Some(opts), Some(dir)) = (options_json, context_dir) {
        register_context_tables(&ctx, dir, opts).await?;
    }

    let sql = if let Some(opts) = options_json {
        format!(
            "SELECT chrom, start, ref, alt, csq, most_severe_consequence FROM annotate_vep('sampled_vcf', '{}', '{}', '{}')",
            sql_literal(cache_source),
            sql_literal(backend),
            sql_literal(opts)
        )
    } else {
        format!(
            "SELECT chrom, start, ref, alt, csq, most_severe_consequence FROM annotate_vep('sampled_vcf', '{}', '{}')",
            sql_literal(cache_source),
            sql_literal(backend)
        )
    };
    let batches = ctx.sql(&sql).await?.collect().await?;

    let mut out = Vec::new();
    for batch in &batches {
        let chrom_col = batch.column_by_name("chrom").ok_or_else(|| {
            DataFusionError::Execution("annotate_vep output missing column chrom".to_string())
        })?;
        let start_col = batch.column_by_name("start").ok_or_else(|| {
            DataFusionError::Execution("annotate_vep output missing column start".to_string())
        })?;
        let ref_col = batch.column_by_name("ref").ok_or_else(|| {
            DataFusionError::Execution("annotate_vep output missing column ref".to_string())
        })?;
        let alt_col = batch.column_by_name("alt").ok_or_else(|| {
            DataFusionError::Execution("annotate_vep output missing column alt".to_string())
        })?;
        let csq_col = batch.column_by_name("csq").ok_or_else(|| {
            DataFusionError::Execution("annotate_vep output missing column csq".to_string())
        })?;
        let most_col = batch
            .column_by_name("most_severe_consequence")
            .ok_or_else(|| {
                DataFusionError::Execution(
                    "annotate_vep output missing column most_severe_consequence".to_string(),
                )
            })?;

        for row in 0..batch.num_rows() {
            let chrom = string_at(chrom_col, row)?;
            let pos = int64_at(start_col, row)?;
            let ref_allele = string_at(ref_col, row)?;
            let alt_alleles = string_at(alt_col, row)?;

            let key = VariantKey {
                chrom: normalize_chrom(&chrom.unwrap_or_default()),
                pos: pos.unwrap_or_default(),
                ref_allele: ref_allele.unwrap_or_default(),
                alt_alleles: alt_alleles.unwrap_or_default(),
            };
            out.push(VariantAnnotation {
                key,
                csq: string_at(csq_col, row)?,
                most_severe_consequence: string_at(most_col, row)?,
            });
        }
    }
    Ok(out)
}

fn string_at(col: &ArrayRef, row: usize) -> Result<Option<String>> {
    if row >= col.len() || col.is_null(row) {
        return Ok(None);
    }
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        return Ok(Some(arr.value(row).to_string()));
    }
    if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        return Ok(Some(arr.value(row).to_string()));
    }
    if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
        return Ok(Some(arr.value(row).to_string()));
    }
    Err(DataFusionError::Execution(format!(
        "expected string-like column, got {:?}",
        col.data_type()
    )))
}

fn int64_at(col: &ArrayRef, row: usize) -> Result<Option<i64>> {
    if row >= col.len() || col.is_null(row) {
        return Ok(None);
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return Ok(Some(arr.value(row)));
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
        return Ok(Some(arr.value(row) as i64));
    }
    if let Some(arr) = col.as_any().downcast_ref::<UInt64Array>() {
        return Ok(Some(arr.value(row) as i64));
    }
    if let Some(arr) = col.as_any().downcast_ref::<UInt32Array>() {
        return Ok(Some(arr.value(row) as i64));
    }
    Err(DataFusionError::Execution(format!(
        "expected integer start column, got {:?}",
        col.data_type()
    )))
}

fn print_report(report: &ComparisonReport, term_report: &TermComparisonReport) {
    println!("\ncomparison report");
    println!("  golden_rows: {}", report.golden_rows);
    println!("  ours_rows: {}", report.ours_rows);
    println!("  intersection_rows: {}", report.intersection_rows);
    println!("  missing_in_ours: {}", report.missing_in_ours);
    println!("  extra_in_ours: {}", report.extra_in_ours);
    println!("  golden_with_csq: {}", report.golden_with_csq);
    println!("  ours_with_csq: {}", report.ours_with_csq);
    println!("  csq_exact_matches: {}", report.csq_exact_matches);
    println!(
        "  most_severe_exact_matches: {}",
        report.most_severe_exact_matches
    );
    if report.intersection_rows > 0 {
        println!(
            "  most_severe_accuracy: {:.1}%",
            report.most_severe_exact_matches as f64 / report.intersection_rows as f64 * 100.0
        );
    }
    println!("  term_comparable_rows: {}", term_report.comparable_rows);
    println!(
        "  term_golden_with_terms: {}",
        term_report.golden_with_terms
    );
    println!("  term_ours_with_terms: {}", term_report.ours_with_terms);
    println!(
        "  term_set_exact_matches: {}",
        term_report.term_set_exact_matches
    );
    if term_report.comparable_rows > 0 {
        println!(
            "  term_set_accuracy: {:.1}%",
            term_report.term_set_exact_matches as f64 / term_report.comparable_rows as f64 * 100.0
        );
    }
    println!(
        "  golden_term_subset_matches: {}",
        term_report.golden_term_subset_matches
    );
}

fn print_csq_field_report(report: &CsqFieldReport) {
    if report.total_entries_compared == 0 {
        println!("\ncsq field report: no entries to compare");
        return;
    }
    println!(
        "\ncsq per-field accuracy ({} entries compared):",
        report.total_entries_compared
    );
    for (i, name) in report.field_names.iter().enumerate() {
        let matches = report.field_match_counts[i];
        let pct = matches as f64 / report.total_entries_compared as f64 * 100.0;
        println!(
            "  {:<20} {}/{} ({:.1}%)",
            name, matches, report.total_entries_compared, pct
        );
    }
}

fn print_unmatched_report(report: &CsqUnmatchedReport) {
    println!("\ncsq entry matching diagnostic:");
    println!("  golden_total_entries: {}", report.golden_total);
    println!("  ours_total_entries: {}", report.ours_total);
    println!("  matched_by_allele_feature: {}", report.matched);
    println!(
        "  golden_only (unmatched): {}",
        report.golden_total.saturating_sub(report.matched)
    );
    println!(
        "  ours_only (extra): {}",
        report.ours_total.saturating_sub(report.matched)
    );

    if !report.golden_only_by_ft.is_empty() {
        println!("  golden_only by Feature_type:");
        let mut sorted: Vec<_> = report.golden_only_by_ft.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));
        for (ft, count) in &sorted {
            let label = if ft.is_empty() { "(empty)" } else { ft.as_str() };
            println!("    {:<25} {}", label, count);
        }
    }
    if !report.ours_only_by_ft.is_empty() {
        println!("  ours_only by Feature_type:");
        let mut sorted: Vec<_> = report.ours_only_by_ft.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));
        for (ft, count) in &sorted {
            let label = if ft.is_empty() { "(empty)" } else { ft.as_str() };
            println!("    {:<25} {}", label, count);
        }
    }

    if !report.golden_only_sample.is_empty() {
        println!("  golden_only samples (allele|feature|feature_type):");
        for (a, f, ft) in &report.golden_only_sample {
            println!("    {}|{}|{}", a, f, ft);
        }
    }
    if !report.ours_only_sample.is_empty() {
        println!("  ours_only samples (allele|feature|feature_type):");
        for (a, f, ft) in &report.ours_only_sample {
            println!("    {}|{}|{}", a, f, ft);
        }
    }
}

fn print_discrepancy_summary(discrepancies: &[VariantDiscrepancy], report: &ComparisonReport) {
    let total = discrepancies.len();
    let most_severe_mismatches = discrepancies
        .iter()
        .filter(|d| !d.most_severe_match)
        .count();
    let term_mismatches = discrepancies.iter().filter(|d| !d.term_set_match).count();
    let missing = discrepancies
        .iter()
        .filter(|d| d.ours_most_severe.is_none())
        .count();

    println!("\ndiscrepancy summary");
    println!("  total_discrepancies: {total}");
    println!("  most_severe_mismatches: {most_severe_mismatches}");
    println!("  term_set_mismatches: {term_mismatches}");
    println!("  missing_in_ours: {missing}");

    // Print first 20 discrepancies as preview.
    let preview_limit = 20;
    let shown = discrepancies.len().min(preview_limit);
    if shown > 0 {
        println!("\nfirst {shown} discrepancies:");
        for d in discrepancies.iter().take(preview_limit) {
            let golden_ms = d.golden_most_severe.as_deref().unwrap_or("(none)");
            let ours_ms = d.ours_most_severe.as_deref().unwrap_or("(none)");
            let ms_marker = if d.most_severe_match { "=" } else { "!" };
            let ts_marker = if d.term_set_match { "=" } else { "!" };

            println!(
                "  {}:{} {}/{} most_severe[{ms_marker}]: {} vs {} terms[{ts_marker}]: {:?} vs {:?}",
                d.key.chrom,
                d.key.pos,
                d.key.ref_allele,
                d.key.alt_alleles,
                golden_ms,
                ours_ms,
                d.golden_terms,
                d.ours_terms,
            );
        }
        if total > preview_limit {
            println!(
                "  ... and {} more (see discrepancies file)",
                total - preview_limit
            );
        }
    }

    // Print accuracy rates.
    if report.intersection_rows > 0 {
        println!(
            "\naccuracy: most_severe={:.1}% term_set={:.1}%",
            (report.intersection_rows - most_severe_mismatches) as f64
                / report.intersection_rows as f64
                * 100.0,
            (report.intersection_rows - term_mismatches) as f64 / report.intersection_rows as f64
                * 100.0,
        );
    }
}

fn write_report(
    report_path: &Path,
    report: &ComparisonReport,
    term_report: &TermComparisonReport,
    csq_field_report: &CsqFieldReport,
    sampled_rows: usize,
    sampling_elapsed_s: f64,
    docker_elapsed_s: f64,
    parse_elapsed_s: f64,
    ours_elapsed_s: f64,
) -> Result<()> {
    let mut f = File::create(report_path).map_err(io_err)?;
    writeln!(f, "annotate_vep vs Ensembl VEP 115 golden benchmark").map_err(io_err)?;
    writeln!(f, "sampled_rows={sampled_rows}").map_err(io_err)?;
    writeln!(f, "sampling_elapsed_s={sampling_elapsed_s:.3}").map_err(io_err)?;
    writeln!(f, "vep_docker_elapsed_s={docker_elapsed_s:.3}").map_err(io_err)?;
    writeln!(f, "golden_parse_elapsed_s={parse_elapsed_s:.3}").map_err(io_err)?;
    writeln!(f, "annotate_vep_elapsed_s={ours_elapsed_s:.3}").map_err(io_err)?;
    writeln!(f, "golden_rows={}", report.golden_rows).map_err(io_err)?;
    writeln!(f, "ours_rows={}", report.ours_rows).map_err(io_err)?;
    writeln!(f, "intersection_rows={}", report.intersection_rows).map_err(io_err)?;
    writeln!(f, "missing_in_ours={}", report.missing_in_ours).map_err(io_err)?;
    writeln!(f, "extra_in_ours={}", report.extra_in_ours).map_err(io_err)?;
    writeln!(f, "golden_with_csq={}", report.golden_with_csq).map_err(io_err)?;
    writeln!(f, "ours_with_csq={}", report.ours_with_csq).map_err(io_err)?;
    writeln!(f, "csq_exact_matches={}", report.csq_exact_matches).map_err(io_err)?;
    writeln!(
        f,
        "most_severe_exact_matches={}",
        report.most_severe_exact_matches
    )
    .map_err(io_err)?;
    if report.intersection_rows > 0 {
        writeln!(
            f,
            "most_severe_accuracy={:.3}",
            report.most_severe_exact_matches as f64 / report.intersection_rows as f64
        )
        .map_err(io_err)?;
    }
    writeln!(f, "term_comparable_rows={}", term_report.comparable_rows).map_err(io_err)?;
    writeln!(
        f,
        "term_golden_with_terms={}",
        term_report.golden_with_terms
    )
    .map_err(io_err)?;
    writeln!(f, "term_ours_with_terms={}", term_report.ours_with_terms).map_err(io_err)?;
    writeln!(
        f,
        "term_set_exact_matches={}",
        term_report.term_set_exact_matches
    )
    .map_err(io_err)?;
    if term_report.comparable_rows > 0 {
        writeln!(
            f,
            "term_set_accuracy={:.3}",
            term_report.term_set_exact_matches as f64 / term_report.comparable_rows as f64
        )
        .map_err(io_err)?;
    }
    writeln!(
        f,
        "golden_term_subset_matches={}",
        term_report.golden_term_subset_matches
    )
    .map_err(io_err)?;
    // CSQ per-field accuracy report.
    writeln!(f, "\n# CSQ per-field accuracy").map_err(io_err)?;
    writeln!(
        f,
        "csq_entries_compared={}",
        csq_field_report.total_entries_compared
    )
    .map_err(io_err)?;
    for (i, name) in csq_field_report.field_names.iter().enumerate() {
        let matches = csq_field_report.field_match_counts[i];
        let pct = if csq_field_report.total_entries_compared > 0 {
            matches as f64 / csq_field_report.total_entries_compared as f64
        } else {
            0.0
        };
        writeln!(
            f,
            "csq_field_{name}={matches}/{} ({pct:.3})",
            csq_field_report.total_entries_compared
        )
        .map_err(io_err)?;
    }
    Ok(())
}

fn write_discrepancies(path: &Path, discrepancies: &[VariantDiscrepancy]) -> Result<()> {
    let mut f = File::create(path).map_err(io_err)?;
    writeln!(f, "# Per-variant discrepancies (golden vs ours)").map_err(io_err)?;
    writeln!(
        f,
        "# chrom:pos ref/alt | most_severe_match | golden_most_severe | ours_most_severe | term_set_match | golden_terms | ours_terms"
    )
    .map_err(io_err)?;
    writeln!(f, "# total: {}", discrepancies.len()).map_err(io_err)?;
    writeln!(f).map_err(io_err)?;

    for d in discrepancies {
        let golden_ms = d.golden_most_severe.as_deref().unwrap_or("(none)");
        let ours_ms = d.ours_most_severe.as_deref().unwrap_or("(none)");
        let golden_terms: Vec<&str> = d.golden_terms.iter().map(|s| s.as_str()).collect();
        let ours_terms: Vec<&str> = d.ours_terms.iter().map(|s| s.as_str()).collect();

        writeln!(
            f,
            "{}:{} {}/{} | ms={} | golden={} | ours={} | ts={} | golden_terms=[{}] | ours_terms=[{}]",
            d.key.chrom,
            d.key.pos,
            d.key.ref_allele,
            d.key.alt_alleles,
            d.most_severe_match,
            golden_ms,
            ours_ms,
            d.term_set_match,
            golden_terms.join(","),
            ours_terms.join(","),
        )
        .map_err(io_err)?;
    }
    Ok(())
}

fn io_err<E: std::fmt::Display>(e: E) -> DataFusionError {
    DataFusionError::Execution(e.to_string())
}
