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
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::golden_benchmark::{
    compare_annotations, ensure_local_copy, normalize_chrom, parse_vep_vcf_annotations,
    sample_gz_vcf_first_n, ComparisonReport, VariantAnnotation, VariantKey,
    DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ, DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ_TBI,
    DEFAULT_EXTERNAL_VEP_CACHE_DIR, DEFAULT_LOCAL_HG002_CHR22_VCF_GZ,
    DEFAULT_LOCAL_HG002_CHR22_VCF_GZ_TBI,
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
}

impl Args {
    fn parse() -> Self {
        let args: Vec<String> = std::env::args().collect();
        Self {
            source_vcf_gz: PathBuf::from(
                args.get(1)
                    .cloned()
                    .unwrap_or_else(|| DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ.to_string()),
            ),
            cache_source: args
                .get(2)
                .cloned()
                .unwrap_or_else(|| DEFAULT_CACHE_SOURCE.to_string()),
            backend: args
                .get(3)
                .cloned()
                .unwrap_or_else(|| DEFAULT_BACKEND.to_string()),
            sample_limit: args
                .get(4)
                .and_then(|s| s.parse::<usize>().ok())
                .unwrap_or(DEFAULT_SAMPLE_LIMIT),
            vep_cache_dir: PathBuf::from(
                args.get(5)
                    .cloned()
                    .unwrap_or_else(|| DEFAULT_EXTERNAL_VEP_CACHE_DIR.to_string()),
            ),
            local_copy_vcf_gz: PathBuf::from(
                args.get(6)
                    .cloned()
                    .unwrap_or_else(|| DEFAULT_LOCAL_HG002_CHR22_VCF_GZ.to_string()),
            ),
            work_dir: PathBuf::from(
                args.get(7)
                    .cloned()
                    .unwrap_or_else(|| DEFAULT_WORK_DIR.to_string()),
            ),
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

    fs::create_dir_all(&args.work_dir).map_err(io_err)?;

    let source_tbi = PathBuf::from(DEFAULT_EXTERNAL_HG002_CHR22_VCF_GZ_TBI);
    let local_tbi = PathBuf::from(DEFAULT_LOCAL_HG002_CHR22_VCF_GZ_TBI);
    ensure_local_copy(
        &args.source_vcf_gz,
        &args.local_copy_vcf_gz,
        Some(&source_tbi),
        Some(&local_tbi),
    )?;

    let sampled_vcf = args.work_dir.join(format!("HG002_chr22_{}.vcf", args.sample_limit));
    let golden_vcf = args
        .work_dir
        .join(format!("HG002_chr22_{}_vep115_golden.vcf", args.sample_limit));
    let report_path = args
        .work_dir
        .join(format!("HG002_chr22_{}_comparison_report.txt", args.sample_limit));

    let sampling_start = Instant::now();
    let sampled_rows = sample_gz_vcf_first_n(&args.local_copy_vcf_gz, &sampled_vcf, args.sample_limit)?;
    let sampling_elapsed = sampling_start.elapsed().as_secs_f64();
    println!(
        "sampled rows: {} (elapsed {:.3}s) -> {}",
        sampled_rows,
        sampling_elapsed,
        sampled_vcf.display()
    );

    let docker_start = Instant::now();
    run_vep_docker(&sampled_vcf, &golden_vcf, &args.work_dir, &args.vep_cache_dir)?;
    let docker_elapsed = docker_start.elapsed().as_secs_f64();
    println!(
        "golden VEP output generated (elapsed {:.3}s) -> {}",
        docker_elapsed,
        golden_vcf.display()
    );

    let parse_start = Instant::now();
    let golden_annotations = parse_vep_vcf_annotations(&golden_vcf)?;
    let parse_elapsed = parse_start.elapsed().as_secs_f64();
    println!(
        "golden parsed rows: {} (elapsed {:.3}s)",
        golden_annotations.len(),
        parse_elapsed
    );

    let ours_start = Instant::now();
    let ours_annotations = run_annotate_vep(&sampled_vcf, &args.cache_source, &args.backend).await?;
    let ours_elapsed = ours_start.elapsed().as_secs_f64();
    println!(
        "annotate_vep rows: {} (elapsed {:.3}s)",
        ours_annotations.len(),
        ours_elapsed
    );

    let report = compare_annotations(&golden_annotations, &ours_annotations);
    print_report(&report);
    write_report(
        &report_path,
        &report,
        sampled_rows,
        sampling_elapsed,
        docker_elapsed,
        parse_elapsed,
        ours_elapsed,
    )?;
    println!("report written to {}", report_path.display());

    Ok(())
}

fn run_vep_docker(
    sampled_vcf: &Path,
    golden_vcf: &Path,
    work_dir: &Path,
    vep_cache_dir: &Path,
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

    let output = Command::new("docker")
        .arg("run")
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
        .arg("--fork")
        .arg("4")
        .arg("--force_overwrite")
        .arg("--no_stats")
        .output()
        .map_err(io_err)?;

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

async fn run_annotate_vep(sampled_vcf: &Path, cache_source: &str, backend: &str) -> Result<Vec<VariantAnnotation>> {
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

    let sql = format!(
        "SELECT chrom, start, ref, alt, csq, most_severe_consequence FROM annotate_vep('sampled_vcf', '{}', '{}')",
        sql_literal(cache_source),
        sql_literal(backend)
    );
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

fn print_report(report: &ComparisonReport) {
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
}

fn write_report(
    report_path: &Path,
    report: &ComparisonReport,
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
    Ok(())
}

fn io_err<E: std::fmt::Display>(e: E) -> DataFusionError {
    DataFusionError::Execution(e.to_string())
}
