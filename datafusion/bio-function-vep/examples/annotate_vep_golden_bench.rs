#![allow(
    clippy::collapsible_if,
    clippy::too_many_arguments,
    unused_variables,
    clippy::manual_saturating_arithmetic
)]
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
    ComparisonReport, CsqFieldReport, CsqMultiplicityReport, CsqUnmatchedReport,
    DEFAULT_LOCAL_HG002_CHR22_VCF_GZ, DEFAULT_LOCAL_VEP_CACHE_DIR, TermComparisonReport,
    VariantAnnotation, VariantDiscrepancy, VariantKey, collect_discrepancies,
    compare_annotation_terms, compare_annotations, compare_csq_fields_with_names,
    csq_field_names_for_mode, diagnose_csq_multiplicity, diagnose_unmatched_csq, ensure_local_copy,
    normalize_chrom, parse_vep_vcf_annotations, sample_gz_vcf_first_n,
};
use datafusion_bio_function_vep::register_vep_functions;

const DEFAULT_BACKEND: &str = "parquet";
const DEFAULT_SAMPLE_LIMIT: usize = 1000;
const DEFAULT_WORK_DIR: &str = "/tmp/annotate_vep_golden_bench";
const ENV_SOURCE_VCF_GZ: &str = "ANNOTATE_VEP_SOURCE_VCF_GZ";
const ENV_CACHE_SOURCE: &str = "ANNOTATE_VEP_CACHE_SOURCE";
const ENV_BACKEND: &str = "ANNOTATE_VEP_BACKEND";
const ENV_SAMPLE_LIMIT: &str = "ANNOTATE_VEP_SAMPLE_LIMIT";
const ENV_VEP_CACHE_DIR: &str = "ANNOTATE_VEP_VEP_CACHE_DIR";
const ENV_LOCAL_COPY_VCF_GZ: &str = "ANNOTATE_VEP_LOCAL_COPY_VCF_GZ";
const ENV_WORK_DIR: &str = "ANNOTATE_VEP_WORK_DIR";
const ENV_CONTEXT_DIR: &str = "ANNOTATE_VEP_CONTEXT_DIR";

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
    /// Whether to use VEP --refseq flag (RefSeq-only cache/transcripts).
    refseq: bool,
    /// Whether to use VEP --merged flag (for merged Ensembl+RefSeq cache).
    merged: bool,
    /// Restrict to GENCODE basic transcripts.
    gencode_basic: bool,
    /// Restrict to GENCODE primary transcripts.
    gencode_primary: bool,
    /// Keep all RefSeq transcripts, including CCDS/EST-style rows.
    all_refseq: bool,
    /// Exclude predicted RefSeq transcripts (XM_/XR_).
    exclude_predicted: bool,
    /// Which steps to run: "ensembl", "datafusion", or "ensembl,datafusion" (default: both).
    /// Comparison only runs when both steps are enabled.
    steps: Vec<String>,
    /// Use interval-overlap fallback for shifted indels (--extended-probes).
    extended_probes: bool,
    /// Enable HGVSc/HGVSp generation in both Ensembl VEP and annotate_vep.
    hgvs: bool,
    /// Mirrors VEP `--shift_hgvs [0|1]` when HGVS output is enabled.
    shift_hgvs: Option<bool>,
    /// Reference FASTA required for offline HGVS generation.
    reference_fasta_path: Option<PathBuf>,
    /// Use VEP --everything flag (enables all features, 80-field CSQ schema).
    everything: bool,
    /// Use fjall KV store for variation lookup + SIFT (--use-fjall).
    /// Cache directory must contain a `variation.fjall/` subdirectory.
    use_fjall: bool,
}

fn parse_cli_bool(raw: &str) -> Option<bool> {
    match raw {
        "1" | "true" => Some(true),
        "0" | "false" => Some(false),
        _ => None,
    }
}

fn env_string(key: &str) -> Option<String> {
    std::env::var(key).ok().filter(|value| !value.is_empty())
}

fn positional_or_env_string(
    positional: &[&String],
    index: usize,
    env_key: &str,
    default: Option<&str>,
) -> Option<String> {
    positional
        .get(index)
        .map(|s| s.to_string())
        .or_else(|| env_string(env_key))
        .or_else(|| default.map(str::to_string))
}

impl Args {
    fn parse() -> Result<Self> {
        let args: Vec<String> = std::env::args().collect();

        let refseq = args.iter().any(|a| a == "--refseq");
        // Check for --merged flag anywhere in args.
        let merged = args.iter().any(|a| a == "--merged");
        // Parse --steps=ensembl,datafusion (default: both).
        let steps: Vec<String> = args
            .iter()
            .find(|a| a.starts_with("--steps="))
            .map(|a| {
                a.strip_prefix("--steps=")
                    .unwrap()
                    .split(',')
                    .map(|s| s.trim().to_lowercase())
                    .collect()
            })
            .unwrap_or_else(|| vec!["ensembl".to_string(), "datafusion".to_string()]);
        let hgvs = args.iter().any(|a| a == "--hgvs");
        let shift_hgvs = args
            .iter()
            .find_map(|a| {
                a.strip_prefix("--shift-hgvs=")
                    .or_else(|| a.strip_prefix("--shift_hgvs="))
            })
            .and_then(parse_cli_bool)
            .or_else(|| {
                args.iter()
                    .any(|a| a == "--shift-hgvs" || a == "--shift_hgvs")
                    .then_some(true)
            })
            .or_else(|| {
                args.iter()
                    .any(|a| a == "--no-shift-hgvs" || a == "--no-shift_hgvs")
                    .then_some(false)
            });
        let reference_fasta_path = args
            .iter()
            .find_map(|a| a.strip_prefix("--reference-fasta-path="))
            .map(PathBuf::from);
        // Filter out flags for positional parsing.
        let positional: Vec<&String> = args.iter().filter(|a| !a.starts_with("--")).collect();
        let cache_source = positional_or_env_string(&positional, 2, ENV_CACHE_SOURCE, None)
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "annotate_vep_golden_bench requires positional arg #2 cache_source or {ENV_CACHE_SOURCE}"
                ))
            })?;

        Ok(Self {
            source_vcf_gz: PathBuf::from(
                positional_or_env_string(
                    &positional,
                    1,
                    ENV_SOURCE_VCF_GZ,
                    Some(DEFAULT_LOCAL_HG002_CHR22_VCF_GZ),
                )
                .unwrap(),
            ),
            cache_source,
            backend: positional_or_env_string(&positional, 3, ENV_BACKEND, Some(DEFAULT_BACKEND))
                .unwrap(),
            sample_limit: positional
                .get(4)
                .and_then(|s| s.parse::<usize>().ok())
                .or_else(|| {
                    env_string(ENV_SAMPLE_LIMIT).and_then(|value| value.parse::<usize>().ok())
                })
                .unwrap_or(DEFAULT_SAMPLE_LIMIT),
            vep_cache_dir: PathBuf::from(
                positional_or_env_string(
                    &positional,
                    5,
                    ENV_VEP_CACHE_DIR,
                    Some(DEFAULT_LOCAL_VEP_CACHE_DIR),
                )
                .unwrap(),
            ),
            local_copy_vcf_gz: PathBuf::from(
                positional_or_env_string(
                    &positional,
                    6,
                    ENV_LOCAL_COPY_VCF_GZ,
                    Some(DEFAULT_LOCAL_HG002_CHR22_VCF_GZ),
                )
                .unwrap(),
            ),
            work_dir: PathBuf::from(
                positional_or_env_string(&positional, 7, ENV_WORK_DIR, Some(DEFAULT_WORK_DIR))
                    .unwrap(),
            ),
            context_dir: positional_or_env_string(&positional, 8, ENV_CONTEXT_DIR, None)
                .map(PathBuf::from),
            refseq,
            merged,
            gencode_basic: args
                .iter()
                .any(|a| a == "--gencode-basic" || a == "--gencode_basic"),
            gencode_primary: args
                .iter()
                .any(|a| a == "--gencode-primary" || a == "--gencode_primary"),
            all_refseq: args
                .iter()
                .any(|a| a == "--all-refseq" || a == "--all_refseq"),
            exclude_predicted: args
                .iter()
                .any(|a| a == "--exclude-predicted" || a == "--exclude_predicted"),
            steps,
            extended_probes: args.iter().any(|a| a == "--extended-probes"),
            hgvs,
            shift_hgvs,
            reference_fasta_path,
            everything: args.iter().any(|a| a == "--everything"),
            use_fjall: args.iter().any(|a| a == "--use-fjall"),
        })
    }

    fn run_ensembl(&self) -> bool {
        self.steps.iter().any(|s| s == "ensembl")
    }

    fn run_datafusion(&self) -> bool {
        self.steps.iter().any(|s| s == "datafusion")
    }
}

fn sql_literal(value: &str) -> String {
    value.replace('\'', "''")
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse()?;

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
    println!("  refseq: {}", args.refseq);
    println!("  merged: {}", args.merged);
    println!("  gencode_basic: {}", args.gencode_basic);
    println!("  gencode_primary: {}", args.gencode_primary);
    println!("  all_refseq: {}", args.all_refseq);
    println!("  exclude_predicted: {}", args.exclude_predicted);
    println!("  steps: {:?}", args.steps);
    println!("  extended_probes: {}", args.extended_probes);
    println!("  hgvs: {}", args.hgvs);
    println!(
        "  shift_hgvs: {}",
        args.shift_hgvs
            .map(|value| value.to_string())
            .unwrap_or_else(|| "(default)".to_string())
    );
    println!(
        "  reference_fasta_path: {}",
        args.reference_fasta_path
            .as_deref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "(none)".to_string())
    );
    println!("  everything: {}", args.everything);
    println!("  use_fjall: {}", args.use_fjall);

    fs::create_dir_all(&args.work_dir).map_err(io_err)?;

    // Derive filename stem from source VCF (e.g., "HG002_chr1" from "HG002_chr1.vcf.gz").
    let vcf_stem = args
        .source_vcf_gz
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or("input")
        .strip_suffix(".vcf.gz")
        .or_else(|| {
            args.source_vcf_gz
                .file_name()
                .and_then(|f| f.to_str())
                .unwrap_or("input")
                .strip_suffix(".vcf")
        })
        .unwrap_or("input")
        .to_string();

    let source_tbi = args.source_vcf_gz.with_extension("gz.tbi");
    let local_tbi = args.local_copy_vcf_gz.with_extension("gz.tbi");
    ensure_local_copy(
        &args.source_vcf_gz,
        &args.local_copy_vcf_gz,
        Some(&source_tbi),
        Some(&local_tbi),
    )?;

    let sampled_vcf = args
        .work_dir
        .join(format!("{}_{}.vcf", vcf_stem, args.sample_limit));
    let golden_vcf = args.work_dir.join(format!(
        "{}_{}_vep115_golden.vcf",
        vcf_stem, args.sample_limit
    ));
    let report_path = args.work_dir.join(format!(
        "{}_{}_comparison_report.txt",
        vcf_stem, args.sample_limit
    ));
    let diff_path = args.work_dir.join(format!(
        "{}_{}_discrepancies.txt",
        vcf_stem, args.sample_limit
    ));

    let sampling_start = Instant::now();
    // sample_limit=0 means "all rows" — use usize::MAX as the limit.
    let effective_limit = if args.sample_limit == 0 {
        usize::MAX
    } else {
        args.sample_limit
    };
    let sampled_rows =
        sample_gz_vcf_first_n(&args.local_copy_vcf_gz, &sampled_vcf, effective_limit)?;
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
        .join(format!("{}_{}_norm.vcf", vcf_stem, args.sample_limit));
    let norm_rows = normalize_vcf(&sampled_vcf, &normalized_vcf)?;
    println!(
        "normalized: {} -> {} rows (multi-allelic decomposed) -> {}",
        sampled_rows,
        norm_rows,
        normalized_vcf.display()
    );

    // --- Step: Ensembl VEP (Docker) ---
    let mut docker_elapsed = 0.0;
    let mut golden_annotations = None;
    if args.run_ensembl() {
        docker_elapsed = if golden_vcf.exists() {
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
                args.refseq,
                args.merged,
                args.gencode_basic,
                args.gencode_primary,
                args.all_refseq,
                args.exclude_predicted,
                args.hgvs,
                args.shift_hgvs,
                args.reference_fasta_path.as_deref(),
                args.everything,
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
        let parsed = parse_vep_vcf_annotations(&golden_vcf)?;
        let parse_elapsed = parse_start.elapsed().as_secs_f64();
        println!(
            "golden parsed rows: {} (elapsed {:.3}s)",
            parsed.len(),
            parse_elapsed
        );
        golden_annotations = Some(parsed);
    }

    // --- Step: DataFusion annotate_vep ---
    let mut ours_annotations = None;
    let mut ours_elapsed = 0.0;
    if args.run_datafusion() {
        let options_json = build_options_json(&args)?;
        let effective_context_dir = args
            .context_dir
            .as_deref()
            .or_else(|| Path::new(&args.cache_source).parent());
        let ours_start = Instant::now();
        let ours = run_annotate_vep(
            &normalized_vcf,
            &args.cache_source,
            &args.backend,
            options_json.as_deref(),
            effective_context_dir,
        )
        .await?;
        ours_elapsed = ours_start.elapsed().as_secs_f64();
        println!(
            "annotate_vep rows: {} (elapsed {:.3}s)",
            ours.len(),
            ours_elapsed
        );
        ours_annotations = Some(ours);
    }

    // --- Comparison (only when both steps ran) ---
    if let (Some(golden), Some(ours)) = (&golden_annotations, &ours_annotations) {
        let report = compare_annotations(golden, ours);
        let term_report = compare_annotation_terms(golden, ours);
        let csq_field_names = csq_field_names_for_mode(args.everything, args.refseq, args.merged);
        let csq_field_report = compare_csq_fields_with_names(golden, ours, &csq_field_names);
        let discrepancies = collect_discrepancies(golden, ours);
        let unmatched_report = diagnose_unmatched_csq(golden, ours);
        let multiplicity_report = diagnose_csq_multiplicity(golden, ours);

        print_report(&report, &term_report);
        print_csq_field_report(&csq_field_report);
        print_unmatched_report(&unmatched_report);
        print_multiplicity_report(&multiplicity_report);
        print_discrepancy_summary(&discrepancies, &report);

        write_report(
            &report_path,
            &report,
            &term_report,
            &csq_field_report,
            sampled_rows,
            sampling_elapsed,
            docker_elapsed,
            0.0,
            ours_elapsed,
        )?;
        write_discrepancies(&diff_path, &discrepancies)?;
        println!("report written to {}", report_path.display());
        println!(
            "discrepancies written to {} ({} variants)",
            diff_path.display(),
            discrepancies.len()
        );

        print_summary_table(&report, &term_report, &csq_field_report, ours_elapsed);
    }

    Ok(())
}

/// Build options_json from context_dir by discovering parquet files.
///
/// Scans context_dir for files matching known patterns and builds the
/// JSON string to pass as the 4th argument to `annotate_vep()`.
fn build_options_json(args: &Args) -> Result<Option<String>> {
    // Detect partitioned per-chromosome cache layout.
    let is_partitioned = Path::new(&args.cache_source).join("variation").is_dir();
    if is_partitioned {
        // Partitioned path: AnnotateProvider::scan() handles everything.
        let mut entries = Vec::new();
        entries.push("\"partitioned\":true".to_string());
        entries.push(format!("\"extended_probes\":{}", args.extended_probes));

        if args.everything {
            entries.push("\"everything\":true".to_string());
            let fasta_path = args.reference_fasta_path.as_ref().ok_or_else(|| {
                DataFusionError::Execution(
                    "--everything requires --reference-fasta-path=/path/to/reference.fa[.gz]"
                        .to_string(),
                )
            })?;
            entries.push(format!(
                "\"reference_fasta_path\":\"{}\"",
                sql_literal(fasta_path.to_str().ok_or_else(|| {
                    DataFusionError::Execution(
                        "reference_fasta_path must be valid UTF-8".to_string(),
                    )
                })?)
            ));
        } else {
            if args.hgvs {
                entries.push("\"hgvs\":true".to_string());
                if let Some(shift_hgvs) = args.shift_hgvs {
                    entries.push(format!("\"shift_hgvs\":{shift_hgvs}"));
                }
                let fasta_path = args.reference_fasta_path.as_ref().ok_or_else(|| {
                    DataFusionError::Execution(
                        "--hgvs requires --reference-fasta-path=/path/to/reference.fa[.gz]"
                            .to_string(),
                    )
                })?;
                entries.push(format!(
                    "\"reference_fasta_path\":\"{}\"",
                    sql_literal(fasta_path.to_str().ok_or_else(|| {
                        DataFusionError::Execution(
                            "reference_fasta_path must be valid UTF-8".to_string(),
                        )
                    })?)
                ));
            }
            entries.push("\"check_existing\":true".to_string());
            entries.push("\"af\":true".to_string());
            entries.push("\"af_1kg\":true".to_string());
            entries.push("\"af_gnomade\":true".to_string());
            entries.push("\"af_gnomadg\":true".to_string());
            entries.push("\"max_af\":true".to_string());
            entries.push("\"pubmed\":true".to_string());
        }
        if args.refseq {
            entries.push("\"refseq\":true".to_string());
        }
        if args.merged {
            entries.push("\"merged\":true".to_string());
        }
        if args.gencode_basic {
            entries.push("\"gencode_basic\":true".to_string());
        }
        if args.gencode_primary {
            entries.push("\"gencode_primary\":true".to_string());
        }
        if args.all_refseq {
            entries.push("\"all_refseq\":true".to_string());
        }
        if args.exclude_predicted {
            entries.push("\"exclude_predicted\":true".to_string());
        }
        if args.use_fjall {
            entries.push("\"use_fjall\":true".to_string());
        }
        return Ok(Some(format!("{{{}}}", entries.join(","))));
    }

    // Use explicit context_dir, or derive from cache_source parent directory.
    let context_dir = args
        .context_dir
        .as_deref()
        .or_else(|| Path::new(&args.cache_source).parent())
        .filter(|p| p.is_dir());

    let Some(context_dir) = context_dir else {
        return Ok(None);
    };

    // Derive base name from cache_source or context_dir.
    // For fjall: cache_source is a directory path, so derive base from context_dir parquet files.
    // For parquet: derive from cache_source filename (e.g. "115_GRCh38_variation_22" -> "115_GRCh38").
    let cache_stem = Path::new(&args.cache_source)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("");

    // For fjall backend: cache_source is a directory, not a parquet file.
    // Don't extract base from directory name — use context_dir discovery.
    let is_fjall_dir = Path::new(&args.cache_source).join("keyspaces").is_dir();

    let base = if !is_fjall_dir {
        if let Some(idx) = cache_stem.find("_variation") {
            &cache_stem[..idx]
        } else if let Some(idx) = cache_stem.find("_variants") {
            &cache_stem[..idx]
        } else {
            cache_stem
        }
    } else {
        // Fjall backend: derive base from a parquet file in context_dir.
        let found = std::fs::read_dir(context_dir).ok().and_then(|entries| {
            entries
                .filter_map(|e| e.ok())
                .find(|e| {
                    let name = e.file_name().to_string_lossy().to_string();
                    name.contains("variation") && name.ends_with(".parquet")
                })
                .map(|e| {
                    let name = e.file_name().to_string_lossy().to_string();
                    let stem = name.trim_end_matches(".parquet");
                    if let Some(idx) = stem.find("_variation") {
                        stem[..idx].to_string()
                    } else {
                        stem.to_string()
                    }
                })
        });
        found
            .map(|s| &*Box::leak(s.into_boxed_str()))
            .unwrap_or(cache_stem)
    };

    // Derive suffix (e.g. "_1_vep" from "115_GRCh38_variation_1_vep").
    // For fjall: derive from the variation parquet in context_dir, not from cache_stem.
    let suffix_source = if is_fjall_dir {
        std::fs::read_dir(context_dir)
            .ok()
            .and_then(|entries| {
                entries
                    .filter_map(|e| e.ok())
                    .find(|e| {
                        let n = e.file_name().to_string_lossy().to_string();
                        n.contains("variation") && n.ends_with(".parquet")
                    })
                    .map(|e| {
                        e.file_name()
                            .to_string_lossy()
                            .trim_end_matches(".parquet")
                            .to_string()
                    })
            })
            .unwrap_or_default()
    } else {
        cache_stem.to_string()
    };
    let suffix = suffix_source
        .strip_prefix(base)
        .and_then(|rest| rest.strip_prefix("_variation"))
        .or_else(|| {
            suffix_source
                .strip_prefix(base)
                .and_then(|rest| rest.strip_prefix("_variants"))
        })
        .unwrap_or("");

    // Build table name -> parquet path mapping for known context tables.
    // For each key, try the primary stem first, then fallback stems.
    let table_specs: Vec<(&str, &[&str])> = vec![
        ("transcripts_table", &["transcript"] as &[&str]),
        ("exons_table", &["exon"]),
        ("translations_table", &["translation_core", "translation"]),
        ("regulatory_table", &["regulatory"]),
        ("motif_table", &["motif"]),
    ];

    let mut entries = Vec::new();
    for (json_key, file_stems) in &table_specs {
        for file_stem in *file_stems {
            let parquet_name = format!("{base}_{file_stem}{suffix}.parquet");
            let parquet_path = context_dir.join(&parquet_name);
            if parquet_path.exists() {
                let table_name = format!("{base}_{file_stem}{suffix}");
                entries.push(format!("\"{json_key}\":\"{table_name}\""));
                break;
            }
        }
    }

    // Discover split translation_sift for sift/polyphen window loading.
    {
        let sift_name = format!("{base}_translation_sift{suffix}.parquet");
        let sift_path = context_dir.join(&sift_name);
        if sift_path.exists() {
            let table_name = format!("{base}_translation_sift{suffix}");
            entries.push(format!("\"translations_sift_table\":\"{table_name}\""));
        }
    }

    // Extended probes: use interval-overlap fallback for shifted indels.
    entries.push(format!("\"extended_probes\":{}", args.extended_probes));

    if args.everything {
        // --everything implies --hgvs and all other flags.
        entries.push("\"everything\":true".to_string());
        // --everything implies HGVS, so we need the reference FASTA.
        let fasta_path = args.reference_fasta_path.as_ref().ok_or_else(|| {
            DataFusionError::Execution(
                "--everything requires --reference-fasta-path=/path/to/reference.fa[.gz]"
                    .to_string(),
            )
        })?;
        entries.push(format!(
            "\"reference_fasta_path\":\"{}\"",
            sql_literal(fasta_path.to_str().ok_or_else(|| {
                DataFusionError::Execution("reference_fasta_path must be valid UTF-8".to_string())
            })?)
        ));
    } else {
        if args.hgvs {
            entries.push("\"hgvs\":true".to_string());
            if let Some(shift_hgvs) = args.shift_hgvs {
                entries.push(format!("\"shift_hgvs\":{shift_hgvs}"));
            }
            let fasta_path = args.reference_fasta_path.as_ref().ok_or_else(|| {
                DataFusionError::Execution(
                    "--hgvs requires --reference-fasta-path=/path/to/reference.fa[.gz]".to_string(),
                )
            })?;
            entries.push(format!(
                "\"reference_fasta_path\":\"{}\"",
                sql_literal(fasta_path.to_str().ok_or_else(|| {
                    DataFusionError::Execution(
                        "reference_fasta_path must be valid UTF-8".to_string(),
                    )
                })?)
            ));
        }

        // Batch 3 flags — enable frequency and clinical CSQ fields.
        entries.push("\"check_existing\":true".to_string());
        entries.push("\"af\":true".to_string());
        entries.push("\"af_1kg\":true".to_string());
        entries.push("\"af_gnomade\":true".to_string());
        entries.push("\"af_gnomadg\":true".to_string());
        entries.push("\"max_af\":true".to_string());
        entries.push("\"pubmed\":true".to_string());
    }

    if args.refseq {
        entries.push("\"refseq\":true".to_string());
    }

    if args.merged {
        entries.push("\"merged\":true".to_string());
    }

    if args.gencode_basic {
        entries.push("\"gencode_basic\":true".to_string());
    }

    if args.gencode_primary {
        entries.push("\"gencode_primary\":true".to_string());
    }

    if args.all_refseq {
        entries.push("\"all_refseq\":true".to_string());
    }

    if args.exclude_predicted {
        entries.push("\"exclude_predicted\":true".to_string());
    }

    if entries.is_empty() {
        return Ok(None);
    }

    Ok(Some(format!("{{{}}}", entries.join(","))))
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
        "translations_sift_table",
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
    refseq: bool,
    merged: bool,
    gencode_basic: bool,
    gencode_primary: bool,
    all_refseq: bool,
    exclude_predicted: bool,
    hgvs: bool,
    shift_hgvs: Option<bool>,
    reference_fasta_path: Option<&Path>,
    everything: bool,
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
        .arg(volume_work);

    if hgvs {
        let fasta_path = reference_fasta_path.ok_or_else(|| {
            DataFusionError::Execution(
                "--hgvs requires --reference-fasta-path=/path/to/reference.fa[.gz]".to_string(),
            )
        })?;
        let fasta_parent = fasta_path.parent().ok_or_else(|| {
            DataFusionError::Execution(format!(
                "reference FASTA has no parent directory: {}",
                fasta_path.display()
            ))
        })?;
        let fasta_name = fasta_path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "invalid reference FASTA file name: {}",
                    fasta_path.display()
                ))
            })?;
        cmd.arg("-v")
            .arg(format!("{}:/fasta:ro", fasta_parent.display()));
        cmd.arg("ensemblorg/ensembl-vep:release_115.2")
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
            .arg("--no_stats")
            .arg("--hgvs")
            .arg("--fasta")
            .arg(format!("/fasta/{fasta_name}"));
        if let Some(value) = shift_hgvs {
            cmd.arg("--shift_hgvs").arg(if value { "1" } else { "0" });
        }
    } else {
        cmd.arg("ensemblorg/ensembl-vep:release_115.2")
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
    }

    if everything {
        // --everything enables all VEP features; VEP determines its own 80-field order.
        cmd.arg("--everything");
    } else {
        // Enable regulatory feature annotations (VEP skips them by default).
        cmd.arg("--regulatory");

        // Batch 1 VEP flags — enable additional CSQ fields for comparison.
        cmd.arg("--variant_class");
        cmd.arg("--canonical");
        cmd.arg("--tsl");
        cmd.arg("--mane");
        cmd.arg("--protein");
        cmd.arg("--gene_phenotype");
        cmd.arg("--ccds");
        cmd.arg("--uniprot");

        // Batch 3 VEP flags — enable frequency and clinical CSQ fields.
        cmd.arg("--check_existing");
        cmd.arg("--af");
        cmd.arg("--af_1kg");
        cmd.arg("--af_gnomade");
        cmd.arg("--af_gnomadg");
        cmd.arg("--max_af");
        cmd.arg("--pubmed");

        let fields = csq_field_names_for_mode(false, refseq, merged).join(",");
        cmd.arg("--fields").arg(fields);
    }

    if refseq {
        cmd.arg("--refseq");
        cmd.arg("--use_given_ref");
    }
    if merged {
        cmd.arg("--merged");
        // Merged caches auto-enable --use_transcript_ref which requires FASTA.
        // Use --use_given_ref to skip the FASTA requirement.
        cmd.arg("--use_given_ref");
    }
    if gencode_basic {
        cmd.arg("--gencode_basic");
    }
    if gencode_primary {
        cmd.arg("--gencode_primary");
    }
    if all_refseq {
        cmd.arg("--all_refseq");
    }
    if exclude_predicted {
        cmd.arg("--exclude_predicted");
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
            "SELECT chrom, start, ref, alt, \"CSQ\", most_severe_consequence FROM annotate_vep('sampled_vcf', '{}', '{}', '{}')",
            sql_literal(cache_source),
            sql_literal(backend),
            sql_literal(opts)
        )
    } else {
        format!(
            "SELECT chrom, start, ref, alt, \"CSQ\", most_severe_consequence FROM annotate_vep('sampled_vcf', '{}', '{}')",
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
        let csq_col = batch.column_by_name("CSQ").ok_or_else(|| {
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
        // Print mismatch samples for fields with mismatches
        if matches < report.total_entries_compared {
            if let Some(samples) = report.field_mismatch_samples.get(&i) {
                for s in samples {
                    println!(
                        "    mismatch: {}:{} {}|{} golden={:?} ours={:?}",
                        s.variant_key.chrom,
                        s.variant_key.pos,
                        s.allele,
                        s.feature,
                        s.golden_val,
                        s.ours_val
                    );
                }
            }
        }
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
            let label = if ft.is_empty() {
                "(empty)"
            } else {
                ft.as_str()
            };
            println!("    {:<25} {}", label, count);
        }
    }
    if !report.ours_only_by_ft.is_empty() {
        println!("  ours_only by Feature_type:");
        let mut sorted: Vec<_> = report.ours_only_by_ft.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));
        for (ft, count) in &sorted {
            let label = if ft.is_empty() {
                "(empty)"
            } else {
                ft.as_str()
            };
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

fn print_multiplicity_report(report: &CsqMultiplicityReport) {
    if report.golden_extra_total == 0 && report.ours_extra_total == 0 {
        return;
    }

    println!("\ncsq multiplicity diagnostic:");
    println!("  golden_extra_duplicates: {}", report.golden_extra_total);
    println!("  ours_extra_duplicates: {}", report.ours_extra_total);
    if !report.samples.is_empty() {
        println!(
            "  samples (chrom:pos ref/alt allele|feature|feature_type golden_count->ours_count):"
        );
        for sample in &report.samples {
            println!(
                "    {}:{} {}/{} {}|{}|{} {}->{}",
                sample.variant_key.chrom,
                sample.variant_key.pos,
                sample.variant_key.ref_allele,
                sample.variant_key.alt_alleles,
                sample.allele,
                sample.feature,
                if sample.feature_type.is_empty() {
                    "(empty)"
                } else {
                    sample.feature_type.as_str()
                },
                sample.golden_count,
                sample.ours_count,
            );
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

fn print_summary_table(
    report: &ComparisonReport,
    term_report: &TermComparisonReport,
    csq: &CsqFieldReport,
    elapsed_s: f64,
) {
    println!("\n{}", "=".repeat(80));
    println!("SUMMARY");
    println!("{}", "=".repeat(80));

    // Variant-level stats
    println!(
        "\nVariants: {} golden, {} ours, {} matched, {} missing, {} extra",
        report.golden_rows,
        report.ours_rows,
        report.intersection_rows,
        report.missing_in_ours,
        report.extra_in_ours,
    );
    if report.intersection_rows > 0 {
        let ms_pct =
            report.most_severe_exact_matches as f64 / report.intersection_rows as f64 * 100.0;
        let ts_pct = term_report.term_set_exact_matches as f64
            / term_report.comparable_rows.max(1) as f64
            * 100.0;
        println!(
            "Accuracy: most_severe={:.2}%, term_set={:.2}%",
            ms_pct, ts_pct,
        );
    }
    println!("Time: {:.1}s", elapsed_s);

    // Per-field table
    if csq.total_entries_compared > 0 {
        let total = csq.total_entries_compared;
        println!(
            "\nCSQ per-field accuracy ({} entries across {} fields):",
            total,
            csq.field_names.len(),
        );
        println!(
            "{:<25} {:>12} {:>12} {:>8}",
            "Field", "Matched", "Mismatches", "Accuracy"
        );
        println!("{}", "-".repeat(60));

        let mut perfect_count = 0usize;
        let mut imperfect: Vec<(&str, usize, f64)> = Vec::new();

        for (i, name) in csq.field_names.iter().enumerate() {
            let matched = csq.field_match_counts[i];
            let mismatches = total - matched;
            let pct = matched as f64 / total as f64 * 100.0;
            let status = if mismatches == 0 { " " } else { "*" };
            println!(
                "{}{:<24} {:>12} {:>12} {:>7.2}%",
                status, name, matched, mismatches, pct,
            );
            if mismatches == 0 {
                perfect_count += 1;
            } else {
                imperfect.push((name, mismatches, pct));
            }
        }

        println!("{}", "-".repeat(60));
        println!(
            "Perfect (0 mismatches): {}/{} fields",
            perfect_count,
            csq.field_names.len(),
        );
        if !imperfect.is_empty() {
            imperfect.sort_by(|a, b| b.1.cmp(&a.1));
            println!(
                "Imperfect: {}/{} fields",
                imperfect.len(),
                csq.field_names.len(),
            );
            for (name, mm, pct) in &imperfect {
                println!("  {:<24} {:>8} mismatches ({:.2}% accuracy)", name, mm, pct);
            }
        }
    }

    println!("{}", "=".repeat(80));
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
