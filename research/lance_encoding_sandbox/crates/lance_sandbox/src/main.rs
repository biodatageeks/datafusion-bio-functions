mod bench;
mod build;
mod config;
mod inspect;
mod keys;
mod payload_sidecar;
mod report;
mod row_sidecar;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use arrow_array::{Array, Int8Array, Int32Array, Int64Array, UInt8Array, UInt32Array, UInt64Array};
use arrow_schema::DataType;
use clap::{Parser, Subcommand};
use flate2::read::MultiGzDecoder;
use futures::TryStreamExt;
use heed::{Database, EnvFlags, EnvOpenOptions};
use lance::Dataset;
use lance::dataset::scanner::MaterializationStyle;

use crate::config::SandboxConfig;

#[derive(Debug, Parser)]
#[command(name = "lance-sandbox")]
#[command(about = "Research-only Lance encoding/layout sandbox")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Build {
        #[arg(long)]
        config: PathBuf,
    },
    Inspect {
        #[arg(long)]
        config: PathBuf,
    },
    InspectPath {
        #[arg(long)]
        dataset_path: PathBuf,
        #[arg(long, default_value = "unknown")]
        lance_version: String,
        #[arg(long = "input-parquet")]
        input_parquet: Vec<PathBuf>,
        #[arg(long, default_value = "position")]
        position_field: String,
        #[arg(long, default_value = "start")]
        position_source_column: String,
        #[arg(long, default_value_t = true)]
        include_physical_columns: bool,
        #[arg(long, default_value_t = true)]
        include_index_sizes: bool,
    },
    Bench {
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        positions_file: Option<PathBuf>,
        #[arg(long)]
        dataset_path: Option<PathBuf>,
        #[arg(long, default_value = "1")]
        lookup_tier: String,
        #[arg(long, value_delimiter = ',')]
        lookup_tiers: Vec<String>,
        #[arg(long)]
        only_index_take_everything: bool,
    },
    Run {
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        positions_file: Option<PathBuf>,
        #[arg(long, default_value = "1")]
        lookup_tier: String,
    },
    ExtractKeys {
        #[arg(long)]
        config: PathBuf,
    },
    SampleCold {
        #[arg(long)]
        config: PathBuf,
        #[arg(long, value_delimiter = ',', default_value = "10000,100000")]
        limits: Vec<usize>,
        #[arg(long)]
        output_dir: Option<PathBuf>,
    },
    /// Read a parquet file and write it as a single-fragment Lance dataset, stamping
    /// `lance-encoding:structural-encoding=fullzip` metadata onto the named columns.
    /// Used to (re)write the bundled-AF candidate via the patched/fixed lance writer.
    WriteLanceBundle {
        #[arg(long)]
        input_parquet: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "2.2")]
        version: String,
        #[arg(long, value_delimiter = ',')]
        fullzip_cols: Vec<String>,
        /// "overwrite" (default) or "append" — append builds a multi-fragment dataset
        /// one source fragment at a time (preserving fragment order/addresses).
        #[arg(long, default_value = "overwrite")]
        mode: String,
        /// "none" (default, lets Lance pick e.g. FSST) or "zstd" — stamps
        /// `lance-encoding:compression` on the fullzip columns to match a zstd baseline.
        #[arg(long, default_value = "none")]
        compression: String,
    },
    /// Time `dataset.take_rows(row_ids, projection)` for an explicit row-address
    /// set and column projection. Used for the AF-bundling encoding A/B
    /// (separate-miniblock columns vs bundled-fullzip lists).
    TakeColsBench {
        #[arg(long)]
        dataset_path: PathBuf,
        /// File of u64 Lance row addresses (fragment<<32 | offset), one per line.
        #[arg(long)]
        rowids_file: PathBuf,
        /// Comma-separated projection columns.
        #[arg(long, value_delimiter = ',')]
        columns: Vec<String>,
        #[arg(long, default_value_t = 5)]
        iterations: usize,
        /// Materialize each taken column to bytes (force full decode).
        #[arg(long, default_value_t = true)]
        materialize: bool,
    },
    /// Replay a batch-delimited row-address file (blank line between batches) as N separate
    /// take_rows calls — the realistic per-batch lookup pattern. Sums disk bytes / IOPS / time.
    TakeBatchesBench {
        #[arg(long)]
        dataset_path: PathBuf,
        #[arg(long)]
        batches_file: PathBuf,
        #[arg(long, value_delimiter = ',')]
        columns: Vec<String>,
        #[arg(long, default_value_t = 2)]
        iterations: usize,
    },
    PayloadPeek {
        #[arg(long)]
        db_path: PathBuf,
        #[arg(long)]
        positions_file: PathBuf,
        #[arg(long, default_value_t = 5)]
        limit: usize,
    },
    PayloadBench {
        #[arg(long)]
        db_path: PathBuf,
        #[arg(long)]
        positions_file: PathBuf,
        #[arg(long)]
        materialize: bool,
        #[arg(long)]
        no_read_ahead: bool,
        #[arg(long)]
        raw_only: bool,
    },
    PayloadExportKeys {
        #[arg(long)]
        db_path: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 100_000)]
        limit: usize,
        #[arg(long, default_value_t = 0)]
        skip: usize,
    },
    VcfPositionMatch {
        #[arg(long)]
        dataset_path: PathBuf,
        #[arg(long)]
        vcf_path: PathBuf,
        #[arg(long, default_value = "chr1")]
        chrom: String,
        #[arg(long, default_value = "position")]
        position_column: String,
        #[arg(long, default_value = "tier")]
        tier_column: String,
        #[arg(long)]
        write_all_positions: Option<PathBuf>,
        #[arg(long)]
        write_matched_cold_positions: Option<PathBuf>,
    },
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Build { config } => {
            let cfg = SandboxConfig::read(&config)?;
            let build = build::build_dataset(&cfg).await?;
            report::write_build(&cfg, &build)?;
        }
        Command::Inspect { config } => {
            let cfg = SandboxConfig::read(&config)?;
            let inspected = inspect::inspect_dataset(&cfg).await?;
            let inspect_path = report::write_inspect(&cfg, &inspected)?;
            report::print_terminal_markdown(&report::render_inspect_stdout(&inspected));
            report::print_terminal_markdown(&format!("Wrote `{}`\n", inspect_path.display()));
        }
        Command::InspectPath {
            dataset_path,
            lance_version,
            input_parquet,
            position_field,
            position_source_column,
            include_physical_columns,
            include_index_sizes,
        } => {
            let inspected = inspect::inspect_dataset_path(
                &dataset_path,
                lance_version,
                &input_parquet,
                Some(position_field.as_str()),
                Some(position_source_column.as_str()),
                include_physical_columns,
                include_index_sizes,
            )
            .await?;
            println!("{}", serde_json::to_string_pretty(&inspected)?);
        }
        Command::Bench {
            config,
            positions_file,
            dataset_path,
            lookup_tier,
            lookup_tiers,
            only_index_take_everything,
        } => {
            let cfg = SandboxConfig::read(&config)?;
            let positions_file = resolve_positions_file(&cfg, &config, positions_file).await?;
            let lookup_tiers = resolve_lookup_tiers(&lookup_tier, lookup_tiers)?;
            let multi_tier = lookup_tiers.len() > 1;
            let benchmark_mode = if only_index_take_everything {
                bench::BenchmarkMode::OnlyIndexTakeEverything
            } else {
                bench::BenchmarkMode::Full
            };
            for lookup_tier in lookup_tiers {
                eprintln!(
                    "Running benchmark for {}",
                    lookup_scope_display(lookup_tier)
                );
                let result = bench::run_benchmark_with_mode(
                    &cfg,
                    &positions_file,
                    dataset_path.as_deref(),
                    lookup_tier,
                    benchmark_mode,
                    |partial| {
                        let path = write_benchmark_progress(
                            &cfg,
                            dataset_path.as_deref(),
                            partial,
                            multi_tier,
                        )?;
                        eprintln!("Updated partial report `{}`", path.display());
                        Ok(())
                    },
                )
                .await?;
                let path =
                    write_benchmark_report(&cfg, dataset_path.as_deref(), &result, multi_tier)?;
                report::print_terminal_markdown(&report::render_benchmark_stdout(&result));
                report::print_terminal_markdown(&format!("Wrote `{}`\n", path.display()));
            }
        }
        Command::Run {
            config,
            positions_file,
            lookup_tier,
        } => {
            let cfg = SandboxConfig::read(&config)?;
            let build = build::build_dataset(&cfg).await?;
            report::write_build(&cfg, &build)?;
            let inspected = inspect::inspect_dataset(&cfg).await?;
            report::write_inspect(&cfg, &inspected)?;
            let positions_file = resolve_positions_file(&cfg, &config, positions_file).await?;
            let benchmark = bench::run_benchmark_with_progress(
                &cfg,
                &positions_file,
                None,
                parse_lookup_tier(&lookup_tier)?,
                |_| Ok(()),
            )
            .await?;
            report::write_benchmark(&cfg, &benchmark)?;
            report::write_summary(&cfg, &inspected, &benchmark)?;
        }
        Command::ExtractKeys { config } => {
            let cfg = SandboxConfig::read(&config)?;
            let positions_file = cfg.resolved_positions_file(&config);
            keys::extract_keys(&cfg, &positions_file).await?;
        }
        Command::SampleCold {
            config,
            limits,
            output_dir,
        } => {
            let cfg = SandboxConfig::read(&config)?;
            let output_dir = output_dir.unwrap_or_else(|| {
                config
                    .parent()
                    .and_then(|p| p.parent())
                    .unwrap_or_else(|| std::path::Path::new("."))
                    .join("inputs")
            });
            let report = keys::sample_cold_positions(&cfg, &output_dir, &limits).await?;
            println!("Unique cold positions: {}", report.unique_cold_positions);
            for file in report.files {
                println!(
                    "Wrote {} rows for limit {} to `{}`",
                    file.rows,
                    file.limit,
                    file.path.display()
                );
            }
        }
        Command::TakeColsBench {
            dataset_path,
            rowids_file,
            columns,
            iterations,
            materialize,
        } => {
            take_cols_bench(&dataset_path, &rowids_file, &columns, iterations, materialize).await?;
        }
        Command::TakeBatchesBench {
            dataset_path,
            batches_file,
            columns,
            iterations,
        } => {
            take_batches_bench(&dataset_path, &batches_file, &columns, iterations).await?;
        }
        Command::WriteLanceBundle {
            input_parquet,
            out,
            version,
            fullzip_cols,
            mode,
            compression,
        } => {
            write_lance_bundle(&input_parquet, &out, &version, &fullzip_cols, &mode, &compression)
                .await?;
        }
        Command::PayloadPeek {
            db_path,
            positions_file,
            limit,
        } => {
            payload_peek(&db_path, &positions_file, limit)?;
        }
        Command::PayloadBench {
            db_path,
            positions_file,
            materialize,
            no_read_ahead,
            raw_only,
        } => {
            payload_bench(
                &db_path,
                &positions_file,
                materialize,
                no_read_ahead,
                raw_only,
            )?;
        }
        Command::PayloadExportKeys {
            db_path,
            output,
            limit,
            skip,
        } => {
            payload_export_keys(&db_path, &output, limit, skip)?;
        }
        Command::VcfPositionMatch {
            dataset_path,
            vcf_path,
            chrom,
            position_column,
            tier_column,
            write_all_positions,
            write_matched_cold_positions,
        } => {
            let report = vcf_position_match(
                &dataset_path,
                &vcf_path,
                &chrom,
                &position_column,
                &tier_column,
            )
            .await?;
            if let Some(path) = write_all_positions {
                write_positions_file(&path, &report.vcf_positions)?;
                eprintln!(
                    "Wrote {} all VCF positions to `{}`",
                    report.vcf_positions.len(),
                    path.display()
                );
            }
            if let Some(path) = write_matched_cold_positions {
                write_positions_file(&path, &report.matched_cold_positions)?;
                eprintln!(
                    "Wrote {} cold-matched VCF positions to `{}`",
                    report.matched_cold_positions.len(),
                    path.display()
                );
            }
            report::print_terminal_markdown(&render_vcf_position_match(&report));
        }
    }
    Ok(())
}

async fn resolve_positions_file(
    cfg: &SandboxConfig,
    config_path: &Path,
    positions_file: Option<PathBuf>,
) -> Result<PathBuf> {
    let using_default_positions = positions_file.is_none();
    let positions_file = positions_file.unwrap_or_else(|| cfg.resolved_positions_file(config_path));
    if !positions_file.exists() {
        if using_default_positions {
            eprintln!(
                "Positions file '{}' is missing; extracting keys first.",
                positions_file.display()
            );
            keys::extract_keys(cfg, &positions_file).await?;
        } else {
            bail!(
                "positions file '{}' does not exist; run sample-cold or extract-keys first",
                positions_file.display()
            );
        }
    }
    Ok(positions_file)
}

fn resolve_lookup_tiers(lookup_tier: &str, lookup_tiers: Vec<String>) -> Result<Vec<Option<u8>>> {
    let mut lookup_tiers = lookup_tiers
        .iter()
        .map(|tier| parse_lookup_tier(tier))
        .collect::<Result<Vec<_>>>()?;
    if lookup_tiers.is_empty() {
        lookup_tiers.push(parse_lookup_tier(lookup_tier)?);
    }
    lookup_tiers.sort_unstable_by_key(|tier| match tier {
        Some(tier) => *tier,
        None => u8::MAX,
    });
    lookup_tiers.dedup();
    Ok(lookup_tiers)
}

fn parse_lookup_tier(raw: &str) -> Result<Option<u8>> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "0" => Ok(Some(0)),
        "1" => Ok(Some(1)),
        "all" | "joint" | "both" | "0+1" | "1+0" => Ok(None),
        other => bail!("lookup tier must be 0, 1, or all; got '{other}'"),
    }
}

fn lookup_scope_display(lookup_tier: Option<u8>) -> String {
    match lookup_tier {
        Some(tier) => format!("lookup_tier={tier}"),
        None => "all tiers without tier filtering".to_string(),
    }
}

fn write_benchmark_report(
    cfg: &SandboxConfig,
    dataset_path: Option<&Path>,
    benchmark: &bench::BenchmarkReport,
    tiered_output: bool,
) -> Result<PathBuf> {
    match (dataset_path, tiered_output) {
        (Some(dataset_path), true) => {
            report::write_benchmark_for_dataset_path_tier(dataset_path, benchmark)
        }
        (Some(dataset_path), false) => {
            report::write_benchmark_for_dataset_path(dataset_path, benchmark)
        }
        (None, true) => report::write_benchmark_for_tier(cfg, benchmark),
        (None, false) => report::write_benchmark(cfg, benchmark),
    }
}

async fn write_lance_bundle(
    input_parquet: &Path,
    out: &Path,
    version: &str,
    fullzip_cols: &[String],
    mode: &str,
    compression: &str,
) -> Result<()> {
    use arrow_schema::{Field, Schema};
    use lance::dataset::{WriteMode, WriteParams};
    use lance_file::version::LanceFileVersion;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use std::fs::File;
    use std::sync::Arc;

    let data_storage_version = match version {
        "2.1" => LanceFileVersion::V2_1,
        "2.2" => LanceFileVersion::V2_2,
        other => anyhow::bail!("unsupported lance version '{other}'"),
    };
    let write_mode = match mode {
        "overwrite" => WriteMode::Overwrite,
        "append" => WriteMode::Append,
        other => anyhow::bail!("unsupported write mode '{other}'"),
    };

    let file = File::open(input_parquet)
        .with_context(|| format!("failed to open '{}'", input_parquet.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)?.build()?;
    let batches = reader.collect::<std::result::Result<Vec<_>, _>>()?;
    anyhow::ensure!(!batches.is_empty(), "parquet file produced no batches");
    let orig_schema = batches[0].schema();

    // Stamp fullzip structural-encoding metadata onto the requested columns.
    let fullzip: std::collections::HashSet<&str> =
        fullzip_cols.iter().map(String::as_str).collect();
    let new_fields = orig_schema
        .fields()
        .iter()
        .map(|f| {
            if fullzip.contains(f.name().as_str()) {
                let mut md = f.metadata().clone();
                md.insert(
                    "lance-encoding:structural-encoding".to_string(),
                    "fullzip".to_string(),
                );
                if compression == "zstd" {
                    md.insert("lance-encoding:compression".to_string(), "zstd".to_string());
                    md.insert(
                        "lance-encoding:compression-level".to_string(),
                        "3".to_string(),
                    );
                }
                Arc::new(Field::clone(f).with_metadata(md))
            } else {
                f.clone()
            }
        })
        .collect::<Vec<_>>();
    let schema = Arc::new(Schema::new_with_metadata(
        new_fields,
        orig_schema.metadata().clone(),
    ));

    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    let restamped = batches
        .into_iter()
        .map(|b| arrow_array::RecordBatch::try_new(schema.clone(), b.columns().to_vec()))
        .collect::<std::result::Result<Vec<_>, _>>()?;

    eprintln!(
        "WriteLanceBundle: rows={total_rows} version={version} fullzip_cols={fullzip_cols:?} -> {}",
        out.display()
    );

    let reader = arrow_array::RecordBatchIterator::new(restamped.into_iter().map(Ok), schema);
    let params = WriteParams {
        mode: write_mode,
        data_storage_version: Some(data_storage_version),
        max_rows_per_file: total_rows + 1,
        ..Default::default()
    };
    lance::Dataset::write(reader, out.to_string_lossy().as_ref(), Some(params))
        .await
        .with_context(|| format!("failed to write '{}'", out.display()))?;
    eprintln!("wrote {} ({total_rows} rows)", out.display());
    Ok(())
}

async fn take_cols_bench(
    dataset_path: &Path,
    rowids_file: &Path,
    columns: &[String],
    iterations: usize,
    materialize: bool,
) -> Result<()> {
    use arrow_array::Array;
    use lance::dataset::ProjectionRequest;

    let row_ids = read_positions_file(rowids_file)?;
    let dataset = Dataset::open(
        dataset_path
            .to_str()
            .context("dataset path is not valid UTF-8")?,
    )
    .await?;

    eprintln!(
        "TakeColsBench: dataset={}, row_ids={}, columns={:?}, iterations={}",
        dataset_path.display(),
        row_ids.len(),
        columns,
        iterations
    );

    let mut wall_secs = Vec::with_capacity(iterations);
    let mut cpu_secs = Vec::with_capacity(iterations);
    let mut decoded_bytes = 0usize;
    let mut taken_rows = 0usize;

    for iter in 0..iterations {
        let projection =
            ProjectionRequest::from_columns(columns.iter().map(String::as_str), dataset.schema());
        let io_bytes_before = lance_io::bytes_read_counter();
        let io_iops_before = lance_io::iops_counter();
        let usage_before = ResourceUsage::snapshot();
        let started = Instant::now();
        let batch = dataset.take_rows(&row_ids, projection).await?;
        // force full materialization of every projected column
        let mut bytes = 0usize;
        if materialize {
            for col in batch.columns() {
                bytes += col.get_array_memory_size();
            }
        }
        let elapsed = started.elapsed().as_secs_f64();
        let usage = ResourceUsage::snapshot().delta(usage_before);
        let cpu = usage.user_seconds + usage.system_seconds;
        let disk_bytes = lance_io::bytes_read_counter().saturating_sub(io_bytes_before);
        let iops = lance_io::iops_counter().saturating_sub(io_iops_before);
        taken_rows = batch.num_rows();
        decoded_bytes = bytes;
        wall_secs.push(elapsed);
        cpu_secs.push(cpu);
        eprintln!(
            "  iter {}: wall={:.4}s cpu={:.4}s rows={} decoded={} disk_bytes_read={} iops={} major_faults={}",
            iter, elapsed, cpu, batch.num_rows(), bytes, disk_bytes, iops, usage.major_faults
        );
    }

    let best_wall = wall_secs.iter().cloned().fold(f64::INFINITY, f64::min);
    let best_cpu = cpu_secs.iter().cloned().fold(f64::INFINITY, f64::min);
    let med = |v: &mut Vec<f64>| {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    let med_wall = med(&mut wall_secs);
    let med_cpu = med(&mut cpu_secs);
    println!(
        "RESULT dataset={} columns={} rows={} decoded_bytes={} best_wall_s={:.4} med_wall_s={:.4} best_cpu_s={:.4} med_cpu_s={:.4} us_per_row_cpu={:.3}",
        dataset_path.display(),
        columns.len(),
        taken_rows,
        decoded_bytes,
        best_wall,
        med_wall,
        best_cpu,
        med_cpu,
        best_cpu * 1_000_000.0 / taken_rows.max(1) as f64
    );
    Ok(())
}

/// Parse a batch-delimited address file: blank line separates batches.
fn read_batches_file(path: &Path) -> Result<Vec<Vec<u64>>> {
    let raw = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read batches file '{}'", path.display()))?;
    let mut batches = Vec::new();
    let mut cur = Vec::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() {
            if !cur.is_empty() {
                batches.push(std::mem::take(&mut cur));
            }
        } else if !t.starts_with('#') {
            cur.push(t.parse::<u64>()?);
        }
    }
    if !cur.is_empty() {
        batches.push(cur);
    }
    Ok(batches)
}

/// Replay the per-batch take pattern: one take_rows() per batch, summing disk bytes/IOPS/time.
async fn take_batches_bench(
    dataset_path: &Path,
    batches_file: &Path,
    columns: &[String],
    iterations: usize,
) -> Result<()> {
    use lance::dataset::ProjectionRequest;

    let batches = read_batches_file(batches_file)?;
    let dataset = Dataset::open(dataset_path.to_str().context("path not utf-8")?).await?;
    let total_ids: usize = batches.iter().map(Vec::len).sum();
    eprintln!(
        "TakeBatchesBench: dataset={}, batches={}, total_ids={}, columns={:?}",
        dataset_path.display(),
        batches.len(),
        total_ids,
        columns
    );

    for iter in 0..iterations {
        let io_b0 = lance_io::bytes_read_counter();
        let io_i0 = lance_io::iops_counter();
        let usage0 = ResourceUsage::snapshot();
        let started = Instant::now();
        let mut rows = 0usize;
        let mut decoded = 0usize;
        for batch_ids in &batches {
            let projection =
                ProjectionRequest::from_columns(columns.iter().map(String::as_str), dataset.schema());
            let b = dataset.take_rows(batch_ids, projection).await?;
            rows += b.num_rows();
            for col in b.columns() {
                decoded += col.get_array_memory_size();
            }
        }
        let wall = started.elapsed().as_secs_f64();
        let usage = ResourceUsage::snapshot().delta(usage0);
        let cpu = usage.user_seconds + usage.system_seconds;
        let disk_bytes = lance_io::bytes_read_counter().saturating_sub(io_b0);
        let iops = lance_io::iops_counter().saturating_sub(io_i0);
        eprintln!(
            "  iter {iter}: batches={} rows={rows} decoded={decoded} wall={wall:.4}s cpu={cpu:.4}s disk_bytes_read={disk_bytes} iops={iops} major_faults={}",
            batches.len(),
            usage.major_faults
        );
    }
    Ok(())
}

fn write_benchmark_progress(
    cfg: &SandboxConfig,
    dataset_path: Option<&Path>,
    benchmark: &bench::BenchmarkReport,
    tiered_output: bool,
) -> Result<PathBuf> {
    match (dataset_path, tiered_output) {
        (Some(dataset_path), true) => {
            report::write_benchmark_progress_for_dataset_path_tier(dataset_path, benchmark)
        }
        (Some(dataset_path), false) => {
            report::write_benchmark_progress_for_dataset_path(dataset_path, benchmark)
        }
        (None, true) => report::write_benchmark_progress_for_tier(cfg, benchmark),
        (None, false) => report::write_benchmark_progress(cfg, benchmark),
    }
}

fn payload_peek(db_path: &Path, positions_file: &Path, limit: usize) -> Result<()> {
    let positions = read_positions_file(positions_file)?;
    let env = unsafe {
        EnvOpenOptions::new()
            .map_size(64 * 1024 * 1024 * 1024)
            .max_dbs(4)
            .open(db_path)
            .with_context(|| format!("failed to open heed env '{}'", db_path.display()))?
    };
    let rtxn = env.read_txn()?;
    let db: Database<payload_sidecar::U32KeyCodec, heed::types::Bytes> = env
        .open_database(&rtxn, Some("payloads"))?
        .context("heed payload database is missing")?;
    let typed_db: Database<
        payload_sidecar::U32KeyCodec,
        payload_sidecar::ZstdEverythingPayloadCodec,
    > = env
        .open_database(&rtxn, Some("payloads"))?
        .context("heed payload database is missing")?;

    for position in positions.iter().take(limit) {
        let position = u32::try_from(*position)
            .with_context(|| format!("position {position} does not fit u32"))?;
        match db.get(&rtxn, &position) {
            Ok(Some(bytes)) => {
                let prefix = &bytes[..bytes.len().min(8)];
                match payload_sidecar::decode_payload_zstd(bytes) {
                    Ok(payload) => {
                        let typed_payload = typed_db.get(&rtxn, &position)?;
                        println!(
                            "position={position} bytes={} prefix={prefix:02x?} rows={} typed_rows={}",
                            bytes.len(),
                            payload.rows.len(),
                            typed_payload
                                .as_ref()
                                .map(|payload| payload.rows.len())
                                .unwrap_or_default()
                        );
                    }
                    Err(error) => {
                        println!(
                            "position={position} bytes={} prefix={prefix:02x?} decode_error={error:#}",
                            bytes.len()
                        );
                    }
                }
            }
            Ok(None) => println!("position={position} missing"),
            Err(error) => println!("position={position} get_error={error:#}"),
        }
    }

    Ok(())
}

fn payload_bench(
    db_path: &Path,
    positions_file: &Path,
    materialize: bool,
    no_read_ahead: bool,
    raw_only: bool,
) -> Result<()> {
    let positions = read_positions_file(positions_file)?;
    let requested_positions = positions.len();
    let env = unsafe {
        let mut options = EnvOpenOptions::new();
        if no_read_ahead {
            options.flags(EnvFlags::NO_READ_AHEAD);
        }
        options
            .map_size(64 * 1024 * 1024 * 1024)
            .max_dbs(4)
            .open(db_path)
            .with_context(|| format!("failed to open heed env '{}'", db_path.display()))?
    };
    let rtxn = env.read_txn()?;

    let usage_before = ResourceUsage::snapshot();
    let started = Instant::now();
    let mut matched_positions = 0usize;
    let mut row_count = 0usize;
    let mut rows = Vec::new();
    let mut raw_bytes = 0usize;

    if raw_only {
        let db: Database<payload_sidecar::U32KeyCodec, heed::types::Bytes> = env
            .open_database(&rtxn, Some("payloads"))?
            .context("heed payload database is missing")?;
        for position in positions {
            let Ok(position) = u32::try_from(position) else {
                continue;
            };
            if let Some(bytes) = db.get(&rtxn, &position)? {
                matched_positions += 1;
                raw_bytes += bytes.len();
            }
        }
    } else {
        let db: Database<
            payload_sidecar::U32KeyCodec,
            payload_sidecar::ZstdEverythingPayloadCodec,
        > = env
            .open_database(&rtxn, Some("payloads"))?
            .context("heed payload database is missing")?;
        for position in positions {
            let Ok(position) = u32::try_from(position) else {
                continue;
            };
            if let Some(payload) = db.get(&rtxn, &position)? {
                matched_positions += 1;
                row_count += payload.rows.len();
                if materialize {
                    rows.extend(payload.rows);
                }
            }
        }
    }

    let elapsed = started.elapsed().as_secs_f64();
    let usage_after = ResourceUsage::snapshot();
    let usage_delta = usage_after.delta(usage_before);
    println!(
        "mode={} positions={} matched={} rows={} materialized={} raw_bytes={} seconds={:.6} minor_faults={} major_faults={} user_seconds={:.6} system_seconds={:.6}",
        if raw_only {
            "raw_bytes"
        } else if materialize {
            "zstd_materialize"
        } else {
            "zstd_count"
        },
        requested_positions,
        matched_positions,
        row_count,
        rows.len(),
        raw_bytes,
        elapsed,
        usage_delta.minor_faults,
        usage_delta.major_faults,
        usage_delta.user_seconds,
        usage_delta.system_seconds
    );
    Ok(())
}

fn payload_export_keys(db_path: &Path, output: &Path, limit: usize, skip: usize) -> Result<()> {
    let env = unsafe {
        EnvOpenOptions::new()
            .map_size(64 * 1024 * 1024 * 1024)
            .max_dbs(4)
            .open(db_path)
            .with_context(|| format!("failed to open heed env '{}'", db_path.display()))?
    };
    let rtxn = env.read_txn()?;
    let db: Database<payload_sidecar::U32KeyCodec, heed::types::Bytes> = env
        .open_database(&rtxn, Some("payloads"))?
        .context("heed payload database is missing")?;
    let mut iter = db.iter(&rtxn)?;
    let mut written = 0usize;
    let mut seen = 0usize;
    let mut text = String::with_capacity(limit.saturating_mul(12));

    while let Some((position, _bytes)) = iter.next().transpose()? {
        if seen < skip {
            seen += 1;
            continue;
        }
        if written >= limit {
            break;
        }
        text.push_str(&position.to_string());
        text.push('\n');
        written += 1;
    }

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory '{}'", parent.display()))?;
    }
    std::fs::write(output, text)
        .with_context(|| format!("failed to write positions file '{}'", output.display()))?;
    println!(
        "wrote={} skip={} limit={} output={}",
        written,
        skip,
        limit,
        output.display()
    );
    Ok(())
}

#[derive(Debug)]
struct VcfPositionMatchReport {
    chrom: String,
    dataset_path: PathBuf,
    vcf_path: PathBuf,
    cache_rows: usize,
    cache_warm_rows: usize,
    cache_cold_rows: usize,
    cache_other_tier_rows: usize,
    cache_unique_positions: usize,
    cache_warm_unique_positions: usize,
    cache_cold_unique_positions: usize,
    cache_warm_cold_overlap_positions: usize,
    vcf_records: usize,
    vcf_unique_positions: usize,
    matched_records_any: usize,
    matched_records_warm: usize,
    matched_records_cold: usize,
    matched_unique_any: usize,
    matched_unique_warm: usize,
    matched_unique_cold: usize,
    matched_unique_both_tiers: usize,
    vcf_positions: Vec<u32>,
    matched_cold_positions: Vec<u32>,
}

#[derive(Debug, Default)]
struct CachePositionSets {
    all_positions: Vec<u32>,
    warm_positions: Vec<u32>,
    cold_positions: Vec<u32>,
    rows: usize,
    warm_rows: usize,
    cold_rows: usize,
    other_tier_rows: usize,
}

async fn vcf_position_match(
    dataset_path: &Path,
    vcf_path: &Path,
    chrom: &str,
    position_column: &str,
    tier_column: &str,
) -> Result<VcfPositionMatchReport> {
    let cache = load_cache_position_sets(dataset_path, position_column, tier_column).await?;
    let vcf = scan_vcf_position_matches(
        vcf_path,
        chrom,
        &cache.all_positions,
        &cache.warm_positions,
        &cache.cold_positions,
    )?;
    let warm_cold_overlap_positions =
        count_sorted_intersection(&cache.warm_positions, &cache.cold_positions);

    Ok(VcfPositionMatchReport {
        chrom: chrom.to_string(),
        dataset_path: dataset_path.to_path_buf(),
        vcf_path: vcf_path.to_path_buf(),
        cache_rows: cache.rows,
        cache_warm_rows: cache.warm_rows,
        cache_cold_rows: cache.cold_rows,
        cache_other_tier_rows: cache.other_tier_rows,
        cache_unique_positions: cache.all_positions.len(),
        cache_warm_unique_positions: cache.warm_positions.len(),
        cache_cold_unique_positions: cache.cold_positions.len(),
        cache_warm_cold_overlap_positions: warm_cold_overlap_positions,
        vcf_records: vcf.records,
        vcf_unique_positions: vcf.unique_positions,
        matched_records_any: vcf.matched_records_any,
        matched_records_warm: vcf.matched_records_warm,
        matched_records_cold: vcf.matched_records_cold,
        matched_unique_any: vcf.matched_unique_any,
        matched_unique_warm: vcf.matched_unique_warm,
        matched_unique_cold: vcf.matched_unique_cold,
        matched_unique_both_tiers: vcf.matched_unique_both_tiers,
        vcf_positions: vcf.positions,
        matched_cold_positions: vcf.matched_cold_positions,
    })
}

async fn load_cache_position_sets(
    dataset_path: &Path,
    position_column: &str,
    tier_column: &str,
) -> Result<CachePositionSets> {
    let dataset = Dataset::open(dataset_path.to_string_lossy().as_ref())
        .await
        .with_context(|| format!("failed to open Lance dataset '{}'", dataset_path.display()))?;
    let mut scanner = dataset.scan();
    scanner
        .project(&[position_column, tier_column])?
        .materialization_style(MaterializationStyle::AllLate);
    let mut stream = scanner.try_into_stream().await?;
    let mut out = CachePositionSets::default();

    while let Some(batch) = stream.try_next().await? {
        let schema = batch.schema();
        let position_index = schema.index_of(position_column)?;
        let tier_index = schema.index_of(tier_column)?;
        let positions = PositionArray::from_array(batch.column(position_index).as_ref())?;
        let tiers = TierArray::from_array(batch.column(tier_index).as_ref())?;

        out.rows += batch.num_rows();
        out.all_positions.reserve(batch.num_rows());
        for row in 0..batch.num_rows() {
            let Some(position) = positions.get(row)? else {
                continue;
            };
            out.all_positions.push(position);
            match tiers.get(row)? {
                Some(0) => {
                    out.warm_rows += 1;
                    out.warm_positions.push(position);
                }
                Some(1) => {
                    out.cold_rows += 1;
                    out.cold_positions.push(position);
                }
                _ => {
                    out.other_tier_rows += 1;
                }
            }
        }
    }

    sort_dedup(&mut out.all_positions);
    sort_dedup(&mut out.warm_positions);
    sort_dedup(&mut out.cold_positions);
    Ok(out)
}

#[derive(Debug, Default)]
struct VcfPositionCounts {
    records: usize,
    unique_positions: usize,
    matched_records_any: usize,
    matched_records_warm: usize,
    matched_records_cold: usize,
    matched_unique_any: usize,
    matched_unique_warm: usize,
    matched_unique_cold: usize,
    matched_unique_both_tiers: usize,
    positions: Vec<u32>,
    matched_cold_positions: Vec<u32>,
}

fn scan_vcf_position_matches(
    vcf_path: &Path,
    chrom: &str,
    all_cache_positions: &[u32],
    warm_cache_positions: &[u32],
    cold_cache_positions: &[u32],
) -> Result<VcfPositionCounts> {
    let file =
        File::open(vcf_path).with_context(|| format!("failed to open '{}'", vcf_path.display()))?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let mut seen_target = false;
    let mut counts = VcfPositionCounts::default();
    let mut unique_positions = HashSet::new();
    let mut matched_any = HashSet::new();
    let mut matched_warm = HashSet::new();
    let mut matched_cold = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let Some(row_chrom) = fields.next() else {
            continue;
        };
        let row_chrom_bare = row_chrom.strip_prefix("chr").unwrap_or(row_chrom);
        if row_chrom_bare != chrom_bare {
            if seen_target {
                break;
            }
            continue;
        }
        seen_target = true;
        let Some(pos_str) = fields.next() else {
            continue;
        };
        let position = pos_str
            .parse::<u32>()
            .with_context(|| format!("invalid VCF POS '{pos_str}'"))?;

        counts.records += 1;
        unique_positions.insert(position);
        let in_warm = warm_cache_positions.binary_search(&position).is_ok();
        let in_cold = cold_cache_positions.binary_search(&position).is_ok();
        let in_any = in_warm || in_cold || all_cache_positions.binary_search(&position).is_ok();

        if in_any {
            counts.matched_records_any += 1;
            matched_any.insert(position);
        }
        if in_warm {
            counts.matched_records_warm += 1;
            matched_warm.insert(position);
        }
        if in_cold {
            counts.matched_records_cold += 1;
            matched_cold.insert(position);
        }
    }

    counts.unique_positions = unique_positions.len();
    counts.matched_unique_any = matched_any.len();
    counts.matched_unique_warm = matched_warm.len();
    counts.matched_unique_cold = matched_cold.len();
    counts.matched_unique_both_tiers = matched_warm.intersection(&matched_cold).count();
    counts.positions = unique_positions.into_iter().collect();
    counts.positions.sort_unstable();
    counts.matched_cold_positions = matched_cold.into_iter().collect();
    counts.matched_cold_positions.sort_unstable();
    Ok(counts)
}

fn write_positions_file(path: &Path, positions: &[u32]) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create '{}'", parent.display()))?;
    }
    let mut out = String::with_capacity(positions.len().saturating_mul(12));
    for position in positions {
        out.push_str(&position.to_string());
        out.push('\n');
    }
    std::fs::write(path, out).with_context(|| format!("failed to write '{}'", path.display()))
}

fn render_vcf_position_match(report: &VcfPositionMatchReport) -> String {
    let mut out = String::new();
    out.push_str("# VCF position match\n\n");
    out.push_str(&format!(
        "- VCF: `{}`\n- Lance dataset: `{}`\n- Chromosome: `{}`\n\n",
        report.vcf_path.display(),
        report.dataset_path.display(),
        report.chrom
    ));
    out.push_str("| Scope | Rows/records | Unique positions | Matched records | Matched unique positions | Match rate unique |\n");
    out.push_str("|---|---:|---:|---:|---:|---:|\n");
    out.push_str(&format!(
        "| Lance union | {} | {} | - | - | - |\n",
        fmt_usize(report.cache_rows),
        fmt_usize(report.cache_unique_positions)
    ));
    out.push_str(&format!(
        "| Lance warm tier | {} | {} | - | - | - |\n",
        fmt_usize(report.cache_warm_rows),
        fmt_usize(report.cache_warm_unique_positions)
    ));
    out.push_str(&format!(
        "| Lance cold tier | {} | {} | - | - | - |\n",
        fmt_usize(report.cache_cold_rows),
        fmt_usize(report.cache_cold_unique_positions)
    ));
    out.push_str(&format!(
        "| VCF chr1 vs union | {} | {} | {} | {} | {:.4}% |\n",
        fmt_usize(report.vcf_records),
        fmt_usize(report.vcf_unique_positions),
        fmt_usize(report.matched_records_any),
        fmt_usize(report.matched_unique_any),
        percentage(report.matched_unique_any, report.vcf_unique_positions)
    ));
    out.push_str(&format!(
        "| VCF chr1 vs warm | {} | {} | {} | {} | {:.4}% |\n",
        fmt_usize(report.vcf_records),
        fmt_usize(report.vcf_unique_positions),
        fmt_usize(report.matched_records_warm),
        fmt_usize(report.matched_unique_warm),
        percentage(report.matched_unique_warm, report.vcf_unique_positions)
    ));
    out.push_str(&format!(
        "| VCF chr1 vs cold | {} | {} | {} | {} | {:.4}% |\n\n",
        fmt_usize(report.vcf_records),
        fmt_usize(report.vcf_unique_positions),
        fmt_usize(report.matched_records_cold),
        fmt_usize(report.matched_unique_cold),
        percentage(report.matched_unique_cold, report.vcf_unique_positions)
    ));
    out.push_str(&format!(
        "- Warm/cold cache overlap positions: `{}`\n",
        fmt_usize(report.cache_warm_cold_overlap_positions)
    ));
    out.push_str(&format!(
        "- Matched VCF unique positions present in both warm and cold: `{}`\n",
        fmt_usize(report.matched_unique_both_tiers)
    ));
    if report.cache_other_tier_rows > 0 {
        out.push_str(&format!(
            "- Lance rows with tier other than 0/1: `{}`\n",
            fmt_usize(report.cache_other_tier_rows)
        ));
    }
    out
}

fn sort_dedup(values: &mut Vec<u32>) {
    values.sort_unstable();
    values.dedup();
}

fn count_sorted_intersection(left: &[u32], right: &[u32]) -> usize {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut count = 0usize;
    while i < left.len() && j < right.len() {
        match left[i].cmp(&right[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
        }
    }
    count
}

fn fmt_usize(value: usize) -> String {
    let raw = value.to_string();
    let mut out = String::with_capacity(raw.len() + raw.len() / 3);
    for (idx, ch) in raw.chars().rev().enumerate() {
        if idx > 0 && idx % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}

fn percentage(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        (numerator as f64 / denominator as f64) * 100.0
    }
}

enum PositionArray<'a> {
    U32(&'a UInt32Array),
    U64(&'a UInt64Array),
    I64(&'a Int64Array),
    I32(&'a Int32Array),
}

impl<'a> PositionArray<'a> {
    fn from_array(array: &'a dyn Array) -> Result<Self> {
        match array.data_type() {
            DataType::UInt32 => Ok(Self::U32(
                array.as_any().downcast_ref::<UInt32Array>().unwrap(),
            )),
            DataType::UInt64 => Ok(Self::U64(
                array.as_any().downcast_ref::<UInt64Array>().unwrap(),
            )),
            DataType::Int64 => Ok(Self::I64(
                array.as_any().downcast_ref::<Int64Array>().unwrap(),
            )),
            DataType::Int32 => Ok(Self::I32(
                array.as_any().downcast_ref::<Int32Array>().unwrap(),
            )),
            other => bail!("unsupported position column type: {other:?}"),
        }
    }

    fn get(&self, row: usize) -> Result<Option<u32>> {
        match self {
            Self::U32(array) => Ok((!array.is_null(row)).then(|| array.value(row))),
            Self::U64(array) => {
                if array.is_null(row) {
                    Ok(None)
                } else {
                    Ok(Some(u32::try_from(array.value(row))?))
                }
            }
            Self::I64(array) => {
                if array.is_null(row) {
                    Ok(None)
                } else {
                    Ok(Some(u32::try_from(array.value(row))?))
                }
            }
            Self::I32(array) => {
                if array.is_null(row) {
                    Ok(None)
                } else {
                    Ok(Some(u32::try_from(array.value(row))?))
                }
            }
        }
    }
}

enum TierArray<'a> {
    U8(&'a UInt8Array),
    I8(&'a Int8Array),
}

impl<'a> TierArray<'a> {
    fn from_array(array: &'a dyn Array) -> Result<Self> {
        match array.data_type() {
            DataType::UInt8 => Ok(Self::U8(
                array.as_any().downcast_ref::<UInt8Array>().unwrap(),
            )),
            DataType::Int8 => Ok(Self::I8(
                array.as_any().downcast_ref::<Int8Array>().unwrap(),
            )),
            other => bail!("unsupported tier column type: {other:?}"),
        }
    }

    fn get(&self, row: usize) -> Result<Option<u8>> {
        match self {
            Self::U8(array) => Ok((!array.is_null(row)).then(|| array.value(row))),
            Self::I8(array) => {
                if array.is_null(row) {
                    Ok(None)
                } else {
                    Ok(Some(u8::try_from(array.value(row))?))
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct ResourceUsage {
    minor_faults: i64,
    major_faults: i64,
    user_seconds: f64,
    system_seconds: f64,
}

impl ResourceUsage {
    fn snapshot() -> Self {
        unsafe {
            let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
            let rc = libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr());
            if rc != 0 {
                return Self {
                    minor_faults: 0,
                    major_faults: 0,
                    user_seconds: 0.0,
                    system_seconds: 0.0,
                };
            }
            let usage = usage.assume_init();
            Self {
                minor_faults: usage.ru_minflt,
                major_faults: usage.ru_majflt,
                user_seconds: timeval_seconds(usage.ru_utime),
                system_seconds: timeval_seconds(usage.ru_stime),
            }
        }
    }

    fn delta(self, before: Self) -> Self {
        Self {
            minor_faults: self.minor_faults - before.minor_faults,
            major_faults: self.major_faults - before.major_faults,
            user_seconds: self.user_seconds - before.user_seconds,
            system_seconds: self.system_seconds - before.system_seconds,
        }
    }
}

fn timeval_seconds(time: libc::timeval) -> f64 {
    time.tv_sec as f64 + (time.tv_usec as f64 / 1_000_000.0)
}

fn read_positions_file(path: &Path) -> Result<Vec<u64>> {
    let raw = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read positions file '{}'", path.display()))?;
    raw.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| Ok(line.parse::<u64>()?))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_accepts_positions_file() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "run",
            "--config",
            "config.toml",
            "--positions-file",
            "positions.txt",
        ])
        .unwrap();

        match cli.command {
            Command::Run {
                config,
                positions_file,
                lookup_tier,
            } => {
                assert_eq!(config, PathBuf::from("config.toml"));
                assert_eq!(positions_file, Some(PathBuf::from("positions.txt")));
                assert_eq!(lookup_tier, "1");
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn inspect_path_accepts_dataset_path_and_lance_version() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "inspect-path",
            "--dataset-path",
            "/tmp/chr1.lance",
            "--lance-version",
            "2.1",
        ])
        .unwrap();

        match cli.command {
            Command::InspectPath {
                dataset_path,
                lance_version,
                input_parquet,
                position_field,
                position_source_column,
                include_physical_columns,
                include_index_sizes,
            } => {
                assert_eq!(dataset_path, PathBuf::from("/tmp/chr1.lance"));
                assert_eq!(lance_version, "2.1");
                assert!(input_parquet.is_empty());
                assert_eq!(position_field, "position");
                assert_eq!(position_source_column, "start");
                assert!(include_physical_columns);
                assert!(include_index_sizes);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn bench_accepts_dataset_path_override() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "bench",
            "--config",
            "config.toml",
            "--positions-file",
            "positions.txt",
            "--dataset-path",
            "/tmp/notebook_chr1.lance",
            "--lookup-tier",
            "0",
        ])
        .unwrap();

        match cli.command {
            Command::Bench {
                config,
                positions_file,
                dataset_path,
                lookup_tier,
                lookup_tiers,
                only_index_take_everything,
            } => {
                assert_eq!(config, PathBuf::from("config.toml"));
                assert_eq!(positions_file, Some(PathBuf::from("positions.txt")));
                assert_eq!(
                    dataset_path,
                    Some(PathBuf::from("/tmp/notebook_chr1.lance"))
                );
                assert_eq!(lookup_tier, "0");
                assert!(lookup_tiers.is_empty());
                assert!(!only_index_take_everything);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn bench_accepts_lookup_tiers() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "bench",
            "--config",
            "config.toml",
            "--positions-file",
            "positions.txt",
            "--lookup-tiers",
            "0,1,all",
        ])
        .unwrap();

        match cli.command {
            Command::Bench {
                config,
                positions_file,
                dataset_path,
                lookup_tier,
                lookup_tiers,
                only_index_take_everything,
            } => {
                assert_eq!(config, PathBuf::from("config.toml"));
                assert_eq!(positions_file, Some(PathBuf::from("positions.txt")));
                assert_eq!(dataset_path, None);
                assert_eq!(lookup_tier, "1");
                assert_eq!(lookup_tiers, vec!["0", "1", "all"]);
                assert!(!only_index_take_everything);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn bench_accepts_only_index_take_everything() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "bench",
            "--config",
            "config.toml",
            "--positions-file",
            "positions.txt",
            "--only-index-take-everything",
        ])
        .unwrap();

        match cli.command {
            Command::Bench {
                config,
                positions_file,
                dataset_path,
                lookup_tier,
                lookup_tiers,
                only_index_take_everything,
            } => {
                assert_eq!(config, PathBuf::from("config.toml"));
                assert_eq!(positions_file, Some(PathBuf::from("positions.txt")));
                assert_eq!(dataset_path, None);
                assert_eq!(lookup_tier, "1");
                assert!(lookup_tiers.is_empty());
                assert!(only_index_take_everything);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn payload_bench_accepts_raw_only() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "payload-bench",
            "--db-path",
            "/tmp/payload.lmdb",
            "--positions-file",
            "positions.txt",
            "--raw-only",
            "--no-read-ahead",
        ])
        .unwrap();

        match cli.command {
            Command::PayloadBench {
                db_path,
                positions_file,
                materialize,
                no_read_ahead,
                raw_only,
            } => {
                assert_eq!(db_path, PathBuf::from("/tmp/payload.lmdb"));
                assert_eq!(positions_file, PathBuf::from("positions.txt"));
                assert!(!materialize);
                assert!(no_read_ahead);
                assert!(raw_only);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn payload_export_keys_accepts_limit_and_skip() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "payload-export-keys",
            "--db-path",
            "/tmp/payload.lmdb",
            "--output",
            "/tmp/keys.txt",
            "--limit",
            "100000",
            "--skip",
            "500000",
        ])
        .unwrap();

        match cli.command {
            Command::PayloadExportKeys {
                db_path,
                output,
                limit,
                skip,
            } => {
                assert_eq!(db_path, PathBuf::from("/tmp/payload.lmdb"));
                assert_eq!(output, PathBuf::from("/tmp/keys.txt"));
                assert_eq!(limit, 100_000);
                assert_eq!(skip, 500_000);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn vcf_position_match_accepts_inputs() {
        let cli = Cli::try_parse_from([
            "lance_sandbox",
            "vcf-position-match",
            "--dataset-path",
            "/tmp/chr1.lance",
            "--vcf-path",
            "/tmp/input.vcf.gz",
            "--chrom",
            "chr1",
            "--write-all-positions",
            "/tmp/all.txt",
            "--write-matched-cold-positions",
            "/tmp/cold.txt",
        ])
        .unwrap();

        match cli.command {
            Command::VcfPositionMatch {
                dataset_path,
                vcf_path,
                chrom,
                position_column,
                tier_column,
                write_all_positions,
                write_matched_cold_positions,
            } => {
                assert_eq!(dataset_path, PathBuf::from("/tmp/chr1.lance"));
                assert_eq!(vcf_path, PathBuf::from("/tmp/input.vcf.gz"));
                assert_eq!(chrom, "chr1");
                assert_eq!(position_column, "position");
                assert_eq!(tier_column, "tier");
                assert_eq!(write_all_positions, Some(PathBuf::from("/tmp/all.txt")));
                assert_eq!(
                    write_matched_cold_positions,
                    Some(PathBuf::from("/tmp/cold.txt"))
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }
}
