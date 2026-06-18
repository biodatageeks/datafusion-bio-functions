use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail};
use arrow_array::{Array, UInt32Array, UInt64Array};
use arrow_schema::DataType;
use datafusion_common::ScalarValue;
use datafusion_physical_plan::displayable;
use futures::TryStreamExt;
use lance::Dataset;
use lance::dataset::ProjectionRequest;
use lance::dataset::index::LanceIndexStoreExt;
use lance::dataset::scanner::{ExecutionSummaryCounts, MaterializationStyle};
use lance::index::{DatasetIndexExt, DatasetIndexInternalExt};
use lance::table::format::IndexMetadata;
use lance_core::ROW_ID;
use lance_core::utils::mask::{NullableRowAddrSet, RowAddrTreeMap, RowSetOps};
use lance_index::IndexCriteria;
use lance_index::metrics::LocalMetricsCollector;
use lance_index::scalar::{IndexStore, SargableQuery, SearchResult, lance_format::LanceIndexStore};
use lance_io::{bytes_read_counter, iops_counter};
use serde::Serialize;

use crate::build::{key_values_from_batch, physical_projection_for_everything};
use crate::config::{KeyDataType, SandboxConfig};
use crate::row_sidecar::{PositionRowIdIndex, PositionRowIdIndexBuilder};

const BTREE_PAGE_DATA_FILE: &str = "page_data.lance";
const BTREE_VALUES_COLUMN: &str = "values";
const BTREE_IDS_COLUMN: &str = "ids";

#[derive(Debug, Default, Clone, Serialize)]
pub struct ScanIoStats {
    pub iops: usize,
    pub requests: usize,
    pub bytes_read: usize,
    pub indices_loaded: usize,
    pub parts_loaded: usize,
    pub index_comparisons: usize,
    pub fragments_scanned: usize,
    pub ranges_scanned: usize,
    pub rows_scanned: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ScanResult {
    pub name: String,
    pub filter_keys: usize,
    pub lookup_batch_size: Option<usize>,
    pub scans: usize,
    pub rows: usize,
    pub result_batches: usize,
    pub selected_positions: usize,
    pub seconds: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub row_id_resolve_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub btree_resolve_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub take_rows_seconds: Option<f64>,
    pub io: ScanIoStats,
}

#[derive(Debug, Clone, Serialize)]
pub struct PhysicalPlanReport {
    pub name: String,
    pub filter: String,
    pub lookup_batch_size: Option<usize>,
    pub projected_column_count: usize,
    pub projection: Vec<String>,
    pub plan: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct BenchmarkReport {
    pub dataset_path: String,
    pub lookup_tier: Option<u8>,
    pub lookup_scope: String,
    pub requested_everything_fields: Vec<String>,
    pub resolved_dataset_projection: Vec<String>,
    pub projected_column_count: usize,
    pub positions_requested: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warm_result: Option<ScanResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cold_prewarm_result: Option<ScanResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sidecar: Option<PositionSidecarReport>,
    pub cold_results: Vec<ScanResult>,
    pub physical_plans: Vec<PhysicalPlanReport>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PositionSidecarReport {
    pub path: String,
    pub built: bool,
    pub build_seconds: f64,
    pub load_seconds: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sidecar_read_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub btree_page_load_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub btree_metadata_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tier_filter_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub btree_page_scan_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub btree_page_stream_batches: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub btree_page_stream_rows: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pair_compare_seconds: Option<f64>,
    pub unique_positions: usize,
    pub row_ids: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BenchmarkMode {
    Full,
    OnlyIndexTakeEverything,
}

impl BenchmarkMode {
    fn run_warm_scan(self) -> bool {
        matches!(self, Self::Full)
    }

    fn run_cold_scanner_paths(self) -> bool {
        matches!(self, Self::Full)
    }

    fn run_unbatched_btree_take(self) -> bool {
        matches!(self, Self::Full)
    }

    fn run_position_only_take(self) -> bool {
        matches!(self, Self::Full)
    }

    fn capture_physical_plans(self) -> bool {
        matches!(self, Self::Full)
    }
}

pub async fn run_benchmark_with_progress<F>(
    config: &SandboxConfig,
    positions_file: &Path,
    dataset_path_override: Option<&Path>,
    lookup_tier: Option<u8>,
    on_update: F,
) -> Result<BenchmarkReport>
where
    F: FnMut(&BenchmarkReport) -> Result<()>,
{
    run_benchmark_with_mode(
        config,
        positions_file,
        dataset_path_override,
        lookup_tier,
        BenchmarkMode::Full,
        on_update,
    )
    .await
}

pub async fn run_benchmark_with_mode<F>(
    config: &SandboxConfig,
    positions_file: &Path,
    dataset_path_override: Option<&Path>,
    lookup_tier: Option<u8>,
    mode: BenchmarkMode,
    mut on_update: F,
) -> Result<BenchmarkReport>
where
    F: FnMut(&BenchmarkReport) -> Result<()>,
{
    let dataset_path = benchmark_dataset_path(config, dataset_path_override);
    let dataset = Dataset::open(dataset_path.to_string_lossy().as_ref())
        .await
        .with_context(|| format!("failed to open Lance dataset '{}'", dataset_path.display()))?;
    let projection = physical_projection_for_everything(config, dataset.schema())?;
    let requested_everything_fields = crate::build::everything_logical_fields(config);
    let mut positions = read_positions(positions_file)?;
    if config.benchmark.sort_position_keys {
        positions.sort_unstable();
    }
    positions.dedup();

    eprintln!(
        "Benchmarking `{}` with {} positions and {} projected physical fields.",
        dataset_path.display(),
        positions.len(),
        projection.len()
    );

    let mut report = BenchmarkReport {
        dataset_path: dataset_path.display().to_string(),
        lookup_tier,
        lookup_scope: lookup_scope_label(lookup_tier),
        requested_everything_fields,
        resolved_dataset_projection: projection.clone(),
        projected_column_count: projection.len(),
        positions_requested: positions.len(),
        warm_result: None,
        cold_prewarm_result: None,
        sidecar: None,
        cold_results: Vec::new(),
        physical_plans: Vec::new(),
    };

    if mode.run_warm_scan() {
        eprintln!("Running warm_full_scan_all_everything_columns...");
        let warm_result = scan_once(
            config,
            &dataset,
            "warm_full_scan_all_everything_columns",
            "tier = 0".to_string(),
            None,
            &projection,
        )
        .await?;
        eprintln!("Finished {}", scan_progress_summary(&warm_result));
        report.warm_result = Some(warm_result);
    }
    on_update(&report)?;

    if mode.run_cold_scanner_paths() {
        if let Some(&lookup_batch_size) = config.benchmark.lookup_batch_sizes.first() {
            if !positions.is_empty() {
                let prewarm_scans = positions.chunks(lookup_batch_size.max(1)).count();
                let key_only_projection = vec![config.key.column.clone()];
                eprintln!(
                    "Prewarming cold index cache: lookup_batch_size={}, scans={}, keys={}, projection={}",
                    lookup_batch_size,
                    prewarm_scans,
                    positions.len(),
                    config.key.column
                );
                let mut prewarm_result = scan_position_chunks(
                    config,
                    &dataset,
                    &positions,
                    lookup_batch_size,
                    &key_only_projection,
                    lookup_tier,
                )
                .await?;
                prewarm_result.name =
                    format!("{}_index_prewarm_key_only", lookup_label(lookup_tier));
                eprintln!("Finished {}", scan_progress_summary(&prewarm_result));
                report.cold_prewarm_result = Some(prewarm_result);
                on_update(&report)?;
            }
        }

        for (idx, &lookup_batch_size) in config.benchmark.lookup_batch_sizes.iter().enumerate() {
            let scans = positions.chunks(lookup_batch_size.max(1)).count();
            eprintln!(
                "Running cold batch {}/{}: lookup_batch_size={}, scans={}, keys={}",
                idx + 1,
                config.benchmark.lookup_batch_sizes.len(),
                lookup_batch_size,
                scans,
                positions.len()
            );
            let cold_result = scan_position_chunks(
                config,
                &dataset,
                &positions,
                lookup_batch_size,
                &projection,
                lookup_tier,
            )
            .await?;
            eprintln!("Finished {}", scan_progress_summary(&cold_result));
            report.cold_results.push(cold_result);
            on_update(&report)?;

            eprintln!(
                "Running row-id scan + take_rows batch {}/{}: lookup_batch_size={}, scans={}, keys={}",
                idx + 1,
                config.benchmark.lookup_batch_sizes.len(),
                lookup_batch_size,
                scans,
                positions.len()
            );
            let row_id_take_result = scan_row_ids_then_take(
                config,
                &dataset,
                &positions,
                lookup_batch_size,
                &projection,
                "everything",
                lookup_tier,
            )
            .await?;
            eprintln!("Finished {}", scan_progress_summary(&row_id_take_result));
            report.cold_results.push(row_id_take_result);
            on_update(&report)?;
        }
    }

    let loaded_sidecar = if !positions.is_empty() {
        eprintln!("Loading or building external position -> row_id sidecar...");
        let sidecar =
            load_or_build_position_sidecar(config, &dataset_path, &dataset, lookup_tier).await?;
        eprintln!(
            "Sidecar ready: unique_positions={}, row_ids={}, built={}, build_seconds={:.3}",
            sidecar.index.unique_positions(),
            sidecar.index.row_ids_len(),
            sidecar.report.built,
            sidecar.report.build_seconds
        );
        report.sidecar = Some(sidecar.report.clone());
        on_update(&report)?;
        Some(sidecar)
    } else {
        None
    };

    if let Some(sidecar) = loaded_sidecar.as_ref() {
        let btree_row_index = sidecar.btree_index.as_ref().unwrap_or(&sidecar.index);
        if mode.run_unbatched_btree_take() {
            eprintln!("Running BTree page-data row map + take_rows everything projection...");
            let btree_everything = btree_direct_take_rows(
                config,
                &dataset,
                btree_row_index,
                &positions,
                &projection,
                "everything",
                lookup_tier,
            )
            .await?;
            eprintln!("Finished {}", scan_progress_summary(&btree_everything));
            report.cold_results.push(btree_everything);
            on_update(&report)?;
        }

        for &lookup_batch_size in &config.benchmark.lookup_batch_sizes {
            let scans = positions.chunks(lookup_batch_size.max(1)).count();
            eprintln!(
                "Running BTree page-data row map + take_rows everything projection batch: lookup_batch_size={}, scans={}, keys={}",
                lookup_batch_size,
                scans,
                positions.len()
            );
            let btree_everything_batched = btree_direct_take_rows_batched(
                config,
                &dataset,
                btree_row_index,
                &positions,
                lookup_batch_size,
                &projection,
                "everything",
                lookup_tier,
            )
            .await?;
            eprintln!(
                "Finished {}",
                scan_progress_summary(&btree_everything_batched)
            );
            report.cold_results.push(btree_everything_batched);
            on_update(&report)?;
        }

        if mode.run_position_only_take() {
            eprintln!("Running BTree page-data row map + take_rows position-only projection...");
            let position_projection = vec![config.key.column.clone()];
            let btree_position = btree_direct_take_rows(
                config,
                &dataset,
                btree_row_index,
                &positions,
                &position_projection,
                "position_only",
                lookup_tier,
            )
            .await?;
            eprintln!("Finished {}", scan_progress_summary(&btree_position));
            report.cold_results.push(btree_position);
            on_update(&report)?;
        }
    }

    if let Some(sidecar) = loaded_sidecar.as_ref() {
        eprintln!("Running sidecar + take_rows everything projection...");
        let sidecar_everything = sidecar_take_rows(
            config,
            &dataset,
            &sidecar.index,
            &positions,
            &projection,
            "everything",
            lookup_tier,
        )
        .await?;
        eprintln!("Finished {}", scan_progress_summary(&sidecar_everything));
        report.cold_results.push(sidecar_everything);
        on_update(&report)?;

        for &lookup_batch_size in &config.benchmark.lookup_batch_sizes {
            let scans = positions.chunks(lookup_batch_size.max(1)).count();
            eprintln!(
                "Running sidecar + take_rows everything projection batch: lookup_batch_size={}, scans={}, keys={}",
                lookup_batch_size,
                scans,
                positions.len()
            );
            let sidecar_everything_batched = sidecar_take_rows_batched(
                config,
                &dataset,
                &sidecar.index,
                &positions,
                lookup_batch_size,
                &projection,
                "everything",
                lookup_tier,
            )
            .await?;
            eprintln!(
                "Finished {}",
                scan_progress_summary(&sidecar_everything_batched)
            );
            report.cold_results.push(sidecar_everything_batched);
            on_update(&report)?;
        }

        if mode.run_position_only_take() {
            eprintln!("Running sidecar + take_rows position-only projection...");
            let position_projection = vec![config.key.column.clone()];
            let sidecar_position = sidecar_take_rows(
                config,
                &dataset,
                &sidecar.index,
                &positions,
                &position_projection,
                "position_only",
                lookup_tier,
            )
            .await?;
            eprintln!("Finished {}", scan_progress_summary(&sidecar_position));
            report.cold_results.push(sidecar_position);
            on_update(&report)?;
        }
    }

    if mode.capture_physical_plans() {
        eprintln!("Capturing representative Lance physical plans...");
        report.physical_plans =
            capture_physical_plans(config, &dataset, &positions, &projection, lookup_tier).await?;
        on_update(&report)?;
    }

    Ok(report)
}

fn benchmark_dataset_path(config: &SandboxConfig, dataset_path_override: Option<&Path>) -> PathBuf {
    dataset_path_override
        .map(Path::to_path_buf)
        .unwrap_or_else(|| config.dataset_path())
}

async fn scan_position_chunks(
    config: &SandboxConfig,
    dataset: &Dataset,
    positions: &[u64],
    lookup_batch_size: usize,
    projection: &[String],
    lookup_tier: Option<u8>,
) -> Result<ScanResult> {
    let started = Instant::now();
    let mut rows = 0usize;
    let mut result_batches = 0usize;
    let mut selected_positions = HashSet::new();
    let mut io = ScanIoStats::default();
    let mut scans = 0usize;
    let mut filter_keys = 0usize;
    let projection_refs = projection.iter().map(String::as_str).collect::<Vec<_>>();

    for chunk in positions.chunks(lookup_batch_size.max(1)) {
        scans += 1;
        filter_keys += chunk.len();
        let filter = position_filter(config, lookup_tier, chunk);
        let stats_slot = Arc::new(Mutex::new(None::<ExecutionSummaryCounts>));
        let stats_sink = Arc::clone(&stats_slot);
        let mut scanner = dataset.scan();
        scanner
            .filter(&filter)?
            .project(&projection_refs)?
            .materialization_style(MaterializationStyle::AllLate)
            .scan_stats_callback(Arc::new(move |stats| {
                *stats_sink.lock().expect("stats lock poisoned") = Some(stats.clone());
            }));
        let mut stream = scanner.try_into_stream().await?;
        while let Some(batch) = stream.try_next().await? {
            rows += batch.num_rows();
            result_batches += 1;
            selected_positions.extend(key_values_from_batch(config, &batch)?);
        }
        if let Some(stats) = stats_slot.lock().expect("stats lock poisoned").as_ref() {
            io.add_execution(stats);
        }
    }

    Ok(ScanResult {
        name: format!(
            "{}_lookup_batch_{lookup_batch_size}",
            lookup_label(lookup_tier)
        ),
        filter_keys,
        lookup_batch_size: Some(lookup_batch_size),
        scans,
        rows,
        result_batches,
        selected_positions: selected_positions.len(),
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: None,
        btree_resolve_seconds: None,
        take_rows_seconds: None,
        io,
    })
}

async fn scan_row_ids_then_take(
    config: &SandboxConfig,
    dataset: &Dataset,
    positions: &[u64],
    lookup_batch_size: usize,
    projection: &[String],
    projection_label: &str,
    lookup_tier: Option<u8>,
) -> Result<ScanResult> {
    let started = Instant::now();
    let resolve_started = Instant::now();
    let mut row_ids = Vec::new();
    let mut row_id_scan_batches = 0usize;
    let mut selected_positions = HashSet::new();
    let mut io = ScanIoStats::default();
    let mut scans = 0usize;
    let mut filter_keys = 0usize;
    let key_projection = [config.key.column.as_str()];

    for chunk in positions.chunks(lookup_batch_size.max(1)) {
        scans += 1;
        filter_keys += chunk.len();
        let filter = position_filter(config, lookup_tier, chunk);
        let stats_slot = Arc::new(Mutex::new(None::<ExecutionSummaryCounts>));
        let stats_sink = Arc::clone(&stats_slot);
        let mut scanner = dataset.scan();
        scanner
            .filter(&filter)?
            .project(&key_projection)?
            .with_row_id()
            .materialization_style(MaterializationStyle::AllLate)
            .scan_stats_callback(Arc::new(move |stats| {
                *stats_sink.lock().expect("stats lock poisoned") = Some(stats.clone());
            }));
        let mut stream = scanner.try_into_stream().await?;
        while let Some(batch) = stream.try_next().await? {
            row_id_scan_batches += 1;
            selected_positions.extend(key_values_from_batch(config, &batch)?);
            row_ids.extend(row_ids_from_batch(&batch)?);
        }
        if let Some(stats) = stats_slot.lock().expect("stats lock poisoned").as_ref() {
            io.add_execution(stats);
        }
    }

    let row_id_resolve_seconds = resolve_started.elapsed().as_secs_f64();
    let (batch, take_io, take_rows_seconds) =
        take_rows_projected(dataset, &row_ids, projection).await?;
    io.add(&take_io);
    selected_positions.extend(key_values_from_batch(config, &batch)?);

    Ok(ScanResult {
        name: format!(
            "{}_rowid_scan_then_take_{projection_label}_batch_{lookup_batch_size}",
            lookup_label(lookup_tier)
        ),
        filter_keys,
        lookup_batch_size: Some(lookup_batch_size),
        scans,
        rows: batch.num_rows(),
        result_batches: row_id_scan_batches + usize::from(batch.num_rows() > 0),
        selected_positions: selected_positions.len(),
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: Some(row_id_resolve_seconds),
        btree_resolve_seconds: None,
        take_rows_seconds: Some(take_rows_seconds),
        io,
    })
}

struct LoadedPositionSidecar {
    report: PositionSidecarReport,
    index: PositionRowIdIndex,
    btree_index: Option<PositionRowIdIndex>,
}

struct BtreePageLoad {
    index: PositionRowIdIndex,
    metadata_seconds: f64,
    tier_filter_seconds: f64,
    page_scan_seconds: f64,
    stream_batches: usize,
    stream_rows: usize,
}

#[derive(Debug, Default, Clone, Copy)]
struct BtreePageStreamStats {
    batches: usize,
    rows: usize,
}

impl BtreePageStreamStats {
    fn add(&mut self, other: Self) {
        self.batches += other.batches;
        self.rows += other.rows;
    }
}

async fn load_or_build_position_sidecar(
    config: &SandboxConfig,
    dataset_path: &Path,
    dataset: &Dataset,
    lookup_tier: Option<u8>,
) -> Result<LoadedPositionSidecar> {
    let path = position_sidecar_path(dataset_path, lookup_tier);
    let started = Instant::now();

    if path.exists() {
        let read_started = Instant::now();
        match fs::File::open(&path)
            .with_context(|| format!("failed to open row sidecar '{}'", path.display()))
            .and_then(|file| {
                PositionRowIdIndex::read_from(file)
                    .with_context(|| format!("failed to read row sidecar '{}'", path.display()))
            }) {
            Ok(index) => {
                let sidecar_read_seconds = read_started.elapsed().as_secs_f64();
                if config.indexes.position_btree {
                    eprintln!(
                        "Validating existing sidecar against BTree page-data row map for {}...",
                        lookup_scope_label(lookup_tier)
                    );
                    let btree_started = Instant::now();
                    let btree_load = build_position_sidecar_from_btree_pages_measured(
                        config,
                        dataset,
                        lookup_tier,
                    )
                    .await
                    .context("failed to build BTree page-data sidecar reference")?;
                    let btree_page_load_seconds = btree_started.elapsed().as_secs_f64();
                    let btree_index = btree_load.index;
                    let compare_started = Instant::now();
                    if let Some(mismatch) = index.pair_set_mismatch(&btree_index) {
                        bail!(
                            "existing sidecar '{}' does not match BTree page-data row map for {}: {mismatch}",
                            path.display(),
                            lookup_scope_label(lookup_tier)
                        );
                    }
                    let pair_compare_seconds = compare_started.elapsed().as_secs_f64();
                    eprintln!("Existing sidecar matches BTree page-data row map.");
                    let report = PositionSidecarReport {
                        path: path.display().to_string(),
                        built: false,
                        build_seconds: 0.0,
                        load_seconds: started.elapsed().as_secs_f64(),
                        sidecar_read_seconds: Some(sidecar_read_seconds),
                        btree_page_load_seconds: Some(btree_page_load_seconds),
                        btree_metadata_seconds: Some(btree_load.metadata_seconds),
                        tier_filter_seconds: Some(btree_load.tier_filter_seconds),
                        btree_page_scan_seconds: Some(btree_load.page_scan_seconds),
                        btree_page_stream_batches: Some(btree_load.stream_batches),
                        btree_page_stream_rows: Some(btree_load.stream_rows),
                        pair_compare_seconds: Some(pair_compare_seconds),
                        unique_positions: index.unique_positions(),
                        row_ids: index.row_ids_len(),
                    };
                    return Ok(LoadedPositionSidecar {
                        report,
                        index,
                        btree_index: Some(btree_index),
                    });
                }
                let report = PositionSidecarReport {
                    path: path.display().to_string(),
                    built: false,
                    build_seconds: 0.0,
                    load_seconds: started.elapsed().as_secs_f64(),
                    sidecar_read_seconds: Some(sidecar_read_seconds),
                    btree_page_load_seconds: None,
                    btree_metadata_seconds: None,
                    tier_filter_seconds: None,
                    btree_page_scan_seconds: None,
                    btree_page_stream_batches: None,
                    btree_page_stream_rows: None,
                    pair_compare_seconds: None,
                    unique_positions: index.unique_positions(),
                    row_ids: index.row_ids_len(),
                };
                return Ok(LoadedPositionSidecar {
                    report,
                    index,
                    btree_index: None,
                });
            }
            Err(error) => {
                eprintln!(
                    "Ignoring unreadable row sidecar '{}': {error:#}",
                    path.display()
                );
                let _ = fs::remove_file(&path);
            }
        }
    }

    let index = build_position_sidecar(config, dataset, lookup_tier).await?;
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create sidecar dir '{}'", parent.display()))?;
    }
    let temp_path = path.with_extension("bin.tmp");
    let mut file = fs::File::create(&temp_path)
        .with_context(|| format!("failed to create row sidecar '{}'", temp_path.display()))?;
    index
        .write_to(&mut file)
        .with_context(|| format!("failed to write row sidecar '{}'", temp_path.display()))?;
    file.sync_all()
        .with_context(|| format!("failed to sync row sidecar '{}'", temp_path.display()))?;
    drop(file);
    fs::rename(&temp_path, &path).with_context(|| {
        format!(
            "failed to move row sidecar '{}' to '{}'",
            temp_path.display(),
            path.display()
        )
    })?;

    let report = PositionSidecarReport {
        path: path.display().to_string(),
        built: true,
        build_seconds: started.elapsed().as_secs_f64(),
        load_seconds: 0.0,
        sidecar_read_seconds: None,
        btree_page_load_seconds: None,
        btree_metadata_seconds: None,
        tier_filter_seconds: None,
        btree_page_scan_seconds: None,
        btree_page_stream_batches: None,
        btree_page_stream_rows: None,
        pair_compare_seconds: None,
        unique_positions: index.unique_positions(),
        row_ids: index.row_ids_len(),
    };
    Ok(LoadedPositionSidecar {
        report,
        index,
        btree_index: None,
    })
}

fn position_sidecar_path(dataset_path: &Path, lookup_tier: Option<u8>) -> PathBuf {
    let file_name = match lookup_tier {
        Some(1) => "position_row_ids.bin".to_string(),
        Some(tier) => format!("position_row_ids_tier{tier}.bin"),
        None => "position_row_ids_all_tiers.bin".to_string(),
    };
    dataset_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("reports")
        .join(file_name)
}

async fn build_position_sidecar(
    config: &SandboxConfig,
    dataset: &Dataset,
    lookup_tier: Option<u8>,
) -> Result<PositionRowIdIndex> {
    if config.indexes.position_btree {
        let index = build_position_sidecar_from_btree_pages(config, dataset, lookup_tier).await?;
        eprintln!(
            "Validating BTree page-data row map against dataset scan for {}...",
            lookup_scope_label(lookup_tier)
        );
        let reference = build_position_sidecar_from_dataset_scan(config, dataset, lookup_tier)
            .await
            .context("failed to build dataset-scan sidecar reference")?;
        if let Some(mismatch) = index.pair_set_mismatch(&reference) {
            bail!(
                "BTree page-data row map does not match dataset-scan sidecar reference for {}: {mismatch}",
                lookup_scope_label(lookup_tier)
            );
        }
        eprintln!("BTree page-data row map matches dataset-scan sidecar reference.");
        return Ok(index);
    }

    build_position_sidecar_from_dataset_scan(config, dataset, lookup_tier).await
}

async fn build_position_sidecar_from_dataset_scan(
    config: &SandboxConfig,
    dataset: &Dataset,
    lookup_tier: Option<u8>,
) -> Result<PositionRowIdIndex> {
    let key_projection = [config.key.column.as_str()];
    let mut scanner = dataset.scan();
    if let Some(tier) = lookup_tier {
        scanner.filter(&format!("tier = {tier}"))?;
    }
    scanner
        .project(&key_projection)?
        .with_row_id()
        .materialization_style(MaterializationStyle::AllLate);

    let mut stream = scanner.try_into_stream().await?;
    if lookup_tier.is_some() {
        let mut builder = PositionRowIdIndexBuilder::new();
        while let Some(batch) = stream.try_next().await? {
            append_position_row_ids_to_builder(config, &batch, &mut builder)?;
        }
        return builder.finish();
    }

    let mut pairs = Vec::new();
    while let Some(batch) = stream.try_next().await? {
        append_position_row_id_pairs(config, &batch, &mut pairs)?;
    }
    pairs.sort_unstable_by_key(|&(position, row_id)| (position, row_id));
    let mut builder = PositionRowIdIndexBuilder::new();
    for (position, row_id) in pairs {
        builder.push_monotonic(position, row_id)?;
    }
    builder.finish()
}

async fn build_position_sidecar_from_btree_pages(
    config: &SandboxConfig,
    dataset: &Dataset,
    lookup_tier: Option<u8>,
) -> Result<PositionRowIdIndex> {
    Ok(
        build_position_sidecar_from_btree_pages_measured(config, dataset, lookup_tier)
            .await?
            .index,
    )
}

async fn build_position_sidecar_from_btree_pages_measured(
    config: &SandboxConfig,
    dataset: &Dataset,
    lookup_tier: Option<u8>,
) -> Result<BtreePageLoad> {
    let metadata_started = Instant::now();
    let index_segments = load_position_btree_segments(config, dataset).await?;
    let metadata_seconds = metadata_started.elapsed().as_secs_f64();

    let tier_started = Instant::now();
    let tier_rows = load_tier_row_filter(dataset, lookup_tier).await?;
    let tier_filter_seconds = tier_started.elapsed().as_secs_f64();

    let page_started = Instant::now();
    if index_segments.len() == 1 {
        let mut builder = PositionRowIdIndexBuilder::new();
        let stream_stats = append_btree_segment_to_builder(
            config,
            dataset,
            &index_segments[0],
            tier_rows.as_ref(),
            &mut builder,
        )
        .await?;
        return Ok(BtreePageLoad {
            index: builder.finish()?,
            metadata_seconds,
            tier_filter_seconds,
            page_scan_seconds: page_started.elapsed().as_secs_f64(),
            stream_batches: stream_stats.batches,
            stream_rows: stream_stats.rows,
        });
    }

    let mut pairs = Vec::new();
    let mut stream_stats = BtreePageStreamStats::default();
    for index_segment in &index_segments {
        stream_stats.add(
            append_btree_segment_pairs(
                config,
                dataset,
                index_segment,
                tier_rows.as_ref(),
                &mut pairs,
            )
            .await?,
        );
    }
    pairs.sort_unstable_by_key(|&(position, row_id)| (position, row_id));
    let mut builder = PositionRowIdIndexBuilder::new();
    for (position, row_id) in pairs {
        builder.push_monotonic(position, row_id)?;
    }
    Ok(BtreePageLoad {
        index: builder.finish()?,
        metadata_seconds,
        tier_filter_seconds,
        page_scan_seconds: page_started.elapsed().as_secs_f64(),
        stream_batches: stream_stats.batches,
        stream_rows: stream_stats.rows,
    })
}

async fn load_position_btree_segments(
    config: &SandboxConfig,
    dataset: &Dataset,
) -> Result<Vec<IndexMetadata>> {
    let index_name = format!("{}_btree_idx", config.key.column);
    let mut index_segments = dataset.load_indices_by_name(&index_name).await?;
    if index_segments.is_empty() {
        if let Some(index) = dataset
            .load_scalar_index(
                IndexCriteria::default()
                    .for_column(&config.key.column)
                    .supports_exact_equality(),
            )
            .await?
        {
            index_segments.push(index);
        }
    }
    if index_segments.is_empty() {
        bail!(
            "missing scalar BTree index '{index_name}' for column '{}'",
            config.key.column
        );
    }
    Ok(index_segments)
}

async fn load_tier_row_filter(
    dataset: &Dataset,
    lookup_tier: Option<u8>,
) -> Result<Option<RowAddrTreeMap>> {
    let Some(tier) = lookup_tier else {
        return Ok(None);
    };

    let metrics = LocalMetricsCollector::default();
    let tier_addrs = search_scalar_index(
        dataset,
        "tier",
        "tier_bitmap_idx",
        SargableQuery::Equals(tier_scalar_value(dataset, tier)?),
        &metrics,
    )
    .await?;
    Ok(Some(tier_addrs.true_rows()))
}

async fn append_btree_segment_to_builder(
    config: &SandboxConfig,
    dataset: &Dataset,
    index_segment: &IndexMetadata,
    tier_rows: Option<&RowAddrTreeMap>,
    builder: &mut PositionRowIdIndexBuilder,
) -> Result<BtreePageStreamStats> {
    let store = LanceIndexStore::from_dataset_for_existing(dataset, index_segment).await?;
    let reader = store.open_index_file(BTREE_PAGE_DATA_FILE).await?;
    let num_rows = reader.num_rows();
    if num_rows == 0 {
        return Ok(BtreePageStreamStats::default());
    }
    builder.reserve(btree_page_capacity_hint(num_rows, tier_rows));
    let mut stream = reader
        .read_range_stream(0..num_rows, Some(&[BTREE_VALUES_COLUMN, BTREE_IDS_COLUMN]))
        .await?;
    let mut stats = BtreePageStreamStats::default();
    while let Some(batch) = stream.try_next().await? {
        stats.batches += 1;
        stats.rows += batch.num_rows();
        append_btree_page_row_ids_to_builder(config, &batch, tier_rows, builder)?;
    }
    Ok(stats)
}

async fn append_btree_segment_pairs(
    config: &SandboxConfig,
    dataset: &Dataset,
    index_segment: &IndexMetadata,
    tier_rows: Option<&RowAddrTreeMap>,
    pairs: &mut Vec<(u32, u64)>,
) -> Result<BtreePageStreamStats> {
    let store = LanceIndexStore::from_dataset_for_existing(dataset, index_segment).await?;
    let reader = store.open_index_file(BTREE_PAGE_DATA_FILE).await?;
    let num_rows = reader.num_rows();
    if num_rows == 0 {
        return Ok(BtreePageStreamStats::default());
    }
    pairs.reserve(btree_page_capacity_hint(num_rows, tier_rows));
    let mut stream = reader
        .read_range_stream(0..num_rows, Some(&[BTREE_VALUES_COLUMN, BTREE_IDS_COLUMN]))
        .await?;
    let mut stats = BtreePageStreamStats::default();
    while let Some(batch) = stream.try_next().await? {
        stats.batches += 1;
        stats.rows += batch.num_rows();
        append_btree_page_row_id_pairs(config, &batch, tier_rows, pairs)?;
    }
    Ok(stats)
}

fn btree_page_capacity_hint(num_rows: usize, tier_rows: Option<&RowAddrTreeMap>) -> usize {
    tier_rows
        .and_then(RowSetOps::len)
        .and_then(|len| usize::try_from(len).ok())
        .unwrap_or(num_rows)
        .min(num_rows)
}

async fn sidecar_take_rows(
    config: &SandboxConfig,
    dataset: &Dataset,
    sidecar: &PositionRowIdIndex,
    positions: &[u64],
    projection: &[String],
    projection_label: &str,
    lookup_tier: Option<u8>,
) -> Result<ScanResult> {
    let started = Instant::now();
    let resolve_started = Instant::now();
    let resolved = sidecar.resolve_sorted_positions(positions);
    let row_id_resolve_seconds = resolve_started.elapsed().as_secs_f64();
    let (batch, io, take_rows_seconds) =
        take_rows_projected(dataset, &resolved.row_ids, projection).await?;
    let selected_positions = key_values_from_batch(config, &batch)?
        .into_iter()
        .collect::<HashSet<_>>()
        .len();

    Ok(ScanResult {
        name: format!(
            "{}_sidecar_take_{projection_label}",
            lookup_label(lookup_tier)
        ),
        filter_keys: positions.len(),
        lookup_batch_size: None,
        scans: 1,
        rows: batch.num_rows(),
        result_batches: usize::from(batch.num_rows() > 0),
        selected_positions,
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: Some(row_id_resolve_seconds),
        btree_resolve_seconds: None,
        take_rows_seconds: Some(take_rows_seconds),
        io,
    })
}

async fn sidecar_take_rows_batched(
    config: &SandboxConfig,
    dataset: &Dataset,
    sidecar: &PositionRowIdIndex,
    positions: &[u64],
    lookup_batch_size: usize,
    projection: &[String],
    projection_label: &str,
    lookup_tier: Option<u8>,
) -> Result<ScanResult> {
    let started = Instant::now();
    let mut rows = 0usize;
    let mut result_batches = 0usize;
    let mut selected_positions = HashSet::new();
    let mut io = ScanIoStats::default();
    let mut scans = 0usize;
    let mut filter_keys = 0usize;
    let mut row_id_resolve_seconds = 0.0f64;
    let mut take_rows_seconds = 0.0f64;
    let mut sidecar_cursor = 0usize;

    for chunk in positions.chunks(lookup_batch_size.max(1)) {
        scans += 1;
        filter_keys += chunk.len();

        let resolve_started = Instant::now();
        let resolved = sidecar.resolve_sorted_positions_from_cursor(chunk, &mut sidecar_cursor);
        row_id_resolve_seconds += resolve_started.elapsed().as_secs_f64();

        if resolved.row_ids.is_empty() {
            continue;
        }

        let (batch, take_io, take_seconds) =
            take_rows_projected(dataset, &resolved.row_ids, projection).await?;
        io.add(&take_io);
        take_rows_seconds += take_seconds;
        rows += batch.num_rows();
        result_batches += usize::from(batch.num_rows() > 0);
        selected_positions.extend(key_values_from_batch(config, &batch)?);
    }

    Ok(ScanResult {
        name: format!(
            "{}_sidecar_take_{projection_label}_batch_{lookup_batch_size}",
            lookup_label(lookup_tier)
        ),
        filter_keys,
        lookup_batch_size: Some(lookup_batch_size),
        scans,
        rows,
        result_batches,
        selected_positions: selected_positions.len(),
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: Some(row_id_resolve_seconds),
        btree_resolve_seconds: None,
        take_rows_seconds: Some(take_rows_seconds),
        io,
    })
}

async fn btree_direct_take_rows(
    config: &SandboxConfig,
    dataset: &Dataset,
    row_index: &PositionRowIdIndex,
    positions: &[u64],
    projection: &[String],
    projection_label: &str,
    lookup_tier: Option<u8>,
) -> Result<ScanResult> {
    let started = Instant::now();
    let resolve_started = Instant::now();
    let resolved = row_index.resolve_sorted_positions(positions);
    let row_ids = resolved.row_ids;
    let mut io = ScanIoStats::default();
    let row_id_resolve_seconds = resolve_started.elapsed().as_secs_f64();

    let (batch, take_io, take_rows_seconds) =
        take_rows_projected(dataset, &row_ids, projection).await?;
    io.add(&take_io);
    let selected_positions = key_values_from_batch(config, &batch)?
        .into_iter()
        .collect::<HashSet<_>>()
        .len();

    Ok(ScanResult {
        name: format!(
            "{}_btree_direct_take_{projection_label}",
            lookup_label(lookup_tier)
        ),
        filter_keys: positions.len(),
        lookup_batch_size: None,
        scans: 1,
        rows: batch.num_rows(),
        result_batches: usize::from(batch.num_rows() > 0),
        selected_positions,
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: Some(row_id_resolve_seconds),
        btree_resolve_seconds: Some(row_id_resolve_seconds),
        take_rows_seconds: Some(take_rows_seconds),
        io,
    })
}

async fn btree_direct_take_rows_batched(
    config: &SandboxConfig,
    dataset: &Dataset,
    row_index: &PositionRowIdIndex,
    positions: &[u64],
    lookup_batch_size: usize,
    projection: &[String],
    projection_label: &str,
    lookup_tier: Option<u8>,
) -> Result<ScanResult> {
    let started = Instant::now();
    let mut rows = 0usize;
    let mut result_batches = 0usize;
    let mut selected_positions = HashSet::new();
    let mut io = ScanIoStats::default();
    let mut scans = 0usize;
    let mut filter_keys = 0usize;
    let mut row_id_resolve_seconds = 0.0f64;
    let mut take_rows_seconds = 0.0f64;
    let mut row_index_cursor = 0usize;

    for chunk in positions.chunks(lookup_batch_size.max(1)) {
        scans += 1;
        filter_keys += chunk.len();

        let resolve_started = Instant::now();
        let resolved = row_index.resolve_sorted_positions_from_cursor(chunk, &mut row_index_cursor);
        let row_ids = resolved.row_ids;
        row_id_resolve_seconds += resolve_started.elapsed().as_secs_f64();

        if row_ids.is_empty() {
            continue;
        }

        let (batch, take_io, take_seconds) =
            take_rows_projected(dataset, &row_ids, projection).await?;
        io.add(&take_io);
        take_rows_seconds += take_seconds;
        rows += batch.num_rows();
        result_batches += usize::from(batch.num_rows() > 0);
        selected_positions.extend(key_values_from_batch(config, &batch)?);
    }

    Ok(ScanResult {
        name: format!(
            "{}_btree_direct_take_{projection_label}_batch_{lookup_batch_size}",
            lookup_label(lookup_tier)
        ),
        filter_keys,
        lookup_batch_size: Some(lookup_batch_size),
        scans,
        rows,
        result_batches,
        selected_positions: selected_positions.len(),
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: Some(row_id_resolve_seconds),
        btree_resolve_seconds: Some(row_id_resolve_seconds),
        take_rows_seconds: Some(take_rows_seconds),
        io,
    })
}

async fn search_scalar_index(
    dataset: &Dataset,
    column: &str,
    index_name: &str,
    query: SargableQuery,
    metrics: &LocalMetricsCollector,
) -> Result<NullableRowAddrSet> {
    let mut index_segments = dataset.load_indices_by_name(index_name).await?;
    if index_segments.is_empty() {
        if let Some(index) = dataset
            .load_scalar_index(
                IndexCriteria::default()
                    .for_column(column)
                    .supports_exact_equality(),
            )
            .await?
        {
            index_segments.push(index);
        }
    }
    if index_segments.is_empty() {
        bail!("missing scalar index '{index_name}' for column '{column}'");
    }

    let mut selections = Vec::with_capacity(index_segments.len());
    for index_segment in index_segments {
        let uuid = index_segment.uuid.to_string();
        let index = dataset.open_scalar_index(column, &uuid, metrics).await?;
        let result = index.search(&query, metrics).await?;
        match result {
            SearchResult::Exact(row_addrs) => selections.push(row_addrs),
            other => {
                bail!(
                    "scalar index '{index_name}' for column '{column}' returned non-exact result: {other:?}"
                );
            }
        }
    }

    Ok(NullableRowAddrSet::union_all(&selections))
}

fn tier_scalar_value(dataset: &Dataset, tier: u8) -> Result<ScalarValue> {
    let field = dataset
        .schema()
        .field("tier")
        .context("tier-scoped BTree direct lookup requires a tier field")?;
    match field.data_type() {
        DataType::UInt8 => Ok(ScalarValue::UInt8(Some(tier))),
        DataType::Int8 => Ok(ScalarValue::Int8(Some(
            i8::try_from(tier).context("tier value does not fit Int8")?,
        ))),
        other => bail!("tier bitmap direct lookup requires Int8 or UInt8 tier field, got {other}"),
    }
}

async fn take_rows_projected(
    dataset: &Dataset,
    row_ids: &[u64],
    projection: &[String],
) -> Result<(arrow_array::RecordBatch, ScanIoStats, f64)> {
    let projection_request =
        ProjectionRequest::from_columns(projection.iter().map(String::as_str), dataset.schema());
    let before = io_counter_snapshot();
    let started = Instant::now();
    let batch = dataset.take_rows(row_ids, projection_request).await?;
    let seconds = started.elapsed().as_secs_f64();
    let io = io_counter_delta(before);
    Ok((batch, io, seconds))
}

fn io_counter_snapshot() -> (u64, u64) {
    (iops_counter(), bytes_read_counter())
}

fn io_counter_delta(before: (u64, u64)) -> ScanIoStats {
    let after = io_counter_snapshot();
    ScanIoStats {
        iops: after.0.saturating_sub(before.0) as usize,
        bytes_read: after.1.saturating_sub(before.1) as usize,
        ..Default::default()
    }
}

fn row_ids_from_batch(batch: &arrow_array::RecordBatch) -> Result<Vec<u64>> {
    let row_ids = batch
        .column_by_name(ROW_ID)
        .context("row-id scan did not return _rowid column")?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .context("_rowid column must be UInt64")?;
    Ok(row_ids.values().to_vec())
}

fn append_position_row_ids_to_builder(
    config: &SandboxConfig,
    batch: &arrow_array::RecordBatch,
    builder: &mut PositionRowIdIndexBuilder,
) -> Result<()> {
    let mut pairs = Vec::with_capacity(batch.num_rows());
    append_position_row_id_pairs(config, batch, &mut pairs)?;
    for (position, row_id) in pairs {
        builder.push_monotonic(position, row_id)?;
    }
    Ok(())
}

fn append_btree_page_row_ids_to_builder(
    config: &SandboxConfig,
    batch: &arrow_array::RecordBatch,
    tier_rows: Option<&RowAddrTreeMap>,
    builder: &mut PositionRowIdIndexBuilder,
) -> Result<()> {
    let value_array = batch
        .column_by_name(BTREE_VALUES_COLUMN)
        .unwrap_or_else(|| batch.column(0));
    let row_ids = batch
        .column_by_name(BTREE_IDS_COLUMN)
        .or_else(|| batch.column_by_name(ROW_ID))
        .context("BTree page data did not return ids column")?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .context("BTree page ids column must be UInt64")?;

    match config.key.data_type {
        KeyDataType::Uint32 => {
            let positions = value_array
                .as_any()
                .downcast_ref::<UInt32Array>()
                .context("BTree page values column must be UInt32")?;
            for row in 0..batch.num_rows() {
                if positions.is_null(row) {
                    continue;
                }
                let row_id = row_ids.value(row);
                if tier_rows.is_none_or(|rows| rows.contains(row_id)) {
                    builder.push_monotonic(positions.value(row), row_id)?;
                }
            }
        }
        KeyDataType::Uint64 => {
            let positions = value_array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .context("BTree page values column must be UInt64")?;
            for row in 0..batch.num_rows() {
                if positions.is_null(row) {
                    continue;
                }
                let row_id = row_ids.value(row);
                if tier_rows.is_none_or(|rows| rows.contains(row_id)) {
                    let position = u32::try_from(positions.value(row))
                        .context("BTree page sidecar only supports UInt32-compatible positions")?;
                    builder.push_monotonic(position, row_id)?;
                }
            }
        }
    }
    Ok(())
}

fn append_btree_page_row_id_pairs(
    config: &SandboxConfig,
    batch: &arrow_array::RecordBatch,
    tier_rows: Option<&RowAddrTreeMap>,
    pairs: &mut Vec<(u32, u64)>,
) -> Result<()> {
    let value_array = batch
        .column_by_name(BTREE_VALUES_COLUMN)
        .unwrap_or_else(|| batch.column(0));
    let row_ids = batch
        .column_by_name(BTREE_IDS_COLUMN)
        .or_else(|| batch.column_by_name(ROW_ID))
        .context("BTree page data did not return ids column")?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .context("BTree page ids column must be UInt64")?;

    match config.key.data_type {
        KeyDataType::Uint32 => {
            let positions = value_array
                .as_any()
                .downcast_ref::<UInt32Array>()
                .context("BTree page values column must be UInt32")?;
            for row in 0..batch.num_rows() {
                if positions.is_null(row) {
                    continue;
                }
                let row_id = row_ids.value(row);
                if tier_rows.is_none_or(|rows| rows.contains(row_id)) {
                    pairs.push((positions.value(row), row_id));
                }
            }
        }
        KeyDataType::Uint64 => {
            let positions = value_array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .context("BTree page values column must be UInt64")?;
            for row in 0..batch.num_rows() {
                if positions.is_null(row) {
                    continue;
                }
                let row_id = row_ids.value(row);
                if tier_rows.is_none_or(|rows| rows.contains(row_id)) {
                    let position = u32::try_from(positions.value(row))
                        .context("BTree page sidecar only supports UInt32-compatible positions")?;
                    pairs.push((position, row_id));
                }
            }
        }
    }

    Ok(())
}

fn append_position_row_id_pairs(
    config: &SandboxConfig,
    batch: &arrow_array::RecordBatch,
    pairs: &mut Vec<(u32, u64)>,
) -> Result<()> {
    let position_idx = batch.schema().index_of(config.key.column.as_str())?;
    let position_array = batch.column(position_idx);
    let row_ids = batch
        .column_by_name(ROW_ID)
        .context("sidecar scan did not return _rowid column")?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .context("_rowid column must be UInt64")?;

    match config.key.data_type {
        KeyDataType::Uint32 => {
            let positions = position_array
                .as_any()
                .downcast_ref::<UInt32Array>()
                .context("position sidecar key column must be UInt32")?;
            for row in 0..batch.num_rows() {
                if !positions.is_null(row) {
                    pairs.push((positions.value(row), row_ids.value(row)));
                }
            }
        }
        KeyDataType::Uint64 => {
            let positions = position_array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .context("position sidecar key column must be UInt64")?;
            for row in 0..batch.num_rows() {
                if !positions.is_null(row) {
                    let position = u32::try_from(positions.value(row))
                        .context("position sidecar only supports UInt32-compatible positions")?;
                    pairs.push((position, row_ids.value(row)));
                }
            }
        }
    }

    Ok(())
}

async fn scan_once(
    config: &SandboxConfig,
    dataset: &Dataset,
    name: &str,
    filter: String,
    lookup_batch_size: Option<usize>,
    projection: &[String],
) -> Result<ScanResult> {
    let started = Instant::now();
    let stats_slot = Arc::new(Mutex::new(None::<ExecutionSummaryCounts>));
    let stats_sink = Arc::clone(&stats_slot);
    let projection_refs = projection.iter().map(String::as_str).collect::<Vec<_>>();
    let mut scanner = dataset.scan();
    scanner
        .filter(&filter)?
        .project(&projection_refs)?
        .materialization_style(MaterializationStyle::AllLate)
        .scan_stats_callback(Arc::new(move |stats| {
            *stats_sink.lock().expect("stats lock poisoned") = Some(stats.clone());
        }));
    let mut stream = scanner.try_into_stream().await?;
    let mut rows = 0usize;
    let mut result_batches = 0usize;
    let mut selected_positions = HashSet::new();
    while let Some(batch) = stream.try_next().await? {
        rows += batch.num_rows();
        result_batches += 1;
        selected_positions.extend(key_values_from_batch(config, &batch)?);
    }
    let mut io = ScanIoStats::default();
    if let Some(stats) = stats_slot.lock().expect("stats lock poisoned").as_ref() {
        io.add_execution(stats);
    }
    Ok(ScanResult {
        name: name.to_string(),
        filter_keys: 0,
        lookup_batch_size,
        scans: 1,
        rows,
        result_batches,
        selected_positions: selected_positions.len(),
        seconds: started.elapsed().as_secs_f64(),
        row_id_resolve_seconds: None,
        btree_resolve_seconds: None,
        take_rows_seconds: None,
        io,
    })
}

async fn capture_physical_plans(
    config: &SandboxConfig,
    dataset: &Dataset,
    positions: &[u64],
    projection: &[String],
    lookup_tier: Option<u8>,
) -> Result<Vec<PhysicalPlanReport>> {
    let mut plans = Vec::new();
    plans.push(
        capture_scan_plan(
            dataset,
            "warm_full_scan_all_everything_columns",
            "tier = 0".to_string(),
            "tier = 0".to_string(),
            None,
            projection,
            None,
        )
        .await?,
    );

    if positions.is_empty() {
        return Ok(plans);
    }

    if let Some(&lookup_batch_size) = config.benchmark.lookup_batch_sizes.first() {
        let chunk = representative_position_chunk(positions, lookup_batch_size);
        let key_only_projection = vec![config.key.column.clone()];
        plans.push(
            capture_scan_plan(
                dataset,
                format!("{}_index_prewarm_key_only", lookup_label(lookup_tier)),
                position_filter(config, lookup_tier, chunk),
                position_filter_summary(config, lookup_tier, chunk),
                Some(lookup_batch_size),
                &key_only_projection,
                Some((&config.key.column, chunk)),
            )
            .await?,
        );
    }

    for &lookup_batch_size in &config.benchmark.lookup_batch_sizes {
        let chunk = representative_position_chunk(positions, lookup_batch_size);
        plans.push(
            capture_scan_plan(
                dataset,
                format!(
                    "{}_lookup_batch_{lookup_batch_size}",
                    lookup_label(lookup_tier)
                ),
                position_filter(config, lookup_tier, chunk),
                position_filter_summary(config, lookup_tier, chunk),
                Some(lookup_batch_size),
                projection,
                Some((&config.key.column, chunk)),
            )
            .await?,
        );
    }

    Ok(plans)
}

fn representative_position_chunk(positions: &[u64], lookup_batch_size: usize) -> &[u64] {
    &positions[..positions.len().min(lookup_batch_size.max(1))]
}

async fn capture_scan_plan(
    dataset: &Dataset,
    name: impl Into<String>,
    filter: String,
    display_filter: String,
    lookup_batch_size: Option<usize>,
    projection: &[String],
    lookup_keys: Option<(&str, &[u64])>,
) -> Result<PhysicalPlanReport> {
    let projection_refs = projection.iter().map(String::as_str).collect::<Vec<_>>();
    let mut scanner = dataset.scan();
    scanner
        .filter(&filter)?
        .project(&projection_refs)?
        .materialization_style(MaterializationStyle::AllLate);
    let plan = scanner.create_plan().await?;
    let plan = format!("{}", displayable(plan.as_ref()).indent(true));
    let plan = lookup_keys
        .map(|(column, keys)| compact_point_lookup_plan_text(&plan, column, keys))
        .unwrap_or(plan);
    Ok(PhysicalPlanReport {
        name: name.into(),
        filter: display_filter,
        lookup_batch_size,
        projected_column_count: projection.len(),
        projection: projection.to_vec(),
        plan,
    })
}

impl ScanIoStats {
    fn add(&mut self, other: &ScanIoStats) {
        self.iops += other.iops;
        self.requests += other.requests;
        self.bytes_read += other.bytes_read;
        self.indices_loaded += other.indices_loaded;
        self.parts_loaded += other.parts_loaded;
        self.index_comparisons += other.index_comparisons;
        self.fragments_scanned += other.fragments_scanned;
        self.ranges_scanned += other.ranges_scanned;
        self.rows_scanned += other.rows_scanned;
    }

    fn add_execution(&mut self, counts: &ExecutionSummaryCounts) {
        self.iops += counts.iops;
        self.requests += counts.requests;
        self.bytes_read += counts.bytes_read;
        self.indices_loaded += counts.indices_loaded;
        self.parts_loaded += counts.parts_loaded;
        self.index_comparisons += counts.index_comparisons;
        self.fragments_scanned += counts
            .all_counts
            .get("fragments_scanned")
            .copied()
            .unwrap_or_default();
        self.ranges_scanned += counts
            .all_counts
            .get("ranges_scanned")
            .copied()
            .unwrap_or_default();
        self.rows_scanned += counts
            .all_counts
            .get("rows_scanned")
            .copied()
            .unwrap_or_default();
    }
}

fn lookup_scope_label(lookup_tier: Option<u8>) -> String {
    match lookup_tier {
        Some(tier) => format!("tier={tier}"),
        None => "all tiers".to_string(),
    }
}

fn lookup_label(lookup_tier: Option<u8>) -> String {
    match lookup_tier {
        Some(1) => "cold".to_string(),
        Some(tier) => format!("tier{tier}"),
        None => "alltiers".to_string(),
    }
}

fn position_filter(config: &SandboxConfig, lookup_tier: Option<u8>, keys: &[u64]) -> String {
    let mut filter = match lookup_tier {
        Some(tier) => format!("tier = {tier} AND {} IN (", config.key.column),
        None => format!("{} IN (", config.key.column),
    };
    for (idx, key) in keys.iter().enumerate() {
        if idx > 0 {
            filter.push_str(", ");
        }
        filter.push_str(&key.to_string());
    }
    filter.push(')');
    filter
}

fn position_filter_summary(
    config: &SandboxConfig,
    lookup_tier: Option<u8>,
    keys: &[u64],
) -> String {
    match lookup_tier {
        Some(tier) => format!(
            "tier = {tier} AND {} IN ({})",
            config.key.column,
            lookup_key_summary(keys)
        ),
        None => format!("{} IN ({})", config.key.column, lookup_key_summary(keys)),
    }
}

fn lookup_key_summary(keys: &[u64]) -> String {
    match (keys.first(), keys.last()) {
        (Some(first), Some(last)) => {
            format!("<{} keys; first={first}; last={last}>", keys.len())
        }
        _ => "<0 keys>".to_string(),
    }
}

fn compact_point_lookup_plan_text(plan: &str, column: &str, keys: &[u64]) -> String {
    let summary = lookup_key_summary(keys);
    let plan = compact_typed_in_list(plan, column, &summary);
    compact_scalar_index_in_list(&plan, column, &summary)
}

fn compact_typed_in_list(plan: &str, column: &str, summary: &str) -> String {
    let needle = format!("{column} IN ([");
    let replacement = format!("{column} IN ([{summary}])");
    let mut compacted = String::with_capacity(plan.len().min(4096));
    let mut rest = plan;

    while let Some(start) = rest.find(&needle) {
        compacted.push_str(&rest[..start]);
        let after_needle = &rest[start + needle.len()..];
        if let Some(end) = after_needle.find("])") {
            compacted.push_str(&replacement);
            rest = &after_needle[end + 2..];
        } else {
            compacted.push_str(&rest[start..]);
            return compacted;
        }
    }

    compacted.push_str(rest);
    compacted
}

fn compact_scalar_index_in_list(plan: &str, column: &str, summary: &str) -> String {
    let needle = format!("{column} IN [");
    let replacement = format!("{column} IN [{summary}]");
    let mut compacted = String::with_capacity(plan.len().min(4096));
    let mut rest = plan;

    while let Some(start) = rest.find(&needle) {
        compacted.push_str(&rest[..start]);
        let after_needle = &rest[start + needle.len()..];
        if let Some(end) = after_needle.find("]@") {
            compacted.push_str(&replacement);
            rest = &after_needle[end..];
        } else {
            compacted.push_str(&rest[start..]);
            return compacted;
        }
    }

    compacted.push_str(rest);
    compacted
}

fn read_positions(path: &Path) -> Result<Vec<u64>> {
    let raw = fs::read_to_string(path)
        .with_context(|| format!("failed to read positions file '{}'", path.display()))?;
    raw.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| Ok(line.parse::<u64>()?))
        .collect()
}

fn scan_progress_summary(result: &ScanResult) -> String {
    format!(
        "{}: rows={}, selected_positions={}, seconds={:.3}, bytes_read={}, iops={}, requests={}, ranges={}",
        result.name,
        result.rows,
        result.selected_positions,
        result.seconds,
        result.io.bytes_read,
        result.io.iops,
        result.io.requests,
        result.io.ranges_scanned
    )
}

#[allow(dead_code)]
fn _duration_seconds(duration: Duration) -> f64 {
    duration.as_secs_f64()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use arrow_array::{ArrayRef, RecordBatch, UInt32Array, UInt64Array};
    use arrow_schema::{DataType, Field, Schema};

    use crate::config::{
        BenchmarkConfig, DatasetConfig, DefaultsConfig, EncodingConfig, IndexConfig, InspectConfig,
        KeyConfig, KeyDataType, KeyMode, SandboxConfig, SidecarConfig,
    };

    fn test_config() -> SandboxConfig {
        SandboxConfig {
            dataset: DatasetConfig {
                name: "test".to_string(),
                cache_root: ".".into(),
                chrom: "chr1".to_string(),
                output_root: ".".into(),
                lance_version: "2.1".to_string(),
                batch_size: 1024,
                overwrite: true,
            },
            key: KeyConfig {
                mode: KeyMode::Position,
                column: "position".to_string(),
                data_type: KeyDataType::Uint32,
            },
            indexes: IndexConfig {
                position_btree: true,
                tier_bitmap: true,
            },
            benchmark: BenchmarkConfig {
                positions_file: "positions.txt".into(),
                vcf_path: "input.vcf.gz".into(),
                lookup_batch_sizes: vec![512],
                sort_position_keys: true,
            },
            sidecar: SidecarConfig::default(),
            inspect: InspectConfig {
                include_physical_columns: true,
                include_index_sizes: true,
            },
            defaults: DefaultsConfig {
                encoding: EncodingConfig {
                    structural: "miniblock".to_string(),
                    compression: "zstd".to_string(),
                    compression_level: 3,
                    dict_values_compression: "zstd".to_string(),
                    dict_values_compression_level: 3,
                    rle_threshold: 0.95,
                    dict_size_ratio: 0.99,
                    dict_divisor: 1,
                    minichunk_size: 4096,
                },
            },
            fields: Default::default(),
            structs: Default::default(),
        }
    }

    #[test]
    fn only_index_take_everything_mode_skips_non_index_paths() {
        let mode = super::BenchmarkMode::OnlyIndexTakeEverything;

        assert!(!mode.run_warm_scan());
        assert!(!mode.run_cold_scanner_paths());
        assert!(!mode.run_unbatched_btree_take());
        assert!(!mode.run_position_only_take());
        assert!(!mode.capture_physical_plans());

        assert!(super::BenchmarkMode::Full.run_unbatched_btree_take());
    }

    #[test]
    fn compacts_large_point_lookup_in_lists_in_display_plans() {
        let keys = (1_000_u64..3_048_u64).collect::<Vec<_>>();
        let literals = keys
            .iter()
            .map(|key| format!("UInt32({key})"))
            .collect::<Vec<_>>()
            .join(", ");
        let raw_plan = format!(
            "LanceRead: full_filter=tier = Int8(1) AND position IN ([{literals}]), row_id=false\n\
             ScalarIndexQuery: query=AND([tier = 1]@tier_bitmap_idx(Bitmap),[position IN [{}]]@position_btree_idx(BTree))",
            keys.iter()
                .map(u64::to_string)
                .collect::<Vec<_>>()
                .join(",")
        );

        let compacted = super::compact_point_lookup_plan_text(&raw_plan, "position", &keys);

        assert!(compacted.contains("position IN ([<2048 keys; first=1000; last=3047>])"));
        assert!(compacted.contains("position IN [<2048 keys; first=1000; last=3047>]"));
        assert!(!compacted.contains("UInt32(1000), UInt32(1001), UInt32(1002)"));
        assert!(!compacted.contains("1000,1001,1002"));
        assert!(compacted.lines().all(|line| line.len() < 1_000));
    }

    #[test]
    fn btree_page_pairs_extract_values_and_ids_columns() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("values", DataType::UInt32, true),
            Field::new("ids", DataType::UInt64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(vec![Some(20), None, Some(10)])) as ArrayRef,
                Arc::new(UInt64Array::from(vec![200, 201, 100])) as ArrayRef,
            ],
        )
        .unwrap();
        let mut pairs = Vec::new();

        super::append_btree_page_row_id_pairs(&test_config(), &batch, None, &mut pairs).unwrap();

        assert_eq!(pairs, vec![(20, 200), (10, 100)]);
    }
}
