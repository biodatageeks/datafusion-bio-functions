//! End-to-end per-chrom plugin cache build: register sources → ingest view →
//! normalization wrapper → variation-frequency join/tier → direct tiered
//! shard write → cache-manifest chrom entry.
//!
//! Both joins run in a bounded, parallel DataFusion context. For each join we
//! first inspect DataFusion's optimized HashJoin build side and use its actual
//! statistics to decide whether the estimated reservation fits the active
//! pool. If it does not, or the statistics are incomplete, only
//! `prefer_hash_join` is changed and DataFusion plans its built-in
//! SortMergeJoin. The final query always has an explicit
//! total `ORDER BY tier, start, ...`, so correctness does not depend on execution-
//! partition or hash-table order and the result can be written directly in
//! the point-lookup shard's physical order.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, Int8Array, UInt32Array};
use datafusion::arrow::compute::cast;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::dataframe::DataFrame;
use datafusion::datasource::MemTable;
use datafusion::execution::memory_pool::{FairSpillPool, MemoryPool};
use datafusion::execution::runtime_env::RuntimeEnvBuilder;
use datafusion::physical_plan::{SendableRecordBatchStream, stream::RecordBatchStreamAdapter};
use datafusion::prelude::{SessionConfig, SessionContext};
use futures::StreamExt;
use log::info;
use std::time::Instant;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::ChromEntry;
use crate::plugin_cache::dedup::{check_assume_unique_sample, dedup_keep_first};
use crate::plugin_cache::join::{
    partition_batches, should_retry_final_hash_plan, tiered_stream_sorted_adaptive,
    tiered_stream_sorted_sort_merge,
};
use crate::plugin_cache::mem_trace;
use crate::plugin_cache::normalize::{
    canonical_contig_str, canonical_contig_udf, wrap_normalization,
};
use crate::plugin_cache::provider::register_sources_for_chrom;
use crate::plugin_cache::source_manifest::SourceManifest;
use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};

/// Reproject `batch` to `schema`'s column order and types (casting where needed,
/// e.g. the normalized view's `Int64` `start`/`end` → the shard's `UInt32`).
fn reproject_cast(batch: &RecordBatch, schema: &SchemaRef) -> Result<RecordBatch> {
    let cols = schema
        .fields()
        .iter()
        .map(|f| {
            let idx = batch.schema().index_of(f.name())?;
            cast(batch.column(idx), f.data_type())
                .map_err(|e| DataFusionError::Execution(format!("cast {}: {e}", f.name())))
        })
        .collect::<Result<Vec<_>>>()?;
    RecordBatch::try_new(Arc::clone(schema), cols)
        .map_err(|e| DataFusionError::Execution(format!("reproject: {e}")))
}

#[derive(Debug, PartialEq, Eq)]
struct TierOrderSummary {
    last_key: Option<(i8, u32)>,
    warm: usize,
    cold: usize,
}

/// Verify `(tier, start)` is lexicographically non-decreasing within `batch`
/// and against the previous batch's final key. Return the batch's final key
/// and its warm/cold row counts.
///
/// The query's explicit total `ORDER BY tier, start, ...` provides this prefix
/// contract. Check
/// it at the storage boundary anyway: a planner/executor regression must abort
/// the write instead of producing a page directory whose binary-search lookup
/// can silently miss rows (`PageDir::resolve_ranges` assumes a sorted run).
fn inspect_tier_start_order(
    batch: &RecordBatch,
    start_idx: usize,
    tier_idx: usize,
    chrom: &str,
    plugin_name: &str,
    last_seen: Option<(i8, u32)>,
) -> Result<TierOrderSummary> {
    let starts = batch
        .column(start_idx)
        .as_any()
        .downcast_ref::<UInt32Array>()
        .ok_or_else(|| DataFusionError::Execution("start column not UInt32".into()))?;
    let tiers = batch
        .column(tier_idx)
        .as_any()
        .downcast_ref::<Int8Array>()
        .ok_or_else(|| DataFusionError::Execution("tier column not Int8".into()))?;
    let mut prev = last_seen;
    let (mut warm, mut cold) = (0usize, 0usize);
    for i in 0..starts.len() {
        if starts.is_null(i) || tiers.is_null(i) {
            return Err(DataFusionError::Execution(format!(
                "plugin_cache[{plugin_name}/{chrom}]: NULL tier/start in final shard row"
            )));
        }
        let key = (tiers.value(i), starts.value(i));
        match key.0 {
            0 => warm += 1,
            1 => cold += 1,
            tier => {
                return Err(DataFusionError::Execution(format!(
                    "plugin_cache[{plugin_name}/{chrom}]: invalid tier {tier}; expected 0 or 1"
                )));
            }
        }
        if let Some(previous) = prev
            && key < previous
        {
            return Err(DataFusionError::Execution(format!(
                "plugin_cache[{plugin_name}/{chrom}]: final shard write is not ordered by \
                 (tier, start): ({}, {}) follows ({}, {}) -- the tier join did not preserve \
                 row order; the \
                 on-disk point-lookup directory requires a sorted shard, so refusing to \
                 write a corrupt one. The adaptive plan's explicit ORDER BY contract was \
                 violated.",
                key.0, key.1, previous.0, previous.1
            )));
        }
        prev = Some(key);
    }
    Ok(TierOrderSummary {
        last_key: prev,
        warm,
        cold,
    })
}

/// Default memory budget for the adaptive tier plan's `FairSpillPool`, in MiB.
///
/// The pool bounds the WHOLE plan, not just its sort -- and a selected
/// `HashJoinExec` cannot spill (DataFusion has no spill path for a hash join's
/// build side), so for the join this budget is a hard ceiling rather than a
/// spill threshold. Only the sort can actually honour it by spilling.
///
/// The previous 2 GiB default was too small for that reality: measured on
/// chr21 (the smallest autosome, 64 GB machine), the join alone reserved
/// 1557 MB for dbNSFP, exhausting the pool before the sort reserved anything
/// and failing the build outright. 8 GiB clears it with room to spare (peak
/// RSS 1.7 GB, builds in 0.2 min) while still bounding a runaway plan.
///
/// CADD also fits this budget once `materialize_probe` copies its build input
/// into execution-sized owned batches. Before that structural fix, DataFusion
/// 53 repeatedly charged shared backing buffers and exhausted every tested
/// ceiling while RSS stayed flat; raising this default only moved the failure.
const DEFAULT_RETRY_SORT_MEMORY_MIB: usize = 8 * 1024;

/// Env override for [`DEFAULT_RETRY_SORT_MEMORY_MIB`], in MiB. The existing
/// variable name is retained for compatibility even though the bounded pool
/// now covers every tier plan rather than only an order-recovery retry.
///
/// A fixed default cannot suit both a 16 GB laptop and a 512 GB node, and the
/// right value depends on the plugin's row width as much as its row count --
/// so make it tunable rather than guessing. Invalid or zero values fall back
/// to the default rather than failing a long build on a typo'd env var.
const RETRY_SORT_MEMORY_ENV: &str = "VEPYR_PLUGIN_RETRY_MEMORY_MIB";
const TARGET_PARTITIONS_ENV: &str = "VEPYR_PLUGIN_TARGET_PARTITIONS";

fn retry_sort_memory_limit_bytes() -> usize {
    let mib = std::env::var(RETRY_SORT_MEMORY_ENV)
        .ok()
        .and_then(|raw| raw.trim().parse::<usize>().ok())
        .filter(|mib| *mib > 0)
        .unwrap_or(DEFAULT_RETRY_SORT_MEMORY_MIB);
    mib * 1024 * 1024
}

fn target_partitions_from(raw: Option<&str>, runtime_default: usize) -> usize {
    raw.and_then(|value| value.trim().parse::<usize>().ok())
        .filter(|partitions| *partitions > 0)
        .unwrap_or_else(|| runtime_default.max(2))
}

/// Session config for the tier-join passes.
///
/// `schema_force_view_types` is DataFusion's default-on parquet behaviour of
/// reading `Utf8` columns as `Utf8View`. It is off here for a specific reason:
/// a `StringViewArray`'s views point into data buffers that are SHARED across
/// the batches a scan emits, and `HashJoinExec` reserves memory per batch via
/// `get_record_batch_memory_size`, which de-duplicates buffers only WITHIN one
/// batch. So every batch re-counts the whole shared buffer, and the build-side
/// reservation grows by (batch count) x (buffer size) while only one copy
/// exists. Measured on CADD chr21: a ~600 MB build side reserved 6.95 GB at an
/// 8 GiB pool and 23.15 GB at 24 GiB -- always the same 592.7 MB increment,
/// once per batch -- while process RSS stayed flat near 16 GB. Reading plain
/// `Utf8` gives each batch its own buffers, so the accounting matches reality.
/// It also drops the `CAST(... AS Utf8View)` nodes the mismatch forced into
/// the join keys.
fn tier_join_config() -> SessionConfig {
    let config = SessionConfig::new();
    // Production follows DataFusion's runtime default with a two-worker
    // minimum. The override exists for controlled scaling measurements and
    // operator debugging; invalid/zero values preserve the automatic default.
    let target_partitions = target_partitions_from(
        std::env::var(TARGET_PARTITIONS_ENV).ok().as_deref(),
        config.target_partitions(),
    );
    config
        .with_target_partitions(target_partitions)
        .set_bool(
            "datafusion.execution.parquet.schema_force_view_types",
            false,
        )
        // Build a CollectLeft HashJoin candidate so the adaptive layer can
        // inspect the exact build side DataFusion selected. It later flips
        // only `prefer_hash_join` when that candidate does not fit, leaving
        // repartitioning, sorting, spilling, and join execution to DataFusion.
        .set_bool("datafusion.optimizer.prefer_hash_join", true)
        .set_usize(
            "datafusion.optimizer.hash_join_single_partition_threshold",
            usize::MAX,
        )
        .set_usize(
            "datafusion.optimizer.hash_join_single_partition_threshold_rows",
            usize::MAX,
        )
}

/// Source-ingest configuration with the same automatic/user-selected
/// parallelism as the downstream joins.
///
/// File-scan repartitioning is left enabled so plain CSV/TSV and Parquet can
/// split large files into ordered byte/row-group ranges. Round-robin
/// repartitioning is disabled because it destroys source order, which the
/// keep-first dedup contract needs. Providers that cannot split an input (for
/// example an unindexed VCF or the current BED provider) still correctly
/// expose one physical partition; no caller-side `target_partitions = 1`
/// restriction is imposed on providers that can parallelize.
fn source_read_config(target_partitions: usize) -> SessionConfig {
    SessionConfig::new()
        .with_target_partitions(target_partitions)
        .with_repartition_file_scans(true)
        .set_bool("datafusion.optimizer.enable_round_robin_repartition", false)
}

/// Execute all normalized source partitions concurrently, then replay their
/// batches in physical partition order.
///
/// `DataFrame::execute_stream()` inserts `CoalescePartitionsExec`, which emits
/// whichever partition produces a batch first and therefore loses source-file
/// order. `collect_partitioned()` also runs every partition concurrently, but
/// DataFusion explicitly sorts the completed results by partition index. File
/// scans and the indexed genomic partition balancer assign ranges in source
/// order, so flattening that outer vector restores the deterministic order
/// required by `dedup_keep_first` without serializing parsing.
async fn ordered_parallel_source_stream(
    df: DataFrame,
) -> Result<(SendableRecordBatchStream, usize)> {
    let schema: SchemaRef = df.schema().inner().clone();
    let partitions = df.collect_partitioned().await?;
    let partition_count = partitions.len();
    let batches = partitions.into_iter().flatten().map(Ok);
    Ok((
        Box::pin(RecordBatchStreamAdapter::new(
            schema,
            futures::stream::iter(batches),
        )),
        partition_count,
    ))
}

/// Deletes the build's scratch shards unless the build reached the point of
/// consuming them.
///
/// `PluginShardWriter::create` truncates, so a leftover `.tmp` is harmless to
/// the *next* build -- but a build that dies partway (a full disk during
/// `writer.write`, or a failed Hash-to-SMJ retry) would otherwise
/// leave up to two chromosome-sized scratch files behind, on the same disk
/// whose exhaustion caused the failure. Every early return from
/// `build_plugin_chrom` past the point the paths are chosen runs this.
struct ScratchGuard(Vec<PathBuf>);

impl ScratchGuard {
    fn new(paths: impl IntoIterator<Item = PathBuf>) -> Self {
        Self(paths.into_iter().collect())
    }

    /// Stop watching: the scratch files have been consumed into the final
    /// shard (or deliberately removed), so dropping must not delete anything.
    fn disarm(mut self) {
        self.0.clear();
    }
}

impl Drop for ScratchGuard {
    fn drop(&mut self) {
        for path in &self.0 {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// Consume the globally, totally sorted stream (whose leading keys are
/// `(tier, start)`) and write it directly to the final temporary shard. The
/// monotonicity check is a final assertion of the adaptive plan's explicit
/// sorted-output contract.
///
/// `PluginShardWriter::create` truncates its target path, so a runtime
/// Hash-to-SMJ fallback safely overwrites any partial first attempt.
async fn write_tiered_shard(
    mut stream: SendableRecordBatchStream,
    out_schema: &SchemaRef,
    path: &Path,
    chrom: &str,
    plugin_name: &str,
) -> Result<(usize, usize, usize)> {
    let mut writer = PluginShardWriter::create(path, Arc::clone(out_schema))?;
    let (mut warm, mut cold) = (0usize, 0usize);
    let start_idx = out_schema.index_of("start")?;
    let tier_idx = out_schema.index_of("tier")?;
    let mut last_key: Option<(i8, u32)> = None;

    while let Some(b) = stream.next().await {
        let batch = b?;
        if batch.num_rows() == 0 {
            continue;
        }
        let reordered = reproject_cast(&batch, out_schema)?;
        let summary = inspect_tier_start_order(
            &reordered,
            start_idx,
            tier_idx,
            chrom,
            plugin_name,
            last_key,
        )?;
        last_key = summary.last_key;
        warm += summary.warm;
        cold += summary.cold;
        writer.write(&reordered)?;
    }
    let rows = writer.finish()?;
    debug_assert_eq!(rows, warm + cold);
    Ok((rows, warm, cold))
}

/// Build one chromosome's plugin shard and make it live. Returns the
/// cache-manifest chrom entry.
pub async fn build_plugin_chrom(
    src: &SourceManifest,
    source_manifest_file: &str,
    variation_shard: &Path,
    output_cache_root: &Path,
    chrom: &str,
) -> Result<ChromEntry> {
    build_plugin_chrom_staged(
        src,
        source_manifest_file,
        variation_shard,
        output_cache_root,
        chrom,
    )
    .await?
    .commit()
}

/// A chromosome shard built to its staging file but not yet live. The
/// previous shard, if any, is untouched until [`StagedShard::commit`]; dropping
/// an uncommitted shard removes the staging file. The builder stages every
/// chromosome of a build and commits them together, so a build that fails on a
/// later chromosome leaves the cache exactly as it found it.
#[derive(Debug)]
pub struct StagedShard {
    pub entry: ChromEntry,
    /// The staging file, or `None` for an empty chromosome (whose commit
    /// removes any previous shard instead).
    staged: Option<PathBuf>,
    live: PathBuf,
}

impl StagedShard {
    /// Replace the live shard with the staged one (or remove it, for an empty
    /// chromosome) and return the manifest entry.
    pub fn commit(mut self) -> Result<ChromEntry> {
        match self.staged.take() {
            Some(staged) => std::fs::rename(&staged, &self.live).map_err(|e| {
                DataFusionError::Execution(format!(
                    "rename {} -> {}: {e}",
                    staged.display(),
                    self.live.display()
                ))
            })?,
            None => {
                let _ = std::fs::remove_file(&self.live);
            }
        }
        Ok(self.entry.clone())
    }
}

impl Drop for StagedShard {
    fn drop(&mut self) {
        if let Some(staged) = self.staged.take() {
            let _ = std::fs::remove_file(staged);
        }
    }
}

/// [`build_plugin_chrom`] up to, but not including, the point where the
/// live shard is touched.
pub async fn build_plugin_chrom_staged(
    src: &SourceManifest,
    _source_manifest_file: &str,
    variation_shard: &Path,
    output_cache_root: &Path,
    chrom: &str,
) -> Result<StagedShard> {
    src.validate()?;
    // Coarse-grained stage timing at `info` level — cheap (a handful of
    // `Instant::now()` calls) and the only thing that would have turned CADD
    // chr6's multi-hour "is it stuck?" investigation (no visibility into
    // which stage was running) into a 30-second log check.
    let t_start = Instant::now();
    let build_config = tier_join_config();
    // Source parsing uses the same adaptive/user-selected parallelism as the
    // downstream joins. Batches are collected concurrently per physical
    // partition and replayed in source-range order below, preserving VEP's
    // first-in-file duplicate semantics without a blanket serial scan.
    let read_ctx =
        SessionContext::new_with_config(source_read_config(build_config.target_partitions()));
    read_ctx.register_udf(canonical_contig_udf());
    // Held until this chrom's stream is fully materialized (dedup) below, then
    // dropped (deleting the decompressed temp) — so build_all keeps at most one temp.
    let src_temps = register_sources_for_chrom(&read_ctx, src, chrom).await?;

    // Ingest view (raw column mapping), then normalized view (contig + coords).
    read_ctx
        .sql(&format!(
            "CREATE OR REPLACE VIEW {} AS {}",
            src.ingest_view_name(),
            src.ingest_sql
        ))
        .await?;
    let value_cols: Vec<String> = src.value_columns.iter().map(|v| v.column.clone()).collect();
    let match_cols = src.match_column_names();
    let norm_sql = wrap_normalization(
        &src.ingest_view_name(),
        src.coordinate_system.clone(),
        &match_cols,
        &value_cols,
    );
    let norm_view = format!("plugin_{}_norm", src.plugin_name);
    // Filter on the same canonicalized contig the normalization applies to the
    // data (`canonical_contig_str` folds `M`/`chrM`/`chrMT` → `MT` and uppercases),
    // so `--chrom M`/`chrM`/`MT` all select the MT rows rather than silently
    // producing a 0-row shard.
    let source_chrom = canonical_contig_str(chrom);
    // Defense-in-depth: `chrom` is a trusted local arg, but reject anything
    // outside a safe contig charset before interpolating it into SQL.
    if source_chrom.is_empty()
        || !source_chrom
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
    {
        return Err(DataFusionError::Execution(format!(
            "invalid contig '{chrom}'"
        )));
    }
    read_ctx
        .sql(&format!(
            "CREATE OR REPLACE VIEW {norm_view} AS SELECT * FROM ({norm_sql}) WHERE chrom = '{source_chrom}'"
        ))
        .await?;

    // Collapse duplicate probe keys to their first source-file occurrence BEFORE
    // the tier join (which reorders rows and would destroy the file-order
    // tiebreak). Two overlapping genes can map the same genomic variant + aa-change
    // to different scores; VEP takes the first-in-file record, so we must too. This
    // executes the source partitions concurrently, restores physical
    // partition order, and keeps the first row per
    // `(start, allele_string, <match cols>)`.
    let norm_df = read_ctx.sql(&format!("SELECT * FROM {norm_view}")).await?;
    let (norm_stream, source_partitions) = ordered_parallel_source_stream(norm_df).await?;
    let norm_schema = norm_stream.schema();
    // `assume_unique` sources are claimed to never repeat a probe key, so the
    // exhaustive keep-first pass (a HashSet<String> with one entry per row —
    // the dominant memory cost on the largest chromosomes) is skipped in
    // favor of a bounded-memory sampled check that still catches a violated
    // claim instead of trusting it blindly (see `check_assume_unique_sample`
    // docs for why this can't be exhaustive without reintroducing the same
    // memory cost the flag exists to avoid).
    let deduped = if src.assume_unique {
        check_assume_unique_sample(norm_stream, &match_cols).await?
    } else {
        dedup_keep_first(norm_stream, &match_cols).await?
    };
    // The source scan is done — drop the decompressed temp before the join leg.
    drop(src_temps);
    let read_rows: usize = deduped.iter().map(|b| b.num_rows()).sum();
    info!(
        "plugin_cache[{}/{chrom}]: read+normalize+dedup done, {read_rows} rows, \
         source_partitions={source_partitions}, {:.1}s elapsed",
        src.plugin_name,
        t_start.elapsed().as_secs_f64()
    );

    // The join and ORDER BY share one bounded spill pool. TracingPool is always
    // installed because it also attributes a rejected reservation to the
    // operator that requested it; detailed peak logging remains opt-in.
    mem_trace::describe_batches(
        "plugin dedup output (join build side when not swapped)",
        &deduped,
    );
    let base_pool: Arc<dyn MemoryPool> =
        Arc::new(FairSpillPool::new(retry_sort_memory_limit_bytes()));
    let tracer = Arc::new(mem_trace::TracingPool::new(Arc::clone(&base_pool)));
    let pool: Arc<dyn MemoryPool> = Arc::clone(&tracer) as _;
    let build_runtime = RuntimeEnvBuilder::new()
        .with_memory_pool(Arc::clone(&pool))
        .build_arc()
        .map_err(|e| DataFusionError::Execution(format!("tier runtime env: {e}")))?;
    // The variation-filter semi join needs only these two columns. Register a
    // narrow zero-copy projection as its own MemTable so DataFusion sees exact
    // key-only statistics and never considers the plugin's wide value columns
    // when choosing the probe-filter join's build side.
    let key_indices = [
        norm_schema.index_of("start")?,
        norm_schema.index_of("allele_string")?,
    ];
    let key_schema = Arc::new(norm_schema.project(&key_indices)?);
    let key_batches = deduped
        .iter()
        .map(|batch| batch.project(&key_indices))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let deduped = partition_batches(deduped, build_config.target_partitions());
    let key_batches = partition_batches(key_batches, build_config.target_partitions());
    let build_ctx = SessionContext::new_with_config_rt(build_config, build_runtime);
    let dedup_view = format!("plugin_{}_dedup", src.plugin_name);
    let key_view = format!("plugin_{}_variation_keys", src.plugin_name);
    let mem = MemTable::try_new(norm_schema.clone(), deduped)
        .map_err(|e| DataFusionError::Execution(format!("dedup memtable: {e}")))?;
    build_ctx.register_table(&dedup_view, Arc::new(mem))?;
    let key_mem = MemTable::try_new(key_schema, key_batches)
        .map_err(|e| DataFusionError::Execution(format!("variation-key memtable: {e}")))?;
    build_ctx.register_table(&key_view, Arc::new(key_mem))?;

    let out_schema = plugin_output_schema(&src.match_columns, &src.value_columns);
    let plugin_dir = output_cache_root.join("plugin").join(&src.plugin_name);
    std::fs::create_dir_all(&plugin_dir)
        .map_err(|e| DataFusionError::Execution(format!("mkdir {}: {e}", plugin_dir.display())))?;
    let file_name = format!("{}.parquet", canonical_chrom_label(chrom));
    let shard_path = plugin_dir.join(&file_name);
    // Write to a sibling temporary path and atomically rename only after the
    // Parquet footer is flushed. A failed/retried hash plan can safely
    // truncate this same path without exposing a partial runtime shard.
    let build_tmp = plugin_dir.join(format!("{file_name}.build.tmp"));
    let scratch = ScratchGuard::new([build_tmp.clone()]);

    // Execute the explicitly sorted parallel plan and stream it directly into
    // one final-order temp shard. If a HashJoin estimate was optimistic, or it
    // fit but starved the downstream external sorter, retry only when the root
    // error and pool trace attribute exhaustion to that hash plan. Other
    // failures propagate as-is.
    let adaptive = match tiered_stream_sorted_adaptive(
        &build_ctx,
        &dedup_view,
        &key_view,
        variation_shard,
        &match_cols,
        pool.as_ref(),
        tracer.as_ref(),
    )
    .await
    {
        Ok(adaptive) => adaptive,
        Err(error) => {
            tracer.report();
            return Err(error);
        }
    };
    let algorithm = adaptive.algorithm;
    let written = write_tiered_shard(
        adaptive.stream,
        &out_schema,
        &build_tmp,
        chrom,
        &src.plugin_name,
    )
    .await;
    let written = match written {
        Err(error) if should_retry_final_hash_plan(&error, algorithm, tracer.as_ref()) => {
            info!(
                "plugin_cache[{}/{chrom}]: final HashJoin plan exhausted the pool \
                 ({read_rows} plugin rows this chrom) -- retrying with DataFusion SortMergeJoin",
                src.plugin_name
            );
            let merge_stream = tiered_stream_sorted_sort_merge(
                &build_ctx,
                &dedup_view,
                &match_cols,
                tracer.as_ref(),
            )
            .await?;
            write_tiered_shard(
                merge_stream,
                &out_schema,
                &build_tmp,
                chrom,
                &src.plugin_name,
            )
            .await
        }
        result => result,
    };
    // Report before propagating: the peak table is most useful on failure.
    tracer.report();
    let (rows, warm, cold) = written?;
    info!(
        "plugin_cache[{}/{chrom}]: tier-join+write done, warm={warm} cold={cold}, {:.1}s elapsed",
        src.plugin_name,
        t_start.elapsed().as_secs_f64()
    );

    // Empty chrom → no shard (matches variation builder cleanup). Its commit
    // removes any stale shard from a previous build so the manifest (rows: 0)
    // matches disk and the runtime never opens a leftover file for an empty
    // chrom.
    if warm + cold == 0 {
        let _ = std::fs::remove_file(&build_tmp);
        return Ok(StagedShard {
            entry: ChromEntry {
                chrom: canonical_chrom_label(chrom),
                file: file_name,
                rows: 0,
                warm: 0,
                cold: 0,
            },
            staged: None,
            live: shard_path,
        });
    }

    scratch.disarm();
    info!(
        "plugin_cache[{}/{chrom}]: staged, rows={rows}, {:.1}s total",
        src.plugin_name,
        t_start.elapsed().as_secs_f64()
    );
    Ok(StagedShard {
        staged: Some(build_tmp),
        live: shard_path,
        entry: ChromEntry {
            chrom: canonical_chrom_label(chrom),
            file: file_name,
            rows,
            warm,
            cold,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::SourceManifest;
    use datafusion::arrow::array::{Float32Array, Int8Array, Int64Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use std::io::Write;
    use std::sync::Arc;

    fn write_gz(path: &std::path::Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    /// Minimal variation-like shard: (chrom, start, allele_string, tier).
    fn write_synthetic_variation(path: &std::path::Path, rows: &[(&str, u32, &str, i8)]) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.0).collect::<Vec<_>>(),
                )),
                Arc::new(UInt32Array::from(
                    rows.iter().map(|r| r.1).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.2).collect::<Vec<_>>(),
                )),
                Arc::new(Int8Array::from(
                    rows.iter().map(|r| r.3).collect::<Vec<_>>(),
                )),
            ],
        )
        .unwrap();
        let file = std::fs::File::create(path).unwrap();
        let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn builds_tiered_shard_with_counts() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("src.tsv.gz");
        // two rows on chr1: 100 (matches warm variation), 300 (miss)
        write_gz(&tsv, "chr1\t100\tA\tG\t0.9\nchr1\t300\tG\tA\t0.7\n");

        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]); // 100 warm (tier 0)

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .unwrap();
        assert_eq!(entry.rows, 2);
        assert_eq!(entry.warm, 1); // start 100 inherited tier 0
        assert_eq!(entry.cold, 1); // start 300 miss -> cold
        assert!(
            out.join("plugin")
                .join("demo")
                .join("chr1.parquet")
                .exists()
        );
    }

    // A staged shard leaves the live one untouched until it is committed;
    // dropping it removes the staging file.
    #[tokio::test(flavor = "multi_thread")]
    async fn a_staged_shard_leaves_the_live_shard_until_committed() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("src.tsv.gz");
        write_gz(&tsv, "chr1\t100\tA\tG\t0.9\n");
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let plugin_dir = out.join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let shard = plugin_dir.join("chr1.parquet");
        let staging = plugin_dir.join("chr1.parquet.build.tmp");
        std::fs::write(&shard, b"previous shard").unwrap();

        let staged = build_plugin_chrom_staged(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .unwrap();
        assert_eq!(staged.entry.rows, 1);
        assert_eq!(std::fs::read(&shard).unwrap(), b"previous shard");
        assert!(staging.exists());
        drop(staged);
        assert!(!staging.exists());
        assert_eq!(std::fs::read(&shard).unwrap(), b"previous shard");

        let entry = build_plugin_chrom_staged(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .unwrap()
            .commit()
            .unwrap();
        assert_eq!(entry.rows, 1);
        assert!(!staging.exists());
        assert_ne!(std::fs::read(&shard).unwrap(), b"previous shard");
    }

    // Two source rows sharing the runtime probe key
    // (start, allele_string, protein_variant) but different scores must collapse
    // to the FIRST in file order (VEP's first-in-file rule) — PR #190 dedup fix.
    #[tokio::test(flavor = "multi_thread")]
    async fn dedups_duplicate_aa_change_keeping_first_in_file() {
        use crate::plugin_cache::lookup::{PluginBufferSlice, PluginLookup, PluginScalar};

        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("am.tsv.gz");
        // Same variant chr3:101 C>T, same aa-change H101Y, two UniProts, two
        // scores. VEP keeps the first (0.0431). A third distinct aa-change at the
        // same position (K55N) must survive.
        write_gz(
            &tsv,
            "chr3\t101\tC\tT\tH101Y\t0.0431\n\
             chr3\t101\tC\tT\tH101Y\t0.0898\n\
             chr3\t101\tC\tT\tK55N\t0.7000\n",
        );

        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("3", 101, "C/T", 0i8)]);

        let toml = format!(
            r##"
plugin_name = "am"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, pv AS protein_variant,
       CAST(score AS FLOAT) AS am_pathogenicity
FROM plugin_am_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "pv",    type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[match_column]]
column = "protein_variant"
template = "{{ref_aa}}{{Protein_position}}{{alt_aa}}"

[[value_columns]]
column = "am_pathogenicity"
csq_field = "am_pathogenicity"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "am.source.toml", &var, &out, "3")
            .await
            .unwrap();
        // H101Y duplicate collapsed → 2 rows survive (H101Y once + K55N).
        assert_eq!(entry.rows, 2, "duplicate H101Y row must be dropped");

        let shard = out.join("plugin").join("am").join("chr3.parquet");
        let lk = PluginLookup::open(
            &shard,
            vec!["protein_variant".into()],
            vec!["am_pathogenicity".into()],
        )
        .await
        .unwrap();
        let batch = lk.take_buffer(&[101]).await.unwrap();
        let slice = PluginBufferSlice::from_batch(&batch, 1, 1).unwrap();
        // The surviving H101Y row carries the FIRST score (0.0431), not 0.0898.
        match slice.probe(101, "C/T", &[Some("H101Y".into())]).unwrap()[0] {
            PluginScalar::F32(v) => assert!((v - 0.0431).abs() < 1e-6, "kept {v}"),
            ref other => panic!("{other:?}"),
        }
        // The distinct aa-change at the same position is untouched.
        match slice.probe(101, "C/T", &[Some("K55N".into())]).unwrap()[0] {
            PluginScalar::F32(v) => assert!((v - 0.7).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
    }

    fn tier_start_batch(tiers: Vec<i8>, starts: Vec<u32>) -> RecordBatch {
        RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("start", DataType::UInt32, false),
                Field::new("tier", DataType::Int8, false),
            ])),
            vec![
                Arc::new(UInt32Array::from(starts)),
                Arc::new(Int8Array::from(tiers)),
            ],
        )
        .unwrap()
    }

    #[test]
    fn tier_start_order_catches_within_tier_disorder() {
        let batch = tier_start_batch(vec![0, 0, 0], vec![10, 20, 15]);
        let err = inspect_tier_start_order(&batch, 0, 1, "1", "demo", None).unwrap_err();
        assert!(err.to_string().contains("not ordered by (tier, start)"));
    }

    #[test]
    fn tier_start_order_catches_cross_batch_tier_regression() {
        let first = tier_start_batch(vec![1], vec![100]);
        let summary = inspect_tier_start_order(&first, 0, 1, "1", "demo", None).unwrap();
        let last_seen = summary.last_key;
        assert_eq!(last_seen, Some((1, 100)));
        let second = tier_start_batch(vec![0], vec![500]);
        let err = inspect_tier_start_order(&second, 0, 1, "1", "demo", last_seen).unwrap_err();
        assert!(err.to_string().contains("not ordered by (tier, start)"));
    }

    #[test]
    fn tier_start_order_accepts_start_reset_at_warm_to_cold_boundary() {
        let batch = tier_start_batch(vec![0, 0, 1, 1], vec![10, 30, 5, 20]);
        let summary = inspect_tier_start_order(&batch, 0, 1, "1", "demo", None).unwrap();
        assert_eq!(summary.last_key, Some((1, 20)));
        assert_eq!((summary.warm, summary.cold), (2, 2));
    }

    // `write_tiered_shard` must surface a genuine reorder as the same error
    // `inspect_tier_start_order` raises directly. The writer must still reject
    // a broken sorted-output contract rather than silently creating a page
    // directory whose lookups can miss rows.
    #[test]
    fn retry_memory_limit_defaults_and_honours_the_env_override() {
        // Serial within one test: these mutate process-wide env.
        unsafe { std::env::remove_var(RETRY_SORT_MEMORY_ENV) };
        assert_eq!(
            retry_sort_memory_limit_bytes(),
            DEFAULT_RETRY_SORT_MEMORY_MIB * 1024 * 1024
        );

        unsafe { std::env::set_var(RETRY_SORT_MEMORY_ENV, "16384") };
        assert_eq!(retry_sort_memory_limit_bytes(), 16384 * 1024 * 1024);

        // A typo must not fail a multi-hour build, and 0 must not create a
        // pool that rejects every allocation.
        for bad in ["", "0", "lots", "-1", "8GiB"] {
            unsafe { std::env::set_var(RETRY_SORT_MEMORY_ENV, bad) };
            assert_eq!(
                retry_sort_memory_limit_bytes(),
                DEFAULT_RETRY_SORT_MEMORY_MIB * 1024 * 1024,
                "{bad:?} should fall back to the default"
            );
        }
        unsafe { std::env::remove_var(RETRY_SORT_MEMORY_ENV) };
    }

    #[test]
    fn retry_memory_default_clears_the_measured_dbnsfp_join_reservation() {
        // chr21, smallest autosome: the join alone reserved 1557 MB for dbNSFP
        // before the sort reserved anything, which the old 2 GiB default could
        // not cover once the rest of the plan is counted. Guard that
        // regression. CADD is deliberately NOT covered -- its reservation
        // scales with whatever pool it is given, so no constant satisfies it.
        let limit = DEFAULT_RETRY_SORT_MEMORY_MIB * 1024 * 1024;
        assert!(
            limit > 2 * 1557 * 1024 * 1024,
            "must clear the measured dbNSFP join with headroom"
        );
    }

    #[test]
    fn target_partitions_default_and_benchmark_override() {
        assert_eq!(target_partitions_from(None, 1), 2);
        assert_eq!(target_partitions_from(None, 16), 16);
        for partitions in [1usize, 2, 4, 8] {
            let raw = partitions.to_string();
            assert_eq!(target_partitions_from(Some(&raw), 16), partitions);
        }
        assert_eq!(target_partitions_from(Some("0"), 16), 16);
        assert_eq!(target_partitions_from(Some("invalid"), 16), 16);
    }

    #[test]
    fn source_reads_use_requested_parallelism_without_round_robin_reordering() {
        let config = source_read_config(8);
        assert_eq!(config.target_partitions(), 8);
        assert!(config.options().optimizer.repartition_file_scans);
        assert!(!config.options().optimizer.enable_round_robin_repartition);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn ordered_parallel_source_replay_preserves_keep_first_across_partitions() {
        fn normalized_batch(score: f32) -> RecordBatch {
            RecordBatch::try_new(
                Arc::new(Schema::new(vec![
                    Field::new("chrom", DataType::Utf8, false),
                    Field::new("start", DataType::Int64, false),
                    Field::new("end", DataType::Int64, false),
                    Field::new("allele_string", DataType::Utf8, false),
                    Field::new("protein_variant", DataType::Utf8, false),
                    Field::new("am_pathogenicity", DataType::Float32, false),
                ])),
                vec![
                    Arc::new(StringArray::from(vec!["1"])),
                    Arc::new(Int64Array::from(vec![100i64])),
                    Arc::new(Int64Array::from(vec![100i64])),
                    Arc::new(StringArray::from(vec!["A/G"])),
                    Arc::new(StringArray::from(vec!["K1R"])),
                    Arc::new(Float32Array::from(vec![score])),
                ],
            )
            .unwrap()
        }

        let first = normalized_batch(0.1);
        let later = normalized_batch(0.9);
        let schema = first.schema();
        let ctx = SessionContext::new_with_config(source_read_config(2));
        ctx.register_table(
            "normalized",
            Arc::new(
                MemTable::try_new(schema, vec![vec![first], vec![later]])
                    .expect("two source partitions"),
            ),
        )
        .unwrap();

        let df = ctx.table("normalized").await.unwrap();
        let (stream, source_partitions) = ordered_parallel_source_stream(df).await.unwrap();
        assert_eq!(source_partitions, 2);
        let deduped = dedup_keep_first(stream, &["protein_variant".into()])
            .await
            .unwrap();
        assert_eq!(deduped.iter().map(RecordBatch::num_rows).sum::<usize>(), 1);
        let score = deduped[0]
            .column(deduped[0].schema().index_of("am_pathogenicity").unwrap())
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap()
            .value(0);
        assert_eq!(score, 0.1, "partition zero's source row must win");
    }

    #[test]
    fn scratch_guard_removes_watched_files_unless_disarmed() {
        let dir = tempfile::tempdir().unwrap();
        let (dropped, kept) = (dir.path().join("a.tmp"), dir.path().join("b.tmp"));
        std::fs::write(&dropped, b"x").unwrap();
        std::fs::write(&kept, b"x").unwrap();

        drop(ScratchGuard::new([dropped.clone()]));
        assert!(
            !dropped.exists(),
            "a failed build must not leave scratch behind"
        );

        ScratchGuard::new([kept.clone()]).disarm();
        assert!(
            kept.exists(),
            "disarm must not delete a consumed scratch file"
        );
    }

    #[test]
    fn scratch_guard_tolerates_files_already_consumed() {
        let dir = tempfile::tempdir().unwrap();
        let gone = dir.path().join("already-removed.tmp");
        drop(ScratchGuard::new([gone.clone()]));
        assert!(!gone.exists());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn write_tiered_shard_flags_disordered_input_as_order_violation() {
        use datafusion::physical_plan::stream::RecordBatchStreamAdapter;

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let make_batch = |starts: &[u32]| {
            let n = starts.len();
            RecordBatch::try_new(
                schema.clone(),
                vec![
                    Arc::new(StringArray::from(vec!["1"; n])),
                    Arc::new(UInt32Array::from(starts.to_vec())),
                    Arc::new(UInt32Array::from(starts.to_vec())),
                    Arc::new(StringArray::from(vec!["A/G"; n])),
                    Arc::new(Int8Array::from(vec![0i8; n])),
                ],
            )
            .unwrap()
        };
        // Two batches, both internally sorted, but the second regresses
        // relative to the first -- exactly what a hash join's build side can
        // produce (each batch may look locally fine, the cross-batch order is
        // what breaks).
        let batches: Vec<Result<RecordBatch>> =
            vec![Ok(make_batch(&[100, 200])), Ok(make_batch(&[50, 300]))];
        let stream: SendableRecordBatchStream = Box::pin(RecordBatchStreamAdapter::new(
            schema.clone(),
            futures::stream::iter(batches),
        ));
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("build.tmp");
        let err = write_tiered_shard(stream, &schema, &path, "1", "demo")
            .await
            .unwrap_err();
        assert!(
            err.to_string().contains("not ordered by (tier, start)"),
            "unexpected error: {err}"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn write_tiered_shard_writes_final_order_directly() {
        use datafusion::physical_plan::stream::RecordBatchStreamAdapter;

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let make_batch = |tiers: Vec<i8>, starts: Vec<u32>| {
            let n = starts.len();
            RecordBatch::try_new(
                schema.clone(),
                vec![
                    Arc::new(StringArray::from(vec!["1"; n])),
                    Arc::new(UInt32Array::from(starts.clone())),
                    Arc::new(UInt32Array::from(starts)),
                    Arc::new(StringArray::from(vec!["A/G"; n])),
                    Arc::new(Int8Array::from(tiers)),
                ],
            )
            .unwrap()
        };
        let batches: Vec<Result<RecordBatch>> = vec![
            Ok(make_batch(vec![0, 0], vec![100, 300])),
            Ok(make_batch(vec![1, 1], vec![50, 200])),
        ];
        let stream: SendableRecordBatchStream = Box::pin(RecordBatchStreamAdapter::new(
            schema.clone(),
            futures::stream::iter(batches),
        ));
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("build.tmp");
        assert_eq!(
            write_tiered_shard(stream, &schema, &path, "1", "demo")
                .await
                .unwrap(),
            (4, 2, 2)
        );

        let file = std::fs::File::open(path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap()
            .build()
            .unwrap();
        let mut keys = Vec::new();
        for batch in reader {
            let batch = batch.unwrap();
            let starts = batch
                .column(batch.schema().index_of("start").unwrap())
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();
            let tiers = batch
                .column(batch.schema().index_of("tier").unwrap())
                .as_any()
                .downcast_ref::<Int8Array>()
                .unwrap();
            keys.extend((0..batch.num_rows()).map(|i| (tiers.value(i), starts.value(i))));
        }
        assert_eq!(keys, vec![(0, 100), (0, 300), (1, 50), (1, 200)]);
    }

    // C1 regression, end-to-end: `HashJoinExec`'s `CollectLeft` mode collects
    // the plugin (LEFT of the LEFT JOIN) into a hash table and streams the
    // variation probe; a MATCHED row's emission order follows the probe's
    // scan (fine, since a real variation shard is itself position-sorted),
    // but an UNMATCHED plugin row is appended afterward by iterating the
    // plugin's own original row order. A plugin whose rows are already
    // ascending never regresses under this mechanism regardless of its size
    // relative to the probe -- but a source that isn't globally sorted (two
    // files concatenated, an upstream ordering bug) surfaces here as soon as
    // any of its rows miss the variation probe. Confirmed empirically (not
    // just asserted) via a temporary EXPLAIN probe against this exact
    // mechanism before writing this test.
    //
    // `build_plugin_chrom` must make the parallel result globally sorted.
    #[tokio::test(flavor = "multi_thread")]
    async fn sparse_plugin_with_disordered_source_still_builds_sorted() {
        let dir = tempfile::tempdir().unwrap();

        // Plugin: 20 rows, DESCENDING -- simulates a source that isn't
        // globally position-sorted (e.g. two per-chrom files concatenated
        // without a merge sort, à la the original CADD SNV+indel bug).
        let mut tsv_body = String::new();
        for i in (1..=20u32).rev() {
            tsv_body.push_str(&format!("chr1\t{}\tA\tG\t0.{i}\n", i * 100));
        }
        let tsv = dir.path().join("sparse.tsv.gz");
        write_gz(&tsv, &tsv_body);

        // Variation shard: none of the plugin's positions have an entry, so
        // every plugin row misses (tier=1/cold). This made the old unsorted
        // HashJoin path emit descending cold rows; the adaptive parallel path
        // must satisfy its explicit ORDER BY contract. The other 5000 rows are
        // unrelated positions, matching production's scale (chr22's variation
        // shard vs. ClinVar's row count is a ~157x skew).
        let mut var_rows: Vec<(&str, u32, String, i8)> = Vec::new();
        for i in 0..5000u32 {
            var_rows.push(("1", 10_000_000 + i, "C/T".to_string(), 1i8));
        }
        let var_path = dir.path().join("var.parquet");
        let var_rows_ref: Vec<(&str, u32, &str, i8)> = var_rows
            .iter()
            .map(|(c, s, a, t)| (*c, *s, a.as_str(), *t))
            .collect();
        write_synthetic_variation(&var_path, &var_rows_ref);

        let toml = format!(
            r##"
plugin_name = "sparse"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_sparse_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "sparse.source.toml", &var_path, &out, "1")
            .await
            .unwrap();
        assert_eq!(entry.rows, 20);
        assert_eq!(entry.warm, 0);
        assert_eq!(
            entry.cold, 20,
            "all 20 plugin rows miss the variation shard"
        );

        // The on-disk shard must actually be position-ascending -- this is
        // the real invariant `PageDir::resolve_ranges` depends on, and the
        // one thing a passing `entry.rows` count alone can't prove.
        let shard_path = out.join("plugin").join("sparse").join("chr1.parquet");
        let file = std::fs::File::open(&shard_path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap()
            .build()
            .unwrap();
        let mut starts = Vec::new();
        for batch in reader {
            let batch = batch.unwrap();
            let col = batch
                .column(batch.schema().index_of("start").unwrap())
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .clone();
            starts.extend((0..col.len()).map(|i| col.value(i)));
        }
        assert_eq!(starts.len(), 20);
        let mut sorted = starts.clone();
        sorted.sort_unstable();
        assert_eq!(starts, sorted, "shard on disk is not position-ascending");
    }

    // --chrom M / chrM / MT must all select the MT rows (data is folded to "MT"
    // by canonical_contig), not silently produce a 0-row shard (PR #190 M1).
    #[tokio::test(flavor = "multi_thread")]
    async fn builds_mt_shard_from_any_mt_alias() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("mt.tsv.gz");
        write_gz(&tsv, "chrM\t100\tA\tG\t0.9\n");
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("MT", 100, "A/G", 0i8)]);
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        for (i, alias) in ["M", "chrM", "MT"].iter().enumerate() {
            let out = dir.path().join(format!("out{i}"));
            let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, alias)
                .await
                .unwrap();
            assert_eq!(entry.rows, 1, "alias '{alias}' should build the MT row");
        }
    }
}
