//! Variation-tier inheritance: within one validated chromosome, LEFT-joins the
//! normalized plugin rows against the variation shard on `(start, allele_string)`
//! and inherits the variation record's `tier` (`COALESCE(v.tier, 1)` — no match
//! → cold). Variation columns drop; only the plugin's value columns plus `tier`
//! survive.

use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::utils::memory::get_record_batch_memory_size;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::MemTable;
use datafusion::execution::memory_pool::MemoryPool;
use datafusion::physical_plan::{
    ExecutionPlan, SendableRecordBatchStream, collect, displayable, execute_stream,
};
use datafusion::prelude::SessionContext;
use log::info;

use crate::plugin_cache::join_strategy::{
    JoinAlgorithm, JoinDecision, choose_for_hash_plan, contains_sort_merge_join,
};
use crate::plugin_cache::mem_trace::TracingPool;

/// SQL that LEFT-joins `normalized_view` to a registered `variation_probe`
/// exposing `(start, allele_string, tier)` and inherits `tier`. Both inputs are
/// already restricted to one chromosome, so carrying that constant string in
/// the physical join key only increases memory. The value columns of
/// `normalized_view` pass through (`p.*`); variation columns drop.
pub fn tier_sql(normalized_view: &str, variation_probe: &str) -> String {
    format!(
        "SELECT p.*, CAST(COALESCE(v.tier, 1) AS TINYINT) AS tier \
         FROM {normalized_view} p \
         LEFT JOIN {variation_probe} v \
         ON p.start = v.start AND p.allele_string = v.allele_string"
    )
}

/// Same join, with the exact physical shard order: warm rows first, then cold,
/// and position-ascending inside each tier. This lets the consumer write one
/// Parquet file directly instead of encoding two tier files and decoding them
/// again for a final merge.
pub fn tier_sql_sorted(normalized_view: &str, variation_probe: &str) -> String {
    format!(
        "SELECT * FROM ({}) tiered ORDER BY tier, start",
        tier_sql(normalized_view, variation_probe)
    )
}

/// Register the variation shard as a `tier`-bearing probe view, then stream the
/// tiered rows (tier inherited from the variation record) in whatever order
/// the join happens to produce.
pub async fn tiered_stream(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
) -> Result<SendableRecordBatchStream> {
    tiered_stream_impl(ctx, normalized_view, variation_shard, false).await
}

/// Same as `tiered_stream`, but output is guaranteed lexicographically sorted
/// by `(tier, start)` (see `tier_sql_sorted`).
pub async fn tiered_stream_sorted(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
) -> Result<SendableRecordBatchStream> {
    tiered_stream_impl(ctx, normalized_view, variation_shard, true).await
}

pub(crate) struct AdaptiveTieredStream {
    pub stream: SendableRecordBatchStream,
    pub algorithm: JoinAlgorithm,
}

fn is_consumer_resources_exhausted(error: &DataFusionError, consumer: &str) -> bool {
    matches!(
        error.find_root(),
        DataFusionError::ResourcesExhausted(message) if message.contains(consumer)
    )
}

pub(crate) fn should_retry_hash_join(
    error: &DataFusionError,
    algorithm: JoinAlgorithm,
    pool: &TracingPool,
) -> bool {
    algorithm == JoinAlgorithm::Hash
        && is_consumer_resources_exhausted(error, "HashJoinInput")
        && pool.hash_join_failed()
}

/// Retry the final sorted query when either its unspillable hash build or the
/// downstream external sorter exhausts the pool. The sorter case matters when
/// the hash estimate technically fits but leaves too little working memory for
/// `ORDER BY tier, start`; dropping the hash plan before replanning as SMJ frees
/// that reservation without relying on a fixed headroom guess.
pub(crate) fn should_retry_final_hash_plan(
    error: &DataFusionError,
    algorithm: JoinAlgorithm,
    pool: &TracingPool,
) -> bool {
    should_retry_hash_join(error, algorithm, pool)
        || (algorithm == JoinAlgorithm::Hash
            && is_consumer_resources_exhausted(error, "ExternalSorter")
            && pool.external_sorter_failed())
}

async fn set_prefer_hash_join(ctx: &SessionContext, prefer: bool) -> Result<()> {
    ctx.sql(&format!(
        "SET datafusion.optimizer.prefer_hash_join = {prefer}"
    ))
    .await?;
    Ok(())
}

fn log_decision(stage: &str, target_partitions: usize, decision: &JoinDecision) {
    info!(
        "plugin_cache: adaptive join stage={stage} target_partitions={target_partitions} algorithm={:?} reason={:?} \
         build_rows={:?} build_data_bytes={:?} hash_build_bytes={:?} available_bytes={:?}",
        decision.algorithm,
        decision.reason,
        decision.build_rows,
        decision.build_data_bytes,
        decision.hash_build_bytes,
        decision.available_bytes,
    );
}

fn log_plan(stage: &str, plan: &dyn ExecutionPlan) {
    if std::env::var("VEPYR_EXPLAIN_TIER_JOIN").is_ok() {
        info!(
            "plugin_cache: adaptive join stage={stage} physical plan:\n{}",
            displayable(plan).indent(true)
        );
    }
}

async fn sort_merge_plan(
    ctx: &SessionContext,
    sql: &str,
    stage: &str,
) -> Result<Arc<dyn ExecutionPlan>> {
    set_prefer_hash_join(ctx, false).await?;
    let plan = ctx.sql(sql).await?.create_physical_plan().await?;
    if !contains_sort_merge_join(plan.as_ref()) {
        return datafusion::common::exec_err!(
            "adaptive join stage {stage} requested SortMergeJoin but DataFusion emitted:\n{}",
            displayable(plan.as_ref()).indent(true)
        );
    }
    Ok(plan)
}

async fn adaptive_plan(
    ctx: &SessionContext,
    sql: &str,
    stage: &str,
    pool: &dyn MemoryPool,
) -> Result<(Arc<dyn ExecutionPlan>, JoinAlgorithm)> {
    set_prefer_hash_join(ctx, true).await?;
    let hash_plan = ctx.sql(sql).await?.create_physical_plan().await?;
    let decision = choose_for_hash_plan(&hash_plan, pool)?;
    log_decision(stage, ctx.copied_config().target_partitions(), &decision);
    if decision.algorithm == JoinAlgorithm::Hash {
        return Ok((hash_plan, JoinAlgorithm::Hash));
    }

    Ok((
        sort_merge_plan(ctx, sql, stage).await?,
        JoinAlgorithm::SortMerge,
    ))
}

/// Concatenate `parts` into a batch that owns its buffers even when there is
/// only one part. Arrow's one-input `concat_batches` fast path is zero-copy, so
/// split that input into two logical parts to force the normal copying path.
fn concat_probe_parts_owned(schema: &SchemaRef, parts: &[RecordBatch]) -> Result<RecordBatch> {
    debug_assert!(!parts.is_empty());
    if parts.len() == 1 {
        let batch = &parts[0];
        let midpoint = batch.num_rows() / 2;
        let first = batch.slice(0, midpoint);
        let second = batch.slice(midpoint, batch.num_rows() - midpoint);
        Ok(concat_batches(schema, [&first, &second])?)
    } else {
        Ok(concat_batches(schema, parts)?)
    }
}

/// Copy input into independently allocated batches no larger than the
/// execution batch size. DataFusion 53 wraps MemTables in `BatchSplitStream`;
/// oversized parents would otherwise be split into zero-copy children and the
/// hash join would reserve the full parent once per child.
fn rechunk_probe_owned(
    schema: &SchemaRef,
    input: Vec<RecordBatch>,
    target_rows: usize,
) -> Result<Vec<RecordBatch>> {
    if target_rows == 0 {
        return datafusion::common::exec_err!("probe materialization batch size must be non-zero");
    }

    let total_rows: usize = input.iter().map(RecordBatch::num_rows).sum();
    let mut output = Vec::with_capacity(total_rows.div_ceil(target_rows));
    let mut parts = Vec::new();
    let mut part_rows = 0usize;

    for batch in input {
        let mut offset = 0usize;
        while offset < batch.num_rows() {
            let take = (target_rows - part_rows).min(batch.num_rows() - offset);
            parts.push(batch.slice(offset, take));
            part_rows += take;
            offset += take;

            if part_rows == target_rows {
                output.push(concat_probe_parts_owned(schema, &parts)?);
                parts.clear();
                part_rows = 0;
            }
        }
    }

    if part_rows != 0 {
        output.push(concat_probe_parts_owned(schema, &parts)?);
    }
    debug_assert_eq!(
        output.iter().map(RecordBatch::num_rows).sum::<usize>(),
        total_rows
    );
    debug_assert!(output.iter().all(|batch| batch.num_rows() <= target_rows));
    Ok(output)
}

/// Distribute already-owned batches across source partitions without copying
/// their Arrow buffers. Greedy byte balancing avoids assigning all large
/// batches to one worker while adapting naturally when a chromosome has only
/// enough batches to occupy a subset of the configured workers.
pub(crate) fn partition_batches(
    batches: Vec<RecordBatch>,
    target_partitions: usize,
) -> Vec<Vec<RecordBatch>> {
    let partition_count = target_partitions.max(1).min(batches.len().max(1));
    let mut partitions = (0..partition_count).map(|_| Vec::new()).collect::<Vec<_>>();
    let mut partition_bytes = vec![0usize; partition_count];

    for batch in batches {
        let partition = partition_bytes
            .iter()
            .enumerate()
            .min_by_key(|(_, bytes)| **bytes)
            .map(|(index, _)| index)
            .expect("at least one batch partition");
        partition_bytes[partition] =
            partition_bytes[partition].saturating_add(get_record_batch_memory_size(&batch));
        partitions[partition].push(batch);
    }
    partitions
}

/// Replace the `plugin_variation_probe` view with an equivalent table whose
/// batches own their buffers.
///
/// The view is a `GROUP BY`, and `GroupedHashAggregateStream` emits its output
/// as batches that SHARE underlying buffers. `HashJoinExec` reserves memory per
/// build-side batch using `get_record_batch_memory_size`, which de-duplicates
/// buffers only WITHIN one batch -- so a shared buffer is charged once per
/// batch that references it. Measured on chr21's variation shard: the grouped
/// output is 181.4 MB of distinct buffers but accounts as 5756.8 MB across its
/// 1791 batches, a 31.7x over-count. The same scan WITHOUT the `GROUP BY`
/// accounts at exactly 1.0x, so this is the aggregate's doing, not the scan's.
///
/// CADD exposes this because the tier join's build side is this probe (its plan
/// swaps to `join_type=Right` because 120M plugin rows dwarf the 10.8M-row
/// probe). An earlier 1,048,576-row materialization removed the aggregate's
/// sharing, but DataFusion's MemTable scan split each parent into 128 zero-copy
/// 8,192-row slices and charged the full parent for every slice. The resulting
/// reservation exceeded 20 GB for ~180 MB of distinct probe buffers.
///
/// Materializing into owned batches no larger than the execution batch size
/// gives each batch fresh, unshared buffers and prevents DataFusion 53 from
/// re-slicing large MemTable parents into zero-copy children. The reservation
/// therefore matches reality. It also lifts the aggregate out of the join plan,
/// so it runs once and the optimizer sees exact statistics for side selection.
fn replace_probe_table(
    ctx: &SessionContext,
    schema: SchemaRef,
    collected: Vec<RecordBatch>,
) -> Result<()> {
    let batches = rechunk_probe_owned(&schema, collected, ctx.copied_config().batch_size())?;
    let partitions = partition_batches(batches, ctx.copied_config().target_partitions());
    ctx.deregister_table("plugin_variation_probe")?;
    ctx.register_table(
        "plugin_variation_probe",
        Arc::new(MemTable::try_new(schema, partitions)?),
    )?;
    Ok(())
}

async fn materialize_probe(ctx: &SessionContext) -> Result<()> {
    let df = ctx.sql("SELECT * FROM plugin_variation_probe").await?;
    let schema: SchemaRef = df.schema().inner().clone();
    let collected = df.collect().await?;
    replace_probe_table(ctx, schema, collected)
}

async fn materialize_probe_adaptive(
    ctx: &SessionContext,
    pool: &dyn MemoryPool,
    tracer: &TracingPool,
) -> Result<()> {
    const SQL: &str = "SELECT * FROM plugin_variation_probe";
    let (hash_or_merge_plan, algorithm) = adaptive_plan(ctx, SQL, "probe", pool).await?;
    let schema = hash_or_merge_plan.schema();
    log_plan("probe", hash_or_merge_plan.as_ref());
    tracer.clear_failures();

    let collected = match collect(hash_or_merge_plan, ctx.task_ctx()).await {
        Ok(batches) => batches,
        Err(error) if should_retry_hash_join(&error, algorithm, tracer) => {
            info!(
                "plugin_cache: probe HashJoin exhausted its reservation; retrying with DataFusion SortMergeJoin"
            );
            let merge_plan = sort_merge_plan(ctx, SQL, "probe_runtime_fallback").await?;
            log_plan("probe_runtime_fallback", merge_plan.as_ref());
            tracer.clear_failures();
            collect(merge_plan, ctx.task_ctx()).await?
        }
        Err(error) => return Err(error),
    };
    replace_probe_table(ctx, schema, collected)
}

async fn register_variation_probe_view(
    ctx: &SessionContext,
    variation_key_view: &str,
    variation_shard: &Path,
) -> Result<()> {
    let shard = variation_shard.to_string_lossy();
    ctx.register_parquet(
        "plugin_variation_raw",
        shard.as_ref(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await?;
    ctx.sql(&format!(
        // GROUP BY (not DISTINCT): the variation shard can carry multiple
        // source rows for the same (start, allele_string) within this
        // already-selected chromosome (e.g.
        // distinct dbSNP/COSMIC-origin entries for one variant). If those
        // duplicates ever disagree on `tier` (one warm, one cold), a plain
        // `SELECT DISTINCT start, allele_string, tier` keeps BOTH rows
        // — the join key is then non-unique on the build side, so the tier
        // LEFT JOIN below fans out the plugin row once per surviving variant,
        // silently inflating the cache with duplicate rows for that key
        // (worse: two rows with the same (start, allele_string) but different
        // tier, which the point-lookup path never expects). `MIN(tier)`
        // collapses to exactly one row per key regardless of how many
        // conflicting-tier duplicates the raw shard has — 0 (warm) wins over
        // 1 (cold), i.e. any evidence of a known/frequent variant marks it
        // warm. This also guarantees `HashJoinExec`'s build side can never
        // fan out, which is a (data-dependent) join-skew risk this GROUP BY
        // removes entirely rather than merely reduces.
        //
        // The LEFT SEMI JOIN against the plugin's own per-chromosome keys
        // restricts the aggregation to keys the plugin table actually has,
        // computed BEFORE
        // the `GROUP BY` rather than after: for a sparse plugin (far fewer
        // rows than the chromosome's full variation shard), grouping the
        // *entire* shard
        // just to inherit tier for a handful of keys is wasted work the LEFT
        // JOIN below never uses (variation rows with no matching plugin key
        // can never survive a `p.*`-projected LEFT JOIN regardless).
        // Semantically a no-op — this narrows which rows get grouped, not
        // which key wins ties or what a hit resolves to. A true semi join
        // cannot multiply variation rows when a plugin key is duplicated, so
        // it removes the DISTINCT aggregate and its shared-buffer output from
        // directly below the join while preserving the previous semantics.
        "CREATE OR REPLACE VIEW plugin_variation_probe AS \
         SELECT v.start, v.allele_string, MIN(v.tier) AS tier \
         FROM plugin_variation_raw v \
         LEFT SEMI JOIN {variation_key_view} p \
         ON v.start = p.start AND v.allele_string = p.allele_string \
         GROUP BY v.start, v.allele_string"
    ))
    .await?;
    Ok(())
}

async fn tiered_stream_impl(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
    sorted: bool,
) -> Result<SendableRecordBatchStream> {
    register_variation_probe_view(ctx, normalized_view, variation_shard).await?;
    materialize_probe(ctx).await?;
    let sql = if sorted {
        tier_sql_sorted(normalized_view, "plugin_variation_probe")
    } else {
        tier_sql(normalized_view, "plugin_variation_probe")
    };
    let df = ctx.sql(&sql).await?;
    // Opt-in plan dump. The tier join's memory behaviour depends on which side
    // of each HashJoinExec the optimizer collects, and that is decided from
    // statistics -- so it cannot be reasoned out from the SQL, it has to be
    // read off the physical plan for the actual data. Costs nothing unset.
    if std::env::var("VEPYR_EXPLAIN_TIER_JOIN").is_ok() {
        let plan = df.clone().create_physical_plan().await?;
        info!(
            "plugin_cache: tier join physical plan (sorted={sorted}):\n{}",
            datafusion::physical_plan::displayable(plan.as_ref()).indent(true)
        );
    }
    df.execute_stream().await
}

/// Build both joins with DataFusion, selecting HashJoin only when its actual
/// optimized build-side statistics fit the active pool. The final SQL always
/// carries an explicit ORDER BY because parallel joins do not preserve a
/// globally useful source order.
pub(crate) async fn tiered_stream_sorted_adaptive(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_key_view: &str,
    variation_shard: &Path,
    pool: &dyn MemoryPool,
    tracer: &TracingPool,
) -> Result<AdaptiveTieredStream> {
    register_variation_probe_view(ctx, variation_key_view, variation_shard).await?;
    materialize_probe_adaptive(ctx, pool, tracer).await?;

    let sql = tier_sql_sorted(normalized_view, "plugin_variation_probe");
    let (plan, algorithm) = adaptive_plan(ctx, &sql, "tier", pool).await?;
    log_plan("tier", plan.as_ref());
    tracer.clear_failures();
    match execute_stream(plan, ctx.task_ctx()) {
        Ok(stream) => Ok(AdaptiveTieredStream { stream, algorithm }),
        Err(error) if should_retry_final_hash_plan(&error, algorithm, tracer) => {
            info!(
                "plugin_cache: final HashJoin plan exhausted the pool before streaming; retrying with DataFusion SortMergeJoin"
            );
            Ok(AdaptiveTieredStream {
                stream: tiered_stream_sorted_sort_merge(ctx, normalized_view, tracer).await?,
                algorithm: JoinAlgorithm::SortMerge,
            })
        }
        Err(error) => Err(error),
    }
}

/// Re-plan only the already-materialized final join as SortMergeJoin. Used
/// when a hash build passed the estimate but DataFusion rejected its runtime
/// reservation while the output stream was being consumed.
pub(crate) async fn tiered_stream_sorted_sort_merge(
    ctx: &SessionContext,
    normalized_view: &str,
    tracer: &TracingPool,
) -> Result<SendableRecordBatchStream> {
    let sql = tier_sql_sorted(normalized_view, "plugin_variation_probe");
    let plan = sort_merge_plan(ctx, &sql, "tier_runtime_fallback").await?;
    log_plan("tier_runtime_fallback", plan.as_ref());
    tracer.clear_failures();
    execute_stream(plan, ctx.task_ctx())
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use datafusion::execution::memory_pool::{
        FairSpillPool, MemoryConsumer, MemoryLimit, MemoryReservation, UnboundedMemoryPool,
    };
    use datafusion::execution::runtime_env::RuntimeEnvBuilder;
    use datafusion::prelude::{SessionConfig, col};
    use futures::{StreamExt, TryStreamExt};
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// Reports a roomy limit to the preflight selector but rejects every hash
    /// build reservation. Other consumers are unbounded so the forced SMJ can
    /// complete; this models optimistic statistics deterministically.
    #[derive(Debug, Default)]
    struct RejectHashPool {
        inner: UnboundedMemoryPool,
        rejected_hash_grows: AtomicUsize,
    }

    impl MemoryPool for RejectHashPool {
        fn grow(&self, reservation: &MemoryReservation, additional: usize) {
            self.inner.grow(reservation, additional);
        }

        fn shrink(&self, reservation: &MemoryReservation, shrink: usize) {
            self.inner.shrink(reservation, shrink);
        }

        fn try_grow(&self, reservation: &MemoryReservation, additional: usize) -> Result<()> {
            if reservation.consumer().name().starts_with("HashJoinInput") {
                self.rejected_hash_grows.fetch_add(1, Ordering::Relaxed);
                return Err(DataFusionError::ResourcesExhausted(format!(
                    "Additional allocation failed for {}",
                    reservation.consumer().name()
                )));
            }
            self.inner.try_grow(reservation, additional)
        }

        fn reserved(&self) -> usize {
            self.inner.reserved()
        }

        fn memory_limit(&self) -> MemoryLimit {
            MemoryLimit::Finite(1024 * 1024)
        }
    }

    fn rejecting_context() -> (SessionContext, Arc<RejectHashPool>, Arc<TracingPool>) {
        let rejecting_pool = Arc::new(RejectHashPool::default());
        let base_pool: Arc<dyn MemoryPool> = Arc::clone(&rejecting_pool) as _;
        let tracer = Arc::new(TracingPool::new(base_pool));
        let runtime_pool: Arc<dyn MemoryPool> = Arc::clone(&tracer) as _;
        let runtime = RuntimeEnvBuilder::new()
            .with_memory_pool(runtime_pool)
            .build_arc()
            .unwrap();
        let ctx = SessionContext::new_with_config_rt(
            SessionConfig::new()
                .with_target_partitions(2)
                .set_bool("datafusion.optimizer.prefer_hash_join", true)
                .set_usize(
                    "datafusion.optimizer.hash_join_single_partition_threshold",
                    usize::MAX,
                )
                .set_usize(
                    "datafusion.optimizer.hash_join_single_partition_threshold_rows",
                    usize::MAX,
                ),
            runtime,
        );
        (ctx, rejecting_pool, tracer)
    }

    #[test]
    fn owned_batches_are_distributed_across_available_workers() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "value",
            DataType::UInt32,
            false,
        )]));
        let batches = (0..6u32)
            .map(|value| {
                RecordBatch::try_new(
                    Arc::clone(&schema),
                    vec![Arc::new(UInt32Array::from(vec![value; 16]))],
                )
                .unwrap()
            })
            .collect::<Vec<_>>();

        let partitions = partition_batches(batches, 3);
        assert_eq!(partitions.len(), 3);
        assert!(partitions.iter().all(|partition| partition.len() == 2));
        assert_eq!(
            partitions
                .iter()
                .flatten()
                .map(RecordBatch::num_rows)
                .sum::<usize>(),
            96
        );
    }

    #[test]
    fn hash_retry_requires_matching_error_and_pool_attribution() {
        let tracer = Arc::new(TracingPool::new(Arc::new(FairSpillPool::new(1))));
        let pool: Arc<dyn MemoryPool> = Arc::clone(&tracer) as _;
        assert!(matches!(pool.memory_limit(), MemoryLimit::Finite(1)));
        let hash_reservation = MemoryConsumer::new("HashJoinInput").register(&pool);
        let hash_error = hash_reservation.try_grow(2).unwrap_err();

        assert!(should_retry_hash_join(
            &hash_error,
            JoinAlgorithm::Hash,
            tracer.as_ref()
        ));
        assert!(!should_retry_hash_join(
            &hash_error,
            JoinAlgorithm::SortMerge,
            tracer.as_ref()
        ));

        tracer.clear_failures();
        let sort_reservation = MemoryConsumer::new("ExternalSorter").register(&pool);
        let sort_error = sort_reservation.try_grow(2).unwrap_err();
        assert!(!should_retry_hash_join(
            &sort_error,
            JoinAlgorithm::Hash,
            tracer.as_ref()
        ));
        assert!(should_retry_final_hash_plan(
            &sort_error,
            JoinAlgorithm::Hash,
            tracer.as_ref()
        ));
        assert!(!should_retry_final_hash_plan(
            &sort_error,
            JoinAlgorithm::SortMerge,
            tracer.as_ref()
        ));
    }

    #[test]
    fn final_hash_plan_retries_when_live_hash_starves_external_sorter() {
        let tracer = Arc::new(TracingPool::new(Arc::new(FairSpillPool::new(100))));
        let pool: Arc<dyn MemoryPool> = Arc::clone(&tracer) as _;
        let hash_reservation = MemoryConsumer::new("HashJoinInput").register(&pool);
        hash_reservation.try_grow(90).unwrap();

        let sort_reservation = MemoryConsumer::new("ExternalSorter").register(&pool);
        let error = sort_reservation.try_grow(20).unwrap_err();
        assert!(should_retry_final_hash_plan(
            &error,
            JoinAlgorithm::Hash,
            tracer.as_ref()
        ));

        // Replanning drops the unspillable build, leaving the same bounded
        // pool able to grant the sorter's working reservation.
        drop(hash_reservation);
        sort_reservation.try_grow(20).unwrap();
    }

    #[tokio::test]
    async fn adaptive_plan_replans_with_datafusion_sort_merge_join() {
        let ctx = SessionContext::new_with_config(
            SessionConfig::new()
                .with_target_partitions(2)
                .set_usize(
                    "datafusion.optimizer.hash_join_single_partition_threshold",
                    usize::MAX,
                )
                .set_usize(
                    "datafusion.optimizer.hash_join_single_partition_threshold_rows",
                    usize::MAX,
                ),
        );
        let schema = Arc::new(Schema::new(vec![Field::new(
            "value",
            DataType::UInt32,
            false,
        )]));
        for (name, rows) in [("build", 128usize), ("probe", 4096usize)] {
            ctx.register_batch(
                name,
                RecordBatch::try_new(
                    Arc::clone(&schema),
                    vec![Arc::new(UInt32Array::from_iter_values(
                        (0..rows).map(|value| value as u32),
                    ))],
                )
                .unwrap(),
            )
            .unwrap();
        }

        let (plan, algorithm) = adaptive_plan(
            &ctx,
            "SELECT p.value FROM probe p INNER JOIN build b ON p.value = b.value",
            "test",
            &FairSpillPool::new(64),
        )
        .await
        .unwrap();
        assert_eq!(algorithm, JoinAlgorithm::SortMerge);
        assert!(contains_sort_merge_join(plan.as_ref()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn probe_runtime_hash_rejection_retries_with_sort_merge() {
        let (ctx, rejecting_pool, tracer) = rejecting_context();
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        for name in ["variation_side", "plugin_side"] {
            ctx.register_batch(
                name,
                RecordBatch::try_new(
                    Arc::clone(&schema),
                    vec![
                        Arc::new(StringArray::from(vec!["1", "1"])),
                        Arc::new(UInt32Array::from(vec![100u32, 200])),
                        Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                        Arc::new(Int8Array::from(vec![0i8, 1])),
                    ],
                )
                .unwrap(),
            )
            .unwrap();
        }
        ctx.sql(
            "CREATE VIEW plugin_variation_probe AS \
             SELECT v.start, v.allele_string, v.tier \
             FROM variation_side v INNER JOIN plugin_side p \
             ON v.start = p.start AND v.allele_string = p.allele_string",
        )
        .await
        .unwrap();

        materialize_probe_adaptive(&ctx, tracer.as_ref(), tracer.as_ref())
            .await
            .unwrap();
        assert_eq!(
            rejecting_pool.rejected_hash_grows.load(Ordering::Relaxed),
            1
        );
        let probe = ctx
            .sql("SELECT * FROM plugin_variation_probe")
            .await
            .unwrap();
        assert_eq!(
            probe
                .schema()
                .fields()
                .iter()
                .map(|field| field.name().as_str())
                .collect::<Vec<_>>(),
            vec!["start", "allele_string", "tier"]
        );
        let rows = probe
            .collect()
            .await
            .unwrap()
            .iter()
            .map(RecordBatch::num_rows)
            .sum::<usize>();
        assert_eq!(rows, 2);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn final_runtime_hash_rejection_retries_with_sort_merge() {
        let (ctx, rejecting_pool, tracer) = rejecting_context();

        let dir = tempfile::tempdir().unwrap();
        let variation_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let variation = RecordBatch::try_new(
            Arc::clone(&variation_schema),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Int8Array::from(vec![0i8, 1])),
            ],
        )
        .unwrap();
        let variation_path = dir.path().join("variation.parquet");
        let file = std::fs::File::create(&variation_path).unwrap();
        let mut writer = ArrowWriter::try_new(file, variation_schema, None).unwrap();
        writer.write(&variation).unwrap();
        writer.close().unwrap();

        let plugin = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("end", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("demo_score", DataType::Float32, true),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 300])),
                Arc::new(UInt32Array::from(vec![100u32, 300])),
                Arc::new(StringArray::from(vec!["A/G", "G/A"])),
                Arc::new(Float32Array::from(vec![Some(0.9f32), Some(0.7)])),
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_demo_norm", plugin).unwrap();

        let adaptive = tiered_stream_sorted_adaptive(
            &ctx,
            "plugin_demo_norm",
            "plugin_demo_norm",
            &variation_path,
            tracer.as_ref(),
            tracer.as_ref(),
        )
        .await
        .unwrap();
        assert_eq!(adaptive.algorithm, JoinAlgorithm::Hash);

        let error = adaptive
            .stream
            .try_collect::<Vec<RecordBatch>>()
            .await
            .unwrap_err();
        assert!(should_retry_hash_join(
            &error,
            adaptive.algorithm,
            tracer.as_ref()
        ));
        assert!(
            rejecting_pool.rejected_hash_grows.load(Ordering::Relaxed) >= 1,
            "the final tier HashJoin must reach the injected rejecting pool"
        );

        let batches = tiered_stream_sorted_sort_merge(&ctx, "plugin_demo_norm", tracer.as_ref())
            .await
            .unwrap()
            .try_collect::<Vec<RecordBatch>>()
            .await
            .unwrap();
        assert_eq!(batches.iter().map(RecordBatch::num_rows).sum::<usize>(), 2);
    }

    // Codex review, PR #196: for a sparse plugin the pre-fix
    // `plugin_variation_probe` unconditionally grouped the ENTIRE variation
    // shard, including keys no plugin row will ever reference. The semi-join
    // narrowing must be a pure performance change, not a semantic one --
    // verify a large-relative-to-plugin variation shard (rows the plugin
    // never touches, PLUS a conflicting-tier duplicate at a key the plugin
    // DOES touch) still tiers exactly as before.
    #[tokio::test(flavor = "multi_thread")]
    async fn tiered_stream_ignores_variation_rows_no_plugin_key_references() {
        let dir = tempfile::tempdir().unwrap();
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        // 100: two conflicting-tier duplicates (must still collapse to warm).
        // 999/1000/1001: irrelevant rows no plugin key ever probes.
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 100, 999, 1000, 1001])),
                Arc::new(StringArray::from(vec!["A/G", "A/G", "G/T", "C/A", "T/G"])),
                Arc::new(Int8Array::from(vec![1i8, 0, 0, 0, 1])),
            ],
        )
        .unwrap();
        let var_path = dir.path().join("var.parquet");
        let file = std::fs::File::create(&var_path).unwrap();
        let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();

        let ctx = SessionContext::new();
        let plug = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("end", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("demo_score", DataType::Float32, true),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 300])),
                Arc::new(UInt32Array::from(vec![100u32, 300])),
                Arc::new(StringArray::from(vec!["A/G", "G/A"])),
                Arc::new(Float32Array::from(vec![Some(0.9f32), Some(0.7)])),
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_demo_norm", plug).unwrap();

        let mut stream = tiered_stream(&ctx, "plugin_demo_norm", &var_path)
            .await
            .unwrap();
        let mut rows: Vec<(u32, i8)> = Vec::new();
        while let Some(b) = stream.next().await {
            let b = b.unwrap();
            let starts = b
                .column(b.schema().index_of("start").unwrap())
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();
            let tiers = b
                .column(b.schema().index_of("tier").unwrap())
                .as_any()
                .downcast_ref::<Int8Array>()
                .unwrap();
            for i in 0..b.num_rows() {
                rows.push((starts.value(i), tiers.value(i)));
            }
        }
        rows.sort();
        // 100: conflicting-tier duplicate collapses to warm (0), matching
        // pre-semi-join behavior exactly. 300: no variation match -> cold (1).
        // The three irrelevant variation rows (999/1000/1001) never surface
        // here since only the plugin's own rows are projected -- their
        // absence from `plugin_variation_raw`'s post-semi-join scan is the
        // thing under test, not something directly observable in the output.
        assert_eq!(rows, vec![(100, 0), (300, 1)]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn variation_probe_uses_true_semi_join_without_distinct_aggregate() {
        let dir = tempfile::tempdir().unwrap();
        let variation_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let variation = RecordBatch::try_new(
            variation_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 100, 999])),
                Arc::new(StringArray::from(vec!["A/G", "A/G", "C/T"])),
                Arc::new(Int8Array::from(vec![1i8, 0, 0])),
            ],
        )
        .unwrap();
        let variation_path = dir.path().join("variation.parquet");
        let file = std::fs::File::create(&variation_path).unwrap();
        let mut writer = ArrowWriter::try_new(file, variation_schema, None).unwrap();
        writer.write(&variation).unwrap();
        writer.close().unwrap();

        // Duplicate right-side keys must not multiply the variation rows.
        let key_schema = Arc::new(Schema::new(vec![
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let keys = RecordBatch::try_new(
            key_schema,
            vec![
                Arc::new(UInt32Array::from(vec![100u32, 100])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
            ],
        )
        .unwrap();
        let ctx = SessionContext::new();
        ctx.register_batch("plugin_keys", keys).unwrap();
        register_variation_probe_view(&ctx, "plugin_keys", &variation_path)
            .await
            .unwrap();

        let df = ctx
            .sql("SELECT * FROM plugin_variation_probe")
            .await
            .unwrap();
        let logical = df.logical_plan().display_indent().to_string();
        assert!(
            logical.contains("LeftSemi Join"),
            "unexpected plan:\n{logical}"
        );
        assert_eq!(
            logical.matches("Aggregate:").count(),
            1,
            "the MIN(tier) aggregation should be the only aggregate:\n{logical}"
        );

        let batches = df.collect().await.unwrap();
        assert_eq!(batches.iter().map(RecordBatch::num_rows).sum::<usize>(), 1);
        let batch = &batches[0];
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
        assert_eq!((starts.value(0), tiers.value(0)), (100, 0));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn inherits_tier_from_variation_and_miss_is_cold() {
        let ctx = SessionContext::new();

        // variation probe now carries `tier` directly (0=warm, 1=cold).
        let var = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("tier", DataType::Int8, false),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Int8Array::from(vec![0i8, 1])), // warm, cold
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_variation_probe", var).unwrap();

        // plugin ingest: warm match, cold match, and a miss.
        let plug = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("end", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("demo_score", DataType::Float32, true),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200, 300])),
                Arc::new(UInt32Array::from(vec![100u32, 200, 300])),
                Arc::new(StringArray::from(vec!["A/G", "C/T", "G/A"])),
                Arc::new(Float32Array::from(vec![Some(0.9f32), Some(0.1), Some(0.7)])),
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_demo_norm", plug).unwrap();

        let df = ctx
            .sql(&tier_sql("plugin_demo_norm", "plugin_variation_probe"))
            .await
            .unwrap()
            .sort(vec![col("start").sort(true, true)])
            .unwrap();
        let batches = df.collect().await.unwrap();
        let b = &batches[0];
        let tier = b
            .column(b.schema().index_of("tier").unwrap())
            .as_any()
            .downcast_ref::<Int8Array>()
            .unwrap();
        assert_eq!(tier.value(0), 0); // 100 -> warm (inherited)
        assert_eq!(tier.value(1), 1); // 200 -> cold (inherited)
        assert_eq!(tier.value(2), 1); // 300 -> no match -> cold
        // value column preserved, variation columns dropped
        assert!(b.schema().index_of("demo_score").is_ok());
        assert!(b.schema().index_of("tier").is_ok());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn owned_probe_chunks_are_not_resplit_or_overcounted() {
        use datafusion::common::utils::memory::get_record_batch_memory_size;
        use datafusion::datasource::MemTable;
        use datafusion::execution::memory_pool::FairSpillPool;
        use datafusion::execution::runtime_env::RuntimeEnvBuilder;
        use datafusion::physical_plan::displayable;

        const EXECUTION_BATCH_ROWS: usize = 8192;
        const PARENT_ROWS: usize = 1024 * 1024;

        let schema = Arc::new(Schema::new(vec![Field::new(
            "value",
            DataType::UInt32,
            false,
        )]));
        let parent = RecordBatch::try_new(
            schema.clone(),
            vec![Arc::new(UInt32Array::from_iter_values(
                (0..PARENT_ROWS).map(|value| value as u32),
            ))],
        )
        .unwrap();
        let parent_size = get_record_batch_memory_size(&parent);

        // Model aggregate output: execution-sized zero-copy slices retaining
        // one large parent. Rechunking must copy even though each output group
        // consists of exactly one input slice.
        let shared_chunks = (0..PARENT_ROWS)
            .step_by(EXECUTION_BATCH_ROWS)
            .map(|offset| parent.slice(offset, EXECUTION_BATCH_ROWS))
            .collect();
        let owned_chunks =
            rechunk_probe_owned(&schema, shared_chunks, EXECUTION_BATCH_ROWS).unwrap();
        assert_eq!(owned_chunks.len(), PARENT_ROWS / EXECUTION_BATCH_ROWS);
        assert!(
            owned_chunks
                .iter()
                .all(|batch| batch.num_rows() <= EXECUTION_BATCH_ROWS)
        );
        assert_eq!(
            owned_chunks
                .iter()
                .map(get_record_batch_memory_size)
                .sum::<usize>(),
            parent_size
        );

        // DataSourceExec must emit the owned batches unchanged and retain 1x
        // accounting, rather than turning large parents into shared slices.
        let ctx = SessionContext::new();
        ctx.register_table(
            "owned_probe",
            Arc::new(MemTable::try_new(schema.clone(), vec![owned_chunks.clone()]).unwrap()),
        )
        .unwrap();
        let executed = ctx
            .sql("SELECT * FROM owned_probe")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        assert_eq!(executed.len(), PARENT_ROWS / EXECUTION_BATCH_ROWS);
        assert_eq!(
            executed
                .iter()
                .map(get_record_batch_memory_size)
                .sum::<usize>(),
            parent_size
        );

        // Confirm the exact CADD plan shape now fits a bounded pool far below
        // the old 128x (512 MiB) reservation.
        let runtime = RuntimeEnvBuilder::new()
            .with_memory_pool(Arc::new(FairSpillPool::new(64 * 1024 * 1024)))
            .build_arc()
            .unwrap();
        let bounded_ctx = SessionContext::new_with_config_rt(
            SessionConfig::new().with_target_partitions(1),
            runtime,
        );
        bounded_ctx
            .register_table(
                "owned_probe",
                Arc::new(MemTable::try_new(schema.clone(), vec![owned_chunks]).unwrap()),
            )
            .unwrap();
        bounded_ctx
            .register_batch(
                "plugin",
                RecordBatch::try_new(
                    schema,
                    vec![Arc::new(UInt32Array::from_iter_values(
                        (0..(2 * PARENT_ROWS)).map(|value| value as u32),
                    ))],
                )
                .unwrap(),
            )
            .unwrap();
        let joined = bounded_ctx
            .sql(
                "SELECT p.value FROM plugin p \
                 LEFT JOIN owned_probe o ON p.value = o.value",
            )
            .await
            .unwrap();
        let plan = joined.clone().create_physical_plan().await.unwrap();
        assert!(
            displayable(plan.as_ref())
                .indent(true)
                .to_string()
                .contains("mode=CollectLeft, join_type=Right")
        );
        assert_eq!(
            joined
                .collect()
                .await
                .unwrap()
                .iter()
                .map(RecordBatch::num_rows)
                .sum::<usize>(),
            2 * PARENT_ROWS
        );
    }
}
