//! Pool-aware choice between DataFusion's built-in hash and sort-merge joins.
//!
//! DataFusion 53 uses statistics to choose a hash-join build side and
//! partition mode, but its HashJoin-vs-SortMergeJoin choice is controlled only
//! by `prefer_hash_join`. This module supplies the missing memory-feasibility
//! check. It does not implement either join.

use std::mem::size_of;
use std::sync::Arc;

use datafusion::common::Result;
use datafusion::common::utils::memory::estimate_memory_size;
use datafusion::execution::memory_pool::{MemoryLimit, MemoryPool};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::joins::join_hash_map::{JoinHashMapU32, JoinHashMapU64};
use datafusion::physical_plan::joins::{HashJoinExec, PartitionMode, SortMergeJoinExec};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum JoinAlgorithm {
    Hash,
    SortMerge,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum JoinDecisionReason {
    BuildFits,
    BuildExceedsPool,
    MissingStatistics,
    UnknownPoolLimit,
    UnexpectedHashMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct JoinDecision {
    pub algorithm: JoinAlgorithm,
    pub reason: JoinDecisionReason,
    pub build_rows: Option<usize>,
    pub build_data_bytes: Option<usize>,
    pub hash_build_bytes: Option<usize>,
    pub available_bytes: Option<usize>,
}

fn find_hash_join(plan: &dyn ExecutionPlan) -> Option<&HashJoinExec> {
    if let Some(join) = plan.as_any().downcast_ref::<HashJoinExec>() {
        return Some(join);
    }
    plan.children()
        .into_iter()
        .find_map(|child| find_hash_join(child.as_ref()))
}

pub(crate) fn contains_sort_merge_join(plan: &dyn ExecutionPlan) -> bool {
    plan.as_any().is::<SortMergeJoinExec>()
        || plan
            .children()
            .into_iter()
            .any(|child| contains_sort_merge_join(child.as_ref()))
}

fn hash_table_bytes(rows: usize) -> Result<usize> {
    if rows > u32::MAX as usize {
        estimate_memory_size::<(u64, u64)>(rows, size_of::<JoinHashMapU64>())
    } else {
        estimate_memory_size::<(u32, u64)>(rows, size_of::<JoinHashMapU32>())
    }
}

/// Decide whether the actual build side in an optimized HashJoin plan fits the
/// active pool. The estimate mirrors the reservations made by DataFusion 53's
/// `collect_left_input`: input batches, hash table, and (conservatively) an
/// outer-join visited bitmap.
pub(crate) fn choose_for_hash_plan(
    plan: &Arc<dyn ExecutionPlan>,
    pool: &dyn MemoryPool,
) -> Result<JoinDecision> {
    let Some(join) = find_hash_join(plan.as_ref()) else {
        return Ok(JoinDecision {
            algorithm: JoinAlgorithm::SortMerge,
            reason: JoinDecisionReason::MissingStatistics,
            build_rows: None,
            build_data_bytes: None,
            hash_build_bytes: None,
            available_bytes: None,
        });
    };

    if join.partition_mode() != &PartitionMode::CollectLeft {
        return Ok(JoinDecision {
            algorithm: JoinAlgorithm::SortMerge,
            reason: JoinDecisionReason::UnexpectedHashMode,
            build_rows: None,
            build_data_bytes: None,
            hash_build_bytes: None,
            available_bytes: None,
        });
    }

    let stats = join.left().partition_statistics(None)?;
    let rows = stats.num_rows.get_value().copied();
    let data_bytes = stats.total_byte_size.get_value().copied();
    let (Some(rows), Some(data_bytes)) = (rows, data_bytes) else {
        return Ok(JoinDecision {
            algorithm: JoinAlgorithm::SortMerge,
            reason: JoinDecisionReason::MissingStatistics,
            build_rows: rows,
            build_data_bytes: data_bytes,
            hash_build_bytes: None,
            available_bytes: None,
        });
    };

    let build_bytes = data_bytes
        .checked_add(hash_table_bytes(rows)?)
        .and_then(|bytes| bytes.checked_add(rows.div_ceil(8)))
        .ok_or_else(|| {
            datafusion::common::DataFusionError::Execution(
                "hash join build-size estimate overflowed usize".into(),
            )
        })?;

    let available = match pool.memory_limit() {
        MemoryLimit::Finite(limit) => Some(limit.saturating_sub(pool.reserved())),
        MemoryLimit::Infinite => Some(usize::MAX),
        MemoryLimit::Unknown => None,
    };
    let Some(available) = available else {
        return Ok(JoinDecision {
            algorithm: JoinAlgorithm::SortMerge,
            reason: JoinDecisionReason::UnknownPoolLimit,
            build_rows: Some(rows),
            build_data_bytes: Some(data_bytes),
            hash_build_bytes: Some(build_bytes),
            available_bytes: None,
        });
    };

    let (algorithm, reason) = if build_bytes <= available {
        (JoinAlgorithm::Hash, JoinDecisionReason::BuildFits)
    } else {
        (
            JoinAlgorithm::SortMerge,
            JoinDecisionReason::BuildExceedsPool,
        )
    };
    Ok(JoinDecision {
        algorithm,
        reason,
        build_rows: Some(rows),
        build_data_bytes: Some(data_bytes),
        hash_build_bytes: Some(build_bytes),
        available_bytes: Some(available),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::UInt32Array;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use datafusion::execution::memory_pool::FairSpillPool;
    use datafusion::prelude::{SessionConfig, SessionContext};

    async fn hash_plan(build_rows: usize) -> Arc<dyn ExecutionPlan> {
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
        for (name, rows) in [("build", build_rows), ("probe", build_rows * 2)] {
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
        ctx.sql("SELECT p.value FROM probe p INNER JOIN build b ON p.value = b.value")
            .await
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap()
    }

    #[tokio::test]
    async fn pool_size_switches_hash_to_sort_merge() {
        let plan = hash_plan(128).await;
        let roomy = FairSpillPool::new(1024 * 1024);
        let tight = FairSpillPool::new(64);
        assert_eq!(
            choose_for_hash_plan(&plan, &roomy).unwrap().algorithm,
            JoinAlgorithm::Hash
        );
        assert_eq!(
            choose_for_hash_plan(&plan, &tight).unwrap().algorithm,
            JoinAlgorithm::SortMerge
        );
    }

    #[tokio::test]
    async fn build_statistics_switch_algorithm_for_the_same_pool() {
        let small_plan = hash_plan(16).await;
        let large_plan = hash_plan(4096).await;
        let measuring_pool = FairSpillPool::new(1024 * 1024);
        let small_bytes = choose_for_hash_plan(&small_plan, &measuring_pool)
            .unwrap()
            .hash_build_bytes
            .unwrap();
        let large_bytes = choose_for_hash_plan(&large_plan, &measuring_pool)
            .unwrap()
            .hash_build_bytes
            .unwrap();
        assert!(small_bytes < large_bytes);

        let same_pool = FairSpillPool::new(small_bytes + (large_bytes - small_bytes) / 2);
        assert_eq!(
            choose_for_hash_plan(&small_plan, &same_pool)
                .unwrap()
                .algorithm,
            JoinAlgorithm::Hash
        );
        assert_eq!(
            choose_for_hash_plan(&large_plan, &same_pool)
                .unwrap()
                .algorithm,
            JoinAlgorithm::SortMerge
        );
    }
}
