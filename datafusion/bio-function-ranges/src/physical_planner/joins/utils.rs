// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

//! Join related functionality used both on logical and physical plans

use std::future::Future;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::common::stats::Precision;
use datafusion::common::{DataFusionError, JoinType, Result, SharedResult};
use datafusion::logical_expr::interval_arithmetic::Interval;
use datafusion::physical_expr::expressions::Column;
use datafusion::physical_plan::execution_plan::Boundedness;
use datafusion::physical_plan::joins::utils::{JoinOn, adjust_right_output_partitioning};
use datafusion::physical_plan::metrics;
use datafusion::physical_plan::metrics::{ExecutionPlanMetricsSet, MetricBuilder};
use datafusion::physical_plan::{ColumnStatistics, ExecutionPlan, Statistics};
use datafusion::physical_plan::{ExecutionPlanProperties, Partitioning};

pub(crate) fn boundedness_from_children<'a>(
    children: impl IntoIterator<Item = &'a Arc<dyn ExecutionPlan>>,
) -> Boundedness {
    let mut unbounded_with_finite_mem = false;
    for child in children {
        match child.boundedness() {
            Boundedness::Unbounded {
                requires_infinite_memory: true,
            } => {
                return Boundedness::Unbounded {
                    requires_infinite_memory: true,
                };
            }
            Boundedness::Unbounded {
                requires_infinite_memory: false,
            } => {
                unbounded_with_finite_mem = true;
            }
            Boundedness::Bounded => {}
        }
    }
    if unbounded_with_finite_mem {
        Boundedness::Unbounded {
            requires_infinite_memory: false,
        }
    } else {
        Boundedness::Bounded
    }
}
use futures::future::{BoxFuture, Shared};
use futures::{FutureExt, ready};
use parking_lot::Mutex;

/// A [`OnceAsync`] can be used to run an async closure once, with subsequent calls
/// to [`OnceAsync::once`] returning a [`OnceFut`] to the same asynchronous computation
///
/// This is useful for joins where the results of one child are buffered in memory
/// and shared across potentially multiple output partitions
pub(crate) struct OnceAsync<T> {
    fut: Mutex<Option<OnceFut<T>>>,
}

impl<T> Default for OnceAsync<T> {
    fn default() -> Self {
        Self {
            fut: Mutex::new(None),
        }
    }
}

impl<T> std::fmt::Debug for OnceAsync<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "OnceAsync")
    }
}

impl<T: 'static> OnceAsync<T> {
    /// If this is the first call to this function on this object, will invoke
    /// `f` to obtain a future and return a [`OnceFut`] referring to this
    ///
    /// If this is not the first call, will return a [`OnceFut`] referring
    /// to the same future as was returned by the first call
    pub(crate) fn once<F, Fut>(&self, f: F) -> OnceFut<T>
    where
        F: FnOnce() -> Fut,
        Fut: Future<Output = Result<T>> + Send + 'static,
    {
        self.fut
            .lock()
            .get_or_insert_with(|| OnceFut::new(f()))
            .clone()
    }
}

/// The shared future type used internally within [`OnceAsync`]
type OnceFutPending<T> = Shared<BoxFuture<'static, SharedResult<Arc<T>>>>;

/// A [`OnceFut`] represents a shared asynchronous computation, that will be evaluated
/// once for all [`Clone`]'s, with [`OnceFut::get`] providing a non-consuming interface
/// to drive the underlying [`Future`] to completion
pub(crate) struct OnceFut<T> {
    state: OnceFutState<T>,
}

impl<T> Clone for OnceFut<T> {
    fn clone(&self) -> Self {
        Self {
            state: self.state.clone(),
        }
    }
}

/// A shared state between statistic aggregators for a join
/// operation.
#[derive(Clone, Debug, Default)]
struct PartialJoinStatistics {
    pub num_rows: usize,
    pub column_statistics: Vec<ColumnStatistics>,
}

/// Estimate the statistics for the given join's output.
pub(crate) fn estimate_join_statistics(
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    on: JoinOn,
    join_type: &JoinType,
) -> Result<Statistics> {
    let left_stats = left.partition_statistics(None)?;
    let right_stats = right.partition_statistics(None)?;
    let schema = left.schema();
    let join_stats = estimate_join_cardinality(join_type, left_stats, right_stats, &on);
    let (num_rows, column_statistics) = match join_stats {
        Some(stats) => (Precision::Inexact(stats.num_rows), stats.column_statistics),
        None => (Precision::Absent, Statistics::unknown_column(&schema)),
    };
    Ok(Statistics {
        num_rows,
        total_byte_size: Precision::Absent,
        column_statistics,
    })
}

// Estimate the cardinality for the given join with input statistics.
fn estimate_join_cardinality(
    join_type: &JoinType,
    left_stats: Statistics,
    right_stats: Statistics,
    on: &JoinOn,
) -> Option<PartialJoinStatistics> {
    match join_type {
        JoinType::Inner
        | JoinType::Left
        | JoinType::Right
        | JoinType::Full
        | JoinType::LeftMark => {
            let (left_col_stats, right_col_stats) = on
                .iter()
                .map(|(left, right)| {
                    match (
                        left.as_any().downcast_ref::<Column>(),
                        right.as_any().downcast_ref::<Column>(),
                    ) {
                        (Some(left), Some(right)) => (
                            left_stats.column_statistics[left.index()].clone(),
                            right_stats.column_statistics[right.index()].clone(),
                        ),
                        _ => (
                            ColumnStatistics::new_unknown(),
                            ColumnStatistics::new_unknown(),
                        ),
                    }
                })
                .unzip::<_, _, Vec<_>, Vec<_>>();

            let ij_cardinality = estimate_inner_join_cardinality(
                Statistics {
                    num_rows: left_stats.num_rows,
                    total_byte_size: Precision::Absent,
                    column_statistics: left_col_stats,
                },
                Statistics {
                    num_rows: right_stats.num_rows,
                    total_byte_size: Precision::Absent,
                    column_statistics: right_col_stats,
                },
            )?;

            let cardinality = match join_type {
                JoinType::Inner => ij_cardinality,
                JoinType::Left => ij_cardinality.max(&left_stats.num_rows),
                JoinType::Right => ij_cardinality.max(&right_stats.num_rows),
                JoinType::Full => ij_cardinality
                    .max(&left_stats.num_rows)
                    .add(&ij_cardinality.max(&right_stats.num_rows))
                    .sub(&ij_cardinality),
                _ => unreachable!(),
            };

            Some(PartialJoinStatistics {
                num_rows: *cardinality.get_value()?,
                column_statistics: left_stats
                    .column_statistics
                    .into_iter()
                    .chain(right_stats.column_statistics)
                    .collect(),
            })
        }

        JoinType::LeftSemi
        | JoinType::RightSemi
        | JoinType::LeftAnti
        | JoinType::RightAnti
        | JoinType::RightMark => None,
    }
}

fn estimate_inner_join_cardinality(
    left_stats: Statistics,
    right_stats: Statistics,
) -> Option<Precision<usize>> {
    let mut join_selectivity = Precision::Absent;
    for (left_stat, right_stat) in left_stats
        .column_statistics
        .iter()
        .zip(right_stats.column_statistics.iter())
    {
        if left_stat.min_value.get_value()? > right_stat.max_value.get_value()? {
            return Some(
                if left_stat.min_value.is_exact().unwrap_or(false)
                    && right_stat.max_value.is_exact().unwrap_or(false)
                {
                    Precision::Exact(0)
                } else {
                    Precision::Inexact(0)
                },
            );
        }
        if left_stat.max_value.get_value()? < right_stat.min_value.get_value()? {
            return Some(
                if left_stat.max_value.is_exact().unwrap_or(false)
                    && right_stat.min_value.is_exact().unwrap_or(false)
                {
                    Precision::Exact(0)
                } else {
                    Precision::Inexact(0)
                },
            );
        }

        let left_max_distinct = max_distinct_count(&left_stats.num_rows, left_stat);
        let right_max_distinct = max_distinct_count(&right_stats.num_rows, right_stat);
        let max_distinct = left_max_distinct.max(&right_max_distinct);
        if max_distinct.get_value().is_some() {
            join_selectivity = max_distinct;
        }
    }

    let left_num_rows = left_stats.num_rows.get_value()?;
    let right_num_rows = right_stats.num_rows.get_value()?;
    match join_selectivity {
        Precision::Exact(value) if value > 0 => {
            Some(Precision::Exact((left_num_rows * right_num_rows) / value))
        }
        Precision::Inexact(value) if value > 0 => {
            Some(Precision::Inexact((left_num_rows * right_num_rows) / value))
        }
        _ => None,
    }
}

fn max_distinct_count(num_rows: &Precision<usize>, stats: &ColumnStatistics) -> Precision<usize> {
    match &stats.distinct_count {
        dc @ (Precision::Exact(_) | Precision::Inexact(_)) => *dc,
        _ => {
            let result = match num_rows {
                Precision::Absent => Precision::Absent,
                Precision::Inexact(count) => {
                    match count.checked_sub(*stats.null_count.get_value().unwrap_or(&0)) {
                        None => Precision::Inexact(0),
                        Some(non_null_count) => Precision::Inexact(non_null_count),
                    }
                }
                Precision::Exact(count) => {
                    let count = count - stats.null_count.get_value().unwrap_or(&0);
                    if stats.null_count.is_exact().unwrap_or(false) {
                        Precision::Exact(count)
                    } else {
                        Precision::Inexact(count)
                    }
                }
            };
            if let (Some(min), Some(max)) =
                (stats.min_value.get_value(), stats.max_value.get_value())
            {
                if let Some(range_dc) = Interval::try_new(min.clone(), max.clone())
                    .ok()
                    .and_then(|e| e.cardinality())
                {
                    let range_dc = range_dc as usize;
                    return if matches!(result, Precision::Absent)
                        || &range_dc < result.get_value().unwrap()
                    {
                        if stats.min_value.is_exact().unwrap()
                            && stats.max_value.is_exact().unwrap()
                        {
                            Precision::Exact(range_dc)
                        } else {
                            Precision::Inexact(range_dc)
                        }
                    } else {
                        result
                    };
                }
            }

            result
        }
    }
}

enum OnceFutState<T> {
    Pending(OnceFutPending<T>),
    Ready(SharedResult<Arc<T>>),
}

impl<T> Clone for OnceFutState<T> {
    fn clone(&self) -> Self {
        match self {
            Self::Pending(p) => Self::Pending(p.clone()),
            Self::Ready(r) => Self::Ready(r.clone()),
        }
    }
}

impl<T: 'static> OnceFut<T> {
    /// Create a new [`OnceFut`] from a [`Future`]
    pub(crate) fn new<Fut>(fut: Fut) -> Self
    where
        Fut: Future<Output = Result<T>> + Send + 'static,
    {
        Self {
            state: OnceFutState::Pending(
                fut.map(|res| res.map(Arc::new).map_err(Arc::new))
                    .boxed()
                    .shared(),
            ),
        }
    }

    /// Get the result of the computation if it is ready, without consuming it
    #[allow(dead_code)]
    pub(crate) fn get(&mut self, cx: &mut Context<'_>) -> Poll<Result<&T>> {
        if let OnceFutState::Pending(fut) = &mut self.state {
            let r = ready!(fut.poll_unpin(cx));
            self.state = OnceFutState::Ready(r);
        }

        match &self.state {
            OnceFutState::Pending(_) => unreachable!(),
            OnceFutState::Ready(r) => Poll::Ready(
                r.as_ref()
                    .map(|r| r.as_ref())
                    .map_err(|e| DataFusionError::External(Box::new(e.clone()))),
            ),
        }
    }

    /// Get shared reference to the result of the computation if it is ready, without consuming it
    #[allow(dead_code)]
    pub(crate) fn get_shared(&mut self, cx: &mut Context<'_>) -> Poll<Result<Arc<T>>> {
        if let OnceFutState::Pending(fut) = &mut self.state {
            let r = ready!(fut.poll_unpin(cx));
            self.state = OnceFutState::Ready(r);
        }

        match &self.state {
            OnceFutState::Pending(_) => unreachable!(),
            OnceFutState::Ready(r) => Poll::Ready(
                r.clone()
                    .map_err(|e| DataFusionError::External(Box::new(e))),
            ),
        }
    }
}

/// Metrics for build & probe joins
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub(crate) struct BuildProbeJoinMetrics {
    /// Total time for collecting build-side of join
    pub(crate) build_time: metrics::Time,
    /// Number of batches consumed by build-side
    pub(crate) build_input_batches: metrics::Count,
    /// Number of rows consumed by build-side
    pub(crate) build_input_rows: metrics::Count,
    /// Memory used by build-side in bytes
    pub(crate) build_mem_used: metrics::Gauge,
    /// Total time for joining probe-side batches to the build-side batches
    pub(crate) join_time: metrics::Time,
    /// Number of batches consumed by probe-side of this operator
    pub(crate) input_batches: metrics::Count,
    /// Number of rows consumed by probe-side this operator
    pub(crate) input_rows: metrics::Count,
    /// Number of batches produced by this operator
    pub(crate) output_batches: metrics::Count,
    /// Number of rows produced by this operator
    pub(crate) output_rows: metrics::Count,
}

impl BuildProbeJoinMetrics {
    pub fn new(partition: usize, metrics: &ExecutionPlanMetricsSet) -> Self {
        let join_time = MetricBuilder::new(metrics).subset_time("join_time", partition);

        let build_time = MetricBuilder::new(metrics).subset_time("build_time", partition);

        let build_input_batches =
            MetricBuilder::new(metrics).counter("build_input_batches", partition);

        let build_input_rows = MetricBuilder::new(metrics).counter("build_input_rows", partition);

        let build_mem_used = MetricBuilder::new(metrics).gauge("build_mem_used", partition);

        let input_batches = MetricBuilder::new(metrics).counter("input_batches", partition);

        let input_rows = MetricBuilder::new(metrics).counter("input_rows", partition);

        let output_batches = MetricBuilder::new(metrics).counter("output_batches", partition);

        let output_rows = MetricBuilder::new(metrics).output_rows(partition);

        Self {
            build_time,
            build_input_batches,
            build_input_rows,
            build_mem_used,
            join_time,
            input_batches,
            input_rows,
            output_batches,
            output_rows,
        }
    }
}
pub(crate) fn symmetric_join_output_partitioning(
    left: &Arc<dyn ExecutionPlan>,
    right: &Arc<dyn ExecutionPlan>,
    join_type: &JoinType,
) -> Result<Partitioning> {
    let left_columns_len = left.schema().fields.len();
    let left_partitioning = left.output_partitioning();
    let right_partitioning = right.output_partitioning();
    match join_type {
        JoinType::Left | JoinType::LeftSemi | JoinType::LeftAnti | JoinType::LeftMark => {
            Ok(left_partitioning.clone())
        }
        JoinType::RightSemi | JoinType::RightAnti | JoinType::RightMark => {
            Ok(right_partitioning.clone())
        }
        JoinType::Inner | JoinType::Right => {
            adjust_right_output_partitioning(right_partitioning, left_columns_len)
        }
        JoinType::Full => Ok(Partitioning::UnknownPartitioning(
            right_partitioning.partition_count(),
        )),
    }
}
