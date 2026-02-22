use crate::physical_planner::intervals::{ColIntervals, parse};
use crate::physical_planner::joins::interval_join::IntervalJoinExec;
use crate::session_context::{Algorithm, BioConfig};
use async_trait::async_trait;
use datafusion::common::tree_node::{Transformed, TransformedResult, TreeNode};
use datafusion::common::{DFSchema, NullEquality, Result};
use datafusion::config::ConfigOptions;
use datafusion::execution::context::SessionState;
use datafusion::logical_expr::{Expr, LogicalPlan};
use datafusion::physical_expr::PhysicalExpr;
use datafusion::physical_expr::expressions::lit;
use datafusion::physical_optimizer::PhysicalOptimizerRule;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::joins::{HashJoinExec, NestedLoopJoinExec, PartitionMode};
use datafusion::physical_planner::{DefaultPhysicalPlanner, PhysicalPlanner};
use log::info;
use std::sync::Arc;

#[derive(Default)]
pub struct BioPhysicalPlanner {
    planner: DefaultPhysicalPlanner,
}

#[derive(Debug, Default)]
pub struct IntervalJoinPhysicalOptimizationRule;

impl PhysicalOptimizerRule for IntervalJoinPhysicalOptimizationRule {
    fn optimize(
        &self,
        plan: Arc<dyn ExecutionPlan>,
        config: &ConfigOptions,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let default = BioConfig::default();
        let bio_config = config.extensions.get::<BioConfig>().unwrap_or(&default);

        if !bio_config.prefer_interval_join {
            info!("bio.prefer_interval_join = false, skipping optimization");
            return Ok(plan);
        }

        let algorithm = bio_config.interval_join_algorithm;

        let low_memory = bio_config.interval_join_low_memory;

        plan.transform_up(|plan| {
            match plan.as_any().downcast_ref::<HashJoinExec>() {
                Some(join_exec) => {
                    info!("HashJoinExec detected");
                    if let Some(intervals) = parse(join_exec.filter()) {
                        info!("Detected HashJoinExec with Range filters. Optimizing into IntervalJoinExec using {algorithm} algorithm...");
                        let new_plan = from_hash_join(
                            join_exec,
                            intervals,
                            algorithm,
                            low_memory,
                        )?;
                        Ok(Transformed::yes(new_plan))
                    } else {
                        info!("Could not build range filter from {}, skipping optimization",
                                    join_exec.filter().map_or(
                                        "".to_string(),
                                        |f| format!("'{}'", f.expression())));
                        Ok(Transformed::no(plan))
                    }
                },
                None => {
                    match plan.as_any().downcast_ref::<NestedLoopJoinExec>() {
                        Some(join_exec) => {
                            info!("NestedLoopJoinExec detected");
                            if let Some(intervals) = parse(join_exec.filter()) {
                                info!("Detected NestedLoopJoinExec with Range filters. Optimizing into IntervalJoinExec using {algorithm} algorithm...");
                                let new_plan = from_nested_loop_join(
                                    join_exec,
                                    intervals,
                                    algorithm,
                                    low_memory,
                                )?;
                                Ok(Transformed::yes(new_plan))
                            } else {
                                info!("Could not build range filter from {}, skipping optimization",
                                    join_exec.filter().map_or(
                                        "".to_string(),
                                        |f| format!("'{}'", f.expression())));
                                Ok(Transformed::no(plan))
                            }
                        },
                        None => Ok(Transformed::no(plan)),
                    }
                }
            }
        }).data()
    }

    fn name(&self) -> &str {
        "IntervalJoinOptimizationRule"
    }

    fn schema_check(&self) -> bool {
        true
    }
}

fn from_hash_join(
    join_exec: &HashJoinExec,
    intervals: ColIntervals,
    algorithm: Algorithm,
    low_memory: bool,
) -> Result<Arc<dyn ExecutionPlan>> {
    let new_plan = IntervalJoinExec::try_new(
        join_exec.left().clone(),
        join_exec.right().clone(),
        join_exec.on.clone(),
        join_exec.filter.clone(),
        intervals,
        &join_exec.join_type,
        join_exec.projection.clone(),
        *join_exec.partition_mode(),
        join_exec.null_equality() == NullEquality::NullEqualsNull,
        algorithm,
        low_memory,
    )?;
    Ok(Arc::new(new_plan))
}

fn from_nested_loop_join(
    join_exec: &NestedLoopJoinExec,
    intervals: ColIntervals,
    algorithm: Algorithm,
    low_memory: bool,
) -> Result<Arc<dyn ExecutionPlan>> {
    let new_plan = IntervalJoinExec::try_new(
        join_exec.left().clone(),
        join_exec.right().clone(),
        vec![(lit(1), lit(1))],
        join_exec.filter().cloned(),
        intervals,
        join_exec.join_type(),
        join_exec.projection().cloned(),
        PartitionMode::CollectLeft,
        true,
        algorithm,
        low_memory,
    )?;

    Ok(Arc::new(new_plan))
}

#[async_trait]
impl PhysicalPlanner for BioPhysicalPlanner {
    async fn create_physical_plan(
        &self,
        logical_plan: &LogicalPlan,
        session_state: &SessionState,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let plan = self
            .planner
            .create_physical_plan(logical_plan, session_state)
            .await?;
        Ok(plan)
    }

    fn create_physical_expr(
        &self,
        expr: &Expr,
        input_dfschema: &DFSchema,
        session_state: &SessionState,
    ) -> Result<Arc<dyn PhysicalExpr>> {
        self.planner
            .create_physical_expr(expr, input_dfschema, session_state)
    }
}
