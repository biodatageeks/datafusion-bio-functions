use async_trait::async_trait;
use datafusion::config::ConfigOptions;
use datafusion::execution::SessionStateBuilder;
use datafusion::execution::runtime_env::RuntimeEnv;
use datafusion::physical_optimizer::optimizer::PhysicalOptimizer;
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_ranges::physical_planner::IntervalJoinPhysicalOptimizationRule;
use datafusion_bio_function_ranges::{BioConfig, register_ranges_functions};
use datafusion_bio_function_vcftools::{
    FusedArrayTransformOptimizerRule, FusedArrayTransformPlanner, enable_fused_array_transform,
};
use log::info;
use std::sync::Arc;

mod query_planner;
use query_planner::BioQueryPlanner;

/// Extension trait for [`SessionContext`] that enables all bio optimizations.
#[async_trait]
pub trait BioSessionExt {
    fn new_with_bio(config: SessionConfig) -> SessionContext;
    fn with_config_rt_bio(config: SessionConfig, runtime: Arc<RuntimeEnv>) -> SessionContext;
}

impl BioSessionExt for SessionContext {
    fn new_with_bio(config: SessionConfig) -> SessionContext {
        info!("Loading BioSessionExt...");
        let runtime = Arc::new(RuntimeEnv::default());
        Self::with_config_rt_bio(config, runtime)
    }

    fn with_config_rt_bio(config: SessionConfig, runtime: Arc<RuntimeEnv>) -> SessionContext {
        let mut rules = PhysicalOptimizer::new().rules;
        rules.retain(|rule| rule.name() != "join_selection");
        rules.push(Arc::new(IntervalJoinPhysicalOptimizationRule));

        let ctx: SessionContext = SessionStateBuilder::new()
            .with_config(config)
            .with_runtime_env(runtime)
            .with_default_features()
            .with_query_planner(Arc::new(BioQueryPlanner::with_extension_planners(vec![
                Arc::new(FusedArrayTransformPlanner::new()),
            ])))
            .with_physical_optimizer_rules(rules)
            .build()
            .into();

        ctx.add_optimizer_rule(Arc::new(FusedArrayTransformOptimizerRule::new()));
        enable_fused_array_transform();

        info!("Initialized unified Bio session context...");

        ctx
    }
}

/// Convenience function: create a [`SessionContext`] with all bio query planning
/// optimizations and range table functions (`coverage`, `count_overlaps`).
pub fn create_bio_session() -> SessionContext {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false);
    let ctx = SessionContext::new_with_bio(config);
    register_ranges_functions(&ctx);
    ctx
}
