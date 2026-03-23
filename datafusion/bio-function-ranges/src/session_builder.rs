use datafusion::config::ConfigOptions;
use datafusion::execution::SessionStateBuilder;
use datafusion::execution::runtime_env::RuntimeEnv;
use datafusion::physical_optimizer::optimizer::PhysicalOptimizer;
use datafusion::prelude::{SessionConfig, SessionContext};
use std::sync::Arc;

use crate::config::BioConfig;
use crate::physical_planner::{IntervalJoinPhysicalOptimizationRule, RangesQueryPlanner};

/// Build a ranges-optimized session from an existing SessionConfig.
pub fn create_ranges_session_with_config(config: SessionConfig) -> SessionContext {
    let runtime = Arc::new(RuntimeEnv::default());
    let mut rules = PhysicalOptimizer::new().rules;
    rules.retain(|rule| rule.name() != "join_selection");
    rules.push(Arc::new(IntervalJoinPhysicalOptimizationRule));

    SessionStateBuilder::new()
        .with_config(config)
        .with_runtime_env(runtime)
        .with_default_features()
        .with_query_planner(Arc::new(RangesQueryPlanner::new()))
        .with_physical_optimizer_rules(rules)
        .build()
        .into()
}

/// Create a ranges-optimized session with default BioConfig and defaults used by tests/examples.
pub fn create_ranges_session() -> SessionContext {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false);
    create_ranges_session_with_config(config)
}
