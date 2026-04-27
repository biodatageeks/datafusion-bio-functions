use crate::physical_planner::BioQueryPlanner;
use crate::physical_planner::IntervalJoinPhysicalOptimizationRule;
use async_trait::async_trait;
use datafusion::common::extensions_options;
use datafusion::config::{ConfigExtension, ConfigField, ConfigOptions};
use datafusion::execution::SessionStateBuilder;
use datafusion::execution::runtime_env::RuntimeEnv;
use datafusion::physical_optimizer::optimizer::PhysicalOptimizer;
use datafusion::prelude::{SessionConfig, SessionContext};
use log::info;
use std::str::FromStr;
use std::sync::Arc;

/// Extension trait for [`SessionContext`] that adds bio-specific functionality.
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
            .with_query_planner(Arc::new(BioQueryPlanner))
            .with_physical_optimizer_rules(rules)
            .build()
            .into();

        info!("Initialized BioQueryPlanner...");

        ctx
    }
}

extensions_options! {
    pub struct BioConfig {
        pub prefer_interval_join: bool, default = true
        pub interval_join_algorithm: Algorithm, default = Algorithm::default()
        pub interval_join_low_memory: bool, default = false
        pub count_overlaps_backend: CountOverlapsBackendMode, default = CountOverlapsBackendMode::default()
    }
}

impl ConfigExtension for BioConfig {
    const PREFIX: &'static str = "bio";
}

#[derive(Debug, Eq, PartialEq, Default, Clone, Copy)]
pub enum CountOverlapsBackendMode {
    #[default]
    Auto,
    Cpu,
    AppleGpu,
}

#[derive(Debug)]
pub struct ParseCountOverlapsBackendModeError(String);

impl std::fmt::Display for ParseCountOverlapsBackendModeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ParseCountOverlapsBackendModeError {}

impl FromStr for CountOverlapsBackendMode {
    type Err = ParseCountOverlapsBackendModeError;

    #[inline]
    fn from_str(s: &str) -> Result<CountOverlapsBackendMode, Self::Err> {
        match s.to_lowercase().replace(['-', '_'], "").as_str() {
            "auto" => Ok(CountOverlapsBackendMode::Auto),
            "cpu" => Ok(CountOverlapsBackendMode::Cpu),
            "applegpu" => Ok(CountOverlapsBackendMode::AppleGpu),
            _ => Err(ParseCountOverlapsBackendModeError(format!(
                "Can't parse '{s}' as CountOverlapsBackendMode"
            ))),
        }
    }
}

impl std::fmt::Display for CountOverlapsBackendMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let val = match self {
            CountOverlapsBackendMode::Auto => "auto",
            CountOverlapsBackendMode::Cpu => "cpu",
            CountOverlapsBackendMode::AppleGpu => "apple_gpu",
        };
        write!(f, "{val}")
    }
}

impl From<ParseCountOverlapsBackendModeError> for datafusion::error::DataFusionError {
    fn from(e: ParseCountOverlapsBackendModeError) -> Self {
        datafusion::error::DataFusionError::External(Box::new(e))
    }
}

impl ConfigField for CountOverlapsBackendMode {
    fn set(&mut self, _key: &str, value: &str) -> datafusion::common::Result<()> {
        *self = value.parse::<CountOverlapsBackendMode>()?;
        Ok(())
    }

    fn visit<V: datafusion::config::Visit>(&self, visitor: &mut V, name: &str, doc: &'static str) {
        visitor.some(name, self, doc)
    }
}

#[derive(Debug, Eq, PartialEq, Default, Clone, Copy)]
pub enum Algorithm {
    #[default]
    Coitrees,
    IntervalTree,
    ArrayIntervalTree,
    Lapper,
    SuperIntervals,
    CoitreesNearest,
    CoitreesCountOverlaps,
}

#[derive(Debug)]
pub struct ParseAlgorithmError(String);

impl std::fmt::Display for ParseAlgorithmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ParseAlgorithmError {}

impl FromStr for Algorithm {
    type Err = ParseAlgorithmError;

    #[inline]
    fn from_str(s: &str) -> Result<Algorithm, Self::Err> {
        match s.to_lowercase().as_str() {
            "coitrees" => Ok(Algorithm::Coitrees),
            "intervaltree" => Ok(Algorithm::IntervalTree),
            "arrayintervaltree" => Ok(Algorithm::ArrayIntervalTree),
            "lapper" => Ok(Algorithm::Lapper),
            "superintervals" => Ok(Algorithm::SuperIntervals),
            "coitreesnearest" => Ok(Algorithm::CoitreesNearest),
            "coitreescountoverlaps" => Ok(Algorithm::CoitreesCountOverlaps),
            _ => Err(ParseAlgorithmError(format!(
                "Can't parse '{s}' as Algorithm"
            ))),
        }
    }
}

impl std::fmt::Display for Algorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let val = match self {
            Algorithm::Coitrees => "Coitrees",
            Algorithm::IntervalTree => "IntervalTree",
            Algorithm::ArrayIntervalTree => "ArrayIntervalTree",
            Algorithm::Lapper => "Lapper",
            Algorithm::SuperIntervals => "SuperIntervals",
            Algorithm::CoitreesNearest => "CoitreesNearest",
            Algorithm::CoitreesCountOverlaps => "CoitreesCountOverlaps",
        };
        write!(f, "{val}")
    }
}

impl From<ParseAlgorithmError> for datafusion::error::DataFusionError {
    fn from(e: ParseAlgorithmError) -> Self {
        datafusion::error::DataFusionError::External(Box::new(e))
    }
}

impl ConfigField for Algorithm {
    fn set(&mut self, _key: &str, value: &str) -> datafusion::common::Result<()> {
        *self = value.parse::<Algorithm>()?;
        Ok(())
    }

    fn visit<V: datafusion::config::Visit>(&self, visitor: &mut V, name: &str, doc: &'static str) {
        visitor.some(name, self, doc)
    }
}

/// Convenience function: create a [`SessionContext`] with interval join
/// optimization and range table functions (`coverage`, `count_overlaps`).
///
/// Equivalent to calling [`BioSessionExt::new_with_bio`] followed by
/// [`register_ranges_functions`](crate::table_function::register_ranges_functions).
pub fn create_bio_session() -> SessionContext {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false);
    let ctx = SessionContext::new_with_bio(config);
    crate::table_function::register_ranges_functions(&ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_overlaps_backend_mode_parses_supported_values() {
        assert_eq!(
            "auto".parse::<CountOverlapsBackendMode>().unwrap(),
            CountOverlapsBackendMode::Auto
        );
        assert_eq!(
            "cpu".parse::<CountOverlapsBackendMode>().unwrap(),
            CountOverlapsBackendMode::Cpu
        );
        assert_eq!(
            "apple_gpu".parse::<CountOverlapsBackendMode>().unwrap(),
            CountOverlapsBackendMode::AppleGpu
        );
        assert_eq!(
            "apple-gpu".parse::<CountOverlapsBackendMode>().unwrap(),
            CountOverlapsBackendMode::AppleGpu
        );
    }

    #[test]
    fn count_overlaps_backend_mode_rejects_unknown_values() {
        assert!("metal".parse::<CountOverlapsBackendMode>().is_err());
    }
}
