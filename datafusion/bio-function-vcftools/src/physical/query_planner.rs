//! QueryPlanner for VCF tools that enables the FusedArrayTransform optimization.
//!
//! This module provides a custom QueryPlanner that uses a physical planner
//! configured with the FusedArrayTransformPlanner extension.

use async_trait::async_trait;
use datafusion::execution::context::{QueryPlanner, SessionState};
use datafusion::logical_expr::LogicalPlan;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_planner::{DefaultPhysicalPlanner, PhysicalPlanner};
use log::info;
use std::sync::Arc;

use super::extension_planner::FusedArrayTransformPlanner;

/// A custom QueryPlanner that uses a physical planner configured with
/// the FusedArrayTransformPlanner extension for handling FusedArrayTransform nodes.
#[derive(Debug, Default)]
pub struct VcfQueryPlanner;

impl VcfQueryPlanner {
    /// Create a new VcfQueryPlanner.
    pub fn new() -> Self {
        Self
    }
}

#[async_trait]
impl QueryPlanner for VcfQueryPlanner {
    async fn create_physical_plan(
        &self,
        logical_plan: &LogicalPlan,
        session_state: &SessionState,
    ) -> datafusion::common::Result<Arc<dyn ExecutionPlan>> {
        // Create a physical planner with our extension planner for FusedArrayTransform
        let physical_planner = DefaultPhysicalPlanner::with_extension_planners(vec![Arc::new(
            FusedArrayTransformPlanner::new(),
        )]);

        let display_string = logical_plan.display();
        info!("VcfQueryPlanner: Creating physical plan for logical plan: {display_string}");

        physical_planner
            .create_physical_plan(logical_plan, session_state)
            .await
    }
}
