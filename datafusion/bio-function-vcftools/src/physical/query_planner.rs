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
///
/// The physical planner is initialized once and reused for all queries,
/// avoiding repeated allocations since both DefaultPhysicalPlanner and
/// FusedArrayTransformPlanner are stateless.
pub struct VcfQueryPlanner {
    physical_planner: DefaultPhysicalPlanner,
}

impl std::fmt::Debug for VcfQueryPlanner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VcfQueryPlanner").finish()
    }
}

impl Default for VcfQueryPlanner {
    fn default() -> Self {
        Self::new()
    }
}

impl VcfQueryPlanner {
    /// Create a new VcfQueryPlanner with a pre-configured physical planner.
    pub fn new() -> Self {
        let physical_planner = DefaultPhysicalPlanner::with_extension_planners(vec![Arc::new(
            FusedArrayTransformPlanner::new(),
        )]);
        Self { physical_planner }
    }
}

#[async_trait]
impl QueryPlanner for VcfQueryPlanner {
    async fn create_physical_plan(
        &self,
        logical_plan: &LogicalPlan,
        session_state: &SessionState,
    ) -> datafusion::common::Result<Arc<dyn ExecutionPlan>> {
        let display_string = logical_plan.display();
        info!("VcfQueryPlanner: Creating physical plan for logical plan: {display_string}");

        self.physical_planner
            .create_physical_plan(logical_plan, session_state)
            .await
    }
}
