use async_trait::async_trait;
use datafusion::execution::context::{QueryPlanner, SessionState};
use datafusion::logical_expr::LogicalPlan;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_planner::{DefaultPhysicalPlanner, ExtensionPlanner, PhysicalPlanner};
use log::info;
use std::sync::Arc;

/// Query planner for the unified bio context.
#[derive(Default)]
pub struct BioQueryPlanner {
    physical_planner: DefaultPhysicalPlanner,
}

impl std::fmt::Debug for BioQueryPlanner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BioQueryPlanner").finish()
    }
}

impl BioQueryPlanner {
    pub fn with_extension_planners(
        extension_planners: Vec<Arc<dyn ExtensionPlanner + Send + Sync>>,
    ) -> Self {
        Self {
            physical_planner: DefaultPhysicalPlanner::with_extension_planners(extension_planners),
        }
    }
}

#[async_trait]
impl QueryPlanner for BioQueryPlanner {
    async fn create_physical_plan(
        &self,
        logical_plan: &LogicalPlan,
        session_state: &SessionState,
    ) -> datafusion::common::Result<Arc<dyn ExecutionPlan>> {
        let display_string = logical_plan.display();
        info!("BioQueryPlanner: Creating physical plan for logical plan: {display_string}");
        self.physical_planner
            .create_physical_plan(logical_plan, session_state)
            .await
    }
}
