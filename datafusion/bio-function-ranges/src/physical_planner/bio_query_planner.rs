use async_trait::async_trait;
use datafusion::execution::context::{QueryPlanner, SessionState};
use datafusion::logical_expr::LogicalPlan;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_planner::PhysicalPlanner;
use log::info;
use std::sync::Arc;

use super::bio_physical_planner::BioPhysicalPlanner;

#[derive(Debug, Default)]
pub struct BioQueryPlanner;

#[async_trait]
impl QueryPlanner for BioQueryPlanner {
    async fn create_physical_plan(
        &self,
        logical_plan: &LogicalPlan,
        session_state: &SessionState,
    ) -> datafusion::common::Result<Arc<dyn ExecutionPlan>> {
        let physical_planner = BioPhysicalPlanner::default();
        let display_string = logical_plan.display();
        info!("Creating physical plan for logical plan: {display_string}");
        physical_planner
            .create_physical_plan(logical_plan, session_state)
            .await
    }
}
