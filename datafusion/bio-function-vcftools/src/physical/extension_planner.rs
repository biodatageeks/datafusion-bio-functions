//! Extension planner for FusedArrayTransform logical node.
//!
//! This planner converts the `FusedArrayTransform` logical node into
//! its physical counterpart `FusedArrayTransformExec`.

use std::sync::Arc;

use async_trait::async_trait;
use datafusion::common::{DFSchema, Result};
use datafusion::execution::SessionState;
use datafusion::logical_expr::{LogicalPlan, UserDefinedLogicalNode};
use datafusion::physical_expr::PhysicalExpr;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_planner::{ExtensionPlanner, PhysicalPlanner};

use crate::common::build_element_schema_from_arrow;
use crate::logical::fused_array_transform::FusedArrayTransform;
use crate::physical::fused_array_transform_exec::FusedArrayTransformExec;

/// Extension planner that converts `FusedArrayTransform` logical nodes
/// to `FusedArrayTransformExec` physical nodes.
#[derive(Debug, Default)]
pub struct FusedArrayTransformPlanner;

impl FusedArrayTransformPlanner {
    /// Create a new FusedArrayTransformPlanner.
    pub fn new() -> Self {
        Self
    }
}

#[async_trait]
impl ExtensionPlanner for FusedArrayTransformPlanner {
    async fn plan_extension(
        &self,
        planner: &dyn PhysicalPlanner,
        node: &dyn UserDefinedLogicalNode,
        _logical_inputs: &[&LogicalPlan],
        physical_inputs: &[Arc<dyn ExecutionPlan>],
        session_state: &SessionState,
    ) -> Result<Option<Arc<dyn ExecutionPlan>>> {
        // Check if this is our node type
        let Some(fused_transform) = node.as_any().downcast_ref::<FusedArrayTransform>() else {
            return Ok(None);
        };

        // We expect exactly one input
        if physical_inputs.len() != 1 {
            return Err(datafusion::common::DataFusionError::Plan(
                "FusedArrayTransform requires exactly one physical input".to_string(),
            ));
        }

        let input = physical_inputs[0].clone();
        let input_schema = input.schema();

        // Build the "element-level" schema for expression conversion
        // This schema represents the columns as if they were unnested (scalar element types)
        let element_schema =
            build_element_schema_from_arrow(&input_schema, fused_transform.array_columns())?;
        let element_df_schema = DFSchema::try_from(element_schema)?;

        // Convert logical expressions to physical expressions using the element schema
        let physical_exprs: Vec<Arc<dyn PhysicalExpr>> = fused_transform
            .transform_exprs()
            .iter()
            .map(|expr| planner.create_physical_expr(expr, &element_df_schema, session_state))
            .collect::<Result<Vec<_>>>()?;

        // Create the physical execution plan
        let exec = FusedArrayTransformExec::try_new(
            input,
            fused_transform.array_columns().to_vec(),
            fused_transform.passthrough_columns().to_vec(),
            fused_transform.output_columns().to_vec(),
            physical_exprs,
        )?;

        Ok(Some(Arc::new(exec)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planner_creation() {
        let planner = FusedArrayTransformPlanner::new();
        assert!(format!("{planner:?}").contains("FusedArrayTransformPlanner"));
    }
}
