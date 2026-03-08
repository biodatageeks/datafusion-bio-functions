//! Physical execution components for fused array transform.
//!
//! This module contains the physical execution plan and extension planner
//! for the fused array transform optimization.

pub mod extension_planner;
pub mod fused_array_transform_exec;
pub mod query_planner;

pub use extension_planner::FusedArrayTransformPlanner;
pub use fused_array_transform_exec::FusedArrayTransformExec;
pub use query_planner::VcfQueryPlanner;
