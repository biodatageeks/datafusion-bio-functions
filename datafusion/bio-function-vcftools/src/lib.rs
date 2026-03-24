//! Fused Array Transform Optimization for Apache DataFusion
//!
//! This crate provides a memory-efficient optimization for queries that follow
//! the pattern: `unnest → transform → array_agg`. Instead of materializing all
//! unnested rows (which can cause memory explosion), this optimization processes
//! rows in a streaming fashion with bounded memory.
//!
//! # Problem
//!
//! Queries like:
//! ```sql
//! SELECT row_idx, array_agg(val * 2)
//! FROM (SELECT unnest(values) as val, ... FROM input)
//! GROUP BY row_idx
//! ```
//!
//! Can consume excessive memory because the GROUP BY materializes all unnested rows
//! before aggregating, even though the grouping perfectly partitions the data back
//! to original rows.
//!
//! # Solution
//!
//! This crate provides the following components:
//!
//! - [`FusedArrayTransform`]: A custom logical node representing the fused operation
//! - [`FusedArrayTransformOptimizerRule`]: A logical optimizer rule that detects the pattern
//! - [`FusedArrayTransformPlanner`]: An extension planner that converts the logical node to physical
//! - [`FusedArrayTransformExec`]: A streaming physical operator
//!
//! # Usage
//!
//! ```ignore
//! use datafusion::prelude::*;
//! use datafusion_bio_function_vcftools::{
//!     FusedArrayTransformOptimizerRule,
//!     FusedArrayTransformPlanner,
//!     enable_fused_array_transform,
//! };
//!
//! // Enable the optimization (thread-safe)
//! enable_fused_array_transform();
//!
//! let ctx = SessionContext::new();
//!
//! // Register the logical optimizer rule
//! ctx.add_optimizer_rule(Arc::new(FusedArrayTransformOptimizerRule::new()));
//!
//! // The optimizer will automatically detect and optimize matching patterns
//! ```

pub mod common;
pub mod logical;
pub mod physical;
pub mod udfs;

// Re-export logical types
pub use logical::fused_array_transform::FusedArrayTransform;
pub use logical::optimizer_rule::{
    FusedArrayTransformOptimizerRule, disable_fused_array_transform, enable_fused_array_transform,
    is_fused_array_transform_enabled, set_fused_array_transform_enabled,
};

// Re-export physical types
pub use physical::extension_planner::FusedArrayTransformPlanner;
pub use physical::fused_array_transform_exec::FusedArrayTransformExec;
pub use physical::query_planner::VcfQueryPlanner;

// Re-export UDF registration
pub use udfs::register_list_udfs;
