//! Logical plan components for fused array transform optimization.
//!
//! This module contains the custom logical node and optimizer rule that
//! detects and rewrites the `unnest → transform → array_agg` pattern.

pub mod fused_array_transform;
pub mod optimizer_rule;

pub use fused_array_transform::FusedArrayTransform;
pub use optimizer_rule::{
    FusedArrayTransformOptimizerRule, disable_fused_array_transform, enable_fused_array_transform,
    is_fused_array_transform_enabled, set_fused_array_transform_enabled,
};
