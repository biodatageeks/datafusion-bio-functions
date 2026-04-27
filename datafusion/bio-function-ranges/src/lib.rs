pub mod array_utils;
pub mod cluster;
pub mod complement;
pub mod count_overlaps;
#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
mod count_overlaps_apple_gpu;
pub mod count_overlaps_rank;
pub mod filter_op;
pub mod grouped_stream;
pub mod interval_tree;
pub mod merge;
pub mod nearest;
pub mod nearest_index;
pub mod overlap;
pub mod physical_planner;
pub mod session_context;
pub mod subtract;
pub mod table_function;

// Re-export key types
pub use cluster::ClusterProvider;
pub use complement::ComplementProvider;
pub use count_overlaps::{
    CountOverlapsProvider, CountOverlapsTimingSnapshot, reset_count_overlaps_timings,
    snapshot_count_overlaps_timings,
};
#[cfg(all(feature = "apple-gpu", target_os = "macos"))]
pub use count_overlaps_apple_gpu::{
    AppleGpuTimingSnapshot, reset_apple_gpu_timings, snapshot_apple_gpu_timings,
};
pub use filter_op::FilterOp;
pub use merge::MergeProvider;
pub use nearest::NearestProvider;
pub use overlap::OverlapProvider;
pub use physical_planner::BioQueryPlanner;
pub use physical_planner::IntervalJoinPhysicalOptimizationRule;
pub use physical_planner::joins::interval_join::IntervalJoinExec;
pub use session_context::{
    Algorithm, BioConfig, BioSessionExt, CountOverlapsBackendMode, create_bio_session,
};
pub use subtract::SubtractProvider;
pub use table_function::register_ranges_functions;
