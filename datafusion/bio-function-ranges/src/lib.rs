pub mod array_utils;
pub mod count_overlaps;
pub mod filter_op;
pub mod interval_tree;
pub mod physical_planner;
pub mod session_context;
pub mod table_function;

// Re-export key types
pub use count_overlaps::CountOverlapsProvider;
pub use filter_op::FilterOp;
pub use physical_planner::BioQueryPlanner;
pub use physical_planner::IntervalJoinPhysicalOptimizationRule;
pub use physical_planner::joins::interval_join::IntervalJoinExec;
pub use session_context::{Algorithm, BioConfig, BioSessionExt, create_bio_session};
pub use table_function::register_ranges_functions;
