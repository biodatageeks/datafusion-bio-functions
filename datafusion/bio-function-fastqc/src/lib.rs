//! Streaming FastQC quality-control modules for Apache DataFusion.
//!
//! Each module implements `QcModule` (update/merge/finalize), mirroring
//! RastQC's process_sequence/merge_from/calculate_results.

pub mod basic_stats;

use std::any::Any;

/// One row of the uniform tidy output schema.
#[derive(Debug, Clone, PartialEq)]
pub struct TidyRow {
    pub module: &'static str,
    pub label: Option<String>,
    pub position: Option<i32>,
    pub metric: String,
    pub value: Option<f64>,
    pub value_str: Option<String>,
}

impl TidyRow {
    pub fn num(module: &'static str, metric: &str, value: f64) -> Self {
        TidyRow {
            module,
            label: None,
            position: None,
            metric: metric.to_string(),
            value: Some(value),
            value_str: None,
        }
    }
    pub fn status(module: &'static str, status: &str) -> Self {
        TidyRow {
            module,
            label: None,
            position: None,
            metric: "status".to_string(),
            value: None,
            value_str: Some(status.to_string()),
        }
    }
}

/// A streaming QC module accumulator.
pub trait QcModule: Send {
    /// Stable identifier, e.g. "basic_stats".
    fn name(&self) -> &'static str;
    /// Fold one read (raw sequence bytes and raw phred+33 quality bytes).
    fn update(&mut self, seq: &[u8], qual: &[u8]);
    /// Merge another accumulator of the SAME concrete type into self.
    fn merge(&mut self, other: &dyn QcModule);
    /// Emit tidy rows for this module's result.
    fn finalize(&self, out: &mut Vec<TidyRow>);
    /// For downcasting in `merge`.
    fn as_any(&self) -> &dyn Any;
}
