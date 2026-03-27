//! VCF Sample QC Table-Valued Function.
//!
//! A high-performance TVF that applies per-site filtering and per-sample genotype
//! QC transformations, writing output directly to VCF. This avoids the UNNEST/GROUP BY
//! overhead by operating directly on Arrow List value arrays.
//!
//! # Usage
//!
//! ```sql
//! SELECT * FROM vcf_sample_qc(
//!   'vcf_table',           -- registered table name with genotypes struct
//!   '/path/output.vcf',    -- output VCF file path
//!   'SAMPLE1,SAMPLE2,...'  -- comma-separated sample names, 'auto', or
//!                          -- path to a text file with one sample name per line
//! )
//! ```

pub mod processing;
pub mod table_function;
pub mod writer;

use datafusion::common::{DataFusionError, Result};
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;
use std::sync::Arc;

/// Configuration for VCF sample QC processing.
#[derive(Debug, Clone)]
pub struct QcConfig {
    /// Minimum site QUAL value.
    pub qual_min: f64,
    /// Minimum average GQ across samples for a site to pass.
    pub site_gq_avg_min: f64,
    /// Minimum average DP across samples for a site to pass.
    pub site_dp_avg_min: f64,
    /// Maximum average DP across samples for a site to pass.
    pub site_dp_avg_max: f64,
    /// Minimum GQ for a sample to be considered "good".
    pub sample_gq_min: f64,
    /// Minimum DP for a sample to be considered "good".
    pub sample_dp_min: f64,
    /// Maximum DP for a sample to be considered "good".
    pub sample_dp_max: f64,
}

impl Default for QcConfig {
    fn default() -> Self {
        Self {
            qual_min: 20.0,
            site_gq_avg_min: 15.0,
            site_dp_avg_min: 15.0,
            site_dp_avg_max: 150.0,
            sample_gq_min: 10.0,
            sample_dp_min: 10.0,
            sample_dp_max: 200.0,
        }
    }
}

/// Register the `vcf_sample_qc` table function on a session context.
pub fn register_sample_qc_function(ctx: &SessionContext) {
    let session = Arc::new(ctx.clone());
    ctx.register_udtf(
        "vcf_sample_qc",
        Arc::new(table_function::VcfSampleQcFunction::new(session)),
    );
}

/// Extract a string literal from an `Expr`.
pub(crate) fn extract_string_arg(arg: &Expr, name: &str) -> Result<String> {
    use datafusion::common::ScalarValue;
    match arg {
        Expr::Literal(ScalarValue::Utf8(Some(val)), _) => Ok(val.clone()),
        Expr::Literal(ScalarValue::LargeUtf8(Some(val)), _) => Ok(val.clone()),
        Expr::Literal(ScalarValue::Utf8View(Some(val)), _) => Ok(val.clone()),
        other => Err(DataFusionError::Plan(format!(
            "vcf_sample_qc() {name} must be a string literal, got: {other}"
        ))),
    }
}

/// Extract a float literal from an `Expr`.
pub(crate) fn extract_f64_arg(arg: &Expr, name: &str) -> Result<f64> {
    use datafusion::common::ScalarValue;
    match arg {
        Expr::Literal(ScalarValue::Float64(Some(v)), _) => Ok(*v),
        Expr::Literal(ScalarValue::Float32(Some(v)), _) => Ok(*v as f64),
        Expr::Literal(ScalarValue::Int64(Some(v)), _) => Ok(*v as f64),
        Expr::Literal(ScalarValue::Int32(Some(v)), _) => Ok(*v as f64),
        other => Err(DataFusionError::Plan(format!(
            "vcf_sample_qc() {name} must be a numeric literal, got: {other}"
        ))),
    }
}
