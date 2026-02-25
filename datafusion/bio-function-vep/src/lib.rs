//! VEP (Variant Effect Predictor) annotation functions for Apache DataFusion.
//!
//! Provides:
//! - `lookup_variants()` table function for known variant lookup via interval join
//! - `match_allele()` scalar UDF for allele matching
//! - `vep_allele()` scalar UDF for VCF→VEP allele conversion

pub mod allele;
pub mod config;
pub mod coordinate;
pub mod lookup_provider;
pub mod schema_contract;
pub mod table_function;

pub use config::AnnotationConfig;
#[cfg(feature = "kv-cache")]
pub use datafusion_bio_function_vep_cache as kv_cache;

use std::sync::Arc;

use datafusion::prelude::SessionContext;

use crate::allele::{match_allele_relaxed_udf, match_allele_udf, vep_allele_udf};
use crate::table_function::LookupFunction;

/// Test-only convenience: create a session with ranges + VEP functions.
#[cfg(test)]
pub(crate) fn create_vep_session() -> SessionContext {
    let ctx = datafusion_bio_function_ranges::create_bio_session();
    register_vep_functions(&ctx);
    ctx
}

/// Register all VEP functions on the given session context.
///
/// Registers:
/// - `match_allele(ref, alt, allele_string)` — scalar UDF
/// - `match_allele_relaxed(ref, alt, allele_string)` — scalar UDF
/// - `vep_allele(ref, alt)` — scalar UDF
/// - `lookup_variants(vcf_table, cache_table [, columns [, prune]])` — table function
pub fn register_vep_functions(ctx: &SessionContext) {
    ctx.register_udf(match_allele_udf());
    ctx.register_udf(match_allele_relaxed_udf());
    ctx.register_udf(vep_allele_udf());

    let session = Arc::new(ctx.clone());
    ctx.register_udtf("lookup_variants", Arc::new(LookupFunction::new(session)));
}
