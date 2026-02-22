//! VEP (Variant Effect Predictor) annotation functions for Apache DataFusion.
//!
//! Provides:
//! - `lookup_variants()` table function for known variant lookup via interval join
//! - `match_allele()` scalar UDF for allele matching
//! - `vep_allele()` scalar UDF for VCF→VEP allele conversion

pub mod allele;
pub mod coordinate;
pub mod lookup_provider;
pub mod schema_contract;
pub mod table_function;

use std::sync::Arc;

use datafusion::prelude::SessionContext;
use datafusion_bio_function_ranges::create_bio_session;

use crate::allele::{match_allele_udf, vep_allele_udf};
use crate::table_function::LookupFunction;

/// Register all VEP functions on the given session context.
///
/// Registers:
/// - `match_allele(ref, alt, allele_string)` — scalar UDF
/// - `vep_allele(ref, alt)` — scalar UDF
/// - `lookup_variants(vcf_table, cache_table [, columns [, prune]])` — table function
pub fn register_vep_functions(ctx: &SessionContext) {
    ctx.register_udf(match_allele_udf());
    ctx.register_udf(vep_allele_udf());

    let session = Arc::new(ctx.clone());
    ctx.register_udtf("lookup_variants", Arc::new(LookupFunction::new(session)));
}

/// Create a session context with both ranges and VEP functions registered.
///
/// This builds on `create_bio_session()` from `bio-function-ranges` which sets up
/// the `BioQueryPlanner` with `IntervalJoinExec` optimization.
pub fn create_vep_session() -> SessionContext {
    let ctx = create_bio_session();
    register_vep_functions(&ctx);
    ctx
}
