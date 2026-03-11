//! VEP (Variant Effect Predictor) annotation functions for Apache DataFusion.
//!
//! Provides:
//! - `lookup_variants()` table function for known variant lookup via equi-join
//! - `match_allele()` scalar UDF for allele matching
//! - `vep_allele()` scalar UDF for VCF→VEP allele conversion

pub mod allele;
pub mod annotate_provider;
pub mod annotate_table_function;
pub mod annotation_store;
pub mod config;
pub mod coordinate;
pub mod golden_benchmark;
pub mod hgvs;
#[cfg(feature = "kv-cache")]
pub mod kv_cache;
pub mod lookup_provider;
pub mod schema_contract;
pub mod so_terms;
pub mod table_function;
pub mod transcript_consequence;

pub use config::AnnotationConfig;

use std::sync::Arc;

use datafusion::prelude::SessionContext;

use crate::allele::{
    match_allele_relaxed_udf, match_allele_udf, vep_allele_udf, vep_norm_end_udf,
    vep_norm_start_udf,
};
use crate::annotate_table_function::AnnotateFunction;
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
/// - `lookup_variants(vcf_table, cache_table [, columns [, match_mode [, extended_probes]]])` — table function
/// - `annotate_vep(vcf_table, cache_source, backend [, options_json])` — table function
pub fn register_vep_functions(ctx: &SessionContext) {
    ctx.register_udf(match_allele_udf());
    ctx.register_udf(match_allele_relaxed_udf());
    ctx.register_udf(vep_allele_udf());
    ctx.register_udf(vep_norm_start_udf());
    ctx.register_udf(vep_norm_end_udf());

    let session = Arc::new(ctx.clone());
    ctx.register_udtf("lookup_variants", Arc::new(LookupFunction::new(session)));
    let session = Arc::new(ctx.clone());
    ctx.register_udtf("annotate_vep", Arc::new(AnnotateFunction::new(session)));
}
