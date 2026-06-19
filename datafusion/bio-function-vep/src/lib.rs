//! VEP (Variant Effect Predictor) annotation functions for Apache DataFusion.
//!
//! Provides:
//! - `annotate_vep()` table function for consequence annotation over a Lance cache
//! - `match_allele()` scalar UDF for allele matching
//! - `vep_allele()` scalar UDF for VCF→VEP allele conversion
#![allow(
    dead_code,
    unused_imports,
    unused_variables,
    clippy::too_many_arguments,
    clippy::collapsible_if,
    clippy::collapsible_else_if,
    clippy::doc_lazy_continuation,
    clippy::doc_overindented_list_items,
    clippy::redundant_guards,
    clippy::manual_div_ceil,
    clippy::field_reassign_with_default,
    clippy::useless_vec,
    clippy::manual_flatten,
    clippy::unnecessary_lazy_evaluations,
    clippy::should_implement_trait,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::or_fun_call,
    clippy::vec_init_then_push,
    clippy::clone_on_copy,
    clippy::single_element_loop,
    clippy::used_underscore_items,
    clippy::empty_line_after_doc_comments,
    clippy::manual_contains,
    clippy::collapsible_str_replace,
    clippy::unnecessary_map_or,
    clippy::assigning_clones,
    clippy::map_entry,
    clippy::cloned_ref_to_slice_refs,
    clippy::unwrap_or_default,
    clippy::needless_option_as_deref
)]

pub mod allele;
pub mod annotate_provider;
pub mod annotate_table_function;
pub mod annotation_store;
#[cfg(feature = "cache-builder")]
pub mod cache_builder;
pub(crate) mod cache_common;
pub(crate) mod cache_source;
pub(crate) mod colocated;
pub mod coordinate;
pub mod golden_benchmark;
pub mod hgvs;
#[cfg(feature = "lance-cache")]
pub mod lance_cache;
pub mod lookup_provider;
pub mod miss_worklist;
pub(crate) mod ordered_drain;
pub mod partitioned_cache;
pub(crate) mod pipeline_trace;
pub mod schema_contract;
pub mod so_terms;
pub mod transcript_consequence;
pub mod vcf_sink;
pub(crate) mod window_planner;

use std::sync::Arc;

use datafusion::prelude::SessionContext;

use crate::allele::{
    match_allele_relaxed_udf, match_allele_udf, vep_allele_udf, vep_norm_end_udf,
    vep_norm_start_udf,
};
use crate::annotate_table_function::AnnotateFunction;

/// Test-only convenience: create a session with VEP functions.
#[cfg(test)]
pub(crate) fn create_vep_session() -> SessionContext {
    let ctx = SessionContext::new();
    register_vep_functions(&ctx);
    ctx
}

/// Register all VEP functions on the given session context.
///
/// Registers:
/// - `match_allele(ref, alt, allele_string)` — scalar UDF
/// - `match_allele_relaxed(ref, alt, allele_string)` — scalar UDF
/// - `vep_allele(ref, alt)` — scalar UDF
/// - `annotate_vep(vcf_table, cache_source, backend [, options_json])` — table function
pub fn register_vep_functions(ctx: &SessionContext) {
    ctx.register_udf(match_allele_udf());
    ctx.register_udf(match_allele_relaxed_udf());
    ctx.register_udf(vep_allele_udf());
    ctx.register_udf(vep_norm_start_udf());
    ctx.register_udf(vep_norm_end_udf());

    let session = Arc::new(ctx.clone());
    ctx.register_udtf("annotate_vep", Arc::new(AnnotateFunction::new(session)));
}
