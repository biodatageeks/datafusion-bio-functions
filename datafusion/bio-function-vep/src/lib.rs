//! VEP (Variant Effect Predictor) annotation functions for Apache DataFusion.
//!
//! Provides:
//! - `lookup_variants()` table function for known variant lookup via equi-join
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
    clippy::unwrap_or_default
)]

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
pub mod miss_worklist;
pub mod partitioned_cache;
pub mod schema_contract;
pub mod so_terms;
pub mod table_function;
pub mod transcript_consequence;
pub mod variant_lookup_exec;
pub mod vcf_sink;

pub use config::AnnotationConfig;

use std::sync::Arc;

use datafusion::prelude::SessionContext;

use crate::allele::{
    match_allele_relaxed_udf, match_allele_udf, vep_allele_udf, vep_norm_end_udf,
    vep_norm_start_udf,
};
use crate::annotate_table_function::AnnotateFunction;
use crate::table_function::LookupFunction;

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
