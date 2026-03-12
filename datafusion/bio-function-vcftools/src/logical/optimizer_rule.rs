//! Logical optimizer rule for fused array transform.
//!
//! This rule detects the pattern `Aggregate → Project* → Unnest` and rewrites
//! it to a `FusedArrayTransform` node for memory-efficient processing.
//!
//! # Supported Patterns
//!
//! Simple pattern (identity):
//! ```text
//! Aggregate → Projection → Unnest
//! ```
//!
//! Transformation CTE pattern:
//! ```text
//! Aggregate → SubqueryAlias → Projection → SubqueryAlias → Projection → Unnest
//! ```

use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use datafusion::common::Result;
use datafusion::common::tree_node::{Transformed, TreeNode};
use datafusion::logical_expr::expr::{AggregateFunction, WindowFunction};
use datafusion::logical_expr::{Expr, Extension, LogicalPlan, SubqueryAlias};
use datafusion::optimizer::{OptimizerConfig, OptimizerRule};
use log::{debug, trace};

use super::fused_array_transform::FusedArrayTransform;

/// Global thread-safe flag to enable/disable the fused array transform optimization.
///
/// This is the preferred way to enable the optimization in production code.
/// Use [`enable_fused_array_transform`] and [`disable_fused_array_transform`]
/// to control this flag, or use [`set_fused_array_transform_enabled`] for
/// conditional enabling.
///
/// The environment variable `BIO_FUSED_ARRAY_TRANSFORM` is still supported
/// for backward compatibility, but the atomic flag takes precedence.
static FUSED_ARRAY_TRANSFORM_ENABLED: AtomicBool = AtomicBool::new(false);

/// Enable the fused array transform optimization globally.
///
/// This is thread-safe and can be called from any thread.
///
/// # Example
///
/// ```
/// use datafusion_bio_function_vcftools::logical::optimizer_rule::enable_fused_array_transform;
///
/// enable_fused_array_transform();
/// ```
pub fn enable_fused_array_transform() {
    FUSED_ARRAY_TRANSFORM_ENABLED.store(true, Ordering::SeqCst);
}

/// Disable the fused array transform optimization globally.
///
/// This is thread-safe and can be called from any thread.
///
/// # Example
///
/// ```
/// use datafusion_bio_function_vcftools::logical::optimizer_rule::disable_fused_array_transform;
///
/// disable_fused_array_transform();
/// ```
pub fn disable_fused_array_transform() {
    FUSED_ARRAY_TRANSFORM_ENABLED.store(false, Ordering::SeqCst);
}

/// Set the fused array transform optimization state.
///
/// This is thread-safe and can be called from any thread.
///
/// # Arguments
///
/// * `enabled` - `true` to enable the optimization, `false` to disable it.
///
/// # Example
///
/// ```
/// use datafusion_bio_function_vcftools::logical::optimizer_rule::set_fused_array_transform_enabled;
///
/// // Enable conditionally based on a configuration
/// set_fused_array_transform_enabled(true);
/// ```
pub fn set_fused_array_transform_enabled(enabled: bool) {
    FUSED_ARRAY_TRANSFORM_ENABLED.store(enabled, Ordering::SeqCst);
}

/// Check if the fused array transform optimization is currently enabled.
///
/// This checks the global atomic flag first, then falls back to the
/// environment variable `BIO_FUSED_ARRAY_TRANSFORM` for backward compatibility.
///
/// # Returns
///
/// `true` if the optimization is enabled, `false` otherwise.
pub fn is_fused_array_transform_enabled() -> bool {
    // Check atomic flag first (takes precedence)
    if FUSED_ARRAY_TRANSFORM_ENABLED.load(Ordering::SeqCst) {
        return true;
    }
    // Fall back to environment variable for backward compatibility
    std::env::var("BIO_FUSED_ARRAY_TRANSFORM")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false)
}

/// Get a human-readable name for a LogicalPlan node type.
fn plan_type_name(plan: &LogicalPlan) -> &'static str {
    match plan {
        LogicalPlan::Projection(_) => "Projection",
        LogicalPlan::Filter(_) => "Filter",
        LogicalPlan::Window(_) => "Window",
        LogicalPlan::Aggregate(_) => "Aggregate",
        LogicalPlan::Sort(_) => "Sort",
        LogicalPlan::Join(_) => "Join",
        LogicalPlan::SubqueryAlias(_) => "SubqueryAlias",
        LogicalPlan::TableScan(_) => "TableScan",
        LogicalPlan::Unnest(_) => "Unnest",
        LogicalPlan::Extension(_) => "Extension",
        _ => "Other",
    }
}

/// Logical optimizer rule that detects `Aggregate → Project → Unnest` patterns
/// and replaces them with a `FusedArrayTransform` logical node.
///
/// # Pattern Detection
///
/// The rule looks for:
/// ```text
/// Aggregate (GROUP BY row_idx with all array_agg functions)
///   └─ Projection (optional transformations)
///        └─ Unnest (unnesting array columns)
///             └─ ...
/// ```
///
/// # Requirements for Optimization
///
/// 1. All aggregate functions must be `array_agg` (or `list_agg`)
/// 2. The grouping columns should include a row identifier
/// 3. The unnest must be on array columns
#[derive(Debug, Default)]
pub struct FusedArrayTransformOptimizerRule;

impl FusedArrayTransformOptimizerRule {
    /// Create a new instance of the optimizer rule.
    pub fn new() -> Self {
        Self
    }

    /// Check if the optimization is enabled.
    ///
    /// This delegates to the global [`is_fused_array_transform_enabled`] function,
    /// which checks the atomic flag first, then falls back to the environment
    /// variable `BIO_FUSED_ARRAY_TRANSFORM` for backward compatibility.
    pub fn is_enabled(&self) -> bool {
        is_fused_array_transform_enabled()
    }
}

impl OptimizerRule for FusedArrayTransformOptimizerRule {
    fn name(&self) -> &str {
        "fused_array_transform"
    }

    fn rewrite(
        &self,
        plan: LogicalPlan,
        _config: &dyn OptimizerConfig,
    ) -> Result<Transformed<LogicalPlan>> {
        if !self.is_enabled() {
            return Ok(Transformed::no(plan));
        }

        debug!("FusedArrayTransformOptimizerRule: rewrite() called");

        plan.transform_up(|node| match try_optimize(&node) {
            Some(Ok(optimized)) => {
                debug!("FusedArrayTransformOptimizerRule: applied optimization");
                Ok(Transformed::yes(optimized))
            }
            Some(Err(e)) => {
                debug!("FusedArrayTransformOptimizerRule: optimization failed: {e}");
                Ok(Transformed::no(node))
            }
            None => Ok(Transformed::no(node)),
        })
    }
}

/// Skip transparent wrapper nodes like SubqueryAlias to find the actual plan node.
fn skip_wrappers(plan: &LogicalPlan) -> &LogicalPlan {
    match plan {
        LogicalPlan::SubqueryAlias(alias) => {
            trace!("Skipping SubqueryAlias: {}", alias.alias);
            skip_wrappers(alias.input.as_ref())
        }
        other => other,
    }
}

/// Result of traversing from Aggregate to Unnest.
struct TraversalResult<'a> {
    /// The Unnest node we found
    unnest: &'a datafusion::logical_expr::Unnest,
    /// Column mapping: output_name -> expression defining it
    /// This is built by traversing projections and composing expressions
    column_definitions: HashMap<String, Expr>,
    /// Column names that uniquely identify each input row (e.g., from ROW_NUMBER() OVER ())
    /// The optimization can only be applied if the GROUP BY includes at least one of these.
    row_identity_columns: HashSet<String>,
}

/// Traverse from the Aggregate's input down to Unnest, collecting column definitions along the way.
///
/// This handles transformation CTEs by traversing multiple Projection layers and
/// composing expressions.
///
/// # Returns
///
/// - `Ok(Some(result))` if the pattern matches and traversal succeeded
/// - `Ok(None)` if the pattern doesn't match (no Unnest found)
/// - `Err(e)` if traversal failed due to expression resolution errors
fn traverse_to_unnest(plan: &LogicalPlan) -> Result<Option<TraversalResult<'_>>> {
    // Skip SubqueryAlias wrappers
    let plan = skip_wrappers(plan);

    match plan {
        LogicalPlan::Unnest(unnest) => {
            // Base case: we found Unnest
            // Map OUTPUT column names (with depth) to INPUT column names (without depth)
            // because FusedArrayTransform operates on unnest's input, not output
            let mut column_definitions = HashMap::new();

            // Build mapping from exec_columns (input placeholder names) to output names
            // exec_columns contains: __unnest_placeholder(indexed.values_a)
            // output schema contains: __unnest_placeholder(indexed.values_a,depth=1)
            for exec_col in &unnest.exec_columns {
                let input_name = exec_col.name().to_string();
                // The output column name has ",depth=N)" appended before the closing paren
                // Input: __unnest_placeholder(col) -> Output: __unnest_placeholder(col,depth=1)
                // Use exact prefix match with ",depth=" to avoid matching columns that share a prefix
                // (e.g., values_a and values_ab should not match each other)
                let input_prefix = input_name.trim_end_matches(')');
                let expected_output_prefix = format!("{input_prefix},depth=");

                // Find the corresponding output column
                let unnest_plan = LogicalPlan::Unnest(unnest.clone());
                let output_schema = unnest_plan.schema();

                for (_qualifier, field) in output_schema.iter() {
                    let output_name = field.name();
                    // Check if this output column corresponds to this input column exactly
                    // by matching the expected pattern: input_prefix + ",depth=" + N + ")"
                    if output_name.starts_with(&expected_output_prefix) {
                        // Map the output name to an expression referencing the INPUT name
                        column_definitions.insert(
                            output_name.to_string(),
                            Expr::Column(input_name.clone().into()),
                        );
                    }
                }

                // Also add the input column directly
                column_definitions.insert(input_name.clone(), Expr::Column(input_name.into()));
            }

            // Also include passthrough columns from the input schema
            let input_schema = unnest.input.schema();
            for (qualifier, field) in input_schema.iter() {
                let name = field.name().to_string();
                column_definitions.entry(name.clone()).or_insert_with(|| {
                    let qualified_col = datafusion::common::Column::new(qualifier.cloned(), name);
                    Expr::Column(qualified_col)
                });
            }

            // Scan the unnest input for row-identity columns (e.g., from ROW_NUMBER() OVER ())
            let row_identity_columns = find_row_identity_columns(unnest.input.as_ref());

            Ok(Some(TraversalResult {
                unnest,
                column_definitions,
                row_identity_columns,
            }))
        }
        LogicalPlan::Projection(projection) => {
            // Recurse into child
            let Some(child_result) = traverse_to_unnest(projection.input.as_ref())? else {
                return Ok(None);
            };

            // Build new column definitions by resolving projection expressions
            let mut new_definitions = HashMap::new();
            // Track row-identity columns through this projection
            let mut new_row_identity_columns = HashSet::new();

            for expr in &projection.expr {
                let (alias, inner_expr) = match expr {
                    Expr::Alias(alias_expr) => {
                        (alias_expr.name.clone(), alias_expr.expr.as_ref().clone())
                    }
                    // For column references, use just the name (not qualified)
                    Expr::Column(col) => (col.name().to_string(), expr.clone()),
                    other => (expr.schema_name().to_string(), other.clone()),
                };

                // Check if this column references a row-identity column
                if let Expr::Column(col) = &inner_expr {
                    if child_result.row_identity_columns.contains(col.name()) {
                        trace!("Row-identity column aliased: {} -> {}", col.name(), alias);
                        new_row_identity_columns.insert(alias.clone());
                    }
                }

                // Resolve the expression against child definitions
                let resolved = resolve_expr(&inner_expr, &child_result.column_definitions)?;
                new_definitions.insert(alias, resolved);
            }

            Ok(Some(TraversalResult {
                unnest: child_result.unnest,
                column_definitions: new_definitions,
                row_identity_columns: new_row_identity_columns,
            }))
        }
        _ => {
            trace!(
                "traverse_to_unnest: expected Projection or Unnest, got {}",
                plan_type_name(plan)
            );
            Ok(None)
        }
    }
}

/// Resolve an expression by substituting column references with their definitions.
/// Uses DataFusion's `transform` to recursively traverse all expression variants.
///
/// # Errors
///
/// Returns an error if the expression tree traversal fails (e.g., due to an
/// unexpected expression variant or internal DataFusion error).
fn resolve_expr(expr: &Expr, definitions: &HashMap<String, Expr>) -> Result<Expr> {
    expr.clone()
        .transform(|e| {
            if let Expr::Column(col) = &e {
                if let Some(definition) = definitions.get(col.name()) {
                    return Ok(Transformed::yes(definition.clone()));
                }
            }
            Ok(Transformed::no(e))
        })
        .map(|t| t.data)
}

/// Find columns that uniquely identify each row by scanning for ROW_NUMBER() OVER () patterns.
///
/// This traverses the plan tree looking for Window nodes that contain row_number()
/// with an empty PARTITION BY clause (which produces unique values per row).
/// It tracks these column names through projections/aliases.
///
/// # Arguments
///
/// * `plan` - The logical plan to scan (typically the input to Unnest)
///
/// # Returns
///
/// A set of column names that are guaranteed to uniquely identify each row.
fn find_row_identity_columns(plan: &LogicalPlan) -> HashSet<String> {
    let mut result = HashSet::new();

    // Recursively scan the plan tree
    fn scan_plan(plan: &LogicalPlan, identity_cols: &mut HashSet<String>) {
        let plan = skip_wrappers(plan);

        match plan {
            LogicalPlan::Window(window) => {
                // Check each window function expression
                for expr in &window.window_expr {
                    if let Some(col_name) = extract_row_identity_column(expr) {
                        identity_cols.insert(col_name);
                    }
                }
                // Continue scanning the input
                scan_plan(window.input.as_ref(), identity_cols);
            }
            LogicalPlan::Projection(projection) => {
                // Scan child first to collect identity columns
                scan_plan(projection.input.as_ref(), identity_cols);

                // Track aliases: if a projection aliases an identity column,
                // the alias is also an identity column
                let mut new_aliases = Vec::new();
                for expr in &projection.expr {
                    if let Expr::Alias(alias) = expr {
                        if let Expr::Column(col) = alias.expr.as_ref() {
                            if identity_cols.contains(col.name()) {
                                new_aliases.push(alias.name.clone());
                            }
                        }
                    }
                }
                for alias in new_aliases {
                    identity_cols.insert(alias);
                }
            }
            other => {
                // For other node types, scan children
                for child in other.inputs() {
                    scan_plan(child, identity_cols);
                }
            }
        }
    }

    scan_plan(plan, &mut result);
    result
}

/// Check if an expression is a ROW_NUMBER() window function with empty PARTITION BY,
/// and extract its output column name.
///
/// ROW_NUMBER() OVER () produces unique sequential values for each row, making it
/// a valid row identity column. However, ROW_NUMBER() OVER (PARTITION BY x) does not
/// produce globally unique values, so we only accept empty partition_by.
fn extract_row_identity_column(expr: &Expr) -> Option<String> {
    // Handle aliased expressions (e.g., ROW_NUMBER() OVER () as row_idx)
    let (output_name, inner_expr) = match expr {
        Expr::Alias(alias) => (alias.name.clone(), alias.expr.as_ref()),
        other => (other.schema_name().to_string(), other),
    };

    // Check if it's a window function (WindowFunction is boxed in Expr::WindowFunction)
    if let Expr::WindowFunction(window_func) = inner_expr {
        let WindowFunction { fun, params } = window_func.as_ref();
        // Check if it's row_number (case-insensitive)
        let func_name = fun.name().to_lowercase();
        if func_name == "row_number" {
            // Only accept if partition_by is empty (ensures global uniqueness)
            if params.partition_by.is_empty() {
                trace!("Found row-identity column: {output_name} (from row_number())");
                return Some(output_name);
            } else {
                trace!(
                    "row_number() with PARTITION BY is not a row-identity column: {output_name}"
                );
            }
        }
    }
    None
}

/// Attempt to detect and optimize the pattern.
///
/// Returns:
/// - `Some(Ok(plan))` if pattern matched and optimization succeeded
/// - `Some(Err(e))` if pattern matched but optimization failed
/// - `None` if pattern doesn't match
fn try_optimize(plan: &LogicalPlan) -> Option<Result<LogicalPlan>> {
    // Check if this is an Aggregate node
    let LogicalPlan::Aggregate(aggregate) = plan else {
        return None;
    };

    trace!("Found Aggregate node, checking for optimization pattern");

    // Check if all aggregate expressions are array_agg
    if !all_array_agg(&aggregate.aggr_expr) {
        trace!("Not all aggregates are array_agg");
        return None;
    }

    trace!("All aggregates are array_agg, traversing to find Unnest");

    // Get the SubqueryAlias name if the immediate child is a SubqueryAlias
    // (we need to wrap our result with the same alias to preserve column qualifiers)
    let subquery_alias_name = if let LogicalPlan::SubqueryAlias(alias) = aggregate.input.as_ref() {
        Some(alias.alias.clone())
    } else {
        None
    };

    // Traverse to find Unnest while collecting column definitions
    let traversal = match traverse_to_unnest(aggregate.input.as_ref()) {
        Ok(Some(t)) => t,
        Ok(None) => {
            trace!("traverse_to_unnest returned None");
            return None;
        }
        Err(e) => {
            // Pattern matched but traversal failed - propagate error
            return Some(Err(e));
        }
    };
    let unnest_plan = traversal.unnest;
    let column_definitions = traversal.column_definitions;
    let row_identity_columns = traversal.row_identity_columns;

    trace!(
        "Found Unnest, column_definitions: {:?}, row_identity_columns: {:?}",
        column_definitions.keys().collect::<Vec<_>>(),
        row_identity_columns
    );

    // Get the input to unnest (this is what we'll use as our input)
    let unnest_input = unnest_plan.input.clone();

    // Extract unnest columns from exec_columns (the INPUT column names, without depth)
    // These are the columns we'll read arrays from in FusedArrayTransform
    let array_columns: Vec<String> = unnest_plan
        .exec_columns
        .iter()
        .map(|col| col.name().to_string())
        .collect();

    if array_columns.is_empty() {
        debug!("No array columns found in Unnest");
        return None;
    }

    // Extract passthrough columns from group-by (all grouping columns that are not aggregated values)
    // These are the columns from the input that we preserve unchanged
    let passthrough_columns: Vec<String> = aggregate
        .group_expr
        .iter()
        .filter_map(|expr| {
            if let Expr::Column(col) = expr {
                Some(col.name().to_string())
            } else {
                None
            }
        })
        .collect();

    // CRITICAL: Verify that GROUP BY includes a row-identity column.
    //
    // FusedArrayTransformExec emits at most one output row per input row and never
    // performs cross-row aggregation. If the GROUP BY doesn't uniquely identify each
    // original input row, the optimization would produce incorrect results.
    //
    // For example, GROUP BY metadata (where multiple rows share the same metadata)
    // should merge values across rows with array_agg, but the optimized plan would
    // return separate arrays per row.
    let has_row_identity = passthrough_columns
        .iter()
        .any(|col| row_identity_columns.contains(col));

    if !has_row_identity {
        debug!(
            "GROUP BY columns {passthrough_columns:?} do not include any row-identity column \
             (available: {row_identity_columns:?}). Cannot safely apply FusedArrayTransform optimization."
        );
        return None;
    }

    // Extract output column names and transform expressions from aggregate expressions
    // For each array_agg(col), we resolve 'col' to get the actual transformation
    let mut output_columns = Vec::new();
    let mut transform_exprs = Vec::new();

    for expr in &aggregate.aggr_expr {
        // Get output column name
        output_columns.push(expr.schema_name().to_string());

        // Extract the argument to array_agg and resolve it
        if let Expr::AggregateFunction(AggregateFunction { params, .. }) = expr {
            if let Some(arg) = params.args.first() {
                // Resolve the argument through the column definitions
                let resolved = match resolve_expr(arg, &column_definitions) {
                    Ok(r) => r,
                    Err(e) => return Some(Err(e)),
                };
                trace!("Resolved array_agg argument: {arg:?} -> {resolved:?}");
                transform_exprs.push(resolved);
            } else {
                // No argument, use identity
                transform_exprs.push(Expr::Literal(datafusion::common::ScalarValue::Null, None));
            }
        }
    }

    debug!(
        "Pattern matched: arrays={array_columns:?}, passthrough={passthrough_columns:?}, outputs={output_columns:?}, transforms={:?}",
        transform_exprs
            .iter()
            .map(|e| e.to_string())
            .collect::<Vec<_>>()
    );

    // Create the FusedArrayTransform node
    let fused_node_result = FusedArrayTransform::try_new(
        unnest_input,
        array_columns,
        passthrough_columns,
        output_columns,
        transform_exprs,
    );

    Some(fused_node_result.and_then(|fused_node| {
        let extension_plan = LogicalPlan::Extension(Extension {
            node: Arc::new(fused_node),
        });

        // Wrap in SubqueryAlias to preserve column qualifiers if original had one
        if let Some(alias_name) = subquery_alias_name {
            trace!("Wrapping FusedArrayTransform in SubqueryAlias: {alias_name}");
            SubqueryAlias::try_new(Arc::new(extension_plan), alias_name)
                .map(LogicalPlan::SubqueryAlias)
        } else {
            Ok(extension_plan)
        }
    }))
}

/// Check if all aggregate expressions are array_agg or list_agg without unsupported modifiers.
///
/// We reject ORDER BY, DISTINCT and FILTER as these require semantics that the fused transform doesn't currently implement.
fn all_array_agg(exprs: &[Expr]) -> bool {
    exprs.iter().all(|expr| {
        if let Expr::AggregateFunction(AggregateFunction { func, params }) = expr {
            let name = func.name().to_lowercase();
            let is_array_agg = name.contains("array_agg") || name.contains("list_agg");

            // Reject unsupported modifiers
            let has_unsupported_modifiers =
                params.distinct || params.filter.is_some() || !params.order_by.is_empty();

            is_array_agg && !has_unsupported_modifiers
        } else {
            false
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::functions_aggregate::array_agg::array_agg;
    use datafusion::logical_expr::expr::{AggregateFunctionParams, Sort};
    use datafusion::logical_expr::{AggregateUDF, col};

    /// Helper to create a simple array_agg(col) expression.
    fn simple_array_agg(column_name: &str) -> Expr {
        array_agg(col(column_name))
    }

    /// Helper to create an array_agg with ORDER BY modifier.
    fn array_agg_with_order_by(agg_col: &str, order_col: &str) -> Expr {
        // Get the array_agg UDF
        let func: Arc<AggregateUDF> = datafusion::functions_aggregate::array_agg::array_agg_udaf();

        Expr::AggregateFunction(AggregateFunction {
            func,
            params: AggregateFunctionParams {
                args: vec![col(agg_col)],
                distinct: false,
                filter: None,
                order_by: vec![Sort::new(col(order_col), true, false)],
                null_treatment: None,
            },
        })
    }

    /// Helper to create an array_agg with DISTINCT modifier.
    fn array_agg_with_distinct(column_name: &str) -> Expr {
        let func: Arc<AggregateUDF> = datafusion::functions_aggregate::array_agg::array_agg_udaf();

        Expr::AggregateFunction(AggregateFunction {
            func,
            params: AggregateFunctionParams {
                args: vec![col(column_name)],
                distinct: true,
                filter: None,
                order_by: vec![],
                null_treatment: None,
            },
        })
    }

    /// Helper to create an array_agg with FILTER modifier.
    fn array_agg_with_filter(agg_col: &str, filter_expr: Expr) -> Expr {
        let func: Arc<AggregateUDF> = datafusion::functions_aggregate::array_agg::array_agg_udaf();

        Expr::AggregateFunction(AggregateFunction {
            func,
            params: AggregateFunctionParams {
                args: vec![col(agg_col)],
                distinct: false,
                filter: Some(Box::new(filter_expr)),
                order_by: vec![],
                null_treatment: None,
            },
        })
    }

    #[test]
    fn test_rule_name() {
        let rule = FusedArrayTransformOptimizerRule::new();
        assert_eq!(rule.name(), "fused_array_transform");
    }

    #[test]
    fn test_rule_disabled_by_default() {
        // Ensure the optimization is disabled
        disable_fused_array_transform();

        let rule = FusedArrayTransformOptimizerRule::new();
        assert!(!rule.is_enabled());
    }

    #[test]
    fn test_rule_enabled_via_api() {
        // Enable the optimization using the thread-safe API
        enable_fused_array_transform();

        let rule = FusedArrayTransformOptimizerRule::new();
        assert!(rule.is_enabled());

        // Clean up
        disable_fused_array_transform();
    }

    #[test]
    fn test_all_array_agg_simple() {
        // Simple array_agg without modifiers should be accepted
        let exprs = vec![simple_array_agg("col_a"), simple_array_agg("col_b")];
        assert!(all_array_agg(&exprs), "Simple array_agg should be accepted");
    }

    #[test]
    fn test_all_array_agg_rejects_order_by() {
        // array_agg with ORDER BY should be REJECTED because the fused transform
        // doesn't preserve element ordering
        let exprs = vec![array_agg_with_order_by("val", "sort_key")];
        assert!(
            !all_array_agg(&exprs),
            "array_agg with ORDER BY should be rejected"
        );
    }

    #[test]
    fn test_all_array_agg_rejects_distinct() {
        // array_agg with DISTINCT should be REJECTED because the fused transform
        // doesn't remove duplicates
        let exprs = vec![array_agg_with_distinct("val")];
        assert!(
            !all_array_agg(&exprs),
            "array_agg with DISTINCT should be rejected"
        );
    }

    #[test]
    fn test_all_array_agg_rejects_filter() {
        // array_agg with FILTER should be REJECTED because the fused transform
        // doesn't apply filtering
        let filter = col("flag").eq(datafusion::logical_expr::lit(true));
        let exprs = vec![array_agg_with_filter("val", filter)];
        assert!(
            !all_array_agg(&exprs),
            "array_agg with FILTER should be rejected"
        );
    }

    #[test]
    fn test_all_array_agg_mixed_simple_and_with_modifiers() {
        // If any aggregate has modifiers, the whole set should be rejected
        let exprs = vec![
            simple_array_agg("col_a"),
            array_agg_with_order_by("col_b", "sort_key"),
        ];
        assert!(
            !all_array_agg(&exprs),
            "Mixed simple and ORDER BY should be rejected"
        );
    }
}
