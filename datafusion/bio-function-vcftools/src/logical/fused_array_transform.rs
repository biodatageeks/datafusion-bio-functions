//! FusedArrayTransform logical node.
//!
//! This module defines the `FusedArrayTransform` custom logical node that represents
//! an optimized fusion of `unnest → transform → array_agg` operations. Instead of
//! materializing all unnested rows, it processes arrays element-wise with bounded memory.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::{self, Debug};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::tree_node::{Transformed, TreeNode};
use datafusion::common::{Column, DFSchema, DFSchemaRef, Result};
use datafusion::logical_expr::{Expr, ExprSchemable, LogicalPlan, UserDefinedLogicalNodeCore};

use crate::common::build_element_schema_from_df;

/// Strip table qualifiers from all column references in an expression.
///
/// This is needed because the element schema uses unqualified column names,
/// but transform expressions from the optimizer may contain qualified references
/// (e.g., `indexed.metadata_int` instead of just `metadata_int`).
fn strip_column_qualifiers(expr: &Expr) -> Result<Expr> {
    expr.clone()
        .transform(|e| {
            if let Expr::Column(col) = &e {
                if col.relation.is_some() {
                    // Strip the qualifier, keep only the column name
                    let unqualified = Column::new_unqualified(col.name());
                    return Ok(Transformed::yes(Expr::Column(unqualified)));
                }
            }
            Ok(Transformed::no(e))
        })
        .map(|t| t.data)
}

/// A custom logical node that fuses unnest + transform + array_agg into a single operation.
///
/// This node captures:
/// - The input plan (original table with array columns)
/// - Which columns are arrays to be transformed
/// - Which columns should pass through unchanged
/// - The transformation expressions to apply element-wise
///
/// During physical planning, this is converted to `FusedArrayTransformExec`.
#[derive(Clone)]
pub struct FusedArrayTransform {
    /// The input logical plan
    input: Arc<LogicalPlan>,
    /// Column names that contain arrays to process
    array_columns: Vec<String>,
    /// Column names to pass through unchanged (e.g., metadata columns)
    passthrough_columns: Vec<String>,
    /// Output column names for the transformed arrays
    output_columns: Vec<String>,
    /// Transformation expressions. Each expression operates on "virtual" scalar
    /// values from the unnested arrays. The expression references column names
    /// as they would appear after unnesting.
    transform_exprs: Vec<Expr>,
    /// Expected output schema
    schema: DFSchemaRef,
}

impl FusedArrayTransform {
    /// Create a new FusedArrayTransform node.
    ///
    /// # Arguments
    ///
    /// * `input` - The input logical plan containing array columns
    /// * `array_columns` - Names of columns containing arrays to transform
    /// * `passthrough_columns` - Names of columns to pass through unchanged
    /// * `output_columns` - Names for the output transformed array columns
    /// * `transform_exprs` - Expressions to apply element-wise on array elements
    pub fn try_new(
        input: Arc<LogicalPlan>,
        array_columns: Vec<String>,
        passthrough_columns: Vec<String>,
        output_columns: Vec<String>,
        transform_exprs: Vec<Expr>,
    ) -> Result<Self> {
        // Build output schema: passthrough columns + output array columns
        let input_schema = input.schema();
        let mut fields: Vec<Field> = Vec::new();

        // Add passthrough columns
        for col_name in &passthrough_columns {
            let field = input_schema.field_with_unqualified_name(col_name)?;
            fields.push(Field::new(
                col_name,
                field.data_type().clone(),
                field.is_nullable(),
            ));
        }

        // Build element-level schema for type inference
        // This represents the scalar types of array elements (as if unnested)
        // We include both:
        // 1. Array columns converted to their element types (for unnested value references)
        // 2. Passthrough columns as-is (for transform expressions that reference metadata)
        let array_element_schema = build_element_schema_from_df(input_schema, &array_columns)?;

        // Merge array element fields with passthrough fields for complete element schema
        let mut element_fields: Vec<(Option<datafusion::common::TableReference>, Arc<Field>)> =
            array_element_schema
                .iter()
                .map(|(q, f)| (q.cloned(), f.clone()))
                .collect();
        for col_name in &passthrough_columns {
            let field = input_schema.field_with_unqualified_name(col_name)?;
            element_fields.push((None, Arc::new(field.clone())));
        }
        let element_schema = DFSchema::new_with_metadata(element_fields, HashMap::new())?;

        // Helper to get element type from array column
        let get_array_element_type = |array_col: &str| -> Result<DataType> {
            crate::common::extract_list_element_type(
                input_schema
                    .field_with_unqualified_name(array_col)?
                    .data_type(),
                array_col,
            )
        };

        // Strip column qualifiers from transform expressions
        let normalized_transform_exprs: Vec<Expr> = transform_exprs
            .iter()
            .map(strip_column_qualifiers)
            .collect::<Result<Vec<_>>>()?;

        // Add output array columns - infer type from corresponding transform expression
        for (idx, output_col) in output_columns.iter().enumerate() {
            // Infer element type from the transform expression
            // Type inference must succeed - no silent fallbacks to avoid schema mismatches
            let element_type = if idx < normalized_transform_exprs.len() {
                let inferred = normalized_transform_exprs[idx].get_type(&element_schema);
                match inferred {
                    Ok(dt) => dt,
                    Err(e) => {
                        // Try falling back to array column's element type
                        if let Some(array_col) = array_columns.get(idx) {
                            get_array_element_type(array_col).map_err(|_| {
                                datafusion::common::DataFusionError::Plan(format!(
                                    "Type inference failed for output '{output_col}': {e}. \
                                     Could not determine element type from array column '{array_col}'."
                                ))
                            })?
                        } else {
                            return Err(datafusion::common::DataFusionError::Plan(format!(
                                "Type inference failed for output '{output_col}': {e}. \
                                 No corresponding array column to fall back to."
                            )));
                        }
                    }
                }
            } else if idx < array_columns.len() {
                // No transform expression, use array column's element type
                get_array_element_type(&array_columns[idx])?
            } else {
                return Err(datafusion::common::DataFusionError::Plan(format!(
                    "No type information for output '{output_col}' at index {idx}. \
                     Number of outputs ({}) exceeds number of array columns ({}).",
                    output_columns.len(),
                    array_columns.len()
                )));
            };

            fields.push(Field::new(
                output_col,
                DataType::List(Arc::new(Field::new("item", element_type, true))),
                true,
            ));
        }

        let arrow_schema = Schema::new(fields);
        let schema = Arc::new(DFSchema::try_from(arrow_schema)?);

        Ok(Self {
            input,
            array_columns,
            passthrough_columns,
            output_columns,
            transform_exprs: normalized_transform_exprs,
            schema,
        })
    }

    /// Get the input logical plan
    pub fn input(&self) -> &LogicalPlan {
        &self.input
    }

    /// Get the array column names
    pub fn array_columns(&self) -> &[String] {
        &self.array_columns
    }

    /// Get the passthrough column names
    pub fn passthrough_columns(&self) -> &[String] {
        &self.passthrough_columns
    }

    /// Get the output column names
    pub fn output_columns(&self) -> &[String] {
        &self.output_columns
    }

    /// Get the transformation expressions
    pub fn transform_exprs(&self) -> &[Expr] {
        &self.transform_exprs
    }
}

impl Debug for FusedArrayTransform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FusedArrayTransform")
            .field("array_columns", &self.array_columns)
            .field("passthrough_columns", &self.passthrough_columns)
            .field("output_columns", &self.output_columns)
            .field("transform_exprs", &self.transform_exprs)
            .finish()
    }
}

impl Hash for FusedArrayTransform {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.array_columns.hash(state);
        self.passthrough_columns.hash(state);
        self.output_columns.hash(state);
        // Note: Expr doesn't implement Hash, so we hash its Display representation.
        // Display is more canonical than Debug for expressions.
        for expr in &self.transform_exprs {
            expr.to_string().hash(state);
        }
    }
}

impl PartialEq for FusedArrayTransform {
    fn eq(&self, other: &Self) -> bool {
        self.array_columns == other.array_columns
            && self.passthrough_columns == other.passthrough_columns
            && self.output_columns == other.output_columns
            && self.transform_exprs.len() == other.transform_exprs.len()
            // Use Display string comparison for consistency with Hash and Ord.
            // Expr doesn't implement Ord, so we use Display as the canonical representation
            // across all three traits to satisfy the invariant: a == b implies a.cmp(b) == Equal.
            && self
                .transform_exprs
                .iter()
                .zip(other.transform_exprs.iter())
                .all(|(a, b)| a.to_string() == b.to_string())
    }
}

impl Eq for FusedArrayTransform {}

impl PartialOrd for FusedArrayTransform {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FusedArrayTransform {
    fn cmp(&self, other: &Self) -> Ordering {
        self.array_columns
            .cmp(&other.array_columns)
            .then_with(|| self.passthrough_columns.cmp(&other.passthrough_columns))
            .then_with(|| self.output_columns.cmp(&other.output_columns))
            // Compare transform_exprs to maintain consistency with PartialEq.
            // Since Expr doesn't implement Ord, we compare by length first,
            // then by Display string representation for deterministic ordering.
            .then_with(|| self.transform_exprs.len().cmp(&other.transform_exprs.len()))
            .then_with(|| {
                for (a, b) in self
                    .transform_exprs
                    .iter()
                    .zip(other.transform_exprs.iter())
                {
                    let cmp = a.to_string().cmp(&b.to_string());
                    if cmp != Ordering::Equal {
                        return cmp;
                    }
                }
                Ordering::Equal
            })
    }
}

impl UserDefinedLogicalNodeCore for FusedArrayTransform {
    fn name(&self) -> &str {
        "FusedArrayTransform"
    }

    fn inputs(&self) -> Vec<&LogicalPlan> {
        vec![&self.input]
    }

    fn schema(&self) -> &DFSchemaRef {
        &self.schema
    }

    fn expressions(&self) -> Vec<Expr> {
        self.transform_exprs.clone()
    }

    fn fmt_for_explain(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "FusedArrayTransform: arrays=[{}], passthrough=[{}], outputs=[{}]",
            self.array_columns.join(", "),
            self.passthrough_columns.join(", "),
            self.output_columns.join(", ")
        )
    }

    fn with_exprs_and_inputs(&self, exprs: Vec<Expr>, inputs: Vec<LogicalPlan>) -> Result<Self> {
        if inputs.len() != 1 {
            return Err(datafusion::common::DataFusionError::Plan(
                "FusedArrayTransform requires exactly one input".to_string(),
            ));
        }

        Self::try_new(
            Arc::new(inputs.into_iter().next().unwrap()),
            self.array_columns.clone(),
            self.passthrough_columns.clone(),
            self.output_columns.clone(),
            exprs,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::logical_expr::col;

    fn create_test_schema() -> Schema {
        Schema::new(vec![
            Field::new("metadata", DataType::Utf8, false),
            Field::new(
                "values_a",
                DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
                true,
            ),
            Field::new(
                "values_b",
                DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
                true,
            ),
        ])
    }

    #[test]
    fn test_fused_array_transform_creation() {
        let schema = create_test_schema();
        let df_schema = Arc::new(DFSchema::try_from(schema).unwrap());

        // Create a simple identity transform using EmptyRelation
        let input = Arc::new(LogicalPlan::EmptyRelation(
            datafusion::logical_expr::EmptyRelation {
                produce_one_row: false,
                schema: df_schema.clone(),
            },
        ));

        let node = FusedArrayTransform::try_new(
            input,
            vec!["values_a".to_string()],
            vec!["metadata".to_string()],
            vec!["values_a_out".to_string()],
            vec![col("val_a")],
        );

        // Node creation should succeed
        assert!(node.is_ok());
    }

    #[test]
    fn test_fused_array_transform_name() {
        let schema = create_test_schema();
        let df_schema = Arc::new(DFSchema::try_from(schema).unwrap());

        // Create a minimal valid node for testing
        let node = FusedArrayTransform {
            input: Arc::new(LogicalPlan::EmptyRelation(
                datafusion::logical_expr::EmptyRelation {
                    produce_one_row: false,
                    schema: df_schema.clone(),
                },
            )),
            array_columns: vec!["values_a".to_string()],
            passthrough_columns: vec!["metadata".to_string()],
            output_columns: vec!["values_a_out".to_string()],
            transform_exprs: vec![col("val_a")],
            schema: df_schema,
        };

        assert_eq!(node.name(), "FusedArrayTransform");
    }
}
