//! FusedArrayTransform logical node.
//!
//! This module defines the `FusedArrayTransform` custom logical node that represents
//! an optimized fusion of `unnest → transform → array_agg` operations. Instead of
//! materializing all unnested rows, it processes arrays element-wise with bounded memory.

use std::cmp::Ordering;
use std::fmt::{self, Debug};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DFSchema, DFSchemaRef, Result};
use datafusion::logical_expr::{Expr, LogicalPlan, UserDefinedLogicalNodeCore};

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

        // Get the element type from the first array column as a fallback
        let fallback_element_type = if let Some(first_array_col) = array_columns.first() {
            let input_field = input_schema.field_with_unqualified_name(first_array_col)?;
            match input_field.data_type() {
                DataType::List(inner) | DataType::LargeList(inner) => inner.data_type().clone(),
                DataType::FixedSizeList(inner, _) => inner.data_type().clone(),
                other => {
                    return Err(datafusion::common::DataFusionError::Plan(format!(
                        "Column '{first_array_col}' has type {other}, expected List type"
                    )));
                }
            }
        } else {
            DataType::Float64 // Default fallback
        };

        // Add output array columns for ALL output columns (not just those with matching array columns)
        // For transformation CTEs, we may have more outputs than input arrays
        for (idx, output_col) in output_columns.iter().enumerate() {
            // Try to get element type from corresponding array column if it exists
            let element_type = if idx < array_columns.len() {
                let array_col = &array_columns[idx];
                let input_field = input_schema.field_with_unqualified_name(array_col)?;
                match input_field.data_type() {
                    DataType::List(inner) | DataType::LargeList(inner) => inner.data_type().clone(),
                    DataType::FixedSizeList(inner, _) => inner.data_type().clone(),
                    other => {
                        return Err(datafusion::common::DataFusionError::Plan(format!(
                            "Column '{array_col}' has type {other}, expected List type"
                        )));
                    }
                }
            } else {
                // For derived columns (e.g., val_sum = val_a + val_b), use fallback type
                fallback_element_type.clone()
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
            transform_exprs,
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
        // Note: Expr doesn't implement Hash, so we hash its debug representation
        for expr in &self.transform_exprs {
            format!("{expr:?}").hash(state);
        }
    }
}

impl PartialEq for FusedArrayTransform {
    fn eq(&self, other: &Self) -> bool {
        self.array_columns == other.array_columns
            && self.passthrough_columns == other.passthrough_columns
            && self.output_columns == other.output_columns
            && self.transform_exprs.len() == other.transform_exprs.len()
            && self
                .transform_exprs
                .iter()
                .zip(other.transform_exprs.iter())
                .all(|(a, b)| format!("{a:?}") == format!("{b:?}"))
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

    fn with_exprs_and_inputs(
        &self,
        exprs: Vec<Expr>,
        inputs: Vec<LogicalPlan>,
    ) -> Result<Self> {
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
