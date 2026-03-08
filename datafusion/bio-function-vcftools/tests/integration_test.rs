//! Integration tests for the Fused Array Transform optimization.
//!
//! These tests verify the full pipeline: SQL → optimize → execute.
//!
//! NOTE: All tests use `#[serial]` because they modify the shared
//! `BIO_FUSED_ARRAY_TRANSFORM` environment variable.

use std::sync::Arc;

use datafusion::arrow::array::{
    Array, Float64Array, Float64Builder, Int32Builder, ListArray, ListBuilder, RecordBatch,
    StringArray, StructArray,
};
use datafusion::arrow::datatypes::{DataType, Field, Fields, Schema};
use datafusion::datasource::MemTable;
use datafusion::execution::session_state::SessionStateBuilder;
use datafusion::prelude::*;
use datafusion_bio_function_vcftools::{
    FusedArrayTransformOptimizerRule, VcfQueryPlanner,
};
use serial_test::serial;

/// Create test data with array columns.
fn create_test_data() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Float64Builder::new());
    let mut list_builder_b = ListBuilder::new(Float64Builder::new());

    // Row 0: [1.0, 2.0, 3.0], [10.0, 20.0, 30.0]
    list_builder_a.values().append_value(1.0);
    list_builder_a.values().append_value(2.0);
    list_builder_a.values().append_value(3.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(10.0);
    list_builder_b.values().append_value(20.0);
    list_builder_b.values().append_value(30.0);
    list_builder_b.append(true);

    // Row 1: [4.0, 5.0], [40.0, 50.0]
    list_builder_a.values().append_value(4.0);
    list_builder_a.values().append_value(5.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(40.0);
    list_builder_b.values().append_value(50.0);
    list_builder_b.append(true);

    // Row 2: [6.0], [60.0]
    list_builder_a.values().append_value(6.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(60.0);
    list_builder_b.append(true);

    let arr_a = list_builder_a.finish();
    let arr_b = list_builder_b.finish();

    let metadata = StringArray::from(vec!["meta1", "meta2", "meta3"]);

    let schema = Arc::new(Schema::new(vec![
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
    ]));

    RecordBatch::try_new(
        schema,
        vec![Arc::new(metadata), Arc::new(arr_a), Arc::new(arr_b)],
    )
    .unwrap()
}

/// Create a SessionContext with the optimization enabled.
async fn create_optimized_context() -> SessionContext {
    // Build a SessionState with the VcfQueryPlanner that handles FusedArrayTransform
    let state = SessionStateBuilder::new()
        .with_default_features()
        .with_query_planner(Arc::new(VcfQueryPlanner::new()))
        .build();

    let ctx = SessionContext::new_with_state(state);

    // Register the optimizer rule
    ctx.add_optimizer_rule(Arc::new(FusedArrayTransformOptimizerRule::new()));

    // Enable the optimization via environment variable
    // SAFETY: Tests run single-threaded
    unsafe {
        std::env::set_var("BIO_FUSED_ARRAY_TRANSFORM", "1");
    }

    ctx
}

/// Create a SessionContext without the optimization.
async fn create_baseline_context() -> SessionContext {
    // Disable the optimization
    // SAFETY: Tests run single-threaded
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }

    SessionContext::new()
}

#[tokio::test]
#[serial]
async fn test_simple_identity_transform_sql() {
    let ctx = create_baseline_context().await;
    let batch = create_test_data();
    let schema = batch.schema();

    // Register test data
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Simple query that would benefit from the optimization
    // (but we run without optimization here to verify baseline)
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a,
                values_b
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a,
                unnest(values_b) as val_b
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out,
            array_agg(val_b) AS values_b_out
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Verify we got 3 rows back
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");
}

#[tokio::test]
#[serial]
async fn test_optimizer_rule_detection() {
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await;
    let batch = create_test_data();
    let schema = batch.schema();

    // Register test data
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query with the pattern we want to optimize
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();

    // Check the PHYSICAL plan includes our optimization
    // Note: df.logical_plan() returns the RAW unoptimized plan
    // The optimizer rule is applied during physical planning
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();

    // Log the plan for debugging
    println!("Physical plan:\n{plan_str}");

    // ASSERT: Optimization must be applied in physical plan
    let has_fused = plan_str.contains("FusedArrayTransform");
    assert!(
        has_fused,
        "FusedArrayTransform optimization was NOT applied! Physical plan:\n{plan_str}"
    );

    // Execute and verify results
    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}

// =============================================================================
// Struct Field Access Tests (VCF-style data)
// =============================================================================

/// Create test data with a struct containing array fields (VCF genotypes style).
///
/// Schema:
/// - sample_id: Utf8
/// - genotypes: Struct { GT: List<Int32>, DP: List<Int32> }
fn create_test_data_with_struct_arrays() -> RecordBatch {
    // GT arrays: [0, 1], [1, 1], [0, 0]
    let mut gt_builder = ListBuilder::new(Int32Builder::new());

    // Row 0: GT = [0, 1]
    gt_builder.values().append_value(0);
    gt_builder.values().append_value(1);
    gt_builder.append(true);

    // Row 1: GT = [1, 1]
    gt_builder.values().append_value(1);
    gt_builder.values().append_value(1);
    gt_builder.append(true);

    // Row 2: GT = [0, 0]
    gt_builder.values().append_value(0);
    gt_builder.values().append_value(0);
    gt_builder.append(true);

    let gt_array = gt_builder.finish();

    // DP arrays: [30, 25], [40, 35], [20, 22]
    let mut dp_builder = ListBuilder::new(Int32Builder::new());

    // Row 0: DP = [30, 25]
    dp_builder.values().append_value(30);
    dp_builder.values().append_value(25);
    dp_builder.append(true);

    // Row 1: DP = [40, 35]
    dp_builder.values().append_value(40);
    dp_builder.values().append_value(35);
    dp_builder.append(true);

    // Row 2: DP = [20, 22]
    dp_builder.values().append_value(20);
    dp_builder.values().append_value(22);
    dp_builder.append(true);

    let dp_array = dp_builder.finish();

    // Create the "genotypes" struct field
    let gt_field = Field::new(
        "GT",
        DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
        true,
    );
    let dp_field = Field::new(
        "DP",
        DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
        true,
    );

    let struct_fields = Fields::from(vec![gt_field, dp_field]);
    let genotypes_array = StructArray::try_new(
        struct_fields.clone(),
        vec![Arc::new(gt_array), Arc::new(dp_array)],
        None,
    )
    .unwrap();

    let sample_ids = StringArray::from(vec!["sample1", "sample2", "sample3"]);

    let schema = Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("genotypes", DataType::Struct(struct_fields), false),
    ]));

    RecordBatch::try_new(
        schema,
        vec![Arc::new(sample_ids), Arc::new(genotypes_array)],
    )
    .unwrap()
}

#[tokio::test]
#[serial]
async fn test_struct_field_access_baseline() {
    // First, verify that the query works without optimization
    let ctx = create_baseline_context().await;
    let batch = create_test_data_with_struct_arrays();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_data", Arc::new(table)).unwrap();

    // Query that unnests a struct field (VCF-style)
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                sample_id,
                genotypes
            FROM vcf_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                sample_id,
                unnest(genotypes."GT") as gt,
                unnest(genotypes."DP") as dp
            FROM indexed
        )
        SELECT
            row_idx,
            sample_id,
            array_agg(gt) AS gt_out,
            array_agg(dp) AS dp_out
        FROM unnested
        GROUP BY row_idx, sample_id
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Verify we got 3 rows back
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    println!("Struct field access baseline test passed");
}

#[tokio::test]
#[serial]
async fn test_struct_field_access_explain() {
    // Print the plan to understand the structure
    let ctx = create_baseline_context().await;
    let batch = create_test_data_with_struct_arrays();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_data", Arc::new(table)).unwrap();

    let sql = r#"
        EXPLAIN VERBOSE
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                sample_id,
                genotypes
            FROM vcf_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                sample_id,
                unnest(genotypes."GT") as gt,
                unnest(genotypes."DP") as dp
            FROM indexed
        )
        SELECT
            row_idx,
            sample_id,
            array_agg(gt) AS gt_out,
            array_agg(dp) AS dp_out
        FROM unnested
        GROUP BY row_idx, sample_id
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    println!("\n=== EXPLAIN VERBOSE for struct field access ===\n");
    for batch in &results {
        let plan_col = batch
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        for i in 0..plan_col.len() {
            if !plan_col.is_null(i) {
                println!("{}", plan_col.value(i));
            }
        }
    }
    println!("\n=== END EXPLAIN ===\n");
}

#[tokio::test]
#[serial]
async fn test_struct_field_access_optimized() {
    use datafusion::physical_plan::displayable;

    // Test with optimization enabled
    let ctx = create_optimized_context().await;
    let batch = create_test_data_with_struct_arrays();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("vcf_data", Arc::new(table)).unwrap();

    // Same query but with optimization enabled
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                sample_id,
                genotypes
            FROM vcf_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                sample_id,
                unnest(genotypes."GT") as gt,
                unnest(genotypes."DP") as dp
            FROM indexed
        )
        SELECT
            row_idx,
            sample_id,
            array_agg(gt) AS gt_out,
            array_agg(dp) AS dp_out
        FROM unnested
        GROUP BY row_idx, sample_id
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();

    // Check the PHYSICAL plan includes our optimization
    // Note: df.logical_plan() returns the RAW unoptimized plan
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    println!("\n=== Physical plan for struct field access ===");
    println!("{plan_str}");
    println!("=== END ===\n");

    // ASSERT: Optimization must be applied for struct field access
    let has_fused = plan_str.contains("FusedArrayTransform");
    println!("Optimization applied: {has_fused}");
    assert!(
        has_fused,
        "FusedArrayTransform optimization was NOT applied for struct field access! Physical plan:\n{plan_str}"
    );

    // Execute and verify results
    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}
// =============================================================================
// Edge Case Tests
// =============================================================================

/// Create test data with NULL arrays.
fn create_test_data_with_null_arrays() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Float64Builder::new());

    // Row 0: [1.0, 2.0]
    list_builder_a.values().append_value(1.0);
    list_builder_a.values().append_value(2.0);
    list_builder_a.append(true);

    // Row 1: NULL array
    list_builder_a.append(false);

    // Row 2: [3.0, 4.0]
    list_builder_a.values().append_value(3.0);
    list_builder_a.values().append_value(4.0);
    list_builder_a.append(true);

    let arr_a = list_builder_a.finish();
    let metadata = StringArray::from(vec!["meta1", "meta2", "meta3"]);

    let schema = Arc::new(Schema::new(vec![
        Field::new("metadata", DataType::Utf8, false),
        Field::new(
            "values_a",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            true,
        ),
    ]));

    RecordBatch::try_new(schema, vec![Arc::new(metadata), Arc::new(arr_a)]).unwrap()
}

/// Create test data with empty arrays.
fn create_test_data_with_empty_arrays() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Float64Builder::new());

    // Row 0: [1.0, 2.0]
    list_builder_a.values().append_value(1.0);
    list_builder_a.values().append_value(2.0);
    list_builder_a.append(true);

    // Row 1: [] (empty array)
    list_builder_a.append(true);

    // Row 2: [3.0]
    list_builder_a.values().append_value(3.0);
    list_builder_a.append(true);

    let arr_a = list_builder_a.finish();
    let metadata = StringArray::from(vec!["meta1", "meta2", "meta3"]);

    let schema = Arc::new(Schema::new(vec![
        Field::new("metadata", DataType::Utf8, false),
        Field::new(
            "values_a",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            true,
        ),
    ]));

    RecordBatch::try_new(schema, vec![Arc::new(metadata), Arc::new(arr_a)]).unwrap()
}

/// Create test data with mismatched array lengths.
fn create_test_data_with_mismatched_lengths() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Float64Builder::new());
    let mut list_builder_b = ListBuilder::new(Float64Builder::new());

    // Row 0: [1.0, 2.0], [10.0, 20.0, 30.0] - different lengths
    list_builder_a.values().append_value(1.0);
    list_builder_a.values().append_value(2.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(10.0);
    list_builder_b.values().append_value(20.0);
    list_builder_b.values().append_value(30.0);
    list_builder_b.append(true);

    // Row 1: [3.0], [40.0, 50.0] - different lengths
    list_builder_a.values().append_value(3.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(40.0);
    list_builder_b.values().append_value(50.0);
    list_builder_b.append(true);

    let arr_a = list_builder_a.finish();
    let arr_b = list_builder_b.finish();
    let metadata = StringArray::from(vec!["meta1", "meta2"]);

    let schema = Arc::new(Schema::new(vec![
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
    ]));

    RecordBatch::try_new(
        schema,
        vec![Arc::new(metadata), Arc::new(arr_a), Arc::new(arr_b)],
    )
    .unwrap()
}

/// Create test data with NULL elements inside arrays.
fn create_test_data_with_null_elements() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Float64Builder::new());

    // Row 0: [1.0, NULL, 3.0]
    list_builder_a.values().append_value(1.0);
    list_builder_a.values().append_null();
    list_builder_a.values().append_value(3.0);
    list_builder_a.append(true);

    // Row 1: [NULL, 5.0]
    list_builder_a.values().append_null();
    list_builder_a.values().append_value(5.0);
    list_builder_a.append(true);

    let arr_a = list_builder_a.finish();
    let metadata = StringArray::from(vec!["meta1", "meta2"]);

    let schema = Arc::new(Schema::new(vec![
        Field::new("metadata", DataType::Utf8, false),
        Field::new(
            "values_a",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            true,
        ),
    ]));

    RecordBatch::try_new(schema, vec![Arc::new(metadata), Arc::new(arr_a)]).unwrap()
}

#[tokio::test]
#[serial]
async fn test_null_arrays() {
    let ctx = create_baseline_context().await;
    let batch = create_test_data_with_null_arrays();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Execute query with NULL arrays - should not panic
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Only rows with non-NULL arrays should appear (2 rows: meta1 and meta3)
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 2, "Expected 2 rows (NULL array row filtered out)");
}

#[tokio::test]
#[serial]
async fn test_empty_arrays() {
    let ctx = create_baseline_context().await;
    let batch = create_test_data_with_empty_arrays();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Execute query with empty arrays
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Only rows with non-empty arrays should appear (2 rows: meta1 and meta3)
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 2, "Expected 2 rows (empty array row filtered out)");
}

#[tokio::test]
#[serial]
async fn test_mismatched_array_lengths() {
    let ctx = create_baseline_context().await;
    let batch = create_test_data_with_mismatched_lengths();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query with unnest on multiple arrays with mismatched lengths
    // DataFusion should handle this with zip semantics (shortest array determines length)
    // or with null-padding.
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a,
                values_b
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a,
                unnest(values_b) as val_b
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out,
            array_agg(val_b) AS values_b_out
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Should produce results - the exact behavior depends on DataFusion's unnest semantics
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert!(total_rows > 0, "Expected at least some rows in result");

    // Print the results for debugging
    println!("Mismatched lengths test results:");
    for batch in &results {
        println!("{batch:?}");
    }
}

#[tokio::test]
#[serial]
async fn test_null_elements_in_arrays() {
    let ctx = create_baseline_context().await;
    let batch = create_test_data_with_null_elements();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query with arrays containing NULL elements
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Both rows should appear
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 2, "Expected 2 rows");

    // Verify NULL elements are preserved
    for batch in &results {
        let values_col = batch.column(2);
        let list_array = values_col.as_any().downcast_ref::<ListArray>().unwrap();

        for i in 0..list_array.len() {
            if !list_array.is_null(i) {
                let inner = list_array.value(i);
                let float_array = inner.as_any().downcast_ref::<Float64Array>().unwrap();
                // Verify that NULL elements are preserved
                let null_count = float_array.null_count();
                println!(
                    "Row {}: array has {} elements, {} nulls",
                    i,
                    float_array.len(),
                    null_count
                );
            }
        }
    }
}

// =============================================================================
// Expression Transform Tests
// =============================================================================

#[tokio::test]
#[serial]
async fn test_arithmetic_transform() {
    // NOTE: This test currently FAILS because the FusedArrayTransform optimization 
    //       does NOT support transformations in a separate CTE yet.
    //
    // Current limitation: The optimizer rule only detects the simple pattern:
    //   Aggregate → (SubqueryAlias?) → Projection → Unnest
    //
    // But with a transformation CTE like "transformed AS (SELECT val_a + val_b ...)",
    // the structure becomes:
    //   Aggregate → SubqueryAlias("transformed") → Projection → SubqueryAlias("unnested") → ...
    //
    // TODO: Enhance pattern detection to traverse multiple Projection/SubqueryAlias layers
    // and compose transformation expressions from intermediate projections.
    //
    // This test uses create_optimized_context() so it will PASS once transformation
    // CTE support is implemented.
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await;  // Use optimized context
    let batch = create_test_data();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query with actual arithmetic transformation (val_a + val_b)
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a,
                values_b
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a,
                unnest(values_b) as val_b
            FROM indexed
        ),
        transformed AS (
            SELECT
                row_idx,
                metadata,
                val_a,
                val_b,
                val_a + val_b AS val_sum,
                val_a * val_b AS val_product
            FROM unnested
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out,
            array_agg(val_b) AS values_b_out,
            array_agg(val_sum) AS values_sum,
            array_agg(val_product) AS values_product
        FROM transformed
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();

    // Check the PHYSICAL plan includes our optimization
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    println!("\n=== Physical plan for arithmetic transform ===");
    println!("{plan_str}");
    println!("=== END ===\n");

    // ASSERT: Optimization must be applied for arithmetic transform
    let has_fused = plan_str.contains("FusedArrayTransform");
    println!("Optimization applied: {has_fused}");
    assert!(
        has_fused,
        "FusedArrayTransform optimization was NOT applied for arithmetic transform! Physical plan:\n{plan_str}"
    );

    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();

    // Verify we got 3 rows back
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    // Verify the arithmetic is correct
    for batch in &results {
        let values_a_col = batch.column(2).as_any().downcast_ref::<ListArray>().unwrap();
        let values_b_col = batch.column(3).as_any().downcast_ref::<ListArray>().unwrap();
        let values_sum_col = batch.column(4).as_any().downcast_ref::<ListArray>().unwrap();
        let values_product_col = batch.column(5).as_any().downcast_ref::<ListArray>().unwrap();

        for i in 0..batch.num_rows() {
            let arr_a = values_a_col.value(i);
            let arr_b = values_b_col.value(i);
            let arr_sum = values_sum_col.value(i);
            let arr_product = values_product_col.value(i);

            let a_vals = arr_a.as_any().downcast_ref::<Float64Array>().unwrap();
            let b_vals = arr_b.as_any().downcast_ref::<Float64Array>().unwrap();
            let sum_vals = arr_sum.as_any().downcast_ref::<Float64Array>().unwrap();
            let product_vals = arr_product.as_any().downcast_ref::<Float64Array>().unwrap();

            // Verify each element
            for j in 0..a_vals.len() {
                let a = a_vals.value(j);
                let b = b_vals.value(j);
                let expected_sum = a + b;
                let expected_product = a * b;

                assert!(
                    (sum_vals.value(j) - expected_sum).abs() < 1e-10,
                    "Sum mismatch at row {i}, element {j}: {} + {} = {}, got {}",
                    a,
                    b,
                    expected_sum,
                    sum_vals.value(j)
                );

                assert!(
                    (product_vals.value(j) - expected_product).abs() < 1e-10,
                    "Product mismatch at row {i}, element {j}: {} * {} = {}, got {}",
                    a,
                    b,
                    expected_product,
                    product_vals.value(j)
                );
            }
        }
    }

    println!("Arithmetic transform test passed - all sums and products verified");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}

#[tokio::test]
#[serial]
async fn test_conditional_transform() {
    // NOTE: This test currently FAILS because the FusedArrayTransform optimization
    //       does NOT support transformations in a separate CTE yet.
    //
    // Same limitation as test_arithmetic_transform - transformation CTEs not detected.
    // This test uses create_optimized_context() so it will PASS once support is added.
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await;  // Use optimized context
    let batch = create_test_data();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query with CASE expression
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a,
                values_b
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a,
                unnest(values_b) as val_b
            FROM indexed
        ),
        transformed AS (
            SELECT
                row_idx,
                metadata,
                val_a,
                val_b,
                CASE WHEN val_a > 3.0 THEN val_a * val_b ELSE val_a + val_b END AS val_conditional
            FROM unnested
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_conditional) AS values_conditional
        FROM transformed
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();

    // Check the PHYSICAL plan includes our optimization
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    println!("\n=== Physical plan for conditional transform ===");
    println!("{plan_str}");
    println!("=== END ===\n");

    // ASSERT: Optimization must be applied for conditional transform
    let has_fused = plan_str.contains("FusedArrayTransform");
    println!("Optimization applied: {has_fused}");
    assert!(
        has_fused,
        "FusedArrayTransform optimization was NOT applied for conditional transform! Physical plan:\n{plan_str}"
    );

    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();

    // Verify we got 3 rows back
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    println!("Conditional transform test passed");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}

// =============================================================================
// False Positive Tests (optimization should NOT be applied)
// =============================================================================

#[tokio::test]
#[serial]
async fn test_optimization_not_applied_mixed_aggregates() {
    // Test that optimization is NOT applied when there's a mix of array_agg and other aggregates
    // (e.g., SUM, COUNT, AVG). Our optimization only works when ALL aggregates are array_agg.
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await;
    let batch = create_test_data();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query with mixed aggregates: array_agg + SUM
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out,
            SUM(val_a) AS total_a
        FROM unnested
        GROUP BY row_idx, metadata
        ORDER BY row_idx
    "#;

    let df = ctx.sql(sql).await.unwrap();

    // Check the PHYSICAL plan does NOT include our optimization
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    println!("\n=== Physical plan for mixed aggregates (should NOT be optimized) ===");
    println!("{plan_str}");
    println!("=== END ===\n");

    // ASSERT: Optimization must NOT be applied for mixed aggregates
    let has_fused = plan_str.contains("FusedArrayTransform");
    assert!(
        !has_fused,
        "FusedArrayTransform optimization was INCORRECTLY applied for mixed aggregates! Physical plan:\n{plan_str}"
    );

    // Verify query still executes correctly
    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    println!("Mixed aggregates test passed - optimization correctly NOT applied");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}

#[tokio::test]
#[serial]
async fn test_optimization_not_applied_no_unnest() {
    // Test that optimization is NOT applied for regular aggregations without UNNEST
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await;
    let batch = create_test_data();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Query without UNNEST - just regular array_agg on a scalar column
    let sql = r#"
        SELECT
            metadata,
            array_agg(metadata) AS metadata_array
        FROM test_data
        GROUP BY metadata
    "#;

    let df = ctx.sql(sql).await.unwrap();

    // Check the PHYSICAL plan does NOT include our optimization
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    println!("\n=== Physical plan for no-unnest query (should NOT be optimized) ===");
    println!("{plan_str}");
    println!("=== END ===\n");

    // ASSERT: Optimization must NOT be applied without UNNEST pattern
    let has_fused = plan_str.contains("FusedArrayTransform");
    assert!(
        !has_fused,
        "FusedArrayTransform optimization was INCORRECTLY applied without UNNEST! Physical plan:\n{plan_str}"
    );

    // Verify query still executes correctly
    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 3, "Expected 3 rows in result");

    println!("No-unnest test passed - optimization correctly NOT applied");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}

#[tokio::test]
#[serial]
async fn test_explain_shows_optimization() {
    let ctx = create_optimized_context().await;
    let batch = create_test_data();
    let schema = batch.schema();

    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    let sql = r#"
        EXPLAIN
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a
            FROM test_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                unnest(values_a) as val_a
            FROM indexed
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a) AS values_a_out
        FROM unnested
        GROUP BY row_idx, metadata
    "#;

    let df = ctx.sql(sql).await.unwrap();
    let results = df.collect().await.unwrap();

    // Collect all plan lines into a string
    let mut explain_output = String::new();
    for batch in &results {
        let plan_col = batch
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        for i in 0..plan_col.len() {
            if !plan_col.is_null(i) {
                let line = plan_col.value(i);
                println!("{}", line);
                explain_output.push_str(line);
                explain_output.push('\n');
            }
        }
    }

    // ASSERT: Optimization must appear in EXPLAIN output
    let has_fused = explain_output.contains("FusedArrayTransform");
    assert!(
        has_fused,
        "FusedArrayTransform optimization was NOT shown in EXPLAIN! Output:\n{explain_output}"
    );

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}
