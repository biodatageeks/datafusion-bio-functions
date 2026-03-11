//! Integration tests for the Fused Array Transform optimization.
//!
//! These tests verify the full pipeline: SQL → optimize → execute.
//!
//! NOTE: All tests use `#[serial]` because they modify the shared
//! `BIO_FUSED_ARRAY_TRANSFORM` environment variable.

use std::sync::Arc;

use datafusion::arrow::array::{
    Array, Float32Array, Float32Builder, Float64Array, Float64Builder, Int32Array, Int32Builder,
    ListArray, ListBuilder, RecordBatch, StringArray, StructArray,
};
use datafusion::arrow::datatypes::{DataType, Field, Fields, Schema};
use datafusion::datasource::MemTable;
use datafusion::execution::session_state::SessionStateBuilder;
use datafusion::prelude::*;
use datafusion_bio_function_vcftools::{FusedArrayTransformOptimizerRule, VcfQueryPlanner};
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
    use datafusion::physical_plan::displayable;

    // Run the query with the optimization enabled
    let ctx = create_optimized_context().await;
    let batch = create_test_data();
    let schema = batch.schema();

    // Register test data
    let table = MemTable::try_new(schema, vec![vec![batch]]).unwrap();
    ctx.register_table("test_data", Arc::new(table)).unwrap();

    // Simple query that should be optimized by FusedArrayTransform
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

    // Verify that the physical plan contains the fused transform
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    assert!(
        plan_str.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for simple identity transform! Physical plan:\n{plan_str}"
    );

    // Execute and verify results are still correct
    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();
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

#[tokio::test]
#[serial]
async fn test_simple_identity_baseline_vs_optimized_same_results() {
    use datafusion::arrow::util::pretty::pretty_format_batches;
    use datafusion::physical_plan::displayable;

    // Reuse the same SQL as the simple identity transform test
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

    // Baseline context: no optimization
    let ctx_baseline = create_baseline_context().await;
    let batch_baseline = create_test_data();
    let schema_baseline = batch_baseline.schema();
    let table_baseline =
        MemTable::try_new(schema_baseline.clone(), vec![vec![batch_baseline]]).unwrap();
    ctx_baseline
        .register_table("test_data", Arc::new(table_baseline))
        .unwrap();

    let df_baseline = ctx_baseline.sql(sql).await.unwrap();
    let physical_plan_baseline = df_baseline.create_physical_plan().await.unwrap();
    let plan_str_baseline = displayable(physical_plan_baseline.as_ref())
        .indent(true)
        .to_string();
    assert!(
        !plan_str_baseline.contains("FusedArrayTransform"),
        "Baseline physical plan unexpectedly contains FusedArrayTransform:\n{plan_str_baseline}"
    );
    let df_baseline2 = ctx_baseline.sql(sql).await.unwrap();
    let baseline_results = df_baseline2.collect().await.unwrap();

    // Optimized context: with FusedArrayTransform
    let ctx_optimized = create_optimized_context().await;
    let batch_optimized = create_test_data();
    let schema_optimized = batch_optimized.schema();
    let table_optimized =
        MemTable::try_new(schema_optimized.clone(), vec![vec![batch_optimized]]).unwrap();
    ctx_optimized
        .register_table("test_data", Arc::new(table_optimized))
        .unwrap();

    let df_optimized = ctx_optimized.sql(sql).await.unwrap();
    let physical_plan_optimized = df_optimized.create_physical_plan().await.unwrap();
    let plan_str_optimized = displayable(physical_plan_optimized.as_ref())
        .indent(true)
        .to_string();
    assert!(
        plan_str_optimized.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for simple identity baseline-vs-optimized test! Physical plan:\n{plan_str_optimized}"
    );
    let df_optimized2 = ctx_optimized.sql(sql).await.unwrap();
    let optimized_results = df_optimized2.collect().await.unwrap();

    // Compare pretty-printed results to ensure outputs are identical
    let baseline_str = pretty_format_batches(&baseline_results)
        .unwrap()
        .to_string();
    let optimized_str = pretty_format_batches(&optimized_results)
        .unwrap()
        .to_string();

    assert_eq!(
        baseline_str, optimized_str,
        "Baseline and optimized results differ for simple identity transform"
    );
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

    // Verify the struct field values are correctly extracted
    // Test data:
    // - Row 0: GT = [0, 1], DP = [30, 25]
    // - Row 1: GT = [1, 1], DP = [40, 35]
    // - Row 2: GT = [0, 0], DP = [20, 22]
    let batch = &results[0];
    let gt_col = batch
        .column_by_name("gt_out")
        .unwrap()
        .as_any()
        .downcast_ref::<ListArray>()
        .unwrap();
    let dp_col = batch
        .column_by_name("dp_out")
        .unwrap()
        .as_any()
        .downcast_ref::<ListArray>()
        .unwrap();

    fn get_i32_list(arr: &ListArray, row: usize) -> Vec<i32> {
        let inner = arr.value(row);
        if let Some(i32_arr) = inner.as_any().downcast_ref::<Int32Array>() {
            (0..i32_arr.len()).map(|i| i32_arr.value(i)).collect()
        } else if let Some(i64_arr) = inner
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Int64Array>()
        {
            // DataFusion may cast to Int64
            (0..i64_arr.len())
                .map(|i| i64_arr.value(i) as i32)
                .collect()
        } else {
            panic!("Unexpected array type for struct field elements");
        }
    }

    let gt0 = get_i32_list(gt_col, 0);
    let gt1 = get_i32_list(gt_col, 1);
    let gt2 = get_i32_list(gt_col, 2);
    let dp0 = get_i32_list(dp_col, 0);
    let dp1 = get_i32_list(dp_col, 1);
    let dp2 = get_i32_list(dp_col, 2);

    assert_eq!(gt0, vec![0, 1], "Row 0 GT should be [0, 1]");
    assert_eq!(gt1, vec![1, 1], "Row 1 GT should be [1, 1]");
    assert_eq!(gt2, vec![0, 0], "Row 2 GT should be [0, 0]");
    assert_eq!(dp0, vec![30, 25], "Row 0 DP should be [30, 25]");
    assert_eq!(dp1, vec![40, 35], "Row 1 DP should be [40, 35]");
    assert_eq!(dp2, vec![20, 22], "Row 2 DP should be [20, 22]");

    println!("Struct field access test passed - verified GT and DP values");

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

/// Create test data that combines NULLs and mismatched array lengths across two value columns.
fn create_test_data_with_nulls_and_mismatched_lengths() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Float64Builder::new());
    let mut list_builder_b = ListBuilder::new(Float64Builder::new());

    // Row 0: values_a = [1.0, NULL, 3.0] (len 3), values_b = [10.0, 20.0] (len 2)
    list_builder_a.values().append_value(1.0);
    list_builder_a.values().append_null();
    list_builder_a.values().append_value(3.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(10.0);
    list_builder_b.values().append_value(20.0);
    list_builder_b.append(true);

    // Row 1: values_a = [NULL, 5.0] (len 2), values_b = [40.0] (len 1)
    list_builder_a.values().append_null();
    list_builder_a.values().append_value(5.0);
    list_builder_a.append(true);

    list_builder_b.values().append_value(40.0);
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

#[tokio::test]
#[serial]
async fn test_null_arrays() {
    use datafusion::arrow::util::pretty::pretty_format_batches;
    use datafusion::physical_plan::displayable;

    let batch = create_test_data_with_null_arrays();
    let schema = batch.schema();

    // Execute query with NULL arrays - baseline and optimized should behave the same
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

    // Baseline context
    let ctx_baseline = create_baseline_context().await;
    let table_baseline = MemTable::try_new(schema.clone(), vec![vec![batch.clone()]]).unwrap();
    ctx_baseline
        .register_table("test_data", Arc::new(table_baseline))
        .unwrap();

    let df_baseline = ctx_baseline.sql(sql).await.unwrap();
    let physical_plan_baseline = df_baseline.create_physical_plan().await.unwrap();
    let plan_str_baseline = displayable(physical_plan_baseline.as_ref())
        .indent(true)
        .to_string();
    assert!(
        !plan_str_baseline.contains("FusedArrayTransform"),
        "Baseline physical plan unexpectedly contains FusedArrayTransform:\n{plan_str_baseline}"
    );
    let df_baseline2 = ctx_baseline.sql(sql).await.unwrap();
    let baseline_results = df_baseline2.collect().await.unwrap();

    // Optimized context
    let ctx_optimized = create_optimized_context().await;
    let table_optimized = MemTable::try_new(schema.clone(), vec![vec![batch]]).unwrap();
    ctx_optimized
        .register_table("test_data", Arc::new(table_optimized))
        .unwrap();

    let df_optimized = ctx_optimized.sql(sql).await.unwrap();
    let physical_plan_optimized = df_optimized.create_physical_plan().await.unwrap();
    let plan_str_optimized = displayable(physical_plan_optimized.as_ref())
        .indent(true)
        .to_string();
    assert!(
        plan_str_optimized.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for null arrays case! Physical plan:\n{plan_str_optimized}"
    );
    let df_optimized2 = ctx_optimized.sql(sql).await.unwrap();
    let optimized_results = df_optimized2.collect().await.unwrap();

    // Both baseline and optimized should filter out the NULL-array row and have identical output
    let baseline_str = pretty_format_batches(&baseline_results)
        .unwrap()
        .to_string();
    let optimized_str = pretty_format_batches(&optimized_results)
        .unwrap()
        .to_string();
    assert_eq!(
        baseline_str, optimized_str,
        "Baseline and optimized results differ for null arrays case"
    );
}

#[tokio::test]
#[serial]
async fn test_empty_arrays() {
    use datafusion::arrow::util::pretty::pretty_format_batches;
    use datafusion::physical_plan::displayable;

    let batch = create_test_data_with_empty_arrays();
    let schema = batch.schema();

    // Execute query with empty arrays - baseline and optimized should behave the same
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

    // Baseline context
    let ctx_baseline = create_baseline_context().await;
    let table_baseline = MemTable::try_new(schema.clone(), vec![vec![batch.clone()]]).unwrap();
    ctx_baseline
        .register_table("test_data", Arc::new(table_baseline))
        .unwrap();

    let df_baseline = ctx_baseline.sql(sql).await.unwrap();
    let physical_plan_baseline = df_baseline.create_physical_plan().await.unwrap();
    let plan_str_baseline = displayable(physical_plan_baseline.as_ref())
        .indent(true)
        .to_string();
    assert!(
        !plan_str_baseline.contains("FusedArrayTransform"),
        "Baseline physical plan unexpectedly contains FusedArrayTransform:\n{plan_str_baseline}"
    );
    let df_baseline2 = ctx_baseline.sql(sql).await.unwrap();
    let baseline_results = df_baseline2.collect().await.unwrap();

    // Optimized context
    let ctx_optimized = create_optimized_context().await;
    let table_optimized = MemTable::try_new(schema.clone(), vec![vec![batch]]).unwrap();
    ctx_optimized
        .register_table("test_data", Arc::new(table_optimized))
        .unwrap();

    let df_optimized = ctx_optimized.sql(sql).await.unwrap();
    let physical_plan_optimized = df_optimized.create_physical_plan().await.unwrap();
    let plan_str_optimized = displayable(physical_plan_optimized.as_ref())
        .indent(true)
        .to_string();
    assert!(
        plan_str_optimized.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for empty arrays case! Physical plan:\n{plan_str_optimized}"
    );
    let df_optimized2 = ctx_optimized.sql(sql).await.unwrap();
    let optimized_results = df_optimized2.collect().await.unwrap();

    // Both baseline and optimized should filter out the empty-array row and have identical output
    let baseline_str = pretty_format_batches(&baseline_results)
        .unwrap()
        .to_string();
    let optimized_str = pretty_format_batches(&optimized_results)
        .unwrap()
        .to_string();
    assert_eq!(
        baseline_str, optimized_str,
        "Baseline and optimized results differ for empty arrays case"
    );
}

#[tokio::test]
#[serial]
async fn test_mismatched_array_lengths() {
    use datafusion::arrow::util::pretty::pretty_format_batches;
    use datafusion::physical_plan::displayable;

    let batch = create_test_data_with_mismatched_lengths();
    let schema = batch.schema();

    // Query with unnest on multiple arrays with mismatched lengths.
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

    // Baseline context
    let ctx_baseline = create_baseline_context().await;
    let table_baseline = MemTable::try_new(schema.clone(), vec![vec![batch.clone()]]).unwrap();
    ctx_baseline
        .register_table("test_data", Arc::new(table_baseline))
        .unwrap();

    let df_baseline = ctx_baseline.sql(sql).await.unwrap();
    let physical_plan_baseline = df_baseline.create_physical_plan().await.unwrap();
    let plan_str_baseline = displayable(physical_plan_baseline.as_ref())
        .indent(true)
        .to_string();
    assert!(
        !plan_str_baseline.contains("FusedArrayTransform"),
        "Baseline physical plan unexpectedly contains FusedArrayTransform:\n{plan_str_baseline}"
    );
    let df_baseline2 = ctx_baseline.sql(sql).await.unwrap();
    let baseline_results = df_baseline2.collect().await.unwrap();

    // Optimized context
    let ctx_optimized = create_optimized_context().await;
    let table_optimized = MemTable::try_new(schema.clone(), vec![vec![batch]]).unwrap();
    ctx_optimized
        .register_table("test_data", Arc::new(table_optimized))
        .unwrap();

    let df_optimized = ctx_optimized.sql(sql).await.unwrap();
    let physical_plan_optimized = df_optimized.create_physical_plan().await.unwrap();
    let plan_str_optimized = displayable(physical_plan_optimized.as_ref())
        .indent(true)
        .to_string();
    assert!(
        plan_str_optimized.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for empty arrays case! Physical plan:\n{plan_str_optimized}"
    );
    let df_optimized2 = ctx_optimized.sql(sql).await.unwrap();
    let optimized_results = df_optimized2.collect().await.unwrap();

    // Both baseline and optimized should filter out the empty-array row and have identical output
    let baseline_str = pretty_format_batches(&baseline_results)
        .unwrap()
        .to_string();
    let optimized_str = pretty_format_batches(&optimized_results)
        .unwrap()
        .to_string();
    assert_eq!(
        baseline_str, optimized_str,
        "Baseline and optimized results differ for empty arrays case"
    );
}

#[tokio::test]
#[serial]
async fn test_null_elements_in_arrays() {
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await;
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

    // Optimization should still be applied when arrays contain NULL elements
    let physical_plan = df.create_physical_plan().await.unwrap();
    let plan_str = displayable(physical_plan.as_ref()).indent(true).to_string();
    assert!(
        plan_str.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for null-elements case! Physical plan:\n{plan_str}"
    );

    let df2 = ctx.sql(sql).await.unwrap();
    let results = df2.collect().await.unwrap();

    // Both rows should appear
    let total_rows: usize = results.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total_rows, 2, "Expected 2 rows");

    // Verify NULL elements are preserved
    // Test data:
    // - Row 0: [1.0, NULL, 3.0] -> 3 elements, 1 null
    // - Row 1: [NULL, 5.0] -> 2 elements, 1 null
    let batch = &results[0];
    let values_col = batch.column(2);
    let list_array = values_col.as_any().downcast_ref::<ListArray>().unwrap();

    // Row 0: [1.0, NULL, 3.0]
    let row0_inner = list_array.value(0);
    let row0_floats = row0_inner.as_any().downcast_ref::<Float64Array>().unwrap();
    assert_eq!(row0_floats.len(), 3, "Row 0 should have 3 elements");
    assert_eq!(row0_floats.null_count(), 1, "Row 0 should have 1 null");
    assert!(!row0_floats.is_null(0), "Row 0[0] should be valid");
    assert!(row0_floats.is_null(1), "Row 0[1] should be null");
    assert!(!row0_floats.is_null(2), "Row 0[2] should be valid");
    assert!(
        (row0_floats.value(0) - 1.0).abs() < 0.001,
        "Row 0[0] should be 1.0"
    );
    assert!(
        (row0_floats.value(2) - 3.0).abs() < 0.001,
        "Row 0[2] should be 3.0"
    );

    // Row 1: [NULL, 5.0]
    let row1_inner = list_array.value(1);
    let row1_floats = row1_inner.as_any().downcast_ref::<Float64Array>().unwrap();
    assert_eq!(row1_floats.len(), 2, "Row 1 should have 2 elements");
    assert_eq!(row1_floats.null_count(), 1, "Row 1 should have 1 null");
    assert!(row1_floats.is_null(0), "Row 1[0] should be null");
    assert!(!row1_floats.is_null(1), "Row 1[1] should be valid");
    assert!(
        (row1_floats.value(1) - 5.0).abs() < 0.001,
        "Row 1[1] should be 5.0"
    );

    println!("Null elements test passed - verified NULL positions are preserved");
}

#[tokio::test]
#[serial]
async fn test_nulls_and_mismatched_lengths_baseline_vs_optimized_same_results() {
    use datafusion::arrow::util::pretty::pretty_format_batches;
    use datafusion::physical_plan::displayable;

    // SQL that unnests both value columns and aggregates them back
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

    // Baseline execution (no optimization)
    let ctx_baseline = create_baseline_context().await;
    let batch_baseline = create_test_data_with_nulls_and_mismatched_lengths();
    let schema_baseline = batch_baseline.schema();
    let table_baseline =
        MemTable::try_new(schema_baseline.clone(), vec![vec![batch_baseline]]).unwrap();
    ctx_baseline
        .register_table("test_data", Arc::new(table_baseline))
        .unwrap();

    let df_baseline = ctx_baseline.sql(sql).await.unwrap();
    let physical_plan_baseline = df_baseline.create_physical_plan().await.unwrap();
    let plan_str_baseline = displayable(physical_plan_baseline.as_ref())
        .indent(true)
        .to_string();
    assert!(
        !plan_str_baseline.contains("FusedArrayTransform"),
        "Baseline physical plan unexpectedly contains FusedArrayTransform:\n{plan_str_baseline}"
    );
    let df_baseline2 = ctx_baseline.sql(sql).await.unwrap();
    let baseline_results = df_baseline2.collect().await.unwrap();

    // Optimized execution (with FusedArrayTransform)
    let ctx_optimized = create_optimized_context().await;
    let batch_optimized = create_test_data_with_nulls_and_mismatched_lengths();
    let schema_optimized = batch_optimized.schema();
    let table_optimized =
        MemTable::try_new(schema_optimized.clone(), vec![vec![batch_optimized]]).unwrap();
    ctx_optimized
        .register_table("test_data", Arc::new(table_optimized))
        .unwrap();

    let df_optimized = ctx_optimized.sql(sql).await.unwrap();
    let physical_plan_optimized = df_optimized.create_physical_plan().await.unwrap();
    let plan_str_optimized = displayable(physical_plan_optimized.as_ref())
        .indent(true)
        .to_string();
    assert!(
        plan_str_optimized.contains("FusedArrayTransform"),
        "FusedArrayTransform optimization was NOT applied for nulls+mismatched-lengths case! Physical plan:\n{plan_str_optimized}"
    );
    let df_optimized2 = ctx_optimized.sql(sql).await.unwrap();
    let optimized_results = df_optimized2.collect().await.unwrap();

    // Compare pretty-printed results to ensure outputs are identical
    let baseline_str = pretty_format_batches(&baseline_results)
        .unwrap()
        .to_string();
    let optimized_str = pretty_format_batches(&optimized_results)
        .unwrap()
        .to_string();

    assert_eq!(
        baseline_str, optimized_str,
        "Baseline and optimized results differ for nulls+mismatched-lengths case"
    );
}

// =============================================================================
// Expression Transform Tests
// =============================================================================

#[tokio::test]
#[serial]
async fn test_arithmetic_transform() {
    use datafusion::physical_plan::displayable;

    let ctx = create_optimized_context().await; // Use optimized context
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
        let values_a_col = batch
            .column(2)
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap();
        let values_b_col = batch
            .column(3)
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap();
        let values_sum_col = batch
            .column(4)
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap();
        let values_product_col = batch
            .column(5)
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap();

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
                let sum_j = sum_vals.value(j);
                let product_j = product_vals.value(j);

                assert!(
                    (sum_j - expected_sum).abs() < 1e-10,
                    "Sum mismatch at row {i}, element {j}: {a} + {b} = {expected_sum}, got {sum_j}"
                );

                assert!(
                    (product_j - expected_product).abs() < 1e-10,
                    "Product mismatch at row {i}, element {j}: {a} * {b} = {expected_product}, got {product_j}"
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

    let ctx = create_optimized_context().await; // Use optimized context
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

    // Verify the CASE expression produces correct values
    // Test data:
    // - Row 0: values_a=[1.0, 2.0, 3.0], values_b=[10.0, 20.0, 30.0]
    // - Row 1: values_a=[4.0, 5.0], values_b=[40.0, 50.0]
    // - Row 2: values_a=[6.0], values_b=[60.0]
    // CASE WHEN val_a > 3.0 THEN val_a * val_b ELSE val_a + val_b END
    // Expected:
    // - Row 0: [1+10, 2+20, 3+30] = [11.0, 22.0, 33.0]
    // - Row 1: [4*40, 5*50] = [160.0, 250.0]
    // - Row 2: [6*60] = [360.0]
    let batch = &results[0];
    let values_cond_col = batch
        .column_by_name("values_conditional")
        .unwrap()
        .as_any()
        .downcast_ref::<ListArray>()
        .unwrap();

    fn get_f64_list(arr: &ListArray, row: usize) -> Vec<f64> {
        let inner = arr.value(row);
        if let Some(f64_arr) = inner.as_any().downcast_ref::<Float64Array>() {
            (0..f64_arr.len()).map(|i| f64_arr.value(i)).collect()
        } else if let Some(f32_arr) = inner.as_any().downcast_ref::<Float32Array>() {
            (0..f32_arr.len())
                .map(|i| f32_arr.value(i) as f64)
                .collect()
        } else {
            panic!("Unexpected array type for values_conditional elements");
        }
    }

    let row0 = get_f64_list(values_cond_col, 0);
    let row1 = get_f64_list(values_cond_col, 1);
    let row2 = get_f64_list(values_cond_col, 2);

    assert_eq!(row0.len(), 3, "Row 0 should have 3 elements");
    assert!(
        (row0[0] - 11.0).abs() < 0.001,
        "Row 0[0] should be 11.0, got {}",
        row0[0]
    );
    assert!(
        (row0[1] - 22.0).abs() < 0.001,
        "Row 0[1] should be 22.0, got {}",
        row0[1]
    );
    assert!(
        (row0[2] - 33.0).abs() < 0.001,
        "Row 0[2] should be 33.0, got {}",
        row0[2]
    );

    assert_eq!(row1.len(), 2, "Row 1 should have 2 elements");
    assert!(
        (row1[0] - 160.0).abs() < 0.001,
        "Row 1[0] should be 160.0, got {}",
        row1[0]
    );
    assert!(
        (row1[1] - 250.0).abs() < 0.001,
        "Row 1[1] should be 250.0, got {}",
        row1[1]
    );

    assert_eq!(row2.len(), 1, "Row 2 should have 1 element");
    assert!(
        (row2[0] - 360.0).abs() < 0.001,
        "Row 2[0] should be 360.0, got {}",
        row2[0]
    );

    println!("Conditional transform test passed - verified CASE expression values");

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
                println!("{line}");
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

// =============================================================================
// Parquet Round-Trip Test with Complex Transform
// =============================================================================

/// Create test data with array columns for the parquet round-trip test.
///
/// Schema:
/// - metadata: Utf8
/// - values_a: List<Int32>
/// - values_b: List<Float32>
/// - values_c: List<Float32>
fn create_parquet_test_data() -> RecordBatch {
    let mut list_builder_a = ListBuilder::new(Int32Builder::new());
    let mut list_builder_b = ListBuilder::new(Float32Builder::new());
    let mut list_builder_c = ListBuilder::new(Float32Builder::new());

    // Row 0: values_a=[1, 2, 3], values_b=[10.0, 20.0, 30.0], values_c=[50.0, 150.0, 80.0]
    // val_c[0]=50 <= 100 -> val_d[0] = 1 + 10.0 = 11.0
    // val_c[1]=150 > 100 -> val_d[1] = 2 * 20.0 = 40.0
    // val_c[2]=80 <= 100 -> val_d[2] = 3 + 30.0 = 33.0
    list_builder_a.values().append_value(1);
    list_builder_a.values().append_value(2);
    list_builder_a.values().append_value(3);
    list_builder_a.append(true);

    list_builder_b.values().append_value(10.0);
    list_builder_b.values().append_value(20.0);
    list_builder_b.values().append_value(30.0);
    list_builder_b.append(true);

    list_builder_c.values().append_value(50.0);
    list_builder_c.values().append_value(150.0);
    list_builder_c.values().append_value(80.0);
    list_builder_c.append(true);

    // Row 1: values_a=[4, 5], values_b=[40.0, 50.0], values_c=[200.0, 60.0]
    // val_c[0]=200 > 100 -> val_d[0] = 4 * 40.0 = 160.0
    // val_c[1]=60 <= 100 -> val_d[1] = 5 + 50.0 = 55.0
    list_builder_a.values().append_value(4);
    list_builder_a.values().append_value(5);
    list_builder_a.append(true);

    list_builder_b.values().append_value(40.0);
    list_builder_b.values().append_value(50.0);
    list_builder_b.append(true);

    list_builder_c.values().append_value(200.0);
    list_builder_c.values().append_value(60.0);
    list_builder_c.append(true);

    // Row 2: values_a=[6], values_b=[60.0], values_c=[101.0]
    // val_c[0]=101 > 100 -> val_d[0] = 6 * 60.0 = 360.0
    list_builder_a.values().append_value(6);
    list_builder_a.append(true);

    list_builder_b.values().append_value(60.0);
    list_builder_b.append(true);

    list_builder_c.values().append_value(101.0);
    list_builder_c.append(true);

    let arr_a = list_builder_a.finish();
    let arr_b = list_builder_b.finish();
    let arr_c = list_builder_c.finish();

    let metadata = StringArray::from(vec!["meta1", "meta2", "meta3"]);

    let schema = Arc::new(Schema::new(vec![
        Field::new("metadata", DataType::Utf8, false),
        Field::new(
            "values_a",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            true,
        ),
        Field::new(
            "values_b",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            true,
        ),
        Field::new(
            "values_c",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            true,
        ),
    ]));

    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(metadata),
            Arc::new(arr_a),
            Arc::new(arr_b),
            Arc::new(arr_c),
        ],
    )
    .unwrap()
}

#[tokio::test]
#[serial]
async fn test_parquet_round_trip_with_transform() {
    use datafusion::dataframe::DataFrameWriteOptions;
    use datafusion::parquet::basic::Compression;
    use datafusion::parquet::file::properties::WriterProperties;
    use datafusion::physical_plan::displayable;

    // Create temp directory for parquet files
    let temp_dir = tempfile::tempdir().unwrap();
    let input_path = temp_dir.path().join("input.parquet");
    let output_path = temp_dir.path().join("output.parquet");

    // Step 1: Create and write input parquet file
    let batch = create_parquet_test_data();
    let schema = batch.schema();

    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let file = std::fs::File::create(&input_path).unwrap();
    let mut writer =
        parquet::arrow::ArrowWriter::try_new(file, schema.clone(), Some(props)).unwrap();
    writer.write(&batch).unwrap();
    writer.close().unwrap();

    // Step 2: Create optimized context and register input
    let ctx = create_optimized_context().await;
    ctx.register_parquet("vcf_data", input_path.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    // Step 3: Run the complex query with CASE transformation
    let sql = r#"
        WITH indexed AS (
            SELECT 
                ROW_NUMBER() OVER () as row_idx,
                metadata,
                values_a,
                values_b,
                values_c
            FROM vcf_data
        ),
        unnested AS (
            SELECT 
                row_idx,
                metadata,
                ROW_NUMBER() OVER (PARTITION BY row_idx ORDER BY row_idx) as val_idx,
                unnest(values_a) as val_a,
                unnest(values_b) as val_b,
                unnest(values_c) as val_c
            FROM indexed
        ),
        transformed AS (
            SELECT
                row_idx,
                val_idx,
                metadata,
                val_a,
                val_b,
                val_c,
                CASE WHEN val_c > 100 THEN val_a * val_b ELSE val_a + val_b END AS val_d
            FROM unnested
        )
        SELECT
            row_idx,
            metadata,
            array_agg(val_a ORDER BY val_idx) AS values_a,
            array_agg(val_b ORDER BY val_idx) AS values_b,
            array_agg(val_c ORDER BY val_idx) AS values_c,
            array_agg(val_d ORDER BY val_idx) AS values_d
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

    // Step 4: Write results to output parquet file
    df2.write_parquet(
        output_path.to_str().unwrap(),
        DataFrameWriteOptions::default(),
        None,
    )
    .await
    .unwrap();

    // Step 5: Read back the output parquet and verify contents
    let ctx2 = SessionContext::new();
    ctx2.register_parquet("output", output_path.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result_df = ctx2
        .sql("SELECT * FROM output ORDER BY row_idx")
        .await
        .unwrap();
    let results = result_df.collect().await.unwrap();

    assert_eq!(results.len(), 1, "Expected 1 batch");
    let batch = &results[0];
    println!("Result schema: {:?}", batch.schema());
    println!(
        "Result columns: {:?}",
        batch
            .schema()
            .fields()
            .iter()
            .map(|f| f.name())
            .collect::<Vec<_>>()
    );
    assert_eq!(batch.num_rows(), 3, "Expected 3 rows");

    // Verify values_d schema type - must be floating-point (result of val_a * val_b or val_a + val_b)
    // This is critical: the type inference bug caused values_d to be Int32 instead of Float32/Float64
    let schema = batch.schema();
    let values_d_field = schema.field_with_name("values_d").unwrap();
    if let DataType::List(inner_field) = values_d_field.data_type() {
        let inner_type = inner_field.data_type();
        assert!(
            matches!(inner_type, DataType::Float32 | DataType::Float64),
            "values_d inner type should be Float32 or Float64, got {inner_type:?}. \
            This indicates a type inference bug in the optimizer!"
        );
    } else {
        panic!(
            "values_d should be a List type, got {:?}",
            values_d_field.data_type()
        );
    }

    // Verify metadata column (can be Utf8 or Utf8View)
    let metadata_col = batch.column_by_name("metadata").unwrap();
    let metadata_values: Vec<String> =
        if let Some(arr) = metadata_col.as_any().downcast_ref::<StringArray>() {
            (0..arr.len()).map(|i| arr.value(i).to_string()).collect()
        } else if let Some(arr) = metadata_col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringViewArray>()
        {
            (0..arr.len()).map(|i| arr.value(i).to_string()).collect()
        } else {
            panic!("Unexpected array type for metadata column");
        };
    assert_eq!(metadata_values[0], "meta1");
    assert_eq!(metadata_values[1], "meta2");
    assert_eq!(metadata_values[2], "meta3");

    // Verify values_d (the transformed column)
    // Row 0: val_c=[50, 150, 80] -> val_d=[1+10, 2*20, 3+30] = [11, 40, 33]
    // Row 1: val_c=[200, 60] -> val_d=[4*40, 5+50] = [160, 55]
    // Row 2: val_c=[101] -> val_d=[6*60] = [360]
    let values_d_col = batch
        .column_by_name("values_d")
        .unwrap()
        .as_any()
        .downcast_ref::<ListArray>()
        .unwrap();

    // Helper to extract Float values from a ListArray at row index (handles Float32 and Float64)
    fn get_f64_list(arr: &ListArray, row: usize) -> Vec<f64> {
        let inner = arr.value(row);
        if let Some(f64_arr) = inner.as_any().downcast_ref::<Float64Array>() {
            (0..f64_arr.len()).map(|i| f64_arr.value(i)).collect()
        } else if let Some(f32_arr) = inner.as_any().downcast_ref::<Float32Array>() {
            (0..f32_arr.len())
                .map(|i| f32_arr.value(i) as f64)
                .collect()
        } else {
            panic!(
                "Unexpected array type for values_d elements: {:?}",
                inner.data_type()
            );
        }
    }

    let row0_d = get_f64_list(values_d_col, 0);
    let row1_d = get_f64_list(values_d_col, 1);
    let row2_d = get_f64_list(values_d_col, 2);

    // Expected: row0=[11.0, 40.0, 33.0], row1=[160.0, 55.0], row2=[360.0]
    assert_eq!(row0_d.len(), 3, "Row 0 values_d should have 3 elements");
    assert!(
        (row0_d[0] - 11.0).abs() < 0.001,
        "Row 0 val_d[0] should be 11.0, got {}",
        row0_d[0]
    );
    assert!(
        (row0_d[1] - 40.0).abs() < 0.001,
        "Row 0 val_d[1] should be 40.0, got {}",
        row0_d[1]
    );
    assert!(
        (row0_d[2] - 33.0).abs() < 0.001,
        "Row 0 val_d[2] should be 33.0, got {}",
        row0_d[2]
    );

    assert_eq!(row1_d.len(), 2, "Row 1 values_d should have 2 elements");
    assert!(
        (row1_d[0] - 160.0).abs() < 0.001,
        "Row 1 val_d[0] should be 160.0, got {}",
        row1_d[0]
    );
    assert!(
        (row1_d[1] - 55.0).abs() < 0.001,
        "Row 1 val_d[1] should be 55.0, got {}",
        row1_d[1]
    );

    assert_eq!(row2_d.len(), 1, "Row 2 values_d should have 1 element");
    assert!(
        (row2_d[0] - 360.0).abs() < 0.001,
        "Row 2 val_d[0] should be 360.0, got {}",
        row2_d[0]
    );

    // Verify values_a is preserved correctly
    let values_a_col = batch
        .column_by_name("values_a")
        .unwrap()
        .as_any()
        .downcast_ref::<ListArray>()
        .unwrap();

    fn get_i64_list(arr: &ListArray, row: usize) -> Vec<i64> {
        let inner = arr.value(row);
        // DataFusion may cast Int32 to Int64 during aggregation
        if let Some(i64_arr) = inner
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Int64Array>()
        {
            (0..i64_arr.len()).map(|i| i64_arr.value(i)).collect()
        } else if let Some(i32_arr) = inner.as_any().downcast_ref::<Int32Array>() {
            (0..i32_arr.len())
                .map(|i| i32_arr.value(i) as i64)
                .collect()
        } else {
            panic!("Unexpected array type for values_a");
        }
    }

    let row0_a = get_i64_list(values_a_col, 0);
    let row1_a = get_i64_list(values_a_col, 1);
    let row2_a = get_i64_list(values_a_col, 2);

    assert_eq!(row0_a, vec![1, 2, 3], "Row 0 values_a should be [1, 2, 3]");
    assert_eq!(row1_a, vec![4, 5], "Row 1 values_a should be [4, 5]");
    assert_eq!(row2_a, vec![6], "Row 2 values_a should be [6]");

    // Verify values_b is preserved correctly
    let values_b_col = batch
        .column_by_name("values_b")
        .unwrap()
        .as_any()
        .downcast_ref::<ListArray>()
        .unwrap();

    let row0_b = get_f64_list(values_b_col, 0);
    let row1_b = get_f64_list(values_b_col, 1);
    let row2_b = get_f64_list(values_b_col, 2);

    assert_eq!(row0_b.len(), 3, "Row 0 values_b should have 3 elements");
    assert!(
        (row0_b[0] - 10.0).abs() < 0.001,
        "Row 0 val_b[0] should be 10.0"
    );
    assert!(
        (row0_b[1] - 20.0).abs() < 0.001,
        "Row 0 val_b[1] should be 20.0"
    );
    assert!(
        (row0_b[2] - 30.0).abs() < 0.001,
        "Row 0 val_b[2] should be 30.0"
    );
    assert_eq!(row1_b.len(), 2, "Row 1 values_b should have 2 elements");
    assert_eq!(row2_b.len(), 1, "Row 2 values_b should have 1 element");

    println!("Parquet round-trip test passed!");
    println!("Row 0 values_d: {row0_d:?}");
    println!("Row 1 values_d: {row1_d:?}");
    println!("Row 2 values_d: {row2_d:?}");

    // Clean up
    unsafe {
        std::env::remove_var("BIO_FUSED_ARRAY_TRANSFORM");
    }
}
