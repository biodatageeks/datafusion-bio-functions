//! Lookup provider for `lookup_variants()` table function.
//!
//! Implements interval join between a VCF table and a variation cache table,
//! reusing the `IntervalJoinExec` from `bio-function-ranges`.

use std::any::Any;
use std::collections::BTreeSet;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::array::{Array, StringArray};
use datafusion::arrow::datatypes::{Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{Expr, SessionContext};

use crate::coordinate::CoordinateNormalizer;
use crate::schema_contract::validate_variation_schema;

/// Table provider that implements variant lookup via interval join.
///
/// Generates an internal SQL plan that joins VCF variants against the variation
/// cache using interval overlap on (chrom, start, end), with allele matching
/// as a post-filter.
pub struct LookupProvider {
    session: Arc<SessionContext>,
    vcf_table: String,
    cache_table: String,
    vcf_schema: Schema,
    /// Columns to select from the cache table.
    cache_columns: Vec<String>,
    /// Whether to auto-prune all-null columns from output.
    #[allow(dead_code)]
    prune_nulls: bool,
    /// Coordinate normalizer for handling different coordinate systems.
    #[allow(dead_code)]
    coord_normalizer: CoordinateNormalizer,
    /// Output schema.
    schema: SchemaRef,
}

impl LookupProvider {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        session: Arc<SessionContext>,
        vcf_table: String,
        cache_table: String,
        vcf_schema: Schema,
        cache_schema: Schema,
        cache_columns: Vec<String>,
        prune_nulls: bool,
    ) -> Result<Self> {
        let cache_schema_ref: SchemaRef = Arc::new(cache_schema.clone());
        validate_variation_schema(&cache_schema_ref)?;

        let vcf_schema_ref: SchemaRef = Arc::new(vcf_schema.clone());
        let coord_normalizer =
            CoordinateNormalizer::from_schemas(&vcf_schema_ref, &cache_schema_ref);

        // Build output schema: all VCF columns + selected cache columns (prefixed with cache_)
        let mut fields: Vec<Arc<Field>> = vcf_schema.fields().iter().cloned().collect();
        for col_name in &cache_columns {
            if let Ok(field) = cache_schema.field_with_name(col_name) {
                fields.push(Arc::new(Field::new(
                    format!("cache_{}", field.name()),
                    field.data_type().clone(),
                    true, // all cache columns are nullable in output (may not match)
                )));
            }
        }
        let schema = Arc::new(Schema::new(fields));

        Ok(Self {
            session,
            vcf_table,
            cache_table,
            vcf_schema,
            cache_columns,
            prune_nulls,
            coord_normalizer,
            schema,
        })
    }
}

/// Check whether the chrom column in the given table uses "chr" prefix (e.g. "chr1").
async fn has_chr_prefix(session: &SessionContext, table: &str) -> Result<bool> {
    let batches = session
        .sql(&format!("SELECT `chrom` FROM `{table}` LIMIT 1"))
        .await?
        .collect()
        .await?;
    if let Some(batch) = batches.first() {
        if batch.num_rows() > 0 {
            if let Some(arr) = batch.column(0).as_any().downcast_ref::<StringArray>() {
                if !arr.is_null(0) {
                    return Ok(arr.value(0).starts_with("chr"));
                }
            }
        }
    }
    Ok(false)
}

impl Debug for LookupProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "LookupProvider {{ vcf: {}, cache: {}, columns: {:?} }}",
            self.vcf_table, self.cache_table, self.cache_columns
        )
    }
}

#[async_trait]
impl TableProvider for LookupProvider {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn table_type(&self) -> TableType {
        TableType::Temporary
    }

    async fn scan(
        &self,
        _state: &dyn Session,
        projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let total_output_fields = self.schema.fields().len();
        let vcf_output_fields = self.vcf_schema.fields().len();
        let projected_indices: Vec<usize> = projection
            .cloned()
            .unwrap_or_else(|| (0..total_output_fields).collect());

        let mut projected_output_exprs = Vec::with_capacity(projected_indices.len());
        let mut required_vcf_columns = BTreeSet::<String>::new();

        for idx in projected_indices {
            if idx >= total_output_fields {
                return Err(DataFusionError::Execution(format!(
                    "lookup_variants(): projection index {idx} out of bounds for schema with {total_output_fields} fields"
                )));
            }

            if idx < vcf_output_fields {
                let name = self.schema.field(idx).name().clone();
                projected_output_exprs.push(format!("a.`{name}` AS `{name}`"));
                required_vcf_columns.insert(name);
            } else {
                // Output schema appends cache columns in the same order as `self.cache_columns`.
                let cache_idx = idx - vcf_output_fields;
                let cache_col_name = self.cache_columns.get(cache_idx).ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "lookup_variants(): cache projection index {cache_idx} out of bounds for {} cache columns",
                        self.cache_columns.len()
                    ))
                })?;
                projected_output_exprs
                    .push(format!("b.`{cache_col_name}` AS `cache_{cache_col_name}`"));
            }
        }

        // Join and allele-match columns are always required from the VCF side.
        for col in ["chrom", "start", "end", "ref", "alt"] {
            required_vcf_columns.insert(col.to_string());
        }

        // Detect whether VCF uses "chr"-prefixed chromosome names.
        // Ensembl VEP cache always uses bare names (e.g. "1", "22").
        let vcf_has_chr = has_chr_prefix(&self.session, &self.vcf_table).await?;

        // Build the VCF FROM source.
        // Always project only required columns so table scans can prune aggressively.
        // Strip "chr" prefix on chrom if needed so join keys match cache naming.
        let vcf_inner_cols = required_vcf_columns
            .iter()
            .map(|name| {
                if vcf_has_chr && name == "chrom" {
                    "replace(`chrom`, 'chr', '') AS `chrom`".to_string()
                } else {
                    format!("`{name}`")
                }
            })
            .collect::<Vec<_>>()
            .join(", ");
        let vcf_from = format!("(SELECT {vcf_inner_cols} FROM `{}`) AS a", self.vcf_table);

        // Build cache-side source with explicit columns only:
        // join keys + allele string + annotation columns involved in output/null-filter.
        let mut required_cache_columns = BTreeSet::<String>::from([
            "chrom".to_string(),
            "start".to_string(),
            "end".to_string(),
            "allele_string".to_string(),
        ]);
        required_cache_columns.extend(self.cache_columns.iter().cloned());

        let cache_inner_cols = required_cache_columns
            .iter()
            .map(|c| format!("`{c}`"))
            .collect::<Vec<_>>()
            .join(", ");

        // Build the ON clause for the interval join.
        //
        // Normalize both tables to 1-based closed coordinates based on
        // `bio.coordinate_system_zero_based` metadata:
        // - zero-based half-open [start, end) -> [start + 1, end]
        // - one-based closed [start, end]     -> [start, end]
        //
        // After normalization, overlap is standard inclusive overlap.
        let vcf_start_expr = if self.coord_normalizer.input_zero_based {
            "CAST(a.`start` AS BIGINT) + 1"
        } else {
            "CAST(a.`start` AS BIGINT)"
        };
        let cache_start_expr = if self.coord_normalizer.cache_zero_based {
            "CAST(b.`start` AS BIGINT) + 1"
        } else {
            "CAST(b.`start` AS BIGINT)"
        };
        let vcf_end_expr = "CAST(a.`end` AS BIGINT)";
        let cache_end_expr = "CAST(b.`end` AS BIGINT)";

        let on_clause = format!(
            "a.`chrom` = b.`chrom` \
             AND {vcf_end_expr} >= {cache_start_expr} \
             AND {vcf_start_expr} <= {cache_end_expr} \
             AND match_allele(a.`ref`, a.`alt`, b.`allele_string`)"
        );

        // Build the LEFT JOIN query.
        // VCF is the left (build) side — small, fully indexed in memory.
        // Cache is the right (probe) side — potentially huge, streamed.
        // LEFT JOIN ensures all VCF variants appear in output.
        //
        // Note: VEP caches may contain duplicate rows at the same
        // (chrom, start, end, allele_string) with different variation_name
        // values (e.g. dbSNP + COSMIC IDs). These produce ~0.03% extra
        // output rows compared to VEP (which deduplicates internally).
        // A post-join GROUP BY would fix this but turns the streaming
        // pipeline into a blocking operation with high memory usage, so
        // we accept the minor discrepancy for performance.
        let select_output = projected_output_exprs.join(", ");
        let query = if self.cache_columns.is_empty() {
            format!(
                "SELECT {select_output} \
                 FROM {vcf_from} LEFT JOIN \
                 (SELECT {cache_inner_cols} FROM `{}`) AS b \
                 ON {on_clause}",
                self.cache_table
            )
        } else {
            // Pre-filter cache rows where all selected annotation columns are NULL
            let null_filter = self
                .cache_columns
                .iter()
                .map(|c| format!("`{c}` IS NOT NULL"))
                .collect::<Vec<_>>()
                .join(" OR ");

            format!(
                "SELECT {select_output} \
                 FROM {vcf_from} LEFT JOIN \
                 (SELECT {cache_inner_cols} FROM `{}` WHERE {null_filter}) AS b \
                 ON {on_clause}",
                self.cache_table,
            )
        };

        let df = self.session.sql(&query).await?;
        df.create_physical_plan().await
    }
}

#[cfg(test)]
mod tests {
    use crate::create_vep_session;
    use datafusion::arrow::array::{
        Array, ArrayRef, Int64Array, RecordBatch, StringArray, StringViewArray,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use datafusion::physical_plan::displayable;
    use std::collections::HashMap;
    use std::sync::Arc;

    /// Extract string values from an array that may be Utf8, Utf8View, or LargeUtf8.
    fn string_values(col: &ArrayRef) -> Vec<Option<String>> {
        use datafusion::arrow::array::LargeStringArray;
        if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else {
            panic!(
                "expected Utf8, Utf8View, or LargeUtf8, got {:?}",
                col.data_type()
            );
        }
    }

    fn schema_with_coord_metadata(fields: Vec<Field>, zero_based: bool) -> Arc<Schema> {
        let mut metadata = HashMap::new();
        metadata.insert(
            "bio.coordinate_system_zero_based".to_string(),
            zero_based.to_string(),
        );
        Arc::new(Schema::new_with_metadata(fields, metadata))
    }

    fn vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                // Three VCF rows: first two match cache, third has no match
                Arc::new(StringArray::from(vec!["chr1", "chr1", "chr2"])),
                Arc::new(Int64Array::from(vec![100, 200, 500])),
                Arc::new(Int64Array::from(vec![101, 201, 501])),
                Arc::new(StringArray::from(vec!["A", "C", "G"])),
                Arc::new(StringArray::from(vec!["G", "T", "A"])),
            ],
        )
        .unwrap();
        MemTable::try_new(schema, vec![vec![batch]]).unwrap()
    }

    fn cache_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![99, 199])),
                Arc::new(Int64Array::from(vec![102, 202])),
                Arc::new(StringArray::from(vec!["rs123", "rs456"])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(StringArray::from(vec!["benign", "pathogenic"])),
            ],
        )
        .unwrap();
        MemTable::try_new(schema, vec![vec![batch]]).unwrap()
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_generates_interval_join() {
        let ctx = create_vep_session();

        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .unwrap();
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_data', 'var_cache')")
            .await
            .unwrap();

        let plan = df.create_physical_plan().await.unwrap();
        let plan_str = displayable(plan.as_ref()).indent(true).to_string();

        assert!(
            plan_str.contains("IntervalJoinExec"),
            "Expected IntervalJoinExec in plan, got:\n{plan_str}"
        );
        assert!(
            plan_str.contains("Left"),
            "Expected LEFT join type in plan, got:\n{plan_str}"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_left_join_preserves_unmatched_vcf_rows() {
        let ctx = create_vep_session();

        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .unwrap();
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_data', 'var_cache', 'clin_sig')")
            .await
            .unwrap();

        // Verify NOT NULL pre-filter is pushed down before the join
        let plan = df.clone().create_physical_plan().await.unwrap();
        let plan_str = displayable(plan.as_ref()).indent(true).to_string();
        let interval_pos = plan_str.find("IntervalJoinExec").unwrap();
        let filter_pos = plan_str.find("IS NOT NULL").unwrap();
        assert!(
            filter_pos > interval_pos,
            "Expected IS NOT NULL filter below IntervalJoinExec (pushed down before probing), got:\n{plan_str}"
        );

        let batches = df.collect().await.unwrap();

        // Count total rows — should be 3 (two matches + one unmatched VCF row)
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 3,
            "Expected 3 rows (2 matched + 1 unmatched), got {total_rows}"
        );

        // Collect all chrom values and cache_clin_sig values
        let mut chroms = Vec::new();
        let mut clin_sigs = Vec::new();
        for batch in &batches {
            let chrom_col = batch
                .column_by_name("chrom")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let clin_sig_col = batch
                .column_by_name("cache_clin_sig")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            for i in 0..batch.num_rows() {
                chroms.push(chrom_col.value(i).to_string());
                clin_sigs.push(if clin_sig_col.is_null(i) {
                    "NULL".to_string()
                } else {
                    clin_sig_col.value(i).to_string()
                });
            }
        }

        // The chr2 variant has no cache match — it should appear with NULL annotation.
        // Chrom is "2" (not "chr2") because the chr prefix was stripped for the join.
        assert!(
            chroms.contains(&"2".to_string()),
            "Expected '2' (unmatched VCF row, chr-stripped) in output, got chroms: {chroms:?}"
        );

        let chr2_idx = chroms.iter().position(|c| c == "2").unwrap();
        assert_eq!(
            clin_sigs[chr2_idx], "NULL",
            "Expected NULL clin_sig for unmatched chr2 row, got: {}",
            clin_sigs[chr2_idx]
        );
    }

    /// Regression test: VCF half-open [start, end) must NOT match VEP cache
    /// entries at the adjacent position (start=end, 1-based closed).
    ///
    /// Before the fix, VCF [100, 101) with weak overlap (>=) would match
    /// cache entries at BOTH position 100 AND position 101, producing
    /// duplicate/false rows.
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_no_false_matches_at_adjacent_positions() {
        let ctx = create_vep_session();

        // VCF: half-open intervals (end = start + 1 for SNVs)
        let vcf_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![101, 201])),
                Arc::new(StringArray::from(vec!["A", "C"])),
                Arc::new(StringArray::from(vec!["G", "T"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        // VEP cache: 1-based closed intervals (start=end for SNVs).
        // Includes TRUE matches at positions 100 and 200, PLUS decoy entries
        // at adjacent positions (101, 201) with the same ALT allele but
        // different REF — these must NOT produce output rows.
        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                // chrom:  100 match, 101 decoy, 200 match, 201 decoy
                Arc::new(StringArray::from(vec!["1", "1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 101, 200, 201])),
                Arc::new(Int64Array::from(vec![100, 101, 200, 201])),
                Arc::new(StringArray::from(vec!["rs100", "rs101", "rs200", "rs201"])),
                // True match alleles + decoy alleles (same ALT, different REF)
                Arc::new(StringArray::from(vec!["A/G", "T/G", "C/T", "G/T"])),
                Arc::new(StringArray::from(vec![
                    "benign",
                    "decoy_101",
                    "pathogenic",
                    "decoy_201",
                ])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_narrow", Arc::new(vcf)).unwrap();
        ctx.register_table("cache_narrow", Arc::new(cache)).unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_narrow', 'cache_narrow', 'variation_name,clin_sig')")
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 2,
            "Expected exactly 2 rows (one match per VCF variant, no decoy matches), got {total_rows}"
        );

        // Verify only the true matches appear (not the decoys)
        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        assert!(
            var_names.contains(&"rs100".to_string()),
            "Expected rs100 in output, got: {var_names:?}"
        );
        assert!(
            var_names.contains(&"rs200".to_string()),
            "Expected rs200 in output, got: {var_names:?}"
        );
        assert!(
            !var_names.contains(&"rs101".to_string()),
            "Decoy rs101 should NOT appear in output, got: {var_names:?}"
        );
        assert!(
            !var_names.contains(&"rs201".to_string()),
            "Decoy rs201 should NOT appear in output, got: {var_names:?}"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_one_based_metadata_matches_point_variants() {
        let ctx = create_vep_session();

        // Explicitly mark both tables as 1-based closed.
        let vcf_schema = schema_with_coord_metadata(
            vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("ref", DataType::Utf8, false),
                Field::new("alt", DataType::Utf8, false),
            ],
            false,
        );
        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        let cache_schema = schema_with_coord_metadata(
            vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
            ],
            false,
        );
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["rs100"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["benign"])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_one_based", Arc::new(vcf)).unwrap();
        ctx.register_table("cache_one_based", Arc::new(cache))
            .unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_one_based', 'cache_one_based', 'variation_name')")
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 1, "Expected a single joined row");

        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        assert_eq!(var_names, vec!["rs100"]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_zero_based_input_to_one_based_cache_uses_metadata() {
        let ctx = create_vep_session();

        // VCF is 0-based half-open: position 100 is [99, 100).
        let vcf_schema = schema_with_coord_metadata(
            vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("ref", DataType::Utf8, false),
                Field::new("alt", DataType::Utf8, false),
            ],
            true,
        );
        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![99])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        // Cache is 1-based closed. Both rows share allele_string, so the
        // interval condition must discriminate the true position.
        let cache_schema = schema_with_coord_metadata(
            vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
            ],
            false,
        );
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 101])),
                Arc::new(Int64Array::from(vec![100, 101])),
                Arc::new(StringArray::from(vec!["rs100", "rs101"])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["benign", "benign"])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_zero_based", Arc::new(vcf)).unwrap();
        ctx.register_table("cache_closed", Arc::new(cache)).unwrap();

        let df = ctx
            .sql(
                "SELECT * FROM lookup_variants('vcf_zero_based', 'cache_closed', 'variation_name')",
            )
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 1,
            "Expected exactly one match after 0-based->1-based normalization"
        );

        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        assert_eq!(var_names, vec!["rs100"]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_one_based_input_to_zero_based_cache_uses_metadata() {
        let ctx = create_vep_session();

        // VCF is 1-based closed: position 100 is [100, 100].
        let vcf_schema = schema_with_coord_metadata(
            vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("ref", DataType::Utf8, false),
                Field::new("alt", DataType::Utf8, false),
            ],
            false,
        );
        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        // Cache is 0-based half-open:
        // pos 100 -> [99, 100), pos 101 -> [100, 101).
        let cache_schema = schema_with_coord_metadata(
            vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::Int64, false),
                Field::new("end", DataType::Int64, false),
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("clin_sig", DataType::Utf8, true),
            ],
            true,
        );
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![99, 100])),
                Arc::new(Int64Array::from(vec![100, 101])),
                Arc::new(StringArray::from(vec!["rs100", "rs101"])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["benign", "benign"])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_closed", Arc::new(vcf)).unwrap();
        ctx.register_table("cache_zero_based", Arc::new(cache))
            .unwrap();

        let df = ctx
            .sql(
                "SELECT * FROM lookup_variants('vcf_closed', 'cache_zero_based', 'variation_name')",
            )
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 1,
            "Expected exactly one match after 1-based->0-based normalization"
        );

        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        assert_eq!(var_names, vec!["rs100"]);
    }

    /// Duplicate cache rows on (chrom, start, end, allele_string) with different
    /// variation_names each produce a separate output row. This is a conscious
    /// trade-off: a post-join GROUP BY would match VEP's one-row-per-variant
    /// semantics but turns the streaming pipeline into a blocking operation.
    /// The ~0.03% extra rows from cache duplicates are acceptable.
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_cache_duplicates_produce_separate_rows() {
        let ctx = create_vep_session();

        let vcf_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let vcf_batch = RecordBatch::try_new(
            vcf_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        // Cache with THREE duplicate rows at the same position + allele.
        // Only variation_name differs (e.g. dbSNP + COSMIC at same position).
        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(StringArray::from(vec!["rs100", "COSM123", "rs100_dup"])),
                Arc::new(StringArray::from(vec!["A/G", "A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["benign", "benign", "benign"])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_dedup", Arc::new(vcf)).unwrap();
        ctx.register_table("cache_dedup", Arc::new(cache)).unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_dedup', 'cache_dedup', 'variation_name,clin_sig')")
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        // Each cache duplicate produces a separate row (no post-join aggregation)
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 3,
            "Expected 3 rows (one per cache duplicate), got {total_rows}"
        );

        // All three variation_names should appear
        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        var_names.sort();
        assert_eq!(
            var_names,
            vec!["COSM123", "rs100", "rs100_dup"],
            "Expected all three variation_names in output"
        );
    }
}
