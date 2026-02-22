//! Lookup provider for `lookup_variants()` table function.
//!
//! Implements interval join between a VCF table and a variation cache table,
//! reusing the `IntervalJoinExec` from `bio-function-ranges`.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::array::{Array, StringArray};
use datafusion::arrow::datatypes::{Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::Result;
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
    cache_schema: Schema,
    /// Columns to select from the cache table.
    cache_columns: Vec<String>,
    /// Whether to auto-prune all-null columns from output.
    #[allow(dead_code)]
    prune_nulls: bool,
    /// Coordinate normalizer for handling different coordinate systems.
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
            cache_schema,
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
        _projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        // Detect whether VCF uses "chr"-prefixed chromosome names.
        // Ensembl VEP cache always uses bare names (e.g. "1", "22").
        let vcf_has_chr = has_chr_prefix(&self.session, &self.vcf_table).await?;

        // Build the VCF FROM source — strip "chr" prefix if needed so the
        // equi-join on chrom matches the cache's bare chromosome names.
        let vcf_from = if vcf_has_chr {
            let inner_cols = self
                .vcf_schema
                .fields()
                .iter()
                .map(|f| {
                    let name = f.name();
                    if name == "chrom" {
                        "replace(`chrom`, 'chr', '') AS `chrom`".to_string()
                    } else {
                        format!("`{name}`")
                    }
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!(
                "(SELECT {inner_cols} FROM `{}`) AS a",
                self.vcf_table
            )
        } else {
            format!("`{}` AS a", self.vcf_table)
        };

        // Build SELECT list for VCF columns
        let select_vcf = self
            .vcf_schema
            .fields()
            .iter()
            .map(|f| {
                let name = f.name();
                format!("a.`{name}` AS `{name}`")
            })
            .collect::<Vec<_>>()
            .join(", ");

        // Build SELECT list for cache columns
        let select_cache = self
            .cache_columns
            .iter()
            .filter_map(|col_name| {
                self.cache_schema
                    .field_with_name(col_name)
                    .ok()
                    .map(|_| format!("b.`{col_name}` AS `cache_{col_name}`"))
            })
            .collect::<Vec<_>>()
            .join(", ");

        // Determine overlap sign based on coordinate systems
        let sign = if self.coord_normalizer.same_system() {
            "="
        } else {
            ""
        };

        // Build the ON clause for the interval join
        let on_clause = format!(
            "a.`chrom` = b.`chrom` \
             AND CAST(a.`end` AS INTEGER) >{sign} CAST(b.`start` AS INTEGER) \
             AND CAST(a.`start` AS INTEGER) <{sign} CAST(b.`end` AS INTEGER) \
             AND match_allele(a.`ref`, a.`alt`, b.`allele_string`)"
        );

        // Build the LEFT JOIN query.
        // VCF is the left (build) side — small, fully indexed in memory.
        // Cache is the right (probe) side — potentially huge, streamed.
        // LEFT JOIN ensures all VCF variants appear in output.
        let query = if select_cache.is_empty() {
            format!(
                "SELECT {select_vcf} \
                 FROM {vcf_from} LEFT JOIN `{}` AS b \
                 ON {on_clause}",
                self.cache_table,
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
                "SELECT {select_vcf}, {select_cache} \
                 FROM {vcf_from} LEFT JOIN \
                 (SELECT * FROM `{}` WHERE {null_filter}) AS b \
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
    use datafusion::arrow::array::{Array, Int64Array, RecordBatch, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use datafusion::physical_plan::displayable;
    use std::sync::Arc;

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
}

