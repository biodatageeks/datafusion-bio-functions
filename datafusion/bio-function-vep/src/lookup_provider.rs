//! Lookup provider for `lookup_variants()` table function.
//!
//! Implements interval join between a VCF table and a variation cache table,
//! reusing the `IntervalJoinExec` from `bio-function-ranges`.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
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

        // Build the interval join query
        let query = if select_cache.is_empty() {
            format!(
                "SELECT {select_vcf} \
                 FROM `{}` AS b, `{}` AS a \
                 WHERE a.`chrom` = b.`chrom` \
                 AND CAST(a.`end` AS INTEGER) >{sign} CAST(b.`start` AS INTEGER) \
                 AND CAST(a.`start` AS INTEGER) <{sign} CAST(b.`end` AS INTEGER) \
                 AND match_allele(a.`ref`, a.`alt`, b.`allele_string`)",
                self.cache_table, self.vcf_table,
            )
        } else {
            format!(
                "SELECT {select_vcf}, {select_cache} \
                 FROM `{}` AS b, `{}` AS a \
                 WHERE a.`chrom` = b.`chrom` \
                 AND CAST(a.`end` AS INTEGER) >{sign} CAST(b.`start` AS INTEGER) \
                 AND CAST(a.`start` AS INTEGER) <{sign} CAST(b.`end` AS INTEGER) \
                 AND match_allele(a.`ref`, a.`alt`, b.`allele_string`)",
                self.cache_table, self.vcf_table,
            )
        };

        let df = self.session.sql(&query).await?;
        df.create_physical_plan().await
    }
}
