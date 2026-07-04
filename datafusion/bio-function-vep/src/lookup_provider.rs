//! Lookup provider for the Parquet variation cache.
//!
//! Builds a `KvLookupExec` (Parquet backend) over the per-chrom variation cache;
//! consumed by the `annotate_vep` annotation path.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::path::PathBuf;
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::array::{Array, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{Expr, SessionContext};

use crate::colocated::ColocatedSink;
use crate::coordinate::CoordinateNormalizer;
use crate::schema_contract::validate_variation_schema;

/// Wrap an execution plan with a `ProjectionExec` if projection is requested.
///
/// DataFusion's physical planner expects `scan()` to return a plan whose schema
/// matches the projected schema. When a projection is provided, we wrap the
/// full-schema plan with a `ProjectionExec` that selects only the projected columns.
fn wrap_with_projection(
    plan: Arc<dyn ExecutionPlan>,
    projection: Option<&Vec<usize>>,
) -> Result<Arc<dyn ExecutionPlan>> {
    use datafusion::physical_expr::expressions::Column;
    use datafusion::physical_plan::projection::ProjectionExec;

    match projection {
        Some(indices) if indices.len() < plan.schema().fields().len() => {
            let schema = plan.schema();
            let exprs: Vec<_> = indices
                .iter()
                .map(|&idx| {
                    let field = schema.field(idx);
                    let expr = Arc::new(Column::new(field.name(), idx))
                        as Arc<dyn datafusion::physical_plan::PhysicalExpr>;
                    (expr, field.name().clone())
                })
                .collect();
            Ok(Arc::new(ProjectionExec::try_new(exprs, plan)?))
        }
        _ => Ok(plan),
    }
}

/// Table provider that implements variant lookup against the Parquet variation
/// cache via `KvLookupExec`.
///
/// VCF variants are streamed and probed against the per-chrom Parquet variation
/// dataset; `match_allele()` is applied as a post-filter and unmatched VCF rows
/// appear with NULL cache columns. `extended_probes = true` uses interval-overlap
/// matching for VEP-style coordinate encodings (insertion-style start > end,
/// prefix-trimmed deletion shifts, tandem-repeat right-shifting).
pub struct LookupProvider {
    session: Arc<SessionContext>,
    vcf_table: String,
    cache_table: String,
    /// Columns to select from the cache table.
    cache_columns: Vec<String>,
    /// Full variation cache schema.
    cache_schema: SchemaRef,
    /// Coordinate normalizer for handling different coordinate systems.
    coord_normalizer: CoordinateNormalizer,
    /// When true, use interval-overlap matching instead of exact coordinate match.
    extended_probes: bool,
    /// Maximum allowed `failed` flag value from the cache.
    allowed_failed: i64,
    /// Output schema.
    schema: SchemaRef,
    /// Optional sink for co-located data collection during probe phase.
    colocated_sink: Option<ColocatedSink>,
    /// Optional partition-local sinks for co-located data collection.
    partition_colocated_sinks: Option<Vec<ColocatedSink>>,
    /// Root of a Parquet variation cache. When set, variation lookup uses
    /// `variation.cache/chrN.cache` plus the sidecar position/bloom indexes.
    #[cfg(feature = "parquet-cache")]
    cache_root: Option<PathBuf>,
    /// When true, the variation lookup uses the Parquet backend
    /// (`variation/chrN.parquet`) instead of Parquet. The cache root is
    /// still carried in `cache_root` (the shared cache base dir).
    #[cfg(feature = "parquet-cache")]
    parquet_backend: bool,
    /// Maximum number of independent cold readers used by lookup.
    target_partitions: usize,
    /// Optional filter to apply to the VCF input (e.g., `chrom = 'chr1'`
    /// for per-contig partitioned annotation).
    vcf_filter: Option<Expr>,
}

fn normalize_cache_output_type(data_type: &DataType) -> DataType {
    match data_type {
        DataType::Utf8View | DataType::LargeUtf8 => DataType::Utf8,
        other => other.clone(),
    }
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
        extended_probes: bool,
        allowed_failed: i64,
    ) -> Result<Self> {
        let cache_schema_ref: SchemaRef = Arc::new(cache_schema.clone());
        validate_variation_schema(&cache_schema_ref)?;

        let vcf_schema_ref: SchemaRef = Arc::new(vcf_schema.clone());
        let coord_normalizer =
            CoordinateNormalizer::from_schemas(&vcf_schema_ref, &cache_schema_ref);

        // Build output schema: all VCF columns + selected cache columns (prefixed with cache_).
        let mut fields: Vec<Arc<Field>> = vcf_schema
            .fields()
            .iter()
            .map(|field| {
                Arc::new(Field::new(
                    field.name(),
                    field.data_type().clone(),
                    field.is_nullable(),
                ))
            })
            .collect();
        for col_name in &cache_columns {
            if let Ok(field) = cache_schema.field_with_name(col_name) {
                fields.push(Arc::new(Field::new(
                    format!("cache_{}", field.name()),
                    normalize_cache_output_type(field.data_type()),
                    true, // all cache columns are nullable in output (may not match)
                )));
            }
        }
        let schema = Arc::new(Schema::new(fields));

        Ok(Self {
            session,
            vcf_table,
            cache_table,
            cache_columns,
            cache_schema: cache_schema_ref,
            coord_normalizer,
            extended_probes,
            allowed_failed,
            schema,
            colocated_sink: None,
            partition_colocated_sinks: None,
            #[cfg(feature = "parquet-cache")]
            cache_root: None,
            #[cfg(feature = "parquet-cache")]
            parquet_backend: false,
            target_partitions: 1,
            vcf_filter: None,
        })
    }

    /// Set the co-located data sink for piggybacked collection during probe.
    pub fn set_colocated_sink(&mut self, sink: ColocatedSink) {
        self.colocated_sink = Some(sink);
    }

    /// Set partition-local co-located data sinks for cache lookup execution.
    pub fn set_partition_colocated_sinks(&mut self, sinks: Vec<ColocatedSink>) {
        self.partition_colocated_sinks = Some(sinks);
    }

    #[cfg(feature = "parquet-cache")]
    pub fn set_cache_root(&mut self, root: impl Into<PathBuf>) {
        self.cache_root = Some(root.into());
    }

    /// Select the Parquet variation backend (default is Parquet). The cache root is
    /// set separately via [`Self::set_cache_root`] (the shared base dir).
    #[cfg(feature = "parquet-cache")]
    pub fn set_parquet_backend(&mut self, parquet: bool) {
        self.parquet_backend = parquet;
    }

    pub fn set_target_partitions(&mut self, target_partitions: usize) {
        self.target_partitions = target_partitions.max(1);
    }

    /// Set an optional filter to apply to VCF input before lookup.
    ///
    /// Used by the partitioned annotation path to scope VCF reading to a
    /// single contig via filter pushdown (`chrom = 'chrN'`), which triggers
    /// tabix-indexed seeks in the VCF table provider.
    pub fn set_vcf_filter(&mut self, filter: Option<Expr>) {
        self.vcf_filter = filter;
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
        // Parquet cache dispatch: when a variation cache root is set, use the
        // KvLookupExec (Parquet backend) instead of the interval join.
        #[cfg(feature = "parquet-cache")]
        if let Some(cache_root) = &self.cache_root {
            use crate::allele::allele_matches;
            use crate::cache::lookup_exec::{KvLookupExec, KvMatchMode};

            let _ = self.parquet_backend;
            let vcf_has_chr = has_chr_prefix(&self.session, &self.vcf_table).await?;
            let vcf_df = self.session.table(&self.vcf_table).await?;
            let vcf_df = if let Some(ref filter) = self.vcf_filter {
                vcf_df.filter(filter.clone())?
            } else {
                vcf_df
            };
            let vcf_plan = vcf_df.create_physical_plan().await?;

            let mut exec = KvLookupExec::new_parquet(
                vcf_plan,
                cache_root.clone(),
                self.cache_schema.clone(),
                self.cache_columns.clone(),
                KvMatchMode::Exact,
                allele_matches as fn(&str, &str, &str) -> bool,
                vcf_has_chr,
                self.coord_normalizer.input_zero_based,
                self.coord_normalizer.cache_zero_based,
                self.extended_probes,
                self.allowed_failed,
            )?;
            exec = exec.with_target_partitions(self.target_partitions);
            if let Some(ref sink) = self.colocated_sink {
                exec = exec.with_colocated_sink(Arc::clone(sink));
            }
            if let Some(ref sinks) = self.partition_colocated_sinks {
                exec = exec.with_colocated_partition_sinks(sinks.clone());
            }
            let plan: Arc<dyn ExecutionPlan> = Arc::new(exec);
            return wrap_with_projection(plan, projection);
        }

        // Parquet is the only supported variation-lookup backend; a LookupProvider
        // is always constructed with a Parquet cache root by the annotation path.
        Err(DataFusionError::Plan(
            "LookupProvider::scan requires a Lance variation cache root".to_string(),
        ))
    }
}
