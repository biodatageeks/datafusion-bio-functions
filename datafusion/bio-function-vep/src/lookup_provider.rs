//! Lookup provider for `lookup_variants()` table function.
//!
//! Implements variant lookup using a self-contained left interval join
//! (`VariantLookupExec`) that builds VCF rows into per-chromosome COITrees
//! and streams cache rows as the probe side. When the cache is a KV store,
//! dispatches to `KvLookupExec` instead.

use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::array::{Array, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::Result;
use datafusion::datasource::{TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{Expr, SessionContext};

use crate::coordinate::CoordinateNormalizer;
use crate::schema_contract::validate_variation_schema;
use crate::variant_lookup_exec::{ColocatedSink, VariantLookupExec};

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

/// Table provider that implements variant lookup via self-contained left
/// interval join (`VariantLookupExec`).
///
/// VCF variants are collected into per-chromosome COITrees (build side).
/// Cache rows are streamed as the probe side with `match_allele()` as a
/// post-filter. Unmatched VCF rows appear with NULL cache columns.
///
/// # Default mode (equi-join)
///
/// Requires VCF and cache coordinates to match exactly after VEP
/// normalization (prefix/suffix stripping). SNVs and simple indels work well.
///
/// # Extended probes mode
///
/// When `extended_probes = true`, uses interval-overlap matching to handle
/// VEP-style coordinate encodings (insertion-style start > end,
/// prefix-trimmed deletion shifts, tandem-repeat right-shifting).
pub struct LookupProvider {
    session: Arc<SessionContext>,
    vcf_table: String,
    cache_table: String,
    /// Columns to select from the cache table.
    cache_columns: Vec<String>,
    /// Coordinate normalizer for handling different coordinate systems.
    coord_normalizer: CoordinateNormalizer,
    /// When true, use interval-overlap matching instead of exact coordinate match.
    extended_probes: bool,
    /// Maximum allowed `failed` flag value from the cache.
    allowed_failed: i64,
    /// Optional indexed reference FASTA used to materialize Ensembl genomic
    /// shift state for colocated existing-variant matching.
    reference_fasta_path: Option<String>,
    /// Output schema.
    schema: SchemaRef,
    /// Optional sink for co-located data collection during probe phase.
    colocated_sink: Option<ColocatedSink>,
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
        reference_fasta_path: Option<String>,
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
            coord_normalizer,
            extended_probes,
            allowed_failed,
            reference_fasta_path,
            schema,
            colocated_sink: None,
        })
    }

    /// Set the co-located data sink for piggybacked collection during probe.
    pub fn set_colocated_sink(&mut self, sink: ColocatedSink) {
        self.colocated_sink = Some(sink);
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
        // KV cache dispatch: if the cache table is a KvCacheTableProvider,
        // use the KvLookupExec instead.
        #[cfg(feature = "kv-cache")]
        {
            use crate::allele::allele_matches;
            use crate::kv_cache::KvCacheTableProvider;
            use crate::kv_cache::cache_exec::{KvLookupExec, KvMatchMode};

            let table_ref = tokio::task::block_in_place(|| {
                tokio::runtime::Handle::current()
                    .block_on(self.session.table_provider(&self.cache_table))
            });
            if let Ok(provider) = table_ref {
                if let Some(kv_provider) = provider.as_any().downcast_ref::<KvCacheTableProvider>()
                {
                    let store = kv_provider.store().clone();

                    let vcf_has_chr = has_chr_prefix(&self.session, &self.vcf_table).await?;

                    let vcf_df = self.session.table(&self.vcf_table).await?;
                    let vcf_plan = vcf_df.create_physical_plan().await?;

                    let mut exec = KvLookupExec::new(
                        vcf_plan,
                        store,
                        self.cache_columns.clone(),
                        KvMatchMode::Exact,
                        allele_matches as fn(&str, &str, &str) -> bool,
                        vcf_has_chr,
                        self.coord_normalizer.input_zero_based,
                        self.coord_normalizer.cache_zero_based,
                        self.extended_probes,
                        self.allowed_failed,
                    )?;
                    exec = exec.with_reference_fasta_path(self.reference_fasta_path.clone());
                    if let Some(ref sink) = self.colocated_sink {
                        exec = exec.with_colocated_sink(Arc::clone(sink));
                    }
                    let plan: Arc<dyn ExecutionPlan> = Arc::new(exec);
                    return wrap_with_projection(plan, projection);
                }
            }
        }

        // Parquet / MemTable path: use VariantLookupExec.
        let vcf_has_chr = has_chr_prefix(&self.session, &self.vcf_table).await?;

        let vcf_df = self.session.table(&self.vcf_table).await?;
        let vcf_plan = vcf_df.create_physical_plan().await?;

        // Project cache to only the columns needed: join keys + requested output columns.
        // This avoids reading all 78 parquet columns when only ~10 are needed.
        let cache_df = self.session.table(&self.cache_table).await?;
        let mut required_cols: Vec<&str> = vec!["chrom", "start", "end", "allele_string"];
        // Only request "failed" if the cache schema actually contains it.
        if cache_df
            .schema()
            .as_arrow()
            .field_with_name("failed")
            .is_ok()
        {
            required_cols.push("failed");
        }
        for col in &self.cache_columns {
            if !required_cols.contains(&col.as_str()) {
                required_cols.push(col.as_str());
            }
        }
        let cache_plan = cache_df
            .select_columns(&required_cols)?
            .create_physical_plan()
            .await?;

        let mut exec = VariantLookupExec::new(
            vcf_plan,
            cache_plan,
            self.cache_columns.clone(),
            vcf_has_chr,
            self.coord_normalizer.clone(),
            self.extended_probes,
            self.allowed_failed,
            self.reference_fasta_path.clone(),
            self.schema.clone(),
        );
        if let Some(ref sink) = self.colocated_sink {
            exec = exec.with_colocated_sink(Arc::clone(sink));
        }
        let plan: Arc<dyn ExecutionPlan> = Arc::new(exec);
        wrap_with_projection(plan, projection)
    }
}

#[cfg(test)]
mod tests {
    use crate::create_vep_session;
    #[cfg(feature = "kv-cache")]
    use crate::kv_cache::{
        KvCacheTableProvider, VepKvStore, position_entry::serialize_position_entry,
    };
    use datafusion::arrow::array::{
        Array, ArrayRef, Int64Array, RecordBatch, StringArray, StringViewArray,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use datafusion::physical_plan::displayable;
    use std::collections::HashMap;
    use std::sync::Arc;
    #[cfg(feature = "kv-cache")]
    use std::time::{SystemTime, UNIX_EPOCH};

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
            Field::new("failed", DataType::Int64, false),
        ]));
        // Cache coordinates must exactly match VCF after normalization.
        // VCF has chr1:(100,101) and chr1:(200,201), cache uses bare chrom names.
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![101, 201])),
                Arc::new(StringArray::from(vec!["rs123", "rs456"])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(StringArray::from(vec!["benign", "pathogenic"])),
                Arc::new(Int64Array::from(vec![0, 0])),
            ],
        )
        .unwrap();
        MemTable::try_new(schema, vec![vec![batch]]).unwrap()
    }

    #[cfg(feature = "kv-cache")]
    fn unique_temp_dir(prefix: &str) -> std::path::PathBuf {
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before UNIX_EPOCH")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("{prefix}-{}-{now}", std::process::id()));
        std::fs::create_dir_all(&path).expect("failed to create temp directory");
        path
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_generates_variant_lookup_exec() {
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
            plan_str.contains("VariantLookupExec"),
            "Expected VariantLookupExec in plan, got:\n{plan_str}"
        );
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_dispatches_to_kv_exec_for_fjall_cache_provider() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));

        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(StringArray::from(vec!["rs123"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["benign"])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .unwrap();

        // col_indices: end(2), variation_name(3), allele_string(4), clin_sig(5), failed(6)
        // (chrom=0 and start=1 are excluded from the entry; end is stored inside)
        let entry = serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-dispatch");

        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 100, &entry).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("var_cache", Arc::new(kv_provider))
            .unwrap();

        let df = ctx
            .sql(
                "SELECT * FROM lookup_variants('vcf_data', 'var_cache', 'variation_name,clin_sig')",
            )
            .await
            .unwrap();

        let plan = df.clone().create_physical_plan().await.unwrap();
        let plan_str = displayable(plan.as_ref()).indent(true).to_string();
        assert!(
            plan_str.contains("KvLookupExec"),
            "Expected KvLookupExec in plan for fjall provider, got:\n{plan_str}"
        );
        assert!(
            !plan_str.contains("IntervalJoinExec"),
            "Did not expect IntervalJoinExec for fjall provider, got:\n{plan_str}"
        );

        let batches = df.collect().await.unwrap();
        let mut annotations = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for value in string_values(col) {
                if let Some(v) = value {
                    annotations.push(v);
                }
            }
        }
        assert!(
            annotations.contains(&"rs123".to_string()),
            "Expected annotation from KV cache, got: {annotations:?}"
        );

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_matches_prefix_shifted_deletion_coordinates() {
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
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["TTA"])),
                Arc::new(StringArray::from(vec!["T"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        ctx.register_table("vcf_data", Arc::new(vcf)).unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["rs_shifted"])),
                Arc::new(StringArray::from(vec!["TA/-"])),
                Arc::new(StringArray::from(vec![Option::<&str>::None])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .unwrap();

        // col_indices: end(2), variation_name(3), allele_string(4), clin_sig(5), failed(6)
        let entry = serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-prefix-shift");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 101, &entry).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("var_cache", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql("SELECT * FROM lookup_variants('vcf_data', 'var_cache', 'variation_name,allele_string', 'exact', true)")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let mut names = Vec::new();
        let mut alleles = Vec::new();
        for batch in &batches {
            names.extend(string_values(
                batch.column_by_name("cache_variation_name").unwrap(),
            ));
            alleles.extend(string_values(
                batch.column_by_name("cache_allele_string").unwrap(),
            ));
        }

        assert!(
            names.iter().any(|v| v.as_deref() == Some("rs_shifted")),
            "Expected rs_shifted from KV lookup, got: {names:?}"
        );
        assert!(
            alleles.iter().any(|v| v.as_deref() == Some("TA/-")),
            "Expected TA/- allele from KV lookup, got: {alleles:?}"
        );

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_filters_failed_rows() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));

        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(Int64Array::from(vec![101, 101])),
                Arc::new(StringArray::from(vec!["rs_keep", "rs_failed"])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["benign", "pathogenic"])),
                Arc::new(Int64Array::from(vec![0, 1])),
            ],
        )
        .unwrap();

        // col_indices: end(2), variation_name(3), allele_string(4), clin_sig(5), failed(6)
        let entry = serialize_position_entry(&[0, 1], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-failed-filter");

        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 100, &entry).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("var_cache", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql(
                "SELECT * FROM lookup_variants('vcf_data', 'var_cache', 'variation_name,clin_sig')",
            )
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let mut names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for value in string_values(col) {
                if let Some(v) = value {
                    names.push(v);
                }
            }
        }

        assert_eq!(names, vec!["rs_keep".to_string()]);

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_matches_repeat_shifted_deletion_coordinates() {
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
                Arc::new(Int64Array::from(vec![55])),
                Arc::new(Int64Array::from(vec![70])),
                Arc::new(StringArray::from(vec!["GAAGAAGAAGAAGAA"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        ctx.register_table("vcf_data", Arc::new(vcf)).unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![66])),
                Arc::new(Int64Array::from(vec![79])),
                Arc::new(StringArray::from(vec!["rs_repeat_shift"])),
                Arc::new(StringArray::from(vec!["AAGAAGAAGAAGAA/-"])),
                Arc::new(StringArray::from(vec![Option::<&str>::None])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .unwrap();

        // col_indices: end(2), variation_name(3), allele_string(4), clin_sig(5), failed(6)
        let entry = serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-repeat-shift");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 66, &entry).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("var_cache", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql("SELECT * FROM lookup_variants('vcf_data', 'var_cache', 'variation_name,allele_string', 'exact', true)")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let mut names = Vec::new();
        let mut alleles = Vec::new();
        for batch in &batches {
            names.extend(string_values(
                batch.column_by_name("cache_variation_name").unwrap(),
            ));
            alleles.extend(string_values(
                batch.column_by_name("cache_allele_string").unwrap(),
            ));
        }

        assert!(
            names
                .iter()
                .any(|v| v.as_deref() == Some("rs_repeat_shift")),
            "Expected rs_repeat_shift from KV lookup, got: {names:?}"
        );
        assert!(
            alleles
                .iter()
                .any(|v| v.as_deref() == Some("AAGAAGAAGAAGAA/-")),
            "Expected AAGAAGAAGAAGAA/- allele from KV lookup, got: {alleles:?}"
        );

        let _ = std::fs::remove_dir_all(&cache_dir);
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

        // Verify the plan uses VariantLookupExec
        let plan = df.clone().create_physical_plan().await.unwrap();
        let plan_str = displayable(plan.as_ref()).indent(true).to_string();
        assert!(
            plan_str.contains("VariantLookupExec"),
            "Expected VariantLookupExec in plan, got:\n{plan_str}"
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
            let chrom_col = batch.column_by_name("chrom").unwrap();
            let clin_sig_col = batch.column_by_name("cache_clin_sig").unwrap();
            let chrom_vals = string_values(chrom_col);
            let clin_sig_vals = string_values(clin_sig_col);
            for i in 0..batch.num_rows() {
                chroms.push(chrom_vals[i].clone().unwrap_or_else(|| "NULL".to_string()));
                clin_sigs.push(
                    clin_sig_vals[i]
                        .clone()
                        .unwrap_or_else(|| "NULL".to_string()),
                );
            }
        }

        // The chr2 variant has no cache match — it should appear with NULL annotation.
        // Output keeps the original VCF chrom value ("chr2").
        assert!(
            chroms.contains(&"chr2".to_string()),
            "Expected 'chr2' (unmatched VCF row) in output, got chroms: {chroms:?}"
        );

        let chr2_idx = chroms.iter().position(|c| c == "chr2").unwrap();
        assert_eq!(
            clin_sigs[chr2_idx], "NULL",
            "Expected NULL clin_sig for unmatched chr2 row, got: {}",
            clin_sigs[chr2_idx]
        );
    }

    /// Regression test: equi-join ensures VCF rows only match cache entries
    /// at exactly the same (chrom, start, end) coordinates. Decoy entries at
    /// adjacent positions must NOT produce output rows.
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_no_false_matches_at_adjacent_positions() {
        let ctx = create_vep_session();

        // VCF: 1-based closed SNVs (start=end for point variants)
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
                Arc::new(Int64Array::from(vec![100, 200])),
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
            Field::new("failed", DataType::Int64, false),
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
                Arc::new(Int64Array::from(vec![0, 0, 0, 0])),
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
    async fn test_lookup_matches_pipe_joined_multi_alt_input() {
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
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G|T"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(StringArray::from(vec!["rs_g_match", "rs_c_miss"])),
                Arc::new(StringArray::from(vec!["A/G", "A/C"])),
                Arc::new(StringArray::from(vec!["benign", "benign"])),
                Arc::new(Int64Array::from(vec![0, 0])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_multi_alt_pipe", Arc::new(vcf))
            .unwrap();
        ctx.register_table("cache_multi_alt_pipe", Arc::new(cache))
            .unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_multi_alt_pipe', 'cache_multi_alt_pipe', 'variation_name')")
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 1, "Expected one row from ALT=G|T matching A/G");

        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        assert_eq!(var_names, vec!["rs_g_match"]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_matches_insertion_style_cache_coordinates() {
        let ctx = create_vep_session();

        // One-based closed insertion at position 100: REF=A ALT=AT.
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
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["AT"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();

        // Cache uses insertion-style coordinates: start=end+1.
        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![101, 102])),
                Arc::new(Int64Array::from(vec![100, 101])),
                Arc::new(StringArray::from(vec!["rs_ins_match", "rs_ins_miss"])),
                Arc::new(StringArray::from(vec!["-/T", "-/T"])),
                Arc::new(StringArray::from(vec!["pathogenic", "pathogenic"])),
                Arc::new(Int64Array::from(vec![0, 0])),
            ],
        )
        .unwrap();
        let cache = MemTable::try_new(cache_schema, vec![vec![cache_batch]]).unwrap();

        ctx.register_table("vcf_insertion", Arc::new(vcf)).unwrap();
        ctx.register_table("cache_insertion_style", Arc::new(cache))
            .unwrap();

        let df = ctx
            .sql("SELECT * FROM lookup_variants('vcf_insertion', 'cache_insertion_style', 'variation_name', 'exact', true)")
            .await
            .unwrap();
        let batches = df.collect().await.unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 1,
            "Expected one match against insertion-style cache coordinates"
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
        assert_eq!(var_names, vec!["rs_ins_match"]);
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
                Field::new("failed", DataType::Int64, false),
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
                Arc::new(Int64Array::from(vec![0])),
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
                Field::new("failed", DataType::Int64, false),
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
                Arc::new(Int64Array::from(vec![0, 0])),
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
                Field::new("failed", DataType::Int64, false),
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
                Arc::new(Int64Array::from(vec![0, 0])),
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
    /// variation_names each produce a separate output row.
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
                Arc::new(Int64Array::from(vec![100])),
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
            Field::new("failed", DataType::Int64, false),
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
                Arc::new(Int64Array::from(vec![0, 0, 0])),
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

    // -----------------------------------------------------------------------
    // KV (fjall) parity tests — mirrors of the Parquet-path tests above
    // -----------------------------------------------------------------------

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_left_join_preserves_unmatched_vcf_rows() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));

        // Two cache entries matching chr1:100 and chr1:200 — nothing for chr2:500.
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![101, 201])),
                Arc::new(StringArray::from(vec!["rs123", "rs456"])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(StringArray::from(vec!["benign", "pathogenic"])),
                Arc::new(Int64Array::from(vec![0, 0])),
            ],
        )
        .unwrap();

        // Store both alleles at their respective positions.
        let entry_100 = serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let entry_200 = serialize_position_entry(&[1], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-left-join");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 100, &entry_100).unwrap();
        store.put_position_entry("1", 200, &entry_200).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("var_cache", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql("SELECT * FROM lookup_variants('vcf_data', 'var_cache', 'clin_sig')")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 3,
            "Expected 3 rows (2 matched + 1 unmatched), got {total_rows}"
        );

        let mut chroms = Vec::new();
        let mut clin_sigs = Vec::new();
        for batch in &batches {
            let chrom_col = batch.column_by_name("chrom").unwrap();
            let clin_sig_col = batch.column_by_name("cache_clin_sig").unwrap();
            let chrom_vals = string_values(chrom_col);
            let clin_sig_vals = string_values(clin_sig_col);
            for i in 0..batch.num_rows() {
                chroms.push(chrom_vals[i].clone().unwrap_or_else(|| "NULL".to_string()));
                clin_sigs.push(
                    clin_sig_vals[i]
                        .clone()
                        .unwrap_or_else(|| "NULL".to_string()),
                );
            }
        }

        assert!(
            chroms.contains(&"chr2".to_string()),
            "Expected 'chr2' (unmatched VCF row) in output, got chroms: {chroms:?}"
        );
        let chr2_idx = chroms.iter().position(|c| c == "chr2").unwrap();
        assert_eq!(
            clin_sigs[chr2_idx], "NULL",
            "Expected NULL clin_sig for unmatched chr2 row, got: {}",
            clin_sigs[chr2_idx]
        );

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_no_false_matches_at_adjacent_positions() {
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
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(StringArray::from(vec!["A", "C"])),
                Arc::new(StringArray::from(vec!["G", "T"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        ctx.register_table("vcf_kv_adj", Arc::new(vcf)).unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));

        // TRUE matches at 100 and 200; DECOYS at 101 and 201.
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 101, 200, 201])),
                Arc::new(Int64Array::from(vec![100, 101, 200, 201])),
                Arc::new(StringArray::from(vec!["rs100", "rs101", "rs200", "rs201"])),
                Arc::new(StringArray::from(vec!["A/G", "T/G", "C/T", "G/T"])),
                Arc::new(StringArray::from(vec![
                    "benign",
                    "decoy_101",
                    "pathogenic",
                    "decoy_201",
                ])),
                Arc::new(Int64Array::from(vec![0, 0, 0, 0])),
            ],
        )
        .unwrap();

        let cache_dir = unique_temp_dir("vep-kv-adj-pos");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();

        // Each position gets its own entry.
        let entry_100 = serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let entry_101 = serialize_position_entry(&[1], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let entry_200 = serialize_position_entry(&[2], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let entry_201 = serialize_position_entry(&[3], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        store.put_position_entry("1", 100, &entry_100).unwrap();
        store.put_position_entry("1", 101, &entry_101).unwrap();
        store.put_position_entry("1", 200, &entry_200).unwrap();
        store.put_position_entry("1", 201, &entry_201).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("cache_kv_adj", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql("SELECT * FROM lookup_variants('vcf_kv_adj', 'cache_kv_adj', 'variation_name,clin_sig')")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 2,
            "Expected exactly 2 rows (one match per VCF variant, no decoys), got {total_rows}"
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
        assert!(var_names.contains(&"rs100".to_string()));
        assert!(var_names.contains(&"rs200".to_string()));
        assert!(!var_names.contains(&"rs101".to_string()));
        assert!(!var_names.contains(&"rs201".to_string()));

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_matches_pipe_joined_multi_alt_input() {
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
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G|T"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        ctx.register_table("vcf_kv_multi", Arc::new(vcf)).unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(Int64Array::from(vec![100, 100])),
                Arc::new(StringArray::from(vec!["rs_g_match", "rs_c_miss"])),
                Arc::new(StringArray::from(vec!["A/G", "A/C"])),
                Arc::new(StringArray::from(vec!["benign", "benign"])),
                Arc::new(Int64Array::from(vec![0, 0])),
            ],
        )
        .unwrap();

        // Both alleles at same position → single position entry.
        let entry = serialize_position_entry(&[0, 1], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-multi-alt");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 100, &entry).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("cache_kv_multi", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql(
                "SELECT * FROM lookup_variants('vcf_kv_multi', 'cache_kv_multi', 'variation_name')",
            )
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 1, "Expected one row from ALT=G|T matching A/G");

        let mut var_names = Vec::new();
        for batch in &batches {
            let col = batch.column_by_name("cache_variation_name").unwrap();
            for val in string_values(col) {
                if let Some(v) = val {
                    var_names.push(v);
                }
            }
        }
        assert_eq!(var_names, vec!["rs_g_match"]);

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_matches_insertion_style_cache_coordinates() {
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
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["AT"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        ctx.register_table("vcf_kv_ins", Arc::new(vcf)).unwrap();

        // Insertion-style cache: start=101, end=100 (start > end).
        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![101, 102])),
                Arc::new(Int64Array::from(vec![100, 101])),
                Arc::new(StringArray::from(vec!["rs_ins_match", "rs_ins_miss"])),
                Arc::new(StringArray::from(vec!["-/T", "-/T"])),
                Arc::new(StringArray::from(vec!["pathogenic", "pathogenic"])),
                Arc::new(Int64Array::from(vec![0, 0])),
            ],
        )
        .unwrap();

        let cache_dir = unique_temp_dir("vep-kv-ins-style");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        // Store at the cache start positions.
        let entry_101 = serialize_position_entry(&[0], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let entry_102 = serialize_position_entry(&[1], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        store.put_position_entry("1", 101, &entry_101).unwrap();
        store.put_position_entry("1", 102, &entry_102).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("cache_kv_ins", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql("SELECT * FROM lookup_variants('vcf_kv_ins', 'cache_kv_ins', 'variation_name', 'exact', true)")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 1,
            "Expected one match against insertion-style cache coordinates"
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
        assert_eq!(var_names, vec!["rs_ins_match"]);

        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_lookup_kv_cache_duplicates_produce_separate_rows() {
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
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .unwrap();
        let vcf = MemTable::try_new(vcf_schema, vec![vec![vcf_batch]]).unwrap();
        ctx.register_table("vcf_kv_dup", Arc::new(vcf)).unwrap();

        let cache_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        // Three alleles at the same position, all A/G with different variation_names.
        let cache_batch = RecordBatch::try_new(
            cache_schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(StringArray::from(vec!["rs100", "COSM123", "rs100_dup"])),
                Arc::new(StringArray::from(vec!["A/G", "A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["benign", "benign", "benign"])),
                Arc::new(Int64Array::from(vec![0, 0, 0])),
            ],
        )
        .unwrap();

        // All three rows go into a single position entry.
        let entry =
            serialize_position_entry(&[0, 1, 2], &cache_batch, &[2, 3, 4, 5, 6], 4).unwrap();
        let cache_dir = unique_temp_dir("vep-kv-dup");
        let store = VepKvStore::create(&cache_dir, cache_schema).unwrap();
        store.put_position_entry("1", 100, &entry).unwrap();
        store.persist().unwrap();
        drop(store);

        let kv_provider = KvCacheTableProvider::open(&cache_dir).unwrap();
        ctx.register_table("cache_kv_dup", Arc::new(kv_provider))
            .unwrap();

        let batches = ctx
            .sql("SELECT * FROM lookup_variants('vcf_kv_dup', 'cache_kv_dup', 'variation_name,clin_sig')")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 3,
            "Expected 3 rows (one per cache duplicate), got {total_rows}"
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
        var_names.sort();
        assert_eq!(
            var_names,
            vec!["COSM123", "rs100", "rs100_dup"],
            "Expected all three variation_names in output"
        );

        let _ = std::fs::remove_dir_all(&cache_dir);
    }
}
