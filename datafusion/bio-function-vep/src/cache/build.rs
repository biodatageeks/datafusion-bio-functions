use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use datafusion::arrow::array::{
    Array, ArrayRef, BinaryArray, BooleanArray, Float64Array, Int8Array, Int64Array, StringArray,
    UInt32Array, UInt64Array,
};
use datafusion::arrow::compute::{
    cast, concat_batches, filter_record_batch, sort_to_indices, take,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{ExecutionPlan, SendableRecordBatchStream};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_ensembl_cache::{
    CacheSourceType as BioFormatsCacheSourceType, EnsemblCacheOptions, EnsemblCacheTableProvider,
    EnsemblEntityKind, VEP_CACHE_REGION_SIZE_BP, build_export_query, translation_core_schema,
    translation_sift_schema,
};
use futures::StreamExt;
use log::info;

use crate::cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest, canonical_chrom_label,
};
use crate::cache::schema::{
    VARIATION_FORBIDDEN_COLUMNS, VARIATION_REQUIRED_COLUMNS, variation_projected_schema,
    with_cache_source_metadata,
};
use crate::cache_builder::EntityStats;
use crate::cache_common::{
    FrequencyFields, PositionFrequency, max_global_af, select_warm_positions,
};
use crate::{
    annotate_provider::{read_compact_predictions, string_at},
    transcript_consequence::CachedPredictions,
};

const DEFAULT_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];
/// Warm/cold tier selection for the single-path variation layout. A genomic
/// `start` is "warm" (tier 0) when its max global allele frequency across
/// `minor_allele_freq`, `AF`, `gnomADg`, `gnomADe` is at least this threshold;
/// the radius extends warm membership to neighboring positions that exist in the
/// data. These must match `warm_cache::split` so the Parquet tier agrees with any
/// other tiered consumer.
const WARM_AF_THRESHOLD: f64 = 0.01;
const WARM_POSITION_RADIUS: i64 = 1;

#[derive(Debug, Clone)]
pub struct CacheBuildOptions {
    pub cache_root: String,
    pub output_dir: String,
    pub partitions: usize,
    pub cache_source_type: BioFormatsCacheSourceType,
    pub overwrite: bool,
    /// Optional allow-list of chromosomes to build (e.g. `["chr1"]`). Matched
    /// after stripping any `chr` prefix, so `"1"` and `"chr1"` are equivalent.
    /// `None` builds every chromosome discovered in the source cache. Used to
    /// scope rebuilds (e.g. a single-chromosome profiling/parity build).
    pub chrom_filter: Option<Vec<String>>,
}

/// Build one chromosome's
/// variation shard as a single no-dictionary, page-indexed `.parquet` file under
/// `variation/<chrom>.parquet`.
///
/// Reuses the exact tiering pipeline — the warm-start pre-pass
/// (`collect_warm_starts`) and `transform_variation_tier_batch` — then applies
/// the Parquet encoding (Boolean flags, 2-array AF, `variation_name` dedup, via
/// [`crate::parquet_cache::write::encode_variation_batch`]) and streams the result
/// to the shard. Rows are written warm tier (0) first then cold (1), each
/// `start`-ordered by the source `ORDER BY start` plan — the two sorted runs the
/// read-side `PageDir` binary search relies on.
pub async fn build_parquet_variation_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    use crate::parquet_cache::write::{
        VariationParquetShardWriter, encode_variation_batch, variation_output_schema,
    };

    let entity_dir = PathBuf::from(&options.output_dir).join("variation");
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let query = build_export_query(
        EnsemblEntityKind::Variation,
        "var",
        Some(source_chrom),
        None,
    );
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    if shard_path.exists() {
        if !options.overwrite {
            return Err(DataFusionError::Execution(format!(
                "parquet variation shard '{}' already exists and overwrite=false",
                shard_path.display()
            )));
        }
        std::fs::remove_file(&shard_path).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to remove existing parquet shard '{}': {err}",
                shard_path.display()
            ))
        })?;
    }

    let source_type = options.cache_source_type.as_str().to_string();
    let warm_starts = Arc::new(collect_warm_starts(options, source_chrom).await?);
    info!(
        "Parquet variation {} warm positions={} (af>={WARM_AF_THRESHOLD}, +/-{WARM_POSITION_RADIUS})",
        shard_path.display(),
        warm_starts.len(),
    );

    // Warm tier (0) first, then cold (1); a chromosome with no warm positions
    // writes a single cold tier. Mirrors the Parquet tiered writer's passes.
    let passes: &[i8] = if warm_starts.is_empty() {
        &[1]
    } else {
        &[0, 1]
    };
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Variation, "var")?;

    // Derive the output schema once from the source plan schema and reuse it for
    // every written batch so the shard carries one stable Arrow schema.
    let schema_plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let out_schema = variation_output_schema(schema_plan.schema().as_ref(), &source_type)?;

    let mut writer = VariationParquetShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
    for &keep_tier in passes {
        let plan = ctx.sql(&query).await?.create_physical_plan().await?;
        let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;
        while let Some(batch) = stream.next().await {
            let tiered =
                transform_variation_tier_batch(batch?, &source_type, &warm_starts, keep_tier)?;
            if tiered.num_rows() == 0 {
                continue;
            }
            let encoded = encode_variation_batch(&tiered, &out_schema)?;
            writer.write(&encoded)?;
        }
    }
    let rows = writer.finish()?;
    if rows == 0 {
        // The writer is created eagerly (before the stream), so a contig with no
        // variants — e.g. a GL*/KI* scaffold present in the source schema but
        // carrying no variant rows — leaves a schema-only shard. Drop it so the
        // on-disk layout matches the manifest (which omits empty contigs), like
        // the lazily-written context/translation shards.
        let _ = std::fs::remove_file(&shard_path);
    }
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
}

/// Run a cheap pre-pass over the source variation table for one chromosome to
/// classify each genomic `start` as warm (common) or cold (rare) using the same
/// allele-frequency selector as the Parquet warm tier.
async fn collect_warm_starts(
    options: &CacheBuildOptions,
    source_chrom: &str,
) -> Result<BTreeSet<i64>> {
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Variation, "var")?;
    // Mixed-case columns must be quoted; DataFusion lowercases bare identifiers.
    let sql = format!(
        "SELECT start, minor_allele_freq, \"AF\", \"gnomADg\", \"gnomADe\" FROM var WHERE chrom = '{}'",
        sql_escape_literal(source_chrom)
    );
    let df = ctx.sql(&sql).await?;
    let mut stream = df.execute_stream().await?;
    let mut rows = Vec::new();
    while let Some(batch) = stream.next().await {
        collect_position_frequencies(&batch?, &mut rows)?;
    }
    Ok(select_warm_positions(
        rows,
        WARM_AF_THRESHOLD,
        WARM_POSITION_RADIUS,
    ))
}

fn collect_position_frequencies(
    batch: &RecordBatch,
    out: &mut Vec<PositionFrequency>,
) -> Result<()> {
    if batch.num_rows() == 0 {
        return Ok(());
    }
    let schema = batch.schema();
    let start_idx = schema.index_of("start")?;
    let start = cast(batch.column(start_idx).as_ref(), &DataType::Int64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation start to Int64: {err}"))
    })?;
    let start = start
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| DataFusionError::Execution("variation start must be Int64".to_string()))?;

    let maf = column_as_f64(batch, "minor_allele_freq")?;
    let af = column_as_str(batch, "AF");
    let gnomadg = column_as_str(batch, "gnomADg");
    let gnomade = column_as_str(batch, "gnomADe");

    out.reserve(batch.num_rows());
    for row in 0..batch.num_rows() {
        if start.is_null(row) {
            continue;
        }
        let af_s = af.as_ref().and_then(|arr| string_at(arr.as_ref(), row));
        let gnomadg_s = gnomadg
            .as_ref()
            .and_then(|arr| string_at(arr.as_ref(), row));
        let gnomade_s = gnomade
            .as_ref()
            .and_then(|arr| string_at(arr.as_ref(), row));
        let fields = FrequencyFields {
            minor_allele_freq: maf
                .as_ref()
                .and_then(|arr| (!arr.is_null(row)).then(|| arr.value(row))),
            af: af_s.as_deref(),
            gnomadg: gnomadg_s.as_deref(),
            gnomade: gnomade_s.as_deref(),
        };
        out.push(PositionFrequency {
            start: start.value(row),
            max_af: max_global_af(&fields),
        });
    }
    Ok(())
}

fn column_as_f64(batch: &RecordBatch, name: &str) -> Result<Option<Float64Array>> {
    let Ok(index) = batch.schema().index_of(name) else {
        return Ok(None);
    };
    let array = cast(batch.column(index).as_ref(), &DataType::Float64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation {name} to Float64: {err}"))
    })?;
    let array = array
        .as_any()
        .downcast_ref::<Float64Array>()
        .ok_or_else(|| DataFusionError::Execution(format!("variation {name} must be Float64")))?
        .clone();
    Ok(Some(array))
}

fn column_as_str(batch: &RecordBatch, name: &str) -> Option<ArrayRef> {
    batch
        .schema()
        .index_of(name)
        .ok()
        .map(|index| batch.column(index).clone())
}

/// Build, for one chromosome,
/// the position-sliced SIFT/PolyPhen shard as a single no-dictionary,
/// page-indexed `.parquet` file under `translation_sift/<chrom>.parquet`,
/// physically sorted by the u64 `key` so the read-side `PageDir` sees one monotone
/// run.
///
/// Reuses the exact row pipeline — the shared `transcript_uid` map and
/// `transform_translation_sift_position_batch` — then globally sorts the
/// accumulated `(key, sift, poly)` rows by `key` (the source dedup plan groups by
/// transcript, not by the derived key) and writes via
/// [`crate::parquet_cache::sift::write_sift_parquet`]. The schema metadata
/// (`SIFT_BLOB_VERSION`) is preserved so the reader decodes blobs unchanged.
pub async fn build_parquet_translation_sift_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let entity_dir = PathBuf::from(&options.output_dir).join("translation_sift");
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let query = build_translation_dedup_query_with_where_clause(
        "tl",
        &format!(" WHERE chrom = '{}'", sql_escape_literal(source_chrom)),
    );
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    if shard_path.exists() {
        if !options.overwrite {
            return Err(DataFusionError::Execution(format!(
                "parquet translation_sift shard '{}' already exists and overwrite=false",
                shard_path.display()
            )));
        }
        std::fs::remove_file(&shard_path).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to remove existing parquet shard '{}': {err}",
                shard_path.display()
            ))
        })?;
    }

    let source_type = options.cache_source_type.as_str().to_string();
    // Same shared uid map the transcript build uses, so keys agree.
    let uid_map = Arc::new(load_transcript_uid_map(options, source_chrom).await?);
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Translation, "tl")?;
    let plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

    let mut batches: Vec<RecordBatch> = Vec::new();
    while let Some(batch) = stream.next().await {
        let batch = drop_row_number_batch(batch?)?;
        let out = transform_translation_sift_position_batch(batch, &source_type, &uid_map)?;
        if out.num_rows() > 0 {
            batches.push(out);
        }
    }
    if batches.is_empty() {
        return Ok(ChromDatasetEntry::new(manifest_chrom, file_name, 0));
    }

    // Global sort by the derived u64 `key` (unique by construction) → one monotone
    // PageDir run, matching the Parquet BTree key order.
    let schema = batches[0].schema();
    let combined = concat_batches(&schema, &batches).map_err(|err| {
        DataFusionError::Execution(format!("failed to concat sift batches: {err}"))
    })?;
    let sorted = sort_record_batch_by_key(&combined)?;
    let rows = crate::parquet_cache::sift::write_sift_parquet(&shard_path, &[sorted])?;
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
}

/// Remove an existing Parquet shard file, honoring `overwrite`. Shared by the
/// Parquet build drivers.
fn remove_existing_parquet_shard(path: &Path, overwrite: bool) -> Result<()> {
    if path.exists() {
        if !overwrite {
            return Err(DataFusionError::Execution(format!(
                "parquet shard '{}' already exists and overwrite=false",
                path.display()
            )));
        }
        std::fs::remove_file(path).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to remove existing parquet shard '{}': {err}",
                path.display()
            ))
        })?;
    }
    Ok(())
}

/// The provider output schema for an entity (used to build the transcript export
/// query, which selects an explicit field list keyed on the provider schema).
fn provider_output_schema(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<SchemaRef> {
    use datafusion::catalog::TableProvider;
    let mut provider_options = EnsemblCacheOptions::new(&options.cache_root)
        .with_cache_source_type(options.cache_source_type);
    provider_options.target_partitions = Some(options.partitions);
    Ok(EnsemblCacheTableProvider::for_entity(kind, provider_options)?.schema())
}

/// Stream a source query to a dictionary-enabled context Parquet shard, applying
/// `transform` to each batch. Preserves source (query) order — no explicit sort —
/// so the shard matches the Parquet dataset row order. Returns the row count
/// (0 → no file written).
async fn write_query_stream_to_parquet<F>(
    ctx: &SessionContext,
    query: &str,
    path: &Path,
    transform: F,
) -> Result<usize>
where
    F: Fn(RecordBatch) -> Result<RecordBatch>,
{
    use crate::parquet_cache::scan::ContextParquetShardWriter;
    let plan = ctx.sql(query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;
    let mut writer: Option<ContextParquetShardWriter> = None;
    let mut rows = 0usize;
    while let Some(batch) = stream.next().await {
        let out = transform(batch?)?;
        if out.num_rows() == 0 {
            continue;
        }
        if writer.is_none() {
            writer = Some(ContextParquetShardWriter::create(path, out.schema())?);
        }
        writer.as_mut().unwrap().write(&out)?;
        rows += out.num_rows();
    }
    if let Some(writer) = writer {
        writer.finish()?;
    }
    Ok(rows)
}

/// Build, for one chromosome, one
/// scan entity (`transcript`/`exon`/`regulatory`/`motif`): dictionary-enabled
/// `parquet.<entity>/<chrom>.parquet`, written verbatim from the same export
/// query + per-batch transform the Parquet path uses (transcript also attaches the
/// shared `transcript_uid`). Read back by a full projected scan.
pub async fn build_parquet_context_entity_chrom(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let entity = entity_output_name(kind);
    let entity_dir = PathBuf::from(&options.output_dir).join(entity);
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let table_name = entity_table_name(kind);
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    remove_existing_parquet_shard(&shard_path, options.overwrite)?;

    // Only transcript's export query needs the provider schema (explicit select list).
    let provider_schema = if kind == EnsemblEntityKind::Transcript {
        Some(provider_output_schema(options, kind)?)
    } else {
        None
    };
    let query = build_export_query(
        kind,
        table_name,
        Some(source_chrom),
        provider_schema.as_deref(),
    );
    let source_type = options.cache_source_type.as_str().to_string();
    let ctx = make_ctx_and_register(options, kind, table_name)?;

    let rows = if kind == EnsemblEntityKind::Transcript {
        let uid_map = Arc::new(load_transcript_uid_map(options, source_chrom).await?);
        write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
            attach_transcript_uid_batch(batch, &source_type, &uid_map)
        })
        .await?
    } else {
        write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
            attach_schema_metadata_to_batch(batch, &source_type)
        })
        .await?
    };
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
}

/// Build `translation_core`,
/// one chromosome: dictionary-enabled `translation_core/<chrom>.parquet`,
/// projected to `translation_core_schema` (big sequences + `list<struct>`
/// `protein_features`). Read back by a full projected scan.
pub async fn build_parquet_translation_core_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let entity_dir = PathBuf::from(&options.output_dir).join("translation_core");
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let query = build_translation_dedup_query_with_where_clause(
        "tl",
        &format!(" WHERE chrom = '{}'", sql_escape_literal(source_chrom)),
    );
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    remove_existing_parquet_shard(&shard_path, options.overwrite)?;

    let target_schema = translation_core_schema(false, options.cache_source_type);
    let source_type = options.cache_source_type.as_str().to_string();
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Translation, "tl")?;
    let rows = write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
        let batch = drop_row_number_batch(batch)?;
        let batch = project_batch_to_schema(batch, Arc::clone(&target_schema))?;
        attach_schema_metadata_to_batch(batch, &source_type)
    })
    .await?;
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
}

/// Output directory for a Parquet entity's shard set: `parquet.<entity>/`.
fn parquet_entity_dir(options: &CacheBuildOptions, entity: &str) -> Result<PathBuf> {
    let path = Path::new(&options.output_dir).join(entity);
    std::fs::create_dir_all(&path).map_err(|err| {
        DataFusionError::Execution(format!("failed to create {}: {err}", path.display()))
    })?;
    Ok(path)
}

/// Skip a Parquet entity build when a manifest and at least one `.parquet` shard
/// already exist and `overwrite` is false. Parquet analogue of
fn should_skip_parquet_entity(entity_dir: &Path, overwrite: bool) -> bool {
    !overwrite
        && entity_dir.join(CHROM_MANIFEST_FILE).exists()
        && std::fs::read_dir(entity_dir)
            .ok()
            .into_iter()
            .flat_map(|entries| entries.flatten())
            .any(|entry| entry.file_name().to_string_lossy().ends_with(".parquet"))
}

/// Build one entity as Parquet shards under `parquet.<entity>/`, writing the
/// per-chrom `chrom_manifest.json`.
/// this is the entry point the [`crate::cache_builder::CacheBuilder`] dispatches
/// to.
pub async fn build_parquet_entity(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<Vec<EntityStats>> {
    match kind {
        EnsemblEntityKind::Variation => build_parquet_variation(options).await,
        EnsemblEntityKind::Translation => build_parquet_translation(options).await,
        EnsemblEntityKind::Transcript
        | EnsemblEntityKind::Exon
        | EnsemblEntityKind::RegulatoryFeature
        | EnsemblEntityKind::MotifFeature => build_parquet_context(options, kind).await,
    }
}

async fn build_parquet_variation(options: &CacheBuildOptions) -> Result<Vec<EntityStats>> {
    let entity_dir = parquet_entity_dir(options, "variation")?;
    if should_skip_parquet_entity(&entity_dir, options.overwrite) {
        info!("variation already exists, skipping (use overwrite to rebuild)");
        return Ok(vec![EntityStats {
            entity: "variation".to_string(),
            parquet_files: vec![],
        }]);
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Variation, "var").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let entry = build_parquet_variation_chrom(options, &chrom).await?;
        if entry.rows > 0 {
            let shard_path = entity_dir.join(&entry.dataset);
            files.push((shard_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: "variation".to_string(),
        parquet_files: files,
    }])
}

async fn build_parquet_context(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<Vec<EntityStats>> {
    let entity = entity_output_name(kind);
    let entity_dir = parquet_entity_dir(options, entity)?;
    if should_skip_parquet_entity(&entity_dir, options.overwrite) {
        info!("{entity} already exists, skipping (use overwrite to rebuild)");
        return Ok(vec![EntityStats {
            entity: entity.to_string(),
            parquet_files: vec![],
        }]);
    }

    let table_name = entity_table_name(kind);
    let chroms = discover_chroms(options, kind, table_name).await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let entry = build_parquet_context_entity_chrom(options, kind, &chrom).await?;
        if entry.rows > 0 {
            let shard_path = entity_dir.join(&entry.dataset);
            files.push((shard_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: entity.to_string(),
        parquet_files: files,
    }])
}

async fn build_parquet_translation(options: &CacheBuildOptions) -> Result<Vec<EntityStats>> {
    let core = build_parquet_translation_core(options).await?;
    let sift = build_parquet_translation_sift(options).await?;
    Ok(vec![core, sift])
}

async fn build_parquet_translation_core(options: &CacheBuildOptions) -> Result<EntityStats> {
    let entity_dir = parquet_entity_dir(options, "translation_core")?;
    if should_skip_parquet_entity(&entity_dir, options.overwrite) {
        info!("translation_core already exists, skipping (use overwrite to rebuild)");
        return Ok(EntityStats {
            entity: "translation_core".to_string(),
            parquet_files: vec![],
        });
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Translation, "tl").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let entry = build_parquet_translation_core_chrom(options, &chrom).await?;
        if entry.rows > 0 {
            let shard_path = entity_dir.join(&entry.dataset);
            files.push((shard_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(EntityStats {
        entity: "translation_core".to_string(),
        parquet_files: files,
    })
}

async fn build_parquet_translation_sift(options: &CacheBuildOptions) -> Result<EntityStats> {
    let entity_dir = parquet_entity_dir(options, "translation_sift")?;
    if should_skip_parquet_entity(&entity_dir, options.overwrite) {
        info!("translation_sift already exists, skipping (use overwrite to rebuild)");
        return Ok(EntityStats {
            entity: "translation_sift".to_string(),
            parquet_files: vec![],
        });
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Translation, "tl").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let entry = build_parquet_translation_sift_chrom(options, &chrom).await?;
        if entry.rows > 0 {
            let shard_path = entity_dir.join(&entry.dataset);
            files.push((shard_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(EntityStats {
        entity: "translation_sift".to_string(),
        parquet_files: files,
    })
}

/// Sort a translation_sift batch ascending by its `key` column.
fn sort_record_batch_by_key(batch: &RecordBatch) -> Result<RecordBatch> {
    let key_idx = batch.schema().index_of("key")?;
    let indices = sort_to_indices(batch.column(key_idx), None, None)
        .map_err(|err| DataFusionError::Execution(format!("failed to sort sift key: {err}")))?;
    let columns = batch
        .columns()
        .iter()
        .map(|c| take(c, &indices, None))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|err| DataFusionError::Execution(format!("failed to take sorted sift: {err}")))?;
    RecordBatch::try_new(batch.schema(), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to rebuild sorted sift batch: {err}"))
    })
}

async fn discover_chroms(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<Vec<String>> {
    discover_chroms_and_schema(options, kind, table_name)
        .await
        .map(|(chroms, _)| chroms)
}

async fn discover_chroms_and_schema(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<(Vec<String>, SchemaRef)> {
    let ctx = make_ctx_and_register(options, kind, table_name)?;
    let table = ctx.table(table_name).await?;
    let schema = Arc::clone(table.schema().inner());
    let chroms = filter_chroms(chroms_from_schema(&schema), options.chrom_filter.as_deref());
    Ok((chroms, schema))
}

/// Restrict discovered chromosomes to the build's allow-list, comparing by the
/// bare label (any `chr` prefix stripped) so `"1"` and `"chr1"` match. `None`
/// passes every chromosome through unchanged.
fn filter_chroms(chroms: Vec<String>, allow: Option<&[String]>) -> Vec<String> {
    let Some(allow) = allow else {
        return chroms;
    };
    let bare = |label: &str| label.strip_prefix("chr").unwrap_or(label).to_string();
    let allowed: std::collections::HashSet<String> = allow.iter().map(|c| bare(c)).collect();
    chroms
        .into_iter()
        .filter(|chrom| allowed.contains(&bare(chrom)))
        .collect()
}

fn make_ctx_and_register(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<SessionContext> {
    let config = SessionConfig::new().with_target_partitions(options.partitions);
    let ctx = SessionContext::new_with_config(config);
    let mut provider_options = EnsemblCacheOptions::new(&options.cache_root)
        .with_cache_source_type(options.cache_source_type);
    provider_options.target_partitions = Some(options.partitions);
    let provider = EnsemblCacheTableProvider::for_entity(kind, provider_options)?;
    ctx.register_table(table_name, provider)?;
    Ok(ctx)
}

/// Project a source variation batch to the Parquet output schema, derive the
/// `tier` column from `warm_starts`, and keep only rows whose tier matches
/// `keep_tier`. The returned batch may be empty (all rows belonged to the other
/// tier); the streaming writer drops empty batches.
fn transform_variation_tier_batch(
    batch: RecordBatch,
    source_type: &str,
    warm_starts: &BTreeSet<i64>,
    keep_tier: i8,
) -> Result<RecordBatch> {
    let schema = batch.schema();
    let target_schema = Arc::new(variation_projected_schema(schema.as_ref(), source_type)?);
    let mut columns = Vec::new();
    let forbidden = VARIATION_FORBIDDEN_COLUMNS
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    for name in VARIATION_REQUIRED_COLUMNS {
        if forbidden.contains(name) {
            continue;
        }
        let (index, _) = schema.column_with_name(name).ok_or_else(|| {
            DataFusionError::Execution(format!("variation source batch missing column {name}"))
        })?;
        if *name == "start" || *name == "end" {
            columns.push(
                cast(batch.column(index).as_ref(), &DataType::UInt32).map_err(|err| {
                    DataFusionError::Execution(format!(
                        "failed to cast variation {name} to UInt32: {err}"
                    ))
                })?,
            );
        } else {
            columns.push(batch.column(index).clone());
        }
    }

    // Derive tier from the source `start` and build the keep mask in one pass.
    let (start_index, _) = schema.column_with_name("start").ok_or_else(|| {
        DataFusionError::Execution("variation source batch missing column start".to_string())
    })?;
    let start = cast(batch.column(start_index).as_ref(), &DataType::Int64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation start to Int64: {err}"))
    })?;
    let start = start
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| DataFusionError::Execution("variation start must be Int64".to_string()))?;

    let mut tier_values = Vec::with_capacity(batch.num_rows());
    let mut keep = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        let tier = if !start.is_null(row) && warm_starts.contains(&start.value(row)) {
            0i8
        } else {
            1i8
        };
        tier_values.push(tier);
        keep.push(tier == keep_tier);
    }
    columns.push(Arc::new(Int8Array::from(tier_values)) as ArrayRef);

    let full = RecordBatch::try_new(target_schema, columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to build Lance variation batch: {err}"))
    })?;
    let mask = BooleanArray::from(keep);
    filter_record_batch(&full, &mask).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to filter Lance variation tier batch: {err}"
        ))
    })
}

fn attach_schema_metadata_to_batch(batch: RecordBatch, source_type: &str) -> Result<RecordBatch> {
    let schema = with_cache_source_metadata(batch.schema().as_ref(), source_type);
    RecordBatch::try_new(Arc::new(schema), batch.columns().to_vec()).map_err(|err| {
        DataFusionError::Execution(format!("failed to attach Lance schema metadata: {err}"))
    })
}

/// Assign a dense `transcript_uid` (0..N) to each id, by the order given. Callers
/// pass a sorted-unique id list, so the uid is a stable function of the id set —
/// the transcript build and the sift build derive identical uids from the same
/// transcript source without depending on each other's artifacts.
fn assign_transcript_uids(sorted_unique_ids: &[String]) -> HashMap<String, u32> {
    sorted_unique_ids
        .iter()
        .enumerate()
        .map(|(index, id)| (id.clone(), index as u32))
        .collect()
}

/// Build the `transcript_id -> transcript_uid` map for one chromosome by ranking
/// the transcript source's sorted-unique `stable_id`s. Shared by the transcript
/// build (writes the `transcript_uid` column) and the sift build (keys rows).
async fn load_transcript_uid_map(
    options: &CacheBuildOptions,
    source_chrom: &str,
) -> Result<HashMap<String, u32>> {
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Transcript, "tx")?;
    let sql = format!(
        "SELECT DISTINCT stable_id FROM tx WHERE chrom = '{}' ORDER BY stable_id",
        sql_escape_literal(source_chrom)
    );
    let df = ctx.sql(&sql).await?;
    let mut stream = df.execute_stream().await?;
    let mut ids = Vec::new();
    while let Some(batch) = stream.next().await {
        let batch = batch?;
        let idx = batch.schema().index_of("stable_id")?;
        let column = batch.column(idx);
        for row in 0..batch.num_rows() {
            if let Some(id) = string_at(column.as_ref(), row) {
                ids.push(id);
            }
        }
    }
    // DISTINCT + ORDER BY already yields sorted-unique ids; normalize defensively
    // so uid assignment is deterministic regardless of plan ordering.
    ids.sort_unstable();
    ids.dedup();
    Ok(assign_transcript_uids(&ids))
}

/// Append a non-null `transcript_uid: UInt32` column to a transcript batch,
/// looked up from `uid_map` by the row's `stable_id`. Errors if any id is absent
/// (the map is built from the same transcript source, so this must not happen).
fn attach_transcript_uid_batch(
    batch: RecordBatch,
    source_type: &str,
    uid_map: &HashMap<String, u32>,
) -> Result<RecordBatch> {
    let schema = batch.schema();
    let sid_idx = schema.index_of("stable_id").map_err(|_| {
        DataFusionError::Execution("transcript batch missing stable_id column".to_string())
    })?;
    let sid = batch.column(sid_idx);

    let mut uids = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        let Some(id) = string_at(sid.as_ref(), row) else {
            return Err(DataFusionError::Execution(
                "transcript stable_id is null".to_string(),
            ));
        };
        let uid = *uid_map.get(&id).ok_or_else(|| {
            DataFusionError::Execution(format!("transcript stable_id {id} missing from uid map"))
        })?;
        uids.push(uid);
    }

    let mut fields = schema
        .fields()
        .iter()
        .map(|field| field.as_ref().clone())
        .collect::<Vec<Field>>();
    fields.push(Field::new("transcript_uid", DataType::UInt32, false));
    let mut columns = batch.columns().to_vec();
    columns.push(Arc::new(UInt32Array::from(uids)) as ArrayRef);

    let new_schema = with_cache_source_metadata(&Schema::new(fields), source_type);
    RecordBatch::try_new(Arc::new(new_schema), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to attach transcript_uid column: {err}"))
    })
}

/// Position-sliced translation_sift schema: one row per `(transcript, position)`.
/// `key = (transcript_uid << 32) | position`; small `sift`/`poly` payloads use
/// miniblock + zstd (v2.1), not FullZip.
fn compact_translation_sift_position_schema(source_type: &str) -> Schema {
    let key = Field::new("key", DataType::UInt64, false);
    let sift = Field::new("sift", DataType::Binary, false);
    let poly = Field::new("poly", DataType::Binary, false);
    let schema = with_cache_source_metadata(&Schema::new(vec![key, sift, poly]), source_type);
    let mut md = schema.metadata().clone();
    md.insert(
        crate::cache_common::SIFT_BLOB_VERSION_KEY.to_string(),
        "2".to_string(),
    );
    Schema::new_with_metadata(schema.fields().clone(), md)
}

/// Explode each transcript's SIFT/PolyPhen matrices into one row per protein
/// position. `(transcript_uid, position)` is unique by construction (uid 1:1,
/// positions grouped), so `key` is a unique BTree key.
fn transform_translation_sift_position_batch(
    batch: RecordBatch,
    source_type: &str,
    uid_map: &HashMap<String, u32>,
) -> Result<RecordBatch> {
    let schema = batch.schema();
    let tid_idx = schema.index_of("transcript_id")?;
    let sift_idx = schema.index_of("sift_predictions").ok();
    let poly_idx = schema.index_of("polyphen_predictions").ok();

    let mut keys: Vec<u64> = Vec::new();
    let mut sift_blobs: Vec<Vec<u8>> = Vec::new();
    let mut poly_blobs: Vec<Vec<u8>> = Vec::new();

    for row in 0..batch.num_rows() {
        let Some(transcript_id) = string_at(batch.column(tid_idx).as_ref(), row) else {
            continue;
        };
        let uid = *uid_map.get(&transcript_id).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "translation transcript_id {transcript_id} missing from transcript uid map"
            ))
        })?;
        let mut preds = CachedPredictions::default();
        if let Some(idx) = sift_idx {
            preds.sift = read_compact_predictions(batch.column(idx).as_ref(), row);
        }
        if let Some(idx) = poly_idx {
            preds.polyphen = read_compact_predictions(batch.column(idx).as_ref(), row);
        }
        preds.sort(); // sorts by (position, amino_acid) → position-major
        append_translation_sift_position_rows(
            uid,
            &preds,
            &mut keys,
            &mut sift_blobs,
            &mut poly_blobs,
        )?;
    }

    let target_schema = Arc::new(compact_translation_sift_position_schema(source_type));
    RecordBatch::try_new(
        target_schema,
        vec![
            Arc::new(UInt64Array::from(keys)) as ArrayRef,
            Arc::new(BinaryArray::from_iter_values(sift_blobs)) as ArrayRef,
            Arc::new(BinaryArray::from_iter_values(poly_blobs)) as ArrayRef,
        ],
    )
    .map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to build position-sliced translation_sift batch: {err}"
        ))
    })
}

/// Linear two-pointer merge over `preds.sift`/`preds.polyphen` (each sorted by
/// position) emitting one `(key, sift, poly)` row per distinct position. Position
/// payloads drop the position field (it lives in `key`).
fn append_translation_sift_position_rows(
    uid: u32,
    preds: &CachedPredictions,
    keys: &mut Vec<u64>,
    sift_blobs: &mut Vec<Vec<u8>>,
    poly_blobs: &mut Vec<Vec<u8>>,
) -> Result<()> {
    let sift = &preds.sift;
    let poly = &preds.polyphen;
    let mut si = 0usize;
    let mut pi = 0usize;
    while si < sift.len() || pi < poly.len() {
        let next_pos = match (sift.get(si), poly.get(pi)) {
            (Some(s), Some(p)) => s.position.min(p.position),
            (Some(s), None) => s.position,
            (None, Some(p)) => p.position,
            (None, None) => break,
        };
        if next_pos < 0 {
            return Err(DataFusionError::Execution(format!(
                "negative protein position {next_pos} in translation_sift"
            )));
        }
        let s_start = si;
        while si < sift.len() && sift[si].position == next_pos {
            si += 1;
        }
        let p_start = pi;
        while pi < poly.len() && poly[pi].position == next_pos {
            pi += 1;
        }
        keys.push(((uid as u64) << 32) | (next_pos as u64));
        sift_blobs.push(crate::cache_common::serialize_position_entries_v2(
            &sift[s_start..si],
            crate::cache_common::SIFT_CODEC,
        )?);
        poly_blobs.push(crate::cache_common::serialize_position_entries_v2(
            &poly[p_start..pi],
            crate::cache_common::POLY_CODEC,
        )?);
    }
    Ok(())
}

fn drop_row_number_batch(batch: RecordBatch) -> Result<RecordBatch> {
    let schema = batch.schema();
    let fields = schema
        .fields()
        .iter()
        .filter(|field| field.name() != "_rn")
        .map(|field| field.as_ref().clone())
        .collect::<Vec<_>>();
    let columns = schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(_, field)| field.name() != "_rn")
        .map(|(index, _)| batch.column(index).clone())
        .collect::<Vec<_>>();
    RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to drop translation row number: {err}"))
    })
}

fn project_batch_to_schema(batch: RecordBatch, target_schema: SchemaRef) -> Result<RecordBatch> {
    let source_schema = batch.schema();
    let mut columns = Vec::<ArrayRef>::with_capacity(target_schema.fields().len());
    for field in target_schema.fields() {
        let (index, _) = source_schema
            .column_with_name(field.name())
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "column '{}' not found in source batch",
                    field.name()
                ))
            })?;
        columns.push(batch.column(index).clone());
    }
    RecordBatch::try_new(target_schema, columns)
        .map_err(|err| DataFusionError::Execution(format!("failed to project Lance batch: {err}")))
}

fn chroms_from_schema(schema: &SchemaRef) -> Vec<String> {
    let mut chroms: Vec<String> = schema
        .metadata()
        .get("bio.vep.chromosomes")
        .and_then(|json| serde_json::from_str(json).ok())
        .unwrap_or_else(|| {
            DEFAULT_CHROMS
                .iter()
                .map(|chrom| (*chrom).to_string())
                .collect()
        });
    chroms.sort_by(|left, right| chrom_cache_order_key(left).cmp(&chrom_cache_order_key(right)));
    chroms
}

fn chrom_cache_order_key(chrom: &str) -> (usize, &str) {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Some(index) = DEFAULT_CHROMS.iter().position(|known| *known == bare) {
        (index, "")
    } else {
        (DEFAULT_CHROMS.len(), chrom)
    }
}

fn entity_output_name(kind: EnsemblEntityKind) -> &'static str {
    match kind {
        EnsemblEntityKind::Variation => "variation",
        EnsemblEntityKind::Transcript => "transcript",
        EnsemblEntityKind::Exon => "exon",
        EnsemblEntityKind::Translation => "translation",
        EnsemblEntityKind::RegulatoryFeature => "regulatory",
        EnsemblEntityKind::MotifFeature => "motif",
    }
}

fn entity_table_name(kind: EnsemblEntityKind) -> &'static str {
    match kind {
        EnsemblEntityKind::Variation => "var",
        EnsemblEntityKind::Transcript => "tx",
        EnsemblEntityKind::Exon => "exon",
        EnsemblEntityKind::Translation => "tl",
        EnsemblEntityKind::RegulatoryFeature => "reg",
        EnsemblEntityKind::MotifFeature => "motif",
    }
}

fn transcript_region_start_expr(start_col: &str) -> String {
    format!(
        "(CAST(FLOOR(({start_col} - 1) / {VEP_CACHE_REGION_SIZE_BP}.0) AS BIGINT) * {VEP_CACHE_REGION_SIZE_BP} + 1)"
    )
}

fn translation_source_region_preference_expr(start_col: &str, source_file_col: &str) -> String {
    let region_start = transcript_region_start_expr(start_col);
    let region_end = format!("({region_start} + {} - 1)", VEP_CACHE_REGION_SIZE_BP);
    format!(
        "CASE WHEN {source_file_col} LIKE CONCAT('%/', CAST({region_start} AS VARCHAR), '-', CAST({region_end} AS VARCHAR), '.gz') THEN 0 ELSE 1 END"
    )
}

fn build_translation_dedup_query_with_where_clause(table_name: &str, where_clause: &str) -> String {
    let source_pref = translation_source_region_preference_expr("start", "source_file");
    format!(
        "SELECT * FROM (\
            SELECT *, ROW_NUMBER() OVER (\
                PARTITION BY chrom, transcript_id \
                ORDER BY {source_pref}, cdna_coding_start NULLS LAST, source_file\
            ) AS _rn \
            FROM {table_name}{where_clause}\
        ) WHERE _rn = 1"
    )
}

fn sql_escape_literal(value: &str) -> String {
    value.replace('\'', "''")
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{
        ArrayRef, BinaryArray, Int8Array, Int32Array, Int64Array, StringArray, UInt32Array,
    };
    use datafusion::arrow::compute::SortOptions;
    use datafusion::arrow::datatypes::{DataType, Field};
    use datafusion::physical_expr::expressions::Column;
    use datafusion::physical_expr_common::sort_expr::{LexOrdering, PhysicalSortExpr};
    #[allow(deprecated)]
    use datafusion::physical_plan::coalesce_batches::CoalesceBatchesExec;
    use datafusion::physical_plan::sorts::sort_preserving_merge::SortPreservingMergeExec;
    use datafusion::physical_plan::test::TestMemoryExec;

    #[test]
    fn variation_transform_drops_legacy_helpers_and_casts_positions() {
        let mut fields = Vec::new();
        let mut columns = Vec::<ArrayRef>::new();
        for name in VARIATION_REQUIRED_COLUMNS {
            match *name {
                "start" | "end" => {
                    fields.push(Field::new(*name, DataType::Int64, false));
                    columns.push(Arc::new(Int64Array::from(vec![10])) as ArrayRef);
                }
                "failed" => {
                    fields.push(Field::new(*name, DataType::Int8, true));
                    columns.push(
                        Arc::new(datafusion::arrow::array::Int8Array::from(vec![Some(0)]))
                            as ArrayRef,
                    );
                }
                _ => {
                    fields.push(Field::new(*name, DataType::Utf8, true));
                    columns.push(Arc::new(StringArray::from(vec![Some("value")])) as ArrayRef);
                }
            }
        }
        fields.push(Field::new("position_key", DataType::UInt64, false));
        columns.push(Arc::new(datafusion::arrow::array::UInt64Array::from(vec![10])) as ArrayRef);
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).unwrap();

        // start=10 is warm; the warm pass (keep_tier=0) keeps it and tags tier=0.
        let warm_starts = BTreeSet::from([10i64]);
        let warm =
            transform_variation_tier_batch(batch.clone(), "ensembl", &warm_starts, 0).unwrap();
        let schema = warm.schema();

        assert_eq!(
            schema.field_with_name("start").unwrap().data_type(),
            &DataType::UInt32
        );
        assert_eq!(
            schema.field_with_name("end").unwrap().data_type(),
            &DataType::UInt32
        );
        assert_eq!(
            schema.field_with_name("tier").unwrap().data_type(),
            &DataType::Int8
        );
        assert!(schema.field_with_name("position_key").is_err());
        assert_eq!(
            schema
                .metadata()
                .get(crate::cache_source::CACHE_SOURCE_METADATA_KEY),
            Some(&"ensembl".to_string())
        );
        assert_eq!(warm.num_rows(), 1);
        let tier = warm
            .column(schema.index_of("tier").unwrap())
            .as_any()
            .downcast_ref::<Int8Array>()
            .unwrap();
        assert_eq!(tier.value(0), 0);

        // The cold pass (keep_tier=1) drops the warm row entirely.
        let cold = transform_variation_tier_batch(batch, "ensembl", &warm_starts, 1).unwrap();
        assert_eq!(cold.num_rows(), 0);
    }

    #[test]
    fn assign_transcript_uids_is_dense_and_unique() {
        let mut ids = vec![
            "ENST2".to_string(),
            "ENST1".to_string(),
            "ENST1".to_string(),
        ];
        ids.sort_unstable();
        ids.dedup();
        let map = assign_transcript_uids(&ids);
        assert_eq!(map.len(), 2);
        assert_eq!(map["ENST1"], 0);
        assert_eq!(map["ENST2"], 1);
    }

    #[test]
    fn attach_transcript_uid_appends_column_and_errors_on_missing() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            Arc::clone(&schema),
            vec![
                Arc::new(StringArray::from(vec!["ENST1", "ENST2"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![10_i64, 20])) as ArrayRef,
            ],
        )
        .unwrap();
        let uid_map: HashMap<String, u32> =
            [("ENST1".to_string(), 0u32), ("ENST2".to_string(), 1u32)]
                .into_iter()
                .collect();

        let out = attach_transcript_uid_batch(batch, "ensembl", &uid_map).unwrap();
        let out_schema = out.schema();
        let uid_idx = out_schema.index_of("transcript_uid").unwrap();
        assert_eq!(out_schema.field(uid_idx).data_type(), &DataType::UInt32);
        let uids = out
            .column(uid_idx)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        assert_eq!(uids.value(0), 0);
        assert_eq!(uids.value(1), 1);
        assert_eq!(
            out_schema
                .metadata()
                .get(crate::cache_source::CACHE_SOURCE_METADATA_KEY),
            Some(&"ensembl".to_string())
        );

        // An id absent from the map is a hard error (consistency guard).
        let orphan = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST_MISSING"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1_i64])) as ArrayRef,
            ],
        )
        .unwrap();
        assert!(attach_transcript_uid_batch(orphan, "ensembl", &uid_map).is_err());
    }

    #[test]
    fn translation_dedup_query_stays_chrom_scoped() {
        let query = build_translation_dedup_query_with_where_clause("tl", " WHERE chrom = 'chr1'");
        assert!(query.contains("PARTITION BY chrom, transcript_id"));
        assert!(query.contains("WHERE chrom = 'chr1'"));
    }

    #[test]
    fn transcript_context_query_uses_provider_schema() {
        let schema = Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("source_file", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
        ]);

        let query = build_export_query(
            EnsemblEntityKind::Transcript,
            "tx",
            Some("1"),
            Some(&schema),
        );

        assert!(query.contains("ROW_NUMBER()"));
        assert!(query.contains("PARTITION BY chrom, stable_id"));
        assert!(query.contains("WHERE chrom = '1'"));
    }

    #[test]
    fn chroms_from_schema_uses_old_cache_order() {
        let metadata = [(
            "bio.vep.chromosomes".to_string(),
            serde_json::json!(["10", "2", "1", "X", "MT", "HSCHR6_MHC_COX_CTG1"]).to_string(),
        )]
        .into_iter()
        .collect();
        let schema = Arc::new(Schema::new_with_metadata(Vec::<Field>::new(), metadata));

        assert_eq!(
            chroms_from_schema(&schema),
            ["1", "2", "10", "X", "MT", "HSCHR6_MHC_COX_CTG1"]
        );
    }

    #[test]
    fn filter_chroms_matches_bare_labels_and_passes_through_when_none() {
        let all = || {
            vec![
                "1".to_string(),
                "2".to_string(),
                "X".to_string(),
                "chr10".to_string(),
            ]
        };
        // None → unchanged.
        assert_eq!(filter_chroms(all(), None), all());
        // `chr1` matches the bare `1`; `2` matches itself; `chr10` matches `10`.
        let allow = vec!["chr1".to_string(), "2".to_string(), "10".to_string()];
        assert_eq!(
            filter_chroms(all(), Some(&allow)),
            ["1".to_string(), "2".to_string(), "chr10".to_string()]
        );
        // No overlap → empty.
        let none = vec!["chrY".to_string()];
        assert!(filter_chroms(all(), Some(&none)).is_empty());
    }

    #[test]
    fn translation_sift_position_rows_group_by_position_with_unique_keys() {
        use crate::cache_common::{SiftBlobVersion, deserialize_position_predictions_versioned};
        use crate::transcript_consequence::CompactPrediction;

        // sift at positions 1 & 2, polyphen only at position 2.
        let mut preds = CachedPredictions {
            sift: vec![
                CompactPrediction {
                    position: 1,
                    amino_acid: 0,
                    prediction: 0,
                    score: 0.3,
                },
                CompactPrediction {
                    position: 2,
                    amino_acid: 2,
                    prediction: 1,
                    score: 0.01,
                },
            ],
            polyphen: vec![CompactPrediction {
                position: 2,
                amino_acid: 2,
                prediction: 4,
                score: 0.9,
            }],
        };
        preds.sort();

        let uid = 7u32;
        let mut keys = Vec::new();
        let mut sift_blobs = Vec::new();
        let mut poly_blobs = Vec::new();
        append_translation_sift_position_rows(
            uid,
            &preds,
            &mut keys,
            &mut sift_blobs,
            &mut poly_blobs,
        )
        .unwrap();

        // One row per distinct position (1 and 2), keyed by (uid<<32)|position.
        assert_eq!(
            keys,
            vec![((uid as u64) << 32) | 1, ((uid as u64) << 32) | 2]
        );
        // Keys are unique.
        let mut sorted = keys.clone();
        sorted.dedup();
        assert_eq!(sorted.len(), keys.len());

        // Position 1: sift only, no polyphen.
        let p1 = deserialize_position_predictions_versioned(
            1,
            &sift_blobs[0],
            &poly_blobs[0],
            SiftBlobVersion::V2DivIndex,
        )
        .unwrap();
        assert_eq!(p1.sift.len(), 1);
        assert!(p1.polyphen.is_empty());
        assert_eq!(p1.sift[0].position, 1);

        // Position 2: both sift and polyphen.
        let p2 = deserialize_position_predictions_versioned(
            2,
            &sift_blobs[1],
            &poly_blobs[1],
            SiftBlobVersion::V2DivIndex,
        )
        .unwrap();
        assert_eq!(p2.sift.len(), 1);
        assert_eq!(p2.polyphen.len(), 1);
        assert_eq!(p2.polyphen[0].prediction, 4);

        // Schema is key/sift/poly.
        let schema = compact_translation_sift_position_schema("ensembl");
        assert_eq!(schema.field(0).name(), "key");
        assert_eq!(schema.field(0).data_type(), &DataType::UInt64);
        assert_eq!(schema.field(1).name(), "sift");
        assert_eq!(schema.field(2).name(), "poly");
    }

    #[test]
    fn append_translation_sift_position_rows_emits_v2() {
        use crate::cache_common::{POLY_CODEC, SIFT_CODEC, deserialize_position_entries_v2};
        use crate::transcript_consequence::CompactPrediction;
        let preds = CachedPredictions {
            sift: vec![
                CompactPrediction {
                    position: 1,
                    amino_acid: 0,
                    prediction: 1,
                    score: 0.03,
                },
                CompactPrediction {
                    position: 1,
                    amino_acid: 5,
                    prediction: 0,
                    score: 1.0,
                },
            ],
            polyphen: vec![CompactPrediction {
                position: 1,
                amino_acid: 0,
                prediction: 2,
                score: 0.998,
            }],
        };
        let mut keys = Vec::new();
        let mut sift_blobs = Vec::new();
        let mut poly_blobs = Vec::new();
        append_translation_sift_position_rows(
            42,
            &preds,
            &mut keys,
            &mut sift_blobs,
            &mut poly_blobs,
        )
        .unwrap();
        assert_eq!(keys, vec![(42u64 << 32) | 1]);
        // v2 sift: 2 entries * 3 bytes
        assert_eq!(sift_blobs[0].len(), 6);
        let s = deserialize_position_entries_v2(1, &sift_blobs[0], SIFT_CODEC).unwrap();
        assert_eq!(s[1].score.to_bits(), 1.0f32.to_bits());
        let p = deserialize_position_entries_v2(1, &poly_blobs[0], POLY_CODEC).unwrap();
        assert_eq!(p[0].score.to_bits(), 0.998f32.to_bits());
    }

    #[test]
    fn sift_position_schema_carries_v2_flag() {
        let schema = compact_translation_sift_position_schema("ensembl");
        assert_eq!(
            schema
                .metadata()
                .get(crate::cache_common::SIFT_BLOB_VERSION_KEY)
                .map(String::as_str),
            Some("2")
        );
    }
}
