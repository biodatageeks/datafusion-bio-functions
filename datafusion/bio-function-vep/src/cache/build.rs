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
};
use futures::StreamExt;
use log::info;

use crate::cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest, canonical_chrom_label,
};
use crate::cache::schema::{
    VARIATION_FORBIDDEN_COLUMNS, VARIATION_OPTIONAL_COLUMNS, VARIATION_REQUIRED_COLUMNS,
    variation_projected_schema, with_cache_identity_metadata,
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
    /// Release independently established from the raw cache and validated
    /// against the compiled vepyr support matrix before any shard is written.
    pub cache_version: String,
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

    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Variation, "var")?;

    // Derive the output schema once from the source plan schema and reuse it for
    // every written batch so the shard carries one stable Arrow schema.
    let schema_plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let projected = variation_output_schema(schema_plan.schema().as_ref(), &source_type)?;
    let out_schema = Arc::new(with_cache_identity_metadata(
        projected.as_ref(),
        &source_type,
        &options.cache_version,
    ));

    // Single-pass execution: stream the DataFusion query once instead of multiple passes.
    // Warm tier batches (rare: ~3-5% of variants, e.g. ~30-40 MB RAM for chr20) are buffered
    // in memory, while cold tier batches are streamed directly to an Arrow IPC temp file.
    // Finally, warm batches are written first to the output Parquet, followed by cold batches.
    let plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

    let rows = if warm_starts.is_empty() {
        // No warm variants: write directly to final parquet shard.
        let mut writer = VariationParquetShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
        while let Some(batch) = stream.next().await {
            let (_, cold_opt) =
                transform_variation_split_batch(batch?, &source_type, &warm_starts)?;
            if let Some(cold) = cold_opt {
                let encoded = encode_variation_batch(&cold, &out_schema)?;
                writer.write(&encoded)?;
            }
        }
        writer.finish()?
    } else {
        use datafusion::arrow::ipc::reader::FileReader;
        use datafusion::arrow::ipc::writer::FileWriter;

        let cold_tmp_path = shard_path.with_extension("cold_ipc.tmp");
        let cold_file = std::fs::File::create(&cold_tmp_path).map_err(|err| {
            DataFusionError::Execution(format!("failed to create '{cold_tmp_path:?}': {err}"))
        })?;
        let mut cold_ipc_writer = FileWriter::try_new(cold_file, &out_schema)
            .map_err(|err| DataFusionError::Execution(format!("open Arrow IPC writer: {err}")))?;

        let mut warm_batches: Vec<RecordBatch> = Vec::new();
        let mut cold_batches_count: usize = 0;

        while let Some(batch) = stream.next().await {
            let (warm_opt, cold_opt) =
                transform_variation_split_batch(batch?, &source_type, &warm_starts)?;
            if let Some(warm) = warm_opt {
                let encoded = encode_variation_batch(&warm, &out_schema)?;
                warm_batches.push(encoded);
            }
            if let Some(cold) = cold_opt {
                let encoded = encode_variation_batch(&cold, &out_schema)?;
                cold_ipc_writer.write(&encoded).map_err(|err| {
                    DataFusionError::Execution(format!("write Arrow IPC batch: {err}"))
                })?;
                cold_batches_count += 1;
            }
        }
        cold_ipc_writer.finish().map_err(|err| {
            DataFusionError::Execution(format!("finish Arrow IPC writer: {err}"))
        })?;
        drop(cold_ipc_writer);

        let mut writer = VariationParquetShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
        for warm in warm_batches {
            writer.write(&warm)?;
        }
        if cold_batches_count > 0 {
            let file = std::fs::File::open(&cold_tmp_path).map_err(|err| {
                DataFusionError::Execution(format!("open '{cold_tmp_path:?}': {err}"))
            })?;
            let reader = FileReader::try_new(file, None)
                .map_err(|err| DataFusionError::Execution(format!("open Arrow IPC reader: {err}")))?;
            for batch in reader {
                let batch = batch.map_err(|err| DataFusionError::Execution(format!("Arrow IPC read: {err}")))?;
                writer.write(&batch)?;
            }
        }
        let _ = std::fs::remove_file(&cold_tmp_path);
        writer.finish()?
    };
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
) -> Result<HashSet<i64>> {
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
    let warm_set = select_warm_positions(
        rows,
        WARM_AF_THRESHOLD,
        WARM_POSITION_RADIUS,
    );
    Ok(warm_set.into_iter().collect())
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
        let out = transform_translation_sift_position_batch(
            batch,
            &source_type,
            &options.cache_version,
            &uid_map,
        )?;
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
        .with_cache_source_type(options.cache_source_type)
        .with_expected_cache_version(&options.cache_version);
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
    let cache_version = options.cache_version.clone();
    let ctx = make_ctx_and_register(options, kind, table_name)?;

    let rows = if kind == EnsemblEntityKind::Transcript {
        let uid_map = Arc::new(load_transcript_uid_map(options, source_chrom).await?);
        write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
            attach_transcript_uid_batch(batch, &source_type, &cache_version, &uid_map)
        })
        .await?
    } else {
        write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
            attach_schema_metadata_to_batch(batch, &source_type, &cache_version)
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

    let target_schema =
        translation_core_schema(false, options.cache_source_type, &options.cache_version);
    let source_type = options.cache_source_type.as_str().to_string();
    let cache_version = options.cache_version.clone();
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Translation, "tl")?;
    let rows = write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
        let batch = drop_row_number_batch(batch)?;
        let batch = project_batch_to_schema(batch, Arc::clone(&target_schema))?;
        attach_schema_metadata_to_batch(batch, &source_type, &cache_version)
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

// ---------------------------------------------------------------------------
// Concurrent build: LPT-ordered, weighted-semaphore-bounded scheduling across
// every (entity, chrom) build primitive. `build_all()` stays the sequential
// default (`build_concurrency == 1`); `CacheBuilder::with_build_concurrency`
// routes here for `> 1`.
// ---------------------------------------------------------------------------

/// Approximate GRCh38 chromosome lengths (bp), used only to order units
/// longest-first (LPT) before scheduling. Exact per-unit cost isn't known
/// upfront; chromosome length is a cheap, good-enough proxy for it.
fn approx_chrom_length_bp(chrom: &str) -> u64 {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    match bare {
        "1" => 248_956_422,
        "2" => 242_193_529,
        "3" => 198_295_559,
        "4" => 190_214_555,
        "5" => 181_538_259,
        "6" => 170_805_979,
        "7" => 159_345_973,
        "8" => 145_138_636,
        "9" => 138_394_717,
        "10" => 133_797_422,
        "11" => 135_086_622,
        "12" => 133_275_309,
        "13" => 114_364_328,
        "14" => 107_043_718,
        "15" => 101_991_189,
        "16" => 90_338_345,
        "17" => 83_257_441,
        "18" => 80_373_285,
        "19" => 58_617_616,
        "20" => 64_444_167,
        "21" => 46_709_983,
        "22" => 50_818_468,
        "X" => 156_040_895,
        "Y" => 57_227_415,
        "MT" | "M" => 16_569,
        _ => 100_000_000, // unknown contig: assume mid-size rather than starve it
    }
}

const VARIATION_REGION_CHUNK_BP: i64 = 16_000_000;

/// One independently-schedulable cache-build task.
enum BuildUnit {
    VariationRegion {
        chrom: String,
        part_idx: usize,
        total_parts: usize,
        start_bp: i64,
        end_bp: i64,
    },
    Context(EnsemblEntityKind, String),
    TranslationCore(String),
    TranslationSift(String),
}

pub struct VariationRegionResult {
    pub chrom: String,
    pub part_idx: usize,
    pub total_parts: usize,
    pub warm_ipc_path: Option<PathBuf>,
    pub cold_ipc_path: Option<PathBuf>,
    pub warm_rows: usize,
    pub cold_rows: usize,
}

enum UnitRunResult {
    Completed(&'static str, ChromDatasetEntry),
    VariationPart(VariationRegionResult),
}

impl BuildUnit {
    fn cost(&self) -> u64 {
        match self {
            BuildUnit::VariationRegion { start_bp, end_bp, .. } => {
                let len = if *end_bp == i64::MAX {
                    VARIATION_REGION_CHUNK_BP as u64
                } else {
                    (*end_bp - *start_bp + 1).max(0) as u64
                };
                len.saturating_mul(100)
            }
            BuildUnit::Context(kind, chrom) => {
                let chrom_len = approx_chrom_length_bp(chrom);
                let weight = match kind {
                    EnsemblEntityKind::Transcript => 10,
                    EnsemblEntityKind::Exon => 5,
                    EnsemblEntityKind::RegulatoryFeature => 2,
                    EnsemblEntityKind::MotifFeature => 1,
                    _ => 1,
                };
                chrom_len.saturating_mul(weight)
            }
            BuildUnit::TranslationCore(chrom) => approx_chrom_length_bp(chrom).saturating_mul(2),
            BuildUnit::TranslationSift(chrom) => approx_chrom_length_bp(chrom).saturating_mul(2),
        }
    }

    /// Semaphore permits this unit consumes.
    fn weight(&self, var_weight: u32) -> u32 {
        match self {
            BuildUnit::VariationRegion { .. } => var_weight,
            _ => 1,
        }
    }

    async fn run(
        self,
        options: Arc<CacheBuildOptions>,
        var_out_schema: Option<Arc<Schema>>,
    ) -> Result<UnitRunResult> {
        match self {
            BuildUnit::VariationRegion {
                chrom,
                part_idx,
                total_parts,
                start_bp,
                end_bp,
            } => {
                let schema = var_out_schema.expect("var_out_schema present for VariationRegion");
                let res = build_parquet_variation_region(
                    options,
                    schema,
                    chrom,
                    part_idx,
                    total_parts,
                    start_bp,
                    end_bp,
                )
                .await?;
                Ok(UnitRunResult::VariationPart(res))
            }
            BuildUnit::Context(kind, chrom) => {
                let entity = entity_output_name(kind);
                let entry = build_parquet_context_entity_chrom(&options, kind, &chrom).await?;
                Ok(UnitRunResult::Completed(entity, entry))
            }
            BuildUnit::TranslationCore(chrom) => {
                let entry = build_parquet_translation_core_chrom(&options, &chrom).await?;
                Ok(UnitRunResult::Completed("translation_core", entry))
            }
            BuildUnit::TranslationSift(chrom) => {
                let entry = build_parquet_translation_sift_chrom(&options, &chrom).await?;
                Ok(UnitRunResult::Completed("translation_sift", entry))
            }
        }
    }
}

async fn collect_warm_starts_region(
    options: &CacheBuildOptions,
    source_chrom: &str,
    start_bp: i64,
    end_bp: i64,
) -> Result<HashSet<i64>> {
    let ctx = make_ctx_and_register_partitions(options, EnsemblEntityKind::Variation, "var", 1)?;
    let scan_start = (start_bp - WARM_POSITION_RADIUS as i64).max(1);
    let scan_end = end_bp.saturating_add(WARM_POSITION_RADIUS as i64);
    let end_filter = if end_bp < i64::MAX - 20000 {
        format!(" AND start <= {scan_end} AND end <= {}", scan_end + 10000)
    } else {
        String::new()
    };
    let sql = format!(
        "SELECT start, minor_allele_freq, \"AF\", \"gnomADg\", \"gnomADe\" FROM var \
         WHERE chrom = '{}' AND start >= {scan_start}{end_filter}",
        sql_escape_literal(source_chrom)
    );
    let df = ctx.sql(&sql).await?;
    let mut stream = df.execute_stream().await?;
    let mut rows = Vec::new();
    while let Some(batch) = stream.next().await {
        collect_position_frequencies(&batch?, &mut rows)?;
    }
    let warm_set = select_warm_positions(
        rows,
        WARM_AF_THRESHOLD,
        WARM_POSITION_RADIUS,
    );
    Ok(warm_set.into_iter().collect())
}

async fn build_parquet_variation_region(
    options: Arc<CacheBuildOptions>,
    out_schema: Arc<Schema>,
    chrom: String,
    part_idx: usize,
    total_parts: usize,
    start_bp: i64,
    end_bp: i64,
) -> Result<VariationRegionResult> {
    use crate::parquet_cache::write::encode_variation_batch;
    use datafusion::arrow::ipc::writer::FileWriter;

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(&chrom);
    let warm_starts = Arc::new(collect_warm_starts_region(&options, source_chrom, start_bp, end_bp).await?);

    let ctx = make_ctx_and_register_partitions(&options, EnsemblEntityKind::Variation, "var", 1)?;
    let end_filter = if end_bp < i64::MAX - 20000 {
        format!(" AND start <= {end_bp} AND end <= {}", end_bp + 10000)
    } else {
        String::new()
    };
    let query = format!(
        "SELECT * FROM var WHERE chrom = '{}' AND start >= {start_bp}{end_filter}",
        sql_escape_literal(source_chrom)
    );
    let plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

    let entity_dir = PathBuf::from(&options.output_dir).join("variation");
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!("failed to create '{:?}': {err}", entity_dir))
    })?;

    let warm_tmp_path = entity_dir.join(format!("{chrom}.part_{part_idx}.warm_ipc.tmp"));
    let cold_tmp_path = entity_dir.join(format!("{chrom}.part_{part_idx}.cold_ipc.tmp"));
    let _ = std::fs::remove_file(&warm_tmp_path);
    let _ = std::fs::remove_file(&cold_tmp_path);

    let warm_file = std::fs::File::create(&warm_tmp_path).map_err(|err| {
        DataFusionError::Execution(format!("create warm IPC {warm_tmp_path:?}: {err}"))
    })?;
    let mut warm_ipc_writer = FileWriter::try_new(warm_file, &out_schema)
        .map_err(|err| DataFusionError::Execution(format!("open warm IPC writer: {err}")))?;

    let cold_file = std::fs::File::create(&cold_tmp_path).map_err(|err| {
        DataFusionError::Execution(format!("create cold IPC {cold_tmp_path:?}: {err}"))
    })?;
    let mut cold_ipc_writer = FileWriter::try_new(cold_file, &out_schema)
        .map_err(|err| DataFusionError::Execution(format!("open cold IPC writer: {err}")))?;

    let source_type = options.cache_source_type.as_str().to_string();
    let mut warm_rows = 0usize;
    let mut cold_rows = 0usize;

    while let Some(batch) = stream.next().await {
        let (warm_opt, cold_opt) =
            transform_variation_split_batch(batch?, &source_type, &warm_starts)?;
        if let Some(warm) = warm_opt {
            warm_rows += warm.num_rows();
            let encoded = encode_variation_batch(&warm, &out_schema)?;
            warm_ipc_writer.write(&encoded).map_err(|err| {
                DataFusionError::Execution(format!("write warm IPC batch: {err}"))
            })?;
        }
        if let Some(cold) = cold_opt {
            cold_rows += cold.num_rows();
            let encoded = encode_variation_batch(&cold, &out_schema)?;
            cold_ipc_writer.write(&encoded).map_err(|err| {
                DataFusionError::Execution(format!("write cold IPC batch: {err}"))
            })?;
        }
    }

    warm_ipc_writer.finish().map_err(|err| {
        DataFusionError::Execution(format!("finish warm IPC writer: {err}"))
    })?;
    cold_ipc_writer.finish().map_err(|err| {
        DataFusionError::Execution(format!("finish cold IPC writer: {err}"))
    })?;

    let warm_ipc_path = if warm_rows > 0 {
        Some(warm_tmp_path)
    } else {
        let _ = std::fs::remove_file(&warm_tmp_path);
        None
    };
    let cold_ipc_path = if cold_rows > 0 {
        Some(cold_tmp_path)
    } else {
        let _ = std::fs::remove_file(&cold_tmp_path);
        None
    };

    Ok(VariationRegionResult {
        chrom,
        part_idx,
        total_parts,
        warm_ipc_path,
        cold_ipc_path,
        warm_rows,
        cold_rows,
    })
}

async fn assemble_variation_shard(
    options: &CacheBuildOptions,
    entity_dir: &Path,
    chrom: &str,
    mut parts: Vec<VariationRegionResult>,
    out_schema: Arc<Schema>,
) -> Result<ChromDatasetEntry> {
    use crate::parquet_cache::write::VariationParquetShardWriter;
    use datafusion::arrow::ipc::reader::FileReader;

    parts.sort_by_key(|p| p.part_idx);
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    remove_existing_parquet_shard(&shard_path, options.overwrite)?;

    let total_rows: usize = parts.iter().map(|p| p.warm_rows + p.cold_rows).sum();
    if total_rows == 0 {
        return Ok(ChromDatasetEntry::new(manifest_chrom, file_name, 0));
    }

    let mut writer = VariationParquetShardWriter::create(&shard_path, Arc::clone(&out_schema))?;

    // Phase 1: Write all warm rows in part order (0, 1, 2, ...)
    for part in &parts {
        if let Some(ref warm_path) = part.warm_ipc_path {
            if warm_path.exists() && part.warm_rows > 0 {
                let file = std::fs::File::open(warm_path).map_err(|err| {
                    DataFusionError::Execution(format!("open warm IPC {warm_path:?}: {err}"))
                })?;
                let reader = FileReader::try_new(file, None).map_err(|err| {
                    DataFusionError::Execution(format!("read warm IPC {warm_path:?}: {err}"))
                })?;
                for batch in reader {
                    let batch = batch.map_err(|err| {
                        DataFusionError::Execution(format!("batch read: {err}"))
                    })?;
                    writer.write(&batch)?;
                }
            }
            let _ = std::fs::remove_file(warm_path);
        }
    }

    // Phase 2: Write all cold rows in part order (0, 1, 2, ...)
    for part in &parts {
        if let Some(ref cold_path) = part.cold_ipc_path {
            if cold_path.exists() && part.cold_rows > 0 {
                let file = std::fs::File::open(cold_path).map_err(|err| {
                    DataFusionError::Execution(format!("open cold IPC {cold_path:?}: {err}"))
                })?;
                let reader = FileReader::try_new(file, None).map_err(|err| {
                    DataFusionError::Execution(format!("read cold IPC {cold_path:?}: {err}"))
                })?;
                for batch in reader {
                    let batch = batch.map_err(|err| {
                        DataFusionError::Execution(format!("batch read: {err}"))
                    })?;
                    writer.write(&batch)?;
                }
            }
            let _ = std::fs::remove_file(cold_path);
        }
    }

    let written = writer.finish()?;
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, written))
}

async fn get_variation_out_schema(options: &CacheBuildOptions) -> Result<SchemaRef> {
    use crate::parquet_cache::write::variation_output_schema;
    let ctx = make_ctx_and_register_partitions(options, EnsemblEntityKind::Variation, "var", 1)?;
    let plan = ctx.sql("SELECT * FROM var LIMIT 0").await?.create_physical_plan().await?;
    let source_type = options.cache_source_type.as_str().to_string();
    let projected = variation_output_schema(plan.schema().as_ref(), &source_type)?;
    Ok(Arc::new(with_cache_identity_metadata(
        projected.as_ref(),
        &source_type,
        &options.cache_version,
    )))
}

/// Read `VEPYR_VAR_WEIGHT` (default 1): extra semaphore permits a variation
/// unit consumes, so a memory-constrained host can bound peak RSS by trading
/// away concurrency for the dominant entity specifically instead of globally.
fn variation_permit_weight() -> u32 {
    std::env::var("VEPYR_VAR_WEIGHT")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|w| *w >= 1)
        .unwrap_or(1)
}

/// Build all entities with cross-unit concurrency: every (entity, chromosome)
/// build primitive becomes one [`BuildUnit`], sorted longest-first (LPT —
/// Longest Processing Time first) so the tail (dominated by the `variation`
/// entity) starts as early as possible, then run through a `JoinSet` bounded
/// by a weighted `Semaphore`. `build_concurrency` caps how many units run at
/// once; `options.partitions` is unrelated — it only controls DataFusion's
/// per-query internal parallelism, not this cross-unit scheduling.
///
/// Each entity's manifest is written once, after all of that entity's units
/// complete, mirroring the sequential path's per-entity
/// `merge_write_to_entity_dir` call.
pub async fn build_all_units(
    options: &CacheBuildOptions,
    build_concurrency: usize,
) -> Result<Vec<EntityStats>> {
    use tokio::sync::Semaphore;
    use tokio::task::JoinSet;

    let build_concurrency = build_concurrency.max(1);
    let options = Arc::new(options.clone());
    let var_weight = variation_permit_weight();

    let entity_kinds = [
        EnsemblEntityKind::Variation,
        EnsemblEntityKind::Transcript,
        EnsemblEntityKind::Exon,
        EnsemblEntityKind::Translation,
        EnsemblEntityKind::RegulatoryFeature,
        EnsemblEntityKind::MotifFeature,
    ];

    // Discovery: cheap, sequential (schema introspection only). Entities
    // already fully built (overwrite=false) short-circuit to a placeholder
    // EntityStats exactly like `build_all()`'s per-entity dispatch does.
    let mut done: Vec<EntityStats> = Vec::new();
    let mut units: Vec<BuildUnit> = Vec::new();
    let mut entity_dirs: HashMap<&'static str, PathBuf> = HashMap::new();

    for kind in entity_kinds {
        match kind {
            EnsemblEntityKind::Variation => {
                let entity_dir = parquet_entity_dir(&options, "variation")?;
                if should_skip_parquet_entity(&entity_dir, options.overwrite) {
                    info!("variation already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: "variation".to_string(),
                        parquet_files: vec![],
                    });
                    continue;
                }
                entity_dirs.insert("variation", entity_dir);
                let chroms =
                    discover_chroms(&options, EnsemblEntityKind::Variation, "var").await?;
                for chrom in chroms {
                    let chrom_len = approx_chrom_length_bp(&chrom) as i64;
                    let chunk_size = VARIATION_REGION_CHUNK_BP;
                    let total_parts = ((chrom_len + chunk_size - 1) / chunk_size).max(1) as usize;
                    for part_idx in 0..total_parts {
                        let start_bp = (part_idx as i64) * chunk_size + 1;
                        let end_bp = if part_idx + 1 == total_parts {
                            i64::MAX
                        } else {
                            (part_idx as i64 + 1) * chunk_size
                        };
                        units.push(BuildUnit::VariationRegion {
                            chrom: chrom.clone(),
                            part_idx,
                            total_parts,
                            start_bp,
                            end_bp,
                        });
                    }
                }
            }
            EnsemblEntityKind::Translation => {
                let core_dir = parquet_entity_dir(&options, "translation_core")?;
                let sift_dir = parquet_entity_dir(&options, "translation_sift")?;
                let core_done = should_skip_parquet_entity(&core_dir, options.overwrite);
                let sift_done = should_skip_parquet_entity(&sift_dir, options.overwrite);
                if core_done {
                    info!("translation_core already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: "translation_core".to_string(),
                        parquet_files: vec![],
                    });
                } else {
                    entity_dirs.insert("translation_core", core_dir);
                }
                if sift_done {
                    info!("translation_sift already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: "translation_sift".to_string(),
                        parquet_files: vec![],
                    });
                } else {
                    entity_dirs.insert("translation_sift", sift_dir);
                }
                if !core_done || !sift_done {
                    let chroms =
                        discover_chroms(&options, EnsemblEntityKind::Translation, "tl").await?;
                    if !core_done {
                        units.extend(chroms.iter().cloned().map(BuildUnit::TranslationCore));
                    }
                    if !sift_done {
                        units.extend(chroms.into_iter().map(BuildUnit::TranslationSift));
                    }
                }
            }
            other => {
                let entity = entity_output_name(other);
                let entity_dir = parquet_entity_dir(&options, entity)?;
                if should_skip_parquet_entity(&entity_dir, options.overwrite) {
                    info!("{entity} already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: entity.to_string(),
                        parquet_files: vec![],
                    });
                    continue;
                }
                entity_dirs.insert(entity, entity_dir);
                let table_name = entity_table_name(other);
                let chroms = discover_chroms(&options, other, table_name).await?;
                units.extend(
                    chroms
                        .into_iter()
                        .map(move |chrom| BuildUnit::Context(other, chrom)),
                );
            }
        }
    }

    // LPT: longest first, so the tail (dominated by `variation`) is in
    // flight as early as possible instead of starting last.
    units.sort_by_key(|b| std::cmp::Reverse(b.cost()));

    let var_out_schema = if units.iter().any(|u| matches!(u, BuildUnit::VariationRegion { .. })) {
        Some(get_variation_out_schema(&options).await?)
    } else {
        None
    };

    let semaphore = Arc::new(Semaphore::new(build_concurrency));
    let mut join_set: JoinSet<Result<UnitRunResult>> = JoinSet::new();
    for unit in units {
        let permits = unit.weight(var_weight).min(build_concurrency as u32);
        let semaphore = Arc::clone(&semaphore);
        let options = Arc::clone(&options);
        let var_out_schema = var_out_schema.clone();
        join_set.spawn(async move {
            let _permit = semaphore
                .acquire_many_owned(permits)
                .await
                .map_err(|err| DataFusionError::Execution(format!("semaphore closed: {err}")))?;
            unit.run(options, var_out_schema).await
        });
    }

    let mut manifest_entries: HashMap<&'static str, Vec<ChromDatasetEntry>> = HashMap::new();
    let mut files: HashMap<&'static str, Vec<(String, usize)>> = HashMap::new();
    let mut variation_parts: HashMap<String, Vec<VariationRegionResult>> = HashMap::new();

    while let Some(joined) = join_set.join_next().await {
        let result = joined
            .map_err(|err| DataFusionError::Execution(format!("build task panicked: {err}")))??;
        match result {
            UnitRunResult::Completed(entity, entry) => {
                if entry.rows > 0 {
                    let entity_dir = entity_dirs
                        .get(entity)
                        .expect("entity_dirs populated for every scheduled unit's entity");
                    let shard_path = entity_dir.join(&entry.dataset);
                    files
                        .entry(entity)
                        .or_default()
                        .push((shard_path.to_string_lossy().to_string(), entry.rows));
                    manifest_entries.entry(entity).or_default().push(entry);
                }
            }
            UnitRunResult::VariationPart(part) => {
                let chrom = part.chrom.clone();
                let total_parts = part.total_parts;
                let list = variation_parts.entry(chrom.clone()).or_default();
                list.push(part);
                if list.len() == total_parts {
                    let completed_parts = variation_parts.remove(&chrom).unwrap();
                    let entity_dir = entity_dirs
                        .get("variation")
                        .expect("variation entity_dir exists");
                    let schema = var_out_schema.as_ref().expect("var_out_schema exists");
                    let entry = assemble_variation_shard(
                        &options,
                        entity_dir,
                        &chrom,
                        completed_parts,
                        Arc::clone(schema),
                    )
                    .await?;
                    if entry.rows > 0 {
                        let shard_path = entity_dir.join(&entry.dataset);
                        files
                            .entry("variation")
                            .or_default()
                            .push((shard_path.to_string_lossy().to_string(), entry.rows));
                        manifest_entries.entry("variation").or_default().push(entry);
                    }
                }
            }
        }
    }

    for (entity, entity_dir) in &entity_dirs {
        let entries = manifest_entries.remove(entity).unwrap_or_default();
        ChromManifest::new(entries).merge_write_to_entity_dir(entity_dir)?;
        done.push(EntityStats {
            entity: entity.to_string(),
            parquet_files: files.remove(*entity).unwrap_or_default(),
        });
    }

    Ok(done)
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

fn make_ctx_and_register_partitions(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
    partitions: usize,
) -> Result<SessionContext> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);
    let mut provider_options = EnsemblCacheOptions::new(&options.cache_root)
        .with_cache_source_type(options.cache_source_type)
        .with_expected_cache_version(&options.cache_version);
    provider_options.target_partitions = Some(partitions);
    let provider = EnsemblCacheTableProvider::for_entity(kind, provider_options)?;
    ctx.register_table(table_name, provider)?;
    Ok(ctx)
}

fn make_ctx_and_register(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<SessionContext> {
    let partitions = match kind {
        EnsemblEntityKind::Variation => options.partitions.max(4),
        _ => 1,
    };
    make_ctx_and_register_partitions(options, kind, table_name, partitions)
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
        if *name == "clin_sig_allele" {
            for optional_name in VARIATION_OPTIONAL_COLUMNS {
                if let Some((optional_index, _)) = schema.column_with_name(optional_name) {
                    columns.push(batch.column(optional_index).clone());
                }
            }
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

fn transform_variation_split_batch(
    batch: RecordBatch,
    source_type: &str,
    warm_starts: &HashSet<i64>,
) -> Result<(Option<RecordBatch>, Option<RecordBatch>)> {
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
        if *name == "clin_sig_allele" {
            for optional_name in VARIATION_OPTIONAL_COLUMNS {
                if let Some((optional_index, _)) = schema.column_with_name(optional_name) {
                    columns.push(batch.column(optional_index).clone());
                }
            }
        }
    }

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
    let mut warm_mask = Vec::with_capacity(batch.num_rows());
    let mut cold_mask = Vec::with_capacity(batch.num_rows());
    let mut warm_count = 0;
    let mut cold_count = 0;
    for row in 0..batch.num_rows() {
        let is_warm = !start.is_null(row) && warm_starts.contains(&start.value(row));
        if is_warm {
            tier_values.push(0i8);
            warm_mask.push(true);
            cold_mask.push(false);
            warm_count += 1;
        } else {
            tier_values.push(1i8);
            warm_mask.push(false);
            cold_mask.push(true);
            cold_count += 1;
        }
    }
    columns.push(Arc::new(Int8Array::from(tier_values)) as ArrayRef);

    let full = RecordBatch::try_new(target_schema, columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to build Lance variation batch: {err}"))
    })?;

    let warm_batch = if warm_count > 0 {
        if warm_count == batch.num_rows() {
            Some(full.clone())
        } else {
            let mask = BooleanArray::from(warm_mask);
            Some(filter_record_batch(&full, &mask).map_err(|err| {
                DataFusionError::Execution(format!("failed to filter variation warm batch: {err}"))
            })?)
        }
    } else {
        None
    };

    let cold_batch = if cold_count > 0 {
        if cold_count == batch.num_rows() {
            Some(full)
        } else {
            let mask = BooleanArray::from(cold_mask);
            Some(filter_record_batch(&full, &mask).map_err(|err| {
                DataFusionError::Execution(format!("failed to filter variation cold batch: {err}"))
            })?)
        }
    } else {
        None
    };

    Ok((warm_batch, cold_batch))
}


fn attach_schema_metadata_to_batch(
    batch: RecordBatch,
    source_type: &str,
    cache_version: &str,
) -> Result<RecordBatch> {
    let schema = with_cache_identity_metadata(batch.schema().as_ref(), source_type, cache_version);
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
    cache_version: &str,
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

    let schema_with_upstream_metadata =
        Schema::new_with_metadata(fields, schema.metadata().clone());
    let new_schema =
        with_cache_identity_metadata(&schema_with_upstream_metadata, source_type, cache_version);
    RecordBatch::try_new(Arc::new(new_schema), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to attach transcript_uid column: {err}"))
    })
}

/// Position-sliced translation_sift schema: one row per `(transcript, position)`.
/// `key = (transcript_uid << 32) | position`; small `sift`/`poly` payloads use
/// miniblock + zstd (v2.1), not FullZip.
fn compact_translation_sift_position_schema(source_type: &str, cache_version: &str) -> Schema {
    let key = Field::new("key", DataType::UInt64, false);
    let sift = Field::new("sift", DataType::Binary, false);
    let poly = Field::new("poly", DataType::Binary, false);
    let schema = with_cache_identity_metadata(
        &Schema::new(vec![key, sift, poly]),
        source_type,
        cache_version,
    );
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
    cache_version: &str,
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

    let target_schema = Arc::new(compact_translation_sift_position_schema(
        source_type,
        cache_version,
    ));
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

/// Numeric start of the 1 Mb cache region a row was read from, parsed out of
/// the `<start>-<end>.gz` source file name.
///
/// Lexicographic ordering on `source_file` is unsafe here: `10000001-11000000`
/// sorts before `9000001-10000000` as text. `TRY_CAST` yields NULL when the
/// name does not match, so a non-region source degrades to the fallback terms
/// rather than erroring.
fn translation_source_region_start_expr(source_file_col: &str) -> String {
    format!(
        "TRY_CAST(REGEXP_REPLACE({source_file_col}, '^.*/([0-9]+)-[0-9]+\\.gz$', '$1') AS BIGINT)"
    )
}

fn build_translation_dedup_query_with_where_clause(table_name: &str, where_clause: &str) -> String {
    // Ensembl VEP resolves a transcript duplicated across cache regions in
    // `merge_features`, keeping the FIRST copy encountered over a feature list
    // ordered by region, so for a whole-chromosome run the lowest region wins.
    //
    // Preferring the region that contains the transcript `start` instead picks
    // the wrong copy whenever a transcript is also present in an earlier
    // region, which changes `protein_features` order and hence the `DOMAINS`
    // CSQ field. The start-region preference is retained as a fallback for
    // rows whose source file name does not parse.
    //
    // Traceability:
    // - Ensembl VEP `AnnotationType/Transcript.pm::merge_features()`
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L246-L280>
    // - Ensembl VEP `AnnotationSource.pm::get_all_features_by_InputBuffer()`
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm#L109-L145>
    let region_start = translation_source_region_start_expr("source_file");
    let source_pref = translation_source_region_preference_expr("start", "source_file");
    format!(
        "SELECT * FROM (\
            SELECT *, ROW_NUMBER() OVER (\
                PARTITION BY chrom, transcript_id \
                ORDER BY {region_start} NULLS LAST, {source_pref}, cdna_coding_start NULLS LAST, source_file\
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
    fn variation_transform_preserves_optional_clin_sig_ref_allele_when_present() {
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
                    columns.push(Arc::new(Int8Array::from(vec![Some(0)])) as ArrayRef);
                }
                _ => {
                    fields.push(Field::new(*name, DataType::Utf8, true));
                    columns.push(Arc::new(StringArray::from(vec![Some("value")])) as ArrayRef);
                }
            }
            if *name == "clin_sig_allele" {
                fields.push(Field::new("clin_sig_ref_allele", DataType::Utf8, false));
                columns.push(Arc::new(StringArray::from(vec!["AGG"])) as ArrayRef);
            }
        }
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).unwrap();
        let transformed =
            transform_variation_tier_batch(batch, "ensembl", &BTreeSet::from([10]), 0).unwrap();
        let schema = transformed.schema();
        let field = schema.field_with_name("clin_sig_ref_allele").unwrap();
        assert_eq!(field.data_type(), &DataType::Utf8);
        assert!(field.is_nullable());
        let values = transformed
            .column(schema.index_of("clin_sig_ref_allele").unwrap())
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(values.value(0), "AGG");
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

        let out = attach_transcript_uid_batch(batch, "ensembl", "115", &uid_map).unwrap();
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
        assert_eq!(
            out_schema
                .metadata()
                .get(crate::cache_identity::CACHE_VERSION_METADATA_KEY),
            Some(&"115".to_string())
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
        assert!(attach_transcript_uid_batch(orphan, "ensembl", "115", &uid_map).is_err());
    }

    #[test]
    fn translation_dedup_query_stays_chrom_scoped() {
        let query = build_translation_dedup_query_with_where_clause("tl", " WHERE chrom = 'chr1'");
        assert!(query.contains("PARTITION BY chrom, transcript_id"));
        assert!(query.contains("WHERE chrom = 'chr1'"));
    }

    /// A transcript present in more than one cache region must resolve to the
    /// copy from the LOWEST region, matching VEP's merge_features first-wins
    /// over a region-ordered feature list. Preferring the region containing the
    /// transcript start picks the later copy and reorders protein_features,
    /// which surfaces as a DOMAINS mismatch (chrX ENST00000381077).
    #[test]
    fn translation_dedup_prefers_the_lowest_source_region() {
        let query = build_translation_dedup_query_with_where_clause("tl", "");
        let order = query
            .split("ORDER BY ")
            .nth(1)
            .expect("ORDER BY present")
            .to_string();
        let region_term = order
            .find("REGEXP_REPLACE")
            .expect("region-start term present");
        let start_pref = order
            .find("CASE WHEN")
            .expect("start-region fallback present");
        assert!(
            region_term < start_pref,
            "lowest-region term must outrank the start-region preference: {order}"
        );
        assert!(
            order.contains("TRY_CAST"),
            "must degrade on unparsable names"
        );
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
        let schema = compact_translation_sift_position_schema("ensembl", "115");
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
        let schema = compact_translation_sift_position_schema("ensembl", "115");
        assert_eq!(
            schema
                .metadata()
                .get(crate::cache_common::SIFT_BLOB_VERSION_KEY)
                .map(String::as_str),
            Some("2")
        );
    }
}
