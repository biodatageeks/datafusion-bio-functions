use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use datafusion::arrow::array::{
    Array, ArrayRef, BinaryArray, BooleanArray, Float64Array, Int8Array, Int64Array, StringArray,
    UInt32Array, UInt64Array,
};
use datafusion::arrow::compute::{cast, concat_batches, filter_record_batch};
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
use lance::dataset::WriteMode;
use lance_file::version::LanceFileVersion;
use log::info;
use tokio::sync::mpsc;

use crate::cache_builder::EntityStats;
use crate::cache_common::{
    FrequencyFields, PositionFrequency, max_global_af, select_warm_positions,
};
use crate::lance_cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest, canonical_chrom_label, dataset_dir_name,
};
use crate::lance_cache::schema::{
    VARIATION_FORBIDDEN_COLUMNS, VARIATION_REQUIRED_COLUMNS, lance_field_metadata,
    validate_variation_schema, with_cache_source_metadata, with_lance_field_metadata,
};
use crate::lance_cache::write::{
    LanceIndexKind, align_batch_to_schema, create_required_index, create_tier_bitmap_index,
    write_record_batch_stream_to_lance_with_mode_and_version,
    write_record_batch_stream_to_lance_with_version,
    write_record_batches_to_lance_with_mode_and_version,
};
use crate::{
    annotate_provider::{read_compact_predictions, string_at},
    transcript_consequence::CachedPredictions,
};

const DEFAULT_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];
const LANCE_CONTEXT_WRITE_CHUNK_ROWS: usize = 4_096;
/// Write chunk for the narrow translation datasets (translation_sift / _core).
/// Large so the dataset lands in a few big fragments that compress well; the
/// rows are small enough that buffering ~1M of them is cheap.
const LANCE_TRANSLATION_WRITE_CHUNK_ROWS: usize = 1_048_576;
const LANCE_VARIATION_STREAM_BATCH_ROWS: usize = 1_048_576;
/// Warm/cold tier selection for the single-path variation layout. A genomic
/// `start` is "warm" (tier 0) when its max global allele frequency across
/// `minor_allele_freq`, `AF`, `gnomADg`, `gnomADe` is at least this threshold;
/// the radius extends warm membership to neighboring positions that exist in the
/// data. These must match `warm_cache::split` so the Lance tier agrees with any
/// other tiered consumer.
const WARM_AF_THRESHOLD: f64 = 0.01;
const WARM_POSITION_RADIUS: i64 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LanceWriteStrategy {
    StreamFullDataset,
    ChunkedAppend { chunk_rows: usize },
}

fn lance_write_strategy(kind: EnsemblEntityKind) -> LanceWriteStrategy {
    match kind {
        EnsemblEntityKind::Variation => LanceWriteStrategy::StreamFullDataset,
        // Translation rows are narrow (translation_sift = key + two small Binary
        // blobs; translation_core is tiny), so a large write chunk is cheap in
        // memory and writes the dataset in a few big fragments. The small 4 096-row
        // chunk used for the wide entities (transcript/exon, which carry long
        // sequence strings) fragments translation_sift into thousands of tiny
        // writes, which Lance's miniblock+zstd compresses far worse (~1.4x vs
        // ~2.3x). Match variation's coalesced write unit instead.
        EnsemblEntityKind::Translation => LanceWriteStrategy::ChunkedAppend {
            chunk_rows: LANCE_TRANSLATION_WRITE_CHUNK_ROWS,
        },
        EnsemblEntityKind::Transcript
        | EnsemblEntityKind::Exon
        | EnsemblEntityKind::RegulatoryFeature
        | EnsemblEntityKind::MotifFeature => LanceWriteStrategy::ChunkedAppend {
            chunk_rows: LANCE_CONTEXT_WRITE_CHUNK_ROWS,
        },
    }
}

fn should_bypass_ordered_plan_root(write_strategy: LanceWriteStrategy) -> bool {
    matches!(write_strategy, LanceWriteStrategy::ChunkedAppend { .. })
}

pub fn lance_entity_dir_name(entity: &str) -> String {
    format!("{entity}.lance")
}

#[derive(Debug, Clone)]
pub struct LanceCacheBuildOptions {
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

pub async fn build_lance_entity(
    options: &LanceCacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<Vec<EntityStats>> {
    match kind {
        EnsemblEntityKind::Variation => build_lance_variation(options).await,
        EnsemblEntityKind::Translation => build_lance_translation(options).await,
        EnsemblEntityKind::Transcript
        | EnsemblEntityKind::Exon
        | EnsemblEntityKind::RegulatoryFeature
        | EnsemblEntityKind::MotifFeature => build_lance_context_entity(options, kind).await,
    }
}

async fn build_lance_variation(options: &LanceCacheBuildOptions) -> Result<Vec<EntityStats>> {
    let entity_dir = entity_dir(options, "variation")?;
    if should_skip_lance_entity(&entity_dir, options.overwrite) {
        info!("variation.lance already exists, skipping (use overwrite to rebuild)");
        return Ok(vec![EntityStats {
            entity: "variation.lance".to_string(),
            parquet_files: vec![],
        }]);
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Variation, "var").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();

    for chrom in chroms {
        let entry = build_lance_variation_chrom(options, &chrom).await?;
        if entry.rows > 0 {
            let dataset_path = entity_dir.join(&entry.dataset);
            files.push((dataset_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: "variation.lance".to_string(),
        parquet_files: files,
    }])
}

pub async fn build_lance_variation_chrom(
    options: &LanceCacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let entity_dir = entity_dir(options, "variation")?;
    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let query = build_export_query(
        EnsemblEntityKind::Variation,
        "var",
        Some(source_chrom),
        None,
    );
    let manifest_chrom = canonical_chrom_label(chrom);
    let dataset_name = dataset_dir_name(&manifest_chrom);
    let dataset_path = entity_dir.join(&dataset_name);
    if dataset_path.exists() {
        if !options.overwrite {
            // The tiered writer's first pass uses `WriteMode::Overwrite`, so a
            // rebuild would clobber the existing chromosome dataset regardless
            // of this flag. Refuse instead of silently overwriting.
            return Err(DataFusionError::Execution(format!(
                "Lance dataset '{}' already exists and overwrite=false",
                dataset_path.display()
            )));
        }
        std::fs::remove_dir_all(&dataset_path).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to remove existing Lance dataset '{}': {err}",
                dataset_path.display()
            ))
        })?;
    }
    let source_type = options.cache_source_type.as_str().to_string();
    let warm_starts = Arc::new(collect_warm_starts(options, source_chrom).await?);
    info!(
        "Lance variation {} warm positions={} (af>={WARM_AF_THRESHOLD}, +/-{WARM_POSITION_RADIUS})",
        dataset_path.display(),
        warm_starts.len(),
    );
    let rows =
        write_variation_tiered_to_lance(options, &query, &dataset_path, warm_starts, source_type)
            .await?;
    Ok(ChromDatasetEntry::new(manifest_chrom, dataset_name, rows))
}

/// Run a cheap pre-pass over the source variation table for one chromosome to
/// classify each genomic `start` as warm (common) or cold (rare) using the same
/// allele-frequency selector as the Parquet warm tier.
async fn collect_warm_starts(
    options: &LanceCacheBuildOptions,
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

/// Write the variation dataset clustered by tier: warm (common) rows first as a
/// dense prefix, cold (rare) rows appended after. Each pass re-runs the
/// `ORDER BY start` source plan and keeps only its tier, so both tiers stay
/// `start`-sorted. The `start` BTree and `tier` bitmap indexes are created once
/// after both passes so they cover all row ids.
async fn write_variation_tiered_to_lance(
    options: &LanceCacheBuildOptions,
    query: &str,
    dataset_path: &Path,
    warm_starts: Arc<BTreeSet<i64>>,
    source_type: String,
) -> Result<usize> {
    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Variation, "var")?;
    write_variation_tiered_with_ctx(&ctx, query, dataset_path, warm_starts, source_type).await
}

async fn write_variation_tiered_with_ctx(
    ctx: &SessionContext,
    query: &str,
    dataset_path: &Path,
    warm_starts: Arc<BTreeSet<i64>>,
    source_type: String,
) -> Result<usize> {
    let passes: &[(i8, WriteMode)] = if warm_starts.is_empty() {
        &[(1, WriteMode::Overwrite)]
    } else {
        &[(0, WriteMode::Overwrite), (1, WriteMode::Append)]
    };

    let mut total_rows = 0usize;
    for &(keep_tier, mode) in passes {
        total_rows += write_variation_tier_pass(
            ctx,
            query,
            dataset_path,
            &warm_starts,
            &source_type,
            keep_tier,
            mode,
        )
        .await?;
    }

    if total_rows > 0 {
        create_required_index(dataset_path, LanceIndexKind::Start).await?;
        create_tier_bitmap_index(dataset_path).await?;
    }
    Ok(total_rows)
}

#[allow(clippy::too_many_arguments)]
async fn write_variation_tier_pass(
    ctx: &SessionContext,
    query: &str,
    dataset_path: &Path,
    warm_starts: &Arc<BTreeSet<i64>>,
    source_type: &str,
    keep_tier: i8,
    mode: WriteMode,
) -> Result<usize> {
    let df = ctx.sql(query).await?;
    let plan = df.create_physical_plan().await?;
    let output_schema = streaming_lance_output_schema(
        EnsemblEntityKind::Variation,
        plan.schema().as_ref(),
        source_type,
    )?;
    let stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

    let rows = Arc::new(AtomicUsize::new(0));
    let warm_starts = Arc::clone(warm_starts);
    let source_type = source_type.to_string();
    let transform = move |batch: RecordBatch| {
        let batch = transform_variation_tier_batch(batch, &source_type, &warm_starts, keep_tier)?;
        crate::lance_cache::af_bundle::bundle_af_columns(&batch)
    };
    let stream = transformed_lance_stream(
        stream,
        Arc::clone(&output_schema),
        Arc::clone(&rows),
        transform,
    );
    let stream =
        coalesce_record_batch_stream(stream, output_schema, LANCE_VARIATION_STREAM_BATCH_ROWS);
    write_record_batch_stream_to_lance_with_mode_and_version(
        dataset_path,
        stream,
        mode,
        // 2.2 is required for the bundled fullzip List<Utf8> AF columns (large-chunk fix).
        LanceFileVersion::V2_2,
    )
    .await?;
    Ok(rows.load(Ordering::Relaxed))
}

async fn build_lance_context_entity(
    options: &LanceCacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<Vec<EntityStats>> {
    let entity = entity_output_name(kind);
    let entity_dir = entity_dir(options, entity)?;
    if should_skip_lance_entity(&entity_dir, options.overwrite) {
        info!("{entity}.lance already exists, skipping (use overwrite to rebuild)");
        return Ok(vec![EntityStats {
            entity: lance_entity_dir_name(entity),
            parquet_files: vec![],
        }]);
    }

    let table_name = entity_table_name(kind);
    let (chroms, provider_schema) = discover_chroms_and_schema(options, kind, table_name).await?;
    let transcript_schema = (kind == EnsemblEntityKind::Transcript).then_some(provider_schema);
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    let index_kind = context_index_kind(kind);

    for chrom in chroms {
        let query = build_export_query(
            kind,
            table_name,
            Some(&chrom),
            transcript_schema.as_ref().map(|schema| schema.as_ref()),
        );
        let manifest_chrom = canonical_chrom_label(&chrom);
        let dataset_name = dataset_dir_name(&manifest_chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        let rows = if kind == EnsemblEntityKind::Transcript {
            // The transcript table carries a dense `transcript_uid` (stable-id rank)
            // used to key the position-sliced translation_sift dataset.
            let uid_map = Arc::new(load_transcript_uid_map(options, &chrom).await?);
            let source_type = options.cache_source_type.as_str().to_string();
            write_query_stream_to_lance(
                options,
                kind,
                table_name,
                &query,
                &dataset_path,
                index_kind,
                LanceFileVersion::V2_2,
                move |batch| attach_transcript_uid_batch(batch, &source_type, &uid_map),
            )
            .await?
        } else {
            let source_type = options.cache_source_type.as_str().to_string();
            write_query_stream_to_lance(
                options,
                kind,
                table_name,
                &query,
                &dataset_path,
                index_kind,
                LanceFileVersion::V2_2,
                move |batch| attach_schema_metadata_to_batch(batch, &source_type),
            )
            .await?
        };
        if rows == 0 {
            continue;
        }
        manifest_entries.push(ChromDatasetEntry::new(manifest_chrom, dataset_name, rows));
        files.push((dataset_path.to_string_lossy().to_string(), rows));
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: lance_entity_dir_name(entity),
        parquet_files: files,
    }])
}

async fn build_lance_translation(options: &LanceCacheBuildOptions) -> Result<Vec<EntityStats>> {
    let core = build_lance_translation_split(
        options,
        "translation_core",
        translation_core_schema(false, options.cache_source_type),
    )
    .await?;
    let sift = build_lance_translation_sift(options).await?;
    Ok(vec![core, sift])
}

async fn build_lance_translation_sift(options: &LanceCacheBuildOptions) -> Result<EntityStats> {
    let entity = "translation_sift";
    let entity_dir = entity_dir(options, entity)?;
    if should_skip_lance_entity(&entity_dir, options.overwrite) {
        info!("{entity}.lance already exists, skipping (use overwrite to rebuild)");
        return Ok(EntityStats {
            entity: lance_entity_dir_name(entity),
            parquet_files: vec![],
        });
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Translation, "tl").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let query = build_translation_dedup_query_with_where_clause(
            "tl",
            &format!(" WHERE chrom = '{}'", sql_escape_literal(&chrom)),
        );
        let manifest_chrom = canonical_chrom_label(&chrom);
        let dataset_name = dataset_dir_name(&manifest_chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        let source_type = options.cache_source_type.as_str().to_string();
        // Same shared uid map the transcript build uses (keyed by stable_id ==
        // translation transcript_id), so position-row keys agree with the
        // transcript table's `transcript_uid`.
        let uid_map = Arc::new(load_transcript_uid_map(options, &chrom).await?);
        let rows = write_query_stream_to_lance(
            options,
            EnsemblEntityKind::Translation,
            "tl",
            &query,
            &dataset_path,
            LanceIndexKind::SiftKey,
            LanceFileVersion::V2_1,
            move |batch| {
                let batch = drop_row_number_batch(batch)?;
                transform_translation_sift_position_batch(batch, &source_type, &uid_map)
            },
        )
        .await?;
        if rows == 0 {
            continue;
        }
        // Collapse the per-batch fragments + stale versions left by the
        // chunked append (reclaims the stale _versions churn).
        crate::lance_cache::write::compact_and_cleanup_sift_dataset(&dataset_path).await?;
        manifest_entries.push(ChromDatasetEntry::new(manifest_chrom, dataset_name, rows));
        files.push((dataset_path.to_string_lossy().to_string(), rows));
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(EntityStats {
        entity: lance_entity_dir_name(entity),
        parquet_files: files,
    })
}

async fn build_lance_translation_split(
    options: &LanceCacheBuildOptions,
    entity: &str,
    target_schema: SchemaRef,
) -> Result<EntityStats> {
    let entity_dir = entity_dir(options, entity)?;
    if should_skip_lance_entity(&entity_dir, options.overwrite) {
        info!("{entity}.lance already exists, skipping (use overwrite to rebuild)");
        return Ok(EntityStats {
            entity: lance_entity_dir_name(entity),
            parquet_files: vec![],
        });
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Translation, "tl").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let query = build_translation_dedup_query_with_where_clause(
            "tl",
            &format!(" WHERE chrom = '{}'", sql_escape_literal(&chrom)),
        );
        let manifest_chrom = canonical_chrom_label(&chrom);
        let dataset_name = dataset_dir_name(&manifest_chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        let source_type = options.cache_source_type.as_str().to_string();
        let rows = write_query_stream_to_lance(
            options,
            EnsemblEntityKind::Translation,
            "tl",
            &query,
            &dataset_path,
            LanceIndexKind::TranscriptId,
            LanceFileVersion::V2_2,
            {
                let target_schema = Arc::clone(&target_schema);
                move |batch| {
                    let batch = drop_row_number_batch(batch)?;
                    let batch = project_batch_to_schema(batch, Arc::clone(&target_schema))?;
                    attach_schema_metadata_to_batch(batch, &source_type)
                }
            },
        )
        .await?;
        if rows == 0 {
            continue;
        }
        manifest_entries.push(ChromDatasetEntry::new(manifest_chrom, dataset_name, rows));
        files.push((dataset_path.to_string_lossy().to_string(), rows));
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(EntityStats {
        entity: lance_entity_dir_name(entity),
        parquet_files: files,
    })
}

fn entity_dir(options: &LanceCacheBuildOptions, entity: &str) -> Result<PathBuf> {
    let path = Path::new(&options.output_dir).join(lance_entity_dir_name(entity));
    std::fs::create_dir_all(&path).map_err(|err| {
        DataFusionError::Execution(format!("failed to create {}: {err}", path.display()))
    })?;
    Ok(path)
}

fn should_skip_lance_entity(entity_dir: &Path, overwrite: bool) -> bool {
    !overwrite
        && entity_dir.join(CHROM_MANIFEST_FILE).exists()
        && std::fs::read_dir(entity_dir)
            .ok()
            .into_iter()
            .flat_map(|entries| entries.flatten())
            .any(|entry| entry.file_name().to_string_lossy().ends_with(".lance"))
}

async fn discover_chroms(
    options: &LanceCacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<Vec<String>> {
    discover_chroms_and_schema(options, kind, table_name)
        .await
        .map(|(chroms, _)| chroms)
}

async fn discover_chroms_and_schema(
    options: &LanceCacheBuildOptions,
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

async fn write_query_stream_to_lance<F>(
    options: &LanceCacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
    query: &str,
    dataset_path: &Path,
    index_kind: LanceIndexKind,
    data_storage_version: LanceFileVersion,
    transform: F,
) -> Result<usize>
where
    F: FnMut(RecordBatch) -> Result<RecordBatch> + Clone + Send + Sync + 'static,
{
    let ctx = make_ctx_and_register(options, kind, table_name)?;
    let df = ctx.sql(query).await?;
    let plan = df.create_physical_plan().await?;
    let write_strategy = lance_write_strategy(kind);
    let stream_output_schema = match write_strategy {
        LanceWriteStrategy::StreamFullDataset => Some(streaming_lance_output_schema(
            kind,
            plan.schema().as_ref(),
            options.cache_source_type.as_str(),
        )?),
        LanceWriteStrategy::ChunkedAppend { .. } => None,
    };
    if should_bypass_ordered_plan_root(write_strategy)
        && let Some(inner) = extract_parallel_source_plan(&plan)
    {
        let partition_count = inner.properties().partitioning.partition_count();
        info!(
            "Lance build {} {} using {partition_count} physical partitions",
            table_name,
            dataset_path.display()
        );
        return match write_strategy {
            LanceWriteStrategy::StreamFullDataset => {
                write_partitioned_plan_to_lance_streaming(
                    &ctx,
                    inner,
                    dataset_path,
                    index_kind,
                    stream_output_schema.expect("streaming output schema"),
                    data_storage_version,
                    transform,
                )
                .await
            }
            LanceWriteStrategy::ChunkedAppend { chunk_rows } => {
                write_partitioned_plan_to_lance(
                    &ctx,
                    inner,
                    dataset_path,
                    index_kind,
                    chunk_rows,
                    data_storage_version,
                    transform,
                )
                .await
            }
        };
    }

    info!(
        "Lance build {} {} using serial plan root={} partitions={}",
        table_name,
        dataset_path.display(),
        plan.name(),
        plan.properties().partitioning.partition_count(),
    );
    let stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;
    match write_strategy {
        LanceWriteStrategy::StreamFullDataset => {
            write_stream_to_lance_streaming(
                stream,
                dataset_path,
                index_kind,
                stream_output_schema.expect("streaming output schema"),
                data_storage_version,
                transform,
            )
            .await
        }
        LanceWriteStrategy::ChunkedAppend { chunk_rows } => {
            write_stream_to_lance(
                stream,
                dataset_path,
                index_kind,
                chunk_rows,
                data_storage_version,
                transform,
            )
            .await
        }
    }
}

async fn write_stream_to_lance_streaming<F>(
    stream: SendableRecordBatchStream,
    dataset_path: &Path,
    index_kind: LanceIndexKind,
    output_schema: SchemaRef,
    data_storage_version: LanceFileVersion,
    transform: F,
) -> Result<usize>
where
    F: FnMut(RecordBatch) -> Result<RecordBatch> + Send + 'static,
{
    let rows = Arc::new(AtomicUsize::new(0));
    let stream = transformed_lance_stream(
        stream,
        Arc::clone(&output_schema),
        Arc::clone(&rows),
        transform,
    );
    let stream =
        coalesce_record_batch_stream(stream, output_schema, LANCE_VARIATION_STREAM_BATCH_ROWS);
    write_record_batch_stream_to_lance_with_version(dataset_path, stream, data_storage_version)
        .await?;
    let written_rows = rows.load(Ordering::Relaxed);
    if written_rows > 0 {
        create_required_index(dataset_path, index_kind).await?;
    }
    Ok(written_rows)
}

async fn write_partitioned_plan_to_lance_streaming<F>(
    ctx: &SessionContext,
    inner: Arc<dyn ExecutionPlan>,
    dataset_path: &Path,
    index_kind: LanceIndexKind,
    output_schema: SchemaRef,
    data_storage_version: LanceFileVersion,
    transform: F,
) -> Result<usize>
where
    F: FnMut(RecordBatch) -> Result<RecordBatch> + Clone + Send + Sync + 'static,
{
    let partition_count = inner.properties().partitioning.partition_count();
    let task_ctx = ctx.task_ctx();
    let rows = Arc::new(AtomicUsize::new(0));
    let mut streams = Vec::with_capacity(partition_count);

    for partition_idx in 0..partition_count {
        let stream = inner.execute(partition_idx, Arc::clone(&task_ctx))?;
        streams.push(transformed_lance_stream(
            stream,
            Arc::clone(&output_schema),
            Arc::clone(&rows),
            transform.clone(),
        ));
    }

    let combined = futures::stream::select_all(streams);
    let combined = Box::pin(RecordBatchStreamAdapter::new(
        Arc::clone(&output_schema),
        combined,
    ));
    let combined =
        coalesce_record_batch_stream(combined, output_schema, LANCE_VARIATION_STREAM_BATCH_ROWS);
    write_record_batch_stream_to_lance_with_version(dataset_path, combined, data_storage_version)
        .await?;
    let written_rows = rows.load(Ordering::Relaxed);
    if written_rows > 0 {
        create_required_index(dataset_path, index_kind).await?;
    }
    Ok(written_rows)
}

fn transformed_lance_stream<F>(
    stream: SendableRecordBatchStream,
    output_schema: SchemaRef,
    rows: Arc<AtomicUsize>,
    mut transform: F,
) -> SendableRecordBatchStream
where
    F: FnMut(RecordBatch) -> Result<RecordBatch> + Send + 'static,
{
    let adapter_schema = Arc::clone(&output_schema);
    let stream = stream.filter_map(move |batch_result| {
        let output_schema = Arc::clone(&output_schema);
        let rows = Arc::clone(&rows);
        futures::future::ready(match batch_result {
            Ok(batch) if batch.num_rows() == 0 => None,
            Ok(batch) => match transform(batch).and_then(|batch| {
                if batch.num_rows() == 0 {
                    Ok(None)
                } else {
                    let row_count = batch.num_rows();
                    align_batch_to_schema(batch, Arc::clone(&output_schema))
                        .map(|batch| Some((batch, row_count)))
                }
            }) {
                Ok(Some((batch, row_count))) => {
                    rows.fetch_add(row_count, Ordering::Relaxed);
                    Some(Ok(batch))
                }
                Ok(None) => None,
                Err(err) => Some(Err(err)),
            },
            Err(err) => Some(Err(err)),
        })
    });
    Box::pin(RecordBatchStreamAdapter::new(adapter_schema, stream))
}

struct CoalesceRecordBatchStreamState {
    input: SendableRecordBatchStream,
    schema: SchemaRef,
    target_rows: usize,
    pending: Vec<RecordBatch>,
    pending_rows: usize,
    carry: Option<RecordBatch>,
    deferred_error: Option<DataFusionError>,
}

impl CoalesceRecordBatchStreamState {
    fn new(input: SendableRecordBatchStream, schema: SchemaRef, target_rows: usize) -> Self {
        Self {
            input,
            schema,
            target_rows: target_rows.max(1),
            pending: Vec::new(),
            pending_rows: 0,
            carry: None,
            deferred_error: None,
        }
    }

    async fn next_batch(&mut self) -> Option<Result<RecordBatch>> {
        loop {
            if let Some(error) = self.deferred_error.take() {
                return Some(Err(error));
            }
            if self.pending_rows >= self.target_rows {
                return Some(self.flush());
            }

            let batch = if let Some(batch) = self.carry.take() {
                batch
            } else {
                match self.input.next().await {
                    Some(Ok(batch)) => batch,
                    Some(Err(error)) => {
                        if self.pending_rows > 0 {
                            self.deferred_error = Some(error);
                            return Some(self.flush());
                        }
                        return Some(Err(error));
                    }
                    None if self.pending_rows > 0 => return Some(self.flush()),
                    None => return None,
                }
            };

            if batch.num_rows() == 0 {
                continue;
            }
            if let Some(flushed) = self.push_batch(batch) {
                return Some(flushed);
            }
        }
    }

    fn push_batch(&mut self, batch: RecordBatch) -> Option<Result<RecordBatch>> {
        let mut offset = 0usize;
        while offset < batch.num_rows() {
            if self.pending_rows >= self.target_rows {
                self.carry = Some(batch.slice(offset, batch.num_rows() - offset));
                return Some(self.flush());
            }
            let remaining_capacity = self.target_rows - self.pending_rows;
            let len = remaining_capacity.min(batch.num_rows() - offset);
            self.pending.push(batch.slice(offset, len));
            self.pending_rows += len;
            offset += len;
        }
        (self.pending_rows >= self.target_rows).then(|| self.flush())
    }

    fn flush(&mut self) -> Result<RecordBatch> {
        let pending = std::mem::take(&mut self.pending);
        self.pending_rows = 0;
        if pending.len() == 1 {
            Ok(pending.into_iter().next().expect("one pending batch"))
        } else {
            concat_batches(&self.schema, pending.iter()).map_err(|err| {
                DataFusionError::Execution(format!(
                    "failed to coalesce Lance output record batches: {err}"
                ))
            })
        }
    }
}

fn coalesce_record_batch_stream(
    input: SendableRecordBatchStream,
    schema: SchemaRef,
    target_rows: usize,
) -> SendableRecordBatchStream {
    let adapter_schema = Arc::clone(&schema);
    let state = CoalesceRecordBatchStreamState::new(input, schema, target_rows);
    // `.fuse()` so the stream keeps returning `None` if the downstream Lance
    // writer polls it again after exhaustion. A bare `unfold` panics on
    // poll-after-`None` ("Unfold must not be polled after it returned
    // `Poll::Ready(None)`"), which crashed the full variation build mid-stream.
    let stream = futures::stream::unfold(state, |mut state| async move {
        state.next_batch().await.map(|batch| (batch, state))
    })
    .fuse();
    Box::pin(RecordBatchStreamAdapter::new(adapter_schema, stream))
}

struct LanceChunkWriter<'a> {
    dataset_path: &'a Path,
    index_kind: LanceIndexKind,
    chunk_rows: usize,
    data_storage_version: LanceFileVersion,
    pending: Vec<RecordBatch>,
    pending_rows: usize,
    total_rows: usize,
    mode: WriteMode,
    wrote_any: bool,
}

impl<'a> LanceChunkWriter<'a> {
    fn new(
        dataset_path: &'a Path,
        index_kind: LanceIndexKind,
        chunk_rows: usize,
        data_storage_version: LanceFileVersion,
    ) -> Self {
        Self {
            dataset_path,
            index_kind,
            chunk_rows,
            data_storage_version,
            pending: Vec::new(),
            pending_rows: 0,
            total_rows: 0,
            mode: WriteMode::Overwrite,
            wrote_any: false,
        }
    }

    async fn push(&mut self, batch: RecordBatch) -> Result<()> {
        if batch.num_rows() == 0 {
            return Ok(());
        }
        let mut offset = 0usize;
        while offset < batch.num_rows() {
            if self.pending_rows >= self.chunk_rows {
                self.flush().await?;
            }
            let remaining_capacity = self.chunk_rows.saturating_sub(self.pending_rows).max(1);
            let len = remaining_capacity.min(batch.num_rows() - offset);
            self.pending_rows += len;
            self.total_rows += len;
            self.pending.push(batch.slice(offset, len));
            offset += len;

            if self.pending_rows >= self.chunk_rows {
                self.flush().await?;
            }
        }
        Ok(())
    }

    async fn flush(&mut self) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        write_record_batches_to_lance_with_mode_and_version(
            self.dataset_path,
            std::mem::take(&mut self.pending),
            self.mode,
            self.data_storage_version,
        )
        .await?;
        self.mode = WriteMode::Append;
        self.pending_rows = 0;
        self.wrote_any = true;
        Ok(())
    }

    async fn finish(mut self) -> Result<usize> {
        self.flush().await?;
        if self.wrote_any {
            create_required_index(self.dataset_path, self.index_kind).await?;
        }
        Ok(self.total_rows)
    }
}

async fn write_stream_to_lance<F>(
    mut stream: SendableRecordBatchStream,
    dataset_path: &Path,
    index_kind: LanceIndexKind,
    chunk_rows: usize,
    data_storage_version: LanceFileVersion,
    mut transform: F,
) -> Result<usize>
where
    F: FnMut(RecordBatch) -> Result<RecordBatch>,
{
    let mut writer =
        LanceChunkWriter::new(dataset_path, index_kind, chunk_rows, data_storage_version);

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }
        let batch = transform(batch)?;
        if batch.num_rows() == 0 {
            continue;
        }
        writer.push(batch).await?;
    }

    writer.finish().await
}

async fn write_partitioned_plan_to_lance<F>(
    ctx: &SessionContext,
    inner: Arc<dyn ExecutionPlan>,
    dataset_path: &Path,
    index_kind: LanceIndexKind,
    chunk_rows: usize,
    data_storage_version: LanceFileVersion,
    transform: F,
) -> Result<usize>
where
    F: FnMut(RecordBatch) -> Result<RecordBatch> + Clone + Send + Sync + 'static,
{
    let partition_count = inner.properties().partitioning.partition_count();
    let task_ctx = ctx.task_ctx();
    let (tx, mut rx) = mpsc::channel::<RecordBatch>((partition_count.max(1) * 2).max(8));
    let mut handles = tokio::task::JoinSet::new();

    for partition_idx in 0..partition_count {
        let mut stream = inner.execute(partition_idx, Arc::clone(&task_ctx))?;
        let tx = tx.clone();
        let mut transform = transform.clone();
        handles.spawn(async move {
            let mut rows = 0usize;
            while let Some(batch_result) = stream.next().await {
                let batch = batch_result?;
                if batch.num_rows() == 0 {
                    continue;
                }
                let batch = transform(batch)?;
                if batch.num_rows() == 0 {
                    continue;
                }
                rows += batch.num_rows();
                tx.send(batch).await.map_err(|_| {
                    DataFusionError::Execution(
                        "Lance partition writer stopped before all batches were sent".to_string(),
                    )
                })?;
            }
            Ok::<usize, DataFusionError>(rows)
        });
    }
    drop(tx);

    let mut writer =
        LanceChunkWriter::new(dataset_path, index_kind, chunk_rows, data_storage_version);
    while let Some(batch) = rx.recv().await {
        writer.push(batch).await?;
    }

    let mut partition_rows = 0usize;
    while let Some(result) = handles.join_next().await {
        partition_rows += result.map_err(|err| {
            DataFusionError::Execution(format!("Lance partition task failed: {err}"))
        })??;
    }

    let written_rows = writer.finish().await?;
    if written_rows != partition_rows {
        return Err(DataFusionError::Execution(format!(
            "Lance partition row count mismatch for {}: wrote {written_rows}, source partitions produced {partition_rows}",
            dataset_path.display()
        )));
    }

    Ok(written_rows)
}

fn extract_parallel_source_plan(plan: &Arc<dyn ExecutionPlan>) -> Option<Arc<dyn ExecutionPlan>> {
    if plan.properties().partitioning.partition_count() > 1 {
        return Some(Arc::clone(plan));
    }

    if is_transparent_single_partition_wrapper(plan.name()) {
        let children = plan.children();
        if children.len() == 1 {
            return extract_parallel_source_plan(children[0]);
        }
    }

    if plan.name() != "SortPreservingMergeExec" {
        return None;
    }
    let children = plan.children();
    if children.len() != 1 {
        return None;
    }
    let inner = Arc::clone(children[0]);
    if inner.properties().partitioning.partition_count() > 1 {
        Some(inner)
    } else {
        None
    }
}

fn is_transparent_single_partition_wrapper(plan_name: &str) -> bool {
    matches!(
        plan_name,
        // These nodes preserve the row set and schema.  They can appear above
        // SortPreservingMergeExec in optimized plans, and bypassing them only
        // drops scheduling or batch-sizing behavior.
        "CoalesceBatchesExec" | "CooperativeExec" | "BufferExec"
    )
}

fn make_ctx_and_register(
    options: &LanceCacheBuildOptions,
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

fn streaming_lance_output_schema(
    kind: EnsemblEntityKind,
    source_schema: &Schema,
    source_type: &str,
) -> Result<SchemaRef> {
    match kind {
        EnsemblEntityKind::Variation => {
            let schema = variation_projected_schema(source_schema, source_type)?;
            let schema = with_lance_field_metadata(&schema);
            // Bundle the 27 AF columns into 3 fullzip List<Utf8> columns for the cache.
            Ok(Arc::new(crate::lance_cache::af_bundle::bundle_schema(
                &schema,
            )))
        }
        _ => Err(DataFusionError::Execution(format!(
            "streaming Lance writer is not configured for {kind:?}"
        ))),
    }
}

fn variation_projected_schema(source_schema: &Schema, source_type: &str) -> Result<Schema> {
    let mut fields = Vec::new();
    let forbidden = VARIATION_FORBIDDEN_COLUMNS
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    for name in VARIATION_REQUIRED_COLUMNS {
        if forbidden.contains(name) {
            continue;
        }
        let (_, field) = source_schema.column_with_name(name).ok_or_else(|| {
            DataFusionError::Execution(format!("variation source batch missing column {name}"))
        })?;
        if *name == "start" || *name == "end" {
            fields.push(Field::new(*name, DataType::UInt32, field.is_nullable()));
        } else {
            fields.push(field.as_ref().clone());
        }
    }

    // Derived warm/cold tier column (0 = warm/common, 1 = cold/rare). Appended
    // here rather than read from the source table.
    fields.push(Field::new("tier", DataType::Int8, false));

    let target_schema = with_cache_source_metadata(&Schema::new(fields), source_type);
    validate_variation_schema(&target_schema)?;
    Ok(target_schema)
}

/// Project a source variation batch to the Lance output schema, derive the
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
    options: &LanceCacheBuildOptions,
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
    let meta = lance_field_metadata();
    let key = Field::new("key", DataType::UInt64, false).with_metadata(meta.clone());
    let sift = Field::new("sift", DataType::Binary, false).with_metadata(meta.clone());
    let poly = Field::new("poly", DataType::Binary, false).with_metadata(meta);
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

fn context_index_kind(kind: EnsemblEntityKind) -> LanceIndexKind {
    match kind {
        EnsemblEntityKind::Exon | EnsemblEntityKind::Translation => LanceIndexKind::TranscriptId,
        EnsemblEntityKind::Variation
        | EnsemblEntityKind::Transcript
        | EnsemblEntityKind::RegulatoryFeature
        | EnsemblEntityKind::MotifFeature => LanceIndexKind::Start,
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
    use lance::index::DatasetIndexExt;

    #[test]
    fn lance_entity_dirs_use_lance_suffix() {
        assert_eq!(lance_entity_dir_name("variation"), "variation.lance");
        assert_eq!(
            lance_entity_dir_name("translation_core"),
            "translation_core.lance"
        );
    }

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

    #[tokio::test]
    async fn tiered_variation_write_clusters_warm_rows_first_with_indexes() {
        use datafusion::arrow::array::Int8Array as ArrowInt8Array;
        use datafusion::datasource::MemTable;
        use futures::TryStreamExt;

        // Synthetic source containing every column the variation transform reads.
        let mut fields = Vec::new();
        let mut columns = Vec::<ArrayRef>::new();
        let starts = [300_i64, 100, 100, 200];
        let n = starts.len();
        for name in VARIATION_REQUIRED_COLUMNS {
            match *name {
                "start" => {
                    fields.push(Field::new("start", DataType::Int64, false));
                    columns.push(Arc::new(Int64Array::from(starts.to_vec())) as ArrayRef);
                }
                "end" => {
                    fields.push(Field::new("end", DataType::Int64, false));
                    columns.push(Arc::new(Int64Array::from(starts.to_vec())) as ArrayRef);
                }
                "failed" => {
                    fields.push(Field::new("failed", DataType::Int8, true));
                    columns.push(Arc::new(ArrowInt8Array::from(vec![Some(0); n])) as ArrayRef);
                }
                "chrom" => {
                    fields.push(Field::new("chrom", DataType::Utf8, false));
                    columns.push(Arc::new(StringArray::from(vec!["chr1"; n])) as ArrayRef);
                }
                other => {
                    fields.push(Field::new(other, DataType::Utf8, true));
                    columns.push(Arc::new(StringArray::from(vec![Some("v"); n])) as ArrayRef);
                }
            }
        }
        let schema = Arc::new(Schema::new(fields));
        let batch = RecordBatch::try_new(Arc::clone(&schema), columns).unwrap();

        let ctx = SessionContext::new();
        ctx.register_table(
            "var",
            Arc::new(MemTable::try_new(schema, vec![vec![batch]]).unwrap()),
        )
        .unwrap();

        // start=100 is warm (two rows); 200/300 are cold.
        let warm_starts = Arc::new(BTreeSet::from([100_i64]));
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("variation/chr1.lance");

        let rows = write_variation_tiered_with_ctx(
            &ctx,
            "SELECT * FROM var ORDER BY chrom, start",
            &path,
            warm_starts,
            "ensembl".to_string(),
        )
        .await
        .unwrap();
        assert_eq!(rows, 4);

        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        assert_eq!(dataset.count_rows(None).await.unwrap(), 4);

        // Storage order must be warm-first: all tier=0 rows precede all tier=1 rows.
        let mut stream = dataset
            .scan()
            .project(&["tier", "start"])
            .unwrap()
            .try_into_stream()
            .await
            .unwrap();
        let mut tiers = Vec::new();
        let mut starts_seen = Vec::new();
        while let Some(batch) = stream.try_next().await.unwrap() {
            let tier = batch
                .column(batch.schema().index_of("tier").unwrap())
                .as_any()
                .downcast_ref::<ArrowInt8Array>()
                .unwrap();
            let start = batch
                .column(batch.schema().index_of("start").unwrap())
                .as_any()
                .downcast_ref::<datafusion::arrow::array::UInt32Array>()
                .unwrap();
            for row in 0..batch.num_rows() {
                tiers.push(tier.value(row));
                starts_seen.push(start.value(row));
            }
        }
        assert_eq!(tiers, vec![0, 0, 1, 1], "warm rows must be written first");
        // Warm block is the two start=100 rows; cold block is 200 then 300.
        assert_eq!(starts_seen, vec![100, 100, 200, 300]);

        let index_names = dataset
            .load_indices()
            .await
            .unwrap()
            .iter()
            .map(|idx| idx.name.clone())
            .collect::<Vec<_>>();
        assert!(
            index_names
                .iter()
                .any(|name| name == crate::lance_cache::write::START_BTREE_INDEX_NAME),
            "missing start btree index: {index_names:?}"
        );
        assert!(
            index_names
                .iter()
                .any(|name| name == crate::lance_cache::write::TIER_BITMAP_INDEX_NAME),
            "missing tier bitmap index: {index_names:?}"
        );
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

    #[tokio::test]
    async fn multi_partition_detector_ignores_non_merge_plan() {
        let ctx = SessionContext::new();
        let df = ctx.sql("SELECT 1 AS x").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        assert!(extract_parallel_source_plan(&plan).is_none());
    }

    #[test]
    #[allow(deprecated)]
    fn parallel_source_detector_unwraps_transparent_merge_wrapper() {
        let schema = Arc::new(Schema::new(vec![Field::new("x", DataType::Int32, false)]));
        let make_batch = |value| {
            RecordBatch::try_new(
                Arc::clone(&schema),
                vec![Arc::new(Int32Array::from(vec![value])) as ArrayRef],
            )
            .unwrap()
        };
        let source = TestMemoryExec::try_new_exec(
            &[
                vec![make_batch(1)],
                vec![make_batch(2)],
                vec![make_batch(3)],
            ],
            Arc::clone(&schema),
            None,
        )
        .unwrap();
        let sort = LexOrdering::new(vec![PhysicalSortExpr {
            expr: Arc::new(Column::new("x", 0)),
            options: SortOptions::default(),
        }])
        .expect("non-empty sort ordering");
        let merge: Arc<dyn ExecutionPlan> = Arc::new(SortPreservingMergeExec::new(sort, source));
        let wrapped: Arc<dyn ExecutionPlan> = Arc::new(CoalesceBatchesExec::new(merge, 8192));

        let parallel = extract_parallel_source_plan(&wrapped)
            .expect("transparent wrapper around merge should expose source partitions");
        assert_eq!(parallel.properties().partitioning.partition_count(), 3);
    }

    #[test]
    fn variation_write_strategy_streams_single_lance_dataset() {
        assert_eq!(
            lance_write_strategy(EnsemblEntityKind::Variation),
            LanceWriteStrategy::StreamFullDataset
        );
        assert_eq!(
            lance_write_strategy(EnsemblEntityKind::Transcript),
            LanceWriteStrategy::ChunkedAppend {
                chunk_rows: LANCE_CONTEXT_WRITE_CHUNK_ROWS
            }
        );
        // Translation (incl. translation_sift) writes in large chunks so the
        // narrow dataset lands in a few big, well-compressed fragments.
        assert_eq!(
            lance_write_strategy(EnsemblEntityKind::Translation),
            LanceWriteStrategy::ChunkedAppend {
                chunk_rows: LANCE_TRANSLATION_WRITE_CHUNK_ROWS
            }
        );
    }

    #[test]
    fn variation_streaming_uses_ordered_plan_root() {
        assert!(!should_bypass_ordered_plan_root(
            LanceWriteStrategy::StreamFullDataset
        ));
        assert!(should_bypass_ordered_plan_root(
            LanceWriteStrategy::ChunkedAppend {
                chunk_rows: LANCE_CONTEXT_WRITE_CHUNK_ROWS
            }
        ));
    }

    #[tokio::test]
    async fn coalesced_stream_combines_small_batches_to_target_rows() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "start",
            DataType::UInt32,
            false,
        )]));
        let make_batch = |values| {
            RecordBatch::try_new(
                Arc::clone(&schema),
                vec![Arc::new(UInt32Array::from(values)) as ArrayRef],
            )
            .unwrap()
        };
        let stream = Box::pin(RecordBatchStreamAdapter::new(
            Arc::clone(&schema),
            futures::stream::iter(vec![
                Ok(make_batch(vec![1, 2])),
                Ok(make_batch(vec![3])),
                Ok(make_batch(vec![4, 5])),
                Ok(make_batch(vec![6])),
            ]),
        ));

        let mut coalesced = coalesce_record_batch_stream(stream, Arc::clone(&schema), 4);
        let first = coalesced.next().await.unwrap().unwrap();
        let second = coalesced.next().await.unwrap().unwrap();

        assert_eq!(first.num_rows(), 4);
        assert_eq!(second.num_rows(), 2);
        assert!(coalesced.next().await.is_none());
        let first_values = first
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let second_values = second
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        assert_eq!(
            (0..first_values.len())
                .map(|idx| first_values.value(idx))
                .collect::<Vec<_>>(),
            [1, 2, 3, 4]
        );
        assert_eq!(
            (0..second_values.len())
                .map(|idx| second_values.value(idx))
                .collect::<Vec<_>>(),
            [5, 6]
        );
    }

    #[tokio::test]
    async fn coalesced_stream_tolerates_poll_after_none() {
        // The Lance stream writer can poll the coalesced output stream again
        // after it has already returned `None`. The unfold backing the stream
        // must tolerate that (fused) rather than panic with
        // "Unfold must not be polled after it returned `Poll::Ready(None)`".
        let schema = Arc::new(Schema::new(vec![Field::new(
            "start",
            DataType::UInt32,
            false,
        )]));
        let batch = RecordBatch::try_new(
            Arc::clone(&schema),
            vec![Arc::new(UInt32Array::from(vec![1u32, 2, 3])) as ArrayRef],
        )
        .unwrap();
        let stream = Box::pin(RecordBatchStreamAdapter::new(
            Arc::clone(&schema),
            futures::stream::iter(vec![Ok(batch)]),
        ));
        let mut coalesced = coalesce_record_batch_stream(stream, Arc::clone(&schema), 4);
        assert_eq!(coalesced.next().await.unwrap().unwrap().num_rows(), 3);
        assert!(coalesced.next().await.is_none());
        // Extra polls after exhaustion must keep yielding None, not panic.
        assert!(coalesced.next().await.is_none());
        assert!(coalesced.next().await.is_none());
    }

    #[tokio::test]
    async fn streaming_writer_writes_lance_dataset_and_start_btree() {
        let dir = tempfile::tempdir().unwrap();
        let source_schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ]));
        let first = RecordBatch::try_new(
            Arc::clone(&source_schema),
            vec![
                Arc::new(StringArray::from(vec!["chr1"])) as ArrayRef,
                Arc::new(UInt32Array::from(vec![10])) as ArrayRef,
                Arc::new(UInt32Array::from(vec![10])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/T"])) as ArrayRef,
                Arc::new(Int8Array::from(vec![Some(0)])) as ArrayRef,
                Arc::new(StringArray::from(vec!["rs1"])) as ArrayRef,
            ],
        )
        .unwrap();
        let second = RecordBatch::try_new(
            Arc::clone(&source_schema),
            vec![
                Arc::new(StringArray::from(vec!["chr1"])) as ArrayRef,
                Arc::new(UInt32Array::from(vec![20])) as ArrayRef,
                Arc::new(UInt32Array::from(vec![20])) as ArrayRef,
                Arc::new(StringArray::from(vec!["G/C"])) as ArrayRef,
                Arc::new(Int8Array::from(vec![Some(0)])) as ArrayRef,
                Arc::new(StringArray::from(vec!["rs2"])) as ArrayRef,
            ],
        )
        .unwrap();
        let stream = Box::pin(RecordBatchStreamAdapter::new(
            Arc::clone(&source_schema),
            futures::stream::iter(vec![Ok(first), Ok(second)]),
        ));
        let output_schema = Arc::new(with_lance_field_metadata(source_schema.as_ref()));
        let dataset_path = dir.path().join("variation/chr1.lance");

        let rows = write_stream_to_lance_streaming(
            stream,
            &dataset_path,
            LanceIndexKind::Start,
            output_schema,
            LanceFileVersion::V2_1,
            Ok,
        )
        .await
        .unwrap();

        assert_eq!(rows, 2);
        let dataset = lance::Dataset::open(dataset_path.to_string_lossy().as_ref())
            .await
            .unwrap();
        assert_eq!(dataset.count_rows(None).await.unwrap(), 2);
        assert!(
            dataset
                .load_indices()
                .await
                .unwrap()
                .iter()
                .any(|idx| idx.name == crate::lance_cache::write::START_BTREE_INDEX_NAME)
        );
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

        // Schema is key/sift/poly with no FullZip (miniblock default).
        let schema = compact_translation_sift_position_schema("ensembl");
        assert_eq!(schema.field(0).name(), "key");
        assert_eq!(schema.field(0).data_type(), &DataType::UInt64);
        assert_eq!(schema.field(1).name(), "sift");
        assert_eq!(schema.field(2).name(), "poly");
        assert_ne!(
            schema
                .field(1)
                .metadata()
                .get("lance-encoding:structural-encoding")
                .map(String::as_str),
            Some("fullzip")
        );
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
