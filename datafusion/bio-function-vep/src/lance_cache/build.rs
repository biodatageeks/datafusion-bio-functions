use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::ArrayRef;
use datafusion::arrow::compute::cast;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
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
use crate::lance_cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest, canonical_chrom_label, dataset_dir_name,
};
use crate::lance_cache::schema::{
    VARIATION_FORBIDDEN_COLUMNS, VARIATION_REQUIRED_COLUMNS, validate_variation_schema,
    with_cache_source_metadata,
};
use crate::lance_cache::write::{
    LanceIndexKind, create_required_index, write_record_batches_to_lance_with_mode_and_version,
};

const DEFAULT_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];
const LANCE_VARIATION_WRITE_CHUNK_ROWS: usize = 1_000_000;
const LANCE_CONTEXT_WRITE_CHUNK_ROWS: usize = 4_096;

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
            fjall_stats: None,
        }]);
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Variation, "var").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();

    for chrom in chroms {
        let query = build_export_query(EnsemblEntityKind::Variation, "var", Some(&chrom), None);
        let manifest_chrom = canonical_chrom_label(&chrom);
        let dataset_name = dataset_dir_name(&manifest_chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        let source_type = options.cache_source_type.as_str().to_string();
        let rows = write_query_stream_to_lance(
            options,
            EnsemblEntityKind::Variation,
            "var",
            &query,
            &dataset_path,
            LanceIndexKind::Start,
            LANCE_VARIATION_WRITE_CHUNK_ROWS,
            LanceFileVersion::V2_1,
            move |batch| transform_variation_batch(batch, &source_type),
        )
        .await?;
        if rows == 0 {
            continue;
        }
        manifest_entries.push(ChromDatasetEntry::new(manifest_chrom, dataset_name, rows));
        files.push((dataset_path.to_string_lossy().to_string(), rows));
    }

    ChromManifest::new(manifest_entries).write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: "variation.lance".to_string(),
        parquet_files: files,
        fjall_stats: None,
    }])
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
            fjall_stats: None,
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
        let source_type = options.cache_source_type.as_str().to_string();
        let rows = write_query_stream_to_lance(
            options,
            kind,
            table_name,
            &query,
            &dataset_path,
            index_kind,
            LANCE_CONTEXT_WRITE_CHUNK_ROWS,
            LanceFileVersion::V2_2,
            move |batch| attach_schema_metadata_to_batch(batch, &source_type),
        )
        .await?;
        if rows == 0 {
            continue;
        }
        manifest_entries.push(ChromDatasetEntry::new(manifest_chrom, dataset_name, rows));
        files.push((dataset_path.to_string_lossy().to_string(), rows));
    }

    ChromManifest::new(manifest_entries).write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: lance_entity_dir_name(entity),
        parquet_files: files,
        fjall_stats: None,
    }])
}

async fn build_lance_translation(options: &LanceCacheBuildOptions) -> Result<Vec<EntityStats>> {
    let core = build_lance_translation_split(
        options,
        "translation_core",
        translation_core_schema(false, options.cache_source_type),
    )
    .await?;
    let sift = build_lance_translation_split(
        options,
        "translation_sift",
        translation_sift_schema(false, options.cache_source_type),
    )
    .await?;
    Ok(vec![core, sift])
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
            fjall_stats: None,
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
            LANCE_CONTEXT_WRITE_CHUNK_ROWS,
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

    ChromManifest::new(manifest_entries).write_to_entity_dir(&entity_dir)?;
    Ok(EntityStats {
        entity: lance_entity_dir_name(entity),
        parquet_files: files,
        fjall_stats: None,
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
    Ok((chroms_from_schema(&schema), schema))
}

async fn write_query_stream_to_lance<F>(
    options: &LanceCacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
    query: &str,
    dataset_path: &Path,
    index_kind: LanceIndexKind,
    chunk_rows: usize,
    data_storage_version: LanceFileVersion,
    transform: F,
) -> Result<usize>
where
    F: FnMut(RecordBatch) -> Result<RecordBatch> + Clone + Send + Sync + 'static,
{
    let ctx = make_ctx_and_register(options, kind, table_name)?;
    let df = ctx.sql(query).await?;
    let plan = df.create_physical_plan().await?;
    if let Some(inner) = extract_multi_partition_inner(&plan) {
        let partition_count = inner.properties().partitioning.partition_count();
        info!(
            "Lance build {} {} using {partition_count} physical partitions",
            table_name,
            dataset_path.display()
        );
        return write_partitioned_plan_to_lance(
            &ctx,
            inner,
            dataset_path,
            index_kind,
            chunk_rows,
            data_storage_version,
            transform,
        )
        .await;
    }

    let stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;
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

fn extract_multi_partition_inner(plan: &Arc<dyn ExecutionPlan>) -> Option<Arc<dyn ExecutionPlan>> {
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

fn transform_variation_batch(batch: RecordBatch, source_type: &str) -> Result<RecordBatch> {
    let schema = batch.schema();
    let mut fields = Vec::new();
    let mut columns = Vec::new();
    let forbidden = VARIATION_FORBIDDEN_COLUMNS
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    for name in VARIATION_REQUIRED_COLUMNS {
        if forbidden.contains(name) {
            continue;
        }
        let (index, field) = schema.column_with_name(name).ok_or_else(|| {
            DataFusionError::Execution(format!("variation source batch missing column {name}"))
        })?;
        if *name == "start" || *name == "end" {
            fields.push(Field::new(*name, DataType::UInt32, field.is_nullable()));
            columns.push(
                cast(batch.column(index).as_ref(), &DataType::UInt32).map_err(|err| {
                    DataFusionError::Execution(format!(
                        "failed to cast variation {name} to UInt32: {err}"
                    ))
                })?,
            );
        } else {
            fields.push(field.as_ref().clone());
            columns.push(batch.column(index).clone());
        }
    }

    let target_schema = with_cache_source_metadata(&Schema::new(fields), source_type);
    validate_variation_schema(&target_schema)?;
    RecordBatch::try_new(Arc::new(target_schema), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to build Lance variation batch: {err}"))
    })
}

#[cfg(test)]
fn transform_variation_batches(
    batches: Vec<RecordBatch>,
    source_type: &str,
) -> Result<Vec<RecordBatch>> {
    batches
        .into_iter()
        .map(|batch| transform_variation_batch(batch, source_type))
        .collect()
}

fn attach_schema_metadata_to_batch(batch: RecordBatch, source_type: &str) -> Result<RecordBatch> {
    let schema = with_cache_source_metadata(batch.schema().as_ref(), source_type);
    RecordBatch::try_new(Arc::new(schema), batch.columns().to_vec()).map_err(|err| {
        DataFusionError::Execution(format!("failed to attach Lance schema metadata: {err}"))
    })
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
    use datafusion::arrow::array::{ArrayRef, Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field};

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

        let transformed = transform_variation_batches(vec![batch], "ensembl").unwrap();
        let schema = transformed[0].schema();

        assert_eq!(
            schema.field_with_name("start").unwrap().data_type(),
            &DataType::UInt32
        );
        assert_eq!(
            schema.field_with_name("end").unwrap().data_type(),
            &DataType::UInt32
        );
        assert!(schema.field_with_name("position_key").is_err());
        assert_eq!(
            schema
                .metadata()
                .get(crate::cache_source::CACHE_SOURCE_METADATA_KEY),
            Some(&"ensembl".to_string())
        );
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

    #[tokio::test]
    async fn multi_partition_detector_ignores_non_merge_plan() {
        let ctx = SessionContext::new();
        let df = ctx.sql("SELECT 1 AS x").await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        assert!(extract_multi_partition_inner(&plan).is_none());
    }
}
