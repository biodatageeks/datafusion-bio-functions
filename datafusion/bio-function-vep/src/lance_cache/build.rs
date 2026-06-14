use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::ArrayRef;
use datafusion::arrow::compute::cast;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_ensembl_cache::{
    CacheSourceType as BioFormatsCacheSourceType, EnsemblCacheOptions, EnsemblCacheTableProvider,
    EnsemblEntityKind, VEP_CACHE_REGION_SIZE_BP, build_export_query, translation_core_schema,
    translation_sift_schema,
};
use futures::StreamExt;
use log::info;

use crate::cache_builder::EntityStats;
use crate::lance_cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest, dataset_dir_name,
};
use crate::lance_cache::schema::{
    VARIATION_FORBIDDEN_COLUMNS, VARIATION_REQUIRED_COLUMNS, validate_variation_schema,
    with_cache_source_metadata,
};
use crate::lance_cache::write::{LanceIndexKind, write_record_batches_to_lance};

const DEFAULT_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

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
        let batches =
            collect_query_batches(options, EnsemblEntityKind::Variation, "var", &query).await?;
        let batches = transform_variation_batches(batches, options.cache_source_type.as_str())?;
        let rows = batch_rows(&batches);
        if rows == 0 {
            continue;
        }
        let dataset_name = dataset_dir_name(&chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        write_record_batches_to_lance(&dataset_path, batches, LanceIndexKind::Start).await?;
        manifest_entries.push(ChromDatasetEntry::new(chrom, dataset_name, rows));
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
    let chroms = discover_chroms(options, kind, table_name).await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    let index_kind = context_index_kind(kind);

    for chrom in chroms {
        let query = build_export_query(kind, table_name, Some(&chrom), None);
        let batches = collect_query_batches(options, kind, table_name, &query).await?;
        let batches =
            attach_schema_metadata_to_batches(batches, options.cache_source_type.as_str())?;
        let rows = batch_rows(&batches);
        if rows == 0 {
            continue;
        }
        let dataset_name = dataset_dir_name(&chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        write_record_batches_to_lance(&dataset_path, batches, index_kind).await?;
        manifest_entries.push(ChromDatasetEntry::new(chrom, dataset_name, rows));
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
        let batches =
            collect_query_batches(options, EnsemblEntityKind::Translation, "tl", &query).await?;
        let batches = drop_row_number_batches(batches)?;
        if batch_rows(&batches) == 0 {
            continue;
        }
        let projected = project_batches_to_schema(batches, Arc::clone(&target_schema))?;
        let projected =
            attach_schema_metadata_to_batches(projected, options.cache_source_type.as_str())?;
        let rows = batch_rows(&projected);
        if rows == 0 {
            continue;
        }
        let dataset_name = dataset_dir_name(&chrom);
        let dataset_path = entity_dir.join(&dataset_name);
        write_record_batches_to_lance(&dataset_path, projected, LanceIndexKind::TranscriptId)
            .await?;
        manifest_entries.push(ChromDatasetEntry::new(chrom, dataset_name, rows));
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
    let ctx = make_ctx_and_register(options, kind, table_name)?;
    let table = ctx.table(table_name).await?;
    Ok(chroms_from_schema(table.schema().inner()))
}

async fn collect_query_batches(
    options: &LanceCacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
    query: &str,
) -> Result<Vec<RecordBatch>> {
    let ctx = make_ctx_and_register(options, kind, table_name)?;
    let df = ctx.sql(query).await?;
    let mut stream = df.execute_stream().await?;
    let mut batches = Vec::new();
    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() > 0 {
            batches.push(batch);
        }
    }
    Ok(batches)
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

fn transform_variation_batches(
    batches: Vec<RecordBatch>,
    source_type: &str,
) -> Result<Vec<RecordBatch>> {
    let mut out = Vec::with_capacity(batches.len());
    for batch in batches {
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
        out.push(
            RecordBatch::try_new(Arc::new(target_schema), columns).map_err(|err| {
                DataFusionError::Execution(format!("failed to build Lance variation batch: {err}"))
            })?,
        );
    }
    Ok(out)
}

fn attach_schema_metadata_to_batches(
    batches: Vec<RecordBatch>,
    source_type: &str,
) -> Result<Vec<RecordBatch>> {
    let mut out = Vec::with_capacity(batches.len());
    for batch in batches {
        let schema = with_cache_source_metadata(batch.schema().as_ref(), source_type);
        out.push(
            RecordBatch::try_new(Arc::new(schema), batch.columns().to_vec()).map_err(|err| {
                DataFusionError::Execution(format!("failed to attach Lance schema metadata: {err}"))
            })?,
        );
    }
    Ok(out)
}

fn drop_row_number_batches(batches: Vec<RecordBatch>) -> Result<Vec<RecordBatch>> {
    let mut out = Vec::with_capacity(batches.len());
    for batch in batches {
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
        out.push(
            RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).map_err(|err| {
                DataFusionError::Execution(format!("failed to drop translation row number: {err}"))
            })?,
        );
    }
    Ok(out)
}

fn project_batches_to_schema(
    batches: Vec<RecordBatch>,
    target_schema: SchemaRef,
) -> Result<Vec<RecordBatch>> {
    batches
        .into_iter()
        .map(|batch| project_batch_to_schema(batch, Arc::clone(&target_schema)))
        .collect()
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
    schema
        .metadata()
        .get("bio.vep.chromosomes")
        .and_then(|json| serde_json::from_str(json).ok())
        .unwrap_or_else(|| {
            DEFAULT_CHROMS
                .iter()
                .map(|chrom| (*chrom).to_string())
                .collect()
        })
}

fn batch_rows(batches: &[RecordBatch]) -> usize {
    batches.iter().map(RecordBatch::num_rows).sum()
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
}
