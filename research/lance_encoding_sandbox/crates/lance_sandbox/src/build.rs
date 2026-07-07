use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use arrow_array::{
    Array, ArrayRef, Int64Array, RecordBatch, RecordBatchIterator, StructArray, UInt8Array,
    UInt32Array, UInt64Array, new_null_array,
};
use arrow_schema::{DataType, Field, Fields, Schema, SchemaRef};
use datafusion_bio_function_vep::warm_cache::reader::{
    WARM_RUNTIME_COLUMNS, projection_for_existing_roots,
};
use lance::dataset::{WriteMode, WriteParams};
use lance::index::DatasetIndexExt;
use lance_arrow::RecordBatchExt;
use lance_index::{
    IndexType,
    scalar::{BuiltinIndexType, ScalarIndexParams},
};
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ParquetRecordBatchReaderBuilder};
use serde::Serialize;

use crate::config::{KeyDataType, KeyMode, SandboxConfig};

const TIER_FIELD: &str = "tier";
const WARM_TIER: u8 = 0;
const COLD_TIER: u8 = 1;

#[derive(Debug, Serialize)]
pub struct BuildReport {
    pub dataset_path: PathBuf,
    pub warm_path: PathBuf,
    pub cold_paths: Vec<PathBuf>,
    pub warm_rows: usize,
    pub cold_rows: usize,
    pub elapsed_seconds: f64,
}

pub async fn build_dataset(config: &SandboxConfig) -> Result<BuildReport> {
    let started = Instant::now();
    let output = config.dataset_path();
    if output.exists() {
        if config.dataset.overwrite {
            fs::remove_dir_all(&output)
                .with_context(|| format!("failed to remove '{}'", output.display()))?;
        } else {
            bail!(
                "dataset '{}' already exists; set dataset.overwrite=true",
                output.display()
            );
        }
    }
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create '{}'", parent.display()))?;
    }

    let input_paths = source_variation_paths(config)?;
    let warm_path = input_paths[0].clone();
    let cold_paths = input_paths[1..].to_vec();

    let projection = source_projection_columns(config);
    let warm_raw = read_parquet_batches(&warm_path, config.dataset.batch_size, &projection)?;
    let mut cold_raw = Vec::new();
    for path in &cold_paths {
        cold_raw.extend(read_parquet_batches(
            path,
            config.dataset.batch_size,
            &projection,
        )?);
    }

    let warm_rows = row_count(&warm_raw);
    let cold_rows = row_count(&cold_raw);
    let warm_batches = transform_batches(config, warm_raw, WARM_TIER)?;
    let cold_batches = transform_batches(config, cold_raw, COLD_TIER)?;
    let all_batches = warm_batches
        .iter()
        .chain(cold_batches.iter())
        .cloned()
        .collect::<Vec<_>>();
    let schema = encoded_schema(config, merged_schema(&all_batches)?)?;

    let warm_batches = align_batches(warm_batches, schema.clone())?;
    let cold_batches = align_batches(cold_batches, schema.clone())?;
    write_tier_batches(
        &output,
        schema.clone(),
        warm_batches,
        WriteMode::Overwrite,
        config,
    )
    .await?;
    write_tier_batches(&output, schema, cold_batches, WriteMode::Append, config).await?;
    create_indexes(config).await?;

    Ok(BuildReport {
        dataset_path: output,
        warm_path,
        cold_paths,
        warm_rows,
        cold_rows,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

fn source_projection_columns(config: &SandboxConfig) -> Vec<String> {
    let mut columns = WARM_RUNTIME_COLUMNS
        .iter()
        .copied()
        .filter(|name| !(config.key.mode == KeyMode::Position && *name == "position_key"))
        .map(str::to_string)
        .collect::<Vec<_>>();
    if !columns.iter().any(|name| name == "start") {
        columns.push("start".to_string());
    }
    columns
}

fn transform_batches(
    config: &SandboxConfig,
    batches: Vec<RecordBatch>,
    tier: u8,
) -> Result<Vec<RecordBatch>> {
    batches
        .into_iter()
        .map(|batch| {
            let batch = add_tier(batch, tier)?;
            let batch = add_or_keep_key(config, batch)?;
            pack_structs(config, batch)
        })
        .collect()
}

fn add_tier(batch: RecordBatch, tier: u8) -> Result<RecordBatch> {
    let values = UInt8Array::from(vec![tier; batch.num_rows()]);
    Ok(batch.try_with_column(
        Field::new(TIER_FIELD, DataType::UInt8, false),
        Arc::new(values),
    )?)
}

fn add_or_keep_key(config: &SandboxConfig, batch: RecordBatch) -> Result<RecordBatch> {
    match config.key.mode {
        KeyMode::Position => {
            let start_idx = batch.schema().index_of("start")?;
            let start = batch
                .column(start_idx)
                .as_any()
                .downcast_ref::<Int64Array>()
                .context("start must be Int64")?;
            let mut values = Vec::with_capacity(start.len());
            for row in 0..start.len() {
                if start.is_null(row) {
                    bail!("cannot derive non-null position from null start at row {row}");
                }
                let value = start.value(row);
                if value < 0 || value > u32::MAX as i64 {
                    bail!("start value {value} cannot fit UInt32 position");
                }
                values.push(value as u32);
            }
            let position = Arc::new(UInt32Array::from(values)) as ArrayRef;
            let batch = if batch.schema().index_of(config.key.column.as_str()).is_ok() {
                batch
            } else {
                batch.try_with_column(
                    Field::new(config.key.column.as_str(), DataType::UInt32, false),
                    position,
                )?
            };
            drop_column(batch, "position_key")
        }
        KeyMode::PositionKey => {
            if batch.schema().index_of(config.key.column.as_str()).is_err() {
                bail!(
                    "position_key compatibility mode requires source column '{}'",
                    config.key.column
                );
            }
            Ok(batch)
        }
    }
}

fn pack_structs(config: &SandboxConfig, batch: RecordBatch) -> Result<RecordBatch> {
    let child_to_struct = config.child_to_struct();
    if child_to_struct.is_empty() {
        return Ok(batch);
    }

    let mut emitted_structs = HashSet::new();
    let mut output_fields = Vec::new();
    let mut output_columns = Vec::new();
    let schema = batch.schema();

    for field in schema.fields() {
        let name = field.name();
        if let Some(struct_name) = child_to_struct.get(name) {
            if emitted_structs.insert(struct_name.clone()) {
                let group = &config.structs[struct_name];
                let mut child_fields = Vec::new();
                let mut child_columns = Vec::new();
                for child in &group.fields {
                    if let Ok(idx) = schema.index_of(child) {
                        child_fields.push(schema.field(idx).as_ref().clone().with_nullable(true));
                        child_columns.push(batch.column(idx).clone());
                    }
                }
                if child_fields.is_empty() {
                    continue;
                }
                let fields = Fields::from(child_fields);
                let struct_array = StructArray::try_new(fields.clone(), child_columns, None)?;
                let metadata = packed_struct_metadata(config, struct_name, group.packed_metadata);
                output_fields.push(
                    Field::new(struct_name, DataType::Struct(fields), true).with_metadata(metadata),
                );
                output_columns.push(Arc::new(struct_array) as ArrayRef);
            }
            continue;
        }
        output_fields.push(field.as_ref().clone());
        output_columns.push(batch.column(schema.index_of(name)?).clone());
    }

    Ok(RecordBatch::try_new(
        Arc::new(Schema::new(output_fields)),
        output_columns,
    )?)
}

fn packed_struct_metadata(
    config: &SandboxConfig,
    struct_name: &str,
    include_structural: bool,
) -> HashMap<String, String> {
    let mut metadata = HashMap::from([("lance-encoding:packed".to_string(), "true".to_string())]);
    if include_structural {
        metadata.insert(
            "lance-encoding:structural-encoding".to_string(),
            config.field_encoding(struct_name).structural,
        );
    }
    metadata
}

fn drop_column(batch: RecordBatch, name: &str) -> Result<RecordBatch> {
    let schema = batch.schema();
    let mut fields = Vec::new();
    let mut columns = Vec::new();
    for (idx, field) in schema.fields().iter().enumerate() {
        if field.name() == name {
            continue;
        }
        fields.push(field.as_ref().clone());
        columns.push(batch.column(idx).clone());
    }
    Ok(RecordBatch::try_new(
        Arc::new(Schema::new(fields)),
        columns,
    )?)
}

fn merged_schema(batches: &[RecordBatch]) -> Result<SchemaRef> {
    let mut fields = Vec::<Field>::new();
    let mut present = BTreeMap::<String, usize>::new();
    for batch in batches {
        for field in batch.schema().fields() {
            if !present.contains_key(field.name()) {
                fields.push(field.as_ref().clone());
            }
            *present.entry(field.name().clone()).or_default() += 1;
        }
    }
    for field in &mut fields {
        if present.get(field.name()).copied().unwrap_or_default() < batches.len()
            && !field.is_nullable()
        {
            *field = field.clone().with_nullable(true);
        }
    }
    Ok(Arc::new(Schema::new(fields)))
}

fn encoded_schema(config: &SandboxConfig, schema: SchemaRef) -> Result<SchemaRef> {
    let fields = schema
        .fields()
        .iter()
        .map(|field| encode_field(config, field.as_ref()))
        .collect::<Vec<_>>();
    Ok(Arc::new(Schema::new(fields)))
}

fn encode_field(config: &SandboxConfig, field: &Field) -> Field {
    let mut metadata = field.metadata().clone();
    for (key, value) in config.field_encoding(field.name()).to_lance_metadata() {
        metadata.insert(key, value);
    }
    let field = if let DataType::Struct(children) = field.data_type() {
        let children = children
            .iter()
            .map(|child| Arc::new(encode_field(config, child.as_ref())))
            .collect::<Vec<_>>();
        Field::new(
            field.name(),
            DataType::Struct(Fields::from(children)),
            field.is_nullable(),
        )
    } else {
        field.clone()
    };
    field.with_metadata(metadata)
}

fn align_batches(batches: Vec<RecordBatch>, schema: SchemaRef) -> Result<Vec<RecordBatch>> {
    batches
        .into_iter()
        .map(|batch| {
            let mut columns = Vec::with_capacity(schema.fields().len());
            for field in schema.fields() {
                match batch.schema().index_of(field.name()) {
                    Ok(idx) => columns.push(retag_array_for_field(
                        batch.column(idx).clone(),
                        field.as_ref(),
                    )?),
                    Err(_) => columns.push(new_null_array(field.data_type(), batch.num_rows())),
                }
            }
            Ok(RecordBatch::try_new(schema.clone(), columns)?)
        })
        .collect()
}

fn retag_array_for_field(array: ArrayRef, field: &Field) -> Result<ArrayRef> {
    let DataType::Struct(fields) = field.data_type() else {
        return Ok(array);
    };
    let struct_array = array
        .as_any()
        .downcast_ref::<StructArray>()
        .with_context(|| format!("field '{}' must be a StructArray", field.name()))?;
    let child_columns = fields
        .iter()
        .zip(struct_array.columns())
        .map(|(child_field, child_array)| {
            retag_array_for_field(child_array.clone(), child_field.as_ref())
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Arc::new(StructArray::try_new(
        fields.clone(),
        child_columns,
        struct_array.nulls().cloned(),
    )?))
}

async fn write_tier_batches(
    path: &Path,
    schema: SchemaRef,
    batches: Vec<RecordBatch>,
    mode: WriteMode,
    config: &SandboxConfig,
) -> Result<()> {
    if batches.is_empty() {
        return Ok(());
    }
    let reader = RecordBatchIterator::new(batches.into_iter().map(Ok), schema);
    let params = WriteParams {
        mode,
        data_storage_version: Some(config.lance_file_version()?),
        ..Default::default()
    };
    lance::Dataset::write(reader, path.to_string_lossy().as_ref(), Some(params))
        .await
        .with_context(|| format!("failed to write Lance dataset '{}'", path.display()))?;
    Ok(())
}

async fn create_indexes(config: &SandboxConfig) -> Result<()> {
    let path = config.dataset_path();
    let mut dataset = lance::Dataset::open(path.to_string_lossy().as_ref()).await?;
    if config.indexes.position_btree {
        dataset
            .create_index(
                &[config.key.column.as_str()],
                IndexType::BTree,
                Some(format!("{}_btree_idx", config.key.column)),
                &ScalarIndexParams::for_builtin(BuiltinIndexType::BTree),
                true,
            )
            .await?;
    }
    if config.indexes.tier_bitmap {
        dataset
            .create_index(
                &[TIER_FIELD],
                IndexType::Bitmap,
                Some("tier_bitmap_idx".to_string()),
                &ScalarIndexParams::for_builtin(BuiltinIndexType::Bitmap),
                true,
            )
            .await?;
    }
    Ok(())
}

fn read_parquet_batches(
    path: &Path,
    batch_size: usize,
    projection_columns: &[String],
) -> Result<Vec<RecordBatch>> {
    let file = File::open(path).with_context(|| format!("failed to open '{}'", path.display()))?;
    let metadata = ArrowReaderMetadata::load(&file, Default::default())?;
    let mask = projection_for_existing_roots(
        metadata.schema(),
        metadata.parquet_schema(),
        projection_columns,
    );
    let reader = ParquetRecordBatchReaderBuilder::new_with_metadata(file, metadata)
        .with_projection(mask)
        .with_batch_size(batch_size.max(1))
        .build()?;
    reader.map(|batch| Ok(batch?)).collect()
}

fn row_count(batches: &[RecordBatch]) -> usize {
    batches.iter().map(RecordBatch::num_rows).sum()
}

pub fn source_variation_paths(config: &SandboxConfig) -> Result<Vec<PathBuf>> {
    let variation_dir = config.dataset.cache_root.join("variation");
    let warm_path = variation_split_file_for_chrom(&variation_dir, &config.dataset.chrom, "warm")
        .with_context(|| {
        format!(
            "missing warm parquet for {} under '{}'",
            config.dataset.chrom,
            variation_dir.display()
        )
    })?;
    let cold_paths = cold_variation_files_for_chrom(&variation_dir, &config.dataset.chrom);
    if cold_paths.is_empty() {
        bail!(
            "missing cold parquet for {} under '{}'",
            config.dataset.chrom,
            variation_dir.display()
        );
    }
    let mut paths = Vec::with_capacity(1 + cold_paths.len());
    paths.push(warm_path);
    paths.extend(cold_paths);
    Ok(paths)
}

fn variation_split_file_for_chrom(
    variation_dir: &Path,
    chrom: &str,
    tier: &str,
) -> Option<PathBuf> {
    let candidates = if let Some(stripped) = chrom.strip_prefix("chr") {
        vec![
            variation_dir.join(format!("{chrom}_{tier}.parquet")),
            variation_dir.join(format!("{stripped}_{tier}.parquet")),
        ]
    } else {
        vec![
            variation_dir.join(format!("{chrom}_{tier}.parquet")),
            variation_dir.join(format!("chr{chrom}_{tier}.parquet")),
        ]
    };
    candidates.into_iter().find(|path| path.is_file())
}

fn cold_variation_files_for_chrom(variation_dir: &Path, chrom: &str) -> Vec<PathBuf> {
    let mut files = Vec::new();
    if let Some(path) = variation_split_file_for_chrom(variation_dir, chrom, "cold") {
        files.push(path);
    }
    files.sort();
    files
}

pub fn everything_logical_fields(config: &SandboxConfig) -> Vec<String> {
    let mut fields = WARM_RUNTIME_COLUMNS
        .iter()
        .copied()
        .filter(|name| !(config.key.mode == KeyMode::Position && *name == "position_key"))
        .map(str::to_string)
        .collect::<Vec<_>>();
    if config.key.mode == KeyMode::Position && !fields.iter().any(|name| name == &config.key.column)
    {
        fields.insert(0, config.key.column.clone());
    }
    fields
}

pub fn physical_projection_for_everything(
    config: &SandboxConfig,
    schema: &lance_core::datatypes::Schema,
) -> Result<Vec<String>> {
    let available = schema
        .fields
        .iter()
        .map(|field| field.name.clone())
        .collect::<HashSet<_>>();
    let configured_child_to_struct = config.child_to_struct();
    let inferred_child_to_struct = schema_child_to_struct(schema);
    let mut projection = Vec::new();
    let mut emitted = HashSet::new();
    let mut missing = Vec::new();

    for logical in everything_logical_fields(config) {
        if available.contains(&logical) {
            projection.push(logical);
            continue;
        }
        if let Some(struct_name) = configured_child_to_struct.get(&logical) {
            if available.contains(struct_name) {
                if emitted.insert(struct_name.clone()) {
                    projection.push(struct_name.clone());
                }
                continue;
            }
        }
        if let Some(struct_name) = inferred_child_to_struct.get(&logical) {
            if emitted.insert(struct_name.clone()) {
                projection.push(struct_name.clone());
            }
            continue;
        }
        missing.push(logical);
    }
    if !missing.is_empty() {
        bail!("missing required everything fields: {}", missing.join(", "));
    }
    Ok(projection)
}

fn schema_child_to_struct(schema: &lance_core::datatypes::Schema) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let mut ambiguous = HashSet::new();

    for field in &schema.fields {
        if !matches!(field.data_type(), DataType::Struct(_)) {
            continue;
        }
        for child in &field.children {
            if map
                .insert(child.name.clone(), field.name.clone())
                .is_some_and(|existing| existing != field.name)
            {
                ambiguous.insert(child.name.clone());
            }
        }
    }

    for child in ambiguous {
        map.remove(&child);
    }
    map
}

pub fn key_values_from_batch(config: &SandboxConfig, batch: &RecordBatch) -> Result<Vec<u64>> {
    let idx = batch.schema().index_of(config.key.column.as_str())?;
    let array = batch.column(idx);
    match config.key.data_type {
        KeyDataType::Uint32 => {
            let array = array
                .as_any()
                .downcast_ref::<UInt32Array>()
                .context("key column must be UInt32")?;
            Ok((0..array.len())
                .filter(|&row| !array.is_null(row))
                .map(|row| array.value(row) as u64)
                .collect())
        }
        KeyDataType::Uint64 => {
            let array = array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .context("key column must be UInt64")?;
            Ok((0..array.len())
                .filter(|&row| !array.is_null(row))
                .map(|row| array.value(row))
                .collect())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lance_core::datatypes::Schema as LanceSchema;

    fn test_config(lance_version: &str, structs_toml: &str) -> SandboxConfig {
        toml::from_str(&format!(
            r#"
[dataset]
name = "test"
cache_root = "/tmp/cache"
chrom = "chr1"
output_root = "/tmp/out"
lance_version = "{lance_version}"
batch_size = 65536
overwrite = true

[key]
mode = "position"
column = "position"
type = "uint32"

[indexes]
position_btree = true
tier_bitmap = true

[benchmark]
positions_file = "inputs/positions.txt"
vcf_path = "/tmp/input.vcf.gz"
lookup_batch_sizes = [512]
sort_position_keys = true

[inspect]
include_physical_columns = true
include_index_sizes = true

[defaults.encoding]
structural = "miniblock"
compression = "zstd"
compression_level = 3
dict_values_compression = "zstd"
dict_values_compression_level = 3
rle_threshold = 0.95
dict_size_ratio = 0.99
dict_divisor = 1
minichunk_size = 4096

{structs_toml}
"#
        ))
        .unwrap()
    }

    fn lance_schema(fields: Vec<Field>) -> LanceSchema {
        let arrow_schema = Schema::new(fields);
        LanceSchema::try_from(&arrow_schema).unwrap()
    }

    fn flat_everything_schema(config: &SandboxConfig) -> LanceSchema {
        let mut fields = everything_logical_fields(config)
            .into_iter()
            .map(|name| {
                let data_type = match name.as_str() {
                    "position" => DataType::UInt32,
                    "tier" => DataType::UInt8,
                    _ => DataType::Utf8,
                };
                Field::new(name, data_type, true)
            })
            .collect::<Vec<_>>();
        fields.push(Field::new("tier", DataType::UInt8, false));
        lance_schema(fields)
    }

    #[test]
    fn everything_projection_keeps_flat_fields_top_level() {
        let config = test_config("2.1", "");
        let schema = flat_everything_schema(&config);

        let projection = physical_projection_for_everything(&config, &schema).unwrap();

        assert!(projection.contains(&"position".to_string()));
        assert!(projection.contains(&"variation_name".to_string()));
        assert!(projection.contains(&"dbsnp_ids".to_string()));
    }

    #[test]
    fn everything_projection_deduplicates_packed_parent_fields() {
        let config = test_config(
            "2.2",
            r#"
[structs.identity_text]
enabled = true
packed_metadata = true
fields = ["variation_name", "dbsnp_ids"]
"#,
        );
        let mut fields = everything_logical_fields(&config)
            .into_iter()
            .filter(|name| name != "variation_name" && name != "dbsnp_ids")
            .map(|name| {
                let data_type = match name.as_str() {
                    "position" => DataType::UInt32,
                    "tier" => DataType::UInt8,
                    _ => DataType::Utf8,
                };
                Field::new(name, data_type, true)
            })
            .collect::<Vec<_>>();
        fields.push(Field::new("tier", DataType::UInt8, false));
        fields.push(Field::new(
            "identity_text",
            DataType::Struct(Fields::from(vec![
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("dbsnp_ids", DataType::Utf8, true),
            ])),
            true,
        ));
        let schema = lance_schema(fields);

        let projection = physical_projection_for_everything(&config, &schema).unwrap();

        assert!(projection.contains(&"position".to_string()));
        assert_eq!(
            projection
                .iter()
                .filter(|name| name.as_str() == "identity_text")
                .count(),
            1
        );
        assert!(!projection.contains(&"variation_name".to_string()));
        assert!(!projection.contains(&"dbsnp_ids".to_string()));
    }

    #[test]
    fn everything_projection_infers_packed_parent_fields_from_schema() {
        let config = test_config("2.1", "");
        let mut fields = everything_logical_fields(&config)
            .into_iter()
            .filter(|name| name != "variation_name" && name != "dbsnp_ids")
            .map(|name| {
                let data_type = match name.as_str() {
                    "position" => DataType::UInt32,
                    "tier" => DataType::UInt8,
                    _ => DataType::Utf8,
                };
                Field::new(name, data_type, true)
            })
            .collect::<Vec<_>>();
        fields.push(Field::new("tier", DataType::UInt8, false));
        fields.push(Field::new(
            "identity_text",
            DataType::Struct(Fields::from(vec![
                Field::new("variation_name", DataType::Utf8, true),
                Field::new("dbsnp_ids", DataType::Utf8, true),
            ])),
            true,
        ));
        let schema = lance_schema(fields);

        let projection = physical_projection_for_everything(&config, &schema).unwrap();

        assert!(projection.contains(&"position".to_string()));
        assert_eq!(
            projection
                .iter()
                .filter(|name| name.as_str() == "identity_text")
                .count(),
            1
        );
        assert!(!projection.contains(&"variation_name".to_string()));
        assert!(!projection.contains(&"dbsnp_ids".to_string()));
    }

    #[test]
    fn align_batches_retains_encoded_struct_child_metadata() {
        let plain_child = Field::new("variation_name", DataType::Utf8, true);
        let plain_fields = Fields::from(vec![plain_child]);
        let struct_array = StructArray::try_new(
            plain_fields.clone(),
            vec![Arc::new(arrow_array::StringArray::from(vec!["rs1"])) as ArrayRef],
            None,
        )
        .unwrap();
        let input_schema = Arc::new(Schema::new(vec![Field::new(
            "identity_text",
            DataType::Struct(plain_fields),
            true,
        )]));
        let input = RecordBatch::try_new(input_schema, vec![Arc::new(struct_array)]).unwrap();
        let encoded_child = Field::new("variation_name", DataType::Utf8, true).with_metadata(
            HashMap::from([("lance-encoding:compression".to_string(), "zstd".to_string())]),
        );
        let encoded_schema = Arc::new(Schema::new(vec![Field::new(
            "identity_text",
            DataType::Struct(Fields::from(vec![encoded_child])),
            true,
        )]));

        let aligned = align_batches(vec![input], encoded_schema.clone()).unwrap();

        assert_eq!(
            aligned[0].column(0).data_type(),
            encoded_schema.field(0).data_type()
        );
    }

    #[test]
    fn pack_structs_keeps_compression_off_packed_parent() {
        let config = test_config(
            "2.2",
            r#"
[structs.identity_text]
enabled = true
packed_metadata = true
fields = ["variation_name", "dbsnp_ids"]
"#,
        );
        let schema = Arc::new(Schema::new(vec![
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("dbsnp_ids", DataType::Utf8, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(arrow_array::StringArray::from(vec!["rs1"])) as ArrayRef,
                Arc::new(arrow_array::StringArray::from(vec!["rs1"])) as ArrayRef,
            ],
        )
        .unwrap();

        let packed = pack_structs(&config, batch).unwrap();
        let parent = packed.schema().field(0).clone();

        assert_eq!(
            parent.metadata().get("lance-encoding:packed"),
            Some(&"true".to_string())
        );
        assert_eq!(
            parent.metadata().get("lance-encoding:structural-encoding"),
            Some(&"miniblock".to_string())
        );
        assert!(!parent.metadata().contains_key("lance-encoding:compression"));
    }
}
