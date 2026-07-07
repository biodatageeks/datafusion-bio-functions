use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::fs::File;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use lance::Dataset;
use lance::index::DatasetIndexExt;
use lance_file::reader::FileReader as LanceFileReader;
use lance_io::scheduler::{ScanScheduler, SchedulerConfig};
use lance_io::utils::CachedFileSize;
use parquet::file::reader::{FileReader as ParquetFileReader, SerializedFileReader};
use serde::Serialize;

use crate::build::{physical_projection_for_everything, source_variation_paths};
use crate::config::SandboxConfig;

#[derive(Debug, Clone, Serialize)]
pub struct InspectReport {
    pub dataset_path: String,
    pub configured_lance_version: String,
    pub observed_lance_file_versions: Vec<String>,
    pub fragments: usize,
    pub rows: usize,
    pub total_size_bytes: u64,
    pub data_size_bytes: u64,
    pub index_size_bytes: u64,
    pub metadata_size_bytes: u64,
    pub other_size_bytes: u64,
    pub file_count: usize,
    pub requested_everything_fields: Vec<String>,
    pub resolved_dataset_projection: Vec<String>,
    pub logical_fields: Vec<FieldInspectRow>,
    pub physical_columns: Vec<PhysicalColumnInspectRow>,
    pub indexes: Vec<IndexInspectRow>,
}

#[derive(Debug, Clone, Serialize)]
pub struct FieldInspectRow {
    pub path: String,
    pub id: i32,
    pub data_type: String,
    pub nullable: bool,
    pub metadata: BTreeMap<String, String>,
    pub compressed_bytes: u64,
    pub pages: u64,
    pub input_parquet: Option<ParquetInputStats>,
    pub encodings_observed: Vec<ObservedEncoding>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PhysicalColumnInspectRow {
    pub field_path: String,
    pub field_id: Option<i32>,
    pub file_column_index: u32,
    pub fragments: usize,
    pub rows: u64,
    pub compressed_bytes: u64,
    pub pages: u64,
    pub encodings_observed: Vec<ObservedEncoding>,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Serialize)]
pub struct ObservedEncoding {
    pub layout: String,
    pub encoding: String,
    pub compression: String,
}

#[derive(Debug, Clone, Eq, PartialEq, Serialize)]
pub struct ParquetInputStats {
    pub compressed_bytes: u64,
    pub files: usize,
    pub column_chunks: usize,
    pub encodings: Vec<String>,
    pub compression: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct IndexInspectRow {
    pub name: String,
    pub fields: Vec<String>,
    pub dataset_version: u64,
    pub index_version: i32,
    pub size_bytes: Option<u64>,
    pub files: BTreeMap<String, u64>,
    pub statistics_json: Option<String>,
}

#[derive(Default)]
struct FieldAgg {
    field_id: Option<i32>,
    file_column_index: u32,
    fragments: usize,
    rows: u64,
    compressed_bytes: u64,
    pages: u64,
    encodings: BTreeSet<ObservedEncoding>,
}

#[derive(Default)]
struct SizeAgg {
    total_size_bytes: u64,
    data_size_bytes: u64,
    index_size_bytes: u64,
    metadata_size_bytes: u64,
    other_size_bytes: u64,
    file_count: usize,
}

#[derive(Default)]
struct ParquetInputAgg {
    compressed_bytes: u64,
    files: BTreeSet<String>,
    column_chunks: usize,
    encodings: BTreeSet<String>,
    compression: BTreeSet<String>,
}

pub async fn inspect_dataset(config: &SandboxConfig) -> Result<InspectReport> {
    let dataset_path = config.dataset_path();
    let dataset = Dataset::open(dataset_path.to_string_lossy().as_ref())
        .await
        .with_context(|| format!("failed to open Lance dataset '{}'", dataset_path.display()))?;

    let requested_everything_fields = crate::build::everything_logical_fields(config);
    let resolved_dataset_projection = physical_projection_for_everything(config, dataset.schema())?;
    let parquet_inputs = inspect_input_parquet(config)?;
    inspect_opened_dataset(
        &dataset_path,
        &dataset,
        config.dataset.lance_version.clone(),
        requested_everything_fields,
        resolved_dataset_projection,
        parquet_inputs,
        config.inspect.include_physical_columns,
        config.inspect.include_index_sizes,
    )
    .await
}

pub async fn inspect_dataset_path(
    dataset_path: &Path,
    configured_lance_version: String,
    input_parquet_paths: &[PathBuf],
    position_field: Option<&str>,
    position_source_column: Option<&str>,
    include_physical_columns: bool,
    include_index_sizes: bool,
) -> Result<InspectReport> {
    let dataset = Dataset::open(dataset_path.to_string_lossy().as_ref())
        .await
        .with_context(|| format!("failed to open Lance dataset '{}'", dataset_path.display()))?;
    let parquet_inputs = if input_parquet_paths.is_empty() {
        BTreeMap::new()
    } else {
        let mut stats = inspect_parquet_paths(input_parquet_paths)?;
        if let (Some(field), Some(source)) = (position_field, position_source_column) {
            if let Some(source_stats) = stats.get(source).cloned() {
                stats.insert(field.to_string(), source_stats);
            }
        }
        stats
    };
    inspect_opened_dataset(
        dataset_path,
        &dataset,
        configured_lance_version,
        Vec::new(),
        Vec::new(),
        parquet_inputs,
        include_physical_columns,
        include_index_sizes,
    )
    .await
}

async fn inspect_opened_dataset(
    dataset_path: &Path,
    dataset: &Dataset,
    configured_lance_version: String,
    requested_everything_fields: Vec<String>,
    resolved_dataset_projection: Vec<String>,
    parquet_inputs: BTreeMap<String, ParquetInputStats>,
    include_physical_columns: bool,
    include_index_sizes: bool,
) -> Result<InspectReport> {
    let size_agg = dataset_size_breakdown(dataset_path)?;
    let observed_lance_file_versions = observed_file_versions(dataset);
    let physical_aggs = if include_physical_columns {
        inspect_physical_columns(dataset).await?
    } else {
        BTreeMap::new()
    };
    let logical_fields = inspect_logical_fields(dataset, &physical_aggs, &parquet_inputs);
    let physical_columns = physical_aggs
        .into_iter()
        .map(|(field_path, agg)| PhysicalColumnInspectRow {
            field_path,
            field_id: agg.field_id,
            file_column_index: agg.file_column_index,
            fragments: agg.fragments,
            rows: agg.rows,
            compressed_bytes: agg.compressed_bytes,
            pages: agg.pages,
            encodings_observed: agg.encodings.into_iter().collect(),
        })
        .collect();
    let indexes = if include_index_sizes {
        inspect_indexes(dataset).await?
    } else {
        Vec::new()
    };
    Ok(InspectReport {
        dataset_path: dataset_path.display().to_string(),
        configured_lance_version,
        observed_lance_file_versions,
        fragments: dataset.count_fragments(),
        rows: dataset
            .iter_fragments()
            .map(|fragment| fragment.num_rows().unwrap_or_default())
            .sum(),
        total_size_bytes: size_agg.total_size_bytes,
        data_size_bytes: size_agg.data_size_bytes,
        index_size_bytes: size_agg.index_size_bytes,
        metadata_size_bytes: size_agg.metadata_size_bytes,
        other_size_bytes: size_agg.other_size_bytes,
        file_count: size_agg.file_count,
        requested_everything_fields,
        resolved_dataset_projection,
        logical_fields,
        physical_columns,
        indexes,
    })
}

async fn inspect_physical_columns(dataset: &Dataset) -> Result<BTreeMap<String, FieldAgg>> {
    let mut by_field = BTreeMap::<String, FieldAgg>::new();
    for fragment in dataset.iter_fragments() {
        for data_file in &fragment.files {
            let object_store = dataset.object_store(None).await?;
            let data_dir = dataset.data_dir();
            let file_path = data_dir.clone().join(data_file.path.as_str());
            let file_size = match data_file.file_size_bytes.get() {
                Some(size) => u64::from(size),
                None => object_store.size(&file_path).await?,
            };
            let scheduler = ScanScheduler::new(
                object_store.clone().into(),
                SchedulerConfig::new(2 * 1024 * 1024 * 1024),
            );
            let file = scheduler
                .open_file(&file_path, &CachedFileSize::new(file_size))
                .await?;
            let metadata = LanceFileReader::read_all_metadata(&file).await?;
            for (column_index, column_metadata) in metadata.column_metadatas.iter().enumerate() {
                let (field_id, field_path) =
                    field_for_file_column(dataset, data_file, column_index as i32)?;
                let pages = column_metadata.pages.len() as u64;
                let compressed_bytes = column_metadata
                    .pages
                    .iter()
                    .map(|page| page.buffer_sizes.iter().sum::<u64>())
                    .sum::<u64>();
                let encoding = metadata
                    .column_infos
                    .get(column_index)
                    .map(|info| {
                        let column_encoding = format!("{:?}", info.encoding);
                        let page_encoding = info
                            .page_infos
                            .first()
                            .map(|page| format!("{:?}", page.encoding));
                        observed_encoding_summary(
                            &column_encoding,
                            page_encoding.as_deref(),
                            info.is_structural(),
                        )
                    })
                    .unwrap_or_else(|| ObservedEncoding {
                        layout: "unknown".to_string(),
                        encoding: "unknown".to_string(),
                        compression: "unknown".to_string(),
                    });
                let entry = by_field.entry(field_path).or_insert_with(|| FieldAgg {
                    field_id,
                    file_column_index: column_index as u32,
                    ..FieldAgg::default()
                });
                entry.fragments += 1;
                entry.rows += metadata.num_rows;
                entry.compressed_bytes += compressed_bytes;
                entry.pages += pages;
                entry.encodings.insert(encoding);
            }
        }
    }
    Ok(by_field)
}

fn observed_encoding_summary(
    column_encoding: &str,
    page_encoding: Option<&str>,
    is_structural: bool,
) -> ObservedEncoding {
    ObservedEncoding {
        layout: layout_summary(page_encoding, is_structural),
        encoding: encoding_summary(column_encoding, page_encoding),
        compression: compression_summary(column_encoding, page_encoding),
    }
}

fn layout_summary(page_encoding: Option<&str>, is_structural: bool) -> String {
    let Some(page) = page_encoding else {
        return "no_pages".to_string();
    };
    if page.contains("MiniBlockLayout") {
        layout_with_suffix("MiniBlock", page)
    } else if page.contains("ConstantLayout") {
        "Constant".to_string()
    } else if page.contains("FullZipLayout") {
        "FullZip".to_string()
    } else if page.contains("BlobLayout") {
        "Blob".to_string()
    } else if page.contains("Legacy(") || !is_structural {
        "Legacy".to_string()
    } else {
        "Structural".to_string()
    }
}

fn layout_with_suffix(layout: &str, page: &str) -> String {
    if page.contains("has_large_chunk: true") {
        format!("{layout}(large_chunk)")
    } else {
        layout.to_string()
    }
}

fn encoding_summary(column_encoding: &str, page_encoding: Option<&str>) -> String {
    let mut parts = Vec::new();
    push_if_contains(&mut parts, column_encoding, "Values", "Values");
    push_if_contains(&mut parts, column_encoding, "ZoneIndex", "ZoneIndex");
    push_if_contains(&mut parts, column_encoding, "Blob", "Blob");

    if let Some(page) = page_encoding {
        push_if_contains(&mut parts, page, "Dictionary", "Dictionary");
        if page.contains("dictionary: Some") {
            push_unique(&mut parts, "Dictionary");
        }
        push_if_contains(&mut parts, page, "Fsst", "FSST");
        push_if_contains(&mut parts, page, "ByteStreamSplit", "ByteStreamSplit");
        if page.contains("InlineBitpacking") || page.contains("OutOfLineBitpacking") {
            push_unique(&mut parts, "Bitpacking");
        }
        push_if_contains(&mut parts, page, "Rle", "RLE");
        push_if_contains(
            &mut parts,
            page,
            "VariablePackedStruct",
            "VariablePackedStruct",
        );
        push_if_contains(&mut parts, page, "PackedStruct", "PackedStruct");
        push_if_contains(&mut parts, page, "Variable(", "Variable");
        push_if_contains(&mut parts, page, "Constant(", "Constant");
        for bits in flat_bits(page) {
            push_unique(&mut parts, &format!("Flat({bits}-bit)"));
        }
    }

    if parts.is_empty() {
        "unknown".to_string()
    } else {
        parts.join("+")
    }
}

fn compression_summary(column_encoding: &str, page_encoding: Option<&str>) -> String {
    let mut parts = Vec::new();
    let mut combined = column_encoding.to_string();
    if let Some(page) = page_encoding {
        combined.push_str(page);
    }
    if combined.contains("CompressionAlgorithmZstd") {
        push_unique(
            &mut parts,
            &compression_with_level("zstd", "CompressionAlgorithmZstd", &combined),
        );
    }
    if combined.contains("CompressionAlgorithmLz4") {
        push_unique(
            &mut parts,
            &compression_with_level("lz4", "CompressionAlgorithmLz4", &combined),
        );
    }
    if parts.is_empty() {
        "none".to_string()
    } else {
        parts.join("+")
    }
}

fn compression_with_level(name: &str, marker: &str, text: &str) -> String {
    text.find(marker)
        .and_then(|start| text[start..].find("level: Some(").map(|idx| start + idx))
        .and_then(|level_start| {
            let value_start = level_start + "level: Some(".len();
            text[value_start..]
                .find(')')
                .map(|value_end| text[value_start..value_start + value_end].to_string())
        })
        .map(|level| format!("{name}(level={level})"))
        .unwrap_or_else(|| name.to_string())
}

fn flat_bits(text: &str) -> Vec<String> {
    let mut bits = Vec::new();
    let mut remaining = text;
    while let Some(start) = remaining.find("bits_per_value: ") {
        let value_start = start + "bits_per_value: ".len();
        let value = remaining[value_start..]
            .chars()
            .take_while(|ch| ch.is_ascii_digit())
            .collect::<String>();
        if !value.is_empty() && !bits.contains(&value) {
            bits.push(value);
        }
        remaining = &remaining[value_start..];
    }
    bits
}

fn push_if_contains(parts: &mut Vec<String>, haystack: &str, needle: &str, label: &str) {
    if haystack.contains(needle) {
        push_unique(parts, label);
    }
}

fn push_unique(parts: &mut Vec<String>, value: &str) {
    if !parts.iter().any(|part| part == value) {
        parts.push(value.to_string());
    }
}

fn field_for_file_column(
    dataset: &Dataset,
    data_file: &lance::table::format::DataFile,
    column_index: i32,
) -> Result<(Option<i32>, String)> {
    for (idx, candidate_column_index) in data_file.column_indices.iter().enumerate() {
        if *candidate_column_index == column_index {
            let field_id = data_file.fields.get(idx).copied();
            if let Some(field_id) = field_id {
                let field_path = dataset
                    .schema()
                    .field_path(field_id)
                    .unwrap_or_else(|_| format!("field_id_{field_id}"));
                return Ok((Some(field_id), field_path));
            }
        }
    }
    Ok((None, format!("file_column_{column_index}")))
}

fn inspect_logical_fields(
    dataset: &Dataset,
    physical_aggs: &BTreeMap<String, FieldAgg>,
    parquet_inputs: &BTreeMap<String, ParquetInputStats>,
) -> Vec<FieldInspectRow> {
    dataset
        .schema()
        .fields_pre_order()
        .map(|field| {
            let path = dataset
                .schema()
                .field_path(field.id)
                .unwrap_or_else(|_| field.name.clone());
            let agg = physical_aggs.get(&path);
            let input_parquet = parquet_inputs.get(&path).cloned();
            FieldInspectRow {
                path,
                id: field.id,
                data_type: field.data_type().to_string(),
                nullable: field.nullable,
                metadata: field
                    .metadata
                    .iter()
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect(),
                compressed_bytes: agg.map(|a| a.compressed_bytes).unwrap_or_default(),
                pages: agg.map(|a| a.pages).unwrap_or_default(),
                input_parquet,
                encodings_observed: agg
                    .map(|a| a.encodings.iter().cloned().collect())
                    .unwrap_or_default(),
            }
        })
        .collect()
}

fn inspect_input_parquet(config: &SandboxConfig) -> Result<BTreeMap<String, ParquetInputStats>> {
    let input_paths = source_variation_paths(config)?;
    let raw_stats = inspect_parquet_paths(&input_paths)?;
    let mut mapped = BTreeMap::new();
    for (source_column, stats) in raw_stats {
        mapped.insert(source_column.clone(), stats.clone());
        if config.key.mode == crate::config::KeyMode::Position && source_column == "start" {
            mapped.insert(config.key.column.clone(), stats);
        }
    }
    Ok(mapped)
}

fn inspect_parquet_paths(paths: &[PathBuf]) -> Result<BTreeMap<String, ParquetInputStats>> {
    let mut aggs = BTreeMap::<String, ParquetInputAgg>::new();
    for path in paths {
        let file =
            File::open(path).with_context(|| format!("failed to open '{}'", path.display()))?;
        let reader = SerializedFileReader::new(file)
            .with_context(|| format!("failed to read Parquet metadata '{}'", path.display()))?;
        let metadata = reader.metadata();
        let path_label = path.display().to_string();
        for row_group in metadata.row_groups() {
            for column in row_group.columns() {
                let root = parquet_root_column(column.column_path().string().as_str());
                let entry = aggs.entry(root).or_default();
                entry.compressed_bytes += column.compressed_size().max(0) as u64;
                entry.files.insert(path_label.clone());
                entry.column_chunks += 1;
                for encoding in column.encodings() {
                    entry.encodings.insert(encoding.to_string());
                }
                entry
                    .compression
                    .insert(parquet_compression_label(column.compression()));
            }
        }
    }
    Ok(aggs
        .into_iter()
        .map(|(column, agg)| {
            (
                column,
                ParquetInputStats {
                    compressed_bytes: agg.compressed_bytes,
                    files: agg.files.len(),
                    column_chunks: agg.column_chunks,
                    encodings: agg.encodings.into_iter().collect(),
                    compression: agg.compression.into_iter().collect(),
                },
            )
        })
        .collect())
}

fn parquet_root_column(path: &str) -> String {
    path.split('.').next().unwrap_or(path).to_string()
}

fn parquet_compression_label(compression: parquet::basic::Compression) -> String {
    compression
        .to_string()
        .split('(')
        .next()
        .unwrap()
        .to_string()
}

async fn inspect_indexes(dataset: &Dataset) -> Result<Vec<IndexInspectRow>> {
    let indices = dataset.load_indices().await?;
    let mut rows = Vec::new();
    for index in indices.iter() {
        let fields = index
            .fields
            .iter()
            .map(|field_id| {
                dataset
                    .schema()
                    .field_path(*field_id)
                    .unwrap_or_else(|_| format!("field_id_{field_id}"))
            })
            .collect::<Vec<_>>();
        let statistics_json = dataset.index_statistics(&index.name).await.ok();
        rows.push(IndexInspectRow {
            name: index.name.clone(),
            fields,
            dataset_version: index.dataset_version,
            index_version: index.index_version,
            size_bytes: index.total_size_bytes(),
            files: index.file_size_map().into_iter().collect(),
            statistics_json,
        });
    }
    rows.sort_by(|left, right| left.name.cmp(&right.name));
    Ok(rows)
}

fn observed_file_versions(dataset: &Dataset) -> Vec<String> {
    let mut versions = BTreeSet::new();
    for fragment in dataset.iter_fragments() {
        for data_file in &fragment.files {
            versions.insert(format!(
                "{}.{}",
                data_file.file_major_version, data_file.file_minor_version
            ));
        }
    }
    versions.into_iter().collect()
}

fn dataset_size_breakdown(path: &Path) -> Result<SizeAgg> {
    let mut agg = SizeAgg::default();
    walk_sizes(path, path, &mut agg)?;
    Ok(agg)
}

fn walk_sizes(root: &Path, path: &Path, agg: &mut SizeAgg) -> Result<()> {
    for entry in
        fs::read_dir(path).with_context(|| format!("failed to read '{}'", path.display()))?
    {
        let entry = entry?;
        let entry_path = entry.path();
        let metadata = entry.metadata()?;
        if metadata.is_dir() {
            walk_sizes(root, &entry_path, agg)?;
            continue;
        }
        if !metadata.is_file() {
            continue;
        }
        let size = metadata.len();
        agg.file_count += 1;
        agg.total_size_bytes += size;
        match classify_dataset_file(root, &entry_path) {
            DatasetFileClass::Data => agg.data_size_bytes += size,
            DatasetFileClass::Index => agg.index_size_bytes += size,
            DatasetFileClass::Metadata => agg.metadata_size_bytes += size,
            DatasetFileClass::Other => agg.other_size_bytes += size,
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
enum DatasetFileClass {
    Data,
    Index,
    Metadata,
    Other,
}

fn classify_dataset_file(root: &Path, path: &Path) -> DatasetFileClass {
    let rel = path.strip_prefix(root).unwrap_or(path);
    let rel_text = rel.to_string_lossy();
    if rel_text.starts_with("_indices/") || rel_text.contains("/_indices/") {
        DatasetFileClass::Index
    } else if rel_text.starts_with("_versions/")
        || rel_text.starts_with("_transactions/")
        || rel_text.starts_with("_deletions/")
    {
        DatasetFileClass::Metadata
    } else if path.extension().and_then(|ext| ext.to_str()) == Some("lance") {
        DatasetFileClass::Data
    } else {
        DatasetFileClass::Other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn observed_encoding_summary_reports_layout_encoding_and_compression() {
        let column = "ColumnEncoding { column_encoding: Some(Values(())) }";
        let page = "Structural(PageLayout { layout: Some(MiniBlockLayout(MiniBlockLayout { value_compression: Some(CompressiveEncoding { compression: Some(General(General { compression: Some(BufferCompression { scheme: CompressionAlgorithmZstd, level: Some(3) }), values: Some(CompressiveEncoding { compression: Some(ByteStreamSplit(ByteStreamSplit { values: Some(CompressiveEncoding { compression: Some(Flat(Flat { bits_per_value: 32, data: None })) }) })) }) })) }), dictionary: Some(CompressiveEncoding { compression: Some(Variable(Variable { offsets: Some(CompressiveEncoding { compression: Some(Flat(Flat { bits_per_value: 32, data: None })) }), values: None })) }), has_large_chunk: false })) })";

        let summary = observed_encoding_summary(column, Some(page), true);

        assert_eq!(summary.layout, "MiniBlock");
        assert_eq!(
            summary.encoding,
            "Values+Dictionary+ByteStreamSplit+Variable+Flat(32-bit)"
        );
        assert_eq!(summary.compression, "zstd(level=3)");
    }

    #[test]
    fn observed_encoding_summary_handles_missing_pages() {
        let summary = observed_encoding_summary(
            "ColumnEncoding { column_encoding: Some(Values(())) }",
            None,
            false,
        );

        assert_eq!(summary.layout, "no_pages");
        assert_eq!(summary.encoding, "Values");
        assert_eq!(summary.compression, "none");
    }
}
