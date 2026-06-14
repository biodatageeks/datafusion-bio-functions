use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use datafusion::arrow::array::{
    Array, ArrayRef, BooleanArray, UInt8Array, UInt64Array, new_null_array,
};
use datafusion::arrow::compute::{concat_batches, filter_record_batch};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::arrow::record_batch::{RecordBatch, RecordBatchIterator};
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use lance::dataset::scanner::ExecutionSummaryCounts;
use lance::dataset::{WriteMode, WriteParams};
use lance::index::DatasetIndexExt;
use lance_file::version::LanceFileVersion;
use lance_index::{
    IndexType,
    scalar::{BuiltinIndexType, ScalarIndexParams},
};

use crate::warm_cache::chrom_cache::WarmProbe;
use crate::warm_cache::chunk::{WarmChunkContext, WarmChunkProbe};
use crate::warm_cache::cold_parquet::ColdProbeResult;

pub const LANCE_VARIATION_DIR: &str = "variation.lance";
pub const WARM_TIER: u8 = 0;
pub const COLD_TIER: u8 = 1;
pub const DEFAULT_LANCE_COLD_LOOKUP_BATCH_SIZE: usize = 2_000;
pub const DEFAULT_LANCE_WARM_SCAN_BATCH_SIZE: usize = 16_384;
pub const DEFAULT_WARM_LANCE_ROWS_PER_FILE: usize = 500_000;
pub const DEFAULT_WARM_LANCE_ROW_GROUP_ROWS: usize = 262_144;
pub const DEFAULT_COLD_LANCE_ROWS_PER_FILE: usize = 1_000_000;
pub const DEFAULT_LANCE_MINICHUNK_SIZE: usize = 16 * 1024;

const LANCE_MINIBLOCK_ZSTD3_METADATA: &[(&str, &str)] = &[
    ("lance-encoding:structural-encoding", "miniblock"),
    ("lance-encoding:compression", "zstd"),
    ("lance-encoding:compression-level", "3"),
    ("lance-encoding:dict-values-compression", "zstd"),
    ("lance-encoding:dict-values-compression-level", "3"),
    ("lance-encoding:rle-threshold", "0.95"),
    ("lance-encoding:dict-size-ratio", "0.99"),
    ("lance-encoding:dict-divisor", "1"),
];

const LANCE_REQUIRED_RUNTIME_COLUMNS: &[&str] = &[
    "position_key",
    "variant_keys",
    "allele_string",
    "end",
    "failed",
];
const LANCE_POSITION_KEY_COLUMN: &str = "position_key";
const LANCE_POSITION_KEY_INDEX_NAME: &str = "position_key_idx";
const LANCE_TIER_COLUMN: &str = "tier";
const LANCE_TIER_BITMAP_INDEX_NAME: &str = "tier_bitmap_idx";

pub struct LanceVariationLookup {
    dataset: lance::Dataset,
    projection_columns: Vec<String>,
    warm_chunk: Option<WarmChunkContext>,
    cold_chunks: std::collections::HashMap<u64, WarmChunkContext>,
    batch_size: usize,
}

pub fn lance_variation_dataset_path(cache_root: &Path, chrom: &str) -> PathBuf {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    cache_root
        .join(LANCE_VARIATION_DIR)
        .join(format!("chr{bare}.lance"))
}

pub async fn read_lance_variation_schema(path: &Path) -> Result<Schema> {
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to open Lance variation dataset '{}': {err}",
                path.display()
            ))
        })?;
    Ok(dataset.schema().into())
}

pub fn lance_position_filter(tier: u8, position_keys: &[u64]) -> String {
    let mut keys = position_keys.to_vec();
    keys.sort_unstable();
    keys.dedup();
    let values = if keys.is_empty() {
        "NULL".to_string()
    } else {
        keys.iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join(",")
    };
    format!("tier = {tier} AND position_key IN ({values})")
}

impl LanceVariationLookup {
    pub async fn open(
        path: &Path,
        projection_columns: Vec<String>,
        batch_size: usize,
    ) -> Result<Self> {
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!(
                    "failed to open Lance variation dataset '{}': {err}",
                    path.display()
                ))
            })?;
        let projection_columns = lance_projection_columns(&dataset, &projection_columns);

        Ok(Self {
            dataset,
            projection_columns,
            warm_chunk: None,
            cold_chunks: std::collections::HashMap::new(),
            batch_size: batch_size.max(1),
        })
    }

    pub async fn load_warm_tier(&mut self) -> Result<()> {
        let batches = self.scan_filter("tier = 0").await?;
        self.warm_chunk = concat_lance_batches(batches)?
            .map(|batch| WarmChunkContext::try_new(0, batch))
            .transpose()?;
        Ok(())
    }

    pub async fn prefetch_cold_positions<I>(&mut self, position_keys: I) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        let mut keys: Vec<u64> = position_keys.into_iter().collect();
        keys.sort_unstable();
        keys.dedup();
        if keys.is_empty() {
            return Ok(());
        }

        for chunk in keys.chunks(self.batch_size) {
            let filter = lance_position_filter(COLD_TIER, chunk);
            let batches = self.scan_filter(&filter).await?;
            let Some(batch) = concat_lance_batches(batches)? else {
                continue;
            };
            let position_idx = batch.schema().index_of("position_key").map_err(|_| {
                DataFusionError::Execution("Lance cold batch missing position_key".into())
            })?;
            for position_key in chunk {
                let Some(position_batch) =
                    filter_batch_for_position(&batch, position_idx, *position_key)?
                else {
                    continue;
                };
                self.cold_chunks.insert(
                    *position_key,
                    WarmChunkContext::try_new_without_variant_index(0, position_batch)?,
                );
            }
        }

        Ok(())
    }

    pub fn warm_probe_exact<V>(
        &self,
        position_key: u64,
        variant_key: u64,
        mut verify_row: V,
    ) -> Result<WarmProbe>
    where
        V: FnMut(&WarmChunkContext, u32) -> Result<bool>,
    {
        let Some(chunk) = self.warm_chunk.as_ref() else {
            return Ok(WarmProbe::NotCovered);
        };

        Ok(
            match chunk.probe_exact(position_key, variant_key, |row, _| verify_row(chunk, row))? {
                WarmChunkProbe::Exact { row, position_rows } => {
                    WarmProbe::Exact { row, position_rows }
                }
                WarmChunkProbe::PositionCoveredNoExact { position_rows } => {
                    WarmProbe::PositionCoveredNoExact { position_rows }
                }
                WarmChunkProbe::NotCovered => WarmProbe::NotCovered,
            },
        )
    }

    pub fn warm_chunk(&self) -> Option<&WarmChunkContext> {
        self.warm_chunk.as_ref()
    }

    pub fn cached_chunk_for_position(&self, position_key: u64) -> Option<&WarmChunkContext> {
        self.cold_chunks.get(&position_key)
    }

    pub fn cold_probe_position_emit_and_visit<P, E, V>(
        &self,
        position_key: u64,
        mut allele_matches: P,
        mut emit_match: E,
        mut visit_row: V,
    ) -> Result<ColdProbeResult>
    where
        P: FnMut(&WarmChunkContext, u32, &str) -> Result<bool>,
        E: FnMut(&WarmChunkContext, u32) -> Result<()>,
        V: FnMut(&WarmChunkContext, u32, &str) -> Result<()>,
    {
        let Some(chunk) = self.cold_chunks.get(&position_key) else {
            return Ok(ColdProbeResult::NotCovered);
        };
        let rows = chunk.rows_for_position(position_key);
        if rows.is_empty() {
            return Ok(ColdProbeResult::NotCovered);
        }

        let mut matched = false;
        for row in rows {
            let Some(allele_string) = chunk.allele_string(row as usize)? else {
                continue;
            };
            visit_row(chunk, row, allele_string)?;
            if !matched && allele_matches(chunk, row, allele_string)? {
                emit_match(chunk, row)?;
                matched = true;
            }
        }

        if matched {
            Ok(ColdProbeResult::Match)
        } else {
            Ok(ColdProbeResult::PositionCoveredNoExact)
        }
    }

    async fn scan_filter(&self, filter: &str) -> Result<Vec<RecordBatch>> {
        let projection = self
            .projection_columns
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        let profile_enabled = lance_profile_enabled();
        let (scan_kind, filter_keys) = classify_lance_scan_filter(filter);
        let scan_batch_size = match scan_kind {
            LanceScanKind::Warm => lance_warm_scan_batch_size(),
            LanceScanKind::Cold | LanceScanKind::Other => self.batch_size,
        };
        let stats_slot = profile_enabled.then(|| Arc::new(Mutex::new(None::<LanceScanStats>)));
        let mut scanner = self.dataset.scan();
        scanner
            .filter(filter)
            .map_err(|err| DataFusionError::Execution(format!("invalid Lance filter: {err}")))?;
        scanner.project(&projection).map_err(|err| {
            DataFusionError::Execution(format!("invalid Lance projection: {err}"))
        })?;
        scanner.batch_size(scan_batch_size);
        if let Some(stats_slot) = stats_slot.as_ref() {
            let stats_sink = Arc::clone(stats_slot);
            scanner.scan_stats_callback(Arc::new(move |stats| {
                *stats_sink
                    .lock()
                    .expect("Lance scan profile stats lock poisoned") =
                    Some(LanceScanStats::from_summary(stats));
            }));
        }

        let scan_started = Instant::now();
        let mut stream = scanner.try_into_stream().await.map_err(|err| {
            DataFusionError::Execution(format!("failed to scan Lance variation dataset: {err}"))
        })?;
        let mut batches = Vec::new();
        let mut rows = 0usize;
        while let Some(batch) = stream.try_next().await.map_err(|err| {
            DataFusionError::Execution(format!("failed to read Lance variation batch: {err}"))
        })? {
            rows += batch.num_rows();
            batches.push(batch);
        }
        if profile_enabled {
            let stats = stats_slot
                .as_ref()
                .and_then(|slot| {
                    slot.lock()
                        .expect("Lance scan profile stats lock poisoned")
                        .clone()
                })
                .unwrap_or_default();
            eprintln!(
                "{}",
                format_lance_scan_profile_line(
                    scan_kind,
                    projection.len(),
                    scan_batch_size,
                    filter_keys,
                    rows,
                    batches.len(),
                    scan_started.elapsed(),
                    &stats,
                )
            );
        }
        Ok(batches)
    }
}

pub async fn write_merged_lance_variation_dataset(
    path: &Path,
    warm_batches: Vec<RecordBatch>,
    cold_batches: Vec<RecordBatch>,
    warm_fragment_rows: usize,
    cold_fragment_rows: usize,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to create Lance variation parent '{}': {err}",
                parent.display()
            ))
        })?;
    }

    let warm_batches = warm_batches
        .into_iter()
        .map(|batch| append_tier_column(batch, WARM_TIER))
        .collect::<Result<Vec<_>>>()?;
    let cold_batches = cold_batches
        .into_iter()
        .map(|batch| append_tier_column(batch, COLD_TIER))
        .collect::<Result<Vec<_>>>()?;
    let batches = warm_batches
        .iter()
        .chain(cold_batches.iter())
        .cloned()
        .collect::<Vec<_>>();

    if batches.is_empty() {
        return Err(DataFusionError::Execution(
            "cannot write empty Lance variation dataset".into(),
        ));
    }
    let schema = lance_2_1_unpacked_schema(merged_lance_schema(&batches)?);
    let warm_batches = warm_batches
        .into_iter()
        .map(|batch| align_batch_to_schema(batch, schema.clone()))
        .collect::<Result<Vec<_>>>()?;
    let cold_batches = cold_batches
        .into_iter()
        .map(|batch| align_batch_to_schema(batch, schema.clone()))
        .collect::<Result<Vec<_>>>()?;

    let mut wrote = false;
    if !warm_batches.is_empty() {
        write_lance_batches(
            path,
            schema.clone(),
            warm_batches,
            WriteMode::Overwrite,
            warm_fragment_rows.max(1),
            DEFAULT_WARM_LANCE_ROW_GROUP_ROWS
                .min(warm_fragment_rows.max(1))
                .max(1),
        )
        .await?;
        wrote = true;
    }
    if !cold_batches.is_empty() {
        write_lance_batches(
            path,
            schema,
            cold_batches,
            if wrote {
                WriteMode::Append
            } else {
                WriteMode::Overwrite
            },
            DEFAULT_COLD_LANCE_ROWS_PER_FILE.max(cold_fragment_rows.max(1)),
            cold_fragment_rows.max(1),
        )
        .await?;
    }
    create_lance_variation_indices(path).await?;
    Ok(())
}

async fn create_lance_variation_indices(path: &Path) -> Result<()> {
    let mut dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to open Lance variation dataset for indexing '{}': {err}",
                path.display()
            ))
        })?;

    dataset
        .create_index(
            &[LANCE_POSITION_KEY_COLUMN],
            IndexType::BTree,
            Some(LANCE_POSITION_KEY_INDEX_NAME.to_string()),
            &ScalarIndexParams::for_builtin(BuiltinIndexType::BTree),
            true,
        )
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to create Lance position_key index '{}': {err}",
                path.display()
            ))
        })?;
    dataset
        .create_index(
            &[LANCE_TIER_COLUMN],
            IndexType::Bitmap,
            Some(LANCE_TIER_BITMAP_INDEX_NAME.to_string()),
            &ScalarIndexParams::for_builtin(BuiltinIndexType::Bitmap),
            true,
        )
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to create Lance tier bitmap index '{}': {err}",
                path.display()
            ))
        })?;
    Ok(())
}

async fn write_lance_batches(
    path: &Path,
    schema: Arc<Schema>,
    batches: Vec<RecordBatch>,
    mode: WriteMode,
    max_rows_per_file: usize,
    max_rows_per_group: usize,
) -> Result<()> {
    let reader = RecordBatchIterator::new(batches.into_iter().map(Ok), schema);
    let params = WriteParams {
        mode,
        max_rows_per_file,
        max_rows_per_group,
        data_storage_version: Some(LanceFileVersion::V2_1),
        ..Default::default()
    };

    lance::Dataset::write(reader, path.to_string_lossy().as_ref(), Some(params))
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to write Lance variation dataset '{}': {err}",
                path.display()
            ))
        })?;
    Ok(())
}

fn merged_lance_schema(batches: &[RecordBatch]) -> Result<Arc<Schema>> {
    let mut fields = Vec::<Field>::new();
    let mut indexes = HashMap::<String, usize>::new();
    let mut present_counts = Vec::<usize>::new();

    for batch in batches {
        for field in batch.schema().fields() {
            let field = field.as_ref();
            if let Some(&idx) = indexes.get(field.name()) {
                let existing = &fields[idx];
                if existing.data_type() != field.data_type() {
                    return Err(DataFusionError::Execution(format!(
                        "cannot merge Lance variation field '{}' with incompatible types: {:?} vs {:?}",
                        field.name(),
                        existing.data_type(),
                        field.data_type()
                    )));
                }
                if existing.metadata() != field.metadata() {
                    return Err(DataFusionError::Execution(format!(
                        "cannot merge Lance variation field '{}' with incompatible metadata",
                        field.name()
                    )));
                }
                present_counts[idx] += 1;
                if field.is_nullable() && !existing.is_nullable() {
                    fields[idx] = existing.clone().with_nullable(true);
                }
            } else {
                indexes.insert(field.name().clone(), fields.len());
                fields.push(field.clone());
                present_counts.push(1);
            }
        }
    }

    for (idx, field) in fields.iter_mut().enumerate() {
        if present_counts[idx] < batches.len() && !field.is_nullable() {
            *field = field.clone().with_nullable(true);
        }
    }

    Ok(Arc::new(Schema::new(fields)))
}

fn lance_2_1_unpacked_schema(schema: Arc<Schema>) -> Arc<Schema> {
    let minichunk_size = lance_minichunk_size().to_string();
    Arc::new(Schema::new(
        schema
            .fields()
            .iter()
            .map(|field| {
                let mut metadata = field.metadata().clone();
                for (key, value) in LANCE_MINIBLOCK_ZSTD3_METADATA {
                    metadata.insert((*key).to_string(), (*value).to_string());
                }
                metadata.insert(
                    "lance-encoding:minichunk-size".to_string(),
                    minichunk_size.clone(),
                );
                field.as_ref().clone().with_metadata(metadata)
            })
            .collect::<Vec<_>>(),
    ))
}

fn align_batch_to_schema(batch: RecordBatch, schema: Arc<Schema>) -> Result<RecordBatch> {
    let mut columns = Vec::with_capacity(schema.fields().len());
    for field in schema.fields() {
        match batch.schema().index_of(field.name()) {
            Ok(idx) => columns.push(batch.column(idx).clone()),
            Err(_) => columns.push(new_null_array(field.data_type(), batch.num_rows())),
        }
    }
    Ok(RecordBatch::try_new(schema, columns)?)
}

fn append_tier_column(batch: RecordBatch, tier: u8) -> Result<RecordBatch> {
    let mut fields = Vec::with_capacity(batch.num_columns() + 1);
    fields.push(Field::new("tier", DataType::UInt8, false));
    fields.extend(
        batch
            .schema()
            .fields()
            .iter()
            .map(|field| field.as_ref().clone()),
    );

    let mut columns = Vec::with_capacity(batch.num_columns() + 1);
    columns.push(Arc::new(UInt8Array::from(vec![tier; batch.num_rows()])) as ArrayRef);
    columns.extend(batch.columns().iter().cloned());

    Ok(RecordBatch::try_new(
        Arc::new(Schema::new(fields)),
        columns,
    )?)
}

fn lance_projection_columns(dataset: &lance::Dataset, requested_columns: &[String]) -> Vec<String> {
    let schema: Schema = dataset.schema().into();
    let mut columns = Vec::new();
    for name in LANCE_REQUIRED_RUNTIME_COLUMNS {
        push_projection_column(&mut columns, &schema, name);
    }
    for name in requested_columns {
        if name != "tier" {
            push_projection_column(&mut columns, &schema, name);
        }
    }
    columns
}

fn lance_profile_enabled() -> bool {
    std::env::var_os("VEP_LANCE_PROFILE").is_some()
        || std::env::var_os("VEP_KV_PROFILE_DETAILED").is_some()
}

fn lance_warm_scan_batch_size() -> usize {
    std::env::var("VEP_LANCE_WARM_SCAN_BATCH_SIZE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|value| *value > 0)
        .unwrap_or(DEFAULT_LANCE_WARM_SCAN_BATCH_SIZE)
}

fn lance_minichunk_size() -> usize {
    std::env::var("VEP_LANCE_MINICHUNK_SIZE")
        .ok()
        .and_then(|value| parse_lance_byte_size(&value))
        .filter(|value| *value > 0)
        .unwrap_or(DEFAULT_LANCE_MINICHUNK_SIZE)
}

fn parse_lance_byte_size(value: &str) -> Option<usize> {
    let trimmed = value.trim();
    let lower = trimmed.to_ascii_lowercase();
    if let Some(number) = lower.strip_suffix("kb").or_else(|| lower.strip_suffix('k')) {
        number.trim().parse::<usize>().ok()?.checked_mul(1024)
    } else {
        trimmed.parse::<usize>().ok()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LanceScanKind {
    Warm,
    Cold,
    Other,
}

impl LanceScanKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Warm => "warm",
            Self::Cold => "cold",
            Self::Other => "other",
        }
    }
}

#[derive(Debug, Default, Clone)]
struct LanceScanStats {
    iops: usize,
    requests: usize,
    bytes_read: usize,
    indices_loaded: usize,
    parts_loaded: usize,
    index_comparisons: usize,
    fragments_scanned: usize,
    ranges_scanned: usize,
    rows_scanned: usize,
}

impl LanceScanStats {
    fn from_summary(summary: &ExecutionSummaryCounts) -> Self {
        Self {
            iops: summary.iops,
            requests: summary.requests,
            bytes_read: summary.bytes_read,
            indices_loaded: summary.indices_loaded,
            parts_loaded: summary.parts_loaded,
            index_comparisons: summary.index_comparisons,
            fragments_scanned: lance_summary_count(summary, "fragments_scanned"),
            ranges_scanned: lance_summary_count(summary, "ranges_scanned"),
            rows_scanned: lance_summary_count(summary, "rows_scanned"),
        }
    }

    fn index_used(&self) -> bool {
        self.indices_loaded > 0 || self.parts_loaded > 0 || self.index_comparisons > 0
    }
}

fn lance_summary_count(summary: &ExecutionSummaryCounts, name: &str) -> usize {
    summary.all_counts.get(name).copied().unwrap_or_default()
}

fn classify_lance_scan_filter(filter: &str) -> (LanceScanKind, usize) {
    let normalized = filter.split_whitespace().collect::<Vec<_>>().join(" ");
    let keys = lance_filter_in_value_count(filter);
    if normalized == "tier = 0" {
        (LanceScanKind::Warm, keys)
    } else if normalized.starts_with("tier = 1") {
        (LanceScanKind::Cold, keys)
    } else {
        (LanceScanKind::Other, keys)
    }
}

fn lance_filter_in_value_count(filter: &str) -> usize {
    let Some(open_idx) = filter.find("IN (") else {
        return 0;
    };
    let value_start = open_idx + "IN (".len();
    let Some(close_idx) = filter[value_start..].find(')') else {
        return 0;
    };
    let values = filter[value_start..value_start + close_idx].trim();
    if values.is_empty() || values.eq_ignore_ascii_case("NULL") {
        0
    } else {
        values.split(',').count()
    }
}

fn format_lance_scan_profile_line(
    kind: LanceScanKind,
    projected_cols: usize,
    batch_size: usize,
    filter_keys: usize,
    rows: usize,
    batches: usize,
    scan_elapsed: Duration,
    stats: &LanceScanStats,
) -> String {
    format!(
        "[vep-lance-profile] scan kind={} projected_cols={} batch_size={} filter_keys={} rows={} batches={} scan_s={:.3} index_used={} indices_loaded={} parts_loaded={} index_comparisons={} fragments_scanned={} ranges_scanned={} rows_scanned={} bytes_read={} requests={} iops={}",
        kind.as_str(),
        projected_cols,
        batch_size,
        filter_keys,
        rows,
        batches,
        scan_elapsed.as_secs_f64(),
        stats.index_used(),
        stats.indices_loaded,
        stats.parts_loaded,
        stats.index_comparisons,
        stats.fragments_scanned,
        stats.ranges_scanned,
        stats.rows_scanned,
        stats.bytes_read,
        stats.requests,
        stats.iops,
    )
}

fn push_projection_column(columns: &mut Vec<String>, schema: &Schema, name: &str) {
    if schema.index_of(name).is_ok() && !columns.iter().any(|existing| existing == name) {
        columns.push(name.to_string());
    }
}

fn concat_lance_batches(batches: Vec<RecordBatch>) -> Result<Option<RecordBatch>> {
    let Some(schema) = batches.first().map(RecordBatch::schema) else {
        return Ok(None);
    };
    Ok(Some(concat_batches(&schema, &batches)?))
}

fn filter_batch_for_position(
    batch: &RecordBatch,
    position_idx: usize,
    position_key: u64,
) -> Result<Option<RecordBatch>> {
    let positions = batch
        .column(position_idx)
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| DataFusionError::Execution("position_key must be UInt64".into()))?;
    let mask = BooleanArray::from(
        (0..positions.len())
            .map(|row| !positions.is_null(row) && positions.value(row) == position_key)
            .collect::<Vec<_>>(),
    );
    let filtered = filter_record_batch(batch, &mask)?;
    if filtered.num_rows() == 0 {
        Ok(None)
    } else {
        Ok(Some(filtered))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use datafusion::arrow::array::{
        ListArray, ListBuilder, StringArray, UInt64Array, UInt64Builder,
    };
    use datafusion::arrow::datatypes::{DataType, Field};
    use datafusion::arrow::record_batch::RecordBatch;
    use futures::TryStreamExt;

    #[test]
    fn lance_variation_dataset_path_normalizes_chr_prefix() {
        let root = Path::new("/cache");
        assert_eq!(
            lance_variation_dataset_path(root, "chr1"),
            Path::new("/cache/variation.lance/chr1.lance")
        );
        assert_eq!(
            lance_variation_dataset_path(root, "1"),
            Path::new("/cache/variation.lance/chr1.lance")
        );
    }

    #[test]
    fn lance_position_filter_sorts_and_deduplicates_keys() {
        assert_eq!(
            lance_position_filter(COLD_TIER, &[30, 10, 30, 20]),
            "tier = 1 AND position_key IN (10,20,30)"
        );
    }

    #[tokio::test]
    async fn materializer_writes_warm_then_cold_tiers() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("variation.lance/chr1.lance");
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, true),
        ]));
        let warm = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![10_u64])),
                Arc::new(StringArray::from(vec![Some("A/T")])),
            ],
        )
        .unwrap();
        let cold = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![20_u64])),
                Arc::new(StringArray::from(vec![Some("G/C")])),
            ],
        )
        .unwrap();

        write_merged_lance_variation_dataset(&path, vec![warm], vec![cold], 65_536, 1_024)
            .await
            .unwrap();

        let ds = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let rows = ds
            .scan()
            .project(&["tier", "position_key"])
            .unwrap()
            .try_into_stream()
            .await
            .unwrap()
            .try_collect::<Vec<_>>()
            .await
            .unwrap();

        assert_eq!(rows.iter().map(|batch| batch.num_rows()).sum::<usize>(), 2);
    }

    #[tokio::test]
    async fn materializer_adds_typed_nulls_for_columns_missing_from_one_tier() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("variation.lance/chr1.lance");

        let mut variant_keys = ListBuilder::new(UInt64Builder::new());
        variant_keys.values().append_value(100);
        variant_keys.append(true);
        let variant_keys = Arc::new(variant_keys.finish()) as ArrayRef;

        let warm_schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("variant_keys", variant_keys.data_type().clone(), true),
            Field::new("allele_string", DataType::Utf8, true),
        ]));
        let warm = RecordBatch::try_new(
            warm_schema,
            vec![
                Arc::new(UInt64Array::from(vec![10_u64])),
                variant_keys,
                Arc::new(StringArray::from(vec![Some("A/T")])),
            ],
        )
        .unwrap();

        let cold_schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, true),
        ]));
        let cold = RecordBatch::try_new(
            cold_schema,
            vec![
                Arc::new(UInt64Array::from(vec![20_u64])),
                Arc::new(StringArray::from(vec![Some("G/C")])),
            ],
        )
        .unwrap();

        write_merged_lance_variation_dataset(&path, vec![warm], vec![cold], 65_536, 1_024)
            .await
            .unwrap();

        let ds = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let rows = ds
            .scan()
            .project(&["tier", "position_key", "variant_keys"])
            .unwrap()
            .try_into_stream()
            .await
            .unwrap()
            .try_collect::<Vec<_>>()
            .await
            .unwrap();
        let schema = rows.first().unwrap().schema();
        let batch = concat_batches(&schema, &rows).unwrap();
        let variant_keys = batch
            .column(batch.schema().index_of("variant_keys").unwrap())
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap();

        assert_eq!(batch.num_rows(), 2);
        assert_eq!(variant_keys.null_count(), 1);
    }

    #[tokio::test]
    async fn materializer_creates_scalar_indices() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("variation.lance/chr1.lance");
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, true),
        ]));
        let warm = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![10_u64, 20])),
                Arc::new(StringArray::from(vec![Some("A/T"), Some("G/C")])),
            ],
        )
        .unwrap();
        let cold = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![30_u64, 40])),
                Arc::new(StringArray::from(vec![Some("C/A"), Some("T/G")])),
            ],
        )
        .unwrap();

        write_merged_lance_variation_dataset(&path, vec![warm], vec![cold], 65_536, 1_024)
            .await
            .unwrap();

        let ds = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let indices = ds.load_indices().await.unwrap();

        let index_names = || indices.iter().map(|index| &index.name).collect::<Vec<_>>();
        let position_index = indices
            .iter()
            .find(|index| index.name == LANCE_POSITION_KEY_INDEX_NAME)
            .unwrap_or_else(|| {
                panic!(
                    "expected a Lance scalar index named {LANCE_POSITION_KEY_INDEX_NAME}, got: {:?}",
                    index_names()
                )
            });
        let index_details = position_index
            .index_details
            .as_ref()
            .expect("expected Lance index details metadata");

        assert!(
            index_details.type_url.ends_with("BTreeIndexDetails"),
            "expected position_key_idx to be a BTREE index, got metadata type URL {:?}",
            index_details.type_url
        );

        let tier_index = indices
            .iter()
            .find(|index| index.name == LANCE_TIER_BITMAP_INDEX_NAME)
            .unwrap_or_else(|| {
                panic!(
                    "expected a Lance scalar index named {LANCE_TIER_BITMAP_INDEX_NAME}, got: {:?}",
                    index_names()
                )
            });
        let index_details = tier_index
            .index_details
            .as_ref()
            .expect("expected Lance index details metadata");

        assert!(
            index_details.type_url.ends_with("BitmapIndexDetails"),
            "expected tier_bitmap_idx to be a bitmap index, got metadata type URL {:?}",
            index_details.type_url
        );
    }

    #[test]
    fn lance_scan_profile_line_reports_index_and_pruning_metrics() {
        let stats = LanceScanStats {
            iops: 12,
            requests: 8,
            bytes_read: 2048,
            indices_loaded: 1,
            parts_loaded: 1,
            index_comparisons: 1,
            fragments_scanned: 8,
            ranges_scanned: 8,
            rows_scanned: 3_628_123,
        };

        let line = format_lance_scan_profile_line(
            LanceScanKind::Warm,
            30,
            DEFAULT_LANCE_WARM_SCAN_BATCH_SIZE,
            0,
            3_628_123,
            58,
            Duration::from_millis(388),
            &stats,
        );

        assert_eq!(
            classify_lance_scan_filter("tier = 1 AND position_key IN (10,20,30)"),
            (LanceScanKind::Cold, 3)
        );
        assert!(line.contains("[vep-lance-profile] scan kind=warm"));
        assert!(line.contains("batch_size=16384"));
        assert!(line.contains("index_used=true"));
        assert!(line.contains("indices_loaded=1"));
        assert!(line.contains("fragments_scanned=8"));
        assert!(line.contains("rows_scanned=3628123"));
        assert!(line.contains("scan_s=0.388"));
    }
}
