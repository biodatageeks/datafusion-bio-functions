use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, ArrayRef, BooleanArray, UInt8Array, UInt64Array};
use datafusion::arrow::compute::{concat_batches, filter_record_batch};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::arrow::record_batch::{RecordBatch, RecordBatchIterator};
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use lance::dataset::{WriteMode, WriteParams};
use lance_file::version::LanceFileVersion;

use crate::warm_cache::chrom_cache::WarmProbe;
use crate::warm_cache::chunk::{WarmChunkContext, WarmChunkProbe};
use crate::warm_cache::cold_parquet::ColdProbeResult;

pub const LANCE_VARIATION_DIR: &str = "variation.lance";
pub const WARM_TIER: u8 = 0;
pub const COLD_TIER: u8 = 1;
pub const DEFAULT_LANCE_COLD_LOOKUP_BATCH_SIZE: usize = 2_000;

const LANCE_REQUIRED_RUNTIME_COLUMNS: &[&str] = &[
    "position_key",
    "variant_keys",
    "allele_string",
    "end",
    "failed",
];

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
        let mut scanner = self.dataset.scan();
        scanner
            .filter(filter)
            .map_err(|err| DataFusionError::Execution(format!("invalid Lance filter: {err}")))?;
        scanner.project(&projection).map_err(|err| {
            DataFusionError::Execution(format!("invalid Lance projection: {err}"))
        })?;
        scanner.batch_size(self.batch_size);

        let mut stream = scanner.try_into_stream().await.map_err(|err| {
            DataFusionError::Execution(format!("failed to scan Lance variation dataset: {err}"))
        })?;
        let mut batches = Vec::new();
        while let Some(batch) = stream.try_next().await.map_err(|err| {
            DataFusionError::Execution(format!("failed to read Lance variation batch: {err}"))
        })? {
            batches.push(batch);
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

    let batches = warm_batches
        .into_iter()
        .map(|batch| append_tier_column(batch, WARM_TIER))
        .chain(
            cold_batches
                .into_iter()
                .map(|batch| append_tier_column(batch, COLD_TIER)),
        )
        .collect::<Result<Vec<_>>>()?;

    let schema = batches
        .first()
        .ok_or_else(|| {
            DataFusionError::Execution("cannot write empty Lance variation dataset".into())
        })?
        .schema();
    let reader = RecordBatchIterator::new(batches.into_iter().map(Ok), schema);

    let params = WriteParams {
        mode: WriteMode::Overwrite,
        max_rows_per_file: cold_fragment_rows.max(1),
        max_rows_per_group: warm_fragment_rows.max(1),
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

    use datafusion::arrow::array::{StringArray, UInt64Array};
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
}
