use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{ArrayRef, UInt8Array};
use datafusion::arrow::datatypes::Schema;
use datafusion::arrow::record_batch::{RecordBatch, RecordBatchIterator};
use datafusion::common::{DataFusionError, Result};
use lance::dataset::{WriteMode, WriteParams};
use lance_file::version::LanceFileVersion;

pub const LANCE_VARIATION_DIR: &str = "variation.lance";
pub const WARM_TIER: u8 = 0;
pub const COLD_TIER: u8 = 1;

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
    fields.push(datafusion::arrow::datatypes::Field::new(
        "tier",
        datafusion::arrow::datatypes::DataType::UInt8,
        false,
    ));
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
