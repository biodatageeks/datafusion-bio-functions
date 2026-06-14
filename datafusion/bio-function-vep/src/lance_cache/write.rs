use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::arrow::record_batch::{RecordBatch, RecordBatchIterator};
use datafusion::common::{DataFusionError, Result};
use lance::dataset::{WriteMode, WriteParams};
use lance::index::DatasetIndexExt;
use lance_file::version::LanceFileVersion;
use lance_index::{
    IndexType,
    scalar::{BuiltinIndexType, ScalarIndexParams},
};

use crate::lance_cache::schema::with_lance_field_metadata;

pub const START_BTREE_INDEX_NAME: &str = "start_btree_idx";
pub const TRANSCRIPT_ID_BTREE_INDEX_NAME: &str = "transcript_id_btree_idx";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LanceIndexKind {
    Start,
    TranscriptId,
}

pub async fn write_record_batches_to_lance(
    path: &Path,
    batches: Vec<RecordBatch>,
    index_kind: LanceIndexKind,
) -> Result<()> {
    write_record_batches_to_lance_with_mode(path, batches, WriteMode::Overwrite).await?;
    create_required_index(path, index_kind).await
}

pub async fn write_record_batches_to_lance_with_mode(
    path: &Path,
    batches: Vec<RecordBatch>,
    mode: WriteMode,
) -> Result<()> {
    write_record_batches_to_lance_with_mode_and_version(path, batches, mode, LanceFileVersion::V2_1)
        .await
}

pub async fn write_record_batches_to_lance_with_mode_and_version(
    path: &Path,
    batches: Vec<RecordBatch>,
    mode: WriteMode,
    data_storage_version: LanceFileVersion,
) -> Result<()> {
    if batches.is_empty() {
        return Err(DataFusionError::Execution(format!(
            "cannot write empty Lance dataset at {}",
            path.display()
        )));
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|err| {
            DataFusionError::Execution(format!("failed to create {}: {err}", parent.display()))
        })?;
    }

    let schema = Arc::new(with_lance_field_metadata(batches[0].schema().as_ref()));
    let aligned = batches
        .into_iter()
        .map(|batch| align_batch_to_schema(batch, Arc::clone(&schema)))
        .collect::<Result<Vec<_>>>()?;
    let reader = RecordBatchIterator::new(aligned.into_iter().map(Ok), Arc::clone(&schema));
    let params = WriteParams {
        mode,
        data_storage_version: Some(data_storage_version),
        ..Default::default()
    };

    lance::Dataset::write(reader, path.to_string_lossy().as_ref(), Some(params))
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to write Lance dataset '{}': {err}",
                path.display()
            ))
        })?;
    Ok(())
}

fn align_batch_to_schema(batch: RecordBatch, schema: Arc<Schema>) -> Result<RecordBatch> {
    let columns = schema
        .fields()
        .iter()
        .map(|field| {
            let index = batch.schema().index_of(field.name()).map_err(|err| {
                DataFusionError::Execution(format!(
                    "batch missing Lance output field '{}': {err}",
                    field.name()
                ))
            })?;
            Ok(batch.column(index).clone())
        })
        .collect::<Result<Vec<_>>>()?;
    RecordBatch::try_new(schema, columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to align batch to Lance schema: {err}"))
    })
}

pub async fn create_required_index(path: &Path, index_kind: LanceIndexKind) -> Result<()> {
    let (column, name) = match index_kind {
        LanceIndexKind::Start => ("start", START_BTREE_INDEX_NAME),
        LanceIndexKind::TranscriptId => ("transcript_id", TRANSCRIPT_ID_BTREE_INDEX_NAME),
    };
    let mut dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to open Lance dataset for indexing '{}': {err}",
                path.display()
            ))
        })?;

    dataset
        .create_index(
            &[column],
            IndexType::BTree,
            Some(name.to_string()),
            &ScalarIndexParams::for_builtin(BuiltinIndexType::BTree),
            true,
        )
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to create Lance BTree index '{name}' on column '{column}' for '{}': {err}",
                path.display()
            ))
        })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field};
    use std::sync::Arc;

    #[tokio::test]
    async fn writes_lance_dataset_and_start_btree() {
        let dir = tempfile::tempdir().unwrap();
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["chr1", "chr1"])),
                Arc::new(UInt32Array::from(vec![10, 20])),
                Arc::new(UInt32Array::from(vec![10, 20])),
                Arc::new(StringArray::from(vec!["A/T", "G/C"])),
                Arc::new(Int8Array::from(vec![Some(0), Some(0)])),
                Arc::new(StringArray::from(vec!["rs1", "rs2"])),
            ],
        )
        .unwrap();
        let dataset_path = dir.path().join("variation/chr1.lance");

        write_record_batches_to_lance(&dataset_path, vec![batch], LanceIndexKind::Start)
            .await
            .unwrap();

        let dataset = lance::Dataset::open(dataset_path.to_string_lossy().as_ref())
            .await
            .unwrap();
        let indexes = dataset.load_indices().await.unwrap();
        let index = indexes
            .iter()
            .find(|idx| idx.name == START_BTREE_INDEX_NAME)
            .unwrap_or_else(|| panic!("missing {START_BTREE_INDEX_NAME}: {indexes:?}"));
        let details = index
            .index_details
            .as_ref()
            .expect("expected Lance index details metadata");
        assert!(
            details.type_url.ends_with("BTreeIndexDetails"),
            "expected BTree index details, got {:?}",
            details.type_url
        );
    }

    #[tokio::test]
    async fn appends_lance_chunks_before_creating_start_btree() {
        let dir = tempfile::tempdir().unwrap();
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ]));
        let first = RecordBatch::try_new(
            Arc::clone(&schema),
            vec![
                Arc::new(StringArray::from(vec!["chr1"])),
                Arc::new(UInt32Array::from(vec![10])),
                Arc::new(UInt32Array::from(vec![10])),
                Arc::new(StringArray::from(vec!["A/T"])),
                Arc::new(Int8Array::from(vec![Some(0)])),
                Arc::new(StringArray::from(vec!["rs1"])),
            ],
        )
        .unwrap();
        let second = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["chr1"])),
                Arc::new(UInt32Array::from(vec![20])),
                Arc::new(UInt32Array::from(vec![20])),
                Arc::new(StringArray::from(vec!["G/C"])),
                Arc::new(Int8Array::from(vec![Some(0)])),
                Arc::new(StringArray::from(vec!["rs2"])),
            ],
        )
        .unwrap();
        let dataset_path = dir.path().join("variation/chr1.lance");

        write_record_batches_to_lance_with_mode(&dataset_path, vec![first], WriteMode::Overwrite)
            .await
            .unwrap();
        write_record_batches_to_lance_with_mode(&dataset_path, vec![second], WriteMode::Append)
            .await
            .unwrap();
        create_required_index(&dataset_path, LanceIndexKind::Start)
            .await
            .unwrap();

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
                .any(|idx| idx.name == START_BTREE_INDEX_NAME)
        );
    }

    #[test]
    fn aligns_batches_to_metadata_enriched_schema() {
        let source_schema = Arc::new(Schema::new(vec![
            Field::new("start", DataType::UInt32, false),
            Field::new("chrom", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            source_schema,
            vec![
                Arc::new(UInt32Array::from(vec![10, 20])),
                Arc::new(StringArray::from(vec!["chr1", "chr1"])),
            ],
        )
        .unwrap();
        let target = Arc::new(with_lance_field_metadata(&Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
        ])));

        let aligned = align_batch_to_schema(batch, target).unwrap();

        assert_eq!(aligned.schema().field(0).name(), "chrom");
        assert_eq!(aligned.schema().field(1).name(), "start");
        assert_eq!(aligned.column(0).len(), 2);
    }
}
