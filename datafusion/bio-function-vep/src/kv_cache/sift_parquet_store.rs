//! Parquet-backed SIFT/PolyPhen prediction store.
//!
//! This store is intended for the compact lookup layout:
//! `transcript_id: Utf8`, `predictions: Binary`, where `predictions` uses the
//! same byte representation as the fjall SIFT keyspace.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, BinaryArray, LargeBinaryArray, LargeStringArray, RecordBatch, StringArray,
    StringViewArray,
};
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReaderBuilder, RowSelection,
};

use super::sift_store::{SiftPredictionStore, deserialize_predictions};
use crate::transcript_consequence::CachedPredictions;

const DEFAULT_BATCH_SIZE: usize = 8192;

#[derive(Clone)]
pub struct SiftParquetStore {
    inner: Arc<SiftParquetStoreInner>,
}

struct SiftParquetStoreInner {
    files: Vec<SiftParquetFile>,
    index: HashMap<String, Vec<RowRef>>,
}

struct SiftParquetFile {
    path: PathBuf,
    metadata: ArrowReaderMetadata,
    id_projection: ProjectionMask,
    predictions_projection: ProjectionMask,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct RowGroupRef {
    file_id: usize,
    row_group_id: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct RowRef {
    file_id: usize,
    row_group_id: usize,
    row_index: usize,
}

impl SiftParquetStore {
    pub fn open_dir(path: impl AsRef<Path>) -> Result<Option<Self>> {
        let path = path.as_ref();
        if !path.exists() {
            return Ok(None);
        }

        let mut paths = if path.is_file() {
            vec![path.to_path_buf()]
        } else {
            let mut paths = std::fs::read_dir(path)
                .map_err(|error| {
                    DataFusionError::Execution(format!(
                        "failed to list sift parquet directory '{}': {error}",
                        path.display()
                    ))
                })?
                .filter_map(|entry| entry.ok().map(|entry| entry.path()))
                .filter(|path| {
                    path.extension()
                        .and_then(|ext| ext.to_str())
                        .is_some_and(|ext| ext.eq_ignore_ascii_case("parquet"))
                })
                .collect::<Vec<_>>();
            paths.sort();
            paths
        };
        paths.sort();

        if paths.is_empty() {
            return Ok(None);
        }

        let mut files = Vec::with_capacity(paths.len());
        let mut index: HashMap<String, Vec<RowRef>> = HashMap::new();
        for path in paths {
            let file_id = files.len();
            let file = Self::open_file(file_id, path)?;
            Self::index_file(file_id, &file, &mut index)?;
            files.push(file);
        }

        Ok(Some(Self {
            inner: Arc::new(SiftParquetStoreInner { files, index }),
        }))
    }

    fn open_file(_file_id: usize, path: PathBuf) -> Result<SiftParquetFile> {
        let file = File::open(&path).map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to open sift lookup parquet '{}': {error}",
                path.display()
            ))
        })?;
        let metadata =
            ArrowReaderMetadata::load(&file, ArrowReaderOptions::default()).map_err(|error| {
                DataFusionError::Execution(format!(
                    "failed to load sift lookup parquet metadata '{}': {error}",
                    path.display()
                ))
            })?;
        let id_projection = required_projection(&metadata, &["transcript_id"])?;
        let predictions_projection = required_projection(&metadata, &["predictions"])?;
        Ok(SiftParquetFile {
            path,
            metadata,
            id_projection,
            predictions_projection,
        })
    }

    fn index_file(
        file_id: usize,
        file: &SiftParquetFile,
        index: &mut HashMap<String, Vec<RowRef>>,
    ) -> Result<()> {
        for row_group_id in 0..file.metadata.metadata().num_row_groups() {
            let mut row_base = 0usize;
            for batch in read_row_group(file, row_group_id, file.id_projection.clone(), None)? {
                let tx_idx = batch.schema().index_of("transcript_id")?;
                let tx_col = batch.column(tx_idx).as_ref();
                for row in 0..batch.num_rows() {
                    let Some(transcript_id) = string_value(tx_col, row)? else {
                        continue;
                    };
                    let row_ref = RowRef {
                        file_id,
                        row_group_id,
                        row_index: row_base + row,
                    };
                    let entry = index.entry(transcript_id.to_string()).or_default();
                    if !entry.contains(&row_ref) {
                        entry.push(row_ref);
                    }
                }
                row_base += batch.num_rows();
            }
        }
        Ok(())
    }
}

impl SiftPredictionStore for SiftParquetStore {
    fn get_many(&self, transcript_ids: &[String]) -> Result<HashMap<String, CachedPredictions>> {
        let mut grouped: HashMap<RowGroupRef, Vec<(usize, String)>> = HashMap::new();
        let mut seen = HashSet::with_capacity(transcript_ids.len());
        for transcript_id in transcript_ids {
            if !seen.insert(transcript_id.as_str()) {
                continue;
            }
            let Some(row_refs) = self.inner.index.get(transcript_id.as_str()) else {
                continue;
            };
            for row_ref in row_refs {
                grouped
                    .entry(RowGroupRef {
                        file_id: row_ref.file_id,
                        row_group_id: row_ref.row_group_id,
                    })
                    .or_default()
                    .push((row_ref.row_index, transcript_id.clone()));
            }
        }

        let mut row_groups = grouped.into_iter().collect::<Vec<_>>();
        row_groups.sort_by_key(|(row_group, _)| (row_group.file_id, row_group.row_group_id));

        let mut out = HashMap::with_capacity(transcript_ids.len());
        for (row_group, mut selected_rows) in row_groups {
            let file = self.inner.files.get(row_group.file_id).ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "sift lookup parquet index referenced missing file {}",
                    row_group.file_id
                ))
            })?;
            selected_rows.sort_by_key(|(row_index, _)| *row_index);
            selected_rows.dedup_by_key(|(row_index, _)| *row_index);
            let row_group_rows = file
                .metadata
                .metadata()
                .row_group(row_group.row_group_id)
                .num_rows() as usize;
            let ranges = consecutive_ranges(&selected_rows);
            let row_selection =
                if ranges.len() == 1 && ranges[0].start == 0 && ranges[0].end == row_group_rows {
                    None
                } else {
                    Some(RowSelection::from_consecutive_ranges(
                        ranges.into_iter(),
                        row_group_rows,
                    ))
                };

            let mut selected_cursor = 0usize;
            for batch in read_row_group(
                file,
                row_group.row_group_id,
                file.predictions_projection.clone(),
                row_selection,
            )? {
                append_selected_predictions(
                    &batch,
                    &selected_rows,
                    &mut selected_cursor,
                    &mut out,
                )?;
                if selected_cursor >= selected_rows.len() {
                    break;
                }
            }
        }
        Ok(out)
    }
}

fn append_selected_predictions(
    batch: &RecordBatch,
    selected_rows: &[(usize, String)],
    selected_cursor: &mut usize,
    out: &mut HashMap<String, CachedPredictions>,
) -> Result<()> {
    let schema = batch.schema();
    let predictions_idx = schema.index_of("predictions")?;
    let predictions_col = batch.column(predictions_idx).as_ref();

    for row in 0..batch.num_rows() {
        let Some((_, transcript_id)) = selected_rows.get(*selected_cursor) else {
            break;
        };
        *selected_cursor += 1;
        if out.contains_key(transcript_id) {
            continue;
        }
        let Some(bytes) = binary_value(predictions_col, row)? else {
            continue;
        };
        out.insert(transcript_id.clone(), deserialize_predictions(bytes)?);
    }
    Ok(())
}

fn consecutive_ranges(selected_rows: &[(usize, String)]) -> Vec<Range<usize>> {
    let Some((first, _)) = selected_rows.first() else {
        return Vec::new();
    };
    let mut ranges = Vec::new();
    let mut start = *first;
    let mut end = start + 1;
    for (row_index, _) in selected_rows.iter().skip(1) {
        if *row_index == end {
            end += 1;
        } else {
            ranges.push(start..end);
            start = *row_index;
            end = start + 1;
        }
    }
    ranges.push(start..end);
    ranges
}

fn read_row_group(
    file: &SiftParquetFile,
    row_group_id: usize,
    projection: ProjectionMask,
    row_selection: Option<RowSelection>,
) -> Result<Vec<RecordBatch>> {
    let parquet_file = File::open(&file.path).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to open sift lookup parquet '{}': {error}",
            file.path.display()
        ))
    })?;
    let mut builder =
        ParquetRecordBatchReaderBuilder::new_with_metadata(parquet_file, file.metadata.clone())
            .with_projection(projection)
            .with_row_groups(vec![row_group_id])
            .with_batch_size(DEFAULT_BATCH_SIZE);
    if let Some(row_selection) = row_selection {
        builder = builder.with_row_selection(row_selection);
    }
    let reader = builder.build().map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to build sift lookup parquet reader '{}': {error}",
            file.path.display()
        ))
    })?;

    reader
        .map(|batch| {
            batch.map_err(|error| {
                DataFusionError::Execution(format!(
                    "failed to read sift lookup parquet row group {row_group_id} from '{}': {error}",
                    file.path.display()
                ))
            })
        })
        .collect()
}

fn required_projection(metadata: &ArrowReaderMetadata, names: &[&str]) -> Result<ProjectionMask> {
    let fields = metadata.schema().fields();
    let mut roots = Vec::with_capacity(names.len());
    for name in names {
        let Some(idx) = fields.iter().position(|field| field.name() == *name) else {
            return Err(DataFusionError::Execution(format!(
                "sift lookup parquet missing required column '{name}'"
            )));
        };
        roots.push(idx);
    }
    Ok(ProjectionMask::roots(metadata.parquet_schema(), roots))
}

fn string_value(array: &dyn Array, row: usize) -> Result<Option<&str>> {
    if row >= array.len() || array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<StringArray>() {
        Ok(Some(array.value(row)))
    } else if let Some(array) = array.as_any().downcast_ref::<StringViewArray>() {
        Ok(Some(array.value(row)))
    } else if let Some(array) = array.as_any().downcast_ref::<LargeStringArray>() {
        Ok(Some(array.value(row)))
    } else {
        Err(DataFusionError::Execution(format!(
            "sift lookup parquet transcript_id expected string array, got {:?}",
            array.data_type()
        )))
    }
}

fn binary_value(array: &dyn Array, row: usize) -> Result<Option<&[u8]>> {
    if row >= array.len() || array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<BinaryArray>() {
        Ok(Some(array.value(row)))
    } else if let Some(array) = array.as_any().downcast_ref::<LargeBinaryArray>() {
        Ok(Some(array.value(row)))
    } else {
        Err(DataFusionError::Execution(format!(
            "sift lookup parquet predictions expected binary array, got {:?}",
            array.data_type()
        )))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use datafusion::arrow::array::{ArrayRef, BinaryArray, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use datafusion::common::Result;
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;

    use super::*;
    use crate::kv_cache::sift_store::{SiftPredictionStore, serialize_predictions};
    use crate::transcript_consequence::{CachedPredictions, CompactPrediction};

    fn make_predictions(position: i32, score: f32) -> CachedPredictions {
        CachedPredictions {
            sift: vec![CompactPrediction {
                position,
                amino_acid: 0,
                prediction: 1,
                score,
            }],
            polyphen: vec![CompactPrediction {
                position: position + 1,
                amino_acid: 1,
                prediction: 2,
                score: score + 0.1,
            }],
        }
    }

    fn write_lookup_parquet(path: &std::path::Path) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("predictions", DataType::Binary, false),
        ]));
        let tx1 = make_predictions(10, 0.2);
        let tx2 = make_predictions(20, 0.3);
        let predictions = vec![serialize_predictions(&tx1), serialize_predictions(&tx2)];
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST0001", "ENST0002"])) as ArrayRef,
                Arc::new(BinaryArray::from_iter_values(predictions)) as ArrayRef,
            ],
        )
        .unwrap();
        let file = std::fs::File::create(path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(1)
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props)).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
    }

    #[test]
    fn parquet_sift_prediction_store_get_many_returns_found_and_skips_missing() -> Result<()> {
        let dir = tempfile::tempdir().unwrap();
        write_lookup_parquet(&dir.path().join("chr4.parquet"));

        let store = SiftParquetStore::open_dir(dir.path())?.expect("parquet store");
        let ids = vec![
            "ENST0002".to_string(),
            "missing".to_string(),
            "ENST0001".to_string(),
            "ENST0002".to_string(),
        ];

        let found = SiftPredictionStore::get_many(&store, &ids)?;

        assert_eq!(found.len(), 2);
        assert_eq!(found["ENST0001"].sift[0].position, 10);
        assert_eq!(found["ENST0001"].polyphen[0].position, 11);
        assert_eq!(found["ENST0002"].sift[0].position, 20);
        Ok(())
    }
}
