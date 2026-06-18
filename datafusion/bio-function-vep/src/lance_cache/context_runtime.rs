use std::path::{Path, PathBuf};
use std::time::Instant;

use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use lance::dataset::ProjectionRequest;

use crate::lance_cache::row_index::{
    StringRowIdIndex, U64RowIdIndex, load_sift_key_index_from_lance_btree,
    load_transcript_id_index_from_lance_btree,
};

#[derive(Debug)]
pub struct TranscriptIdLanceLookup {
    dataset: lance::Dataset,
    projection: Vec<String>,
    index: StringRowIdIndex,
    path: PathBuf,
}

impl TranscriptIdLanceLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let profile_enabled = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        let opened_at = profile_enabled.then(Instant::now);
        let dataset = open_lance_dataset(path).await?;
        let index = load_transcript_id_index_from_lance_btree(&dataset).await?;
        let projection = existing_projection_columns(
            dataset_schema_field_names(&dataset),
            projection.iter().map(String::as_str),
        );
        if let Some(started) = opened_at {
            eprintln!(
                "[vep-lance-profile] transcript_lookup_open path={} projection_cols={} row_ids={} unique_values={} seconds={:.3}",
                path.display(),
                projection.len(),
                index.row_ids_len(),
                index.unique_values(),
                started.elapsed().as_secs_f64(),
            );
        }
        Ok(Self {
            dataset,
            projection,
            index,
            path: path.to_path_buf(),
        })
    }

    pub async fn take_transcript_ids(&self, transcript_ids: &[String]) -> Result<RecordBatch> {
        let mut ids = transcript_ids.to_vec();
        ids.sort_unstable();
        ids.dedup();
        let row_ids = self.index.resolve_sorted_values(&ids);
        let profile_enabled = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        let take_started = profile_enabled.then(Instant::now);
        if profile_enabled {
            eprintln!(
                "[vep-lance-profile] transcript_take_start path={} requested_ids={} unique_ids={} row_ids={} projection_cols={}",
                self.path.display(),
                transcript_ids.len(),
                ids.len(),
                row_ids.len(),
                self.projection.len(),
            );
        }
        let projection_request = ProjectionRequest::from_columns(
            self.projection.iter().map(String::as_str),
            self.dataset.schema(),
        );
        let batch = self
            .dataset
            .take_rows(&row_ids, projection_request)
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!("Lance transcript_id take_rows failed: {err}"))
            })?;
        if let Some(started) = take_started {
            eprintln!(
                "[vep-lance-profile] transcript_take path={} row_ids={} batch_rows={} seconds={:.3}",
                self.path.display(),
                row_ids.len(),
                batch.num_rows(),
                started.elapsed().as_secs_f64(),
            );
        }
        Ok(batch)
    }

    pub fn row_ids_len(&self) -> usize {
        self.index.row_ids_len()
    }

    pub fn unique_values(&self) -> usize {
        self.index.unique_values()
    }
}

/// Position-sliced SIFT lookup: loads the `key` BTree index into an in-memory
/// `u64 key -> row_id` map at open, then resolves queried keys to a batched
/// `take_rows([key, sift, poly])`. Mirrors [`TranscriptIdLanceLookup`] but for
/// the bit-packed `(transcript_uid, position)` key.
#[derive(Debug)]
pub struct KeyU64LanceLookup {
    dataset: lance::Dataset,
    projection: Vec<String>,
    index: U64RowIdIndex,
    path: PathBuf,
}

impl KeyU64LanceLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let profile_enabled = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        let opened_at = profile_enabled.then(Instant::now);
        let dataset = open_lance_dataset(path).await?;
        let index = load_sift_key_index_from_lance_btree(&dataset).await?;
        let projection = existing_projection_columns(
            dataset_schema_field_names(&dataset),
            projection.iter().map(String::as_str),
        );
        if let Some(started) = opened_at {
            eprintln!(
                "[vep-lance-profile] sift_key_lookup_open path={} projection_cols={} keys={} seconds={:.3}",
                path.display(),
                projection.len(),
                index.len(),
                started.elapsed().as_secs_f64(),
            );
        }
        Ok(Self {
            dataset,
            projection,
            index,
            path: path.to_path_buf(),
        })
    }

    /// Resolve `keys` to row ids (skipping keys absent from the index) and take
    /// the projected columns. Returns the fetched batch plus the keys that were
    /// present, aligned to the batch's rows in row-id order.
    pub async fn take_keys(&self, keys: &[u64]) -> Result<(RecordBatch, Vec<u64>)> {
        let mut sorted = keys.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        let (row_ids, present) = self.index.resolve(&sorted);
        let profile_enabled = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        let take_started = profile_enabled.then(Instant::now);
        let projection_request = ProjectionRequest::from_columns(
            self.projection.iter().map(String::as_str),
            self.dataset.schema(),
        );
        let batch = self
            .dataset
            .take_rows(&row_ids, projection_request)
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!("Lance sift key take_rows failed: {err}"))
            })?;
        if let Some(started) = take_started {
            eprintln!(
                "[vep-lance-profile] sift_key_take path={} requested_keys={} unique_keys={} row_ids={} batch_rows={} seconds={:.3}",
                self.path.display(),
                keys.len(),
                sorted.len(),
                row_ids.len(),
                batch.num_rows(),
                started.elapsed().as_secs_f64(),
            );
        }
        Ok((batch, present))
    }

    pub fn keys_len(&self) -> usize {
        self.index.len()
    }
}

pub async fn scan_projected_existing_columns(
    path: &Path,
    requested_columns: &[&str],
) -> Result<Vec<RecordBatch>> {
    let dataset = open_lance_dataset(path).await?;
    let projection = existing_projection_columns(
        dataset_schema_field_names(&dataset),
        requested_columns.iter().copied(),
    );
    let projection_refs = projection.iter().map(String::as_str).collect::<Vec<_>>();

    let mut scanner = dataset.scan();
    scanner.project(&projection_refs).map_err(|err| {
        DataFusionError::Execution(format!(
            "invalid Lance projection for '{}': {err}",
            path.display()
        ))
    })?;
    let mut stream = scanner.try_into_stream().await.map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to scan Lance context dataset '{}': {err}",
            path.display()
        ))
    })?;

    let mut batches = Vec::new();
    while let Some(batch) = stream.try_next().await.map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to read Lance context dataset '{}': {err}",
            path.display()
        ))
    })? {
        batches.push(batch);
    }
    Ok(batches)
}

async fn open_lance_dataset(path: &Path) -> Result<lance::Dataset> {
    lance::Dataset::open(path.to_string_lossy().as_ref())
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to open Lance context dataset '{}': {err}",
                path.display()
            ))
        })
}

fn dataset_schema_field_names(dataset: &lance::Dataset) -> Vec<String> {
    let schema: datafusion::arrow::datatypes::Schema = dataset.schema().into();
    schema
        .fields()
        .iter()
        .map(|field| field.name().clone())
        .collect()
}

fn existing_projection_columns<'a, I>(field_names: Vec<String>, requested: I) -> Vec<String>
where
    I: IntoIterator<Item = &'a str>,
{
    let available = field_names
        .into_iter()
        .collect::<std::collections::HashSet<_>>();
    let mut projection = Vec::new();
    for requested in requested {
        let name = normalize_projection_name(requested);
        if available.contains(name) && !projection.iter().any(|column| column == name) {
            projection.push(name.to_string());
        }
    }
    projection
}

fn normalize_projection_name(name: &str) -> &str {
    name.strip_prefix('"')
        .and_then(|rest| rest.strip_suffix('"'))
        .unwrap_or(name)
}

#[cfg(test)]
mod tests {
    use datafusion::arrow::array::{
        BinaryArray, Int32Array, StringArray, UInt32Array, UInt64Array,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[tokio::test]
    async fn sift_key_lookup_resolves_keys_via_take_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("translation_sift.lance");
        let schema = Arc::new(Schema::new(vec![
            Field::new("key", DataType::UInt64, false),
            Field::new("sift", DataType::Binary, true),
            Field::new("poly", DataType::Binary, true),
        ]));
        // Three position rows for two transcript uids; keys must be ascending.
        let k0 = (1u64 << 32) | 10;
        let k1 = (1u64 << 32) | 11;
        let k2 = (2u64 << 32) | 5;
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![k0, k1, k2])),
                Arc::new(BinaryArray::from_opt_vec(vec![
                    Some(b"a".as_ref()),
                    Some(b"b".as_ref()),
                    Some(b"c".as_ref()),
                ])),
                Arc::new(BinaryArray::from_opt_vec(vec![
                    Some(b"x".as_ref()),
                    None,
                    Some(b"z".as_ref()),
                ])),
            ],
        )
        .unwrap();
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::SiftKey,
        )
        .await
        .unwrap();

        let lookup =
            super::KeyU64LanceLookup::open(&path, vec!["key".into(), "sift".into(), "poly".into()])
                .await
                .unwrap();
        assert_eq!(lookup.keys_len(), 3);

        // Query a present key plus a missing one; only the present key resolves.
        let (batch, present) = lookup.take_keys(&[k1, (9u64 << 32) | 1]).await.unwrap();
        assert_eq!(present, vec![k1]);
        assert_eq!(batch.num_rows(), 1);
        let keys = batch
            .column_by_name("key")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();
        assert_eq!(keys.value(0), k1);
    }

    #[tokio::test]
    async fn transcript_id_payload_lookup_uses_take_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("translation_core.lance");
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int32, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST0001", "ENST0002", "ENST0003"])),
                Arc::new(Int32Array::from(vec![100, 200, 300])),
            ],
        )
        .unwrap();
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::TranscriptId,
        )
        .await
        .unwrap();

        let lookup = super::TranscriptIdLanceLookup::open(
            &path,
            vec!["transcript_id".into(), "cds_len".into()],
        )
        .await
        .unwrap();
        let batch = lookup
            .take_transcript_ids(&["ENST0002".to_string()])
            .await
            .unwrap();

        assert_eq!(batch.num_rows(), 1);
        assert_eq!(batch.schema().field(0).name(), "transcript_id");
    }

    #[tokio::test]
    async fn projected_scan_reads_existing_lance_columns_without_sql() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("transcript.lance");
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST0001"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(UInt32Array::from(vec![10])),
            ],
        )
        .unwrap();
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();

        let batches = super::scan_projected_existing_columns(
            &path,
            &["transcript_id", "missing_column", "\"start\""],
        )
        .await
        .unwrap();

        assert_eq!(batches.iter().map(RecordBatch::num_rows).sum::<usize>(), 1);
        assert_eq!(
            batches[0]
                .schema()
                .fields()
                .iter()
                .map(|field| field.name().as_str())
                .collect::<Vec<_>>(),
            ["transcript_id", "start"]
        );
    }
}
