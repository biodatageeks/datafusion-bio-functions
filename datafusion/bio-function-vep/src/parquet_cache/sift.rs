//! Parquet point-lookup for the `translation_sift` cache entity.
//!
//! Position-sliced SIFT/PolyPhen blobs keyed by a `u64` `(transcript_uid<<32) |
//! position`. Mirrors the Lance `KeyU64LanceLookup` contract — `open` + a
//! `take_keys(&[u64]) -> (RecordBatch, present)` — but resolves keys through the
//! footer [`PageDir`] (u64 key column) + [`CoalescingAsyncReader`] instead of an
//! in-memory BTree. The shard is written no-dictionary, page-indexed, and
//! physically sorted by `key` (a single monotone run — SIFT keys are unique), so
//! the reader is the same three-phase loop as the variation lookup: resolve
//! candidate pages → exact key offsets → projected payload take.
//!
//! The returned batch (`key: UInt64`, `sift: Binary`, `poly: Binary`) is the same
//! shape the Lance path yields, so `position_predictions_from_batch` decodes it
//! unchanged.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, RecordBatch, UInt64Array};
use datafusion::arrow::compute::concat_batches;
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use parquet::arrow::ArrowWriter;
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};
use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;

use crate::parquet_cache::page_dir::{
    CoalescingAsyncReader, IoCounters, PageDir, selection_from_offsets, selection_from_ranges,
};
use crate::parquet_cache::write::point_lookup_writer_properties;

/// Coalescing-reader window (bytes) — matches the variation lookup.
const COALESCE_GAP_BYTES: u64 = 64 * 1024;

/// Write already-sorted (`key` ascending) translation_sift batches to a single
/// no-dictionary, page-indexed Parquet shard. Rows MUST be globally sorted by
/// `key` so the footer `PageDir` sees one monotone run.
pub fn write_sift_parquet(path: &Path, batches: &[RecordBatch]) -> Result<usize> {
    if batches.is_empty() {
        return Ok(0);
    }
    let schema = batches[0].schema();
    let props = point_lookup_writer_properties(&schema, &["key"]);
    let file = std::fs::File::create(path).map_err(|e| {
        DataFusionError::Execution(format!("create parquet '{}': {e}", path.display()))
    })?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| DataFusionError::Execution(format!("open ArrowWriter: {e}")))?;
    let mut rows = 0usize;
    for batch in batches {
        writer
            .write(batch)
            .map_err(|e| DataFusionError::Execution(format!("write batch: {e}")))?;
        rows += batch.num_rows();
    }
    writer
        .close()
        .map_err(|e| DataFusionError::Execution(format!("close ArrowWriter: {e}")))?;
    Ok(rows)
}

/// Per-chromosome Parquet translation_sift point-lookup.
pub struct SinglePathParquetSiftLookup {
    path: PathBuf,
    meta: ArrowReaderMetadata,
    page_dir: PageDir,
    key_leaf: usize,
    /// Physical top-level columns to read for the payload take (always incl. `key`).
    projection: Vec<String>,
}

impl SinglePathParquetSiftLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let file = std::fs::File::open(path).map_err(|e| {
            DataFusionError::Execution(format!("open parquet sift '{}': {e}", path.display()))
        })?;
        // `PageIndexPolicy` is private in parquet 58, so the deprecated setter is
        // the only usable way to force the page index on.
        #[allow(deprecated)]
        let meta =
            ArrowReaderMetadata::load(&file, ArrowReaderOptions::new().with_page_index(true))
                .map_err(|e| DataFusionError::Execution(format!("load parquet metadata: {e}")))?;
        let key_leaf = meta
            .parquet_schema()
            .columns()
            .iter()
            .position(|c| c.name() == "key")
            .ok_or_else(|| {
                DataFusionError::Execution("sift parquet has no 'key' column".to_string())
            })?;
        let page_dir = PageDir::build(&meta, key_leaf)?;
        let mut projection = projection;
        if !projection.iter().any(|c| c == "key") {
            projection.push("key".to_string());
        }
        Ok(Self {
            path: path.to_path_buf(),
            meta,
            page_dir,
            key_leaf,
            projection,
        })
    }

    /// Resolve `keys` to their rows. Returns the payload batch (`key`/`sift`/`poly`
    /// for the matched keys, in `key` order) and the distinct keys found — the
    /// same `(RecordBatch, present)` contract as the Lance `take_keys`.
    pub async fn take_keys(&self, keys: &[u64]) -> Result<(RecordBatch, Vec<u64>)> {
        let mut sorted: Vec<u64> = keys.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        let probe_set: HashSet<u64> = sorted.iter().copied().collect();
        let counters = IoCounters::new();

        let ranges = self.page_dir.resolve_ranges(&sorted);
        let offsets = self.exact_offsets(&ranges, &probe_set, &counters).await?;
        let batch = self.take_payload(&offsets, &counters).await?;
        let present = distinct_present(&batch, &probe_set)?;
        Ok((batch, present))
    }

    async fn exact_offsets(
        &self,
        ranges: &[(u64, u64)],
        probe_set: &HashSet<u64>,
        counters: &IoCounters,
    ) -> Result<Vec<u64>> {
        if ranges.is_empty() {
            return Ok(Vec::new());
        }
        let pq_schema = self.meta.parquet_schema();
        let key_mask = ProjectionMask::leaves(pq_schema, [self.key_leaf]);
        let reader = CoalescingAsyncReader::new(
            self.open_async().await?,
            counters.clone(),
            COALESCE_GAP_BYTES,
        );
        let mut stream =
            ParquetRecordBatchStreamBuilder::new_with_metadata(reader, self.meta.clone())
                .with_projection(key_mask)
                .with_row_selection(selection_from_ranges(ranges))
                .with_batch_size(8192)
                .build()
                .map_err(|e| DataFusionError::Execution(format!("build key stream: {e}")))?;
        let mut offsets = Vec::new();
        let mut cursor = ranges.iter().flat_map(|&(a, b)| a..b);
        while let Some(b) = stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read key batch: {e}")))?
        {
            let ka = b
                .column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| DataFusionError::Execution("key column must be UInt64".into()))?;
            for &v in ka.values() {
                let off = cursor.next().ok_or_else(|| {
                    DataFusionError::Execution("row range cursor underflow".into())
                })?;
                if probe_set.contains(&v) {
                    offsets.push(off);
                }
            }
        }
        Ok(offsets)
    }

    async fn take_payload(&self, offsets: &[u64], counters: &IoCounters) -> Result<RecordBatch> {
        let pq_schema = self.meta.parquet_schema();
        let arrow_schema = self.meta.schema();
        let root_indices: Vec<usize> = self
            .projection
            .iter()
            .filter_map(|n| arrow_schema.index_of(n).ok())
            .collect();
        let mask = ProjectionMask::roots(pq_schema, root_indices);
        let reader = CoalescingAsyncReader::new(
            self.open_async().await?,
            counters.clone(),
            COALESCE_GAP_BYTES,
        );
        let builder = ParquetRecordBatchStreamBuilder::new_with_metadata(reader, self.meta.clone())
            .with_projection(mask)
            .with_row_selection(selection_from_offsets(offsets))
            .with_batch_size(8192);
        let proj_schema = builder.schema().clone();
        let mut stream = builder
            .build()
            .map_err(|e| DataFusionError::Execution(format!("build payload stream: {e}")))?;
        let mut taken: Vec<RecordBatch> = Vec::new();
        while let Some(b) = stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read payload batch: {e}")))?
        {
            taken.push(b);
        }
        if taken.is_empty() {
            Ok(RecordBatch::new_empty(proj_schema))
        } else {
            concat_batches(&taken[0].schema(), &taken)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
        }
    }

    async fn open_async(&self) -> Result<tokio::fs::File> {
        tokio::fs::File::open(&self.path).await.map_err(|e| {
            DataFusionError::Execution(format!("reopen parquet '{}': {e}", self.path.display()))
        })
    }
}

/// Distinct requested keys present in the result batch.
fn distinct_present(batch: &RecordBatch, probe_set: &HashSet<u64>) -> Result<Vec<u64>> {
    let key = batch
        .column_by_name("key")
        .and_then(|c| c.as_any().downcast_ref::<UInt64Array>())
        .ok_or_else(|| DataFusionError::Execution("sift batch missing UInt64 key".into()))?;
    let mut seen: HashSet<u64> = HashSet::new();
    let mut out = Vec::new();
    for &v in key.values() {
        if probe_set.contains(&v) && seen.insert(v) {
            out.push(v);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{BinaryArray, UInt64Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    fn sift_batch(keys: &[u64]) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("key", DataType::UInt64, false),
            Field::new("sift", DataType::Binary, false),
            Field::new("poly", DataType::Binary, false),
        ]));
        let sift: Vec<Vec<u8>> = keys.iter().map(|k| vec![(k & 0xff) as u8, 0x01]).collect();
        let poly: Vec<Vec<u8>> = keys.iter().map(|k| vec![(k & 0xff) as u8, 0x02]).collect();
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(keys.to_vec())),
                Arc::new(BinaryArray::from_iter_values(
                    sift.iter().map(|v| v.as_slice()),
                )),
                Arc::new(BinaryArray::from_iter_values(
                    poly.iter().map(|v| v.as_slice()),
                )),
            ],
        )
        .unwrap()
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn sift_lookup_round_trips_keys_across_pages() {
        // 3000 strictly-increasing u64 keys spanning many 512-row pages, incl. a
        // large `uid` so the high 32 bits are exercised.
        let n = 3000u64;
        let keys: Vec<u64> = (0..n).map(|i| ((i + 100) << 32) | (i * 3 + 7)).collect();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        assert_eq!(
            write_sift_parquet(tmp.path(), &[sift_batch(&keys)]).unwrap(),
            n as usize
        );

        let lookup = SinglePathParquetSiftLookup::open(
            tmp.path(),
            vec!["key".to_string(), "sift".to_string(), "poly".to_string()],
        )
        .await
        .unwrap();

        // Probe a mix incl. first/last/middle and one absent key.
        let probes = vec![keys[0], keys[1], keys[1500], keys[2999], u64::MAX];
        let (batch, present) = lookup.take_keys(&probes).await.unwrap();
        assert_eq!(batch.num_rows(), 4);
        assert_eq!(present.len(), 4);

        let out_key = batch
            .column_by_name("key")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();
        let out_sift = batch
            .column_by_name("sift")
            .unwrap()
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap();
        for r in 0..batch.num_rows() {
            let k = out_key.value(r);
            assert!(probes.contains(&k));
            // sift blob encodes (k & 0xff, 0x01) per the fixture.
            assert_eq!(out_sift.value(r), &[(k & 0xff) as u8, 0x01]);
        }
    }
}
