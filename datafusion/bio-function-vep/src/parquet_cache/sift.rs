//! Parquet point-lookup for the `translation_sift` cache entity.
//!
//! Position-sliced SIFT/PolyPhen blobs keyed by a `u64` `(transcript_uid<<32) |
//! position`. Exposes `open` + a `take_keys(&[u64]) -> (RecordBatch, present)`
//! and resolves keys through the
//! footer [`PageDir`] (u64 key column) + [`CoalescingAsyncReader`] instead of an
//! in-memory BTree. The shard is written no-dictionary, page-indexed, and
//! physically sorted by `key` (a single monotone run — SIFT keys are unique), so
//! the reader is the same three-phase loop as the variation lookup: resolve
//! candidate pages → exact key offsets → projected payload take.
//!
//! The returned batch (`key: UInt64`, `sift: Binary`, `poly: Binary`) is the same
//! shape the Parquet path yields, so `position_predictions_from_batch` decodes it
//! unchanged.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

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
pub(crate) mod decoded_pages;

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
    decoded_cache: Option<decoded_pages::DecodedPageCache>,
    profile: bool,
    point_loads: AtomicU64,
    requested_keys: AtomicU64,
    decoded_key_rows: AtomicU64,
    physical_bytes: AtomicU64,
    read_ops: AtomicU64,
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
            decoded_cache: None,
            profile: std::env::var_os("VEP_ENGINE_PROFILE").is_some_and(|value| value == "1"),
            point_loads: AtomicU64::new(0),
            requested_keys: AtomicU64::new(0),
            decoded_key_rows: AtomicU64::new(0),
            physical_bytes: AtomicU64::new(0),
            read_ops: AtomicU64::new(0),
            path: path.to_path_buf(),
            meta,
            page_dir,
            key_leaf,
            projection,
        })
    }

    /// Resolve `keys` to their rows. Returns the payload batch (`key`/`sift`/`poly`
    /// for the matched keys, in `key` order) and the distinct keys found — the
    /// same `(RecordBatch, present)` contract as the Parquet `take_keys`.
    pub async fn take_keys(&self, keys: &[u64]) -> Result<(RecordBatch, Vec<u64>)> {
        if let Some(result) = self.take_keys_cached(keys).await? {
            return Ok(result);
        }
        self.take_keys_uncached(keys).await
    }

    async fn take_keys_uncached(&self, keys: &[u64]) -> Result<(RecordBatch, Vec<u64>)> {
        let mut sorted: Vec<u64> = keys.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        let probe_set: HashSet<u64> = sorted.iter().copied().collect();
        let counters = IoCounters::new();

        let ranges = self.page_dir.resolve_ranges(&sorted);
        let offsets = self.exact_offsets(&ranges, &probe_set, &counters).await?;
        let batch = self.take_payload(&offsets, &counters).await?;
        let present = distinct_present(&batch, &probe_set)?;
        if self.profile {
            let (bytes, ops) = counters.snapshot();
            self.point_loads.fetch_add(1, Ordering::Relaxed);
            self.requested_keys
                .fetch_add(sorted.len() as u64, Ordering::Relaxed);
            self.physical_bytes.fetch_add(bytes, Ordering::Relaxed);
            self.read_ops.fetch_add(ops, Ordering::Relaxed);
        }
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
            if self.profile {
                self.decoded_key_rows
                    .fetch_add(b.num_rows() as u64, Ordering::Relaxed);
            }
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

impl Drop for SinglePathParquetSiftLookup {
    fn drop(&mut self) {
        if self.profile {
            eprintln!(
                "[VEP_SIFT_LOOKUP] point_loads={} requested_keys={} decoded_key_rows={} physical_bytes={} read_ops={}",
                self.point_loads.load(Ordering::Relaxed),
                self.requested_keys.load(Ordering::Relaxed),
                self.decoded_key_rows.load(Ordering::Relaxed),
                self.physical_bytes.load(Ordering::Relaxed),
                self.read_ops.load(Ordering::Relaxed)
            );
        }
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

impl SinglePathParquetSiftLookup {
    pub(crate) fn with_page_cache(mut self, budget: Arc<decoded_pages::PageCacheBudget>) -> Self {
        self.decoded_cache = Some(decoded_pages::DecodedPageCache::new(budget));
        self
    }

    async fn take_keys_cached(&self, keys: &[u64]) -> Result<Option<(RecordBatch, Vec<u64>)>> {
        let Some(cache) = &self.decoded_cache else {
            return Ok(None);
        };
        let mut sorted = keys.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        let ranges = self.page_dir.resolve_page_ranges(&sorted);
        if ranges.is_empty() {
            return Ok(None);
        }
        let requested: HashSet<u64> = sorted.iter().copied().collect();
        let mut pages = Vec::new();
        let mut selected: Vec<(usize, usize)> = Vec::new();
        for range in ranges {
            let Some(page) = cache.get(self, range).await? else {
                // Admission failure is not a negative lookup. Re-run the original
                // exact point reader, including any overlapping key pages.
                return Ok(None);
            };
            let ka = page
                .batch
                .column_by_name("key")
                .and_then(|a| a.as_any().downcast_ref::<UInt64Array>())
                .ok_or_else(|| {
                    DataFusionError::Execution("cached SIFT page has no UInt64 key".into())
                })?;
            let page_index = pages.len();
            let values = ka.values();
            if page.keys_sorted && !values.is_empty() {
                // Sorted pages need only the requested keys in their extent,
                // not a 512-row hash probe for every repeated single-position hit.
                let first = sorted.partition_point(|key| *key < values[0]);
                let end = sorted.partition_point(|key| *key <= values[values.len() - 1]);
                for key in &sorted[first..end] {
                    let lo = values.partition_point(|value| value < key);
                    let hi = values.partition_point(|value| value <= key);
                    selected.extend((lo..hi).map(|row| (page_index, row)));
                }
            } else {
                // Legacy unsorted pages preserve physical order and duplicates.
                selected.extend(
                    values.iter().enumerate().filter_map(|(row, key)| {
                        requested.contains(key).then_some((page_index, row))
                    }),
                );
            }
            // Keep leases pinned until the independently owned result is built.
            pages.push(page);
        }
        let schema = self.payload_schema()?;
        let batch = if selected.is_empty() {
            RecordBatch::new_empty(schema)
        } else {
            use datafusion::arrow::array::{MutableArrayData, make_array};
            let columns = (0..schema.fields().len())
                .map(|column| {
                    let data: Vec<_> = pages
                        .iter()
                        .map(|page| page.batch.column(column).to_data())
                        .collect();
                    let mut output =
                        MutableArrayData::new(data.iter().collect(), false, selected.len());
                    for &(page, row) in &selected {
                        output.try_extend(page, row, row + 1)?;
                    }
                    Ok(make_array(output.freeze()))
                })
                .collect::<Result<Vec<_>, datafusion::arrow::error::ArrowError>>()
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            RecordBatch::try_new(schema, columns)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
        };
        let present = distinct_present(&batch, &requested)?;
        Ok(Some((batch, present)))
    }

    fn payload_schema(&self) -> Result<Arc<datafusion::arrow::datatypes::Schema>> {
        let schema = self.meta.schema();
        let mut roots: Vec<_> = self
            .projection
            .iter()
            .filter_map(|n| schema.index_of(n).ok())
            .collect();
        roots.sort_unstable();
        roots.dedup();
        schema
            .project(&roots)
            .map(Arc::new)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }

    /// Upper bounds from the footer's per-page *unencoded* byte counts. If an
    /// older file lacks these counts, or uses dictionaries/nested payloads,
    /// keep the original bounded point reader instead of guessing an expansion.
    fn decoded_page_bound(&self, range: (u64, u64)) -> Option<(usize, usize)> {
        use datafusion::arrow::datatypes::DataType;
        let metadata = self.meta.metadata();
        let offsets = metadata.offset_index()?;
        let mut raw = 0usize;
        let mut rows = 0usize;
        let mut compressed = 0usize;
        let mut rg_start = 0u64;
        for rg_index in 0..metadata.num_row_groups() {
            let rg = metadata.row_group(rg_index);
            let rg_rows = u64::try_from(rg.num_rows()).ok()?;
            let rg_end = rg_start.checked_add(rg_rows)?;
            if range.0 < rg_end && range.1 > rg_start {
                for column in &self.projection {
                    let index = self.meta.schema().index_of(column).ok()?;
                    let column_meta = rg.column(index);
                    if column_meta.dictionary_page_offset().is_some() {
                        return None;
                    }
                    let offset = offsets.get(rg_index)?.get(index)?;
                    let locations = offset.page_locations();
                    for (page, location) in locations.iter().enumerate() {
                        let first =
                            rg_start.checked_add(u64::try_from(location.first_row_index).ok()?)?;
                        let end = match locations.get(page + 1) {
                            Some(next) => {
                                rg_start.checked_add(u64::try_from(next.first_row_index).ok()?)?
                            }
                            None => rg_end,
                        };
                        if first >= range.1 || end <= range.0 {
                            continue;
                        }
                        let page_rows = usize::try_from(end.checked_sub(first)?).ok()?;
                        let bytes = match self.meta.schema().field(index).data_type() {
                            DataType::UInt64 => page_rows.checked_mul(8)?,
                            DataType::Binary | DataType::LargeBinary => usize::try_from(
                                *offset.unencoded_byte_array_data_bytes()?.get(page)?,
                            )
                            .ok()?,
                            _ => return None,
                        };
                        raw = raw.checked_add(bytes)?;
                        rows = rows.checked_add(page_rows)?;
                        compressed = compressed
                            .checked_add(usize::try_from(location.compressed_page_size).ok()?)?;
                    }
                }
            }
            rg_start = rg_end;
        }
        if rows == 0 {
            return None;
        }
        // Cover Arrow offsets/validity/alignment and a concatenation copy. The
        // scratch reservation additionally covers compressed reads and decoder
        // working buffers, and is checked before any bulk payload read.
        let resident = raw
            .checked_mul(2)?
            .checked_add(rows.checked_mul(64)?)?
            .checked_add(64 * 1024)?;
        let scratch = raw
            .checked_mul(4)?
            .checked_add(rows.checked_mul(128)?)?
            .checked_add(compressed)?
            .checked_add(4 * 1024 * 1024)?;
        Some((resident, scratch))
    }

    async fn read_page_range(
        &self,
        range: (u64, u64),
        counters: &IoCounters,
    ) -> Result<RecordBatch> {
        let roots: Vec<_> = self
            .projection
            .iter()
            .filter_map(|n| self.meta.schema().index_of(n).ok())
            .collect();
        let reader = CoalescingAsyncReader::new(
            self.open_async().await?,
            counters.clone(),
            COALESCE_GAP_BYTES,
        );
        let mut stream =
            ParquetRecordBatchStreamBuilder::new_with_metadata(reader, self.meta.clone())
                .with_projection(ProjectionMask::roots(self.meta.parquet_schema(), roots))
                .with_row_selection(selection_from_ranges(&[range]))
                .with_batch_size(8192)
                .build()
                .map_err(|e| DataFusionError::Execution(format!("build SIFT page stream: {e}")))?;
        let mut batches = Vec::new();
        while let Some(batch) = stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read SIFT page: {e}")))?
        {
            batches.push(batch);
        }
        let schema = self.payload_schema()?;
        if batches.len() == 1 {
            return Ok(batches.pop().unwrap());
        }
        concat_batches(&schema, &batches)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
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
