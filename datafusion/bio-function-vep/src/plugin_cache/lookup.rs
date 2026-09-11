//! Buffer-batched, page-scoped plugin lookup (spec §5).
//!
//! Per annotation buffer, one PageDir take reads **only** the candidate pages for
//! that buffer's positions — never the whole shard ([`PluginLookup::take_buffer`]).
//! The resulting compact batch becomes a [`PluginBufferSlice`], against which the
//! transcript engine runs synchronous per-transcript [`PluginBufferSlice::probe`]s
//! (filter by `allele_string` + per-transcript match discriminator). Built on the
//! shared `parquet_cache::page_dir` primitives; `SinglePathParquetVariationLookup`
//! is untouched.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use std::sync::Arc;

use datafusion::arrow::array::{
    Array, Float32Array, Int32Array, LargeStringArray, StringArray, StringViewArray, UInt32Array,
};
use datafusion::arrow::datatypes::{Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};
use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;

use coitrees::{COITree, GenericInterval, Interval, IntervalTree};

use crate::parquet_cache::page_dir::{
    CoalescingAsyncReader, IoCounters, PageDir, selection_from_offsets, selection_from_ranges,
};

/// Byte gap under which nearby page ranges coalesce into one read.
const COALESCE_GAP_BYTES: u64 = 512 * 1024;

/// One decoded plugin value.
#[derive(Debug, Clone, PartialEq)]
pub enum PluginScalar {
    Str(String),
    F32(f32),
    I32(i32),
    Null,
}

/// A per-chrom plugin shard handle: footer metadata + a [`PageDir`] over the
/// `start` leaf. Opening reads no data; each buffer's rows are read on demand by
/// [`Self::take_buffer`].
pub struct PluginLookup {
    path: PathBuf,
    meta: ArrowReaderMetadata,
    page_dir: PageDir,
    start_leaf: usize,
    /// Payload columns read per take: `allele_string`, `<match cols…>`,
    /// `<value cols…>` (`start` is always prepended at take time).
    projection: Vec<String>,
    n_match: usize,
    n_values: usize,
}

impl PluginLookup {
    /// Open a plugin shard: load page-indexed metadata and build the `start`
    /// [`PageDir`]. `match_columns` then `value_columns` form the payload
    /// projection.
    pub async fn open(
        shard: &Path,
        match_columns: Vec<String>,
        value_columns: Vec<String>,
    ) -> Result<Self> {
        let file = std::fs::File::open(shard).map_err(|e| {
            DataFusionError::Execution(format!("open plugin shard '{}': {e}", shard.display()))
        })?;
        // `with_page_index_policy` is the non-deprecated API, but `PageIndexPolicy`
        // is private in parquet 58, so the deprecated setter is the only usable one.
        #[allow(deprecated)]
        let meta =
            ArrowReaderMetadata::load(&file, ArrowReaderOptions::new().with_page_index(true))
                .map_err(|e| DataFusionError::Execution(format!("load plugin metadata: {e}")))?;
        let start_leaf = meta
            .parquet_schema()
            .columns()
            .iter()
            .position(|c| c.name() == "start")
            .ok_or_else(|| {
                DataFusionError::Execution("plugin shard has no 'start' column".into())
            })?;
        let page_dir = PageDir::build(&meta, start_leaf)?;
        let n_match = match_columns.len();
        let n_values = value_columns.len();
        let mut projection = vec!["allele_string".to_string()];
        projection.extend(match_columns);
        projection.extend(value_columns);
        Ok(Self {
            path: shard.to_path_buf(),
            meta,
            page_dir,
            start_leaf,
            projection,
            n_match,
            n_values,
        })
    }

    /// Number of match discriminator columns / value columns in the payload.
    pub fn n_match(&self) -> usize {
        self.n_match
    }
    pub fn n_values(&self) -> usize {
        self.n_values
    }

    async fn open_async(&self) -> Result<tokio::fs::File> {
        tokio::fs::File::open(&self.path).await.map_err(|e| {
            DataFusionError::Execution(format!(
                "reopen plugin shard '{}': {e}",
                self.path.display()
            ))
        })
    }

    /// The projected payload column order: `start` + [`Self::projection`].
    fn payload_columns(&self) -> Vec<String> {
        let mut cols = vec!["start".to_string()];
        cols.extend(self.projection.iter().cloned());
        cols
    }

    fn payload_mask(&self) -> ProjectionMask {
        let arrow_schema = self.meta.schema();
        let root_indices: Vec<usize> = self
            .payload_columns()
            .iter()
            .filter_map(|n| arrow_schema.index_of(n).ok())
            .collect();
        ProjectionMask::roots(self.meta.parquet_schema(), root_indices)
    }

    /// The Arrow schema of a take batch: the payload columns in ascending file
    /// order (which is `[start, allele_string, match…, value…]` by shard layout).
    /// Derived from the file schema, not `builder.schema()` (which returns the
    /// full, unprojected schema in parquet 58).
    fn projected_schema(&self) -> SchemaRef {
        let full = self.meta.schema();
        let mut fields: Vec<_> = self
            .payload_columns()
            .iter()
            .filter_map(|n| full.index_of(n).ok())
            .collect();
        fields.sort_unstable();
        Arc::new(Schema::new(
            fields
                .into_iter()
                .map(|i| full.field(i).clone())
                .collect::<Vec<_>>(),
        ))
    }

    /// Read the candidate rows for `sorted_unique_starts` — only the pages those
    /// positions fall on. Returns a batch with columns
    /// `[start, allele_string, <match cols…>, <value cols…>]`.
    pub async fn take_buffer(&self, sorted_unique_starts: &[u32]) -> Result<RecordBatch> {
        if sorted_unique_starts.is_empty() {
            return Ok(RecordBatch::new_empty(self.projected_schema()));
        }

        let counters = IoCounters::new();
        let probe_set: HashSet<u32> = sorted_unique_starts.iter().copied().collect();
        let probes64: Vec<u64> = sorted_unique_starts.iter().map(|&s| s as u64).collect();
        let ranges = self.page_dir.resolve_ranges(&probes64);

        // Phase 2: start-only read over candidate pages → exact row offsets.
        let start_mask = ProjectionMask::leaves(self.meta.parquet_schema(), [self.start_leaf]);
        let mut start_stream = ParquetRecordBatchStreamBuilder::new_with_metadata(
            CoalescingAsyncReader::new(
                self.open_async().await?,
                counters.clone(),
                COALESCE_GAP_BYTES,
            ),
            self.meta.clone(),
        )
        .with_projection(start_mask)
        .with_row_selection(selection_from_ranges(&ranges))
        .with_batch_size(8192)
        .build()
        .map_err(|e| DataFusionError::Execution(format!("build start stream: {e}")))?;
        let mut offsets = Vec::new();
        let mut cursor = ranges.iter().flat_map(|&(a, b)| a..b);
        while let Some(b) = start_stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read start batch: {e}")))?
        {
            let sa = b
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("start column must be UInt32".into()))?;
            for &v in sa.values() {
                let off = cursor.next().ok_or_else(|| {
                    DataFusionError::Execution("row range cursor underflow".into())
                })?;
                if probe_set.contains(&v) {
                    offsets.push(off);
                }
            }
        }

        // Phase 3: projected payload take at the exact offsets.
        let builder = ParquetRecordBatchStreamBuilder::new_with_metadata(
            CoalescingAsyncReader::new(
                self.open_async().await?,
                counters.clone(),
                COALESCE_GAP_BYTES,
            ),
            self.meta.clone(),
        )
        .with_projection(self.payload_mask())
        .with_row_selection(selection_from_offsets(&offsets))
        .with_batch_size(8192);
        let proj_schema = self.projected_schema();
        let mut payload_stream = builder
            .build()
            .map_err(|e| DataFusionError::Execution(format!("build payload stream: {e}")))?;
        let mut taken: Vec<RecordBatch> = Vec::new();
        while let Some(b) = payload_stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read payload batch: {e}")))?
        {
            taken.push(b);
        }
        if taken.is_empty() {
            return Ok(RecordBatch::new_empty(proj_schema));
        }
        datafusion::arrow::compute::concat_batches(&taken[0].schema(), &taken)
            .map_err(|e| DataFusionError::Execution(format!("concat payload: {e}")))
    }
}

/// Lookup key: `(start, allele_string, <match discriminator values…>)`.
type SliceKey = (u32, String, Vec<Option<String>>);

/// The in-memory working set for one buffer × one plugin, built from a
/// [`PluginLookup::take_buffer`] batch. Sync-probed per transcript.
pub struct PluginBufferSlice {
    rows: HashMap<SliceKey, Vec<PluginScalar>>,
}

impl PluginBufferSlice {
    /// Build the slice from a take batch. `n_match`/`n_values` give the payload
    /// layout (`[start, allele_string, match×n_match, value×n_values]`).
    pub fn from_batch(batch: &RecordBatch, n_match: usize, n_values: usize) -> Result<Self> {
        let mut rows = HashMap::new();
        if batch.num_rows() == 0 {
            return Ok(Self { rows });
        }
        let start = batch
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| DataFusionError::Execution("start column must be UInt32".into()))?;
        for r in 0..batch.num_rows() {
            let allele = string_value(batch.column(1).as_ref(), r)?
                .ok_or_else(|| DataFusionError::Execution("null allele_string".into()))?;
            let mut match_vals = Vec::with_capacity(n_match);
            for c in 0..n_match {
                match_vals.push(string_value(batch.column(2 + c).as_ref(), r)?);
            }
            let mut values = Vec::with_capacity(n_values);
            for c in 0..n_values {
                values.push(decode_scalar(batch.column(2 + n_match + c).as_ref(), r)?);
            }
            rows.insert((start.value(r), allele, match_vals), values);
        }
        Ok(Self { rows })
    }

    /// Probe for `(start, allele_string, match_values)`; `None` on any miss.
    /// A `None` element of `match_values` (e.g. a non-missense consequence with no
    /// amino-acid change) will not match a stored `Some` discriminator — the
    /// per-transcript gate.
    pub fn probe(
        &self,
        start: u32,
        allele_string: &str,
        match_values: &[Option<String>],
    ) -> Option<Vec<PluginScalar>> {
        self.rows
            .get(&(start, allele_string.to_string(), match_values.to_vec()))
            .cloned()
    }
}

/// One interval row's values, in shard column order. Its position in
/// `IntervalLookup::rows` is the file ordinal, which is also the COITree
/// metadata, so the tie-break "first record in file order" is a `min` over
/// the ordinals the tree reports.
struct IntervalRow {
    values: Vec<PluginScalar>,
}

/// Whole-shard interval lookup for `LookupKind::Interval` plugins: one
/// `COITree` per discriminator tuple over closed 1-based spans, read once at
/// open and probed synchronously. `probe` returns the overlapping row with the
/// smallest shard ordinal. The shard is written in `(start, arrival)` order
/// (`build::write_interval_shard`), which is the iteration order of the
/// position-sorted tabix file an Ensembl plugin reads, so this is the record
/// such a plugin takes. A plugin with no match columns has a single tree.
pub struct IntervalLookup {
    rows: Vec<IntervalRow>,
    trees: HashMap<Vec<Option<String>>, COITree<usize, u32>>,
    n_values: usize,
}

impl IntervalLookup {
    /// Read the whole shard (`start, end, <match…>, <value…>`) and index it.
    pub async fn open(
        shard: &Path,
        match_columns: Vec<String>,
        value_columns: Vec<String>,
    ) -> Result<Self> {
        let file = tokio::fs::File::open(shard).await.map_err(|e| {
            DataFusionError::Execution(format!("open plugin shard '{}': {e}", shard.display()))
        })?;
        let builder = ParquetRecordBatchStreamBuilder::new(file)
            .await
            .map_err(|e| DataFusionError::Execution(format!("open interval shard: {e}")))?;
        let arrow_schema = builder.schema().clone();
        let mut wanted = vec!["start".to_string(), "end".to_string()];
        wanted.extend(match_columns.iter().cloned());
        wanted.extend(value_columns.iter().cloned());
        let roots: Vec<usize> = wanted
            .iter()
            .map(|n| {
                arrow_schema.index_of(n).map_err(|e| {
                    DataFusionError::Execution(format!("interval shard column {n}: {e}"))
                })
            })
            .collect::<Result<_>>()?;
        let mask = ProjectionMask::roots(builder.parquet_schema(), roots);
        let mut stream = builder
            .with_projection(mask)
            .with_batch_size(8192)
            .build()
            .map_err(|e| DataFusionError::Execution(format!("build interval stream: {e}")))?;
        let n_match = match_columns.len();
        let n_values = value_columns.len();
        let mut rows: Vec<IntervalRow> = Vec::new();
        let mut intervals: HashMap<Vec<Option<String>>, Vec<Interval<usize>>> = HashMap::new();
        while let Some(b) = stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read interval batch: {e}")))?
        {
            let schema = b.schema();
            let col = |name: &str| -> Result<usize> {
                schema
                    .index_of(name)
                    .map_err(|e| DataFusionError::Execution(e.to_string()))
            };
            let start = b
                .column(col("start")?)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("start column must be UInt32".into()))?;
            let end = b
                .column(col("end")?)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("end column must be UInt32".into()))?;
            let match_idx: Vec<usize> = match_columns
                .iter()
                .map(|n| col(n))
                .collect::<Result<_>>()?;
            let value_idx: Vec<usize> = value_columns
                .iter()
                .map(|n| col(n))
                .collect::<Result<_>>()?;
            for r in 0..b.num_rows() {
                let mut key = Vec::with_capacity(n_match);
                for &i in &match_idx {
                    key.push(string_value(b.column(i).as_ref(), r)?);
                }
                let mut values = Vec::with_capacity(n_values);
                for &i in &value_idx {
                    values.push(decode_scalar(b.column(i).as_ref(), r)?);
                }
                let ordinal = rows.len();
                rows.push(IntervalRow { values });
                // Closed 1-based span; positions fit i32 for every contig.
                intervals.entry(key).or_default().push(Interval::new(
                    i32::try_from(start.value(r)).unwrap_or(i32::MAX),
                    i32::try_from(end.value(r)).unwrap_or(i32::MAX),
                    ordinal,
                ));
            }
        }
        let trees = intervals
            .into_iter()
            .map(|(key, ivs)| (key, COITree::new(&ivs)))
            .collect();
        Ok(Self {
            rows,
            trees,
            n_values,
        })
    }

    pub fn n_values(&self) -> usize {
        self.n_values
    }

    /// First row (file order) whose inclusive span overlaps `[span_start, span_end]`
    /// under the given discriminators; `None` on a miss or a `None` discriminator
    /// that the shard stores as `Some`.
    pub fn probe(
        &self,
        span_start: u32,
        span_end: u32,
        match_values: &[Option<String>],
    ) -> Option<&[PluginScalar]> {
        let tree = self.trees.get(match_values)?;
        let (first, last) = (
            i32::try_from(span_start).unwrap_or(i32::MAX),
            i32::try_from(span_end).unwrap_or(i32::MAX),
        );
        let mut best: Option<usize> = None;
        tree.query(first, last, |node| {
            let ordinal = *GenericInterval::<usize>::metadata(node);
            best = Some(best.map_or(ordinal, |b| b.min(ordinal)));
        });
        best.map(|i| self.rows[i].values.as_slice())
    }
}

/// Test fixtures shared with the registry tests.
#[cfg(test)]
pub(crate) mod test_support {
    use super::*;
    use crate::plugin_cache::cache_manifest::LookupKind;
    use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};
    use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};

    /// Three rows for ENSG1 in file order wide, narrow, late (all cover 160;
    /// a tree visits them in its own order, so the ordinal tie-break is what
    /// makes "wide" win) and one row for ENSG2 with no value.
    pub(crate) fn write_interval_shard(path: &std::path::Path) {
        let matches = vec![MatchColumn {
            column: "gene_id".into(),
            template: "{Gene}".into(),
        }];
        let vals = vec![ValueColumn {
            column: "rat".into(),
            csq_field: "PO_Rat".into(),
            ty: ValueType::Utf8,
            description: None,
        }];
        let schema = plugin_output_schema(LookupKind::Interval, &matches, &vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 150, 155, 400])),
                Arc::new(UInt32Array::from(vec![300u32, 200, 350, 500])),
                Arc::new(StringArray::from(vec!["ENSG1", "ENSG1", "ENSG1", "ENSG2"])),
                Arc::new(StringArray::from(vec![
                    Some("wide"),
                    Some("narrow"),
                    Some("late"),
                    None,
                ])),
                Arc::new(Int8Array::from(vec![1i8, 1, 1, 1])),
            ],
        )
        .unwrap();
        let mut w = PluginShardWriter::create(path, schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();
    }
}

/// Read a string cell across the Utf8 / LargeUtf8 / Utf8View encodings
/// DataFusion may materialize.
fn string_value(col: &dyn Array, row: usize) -> Result<Option<String>> {
    if col.is_null(row) {
        return Ok(None);
    }
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return Ok(Some(a.value(row).to_string()));
    }
    if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
        return Ok(Some(a.value(row).to_string()));
    }
    if let Some(a) = col.as_any().downcast_ref::<LargeStringArray>() {
        return Ok(Some(a.value(row).to_string()));
    }
    Err(DataFusionError::Execution(format!(
        "expected a string column, got {}",
        col.data_type()
    )))
}

/// Decode one cell into a [`PluginScalar`] based on the column's Arrow type.
fn decode_scalar(col: &dyn Array, row: usize) -> Result<PluginScalar> {
    if col.is_null(row) {
        return Ok(PluginScalar::Null);
    }
    if let Some(a) = col.as_any().downcast_ref::<Float32Array>() {
        return Ok(PluginScalar::F32(a.value(row)));
    }
    if let Some(a) = col.as_any().downcast_ref::<Int32Array>() {
        return Ok(PluginScalar::I32(a.value(row)));
    }
    if let Some(s) = string_value(col, row)? {
        return Ok(PluginScalar::Str(s));
    }
    Ok(PluginScalar::Null)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::cache_manifest::LookupKind;
    use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};
    use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    /// Three positions; 100 has two protein_variants (multi-isoform).
    fn write_shard(path: &std::path::Path) {
        let matches = vec![MatchColumn {
            column: "protein_variant".into(),
            template: "{ref_aa}{Protein_position}{alt_aa}".into(),
        }];
        let vals = vec![ValueColumn {
            column: "am".into(),
            csq_field: "am".into(),
            ty: ValueType::Float32,
            description: None,
        }];
        let schema = plugin_output_schema(LookupKind::Point, &matches, &vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["22", "22", "22"])),
                Arc::new(UInt32Array::from(vec![100u32, 100, 200])),
                Arc::new(UInt32Array::from(vec![100u32, 100, 200])),
                Arc::new(StringArray::from(vec!["A/G", "A/G", "C/T"])),
                Arc::new(StringArray::from(vec!["R12G", "R78G", "S9L"])),
                Arc::new(Float32Array::from(vec![0.0392f32, 0.0427, 0.5])),
                Arc::new(Int8Array::from(vec![1i8, 1, 1])),
            ],
        )
        .unwrap();
        let mut w = PluginShardWriter::create(path, schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn take_buffer_reads_candidate_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr22.parquet");
        write_shard(&path);
        let lk = PluginLookup::open(&path, vec!["protein_variant".into()], vec!["am".into()])
            .await
            .unwrap();
        // buffer requests only position 100 → both isoform rows; 200 excluded
        let batch = lk.take_buffer(&[100]).await.unwrap();
        assert_eq!(batch.num_rows(), 2);
        let names: Vec<_> = batch
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();
        assert_eq!(
            names,
            vec!["start", "allele_string", "protein_variant", "am"]
        );
        // empty request → empty batch
        assert_eq!(lk.take_buffer(&[]).await.unwrap().num_rows(), 0);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn interval_probe_overlaps_and_gates_on_discriminator() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.parquet");
        super::test_support::write_interval_shard(&path);
        let lk = IntervalLookup::open(&path, vec!["gene_id".into()], vec!["rat".into()])
            .await
            .unwrap();
        let g1 = [Some("ENSG1".to_string())];
        let g2 = [Some("ENSG2".to_string())];
        // inside all three ENSG1 spans → first in FILE order, whatever the tree visits first
        assert_eq!(
            lk.probe(160, 160, &g1).unwrap(),
            &[PluginScalar::Str("wide".into())]
        );
        // only "late" covers 320..350
        assert_eq!(
            lk.probe(320, 320, &g1).unwrap(),
            &[PluginScalar::Str("late".into())]
        );
        // boundaries are inclusive
        assert_eq!(
            lk.probe(300, 300, &g1).unwrap(),
            &[PluginScalar::Str("wide".into())]
        );
        assert_eq!(
            lk.probe(99, 100, &g1).unwrap(),
            &[PluginScalar::Str("wide".into())]
        );
        assert!(lk.probe(351, 351, &g1).is_none());
        // an insertion span [pos, pos+1] straddling a start boundary
        assert_eq!(lk.probe(399, 400, &g2).unwrap(), &[PluginScalar::Null]);
        // discriminator gate
        assert!(lk.probe(160, 160, &[Some("ENSG9".to_string())]).is_none());
        assert!(lk.probe(160, 160, &[None]).is_none());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn slice_probe_gates_on_discriminator() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr22.parquet");
        write_shard(&path);
        let lk = PluginLookup::open(&path, vec!["protein_variant".into()], vec!["am".into()])
            .await
            .unwrap();
        let batch = lk.take_buffer(&[100]).await.unwrap();
        let slice = PluginBufferSlice::from_batch(&batch, 1, 1).unwrap();
        // isoform-specific hits
        match slice.probe(100, "A/G", &[Some("R78G".into())]).unwrap()[0] {
            PluginScalar::F32(v) => assert!((v - 0.0427).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        match slice.probe(100, "A/G", &[Some("R12G".into())]).unwrap()[0] {
            PluginScalar::F32(v) => assert!((v - 0.0392).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        // gate: no aa-change (non-missense) → miss
        assert!(slice.probe(100, "A/G", &[None]).is_none());
        // wrong allele → miss
        assert!(slice.probe(100, "C/T", &[Some("R12G".into())]).is_none());
        // unknown aa-change → miss
        assert!(slice.probe(100, "A/G", &[Some("X9Y".into())]).is_none());
    }
}
