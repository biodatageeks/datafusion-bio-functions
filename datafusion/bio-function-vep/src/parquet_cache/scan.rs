//! Dictionary-enabled Parquet for the context (scan) cache entities:
//! `transcript`, `exon`, `regulatory`, `motif`, `translation_core`.
//!
//! Unlike the point-lookup entities (variation, translation_sift), these are read
//! by a **full projected scan** — every row, a subset of columns — so they use the
//! opposite writer profile: **dictionary enabled** (free for scans, and needed for
//! `translation_core`'s big sequences + `list<struct>` to match Lance). No encode
//! / decode step: the rows are written verbatim (same Arrow schema as the Lance
//! path) and read back unchanged, so the format-agnostic per-entity parsers are
//! untouched.
//!
//! [`read_context_parquet`] mirrors `scan_projected_existing_columns`: a full
//! read, projected to the requested columns that exist (missing columns tolerated,
//! quotes stripped), returning `Vec<RecordBatch>` in file (write) order.

use std::collections::HashSet;
use std::path::Path;

use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use parquet::arrow::ArrowWriter;
use parquet::arrow::ProjectionMask;
use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;

/// Row-group row cap for context shards. Bounded so a row group of big
/// `translation_core`/`transcript` sequence rows stays a sane size in memory
/// during the write (the writer buffers a full row group before flushing).
const CONTEXT_ROW_GROUP_ROWS: usize = 50_000;

/// Dictionary-enabled `WriterProperties` for a context (scan) shard: zstd(3),
/// dictionaries on (free for full scans; closes the `translation_core` gap vs
/// Lance), bounded row groups. No declared sort order — the read path is a full
/// scan, and the physical order already follows the source query.
pub fn context_writer_properties() -> WriterProperties {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(3).unwrap()))
        .set_dictionary_enabled(true)
        .set_max_row_group_row_count(Some(CONTEXT_ROW_GROUP_ROWS))
        .build()
}

/// Streaming writer for one context shard. Opens the `ArrowWriter` with
/// [`context_writer_properties`] once and appends batches in source order.
pub struct ContextParquetShardWriter {
    writer: ArrowWriter<std::fs::File>,
    rows: usize,
}

impl ContextParquetShardWriter {
    pub fn create(path: &Path, schema: SchemaRef) -> Result<Self> {
        let file = std::fs::File::create(path).map_err(|e| {
            DataFusionError::Execution(format!("create parquet '{}': {e}", path.display()))
        })?;
        let writer = ArrowWriter::try_new(file, schema, Some(context_writer_properties()))
            .map_err(|e| DataFusionError::Execution(format!("open ArrowWriter: {e}")))?;
        Ok(Self { writer, rows: 0 })
    }

    pub fn write(&mut self, batch: &RecordBatch) -> Result<()> {
        self.writer
            .write(batch)
            .map_err(|e| DataFusionError::Execution(format!("write batch: {e}")))?;
        self.rows += batch.num_rows();
        Ok(())
    }

    pub fn finish(self) -> Result<usize> {
        self.writer
            .close()
            .map_err(|e| DataFusionError::Execution(format!("close ArrowWriter: {e}")))?;
        Ok(self.rows)
    }
}

/// Read a context shard, projected to the requested columns that exist in the
/// schema. Mirrors `scan_projected_existing_columns`: full scan (all rows),
/// missing columns silently dropped, surrounding quotes on a requested name
/// stripped so `"\"end\""` matches `end`. Returns batches in file order.
pub async fn read_context_parquet(
    path: &Path,
    requested_columns: &[&str],
) -> Result<Vec<RecordBatch>> {
    let file = tokio::fs::File::open(path).await.map_err(|e| {
        DataFusionError::Execution(format!("open context parquet '{}': {e}", path.display()))
    })?;
    let builder = ParquetRecordBatchStreamBuilder::new(file)
        .await
        .map_err(|e| DataFusionError::Execution(format!("read context parquet metadata: {e}")))?;
    let arrow_schema = builder.schema().clone();
    let root_indices = existing_projection_indices(&arrow_schema, requested_columns);
    let mask = ProjectionMask::roots(builder.parquet_schema(), root_indices);
    let mut stream = builder
        .with_projection(mask)
        .build()
        .map_err(|e| DataFusionError::Execution(format!("build context scan: {e}")))?;
    let mut batches = Vec::new();
    while let Some(batch) = stream
        .try_next()
        .await
        .map_err(|e| DataFusionError::Execution(format!("read context batch: {e}")))?
    {
        batches.push(batch);
    }
    Ok(batches)
}

/// Top-level field indices for the requested columns that exist (deduped,
/// quote-stripped), preserving request order.
fn existing_projection_indices(schema: &SchemaRef, requested: &[&str]) -> Vec<usize> {
    let mut seen: HashSet<usize> = HashSet::new();
    let mut out = Vec::new();
    for name in requested {
        let n = name.trim_matches('"');
        if let Ok(idx) = schema.index_of(n) {
            if seen.insert(idx) {
                out.push(idx);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Array, Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn ctx_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("peptide_seq", DataType::Utf8, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["t1", "t2", "t3"])),
                Arc::new(StringArray::from(vec!["chr1", "chr1", "chr1"])),
                Arc::new(Int64Array::from(vec![100, 200, 300])),
                Arc::new(StringArray::from(vec![Some("MAAA"), None, Some("MBBB")])),
            ],
        )
        .unwrap()
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn context_round_trips_projected_scan() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = ContextParquetShardWriter::create(tmp.path(), ctx_batch().schema()).unwrap();
        w.write(&ctx_batch()).unwrap();
        assert_eq!(w.finish().unwrap(), 3);

        // Request incl. a quoted name, a missing column, and a real subset.
        let batches = read_context_parquet(
            tmp.path(),
            &[
                "transcript_id",
                "\"start\"",
                "peptide_seq",
                "does_not_exist",
            ],
        )
        .await
        .unwrap();
        assert_eq!(batches.iter().map(RecordBatch::num_rows).sum::<usize>(), 3);
        let b = &batches[0];
        // Missing column dropped; the three present columns are readable by name.
        assert!(b.column_by_name("does_not_exist").is_none());
        let tid = b
            .column_by_name("transcript_id")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let start = b
            .column_by_name("start")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        let pep = b
            .column_by_name("peptide_seq")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!((tid.value(0), start.value(2)), ("t1", 300));
        assert!(pep.is_null(1) && pep.value(0) == "MAAA");
    }
}
