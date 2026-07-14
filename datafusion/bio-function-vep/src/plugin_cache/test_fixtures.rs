//! Test-only fixture builders shared by the plugin-cache build tests: a gzipped
//! raw source and a minimal variation shard to join against.

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;

/// Write `body` to `path` as gzip (the compression the CSV/TSV provider expects).
pub fn write_gz(path: &Path, body: &str) {
    let f = std::fs::File::create(path).unwrap();
    let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
    enc.write_all(body.as_bytes()).unwrap();
    enc.finish().unwrap();
}

/// Write `body` to `path` as **BGZF** — the block-gzip every real `.vcf.gz` uses.
///
/// Not interchangeable with [`write_gz`] here: bio-formats' `get_compression_type`
/// sniffs the first 18 bytes for BGZF's `BC` extra subfield and routes BGZF and plain
/// GZIP down *different* header readers. A flate2 fixture would silently exercise the
/// GZIP path and leave the BGZF one — the one production actually takes — untested.
pub fn write_bgzf(path: &Path, body: &str) {
    let f = std::fs::File::create(path).unwrap();
    let mut w = noodles_bgzf::io::Writer::new(f);
    w.write_all(body.as_bytes()).unwrap();
    w.finish().unwrap();
}

/// Write a minimal `variation/<chrom>.parquet` shard: the columns the plugin build
/// joins against (`chrom`, `start`, `allele_string`) plus the inherited `tier`.
pub fn write_variation(path: &Path, rows: &[(&str, u32, &str, i8)]) {
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("allele_string", DataType::Utf8, false),
        Field::new("tier", DataType::Int8, false),
    ]));
    let batch = RecordBatch::try_new(
        Arc::clone(&schema),
        vec![
            Arc::new(StringArray::from(
                rows.iter().map(|r| r.0).collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter().map(|r| r.1).collect::<Vec<_>>(),
            )),
            Arc::new(StringArray::from(
                rows.iter().map(|r| r.2).collect::<Vec<_>>(),
            )),
            Arc::new(Int8Array::from(
                rows.iter().map(|r| r.3).collect::<Vec<_>>(),
            )),
        ],
    )
    .unwrap();
    let file = std::fs::File::create(path).unwrap();
    let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
}
