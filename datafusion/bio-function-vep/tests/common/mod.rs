use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

const CACHE_SOURCE_METADATA_KEY: &str = "bio.vep.cache_source_type";

pub fn cache_with_source_metadata(cache_dir: &Path, source: &str) -> tempfile::TempDir {
    let tempdir = tempfile::TempDir::new().expect("create metadata cache tempdir");
    copy_cache_with_metadata(cache_dir, tempdir.path(), source)
        .expect("copy golden cache with source metadata");
    tempdir
}

fn copy_cache_with_metadata(
    source_dir: &Path,
    target_dir: &Path,
    source: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(target_dir)?;
    for entry in std::fs::read_dir(source_dir)? {
        let entry = entry?;
        let source_path = entry.path();
        let target_path = target_dir.join(entry.file_name());
        if source_path.is_dir() {
            copy_cache_with_metadata(&source_path, &target_path, source)?;
        } else if source_path
            .extension()
            .is_some_and(|extension| extension == "parquet")
        {
            rewrite_parquet_with_metadata(&source_path, &target_path, source)?;
        } else {
            std::fs::copy(&source_path, &target_path)?;
        }
    }
    Ok(())
}

fn rewrite_parquet_with_metadata(
    source_path: &Path,
    target_path: &Path,
    source: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let source_file = File::open(source_path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(source_file)?;
    let mut schema = builder.schema().as_ref().clone();
    let mut metadata = schema.metadata().clone();
    metadata.insert(CACHE_SOURCE_METADATA_KEY.to_string(), source.to_string());
    schema = schema.with_metadata(metadata);
    let schema = Arc::new(schema);

    if let Some(parent) = target_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let target_file = File::create(target_path)?;
    let mut writer = ArrowWriter::try_new(target_file, schema.clone(), None)?;
    for batch in builder.build()? {
        let batch = batch?;
        let batch = RecordBatch::try_new(schema.clone(), batch.columns().to_vec())?;
        writer.write(&batch)?;
    }
    writer.close()?;
    Ok(())
}
