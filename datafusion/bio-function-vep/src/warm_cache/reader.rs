use std::fs::File;
use std::path::Path;

use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::schema::types::SchemaDescriptor;

use crate::warm_cache::chunk::WarmChunkContext;

pub const WARM_RUNTIME_COLUMNS: &[&str] = &[
    "position_key",
    "variant_keys",
    "chrom",
    "start",
    "end",
    "allele_string",
    "variation_name",
    "failed",
    "somatic",
    "strand",
    "minor_allele",
    "minor_allele_freq",
    "clin_sig",
    "phenotype_or_disease",
    "clinical_impact",
    "clin_sig_allele",
    "AF",
    "gnomADg",
    "gnomADe",
];

pub fn load_warm_chunk_row_group(
    path: impl AsRef<Path>,
    row_group_id: usize,
    batch_size: usize,
) -> Result<WarmChunkContext> {
    let file = File::open(path.as_ref()).map_err(|e| {
        DataFusionError::Execution(format!("failed to open warm parquet chunk: {e}"))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    if row_group_id >= builder.metadata().num_row_groups() {
        return Err(DataFusionError::Execution(format!(
            "warm row group {row_group_id} out of range; file has {} row groups",
            builder.metadata().num_row_groups()
        )));
    }

    let mask = projection_for_existing_roots(
        builder.schema(),
        builder.parquet_schema(),
        WARM_RUNTIME_COLUMNS,
    );
    let reader = builder
        .with_projection(mask)
        .with_row_groups(vec![row_group_id])
        .with_batch_size(batch_size)
        .build()?;

    let batches = reader.collect::<std::result::Result<Vec<_>, _>>()?;
    let batch = match batches.as_slice() {
        [] => {
            return Err(DataFusionError::Execution(format!(
                "warm row group {row_group_id} produced no batches"
            )));
        }
        [single] => single.clone(),
        _ => concat_batches(&batches[0].schema(), batches.iter())
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?,
    };

    WarmChunkContext::try_new(row_group_id, batch)
}

pub fn projection_for_existing_roots(
    arrow_schema: &SchemaRef,
    parquet_schema: &SchemaDescriptor,
    names: &[&str],
) -> ProjectionMask {
    let root_indices = arrow_schema
        .fields()
        .iter()
        .enumerate()
        .filter_map(|(idx, field)| {
            names
                .iter()
                .any(|name| field.name().as_str() == *name)
                .then_some(idx)
        });
    ProjectionMask::roots(parquet_schema, root_indices)
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::sync::Arc;

    use datafusion::arrow::array::{
        ArrayRef, Int64Array, ListBuilder, RecordBatch, StringArray, UInt64Array, UInt64Builder,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;

    use super::*;
    use crate::warm_cache::key::{position_key, variant_key};

    fn warm_batch(starts: &[i64], alts: &[&str]) -> RecordBatch {
        let mut variant_keys = ListBuilder::new(UInt64Builder::new());
        for (start, allele) in starts.iter().zip(alts.iter()) {
            let (reference, alternate) = allele.split_once('/').unwrap();
            variant_keys
                .values()
                .append_value(variant_key("1", *start, reference, alternate).unwrap());
            variant_keys.append(true);
        }

        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new_list(
                "variant_keys",
                Arc::new(Field::new_list_field(DataType::UInt64, true)),
                false,
            ),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("unused_column", DataType::Utf8, false),
        ]));

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(
                    starts
                        .iter()
                        .map(|start| position_key("1", *start).unwrap())
                        .collect::<Vec<_>>(),
                )) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["1"; starts.len()])) as ArrayRef,
                Arc::new(Int64Array::from(starts.to_vec())) as ArrayRef,
                Arc::new(StringArray::from(alts.to_vec())) as ArrayRef,
                Arc::new(StringArray::from(vec!["ignored"; starts.len()])) as ArrayRef,
            ],
        )
        .unwrap()
    }

    #[test]
    fn reader_loads_requested_row_group_as_warm_chunk() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_warm.parquet");
        let file = File::create(&path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(2)
            .build();
        let batch1 = warm_batch(&[101, 102], &["A/G", "C/T"]);
        let batch2 = warm_batch(&[205], &["G/A"]);
        let mut writer = ArrowWriter::try_new(file, batch1.schema(), Some(props)).unwrap();
        writer.write(&batch1).unwrap();
        writer.write(&batch2).unwrap();
        writer.close().unwrap();

        let chunk = load_warm_chunk_row_group(&path, 1, 1024).unwrap();

        assert_eq!(chunk.row_group_id, 1);
        assert!(chunk.contains_position(position_key("1", 205).unwrap()));
        assert!(!chunk.contains_position(position_key("1", 101).unwrap()));
        assert!(chunk.batch.schema().index_of("unused_column").is_err());
    }
}
