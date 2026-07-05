//! Tiered plugin shard writer — generalizes the variation writer to an
//! arbitrary value schema. Reuses `parquet_cache::write::point_lookup_writer_properties`
//! so plugin shards inherit the variation shard's no-dictionary, small-page,
//! page-indexed, `(tier, start)`-sorted point-lookup layout verbatim.

use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ArrowWriter;

use crate::parquet_cache::write::point_lookup_writer_properties;
use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};

fn arrow_type(ty: ValueType) -> DataType {
    match ty {
        ValueType::Utf8 => DataType::Utf8,
        ValueType::Float32 => DataType::Float32,
        ValueType::Int32 => DataType::Int32,
    }
}

/// Physical output schema for a plugin shard: shared key columns, the optional
/// per-transcript match columns (Utf8 discriminators, §3.4), the plugin's value
/// columns (nullable), then the derived `tier`.
pub fn plugin_output_schema(matches: &[MatchColumn], values: &[ValueColumn]) -> SchemaRef {
    let mut fields = vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new("allele_string", DataType::Utf8, false),
    ];
    for m in matches {
        fields.push(Field::new(&m.column, DataType::Utf8, true));
    }
    for v in values {
        fields.push(Field::new(&v.column, arrow_type(v.ty), true));
    }
    fields.push(Field::new("tier", DataType::Int8, false));
    Arc::new(Schema::new(fields))
}

/// Streaming Parquet writer for a plugin shard.
pub struct PluginShardWriter {
    writer: ArrowWriter<std::fs::File>,
    rows: usize,
}

impl PluginShardWriter {
    pub fn create(path: &Path, schema: SchemaRef) -> Result<Self> {
        let props = point_lookup_writer_properties(&schema, &["tier", "start"]);
        let file = std::fs::File::create(path).map_err(|e| {
            DataFusionError::Execution(format!("create plugin shard '{}': {e}", path.display()))
        })?;
        let writer = ArrowWriter::try_new(file, schema, Some(props))
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

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    fn f32_value_col() -> Vec<ValueColumn> {
        vec![ValueColumn {
            column: "am_pathogenicity".into(),
            csq_field: "am_pathogenicity".into(),
            ty: ValueType::Float32,
        }]
    }

    #[test]
    fn output_schema_has_key_values_tier_in_order() {
        let s = plugin_output_schema(&[], &f32_value_col());
        let names: Vec<_> = s.fields().iter().map(|f| f.name().clone()).collect();
        assert_eq!(
            names,
            vec![
                "chrom",
                "start",
                "end",
                "allele_string",
                "am_pathogenicity",
                "tier"
            ]
        );
        assert_eq!(
            s.field_with_name("tier").unwrap().data_type(),
            &DataType::Int8
        );
        assert_eq!(
            s.field_with_name("start").unwrap().data_type(),
            &DataType::UInt32
        );
    }

    #[test]
    fn match_columns_sit_between_allele_and_values() {
        let matches = vec![MatchColumn {
            column: "protein_variant".into(),
            engine_attr: "amino_acid_change".into(),
        }];
        let s = plugin_output_schema(&matches, &f32_value_col());
        let names: Vec<_> = s.fields().iter().map(|f| f.name().clone()).collect();
        assert_eq!(
            names,
            vec![
                "chrom",
                "start",
                "end",
                "allele_string",
                "protein_variant",
                "am_pathogenicity",
                "tier"
            ]
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn writes_readable_shard() {
        let schema = plugin_output_schema(&[], &f32_value_col());
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(UInt32Array::from(vec![100u32])),
                Arc::new(UInt32Array::from(vec![100u32])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(Float32Array::from(vec![0.9f32])),
                Arc::new(Int8Array::from(vec![0i8])),
            ],
        )
        .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.parquet");
        let mut w = PluginShardWriter::create(&path, schema).unwrap();
        w.write(&batch).unwrap();
        assert_eq!(w.finish().unwrap(), 1);
        assert!(path.exists());
    }
}
