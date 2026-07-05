//! Runtime per-chrom plugin probe.
//!
//! The annotation engine must probe plugins **synchronously on the annotating
//! thread** (no async/Rayon pool — the aux-pool oversubscription constraint).
//! The variation `PageDir` take path is async, so it cannot be called
//! synchronously per variant without blocking. For the prototype `open()`
//! (async) loads the per-chrom shard into an in-memory map keyed by
//! `(start, allele_string)`, and `probe()` is a synchronous map lookup. A
//! sync `PageDir` reader that keeps the shard on disk is a later optimization
//! if plugin-cache memory becomes a concern.

use std::collections::HashMap;
use std::path::Path;

use datafusion::arrow::array::{
    Array, Float32Array, Int32Array, LargeStringArray, StringArray, StringViewArray, UInt32Array,
};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{ParquetReadOptions, SessionContext};

/// One decoded plugin value.
#[derive(Debug, Clone, PartialEq)]
pub enum PluginScalar {
    Str(String),
    F32(f32),
    I32(i32),
    Null,
}

/// A plugin's values for one variant, positionally aligned to `value_columns`.
#[derive(Debug, Clone)]
pub struct PluginValueRow {
    pub values: Vec<PluginScalar>,
}

/// In-memory point-lookup over one plugin's per-chrom shard.
pub struct PluginLookup {
    value_columns: Vec<String>,
    rows: HashMap<(u32, String), Vec<PluginScalar>>,
}

impl PluginLookup {
    /// Load a plugin shard into memory, projecting `(start, allele_string)` plus
    /// the requested value columns.
    pub async fn open(shard: &Path, value_columns: Vec<String>) -> Result<Self> {
        let ctx = SessionContext::new();
        let path = shard.to_string_lossy();
        ctx.register_parquet("plugin_shard", path.as_ref(), ParquetReadOptions::default())
            .await?;
        let projection = std::iter::once("start".to_string())
            .chain(std::iter::once("allele_string".to_string()))
            .chain(value_columns.iter().cloned())
            .collect::<Vec<_>>()
            .join(", ");
        let batches = ctx
            .sql(&format!("SELECT {projection} FROM plugin_shard"))
            .await?
            .collect()
            .await?;

        let mut rows: HashMap<(u32, String), Vec<PluginScalar>> = HashMap::new();
        for batch in &batches {
            let start = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("start column not UInt32".into()))?;
            let allele = batch.column(1);
            for r in 0..batch.num_rows() {
                let mut values = Vec::with_capacity(value_columns.len());
                for c in 0..value_columns.len() {
                    values.push(decode_scalar(batch.column(2 + c), r)?);
                }
                let allele_str = string_value(allele.as_ref(), r)?.ok_or_else(|| {
                    DataFusionError::Execution("null allele_string in plugin shard".into())
                })?;
                rows.insert((start.value(r), allele_str), values);
            }
        }
        Ok(Self {
            value_columns,
            rows,
        })
    }

    /// The value columns this lookup projects, in order.
    pub fn value_columns(&self) -> &[String] {
        &self.value_columns
    }

    /// Probe for `(start, allele_string)`; `None` on a position/allele miss.
    pub fn probe(&self, start: u32, allele_string: &str) -> Result<Option<PluginValueRow>> {
        Ok(self
            .rows
            .get(&(start, allele_string.to_string()))
            .map(|values| PluginValueRow {
                values: values.clone(),
            }))
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
    // Strings (any encoding) fall through to the shared reader.
    if let Some(s) = string_value(col, row)? {
        return Ok(PluginScalar::Str(s));
    }
    Ok(PluginScalar::Null)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::{ValueColumn, ValueType};
    use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[tokio::test(flavor = "multi_thread")]
    async fn probe_hit_and_miss() {
        let vals = vec![ValueColumn {
            column: "am_pathogenicity".into(),
            csq_field: "am_pathogenicity".into(),
            ty: ValueType::Float32,
        }];
        let schema = plugin_output_schema(&vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Float32Array::from(vec![0.9f32, 0.1])),
                Arc::new(Int8Array::from(vec![0i8, 0])),
            ],
        )
        .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.parquet");
        let mut w = PluginShardWriter::create(&path, schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

        let lookup = PluginLookup::open(&path, vec!["am_pathogenicity".into()])
            .await
            .unwrap();
        let hit = lookup.probe(100, "A/G").unwrap().unwrap();
        match &hit.values[0] {
            PluginScalar::F32(v) => assert!((*v - 0.9).abs() < 1e-6),
            other => panic!("{other:?}"),
        }
        assert!(lookup.probe(100, "C/T").unwrap().is_none()); // wrong allele
        assert!(lookup.probe(999, "A/G").unwrap().is_none()); // no position
    }
}
