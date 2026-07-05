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

/// Lookup key: `(start, allele_string, <match discriminator values…>)`. The
/// match values (§3.4) are `None` when the source row has none; a runtime probe
/// only hits when its supplied discriminator equals the stored one.
type LookupKey = (u32, String, Vec<Option<String>>);

/// In-memory point-lookup over one plugin's per-chrom shard.
pub struct PluginLookup {
    value_columns: Vec<String>,
    match_columns: Vec<String>,
    rows: HashMap<LookupKey, Vec<PluginScalar>>,
}

impl PluginLookup {
    /// Load a plugin shard into memory, projecting `(start, allele_string)`, the
    /// match discriminator columns, then the requested value columns.
    pub async fn open(
        shard: &Path,
        match_columns: Vec<String>,
        value_columns: Vec<String>,
    ) -> Result<Self> {
        let ctx = SessionContext::new();
        let path = shard.to_string_lossy();
        ctx.register_parquet("plugin_shard", path.as_ref(), ParquetReadOptions::default())
            .await?;
        let projection = ["start".to_string(), "allele_string".to_string()]
            .into_iter()
            .chain(match_columns.iter().cloned())
            .chain(value_columns.iter().cloned())
            .collect::<Vec<_>>()
            .join(", ");
        let batches = ctx
            .sql(&format!("SELECT {projection} FROM plugin_shard"))
            .await?
            .collect()
            .await?;

        let m = match_columns.len();
        let mut rows: HashMap<LookupKey, Vec<PluginScalar>> = HashMap::new();
        for batch in &batches {
            let start = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("start column not UInt32".into()))?;
            let allele = batch.column(1);
            for r in 0..batch.num_rows() {
                let mut match_vals = Vec::with_capacity(m);
                for c in 0..m {
                    match_vals.push(string_value(batch.column(2 + c).as_ref(), r)?);
                }
                let mut values = Vec::with_capacity(value_columns.len());
                for c in 0..value_columns.len() {
                    values.push(decode_scalar(batch.column(2 + m + c), r)?);
                }
                let allele_str = string_value(allele.as_ref(), r)?.ok_or_else(|| {
                    DataFusionError::Execution("null allele_string in plugin shard".into())
                })?;
                rows.insert((start.value(r), allele_str, match_vals), values);
            }
        }
        Ok(Self {
            value_columns,
            match_columns,
            rows,
        })
    }

    /// The value columns this lookup projects, in order.
    pub fn value_columns(&self) -> &[String] {
        &self.value_columns
    }

    /// The match discriminator columns this lookup keys on, in order.
    pub fn match_columns(&self) -> &[String] {
        &self.match_columns
    }

    /// Probe for `(start, allele_string, match_values)`; `None` on any miss.
    /// `match_values` must align with [`Self::match_columns`]. A `None` element
    /// (e.g. a non-missense transcript with no amino-acid change) will not match a
    /// stored `Some` discriminator, which is exactly the per-transcript gate.
    pub fn probe(
        &self,
        start: u32,
        allele_string: &str,
        match_values: &[Option<String>],
    ) -> Result<Option<PluginValueRow>> {
        let key = (start, allele_string.to_string(), match_values.to_vec());
        Ok(self.rows.get(&key).map(|values| PluginValueRow {
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
    use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};
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
        let schema = plugin_output_schema(&[], &vals);
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

        let lookup = PluginLookup::open(&path, vec![], vec!["am_pathogenicity".into()])
            .await
            .unwrap();
        let hit = lookup.probe(100, "A/G", &[]).unwrap().unwrap();
        match &hit.values[0] {
            PluginScalar::F32(v) => assert!((*v - 0.9).abs() < 1e-6),
            other => panic!("{other:?}"),
        }
        assert!(lookup.probe(100, "C/T", &[]).unwrap().is_none()); // wrong allele
        assert!(lookup.probe(999, "A/G", &[]).unwrap().is_none()); // no position
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn probe_with_match_discriminator_gates() {
        // Two rows at the same (start, allele) with distinct protein_variant —
        // the multi-isoform case. Only the matching discriminator hits.
        let matches = vec![MatchColumn {
            column: "protein_variant".into(),
            engine_attr: "amino_acid_change".into(),
        }];
        let vals = vec![ValueColumn {
            column: "am_pathogenicity".into(),
            csq_field: "am_pathogenicity".into(),
            ty: ValueType::Float32,
        }];
        let schema = plugin_output_schema(&matches, &vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 100])),
                Arc::new(UInt32Array::from(vec![100u32, 100])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["R12G", "R78G"])),
                Arc::new(Float32Array::from(vec![0.0392f32, 0.0427])),
                Arc::new(Int8Array::from(vec![1i8, 1])),
            ],
        )
        .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.parquet");
        let mut w = PluginShardWriter::create(&path, schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

        let lookup = PluginLookup::open(
            &path,
            vec!["protein_variant".into()],
            vec!["am_pathogenicity".into()],
        )
        .await
        .unwrap();

        let hit = lookup
            .probe(100, "A/G", &[Some("R78G".into())])
            .unwrap()
            .unwrap();
        match &hit.values[0] {
            PluginScalar::F32(v) => assert!((*v - 0.0427).abs() < 1e-6),
            other => panic!("{other:?}"),
        }
        // matching aa-change of the other isoform hits its own score
        let hit2 = lookup
            .probe(100, "A/G", &[Some("R12G".into())])
            .unwrap()
            .unwrap();
        match &hit2.values[0] {
            PluginScalar::F32(v) => assert!((*v - 0.0392).abs() < 1e-6),
            other => panic!("{other:?}"),
        }
        // non-missense line (no aa-change) → miss → empty (the gate)
        assert!(lookup.probe(100, "A/G", &[None]).unwrap().is_none());
        // unknown aa-change → miss
        assert!(
            lookup
                .probe(100, "A/G", &[Some("X99Y".into())])
                .unwrap()
                .is_none()
        );
    }
}
