use datafusion::arrow::array::{Array, ListArray, RecordBatch, UInt64Array};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use smallvec::SmallVec;

use crate::warm_cache::key::{VariantIndex, new_variant_index};

#[derive(Debug, Clone, Copy)]
pub struct WarmColumnIndices {
    pub position_key: usize,
    pub variant_keys: usize,
}

#[derive(Debug)]
pub struct WarmChunkContext {
    pub row_group_id: usize,
    pub min_position_key: u64,
    pub max_position_key: u64,
    pub position_keys: Vec<u64>,
    pub variant_index: VariantIndex,
    pub batch: RecordBatch,
    pub columns: WarmColumnIndices,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WarmProbeResult {
    Hit(SmallVec<[u32; 1]>),
    DefinitiveMiss,
    NotCovered,
}

impl WarmChunkContext {
    pub fn try_new(row_group_id: usize, batch: RecordBatch) -> Result<Self> {
        let columns = WarmColumnIndices::new(&batch.schema())?;
        let position_array = batch
            .column(columns.position_key)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| DataFusionError::Execution("position_key must be UInt64".into()))?;
        let variant_array = batch
            .column(columns.variant_keys)
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| {
                DataFusionError::Execution("variant_keys must be List<UInt64>".into())
            })?;
        let variant_values = variant_array
            .values()
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| {
                DataFusionError::Execution("variant_keys values must be UInt64".into())
            })?;

        let mut position_keys = Vec::with_capacity(batch.num_rows());
        let mut variant_index = new_variant_index(batch.num_rows());
        let offsets = variant_array.offsets();

        for row in 0..batch.num_rows() {
            if position_array.is_null(row) {
                continue;
            }

            position_keys.push(position_array.value(row));

            if variant_array.is_null(row) {
                continue;
            }

            let start = offsets[row] as usize;
            let end = offsets[row + 1] as usize;
            for value_idx in start..end {
                if variant_values.is_null(value_idx) {
                    continue;
                }
                variant_index
                    .entry(variant_values.value(value_idx))
                    .or_default()
                    .push(row as u32);
            }
        }

        position_keys.sort_unstable();
        position_keys.dedup();
        let min_position_key = position_keys.first().copied().unwrap_or(0);
        let max_position_key = position_keys.last().copied().unwrap_or(0);

        Ok(Self {
            row_group_id,
            min_position_key,
            max_position_key,
            position_keys,
            variant_index,
            batch,
            columns,
        })
    }

    pub fn contains_position(&self, position_key: u64) -> bool {
        position_key >= self.min_position_key
            && position_key <= self.max_position_key
            && self.position_keys.binary_search(&position_key).is_ok()
    }

    pub fn lookup_variant(&self, variant_key: u64) -> SmallVec<[u32; 1]> {
        self.variant_index
            .get(&variant_key)
            .cloned()
            .unwrap_or_default()
    }

    pub fn probe(&self, position_key: u64, variant_key: u64) -> WarmProbeResult {
        let rows = self.lookup_variant(variant_key);
        if !rows.is_empty() {
            return WarmProbeResult::Hit(rows);
        }
        if self.contains_position(position_key) {
            WarmProbeResult::DefinitiveMiss
        } else {
            WarmProbeResult::NotCovered
        }
    }
}

impl WarmColumnIndices {
    fn new(schema: &SchemaRef) -> Result<Self> {
        Ok(Self {
            position_key: schema.index_of("position_key").map_err(|_| {
                DataFusionError::Execution("required warm column position_key not found".into())
            })?,
            variant_keys: schema.index_of("variant_keys").map_err(|_| {
                DataFusionError::Execution("required warm column variant_keys not found".into())
            })?,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use datafusion::arrow::array::{
        ArrayRef, Int64Array, ListBuilder, RecordBatch, StringArray, UInt64Array, UInt64Builder,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    use super::*;
    use crate::warm_cache::key::{position_key, variant_key};

    fn test_batch() -> RecordBatch {
        let k1 = variant_key("1", 101, "A", "G").unwrap();
        let k2 = variant_key("1", 101, "A", "T").unwrap();
        let k3 = variant_key("1", 205, "C", "A").unwrap();

        let mut variant_keys = ListBuilder::new(UInt64Builder::new());
        variant_keys.values().append_value(k1);
        variant_keys.append(true);
        variant_keys.values().append_value(k2);
        variant_keys.append(true);
        variant_keys.values().append_value(k3);
        variant_keys.append(true);

        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new_list(
                "variant_keys",
                Arc::new(Field::new_list_field(DataType::UInt64, true)),
                false,
            ),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
        ]));

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", 101).unwrap(),
                    position_key("1", 101).unwrap(),
                    position_key("1", 205).unwrap(),
                ])) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G", "A/T", "C/A"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![101, 101, 205])) as ArrayRef,
            ],
        )
        .unwrap()
    }

    #[test]
    fn chunk_context_indexes_positions_and_variant_keys() {
        let batch = test_batch();
        let k2 = variant_key("1", 101, "A", "T").unwrap();
        let chunk = WarmChunkContext::try_new(7, batch).unwrap();

        assert_eq!(chunk.row_group_id, 7);
        assert!(chunk.contains_position(position_key("1", 101).unwrap()));
        assert!(chunk.contains_position(position_key("1", 205).unwrap()));
        assert_eq!(chunk.lookup_variant(k2).as_slice(), &[1]);
    }

    #[test]
    fn chunk_probe_distinguishes_hit_definitive_miss_and_not_covered() {
        let batch = test_batch();
        let chunk = WarmChunkContext::try_new(0, batch).unwrap();

        let hit_key = variant_key("1", 101, "A", "T").unwrap();
        assert_eq!(
            chunk.probe(position_key("1", 101).unwrap(), hit_key),
            WarmProbeResult::Hit(vec![1].into())
        );

        assert_eq!(
            chunk.probe(position_key("1", 101).unwrap(), 999),
            WarmProbeResult::DefinitiveMiss
        );

        assert_eq!(
            chunk.probe(position_key("1", 999).unwrap(), 999),
            WarmProbeResult::NotCovered
        );
    }
}
