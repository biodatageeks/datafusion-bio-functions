use std::ops::Range;
use std::sync::OnceLock;

use datafusion::arrow::array::{
    Array, Int8Array, Int16Array, Int32Array, Int64Array, ListArray, RecordBatch, StringArray,
    UInt8Array, UInt16Array, UInt32Array, UInt64Array,
};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use smallvec::SmallVec;

use crate::warm_cache::key::{VariantIndex, new_variant_index};

#[derive(Debug, Clone, Copy)]
pub struct WarmColumnIndices {
    pub position_key: usize,
    pub variant_keys: Option<usize>,
    pub allele_string: Option<usize>,
    pub end: Option<usize>,
    pub failed: Option<usize>,
}

#[derive(Debug)]
pub struct WarmChunkContext {
    pub row_group_id: usize,
    pub min_position_key: u64,
    pub max_position_key: u64,
    position_rows: Vec<PositionRowRange>,
    pub variant_index: VariantIndex,
    allele_strings: Option<StringArray>,
    output_indices: OnceLock<Option<Vec<Option<usize>>>>,
    pub batch: RecordBatch,
    pub columns: WarmColumnIndices,
}

#[derive(Debug, Clone, Copy)]
struct PositionRowRange {
    position_key: u64,
    start: u32,
    end: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WarmProbeResult {
    Hit(SmallVec<[u32; 1]>),
    DefinitiveMiss,
    NotCovered,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WarmChunkProbe {
    Exact { row: u32, position_rows: Range<u32> },
    PositionCoveredNoExact { position_rows: Range<u32> },
    NotCovered,
}

impl WarmChunkContext {
    pub fn try_new(row_group_id: usize, batch: RecordBatch) -> Result<Self> {
        Self::try_new_with_variant_index(row_group_id, batch, true)
    }

    pub fn try_new_without_variant_index(row_group_id: usize, batch: RecordBatch) -> Result<Self> {
        Self::try_new_with_variant_index(row_group_id, batch, false)
    }

    fn try_new_with_variant_index(
        row_group_id: usize,
        batch: RecordBatch,
        build_variant_index: bool,
    ) -> Result<Self> {
        let columns = WarmColumnIndices::new(&batch.schema())?;
        let position_array = batch
            .column(columns.position_key)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| DataFusionError::Execution("position_key must be UInt64".into()))?;
        let variant_arrays = if build_variant_index {
            let variant_keys_idx = columns.variant_keys.ok_or_else(|| {
                DataFusionError::Execution("required warm column variant_keys not found".into())
            })?;
            let variant_array = batch
                .column(variant_keys_idx)
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
            Some((variant_array, variant_values, variant_array.offsets()))
        } else {
            None
        };
        let allele_strings = columns
            .allele_string
            .map(|idx| {
                batch
                    .column(idx)
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .cloned()
                    .ok_or_else(|| {
                        DataFusionError::Execution(
                            "warm allele_string column must be Utf8".to_string(),
                        )
                    })
            })
            .transpose()?;

        let mut min_position_key = u64::MAX;
        let mut max_position_key = 0_u64;
        let mut saw_position = false;
        let mut position_rows: Vec<PositionRowRange> = Vec::new();
        let mut variant_index = if build_variant_index {
            new_variant_index(batch.num_rows())
        } else {
            new_variant_index(0)
        };

        for row in 0..batch.num_rows() {
            if position_array.is_null(row) {
                continue;
            }

            let position_key = position_array.value(row);
            saw_position = true;
            min_position_key = min_position_key.min(position_key);
            max_position_key = max_position_key.max(position_key);
            match position_rows.last_mut() {
                Some(range) if range.position_key == position_key => {
                    range.end = row as u32 + 1;
                }
                _ => position_rows.push(PositionRowRange {
                    position_key,
                    start: row as u32,
                    end: row as u32 + 1,
                }),
            }

            if let Some((variant_array, variant_values, offsets)) = variant_arrays.as_ref() {
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
        }

        if !saw_position {
            min_position_key = 0;
        }

        Ok(Self {
            row_group_id,
            min_position_key,
            max_position_key,
            position_rows,
            variant_index,
            allele_strings,
            output_indices: OnceLock::new(),
            batch,
            columns,
        })
    }

    pub fn contains_position(&self, position_key: u64) -> bool {
        position_key >= self.min_position_key
            && position_key <= self.max_position_key
            && !self.rows_for_position(position_key).is_empty()
    }

    pub fn lookup_variant(&self, variant_key: u64) -> SmallVec<[u32; 1]> {
        self.variant_index
            .get(&variant_key)
            .cloned()
            .unwrap_or_default()
    }

    pub fn lookup_variant_rows(&self, variant_key: u64) -> &[u32] {
        self.variant_index
            .get(&variant_key)
            .map(SmallVec::as_slice)
            .unwrap_or(&[])
    }

    pub fn rows_for_position(&self, position_key: u64) -> Range<u32> {
        self.position_rows
            .binary_search_by_key(&position_key, |range| range.position_key)
            .map(|idx| self.position_rows[idx].start..self.position_rows[idx].end)
            .unwrap_or(0..0)
    }

    pub fn is_boundary_position(&self, position_key: u64) -> bool {
        position_key == self.min_position_key || position_key == self.max_position_key
    }

    pub fn allele_string(&self, row: usize) -> Result<Option<&str>> {
        let Some(alleles) = self.allele_strings.as_ref() else {
            return Ok(None);
        };
        if row >= alleles.len() || alleles.is_null(row) {
            Ok(None)
        } else {
            Ok(Some(alleles.value(row)))
        }
    }

    pub fn i64_value(&self, column_idx: Option<usize>, row: usize) -> Option<i64> {
        let array = self.batch.column(column_idx?);
        if row >= array.len() || array.is_null(row) {
            return None;
        }

        if let Some(array) = array.as_any().downcast_ref::<Int8Array>() {
            Some(array.value(row) as i64)
        } else if let Some(array) = array.as_any().downcast_ref::<Int16Array>() {
            Some(array.value(row) as i64)
        } else if let Some(array) = array.as_any().downcast_ref::<Int32Array>() {
            Some(array.value(row) as i64)
        } else if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
            Some(array.value(row))
        } else if let Some(array) = array.as_any().downcast_ref::<UInt8Array>() {
            Some(array.value(row) as i64)
        } else if let Some(array) = array.as_any().downcast_ref::<UInt16Array>() {
            Some(array.value(row) as i64)
        } else if let Some(array) = array.as_any().downcast_ref::<UInt32Array>() {
            Some(array.value(row) as i64)
        } else {
            array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .map(|array| array.value(row) as i64)
        }
    }

    pub fn output_indices(
        &self,
        cache_columns: &[String],
        col_map: &[usize],
    ) -> Option<&[Option<usize>]> {
        self.output_indices
            .get_or_init(|| build_output_indices(&self.batch, cache_columns, col_map))
            .as_deref()
    }

    pub fn probe(&self, position_key: u64, variant_key: u64) -> WarmProbeResult {
        let rows = self.lookup_variant_rows(variant_key);
        if !rows.is_empty() {
            return WarmProbeResult::Hit(rows.iter().copied().collect());
        }
        if self.contains_position(position_key) {
            WarmProbeResult::DefinitiveMiss
        } else {
            WarmProbeResult::NotCovered
        }
    }

    pub fn probe_exact<F>(
        &self,
        position_key: u64,
        variant_key: u64,
        mut verify_row: F,
    ) -> Result<WarmChunkProbe>
    where
        F: FnMut(u32, &str) -> Result<bool>,
    {
        let position_rows = self.rows_for_position(position_key);
        if position_rows.is_empty() {
            return Ok(WarmChunkProbe::NotCovered);
        }

        for &row in self.lookup_variant_rows(variant_key) {
            if position_rows.contains(&row)
                && let Some(allele_string) = self.allele_string(row as usize)?
                && verify_row(row, allele_string)?
            {
                return Ok(WarmChunkProbe::Exact { row, position_rows });
            }
        }

        Ok(WarmChunkProbe::PositionCoveredNoExact { position_rows })
    }
}

fn build_output_indices(
    batch: &RecordBatch,
    cache_columns: &[String],
    col_map: &[usize],
) -> Option<Vec<Option<usize>>> {
    cache_columns
        .iter()
        .zip(col_map.iter())
        .map(|(name, entry_idx)| {
            if *entry_idx == usize::MAX {
                Some(None)
            } else {
                batch.schema().index_of(name).ok().map(Some)
            }
        })
        .collect()
}

impl WarmColumnIndices {
    fn new(schema: &SchemaRef) -> Result<Self> {
        Ok(Self {
            position_key: schema.index_of("position_key").map_err(|_| {
                DataFusionError::Execution("required warm column position_key not found".into())
            })?,
            variant_keys: schema.index_of("variant_keys").ok(),
            allele_string: schema.index_of("allele_string").ok(),
            end: schema.index_of("end").ok(),
            failed: schema.index_of("failed").ok(),
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
            Field::new("end", DataType::Int64, false),
            Field::new("failed", DataType::Int64, false),
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
                Arc::new(Int64Array::from(vec![101, 101, 205])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0, 1, 0])) as ArrayRef,
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
        assert_eq!(
            chunk.rows_for_position(position_key("1", 101).unwrap()),
            0..2
        );
        assert_eq!(chunk.lookup_variant(k2).as_slice(), &[1]);
        assert_eq!(chunk.lookup_variant_rows(k2), &[1]);
        assert_eq!(chunk.allele_string(1).unwrap(), Some("A/T"));
        assert_eq!(chunk.i64_value(chunk.columns.end, 2), Some(205));
        assert_eq!(chunk.i64_value(chunk.columns.failed, 1), Some(1));
    }

    #[test]
    fn chunk_context_can_skip_variant_index_for_position_lookup() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", 101).unwrap(),
                    position_key("1", 101).unwrap(),
                ])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G", "A/T"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![101, 101])) as ArrayRef,
            ],
        )
        .unwrap();

        let chunk = WarmChunkContext::try_new_without_variant_index(3, batch).unwrap();

        assert!(chunk.contains_position(position_key("1", 101).unwrap()));
        assert_eq!(
            chunk.rows_for_position(position_key("1", 101).unwrap()),
            0..2
        );
        assert_eq!(chunk.allele_string(1).unwrap(), Some("A/T"));
        assert!(chunk.lookup_variant(42).is_empty());
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

    #[test]
    fn chunk_probe_exact_verifies_candidate_rows() {
        let batch = test_batch();
        let chunk = WarmChunkContext::try_new(0, batch).unwrap();
        let hit_key = variant_key("1", 101, "A", "T").unwrap();

        let probe = chunk
            .probe_exact(position_key("1", 101).unwrap(), hit_key, |_, allele| {
                Ok(allele == "A/T")
            })
            .unwrap();

        assert_eq!(
            probe,
            WarmChunkProbe::Exact {
                row: 1,
                position_rows: 0..2
            }
        );

        let miss = chunk
            .probe_exact(position_key("1", 101).unwrap(), 999, |_, _| Ok(true))
            .unwrap();
        assert_eq!(
            miss,
            WarmChunkProbe::PositionCoveredNoExact {
                position_rows: 0..2
            }
        );

        let not_covered = chunk
            .probe_exact(position_key("1", 999).unwrap(), 999, |_, _| Ok(true))
            .unwrap();
        assert_eq!(not_covered, WarmChunkProbe::NotCovered);
    }
}
