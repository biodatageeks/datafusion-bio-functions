//! Compact binary position index for probing cache windows.
//!
//! Replaces full Arrow IPC deserialization for the probing phase.
//! Only `start`, `end`, and `allele_string` are stored, in a flat
//! binary format that avoids Arrow/IPC overhead entirely.
//!
//! Binary layout:
//! ```text
//! [4B num_entries (LE)]
//! [entries: num_entries * 22 bytes each]
//!   [8B start_i64_le][8B end_i64_le][4B allele_offset_u32_le][2B allele_len_u16_le]
//! [allele_pool: remaining bytes, UTF-8]
//! ```

use std::collections::HashMap;

use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::common::{DataFusionError, Result};

use crate::allele_index::AlleleMatcher;

const ENTRY_SIZE: usize = 22; // 8 + 8 + 4 + 2

/// A single entry in the position index.
#[derive(Debug, Clone, Copy)]
struct PosEntry {
    start: i64,
    end: i64,
    allele_offset: u32,
    allele_len: u16,
}

/// Compact binary position index for a cache window.
///
/// Supports fast position-based probing and allele matching without
/// deserializing the full Arrow IPC batch.
pub struct PositionIndex {
    entries: Vec<PosEntry>,
    allele_pool: String,
    position_map: HashMap<(i64, i64), Vec<usize>>,
}

impl PositionIndex {
    /// Build from a RecordBatch (extracts start, end, allele_string columns).
    pub fn from_batch(batch: &RecordBatch) -> Result<Self> {
        let schema = batch.schema();
        let start_col = schema.index_of("start").map_err(|e| {
            DataFusionError::Execution(format!("Cache batch missing 'start' column: {e}"))
        })?;
        let end_col = schema.index_of("end").map_err(|e| {
            DataFusionError::Execution(format!("Cache batch missing 'end' column: {e}"))
        })?;
        let allele_col = schema.index_of("allele_string").map_err(|e| {
            DataFusionError::Execution(format!("Cache batch missing 'allele_string' column: {e}"))
        })?;

        let starts = batch
            .column(start_col)
            .as_primitive_opt::<datafusion::arrow::datatypes::Int64Type>()
            .ok_or_else(|| DataFusionError::Execution("'start' column is not Int64".into()))?;
        let ends = batch
            .column(end_col)
            .as_primitive_opt::<datafusion::arrow::datatypes::Int64Type>()
            .ok_or_else(|| DataFusionError::Execution("'end' column is not Int64".into()))?;

        let allele_array = batch.column(allele_col);

        let mut entries = Vec::with_capacity(batch.num_rows());
        let mut allele_pool = String::new();

        for i in 0..batch.num_rows() {
            let start = if starts.is_null(i) {
                0
            } else {
                starts.value(i)
            };
            let end = if ends.is_null(i) { 0 } else { ends.value(i) };

            let allele_str = get_string_value(allele_array, i);
            let allele_offset = allele_pool.len() as u32;
            let allele_len = allele_str.len() as u16;
            allele_pool.push_str(allele_str);

            entries.push(PosEntry {
                start,
                end,
                allele_offset,
                allele_len,
            });
        }

        let position_map = build_position_map(&entries);
        Ok(Self {
            entries,
            allele_pool,
            position_map,
        })
    }

    /// Serialize to compact binary format for fjall storage.
    pub fn to_bytes(&self) -> Vec<u8> {
        let num = self.entries.len() as u32;
        let pool_bytes = self.allele_pool.as_bytes();
        let total = 4 + self.entries.len() * ENTRY_SIZE + pool_bytes.len();
        let mut buf = Vec::with_capacity(total);

        buf.extend_from_slice(&num.to_le_bytes());
        for e in &self.entries {
            buf.extend_from_slice(&e.start.to_le_bytes());
            buf.extend_from_slice(&e.end.to_le_bytes());
            buf.extend_from_slice(&e.allele_offset.to_le_bytes());
            buf.extend_from_slice(&e.allele_len.to_le_bytes());
        }
        buf.extend_from_slice(pool_bytes);
        buf
    }

    /// Deserialize from compact binary format.
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        if data.len() < 4 {
            return Err(DataFusionError::Execution(
                "PositionIndex data too short".into(),
            ));
        }

        let num = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
        let entries_end = 4 + num * ENTRY_SIZE;
        if data.len() < entries_end {
            return Err(DataFusionError::Execution(format!(
                "PositionIndex data truncated: need {} bytes for {} entries, got {}",
                entries_end,
                num,
                data.len()
            )));
        }

        let mut entries = Vec::with_capacity(num);
        let mut offset = 4;
        for _ in 0..num {
            let start = i64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
            let end = i64::from_le_bytes(data[offset + 8..offset + 16].try_into().unwrap());
            let allele_offset =
                u32::from_le_bytes(data[offset + 16..offset + 20].try_into().unwrap());
            let allele_len = u16::from_le_bytes(data[offset + 20..offset + 22].try_into().unwrap());
            entries.push(PosEntry {
                start,
                end,
                allele_offset,
                allele_len,
            });
            offset += ENTRY_SIZE;
        }

        let allele_pool = std::str::from_utf8(&data[entries_end..])
            .map_err(|e| DataFusionError::Execution(format!("Invalid UTF-8 in allele pool: {e}")))?
            .to_string();

        let position_map = build_position_map(&entries);
        Ok(Self {
            entries,
            allele_pool,
            position_map,
        })
    }

    /// Find matching row indices by position + allele.
    pub fn find_matches(
        &self,
        start: i64,
        end: i64,
        vcf_ref: &str,
        vcf_alt: &str,
        matcher: AlleleMatcher,
    ) -> Vec<usize> {
        let mut out = Vec::new();
        self.append_matches(start, end, vcf_ref, vcf_alt, matcher, &mut out);
        out
    }

    /// Append matching row indices by position + allele to an existing buffer.
    pub fn append_matches(
        &self,
        start: i64,
        end: i64,
        vcf_ref: &str,
        vcf_alt: &str,
        matcher: AlleleMatcher,
        out: &mut Vec<usize>,
    ) {
        let Some(row_indices) = self.position_map.get(&(start, end)) else {
            return;
        };
        for &i in row_indices {
            let allele = self.allele_string(i);
            if matcher(vcf_ref, vcf_alt, allele) {
                out.push(i);
            }
        }
    }

    /// Return true if at least one row matches by position + allele.
    pub fn has_match(
        &self,
        start: i64,
        end: i64,
        vcf_ref: &str,
        vcf_alt: &str,
        matcher: AlleleMatcher,
    ) -> bool {
        let Some(row_indices) = self.position_map.get(&(start, end)) else {
            return false;
        };
        row_indices.iter().copied().any(|i| {
            let allele = self.allele_string(i);
            matcher(vcf_ref, vcf_alt, allele)
        })
    }

    /// Find all co-located row indices (position match only, no allele check).
    pub fn find_colocated(&self, start: i64, end: i64) -> Vec<usize> {
        self.position_map
            .get(&(start, end))
            .cloned()
            .unwrap_or_default()
    }

    /// Number of rows in the index.
    pub fn num_rows(&self) -> usize {
        self.entries.len()
    }

    /// Get the allele string for a row.
    fn allele_string(&self, row: usize) -> &str {
        let e = &self.entries[row];
        let start = e.allele_offset as usize;
        let end = start + e.allele_len as usize;
        &self.allele_pool[start..end]
    }
}

fn build_position_map(entries: &[PosEntry]) -> HashMap<(i64, i64), Vec<usize>> {
    let mut map: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, e) in entries.iter().enumerate() {
        map.entry((e.start, e.end)).or_default().push(i);
    }
    map
}

fn get_string_value(col: &dyn Array, i: usize) -> &str {
    if col.is_null(i) {
        return "";
    }
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        arr.value(i)
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        arr.value(i)
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::LargeStringArray>()
    {
        arr.value(i)
    } else {
        ""
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn test_matcher(vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
        let mut parts = allele_string.split('/');
        let Some(cache_ref) = parts.next() else {
            return false;
        };
        if cache_ref != vcf_ref {
            return false;
        }
        parts.any(|a| a == vcf_alt)
    }

    fn make_test_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100, 100, 200])),
                Arc::new(Int64Array::from(vec![100, 100, 100, 200])),
                Arc::new(StringArray::from(vec!["rs1", "rs2", "rs3", "rs4"])),
                Arc::new(StringArray::from(vec!["A/G", "A/T", "C/T", "G/A"])),
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_from_batch_roundtrip() {
        let batch = make_test_batch();
        let index = PositionIndex::from_batch(&batch).unwrap();
        assert_eq!(index.num_rows(), 4);

        let bytes = index.to_bytes();
        let restored = PositionIndex::from_bytes(&bytes).unwrap();
        assert_eq!(restored.num_rows(), 4);

        // Verify allele strings survived roundtrip
        assert_eq!(restored.allele_string(0), "A/G");
        assert_eq!(restored.allele_string(1), "A/T");
        assert_eq!(restored.allele_string(2), "C/T");
        assert_eq!(restored.allele_string(3), "G/A");
    }

    #[test]
    fn test_find_matches_after_roundtrip() {
        let batch = make_test_batch();
        let index = PositionIndex::from_batch(&batch).unwrap();
        let bytes = index.to_bytes();
        let index = PositionIndex::from_bytes(&bytes).unwrap();

        let matches = index.find_matches(100, 100, "A", "G", test_matcher);
        assert_eq!(matches, vec![0]);

        let matches = index.find_matches(100, 100, "A", "T", test_matcher);
        assert_eq!(matches, vec![1]);

        let matches = index.find_matches(200, 200, "G", "A", test_matcher);
        assert_eq!(matches, vec![3]);

        let matches = index.find_matches(300, 300, "A", "G", test_matcher);
        assert!(matches.is_empty());
    }

    #[test]
    fn test_find_colocated() {
        let batch = make_test_batch();
        let index = PositionIndex::from_batch(&batch).unwrap();

        assert_eq!(index.find_colocated(100, 100).len(), 3);
        assert_eq!(index.find_colocated(200, 200).len(), 1);
        assert!(index.find_colocated(300, 300).is_empty());
    }

    #[test]
    fn test_empty_batch() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::new_empty(schema);
        let index = PositionIndex::from_batch(&batch).unwrap();
        assert_eq!(index.num_rows(), 0);

        let bytes = index.to_bytes();
        let restored = PositionIndex::from_bytes(&bytes).unwrap();
        assert_eq!(restored.num_rows(), 0);
    }
}
