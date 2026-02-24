//! Per-window in-memory position index for fast allele matching.
//!
//! Given a cache window (RecordBatch or PositionIndex), builds a position-based
//! index. Allele matching functions are injected as function pointers to avoid
//! a cyclic dependency on bio-function-vep.

use datafusion::arrow::array::RecordBatch;
use datafusion::common::Result;

use crate::position_index::PositionIndex;

/// Function signature for allele matching: `(vcf_ref, vcf_alt, allele_string) -> bool`.
pub type AlleleMatcher = fn(&str, &str, &str) -> bool;

/// Index over a cache window for fast position-based lookup.
///
/// Maps `(start, end)` pairs to row indices, enabling O(1) lookup by
/// position before applying allele matching.
///
/// Supports two construction paths:
/// - `from_position_index`: pre-built position index (fast, no Arrow IPC decode)
/// - `from_batch`: utility path that builds a position index from a RecordBatch
pub struct WindowAlleleIndex {
    pos_index: PositionIndex,
}

impl WindowAlleleIndex {
    /// Construct from a pre-built PositionIndex (v1 format — fast, no IPC decode).
    pub fn from_position_index(index: PositionIndex) -> Self {
        Self { pos_index: index }
    }

    /// Build an index from a cache window batch.
    ///
    /// The batch must have `start` (Int64), `end` (Int64), and `allele_string` (Utf8) columns.
    pub fn from_batch(batch: RecordBatch) -> Result<Self> {
        let pos_index = PositionIndex::from_batch(&batch)?;
        Ok(Self { pos_index })
    }

    /// Find rows matching with the given matcher function at a given position.
    pub fn find_matches(
        &self,
        start: i64,
        end: i64,
        vcf_ref: &str,
        vcf_alt: &str,
        matcher: AlleleMatcher,
    ) -> Vec<usize> {
        self.pos_index
            .find_matches(start, end, vcf_ref, vcf_alt, matcher)
    }

    /// Find all co-located rows (position match only, no allele check).
    pub fn find_colocated(&self, start: i64, end: i64) -> Vec<usize> {
        self.pos_index.find_colocated(start, end)
    }

    /// Get the number of rows in the index.
    pub fn num_rows(&self) -> usize {
        self.pos_index.num_rows()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    /// Simple exact matcher for testing (no vep dependency needed).
    fn test_exact_matcher(vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
        let mut parts = allele_string.split('/');
        let Some(cache_ref) = parts.next() else {
            return false;
        };
        if cache_ref != vcf_ref {
            return false;
        }
        parts.any(|a| a == vcf_alt)
    }

    fn make_cache_batch() -> RecordBatch {
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
    fn test_exact_match() {
        let batch = make_cache_batch();
        let index = WindowAlleleIndex::from_batch(batch).unwrap();

        // A->G should match only "A/G" (row 0)
        let matches = index.find_matches(100, 100, "A", "G", test_exact_matcher);
        assert_eq!(matches, vec![0]);

        // A->T should match only "A/T" (row 1)
        let matches = index.find_matches(100, 100, "A", "T", test_exact_matcher);
        assert_eq!(matches, vec![1]);
    }

    #[test]
    fn test_colocated() {
        let batch = make_cache_batch();
        let index = WindowAlleleIndex::from_batch(batch).unwrap();

        assert_eq!(index.find_colocated(100, 100).len(), 3);
        assert_eq!(index.find_colocated(200, 200).len(), 1);
        assert!(index.find_colocated(300, 300).is_empty());
    }

    #[test]
    fn test_from_position_index() {
        let batch = make_cache_batch();
        let pos_index = PositionIndex::from_batch(&batch).unwrap();
        let bytes = pos_index.to_bytes();
        let restored = PositionIndex::from_bytes(&bytes).unwrap();
        let index = WindowAlleleIndex::from_position_index(restored);

        let matches = index.find_matches(100, 100, "A", "G", test_exact_matcher);
        assert_eq!(matches, vec![0]);

        assert_eq!(index.find_colocated(100, 100).len(), 3);
        assert_eq!(index.num_rows(), 4);
    }
}
