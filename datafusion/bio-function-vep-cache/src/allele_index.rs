//! Per-window in-memory position index for fast allele matching.
//!
//! Given a cache window (RecordBatch), builds a position-based index.
//! Allele matching functions are injected as function pointers to avoid
//! a cyclic dependency on bio-function-vep.

use std::collections::HashMap;

use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::common::{DataFusionError, Result};

/// Function signature for allele matching: `(vcf_ref, vcf_alt, allele_string) -> bool`.
pub type AlleleMatcher = fn(&str, &str, &str) -> bool;

/// Index over a cache window's RecordBatch for fast position-based lookup.
///
/// Maps `(start, end)` pairs to row indices in the batch, enabling
/// O(1) lookup by position before applying allele matching.
pub struct WindowAlleleIndex {
    batch: RecordBatch,
    position_index: HashMap<(i64, i64), Vec<usize>>,
    allele_string_col: usize,
}

impl WindowAlleleIndex {
    /// Build an index from a cache window batch.
    ///
    /// The batch must have `start` (Int64), `end` (Int64), and `allele_string` (Utf8) columns.
    pub fn from_batch(batch: RecordBatch) -> Result<Self> {
        let schema = batch.schema();
        let start_col = schema.index_of("start").map_err(|e| {
            DataFusionError::Execution(format!("Cache batch missing 'start' column: {e}"))
        })?;
        let end_col = schema.index_of("end").map_err(|e| {
            DataFusionError::Execution(format!("Cache batch missing 'end' column: {e}"))
        })?;
        let allele_string_col = schema.index_of("allele_string").map_err(|e| {
            DataFusionError::Execution(format!("Cache batch missing 'allele_string' column: {e}"))
        })?;

        let starts = batch
            .column(start_col)
            .as_primitive_opt::<datafusion::arrow::datatypes::Int64Type>()
            .ok_or_else(|| DataFusionError::Execution("'start' column is not Int64".to_string()))?;
        let ends = batch
            .column(end_col)
            .as_primitive_opt::<datafusion::arrow::datatypes::Int64Type>()
            .ok_or_else(|| DataFusionError::Execution("'end' column is not Int64".to_string()))?;

        let mut position_index: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for i in 0..batch.num_rows() {
            if !starts.is_null(i) && !ends.is_null(i) {
                let start = starts.value(i);
                let end = ends.value(i);
                position_index.entry((start, end)).or_default().push(i);
            }
        }

        Ok(Self {
            batch,
            position_index,
            allele_string_col,
        })
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
        let Some(row_indices) = self.position_index.get(&(start, end)) else {
            return Vec::new();
        };

        let allele_col = self.batch.column(self.allele_string_col);
        if let Some(arr) = allele_col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringArray>()
        {
            row_indices
                .iter()
                .copied()
                .filter(|&i| !arr.is_null(i) && matcher(vcf_ref, vcf_alt, arr.value(i)))
                .collect()
        } else if let Some(arr) = allele_col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringViewArray>()
        {
            row_indices
                .iter()
                .copied()
                .filter(|&i| !arr.is_null(i) && matcher(vcf_ref, vcf_alt, arr.value(i)))
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Find all co-located rows (position match only, no allele check).
    pub fn find_colocated(&self, start: i64, end: i64) -> Vec<usize> {
        self.position_index
            .get(&(start, end))
            .cloned()
            .unwrap_or_default()
    }

    /// Get the underlying batch.
    pub fn batch(&self) -> &RecordBatch {
        &self.batch
    }

    /// Get the number of rows in the index.
    pub fn num_rows(&self) -> usize {
        self.batch.num_rows()
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
}
