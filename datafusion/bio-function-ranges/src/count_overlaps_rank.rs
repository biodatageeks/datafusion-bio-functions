use std::sync::Arc;

use ahash::AHashMap;
use datafusion::arrow::array::{ArrayRef, Int64Array, RecordBatch, RecordBatchOptions};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};

use crate::array_utils::get_join_col_arrays;
use crate::filter_op::FilterOp;

pub const MISSING_CONTIG_ID: u32 = u32::MAX;

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ContigRange {
    pub start_offset: u32,
    pub end_offset: u32,
    pub len: u32,
    pub _pad: u32,
}

#[derive(Debug)]
pub struct CountOverlapsRankIndex {
    chrom_to_id: AHashMap<String, u32>,
    contigs: Vec<ContigRange>,
    starts_sorted: Vec<i64>,
    ends_sorted: Vec<i64>,
    interval_count: usize,
}

#[derive(Debug)]
pub struct EncodedQueryBatch {
    pub query_contigs: Vec<u32>,
    pub query_starts: Vec<i64>,
    pub query_ends: Vec<i64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OutputColumnSource {
    Input(usize),
    Count,
}

impl CountOverlapsRankIndex {
    pub fn build_from_batches(
        batches: &[RecordBatch],
        columns: (&str, &str, &str),
    ) -> Result<Self> {
        let mut by_contig = AHashMap::<String, (Vec<i64>, Vec<i64>)>::default();
        let mut interval_count = 0usize;

        for batch in batches {
            let (contig_arr, start_arr, end_arr) =
                get_join_col_arrays(batch, (columns.0, columns.1, columns.2))?;

            for i in 0..batch.num_rows() {
                let contig = contig_arr.value(i);
                let start = i64::from(start_arr.value(i)?);
                let end = i64::from(end_arr.value(i)?);
                if let Some((starts, ends)) = by_contig.get_mut(contig) {
                    starts.push(start);
                    ends.push(end);
                } else {
                    by_contig.insert(contig.to_owned(), (vec![start], vec![end]));
                }
                interval_count += 1;
            }
        }

        let mut names = by_contig.keys().cloned().collect::<Vec<_>>();
        names.sort_unstable();

        let mut chrom_to_id = AHashMap::with_capacity(names.len());
        let mut contigs = Vec::with_capacity(names.len());
        let mut starts_sorted = Vec::with_capacity(interval_count);
        let mut ends_sorted = Vec::with_capacity(interval_count);

        for (chrom_id, name) in names.into_iter().enumerate() {
            let chrom_id = u32::try_from(chrom_id).map_err(|_| {
                DataFusionError::Execution(
                    "count_overlaps GPU index supports at most u32::MAX contigs".to_string(),
                )
            })?;
            let (mut starts, mut ends) = by_contig.remove(&name).ok_or_else(|| {
                DataFusionError::Internal(format!(
                    "missing contig while building rank index: {name}"
                ))
            })?;
            starts.sort_unstable();
            ends.sort_unstable();

            let start_offset = u32::try_from(starts_sorted.len()).map_err(|_| {
                DataFusionError::Execution(
                    "count_overlaps GPU index start offset exceeds u32::MAX".to_string(),
                )
            })?;
            let end_offset = u32::try_from(ends_sorted.len()).map_err(|_| {
                DataFusionError::Execution(
                    "count_overlaps GPU index end offset exceeds u32::MAX".to_string(),
                )
            })?;
            let len = u32::try_from(starts.len()).map_err(|_| {
                DataFusionError::Execution(
                    "count_overlaps GPU index contig length exceeds u32::MAX".to_string(),
                )
            })?;

            starts_sorted.extend_from_slice(&starts);
            ends_sorted.extend_from_slice(&ends);
            chrom_to_id.insert(name, chrom_id);
            contigs.push(ContigRange {
                start_offset,
                end_offset,
                len,
                _pad: 0,
            });
        }

        Ok(Self {
            chrom_to_id,
            contigs,
            starts_sorted,
            ends_sorted,
            interval_count,
        })
    }

    pub fn interval_count(&self) -> usize {
        self.interval_count
    }

    pub fn contigs(&self) -> &[ContigRange] {
        &self.contigs
    }

    pub fn starts_sorted(&self) -> &[i64] {
        &self.starts_sorted
    }

    pub fn ends_sorted(&self) -> &[i64] {
        &self.ends_sorted
    }

    pub fn encode_query_batch(
        &self,
        batch: &RecordBatch,
        columns: (&str, &str, &str),
        filter_op: FilterOp,
    ) -> Result<EncodedQueryBatch> {
        let (contig_arr, start_arr, end_arr) =
            get_join_col_arrays(batch, (columns.0, columns.1, columns.2))?;
        let mut query_contigs = Vec::with_capacity(batch.num_rows());
        let mut query_starts = Vec::with_capacity(batch.num_rows());
        let mut query_ends = Vec::with_capacity(batch.num_rows());
        let strict = filter_op == FilterOp::Strict;

        for i in 0..batch.num_rows() {
            let contig = contig_arr.value(i);
            let mut start = i64::from(start_arr.value(i)?);
            let mut end = i64::from(end_arr.value(i)?);
            if strict {
                start += 1;
                end -= 1;
            }
            query_contigs.push(
                self.chrom_to_id
                    .get(contig)
                    .copied()
                    .unwrap_or(MISSING_CONTIG_ID),
            );
            query_starts.push(start);
            query_ends.push(end);
        }

        Ok(EncodedQueryBatch {
            query_contigs,
            query_starts,
            query_ends,
        })
    }

    pub fn count_encoded_queries(&self, queries: &EncodedQueryBatch) -> Vec<i64> {
        queries
            .query_contigs
            .iter()
            .zip(queries.query_starts.iter())
            .zip(queries.query_ends.iter())
            .map(|((&contig_id, &start), &end)| self.count_one(contig_id, start, end))
            .collect()
    }

    fn count_one(&self, contig_id: u32, start: i64, end: i64) -> i64 {
        if contig_id == MISSING_CONTIG_ID || start > end {
            return 0;
        }
        let Some(range) = self.contigs.get(contig_id as usize) else {
            return 0;
        };
        let start_offset = range.start_offset as usize;
        let end_offset = range.end_offset as usize;
        let len = range.len as usize;
        let starts = &self.starts_sorted[start_offset..start_offset + len];
        let ends = &self.ends_sorted[end_offset..end_offset + len];
        let started = upper_bound(starts, end);
        let ended_prior = lower_bound(ends, start);
        i64::try_from(started - ended_prior).unwrap_or(i64::MAX)
    }
}

pub fn append_count_column(
    batch: &RecordBatch,
    schema: SchemaRef,
    counts: Vec<i64>,
) -> Result<RecordBatch> {
    let count_arr = Arc::new(Int64Array::from(counts));
    let mut columns = Vec::with_capacity(batch.num_columns() + 1);
    columns.extend_from_slice(batch.columns());
    columns.push(count_arr);
    RecordBatch::try_new(schema, columns).map_err(Into::into)
}

pub fn build_output_batch(
    batch: &RecordBatch,
    output_schema: SchemaRef,
    sources: &[OutputColumnSource],
    counts: Vec<i64>,
) -> Result<RecordBatch> {
    let mut count_array = None::<ArrayRef>;
    let mut counts = Some(counts);
    let columns = sources
        .iter()
        .map(|source| match source {
            OutputColumnSource::Input(idx) => Ok(batch.column(*idx).clone()),
            OutputColumnSource::Count => Ok(count_array
                .get_or_insert_with(|| {
                    Arc::new(Int64Array::from(counts.take().unwrap_or_default())) as ArrayRef
                })
                .clone()),
        })
        .collect::<Result<Vec<_>>>()?;

    if columns.is_empty() {
        RecordBatch::try_new_with_options(
            output_schema,
            columns,
            &RecordBatchOptions::new().with_row_count(Some(batch.num_rows())),
        )
        .map_err(Into::into)
    } else {
        RecordBatch::try_new(output_schema, columns).map_err(Into::into)
    }
}

fn lower_bound(values: &[i64], needle: i64) -> usize {
    let mut lo = 0usize;
    let mut hi = values.len();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if values[mid] < needle {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

fn upper_bound(values: &[i64], needle: i64) -> usize {
    let mut lo = 0usize;
    let mut hi = values.len();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if values[mid] <= needle {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use datafusion::arrow::array::{Int32Array, RecordBatch, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    use super::*;

    fn batch(rows: &[(&str, i32, i32)]) -> Result<RecordBatch> {
        RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("contig", DataType::Utf8, false),
                Field::new("pos_start", DataType::Int32, false),
                Field::new("pos_end", DataType::Int32, false),
            ])),
            vec![
                Arc::new(StringArray::from_iter_values(rows.iter().map(|r| r.0))),
                Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.1))),
                Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.2))),
            ],
        )
        .map_err(Into::into)
    }

    #[test]
    fn rank_index_counts_duplicates_points_and_missing_contigs() -> Result<()> {
        let left = batch(&[
            ("chr1", 10, 100),
            ("chr1", 20, 30),
            ("chr1", 20, 30),
            ("chr1", 50, 50),
            ("chr1", 101, 110),
            ("chr2", 1, 5),
        ])?;
        let right = batch(&[
            ("chr1", 25, 25),
            ("chr1", 50, 50),
            ("chr1", 100, 101),
            ("chr3", 1, 100),
        ])?;

        let index = CountOverlapsRankIndex::build_from_batches(
            &[left],
            ("contig", "pos_start", "pos_end"),
        )?;
        let queries =
            index.encode_query_batch(&right, ("contig", "pos_start", "pos_end"), FilterOp::Weak)?;
        assert_eq!(index.count_encoded_queries(&queries), vec![3, 2, 2, 0]);
        Ok(())
    }

    #[test]
    fn rank_index_strict_empty_query_returns_zero() -> Result<()> {
        let left = batch(&[("a", 10, 20), ("a", 20, 30), ("a", 15, 25)])?;
        let right = batch(&[("a", 20, 20), ("a", 19, 21)])?;

        let index = CountOverlapsRankIndex::build_from_batches(
            &[left],
            ("contig", "pos_start", "pos_end"),
        )?;
        let queries = index.encode_query_batch(
            &right,
            ("contig", "pos_start", "pos_end"),
            FilterOp::Strict,
        )?;
        assert_eq!(index.count_encoded_queries(&queries), vec![0, 3]);
        Ok(())
    }
}
