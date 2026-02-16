use std::collections::HashMap;

use std::collections::HashMap as StdHashMap;

use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::arrow::datatypes::{DataType, SchemaRef, UInt32Type};
use serde::Deserialize;

/// Metadata key for BAM reference sequences stored as JSON in Arrow schema metadata.
/// Inlined from `datafusion-bio-format-core` to avoid a runtime dependency.
const BAM_REFERENCE_SEQUENCES_KEY: &str = "bio.bam.reference_sequences";

/// Reference sequence metadata (contig name + length) from BAM header.
/// Inlined from `datafusion-bio-format-core` to avoid a runtime dependency.
#[derive(Deserialize)]
struct ReferenceSequenceMetadata {
    name: String,
    length: usize,
    #[serde(default)]
    #[allow(dead_code)]
    other_fields: StdHashMap<String, String>,
}

use crate::cigar;
use crate::filter::ReadFilter;
use crate::schema;

/// Pre-computed column indices for the BAM schema.
///
/// Avoids repeated linear `schema.index_of()` lookups per batch.
pub struct ColumnIndices {
    pub chrom: usize,
    pub start: usize,
    pub flags: usize,
    pub cigar: usize,
    pub mapq: usize,
    /// Whether the cigar column is binary (DataType::Binary) or string (DataType::Utf8).
    pub binary_cigar: bool,
}

impl ColumnIndices {
    pub fn from_schema(s: &SchemaRef) -> Self {
        let cigar = s.index_of(schema::COL_CIGAR).expect("missing cigar column");
        let binary_cigar = matches!(s.field(cigar).data_type(), DataType::Binary);
        Self {
            chrom: s.index_of(schema::COL_CHROM).expect("missing chrom column"),
            start: s.index_of(schema::COL_START).expect("missing start column"),
            flags: s.index_of(schema::COL_FLAGS).expect("missing flags column"),
            cigar,
            mapq: s
                .index_of(schema::COL_MAPPING_QUALITY)
                .expect("missing mapping_quality column"),
            binary_cigar,
        }
    }
}

/// Events for a single contig, stored as an unsorted flat list of (position, delta).
///
/// Uses a Vec for O(1) appends and cache-friendly sequential writes.
/// For sorted BAM input, the list is nearly sorted, so the final
/// `sort_unstable` (timsort) runs in ~O(n).
#[derive(Default)]
pub struct ContigEvents {
    pub events: Vec<(u32, i32)>,
    pub min_pos: Option<u32>,
    pub max_pos: Option<u32>,
}

impl ContigEvents {
    pub fn push(&mut self, pos: u32, delta: i32) {
        self.events.push((pos, delta));
        self.min_pos = Some(self.min_pos.map_or(pos, |m| m.min(pos)));
        self.max_pos = Some(self.max_pos.map_or(pos, |m| m.max(pos)));
    }

    /// Update min/max bounds from events added since `prev_len`.
    pub fn update_bounds_from(&mut self, prev_len: usize) {
        for &(pos, _) in &self.events[prev_len..] {
            self.min_pos = Some(self.min_pos.map_or(pos, |m| m.min(pos)));
            self.max_pos = Some(self.max_pos.map_or(pos, |m| m.max(pos)));
        }
    }

    /// Merge bounds from another `ContigEvents`.
    pub fn merge_bounds(&mut self, other: &ContigEvents) {
        if let Some(other_min) = other.min_pos {
            self.min_pos = Some(self.min_pos.map_or(other_min, |m| m.min(other_min)));
        }
        if let Some(other_max) = other.max_pos {
            self.max_pos = Some(self.max_pos.map_or(other_max, |m| m.max(other_max)));
        }
    }
}

/// Accumulates coverage events from a RecordBatch into per-contig event maps.
///
/// Processes rows from the input batch, applying the read filter, parsing CIGAR
/// strings, and accumulating the resulting events.
pub fn process_batch(
    batch: &RecordBatch,
    filter: &ReadFilter,
    contig_events: &mut HashMap<String, ContigEvents>,
    col_idx: &ColumnIndices,
) {
    let chrom_arr = batch.column(col_idx.chrom).as_string::<i32>();
    let start_arr = batch.column(col_idx.start).as_primitive::<UInt32Type>();
    let flags_arr = batch.column(col_idx.flags).as_primitive::<UInt32Type>();
    let mapq_arr = batch.column(col_idx.mapq).as_primitive::<UInt32Type>();

    if col_idx.binary_cigar {
        let cigar_arr = batch.column(col_idx.cigar).as_binary::<i32>();
        for row in 0..batch.num_rows() {
            if chrom_arr.is_null(row) || start_arr.is_null(row) {
                continue;
            }
            let chrom = chrom_arr.value(row);
            let start = start_arr.value(row);
            let flags = flags_arr.value(row);
            let cigar_bytes = cigar_arr.value(row);
            let mapq = mapq_arr.value(row);

            if cigar_bytes.is_empty() {
                continue;
            }
            if !filter.passes(flags, mapq) {
                continue;
            }

            if !contig_events.contains_key(chrom) {
                contig_events.insert(chrom.to_string(), ContigEvents::default());
            }
            let contig_entry = contig_events.get_mut(chrom).unwrap();
            let prev_len = contig_entry.events.len();
            cigar::apply_binary_cigar_to_event_list(start, cigar_bytes, &mut contig_entry.events);
            contig_entry.update_bounds_from(prev_len);
        }
    } else {
        let cigar_arr = batch.column(col_idx.cigar).as_string::<i32>();
        for row in 0..batch.num_rows() {
            if chrom_arr.is_null(row) || start_arr.is_null(row) {
                continue;
            }
            let chrom = chrom_arr.value(row);
            let start = start_arr.value(row);
            let flags = flags_arr.value(row);
            let cigar_str = cigar_arr.value(row);
            let mapq = mapq_arr.value(row);

            if cigar_str == "*" {
                continue;
            }
            if !filter.passes(flags, mapq) {
                continue;
            }

            if !contig_events.contains_key(chrom) {
                contig_events.insert(chrom.to_string(), ContigEvents::default());
            }
            let contig_entry = contig_events.get_mut(chrom).unwrap();
            let prev_len = contig_entry.events.len();
            cigar::apply_cigar_to_event_list(start, cigar_str, &mut contig_entry.events);
            contig_entry.update_bounds_from(prev_len);
        }
    }
}

/// Dense depth accumulator for a single contig.
///
/// Uses a flat `Vec<i32>` indexed by genomic position for O(1) writes
/// and excellent cache locality (mosdepth-style).
///
/// Tracks the range of positions actually written to, so finalization
/// can scan only the touched region instead of the entire contig length.
pub struct DenseContigDepth {
    pub depth: Vec<i32>,
    min_touched: Option<usize>,
    max_touched: Option<usize>,
}

impl DenseContigDepth {
    pub fn new(length: usize) -> Self {
        Self {
            depth: vec![0i32; length + 1], // +1 for end sentinel
            min_touched: None,
            max_touched: None,
        }
    }

    /// Widen the tracked bounds to include `[lo, hi]`.
    pub fn update_bounds(&mut self, lo: usize, hi: usize) {
        self.min_touched = Some(self.min_touched.map_or(lo, |m| m.min(lo)));
        self.max_touched = Some(self.max_touched.map_or(hi, |m| m.max(hi)));
    }

    /// Return the slice of the depth array that was actually modified,
    /// or an empty slice if nothing was touched.
    pub fn touched_range(&self) -> &[i32] {
        match (self.min_touched, self.max_touched) {
            (Some(lo), Some(hi)) => {
                let end = (hi + 1).min(self.depth.len());
                &self.depth[lo..end]
            }
            _ => &[],
        }
    }

    /// Return the start offset of the touched region (0 if untouched).
    pub fn touched_start(&self) -> usize {
        self.min_touched.unwrap_or(0)
    }
}

/// Extract contig lengths from Arrow schema metadata.
///
/// Looks for the `bio.bam.reference_sequences` key in schema metadata,
/// parses the JSON array of `ReferenceSequenceMetadata`, and returns
/// a map of contig name -> length.
pub fn extract_contig_lengths(schema: &SchemaRef) -> Option<HashMap<String, usize>> {
    let metadata = schema.metadata();
    let json_str = metadata.get(BAM_REFERENCE_SEQUENCES_KEY)?;
    let refs: Vec<ReferenceSequenceMetadata> = serde_json::from_str(json_str).ok()?;
    let map = refs.into_iter().map(|r| (r.name, r.length)).collect();
    Some(map)
}

/// Process a RecordBatch, writing CIGAR events directly to dense depth arrays.
///
/// For each row, applies the CIGAR string directly to the contig's depth array
/// via `apply_cigar_to_depth()`. Returns the contigs seen in this batch in order
/// of first appearance.
pub fn process_batch_dense(
    batch: &RecordBatch,
    filter: &ReadFilter,
    contig_depths: &mut HashMap<String, DenseContigDepth>,
    contig_lengths: &HashMap<String, usize>,
    col_idx: &ColumnIndices,
) -> Vec<String> {
    let chrom_arr = batch.column(col_idx.chrom).as_string::<i32>();
    let start_arr = batch.column(col_idx.start).as_primitive::<UInt32Type>();
    let flags_arr = batch.column(col_idx.flags).as_primitive::<UInt32Type>();
    let mapq_arr = batch.column(col_idx.mapq).as_primitive::<UInt32Type>();

    let mut seen_contigs: Vec<String> = Vec::new();

    /// Ensure we have a depth array for this contig, returning whether it's available.
    fn ensure_contig(
        chrom: &str,
        contig_depths: &mut HashMap<String, DenseContigDepth>,
        contig_lengths: &HashMap<String, usize>,
        seen_contigs: &mut Vec<String>,
    ) -> bool {
        if !contig_depths.contains_key(chrom) {
            if let Some(&len) = contig_lengths.get(chrom) {
                let chrom_owned = chrom.to_string();
                contig_depths.insert(chrom_owned.clone(), DenseContigDepth::new(len));
                seen_contigs.push(chrom_owned);
            } else {
                return false;
            }
        } else if !seen_contigs.iter().any(|s| s.as_str() == chrom) {
            seen_contigs.push(chrom.to_string());
        }
        true
    }

    if col_idx.binary_cigar {
        let cigar_arr = batch.column(col_idx.cigar).as_binary::<i32>();
        for row in 0..batch.num_rows() {
            if chrom_arr.is_null(row) || start_arr.is_null(row) {
                continue;
            }
            let chrom = chrom_arr.value(row);
            let start = start_arr.value(row);
            let flags = flags_arr.value(row);
            let cigar_bytes = cigar_arr.value(row);
            let mapq = mapq_arr.value(row);

            if cigar_bytes.is_empty() {
                continue;
            }
            if !filter.passes(flags, mapq) {
                continue;
            }
            if !ensure_contig(chrom, contig_depths, contig_lengths, &mut seen_contigs) {
                continue;
            }
            let entry = contig_depths.get_mut(chrom).unwrap();
            if let Some((lo, hi)) =
                cigar::apply_binary_cigar_to_depth(start, cigar_bytes, &mut entry.depth)
            {
                entry.update_bounds(lo, hi);
            }
        }
    } else {
        let cigar_arr = batch.column(col_idx.cigar).as_string::<i32>();
        for row in 0..batch.num_rows() {
            if chrom_arr.is_null(row) || start_arr.is_null(row) {
                continue;
            }
            let chrom = chrom_arr.value(row);
            let start = start_arr.value(row);
            let flags = flags_arr.value(row);
            let cigar_str = cigar_arr.value(row);
            let mapq = mapq_arr.value(row);

            if cigar_str == "*" {
                continue;
            }
            if !filter.passes(flags, mapq) {
                continue;
            }
            if !ensure_contig(chrom, contig_depths, contig_lengths, &mut seen_contigs) {
                continue;
            }
            let entry = contig_depths.get_mut(chrom).unwrap();
            if let Some((lo, hi)) = cigar::apply_cigar_to_depth(start, cigar_str, &mut entry.depth)
            {
                entry.update_bounds(lo, hi);
            }
        }
    }

    seen_contigs
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{BinaryArray, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    /// Encode a CIGAR op as packed LE u32 bytes.
    fn encode_op(len: u32, code: u32) -> [u8; 4] {
        ((len << 4) | code).to_le_bytes()
    }

    fn make_batch_binary(
        chroms: Vec<Option<&str>>,
        starts: Vec<Option<u32>>,
        ends: Vec<Option<u32>>,
        flags: Vec<u32>,
        cigars: Vec<&[u8]>,
        mapqs: Vec<u32>,
    ) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("start", DataType::UInt32, true),
            Field::new("end", DataType::UInt32, true),
            Field::new("flags", DataType::UInt32, false),
            Field::new("cigar", DataType::Binary, false),
            Field::new("mapping_quality", DataType::UInt32, false),
        ]));

        let chrom_arr: StringArray = chroms.into_iter().collect();
        let start_arr: UInt32Array = starts.into_iter().collect();
        let end_arr: UInt32Array = ends.into_iter().collect();
        let flags_arr = UInt32Array::from(flags);
        let cigar_arr = BinaryArray::from(cigars);
        let mapq_arr = UInt32Array::from(mapqs);

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(chrom_arr),
                Arc::new(start_arr),
                Arc::new(end_arr),
                Arc::new(flags_arr),
                Arc::new(cigar_arr),
                Arc::new(mapq_arr),
            ],
        )
        .unwrap()
    }

    fn make_batch(
        chroms: Vec<Option<&str>>,
        starts: Vec<Option<u32>>,
        ends: Vec<Option<u32>>,
        flags: Vec<u32>,
        cigars: Vec<&str>,
        mapqs: Vec<u32>,
    ) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("start", DataType::UInt32, true),
            Field::new("end", DataType::UInt32, true),
            Field::new("flags", DataType::UInt32, false),
            Field::new("cigar", DataType::Utf8, false),
            Field::new("mapping_quality", DataType::UInt32, false),
        ]));

        let chrom_arr: StringArray = chroms.into_iter().collect();
        let start_arr: UInt32Array = starts.into_iter().collect();
        let end_arr: UInt32Array = ends.into_iter().collect();
        let flags_arr = UInt32Array::from(flags);
        let cigar_arr = StringArray::from(cigars);
        let mapq_arr = UInt32Array::from(mapqs);

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(chrom_arr),
                Arc::new(start_arr),
                Arc::new(end_arr),
                Arc::new(flags_arr),
                Arc::new(cigar_arr),
                Arc::new(mapq_arr),
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_single_read() {
        let batch = make_batch(
            vec![Some("chr1")],
            vec![Some(100)],
            vec![Some(110)],
            vec![0],
            vec!["10M"],
            vec![60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 1);
        let events = &contig_events["chr1"].events;
        assert!(events.contains(&(100, 1)));
        assert!(events.contains(&(110, -1)));
    }

    #[test]
    fn test_overlapping_reads() {
        let batch = make_batch(
            vec![Some("chr1"), Some("chr1")],
            vec![Some(0), Some(5)],
            vec![Some(10), Some(15)],
            vec![0, 0],
            vec!["10M", "10M"],
            vec![60, 60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        let events = &contig_events["chr1"].events;
        // 4 events: +1@0, -1@10 from read1; +1@5, -1@15 from read2
        assert_eq!(events.len(), 4);
    }

    #[test]
    fn test_multi_contig() {
        let batch = make_batch(
            vec![Some("chr1"), Some("chr2")],
            vec![Some(0), Some(100)],
            vec![Some(10), Some(110)],
            vec![0, 0],
            vec!["10M", "10M"],
            vec![60, 60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 2);
        assert!(contig_events.contains_key("chr1"));
        assert!(contig_events.contains_key("chr2"));
    }

    #[test]
    fn test_filtered_reads_excluded() {
        let batch = make_batch(
            vec![Some("chr1"), Some("chr1")],
            vec![Some(0), Some(100)],
            vec![Some(10), Some(110)],
            vec![0, 4], // second read is unmapped
            vec!["10M", "10M"],
            vec![60, 60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 1);
        let events = &contig_events["chr1"].events;
        assert_eq!(events.len(), 2); // only events from first read
    }

    #[test]
    fn test_null_chroms_skipped() {
        let batch = make_batch(
            vec![None, Some("chr1")],
            vec![Some(0), Some(100)],
            vec![Some(10), Some(110)],
            vec![0, 0],
            vec!["10M", "10M"],
            vec![60, 60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 1);
        assert!(contig_events.contains_key("chr1"));
    }

    #[test]
    fn test_star_cigar_skipped() {
        let batch = make_batch(
            vec![Some("chr1"), Some("chr1")],
            vec![Some(0), Some(100)],
            vec![Some(0), Some(110)],
            vec![0, 0],
            vec!["*", "10M"],
            vec![0, 60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 1);
        let events = &contig_events["chr1"].events;
        assert!(events.contains(&(100, 1)));
        assert!(events.contains(&(110, -1)));
    }

    // --- Tests for dense path ---

    #[test]
    fn test_extract_contig_lengths() {
        let metadata = HashMap::from([(
            BAM_REFERENCE_SEQUENCES_KEY.to_string(),
            r#"[{"name":"chr1","length":1000},{"name":"chr2","length":500}]"#.to_string(),
        )]);
        let schema = Arc::new(
            Schema::new(vec![Field::new("dummy", DataType::Utf8, false)]).with_metadata(metadata),
        );
        let lengths = extract_contig_lengths(&schema).unwrap();
        assert_eq!(lengths["chr1"], 1000);
        assert_eq!(lengths["chr2"], 500);
    }

    #[test]
    fn test_extract_contig_lengths_missing_key() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "dummy",
            DataType::Utf8,
            false,
        )]));
        assert!(extract_contig_lengths(&schema).is_none());
    }

    #[test]
    fn test_process_batch_dense_single_read() {
        let batch = make_batch(
            vec![Some("chr1")],
            vec![Some(100)],
            vec![Some(110)],
            vec![0],
            vec!["10M"],
            vec![60],
        );
        let contig_lengths = HashMap::from([("chr1".to_string(), 1000usize)]);
        let filter = ReadFilter::default();
        let mut contig_depths = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        let seen = process_batch_dense(
            &batch,
            &filter,
            &mut contig_depths,
            &contig_lengths,
            &col_idx,
        );

        assert_eq!(seen, vec!["chr1"]);
        let depth = &contig_depths["chr1"].depth;
        assert_eq!(depth[100], 1);
        assert_eq!(depth[110], -1);
    }

    #[test]
    fn test_process_batch_dense_filtered_reads() {
        let batch = make_batch(
            vec![Some("chr1"), Some("chr1")],
            vec![Some(0), Some(100)],
            vec![Some(10), Some(110)],
            vec![0, 4], // second is unmapped
            vec!["10M", "10M"],
            vec![60, 60],
        );
        let contig_lengths = HashMap::from([("chr1".to_string(), 1000usize)]);
        let filter = ReadFilter::default();
        let mut contig_depths = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch_dense(
            &batch,
            &filter,
            &mut contig_depths,
            &contig_lengths,
            &col_idx,
        );

        let depth = &contig_depths["chr1"].depth;
        assert_eq!(depth[0], 1);
        assert_eq!(depth[10], -1);
        // Filtered read at 100 should NOT be present
        assert_eq!(depth[100], 0);
    }

    // --- Tests for binary CIGAR path ---

    #[test]
    fn test_binary_single_read() {
        let cigar_10m = encode_op(10, 0); // 10M
        let batch = make_batch_binary(
            vec![Some("chr1")],
            vec![Some(100)],
            vec![Some(110)],
            vec![0],
            vec![&cigar_10m],
            vec![60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        assert!(col_idx.binary_cigar);
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 1);
        let events = &contig_events["chr1"].events;
        assert!(events.contains(&(100, 1)));
        assert!(events.contains(&(110, -1)));
    }

    #[test]
    fn test_binary_empty_cigar_skipped() {
        let cigar_10m = encode_op(10, 0); // 10M
        let batch = make_batch_binary(
            vec![Some("chr1"), Some("chr1")],
            vec![Some(0), Some(100)],
            vec![Some(0), Some(110)],
            vec![0, 0],
            vec![&[] as &[u8], &cigar_10m],
            vec![0, 60],
        );
        let filter = ReadFilter::default();
        let mut contig_events = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        process_batch(&batch, &filter, &mut contig_events, &col_idx);

        assert_eq!(contig_events.len(), 1);
        let events = &contig_events["chr1"].events;
        assert!(events.contains(&(100, 1)));
        assert!(events.contains(&(110, -1)));
    }

    #[test]
    fn test_binary_dense_single_read() {
        let cigar_10m = encode_op(10, 0); // 10M
        let batch = make_batch_binary(
            vec![Some("chr1")],
            vec![Some(100)],
            vec![Some(110)],
            vec![0],
            vec![&cigar_10m],
            vec![60],
        );
        let contig_lengths = HashMap::from([("chr1".to_string(), 1000usize)]);
        let filter = ReadFilter::default();
        let mut contig_depths = HashMap::new();
        let col_idx = ColumnIndices::from_schema(&batch.schema());
        assert!(col_idx.binary_cigar);
        let seen = process_batch_dense(
            &batch,
            &filter,
            &mut contig_depths,
            &contig_lengths,
            &col_idx,
        );

        assert_eq!(seen, vec!["chr1"]);
        let depth = &contig_depths["chr1"].depth;
        assert_eq!(depth[100], 1);
        assert_eq!(depth[110], -1);
    }
}
