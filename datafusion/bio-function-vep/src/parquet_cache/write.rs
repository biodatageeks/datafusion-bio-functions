//! Parquet writer for the variation cache shard.
//!
//! The point-lookup entities (variation, translation_sift) are written with a
//! lookup-optimized `WriterProperties`: **no dictionary** (a dictionary forces
//! loading the whole row-group dictionary to read any single row — the take-time
//! "dictionary tax"), zstd(3) to recover size, small ~4 KiB / 512-row data pages,
//! and page-level statistics so the footer carries the `ColumnIndex` +
//! `OffsetIndex` the read-side `PageDir` resolves against. Rows are declared
//! sorted by `(tier, start)` so `start` is monotonic within each tier run — the
//! two-run structure the `PageDir` binary search relies on.

use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::metadata::SortingColumn;
use parquet::file::properties::{EnabledStatistics, WriterProperties};

/// Row-group size for variation shards (rows). Matches the Lance/Parquet
/// reference build; large enough to amortize footer metadata, small enough to
/// keep row-group pruning useful.
const VARIATION_ROW_GROUP_ROWS: usize = 1_000_000;

/// Lookup-optimized `WriterProperties` for a variation Parquet shard.
///
/// `schema` is the *output* schema (post projection/tier/Boolean/AF encoding);
/// the `tier` and `start` columns are used as the declared sort order when
/// present.
pub fn variation_writer_properties(schema: &SchemaRef) -> WriterProperties {
    // Declared sort: (tier, start). `SortingColumn` indexes are into the leaf
    // column order; `tier` and `start` are top-level primitives so their leaf
    // index equals their field index.
    let sorting: Vec<SortingColumn> = ["tier", "start"]
        .iter()
        .filter_map(|name| schema.index_of(name).ok())
        .map(|idx| SortingColumn {
            column_idx: idx as i32,
            descending: false,
            nulls_first: false,
        })
        .collect();

    WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(3).unwrap()))
        // Artifact #1: no dictionary → no per-take dictionary load. zstd recovers
        // the compression (the no-dict file is actually smaller).
        .set_dictionary_enabled(false)
        // Small pages = fine-grained page index → cheap point-lookup resolution.
        .set_data_page_size_limit(4 * 1024)
        .set_data_page_row_count_limit(512)
        // Page-level statistics emit ColumnIndex + OffsetIndex in the footer,
        // which the read-side `PageDir` uses as the position→page directory.
        .set_statistics_enabled(EnabledStatistics::Page)
        .set_max_row_group_row_count(Some(VARIATION_ROW_GROUP_ROWS))
        .set_sorting_columns(Some(sorting))
        .build()
}

/// Write already-encoded variation batches to a single Parquet shard at `path`
/// using [`variation_writer_properties`]. Returns the total row count. This is a
/// thin, generic writer over batches that have already been projected/encoded to
/// the output schema (Boolean flags, 2-array AF, etc.); the per-contig
/// source→tier→encode pipeline is built on top of it.
pub fn write_variation_parquet(path: &std::path::Path, batches: &[RecordBatch]) -> Result<usize> {
    if batches.is_empty() {
        return Ok(0);
    }
    let schema = batches[0].schema();
    let props = variation_writer_properties(&schema);
    let file = std::fs::File::create(path).map_err(|e| {
        DataFusionError::Execution(format!("create parquet '{}': {e}", path.display()))
    })?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| DataFusionError::Execution(format!("open ArrowWriter: {e}")))?;
    let mut rows = 0usize;
    for batch in batches {
        writer
            .write(batch)
            .map_err(|e| DataFusionError::Execution(format!("write batch: {e}")))?;
        rows += batch.num_rows();
    }
    writer
        .close()
        .map_err(|e| DataFusionError::Execution(format!("close ArrowWriter: {e}")))?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::schema::types::ColumnPath;
    use std::sync::Arc;

    #[test]
    fn variation_writer_properties_are_no_dict_small_page_indexed() {
        let schema: SchemaRef = Arc::new(Schema::new(vec![
            Field::new("tier", DataType::Int8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let props = variation_writer_properties(&schema);

        // No dictionary anywhere (the take-time tax the whole design avoids).
        assert!(!props.dictionary_enabled(&ColumnPath::from("start")));
        assert!(!props.dictionary_enabled(&ColumnPath::from("allele_string")));
        // Small, lookup-friendly pages.
        assert_eq!(props.data_page_row_count_limit(), 512);
        // Page index present (ColumnIndex + OffsetIndex) for the footer PageDir.
        assert!(matches!(
            props.statistics_enabled(&ColumnPath::from("start")),
            EnabledStatistics::Page
        ));
    }

    #[test]
    fn variation_writer_properties_declare_tier_start_sort() {
        let schema: SchemaRef = Arc::new(Schema::new(vec![
            Field::new("tier", DataType::Int8, false),
            Field::new("start", DataType::UInt32, false),
        ]));
        let props = variation_writer_properties(&schema);
        let sorting = props.sorting_columns().expect("sorting columns set");
        assert_eq!(sorting.len(), 2);
        assert_eq!(sorting[0].column_idx, 0); // tier
        assert_eq!(sorting[1].column_idx, 1); // start
    }

    /// Full core mechanism end-to-end on a synthetic fixture: encode (Boolean
    /// flag + 2-array AF) → write no-dict/page-indexed Parquet → footer PageDir
    /// resolve → exact-offset take via the coalescing async reader → reconstruct
    /// the AF string via `%.4g`. Proves the pipeline is lossless.
    #[tokio::test(flavor = "multi_thread")]
    async fn variation_roundtrip_encode_write_pageindex_take() {
        use crate::parquet_cache::encode::{encode_af_2array, format_g4, presence_boolean};
        use crate::parquet_cache::page_dir::{
            CoalescingAsyncReader, IoCounters, PageDir, selection_from_offsets,
            selection_from_ranges,
        };
        use datafusion::arrow::array::{
            Array, BooleanArray, Float32Array, Int8Array, ListArray, StringArray, UInt32Array,
        };
        use datafusion::arrow::datatypes::{DataType, Field, Schema};
        use futures::TryStreamExt;
        use parquet::arrow::ProjectionMask;
        use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};
        use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;
        use std::collections::HashSet;
        use std::sync::Arc;

        const N_POPS: usize = 11;
        // 11 populations, two intentionally missing (empty), all tokens %.4g-clean.
        const AF: &str =
            "A:0.08806|A:0.0367|A:2.682e-05||A:0.9969|A:0.5|A:0.006912||A:9.911e-05|A:0.4253|A:1";

        // Fixture: 2000 strictly-ascending starts -> many small pages.
        let n = 2000usize;
        let starts: Vec<u32> = (0..n).map(|i| (i as u32) * 1000 + 100).collect();
        let failed_src: Vec<Option<i8>> = (0..n)
            .map(|i| if i % 7 == 0 { Some(1) } else { None })
            .collect();

        let failed_int = Int8Array::from(failed_src.clone());
        let failed_bool = presence_boolean(&failed_int, "failed").unwrap();
        let af = encode_af_2array(&StringArray::from(vec![AF; n]), N_POPS).unwrap();

        let schema = Arc::new(Schema::new(vec![
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Boolean, false),
            Field::new("af_gnomadg_alleles", af.alleles.data_type().clone(), true),
            Field::new("af_gnomadg_freqs", af.freqs.data_type().clone(), true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(starts.clone())),
                Arc::new(StringArray::from(vec!["A"; n])),
                Arc::new(failed_bool),
                Arc::new(af.alleles),
                Arc::new(af.freqs),
            ],
        )
        .unwrap();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_path_buf();
        assert_eq!(write_variation_parquet(&path, &[batch]).unwrap(), n);

        // Reopen with the page index; build the PageDir from the footer.
        let stdf = std::fs::File::open(&path).unwrap();
        // `with_page_index_policy` is the non-deprecated API, but `PageIndexPolicy`
        // is private in parquet 58 — the deprecated setter is the only usable one.
        #[allow(deprecated)]
        let meta =
            ArrowReaderMetadata::load(&stdf, ArrowReaderOptions::new().with_page_index(true))
                .unwrap();
        let pq_schema = meta.parquet_schema();
        let start_leaf = pq_schema
            .columns()
            .iter()
            .position(|c| c.name() == "start")
            .unwrap();
        let pd = PageDir::build(&meta, start_leaf).unwrap();

        // Probe a few known positions.
        let probes: Vec<u32> = vec![starts[3], starts[100], starts[1999]];
        let ranges = pd.resolve_ranges(&probes);
        assert!(!ranges.is_empty());

        // Phase 2: start-only read of candidate pages -> exact offsets.
        let start_mask = ProjectionMask::leaves(pq_schema, [start_leaf]);
        let counters = IoCounters::new();
        let reader = CoalescingAsyncReader::new(
            tokio::fs::File::open(&path).await.unwrap(),
            counters.clone(),
            64 * 1024,
        );
        let mut stream = ParquetRecordBatchStreamBuilder::new_with_metadata(reader, meta.clone())
            .with_projection(start_mask)
            .with_row_selection(selection_from_ranges(&ranges))
            .with_batch_size(8192)
            .build()
            .unwrap();
        let probe_set: HashSet<u32> = probes.iter().copied().collect();
        let mut offsets = Vec::new();
        let mut cursor = ranges.iter().flat_map(|&(a, b)| a..b);
        while let Some(b) = stream.try_next().await.unwrap() {
            let sa = b.column(0).as_any().downcast_ref::<UInt32Array>().unwrap();
            for &v in sa.values() {
                let off = cursor.next().unwrap();
                if probe_set.contains(&v) {
                    offsets.push(off);
                }
            }
        }
        assert_eq!(offsets.len(), probes.len());

        // Phase 3: full-projection payload take at the exact offsets.
        let reader = CoalescingAsyncReader::new(
            tokio::fs::File::open(&path).await.unwrap(),
            counters,
            64 * 1024,
        );
        let mut stream = ParquetRecordBatchStreamBuilder::new_with_metadata(reader, meta.clone())
            .with_row_selection(selection_from_offsets(&offsets))
            .with_batch_size(8192)
            .build()
            .unwrap();
        let mut taken: Vec<RecordBatch> = Vec::new();
        while let Some(b) = stream.try_next().await.unwrap() {
            taken.push(b);
        }
        let out = datafusion::arrow::compute::concat_batches(&taken[0].schema(), &taken).unwrap();
        assert_eq!(out.num_rows(), probes.len());

        let out_start = out
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let out_failed = out
            .column(2)
            .as_any()
            .downcast_ref::<BooleanArray>()
            .unwrap();
        let out_alleles = out.column(3).as_any().downcast_ref::<ListArray>().unwrap();
        let out_freqs = out.column(4).as_any().downcast_ref::<ListArray>().unwrap();

        for r in 0..out.num_rows() {
            let orig_idx = starts
                .iter()
                .position(|&x| x == out_start.value(r))
                .unwrap();
            assert_eq!(out_failed.value(r), failed_src[orig_idx].is_some());

            // Reconstruct the AF string (allele-major, N_POPS positional freqs).
            let al = out_alleles.value(r);
            let al = al.as_any().downcast_ref::<StringArray>().unwrap();
            let fr = out_freqs.value(r);
            let fr = fr.as_any().downcast_ref::<Float32Array>().unwrap();
            let mut pops_out: Vec<String> = Vec::with_capacity(N_POPS);
            for p in 0..N_POPS {
                let mut parts: Vec<String> = Vec::new();
                for j in 0..al.len() {
                    let idx = j * N_POPS + p;
                    if !fr.is_null(idx) {
                        parts.push(format!("{}:{}", al.value(j), format_g4(fr.value(idx))));
                    }
                }
                pops_out.push(parts.join(","));
            }
            assert_eq!(pops_out.join("|"), AF, "AF round-trip mismatch at row {r}");
        }
    }
}
