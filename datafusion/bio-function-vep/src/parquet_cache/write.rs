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

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{Array, ArrayRef, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::metadata::SortingColumn;
use parquet::file::properties::{EnabledStatistics, WriterProperties};

use crate::lance_cache::af_bundle::{AF_GROUPS, af_column_order, concat_group};
use crate::lance_cache::build::variation_projected_schema;
use crate::parquet_cache::encode::{
    AfArrays, dedup_variation_name, encode_af_2array, presence_boolean,
};

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
pub fn write_variation_parquet(path: &Path, batches: &[RecordBatch]) -> Result<usize> {
    if batches.is_empty() {
        return Ok(0);
    }
    let mut writer = VariationParquetShardWriter::create(path, batches[0].schema())?;
    for batch in batches {
        writer.write(batch)?;
    }
    writer.finish()
}

/// The three binary variation flags stored as non-nullable presence `Boolean`.
pub const VARIATION_FLAG_COLUMNS: [&str; 3] = ["failed", "somatic", "phenotype_or_disease"];

/// Physical output schema for a Parquet variation shard.
///
/// Starts from the Lance projected schema (non-AF required columns with
/// `start`/`end` as `UInt32`, the 27 AF `Utf8` columns, and the derived `tier`),
/// then: re-types the binary flags to non-nullable `Boolean`, drops the 27 AF
/// string columns, and appends the 3 struct-of-arrays AF pairs
/// (`<grp>_alleles: List<Utf8>`, `<grp>_freqs: List<Float32>`). This is the exact
/// physical shape [`crate::parquet_cache::variation_lookup`] reconstructs from.
pub fn variation_output_schema(source_schema: &Schema, source_type: &str) -> Result<SchemaRef> {
    let projected = variation_projected_schema(source_schema, source_type)?;
    let af_cols: HashSet<&str> = af_column_order().into_iter().collect();
    let mut fields: Vec<Field> = Vec::new();
    for field in projected.fields() {
        let name = field.name().as_str();
        if af_cols.contains(name) {
            continue; // 27 AF Utf8 columns → replaced by the 2-array pairs below.
        }
        if VARIATION_FLAG_COLUMNS.contains(&name) {
            fields.push(Field::new(name, DataType::Boolean, false));
        } else {
            fields.push(field.as_ref().clone());
        }
    }
    let (alleles_ty, freqs_ty) = af_array_datatypes()?;
    for (grp, _) in AF_GROUPS {
        fields.push(Field::new(
            format!("{grp}_alleles"),
            alleles_ty.clone(),
            true,
        ));
        fields.push(Field::new(format!("{grp}_freqs"), freqs_ty.clone(), true));
    }
    Ok(Arc::new(Schema::new_with_metadata(
        fields,
        projected.metadata().clone(),
    )))
}

/// The `List` datatypes [`encode_af_2array`] produces, derived from an empty
/// encode so the declared shard schema matches the encoded arrays exactly.
fn af_array_datatypes() -> Result<(DataType, DataType)> {
    let empty = StringArray::from(Vec::<Option<&str>>::new());
    let af = encode_af_2array(&empty, 1)?;
    Ok((af.alleles.data_type().clone(), af.freqs.data_type().clone()))
}

/// Encode one tiered variation batch (Lance projected schema — `Int8` flags, 27
/// AF `Utf8` columns) into the physical Parquet layout `out_schema` describes:
/// Boolean presence flags, 2-array AF (per [`AF_GROUPS`]), and `variation_name`
/// nulled where it equals `dbsnp_ids`. Non-AF, non-flag columns pass through by
/// name, so the caller controls their order via `out_schema`.
pub fn encode_variation_batch(
    tier_batch: &RecordBatch,
    out_schema: &SchemaRef,
) -> Result<RecordBatch> {
    // Bundle each AF group's members into a group string, then to struct-of-arrays.
    let mut groups: HashMap<&str, AfArrays> = HashMap::new();
    for (grp, members) in AF_GROUPS {
        let member_arrays = members
            .iter()
            .map(|m| string_col(tier_batch, m))
            .collect::<Result<Vec<_>>>()?;
        let concat = concat_group(&member_arrays)?;
        groups.insert(*grp, encode_af_2array(&concat, members.len())?);
    }
    let dbsnp = string_col(tier_batch, "dbsnp_ids")?;

    let mut columns: Vec<ArrayRef> = Vec::with_capacity(out_schema.fields().len());
    for field in out_schema.fields() {
        let name = field.name().as_str();
        if let Some(grp) = name.strip_suffix("_alleles") {
            columns.push(Arc::new(af_group(&groups, grp)?.alleles.clone()) as ArrayRef);
        } else if let Some(grp) = name.strip_suffix("_freqs") {
            columns.push(Arc::new(af_group(&groups, grp)?.freqs.clone()) as ArrayRef);
        } else if VARIATION_FLAG_COLUMNS.contains(&name) {
            let src = column(tier_batch, name)?;
            columns.push(Arc::new(presence_boolean(src.as_ref(), name)?) as ArrayRef);
        } else if name == "variation_name" {
            let vn = string_col(tier_batch, name)?;
            columns.push(Arc::new(dedup_variation_name(vn, dbsnp)) as ArrayRef);
        } else {
            columns.push(column(tier_batch, name)?);
        }
    }
    RecordBatch::try_new(Arc::clone(out_schema), columns).map_err(|e| {
        DataFusionError::Execution(format!("failed to assemble encoded variation batch: {e}"))
    })
}

fn af_group<'a>(groups: &'a HashMap<&str, AfArrays>, grp: &str) -> Result<&'a AfArrays> {
    groups
        .get(grp)
        .ok_or_else(|| DataFusionError::Execution(format!("unknown AF group '{grp}'")))
}

fn column(batch: &RecordBatch, name: &str) -> Result<ArrayRef> {
    let idx = batch.schema().index_of(name).map_err(|_| {
        DataFusionError::Execution(format!("variation batch missing column '{name}'"))
    })?;
    Ok(batch.column(idx).clone())
}

fn string_col<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a StringArray> {
    let idx = batch.schema().index_of(name).map_err(|_| {
        DataFusionError::Execution(format!("variation batch missing column '{name}'"))
    })?;
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            DataFusionError::Execution(format!("variation column '{name}' must be Utf8"))
        })
}

/// Streaming writer for a single variation shard. Opens the `ArrowWriter` with
/// [`variation_writer_properties`] once, then accepts already-encoded batches
/// (warm tier first, cold second) so a whole chromosome writes with bounded
/// memory while keeping the `(tier, start)` run order the read-side `PageDir`
/// depends on.
pub struct VariationParquetShardWriter {
    writer: ArrowWriter<std::fs::File>,
    rows: usize,
}

impl VariationParquetShardWriter {
    pub fn create(path: &Path, schema: SchemaRef) -> Result<Self> {
        let props = variation_writer_properties(&schema);
        let file = std::fs::File::create(path).map_err(|e| {
            DataFusionError::Execution(format!("create parquet '{}': {e}", path.display()))
        })?;
        let writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| DataFusionError::Execution(format!("open ArrowWriter: {e}")))?;
        Ok(Self { writer, rows: 0 })
    }

    pub fn write(&mut self, batch: &RecordBatch) -> Result<()> {
        self.writer
            .write(batch)
            .map_err(|e| DataFusionError::Execution(format!("write batch: {e}")))?;
        self.rows += batch.num_rows();
        Ok(())
    }

    pub fn finish(self) -> Result<usize> {
        self.writer
            .close()
            .map_err(|e| DataFusionError::Execution(format!("close ArrowWriter: {e}")))?;
        Ok(self.rows)
    }
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

    /// A synthetic Ensembl-cache source schema covering every required column,
    /// with the types the source provider yields (coords Int64, flags Int8,
    /// `minor_allele_freq` Float64, everything else Utf8).
    fn synthetic_source_schema() -> Schema {
        use crate::lance_cache::schema::VARIATION_REQUIRED_COLUMNS;
        let fields: Vec<Field> = VARIATION_REQUIRED_COLUMNS
            .iter()
            .map(|&n| {
                let dt = match n {
                    "start" | "end" => DataType::Int64,
                    "failed" | "somatic" | "phenotype_or_disease" => DataType::Int8,
                    "minor_allele_freq" => DataType::Float64,
                    _ => DataType::Utf8,
                };
                Field::new(n, dt, true)
            })
            .collect();
        Schema::new(fields)
    }

    #[test]
    fn variation_output_schema_has_boolean_flags_and_af_arrays() {
        let schema = variation_output_schema(&synthetic_source_schema(), "merged").unwrap();

        // Binary flags become non-nullable presence Booleans.
        for flag in VARIATION_FLAG_COLUMNS {
            let f = schema.field_with_name(flag).unwrap();
            assert_eq!(f.data_type(), &DataType::Boolean, "{flag}");
            assert!(!f.is_nullable(), "{flag} must be non-nullable");
        }
        // Lookup key + tier keep their types.
        assert_eq!(
            schema.field_with_name("start").unwrap().data_type(),
            &DataType::UInt32
        );
        assert_eq!(
            schema.field_with_name("tier").unwrap().data_type(),
            &DataType::Int8
        );
        // The 27 AF Utf8 columns are dropped; the 3 group array-pairs replace them.
        assert!(schema.field_with_name("AF").is_err());
        assert!(schema.field_with_name("gnomADg").is_err());
        for (grp, _) in AF_GROUPS {
            let alleles = schema.field_with_name(&format!("{grp}_alleles")).unwrap();
            let freqs = schema.field_with_name(&format!("{grp}_freqs")).unwrap();
            assert!(
                matches!(alleles.data_type(), DataType::List(_)),
                "{grp}_alleles"
            );
            assert!(
                matches!(freqs.data_type(), DataType::List(_)),
                "{grp}_freqs"
            );
        }
    }

    /// End-to-end for the encode step: build a projected-schema batch (as
    /// `transform_variation_tier_batch` would), `encode_variation_batch` it, write
    /// the shard, and read it back through the production lookup — the logical AF
    /// column and coalesced `variation_name` must reconstruct.
    #[tokio::test(flavor = "multi_thread")]
    async fn encode_variation_batch_round_trips_via_reader() {
        use crate::parquet_cache::variation_lookup::SinglePathParquetVariationLookup;
        use datafusion::arrow::array::{Float64Array, Int8Array, StringArray, UInt32Array};

        let source_type = "merged";
        let source_schema = synthetic_source_schema();
        let out_schema = variation_output_schema(&source_schema, source_type).unwrap();
        let projected = Arc::new(variation_projected_schema(&source_schema, source_type).unwrap());

        let n = 300usize;
        let starts: Vec<u32> = (0..n).map(|i| (i as u32) * 7 + 11).collect();
        let af_members: HashSet<&str> = af_column_order().into_iter().collect();

        // Build a tiered batch matching the projected schema field-by-field. AF is
        // set only on gnomADg (group `af_gnomadg`, population 0); other AF members
        // are empty. `variation_name` is null on even rows (to exercise coalesce).
        let mut columns: Vec<ArrayRef> = Vec::new();
        for field in projected.fields() {
            let name = field.name().as_str();
            let arr: ArrayRef = match field.data_type() {
                DataType::UInt32 => Arc::new(UInt32Array::from(starts.clone())),
                DataType::Int8 => {
                    if name == "failed" {
                        Arc::new(Int8Array::from(
                            (0..n)
                                .map(|i| if i % 5 == 0 { Some(1) } else { None })
                                .collect::<Vec<_>>(),
                        ))
                    } else {
                        let v = if name == "tier" { Some(0) } else { None };
                        Arc::new(Int8Array::from(vec![v; n]))
                    }
                }
                DataType::Float64 => Arc::new(Float64Array::from(vec![None::<f64>; n])),
                DataType::Utf8 => {
                    let vals: Vec<Option<String>> = (0..n)
                        .map(|i| match name {
                            "chrom" => Some("chr1".to_string()),
                            "allele_string" => Some("A/G".to_string()),
                            "dbsnp_ids" => Some(format!("rs{i}")),
                            "variation_name" => (i % 2 != 0).then(|| format!("rs{i}")),
                            "gnomADg" => Some("A:0.5".to_string()),
                            _ if af_members.contains(name) => Some(String::new()),
                            _ => None,
                        })
                        .collect();
                    Arc::new(StringArray::from(vals))
                }
                other => panic!("unexpected projected field type {other:?} for {name}"),
            };
            columns.push(arr);
        }
        let tier_batch = RecordBatch::try_new(Arc::clone(&projected), columns).unwrap();

        let encoded = encode_variation_batch(&tier_batch, &out_schema).unwrap();
        assert_eq!(encoded.schema(), out_schema);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_variation_parquet(tmp.path(), &[encoded]).unwrap();

        let lookup = SinglePathParquetVariationLookup::open(
            tmp.path(),
            vec!["gnomADg".to_string(), "variation_name".to_string()],
        )
        .await
        .unwrap();
        let mut cursor = lookup.new_cursor();
        let probes = vec![starts[0], starts[41], starts[299]];
        let taken = lookup.resolve_and_take(&probes, &mut cursor).await.unwrap();
        assert_eq!(taken.resolved.matched_positions, probes.len());

        let out = &taken.batch;
        let gnomadg = out
            .column_by_name("gnomADg")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let out_start = out
            .column_by_name("start")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let out_vn = out
            .column_by_name("variation_name")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        for r in 0..out.num_rows() {
            assert_eq!(gnomadg.value(r), "A:0.5");
            let orig = starts
                .iter()
                .position(|&x| x == out_start.value(r))
                .unwrap();
            // Even rows had a null variation_name → coalesced from dbsnp (`rs{orig}`).
            assert_eq!(out_vn.value(r), format!("rs{orig}"));
        }
    }
}
