//! Parquet point-lookup reader for the variation cache, producing the same
//! [`TakenVariationRows`] contract as the Parquet reader.
//!
//! Three phases (the proven `pq_take` loop): (1) resolve candidate page
//! row-ranges from the footer [`PageDir`]; (2) a `start`-only read of those pages
//! to find the exact row offsets whose `start` is a requested position; (3) a
//! projected payload take at those offsets via the [`CoalescingAsyncReader`].
//!
//! The physical AF columns are the 2-array struct-of-arrays
//! (`<grp>_alleles`/`<grp>_freqs`); on read they are reconstructed to the 3
//! pipe-joined group strings and then expanded to the 27 logical AF columns via
//! [`unbundle_af_columns`] — byte-identical to the Parquet path. `variation_name`
//! is reconstructed via `coalesce(variation_name, dbsnp_ids)`. Binary flags stay
//! `Boolean` (downstream reads them through the `batch_i64_value` Boolean arm).

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use datafusion::arrow::array::{Array, ArrayRef, ListArray, RecordBatch, StringArray, UInt32Array};
use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};
use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;
use std::sync::Arc;

use crate::cache::af_bundle::{AF_GROUPS, unbundle_af_columns};
use crate::cache::variation_runtime::{
    ResolvedRowIds, TakenVariationRows, ensure_runtime_projection,
};
use crate::parquet_cache::encode::reconstruct_af_group_string;
use crate::parquet_cache::page_dir::{
    CoalescingAsyncReader, IoCounters, PageDir, selection_from_offsets, selection_from_ranges,
};

/// Coalescing-reader window (bytes). 64 KiB is the measured sweet spot.
const COALESCE_GAP_BYTES: u64 = 64 * 1024;

/// Per-contig Parquet variation lookup.
pub struct SinglePathParquetVariationLookup {
    path: PathBuf,
    meta: ArrowReaderMetadata,
    page_dir: PageDir,
    start_leaf: usize,
    /// Sanitized logical projection (post [`ensure_runtime_projection`]).
    projection: Vec<String>,
}

impl SinglePathParquetVariationLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let file = std::fs::File::open(path).map_err(|e| {
            DataFusionError::Execution(format!("open parquet variation '{}': {e}", path.display()))
        })?;
        // `with_page_index_policy` is the non-deprecated API, but `PageIndexPolicy`
        // is private in parquet 58, so the deprecated setter is the only usable one.
        #[allow(deprecated)]
        let meta =
            ArrowReaderMetadata::load(&file, ArrowReaderOptions::new().with_page_index(true))
                .map_err(|e| DataFusionError::Execution(format!("load parquet metadata: {e}")))?;
        let start_leaf = meta
            .parquet_schema()
            .columns()
            .iter()
            .position(|c| c.name() == "start")
            .ok_or_else(|| {
                DataFusionError::Execution("variation parquet has no 'start' column".to_string())
            })?;
        let page_dir = PageDir::build(&meta, start_leaf)?;
        Ok(Self {
            path: path.to_path_buf(),
            meta,
            page_dir,
            start_leaf,
            projection: ensure_runtime_projection(projection),
        })
    }

    pub fn projection(&self) -> &[String] {
        &self.projection
    }

    /// Physical top-level column names to read for the payload take, derived from
    /// the logical projection: logical AF columns map to their group's 2-array
    /// pair, `variation_name` also pulls `dbsnp_ids`, everything else is physical.
    /// `start` is always included.
    fn physical_columns(&self) -> Vec<String> {
        let arrow_schema = self.meta.schema();
        let mut names: Vec<String> = Vec::new();
        let push = |n: &str, names: &mut Vec<String>| {
            if arrow_schema.index_of(n).is_ok() && !names.iter().any(|x| x == n) {
                names.push(n.to_string());
            }
        };
        for col in &self.projection {
            if let Some((grp, _)) = AF_GROUPS
                .iter()
                .find(|(_, members)| members.contains(&col.as_str()))
            {
                push(&format!("{grp}_alleles"), &mut names);
                push(&format!("{grp}_freqs"), &mut names);
            } else if col == "variation_name" {
                push("variation_name", &mut names);
                push("dbsnp_ids", &mut names);
            } else {
                push(col, &mut names);
            }
        }
        push("start", &mut names);
        names
    }

    pub async fn resolve_and_take(
        &self,
        sorted_unique_starts: &[u32],
    ) -> Result<TakenVariationRows> {
        let probe_set: HashSet<u32> = sorted_unique_starts.iter().copied().collect();
        let counters = IoCounters::new();

        // Phase 1: candidate page row-ranges from the footer PageDir. The PageDir
        // holds keys as u64; variation's `start` probes widen losslessly.
        let probes64: Vec<u64> = sorted_unique_starts.iter().map(|&s| s as u64).collect();
        let ranges = self.page_dir.resolve_ranges(&probes64);

        // Phase 2: start-only read of the candidate pages -> exact row offsets.
        let offsets = self.exact_offsets(&ranges, &probe_set, &counters).await?;

        // Phase 3: projected payload take at the exact offsets.
        let phys = self.take_payload(&offsets, &counters).await?;

        // Post-process to the Parquet-equivalent logical batch.
        let batch = self.to_logical_batch(&phys)?;

        // matched positions = distinct requested starts present in the result.
        let matched = distinct_matched(&batch, &probe_set)?;
        let resolved = ResolvedRowIds {
            requested_positions: sorted_unique_starts.len(),
            matched_positions: matched,
            row_ids: offsets,
        };
        Ok(TakenVariationRows { resolved, batch })
    }

    async fn exact_offsets(
        &self,
        ranges: &[(u64, u64)],
        probe_set: &HashSet<u32>,
        counters: &IoCounters,
    ) -> Result<Vec<u64>> {
        if ranges.is_empty() {
            return Ok(Vec::new());
        }
        let pq_schema = self.meta.parquet_schema();
        let start_mask = ProjectionMask::leaves(pq_schema, [self.start_leaf]);
        let reader = CoalescingAsyncReader::new(
            self.open_async().await?,
            counters.clone(),
            COALESCE_GAP_BYTES,
        );
        let mut stream =
            ParquetRecordBatchStreamBuilder::new_with_metadata(reader, self.meta.clone())
                .with_projection(start_mask)
                .with_row_selection(selection_from_ranges(ranges))
                .with_batch_size(8192)
                .build()
                .map_err(|e| DataFusionError::Execution(format!("build start stream: {e}")))?;
        let mut offsets = Vec::new();
        let mut cursor = ranges.iter().flat_map(|&(a, b)| a..b);
        while let Some(b) = stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read start batch: {e}")))?
        {
            let sa = b
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| DataFusionError::Execution("start column must be UInt32".into()))?;
            for &v in sa.values() {
                let off = cursor.next().ok_or_else(|| {
                    DataFusionError::Execution("row range cursor underflow".into())
                })?;
                if probe_set.contains(&v) {
                    offsets.push(off);
                }
            }
        }
        Ok(offsets)
    }

    async fn take_payload(&self, offsets: &[u64], counters: &IoCounters) -> Result<RecordBatch> {
        let pq_schema = self.meta.parquet_schema();
        let arrow_schema = self.meta.schema();
        let root_indices: Vec<usize> = self
            .physical_columns()
            .iter()
            .filter_map(|n| arrow_schema.index_of(n).ok())
            .collect();
        let mask = ProjectionMask::roots(pq_schema, root_indices);
        let reader = CoalescingAsyncReader::new(
            self.open_async().await?,
            counters.clone(),
            COALESCE_GAP_BYTES,
        );
        let builder = ParquetRecordBatchStreamBuilder::new_with_metadata(reader, self.meta.clone())
            .with_projection(mask)
            .with_row_selection(selection_from_offsets(offsets))
            .with_batch_size(8192);
        let proj_schema = builder.schema().clone();
        let mut stream = builder
            .build()
            .map_err(|e| DataFusionError::Execution(format!("build payload stream: {e}")))?;
        let mut taken: Vec<RecordBatch> = Vec::new();
        while let Some(b) = stream
            .try_next()
            .await
            .map_err(|e| DataFusionError::Execution(format!("read payload batch: {e}")))?
        {
            taken.push(b);
        }
        if taken.is_empty() {
            Ok(RecordBatch::new_empty(proj_schema))
        } else {
            concat_batches(&taken[0].schema(), &taken)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
        }
    }

    async fn open_async(&self) -> Result<tokio::fs::File> {
        tokio::fs::File::open(&self.path).await.map_err(|e| {
            DataFusionError::Execution(format!("reopen parquet '{}': {e}", self.path.display()))
        })
    }

    /// Reconstruct the group AF strings + unbundle to 27 logical columns, and
    /// coalesce `variation_name`. Non-AF columns pass through unchanged.
    fn to_logical_batch(&self, phys: &RecordBatch) -> Result<RecordBatch> {
        let schema = phys.schema();
        let dbsnp = phys
            .column_by_name("dbsnp_ids")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>());

        let mut fields: Vec<Arc<Field>> = Vec::new();
        let mut cols: Vec<ArrayRef> = Vec::new();

        for (i, f) in schema.fields().iter().enumerate() {
            let name = f.name();
            // The 2-array AF columns are replaced by reconstructed group strings.
            if name.ends_with("_alleles") || name.ends_with("_freqs") {
                let base = name.trim_end_matches("_alleles").trim_end_matches("_freqs");
                if AF_GROUPS.iter().any(|(g, _)| *g == base) {
                    continue;
                }
            }
            if name == "variation_name" {
                let vn = phys
                    .column(i)
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .ok_or_else(|| {
                        DataFusionError::Execution("variation_name must be Utf8".into())
                    })?;
                fields.push(Arc::new(Field::new("variation_name", DataType::Utf8, true)));
                cols.push(Arc::new(coalesce_variation_name(vn, dbsnp)) as ArrayRef);
            } else {
                fields.push(f.clone());
                cols.push(phys.column(i).clone());
            }
        }

        // Append reconstructed group string columns for each present group.
        for (grp, members) in AF_GROUPS {
            let alleles = phys
                .column_by_name(&format!("{grp}_alleles"))
                .and_then(|c| c.as_any().downcast_ref::<ListArray>());
            let freqs = phys
                .column_by_name(&format!("{grp}_freqs"))
                .and_then(|c| c.as_any().downcast_ref::<ListArray>());
            if let (Some(alleles), Some(freqs)) = (alleles, freqs) {
                let s = reconstruct_af_group_string(alleles, freqs, members.len())?;
                fields.push(Arc::new(Field::new(*grp, DataType::Utf8, true)));
                cols.push(Arc::new(s) as ArrayRef);
            }
        }

        let bundled = RecordBatch::try_new(
            Arc::new(Schema::new_with_metadata(fields, schema.metadata().clone())),
            cols,
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        unbundle_af_columns(&bundled)
    }
}

/// `coalesce(vn, dbsnp)`: keep `vn` where present, else fall back to `dbsnp`.
fn coalesce_variation_name(vn: &StringArray, dbsnp: Option<&StringArray>) -> StringArray {
    match dbsnp {
        None => vn.clone(),
        Some(db) => (0..vn.len())
            .map(|i| {
                if !vn.is_null(i) {
                    Some(vn.value(i))
                } else if i < db.len() && !db.is_null(i) {
                    Some(db.value(i))
                } else {
                    None
                }
            })
            .collect(),
    }
}

fn distinct_matched(batch: &RecordBatch, probe_set: &HashSet<u32>) -> Result<usize> {
    let start = batch
        .column_by_name("start")
        .and_then(|c| c.as_any().downcast_ref::<UInt32Array>())
        .ok_or_else(|| DataFusionError::Execution("logical batch missing UInt32 start".into()))?;
    let mut seen: HashSet<u32> = HashSet::new();
    for &v in start.values() {
        if probe_set.contains(&v) {
            seen.insert(v);
        }
    }
    Ok(seen.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parquet_cache::encode::{encode_af_2array, presence_boolean};
    use crate::parquet_cache::write::write_variation_parquet;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    #[tokio::test(flavor = "multi_thread")]
    async fn parquet_lookup_reconstructs_logical_af_and_vn() {
        const N_POPS: usize = 11; // gnomADg group
        const AF: &str =
            "A:0.08806|A:0.0367|A:2.682e-05||A:0.9969|A:0.5|A:0.006912||A:9.911e-05|A:0.4253|A:1";
        let n = 500usize;
        let starts: Vec<u32> = (0..n).map(|i| (i as u32) * 10 + 5).collect();
        let start_arr = UInt32Array::from(starts.clone());
        let allele = StringArray::from((0..n).map(|_| "A/G").collect::<Vec<_>>());
        let failed = presence_boolean(
            &Int8Array::from(
                (0..n)
                    .map(|i| if i % 5 == 0 { Some(1) } else { None })
                    .collect::<Vec<_>>(),
            ),
            "failed",
        )
        .unwrap();
        // variation_name: null on even rows (should coalesce to dbsnp), set on odd.
        let vn = StringArray::from(
            (0..n)
                .map(|i| {
                    if i % 2 == 0 {
                        None
                    } else {
                        Some(format!("rs{i}"))
                    }
                })
                .collect::<Vec<_>>(),
        );
        let dbsnp = StringArray::from((0..n).map(|i| Some(format!("rs{i}"))).collect::<Vec<_>>());
        let af = encode_af_2array(
            &StringArray::from((0..n).map(|_| AF).collect::<Vec<_>>()),
            N_POPS,
        )
        .unwrap();

        let schema = Arc::new(Schema::new(vec![
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Boolean, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("dbsnp_ids", DataType::Utf8, true),
            Field::new("af_gnomadg_alleles", af.alleles.data_type().clone(), true),
            Field::new("af_gnomadg_freqs", af.freqs.data_type().clone(), true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(start_arr),
                Arc::new(allele),
                Arc::new(failed),
                Arc::new(vn),
                Arc::new(dbsnp),
                Arc::new(af.alleles),
                Arc::new(af.freqs),
            ],
        )
        .unwrap();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_variation_parquet(tmp.path(), &[batch]).unwrap();

        let lookup = SinglePathParquetVariationLookup::open(
            tmp.path(),
            vec![
                "gnomADg".to_string(),
                "gnomADg_NFE".to_string(),
                "variation_name".to_string(),
            ],
        )
        .await
        .unwrap();
        let probes = vec![starts[4], starts[101], starts[400]]; // rows 4,101,400
        let taken = lookup.resolve_and_take(&probes).await.unwrap();

        assert_eq!(taken.resolved.matched_positions, probes.len());
        let out = &taken.batch;
        // The 27 logical AF columns must be present (via unbundle).
        let gnomadg = out
            .column_by_name("gnomADg")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let gnomadg_nfe = out
            .column_by_name("gnomADg_NFE")
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
        // AF string groups: gnomADg = population 0 = "A:0.08806"; gnomADg_NFE is
        // population index 8 (per AF_GROUPS gnomADg order) = "A:9.911e-05".
        for r in 0..out.num_rows() {
            assert_eq!(gnomadg.value(r), "A:0.08806");
            assert_eq!(gnomadg_nfe.value(r), "A:9.911e-05");
            // variation_name coalesced: every row has a value (even rows from dbsnp).
            let orig = starts
                .iter()
                .position(|&x| x == out_start.value(r))
                .unwrap();
            assert_eq!(out_vn.value(r), format!("rs{orig}"));
        }
    }

    /// Build a small shard whose three AF groups cover the shapes the read path
    /// has to reproduce byte for byte, and return it with the source group text.
    /// `None` is an absent group; every `Some` has exactly one `|`-separated
    /// segment per population, so segment `p` IS the expected member value.
    #[allow(clippy::type_complexity)]
    fn af_shape_fixture() -> (
        tempfile::NamedTempFile,
        Vec<u32>,
        Vec<(&'static str, Vec<Option<&'static str>>)>,
    ) {
        let groups: Vec<(&'static str, Vec<Option<&'static str>>)> = vec![
            (
                "af_global",
                vec![
                    Some("A:0.1|A:0.2|A:0.3|A:0.4|A:0.5|A:0.6"),
                    // multi-allelic; population 1 lacks G, population 2 is absent
                    Some("A:0.1,G:0.2|A:0.3||A:0.5,G:0.6|A:0,G:1|A:0.25,G:0.75"),
                    None,
                    // scientific notation, then five absent populations
                    Some("T:2.682e-05|||||"),
                    None,
                    Some("C:0.5|C:0.5|C:0.5|C:0.5|C:0.5|C:0.5"),
                ],
            ),
            (
                "af_gnomade",
                vec![
                    None,
                    Some(
                        "A:9.911e-05|A:0.4253|A:0.006912|A:0.9969|A:0.0367|A:0.08806|A:4.248e-06|A:5.896e-05|A:1|A:0",
                    ),
                    None,
                    // deletion allele and a long insertion allele
                    Some("-:9.47e-05,AAAAAAAAAA:0.0003568|-:0.5,AAAAAAAAAA:0.25||||||||"),
                    None,
                    None,
                ],
            ),
            (
                "af_gnomadg",
                vec![
                    Some(
                        "A:0.08806|A:0.0367|A:2.682e-05||A:0.9969|A:0.5|A:0.006912||A:9.911e-05|A:0.4253|A:1",
                    ),
                    Some("C:0.0003994,G:0.0001997|C:0.001,G:0.002|||||||||C:1.47e-05,G:2.412e-05"),
                    None,
                    Some("||||||||||T:0.01"),
                    Some("G:0|G:0|G:0|G:0|G:0|G:0|G:0|G:0|G:0|G:0|G:0"),
                    None,
                ],
            ),
        ];
        let n = 6usize;
        let starts: Vec<u32> = (0..n).map(|i| (i as u32) * 10 + 5).collect();
        let mut fields = vec![
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Boolean, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("dbsnp_ids", DataType::Utf8, true),
        ];
        let failed = presence_boolean(
            &Int8Array::from((0..n).map(|_| None::<i8>).collect::<Vec<_>>()),
            "failed",
        )
        .unwrap();
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(UInt32Array::from(starts.clone())),
            Arc::new(StringArray::from(vec!["A/G"; n])),
            Arc::new(failed),
            Arc::new(StringArray::from(vec![None::<&str>; n])),
            Arc::new(StringArray::from(
                (0..n).map(|i| Some(format!("rs{i}"))).collect::<Vec<_>>(),
            )),
        ];
        for (grp, rows) in &groups {
            let width = AF_GROUPS.iter().find(|(g, _)| g == grp).unwrap().1.len();
            let af = encode_af_2array(&StringArray::from(rows.clone()), width).unwrap();
            fields.push(Field::new(
                format!("{grp}_alleles"),
                af.alleles.data_type().clone(),
                true,
            ));
            fields.push(Field::new(
                format!("{grp}_freqs"),
                af.freqs.data_type().clone(),
                true,
            ));
            cols.push(Arc::new(af.alleles));
            cols.push(Arc::new(af.freqs));
        }
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), cols).unwrap();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_variation_parquet(tmp.path(), &[batch]).unwrap();
        (tmp, starts, groups)
    }

    /// Pins the whole logical batch, not one or two columns of it: the 27 AF
    /// members must come out in `AF_GROUPS` order after the non-AF columns, as
    /// nullable `Utf8` with NO nulls (an absent group or population is `""` --
    /// a NULL here would surface as NULL in `lookup_variants()` output), and
    /// every cell must equal its segment of the source text exactly.
    #[tokio::test(flavor = "multi_thread")]
    async fn logical_batch_pins_all_27_af_members_byte_for_byte() {
        let (tmp, starts, groups) = af_shape_fixture();
        let mut projection: Vec<String> = crate::cache::af_bundle::af_column_order()
            .into_iter()
            .map(str::to_string)
            .collect();
        projection.extend(["variation_name", "allele_string", "failed"].map(str::to_string));
        let lookup = SinglePathParquetVariationLookup::open(tmp.path(), projection)
            .await
            .unwrap();
        let taken = lookup.resolve_and_take(&starts).await.unwrap();
        let out = &taken.batch;
        assert_eq!(out.num_rows(), starts.len());

        let names: Vec<String> = out
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();
        let af_order = crate::cache::af_bundle::af_column_order();
        let (non_af, af) = names.split_at(names.len() - af_order.len());
        assert_eq!(
            non_af,
            [
                "start",
                "allele_string",
                "failed",
                "variation_name",
                "dbsnp_ids"
            ],
            "non-AF columns keep the file's order"
        );
        assert_eq!(
            af,
            af_order.as_slice(),
            "AF members in AF_GROUPS order, last"
        );
        assert_eq!(out.schema().metadata(), lookup.meta.schema().metadata());

        for (grp, rows) in &groups {
            let members = AF_GROUPS.iter().find(|(g, _)| g == grp).unwrap().1;
            for (p, member) in members.iter().enumerate() {
                let field = out.schema().field_with_name(member).unwrap().clone();
                assert_eq!(field.data_type(), &DataType::Utf8, "{member}");
                assert!(field.is_nullable(), "{member}");
                assert!(field.metadata().is_empty(), "{member}");
                let col = out
                    .column_by_name(member)
                    .unwrap()
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .unwrap();
                assert_eq!(col.null_count(), 0, "{member} must hold \"\", never NULL");
                for (r, src) in rows.iter().enumerate() {
                    let expected = src.map_or("", |s| s.split('|').nth(p).unwrap());
                    assert_eq!(col.value(r), expected, "{member} row {r}");
                }
            }
        }
    }

    /// Asking for one member still materialises its whole group and nothing
    /// else: `physical_columns` pulls the group's list pair, and downstream
    /// resolves the members by name.
    #[tokio::test(flavor = "multi_thread")]
    async fn logical_batch_emits_the_full_group_for_a_single_projected_member() {
        let (tmp, starts, _) = af_shape_fixture();
        let lookup =
            SinglePathParquetVariationLookup::open(tmp.path(), vec!["gnomADg_NFE".to_string()])
                .await
                .unwrap();
        let taken = lookup.resolve_and_take(&starts).await.unwrap();
        let schema = taken.batch.schema();
        let af_order = crate::cache::af_bundle::af_column_order();
        let af: Vec<&str> = schema
            .fields()
            .iter()
            .map(|f| f.name().as_str())
            .filter(|n| af_order.contains(n))
            .collect();
        let gnomadg = AF_GROUPS
            .iter()
            .find(|(g, _)| *g == "af_gnomadg")
            .unwrap()
            .1;
        assert_eq!(af, gnomadg);
    }

    /// Probes that hit nothing still return every AF member, in order. (The
    /// non-AF part of an empty batch is the file's unprojected schema today --
    /// `take_payload` builds it from the reader builder -- so only the AF tail is
    /// pinned here; consumers resolve columns by name.)
    #[tokio::test(flavor = "multi_thread")]
    async fn logical_batch_keeps_the_af_members_when_nothing_matches() {
        let (tmp, starts, _) = af_shape_fixture();
        let af_order = crate::cache::af_bundle::af_column_order();
        let projection: Vec<String> = af_order.iter().map(|c| c.to_string()).collect();
        let lookup = SinglePathParquetVariationLookup::open(tmp.path(), projection)
            .await
            .unwrap();
        let miss = lookup.resolve_and_take(&[starts[0] + 1]).await.unwrap();
        assert_eq!(miss.batch.num_rows(), 0);
        assert_eq!(miss.resolved.matched_positions, 0);
        let schema = miss.batch.schema();
        let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert_eq!(&names[names.len() - af_order.len()..], af_order.as_slice());
        for member in &af_order {
            assert_eq!(
                schema.field_with_name(member).unwrap().data_type(),
                &DataType::Utf8
            );
        }
    }
}
