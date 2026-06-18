use std::path::Path;
use std::time::Instant;

use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use lance::dataset::ProjectionRequest;

use crate::lance_cache::row_index::{
    PositionRowIdIndex, ResolvedRowIds, load_start_index_from_lance_btree,
};
use crate::lance_cache::schema::VARIATION_FORBIDDEN_COLUMNS;

#[derive(Debug)]
pub struct TakenVariationRows {
    pub resolved: ResolvedRowIds,
    pub batch: RecordBatch,
}

#[derive(Debug)]
pub struct SinglePathLanceVariationLookup {
    dataset: lance::Dataset,
    projection: Vec<String>,
    index: PositionRowIdIndex,
}

impl SinglePathLanceVariationLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self> {
        let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!(
                    "failed to open Lance variation dataset '{}': {err}",
                    path.display()
                ))
            })?;
        let index = load_start_index_from_lance_btree(&dataset).await?;
        Ok(Self {
            dataset,
            projection: ensure_runtime_projection(projection),
            index,
        })
    }

    pub async fn resolve_and_take(
        &self,
        sorted_unique_starts: &[u32],
        cursor: &mut usize,
    ) -> Result<TakenVariationRows> {
        let resolved = self
            .index
            .resolve_sorted_positions_from_cursor(sorted_unique_starts, cursor);
        let profile_enabled = std::env::var_os("VEP_LANCE_PROFILE").is_some();
        if profile_enabled {
            eprintln!(
                "[vep-lance-profile] variation_resolve requested_starts={} matched_positions={} row_ids={} cursor={}",
                resolved.requested_positions,
                resolved.matched_positions,
                resolved.row_ids.len(),
                *cursor,
            );
        }
        // The cache stores the 27 AF columns bundled into 3 fullzip List<Utf8> columns.
        // Project the bundled columns for the take, then expand them back so downstream
        // sees the 27 logical AF columns unchanged.
        let physical_projection =
            crate::lance_cache::af_bundle::bundle_projection(&self.projection);
        let projection_request = ProjectionRequest::from_columns(
            physical_projection.iter().map(String::as_str),
            self.dataset.schema(),
        );
        let take_started = profile_enabled.then(Instant::now);
        let batch = self
            .dataset
            .take_rows(&resolved.row_ids, projection_request)
            .await
            .map_err(|err| DataFusionError::Execution(format!("Lance take_rows failed: {err}")))?;
        let batch = crate::lance_cache::af_bundle::unbundle_af_columns(&batch)?;
        if let Some(started) = take_started {
            eprintln!(
                "[vep-lance-profile] variation_take row_ids={} batch_rows={} seconds={:.3}",
                resolved.row_ids.len(),
                batch.num_rows(),
                started.elapsed().as_secs_f64(),
            );
        }
        Ok(TakenVariationRows { resolved, batch })
    }

    pub fn projection(&self) -> &[String] {
        &self.projection
    }

    pub fn row_ids_len(&self) -> usize {
        self.index.row_ids_len()
    }

    pub fn unique_positions(&self) -> usize {
        self.index.unique_positions()
    }
}

pub fn ensure_runtime_projection(projection: Vec<String>) -> Vec<String> {
    let mut sanitized = Vec::with_capacity(projection.len() + 4);
    for column in projection {
        // `tier` is a build-only clustering column: it is stored in the dataset
        // (so it is no longer in VARIATION_FORBIDDEN_COLUMNS) but must never be
        // materialized into annotation output.
        let excluded = column == "tier"
            || VARIATION_FORBIDDEN_COLUMNS
                .iter()
                .any(|forbidden| column == *forbidden);
        if !excluded && !sanitized.iter().any(|existing| existing == &column) {
            sanitized.push(column);
        }
    }
    for required in ["start", "end", "allele_string", "failed"] {
        if !sanitized.iter().any(|column| column == required) {
            sanitized.push(required.to_string());
        }
    }
    sanitized
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    #[tokio::test]
    async fn lookup_resolves_start_rows_with_take_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.lance");
        let batch = test_variation_batch(vec![
            ("chr1", 10, 10, "A/T", 0, "rs10a"),
            ("chr1", 10, 10, "A/G", 0, "rs10b"),
            ("chr1", 20, 20, "C/G", 0, "rs20"),
        ]);
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![batch],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();
        let lookup = SinglePathLanceVariationLookup::open(
            &path,
            vec![
                "start".into(),
                "end".into(),
                "allele_string".into(),
                "failed".into(),
                "variation_name".into(),
            ],
        )
        .await
        .unwrap();

        let mut cursor = 0;
        let result = lookup
            .resolve_and_take(&[10, 20], &mut cursor)
            .await
            .unwrap();

        assert_eq!(result.batch.num_rows(), 3);
        assert_eq!(result.resolved.matched_positions, 2);
        assert_eq!(lookup.row_ids_len(), 3);
        assert_eq!(lookup.unique_positions(), 2);
    }

    #[test]
    fn runtime_projection_includes_required_lookup_columns_once() {
        let projection = ensure_runtime_projection(vec![
            "variation_name".into(),
            "start".into(),
            "allele_string".into(),
        ]);

        assert_eq!(
            projection,
            vec!["variation_name", "start", "allele_string", "end", "failed"]
        );
    }

    #[test]
    fn runtime_projection_drops_legacy_warm_cold_columns() {
        let projection = ensure_runtime_projection(vec![
            "position_key".into(),
            "variant_keys".into(),
            "tier".into(),
            "variation_name".into(),
            "position_key".into(),
        ]);

        assert_eq!(
            projection,
            vec!["variation_name", "start", "end", "allele_string", "failed"]
        );
    }

    #[tokio::test]
    async fn lance_lookup_round_trips_bundled_af_columns() {
        use crate::lance_cache::af_bundle::{af_column_order, bundle_af_columns};
        use datafusion::arrow::array::ArrayRef;

        let order = af_column_order();
        // 3 rows: all-present, gnomADg-only, all-absent.
        let af_for = |col: &str, row: usize| -> Option<String> {
            let is_gnomadg = col.starts_with("gnomADg");
            match row {
                0 => Some(format!("v_{col}")),
                1 => {
                    if is_gnomadg {
                        Some(format!("g_{col}"))
                    } else {
                        Some(String::new())
                    }
                }
                _ => Some(String::new()),
            }
        };
        let mut fields = vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ];
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(vec!["chr1", "chr1", "chr1"])),
            Arc::new(UInt32Array::from(vec![10u32, 20, 30])),
            Arc::new(UInt32Array::from(vec![10u32, 20, 30])),
            Arc::new(StringArray::from(vec!["A/T", "C/G", "G/A"])),
            Arc::new(Int8Array::from(vec![Some(0i8), Some(0), Some(0)])),
            Arc::new(StringArray::from(vec!["rs10", "rs20", "rs30"])),
        ];
        for col in &order {
            fields.push(Field::new(*col, DataType::Utf8, true));
            let vals: Vec<Option<String>> = (0..3).map(|r| af_for(col, r)).collect();
            cols.push(Arc::new(StringArray::from(vals)) as ArrayRef);
        }
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), cols).unwrap();
        let bundled = bundle_af_columns(&batch).unwrap();

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("var.lance");
        crate::lance_cache::write::write_record_batches_to_lance(
            &path,
            vec![bundled],
            crate::lance_cache::write::LanceIndexKind::Start,
        )
        .await
        .unwrap();

        // Projection includes the 27 AF column names (as the cold-parquet projection would).
        let mut projection: Vec<String> = vec![
            "variation_name".into(),
            "start".into(),
            "allele_string".into(),
        ];
        projection.extend(order.iter().map(|s| s.to_string()));
        let lookup = SinglePathLanceVariationLookup::open(&path, projection)
            .await
            .unwrap();
        let mut cursor = 0;
        let result = lookup
            .resolve_and_take(&[10, 20, 30], &mut cursor)
            .await
            .unwrap();
        let out = result.batch;
        assert_eq!(out.num_rows(), 3);

        // Every one of the 27 AF columns is present (unbundled) with exact values.
        for col in &order {
            let idx = out
                .schema()
                .index_of(col)
                .expect("AF column present after unbundle");
            let arr = out
                .column(idx)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            assert_eq!(arr.value(0), format!("v_{col}"), "row0 {col}");
            let expected1 = if col.starts_with("gnomADg") {
                format!("g_{col}")
            } else {
                String::new()
            };
            assert_eq!(arr.value(1), expected1, "row1 {col}");
            assert_eq!(arr.value(2), "", "row2 {col}");
        }
    }

    fn test_variation_batch(rows: Vec<(&str, u32, u32, &str, i8, &str)>) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int8, true),
            Field::new("variation_name", DataType::Utf8, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(
                    rows.iter().map(|row| row.0).collect::<Vec<_>>(),
                )),
                Arc::new(UInt32Array::from(
                    rows.iter().map(|row| row.1).collect::<Vec<_>>(),
                )),
                Arc::new(UInt32Array::from(
                    rows.iter().map(|row| row.2).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|row| row.3).collect::<Vec<_>>(),
                )),
                Arc::new(Int8Array::from(
                    rows.iter().map(|row| Some(row.4)).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|row| row.5).collect::<Vec<_>>(),
                )),
            ],
        )
        .unwrap()
    }
}
