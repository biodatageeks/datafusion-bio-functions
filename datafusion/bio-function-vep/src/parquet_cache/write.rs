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
}
