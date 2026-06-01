use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use datafusion::arrow::compute::concat_batches;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReaderBuilder, RowSelection,
};
use parquet::file::page_index::index::Index;
use parquet::file::statistics::Statistics;

use crate::warm_cache::chunk::WarmChunkContext;
use crate::warm_cache::reader::projection_for_existing_roots;

const DEFAULT_COLD_PARQUET_BATCH_SIZE: usize = 262_144;
const DEFAULT_COLD_PARQUET_ROW_GROUP_CACHE: usize = 2;

#[derive(Debug)]
pub struct ColdParquetLookup {
    path: PathBuf,
    metadata: ArrowReaderMetadata,
    projection_columns: Vec<String>,
    batch_size: usize,
    max_cached_row_groups: usize,
    position_leaf: usize,
    page_layout: ColdPageLayoutStats,
    cursor: ColdRowGroupCursor,
    cache: VecDeque<CachedColdRowGroup>,
    row_groups_touched: Vec<bool>,
    row_groups_unique_touched: u64,
    pub row_groups_loaded: u64,
    pub row_group_cache_hits: u64,
    pub row_group_cache_misses: u64,
    pub rows_loaded: u64,
    pub load_time: Duration,
    page_index_probes: u64,
    page_index_available_probes: u64,
    page_index_unavailable_probes: u64,
    page_index_pages_in_probed_row_groups: u64,
    page_index_candidate_pages: u64,
    page_index_candidate_rows: u64,
    page_index_candidate_misses: u64,
    page_index_unique_candidate_pages: HashSet<(usize, usize)>,
    page_index_unique_candidate_rows: u64,
}

#[derive(Debug)]
pub struct ColdParquetLookupSet {
    lookups: Vec<ColdParquetLookup>,
    cursor: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ColdParquetPositionGroup {
    lookup_idx: usize,
    row_group_id: usize,
}

#[derive(Debug)]
struct CachedColdRowGroup {
    row_group_id: usize,
    row_ranges: Option<Vec<Range<usize>>>,
    chunk: WarmChunkContext,
}

#[derive(Debug, Default)]
struct PageProbeStats {
    probes: u64,
    available_probes: u64,
    unavailable_probes: u64,
    pages_in_probed_row_groups: u64,
    candidate_pages: u64,
    candidate_rows: u64,
    candidate_misses: u64,
    unique_candidate_pages: Vec<(usize, usize, usize)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColdProbeResult {
    Match,
    PositionCoveredNoExact,
    NotCovered,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColdRowGroupMeta {
    pub rows: usize,
    pub min_position_key: u64,
    pub max_position_key: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ColdRowGroupCursorStats {
    pub metadata_probes: u64,
    pub current_hits: u64,
    pub previous_hits: u64,
    pub advanced_hits: u64,
    pub binary_search_hits: u64,
    pub misses: u64,
    pub skipped_ahead: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ColdPageLayoutStats {
    pub position_page_index_loaded: bool,
    pub position_column_index_loaded: bool,
    pub position_pages_total: u64,
    pub position_bloom_filter_row_groups: u64,
    pub position_bloom_filter_bytes: u64,
}

#[derive(Debug, Clone)]
pub struct ColdRowGroupCursor {
    row_groups: Vec<ColdRowGroupMeta>,
    cursor: usize,
    stats: ColdRowGroupCursorStats,
}

impl ColdParquetLookup {
    pub fn open(
        path: impl AsRef<Path>,
        projection_columns: Vec<String>,
        batch_size: usize,
        max_cached_row_groups: usize,
    ) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path).map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to open cold variation parquet '{}': {error}",
                path.display()
            ))
        })?;
        let metadata_options =
            ArrowReaderOptions::new().with_page_index(cold_parquet_load_page_index());
        let metadata = ArrowReaderMetadata::load(&file, metadata_options)?;
        reject_deprecated_cold_variant_keys(&metadata, &path)?;
        let position_leaf = cold_position_leaf_index(&metadata)?;
        let row_groups = cold_row_group_metadata(&metadata, position_leaf)?;
        validate_cold_row_group_order(&row_groups)?;
        let page_layout = cold_page_layout_stats(&metadata, position_leaf);
        let row_group_count = metadata.metadata().num_row_groups();

        Ok(Self {
            path,
            metadata,
            projection_columns,
            batch_size: batch_size.max(1),
            max_cached_row_groups: max_cached_row_groups.max(1),
            position_leaf,
            page_layout,
            cursor: ColdRowGroupCursor::new(row_groups),
            cache: VecDeque::new(),
            row_groups_touched: vec![false; row_group_count],
            row_groups_unique_touched: 0,
            row_groups_loaded: 0,
            row_group_cache_hits: 0,
            row_group_cache_misses: 0,
            rows_loaded: 0,
            load_time: Duration::ZERO,
            page_index_probes: 0,
            page_index_available_probes: 0,
            page_index_unavailable_probes: 0,
            page_index_pages_in_probed_row_groups: 0,
            page_index_candidate_pages: 0,
            page_index_candidate_rows: 0,
            page_index_candidate_misses: 0,
            page_index_unique_candidate_pages: HashSet::new(),
            page_index_unique_candidate_rows: 0,
        })
    }

    pub fn from_env(path: impl AsRef<Path>, projection_columns: Vec<String>) -> Result<Self> {
        Self::open(
            path,
            projection_columns,
            cold_parquet_batch_size(),
            cold_parquet_row_group_cache(),
        )
    }

    fn position_range(&self) -> Option<(u64, u64)> {
        let first = self.cursor.row_groups.first()?;
        let last = self.cursor.row_groups.last()?;
        Some((first.min_position_key, last.max_position_key))
    }

    fn contains_position(&self, position_key: u64) -> bool {
        self.position_range()
            .is_some_and(|(min, max)| position_key >= min && position_key <= max)
    }

    pub fn probe_position_emit_and_visit<P, E, V>(
        &mut self,
        position_key: u64,
        mut allele_matches: P,
        mut emit_match: E,
        mut visit_row: V,
    ) -> Result<ColdProbeResult>
    where
        P: FnMut(&WarmChunkContext, u32, &str) -> Result<bool>,
        E: FnMut(&WarmChunkContext, u32) -> Result<()>,
        V: FnMut(&WarmChunkContext, u32, &str) -> Result<()>,
    {
        let Some(row_group_id) = self.cursor.find_row_group(position_key) else {
            return Ok(ColdProbeResult::NotCovered);
        };
        let row_ranges = self.candidate_page_row_ranges(row_group_id, position_key);
        if row_ranges.as_ref().is_some_and(Vec::is_empty) {
            return Ok(ColdProbeResult::NotCovered);
        }
        let cache_idx = self.ensure_row_group(row_group_id, row_ranges)?;
        let chunk = &self.cache[cache_idx].chunk;
        let rows = chunk.rows_for_position(position_key);
        if rows.is_empty() {
            return Ok(ColdProbeResult::NotCovered);
        }

        let mut matched = false;
        for row in rows {
            let Some(allele_string) = chunk.allele_string(row as usize)? else {
                continue;
            };
            visit_row(chunk, row, allele_string)?;
            if !matched && allele_matches(chunk, row, allele_string)? {
                emit_match(chunk, row)?;
                matched = true;
            }
        }

        if matched {
            Ok(ColdProbeResult::Match)
        } else {
            Ok(ColdProbeResult::PositionCoveredNoExact)
        }
    }

    pub fn prefetch_positions<I>(&mut self, position_keys: I) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        self.prefetch_positions_inner(position_keys, false)
    }

    pub fn prefetch_positions_retaining<I>(&mut self, position_keys: I) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        self.prefetch_positions_inner(position_keys, true)
    }

    fn prefetch_positions_inner<I>(
        &mut self,
        position_keys: I,
        retain_prefetched: bool,
    ) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        let mut full_row_groups = HashSet::new();
        let mut ranges_by_row_group: HashMap<usize, Vec<Range<usize>>> = HashMap::new();
        let mut probe_stats = PageProbeStats::default();
        for position_key in position_keys {
            let Some(row_group_id) = self.cursor.find_row_group(position_key) else {
                continue;
            };
            let Some(row_ranges) = self.candidate_page_row_ranges_for_probe(
                row_group_id,
                position_key,
                &mut probe_stats,
            ) else {
                full_row_groups.insert(row_group_id);
                continue;
            };
            if row_ranges.is_empty() {
                continue;
            }
            ranges_by_row_group
                .entry(row_group_id)
                .or_default()
                .extend(row_ranges);
        }

        let ranged_row_groups = ranges_by_row_group
            .keys()
            .filter(|row_group_id| !full_row_groups.contains(row_group_id))
            .count();
        let old_max_cached_row_groups = self.max_cached_row_groups;
        if retain_prefetched {
            let required_capacity = self
                .cache
                .len()
                .saturating_add(full_row_groups.len())
                .saturating_add(ranged_row_groups);
            self.max_cached_row_groups = self.max_cached_row_groups.max(required_capacity);
        }

        let result = (|| {
            let mut full_row_group_ids = full_row_groups.iter().copied().collect::<Vec<_>>();
            full_row_group_ids.sort_unstable();
            for row_group_id in full_row_group_ids {
                self.ensure_row_group(row_group_id, None)?;
            }

            let mut ranged_row_groups = ranges_by_row_group.into_iter().collect::<Vec<_>>();
            ranged_row_groups.sort_by_key(|(row_group_id, _)| *row_group_id);
            for (row_group_id, mut row_ranges) in ranged_row_groups {
                if full_row_groups.contains(&row_group_id) {
                    continue;
                }
                merge_row_ranges(&mut row_ranges);
                self.ensure_row_group(row_group_id, Some(row_ranges))?;
            }

            self.record_page_probe_stats(probe_stats);
            Ok(())
        })();

        self.max_cached_row_groups = old_max_cached_row_groups;
        result
    }

    pub fn trim_cache_to_capacity(&mut self) {
        while self.cache.len() > self.max_cached_row_groups {
            self.cache.pop_front();
        }
    }

    pub fn cached_chunk_for_position(&self, position_key: u64) -> Option<&WarmChunkContext> {
        let row_group_id = self.find_row_group_without_stats(position_key)?;
        self.cache
            .iter()
            .rev()
            .find(|cached| cached.row_group_id == row_group_id)
            .map(|cached| &cached.chunk)
    }

    fn ensure_row_group(
        &mut self,
        row_group_id: usize,
        row_ranges: Option<Vec<Range<usize>>>,
    ) -> Result<usize> {
        if let Some(idx) = self.cache.iter().position(|cached| {
            cached.row_group_id == row_group_id
                && cached_row_ranges_cover(cached.row_ranges.as_deref(), row_ranges.as_deref())
        }) {
            self.row_group_cache_hits += 1;
            return Ok(idx);
        }

        self.row_group_cache_misses += 1;
        let started = Instant::now();
        let chunk = self.load_row_group(row_group_id, row_ranges.as_deref())?;
        self.load_time += started.elapsed();
        self.rows_loaded += chunk.batch.num_rows() as u64;
        self.row_groups_loaded += 1;
        if !self.row_groups_touched[row_group_id] {
            self.row_groups_touched[row_group_id] = true;
            self.row_groups_unique_touched += 1;
        }

        self.cache.push_back(CachedColdRowGroup {
            row_group_id,
            row_ranges,
            chunk,
        });
        while self.cache.len() > self.max_cached_row_groups {
            self.cache.pop_front();
        }
        Ok(self.cache.len() - 1)
    }

    fn load_row_group(
        &self,
        row_group_id: usize,
        row_ranges: Option<&[Range<usize>]>,
    ) -> Result<WarmChunkContext> {
        let mask = projection_for_existing_roots(
            self.metadata.schema(),
            self.metadata.parquet_schema(),
            &self.projection_columns,
        );
        let file = File::open(&self.path).map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to open cold variation parquet '{}': {error}",
                self.path.display()
            ))
        })?;
        let mut builder =
            ParquetRecordBatchReaderBuilder::new_with_metadata(file, self.metadata.clone())
                .with_projection(mask)
                .with_row_groups(vec![row_group_id])
                .with_batch_size(self.batch_size);
        if let Some(row_ranges) = row_ranges
            && !row_ranges.is_empty()
        {
            let row_group_rows = self.cursor.row_groups[row_group_id].rows;
            let selection =
                RowSelection::from_consecutive_ranges(row_ranges.iter().cloned(), row_group_rows);
            builder = builder.with_row_selection(selection);
        }
        let reader = builder.build()?;
        let batches = reader.collect::<std::result::Result<Vec<_>, _>>()?;
        let batch = match batches.as_slice() {
            [] => {
                return Err(DataFusionError::Execution(format!(
                    "cold row group {row_group_id} produced no batches"
                )));
            }
            [single] => single.clone(),
            _ => concat_batches(&batches[0].schema(), batches.iter())
                .map_err(|error| DataFusionError::ArrowError(Box::new(error), None))?,
        };

        WarmChunkContext::try_new_without_variant_index(row_group_id, batch)
    }

    fn candidate_page_row_ranges(
        &mut self,
        row_group_id: usize,
        position_key: u64,
    ) -> Option<Vec<Range<usize>>> {
        let mut probe_stats = PageProbeStats::default();
        let row_ranges =
            self.candidate_page_row_ranges_for_probe(row_group_id, position_key, &mut probe_stats);
        self.record_page_probe_stats(probe_stats);
        row_ranges
    }

    fn candidate_page_row_ranges_without_recording(
        &self,
        row_group_id: usize,
        position_key: u64,
    ) -> Option<Vec<Range<usize>>> {
        let mut probe_stats = PageProbeStats::default();
        self.candidate_page_row_ranges_for_probe(row_group_id, position_key, &mut probe_stats)
    }

    fn candidate_page_row_ranges_for_probe(
        &self,
        row_group_id: usize,
        position_key: u64,
        stats: &mut PageProbeStats,
    ) -> Option<Vec<Range<usize>>> {
        stats.probes += 1;
        let metadata = self.metadata.metadata();
        let Some(column_index) = metadata.column_index() else {
            stats.unavailable_probes += 1;
            return None;
        };
        let Some(offset_index) = metadata.offset_index() else {
            stats.unavailable_probes += 1;
            return None;
        };
        let Some(row_group_offsets) = offset_index.get(row_group_id) else {
            stats.unavailable_probes += 1;
            return None;
        };
        let Some(position_offsets) = row_group_offsets.get(self.position_leaf) else {
            stats.unavailable_probes += 1;
            return None;
        };
        let Some(row_group_indexes) = column_index.get(row_group_id) else {
            stats.unavailable_probes += 1;
            return None;
        };
        let Some(position_index) = row_group_indexes.get(self.position_leaf) else {
            stats.unavailable_probes += 1;
            return None;
        };
        let Index::INT64(position_index) = position_index else {
            stats.unavailable_probes += 1;
            return None;
        };

        let page_locations = position_offsets.page_locations();
        if page_locations.is_empty() {
            stats.unavailable_probes += 1;
            return None;
        }

        stats.available_probes += 1;
        stats.pages_in_probed_row_groups += page_locations.len() as u64;

        let mut row_ranges = Vec::new();
        for (page_idx, page_index) in position_index.indexes.iter().enumerate() {
            if page_idx >= page_locations.len() {
                continue;
            }
            let Some(min) = page_index.min.and_then(|value| u64::try_from(value).ok()) else {
                continue;
            };
            let Some(max) = page_index.max.and_then(|value| u64::try_from(value).ok()) else {
                continue;
            };
            if position_key < min || position_key > max {
                continue;
            }

            let rows = page_row_count(
                page_locations,
                page_idx,
                self.cursor.row_groups[row_group_id].rows,
            );
            let start = page_locations[page_idx].first_row_index as usize;
            let end = start.saturating_add(rows);
            stats.candidate_pages += 1;
            stats.candidate_rows += rows as u64;
            stats
                .unique_candidate_pages
                .push((row_group_id, page_idx, rows));
            row_ranges.push(start..end);
        }

        if stats.candidate_pages == 0 {
            stats.candidate_misses += 1;
        }

        if row_ranges.len() == 1
            && row_ranges[0].start == 0
            && row_ranges[0].end == self.cursor.row_groups[row_group_id].rows
        {
            None
        } else {
            Some(row_ranges)
        }
    }

    fn record_page_probe_stats(&mut self, stats: PageProbeStats) {
        self.page_index_probes += stats.probes;
        self.page_index_available_probes += stats.available_probes;
        self.page_index_unavailable_probes += stats.unavailable_probes;
        self.page_index_pages_in_probed_row_groups += stats.pages_in_probed_row_groups;
        self.page_index_candidate_pages += stats.candidate_pages;
        self.page_index_candidate_rows += stats.candidate_rows;
        self.page_index_candidate_misses += stats.candidate_misses;
        for (row_group_id, page_idx, rows) in stats.unique_candidate_pages {
            if self
                .page_index_unique_candidate_pages
                .insert((row_group_id, page_idx))
            {
                self.page_index_unique_candidate_rows += rows as u64;
            }
        }
    }

    fn find_row_group_without_stats(&self, position_key: u64) -> Option<usize> {
        self.cursor
            .row_groups
            .binary_search_by(|row_group| {
                if position_key < row_group.min_position_key {
                    std::cmp::Ordering::Greater
                } else if position_key > row_group.max_position_key {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .ok()
    }

    pub fn stats_snapshot(&self) -> ColdParquetStats {
        let cursor_stats = self.cursor.stats_snapshot();
        ColdParquetStats {
            row_groups_total: self.cursor.row_groups.len() as u64,
            row_groups_unique_touched: self.row_groups_unique_touched,
            row_groups_loaded: self.row_groups_loaded,
            row_group_cache_hits: self.row_group_cache_hits,
            row_group_cache_misses: self.row_group_cache_misses,
            rows_loaded: self.rows_loaded,
            load_time: self.load_time,
            row_group_metadata_probes: cursor_stats.metadata_probes,
            row_group_current_hits: cursor_stats.current_hits,
            row_group_previous_hits: cursor_stats.previous_hits,
            row_group_advanced_hits: cursor_stats.advanced_hits,
            row_group_binary_search_hits: cursor_stats.binary_search_hits,
            row_group_metadata_misses: cursor_stats.misses,
            row_group_skipped_ahead: cursor_stats.skipped_ahead,
            position_page_index_loaded: self.page_layout.position_page_index_loaded,
            position_column_index_loaded: self.page_layout.position_column_index_loaded,
            position_pages_total: self.page_layout.position_pages_total,
            position_bloom_filter_row_groups: self.page_layout.position_bloom_filter_row_groups,
            position_bloom_filter_bytes: self.page_layout.position_bloom_filter_bytes,
            page_index_probes: self.page_index_probes,
            page_index_available_probes: self.page_index_available_probes,
            page_index_unavailable_probes: self.page_index_unavailable_probes,
            page_index_pages_in_probed_row_groups: self.page_index_pages_in_probed_row_groups,
            page_index_candidate_pages: self.page_index_candidate_pages,
            page_index_candidate_rows: self.page_index_candidate_rows,
            page_index_candidate_misses: self.page_index_candidate_misses,
            page_index_unique_candidate_pages: self.page_index_unique_candidate_pages.len() as u64,
            page_index_unique_candidate_rows: self.page_index_unique_candidate_rows,
        }
    }
}

impl ColdParquetLookupSet {
    pub fn open<I, P>(
        paths: I,
        projection_columns: Vec<String>,
        batch_size: usize,
        max_cached_row_groups: usize,
    ) -> Result<Self>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
    {
        let mut lookups = paths
            .into_iter()
            .map(|path| {
                ColdParquetLookup::open(
                    path,
                    projection_columns.clone(),
                    batch_size,
                    max_cached_row_groups,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        if lookups.is_empty() {
            return Err(DataFusionError::Execution(
                "no cold variation parquet files found".into(),
            ));
        }
        lookups.sort_by_key(|lookup| {
            lookup
                .position_range()
                .map(|(min, _)| min)
                .unwrap_or_default()
        });
        validate_cold_file_order(&lookups)?;
        Ok(Self { lookups, cursor: 0 })
    }

    pub fn from_env<I, P>(paths: I, projection_columns: Vec<String>) -> Result<Self>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
    {
        Self::open(
            paths,
            projection_columns,
            cold_parquet_batch_size(),
            cold_parquet_row_group_cache(),
        )
    }

    pub fn probe_position_emit_and_visit<P, E, V>(
        &mut self,
        position_key: u64,
        allele_matches: P,
        emit_match: E,
        visit_row: V,
    ) -> Result<ColdProbeResult>
    where
        P: FnMut(&WarmChunkContext, u32, &str) -> Result<bool>,
        E: FnMut(&WarmChunkContext, u32) -> Result<()>,
        V: FnMut(&WarmChunkContext, u32, &str) -> Result<()>,
    {
        let Some(lookup_idx) = self.find_lookup(position_key) else {
            return Ok(ColdProbeResult::NotCovered);
        };
        self.lookups[lookup_idx].probe_position_emit_and_visit(
            position_key,
            allele_matches,
            emit_match,
            visit_row,
        )
    }

    pub fn prefetch_positions<I>(&mut self, position_keys: I) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        self.prefetch_positions_inner(position_keys, false)
    }

    pub fn prefetch_positions_retaining<I>(&mut self, position_keys: I) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        self.prefetch_positions_inner(position_keys, true)
    }

    fn prefetch_positions_inner<I>(
        &mut self,
        position_keys: I,
        retain_prefetched: bool,
    ) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        let mut keys_by_lookup: HashMap<usize, Vec<u64>> = HashMap::new();
        for position_key in position_keys {
            if let Some(lookup_idx) = self.find_lookup_without_stats(position_key) {
                keys_by_lookup
                    .entry(lookup_idx)
                    .or_default()
                    .push(position_key);
            }
        }

        let mut keys_by_lookup = keys_by_lookup.into_iter().collect::<Vec<_>>();
        keys_by_lookup.sort_by_key(|(lookup_idx, _)| *lookup_idx);
        for (lookup_idx, position_keys) in keys_by_lookup {
            if retain_prefetched {
                self.lookups[lookup_idx].prefetch_positions_retaining(position_keys)?;
            } else {
                self.lookups[lookup_idx].prefetch_positions(position_keys)?;
            }
        }

        Ok(())
    }

    pub fn cached_chunk_for_position(&self, position_key: u64) -> Option<&WarmChunkContext> {
        let lookup_idx = self.find_lookup_without_stats(position_key)?;
        self.lookups[lookup_idx].cached_chunk_for_position(position_key)
    }

    pub fn trim_cache_to_capacity(&mut self) {
        for lookup in &mut self.lookups {
            lookup.trim_cache_to_capacity();
        }
    }

    pub fn position_group(&self, position_key: u64) -> Option<ColdParquetPositionGroup> {
        let lookup_idx = self.find_lookup_without_stats(position_key)?;
        let row_group_id = self.lookups[lookup_idx].find_row_group_without_stats(position_key)?;
        Some(ColdParquetPositionGroup {
            lookup_idx,
            row_group_id,
        })
    }

    pub fn prefetch_position_group<I>(
        &mut self,
        group: ColdParquetPositionGroup,
        position_keys: I,
    ) -> Result<()>
    where
        I: IntoIterator<Item = u64>,
    {
        let Some(lookup) = self.lookups.get_mut(group.lookup_idx) else {
            return Ok(());
        };

        let row_group_position_keys = position_keys
            .into_iter()
            .filter(|position_key| {
                lookup
                    .find_row_group_without_stats(*position_key)
                    .is_some_and(|row_group_id| row_group_id == group.row_group_id)
            })
            .collect::<Vec<_>>();
        lookup.prefetch_positions(row_group_position_keys)
    }

    pub fn stats_snapshot(&self) -> ColdParquetStats {
        let mut aggregate = ColdParquetStats::default();
        for lookup in &self.lookups {
            aggregate += lookup.stats_snapshot();
        }
        aggregate
    }

    fn find_lookup(&mut self, position_key: u64) -> Option<usize> {
        if self
            .lookups
            .get(self.cursor)
            .is_some_and(|lookup| lookup.contains_position(position_key))
        {
            return Some(self.cursor);
        }
        if let Some(previous) = self.cursor.checked_sub(1)
            && self.lookups[previous].contains_position(position_key)
        {
            self.cursor = previous;
            return Some(previous);
        }
        while self.cursor < self.lookups.len()
            && self.lookups[self.cursor]
                .position_range()
                .is_some_and(|(_, max)| position_key > max)
        {
            self.cursor += 1;
        }
        if self
            .lookups
            .get(self.cursor)
            .is_some_and(|lookup| lookup.contains_position(position_key))
        {
            return Some(self.cursor);
        }
        let result = self
            .lookups
            .binary_search_by(|lookup| {
                let Some((min, max)) = lookup.position_range() else {
                    return std::cmp::Ordering::Less;
                };
                if position_key < min {
                    std::cmp::Ordering::Greater
                } else if position_key > max {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .ok();
        if let Some(idx) = result {
            self.cursor = idx;
        }
        result
    }

    fn find_lookup_without_stats(&self, position_key: u64) -> Option<usize> {
        self.lookups
            .binary_search_by(|lookup| {
                let Some((min, max)) = lookup.position_range() else {
                    return std::cmp::Ordering::Less;
                };
                if position_key < min {
                    std::cmp::Ordering::Greater
                } else if position_key > max {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .ok()
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ColdParquetStats {
    pub row_groups_total: u64,
    pub row_groups_unique_touched: u64,
    pub row_groups_loaded: u64,
    pub row_group_cache_hits: u64,
    pub row_group_cache_misses: u64,
    pub rows_loaded: u64,
    pub load_time: Duration,
    pub row_group_metadata_probes: u64,
    pub row_group_current_hits: u64,
    pub row_group_previous_hits: u64,
    pub row_group_advanced_hits: u64,
    pub row_group_binary_search_hits: u64,
    pub row_group_metadata_misses: u64,
    pub row_group_skipped_ahead: u64,
    pub position_page_index_loaded: bool,
    pub position_column_index_loaded: bool,
    pub position_pages_total: u64,
    pub position_bloom_filter_row_groups: u64,
    pub position_bloom_filter_bytes: u64,
    pub page_index_probes: u64,
    pub page_index_available_probes: u64,
    pub page_index_unavailable_probes: u64,
    pub page_index_pages_in_probed_row_groups: u64,
    pub page_index_candidate_pages: u64,
    pub page_index_candidate_rows: u64,
    pub page_index_candidate_misses: u64,
    pub page_index_unique_candidate_pages: u64,
    pub page_index_unique_candidate_rows: u64,
}

impl std::ops::AddAssign for ColdParquetStats {
    fn add_assign(&mut self, rhs: Self) {
        self.row_groups_total += rhs.row_groups_total;
        self.row_groups_unique_touched += rhs.row_groups_unique_touched;
        self.row_groups_loaded += rhs.row_groups_loaded;
        self.row_group_cache_hits += rhs.row_group_cache_hits;
        self.row_group_cache_misses += rhs.row_group_cache_misses;
        self.rows_loaded += rhs.rows_loaded;
        self.load_time += rhs.load_time;
        self.row_group_metadata_probes += rhs.row_group_metadata_probes;
        self.row_group_current_hits += rhs.row_group_current_hits;
        self.row_group_previous_hits += rhs.row_group_previous_hits;
        self.row_group_advanced_hits += rhs.row_group_advanced_hits;
        self.row_group_binary_search_hits += rhs.row_group_binary_search_hits;
        self.row_group_metadata_misses += rhs.row_group_metadata_misses;
        self.row_group_skipped_ahead += rhs.row_group_skipped_ahead;
        self.position_page_index_loaded |= rhs.position_page_index_loaded;
        self.position_column_index_loaded |= rhs.position_column_index_loaded;
        self.position_pages_total += rhs.position_pages_total;
        self.position_bloom_filter_row_groups += rhs.position_bloom_filter_row_groups;
        self.position_bloom_filter_bytes += rhs.position_bloom_filter_bytes;
        self.page_index_probes += rhs.page_index_probes;
        self.page_index_available_probes += rhs.page_index_available_probes;
        self.page_index_unavailable_probes += rhs.page_index_unavailable_probes;
        self.page_index_pages_in_probed_row_groups += rhs.page_index_pages_in_probed_row_groups;
        self.page_index_candidate_pages += rhs.page_index_candidate_pages;
        self.page_index_candidate_rows += rhs.page_index_candidate_rows;
        self.page_index_candidate_misses += rhs.page_index_candidate_misses;
        self.page_index_unique_candidate_pages += rhs.page_index_unique_candidate_pages;
        self.page_index_unique_candidate_rows += rhs.page_index_unique_candidate_rows;
    }
}

impl ColdRowGroupCursor {
    pub fn new(row_groups: Vec<ColdRowGroupMeta>) -> Self {
        Self {
            row_groups,
            cursor: 0,
            stats: ColdRowGroupCursorStats::default(),
        }
    }

    pub fn find_row_group(&mut self, position_key: u64) -> Option<usize> {
        self.stats.metadata_probes += 1;
        if self
            .row_groups
            .get(self.cursor)
            .is_some_and(|row_group| row_group.contains(position_key))
        {
            self.stats.current_hits += 1;
            return Some(self.cursor);
        }

        if let Some(previous) = self.cursor.checked_sub(1)
            && self.row_groups[previous].contains(position_key)
        {
            self.stats.previous_hits += 1;
            return Some(previous);
        }

        let start_cursor = self.cursor;
        while self.cursor < self.row_groups.len()
            && position_key > self.row_groups[self.cursor].max_position_key
        {
            self.cursor += 1;
        }
        self.stats.skipped_ahead += (self.cursor - start_cursor) as u64;

        if self
            .row_groups
            .get(self.cursor)
            .is_some_and(|row_group| row_group.contains(position_key))
        {
            self.stats.advanced_hits += 1;
            return Some(self.cursor);
        }

        let result = self
            .row_groups
            .binary_search_by(|row_group| {
                if position_key < row_group.min_position_key {
                    std::cmp::Ordering::Greater
                } else if position_key > row_group.max_position_key {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .inspect(|idx| {
                self.cursor = *idx;
            })
            .ok();
        if result.is_some() {
            self.stats.binary_search_hits += 1;
        } else {
            self.stats.misses += 1;
        }
        result
    }

    pub fn stats_snapshot(&self) -> ColdRowGroupCursorStats {
        self.stats
    }
}

impl ColdRowGroupMeta {
    fn contains(&self, position_key: u64) -> bool {
        position_key >= self.min_position_key && position_key <= self.max_position_key
    }
}

pub fn validate_cold_row_group_order(row_groups: &[ColdRowGroupMeta]) -> Result<()> {
    let mut previous_max = None;
    for (idx, row_group) in row_groups.iter().enumerate() {
        if previous_max.is_some_and(|max| max >= row_group.min_position_key) {
            return Err(DataFusionError::Execution(format!(
                "position_key split across cold row groups: row_group={} previous_max={} current_min={}",
                idx,
                previous_max.unwrap_or_default(),
                row_group.min_position_key
            )));
        }
        previous_max = Some(row_group.max_position_key);
    }
    Ok(())
}

fn validate_cold_file_order(lookups: &[ColdParquetLookup]) -> Result<()> {
    let mut previous_max = None;
    for (idx, lookup) in lookups.iter().enumerate() {
        let Some((min, max)) = lookup.position_range() else {
            continue;
        };
        if previous_max.is_some_and(|previous| previous >= min) {
            return Err(DataFusionError::Execution(format!(
                "position_key split across cold parquet files: file={} previous_max={} current_min={}",
                idx,
                previous_max.unwrap_or_default(),
                min
            )));
        }
        previous_max = Some(max);
    }
    Ok(())
}

pub fn cold_parquet_projection_columns(
    cache_columns: &[String],
    include_colocated: bool,
) -> Vec<String> {
    let mut columns = Vec::with_capacity(cache_columns.len() + 16);
    for name in ["position_key", "allele_string", "end", "failed"] {
        push_unique_column(&mut columns, name);
    }
    for name in cache_columns {
        push_unique_column(&mut columns, name);
    }
    if include_colocated {
        for name in [
            "variation_name",
            "end",
            "failed",
            "somatic",
            "phenotype_or_disease",
            "clin_sig",
            "clin_sig_allele",
            "pubmed",
        ] {
            push_unique_column(&mut columns, name);
        }
        for name in crate::variant_lookup_exec::AF_COL_NAMES {
            push_unique_column(&mut columns, name);
        }
    }
    columns
}

fn cold_position_leaf_index(metadata: &ArrowReaderMetadata) -> Result<usize> {
    metadata
        .parquet_schema()
        .columns()
        .iter()
        .position(|column| column.name() == "position_key")
        .ok_or_else(|| {
            DataFusionError::Execution("cold parquet missing position_key column".into())
        })
}

fn reject_deprecated_cold_variant_keys(metadata: &ArrowReaderMetadata, path: &Path) -> Result<()> {
    if metadata.schema().index_of("variant_keys").is_ok() {
        return Err(DataFusionError::Execution(format!(
            "cold variation parquet '{}' contains deprecated variant_keys; rebuild indexed_parquet cache",
            path.display()
        )));
    }
    Ok(())
}

fn cold_row_group_metadata(
    metadata: &ArrowReaderMetadata,
    position_leaf: usize,
) -> Result<Vec<ColdRowGroupMeta>> {
    (0..metadata.metadata().num_row_groups())
        .map(|row_group| {
            let metadata_row_group = metadata.metadata().row_group(row_group);
            let stats = metadata_row_group
                .column(position_leaf)
                .statistics()
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "cold row group {row_group} missing position_key statistics"
                    ))
                })?;
            let (min_position_key, max_position_key) = match stats {
                Statistics::Int64(stats) => {
                    let min = *stats.min_opt().ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "cold row group {row_group} missing position_key min"
                        ))
                    })?;
                    let max = *stats.max_opt().ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "cold row group {row_group} missing position_key max"
                        ))
                    })?;
                    (min as u64, max as u64)
                }
                other => {
                    return Err(DataFusionError::Execution(format!(
                        "unsupported position_key statistics for cold row group {row_group}: {other:?}"
                    )));
                }
            };

            Ok(ColdRowGroupMeta {
                rows: metadata_row_group.num_rows() as usize,
                min_position_key,
                max_position_key,
            })
        })
        .collect()
}

fn cold_page_layout_stats(
    metadata: &ArrowReaderMetadata,
    position_leaf: usize,
) -> ColdPageLayoutStats {
    let parquet_metadata = metadata.metadata();
    let mut stats = ColdPageLayoutStats {
        position_page_index_loaded: parquet_metadata.offset_index().is_some(),
        position_column_index_loaded: parquet_metadata.column_index().is_some(),
        ..ColdPageLayoutStats::default()
    };

    if let Some(offset_index) = parquet_metadata.offset_index() {
        for row_group in 0..parquet_metadata.num_row_groups() {
            if let Some(position_offsets) = offset_index
                .get(row_group)
                .and_then(|columns| columns.get(position_leaf))
            {
                stats.position_pages_total += position_offsets.page_locations().len() as u64;
            }
        }
    }

    for row_group in 0..parquet_metadata.num_row_groups() {
        let column = parquet_metadata.row_group(row_group).column(position_leaf);
        if column.bloom_filter_offset().is_some() || column.bloom_filter_length().is_some() {
            stats.position_bloom_filter_row_groups += 1;
            stats.position_bloom_filter_bytes += column
                .bloom_filter_length()
                .and_then(|length| u64::try_from(length).ok())
                .unwrap_or(0);
        }
    }

    stats
}

fn page_row_count(
    page_locations: &[parquet::format::PageLocation],
    page_idx: usize,
    row_group_rows: usize,
) -> usize {
    let start = page_locations[page_idx].first_row_index as usize;
    let end = page_locations
        .get(page_idx + 1)
        .map(|location| location.first_row_index as usize)
        .unwrap_or(row_group_rows);
    end.saturating_sub(start)
}

fn merge_row_ranges(ranges: &mut Vec<Range<usize>>) {
    ranges.retain(|range| range.start < range.end);
    ranges.sort_by_key(|range| range.start);

    let mut merged: Vec<Range<usize>> = Vec::with_capacity(ranges.len());
    for range in ranges.drain(..) {
        if let Some(last) = merged.last_mut()
            && range.start <= last.end
        {
            last.end = last.end.max(range.end);
            continue;
        }
        merged.push(range);
    }
    *ranges = merged;
}

fn cached_row_ranges_cover(
    cached: Option<&[Range<usize>]>,
    requested: Option<&[Range<usize>]>,
) -> bool {
    match (cached, requested) {
        // A full-row-group cache satisfies any page-level request.
        (None, _) => true,
        // A page-level cache does not satisfy a full-row-group request.
        (Some(_), None) => false,
        (Some(cached), Some(requested)) => ranges_cover(cached, requested),
    }
}

fn ranges_cover(cached: &[Range<usize>], requested: &[Range<usize>]) -> bool {
    requested.iter().all(|requested| {
        cached
            .iter()
            .any(|cached| cached.start <= requested.start && cached.end >= requested.end)
    })
}

pub fn cold_parquet_batch_size() -> usize {
    std::env::var("VEP_COLD_PARQUET_BATCH_SIZE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_COLD_PARQUET_BATCH_SIZE)
        .max(1)
}

pub fn cold_parquet_row_group_cache() -> usize {
    std::env::var("VEP_COLD_PARQUET_ROW_GROUP_CACHE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_COLD_PARQUET_ROW_GROUP_CACHE)
        .max(1)
}

pub fn cold_parquet_load_page_index() -> bool {
    std::env::var("VEP_COLD_PARQUET_LOAD_PAGE_INDEX")
        .ok()
        .map(|value| matches!(value.as_str(), "1" | "true" | "TRUE" | "yes" | "YES"))
        .unwrap_or(true)
}

fn push_unique_column(columns: &mut Vec<String>, name: &str) {
    if !columns.iter().any(|existing| existing == name) {
        columns.push(name.to_string());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use datafusion::arrow::array::{
        ArrayRef, Int64Array, ListBuilder, RecordBatch, StringArray, UInt64Array, UInt64Builder,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;

    use crate::warm_cache::key::position_key;

    #[test]
    fn cold_cursor_advances_monotonically_and_can_step_back_to_cached_group() {
        let mut cursor = ColdRowGroupCursor::new(vec![
            ColdRowGroupMeta {
                rows: 10,
                min_position_key: 10,
                max_position_key: 20,
            },
            ColdRowGroupMeta {
                rows: 10,
                min_position_key: 21,
                max_position_key: 30,
            },
            ColdRowGroupMeta {
                rows: 10,
                min_position_key: 31,
                max_position_key: 40,
            },
        ]);

        assert_eq!(cursor.find_row_group(12), Some(0));
        assert_eq!(cursor.find_row_group(25), Some(1));
        assert_eq!(cursor.find_row_group(39), Some(2));
        assert_eq!(cursor.find_row_group(26), Some(1));
        assert_eq!(cursor.find_row_group(41), None);

        let stats = cursor.stats_snapshot();
        assert_eq!(stats.metadata_probes, 5);
        assert_eq!(stats.current_hits, 1);
        assert_eq!(stats.advanced_hits, 2);
        assert_eq!(stats.previous_hits, 1);
        assert_eq!(stats.binary_search_hits, 0);
        assert_eq!(stats.misses, 1);
        assert_eq!(stats.skipped_ahead, 3);
    }

    #[test]
    fn cold_cursor_rejects_overlapping_row_group_ranges() {
        let err = validate_cold_row_group_order(&[
            ColdRowGroupMeta {
                rows: 10,
                min_position_key: 10,
                max_position_key: 20,
            },
            ColdRowGroupMeta {
                rows: 10,
                min_position_key: 20,
                max_position_key: 30,
            },
        ])
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("position_key split across cold row groups")
        );
    }

    #[test]
    fn cold_lookup_rejects_deprecated_variant_keys_column() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("cold_with_variant_keys.parquet");
        let mut variant_keys = ListBuilder::new(UInt64Builder::new());
        variant_keys.values().append_value(42);
        variant_keys.append(true);
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new_list(
                "variant_keys",
                Arc::new(Field::new_list_field(DataType::UInt64, true)),
                false,
            ),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("end", DataType::Int64, false),
            Field::new("failed", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![position_key("1", 10).unwrap()])) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![10])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0])) as ArrayRef,
            ],
        )
        .unwrap();
        let file = File::create(&path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let err =
            ColdParquetLookup::open(&path, cold_parquet_projection_columns(&[], false), 64, 2)
                .unwrap_err();

        assert!(err.to_string().contains("deprecated variant_keys"));
    }

    fn open_test_lookup_with_page_index(path: &Path) -> ColdParquetLookup {
        let file = File::open(path).unwrap();
        let metadata_options = ArrowReaderOptions::new().with_page_index(true);
        let metadata = ArrowReaderMetadata::load(&file, metadata_options).unwrap();
        let position_leaf = cold_position_leaf_index(&metadata).unwrap();
        let row_groups = cold_row_group_metadata(&metadata, position_leaf).unwrap();
        let page_layout = cold_page_layout_stats(&metadata, position_leaf);
        let row_group_count = metadata.metadata().num_row_groups();

        ColdParquetLookup {
            path: path.to_path_buf(),
            metadata,
            projection_columns: cold_parquet_projection_columns(&[], false),
            batch_size: 64,
            max_cached_row_groups: 2,
            position_leaf,
            page_layout,
            cursor: ColdRowGroupCursor::new(row_groups),
            cache: VecDeque::new(),
            row_groups_touched: vec![false; row_group_count],
            row_groups_unique_touched: 0,
            row_groups_loaded: 0,
            row_group_cache_hits: 0,
            row_group_cache_misses: 0,
            rows_loaded: 0,
            load_time: Duration::ZERO,
            page_index_probes: 0,
            page_index_available_probes: 0,
            page_index_unavailable_probes: 0,
            page_index_pages_in_probed_row_groups: 0,
            page_index_candidate_pages: 0,
            page_index_candidate_rows: 0,
            page_index_candidate_misses: 0,
            page_index_unique_candidate_pages: HashSet::new(),
            page_index_unique_candidate_rows: 0,
        }
    }

    #[test]
    fn cold_lookup_uses_page_index_row_selection_to_load_candidate_page_only() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("cold.parquet");
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("end", DataType::Int64, false),
            Field::new("failed", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", 10).unwrap(),
                    position_key("1", 10).unwrap(),
                    position_key("1", 20).unwrap(),
                    position_key("1", 20).unwrap(),
                    position_key("1", 30).unwrap(),
                    position_key("1", 30).unwrap(),
                    position_key("1", 40).unwrap(),
                    position_key("1", 40).unwrap(),
                ])) as ArrayRef,
                Arc::new(StringArray::from(vec![
                    "A/G", "A/T", "A/G", "A/T", "A/G", "A/T", "A/G", "A/T",
                ])) as ArrayRef,
                Arc::new(Int64Array::from(vec![10, 10, 20, 20, 30, 30, 40, 40])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0, 0, 0, 0, 0, 0, 0, 0])) as ArrayRef,
            ],
        )
        .unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(batch.num_rows())
            .set_write_batch_size(2)
            .set_data_page_row_count_limit(2)
            .build();
        let file = File::create(&path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props)).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let mut lookup = open_test_lookup_with_page_index(&path);
        assert!(lookup.page_layout.position_page_index_loaded);
        assert!(lookup.page_layout.position_column_index_loaded);
        assert_eq!(lookup.page_layout.position_pages_total, 4);

        let mut emitted = 0;
        let mut visited = Vec::new();
        let result = lookup
            .probe_position_emit_and_visit(
                position_key("1", 30).unwrap(),
                |_, _, allele_string| Ok(allele_string == "A/G"),
                |_, _| {
                    emitted += 1;
                    Ok(())
                },
                |_, row, allele_string| {
                    visited.push((row, allele_string.to_string()));
                    Ok(())
                },
            )
            .unwrap();

        assert_eq!(result, ColdProbeResult::Match);
        assert_eq!(emitted, 1);
        assert_eq!(
            visited,
            vec![(0, "A/G".to_string()), (1, "A/T".to_string())]
        );
        assert_eq!(lookup.page_index_available_probes, 1);
        assert_eq!(lookup.page_index_candidate_pages, 1);
        assert_eq!(lookup.page_index_candidate_rows, 2);
        assert_eq!(lookup.rows_loaded, 2);
    }

    #[test]
    fn cold_lookup_prefetch_merges_page_ranges_for_same_row_group() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("cold.parquet");
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("end", DataType::Int64, false),
            Field::new("failed", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", 10).unwrap(),
                    position_key("1", 10).unwrap(),
                    position_key("1", 20).unwrap(),
                    position_key("1", 20).unwrap(),
                    position_key("1", 30).unwrap(),
                    position_key("1", 30).unwrap(),
                    position_key("1", 40).unwrap(),
                    position_key("1", 40).unwrap(),
                ])) as ArrayRef,
                Arc::new(StringArray::from(vec![
                    "A/G", "A/T", "A/G", "A/T", "A/G", "A/T", "A/G", "A/T",
                ])) as ArrayRef,
                Arc::new(Int64Array::from(vec![10, 10, 20, 20, 30, 30, 40, 40])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0, 0, 0, 0, 0, 0, 0, 0])) as ArrayRef,
            ],
        )
        .unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(batch.num_rows())
            .set_write_batch_size(2)
            .set_data_page_row_count_limit(2)
            .build();
        let file = File::create(&path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props)).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let mut lookup = open_test_lookup_with_page_index(&path);
        lookup
            .prefetch_positions([
                position_key("1", 10).unwrap(),
                position_key("1", 30).unwrap(),
            ])
            .unwrap();

        assert_eq!(lookup.row_groups_loaded, 1);
        assert_eq!(lookup.rows_loaded, 4);

        for key in [
            position_key("1", 10).unwrap(),
            position_key("1", 30).unwrap(),
        ] {
            let result = lookup
                .probe_position_emit_and_visit(
                    key,
                    |_, _, allele_string| Ok(allele_string == "A/G"),
                    |_, _| Ok(()),
                    |_, _, _| Ok(()),
                )
                .unwrap();
            assert_eq!(result, ColdProbeResult::Match);
        }

        assert_eq!(lookup.row_groups_loaded, 1);
        assert_eq!(lookup.row_group_cache_hits, 2);
    }

    #[test]
    fn cold_lookup_prefetch_records_page_index_and_exposes_cached_chunk() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("cold.parquet");
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("end", DataType::Int64, false),
            Field::new("failed", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![
                    position_key("1", 10).unwrap(),
                    position_key("1", 10).unwrap(),
                    position_key("1", 20).unwrap(),
                    position_key("1", 20).unwrap(),
                ])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G", "A/T", "A/G", "A/T"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![10, 10, 20, 20])) as ArrayRef,
                Arc::new(Int64Array::from(vec![0, 0, 0, 0])) as ArrayRef,
            ],
        )
        .unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(batch.num_rows())
            .set_write_batch_size(2)
            .set_data_page_row_count_limit(2)
            .build();
        let file = File::create(&path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props)).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let mut lookup = open_test_lookup_with_page_index(&path);
        let key = position_key("1", 20).unwrap();
        lookup.prefetch_positions_retaining([key]).unwrap();

        assert_eq!(lookup.page_index_probes, 1);
        assert_eq!(lookup.page_index_available_probes, 1);
        assert_eq!(lookup.row_groups_loaded, 1);
        assert_eq!(lookup.rows_loaded, 2);

        let chunk = lookup.cached_chunk_for_position(key).unwrap();
        let rows = chunk.rows_for_position(key);
        assert_eq!(rows, 0..2);
        assert_eq!(chunk.allele_string(0).unwrap(), Some("A/G"));
        assert_eq!(chunk.allele_string(1).unwrap(), Some("A/T"));
        assert_eq!(lookup.page_index_probes, 1);
    }
}
