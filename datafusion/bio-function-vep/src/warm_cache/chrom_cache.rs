use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Duration, Instant};

use datafusion::arrow::compute::concat_batches;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ParquetRecordBatchReaderBuilder};
use parquet::file::statistics::Statistics;

use crate::kv_cache::position_index::{PositionIndex, PositionIndexSource};
use crate::variant_lookup_exec::AF_COL_NAMES;
use crate::warm_cache::chunk::{WarmChunkContext, WarmChunkProbe};
use crate::warm_cache::reader::projection_for_existing_roots;

#[derive(Debug)]
pub enum WarmChromOpen {
    Available(Box<WarmChromCache>),
    Unavailable(String),
}

#[derive(Debug)]
pub struct WarmChromCache {
    path: PathBuf,
    metadata: ArrowReaderMetadata,
    projection_columns: Vec<String>,
    plans: Vec<WarmChunkPlan>,
    current_plan_idx: Option<usize>,
    current: Option<WarmChunkContext>,
    cold_positions: Option<Arc<PositionIndex>>,
    cold_position_source: Option<PositionIndexSource>,
    batch_size: usize,
    pub chunks_loaded: u64,
    pub chunk_rows: u64,
    pub chunk_load: Duration,
    pub chunk_index: Duration,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WarmChunkPlan {
    pub row_groups: Range<usize>,
    pub min_position_key: u64,
    pub max_position_key: u64,
    pub rows: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WarmRowGroupMeta {
    pub rows: usize,
    pub min_position_key: u64,
    pub max_position_key: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WarmProbe {
    Exact { row: u32, position_rows: Range<u32> },
    PositionCoveredNoExact { position_rows: Range<u32> },
    NotCovered,
}

impl WarmChromCache {
    pub fn open(
        path: impl AsRef<Path>,
        chrom: &str,
        index_dir: impl AsRef<Path>,
        cold_parquet: Option<&Path>,
        batch_size: usize,
        projection_columns: Vec<String>,
    ) -> Result<WarmChromOpen> {
        Self::open_with_optional_position_index(
            path,
            chrom,
            Some(index_dir.as_ref()),
            cold_parquet,
            batch_size,
            projection_columns,
        )
    }

    pub fn open_with_optional_position_index(
        path: impl AsRef<Path>,
        chrom: &str,
        index_dir: Option<&Path>,
        cold_parquet: Option<&Path>,
        batch_size: usize,
        projection_columns: Vec<String>,
    ) -> Result<WarmChromOpen> {
        let path = path.as_ref().to_path_buf();
        let file = match File::open(&path) {
            Ok(file) => file,
            Err(error) => return Ok(WarmChromOpen::Unavailable(error.to_string())),
        };
        let metadata = ArrowReaderMetadata::load(&file, Default::default())?;
        let row_groups = warm_row_group_metadata(&metadata)?;
        let plans = plan_warm_chunks_from_row_groups(&row_groups)?;
        let (cold_positions, cold_position_source) = if let Some(index_dir) = index_dir {
            match PositionIndex::shared_for_chrom(index_dir, chrom, cold_parquet, batch_size) {
                Ok((index, source)) => (Some(index), Some(source)),
                Err(error) => return Ok(WarmChromOpen::Unavailable(error.to_string())),
            }
        } else {
            (None, None)
        };

        Ok(WarmChromOpen::Available(Box::new(Self {
            path,
            metadata,
            projection_columns,
            plans,
            current_plan_idx: None,
            current: None,
            cold_positions,
            cold_position_source,
            batch_size,
            chunks_loaded: 0,
            chunk_rows: 0,
            chunk_load: Duration::ZERO,
            chunk_index: Duration::ZERO,
        })))
    }

    pub fn probe<F>(
        &mut self,
        position_key: u64,
        variant_key: u64,
        verify_row: F,
    ) -> Result<WarmProbe>
    where
        F: FnMut(u32, &str) -> Result<bool>,
    {
        let Some(plan_idx) = self.find_plan(position_key) else {
            return Ok(WarmProbe::NotCovered);
        };
        self.ensure_loaded(plan_idx)?;
        let chunk = self
            .current
            .as_ref()
            .ok_or_else(|| DataFusionError::Execution("warm chunk not loaded".into()))?;

        Ok(
            match chunk.probe_exact(position_key, variant_key, verify_row)? {
                WarmChunkProbe::Exact { row, position_rows } => {
                    WarmProbe::Exact { row, position_rows }
                }
                WarmChunkProbe::PositionCoveredNoExact { position_rows } => {
                    WarmProbe::PositionCoveredNoExact { position_rows }
                }
                WarmChunkProbe::NotCovered => WarmProbe::NotCovered,
            },
        )
    }

    pub fn current_chunk(&self) -> Option<&WarmChunkContext> {
        self.current.as_ref()
    }

    #[inline]
    pub fn cold_may_contain(&self, position_key: u64) -> bool {
        self.cold_positions
            .as_ref()
            .is_none_or(|index| index.contains_position_key(position_key))
    }

    pub fn cold_position_source(&self) -> Option<PositionIndexSource> {
        self.cold_position_source
    }

    pub fn cold_position_len(&self) -> usize {
        self.cold_positions
            .as_ref()
            .map_or(0, |positions| positions.len())
    }

    pub fn cold_position_storage_bytes(&self) -> usize {
        self.cold_positions
            .as_ref()
            .map_or(0, |positions| positions.storage_bytes())
    }

    fn find_plan(&self, position_key: u64) -> Option<usize> {
        self.plans
            .binary_search_by(|plan| {
                if position_key < plan.min_position_key {
                    std::cmp::Ordering::Greater
                } else if position_key > plan.max_position_key {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .ok()
    }

    fn ensure_loaded(&mut self, plan_idx: usize) -> Result<()> {
        if self.current_plan_idx == Some(plan_idx) {
            return Ok(());
        }

        let started = Instant::now();
        let chunk = self.load_plan(&self.plans[plan_idx])?;
        self.chunk_load += started.elapsed();
        self.chunk_rows += chunk.batch.num_rows() as u64;
        self.chunks_loaded += 1;
        self.current_plan_idx = Some(plan_idx);
        self.current = Some(chunk);
        Ok(())
    }

    fn load_plan(&self, plan: &WarmChunkPlan) -> Result<WarmChunkContext> {
        debug_assert_eq!(plan.row_groups.end, plan.row_groups.start + 1);
        let mask = ProjectionMask::roots(
            self.metadata.parquet_schema(),
            self.metadata
                .schema()
                .fields()
                .iter()
                .enumerate()
                .filter_map(|(idx, field)| {
                    self.projection_columns
                        .iter()
                        .any(|name| field.name() == name)
                        .then_some(idx)
                }),
        );
        let file = File::open(&self.path).map_err(|error| {
            DataFusionError::Execution(format!("failed to open warm parquet chunk: {error}"))
        })?;
        let reader =
            ParquetRecordBatchReaderBuilder::new_with_metadata(file, self.metadata.clone())
                .with_projection(mask)
                .with_row_groups(vec![plan.row_groups.start])
                .with_batch_size(plan.rows.max(1))
                .build()?;

        let batches = reader.collect::<std::result::Result<Vec<_>, _>>()?;
        let batch = match batches.as_slice() {
            [] => {
                return Err(DataFusionError::Execution(format!(
                    "warm row group {} produced no batches",
                    plan.row_groups.start
                )));
            }
            [single] => single.clone(),
            _ => concat_batches(&batches[0].schema(), batches.iter())
                .map_err(|error| DataFusionError::ArrowError(Box::new(error), None))?,
        };

        WarmChunkContext::try_new(plan.row_groups.start, batch)
    }
}

pub fn plan_warm_chunks_from_row_groups(
    row_groups: &[WarmRowGroupMeta],
) -> Result<Vec<WarmChunkPlan>> {
    let mut plans = Vec::with_capacity(row_groups.len());
    let mut previous_max = None;
    for (idx, row_group) in row_groups.iter().enumerate() {
        if previous_max.is_some_and(|max| max >= row_group.min_position_key) {
            return Err(DataFusionError::Execution(format!(
                "position_key split across row groups: row_group={} previous_max={} current_min={}",
                idx,
                previous_max.unwrap_or_default(),
                row_group.min_position_key
            )));
        }
        plans.push(WarmChunkPlan {
            row_groups: idx..idx + 1,
            min_position_key: row_group.min_position_key,
            max_position_key: row_group.max_position_key,
            rows: row_group.rows,
        });
        previous_max = Some(row_group.max_position_key);
    }
    Ok(plans)
}

fn warm_row_group_metadata(metadata: &ArrowReaderMetadata) -> Result<Vec<WarmRowGroupMeta>> {
    let position_leaf = metadata
        .parquet_schema()
        .columns()
        .iter()
        .position(|column| column.name() == "position_key")
        .ok_or_else(|| {
            DataFusionError::Execution("warm parquet missing position_key column".into())
        })?;

    (0..metadata.metadata().num_row_groups())
        .map(|row_group| {
            let metadata_row_group = metadata.metadata().row_group(row_group);
            let stats = metadata_row_group
                .column(position_leaf)
                .statistics()
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "warm row group {row_group} missing position_key statistics"
                    ))
                })?;
            let (min_position_key, max_position_key) = match stats {
                Statistics::Int64(stats) => {
                    let min = *stats.min_opt().ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "warm row group {row_group} missing position_key min"
                        ))
                    })?;
                    let max = *stats.max_opt().ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "warm row group {row_group} missing position_key max"
                        ))
                    })?;
                    (u64::try_from(min).map_err(|_| {
                        DataFusionError::Execution(format!(
                            "warm row group {row_group} has negative position_key min: {min}"
                        ))
                    })?, u64::try_from(max).map_err(|_| {
                        DataFusionError::Execution(format!(
                            "warm row group {row_group} has negative position_key max: {max}"
                        ))
                    })?)
                }
                other => {
                    return Err(DataFusionError::Execution(format!(
                        "unsupported position_key statistics for warm row group {row_group}: {other:?}"
                    )));
                }
            };

            Ok(WarmRowGroupMeta {
                rows: metadata_row_group.num_rows() as usize,
                min_position_key,
                max_position_key,
            })
        })
        .collect()
}

pub fn projection_columns_for_cache(
    cache_columns: &[String],
    include_colocated: bool,
) -> Vec<String> {
    let mut columns = Vec::with_capacity(cache_columns.len() + 16);
    for name in [
        "position_key",
        "variant_keys",
        "allele_string",
        "end",
        "failed",
    ] {
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
        for name in AF_COL_NAMES {
            push_unique_column(&mut columns, name);
        }
    }
    columns
}

fn push_unique_column(columns: &mut Vec<String>, name: &str) {
    if !columns.iter().any(|existing| existing == name) {
        columns.push(name.to_string());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_plan_uses_one_row_group_per_chunk() {
        let row_groups = vec![
            WarmRowGroupMeta {
                rows: 500_000,
                min_position_key: 10,
                max_position_key: 20,
            },
            WarmRowGroupMeta {
                rows: 500_000,
                min_position_key: 21,
                max_position_key: 30,
            },
            WarmRowGroupMeta {
                rows: 500_000,
                min_position_key: 31,
                max_position_key: 40,
            },
        ];

        let plan = plan_warm_chunks_from_row_groups(&row_groups).unwrap();

        assert_eq!(plan.len(), 3);
        assert_eq!(plan[0].row_groups, 0..1);
        assert_eq!(plan[0].min_position_key, 10);
        assert_eq!(plan[0].max_position_key, 20);
        assert_eq!(plan[1].row_groups, 1..2);
        assert_eq!(plan[2].row_groups, 2..3);
    }

    #[test]
    fn chunk_plan_rejects_position_split_across_row_groups() {
        let row_groups = vec![
            WarmRowGroupMeta {
                rows: 500_000,
                min_position_key: 10,
                max_position_key: 20,
            },
            WarmRowGroupMeta {
                rows: 500_000,
                min_position_key: 20,
                max_position_key: 30,
            },
        ];

        let err = plan_warm_chunks_from_row_groups(&row_groups).unwrap_err();
        assert!(
            err.to_string()
                .contains("position_key split across row groups")
        );
    }
}
