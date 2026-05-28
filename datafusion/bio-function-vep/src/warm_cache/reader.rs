use std::fs::File;
use std::path::{Path, PathBuf};

use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ParquetRecordBatchReaderBuilder};
use parquet::file::statistics::Statistics;
use parquet::schema::types::SchemaDescriptor;

use crate::warm_cache::chunk::WarmChunkContext;

pub const WARM_RUNTIME_COLUMNS: &[&str] = &[
    "position_key",
    "variant_keys",
    "chrom",
    "start",
    "end",
    "allele_string",
    "variation_name",
    "failed",
    "somatic",
    "strand",
    "minor_allele",
    "minor_allele_freq",
    "clin_sig",
    "phenotype_or_disease",
    "clinical_impact",
    "clin_sig_allele",
    "AF",
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_SAS",
    "gnomADg_REMAINING",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_NFE",
    "gnomADe_SAS",
    "gnomADe_MID",
    "gnomADe_REMAINING",
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
    "pubmed",
];

pub const WARM_DIRECT_RUNTIME_COLUMNS: &[&str] = &[
    "position_key",
    "chrom",
    "start",
    "end",
    "allele_string",
    "variation_name",
    "failed",
    "somatic",
    "strand",
    "minor_allele",
    "minor_allele_freq",
    "clin_sig",
    "phenotype_or_disease",
    "clinical_impact",
    "clin_sig_allele",
    "AF",
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_SAS",
    "gnomADg_REMAINING",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_NFE",
    "gnomADe_SAS",
    "gnomADe_MID",
    "gnomADe_REMAINING",
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
    "pubmed",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WarmPositionProbe {
    Match,
    DefinitiveMiss,
    NotCovered,
    Boundary,
}

#[derive(Debug)]
pub struct WarmChunkWindow {
    path: PathBuf,
    metadata: ArrowReaderMetadata,
    batch_size: usize,
    total_row_groups: usize,
    next_row_group: usize,
    projection_columns: Vec<String>,
    build_variant_index: bool,
    row_group_ranges: Vec<Option<(u64, u64)>>,
    previous: Option<WarmChunkContext>,
    current: Option<WarmChunkContext>,
    pub chunks_loaded: u64,
    pub chunk_rows: u64,
    pub chunk_load: std::time::Duration,
}

impl WarmChunkWindow {
    pub fn open(path: impl AsRef<Path>, batch_size: usize) -> Result<Self> {
        Self::open_with_options(
            path,
            batch_size,
            static_projection_columns(WARM_RUNTIME_COLUMNS),
            true,
        )
    }

    pub fn open_direct(
        path: impl AsRef<Path>,
        batch_size: usize,
        cache_columns: &[String],
    ) -> Result<Self> {
        Self::open_with_options(
            path,
            batch_size,
            warm_direct_projection_columns(cache_columns),
            false,
        )
    }

    fn open_with_options(
        path: impl AsRef<Path>,
        batch_size: usize,
        projection_columns: Vec<String>,
        build_variant_index: bool,
    ) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path).map_err(|e| {
            DataFusionError::Execution(format!("failed to open warm parquet file: {e}"))
        })?;
        let metadata = ArrowReaderMetadata::load(&file, Default::default())?;
        let total_row_groups = metadata.metadata().num_row_groups();
        let row_group_ranges = warm_position_key_row_group_ranges(&metadata);
        Ok(Self {
            path,
            metadata,
            batch_size,
            total_row_groups,
            next_row_group: 0,
            projection_columns,
            build_variant_index,
            row_group_ranges,
            previous: None,
            current: None,
            chunks_loaded: 0,
            chunk_rows: 0,
            chunk_load: std::time::Duration::ZERO,
        })
    }

    pub fn probe_position<F>(
        &mut self,
        position_key: u64,
        mut allele_matches: F,
    ) -> Result<WarmPositionProbe>
    where
        F: FnMut(&str) -> bool,
    {
        self.probe_position_emit(
            position_key,
            |_, _, allele_string| Ok(allele_matches(allele_string)),
            |_, _| Ok(()),
        )
    }

    pub fn probe_position_emit<P, E>(
        &mut self,
        position_key: u64,
        mut allele_matches: P,
        mut emit_match: E,
    ) -> Result<WarmPositionProbe>
    where
        P: FnMut(&WarmChunkContext, u32, &str) -> Result<bool>,
        E: FnMut(&WarmChunkContext, u32) -> Result<()>,
    {
        if self.current.is_none() {
            self.seek_to_position(position_key);
            self.load_next()?;
        }

        while self
            .current
            .as_ref()
            .is_some_and(|chunk| position_key > chunk.max_position_key)
        {
            self.seek_to_position(position_key);
            if !self.load_next()? {
                break;
            }
        }

        let chunks = [self.current.as_ref(), self.previous.as_ref()];
        let mut covered = false;
        let mut boundary = false;

        for chunk in chunks.into_iter().flatten() {
            let rows = chunk.rows_for_position(position_key);
            if rows.is_empty() {
                continue;
            }

            covered = true;
            boundary |= self.boundary_may_be_incomplete(chunk, position_key);
        }

        if boundary {
            return Ok(WarmPositionProbe::Boundary);
        }

        for chunk in chunks.into_iter().flatten() {
            let rows = chunk.rows_for_position(position_key);
            if rows.is_empty() {
                continue;
            }

            for row in rows {
                if let Some(allele_string) = chunk.allele_string(row as usize)?
                    && allele_matches(chunk, row, allele_string)?
                {
                    emit_match(chunk, row)?;
                    return Ok(WarmPositionProbe::Match);
                }
            }
        }

        if covered {
            Ok(WarmPositionProbe::DefinitiveMiss)
        } else {
            Ok(WarmPositionProbe::NotCovered)
        }
    }

    pub fn probe_position_emit_and_visit<P, E, V>(
        &mut self,
        position_key: u64,
        mut allele_matches: P,
        mut emit_match: E,
        mut visit_row: V,
    ) -> Result<WarmPositionProbe>
    where
        P: FnMut(&WarmChunkContext, u32, &str) -> Result<bool>,
        E: FnMut(&WarmChunkContext, u32) -> Result<()>,
        V: FnMut(&WarmChunkContext, u32, &str) -> Result<()>,
    {
        if self.current.is_none() {
            self.seek_to_position(position_key);
            self.load_next()?;
        }

        while self
            .current
            .as_ref()
            .is_some_and(|chunk| position_key > chunk.max_position_key)
        {
            self.seek_to_position(position_key);
            if !self.load_next()? {
                break;
            }
        }

        let chunks = [self.current.as_ref(), self.previous.as_ref()];
        let mut covered = false;
        let mut boundary = false;
        let mut matched = false;

        for chunk in chunks.into_iter().flatten() {
            let rows = chunk.rows_for_position(position_key);
            if rows.is_empty() {
                continue;
            }

            covered = true;
            boundary |= self.boundary_may_be_incomplete(chunk, position_key);
        }

        if boundary {
            return Ok(WarmPositionProbe::Boundary);
        }

        for chunk in chunks.into_iter().flatten() {
            let rows = chunk.rows_for_position(position_key);
            if rows.is_empty() {
                continue;
            }

            for row in rows {
                if let Some(allele_string) = chunk.allele_string(row as usize)? {
                    visit_row(chunk, row, allele_string)?;
                    if !matched && allele_matches(chunk, row, allele_string)? {
                        emit_match(chunk, row)?;
                        matched = true;
                    }
                }
            }
        }

        if matched {
            Ok(WarmPositionProbe::Match)
        } else if covered {
            Ok(WarmPositionProbe::DefinitiveMiss)
        } else {
            Ok(WarmPositionProbe::NotCovered)
        }
    }

    fn load_next(&mut self) -> Result<bool> {
        if self.next_row_group >= self.total_row_groups {
            return Ok(false);
        }

        let row_group_id = self.next_row_group;
        self.next_row_group += 1;
        let started = std::time::Instant::now();
        let chunk = load_warm_chunk_row_group_from_metadata(
            &self.path,
            &self.metadata,
            row_group_id,
            self.batch_size,
            &self.projection_columns,
            self.build_variant_index,
        )?;
        self.chunk_load += started.elapsed();
        self.chunk_rows += chunk.batch.num_rows() as u64;
        self.chunks_loaded += 1;
        self.previous = self.current.take();
        self.current = Some(chunk);
        Ok(true)
    }

    fn seek_to_position(&mut self, position_key: u64) {
        while self.next_row_group < self.total_row_groups {
            let Some((_, max_position_key)) = self.row_group_ranges[self.next_row_group] else {
                break;
            };
            if max_position_key >= position_key {
                break;
            }
            self.next_row_group += 1;
        }
    }

    fn boundary_may_be_incomplete(&self, chunk: &WarmChunkContext, position_key: u64) -> bool {
        (position_key == chunk.max_position_key && self.next_row_group < self.total_row_groups)
            || (position_key == chunk.min_position_key && self.previous.is_some())
    }
}

fn warm_position_key_row_group_ranges(metadata: &ArrowReaderMetadata) -> Vec<Option<(u64, u64)>> {
    let position_leaf = metadata
        .parquet_schema()
        .columns()
        .iter()
        .position(|column| column.name() == "position_key");
    let Some(position_leaf) = position_leaf else {
        return vec![None; metadata.metadata().num_row_groups()];
    };

    (0..metadata.metadata().num_row_groups())
        .map(|row_group| {
            let stats = metadata
                .metadata()
                .row_group(row_group)
                .column(position_leaf)
                .statistics()?;
            match stats {
                Statistics::Int64(stats) => Some((
                    (*stats.min_opt()? as i128).try_into().ok()?,
                    (*stats.max_opt()? as i128).try_into().ok()?,
                )),
                _ => None,
            }
        })
        .collect()
}

pub fn load_warm_chunk_row_group(
    path: impl AsRef<Path>,
    row_group_id: usize,
    batch_size: usize,
) -> Result<WarmChunkContext> {
    let projection_columns = static_projection_columns(WARM_RUNTIME_COLUMNS);
    load_warm_chunk_row_group_with_options(
        path,
        row_group_id,
        batch_size,
        &projection_columns,
        true,
    )
}

fn load_warm_chunk_row_group_with_options(
    path: impl AsRef<Path>,
    row_group_id: usize,
    batch_size: usize,
    projection_columns: &[String],
    build_variant_index: bool,
) -> Result<WarmChunkContext> {
    let file = File::open(path.as_ref()).map_err(|e| {
        DataFusionError::Execution(format!("failed to open warm parquet chunk: {e}"))
    })?;
    let metadata = ArrowReaderMetadata::load(&file, Default::default())?;
    load_warm_chunk_row_group_from_metadata(
        path,
        &metadata,
        row_group_id,
        batch_size,
        projection_columns,
        build_variant_index,
    )
}

fn load_warm_chunk_row_group_from_metadata(
    path: impl AsRef<Path>,
    metadata: &ArrowReaderMetadata,
    row_group_id: usize,
    batch_size: usize,
    projection_columns: &[String],
    build_variant_index: bool,
) -> Result<WarmChunkContext> {
    if row_group_id >= metadata.metadata().num_row_groups() {
        return Err(DataFusionError::Execution(format!(
            "warm row group {row_group_id} out of range; file has {} row groups",
            metadata.metadata().num_row_groups()
        )));
    }

    let mask = projection_for_existing_roots(
        metadata.schema(),
        metadata.parquet_schema(),
        projection_columns,
    );
    let file = File::open(path.as_ref()).map_err(|e| {
        DataFusionError::Execution(format!("failed to open warm parquet chunk: {e}"))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::new_with_metadata(file, metadata.clone());
    let reader = builder
        .with_projection(mask)
        .with_row_groups(vec![row_group_id])
        .with_batch_size(batch_size)
        .build()?;

    let batches = reader.collect::<std::result::Result<Vec<_>, _>>()?;
    let batch = match batches.as_slice() {
        [] => {
            return Err(DataFusionError::Execution(format!(
                "warm row group {row_group_id} produced no batches"
            )));
        }
        [single] => single.clone(),
        _ => concat_batches(&batches[0].schema(), batches.iter())
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?,
    };

    if build_variant_index {
        WarmChunkContext::try_new(row_group_id, batch)
    } else {
        WarmChunkContext::try_new_without_variant_index(row_group_id, batch)
    }
}

pub fn projection_for_existing_roots<S: AsRef<str>>(
    arrow_schema: &SchemaRef,
    parquet_schema: &SchemaDescriptor,
    names: &[S],
) -> ProjectionMask {
    let root_indices = arrow_schema
        .fields()
        .iter()
        .enumerate()
        .filter_map(|(idx, field)| {
            names
                .iter()
                .any(|name| field.name().as_str() == name.as_ref())
                .then_some(idx)
        });
    ProjectionMask::roots(parquet_schema, root_indices)
}

fn static_projection_columns(names: &[&str]) -> Vec<String> {
    names.iter().map(|name| (*name).to_string()).collect()
}

fn warm_direct_projection_columns(cache_columns: &[String]) -> Vec<String> {
    let mut columns = Vec::with_capacity(4 + cache_columns.len());
    for name in ["position_key", "allele_string", "end", "failed"] {
        push_unique_column(&mut columns, name);
    }
    for name in cache_columns {
        push_unique_column(&mut columns, name);
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
    use std::fs::File;
    use std::sync::Arc;

    use datafusion::arrow::array::{
        ArrayRef, Int64Array, ListBuilder, RecordBatch, StringArray, UInt64Array, UInt64Builder,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;

    use super::*;
    use crate::warm_cache::key::{position_key, variant_key};

    fn warm_batch(starts: &[i64], alts: &[&str]) -> RecordBatch {
        let mut variant_keys = ListBuilder::new(UInt64Builder::new());
        for (start, allele) in starts.iter().zip(alts.iter()) {
            let (reference, alternate) = allele.split_once('/').unwrap();
            variant_keys
                .values()
                .append_value(variant_key("1", *start, reference, alternate).unwrap());
            variant_keys.append(true);
        }

        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new_list(
                "variant_keys",
                Arc::new(Field::new_list_field(DataType::UInt64, true)),
                false,
            ),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("unused_column", DataType::Utf8, false),
        ]));

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(
                    starts
                        .iter()
                        .map(|start| position_key("1", *start).unwrap())
                        .collect::<Vec<_>>(),
                )) as ArrayRef,
                Arc::new(variant_keys.finish()) as ArrayRef,
                Arc::new(StringArray::from(vec!["1"; starts.len()])) as ArrayRef,
                Arc::new(Int64Array::from(starts.to_vec())) as ArrayRef,
                Arc::new(StringArray::from(alts.to_vec())) as ArrayRef,
                Arc::new(StringArray::from(vec!["ignored"; starts.len()])) as ArrayRef,
            ],
        )
        .unwrap()
    }

    #[test]
    fn reader_loads_requested_row_group_as_warm_chunk() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_warm.parquet");
        let file = File::create(&path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(2)
            .build();
        let batch1 = warm_batch(&[101, 102], &["A/G", "C/T"]);
        let batch2 = warm_batch(&[205], &["G/A"]);
        let mut writer = ArrowWriter::try_new(file, batch1.schema(), Some(props)).unwrap();
        writer.write(&batch1).unwrap();
        writer.write(&batch2).unwrap();
        writer.close().unwrap();

        let chunk = load_warm_chunk_row_group(&path, 1, 1024).unwrap();

        assert_eq!(chunk.row_group_id, 1);
        assert!(chunk.contains_position(position_key("1", 205).unwrap()));
        assert!(!chunk.contains_position(position_key("1", 101).unwrap()));
        assert!(chunk.batch.schema().index_of("unused_column").is_err());
    }

    #[test]
    fn direct_window_omits_variant_keys_but_still_probes_positions() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_warm.parquet");
        let file = File::create(&path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(2)
            .build();
        let batch1 = warm_batch(&[101, 102], &["A/G", "C/T"]);
        let batch2 = warm_batch(&[205], &["G/A"]);
        let mut writer = ArrowWriter::try_new(file, batch1.schema(), Some(props)).unwrap();
        writer.write(&batch1).unwrap();
        writer.write(&batch2).unwrap();
        writer.close().unwrap();

        let mut window = WarmChunkWindow::open_direct(&path, 1024, &[]).unwrap();

        assert_eq!(
            window
                .probe_position(position_key("1", 101).unwrap(), |allele| allele == "A/G")
                .unwrap(),
            WarmPositionProbe::Match
        );
        let current = window.current.as_ref().unwrap();
        assert!(current.batch.schema().index_of("variant_keys").is_err());
        assert!(
            current
                .lookup_variant(variant_key("1", 101, "A", "G").unwrap())
                .is_empty()
        );
    }

    #[test]
    fn warm_window_probes_position_rows_with_boundary_fallback() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_warm.parquet");
        let file = File::create(&path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(2)
            .build();
        let batch1 = warm_batch(&[101, 102], &["A/G", "C/T"]);
        let batch2 = warm_batch(&[205], &["G/A"]);
        let mut writer = ArrowWriter::try_new(file, batch1.schema(), Some(props)).unwrap();
        writer.write(&batch1).unwrap();
        writer.write(&batch2).unwrap();
        writer.close().unwrap();

        let mut window = WarmChunkWindow::open(&path, 1024).unwrap();

        assert_eq!(
            window
                .probe_position(position_key("1", 101).unwrap(), |allele| allele == "A/G")
                .unwrap(),
            WarmPositionProbe::Match
        );
        assert_eq!(
            window
                .probe_position(position_key("1", 101).unwrap(), |_| false)
                .unwrap(),
            WarmPositionProbe::DefinitiveMiss
        );
        assert_eq!(
            window
                .probe_position(position_key("1", 102).unwrap(), |_| false)
                .unwrap(),
            WarmPositionProbe::Boundary
        );
        let mut emitted = false;
        let mut visited = false;
        assert_eq!(
            window
                .probe_position_emit_and_visit(
                    position_key("1", 102).unwrap(),
                    |_, _, allele| Ok(allele == "C/T"),
                    |_, _| {
                        emitted = true;
                        Ok(())
                    },
                    |_, _, _| {
                        visited = true;
                        Ok(())
                    },
                )
                .unwrap(),
            WarmPositionProbe::Boundary
        );
        assert!(!emitted);
        assert!(!visited);
        assert_eq!(
            window
                .probe_position(position_key("1", 999).unwrap(), |_| false)
                .unwrap(),
            WarmPositionProbe::NotCovered
        );
    }

    #[test]
    fn warm_window_seeks_to_first_relevant_row_group() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_warm.parquet");
        let file = File::create(&path).unwrap();
        let props = WriterProperties::builder()
            .set_max_row_group_size(2)
            .build();
        let batch1 = warm_batch(&[101, 102], &["A/G", "C/T"]);
        let batch2 = warm_batch(&[205], &["G/A"]);
        let mut writer = ArrowWriter::try_new(file, batch1.schema(), Some(props)).unwrap();
        writer.write(&batch1).unwrap();
        writer.write(&batch2).unwrap();
        writer.close().unwrap();

        let mut window = WarmChunkWindow::open(&path, 1024).unwrap();

        assert_eq!(
            window
                .probe_position(position_key("1", 205).unwrap(), |allele| allele == "G/A")
                .unwrap(),
            WarmPositionProbe::Match
        );
        assert_eq!(window.chunks_loaded, 1);
        assert_eq!(window.chunk_rows, 1);
        assert_eq!(window.current.as_ref().unwrap().row_group_id, 1);
    }
}
