use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, LazyLock, Mutex};

use datafusion::arrow::array::{
    Array, Int32Array, Int64Array, RecordBatch, UInt32Array, UInt64Array,
};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::schema::types::SchemaDescriptor;

const MAGIC: &[u8; 8] = b"VPBIT01\0";
const POSITION_KEY_BITS: u64 = 48;
const POSITION_KEY_MASK: u64 = (1_u64 << POSITION_KEY_BITS) - 1;
const BITS_PER_WORD: usize = u64::BITS as usize;

static SHARED_PARQUET_INDICES: LazyLock<
    Mutex<std::collections::HashMap<PathBuf, Arc<SharedPositionIndexSlot>>>,
> = LazyLock::new(|| Mutex::new(std::collections::HashMap::new()));

#[derive(Debug, Default)]
struct SharedPositionIndexSlot {
    value: Mutex<Option<Arc<PositionIndex>>>,
}

#[derive(Debug, Clone, Default)]
pub struct PositionIndex {
    bits: Vec<u64>,
    len: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PositionIndexSource {
    Persisted,
    ParquetFallback,
}

impl PositionIndex {
    pub fn from_positions<I>(positions: I) -> Result<Self>
    where
        I: IntoIterator<Item = i64>,
    {
        let mut index = Self::default();
        for position in positions {
            index.insert_position(position)?;
        }
        Ok(index)
    }

    pub fn from_position_keys<I>(position_keys: I) -> Result<Self>
    where
        I: IntoIterator<Item = u64>,
    {
        let mut index = Self::default();
        for position_key in position_keys {
            index.insert_position_key(position_key)?;
        }
        Ok(index)
    }

    pub fn from_parquet(path: impl AsRef<Path>, batch_size: usize) -> Result<Self> {
        let file = File::open(path.as_ref()).map_err(io_err)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let projection_columns = if builder.schema().index_of("position_key").is_ok() {
            vec!["position_key"]
        } else {
            vec!["start"]
        };
        let mask = projection_for_existing_roots(
            builder.schema(),
            builder.parquet_schema(),
            &projection_columns,
        );
        let reader = builder
            .with_projection(mask)
            .with_batch_size(batch_size)
            .build()?;

        let mut index = Self::default();
        for batch in reader {
            index.append_batch(&batch?)?;
        }
        Ok(index)
    }

    pub fn shared_from_parquet(
        path: impl AsRef<Path>,
        batch_size: usize,
    ) -> Result<(Arc<Self>, bool)> {
        let path = path.as_ref().to_path_buf();
        let slot = {
            let mut indices = SHARED_PARQUET_INDICES
                .lock()
                .map_err(|_| DataFusionError::Execution("position index cache poisoned".into()))?;
            indices
                .entry(path.clone())
                .or_insert_with(|| Arc::new(SharedPositionIndexSlot::default()))
                .clone()
        };

        let mut value = slot
            .value
            .lock()
            .map_err(|_| DataFusionError::Execution("position index slot poisoned".into()))?;
        if let Some(index) = value.as_ref() {
            return Ok((Arc::clone(index), false));
        }

        let index = Arc::new(Self::from_parquet(&path, batch_size)?);
        *value = Some(Arc::clone(&index));
        Ok((index, true))
    }

    pub fn shared_for_chrom(
        index_dir: impl AsRef<Path>,
        chrom: &str,
        cold_parquet: Option<&Path>,
        batch_size: usize,
    ) -> Result<(Arc<Self>, PositionIndexSource)> {
        if let Some(path) = find_position_index_file(index_dir, chrom) {
            return Ok((
                Arc::new(Self::read_from_path(path)?),
                PositionIndexSource::Persisted,
            ));
        }

        if let Some(path) = cold_parquet {
            let (index, _) = Self::shared_from_parquet(path, batch_size)?;
            return Ok((index, PositionIndexSource::ParquetFallback));
        }

        Err(DataFusionError::Execution(format!(
            "missing cold position source for {chrom}: no persisted .posidx and no cold parquet fallback"
        )))
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn storage_bytes(&self) -> usize {
        self.bits.len() * std::mem::size_of::<u64>()
    }

    pub fn contains(&self, position: i64) -> bool {
        let Ok(position) = usize::try_from(position) else {
            return false;
        };
        let word = position / BITS_PER_WORD;
        let Some(bits) = self.bits.get(word) else {
            return false;
        };
        let bit = position % BITS_PER_WORD;
        bits & (1_u64 << bit) != 0
    }

    pub fn contains_position_key(&self, position_key: u64) -> bool {
        let Ok(position) = i64::try_from(position_key & POSITION_KEY_MASK) else {
            return false;
        };
        self.contains(position)
    }

    pub fn write_to_path(&self, path: impl AsRef<Path>) -> Result<()> {
        if let Some(parent) = path.as_ref().parent() {
            std::fs::create_dir_all(parent).map_err(io_err)?;
        }
        let mut file = File::create(path.as_ref()).map_err(io_err)?;
        file.write_all(MAGIC).map_err(io_err)?;
        file.write_all(&(self.len as u64).to_le_bytes())
            .map_err(io_err)?;
        file.write_all(&(self.bits.len() as u64).to_le_bytes())
            .map_err(io_err)?;
        for word in &self.bits {
            file.write_all(&word.to_le_bytes()).map_err(io_err)?;
        }
        Ok(())
    }

    pub fn read_from_path(path: impl AsRef<Path>) -> Result<Self> {
        let mut file = File::open(path.as_ref()).map_err(io_err)?;
        let mut magic = [0_u8; 8];
        file.read_exact(&mut magic).map_err(io_err)?;
        if &magic != MAGIC {
            return Err(DataFusionError::Execution(format!(
                "invalid variation position index magic in {}",
                path.as_ref().display()
            )));
        }

        let mut count_buf = [0_u8; 8];
        file.read_exact(&mut count_buf).map_err(io_err)?;
        let count = u64::from_le_bytes(count_buf);
        let len = usize::try_from(count).map_err(|_| {
            DataFusionError::Execution(format!(
                "variation position index too large: {} entries",
                count
            ))
        })?;

        let mut word_count_buf = [0_u8; 8];
        file.read_exact(&mut word_count_buf).map_err(io_err)?;
        let word_count = u64::from_le_bytes(word_count_buf);
        let word_count = usize::try_from(word_count).map_err(|_| {
            DataFusionError::Execution(format!(
                "variation position index too large: {} words",
                word_count
            ))
        })?;
        let byte_len = word_count.checked_mul(8).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "variation position index byte length overflow: {} words",
                word_count
            ))
        })?;
        let mut bytes = vec![0_u8; byte_len];
        file.read_exact(&mut bytes).map_err(io_err)?;

        let mut bits = Vec::with_capacity(word_count);
        for chunk in bytes.chunks_exact(8) {
            bits.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        Ok(Self { bits, len })
    }

    fn append_batch(&mut self, batch: &RecordBatch) -> Result<()> {
        let schema = batch.schema();
        if let Ok(position_key_idx) = schema.index_of("position_key") {
            append_position_keys(self, batch.column(position_key_idx).as_ref())
        } else {
            let start_idx = schema.index_of("start")?;
            append_positions(self, batch.column(start_idx).as_ref())
        }
    }

    fn insert_position_key(&mut self, position_key: u64) -> Result<bool> {
        let position = i64::try_from(position_key & POSITION_KEY_MASK).map_err(|_| {
            DataFusionError::Execution(format!("position key out of range: {position_key}"))
        })?;
        self.insert_position(position)
    }

    fn insert_position(&mut self, position: i64) -> Result<bool> {
        let position = u32::try_from(position).map_err(|_| {
            DataFusionError::Execution(format!("position index value out of u32 range: {position}"))
        })? as usize;
        let word = position / BITS_PER_WORD;
        if word >= self.bits.len() {
            self.bits.resize(word + 1, 0);
        }
        let mask = 1_u64 << (position % BITS_PER_WORD);
        let was_absent = self.bits[word] & mask == 0;
        if was_absent {
            self.bits[word] |= mask;
            self.len += 1;
        }
        Ok(was_absent)
    }
}

pub fn position_index_file(dir: impl AsRef<Path>, chrom: &str) -> std::path::PathBuf {
    dir.as_ref().join(format!("{chrom}.posidx"))
}

pub fn find_position_index_file(dir: impl AsRef<Path>, chrom: &str) -> Option<std::path::PathBuf> {
    let direct = position_index_file(&dir, chrom);
    if direct.is_file() {
        return Some(direct);
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        let stripped = position_index_file(&dir, stripped);
        if stripped.is_file() {
            return Some(stripped);
        }
    } else {
        let prefixed = position_index_file(&dir, &format!("chr{chrom}"));
        if prefixed.is_file() {
            return Some(prefixed);
        }
    }

    None
}

pub fn cold_variation_file_for_chrom(
    variation_dir: impl AsRef<Path>,
    chrom: &str,
) -> Option<std::path::PathBuf> {
    variation_split_file_for_chrom(variation_dir, chrom, "cold")
}

fn variation_split_file_for_chrom(
    variation_dir: impl AsRef<Path>,
    chrom: &str,
    suffix: &str,
) -> Option<std::path::PathBuf> {
    let direct = variation_dir
        .as_ref()
        .join(format!("{chrom}_{suffix}.parquet"));
    if direct.is_file() {
        return Some(direct);
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        let stripped = variation_dir
            .as_ref()
            .join(format!("{stripped}_{suffix}.parquet"));
        if stripped.is_file() {
            return Some(stripped);
        }
    } else {
        let prefixed = variation_dir
            .as_ref()
            .join(format!("chr{chrom}_{suffix}.parquet"));
        if prefixed.is_file() {
            return Some(prefixed);
        }
    }

    None
}

fn append_position_keys(index: &mut PositionIndex, array: &dyn Array) -> Result<()> {
    if let Some(array) = array.as_any().downcast_ref::<UInt64Array>() {
        for row in 0..array.len() {
            if !array.is_null(row) {
                index.insert_position_key(array.value(row))?;
            }
        }
    } else if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        for row in 0..array.len() {
            if !array.is_null(row) {
                let value = u64::try_from(array.value(row)).map_err(|_| {
                    DataFusionError::Execution(format!(
                        "negative position_key in cold variation parquet: {}",
                        array.value(row)
                    ))
                })?;
                index.insert_position_key(value)?;
            }
        }
    } else {
        return Err(DataFusionError::Execution(
            "position_key must be UInt64 or Int64".into(),
        ));
    }
    Ok(())
}

fn append_positions(index: &mut PositionIndex, array: &dyn Array) -> Result<()> {
    for row in 0..array.len() {
        if let Some(position) = int64_value(array, row) {
            index.insert_position(position)?;
        }
    }
    Ok(())
}

fn int64_value(array: &dyn Array, row: usize) -> Option<i64> {
    if array.is_null(row) {
        return None;
    }
    if let Some(array) = array.as_any().downcast_ref::<Int32Array>() {
        Some(array.value(row) as i64)
    } else if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        Some(array.value(row))
    } else if let Some(array) = array.as_any().downcast_ref::<UInt32Array>() {
        Some(array.value(row) as i64)
    } else {
        array
            .as_any()
            .downcast_ref::<UInt64Array>()
            .and_then(|array| i64::try_from(array.value(row)).ok())
    }
}

fn projection_for_existing_roots<S: AsRef<str>>(
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

fn io_err(error: std::io::Error) -> DataFusionError {
    DataFusionError::Execution(error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use datafusion::arrow::array::{ArrayRef, UInt64Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;

    #[test]
    fn position_index_sorts_deduplicates_and_checks_bounds() {
        let index = PositionIndex::from_positions([42, 7, 42, 100]).unwrap();

        assert!(index.contains(7));
        assert!(index.contains(42));
        assert!(index.contains(100));
        assert!(!index.contains(8));
        assert!(!index.contains(-1));
        assert!(!index.contains(i64::from(u32::MAX) + 1));
        assert_eq!(index.len(), 3);
    }

    #[test]
    fn position_index_uses_compact_bitset_storage() {
        let index = PositionIndex::from_positions([0, 1, 64, 64, 130]).unwrap();

        assert_eq!(index.len(), 4);
        assert_eq!(index.storage_bytes(), 24);
        assert!(index.contains(130));
        assert!(!index.contains(129));
    }

    #[test]
    fn position_index_round_trips_binary_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.posidx");
        let index = PositionIndex::from_positions([2, 1, 1, 3]).unwrap();

        index.write_to_path(&path).unwrap();
        let loaded = PositionIndex::read_from_path(&path).unwrap();

        assert_eq!(loaded.len(), 3);
        assert!(loaded.contains(1));
        assert!(loaded.contains(2));
        assert!(loaded.contains(3));
        assert!(!loaded.contains(4));
    }

    #[test]
    fn position_index_checks_position_key_membership() {
        let index = PositionIndex::from_position_keys([
            (1_u64 << POSITION_KEY_BITS) | 20,
            (1_u64 << POSITION_KEY_BITS) | 40,
        ])
        .unwrap();

        assert!(index.contains_position_key((1_u64 << POSITION_KEY_BITS) | 20));
        assert!(!index.contains_position_key((1_u64 << POSITION_KEY_BITS) | 21));
    }

    #[test]
    fn shared_loader_prefers_persisted_posidx() {
        let dir = tempfile::tempdir().unwrap();
        let index_dir = dir.path().join("variation.position_index");
        let path = position_index_file(&index_dir, "chr1");
        PositionIndex::from_positions([10, 20, 30])
            .unwrap()
            .write_to_path(&path)
            .unwrap();

        let (index, source) =
            PositionIndex::shared_for_chrom(&index_dir, "chr1", None, 1024).unwrap();

        assert_eq!(source, PositionIndexSource::Persisted);
        assert!(index.contains(20));
        assert!(!index.contains(21));
    }

    #[test]
    fn shared_loader_builds_from_cold_parquet_when_posidx_missing() {
        let dir = tempfile::tempdir().unwrap();
        let index_dir = dir.path().join("variation.position_index");
        let cold_parquet = write_test_position_parquet(
            &dir.path().join("chr1_cold.parquet"),
            &[
                (1_u64 << POSITION_KEY_BITS) | 10,
                (1_u64 << POSITION_KEY_BITS) | 20,
                (1_u64 << POSITION_KEY_BITS) | 20,
                (1_u64 << POSITION_KEY_BITS) | 30,
            ],
        );

        let (index, source) =
            PositionIndex::shared_for_chrom(&index_dir, "chr1", Some(&cold_parquet), 1024).unwrap();

        assert_eq!(source, PositionIndexSource::ParquetFallback);
        assert!(index.contains(20));
        assert!(!index.contains(21));
    }

    #[test]
    fn shared_loader_errors_when_no_position_source_exists() {
        let dir = tempfile::tempdir().unwrap();
        let err = PositionIndex::shared_for_chrom(
            dir.path().join("variation.position_index"),
            "chr1",
            None,
            1024,
        )
        .unwrap_err();

        assert!(err.to_string().contains("missing cold position source"));
    }

    fn write_test_position_parquet(path: &Path, positions: &[u64]) -> PathBuf {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "position_key",
            DataType::UInt64,
            false,
        )]));
        let batch = RecordBatch::try_new(
            Arc::clone(&schema),
            vec![Arc::new(UInt64Array::from(positions.to_vec())) as ArrayRef],
        )
        .unwrap();
        let file = File::create(path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
        path.to_path_buf()
    }
}
