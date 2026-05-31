use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use datafusion::arrow::array::{
    Array, Int64Array, LargeStringArray, ListArray, RecordBatch, StringArray, StringViewArray,
    UInt64Array,
};
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::schema::types::SchemaDescriptor;

use crate::allele::vcf_to_vep_allele;
use crate::warm_cache::key::variant_key_from_position;

const MAGIC: &[u8; 8] = b"VPVBF01\0";
const DEFAULT_HASH_COUNT: u8 = 7;
const MIN_BITS: u64 = 64;
const BITS_PER_WORD: u64 = u64::BITS as u64;

#[derive(Debug, Clone)]
pub struct VariantBloomIndex {
    bits: Vec<u64>,
    bit_count: u64,
    hash_count: u8,
    inserted: u64,
}

impl VariantBloomIndex {
    pub fn with_expected_items(expected_items: u64, bits_per_key: u32) -> Result<Self> {
        let expected_items = expected_items.max(1);
        let bit_count = expected_items
            .checked_mul(u64::from(bits_per_key.max(1)))
            .ok_or_else(|| DataFusionError::Execution("variant bloom bit count overflow".into()))?
            .max(MIN_BITS);
        let word_count = bit_count.div_ceil(BITS_PER_WORD);
        let word_count = usize::try_from(word_count).map_err(|_| {
            DataFusionError::Execution(format!("variant bloom index too large: {word_count} words"))
        })?;
        let hash_count = hash_count_for_bits_per_key(bits_per_key);
        Ok(Self {
            bits: vec![0; word_count],
            bit_count,
            hash_count,
            inserted: 0,
        })
    }

    pub fn from_parquet(
        path: impl AsRef<Path>,
        batch_size: usize,
        bits_per_key: u32,
    ) -> Result<Self> {
        let file = File::open(path.as_ref()).map_err(io_err)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let expected_items = u64::try_from(builder.metadata().file_metadata().num_rows())
            .unwrap_or(1)
            .saturating_mul(4)
            .max(1);
        let mask = projection_for_existing_roots(
            builder.schema(),
            builder.parquet_schema(),
            &["position_key", "allele_string", "variant_keys"],
        );
        let reader = builder
            .with_projection(mask)
            .with_batch_size(batch_size.max(1))
            .build()?;

        let mut index = Self::with_expected_items(expected_items, bits_per_key)?;
        for batch in reader {
            index.append_batch(&batch?)?;
        }
        Ok(index)
    }

    pub fn insert(&mut self, key: u64) {
        let (h1, h2) = hash_pair(key);
        for i in 0..self.hash_count {
            let bit = h1.wrapping_add(u64::from(i).wrapping_mul(h2)) % self.bit_count;
            let word = (bit / BITS_PER_WORD) as usize;
            let offset = bit % BITS_PER_WORD;
            self.bits[word] |= 1_u64 << offset;
        }
        self.inserted += 1;
    }

    pub fn contains(&self, key: u64) -> bool {
        let (h1, h2) = hash_pair(key);
        for i in 0..self.hash_count {
            let bit = h1.wrapping_add(u64::from(i).wrapping_mul(h2)) % self.bit_count;
            let word = (bit / BITS_PER_WORD) as usize;
            let offset = bit % BITS_PER_WORD;
            if self.bits.get(word).copied().unwrap_or_default() & (1_u64 << offset) == 0 {
                return false;
            }
        }
        true
    }

    pub fn contains_any<I>(&self, keys: I) -> bool
    where
        I: IntoIterator<Item = u64>,
    {
        keys.into_iter().any(|key| self.contains(key))
    }

    pub fn inserted(&self) -> u64 {
        self.inserted
    }

    pub fn bit_count(&self) -> u64 {
        self.bit_count
    }

    pub fn hash_count(&self) -> u8 {
        self.hash_count
    }

    pub fn storage_bytes(&self) -> usize {
        self.bits.len() * std::mem::size_of::<u64>()
    }

    pub fn write_to_path(&self, path: impl AsRef<Path>) -> Result<()> {
        if let Some(parent) = path.as_ref().parent() {
            std::fs::create_dir_all(parent).map_err(io_err)?;
        }
        let mut file = File::create(path.as_ref()).map_err(io_err)?;
        file.write_all(MAGIC).map_err(io_err)?;
        file.write_all(&self.bit_count.to_le_bytes())
            .map_err(io_err)?;
        file.write_all(&[self.hash_count]).map_err(io_err)?;
        file.write_all(&[0_u8; 7]).map_err(io_err)?;
        file.write_all(&self.inserted.to_le_bytes())
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
                "invalid variant bloom index magic in {}",
                path.as_ref().display()
            )));
        }

        let bit_count = read_u64(&mut file)?;
        let mut hash_count = [0_u8; 1];
        file.read_exact(&mut hash_count).map_err(io_err)?;
        let mut reserved = [0_u8; 7];
        file.read_exact(&mut reserved).map_err(io_err)?;
        let inserted = read_u64(&mut file)?;
        let word_count = read_u64(&mut file)?;
        let word_count = usize::try_from(word_count).map_err(|_| {
            DataFusionError::Execution(format!("variant bloom index too large: {word_count} words"))
        })?;
        let byte_len = word_count.checked_mul(8).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "variant bloom index byte length overflow: {word_count} words"
            ))
        })?;
        let mut bytes = vec![0_u8; byte_len];
        file.read_exact(&mut bytes).map_err(io_err)?;
        let mut bits = Vec::with_capacity(word_count);
        for chunk in bytes.chunks_exact(8) {
            bits.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        if bit_count == 0 || hash_count[0] == 0 {
            return Err(DataFusionError::Execution(format!(
                "invalid variant bloom index dimensions in {}",
                path.as_ref().display()
            )));
        }
        Ok(Self {
            bits,
            bit_count,
            hash_count: hash_count[0],
            inserted,
        })
    }

    fn append_batch(&mut self, batch: &RecordBatch) -> Result<()> {
        let schema = batch.schema();
        if let (Ok(position_idx), Ok(allele_idx)) = (
            schema.index_of("position_key"),
            schema.index_of("allele_string"),
        ) {
            append_allele_match_keys(
                self,
                batch.column(position_idx).as_ref(),
                batch.column(allele_idx).as_ref(),
            )?;
        }
        if let Ok(variant_idx) = schema.index_of("variant_keys") {
            append_variant_keys(self, batch.column(variant_idx).as_ref())?;
        }
        Ok(())
    }
}

pub fn variant_bloom_index_file(dir: impl AsRef<Path>, chrom: &str) -> PathBuf {
    dir.as_ref().join(format!("{chrom}.varbf"))
}

pub fn find_variant_bloom_index_file(dir: impl AsRef<Path>, chrom: &str) -> Option<PathBuf> {
    let direct = variant_bloom_index_file(&dir, chrom);
    if direct.is_file() {
        return Some(direct);
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        let stripped = variant_bloom_index_file(&dir, stripped);
        if stripped.is_file() {
            return Some(stripped);
        }
    } else {
        let prefixed = variant_bloom_index_file(&dir, &format!("chr{chrom}"));
        if prefixed.is_file() {
            return Some(prefixed);
        }
    }

    None
}

fn append_variant_keys(index: &mut VariantBloomIndex, array: &dyn Array) -> Result<()> {
    let variant_array = array
        .as_any()
        .downcast_ref::<ListArray>()
        .ok_or_else(|| DataFusionError::Execution("variant_keys must be List<UInt64>".into()))?;
    let variant_values = variant_array
        .values()
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| DataFusionError::Execution("variant_keys values must be UInt64".into()))?;
    let offsets = variant_array.offsets();
    for row in 0..variant_array.len() {
        if variant_array.is_null(row) {
            continue;
        }
        let start = offsets[row] as usize;
        let end = offsets[row + 1] as usize;
        for value_idx in start..end {
            if !variant_values.is_null(value_idx) {
                index.insert(variant_values.value(value_idx));
            }
        }
    }
    Ok(())
}

fn append_allele_match_keys(
    index: &mut VariantBloomIndex,
    position_array: &dyn Array,
    allele_array: &dyn Array,
) -> Result<()> {
    for row in 0..position_array.len() {
        let Some(position_key) = position_key_value(position_array, row)? else {
            continue;
        };
        let Some(allele_string) = string_value(allele_array, row)? else {
            continue;
        };
        let Some((reference, alternates)) = allele_string.split_once('/') else {
            continue;
        };

        for alternate in alternates
            .split('/')
            .filter(|alternate| !alternate.is_empty())
        {
            insert_match_key(index, position_key, reference, alternate);

            let (left_ref, left_alt) = vcf_to_vep_allele(reference, alternate);
            insert_match_key(index, position_key, &left_ref, &left_alt);

            let (right_ref, right_alt) = trim_right_first(reference, alternate);
            insert_match_key(index, position_key, &right_ref, &right_alt);
        }
    }
    Ok(())
}

fn insert_match_key(index: &mut VariantBloomIndex, position_key: u64, reference: &str, alt: &str) {
    index.insert(variant_key_from_position(position_key, reference, alt));
}

fn position_key_value(array: &dyn Array, row: usize) -> Result<Option<u64>> {
    if array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<UInt64Array>() {
        return Ok(Some(array.value(row)));
    }
    if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        return u64::try_from(array.value(row)).map(Some).map_err(|_| {
            DataFusionError::Execution(format!(
                "negative position_key in cold variation parquet: {}",
                array.value(row)
            ))
        });
    }
    Err(DataFusionError::Execution(
        "position_key must be UInt64 or Int64".into(),
    ))
}

fn string_value(array: &dyn Array, row: usize) -> Result<Option<&str>> {
    if array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<StringArray>() {
        return Ok(Some(array.value(row)));
    }
    if let Some(array) = array.as_any().downcast_ref::<StringViewArray>() {
        return Ok(Some(array.value(row)));
    }
    if let Some(array) = array.as_any().downcast_ref::<LargeStringArray>() {
        return Ok(Some(array.value(row)));
    }
    Err(DataFusionError::Execution(
        "allele_string must be Utf8, Utf8View, or LargeUtf8".into(),
    ))
}

fn trim_right_first(ref_allele: &str, alt_allele: &str) -> (String, String) {
    let mut r = ref_allele.as_bytes().to_vec();
    let mut a = alt_allele.as_bytes().to_vec();

    while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
        r.pop();
        a.pop();
    }

    while !r.is_empty() && !a.is_empty() && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
    }

    let ref_trimmed = if r.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8_lossy(&r).into_owned()
    };
    let alt_trimmed = if a.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8_lossy(&a).into_owned()
    };

    (ref_trimmed, alt_trimmed)
}

fn hash_count_for_bits_per_key(bits_per_key: u32) -> u8 {
    if bits_per_key == 0 {
        return DEFAULT_HASH_COUNT;
    }
    ((f64::from(bits_per_key) * std::f64::consts::LN_2).round() as u8).clamp(1, 32)
}

fn hash_pair(key: u64) -> (u64, u64) {
    let key_bytes = key.to_le_bytes();
    let h1 = rapidhash::v3::rapidhash_v3(&key_bytes);
    let mut salted = [0_u8; 16];
    salted[..8].copy_from_slice(&key_bytes);
    salted[8..].copy_from_slice(&0x9E37_79B9_7F4A_7C15_u64.to_le_bytes());
    let h2 = rapidhash::v3::rapidhash_v3(&salted) | 1;
    (h1, h2)
}

fn read_u64(file: &mut File) -> Result<u64> {
    let mut buf = [0_u8; 8];
    file.read_exact(&mut buf).map_err(io_err)?;
    Ok(u64::from_le_bytes(buf))
}

fn projection_for_existing_roots<S: AsRef<str>>(
    arrow_schema: &datafusion::arrow::datatypes::SchemaRef,
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

    #[test]
    fn variant_bloom_index_round_trips_binary_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1.varbf");
        let mut index = VariantBloomIndex::with_expected_items(4, 10).unwrap();

        index.insert(11);
        index.insert(42);
        index.write_to_path(&path).unwrap();

        let loaded = VariantBloomIndex::read_from_path(&path).unwrap();

        assert_eq!(loaded.inserted(), 2);
        assert!(loaded.contains(11));
        assert!(loaded.contains(42));
        assert!(!loaded.contains(99));
    }

    #[test]
    fn variant_bloom_index_checks_any_candidate() {
        let mut index = VariantBloomIndex::with_expected_items(4, 10).unwrap();
        index.insert(42);

        assert!(index.contains_any([7, 42, 99]));
        assert!(!index.contains_any([7, 99]));
    }

    #[test]
    fn variant_bloom_index_file_matches_chrom_variants() {
        let dir = tempfile::tempdir().unwrap();
        let path = variant_bloom_index_file(dir.path(), "chr4");

        assert_eq!(path, dir.path().join("chr4.varbf"));
    }

    #[test]
    fn variant_bloom_index_builds_from_parquet_variant_keys() {
        use std::sync::Arc;

        use datafusion::arrow::array::{ArrayRef, ListBuilder, UInt64Builder};
        use datafusion::arrow::datatypes::{DataType, Field, Schema};
        use parquet::arrow::ArrowWriter;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_cold.parquet");
        let mut keys = ListBuilder::new(UInt64Builder::new());
        keys.values().append_value(11);
        keys.values().append_value(42);
        keys.append(true);
        keys.values().append_value(99);
        keys.append(true);
        let schema = Arc::new(Schema::new(vec![Field::new_list(
            "variant_keys",
            Arc::new(Field::new_list_field(DataType::UInt64, true)),
            false,
        )]));
        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(keys.finish()) as ArrayRef])
            .unwrap();
        let file = File::create(&path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let index = VariantBloomIndex::from_parquet(&path, 1024, 10).unwrap();

        assert_eq!(index.inserted(), 3);
        assert!(index.contains(11));
        assert!(index.contains(42));
        assert!(index.contains(99));
    }

    #[test]
    fn variant_bloom_index_builds_augmented_keys_from_allele_string() {
        use std::sync::Arc;

        use datafusion::arrow::array::{ArrayRef, StringArray};
        use datafusion::arrow::datatypes::{DataType, Field, Schema};
        use parquet::arrow::ArrowWriter;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr1_cold.parquet");
        let position_key = 101;
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(UInt64Array::from(vec![position_key])) as ArrayRef,
                Arc::new(StringArray::from(vec!["GCC/GCCCAGCC"])) as ArrayRef,
            ],
        )
        .unwrap();
        let file = File::create(&path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let index = VariantBloomIndex::from_parquet(&path, 1024, 10).unwrap();

        assert!(index.contains(variant_key_from_position(position_key, "-", "GCCCA")));
    }
}
