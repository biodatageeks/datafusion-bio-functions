//! Per-position entry serialization for V5 format.
//!
//! Each genomic position `(chrom, start, end)` stores one binary value in fjall
//! containing an allele table and column-major data for all cache columns.
//!
//! Binary layout:
//! ```text
//! HEADER (4 bytes):
//!   [2B num_alleles u16 LE]
//!   [2B num_cols u16 LE]
//!
//! ALLELE TABLE (variable):
//!   For each allele (num_alleles entries):
//!     [2B allele_len u16 LE]
//!     [allele_len bytes UTF-8]
//!
//! COLUMN DATA (num_cols entries):
//!   For each column:
//!     [1B type_code]
//!     [Packed values for num_alleles rows]
//! ```

use datafusion::arrow::array::{
    Array, ArrayBuilder, BooleanArray, BooleanBuilder, Float32Array, Float32Builder, Float64Array,
    Float64Builder, Int8Array, Int8Builder, Int16Array, Int16Builder, Int32Array, Int32Builder,
    Int64Array, Int64Builder, RecordBatch, StringBuilder, UInt16Array, UInt16Builder, UInt32Array,
    UInt32Builder, UInt64Array, UInt64Builder,
};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::{DataFusionError, Result};

// Reuse the same type codes as V4 (from mmap_block_store.rs).
fn encode_data_type(dt: &DataType) -> Result<u8> {
    match dt {
        DataType::Int32 => Ok(0),
        DataType::Int64 => Ok(1),
        DataType::Utf8 | DataType::Utf8View => Ok(2),
        DataType::Boolean => Ok(3),
        DataType::Float32 => Ok(4),
        DataType::Float64 => Ok(5),
        DataType::UInt32 => Ok(6),
        DataType::UInt64 => Ok(7),
        DataType::LargeUtf8 => Ok(8),
        DataType::Int8 => Ok(9),
        DataType::Int16 => Ok(10),
        DataType::UInt8 => Ok(11),
        DataType::UInt16 => Ok(12),
        other => Err(DataFusionError::Execution(format!(
            "V5 format does not support data type: {other}"
        ))),
    }
}

fn decode_data_type(code: u8) -> Result<DataType> {
    match code {
        0 => Ok(DataType::Int32),
        1 => Ok(DataType::Int64),
        2 => Ok(DataType::Utf8),
        3 => Ok(DataType::Boolean),
        4 => Ok(DataType::Float32),
        5 => Ok(DataType::Float64),
        6 => Ok(DataType::UInt32),
        7 => Ok(DataType::UInt64),
        8 => Ok(DataType::LargeUtf8),
        9 => Ok(DataType::Int8),
        10 => Ok(DataType::Int16),
        11 => Ok(DataType::UInt8),
        12 => Ok(DataType::UInt16),
        other => Err(DataFusionError::Execution(format!(
            "V5 unknown data type code: {other}"
        ))),
    }
}

/// Serialize a position entry from selected rows of a RecordBatch.
///
/// `rows`: row indices in the batch belonging to this position group.
/// `batch`: source RecordBatch.
/// `col_indices`: which columns to serialize (indices into `batch`).
/// `allele_col_idx`: index of the `allele_string` column in `batch`.
pub fn serialize_position_entry(
    rows: &[usize],
    batch: &RecordBatch,
    col_indices: &[usize],
    allele_col_idx: usize,
) -> Result<Vec<u8>> {
    let num_alleles = rows.len() as u16;
    let num_cols = col_indices.len() as u16;

    let mut buf = Vec::with_capacity(256);

    // Header.
    buf.extend_from_slice(&num_alleles.to_le_bytes());
    buf.extend_from_slice(&num_cols.to_le_bytes());

    // Allele table.
    let allele_col = batch.column(allele_col_idx);
    for &row in rows {
        let allele_str = get_string_value(allele_col.as_ref(), row);
        let len = allele_str.len() as u16;
        buf.extend_from_slice(&len.to_le_bytes());
        buf.extend_from_slice(allele_str.as_bytes());
    }

    // Column data.
    for &col_idx in col_indices {
        let col = batch.column(col_idx);
        let dt = col.data_type();
        buf.push(encode_data_type(dt)?);
        serialize_column_values(&mut buf, col.as_ref(), rows)?;
    }

    Ok(buf)
}

/// Serialize column values for a set of rows into the buffer.
fn serialize_column_values(buf: &mut Vec<u8>, col: &dyn Array, rows: &[usize]) -> Result<()> {
    let n = rows.len();

    // Null bitmap: ceil(n/8) bytes, bit set = non-null.
    let null_bytes = n.div_ceil(8);
    let null_start = buf.len();
    buf.resize(null_start + null_bytes, 0);
    for (bit_idx, &row) in rows.iter().enumerate() {
        if !col.is_null(row) {
            buf[null_start + bit_idx / 8] |= 1 << (bit_idx % 8);
        }
    }

    match col.data_type() {
        DataType::Int32 => {
            let arr = col.as_any().downcast_ref::<Int32Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::Int64 => {
            let arr = col.as_any().downcast_ref::<Int64Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::Float32 => {
            let arr = col.as_any().downcast_ref::<Float32Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) {
                    0.0f32
                } else {
                    arr.value(row)
                };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::Float64 => {
            let arr = col.as_any().downcast_ref::<Float64Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) {
                    0.0f64
                } else {
                    arr.value(row)
                };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::UInt32 => {
            let arr = col.as_any().downcast_ref::<UInt32Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::UInt64 => {
            let arr = col.as_any().downcast_ref::<UInt64Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::Int8 => {
            let arr = col.as_any().downcast_ref::<Int8Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.push(v as u8);
            }
        }
        DataType::Int16 => {
            let arr = col.as_any().downcast_ref::<Int16Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::UInt8 => {
            let arr = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::UInt8Array>()
                .unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.push(v);
            }
        }
        DataType::UInt16 => {
            let arr = col.as_any().downcast_ref::<UInt16Array>().unwrap();
            for &row in rows {
                let v = if arr.is_null(row) { 0 } else { arr.value(row) };
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        DataType::Boolean => {
            let arr = col.as_any().downcast_ref::<BooleanArray>().unwrap();
            let value_bytes = n.div_ceil(8);
            let value_start = buf.len();
            buf.resize(value_start + value_bytes, 0);
            for (bit_idx, &row) in rows.iter().enumerate() {
                if !arr.is_null(row) && arr.value(row) {
                    buf[value_start + bit_idx / 8] |= 1 << (bit_idx % 8);
                }
            }
        }
        DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => {
            serialize_string_values(buf, col, rows)?;
        }
        other => {
            return Err(DataFusionError::Execution(format!(
                "V5 serialize: unsupported data type {other}"
            )));
        }
    }
    Ok(())
}

/// Serialize string column values.
///
/// Layout: `[4B total_string_bytes][N × 4B cumulative offset][string bytes]`
fn serialize_string_values(buf: &mut Vec<u8>, col: &dyn Array, rows: &[usize]) -> Result<()> {
    // Collect strings first.
    let mut strings: Vec<&str> = Vec::with_capacity(rows.len());
    let mut total_bytes = 0u32;
    for &row in rows {
        let s = get_string_value(col, row);
        total_bytes += s.len() as u32;
        strings.push(s);
    }

    buf.extend_from_slice(&total_bytes.to_le_bytes());

    // Cumulative offsets.
    let mut cumulative = 0u32;
    for s in &strings {
        cumulative += s.len() as u32;
        buf.extend_from_slice(&cumulative.to_le_bytes());
    }

    // String bytes.
    for s in &strings {
        buf.extend_from_slice(s.as_bytes());
    }

    Ok(())
}

/// Zero-copy reader for a serialized position entry.
pub struct PositionEntryReader<'a> {
    data: &'a [u8],
    num_alleles: usize,
    num_cols: usize,
    /// Precomputed byte offsets into the allele table: (offset, len) for each allele.
    allele_offsets: Vec<(usize, usize)>,
    /// Precomputed (byte_offset_after_type_byte, type_code) for each column.
    col_offsets: Vec<(usize, u8)>,
}

impl<'a> PositionEntryReader<'a> {
    /// Parse header and build allele offset table.
    pub fn new(data: &'a [u8]) -> Result<Self> {
        if data.len() < 4 {
            return Err(DataFusionError::Execution(
                "V5 position entry too short for header".into(),
            ));
        }
        let num_alleles = u16::from_le_bytes([data[0], data[1]]) as usize;
        let num_cols = u16::from_le_bytes([data[2], data[3]]) as usize;

        let mut offset = 4;
        let mut allele_offsets = Vec::with_capacity(num_alleles);
        for _ in 0..num_alleles {
            if offset + 2 > data.len() {
                return Err(DataFusionError::Execution(
                    "V5 position entry: truncated allele table".into(),
                ));
            }
            let len = u16::from_le_bytes([data[offset], data[offset + 1]]) as usize;
            offset += 2;
            allele_offsets.push((offset, len));
            offset += len;
        }

        // Precompute column offsets to avoid O(N^2) sequential walks.
        let col_data_offset = offset;
        let null_bytes = num_alleles.div_ceil(8);
        let mut col_offsets = Vec::with_capacity(num_cols);
        let mut col_off = col_data_offset;
        for col in 0..num_cols {
            if col_off >= data.len() {
                return Err(DataFusionError::Execution(format!(
                    "V5 position entry: truncated column data at col {col}"
                )));
            }
            let type_code = data[col_off];
            col_offsets.push((col_off + 1, type_code)); // +1 to skip type byte
            col_off += 1; // type byte
            let dt = decode_data_type(type_code)?;
            col_off += column_packed_size(&dt, num_alleles, null_bytes, &data[col_off..])?;
        }

        Ok(Self {
            data,
            num_alleles,
            num_cols,
            allele_offsets,
            col_offsets,
        })
    }

    pub fn num_alleles(&self) -> usize {
        self.num_alleles
    }

    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    /// Get the allele string for a given allele index.
    pub fn allele_string(&self, idx: usize) -> &str {
        let (off, len) = self.allele_offsets[idx];
        std::str::from_utf8(&self.data[off..off + len]).unwrap_or("")
    }

    /// Append selected allele rows for a given column into an ArrayBuilder.
    ///
    /// `col_idx`: which column (0-based) in the entry.
    /// `allele_rows`: which allele indices to extract.
    /// `builder`: target builder matching the column type.
    pub fn append_column_values(
        &self,
        col_idx: usize,
        allele_rows: &[usize],
        builder: &mut dyn ArrayBuilder,
    ) -> Result<()> {
        if col_idx >= self.col_offsets.len() {
            return Err(DataFusionError::Execution(format!(
                "V5: column index {col_idx} out of range (have {})",
                self.col_offsets.len()
            )));
        }
        let (col_offset, type_code) = self.col_offsets[col_idx];
        let dt = decode_data_type(type_code)?;
        let data = &self.data[col_offset..];
        let n = self.num_alleles;

        match dt {
            DataType::Int32 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<Int32Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected Int32Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 4;
                        let v = i32::from_le_bytes(data[off..off + 4].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::Int64 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<Int64Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected Int64Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 8;
                        let v = i64::from_le_bytes(data[off..off + 8].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::Float32 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<Float32Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected Float32Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 4;
                        let v = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::Float64 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<Float64Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected Float64Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 8;
                        let v = f64::from_le_bytes(data[off..off + 8].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::UInt32 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<UInt32Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected UInt32Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 4;
                        let v = u32::from_le_bytes(data[off..off + 4].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::UInt64 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<UInt64Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected UInt64Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 8;
                        let v = u64::from_le_bytes(data[off..off + 8].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::Int8 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<Int8Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected Int8Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        b.append_value(data[values_off + row] as i8);
                    }
                }
            }
            DataType::Int16 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<Int16Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected Int16Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 2;
                        let v = i16::from_le_bytes(data[off..off + 2].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::UInt8 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<datafusion::arrow::array::UInt8Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected UInt8Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        b.append_value(data[values_off + row]);
                    }
                }
            }
            DataType::UInt16 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<UInt16Builder>()
                    .ok_or_else(|| DataFusionError::Execution("expected UInt16Builder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let off = values_off + row * 2;
                        let v = u16::from_le_bytes(data[off..off + 2].try_into().unwrap());
                        b.append_value(v);
                    }
                }
            }
            DataType::Boolean => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<BooleanBuilder>()
                    .ok_or_else(|| DataFusionError::Execution("expected BooleanBuilder".into()))?;
                let values_off = null_bitmap_size(n);
                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let bit = (data[values_off + row / 8] >> (row % 8)) & 1;
                        b.append_value(bit != 0);
                    }
                }
            }
            DataType::Utf8 | DataType::LargeUtf8 => {
                let b = builder
                    .as_any_mut()
                    .downcast_mut::<StringBuilder>()
                    .ok_or_else(|| DataFusionError::Execution("expected StringBuilder".into()))?;
                let str_data_off = null_bitmap_size(n);
                // String layout: [4B total_bytes][N*4B cumulative offsets][string data]
                // Skip total_str_bytes (4 bytes) — we use per-row offsets directly.
                let offsets_start = str_data_off + 4;
                let string_data_start = offsets_start + n * 4;

                for &row in allele_rows {
                    if !is_non_null(data, row) {
                        b.append_null();
                    } else {
                        let start = if row == 0 {
                            0
                        } else {
                            let off = offsets_start + (row - 1) * 4;
                            u32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as usize
                        };
                        let end_off = offsets_start + row * 4;
                        let end = u32::from_le_bytes(data[end_off..end_off + 4].try_into().unwrap())
                            as usize;
                        let s = std::str::from_utf8(
                            &data[string_data_start + start..string_data_start + end],
                        )
                        .unwrap_or("");
                        b.append_value(s);
                    }
                }
            }
            other => {
                return Err(DataFusionError::Execution(format!(
                    "V5 read: unsupported data type {other}"
                )));
            }
        }
        Ok(())
    }

}

/// Compute the packed size of a single column's data (after the type byte).
fn column_packed_size(dt: &DataType, n: usize, null_bytes: usize, data: &[u8]) -> Result<usize> {
    let fixed = match dt {
        DataType::Int32 | DataType::Float32 | DataType::UInt32 => null_bytes + n * 4,
        DataType::Int64 | DataType::Float64 | DataType::UInt64 => null_bytes + n * 8,
        DataType::Int8 | DataType::UInt8 => null_bytes + n,
        DataType::Int16 | DataType::UInt16 => null_bytes + n * 2,
        DataType::Boolean => null_bytes + n.div_ceil(8),
        DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View => {
            // null_bitmap + 4B total_string_bytes + N*4B offsets + string data
            if data.len() < null_bytes + 4 {
                return Err(DataFusionError::Execution(
                    "V5: truncated string column header".into(),
                ));
            }
            let total_str =
                u32::from_le_bytes(data[null_bytes..null_bytes + 4].try_into().unwrap()) as usize;
            null_bytes + 4 + n * 4 + total_str
        }
        other => {
            return Err(DataFusionError::Execution(format!(
                "V5 column_packed_size: unsupported type {other}"
            )));
        }
    };
    Ok(fixed)
}

/// Check if a row is non-null in the null bitmap at the start of `data`.
/// Returns `true` if the bit is set (non-null).
#[inline(always)]
fn is_non_null(data: &[u8], row: usize) -> bool {
    (data[row / 8] >> (row % 8)) & 1 != 0
}

/// Byte offset where values start after the null bitmap.
#[inline(always)]
fn null_bitmap_size(n: usize) -> usize {
    n.div_ceil(8)
}

/// Extract a string value from Utf8, Utf8View, or LargeUtf8 array.
fn get_string_value(col: &dyn Array, i: usize) -> &str {
    if col.is_null(i) {
        return "";
    }
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        arr.value(i)
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        arr.value(i)
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::LargeStringArray>()
    {
        arr.value(i)
    } else {
        ""
    }
}

/// Create an appropriate ArrayBuilder for a DataType.
pub fn make_builder(dt: &DataType, capacity: usize) -> Result<Box<dyn ArrayBuilder>> {
    match dt {
        DataType::Int32 => Ok(Box::new(Int32Builder::with_capacity(capacity))),
        DataType::Int64 => Ok(Box::new(Int64Builder::with_capacity(capacity))),
        DataType::Float32 => Ok(Box::new(Float32Builder::with_capacity(capacity))),
        DataType::Float64 => Ok(Box::new(Float64Builder::with_capacity(capacity))),
        DataType::UInt32 => Ok(Box::new(UInt32Builder::with_capacity(capacity))),
        DataType::UInt64 => Ok(Box::new(UInt64Builder::with_capacity(capacity))),
        DataType::Int8 => Ok(Box::new(Int8Builder::with_capacity(capacity))),
        DataType::Int16 => Ok(Box::new(Int16Builder::with_capacity(capacity))),
        DataType::UInt8 => Ok(Box::new(
            datafusion::arrow::array::UInt8Builder::with_capacity(capacity),
        )),
        DataType::UInt16 => Ok(Box::new(UInt16Builder::with_capacity(capacity))),
        DataType::Boolean => Ok(Box::new(BooleanBuilder::with_capacity(capacity))),
        // For output, we always use Utf8 (normalized from Utf8View/LargeUtf8).
        DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => Ok(Box::new(
            StringBuilder::with_capacity(capacity, capacity * 16),
        )),
        other => Err(DataFusionError::Execution(format!(
            "V5 make_builder: unsupported type {other}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::StringArray;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn make_test_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("score", DataType::Int32, true),
            Field::new("flag", DataType::Boolean, true),
            Field::new("ratio", DataType::Float64, true),
        ]));
        RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(StringArray::from(vec!["A/G", "A/T", "C/G"])),
                Arc::new(StringArray::from(vec![Some("rs123"), None, Some("rs789")])),
                Arc::new(Int32Array::from(vec![Some(42), Some(99), None])),
                Arc::new(BooleanArray::from(vec![Some(true), Some(false), None])),
                Arc::new(Float64Array::from(vec![Some(0.5), None, Some(1.0)])),
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_position_entry_roundtrip() {
        let batch = make_test_batch();
        let rows = vec![0, 1, 2];
        let col_indices = vec![4, 5, 6, 7]; // variation_name, score, flag, ratio
        let allele_col_idx = 3; // allele_string

        let bytes = serialize_position_entry(&rows, &batch, &col_indices, allele_col_idx).unwrap();

        let reader = PositionEntryReader::new(&bytes).unwrap();
        assert_eq!(reader.num_alleles(), 3);
        assert_eq!(reader.num_cols(), 4);
        assert_eq!(reader.allele_string(0), "A/G");
        assert_eq!(reader.allele_string(1), "A/T");
        assert_eq!(reader.allele_string(2), "C/G");

        // Read variation_name (Utf8, col 0) for all alleles.
        let mut str_builder = StringBuilder::new();
        reader
            .append_column_values(0, &[0, 1, 2], &mut str_builder)
            .unwrap();
        let arr = str_builder.finish();
        assert_eq!(arr.value(0), "rs123");
        assert!(arr.is_null(1));
        assert_eq!(arr.value(2), "rs789");

        // Read score (Int32, col 1).
        let mut int_builder = Int32Builder::new();
        reader
            .append_column_values(1, &[0, 1, 2], &mut int_builder)
            .unwrap();
        let arr = int_builder.finish();
        assert_eq!(arr.value(0), 42);
        assert_eq!(arr.value(1), 99);
        assert!(arr.is_null(2));

        // Read flag (Boolean, col 2).
        let mut bool_builder = BooleanBuilder::new();
        reader
            .append_column_values(2, &[0, 1, 2], &mut bool_builder)
            .unwrap();
        let arr = bool_builder.finish();
        assert!(arr.value(0));
        assert!(!arr.value(1));
        assert!(arr.is_null(2));

        // Read ratio (Float64, col 3).
        let mut f64_builder = Float64Builder::new();
        reader
            .append_column_values(3, &[0, 1, 2], &mut f64_builder)
            .unwrap();
        let arr = f64_builder.finish();
        assert!((arr.value(0) - 0.5).abs() < f64::EPSILON);
        assert!(arr.is_null(1));
        assert!((arr.value(2) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_position_entry_partial_rows() {
        let batch = make_test_batch();
        let rows = vec![0, 2]; // Only first and third alleles
        let col_indices = vec![4]; // variation_name only
        let allele_col_idx = 3;

        let bytes = serialize_position_entry(&rows, &batch, &col_indices, allele_col_idx).unwrap();

        let reader = PositionEntryReader::new(&bytes).unwrap();
        assert_eq!(reader.num_alleles(), 2);
        assert_eq!(reader.allele_string(0), "A/G");
        assert_eq!(reader.allele_string(1), "C/G");

        // Read only allele row 1 from the entry.
        let mut str_builder = StringBuilder::new();
        reader
            .append_column_values(0, &[1], &mut str_builder)
            .unwrap();
        let arr = str_builder.finish();
        assert_eq!(arr.value(0), "rs789");
    }

    #[test]
    fn test_position_entry_empty() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("val", DataType::Int64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(Int64Array::from(vec![Some(42)])),
            ],
        )
        .unwrap();

        let bytes = serialize_position_entry(&[0], &batch, &[1], 0).unwrap();
        let reader = PositionEntryReader::new(&bytes).unwrap();
        assert_eq!(reader.num_alleles(), 1);
        assert_eq!(reader.num_cols(), 1);
        assert_eq!(reader.allele_string(0), "A/G");

        let mut builder = Int64Builder::new();
        reader.append_column_values(0, &[0], &mut builder).unwrap();
        let arr = builder.finish();
        assert_eq!(arr.value(0), 42);
    }
}
