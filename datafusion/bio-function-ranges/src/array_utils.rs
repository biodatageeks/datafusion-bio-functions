use datafusion::arrow::array::{
    GenericStringArray, Int32Array, Int64Array, RecordBatch, StringViewArray, UInt32Array,
    UInt64Array,
};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::{DataFusionError, Result};

pub enum ContigArray<'a> {
    GenericString(&'a GenericStringArray<i64>),
    Utf8View(&'a StringViewArray),
    Utf8(&'a GenericStringArray<i32>),
}

impl ContigArray<'_> {
    pub fn value(&self, i: usize) -> &str {
        match self {
            ContigArray::GenericString(arr) => arr.value(i),
            ContigArray::Utf8View(arr) => arr.value(i),
            ContigArray::Utf8(arr) => arr.value(i),
        }
    }
}

pub enum PosArray<'a> {
    Int32(&'a Int32Array),
    Int64(&'a Int64Array),
    UInt32(&'a UInt32Array),
    UInt64(&'a UInt64Array),
}

impl PosArray<'_> {
    pub fn value(&self, i: usize) -> Result<i32> {
        match self {
            PosArray::Int32(arr) => Ok(arr.value(i)),
            PosArray::Int64(arr) => {
                let v = arr.value(i);
                i32::try_from(v).map_err(|_| {
                    DataFusionError::Execution(format!(
                        "coordinate value {v} at row {i} overflows i32 (max {})",
                        i32::MAX
                    ))
                })
            }
            PosArray::UInt32(arr) => {
                let v = arr.value(i);
                i32::try_from(v).map_err(|_| {
                    DataFusionError::Execution(format!(
                        "coordinate value {v} at row {i} overflows i32 (max {})",
                        i32::MAX
                    ))
                })
            }
            PosArray::UInt64(arr) => {
                let v = arr.value(i);
                i32::try_from(v).map_err(|_| {
                    DataFusionError::Execution(format!(
                        "coordinate value {v} at row {i} overflows i32 (max {})",
                        i32::MAX
                    ))
                })
            }
        }
    }
}

/// Extract contig, start, and end column arrays from a [`RecordBatch`].
///
/// Returns an error if a column is missing or has an unsupported data type.
pub fn get_join_col_arrays<'a>(
    batch: &'a RecordBatch,
    columns: (&str, &str, &str),
) -> Result<(ContigArray<'a>, PosArray<'a>, PosArray<'a>)> {
    let contig_col = batch.column_by_name(columns.0).ok_or_else(|| {
        DataFusionError::Plan(format!(
            "contig column '{}' not found in batch with columns: {:?}",
            columns.0,
            batch
                .schema()
                .fields()
                .iter()
                .map(|f| f.name())
                .collect::<Vec<_>>()
        ))
    })?;

    let contig_any = contig_col.as_any();
    let contig_arr = match contig_col.data_type() {
        DataType::LargeUtf8 => {
            let arr = contig_any
                .downcast_ref::<GenericStringArray<i64>>()
                .ok_or_else(|| {
                    DataFusionError::Internal(format!(
                        "failed to downcast contig column '{}' to LargeUtf8",
                        columns.0
                    ))
                })?;
            ContigArray::GenericString(arr)
        }
        DataType::Utf8View => {
            let arr = contig_any
                .downcast_ref::<StringViewArray>()
                .ok_or_else(|| {
                    DataFusionError::Internal(format!(
                        "failed to downcast contig column '{}' to Utf8View",
                        columns.0
                    ))
                })?;
            ContigArray::Utf8View(arr)
        }
        DataType::Utf8 => {
            let arr = contig_any
                .downcast_ref::<GenericStringArray<i32>>()
                .ok_or_else(|| {
                    DataFusionError::Internal(format!(
                        "failed to downcast contig column '{}' to Utf8",
                        columns.0
                    ))
                })?;
            ContigArray::Utf8(arr)
        }
        dt => {
            return Err(DataFusionError::NotImplemented(format!(
                "unsupported data type {dt:?} for contig column '{}'; expected Utf8, LargeUtf8, or Utf8View",
                columns.0
            )));
        }
    };

    let start_arr = extract_pos_array(batch, columns.1, "start")?;
    let end_arr = extract_pos_array(batch, columns.2, "end")?;

    Ok((contig_arr, start_arr, end_arr))
}

/// Extract a position (start/end) column as a [`PosArray`].
fn extract_pos_array<'a>(
    batch: &'a RecordBatch,
    col_name: &str,
    label: &str,
) -> Result<PosArray<'a>> {
    let col = batch.column_by_name(col_name).ok_or_else(|| {
        DataFusionError::Plan(format!(
            "{label} column '{col_name}' not found in batch with columns: {:?}",
            batch
                .schema()
                .fields()
                .iter()
                .map(|f| f.name())
                .collect::<Vec<_>>()
        ))
    })?;

    let any = col.as_any();
    match col.data_type() {
        DataType::Int32 => Ok(PosArray::Int32(
            any.downcast_ref::<Int32Array>().ok_or_else(|| {
                DataFusionError::Internal(format!(
                    "failed to downcast {label} column '{col_name}' to Int32"
                ))
            })?,
        )),
        DataType::Int64 => Ok(PosArray::Int64(
            any.downcast_ref::<Int64Array>().ok_or_else(|| {
                DataFusionError::Internal(format!(
                    "failed to downcast {label} column '{col_name}' to Int64"
                ))
            })?,
        )),
        DataType::UInt32 => Ok(PosArray::UInt32(
            any.downcast_ref::<UInt32Array>().ok_or_else(|| {
                DataFusionError::Internal(format!(
                    "failed to downcast {label} column '{col_name}' to UInt32"
                ))
            })?,
        )),
        DataType::UInt64 => Ok(PosArray::UInt64(
            any.downcast_ref::<UInt64Array>().ok_or_else(|| {
                DataFusionError::Internal(format!(
                    "failed to downcast {label} column '{col_name}' to UInt64"
                ))
            })?,
        )),
        dt => Err(DataFusionError::NotImplemented(format!(
            "unsupported data type {dt:?} for {label} column '{col_name}'; expected Int32, Int64, UInt32, or UInt64"
        ))),
    }
}
