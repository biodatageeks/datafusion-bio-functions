use datafusion::arrow::array::{
    GenericStringArray, Int32Array, Int64Array, RecordBatch, StringViewArray, UInt32Array,
    UInt64Array,
};
use datafusion::arrow::datatypes::DataType;

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
    pub fn value(&self, i: usize) -> i32 {
        match self {
            PosArray::Int32(arr) => arr.value(i),
            PosArray::Int64(arr) => arr.value(i) as i32,
            PosArray::UInt32(arr) => arr.value(i) as i32,
            PosArray::UInt64(arr) => arr.value(i) as i32,
        }
    }
}

pub fn get_join_col_arrays(
    batch: &RecordBatch,
    columns: (String, String, String),
) -> (ContigArray<'_>, PosArray<'_>, PosArray<'_>) {
    let contig_arr = match batch.column_by_name(&columns.0).unwrap().data_type() {
        DataType::LargeUtf8 => {
            let contig_arr = batch
                .column_by_name(&columns.0)
                .unwrap()
                .as_any()
                .downcast_ref::<GenericStringArray<i64>>()
                .unwrap();
            ContigArray::GenericString(contig_arr)
        }
        DataType::Utf8View => {
            let contig_arr = batch
                .column_by_name(&columns.0)
                .unwrap()
                .as_any()
                .downcast_ref::<StringViewArray>()
                .unwrap();
            ContigArray::Utf8View(contig_arr)
        }
        DataType::Utf8 => {
            let contig_arr = batch
                .column_by_name(&columns.0)
                .unwrap()
                .as_any()
                .downcast_ref::<GenericStringArray<i32>>()
                .unwrap();
            ContigArray::Utf8(contig_arr)
        }
        _ => todo!(),
    };

    let start_arr = match batch.column_by_name(&columns.1).unwrap().data_type() {
        DataType::Int32 => {
            let start_arr = batch
                .column_by_name(&columns.1)
                .unwrap()
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap();
            PosArray::Int32(start_arr)
        }
        DataType::Int64 => {
            let start_arr = batch
                .column_by_name(&columns.1)
                .unwrap()
                .as_any()
                .downcast_ref::<Int64Array>()
                .unwrap();
            PosArray::Int64(start_arr)
        }
        DataType::UInt32 => {
            let start_arr = batch
                .column_by_name(&columns.1)
                .unwrap()
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();
            PosArray::UInt32(start_arr)
        }
        DataType::UInt64 => {
            let start_arr = batch
                .column_by_name(&columns.1)
                .unwrap()
                .as_any()
                .downcast_ref::<UInt64Array>()
                .unwrap();
            PosArray::UInt64(start_arr)
        }
        dt => panic!("Unsupported data type for start column: {dt:?}"),
    };

    let end_arr = match batch.column_by_name(&columns.2).unwrap().data_type() {
        DataType::Int32 => {
            let end_arr = batch
                .column_by_name(&columns.2)
                .unwrap()
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap();
            PosArray::Int32(end_arr)
        }
        DataType::Int64 => {
            let end_arr = batch
                .column_by_name(&columns.2)
                .unwrap()
                .as_any()
                .downcast_ref::<Int64Array>()
                .unwrap();
            PosArray::Int64(end_arr)
        }
        DataType::UInt32 => {
            let end_arr = batch
                .column_by_name(&columns.2)
                .unwrap()
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();
            PosArray::UInt32(end_arr)
        }
        DataType::UInt64 => {
            let end_arr = batch
                .column_by_name(&columns.2)
                .unwrap()
                .as_any()
                .downcast_ref::<UInt64Array>()
                .unwrap();
            PosArray::UInt64(end_arr)
        }
        dt => panic!("Unsupported data type for end column: {dt:?}"),
    };

    (contig_arr, start_arr, end_arr)
}
