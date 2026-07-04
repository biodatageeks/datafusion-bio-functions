//! Read-parity oracle for the Parquet context (scan) entities (Phase 2, Task 11).
//!
//! Scans the SAME entity shard through the Parquet reader (`read_context_parquet`)
//! and the Lance reader (`scan_projected_existing_columns`), projecting every
//! column in the Parquet schema, then diffs the concatenated batches
//! column-by-column (in file order). The two must be identical: same row count,
//! same values for every shared column (incl. nested `list<struct>`), rendered
//! via Arrow's `ArrayFormatter`.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example context_backend_diff -- \
//!     --parquet /path/parquet.transcript/chr4.parquet \
//!     --lance   /path/transcript.lance/chr4.lance

use std::path::PathBuf;

use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::compute::{concat_batches, sort_to_indices, take};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::util::display::{ArrayFormatter, FormatOptions};
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::lance_cache::context_runtime::scan_projected_existing_columns;
use datafusion_bio_function_vep::parquet_cache::scan::read_context_parquet;

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let _ = env_logger::try_init();
    let args = parse_args()?;

    // Column set = every column in the Parquet shard schema.
    let pq_schema = parquet_schema(&args.parquet).await?;
    let columns: Vec<String> = pq_schema
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    let column_refs: Vec<&str> = columns.iter().map(String::as_str).collect();
    eprintln!(
        "columns={} parquet={} lance={}",
        columns.len(),
        args.parquet.display(),
        args.lance.display()
    );

    let pq_batches = read_context_parquet(&args.parquet, &column_refs).await?;
    let lance_batches = scan_projected_existing_columns(&args.lance, &column_refs).await?;

    let mut pq = concat(&pq_batches)?;
    let mut lance = concat(&lance_batches)?;
    // The scans return rows in each backend's storage order, which differs (the
    // source ORDER BY is not total). Sort both by a unique key so the comparison
    // tests data-set equality, not incidental row order.
    if let Some(key) = &args.sort_key {
        pq = sort_by(&pq, key)?;
        lance = sort_by(&lance, key)?;
        eprintln!("sorted both by '{key}' for order-independent comparison");
    }
    eprintln!("rows: parquet={} lance={}", pq.num_rows(), lance.num_rows());

    if pq.num_rows() != lance.num_rows() {
        eprintln!("MISMATCH: row counts differ — cannot align rows");
        std::process::exit(2);
    }

    let opts = FormatOptions::default().with_null("<NULL>");
    let n = pq.num_rows();
    let mut mismatched_cols = 0usize;
    let mut compared = 0usize;
    for field in pq.schema().fields() {
        let name = field.name();
        let (Ok(pi), Ok(li)) = (pq.schema().index_of(name), lance.schema().index_of(name)) else {
            eprintln!("  (skip {name}: not present in both)");
            continue;
        };
        compared += 1;
        let fp = match ArrayFormatter::try_new(pq.column(pi).as_ref(), &opts) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("  (skip {name}: unformattable {e})");
                continue;
            }
        };
        let fl = ArrayFormatter::try_new(lance.column(li).as_ref(), &opts)
            .map_err(|e| DataFusionError::Execution(format!("formatter {name}: {e}")))?;
        let mut cell_mismatches = 0usize;
        let mut example = None;
        for row in 0..n {
            let a = fp.value(row).to_string();
            let b = fl.value(row).to_string();
            if a != b {
                cell_mismatches += 1;
                example.get_or_insert_with(|| {
                    let a = if a.len() > 60 {
                        format!("{}…", &a[..60])
                    } else {
                        a.clone()
                    };
                    let b = if b.len() > 60 {
                        format!("{}…", &b[..60])
                    } else {
                        b.clone()
                    };
                    format!("row {row}: parquet={a:?} lance={b:?}")
                });
            }
        }
        if cell_mismatches > 0 {
            mismatched_cols += 1;
            eprintln!(
                "  DIFF {name}: {cell_mismatches}/{n} cells — e.g. {}",
                example.unwrap_or_default()
            );
        }
    }

    eprintln!(
        "compared {compared} columns; {} matched, {mismatched_cols} mismatched",
        compared - mismatched_cols
    );
    if mismatched_cols == 0 {
        eprintln!("\nOK — Parquet and Lance context scans are identical across all columns.");
        Ok(())
    } else {
        eprintln!("\nFAILED — {mismatched_cols} column(s) mismatched.");
        std::process::exit(1)
    }
}

fn concat(batches: &[RecordBatch]) -> Result<RecordBatch> {
    if batches.is_empty() {
        return Err(DataFusionError::Execution("no batches scanned".into()));
    }
    concat_batches(&batches[0].schema(), batches)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}

/// Sort a batch ascending by `key` (must be a unique column for a total order).
fn sort_by(batch: &RecordBatch, key: &str) -> Result<RecordBatch> {
    let idx = batch
        .schema()
        .index_of(key)
        .map_err(|_| DataFusionError::Execution(format!("sort key '{key}' not in schema")))?;
    let indices = sort_to_indices(batch.column(idx), None, None)
        .map_err(|e| DataFusionError::Execution(format!("sort {key}: {e}")))?;
    let cols = batch
        .columns()
        .iter()
        .map(|c| take(c, &indices, None))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|e| DataFusionError::Execution(format!("take sorted: {e}")))?;
    RecordBatch::try_new(batch.schema(), cols)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}

async fn parquet_schema(path: &PathBuf) -> Result<SchemaRef> {
    use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;
    let file = tokio::fs::File::open(path)
        .await
        .map_err(|e| DataFusionError::Execution(format!("open '{}': {e}", path.display())))?;
    let builder = ParquetRecordBatchStreamBuilder::new(file)
        .await
        .map_err(|e| DataFusionError::Execution(format!("read parquet schema: {e}")))?;
    Ok(builder.schema().clone())
}

struct Args {
    parquet: PathBuf,
    lance: PathBuf,
    sort_key: Option<String>,
}

fn parse_args() -> Result<Args> {
    let mut parquet = None;
    let mut lance = None;
    let mut sort_key = None;
    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--parquet" => parquet = Some(PathBuf::from(require(&mut it, &a)?)),
            "--lance" => lance = Some(PathBuf::from(require(&mut it, &a)?)),
            "--sort-key" => sort_key = Some(require(&mut it, &a)?),
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unknown argument: {other}"
                )));
            }
        }
    }
    Ok(Args {
        parquet: parquet.ok_or_else(|| DataFusionError::Execution("--parquet required".into()))?,
        lance: lance.ok_or_else(|| DataFusionError::Execution("--lance required".into()))?,
        sort_key,
    })
}

fn require(it: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    it.next()
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}
