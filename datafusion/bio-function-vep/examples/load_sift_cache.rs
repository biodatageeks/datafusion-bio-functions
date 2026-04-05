//! Convert translation_sift parquet → fjall sift keyspace.
//!
//! Usage:
//!   cargo run --release --example load_sift_cache --features kv-cache -- \
//!     <translation_sift_parquet> <fjall_db_path>
//!
//! The sift keyspace is added to an existing fjall database (created by load_cache).
//! If the database doesn't exist, it is created.

use std::time::Instant;

use datafusion::arrow::array::{
    Array, Float32Array, Int32Array, ListArray, StringArray, StringViewArray,
};
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::kv_cache::sift_store::SiftKvStore;
use datafusion_bio_function_vep::transcript_consequence::{CachedPredictions, CompactPrediction};

fn read_compact_predictions_from_batch(col: &dyn Array, row: usize) -> Vec<CompactPrediction> {
    if col.is_null(row) {
        return Vec::new();
    }
    let Some(list_arr) = col.as_any().downcast_ref::<ListArray>() else {
        return Vec::new();
    };
    let offsets = list_arr.offsets();
    let start_off = offsets[row] as usize;
    let end_off = offsets[row + 1] as usize;
    if start_off == end_off {
        return Vec::new();
    }
    let values = list_arr.values();
    let struct_arr = values
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StructArray>()
        .unwrap();
    let Some(positions) = struct_arr.column_by_name("position") else {
        return Vec::new();
    };
    let Some(amino_acids) = struct_arr.column_by_name("amino_acid") else {
        return Vec::new();
    };
    let Some(predictions) = struct_arr.column_by_name("prediction") else {
        return Vec::new();
    };
    let Some(scores) = struct_arr.column_by_name("score") else {
        return Vec::new();
    };
    let pos_arr = positions.as_any().downcast_ref::<Int32Array>();
    let aa_arr = amino_acids.as_any().downcast_ref::<StringArray>();
    let aa_view = amino_acids.as_any().downcast_ref::<StringViewArray>();
    let pred_arr = predictions.as_any().downcast_ref::<StringArray>();
    let pred_view = predictions.as_any().downcast_ref::<StringViewArray>();
    let score_arr = scores.as_any().downcast_ref::<Float32Array>();

    let mut out = Vec::with_capacity(end_off - start_off);
    for i in start_off..end_off {
        let Some(pos) = pos_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
        else {
            continue;
        };
        let aa_str = aa_arr
            .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
            .or_else(|| aa_view.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }));
        let pred_str = pred_arr
            .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
            .or_else(|| pred_view.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }));
        let score = score_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
        if let (Some(aa), Some(pred), Some(sc)) = (aa_str, pred_str, score)
            && let Some(aa_idx) = CompactPrediction::encode_amino_acid(aa)
        {
            out.push(CompactPrediction {
                position: pos,
                amino_acid: aa_idx,
                prediction: CompactPrediction::encode_prediction(pred),
                score: sc,
            });
        }
    }
    out
}

fn string_at(array: &dyn Array, row: usize) -> Option<&str> {
    if array.is_null(row) {
        return None;
    }
    if let Some(a) = array.as_any().downcast_ref::<StringArray>() {
        return Some(a.value(row));
    }
    if let Some(a) = array.as_any().downcast_ref::<StringViewArray>() {
        return Some(a.value(row));
    }
    None
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <translation_sift_parquet> <fjall_db_path>",
            args[0]
        );
        std::process::exit(1);
    }

    let parquet_path = &args[1];
    let fjall_path = &args[2];

    eprintln!("Loading sift predictions from: {parquet_path}");
    eprintln!("Into fjall database: {fjall_path}");

    let t_start = Instant::now();

    // Open or create the fjall database.
    let db = fjall::Database::builder(fjall_path)
        .cache_size(256 * 1024 * 1024)
        .manual_journal_persist(true)
        .open()
        .map_err(|e| DataFusionError::External(Box::new(e)))?;

    let sift_store = SiftKvStore::create(&db)?;

    // Read parquet with parquet-rs directly (not DataFusion).
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(parquet_path)
        .map_err(|e| DataFusionError::Execution(format!("cannot open {parquet_path}: {e}")))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| DataFusionError::Execution(format!("parquet error: {e}")))?
        .build()
        .map_err(|e| DataFusionError::Execution(format!("parquet error: {e}")))?;

    let mut total_transcripts = 0u64;
    let mut total_sift_entries = 0u64;
    let mut total_polyphen_entries = 0u64;

    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| DataFusionError::Execution(format!("batch read error: {e}")))?;
        let schema = batch.schema();
        let tid_idx = schema.index_of("transcript_id").ok();
        let sift_idx = schema.index_of("sift_predictions").ok();
        let poly_idx = schema.index_of("polyphen_predictions").ok();

        let Some(tid_idx) = tid_idx else { continue };

        for row in 0..batch.num_rows() {
            let Some(transcript_id) = string_at(batch.column(tid_idx).as_ref(), row) else {
                continue;
            };

            let mut preds = CachedPredictions::default();
            if let Some(idx) = sift_idx {
                preds.sift = read_compact_predictions_from_batch(batch.column(idx).as_ref(), row);
            }
            if let Some(idx) = poly_idx {
                preds.polyphen =
                    read_compact_predictions_from_batch(batch.column(idx).as_ref(), row);
            }

            if preds.sift.is_empty() && preds.polyphen.is_empty() {
                continue;
            }

            preds.sort();
            total_sift_entries += preds.sift.len() as u64;
            total_polyphen_entries += preds.polyphen.len() as u64;
            total_transcripts += 1;

            sift_store.put(transcript_id, &preds)?;
        }
    }

    db.persist(fjall::PersistMode::SyncAll)
        .map_err(|e| DataFusionError::External(Box::new(e)))?;

    // Major-compact the sift keyspace so the LSM tree is fully optimized
    // for reads (bloom filters built, levels merged).
    // Without this, every lookup probes hundreds of L0 SSTs.
    eprintln!("Running major compaction on sift keyspace...");
    let compact_start = Instant::now();
    sift_store
        .keyspace()
        .major_compact()
        .map_err(|e| DataFusionError::External(Box::new(e)))?;
    eprintln!(
        "Major compaction completed in {:.1}s",
        compact_start.elapsed().as_secs_f64()
    );

    // Skip fjall's Drop which deadlocks when background worker threads are
    // busy compacting (issue #86). All data is persisted and compacted above.
    std::mem::forget(sift_store);
    std::mem::forget(db);

    let elapsed = t_start.elapsed().as_secs_f64();
    eprintln!(
        "Done in {elapsed:.1}s: {total_transcripts} transcripts, {total_sift_entries} sift entries, {total_polyphen_entries} polyphen entries"
    );

    Ok(())
}
