//! Read-parity oracle for the Parquet translation_sift backend (Phase 2, Task 11).
//!
//! Runs the SAME u64 keys through the Parquet point-lookup
//! (`SinglePathParquetSiftLookup`) and the Lance point-lookup
//! (`KeyU64LanceLookup`) against real per-chromosome shards, then compares the
//! returned `key -> (sift, poly)` blob maps. The two must be identical: same
//! matched keys, byte-identical SIFT/PolyPhen blobs.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example sift_backend_diff -- \
//!     --parquet /path/parquet.translation_sift/chr4.parquet \
//!     --lance   /path/translation_sift.lance/chr4.lance \
//!     --keys /path/sift_keys.txt

use std::collections::BTreeMap;
use std::path::PathBuf;

use datafusion::arrow::array::{Array, BinaryArray, RecordBatch, UInt64Array};
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::lance_cache::context_runtime::KeyU64LanceLookup;
use datafusion_bio_function_vep::parquet_cache::sift::SinglePathParquetSiftLookup;

const PROJECTION: [&str; 3] = ["key", "sift", "poly"];

/// `key -> (sift_blob, poly_blob)`.
type BlobMap = BTreeMap<u64, (Vec<u8>, Vec<u8>)>;

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let _ = env_logger::try_init();
    let args = parse_args()?;
    let keys = read_keys(&args.keys)?;
    eprintln!(
        "keys={} parquet={} lance={}",
        keys.len(),
        args.parquet.display(),
        args.lance.display()
    );

    let pq = SinglePathParquetSiftLookup::open(&args.parquet, projection()).await?;
    let (pq_batch, pq_present) = pq.take_keys(&keys).await?;

    let lance = KeyU64LanceLookup::open(&args.lance, projection()).await?;
    let (lance_batch, lance_present) = lance.take_keys(&keys).await?;

    eprintln!(
        "present: parquet={} lance={}   rows: parquet={} lance={}",
        pq_present.len(),
        lance_present.len(),
        pq_batch.num_rows(),
        lance_batch.num_rows()
    );

    let pq_map = key_blob_map(&pq_batch)?;
    let lance_map = key_blob_map(&lance_batch)?;

    let mut failures = 0usize;
    if pq_map.len() != lance_map.len() {
        eprintln!(
            "MISMATCH: distinct-key counts differ (parquet={} lance={})",
            pq_map.len(),
            lance_map.len()
        );
        failures += 1;
    }
    let mut sift_mismatch = 0usize;
    let mut poly_mismatch = 0usize;
    let mut missing = 0usize;
    for (k, (l_sift, l_poly)) in &lance_map {
        match pq_map.get(k) {
            None => {
                missing += 1;
                if missing <= 5 {
                    eprintln!("  key {k} present in Lance but missing in Parquet");
                }
            }
            Some((p_sift, p_poly)) => {
                if p_sift != l_sift {
                    sift_mismatch += 1;
                    if sift_mismatch <= 5 {
                        eprintln!("  key {k} SIFT blob differs");
                    }
                }
                if p_poly != l_poly {
                    poly_mismatch += 1;
                    if poly_mismatch <= 5 {
                        eprintln!("  key {k} PolyPhen blob differs");
                    }
                }
            }
        }
    }
    // Keys present in Parquet but not Lance.
    let extra = pq_map.keys().filter(|k| !lance_map.contains_key(k)).count();

    eprintln!(
        "compared {} keys: sift_mismatch={sift_mismatch} poly_mismatch={poly_mismatch} missing_in_parquet={missing} extra_in_parquet={extra}",
        lance_map.len()
    );
    failures += sift_mismatch + poly_mismatch + missing + extra;

    if failures == 0 {
        eprintln!(
            "\nOK — Parquet and Lance translation_sift reads are byte-identical (key + sift + poly)."
        );
        Ok(())
    } else {
        eprintln!("\nFAILED — {failures} discrepancies.");
        std::process::exit(1)
    }
}

fn projection() -> Vec<String> {
    PROJECTION.iter().map(|s| s.to_string()).collect()
}

/// `key -> (sift_blob, poly_blob)` for every row in a sift take batch.
fn key_blob_map(batch: &RecordBatch) -> Result<BlobMap> {
    let key = batch
        .column_by_name("key")
        .and_then(|c| c.as_any().downcast_ref::<UInt64Array>())
        .ok_or_else(|| DataFusionError::Execution("batch missing UInt64 key".into()))?;
    let sift = batch
        .column_by_name("sift")
        .and_then(|c| c.as_any().downcast_ref::<BinaryArray>());
    let poly = batch
        .column_by_name("poly")
        .and_then(|c| c.as_any().downcast_ref::<BinaryArray>());
    let blob = |arr: Option<&BinaryArray>, i: usize| -> Vec<u8> {
        match arr {
            Some(a) if a.is_valid(i) => a.value(i).to_vec(),
            _ => Vec::new(),
        }
    };
    let mut m = BTreeMap::new();
    for i in 0..batch.num_rows() {
        if key.is_null(i) {
            continue;
        }
        m.insert(key.value(i), (blob(sift, i), blob(poly, i)));
    }
    Ok(m)
}

struct Args {
    parquet: PathBuf,
    lance: PathBuf,
    keys: PathBuf,
}

fn parse_args() -> Result<Args> {
    let mut parquet = None;
    let mut lance = None;
    let mut keys = None;
    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--parquet" => parquet = Some(PathBuf::from(require(&mut it, &a)?)),
            "--lance" => lance = Some(PathBuf::from(require(&mut it, &a)?)),
            "--keys" => keys = Some(PathBuf::from(require(&mut it, &a)?)),
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
        keys: keys.ok_or_else(|| DataFusionError::Execution("--keys required".into()))?,
    })
}

fn require(it: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    it.next()
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}

fn read_keys(path: &PathBuf) -> Result<Vec<u64>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| DataFusionError::Execution(format!("read keys '{}': {e}", path.display())))?;
    let mut v: Vec<u64> = text
        .split_whitespace()
        .map(|t| {
            t.parse::<u64>()
                .map_err(|e| DataFusionError::Execution(format!("bad key '{t}': {e}")))
        })
        .collect::<Result<Vec<_>>>()?;
    v.sort_unstable();
    v.dedup();
    Ok(v)
}
