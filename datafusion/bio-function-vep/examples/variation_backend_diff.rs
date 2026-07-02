//! Read-parity oracle for the Parquet variation backend (Phase 2, Task 8).
//!
//! Runs the SAME sorted-unique positions through BOTH the Parquet point-lookup
//! (`SinglePathParquetVariationLookup`) and the Lance point-lookup
//! (`SinglePathLanceVariationLookup`) against real per-chromosome shards, then
//! diffs the produced logical batches column-by-column. The two must be
//! byte-identical: same matched positions, same rows, same values for all 27
//! logical AF columns + scalars. Binary flags are compared by presence
//! (Parquet `Boolean` vs Lance `Int8` null/1).
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example variation_backend_diff -- \
//!     --parquet /path/parquet.variation/chr1.parquet \
//!     --lance   /path/variation.lance/chr1.lance \
//!     --positions /path/probe_positions.txt

use std::path::PathBuf;

use datafusion::arrow::array::{Array, BooleanArray, Int8Array, StringArray, UInt32Array};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::arrow::util::display::{ArrayFormatter, FormatOptions};
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::lance_cache::af_bundle::af_column_order;
use datafusion_bio_function_vep::lance_cache::variation_runtime::SinglePathLanceVariationLookup;
use datafusion_bio_function_vep::parquet_cache::variation_lookup::SinglePathParquetVariationLookup;

const FLAG_COLUMNS: [&str; 3] = ["failed", "somatic", "phenotype_or_disease"];

/// A rich projection: every scalar the annotation reads plus all 27 logical AF
/// columns, so the diff covers the full row. `start`/`end`/`allele_string`/
/// `failed` are appended by `ensure_runtime_projection` on both backends.
fn projection() -> Vec<String> {
    let mut p: Vec<String> = [
        "variation_name",
        "allele_string",
        "failed",
        "somatic",
        "phenotype_or_disease",
        "clin_sig",
        "clin_sig_allele",
        "clinical_impact",
        "pubmed",
        "minor_allele",
        "minor_allele_freq",
        "clinvar_ids",
        "cosmic_ids",
    ]
    .into_iter()
    .map(String::from)
    .collect();
    p.extend(af_column_order().iter().map(|s| s.to_string()));
    p
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let _ = env_logger::try_init();
    let args = parse_args()?;

    let positions = read_positions(&args.positions)?;
    eprintln!(
        "probes={} parquet={} lance={}",
        positions.len(),
        args.parquet.display(),
        args.lance.display()
    );

    // Parquet lookup.
    let pq = SinglePathParquetVariationLookup::open(&args.parquet, projection()).await?;
    let mut pq_cursor = pq.new_cursor();
    let pq_taken = pq.resolve_and_take(&positions, &mut pq_cursor).await?;

    // Lance lookup.
    let lance = SinglePathLanceVariationLookup::open(&args.lance, projection()).await?;
    let mut lance_cursor = lance.new_cursor();
    let lance_taken = lance
        .resolve_and_take(&positions, &mut lance_cursor)
        .await?;

    eprintln!(
        "matched_positions: parquet={} lance={}",
        pq_taken.resolved.matched_positions, lance_taken.resolved.matched_positions
    );
    eprintln!(
        "rows: parquet={} lance={}",
        pq_taken.batch.num_rows(),
        lance_taken.batch.num_rows()
    );

    let mut failures = 0usize;
    if pq_taken.resolved.matched_positions != lance_taken.resolved.matched_positions {
        eprintln!("MISMATCH: matched_positions differ");
        failures += 1;
    }
    if pq_taken.batch.num_rows() != lance_taken.batch.num_rows() {
        eprintln!("MISMATCH: row counts differ — reporting per-position count diffs");
        report_count_diffs(&pq_taken.batch, &lance_taken.batch)?;
        std::process::exit(2);
    }

    failures += diff_batches(&pq_taken.batch, &lance_taken.batch)?;

    if failures == 0 {
        eprintln!(
            "\nOK — Parquet and Lance variation reads are byte-identical across all columns."
        );
        Ok(())
    } else {
        eprintln!("\nFAILED — {failures} column(s) or checks mismatched.");
        std::process::exit(1)
    }
}

/// Compare two logical variation batches column-by-column after aligning rows by
/// `(start, allele_string, variation_name)`. Returns the number of columns with
/// at least one mismatched cell (also prints a per-column tally).
fn diff_batches(pq: &RecordBatch, lance: &RecordBatch) -> Result<usize> {
    let perm_pq = sort_perm(pq)?;
    let perm_lance = sort_perm(lance)?;
    let n = pq.num_rows();

    // Columns present in BOTH batches (compare by name; order-independent).
    let pq_names: Vec<String> = pq
        .schema()
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();

    let opts = FormatOptions::default().with_null("<NULL>");
    let mut mismatched_cols = 0usize;
    let mut compared = 0usize;

    for name in &pq_names {
        let (Ok(pi), Ok(li)) = (pq.schema().index_of(name), lance.schema().index_of(name)) else {
            continue; // not in both → skip
        };
        compared += 1;
        let mut cell_mismatches = 0usize;
        let mut first_example: Option<String> = None;

        if FLAG_COLUMNS.contains(&name.as_str()) {
            // Presence compare: Parquet Boolean vs Lance Int8 (null | 1).
            let pq_bool = pq.column(pi).as_any().downcast_ref::<BooleanArray>();
            let lance_i8 = lance.column(li).as_any().downcast_ref::<Int8Array>();
            for row in 0..n {
                let rp = perm_pq[row];
                let rl = perm_lance[row];
                let p_present = match pq_bool {
                    Some(a) => a.value(rp),
                    None => !pq.column(pi).is_null(rp),
                };
                let l_present = match lance_i8 {
                    Some(a) => !a.is_null(rl),
                    None => !lance.column(li).is_null(rl),
                };
                if p_present != l_present {
                    cell_mismatches += 1;
                    first_example.get_or_insert_with(|| {
                        format!("row {row}: parquet_present={p_present} lance_present={l_present}")
                    });
                }
            }
        } else {
            let fp = ArrayFormatter::try_new(pq.column(pi).as_ref(), &opts)
                .map_err(|e| DataFusionError::Execution(format!("formatter {name}: {e}")))?;
            let fl = ArrayFormatter::try_new(lance.column(li).as_ref(), &opts)
                .map_err(|e| DataFusionError::Execution(format!("formatter {name}: {e}")))?;
            for row in 0..n {
                let a = fp.value(perm_pq[row]).to_string();
                let b = fl.value(perm_lance[row]).to_string();
                if a != b {
                    cell_mismatches += 1;
                    first_example
                        .get_or_insert_with(|| format!("row {row}: parquet={a:?} lance={b:?}"));
                }
            }
        }

        if cell_mismatches > 0 {
            mismatched_cols += 1;
            eprintln!(
                "  DIFF {name}: {cell_mismatches}/{n} cells — e.g. {}",
                first_example.unwrap_or_default()
            );
        }
    }

    eprintln!(
        "compared {compared} shared columns; {} matched, {mismatched_cols} mismatched",
        compared - mismatched_cols
    );
    Ok(mismatched_cols)
}

/// Print positions whose per-position row count differs between the two result
/// batches (diagnostic for a row-count mismatch).
fn report_count_diffs(pq: &RecordBatch, lance: &RecordBatch) -> Result<()> {
    use std::collections::BTreeMap;
    let counts = |b: &RecordBatch| -> Result<BTreeMap<u32, usize>> {
        let start = b
            .column_by_name("start")
            .and_then(|c| c.as_any().downcast_ref::<UInt32Array>())
            .ok_or_else(|| DataFusionError::Execution("batch missing UInt32 start".into()))?;
        let mut m = BTreeMap::new();
        for i in 0..b.num_rows() {
            *m.entry(start.value(i)).or_insert(0) += 1;
        }
        Ok(m)
    };
    let pc = counts(pq)?;
    let lc = counts(lance)?;
    let mut keys: Vec<u32> = pc.keys().chain(lc.keys()).copied().collect();
    keys.sort_unstable();
    keys.dedup();
    let mut n = 0;
    for k in keys {
        let a = *pc.get(&k).unwrap_or(&0);
        let b = *lc.get(&k).unwrap_or(&0);
        if a != b {
            eprintln!("  pos={k}: parquet={a} lance={b}");
            n += 1;
        }
    }
    eprintln!("  {n} position(s) with differing row counts");
    Ok(())
}

/// Row permutation sorting a batch by `(start, allele_string, variation_name)`.
fn sort_perm(batch: &RecordBatch) -> Result<Vec<usize>> {
    let start = batch
        .column_by_name("start")
        .and_then(|c| c.as_any().downcast_ref::<UInt32Array>())
        .ok_or_else(|| DataFusionError::Execution("batch missing UInt32 start".into()))?;
    let allele = batch
        .column_by_name("allele_string")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let vn = batch
        .column_by_name("variation_name")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let at = |arr: Option<&StringArray>, i: usize| -> String {
        match arr {
            Some(a) if a.is_valid(i) => a.value(i).to_string(),
            _ => String::new(),
        }
    };
    let mut perm: Vec<usize> = (0..batch.num_rows()).collect();
    perm.sort_by(|&x, &y| {
        (start.value(x), at(allele, x), at(vn, x)).cmp(&(start.value(y), at(allele, y), at(vn, y)))
    });
    Ok(perm)
}

struct Args {
    parquet: PathBuf,
    lance: PathBuf,
    positions: PathBuf,
}

fn parse_args() -> Result<Args> {
    let mut parquet = None;
    let mut lance = None;
    let mut positions = None;
    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--parquet" => parquet = Some(PathBuf::from(require(&mut it, &a)?)),
            "--lance" => lance = Some(PathBuf::from(require(&mut it, &a)?)),
            "--positions" => positions = Some(PathBuf::from(require(&mut it, &a)?)),
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
        positions: positions
            .ok_or_else(|| DataFusionError::Execution("--positions required".into()))?,
    })
}

fn require(it: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    it.next()
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}

/// Read whitespace/newline-separated u32 positions, sorted-unique ascending.
fn read_positions(path: &PathBuf) -> Result<Vec<u32>> {
    let text = std::fs::read_to_string(path).map_err(|e| {
        DataFusionError::Execution(format!("read positions '{}': {e}", path.display()))
    })?;
    let mut v: Vec<u32> = text
        .split_whitespace()
        .map(|t| {
            t.parse::<u32>()
                .map_err(|e| DataFusionError::Execution(format!("bad position '{t}': {e}")))
        })
        .collect::<Result<Vec<_>>>()?;
    v.sort_unstable();
    v.dedup();
    Ok(v)
}
