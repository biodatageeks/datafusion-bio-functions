//! Build-time keep-first dedup over the normalized plugin stream.
//!
//! AlphaMissense (and any per-transcript plugin) can carry two source rows for
//! the **same** runtime probe key `(start, allele_string, <match cols…>)` with
//! **different** values — e.g. two overlapping UniProts sharing a genomic variant
//! and amino-acid change but scoring it differently. VEP resolves this by taking
//! the **first record in file order** (`AlphaMissense.pm` `return`s on the first
//! aa-change match). Our runtime slice, by contrast, inserts each shard row into a
//! `HashMap` keyed on the probe key, so a duplicate key silently keeps whichever
//! row the shard emits *last* — the opposite of VEP.
//!
//! This pass collapses each key to its **first** occurrence in stream order. It
//! MUST run over a single-partition, file-ordered stream (the builder sets
//! `target_partitions = 1`) and BEFORE the tier LEFT-JOIN, which reorders rows and
//! would otherwise destroy the file-order tiebreak.

use std::collections::HashSet;

use datafusion::arrow::array::BooleanArray;
use datafusion::arrow::compute::filter_record_batch;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::arrow::util::display::{ArrayFormatter, FormatOptions};
use datafusion::common::{DataFusionError, Result};
use datafusion::physical_plan::SendableRecordBatchStream;
use futures::StreamExt;

/// ASCII unit separator between key fields — will not occur in genomic contig /
/// allele / amino-acid-change strings, so it can't collide two distinct keys.
const KEY_SEP: char = '\u{1f}';

/// Consume `stream` (single-partition, source-file order) and return its batches
/// with only the **first** row per `(start, allele_string, <match col values…>)`
/// key retained — matching VEP's first-in-file rule. Later duplicates are dropped.
///
/// The key is exactly the runtime probe key (`PluginBufferSlice::probe`): `start`,
/// `allele_string`, then each match-discriminator column in manifest order. `chrom`
/// is constant within a per-chrom build and `end` is not part of the probe key, so
/// neither participates.
pub async fn dedup_keep_first(
    mut stream: SendableRecordBatchStream,
    match_columns: &[String],
) -> Result<Vec<RecordBatch>> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut out: Vec<RecordBatch> = Vec::new();
    // Render nulls to a token distinct from an empty string so `Some("")` and
    // `None` never collapse into the same key.
    let opts = FormatOptions::default().with_null("\u{0}NULL");

    let mut key_names: Vec<&str> = vec!["start", "allele_string"];
    key_names.extend(match_columns.iter().map(|s| s.as_str()));

    while let Some(batch) = stream.next().await {
        let batch = batch?;
        let schema = batch.schema();
        // Build the keep mask in an inner scope so the `formatters` (which borrow
        // `batch`) drop before we move/filter `batch` below.
        let (keep, all_kept) = {
            let formatters = key_names
                .iter()
                .map(|name| {
                    let idx = schema.index_of(name)?;
                    ArrayFormatter::try_new(batch.column(idx).as_ref(), &opts).map_err(|e| {
                        DataFusionError::Execution(format!("key formatter {name}: {e}"))
                    })
                })
                .collect::<Result<Vec<_>>>()?;

            let mut keep = Vec::with_capacity(batch.num_rows());
            let mut all_kept = true;
            let mut key = String::new();
            for r in 0..batch.num_rows() {
                use std::fmt::Write;
                key.clear();
                for f in &formatters {
                    let _ = write!(key, "{}{KEY_SEP}", f.value(r));
                }
                // `insert` returns true iff the key was not already present.
                let first = seen.insert(std::mem::take(&mut key));
                all_kept &= first;
                keep.push(first);
            }
            (keep, all_kept)
        };

        if all_kept {
            out.push(batch);
        } else if keep.iter().any(|&k| k) {
            let mask = BooleanArray::from(keep);
            out.push(
                filter_record_batch(&batch, &mask)
                    .map_err(|e| DataFusionError::Execution(format!("dedup filter: {e}")))?,
            );
        }
        // else: every row was a duplicate — drop the whole batch.
    }
    Ok(out)
}

/// Bound on how many distinct keys `check_assume_unique_sample` tracks before
/// it stops checking. An `assume_unique` source skips `dedup_keep_first`
/// entirely to avoid materializing a `HashSet<String>` sized to the whole
/// chromosome (the dominant remaining memory cost on CADD/SpliceAI's largest
/// chromosomes) -- so this check trades exhaustive coverage for a hard,
/// small memory ceiling: it validates the FIRST `SAMPLE_CAP` distinct keys
/// seen (in source-file order) and then stops, rather than either costing
/// nothing (leaving the assumption fully untested) or costing exactly what
/// the flag exists to avoid (an exhaustive `HashSet`).
const SAMPLE_CAP: usize = 2_000_000;

/// Validate (a bounded prefix of) an `assume_unique` source's claim that it
/// never repeats a runtime probe key, without filtering anything.
///
/// A violated assumption is caught here rather than left to the runtime
/// lookup `HashMap`, which keeps the *last* row of a duplicate key --
/// silently inverting VEP's first-in-file rule (see module docs). This
/// checks only the first `SAMPLE_CAP` distinct keys encountered so memory
/// stays bounded regardless of chromosome size; a source whose duplicates
/// only occur past that prefix will not be caught by this pass.
pub async fn check_assume_unique_sample(
    mut stream: SendableRecordBatchStream,
    match_columns: &[String],
) -> Result<Vec<RecordBatch>> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut out: Vec<RecordBatch> = Vec::new();
    let opts = FormatOptions::default().with_null("\u{0}NULL");

    let mut key_names: Vec<&str> = vec!["start", "allele_string"];
    key_names.extend(match_columns.iter().map(|s| s.as_str()));

    while let Some(batch) = stream.next().await {
        let batch = batch?;
        if seen.len() < SAMPLE_CAP {
            let schema = batch.schema();
            let formatters = key_names
                .iter()
                .map(|name| {
                    let idx = schema.index_of(name)?;
                    ArrayFormatter::try_new(batch.column(idx).as_ref(), &opts).map_err(|e| {
                        DataFusionError::Execution(format!("key formatter {name}: {e}"))
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            let mut key = String::new();
            for r in 0..batch.num_rows() {
                if seen.len() >= SAMPLE_CAP {
                    break;
                }
                use std::fmt::Write;
                key.clear();
                for f in &formatters {
                    let _ = write!(key, "{}{KEY_SEP}", f.value(r));
                }
                if !seen.insert(std::mem::take(&mut key)) {
                    return Err(DataFusionError::Execution(
                        "assume_unique source violated its own claim: a duplicate runtime \
                         probe key was found in the first sampled rows. The manifest's \
                         `assume_unique = true` is unsafe for this source -- remove it so \
                         the keep-first dedup pass runs (VEP's first-in-file rule), or fix \
                         the source if the duplicate is unexpected."
                            .into(),
                    ));
                }
            }
        }
        out.push(batch);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float32Array, Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::prelude::SessionContext;
    use std::sync::Arc;

    /// Normalized-shaped batch: chrom, start, end, allele_string, protein_variant,
    /// am_pathogenicity. Two rows share the runtime key (start 101, A/G, H101Y)
    /// but score differently; the first must win.
    fn norm_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("protein_variant", DataType::Utf8, true),
            Field::new("am_pathogenicity", DataType::Float32, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["3", "3", "3", "3"])),
                Arc::new(Int64Array::from(vec![101i64, 101, 101, 202])),
                Arc::new(Int64Array::from(vec![101i64, 101, 101, 202])),
                Arc::new(StringArray::from(vec!["C/T", "C/T", "C/T", "G/A"])),
                // row0/row1 share H101Y (dup key, keep first=0.0431);
                // row2 is a different aa-change at same pos → kept;
                // row3 different pos → kept.
                Arc::new(StringArray::from(vec![
                    Some("H101Y"),
                    Some("H101Y"),
                    Some("K55N"),
                    Some("D43Y"),
                ])),
                Arc::new(Float32Array::from(vec![
                    Some(0.0431f32),
                    Some(0.0898),
                    Some(0.7),
                    Some(0.12),
                ])),
            ],
        )
        .unwrap()
    }

    async fn stream_of(batch: RecordBatch) -> SendableRecordBatchStream {
        // Single-partition MemTable → single-partition, in-order stream.
        let ctx = SessionContext::new();
        ctx.register_batch("t", batch).unwrap();
        ctx.sql("SELECT * FROM t")
            .await
            .unwrap()
            .execute_stream()
            .await
            .unwrap()
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn keeps_first_of_duplicate_key() {
        let stream = stream_of(norm_batch()).await;
        let out = dedup_keep_first(stream, &["protein_variant".to_string()])
            .await
            .unwrap();
        let rows: usize = out.iter().map(|b| b.num_rows()).sum();
        assert_eq!(rows, 3, "one duplicate H101Y row must be dropped");

        // Collect (protein_variant, score) survivors.
        let mut survivors: Vec<(String, f32)> = Vec::new();
        for b in &out {
            let pv = b
                .column(b.schema().index_of("protein_variant").unwrap())
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let sc = b
                .column(b.schema().index_of("am_pathogenicity").unwrap())
                .as_any()
                .downcast_ref::<Float32Array>()
                .unwrap();
            for r in 0..b.num_rows() {
                survivors.push((pv.value(r).to_string(), sc.value(r)));
            }
        }
        // H101Y appears once, with the FIRST score (0.0431), not 0.0898.
        let h101y: Vec<_> = survivors.iter().filter(|(pv, _)| pv == "H101Y").collect();
        assert_eq!(h101y.len(), 1);
        assert!((h101y[0].1 - 0.0431).abs() < 1e-6, "kept {:?}", h101y[0]);
        // The distinct aa-change at the same position survives (not over-deduped).
        assert!(survivors.iter().any(|(pv, _)| pv == "K55N"));
        assert!(survivors.iter().any(|(pv, _)| pv == "D43Y"));
    }

    // I2 fix: an `assume_unique` source that actually repeats a probe key must
    // be caught, not silently accepted (the runtime lookup would then keep the
    // *last* duplicate, inverting VEP's first-in-file rule).
    #[tokio::test(flavor = "multi_thread")]
    async fn check_assume_unique_sample_rejects_a_real_duplicate() {
        let stream = stream_of(norm_batch()).await;
        let err = check_assume_unique_sample(stream, &["protein_variant".to_string()])
            .await
            .unwrap_err();
        assert!(err.to_string().contains("violated its own claim"));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn check_assume_unique_sample_passes_genuinely_unique_input() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("score", DataType::Float32, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![50i64, 51])),
                Arc::new(Int64Array::from(vec![50i64, 51])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Float32Array::from(vec![Some(1.0f32), Some(2.0)])),
            ],
        )
        .unwrap();
        let stream = stream_of(batch).await;
        let out = check_assume_unique_sample(stream, &[]).await.unwrap();
        let rows: usize = out.iter().map(|b| b.num_rows()).sum();
        assert_eq!(rows, 2, "unique input must pass through untouched");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn per_variant_plugin_dedups_on_start_allele_only() {
        // No match columns: two rows same (start, allele) → keep first.
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("score", DataType::Float32, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![50i64, 50])),
                Arc::new(Int64Array::from(vec![50i64, 50])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
                Arc::new(Float32Array::from(vec![Some(1.0f32), Some(2.0)])),
            ],
        )
        .unwrap();
        let stream = stream_of(batch).await;
        let out = dedup_keep_first(stream, &[]).await.unwrap();
        let rows: usize = out.iter().map(|b| b.num_rows()).sum();
        assert_eq!(rows, 1);
        let sc = out[0]
            .column(out[0].schema().index_of("score").unwrap())
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        assert!((sc.value(0) - 1.0).abs() < 1e-6);
    }
}
