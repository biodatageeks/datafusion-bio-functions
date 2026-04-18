//! Populate exon and translation keyspaces in a **separate** context fjall
//! database from Parquet files.
//!
//! IMPORTANT: Context keyspaces are written to `<fjall_db_path>/context.fjall`,
//! NOT into `variation.fjall` directly. Adding extra keyspaces to variation.fjall
//! degrades point-lookup performance by ~10% due to block cache pressure from
//! unused keyspace metadata.
//!
//! Usage:
//!   load_context_kv <cache_dir> <exon_parquet> <translation_parquet>

use std::time::Instant;

use datafusion::arrow::array::{Array, Int64Array, ListArray, StringArray, StringViewArray};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::SessionContext;
use datafusion_bio_function_vep::kv_cache::loader::{
    load_exons_into_kv, load_translations_into_kv,
};
use datafusion_bio_function_vep::transcript_consequence::{
    ExonFeature, ProteinDomainFeature, TranslationFeature,
};

fn string_at(col: &dyn Array, row: usize) -> Option<String> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
        return Some(a.value(row).to_string());
    }
    None
}

fn int64_at(col: &dyn Array, row: usize) -> Option<i64> {
    use datafusion::arrow::array::{Int8Array, Int16Array, Int32Array};
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<Int64Array>() {
        return Some(a.value(row));
    }
    if let Some(a) = col.as_any().downcast_ref::<Int32Array>() {
        return Some(a.value(row) as i64);
    }
    if let Some(a) = col.as_any().downcast_ref::<Int16Array>() {
        return Some(a.value(row) as i64);
    }
    if let Some(a) = col.as_any().downcast_ref::<Int8Array>() {
        return Some(a.value(row) as i64);
    }
    None
}

fn read_protein_features(col: &dyn Array, row: usize) -> Vec<ProteinDomainFeature> {
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
        .downcast_ref::<datafusion::arrow::array::StructArray>();
    let Some(struct_arr) = struct_arr else {
        return Vec::new();
    };
    let analysis_col = struct_arr.column_by_name("analysis");
    let hseqname_col = struct_arr.column_by_name("hseqname");
    let Some(start_col) = struct_arr.column_by_name("start") else {
        return Vec::new();
    };
    let Some(end_col) = struct_arr.column_by_name("end") else {
        return Vec::new();
    };
    let start_arr = start_col.as_any().downcast_ref::<Int64Array>();
    let end_arr = end_col.as_any().downcast_ref::<Int64Array>();

    let mut out = Vec::with_capacity(end_off - start_off);
    for i in start_off..end_off {
        let s = start_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
        let e = end_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
        let analysis = analysis_col.and_then(|c| string_at(c.as_ref(), i));
        let hseqname = hseqname_col.and_then(|c| string_at(c.as_ref(), i));
        out.push(ProteinDomainFeature {
            analysis,
            hseqname,
            start: s.unwrap_or(0),
            end: e.unwrap_or(0),
        });
    }
    out
}

async fn load_exons_from_parquet(ctx: &SessionContext, table: &str) -> Result<Vec<ExonFeature>> {
    let batches = ctx
        .sql(&format!("SELECT * FROM `{table}`"))
        .await?
        .collect()
        .await?;
    let mut out = Vec::new();
    for batch in &batches {
        let schema = batch.schema();
        let tx_idx = schema
            .index_of("transcript_id")
            .or_else(|_| schema.index_of("stable_id"))?;
        let exon_idx = schema.index_of("exon_number")?;
        let start_idx = schema.index_of("start")?;
        let end_idx = schema.index_of("end")?;
        for row in 0..batch.num_rows() {
            let Some(transcript_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                continue;
            };
            let Some(exon_number) = int64_at(batch.column(exon_idx).as_ref(), row) else {
                continue;
            };
            let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                continue;
            };
            let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                continue;
            };
            out.push(ExonFeature {
                transcript_id,
                exon_number: exon_number as i32,
                start,
                end,
            });
        }
    }
    Ok(out)
}

async fn load_translations_from_parquet(
    ctx: &SessionContext,
    table: &str,
) -> Result<Vec<TranslationFeature>> {
    let batches = ctx
        .sql(&format!("SELECT * FROM `{table}`"))
        .await?
        .collect()
        .await?;
    let mut out = Vec::new();
    for batch in &batches {
        let schema = batch.schema();
        let tx_idx = schema
            .index_of("transcript_id")
            .or_else(|_| schema.index_of("stable_id"))?;
        let cds_len_idx = schema
            .index_of("cds_len")
            .or_else(|_| schema.index_of("cds_length"))
            .ok();
        let protein_len_idx = schema.index_of("protein_len").ok();
        let translation_seq_idx = schema.index_of("translation_seq").ok();
        let cds_seq_idx = schema
            .index_of("cds_sequence")
            .or_else(|_| schema.index_of("cds_seq"))
            .or_else(|_| schema.index_of("coding_sequence"))
            .ok();
        // Canonical (pre-BAM-edit) translation/CDS columns added upstream
        // in d26e370. Strict: absent columns propagate as `None` — legacy
        // parquet caches must be regenerated.
        let translation_seq_canonical_idx = schema.index_of("translation_seq_canonical").ok();
        let cds_seq_canonical_idx = schema.index_of("cds_sequence_canonical").ok();
        let stable_id_idx = schema.index_of("stable_id").ok();
        let version_idx = schema.index_of("version").ok();
        let pf_idx = schema.index_of("protein_features").ok();

        for row in 0..batch.num_rows() {
            let Some(transcript_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                continue;
            };
            let cds_len = cds_len_idx
                .and_then(|i| int64_at(batch.column(i).as_ref(), row))
                .and_then(|v| usize::try_from(v).ok());
            let protein_len = protein_len_idx
                .and_then(|i| int64_at(batch.column(i).as_ref(), row))
                .and_then(|v| usize::try_from(v).ok());
            let translation_seq =
                translation_seq_idx.and_then(|i| string_at(batch.column(i).as_ref(), row));
            let cds_sequence = cds_seq_idx.and_then(|i| string_at(batch.column(i).as_ref(), row));
            let stable_id = stable_id_idx.and_then(|i| string_at(batch.column(i).as_ref(), row));
            let version = version_idx
                .and_then(|i| int64_at(batch.column(i).as_ref(), row))
                .and_then(|v| i32::try_from(v).ok());
            let protein_features = pf_idx
                .map(|i| read_protein_features(batch.column(i).as_ref(), row))
                .unwrap_or_default();

            let translation_seq_canonical = translation_seq_canonical_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row));
            let cds_sequence_canonical =
                cds_seq_canonical_idx.and_then(|i| string_at(batch.column(i).as_ref(), row));
            out.push(TranslationFeature {
                transcript_id,
                cds_len,
                protein_len,
                translation_seq,
                cds_sequence,
                translation_seq_canonical,
                cds_sequence_canonical,
                stable_id,
                version,
                protein_features,
            });
        }
    }
    Ok(out)
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <cache_dir> <exon_parquet> <translation_parquet>",
            args[0]
        );
        std::process::exit(1);
    }
    let cache_dir = &args[1];
    let exon_parquet = &args[2];
    let translation_parquet = &args[3];

    // Write to a separate context.fjall database, NOT into variation.fjall.
    // Extra keyspaces in variation.fjall cause ~10% performance degradation.
    let fjall_path = format!("{}/context.fjall", cache_dir);

    let ctx = SessionContext::new();
    ctx.register_parquet("exons", exon_parquet, Default::default())
        .await?;
    ctx.register_parquet("translations", translation_parquet, Default::default())
        .await?;

    println!("Loading exons from {exon_parquet}...");
    let t0 = Instant::now();
    let exons = load_exons_from_parquet(&ctx, "exons").await?;
    println!(
        "  parsed {} exons in {:.3}s",
        exons.len(),
        t0.elapsed().as_secs_f64()
    );

    println!("Loading translations from {translation_parquet}...");
    let t0 = Instant::now();
    let translations = load_translations_from_parquet(&ctx, "translations").await?;
    println!(
        "  parsed {} translations in {:.3}s",
        translations.len(),
        t0.elapsed().as_secs_f64()
    );

    println!("Opening fjall database at {fjall_path}...");
    let db = fjall::Database::builder(&fjall_path)
        .cache_size(256 * 1024 * 1024)
        .open()
        .map_err(|e| DataFusionError::External(Box::new(e)))?;

    let t0 = Instant::now();
    let n_ex = load_exons_into_kv(&db, &exons)?;
    println!("  wrote {n_ex} exons in {:.3}s", t0.elapsed().as_secs_f64());

    let t0 = Instant::now();
    let n_tl = load_translations_into_kv(&db, &translations)?;
    println!(
        "  wrote {n_tl} translations in {:.3}s",
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    db.persist(fjall::PersistMode::SyncAll)
        .map_err(|e| DataFusionError::External(Box::new(e)))?;
    println!("  persisted in {:.3}s", t0.elapsed().as_secs_f64());

    // Let Drop run to trigger GC of old SSTs.
    drop(db);

    println!("Done.");
    Ok(())
}
