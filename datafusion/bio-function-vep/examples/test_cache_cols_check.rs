use datafusion::arrow::array::{Array, LargeStringArray, StringArray, StringViewArray};
use datafusion::prelude::{ParquetReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::register_vep_functions;
use std::sync::Arc;

fn str_at(col: &dyn Array, row: usize) -> Option<String> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = col.as_any().downcast_ref::<LargeStringArray>() {
        return Some(a.value(row).to_string());
    }
    None
}

#[tokio::main]
async fn main() -> datafusion::common::Result<()> {
    let ctx = SessionContext::new_with_config(SessionConfig::new().with_target_partitions(1));
    register_vep_functions(&ctx);

    let vcf = VcfTableProvider::new(
        "/tmp/test_cache_cols.vcf".to_string(),
        Some(vec![]),
        Some(vec![]),
        None,
        false,
    )?;
    ctx.register_table("test_vcf", Arc::new(vcf))?;

    let base = "/Users/mwiewior/research/data/vep";
    for (name, file) in &[
        (
            "115_GRCh38_transcript_22",
            "115_GRCh38_transcript_22.parquet",
        ),
        ("115_GRCh38_exon_22", "115_GRCh38_exon_22.parquet"),
        (
            "115_GRCh38_translation_22",
            "115_GRCh38_translation_22.parquet",
        ),
        (
            "115_GRCh38_regulatory_22",
            "115_GRCh38_regulatory_22.parquet",
        ),
        ("115_GRCh38_motif_22", "115_GRCh38_motif_22.parquet"),
    ] {
        ctx.register_parquet(
            *name,
            &format!("{base}/{file}"),
            ParquetReadOptions::default(),
        )
        .await?;
    }

    let cols = vec![
        "variation_name",
        "clin_sig",
        "clin_sig_allele",
        "clinical_impact",
        "phenotype_or_disease",
        "pubmed",
        "somatic",
        "minor_allele",
        "minor_allele_freq",
        "AF",
        "AFR",
        "AMR",
        "EAS",
        "EUR",
        "SAS",
        "gnomADe",
        "gnomADe_AFR",
        "gnomADe_AMR",
        "gnomADe_ASJ",
        "gnomADe_EAS",
        "gnomADe_FIN",
        "gnomADe_NFE",
        "gnomADe_SAS",
        "gnomADe_MID",
        "gnomADe_REMAINING",
        "gnomADg",
        "gnomADg_AFR",
        "gnomADg_AMI",
        "gnomADg_AMR",
        "gnomADg_ASJ",
        "gnomADg_EAS",
        "gnomADg_FIN",
        "gnomADg_MID",
        "gnomADg_NFE",
        "gnomADg_SAS",
        "gnomADg_REMAINING",
        "clinvar_ids",
        "cosmic_ids",
        "dbsnp_ids",
    ];
    let col_list = cols
        .iter()
        .map(|c| format!("`{c}`"))
        .collect::<Vec<_>>()
        .join(", ");

    let sql = format!(
        "SELECT chrom, start, alt, most_severe_consequence, {col_list} \
         FROM annotate_vep('test_vcf', \
           '{base}/115_GRCh38_variation_22.parquet', 'parquet', \
           '{{\"transcripts_table\":\"115_GRCh38_transcript_22\",\"exons_table\":\"115_GRCh38_exon_22\",\"translations_table\":\"115_GRCh38_translation_22\",\"regulatory_table\":\"115_GRCh38_regulatory_22\",\"motif_table\":\"115_GRCh38_motif_22\",\"extended_probes\":false}}')"
    );

    let batches = ctx.sql(&sql).await?.collect().await?;
    for batch in &batches {
        let schema = batch.schema();
        for row in 0..batch.num_rows() {
            let chrom = str_at(batch.column(0).as_ref(), row).unwrap_or_default();
            let start = str_at(batch.column(1).as_ref(), row).unwrap_or_default();
            let alt = str_at(batch.column(2).as_ref(), row).unwrap_or_default();
            let most = str_at(batch.column(3).as_ref(), row).unwrap_or_default();
            println!(
                "=== {}:{} alt={} most_severe={} ===",
                chrom, start, alt, most
            );
            for (i, &col_name) in cols.iter().enumerate() {
                let col_idx = 4 + i;
                let val = str_at(batch.column(col_idx).as_ref(), row);
                let null_count = batch.column(col_idx).null_count();
                let dtype = batch.column(col_idx).data_type().clone();
                println!(
                    "  {:30} = {:?}  (type={:?}, nulls={})",
                    col_name, val, dtype, null_count
                );
            }
            println!();
        }
    }
    Ok(())
}
