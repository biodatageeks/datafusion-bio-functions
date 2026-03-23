//! VCF output sink for annotate_vep() results.
//!
//! Provides a helper function to annotate variants and write the results
//! directly to a VCF file, preserving all original VCF columns (INFO, FORMAT,
//! sample genotypes) plus the annotation columns added by the pipeline.

use std::path::Path;
use std::sync::Arc;

use datafusion::common::Result;
use datafusion::datasource::TableProvider;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::serializer::batch_to_vcf_lines;
use datafusion_bio_format_vcf::{VcfCompressionType, VcfLocalWriter};
use indicatif::{ProgressBar, ProgressStyle};

/// Configuration for VCF annotation output.
#[derive(Debug, Clone, Default)]
pub struct AnnotateVcfConfig {
    /// Enable all annotation features (80-field CSQ, SIFT, PolyPhen, etc.).
    pub everything: bool,
    /// Use interval-overlap fallback for shifted indels.
    pub extended_probes: bool,
    /// Path to indexed reference FASTA (required for `everything` / `hgvs`).
    pub reference_fasta_path: Option<String>,
    /// Use fjall KV store for variation lookup + SIFT.
    pub use_fjall: bool,
    /// Enable HGVS notation.
    pub hgvs: bool,
    /// Use merged Ensembl+RefSeq cache.
    pub merged: bool,
    /// Maximum allowed `failed` flag value from cache.
    pub failed: Option<i64>,
    /// Upstream/downstream distance for transcript overlap.
    pub distance: Option<String>,
    /// Output compression type.
    pub compression: VcfCompressionType,
    /// Show a progress bar on stderr during annotation + VCF write.
    /// When true, counts input variants and displays an indicatif progress bar.
    /// Note: enabling this causes an extra `COUNT(*)` scan of the input table
    /// before annotation begins, which adds a small overhead.
    pub show_progress: bool,
}

impl AnnotateVcfConfig {
    /// Serialize to JSON options string for annotate_vep() SQL.
    fn to_options_json(&self) -> String {
        let mut opts = serde_json::Map::new();
        opts.insert("partitioned".into(), serde_json::Value::Bool(true));
        if self.everything {
            opts.insert("everything".into(), serde_json::Value::Bool(true));
        }
        if self.extended_probes {
            opts.insert("extended_probes".into(), serde_json::Value::Bool(true));
        }
        if let Some(ref fasta) = self.reference_fasta_path {
            opts.insert(
                "reference_fasta_path".into(),
                serde_json::Value::String(fasta.clone()),
            );
        }
        if self.use_fjall {
            opts.insert("use_fjall".into(), serde_json::Value::Bool(true));
        }
        if self.hgvs {
            opts.insert("hgvs".into(), serde_json::Value::Bool(true));
        }
        if self.merged {
            opts.insert("merged".into(), serde_json::Value::Bool(true));
        }
        if let Some(failed) = self.failed {
            opts.insert(
                "failed".into(),
                serde_json::Value::Number(serde_json::Number::from(failed)),
            );
        }
        if let Some(ref dist) = self.distance {
            opts.insert("distance".into(), serde_json::Value::String(dist.clone()));
        }
        serde_json::to_string(&serde_json::Value::Object(opts)).unwrap()
    }
}

/// Annotate variants from a registered VCF table and write results to a VCF file.
///
/// This function:
/// 1. Looks up the VCF input table's schema to recover INFO/FORMAT field metadata
/// 2. Optionally counts input rows for a progress bar (`config.show_progress`)
/// 3. Runs `annotate_vep()` via SQL and streams annotated batches to a VCF file
///
/// All original VCF columns are preserved in the output. The `csq` column from
/// annotation is added as an INFO field. The 87 structured annotation columns
/// plus `most_severe_consequence` (Allele, Consequence, SYMBOL, AF, etc.) are
/// NOT written to the VCF.
///
/// # Returns
///
/// The number of rows written.
pub async fn annotate_to_vcf(
    ctx: &SessionContext,
    vcf_table: &str,
    cache_source: &str,
    backend: &str,
    output_path: &Path,
    config: &AnnotateVcfConfig,
) -> Result<usize> {
    // 1. Get VCF input schema for field classification (INFO vs FORMAT metadata).
    let vcf_provider = ctx.table_provider(vcf_table).await?;
    let vcf_schema = vcf_provider.schema();

    // Classify columns using "bio.vcf.field.field_type" metadata from the input schema.
    let core_vcf = [
        "chrom", "start", "end", "id", "ref", "alt", "qual", "filter",
    ];
    let mut info_fields: Vec<String> = Vec::new();
    let mut format_fields: Vec<String> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();

    for field in vcf_schema.fields() {
        let name = field.name().as_str();
        if core_vcf.contains(&name) {
            continue;
        }
        match field
            .metadata()
            .get("bio.vcf.field.field_type")
            .map(|s| s.as_str())
        {
            Some("INFO") => {
                info_fields.push(name.to_string());
            }
            Some("FORMAT") => {
                if let Some(format_id) = field.metadata().get("bio.vcf.field.format_id") {
                    if name.len() > format_id.len() + 1 && name.ends_with(format_id.as_str()) {
                        let sample = name[..name.len() - format_id.len() - 1].to_string();
                        if !sample.is_empty() && !sample_names.contains(&sample) {
                            sample_names.push(sample);
                        }
                    }
                }
                format_fields.push(name.to_string());
            }
            _ => {}
        }
    }

    if sample_names.is_empty() {
        if let Some(json) = vcf_schema.metadata().get("bio.vcf.samples") {
            if let Ok(names) = serde_json::from_str::<Vec<String>>(json) {
                sample_names = names;
            }
        }
    }
    if sample_names.is_empty() && !format_fields.is_empty() {
        sample_names.push("SAMPLE".to_string());
    }

    let unique_format_tags: Vec<String> = if sample_names.len() <= 1 {
        format_fields.clone()
    } else {
        let mut tags = Vec::new();
        for name in &format_fields {
            if let Some(tag) = name.rsplit('_').next() {
                let tag_str = tag.to_string();
                if !tags.contains(&tag_str) {
                    tags.push(tag_str);
                }
            }
        }
        tags
    };

    // 2. Set up progress bar if requested.
    let pb = if config.show_progress {
        let total = ctx
            .sql(&format!("SELECT COUNT(*) AS n FROM `{vcf_table}`"))
            .await?
            .collect()
            .await?[0]
            .column(0)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Int64Array>()
            .map(|a| a.value(0) as u64)
            .unwrap_or(0);
        let pb = ProgressBar::new(total);
        pb.set_style(
            ProgressStyle::with_template(
                "  {spinner:.green} {bar:40.cyan/blue} {pos}/{len} [{elapsed_precise}] (eta {eta})",
            )
            .unwrap()
            .progress_chars("##-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(200));
        pb
    } else {
        ProgressBar::hidden()
    };

    // 3. Build annotation SQL — select only VCF-relevant columns + csq.
    let mut select_cols: Vec<String> = Vec::new();
    for name in &core_vcf {
        select_cols.push(format!("`{name}`"));
    }
    for name in &info_fields {
        select_cols.push(format!("`{name}`"));
    }
    select_cols.push("`csq`".to_string());
    for name in &format_fields {
        select_cols.push(format!("`{name}`"));
    }
    let select_list = select_cols.join(", ");

    let options_json = config.to_options_json();
    let opts_clause = format!(", '{}'", options_json.replace('\'', "''"));
    let sql = format!(
        "SELECT {select_list} FROM annotate_vep('{vcf_table}', '{}', '{}'{opts_clause})",
        cache_source.replace('\'', "''"),
        backend.replace('\'', "''"),
    );

    let mut vcf_info_fields = info_fields;
    vcf_info_fields.push("csq".to_string());

    // 4. Build the output schema for the VCF header.
    let df = ctx.sql(&sql).await?;
    let df_schema = df.schema();
    let output_fields: Vec<datafusion::arrow::datatypes::Field> = df_schema
        .fields()
        .iter()
        .map(|df_field| {
            let name = df_field.name();
            let arrow_field = datafusion::arrow::datatypes::Field::new(
                name,
                df_field.data_type().clone(),
                df_field.is_nullable(),
            );
            if let Ok(input_field) = vcf_schema.field_with_name(name) {
                let mut merged_metadata = input_field.metadata().clone();
                for (k, v) in arrow_field.metadata() {
                    merged_metadata.insert(k.clone(), v.clone());
                }
                arrow_field.with_metadata(merged_metadata)
            } else if name == "csq" {
                let mut meta = std::collections::HashMap::new();
                meta.insert("bio.vcf.field.field_type".to_string(), "INFO".to_string());
                meta.insert(
                    "bio.vcf.field.description".to_string(),
                    "Consequence annotations from annotate_vep".to_string(),
                );
                meta.insert("bio.vcf.field.number".to_string(), ".".to_string());
                meta.insert("bio.vcf.field.type".to_string(), "String".to_string());
                arrow_field.with_metadata(meta)
            } else {
                arrow_field
            }
        })
        .collect();

    let write_schema = Arc::new(
        datafusion::arrow::datatypes::Schema::new(output_fields)
            .with_metadata(vcf_schema.metadata().clone()),
    );

    // 5. Write VCF header, then stream batches directly to file.
    let mut writer = VcfLocalWriter::with_compression(output_path, config.compression)?;
    writer.write_header(
        &write_schema,
        &vcf_info_fields,
        &unique_format_tags,
        &sample_names,
    )?;

    let coordinate_zero_based = vcf_schema
        .metadata()
        .get("bio.coordinate_system_zero_based")
        .is_some_and(|v| v == "true");

    use futures::StreamExt;
    let mut stream = df.execute_stream().await?;
    let mut total_rows = 0;
    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        let lines = batch_to_vcf_lines(
            &batch,
            &vcf_info_fields,
            &unique_format_tags,
            &sample_names,
            coordinate_zero_based,
        )?;
        total_rows += lines.len();
        writer.write_records(&lines)?;
        pb.inc(lines.len() as u64);
    }

    writer.finish()?;
    pb.finish_and_clear();

    Ok(total_rows)
}
