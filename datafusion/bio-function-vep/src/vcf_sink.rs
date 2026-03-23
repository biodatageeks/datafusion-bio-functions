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
/// 2. Runs `annotate_vep()` via SQL to produce annotated batches
/// 3. Writes the annotated batches to a VCF file using `VcfLocalWriter`
///
/// All original VCF columns are preserved in the output. The `csq` column from
/// annotation is added as an INFO field. The 87 typed annotation columns
/// (Allele, Consequence, SYMBOL, AF, etc.) are NOT written to the VCF.
///
/// # Arguments
///
/// * `ctx` - Session context with VEP functions registered and VCF table available
/// * `vcf_table` - Name of the registered VCF table to annotate
/// * `cache_source` - Path to the VEP cache (parquet directory or fjall store)
/// * `backend` - Cache backend type: `"parquet"`
/// * `output_path` - Path to write the output VCF file
/// * `config` - Annotation and output configuration
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
    let options_json = config.to_options_json();
    annotate_to_vcf_with_options(
        ctx,
        vcf_table,
        cache_source,
        backend,
        Some(&options_json),
        output_path,
        config.compression,
    )
    .await
}

/// Lower-level annotation-to-VCF function that takes a raw options JSON string.
///
/// Prefer [`annotate_to_vcf`] with [`AnnotateVcfConfig`] for typed parameters.
pub async fn annotate_to_vcf_with_options(
    ctx: &SessionContext,
    vcf_table: &str,
    cache_source: &str,
    backend: &str,
    options_json: Option<&str>,
    output_path: &Path,
    compression: VcfCompressionType,
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
                // Extract sample name from "{sample}_{format_id}" column naming.
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

    // Fallback: get sample names from schema-level metadata.
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

    // Deduplicate format fields to unique FORMAT tags for the VCF header.
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

    // 2. Build annotation SQL — select only VCF-relevant columns + csq.
    //    The 87 typed annotation columns (Allele, Consequence, SYMBOL, AF, etc.)
    //    are Arrow-only and should NOT appear in the VCF output.
    let mut select_cols: Vec<String> = Vec::new();
    for name in &core_vcf {
        select_cols.push(format!("`{name}`"));
    }
    for name in &info_fields {
        select_cols.push(format!("`{name}`"));
    }
    // CSQ must be explicitly projected so it gets computed (skip_csq default).
    select_cols.push("`csq`".to_string());
    for name in &format_fields {
        select_cols.push(format!("`{name}`"));
    }
    let select_list = select_cols.join(", ");

    let opts_clause = options_json
        .map(|o| format!(", '{}'", o.replace('\'', "''")))
        .unwrap_or_default();
    let sql = format!(
        "SELECT {select_list} FROM annotate_vep('{vcf_table}', '{}', '{}'{opts_clause})",
        cache_source.replace('\'', "''"),
        backend.replace('\'', "''"),
    );

    let batches = ctx.sql(&sql).await?.collect().await?;

    // 3. Add CSQ to INFO fields for VCF output (annotation adds a `csq` column).
    let mut vcf_info_fields = info_fields;
    if batches
        .first()
        .is_some_and(|b| b.schema().index_of("csq").is_ok())
    {
        vcf_info_fields.push("csq".to_string());
    }

    // 4. Build the output schema reference for the header builder.
    //    We need the schema that has the field metadata — the annotation output
    //    schema loses it, so we merge: take the output schema fields but attach
    //    metadata from the VCF input schema where available, plus add `csq` with
    //    proper INFO metadata.
    let output_schema = if let Some(first) = batches.first() {
        first.schema()
    } else {
        vcf_schema.clone()
    };

    // Build a schema that carries the original VCF metadata on each field.
    let merged_fields: Vec<_> = output_schema
        .fields()
        .iter()
        .map(|f| {
            if let Ok(input_field) = vcf_schema.field_with_name(f.name()) {
                // Carry over metadata from the input VCF field.
                let mut merged_metadata = input_field.metadata().clone();
                // The output field metadata takes precedence.
                for (k, v) in f.metadata() {
                    merged_metadata.insert(k.clone(), v.clone());
                }
                f.as_ref().clone().with_metadata(merged_metadata)
            } else if f.name() == "csq" {
                // Add INFO metadata for the CSQ annotation field.
                let mut meta = f.metadata().clone();
                meta.insert("bio.vcf.field.field_type".to_string(), "INFO".to_string());
                meta.insert(
                    "bio.vcf.field.description".to_string(),
                    "Consequence annotations from annotate_vep".to_string(),
                );
                meta.insert("bio.vcf.field.number".to_string(), ".".to_string());
                meta.insert("bio.vcf.field.type".to_string(), "String".to_string());
                f.as_ref().clone().with_metadata(meta)
            } else {
                f.as_ref().clone()
            }
        })
        .collect();

    // Merge schema-level metadata from the VCF input (file format, contigs, filters, etc.)
    let mut merged_schema_metadata = vcf_schema.metadata().clone();
    for (k, v) in output_schema.metadata() {
        merged_schema_metadata.insert(k.clone(), v.clone());
    }

    let write_schema = Arc::new(
        datafusion::arrow::datatypes::Schema::new(merged_fields)
            .with_metadata(merged_schema_metadata),
    );

    // 5. Write VCF file.
    let mut writer = VcfLocalWriter::with_compression(output_path, compression)?;
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

    let mut total_rows = 0;
    for batch in &batches {
        let lines = batch_to_vcf_lines(
            batch,
            &vcf_info_fields,
            &unique_format_tags,
            &sample_names,
            coordinate_zero_based,
        )?;
        total_rows += lines.len();
        writer.write_records(&lines)?;
    }

    writer.finish()?;

    Ok(total_rows)
}
