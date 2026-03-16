//! Profile the VEP annotation pipeline with fine-grained timing.
//!
//! Usage:
//!   VEP_PROFILE=1 cargo run --release --example profile_annotation -- \
//!     /path/to/input.vcf.gz \
//!     /path/to/cache_dir \
//!     [sample_limit]
//!
//! The cache_dir should contain:
//!   115_GRCh38_variation_1_vep.parquet
//!   115_GRCh38_transcript_1_vep.parquet
//!   115_GRCh38_exon_1_vep.parquet
//!   115_GRCh38_translation_1_vep.parquet
//!   115_GRCh38_regulatory_1_vep.parquet
//!   115_GRCh38_motif_1_vep.parquet
//!
//! Example (1000 variants from chr1, --everything mode):
//!   VEP_PROFILE=1 cargo run --release --example profile_annotation -- \
//!     vep-benchmark/data/HG002_chr1.vcf.gz \
//!     /Users/mwiewior/research/data/vep/chr1-vep \
//!     1000 \
//!     --everything \
//!     --reference-fasta-path=/Users/mwiewior/research/data/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::Array;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::TableProvider;
use datafusion::prelude::{ParquetReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::golden_benchmark::sample_gz_vcf_first_n;
use datafusion_bio_function_vep::register_vep_functions;

fn sql_literal(s: &str) -> String {
    s.replace('\'', "''")
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    // Ensure profiling is on.
    if std::env::var("VEP_PROFILE").is_err() {
        // SAFETY: called before any threads are spawned.
        unsafe { std::env::set_var("VEP_PROFILE", "1") };
    }

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: VEP_PROFILE=1 cargo run --release --example profile_annotation -- \
             <vcf.gz> <cache_dir> [sample_limit] [--everything] \
             [--reference-fasta-path=<path>] [--hgvs]"
        );
        std::process::exit(1);
    }

    let vcf_gz = PathBuf::from(&args[1]);
    let cache_dir = PathBuf::from(&args[2]);
    let sample_limit: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1000);

    let everything = args.iter().any(|a| a == "--everything");
    let hgvs = args.iter().any(|a| a == "--hgvs") || everything;
    let reference_fasta_path = args
        .iter()
        .find_map(|a| a.strip_prefix("--reference-fasta-path="))
        .map(|s| s.to_string());
    let output_vcf = args
        .iter()
        .find_map(|a| a.strip_prefix("--output="))
        .map(|s| PathBuf::from(s));

    eprintln!("=== VEP Annotation Profiler ===");
    eprintln!("VCF:          {}", vcf_gz.display());
    eprintln!("Cache dir:    {}", cache_dir.display());
    eprintln!("Sample limit: {}", sample_limit);
    eprintln!("Everything:   {}", everything);
    eprintln!("HGVS:         {}", hgvs);
    eprintln!(
        "Ref FASTA:    {}",
        reference_fasta_path.as_deref().unwrap_or("none")
    );
    eprintln!();

    // Sample VCF.
    let work_dir = tempfile::tempdir()
        .map_err(|e| DataFusionError::Execution(format!("failed to create temp dir: {e}")))?;
    let sampled_vcf = work_dir.path().join("sampled.vcf");
    let t_sample = Instant::now();
    let n_sampled = sample_gz_vcf_first_n(&vcf_gz, &sampled_vcf, sample_limit)?;
    eprintln!(
        "[PREP] Sampled {} variants in {:.1}ms",
        n_sampled,
        t_sample.elapsed().as_secs_f64() * 1000.0
    );

    // Detect cache file naming (e.g. "115_GRCh38_variation_1_vep").
    let variation_file = find_parquet(&cache_dir, "variation")?;
    let base = extract_base(&variation_file);
    eprintln!("[PREP] Base name: {}", base);

    // Set up session.
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    register_vep_functions(&ctx);

    // Register VCF.
    // None = include all INFO and FORMAT/sample fields from the VCF header.
    // This ensures original VCF fields flow through the annotation pipeline.
    let vcf = VcfTableProvider::new(sampled_vcf.display().to_string(), None, None, None, false)?;
    // Capture VCF field metadata before registration (metadata is lost through the pipeline).
    let vcf_schema = vcf.schema();
    ctx.register_table("sampled_vcf", Arc::new(vcf))?;

    // Register cache parquet tables.
    let variation_path = cache_dir.join(&variation_file);
    ctx.register_parquet(
        "variation_cache",
        variation_path.display().to_string().as_str(),
        ParquetReadOptions::default(),
    )
    .await?;

    let mut options = Vec::new();

    // Register context tables.
    // Try translation_core first, then fall back to translation (unsplit layout).
    let context_tables = [
        ("transcript", "transcripts_table"),
        ("exon", "exons_table"),
        ("translation_core", "translations_table"),
        ("regulatory", "regulatory_table"),
        ("motif", "motif_table"),
    ];
    for (file_stem, json_key) in &context_tables {
        let stem_result = find_parquet_stem(&cache_dir, file_stem, &base).or_else(|_| {
            // Fallback: "translation_core" -> "translation" for unsplit layout
            if *file_stem == "translation_core" {
                find_parquet_stem(&cache_dir, "translation", &base)
            } else {
                Err(DataFusionError::Execution("not found".to_string()))
            }
        });
        if let Ok(filename) = stem_result {
            let table_name = filename.trim_end_matches(".parquet");
            let path = cache_dir.join(&filename);
            ctx.register_parquet(
                table_name,
                path.display().to_string().as_str(),
                ParquetReadOptions::default(),
            )
            .await?;
            options.push(format!("\"{}\":\"{}\"", json_key, sql_literal(table_name)));
            eprintln!("[PREP] Registered {} -> {}", json_key, filename);
        }
    }

    // Register translation_sift separately for sift/polyphen window loading.
    if let Ok(filename) = find_parquet_stem(&cache_dir, "translation_sift", &base) {
        let table_name = filename.trim_end_matches(".parquet");
        let path = cache_dir.join(&filename);
        ctx.register_parquet(
            table_name,
            path.display().to_string().as_str(),
            ParquetReadOptions::default(),
        )
        .await?;
        options.push(format!(
            "\"translations_sift_table\":\"{}\"",
            sql_literal(table_name)
        ));
        eprintln!("[PREP] Registered translations_sift_table -> {}", filename);
    }

    options.push("\"extended_probes\":true".to_string());
    options.push("\"check_existing\":true".to_string());
    options.push("\"af\":true".to_string());
    options.push("\"af_1kg\":true".to_string());
    options.push("\"af_gnomade\":true".to_string());
    options.push("\"af_gnomadg\":true".to_string());
    options.push("\"max_af\":true".to_string());
    options.push("\"pubmed\":true".to_string());

    if everything {
        options.push("\"everything\":true".to_string());
    }
    if hgvs || everything {
        let fasta = reference_fasta_path.as_deref().ok_or_else(|| {
            DataFusionError::Execution(
                "--hgvs/--everything requires --reference-fasta-path=<path>".to_string(),
            )
        })?;
        if !everything {
            options.push("\"hgvs\":true".to_string());
        }
        options.push(format!(
            "\"reference_fasta_path\":\"{}\"",
            sql_literal(fasta)
        ));
    }

    let options_json = format!("{{{}}}", options.join(","));
    eprintln!("[PREP] options_json: {}", options_json);
    eprintln!();

    // Run annotation with profiling.
    let sql = format!(
        "SELECT * FROM annotate_vep('sampled_vcf', '{}', '{}', '{}')",
        sql_literal(variation_path.display().to_string().as_str()),
        sql_literal("parquet"),
        sql_literal(&options_json)
    );

    eprintln!("=== Running annotation ===");
    let t_total = Instant::now();
    let batches = ctx.sql(&sql).await?.collect().await?;
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    let total_ms = t_total.elapsed().as_secs_f64() * 1000.0;
    eprintln!();
    eprintln!("=== Results ===");
    eprintln!("Output rows:  {}", total_rows);
    eprintln!(
        "Total time:   {:.1}ms ({:.2}s)",
        total_ms,
        total_ms / 1000.0
    );
    eprintln!(
        "Throughput:   {:.0} variants/sec",
        total_rows as f64 / (total_ms / 1000.0)
    );

    // Print schema.
    if let Some(batch) = batches.first() {
        eprintln!("Output cols:  {}", batch.schema().fields().len());
    }

    // Write output VCF if --output=<path> was specified.
    if let Some(output_path) = &output_vcf {
        let t_write = Instant::now();
        write_vcf_output(&batches, output_path, &vcf_schema)?;
        let write_ms = t_write.elapsed().as_secs_f64() * 1000.0;
        eprintln!(
            "VCF output:   {} ({:.1}ms)",
            output_path.display(),
            write_ms
        );
    }

    Ok(())
}

fn find_parquet(dir: &Path, stem_contains: &str) -> Result<String> {
    for entry in std::fs::read_dir(dir).map_err(|e| {
        DataFusionError::Execution(format!("cannot read dir {}: {e}", dir.display()))
    })? {
        let entry =
            entry.map_err(|e| DataFusionError::Execution(format!("dir entry error: {e}")))?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.contains(stem_contains) && name.ends_with(".parquet") {
            return Ok(name);
        }
    }
    Err(DataFusionError::Execution(format!(
        "no parquet file containing '{}' found in {}",
        stem_contains,
        dir.display()
    )))
}

fn find_parquet_stem(dir: &Path, stem: &str, base: &str) -> Result<String> {
    // Try "{base}_{stem}_*_vep.parquet" pattern first, then "{base}_{stem}.parquet".
    for entry in std::fs::read_dir(dir).map_err(|e| {
        DataFusionError::Execution(format!("cannot read dir {}: {e}", dir.display()))
    })? {
        let entry =
            entry.map_err(|e| DataFusionError::Execution(format!("dir entry error: {e}")))?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with(base) && name.contains(stem) && name.ends_with(".parquet") {
            return Ok(name);
        }
    }
    Err(DataFusionError::Execution(format!(
        "no parquet file for '{}' with base '{}' in {}",
        stem,
        base,
        dir.display()
    )))
}

fn string_at(array: &dyn Array, row: usize) -> Option<String> {
    use datafusion::arrow::array::{LargeStringArray, StringArray, StringViewArray};
    if array.is_null(row) {
        return None;
    }
    if let Some(a) = array.as_any().downcast_ref::<StringArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = array.as_any().downcast_ref::<LargeStringArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = array.as_any().downcast_ref::<StringViewArray>() {
        return Some(a.value(row).to_string());
    }
    None
}

fn int64_at(array: &dyn Array, row: usize) -> Option<i64> {
    use datafusion::arrow::array::{Int32Array, Int64Array, UInt32Array, UInt64Array};
    if array.is_null(row) {
        return None;
    }
    if let Some(a) = array.as_any().downcast_ref::<Int64Array>() {
        return Some(a.value(row));
    }
    if let Some(a) = array.as_any().downcast_ref::<Int32Array>() {
        return Some(a.value(row) as i64);
    }
    if let Some(a) = array.as_any().downcast_ref::<UInt64Array>() {
        return Some(a.value(row) as i64);
    }
    if let Some(a) = array.as_any().downcast_ref::<UInt32Array>() {
        return Some(a.value(row) as i64);
    }
    // Try string (some VCF providers store POS as string)
    if let Some(s) = string_at(array, row) {
        return s.parse::<i64>().ok();
    }
    None
}

fn write_vcf_output(
    batches: &[datafusion::arrow::array::RecordBatch],
    path: &Path,
    vcf_input_schema: &datafusion::arrow::datatypes::SchemaRef,
) -> Result<()> {
    use std::io::BufWriter;

    let mut file = BufWriter::new(std::fs::File::create(path).map_err(|e| {
        DataFusionError::Execution(format!("cannot create {}: {e}", path.display()))
    })?);

    let Some(first_batch) = batches.first() else {
        return Ok(());
    };
    let schema = first_batch.schema();

    // Classify columns using VCF field metadata from the ORIGINAL VCF input schema.
    // The annotation pipeline output loses metadata, so we look up each column name
    // in the input schema to find its "bio.vcf.field.field_type" = "INFO" or "FORMAT".
    let core_vcf = [
        "chrom", "start", "end", "id", "ref", "alt", "qual", "filter",
    ];
    let mut info_col_indices: Vec<(usize, String)> = Vec::new();
    let mut format_col_indices: Vec<(usize, String)> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();

    // Build a lookup from the VCF input schema.
    let vcf_field_types: std::collections::HashMap<&str, &str> = vcf_input_schema
        .fields()
        .iter()
        .filter_map(|f| {
            f.metadata()
                .get("bio.vcf.field.field_type")
                .map(|ft| (f.name().as_str(), ft.as_str()))
        })
        .collect();

    for (idx, field) in schema.fields().iter().enumerate() {
        let name = field.name().as_str();
        if core_vcf.contains(&name) {
            continue;
        }
        match vcf_field_types.get(name).copied() {
            Some("INFO") => {
                info_col_indices.push((idx, name.to_string()));
            }
            Some("FORMAT") => {
                let format_id = vcf_input_schema
                    .field_with_name(name)
                    .ok()
                    .and_then(|f| f.metadata().get("bio.vcf.field.format_id").cloned())
                    .unwrap_or_else(|| name.to_string());
                let sample = if name.ends_with(&format_id) && name.len() > format_id.len() {
                    name[..name.len() - format_id.len() - 1].to_string()
                } else {
                    String::new()
                };
                if !sample.is_empty() && !sample_names.contains(&sample) {
                    sample_names.push(sample);
                }
                format_col_indices.push((idx, name.to_string()));
            }
            _ => {}
        }
    }
    if sample_names.is_empty() && !format_col_indices.is_empty() {
        if let Some(samples_json) = vcf_input_schema.metadata().get("bio.vcf.samples") {
            if let Ok(names) = serde_json::from_str::<Vec<String>>(samples_json) {
                sample_names = names;
            }
        }
        if sample_names.is_empty() {
            sample_names.push("SAMPLE".to_string());
        }
    }

    // Write VCF header.
    let wr = |e: std::io::Error| DataFusionError::Execution(format!("write error: {e}"));
    writeln!(file, "##fileformat=VCFv4.2").map_err(wr)?;
    writeln!(
        file,
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from annotate_vep\">"
    )
    .map_err(wr)?;
    write!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").map_err(wr)?;
    if !format_col_indices.is_empty() {
        write!(file, "\tFORMAT").map_err(wr)?;
        for name in &sample_names {
            write!(file, "\t{name}").map_err(wr)?;
        }
    }
    writeln!(file).map_err(wr)?;

    // Pre-compute FORMAT tag order for each sample.
    // For single-sample: format tags are just column names (GT, GQ, DP, ...).
    // For multi-sample: extract tag from "sampleName_tag" pattern.
    let format_tags: Vec<String> = if sample_names.len() <= 1 {
        format_col_indices
            .iter()
            .map(|(_, name)| name.clone())
            .collect()
    } else {
        // Deduplicate tags from "{sample}_{tag}" column names.
        let mut tags = Vec::new();
        for (_, name) in &format_col_indices {
            if let Some(tag) = name.split('_').last() {
                if !tags.contains(&tag.to_string()) {
                    tags.push(tag.to_string());
                }
            }
        }
        tags
    };

    for batch in batches {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").ok();
        let start_idx = schema.index_of("start").ok();
        let id_idx = schema.index_of("id").ok();
        let ref_idx = schema.index_of("ref").ok();
        let alt_idx = schema.index_of("alt").ok();
        let qual_idx = schema.index_of("qual").ok();
        let filter_idx = schema.index_of("filter").ok();
        let csq_idx = schema.index_of("csq").ok();

        for row in 0..batch.num_rows() {
            let chrom = chrom_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .unwrap_or_default();
            let pos = start_idx
                .and_then(|i| int64_at(batch.column(i).as_ref(), row))
                .unwrap_or(0);
            let id = id_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| ".".to_string());
            let ref_al = ref_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .unwrap_or_default();
            let alt_al = alt_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .unwrap_or_default();
            let qual = qual_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| ".".to_string());
            let filter = filter_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| ".".to_string());
            let csq = csq_idx
                .and_then(|i| string_at(batch.column(i).as_ref(), row))
                .unwrap_or_default();

            // Build INFO field: original INFO fields + CSQ.
            let mut info_parts: Vec<String> = Vec::new();
            for (idx, name) in &info_col_indices {
                if let Some(val) = string_at(batch.column(*idx).as_ref(), row) {
                    if !val.is_empty() {
                        info_parts.push(format!("{name}={val}"));
                    }
                }
            }
            if !csq.is_empty() {
                info_parts.push(format!("CSQ={csq}"));
            }
            let info = if info_parts.is_empty() {
                ".".to_string()
            } else {
                info_parts.join(";")
            };

            write!(
                file,
                "{chrom}\t{pos}\t{id}\t{ref_al}\t{alt_al}\t{qual}\t{filter}\t{info}"
            )
            .map_err(wr)?;

            // Write FORMAT and sample columns.
            if !format_col_indices.is_empty() {
                write!(file, "\t{}", format_tags.join(":")).map_err(wr)?;
                for sample in &sample_names {
                    let mut sample_values: Vec<String> = Vec::new();
                    for tag in &format_tags {
                        // Column name: "tag" (single sample) or "sample_tag" (multi sample)
                        let col_name = if sample_names.len() <= 1 {
                            tag.clone()
                        } else {
                            format!("{sample}_{tag}")
                        };
                        let val = schema
                            .index_of(&col_name)
                            .ok()
                            .and_then(|i| string_at(batch.column(i).as_ref(), row))
                            .unwrap_or_else(|| ".".to_string());
                        sample_values.push(val);
                    }
                    write!(file, "\t{}", sample_values.join(":")).map_err(wr)?;
                }
            }
            writeln!(file).map_err(wr)?;
        }
    }

    file.flush()
        .map_err(|e| DataFusionError::Execution(format!("flush error: {e}")))?;
    Ok(())
}

fn extract_base(variation_filename: &str) -> String {
    // "115_GRCh38_variation_1_vep.parquet" -> "115_GRCh38"
    if let Some(idx) = variation_filename.find("_variation") {
        variation_filename[..idx].to_string()
    } else {
        variation_filename.trim_end_matches(".parquet").to_string()
    }
}
