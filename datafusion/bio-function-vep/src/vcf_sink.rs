//! VCF output sink for annotate_vep() results.
//!
//! Provides [`annotate_to_vcf`] — a single-call function that reads a VCF,
//! annotates it, and streams results to an output VCF file.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use std::time::{Duration, Instant};

use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::TableProvider;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::serializer::{VcfRecordLine, batch_to_vcf_lines};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_format_vcf::{VcfCompressionType, VcfLocalWriter};
use indicatif::{ProgressBar, ProgressStyle};

use crate::cache_source::CacheSourceType;

/// Callback invoked after each batch is written to VCF.
/// Arguments: (rows_in_batch, total_rows_written_so_far, total_input_rows).
/// `total_input_rows` is 0 if the count was not computed (show_progress=false).
/// Used by Python wrappers (vepyr) to drive tqdm progress bars in Jupyter.
pub type OnBatchWritten = Box<dyn Fn(usize, usize, usize) + Send + Sync>;

/// Ensembl VEP release/115 default `--buffer_size`.
pub const VEP_DEFAULT_BUFFER_SIZE: usize = 5000;

fn sink_profile_enabled() -> bool {
    std::env::var("VEP_PROFILE").is_ok() || std::env::var("VEP_VCF_PROFILE").is_ok()
}

#[derive(Debug, Default)]
struct VcfSinkProfile {
    stream_next: Duration,
    batch_to_lines: Duration,
    format_wait: Duration,
    write_records: Duration,
    writer_finish: Duration,
    batches: usize,
    rows: usize,
    lines: usize,
    body_chunk_bytes: usize,
    format_jobs: usize,
    format_inflight_max: usize,
}

impl VcfSinkProfile {
    fn summary_line(&self) -> String {
        format!(
            "[VEP_PROFILE] vcf_sink_profile batches={} rows={} lines={} body_chunk_bytes={} format_jobs={} format_inflight_max={} stream_next={:.3}s batch_to_lines={:.3}s format_wait={:.3}s write_records={:.3}s writer_finish={:.3}s",
            self.batches,
            self.rows,
            self.lines,
            self.body_chunk_bytes,
            self.format_jobs,
            self.format_inflight_max,
            self.stream_next.as_secs_f64(),
            self.batch_to_lines.as_secs_f64(),
            self.format_wait.as_secs_f64(),
            self.write_records.as_secs_f64(),
            self.writer_finish.as_secs_f64(),
        )
    }
}

#[derive(Debug)]
struct FormattedVcfBatch {
    batch_id: usize,
    input_rows: usize,
    lines: usize,
    bytes: Vec<u8>,
    format_duration: Duration,
}

fn parallel_vcf_format_enabled(config: &AnnotateVcfConfig) -> bool {
    config.target_partitions > 1
        && std::env::var("VEP_VCF_PARALLEL_FORMAT")
            .map(|value| value != "0")
            .unwrap_or(true)
        && std::env::var("VEP_VCF_BYTE_CHUNKS")
            .map(|value| value != "0")
            .unwrap_or(true)
}

fn parallel_vcf_format_jobs(config: &AnnotateVcfConfig) -> usize {
    config.target_partitions.max(1)
}

fn vcf_lines_to_body_chunk(lines: Vec<VcfRecordLine>) -> (Vec<u8>, usize) {
    let line_count = lines.len();
    let byte_len = lines.iter().map(|record| record.line.len() + 1).sum();
    let mut bytes = Vec::with_capacity(byte_len);
    for record in lines {
        bytes.extend_from_slice(record.line.as_bytes());
        bytes.push(b'\n');
    }
    (bytes, line_count)
}

fn format_vcf_body_chunk(
    batch_id: usize,
    batch: RecordBatch,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
) -> Result<FormattedVcfBatch> {
    let input_rows = batch.num_rows();
    let format_started = Instant::now();
    let lines = batch_to_vcf_lines(
        &batch,
        vcf_info_fields.as_slice(),
        unique_format_tags.as_slice(),
        sample_names.as_slice(),
        coordinate_zero_based,
    )?;
    let (bytes, line_count) = vcf_lines_to_body_chunk(lines);
    Ok(FormattedVcfBatch {
        batch_id,
        input_rows,
        lines: line_count,
        bytes,
        format_duration: format_started.elapsed(),
    })
}

fn write_vcf_body_chunk(writer: &mut VcfLocalWriter, bytes: &[u8]) -> Result<()> {
    if bytes.is_empty() {
        return Ok(());
    }

    let write_result = match writer {
        VcfLocalWriter::Plain(writer) => writer.write_all(bytes),
        VcfLocalWriter::Gzip(writer) => writer.write_all(bytes),
        VcfLocalWriter::Bgzf(writer) => writer.write_all(bytes),
    };

    write_result.map_err(|e| DataFusionError::Execution(format!("Failed to write VCF chunk: {e}")))
}

fn drain_ready_vcf_chunks(
    ready: &mut BTreeMap<usize, FormattedVcfBatch>,
    next_write_batch_id: &mut usize,
    writer: &mut VcfLocalWriter,
    pb: &ProgressBar,
    config: &AnnotateVcfConfig,
    total_input: usize,
    total_rows: &mut usize,
    sink_profile: &mut Option<VcfSinkProfile>,
) -> Result<()> {
    while let Some(chunk) = ready.remove(next_write_batch_id) {
        let write_started = Instant::now();
        write_vcf_body_chunk(writer, &chunk.bytes)?;
        if let Some(profile) = sink_profile.as_mut() {
            profile.write_records += write_started.elapsed();
            profile.batch_to_lines += chunk.format_duration;
            profile.batches += 1;
            profile.rows += chunk.input_rows;
            profile.lines += chunk.lines;
            profile.body_chunk_bytes += chunk.bytes.len();
        }
        *total_rows += chunk.lines;
        pb.inc(chunk.lines as u64);
        if let Some(ref cb) = config.on_batch_written {
            cb(chunk.lines, *total_rows, total_input);
        }
        *next_write_batch_id += 1;
    }

    Ok(())
}

/// Configuration for VCF annotation output.
pub struct AnnotateVcfConfig {
    /// Enable all annotation features (80-field CSQ, SIFT, PolyPhen, etc.).
    pub everything: bool,
    /// Emit one consequence per variant (`--pick`).
    pub pick: bool,
    /// Emit one consequence per allele (`--pick_allele`).
    pub pick_allele: bool,
    /// Emit one consequence per gene and retain non-transcript rows (`--per_gene`).
    pub per_gene: bool,
    /// Emit one consequence per allele and gene (`--pick_allele_gene`).
    pub pick_allele_gene: bool,
    /// Add standalone `PICK=1` marker for one consequence per variant.
    pub flag_pick: bool,
    /// Add standalone `PICK=1` marker for one consequence per allele.
    pub flag_pick_allele: bool,
    /// Add standalone `PICK=1` markers for VEP `--flag_pick_allele_gene`.
    /// VEP also marks retained non-transcript regulatory/motif/intergenic rows.
    pub flag_pick_allele_gene: bool,
    /// Override Ensembl VEP's default pick-order ranking.
    pub pick_order: Option<String>,
    /// Use interval-overlap fallback for shifted indels.
    pub extended_probes: bool,
    /// Path to indexed reference FASTA (required for `everything` / `hgvs`).
    pub reference_fasta_path: Option<String>,
    /// Use fjall KV store for variation lookup + SIFT.
    pub use_fjall: bool,
    /// Enable HGVS notation.
    pub hgvs: bool,
    /// Enable transcript HGVS notation explicitly.
    pub hgvsc: bool,
    /// Enable protein HGVS notation explicitly.
    pub hgvsp: bool,
    /// Enable 3' HGVS shifting explicitly.
    pub shift_hgvs: Option<bool>,
    /// Don't URI-escape HGVS output.
    pub no_escape: bool,
    /// Remove version from HGVSp IDs.
    pub remove_hgvsp_version: bool,
    /// Format HGVSp using prediction-style parentheses.
    pub hgvsp_use_prediction: bool,
    /// Use RefSeq cache/transcripts in place of Ensembl transcripts.
    pub refseq: bool,
    /// Use merged Ensembl+RefSeq cache.
    pub merged: bool,
    /// Restrict to GENCODE basic transcripts.
    pub gencode_basic: bool,
    /// Restrict to GENCODE primary transcripts.
    pub gencode_primary: bool,
    /// Keep all RefSeq transcripts, including CCDS/EST-style rows.
    pub all_refseq: bool,
    /// Exclude predicted RefSeq transcripts (XM_/XR_).
    pub exclude_predicted: bool,
    /// Maximum allowed `failed` flag value from cache.
    pub failed: Option<i64>,
    /// Upstream/downstream distance for transcript overlap.
    pub distance: Option<String>,
    /// Number of input variants per VEP-style annotation buffer.
    pub buffer_size: usize,
    /// DataFusion target partitions for the annotation session.
    pub target_partitions: usize,
    /// Output compression type.
    pub compression: VcfCompressionType,
    /// Show an indicatif progress bar on stderr (for Rust CLI).
    /// For Python/Jupyter, use `on_batch_written` with tqdm instead.
    pub show_progress: bool,
    /// Optional callback invoked after each batch is written.
    /// Used by Python wrappers to drive tqdm progress bars that work in Jupyter.
    pub on_batch_written: Option<OnBatchWritten>,
}

impl Default for AnnotateVcfConfig {
    fn default() -> Self {
        Self {
            everything: false,
            pick: false,
            pick_allele: false,
            per_gene: false,
            pick_allele_gene: false,
            flag_pick: false,
            flag_pick_allele: false,
            flag_pick_allele_gene: false,
            pick_order: None,
            extended_probes: false,
            reference_fasta_path: None,
            use_fjall: false,
            hgvs: false,
            hgvsc: false,
            hgvsp: false,
            shift_hgvs: None,
            no_escape: false,
            remove_hgvsp_version: false,
            hgvsp_use_prediction: false,
            refseq: false,
            merged: false,
            gencode_basic: false,
            gencode_primary: false,
            all_refseq: false,
            exclude_predicted: false,
            failed: None,
            distance: None,
            buffer_size: VEP_DEFAULT_BUFFER_SIZE,
            target_partitions: 1,
            compression: VcfCompressionType::Plain,
            show_progress: false,
            on_batch_written: None,
        }
    }
}

impl std::fmt::Debug for AnnotateVcfConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AnnotateVcfConfig")
            .field("everything", &self.everything)
            .field("buffer_size", &self.buffer_size)
            .field("target_partitions", &self.target_partitions)
            .field("compression", &self.compression)
            .field("show_progress", &self.show_progress)
            .field("on_batch_written", &self.on_batch_written.is_some())
            .finish()
    }
}

impl AnnotateVcfConfig {
    fn to_options_json(&self) -> String {
        let mut opts = serde_json::Map::new();
        opts.insert("partitioned".into(), serde_json::Value::Bool(true));
        if self.everything {
            opts.insert("everything".into(), serde_json::Value::Bool(true));
        }
        for (key, enabled) in [
            ("pick", self.pick),
            ("pick_allele", self.pick_allele),
            ("per_gene", self.per_gene),
            ("pick_allele_gene", self.pick_allele_gene),
            ("flag_pick", self.flag_pick),
            ("flag_pick_allele", self.flag_pick_allele),
            ("flag_pick_allele_gene", self.flag_pick_allele_gene),
        ] {
            if enabled {
                opts.insert(key.into(), serde_json::Value::Bool(true));
            }
        }
        if let Some(ref pick_order) = self.pick_order {
            opts.insert(
                "pick_order".into(),
                serde_json::Value::String(pick_order.clone()),
            );
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
        if self.hgvsc {
            opts.insert("hgvsc".into(), serde_json::Value::Bool(true));
        }
        if self.hgvsp {
            opts.insert("hgvsp".into(), serde_json::Value::Bool(true));
        }
        if let Some(shift_hgvs) = self.shift_hgvs {
            opts.insert("shift_hgvs".into(), serde_json::Value::Bool(shift_hgvs));
        }
        if self.no_escape {
            opts.insert("no_escape".into(), serde_json::Value::Bool(true));
        }
        if self.remove_hgvsp_version {
            opts.insert("remove_hgvsp_version".into(), serde_json::Value::Bool(true));
        }
        if self.hgvsp_use_prediction {
            opts.insert("hgvsp_use_prediction".into(), serde_json::Value::Bool(true));
        }
        if self.gencode_basic {
            opts.insert("gencode_basic".into(), serde_json::Value::Bool(true));
        }
        if self.gencode_primary {
            opts.insert("gencode_primary".into(), serde_json::Value::Bool(true));
        }
        if self.all_refseq {
            opts.insert("all_refseq".into(), serde_json::Value::Bool(true));
        }
        if self.exclude_predicted {
            opts.insert("exclude_predicted".into(), serde_json::Value::Bool(true));
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
        opts.insert(
            "buffer_size".into(),
            serde_json::Value::Number(serde_json::Number::from(self.buffer_size)),
        );
        serde_json::to_string(&serde_json::Value::Object(opts)).unwrap()
    }

    fn include_pick_output(&self) -> bool {
        self.flag_pick || self.flag_pick_allele || self.flag_pick_allele_gene
    }
}

fn csq_header_description(
    config: &AnnotateVcfConfig,
    cache_source_type: CacheSourceType,
) -> String {
    let field_names = crate::golden_benchmark::csq_field_names_for_mode_with_pick(
        config.everything,
        cache_source_type == CacheSourceType::RefSeq,
        cache_source_type == CacheSourceType::Merged,
        config.include_pick_output(),
    );
    let format_list = field_names.join("|");
    format!("Consequence annotations from annotate_vep. Format: {format_list}")
}

/// Annotate a VCF file and write results to an output VCF.
///
/// Handles everything in a single call:
/// 1. Reads the input VCF (all INFO/FORMAT fields preserved)
/// 2. Creates a session, registers VEP functions and the VCF table
/// 3. Runs annotation and streams results to the output VCF
///
/// The 87 structured annotation columns plus `most_severe_consequence`
/// are NOT written to the VCF — only core VCF columns, original INFO/FORMAT
/// fields, and the `csq` annotation are included.
///
/// Cache source mode is read from Arrow schema metadata on a parquet file under
/// `{cache_source}/variation`. That directory must be readable even when
/// `backend` selects a non-parquet annotation store such as `fjall`.
///
/// # Returns
///
/// The number of rows written.
pub async fn annotate_to_vcf(
    input_vcf: &str,
    cache_source: &str,
    backend: &str,
    output_vcf: &str,
    config: &AnnotateVcfConfig,
) -> Result<usize> {
    if config.refseq || config.merged {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): refseq and merged config fields are unsupported; cache source mode must come from cache schema metadata bio.vep.cache_source_type".to_string(),
        ));
    }
    let cache_source_type = CacheSourceType::from_partitioned_cache_source(cache_source)?;

    // 1. Create session and register VCF table.
    let session_config = datafusion::prelude::SessionConfig::new()
        .with_target_partitions(config.target_partitions.max(1));
    let ctx = SessionContext::new_with_config(session_config);
    crate::register_vep_functions(&ctx);

    let vcf_path = input_vcf.to_string();
    let vcf_provider = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(vcf_path, None, None, None, false)
    })
    .await
    .map_err(|e| datafusion::common::DataFusionError::External(Box::new(e)))??;

    let vcf_schema = vcf_provider.schema();
    ctx.register_table("__vep_vcf", Arc::new(vcf_provider))?;
    let vcf_table = "__vep_vcf";

    // 2. Classify columns from VCF schema metadata.
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

    // 3. Count input rows (for progress bar and/or callback total).
    let need_count = config.show_progress || config.on_batch_written.is_some();
    let total_input: usize = if need_count {
        ctx.sql(&format!("SELECT COUNT(*) AS n FROM `{vcf_table}`"))
            .await?
            .collect()
            .await?[0]
            .column(0)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Int64Array>()
            .map(|a| a.value(0) as usize)
            .unwrap_or(0)
    } else {
        0
    };

    let pb = if config.show_progress {
        let pb = ProgressBar::new(total_input as u64);
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

    // 4. Build annotation SQL — only VCF-relevant columns + csq.
    let mut select_cols: Vec<String> = Vec::new();
    for name in &core_vcf {
        select_cols.push(format!("`{name}`"));
    }
    for name in &info_fields {
        select_cols.push(format!("`{name}`"));
    }
    select_cols.push("\"CSQ\"".to_string());
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
    vcf_info_fields.push("CSQ".to_string());

    // 5. Build output schema with merged metadata for VCF header.
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
            } else if name == "CSQ" {
                let description = csq_header_description(config, cache_source_type);
                let mut meta = std::collections::HashMap::new();
                meta.insert("bio.vcf.field.field_type".to_string(), "INFO".to_string());
                meta.insert("bio.vcf.field.description".to_string(), description);
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

    // 6. Stream annotated batches to VCF file.
    let output_path = Path::new(output_vcf);
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
    let mut sink_profile = sink_profile_enabled().then(VcfSinkProfile::default);

    if parallel_vcf_format_enabled(config) {
        let vcf_info_fields = Arc::new(vcf_info_fields);
        let unique_format_tags = Arc::new(unique_format_tags);
        let sample_names = Arc::new(sample_names);
        let max_format_jobs = parallel_vcf_format_jobs(config);
        let mut format_jobs = tokio::task::JoinSet::new();
        let mut ready_chunks = BTreeMap::new();
        let mut next_input_batch_id = 0;
        let mut next_write_batch_id = 0;
        let mut input_done = false;

        loop {
            while !input_done && format_jobs.len() < max_format_jobs {
                let stream_started = Instant::now();
                let next_batch = stream.next().await;
                if let Some(profile) = sink_profile.as_mut() {
                    profile.stream_next += stream_started.elapsed();
                }

                let Some(batch_result) = next_batch else {
                    input_done = true;
                    break;
                };

                let batch = batch_result?;
                let batch_id = next_input_batch_id;
                next_input_batch_id += 1;
                let vcf_info_fields = Arc::clone(&vcf_info_fields);
                let unique_format_tags = Arc::clone(&unique_format_tags);
                let sample_names = Arc::clone(&sample_names);
                format_jobs.spawn_blocking(move || {
                    format_vcf_body_chunk(
                        batch_id,
                        batch,
                        vcf_info_fields,
                        unique_format_tags,
                        sample_names,
                        coordinate_zero_based,
                    )
                });
                if let Some(profile) = sink_profile.as_mut() {
                    profile.format_jobs += 1;
                    profile.format_inflight_max =
                        profile.format_inflight_max.max(format_jobs.len());
                }
            }

            drain_ready_vcf_chunks(
                &mut ready_chunks,
                &mut next_write_batch_id,
                &mut writer,
                &pb,
                config,
                total_input,
                &mut total_rows,
                &mut sink_profile,
            )?;

            if input_done && format_jobs.is_empty() {
                break;
            }

            if format_jobs.is_empty() {
                continue;
            }

            let wait_started = Instant::now();
            let join_result = format_jobs.join_next().await;
            if let Some(profile) = sink_profile.as_mut() {
                profile.format_wait += wait_started.elapsed();
            }
            let Some(join_result) = join_result else {
                continue;
            };
            let formatted = join_result.map_err(|e| {
                DataFusionError::Execution(format!("VCF formatter task failed: {e}"))
            })??;
            ready_chunks.insert(formatted.batch_id, formatted);
        }

        drain_ready_vcf_chunks(
            &mut ready_chunks,
            &mut next_write_batch_id,
            &mut writer,
            &pb,
            config,
            total_input,
            &mut total_rows,
            &mut sink_profile,
        )?;
    } else {
        loop {
            let stream_started = Instant::now();
            let next_batch = stream.next().await;
            if let Some(profile) = sink_profile.as_mut() {
                profile.stream_next += stream_started.elapsed();
            }
            let Some(batch_result) = next_batch else {
                break;
            };
            let batch = batch_result?;
            let input_rows = batch.num_rows();
            let lines_started = Instant::now();
            let lines = batch_to_vcf_lines(
                &batch,
                &vcf_info_fields,
                &unique_format_tags,
                &sample_names,
                coordinate_zero_based,
            )?;
            if let Some(profile) = sink_profile.as_mut() {
                profile.batch_to_lines += lines_started.elapsed();
                profile.batches += 1;
                profile.rows += input_rows;
                profile.lines += lines.len();
            }
            total_rows += lines.len();
            let write_started = Instant::now();
            writer.write_records(&lines)?;
            if let Some(profile) = sink_profile.as_mut() {
                profile.write_records += write_started.elapsed();
            }
            pb.inc(lines.len() as u64);
            if let Some(ref cb) = config.on_batch_written {
                cb(lines.len(), total_rows, total_input);
            }
        }
    }

    let finish_started = Instant::now();
    writer.finish()?;
    if let Some(profile) = sink_profile.as_mut() {
        profile.writer_finish += finish_started.elapsed();
    }
    pb.finish_and_clear();
    if let Some(profile) = sink_profile {
        eprintln!("{}", profile.summary_line());
    }

    Ok(total_rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_options_json_emits_pick_flags() {
        let config = AnnotateVcfConfig {
            everything: true,
            pick: true,
            pick_allele: true,
            per_gene: true,
            pick_allele_gene: true,
            flag_pick: true,
            flag_pick_allele: true,
            flag_pick_allele_gene: true,
            pick_order: Some("mane_select,tsl,canonical".to_string()),
            ..Default::default()
        };

        let json = config.to_options_json();
        assert!(json.contains("\"everything\":true"));
        assert!(json.contains("\"pick\":true"));
        assert!(json.contains("\"pick_allele\":true"));
        assert!(json.contains("\"per_gene\":true"));
        assert!(json.contains("\"pick_allele_gene\":true"));
        assert!(json.contains("\"flag_pick\":true"));
        assert!(json.contains("\"flag_pick_allele\":true"));
        assert!(json.contains("\"flag_pick_allele_gene\":true"));
        assert!(json.contains("\"pick_order\":\"mane_select,tsl,canonical\""));
    }

    #[test]
    fn test_to_options_json_emits_buffer_size() {
        let config = AnnotateVcfConfig {
            buffer_size: 1234,
            ..Default::default()
        };

        let json = config.to_options_json();
        assert!(json.contains("\"buffer_size\":1234"));
    }

    #[test]
    fn test_default_target_partitions_is_one() {
        assert_eq!(AnnotateVcfConfig::default().target_partitions, 1);
    }

    #[test]
    fn test_to_options_json_does_not_emit_source_selectors() {
        let config = AnnotateVcfConfig {
            refseq: true,
            merged: true,
            ..Default::default()
        };

        let json = config.to_options_json();
        assert!(!json.contains("\"refseq\""));
        assert!(!json.contains("\"merged\""));
    }

    #[test]
    fn test_csq_header_description_matches_vep_pick_layout() {
        let config = AnnotateVcfConfig {
            everything: true,
            flag_pick_allele_gene: true,
            ..Default::default()
        };

        let description = csq_header_description(&config, CacheSourceType::Ensembl);
        assert!(description.starts_with("Consequence annotations from annotate_vep. Format: "));
        assert!(description.contains("|FLAGS|PICK|VARIANT_CLASS|"));
    }

    #[test]
    fn test_csq_header_description_omits_pick_for_filter_modes() {
        let config = AnnotateVcfConfig {
            everything: true,
            pick_allele_gene: true,
            ..Default::default()
        };

        let description = csq_header_description(&config, CacheSourceType::Ensembl);
        assert!(!description.contains("|FLAGS|PICK|VARIANT_CLASS|"));
    }

    #[test]
    fn test_vcf_sink_profile_summary_formats_stage_timings() {
        let mut profile = VcfSinkProfile::default();
        profile.stream_next += std::time::Duration::from_millis(10);
        profile.batch_to_lines += std::time::Duration::from_millis(20);
        profile.format_wait += std::time::Duration::from_millis(25);
        profile.write_records += std::time::Duration::from_millis(30);
        profile.writer_finish += std::time::Duration::from_millis(40);
        profile.batches = 2;
        profile.rows = 3;
        profile.lines = 4;
        profile.body_chunk_bytes = 5;
        profile.format_jobs = 6;
        profile.format_inflight_max = 7;

        let line = profile.summary_line();

        assert!(line.contains("[VEP_PROFILE] vcf_sink_profile"));
        assert!(line.contains("batches=2"));
        assert!(line.contains("rows=3"));
        assert!(line.contains("lines=4"));
        assert!(line.contains("body_chunk_bytes=5"));
        assert!(line.contains("format_jobs=6"));
        assert!(line.contains("format_inflight_max=7"));
        assert!(line.contains("stream_next=0.010s"));
        assert!(line.contains("batch_to_lines=0.020s"));
        assert!(line.contains("format_wait=0.025s"));
        assert!(line.contains("write_records=0.030s"));
        assert!(line.contains("writer_finish=0.040s"));
    }

    #[test]
    fn test_vcf_lines_to_body_chunk_appends_record_newlines() {
        let (chunk, lines) = vcf_lines_to_body_chunk(vec![
            VcfRecordLine {
                line: "chr1\t1\t.\tA\tC\t.\t.\t.".to_string(),
            },
            VcfRecordLine {
                line: "chr1\t2\t.\tG\tT\t.\t.\t.".to_string(),
            },
        ]);

        assert_eq!(lines, 2);
        assert_eq!(
            chunk,
            b"chr1\t1\t.\tA\tC\t.\t.\t.\nchr1\t2\t.\tG\tT\t.\t.\t.\n"
        );
    }
}
