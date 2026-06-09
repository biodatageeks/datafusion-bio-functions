//! VCF output sink for annotate_vep() results.
//!
//! Provides [`annotate_to_vcf`] — a single-call function that reads a VCF,
//! annotates it, and streams results to an output VCF file.

use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
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
use crate::pipeline_trace::{self, PipelineTraceValue as TraceValue};

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
    contig_partitions: usize,
    contig_inflight_max: usize,
}

impl VcfSinkProfile {
    fn summary_line(&self) -> String {
        format!(
            "[VEP_PROFILE] vcf_sink_profile batches={} rows={} lines={} body_chunk_bytes={} format_jobs={} format_inflight_max={} contig_partitions={} contig_inflight_max={} stream_next={:.3}s batch_to_lines={:.3}s format_wait={:.3}s write_records={:.3}s writer_finish={:.3}s",
            self.batches,
            self.rows,
            self.lines,
            self.body_chunk_bytes,
            self.format_jobs,
            self.format_inflight_max,
            self.contig_partitions,
            self.contig_inflight_max,
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

#[derive(Debug)]
struct PartitionedVcfBody {
    partition_id: usize,
    path: PathBuf,
    batches: usize,
    input_rows: usize,
    lines: usize,
    bytes: usize,
    stream_next: Duration,
    format_duration: Duration,
    temp_write_duration: Duration,
}

struct VcfBodyTempDir {
    path: PathBuf,
}

impl VcfBodyTempDir {
    fn new() -> Result<Self> {
        let unique = format!(
            "vep_vcf_body_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|duration| duration.as_nanos())
                .unwrap_or(0)
        );
        let path = std::env::temp_dir().join(unique);
        std::fs::create_dir(&path).map_err(|e| {
            DataFusionError::Execution(format!(
                "Failed to create VCF body shard tempdir {}: {e}",
                path.display()
            ))
        })?;
        Ok(Self { path })
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for VcfBodyTempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct VepConcurrencyPlan {
    lookup_partitions: usize,
    chromosome_lanes: usize,
    inline_lookup: bool,
    spawn_vcf_provider_open: bool,
}

impl VepConcurrencyPlan {
    fn from_config(config: &AnnotateVcfConfig) -> Self {
        Self::from_forks(config.forks.unwrap_or(0), config.lookup_partitions.max(1))
    }

    fn from_forks(forks: usize, lookup_partitions: usize) -> Self {
        if forks == 0 {
            return Self {
                lookup_partitions: 1,
                chromosome_lanes: 1,
                inline_lookup: true,
                spawn_vcf_provider_open: false,
            };
        }

        let chromosome_lanes = forks.max(1);
        Self {
            lookup_partitions,
            chromosome_lanes,
            inline_lookup: false,
            spawn_vcf_provider_open: true,
        }
    }
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
    pipeline_trace::emit(
        "vcf_format",
        "start",
        &[
            ("batch_id", TraceValue::Usize(batch_id)),
            ("rows", TraceValue::Usize(input_rows)),
        ],
    );
    let lines = batch_to_vcf_lines(
        &batch,
        vcf_info_fields.as_slice(),
        unique_format_tags.as_slice(),
        sample_names.as_slice(),
        coordinate_zero_based,
    )?;
    let (bytes, line_count) = vcf_lines_to_body_chunk(lines);
    pipeline_trace::emit(
        "vcf_format",
        "done",
        &[
            ("batch_id", TraceValue::Usize(batch_id)),
            ("rows", TraceValue::Usize(input_rows)),
            ("lines", TraceValue::Usize(line_count)),
            ("bytes", TraceValue::Usize(bytes.len())),
            ("elapsed", TraceValue::Duration(format_started.elapsed())),
        ],
    );
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
        let write_elapsed = write_started.elapsed();
        pipeline_trace::emit(
            "vcf_write",
            "done",
            &[
                ("batch_id", TraceValue::Usize(chunk.batch_id)),
                ("rows", TraceValue::Usize(chunk.input_rows)),
                ("lines", TraceValue::Usize(chunk.lines)),
                ("bytes", TraceValue::Usize(chunk.bytes.len())),
                ("elapsed", TraceValue::Duration(write_elapsed)),
            ],
        );
        if let Some(profile) = sink_profile.as_mut() {
            profile.write_records += write_elapsed;
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

fn copy_body_file_to_writer(path: &Path, writer: &mut VcfLocalWriter) -> Result<Duration> {
    let started = Instant::now();
    let mut file = std::fs::File::open(path).map_err(|e| {
        DataFusionError::Execution(format!(
            "Failed to open VCF body shard {}: {e}",
            path.display()
        ))
    })?;
    let mut buffer = vec![0u8; 8 * 1024 * 1024];
    loop {
        let bytes = file.read(&mut buffer).map_err(|e| {
            DataFusionError::Execution(format!(
                "Failed to read VCF body shard {}: {e}",
                path.display()
            ))
        })?;
        if bytes == 0 {
            break;
        }
        write_vcf_body_chunk(writer, &buffer[..bytes])?;
    }
    Ok(started.elapsed())
}

fn partition_body_job_limit(partition_count: usize, max_contig_jobs: usize) -> usize {
    max_contig_jobs.max(1).min(partition_count.max(1))
}

async fn write_vcf_partition_body(
    partition_id: usize,
    mut stream: datafusion::physical_plan::SendableRecordBatchStream,
    path: PathBuf,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
) -> Result<PartitionedVcfBody> {
    let file = std::fs::File::create(&path).map_err(|e| {
        DataFusionError::Execution(format!(
            "Failed to create VCF body shard {}: {e}",
            path.display()
        ))
    })?;
    let mut writer = std::io::BufWriter::new(file);
    let mut batch_id = 0usize;
    let mut result = PartitionedVcfBody {
        partition_id,
        path,
        batches: 0,
        input_rows: 0,
        lines: 0,
        bytes: 0,
        stream_next: Duration::ZERO,
        format_duration: Duration::ZERO,
        temp_write_duration: Duration::ZERO,
    };

    use futures::StreamExt;
    loop {
        let stream_started = Instant::now();
        let next_batch = stream.next().await;
        result.stream_next += stream_started.elapsed();
        let Some(batch_result) = next_batch else {
            break;
        };
        let batch = batch_result?;
        let input_rows = batch.num_rows();
        let vcf_info_fields = Arc::clone(&vcf_info_fields);
        let unique_format_tags = Arc::clone(&unique_format_tags);
        let sample_names = Arc::clone(&sample_names);
        let formatted = tokio::task::spawn_blocking(move || {
            format_vcf_body_chunk(
                batch_id,
                batch,
                vcf_info_fields,
                unique_format_tags,
                sample_names,
                coordinate_zero_based,
            )
        })
        .await
        .map_err(|e| DataFusionError::Execution(format!("VCF formatter task failed: {e}")))??;
        batch_id += 1;

        let write_started = Instant::now();
        writer.write_all(&formatted.bytes).map_err(|e| {
            DataFusionError::Execution(format!(
                "Failed to write VCF body shard {}: {e}",
                result.path.display()
            ))
        })?;
        result.temp_write_duration += write_started.elapsed();
        result.batches += 1;
        result.input_rows += input_rows;
        result.lines += formatted.lines;
        result.bytes += formatted.bytes.len();
        result.format_duration += formatted.format_duration;
    }

    let write_started = Instant::now();
    writer.flush().map_err(|e| {
        DataFusionError::Execution(format!(
            "Failed to flush VCF body shard {}: {e}",
            result.path.display()
        ))
    })?;
    result.temp_write_duration += write_started.elapsed();
    Ok(result)
}

async fn schedule_partition_body_jobs(
    streams: Vec<(usize, datafusion::physical_plan::SendableRecordBatchStream)>,
    tempdir_path: PathBuf,
    max_jobs: usize,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
    completed: tokio::sync::mpsc::UnboundedSender<PartitionedVcfBody>,
) -> Result<usize> {
    let mut stream_iter = streams.into_iter();
    let mut jobs = tokio::task::JoinSet::new();
    let mut inflight_max = 0usize;
    let max_jobs = max_jobs.max(1);

    loop {
        while jobs.len() < max_jobs {
            let Some((partition_id, stream)) = stream_iter.next() else {
                break;
            };
            let path = tempdir_path.join(format!("partition_{partition_id:04}.vcf.body"));
            let vcf_info_fields = Arc::clone(&vcf_info_fields);
            let unique_format_tags = Arc::clone(&unique_format_tags);
            let sample_names = Arc::clone(&sample_names);
            jobs.spawn(write_vcf_partition_body(
                partition_id,
                stream,
                path,
                vcf_info_fields,
                unique_format_tags,
                sample_names,
                coordinate_zero_based,
            ));
            inflight_max = inflight_max.max(jobs.len());
        }

        if jobs.is_empty() {
            break;
        }

        let Some(joined) = jobs.join_next().await else {
            continue;
        };
        let body = joined.map_err(|e| {
            DataFusionError::Execution(format!("VCF body shard task failed: {e}"))
        })??;
        if completed.send(body).is_err() {
            break;
        }
    }

    Ok(inflight_max)
}

fn spawn_partition_body_scheduler(
    streams: Vec<(usize, datafusion::physical_plan::SendableRecordBatchStream)>,
    tempdir_path: PathBuf,
    max_jobs: usize,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
) -> (
    tokio::sync::mpsc::UnboundedReceiver<PartitionedVcfBody>,
    tokio::task::JoinHandle<Result<usize>>,
) {
    let (completed_tx, completed_rx) = tokio::sync::mpsc::unbounded_channel();
    let handle = tokio::spawn(schedule_partition_body_jobs(
        streams,
        tempdir_path,
        max_jobs,
        vcf_info_fields,
        unique_format_tags,
        sample_names,
        coordinate_zero_based,
        completed_tx,
    ));
    (completed_rx, handle)
}

#[allow(clippy::too_many_arguments)]
async fn stream_vcf_partition_to_writer(
    partition_id: usize,
    mut stream: datafusion::physical_plan::SendableRecordBatchStream,
    writer: &mut VcfLocalWriter,
    pb: &ProgressBar,
    config: &AnnotateVcfConfig,
    total_input: usize,
    total_rows: &mut usize,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
    sink_profile: &mut Option<VcfSinkProfile>,
) -> Result<()> {
    let mut batch_id = 0usize;

    use futures::StreamExt;
    loop {
        let stream_started = Instant::now();
        let next_batch = stream.next().await;
        let stream_elapsed = stream_started.elapsed();
        if let Some(profile) = sink_profile.as_mut() {
            profile.stream_next += stream_elapsed;
        }
        let Some(batch_result) = next_batch else {
            break;
        };

        let batch = batch_result?;
        let input_rows = batch.num_rows();
        let vcf_info_fields = Arc::clone(&vcf_info_fields);
        let unique_format_tags = Arc::clone(&unique_format_tags);
        let sample_names = Arc::clone(&sample_names);
        let formatted = tokio::task::spawn_blocking(move || {
            format_vcf_body_chunk(
                batch_id,
                batch,
                vcf_info_fields,
                unique_format_tags,
                sample_names,
                coordinate_zero_based,
            )
        })
        .await
        .map_err(|e| DataFusionError::Execution(format!("VCF formatter task failed: {e}")))??;
        batch_id += 1;

        let write_started = Instant::now();
        write_vcf_body_chunk(writer, &formatted.bytes)?;
        let write_elapsed = write_started.elapsed();
        pipeline_trace::emit(
            "vcf_write",
            "done",
            &[
                ("partition_id", TraceValue::Usize(partition_id)),
                ("batch_id", TraceValue::Usize(formatted.batch_id)),
                ("rows", TraceValue::Usize(input_rows)),
                ("lines", TraceValue::Usize(formatted.lines)),
                ("bytes", TraceValue::Usize(formatted.bytes.len())),
                ("elapsed", TraceValue::Duration(write_elapsed)),
            ],
        );

        if let Some(profile) = sink_profile.as_mut() {
            profile.batch_to_lines += formatted.format_duration;
            profile.write_records += write_elapsed;
            profile.batches += 1;
            profile.rows += input_rows;
            profile.lines += formatted.lines;
            profile.body_chunk_bytes += formatted.bytes.len();
            profile.format_jobs += 1;
        }
        *total_rows += formatted.lines;
        pb.inc(formatted.lines as u64);
        if let Some(ref cb) = config.on_batch_written {
            cb(formatted.lines, *total_rows, total_input);
        }
    }

    Ok(())
}

async fn stream_partitioned_vcf_body(
    df: datafusion::dataframe::DataFrame,
    writer: &mut VcfLocalWriter,
    pb: &ProgressBar,
    config: &AnnotateVcfConfig,
    total_input: usize,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
    max_contig_jobs: usize,
    sink_profile: &mut Option<VcfSinkProfile>,
) -> Result<usize> {
    let streams = df.execute_stream_partitioned().await?;
    let partition_count = streams.len();
    if partition_count <= 1 {
        if let Some(profile) = sink_profile.as_mut() {
            profile.contig_partitions = partition_count;
        }

        let Some(stream) = streams.into_iter().next() else {
            return Ok(0);
        };
        let mut total_rows = 0usize;
        stream_vcf_partition_to_writer(
            0,
            stream,
            writer,
            pb,
            config,
            total_input,
            &mut total_rows,
            vcf_info_fields,
            unique_format_tags,
            sample_names,
            coordinate_zero_based,
            sink_profile,
        )
        .await?;
        return Ok(total_rows);
    }

    let tempdir = VcfBodyTempDir::new()?;
    let max_contig_jobs = partition_body_job_limit(partition_count, max_contig_jobs);
    if let Some(profile) = sink_profile.as_mut() {
        profile.contig_partitions = partition_count;
    }

    let streams = streams.into_iter().enumerate().collect::<Vec<_>>();
    let mut total_rows = 0usize;
    let tempdir_path = tempdir.path().to_path_buf();

    let (mut body_rx, scheduler_handle) = spawn_partition_body_scheduler(
        streams,
        tempdir_path,
        max_contig_jobs,
        Arc::clone(&vcf_info_fields),
        Arc::clone(&unique_format_tags),
        Arc::clone(&sample_names),
        coordinate_zero_based,
    );

    let mut scheduler_handle = Some(scheduler_handle);
    let mut ready: BTreeMap<usize, PartitionedVcfBody> = BTreeMap::new();
    let mut next_write_partition_id = 0usize;
    while next_write_partition_id < partition_count {
        if let Some(body) = ready.remove(&next_write_partition_id) {
            let copy_duration = copy_body_file_to_writer(&body.path, writer)?;
            let _ = std::fs::remove_file(&body.path);
            if let Some(profile) = sink_profile.as_mut() {
                profile.stream_next += body.stream_next;
                profile.batch_to_lines += body.format_duration;
                profile.write_records += body.temp_write_duration + copy_duration;
                profile.batches += body.batches;
                profile.rows += body.input_rows;
                profile.lines += body.lines;
                profile.body_chunk_bytes += body.bytes;
                profile.format_jobs += body.batches;
            }
            total_rows += body.lines;
            pb.inc(body.lines as u64);
            if let Some(ref cb) = config.on_batch_written {
                cb(body.lines, total_rows, total_input);
            }
            next_write_partition_id += 1;
            continue;
        }

        let Some(body) = body_rx.recv().await else {
            if let Some(handle) = scheduler_handle.take() {
                handle.await.map_err(|e| {
                    DataFusionError::Execution(format!("VCF body shard scheduler task failed: {e}"))
                })??;
            }
            return Err(DataFusionError::Execution(format!(
                "partitioned VCF body writer missing partition {next_write_partition_id}"
            )));
        };
        ready.insert(body.partition_id, body);
    }

    let contig_inflight_max = if let Some(handle) = scheduler_handle.take() {
        handle.await.map_err(|e| {
            DataFusionError::Execution(format!("VCF body shard scheduler task failed: {e}"))
        })??
    } else {
        0
    };

    if let Some(profile) = sink_profile.as_mut() {
        profile.contig_inflight_max = profile.contig_inflight_max.max(contig_inflight_max);
    }

    Ok(total_rows)
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
    /// User-facing chromosome lane count. Missing or `Some(0)` is a strict
    /// single-lane path with no helper lookup task; `Some(n)` runs up to `n`
    /// contigs concurrently.
    pub forks: Option<usize>,
    /// DataFusion partitions used by the per-contig variation lookup.
    pub lookup_partitions: usize,
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
            forks: Some(0),
            lookup_partitions: 1,
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
            .field("forks", &self.forks)
            .field("lookup_partitions", &self.lookup_partitions)
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
        if let Some(forks) = self.forks {
            let worker_forks = if forks == 0 {
                0
            } else {
                self.lookup_partitions.max(1)
            };
            opts.insert(
                "forks".into(),
                serde_json::Value::Number(serde_json::Number::from(worker_forks)),
            );
            opts.insert("inline_lookup".into(), serde_json::Value::Bool(forks == 0));
        }
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
    let concurrency_plan = VepConcurrencyPlan::from_config(config);
    if sink_profile_enabled() {
        eprintln!(
            "[VEP_PROFILE] concurrency_plan lookup_partitions={} chromosome_lanes={} inline_lookup={} spawn_vcf_provider_open={}",
            concurrency_plan.lookup_partitions,
            concurrency_plan.chromosome_lanes,
            concurrency_plan.inline_lookup,
            concurrency_plan.spawn_vcf_provider_open
        );
    }

    // 1. Create session and register VCF table.
    let session_config = datafusion::prelude::SessionConfig::new()
        .with_target_partitions(concurrency_plan.lookup_partitions);
    let ctx = SessionContext::new_with_config(session_config);
    crate::register_vep_functions(&ctx);

    let vcf_path = input_vcf.to_string();
    let vcf_provider = if concurrency_plan.spawn_vcf_provider_open {
        tokio::task::spawn_blocking(move || {
            VcfTableProvider::new(vcf_path, None, None, None, false)
        })
        .await
        .map_err(|e| datafusion::common::DataFusionError::External(Box::new(e)))??
    } else {
        VcfTableProvider::new(vcf_path, None, None, None, false)?
    };

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

    let mut total_rows = 0;
    let mut sink_profile = sink_profile_enabled().then(VcfSinkProfile::default);

    if concurrency_plan.chromosome_lanes > 1 {
        total_rows = stream_partitioned_vcf_body(
            df,
            &mut writer,
            &pb,
            config,
            total_input,
            Arc::new(vcf_info_fields),
            Arc::new(unique_format_tags),
            Arc::new(sample_names),
            coordinate_zero_based,
            concurrency_plan.chromosome_lanes,
            &mut sink_profile,
        )
        .await?;
    } else {
        use futures::StreamExt;
        let mut stream = df.execute_stream().await?;
        let mut next_serial_batch_id = 0usize;
        loop {
            let stream_started = Instant::now();
            let next_batch = stream.next().await;
            let stream_elapsed = stream_started.elapsed();
            if let Some(profile) = sink_profile.as_mut() {
                profile.stream_next += stream_elapsed;
            }
            let Some(batch_result) = next_batch else {
                break;
            };
            let batch = batch_result?;
            let input_rows = batch.num_rows();
            let batch_id = next_serial_batch_id;
            next_serial_batch_id += 1;
            pipeline_trace::emit(
                "vcf_stream",
                "batch_ready",
                &[
                    ("batch_id", TraceValue::Usize(batch_id)),
                    ("rows", TraceValue::Usize(input_rows)),
                    ("stream_wait", TraceValue::Duration(stream_elapsed)),
                ],
            );
            let lines_started = Instant::now();
            pipeline_trace::emit(
                "vcf_format",
                "start",
                &[
                    ("batch_id", TraceValue::Usize(batch_id)),
                    ("rows", TraceValue::Usize(input_rows)),
                ],
            );
            let lines = batch_to_vcf_lines(
                &batch,
                &vcf_info_fields,
                &unique_format_tags,
                &sample_names,
                coordinate_zero_based,
            )?;
            let format_elapsed = lines_started.elapsed();
            pipeline_trace::emit(
                "vcf_format",
                "done",
                &[
                    ("batch_id", TraceValue::Usize(batch_id)),
                    ("rows", TraceValue::Usize(input_rows)),
                    ("lines", TraceValue::Usize(lines.len())),
                    ("elapsed", TraceValue::Duration(format_elapsed)),
                ],
            );
            if let Some(profile) = sink_profile.as_mut() {
                profile.batch_to_lines += format_elapsed;
                profile.batches += 1;
                profile.rows += input_rows;
                profile.lines += lines.len();
            }
            total_rows += lines.len();
            let write_started = Instant::now();
            writer.write_records(&lines)?;
            let write_elapsed = write_started.elapsed();
            pipeline_trace::emit(
                "vcf_write",
                "done",
                &[
                    ("batch_id", TraceValue::Usize(batch_id)),
                    ("rows", TraceValue::Usize(input_rows)),
                    ("lines", TraceValue::Usize(lines.len())),
                    ("elapsed", TraceValue::Duration(write_elapsed)),
                ],
            );
            if let Some(profile) = sink_profile.as_mut() {
                profile.write_records += write_elapsed;
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
    use std::pin::Pin;
    use std::task::{Context, Poll};

    use datafusion::arrow::datatypes::{Schema, SchemaRef};
    use datafusion::execution::RecordBatchStream;
    use datafusion::physical_plan::SendableRecordBatchStream;

    struct EmptyBatchStream {
        schema: SchemaRef,
    }

    impl futures::Stream for EmptyBatchStream {
        type Item = Result<RecordBatch>;

        fn poll_next(self: Pin<&mut Self>, _cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
            Poll::Ready(None)
        }
    }

    impl RecordBatchStream for EmptyBatchStream {
        fn schema(&self) -> SchemaRef {
            self.schema.clone()
        }
    }

    fn empty_batch_stream() -> SendableRecordBatchStream {
        Box::pin(EmptyBatchStream {
            schema: Arc::new(Schema::empty()),
        })
    }

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
    fn test_forks_zero_derives_strict_inline_concurrency_plan() {
        let config = AnnotateVcfConfig {
            forks: Some(0),
            ..Default::default()
        };

        let plan = VepConcurrencyPlan::from_config(&config);

        assert_eq!(plan.lookup_partitions, 1);
        assert_eq!(plan.chromosome_lanes, 1);
        assert!(plan.inline_lookup);
        assert!(!plan.spawn_vcf_provider_open);
    }

    #[test]
    fn test_forks_one_preserves_legacy_single_partition_plan() {
        let config = AnnotateVcfConfig {
            forks: Some(1),
            ..Default::default()
        };

        let plan = VepConcurrencyPlan::from_config(&config);

        assert_eq!(plan.lookup_partitions, 1);
        assert_eq!(plan.chromosome_lanes, 1);
        assert!(!plan.inline_lookup);
        assert!(plan.spawn_vcf_provider_open);
    }

    #[test]
    fn test_forks_nonzero_uses_fork_count_for_chromosome_lanes() {
        let config = AnnotateVcfConfig {
            forks: Some(4),
            lookup_partitions: 3,
            ..Default::default()
        };

        let plan = VepConcurrencyPlan::from_config(&config);

        assert_eq!(plan.lookup_partitions, 3);
        assert_eq!(plan.chromosome_lanes, 4);
        assert!(!plan.inline_lookup);
        assert!(plan.spawn_vcf_provider_open);
    }

    #[test]
    fn test_to_options_json_emits_forks_plan_when_set() {
        let config = AnnotateVcfConfig {
            forks: Some(0),
            ..Default::default()
        };

        let json = config.to_options_json();
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(value["forks"], 0);
        assert_eq!(value["inline_lookup"], true);
    }

    #[test]
    fn test_to_options_json_emits_only_forks_for_chromosome_lanes() {
        let config = AnnotateVcfConfig {
            forks: Some(8),
            lookup_partitions: 4,
            ..Default::default()
        };

        let json = config.to_options_json();
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(value["forks"], 4);
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
        profile.contig_partitions = 8;
        profile.contig_inflight_max = 9;

        let line = profile.summary_line();

        assert!(line.contains("[VEP_PROFILE] vcf_sink_profile"));
        assert!(line.contains("batches=2"));
        assert!(line.contains("rows=3"));
        assert!(line.contains("lines=4"));
        assert!(line.contains("body_chunk_bytes=5"));
        assert!(line.contains("format_jobs=6"));
        assert!(line.contains("format_inflight_max=7"));
        assert!(line.contains("contig_partitions=8"));
        assert!(line.contains("contig_inflight_max=9"));
        assert!(line.contains("stream_next=0.010s"));
        assert!(line.contains("batch_to_lines=0.020s"));
        assert!(line.contains("format_wait=0.025s"));
        assert!(line.contains("write_records=0.030s"));
        assert!(line.contains("writer_finish=0.040s"));
    }

    #[test]
    fn test_partitioned_body_writer_schedules_partition_zero_as_a_job() {
        assert_eq!(partition_body_job_limit(22, 12), 12);
        assert_eq!(partition_body_job_limit(3, 12), 3);
        assert_eq!(partition_body_job_limit(22, 1), 1);
        assert_eq!(partition_body_job_limit(0, 12), 1);
    }

    #[tokio::test]
    async fn test_partition_body_scheduler_uses_full_job_limit() {
        let tempdir = tempfile::tempdir().unwrap();
        let streams = (0..5)
            .map(|partition_id| (partition_id, empty_batch_stream()))
            .collect::<Vec<_>>();
        let (completed_tx, mut completed_rx) = tokio::sync::mpsc::unbounded_channel();

        let inflight_max = schedule_partition_body_jobs(
            streams,
            tempdir.path().to_path_buf(),
            2,
            Arc::new(Vec::new()),
            Arc::new(Vec::new()),
            Arc::new(Vec::new()),
            false,
            completed_tx,
        )
        .await
        .unwrap();

        let mut completed_ids = Vec::new();
        while let Some(body) = completed_rx.recv().await {
            completed_ids.push(body.partition_id);
        }
        completed_ids.sort_unstable();

        assert_eq!(inflight_max, 2);
        assert_eq!(completed_ids, vec![0, 1, 2, 3, 4]);
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
