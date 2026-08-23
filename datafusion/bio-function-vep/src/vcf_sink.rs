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
use datafusion_bio_format_core::metadata::VCF_HEADER_RAW_LINES_KEY;
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
    spawn_vcf_provider_open: bool,
}

impl VepConcurrencyPlan {
    fn from_config(config: &AnnotateVcfConfig) -> Self {
        // Single `workers` knob: N position-ordered lookup partitions feeding
        // N independent annotation range pipelines, one contig at a time.
        // Correctness for workers>1 relies on the indexed-input guard in
        // annotate_to_vcf (contiguous, position-ordered partitions).
        let lookup_partitions = config.workers.max(1);
        Self {
            lookup_partitions,
            spawn_vcf_provider_open: lookup_partitions > 1,
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

/// Published by `ContigAnnotationStream` the instant a contig's shard workers
/// all finish (shards flushed). The assembler thread appends the id range
/// `[first_shard_id, end_shard_id)` — which is contiguous and equals final-VCF
/// position order — and adds `output_lines` to the running body-line total.
pub(crate) struct ContigShardRange {
    pub(crate) first_shard_id: usize,
    pub(crate) end_shard_id: usize,
    pub(crate) output_lines: usize,
}

/// Formatter context + tempdir threaded into the `threads>1` engine workers so
/// each fused worker can serialize its annotated batches directly into its own
/// VCF body shard (`format_vcf_body_chunk` + `BufWriter<File>`), removing the
/// ordered drain / output channel. Derived in `annotate_to_vcf` from the VCF
/// schema (independent of the SQL path) and injected via
/// `AnnotateProvider` → `ContigAnnotationConfig`.
#[derive(Clone)]
pub(crate) struct VcfShardContext {
    pub(crate) vcf_info_fields: Arc<Vec<String>>,
    pub(crate) unique_format_tags: Arc<Vec<String>>,
    pub(crate) sample_names: Arc<Vec<String>>,
    pub(crate) coordinate_zero_based: bool,
    pub(crate) tempdir: PathBuf,
    /// Output compression for the per-worker shard files. Only ever `Bgzf`
    /// (bgzf output → parallel per-worker compression + raw block concat) or
    /// `Plain` (plain/gzip output → text shards).
    pub(crate) shard_compression: VcfCompressionType,
    /// Running count of input rows annotated across all shard workers. Polled by
    /// `drive_sharded_vcf_annotation` to drive the progress bar / callback
    /// incrementally during annotation (shards are written concurrently and only
    /// concatenated at the end, so this is the only live progress signal).
    pub(crate) rows_done: Arc<std::sync::atomic::AtomicUsize>,
    /// One message per completed contig, in completion (= ascending id) order.
    /// Dropped when annotation finishes, which signals the assembler to finalize.
    pub(crate) contig_done_tx: std::sync::mpsc::Sender<ContigShardRange>,
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

/// Reusable streaming VCF body shard writer: owns a `BufWriter<File>` + the
/// per-batch formatter context, and streams formatted batches into one shard
/// file with no buffering of the batch. Used by the `workers>1` engine workers
/// (which call the sync [`write_batch`](Self::write_batch) inside
/// `block_in_place`).
pub(crate) struct VcfBodyShardWriter {
    /// Owns the shard's output sink. `Plain` for text shards (plain/gzip output),
    /// `Bgzf` for bgzf output — workers compress their own shard in parallel; the
    /// assembler then does a raw block concat. Reuses the bio-formats writer enum
    /// so no direct `noodles-bgzf` dependency is needed here.
    writer: VcfLocalWriter,
    vcf_info_fields: Arc<Vec<String>>,
    unique_format_tags: Arc<Vec<String>>,
    sample_names: Arc<Vec<String>>,
    coordinate_zero_based: bool,
    batch_id: usize,
    pub(crate) lines: usize,
    pub(crate) input_rows: usize,
    pub(crate) bytes: usize,
}

impl VcfBodyShardWriter {
    pub(crate) fn create(
        path: &Path,
        vcf_info_fields: Arc<Vec<String>>,
        unique_format_tags: Arc<Vec<String>>,
        sample_names: Arc<Vec<String>>,
        coordinate_zero_based: bool,
        compression: VcfCompressionType,
    ) -> Result<Self> {
        let writer = VcfLocalWriter::with_compression(path, compression)?;
        Ok(Self {
            writer,
            vcf_info_fields,
            unique_format_tags,
            sample_names,
            coordinate_zero_based,
            batch_id: 0,
            lines: 0,
            input_rows: 0,
            bytes: 0,
        })
    }

    /// Format one batch with the owned context and stream it to the shard file.
    /// Synchronous: callers must be on a thread where blocking is acceptable
    /// (the engine workers run inside `block_in_place`). No batch buffering.
    pub(crate) fn write_batch(&mut self, batch: RecordBatch) -> Result<()> {
        let formatted = format_vcf_body_chunk(
            self.batch_id,
            batch,
            Arc::clone(&self.vcf_info_fields),
            Arc::clone(&self.unique_format_tags),
            Arc::clone(&self.sample_names),
            self.coordinate_zero_based,
        )?;
        self.batch_id += 1;
        self.write_formatted(&formatted)?;
        Ok(())
    }

    /// Allocate the next batch id (for callers that format off-thread).
    fn next_batch_id(&mut self) -> usize {
        let id = self.batch_id;
        self.batch_id += 1;
        id
    }

    /// Write already-formatted bytes and accumulate counters. Returns the IO
    /// duration so async callers can attribute write time. Used by
    /// [`write_vcf_partition_body`] which formats via `spawn_blocking`.
    fn write_formatted(&mut self, formatted: &FormattedVcfBatch) -> Result<Duration> {
        let started = Instant::now();
        write_vcf_body_chunk(&mut self.writer, &formatted.bytes)?;
        let elapsed = started.elapsed();
        self.lines += formatted.lines;
        self.input_rows += formatted.input_rows;
        self.bytes += formatted.bytes.len();
        Ok(elapsed)
    }

    pub(crate) fn finish(self) -> Result<()> {
        // Consumes the writer: flushes (Plain) or writes the final BGZF block +
        // 28-byte EOF marker (Bgzf). The assembler strips that trailing EOF when
        // concatenating bgzf shards.
        self.writer.finish()
    }
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
    /// Optional assertion against embedded per-shard cache version metadata.
    pub expected_cache_version: Option<String>,
    /// Path to indexed reference FASTA (required for `everything` / `hgvs`).
    pub reference_fasta_path: Option<String>,
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
    /// Maximum number of independent cold-Parquet lookup shards per contig.
    pub target_partitions: usize,
    /// Number of parallel fused window-annotation pipelines within a contig
    /// (`1` = serial). The single annotation-concurrency knob.
    pub workers: usize,
    /// Output compression type.
    pub compression: VcfCompressionType,
    /// Show an indicatif progress bar on stderr (for Rust CLI).
    /// For Python/Jupyter, use `on_batch_written` with tqdm instead.
    pub show_progress: bool,
    /// Optional callback invoked after each batch is written.
    /// Used by Python wrappers to drive tqdm progress bars that work in Jupyter.
    pub on_batch_written: Option<OnBatchWritten>,
    /// Root of a custom-plugin cache (`<root>/plugin/*/manifest.json`). When
    /// `Some`, plugin CSQ fields are appended to output. `None` = disabled
    /// (byte-identical to no-plugin output).
    pub plugin_cache_root: Option<std::path::PathBuf>,
    /// Tool name recorded in the output header's provenance lines, e.g.
    /// `"vepyr"`. `None` records this engine's own crate name. The value only
    /// labels provenance; it never claims the output came from Ensembl VEP.
    pub provenance_tool_name: Option<String>,
    /// Version recorded next to [`Self::provenance_tool_name`]. `None` records
    /// this engine's crate version.
    pub provenance_tool_version: Option<String>,
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
            expected_cache_version: None,
            reference_fasta_path: None,
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
            workers: 1,
            compression: VcfCompressionType::Plain,
            show_progress: false,
            on_batch_written: None,
            plugin_cache_root: None,
            provenance_tool_name: None,
            provenance_tool_version: None,
        }
    }
}

impl std::fmt::Debug for AnnotateVcfConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AnnotateVcfConfig")
            .field("everything", &self.everything)
            .field("buffer_size", &self.buffer_size)
            .field("workers", &self.workers)
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
        if let Some(ref expected_cache_version) = self.expected_cache_version {
            opts.insert(
                "expected_cache_version".into(),
                serde_json::Value::String(expected_cache_version.clone()),
            );
        }
        if let Some(ref fasta) = self.reference_fasta_path {
            opts.insert(
                "reference_fasta_path".into(),
                serde_json::Value::String(fasta.clone()),
            );
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
        opts.insert(
            "target_partitions".into(),
            serde_json::Value::Number(serde_json::Number::from(self.target_partitions.max(1))),
        );
        opts.insert(
            "workers".into(),
            serde_json::Value::Number(serde_json::Number::from(self.workers.max(1))),
        );
        // Custom-plugin cache root flows to the annotate_vep UDTF (workers=1 / SQL
        // path) via options_json, mirroring every other config field. The sharded
        // (workers>1) path passes it via `with_plugin_cache_root` instead.
        if let Some(ref root) = self.plugin_cache_root {
            opts.insert(
                "plugin_cache_root".into(),
                serde_json::Value::String(root.to_string_lossy().into_owned()),
            );
        }
        serde_json::to_string(&serde_json::Value::Object(opts)).unwrap()
    }

    fn to_options_json_with_backend(&self, backend: &str) -> String {
        let mut value: serde_json::Value =
            serde_json::from_str(&self.to_options_json()).expect("generated options JSON is valid");
        if let Some(opts) = value.as_object_mut() {
            opts.insert(
                "cache_format".into(),
                serde_json::Value::String(cache_format_for_backend(backend).to_string()),
            );
        }
        value.to_string()
    }

    fn include_pick_output(&self) -> bool {
        self.flag_pick || self.flag_pick_allele || self.flag_pick_allele_gene
    }
}

fn cache_format_for_backend(_backend: &str) -> &str {
    // Parquet is the only supported cache format.
    "parquet"
}

fn cache_source_type_from_cache_source_for_backend(
    cache_source: &str,
    _backend: &str,
) -> Result<CacheSourceType> {
    #[cfg(feature = "parquet-cache")]
    {
        CacheSourceType::from_partitioned_parquet_cache_source(cache_source)
    }
    #[cfg(not(feature = "parquet-cache"))]
    {
        let _ = cache_source;
        Err(DataFusionError::Plan(
            "annotate_to_vcf(): Parquet cache source metadata requires the parquet-cache feature"
                .to_string(),
        ))
    }
}

/// Header lines recording how this run was configured, so its output can be
/// reproduced.
///
/// Deliberately excludes a wall-clock timestamp. Ensembl VEP writes one into
/// its `##VEP=` line, which makes VEP's own output differ byte for byte between
/// two identical runs; omitting it keeps annotation output byte-reproducible.
///
/// Deliberately does not write a `##VEP=` line. That would assert Ensembl VEP
/// produced the file, at a time it did not run, from a path that may not exist
/// — a false entry in an audit trail. It would not buy byte parity either,
/// because VEP's timestamp differs every run.
///
/// The Ensembl cache's data-source versions (gnomAD, ClinVar, dbSNP, GENCODE,
/// SIFT, PolyPhen …) determine what is inside `CSQ` and belong here, but are
/// dropped when a raw cache is converted, so they are not yet available at
/// annotation time. Tracked separately.
/// Marks this engine's own identity inside a provenance line, independent of the
/// configurable tool name, so a run can recognise provenance an earlier run wrote
/// under a different label.
const ENGINE_MARKER_ATTR: &str = concat!("engine=\"", env!("CARGO_PKG_NAME"), " ");
const ENGINE_MARKER_JSON: &str = concat!("\"engine\":\"", env!("CARGO_PKG_NAME"), "\"");

/// True when `line` is provenance describing a `CSQ` field that is about to be
/// replaced, and so must not be carried into the new header.
///
/// Matches three things:
///
/// - Ensembl VEP's own two lines, by exact key. `##VEP` as a bare prefix would
///   also swallow unrelated meta-information keys such as `##VEPStatus`, which
///   are valid VCF and belong to the source header.
/// - Provenance written under the currently configured tool name.
/// - Provenance written under *any* earlier tool name. A file annotated with the
///   engine default and then re-annotated by a wrapper that sets
///   `provenance_tool_name` would otherwise keep the first run's lines.
///
/// The third case is recognised by the *shape* [`provenance_header_lines`]
/// emits, not by a bare substring search: an identity line is `##<key>="…"`
/// whose value carries the engine attribute, and a command line is
/// `##<key>-command-line='…'` whose JSON carries the engine field. An unrelated
/// metadata line that merely mentions the engine — say
/// `##audit={"engine":"…","status":"reviewed"}` — is valid source metadata and
/// survives, because its value opens with `{` rather than a quote.
fn is_stale_provenance_line(line: &str, tool: &str) -> bool {
    if line.starts_with("##VEP=") || line.starts_with("##VEP-command-line=") {
        return true;
    }
    if line.starts_with(&format!("##{tool}="))
        || line.starts_with(&format!("##{tool}-command-line="))
    {
        return true;
    }

    let Some((key, value)) = line
        .strip_prefix("##")
        .and_then(|rest| rest.split_once('='))
    else {
        return false;
    };
    if key.is_empty() || key.contains(char::is_whitespace) {
        return false;
    }

    match key.strip_suffix("-command-line") {
        Some(name) => {
            !name.is_empty() && value.starts_with('\'') && value.contains(ENGINE_MARKER_JSON)
        }
        None => value.starts_with('"') && value.contains(ENGINE_MARKER_ATTR),
    }
}

fn provenance_header_lines(
    input_vcf: &str,
    cache_source: &str,
    backend: &str,
    output_vcf: &str,
    config: &AnnotateVcfConfig,
    cache_source_type: CacheSourceType,
) -> Vec<String> {
    let tool = config
        .provenance_tool_name
        .as_deref()
        .unwrap_or(env!("CARGO_PKG_NAME"));
    let version = config
        .provenance_tool_version
        .as_deref()
        .unwrap_or(env!("CARGO_PKG_VERSION"));

    let mut identity = format!(
        "##{tool}=\"{version}\" engine=\"{} {}\" cache=\"{cache_source}\" \
         cache_type=\"{}\" backend=\"{backend}\"",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION"),
        cache_source_type.as_str(),
    );

    // The pinned VEP code/cache pair this run reproduces. Only stated when the
    // cache version is known up front; an unpinned run says nothing rather than
    // guessing.
    if let Some(cache_version) = config.expected_cache_version.as_deref()
        && let Ok(target) = crate::vep_semantics::target_for_cache_version(cache_version)
    {
        identity.push_str(&format!(
            " cache_version=\"{}\" vep_codebase=\"{}\" api=\"{}\" ensembl={}.{} \
             ensembl-variation={}.{}",
            target.cache_version,
            target.vep_codebase_version,
            target.api_version,
            target.cache_version,
            target.ensembl_core_revision,
            target.cache_version,
            target.ensembl_variation_revision,
        ));
    }

    let invocation = serde_json::json!({
        "engine": env!("CARGO_PKG_NAME"),
        "input": input_vcf,
        "output": output_vcf,
        "cache": cache_source,
        "backend": backend,
        "reference_fasta": config.reference_fasta_path,
        "plugin_cache_root": config.plugin_cache_root
            .as_ref()
            .map(|p| p.display().to_string()),
        "options": serde_json::from_str::<serde_json::Value>(
            &config.to_options_json_with_backend(backend),
        )
        .unwrap_or(serde_json::Value::Null),
    });

    vec![identity, format!("##{tool}-command-line='{invocation}'")]
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
    let mut format_list = field_names.join("|");
    // Trailing custom-plugin CSQ fields (spec §5): read cheaply from manifests.
    #[cfg(feature = "parquet-cache")]
    if let Some(root) = &config.plugin_cache_root {
        for name in crate::plugin_cache::registry::PluginRegistry::field_names(root) {
            format_list.push('|');
            format_list.push_str(&name);
        }
    }
    // Byte-identical to the line Ensembl VEP writes in OutputFactory/VCF.pm.
    // Downstream tooling matches on this prefix, so it must not be reworded.
    format!("Consequence annotations from Ensembl VEP. Format: {format_list}")
}

/// Canonical 28-byte BGZF end-of-file marker (an empty BGZF block, SAM spec).
/// Every finished bgzf file (`VcfLocalWriter::Bgzf::finish()`) ends with exactly
/// these bytes. Written once, at the very end of an assembled BGZF file; stripped
/// from each shard/header piece before concatenation.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Read a finished bgzf file and return its bytes with the trailing 28-byte EOF
/// marker removed, ready for raw concatenation into a larger bgzf stream.
/// Read a small finished bgzf file (the header) whole and return its bytes with
/// the trailing 28-byte EOF marker removed. Used only for the header, which is
/// tiny; shards stream via [`append_bgzf_shard`] instead. The EOF check is a real
/// runtime error (not a `debug_assert`) so it holds in release builds.
fn bgzf_blocks_no_eof(path: &Path) -> Result<Vec<u8>> {
    let mut bytes = std::fs::read(path).map_err(|e| {
        DataFusionError::Execution(format!("read bgzf file {}: {e}", path.display()))
    })?;
    if bytes.len() < BGZF_EOF.len() || bytes[bytes.len() - BGZF_EOF.len()..] != BGZF_EOF {
        return Err(DataFusionError::Execution(format!(
            "bgzf file {} does not end with the canonical EOF marker",
            path.display()
        )));
    }
    let keep = bytes.len() - BGZF_EOF.len();
    bytes.truncate(keep);
    Ok(bytes)
}

/// Append a finished bgzf shard file to `out`, stripping its trailing 28-byte EOF
/// marker, streaming the body in `buf`-sized chunks so a large shard is never
/// loaded whole into memory (mirrors the plain assembler's bounded-buffer copy).
/// Uses the file length to locate the marker, then reads + validates it. Returns
/// a real runtime error (holds in release) if the shard is shorter than the
/// marker or does not end with it — never silently truncating corrupt input.
fn append_bgzf_shard<W: Write>(path: &Path, out: &mut W, buf: &mut [u8]) -> Result<()> {
    let mut f = std::fs::File::open(path).map_err(|e| {
        DataFusionError::Execution(format!("assembler open bgzf shard {}: {e}", path.display()))
    })?;
    let total = f
        .metadata()
        .map_err(|e| {
            DataFusionError::Execution(format!("assembler stat bgzf shard {}: {e}", path.display()))
        })?
        .len();
    let eof_len = BGZF_EOF.len() as u64;
    if total < eof_len {
        return Err(DataFusionError::Execution(format!(
            "bgzf shard {} is {total} bytes, shorter than the canonical EOF marker",
            path.display()
        )));
    }
    // Stream the body (everything except the trailing EOF marker) in bounded chunks.
    let mut remaining = total - eof_len;
    while remaining > 0 {
        let want = remaining.min(buf.len() as u64) as usize;
        f.read_exact(&mut buf[..want]).map_err(|e| {
            DataFusionError::Execution(format!("assembler read bgzf shard {}: {e}", path.display()))
        })?;
        out.write_all(&buf[..want])
            .map_err(|e| DataFusionError::Execution(format!("assembler write bgzf shard: {e}")))?;
        remaining -= want as u64;
    }
    // Read + validate the trailing EOF marker; never write it (the assembler emits
    // exactly one canonical EOF at the very end of the assembled file).
    let mut marker = [0u8; 28]; // == BGZF_EOF.len()
    f.read_exact(&mut marker).map_err(|e| {
        DataFusionError::Execution(format!(
            "assembler read bgzf EOF of {}: {e}",
            path.display()
        ))
    })?;
    if marker != BGZF_EOF {
        return Err(DataFusionError::Execution(format!(
            "bgzf shard {} does not end with the canonical EOF marker",
            path.display()
        )));
    }
    Ok(())
}

/// Dedicated assembler thread for Plain/Gzip output. Owns the header-primed
/// `VcfLocalWriter` (targeting `temp_output`) and appends each completed contig's
/// shard range — in arrival (= ascending id = position) order — as raw body
/// bytes (Plain = byte copy; Gzip = compress-on-append). Runs concurrently with
/// the next contig's prepare + annotation, so the copy overlaps rather than
/// tailing after 100%. On channel close (annotation finished, no error) it
/// finalizes the writer and atomically renames `temp_output` → `final_output`.
/// Returns the total body-line count (summed from worker-reported `output_lines`,
/// so no re-scan). Never advances the progress bar (already driven by rows_done).
fn run_assembler_thread(
    mut writer: VcfLocalWriter,
    final_output: PathBuf,
    temp_output: PathBuf,
    tempdir: PathBuf,
    rx: std::sync::mpsc::Receiver<ContigShardRange>,
    finalize: Arc<std::sync::atomic::AtomicBool>,
) -> std::thread::JoinHandle<Result<usize>> {
    std::thread::spawn(move || -> Result<usize> {
        let mut total_lines = 0usize;
        let mut buffer = vec![0u8; 8 * 1024 * 1024];
        while let Ok(range) = rx.recv() {
            for id in range.first_shard_id..range.end_shard_id {
                let path = tempdir.join(format!("partition_{id:04}.vcf.body"));
                let mut f = std::fs::File::open(&path).map_err(|e| {
                    DataFusionError::Execution(format!(
                        "assembler open shard {}: {e}",
                        path.display()
                    ))
                })?;
                loop {
                    let n = f.read(&mut buffer).map_err(|e| {
                        DataFusionError::Execution(format!(
                            "assembler read shard {}: {e}",
                            path.display()
                        ))
                    })?;
                    if n == 0 {
                        break;
                    }
                    write_vcf_body_chunk(&mut writer, &buffer[..n])?;
                }
                let _ = std::fs::remove_file(&path);
            }
            total_lines += range.output_lines;
        }
        // Channel closed. Finalize (finish + atomic rename) ONLY if the driver
        // signalled success; otherwise the run failed mid-annotation and this temp
        // holds only the contigs completed before the error — discard it so no
        // truncated-but-valid-looking VCF is ever renamed into the output path.
        if !finalize.load(std::sync::atomic::Ordering::Acquire) {
            drop(writer);
            let _ = std::fs::remove_file(&temp_output);
            return Ok(0);
        }
        writer.finish()?;
        std::fs::rename(&temp_output, &final_output).map_err(|e| {
            DataFusionError::Execution(format!(
                "assembler rename {} -> {}: {e}",
                temp_output.display(),
                final_output.display()
            ))
        })?;
        Ok(total_lines)
    })
}

/// Dedicated assembler thread for BGZF output. Shards are already bgzf-compressed
/// (by the workers, in parallel), so this owns a RAW `BufWriter<File>` — it must
/// NOT re-compress. It writes the pre-compressed `header_blocks`, then for each
/// completed contig range raw-appends every shard minus its trailing 28-byte EOF
/// marker, then writes one canonical `BGZF_EOF`, flushes, and atomically renames.
/// Returns the total body-line count (worker-reported). Concurrent with the next
/// contig's prepare + annotation, exactly like the Plain assembler.
fn run_assembler_thread_bgzf(
    temp_output: PathBuf,
    final_output: PathBuf,
    tempdir: PathBuf,
    header_blocks: Vec<u8>,
    rx: std::sync::mpsc::Receiver<ContigShardRange>,
    finalize: Arc<std::sync::atomic::AtomicBool>,
) -> std::thread::JoinHandle<Result<usize>> {
    std::thread::spawn(move || -> Result<usize> {
        let file = std::fs::File::create(&temp_output).map_err(|e| {
            DataFusionError::Execution(format!("assembler create {}: {e}", temp_output.display()))
        })?;
        let mut out = std::io::BufWriter::new(file);
        out.write_all(&header_blocks)
            .map_err(|e| DataFusionError::Execution(format!("assembler write header: {e}")))?;
        let mut total_lines = 0usize;
        let mut buffer = vec![0u8; 8 * 1024 * 1024];
        while let Ok(range) = rx.recv() {
            for id in range.first_shard_id..range.end_shard_id {
                let path = tempdir.join(format!("partition_{id:04}.vcf.body"));
                append_bgzf_shard(&path, &mut out, &mut buffer)?;
                let _ = std::fs::remove_file(&path);
            }
            total_lines += range.output_lines;
        }
        // Channel closed. Discard the partial temp on a mid-run failure (see the
        // plain assembler) — do NOT write the EOF or rename a truncated file.
        if !finalize.load(std::sync::atomic::Ordering::Acquire) {
            drop(out);
            let _ = std::fs::remove_file(&temp_output);
            return Ok(0);
        }
        out.write_all(&BGZF_EOF)
            .map_err(|e| DataFusionError::Execution(format!("assembler write EOF: {e}")))?;
        out.flush()
            .map_err(|e| DataFusionError::Execution(format!("assembler flush: {e}")))?;
        drop(out);
        std::fs::rename(&temp_output, &final_output).map_err(|e| {
            DataFusionError::Execution(format!(
                "assembler rename {} -> {}: {e}",
                temp_output.display(),
                final_output.display()
            ))
        })?;
        Ok(total_lines)
    })
}

/// Drive the `threads>1` sharded-VCF-output path: build the annotation plan
/// directly (bypassing SQL) with a [`VcfShardContext`] so each fused worker
/// streams its own position-ordered VCF body shard into `shard_ctx.tempdir`,
/// then concatenate the shards in ascending id (= position) order into `writer`.
/// Returns the total number of body lines written.
#[allow(clippy::too_many_arguments)]
async fn drive_sharded_vcf_annotation(
    ctx: &SessionContext,
    vcf_table: &str,
    cache_source: &str,
    backend: &str,
    cache_source_type: CacheSourceType,
    options_json: String,
    vcf_schema: &Arc<datafusion::arrow::datatypes::Schema>,
    projection_names: &[String],
    shard_ctx: Arc<VcfShardContext>,
    pb: &ProgressBar,
    config: &AnnotateVcfConfig,
    total_input: usize,
) -> Result<()> {
    use crate::annotate_provider::AnnotateProvider;
    use crate::annotation_store::AnnotationBackend;
    use datafusion::physical_plan::ExecutionPlan;
    use futures::StreamExt;

    let backend_enum = AnnotationBackend::parse(backend)?;
    // Reuse AnnotateProvider::new + TableProvider::scan (crate-reachable, all
    // inputs derived from vcf_schema, not SQL) and inject the shard context.
    let provider = AnnotateProvider::new(
        Arc::new(ctx.clone()),
        vcf_table.to_string(),
        cache_source.to_string(),
        backend_enum,
        cache_source_type,
        Some(options_json),
        (**vcf_schema).clone(),
    )?
    .with_vcf_shard_ctx(Arc::clone(&shard_ctx))
    .with_plugin_cache_root(config.plugin_cache_root.clone());

    let full_schema = TableProvider::schema(&provider);
    let projection: Vec<usize> = projection_names
        .iter()
        .map(|name| full_schema.index_of(name).map_err(DataFusionError::from))
        .collect::<Result<_>>()?;

    let state = ctx.state();
    let plan = provider.scan(&state, Some(&projection), &[], None).await?;
    let partition_count = plan.properties().output_partitioning().partition_count();
    if partition_count != 1 {
        return Err(DataFusionError::Internal(format!(
            "sharded VCF output expects a single output partition, got {partition_count}; \
             threads>1 must not be combined with forks/contig_parallelism"
        )));
    }

    // Drive the plan to completion: the workers write their shards; this stream
    // yields no row batches. `?` propagates the first worker error (no concat).
    // Workers bump `shard_ctx.rows_done` as they annotate each window; a sibling
    // OS thread polls it so the progress bar / callback advance DURING annotation
    // rather than jumping only at concat. The poller uses `std::thread::sleep`
    // (not a tokio timer) so it works regardless of whether the caller's runtime
    // enabled the time driver.
    use std::sync::atomic::Ordering::Relaxed;
    let mut stream = plan.execute(0, ctx.task_ctx())?;
    let progress_started = Instant::now();
    let progress_trace = std::env::var_os("VEP_PROFILE").is_some();
    let mut last = 0usize;
    let mut report = |now: usize| {
        if now > last {
            let prev = last;
            last = now;
            pb.set_position(now as u64);
            if progress_trace {
                eprintln!(
                    "[VEP_PROGRESS] annotated {now}/{total_input} (+{}) at {:.1}s",
                    now - prev,
                    progress_started.elapsed().as_secs_f64(),
                );
            }
            if let Some(ref cb) = config.on_batch_written {
                cb(now - prev, now, total_input);
            }
        }
    };

    // The sharded stream yields no row batches — workers write shards directly
    // and bump `shard_ctx.rows_done` as they annotate. A dedicated OS thread
    // samples that counter every 150ms and forwards snapshots over a channel, so
    // the progress bar / callback advance DURING annotation. We AWAIT the stream
    // on the caller's runtime (instead of `block_in_place` + `block_on`, which
    // panics on a current-thread runtime), keeping the worker tasks driven on
    // whatever runtime spawned them. The poller uses `std::thread::sleep` (not a
    // tokio timer) and an unbounded `mpsc`, so it needs no time driver.
    let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let (tick_tx, mut tick_rx) = tokio::sync::mpsc::unbounded_channel::<usize>();
    let poller = {
        let stop = Arc::clone(&stop);
        let shard_ctx = Arc::clone(&shard_ctx);
        std::thread::spawn(move || {
            while !stop.load(Relaxed) {
                std::thread::sleep(std::time::Duration::from_millis(150));
                let _ = tick_tx.send(shard_ctx.rows_done.load(Relaxed));
            }
        })
    };

    let drive_result: Result<()> = loop {
        tokio::select! {
            item = stream.next() => match item {
                Some(batch) => {
                    if let Err(e) = batch {
                        break Err(e);
                    }
                }
                None => break Ok(()),
            },
            Some(now) = tick_rx.recv() => report(now),
        }
    };
    stop.store(true, Relaxed);
    let _ = poller.join();
    drive_result?;
    // Final flush: report any rows annotated since the poller's last tick.
    report(shard_ctx.rows_done.load(Relaxed));

    // Shard assembly no longer happens here: each contig's shard range is
    // published (annotate_provider) and concatenated by the dedicated assembler
    // thread spawned in `annotate_to_vcf`, concurrently with the next contig's
    // prepare + annotation. This function only drives annotation + progress.
    Ok(())
}

pub async fn annotate_to_vcf(
    input_vcf: &str,
    cache_source: &str,
    backend: &str,
    output_vcf: &str,
    config: &AnnotateVcfConfig,
) -> Result<usize> {
    if config.workers == 0 {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): workers must be a positive integer".to_string(),
        ));
    }
    if config.target_partitions == 0 {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): target_partitions must be a positive integer".to_string(),
        ));
    }
    if config.refseq || config.merged {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): refseq and merged config fields are unsupported; cache source mode must come from cache schema metadata bio.vep.cache_source_type".to_string(),
        ));
    }
    // Parallel annotation (workers>1) splits the contig into N position-range
    // lookup partitions. That requires a bgzipped + tabix-indexed input so the
    // VCF scan yields contiguous, position-ordered partitions; an unindexed
    // input would force a round-robin repartition that scrambles output order.
    if config.workers > 1
        && !std::path::Path::new(&format!("{input_vcf}.tbi")).exists()
        && !std::path::Path::new(&format!("{input_vcf}.csi")).exists()
    {
        return Err(DataFusionError::Plan(format!(
            "annotate_to_vcf(): workers>1 requires a tabix-indexed input (`{input_vcf}.tbi` or `.csi`); \
             bgzip + tabix the VCF, or run with workers=1"
        )));
    }
    let cache_source_type = cache_source_type_from_cache_source_for_backend(cache_source, backend)?;
    let concurrency_plan = VepConcurrencyPlan::from_config(config);
    if sink_profile_enabled() {
        eprintln!(
            "[VEP_PROFILE] concurrency_plan lookup_partitions={} workers={} cold_parquet_target_partitions={} spawn_vcf_provider_open={}",
            concurrency_plan.lookup_partitions,
            config.workers,
            config.target_partitions,
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

    let options_json = config.to_options_json_with_backend(backend);
    let opts_clause = format!(", '{}'", options_json.replace('\'', "''"));
    let sql = format!(
        "SELECT {select_list} FROM annotate_vep('{vcf_table}', '{}', '{}'{opts_clause})",
        cache_source.replace('\'', "''"),
        backend.replace('\'', "''"),
    );

    // Column set the formatter consumes (same as the serial SELECT list): core
    // VCF columns + INFO fields + CSQ + per-sample FORMAT columns. Used by the
    // threads>1 sharded path to project the annotation plan identically.
    let projection_names: Vec<String> = core_vcf
        .iter()
        .map(|name| (*name).to_string())
        .chain(info_fields.iter().cloned())
        .chain(std::iter::once("CSQ".to_string()))
        .chain(format_fields.iter().cloned())
        .collect();

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

    // Append this run's provenance to the captured source header, so the writer
    // emits it verbatim after the input's own lines and before the CSQ
    // declaration. Any provenance from an earlier annotation is stripped first:
    // it describes a CSQ field that is being replaced. Ensembl VEP does the same
    // with its own `##VEP` lines (OutputFactory/VCF.pm).
    let mut write_metadata = vcf_schema.metadata().clone();
    if let Some(raw_json) = write_metadata.get(VCF_HEADER_RAW_LINES_KEY)
        && let Ok(existing) = serde_json::from_str::<Vec<String>>(raw_json)
    {
        let tool = config
            .provenance_tool_name
            .as_deref()
            .unwrap_or(env!("CARGO_PKG_NAME"));
        let mut lines: Vec<String> = existing
            .into_iter()
            .filter(|line| !is_stale_provenance_line(line, tool))
            .collect();
        lines.extend(provenance_header_lines(
            input_vcf,
            cache_source,
            backend,
            output_vcf,
            config,
            cache_source_type,
        ));
        if let Ok(json) = serde_json::to_string(&lines) {
            write_metadata.insert(VCF_HEADER_RAW_LINES_KEY.to_string(), json);
        }
    }

    let write_schema = Arc::new(
        datafusion::arrow::datatypes::Schema::new(output_fields).with_metadata(write_metadata),
    );

    // 6. Stream annotated batches to VCF file.
    let output_path = Path::new(output_vcf);

    let coordinate_zero_based = vcf_schema
        .metadata()
        .get("bio.coordinate_system_zero_based")
        .is_some_and(|v| v == "true");

    let mut total_rows = 0;
    let mut sink_profile = sink_profile_enabled().then(VcfSinkProfile::default);

    if config.workers > 1 {
        // Sharded VCF output: bypass the SQL/DataFrame execution and drive the
        // annotation plan directly with a VcfShardContext, so each fused worker
        // streams its own position-ordered VCF body shard. A dedicated assembler
        // thread concatenates each contig's shards the instant that contig
        // finishes — concurrently with the next contig's prepare + annotation —
        // then atomically renames the temp output into place.
        let tempdir = VcfBodyTempDir::new()?;
        // Temp assembly path in the SAME directory as the output (so the final
        // rename is atomic on one filesystem).
        let mut temp_os = output_path.as_os_str().to_owned();
        temp_os.push(".part");
        let temp_output = PathBuf::from(temp_os);

        let (contig_done_tx, contig_done_rx) = std::sync::mpsc::channel::<ContigShardRange>();
        // bgzf output → parallel per-worker bgzf shards + raw block concat;
        // plain/gzip → plain text shards + re-copy/compress at the assembler.
        let shard_compression = if config.compression == VcfCompressionType::Bgzf {
            VcfCompressionType::Bgzf
        } else {
            VcfCompressionType::Plain
        };

        // Spawn the assembler (owns the output writer); its mode matches the
        // output compression. Header is written BEFORE the field vecs are moved
        // into the shard context.
        // Set to true only if annotation completes successfully; the assembler
        // reads it after the channel closes and finalizes (finish + atomic rename)
        // ONLY when true. On a mid-run error the assembler discards its partial
        // temp file instead of renaming a truncated-but-valid-looking VCF into the
        // output path.
        let finalize = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let assembler = if config.compression == VcfCompressionType::Bgzf {
            // Pre-compress the header into bgzf blocks (strip its EOF); the raw
            // assembler writes those, then the already-compressed shards.
            let header_tmp = tempdir.path().join("header.bgz");
            let mut hw = VcfLocalWriter::with_compression(&header_tmp, VcfCompressionType::Bgzf)?;
            hw.write_header(
                &write_schema,
                &vcf_info_fields,
                &unique_format_tags,
                &sample_names,
            )?;
            hw.finish()?;
            let header_blocks = bgzf_blocks_no_eof(&header_tmp)?;
            let _ = std::fs::remove_file(&header_tmp);
            run_assembler_thread_bgzf(
                temp_output.clone(),
                output_path.to_path_buf(),
                tempdir.path().to_path_buf(),
                header_blocks,
                contig_done_rx,
                Arc::clone(&finalize),
            )
        } else {
            // Plain/Gzip: write the header to the temp writer; the assembler
            // owns it, appends shard bodies, finishes, and renames.
            let mut writer = VcfLocalWriter::with_compression(&temp_output, config.compression)?;
            writer.write_header(
                &write_schema,
                &vcf_info_fields,
                &unique_format_tags,
                &sample_names,
            )?;
            run_assembler_thread(
                writer,
                output_path.to_path_buf(),
                temp_output.clone(),
                tempdir.path().to_path_buf(),
                contig_done_rx,
                Arc::clone(&finalize),
            )
        };

        let shard_ctx = Arc::new(VcfShardContext {
            vcf_info_fields: Arc::new(vcf_info_fields),
            unique_format_tags: Arc::new(unique_format_tags),
            sample_names: Arc::new(sample_names),
            coordinate_zero_based,
            shard_compression,
            tempdir: tempdir.path().to_path_buf(),
            rows_done: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            contig_done_tx,
        });
        let drive_result = drive_sharded_vcf_annotation(
            &ctx,
            vcf_table,
            cache_source,
            backend,
            cache_source_type,
            options_json.clone(),
            &vcf_schema,
            &projection_names,
            Arc::clone(&shard_ctx),
            &pb,
            config,
            total_input,
        )
        .await;
        // Signal the assembler whether to finalize BEFORE closing the channel, so
        // a mid-run error never renames a partial file into the output path.
        if drive_result.is_ok() {
            finalize.store(true, std::sync::atomic::Ordering::Release);
        }
        // Drop the last sender (the one held on shard_ctx) so the assembler sees
        // the channel close. Any per-contig senders were already sent-and-dropped
        // by the stream.
        drop(shard_ctx);
        // ALWAYS join the assembler (even on error) so the thread is never detached
        // and its cleanup/result is observed. Propagate the drive error first (the
        // assembler has already discarded its partial temp on failure), then any
        // assembler error.
        let assembled = assembler
            .join()
            .map_err(|_| DataFusionError::Execution("assembler thread panicked".to_string()))?;
        drive_result?;
        total_rows = assembled?;
        // tempdir (and any leftover shard files) removed on drop here.
        drop(tempdir);
    } else {
        let mut writer = VcfLocalWriter::with_compression(output_path, config.compression)?;
        writer.write_header(
            &write_schema,
            &vcf_info_fields,
            &unique_format_tags,
            &sample_names,
        )?;
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
        let finish_started = Instant::now();
        writer.finish()?;
        if let Some(profile) = sink_profile.as_mut() {
            profile.writer_finish += finish_started.elapsed();
        }
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

    // bgzf is a stream of concatenated gzip members; flate2's MultiGzDecoder
    // decodes it without needing an external bgzip binary.
    fn gunzip_to_string(path: &Path) -> String {
        use std::io::Read;
        let f = std::fs::File::open(path).unwrap();
        let mut d = flate2::read::MultiGzDecoder::new(f);
        let mut s = String::new();
        d.read_to_string(&mut s).unwrap();
        s
    }

    fn unique_dir(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn contig_shard_range_is_contiguous_and_sums() {
        let a = ContigShardRange {
            first_shard_id: 0,
            end_shard_id: 3,
            output_lines: 10,
        };
        let b = ContigShardRange {
            first_shard_id: 3,
            end_shard_id: 5,
            output_lines: 4,
        };
        assert_eq!(a.end_shard_id, b.first_shard_id); // contiguous, no gaps/overlap
        assert_eq!(a.output_lines + b.output_lines, 14);
    }

    #[test]
    fn assembler_plain_appends_ranges_in_order_and_counts_lines() {
        use std::io::Write;
        let dir = unique_dir("asm_plain");
        for (id, body) in [(0usize, "a\nb\n"), (1, "c\n"), (2, "d\ne\nf\n")] {
            let mut f =
                std::fs::File::create(dir.join(format!("partition_{id:04}.vcf.body"))).unwrap();
            f.write_all(body.as_bytes()).unwrap();
        }
        let out = dir.join("out.vcf");
        let tmp = dir.join("out.vcf.part");
        let writer = VcfLocalWriter::with_compression(&tmp, VcfCompressionType::Plain).unwrap();
        let (tx, rx) = std::sync::mpsc::channel();
        let h = run_assembler_thread(
            writer,
            out.clone(),
            tmp,
            dir.clone(),
            rx,
            Arc::new(std::sync::atomic::AtomicBool::new(true)),
        );
        tx.send(ContigShardRange {
            first_shard_id: 0,
            end_shard_id: 2,
            output_lines: 3,
        })
        .unwrap();
        tx.send(ContigShardRange {
            first_shard_id: 2,
            end_shard_id: 3,
            output_lines: 3,
        })
        .unwrap();
        drop(tx);
        let total = h.join().unwrap().unwrap();
        assert_eq!(total, 6);
        assert_eq!(std::fs::read_to_string(&out).unwrap(), "a\nb\nc\nd\ne\nf\n");
        assert!(
            !dir.join("partition_0000.vcf.body").exists(),
            "shards removed after copy"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn shard_writer_bgzf_empty_shard_is_decodable() {
        let dir = unique_dir("shard_bgzf");
        let path = dir.join("partition_0000.vcf.body");
        let w = VcfBodyShardWriter::create(
            &path,
            Arc::new(vec![]),
            Arc::new(vec![]),
            Arc::new(vec![]),
            false,
            VcfCompressionType::Bgzf,
        )
        .unwrap();
        let lines = w.lines;
        w.finish().unwrap(); // flushes bgzf blocks + EOF
        assert_eq!(lines, 0);
        assert_eq!(gunzip_to_string(&path), "");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn assembler_bgzf_raw_concat_decodes_in_order() {
        let dir = unique_dir("asm_bgzf");
        // Produce bgzf shards exactly as VcfBodyShardWriter(Bgzf) does: body
        // through a VcfLocalWriter::Bgzf, finish() (appends EOF).
        let write_bgzf = |path: &Path, body: &str| {
            let mut w = VcfLocalWriter::with_compression(path, VcfCompressionType::Bgzf).unwrap();
            write_vcf_body_chunk(&mut w, body.as_bytes()).unwrap();
            w.finish().unwrap();
        };
        for (id, body) in [(0usize, "a\nb\n"), (1usize, "c\n")] {
            write_bgzf(&dir.join(format!("partition_{id:04}.vcf.body")), body);
        }
        // Header blocks: compress "##header\n#CHROM\n" and strip its EOF.
        let hpath = dir.join("hdr.bgz");
        write_bgzf(&hpath, "##header\n#CHROM\n");
        let header_blocks = bgzf_blocks_no_eof(&hpath).unwrap();
        std::fs::remove_file(&hpath).ok();

        let out = dir.join("out.vcf.gz");
        let tmp = dir.join("out.vcf.gz.part");
        let (tx, rx) = std::sync::mpsc::channel();
        let h = run_assembler_thread_bgzf(
            tmp.clone(),
            out.clone(),
            dir.clone(),
            header_blocks,
            rx,
            Arc::new(std::sync::atomic::AtomicBool::new(true)),
        );
        tx.send(ContigShardRange {
            first_shard_id: 0,
            end_shard_id: 1,
            output_lines: 2,
        })
        .unwrap();
        tx.send(ContigShardRange {
            first_shard_id: 1,
            end_shard_id: 2,
            output_lines: 1,
        })
        .unwrap();
        drop(tx);
        let total = h.join().unwrap().unwrap();
        assert_eq!(total, 3);
        assert_eq!(gunzip_to_string(&out), "##header\n#CHROM\na\nb\nc\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn assembler_discards_partial_output_when_not_finalized() {
        // On a mid-run failure the driver never sets `finalize`; the assembler must
        // NOT finish/rename its partial temp into the final output path — it must
        // discard the temp and produce no output file.
        use std::io::Write;
        let dir = unique_dir("asm_nofinal");
        let mut f = std::fs::File::create(dir.join("partition_0000.vcf.body")).unwrap();
        f.write_all(b"a\nb\n").unwrap();
        let out = dir.join("out.vcf");
        let tmp = dir.join("out.vcf.part");
        let writer = VcfLocalWriter::with_compression(&tmp, VcfCompressionType::Plain).unwrap();
        let (tx, rx) = std::sync::mpsc::channel();
        let h = run_assembler_thread(
            writer,
            out.clone(),
            tmp.clone(),
            dir.clone(),
            rx,
            Arc::new(std::sync::atomic::AtomicBool::new(false)), // driver "failed"
        );
        tx.send(ContigShardRange {
            first_shard_id: 0,
            end_shard_id: 1,
            output_lines: 2,
        })
        .unwrap();
        drop(tx); // channel closes, but finalize == false
        let total = h.join().unwrap().unwrap();
        assert_eq!(total, 0, "no rows finalized on failure");
        assert!(!out.exists(), "final output must NOT be created on failure");
        assert!(!tmp.exists(), "temp .part must be cleaned up on failure");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn assembler_bgzf_errors_on_shard_missing_eof_marker() {
        // A shard that does not end with the canonical 28-byte BGZF EOF marker must
        // produce a runtime error (not silent truncation) — the check must hold in
        // release builds too, so it cannot be a debug_assert.
        use std::io::Write;
        let dir = unique_dir("asm_bgzf_badeof");
        let mut f = std::fs::File::create(dir.join("partition_0000.vcf.body")).unwrap();
        f.write_all(b"this shard does not end with the bgzf eof marker")
            .unwrap();
        let out = dir.join("out.vcf.gz");
        let tmp = dir.join("out.vcf.gz.part");
        let (tx, rx) = std::sync::mpsc::channel();
        let h = run_assembler_thread_bgzf(
            tmp,
            out.clone(),
            dir.clone(),
            Vec::new(),
            rx,
            Arc::new(std::sync::atomic::AtomicBool::new(true)),
        );
        tx.send(ContigShardRange {
            first_shard_id: 0,
            end_shard_id: 1,
            output_lines: 1,
        })
        .unwrap();
        drop(tx);
        let res = h.join().unwrap();
        assert!(
            res.is_err(),
            "must error on a shard without the canonical EOF marker"
        );
        assert!(!out.exists(), "no output file on error");
        std::fs::remove_dir_all(&dir).ok();
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

    // Guards the I1 wiring bug: the workers=1 SQL path builds the AnnotateProvider
    // from options_json, so plugin_cache_root MUST round-trip through it (else the
    // plugin registry is never opened and CSQ plugin fields emit empty while the
    // header still lists them — a header/body width divergence).
    #[test]
    fn test_to_options_json_emits_plugin_cache_root() {
        let none = AnnotateVcfConfig::default().to_options_json();
        assert!(!none.contains("plugin_cache_root"));
        let config = AnnotateVcfConfig {
            plugin_cache_root: Some(std::path::PathBuf::from("/tmp/plugin_cache")),
            ..Default::default()
        };
        let json = config.to_options_json();
        assert!(json.contains("\"plugin_cache_root\":\"/tmp/plugin_cache\""));
    }

    #[test]
    fn test_to_options_json_with_backend_emits_cache_format() {
        let config = AnnotateVcfConfig::default();

        let json = config.to_options_json_with_backend("lance");
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(value["cache_format"], "parquet");
    }

    #[test]
    fn concurrency_plan_serial_uses_single_lookup_partition() {
        let config = AnnotateVcfConfig::default();
        let plan = VepConcurrencyPlan::from_config(&config);
        assert_eq!(plan.lookup_partitions, 1);
        assert!(!plan.spawn_vcf_provider_open);
    }

    #[test]
    fn concurrency_plan_workers_drive_lookup_partitions() {
        let config = AnnotateVcfConfig {
            workers: 8,
            ..AnnotateVcfConfig::default()
        };
        let plan = VepConcurrencyPlan::from_config(&config);
        assert_eq!(plan.lookup_partitions, 8);
        assert!(plan.spawn_vcf_provider_open);
    }

    #[test]
    fn to_options_json_emits_workers_not_threads_or_forks() {
        let config = AnnotateVcfConfig {
            workers: 4,
            ..AnnotateVcfConfig::default()
        };
        let json: serde_json::Value = serde_json::from_str(&config.to_options_json()).unwrap();
        assert_eq!(json["workers"], 4);
        assert!(json.get("threads").is_none());
        assert!(json.get("forks").is_none());
        assert!(json.get("contig_parallelism").is_none());
        assert!(json.get("inline_lookup").is_none());
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

    /// Ensembl VEP writes this header line from `OutputFactory/VCF.pm`:
    ///
    ///   '##INFO=<ID=%s,Number=.,Type=String,Description="Consequence
    ///    annotations from Ensembl VEP. Format: %s">'
    ///
    /// Byte parity with a VEP run requires the prefix to match exactly, so it is
    /// pinned here rather than described loosely.
    #[test]
    fn stale_filter_removes_ensembl_vep_provenance() {
        assert!(is_stale_provenance_line(
            "##VEP=\"v116.0\" API=\"v116\"",
            "vepyr"
        ));
        assert!(is_stale_provenance_line(
            "##VEP-command-line='vep --everything'",
            "vepyr"
        ));
    }

    #[test]
    fn stale_filter_keeps_unrelated_keys_that_merely_start_with_vep() {
        // Arbitrary meta-information keys are valid VCF. `##VEPStatus` is not
        // Ensembl provenance and the source header must keep it.
        for line in [
            "##VEPStatus=reviewed",
            "##VEPTools=custom",
            "##VEP_NOTES=see methods",
        ] {
            assert!(!is_stale_provenance_line(line, "vepyr"), "{line}");
        }
    }

    #[test]
    fn stale_filter_removes_provenance_written_under_another_tool_name() {
        // A file annotated by an engine-default run and then re-annotated by a
        // wrapper that sets `provenance_tool_name` must not keep the first run's
        // lines: they describe the CSQ being replaced. Recognition keys off this
        // engine's identity, which is present whatever the configured name.
        let engine_default = provenance_header_lines(
            "/in.vcf",
            "/cache",
            "parquet",
            "/out.vcf",
            &AnnotateVcfConfig::default(),
            CacheSourceType::Ensembl,
        );
        for line in &engine_default {
            assert!(
                is_stale_provenance_line(line, "vepyr"),
                "not recognised under a new tool name: {line}"
            );
        }
    }

    #[test]
    fn stale_filter_removes_its_own_provenance_under_the_same_name() {
        let config = AnnotateVcfConfig {
            provenance_tool_name: Some("vepyr".to_string()),
            ..Default::default()
        };
        for line in provenance_for(&config) {
            assert!(is_stale_provenance_line(&line, "vepyr"), "{line}");
        }
    }

    #[test]
    fn stale_filter_keeps_generic_metadata_that_merely_mentions_the_engine() {
        // Recognition must key off the *shape* this module emits, not a bare
        // substring. An unrelated metadata line that happens to name the engine
        // is valid source metadata and must survive.
        for line in [
            r#"##audit={"engine":"datafusion-bio-function-vep","status":"reviewed"}"#,
            r#"##pipeline={"steps":[{"engine":"datafusion-bio-function-vep"}]}"#,
            "##INFO=<ID=X,Number=1,Type=String,Description=\"engine=\\\"datafusion-bio-function-vep 1.0\\\"\">",
        ] {
            assert!(!is_stale_provenance_line(line, "vepyr"), "{line}");
        }
    }

    #[test]
    fn stale_filter_keeps_ordinary_header_lines() {
        for line in [
            "##fileformat=VCFv4.2",
            "##fileDate=20160824",
            "##contig=<ID=chr1,length=248956422,assembly=GRCh38>",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##bcftools_normVersion=1.21+htslib-1.21",
        ] {
            assert!(!is_stale_provenance_line(line, "vepyr"), "{line}");
        }
    }

    fn provenance_for(config: &AnnotateVcfConfig) -> Vec<String> {
        provenance_header_lines(
            "/in/sample.vcf.gz",
            "/caches/116_GRCh38_merged",
            "parquet",
            "/out/sample.annotated.vcf",
            config,
            CacheSourceType::Merged,
        )
    }

    #[test]
    fn provenance_records_tool_identity_and_cache() {
        let config = AnnotateVcfConfig {
            provenance_tool_name: Some("vepyr".to_string()),
            provenance_tool_version: Some("0.3.0".to_string()),
            ..Default::default()
        };
        let lines = provenance_for(&config);
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("##vepyr=\"0.3.0\""), "{}", lines[0]);
        assert!(
            lines[0].contains("cache=\"/caches/116_GRCh38_merged\""),
            "{}",
            lines[0]
        );
        assert!(lines[0].contains("cache_type=\"merged\""), "{}", lines[0]);
        assert!(lines[0].contains("backend=\"parquet\""), "{}", lines[0]);
        assert!(
            lines[1].starts_with("##vepyr-command-line='"),
            "{}",
            lines[1]
        );
    }

    #[test]
    fn provenance_never_records_a_timestamp() {
        // A wall-clock stamp is what makes Ensembl VEP's own output differ
        // between two identical runs. Omitting it keeps output byte-reproducible.
        let config = AnnotateVcfConfig::default();
        for line in provenance_for(&config) {
            assert!(!line.contains("time="), "{line}");
            for year in ["202", "203"] {
                assert!(!line.contains(year), "looks like a date: {line}");
            }
        }
    }

    #[test]
    fn provenance_is_identical_across_runs() {
        let config = AnnotateVcfConfig::default();
        assert_eq!(provenance_for(&config), provenance_for(&config));
    }

    #[test]
    fn provenance_never_claims_to_be_ensembl_vep() {
        // Writing a `##VEP=` line would assert Ensembl VEP produced the file.
        let config = AnnotateVcfConfig {
            provenance_tool_name: Some("vepyr".to_string()),
            ..Default::default()
        };
        for line in provenance_for(&config) {
            assert!(!line.starts_with("##VEP"), "{line}");
        }
    }

    #[test]
    fn provenance_states_the_pinned_vep_target_when_the_cache_version_is_known() {
        let config = AnnotateVcfConfig {
            provenance_tool_name: Some("vepyr".to_string()),
            expected_cache_version: Some("116".to_string()),
            ..Default::default()
        };
        let identity = &provenance_for(&config)[0];
        assert!(identity.contains("cache_version=\"116\""), "{identity}");
        assert!(identity.contains("vep_codebase=\"116.0\""), "{identity}");
        assert!(identity.contains("api=\"116\""), "{identity}");
        assert!(identity.contains("ensembl=116.c0cf13d"), "{identity}");
        assert!(
            identity.contains("ensembl-variation=116.2fb834b"),
            "{identity}"
        );
    }

    #[test]
    fn provenance_omits_the_vep_target_when_the_cache_version_is_unknown() {
        // An unpinned run says nothing rather than guessing.
        let config = AnnotateVcfConfig::default();
        let identity = &provenance_for(&config)[0];
        assert!(!identity.contains("vep_codebase="), "{identity}");
        assert!(!identity.contains("cache_version="), "{identity}");
    }

    #[test]
    fn provenance_command_line_carries_paths_and_options() {
        let config = AnnotateVcfConfig {
            provenance_tool_name: Some("vepyr".to_string()),
            everything: true,
            hgvs: true,
            merged: true,
            reference_fasta_path: Some("/ref/GRCh38.fa".to_string()),
            ..Default::default()
        };
        let command_line = &provenance_for(&config)[1];
        for expected in [
            "/in/sample.vcf.gz",
            "/out/sample.annotated.vcf",
            "/caches/116_GRCh38_merged",
            "/ref/GRCh38.fa",
            "\"everything\":true",
            "\"hgvs\":true",
        ] {
            assert!(
                command_line.contains(expected),
                "missing {expected}: {command_line}"
            );
        }
    }

    #[test]
    fn provenance_defaults_to_the_engine_identity() {
        let config = AnnotateVcfConfig::default();
        let identity = &provenance_for(&config)[0];
        assert!(
            identity.starts_with(&format!("##{}=", env!("CARGO_PKG_NAME"))),
            "{identity}"
        );
    }

    #[test]
    fn test_csq_header_description_uses_the_ensembl_vep_prefix() {
        let config = AnnotateVcfConfig {
            everything: true,
            ..Default::default()
        };

        let description = csq_header_description(&config, CacheSourceType::Ensembl);

        assert!(
            description.starts_with("Consequence annotations from Ensembl VEP. Format: "),
            "CSQ description must match Ensembl VEP byte for byte, got: {description}"
        );
    }

    /// The prefix is fixed regardless of cache source or PICK layout — only the
    /// field list after `Format: ` varies.
    #[test]
    fn test_csq_header_description_prefix_is_stable_across_cache_sources() {
        const PREFIX: &str = "Consequence annotations from Ensembl VEP. Format: ";
        for source in [
            CacheSourceType::Ensembl,
            CacheSourceType::RefSeq,
            CacheSourceType::Merged,
        ] {
            for flag_pick_allele_gene in [false, true] {
                let config = AnnotateVcfConfig {
                    everything: true,
                    flag_pick_allele_gene,
                    ..Default::default()
                };
                let description = csq_header_description(&config, source);
                assert!(
                    description.starts_with(PREFIX),
                    "{source:?} (pick={flag_pick_allele_gene}) produced: {description}"
                );
            }
        }
    }

    #[test]
    fn test_csq_header_description_matches_vep_pick_layout() {
        let config = AnnotateVcfConfig {
            everything: true,
            flag_pick_allele_gene: true,
            ..Default::default()
        };

        let description = csq_header_description(&config, CacheSourceType::Ensembl);
        assert!(description.starts_with("Consequence annotations from Ensembl VEP. Format: "));
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
