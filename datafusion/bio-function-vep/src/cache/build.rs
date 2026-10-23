use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use datafusion::arrow::array::{
    Array, ArrayRef, BinaryArray, BooleanArray, Float64Array, Float64Builder, Int8Array,
    Int8Builder, Int64Array, StringArray, StringBuilder, UInt32Array, UInt32Builder, UInt64Array,
};
use datafusion::arrow::compute::{
    cast, concat_batches, filter_record_batch, sort_to_indices, take,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use rayon::prelude::*;
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::memory_pool::FairSpillPool;
use datafusion::execution::runtime_env::RuntimeEnvBuilder;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{ExecutionPlan, SendableRecordBatchStream};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_ensembl_cache::{
    CacheSourceType as BioFormatsCacheSourceType, EnsemblCacheOptions, EnsemblCacheTableProvider,
    EnsemblEntityKind, VEP_CACHE_REGION_SIZE_BP, build_export_query, translation_core_schema,
};
use futures::StreamExt;
use log::info;
use parquet::arrow::ArrowWriter;

const MEMORY_TOKEN_BYTES: u64 = 256 * 1024 * 1024;
const MEMORY_BUDGET_FRACTION_NUM: u64 = 80;
const MEMORY_BUDGET_FRACTION_DEN: u64 = 100;

use crate::cache::af_bundle::af_column_order;
use crate::cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest, canonical_chrom_label,
};
use crate::cache::schema::{
    VARIATION_FORBIDDEN_COLUMNS, VARIATION_OPTIONAL_COLUMNS, VARIATION_REQUIRED_COLUMNS,
    variation_projected_schema, with_cache_identity_metadata,
};
use crate::cache_builder::EntityStats;
use crate::cache_common::{
    FrequencyFields, PositionFrequency, max_global_af, select_warm_positions,
};
use crate::{
    annotate_provider::{read_compact_predictions, string_at},
    transcript_consequence::CachedPredictions,
};

const DEFAULT_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];
/// Warm/cold tier selection for the single-path variation layout. A genomic
/// `start` is "warm" (tier 0) when its max global allele frequency across
/// `minor_allele_freq`, `AF`, `gnomADg`, `gnomADe` is at least this threshold;
/// the radius extends warm membership to neighboring positions that exist in the
/// data. These must match `warm_cache::split` so the Parquet tier agrees with any
/// other tiered consumer.
const WARM_AF_THRESHOLD: f64 = 0.01;
const WARM_POSITION_RADIUS: i64 = 1;

#[derive(Debug, Clone)]
pub struct CacheBuildOptions {
    pub cache_root: String,
    pub output_dir: String,
    pub partitions: usize,
    pub cache_source_type: BioFormatsCacheSourceType,
    /// Release independently established from the raw cache and validated
    /// against the compiled vepyr support matrix before any shard is written.
    pub cache_version: String,
    pub overwrite: bool,
    /// Optional allow-list of chromosomes to build (e.g. `["chr1"]`). Matched
    /// after stripping any `chr` prefix, so `"1"` and `"chr1"` are equivalent.
    /// `None` builds every chromosome discovered in the source cache. Used to
    /// scope rebuilds (e.g. a single-chromosome profiling/parity build).
    pub chrom_filter: Option<Vec<String>>,
}

/// Build one chromosome's
/// variation shard as a single no-dictionary, page-indexed `.parquet` file under
/// `variation/<chrom>.parquet`.
///
/// Reuses the exact tiering pipeline — the warm-start pre-pass
/// (`collect_warm_starts`) and `transform_variation_tier_batch` — then applies
/// the Parquet encoding (Boolean flags, 2-array AF, `variation_name` dedup, via
/// [`crate::parquet_cache::write::encode_variation_batch`]) and streams the result
struct NativeSourceIds {
    clinvar: String,
    cosmic: String,
    dbsnp: String,
}

impl NativeSourceIds {
    fn new() -> Self {
        Self {
            clinvar: String::with_capacity(32),
            cosmic: String::with_capacity(32),
            dbsnp: String::with_capacity(32),
        }
    }

    fn clear(&mut self) {
        self.clinvar.clear();
        self.cosmic.clear();
        self.dbsnp.clear();
    }

    fn populate(&mut self, var_name: Option<&str>, var_synonyms: Option<&str>) {
        if let Some(name) = var_name {
            self.assign_by_pattern(name);
        }
        if let Some(synonyms) = var_synonyms {
            self.parse_var_synonyms(synonyms);
        }
    }

    fn assign_by_pattern(&mut self, value: &str) {
        let value = value.trim();
        if value.is_empty() || value == "." {
            return;
        }
        if value.starts_with("rs") && value[2..].chars().all(|c| c.is_ascii_digit()) {
            Self::push_value(&mut self.dbsnp, value);
        } else if value.starts_with("COSM") || value.starts_with("COSV") {
            Self::push_value(&mut self.cosmic, value);
        } else if value.starts_with("RCV") || value.starts_with("VCV") {
            Self::push_value(&mut self.clinvar, value);
        }
    }

    fn parse_var_synonyms(&mut self, synonyms: &str) {
        let mut current_slot: Option<u8> = None;
        for token in synonyms.split(',') {
            let token = token.trim();
            if token.is_empty() {
                continue;
            }

            if let Some((source, value)) = token.split_once("::") {
                current_slot = Self::lookup_source(source);
                if let Some(slot) = current_slot {
                    self.push_slot(slot, value.trim());
                } else {
                    self.assign_by_pattern(token);
                }
                continue;
            }

            if let Some((source, value)) = token.split_once(':') {
                let slot = Self::lookup_source(source);
                if let Some(s) = slot {
                    current_slot = Some(s);
                    self.push_slot(s, value.trim());
                    continue;
                }
            }

            if let Some(slot) = current_slot {
                self.push_slot(slot, token);
            } else {
                self.assign_by_pattern(token);
            }
        }
    }

    fn lookup_source(source: &str) -> Option<u8> {
        let s = source.to_ascii_lowercase();
        if s.contains("dbsnp") {
            Some(0)
        } else if s.contains("cosmic") {
            Some(1)
        } else if s.contains("clinvar") {
            Some(2)
        } else {
            None
        }
    }

    fn push_slot(&mut self, slot: u8, val: &str) {
        match slot {
            0 => Self::push_value(&mut self.dbsnp, val),
            1 => Self::push_value(&mut self.cosmic, val),
            2 => Self::push_value(&mut self.clinvar, val),
            _ => {}
        }
    }

    fn push_value(buf: &mut String, val: &str) {
        let val = val.trim();
        if val.is_empty() || val == "." {
            return;
        }
        for existing in buf.split(',') {
            if existing == val {
                return;
            }
        }
        if !buf.is_empty() {
            buf.push(',');
        }
        buf.push_str(val);
    }
}

fn intermediate_variation_schema(has_clin_sig_ref_allele: bool) -> SchemaRef {
    let mut fields = vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new("allele_string", DataType::Utf8, false),
        Field::new("failed", DataType::Int8, true),
        Field::new("variation_name", DataType::Utf8, false),
        Field::new("clin_sig", DataType::Utf8, true),
        Field::new("clin_sig_allele", DataType::Utf8, true),
    ];
    if has_clin_sig_ref_allele {
        fields.push(Field::new("clin_sig_ref_allele", DataType::Utf8, true));
    }
    fields.extend(vec![
        Field::new("clinical_impact", DataType::Utf8, true),
        Field::new("phenotype_or_disease", DataType::Int8, true),
        Field::new("pubmed", DataType::Utf8, true),
        Field::new("somatic", DataType::Int8, true),
        Field::new("minor_allele", DataType::Utf8, true),
        Field::new("minor_allele_freq", DataType::Float64, true),
    ]);
    for name in af_column_order() {
        fields.push(Field::new(name, DataType::Utf8, true));
    }
    fields.extend(vec![
        Field::new("clinvar_ids", DataType::Utf8, true),
        Field::new("cosmic_ids", DataType::Utf8, true),
        Field::new("dbsnp_ids", DataType::Utf8, true),
        Field::new("tier", DataType::Int8, false),
    ]);
    Arc::new(Schema::new(fields))
}

struct NativeVariationBatchBuilder {
    chrom: StringBuilder,
    start: UInt32Builder,
    end: UInt32Builder,
    allele_string: StringBuilder,
    failed: Int8Builder,
    variation_name: StringBuilder,
    clin_sig: StringBuilder,
    clin_sig_allele: StringBuilder,
    clin_sig_ref_allele: StringBuilder,
    clinical_impact: StringBuilder,
    phenotype_or_disease: Int8Builder,
    pubmed: StringBuilder,
    somatic: Int8Builder,
    minor_allele: StringBuilder,
    minor_allele_freq: Float64Builder,
    af_builders: Vec<StringBuilder>,
    clinvar_ids: StringBuilder,
    cosmic_ids: StringBuilder,
    dbsnp_ids: StringBuilder,
    tier: Int8Builder,
    pub num_rows: usize,
}

impl NativeVariationBatchBuilder {
    fn new(capacity: usize) -> Self {
        let mut af_builders = Vec::with_capacity(27);
        for _ in 0..27 {
            af_builders.push(StringBuilder::with_capacity(capacity, capacity * 2));
        }
        Self {
            chrom: StringBuilder::with_capacity(capacity, capacity * 2),
            start: UInt32Builder::with_capacity(capacity),
            end: UInt32Builder::with_capacity(capacity),
            allele_string: StringBuilder::with_capacity(capacity, capacity * 4),
            failed: Int8Builder::with_capacity(capacity),
            variation_name: StringBuilder::with_capacity(capacity, capacity * 12),
            clin_sig: StringBuilder::with_capacity(capacity, capacity * 2),
            clin_sig_allele: StringBuilder::with_capacity(capacity, capacity * 2),
            clin_sig_ref_allele: StringBuilder::with_capacity(capacity, capacity * 2),
            clinical_impact: StringBuilder::with_capacity(capacity, capacity * 2),
            phenotype_or_disease: Int8Builder::with_capacity(capacity),
            pubmed: StringBuilder::with_capacity(capacity, capacity * 4),
            somatic: Int8Builder::with_capacity(capacity),
            minor_allele: StringBuilder::with_capacity(capacity, capacity * 2),
            minor_allele_freq: Float64Builder::with_capacity(capacity),
            af_builders,
            clinvar_ids: StringBuilder::with_capacity(capacity, capacity * 4),
            cosmic_ids: StringBuilder::with_capacity(capacity, capacity * 4),
            dbsnp_ids: StringBuilder::with_capacity(capacity, capacity * 12),
            tier: Int8Builder::with_capacity(capacity),
            num_rows: 0,
        }
    }

    #[inline]
    fn add_row(
        &mut self,
        chrom: &str,
        fields: &[&str],
        tier_val: i8,
        ids: &NativeSourceIds,
        has_clin_sig_ref: bool,
    ) {
        self.chrom.append_value(chrom);

        let start_val: u32 = fields.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
        let end_val: u32 = fields
            .get(5)
            .and_then(|s| if *s == "." || s.is_empty() { None } else { s.parse().ok() })
            .unwrap_or(start_val);
        self.start.append_value(start_val);
        self.end.append_value(end_val);

        self.allele_string.append_value(fields.get(6).copied().unwrap_or("."));

        if fields.get(2) == Some(&"1") {
            self.failed.append_value(1);
        } else {
            self.failed.append_null();
        }

        self.variation_name.append_value(fields.get(1).copied().unwrap_or("."));

        append_opt_str(&mut self.clin_sig, fields.get(8).copied());
        append_opt_str(&mut self.clin_sig_allele, fields.get(10).copied());
        if has_clin_sig_ref {
            append_opt_str(&mut self.clin_sig_ref_allele, fields.get(11).copied());
        }
        append_opt_str(&mut self.clinical_impact, fields.get(12).copied());

        if fields.get(9) == Some(&"1") {
            self.phenotype_or_disease.append_value(1);
        } else {
            self.phenotype_or_disease.append_null();
        }

        append_opt_str(&mut self.pubmed, fields.get(13).copied());

        if fields.get(3) == Some(&"1") {
            self.somatic.append_value(1);
        } else {
            self.somatic.append_null();
        }

        self.minor_allele.append_null();
        self.minor_allele_freq.append_null();

        const AF_INDICES: [usize; 27] = [
            15, 16, 17, 18, 19, 20,
            21, 22, 23, 24, 25, 26, 27, 29, 30, 28,
            31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 40,
        ];
        for (i, &tab_idx) in AF_INDICES.iter().enumerate() {
            append_opt_str(&mut self.af_builders[i], fields.get(tab_idx).copied());
        }

        if ids.clinvar.is_empty() {
            self.clinvar_ids.append_null();
        } else {
            self.clinvar_ids.append_value(&ids.clinvar);
        }
        if ids.cosmic.is_empty() {
            self.cosmic_ids.append_null();
        } else {
            self.cosmic_ids.append_value(&ids.cosmic);
        }
        if ids.dbsnp.is_empty() {
            self.dbsnp_ids.append_null();
        } else {
            self.dbsnp_ids.append_value(&ids.dbsnp);
        }

        self.tier.append_value(tier_val);
        self.num_rows += 1;
    }

    fn finish(mut self, schema: &SchemaRef, has_clin_sig_ref: bool) -> RecordBatch {
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(self.chrom.finish()),
            Arc::new(self.start.finish()),
            Arc::new(self.end.finish()),
            Arc::new(self.allele_string.finish()),
            Arc::new(self.failed.finish()),
            Arc::new(self.variation_name.finish()),
            Arc::new(self.clin_sig.finish()),
            Arc::new(self.clin_sig_allele.finish()),
        ];
        if has_clin_sig_ref {
            columns.push(Arc::new(self.clin_sig_ref_allele.finish()));
        }
        columns.extend(vec![
            Arc::new(self.clinical_impact.finish()) as ArrayRef,
            Arc::new(self.phenotype_or_disease.finish()),
            Arc::new(self.pubmed.finish()),
            Arc::new(self.somatic.finish()),
            Arc::new(self.minor_allele.finish()),
            Arc::new(self.minor_allele_freq.finish()),
        ]);
        for mut b in self.af_builders {
            columns.push(Arc::new(b.finish()));
        }
        columns.extend(vec![
            Arc::new(self.clinvar_ids.finish()) as ArrayRef,
            Arc::new(self.cosmic_ids.finish()),
            Arc::new(self.dbsnp_ids.finish()),
            Arc::new(self.tier.finish()),
        ]);
        RecordBatch::try_new(Arc::clone(schema), columns).expect("valid batch")
    }
}

#[inline]
fn append_opt_str(builder: &mut StringBuilder, val: Option<&str>) {
    match val {
        Some(s) if !s.is_empty() && s != "." => builder.append_value(s),
        _ => builder.append_null(),
    }
}

async fn build_parquet_variation_chrom_native(
    options: &CacheBuildOptions,
    chrom: &str,
    all_vars_path: &Path,
) -> Result<ChromDatasetEntry> {
    use crate::parquet_cache::write::variation_output_schema;

    let entity_dir = PathBuf::from(&options.output_dir).join("variation");
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    remove_existing_parquet_shard(&shard_path, options.overwrite)?;

    let source_type = options.cache_source_type.as_str().to_string();
    let provider_opts = EnsemblCacheOptions::new(&options.cache_root)
        .with_cache_source_type(options.cache_source_type)
        .with_expected_cache_version(&options.cache_version);
    let provider =
        EnsemblCacheTableProvider::for_entity(EnsemblEntityKind::Variation, provider_opts)?;
    let source_schema = provider.schema();
    let projected = variation_output_schema(source_schema.as_ref(), &source_type)?;
    let out_schema = Arc::new(with_cache_identity_metadata(
        projected.as_ref(),
        &source_type,
        &options.cache_version,
    ));
    let has_clin_sig_ref = source_schema
        .column_with_name("clin_sig_ref_allele")
        .is_some();
    let intermediate_schema = intermediate_variation_schema(has_clin_sig_ref);

    let all_vars_path_buf = all_vars_path.to_path_buf();
    let chrom_string = chrom.to_string();

    tokio::task::spawn_blocking(move || {
        build_parquet_variation_chrom_native_blocking(
            &chrom_string,
            &all_vars_path_buf,
            &shard_path,
            &manifest_chrom,
            &file_name,
            out_schema,
            intermediate_schema,
            has_clin_sig_ref,
        )
    })
    .await
    .map_err(|err| DataFusionError::Execution(format!("spawn_blocking join: {err}")))?
}

fn build_parquet_variation_chrom_native_blocking(
    chrom: &str,
    all_vars_path: &Path,
    shard_path: &Path,
    manifest_chrom: &str,
    file_name: &str,
    out_schema: SchemaRef,
    intermediate_schema: SchemaRef,
    has_clin_sig_ref: bool,
) -> Result<ChromDatasetEntry> {
    use crate::parquet_cache::write::{VariationParquetShardWriter, encode_variation_batch};
    use std::io::BufRead;

    let t_start = std::time::Instant::now();
    // Experimental low-memory mode for distributed workers.  It keeps the
    // native single-pass scanner but disables warm-tier classification, so no
    // encoded warm batches accumulate in RAM.  This is benchmark-only until
    // read-side tier semantics are validated for production.
    let warm_starts = if std::env::var("VEPYR_DISABLE_WARM").as_deref() == Ok("1") {
        Arc::new(std::collections::HashSet::new())
    } else {
        Arc::new(collect_warm_starts_native(all_vars_path)?)
    };

    let file_size = std::fs::metadata(all_vars_path)
        .map_err(|e| DataFusionError::Execution(format!("stat {}: {e}", all_vars_path.display())))?
        .len();
    let requested_external_count = std::env::var("VEPYR_VAR_RANGE_COUNT")
        .ok()
        .map(|value| {
            value.parse::<usize>().map_err(|err| {
                DataFusionError::Execution(format!("invalid VEPYR_VAR_RANGE_COUNT: {err}"))
            })
        })
        .transpose()?;
    let num_partitions = if let Some(count) = requested_external_count {
        if count == 0 {
            return Err(DataFusionError::Execution(
                "VEPYR_VAR_RANGE_COUNT must be positive".to_string(),
            ));
        }
        count
    } else if warm_starts.is_empty() {
        // Low-memory mode must consume the cold stream while the scanner runs;
        // one partition gives a bounded producer/consumer path without a large
        // cross-partition queue.
        1
    } else {
        std::env::var("VEPYR_VAR_PARTITIONS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or_else(|| ((file_size / (40 * 1024 * 1024)).max(2).min(6)) as usize)
    };
    let partitions = compute_bgzf_byte_partitions(all_vars_path, num_partitions)?;
    // Experimental external distributed mode: one container handles exactly
    // one BGZF range. The caller launches separate synchronous workers, so
    // this does not add parallelism inside a worker. Each worker writes its
    // shard into an isolated output volume; the orchestrator may assemble the
    // range shards afterwards.
    let partitions = match (
        std::env::var("VEPYR_VAR_RANGE_INDEX").ok(),
        std::env::var("VEPYR_VAR_RANGE_COUNT").ok(),
    ) {
        (Some(index), Some(count)) => {
            let index = index.parse::<usize>().map_err(|err| {
                DataFusionError::Execution(format!("invalid VEPYR_VAR_RANGE_INDEX: {err}"))
            })?;
            let count = count.parse::<usize>().map_err(|err| {
                DataFusionError::Execution(format!("invalid VEPYR_VAR_RANGE_COUNT: {err}"))
            })?;
            if count == 0 || index >= count {
                return Err(DataFusionError::Execution(format!(
                    "invalid external BGZF range {index}/{count}"
                )));
            }
            if count != num_partitions {
                return Err(DataFusionError::Execution(format!(
                    "VEPYR_VAR_RANGE_COUNT={count} must equal computed partition count {num_partitions}"
                )));
            }
            vec![partitions[index].clone()]
        }
        (None, None) => partitions,
        _ => {
            return Err(DataFusionError::Execution(
                "VEPYR_VAR_RANGE_INDEX and VEPYR_VAR_RANGE_COUNT must be set together"
                    .to_string(),
            ));
        }
    };

    let mut writer = VariationParquetShardWriter::create(shard_path, Arc::clone(&out_schema))?;

    // Streaming channels for cold tier batches (one per partition)
    let (senders, receivers): (Vec<_>, Vec<_>) = (0..partitions.len())
        .map(|_| std::sync::mpsc::channel::<RecordBatch>())
        .unzip();

    let all_vars_path_buf = all_vars_path.to_path_buf();
    let warm_starts_clone = Arc::clone(&warm_starts);
    let intermediate_schema_clone = Arc::clone(&intermediate_schema);
    let out_schema_clone = Arc::clone(&out_schema);
    let chrom_string = chrom.to_string();

    // Single-Pass Rayon Parallel Scan:
    // Decompresses BGZF once. Rows with start in warm_starts are buffered as warm batches
    // in memory; rows with start NOT in warm_starts are streamed across channels as cold batches.
    let rayon_handle = std::thread::spawn(move || {
        partitions
            .into_par_iter()
            .zip(senders.into_par_iter())
            .enumerate()
            .map(|(part_idx, (part, tx))| {
                let file = match std::fs::File::open(&all_vars_path_buf) {
                    Ok(f) => f,
                    Err(_) => return (part_idx, Vec::new()),
                };
                let mut reader = noodles_bgzf::io::Reader::new(
                    std::io::BufReader::with_capacity(BGZF_IO_BUFFER_SIZE, file),
                );
                if part.start_compressed > 0 {
                    if let Ok(vpos) =
                        noodles_bgzf::VirtualPosition::try_from((part.start_compressed, 0))
                    {
                        let _ = reader.seek(vpos);
                        if part.skip_first_line {
                            let mut discard = String::new();
                            let _ = reader.read_line(&mut discard);
                        }
                    }
                }

                let mut warm_builder = NativeVariationBatchBuilder::new(50_000);
                let mut cold_builder = NativeVariationBatchBuilder::new(100_000);
                let mut warm_batches: Vec<RecordBatch> = Vec::new();
                let mut ids = NativeSourceIds::new();
                let mut line = String::new();

                loop {
                    if reader.virtual_position().compressed() >= part.end_compressed {
                        break;
                    }
                    line.clear();
                    match reader.read_line(&mut line) {
                        Ok(0) => break,
                        Ok(_) => {
                            let trimmed = line.trim();
                            if trimmed.is_empty() || trimmed.starts_with('#') {
                                continue;
                            }
                            let mut fields = [""; 48];
                            let mut count = 0;
                            for f in trimmed.split('\t') {
                                if count < 48 {
                                    fields[count] = f;
                                    count += 1;
                                }
                            }
                            let start_val: i64 = fields[4].parse().unwrap_or(0);
                            ids.clear();
                            ids.populate(fields.get(1).copied(), fields.get(14).copied());

                            if warm_starts_clone.contains(&start_val) {
                                warm_builder.add_row(&chrom_string, &fields[..count], 0, &ids, has_clin_sig_ref);
                                if warm_builder.num_rows >= 100_000 {
                                    let batch = warm_builder.finish(&intermediate_schema_clone, has_clin_sig_ref);
                                    if let Ok(encoded) = encode_variation_batch(&batch, &out_schema_clone) {
                                        warm_batches.push(encoded);
                                    }
                                    warm_builder = NativeVariationBatchBuilder::new(100_000);
                                }
                            } else {
                                cold_builder.add_row(&chrom_string, &fields[..count], 1, &ids, has_clin_sig_ref);
                                if cold_builder.num_rows >= 100_000 {
                                    let batch = cold_builder.finish(&intermediate_schema_clone, has_clin_sig_ref);
                                    if let Ok(encoded) = encode_variation_batch(&batch, &out_schema_clone) {
                                        let _ = tx.send(encoded);
                                    }
                                    cold_builder = NativeVariationBatchBuilder::new(100_000);
                                }
                            }
                        }
                        Err(_) => break,
                    }
                }

                if warm_builder.num_rows > 0 {
                    let batch = warm_builder.finish(&intermediate_schema_clone, has_clin_sig_ref);
                    if let Ok(encoded) = encode_variation_batch(&batch, &out_schema_clone) {
                        warm_batches.push(encoded);
                    }
                }

                if cold_builder.num_rows > 0 {
                    let batch = cold_builder.finish(&intermediate_schema_clone, has_clin_sig_ref);
                    if let Ok(encoded) = encode_variation_batch(&batch, &out_schema_clone) {
                        let _ = tx.send(encoded);
                    }
                }

                (part_idx, warm_batches)
            })
            .collect::<Vec<_>>()
    });

    if warm_starts.is_empty() {
        // In low-memory mode the scanner sends cold batches through one channel.
        // Consume it before joining the scanner; joining first would retain every
        // batch in an unbounded channel and recreate the OOM condition.
        while let Ok(batch) = receivers[0].recv() {
            writer.write(&batch)?;
        }
        rayon_handle.join().map_err(|_| {
            DataFusionError::Execution("rayon partition scan panicked".to_string())
        })?;
        let finished_rows = writer.finish()?;
        eprintln!(
            "[VAR_TIMING] chrom {} native low-memory build finished in {:.2}s ({} rows)",
            chrom,
            t_start.elapsed().as_secs_f64(),
            finished_rows
        );
        return Ok(ChromDatasetEntry::new(
            manifest_chrom,
            file_name,
            finished_rows,
        ));
    }

    let mut warm_results = rayon_handle.join().map_err(|_| {
        DataFusionError::Execution("rayon partition scan panicked".to_string())
    })?;

    // 1. Write warm tier batches in partition order
    warm_results.sort_by_key(|(idx, _)| *idx);
    for (_, batches) in warm_results {
        for b in batches {
            writer.write(&b)?;
        }
    }

    // 2. Write cold tier batches in partition order
    for rx in receivers {
        while let Ok(b) = rx.recv() {
            writer.write(&b)?;
        }
    }

    let finished_rows = writer.finish()?;
    eprintln!(
        "[VAR_TIMING] chrom {} native single-pass build finished in {:.2}s ({} rows)",
        chrom,
        t_start.elapsed().as_secs_f64(),
        finished_rows
    );
    Ok(ChromDatasetEntry::new(
        manifest_chrom,
        file_name,
        finished_rows,
    ))
}

/// to the shard. Rows are written warm tier (0) first then cold (1), each
/// `start`-ordered by the source `ORDER BY start` plan — the two sorted runs the
/// read-side `PageDir` binary search relies on.
pub async fn build_parquet_variation_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    use crate::parquet_cache::write::{
        VariationParquetShardWriter, encode_variation_batch, variation_output_schema,
    };

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let all_vars_path = Path::new(&options.cache_root).join(source_chrom).join("all_vars.gz");
    let use_df = std::env::var("VEPYR_USE_DATAFUSION_VARIATION").as_deref() == Ok("1");
    if all_vars_path.exists() && !use_df {
        return build_parquet_variation_chrom_native(options, chrom, &all_vars_path).await;
    }

    let entity_dir = PathBuf::from(&options.output_dir).join("variation");
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let query = build_export_query(
        EnsemblEntityKind::Variation,
        "var",
        Some(source_chrom),
        None,
    );
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    if shard_path.exists() {
        if !options.overwrite {
            return Err(DataFusionError::Execution(format!(
                "parquet variation shard '{}' already exists and overwrite=false",
                shard_path.display()
            )));
        }
        std::fs::remove_file(&shard_path).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to remove existing parquet shard '{}': {err}",
                shard_path.display()
            ))
        })?;
    }

    let source_type = options.cache_source_type.as_str().to_string();
    let t_warm_start = std::time::Instant::now();
    let warm_starts = Arc::new(collect_warm_starts(options, source_chrom).await?);
    eprintln!(
        "[VAR_TIMING] chrom {} warm_starts took {:.2}s ({} positions)",
        chrom,
        t_warm_start.elapsed().as_secs_f64(),
        warm_starts.len()
    );
    info!(
        "Parquet variation {} warm positions={} (af>={WARM_AF_THRESHOLD}, +/-{WARM_POSITION_RADIUS})",
        shard_path.display(),
        warm_starts.len(),
    );

    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Variation, "var")?;

    // Derive the output schema once from the source plan schema and reuse it for
    // every written batch so the shard carries one stable Arrow schema.
    let schema_plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let projected = variation_output_schema(schema_plan.schema().as_ref(), &source_type)?;
    let out_schema = Arc::new(with_cache_identity_metadata(
        projected.as_ref(),
        &source_type,
        &options.cache_version,
    ));

    // Single-pass execution: stream the DataFusion query once instead of multiple passes.
    // Warm tier batches (rare: ~3-5% of variants, e.g. ~30-40 MB RAM for chr20) are buffered
    // in memory, while cold tier batches are streamed directly to an Arrow IPC temp file.
    // Finally, warm batches are written first to the output Parquet, followed by cold batches.
    let t_stream_start = std::time::Instant::now();
    let plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

    let rows = if warm_starts.is_empty() {
        // No warm variants: write directly to final parquet shard.
        let mut writer = VariationParquetShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
        while let Some(batch) = stream.next().await {
            let (_, cold_opt) =
                transform_variation_split_batch(batch?, &source_type, &warm_starts)?;
            if let Some(cold) = cold_opt {
                let encoded = encode_variation_batch(&cold, &out_schema)?;
                writer.write(&encoded)?;
            }
        }
        writer.finish()?
    } else {
        let mut warm_batches: Vec<RecordBatch> = Vec::new();
        let mut cold_batches: Vec<RecordBatch> = Vec::new();

        while let Some(batch) = stream.next().await {
            let (warm_opt, cold_opt) =
                transform_variation_split_batch(batch?, &source_type, &warm_starts)?;
            if let Some(warm) = warm_opt {
                let encoded = encode_variation_batch(&warm, &out_schema)?;
                warm_batches.push(encoded);
            }
            if let Some(cold) = cold_opt {
                let encoded = encode_variation_batch(&cold, &out_schema)?;
                cold_batches.push(encoded);
            }
        }
        eprintln!(
            "[VAR_TIMING] chrom {} DataFusion stream finished in {:.2}s (warm batches: {}, cold batches: {})",
            chrom,
            t_stream_start.elapsed().as_secs_f64(),
            warm_batches.len(),
            cold_batches.len()
        );

        let t_parquet_start = std::time::Instant::now();
        let mut writer = VariationParquetShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
        for warm in warm_batches {
            writer.write(&warm)?;
        }
        for cold in cold_batches {
            writer.write(&cold)?;
        }
        let finished_rows = writer.finish()?;
        eprintln!(
            "[VAR_TIMING] chrom {} in-memory Parquet write took {:.2}s (total {} rows)",
            chrom,
            t_parquet_start.elapsed().as_secs_f64(),
            finished_rows
        );
        finished_rows
    };
    if rows == 0 {
        // The writer is created eagerly (before the stream), so a contig with no
        // variants — e.g. a GL*/KI* scaffold present in the source schema but
        // carrying no variant rows — leaves a schema-only shard. Drop it so the
        // on-disk layout matches the manifest (which omits empty contigs), like
        // the lazily-written context/translation shards.
        let _ = std::fs::remove_file(&shard_path);
    }
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
}

const BGZF_IO_BUFFER_SIZE: usize = 64 * 1024;

#[derive(Debug, Clone)]
struct BgzfBytePartition {
    start_compressed: u64,
    end_compressed: u64,
    skip_first_line: bool,
}

fn scan_bgzf_block_offsets(file_path: &Path, file_size: u64) -> Result<Vec<u64>> {
    use std::io::{Read, Seek, SeekFrom};
    let mut file = std::fs::File::open(file_path)
        .map_err(|e| DataFusionError::Execution(format!("open {}: {e}", file_path.display())))?;
    let mut offsets = Vec::new();
    let mut pos = 0u64;
    let mut header = [0u8; 18];

    while pos + 28 <= file_size {
        file.seek(SeekFrom::Start(pos))
            .map_err(|e| DataFusionError::Execution(format!("seek: {e}")))?;
        if file.read_exact(&mut header).is_err() {
            break;
        }
        if header[0] != 0x1f || header[1] != 0x8b || header[2] != 0x08 || header[3] != 0x04 {
            break;
        }
        let bsize = u16::from_le_bytes([header[16], header[17]]) as u64 + 1;
        let is_eof_marker = bsize == 28 && pos + bsize >= file_size;
        if !is_eof_marker {
            offsets.push(pos);
        }
        pos += bsize;
    }
    Ok(offsets)
}

fn compute_bgzf_byte_partitions(file_path: &Path, num_partitions: usize) -> Result<Vec<BgzfBytePartition>> {
    let file_size = std::fs::metadata(file_path)
        .map_err(|e| DataFusionError::Execution(format!("stat {}: {e}", file_path.display())))?
        .len();
    if file_size == 0 || num_partitions <= 1 {
        return Ok(vec![BgzfBytePartition {
            start_compressed: 0,
            end_compressed: file_size,
            skip_first_line: false,
        }]);
    }

    let block_offsets = scan_bgzf_block_offsets(file_path, file_size)?;
    if block_offsets.len() <= 1 {
        return Ok(vec![BgzfBytePartition {
            start_compressed: 0,
            end_compressed: file_size,
            skip_first_line: false,
        }]);
    }

    let target_chunk = file_size / num_partitions as u64;
    let mut boundaries = vec![0u64];
    let mut current_start = 0u64;

    for i in 1..num_partitions {
        let target_offset = target_chunk * i as u64;
        let best = block_offsets
            .iter()
            .copied()
            .filter(|&off| off > current_start)
            .min_by_key(|&off| (off as i64 - target_offset as i64).unsigned_abs())
            .unwrap_or(file_size);

        if best < file_size && best > current_start {
            boundaries.push(best);
            current_start = best;
        }
    }

    let mut partitions = Vec::with_capacity(boundaries.len());
    for i in 0..boundaries.len() {
        let start = boundaries[i];
        let end = if i + 1 < boundaries.len() {
            boundaries[i + 1]
        } else {
            file_size
        };

        let skip_first_line = if start > 0 {
            let prev_block_offset = block_offsets
                .iter()
                .copied()
                .filter(|&off| off < start)
                .last()
                .unwrap_or(0);
            !bgzf_block_ends_with_newline(file_path, prev_block_offset)?
        } else {
            false
        };

        partitions.push(BgzfBytePartition {
            start_compressed: start,
            end_compressed: end,
            skip_first_line,
        });
    }

    Ok(partitions)
}

fn bgzf_block_ends_with_newline(file_path: &Path, block_offset: u64) -> Result<bool> {
    use std::io::Read;
    let file = std::fs::File::open(file_path)
        .map_err(|e| DataFusionError::Execution(format!("open: {e}")))?;
    let mut reader = noodles_bgzf::io::Reader::new(std::io::BufReader::with_capacity(BGZF_IO_BUFFER_SIZE, file));
    let vpos = noodles_bgzf::VirtualPosition::try_from((block_offset, 0))
        .map_err(|e| DataFusionError::Execution(format!("vpos: {e}")))?;
    reader.seek(vpos).map_err(|e| DataFusionError::Execution(format!("seek: {e}")))?;

    let mut buf = [0u8; 64 * 1024];
    let mut last_byte = None;
    loop {
        let n = reader.read(&mut buf).map_err(|e| DataFusionError::Execution(format!("read: {e}")))?;
        if n == 0 {
            break;
        }
        last_byte = Some(buf[n - 1]);
        if reader.virtual_position().compressed() != block_offset {
            break;
        }
    }
    Ok(last_byte == Some(b'\n'))
}

fn collect_warm_starts_native(path: &Path) -> Result<HashSet<i64>> {
    use rayon::prelude::*;
    use std::io::BufRead;

    let num_threads = std::env::var("VEPYR_VAR_PARTITIONS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(6);
    let partitions = compute_bgzf_byte_partitions(path, num_threads)?;

    let results: Vec<(Vec<i64>, Vec<i64>)> = partitions
        .into_par_iter()
        .map(|part| {
            let file = std::fs::File::open(path).expect("open file");
            let mut reader = noodles_bgzf::io::Reader::new(std::io::BufReader::with_capacity(BGZF_IO_BUFFER_SIZE, file));
            if part.start_compressed > 0 {
                let vpos = noodles_bgzf::VirtualPosition::try_from((part.start_compressed, 0)).expect("vpos");
                reader.seek(vpos).expect("seek");
                if part.skip_first_line {
                    let mut discard = String::new();
                    let _ = reader.read_line(&mut discard);
                }
            }

            let mut all_starts = Vec::with_capacity(3_000_000);
            let mut high_af_starts = Vec::with_capacity(100_000);
            let mut line = String::new();

            loop {
                if reader.virtual_position().compressed() >= part.end_compressed {
                    break;
                }
                line.clear();
                match reader.read_line(&mut line) {
                    Ok(0) => break,
                    Ok(_) => {
                        let trimmed = line.trim();
                        if trimmed.is_empty() || trimmed.starts_with('#') {
                            continue;
                        }
                        let mut tab_idx = 0;
                        let mut start_val: i64 = 0;
                        let mut max_af = 0.0f64;

                        for field in trimmed.split('\t') {
                            match tab_idx {
                                4 => {
                                    if let Ok(v) = field.parse::<i64>() {
                                        start_val = v;
                                    }
                                }
                                10 => {
                                    if let Ok(v) = field.parse::<f64>() {
                                        if v.is_finite() && v > max_af {
                                            max_af = v;
                                        }
                                    }
                                }
                                15 | 21 | 31 => {
                                    let v = crate::cache_common::max_af_from_pairs(Some(field));
                                    if v > max_af {
                                        max_af = v;
                                    }
                                }
                                _ => {}
                            }
                            tab_idx += 1;
                            if tab_idx > 31 {
                                break;
                            }
                        }

                        if start_val > 0 {
                            all_starts.push(start_val);
                            if max_af >= WARM_AF_THRESHOLD {
                                high_af_starts.push(start_val);
                            }
                        }
                    }
                    Err(_) => break,
                }
            }
            (all_starts, high_af_starts)
        })
        .collect();

    let mut source_positions: HashSet<i64> = HashSet::new();
    let mut high_af: Vec<i64> = Vec::new();

    for (starts, highs) in results {
        source_positions.extend(starts);
        high_af.extend(highs);
    }

    let mut warm_set: HashSet<i64> = HashSet::new();
    for start in high_af {
        for offset in -WARM_POSITION_RADIUS..=WARM_POSITION_RADIUS {
            let cand = start + offset;
            if source_positions.contains(&cand) {
                warm_set.insert(cand);
            }
        }
    }

    Ok(warm_set)
}

/// Run a cheap pre-pass over the source variation table for one chromosome to
/// classify each genomic `start` as warm (common) or cold (rare) using the same
/// allele-frequency selector as the Parquet warm tier.
async fn collect_warm_starts(
    options: &CacheBuildOptions,
    source_chrom: &str,
) -> Result<HashSet<i64>> {
    let all_vars_path = Path::new(&options.cache_root).join(source_chrom).join("all_vars.gz");
    if all_vars_path.exists() {
        return collect_warm_starts_native(&all_vars_path);
    }

    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Variation, "var")?;
    // Mixed-case columns must be quoted; DataFusion lowercases bare identifiers.
    let sql = format!(
        "SELECT start, minor_allele_freq, \"AF\", \"gnomADg\", \"gnomADe\" FROM var WHERE chrom = '{}'",
        sql_escape_literal(source_chrom)
    );
    let df = ctx.sql(&sql).await?;
    let mut stream = df.execute_stream().await?;
    let mut rows = Vec::new();
    while let Some(batch) = stream.next().await {
        collect_position_frequencies(&batch?, &mut rows)?;
    }
    let warm_set = select_warm_positions(
        rows,
        WARM_AF_THRESHOLD,
        WARM_POSITION_RADIUS,
    );
    Ok(warm_set.into_iter().collect())
}

fn collect_position_frequencies(
    batch: &RecordBatch,
    out: &mut Vec<PositionFrequency>,
) -> Result<()> {
    if batch.num_rows() == 0 {
        return Ok(());
    }
    let schema = batch.schema();
    let start_idx = schema.index_of("start")?;
    let start = cast(batch.column(start_idx).as_ref(), &DataType::Int64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation start to Int64: {err}"))
    })?;
    let start = start
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| DataFusionError::Execution("variation start must be Int64".to_string()))?;

    let maf = column_as_f64(batch, "minor_allele_freq")?;
    let af = column_as_str(batch, "AF");
    let gnomadg = column_as_str(batch, "gnomADg");
    let gnomade = column_as_str(batch, "gnomADe");

    out.reserve(batch.num_rows());
    for row in 0..batch.num_rows() {
        if start.is_null(row) {
            continue;
        }
        let af_s = af.as_ref().and_then(|arr| string_at(arr.as_ref(), row));
        let gnomadg_s = gnomadg
            .as_ref()
            .and_then(|arr| string_at(arr.as_ref(), row));
        let gnomade_s = gnomade
            .as_ref()
            .and_then(|arr| string_at(arr.as_ref(), row));
        let fields = FrequencyFields {
            minor_allele_freq: maf
                .as_ref()
                .and_then(|arr| (!arr.is_null(row)).then(|| arr.value(row))),
            af: af_s.as_deref(),
            gnomadg: gnomadg_s.as_deref(),
            gnomade: gnomade_s.as_deref(),
        };
        out.push(PositionFrequency {
            start: start.value(row),
            max_af: max_global_af(&fields),
        });
    }
    Ok(())
}

fn column_as_f64(batch: &RecordBatch, name: &str) -> Result<Option<Float64Array>> {
    let Ok(index) = batch.schema().index_of(name) else {
        return Ok(None);
    };
    let array = cast(batch.column(index).as_ref(), &DataType::Float64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation {name} to Float64: {err}"))
    })?;
    let array = array
        .as_any()
        .downcast_ref::<Float64Array>()
        .ok_or_else(|| DataFusionError::Execution(format!("variation {name} must be Float64")))?
        .clone();
    Ok(Some(array))
}

fn column_as_str(batch: &RecordBatch, name: &str) -> Option<ArrayRef> {
    batch
        .schema()
        .index_of(name)
        .ok()
        .map(|index| batch.column(index).clone())
}

/// Build, for one chromosome,
/// the position-sliced SIFT/PolyPhen shard as a single no-dictionary,
/// page-indexed `.parquet` file under `translation_sift/<chrom>.parquet`,
/// physically sorted by the u64 `key` so the read-side `PageDir` sees one monotone
/// run.
///
/// Reuses the exact row pipeline — the shared `transcript_uid` map and
/// `transform_translation_sift_position_batch` — then globally sorts the
/// accumulated `(key, sift, poly)` rows by `key` (the source dedup plan groups by
/// transcript, not by the derived key) and writes via
/// [`crate::parquet_cache::sift::write_sift_parquet`]. The schema metadata
/// (`SIFT_BLOB_VERSION`) is preserved so the reader decodes blobs unchanged.
/// Fused translation build for one chromosome: executes the DataFusion deduplication
/// window query once over the region files and streams batches simultaneously to
/// translation_core (projected context shard) and translation_sift (position-sliced
/// SIFT/PolyPhen blobs).
pub async fn build_parquet_translation_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
    build_core: bool,
    build_sift: bool,
) -> Result<(Option<ChromDatasetEntry>, Option<ChromDatasetEntry>)> {
    if !build_core && !build_sift {
        return Ok((None, None));
    }

    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let source_type = options.cache_source_type.as_str().to_string();

    let (core_shard_path, core_target_schema) = if build_core {
        let core_dir = PathBuf::from(&options.output_dir).join("translation_core");
        std::fs::create_dir_all(&core_dir).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to create '{}': {err}",
                core_dir.display()
            ))
        })?;
        let shard_path = core_dir.join(&file_name);
        remove_existing_parquet_shard(&shard_path, options.overwrite)?;
        let schema = translation_core_schema(false, options.cache_source_type, &options.cache_version);
        (Some(shard_path), Some(schema))
    } else {
        (None, None)
    };

    let (sift_shard_path, uid_map) = if build_sift {
        let sift_dir = PathBuf::from(&options.output_dir).join("translation_sift");
        std::fs::create_dir_all(&sift_dir).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to create '{}': {err}",
                sift_dir.display()
            ))
        })?;
        let shard_path = sift_dir.join(&file_name);
        remove_existing_parquet_shard(&shard_path, options.overwrite)?;
        let map = Arc::new(load_transcript_uid_map(options, source_chrom).await?);
        (Some(shard_path), Some(map))
    } else {
        (None, None)
    };

    let base_query = build_translation_dedup_query_with_where_clause(
        "tl",
        &format!(" WHERE chrom = '{}'", sql_escape_literal(source_chrom)),
    );
    // The legacy path accumulates every transformed SIFT batch before sorting.
    // The streaming path keeps the same key order by sorting the deduplicated
    // transcript rows upstream, then writes each transformed batch immediately.
    let stream_sift = std::env::var("VEPYR_TRANSLATION_STREAM_SORT").as_deref() == Ok("1");
    let query = if stream_sift {
        format!("{base_query} ORDER BY transcript_id")
    } else {
        base_query
    };
    let mut translation_options = options.clone();
    if let Ok(value) = std::env::var("VEPYR_TRANSLATION_PARTITIONS") {
        if let Ok(partitions) = value.parse::<usize>() {
            translation_options.partitions = partitions.max(1);
        }
    }
    eprintln!(
        "[TRANSLATION_CONFIG] chrom={} partitions={}",
        chrom, translation_options.partitions
    );
    let ctx = make_ctx_and_register(
        &translation_options,
        EnsemblEntityKind::Translation,
        "tl",
    )?;
    let plan = ctx.sql(&query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

    use crate::parquet_cache::scan::ContextParquetShardWriter;
    use crate::parquet_cache::write::point_lookup_writer_properties;
    let mut core_writer: Option<ContextParquetShardWriter> = None;
    let mut sift_batches: Vec<RecordBatch> = Vec::new();
    let mut sift_writer: Option<ArrowWriter<std::fs::File>> = None;
    let mut sift_rows = 0usize;

    while let Some(batch) = stream.next().await {
        let batch = drop_row_number_batch(batch?)?;

        if build_core {
            let target_schema = core_target_schema.as_ref().unwrap();
            let out = project_batch_to_schema(batch.clone(), Arc::clone(target_schema))?;
            let out = attach_schema_metadata_to_batch(out, &source_type, &options.cache_version)?;
            if out.num_rows() > 0 {
                if core_writer.is_none() {
                    let path = core_shard_path.as_ref().unwrap();
                    core_writer = Some(ContextParquetShardWriter::create(path, out.schema())?);
                }
                core_writer.as_mut().unwrap().write(&out)?;
            }
        }

        if build_sift {
            let map = uid_map.as_ref().unwrap();
            let out = transform_translation_sift_position_batch(
                batch,
                &source_type,
                &options.cache_version,
                map,
            )?;
            if out.num_rows() > 0 {
                if stream_sift {
                    if sift_writer.is_none() {
                        let path = sift_shard_path.as_ref().unwrap();
                        let props = point_lookup_writer_properties(&out.schema(), &["key"]);
                        let file = std::fs::File::create(path).map_err(|err| {
                            DataFusionError::Execution(format!(
                                "failed to create '{}': {err}", path.display()
                            ))
                        })?;
                        sift_writer = Some(ArrowWriter::try_new(file, out.schema(), Some(props)).map_err(|err| {
                            DataFusionError::Execution(format!("failed to open sift writer: {err}"))
                        })?);
                    }
                    sift_rows += out.num_rows();
                    sift_writer.as_mut().unwrap().write(&out).map_err(|err| {
                        DataFusionError::Execution(format!("failed to write sift batch: {err}"))
                    })?;
                } else {
                    sift_batches.push(out);
                }
            }
        }
    }

    let core_entry = if build_core {
        let rows = if let Some(writer) = core_writer {
            writer.finish()?
        } else {
            0
        };
        Some(ChromDatasetEntry::new(manifest_chrom.clone(), file_name.clone(), rows))
    } else {
        None
    };

    let sift_entry = if build_sift {
        let shard_path = sift_shard_path.as_ref().unwrap();
        let rows = if stream_sift {
            if let Some(writer) = sift_writer {
                writer.close().map_err(|err| {
                    DataFusionError::Execution(format!("failed to close sift writer: {err}"))
                })?;
            }
            sift_rows
        } else if sift_batches.is_empty() {
            0
        } else {
            let schema = sift_batches[0].schema();
            let combined = concat_batches(&schema, &sift_batches).map_err(|err| {
                DataFusionError::Execution(format!("failed to concat sift batches: {err}"))
            })?;
            let sorted = sort_record_batch_by_key(&combined)?;
            crate::parquet_cache::sift::write_sift_parquet(shard_path, &[sorted])?
        };
        Some(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
    } else {
        None
    };

    Ok((core_entry, sift_entry))
}

pub async fn build_parquet_translation_sift_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let (_, sift) = build_parquet_translation_chrom(options, chrom, false, true).await?;
    sift.ok_or_else(|| DataFusionError::Execution("translation_sift build failed".into()))
}

/// Remove an existing Parquet shard file, honoring `overwrite`. Shared by the
/// Parquet build drivers.
fn remove_existing_parquet_shard(path: &Path, overwrite: bool) -> Result<()> {
    if path.exists() {
        if !overwrite {
            return Err(DataFusionError::Execution(format!(
                "parquet shard '{}' already exists and overwrite=false",
                path.display()
            )));
        }
        std::fs::remove_file(path).map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to remove existing parquet shard '{}': {err}",
                path.display()
            ))
        })?;
    }
    Ok(())
}

/// The provider output schema for an entity (used to build the transcript export
/// query, which selects an explicit field list keyed on the provider schema).
fn provider_output_schema(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<SchemaRef> {
    use datafusion::catalog::TableProvider;
    let mut provider_options = EnsemblCacheOptions::new(&options.cache_root)
        .with_cache_source_type(options.cache_source_type)
        .with_expected_cache_version(&options.cache_version);
    provider_options.target_partitions = Some(options.partitions);
    Ok(EnsemblCacheTableProvider::for_entity(kind, provider_options)?.schema())
}

/// Stream a source query to a dictionary-enabled context Parquet shard, applying
/// `transform` to each batch. Preserves source (query) order — no explicit sort —
/// so the shard matches the Parquet dataset row order. Returns the row count
/// (0 → no file written).
async fn write_query_stream_to_parquet<F>(
    ctx: &SessionContext,
    query: &str,
    path: &Path,
    transform: F,
) -> Result<usize>
where
    F: Fn(RecordBatch) -> Result<RecordBatch>,
{
    use crate::parquet_cache::scan::ContextParquetShardWriter;
    let plan = ctx.sql(query).await?.create_physical_plan().await?;
    let mut stream = datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;
    let mut writer: Option<ContextParquetShardWriter> = None;
    let mut rows = 0usize;
    while let Some(batch) = stream.next().await {
        let out = transform(batch?)?;
        if out.num_rows() == 0 {
            continue;
        }
        if writer.is_none() {
            writer = Some(ContextParquetShardWriter::create(path, out.schema())?);
        }
        writer.as_mut().unwrap().write(&out)?;
        rows += out.num_rows();
    }
    if let Some(writer) = writer {
        writer.finish()?;
    }
    Ok(rows)
}

/// Build, for one chromosome, one
/// scan entity (`transcript`/`exon`/`regulatory`/`motif`): dictionary-enabled
/// `parquet.<entity>/<chrom>.parquet`, written verbatim from the same export
/// query + per-batch transform the Parquet path uses (transcript also attaches the
/// shared `transcript_uid`). Read back by a full projected scan.
pub async fn build_parquet_context_entity_chrom(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let entity = entity_output_name(kind);
    let entity_dir = PathBuf::from(&options.output_dir).join(entity);
    std::fs::create_dir_all(&entity_dir).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to create '{}': {err}",
            entity_dir.display()
        ))
    })?;

    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    let table_name = entity_table_name(kind);
    let manifest_chrom = canonical_chrom_label(chrom);
    let file_name = format!("{manifest_chrom}.parquet");
    let shard_path = entity_dir.join(&file_name);
    remove_existing_parquet_shard(&shard_path, options.overwrite)?;

    // Only transcript's export query needs the provider schema (explicit select list).
    let provider_schema = if kind == EnsemblEntityKind::Transcript {
        Some(provider_output_schema(options, kind)?)
    } else {
        None
    };
    let mut query = build_export_query(
        kind,
        table_name,
        Some(source_chrom),
        provider_schema.as_deref(),
    );
    // Optional external genomic-range partition. Each synchronous worker can
    // process a disjoint source range without materialising the full large
    // chromosome. Duplicate transcript/exon rows share the same `start`, so
    // filtering before the per-shard dedup keeps them in one range.
    if let (Ok(start), Ok(end)) = (
        std::env::var("VEPYR_CONTEXT_START").ok().and_then(|v| v.parse::<i64>().ok()).ok_or(()),
        std::env::var("VEPYR_CONTEXT_END").ok().and_then(|v| v.parse::<i64>().ok()).ok_or(()),
    ) {
        let needle = format!(" WHERE chrom = '{}'", sql_escape_literal(source_chrom));
        let replacement = format!(
            " WHERE chrom = '{}' AND start >= {} AND start < {}",
            sql_escape_literal(source_chrom), start, end
        );
        query = query.replacen(&needle, &replacement, 1);
        eprintln!("[CONTEXT_RANGE] chrom={} start={} end={}", chrom, start, end);
    }
    let source_type = options.cache_source_type.as_str().to_string();
    let cache_version = options.cache_version.clone();
    let ctx = make_ctx_and_register(options, kind, table_name)?;

    let rows = if kind == EnsemblEntityKind::Transcript {
        let uid_map = Arc::new(load_transcript_uid_map(options, source_chrom).await?);
        write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
            attach_transcript_uid_batch(batch, &source_type, &cache_version, &uid_map)
        })
        .await?
    } else {
        write_query_stream_to_parquet(&ctx, &query, &shard_path, move |batch| {
            attach_schema_metadata_to_batch(batch, &source_type, &cache_version)
        })
        .await?
    };
    Ok(ChromDatasetEntry::new(manifest_chrom, file_name, rows))
}

/// Build `translation_core`,
/// one chromosome: dictionary-enabled `translation_core/<chrom>.parquet`,
/// projected to `translation_core_schema` (big sequences + `list<struct>`
/// `protein_features`). Read back by a full projected scan.
pub async fn build_parquet_translation_core_chrom(
    options: &CacheBuildOptions,
    chrom: &str,
) -> Result<ChromDatasetEntry> {
    let (core, _) = build_parquet_translation_chrom(options, chrom, true, false).await?;
    core.ok_or_else(|| DataFusionError::Execution("translation_core build failed".into()))
}

/// Output directory for a Parquet entity's shard set: `parquet.<entity>/`.
fn parquet_entity_dir(options: &CacheBuildOptions, entity: &str) -> Result<PathBuf> {
    let path = Path::new(&options.output_dir).join(entity);
    std::fs::create_dir_all(&path).map_err(|err| {
        DataFusionError::Execution(format!("failed to create {}: {err}", path.display()))
    })?;
    Ok(path)
}

/// Skip a Parquet entity build when a manifest and at least one `.parquet` shard
/// already exist and `overwrite` is false. Parquet analogue of
fn should_skip_parquet_entity(entity_dir: &Path, overwrite: bool) -> bool {
    !overwrite
        && entity_dir.join(CHROM_MANIFEST_FILE).exists()
        && std::fs::read_dir(entity_dir)
            .ok()
            .into_iter()
            .flat_map(|entries| entries.flatten())
            .any(|entry| entry.file_name().to_string_lossy().ends_with(".parquet"))
}

/// Build one entity as Parquet shards under `parquet.<entity>/`, writing the
/// per-chrom `chrom_manifest.json`.
/// this is the entry point the [`crate::cache_builder::CacheBuilder`] dispatches
/// to.
pub async fn build_parquet_entity(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<Vec<EntityStats>> {
    match kind {
        EnsemblEntityKind::Variation => build_parquet_variation(options).await,
        EnsemblEntityKind::Translation => build_parquet_translation(options).await,
        EnsemblEntityKind::Transcript
        | EnsemblEntityKind::Exon
        | EnsemblEntityKind::RegulatoryFeature
        | EnsemblEntityKind::MotifFeature => build_parquet_context(options, kind).await,
    }
}

async fn build_parquet_variation(options: &CacheBuildOptions) -> Result<Vec<EntityStats>> {
    let entity_dir = parquet_entity_dir(options, "variation")?;
    if should_skip_parquet_entity(&entity_dir, options.overwrite) {
        info!("variation already exists, skipping (use overwrite to rebuild)");
        return Ok(vec![EntityStats {
            entity: "variation".to_string(),
            parquet_files: vec![],
        }]);
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Variation, "var").await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let entry = build_parquet_variation_chrom(options, &chrom).await?;
        if entry.rows > 0 {
            let shard_path = entity_dir.join(&entry.dataset);
            files.push((shard_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: "variation".to_string(),
        parquet_files: files,
    }])
}

async fn build_parquet_context(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
) -> Result<Vec<EntityStats>> {
    let entity = entity_output_name(kind);
    let entity_dir = parquet_entity_dir(options, entity)?;
    if should_skip_parquet_entity(&entity_dir, options.overwrite) {
        info!("{entity} already exists, skipping (use overwrite to rebuild)");
        return Ok(vec![EntityStats {
            entity: entity.to_string(),
            parquet_files: vec![],
        }]);
    }

    let table_name = entity_table_name(kind);
    let chroms = discover_chroms(options, kind, table_name).await?;
    let mut manifest_entries = Vec::new();
    let mut files = Vec::new();
    for chrom in chroms {
        let entry = build_parquet_context_entity_chrom(options, kind, &chrom).await?;
        if entry.rows > 0 {
            let shard_path = entity_dir.join(&entry.dataset);
            files.push((shard_path.to_string_lossy().to_string(), entry.rows));
            manifest_entries.push(entry);
        }
    }

    ChromManifest::new(manifest_entries).merge_write_to_entity_dir(&entity_dir)?;
    Ok(vec![EntityStats {
        entity: entity.to_string(),
        parquet_files: files,
    }])
}

async fn build_parquet_translation(options: &CacheBuildOptions) -> Result<Vec<EntityStats>> {
    let core_dir = parquet_entity_dir(options, "translation_core")?;
    let sift_dir = parquet_entity_dir(options, "translation_sift")?;
    let core_done = should_skip_parquet_entity(&core_dir, options.overwrite);
    let sift_done = should_skip_parquet_entity(&sift_dir, options.overwrite);

    if core_done && sift_done {
        info!("translation_core and translation_sift already exist, skipping (use overwrite to rebuild)");
        return Ok(vec![
            EntityStats {
                entity: "translation_core".to_string(),
                parquet_files: vec![],
            },
            EntityStats {
                entity: "translation_sift".to_string(),
                parquet_files: vec![],
            },
        ]);
    }

    let chroms = discover_chroms(options, EnsemblEntityKind::Translation, "tl").await?;
    let mut core_manifest = Vec::new();
    let mut sift_manifest = Vec::new();
    let mut core_files = Vec::new();
    let mut sift_files = Vec::new();

    for chrom in chroms {
        let (core_opt, sift_opt) = build_parquet_translation_chrom(
            options,
            &chrom,
            !core_done,
            !sift_done,
        ).await?;
        if let Some(entry) = core_opt {
            if entry.rows > 0 {
                let shard_path = core_dir.join(&entry.dataset);
                core_files.push((shard_path.to_string_lossy().to_string(), entry.rows));
                core_manifest.push(entry);
            }
        }
        if let Some(entry) = sift_opt {
            if entry.rows > 0 {
                let shard_path = sift_dir.join(&entry.dataset);
                sift_files.push((shard_path.to_string_lossy().to_string(), entry.rows));
                sift_manifest.push(entry);
            }
        }
    }

    if !core_done {
        ChromManifest::new(core_manifest).merge_write_to_entity_dir(&core_dir)?;
    }
    if !sift_done {
        ChromManifest::new(sift_manifest).merge_write_to_entity_dir(&sift_dir)?;
    }

    let mut results = Vec::new();
    results.push(EntityStats {
        entity: "translation_core".to_string(),
        parquet_files: if !core_done { core_files } else { vec![] },
    });
    results.push(EntityStats {
        entity: "translation_sift".to_string(),
        parquet_files: if !sift_done { sift_files } else { vec![] },
    });
    Ok(results)
}

// ---------------------------------------------------------------------------
// Concurrent build: LPT-ordered, weighted-semaphore-bounded scheduling across
// every (entity, chrom) build primitive. `build_all()` stays the sequential
// default (`build_concurrency == 1`); `CacheBuilder::with_build_concurrency`
// routes here for `> 1`.
// ---------------------------------------------------------------------------

/// Approximate GRCh38 chromosome lengths (bp), used only to order units
/// longest-first (LPT) before scheduling. Exact per-unit cost isn't known
/// upfront; chromosome length is a cheap, good-enough proxy for it.
fn approx_chrom_length_bp(chrom: &str) -> u64 {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    match bare {
        "1" => 248_956_422,
        "2" => 242_193_529,
        "3" => 198_295_559,
        "4" => 190_214_555,
        "5" => 181_538_259,
        "6" => 170_805_979,
        "7" => 159_345_973,
        "8" => 145_138_636,
        "9" => 138_394_717,
        "10" => 133_797_422,
        "11" => 135_086_622,
        "12" => 133_275_309,
        "13" => 114_364_328,
        "14" => 107_043_718,
        "15" => 101_991_189,
        "16" => 90_338_345,
        "17" => 83_257_441,
        "18" => 80_373_285,
        "19" => 58_617_616,
        "20" => 64_444_167,
        "21" => 46_709_983,
        "22" => 50_818_468,
        "X" => 156_040_895,
        "Y" => 57_227_415,
        "MT" | "M" => 16_569,
        _ => 100_000_000, // unknown contig: assume mid-size rather than starve it
    }
}

/// One independently-schedulable cache-build task.
#[derive(Debug)]
enum BuildUnit {
    Variation(String),
    Context(EnsemblEntityKind, String),
    Translation {
        chrom: String,
        build_core: bool,
        build_sift: bool,
    },
}

enum UnitRunResult {
    Completed(&'static str, ChromDatasetEntry),
    TranslationCompleted(Option<ChromDatasetEntry>, Option<ChromDatasetEntry>),
}

impl BuildUnit {
    fn cost(&self) -> u64 {
        match self {
            BuildUnit::Variation(chrom) => {
                let chrom_len = approx_chrom_length_bp(chrom);
                chrom_len.saturating_mul(100)
            }
            BuildUnit::Context(kind, chrom) => {
                let chrom_len = approx_chrom_length_bp(chrom);
                let weight = match kind {
                    EnsemblEntityKind::Transcript => 10,
                    EnsemblEntityKind::Exon => 5,
                    EnsemblEntityKind::RegulatoryFeature => 2,
                    EnsemblEntityKind::MotifFeature => 1,
                    _ => 1,
                };
                chrom_len.saturating_mul(weight)
            }
            BuildUnit::Translation { chrom, build_core, build_sift } => {
                let mut factor = 0u64;
                if *build_core { factor += 20; }
                if *build_sift { factor += 30; }
                approx_chrom_length_bp(chrom).saturating_mul(factor)
            }
        }
    }

    /// Semaphore permits this unit consumes.
    fn weight(&self, var_weight: u32) -> u32 {
        match self {
            BuildUnit::Variation(_) => var_weight,
            _ => 1,
        }
    }

    /// Conservative resident-memory estimate used by the adaptive scheduler.
    /// Variation is the dominant unit; estimates are deliberately capped so
    /// larger-memory hosts can run multiple chromosomes while small hosts keep
    /// one heavy unit inside the OS budget.
    fn memory_weight_tokens(&self, capacity: u32) -> u32 {
        let scale = |chrom: &str| {
            (approx_chrom_length_bp(chrom) as f64 / approx_chrom_length_bp("1") as f64)
                .clamp(0.20, 1.0)
        };
        let estimate = match self {
            BuildUnit::Variation(chrom) => (96.0 * scale(chrom)).ceil() as u32, // chr1: 24 GiB
            BuildUnit::Translation { chrom, .. } => (48.0 * scale(chrom)).ceil() as u32,
            BuildUnit::Context(_, chrom) => (8.0 * scale(chrom)).ceil() as u32,
        };
        estimate.min(capacity).max(1)
    }

    async fn run(
        self,
        options: Arc<CacheBuildOptions>,
    ) -> Result<UnitRunResult> {
        match self {
            BuildUnit::Variation(chrom) => {
                let entry = build_parquet_variation_chrom(&options, &chrom).await?;
                Ok(UnitRunResult::Completed("variation", entry))
            }
            BuildUnit::Context(kind, chrom) => {
                let entity = entity_output_name(kind);
                let entry = build_parquet_context_entity_chrom(&options, kind, &chrom).await?;
                Ok(UnitRunResult::Completed(entity, entry))
            }
            BuildUnit::Translation { chrom, build_core, build_sift } => {
                let (core, sift) = build_parquet_translation_chrom(&options, &chrom, build_core, build_sift).await?;
                Ok(UnitRunResult::TranslationCompleted(core, sift))
            }
        }
    }
}

fn host_memory_bytes() -> Option<u64> {
    #[cfg(target_os = "macos")]
    {
        let output = std::process::Command::new("sysctl")
            .args(["-n", "hw.memsize"])
            .output()
            .ok()?;
        return String::from_utf8_lossy(&output.stdout).trim().parse().ok();
    }
    #[cfg(target_os = "linux")]
    {
        let text = std::fs::read_to_string("/proc/meminfo").ok()?;
        let kb = text.lines().find_map(|line| {
            line.strip_prefix("MemTotal:")?.split_whitespace().next()?.parse::<u64>().ok()
        })?;
        return Some(kb * 1024);
    }
    #[allow(unreachable_code)]
    None
}

fn adaptive_memory_capacity_tokens() -> u32 {
    let Some(total) = host_memory_bytes() else { return 96 };
    let budget = total.saturating_mul(MEMORY_BUDGET_FRACTION_NUM) / MEMORY_BUDGET_FRACTION_DEN;
    (budget / MEMORY_TOKEN_BYTES).max(1).min(u32::MAX as u64) as u32
}

async fn wait_for_memory_headroom() {
    let minimum_free = std::env::var("VEPYR_MIN_FREE_MEMORY_PERCENT")
        .ok().and_then(|v| v.parse::<u8>().ok()).unwrap_or(20);
    loop {
        let free = std::process::Command::new("memory_pressure")
            .args(["-Q"])
            .output()
            .ok()
            .and_then(|o| String::from_utf8(o.stdout).ok())
            .and_then(|s| s.lines().find_map(|line| {
                let (key, value) = line.split_once(':')?;
                (key.contains("free percentage"))
                    .then(|| value.trim().trim_end_matches('%').parse::<u8>().ok())
                    .flatten()
            }));
        if free.map(|v| v >= minimum_free).unwrap_or(true) { return; }
        tokio::time::sleep(std::time::Duration::from_millis(500)).await;
    }
}







/// Read `VEPYR_VAR_WEIGHT` (default 2): extra semaphore permits a variation
/// unit consumes, so a memory-constrained host can bound peak RSS by trading
/// away concurrency for the dominant entity specifically instead of globally.
fn variation_permit_weight() -> u32 {
    std::env::var("VEPYR_VAR_WEIGHT")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|w| *w >= 1)
        .unwrap_or(2)
}

/// Build all entities with cross-unit concurrency: every (entity, chromosome)
/// build primitive becomes one [`BuildUnit`], sorted longest-first (LPT —
/// Longest Processing Time first) so the tail (dominated by the `variation`
/// entity) starts as early as possible, then run through a `JoinSet` bounded
/// by a weighted `Semaphore`. `build_concurrency` caps how many units run at
/// once; `options.partitions` is unrelated — it only controls DataFusion's
/// per-query internal parallelism, not this cross-unit scheduling.
///
/// Each entity's manifest is written once, after all of that entity's units
/// complete, mirroring the sequential path's per-entity
/// `merge_write_to_entity_dir` call.
pub async fn build_all_units(
    options: &CacheBuildOptions,
    build_concurrency: usize,
) -> Result<Vec<EntityStats>> {
    use tokio::sync::Semaphore;
    use tokio::task::JoinSet;

    let build_concurrency = build_concurrency.max(1);
    let options = Arc::new(options.clone());
    let _var_weight = variation_permit_weight();

    let entity_kinds = [
        EnsemblEntityKind::Variation,
        EnsemblEntityKind::Transcript,
        EnsemblEntityKind::Exon,
        EnsemblEntityKind::Translation,
        EnsemblEntityKind::RegulatoryFeature,
        EnsemblEntityKind::MotifFeature,
    ];

    // Discovery: cheap, sequential (schema introspection only). Entities
    // already fully built (overwrite=false) short-circuit to a placeholder
    // EntityStats exactly like `build_all()`'s per-entity dispatch does.
    let mut done: Vec<EntityStats> = Vec::new();
    let mut units: Vec<BuildUnit> = Vec::new();
    let mut entity_dirs: HashMap<&'static str, PathBuf> = HashMap::new();

    for kind in entity_kinds {
        match kind {
            EnsemblEntityKind::Variation => {
                let entity_dir = parquet_entity_dir(&options, "variation")?;
                if should_skip_parquet_entity(&entity_dir, options.overwrite) {
                    info!("variation already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: "variation".to_string(),
                        parquet_files: vec![],
                    });
                    continue;
                }
                entity_dirs.insert("variation", entity_dir);
                let chroms =
                    discover_chroms(&options, EnsemblEntityKind::Variation, "var").await?;
                units.extend(chroms.into_iter().map(BuildUnit::Variation));
            }
            EnsemblEntityKind::Translation => {
                let core_dir = parquet_entity_dir(&options, "translation_core")?;
                let sift_dir = parquet_entity_dir(&options, "translation_sift")?;
                let core_done = should_skip_parquet_entity(&core_dir, options.overwrite);
                let sift_done = should_skip_parquet_entity(&sift_dir, options.overwrite);
                if core_done {
                    info!("translation_core already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: "translation_core".to_string(),
                        parquet_files: vec![],
                    });
                } else {
                    entity_dirs.insert("translation_core", core_dir);
                }
                if sift_done {
                    info!("translation_sift already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: "translation_sift".to_string(),
                        parquet_files: vec![],
                    });
                } else {
                    entity_dirs.insert("translation_sift", sift_dir);
                }
                if !core_done || !sift_done {
                    let chroms =
                        discover_chroms(&options, EnsemblEntityKind::Translation, "tl").await?;
                    units.extend(chroms.into_iter().map(|chrom| BuildUnit::Translation {
                        chrom,
                        build_core: !core_done,
                        build_sift: !sift_done,
                    }));
                }
            }
            other => {
                let entity = entity_output_name(other);
                let entity_dir = parquet_entity_dir(&options, entity)?;
                if should_skip_parquet_entity(&entity_dir, options.overwrite) {
                    info!("{entity} already exists, skipping (use overwrite to rebuild)");
                    done.push(EntityStats {
                        entity: entity.to_string(),
                        parquet_files: vec![],
                    });
                    continue;
                }
                entity_dirs.insert(entity, entity_dir);
                let table_name = entity_table_name(other);
                let chroms = discover_chroms(&options, other, table_name).await?;
                units.extend(
                    chroms
                        .into_iter()
                        .map(move |chrom| BuildUnit::Context(other, chrom)),
                );
            }
        }
    }

    // LPT: longest first, so the tail (dominated by `variation`) is in
    // flight as early as possible instead of starting last.
    units.sort_by_key(|b| std::cmp::Reverse(b.cost()));

    let memory_capacity = adaptive_memory_capacity_tokens();
    eprintln!(
        "[MEMORY_GATE] capacity_tokens={} budget_gib={:.1}",
        memory_capacity,
        memory_capacity as f64 * MEMORY_TOKEN_BYTES as f64 / 1024f64.powi(3)
    );
    let var_concurrency = std::env::var("VEPYR_VAR_CONCURRENCY")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&c| c >= 1)
        .unwrap_or(build_concurrency);
    let memory_semaphore = Arc::new(Semaphore::new(memory_capacity as usize));
    let var_semaphore = Arc::new(Semaphore::new(var_concurrency));
    let context_semaphore = Arc::new(Semaphore::new(build_concurrency));
    let mut join_set: JoinSet<Result<UnitRunResult>> = JoinSet::new();
    for unit in units {
        let is_var = matches!(unit, BuildUnit::Variation(_));
        let semaphore = if is_var {
            Arc::clone(&var_semaphore)
        } else {
            Arc::clone(&context_semaphore)
        };
        let memory_semaphore = Arc::clone(&memory_semaphore);
        let memory_weight = unit.memory_weight_tokens(memory_capacity);
        let options = Arc::clone(&options);
        let unit_desc = format!("{:?}", unit);
        join_set.spawn(async move {
            wait_for_memory_headroom().await;
            let _memory_permits = memory_semaphore
                .acquire_many_owned(memory_weight)
                .await
                .map_err(|err| DataFusionError::Execution(format!("memory semaphore closed: {err}")))?;
            let _permit = semaphore
                .acquire_owned()
                .await
                .map_err(|err| DataFusionError::Execution(format!("semaphore closed: {err}")))?;
            let t0 = std::time::Instant::now();
            let res = unit.run(options).await;
            eprintln!("[TIMING] Unit {} finished in {:.2}s", unit_desc, t0.elapsed().as_secs_f64());
            res
        });
    }

    let mut manifest_entries: HashMap<&'static str, Vec<ChromDatasetEntry>> = HashMap::new();
    let mut files: HashMap<&'static str, Vec<(String, usize)>> = HashMap::new();

    while let Some(joined) = join_set.join_next().await {
        let result = joined
            .map_err(|err| DataFusionError::Execution(format!("build task panicked: {err}")))??;
        match result {
            UnitRunResult::Completed(entity, entry) => {
                if entry.rows > 0 {
                    let entity_dir = entity_dirs
                        .get(entity)
                        .expect("entity_dirs populated for every scheduled unit's entity");
                    let shard_path = entity_dir.join(&entry.dataset);
                    files
                        .entry(entity)
                        .or_default()
                        .push((shard_path.to_string_lossy().to_string(), entry.rows));
                    manifest_entries.entry(entity).or_default().push(entry);
                }
            }
            UnitRunResult::TranslationCompleted(core_opt, sift_opt) => {
                if let Some(entry) = core_opt {
                    if entry.rows > 0 {
                        let entity_dir = entity_dirs
                            .get("translation_core")
                            .expect("entity_dirs populated for translation_core");
                        let shard_path = entity_dir.join(&entry.dataset);
                        files
                            .entry("translation_core")
                            .or_default()
                            .push((shard_path.to_string_lossy().to_string(), entry.rows));
                        manifest_entries.entry("translation_core").or_default().push(entry);
                    }
                }
                if let Some(entry) = sift_opt {
                    if entry.rows > 0 {
                        let entity_dir = entity_dirs
                            .get("translation_sift")
                            .expect("entity_dirs populated for translation_sift");
                        let shard_path = entity_dir.join(&entry.dataset);
                        files
                            .entry("translation_sift")
                            .or_default()
                            .push((shard_path.to_string_lossy().to_string(), entry.rows));
                        manifest_entries.entry("translation_sift").or_default().push(entry);
                    }
                }
            }
        }
    }

    for (entity, entity_dir) in &entity_dirs {
        let entries = manifest_entries.remove(entity).unwrap_or_default();
        ChromManifest::new(entries).merge_write_to_entity_dir(entity_dir)?;
        done.push(EntityStats {
            entity: entity.to_string(),
            parquet_files: files.remove(*entity).unwrap_or_default(),
        });
    }

    Ok(done)
}

/// Sort a translation_sift batch ascending by its `key` column.
fn sort_record_batch_by_key(batch: &RecordBatch) -> Result<RecordBatch> {
    let key_idx = batch.schema().index_of("key")?;
    let indices = sort_to_indices(batch.column(key_idx), None, None)
        .map_err(|err| DataFusionError::Execution(format!("failed to sort sift key: {err}")))?;
    let columns = batch
        .columns()
        .iter()
        .map(|c| take(c, &indices, None))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|err| DataFusionError::Execution(format!("failed to take sorted sift: {err}")))?;
    RecordBatch::try_new(batch.schema(), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to rebuild sorted sift batch: {err}"))
    })
}

async fn discover_chroms(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<Vec<String>> {
    discover_chroms_and_schema(options, kind, table_name)
        .await
        .map(|(chroms, _)| chroms)
}

async fn discover_chroms_and_schema(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<(Vec<String>, SchemaRef)> {
    let ctx = make_ctx_and_register(options, kind, table_name)?;
    let table = ctx.table(table_name).await?;
    let schema = Arc::clone(table.schema().inner());
    let chroms = filter_chroms(chroms_from_schema(&schema), options.chrom_filter.as_deref());
    Ok((chroms, schema))
}

/// Restrict discovered chromosomes to the build's allow-list, comparing by the
/// bare label (any `chr` prefix stripped) so `"1"` and `"chr1"` match. `None`
/// passes every chromosome through unchanged.
fn filter_chroms(chroms: Vec<String>, allow: Option<&[String]>) -> Vec<String> {
    let Some(allow) = allow else {
        return chroms;
    };
    let bare = |label: &str| label.strip_prefix("chr").unwrap_or(label).to_string();
    let allowed: std::collections::HashSet<String> = allow.iter().map(|c| bare(c)).collect();
    chroms
        .into_iter()
        .filter(|chrom| allowed.contains(&bare(chrom)))
        .collect()
}

fn make_ctx_and_register_partitions(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
    partitions: usize,
) -> Result<SessionContext> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let spill_translation = kind == EnsemblEntityKind::Translation
        && std::env::var("VEPYR_TRANSLATION_SPILL").as_deref() == Ok("1");
    let ctx = if spill_translation {
        let gib = std::env::var("VEPYR_TRANSLATION_SPILL_GIB")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(6)
            .max(1);
        let runtime = RuntimeEnvBuilder::new()
            .with_memory_pool(Arc::new(FairSpillPool::new(
                (gib * 1024 * 1024 * 1024) as usize,
            )))
            .build_arc()
            .map_err(|err| DataFusionError::Execution(format!("translation spill runtime: {err}")))?;
        SessionContext::new_with_config_rt(config, runtime)
    } else {
        SessionContext::new_with_config(config)
    };
    let mut provider_options = EnsemblCacheOptions::new(&options.cache_root)
        .with_cache_source_type(options.cache_source_type)
        .with_expected_cache_version(&options.cache_version);
    provider_options.target_partitions = Some(partitions);
    let provider = EnsemblCacheTableProvider::for_entity(kind, provider_options)?;
    ctx.register_table(table_name, provider)?;
    Ok(ctx)
}

fn make_ctx_and_register(
    options: &CacheBuildOptions,
    kind: EnsemblEntityKind,
    table_name: &str,
) -> Result<SessionContext> {
    let partitions = match kind {
        EnsemblEntityKind::Variation => options.partitions.max(4),
        EnsemblEntityKind::Translation | EnsemblEntityKind::Transcript => options.partitions.max(4),
        _ => options.partitions.max(2),
    };
    make_ctx_and_register_partitions(options, kind, table_name, partitions)
}

/// Project a source variation batch to the Parquet output schema, derive the
/// `tier` column from `warm_starts`, and keep only rows whose tier matches
/// `keep_tier`. The returned batch may be empty (all rows belonged to the other
/// tier); the streaming writer drops empty batches.
fn transform_variation_tier_batch(
    batch: RecordBatch,
    source_type: &str,
    warm_starts: &BTreeSet<i64>,
    keep_tier: i8,
) -> Result<RecordBatch> {
    let schema = batch.schema();
    let target_schema = Arc::new(variation_projected_schema(schema.as_ref(), source_type)?);
    let mut columns = Vec::new();
    let forbidden = VARIATION_FORBIDDEN_COLUMNS
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    for name in VARIATION_REQUIRED_COLUMNS {
        if forbidden.contains(name) {
            continue;
        }
        let (index, _) = schema.column_with_name(name).ok_or_else(|| {
            DataFusionError::Execution(format!("variation source batch missing column {name}"))
        })?;
        if *name == "start" || *name == "end" {
            columns.push(
                cast(batch.column(index).as_ref(), &DataType::UInt32).map_err(|err| {
                    DataFusionError::Execution(format!(
                        "failed to cast variation {name} to UInt32: {err}"
                    ))
                })?,
            );
        } else {
            columns.push(batch.column(index).clone());
        }
        if *name == "clin_sig_allele" {
            for optional_name in VARIATION_OPTIONAL_COLUMNS {
                if let Some((optional_index, _)) = schema.column_with_name(optional_name) {
                    columns.push(batch.column(optional_index).clone());
                }
            }
        }
    }

    // Derive tier from the source `start` and build the keep mask in one pass.
    let (start_index, _) = schema.column_with_name("start").ok_or_else(|| {
        DataFusionError::Execution("variation source batch missing column start".to_string())
    })?;
    let start = cast(batch.column(start_index).as_ref(), &DataType::Int64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation start to Int64: {err}"))
    })?;
    let start = start
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| DataFusionError::Execution("variation start must be Int64".to_string()))?;

    let mut tier_values = Vec::with_capacity(batch.num_rows());
    let mut keep = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        let tier = if !start.is_null(row) && warm_starts.contains(&start.value(row)) {
            0i8
        } else {
            1i8
        };
        tier_values.push(tier);
        keep.push(tier == keep_tier);
    }
    columns.push(Arc::new(Int8Array::from(tier_values)) as ArrayRef);

    let full = RecordBatch::try_new(target_schema, columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to build Lance variation batch: {err}"))
    })?;
    let mask = BooleanArray::from(keep);
    filter_record_batch(&full, &mask).map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to filter Lance variation tier batch: {err}"
        ))
    })
}

fn transform_variation_split_batch(
    batch: RecordBatch,
    source_type: &str,
    warm_starts: &HashSet<i64>,
) -> Result<(Option<RecordBatch>, Option<RecordBatch>)> {
    let schema = batch.schema();
    let target_schema = Arc::new(variation_projected_schema(schema.as_ref(), source_type)?);
    let mut columns = Vec::new();
    let forbidden = VARIATION_FORBIDDEN_COLUMNS
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    for name in VARIATION_REQUIRED_COLUMNS {
        if forbidden.contains(name) {
            continue;
        }
        let (index, _) = schema.column_with_name(name).ok_or_else(|| {
            DataFusionError::Execution(format!("variation source batch missing column {name}"))
        })?;
        if *name == "start" || *name == "end" {
            columns.push(
                cast(batch.column(index).as_ref(), &DataType::UInt32).map_err(|err| {
                    DataFusionError::Execution(format!(
                        "failed to cast variation {name} to UInt32: {err}"
                    ))
                })?,
            );
        } else {
            columns.push(batch.column(index).clone());
        }
        if *name == "clin_sig_allele" {
            for optional_name in VARIATION_OPTIONAL_COLUMNS {
                if let Some((optional_index, _)) = schema.column_with_name(optional_name) {
                    columns.push(batch.column(optional_index).clone());
                }
            }
        }
    }

    let (start_index, _) = schema.column_with_name("start").ok_or_else(|| {
        DataFusionError::Execution("variation source batch missing column start".to_string())
    })?;
    let start = cast(batch.column(start_index).as_ref(), &DataType::Int64).map_err(|err| {
        DataFusionError::Execution(format!("failed to cast variation start to Int64: {err}"))
    })?;
    let start = start
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| DataFusionError::Execution("variation start must be Int64".to_string()))?;

    let mut tier_values = Vec::with_capacity(batch.num_rows());
    let mut warm_mask = Vec::with_capacity(batch.num_rows());
    let mut cold_mask = Vec::with_capacity(batch.num_rows());
    let mut warm_count = 0;
    let mut cold_count = 0;
    for row in 0..batch.num_rows() {
        let is_warm = !start.is_null(row) && warm_starts.contains(&start.value(row));
        if is_warm {
            tier_values.push(0i8);
            warm_mask.push(true);
            cold_mask.push(false);
            warm_count += 1;
        } else {
            tier_values.push(1i8);
            warm_mask.push(false);
            cold_mask.push(true);
            cold_count += 1;
        }
    }
    columns.push(Arc::new(Int8Array::from(tier_values)) as ArrayRef);

    let full = RecordBatch::try_new(target_schema, columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to build Lance variation batch: {err}"))
    })?;

    let warm_batch = if warm_count > 0 {
        if warm_count == batch.num_rows() {
            Some(full.clone())
        } else {
            let mask = BooleanArray::from(warm_mask);
            Some(filter_record_batch(&full, &mask).map_err(|err| {
                DataFusionError::Execution(format!("failed to filter variation warm batch: {err}"))
            })?)
        }
    } else {
        None
    };

    let cold_batch = if cold_count > 0 {
        if cold_count == batch.num_rows() {
            Some(full)
        } else {
            let mask = BooleanArray::from(cold_mask);
            Some(filter_record_batch(&full, &mask).map_err(|err| {
                DataFusionError::Execution(format!("failed to filter variation cold batch: {err}"))
            })?)
        }
    } else {
        None
    };

    Ok((warm_batch, cold_batch))
}


fn attach_schema_metadata_to_batch(
    batch: RecordBatch,
    source_type: &str,
    cache_version: &str,
) -> Result<RecordBatch> {
    let schema = with_cache_identity_metadata(batch.schema().as_ref(), source_type, cache_version);
    RecordBatch::try_new(Arc::new(schema), batch.columns().to_vec()).map_err(|err| {
        DataFusionError::Execution(format!("failed to attach Lance schema metadata: {err}"))
    })
}

/// Assign a dense `transcript_uid` (0..N) to each id, by the order given. Callers
/// pass a sorted-unique id list, so the uid is a stable function of the id set —
/// the transcript build and the sift build derive identical uids from the same
/// transcript source without depending on each other's artifacts.
fn assign_transcript_uids(sorted_unique_ids: &[String]) -> HashMap<String, u32> {
    sorted_unique_ids
        .iter()
        .enumerate()
        .map(|(index, id)| (id.clone(), index as u32))
        .collect()
}

/// Build the `transcript_id -> transcript_uid` map for one chromosome by ranking
static TRANSCRIPT_UID_CACHE: std::sync::LazyLock<
    std::sync::Mutex<HashMap<String, HashMap<String, u32>>>,
> = std::sync::LazyLock::new(|| std::sync::Mutex::new(HashMap::new()));

/// Load the chromosome's `transcript_uid` map, an ascending dense `0..N` rank of
/// the transcript source's sorted-unique `stable_id`s. Shared by the transcript
/// build (writes the `transcript_uid` column) and the sift build (keys rows).
async fn load_transcript_uid_map(
    options: &CacheBuildOptions,
    source_chrom: &str,
) -> Result<HashMap<String, u32>> {
    let cache_key = format!("{}:{}:{}", options.cache_root, options.cache_version, source_chrom);
    if let Ok(guard) = TRANSCRIPT_UID_CACHE.lock() {
        if let Some(cached) = guard.get(&cache_key) {
            return Ok(cached.clone());
        }
    }

    let ctx = make_ctx_and_register(options, EnsemblEntityKind::Transcript, "tx")?;
    let sql = format!(
        "SELECT DISTINCT stable_id FROM tx WHERE chrom = '{}' ORDER BY stable_id",
        sql_escape_literal(source_chrom)
    );
    let df = ctx.sql(&sql).await?;
    let mut stream = df.execute_stream().await?;
    let mut ids = Vec::new();
    while let Some(batch) = stream.next().await {
        let batch = batch?;
        let idx = batch.schema().index_of("stable_id")?;
        let column = batch.column(idx);
        for row in 0..batch.num_rows() {
            if let Some(id) = string_at(column.as_ref(), row) {
                ids.push(id);
            }
        }
    }
    // DISTINCT + ORDER BY already yields sorted-unique ids; normalize defensively
    // so uid assignment is deterministic regardless of plan ordering.
    ids.sort_unstable();
    ids.dedup();
    let map = assign_transcript_uids(&ids);
    if let Ok(mut guard) = TRANSCRIPT_UID_CACHE.lock() {
        guard.insert(cache_key, map.clone());
    }
    Ok(map)
}

/// Append a non-null `transcript_uid: UInt32` column to a transcript batch,
/// looked up from `uid_map` by the row's `stable_id`. Errors if any id is absent
/// (the map is built from the same transcript source, so this must not happen).
fn attach_transcript_uid_batch(
    batch: RecordBatch,
    source_type: &str,
    cache_version: &str,
    uid_map: &HashMap<String, u32>,
) -> Result<RecordBatch> {
    let schema = batch.schema();
    let sid_idx = schema.index_of("stable_id").map_err(|_| {
        DataFusionError::Execution("transcript batch missing stable_id column".to_string())
    })?;
    let sid = batch.column(sid_idx);

    let mut uids = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        let Some(id) = string_at(sid.as_ref(), row) else {
            return Err(DataFusionError::Execution(
                "transcript stable_id is null".to_string(),
            ));
        };
        let uid = *uid_map.get(&id).ok_or_else(|| {
            DataFusionError::Execution(format!("transcript stable_id {id} missing from uid map"))
        })?;
        uids.push(uid);
    }

    let mut fields = schema
        .fields()
        .iter()
        .map(|field| field.as_ref().clone())
        .collect::<Vec<Field>>();
    fields.push(Field::new("transcript_uid", DataType::UInt32, false));
    let mut columns = batch.columns().to_vec();
    columns.push(Arc::new(UInt32Array::from(uids)) as ArrayRef);

    let schema_with_upstream_metadata =
        Schema::new_with_metadata(fields, schema.metadata().clone());
    let new_schema =
        with_cache_identity_metadata(&schema_with_upstream_metadata, source_type, cache_version);
    RecordBatch::try_new(Arc::new(new_schema), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to attach transcript_uid column: {err}"))
    })
}

/// Position-sliced translation_sift schema: one row per `(transcript, position)`.
/// `key = (transcript_uid << 32) | position`; small `sift`/`poly` payloads use
/// miniblock + zstd (v2.1), not FullZip.
fn compact_translation_sift_position_schema(source_type: &str, cache_version: &str) -> Schema {
    let key = Field::new("key", DataType::UInt64, false);
    let sift = Field::new("sift", DataType::Binary, false);
    let poly = Field::new("poly", DataType::Binary, false);
    let schema = with_cache_identity_metadata(
        &Schema::new(vec![key, sift, poly]),
        source_type,
        cache_version,
    );
    let mut md = schema.metadata().clone();
    md.insert(
        crate::cache_common::SIFT_BLOB_VERSION_KEY.to_string(),
        "2".to_string(),
    );
    Schema::new_with_metadata(schema.fields().clone(), md)
}

/// Explode each transcript's SIFT/PolyPhen matrices into one row per protein
/// position. `(transcript_uid, position)` is unique by construction (uid 1:1,
/// positions grouped), so `key` is a unique BTree key.
fn transform_translation_sift_position_batch(
    batch: RecordBatch,
    source_type: &str,
    cache_version: &str,
    uid_map: &HashMap<String, u32>,
) -> Result<RecordBatch> {
    let schema = batch.schema();
    let tid_idx = schema.index_of("transcript_id")?;
    let sift_idx = schema.index_of("sift_predictions").ok();
    let poly_idx = schema.index_of("polyphen_predictions").ok();

    let mut keys: Vec<u64> = Vec::new();
    let mut sift_blobs: Vec<Vec<u8>> = Vec::new();
    let mut poly_blobs: Vec<Vec<u8>> = Vec::new();

    for row in 0..batch.num_rows() {
        let Some(transcript_id) = string_at(batch.column(tid_idx).as_ref(), row) else {
            continue;
        };
        let uid = *uid_map.get(&transcript_id).ok_or_else(|| {
            DataFusionError::Execution(format!(
                "translation transcript_id {transcript_id} missing from transcript uid map"
            ))
        })?;
        let mut preds = CachedPredictions::default();
        if let Some(idx) = sift_idx {
            preds.sift = read_compact_predictions(batch.column(idx).as_ref(), row);
        }
        if let Some(idx) = poly_idx {
            preds.polyphen = read_compact_predictions(batch.column(idx).as_ref(), row);
        }
        preds.sort(); // sorts by (position, amino_acid) → position-major
        append_translation_sift_position_rows(
            uid,
            &preds,
            &mut keys,
            &mut sift_blobs,
            &mut poly_blobs,
        )?;
    }

    let target_schema = Arc::new(compact_translation_sift_position_schema(
        source_type,
        cache_version,
    ));
    RecordBatch::try_new(
        target_schema,
        vec![
            Arc::new(UInt64Array::from(keys)) as ArrayRef,
            Arc::new(BinaryArray::from_iter_values(sift_blobs)) as ArrayRef,
            Arc::new(BinaryArray::from_iter_values(poly_blobs)) as ArrayRef,
        ],
    )
    .map_err(|err| {
        DataFusionError::Execution(format!(
            "failed to build position-sliced translation_sift batch: {err}"
        ))
    })
}

/// Linear two-pointer merge over `preds.sift`/`preds.polyphen` (each sorted by
/// position) emitting one `(key, sift, poly)` row per distinct position. Position
/// payloads drop the position field (it lives in `key`).
fn append_translation_sift_position_rows(
    uid: u32,
    preds: &CachedPredictions,
    keys: &mut Vec<u64>,
    sift_blobs: &mut Vec<Vec<u8>>,
    poly_blobs: &mut Vec<Vec<u8>>,
) -> Result<()> {
    let sift = &preds.sift;
    let poly = &preds.polyphen;
    let mut si = 0usize;
    let mut pi = 0usize;
    while si < sift.len() || pi < poly.len() {
        let next_pos = match (sift.get(si), poly.get(pi)) {
            (Some(s), Some(p)) => s.position.min(p.position),
            (Some(s), None) => s.position,
            (None, Some(p)) => p.position,
            (None, None) => break,
        };
        if next_pos < 0 {
            return Err(DataFusionError::Execution(format!(
                "negative protein position {next_pos} in translation_sift"
            )));
        }
        let s_start = si;
        while si < sift.len() && sift[si].position == next_pos {
            si += 1;
        }
        let p_start = pi;
        while pi < poly.len() && poly[pi].position == next_pos {
            pi += 1;
        }
        keys.push(((uid as u64) << 32) | (next_pos as u64));
        sift_blobs.push(crate::cache_common::serialize_position_entries_v2(
            &sift[s_start..si],
            crate::cache_common::SIFT_CODEC,
        )?);
        poly_blobs.push(crate::cache_common::serialize_position_entries_v2(
            &poly[p_start..pi],
            crate::cache_common::POLY_CODEC,
        )?);
    }
    Ok(())
}

fn drop_row_number_batch(batch: RecordBatch) -> Result<RecordBatch> {
    let schema = batch.schema();
    let fields = schema
        .fields()
        .iter()
        .filter(|field| field.name() != "_rn")
        .map(|field| field.as_ref().clone())
        .collect::<Vec<_>>();
    let columns = schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(_, field)| field.name() != "_rn")
        .map(|(index, _)| batch.column(index).clone())
        .collect::<Vec<_>>();
    RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).map_err(|err| {
        DataFusionError::Execution(format!("failed to drop translation row number: {err}"))
    })
}

fn project_batch_to_schema(batch: RecordBatch, target_schema: SchemaRef) -> Result<RecordBatch> {
    let source_schema = batch.schema();
    let mut columns = Vec::<ArrayRef>::with_capacity(target_schema.fields().len());
    for field in target_schema.fields() {
        let (index, _) = source_schema
            .column_with_name(field.name())
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "column '{}' not found in source batch",
                    field.name()
                ))
            })?;
        columns.push(batch.column(index).clone());
    }
    RecordBatch::try_new(target_schema, columns)
        .map_err(|err| DataFusionError::Execution(format!("failed to project Lance batch: {err}")))
}

fn chroms_from_schema(schema: &SchemaRef) -> Vec<String> {
    let mut chroms: Vec<String> = schema
        .metadata()
        .get("bio.vep.chromosomes")
        .and_then(|json| serde_json::from_str(json).ok())
        .unwrap_or_else(|| {
            DEFAULT_CHROMS
                .iter()
                .map(|chrom| (*chrom).to_string())
                .collect()
        });
    chroms.sort_by(|left, right| chrom_cache_order_key(left).cmp(&chrom_cache_order_key(right)));
    chroms
}

fn chrom_cache_order_key(chrom: &str) -> (usize, &str) {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Some(index) = DEFAULT_CHROMS.iter().position(|known| *known == bare) {
        (index, "")
    } else {
        (DEFAULT_CHROMS.len(), chrom)
    }
}

fn entity_output_name(kind: EnsemblEntityKind) -> &'static str {
    match kind {
        EnsemblEntityKind::Variation => "variation",
        EnsemblEntityKind::Transcript => "transcript",
        EnsemblEntityKind::Exon => "exon",
        EnsemblEntityKind::Translation => "translation",
        EnsemblEntityKind::RegulatoryFeature => "regulatory",
        EnsemblEntityKind::MotifFeature => "motif",
    }
}

fn entity_table_name(kind: EnsemblEntityKind) -> &'static str {
    match kind {
        EnsemblEntityKind::Variation => "var",
        EnsemblEntityKind::Transcript => "tx",
        EnsemblEntityKind::Exon => "exon",
        EnsemblEntityKind::Translation => "tl",
        EnsemblEntityKind::RegulatoryFeature => "reg",
        EnsemblEntityKind::MotifFeature => "motif",
    }
}

fn transcript_region_start_expr(start_col: &str) -> String {
    format!(
        "(CAST(FLOOR(({start_col} - 1) / {VEP_CACHE_REGION_SIZE_BP}.0) AS BIGINT) * {VEP_CACHE_REGION_SIZE_BP} + 1)"
    )
}

fn translation_source_region_preference_expr(start_col: &str, source_file_col: &str) -> String {
    let region_start = transcript_region_start_expr(start_col);
    let region_end = format!("({region_start} + {} - 1)", VEP_CACHE_REGION_SIZE_BP);
    format!(
        "CASE WHEN {source_file_col} LIKE CONCAT('%/', CAST({region_start} AS VARCHAR), '-', CAST({region_end} AS VARCHAR), '.gz') THEN 0 ELSE 1 END"
    )
}

/// Numeric start of the 1 Mb cache region a row was read from, parsed out of
/// the `<start>-<end>.gz` source file name.
///
/// Lexicographic ordering on `source_file` is unsafe here: `10000001-11000000`
/// sorts before `9000001-10000000` as text. `TRY_CAST` yields NULL when the
/// name does not match, so a non-region source degrades to the fallback terms
/// rather than erroring.
fn translation_source_region_start_expr(source_file_col: &str) -> String {
    format!(
        "TRY_CAST(REGEXP_REPLACE({source_file_col}, '^.*/([0-9]+)-[0-9]+\\.gz$', '$1') AS BIGINT)"
    )
}

fn build_translation_dedup_query_with_where_clause(table_name: &str, where_clause: &str) -> String {
    // Ensembl VEP resolves a transcript duplicated across cache regions in
    // `merge_features`, keeping the FIRST copy encountered over a feature list
    // ordered by region, so for a whole-chromosome run the lowest region wins.
    //
    // Preferring the region that contains the transcript `start` instead picks
    // the wrong copy whenever a transcript is also present in an earlier
    // region, which changes `protein_features` order and hence the `DOMAINS`
    // CSQ field. The start-region preference is retained as a fallback for
    // rows whose source file name does not parse.
    //
    // Traceability:
    // - Ensembl VEP `AnnotationType/Transcript.pm::merge_features()`
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L246-L280>
    // - Ensembl VEP `AnnotationSource.pm::get_all_features_by_InputBuffer()`
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm#L109-L145>
    let region_start = translation_source_region_start_expr("source_file");
    let source_pref = translation_source_region_preference_expr("start", "source_file");
    format!(
        "SELECT * FROM (\
            SELECT *, ROW_NUMBER() OVER (\
                PARTITION BY chrom, transcript_id \
                ORDER BY {region_start} NULLS LAST, {source_pref}, cdna_coding_start NULLS LAST, source_file\
            ) AS _rn \
            FROM {table_name}{where_clause}\
        ) WHERE _rn = 1"
    )
}

fn sql_escape_literal(value: &str) -> String {
    value.replace('\'', "''")
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{
        ArrayRef, BinaryArray, Int8Array, Int32Array, Int64Array, StringArray, UInt32Array,
    };
    use datafusion::arrow::compute::SortOptions;
    use datafusion::arrow::datatypes::{DataType, Field};
    use datafusion::physical_expr::expressions::Column;
    use datafusion::physical_expr_common::sort_expr::{LexOrdering, PhysicalSortExpr};
    #[allow(deprecated)]
    use datafusion::physical_plan::coalesce_batches::CoalesceBatchesExec;
    use datafusion::physical_plan::sorts::sort_preserving_merge::SortPreservingMergeExec;
    use datafusion::physical_plan::test::TestMemoryExec;

    #[test]
    fn variation_transform_drops_legacy_helpers_and_casts_positions() {
        let mut fields = Vec::new();
        let mut columns = Vec::<ArrayRef>::new();
        for name in VARIATION_REQUIRED_COLUMNS {
            match *name {
                "start" | "end" => {
                    fields.push(Field::new(*name, DataType::Int64, false));
                    columns.push(Arc::new(Int64Array::from(vec![10])) as ArrayRef);
                }
                "failed" => {
                    fields.push(Field::new(*name, DataType::Int8, true));
                    columns.push(
                        Arc::new(datafusion::arrow::array::Int8Array::from(vec![Some(0)]))
                            as ArrayRef,
                    );
                }
                _ => {
                    fields.push(Field::new(*name, DataType::Utf8, true));
                    columns.push(Arc::new(StringArray::from(vec![Some("value")])) as ArrayRef);
                }
            }
        }
        fields.push(Field::new("position_key", DataType::UInt64, false));
        columns.push(Arc::new(datafusion::arrow::array::UInt64Array::from(vec![10])) as ArrayRef);
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).unwrap();

        // start=10 is warm; the warm pass (keep_tier=0) keeps it and tags tier=0.
        let warm_starts = BTreeSet::from([10i64]);
        let warm =
            transform_variation_tier_batch(batch.clone(), "ensembl", &warm_starts, 0).unwrap();
        let schema = warm.schema();

        assert_eq!(
            schema.field_with_name("start").unwrap().data_type(),
            &DataType::UInt32
        );
        assert_eq!(
            schema.field_with_name("end").unwrap().data_type(),
            &DataType::UInt32
        );
        assert_eq!(
            schema.field_with_name("tier").unwrap().data_type(),
            &DataType::Int8
        );
        assert!(schema.field_with_name("position_key").is_err());
        assert_eq!(
            schema
                .metadata()
                .get(crate::cache_source::CACHE_SOURCE_METADATA_KEY),
            Some(&"ensembl".to_string())
        );
        assert_eq!(warm.num_rows(), 1);
        let tier = warm
            .column(schema.index_of("tier").unwrap())
            .as_any()
            .downcast_ref::<Int8Array>()
            .unwrap();
        assert_eq!(tier.value(0), 0);

        // The cold pass (keep_tier=1) drops the warm row entirely.
        let cold = transform_variation_tier_batch(batch, "ensembl", &warm_starts, 1).unwrap();
        assert_eq!(cold.num_rows(), 0);
    }

    #[test]
    fn variation_transform_preserves_optional_clin_sig_ref_allele_when_present() {
        let mut fields = Vec::new();
        let mut columns = Vec::<ArrayRef>::new();
        for name in VARIATION_REQUIRED_COLUMNS {
            match *name {
                "start" | "end" => {
                    fields.push(Field::new(*name, DataType::Int64, false));
                    columns.push(Arc::new(Int64Array::from(vec![10])) as ArrayRef);
                }
                "failed" => {
                    fields.push(Field::new(*name, DataType::Int8, true));
                    columns.push(Arc::new(Int8Array::from(vec![Some(0)])) as ArrayRef);
                }
                _ => {
                    fields.push(Field::new(*name, DataType::Utf8, true));
                    columns.push(Arc::new(StringArray::from(vec![Some("value")])) as ArrayRef);
                }
            }
            if *name == "clin_sig_allele" {
                fields.push(Field::new("clin_sig_ref_allele", DataType::Utf8, false));
                columns.push(Arc::new(StringArray::from(vec!["AGG"])) as ArrayRef);
            }
        }
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), columns).unwrap();
        let transformed =
            transform_variation_tier_batch(batch, "ensembl", &BTreeSet::from([10]), 0).unwrap();
        let schema = transformed.schema();
        let field = schema.field_with_name("clin_sig_ref_allele").unwrap();
        assert_eq!(field.data_type(), &DataType::Utf8);
        assert!(field.is_nullable());
        let values = transformed
            .column(schema.index_of("clin_sig_ref_allele").unwrap())
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(values.value(0), "AGG");
    }

    #[test]
    fn assign_transcript_uids_is_dense_and_unique() {
        let mut ids = vec![
            "ENST2".to_string(),
            "ENST1".to_string(),
            "ENST1".to_string(),
        ];
        ids.sort_unstable();
        ids.dedup();
        let map = assign_transcript_uids(&ids);
        assert_eq!(map.len(), 2);
        assert_eq!(map["ENST1"], 0);
        assert_eq!(map["ENST2"], 1);
    }

    #[test]
    fn attach_transcript_uid_appends_column_and_errors_on_missing() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            Arc::clone(&schema),
            vec![
                Arc::new(StringArray::from(vec!["ENST1", "ENST2"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![10_i64, 20])) as ArrayRef,
            ],
        )
        .unwrap();
        let uid_map: HashMap<String, u32> =
            [("ENST1".to_string(), 0u32), ("ENST2".to_string(), 1u32)]
                .into_iter()
                .collect();

        let out = attach_transcript_uid_batch(batch, "ensembl", "115", &uid_map).unwrap();
        let out_schema = out.schema();
        let uid_idx = out_schema.index_of("transcript_uid").unwrap();
        assert_eq!(out_schema.field(uid_idx).data_type(), &DataType::UInt32);
        let uids = out
            .column(uid_idx)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        assert_eq!(uids.value(0), 0);
        assert_eq!(uids.value(1), 1);
        assert_eq!(
            out_schema
                .metadata()
                .get(crate::cache_source::CACHE_SOURCE_METADATA_KEY),
            Some(&"ensembl".to_string())
        );
        assert_eq!(
            out_schema
                .metadata()
                .get(crate::cache_identity::CACHE_VERSION_METADATA_KEY),
            Some(&"115".to_string())
        );

        // An id absent from the map is a hard error (consistency guard).
        let orphan = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST_MISSING"])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1_i64])) as ArrayRef,
            ],
        )
        .unwrap();
        assert!(attach_transcript_uid_batch(orphan, "ensembl", "115", &uid_map).is_err());
    }

    #[test]
    fn translation_dedup_query_stays_chrom_scoped() {
        let query = build_translation_dedup_query_with_where_clause("tl", " WHERE chrom = 'chr1'");
        assert!(query.contains("PARTITION BY chrom, transcript_id"));
        assert!(query.contains("WHERE chrom = 'chr1'"));
    }

    /// A transcript present in more than one cache region must resolve to the
    /// copy from the LOWEST region, matching VEP's merge_features first-wins
    /// over a region-ordered feature list. Preferring the region containing the
    /// transcript start picks the later copy and reorders protein_features,
    /// which surfaces as a DOMAINS mismatch (chrX ENST00000381077).
    #[test]
    fn translation_dedup_prefers_the_lowest_source_region() {
        let query = build_translation_dedup_query_with_where_clause("tl", "");
        let order = query
            .split("ORDER BY ")
            .nth(1)
            .expect("ORDER BY present")
            .to_string();
        let region_term = order
            .find("REGEXP_REPLACE")
            .expect("region-start term present");
        let start_pref = order
            .find("CASE WHEN")
            .expect("start-region fallback present");
        assert!(
            region_term < start_pref,
            "lowest-region term must outrank the start-region preference: {order}"
        );
        assert!(
            order.contains("TRY_CAST"),
            "must degrade on unparsable names"
        );
    }

    #[test]
    fn transcript_context_query_uses_provider_schema() {
        let schema = Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("source_file", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
        ]);

        let query = build_export_query(
            EnsemblEntityKind::Transcript,
            "tx",
            Some("1"),
            Some(&schema),
        );

        assert!(query.contains("ROW_NUMBER()"));
        assert!(query.contains("PARTITION BY chrom, stable_id"));
        assert!(query.contains("WHERE chrom = '1'"));
    }

    #[test]
    fn chroms_from_schema_uses_old_cache_order() {
        let metadata = [(
            "bio.vep.chromosomes".to_string(),
            serde_json::json!(["10", "2", "1", "X", "MT", "HSCHR6_MHC_COX_CTG1"]).to_string(),
        )]
        .into_iter()
        .collect();
        let schema = Arc::new(Schema::new_with_metadata(Vec::<Field>::new(), metadata));

        assert_eq!(
            chroms_from_schema(&schema),
            ["1", "2", "10", "X", "MT", "HSCHR6_MHC_COX_CTG1"]
        );
    }

    #[test]
    fn filter_chroms_matches_bare_labels_and_passes_through_when_none() {
        let all = || {
            vec![
                "1".to_string(),
                "2".to_string(),
                "X".to_string(),
                "chr10".to_string(),
            ]
        };
        // None → unchanged.
        assert_eq!(filter_chroms(all(), None), all());
        // `chr1` matches the bare `1`; `2` matches itself; `chr10` matches `10`.
        let allow = vec!["chr1".to_string(), "2".to_string(), "10".to_string()];
        assert_eq!(
            filter_chroms(all(), Some(&allow)),
            ["1".to_string(), "2".to_string(), "chr10".to_string()]
        );
        // No overlap → empty.
        let none = vec!["chrY".to_string()];
        assert!(filter_chroms(all(), Some(&none)).is_empty());
    }

    #[test]
    fn translation_sift_position_rows_group_by_position_with_unique_keys() {
        use crate::cache_common::{SiftBlobVersion, deserialize_position_predictions_versioned};
        use crate::transcript_consequence::CompactPrediction;

        // sift at positions 1 & 2, polyphen only at position 2.
        let mut preds = CachedPredictions {
            sift: vec![
                CompactPrediction {
                    position: 1,
                    amino_acid: 0,
                    prediction: 0,
                    score: 0.3,
                },
                CompactPrediction {
                    position: 2,
                    amino_acid: 2,
                    prediction: 1,
                    score: 0.01,
                },
            ],
            polyphen: vec![CompactPrediction {
                position: 2,
                amino_acid: 2,
                prediction: 4,
                score: 0.9,
            }],
        };
        preds.sort();

        let uid = 7u32;
        let mut keys = Vec::new();
        let mut sift_blobs = Vec::new();
        let mut poly_blobs = Vec::new();
        append_translation_sift_position_rows(
            uid,
            &preds,
            &mut keys,
            &mut sift_blobs,
            &mut poly_blobs,
        )
        .unwrap();

        // One row per distinct position (1 and 2), keyed by (uid<<32)|position.
        assert_eq!(
            keys,
            vec![((uid as u64) << 32) | 1, ((uid as u64) << 32) | 2]
        );
        // Keys are unique.
        let mut sorted = keys.clone();
        sorted.dedup();
        assert_eq!(sorted.len(), keys.len());

        // Position 1: sift only, no polyphen.
        let p1 = deserialize_position_predictions_versioned(
            1,
            &sift_blobs[0],
            &poly_blobs[0],
            SiftBlobVersion::V2DivIndex,
        )
        .unwrap();
        assert_eq!(p1.sift.len(), 1);
        assert!(p1.polyphen.is_empty());
        assert_eq!(p1.sift[0].position, 1);

        // Position 2: both sift and polyphen.
        let p2 = deserialize_position_predictions_versioned(
            2,
            &sift_blobs[1],
            &poly_blobs[1],
            SiftBlobVersion::V2DivIndex,
        )
        .unwrap();
        assert_eq!(p2.sift.len(), 1);
        assert_eq!(p2.polyphen.len(), 1);
        assert_eq!(p2.polyphen[0].prediction, 4);

        // Schema is key/sift/poly.
        let schema = compact_translation_sift_position_schema("ensembl", "115");
        assert_eq!(schema.field(0).name(), "key");
        assert_eq!(schema.field(0).data_type(), &DataType::UInt64);
        assert_eq!(schema.field(1).name(), "sift");
        assert_eq!(schema.field(2).name(), "poly");
    }

    #[test]
    fn append_translation_sift_position_rows_emits_v2() {
        use crate::cache_common::{POLY_CODEC, SIFT_CODEC, deserialize_position_entries_v2};
        use crate::transcript_consequence::CompactPrediction;
        let preds = CachedPredictions {
            sift: vec![
                CompactPrediction {
                    position: 1,
                    amino_acid: 0,
                    prediction: 1,
                    score: 0.03,
                },
                CompactPrediction {
                    position: 1,
                    amino_acid: 5,
                    prediction: 0,
                    score: 1.0,
                },
            ],
            polyphen: vec![CompactPrediction {
                position: 1,
                amino_acid: 0,
                prediction: 2,
                score: 0.998,
            }],
        };
        let mut keys = Vec::new();
        let mut sift_blobs = Vec::new();
        let mut poly_blobs = Vec::new();
        append_translation_sift_position_rows(
            42,
            &preds,
            &mut keys,
            &mut sift_blobs,
            &mut poly_blobs,
        )
        .unwrap();
        assert_eq!(keys, vec![(42u64 << 32) | 1]);
        // v2 sift: 2 entries * 3 bytes
        assert_eq!(sift_blobs[0].len(), 6);
        let s = deserialize_position_entries_v2(1, &sift_blobs[0], SIFT_CODEC).unwrap();
        assert_eq!(s[1].score.to_bits(), 1.0f32.to_bits());
        let p = deserialize_position_entries_v2(1, &poly_blobs[0], POLY_CODEC).unwrap();
        assert_eq!(p[0].score.to_bits(), 0.998f32.to_bits());
    }

    #[test]
    fn sift_position_schema_carries_v2_flag() {
        let schema = compact_translation_sift_position_schema("ensembl", "115");
        assert_eq!(
            schema
                .metadata()
                .get(crate::cache_common::SIFT_BLOB_VERSION_KEY)
                .map(String::as_str),
            Some("2")
        );
    }
}
