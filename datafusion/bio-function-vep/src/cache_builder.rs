//! Cache builder: converts raw Ensembl VEP cache to parquet + fjall.
//!
//! Ported from vepyr `convert.rs` and extended with:
//! - Fjall dual-sink for `variation` (single pass via `start_ingestion()`)
//! - Fjall second pass for `translation_sift` (re-sorted by `transcript_id`)
//! - Progress callback for driving tqdm bars in Python wrappers

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{
    Array, AsArray, LargeStringArray, RecordBatch, StringArray, StringViewArray,
};
use datafusion::arrow::datatypes::{DataType, Int64Type, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::parquet::arrow::ArrowWriter;
use datafusion::parquet::basic::Compression;
use datafusion::parquet::file::properties::WriterProperties;
use datafusion::parquet::format::SortingColumn;
use datafusion::parquet::schema::types::ColumnPath;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{SessionConfig, SessionContext};
use futures::StreamExt;
use log::info;

use datafusion_bio_format_ensembl_cache::{
    EnsemblCacheOptions, EnsemblCacheTableProvider, EnsemblEntityKind,
};

use crate::annotate_provider::read_compact_predictions;
use crate::kv_cache::LoadStats;
use crate::kv_cache::key_encoding::{chrom_to_code, encode_position_key};
use crate::kv_cache::kv_store::VepKvStore;
use crate::kv_cache::position_entry::serialize_position_entry;
use crate::kv_cache::sift_store::SiftKvStore;
use crate::transcript_consequence::CachedPredictions;

/// Progress callback: `(entity, format, batch_rows, total_rows, total_expected)`.
///
/// - `entity`: "variation", "transcript", "exon", etc.
/// - `format`: "parquet" or "fjall"
/// - `batch_rows`: rows processed in this batch
/// - `total_rows`: cumulative rows processed so far for this (entity, format)
/// - `total_expected`: total rows expected (0 if unknown)
pub type OnProgress = Box<dyn Fn(&str, &str, usize, usize, usize) + Send + Sync>;

/// Main chromosomes that get their own parquet file.
const MAIN_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y",
];

/// Chromosomes in fjall key encoding order (chrom_code ascending).
const CHROM_CODE_ORDER: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// Statistics returned after building one entity.
#[derive(Debug, Clone)]
pub struct EntityStats {
    pub entity: String,
    pub parquet_files: Vec<(String, usize)>,
    pub fjall_stats: Option<LoadStats>,
}

/// Builder for converting Ensembl VEP cache to parquet and fjall.
pub struct CacheBuilder {
    cache_root: String,
    output_dir: String,
    partitions: usize,
    build_fjall: bool,
    overwrite: bool,
    zstd_level: i32,
    dict_size_kb: u32,
    on_progress: Option<Arc<OnProgress>>,
}

impl CacheBuilder {
    pub fn new(cache_root: impl Into<String>, output_dir: impl Into<String>) -> Self {
        Self {
            cache_root: cache_root.into(),
            output_dir: output_dir.into(),
            partitions: 8,
            build_fjall: true,
            overwrite: false,
            zstd_level: 3,
            dict_size_kb: 112,
            on_progress: None,
        }
    }

    pub fn with_partitions(mut self, n: usize) -> Self {
        self.partitions = n;
        self
    }

    pub fn with_build_fjall(mut self, enabled: bool) -> Self {
        self.build_fjall = enabled;
        self
    }

    pub fn with_overwrite(mut self, overwrite: bool) -> Self {
        self.overwrite = overwrite;
        self
    }

    pub fn with_zstd_level(mut self, level: i32) -> Self {
        self.zstd_level = level;
        self
    }

    pub fn with_dict_size_kb(mut self, size: u32) -> Self {
        self.dict_size_kb = size;
        self
    }

    pub fn with_on_progress(mut self, cb: OnProgress) -> Self {
        self.on_progress = Some(Arc::new(cb));
        self
    }

    /// Build all entities (parquet + optional fjall).
    pub async fn build_all(&self) -> Result<Vec<EntityStats>> {
        let entities = [
            "variation",
            "transcript",
            "exon",
            "translation",
            "regulatory",
            "motif",
        ];
        let mut results = Vec::new();
        for entity in entities {
            match self.build_entity(entity).await {
                Ok(stats) => results.extend(stats),
                Err(e) => {
                    let msg = e.to_string();
                    if msg.contains("No source files discovered") || msg.contains("skipped") {
                        info!("{entity}: skipped (no source files)");
                    } else {
                        return Err(e);
                    }
                }
            }
        }
        Ok(results)
    }

    /// Build a single entity. Returns one or more EntityStats (translation splits into two).
    ///
    /// When `overwrite` is false (default), skips entities whose output
    /// already exists:
    /// - For parquet-only entities: skips if the parquet directory contains `.parquet` files
    /// - For variation: skips parquet if files exist, skips fjall if `variation.fjall` exists
    /// - For translation: skips parquet if files exist, skips sift fjall if `translation_sift.fjall` exists
    pub async fn build_entity(&self, entity: &str) -> Result<Vec<EntityStats>> {
        let kind = parse_entity(entity)
            .ok_or_else(|| DataFusionError::Execution(format!("Unknown entity: {entity}")))?;

        // Create output directories
        let subdir = entity_subdir(kind);
        let entity_dir = format!("{}/{}", self.output_dir, subdir);
        std::fs::create_dir_all(&entity_dir).map_err(|e| {
            DataFusionError::Execution(format!("Failed to create dir {entity_dir}: {e}"))
        })?;

        if kind == EnsemblEntityKind::Translation {
            std::fs::create_dir_all(format!("{}/translation_core", self.output_dir))
                .map_err(|e| DataFusionError::Execution(format!("Failed to create dir: {e}")))?;
            std::fs::create_dir_all(format!("{}/translation_sift", self.output_dir))
                .map_err(|e| DataFusionError::Execution(format!("Failed to create dir: {e}")))?;
        }

        // Skip if all outputs exist (unless overwrite is set).
        // Variation and translation have partial-skip logic in their
        // own build methods (e.g. parquet exists but fjall missing).
        if !self.overwrite && kind != EnsemblEntityKind::Variation && kind != EnsemblEntityKind::Translation {
            if dir_has_parquet_files(&entity_dir) {
                info!("{entity}: parquet already exists, skipping (use overwrite to rebuild)");
                return Ok(vec![EntityStats {
                    entity: subdir.to_string(),
                    parquet_files: vec![],
                    fjall_stats: None,
                }]);
            }
        }

        if kind == EnsemblEntityKind::Variation {
            return self.build_variation().await;
        }

        if kind == EnsemblEntityKind::Translation {
            return self.build_translation().await;
        }

        // All other entities: parquet only
        let results = self.build_parquet_entity(kind).await?;
        Ok(vec![EntityStats {
            entity: subdir.to_string(),
            parquet_files: results,
            fjall_stats: None,
        }])
    }

    /// Build variation: parquet + optional fjall.
    ///
    /// When parquet files already exist and fjall is missing, reads from the
    /// existing parquet files to populate fjall without re-parsing the source.
    async fn build_variation(&self) -> Result<Vec<EntityStats>> {
        let kind = EnsemblEntityKind::Variation;
        let table_name = "var";
        let parquet_dir = format!("{}/variation", self.output_dir);
        let parquet_exists = dir_has_parquet_files(&parquet_dir);
        let fjall_dir_path = format!("{}/variation.fjall", self.output_dir);
        let fjall_exists = Path::new(&fjall_dir_path).exists();

        if !self.overwrite {
            let need_parquet = !parquet_exists;
            let need_fjall = self.build_fjall && !fjall_exists;

            if !need_parquet && !need_fjall {
                info!("variation: all outputs exist, skipping (use overwrite to rebuild)");
                return Ok(vec![EntityStats {
                    entity: "variation".to_string(),
                    parquet_files: vec![],
                    fjall_stats: None,
                }]);
            }

            // Parquet exists but fjall missing → rebuild fjall from parquet
            if !need_parquet && need_fjall {
                info!("variation: parquet exists, building fjall from existing parquet files");
                return self.build_variation_fjall_from_parquet().await;
            }
        }

        // Discover chromosomes from schema metadata
        let init_ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
        let provider_schema = {
            let table = init_ctx.table(table_name).await?;
            table.schema().inner().clone()
        };
        let chroms = chroms_from_schema(&provider_schema);
        drop(init_ctx);

        let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
        let (main_chroms, other_chroms) = split_chroms(&chroms, &main_set);

        info!(
            "variation: {} main chroms, {} other contigs",
            main_chroms.len(),
            other_chroms.len()
        );

        let mut parquet_results: Vec<(String, usize)> = Vec::new();

        // --- Fjall setup ---
        let fjall_state = if self.build_fjall {
            let fjall_dir = format!("{}/variation.fjall", self.output_dir);
            if Path::new(&fjall_dir).exists() && !self.overwrite {
                info!("variation.fjall already exists, skipping fjall build (use overwrite to rebuild)");
                None
            } else {
                // We need a schema for the VepKvStore — get it from a probe query
                let probe_ctx =
                    make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
                let probe_df = probe_ctx
                    .sql(&format!("SELECT * FROM {table_name} LIMIT 0"))
                    .await?;
                let schema: SchemaRef = Arc::new(probe_df.schema().as_arrow().clone());
                drop(probe_ctx);

                let store = VepKvStore::create(Path::new(&fjall_dir), schema.clone())?;

                // Train zstd dictionary from sample
                let dict = self
                    .train_variation_dict(&store, &schema, kind, table_name)
                    .await?;

                Some(FjallVariationState {
                    store,
                    schema,
                    dict,
                    fjall_dir,
                })
            }
        } else {
            None
        };

        let fjall_start_time = Instant::now();

        // Start ingestion if fjall is active
        let mut ingestion = if let Some(ref state) = fjall_state {
            Some(
                state
                    .store
                    .data_partition()
                    .start_ingestion()
                    .map_err(|e| DataFusionError::External(Box::new(e)))?,
            )
        } else {
            None
        };

        let mut compressor = if let Some(ref state) = fjall_state {
            if let Some(ref dict) = state.dict {
                Some(
                    zstd::bulk::Compressor::with_dictionary(self.zstd_level, dict).map_err(
                        |e| {
                            DataFusionError::Execution(format!(
                                "failed to create zstd compressor: {e}"
                            ))
                        },
                    )?,
                )
            } else {
                None
            }
        } else {
            None
        };

        let mut fjall_total_positions = 0u64;
        let mut fjall_total_variants = 0u64;
        let mut fjall_total_bytes = 0u64;

        // Column indices for fjall (cached once since schema is constant)
        let fjall_col_info = fjall_state.as_ref().map(|state| {
            let chrom_col_idx = state.schema.index_of("chrom").unwrap();
            let start_col_idx = state.schema.index_of("start").unwrap();
            let allele_col_idx = state.schema.index_of("allele_string").unwrap();
            let col_indices: Vec<usize> = (0..state.schema.fields().len())
                .filter(|&i| i != chrom_col_idx && i != start_col_idx)
                .collect();
            (chrom_col_idx, start_col_idx, allele_col_idx, col_indices)
        });

        // Build a unified processing list: main chroms (in chrom_code order)
        // then non-main contigs (sorted by hash-based chrom_code for ascending
        // fjall keys). Main chroms get individual parquet files; other contigs
        // share other.parquet. All contigs feed fjall.
        //
        // Each entry: (chrom, output_path, is_other)
        let mut chrom_batches: Vec<(String, String, bool)> = Vec::new();
        for chrom in CHROM_CODE_ORDER {
            if main_chroms.iter().any(|c| c == chrom) {
                let out = format!("{}/variation/chr{chrom}.parquet", self.output_dir);
                chrom_batches.push((chrom.to_string(), out, false));
            }
        }

        // Register non-canonical contigs with collision-free sequential codes,
        // then sort by the assigned code for ascending fjall key order.
        if !other_chroms.is_empty() {
            use crate::kv_cache::key_encoding::register_non_canonical_contigs;

            let refs: Vec<&str> = other_chroms.iter().map(|s| s.as_str()).collect();
            let mapping = register_non_canonical_contigs(&refs);

            // Store the mapping in fjall for read-time lookup
            if let Some(ref state) = fjall_state {
                state.store.store_contig_codes(&mapping)?;
            }

            info!(
                "variation: {} non-main contigs registered ({} codes) → other.parquet + fjall",
                other_chroms.len(),
                mapping.len()
            );
        }

        // Now sort other_chroms by their (collision-free) code
        let mut sorted_other = other_chroms.clone();
        sorted_other.sort_by_key(|c| chrom_to_code(c));
        if !sorted_other.is_empty() {
            let other_out = format!("{}/variation/other.parquet", self.output_dir);
            for chrom in &sorted_other {
                chrom_batches.push((chrom.clone(), other_out.clone(), true));
            }
        }

        let mut accum: Option<PositionAccumulator> = None;
        let mut total_parquet_rows: usize = 0;

        // Shared writer for other.parquet (opened lazily)
        let mut other_writer: Option<ArrowWriter<File>> = None;
        let mut other_total_rows = 0usize;

        for (chrom, output_file, is_other) in &chrom_batches {
            let ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
            let query = build_query(kind, table_name, Some(chrom));

            let df = ctx.sql(&query).await?;
            let plan = df.create_physical_plan().await?;

            // Check if the plan has a SortPreservingMergeExec wrapping a
            // multi-partition inner plan.  If so, use parallel execution.
            let inner_plan = extract_multi_partition_inner(&plan);
            let num_partitions = inner_plan
                .as_ref()
                .map(|p| p.properties().partitioning.partition_count())
                .unwrap_or(1);

            info!(
                "variation: [{chrom}] {} partitions (parallel={})",
                num_partitions,
                inner_plan.is_some()
            );

            let schema = plan.schema();
            let sk = sort_key(kind);

            // For main chroms: create a new writer per chrom.
            // For other contigs: use shared other_writer (opened lazily).
            let mut main_writer: Option<ArrowWriter<File>> = None;
            let writer = if *is_other {
                if other_writer.is_none() {
                    other_writer = Some(create_writer(output_file, &schema, kind, sk, None)?);
                }
                other_writer.as_mut().unwrap()
            } else {
                main_writer = Some(create_writer(output_file, &schema, kind, sk, None)?);
                main_writer.as_mut().unwrap()
            };

            let mut chrom_rows = 0usize;

            if let Some(inner) = inner_plan.filter(|_| num_partitions > 1) {
                // --- Parallel path: write N temp parquet files, then merge ---
                let temp_dir = format!("{}/variation/_tmp_{chrom}", self.output_dir);
                std::fs::create_dir_all(&temp_dir).map_err(|e| {
                    DataFusionError::Execution(format!("Failed to create temp dir: {e}"))
                })?;

                // Phase 1: Execute all partitions in parallel, each writes a temp parquet
                let task_ctx = ctx.task_ctx();
                let mut handles = tokio::task::JoinSet::new();
                for partition_idx in 0..num_partitions {
                    let stream = inner.execute(partition_idx, Arc::clone(&task_ctx))?;
                    let part_schema = stream.schema();
                    let part_file = format!("{temp_dir}/part_{partition_idx}.parquet");
                    let part_sk = sort_key(kind);
                    handles.spawn(async move {
                        let mut pw = create_writer(
                            &part_file,
                            &part_schema,
                            kind,
                            part_sk,
                            None,
                        )?;
                        let mut rows = 0usize;
                        let mut stream = stream;
                        while let Some(batch_result) = stream.next().await {
                            let batch = batch_result?;
                            if batch.num_rows() == 0 {
                                continue;
                            }
                            pw.write(&batch)?;
                            rows += batch.num_rows();
                        }
                        pw.close().map_err(|e| {
                            DataFusionError::Execution(format!(
                                "Failed to close temp parquet: {e}"
                            ))
                        })?;
                        // Return partition_idx so we can sort by it after JoinSet
                        Ok::<(usize, String, usize), DataFusionError>((
                            partition_idx, part_file, rows,
                        ))
                    });
                }

                // Collect results — JoinSet returns in completion order
                let mut part_files: Vec<(usize, String, usize)> = Vec::new();
                while let Some(result) = handles.join_next().await {
                    let (idx, file, rows) = result
                        .map_err(|e| DataFusionError::Execution(format!("Task join error: {e}")))?
                        .map_err(|e: DataFusionError| e)?;
                    part_files.push((idx, file, rows));
                }
                // Sort by partition index for ascending position order
                part_files.sort_by_key(|(idx, _, _)| *idx);

                let parallel_rows: usize = part_files.iter().map(|(_, _, r)| *r).sum();
                info!(
                    "variation: [{chrom}] phase 1 done: {} rows in {} partitions",
                    format_rows(parallel_rows),
                    num_partitions
                );

                // Phase 2: Read temp parquet files in partition order → final parquet + fjall
                // Use target_partitions=1 to preserve row order within each file.
                let read_config = SessionConfig::new().with_target_partitions(1);
                let read_ctx = SessionContext::new_with_config(read_config);
                for (_, part_file, _) in &part_files {
                    read_ctx
                        .register_parquet("_part", part_file, Default::default())
                        .await?;
                    let part_batches = read_ctx.sql("SELECT * FROM _part").await?.collect().await?;
                    read_ctx.deregister_table("_part")?;

                    for batch in &part_batches {
                        if batch.num_rows() == 0 {
                            continue;
                        }

                        writer.write(batch)?;
                        chrom_rows += batch.num_rows();
                        total_parquet_rows += batch.num_rows();

                        if let Some(ref cb) = self.on_progress {
                            cb(
                                "variation",
                                "parquet",
                                batch.num_rows(),
                                total_parquet_rows,
                                0,
                            );
                        }

                        // Feed fjall ingestion
                        if let (
                            Some(ing),
                            Some((chrom_col_idx, start_col_idx, allele_col_idx, col_indices)),
                        ) = (&mut ingestion, &fjall_col_info)
                        {
                            let starts =
                                batch.column(*start_col_idx).as_primitive::<Int64Type>();
                            let chrom_col = batch.column(*chrom_col_idx);

                            for row in 0..batch.num_rows() {
                                let start = starts.value(row);
                                let row_chrom = string_value(chrom_col.as_ref(), row);
                                let chrom_code = chrom_to_code(row_chrom);

                                let should_flush = accum.as_ref().is_some_and(|a| {
                                    a.chrom_code != chrom_code || a.start != start
                                });

                                if should_flush {
                                    let a = accum.take().unwrap();
                                    let (key, value) = a.finish_entry(
                                        col_indices,
                                        *allele_col_idx,
                                        &mut compressor,
                                    )?;
                                    ing.write(&key, &value)
                                        .map_err(|e| DataFusionError::External(Box::new(e)))?;
                                    fjall_total_positions += 1;
                                    fjall_total_bytes += value.len() as u64;
                                }

                                fjall_total_variants += 1;

                                match &mut accum {
                                    Some(a) => a.add_row(row, batch),
                                    None => {
                                        accum = Some(PositionAccumulator::new(
                                            chrom_code,
                                            row_chrom.to_string(),
                                            start,
                                            row,
                                            batch,
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }

                // Cleanup temp files
                if let Err(e) = std::fs::remove_dir_all(&temp_dir) {
                    log::warn!("Failed to clean up temp dir {temp_dir}: {e}");
                }
            } else {
                // --- Sequential path: single partition, stream directly ---
                let mut stream =
                    datafusion::physical_plan::execute_stream(plan, ctx.task_ctx())?;

                while let Some(batch_result) = stream.next().await {
                    let batch = batch_result?;
                    if batch.num_rows() == 0 {
                        continue;
                    }

                    writer.write(&batch)?;
                    chrom_rows += batch.num_rows();
                    total_parquet_rows += batch.num_rows();

                    if let Some(ref cb) = self.on_progress {
                        cb(
                            "variation",
                            "parquet",
                            batch.num_rows(),
                            total_parquet_rows,
                            0,
                        );
                    }

                    // Feed fjall ingestion
                    if let (
                        Some(ing),
                        Some((chrom_col_idx, start_col_idx, allele_col_idx, col_indices)),
                    ) = (&mut ingestion, &fjall_col_info)
                    {
                        let starts =
                            batch.column(*start_col_idx).as_primitive::<Int64Type>();
                        let chrom_col = batch.column(*chrom_col_idx);

                        for row in 0..batch.num_rows() {
                            let start = starts.value(row);
                            let row_chrom = string_value(chrom_col.as_ref(), row);
                            let chrom_code = chrom_to_code(row_chrom);

                            let should_flush = accum
                                .as_ref()
                                .is_some_and(|a| {
                                    a.chrom_code != chrom_code || a.start != start
                                });

                            if should_flush {
                                let a = accum.take().unwrap();
                                let (key, value) = a.finish_entry(
                                    col_indices,
                                    *allele_col_idx,
                                    &mut compressor,
                                )?;
                                ing.write(&key, &value)
                                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                                fjall_total_positions += 1;
                                fjall_total_bytes += value.len() as u64;
                            }

                            fjall_total_variants += 1;

                            match &mut accum {
                                Some(a) => a.add_row(row, &batch),
                                None => {
                                    accum = Some(PositionAccumulator::new(
                                        chrom_code,
                                        row_chrom.to_string(),
                                        start,
                                        row,
                                        &batch,
                                    ));
                                }
                            }
                        }
                    }
                }
            }

            // Flush fjall accumulator between chromosomes
            if let (Some(ing), Some((_, _, allele_col_idx, col_indices))) =
                (&mut ingestion, &fjall_col_info)
            {
                if let Some(a) = accum.take() {
                    let (key, value) =
                        a.finish_entry(col_indices, *allele_col_idx, &mut compressor)?;
                    ing.write(&key, &value)
                        .map_err(|e| DataFusionError::External(Box::new(e)))?;
                    fjall_total_positions += 1;
                    fjall_total_bytes += value.len() as u64;
                }
                if let Some(ref cb) = self.on_progress {
                    cb("variation", "fjall", 0, fjall_total_positions as usize, 0);
                }
            }

            // Close main-chrom writer; other_writer stays open across contigs
            if let Some(w) = main_writer {
                w.close().map_err(|e| {
                    DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
                })?;
            }

            if *is_other {
                other_total_rows += chrom_rows;
            } else if chrom_rows > 0 {
                parquet_results.push((output_file.clone(), chrom_rows));
                info!("variation: {chrom} {} rows", format_rows(chrom_rows));
            }
        }

        // Close shared other.parquet writer
        if let Some(w) = other_writer {
            w.close().map_err(|e| {
                DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
            })?;
        }
        if other_total_rows > 0 {
            let other_out = format!("{}/variation/other.parquet", self.output_dir);
            parquet_results.push((other_out, other_total_rows));
            info!(
                "variation: other {} rows ({} contigs)",
                format_rows(other_total_rows),
                sorted_other.len()
            );
        }

        // Finalize fjall
        let fjall_stats = if let Some(ref mut ing) = ingestion {
            // Take ingestion out to call finish()
            let ingestion = ingestion.take().unwrap();
            ingestion
                .finish()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;

            if let Some(ref state) = fjall_state {
                state.store.persist()?;

                info!("Running major compaction on variation.fjall...");
                let compact_start = Instant::now();
                state
                    .store
                    .data_partition()
                    .major_compact()
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                info!(
                    "Major compaction completed in {:.1}s",
                    compact_start.elapsed().as_secs_f64()
                );
            }

            Some(LoadStats {
                total_variants: fjall_total_variants,
                total_positions: fjall_total_positions,
                total_bytes: fjall_total_bytes,
                elapsed_secs: fjall_start_time.elapsed().as_secs_f64(),
            })
        } else {
            None
        };

        Ok(vec![EntityStats {
            entity: "variation".to_string(),
            parquet_files: parquet_results,
            fjall_stats,
        }])
    }

    /// Train zstd dictionary from a sample of variation data.
    async fn train_variation_dict(
        &self,
        store: &VepKvStore,
        schema: &SchemaRef,
        kind: EnsemblEntityKind,
        table_name: &str,
    ) -> Result<Option<Arc<Vec<u8>>>> {
        let ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
        let sample_df = ctx
            .sql(&format!("SELECT * FROM {table_name} LIMIT 10000"))
            .await?;
        let batches = sample_df.collect().await?;

        let chrom_col_idx = schema.index_of("chrom")?;
        let start_col_idx = schema.index_of("start")?;
        let allele_col_idx = schema.index_of("allele_string")?;
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_col_idx && i != start_col_idx)
            .collect();

        let mut samples: Vec<Vec<u8>> = Vec::new();
        for batch in &batches {
            if batch.num_rows() == 0 {
                continue;
            }
            let chrom_col = batch.column(chrom_col_idx);
            let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();

            let mut groups: HashMap<(String, i64), Vec<usize>> = HashMap::new();
            for row in 0..batch.num_rows() {
                let chrom = string_value(chrom_col.as_ref(), row);
                let start = starts.value(row);
                groups
                    .entry((chrom.to_string(), start))
                    .or_default()
                    .push(row);
            }

            for rows in groups.values() {
                let value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
                samples.push(value);
            }
        }

        let min_training_samples = 100;
        if samples.len() < min_training_samples {
            info!(
                "only {} sample entries, skipping dictionary training (need >= {})",
                samples.len(),
                min_training_samples
            );
            return Ok(None);
        }

        let max_samples = 10_000.min(samples.len());
        let sample_refs: Vec<&[u8]> = samples[..max_samples]
            .iter()
            .map(|v| v.as_slice())
            .collect();

        let dict_size = self.dict_size_kb as usize * 1024;
        match zstd::dict::from_continuous(
            &sample_refs.concat(),
            &sample_refs.iter().map(|s| s.len()).collect::<Vec<_>>(),
            dict_size,
        ) {
            Ok(dict) => {
                info!(
                    "zstd dict trained: {} samples, dict {} bytes",
                    max_samples,
                    dict.len()
                );
                store.store_dict(&dict)?;
                store.store_zstd_level(self.zstd_level)?;
                Ok(Some(Arc::new(dict)))
            }
            Err(e) => {
                info!("zstd dict training failed (falling back to uncompressed): {e}");
                Ok(None)
            }
        }
    }

    /// Build variation fjall from existing parquet files.
    ///
    /// Reads chr*.parquet and other.parquet in chrom_code order and feeds fjall.
    async fn build_variation_fjall_from_parquet(&self) -> Result<Vec<EntityStats>> {
        let parquet_dir = format!("{}/variation", self.output_dir);
        let fjall_dir = format!("{}/variation.fjall", self.output_dir);

        // Discover parquet files and sort in chrom_code order
        let mut parquet_files: Vec<(u16, String)> = Vec::new();
        for entry in std::fs::read_dir(&parquet_dir).map_err(|e| {
            DataFusionError::Execution(format!("Failed to read dir {parquet_dir}: {e}"))
        })? {
            let entry = entry.map_err(|e| {
                DataFusionError::Execution(format!("Failed to read dir entry: {e}"))
            })?;
            let path = entry.path();
            if path.extension().and_then(|e| e.to_str()) != Some("parquet") {
                continue;
            }
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or_default();
            // chr1.parquet → "1", other.parquet → sort last
            let chrom = stem.strip_prefix("chr").unwrap_or(stem);
            let code = if chrom == "other" {
                u16::MAX
            } else {
                chrom_to_code(chrom)
            };
            parquet_files.push((code, path.to_string_lossy().to_string()));
        }
        parquet_files.sort_by_key(|(code, _)| *code);

        if parquet_files.is_empty() {
            return Err(DataFusionError::Execution(
                "No parquet files found for variation fjall rebuild".to_string(),
            ));
        }

        info!(
            "variation: rebuilding fjall from {} parquet files",
            parquet_files.len()
        );

        // Create fjall store
        let read_ctx = SessionContext::new_with_config(
            SessionConfig::new().with_target_partitions(1),
        );
        // Probe schema from first file
        read_ctx
            .register_parquet("_probe", &parquet_files[0].1, Default::default())
            .await?;
        let probe_df = read_ctx.sql("SELECT * FROM _probe LIMIT 0").await?;
        let schema: SchemaRef = Arc::new(probe_df.schema().as_arrow().clone());
        read_ctx.deregister_table("_probe")?;

        let store = VepKvStore::create(Path::new(&fjall_dir), schema.clone())?;

        // Register non-canonical contigs from the parquet files
        let non_canonical: Vec<String> = parquet_files
            .iter()
            .filter_map(|(_, path)| {
                let stem = Path::new(path)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or_default();
                let chrom = stem.strip_prefix("chr").unwrap_or(stem);
                if chrom == "other" || CHROM_TO_CODE_SET.contains(&chrom) {
                    None // skip canonical and "other" (handled separately)
                } else {
                    Some(chrom.to_string())
                }
            })
            .collect();
        if !non_canonical.is_empty() {
            let refs: Vec<&str> = non_canonical.iter().map(|s| s.as_str()).collect();
            let mapping =
                crate::kv_cache::key_encoding::register_non_canonical_contigs(&refs);
            store.store_contig_codes(&mapping)?;
        }

        // Train dict from sample
        let dict = {
            read_ctx
                .register_parquet("_sample", &parquet_files[0].1, Default::default())
                .await?;
            let sample_df = read_ctx.sql("SELECT * FROM _sample LIMIT 10000").await?;
            let batches = sample_df.collect().await?;
            read_ctx.deregister_table("_sample")?;
            self.train_dict_from_batches(&store, &schema, &batches)
                .await?
        };

        let mut compressor = dict.as_ref().map(|d| {
            zstd::bulk::Compressor::with_dictionary(self.zstd_level, d)
                .expect("create zstd compressor")
        });

        // Start ingestion
        let mut ing = store
            .data_partition()
            .start_ingestion()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let chrom_col_idx = schema.index_of("chrom")?;
        let start_col_idx = schema.index_of("start")?;
        let allele_col_idx = schema.index_of("allele_string")?;
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_col_idx && i != start_col_idx)
            .collect();

        let mut accum: Option<PositionAccumulator> = None;
        let mut total_positions = 0u64;
        let mut total_variants = 0u64;
        let mut total_bytes = 0u64;
        let start_time = Instant::now();

        let mut time_parquet_read = std::time::Duration::ZERO;
        let mut time_serialize_compress = std::time::Duration::ZERO;
        let mut time_fjall_write = std::time::Duration::ZERO;

        for (_, parquet_path) in &parquet_files {
            read_ctx
                .register_parquet("_part", parquet_path, Default::default())
                .await?;
            let df = read_ctx.sql("SELECT * FROM _part").await?;
            let mut stream = df.execute_stream().await?;
            read_ctx.deregister_table("_part")?;

            while let Some(batch_result) = stream.next().await {
                let t0 = Instant::now();
                let batch = batch_result?;
                time_parquet_read += t0.elapsed();

                if batch.num_rows() == 0 {
                    continue;
                }
                let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
                let chrom_col = batch.column(chrom_col_idx);

                for row in 0..batch.num_rows() {
                    let start = starts.value(row);
                    let row_chrom = string_value(chrom_col.as_ref(), row);
                    let chrom_code = chrom_to_code(row_chrom);

                    let should_flush = accum
                        .as_ref()
                        .is_some_and(|a| a.chrom_code != chrom_code || a.start != start);

                    if should_flush {
                        let t1 = Instant::now();
                        let a = accum.take().unwrap();
                        let (key, value) =
                            a.finish_entry(&col_indices, allele_col_idx, &mut compressor)?;
                        time_serialize_compress += t1.elapsed();

                        let t2 = Instant::now();
                        ing.write(&key, &value)
                            .map_err(|e| DataFusionError::External(Box::new(e)))?;
                        time_fjall_write += t2.elapsed();

                        total_positions += 1;
                        total_bytes += value.len() as u64;
                    }

                    total_variants += 1;

                    match &mut accum {
                        Some(a) => a.add_row(row, &batch),
                        None => {
                            accum = Some(PositionAccumulator::new(
                                chrom_code,
                                row_chrom.to_string(),
                                start,
                                row,
                                &batch,
                            ));
                        }
                    }
                }
            }

            // Flush between files
            if let Some(a) = accum.take() {
                let (key, value) =
                    a.finish_entry(&col_indices, allele_col_idx, &mut compressor)?;
                ing.write(&key, &value)
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                total_positions += 1;
                total_bytes += value.len() as u64;
            }
        }

        ing.finish()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        store.persist()?;

        let total_accounted = time_parquet_read + time_serialize_compress + time_fjall_write;
        let total_wall = start_time.elapsed();
        let time_other = total_wall.saturating_sub(total_accounted);
        info!(
            "variation.fjall timing breakdown:\n\
             \x20 parquet read:         {:.1}s ({:.0}%)\n\
             \x20 serialize+compress:   {:.1}s ({:.0}%)\n\
             \x20 fjall write:          {:.1}s ({:.0}%)\n\
             \x20 other (row iteration):{:.1}s ({:.0}%)\n\
             \x20 total:                {:.1}s",
            time_parquet_read.as_secs_f64(),
            100.0 * time_parquet_read.as_secs_f64() / total_wall.as_secs_f64(),
            time_serialize_compress.as_secs_f64(),
            100.0 * time_serialize_compress.as_secs_f64() / total_wall.as_secs_f64(),
            time_fjall_write.as_secs_f64(),
            100.0 * time_fjall_write.as_secs_f64() / total_wall.as_secs_f64(),
            time_other.as_secs_f64(),
            100.0 * time_other.as_secs_f64() / total_wall.as_secs_f64(),
            total_wall.as_secs_f64(),
        );

        info!(
            "Running major compaction on variation.fjall...",
        );
        let compact_start = Instant::now();
        store
            .data_partition()
            .major_compact()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        info!(
            "Major compaction completed in {:.1}s",
            compact_start.elapsed().as_secs_f64()
        );

        let elapsed = start_time.elapsed().as_secs_f64();
        info!(
            "variation.fjall rebuilt: {} positions, {} variants, {:.1} MB in {:.1}s",
            total_positions,
            total_variants,
            total_bytes as f64 / 1_048_576.0,
            elapsed
        );

        Ok(vec![EntityStats {
            entity: "variation".to_string(),
            parquet_files: vec![],
            fjall_stats: Some(LoadStats {
                total_variants,
                total_positions,
                total_bytes,
                elapsed_secs: elapsed,
            }),
        }])
    }

    /// Train zstd dictionary from pre-collected batches.
    async fn train_dict_from_batches(
        &self,
        store: &VepKvStore,
        schema: &SchemaRef,
        batches: &[RecordBatch],
    ) -> Result<Option<Arc<Vec<u8>>>> {
        let chrom_col_idx = schema.index_of("chrom")?;
        let start_col_idx = schema.index_of("start")?;
        let allele_col_idx = schema.index_of("allele_string")?;
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_col_idx && i != start_col_idx)
            .collect();

        let mut samples: Vec<Vec<u8>> = Vec::new();
        for batch in batches {
            if batch.num_rows() == 0 {
                continue;
            }
            let chrom_col = batch.column(chrom_col_idx);
            let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();

            let mut groups: HashMap<(String, i64), Vec<usize>> = HashMap::new();
            for row in 0..batch.num_rows() {
                let chrom = string_value(chrom_col.as_ref(), row);
                let start = starts.value(row);
                groups
                    .entry((chrom.to_string(), start))
                    .or_default()
                    .push(row);
            }

            for rows in groups.values() {
                let value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
                samples.push(value);
            }
        }

        let min_training_samples = 100;
        if samples.len() < min_training_samples {
            info!(
                "only {} sample entries, skipping dictionary training (need >= {})",
                samples.len(),
                min_training_samples
            );
            return Ok(None);
        }

        let max_samples = 10_000.min(samples.len());
        let sample_refs: Vec<&[u8]> = samples[..max_samples]
            .iter()
            .map(|v| v.as_slice())
            .collect();

        let dict_size = self.dict_size_kb as usize * 1024;
        match zstd::dict::from_continuous(
            &sample_refs.concat(),
            &sample_refs.iter().map(|s| s.len()).collect::<Vec<_>>(),
            dict_size,
        ) {
            Ok(dict) => {
                info!(
                    "zstd dict trained: {} samples, dict {} bytes",
                    max_samples,
                    dict.len()
                );
                store.store_dict(&dict)?;
                store.store_zstd_level(self.zstd_level)?;
                Ok(Some(Arc::new(dict)))
            }
            Err(e) => {
                info!("zstd dict training failed (falling back to uncompressed): {e}");
                Ok(None)
            }
        }
    }

    /// Build translation: split into core + sift parquet, then fjall for sift.
    async fn build_translation(&self) -> Result<Vec<EntityStats>> {
        let table_name = "tl";
        let kind = EnsemblEntityKind::Translation;

        // Skip if all outputs exist
        if !self.overwrite {
            let core_exists =
                dir_has_parquet_files(&format!("{}/translation_core", self.output_dir));
            let sift_parquet_exists =
                dir_has_parquet_files(&format!("{}/translation_sift", self.output_dir));
            let sift_fjall_exists =
                Path::new(&format!("{}/translation_sift.fjall", self.output_dir)).exists();

            if core_exists && sift_parquet_exists && (sift_fjall_exists || !self.build_fjall) {
                info!("translation: all outputs exist, skipping (use overwrite to rebuild)");
                return Ok(vec![EntityStats {
                    entity: "translation".to_string(),
                    parquet_files: vec![],
                    fjall_stats: None,
                }]);
            }
        }

        // Discover chroms
        let init_ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
        let provider_schema = {
            let table = init_ctx.table(table_name).await?;
            table.schema().inner().clone()
        };
        let chroms = chroms_from_schema(&provider_schema);
        drop(init_ctx);

        let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
        let (main_chroms, other_chroms) = split_chroms(&chroms, &main_set);

        let mut core_results: Vec<(String, usize)> = Vec::new();
        let mut sift_results: Vec<(String, usize)> = Vec::new();

        // Process each main chromosome
        for chrom in &main_chroms {
            let (core_res, sift_res) = self
                .build_translation_chrom(chrom, kind, table_name)
                .await?;
            if let Some(r) = core_res {
                core_results.push(r);
            }
            if let Some(r) = sift_res {
                sift_results.push(r);
            }
        }

        // Process other contigs
        if !other_chroms.is_empty() {
            let (core_res, sift_res) = self
                .build_translation_multi_chrom(&other_chroms, kind, table_name)
                .await?;
            if let Some(r) = core_res {
                core_results.push(r);
            }
            if let Some(r) = sift_res {
                sift_results.push(r);
            }
        }

        // Build fjall for translation_sift (second pass, re-sorted by transcript_id)
        let sift_fjall_stats = if self.build_fjall && !sift_results.is_empty() {
            self.build_sift_fjall(&sift_results).await?
        } else {
            None
        };

        Ok(vec![
            EntityStats {
                entity: "translation_core".to_string(),
                parquet_files: core_results,
                fjall_stats: None,
            },
            EntityStats {
                entity: "translation_sift".to_string(),
                parquet_files: sift_results,
                fjall_stats: sift_fjall_stats,
            },
        ])
    }

    async fn build_translation_chrom(
        &self,
        chrom: &str,
        kind: EnsemblEntityKind,
        table_name: &str,
    ) -> Result<(Option<(String, usize)>, Option<(String, usize)>)> {
        let ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;

        let dedup_query = format!(
            "SELECT * FROM (\
                SELECT *, ROW_NUMBER() OVER (\
                    PARTITION BY transcript_id \
                    ORDER BY cdna_coding_start NULLS LAST\
                ) AS _rn \
                FROM {table_name} WHERE chrom = '{chrom}'\
            ) WHERE _rn = 1"
        );
        let df = ctx.sql(&dedup_query).await?;
        let schema = df.schema().clone();
        let cols: Vec<_> = schema
            .columns()
            .into_iter()
            .filter(|c| c.name() != "_rn")
            .collect();
        let df = df.select_columns(&cols.iter().map(|c| c.name()).collect::<Vec<_>>())?;
        let deduped = df.collect().await?;

        if deduped.is_empty() || deduped.iter().all(|b| b.num_rows() == 0) {
            return Ok((None, None));
        }

        let mem_table =
            datafusion::datasource::MemTable::try_new(deduped[0].schema(), vec![deduped])?;
        let split_ctx = SessionContext::new();
        split_ctx.register_table("_tl_deduped", Arc::new(mem_table))?;

        // translation_core
        let core_schema = datafusion_bio_format_ensembl_cache::translation_core_schema(false);
        let core_select = core_schema
            .fields()
            .iter()
            .map(|f| format!("\"{}\"", f.name()))
            .collect::<Vec<_>>()
            .join(", ");
        let core_file = format!("{}/translation_core/chr{chrom}.parquet", self.output_dir);
        let core_query = format!("SELECT {core_select} FROM _tl_deduped ORDER BY transcript_id");

        let mut w = create_writer(&core_file, &core_schema, kind, &["transcript_id"], None)?;
        let core_df = split_ctx.sql(&core_query).await?;
        let mut stream = core_df.execute_stream().await?;
        let mut core_rows = 0usize;
        while let Some(batch_result) = stream.next().await {
            let batch = batch_result?;
            if batch.num_rows() == 0 {
                continue;
            }
            let batch = project_batch(&batch, &core_schema)?;
            core_rows += batch.num_rows();
            w.write(&batch)?;

            if let Some(ref cb) = self.on_progress {
                cb(
                    "translation_core",
                    "parquet",
                    batch.num_rows(),
                    core_rows,
                    0,
                );
            }
        }
        w.close().map_err(|e| {
            DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
        })?;

        // translation_sift
        let sift_schema = datafusion_bio_format_ensembl_cache::translation_sift_schema(false);
        let sift_select = sift_schema
            .fields()
            .iter()
            .map(|f| format!("\"{}\"", f.name()))
            .collect::<Vec<_>>()
            .join(", ");
        let sift_file = format!("{}/translation_sift/chr{chrom}.parquet", self.output_dir);
        let sift_query = format!("SELECT {sift_select} FROM _tl_deduped ORDER BY chrom, start");

        let mut w = create_writer(
            &sift_file,
            &sift_schema,
            kind,
            &["chrom", "start"],
            Some(256),
        )?;
        let sift_df = split_ctx.sql(&sift_query).await?;
        let mut stream = sift_df.execute_stream().await?;
        let mut sift_rows = 0usize;
        while let Some(batch_result) = stream.next().await {
            let batch = batch_result?;
            if batch.num_rows() == 0 {
                continue;
            }
            let batch = project_batch(&batch, &sift_schema)?;
            sift_rows += batch.num_rows();
            w.write(&batch)?;

            if let Some(ref cb) = self.on_progress {
                cb(
                    "translation_sift",
                    "parquet",
                    batch.num_rows(),
                    sift_rows,
                    0,
                );
            }
        }
        w.close().map_err(|e| {
            DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
        })?;

        info!(
            "translation: chr{chrom} core={} sift={}",
            format_rows(core_rows),
            format_rows(sift_rows)
        );

        let core_result = if core_rows > 0 {
            Some((core_file, core_rows))
        } else {
            None
        };
        let sift_result = if sift_rows > 0 {
            Some((sift_file, sift_rows))
        } else {
            None
        };
        Ok((core_result, sift_result))
    }

    async fn build_translation_multi_chrom(
        &self,
        other_chroms: &[String],
        kind: EnsemblEntityKind,
        table_name: &str,
    ) -> Result<(Option<(String, usize)>, Option<(String, usize)>)> {
        let ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
        let in_list = other_chroms
            .iter()
            .map(|c| format!("'{c}'"))
            .collect::<Vec<_>>()
            .join(", ");
        let dedup_query = format!(
            "SELECT * FROM (\
                SELECT *, ROW_NUMBER() OVER (\
                    PARTITION BY transcript_id \
                    ORDER BY cdna_coding_start NULLS LAST\
                ) AS _rn \
                FROM {table_name} WHERE chrom IN ({in_list})\
            ) WHERE _rn = 1"
        );
        let df = ctx.sql(&dedup_query).await?;
        let schema = df.schema().clone();
        let cols: Vec<_> = schema
            .columns()
            .into_iter()
            .filter(|c| c.name() != "_rn")
            .collect();
        let df = df.select_columns(&cols.iter().map(|c| c.name()).collect::<Vec<_>>())?;
        let deduped = df.collect().await?;

        if deduped.is_empty() || deduped.iter().all(|b| b.num_rows() == 0) {
            return Ok((None, None));
        }

        let mem_table =
            datafusion::datasource::MemTable::try_new(deduped[0].schema(), vec![deduped])?;
        let split_ctx = SessionContext::new();
        split_ctx.register_table("_tl_deduped", Arc::new(mem_table))?;

        let core_schema = datafusion_bio_format_ensembl_cache::translation_core_schema(false);
        let core_select = core_schema
            .fields()
            .iter()
            .map(|f| format!("\"{}\"", f.name()))
            .collect::<Vec<_>>()
            .join(", ");
        let core_file = format!("{}/translation_core/other.parquet", self.output_dir);
        let mut w = create_writer(&core_file, &core_schema, kind, &["transcript_id"], None)?;
        let core_rows = stream_to_writer_with_progress(
            &split_ctx,
            &format!("SELECT {core_select} FROM _tl_deduped ORDER BY transcript_id"),
            &mut w,
            false,
            self.on_progress.as_deref(),
            "translation_core",
        )
        .await?;
        w.close().map_err(|e| {
            DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
        })?;

        let sift_schema = datafusion_bio_format_ensembl_cache::translation_sift_schema(false);
        let sift_select = sift_schema
            .fields()
            .iter()
            .map(|f| format!("\"{}\"", f.name()))
            .collect::<Vec<_>>()
            .join(", ");
        let sift_file = format!("{}/translation_sift/other.parquet", self.output_dir);
        let mut w = create_writer(
            &sift_file,
            &sift_schema,
            kind,
            &["chrom", "start"],
            Some(256),
        )?;
        let sift_rows = stream_to_writer_with_progress(
            &split_ctx,
            &format!("SELECT {sift_select} FROM _tl_deduped ORDER BY chrom, start"),
            &mut w,
            false,
            self.on_progress.as_deref(),
            "translation_sift",
        )
        .await?;
        w.close().map_err(|e| {
            DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
        })?;

        info!(
            "translation: other ({} contigs) core={} sift={}",
            other_chroms.len(),
            format_rows(core_rows),
            format_rows(sift_rows)
        );

        let core_result = if core_rows > 0 {
            Some((core_file, core_rows))
        } else {
            None
        };
        let sift_result = if sift_rows > 0 {
            Some((sift_file, sift_rows))
        } else {
            None
        };
        Ok((core_result, sift_result))
    }

    /// Build fjall for translation_sift: re-read parquet sorted by transcript_id,
    /// parse predictions, and ingest into SiftKvStore.
    async fn build_sift_fjall(
        &self,
        sift_parquet_files: &[(String, usize)],
    ) -> Result<Option<LoadStats>> {
        let fjall_dir = format!("{}/translation_sift.fjall", self.output_dir);
        if Path::new(&fjall_dir).exists() && !self.overwrite {
            info!("translation_sift.fjall already exists, skipping (use overwrite to rebuild)");
            return Ok(None);
        }

        let start_time = Instant::now();

        // Register all sift parquet files as a single table
        let ctx = SessionContext::new();
        let parquet_paths: Vec<&str> = sift_parquet_files.iter().map(|(p, _)| p.as_str()).collect();

        // Register each parquet file and UNION ALL them
        for (i, path) in parquet_paths.iter().enumerate() {
            ctx.register_parquet(&format!("_sift_{i}"), path, Default::default())
                .await?;
        }

        let union_query = if parquet_paths.len() == 1 {
            "SELECT * FROM _sift_0 ORDER BY transcript_id".to_string()
        } else {
            let unions: Vec<String> = (0..parquet_paths.len())
                .map(|i| format!("SELECT * FROM _sift_{i}"))
                .collect();
            format!(
                "SELECT * FROM ({}) ORDER BY transcript_id",
                unions.join(" UNION ALL ")
            )
        };

        let df = ctx.sql(&union_query).await?;
        let mut stream = df.execute_stream().await?;

        // Open fjall database
        let db = fjall::Database::builder(&fjall_dir)
            .cache_size(64 * 1024 * 1024)
            .open()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        let sift_store = SiftKvStore::create(&db)?;

        // Stream predictions directly into fjall, grouping by transcript_id
        // on-the-fly (data arrives sorted). Only one transcript is buffered.
        let mut current_transcript: Option<String> = None;
        let mut current_preds = CachedPredictions::default();
        let mut total_rows = 0usize;
        let mut transcript_count = 0usize;

        while let Some(batch_result) = stream.next().await {
            let batch = batch_result?;
            if batch.num_rows() == 0 {
                continue;
            }

            let schema = batch.schema();
            let tid_idx = schema.index_of("transcript_id")?;
            let sift_idx = schema.index_of("sift_predictions").ok();
            let poly_idx = schema.index_of("polyphen_predictions").ok();

            for row in 0..batch.num_rows() {
                let transcript_id = string_value(batch.column(tid_idx).as_ref(), row).to_string();

                let mut row_preds = CachedPredictions::default();
                if let Some(idx) = sift_idx {
                    row_preds.sift = read_compact_predictions(batch.column(idx).as_ref(), row);
                }
                if let Some(idx) = poly_idx {
                    row_preds.polyphen = read_compact_predictions(batch.column(idx).as_ref(), row);
                }

                if current_transcript.as_deref() != Some(&transcript_id) {
                    // Flush previous transcript directly to fjall
                    if let Some(tid) = current_transcript.take() {
                        current_preds.sort();
                        sift_store.put(&tid, &current_preds)?;
                        current_preds = CachedPredictions::default();
                        transcript_count += 1;
                    }
                    current_transcript = Some(transcript_id);
                }
                current_preds.sift.extend(row_preds.sift);
                current_preds.polyphen.extend(row_preds.polyphen);
                total_rows += 1;
            }

            if let Some(ref cb) = self.on_progress {
                cb("translation_sift", "fjall", batch.num_rows(), total_rows, 0);
            }
        }

        // Flush last transcript
        if let Some(tid) = current_transcript.take() {
            current_preds.sort();
            sift_store.put(&tid, &current_preds)?;
            transcript_count += 1;
        }

        db.persist(fjall::PersistMode::SyncAll)
            .map_err(|e| DataFusionError::External(Box::new(e)))?;

        let elapsed = start_time.elapsed().as_secs_f64();
        info!(
            "translation_sift.fjall: {} transcripts from {} rows in {:.1}s",
            transcript_count, total_rows, elapsed
        );

        Ok(Some(LoadStats {
            total_variants: total_rows as u64,
            total_positions: transcript_count as u64,
            total_bytes: 0, // fjall manages its own storage
            elapsed_secs: elapsed,
        }))
    }

    /// Build a parquet-only entity (transcript, exon, regulatory, motif).
    async fn build_parquet_entity(&self, kind: EnsemblEntityKind) -> Result<Vec<(String, usize)>> {
        let table_name = entity_table_name(kind);
        let subdir = entity_subdir(kind);
        let needs_rn_drop = matches!(
            kind,
            EnsemblEntityKind::Transcript | EnsemblEntityKind::Exon
        );

        // Discover chroms
        let init_ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
        let provider_schema = {
            let table = init_ctx.table(table_name).await?;
            table.schema().inner().clone()
        };
        let chroms = chroms_from_schema(&provider_schema);
        drop(init_ctx);

        let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
        let (main_chroms, other_chroms) = split_chroms(&chroms, &main_set);

        info!(
            "{subdir}: {} main chroms, {} other contigs",
            main_chroms.len(),
            other_chroms.len()
        );

        let mut all_results = Vec::new();
        let global_start = Instant::now();
        let mut total_rows: usize = 0;

        for chrom in &main_chroms {
            let ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
            let query = build_query(kind, table_name, Some(chrom));
            let output_file = format!("{}/{subdir}/chr{chrom}.parquet", self.output_dir);

            let df = ctx.sql(&query).await?;
            let df = if needs_rn_drop {
                let schema = df.schema().clone();
                let cols: Vec<_> = schema
                    .columns()
                    .into_iter()
                    .filter(|c| c.name() != "_rn")
                    .collect();
                df.select_columns(&cols.iter().map(|c| c.name()).collect::<Vec<_>>())?
            } else {
                df
            };
            let mut stream = df.execute_stream().await?;
            let schema = stream.schema();
            let sk = sort_key(kind);
            let mut writer = create_writer(&output_file, &schema, kind, sk, None)?;

            let mut chrom_rows = 0usize;
            while let Some(batch_result) = stream.next().await {
                let batch = batch_result?;
                if batch.num_rows() == 0 {
                    continue;
                }
                chrom_rows += batch.num_rows();
                total_rows += batch.num_rows();
                writer.write(&batch)?;

                if let Some(ref cb) = self.on_progress {
                    cb(subdir, "parquet", batch.num_rows(), total_rows, 0);
                }
            }
            writer.close().map_err(|e| {
                DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
            })?;

            if chrom_rows > 0 {
                all_results.push((output_file, chrom_rows));
            }
        }

        // Process remaining contigs as "other.parquet"
        if !other_chroms.is_empty() {
            let ctx = make_ctx_and_register(&self.cache_root, kind, table_name, self.partitions)?;
            let other_refs: Vec<&str> = other_chroms.iter().map(|s| s.as_str()).collect();
            let query = build_query_multi_chrom(kind, table_name, &other_refs);
            let output_file = format!("{}/{subdir}/other.parquet", self.output_dir);

            let df = ctx.sql(&query).await?;
            let df = if needs_rn_drop {
                let schema = df.schema().clone();
                let cols: Vec<_> = schema
                    .columns()
                    .into_iter()
                    .filter(|c| c.name() != "_rn")
                    .collect();
                df.select_columns(&cols.iter().map(|c| c.name()).collect::<Vec<_>>())?
            } else {
                df
            };
            let mut stream = df.execute_stream().await?;
            let schema = stream.schema();
            let sk = sort_key(kind);
            let mut writer = create_writer(&output_file, &schema, kind, sk, None)?;

            let mut other_rows = 0usize;
            while let Some(batch_result) = stream.next().await {
                let batch = batch_result?;
                if batch.num_rows() == 0 {
                    continue;
                }
                other_rows += batch.num_rows();
                total_rows += batch.num_rows();
                writer.write(&batch)?;

                if let Some(ref cb) = self.on_progress {
                    cb(subdir, "parquet", batch.num_rows(), total_rows, 0);
                }
            }
            writer.close().map_err(|e| {
                DataFusionError::Execution(format!("Failed to close parquet writer: {e}"))
            })?;

            if other_rows > 0 {
                all_results.push((output_file, other_rows));
            }
        }

        let elapsed = global_start.elapsed().as_secs_f64();
        print_progress(subdir, total_rows, elapsed);
        Ok(all_results)
    }
}

// ---------------------------------------------------------------------------
// Position accumulator for fjall dual-sink
// ---------------------------------------------------------------------------

/// Accumulates rows for a single genomic position during sorted streaming.
/// Handles rows spanning multiple record batches by storing (batch, row_indices) pairs.
struct PositionAccumulator {
    chrom_code: u16,
    chrom: String,
    start: i64,
    /// Each entry is a (batch, row_indices) pair. A new entry is added when
    /// the stream advances to a new RecordBatch while the position key stays the same.
    segments: Vec<(RecordBatch, Vec<usize>)>,
}

impl PositionAccumulator {
    fn new(chrom_code: u16, chrom: String, start: i64, row: usize, batch: &RecordBatch) -> Self {
        Self {
            chrom_code,
            chrom,
            start,
            segments: vec![(batch.clone(), vec![row])],
        }
    }

    fn add_row(&mut self, row: usize, batch: &RecordBatch) {
        // Always push a new segment per batch reference. The multi-segment
        // path in `finish_entry` handles merging via take + concat_batches.
        // This avoids an unreliable same-batch heuristic (DataFusion reuses
        // the same SchemaRef Arc across all batches in a stream, so pointer
        // comparison doesn't distinguish distinct batches).
        self.segments.push((batch.clone(), vec![row]));
    }

    fn finish_entry(
        self,
        col_indices: &[usize],
        allele_col_idx: usize,
        compressor: &mut Option<zstd::bulk::Compressor<'_>>,
    ) -> Result<(Vec<u8>, Vec<u8>)> {
        let key = encode_position_key(&self.chrom, self.start);

        // If all rows are in one batch (common case), serialize directly.
        // Otherwise, merge segments by slicing each batch to its rows and concatenating.
        let raw_value = if self.segments.len() == 1 {
            let (batch, rows) = &self.segments[0];
            serialize_position_entry(rows, batch, col_indices, allele_col_idx)?
        } else {
            // Collect all rows across batch boundaries into a single batch
            let mut all_batches: Vec<RecordBatch> = Vec::new();
            for (batch, rows) in &self.segments {
                let indices: Vec<u32> = rows.iter().map(|&r| r as u32).collect();
                let idx_array = datafusion::arrow::array::UInt32Array::from(indices);
                let taken = datafusion::arrow::compute::take_record_batch(batch, &idx_array)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                all_batches.push(taken);
            }
            let merged = datafusion::arrow::compute::concat_batches(
                &all_batches[0].schema(),
                all_batches.iter(),
            )
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            let row_indices: Vec<usize> = (0..merged.num_rows()).collect();
            serialize_position_entry(&row_indices, &merged, col_indices, allele_col_idx)?
        };

        let value = if let Some(comp) = compressor {
            comp.compress(&raw_value)
                .map_err(|e| DataFusionError::Execution(format!("zstd compression failed: {e}")))?
        } else {
            raw_value
        };

        Ok((key, value))
    }
}

// ---------------------------------------------------------------------------
// Fjall state for variation dual-sink
// ---------------------------------------------------------------------------

struct FjallVariationState {
    store: VepKvStore,
    schema: SchemaRef,
    dict: Option<Arc<Vec<u8>>>,
    #[allow(dead_code)]
    fjall_dir: String,
}

// ---------------------------------------------------------------------------
// Helper functions (ported from vepyr convert.rs)
// ---------------------------------------------------------------------------

fn parse_entity(name: &str) -> Option<EnsemblEntityKind> {
    match name {
        "variation" => Some(EnsemblEntityKind::Variation),
        "transcript" => Some(EnsemblEntityKind::Transcript),
        "exon" => Some(EnsemblEntityKind::Exon),
        "translation" => Some(EnsemblEntityKind::Translation),
        "regulatory" => Some(EnsemblEntityKind::RegulatoryFeature),
        "motif" => Some(EnsemblEntityKind::MotifFeature),
        _ => None,
    }
}

fn entity_subdir(kind: EnsemblEntityKind) -> &'static str {
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

fn row_group_size(kind: EnsemblEntityKind) -> usize {
    match kind {
        EnsemblEntityKind::Variation => 100_000,
        EnsemblEntityKind::Transcript => 8_000,
        EnsemblEntityKind::Exon => 45_000,
        EnsemblEntityKind::Translation => 6_000,
        EnsemblEntityKind::RegulatoryFeature => 9_000,
        EnsemblEntityKind::MotifFeature => 10_000,
    }
}

fn sort_key(kind: EnsemblEntityKind) -> &'static [&'static str] {
    match kind {
        EnsemblEntityKind::Exon => &["transcript_id", "start"],
        _ => &["chrom", "start"],
    }
}

fn sorting_columns_for(schema: &SchemaRef, sort_columns: &[&str]) -> Option<Vec<SortingColumn>> {
    let cols: Vec<SortingColumn> = sort_columns
        .iter()
        .filter_map(|name| {
            schema
                .column_with_name(name)
                .map(|(idx, _)| SortingColumn::new(idx as i32, false, false))
        })
        .collect();
    if cols.len() == sort_columns.len() {
        Some(cols)
    } else {
        None
    }
}

fn writer_properties(
    kind: EnsemblEntityKind,
    schema: &SchemaRef,
    sort_columns: &[&str],
    rg_size_override: Option<usize>,
) -> WriterProperties {
    let rg_size = rg_size_override.unwrap_or_else(|| row_group_size(kind));
    let sorting = sorting_columns_for(schema, sort_columns);

    let mut builder = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(rg_size)
        .set_sorting_columns(sorting);

    if matches!(
        kind,
        EnsemblEntityKind::Translation | EnsemblEntityKind::Exon
    ) {
        builder = builder.set_column_bloom_filter_enabled(ColumnPath::from("transcript_id"), true);
    }

    builder.build()
}

fn build_query(kind: EnsemblEntityKind, table_name: &str, chrom_filter: Option<&str>) -> String {
    let where_clause = chrom_filter
        .map(|c| format!(" WHERE chrom = '{c}'"))
        .unwrap_or_default();

    match kind {
        EnsemblEntityKind::Transcript => {
            format!(
                "SELECT * FROM (\
                    SELECT *, ROW_NUMBER() OVER (\
                        PARTITION BY stable_id \
                        ORDER BY cds_start NULLS LAST\
                    ) AS _rn \
                    FROM {table_name}{where_clause}\
                ) WHERE _rn = 1 \
                ORDER BY chrom, start"
            )
        }
        EnsemblEntityKind::Translation => unreachable!("use build_translation() instead"),
        EnsemblEntityKind::Exon => {
            format!(
                "SELECT * FROM (\
                    SELECT *, ROW_NUMBER() OVER (\
                        PARTITION BY transcript_id, exon_number \
                        ORDER BY stable_id NULLS LAST\
                    ) AS _rn \
                    FROM {table_name}{where_clause}\
                ) WHERE _rn = 1 \
                ORDER BY transcript_id, start"
            )
        }
        _ => {
            format!("SELECT * FROM {table_name}{where_clause} ORDER BY chrom, start")
        }
    }
}

fn build_query_multi_chrom(kind: EnsemblEntityKind, table_name: &str, chroms: &[&str]) -> String {
    let list = chroms
        .iter()
        .map(|c| format!("'{c}'"))
        .collect::<Vec<_>>()
        .join(", ");
    let where_clause = format!(" WHERE chrom IN ({list})");

    match kind {
        EnsemblEntityKind::Transcript => {
            format!(
                "SELECT * FROM (\
                    SELECT *, ROW_NUMBER() OVER (\
                        PARTITION BY stable_id \
                        ORDER BY cds_start NULLS LAST\
                    ) AS _rn \
                    FROM {table_name}{where_clause}\
                ) WHERE _rn = 1 \
                ORDER BY chrom, start"
            )
        }
        EnsemblEntityKind::Translation => unreachable!(),
        EnsemblEntityKind::Exon => {
            format!(
                "SELECT * FROM (\
                    SELECT *, ROW_NUMBER() OVER (\
                        PARTITION BY transcript_id, exon_number \
                        ORDER BY stable_id NULLS LAST\
                    ) AS _rn \
                    FROM {table_name}{where_clause}\
                ) WHERE _rn = 1 \
                ORDER BY transcript_id, start"
            )
        }
        _ => {
            format!("SELECT * FROM {table_name}{where_clause} ORDER BY chrom, start")
        }
    }
}

fn project_batch(batch: &RecordBatch, target_schema: &SchemaRef) -> Result<RecordBatch> {
    let source_schema = batch.schema();
    let mut columns = Vec::with_capacity(target_schema.fields().len());
    for field in target_schema.fields() {
        let (idx, _) = source_schema
            .column_with_name(field.name())
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "Column '{}' not found in source batch",
                    field.name()
                ))
            })?;
        columns.push(batch.column(idx).clone());
    }
    RecordBatch::try_new(target_schema.clone(), columns)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}

fn chroms_from_schema(schema: &SchemaRef) -> Option<Vec<String>> {
    let result = schema
        .metadata()
        .get("bio.vep.chromosomes")
        .and_then(|json| serde_json::from_str(json).ok());
    if result.is_none() {
        info!(
            "Schema metadata key 'bio.vep.chromosomes' not found; \
             falling back to default MAIN_CHROMS list"
        );
    }
    result
}

fn format_rows(n: usize) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}k", n as f64 / 1_000.0)
    } else {
        format!("{n}")
    }
}

fn print_progress(label: &str, rows: usize, elapsed: f64) {
    let rate = if elapsed > 0.0 {
        format!("{}/s", format_rows((rows as f64 / elapsed) as usize))
    } else {
        "? rows/s".to_string()
    };
    eprintln!(
        "  {label}: {} rows [{:.1}s, {rate}]",
        format_rows(rows),
        elapsed
    );
}

/// If the plan is a `SortPreservingMergeExec` wrapping a multi-partition inner
/// plan, return the inner plan.  This allows bypassing the merge and executing
/// partitions in parallel with larger buffers.
fn extract_multi_partition_inner(
    plan: &Arc<dyn ExecutionPlan>,
) -> Option<Arc<dyn ExecutionPlan>> {
    if plan.name() != "SortPreservingMergeExec" {
        return None;
    }
    let children = plan.children();
    if children.len() != 1 {
        return None;
    }
    let inner = Arc::clone(children[0]);
    if inner.properties().partitioning.partition_count() > 1 {
        Some(inner)
    } else {
        None
    }
}

/// Set of canonical chromosome names for quick lookup.
const CHROM_TO_CODE_SET: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// Check if a directory contains at least one `.parquet` file.
fn dir_has_parquet_files(dir: &str) -> bool {
    Path::new(dir).is_dir()
        && std::fs::read_dir(dir)
            .ok()
            .map(|entries| {
                entries.flatten().any(|e| {
                    e.path()
                        .extension()
                        .and_then(|ext| ext.to_str())
                        == Some("parquet")
                })
            })
            .unwrap_or(false)
}

fn make_ctx_and_register(
    cache_root: &str,
    kind: EnsemblEntityKind,
    table_name: &str,
    partitions: usize,
) -> Result<SessionContext> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);
    let mut options = EnsemblCacheOptions::new(cache_root);
    options.target_partitions = Some(partitions);
    let provider = EnsemblCacheTableProvider::for_entity(kind, options)?;
    ctx.register_table(table_name, provider)?;
    Ok(ctx)
}

fn create_writer(
    path: &str,
    schema: &SchemaRef,
    kind: EnsemblEntityKind,
    sort_cols: &[&str],
    rg_override: Option<usize>,
) -> Result<ArrowWriter<File>> {
    let props = writer_properties(kind, schema, sort_cols, rg_override);
    let file = File::create(path)
        .map_err(|e| DataFusionError::Execution(format!("Failed to create file {path}: {e}")))?;
    ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| DataFusionError::Execution(format!("Failed to create ArrowWriter: {e}")))
}

async fn stream_to_writer_with_progress(
    ctx: &SessionContext,
    query: &str,
    writer: &mut ArrowWriter<File>,
    needs_rn_drop: bool,
    on_progress: Option<&OnProgress>,
    entity: &str,
) -> Result<usize> {
    let df = ctx.sql(query).await?;
    let df = if needs_rn_drop {
        let schema = df.schema().clone();
        let cols: Vec<_> = schema
            .columns()
            .into_iter()
            .filter(|c| c.name() != "_rn")
            .collect();
        df.select_columns(&cols.iter().map(|c| c.name()).collect::<Vec<_>>())?
    } else {
        df
    };
    let mut stream = df.execute_stream().await?;
    let mut rows = 0usize;
    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }
        rows += batch.num_rows();
        writer.write(&batch)?;

        if let Some(cb) = on_progress {
            cb(entity, "parquet", batch.num_rows(), rows, 0);
        }
    }
    Ok(rows)
}

fn split_chroms(
    chroms: &Option<Vec<String>>,
    main_set: &HashSet<&str>,
) -> (Vec<String>, Vec<String>) {
    match chroms {
        Some(all) => {
            let main: Vec<String> = all
                .iter()
                .filter(|c| main_set.contains(c.as_str()))
                .cloned()
                .collect();
            let other: Vec<String> = all
                .iter()
                .filter(|c| !main_set.contains(c.as_str()))
                .cloned()
                .collect();
            (main, other)
        }
        None => (MAIN_CHROMS.iter().map(|s| s.to_string()).collect(), vec![]),
    }
}

/// Extract a string value from a Utf8, LargeUtf8, or Utf8View array.
fn string_value(col: &dyn Array, row: usize) -> &str {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        arr.value(row)
    } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        arr.value(row)
    } else if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
        arr.value(row)
    } else {
        panic!("Expected Utf8, LargeUtf8, or Utf8View column")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::collections::HashMap;
    use std::sync::Arc;

    // -----------------------------------------------------------------------
    // Helper: create a test schema and RecordBatch for variation-like data
    // -----------------------------------------------------------------------
    fn variation_schema() -> SchemaRef {
        Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("variation_name", DataType::Utf8, true),
        ]))
    }

    fn make_batch(
        chroms: Vec<&str>,
        starts: Vec<i64>,
        ends: Vec<i64>,
        alleles: Vec<&str>,
        names: Vec<&str>,
    ) -> RecordBatch {
        RecordBatch::try_new(
            variation_schema(),
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts)),
                Arc::new(Int64Array::from(ends)),
                Arc::new(StringArray::from(alleles)),
                Arc::new(StringArray::from(names)),
            ],
        )
        .unwrap()
    }

    // -----------------------------------------------------------------------
    // parse_entity
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_entity_known() {
        assert!(matches!(
            parse_entity("variation"),
            Some(EnsemblEntityKind::Variation)
        ));
        assert!(matches!(
            parse_entity("transcript"),
            Some(EnsemblEntityKind::Transcript)
        ));
        assert!(matches!(
            parse_entity("exon"),
            Some(EnsemblEntityKind::Exon)
        ));
        assert!(matches!(
            parse_entity("translation"),
            Some(EnsemblEntityKind::Translation)
        ));
        assert!(matches!(
            parse_entity("regulatory"),
            Some(EnsemblEntityKind::RegulatoryFeature)
        ));
        assert!(matches!(
            parse_entity("motif"),
            Some(EnsemblEntityKind::MotifFeature)
        ));
    }

    #[test]
    fn test_parse_entity_unknown() {
        assert!(parse_entity("unknown").is_none());
        assert!(parse_entity("").is_none());
        assert!(parse_entity("Variation").is_none()); // case-sensitive
    }

    // -----------------------------------------------------------------------
    // entity_subdir
    // -----------------------------------------------------------------------
    #[test]
    fn test_entity_subdir() {
        assert_eq!(entity_subdir(EnsemblEntityKind::Variation), "variation");
        assert_eq!(entity_subdir(EnsemblEntityKind::Transcript), "transcript");
        assert_eq!(entity_subdir(EnsemblEntityKind::Exon), "exon");
        assert_eq!(entity_subdir(EnsemblEntityKind::Translation), "translation");
        assert_eq!(
            entity_subdir(EnsemblEntityKind::RegulatoryFeature),
            "regulatory"
        );
        assert_eq!(entity_subdir(EnsemblEntityKind::MotifFeature), "motif");
    }

    // -----------------------------------------------------------------------
    // entity_table_name
    // -----------------------------------------------------------------------
    #[test]
    fn test_entity_table_name() {
        assert_eq!(entity_table_name(EnsemblEntityKind::Variation), "var");
        assert_eq!(entity_table_name(EnsemblEntityKind::Transcript), "tx");
        assert_eq!(entity_table_name(EnsemblEntityKind::Exon), "exon");
        assert_eq!(entity_table_name(EnsemblEntityKind::Translation), "tl");
        assert_eq!(
            entity_table_name(EnsemblEntityKind::RegulatoryFeature),
            "reg"
        );
        assert_eq!(entity_table_name(EnsemblEntityKind::MotifFeature), "motif");
    }

    // -----------------------------------------------------------------------
    // row_group_size
    // -----------------------------------------------------------------------
    #[test]
    fn test_row_group_size() {
        assert_eq!(row_group_size(EnsemblEntityKind::Variation), 100_000);
        assert_eq!(row_group_size(EnsemblEntityKind::Transcript), 8_000);
        assert_eq!(row_group_size(EnsemblEntityKind::Exon), 45_000);
        assert_eq!(row_group_size(EnsemblEntityKind::Translation), 6_000);
        assert_eq!(row_group_size(EnsemblEntityKind::RegulatoryFeature), 9_000);
        assert_eq!(row_group_size(EnsemblEntityKind::MotifFeature), 10_000);
    }

    // -----------------------------------------------------------------------
    // sort_key
    // -----------------------------------------------------------------------
    #[test]
    fn test_sort_key() {
        assert_eq!(
            sort_key(EnsemblEntityKind::Exon),
            &["transcript_id", "start"]
        );
        assert_eq!(sort_key(EnsemblEntityKind::Variation), &["chrom", "start"]);
        assert_eq!(sort_key(EnsemblEntityKind::Transcript), &["chrom", "start"]);
        assert_eq!(
            sort_key(EnsemblEntityKind::RegulatoryFeature),
            &["chrom", "start"]
        );
    }

    // -----------------------------------------------------------------------
    // format_rows
    // -----------------------------------------------------------------------
    #[test]
    fn test_format_rows() {
        assert_eq!(format_rows(0), "0");
        assert_eq!(format_rows(999), "999");
        assert_eq!(format_rows(1_000), "1.0k");
        assert_eq!(format_rows(1_500), "1.5k");
        assert_eq!(format_rows(999_999), "1000.0k");
        assert_eq!(format_rows(1_000_000), "1.0M");
        assert_eq!(format_rows(1_170_000_000), "1170.0M");
    }

    // -----------------------------------------------------------------------
    // split_chroms
    // -----------------------------------------------------------------------
    #[test]
    fn test_split_chroms_with_known_chroms() {
        let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
        let chroms = Some(vec![
            "1".to_string(),
            "2".to_string(),
            "MT".to_string(),
            "GL000220.1".to_string(),
        ]);
        let (main, other) = split_chroms(&chroms, &main_set);
        assert_eq!(main, vec!["1", "2"]);
        assert_eq!(other, vec!["MT", "GL000220.1"]);
    }

    #[test]
    fn test_split_chroms_none_defaults_to_main() {
        let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
        let (main, other) = split_chroms(&None, &main_set);
        assert_eq!(main.len(), MAIN_CHROMS.len());
        assert!(other.is_empty());
    }

    #[test]
    fn test_split_chroms_all_non_canonical() {
        let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
        let chroms = Some(vec!["KI270442.1".to_string(), "GL000220.1".to_string()]);
        let (main, other) = split_chroms(&chroms, &main_set);
        assert!(main.is_empty());
        assert_eq!(other.len(), 2);
    }

    // -----------------------------------------------------------------------
    // chroms_from_schema
    // -----------------------------------------------------------------------
    #[test]
    fn test_chroms_from_schema_present() {
        let mut metadata = HashMap::new();
        metadata.insert(
            "bio.vep.chromosomes".to_string(),
            r#"["1","2","X"]"#.to_string(),
        );
        let schema = Arc::new(Schema::new_with_metadata(
            vec![Field::new("chrom", DataType::Utf8, false)],
            metadata,
        ));
        let chroms = chroms_from_schema(&schema);
        assert_eq!(
            chroms,
            Some(vec!["1".to_string(), "2".to_string(), "X".to_string()])
        );
    }

    #[test]
    fn test_chroms_from_schema_absent() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "chrom",
            DataType::Utf8,
            false,
        )]));
        assert!(chroms_from_schema(&schema).is_none());
    }

    #[test]
    fn test_chroms_from_schema_malformed_json() {
        let mut metadata = HashMap::new();
        metadata.insert("bio.vep.chromosomes".to_string(), "not json".to_string());
        let schema = Arc::new(Schema::new_with_metadata(
            vec![Field::new("chrom", DataType::Utf8, false)],
            metadata,
        ));
        assert!(chroms_from_schema(&schema).is_none());
    }

    // -----------------------------------------------------------------------
    // build_query
    // -----------------------------------------------------------------------
    #[test]
    fn test_build_query_variation_no_filter() {
        let q = build_query(EnsemblEntityKind::Variation, "var", None);
        assert_eq!(q, "SELECT * FROM var ORDER BY chrom, start");
    }

    #[test]
    fn test_build_query_variation_with_filter() {
        let q = build_query(EnsemblEntityKind::Variation, "var", Some("1"));
        assert_eq!(
            q,
            "SELECT * FROM var WHERE chrom = '1' ORDER BY chrom, start"
        );
    }

    #[test]
    fn test_build_query_transcript_dedup() {
        let q = build_query(EnsemblEntityKind::Transcript, "tx", Some("X"));
        assert!(q.contains("ROW_NUMBER()"));
        assert!(q.contains("PARTITION BY stable_id"));
        assert!(q.contains("WHERE _rn = 1"));
        assert!(q.contains("ORDER BY chrom, start"));
        assert!(q.contains("WHERE chrom = 'X'"));
    }

    #[test]
    fn test_build_query_exon_dedup() {
        let q = build_query(EnsemblEntityKind::Exon, "exon", None);
        assert!(q.contains("PARTITION BY transcript_id, exon_number"));
        assert!(q.contains("ORDER BY transcript_id, start"));
    }

    // -----------------------------------------------------------------------
    // build_query_multi_chrom
    // -----------------------------------------------------------------------
    #[test]
    fn test_build_query_multi_chrom() {
        let q = build_query_multi_chrom(EnsemblEntityKind::Variation, "var", &["MT", "GL000220"]);
        assert!(q.contains("WHERE chrom IN ('MT', 'GL000220')"));
        assert!(q.contains("ORDER BY chrom, start"));
    }

    #[test]
    fn test_build_query_multi_chrom_transcript() {
        let q = build_query_multi_chrom(EnsemblEntityKind::Transcript, "tx", &["1", "2"]);
        assert!(q.contains("WHERE chrom IN ('1', '2')"));
        assert!(q.contains("ROW_NUMBER()"));
        assert!(q.contains("WHERE _rn = 1"));
    }

    // -----------------------------------------------------------------------
    // sorting_columns_for
    // -----------------------------------------------------------------------
    #[test]
    fn test_sorting_columns_for_all_present() {
        let schema = variation_schema();
        let result = sorting_columns_for(&schema, &["chrom", "start"]);
        assert!(result.is_some());
        let cols = result.unwrap();
        assert_eq!(cols.len(), 2);
    }

    #[test]
    fn test_sorting_columns_for_missing_column() {
        let schema = variation_schema();
        let result = sorting_columns_for(&schema, &["chrom", "nonexistent"]);
        assert!(result.is_none());
    }

    // -----------------------------------------------------------------------
    // project_batch
    // -----------------------------------------------------------------------
    #[test]
    fn test_project_batch_subset() {
        let batch = make_batch(vec!["1"], vec![100], vec![100], vec!["A/G"], vec!["rs1"]);
        let target = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let projected = project_batch(&batch, &target).unwrap();
        assert_eq!(projected.num_columns(), 2);
        assert_eq!(projected.num_rows(), 1);
        assert_eq!(projected.schema().field(0).name(), "chrom");
        assert_eq!(projected.schema().field(1).name(), "allele_string");
    }

    #[test]
    fn test_project_batch_missing_column() {
        let batch = make_batch(vec!["1"], vec![100], vec![100], vec!["A/G"], vec!["rs1"]);
        let target = Arc::new(Schema::new(vec![Field::new(
            "nonexistent",
            DataType::Utf8,
            false,
        )]));
        assert!(project_batch(&batch, &target).is_err());
    }

    // -----------------------------------------------------------------------
    // string_value
    // -----------------------------------------------------------------------
    #[test]
    fn test_string_value_utf8() {
        let arr = StringArray::from(vec!["hello", "world"]);
        assert_eq!(string_value(&arr, 0), "hello");
        assert_eq!(string_value(&arr, 1), "world");
    }

    #[test]
    fn test_string_value_large_utf8() {
        let arr = LargeStringArray::from(vec!["large", "string"]);
        assert_eq!(string_value(&arr, 0), "large");
        assert_eq!(string_value(&arr, 1), "string");
    }

    // -----------------------------------------------------------------------
    // PositionAccumulator — single batch
    // -----------------------------------------------------------------------
    #[test]
    fn test_accumulator_single_batch_single_row() {
        let batch = make_batch(vec!["1"], vec![100], vec![100], vec!["A/G"], vec!["rs1"]);
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").unwrap();
        let start_idx = schema.index_of("start").unwrap();
        let allele_idx = schema.index_of("allele_string").unwrap();
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_idx && i != start_idx)
            .collect();

        let accum = PositionAccumulator::new(1, "1".to_string(), 100, 0, &batch);
        let (key, value) = accum
            .finish_entry(&col_indices, allele_idx, &mut None)
            .unwrap();

        assert_eq!(key.len(), 10); // 2-byte chrom + 8-byte start
        assert!(!value.is_empty());
    }

    #[test]
    fn test_accumulator_single_batch_multiple_rows() {
        let batch = make_batch(
            vec!["1", "1"],
            vec![100, 100],
            vec![100, 100],
            vec!["A/G", "A/T"],
            vec!["rs1", "rs2"],
        );
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").unwrap();
        let start_idx = schema.index_of("start").unwrap();
        let allele_idx = schema.index_of("allele_string").unwrap();
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_idx && i != start_idx)
            .collect();

        let mut accum = PositionAccumulator::new(1, "1".to_string(), 100, 0, &batch);
        accum.add_row(1, &batch);

        // new() creates one segment, add_row() always pushes another
        assert_eq!(accum.segments.len(), 2);

        let (key, value) = accum
            .finish_entry(&col_indices, allele_idx, &mut None)
            .unwrap();
        assert!(!value.is_empty());

        // Verify the serialized entry has 2 alleles
        let reader = crate::kv_cache::position_entry::PositionEntryReader::new(&value).unwrap();
        assert_eq!(reader.num_alleles(), 2);
    }

    // -----------------------------------------------------------------------
    // PositionAccumulator — cross-batch
    // -----------------------------------------------------------------------
    #[test]
    fn test_accumulator_cross_batch_preserves_data() {
        // Create two different batches with the same position
        let batch1 = make_batch(vec!["1"], vec![100], vec![100], vec!["A/G"], vec!["rs1"]);
        let batch2 = make_batch(vec!["1"], vec![100], vec![100], vec!["A/T"], vec!["rs2"]);

        let schema = batch1.schema();
        let chrom_idx = schema.index_of("chrom").unwrap();
        let start_idx = schema.index_of("start").unwrap();
        let allele_idx = schema.index_of("allele_string").unwrap();
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_idx && i != start_idx)
            .collect();

        let mut accum = PositionAccumulator::new(1, "1".to_string(), 100, 0, &batch1);
        accum.add_row(0, &batch2); // Different batch!

        assert_eq!(
            accum.segments.len(),
            2,
            "Should have 2 segments for 2 batches"
        );

        let (key, value) = accum
            .finish_entry(&col_indices, allele_idx, &mut None)
            .unwrap();
        assert!(!value.is_empty());

        // Verify both alleles are preserved
        let reader = crate::kv_cache::position_entry::PositionEntryReader::new(&value).unwrap();
        assert_eq!(
            reader.num_alleles(),
            2,
            "Cross-batch accumulator must preserve both alleles"
        );
    }

    // -----------------------------------------------------------------------
    // PositionAccumulator — with compression
    // -----------------------------------------------------------------------
    #[test]
    fn test_accumulator_with_compression() {
        let batch = make_batch(vec!["1"], vec![100], vec![100], vec!["A/G"], vec!["rs1"]);
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").unwrap();
        let start_idx = schema.index_of("start").unwrap();
        let allele_idx = schema.index_of("allele_string").unwrap();
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_idx && i != start_idx)
            .collect();

        // Without compression
        let accum_raw = PositionAccumulator::new(1, "1".to_string(), 100, 0, &batch);
        let (_, raw_value) = accum_raw
            .finish_entry(&col_indices, allele_idx, &mut None)
            .unwrap();

        // With compression (no dictionary — plain zstd)
        let accum_comp = PositionAccumulator::new(1, "1".to_string(), 100, 0, &batch);
        let mut compressor = Some(zstd::bulk::Compressor::new(3).unwrap());
        let (_, comp_value) = accum_comp
            .finish_entry(&col_indices, allele_idx, &mut compressor)
            .unwrap();

        assert!(!raw_value.is_empty());
        assert!(!comp_value.is_empty());
        // Compressed may be larger for tiny payloads, but should be different bytes
        assert_ne!(raw_value, comp_value);
    }

    // -----------------------------------------------------------------------
    // CacheBuilder — builder pattern defaults
    // -----------------------------------------------------------------------
    #[test]
    fn test_cache_builder_defaults() {
        let builder = CacheBuilder::new("/cache", "/output");
        assert_eq!(builder.cache_root, "/cache");
        assert_eq!(builder.output_dir, "/output");
        assert_eq!(builder.partitions, 8);
        assert!(builder.build_fjall);
        assert_eq!(builder.zstd_level, 3);
        assert_eq!(builder.dict_size_kb, 112);
        assert!(builder.on_progress.is_none());
    }

    #[test]
    fn test_cache_builder_with_overrides() {
        let builder = CacheBuilder::new("/cache", "/output")
            .with_partitions(4)
            .with_build_fjall(false)
            .with_zstd_level(9)
            .with_dict_size_kb(256);
        assert_eq!(builder.partitions, 4);
        assert!(!builder.build_fjall);
        assert_eq!(builder.zstd_level, 9);
        assert_eq!(builder.dict_size_kb, 256);
    }

    #[test]
    fn test_cache_builder_with_progress() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let calls = Arc::new(AtomicUsize::new(0));
        let calls_clone = Arc::clone(&calls);
        let builder = CacheBuilder::new("/cache", "/output").with_on_progress(Box::new(
            move |_entity, _format, _batch, _total, _expected| {
                calls_clone.fetch_add(1, Ordering::SeqCst);
            },
        ));
        assert!(builder.on_progress.is_some());

        // Invoke callback to verify it works
        if let Some(ref cb) = builder.on_progress {
            cb("variation", "parquet", 100, 100, 0);
        }
        assert_eq!(calls.load(Ordering::SeqCst), 1);
    }

    // -----------------------------------------------------------------------
    // build_entity — unknown entity returns error
    // -----------------------------------------------------------------------
    #[tokio::test(flavor = "multi_thread")]
    async fn test_build_entity_unknown() {
        let builder = CacheBuilder::new("/nonexistent", "/nonexistent");
        let result = builder.build_entity("foobar").await;
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Unknown entity"));
    }

    // -----------------------------------------------------------------------
    // MAIN_CHROMS and CHROM_CODE_ORDER consistency
    // -----------------------------------------------------------------------
    #[test]
    fn test_main_chroms_in_code_order() {
        // Every MAIN_CHROMS entry must appear in CHROM_CODE_ORDER
        for chrom in MAIN_CHROMS {
            assert!(
                CHROM_CODE_ORDER.contains(chrom),
                "{chrom} in MAIN_CHROMS but not in CHROM_CODE_ORDER"
            );
        }
    }

    #[test]
    fn test_chrom_code_order_is_ascending() {
        // Verify CHROM_CODE_ORDER produces ascending chrom_codes
        let codes: Vec<u16> = CHROM_CODE_ORDER.iter().map(|c| chrom_to_code(c)).collect();
        for i in 1..codes.len() {
            assert!(
                codes[i] > codes[i - 1],
                "CHROM_CODE_ORDER not ascending at index {i}: {} (code {}) should be > {} (code {})",
                CHROM_CODE_ORDER[i],
                codes[i],
                CHROM_CODE_ORDER[i - 1],
                codes[i - 1]
            );
        }
    }

    // -----------------------------------------------------------------------
    // writer_properties
    // -----------------------------------------------------------------------
    #[test]
    fn test_writer_properties_variation() {
        let schema = variation_schema();
        let props = writer_properties(
            EnsemblEntityKind::Variation,
            &schema,
            &["chrom", "start"],
            None,
        );
        assert_eq!(props.max_row_group_size(), 100_000);
    }

    #[test]
    fn test_writer_properties_rg_override() {
        let schema = variation_schema();
        let props = writer_properties(
            EnsemblEntityKind::Variation,
            &schema,
            &["chrom", "start"],
            Some(256),
        );
        assert_eq!(props.max_row_group_size(), 256);
    }

    // -----------------------------------------------------------------------
    // SiftKvStore::ingest_sorted
    // -----------------------------------------------------------------------
    #[test]
    fn test_sift_ingest_sorted_roundtrip() {
        use crate::transcript_consequence::CompactPrediction;
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let entries = vec![
            (
                "ENST00000111111".to_string(),
                CachedPredictions {
                    sift: vec![CompactPrediction {
                        position: 10,
                        amino_acid: 1,
                        prediction: 0,
                        score: 0.05,
                    }],
                    polyphen: vec![],
                },
            ),
            (
                "ENST00000222222".to_string(),
                CachedPredictions {
                    sift: vec![],
                    polyphen: vec![CompactPrediction {
                        position: 20,
                        amino_acid: 2,
                        prediction: 2,
                        score: 0.88,
                    }],
                },
            ),
        ];

        let store = SiftKvStore::ingest_sorted(&db, entries.into_iter()).unwrap();

        let preds1 = store.get("ENST00000111111").unwrap().unwrap();
        assert_eq!(preds1.sift.len(), 1);
        assert_eq!(preds1.sift[0].position, 10);

        let preds2 = store.get("ENST00000222222").unwrap().unwrap();
        assert_eq!(preds2.polyphen.len(), 1);
        assert!((preds2.polyphen[0].score - 0.88).abs() < f32::EPSILON);

        assert!(store.get("MISSING").unwrap().is_none());
    }

    #[test]
    fn test_sift_ingest_sorted_empty() {
        let dir = tempfile::tempdir().unwrap();
        let db = fjall::Database::builder(dir.path())
            .cache_size(64 * 1024 * 1024)
            .open()
            .unwrap();

        let entries: Vec<(String, CachedPredictions)> = vec![];
        let _store = SiftKvStore::ingest_sorted(&db, entries.into_iter()).unwrap();
    }

    #[test]
    fn test_sift_ingest_sorted_persistence() {
        use crate::transcript_consequence::CompactPrediction;
        let dir = tempfile::tempdir().unwrap();

        // Write phase
        {
            let db = fjall::Database::builder(dir.path())
                .cache_size(64 * 1024 * 1024)
                .open()
                .unwrap();
            let entries = vec![(
                "ENST00000333333".to_string(),
                CachedPredictions {
                    sift: vec![CompactPrediction {
                        position: 42,
                        amino_acid: 3,
                        prediction: 1,
                        score: 0.5,
                    }],
                    polyphen: vec![],
                },
            )];
            SiftKvStore::ingest_sorted(&db, entries.into_iter()).unwrap();
            db.persist(fjall::PersistMode::SyncAll).unwrap();
        }

        // Read phase — reopen
        let store = SiftKvStore::open_path(dir.path())
            .unwrap()
            .expect("should find persisted sift data");
        let preds = store.get("ENST00000333333").unwrap().unwrap();
        assert_eq!(preds.sift.len(), 1);
        assert_eq!(preds.sift[0].position, 42);
    }

    /// Verify that the parallel execution path (Option B: temp parquet files)
    /// works against the real VEP 115 tabix cache.
    #[tokio::test(flavor = "multi_thread")]
    #[ignore] // requires real VEP cache at well-known path
    async fn test_tabix_variation_plan_is_parallel() {
        use datafusion::physical_plan::displayable;
        use std::time::Instant;

        let cache_root = "/Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38";
        if !std::path::Path::new(cache_root).exists() {
            eprintln!("Skipping: VEP cache not found at {cache_root}");
            return;
        }

        let partitions = 8;
        let ctx = make_ctx_and_register(
            cache_root,
            EnsemblEntityKind::Variation,
            "var",
            partitions,
        )
        .unwrap();

        let query = "SELECT * FROM var WHERE chrom = '1' ORDER BY chrom, start";
        let df = ctx.sql(query).await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let plan_str = format!("{}", displayable(plan.as_ref()).indent(true));
        eprintln!("Plan:\n{plan_str}");

        // Extract inner multi-partition plan (same as cache builder does)
        let inner = extract_multi_partition_inner(&plan);
        assert!(inner.is_some(), "should extract inner plan from SortPreservingMergeExec");
        let inner = inner.unwrap();
        let n = inner.properties().partitioning.partition_count();
        eprintln!("Inner plan: {n} partitions");
        assert!(n > 1, "expected >1 partitions, got {n}");

        // Phase 1: Execute all partitions in parallel via JoinSet (like cache builder)
        let task_ctx = ctx.task_ctx();
        let mut handles = tokio::task::JoinSet::new();
        let phase1_start = Instant::now();

        for i in 0..n {
            let stream = inner.execute(i, Arc::clone(&task_ctx)).unwrap();
            handles.spawn(async move {
                let mut stream = stream;
                let mut rows = 0usize;
                let mut first_start: Option<i64> = None;
                let mut last_start: Option<i64> = None;
                let mut prev_start: Option<i64> = None;
                while let Some(batch) = stream.next().await {
                    let batch = batch.unwrap();
                    if batch.num_rows() == 0 {
                        continue;
                    }
                    let starts = batch
                        .column_by_name("start")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<datafusion::arrow::array::Int64Array>()
                        .unwrap();
                    for j in 0..starts.len() {
                        let s = starts.value(j);
                        if first_start.is_none() {
                            first_start = Some(s);
                        }
                        if let Some(prev) = prev_start {
                            assert!(
                                s >= prev,
                                "partition {i}: rows not sorted: {s} < {prev}"
                            );
                        }
                        prev_start = Some(s);
                        last_start = Some(s);
                    }
                    rows += batch.num_rows();
                }
                (i, rows, first_start, last_start)
            });
        }

        let mut partition_results: Vec<(usize, usize, Option<i64>, Option<i64>)> = Vec::new();
        while let Some(result) = handles.join_next().await {
            let (idx, rows, first_start, last_start) = result.unwrap();
            partition_results.push((idx, rows, first_start, last_start));
        }
        partition_results.sort_by_key(|(idx, _, _, _)| *idx);
        let partition_rows: Vec<(usize, usize)> = partition_results
            .iter()
            .map(|(idx, rows, _, _)| (*idx, *rows))
            .collect();

        let phase1_elapsed = phase1_start.elapsed();
        let total_rows: usize = partition_rows.iter().map(|(_, r)| *r).sum();
        eprintln!("Phase 1 done: {total_rows} rows in {:.1}s", phase1_elapsed.as_secs_f64());
        for (idx, rows) in &partition_rows {
            eprintln!("  partition {idx}: {rows} rows");
        }

        assert!(total_rows > 100_000, "expected >100K rows for chr1, got {total_rows}");

        // Verify all partitions did useful work (not just 1)
        let active_partitions = partition_rows.iter().filter(|(_, r)| *r > 0).count();
        eprintln!("Active partitions (>0 rows): {active_partitions}/{n}");
        assert!(
            active_partitions > 1,
            "expected multiple active partitions, got {active_partitions}"
        );

        // Verify cross-partition ordering: last_start of partition i <= first_start of partition i+1
        eprintln!("Cross-partition ordering check:");
        for w in partition_results.windows(2) {
            let (idx_a, _, _, last_a) = &w[0];
            let (idx_b, _, first_b, _) = &w[1];
            if let (Some(last), Some(first)) = (last_a, first_b) {
                eprintln!(
                    "  partition {idx_a} last_start={last} <= partition {idx_b} first_start={first} : {}",
                    last <= first
                );
                assert!(
                    last <= first,
                    "cross-partition order violation: partition {idx_a} last={last} > partition {idx_b} first={first}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // other-chroms sorting and fjall key ordering
    // -----------------------------------------------------------------------

    #[test]
    fn test_other_chroms_sorted_by_chrom_code() {
        let mut contigs = vec![
            "GL000220.1".to_string(),
            "HG1012_PATCH".to_string(),
            "KI270733.1".to_string(),
        ];
        contigs.sort_by_key(|c| chrom_to_code(c));

        // Verify codes are ascending
        let codes: Vec<u16> = contigs.iter().map(|c| chrom_to_code(c)).collect();
        for w in codes.windows(2) {
            assert!(
                w[0] <= w[1],
                "chrom codes not ascending: {} (code {}) > {} (code {})",
                contigs[0],
                w[0],
                contigs[1],
                w[1]
            );
        }
    }

    #[test]
    fn test_other_chroms_after_main_chroms_in_chrom_code() {
        // All main chrom codes (1-25) must be less than non-canonical codes (>= 26)
        let main_max = CHROM_CODE_ORDER
            .iter()
            .map(|c| chrom_to_code(c))
            .max()
            .unwrap();
        let non_canonical_code = chrom_to_code("GL000220.1");
        assert!(
            non_canonical_code > main_max,
            "non-canonical code {non_canonical_code} should be > main max {main_max}"
        );
    }

    #[test]
    fn test_chrom_batches_ordering_main_then_other() {
        // Simulate the batch building logic
        let main_chroms = vec!["1".to_string(), "22".to_string(), "X".to_string()];
        let other_chroms = vec![
            "GL000220.1".to_string(),
            "HG1012_PATCH".to_string(),
        ];

        let mut batches: Vec<(String, bool)> = Vec::new();
        for chrom in CHROM_CODE_ORDER {
            if main_chroms.iter().any(|c| c == chrom) {
                batches.push((chrom.to_string(), false));
            }
        }
        let mut sorted_other = other_chroms;
        sorted_other.sort_by_key(|c| chrom_to_code(c));
        for chrom in &sorted_other {
            batches.push((chrom.clone(), true));
        }

        // Verify all main chroms come first
        let main_end = batches.iter().position(|(_, is_other)| *is_other);
        if let Some(idx) = main_end {
            assert!(
                batches[..idx].iter().all(|(_, is_other)| !is_other),
                "all entries before first other should be main chroms"
            );
            assert!(
                batches[idx..].iter().all(|(_, is_other)| *is_other),
                "all entries after first other should be other chroms"
            );
        }

        // Verify chrom_codes are globally ascending
        let codes: Vec<u16> = batches.iter().map(|(c, _)| chrom_to_code(c)).collect();
        for w in codes.windows(2) {
            assert!(
                w[0] < w[1],
                "chrom codes not strictly ascending: {} >= {}",
                w[0],
                w[1]
            );
        }
    }

    #[test]
    fn test_build_query_single_chrom_for_other() {
        // Verify that per-contig queries use single WHERE chrom = '...'
        // instead of massive IN clause
        let q = build_query(EnsemblEntityKind::Variation, "var", Some("GL000220.1"));
        assert_eq!(
            q,
            "SELECT * FROM var WHERE chrom = 'GL000220.1' ORDER BY chrom, start"
        );
        // Should NOT contain IN clause
        assert!(!q.contains("IN ("));
    }

    #[tokio::test(flavor = "multi_thread")]
    #[ignore] // requires real VEP cache
    async fn test_real_cache_chroms_from_schema() {
        let cache_root = "/Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38";
        if !std::path::Path::new(cache_root).exists() {
            return;
        }
        let ctx = make_ctx_and_register(
            cache_root,
            EnsemblEntityKind::Variation,
            "var",
            1,
        )
        .unwrap();
        let table = ctx.table("var").await.unwrap();
        let schema = table.schema().inner().clone();
        let chroms = chroms_from_schema(&schema);
        eprintln!("chroms_from_schema: {:?}", chroms.as_ref().map(|c| c.len()));
        if let Some(ref all) = chroms {
            let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
            let main_count = all.iter().filter(|c| main_set.contains(c.as_str())).count();
            let other_count = all.len() - main_count;
            eprintln!("main: {main_count}, other: {other_count}");
            eprintln!("first 10: {:?}", &all[..all.len().min(10)]);
            eprintln!("contains MT: {}", all.iter().any(|c| c == "MT"));
        } else {
            eprintln!("WARNING: chroms_from_schema returned None!");
            eprintln!("Schema metadata keys: {:?}", schema.metadata().keys().collect::<Vec<_>>());
        }
        assert!(chroms.is_some(), "chroms_from_schema should not be None for real VEP cache");
        let all = chroms.unwrap();
        assert!(all.len() > 25, "expected >25 chroms, got {}", all.len());
    }

    #[test]
    fn test_mt_included_in_main_chrom_code_order() {
        // MT is in CHROM_CODE_ORDER but not in MAIN_CHROMS — verify it's
        // processed via the main chrom loop (it has a canonical code)
        assert!(CHROM_CODE_ORDER.contains(&"MT"));
        assert!(!MAIN_CHROMS.contains(&"MT"));
        // MT gets processed in the main loop IF present in the discovered chroms.
        // Its code (25) is less than NON_CANONICAL_START (26), so fjall ordering
        // is correct: main codes 1-25 < non-canonical codes >= 26.
        assert!(chrom_to_code("MT") < 26);
    }

    // -----------------------------------------------------------------------
    // extract_multi_partition_inner
    // -----------------------------------------------------------------------

    #[tokio::test(flavor = "multi_thread")]
    async fn test_extract_multi_partition_inner_non_merge_returns_none() {
        // A simple plan (no SortPreservingMergeExec) should return None
        let ctx = SessionContext::new();
        let df = ctx
            .sql("SELECT 1 as x")
            .await
            .unwrap();
        let plan = df.create_physical_plan().await.unwrap();
        assert!(
            extract_multi_partition_inner(&plan).is_none(),
            "non-merge plan should return None"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    #[ignore] // requires real VEP cache
    async fn test_parallel_execution_produces_same_results() {
        let cache_root = "/Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38";
        if !std::path::Path::new(cache_root).exists() {
            return;
        }

        // Run with 1 partition (sequential) and 8 partitions (parallel)
        // on a small chromosome to compare row counts.
        for partitions in [1, 8] {
            let ctx = make_ctx_and_register(
                cache_root,
                EnsemblEntityKind::Variation,
                "var",
                partitions,
            )
            .unwrap();

            let query = "SELECT COUNT(*) as cnt FROM var WHERE chrom = '22'";
            let batches = ctx.sql(query).await.unwrap().collect().await.unwrap();
            let count = batches[0]
                .column(0)
                .as_any()
                .downcast_ref::<datafusion::arrow::array::Int64Array>()
                .unwrap()
                .value(0);
            eprintln!("partitions={partitions}: chr22 count={count}");

            // chr22 should have a substantial number of rows
            assert!(count > 100_000, "chr22 too few rows: {count}");
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    #[ignore] // requires real VEP cache
    async fn test_parallel_plan_extracts_inner() {
        let cache_root = "/Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38";
        if !std::path::Path::new(cache_root).exists() {
            return;
        }

        let ctx = make_ctx_and_register(
            cache_root,
            EnsemblEntityKind::Variation,
            "var",
            8,
        )
        .unwrap();

        let query = "SELECT * FROM var WHERE chrom = '22' ORDER BY chrom, start";
        let df = ctx.sql(query).await.unwrap();
        let plan = df.create_physical_plan().await.unwrap();

        let inner = extract_multi_partition_inner(&plan);
        assert!(inner.is_some(), "should extract inner plan from SortPreservingMergeExec");

        let inner = inner.unwrap();
        let n = inner.properties().partitioning.partition_count();
        eprintln!("Inner plan: {} partitions", n);
        assert!(n > 1, "inner plan should have >1 partitions, got {n}");
    }

    // -----------------------------------------------------------------------
    // dir_has_parquet_files
    // -----------------------------------------------------------------------

    #[test]
    fn test_dir_has_parquet_files_with_files() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("chr1.parquet"), b"PAR1").unwrap();
        assert!(dir_has_parquet_files(dir.path().to_str().unwrap()));
    }

    #[test]
    fn test_dir_has_parquet_files_empty_dir() {
        let dir = tempfile::tempdir().unwrap();
        assert!(!dir_has_parquet_files(dir.path().to_str().unwrap()));
    }

    #[test]
    fn test_dir_has_parquet_files_nonexistent() {
        assert!(!dir_has_parquet_files("/nonexistent/path"));
    }

    #[test]
    fn test_dir_has_parquet_files_only_non_parquet() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("data.csv"), b"a,b").unwrap();
        assert!(!dir_has_parquet_files(dir.path().to_str().unwrap()));
    }

    // -----------------------------------------------------------------------
    // overwrite flag
    // -----------------------------------------------------------------------

    #[test]
    fn test_builder_overwrite_default_false() {
        let builder = CacheBuilder::new("/cache", "/output");
        assert!(!builder.overwrite);
    }

    #[test]
    fn test_builder_with_overwrite() {
        let builder = CacheBuilder::new("/cache", "/output").with_overwrite(true);
        assert!(builder.overwrite);
    }

    // -----------------------------------------------------------------------
    // skip logic integration
    // -----------------------------------------------------------------------

    #[tokio::test(flavor = "multi_thread")]
    async fn test_build_entity_skips_existing_parquet() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().to_str().unwrap();

        // Create a fake parquet file for transcript
        let transcript_dir = dir.path().join("transcript");
        std::fs::create_dir_all(&transcript_dir).unwrap();
        std::fs::write(transcript_dir.join("chr1.parquet"), b"PAR1").unwrap();

        let builder = CacheBuilder::new("/nonexistent_cache", output);
        // Should skip because parquet exists and overwrite=false
        let result = builder.build_entity("transcript").await;
        // Will succeed (skip) even though cache doesn't exist
        assert!(result.is_ok());
        let stats = result.unwrap();
        assert_eq!(stats[0].entity, "transcript");
        assert!(stats[0].parquet_files.is_empty(), "should return empty (skipped)");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_build_entity_does_not_skip_with_overwrite() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().to_str().unwrap();

        // Create a fake parquet file for transcript
        let transcript_dir = dir.path().join("transcript");
        std::fs::create_dir_all(&transcript_dir).unwrap();
        std::fs::write(transcript_dir.join("chr1.parquet"), b"PAR1").unwrap();

        let builder = CacheBuilder::new("/nonexistent_cache", output).with_overwrite(true);
        // Should NOT skip because overwrite=true — will fail because cache doesn't exist
        let result = builder.build_entity("transcript").await;
        assert!(result.is_err(), "should fail (not skip) with overwrite=true");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_variation_skips_when_both_exist() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().to_str().unwrap();

        // Create fake variation parquet + fjall
        let var_dir = dir.path().join("variation");
        std::fs::create_dir_all(&var_dir).unwrap();
        std::fs::write(var_dir.join("chr1.parquet"), b"PAR1").unwrap();
        let fjall_dir = dir.path().join("variation.fjall");
        std::fs::create_dir_all(&fjall_dir).unwrap();

        let builder = CacheBuilder::new("/nonexistent_cache", output);
        let result = builder.build_entity("variation").await;
        assert!(result.is_ok());
        let stats = result.unwrap();
        assert!(stats[0].parquet_files.is_empty(), "should skip");
        assert!(stats[0].fjall_stats.is_none(), "should skip fjall too");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_variation_rebuilds_fjall_when_parquet_exists() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().to_str().unwrap();

        // Create fake variation parquet but NO fjall
        let var_dir = dir.path().join("variation");
        std::fs::create_dir_all(&var_dir).unwrap();
        std::fs::write(var_dir.join("chr1.parquet"), b"PAR1").unwrap();

        let builder = CacheBuilder::new("/nonexistent_cache", output);
        let result = builder.build_entity("variation").await;
        // Will try to rebuild fjall from parquet but fail because the
        // parquet file is not a real parquet file
        assert!(result.is_err(), "should attempt fjall rebuild and fail on fake parquet");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_translation_skips_when_all_exist() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().to_str().unwrap();

        // Create all translation outputs
        let core_dir = dir.path().join("translation_core");
        std::fs::create_dir_all(&core_dir).unwrap();
        std::fs::write(core_dir.join("chr1.parquet"), b"PAR1").unwrap();
        let sift_dir = dir.path().join("translation_sift");
        std::fs::create_dir_all(&sift_dir).unwrap();
        std::fs::write(sift_dir.join("chr1.parquet"), b"PAR1").unwrap();
        let fjall_dir = dir.path().join("translation_sift.fjall");
        std::fs::create_dir_all(&fjall_dir).unwrap();

        let builder = CacheBuilder::new("/nonexistent_cache", output);
        let result = builder.build_entity("translation").await;
        assert!(result.is_ok());
        let stats = result.unwrap();
        assert!(stats[0].parquet_files.is_empty(), "should skip");
    }
}
