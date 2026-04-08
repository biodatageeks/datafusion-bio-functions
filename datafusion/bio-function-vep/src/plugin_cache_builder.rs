//! Plugin cache builder: converts external plugin source files to parquet and fjall.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use datafusion::arrow::array::{
    Array, ArrayRef, Float32Array, Float32Builder, Int32Array, Int32Builder, StringArray,
    StringBuilder, UInt32Array, UInt32Builder,
};
use datafusion::arrow::compute::cast;
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::file_format::file_compression_type::FileCompressionType;
use datafusion::parquet::arrow::ArrowWriter;
use datafusion::parquet::basic::Compression as ParquetCompression;
use datafusion::parquet::file::properties::WriterProperties;
use datafusion::parquet::format::SortingColumn;
use datafusion::parquet::schema::types::ColumnPath;
use datafusion::prelude::{CsvReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use flate2::Compression as GzCompression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use futures::StreamExt;
use once_cell::sync::Lazy;
use tempfile::{Builder, NamedTempFile};

use crate::kv_cache::CacheLoader;
use crate::plugin::PluginKind;

static SANITIZED_TSVS: Lazy<Mutex<Vec<NamedTempFile>>> = Lazy::new(|| Mutex::new(Vec::new()));

const MAIN_CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y",
];

const PLUGIN_ROW_GROUP_SIZE: usize = 100_000;

pub fn convert_plugin(
    source_path: &str,
    output_dir: &str,
    plugin: PluginKind,
    partitions: usize,
    memory_limit_gb: usize,
) -> Result<Vec<(String, usize)>> {
    match plugin {
        PluginKind::ClinVar => convert_plugin_to_parquet(
            "clinvar",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
        ),
        PluginKind::Cadd => {
            let (snv_source, indel_source) = resolve_cadd_sources(source_path)?;
            convert_cadd_sources_to_parquet(
                &snv_source,
                &indel_source,
                output_dir,
                partitions,
                memory_limit_gb,
                None,
            )
        }
        PluginKind::SpliceAI => convert_plugin_to_parquet(
            "spliceai",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
        ),
        PluginKind::AlphaMissense => convert_plugin_to_parquet(
            "alphamissense",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
        ),
        PluginKind::DbNSFP => convert_plugin_to_parquet(
            "dbnsfp",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
        ),
    }
}

pub fn convert_plugin_to_parquet(
    plugin_name: &str,
    source_path: &str,
    output_dir: &str,
    partitions: usize,
    memory_limit_gb: usize,
    chromosomes: Option<Vec<String>>,
) -> Result<Vec<(String, usize)>> {
    let rt = tokio::runtime::Runtime::new()
        .map_err(|e| DataFusionError::Execution(format!("Failed to create runtime: {e}")))?;
    rt.block_on(convert_plugin_to_parquet_async(
        plugin_name,
        source_path,
        output_dir,
        partitions,
        memory_limit_gb,
        chromosomes,
    ))
}

pub fn convert_cadd_sources_to_parquet(
    snv_source_path: &str,
    indel_source_path: &str,
    output_dir: &str,
    partitions: usize,
    memory_limit_gb: usize,
    chromosomes: Option<Vec<String>>,
) -> Result<Vec<(String, usize)>> {
    let rt = tokio::runtime::Runtime::new()
        .map_err(|e| DataFusionError::Execution(format!("Failed to create runtime: {e}")))?;
    rt.block_on(convert_cadd_sources_to_parquet_async(
        snv_source_path,
        indel_source_path,
        output_dir,
        partitions,
        memory_limit_gb,
        chromosomes,
    ))
}

pub fn build_plugin_fjall_from_parquet(
    plugin_name: &str,
    parquet_dir: &str,
    output_path: &str,
    partitions: usize,
    chromosomes: Option<Vec<String>>,
) -> Result<(String, usize)> {
    let rt = tokio::runtime::Runtime::new()
        .map_err(|e| DataFusionError::Execution(format!("Failed to create runtime: {e}")))?;
    rt.block_on(build_plugin_fjall_from_parquet_async(
        plugin_name,
        parquet_dir,
        output_path,
        partitions,
        chromosomes,
    ))
}

async fn convert_plugin_to_parquet_async(
    plugin_name: &str,
    source_path: &str,
    output_dir: &str,
    partitions: usize,
    _memory_limit_gb: usize,
    chromosomes: Option<Vec<String>>,
) -> Result<Vec<(String, usize)>> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);

    match plugin_name {
        "clinvar" | "spliceai" => register_vcf_source(&ctx, source_path, "source").await?,
        "cadd" | "alphamissense" | "dbnsfp" => {
            register_tsv_source(&ctx, source_path, "source").await?
        }
        _ => {
            return Err(DataFusionError::Execution(format!(
                "Unknown plugin: {plugin_name}"
            )));
        }
    };

    let select_query = plugin_select_query(plugin_name)?;
    let chrom_col = plugin_chrom_column(plugin_name);
    let chroms = discover_chromosomes(&ctx, "source", chrom_col).await?;
    let requested: Option<HashSet<String>> = chromosomes.map(|values| {
        values
            .into_iter()
            .map(|value| normalize_chrom(&value))
            .collect()
    });
    let grouped = group_chromosomes(chroms, requested.as_ref());
    let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
    let (main_chroms, other_chroms): (Vec<(String, Vec<String>)>, Vec<(String, Vec<String>)>) =
        grouped
            .into_iter()
            .partition(|(canonical, _)| main_set.contains(canonical.as_str()));

    eprintln!(
        "  {plugin_name}: {} main chroms, {} other contigs",
        main_chroms.len(),
        other_chroms.len()
    );

    let mut all_results = Vec::new();
    let global_start = Instant::now();
    let mut total_rows = 0usize;

    for (norm, raw_chroms) in &main_chroms {
        let where_clause = build_where_clause(chrom_col, raw_chroms);
        let query = format!("{select_query}{where_clause} ORDER BY chrom, pos, ref, alt");
        let output_file = format!("{output_dir}/chr{norm}.parquet");

        let rows = write_query_to_parquet(&ctx, plugin_name, &query, &output_file).await?;
        total_rows += rows;
        if rows > 0 {
            eprintln!("  {plugin_name}: chr{norm} {rows} rows (total: {total_rows})");
        }
        all_results.push((output_file, rows));
    }

    if !other_chroms.is_empty() {
        let all_other_raw: Vec<String> = other_chroms
            .iter()
            .flat_map(|(_, raws)| raws.iter().cloned())
            .collect();
        let where_clause = build_where_clause(chrom_col, &all_other_raw);
        let query = format!("{select_query}{where_clause} ORDER BY chrom, pos, ref, alt");
        let output_file = format!("{output_dir}/other.parquet");

        let rows = write_query_to_parquet(&ctx, plugin_name, &query, &output_file).await?;
        total_rows += rows;
        eprintln!(
            "  {plugin_name}: other ({} contigs) {rows} rows",
            other_chroms.len()
        );
        all_results.push((output_file, rows));
    }

    eprintln!(
        "  {plugin_name}: {total_rows} total rows [{:.1}s]",
        global_start.elapsed().as_secs_f64()
    );
    Ok(all_results)
}

async fn convert_cadd_sources_to_parquet_async(
    snv_source_path: &str,
    indel_source_path: &str,
    output_dir: &str,
    partitions: usize,
    _memory_limit_gb: usize,
    chromosomes: Option<Vec<String>>,
) -> Result<Vec<(String, usize)>> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);
    register_tsv_source(&ctx, snv_source_path, "source_snv").await?;
    register_tsv_source(&ctx, indel_source_path, "source_indel").await?;

    let mut chroms = discover_chromosomes(&ctx, "source_snv", "chrom").await?;
    chroms.extend(discover_chromosomes(&ctx, "source_indel", "chrom").await?);
    let requested: Option<HashSet<String>> = chromosomes.map(|values| {
        values
            .into_iter()
            .map(|value| normalize_chrom(&value))
            .collect()
    });
    let grouped = group_chromosomes(chroms, requested.as_ref());
    let main_set: HashSet<&str> = MAIN_CHROMS.iter().copied().collect();
    let (main_chroms, other_chroms): (Vec<(String, Vec<String>)>, Vec<(String, Vec<String>)>) =
        grouped
            .into_iter()
            .partition(|(canonical, _)| main_set.contains(canonical.as_str()));

    eprintln!(
        "  cadd: {} main chroms, {} other contigs",
        main_chroms.len(),
        other_chroms.len()
    );

    let mut all_results = Vec::new();
    let global_start = Instant::now();
    let mut total_rows = 0usize;

    for (norm, raw_chroms) in &main_chroms {
        let query = cadd_union_query(raw_chroms);
        let output_file = format!("{output_dir}/chr{norm}.parquet");
        let rows = write_query_to_parquet(&ctx, "cadd", &query, &output_file).await?;
        total_rows += rows;
        if rows > 0 {
            eprintln!("  cadd: chr{norm} {rows} rows (total: {total_rows})");
        }
        all_results.push((output_file, rows));
    }

    if !other_chroms.is_empty() {
        let all_other_raw: Vec<String> = other_chroms
            .iter()
            .flat_map(|(_, raws)| raws.iter().cloned())
            .collect();
        let query = cadd_union_query(&all_other_raw);
        let output_file = format!("{output_dir}/other.parquet");
        let rows = write_query_to_parquet(&ctx, "cadd", &query, &output_file).await?;
        total_rows += rows;
        eprintln!("  cadd: other ({} contigs) {rows} rows", other_chroms.len());
        all_results.push((output_file, rows));
    }

    eprintln!(
        "  cadd: {total_rows} total rows [{:.1}s]",
        global_start.elapsed().as_secs_f64()
    );
    Ok(all_results)
}

async fn build_plugin_fjall_from_parquet_async(
    plugin_name: &str,
    parquet_dir: &str,
    output_path: &str,
    partitions: usize,
    chromosomes: Option<Vec<String>>,
) -> Result<(String, usize)> {
    if !Path::new(parquet_dir).is_dir() {
        return Err(DataFusionError::Execution(format!(
            "Plugin parquet directory not found: {parquet_dir}"
        )));
    }

    if Path::new(output_path).exists() {
        std::fs::remove_dir_all(output_path).map_err(|e| {
            DataFusionError::Execution(format!(
                "Failed to remove existing plugin fjall store {output_path}: {e}"
            ))
        })?;
    }

    let ctx =
        SessionContext::new_with_config(SessionConfig::new().with_target_partitions(partitions));
    ctx.register_parquet("plugin_parquet", parquet_dir, Default::default())
        .await?;

    let select_expr = plugin_fjall_projection(plugin_name)?;
    let view_name = if chromosomes
        .as_ref()
        .is_some_and(|values| !values.is_empty())
    {
        let chrom_filter = chromosomes
            .unwrap_or_default()
            .into_iter()
            .flat_map(|chrom| chrom_aliases_for_sql(&chrom))
            .collect::<Vec<_>>()
            .join(", ");
        let sql = format!(
            "CREATE VIEW plugin_source AS \
             SELECT {select_expr} FROM plugin_parquet \
             WHERE chrom IN ({chrom_filter})"
        );
        ctx.sql(&sql).await?;
        "plugin_source"
    } else {
        let sql = format!("CREATE VIEW plugin_source AS SELECT {select_expr} FROM plugin_parquet");
        ctx.sql(&sql).await?;
        "plugin_source"
    };

    let loader = CacheLoader::new(view_name, output_path).with_parallelism(partitions);
    let stats = loader.load(&ctx).await?;
    Ok((output_path.to_string(), stats.total_variants as usize))
}

fn group_chromosomes(
    chroms: Vec<String>,
    requested: Option<&HashSet<String>>,
) -> Vec<(String, Vec<String>)> {
    let mut grouped: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for chrom in chroms {
        let canonical = normalize_chrom(&chrom);
        if requested.is_some_and(|wanted| !wanted.contains(&canonical)) {
            continue;
        }
        grouped.entry(canonical).or_default().insert(chrom);
    }
    grouped
        .into_iter()
        .map(|(canonical, raw_values)| (canonical, raw_values.into_iter().collect()))
        .collect()
}

fn build_where_clause(chrom_col: &str, raw_chroms: &[String]) -> String {
    if raw_chroms.len() == 1 {
        format!(
            " WHERE {chrom_col} = '{}'",
            raw_chroms[0].replace('\'', "''")
        )
    } else {
        let in_list = raw_chroms
            .iter()
            .map(|value| format!("'{}'", value.replace('\'', "''")))
            .collect::<Vec<_>>()
            .join(", ");
        format!(" WHERE {chrom_col} IN ({in_list})")
    }
}

fn cadd_union_query(raw_chroms: &[String]) -> String {
    let snv_where = build_where_clause("chrom", raw_chroms);
    let indel_where = build_where_clause("chrom", raw_chroms);
    format!(
        "SELECT * FROM (\
         SELECT chrom AS chrom, CAST(pos AS INTEGER) AS pos, ref AS ref, alt AS alt, CAST(rawscore AS FLOAT) AS raw_score, CAST(phred AS FLOAT) AS phred_score FROM source_snv{snv_where} \
         UNION ALL \
         SELECT chrom AS chrom, CAST(pos AS INTEGER) AS pos, ref AS ref, alt AS alt, CAST(rawscore AS FLOAT) AS raw_score, CAST(phred AS FLOAT) AS phred_score FROM source_indel{indel_where}\
         ) ORDER BY chrom, pos, ref, alt"
    )
}

fn resolve_cadd_sources(source_path: &str) -> Result<(String, String)> {
    let path = Path::new(source_path);
    if path.is_file()
        && path.file_name().and_then(|name| name.to_str()) == Some("whole_genome_SNVs.tsv.gz")
    {
        let indel = path.with_file_name("gnomad.genomes.r4.0.indel.tsv.gz");
        if indel.is_file() {
            return Ok((
                path.to_string_lossy().into_owned(),
                indel.to_string_lossy().into_owned(),
            ));
        }
    }
    if path.is_dir() {
        let snv = path.join("whole_genome_SNVs.tsv.gz");
        let indel = path.join("gnomad.genomes.r4.0.indel.tsv.gz");
        if snv.is_file() && indel.is_file() {
            return Ok((
                snv.to_string_lossy().into_owned(),
                indel.to_string_lossy().into_owned(),
            ));
        }
        return Err(DataFusionError::Execution(format!(
            "CADD source directory must contain whole_genome_SNVs.tsv.gz and gnomad.genomes.r4.0.indel.tsv.gz: {source_path}"
        )));
    }

    let manifest_text = std::fs::read_to_string(path).map_err(|e| {
        DataFusionError::Execution(format!(
            "Failed to read CADD source manifest {source_path}: {e}"
        ))
    })?;
    let manifest: serde_json::Value = serde_json::from_str(&manifest_text).map_err(|e| {
        DataFusionError::Execution(format!(
            "CADD source manifest must be valid JSON with 'snv' and 'indel' keys: {e}"
        ))
    })?;
    let snv = manifest
        .get("snv")
        .and_then(|value| value.as_str())
        .ok_or_else(|| {
            DataFusionError::Execution("CADD source manifest is missing string key 'snv'".into())
        })?;
    let indel = manifest
        .get("indel")
        .and_then(|value| value.as_str())
        .ok_or_else(|| {
            DataFusionError::Execution("CADD source manifest is missing string key 'indel'".into())
        })?;
    Ok((snv.to_string(), indel.to_string()))
}

async fn register_vcf_source(
    ctx: &SessionContext,
    source_path: &str,
    table_name: &str,
) -> Result<()> {
    let provider = VcfTableProvider::new(source_path.to_string(), None, None, None, false)?;
    ctx.register_table(table_name, Arc::new(provider))?;
    Ok(())
}

async fn register_tsv_source(
    ctx: &SessionContext,
    source_path: &str,
    table_name: &str,
) -> Result<()> {
    let (table_path, is_gzip) = sanitize_plugin_source(source_path)?;
    let mut options = CsvReadOptions::new()
        .delimiter(b'\t')
        .has_header(true)
        .file_extension(".gz");
    if is_gzip {
        options = options.file_compression_type(FileCompressionType::GZIP);
    }
    ctx.register_csv(table_name, &table_path, options).await?;
    Ok(())
}

fn sanitize_plugin_source(source_path: &str) -> Result<(String, bool)> {
    let mut source_file = File::open(source_path).map_err(|e| {
        DataFusionError::Execution(format!("Failed to open plugin source {source_path}: {e}"))
    })?;

    let mut signature = [0u8; 2];
    let is_gzip = match source_file.read_exact(&mut signature) {
        Ok(_) => {
            source_file.seek(SeekFrom::Start(0)).map_err(|e| {
                DataFusionError::Execution(format!("Failed to seek plugin source: {e}"))
            })?;
            signature == [0x1f, 0x8b]
        }
        Err(_) => {
            source_file.seek(SeekFrom::Start(0)).map_err(|e| {
                DataFusionError::Execution(format!("Failed to seek plugin source: {e}"))
            })?;
            false
        }
    };

    if !is_gzip {
        return Ok((source_path.to_string(), false));
    }

    let temp = Builder::new()
        .prefix("vepyr_plugin_tsv_")
        .suffix(".sanitized.tsv.gz")
        .tempfile()
        .map_err(|e| {
            DataFusionError::Execution(format!("Failed to create sanitized plugin source: {e}"))
        })?;
    let mut reader = BufReader::new(MultiGzDecoder::new(source_file));
    let writer = temp.reopen().map_err(|e| {
        DataFusionError::Execution(format!("Failed to open temp file for writing: {e}"))
    })?;
    let mut encoder = GzEncoder::new(writer, GzCompression::default());

    let mut buffer = String::new();
    let mut header_written = false;

    while reader.read_line(&mut buffer)? != 0 {
        let line = buffer.trim_end_matches(&['\r', '\n'][..]);
        if line.is_empty() {
            buffer.clear();
            continue;
        }
        if !header_written {
            if line.starts_with("##") || !line.contains('\t') {
                buffer.clear();
                continue;
            }
            let header_line = line.trim_start_matches('#').to_lowercase();
            encoder.write_all(header_line.as_bytes())?;
            encoder.write_all(b"\n")?;
            header_written = true;
        } else {
            encoder.write_all(line.as_bytes())?;
            encoder.write_all(b"\n")?;
        }
        buffer.clear();
    }

    if !header_written {
        return Err(DataFusionError::Execution(format!(
            "Plugin source {source_path} missing header row"
        )));
    }

    encoder.finish().map_err(|e| {
        DataFusionError::Execution(format!("Failed to finish sanitized source: {e}"))
    })?;

    let sanitized_path = temp.path().to_string_lossy().into_owned();
    let mut storage = SANITIZED_TSVS
        .lock()
        .map_err(|e| DataFusionError::Execution(format!("Sanitized path cache poisoned: {e}")))?;
    storage.push(temp);
    Ok((sanitized_path, true))
}

fn plugin_chrom_column(plugin_name: &str) -> &'static str {
    match plugin_name {
        "clinvar" | "spliceai" => "chrom",
        "cadd" => "chrom",
        "alphamissense" => "chrom",
        "dbnsfp" => "chr",
        _ => "chrom",
    }
}

fn plugin_select_query(plugin_name: &str) -> Result<String> {
    let query = match plugin_name {
        "clinvar" => {
            "SELECT chrom, CAST(start AS INTEGER) AS pos, ref AS ref, alt AS alt, array_to_string(\"CLNSIG\", '|') AS clnsig, array_to_string(\"CLNREVSTAT\", '|') AS clnrevstat, array_to_string(\"CLNDN\", '|') AS clndn, \"CLNVC\" AS clnvc, array_to_string(\"CLNVI\", '|') AS clnvi, CAST(\"AF_ESP\" AS FLOAT) AS af_esp, CAST(\"AF_EXAC\" AS FLOAT) AS af_exac, CAST(\"AF_TGP\" AS FLOAT) AS af_tgp FROM source"
        }
        "cadd" => {
            "SELECT chrom AS chrom, CAST(pos AS INTEGER) AS pos, ref AS ref, alt AS alt, CAST(rawscore AS FLOAT) AS raw_score, CAST(phred AS FLOAT) AS phred_score FROM source"
        }
        "alphamissense" => {
            "SELECT chrom AS chrom, CAST(pos AS INTEGER) AS pos, ref AS ref, alt AS alt, genome, uniprot_id, transcript_id, protein_variant, CAST(am_pathogenicity AS FLOAT) AS am_pathogenicity, am_class FROM source"
        }
        "spliceai" => {
            "SELECT chrom, CAST(start AS INTEGER) AS pos, ref AS ref, alt AS alt, split_part(\"SpliceAI\", '|', 2) AS symbol, CAST(split_part(\"SpliceAI\", '|', 3) AS FLOAT) AS ds_ag, CAST(split_part(\"SpliceAI\", '|', 4) AS FLOAT) AS ds_al, CAST(split_part(\"SpliceAI\", '|', 5) AS FLOAT) AS ds_dg, CAST(split_part(\"SpliceAI\", '|', 6) AS FLOAT) AS ds_dl, CAST(NULLIF(regexp_replace(split_part(\"SpliceAI\", '|', 7), '[^0-9-]', ''), '') AS INTEGER) AS dp_ag, CAST(NULLIF(regexp_replace(split_part(\"SpliceAI\", '|', 8), '[^0-9-]', ''), '') AS INTEGER) AS dp_al, CAST(NULLIF(regexp_replace(split_part(\"SpliceAI\", '|', 9), '[^0-9-]', ''), '') AS INTEGER) AS dp_dg, CAST(NULLIF(regexp_replace(split_part(\"SpliceAI\", '|', 10), '[^0-9-]', ''), '') AS INTEGER) AS dp_dl FROM source"
        }
        "dbnsfp" => {
            "SELECT chr AS chrom, CAST(\"pos(1-based)\" AS INTEGER) AS pos, ref AS ref, alt AS alt, sift4g_score, sift4g_pred, polyphen2_hdiv_score, polyphen2_hvar_score, lrt_score, lrt_pred, mutationtaster_score, mutationtaster_pred, fathmm_score, fathmm_pred, provean_score, provean_pred, vest4_score, metasvm_score, metasvm_pred, metalr_score, metalr_pred, CAST(revel_score AS FLOAT) AS revel_score, CAST(\"gerp++_rs\" AS FLOAT) AS gerp_rs, CAST(phyloP100way_vertebrate AS FLOAT) AS phylop100way, CAST(phyloP30way_mammalian AS FLOAT) AS phylop30way, CAST(phastCons100way_vertebrate AS FLOAT) AS phastcons100way, CAST(phastCons30way_mammalian AS FLOAT) AS phastcons30way, CAST(\"siphy_29way_logodds\" AS FLOAT) AS siphy_29way, CAST(cadd_raw AS FLOAT) AS cadd_raw, CAST(cadd_phred AS FLOAT) AS cadd_phred FROM source"
        }
        _ => {
            return Err(DataFusionError::Execution(format!(
                "Unknown plugin: {plugin_name}"
            )));
        }
    };
    Ok(query.to_string())
}

fn plugin_fjall_projection(plugin_name: &str) -> Result<&'static str> {
    match plugin_name {
        "clinvar" => Ok(
            "chrom, CAST(pos AS BIGINT) AS start, CAST(pos AS BIGINT) AS end, concat(ref, '/', alt) AS allele_string, clnsig, clnrevstat, clndn, clnvc, clnvi, af_esp, af_exac, af_tgp",
        ),
        "cadd" => Ok(
            "chrom, CAST(pos AS BIGINT) AS start, CAST(pos AS BIGINT) AS end, concat(ref, '/', alt) AS allele_string, raw_score, phred_score",
        ),
        "spliceai" => Ok(
            "chrom, CAST(pos AS BIGINT) AS start, CAST(pos AS BIGINT) AS end, concat(ref, '/', alt) AS allele_string, symbol, ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl",
        ),
        "alphamissense" => Ok(
            "chrom, CAST(pos AS BIGINT) AS start, CAST(pos AS BIGINT) AS end, concat(ref, '/', alt) AS allele_string, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class",
        ),
        "dbnsfp" => Ok(
            "chrom, CAST(pos AS BIGINT) AS start, CAST(pos AS BIGINT) AS end, concat(ref, '/', alt) AS allele_string, sift4g_score, sift4g_pred, polyphen2_hdiv_score, polyphen2_hvar_score, lrt_score, lrt_pred, mutationtaster_score, mutationtaster_pred, fathmm_score, fathmm_pred, provean_score, provean_pred, vest4_score, metasvm_score, metasvm_pred, metalr_score, metalr_pred, revel_score, gerp_rs, phylop100way, phylop30way, phastcons100way, phastcons30way, siphy_29way, cadd_raw, cadd_phred",
        ),
        _ => Err(DataFusionError::Execution(format!(
            "Unknown plugin for fjall build: {plugin_name}"
        ))),
    }
}

async fn discover_chromosomes(
    ctx: &SessionContext,
    table_name: &str,
    chrom_col: &str,
) -> Result<Vec<String>> {
    let query = format!("SELECT DISTINCT {chrom_col} AS chrom FROM {table_name} ORDER BY chrom");
    let df = ctx.sql(&query).await?;
    let batches = df.collect().await?;

    let mut chroms = Vec::new();
    for batch in &batches {
        let col = batch
            .column(0)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringArray>()
            .ok_or_else(|| {
                DataFusionError::Execution("Chromosome column is not StringArray".to_string())
            })?;
        for i in 0..col.len() {
            if !col.is_null(i) {
                chroms.push(col.value(i).to_string());
            }
        }
    }
    Ok(chroms)
}

fn normalize_chrom(chrom: &str) -> String {
    chrom.strip_prefix("chr").unwrap_or(chrom).to_string()
}

async fn write_query_to_parquet(
    ctx: &SessionContext,
    plugin_name: &str,
    query: &str,
    output_path: &str,
) -> Result<usize> {
    let df = ctx.sql(query).await?;
    let mut stream = df.execute_stream().await?;
    let schema = normalize_plugin_schema(&stream.schema());

    let props = plugin_writer_properties(&schema);
    let file = File::create(output_path).map_err(|e| DataFusionError::Execution(format!("{e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| DataFusionError::Execution(format!("{e}")))?;

    let mut rows = 0usize;
    while let Some(batch_result) = stream.next().await {
        let batch = normalize_plugin_batch(plugin_name, &batch_result?)?;
        if batch.num_rows() == 0 {
            continue;
        }
        rows += batch.num_rows();
        writer.write(&batch)?;
    }
    writer.close()?;
    Ok(rows)
}

fn normalize_plugin_batch(plugin_name: &str, batch: &RecordBatch) -> Result<RecordBatch> {
    let batch = if plugin_name == "clinvar" {
        decompose_multi_allelic_batch(batch)?
    } else {
        batch.clone()
    };

    coerce_plugin_batch_types(&batch)
}

fn coerce_plugin_batch_types(batch: &RecordBatch) -> Result<RecordBatch> {
    let Some((pos_idx, _)) = batch.schema().column_with_name("pos") else {
        return Ok(batch.clone());
    };
    if batch.column(pos_idx).data_type() == &DataType::UInt32 {
        return Ok(batch.clone());
    }

    let casted = cast(batch.column(pos_idx), &DataType::UInt32).map_err(|e| {
        DataFusionError::Execution(format!("Failed to cast plugin column 'pos' to UInt32: {e}"))
    })?;

    let mut fields = batch.schema().fields().iter().cloned().collect::<Vec<_>>();
    fields[pos_idx] = Arc::new(Field::new("pos", DataType::UInt32, true));
    let schema = Arc::new(Schema::new(fields));

    let mut columns = batch.columns().to_vec();
    columns[pos_idx] = casted;
    RecordBatch::try_new(schema, columns)
        .map_err(|e| DataFusionError::Execution(format!("Failed to rebuild plugin batch: {e}")))
}

fn normalize_plugin_schema(schema: &SchemaRef) -> SchemaRef {
    let fields = schema
        .fields()
        .iter()
        .map(|field| {
            if field.name() == "pos" && field.data_type() != &DataType::UInt32 {
                Arc::new(Field::new("pos", DataType::UInt32, field.is_nullable()))
            } else {
                field.clone()
            }
        })
        .collect::<Vec<_>>();
    Arc::new(Schema::new(fields))
}

fn decompose_multi_allelic_batch(batch: &RecordBatch) -> Result<RecordBatch> {
    let Some((alt_idx, _)) = batch.schema().column_with_name("alt") else {
        return Ok(batch.clone());
    };
    let alt_array = batch
        .column(alt_idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            DataFusionError::Execution("ClinVar alt column must be StringArray".to_string())
        })?;
    if !(0..alt_array.len())
        .any(|row| !alt_array.is_null(row) && alt_array.value(row).contains([',', '|']))
    {
        return Ok(batch.clone());
    }

    let mut expanded_indices = Vec::new();
    let mut expanded_alts = Vec::new();
    for row in 0..batch.num_rows() {
        if alt_array.is_null(row) {
            expanded_indices.push(row);
            expanded_alts.push(None);
            continue;
        }
        let value = alt_array.value(row);
        let mut seen_split = false;
        for alt in value
            .split([',', '|'])
            .map(str::trim)
            .filter(|alt| !alt.is_empty())
        {
            expanded_indices.push(row);
            expanded_alts.push(Some(alt.to_string()));
            seen_split = true;
        }
        if !seen_split {
            expanded_indices.push(row);
            expanded_alts.push(Some(value.to_string()));
        }
    }

    let schema = batch.schema();
    let mut columns: Vec<ArrayRef> = Vec::with_capacity(batch.num_columns());
    for (col_idx, field) in schema.fields().iter().enumerate() {
        if col_idx == alt_idx {
            let mut builder =
                StringBuilder::with_capacity(expanded_alts.len(), expanded_alts.len() * 8);
            for value in &expanded_alts {
                if let Some(value) = value {
                    builder.append_value(value);
                } else {
                    builder.append_null();
                }
            }
            columns.push(Arc::new(builder.finish()));
            continue;
        }
        columns.push(expand_column(
            batch.column(col_idx),
            field.data_type(),
            &expanded_indices,
        )?);
    }

    RecordBatch::try_new(schema.clone(), columns).map_err(|e| {
        DataFusionError::Execution(format!(
            "Failed to rebuild ClinVar batch after multi-allelic decomposition: {e}"
        ))
    })
}

fn expand_column(
    source: &ArrayRef,
    data_type: &DataType,
    row_indices: &[usize],
) -> Result<ArrayRef> {
    match data_type {
        DataType::Utf8 => {
            let array = source
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "Expected Utf8 column during ClinVar decomposition".into(),
                    )
                })?;
            let mut builder =
                StringBuilder::with_capacity(row_indices.len(), row_indices.len() * 16);
            for &row in row_indices {
                if array.is_null(row) {
                    builder.append_null();
                } else {
                    builder.append_value(array.value(row));
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        DataType::UInt32 => {
            let array = source
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "Expected UInt32 column during ClinVar decomposition".into(),
                    )
                })?;
            let mut builder = UInt32Builder::with_capacity(row_indices.len());
            for &row in row_indices {
                if array.is_null(row) {
                    builder.append_null();
                } else {
                    builder.append_value(array.value(row));
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        DataType::Int32 => {
            let array = source
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "Expected Int32 column during ClinVar decomposition".into(),
                    )
                })?;
            let mut builder = Int32Builder::with_capacity(row_indices.len());
            for &row in row_indices {
                if array.is_null(row) {
                    builder.append_null();
                } else {
                    builder.append_value(array.value(row));
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        DataType::Float32 => {
            let array = source
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "Expected Float32 column during ClinVar decomposition".into(),
                    )
                })?;
            let mut builder = Float32Builder::with_capacity(row_indices.len());
            for &row in row_indices {
                if array.is_null(row) {
                    builder.append_null();
                } else {
                    builder.append_value(array.value(row));
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        other => Err(DataFusionError::Execution(format!(
            "Unsupported ClinVar decomposition column type: {other}"
        ))),
    }
}

fn plugin_writer_properties(schema: &SchemaRef) -> WriterProperties {
    let sort_cols = ["chrom", "pos", "ref", "alt"];
    let sorting = {
        let cols: Vec<SortingColumn> = sort_cols
            .iter()
            .filter_map(|name| {
                schema
                    .column_with_name(name)
                    .map(|(idx, _)| SortingColumn::new(idx as i32, false, false))
            })
            .collect();
        if cols.len() == sort_cols.len() {
            Some(cols)
        } else {
            None
        }
    };

    let mut builder = WriterProperties::builder()
        .set_compression(ParquetCompression::ZSTD(Default::default()))
        .set_max_row_group_size(PLUGIN_ROW_GROUP_SIZE)
        .set_sorting_columns(sorting);

    if schema.column_with_name("chrom").is_some() {
        builder = builder.set_column_bloom_filter_enabled(ColumnPath::from("chrom"), true);
    }
    if schema.column_with_name("pos").is_some() {
        builder = builder.set_column_bloom_filter_enabled(ColumnPath::from("pos"), true);
    }

    builder.build()
}

fn chrom_aliases_for_sql(chrom: &str) -> Vec<String> {
    let normalized = chrom.strip_prefix("chr").unwrap_or(chrom);
    if normalized == chrom {
        vec![
            format!("'{}'", normalized.replace('\'', "''")),
            format!("'chr{}'", normalized.replace('\'', "''")),
        ]
    } else {
        vec![
            format!("'{}'", normalized.replace('\'', "''")),
            format!("'{}'", chrom.replace('\'', "''")),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::{
        PLUGIN_ROW_GROUP_SIZE, convert_plugin_to_parquet, normalize_chrom, plugin_select_query,
        plugin_writer_properties,
    };
    use datafusion::arrow::array::{Float32Array, Int32Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::parquet::basic::Compression;
    use datafusion::parquet::schema::types::ColumnPath;
    use flate2::Compression as GzipCompression;
    use flate2::write::GzEncoder;
    use std::fs::File;
    use std::io::Write;
    use std::sync::Arc;

    fn write_gzip_text(path: &std::path::Path, text: &str) {
        let file = File::create(path).expect("create gzip");
        let mut encoder = GzEncoder::new(file, GzipCompression::default());
        encoder.write_all(text.as_bytes()).expect("write gzip");
        encoder.finish().expect("finish gzip");
    }

    fn write_text(path: &std::path::Path, text: &str) {
        std::fs::write(path, text).expect("write text");
    }

    #[test]
    fn normalize_chrom_strips_prefix_only_when_present() {
        assert_eq!(normalize_chrom("chr7"), "7");
        assert_eq!(normalize_chrom("X"), "X");
    }

    #[test]
    fn spliceai_query_projects_expected_columns() {
        let query = plugin_select_query("spliceai").expect("query");
        assert!(query.contains(" AS ds_ag"));
        assert!(query.contains(" AS symbol"));
    }

    #[test]
    fn plugin_queries_use_datafusion_sql_type_names() {
        for plugin in ["clinvar", "cadd", "alphamissense", "spliceai", "dbnsfp"] {
            let query = plugin_select_query(plugin).expect("query");
            assert!(!query.contains("INT64"));
            assert!(!query.contains(" AS INT)"));
        }
    }

    #[test]
    fn clinvar_query_uses_quoted_vcf_info_fields() {
        let query = plugin_select_query("clinvar").expect("query");
        assert!(query.contains("array_to_string(\"CLNSIG\", '|') AS clnsig"));
        assert!(query.contains("array_to_string(\"CLNDN\", '|') AS clndn"));
        assert!(query.contains("array_to_string(\"CLNREVSTAT\", '|') AS clnrevstat"));
        assert!(query.contains("\"CLNVC\" AS clnvc"));
        assert!(query.contains(" AS clnvi"));
        assert!(query.contains(" AS af_esp"));
        assert!(query.contains(" AS af_exac"));
        assert!(query.contains(" AS af_tgp"));
    }

    #[test]
    fn alphamissense_query_projects_extended_columns() {
        let query = plugin_select_query("alphamissense").expect("query");
        assert!(query.contains("genome"));
        assert!(query.contains("uniprot_id"));
        assert!(query.contains("transcript_id"));
        assert!(query.contains("protein_variant"));
    }

    #[test]
    fn dbnsfp_query_projects_extended_columns() {
        let query = plugin_select_query("dbnsfp").expect("query");
        assert!(query.contains("sift4g_score"));
        assert!(query.contains("polyphen2_hdiv_score"));
        assert!(query.contains("phylop100way"));
        assert!(query.contains("cadd_phred"));
    }

    #[test]
    fn plugin_writer_properties_enable_expected_layout() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("pos", DataType::UInt32, true),
            Field::new("ref", DataType::Utf8, true),
            Field::new("alt", DataType::Utf8, true),
        ]));
        let props = plugin_writer_properties(&schema);
        assert_eq!(props.max_row_group_size(), PLUGIN_ROW_GROUP_SIZE);
        assert_eq!(
            props.compression(&ColumnPath::from("chrom")),
            Compression::ZSTD(Default::default())
        );
        assert!(props.sorting_columns().is_some());
    }

    #[test]
    fn alphamissense_round_trip_preserves_extended_schema() {
        let temp = tempfile::tempdir().expect("tempdir");
        let source = temp.path().join("alphamissense.tsv.gz");
        write_gzip_text(
            &source,
            concat!(
                "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n",
                "X\t200\tC\tT\tGRCh38\tP2\tENST2\tp.C2T\t0.10\tbenign\n",
                "X\t100\tA\tG\tGRCh38\tP1\tENST1\tp.A1G\t0.42\tlikely_pathogenic\n",
            ),
        );

        let output = temp.path().join("out");
        std::fs::create_dir_all(&output).expect("output dir");
        let results = convert_plugin_to_parquet(
            "alphamissense",
            source.to_str().unwrap(),
            output.to_str().unwrap(),
            1,
            32,
            None,
        )
        .expect("convert");
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1, 2);

        let runtime = tokio::runtime::Runtime::new().expect("runtime");
        let batches = runtime.block_on(async {
            let ctx = datafusion::prelude::SessionContext::new();
            let df = ctx
                .read_parquet(results[0].0.as_str(), Default::default())
                .await
                .expect("read parquet");
            df.collect().await.expect("collect")
        });
        let batch = &batches[0];
        assert!(batch.schema().index_of("genome").is_ok());
        assert!(batch.schema().index_of("uniprot_id").is_ok());
        assert!(batch.schema().index_of("transcript_id").is_ok());
        assert!(batch.schema().index_of("protein_variant").is_ok());
        let pos = batch
            .column(batch.schema().index_of("pos").unwrap())
            .as_any()
            .downcast_ref::<UInt32Array>()
            .expect("pos");
        assert_eq!(pos.values(), &[100, 200]);
    }

    #[test]
    fn spliceai_round_trip_preserves_columns() {
        let temp = tempfile::tempdir().expect("tempdir");
        let source = temp.path().join("spliceai.vcf.gz");
        write_gzip_text(
            &source,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI composite field\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "1\t100\t.\tA\tG\t.\t.\tSpliceAI=A|GENE1|0.71|0.12|0.33|0.01|1|2|3|4\n",
            ),
        );

        let output = temp.path().join("out");
        std::fs::create_dir_all(&output).expect("output dir");
        let results = convert_plugin_to_parquet(
            "spliceai",
            source.to_str().unwrap(),
            output.to_str().unwrap(),
            1,
            32,
            None,
        )
        .expect("convert");
        assert_eq!(results[0].1, 1);

        let runtime = tokio::runtime::Runtime::new().expect("runtime");
        let batches = runtime.block_on(async {
            let ctx = datafusion::prelude::SessionContext::new();
            let df = ctx
                .read_parquet(results[0].0.as_str(), Default::default())
                .await
                .expect("read parquet");
            df.collect().await.expect("collect")
        });
        let batch = &batches[0];
        for column in [
            "symbol",
            "ds_ag",
            "ds_al",
            "ds_dg",
            "ds_dl",
            "dp_ag",
            "dp_al",
            "dp_dg",
            "dp_dl",
        ] {
            assert!(batch.schema().index_of(column).is_ok(), "{column}");
        }

        let ds_ag = batch
            .column(batch.schema().index_of("ds_ag").unwrap())
            .as_any()
            .downcast_ref::<Float32Array>()
            .expect("ds_ag array");
        assert_eq!(ds_ag.value(0), 0.71);

        let dp_ag = batch
            .column(batch.schema().index_of("dp_ag").unwrap())
            .as_any()
            .downcast_ref::<Int32Array>()
            .expect("dp_ag array");
        assert_eq!(dp_ag.value(0), 1);
    }

    #[test]
    fn dbnsfp_round_trip_preserves_extended_schema() {
        let temp = tempfile::tempdir().expect("tempdir");
        let source = temp.path().join("dbnsfp.tsv.gz");
        write_gzip_text(
            &source,
            concat!(
                "#chr\tpos(1-based)\tref\talt\tSIFT4G_score\tSIFT4G_pred\tPolyphen2_HDIV_score\tPolyphen2_HVAR_score\tLRT_score\tLRT_pred\tMutationTaster_score\tREVEL_score\tMetaSVM_score\tMetaSVM_pred\tMetaLR_score\tMetaLR_pred\tGERP++_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphastCons100way_vertebrate\tphastCons30way_mammalian\tSiPhy_29way_logOdds\tCADD_raw\tCADD_phred\tFATHMM_score\tFATHMM_pred\tPROVEAN_score\tPROVEAN_pred\tVEST4_score\tBayesDel_addAF_score\tBayesDel_noAF_score\tMutationTaster_pred\n",
                "X\t100\tA\tG\t0.1\tT\t0.2\t0.3\t0.4\tD\t0.5\t0.6\t0.7\tD\t0.8\tD\t1.0\t1.1\t1.2\t1.3\t1.4\t1.5\t1.6\t10.0\t0.9\tD\tN\tD\t0.7\t0.11\t0.12\tD\n",
            ),
        );

        let output = temp.path().join("out");
        std::fs::create_dir_all(&output).expect("output dir");
        let results = convert_plugin_to_parquet(
            "dbnsfp",
            source.to_str().unwrap(),
            output.to_str().unwrap(),
            1,
            32,
            None,
        )
        .expect("convert");
        assert_eq!(results[0].1, 1);

        let runtime = tokio::runtime::Runtime::new().expect("runtime");
        let batches = runtime.block_on(async {
            let ctx = datafusion::prelude::SessionContext::new();
            let df = ctx
                .read_parquet(results[0].0.as_str(), Default::default())
                .await
                .expect("read parquet");
            df.collect().await.expect("collect")
        });
        let batch = &batches[0];
        for column in [
            "sift4g_score",
            "polyphen2_hdiv_score",
            "revel_score",
            "phylop100way",
            "cadd_phred",
        ] {
            assert!(batch.schema().index_of(column).is_ok(), "{column}");
        }
    }

    #[test]
    fn clinvar_round_trip_includes_extended_columns() {
        let temp = tempfile::tempdir().expect("tempdir");
        let source = temp.path().join("clinvar.vcf.gz");
        write_gzip_text(
            &source,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">\n",
                "##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"Disease names\">\n",
                "##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description=\"Review status\">\n",
                "##INFO=<ID=CLNVC,Number=1,Type=String,Description=\"Variant type\">\n",
                "##INFO=<ID=CLNVI,Number=.,Type=String,Description=\"Variant identifiers\">\n",
                "##INFO=<ID=AF_ESP,Number=1,Type=Float,Description=\"ESP AF\">\n",
                "##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description=\"ExAC AF\">\n",
                "##INFO=<ID=AF_TGP,Number=1,Type=Float,Description=\"TGP AF\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "1\t100\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNDN=Disease;CLNREVSTAT=criteria_provided;CLNVC=single_nucleotide_variant;CLNVI=RCV:1;AF_ESP=0.1;AF_EXAC=0.2;AF_TGP=0.3\n",
            ),
        );

        let output = temp.path().join("out");
        std::fs::create_dir_all(&output).expect("output dir");
        let results = convert_plugin_to_parquet(
            "clinvar",
            source.to_str().unwrap(),
            output.to_str().unwrap(),
            1,
            32,
            None,
        )
        .expect("convert");
        assert_eq!(results[0].1, 1);

        let runtime = tokio::runtime::Runtime::new().expect("runtime");
        let batches = runtime.block_on(async {
            let ctx = datafusion::prelude::SessionContext::new();
            let df = ctx
                .read_parquet(results[0].0.as_str(), Default::default())
                .await
                .expect("read parquet");
            df.collect().await.expect("collect")
        });
        let batch = &batches[0];
        for column in [
            "clnsig",
            "clndn",
            "clnrevstat",
            "clnvc",
            "clnvi",
            "af_esp",
            "af_exac",
            "af_tgp",
        ] {
            assert!(batch.schema().index_of(column).is_ok(), "{column}");
        }
    }

    #[test]
    fn clinvar_round_trip_decomposes_multi_allelic_rows() {
        let temp = tempfile::tempdir().expect("tempdir");
        let source = temp.path().join("clinvar_multi.vcf.gz");
        write_gzip_text(
            &source,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">\n",
                "##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"Disease names\">\n",
                "##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description=\"Review status\">\n",
                "##INFO=<ID=CLNVC,Number=1,Type=String,Description=\"Variant type\">\n",
                "##INFO=<ID=CLNVI,Number=.,Type=String,Description=\"Variant identifiers\">\n",
                "##INFO=<ID=AF_ESP,Number=1,Type=Float,Description=\"ESP AF\">\n",
                "##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description=\"ExAC AF\">\n",
                "##INFO=<ID=AF_TGP,Number=1,Type=Float,Description=\"TGP AF\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "1\t100\t.\tA\tG,T\t.\t.\tCLNSIG=Pathogenic;CLNDN=Disease;CLNREVSTAT=criteria_provided;CLNVC=single_nucleotide_variant;CLNVI=RCV:1;AF_ESP=0.1;AF_EXAC=0.2;AF_TGP=0.3\n",
            ),
        );

        let output = temp.path().join("out");
        std::fs::create_dir_all(&output).expect("output dir");
        let results = convert_plugin_to_parquet(
            "clinvar",
            source.to_str().unwrap(),
            output.to_str().unwrap(),
            1,
            32,
            None,
        )
        .expect("convert");
        assert_eq!(results[0].1, 2);

        let runtime = tokio::runtime::Runtime::new().expect("runtime");
        let batches = runtime.block_on(async {
            let ctx = datafusion::prelude::SessionContext::new();
            let df = ctx
                .read_parquet(results[0].0.as_str(), Default::default())
                .await
                .expect("read parquet");
            df.collect().await.expect("collect")
        });
        let batch = &batches[0];
        let alt_idx = batch.schema().index_of("alt").unwrap();
        let alt_col = batch.column(alt_idx);
        let alt0 = datafusion::arrow::util::display::array_value_to_string(alt_col.as_ref(), 0)
            .expect("alt0");
        let alt1 = datafusion::arrow::util::display::array_value_to_string(alt_col.as_ref(), 1)
            .expect("alt1");
        assert_eq!(alt0, "G");
        assert_eq!(alt1, "T");
    }
}
