//! Plugin cache builder: converts external plugin source files to parquet and fjall.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{Array, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::parquet::arrow::ArrowWriter;
use datafusion::parquet::basic::Compression as ParquetCompression;
use datafusion::parquet::file::properties::WriterProperties;
use datafusion::parquet::format::SortingColumn;
use datafusion::parquet::schema::types::ColumnPath;
use datafusion::prelude::{CsvReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vep_plugin::{
    PluginSourceKind, build_where_clause as plugin_build_where_clause, cadd_union_query,
    normalize_plugin_batch as format_normalize_plugin_batch,
    normalize_plugin_schema as format_normalize_plugin_schema, register_plugin_source,
};
use futures::StreamExt;

use crate::kv_cache::CacheLoader;
use crate::plugin::PluginKind;

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
    assume_sorted_input: bool,
    preview_rows: Option<usize>,
) -> Result<Vec<(String, usize)>> {
    match plugin {
        PluginKind::ClinVar => convert_plugin_to_parquet(
            "clinvar",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
            assume_sorted_input,
            preview_rows,
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
                assume_sorted_input,
                preview_rows,
            )
        }
        PluginKind::SpliceAI => convert_plugin_to_parquet(
            "spliceai",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
            assume_sorted_input,
            preview_rows,
        ),
        PluginKind::AlphaMissense => convert_plugin_to_parquet(
            "alphamissense",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
            assume_sorted_input,
            preview_rows,
        ),
        PluginKind::DbNSFP => convert_plugin_to_parquet(
            "dbnsfp",
            source_path,
            output_dir,
            partitions,
            memory_limit_gb,
            None,
            assume_sorted_input,
            preview_rows,
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
    assume_sorted_input: bool,
    preview_rows: Option<usize>,
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
        assume_sorted_input,
        preview_rows,
    ))
}

pub fn convert_cadd_sources_to_parquet(
    snv_source_path: &str,
    indel_source_path: &str,
    output_dir: &str,
    partitions: usize,
    memory_limit_gb: usize,
    chromosomes: Option<Vec<String>>,
    assume_sorted_input: bool,
    preview_rows: Option<usize>,
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
        assume_sorted_input,
        preview_rows,
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
    assume_sorted_input: bool,
    preview_rows: Option<usize>,
) -> Result<Vec<(String, usize)>> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);
    let plugin_kind = PluginSourceKind::from_name(plugin_name)?;
    register_plugin_source(&ctx, plugin_kind, source_path, "source").await?;

    let select_query = plugin_kind.select_query();
    let chrom_col = plugin_kind.chrom_column();
    let requested: Option<HashSet<String>> = chromosomes.map(|values| {
        values
            .into_iter()
            .map(|value| normalize_chrom(&value))
            .collect()
    });
    let grouped = if let Some(requested_values) = requested.as_ref() {
        group_requested_chromosomes(requested_values)
    } else {
        let chroms = discover_chromosomes(&ctx, "source", chrom_col).await?;
        group_chromosomes(chroms, None)
    };
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
        let where_clause = plugin_build_where_clause(chrom_col, raw_chroms);
        let query = format!(
            "{select_query}{where_clause}{}{}",
            plugin_order_clause(plugin_name, assume_sorted_input),
            plugin_limit_clause(preview_rows)
        );
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
        let where_clause = plugin_build_where_clause(chrom_col, &all_other_raw);
        let query = format!(
            "{select_query}{where_clause}{}{}",
            plugin_order_clause(plugin_name, assume_sorted_input),
            plugin_limit_clause(preview_rows)
        );
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
    assume_sorted_input: bool,
    preview_rows: Option<usize>,
) -> Result<Vec<(String, usize)>> {
    let config = SessionConfig::new().with_target_partitions(partitions);
    let ctx = SessionContext::new_with_config(config);
    register_plugin_source(&ctx, PluginSourceKind::Cadd, snv_source_path, "source_snv").await?;
    register_plugin_source(
        &ctx,
        PluginSourceKind::Cadd,
        indel_source_path,
        "source_indel",
    )
    .await?;

    let requested: Option<HashSet<String>> = chromosomes.map(|values| {
        values
            .into_iter()
            .map(|value| normalize_chrom(&value))
            .collect()
    });
    let grouped = if let Some(requested_values) = requested.as_ref() {
        group_requested_chromosomes(requested_values)
    } else {
        let mut chroms = discover_chromosomes(&ctx, "source_snv", "chrom").await?;
        chroms.extend(discover_chromosomes(&ctx, "source_indel", "chrom").await?);
        group_chromosomes(chroms, None)
    };
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
    if assume_sorted_input {
        eprintln!("  cadd: assume_sorted_input ignored because CADD merges SNV + indel sources");
    }

    let mut all_results = Vec::new();
    let global_start = Instant::now();
    let mut total_rows = 0usize;

    for (norm, raw_chroms) in &main_chroms {
        let query = cadd_union_query(raw_chroms, preview_rows);
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
        let query = cadd_union_query(&all_other_raw, preview_rows);
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

    let select_expr = PluginSourceKind::from_name(plugin_name)?.fjall_projection();
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

fn group_requested_chromosomes(requested: &HashSet<String>) -> Vec<(String, Vec<String>)> {
    let mut grouped: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for chrom in requested {
        let canonical = normalize_chrom(chrom);
        for alias in chrom_aliases(&canonical) {
            grouped.entry(canonical.clone()).or_default().insert(alias);
        }
    }
    grouped
        .into_iter()
        .map(|(canonical, raw_values)| (canonical, raw_values.into_iter().collect()))
        .collect()
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

fn plugin_chrom_column(plugin_name: &str) -> &'static str {
    PluginSourceKind::from_name(plugin_name)
        .map(PluginSourceKind::chrom_column)
        .unwrap_or("chrom")
}

fn plugin_select_query(plugin_name: &str) -> Result<String> {
    Ok(PluginSourceKind::from_name(plugin_name)?
        .select_query()
        .to_string())
}

fn plugin_fjall_projection(plugin_name: &str) -> Result<&'static str> {
    Ok(PluginSourceKind::from_name(plugin_name)?.fjall_projection())
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

fn plugin_order_clause(plugin_name: &str, assume_sorted_input: bool) -> &'static str {
    if assume_sorted_input && plugin_name != "cadd" {
        ""
    } else {
        " ORDER BY chrom, pos, ref, alt"
    }
}

fn plugin_limit_clause(preview_rows: Option<usize>) -> String {
    preview_rows
        .map(|limit| format!(" LIMIT {limit}"))
        .unwrap_or_default()
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
    format_normalize_plugin_batch(PluginSourceKind::from_name(plugin_name)?, batch)
}

fn coerce_plugin_batch_types(batch: &RecordBatch) -> Result<RecordBatch> {
    let kind = if batch.schema().column_with_name("clnsig").is_some() {
        PluginSourceKind::ClinVar
    } else {
        return Ok(batch.clone());
    };
    format_normalize_plugin_batch(kind, batch)
}

fn normalize_plugin_schema(schema: &SchemaRef) -> SchemaRef {
    format_normalize_plugin_schema(schema)
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

fn chrom_aliases(chrom: &str) -> Vec<String> {
    let normalized = chrom.strip_prefix("chr").unwrap_or(chrom);
    if normalized == chrom {
        vec![normalized.to_string(), format!("chr{normalized}")]
    } else {
        vec![normalized.to_string(), chrom.to_string()]
    }
}

#[cfg(test)]
mod tests {
    use super::{
        PLUGIN_ROW_GROUP_SIZE, convert_plugin_to_parquet, group_requested_chromosomes,
        normalize_chrom, plugin_select_query, plugin_writer_properties,
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
    fn requested_chromosomes_expand_chr_aliases_without_scan() {
        let grouped = group_requested_chromosomes(
            &["Y".to_string(), "MT".to_string()].into_iter().collect(),
        );
        assert!(grouped.contains(&(String::from("MT"), vec![String::from("MT"), String::from("chrMT")])));
        assert!(grouped.contains(&(String::from("Y"), vec![String::from("Y"), String::from("chrY")])));
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
            false,
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
            false,
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
            "symbol", "ds_ag", "ds_al", "ds_dg", "ds_dl", "dp_ag", "dp_al", "dp_dg", "dp_dl",
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
                "#chr\tpos(1-based)\tref\talt\tSIFT4G_score\tSIFT4G_pred\tPolyphen2_HDIV_score\tPolyphen2_HVAR_score\tMutationTaster_score\tMutationTaster_pred\tPROVEAN_score\tPROVEAN_pred\tVEST4_score\tMetaSVM_score\tMetaSVM_pred\tMetaLR_score\tMetaLR_pred\tREVEL_score\tGERP++_RS\tphyloP100way_vertebrate\tphastCons100way_vertebrate\tCADD_raw\tCADD_phred\n",
                "X\t100\tA\tG\t0.1\tT\t0.2\t0.3\t0.5\tD\tN\tD\t0.7\t0.8\tD\t0.9\tD\t0.6\t1.0\t1.1\t1.3\t1.6\t10.0\n",
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
            false,
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
            "lrt_score",
            "fathmm_score",
            "phylop30way",
            "phastcons30way",
            "siphy_29way",
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
            false,
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
            false,
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
