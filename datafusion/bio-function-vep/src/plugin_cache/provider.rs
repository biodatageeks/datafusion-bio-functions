//! Provider factory: register a source manifest's raw tables under their
//! `plugin_<name>_src[_<part>]` names. CSV/TSV/Parquet use builtin DataFusion
//! providers; VCF uses bio-formats' `VcfTableProvider`. BED is not wired yet
//! (no plugin needs it).

use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{CsvReadOptions, ParquetReadOptions, SessionContext};
use datafusion_bio_format_core::object_storage::ObjectStorageOptions;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;

use crate::plugin_cache::source_manifest::{CsvParams, ProviderKind, SourceManifest, ValueType};

fn arrow_type(ty: ValueType) -> DataType {
    match ty {
        ValueType::Utf8 => DataType::Utf8,
        ValueType::Float32 => DataType::Float32,
        ValueType::Int32 => DataType::Int32,
    }
}

fn csv_schema(csv: &CsvParams) -> Schema {
    Schema::new(
        csv.schema
            .iter()
            .map(|f| Field::new(&f.name, arrow_type(f.ty), true))
            .collect::<Vec<_>>(),
    )
}

/// Materialize a plain (uncompressed) CSV/TSV path DataFusion can read.
///
/// The workspace builds DataFusion with `default-features = false` and no
/// `compression` feature (the `xz2`/liblzma link collision with noodles-cram —
/// see the root `Cargo.toml`), so `register_csv` cannot decompress input files.
/// For gzip sources we stream-decompress to a temp file and register that; plain
/// sources pass through unchanged. Returns `(path, Some(temp_path))` for gzip:
/// the caller must keep the `TempPath` alive for the duration of query execution
/// and drop it afterwards, which deletes the temp — so a multi-chrom `build_all`
/// keeps at most one decompressed copy on disk at a time (rather than leaking one
/// per chromosome). Plain sources return `(path, None)`.
fn materialize_plain(path: &str, gzip: bool) -> Result<(String, Option<tempfile::TempPath>)> {
    if !gzip {
        return Ok((path.to_string(), None));
    }
    let src = std::fs::File::open(path)
        .map_err(|e| DataFusionError::Execution(format!("open gzip source '{path}': {e}")))?;
    let mut decoder = flate2::read::MultiGzDecoder::new(src);
    let mut tmp = tempfile::Builder::new()
        .prefix("plugin_src_")
        .suffix(".tsv")
        .tempfile()
        .map_err(|e| DataFusionError::Execution(format!("create temp for '{path}': {e}")))?;
    // Stream-decompress straight to the temp file. These sources can be tens of
    // GB uncompressed (CADD/dbNSFP), so buffering the whole thing in memory
    // (`read_to_end`) would OOM the build.
    std::io::copy(&mut decoder, tmp.as_file_mut())
        .map_err(|e| DataFusionError::Execution(format!("decompress '{path}': {e}")))?;
    // Return the temp handle so the caller can keep it alive only for this
    // chrom's build and drop (delete) it afterwards — leaking it (`keep()`)
    // would accumulate one tens-of-GB temp per chromosome across `build_all`.
    let temp_path = tmp.into_temp_path();
    let plain = temp_path.to_string_lossy().into_owned();
    Ok((plain, Some(temp_path)))
}

/// Register every source in `manifest` as a DataFusion table. Returns any
/// decompressed temp files the caller must keep alive for the duration of query
/// execution (dropping them deletes the temp — see `materialize_plain`).
pub async fn register_sources(
    ctx: &SessionContext,
    manifest: &SourceManifest,
) -> Result<Vec<tempfile::TempPath>> {
    let mut temps: Vec<tempfile::TempPath> = Vec::new();
    for spec in &manifest.sources {
        let table = spec.table_name(&manifest.plugin_name);
        match spec.provider {
            ProviderKind::Csv | ProviderKind::Tsv => {
                let csv = spec.csv.as_ref().ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "csv/tsv source '{table}' missing [source.csv]"
                    ))
                })?;
                let schema = csv_schema(csv);
                let delim = csv.delimiter.as_bytes().first().copied().unwrap_or(b'\t');
                let gz = csv.compression.as_deref() == Some("gzip");
                let (plain_path, temp) = materialize_plain(&spec.path, gz)?;
                if let Some(t) = temp {
                    temps.push(t);
                }
                let ext = std::path::Path::new(&plain_path)
                    .extension()
                    .and_then(|e| e.to_str())
                    .map(|e| format!(".{e}"))
                    .unwrap_or_default();
                let mut opts = CsvReadOptions::new()
                    .delimiter(delim)
                    .has_header(csv.has_header)
                    .schema(&schema)
                    .file_extension(&ext);
                if let Some(c) = csv
                    .comment
                    .as_ref()
                    .and_then(|s| s.as_bytes().first().copied())
                {
                    opts = opts.comment(c);
                }
                ctx.register_csv(&table, &plain_path, opts).await?;
            }
            ProviderKind::Parquet => {
                ctx.register_parquet(&table, &spec.path, ParquetReadOptions::default())
                    .await?;
            }
            ProviderKind::Vcf => {
                let info_fields = spec.vcf.as_ref().and_then(|v| v.info_fields.clone());
                // ObjectStorageOptions::default() sets compression_type = AUTO, which is what
                // lets one code path read both plain `.vcf` and BGZF `.vcf.gz`. `None` would
                // behave identically today (the reader `unwrap_or_default()`s it at every use),
                // but passing the value explicitly keeps the compression policy visible here.
                //
                // coordinate_system_zero_based = false: VCF POS is 1-based and the plugin cache
                // stores 1-based start/end, so the reader must not shift. The manifest's
                // `coordinate_system` remains the single source of truth for any shift the
                // builder applies (see plugin_cache::build::wrap_normalization).
                let path = &spec.path;
                let vcf_table = VcfTableProvider::new(
                    spec.path.clone(),
                    info_fields,
                    None,
                    Some(ObjectStorageOptions::default()),
                    false,
                )
                .map_err(|e| {
                    DataFusionError::Execution(format!(
                        "open vcf source '{path}' for table '{table}': {e}"
                    ))
                })?;
                ctx.register_table(table.as_str(), Arc::new(vcf_table))
                    .map_err(|e| {
                        DataFusionError::Execution(format!(
                            "register vcf table '{table}' ('{path}'): {e}"
                        ))
                    })?;
            }
            ProviderKind::Bed => {
                return Err(DataFusionError::NotImplemented(format!(
                    "provider 'bed' is not wired yet (table '{table}'); \
                     no plugin needs it — wire datafusion-bio-format-bed the way \
                     ProviderKind::Vcf is wired when one does"
                )));
            }
        }
    }
    Ok(temps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::SourceManifest;
    use datafusion::arrow::array::{Array, Int64Array};
    use std::io::Write;

    fn write_gz(path: &std::path::Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn registers_csv_source_with_declared_schema() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("snv.tsv.gz");
        // headerless, comment line, two data rows
        write_gz(
            &tsv,
            "#comment\nchr1\t100\tA\tG\t0.9\nchr1\t200\tC\tT\t0.1\n",
        );

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  comment = "#"
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
"##,
            tsv.display()
        );

        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let ctx = SessionContext::new();
        // Keep the temp alive for the query below (drop deletes it).
        let _temps = register_sources(&ctx, &manifest).await.unwrap();
        let n = ctx
            .sql("SELECT count(*) AS c FROM plugin_demo_src")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        let c = n[0]
            .column(0)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap()
            .value(0);
        assert_eq!(c, 2);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn registers_vcf_source_and_projects_info_fields() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("demo.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr22\t22893742\t.\tC\tG\t.\t.\tSCORE=0.9\n\
             chr22\t22893800\t.\tA\tT\t.\t.\tSCORE=0.1\n",
        )
        .unwrap();

        let toml_src = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "{}"
  [source.vcf]
  info_fields = ["SCORE"]

[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
"##,
            vcf.display()
        );

        let manifest: SourceManifest = toml::from_str(&toml_src).unwrap();
        let ctx = SessionContext::new();
        let _temps = register_sources(&ctx, &manifest).await.unwrap();

        // The reader exposes VCF POS as `start` (plus a matching `end`), NOT `pos`.
        // INFO columns are bare, case-sensitive keys -> must be backticked.
        let batches = ctx
            .sql("SELECT chrom, `start`, `SCORE` FROM plugin_demo_src ORDER BY `start`")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();

        let rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(rows, 2, "both VCF records must be visible");

        // The reader must report 1-based POS (22893742), matching the cache's 1-based start.
        let pos = batches[0]
            .column(1)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::UInt32Array>()
            .expect("start is UInt32");
        assert_eq!(pos.value(0), 22_893_742);

        // The selected INFO field must actually carry its parsed value — a column that
        // resolved by name but came back all-null would otherwise pass unnoticed.
        let score = batches[0]
            .column(2)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Float32Array>()
            .expect("SCORE is Float32");
        assert!(score.is_valid(0), "SCORE must be parsed, not null");
        assert!(
            (score.value(0) - 0.9).abs() < 1e-6,
            "SCORE = {}, want 0.9",
            score.value(0)
        );
    }
}
