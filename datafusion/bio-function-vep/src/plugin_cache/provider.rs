//! Provider factory: register a source manifest's raw tables under their
//! `plugin_<name>_src[_<part>]` names. CSV/TSV/Parquet use builtin DataFusion
//! providers; VCF uses `datafusion-bio-format-vcf`'s `VcfTableProvider` (INFO
//! fields come back as typed top-level columns named after the raw INFO tag —
//! query them directly in `ingest_sql`). BED is not wired in the prototype.

use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{CsvReadOptions, ParquetReadOptions, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;

use crate::plugin_cache::source_manifest::{
    CoordinateSystem, CsvParams, ProviderKind, SourceManifest, ValueType,
};

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
                // `coordinate_system` describes how the manifest's `ingest_sql`
                // should interpret positions; the VCF provider's own
                // `coordinate_system_zero_based` flag controls how IT reports
                // `start`/`end` from the file's native 1-based VCF `POS` — keep
                // them in lock-step so `ingest_sql` always sees the coordinate
                // system the manifest declares, regardless of source format.
                let zero_based = matches!(manifest.coordinate_system, CoordinateSystem::ZeroBasedHalfOpen);
                // No manifest knob yet for selecting a INFO/FORMAT subset — read
                // every INFO field the VCF header declares (`None`) and let
                // `ingest_sql` project down to what it needs, same as the CSV
                // path relies on `ingest_sql` rather than a narrower provider
                // schema. Revisit if a source has so many INFO fields that
                // schema inference itself becomes the bottleneck.
                let provider = VcfTableProvider::new(spec.path.clone(), None, None, None, zero_based)
                    .map_err(|e| {
                        DataFusionError::Execution(format!("open VCF source '{table}': {e}"))
                    })?;
                ctx.register_table(&table, Arc::new(provider))?;
            }
            ProviderKind::Bed => {
                return Err(DataFusionError::NotImplemented(format!(
                    "provider {:?} not wired in prototype (table '{table}')",
                    spec.provider
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
    use datafusion::arrow::array::Int64Array;
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

[tier]
threshold = 0.01
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
}
