//! Provider factory: register a source manifest's raw tables under their
//! `plugin_<name>_src[_<part>]` names. CSV/TSV/Parquet use builtin DataFusion
//! providers; VCF uses `datafusion-bio-format-vcf`'s `VcfTableProvider` (INFO
//! fields come back as typed top-level columns named after the raw INFO tag —
//! query them directly in `ingest_sql`); BED uses `datafusion-bio-format-bed`'s
//! `BedTableProvider`, whose schema is only ever `chrom, start, end, name`
//! regardless of the declared column count (`BEDFields` selects how many
//! columns the reader parses off each line, not how many it exposes) — a
//! source needing more than one extra field packs it into `name` and splits
//! it back out in `ingest_sql`, the same trick SpliceAI's flattened INFO tag
//! uses.

use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{CsvReadOptions, ParquetReadOptions, SessionContext};
use datafusion_bio_format_bed::table_provider::{BEDFields, BedTableProvider};
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
                let zero_based = matches!(
                    manifest.coordinate_system,
                    CoordinateSystem::ZeroBasedHalfOpen
                );
                // No manifest knob yet for selecting a INFO/FORMAT subset — read
                // every INFO field the VCF header declares (`None`) and let
                // `ingest_sql` project down to what it needs, same as the CSV
                // path relies on `ingest_sql` rather than a narrower provider
                // schema. Revisit if a source has so many INFO fields that
                // schema inference itself becomes the bottleneck.
                let provider =
                    VcfTableProvider::new(spec.path.clone(), None, None, None, zero_based)
                        .map_err(|e| {
                            DataFusionError::Execution(format!("open VCF source '{table}': {e}"))
                        })?;
                ctx.register_table(&table, Arc::new(provider))?;
            }
            ProviderKind::Bed => {
                // Same reasoning as the VCF branch above: the manifest's
                // `coordinate_system` drives what `ingest_sql` expects, so
                // keep the provider's own zero-based flag in lock-step with
                // it rather than trusting the file's own convention.
                let zero_based = matches!(
                    manifest.coordinate_system,
                    CoordinateSystem::ZeroBasedHalfOpen
                );
                // No manifest knob yet for BED4/5/6 -- `determine_schema` in
                // `datafusion-bio-format-bed` only ever exposes `chrom, start,
                // end, name` regardless of variant, so BED4 is the only
                // choice that matters until a source needs more raw columns
                // parsed off each line than that.
                let provider =
                    BedTableProvider::new(spec.path.clone(), BEDFields::BED4, None, zero_based)
                        .map_err(|e| {
                            DataFusionError::Execution(format!("open BED source '{table}': {e}"))
                        })?;
                ctx.register_table(&table, Arc::new(provider))?;
            }
        }
    }
    Ok(temps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::SourceManifest;
    use datafusion::arrow::array::{Int64Array, StringArray};
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

    // I4: document (and pin, so a `datafusion-bio-format-vcf` version bump that
    // changes this can't slip by silently) how the VCF provider actually
    // reports a multi-allelic ALT. It comes back as one string joined per-
    // allele (confirmed below: '|', not the VCF spec's own ','), not split
    // into separate rows -- an `ingest_sql` that assumes the wrong separator
    // (or assumes it's split already) builds an `allele_string` that never
    // matches the per-allele variation key, i.e. a silent miss on every
    // multi-allelic site. This test exists so that would fail loudly here
    // instead of being discovered downstream.
    #[tokio::test(flavor = "multi_thread")]
    async fn vcf_source_multiallelic_alt_shape_is_pinned() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("multiallelic.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t100\t.\tG\tA,C\t.\t.\t.\n",
        )
        .unwrap();

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"
value_columns = []

[[source]]
provider = "vcf"
path = "{}"
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let ctx = SessionContext::new();
        let _temps = register_sources(&ctx, &manifest).await.unwrap();
        let rows = ctx
            .sql("SELECT alt FROM plugin_demo_src")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        assert_eq!(rows.len(), 1, "one row per VCF line, not per allele");
        let alt = rows[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("alt column is a plain string, not a list -- if this ever fails, the ALT shape changed and every manifest built on this provider needs re-auditing");
        // Pinned current behavior: `datafusion-bio-format-vcf` v1.8.8 joins a
        // multi-allelic ALT with '|' (NOT ',' -- the naive guess from the VCF
        // spec's own comma-separated ALT column would be wrong here). Any
        // manifest built directly on `ProviderKind::Vcf` for a source that can
        // be multi-allelic MUST split `alt` on '|' in its own `ingest_sql`
        // (the way the CADD/ClinVar bcftools-flatten preprocessing already
        // splits multi-value INFO tags for VCF sources that go through the
        // CSV path instead).
        assert_eq!(
            alt.value(0),
            "A|C",
            "ALT shape changed -- re-audit every VCF-provider manifest's \
             ingest_sql for correct multi-allelic splitting"
        );
    }
}
