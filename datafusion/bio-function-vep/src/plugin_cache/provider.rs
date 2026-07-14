//! Provider factory: register a source manifest's raw tables under their
//! `plugin_<name>_src[_<part>]` names. CSV/TSV/Parquet use builtin DataFusion
//! providers; VCF uses bio-formats' `VcfTableProvider`. BED is not wired yet
//! (no plugin needs it).

use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{CsvReadOptions, ParquetReadOptions, SessionContext};
use datafusion_bio_format_core::object_storage::ObjectStorageOptions;
use datafusion_bio_format_vcf::storage::get_header;
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

/// Reject `info_fields` keys the VCF header does not declare, BEFORE they reach
/// `VcfTableProvider::new`.
///
/// The reader resolves each requested INFO key with `header.infos().get(tag).unwrap()`
/// (`table_provider.rs:174`), so an unknown key does not return an `Err` — it PANICS:
///
/// ```text
/// thread panicked at datafusion-bio-format-vcf/src/table_provider.rs:174:51:
/// called `Option::unwrap()` on a `None` value
/// ```
///
/// naming no plugin, no manifest, no key and no file. (Which is why the `.map_err(…)`
/// wrapped around `VcfTableProvider::new` below cannot help: a panic is not an `Err`.)
/// And this is the single most likely authoring error, because the reader exposes INFO
/// keys as bare, CASE-SENSITIVE column names: `info_fields = ["score"]` for a header
/// that declares `SCORE` is enough.
///
/// So diff the request against the header first and fail with a diagnostic that names
/// the plugin, the file, every offending key, and the keys the header actually declares.
/// `get_header` is the same entry point the reader itself uses, so this handles plain
/// `.vcf`, BGZF `.vcf.gz` and remote object stores identically to the read path.
///
/// `info_fields = None` (take every declared key) and `info_fields = []` (take none)
/// both request nothing that could be absent, and are left alone.
async fn reject_undeclared_info_fields(
    plugin: &str,
    table: &str,
    path: &str,
    requested: &[String],
) -> Result<()> {
    if requested.is_empty() {
        return Ok(());
    }
    let header = get_header(path.to_string(), Some(ObjectStorageOptions::default()))
        .await
        .map_err(|e| {
            DataFusionError::Execution(format!(
                "plugin '{plugin}': read the VCF header of '{path}' (table '{table}'): {e}"
            ))
        })?;
    let declared: Vec<String> = header.infos().keys().map(ToString::to_string).collect();
    let missing: Vec<&str> = requested
        .iter()
        .map(String::as_str)
        .filter(|key| !declared.iter().any(|d| d == key))
        .collect();
    if missing.is_empty() {
        return Ok(());
    }
    Err(DataFusionError::Execution(format!(
        "plugin '{plugin}': [source.vcf] info_fields names {} the VCF header of '{path}' does \
         not declare: {}. INFO keys are bare and CASE-SENSITIVE — check the spelling and the \
         case. The header declares: [{}].",
        if missing.len() == 1 { "a key" } else { "keys" },
        missing
            .iter()
            .map(|k| format!("'{k}'"))
            .collect::<Vec<_>>()
            .join(", "),
        declared.join(", ")
    )))
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
                // A key the header does not declare would PANIC inside the reader with no
                // context at all (see `reject_undeclared_info_fields`) — diff it first.
                if let Some(requested) = info_fields.as_deref() {
                    reject_undeclared_info_fields(
                        &manifest.plugin_name,
                        &table,
                        &spec.path,
                        requested,
                    )
                    .await?;
                }
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
    use crate::plugin_cache::test_fixtures::{write_bgzf, write_gz};
    use datafusion::arrow::array::{Array, Int64Array};

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
        // TWO INFO keys, of which the manifest selects one. With a single-key header
        // every assertion below would hold identically for `info_fields = None`, and
        // the *selection* half of this test would be vacuous.
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             ##INFO=<ID=NOISE,Number=1,Type=Float,Description=\"an INFO key the manifest omits\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr22\t22893742\t.\tC\tG\t.\t.\tSCORE=0.9;NOISE=1.5\n\
             chr22\t22893800\t.\tA\tT\t.\t.\tSCORE=0.1;NOISE=2.5\n",
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

        // `info_fields` is a SELECTION: only the requested INFO keys become columns.
        // NOISE is declared in the header and present in every record, so if the
        // selection were ignored (or `info_fields` silently dropped, as a typo'd
        // `[source.vcf]` key once caused) it would show up here.
        let table = ctx.table("plugin_demo_src").await.unwrap();
        let columns: Vec<&str> = table
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().as_str())
            .collect();
        assert!(
            columns.contains(&"SCORE"),
            "the selected INFO key must be materialized, got {columns:?}"
        );
        assert!(
            !columns.contains(&"NOISE"),
            "an INFO key absent from info_fields must NOT be materialized, got {columns:?}"
        );

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

    /// A manifest that names an INFO key the VCF header does not declare — a typo, or
    /// (far likelier) the wrong case, since INFO keys are bare and CASE-SENSITIVE — used
    /// to reach `VcfTableProvider::new`, where the reader does
    /// `header_infos.get(tag).unwrap()` and PANICS:
    ///
    /// ```text
    /// thread panicked at datafusion-bio-format-vcf/src/table_provider.rs:174:51:
    /// called `Option::unwrap()` on a `None` value
    /// ```
    ///
    /// No plugin, no manifest, no key, no file. (The careful `.map_err(...)` around
    /// `VcfTableProvider::new` catches nothing: a panic is not an `Err`.) Diff the
    /// requested keys against the header BEFORE handing them to the reader instead.
    #[tokio::test(flavor = "multi_thread")]
    async fn miscased_info_field_is_rejected_with_the_header_keys() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("demo.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             ##INFO=<ID=NOISE,Number=1,Type=Float,Description=\"noise\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr22\t22893742\t.\tC\tG\t.\t.\tSCORE=0.9;NOISE=1.5\n",
        )
        .unwrap();

        // `score` — right key, wrong case. The single most likely authoring error.
        let toml_src = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "{}"
  [source.vcf]
  info_fields = ["score"]

[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml_src).unwrap();
        let ctx = SessionContext::new();
        let err = register_sources(&ctx, &manifest)
            .await
            .expect_err("an INFO key the header does not declare must be an error, not a panic");
        let msg = err.to_string();
        assert!(msg.contains("demo"), "must name the plugin: {msg}");
        assert!(msg.contains("score"), "must name the offending key: {msg}");
        assert!(
            msg.contains("SCORE") && msg.contains("NOISE"),
            "must list the INFO keys the header actually declares: {msg}"
        );
        assert!(
            msg.contains(&vcf.display().to_string()),
            "must name the file: {msg}"
        );
    }

    /// The same guard must work on BGZF (`.vcf.gz`) — the form every real plugin source
    /// ships in — not just plain `.vcf`.
    #[tokio::test(flavor = "multi_thread")]
    async fn unknown_info_field_is_rejected_in_a_bgzf_vcf() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("demo.vcf.gz");
        write_bgzf(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=MMCNT1,Number=1,Type=Integer,Description=\"mm count\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG\t.\t.\tMMCNT1=3\n",
        );

        let toml_src = format!(
            r##"
plugin_name = "mastermind"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "{}"
  [source.vcf]
  info_fields = ["MMCNT1", "BOGUS"]

[[value_columns]]
column = "mmcnt1"
csq_field = "MM"
type = "Int32"
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml_src).unwrap();
        let ctx = SessionContext::new();
        let err = register_sources(&ctx, &manifest)
            .await
            .expect_err("BOGUS is not declared in the BGZF header");
        let msg = err.to_string();
        assert!(msg.contains("BOGUS"), "must name the offending key: {msg}");
        assert!(
            msg.contains("MMCNT1"),
            "must list the header's declared keys: {msg}"
        );
        assert!(msg.contains("mastermind"), "must name the plugin: {msg}");
    }

    /// The valid cases must keep working through the new header check:
    /// - every requested key is declared -> registers;
    /// - `info_fields = []` -> no INFO columns, no header complaint;
    /// - an omitted `[source.vcf]` -> all INFO keys, header never diffed.
    #[tokio::test(flavor = "multi_thread")]
    async fn declared_empty_and_omitted_info_fields_all_register() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("demo.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG\t.\t.\tSCORE=0.9\n",
        )
        .unwrap();

        let mk = |params: &str| {
            format!(
                r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "{}"
{params}

[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
"##,
                vcf.display()
            )
        };

        for (label, params, wants_score) in [
            (
                "exact key",
                "  [source.vcf]\n  info_fields = [\"SCORE\"]",
                true,
            ),
            (
                "empty selection",
                "  [source.vcf]\n  info_fields = []",
                false,
            ),
            ("omitted [source.vcf]", "", true),
        ] {
            let manifest: SourceManifest = toml::from_str(&mk(params)).unwrap();
            let ctx = SessionContext::new();
            let _temps = register_sources(&ctx, &manifest)
                .await
                .unwrap_or_else(|e| panic!("{label} must register: {e}"));
            let table = ctx.table("plugin_demo_src").await.unwrap();
            let has_score = table.schema().field_with_name(None, "SCORE").is_ok();
            assert_eq!(has_score, wants_score, "{label}: SCORE column presence");
        }
    }

    /// Pins the two reader behaviours `VcfParams`' docs warn manifest authors about,
    /// so a bio-formats bump that changes either one fails here instead of silently
    /// invalidating the documented `ingest_sql` guidance:
    ///
    /// 1. a multi-allelic record is ONE row whose `alt` is the ALTs joined by `|`
    ///    (so `concat(ref,'/',alt)` yields the un-matchable `A/G|T`);
    /// 2. `end` is `POS + len(REF) - 1` for an indel, not `start`.
    ///
    /// This is WHY `build::reject_pipe_joined_alleles` exists — and the tail of this
    /// test asserts that guard actually fires on exactly the value pinned here, so the
    /// characterization and the guard can never drift apart.
    #[tokio::test(flavor = "multi_thread")]
    async fn multiallelic_alt_is_pipe_joined_and_indel_end_extends_past_start() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("multi.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG,T\t.\t.\t.\n\
             chr1\t200\t.\tACGT\tA\t.\t.\t.\n",
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

        let batches = ctx
            .sql(
                "SELECT `start`, `end`, `ref`, alt, concat(`ref`, '/', alt) AS allele_string \
                 FROM plugin_demo_src ORDER BY `start`",
            )
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        let b = &batches[0];
        let col = |i: usize| {
            b.column(i)
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringArray>()
                .expect("Utf8 column")
        };
        let starts = b
            .column(0)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::UInt32Array>()
            .expect("start is UInt32");
        let ends = b
            .column(1)
            .as_any()
            .downcast_ref::<datafusion::arrow::array::UInt32Array>()
            .expect("end is UInt32");

        // Row 0: A -> G,T collapses into ONE row, alt = "G|T".
        assert_eq!(
            b.num_rows(),
            2,
            "the reader does NOT split multi-allelic records into one row per ALT"
        );
        assert_eq!(col(3).value(0), "G|T", "ALTs are joined with '|'");
        assert_eq!(
            col(4).value(0),
            "A/G|T",
            "so the naive concat(ref,'/',alt) yields a key that matches nothing"
        );

        // Row 1: ACGT -> A is a deletion; end extends past start (200 + 4 - 1 = 203).
        assert_eq!(starts.value(1), 200);
        assert_eq!(
            ends.value(1),
            203,
            "an indel's end is POS + len(REF) - 1, not start"
        );

        // …and the un-matchable "A/G|T" pinned above must now FAIL the build rather
        // than silently produce a shard that annotates nothing. (Before the guard this
        // returned ChromEntry { rows: 1, warm: 0, cold: 1 } and no error.)
        let var = dir.path().join("var.parquet");
        crate::plugin_cache::test_fixtures::write_variation(&var, &[("1", 100, "A/G", 0)]);
        let naive = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, `start`, `end`, concat(`ref`, '/', alt) AS allele_string,
       CAST(qual AS FLOAT) AS score
FROM plugin_demo_src
"""

[[source]]
provider = "vcf"
path = "{}"

[[value_columns]]
column = "score"
csq_field = "DEMO"
type = "Float32"
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&naive).unwrap();
        let err = crate::plugin_cache::build::build_plugin_chrom(
            &manifest,
            "demo.source.toml",
            &var,
            &dir.path().join("out"),
            "1",
        )
        .await
        .expect_err("the pipe-joined allele_string pinned above must fail the build");
        assert!(
            err.to_string().contains("A/G|T"),
            "the guard must quote the very value this test pins: {err}"
        );
    }
}
