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

use std::io::{BufWriter, Write};
use std::sync::Arc;

use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{CsvReadOptions, ParquetReadOptions, SessionContext};
use datafusion_bio_format_bed::table_provider::{BEDFields, BedTableProvider};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use noodles_core_tabix::Region;
use noodles_csi::BinningIndex;

use crate::plugin_cache::normalize::canonical_contig_str;
use crate::plugin_cache::source_manifest::{
    CoordinateSystem, CsvParams, ProviderKind, SourceIndex, SourceManifest, ValueType,
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

/// Query one chromosome from a BGZF source through its sibling tabix index and
/// materialize only those records to a plain file DataFusion can split-scan.
fn materialize_tabix_chrom(path: &str, chrom: &str) -> Result<(String, tempfile::TempPath)> {
    let mut reader = noodles_tabix::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|e| {
            DataFusionError::Execution(format!(
                "open BGZF/tabix source '{path}' (expected sibling '{path}.tbi'): {e}"
            ))
        })?;

    let header = reader.index().header().ok_or_else(|| {
        DataFusionError::Execution(format!("tabix index for '{path}' has no header"))
    })?;
    let canonical = canonical_contig_str(chrom);
    let mut matching_names = Vec::new();
    for name in header.reference_sequence_names() {
        let raw = std::str::from_utf8(name.as_ref()).map_err(|e| {
            DataFusionError::Execution(format!(
                "tabix index for '{path}' contains a non-UTF-8 contig name: {e}"
            ))
        })?;
        if canonical_contig_str(raw) == canonical {
            matching_names.push(raw.to_owned());
        }
    }
    if matching_names.len() > 1 {
        return Err(DataFusionError::Execution(format!(
            "tabix index for '{path}' has multiple contigs equivalent to '{chrom}', including \
             '{}' and '{}'",
            matching_names[0], matching_names[1]
        )));
    }
    let source_contig = matching_names.pop();

    let mut tmp = tempfile::Builder::new()
        .prefix("plugin_src_")
        .suffix(".tsv")
        .tempfile()
        .map_err(|e| DataFusionError::Execution(format!("create temp for '{path}': {e}")))?;

    // A source can legitimately omit a chromosome. An empty explicit-schema
    // CSV table produces a zero-row shard without scanning any other contig.
    if let Some(source_contig) = source_contig {
        let region = Region::new(source_contig.as_str(), ..);
        let records = reader.query(&region).map_err(|e| {
            DataFusionError::Execution(format!(
                "query BGZF/tabix source '{path}' for contig '{source_contig}': {e}"
            ))
        })?;
        let mut writer = BufWriter::new(tmp.as_file_mut());
        for result in records {
            let record = result.map_err(|e| {
                DataFusionError::Execution(format!(
                    "read BGZF/tabix source '{path}' for contig '{source_contig}': {e}"
                ))
            })?;
            writeln!(writer, "{}", record.as_ref()).map_err(|e| {
                DataFusionError::Execution(format!(
                    "write tabix slice for '{path}' contig '{source_contig}': {e}"
                ))
            })?;
        }
        writer.flush().map_err(|e| {
            DataFusionError::Execution(format!(
                "flush tabix slice for '{path}' contig '{source_contig}': {e}"
            ))
        })?;
    }

    let temp_path = tmp.into_temp_path();
    let plain = temp_path.to_string_lossy().into_owned();
    Ok((plain, temp_path))
}

/// Materialize a plain CSV/TSV path DataFusion can read.
///
/// Tabix sources use indexed chromosome selection above. Non-indexed gzip
/// remains supported for backwards compatibility and is decompressed in full.
/// The workspace cannot enable DataFusion's generic compression feature because
/// its `xz2` dependency collides with noodles-cram's liblzma linkage. Returned
/// temp handles keep the staged file alive only until source parsing finishes.
fn materialize_plain(
    path: &str,
    gzip: bool,
    index: Option<SourceIndex>,
    chrom: Option<&str>,
) -> Result<(String, Option<tempfile::TempPath>)> {
    if index == Some(SourceIndex::Tabix) {
        let chrom = chrom.ok_or_else(|| {
            DataFusionError::Execution(format!(
                "BGZF/tabix source '{path}' requires a chromosome-scoped plugin-cache build"
            ))
        })?;
        let (plain, temp) = materialize_tabix_chrom(path, chrom)?;
        return Ok((plain, Some(temp)));
    }
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

async fn register_sources_impl(
    ctx: &SessionContext,
    manifest: &SourceManifest,
    chrom: Option<&str>,
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
                let (plain_path, temp) = materialize_plain(&spec.path, gz, spec.index, chrom)?;
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
                // `index` is source-level metadata because VCF and delimited
                // text can both be tabix-indexed. VcfTableProvider already
                // owns VCF index discovery/usage, so no extra plumbing is
                // needed here yet; retaining the declaration keeps manifests
                // accurate and gives future planning code an explicit signal.
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
                let provider = if spec.record_layout {
                    provider.with_record_layout().map_err(|e| {
                        DataFusionError::Execution(format!(
                            "carry VCF record layout for source '{table}': {e}"
                        ))
                    })?
                } else {
                    provider
                };
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

/// Register unscoped sources. Indexed text inputs deliberately fail here: use
/// [`register_sources_for_chrom`] so a whole BGZF file is never scanned by
/// accident.
pub async fn register_sources(
    ctx: &SessionContext,
    manifest: &SourceManifest,
) -> Result<Vec<tempfile::TempPath>> {
    register_sources_impl(ctx, manifest, None).await
}

/// Register every source for one chromosome. Tabix-indexed text sources query
/// only that contig; other providers retain their existing behavior. Returned
/// temp handles keep staged records alive for query execution and delete them
/// on drop.
pub async fn register_sources_for_chrom(
    ctx: &SessionContext,
    manifest: &SourceManifest,
    chrom: &str,
) -> Result<Vec<tempfile::TempPath>> {
    register_sources_impl(ctx, manifest, Some(chrom)).await
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::SourceManifest;
    use datafusion::arrow::array::{Array, Int64Array, StringArray};
    use noodles_core_tabix::Position;
    use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
    use std::io::Write;

    fn write_gz(path: &std::path::Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    fn write_bgzf_tabix_bed(path: &std::path::Path, rows: &[(&str, usize, usize, &str)]) {
        let file = std::fs::File::create(path).unwrap();
        let mut writer = noodles_bgzf::io::Writer::new(file);
        let mut indexer = noodles_tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());
        let mut chunk_start = writer.virtual_position();

        for &(chrom, start, end, value) in rows {
            writeln!(writer, "{chrom}\t{start}\t{end}\t{value}").unwrap();
            let chunk_end = writer.virtual_position();
            indexer
                .add_record(
                    chrom,
                    Position::try_from(start + 1).unwrap(),
                    Position::try_from(end).unwrap(),
                    Chunk::new(chunk_start, chunk_end),
                )
                .unwrap();
            chunk_start = chunk_end;
        }
        writer.finish().unwrap();

        let index_path = format!("{}.tbi", path.display());
        let index_file = std::fs::File::create(index_path).unwrap();
        let mut index_writer = noodles_tabix::io::Writer::new(index_file);
        index_writer.write_index(&indexer.build()).unwrap();
    }

    fn indexed_bed_manifest(path: &std::path::Path) -> SourceManifest {
        toml::from_str(&format!(
            r##"
plugin_name = "demo"
coordinate_system = "0-based-half-open"
ingest_sql = "SELECT 1"

[[source]]
provider = "tsv"
path = "{}"
index = "tabix"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "start", type = "Utf8" }},
    {{ name = "end",   type = "Utf8" }},
    {{ name = "value", type = "Utf8" }},
  ]

[[value_columns]]
column = "value"
csq_field = "DEMO"
type = "Utf8"
"##,
            path.display()
        ))
        .unwrap()
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn tabix_source_materializes_only_requested_canonical_contig() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("whole.tsv.gz");
        write_bgzf_tabix_bed(
            &tsv,
            &[
                ("chr1", 99, 100, "one-a"),
                ("chr1", 199, 200, "one-b"),
                ("chr2", 299, 300, "two"),
            ],
        );
        let manifest = indexed_bed_manifest(&tsv);
        let ctx = SessionContext::new();

        // The caller uses canonical "1" while the source index stores "chr1".
        let temps = register_sources_for_chrom(&ctx, &manifest, "1")
            .await
            .unwrap();
        assert_eq!(temps.len(), 1);
        let rows = ctx
            .sql("SELECT chrom, value FROM plugin_demo_src ORDER BY CAST(start AS INT)")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        assert_eq!(rows.iter().map(|b| b.num_rows()).sum::<usize>(), 2);
        let chrom = rows[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let value = rows[0]
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(chrom.value(0), "chr1");
        assert_eq!(chrom.value(1), "chr1");
        assert_eq!(value.value(0), "one-a");
        assert_eq!(value.value(1), "one-b");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn tabix_source_requires_a_sibling_index() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("whole.tsv.gz");
        let file = std::fs::File::create(&tsv).unwrap();
        noodles_bgzf::io::Writer::new(file).finish().unwrap();
        let manifest = indexed_bed_manifest(&tsv);
        let ctx = SessionContext::new();

        let error = register_sources_for_chrom(&ctx, &manifest, "1")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("expected sibling"), "{error}");
        assert!(error.contains(".tbi"), "{error}");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn tabix_subsets_every_part_of_a_multi_source_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let snv = dir.path().join("snv.tsv.gz");
        let indel = dir.path().join("indel.tsv.gz");
        write_bgzf_tabix_bed(&snv, &[("1", 99, 100, "snv-1"), ("2", 199, 200, "snv-2")]);
        write_bgzf_tabix_bed(
            &indel,
            &[("1", 299, 300, "indel-1"), ("2", 399, 400, "indel-2")],
        );
        let source = |part: &str, path: &std::path::Path| {
            format!(
                r##"
[[source]]
part = "{part}"
provider = "csv"
path = "{}"
index = "tabix"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "start", type = "Utf8" }},
    {{ name = "end",   type = "Utf8" }},
    {{ name = "value", type = "Utf8" }},
  ]
"##,
                path.display()
            )
        };
        let manifest: SourceManifest = toml::from_str(&format!(
            r##"
plugin_name = "demo"
coordinate_system = "0-based-half-open"
ingest_sql = """
SELECT * FROM plugin_demo_src_snv
UNION ALL
SELECT * FROM plugin_demo_src_indel
"""
{}
{}
[[value_columns]]
column = "value"
csq_field = "DEMO"
type = "Utf8"
"##,
            source("snv", &snv),
            source("indel", &indel)
        ))
        .unwrap();
        let ctx = SessionContext::new();

        let temps = register_sources_for_chrom(&ctx, &manifest, "1")
            .await
            .unwrap();
        assert_eq!(temps.len(), 2);
        let rows = ctx
            .sql(
                "SELECT value FROM (\
                 SELECT value FROM plugin_demo_src_snv \
                 UNION ALL \
                 SELECT value FROM plugin_demo_src_indel) ORDER BY value",
            )
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        let values = rows[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(values.len(), 2);
        assert_eq!(values.value(0), "indel-1");
        assert_eq!(values.value(1), "snv-1");
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

    #[tokio::test(flavor = "multi_thread")]
    async fn vcf_record_layout_distinguishes_explicit_null_from_absent_key() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("missing-info.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n\
             ##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"Disease name\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t100\t.\tG\tA\t.\t.\tCLNDN=.\n\
             1\t101\t.\tG\tA\t.\t.\t.\n",
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
record_layout = true
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let ctx = SessionContext::new();
        let _temps = register_sources(&ctx, &manifest).await.unwrap();
        let batches = ctx
            .sql(
                "SELECT start, \"CLNDN\", \"_vcf_info_keys\" \
                 FROM plugin_demo_src ORDER BY start",
            )
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        assert_eq!(batch.num_rows(), 2);
        let clndn = batch.column(1);
        assert!(clndn.is_null(0));
        assert!(clndn.is_null(1));
        let info_keys = batch
            .column(2)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(info_keys.value(0), "CLNDN");
        assert_eq!(info_keys.value(1), "");
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

    /// The BED branch hardcodes `BEDFields::BED4` on the claim that
    /// `determine_schema` exposes exactly `chrom, start, end, name` whatever
    /// the variant, so a source needing more than one extra field packs it
    /// into `name` and splits it back out in `ingest_sql`. That claim is load
    /// bearing for every BED manifest and was only ever checked by hand
    /// against one bio-formats release. Pin it the way the VCF ALT shape above
    /// is pinned: a bump that widens or renames the BED schema fails here
    /// rather than silently producing an `ingest_sql` that selects the wrong
    /// column.
    #[tokio::test(flavor = "multi_thread")]
    async fn bed_source_schema_is_pinned_to_chrom_start_end_name() {
        let dir = tempfile::tempdir().unwrap();
        let bed = dir.path().join("packed.bed");
        // BED is 0-based half-open; `name` carries the packed extra fields.
        std::fs::write(&bed, "chr1\t99\t100\trs1|A/G|Pathogenic\n").unwrap();

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "0-based-half-open"
ingest_sql = "SELECT 1"
value_columns = []

[[source]]
provider = "bed"
path = "{}"
"##,
            bed.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let ctx = SessionContext::new();
        let _temps = register_sources(&ctx, &manifest).await.unwrap();

        let schema = ctx
            .table_provider("plugin_demo_src")
            .await
            .unwrap()
            .schema();
        let columns: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert_eq!(
            columns,
            vec!["chrom", "start", "end", "name"],
            "BED4 schema changed -- the `name`-packing trick every BED manifest \
             relies on needs re-auditing, and the hardcoded BEDFields::BED4 in \
             the ProviderKind::Bed branch may no longer be the right choice"
        );

        let rows = ctx
            .sql("SELECT name FROM plugin_demo_src")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        let name = rows[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("name is a plain string; packed extra fields are split in ingest_sql");
        assert_eq!(
            name.value(0),
            "rs1|A/G|Pathogenic",
            "BED `name` must arrive verbatim -- any splitting or trimming here \
             would silently corrupt every packed BED manifest"
        );
    }
}
