use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use termimad::crossterm::style::{Attribute, Color};
use termimad::{Alignment, MadSkin, ROUNDED_TABLE_BORDER_CHARS};

use crate::bench::{BenchmarkReport, PhysicalPlanReport, PositionSidecarReport, ScanResult};
use crate::build::BuildReport;
use crate::config::SandboxConfig;
use crate::inspect::InspectReport;

pub fn write_inspect(config: &SandboxConfig, report: &InspectReport) -> Result<PathBuf> {
    let path = reports_dir(config).join("inspect.json");
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_build(config: &SandboxConfig, report: &BuildReport) -> Result<()> {
    let path = reports_dir(config).join("build.json");
    write_json(&path, report)
}

pub fn write_benchmark(config: &SandboxConfig, report: &BenchmarkReport) -> Result<PathBuf> {
    let path = reports_dir(config).join("benchmark.json");
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_for_tier(
    config: &SandboxConfig,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path =
        reports_dir(config).join(format!("benchmark_{}.json", benchmark_scope_suffix(report)));
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_for_dataset_path(
    dataset_path: &Path,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path = reports_dir_for_dataset_path(dataset_path).join("benchmark.json");
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_for_dataset_path_tier(
    dataset_path: &Path,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path = reports_dir_for_dataset_path(dataset_path)
        .join(format!("benchmark_{}.json", benchmark_scope_suffix(report)));
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_progress(
    config: &SandboxConfig,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path = reports_dir(config).join("benchmark.partial.json");
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_progress_for_tier(
    config: &SandboxConfig,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path = reports_dir(config).join(format!(
        "benchmark_{}.partial.json",
        benchmark_scope_suffix(report)
    ));
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_progress_for_dataset_path(
    dataset_path: &Path,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path = reports_dir_for_dataset_path(dataset_path).join("benchmark.partial.json");
    write_json(&path, report)?;
    Ok(path)
}

pub fn write_benchmark_progress_for_dataset_path_tier(
    dataset_path: &Path,
    report: &BenchmarkReport,
) -> Result<PathBuf> {
    let path = reports_dir_for_dataset_path(dataset_path).join(format!(
        "benchmark_{}.partial.json",
        benchmark_scope_suffix(report)
    ));
    write_json(&path, report)?;
    Ok(path)
}

fn benchmark_scope_suffix(report: &BenchmarkReport) -> String {
    match report.lookup_tier {
        Some(tier) => format!("tier{tier}"),
        None => "all_tiers".to_string(),
    }
}

pub fn write_summary(
    config: &SandboxConfig,
    inspect: &InspectReport,
    benchmark: &BenchmarkReport,
) -> Result<()> {
    let dir = reports_dir(config);
    fs::create_dir_all(&dir).with_context(|| format!("failed to create '{}'", dir.display()))?;
    let effective_config = toml::to_string_pretty(config)?;
    fs::write(dir.join("config.effective.toml"), effective_config)?;
    fs::write(
        dir.join("summary.md"),
        render_summary(config, inspect, benchmark)?,
    )?;
    Ok(())
}

fn write_json<T: serde::Serialize>(path: &PathBuf, value: &T) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create '{}'", parent.display()))?;
    }
    let json = serde_json::to_string_pretty(value)?;
    fs::write(path, json).with_context(|| format!("failed to write '{}'", path.display()))?;
    Ok(())
}

fn reports_dir(config: &SandboxConfig) -> PathBuf {
    config.run_root().join("reports")
}

pub fn reports_dir_for_dataset_path(dataset_path: &Path) -> PathBuf {
    dataset_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("reports")
}

fn render_summary(
    config: &SandboxConfig,
    inspect: &InspectReport,
    benchmark: &BenchmarkReport,
) -> Result<String> {
    let mut out = String::new();
    out.push_str("# Lance Encoding Sandbox Report\n\n");
    out.push_str(&format!("Dataset: `{}`\n\n", inspect.dataset_path));
    out.push_str(&format!(
        "Config: Lance {}, key `{}` as {:?}, projection `everything`.\n\n",
        config.dataset.lance_version, config.key.column, config.key.data_type,
    ));

    out.push_str("## Dataset Size\n\n");
    out.push_str("| Metric | Value |\n|---|---:|\n");
    out.push_str(&metric("Total size", inspect.total_size_bytes));
    out.push_str(&metric("Data size", inspect.data_size_bytes));
    out.push_str(&metric("Index size", inspect.index_size_bytes));
    out.push_str(&metric("Metadata size", inspect.metadata_size_bytes));
    out.push_str(&metric("Other size", inspect.other_size_bytes));
    out.push_str(&format!("| Files | {} |\n", inspect.file_count));
    out.push_str(&format!("| Fragments | {} |\n", inspect.fragments));
    out.push_str(&format!("| Rows | {} |\n", inspect.rows));
    out.push('\n');

    out.push_str("## Benchmark\n\n");
    if let Some(sidecar) = &benchmark.sidecar {
        out.push_str("### Position Sidecar\n\n");
        out.push_str(&render_sidecar_table(sidecar));
    }
    if let Some(prewarm) = &benchmark.cold_prewarm_result {
        out.push_str("### Cold Index Prewarm\n\n");
        out.push_str("| Path | Batch keys | Scans | Rows | Selected positions | Seconds | Resolve s | Take s | Rows/s | Bytes read | IOPs | Requests | Ranges scanned | Fragments scanned |\n");
        out.push_str("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
        out.push_str(&scan_row(prewarm));
        out.push('\n');
    }
    out.push_str("### Measured Scans\n\n");
    out.push_str("| Path | Batch keys | Scans | Rows | Selected positions | Seconds | Resolve s | Take s | Rows/s | Bytes read | IOPs | Requests | Ranges scanned | Fragments scanned |\n");
    out.push_str("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    if let Some(warm_result) = &benchmark.warm_result {
        out.push_str(&scan_row(warm_result));
    }
    for cold in &benchmark.cold_results {
        out.push_str(&scan_row(cold));
    }
    out.push('\n');

    out.push_str(&render_physical_plans(&benchmark.physical_plans));

    out.push_str("## Everything Projection\n\n");
    out.push_str(&format!(
        "- Requested logical fields: {}\n",
        benchmark.requested_everything_fields.len()
    ));
    out.push_str(&format!(
        "- Resolved physical projection fields: {}\n\n",
        benchmark.projected_column_count
    ));
    out.push_str("| Physical field | Logical mapping |\n|---|---|\n");
    for physical in &benchmark.resolved_dataset_projection {
        let children = config
            .structs
            .get(physical)
            .filter(|group| group.enabled)
            .map(|group| group.fields.join(", "))
            .unwrap_or_else(|| physical.clone());
        out.push_str(&format!("| `{}` | `{}` |\n", physical, children));
    }
    out.push('\n');

    out.push_str("## Schema And Encodings\n\n");
    out.push_str("| Field | Type | Nullable | Lance compressed bytes | Input parquet bytes | Input parquet encodings | Input parquet compression | Pages | Metadata encoding | Observed encodings |\n");
    out.push_str("|---|---|---:|---:|---:|---|---|---:|---|---|\n");
    for field in &inspect.logical_fields {
        let requested = requested_metadata(field);
        let observed = if field.encodings_observed.is_empty() {
            String::new()
        } else {
            field
                .encodings_observed
                .iter()
                .map(observed_encoding_label)
                .collect::<Vec<_>>()
                .join("<br>")
        };
        let (parquet_bytes, parquet_encodings, parquet_compression) = input_parquet_columns(field);
        out.push_str(&format!(
            "| `{}` | `{}` | {} | {} | {} | {} | {} | {} | `{}` | {} |\n",
            field.path,
            field.data_type.replace('|', "\\|"),
            field.nullable,
            field.compressed_bytes,
            parquet_bytes,
            parquet_encodings.replace('|', "\\|"),
            parquet_compression.replace('|', "\\|"),
            field.pages,
            requested.replace('|', "\\|"),
            observed.replace('|', "\\|")
        ));
    }
    out.push('\n');

    out.push_str("## Indexes\n\n");
    out.push_str("| Index | Fields | Size | Files |\n|---|---|---:|---:|\n");
    for index in &inspect.indexes {
        out.push_str(&format!(
            "| `{}` | `{}` | {} | {} |\n",
            index.name,
            index.fields.join(", "),
            index
                .size_bytes
                .map(human_bytes)
                .unwrap_or_else(|| "unknown".to_string()),
            index.files.len()
        ));
    }
    Ok(out)
}

pub fn render_inspect_stdout(inspect: &InspectReport) -> String {
    let mut out = String::new();
    out.push_str("# Lance Inspect\n\n");
    out.push_str(&format!("Dataset: `{}`\n\n", inspect.dataset_path));
    out.push_str("## Dataset\n\n");
    out.push_str("|-|-:|\n| **Metric** | **Value** |\n|-|-:|\n");
    out.push_str(&format!(
        "| Configured Lance version | {} |\n",
        inspect.configured_lance_version
    ));
    out.push_str(&format!(
        "| Observed Lance file versions | {} |\n",
        inspect.observed_lance_file_versions.join(", ")
    ));
    out.push_str(&format!("| Rows | {} |\n", inspect.rows));
    out.push_str(&format!("| Fragments | {} |\n", inspect.fragments));
    out.push_str(&format!("| Files | {} |\n", inspect.file_count));
    out.push_str(&metric("Total size", inspect.total_size_bytes));
    out.push_str(&metric("Data size", inspect.data_size_bytes));
    out.push_str(&metric("Index size", inspect.index_size_bytes));
    out.push_str(&metric("Metadata size", inspect.metadata_size_bytes));
    out.push_str(&metric("Other size", inspect.other_size_bytes));
    out.push('\n');

    out.push_str("## Everything Projection\n\n");
    out.push_str(&format!(
        "- Requested logical fields: {}\n",
        inspect.requested_everything_fields.len()
    ));
    out.push_str(&format!(
        "- Resolved physical projection fields: {}\n\n",
        inspect.resolved_dataset_projection.len()
    ));

    out.push_str("## Schema And Encodings\n\n");
    out.push_str("|-|-|:-:|-:|-:|-|-|-:|-|-|\n");
    out.push_str("| **Field** | **Type** | **Nullable** | **Lance compressed bytes** | **Input parquet bytes** | **Input parquet encodings** | **Input parquet compression** | **Pages** | **Metadata encoding** | **Observed encodings** |\n");
    out.push_str("|-|-|:-:|-:|-:|-|-|-:|-|-|\n");
    for field in &inspect.logical_fields {
        let observed = if field.encodings_observed.is_empty() {
            String::new()
        } else {
            field
                .encodings_observed
                .iter()
                .map(observed_encoding_label)
                .collect::<Vec<_>>()
                .join("<br>")
        };
        let (parquet_bytes, parquet_encodings, parquet_compression) = input_parquet_columns(field);
        out.push_str(&format!(
            "| `{}` | `{}` | {} | {} | {} | {} | {} | {} | `{}` | {} |\n",
            field.path,
            field.data_type.replace('|', "\\|"),
            field.nullable,
            field.compressed_bytes,
            parquet_bytes,
            parquet_encodings.replace('|', "\\|"),
            parquet_compression.replace('|', "\\|"),
            field.pages,
            requested_metadata(field).replace('|', "\\|"),
            observed.replace('|', "\\|")
        ));
    }
    out.push('\n');

    out.push_str("## Indexes\n\n");
    out.push_str("|-|-|-:|-:|\n| **Index** | **Fields** | **Size** | **Files** |\n|-|-|-:|-:|\n");
    for index in &inspect.indexes {
        out.push_str(&format!(
            "| `{}` | `{}` | {} | {} |\n",
            index.name,
            index.fields.join(", "),
            index
                .size_bytes
                .map(human_bytes)
                .unwrap_or_else(|| "unknown".to_string()),
            index.files.len()
        ));
    }
    out
}

pub fn render_benchmark_stdout(benchmark: &BenchmarkReport) -> String {
    let mut out = String::new();
    out.push_str("# Lance Benchmark\n\n");
    out.push_str(&format!("Dataset: `{}`\n\n", benchmark.dataset_path));
    out.push_str(&format!(
        "- Positions requested: {}\n",
        benchmark.positions_requested
    ));
    out.push_str(&format!("- Lookup scope: {}\n", benchmark.lookup_scope));
    out.push_str(&format!(
        "- Requested logical fields: {}\n",
        benchmark.requested_everything_fields.len()
    ));
    out.push_str(&format!(
        "- Resolved physical projection fields: {}\n\n",
        benchmark.projected_column_count
    ));
    if let Some(sidecar) = &benchmark.sidecar {
        out.push_str("## Position Sidecar\n\n");
        out.push_str(&render_sidecar_table(sidecar));
    }
    if let Some(prewarm) = &benchmark.cold_prewarm_result {
        out.push_str("## Cold Index Prewarm\n\n");
        out.push_str("|-|:-:|:-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|\n");
        out.push_str("| **Path** | **Batch keys** | **Scans** | **Rows** | **Selected positions** | **Seconds** | **Resolve s** | **Take s** | **Rows/s** | **Bytes read** | **IOPs** | **Requests** | **Ranges scanned** | **Fragments scanned** |\n");
        out.push_str("|-|:-:|:-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|\n");
        out.push_str(&scan_row(prewarm));
        out.push('\n');
    }

    out.push_str("## Measured Scans\n\n");
    out.push_str("|-|:-:|:-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|\n");
    out.push_str("| **Path** | **Batch keys** | **Scans** | **Rows** | **Selected positions** | **Seconds** | **Resolve s** | **Take s** | **Rows/s** | **Bytes read** | **IOPs** | **Requests** | **Ranges scanned** | **Fragments scanned** |\n");
    out.push_str("|-|:-:|:-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|-:|\n");
    if let Some(warm_result) = &benchmark.warm_result {
        out.push_str(&scan_row(warm_result));
    }
    for cold in &benchmark.cold_results {
        out.push_str(&scan_row(cold));
    }
    out.push('\n');
    out.push_str(&render_physical_plans(&benchmark.physical_plans));
    out
}

pub fn print_terminal_markdown(markdown: &str) {
    let (width, _) = termimad::terminal_size();
    print!("{}", render_terminal_markdown(markdown, width as usize));
}

pub fn render_terminal_markdown(markdown: &str, width: usize) -> String {
    terminal_skin()
        .text(markdown, Some(width.max(40)))
        .to_string()
}

fn terminal_skin() -> MadSkin {
    let mut skin = MadSkin::default();
    skin.table_border_chars = ROUNDED_TABLE_BORDER_CHARS;
    skin.set_headers_fg(Color::Cyan);
    skin.headers[0].align = Alignment::Left;
    skin.headers[0].add_attr(Attribute::Bold);
    skin.headers[1].set_fg(Color::Blue);
    skin.headers[1].add_attr(Attribute::Bold);
    skin.bold.set_fg(Color::Yellow);
    skin.inline_code.set_fg(Color::Green);
    skin.table.align = Alignment::Left;
    skin
}

fn metric(name: &str, bytes: u64) -> String {
    format!("| {name} | {} |\n", human_bytes(bytes))
}

fn scan_row(result: &ScanResult) -> String {
    let rows_per_second = if result.seconds > 0.0 {
        result.rows as f64 / result.seconds
    } else {
        0.0
    };
    format!(
        "| `{}` | {} | {} | {} | {} | {:.6} | {} | {} | {:.0} | {} | {} | {} | {} | {} |\n",
        result.name,
        result.lookup_batch_size.unwrap_or_default(),
        result.scans,
        result.rows,
        result.selected_positions,
        result.seconds,
        optional_seconds(result.row_id_resolve_seconds),
        optional_seconds(result.take_rows_seconds),
        rows_per_second,
        human_bytes(result.io.bytes_read as u64),
        result.io.iops,
        result.io.requests,
        result.io.ranges_scanned,
        result.io.fragments_scanned,
    )
}

fn render_sidecar_table(sidecar: &PositionSidecarReport) -> String {
    format!(
        "Sidecar path: `{}`\n\n\
         | Built | Build s | Load s | Sidecar read s | BTree page load s | BTree stream batches | BTree stream rows | Pair compare s | Unique positions | Row ids |\n\
         |---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n\
         | {} | {:.6} | {:.6} | {} | {} | {} | {} | {} | {} | {} |\n\n",
        sidecar.path,
        sidecar.built,
        sidecar.build_seconds,
        sidecar.load_seconds,
        optional_seconds(sidecar.sidecar_read_seconds),
        optional_seconds(sidecar.btree_page_load_seconds),
        optional_usize(sidecar.btree_page_stream_batches),
        optional_usize(sidecar.btree_page_stream_rows),
        optional_seconds(sidecar.pair_compare_seconds),
        sidecar.unique_positions,
        sidecar.row_ids
    )
}

fn optional_seconds(seconds: Option<f64>) -> String {
    seconds
        .map(|value| format!("{value:.6}"))
        .unwrap_or_else(|| "-".to_string())
}

fn optional_usize(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "-".to_string())
}

fn render_physical_plans(plans: &[PhysicalPlanReport]) -> String {
    const MAX_INLINE_CHARS: usize = 220;
    const MAX_PLAN_LINE_CHARS: usize = 240;

    if plans.is_empty() {
        return String::new();
    }

    let mut out = String::new();
    out.push_str("## Physical Plans\n\n");
    for plan in plans {
        out.push_str(&format!("### {}\n\n", plan.name));
        out.push_str(&format!(
            "- Batch keys: {}\n",
            plan.lookup_batch_size.unwrap_or_default()
        ));
        out.push_str(&format!(
            "- Projected fields: {}\n",
            plan.projected_column_count
        ));
        let filter = truncate_display_line(&plan.filter.replace('`', "\\`"), MAX_INLINE_CHARS);
        out.push_str(&format!("- Filter: `{}`\n\n", filter));
        out.push_str("```text\n");
        let display_plan = render_plan_text_for_terminal(&plan.plan, MAX_PLAN_LINE_CHARS);
        out.push_str(&display_plan);
        if !display_plan.ends_with('\n') {
            out.push('\n');
        }
        out.push_str("```\n\n");
    }
    out
}

fn render_plan_text_for_terminal(plan: &str, max_line_chars: usize) -> String {
    plan.lines()
        .map(|line| truncate_display_line(line, max_line_chars))
        .collect::<Vec<_>>()
        .join("\n")
}

fn truncate_display_line(line: &str, max_chars: usize) -> String {
    let len = line.chars().count();
    if len <= max_chars {
        return line.to_string();
    }

    let keep = max_chars.saturating_sub(40).max(16);
    let omitted = len.saturating_sub(keep);
    let marker = format!(" ... <truncated {omitted} chars>");
    let keep = max_chars.saturating_sub(marker.chars().count()).max(16);
    let omitted = len.saturating_sub(keep);
    let marker = format!(" ... <truncated {omitted} chars>");

    let mut truncated = line.chars().take(keep).collect::<String>();
    truncated.push_str(&marker);
    truncated
}

fn requested_metadata(field: &crate::inspect::FieldInspectRow) -> String {
    let mut parts = Vec::new();
    for key in [
        "lance-encoding:structural-encoding",
        "lance-encoding:compression",
        "lance-encoding:compression-level",
        "lance-encoding:dict-values-compression",
        "lance-encoding:dict-values-compression-level",
        "lance-encoding:minichunk-size",
        "lance-encoding:packed",
    ] {
        if let Some(value) = field.metadata.get(key) {
            parts.push(format!("{key}={value}"));
        }
    }
    parts.join("; ")
}

fn observed_encoding_label(encoding: &crate::inspect::ObservedEncoding) -> String {
    format!(
        "layout={}; encoding={}; compression={}",
        encoding.layout, encoding.encoding, encoding.compression
    )
}

fn input_parquet_columns(field: &crate::inspect::FieldInspectRow) -> (String, String, String) {
    field
        .input_parquet
        .as_ref()
        .map(|input| {
            (
                input.compressed_bytes.to_string(),
                input.encodings.join(","),
                input.compression.join(","),
            )
        })
        .unwrap_or_else(|| (String::new(), String::new(), String::new()))
}

fn human_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0usize;
    while value >= 1024.0 && unit + 1 < UNITS.len() {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} {}", UNITS[unit])
    } else {
        format!("{value:.2} {}", UNITS[unit])
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::bench::{
        BenchmarkReport, PhysicalPlanReport, PositionSidecarReport, ScanIoStats, ScanResult,
    };
    use crate::inspect::{FieldInspectRow, InspectReport, ObservedEncoding, ParquetInputStats};

    #[test]
    fn inspect_stdout_includes_dataset_schema_and_indexes() {
        let report = InspectReport {
            dataset_path: "/tmp/chr1.lance".to_string(),
            configured_lance_version: "2.1".to_string(),
            observed_lance_file_versions: vec!["2.1".to_string()],
            fragments: 2,
            rows: 42,
            total_size_bytes: 1024,
            data_size_bytes: 900,
            index_size_bytes: 100,
            metadata_size_bytes: 24,
            other_size_bytes: 0,
            file_count: 4,
            requested_everything_fields: vec!["start".to_string()],
            resolved_dataset_projection: vec!["start".to_string()],
            logical_fields: vec![FieldInspectRow {
                path: "start".to_string(),
                id: 1,
                data_type: "Int64".to_string(),
                nullable: false,
                metadata: BTreeMap::from([(
                    "lance-encoding:compression".to_string(),
                    "zstd".to_string(),
                )]),
                compressed_bytes: 128,
                pages: 1,
                input_parquet: Some(ParquetInputStats {
                    compressed_bytes: 256,
                    files: 2,
                    column_chunks: 4,
                    encodings: vec!["PLAIN".to_string(), "RLE_DICTIONARY".to_string()],
                    compression: vec!["ZSTD".to_string()],
                }),
                encodings_observed: vec![ObservedEncoding {
                    layout: "MiniBlock".to_string(),
                    encoding: "Values+Flat(64-bit)".to_string(),
                    compression: "zstd(level=3)".to_string(),
                }],
            }],
            physical_columns: Vec::new(),
            indexes: Vec::new(),
        };

        let rendered = super::render_inspect_stdout(&report);
        assert!(rendered.contains("# Lance Inspect"));
        assert!(rendered.contains("Dataset: `/tmp/chr1.lance`"));
        assert!(rendered.contains(
            "| `start` | `Int64` | false | 128 | 256 | PLAIN,RLE_DICTIONARY | ZSTD | 1 |"
        ));
        assert!(rendered.contains("PLAIN,RLE_DICTIONARY"));
        assert!(rendered.contains("ZSTD"));
        assert!(rendered.contains("lance-encoding:compression=zstd"));
        assert!(
            rendered.contains(
                "layout=MiniBlock; encoding=Values+Flat(64-bit); compression=zstd(level=3)"
            )
        );
        assert!(rendered.contains("## Indexes"));
    }

    #[test]
    fn benchmark_stdout_renders_markdown_table() {
        let report = BenchmarkReport {
            dataset_path: "/tmp/chr1.lance".to_string(),
            lookup_tier: Some(1),
            lookup_scope: "tier=1".to_string(),
            requested_everything_fields: vec!["start".to_string(), "allele_string".to_string()],
            resolved_dataset_projection: vec!["start".to_string(), "allele_string".to_string()],
            projected_column_count: 2,
            positions_requested: 10_000,
            warm_result: Some(scan_result(
                "warm_full_scan_all_everything_columns",
                None,
                1,
                100,
                0.5,
            )),
            cold_prewarm_result: Some(scan_result(
                "cold_index_prewarm_key_only",
                Some(238),
                1,
                2,
                0.25,
            )),
            sidecar: Some(PositionSidecarReport {
                path: "/tmp/position_row_ids.bin".to_string(),
                built: false,
                build_seconds: 0.0,
                load_seconds: 0.125,
                sidecar_read_seconds: Some(0.025),
                btree_page_load_seconds: Some(0.09),
                btree_metadata_seconds: Some(0.005),
                tier_filter_seconds: Some(0.015),
                btree_page_scan_seconds: Some(0.07),
                btree_page_stream_batches: Some(20),
                btree_page_stream_rows: Some(100),
                pair_compare_seconds: Some(0.01),
                unique_positions: 80,
                row_ids: 100,
            }),
            cold_results: vec![scan_result("cold_lookup_batch_512", Some(512), 20, 12, 4.0)],
            physical_plans: vec![PhysicalPlanReport {
                name: "warm_full_scan_all_everything_columns".to_string(),
                filter: "tier = 0".to_string(),
                lookup_batch_size: None,
                projected_column_count: 2,
                projection: vec!["start".to_string(), "allele_string".to_string()],
                plan: "LanceScanExec: uri=/tmp/chr1.lance".to_string(),
            }],
        };

        let rendered = super::render_benchmark_stdout(&report);

        assert!(rendered.contains("# Lance Benchmark"));
        assert!(rendered.contains("- Positions requested: 10000"));
        assert!(rendered.contains("- Lookup scope: tier=1"));
        assert!(rendered.contains("## Position Sidecar"));
        assert!(rendered.contains("Sidecar path: `/tmp/position_row_ids.bin`"));
        assert!(rendered.contains(
            "| false | 0.000000 | 0.125000 | 0.025000 | 0.090000 | 20 | 100 | 0.010000 | 80 | 100 |"
        ));
        assert!(rendered.contains("## Cold Index Prewarm"));
        assert!(rendered.contains("| `cold_index_prewarm_key_only` | 238 | 1 | 2 | 2 | 0.250000 | - | - | 8 | 1.00 KiB | 7 | 3 | 5 | 2 |"));
        assert!(rendered.contains("## Measured Scans"));
        assert!(rendered.contains("| **Path** | **Batch keys** | **Scans** | **Rows** | **Selected positions** | **Seconds** | **Resolve s** | **Take s** | **Rows/s** | **Bytes read** | **IOPs** | **Requests** | **Ranges scanned** | **Fragments scanned** |"));
        assert!(rendered.contains("| `warm_full_scan_all_everything_columns` | 0 | 1 | 100 | 100 | 0.500000 | - | - | 200 | 1.00 KiB | 7 | 3 | 5 | 2 |"));
        assert!(rendered.contains(
            "| `cold_lookup_batch_512` | 512 | 20 | 12 | 12 | 4.000000 | - | - | 3 | 1.00 KiB | 7 | 3 | 5 | 2 |"
        ));
        assert!(rendered.contains("## Physical Plans"));
        assert!(rendered.contains("### warm_full_scan_all_everything_columns"));
        assert!(rendered.contains("```text\nLanceScanExec: uri=/tmp/chr1.lance\n```"));
        assert!(super::render_terminal_markdown(&rendered, 120).contains('╭'));
    }

    #[test]
    fn terminal_markdown_renderer_formats_tables() {
        let rendered = super::render_terminal_markdown(
            "# Title\n\n|-|-:|\n| **Metric** | **Value** |\n|-|-:|\n| Rows | 42 |\n",
            80,
        );

        assert!(rendered.contains("Title"));
        assert!(rendered.contains("Rows"));
        assert!(rendered.contains("42"));
        assert!(rendered.contains('╭'));
    }

    #[test]
    fn benchmark_stdout_bounds_physical_plan_lines_for_terminal_rendering() {
        let mut report = BenchmarkReport {
            dataset_path: "/tmp/chr1.lance".to_string(),
            lookup_tier: Some(1),
            lookup_scope: "tier=1".to_string(),
            requested_everything_fields: vec!["start".to_string()],
            resolved_dataset_projection: vec!["start".to_string()],
            projected_column_count: 1,
            positions_requested: 10_000,
            warm_result: Some(scan_result(
                "warm_full_scan_all_everything_columns",
                None,
                1,
                100,
                0.5,
            )),
            cold_prewarm_result: None,
            sidecar: None,
            cold_results: Vec::new(),
            physical_plans: Vec::new(),
        };
        let long_line = format!("ProjectionExec: expr=[{}]", "x".repeat(5_000));
        report.physical_plans.push(PhysicalPlanReport {
            name: "cold_lookup_batch_2048".to_string(),
            filter: format!("tier = 1 AND position IN ({})", "1, ".repeat(2_048)),
            lookup_batch_size: Some(2_048),
            projected_column_count: 47,
            projection: vec!["position".to_string()],
            plan: long_line,
        });

        let rendered = super::render_benchmark_stdout(&report);

        assert!(rendered.contains("<truncated"));
        assert!(!rendered.contains(&"x".repeat(1_000)));
        assert!(rendered.lines().all(|line| line.len() < 512));
        assert!(super::render_terminal_markdown(&rendered, 120).contains("<truncated"));
    }

    #[test]
    fn btree_resolve_seconds_is_serialized_only_for_btree_results() {
        let mut btree = scan_result("cold_btree_direct_take_everything", None, 1, 10, 0.25);
        btree.row_id_resolve_seconds = Some(0.05);
        btree.btree_resolve_seconds = Some(0.05);
        let btree_json = serde_json::to_value(&btree).expect("scan result serializes");
        assert_eq!(btree_json["row_id_resolve_seconds"], 0.05);
        assert_eq!(btree_json["btree_resolve_seconds"], 0.05);

        let sidecar = scan_result("cold_sidecar_take_everything", None, 1, 10, 0.25);
        let sidecar_json = serde_json::to_value(&sidecar).expect("scan result serializes");
        assert!(sidecar_json.get("btree_resolve_seconds").is_none());
    }

    fn scan_result(
        name: &str,
        lookup_batch_size: Option<usize>,
        scans: usize,
        rows: usize,
        seconds: f64,
    ) -> ScanResult {
        ScanResult {
            name: name.to_string(),
            filter_keys: lookup_batch_size.unwrap_or_default(),
            lookup_batch_size,
            scans,
            rows,
            result_batches: 1,
            selected_positions: rows,
            seconds,
            row_id_resolve_seconds: None,
            btree_resolve_seconds: None,
            take_rows_seconds: None,
            io: ScanIoStats {
                iops: 7,
                requests: 3,
                bytes_read: 1024,
                indices_loaded: 0,
                parts_loaded: 0,
                index_comparisons: 0,
                fragments_scanned: 2,
                ranges_scanned: 5,
                rows_scanned: rows,
            },
        }
    }
}
