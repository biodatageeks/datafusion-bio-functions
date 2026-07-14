//! Shared harness for Perl→Rust ported golden tests (v115).
//!
//! Per-port tests call one of three entry points:
//! - [`run_and_compare_csq`]                       — single golden per case (C1)
//! - [`run_and_compare_csq_with_flags`]            — per-flag-set golden    (C2)
//! - [`run_and_compare_csq_filtered`]              — row-predicate filter   (C3)
//!
//! Each helper annotates `input.vcf` against the shared v115 parquet cache
//! (`vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep/`), reads
//! the CSQ field list from the golden's `##INFO=<ID=CSQ,...,Format: …>` header,
//! and compares per-row, per-allele-group, per-field. Mismatches on
//! `hard_fields` fail the test; mismatches on the rest are logged but
//! non-fatal (per verification.md Layer 5).
//!
//! Reference pattern: `tests/vcf_roundtrip_golden.rs` (the workspace's
//! source-of-truth for the annotate→read-back→compare flow).
//!
//! Spec sources:
//!   - porting-tests/2026-05-25-implementation-design.md §2 Batch 0
//!   - porting-tests/REVIEW_NOTES.md §C1 / §C2 / §C3
//!   - porting-tests/verification.md Layer 5 (hard/soft separation)
#![allow(dead_code)]

use std::collections::HashSet;
use std::sync::Arc;

use datafusion::arrow::array::{Array, LargeListArray, ListArray, StringArray, StringViewArray};
use datafusion::common::Result;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::vcf_sink;

/// Predicate applied to a golden CSQ group's tokens to decide whether the group
/// takes part in the comparison (REVIEW_NOTES.md §C3). The lifetime parameter
/// keeps non-`'static` closures usable, which the `dyn` object-lifetime default
/// would otherwise forbid.
type RowFilter<'a> = dyn Fn(&[String]) -> bool + Sync + 'a;

// ───────────────────────── path helpers ─────────────────────────

/// Resolve a path relative to the workspace root.
fn workspace_path(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

/// Detect a Git LFS pointer file (so a partially-pulled checkout SKIPs
/// rather than panicking on what looks like a fixture).
fn is_lfs_pointer(path: &std::path::Path) -> bool {
    std::fs::read_to_string(path)
        .map(|s| s.starts_with("version https://git-lfs.github.com"))
        .unwrap_or(false)
}

// ───────────────────────── CSQ extraction ─────────────────────────

/// Join the string elements of one list cell (a per-row `List<Utf8>`),
/// reconstructing the comma-separated CSQ string. VEP encodes CSQ as a
/// multi-valued (`Number=.`) INFO field, so `VcfTableProvider` surfaces
/// it as a List and each element is one per-allele-group CSQ entry.
fn join_string_elements(elems: &dyn Array) -> String {
    if let Some(s) = elems.as_any().downcast_ref::<StringArray>() {
        return (0..s.len())
            .filter(|&i| !s.is_null(i))
            .map(|i| s.value(i))
            .collect::<Vec<_>>()
            .join(",");
    }
    if let Some(s) = elems.as_any().downcast_ref::<StringViewArray>() {
        return (0..s.len())
            .filter(|&i| !s.is_null(i))
            .map(|i| s.value(i))
            .collect::<Vec<_>>()
            .join(",");
    }
    panic!(
        "port_common: unhandled CSQ list-element type {:?} — refusing to compare vacuously",
        elems.data_type()
    );
}

/// Extract the `CSQ` value for one row as a String. Handles scalar string
/// representations AND the `List<Utf8>` representation the VCF provider
/// emits for multi-valued CSQ. PANICS on any unhandled array type rather
/// than silently returning `""` — a silent empty string here would make
/// the whole comparison vacuous (empty == empty → always 0 mismatches).
fn csq_at(col: &dyn Array, row: usize) -> String {
    if col.is_null(row) {
        return String::new();
    }
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return a.value(row).to_string();
    }
    if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
        return a.value(row).to_string();
    }
    if let Some(a) = col.as_any().downcast_ref::<ListArray>() {
        return join_string_elements(a.value(row).as_ref());
    }
    if let Some(a) = col.as_any().downcast_ref::<LargeListArray>() {
        return join_string_elements(a.value(row).as_ref());
    }
    panic!(
        "port_common: unhandled CSQ array type {:?} — refusing to compare vacuously",
        col.data_type()
    );
}

/// Read the `CSQ` column out of a record batch into per-row strings.
fn read_csq_column(batch: &datafusion::arrow::record_batch::RecordBatch) -> Vec<String> {
    let col = batch.column(
        batch
            .schema()
            .index_of("CSQ")
            .expect("CSQ column not found in batch"),
    );
    (0..batch.num_rows())
        .map(|row| csq_at(col.as_ref(), row))
        .collect()
}

/// Parse one row's CSQ string into per-allele-group token lists. Each
/// group is one transcript/feature annotation (pipe-delimited tokens).
fn parse_csq_row(csq: &str) -> Vec<Vec<String>> {
    if csq.is_empty() {
        return Vec::new();
    }
    csq.split(',')
        .map(|group| group.split('|').map(str::to_string).collect())
        .collect()
}

/// Extract the CSQ `Format: ` field-name list from a golden VCF header
/// line such as `##INFO=<ID=CSQ,...,Description="...Format: A|B|C|...">`.
fn extract_csq_format(header_text: &str) -> Vec<String> {
    let Some(line) = header_text
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
    else {
        return Vec::new();
    };
    let Some(start) = line.find("Format: ") else {
        return Vec::new();
    };
    let rest = &line[start + "Format: ".len()..];
    let end = rest.find('"').unwrap_or(rest.len());
    rest[..end].split('|').map(str::to_string).collect()
}

// ───────────────────────── config helpers ─────────────────────────

/// The standard `--everything` config used by every port. Returns a
/// fresh `AnnotateVcfConfig` so callers can mutate it (e.g. apply
/// per-test CLI flags for C2).
fn base_config(ref_fasta: &str) -> vcf_sink::AnnotateVcfConfig {
    vcf_sink::AnnotateVcfConfig {
        everything: true,
        extended_probes: true,
        reference_fasta_path: Some(ref_fasta.to_string()),
        ..Default::default()
    }
}

/// Apply one CLI flag override to an existing config. Supported set is
/// intentionally narrow — extend per-port-need rather than speculatively.
fn apply_cli_flag(config: &mut vcf_sink::AnnotateVcfConfig, key: &str, value: &str) {
    match (key, value) {
        ("pick", "1") => config.flag_pick = true,
        ("pick_allele", "1") => config.flag_pick_allele = true,
        ("pick_allele_gene", "1") => config.flag_pick_allele_gene = true,
        ("per_gene", "1") => config.per_gene = true,
        ("gencode_basic", "1") => config.gencode_basic = true,
        ("gencode_primary", "1") => config.gencode_primary = true,
        ("exclude_predicted", "1") => config.exclude_predicted = true,
        ("distance", v) => config.distance = Some(v.to_string()),
        (k, v) => eprintln!("port_common: unsupported cli_flag ({k:?},{v:?}) — ignored"),
    }
}

// ───────────────────────── shared body ─────────────────────────

/// All three entry points delegate to this. `row_filter` (when `Some`)
/// is invoked on each golden CSQ allele-group's parsed token list; the
/// group is included in comparison iff the filter returns true.
async fn compare_csq_inner(
    input_vcf: &str,
    golden_vcf: &str,
    cache_path: &str,
    config: &vcf_sink::AnnotateVcfConfig,
    hard_fields: &[&str],
    row_filter: Option<&RowFilter<'_>>,
    case_label: &str,
) -> Result<bool> {
    // 1. Annotate input.vcf to a tmpfile.
    // NOTE: when PORT_DEBUG_OUTPUT_DIR is set, write to that dir instead of a
    // tempdir so the failed-test output VCF persists for inspection.
    let debug_dir = std::env::var("PORT_DEBUG_OUTPUT_DIR").ok();
    let _tmp_dir_keepalive;
    let output_path = if let Some(d) = debug_dir.as_deref() {
        std::fs::create_dir_all(d).ok();
        std::path::PathBuf::from(d).join(format!("{case_label}.annotated.vcf"))
    } else {
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("annotated.vcf");
        _tmp_dir_keepalive = td;
        p
    };
    let rows_written = vcf_sink::annotate_to_vcf(
        input_vcf,
        cache_path,
        "parquet",
        output_path.to_str().unwrap(),
        config,
    )
    .await?;
    eprintln!(
        "port_common[{case_label}]: annotate_to_vcf returned {rows_written} rows; output at {}",
        output_path.display()
    );

    // 2. Read engine output + golden via VcfTableProvider. The provider
    //    uses futures::executor::block_on() internally; wrap in
    //    spawn_blocking to avoid the deadlock inside the test runtime.
    let ctx = SessionContext::new();
    let out_str = output_path.display().to_string();
    let output_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(out_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    let golden_str = golden_vcf.to_string();
    let golden_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(golden_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    ctx.register_table("output_vcf", Arc::new(output_prov))?;
    ctx.register_table("golden_vcf", Arc::new(golden_prov))?;

    let output_batches = ctx
        .sql("SELECT * FROM output_vcf ORDER BY start")
        .await?
        .collect()
        .await?;
    let golden_batches = ctx
        .sql("SELECT * FROM golden_vcf ORDER BY start")
        .await?
        .collect()
        .await?;
    let output_batch =
        datafusion::arrow::compute::concat_batches(&output_batches[0].schema(), &output_batches)
            .unwrap();
    let golden_batch =
        datafusion::arrow::compute::concat_batches(&golden_batches[0].schema(), &golden_batches)
            .unwrap();

    let ours = read_csq_column(&output_batch);
    let golden = read_csq_column(&golden_batch);

    // 3. Vacuity guard + row-count match.
    assert!(
        golden.iter().any(|s| !s.is_empty()),
        "port case '{case_label}': golden CSQ extraction produced only empty strings — \
         comparator is reading nothing (vacuity guard tripped)"
    );
    assert_eq!(
        ours.len(),
        golden.len(),
        "port case '{case_label}': output has {} rows, golden has {} rows",
        ours.len(),
        golden.len()
    );

    // 4. Parse Format field list from the golden header.
    let golden_header_text = std::fs::read_to_string(golden_vcf).unwrap_or_default();
    let csq_format = extract_csq_format(&golden_header_text);
    assert!(
        !csq_format.is_empty(),
        "port case '{case_label}': could not extract CSQ Format declaration \
         (looked for '##INFO=<ID=CSQ,...Format: ...')"
    );

    // 5. Field-by-field compare, separating hard from soft.
    let hard_set: HashSet<&str> = hard_fields.iter().copied().collect();
    let mut hard_mismatches: Vec<String> = Vec::new();
    let mut soft_mismatches: Vec<String> = Vec::new();

    for (row_idx, (ours_csq, golden_csq)) in ours.iter().zip(golden.iter()).enumerate() {
        let ours_groups = parse_csq_row(ours_csq);
        let golden_groups = parse_csq_row(golden_csq);

        // Compare group-by-group. The number of groups (one per
        // transcript/feature/allele) must match in steady state; if they
        // don't, that's itself a structural failure recorded as hard
        // under a synthetic "_group_count" sentinel.
        if ours_groups.len() != golden_groups.len() {
            hard_mismatches.push(format!(
                "row {row_idx} group-count: ours={} groups, golden={} groups",
                ours_groups.len(),
                golden_groups.len()
            ));
            continue;
        }

        for (g_idx, (og, gg)) in ours_groups.iter().zip(golden_groups.iter()).enumerate() {
            // Apply row filter (C3) on the GOLDEN group's tokens.
            if let Some(filter) = row_filter
                && !filter(gg)
            {
                continue;
            }
            for (fld_idx, fld_name) in csq_format.iter().enumerate() {
                let ov = og.get(fld_idx).map(String::as_str).unwrap_or("");
                let gv = gg.get(fld_idx).map(String::as_str).unwrap_or("");
                if ov != gv {
                    let msg = format!(
                        "row {row_idx} group {g_idx} field '{fld_name}': ours={ov:?} golden={gv:?}"
                    );
                    if hard_set.contains(fld_name.as_str()) {
                        hard_mismatches.push(msg);
                    } else {
                        soft_mismatches.push(msg);
                    }
                }
            }
        }
    }

    if !soft_mismatches.is_empty() {
        eprintln!(
            "port case '{case_label}': {} soft-field mismatch(es) (non-fatal):",
            soft_mismatches.len()
        );
        for m in soft_mismatches.iter().take(20) {
            eprintln!("  {m}");
        }
        if soft_mismatches.len() > 20 {
            eprintln!("  … and {} more", soft_mismatches.len() - 20);
        }
    }

    assert!(
        hard_mismatches.is_empty(),
        "port case '{case_label}': {} hard-field mismatch(es):\n{}",
        hard_mismatches.len(),
        hard_mismatches
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );

    Ok(true)
}

// ───────────────────────── public API ─────────────────────────

/// (C1) Annotate `tests/data/port/<name>/input.vcf` against the shared
/// v115 parquet cache and compare the produced CSQ column against
/// `tests/data/port/<name>/golden.vcf`. Mismatches on any field in
/// `hard_fields` fail the test; other CSQ fields are reported but soft.
///
/// Returns `Ok(true)` on a successful comparison and `Ok(false)` when
/// fixtures are missing (callers should `assert!(ran)` so a clean
/// checkout fails informatively rather than silently passing).
pub async fn run_and_compare_csq(name: &str, hard_fields: &[&str]) -> Result<bool> {
    let input_vcf = workspace_path(&format!(
        "datafusion/bio-function-vep/tests/data/port/{name}/input.vcf"
    ));
    let golden_vcf = workspace_path(&format!(
        "datafusion/bio-function-vep/tests/data/port/{name}/golden.vcf"
    ));
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");

    if !input_vcf.exists()
        || !golden_vcf.exists()
        || !cache_path.exists()
        || !ref_fasta.exists()
        || is_lfs_pointer(&input_vcf)
        || is_lfs_pointer(&golden_vcf)
    {
        eprintln!(
            "Skipping port case '{name}': fixture(s) missing or LFS-stubbed \
             (input: {} | golden: {} | cache: {} | ref: {})",
            input_vcf.display(),
            golden_vcf.display(),
            cache_path.display(),
            ref_fasta.display()
        );
        return Ok(false);
    }

    let ref_fasta_str = ref_fasta.to_str().unwrap().to_string();
    let config = base_config(&ref_fasta_str);
    compare_csq_inner(
        input_vcf.to_str().unwrap(),
        golden_vcf.to_str().unwrap(),
        cache_path.to_str().unwrap(),
        &config,
        hard_fields,
        None,
        name,
    )
    .await
}

/// (C2) Same as [`run_and_compare_csq`] but reads the golden from a
/// per-flag-set subfolder (`tests/data/port/<name>/<flag_subdir>/golden.vcf`)
/// and applies `cli_flags` to the engine config. Used by tests where one
/// fixture has multiple goldens (e.g., `--pick`, `--per_gene`).
///
/// Supported flag keys (see [`apply_cli_flag`] for the canonical list):
/// `pick`, `pick_allele`, `pick_allele_gene`, `per_gene`, `gencode_basic`,
/// `gencode_primary`, `exclude_predicted`, `distance`. Unknown flags are
/// reported via `eprintln!` and skipped — extending the supported set is
/// a per-port decision, not Batch 0 work.
pub async fn run_and_compare_csq_with_flags(
    name: &str,
    hard_fields: &[&str],
    flag_subdir: &str,
    cli_flags: &[(&str, &str)],
) -> Result<bool> {
    let input_vcf = workspace_path(&format!(
        "datafusion/bio-function-vep/tests/data/port/{name}/input.vcf"
    ));
    let golden_vcf = workspace_path(&format!(
        "datafusion/bio-function-vep/tests/data/port/{name}/{flag_subdir}/golden.vcf"
    ));
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");

    if !input_vcf.exists()
        || !golden_vcf.exists()
        || !cache_path.exists()
        || !ref_fasta.exists()
        || is_lfs_pointer(&input_vcf)
        || is_lfs_pointer(&golden_vcf)
    {
        eprintln!("Skipping port case '{name}/{flag_subdir}': fixture(s) missing or LFS-stubbed");
        return Ok(false);
    }

    let ref_fasta_str = ref_fasta.to_str().unwrap().to_string();
    let mut config = base_config(&ref_fasta_str);
    for (k, v) in cli_flags {
        apply_cli_flag(&mut config, k, v);
    }

    let label = format!("{name}/{flag_subdir}");
    compare_csq_inner(
        input_vcf.to_str().unwrap(),
        golden_vcf.to_str().unwrap(),
        cache_path.to_str().unwrap(),
        &config,
        hard_fields,
        None,
        &label,
    )
    .await
}

/// (C3) Same as [`run_and_compare_csq`] but filters CSQ allele-groups by
/// a predicate before comparison. Used by `Cache_RegFeat` to limit the
/// comparison to regulatory/motif rows (ignoring transcript rows that
/// belong to other ports' responsibility).
///
/// The predicate receives the per-token list of one CSQ allele-group
/// (e.g. `["G", "regulatory_region_variant", "MODIFIER", "", ...]`) and
/// returns `true` to include the group, `false` to skip it. The token
/// order is the CSQ Format declared in the golden header.
pub async fn run_and_compare_csq_filtered<F>(
    name: &str,
    hard_fields: &[&str],
    row_filter: F,
) -> Result<bool>
where
    F: Fn(&[String]) -> bool + Sync,
{
    let input_vcf = workspace_path(&format!(
        "datafusion/bio-function-vep/tests/data/port/{name}/input.vcf"
    ));
    let golden_vcf = workspace_path(&format!(
        "datafusion/bio-function-vep/tests/data/port/{name}/golden.vcf"
    ));
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");

    if !input_vcf.exists()
        || !golden_vcf.exists()
        || !cache_path.exists()
        || !ref_fasta.exists()
        || is_lfs_pointer(&input_vcf)
        || is_lfs_pointer(&golden_vcf)
    {
        eprintln!("Skipping port case '{name}' (filtered): fixture(s) missing or LFS-stubbed");
        return Ok(false);
    }

    let ref_fasta_str = ref_fasta.to_str().unwrap().to_string();
    let config = base_config(&ref_fasta_str);
    compare_csq_inner(
        input_vcf.to_str().unwrap(),
        golden_vcf.to_str().unwrap(),
        cache_path.to_str().unwrap(),
        &config,
        hard_fields,
        Some(&row_filter),
        name,
    )
    .await
}
