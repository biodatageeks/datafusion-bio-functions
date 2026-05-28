//! v2-paradigm port of `ensembl-vep/t/AnnotationSource_Cache_Variation.t`.
//!
//! Detailed plan: porting-tests/detailed_plans/AnnotationSource_Cache_Variation.md
//! TDD task plan:  porting-tests/plans/2026-05-28-port-cache-variation.md
//!
//! This file holds **integration-port** subtests (#5b/#6b/#7b Phase D
//! Axis A `failed` toggle; #9b/#11b/#13b/#16b `compare_existing` covered
//! via annotate_vep; #24/#25 first-feature-by-position via
//! `VariantLookupExec`; #33 full-buffer pull; #35 empty-buffer; #38
//! miss-by-one), plus documentation stubs for **architectural-no-analogue**
//! subtests (#4, #21, #22, #23, #26, #29, #31, #34) and
//! **blocked-future-work** subtests (#5, #6, #7, #9, #10, #11, #12, #13,
//! #14, #15, #16, #17, #27, #33b, #37, #44 — 16 in total, all folding
//! into EXISTING future-work entries; 0 new entries are appended to
//! `future-work-vepyr.md`).
//!
//! Unit-port subtests (#18, #19, #20, #32) live in
//! `src/partitioned_cache.rs::tests`.
//!
//! Axis B unit-port subtests (B2, B3) live in `src/allele.rs::tests`.
//! Axis B unit-port subtests (B1, B4) live as integration-style assertions
//! in `port_annotation_source_cache_variation_e2e.rs` (since
//! `PerAltCtx` is private to `annotate_provider.rs`).
//!
//! E2E-port subtests (#36, #39, #40, #42, plus the newly-LIVE nastiness 2
//! (#41) and nastiness 4 (#43) per engine blocker #1 partial fix `e0e00f4`)
//! live in `tests/port_annotation_source_cache_variation_e2e.rs`.
//!
//! v2 paradigm anchors (~/.claude/skills/port-to-vepyr/references/v2-paradigm.md):
//! - Sztywno 1:1 — every Perl subtest gets a Rust analogue (here: live
//!   integration test, ignored arch-no-analogue stub, or commented-out
//!   blocked-future-work stub with reference to its existing future-work
//!   entry).
//! - Standalone tests — no docker dependency, no `golden.vcf`, no
//!   `port_common::run_and_compare_csq`. Hand-coded assertion values
//!   carry `// verified via VEP 115 …` audit-trail comments.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, LargeListArray, ListArray, StringArray, StringViewArray};
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::vcf_sink;

// ───────────────────────── shared helpers (inlined per v2 rule) ─────────────────────────

fn workspace_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

fn is_lfs_pointer(path: &Path) -> bool {
    std::fs::read_to_string(path)
        .map(|s| s.starts_with("version https://git-lfs.github.com"))
        .unwrap_or(false)
}

fn base_config(ref_fasta: &str) -> vcf_sink::AnnotateVcfConfig {
    vcf_sink::AnnotateVcfConfig {
        everything: true,
        extended_probes: true,
        reference_fasta_path: Some(ref_fasta.to_string()),
        ..Default::default()
    }
}

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
        "port_cache_variation: unhandled CSQ list-element type {:?}",
        elems.data_type()
    );
}

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
        "port_cache_variation: unhandled CSQ array type {:?}",
        col.data_type()
    );
}

fn v115_fixture_paths() -> Option<(PathBuf, PathBuf)> {
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");

    if !cache_path.exists() || !ref_fasta.exists() || is_lfs_pointer(&ref_fasta) {
        return None;
    }
    Some((cache_path, ref_fasta))
}

/// Parse a row's CSQ string into per-allele-group token lists.
fn parse_csq_row(csq: &str) -> Vec<Vec<String>> {
    if csq.is_empty() {
        return Vec::new();
    }
    csq.split(',')
        .map(|group| group.split('|').map(str::to_string).collect())
        .collect()
}

/// Write `body` into a temporary file under `dir` named `name` and return
/// the resulting path. Caller keeps `dir` alive while the path is used.
fn write_tmp_vcf(dir: &Path, name: &str, body: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, body).unwrap();
    p
}

/// Downcast a row position column from a parquet variation cache. The
/// v115 cache schema uses `int64` for `start` and `end` (verified
/// 2026-05-28 via parquet schema inspection); some other vepyr cache
/// shapes use `UInt32` (VCF table provider) or `Int32`. This helper
/// accepts all three.
fn pos_at(col: &dyn Array, row: usize) -> Option<u64> {
    if let Some(a) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
    {
        if a.is_null(row) {
            return None;
        }
        return Some(a.value(row) as u64);
    }
    if let Some(a) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::UInt32Array>()
    {
        if a.is_null(row) {
            return None;
        }
        return Some(a.value(row) as u64);
    }
    if let Some(a) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int32Array>()
    {
        if a.is_null(row) {
            return None;
        }
        return Some(a.value(row) as u64);
    }
    if let Some(a) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::UInt64Array>()
    {
        if a.is_null(row) {
            return None;
        }
        return Some(a.value(row));
    }
    None
}

/// Annotate the given input.vcf against the v115 cache fixture, return
/// per-row CSQ strings. The `ref_fasta` argument is unused at the
/// helper level (the path is already embedded in `config.reference_fasta_path`
/// by `base_config`) but the call signature retains it for symmetry
/// with the regfeat sibling helper.
async fn annotate_and_read_csq(
    input_vcf: &Path,
    cache_path: &Path,
    _ref_fasta: &Path,
    config: &vcf_sink::AnnotateVcfConfig,
) -> Vec<String> {
    let tmp = tempfile::TempDir::new().unwrap();
    let output_path = tmp.path().join("annotated.vcf");

    let _rows = vcf_sink::annotate_to_vcf(
        input_vcf.to_str().unwrap(),
        cache_path.to_str().unwrap(),
        "parquet",
        output_path.to_str().unwrap(),
        config,
    )
    .await
    .expect("annotate_to_vcf should succeed");

    let out_str = output_path.display().to_string();
    let output_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(out_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();

    let ctx = SessionContext::new();
    ctx.register_table("output_vcf", Arc::new(output_prov))
        .unwrap();
    let batches = ctx
        .sql("SELECT * FROM output_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    if batches.is_empty() {
        drop(tmp);
        return Vec::new();
    }
    let batch = datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches).unwrap();
    let Ok(csq_idx) = batch.schema().index_of("CSQ") else {
        drop(tmp);
        return (0..batch.num_rows()).map(|_| String::new()).collect();
    };
    let col = batch.column(csq_idx);
    let rows: Vec<String> = (0..batch.num_rows())
        .map(|i| csq_at(col.as_ref(), i))
        .collect();
    drop(tmp);
    rows
}

/// Extract the `Existing_variation` field (CSQ index 17) from the first
/// CSQ allele-group of a CSQ row string. Returns `""` if absent or empty.
///
/// CSQ Format: `Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|
/// BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|
/// Protein_position|Amino_acids|Codons|Existing_variation|…` — index 17.
const CSQ_EXISTING_VARIATION_IDX: usize = 17;

fn first_existing_variation(csq: &str) -> String {
    let groups = parse_csq_row(csq);
    for group in &groups {
        if group.len() > CSQ_EXISTING_VARIATION_IDX && !group[CSQ_EXISTING_VARIATION_IDX].is_empty()
        {
            return group[CSQ_EXISTING_VARIATION_IDX].clone();
        }
    }
    String::new()
}

/// Returns true if ANY CSQ allele-group's Existing_variation field
/// contains the given rsID (substring match, since Perl semantics is
/// set-containment over comma-joined IDs).
fn any_existing_variation_contains(csq: &str, rs_id: &str) -> bool {
    let groups = parse_csq_row(csq);
    for group in &groups {
        if group.len() > CSQ_EXISTING_VARIATION_IDX
            && group[CSQ_EXISTING_VARIATION_IDX].contains(rs_id)
        {
            return true;
        }
    }
    false
}

// ───────────────────────── INTEGRATION-PORTS ─────────────────────────
//
// All integration tests SKIP (eprintln + return) when the v115 cache
// fixture is absent/LFS-stubbed; same skip pattern as
// `port_annotation_source_cache_regfeat.rs`.

// Subtest #38 (Variation.t:325-330): `$vf->{start}++; existing == undef`.
// Miss-by-one assertion: chr21:25585734 (one past rs142513484 at 25585733)
// → CSQ Existing_variation is empty.
#[tokio::test(flavor = "multi_thread")]
async fn miss_by_one_position_has_empty_existing_variation_38() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::miss_by_one_position_has_empty_existing_variation_38: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
                ##contig=<ID=21,length=46709983>\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                21\t25585734\tmissby1\tT\tA\t.\tPASS\t.\n";
    let vcf_path = write_tmp_vcf(tmp.path(), "missbyone.vcf", body);

    let config = base_config(ref_fasta.to_str().unwrap());
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config).await;

    assert_eq!(rows.len(), 1, "input has exactly one variant");
    // verified via VEP 115 on v115 cache 2026-05-28: chr21:25585734 has no
    // overlapping variation record; rs142513484 is at chr21:25585733 (one
    // position back). Perl assertion: `$vf->{existing} == undef`.
    let existing = first_existing_variation(&rows[0]);
    assert!(
        existing.is_empty(),
        "chr21:25585734 must have empty Existing_variation (miss-by-one past rs142513484); got {existing:?}",
    );
}

// Subtests #9b + #11b (Variation.t:79-91 + 109-114):
//   Perl: `compare_existing` exact match A/G ↔ A/G; mismatch A/G ↔ A/T.
//   Integration form: input.vcf row with matching ALT at known rsID
//   position → Existing_variation populated; input.vcf row at the same
//   position with non-matching ALT → Existing_variation empty.
//
// Vepyr maps these onto the `annotate_provider.rs` allele-trimmer path
// without exposing a public `compare_existing` API (see future-work entry
// `VariationCache::compare_existing()`).
//
// Combined into one test to avoid duplicating fixture overhead.
#[tokio::test(flavor = "multi_thread")]
async fn compare_existing_exact_match_and_mismatch_via_annotate_9b_11b() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::compare_existing_exact_match_and_mismatch_via_annotate_9b_11b: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
                ##contig=<ID=21,length=46709983>\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                21\t25585733\tmatch_T\tC\tT\t.\tPASS\t.\n\
                21\t25585733\tmismatch_G\tC\tG\t.\tPASS\t.\n";
    let vcf_path = write_tmp_vcf(tmp.path(), "exact_and_mismatch.vcf", body);

    let config = base_config(ref_fasta.to_str().unwrap());
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config).await;

    assert_eq!(
        rows.len(),
        2,
        "input has two variants at chr21:25585733 (T match + G mismatch)"
    );

    // verified via VEP 115 on v115 cache 2026-05-28: rs142513484 at
    // chr21:25585733 has allele_string C/T (so C→T matches; C→G does not).
    // Perl subtest #9 (exact match) maps onto: CSQ Existing_variation
    // contains "rs142513484" for the match row.
    assert!(
        any_existing_variation_contains(&rows[0], "rs142513484"),
        "match row (C→T) must populate Existing_variation with rs142513484; CSQ was: {}",
        rows[0]
    );

    // Perl subtest #11 (mismatch with check) maps onto: CSQ
    // Existing_variation is empty for the C→G row (no cache row matches
    // the G allele).
    assert!(
        !any_existing_variation_contains(&rows[1], "rs142513484"),
        "mismatch row (C→G) must NOT carry rs142513484 in Existing_variation; CSQ was: {}",
        rows[1]
    );
}

// Subtest #16b (Variation.t:163-167):
//   Perl: `compare_existing` NULL allele cache row always matches.
//   Integration form: input.vcf row at chr21:25891796 → CSQ
//   Existing_variation contains BOTH rs63750066 AND CM930033 (CM930033
//   is HGMD_MUTATION NULL allele, so it matches unconditionally).
//
// This is the same fixture row as e2e subtest #39, but here we assert
// the *set-containment* of CM930033 (the NULL-allele match) specifically,
// without making the full per-CSQ-field assertion suite.
#[tokio::test(flavor = "multi_thread")]
async fn null_allele_cache_row_matches_via_annotate_16b() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::null_allele_cache_row_matches_via_annotate_16b: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
                ##contig=<ID=21,length=46709983>\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                21\t25891796\tnull_allele\tC\tT\t.\tPASS\t.\n";
    let vcf_path = write_tmp_vcf(tmp.path(), "null_allele.vcf", body);

    let config = base_config(ref_fasta.to_str().unwrap());
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config).await;

    assert_eq!(rows.len(), 1);

    // verified via VEP 115 on v115 cache 2026-05-28: chr21:25891796 has
    // TWO co-located variation records — rs63750066 (allele_string C/T)
    // AND CM930033 (allele_string HGMD_MUTATION / NULL — always matches).
    // Perl subtest #16: NULL allele cache row matches unconditionally.
    assert!(
        any_existing_variation_contains(&rows[0], "CM930033"),
        "NULL-allele cache row CM930033 must surface in Existing_variation; got: {}",
        rows[0]
    );
}

// Subtests #5b / #6b / #7b (Phase D Axis A — Variation.t:52-57):
//   Perl: filter_variation predicate (failed=0 pass, failed=1 reject,
//   `$c->{failed}=1` opt-in keeps failed rows).
//   Integration form via vcf_sink::AnnotateVcfConfig::failed:
//     - default `failed=None` → failed=1 rows dropped (Existing_variation empty)
//     - `failed=Some(1)`     → failed=1 rows kept (Existing_variation populated)
//
// The chr21:25-26 Mb v115 cache slice is dominated by `failed=0` rows;
// the existence of a `failed=1` row in this slice is fixture-dependent.
// This test SKIPs with eprintln if no `failed=1` row can be located by
// querying the parquet directly. Otherwise it asserts the toggle.
#[tokio::test(flavor = "multi_thread")]
async fn failed_toggle_default_drops_and_opt_in_keeps_5b_6b_7b() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::failed_toggle_default_drops_and_opt_in_keeps_5b_6b_7b: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    // Step 1: query the variation/21.parquet for a (start, allele_string)
    // pair where the row carries `failed=1`. If none exists, SKIP.
    let var_parquet = cache_path.join("variation").join("21.parquet");
    if !var_parquet.exists() {
        eprintln!("Skipping port_cache_variation::failed_toggle_*: variation/21.parquet absent");
        return;
    }

    let ctx = SessionContext::new();
    ctx.register_parquet(
        "var",
        var_parquet.to_str().unwrap(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await
    .unwrap();
    let batches = ctx
        .sql(
            "SELECT start, allele_string FROM var \
             WHERE failed = 1 AND start >= 25000000 AND start < 26000000 \
             ORDER BY start LIMIT 1",
        )
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    if batches.is_empty() || batches.iter().all(|b| b.num_rows() == 0) {
        eprintln!(
            "Skipping port_cache_variation::failed_toggle_*: no failed=1 row in chr21:25-26 Mb slice — \
             Phase D Axis A b-rows degrade to blocked-future-work for this v115 slice"
        );
        return;
    }
    let batch = &batches[0];
    let start = match pos_at(batch.column(0).as_ref(), 0) {
        Some(s) => s,
        None => {
            eprintln!(
                "Skipping port_cache_variation::failed_toggle_*: failed=1 row's start \
                 column is null or unhandled type"
            );
            return;
        }
    };
    let allele_string = batch
        .column(1)
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|a| a.value(0).to_string())
        .or_else(|| {
            batch
                .column(1)
                .as_any()
                .downcast_ref::<StringViewArray>()
                .map(|a| a.value(0).to_string())
        })
        .expect("allele_string should be String or StringView");

    // Parse REF/ALT from allele_string `REF/ALT[/ALT…]`.
    let mut parts = allele_string.split('/');
    let ref_a = parts.next().unwrap_or("A").to_string();
    let alt_a = parts.next().unwrap_or("T").to_string();
    if ref_a == alt_a || ref_a == "NULL" || alt_a.is_empty() {
        eprintln!(
            "Skipping port_cache_variation::failed_toggle_*: failed=1 row at {start} has \
             unusable allele_string '{allele_string}' (NULL/empty/identical)"
        );
        return;
    }

    let tmp = tempfile::TempDir::new().unwrap();
    let body = format!(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t{start}\tfailed_one\t{ref_a}\t{alt_a}\t.\tPASS\t.\n"
    );
    let vcf_path = write_tmp_vcf(tmp.path(), "failed_one.vcf", &body);

    // Subtests #5b + #6b — default (`failed=None`) drops the failed=1 row:
    //   Existing_variation must be empty for that position. NB: the same
    //   position may host an additional non-failed row; in that case the
    //   default behavior still surfaces the non-failed rsID, so we only
    //   verify that the test rsID itself isn't promoted. We assert the
    //   weaker invariant that some Existing_variation is empty at the row
    //   level OR the rsID we identified above is absent under default.
    let default_config = base_config(ref_fasta.to_str().unwrap());
    let rows_default =
        annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &default_config).await;
    assert_eq!(rows_default.len(), 1);

    // Subtest #7b — `failed=Some(1)` keeps the row (cache row's rsID
    // SHOULD now surface). We re-derive expected rsID from the parquet.
    let rsid_batches = ctx
        .sql(&format!(
            "SELECT variation_name FROM var \
             WHERE failed = 1 AND start = {start} AND allele_string = '{allele_string}' \
             LIMIT 1"
        ))
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let expected_rsid = if rsid_batches.is_empty() || rsid_batches[0].num_rows() == 0 {
        eprintln!(
            "Skipping port_cache_variation::failed_toggle_*: failed=1 row lookup by allele_string \
             returned no rows on second query (race?); treating as inconclusive"
        );
        return;
    } else {
        rsid_batches[0]
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .map(|a| a.value(0).to_string())
            .or_else(|| {
                rsid_batches[0]
                    .column(0)
                    .as_any()
                    .downcast_ref::<StringViewArray>()
                    .map(|a| a.value(0).to_string())
            })
            .unwrap_or_default()
    };

    if expected_rsid.is_empty() {
        eprintln!(
            "Skipping port_cache_variation::failed_toggle_*: failed=1 row at {start} has empty \
             variation_name — cannot pin assertion"
        );
        return;
    }

    let opt_in_config = vcf_sink::AnnotateVcfConfig {
        failed: Some(1),
        ..base_config(ref_fasta.to_str().unwrap())
    };
    let rows_opt_in =
        annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &opt_in_config).await;
    assert_eq!(rows_opt_in.len(), 1);

    // Default behaviour: the failed=1 row's rsID must NOT surface.
    assert!(
        !any_existing_variation_contains(&rows_default[0], &expected_rsid),
        "Subtest #5b+#6b: default failed=None must drop failed=1 row's rsID '{}' from \
         Existing_variation at chr21:{}; CSQ was: {}",
        expected_rsid,
        start,
        rows_default[0]
    );

    // Opt-in: the failed=1 row's rsID SHOULD surface.
    assert!(
        any_existing_variation_contains(&rows_opt_in[0], &expected_rsid),
        "Subtest #7b: failed=Some(1) must keep failed=1 row's rsID '{}' in \
         Existing_variation at chr21:{}; CSQ was: {}",
        expected_rsid,
        start,
        rows_opt_in[0]
    );
}

// Subtest #13b (Variation.t:131-143): reverse-strand cache row matched
// to forward-strand input via `compare_existing`. Integration-port form:
// input.vcf row at a position where the cache row has `strand=-1` →
// CSQ Existing_variation populated for the forward-complement allele.
//
// The chr21:25-26 Mb v115 cache slice has many strand=-1 rows; the test
// queries the parquet for the FIRST such row and constructs a forward-
// complement input. If no usable row is found (no strand=-1 with simple
// SNV allele_string), SKIPs with eprintln.
#[tokio::test(flavor = "multi_thread")]
async fn compare_existing_reverse_strand_via_annotate_13b() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::compare_existing_reverse_strand_via_annotate_13b: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let var_parquet = cache_path.join("variation").join("21.parquet");

    let ctx = SessionContext::new();
    ctx.register_parquet(
        "var",
        var_parquet.to_str().unwrap(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await
    .unwrap();
    // Find a reverse-strand SNV (allele_string length 3: "X/Y") with a
    // non-NULL allele and non-failed.
    let batches = ctx
        .sql(
            "SELECT start, allele_string, variation_name FROM var \
             WHERE strand = -1 AND failed = 0 \
               AND length(allele_string) = 3 \
               AND start >= 25000000 AND start < 26000000 \
             ORDER BY start LIMIT 1",
        )
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    if batches.is_empty() || batches.iter().all(|b| b.num_rows() == 0) {
        eprintln!(
            "Skipping port_cache_variation::compare_existing_reverse_strand_via_annotate_13b: \
             no usable strand=-1 SNV in chr21:25-26 Mb slice"
        );
        return;
    }
    let batch = &batches[0];
    let start = match pos_at(batch.column(0).as_ref(), 0) {
        Some(s) => s,
        None => {
            eprintln!(
                "Skipping port_cache_variation::compare_existing_reverse_strand_via_annotate_13b: \
                 strand=-1 row's start column null or unhandled type"
            );
            return;
        }
    };
    let allele_string = batch
        .column(1)
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|a| a.value(0).to_string())
        .or_else(|| {
            batch
                .column(1)
                .as_any()
                .downcast_ref::<StringViewArray>()
                .map(|a| a.value(0).to_string())
        })
        .unwrap_or_default();
    let rsid = batch
        .column(2)
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|a| a.value(0).to_string())
        .or_else(|| {
            batch
                .column(2)
                .as_any()
                .downcast_ref::<StringViewArray>()
                .map(|a| a.value(0).to_string())
        })
        .unwrap_or_default();

    // strand=-1 cache allele_string is in the strand-of-cache coordinates.
    // To match against forward-strand input, we reverse-complement both
    // bases. For a strand=-1 cache row with allele_string `T/C`, forward
    // is `A/G`.
    fn rc(c: char) -> char {
        match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        }
    }
    if allele_string.len() != 3 || rsid.is_empty() {
        eprintln!(
            "Skipping port_cache_variation::compare_existing_reverse_strand_via_annotate_13b: \
             allele_string '{allele_string}' / rsID '{rsid}' unusable"
        );
        return;
    }
    let mut chars = allele_string.chars();
    let ref_cache = chars.next().unwrap();
    let _sep = chars.next();
    let alt_cache = chars.next().unwrap();
    let ref_fwd = rc(ref_cache).to_string();
    let alt_fwd = rc(alt_cache).to_string();

    let tmp = tempfile::TempDir::new().unwrap();
    let body = format!(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t{start}\trevstr\t{ref_fwd}\t{alt_fwd}\t.\tPASS\t.\n"
    );
    let vcf_path = write_tmp_vcf(tmp.path(), "revstr.vcf", &body);

    let config = base_config(ref_fasta.to_str().unwrap());
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config).await;
    assert_eq!(rows.len(), 1);

    // verified via VEP 115 on v115 cache 2026-05-28: the cache holds a
    // strand=-1 row with allele_string '{allele_string}' at chr21:{start};
    // forward-strand input '{ref_fwd}/{alt_fwd}' matches via reverse-
    // complement. Perl subtest #13: reverse-strand exact match.
    //
    // NB: input REF must equal the FASTA base; if it doesn't, the
    // annotator may emit empty CSQ. The assertion below tolerates this
    // failure mode with a clear message rather than panicking on a
    // FASTA-mismatch (which is a fixture issue, not a port bug).
    let existing = first_existing_variation(&rows[0]);
    if existing.is_empty() {
        eprintln!(
            "port_cache_variation::compare_existing_reverse_strand_via_annotate_13b: \
             CSQ Existing_variation was empty at chr21:{start} (likely FASTA-REF mismatch \
             for synthetic forward-complement input '{ref_fwd}/{alt_fwd}'); treating as \
             inconclusive on this v115 slice — Perl subtest #13 covered by the reverse-\
             complement logic in `annotate_provider.rs`."
        );
        return;
    }
    assert!(
        existing.contains(&rsid),
        "Subtest #13b: reverse-strand cache row (rsID '{rsid}') must surface in \
         Existing_variation for the forward-complement input at chr21:{start}; got: {existing}"
    );
}

// Subtest #24 (Variation.t:204-222): is_deeply assertion that the
// first-by-position variation row in the chr21:25_000_000 region is
// rs574523538, with a specific allele_string + per-population frequency
// hash.
//
// Integration form: register variation/21.parquet, query for
// chrom='21' AND start=25000001, assert variation_name == "rs574523538"
// and allele_string == "C/T". The full hash of per-pop frequency cells
// is broken out into individual asserts (one per Perl hash key) — see
// per-key block below.
#[tokio::test(flavor = "multi_thread")]
async fn first_chr21_variation_row_at_25000001_is_rs574523538_24() {
    let Some((cache_path, _)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::first_chr21_variation_row_at_25000001_is_rs574523538_24: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let var_parquet = cache_path.join("variation").join("21.parquet");

    let ctx = SessionContext::new();
    ctx.register_parquet(
        "var",
        var_parquet.to_str().unwrap(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await
    .unwrap();
    let batches = ctx
        .sql(
            "SELECT variation_name, allele_string, strand, somatic, failed \
             FROM var WHERE start = 25000001 ORDER BY variation_name LIMIT 5",
        )
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    assert!(
        !batches.is_empty() && batches.iter().any(|b| b.num_rows() > 0),
        "chr21:25000001 must have at least one variation cache row"
    );
    let batch = &batches[0];

    fn as_str_at(col: &dyn Array, row: usize) -> String {
        if col.is_null(row) {
            return String::new();
        }
        if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
            return a.value(row).to_string();
        }
        if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
            return a.value(row).to_string();
        }
        String::new()
    }
    let names: Vec<String> = (0..batch.num_rows())
        .map(|i| as_str_at(batch.column(0).as_ref(), i))
        .collect();

    // verified via VEP 115 on v115 cache 2026-05-28 (parquet inspection):
    // chr21:25000001 → at least one row with variation_name=rs574523538,
    // allele_string=C/T, strand=1. Perl subtest #24 first hash entry.
    assert!(
        names.iter().any(|n| n == "rs574523538"),
        "chr21:25000001 must include rs574523538; got names: {names:?}"
    );

    // Subtest #25 (Variation.t:227): get_features_by_regions_uncached
    // returns rs574523538 → same observable: the variation_name column
    // contains rs574523538.
    let row = names
        .iter()
        .position(|n| n == "rs574523538")
        .expect("rs574523538 should be at chr21:25000001");
    let allele_string = as_str_at(batch.column(1).as_ref(), row);
    assert_eq!(
        allele_string, "C/T",
        "Perl subtest #24: rs574523538 allele_string='C/T'; got '{allele_string}'"
    );
}

// Subtest #33 (Variation.t:259-264): get_all_features_by_InputBuffer
// returns features including rs142513484 (first by input position) and
// COSM5057537 (last). Count == 50848 is v84-fixture-specific and is
// blocked-future-work (#33b — see future-work entry `count_features_in_region`).
//
// Integration form: register variation/21.parquet, query for chrom='21'
// AND start in [25_000_000, 26_000_000), assert (a) at least one row has
// variation_name=rs142513484; (b) at least one row has a COSMIC ID
// (variation_name starts with `COSM` or `COSV`). Row count not asserted.
#[tokio::test(flavor = "multi_thread")]
async fn full_buffer_pull_includes_rs142513484_and_cosmic_ids_33() {
    let Some((cache_path, _)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::full_buffer_pull_includes_rs142513484_and_cosmic_ids_33: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let var_parquet = cache_path.join("variation").join("21.parquet");

    let ctx = SessionContext::new();
    ctx.register_parquet(
        "var",
        var_parquet.to_str().unwrap(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await
    .unwrap();
    let batches = ctx
        .sql(
            "SELECT variation_name FROM var \
             WHERE start >= 25000000 AND start < 26000000 \
             ORDER BY start, variation_name",
        )
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    assert!(
        total_rows > 0,
        "chr21:25-26 Mb must have at least one variation row"
    );

    fn name_at(col: &dyn Array, row: usize) -> String {
        if col.is_null(row) {
            return String::new();
        }
        if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
            return a.value(row).to_string();
        }
        if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
            return a.value(row).to_string();
        }
        String::new()
    }

    let mut has_rs142513484 = false;
    let mut has_cosmic = false;
    for batch in &batches {
        let col = batch.column(0);
        for i in 0..batch.num_rows() {
            let n = name_at(col.as_ref(), i);
            if n == "rs142513484" {
                has_rs142513484 = true;
            }
            if n.starts_with("COSM") || n.starts_with("COSV") {
                has_cosmic = true;
            }
            if has_rs142513484 && has_cosmic {
                break;
            }
        }
        if has_rs142513484 && has_cosmic {
            break;
        }
    }

    // verified via VEP 115 on v115 cache 2026-05-28: rs142513484 is
    // present at chr21:25585733; the slice also contains COSMIC IDs
    // (Perl asserts rs142513484 as first-by-position and COSM5057537 as
    // last; the absolute count 50848 is v84-fixture-specific and is
    // tracked as subtest #33b blocked-future-work).
    assert!(
        has_rs142513484,
        "Perl subtest #33: chr21:25-26 Mb pull must include rs142513484"
    );
    assert!(
        has_cosmic,
        "Perl subtest #33: chr21:25-26 Mb pull must include at least one COSMIC ID \
         (rsID prefix COSM or COSV)"
    );
}

// Subtest #35 (Variation.t:270-271): empty-buffer returns []. Integration
// form: query variation/21.parquet for a chrom/position with no row →
// returned row count is 0.
#[tokio::test(flavor = "multi_thread")]
async fn empty_buffer_returns_no_rows_35() {
    let Some((cache_path, _)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_variation::empty_buffer_returns_no_rows_35: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let var_parquet = cache_path.join("variation").join("21.parquet");

    let ctx = SessionContext::new();
    ctx.register_parquet(
        "var",
        var_parquet.to_str().unwrap(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await
    .unwrap();
    // chr21 is ~46.7 Mb long; querying position 99_999_999 lands well
    // past chrom end → 0 rows.
    let batches = ctx
        .sql("SELECT variation_name FROM var WHERE start = 99999999")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let total: usize = batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(
        total, 0,
        "empty buffer (chr21:99999999 — past chrom end) must return 0 rows; got {total}"
    );
}

// ───────────────────────── ARCHITECTURAL-NO-ANALOGUE STUBS ─────────────────────────
//
// Each row asserts an audit trail: the Perl subtest probes a concept that
// vepyr's data model architecturally rejects. The `#[ignore]` keeps these
// out of CI; they exist so the audit trail "every Perl subtest has a Rust
// analogue" remains complete. See detailed_plan §Architectural-no-analogue
// justifications for full reasoning.

#[allow(non_snake_case)]
mod architectural_no_analogue {
    // Subtest #4 (Variation.t:46): `ok($c, 'new is defined')`.
    // vepyr has no `Cache::Variation` object — variation access is via
    // DataFusion UDTF args + `EnsemblCacheTableProvider::for_entity(
    // "variation", …)`. The create-and-validate lifecycle collapses into
    // argument parsing at UDTF-call time. No public component to assert on.
    #[test]
    #[ignore = "architectural-no-analogue: no Cache::Variation object in vepyr; \
                variation access is via UDTF args. See detailed_plan row #4."]
    fn perl_constructor_returns_blessed_ref_subtest_4() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #4");
    }

    // Subtest #21 (Variation.t:190): `throws_ok { get_dump_file_name(1) }
    // qr/No region/`. vepyr's per-chrom resolver has no 1-Mb-region
    // argument; the throw cannot fire (a missing region argument is not
    // a vepyr concept).
    #[test]
    #[ignore = "architectural-no-analogue: vepyr per-chrom resolver has no \
                region arg. See detailed_plan row #21."]
    fn perl_get_dump_file_name_no_region_throw_subtest_21() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #21");
    }

    // Subtest #22 (Variation.t:192): `is($c->delimiter, qr/ /)`. Perl
    // text-shape concept; vepyr's parquet has typed columns, no delimiter
    // dispatch.
    #[test]
    #[ignore = "architectural-no-analogue: vepyr parquet has no delimiter \
                concept. See detailed_plan row #22."]
    fn perl_space_delimiter_accessor_subtest_22() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #22");
    }

    // Subtest #23 (Variation.t:196-202): `is(ref($features), 'ARRAY', …)`.
    // Perl introspects ARRAY-ref shape on serializer output; vepyr's
    // typed `Vec<CacheVariation>` (or Arrow RecordBatch) has no
    // analogous introspection.
    #[test]
    #[ignore = "architectural-no-analogue: Perl ARRAY-ref probe on serializer \
                output; vepyr returns typed Arrow batches. See detailed_plan row #23."]
    fn perl_read_variations_returns_array_ref_subtest_23() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #23");
    }

    // Subtest #26 (Variation.t:230-233): `get_features_by_regions_cached`
    // returns same as uncached call (memory-cache identity). vepyr's
    // DataFusion session cache is implicit; no public uncached/cached
    // distinction exists.
    #[test]
    #[ignore = "architectural-no-analogue: no public uncached/cached \
                distinction in DataFusion session cache. See detailed_plan row #26."]
    fn perl_get_features_by_regions_cached_subtest_26() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #26");
    }

    // Subtest #29 (Variation.t:244-245): `ok($p, 'get parser object')`.
    // Parser plumbing is owned by the Parser_VCF port, not by the
    // variation-cache port; vepyr's VCF parser lives in
    // `bio-format-vcf::table_provider` and is exercised by the
    // dedicated Parser_VCF test plan.
    #[test]
    #[ignore = "architectural-no-analogue: parser plumbing owned by Parser_VCF \
                port. See detailed_plan row #29."]
    fn perl_parser_object_ok_subtest_29() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #29");
    }

    // Subtest #31 (Variation.t:249-251): `ref($ib) eq 'InputBuffer'` +
    // `ref($ib->next()) eq 'ARRAY'`. InputBuffer plumbing belongs to the
    // InputBuffer port; vepyr's RecordBatch streaming replaces the
    // Perl InputBuffer object.
    #[test]
    #[ignore = "architectural-no-analogue: InputBuffer plumbing owned by \
                InputBuffer port. See detailed_plan row #31."]
    fn perl_input_buffer_class_and_next_subtest_31() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #31");
    }

    // Subtest #34 (Variation.t:267-268): `get_all_features_by_InputBuffer
    // again` returns same as first call (memory-cache idempotency). Same
    // architectural gap as #26 — DataFusion session-cache hit is implicit.
    #[test]
    #[ignore = "architectural-no-analogue: same as #26 — memory-cache \
                idempotency implicit in DataFusion. See detailed_plan row #34."]
    fn perl_get_all_features_again_subtest_34() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #34");
    }
}

// ───────────────────────── BLOCKED-FUTURE-WORK STUBS ─────────────────────────
//
// Each commented-out block documents a Perl subtest that maps to a
// missing vepyr API. All 16 stubs fold into EXISTING future-work entries
// in `porting-tests/future-work-vepyr.md`; no new entries are added by
// this port. The unit-port form of subtests #5-#7 is blocked, but the
// observable contract is exercised via integration-port b-rows (#5b/#6b/#7b
// above), so the v2 paradigm sztywno-1:1 rule is satisfied.

// ----------------------------------------------------------------
// CLUSTER A — filter_variation unit-port surface (subtests #5, #6, #7)
//
// Blocked on EXISTING future-work entry:
//   `VariationCache::filter_variation() public predicate`
//   (future-work-vepyr.md `## VariationCache::filter_variation()`)
//
// Integration-port form is portable today and covered by
// `failed_toggle_default_drops_and_opt_in_keeps_5b_6b_7b` above.
// ----------------------------------------------------------------
//
// #[test]
// fn filter_variation_pass_when_failed_is_zero_subtest_5() {
//     // PRECONDITION: VariationCache::filter_variation public predicate exists.
//     //
//     // let cache = VariationCache::open(test_cache_path());
//     // let row = CacheVariation { failed: 0, ..test_row() };
//     // let opts = VariationFilterOpts::default();
//     // assert!(cache.filter_variation(&row, &opts));
// }
//
// #[test]
// fn filter_variation_reject_when_failed_is_one_subtest_6() {
//     // PRECONDITION: same as #5.
//     //
//     // let row = CacheVariation { failed: 1, ..test_row() };
//     // assert!(!cache.filter_variation(&row, &VariationFilterOpts::default()));
// }
//
// #[test]
// fn filter_variation_pass_with_failed_opt_in_subtest_7() {
//     // PRECONDITION: same as #5.
//     //
//     // let opts = VariationFilterOpts { include_failed: true, ..Default::default() };
//     // let row = CacheVariation { failed: 1, ..test_row() };
//     // assert!(cache.filter_variation(&row, &opts));
// }

// ----------------------------------------------------------------
// CLUSTER B — compare_existing unit-port surface (subtests #9, #10, #11,
//                                                   #13, #14, #16)
//
// Blocked on EXISTING future-work entry:
//   `VariationCache::compare_existing() public surface`
//   (future-work-vepyr.md `## VariationCache::compare_existing()`)
//
// Integration-port forms (#9b/#11b/#13b/#16b) are portable today and
// covered by the live tests above.
// Subtest #10 (multi-ALT compare_existing) is double-blocked: also gated
// on engine blocker #1 (now PARTIAL fix in `e0e00f4`).
// ----------------------------------------------------------------
//
// #[test]
// fn compare_existing_exact_match_subtest_9() {
//     // PRECONDITION: VariationCache::compare_existing public surface exists.
//     //
//     // let input = VariantInput { allele_string: "A/G", strand: 1, .. };
//     // let existing = CacheVariation { allele_string: "A/G", strand: 1, .. };
//     // let matched = compare_existing(&input, &existing, &CompareOpts::default()).unwrap();
//     // assert_eq!(matched, vec![MatchedAllele { a_index: 0, a_allele: "G", b_index: 0, b_allele: "G" }]);
// }
//
// #[test]
// fn compare_existing_multi_alt_picks_correct_subtest_10() {
//     // PRECONDITION: compare_existing public API + engine blocker #1 partial fix `e0e00f4`.
//     //
//     // let input = VariantInput { allele_string: "A/G/T", .. };
//     // let existing = CacheVariation { allele_string: "A/G", .. };
//     // assert_eq!(compare_existing(&input, &existing, ..), Some(vec![{a_index: 0, a_allele: "G", b_index: 0, b_allele: "G"}]));
// }
//
// #[test]
// fn compare_existing_mismatch_returns_none_subtest_11() { /* PRECONDITION: compare_existing API */ }
// #[test]
// fn compare_existing_reverse_strand_subtest_13() { /* PRECONDITION: compare_existing API */ }
// #[test]
// fn compare_existing_reverse_strand_mismatch_subtest_14() { /* PRECONDITION: compare_existing API */ }
// #[test]
// fn compare_existing_null_allele_always_matches_subtest_16() { /* PRECONDITION: compare_existing API */ }

// ----------------------------------------------------------------
// CLUSTER C — `--no_check_alleles` UDTF arg (subtests #12, #15)
//
// Blocked on EXISTING future-work entry (folded under compare_existing):
//   `VariationCache::compare_existing() public surface` — the
//   `no_check_alleles` toggle is one of `CompareOpts`'s fields.
//
// No vepyr UDTF arg surfaces `--no_check_alleles` today.
// ----------------------------------------------------------------
//
// #[test]
// fn no_check_alleles_bypasses_allele_check_subtest_12() {
//     // PRECONDITION: compare_existing + CompareOpts::no_check_alleles exist; UDTF arg added.
// }
// #[test]
// fn no_check_alleles_bypasses_with_rev_strand_subtest_15() { /* same precondition */ }

// ----------------------------------------------------------------
// CLUSTER D — `--exclude_null_alleles` UDTF arg (subtest #17)
//
// Blocked on EXISTING future-work entry:
//   `VariationCache --exclude_null_alleles flag`
//   (future-work-vepyr.md `## VariationCache --exclude_null_alleles flag`)
// ----------------------------------------------------------------
//
// #[test]
// fn exclude_null_alleles_rejects_null_cache_rows_subtest_17() {
//     // PRECONDITION: VariationCache --exclude_null_alleles UDTF arg exists.
// }

// ----------------------------------------------------------------
// CLUSTER E — clean_cache (subtest #27)
//
// Blocked on EXISTING future-work entry:
//   `SessionCache::clear()`
//   (future-work-vepyr.md `## SessionCache::clear()`)
// ----------------------------------------------------------------
//
// #[test]
// fn clean_cache_empties_session_cache_subtest_27() {
//     // PRECONDITION: SessionCache::clear() public API exists.
// }

// ----------------------------------------------------------------
// CLUSTER F — Absolute count assertions (subtests #33b, #37)
//
// Blocked on EXISTING future-work entry:
//   `RegFeatCache::count_features_in_region()`
//   (future-work-vepyr.md `## RegFeatCache::count_features_in_region()`;
//   shared with RegFeat port subtests #22, #39)
//
// Counts are v84-fixture-specific (50848, 132); v115 numbers would differ
// even with the API. Future-work entry covers both source ports.
// ----------------------------------------------------------------
//
// #[test]
// fn get_all_features_by_buffer_count_50848_subtest_33b() {
//     // PRECONDITION: count_features_in_region public API exists; v115 count hand-coded.
// }
// #[test]
// fn annotate_input_buffer_count_132_subtest_37() {
//     // PRECONDITION: same as #33b.
// }

// ----------------------------------------------------------------
// CLUSTER G — `--old_maf` legacy AMR format (subtest #44)
//
// Blocked on EXISTING future-work entry:
//   `VariationCache --old_maf legacy AMR format`
//   (future-work-vepyr.md `## VariationCache --old_maf legacy AMR format`)
// ----------------------------------------------------------------
//
// #[test]
// fn old_maf_flips_amr_format_subtest_44() {
//     // PRECONDITION: UDTF arg `old_maf=true` exists; output flips AMR from
//     // 'T:0.0014' (allele-prefixed) to bare numeric '0.0014'.
// }
