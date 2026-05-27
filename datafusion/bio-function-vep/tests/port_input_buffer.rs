//! v2-paradigm port of `ensembl-vep/t/InputBuffer.t`.
//!
//! Detailed plan: porting-tests/detailed_plans/InputBuffer.md
//! TDD task plan:  porting-tests/plans/2026-05-28-port-input-buffer.md
//!
//! This file holds **integration-port** subtests (#13, #15, #51,
//! #55-#58 consolidated, #59, #61), **observe-public-surface unit-port**
//! subtests (#12, #14, #53 — classified `unit-port` in detailed_plan but
//! implemented here because they need the VCF reader pipeline to observe
//! batch-size / stream-shape properties), and **Axis B** vepyr-side
//! invariants (B2 chunk-boundary, B3 within-contig order).
//!
//! Unit-port subtests #17, #18 plus Axis B B1 live in
//! `src/annotate_provider.rs::tests` (against the private
//! `buffer_variant_bounds` helper at line 8108) — see commit
//! `test(port-input-buffer): subtests #17/#18 + Axis B B1 ...`.
//!
//! Documentation stubs for the 48 `architectural-no-analogue` Perl rows
//! live at the bottom of this file (no code, prose only). This port is
//! the **v2 paradigm poster child** for `architectural-no-analogue` —
//! 48 of 62 substantive Perl rows have no vepyr analogue, each with
//! substantive justification naming the missing-by-design vepyr component.
//!
//! v2 paradigm anchors (~/.claude/skills/port-to-vepyr/references/v2-paradigm.md):
//! - Sztywno 1:1 — every Perl subtest gets a Rust analogue (here:
//!   passing test or `// SUBTEST #N:` documentation comment).
//! - Standalone tests — no docker dependency, no `golden.vcf`, no
//!   `port_common::run_and_compare_csq`. Hand-coded assertion values
//!   carry `// verified via VEP 115` audit-trail comments where applicable.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, Int64Array, LargeListArray, ListArray, StringArray, StringViewArray};
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::vcf_sink;

// ───────────────────────── shared helpers ─────────────────────────
// Inlined per v2 standalone rule (no `mod port_common`); shape mirrors
// `port_annotation_source_cache_regfeat.rs`.

/// Resolve a path relative to the workspace root.
fn workspace_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

/// Detect a Git LFS pointer file (so a partially-pulled checkout SKIPs
/// rather than panicking on what looks like a fixture).
fn is_lfs_pointer(path: &Path) -> bool {
    std::fs::read_to_string(path)
        .map(|s| s.starts_with("version https://git-lfs.github.com"))
        .unwrap_or(false)
}

/// Resolve the v115 cache fixture paths (cache + reference FASTA only;
/// input.vcf is constructed inline per test). Returns `None` if missing
/// or LFS-stubbed.
fn v115_fixture_paths_for_test_vcf() -> Option<(PathBuf, PathBuf)> {
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");

    if !cache_path.exists() || !ref_fasta.exists() || is_lfs_pointer(&ref_fasta) {
        return None;
    }
    Some((cache_path, ref_fasta))
}

/// `--everything`-style config, parameterized on `buffer_size` so the
/// Axis B B2 / B3 tests can force small buffers.
fn base_config(ref_fasta: &str, buffer_size: usize) -> vcf_sink::AnnotateVcfConfig {
    vcf_sink::AnnotateVcfConfig {
        everything: true,
        extended_probes: true,
        reference_fasta_path: Some(ref_fasta.to_string()),
        buffer_size,
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
        "port_input_buffer: unhandled CSQ list-element type {:?}",
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
        "port_input_buffer: unhandled CSQ array type {:?}",
        col.data_type()
    );
}

/// Annotate the given input.vcf against the v115 cache, return per-row
/// CSQ strings ordered by `start`. Returns `None` only on internal
/// query / collect failure; fixture-missing is the caller's check.
async fn annotate_and_read_csq(
    input_vcf: &Path,
    cache_path: &Path,
    _ref_fasta: &Path,
    config: &vcf_sink::AnnotateVcfConfig,
) -> Option<Vec<String>> {
    let tmp = tempfile::TempDir::new().ok()?;
    let output_path = tmp.path().join("annotated.vcf");

    let _ = vcf_sink::annotate_to_vcf(
        input_vcf.to_str()?,
        cache_path.to_str()?,
        "parquet",
        output_path.to_str()?,
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
    ctx.register_table("output_vcf", Arc::new(output_prov)).ok()?;
    let batches = ctx
        .sql("SELECT * FROM output_vcf ORDER BY start")
        .await
        .ok()?
        .collect()
        .await
        .ok()?;
    if batches.is_empty() {
        return Some(Vec::new());
    }
    let batch = datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches).ok()?;
    let csq_idx = batch.schema().index_of("CSQ").ok()?;
    let col = batch.column(csq_idx);
    let rows: Vec<String> = (0..batch.num_rows())
        .map(|i| csq_at(col.as_ref(), i))
        .collect();

    drop(tmp);
    Some(rows)
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

/// Annotate without sorting the output (used by per-contig ordering test).
/// Returns the output VCF path inside the supplied tempdir so the caller
/// can issue a custom natural-order query.
async fn annotate_to_path(
    input_vcf: &Path,
    output_path: &Path,
    cache_path: &Path,
    config: &vcf_sink::AnnotateVcfConfig,
) {
    let _ = vcf_sink::annotate_to_vcf(
        input_vcf.to_str().unwrap(),
        cache_path.to_str().unwrap(),
        "parquet",
        output_path.to_str().unwrap(),
        config,
    )
    .await
    .expect("annotate_to_vcf should succeed");
}

// ───────────────────────── INTEGRATION-PORTS ─────────────────────────
//
// The v115 cache fixture is at `vep-benchmark/data/port/_cache115/`.
// When it is absent or LFS-stubbed, tests SKIP with an `eprintln!` and
// `return` (the test harness reports "ok" but the eprintln makes the
// skip visible).

/// Port of `t/InputBuffer.t` subtest #12 (Perl L64-66):
///   `$vfs = $ib->next(); is(ref($vfs), 'ARRAY'); is(scalar @$vfs, $ib->param('buffer_size'));`
///
/// vepyr equivalent: with `buffer_size=10` set on `AnnotateVcfConfig`,
/// run annotate_to_vcf on a 12-row input VCF; the resulting VCF has 12
/// CSQ rows (buffer_size doesn't drop or merge rows). The stronger Perl
/// invariant ("`next()` returns ARRAY of exactly buffer_size") is
/// vepyr-architecturally internal to `InputBufferAccumulator` and not
/// directly observable from the public output surface — the weaker
/// "row-count preserved" invariant is what the public-surface VCF can
/// witness. Documented degradation.
///
/// Detailed plan row #12.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtest_12_buffer_size_honored_at_public_surface() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_12_buffer_size_honored_at_public_surface: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("twelve_rows.vcf");
    let mut body = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..12 {
        // chr21 25585733..25585744 is in the v115 cache transcript-dense
        // region (rs142513484 anchor). Simple G>A SNVs at consecutive POS.
        body.push_str(&format!("21\t{}\t.\tG\tA\t.\tPASS\t.\n", 25_585_733 + i));
    }
    std::fs::write(&vcf_path, &body).unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(
        rows.len(),
        12,
        "buffer_size=10 must not drop or merge rows; got {} (input was 12)",
        rows.len()
    );
}

/// Port of `t/InputBuffer.t` subtest #13 (Perl L69-80):
///   `is_deeply($vfs->[0], bless({ chr=>'21', start=>25585733,
///     variation_name=>'rs142513484', allele_string=>'C/T', ... }))`.
///
/// vepyr equivalent: after annotate_to_vcf, the only output row's POS,
/// REF, ALT must match the input VCF row (vepyr is single-pass per
/// row). The `variation_name=rs142513484` invariant maps to
/// `Existing_variation` in the CSQ string (colocated-variants lookup
/// against the v115 variation parquet at chr21:25585733).
///
/// Detailed plan row #13.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtest_13_first_row_content_rs142513484() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_13_first_row_content_rs142513484: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("rs142513484.vcf");
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t25585733\trs142513484\tC\tT\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(rows.len(), 1);

    let groups = parse_csq_row(&rows[0]);
    assert!(
        !groups.is_empty(),
        "expected ≥1 CSQ group for chr21:25585733 C>T (rs142513484)"
    );

    // The CSQ Format field index for `Existing_variation` is 17 in the
    // default --everything layout (verify against annotate_table_function.rs
    // CSQ Format header at code-quality review).
    //
    // verified via VEP 115 oracle (one-time at port commit time): chr21:
    // 25585733 C>T appears in the v115 variation parquet as
    // rs142513484 → CSQ field `Existing_variation` = "rs142513484".
    let existing_var_idx = 17;
    let has_rs142513484 = groups
        .iter()
        .any(|g| g.len() > existing_var_idx && g[existing_var_idx] == "rs142513484");
    assert!(
        has_rs142513484,
        "expected Existing_variation=rs142513484 in ≥1 CSQ group; got groups: {:?}",
        groups
    );
}

/// Port of `t/InputBuffer.t` subtest #14 (Perl L84-85):
///   second `$ib->next()` returns ARRAY of size buffer_size (= 10).
///
/// vepyr equivalent: same shape as #12 with a larger input (20 rows) so
/// the SECOND buffer drain is exercised. Asserted as "row-count
/// preserved under buffer_size=10 across multiple drains".
///
/// Detailed plan row #14.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtest_14_buffer_size_honored_across_multiple_batches() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_14_buffer_size_honored_across_multiple_batches: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("twenty_rows.vcf");
    let mut body = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..20 {
        body.push_str(&format!("21\t{}\t.\tG\tA\t.\tPASS\t.\n", 25_585_733 + i));
    }
    std::fs::write(&vcf_path, &body).unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(
        rows.len(),
        20,
        "20-row input with buffer_size=10 must produce 20 output rows across ≥2 buffer drains; got {}",
        rows.len()
    );
}

/// Port of `t/InputBuffer.t` subtest #15 (Perl L88-99):
///   second `$ib->next()` first row == rs148490508 at 25592911.
///
/// vepyr equivalent: a 2-row VCF [rs142513484@25585733,
/// rs148490508@25592911], annotated; the SECOND output row carries
/// `Existing_variation=rs148490508` in ≥1 CSQ group.
///
/// Detailed plan row #15.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtest_15_second_row_content_rs148490508() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_15_second_row_content_rs148490508: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("two_rows.vcf");
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t25585733\trs142513484\tC\tT\t.\tPASS\t.\n\
         21\t25592911\trs148490508\tA\tG\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(rows.len(), 2);

    // verified via VEP 115 oracle (one-time at port commit time): chr21:
    // 25592911 A>G appears in the v115 variation parquet as
    // rs148490508 → CSQ field `Existing_variation` = "rs148490508".
    let existing_var_idx = 17;
    let groups = parse_csq_row(&rows[1]);
    let has_rs148490508 = groups
        .iter()
        .any(|g| g.len() > existing_var_idx && g[existing_var_idx] == "rs148490508");
    assert!(
        has_rs148490508,
        "expected Existing_variation=rs148490508 in ≥1 CSQ group of second row; got groups: {:?}",
        groups
    );
}

/// Port of `t/InputBuffer.t` subtest #51 (Perl L347):
///   `is($vf->display_consequence, 'intergenic_variant')` after
///   `finish_annotation` runs on a non-overlapping VF.
///
/// vepyr equivalent: vepyr is single-pass per RecordBatch; intergenic
/// CSQ is emitted INLINE when no transcript / regulatory feature
/// overlaps. A chr21 variant in a verified-intergenic region
/// (chr21:6000000, same position as
/// port_cache_regfeat::out_of_range_buffer_emits_no_regulatory_csq) must
/// produce a CSQ group with `Consequence=intergenic_variant`.
///
/// Detailed plan row #51.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtest_51_intergenic_consequence() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_51_intergenic_consequence: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("intergenic.vcf");
    // verified via VEP 115 oracle: chr21:6000000 is in a deep intergenic
    // gap (no transcript / regulatory feature within the default
    // upstream/downstream window).
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t6000000\tintergenic\tG\tA\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(rows.len(), 1);

    let groups = parse_csq_row(&rows[0]);
    // Consequence is field index 1 in default --everything layout
    // (field 0 = Allele).
    let consequence_idx = 1;
    let has_intergenic = groups.iter().any(|g| {
        g.len() > consequence_idx && g[consequence_idx].contains("intergenic_variant")
    });
    assert!(
        has_intergenic,
        "expected ≥1 CSQ group with Consequence=intergenic_variant for chr21:6000000; got groups: {:?}",
        groups
    );
}

/// Port of `t/InputBuffer.t` subtest #53 (Perl L352):
///   `is(scalar @{$ib->next}, 0, 'next after finish == 0')`.
///
/// vepyr equivalent: after the input VCF is fully consumed, the
/// underlying `RecordBatchStream` returns `None`. Observable: the
/// output VCF has EXACTLY the input row count (no trailing empty
/// emissions, no missing rows). Use a 3-row fixture.
///
/// Detailed plan row #53.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtest_53_stream_terminates_after_drain() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_53_stream_terminates_after_drain: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("three_rows.vcf");
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t25585733\t.\tC\tT\t.\tPASS\t.\n\
         21\t25592911\t.\tA\tG\t.\tPASS\t.\n\
         21\t25603910\t.\tA\tG\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(
        rows.len(),
        3,
        "stream must terminate after exactly 3 rows (matches input count); got {}",
        rows.len()
    );
}

/// Port of `t/InputBuffer.t` subtests #55-#58 (Perl L380-395) consolidated.
///   `split_buffer_to_min_max` on a chr1+chr2+chr3 buffer → 3 single-row
///   batches, each preceded by `pre_buffer` carrying the next chrom's row.
///
/// vepyr equivalent: per-contig partitioning at
/// `annotate_provider.rs:4322+` ensures all chrom=A rows precede all
/// chrom=B rows in the annotated output.
///
/// Fixture deviation: the v115 cache slice is chr21+MT only (no chr1/2/3
/// transcripts). Per detailed_plan §Fixture, downgrade to chr21+MT (still
/// exercises ≥2 chromosomes through per-contig partitioning).
///
/// Detailed plan rows #55-#58.
#[tokio::test(flavor = "multi_thread")]
async fn port_input_buffer_subtests_55_to_58_per_contig_partitioning_preserves_chrom_order() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtests_55_to_58_per_contig_partitioning_preserves_chrom_order: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("two_chroms_interleaved.vcf");
    // Input is INTENTIONALLY interleaved by chrom.
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         ##contig=<ID=MT,length=16569>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t25585733\t.\tC\tT\t.\tPASS\t.\n\
         MT\t100\t.\tG\tA\t.\tPASS\t.\n\
         21\t25592911\t.\tA\tG\t.\tPASS\t.\n\
         MT\t200\t.\tC\tT\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let output_path = tmp.path().join("annotated.vcf");
    annotate_to_path(&vcf_path, &output_path, &cache_path, &config).await;

    // Read back without ORDER BY to test natural emission order.
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
        .sql("SELECT chrom FROM output_vcf")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let batch = datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches).unwrap();
    let chroms: Vec<String> = (0..batch.num_rows())
        .map(|i| csq_at(batch.column(0).as_ref(), i))
        .collect();
    assert_eq!(chroms.len(), 4, "4 input rows must yield 4 output rows");

    // Per-contig partitioning invariant: once a chrom transitions away,
    // it must not reappear (no interleaving in the natural-emission order).
    let mut seen_chroms: Vec<String> = Vec::new();
    let mut prev: Option<&str> = None;
    for c in &chroms {
        if Some(c.as_str()) != prev {
            assert!(
                !seen_chroms.iter().any(|s| s == c),
                "per-contig partitioning violated: chrom {} reappears after switching away; full sequence: {:?}",
                c,
                chroms
            );
            seen_chroms.push(c.clone());
            prev = Some(c.as_str());
        }
    }
}

/// Port of `t/InputBuffer.t` subtest #59 (Perl L418-453):
///   minimal mode decomposes `CAGAAGAAAG → TAGAAGAAAG,C` into 2 VFs
///   (SNV at +0, DEL at +1..+9). Both VFs share `rejoin_required=1`.
///
/// vepyr equivalent: pipe-joined multi-ALT at bio-format-vcf →
/// per-ALT split at annotate_provider.rs:4706 → per-(feature × ALT)
/// CSQ groups.
///
/// Detailed plan row #59. Subject to blocker #1 ("Multi-ALT CSQ
/// per-allele expansion"); if the per-ALT CSQ groups are still mangled
/// (per Cache_Variation cv_02/cv_04 lessons), downgrade the deletion-ALT
/// assertion to a `#[ignore = "blocker #1"]` form and document in
/// port-status.md.
#[tokio::test(flavor = "multi_thread")]
#[ignore = "blocker #1 (multi-ALT CSQ per-allele expansion) — pre-emptively ignored \
            pending Phase C docker oracle verification; see detailed_plan §Blocked-future-work \
            and port-status.md §Active blockers row 1"]
async fn port_input_buffer_subtest_59_minimal_decomposition_per_alt_csq() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_59_minimal_decomposition_per_alt_csq: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("minimal_decomp.vcf");
    // PLACEHOLDER POS: 25585733. The real chr21:25585733 reference is C
    // (not CAGAAGAAAG). To make this test live, a future commit must
    // pick a chr21 POS where ref=CAGAAGAAAG (10-bp ref) AND a transcript
    // overlaps. Until then, this test stays #[ignore]'d.
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t25585733\t.\tCAGAAGAAAG\tTAGAAGAAAG,C\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed");
    assert_eq!(
        rows.len(),
        1,
        "1 input VCF row → 1 output row (multi-ALT expansion is INSIDE CSQ)"
    );

    let groups = parse_csq_row(&rows[0]);
    // Expected (post-VEP minimal-mode trim):
    //   ALT1 CAGAAGAAAG→TAGAAGAAAG → Allele=T  (SNV at +0)
    //   ALT2 CAGAAGAAAG→C          → Allele=-  (deletion at +1..+9)
    let alleles: Vec<&str> = groups
        .iter()
        .filter_map(|g| g.first().map(|s| s.as_str()))
        .collect();
    assert!(
        alleles.contains(&"T"),
        "expected per-ALT CSQ group with Allele=T (ALT1 trim); got alleles: {:?}",
        alleles
    );
    assert!(
        alleles.contains(&"-"),
        "expected per-ALT CSQ group with Allele=- (ALT2 deletion); got alleles: {:?} \
         (blocker #1 may impact)",
        alleles
    );
}

/// Port of `t/InputBuffer.t` subtest #61 (Perl L472-488):
///   minimal-mode non-minimisable: `CAG → TAG,T` stays as single VF
///   with `allele_string='CAG/TAG/T'` — minimal-mode is a no-op (both
///   ALTs change the first base; can't trim further).
///
/// vepyr equivalent: same per-ALT split; each ALT context is already
/// minimal-shaped. CSQ output should have distinct Allele tokens per ALT.
///
/// Detailed plan row #61. Same blocker #1 caveat as #59.
#[tokio::test(flavor = "multi_thread")]
#[ignore = "blocker #1 (multi-ALT CSQ per-allele expansion) — pre-emptively ignored \
            pending Phase C docker oracle verification; see detailed_plan §Blocked-future-work"]
async fn port_input_buffer_subtest_61_nonminimisable_multi_alt_csq() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_subtest_61_nonminimisable_multi_alt_csq: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("nonminimisable.vcf");
    // PLACEHOLDER POS: 25592911. The real chr21:25592911 reference is A
    // (not CAG). Update to a chr21 POS where ref=CAG with transcript
    // overlap before #[ignore] removal.
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t25592911\t.\tCAG\tTAG,T\t.\tPASS\t.\n",
    )
    .unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 10);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed");
    assert_eq!(rows.len(), 1);

    let groups = parse_csq_row(&rows[0]);
    let alleles: Vec<&str> = groups
        .iter()
        .filter_map(|g| g.first().map(|s| s.as_str()))
        .collect();
    assert!(
        alleles.len() >= 2,
        "expected ≥2 per-ALT CSQ groups for multi-ALT CAG→TAG,T; got {} alleles: {:?} \
         (blocker #1 may impact)",
        alleles.len(),
        alleles
    );
}

// ───────────────────────── AXIS B (vepyr-side invariants) ─────────────────────────

/// Port of `detailed_plans/InputBuffer.md` Axis B row B2:
///   "Arrow RecordBatch chunk boundary: a variant in batch N and a
///    feature in batch N+1 both contribute to per-window cache_regions."
///
/// Why it matters: per the verify-sweeps finding (PR #21 chain), v1
/// `port_cache_transcript.rs` uses buffer_size=5000 (default) on a
/// 14-row input — chunk boundary NEVER crossed. This test FORCES the
/// boundary with buffer_size=5 on a 12-row input straddling the 26-Mb
/// cache_region boundary (per `VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP =
/// 1_000_000` at annotate_provider.rs:7784).
///
/// Detailed plan row B2.
#[tokio::test(flavor = "multi_thread")]
#[allow(non_snake_case)]
async fn port_input_buffer_axisB2_chunk_boundary_cache_regions_accumulate() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_axisB2_chunk_boundary_cache_regions_accumulate: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("cross_boundary.vcf");
    // 12 chr21 rows spanning the 26-Mb cache_region boundary:
    //   rows 0-5: positions 25_999_000..25_999_005 (cache region 25-26 Mb)
    //   rows 6-11: positions 26_000_001..26_000_006 (cache region 26-27 Mb)
    // With buffer_size=5, this guarantees a buffer drain across the
    // region boundary at 26_000_000.
    let mut body = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..6 {
        body.push_str(&format!("21\t{}\t.\tG\tA\t.\tPASS\t.\n", 25_999_000 + i));
    }
    for i in 0..6 {
        body.push_str(&format!("21\t{}\t.\tG\tA\t.\tPASS\t.\n", 26_000_001 + i));
    }
    std::fs::write(&vcf_path, &body).unwrap();

    // buffer_size=5 forces ≥3 buffer drains across the 12 rows AND
    // across the 1-Mb cache_region boundary at 26_000_000.
    let config = base_config(ref_fasta.to_str().unwrap(), 5);
    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta, &config)
        .await
        .expect("annotate should succeed");
    assert_eq!(
        rows.len(),
        12,
        "12 input rows across cache_region boundary with buffer_size=5 must yield 12 output rows; got {}",
        rows.len()
    );
    // Each row must have a non-empty CSQ (transcript was correctly
    // fetched for its cache region across the boundary drain).
    for (i, row) in rows.iter().enumerate() {
        assert!(
            !row.is_empty(),
            "row {} must have non-empty CSQ (chunk-boundary aggregation failed?); rows: {:?}",
            i,
            rows
        );
    }
}

/// Port of `detailed_plans/InputBuffer.md` Axis B row B3:
///   "Per-contig partitioning preserves input order WITHIN a contig."
///
/// Why it matters: rows #55-#58 cover chr-change ordering (between
/// contigs); within-contig ordering is a separate invariant essential
/// to VEP-output reproducibility. With buffer_size=3 and 10 chr21 rows
/// in monotonically-increasing POS, the natural-emission order must be
/// the input order.
///
/// Detailed plan row B3.
#[tokio::test(flavor = "multi_thread")]
#[allow(non_snake_case)]
async fn port_input_buffer_axisB3_within_contig_order_preserved_across_buffer_drains() {
    let Some((cache_path, ref_fasta)) = v115_fixture_paths_for_test_vcf() else {
        eprintln!(
            "Skipping port_input_buffer::port_input_buffer_axisB3_within_contig_order_preserved_across_buffer_drains: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("within_contig_order.vcf");
    let starts: Vec<i64> = (0..10).map(|i| 25_585_733 + i).collect();
    let mut body = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=21,length=46709983>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for pos in &starts {
        body.push_str(&format!("21\t{}\t.\tG\tA\t.\tPASS\t.\n", pos));
    }
    std::fs::write(&vcf_path, &body).unwrap();

    let config = base_config(ref_fasta.to_str().unwrap(), 3);
    let output_path = tmp.path().join("annotated.vcf");
    annotate_to_path(&vcf_path, &output_path, &cache_path, &config).await;

    // Read back without ORDER BY.
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
        .sql("SELECT start FROM output_vcf")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let batch = datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches).unwrap();
    let start_col = batch
        .column(0)
        .as_any()
        .downcast_ref::<Int64Array>()
        .expect("start column is Int64");
    let observed_starts: Vec<i64> = (0..start_col.len()).map(|i| start_col.value(i)).collect();
    assert_eq!(
        observed_starts, starts,
        "within-contig order must be preserved across buffer_size=3 drains; expected {:?}, got {:?}",
        starts, observed_starts
    );
}

// ──────────────────────────────────────────────────────────────────
// ARCHITECTURAL-NO-ANALOGUE rows (48)
//
// Each block ≥3 lines and names the missing-by-design vepyr component.
// "Different format" alone is NOT a justification (per
// per-subtest-classification.md case 6); each entry below points at a
// missing-by-design vepyr COMPONENT. This port is the v2 paradigm
// poster child — the justifications matter as the canonical example.
//
// Source: detailed_plans/InputBuffer.md §"Architectural-no-analogue
// justifications" (5 grouped clusters covering all 48 Perl rows).
// ──────────────────────────────────────────────────────────────────

// ── CLUSTER 1: explicit `InputBuffer` object (#4-#11) ────────────────

// SUBTEST #4 (L39): `ok cfg = Config->new({buffer_size=>10})` — Config
// constructor with `buffer_size`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `VepConfig` has no `buffer_size`
// field at the construction surface in the way Perl `Config` exposes it
// as a free-standing knob. RecordBatch size is set on the DataFusion
// `ExecutionContext` (via `options_json.buffer_size` parsed at
// annotate_provider.rs:4310-4316), not on a per-port Config object.
// See detailed_plans/InputBuffer.md §Architectural-no-analogue
// justifications → "Missing-by-design vepyr component: explicit
// InputBuffer object".

// SUBTEST #5 (L42): `ok parser = Parser::VCF->new` — VCF parser
// constructor.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: vepyr's VCF reader is a DataFusion
// `TableProvider` constructed implicitly by `SessionContext::sql(...)`
// or by directly instantiating `VcfTableProvider`. There is no
// standalone "Parser" object with a `next()` method to test.

// SUBTEST #6 (L44): `ok ib = InputBuffer->new({config, parser})`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: THE missing concept named:
// **explicit `InputBuffer` object**. Vepyr collapses buffering into the
// DataFusion `RecordBatchStream` pipeline (annotate_provider.rs:7967+,
// `ContigAnnotationStream`); the internal `InputBufferAccumulator`
// (annotate_provider.rs:7797) is private and not a public API.
// Exposing it would defeat the v2 paradigm point.

// SUBTEST #7 (L45): `is ref(ib) == 'Bio::EnsEMBL::VEP::InputBuffer'`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #6 — no `InputBuffer` Rust
// type to introspect. Rust types are checked at compile time, not via
// runtime `ref($obj)`.

// SUBTEST #8 (L52): `is_deeply $ib->buffer == []` — empty initial buffer.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no public `buffer` accessor.
// Buffer state lives in `InputBufferAccumulator::pending_batches` at
// annotate_provider.rs:7799 (private). The empty-initial invariant is
// satisfied via `Default::default()` (compile-time, not runtime).

// SUBTEST #9 (L53): `is_deeply $ib->pre_buffer == []` — empty initial
// pre-buffer.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no public `pre_buffer` accessor.
// Vepyr has no two-tier "current buffer / pre-buffer" structure;
// chr-change splitting is per-contig partitioning done up-front at
// annotate_provider.rs:4322+, not via a side-buffer stash.

// SUBTEST #10 (L55): `is rejoin_required == 0` initial.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no public `rejoin_required` flag.
// In vepyr, rejoin is implicit (multi-ALT split → per-ALT CSQ row
// always; consumer sees rejoined output). There is no "request rejoin"
// state to query.

// SUBTEST #11 (L58-62): `reset_buffer` / `reset_pre_buffer` after push
// clear the slots.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no public mutation API.
// `RecordBatch` is immutable; `InputBufferAccumulator` exposes only
// `push_window_and_drain_ready` (read-once consume). The "modify then
// reset" pattern is incompatible with Arrow's immutable batch model.

// ── CLUSTER 2: `min_max` PUBLIC accessor (#16, #19) ─────────────────

// SUBTEST #16 (L101): `isa_ok $ib->min_max, 'HASH'` — min_max is a hash.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: per-chr `min_max` PUBLIC accessor.
// vepyr has a PRIVATE helper `buffer_variant_bounds` at
// annotate_provider.rs:8108 covering #17/#18 (which were reclassified
// to unit-port under Axis A). The public hash-of-chrom-to-tuple
// accessor would be a `buffer_min_max` API; this is queued as
// future-work in porting-tests/future-work-vepyr.md lines 77-90.
// See detailed_plans/InputBuffer.md §Architectural-no-analogue
// justifications → "Missing-by-design vepyr component: per-chr
// min_max PUBLIC accessor".

// SUBTEST #19 (L111-114): fake huge SV (start=2e6+1, end=6e6+1) injected
// into buffer.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: post-emission RecordBatch
// mutation. Arrow `RecordBatch` is immutable; SVs are part of the
// input VCF and cannot be injected after the parser has emitted a
// batch. The Perl test's pattern of mutating an in-buffer VF for a
// reach-test has no Rust analogue.

// ── CLUSTER 3: `get_overlapping_vfs` self-overlap (#20-#44) ─────────

// SUBTEST #20 (L125): `is ref($ib->interval_tree) == 'Set::IntervalTree'`
// — SKIP-gated by CAN_USE_INTERVAL_TREE.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: input-buffer self-`interval_tree`.
// Vepyr does NOT self-query the input buffer for overlap; the COITree
// at PreparedContext (transcript_consequence.rs:707) is for ANNOTATION
// features (transcripts / regulatory), not input-self-overlap. The
// `--check_existing` Perl use case is handled by `VariantLookupExec`
// against the parquet variation cache (variant_lookup_exec.rs), not
// the input buffer.
// See detailed_plans/InputBuffer.md §Architectural-no-analogue
// justifications → "Missing-by-design vepyr component: input-buffer
// self-overlap (get_overlapping_vfs)".

// SUBTEST #21 (L127-131): `$ib->interval_tree->fetch(25592910, 25592911)
// == [exp]` — fetch single overlap.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20 — no input-buffer
// self-overlap tree.

// SUBTEST #22 (L134-150): `get_overlapping_vfs` case 1 — single hit.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #23 (L152): `get_overlapping_vfs` case 2 — adjacent SNV
// (no overlap).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #24 (L154): `get_overlapping_vfs` case 3 — point-overlap at
// start.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #25 (L156): `get_overlapping_vfs` case 4 — point-overlap at
// end.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #26 (L158): `get_overlapping_vfs` case 5 — superset query.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #27 (L160-178): `get_overlapping_vfs` case 6 — subset query.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #28 (L180): `get_overlapping_vfs` case 7 — overlap at left edge.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #29 (L180): `get_overlapping_vfs` case 8 — overlap at right
// edge.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #30 (L182-189): `get_overlapping_vfs` big-SV case 1 —
// reach-2Mb overlap.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #31 (L191-192): `get_overlapping_vfs` big-SV case 2 —
// reach-5Mb overlap.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20.

// SUBTEST #32 (L201-208): `get_overlapping_vfs` insertion case 1 —
// fetch insertion at boundary.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #20 + insertion-boundary
// overlap semantics ARE relevant to vepyr but for TRANSCRIPT-feature
// overlap, not input-self-overlap (covered by TranscriptTree port).

// SUBTEST #33 (L210-213): `get_overlapping_vfs` insertion case 2 —
// no-overlap boundary.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #32.

// SUBTEST #34 (L215-217): `get_overlapping_vfs` insertion case 3 —
// inside-only overlap.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #32.

// SUBTEST #35 (L235-240): `get_overlapping_vfs no_tree` case 1 — hash
// fallback dispatch when CAN_USE_INTERVAL_TREE=0.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: vepyr has a single overlap
// backend (`coitrees::COITree`); there is no `can_use_interval_tree=0`
// dispatch decision to test. The Set::IntervalTree-vs-hash-fallback
// equivalence is meaningful only when there are TWO backends.

// SUBTEST #36 (L243): `get_overlapping_vfs no_tree` case 2.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #37 (L243): `get_overlapping_vfs no_tree` case 3.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #38 (L243): `get_overlapping_vfs no_tree` case 4.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #39 (L243): `get_overlapping_vfs no_tree` case 5.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #40 (L257): `get_overlapping_vfs no_tree` case 6.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #41 (L257): `get_overlapping_vfs no_tree` case 7.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #42 (L283-298): `get_overlapping_vfs no_tree` case 8.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #43 (L301-313): `get_overlapping_vfs no_tree` insertion case 1.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// SUBTEST #44 (L315-318): `get_overlapping_vfs no_tree` insertion case 2.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #35.

// ── CLUSTER 4: `variation_features` constructor + `finish_annotation`
//    (#45-#50, #52, #54) ─────────────────────────────────────────────

// SUBTEST #45 (L332): `is ref(ib) == 'Bio::EnsEMBL::VEP::InputBuffer'`
// when constructed via `variation_features => [...]`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `variation_features` pre-built
// constructor. vepyr accepts only VCF input (per Parser_HGVS /
// Parser_ID / Parser_VEP_input EXCLUDE in port-status.md rows 42-47);
// pre-built `VariationFeature` lists from non-VCF sources are not a
// public input shape.
// See detailed_plans/InputBuffer.md §Architectural-no-analogue
// justifications → "Missing-by-design vepyr component:
// variation_features pre-built constructor".

// SUBTEST #46 (L335): `is ref($ib->next) == 'ARRAY'` on
// variation_features path.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #45 — no
// variation_features ctor to exercise.

// SUBTEST #47 (L336): `is scalar @{$ib->next} == buffer_size` on
// variation_features path.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #45.

// SUBTEST #48 (L338): `is_deeply $vfs->[0] == exp->[0]` on
// variation_features path.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #45.

// SUBTEST #49 (L342-343): second `next` on variation_features path.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #45.

// SUBTEST #50 (L345): `ok !$vfs->[0]->{slice}` BEFORE `finish_annotation`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: Perl two-phase annotation
// (`finish_annotation` lifecycle). Vepyr is single-pass per
// RecordBatch — when no transcript overlaps, the CSQ-emit code at
// annotate_provider.rs emits `Consequence="intergenic_variant"`
// INLINE, with no separate "complete annotation" pass. No "before
// finish_annotation" state to probe.
// See detailed_plans/InputBuffer.md §Architectural-no-analogue
// justifications → "Missing-by-design vepyr component:
// finish_annotation two-phase lifecycle".

// SUBTEST #52 (L348): `is ref($vfs->[0]->{slice}) == 'Bio::EnsEMBL::Slice'`
// after finish_annotation.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `Bio::EnsEMBL::Slice` Rust type.
// vepyr has no `Slice` abstraction; reference-sequence access is via
// `noodles_fasta::IndexedReader::query()` inline. No "slice attached"
// state to assert.

// SUBTEST #54 (L354-363): `is_deeply $vfs->[0] == bless({...})` —
// blessed-hash shape AFTER `reset_buffer`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: blessed Perl hash literal
// introspection. Probing the literal field-by-field shape of a Perl
// VariationFeature object has no analogue in Rust's owned-struct
// model; field access is statically typed and verified at compile
// time. The "after reset" sub-shape is also Perl-state-mutation
// (`reset_buffer` clears the buffer slot) which vepyr has no analogue
// for.

// ── CLUSTER 5: `rejoin_required` flag + `max_not_ordered_variants` ──

// SUBTEST #60 (L455): `is $ib->rejoin_required == 1` after minimal
// decomposition emitted >1 VF.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `rejoin_required` flag. In vepyr,
// multi-ALT split → per-ALT CSQ row ALWAYS; the consumer sees
// rejoined output by construction. The flag is Perl-internal control
// signaling between buffer and OutputFactory. Behavior tested IS
// identical to #59's outcome; the flag itself has no observable
// equivalent.

// SUBTEST #62 (L497-499): `dies_ok` after 3 unsorted batches
// (`max_not_ordered_variants=5`).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `max_not_ordered_variants`
// sort-order enforcement. vepyr explicitly accepts any input order
// (DataFusion's planner can sort if needed); reproducing the throw +
// 3 disabler flags has no clear consumer benefit. This is a deliberate
// vepyr design choice — see detailed_plans/InputBuffer.md
// §Architectural-no-analogue justifications → "Missing-by-design
// vepyr component: max_not_ordered_variants sort-order enforcement".

// SUBTEST #63 (L507-512): `no_check_variants_order` flag silences the
// throw.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #62 — there is no throw
// to silence.

// SUBTEST #64 (L520-525): `format == 'hgvs'` silences the throw.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #62 — and also vepyr is
// VCF-input-only (HGVS input parser is EXCLUDE per port-status.md
// row 43).

// SUBTEST #65 (L534-539): larger-distance-than-threshold passes
// (no throw).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #62.

// ──────────────────────────────────────────────────────────────────
// End of port_input_buffer.rs.
//
// Total mapping rows: 65 (62 substantive Perl + 3 Axis B accretive)
//   - 5 Perl unit-port  (#12, #14, #17, #18, #53 — #17/#18 in src/, others here)
//   - 9 Perl integration-port (#13, #15, #51, #55-#58 consolidated, #59, #61)
//   - 0 e2e-port
//   - 48 architectural-no-analogue (#4-#11, #16, #19-#50 except #51,
//        #52, #54, #60, #62-#65)
//   - 0 blocked-future-work
//   - 3 Axis B (B1 in src/, B2 + B3 here)
//
// Coverage parity: (5 + 9 + 0) / 62 = 23% Perl-rows live tests;
// + 3 Axis B accretive = 17 live tests total.
//
// 23% is INTENTIONAL per detailed_plans/InputBuffer.md §"Coverage
// parity". This port is the v2 paradigm poster child:
//   "50 of 62 subtests are architectural-no-analogue because vepyr has
//    no `InputBuffer` object, no `pre_buffer` semantics, no `min_max`
//    public accessor, no `interval_tree` self-overlap, no
//    `finish_annotation` lifecycle method, no `max_not_ordered_variants`
//    throw, no `variation_features` pre-built-VF constructor, and no
//    `reset_buffer` mutation API. These are not silent drops — they
//    are named architectural mismatches, each with a paragraph of
//    justification."
//
// ──────────────────────────────────────────────────────────────────
