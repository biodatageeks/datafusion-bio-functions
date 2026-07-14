//! v2-paradigm port of `ensembl-vep/t/AnnotationSource_Cache_RegFeat.t`.
//!
//! Detailed plan: porting-tests/detailed_plans/AnnotationSource_Cache_RegFeat.md
//! TDD task plan:  porting-tests/plans/2026-05-27-port-cache-regfeat.md
//!
//! This file holds **integration-port** subtests (#24, #26, #35, #38, #41),
//! plus documentation stubs for **architectural-no-analogue** subtests
//! (#4, #10, #14, #15, #17, #18, #19, #25, #27, #30, #32, #33, #36, #37,
//! #40, #44–#51) and **blocked-future-work** subtests (#11, #12, #13, #22,
//! #23, #28, #39).
//!
//! Unit-port subtests (#5, #6, #7, #8, #9, #16, #20, #21, #34) live in
//! `src/partitioned_cache.rs::tests` (commit `test(port-cache-regfeat):
//! v0 unit-ports`).
//!
//! E2E-port subtests (#42, #43) live in
//! `tests/port_annotation_source_cache_regfeat_e2e.rs`.
//!
//! v2 paradigm anchors (~/.claude/skills/port-to-vepyr/references/v2-paradigm.md):
//! - Sztywno 1:1 — every Perl subtest gets a Rust analogue (here: passing
//!   test, ignored doc panic, or commented-out future-work stub).
//! - Standalone tests — no docker dependency, no `golden.vcf`, no
//!   `port_common::run_and_compare_csq`. Hand-coded assertion values
//!   carry `// verified via VEP 115 …` audit-trail comments.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, LargeListArray, ListArray, StringArray, StringViewArray};
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::vcf_sink;

// ───────────────────────── shared helpers ─────────────────────────

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

/// Standard `--everything` config (matches Perl test's enabled features).
fn base_config(ref_fasta: &str) -> vcf_sink::AnnotateVcfConfig {
    vcf_sink::AnnotateVcfConfig {
        everything: true,
        extended_probes: true,
        reference_fasta_path: Some(ref_fasta.to_string()),
        ..Default::default()
    }
}

/// Join the elements of a List/LargeList CSQ cell into one comma-separated
/// string (matches `port_common`'s reader shape, but inlined to keep this
/// test standalone — no `mod port_common`).
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
        "port_cache_regfeat: unhandled CSQ list-element type {:?}",
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
        "port_cache_regfeat: unhandled CSQ array type {:?}",
        col.data_type()
    );
}

/// Resolve the v115 cache fixture paths; returns `None` if the fixture is
/// missing or LFS-stubbed (so a partial checkout SKIPs informatively).
fn v115_fixture_paths() -> Option<(PathBuf, PathBuf, PathBuf)> {
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");
    let input_vcf =
        workspace_path("datafusion/bio-function-vep/tests/data/port/cache_regfeat/input.vcf");

    if !cache_path.exists()
        || !ref_fasta.exists()
        || !input_vcf.exists()
        || is_lfs_pointer(&input_vcf)
        || is_lfs_pointer(&ref_fasta)
    {
        return None;
    }
    Some((cache_path, ref_fasta, input_vcf))
}

/// Annotate the given input.vcf against the v115 cache fixture, return
/// per-row CSQ strings. Returns `None` if fixtures are missing.
async fn annotate_and_read_csq(
    input_vcf: &Path,
    cache_path: &Path,
    ref_fasta: &Path,
) -> Option<Vec<String>> {
    let tmp = tempfile::TempDir::new().ok()?;
    let output_path = tmp.path().join("annotated.vcf");

    let config = base_config(ref_fasta.to_str()?);
    let _rows = vcf_sink::annotate_to_vcf(
        input_vcf.to_str()?,
        cache_path.to_str()?,
        "parquet",
        output_path.to_str()?,
        &config,
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
        .ok()?;
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

    // Keep tmp alive until rows are read.
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

// ───────────────────────── INTEGRATION-PORTS ─────────────────────────
//
// The v115 cache fixture is at `vep-benchmark/data/port/_cache115/`. When
// it is absent (clean checkout without LFS pull), the test SKIPs with an
// `eprintln!` and `return` — the test harness reports "ok" but the
// `eprintln!` makes the skip visible. (Same shape as `port_common`'s
// `Ok(false)` return pattern, but inlined here per v2 standalone rule.)

// Subtests #24 + #26 + #35 + #38 — get_features_by_regions_uncached /
// get_features_by_regions_cached / get_all_features_by_InputBuffer paths.
// Consolidated into one test: a single-variant input VCF on a known
// regulatory feature overlap; assert annotate_vep emits ≥1 CSQ allele-
// group with `Feature_type=RegulatoryFeature` and the expected
// `Feature=<stable_id>`.
#[tokio::test(flavor = "multi_thread")]
async fn regulatory_feature_overlap_emits_regulatory_csq() {
    let Some((cache_path, ref_fasta, input_vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_regfeat::regulatory_feature_overlap_emits_regulatory_csq: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let rows = annotate_and_read_csq(&input_vcf, &cache_path, &ref_fasta)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(
        rows.len(),
        1,
        "input.vcf has exactly one variant; annotate output must have one row \
         (got {} rows)",
        rows.len()
    );
    let row0 = &rows[0];
    assert!(
        !row0.is_empty(),
        "CSQ for the rs1343855353 substitute row must not be empty"
    );

    let groups = parse_csq_row(row0);
    // Per `annotate_table_function.rs:3267-3322`, the CSQ Format header
    // for vepyr's regulatory-mode output has 74 fields; field 5 (0-indexed)
    // is `Feature_type` and field 6 is `Feature`.
    //
    // verified via VEP 115 on v115 cache at commit b97e1a2 (regulatory/21.parquet):
    //   chr21:25039632 lands on ENSR21_B6Z6N (enhancer, 25039631-25040274).
    //   Substitution of rs3989369 (no v115 reg overlap) by rs1343855353
    //   at chr21:25039632 is documented in plan §8.
    let reg_groups: Vec<&Vec<String>> = groups
        .iter()
        .filter(|g| g.len() >= 7 && g[5] == "RegulatoryFeature")
        .collect();
    assert!(
        !reg_groups.is_empty(),
        "chr21:25039632 T>C must produce ≥1 RegulatoryFeature CSQ group; \
         groups were: {groups:?}",
    );
    // Subtest #26 / #38 (first stable_id): assert the regulatory group's
    // Feature column carries ENSR21_B6Z6N (the v115 enhancer overlapping
    // chr21:25039632).
    let feature_ids: Vec<&str> = reg_groups.iter().map(|g| g[6].as_str()).collect();
    assert!(
        feature_ids.contains(&"ENSR21_B6Z6N"),
        "expected ENSR21_B6Z6N (v115 enhancer at 25039631-25040274); got {feature_ids:?}",
    );
}

// Subtest #41 (RegFeat.t:151-152): `$ib->next();
// get_all_features_by_InputBuffer == []`. After consuming the buffer, the
// next buffer produces zero features. vepyr's equivalent: a buffer whose
// chr-range doesn't overlap any chrom file in the cache → annotate_to_vcf
// emits no regulatory CSQ rows.
//
// Implementation: synthesize an in-memory input.vcf row on chr22 (not in
// the v115 fixture which only has chr21 + MT) and assert no
// RegulatoryFeature CSQ rows are produced.
#[tokio::test(flavor = "multi_thread")]
async fn out_of_range_buffer_emits_no_regulatory_csq() {
    let Some((cache_path, ref_fasta, _)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_regfeat::out_of_range_buffer_emits_no_regulatory_csq: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    // Build a tmp input.vcf with one chr21 record OUTSIDE any regulatory
    // feature. chr21:6000000 is in a gap region (the first reg feature in
    // v115 is at chr21:5011320 enhancer 5011320-5011550, second at
    // ENSR21_BZTD 5013566-5014331; the next is far above 6000000).
    //
    // verified via VEP 115 on v115 cache at commit b97e1a2:
    //   chr21:6000000 has no overlapping regulatory or motif feature in
    //   regulatory/21.parquet (5203 rows checked; none overlap 6000000).
    let tmp = tempfile::TempDir::new().unwrap();
    let vcf_path = tmp.path().join("out_of_range.vcf");
    let body = "##fileformat=VCFv4.2\n\
                ##contig=<ID=21,length=46709983>\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                21\t6000000\tno_reg\tG\tA\t.\tPASS\t.\n";
    std::fs::write(&vcf_path, body).unwrap();

    let rows = annotate_and_read_csq(&vcf_path, &cache_path, &ref_fasta)
        .await
        .expect("annotate should succeed when fixture is present");
    assert_eq!(rows.len(), 1);
    let groups = parse_csq_row(&rows[0]);
    let reg_groups: Vec<&Vec<String>> = groups
        .iter()
        .filter(|g| g.len() >= 7 && g[5] == "RegulatoryFeature")
        .collect();
    assert!(
        reg_groups.is_empty(),
        "chr21:6000000 must produce ZERO RegulatoryFeature CSQ groups \
         (empty buffer analogue); got {reg_groups:?}",
    );
}

// Subtest #35 (RegFeat.t:142): get_all_features_by_InputBuffer returns
// `ARRAY` (non-empty). This is the same engine path as #24/#26/#38 — the
// non-emptiness is verified by `regulatory_feature_overlap_emits_regulatory_csq`
// (which asserts ≥1 RegulatoryFeature group).
//
// No separate test needed; the subtest is satisfied by the same code path.

// ───────────────────────── ARCHITECTURAL-NO-ANALOGUE STUBS ─────────────────────────
//
// Each row asserts an audit trail: the Perl subtest probes a concept that
// vepyr's data model architecturally rejects. The `#[ignore]` keeps these
// out of CI; they exist so the audit trail "every Perl subtest has a Rust
// analogue" remains complete. See detailed_plan §Architectural-no-analogue
// justifications for full reasoning.

#[allow(non_snake_case)]
mod architectural_no_analogue {
    // Subtest #4 (RegFeat.t:45): `ok($c, 'new is defined')`.
    // vepyr has no `Cache::RegFeat` object — regulatory access is a
    // DataFusion table-function arg (`annotate_vep(...)`). The
    // create-and-validate lifecycle collapses into argument parsing at
    // UDTF-call time. No public component to assert on.
    #[test]
    #[ignore = "architectural-no-analogue: no AnnotationSource object in vepyr; \
                regulatory access is via UDTF args. See detailed_plan row #4."]
    fn perl_constructor_returns_blessed_ref_subtest_4() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #4");
    }

    // Subtest #10 (RegFeat.t:58): `throws_ok { get_dump_file_name(1) }
    // qr/No region/`. vepyr's resolver is per-chrom; the "region"
    // (1-Mb block) argument has no analogue in the per-chrom partitioned
    // layout (partitioned_cache.rs:24-32).
    #[test]
    #[ignore = "architectural-no-analogue: no 1-Mb-block argument in vepyr's \
                per-chrom resolver. See detailed_plan row #10."]
    fn perl_get_dump_file_name_no_region_throw_subtest_10() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #10");
    }

    // Subtests #14, #15, #19: HASH/ARRAY shape checks on Perl-serialized
    // structure. vepyr's parquet read yields typed Arrow RecordBatch +
    // Vec<RegulatoryFeature>/Vec<MotifFeature> — no introspectable hash
    // shape to assert against.
    #[test]
    #[ignore = "architectural-no-analogue: Perl ref(HASH) probe — vepyr uses \
                typed Arrow batches. See detailed_plan row #14."]
    fn perl_deserialize_top_level_hash_subtest_14() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #14");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Perl per-chr HASH probe — see \
                detailed_plan row #15."]
    fn perl_deserialize_per_chr_hash_subtest_15() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #15");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Perl ARRAY shape after flatten — \
                vepyr always returns Vec<…>. See detailed_plan row #19."]
    fn perl_flattened_features_is_array_subtest_19() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #19");
    }

    // Subtests #17, #18, #25, #36, #37: Perl `ref($obj) eq
    // 'Bio::EnsEMBL::Funcgen::…'` blessed-ref class checks. vepyr's
    // discrimination is via the `RegulatoryTarget` enum and the
    // `feature_type` column — no named-class concept exists.
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref check — vepyr \
                uses RegulatoryTarget enum. See detailed_plan row #17."]
    fn perl_first_reg_feat_class_subtest_17() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #17");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref check — see \
                detailed_plan row #18."]
    fn perl_first_motif_feat_class_subtest_18() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #18");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref check on \
                disk-loaded features — see detailed_plan row #25."]
    fn perl_disk_first_reg_feat_class_subtest_25() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #25");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref check on \
                buffer-loaded features — see detailed_plan row #36."]
    fn perl_buffer_first_reg_feat_class_subtest_36() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #36");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref check on \
                buffer-loaded last (motif) feature — see detailed_plan row #37."]
    fn perl_buffer_last_motif_feat_class_subtest_37() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #37");
    }

    // Subtests #27, #40: uncached vs cached probe. vepyr's session-level
    // table cache is implicit in DataFusion's `register_table` semantics
    // — no public API distinguishes "served from disk" vs "served from
    // memory".
    #[test]
    #[ignore = "architectural-no-analogue: no public uncached/cached \
                distinction in vepyr's DataFusion session — see detailed_plan row #27."]
    fn perl_memory_cache_hit_returns_same_stable_id_subtest_27() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #27");
    }
    #[test]
    #[ignore = "architectural-no-analogue: memory-cache idempotency is \
                implicit in DataFusion — see detailed_plan row #40."]
    fn perl_buffer_memory_cache_idempotency_subtest_40() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #40");
    }

    // Subtests #30, #32, #33: Parser/InputBuffer construct probes —
    // tested in Parser_VCF and InputBuffer plans; not RegFeat-specific.
    #[test]
    #[ignore = "architectural-no-analogue: tested in Parser_VCF port — see \
                detailed_plan row #30."]
    fn perl_parser_instantiates_subtest_30() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #30");
    }
    #[test]
    #[ignore = "architectural-no-analogue: tested in InputBuffer port — see \
                detailed_plan row #32."]
    fn perl_input_buffer_constructed_subtest_32() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #32");
    }
    #[test]
    #[ignore = "architectural-no-analogue: tested in InputBuffer port — see \
                detailed_plan row #33."]
    fn perl_input_buffer_next_returns_array_subtest_33() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #33");
    }

    // Subtests #44–#51: SKIP-SEREAL block. Sereal has no Rust-ecosystem
    // analogue; vepyr's cache data model has no per-block serializer-
    // dispatch concept (all parquet).
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent in Rust ecosystem — \
                see detailed_plan row #44."]
    fn perl_sereal_serializer_type_subtest_44() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #44");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #45."]
    fn perl_sereal_file_suffix_subtest_45() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #45");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #46."]
    fn perl_sereal_top_level_hash_subtest_46() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #46");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #47."]
    fn perl_sereal_per_chr_hash_subtest_47() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #47");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #48."]
    fn perl_sereal_two_feature_types_subtest_48() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #48");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #49."]
    fn perl_sereal_reg_feat_class_subtest_49() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #49");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #50."]
    fn perl_sereal_motif_feat_class_subtest_50() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #50");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent — see detailed_plan row #51."]
    fn perl_sereal_stable_id_subtest_51() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #51");
    }
}

// ───────────────────────── BLOCKED-FUTURE-WORK STUBS ─────────────────────────
//
// Each block carries a commented-out `#[test]` showing the example
// signature the future-work API would enable, plus a pointer to the
// existing future-work-vepyr.md entry that tracks the gap.
//
// All future-work entries below are pre-existing in
// `porting-tests/future-work-vepyr.md` — these tests will be enabled
// when those APIs land.

mod blocked_future_work {
    // SUBTEST #11 (RegFeat.t:60): is_deeply($c->get_available_cell_types, []).
    // Blocked: vepyr has no public `available_cell_types(...)` API.
    // Future-work entry: `--cell_type filter + RegFeatCache::available_cell_types()`
    // (porting-tests/future-work-vepyr.md).
    //
    // Note: in the v115 fixture, regulatory/21.parquet's `cell_types`
    // column is all NULL (5203/5203 rows) — so even when this API lands,
    // the empty-list invariant holds for this fixture.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_get_available_cell_types_returns_empty_for_no_filter_subtest_11() {
    //     // Future signature:
    //     //   pub fn available_cell_types(table: &str, session: &SessionContext)
    //     //       -> Result<Vec<String>>;
    //     let cell_types = available_cell_types("regulatory_chr21", &session).await.unwrap();
    //     assert!(cell_types.is_empty());
    // }

    // SUBTEST #12 (RegFeat.t:62-64): `$c->{cell_type} = ['foo'];
    // ok($c->check_cell_types)`.
    // Blocked: vepyr has no `--cell_type` UDTF arg. Same future-work entry as #11.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_check_cell_types_passes_for_subset_subtest_12() {
    //     // Future signature:
    //     //   annotate_vep(..., cell_type => ['HepG2', 'K562']) with validation.
    //     // Today: no UDTF surface.
    // }

    // SUBTEST #13 (RegFeat.t:66-67): `throws_ok check_cell_types
    // qr/unavailable/`. Same gap as #12.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_check_cell_types_throws_for_unknown_subtest_13() {
    //     // Future: annotate_vep(..., cell_type => ['nonexistent']) should error.
    // }

    // SUBTEST #22 (RegFeat.t:102): `is(scalar @$features, 1771)`.
    // Blocked: vepyr has no public `count_features_in_region(...)` API.
    // Future-work entry: `RegFeatCache::count_features_in_region()`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_flatten_count_is_1771_subtest_22() {
    //     // Future signature:
    //     //   pub fn count_features_in_region(table: &str, chrom: &str,
    //     //                                    start: i64, end: i64)
    //     //       -> Result<usize>;
    //     // Note: the 1771 value is v84-fixture-specific; the v115 count
    //     // would differ. Hand-coded expected derived from v115 cache at
    //     // commit time when the API lands.
    // }

    // SUBTEST #23 (RegFeat.t:104): merge_features dedup count 1605 / 1771.
    // Blocked: vepyr has no public `merge_features(…)` dedup helper.
    // Future-work entry: `merge_features dedup helper`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_merge_features_dedup_count_subtest_23() {
    //     // Future signature:
    //     //   pub fn merge_features<T: HasStableId>(feats: Vec<T>) -> Vec<T>;
    //     // Dedup count is v84-fixture-specific (1605/1771); v115 count
    //     // derived at commit time.
    // }

    // SUBTEST #28 (RegFeat.t:119-120): `$c->clean_cache(); is_deeply($c->cache, {})`.
    // Blocked: vepyr has no public eviction API.
    // Future-work entry: `clear_session_cache()`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_clean_cache_empties_subtest_28() {
    //     // Future signature:
    //     //   pub fn clear_session_cache(session: &SessionContext) -> usize;
    // }

    // SUBTEST #39 (RegFeat.t:146): `is(scalar @$features, 532)`.
    // Same future-work entry as #22 (count_features_in_region).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_buffer_feature_count_is_532_subtest_39() {
    //     // Future: count_features_in_region(table, chr21, buffer_start, buffer_end).
    //     // 532 is v84-fixture-specific; v115 count derived at commit time.
    // }
}
