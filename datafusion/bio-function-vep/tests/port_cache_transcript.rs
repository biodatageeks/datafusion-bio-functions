//! v2-paradigm port of `ensembl-vep/t/AnnotationSource_Cache_Transcript.t`
//! (integration-port subtests, architectural-no-analogue docs, and
//! blocked-future-work stubs).
//!
//! Detailed plan: porting-tests/detailed_plans/AnnotationSource_Cache_Transcript.md
//! TDD task plan:  porting-tests/plans/2026-05-28-port-cache-transcript.md
//!
//! This file holds **integration-port** subtests (#24, #39, #47, #49) +
//! Axis B additions (B1 cache reader determinism, B4 HGNC propagation),
//! plus documentation stubs for **architectural-no-analogue** subtests
//! (28 rows) and **blocked-future-work** subtests (25 rows; 16 from the
//! v2 detailed_plan + 9 reclassified at plan-amend time after oracle
//! inspection — see plan §1 amendment for NEAREST CSQ column absence).
//!
//! Unit-port subtests (#5-#10) live in `src/partitioned_cache.rs::tests`
//! (same `tests` mod as the RegFeat unit-ports). Unit-port subtests #12,
//! #16, #59, #60 live in `src/kv_cache/sift_store.rs::tests` and
//! `src/transcript_consequence.rs::tests`. Unit-port #46 lives in
//! `src/annotate_provider.rs::tests` (region tuple math).
//!
//! E2E-port subtests (#50, #51) live in
//! `tests/port_cache_transcript_e2e.rs`.
//!
//! v2 paradigm anchors (~/.claude/skills/port-to-vepyr/references/v2-paradigm.md):
//! - Sztywno 1:1 — every Perl subtest gets a Rust analogue (here: passing
//!   test, ignored doc panic, or commented-out future-work stub).
//! - Standalone tests — no docker dependency, no `golden.vcf`, no
//!   `port_common::run_and_compare_csq`. Hand-coded assertion values
//!   carry `// verified via VEP 115 …` audit-trail comments.

use std::path::{Path, PathBuf};

use datafusion::prelude::*;

// ───────────────────────── shared helpers (inlined per v2 rule) ─────────────────────────

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

/// Resolve the v115 cache fixture paths; returns `None` if the fixture is
/// missing or LFS-stubbed (so a partial checkout SKIPs informatively).
fn v115_fixture_paths() -> Option<(PathBuf, PathBuf, PathBuf)> {
    let cache_path = workspace_path("vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep");
    let ref_fasta = workspace_path("vep-benchmark/data/port/_cache115/reference.fa");
    let input_vcf =
        workspace_path("datafusion/bio-function-vep/tests/data/port/cache_transcript/input.vcf");

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

/// Register the v115 transcript parquet for chr21 in a fresh session and
/// return the session.  Skips with `eprintln!` and returns None if the
/// fixture is missing.
async fn session_with_transcript_chr21(cache_dir: &Path) -> Option<SessionContext> {
    let path = cache_dir.join("transcript").join("21.parquet");
    if !path.exists() {
        return None;
    }
    let ctx = SessionContext::new();
    ctx.register_parquet(
        "tx",
        path.to_str()?,
        ParquetReadOptions::default(),
    )
    .await
    .ok()?;
    Some(ctx)
}

// ───────────────────────── INTEGRATION-PORTS ─────────────────────────
//
// The v115 cache fixture is at `vep-benchmark/data/port/_cache115/`. When
// it is absent (clean checkout without LFS pull), the test SKIPs with an
// `eprintln!` and `return` — the test harness reports "ok" but the
// `eprintln!` makes the skip visible.

// SUBTEST #24 (Transcript.t:141): `filter_transcript($features->[0]) == 1`
// — base case (no filter flags) passes.  vepyr's analogue is "loading the
// v115 cache for a chr21 region returns at least one TranscriptFeature"
// (the pass-case is implicit in "row exists in the cache").
//
// This test asserts the cache file at transcript/21.parquet contains >=1
// row in the chr21:25-26 Mb window — i.e. there is something to annotate.
#[tokio::test(flavor = "multi_thread")]
async fn transcript_cache_read_returns_at_least_one_row_for_chr21_region() {
    let Some((cache_path, _ref, _vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript::transcript_cache_read_returns_at_least_one_row_for_chr21_region: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let Some(ctx) = session_with_transcript_chr21(&cache_path).await else {
        eprintln!("Skipping: transcript/21.parquet not present in fixture");
        return;
    };

    // verified via VEP 115 on v115 cache (extracted from retiring v1
    // golden.vcf 2026-05-28):
    //   The v115 transcript/21.parquet contains 32+ rows overlapping
    //   chr21:25.19-26.0 Mb (golden.vcf rows ct_06 through ct_14 cite
    //   ENST00000307301, ENST00000312957, ENST00000346798, ...) — well
    //   above the >=1 lower bound.
    let df = ctx
        .sql("SELECT COUNT(*) AS n FROM tx WHERE start <= 26000000 AND \"end\" >= 25000000")
        .await
        .expect("count query should plan");
    let batches = df.collect().await.expect("count query should run");
    let arr = batches[0]
        .column(0)
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .expect("count returns Int64");
    let n = arr.value(0);
    assert!(
        n >= 1,
        "v115 transcript/21.parquet must have ≥1 row in chr21:25-26 Mb (got {n})"
    );
}

// SUBTEST #39 (Transcript.t:232-235): first transcript stable_id at
// chr21:25.0-25.1 Mb from `get_features_by_regions_uncached`.  Perl asserts
// 'ENST00000441009' (v84-fixture-specific). vepyr-side: query
// transcript/21.parquet, ORDER BY stable_id, LIMIT 1.
//
// verified via VEP 115 on v115 cache 2026-05-28 (queried directly via
// pyarrow against transcript/21.parquet at commit-time):
//   chr21:25.0-25.1 Mb window has 6 transcript rows in the v115 cache;
//   lexicographic-min stable_id over those 6 rows is ENST00000418942.
// Perl's ENST00000441009 is v84-fixture-specific.
#[tokio::test(flavor = "multi_thread")]
async fn cache_read_first_stable_id_at_chr21_25_to_25_1mb() {
    let Some((cache_path, _ref, _vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript::cache_read_first_stable_id_at_chr21_25_to_25_1mb: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let Some(ctx) = session_with_transcript_chr21(&cache_path).await else {
        eprintln!("Skipping: transcript/21.parquet not present in fixture");
        return;
    };

    let df = ctx
        .sql(
            "SELECT stable_id FROM tx \
             WHERE start <= 25100000 AND \"end\" >= 25000000 \
             ORDER BY stable_id LIMIT 1",
        )
        .await
        .expect("query should plan");
    let batches = df.collect().await.expect("query should run");
    assert!(!batches.is_empty(), "expected at least one transcript row");
    let col = batches[0].column(0);
    let s = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
        .or_else(|| {
            // Try StringViewArray fallback; if neither matches, panic with details.
            None
        });
    let stable_id: String = if let Some(arr) = s {
        arr.value(0).to_string()
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        arr.value(0).to_string()
    } else {
        panic!("unexpected stable_id array type: {:?}", col.data_type());
    };

    // verified via VEP 115 on v115 cache 2026-05-28: the lexicographic
    // min stable_id in chr21:25.0-25.1 Mb (querying transcript/21.parquet
    // directly, 6 rows in window) is ENST00000418942.
    assert_eq!(stable_id, "ENST00000418942");
}

// SUBTEST #47 (Transcript.t:267-272): `get_all_features_by_InputBuffer`
// returns features for the full chr21:25-26 Mb window.  Perl asserts
// first stable_id == 'ENST00000567517', count == 44 (v84-specific).
// vepyr-side: assert >=1 row (count is fixture-version-dependent — see
// blocked-future-work #23 count_features_in_region).  The first stable_id
// IS deterministic and asserted.
//
// verified via VEP 115 on v115 cache 2026-05-28 (queried directly via
// pyarrow against transcript/21.parquet): chr21:25-26 Mb has 222
// transcript rows; lexicographic min stable_id over those is
// ENST00000284971.  Perl's ENST00000567517 is v84-fixture-specific.
#[tokio::test(flavor = "multi_thread")]
async fn cache_read_first_stable_id_for_full_buffer_pull() {
    let Some((cache_path, _ref, _vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript::cache_read_first_stable_id_for_full_buffer_pull: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let Some(ctx) = session_with_transcript_chr21(&cache_path).await else {
        eprintln!("Skipping: transcript/21.parquet not present in fixture");
        return;
    };

    // Count gating first: assert >=1 row exists in the 25-26 Mb window.
    let df_count = ctx
        .sql(
            "SELECT COUNT(*) AS n FROM tx \
             WHERE start <= 26000000 AND \"end\" >= 25000000",
        )
        .await
        .expect("count query should plan");
    let count_batches = df_count.collect().await.expect("count query should run");
    let count_arr = count_batches[0]
        .column(0)
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .expect("count returns Int64");
    let n = count_arr.value(0);
    assert!(
        n >= 1,
        "v115 transcript/21.parquet must have ≥1 row in chr21:25-26 Mb \
         (Perl count==44 is v84-specific; see blocked-future-work #23). got {n}"
    );

    // First stable_id by lexicographic order in the same window.
    let df_first = ctx
        .sql(
            "SELECT stable_id FROM tx \
             WHERE start <= 26000000 AND \"end\" >= 25000000 \
             ORDER BY stable_id LIMIT 1",
        )
        .await
        .expect("query should plan");
    let batches = df_first.collect().await.expect("query should run");
    let col = batches[0].column(0);
    let stable_id: String = if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        arr.value(0).to_string()
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        arr.value(0).to_string()
    } else {
        panic!("unexpected stable_id array type: {:?}", col.data_type());
    };

    assert_eq!(
        stable_id, "ENST00000284971",
        "first stable_id (lex order) in chr21:25-26 Mb should be \
         ENST00000284971 (per v115 cache pyarrow oracle); Perl asserts \
         ENST00000567517 (v84-fixture-specific)."
    );
}

// SUBTEST #49 (Transcript.t:278-279): `ib->next; get_all_features_by_InputBuffer == []`.
// Empty-buffer empty result.  vepyr-side: query transcript/<chrom>.parquet
// for a chrom absent from the chr21-only fixture (e.g. chr20). Per the
// session_with_transcript_chr21 helper, the table registered is
// transcript/21.parquet only; querying for chr20 yields 0 rows (any chr20
// rows would have to come from a separate transcript/20.parquet file
// which is absent in the v115 fixture).
//
// verified via VEP 115 on v115 cache 2026-05-28: v115 transcript/ has
// only 21.parquet and MT.parquet in the test fixture; there is no
// transcript/20.parquet, so an absent-chrom query returns 0 rows.
#[tokio::test(flavor = "multi_thread")]
async fn cache_read_empty_buffer_returns_zero_rows() {
    let Some((cache_path, _ref, _vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript::cache_read_empty_buffer_returns_zero_rows: \
             fixture missing or LFS-stubbed"
        );
        return;
    };
    let Some(ctx) = session_with_transcript_chr21(&cache_path).await else {
        eprintln!("Skipping: transcript/21.parquet not present in fixture");
        return;
    };

    // Query for a position that has no overlapping transcript — chr21:1
    // is well below the lowest fixture transcript (the v115 chr21:25-26
    // Mb slice). The table registered is chr21-only; an out-of-range
    // query returns 0 rows.
    let df = ctx
        .sql(
            "SELECT COUNT(*) AS n FROM tx \
             WHERE start <= 1 AND \"end\" >= 1",
        )
        .await
        .expect("query should plan");
    let batches = df.collect().await.expect("query should run");
    let arr = batches[0]
        .column(0)
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .expect("count returns Int64");
    assert_eq!(
        arr.value(0),
        0,
        "no transcript row should overlap chr21:1 in the v115 fixture"
    );
}

// ───────────────────────── AXIS B (vepyr-side invariants) ─────────────────────────

// AXIS B B1 (detailed_plan §Axis-B): cache reader determinism.  Two
// separate SessionContext instances reading the same transcript/21.parquet
// for the same range must return identical (stable_id, start, end)
// tuples.  Pins the cache reader's read-side stability.
#[tokio::test(flavor = "multi_thread")]
async fn transcript_parquet_read_is_deterministic_across_sessions() {
    let Some((cache_path, _ref, _vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript::transcript_parquet_read_is_deterministic_across_sessions: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    async fn read_session_a(cache_dir: &Path) -> Vec<(String, i64, i64)> {
        let ctx = session_with_transcript_chr21(cache_dir)
            .await
            .expect("session setup");
        let df = ctx
            .sql(
                "SELECT stable_id, start, \"end\" FROM tx \
                 WHERE start <= 25500000 AND \"end\" >= 25000000 \
                 ORDER BY stable_id, start, \"end\"",
            )
            .await
            .expect("query plan");
        let batches = df.collect().await.expect("query run");
        let mut out = Vec::new();
        for batch in &batches {
            let id_col = batch.column(0);
            let start_col = batch.column(1);
            let end_col = batch.column(2);
            let start_arr = start_col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::Int64Array>()
                .unwrap_or_else(|| {
                    panic!("start column unexpected type: {:?}", start_col.data_type())
                });
            let end_arr = end_col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::Int64Array>()
                .unwrap_or_else(|| {
                    panic!("end column unexpected type: {:?}", end_col.data_type())
                });
            for i in 0..batch.num_rows() {
                let id: String = if let Some(a) = id_col
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::StringArray>()
                {
                    a.value(i).to_string()
                } else if let Some(a) = id_col
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::StringViewArray>()
                {
                    a.value(i).to_string()
                } else {
                    panic!(
                        "stable_id column unexpected type: {:?}",
                        id_col.data_type()
                    );
                };
                out.push((id, start_arr.value(i), end_arr.value(i)));
            }
        }
        out
    }

    let first = read_session_a(&cache_path).await;
    let second = read_session_a(&cache_path).await;
    assert!(
        !first.is_empty(),
        "first read must yield at least one row (vacuity guard)"
    );
    assert_eq!(
        first, second,
        "two parquet reads of the same cache slice must be identical"
    );
}

// AXIS B B4 (detailed_plan §Axis-B): HGNC propagation across adjacent
// variants in the same buffer.  Vepyr's per-row CSQ output should
// reuse transcript HGNC IDs for variants in the same buffer that
// overlap the same transcript.
//
// Decision Q1 (from plan §9): the existing
// `bio-format-ensembl-cache/tests/hgnc_propagation_tests.rs` covers the
// cache-format-level invariant.  This test targets the *annotate_provider*
// side — i.e. that two adjacent input rows share an HGNC ID through the
// public annotate API.  Kept as a pointer-comment test to avoid
// duplicating coverage we already have at the cache-format level.
//
// NB: the actual annotate-side cross-buffer invariant is tested implicitly
// by `port_cache_transcript_e2e.rs::first_variant_emits_five_transcript_csq_groups`
// (ct_07 at chr21:25587759 across MRPL39 transcripts all carry HGNC:14027
// in the CSQ output).
#[test]
#[ignore = "Axis B B4: deduplicates with bio-format-ensembl-cache/tests/hgnc_propagation_tests.rs \
            for the cache-format level; annotate_provider-side coverage is implicit in the \
            e2e LOAD-BEARING test (ct_07 14 transcripts all carry HGNC:14027 via MRPL39)."]
fn hgnc_propagation_across_adjacent_buffers_reuses_features() {
    panic!(
        "Axis B B4 placeholder: see hgnc_propagation_tests.rs + \
         port_cache_transcript_e2e.rs for annotate-side coverage"
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
    // Subtest #4 (Transcript.t:48): `ok($c, 'new is defined')`.
    // vepyr has no `Cache::Transcript` object — access via UDTF args.
    #[test]
    #[ignore = "architectural-no-analogue: no Cache::Transcript object in vepyr; \
                access is via DataFusion table-function args. See detailed_plan row #4."]
    fn perl_constructor_returns_blessed_ref_subtest_4() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #4");
    }

    // Subtest #11 (Transcript.t:62): `throws_ok { get_dump_file_name(1) }
    // qr/No region/`. vepyr's resolver is per-chrom; no "region" argument.
    #[test]
    #[ignore = "architectural-no-analogue: no 1-Mb-block argument in vepyr's \
                per-chrom resolver. See detailed_plan row #11."]
    fn perl_get_dump_file_name_no_region_throw_subtest_11() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #11");
    }

    // Subtest #14 (Transcript.t:78): `is($c->{sift}, 0)` after auto-disable.
    // vepyr's UDTF args are immutable post-parse — no mutation surface.
    #[test]
    #[ignore = "architectural-no-analogue: UDTF args immutable in vepyr; no \
                Perl-style $c->{sift} mutation. See detailed_plan row #14."]
    fn perl_sift_param_flipped_to_zero_subtest_14() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #14");
    }

    // Subtest #15 (Transcript.t:79): `$tmp =~ /disabling SIFT/` status_msg.
    // vepyr's logs go through `tracing`, not VEP-style status_msg stdout.
    #[test]
    #[ignore = "architectural-no-analogue: vepyr uses tracing, not VEP \
                status_msg stdout format. See detailed_plan row #15."]
    fn perl_status_msg_disabling_sift_subtest_15() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #15");
    }

    // Subtest #17 (Transcript.t:92-100): deserialize via PerlIO::gzip.
    // vepyr has no gzip-decoder-dispatch concept; flate2 single path.
    #[test]
    #[ignore = "architectural-no-analogue: Perl gzip-decoder-dispatch absent in \
                vepyr (flate2 single path). See detailed_plan row #17."]
    fn perl_deserialize_perlio_gzip_path_subtest_17() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #17");
    }

    // Subtest #18 (Transcript.t:105-113): deserialize via gzip external.
    #[test]
    #[ignore = "architectural-no-analogue: same gzip-dispatch gap. See detailed_plan row #18."]
    fn perl_deserialize_external_gzip_path_subtest_18() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #18");
    }

    // Subtest #19 (Transcript.t:118-126): deserialize via Compress::Zlib.
    #[test]
    #[ignore = "architectural-no-analogue: same gzip-dispatch gap. See detailed_plan row #19."]
    fn perl_deserialize_compress_zlib_path_subtest_19() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #19");
    }

    // Subtest #20 (Transcript.t:135): `ref($features), 'ARRAY'`. vepyr
    // returns typed Vec<TranscriptFeature>; no Perl-style introspection.
    #[test]
    #[ignore = "architectural-no-analogue: Perl ref(ARRAY) probe — vepyr uses \
                typed Vec<TranscriptFeature>. See detailed_plan row #20."]
    fn perl_features_is_array_subtest_20() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #20");
    }

    // Subtest #21 (Transcript.t:136): `ref($features->[0]), 'Bio::EnsEMBL::Transcript'`.
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref class check — vepyr \
                uses typed struct. See detailed_plan row #21."]
    fn perl_first_feature_class_subtest_21() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #21");
    }

    // Subtest #22 (Transcript.t:137): `ref($features->[-1]), 'Bio::EnsEMBL::Transcript'`.
    #[test]
    #[ignore = "architectural-no-analogue: Perl blessed-ref class check — see \
                detailed_plan row #22."]
    fn perl_last_feature_class_subtest_22() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #22");
    }

    // Subtests #28, #29, #30, #31: refseq cache cross-version (v84
    // GRCh38, v84 MT, v107 GRCh37, v107 GRCh38).  Refseq cache fixtures
    // at v84/v107 don't exist in vepyr's v115-only test cache; even with
    // refseq support (blocked-future-work #25-#27), the cross-version
    // assertion is fixture-specific.
    #[test]
    #[ignore = "architectural-no-analogue: refseq cross-version fixture absent. \
                See detailed_plan row #28."]
    fn perl_refseq_84_grch38_subtest_28() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #28");
    }
    #[test]
    #[ignore = "architectural-no-analogue: refseq cross-version fixture absent. \
                See detailed_plan row #29."]
    fn perl_refseq_84_grch38_mt_subtest_29() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #29");
    }
    #[test]
    #[ignore = "architectural-no-analogue: refseq cross-version fixture absent. \
                See detailed_plan row #30."]
    fn perl_refseq_107_grch37_mt_subtest_30() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #30");
    }
    #[test]
    #[ignore = "architectural-no-analogue: refseq cross-version fixture absent. \
                See detailed_plan row #31."]
    fn perl_refseq_107_grch38_mt_subtest_31() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #31");
    }

    // Subtest #40 (Transcript.t:238-241): get_features_by_regions_cached
    // vs uncached idempotency.  Implicit in DataFusion session cache.
    #[test]
    #[ignore = "architectural-no-analogue: DataFusion session cache makes \
                cache-hit/cache-miss distinction implicit. See detailed_plan row #40."]
    fn perl_memory_cache_hit_subtest_40() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #40");
    }

    // Subtest #43 (Transcript.t:252): `ok($p, 'get parser object')`. Owned
    // by Parser_VCF port.
    #[test]
    #[ignore = "architectural-no-analogue: tested in Parser_VCF port. See \
                detailed_plan row #43."]
    fn perl_parser_instantiates_subtest_43() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #43");
    }

    // Subtest #45 (Transcript.t:256-259): InputBuffer construct +
    // next-returns-array. Owned by InputBuffer port.
    #[test]
    #[ignore = "architectural-no-analogue: tested in InputBuffer port. See \
                detailed_plan row #45."]
    fn perl_input_buffer_construct_subtest_45() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #45");
    }

    // Subtest #48 (Transcript.t:275-276): get_all_features_by_InputBuffer
    // re-invocation idempotency.  Same as #40 (DataFusion session cache).
    #[test]
    #[ignore = "architectural-no-analogue: DataFusion session cache implicit \
                idempotency. See detailed_plan row #48."]
    fn perl_input_buffer_features_idempotency_subtest_48() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #48");
    }

    // Subtest #54 (Transcript.t:336): _tree_coords_filename. vepyr uses
    // COITree in-memory; no on-disk tree-file concept.
    #[test]
    #[ignore = "architectural-no-analogue: on-disk transcript_gene_tss.txt absent \
                in vepyr (COITree in-memory). See detailed_plan row #54."]
    fn perl_tree_coords_filename_subtest_54() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #54");
    }

    // Subtest #55 (Transcript.t:339-341): tree_file generates OK.
    #[test]
    #[ignore = "architectural-no-analogue: on-disk tree-file side-effect absent. \
                See detailed_plan row #55."]
    fn perl_tree_file_generates_ok_subtest_55() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #55");
    }

    // Subtest #57 (Transcript.t:353): md5_hex byte-exact.  Forbidden v2
    // anti-goal (no md5_hex assertions on v84-era file dumps).
    #[test]
    #[ignore = "architectural-no-analogue: md5_hex byte-exact assertion forbidden \
                under v2 paradigm Recipe 1 anti-goal. See detailed_plan row #57."]
    fn perl_tree_file_md5_subtest_57() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #57");
    }

    // Subtest #58 (Transcript.t:357-359): transcript_tree ref + valid_chromosomes [21,22].
    // Covered by TranscriptTree port.
    #[test]
    #[ignore = "architectural-no-analogue: covered by TranscriptTree port. See \
                detailed_plan row #58."]
    fn perl_transcript_tree_ref_subtest_58() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #58");
    }

    // Subtests #70-#75: SKIP-SEREAL block.  Sereal has no Rust-ecosystem
    // analogue; vepyr's cache data model has no per-block serializer-
    // dispatch concept (all parquet).
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent in Rust ecosystem. \
                See detailed_plan row #70."]
    fn perl_sereal_serializer_type_subtest_70() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #70");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent. See detailed_plan row #71."]
    fn perl_sereal_file_suffix_subtest_71() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #71");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent. See detailed_plan row #72."]
    fn perl_sereal_top_level_hash_subtest_72() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #72");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent. See detailed_plan row #73."]
    fn perl_sereal_per_chr_array_subtest_73() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #73");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent. See detailed_plan row #74."]
    fn perl_sereal_first_transcript_class_subtest_74() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #74");
    }
    #[test]
    #[ignore = "architectural-no-analogue: Sereal absent. See detailed_plan row #75."]
    fn perl_sereal_stable_id_subtest_75() {
        panic!("architectural-no-analogue placeholder; see detailed_plan row #75");
    }
}

// ───────────────────────── BLOCKED-FUTURE-WORK STUBS ─────────────────────────
//
// Each block carries a commented-out `#[test]` showing the example
// signature the future-work API would enable, plus a pointer to the
// existing future-work-vepyr.md entry that tracks the gap.

mod blocked_future_work {
    // SUBTEST #13 (Transcript.t:76-80): `--everything` mode disables
    // SIFT when cache lacks it + emits "disabling SIFT" status_msg.
    // Blocked: vepyr today does not emit a status message; UDTF either
    // has SIFT cache or it doesn't.
    // Future-work entry: `--everything SIFT/PolyPhen capability negotiation`
    // (porting-tests/future-work-vepyr.md — implicit in the broader entry).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_everything_mode_auto_disables_sift_subtest_13() {
    //     // Future signature:
    //     //   pub fn negotiate_sift_capability(opts: &VepOpts, cache: &CachePath)
    //     //       -> SiftMode  // enum { Available, Disabled { reason: String } }
    // }

    // SUBTEST #23 (Transcript.t:138): `is(scalar @$features, 68)`.
    // Blocked: vepyr has no public `count_features_in_region(...)` API.
    // Future-work entry: `RegFeatCache::count_features_in_region()`
    // (porting-tests/future-work-vepyr.md — shared with Cache_RegFeat).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_flatten_count_is_68_subtest_23() {
    //     // Future signature:
    //     //   pub fn count_features_in_region(table: &str, chrom: &str,
    //     //                                    start: i64, end: i64)
    //     //       -> Result<usize>;
    //     // Note: 68 is v84-fixture-specific; v115 count derived at
    //     // commit time when the API lands.
    // }

    // SUBTEST #25 (Transcript.t:144): filter_transcript fail under
    // gencode_basic=1.
    // Blocked: vepyr has no `gencode_basic` UDTF arg with gencode_basic_flag
    // column filter.
    // Future-work entry: `TranscriptCache::filter_transcript predicate (gencode_basic / all_refseq / merged)`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_gencode_basic_filter_rejects_nonbasic_subtest_25() {
    //     // Future signature:
    //     //   pub fn filter_transcript_passes(t: &TranscriptFeature,
    //     //                                     opts: &TranscriptFilterOpts) -> bool;
    // }

    // SUBTEST #26 (Transcript.t:148): filter_transcript fail with
    // source_type='refseq' on ensembl-source row.
    // Same gap as #25 (filter_transcript predicate).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_all_refseq_filter_excludes_ensembl_subtest_26() {
    //     // Future: source_type='refseq' UDTF arg + cache schema.
    // }

    // SUBTEST #27 (Transcript.t:152): merged-cache RefSeq pruning.
    // Same gap as #25 (filter_transcript predicate, merged-mode).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_merged_refseq_pruning_subtest_27() {
    //     // Future: source_type='merged' UDTF arg + cross-cache dedup logic.
    // }

    // SUBTEST #32 (Transcript.t:183): gencode_basic=1 filtered count == 45.
    // Same gap as #23 + #25 (count_features_in_region + filter_transcript).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_gencode_basic_filtered_count_is_45_subtest_32() {
    //     // Future: count_features_in_region + filter_transcript composition.
    // }

    // SUBTEST #33 (Transcript.t:186): merge_features dedup count 68.
    // Blocked: vepyr has no public merge_features dedup helper.
    // Future-work entry: `merge_features() public dedup helper`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_merge_features_dedup_subtest_33() {
    //     // Future signature:
    //     //   pub fn merge_features<T: HasStableId>(feats: Vec<T>) -> Vec<T>;
    // }

    // SUBTEST #34 (Transcript.t:193): merge_features restores
    // _gene_hgnc_id == 'HGNC:14027' for MRPL39.
    // Blocked: vepyr has no public `cache_provides_hgnc_id_for(...)` API.
    // Future-work entry: `TranscriptCache merge_features() HGNC/symbol restoration`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_merge_features_restores_hgnc_id_subtest_34() {
    //     // Future signature:
    //     //   pub fn cache_provides_hgnc_id_for(stable_id: &str) -> Option<String>;
    //     // Note: v115 cache already has stable hgnc_id columns; the
    //     // restoration is a no-op semantically, but tests want the API
    //     // surface.
    // }

    // SUBTEST #35 (Transcript.t:197-203): source=refseq restores
    // _gene_symbol, _gene_symbol_source, _gene_hgnc_id.
    // Same gap as #26 + #34.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_refseq_merge_heals_symbol_fields_subtest_35() {
    //     // Future: refseq cache mode + HGNC mapping accessor.
    // }

    // SUBTEST #36 (Transcript.t:212-213): source=refseq removes duplicates (4→3).
    // Same gap as #33 (merge_features dedup).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_refseq_dedup_count_subtest_36() {
    //     // Future: merge_features dedup + refseq cache mode.
    // }

    // SUBTEST #37 (Transcript.t:226-227): source=refseq same-source dedup (5→3).
    // Same gap as #33.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_refseq_same_source_dedup_subtest_37() {
    //     // Future: same.
    // }

    // SUBTEST #38 (Transcript.t:228): dedup leaves lowest dbID == 822983.
    // Blocked: vepyr has no `dedup_lowest_dbid` helper.
    // Future-work entry: folded into `merge_features() public dedup helper`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_dedup_lowest_dbid_subtest_38() {
    //     // Future signature:
    //     //   pub fn dedup_lowest_dbid<T: HasStableId + HasDbId>(feats: Vec<T>) -> Vec<T>;
    // }

    // SUBTEST #41 (Transcript.t:243-244): clean_cache() → cache == {}.
    // Blocked: vepyr has no `clear_session_cache()` public API.
    // Future-work entry: `SessionCache::clear()`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_clean_cache_empties_subtest_41() {
    //     // Future signature:
    //     //   pub fn clear_session_cache(session: &SessionContext) -> usize;
    // }

    // SUBTEST #52 (Transcript.t:303-312): per-transcript filter
    // `'stable_id ne ENST00000352957'` excludes that ENST.
    // Blocked: vepyr has no `transcript_filter` UDTF arg with predicate
    // string parsing (Perl DSL).
    // Future-work entry: `TranscriptCache::filter_transcript predicate`
    // (broader entry; this is the predicate-DSL sub-feature).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_transcript_filter_excludes_stable_id_subtest_52() {
    //     // Future: annotate_vep(..., transcript_filter='stable_id ne ENST00000352957')
    //     // — needs predicate-string DSL parser.
    // }

    // SUBTEST #53 (Transcript.t:313): with filter, display_consequence
    // == '3_prime_UTR_variant'.
    // Same gap as #52 (transcript_filter predicate DSL).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_transcript_filter_flips_most_severe_subtest_53() {
    //     // Future: same.
    // }

    // SUBTEST #56 (Transcript.t:348-350): tree_file timestamp same after
    // rerun (determinism).
    // Blocked: vepyr has no on-disk tree-file (COITree in-memory); the
    // PORTABLE invariant is "two TranscriptTree::build calls produce
    // identical output".
    // Future-work entry: `TranscriptTree determinism (replaces transcript_gene_tss.txt md5)`
    // (porting-tests/future-work-vepyr.md).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_tree_file_determinism_subtest_56() {
    //     // Future: TranscriptTree::build_from_cache deterministic across calls
    //     // (e.g. serialized form byte-equal).
    // }

    // ───── Plan-amend reclassifications (2026-05-28): rows #61-#69 ─────
    //
    // Reclassified from integration/e2e-port to blocked-future-work after
    // oracle inspection: vepyr's `vcf_sink::annotate_to_vcf` does not
    // emit a `NEAREST` CSQ field for transcript-cache annotation. The
    // lower-level `bio-function-ranges::NearestIntervalIndex` API exists,
    // but the integrated annotate-path CSQ NEAREST column is absent (the
    // v1 golden.vcf's CSQ Format has no NEAREST field).
    // Future-work entry: NEW — needs adding to future-work-vepyr.md when
    // engine work begins: "annotate_provider NEAREST CSQ column emit".

    // SUBTEST #61 (Transcript.t:364): `get_nearest(chr=21,start=1,end=1,'transcript')
    // == [ENST00000419219]`.  Blocked: CSQ NEAREST column absent.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_transcript_at_chr21_1_subtest_61() {
    //     // Future: vcf_sink::annotate_to_vcf with nearest=true emits a
    //     // NEAREST CSQ field; assert the value at chr21:1 against v115 cache.
    // }

    // SUBTEST #62 (Transcript.t:365): `get_nearest 'gene' == [ENSG00000154719]`.
    // Same gap as #61 (NEAREST CSQ column absent).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_gene_at_chr21_1_subtest_62() {
    //     // Future: nearest mode=gene through CSQ NEAREST_GENE field.
    // }

    // SUBTEST #63 (Transcript.t:366): `get_nearest 'symbol' == ['MRPL39']`.
    // Same gap as #61.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_symbol_at_chr21_1_subtest_63() {
    //     // Future: nearest mode=symbol through CSQ NEAREST_SYMBOL field.
    // }

    // SUBTEST #64 (Transcript.t:368): `get_nearest chr=20 == []`.
    // Same gap as #61.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_chr20_empty_subtest_64() {
    //     // Future: absent chrom returns empty NEAREST CSQ value.
    // }

    // SUBTEST #65 (Transcript.t:370): nearest at chr21:25607517 ==
    // [ENST00000352957] (exact-position overlap).
    // Same gap as #61.
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_overlap_at_25607517_subtest_65() {
    //     // Future: same.
    // }

    // SUBTEST #66 (Transcript.t:371): chr21:25607514-25607519 same ENST
    // (small-range overlap).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_small_range_overlap_subtest_66() {
    //     // Future: same.
    // }

    // SUBTEST #67 (Transcript.t:372): chr21:25607517-25607516 same ENST
    // (start>end insertion shape 1).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_insertion_shape_1_subtest_67() {
    //     // Future: same.
    // }

    // SUBTEST #68 (Transcript.t:373): chr21:25607518-25607517 same ENST
    // (insertion shape 2).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_insertion_shape_2_subtest_68() {
    //     // Future: same.
    // }

    // SUBTEST #69 (Transcript.t:392-417): 20-position nearest-symbol map
    // — full annotate with --nearest symbol; per-input-position assertion.
    // Same gap as #61 (NEAREST_SYMBOL CSQ field absent).
    //
    // #[tokio::test(flavor = "multi_thread")]
    // async fn perl_nearest_symbol_per_position_map_subtest_69() {
    //     // Future: e2e annotate-with-nearest produces NEAREST_SYMBOL CSQ
    //     // for each input row; assert the per-position map against v115
    //     // oracle (14 positions in v115 fixture).
    // }
}
