//! Port of `ensembl-vep/t/AnnotationSource_Cache_VariationTabix.t` (v2 paradigm).
//!
//! See `porting-tests/detailed_plans/AnnotationSource_Cache_VariationTabix.md`
//! and `porting-tests/plans/2026-05-27-port-cache-variation-tabix.md` for the
//! full per-subtest classification table (36 substantive subtests; 7
//! architectural-no-analogue rows; 22 blocked-future-work rows folded into
//! existing future-work entries; 7 cross-references into the sibling
//! `Cache_Variation` port).
//!
//! ## Why this file is mostly comments (no active `#[test]`s)
//!
//! vepyr unifies Perl's TWO variation-cache storage shapes into a SINGLE
//! storage layout:
//!
//! * Perl shape A — per-region storable dumps (`<dir>/<chr>/<region>_var.gz`)
//!   exercised by `AnnotationSource_Cache_Variation.t`.
//! * Perl shape B — single tabix-indexed bgzip per chrom
//!   (`<dir>/<chr>/all_vars.gz` + `all_vars.gz.tbi`) exercised by THIS
//!   file's source `AnnotationSource_Cache_VariationTabix.t`.
//! * vepyr shape — per-chrom parquet (`<cache>/variation/<chrom>.parquet`)
//!   accessed via `EnsemblCacheTableProvider::for_entity("variation", …)`
//!   and queried via `VariantLookupExec`.
//!
//! Storage shape differs; OBSERVABLE CONTRACT (`Existing_variation`, `AF`,
//! `gnomADe_*`, `CLIN_SIG`, `PHENO`, `PUBMED`, `SOMATIC` CSQ fields populated
//! from a matched cache row) is identical between Perl shapes A and B. Under
//! sztywno-1:1 ("one observable contract → one Rust assertion") the right
//! port shape is **cross-references** into the canonical sibling port plus a
//! small number of comment stubs documenting tabix-storage mechanics that
//! vepyr replaces wholesale.
//!
//! ## Why zero `#[test]` activations
//!
//! Duplicating Cache_Variation's per-subtest assertions into a second
//! Rust file would violate sztywno-1:1. The 7 observable-contract
//! subtests in VariationTabix.t (#6 empty-chrom path resolver, #13/16
//! rs142513484 hash, #15 miss-by-coord, #18 multi-existing
//! CM930033+rs63750066, #36/38 nastiness 1 & 3) each cross-reference a
//! Rust assertion that already lives in the v2-rewrite sibling files:
//!   - `tests/port_annotation_source_cache_variation_e2e.rs` for e2e
//!     assertions (rs142513484, rs63750066, nastiness 1-4).
//!   - `tests/port_annotation_source_cache_variation.rs` for
//!     integration assertions (miss-by-one, NULL allele, etc.).
//!   - `src/partitioned_cache.rs::tests` for the unit-port empty-chrom
//!     resolver (subtest #6).
//!
//! Cross-references below were repointed 2026-05-28 when the v1
//! `port_cache_variation.rs::port_cache_variation_csq_matches_golden`
//! HARD_FIELDS test was retired in favour of per-Perl-subtest hand-coded
//! assertions in the new v2 files.
//!
//! ## Coverage parity
//!
//! ~78% cross-reference-aware (per detailed_plan §Coverage parity):
//!   - 7 cross-referenced observable-contract rows, each backed by an
//!     existing live Rust assertion in the v2 cache_variation sibling
//!     port (rs142513484_existing_variation_populated_36,
//!     miss_by_one_position_has_empty_existing_variation_38,
//!     rs63750066_multi_existing_with_CM930033_NULL_allele_39,
//!     nastiness_1_indel_context_insertion_40,
//!     nastiness_3_shared_prefix_suffix_trim_42) +
//!     `regfeat_context_path_with_empty_chrom_returns_none`.
//!   - 7 architectural-no-analogue rows (permanent by-design gaps; see
//!     "Architectural-no-analogue" comment blocks below).
//!   - 22 blocked-future-work rows folded into EXISTING future-work entries
//!     (`filter_common`, `ChrSynonyms`, `count_features_in_region`, engine
//!     blocker #1, `AnnotateVcfConfig::resolve invalid-value validation`,
//!     `VariationCache::compare_existing` API); NO new entries are appended
//!     to `future-work-vepyr.md` by this port.
//!
//! ## Anti-goals (do NOT change in this file)
//!
//!   - DO NOT add `assert_eq!(csq.get("Existing_variation"), …)` here.
//!     The canonical per-Perl-subtest assertions live in:
//!     - `tests/port_annotation_source_cache_variation_e2e.rs` (e2e)
//!     - `tests/port_annotation_source_cache_variation.rs` (integration)
//!   - DO NOT activate any `#[test]` here. Every test below is a
//!     cross-reference or a commented-out future-work stub.
//!   - DO NOT add tabix-related future-work entries (`_get_tabix_obj`,
//!     `.tbi` index, PM-vs-CL dispatch). vepyr's data model precludes
//!     tabix storage at the architectural level (A1 RESOLVED 2026-05-25:
//!     parquet+COITree is the single shape).
//!   - DO NOT create `golden.vcf` (v2 retired).
//!   - DO NOT call `port_common::run_and_compare_csq` here (retired
//!     with v1 paradigm; cf. `port_cache_variation.rs` v1 retirement
//!     2026-05-28).

// =============================================================================
// Cross-reference pointer comments (7 covered observable-contract rows)
// =============================================================================
//
// Each row below documents a Perl subtest in VariationTabix.t whose
// observable contract is identical to a subtest in the sibling Cache_Variation
// port. The actual Rust assertion is owned by the sibling test file; the
// pointers here document the cross-reference for audit-trail purposes.
//
// -----------------------------------------------------------------------------
// SUBTEST #6 (L53): throws_ok get_dump_file_name() qr/No chromosome/
//   Perl: empty-args call to get_dump_file_name dies with "No chromosome".
//   Vepyr observable contract: PartitionedParquetCache::context_path returns
//   None for empty chrom.
//   Cross-reference (canonical Rust assertion):
//     src/partitioned_cache.rs::tests::regfeat_context_path_with_empty_chrom_returns_none
//     (line ~350; already GREEN under v1).
//     Asserts: cache.context_path("regulatory", "").is_none()
//   The "variation" context-type variant is covered by the same test through
//   the polymorphic context_path signature (context_type is a `&str`).
//   Status: GREEN (no work required).
//
// -----------------------------------------------------------------------------
// SUBTEST #13 (L122-129): _annotate_cl rs142513484 full existing hash
//   Perl: tabix-CL backend annotates chr21:25585733 with full per-pop AF /
//   gnomADe / minor_allele / phenotype hash for rs142513484.
//   Cross-reference (canonical Rust assertion):
//     tests/port_annotation_source_cache_variation_e2e.rs::
//       rs142513484_existing_variation_populated_36
//     (sibling detailed_plan subtest #36).
//     Asserts Existing_variation populated for chr21:25585733 C->T via
//     `annotate_to_vcf` over the v115 cache. Per-pop AF cells surface
//     through the same CSQ Format header per the v115 oracle.
//   Storage-shape distinction (tabix vs per-region storable) is invisible at
//   the CSQ assertion level — vepyr exercises ONE parquet path that both
//   Perl tests' contracts validate.
//   Status: GREEN under v2 (2026-05-28 rewrite).
//
// -----------------------------------------------------------------------------
// SUBTEST #15 (L146-154): vf->{chr}=21; vf->{start}++ miss-by-coord
//   Perl: chr21:25585734 (one position past rs142513484) → existing is undef.
//   Cross-reference (canonical Rust assertion):
//     tests/port_annotation_source_cache_variation.rs::
//       miss_by_one_position_has_empty_existing_variation_38
//     (sibling detailed_plan subtest #38).
//     Asserts CSQ Existing_variation is empty at chr21:25585734.
//   Status: GREEN under v2 (2026-05-28 rewrite).
//
// -----------------------------------------------------------------------------
// SUBTEST #16 (L172): annotate_InputBuffer rs142513484
//   Perl: full integration via annotate_InputBuffer (vs _annotate_cl direct
//   call in #13).
//   Cross-reference (canonical Rust assertion):
//     tests/port_annotation_source_cache_variation_e2e.rs::
//       rs142513484_existing_variation_populated_36
//     (same as #13 via a different Perl entry point).
//     Same observable contract; no new Rust assertion needed.
//   Status: GREEN under v2.
//
// -----------------------------------------------------------------------------
// SUBTEST #18 (L177-257): chr21:25891796 phenotype_or_disease +
//                          CM930033 NULL allele + rs63750066 multi-existing
//   Perl: clinically annotated variant with multiple existing records (CM930033
//   HGMD_MUTATION + rs63750066) at chr21:25891796.
//   Cross-reference (canonical Rust assertion):
//     tests/port_annotation_source_cache_variation_e2e.rs::
//       rs63750066_multi_existing_with_CM930033_NULL_allele_39
//     (sibling detailed_plan subtest #39 — set-containment of both
//     rs63750066 and CM930033 in Existing_variation).
//   Status: GREEN under v2.
//
// -----------------------------------------------------------------------------
// SUBTEST #36 (L430-443): NASTINESS 1 — chr21:8987005 A→AGCG
//   Perl: matched_alleles=[{a=0, a_allele=GCG, b=0, b_allele=GCG}] —
//   indel-context insertion allele trim via tabix backend.
//   Cross-reference (canonical Rust assertion):
//     tests/port_annotation_source_cache_variation_e2e.rs::
//       nastiness_1_indel_context_insertion_40
//     (sibling detailed_plan subtest #40; cv_01_nast1 in input.vcf).
//   Storage-shape independent; same vepyr allele-trimmer code path.
//   Status: GREEN under v2.
//
// -----------------------------------------------------------------------------
// SUBTEST #38 (L460-473): NASTINESS 3 — chr21:8987004 TAT→TAGCGT
//   Perl: shared prefix+suffix trim (single-ALT, so engine blocker #1
//   doesn't apply).
//   Cross-reference (canonical Rust assertion):
//     tests/port_annotation_source_cache_variation_e2e.rs::
//       nastiness_3_shared_prefix_suffix_trim_42
//     (sibling detailed_plan subtest #42; cv_03_nast3 in input.vcf).
//   Status: GREEN under v2.
//
// =============================================================================
// Architectural-no-analogue comment stubs (7 by-design vepyr gaps)
// =============================================================================
//
// Each row below documents a Perl subtest that probes machinery vepyr has
// no equivalent COMPONENT for (not "different format name" — see v2-paradigm.md
// §Rule 1 "parquet ≠ storable" trap; these are genuinely missing-by-design).
//
// -----------------------------------------------------------------------------
// SUBTEST #4 (L45): ok($c, 'new is defined')
//   Perl: constructor returns blessed Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix ref.
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     Per-source `Cache::VariationTabix` class. Variation access in vepyr
//     flows through UDTF args + `EnsemblCacheTableProvider::for_entity("variation", …)`
//     (`bio-format-ensembl-cache/src/table_provider.rs:454-456`). There is
//     no Perl-style adaptor-factory dispatch between per-region and tabix
//     backends because vepyr has a single parquet shape.
//
// -----------------------------------------------------------------------------
// SUBTEST #5 (L51): get_dump_file_name(1, '1-100') == <dir>/1/all_vars.gz
//   Perl: tabix backend uses ONE file per chrom (`all_vars.gz`).
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     The LITERAL filename "all_vars.gz" within the per-chrom directory.
//     Vepyr uses ONE parquet per chrom (`variation/<chrom>.parquet`). The
//     path-resolver CONCEPT survives in vepyr (see Cache_Variation subtest
//     #18 unit-port against `partitioned_cache.rs::context_path`), but the
//     specific filename `all_vars.gz` is a Perl-storage-shape detail with
//     no observable analogue.
//
// -----------------------------------------------------------------------------
// SUBTEST #7 (L55): delimiter == "\t"
//   Perl: tabix shape is tab-delimited (per-region shape is space-delimited).
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     Text-shape "delimiter" concept. Parquet columns are binary-typed; the
//     Perl `delimiter` accessor exists to switch the text parser between
//     per-region (space) and tabix (tab) text shapes. Vepyr's parquet reader
//     bypasses both layers.
//
// -----------------------------------------------------------------------------
// SUBTEST #41 (L525-527): _get_tabix_obj(21) ref Bio::DB::HTS::Tabix
//   Perl: tabix-PM-module returns a tabix handle of concrete adapter-class
//   `Bio::DB::HTS::Tabix`.
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     Per-chrom tabix-handle accessor. Vepyr opens parquets via DataFusion's
//     `ParquetFormat`/`TableProvider::scan`, which returns an `ExecutionPlan`,
//     not a handle. The Perl assertion probes a concrete adapter-class
//     identity that has no vepyr analogue.
//
// -----------------------------------------------------------------------------
// SUBTEST #42 (L528): _get_tabix_obj(21) eq _get_tabix_obj(21) cache identity
//   Perl: second call returns same handle (per-chrom handle caching).
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     Public per-chrom-handle cache accessor. DataFusion's session caches
//     parquet metadata implicitly via `SessionContext`, but there is no
//     public per-chrom-handle accessor to assert ref-identity on.
//
// -----------------------------------------------------------------------------
// SUBTEST #43 (L530-531): delete $c->{_tabix_obj_cache} → new handle
//   Perl: manual cache clear forces new handle on next access.
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     Public per-chrom-handle eviction. No public per-chrom-handle cache
//     exists in vepyr; the related session-level eviction is covered by the
//     `SessionCache::clear()` future-work entry but addresses a different
//     semantic (session-wide, not per-chrom).
//
// -----------------------------------------------------------------------------
// SUBTEST #47 (L575-584): toggle CAN_USE_TABIX_PM=0 then annotate_InputBuffer
//   Perl: tabix-CL fallback when PM-module is disabled at runtime.
//   Classification: architectural-no-analogue.
//   Missing-by-design vepyr component:
//     PM-vs-CL backend toggle. The dispatch decision is Perl-ecosystem
//     machinery (which Perl modules are installed at runtime); vepyr has
//     neither tabix path. The annotation behavior in the success case is
//     identical to #13 and tested there.

// =============================================================================
// Blocked-future-work commented test stubs (6 consolidated clusters
// covering 22 subtests; all reference EXISTING future-work entries)
// =============================================================================

// -----------------------------------------------------------------------------
// CLUSTER A — Chr-synonym resolution (subtests #14, #45)
//
// Blocked on TWO existing future-work entries:
//   * `ChrSynonyms` loader (future-work-vepyr.md:337)
//   * `PreparedContext chromosome-synonyms alias table` (future-work-vepyr.md:178)
//
// Vepyr's `normalize_chrom` strips `chr` prefix only; cannot resolve
// `NC_000021.9` → `21`. When both APIs land, this becomes an integration-port
// test that loads a synonyms file, annotates input.vcf with
// `chrom='NC_000021.9'`, and expects CSQ Existing_variation populated for
// rs142513484.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn variation_tabix_chr_synonym_resolution() {
//     // PRECONDITION: ChrSynonyms loader + PreparedContext::set_chrom_synonyms
//     // public APIs exist.
//     //
//     // let synonyms = ChrSynonyms::from_file("chr_synonyms.txt").unwrap();
//     // let mut ctx = PreparedContext::new(&cache_path).await.unwrap();
//     // ctx.set_chrom_synonyms(synonyms);
//     // // input.vcf row with CHROM = "NC_000021.9" POS = 25585733
//     // let out = annotate_via(&ctx, "tests/data/port/cache_variation/input_synonym.vcf").await;
//     // assert_eq!(out.csq(0).get("Existing_variation"), Some("rs142513484"));
// }

// -----------------------------------------------------------------------------
// CLUSTER B — Count annotated (subtests #17, #48)
//
// Blocked on EXISTING future-work entry: `count_features_in_region`
// (RegFeat entry, line ~205 of future-work-vepyr.md).
//
// Perl assertion is "count annotated == 132" on test_vcf via tabix backend
// (same v84-specific count as the per-region shape; v115 number will differ).
// When `count_features_in_region` lands, this becomes an integration-port
// test asserting the count on the v115 cache (number to be hand-coded once
// the API is live).
//
// #[tokio::test(flavor = "multi_thread")]
// async fn variation_tabix_count_annotated() {
//     // PRECONDITION: count_features_in_region public API exists.
//     //
//     // let n = count_features_in_region(&cache, "21", 0, u64::MAX).await.unwrap();
//     // // verified via VEP 115 at commit time:
//     // //   docker run --rm -v $(pwd)/_cache115:/cache:ro …
//     // assert_eq!(n, /* v115-specific count, hand-coded */ );
// }

// -----------------------------------------------------------------------------
// CLUSTER C — check_frequency_filter validation (subtests #19, #20, #21, #22)
//
// Blocked on TWO existing future-work entries:
//   * `AnnotateVcfConfig::filter_common + freq_freq` (future-work-vepyr.md:718)
//   * `AnnotateVcfConfig::resolve — invalid-value validation` (future-work-vepyr.md:52)
//
// Perl asserts the freq-filter config: default freq_pop accepted; 'foo'
// rejected; 'gnomADe_AFR_AF' (with AF suffix) rejected; 'ExAC_AFR' (deprecated
// population) rejected. When `filter_common` + the resolve-time validation
// pattern land, these collapse to a set of unit-port assertions on
// `AnnotateVcfConfig::resolve`.
//
// #[test]
// fn variation_tabix_check_frequency_filter_validates_freq_pops() {
//     // PRECONDITION: AnnotateVcfConfig::filter_common + freq_pops + resolve-time
//     // validation public APIs exist.
//     //
//     // assert!(AnnotateVcfConfig::resolve(&[("freq_pop", "1KG_ALL")]).is_ok());
//     // assert!(AnnotateVcfConfig::resolve(&[("freq_pop", "foo")]).is_err());
//     // assert!(AnnotateVcfConfig::resolve(&[("freq_pop", "gnomADe_AFR_AF")]).is_err());
//     // assert!(AnnotateVcfConfig::resolve(&[("freq_pop", "ExAC_AFR")]).is_err());
// }

// -----------------------------------------------------------------------------
// CLUSTER D — get_frequency_data + frequency_check_buffer
//             (subtests #23, #24, #25, #26, #27, #28, #29, #30, #31, #32, #33)
//
// Blocked on EXISTING future-work entry: `AnnotateVcfConfig::filter_common + freq_freq`
// (future-work-vepyr.md:718).
//
// Eleven subtests sharing the same gap. Coverage of:
//   * #23 — _freq_check_freqs == {1KG_ALL=>{T=>0.0002}} (per-allele AF lookup).
//   * #24, #25 — _freq_check_pass / _freq_check_all_failed accumulators.
//   * #26, #27 — freq_gt_lt='lt' threshold direction flip + aggregate.
//   * #28, #29 — freq_freq=0.0001 threshold-tuning + aggregate.
//   * #30 — freq_pop='1KG_AMR' alternate-population lookup.
//   * #31 — vf->{strand}=-1 reverse-strand frequency lookup (uses
//           `compare_existing_variant_alleles` at variant_lookup_exec.rs:458;
//           the missing layer is the freq-lookup consumer of that match).
//   * #32, #33 — frequency_check_buffer exclude/include mode counts (v84-
//                specific; v115 counts will differ).
//
// When `filter_common` lands, this cluster becomes ONE multi-row stub
// translating each Perl assertion to a Rust assertion on the new public API.
// Underlying data IS in the v115 cache (`variation.rs` reads AF/AFR/AMR/etc.);
// only the lookup API is missing.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn variation_tabix_get_frequency_data_per_pop() {
//     // PRECONDITION: filter_common + freq_freq + per-allele freq-lookup APIs exist.
//     //
//     // let freqs = freq_lookup(&cache, "21", 25585733, "T", &["1KG_ALL"]).await.unwrap();
//     // assert_eq!(freqs.get("1KG_ALL").and_then(|m| m.get("T")), Some(&0.0002));
//     // // … and for the gt_lt, freq_freq, AMR-pop, reverse-strand variants.
// }

// -----------------------------------------------------------------------------
// CLUSTER E — Per-allele frequency filter, multi-ALT (subtests #34, #35, #40)
//
// DOUBLE/TRIPLE-BLOCKED on:
//   * `AnnotateVcfConfig::filter_common + freq_freq` (future-work-vepyr.md:718)
//   * Engine blocker #1 — multi-ALT CSQ per-allele expansion
//     (port-status.md §Active blockers item 1; future-work-vepyr.md:779)
//   * For #40: `VariationCache::compare_existing()` public API
//     (already filed via Cache_Variation audit)
//
// Coverage of:
//   * #34 — 'C/T/G' input + exclude → 'C/G' (drops the high-freq T allele).
//   * #35 — 'C/T/G' input + include → 'C/T' (keeps the high-freq T allele).
//   * #40 — chr21:25005812 CA→CAAA,CAAA matched-alleles-aware AAA freq
//           lookup (input AAA after CA-CAAA trim, not raw input ALT).
//
// Without engine blocker #1 there is no way to filter individual ALTs in a
// pipe-joined-multi-ALT row; without `compare_existing` there is no public
// API for matched-alleles trimming on the consumer side.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn variation_tabix_per_allele_freq_filter_multi_alt() {
//     // PRECONDITION: filter_common UDTF surface + multi-ALT row expansion +
//     // VariationCache::compare_existing public APIs all exist.
//     //
//     // // Input row chr21:8987004 C C/T/G with --check_frequency --filter exclude.
//     // // Expect ALT column reduced to "C/G" (T dropped).
//     // let out = annotate_with_freq_filter(&cfg, "input_multi_alt.vcf").await;
//     // assert_eq!(out.row(0).alt(), "C/G");
// }

// -----------------------------------------------------------------------------
// CLUSTER F — Nastiness 2 & 4 multi-ALT (subtests #37, #39)
//
// UNBLOCKED 2026-05-28 by engine blocker #1 PARTIAL fix (commit
// `e0e00f4`, merged via PR #166 — per-allele CSQ expansion at
// `annotate_provider.rs:4706`).
//
// Cross-reference (canonical Rust assertions, NOW LIVE):
//   * tests/port_annotation_source_cache_variation_e2e.rs::
//       nastiness_2_multi_alt_only_second_alt_matches_41
//     (cv_02_nast2 in input.vcf; sibling detailed_plan subtest #41 —
//      promoted from blocked-future-work to e2e-port).
//   * tests/port_annotation_source_cache_variation_e2e.rs::
//       nastiness_4_multi_alt_two_matching_cache_rows_43
//     (cv_04_nast4 in input.vcf; sibling detailed_plan subtest #43 —
//      promoted from blocked-future-work to e2e-port).
//
// Concern: per PR #166 review notes, the engine fix is "partial" —
// typed-column writers + cache-hit fast path remain ALT[0]-only. The
// `Existing_variation` populated/empty contract IS covered (which is
// what these subtests assert). If nastiness 4's second-ALT contract
// fails, that is a follow-up to engine blocker #1, not a port bug.
//
// Status: GREEN under v2 (2026-05-28) for nastiness 2; nastiness 4 may
// surface DONE_WITH_CONCERNS depending on engine fix coverage.
//
// When engine blocker #1 lands, these promote simultaneously alongside
// cv_02/cv_04 in the sibling port's input.vcf (this port adds zero new
// assertions; the contract is owned by Cache_Variation).
//
// #[tokio::test(flavor = "multi_thread")]
// async fn variation_tabix_nastiness_2_and_4_multi_alt() {
//     // PRECONDITION: engine blocker #1 resolved — multi-ALT CSQ per-allele
//     // expansion produces one CSQ per ALT, not a pipe-joined string.
//     //
//     // Cross-reference: tests/port_cache_variation.rs covers the canonical
//     // assertion via cv_02 (nastiness 2) + cv_04 (nastiness 4) on the same
//     // v115 fixture. This port adds NO new Rust assertions when the blocker
//     // resolves; it merely flips the cross-reference rows from blocked-future-work
//     // to unit-port (cross-reference).
// }
