//! v2-paradigm port of `ensembl-vep/t/AnnotationSourceAdaptor.t`.
//!
//! Detailed plan: `porting-tests/detailed_plans/AnnotationSourceAdaptor.md`
//! (AUDITED 2026-05-27, Phase-D Axis B addition same date).
//! TDD task plan: `porting-tests/plans/2026-05-28-port-annotation-source-adaptor.md`.
//!
//! ## Why this file is mostly inline documentation
//!
//! `Bio::EnsEMBL::VEP::AnnotationSourceAdaptor` is a Perl factory class that
//! returns an ordered list of `AnnotationSource` blessed-hashes — one
//! `Cache::Transcript`, optionally a `Cache::Variation`, optionally a
//! `Cache::RegFeat`, plus any number of `--custom` file sources and the
//! whole DB-mode parallel hierarchy (`Database::Transcript`,
//! `Database::Variation`, `Database::RegFeat`, `Database::StructuralVariation`).
//!
//! Vepyr replaces the entire factory layer with **stateless typed
//! constructors**:
//!
//! - `EnsemblCacheTableProvider::for_entity(kind, options)` (and the per-kind
//!   `VariationTableProvider::new` / `TranscriptTableProvider::new` / ...)
//!   replace `get_all_from_cache()`. Each caller picks the entity kinds it
//!   wants; there is no list-of-sources object and no first-source-wins
//!   ordering.
//! - There is no `get_all_custom()` analogue. `--custom` is out of vepyr
//!   scope (no `AnnotationSource::File::VCF` / `BED` / `BigWig` / `GFF` /
//!   `GTF` reader).
//! - There is no `get_all_from_database()` analogue. vepyr is
//!   offline-cache-only by architectural decision.
//!
//! ## Coverage parity
//!
//! - Perl denominator (29 raw subtests; 26 substantive after dropping 3
//!   boilerplate `use_ok` / `done_testing` rows): **1 / 26 = 3.8 %**
//!   substantive subtests classified as `unit-port` / `integration-port` /
//!   `e2e-port`.
//! - Vepyr-portable denominator (1 in-scope row + 1 Axis B addition): **2 / 2 = 100 %**.
//! - The 25 `architectural-no-analogue` rows are documented inline below
//!   (one comment block each). 0 `blocked-future-work` rows: every
//!   out-of-scope row here is a deliberate architectural exclusion
//!   (`--custom` and DB-mode), not a missing vepyr public API.
//!
//! See `detailed_plans/AnnotationSourceAdaptor.md` §"Coverage parity"
//! and §"Architectural-no-analogue justifications" for the substantive
//! per-row rationale.
//!
//! ## v2 paradigm anchors
//!
//! `~/.claude/skills/port-to-vepyr/references/v2-paradigm.md`:
//!
//! - **Sztywno 1:1** — every Perl subtest gets a Rust analogue (here:
//!   passing test or inline `// SUBTEST #N …` comment naming the
//!   missing-by-design vepyr component).
//! - **Standalone tests** — no docker, no `golden.vcf`, no
//!   `port_common::run_and_compare_csq`. The two unit-port tests SKIP
//!   gracefully when the v115 fixture is absent (LFS pointer or missing
//!   path); SKIP shape mirrors `tests/port_annotation_source_cache_regfeat.rs:107-124`.

// ──────────────────────────────────────────────────────────────────────────
// Architectural-no-analogue rows (25)
//
// One comment block per row in the detailed_plan’s per-subtest table.
// Each block names a specific missing-by-design vepyr component. “Different
// format” / “Perl-internal” alone is NOT a justification (per
// per-subtest-classification.md case 6).
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #2 (L32-33 of AnnotationSourceAdaptor.t): `AnnotationSourceAdaptor->new()`
// constructs an object with no args.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: there is no `AnnotationSourceAdaptor` Rust
// struct. vepyr's annotation-source factory pattern is replaced by stateless
// `EnsemblCacheTableProvider::for_entity(kind, options)` construction sites
// invoked directly by callers (`AnnotateProvider` and friends). The Perl
// zero-arg constructor's "returns defined blessed ref" semantic has no Rust
// equivalent (Rust types are statically typed, not bless-tagged).

// SUBTEST #3 (L35): `ref($asa) == 'Bio::EnsEMBL::VEP::AnnotationSourceAdaptor'`
// — class-identity check.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #2 — no Rust class shape to
// introspect at runtime. Rust types are compile-time, not via `ref($obj)`
// string comparison.

// SUBTEST #5 (L42-43): `Config->new(...)` returns a defined config.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: covered conceptually by
// `detailed_plans/Config.md` row 3 — `AnnotateVcfConfig` is a clap-derived
// struct with no separate "construct & validate" object-instantiation step.
// There is no `Config->new`-returns-defined assertion to mirror.

// SUBTEST #6 (L45): `AnnotationSourceAdaptor->new({config => $cfg})` constructs
// with config.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #2 — no `AnnotationSourceAdaptor`
// class; no per-instance `_config` slot. `AnnotateVcfConfig` is passed directly
// to engine functions, not stored on a factory adaptor object.

// SUBTEST #8 (L90): `$asa->get_all()` mirrors `$asa->get_all_from_cache()`
// when offline (the default).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: there is no `get_all` aggregator in vepyr
// (the list-of-sources concept does not exist). The offline-vs-non-offline
// branching that `get_all` collapses in Perl is moot in vepyr — vepyr is
// offline-cache-only by design.

// SUBTEST #9 (L92-94): with `check_existing=1`, `$asa->get_all()->[0]` is a
// `Cache::Variation` blessed hash — the variation source comes first.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: there is no ordered `Vec<AnnotationSource>`
// in vepyr's annotation pipeline. Perl `Runner.pm` iterates sources in order,
// and the first match "wins"; vepyr's `AnnotateProvider` has explicit per-field
// merge logic in `transcript_consequence.rs` + `annotate_provider.rs`. The
// observable CSQ-side equivalent (`Existing_variation` populated under
// `--check_existing`) is owned by `AnnotationSource_Cache_Variation.md`'s
// integration-port rows — duplicating it here would split that contract.

// SUBTEST #10 (L96-97): `custom='format=vcf,...'` (no `file=`) throws
// `/No file was added/`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `--custom` is out of vepyr scope. There
// is no `AnnotationSource::File::VCF` / `BED` / `BigWig` / `GFF` / `GTF`
// reader, and therefore no validation surface to mirror.

// SUBTEST #11 (L99-100): `custom='file=...'` (no `format=`) throws
// `/No format specified/`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — `--custom` validation has
// no vepyr-side analogue.

// SUBTEST #12 (L102-103): `custom='...,format=foo,...'` throws
// `/Unknown or unsupported format/`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — no custom-format dispatch
// table on the vepyr side.

// SUBTEST #13 (L105-106): `custom='...,random=2,...'` throws
// `/options are not supported: random/`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — no allowed-option allowlist
// for `--custom` on the vepyr side.

// SUBTEST #14 (L108-111): `no_remote=1` + `file=http://...` throws
// `/remote data files disabled/`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: vepyr never reaches remote URLs — it is
// offline-cache-only by design. There is no `--no_remote` flag and no
// remote-URL fetch path to gate.

// SUBTEST #16 (L120-152, SKIP-gated on Bio::DB::HTS::Tabix): well-formed
// `file=<vcf>,format=vcf,short_name=test,type=exact` builds an `AnnotationSource::File::VCF`
// blessed hash with the expected slot shape.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — `--custom` is out of vepyr
// scope; there is no `AnnotationSource::File::VCF` constructor to mirror.

// SUBTEST #17 (L154-188, SKIP): the same spec with `type=` omitted defaults
// `report_coords` to undef / `type` to "overlap".
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — no custom-source `type`
// option exists on the vepyr side.

// SUBTEST #18 (L190-225, SKIP): the same spec with `coords=1` sets
// `report_coords=1`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — no `report_coords` option
// because no custom-source path exists.

// SUBTEST #19 (L227-262, SKIP): the same spec with `fields=FOO%BAR` parses
// `fields => ['FOO', 'BAR']`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — no custom-source `fields`
// list option exists.

// SUBTEST #20 (L264-365, SKIP): the `###CHR###` substitution token in a
// file-path template expands into one source per chromosome.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — vepyr's per-chromosome
// cache files (`discover_*_files` in `bio-format-ensembl-cache::discovery`)
// have a *cache-internal* per-chrom layout, but there is no user-templated
// `###CHR###` substitution path in `--custom` mode (because `--custom` is
// out of scope).

// SUBTEST #21 (L370-395, SKIP unless DB available): DB-mode
// `get_all_from_database()` returns a `Database::Transcript` blessed hash.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: MySQL DB-mode is out of vepyr scope.
// vepyr is offline-cache-only; there is no `AnnotationSource::Database::Transcript`
// analogue. Cache-mode (row 7) covers the equivalent observable contract.

// SUBTEST #22 (L397-422, SKIP): DB-mode `get_all()` mirrors
// `get_all_from_database()` when `database=1`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no `get_all` aggregator in vepyr
// (see #8); additionally, no DB-mode branching to aggregate over.

// SUBTEST #23 (L424-432, SKIP): with `check_existing=1`, `get_all()->[1]`
// is `Database::Variation`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 — no
// `AnnotationSource::Database::Variation` analogue. The Perl ordering claim
// (same as #9) has no vepyr equivalent either.

// SUBTEST #24 (L434-442, SKIP): with `regulatory=1`, `get_all()->[1]` is
// `Database::RegFeat`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 — no
// `AnnotationSource::Database::RegFeat` analogue. Cache-mode regulatory
// (row 3 of `categorization.md`) covers the equivalent observable contract.

// SUBTEST #25 (L444-452, SKIP): with `check_svs=1`, `get_all()->[1]` is
// `Database::StructuralVariation`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 — no DB-mode SV analogue.
// Additionally, SV/CNV/BND classifier work is tracked as an independent
// engine blocker in `port-status.md`; this subtest cannot be mirrored even
// at the cache-mode level today.

// SUBTEST #26 (L454-466, SKIP): on the `homo_vepiens_coreonly` MultiTestDB
// (no var/reg DB), DB-mode `get_all_from_database()` returns exactly one
// `Database::Transcript`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 — vepyr has no concept of a
// "core-only" species DB because it has no DB-mode at all.

// SUBTEST #27 (L468-476, SKIP): same `homo_vepiens_coreonly` species with
// `check_existing=1` still returns 1 source (no var DB available).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 / #26.

// SUBTEST #28 (L478-486, SKIP): same species with `regulatory=1` still
// returns 1 source (no reg DB available).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 / #26.

// SUBTEST #29 (L488-496, SKIP): same species with `check_svs=1` still
// returns 1 source (no var DB available, hence no SV either).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #21 / #26.

// ──────────────────────────────────────────────────────────────────────────
// Active unit-port `#[test]`s (rows 7 + Axis B B1)
//
// The two substantive vepyr-portable rows. Both are gated on the
// `cache-builder` Cargo feature because that feature is what pulls in the
// `datafusion-bio-format-ensembl-cache` crate (see Cargo.toml). When the
// feature is off, the tests are inert and don't reach the test runner.
//
// Both tests SKIP gracefully when the v115 fixture is missing or LFS-stubbed
// (matches `tests/port_annotation_source_cache_regfeat.rs:107-124`).
// ──────────────────────────────────────────────────────────────────────────

#[cfg(feature = "cache-builder")]
mod active {
    use std::path::{Path, PathBuf};

    use datafusion_bio_format_ensembl_cache::{
        CacheSourceType, EnsemblCacheOptions, EnsemblCacheTableProvider, EnsemblEntityKind,
        VariationTableProvider,
    };

    /// Resolve a path relative to the workspace root (mirrors the helper in
    /// `tests/port_annotation_source_cache_regfeat.rs`).
    fn workspace_path(rel: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(rel)
    }

    /// Detect a Git LFS pointer file so a partially-pulled checkout SKIPs
    /// rather than panicking on a sentinel.
    fn is_lfs_pointer(path: &Path) -> bool {
        std::fs::read_to_string(path)
            .map(|s| s.starts_with("version https://git-lfs.github.com"))
            .unwrap_or(false)
    }

    /// Resolve the v115 storable native cache root used by the Perl AnnotationSourceAdaptor
    /// equivalent path. `EnsemblCacheTableProvider::for_entity` reads raw VEP storable
    /// caches (with `info.txt` + per-chrom region-named files like
    /// `1-1000000_transcript.storable.gz`), NOT the vepyr-output parquet partition tree
    /// (which lives under `parquet/115_GRCh38_vep/<entity>/<chrom>.parquet` and is
    /// produced by `cache_builder` from the raw cache via the same provider, with the
    /// raw cache as input).
    ///
    /// The v115 fixture stores the raw cache under
    /// `vep-benchmark/data/port/_cache115/native_cache/homo_sapiens/115_GRCh38/`.
    ///
    /// Returns `None` (with an `eprintln!` SKIP reason) when the native cache or its
    /// `info.txt` is missing or LFS-stubbed.
    fn v115_storable_cache_root() -> Option<PathBuf> {
        let cache_root = workspace_path(
            "vep-benchmark/data/port/_cache115/native_cache/homo_sapiens/115_GRCh38",
        );
        let info_txt = cache_root.join("info.txt");

        if !cache_root.exists()
            || !info_txt.exists()
            || is_lfs_pointer(&info_txt)
        {
            return None;
        }

        // Sanity-check the per-chrom dir layout is intact (not LFS-stubbed).
        // The v115 fixture stores raw storable files in `21/` and `MT/`.
        let chr21_dir = cache_root.join("21");
        let mt_dir = cache_root.join("MT");
        if !chr21_dir.exists() || !mt_dir.exists() {
            return None;
        }

        Some(cache_root)
    }


    // SUBTEST #7 (L52-88 of AnnotationSourceAdaptor.t): default-config
    // `get_all_from_cache()` returns `[Cache::Transcript blessed-hash]`.
    // Classification: unit-port.
    // Vepyr analogue: `EnsemblCacheTableProvider::for_entity(
    //     EnsemblEntityKind::Transcript, options)` builds a
    // `TranscriptTableProvider` from the v115 fixture. The Perl blessed-hash
    // *shape* (`_config`, `serializer_type`, `dir`, ...) is Perl-internal;
    // vepyr's analogue surface is "the right provider was constructed with
    // the right schema" — which is what we assert.
    #[test]
    fn row7_transcript_provider_from_v115_cache() {
        let Some(cache_root) = v115_storable_cache_root() else {
            eprintln!(
                "SKIP row7_transcript_provider_from_v115_cache: v115 native_cache fixture or its info.txt missing/LFS-stubbed"
            );
            return;
        };

        let options = EnsemblCacheOptions::new(cache_root.to_string_lossy().into_owned())
            .with_cache_source_type(CacheSourceType::Ensembl);

        let provider =
            EnsemblCacheTableProvider::for_entity(EnsemblEntityKind::Transcript, options)
                .expect("transcript provider should be constructible from v115 native cache");

        // The schema contract is set by `transcript_schema()` in
        // `bio-format-ensembl-cache::schema`. Asserting the transcript-shaped
        // columns are present pins "the right provider was constructed with
        // the right schema" without coupling to the full field list (which
        // evolves with cache version). This is the vepyr-portable equivalent
        // of Perl row 7's blessed-hash shape assertion (`_config`, `dir`,
        // `serializer_type`, ...).
        //
        // Field names verified at commit time (2026-05-28) against the
        // observed schema for the v115 storable cache via this very test's
        // earlier failure message: `chrom`, `start`, `end`, `strand`,
        // `stable_id` (= transcript stable ID), `biotype`, `source`,
        // `gene_stable_id`, `gene_symbol`, `exons`, `cdna_seq`,
        // `peptide_seq`, ... `mane_select`, `appris`.
        let schema = provider.schema();
        for required in ["stable_id", "biotype", "gene_stable_id", "exons"] {
            assert!(
                schema.column_with_name(required).is_some(),
                "transcript schema must contain {required}; got fields: {:?}",
                schema
                    .fields()
                    .iter()
                    .map(|f| f.name().clone())
                    .collect::<Vec<_>>()
            );
        }
    }

    // AXIS B ROW B1: `VariationTableProvider::chromosomes()` on the v115
    // fixture returns `Some(["21", "MT"])` — symmetry with row 7's per-entity
    // construction for the Variation entity (separate file-discovery path:
    // `discover_variation_files` instead of `discover_transcript_files`).
    // Classification: unit-port.
    // Vepyr analogue: `VariationTableProvider::new` + `chromosomes()` at
    // `bio-format-ensembl-cache::table_provider:481-494`.
    #[test]
    fn row_b1_variation_provider_v115_chromosomes() {
        let Some(cache_root) = v115_storable_cache_root() else {
            eprintln!(
                "SKIP row_b1_variation_provider_v115_chromosomes: v115 native_cache fixture or its info.txt missing/LFS-stubbed"
            );
            return;
        };

        let options = EnsemblCacheOptions::new(cache_root.to_string_lossy().into_owned())
            .with_cache_source_type(CacheSourceType::Ensembl);

        let provider = VariationTableProvider::new(options)
            .expect("variation provider should be constructible from v115 native cache");

        let chroms = provider
            .chromosomes()
            .expect("variation provider must expose chromosomes from per-chr file layout");

        // Verified at commit time (2026-05-28) by:
        //   ls vep-benchmark/data/port/_cache115/native_cache/homo_sapiens/115_GRCh38/
        // → 21/, MT/, info.txt
        // The per-chrom dirs each contain `all_vars.gz` plus storable region
        // files, which is the v115 tabix-variation layout the provider's
        // `discover_variation_files` consumes.
        let mut sorted: Vec<&str> = chroms.iter().map(String::as_str).collect();
        sorted.sort();
        assert_eq!(sorted, vec!["21", "MT"]);
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Stub `#[test]` when the `cache-builder` feature is OFF.
//
// Keeps the test target non-empty under `cargo test --no-default-features`
// (or any configuration that drops `cache-builder`), and surfaces the SKIP
// reason so reviewers know the feature gate is the only thing preventing
// active coverage.
// ──────────────────────────────────────────────────────────────────────────

#[cfg(not(feature = "cache-builder"))]
#[test]
fn cache_builder_feature_required_for_active_rows() {
    eprintln!(
        "SKIP: rows 7 + B1 require the `cache-builder` feature \
         (pulls in datafusion-bio-format-ensembl-cache). \
         Run with `cargo test -p datafusion-bio-function-vep \
         --features cache-builder --test port_annotation_source_adaptor`."
    );
}
