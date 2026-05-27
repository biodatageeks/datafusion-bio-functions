//! v2-paradigm port of `ensembl-vep/t/BaseVEP.t`.
//!
//! See `porting-tests/detailed_plans/BaseVEP.md` (AUDITED 2026-05-27) for the
//! authoritative per-subtest classification table.
//!
//! ## Why this file has 0 active `#[test]` blocks
//!
//! `Bio::EnsEMBL::VEP::BaseVEP` is the abstract base class for every
//! `Bio::EnsEMBL::VEP::*` object. Its 924 lines centralise Perl-ecosystem
//! plumbing — `AUTOLOAD` shortcuts, dynamic `param(key)` dispatch, the
//! Perl-API adaptor factory, lazy-init Stats blessed-hashes, fake-`Slice`
//! placeholders for unknown chromosomes, the `warning_file` sidecar with
//! `WARNING: $msg\n` prefix, the `skipped_variants` accumulator, and the
//! DB-mode multi-species adaptor fallback.
//!
//! Vepyr replaces this whole plumbing layer with:
//!
//! - typed Rust structs (`AnnotateVcfConfig` in `vcf_sink.rs:29` for config;
//!   no AUTOLOAD; no dynamic field-name dispatch);
//! - DataFusion table providers (`EnsemblCacheTableProvider` registered on a
//!   `SessionContext`; no Perl-style adaptor factory);
//! - the `log` crate (no `warning_file` sidecar; no `WARNING:` prefix);
//! - single-species cache-only operation (`homo_sapiens` + `GRCh38`; no
//!   DB-mode, no multi-species fallback, no fake-Slice placeholder).
//!
//! Coverage parity is therefore **0 / 47 = 0%** of substantive subtests
//! classified as `unit-port` / `integration-port` / `e2e-port`. The 47
//! substantive subtests split as:
//!
//! - **38 architectural-no-analogue** rows — each documented inline below
//!   naming the specific Perl construct that has no vepyr analogue.
//! - **9 blocked-future-work** rows — each present as a commented-out
//!   `#[test]` block pointing at a named entry in
//!   `porting-tests/future-work-vepyr.md`.
//!
//! See `detailed_plans/BaseVEP.md` §"Coverage parity" for the full
//! justification of the 0% target. The 9 blocked-future-work rows surface
//! real public-API gaps (`BaseContext::fasta_reader`, `chr_lengths`,
//! `ChrSynonyms`, `fasta_dir_autodetect`); promoting them to `unit-port`
//! requires a small refactor extracting the FASTA reader currently opened
//! inline at `annotate_provider.rs:4424` and `:8060` (type alias at
//! `annotate_provider.rs:7795`). That refactor is intentionally out of
//! scope for this port.

// ──────────────────────────────────────────────────────────────────────────
// Architectural-no-analogue rows (38)
//
// Each comment block names the specific Perl construct that has no vepyr
// analogue. "Different format" alone is NOT a justification (per
// per-subtest-classification.md case 6); each entry below points at a
// missing-by-design vepyr COMPONENT.
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #2 (L34-35 of BaseVEP.t): `BaseVEP->new()` defined with no args.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: there is no `BaseVEP` Rust struct.
// Vepyr's components are flat functions or single-purpose structs; no abstract
// base class wraps Config + FASTA + chr-synonyms. The Perl `bless {}, $class`
// constructor has no runtime equivalent.

// SUBTEST #3 (L37): `ref(bv) == 'Bio::EnsEMBL::VEP::BaseVEP'`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #2 — no `BaseVEP` class shape to
// introspect. Rust types are checked at compile time, not via `ref($obj)`.

// SUBTEST #4 (L39): `new('not a hashref')` throws.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: Perl's runtime arg-type-validation
// (`assert_ref($_[0], 'HASH')`) is compile-time in Rust. There is no analogous
// runtime "wrong-type-arg" assertion to test.

// SUBTEST #5 (L41): `new({config => {}})` throws when config is wrong type.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: Perl's `assert_ref($config,
// 'Bio::EnsEMBL::VEP::Config')` is compile-time-enforced in Rust (the
// `AnnotateVcfConfig` field types are statically checked). No runtime
// type-assertion path exists.

// SUBTEST #7 (L51-52): `Config->new` instantiates a Config object.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: covered conceptually by
// `detailed_plans/Config.md` row 2 (`AnnotateVcfConfig::default()`); BaseVEP's
// blessed-hash `_config` slot has no Rust analogue (no `bv->config` accessor;
// `AnnotateVcfConfig` IS the config).

// SUBTEST #8 (L54-55): `bv = BaseVEP->new({config => $cfg})` — construct with
// config.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no `BaseVEP` class; no `_config` slot
// pattern. `AnnotateVcfConfig` is passed directly to engine functions, not
// composed into a base-class object.

// SUBTEST #9 (L57): `bv->config` returns ref to the Config.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: no `config()` accessor because there is
// no wrapping object to accessor-fy. `AnnotateVcfConfig` is the value, not a
// field on something else.

// SUBTEST #10 (L59): `bv->param('species')` get.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: dynamic `param(key)` getter does not
// exist. `AnnotateVcfConfig` field access is direct typed Rust (e.g.,
// `cfg.reference_fasta_path`); there is no string-keyed dispatch table.

// SUBTEST #11 (L60): `bv->param('species','mouse')` set.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: dynamic `param(key, val)` setter. The
// `AnnotateVcfConfig` struct is constructed once (via builder/`Default`); no
// in-place mutation API exists.

// SUBTEST #12 (L61): `bv->param()` with no key throws.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #10 — no `param()` API to call
// without args.

// SUBTEST #13 (L64): `bv->species` get.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `AnnotateVcfConfig` has no `species`
// field. Vepyr is single-species (`homo_sapiens` / `GRCh38`) per
// `vepyr/CLAUDE.md` constraints; the multi-species accessor pattern has no
// equivalent.

// SUBTEST #14 (L65): `bv->species('human')` set.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #13 — species is not a tunable
// field in vepyr's data model.

// SUBTEST #15 (L68-78): `bv->stats == Bio::EnsEMBL::VEP::Stats blessed`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `Bio::EnsEMBL::VEP::Stats` lazy-init
// blessed-hash accessor. Vepyr has no per-run statistics writer at all —
// see `detailed_plans/Stats.md` (entire Stats port is blocked-future-work
// against engine blocker #4 in `port-status.md`).

// SUBTEST #16 (L83): `bv->get_adaptor('variation', 'VariationFeature')` in
// offline mode returns a fake adaptor.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: Perl-API adaptor factory. Vepyr cache
// readers are `EnsemblCacheTableProvider`s registered directly on a
// DataFusion `SessionContext`; the lookup-by-`(group, type)` string-keyed
// factory has no analogue, and the "fabricate a fake VariationFeatureAdaptor
// in offline mode" pattern has no analogue.

// SUBTEST #17 (L86): `bv->get_adaptor()` with no group throws.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #16 — no `get_adaptor` API to
// test the group-required error on.

// SUBTEST #18 (L87): `bv->get_adaptor('core')` with no type throws.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #16 — no type-required error to
// test.

// SUBTEST #19 (L95-96): `add_shortcuts('test')` copies the Config value into
// `$self->{test}`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: AUTOLOAD-style accelerator. Rust field
// access is direct; there is no `$self->{key}` hash-cache to populate, and
// no AUTOLOAD method-dispatch fallback to short-circuit.

// SUBTEST #20 (L103-104): `add_shortcuts(['test'])` accepts arrayref form.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #19 — AUTOLOAD shortcut
// machinery is Perl-runtime-specific.

// SUBTEST #21 (L106): `add_shortcuts('test')` second time throws (refuses
// overwrite).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #19. The "refuses overwrite"
// invariant lives in the shortcut-cache pattern; no cache → no invariant.

// SUBTEST #27 (L131-155): `bv->get_slice('1')` returns a fake-Slice with
// `seq_region_length=1` for unknown chromosomes.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: fake-`Slice` placeholder fabrication.
// Perl fabricates a `Bio::EnsEMBL::Slice` object with `is_fake => 1,
// seq_region_length => 1` so downstream `defined($slice)` checks pass even
// for unmapped chromosomes. Vepyr returns `None` (no placeholder object);
// downstream Rust callers use `Option<...>` and `match`/`?` instead of
// "defined-but-fake" sentinels. This is a deliberate semantic divergence,
// not a missing API.

// SUBTEST #33 (L194-195): `status_msg('Hello')` prints to STDOUT (skipped
// when `--quiet`).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: per-instance `status_msg` STDOUT
// writer. Vepyr uses `log::info!` macros routed through `env_logger`;
// testing STDOUT byte output for an info-level log is brittle and
// orthogonal to CSQ correctness.

// SUBTEST #34 (L217): `bv->warning_fh == FileHandle`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `warning_fh` per-instance file-handle
// accessor backed by a `warning_file` config option. Vepyr's `log` macros
// target an `env_logger` subscriber configured by the caller; there is no
// per-config "warning file handle" public API.

// SUBTEST #35 (L219): `warning_msg('test warning message')` returns true.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `warning_msg` that prefixes
// `"WARNING: "` and writes to the per-instance warning file. Vepyr emits
// `log::warn!` macros; there is no formatted-prefix-and-file convention.

// SUBTEST #36 (L222-227): warning file content contains `'WARNING: test
// warning message'`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `warning_file` sidecar with newline-
// delimited `WARNING: $msg\n` entries. Vepyr's log output is structured
// (`env_logger` line format, JSON if configured), not a flat sidecar.

// SUBTEST #37 (L237): `warning_msg('Hello')` with `warning_file => 'STDERR'`
// writes to STDERR.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `warning_file => 'STDERR'` config-
// driven re-routing of warnings to STDERR. `env_logger` defaults to STDERR
// already; testing byte-level output of a `log::warn!` call is brittle.

// SUBTEST #38 (L241): `get_database_assembly == undef` when DB not
// configured.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `get_database_assembly`. DB mode is
// out of scope by architectural decision (offline-cache only). The function
// is `undef`-only when DB is unavailable, which is vepyr's permanent state;
// there is no DB-mode-aware accessor to test.

// SUBTEST #39 (L252-253): `skipped_variant_msg("Invalid variant", 3)`
// records into `skipped_variants` and emits a warning.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `skipped_variants` accumulator
// (`arrayref` of `{description, line_number, line}` records) plus the
// `*_warnings.txt` sidecar concept. Vepyr's "skip invalid variant" path
// emits a single `log::warn!` and continues; there is no public accumulator
// and no sidecar product decision behind it.

// SUBTEST #40 (L256): `skipped[0]->{description} =~ /Invalid variant/`
// preserved.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #39 — no accumulator entries
// to introspect.

// SUBTEST #41 (L257): `skipped[0]->{line_number} == 3` preserved.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #39.

// SUBTEST #42 (L258): `skipped[0]->{line} == undef` default.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #39.

// SUBTEST #43 (L291, SKIP'd unless DB configured): `bv->registry` after
// DB-mode setup returns a loaded registry.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: `Bio::EnsEMBL::Registry` loading
// mechanism. DB mode is out of scope; the whole block is SKIP'd when DB
// unavailable in Perl too.

// SUBTEST #44 (L293): `get_adaptor('core','slice')` in DB-mode returns a
// `SliceAdaptor`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: Perl-API `SliceAdaptor` (DB-backed
// `Bio::EnsEMBL::DBSQL::SliceAdaptor`). DB mode out of scope; no adaptor
// factory in vepyr.

// SUBTEST #45 (L305): `bv->registry` populated from a registry file path.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: registry-file loader. DB mode out of
// scope.

// SUBTEST #46 (L307): `get_adaptor('core','slice')` from registry-file
// variant.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #44 — DB-mode registry-file
// path.

// SUBTEST #47 (L309): `get_slice('21')` in DB mode returns a live Ensembl
// `Slice`.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: DB-backed `Slice` fetcher. Vepyr
// reads reference sequence inline via `noodles_fasta::IndexedReader::query`
// at HGVS / shift-3' callsites; no `Slice` abstraction.

// SUBTEST #48 (L311): `get_database_assembly == 'GRCh38'` in DB mode.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: DB-derived assembly accessor. DB mode
// out of scope; assembly is configured on `AnnotateVcfConfig` at
// construction.

// SUBTEST #49 (L325-331): `get_adaptor('variation','VariationFeature')` for
// a species without a variation DB falls back to a fake VF adaptor.
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: multi-species adaptor-fallback policy.
// Vepyr is single-species; the policy "no variation DB for this species →
// fabricate a fake VF adaptor so downstream calls no-op" has no analogue.

// SUBTEST #50 (L333): `get_adaptor('variation','Variation')` for a species
// without a variation DB returns `undef` (no fake fallback for non-VF
// types).
// Classification: architectural-no-analogue.
// Missing-by-design vepyr component: same as #49 — multi-species + adaptor
// factory both out of scope.

// ──────────────────────────────────────────────────────────────────────────
// Blocked-future-work rows (9)
//
// Each commented `#[test]` block names the missing vepyr public API and
// points at the corresponding entry in `future-work-vepyr.md`. When the
// API lands, the corresponding stub is uncommented, the assertion is
// filled in, and the coverage parity counter in
// `detailed_plans/BaseVEP.md` is updated.
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #22 (L108): `bv->fasta_db == undef` when no FASTA configured.
// Blocked: vepyr's FASTA reader is opened inline inside
// `annotate_provider.rs:4424` and `:8060` (type alias at line 7795). There
// is no public `fasta_reader(cfg) -> Option<...>` accessor returning
// `None` when `cfg.reference_fasta_path.is_none()`.
// See `future-work-vepyr.md`::`BaseContext::fasta_reader / get_slice —
// public FASTA accessors`.
//
// #[test]
// fn fasta_reader_returns_none_when_no_fasta_configured() {
//     let cfg = AnnotateVcfConfig::default(); // reference_fasta_path = None
//     let reader = base_context::fasta_reader(&cfg);
//     assert!(reader.is_none());
// }

// SUBTEST #23 (L111): `bv->get_slice('1') == undef` when no FASTA configured.
// Blocked: vepyr's slice fetch is inline via `IndexedReader::query()`;
// there is no `get_slice(chr, start, end) -> Option<Vec<u8>>` public
// accessor that returns `None` when the reader is not configured.
// See `future-work-vepyr.md`::`BaseContext::fasta_reader / get_slice`.
//
// #[test]
// fn get_slice_returns_none_when_no_fasta_configured() {
//     let cfg = AnnotateVcfConfig::default();
//     let slice = base_context::get_slice(&cfg, "1", 1, 100);
//     assert!(slice.is_none());
// }

// SUBTEST #25 (L118-122): `bv->fasta_db == HTS::Faidx || Bio::DB::Fasta`
// after configuring a FASTA path.
// Blocked: same missing API as #22. Once the public `fasta_reader(cfg)`
// accessor exists, this test asserts `Some(_)` after setting
// `reference_fasta_path = Some(<test_fasta>)`.
// See `future-work-vepyr.md`::`BaseContext::fasta_reader / get_slice`.
//
// #[test]
// fn fasta_reader_returns_some_when_fasta_configured() {
//     let cfg = AnnotateVcfConfig {
//         reference_fasta_path: Some(test_fasta_path()),
//         ..Default::default()
//     };
//     let reader = base_context::fasta_reader(&cfg);
//     assert!(reader.is_some());
// }

// SUBTEST #26 (L124-129): `bv->chr_lengths == {21 => 46709983, MT => 16569}`
// FAI-derived chromosome-length map.
// Blocked: vepyr does not surface chromosome lengths anywhere — the FAI is
// consumed only by noodles internally. Adding a public
// `chr_lengths(fasta_index_path) -> Result<HashMap<String, u64>>` wrapper
// over `noodles_fasta::fai::Index` is a small extraction.
// See `future-work-vepyr.md`::`BaseContext::chr_lengths — FAI-derived
// chromosome length map`.
//
// #[test]
// fn chr_lengths_from_fai() {
//     let fai_path = test_fasta_fai_path();
//     let lengths = base_context::chr_lengths(&fai_path).unwrap();
//     assert_eq!(lengths.get("21"), Some(&46_709_983u64));
//     assert_eq!(lengths.get("MT"), Some(&16_569u64));
// }

// SUBTEST #28 (L157): `bv->chromosome_synonyms == {}` initial.
// Blocked: vepyr has no `ChrSynonyms` struct. Adding `ChrSynonyms::default()`
// plus `from_file(path)` plus `synonyms_of(chr)` plus `is_empty()` is a
// small additive feature.
// See `future-work-vepyr.md`::`ChrSynonyms — chr-synonym loader and map`.
//
// #[test]
// fn chr_synonyms_default_is_empty() {
//     let syns = ChrSynonyms::default();
//     assert!(syns.is_empty());
// }

// SUBTEST #29 (L158-159): `chromosome_synonyms(file)` loads synonyms,
// returns the resulting hashref.
// Blocked: same missing API as #28 — no `ChrSynonyms::from_file` loader
// for the 2-column TSV `canonical\tsynonym` format.
// See `future-work-vepyr.md`::`ChrSynonyms — chr-synonym loader and map`.
//
// #[test]
// fn chr_synonyms_load_from_file() {
//     let path = test_chr_synonyms_path();
//     let syns = ChrSynonyms::from_file(&path).unwrap();
//     assert!(!syns.is_empty());
// }

// SUBTEST #30 (L160-167): loaded synonym map equals
// `{21 => {CM000683.2 => 1, NC_000021.9 => 1, chr21 => 1}}`.
// Blocked: same missing API as #28 — the `synonyms_of(canonical)` getter
// must surface the membership set for canonical chrom `21`.
// See `future-work-vepyr.md`::`ChrSynonyms — chr-synonym loader and map`.
//
// #[test]
// fn chr_synonyms_for_21_contains_canonical_aliases() {
//     let path = test_chr_synonyms_path();
//     let syns = ChrSynonyms::from_file(&path).unwrap();
//     let names = syns.synonyms_of("21").expect("21 in map");
//     assert!(names.contains("CM000683.2"));
//     assert!(names.contains("NC_000021.9"));
//     assert!(names.contains("chr21"));
// }

// SUBTEST #31 (L179): `bv->fasta_db` after `fasta_dir` auto-detect opens
// the single `.fa` inside the directory.
// Blocked: vepyr requires explicit `reference_fasta=<path>` Python kwarg;
// there is no `fasta_dir_autodetect(dir) -> Result<PathBuf>` helper that
// walks a directory and returns the unique `.fa`/`.fa.gz` match.
// See `future-work-vepyr.md`::`fasta_dir_autodetect — directory walk for
// a single `.fa``.
//
// #[test]
// fn fasta_dir_autodetect_finds_unique_fa() {
//     let dir = test_fasta_dir_path();
//     let fasta_path = base_context::fasta_dir_autodetect(&dir).unwrap();
//     assert!(fasta_path.ends_with(".fa") || fasta_path.ends_with(".fa.gz"));
// }

// SUBTEST #32 (L180): `bv->param('fasta')` returns the auto-detected path
// (side-effect on Config).
// Blocked: same missing API as #31. The Perl semantics also imply that the
// helper mutates the Config; vepyr's Config is immutable post-construction,
// so the Rust shape would return the detected path for the caller to set
// (not mutate in place).
// See `future-work-vepyr.md`::`fasta_dir_autodetect — directory walk for
// a single `.fa``.
//
// #[test]
// fn fasta_dir_autodetect_returns_path_for_caller_to_set() {
//     let dir = test_fasta_dir_path();
//     let fasta_path = base_context::fasta_dir_autodetect(&dir).unwrap();
//     let cfg = AnnotateVcfConfig {
//         reference_fasta_path: Some(fasta_path.clone()),
//         ..Default::default()
//     };
//     assert_eq!(cfg.reference_fasta_path.as_deref(), Some(fasta_path.as_path()));
// }

// ──────────────────────────────────────────────────────────────────────────
// End of port_base_vep.rs.
//
// Total substantive subtests covered: 47 (3 use_ok boilerplate rows omitted
// from coverage parity per per-subtest-classification.md case 1).
//   - 38 architectural-no-analogue (rows 2-5, 7-21, 27, 33-50).
//   - 9 blocked-future-work (rows 22, 23, 25, 26, 28, 29, 30, 31, 32).
//   - 0 unit-port / integration-port / e2e-port.
//
// Coverage parity: 0 / 47 = 0% — justified in
// `detailed_plans/BaseVEP.md` §"Coverage parity".
// ──────────────────────────────────────────────────────────────────────────
