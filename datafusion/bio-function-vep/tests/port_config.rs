//! v115 port of `ensembl-vep/t/Config.t` (v2 paradigm — sztywno 1:1).
//!
//! See `porting-tests/detailed_plans/Config.md` and the per-port plan
//! `porting-tests/plans/2026-05-27-port-config.md`.
//!
//! ## What this file covers
//!
//! Perl `Bio::EnsEMBL::VEP::Config` is a 1178-line INI/option-set/DB-mode/
//! validator/deprecation engine. Its CLI surface is `--*` flags, INI files,
//! and a hash-shaped `param(key, value)` runtime accessor. Vepyr's
//! "config" surface is fundamentally different and ~10x smaller:
//!
//! 1. `AnnotateVcfConfig` (this crate, `src/vcf_sink.rs:29-95`) — a plain
//!    Rust struct with ~30 boolean/option fields. `Default` impl sets every
//!    field to `false`/`None`/library-default. There is no internal
//!    option-set expansion, no INI loader, no validator, no plugin system,
//!    no DB mode, no `--safe`, no `--custom`.
//! 2. `AnnotationConfig` (this crate, `src/config.rs:21`) — DataFusion
//!    session-level config exposing only `cache_size_mb`, `zstd_level`,
//!    `dict_size_kb`. Not VEP-related; not exercised by Config.t.
//! 3. Python `annotate()` kwargs (`vepyr/src/vepyr/__init__.py:417-621`) —
//!    the user-facing surface. Validation (`raise ValueError(...)` for
//!    `everything`-requires-fasta, refseq/merged mutual exclusion,
//!    gencode_basic/gencode_primary mutual exclusion) lives here at lines
//!    608-621.
//!
//! Vepyr has NO clap CLI. This is a deliberate architectural decision per
//! `vepyr/CLAUDE.md` Constraints ("Product surface: Python-first library
//! API"). Do not add one to "increase parity" with Perl.
//!
//! ## Coverage parity
//!
//! - Total substantive subtests in Config.t: 37 (after stripping #1 use_ok).
//! - unit-port (substantive Rust test): 1 (subtest #2 — see `t1_*` below).
//! - architectural-no-analogue: 32 (documented inline).
//! - blocked-future-work: 4 (commented-out test stubs + entries in
//!   `porting-tests/future-work-vepyr.md` lines 36-66).
//! - Parity: 1 / 37 = 3% — well below 90% target; this is a
//!   mostly-documentation port. See detailed_plan §"Coverage parity" for
//!   the full justification.
//!
//! ## Anti-goals
//!
//! - Do not add clap. Vepyr is Python-first.
//! - Do not add an INI parser.
//! - Do not add a `param(key, value)` dynamic getter/setter.
//! - Do not add a `species` field to `AnnotateVcfConfig`.
//! - Do not add database-mode / `--custom` / `--safe` / plugin-system flags.

use datafusion_bio_format_vcf::VcfCompressionType;
use datafusion_bio_function_vep::vcf_sink::{AnnotateVcfConfig, VEP_DEFAULT_BUFFER_SIZE};

// =========================================================================
// SUBSTANTIVE PORT — subtest #2 (constructor returns object with defaults)
// =========================================================================

/// Subtest #2 of Config.t: `Bio::EnsEMBL::VEP::Config->new()` returns an
/// object; the vepyr analogue is `AnnotateVcfConfig::default()`. This test
/// pins every field's default value, so a regression in `Default` impl is
/// caught at compile-or-test time.
///
/// Hand-coded against `vcf_sink.rs:97-133` (Default impl).
#[test]
fn t1_annotate_vcf_config_default_all_off() {
    let cfg = AnnotateVcfConfig::default();

    // 24 bool fields — all default to false.
    assert!(!cfg.everything, "everything default false");
    assert!(!cfg.pick, "pick default false");
    assert!(!cfg.pick_allele, "pick_allele default false");
    assert!(!cfg.per_gene, "per_gene default false");
    assert!(!cfg.pick_allele_gene, "pick_allele_gene default false");
    assert!(!cfg.flag_pick, "flag_pick default false");
    assert!(!cfg.flag_pick_allele, "flag_pick_allele default false");
    assert!(!cfg.flag_pick_allele_gene, "flag_pick_allele_gene default false");
    assert!(!cfg.extended_probes, "extended_probes default false");
    assert!(!cfg.hgvs, "hgvs default false");
    assert!(!cfg.hgvsc, "hgvsc default false");
    assert!(!cfg.hgvsp, "hgvsp default false");
    assert!(!cfg.no_escape, "no_escape default false");
    assert!(!cfg.remove_hgvsp_version, "remove_hgvsp_version default false");
    assert!(!cfg.hgvsp_use_prediction, "hgvsp_use_prediction default false");
    assert!(!cfg.refseq, "refseq default false");
    assert!(!cfg.merged, "merged default false");
    assert!(!cfg.gencode_basic, "gencode_basic default false");
    assert!(!cfg.gencode_primary, "gencode_primary default false");
    assert!(!cfg.all_refseq, "all_refseq default false");
    assert!(!cfg.exclude_predicted, "exclude_predicted default false");
    assert!(!cfg.show_progress, "show_progress default false");

    // Option<_> fields — all default to None.
    assert!(cfg.pick_order.is_none(), "pick_order default None");
    assert!(
        cfg.reference_fasta_path.is_none(),
        "reference_fasta_path default None"
    );
    assert!(cfg.shift_hgvs.is_none(), "shift_hgvs default None");
    assert!(cfg.failed.is_none(), "failed default None");
    assert!(cfg.distance.is_none(), "distance default None");
    assert!(
        cfg.on_batch_written.is_none(),
        "on_batch_written default None"
    );

    // Library-default fields.
    assert_eq!(
        cfg.buffer_size, VEP_DEFAULT_BUFFER_SIZE,
        "buffer_size default == VEP_DEFAULT_BUFFER_SIZE (5000)"
    );
    assert!(
        matches!(cfg.compression, VcfCompressionType::Plain),
        "compression default Plain"
    );
}

// =========================================================================
// ARCHITECTURAL-NO-ANALOGUE — 32 rows
// =========================================================================
//
// The following subtests of Config.t have no Rust analogue because the
// underlying Perl mechanism is missing-by-design in vepyr. Each row names
// the Perl construct and the corresponding vepyr-side absence. These are
// NOT API gaps (which would be `blocked-future-work`); they are design
// invariants of vepyr's scope.
//
// --- Rust type system replaces runtime introspection (rows 3-7) ---
//
// #3  Perl `ref($cfg) == 'Bio::EnsEMBL::VEP::Config'` — Rust types are
//     compile-time. Runtime class-name introspection has no equivalent.
// #4  Perl `Config->new('not a hashref')` throws — Rust's type system
//     makes "wrong type" a compile error, not a runtime failure.
// #5  Perl `$cfg->param('testing', 'goodbye')` set/get — vepyr's
//     `AnnotateVcfConfig` is a plain struct accessed by field
//     (`cfg.everything`). There is no Perl-style dynamic `param(key)`
//     getter/setter; key names are compile-time, not runtime.
// #6  Perl `param('testing', 'goodbye')` returns the set value — same
//     justification as #5.
// #7  Perl `param()` with no key throws — same; no dynamic key dispatch.
//
// --- No INI parser (rows 8-16, 24-25) ---
//
// #8  Perl `read_config_from_file()` no-arg throws — vepyr has no INI
//     loader (config surface is Python kwargs + Rust struct, by design).
// #9  Perl `read_config_from_file('does_not_exist')` throws — same.
// #10 Perl `read_config_from_file(test_ini)` ok — same.
// #11 Perl INI `format == ensembl` — no INI loader.
// #12 Perl INI `dir == '/this/is/a path with spaces'` (quoted-string
//     parsing) — no INI loader.
// #13 Perl INI `terms == display` (flag-not-allowed-multiple
//     enforcement) — no INI loader.
// #14 Perl INI `fields == 'Allele,Consequence,...'` (comma-list as
//     string) — no INI loader. Vepyr emits a fixed CSQ field list
//     driven by `everything` + per-feature kwargs; no `--fields`
//     whitelist projection.
// #15 Perl INI `plugin == ['test.pm,...', 'anotherPlugin.pm']`
//     (multi-allowed appendable arrayref) — vepyr has NO plugin system;
//     `--plugin` has no vepyr-side analogue.
// #16 Perl INI re-read appends plugin (4 entries) — same.
// #24 Perl `--config test_ini` sets format — no INI loader.
// #25 Perl `--dir <path>` auto-detect `dir/vep.ini` — no auto-detect INI.
//
// --- Single-mode / single-species / cache-only by design
//     (rows 17, 30, 33) ---
//
// #17 Perl `Config->new` defaults `species == homo_sapiens` — vepyr's
//     scope is fixed `homo_sapiens` / `GRCh38` (per `vepyr/CLAUDE.md`
//     Constraints — single species and assembly). There is no `species`
//     field on `AnnotateVcfConfig`; vepyr cannot annotate mouse.
// #30 Perl `{database => 1, cache => 1}` incompatible — vepyr is
//     cache-only; `--database` flag has no analogue.
// #33 Perl `{database => 0}` "no source of gene data" throws — vepyr
//     unconditionally requires `cache_dir` at the Python boundary
//     (`annotate(vcf, cache_dir, ...)`); the "you must pick one of three
//     modes" check has no equivalent because vepyr has one mode.
//
// --- No option-set indirection / DB-mode collapse (rows 18, 19, 36-38) ---
//
// #18 Perl `regulatory == undef` for empty `cell_type` arrayref —
//     Perl-Config-internal option-set indirection; vepyr has no
//     `cell_type` flag and no option-set indirection.
// #19 Perl `{genomes => 1}` sets host/port — DB-mode option-set; vepyr
//     is offline-only, no `--genomes`.
// #36 Perl `{everything, database}` turns off `af_1kg` — DB mode out of
//     scope (the "everything + database collapse" rule doesn't apply).
// #37 Same: turns off `af_gnomad` — same.
// #38 Same: turns off `pubmed` — same.
//
// --- Out-of-scope deprecated / custom / per-sample flags
//     (rows 21, 22, 23, 26, 27, 31, 32) ---
//
// #21 Perl `{check_frequency => 1}` sets `check_existing 1` —
//     `--check_frequency` is a deprecated alias for `--af` in VEP 115.
//     Vepyr exposes `af` directly; deprecated aliases have no analogue
//     (would be out-of-scope as deprecated synonym anyway).
// #22 Perl `{gff => 'test'}` → custom annotation — `--custom` and
//     `--gff` are custom-annotation flags. Vepyr does not support
//     custom annotation sources (per `vepyr/CLAUDE.md` constraint:
//     outputs are cache-driven CSQ; no `--custom` in engine).
// #23 Perl `{ucsc_assembly, phyloP, custom}` template substitution
//     (URL templating in option-set) — same; custom annotation out of
//     scope.
// #26 Perl `individual => 'dave,barry,keith'` comma-string to list —
//     vepyr has no `--individual` / per-sample selection flag
//     (single-sample annotation is the v1 scope).
// #27 Perl `individual_zyg` list conversion — same.
// #31 Perl `{safe, most_severe, symbol}` passes with safe-on —
//     vepyr has no `--safe` mode (safe-mode bypasses conflict checks;
//     not in vepyr's data model).
// #32 Perl `{phyloP => 1}` no-custom throws — `--phyloP` is a
//     custom-annotation companion flag; not in vepyr scope.
//
// --- Boilerplate (row 1, omitted from coverage parity denominator) ---
//
// #1  Perl `use_ok Config` — cargo handles module-load at compile time;
//     this row is omitted from the coverage-parity denominator (37
//     substantive subtests, not 38).

// =========================================================================
// BLOCKED-FUTURE-WORK — 4 entries covering 5 subtests
// =========================================================================
//
// These subtests probe behaviors that vepyr COULD support with reasonable
// effort (small additions to the public API surface), but does not today.
// Each row points to its entry in `porting-tests/future-work-vepyr.md`.

// SUBTEST #20 — Perl `{af => 1}` sets `check_existing 1` (option-set fan-out).
// Blocked-future-work: vepyr has no `AnnotateVcfConfig::resolve()` method
// performing option-set fan-out. Today the fan-out lives in the Python
// `annotate()` wrapper (`vepyr/src/vepyr/__init__.py:649-674`). Surfacing
// it on the Rust side lets non-Python callers verify the rules.
// Future-work entry:
//   porting-tests/future-work-vepyr.md
//     :: "AnnotateVcfConfig::resolve — option-set fan-out"
//
// #[test]
// fn t20_af_implies_check_existing() {
//     let cfg = AnnotateVcfConfig { af: true, ..Default::default() }
//         .resolve()
//         .unwrap();
//     assert!(cfg.check_existing, "--af must imply --check_existing");
// }

// SUBTEST #28 — Perl `{convert => 1}` throws "deprecated".
// Blocked-future-work: vepyr has no `ConfigError::DeprecatedFlag` and no
// deprecated-flag list. If a user passed `--convert` through Python
// **kwargs, it would silently no-op rather than erroring.
// Future-work entry:
//   porting-tests/future-work-vepyr.md
//     :: "AnnotateVcfConfig::resolve — deprecated-flag rejection"
//
// #[test]
// fn t28_convert_flag_is_deprecated() {
//     // Hypothetical: a future `convert: bool` field on the struct.
//     // let cfg = AnnotateVcfConfig { convert: true, ..Default::default() };
//     // let err = cfg.resolve().unwrap_err();
//     // assert!(matches!(err, ConfigError::DeprecatedFlag(name) if name == "convert"));
// }

// SUBTEST #29 — Perl `{format => 'gobbledegook'}` throws "not a valid value".
// Blocked-future-work: vepyr's `pick_order: Option<String>` is unvalidated
// today. A `ConfigError::InvalidValue { field, value }` returned from
// `resolve()` would mirror Perl's behavior.
// Future-work entry:
//   porting-tests/future-work-vepyr.md
//     :: "AnnotateVcfConfig::resolve — invalid-value validation"
//
// #[test]
// fn t29_invalid_pick_order_rejected() {
//     let cfg = AnnotateVcfConfig {
//         pick_order: Some("gobbledegook".to_string()),
//         ..Default::default()
//     };
//     let err = cfg.resolve().unwrap_err();
//     assert!(matches!(err, ConfigError::InvalidValue { field, value }
//         if field == "pick_order" && value == "gobbledegook"));
// }

// SUBTESTS #34 + #35 — Perl `{output_file => STDOUT, verbose => 1}` flips
// `quiet => 1, verbose => undef` (STDOUT-implies-quiet side effect).
// Blocked-future-work: today `AnnotateVcfConfig` has no `output_file`
// field; the output path is passed as a separate argument to
// `annotate_to_vcf_file()`. A future `OutputTarget::Stdout` enum + a
// `resolve()` rule (downgrade logger when stdout) would mirror Perl.
// Counted as 2 subtests but 1 future-work entry (shared API).
// Future-work entry:
//   porting-tests/future-work-vepyr.md
//     :: "AnnotateVcfConfig::output_target + STDOUT-implies-quiet"
//
// #[test]
// fn t34_t35_stdout_output_implies_quiet() {
//     let cfg = AnnotateVcfConfig {
//         // output_target: OutputTarget::Stdout,
//         show_progress: true,
//         ..Default::default()
//     }
//     .resolve()
//     .unwrap();
//     assert!(!cfg.show_progress, "STDOUT output must downgrade progress reporting");
// }
