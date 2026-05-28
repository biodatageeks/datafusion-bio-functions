//! Port of `ensembl-vep/t/Stats.t` (v2 paradigm).
//!
//! See `porting-tests/detailed_plans/Stats.md` (AUDITED 2026-05-27) for the
//! authoritative per-subtest classification table (27 behavioral subtests,
//! all `blocked-future-work` against the missing `Stats` writer module).
//!
//! **This file deliberately contains 0 `#[test]` functions.** Every subtest
//! is gated on a `Stats` Rust module + `FinishedStats::write_text` /
//! `write_html` serialisers that do not yet exist; once they land, the 27
//! commented-out test stubs below promote to live tests simultaneously.
//!
//! ## Why this file has 0 active `#[test]` blocks
//!
//! `Bio::EnsEMBL::VEP::Stats` is VEP's per-run summary-statistics writer.
//! It emits two sidecar files next to the main output VCF:
//!
//! - `*_summary.txt` — `[VEP run statistics]` header + ASCII tables.
//! - `*_summary.html` — HTML page with embedded pie/bar charts.
//!
//! The Perl `Stats` object is constructed empty (`Stats->new` returns a
//! blessed hash with `{stats => {counters => {}}}`), accumulates 12 counter
//! classes as variants/VFOAs are logged (`var_count`, `classes`, `var_cons`,
//! `chr` × 1Mb-bin, `allele_changes`, `consequences`, `protein_pos`, `gene`,
//! `transcript`, `SIFT`/`PolyPhen`, `filtered_variants`), and at end-of-run
//! produces a `finished_stats` structure that `dump_text` / `dump_html`
//! serialize.
//!
//! Vepyr has no per-run summary-statistics writer at all today. The closest
//! existing concept is `EntityStats` in `cache_builder.rs:93`, but that
//! tracks cache-build progress, not per-annotation-run variant counters.
//! There is no `--stats-file` flag, no `*_summary.txt` sidecar, no HTML
//! chart renderer, and no per-VFOA counter accumulator.
//!
//! A3 RESOLVED 2026-05-25 ("should vepyr emit per-run summary stats?") =
//! "yes, eventually" — see `porting-tests/port-status.md §"Active blockers"`
//! item 4. Per the audit, the Stats writer shares architectural status with
//! the three other planned-but-unbuilt sinks (vep_output_sink, tab_sink,
//! json_sink): different output format, same input data model, no
//! architectural precludence.
//!
//! Missing API surface (see `porting-tests/future-work-vepyr.md` entries):
//!   1. `Stats::new + Counters` — base struct + initial empty state.
//!   2. `Stats::info / set_info` — CLI args / file path map.
//!   3. `Stats::start_time / run_time / end_time` — chrono wall-clock.
//!   4. `Stats::log_lines_read / increment_filtered_variants` — counters.
//!   5. `Stats::log_variant` — per-variant counter update.
//!   6. `Stats::log_csq_row` — per-VFOA counter update.
//!   7. `Stats::log_sift_polyphen + Stats::finish() -> FinishedStats`.
//!   8. `FinishedStats::write_text / write_html` — serialisers.
//!
//! Coverage parity: 0 / 27 = 0% — by design until the Stats module lands.
//! Sibling Tier 5 ports (same A3-resolved pattern):
//!   - `tests/port_output_factory_vep_output.rs` (legacy VEP-tab sink),
//!   - `tests/port_output_factory_tab.rs` (flat-tab sink),
//!   - `tests/port_output_factory_json.rs` (JSON sink).

// ──────────────────────────────────────────────────────────────────────────
// Blocked-future-work rows (27)
//
// Subtests numbered per detailed_plans/Stats.md (rows #1, #2 are `use_ok`
// boilerplate omitted from coverage parity). Each commented `#[test]` block
// names the missing vepyr API and points at the corresponding entry in
// `future-work-vepyr.md`. When the API lands, the corresponding stub is
// uncommented, the assertion filled in (with hand-coded oracle values per
// v2 paradigm Rule 2), and the coverage parity counter updated.
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #3 (L41 of Stats.t): `Stats->new` is defined.
// Blocked: `Stats::new() -> Self` does not exist.
// See `future-work-vepyr.md`::`Stats::new + Counters struct`.
//
// #[test]
// fn stats_new_is_defined() {
//     let s = Stats::new();
//     // constructor returns a value; no panic.
//     drop(s);
// }

// SUBTEST #4 (L42-50): initial shape `{stats => {counters => {}}, …}`.
// Blocked: same as #3 — `Stats::counters() -> &Counters` (or
// `Counters::is_empty(&self) -> bool`) accessor.
// See `future-work-vepyr.md`::`Stats::new + Counters struct`.
//
// #[test]
// fn stats_new_initial_counters_empty() {
//     let s = Stats::new();
//     assert!(s.counters().is_empty());
// }

// SUBTEST #5 (L57): `info` returns hashref.
// Blocked: `Stats::info() -> &HashMap<String, String>` + `set_info(k, v)`.
// See `future-work-vepyr.md`::`Stats::info / set_info`.
//
// #[test]
// fn stats_info_round_trip() {
//     let mut s = Stats::new();
//     s.set_info("foo", "bar");
//     assert_eq!(s.info().get("foo"), Some(&"bar".to_string()));
// }

// SUBTEST #6 (L59): `start_time == get_time()` formatted timestamp.
// Blocked: `Stats::start_time() -> &str` (RFC 3339 or VEP-format).
// See `future-work-vepyr.md`::`Stats::start_time / run_time / end_time`.
//
// #[test]
// fn stats_start_time_set_at_new() {
//     let s = Stats::new();
//     assert!(!s.start_time().is_empty());
// }

// SUBTEST #7 (L60): `run_time_start == time()` epoch seconds.
// Blocked: internal `Instant::now()` capture (folds into Stats::new entry).
// See `future-work-vepyr.md`::`Stats::start_time / run_time / end_time`.
//
// #[test]
// fn stats_run_time_start_captured() {
//     let s = Stats::new();
//     // implementation detail: vepyr stores Instant; no direct accessor needed
//     // beyond `run_time()` (#8).
//     drop(s);
// }

// SUBTEST #8 (L63): `run_time =~ /^\d+$/` integer elapsed seconds.
// Blocked: `Stats::run_time() -> u64`.
// See `future-work-vepyr.md`::`Stats::start_time / run_time / end_time`.
//
// #[test]
// fn stats_run_time_returns_elapsed_seconds() {
//     let s = Stats::new();
//     std::thread::sleep(std::time::Duration::from_secs(1));
//     assert!(s.run_time() >= 1);
// }

// SUBTEST #9 (L64): `end_time == get_time()` formatted at finalize.
// Blocked: `Stats::end_time() -> Option<&str>` set when `finish()` is called.
// See `future-work-vepyr.md`::`Stats::start_time / run_time / end_time`.
//
// #[test]
// fn stats_end_time_set_at_finish() {
//     let s = Stats::new();
//     let finished = s.finish();
//     assert!(!finished.end_time().is_empty());
// }

// SUBTEST #10 (L66-67): `log_lines_read(10) → lines_read == 10`.
// Blocked: `Stats::log_lines_read(n)` + `Stats::lines_read()` getter.
// See `future-work-vepyr.md`::`Stats::log_lines_read / increment_filtered_variants`.
//
// #[test]
// fn stats_log_lines_read_accumulates() {
//     let mut s = Stats::new();
//     s.log_lines_read(10);
//     assert_eq!(s.lines_read(), 10);
// }

// SUBTEST #11 (L69-76): `increment_filtered_variants(5)` counter.
// Blocked: `Stats::increment_filtered_variants(n)` + getter. Companion to
// future `--filter` flag (also blocked).
// See `future-work-vepyr.md`::`Stats::log_lines_read / increment_filtered_variants`.
//
// #[test]
// fn stats_filtered_variants_accumulates() {
//     let mut s = Stats::new();
//     s.increment_filtered_variants(5);
//     assert_eq!(s.filtered_variants(), 5);
// }

// SUBTEST #12 (L85-106): first SNV `log_VariationFeature(vf1)` — classes=SNV,
// var_cons=missense, chr=21/25M, allele_changes=C/T, var_count=1.
// Blocked: `Stats::log_variant(vcf_row: &VcfRow)` updating 5 counters.
// See `future-work-vepyr.md`::`Stats::log_variant`.
//
// #[test]
// fn stats_log_variant_first_snv() {
//     let mut s = Stats::new();
//     let row: VcfRow = /* SNV: 21 25585733 C T missense */ todo!();
//     s.log_variant(&row);
//     let c = s.counters();
//     assert_eq!(c.classes.get("SNV"), Some(&1));
//     assert_eq!(c.var_cons.get("missense_variant"), Some(&1));
//     assert_eq!(c.var_count, 1);
//     // chr_bins[21][25_000_000] == 1; allele_changes[C/T] == 1
// }

// SUBTEST #13 (L108-130): second `log_VariationFeature` — counters increment.
// Blocked: same as #12.
// See `future-work-vepyr.md`::`Stats::log_variant`.
//
// #[test]
// fn stats_log_variant_second_call_increments() {
//     let mut s = Stats::new();
//     let row1: VcfRow = todo!();
//     let row2: VcfRow = todo!();
//     s.log_variant(&row1);
//     s.log_variant(&row2);
//     assert_eq!(s.counters().var_count, 2);
//     assert_eq!(s.counters().classes.get("SNV"), Some(&2));
// }

// SUBTEST #14 (L132-156): multi-allele 'T/C/G' — extra ALT tallied as
// separate `T/G` entry. One-row-per-ALT in `allele_changes`.
// Blocked: `Stats::log_variant` multi-ALT handling.
// See `future-work-vepyr.md`::`Stats::log_variant`.
//
// #[test]
// fn stats_log_variant_multi_allele() {
//     let mut s = Stats::new();
//     let row: VcfRow = /* 21 X T C,G */ todo!();
//     s.log_variant(&row);
//     let c = s.counters();
//     assert_eq!(c.allele_changes.get("T/C"), Some(&1));
//     assert_eq!(c.allele_changes.get("T/G"), Some(&1));
// }

// SUBTEST #15 (L158-174): `log_TranscriptVariationAllele` — per-VFOA
// counters: protein_pos=9, gene=ENSG00000154719. Does NOT bump
// `transcript` or `consequences`.
// Blocked: `Stats::log_csq_row(csq: &CsqRow)` for transcript rows.
// See `future-work-vepyr.md`::`Stats::log_csq_row`.
//
// #[test]
// fn stats_log_csq_row_transcript_only() {
//     let mut s = Stats::new();
//     let csq: CsqRow = /* transcript-level row */ todo!();
//     s.log_csq_row(&csq);
//     let c = s.counters();
//     assert_eq!(c.protein_pos.get(&9), Some(&1));
//     assert_eq!(c.gene.get("ENSG00000154719"), Some(&1));
// }

// SUBTEST #16 (L176-196): `log_VariationFeatureOverlapAllele` — full VFOA
// counters: protein_pos, gene, transcript=ENST00000352957,
// consequences=missense. Superset of #15.
// Blocked: `Stats::log_csq_row` full VFOA branch.
// See `future-work-vepyr.md`::`Stats::log_csq_row`.
//
// #[test]
// fn stats_log_csq_row_full_vfoa() {
//     let mut s = Stats::new();
//     let csq: CsqRow = /* full VFOA */ todo!();
//     s.log_csq_row(&csq);
//     let c = s.counters();
//     assert_eq!(c.transcript.get("ENST00000352957"), Some(&1));
//     assert_eq!(c.consequences.get("missense_variant"), Some(&1));
// }

// SUBTEST #17 (L199-219): SV `log_TranscriptStructuralVariationAllele` —
// SV-variant counters. Vepyr SV support is itself blocked (port-status.md
// active blocker #2).
// Blocked: `Stats::log_csq_row` SV branch.
// See `future-work-vepyr.md`::`Stats::log_csq_row`.
//
// #[test]
// fn stats_log_csq_row_sv_transcript() {
//     let mut s = Stats::new();
//     let csq: CsqRow = /* SV-transcript row */ todo!();
//     s.log_csq_row(&csq);
//     assert!(s.counters().consequences.contains_key("feature_elongation"));
// }

// SUBTEST #18 (L221-242): SV `log_VariationFeatureOverlapAllele` DUP —
// consequences=feature_elongation + coding_sequence_variant.
// Blocked: same as #17.
// See `future-work-vepyr.md`::`Stats::log_csq_row`.
//
// #[test]
// fn stats_log_csq_row_sv_dup_vfoa() {
//     let mut s = Stats::new();
//     let csq: CsqRow = /* SV DUP VFOA */ todo!();
//     s.log_csq_row(&csq);
//     let c = s.counters();
//     assert!(c.consequences.contains_key("feature_elongation"));
//     assert!(c.consequences.contains_key("coding_sequence_variant"));
// }

// SUBTEST #19 (L245-263): intergenic `log_VariationFeatureOverlapAllele` —
// consequences=intergenic_variant only (no gene/transcript/protein_pos).
// Blocked: `Stats::log_csq_row` intergenic branch.
// See `future-work-vepyr.md`::`Stats::log_csq_row`.
//
// #[test]
// fn stats_log_csq_row_intergenic() {
//     let mut s = Stats::new();
//     let csq: CsqRow = /* intergenic row */ todo!();
//     s.log_csq_row(&csq);
//     let c = s.counters();
//     assert_eq!(c.consequences.get("intergenic_variant"), Some(&1));
//     assert!(c.gene.is_empty());
//     assert!(c.transcript.is_empty());
// }

// SUBTEST #20 (L266-276): `log_sift_polyphen('SIFT', 'deleterious')` —
// per-tool prediction tally. SIFT/PolyPhen are themselves blocked engine
// features.
// Blocked: `Stats::log_sift_polyphen(tool, class)`.
// See `future-work-vepyr.md`::`Stats::log_sift_polyphen + Stats::finish() → FinishedStats`.
//
// #[test]
// fn stats_log_sift_polyphen_tally() {
//     let mut s = Stats::new();
//     s.log_sift_polyphen("SIFT", "deleterious");
//     assert_eq!(s.counters().sift.get("deleterious"), Some(&1));
// }

// SUBTEST #21 (L287-320): `finished_stats general_stats` — 7-row table
// (Lines read, Variants processed, Variants filtered out, Novel/existing,
// Overlapped genes/transcripts/regulatory features).
// Blocked: `FinishedStats::general -> Vec<(String, String)>`.
// See `future-work-vepyr.md`::`Stats::log_sift_polyphen + Stats::finish() → FinishedStats`.
//
// #[test]
// fn stats_finished_general_has_seven_rows() {
//     let s = Stats::new();
//     let finished = s.finish();
//     assert_eq!(finished.general.len(), 7);
// }

// SUBTEST #22 (L322-405): `finished_stats charts data` — 7 chart buckets
// (class, MSC, consequences, var-cons, chr, position-on-chr histogram,
// AF-bin histogram).
// Blocked: `FinishedStats::charts -> Vec<ChartData>`.
// See `future-work-vepyr.md`::`Stats::log_sift_polyphen + Stats::finish() → FinishedStats`.
//
// #[test]
// fn stats_finished_charts_has_seven_buckets() {
//     let s = Stats::new();
//     let finished = s.finish();
//     assert_eq!(finished.charts.len(), 7);
// }

// SUBTEST #23 (L408): `run_stats version =~ /\d+ \(\d+\)/` — version
// pattern `<release> (<sub>)`. Vepyr's version is CARGO_PKG_VERSION
// (semver); test would adapt to vepyr's format.
// Blocked: `FinishedStats::run_info` first row.
// See `future-work-vepyr.md`::`Stats::log_sift_polyphen + Stats::finish() → FinishedStats`.
//
// #[test]
// fn stats_finished_run_info_includes_version() {
//     let s = Stats::new();
//     let finished = s.finish();
//     let (label, _value) = &finished.run_info[0];
//     assert!(label.contains("version"));
// }

// SUBTEST #24 (L409-423): `run_stats remaining headers` — 9-header list
// (VEP version, Annotation sources, Species, Command line options,
// Start time, End time, Run time, Input file, Output file).
// Blocked: `FinishedStats::run_info` row labels.
// See `future-work-vepyr.md`::`Stats::log_sift_polyphen + Stats::finish() → FinishedStats`.
//
// #[test]
// fn stats_finished_run_info_has_nine_rows() {
//     let s = Stats::new();
//     let finished = s.finish();
//     assert_eq!(finished.run_info.len(), 9);
// }

// SUBTEST #25 (L427-429): `dump_text` starts with `[VEP run statistics]`.
// Blocked: `FinishedStats::write_text(w)` emits header as first line.
// See `future-work-vepyr.md`::`FinishedStats::write_text / write_html`.
//
// #[test]
// fn stats_dump_text_starts_with_header() {
//     let s = Stats::new();
//     let finished = s.finish();
//     let mut buf = Vec::new();
//     finished.write_text(&mut buf).unwrap();
//     let text = String::from_utf8(buf).unwrap();
//     assert!(text.starts_with("[VEP run statistics]"));
// }

// SUBTEST #26 (L428): `dump_text` returns true (non-error write).
// Blocked: `Result<(), io::Error>` from `write_text`.
// See `future-work-vepyr.md`::`FinishedStats::write_text / write_html`.
//
// #[test]
// fn stats_dump_text_succeeds() {
//     let s = Stats::new();
//     let finished = s.finish();
//     let mut buf = Vec::new();
//     assert!(finished.write_text(&mut buf).is_ok());
// }

// SUBTEST #27 (L433-435): `dump_html` starts with `<html>` tag.
// Blocked: `FinishedStats::write_html(w)` emits HTML doc.
// See `future-work-vepyr.md`::`FinishedStats::write_text / write_html`.
//
// #[test]
// fn stats_dump_html_starts_with_html_tag() {
//     let s = Stats::new();
//     let finished = s.finish();
//     let mut buf = Vec::new();
//     finished.write_html(&mut buf).unwrap();
//     let html = String::from_utf8(buf).unwrap();
//     assert!(html.starts_with("<html>"));
// }

// SUBTEST #28 (L434): `dump_html` returns true (non-error write).
// Blocked: `Result<(), io::Error>` from `write_html`.
// See `future-work-vepyr.md`::`FinishedStats::write_text / write_html`.
//
// #[test]
// fn stats_dump_html_succeeds() {
//     let s = Stats::new();
//     let finished = s.finish();
//     let mut buf = Vec::new();
//     assert!(finished.write_html(&mut buf).is_ok());
// }

// SUBTEST #29 (L444-451): `HTML::Lint` clean — validates HTML tree. Rust
// analogue is structural-parse-success (e.g., `tl::parse(&html, ..).is_ok()`),
// not byte-for-byte HTML 4 validation.
// Blocked: same writer as #27.
// See `future-work-vepyr.md`::`FinishedStats::write_text / write_html`.
//
// #[test]
// fn stats_dump_html_is_structurally_valid() {
//     let s = Stats::new();
//     let finished = s.finish();
//     let mut buf = Vec::new();
//     finished.write_html(&mut buf).unwrap();
//     let html = String::from_utf8(buf).unwrap();
//     // structural parse OK (not full HTML 4 validation):
//     // assert!(tl::parse(&html, tl::ParserOptions::default()).is_ok());
//     assert!(html.contains("</html>"));
// }

// ──────────────────────────────────────────────────────────────────────────
// End of port_stats.rs.
//
// Total substantive subtests covered: 27 (rows #1, #2 are `use_ok`
// boilerplate omitted from coverage parity).
//   - 27 blocked-future-work (rows #3 through #29).
//   - 0 architectural-no-analogue.
//   - 0 unit-port / integration-port / e2e-port.
//
// Coverage parity: 0 / 27 = 0% — justified in
// `detailed_plans/Stats.md §"Coverage parity"`. Tier 5 paperwork-only port.
// ──────────────────────────────────────────────────────────────────────────
