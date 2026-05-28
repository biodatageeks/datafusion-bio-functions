//! Port of `ensembl-vep/t/OutputFactory_Tab.t` (v2 paradigm).
//!
//! See `porting-tests/detailed_plans/OutputFactory_Tab.md` (AUDITED
//! 2026-05-27) for the per-subtest classification table (13 behavioral
//! subtests, all `blocked-future-work` against the missing `tab_sink`
//! module in vepyr).
//!
//! **This file deliberately contains 0 `#[test]` functions.** Every subtest
//! is gated on a `tab_sink` Rust module that does not yet exist; once it
//! lands, the 13 commented-out test stubs below promote to live tests
//! simultaneously.
//!
//! ## Why this file has 0 active `#[test]` blocks
//!
//! `Bio::EnsEMBL::VEP::OutputFactory::Tab` emits VEP's `--tab` flat-tab
//! output format. Unlike VEP-tab (which packs everything into a single
//! `Extra` column), this format expands every CSQ field into its own
//! dedicated tab-separated column. The default 18-column layout is:
//! Uploaded_variation, Location, Allele, Gene, Feature, Feature_type,
//! Consequence, cDNA_position, CDS_position, Protein_position, Amino_acids,
//! Codons, Existing_variation, IMPACT, DISTANCE, STRAND, FLAGS, custom_test.
//!
//! Vepyr's only output sink today is `vcf_sink::annotate_to_vcf` — see
//! `/Users/wojtek/Documents/vepyr/datafusion-bio-functions/datafusion/bio-function-vep/src/vcf_sink.rs`.
//! There is no `FlatTabSink`, no `--output-format tab` flag, no `--fields`
//! projection knob, no per-CSQ-field column expander, and no
//! `output_hash_to_line` analogue.
//!
//! A3 RESOLVED 2026-05-25 "yes, eventually" — see
//! `porting-tests/port-status.md §"Active blockers"` item 4. Per the audit,
//! the `tab_sink` would mirror the shape of the existing `vcf_sink` and
//! consume the same per-VFOA rows that vepyr already computes; only the
//! row formatter, header generator, and `FieldSelector` projection are
//! missing.
//!
//! Missing API surface (see `future-work-vepyr.md::"Module: tab_sink"`):
//!   1. `tab_sink::FlatTabSink::new(cfg)`.
//!   2. `FlatTabSink::fields() -> &[CsqField]` (default 18-column list).
//!   3. `FlatTabSink::headers() -> Vec<String>` (24-line preamble + flat
//!      column header line).
//!   4. `FlatTabSink::format_row(row) -> String` (N-column tab-separated).
//!   5. `FlatTabSink::with_field_selector(sel: FieldSelector)` (--fields
//!      Location,HGVSc projection).
//!   6. `tab_sink::annotate_to_flat_tab(input_vcf, cache, sel) -> Vec<String>`.
//!
//! Secondary blockers (sub-test #8 plugin loading, #14 `--custom`
//! annotation source) are documented inline; both are A3-deferred.
//!
//! Coverage parity: 0 / 13 = 0% — by design until the sink lands.
//! Sibling Tier 5 ports (same A3-resolved pattern):
//!   - `tests/port_stats.rs` (Stats writer),
//!   - `tests/port_output_factory_vep_output.rs` (legacy VEP-tab sink),
//!   - `tests/port_output_factory_json.rs` (JSON sink).

// ──────────────────────────────────────────────────────────────────────────
// Blocked-future-work rows (13)
//
// Subtests numbered per detailed_plans/OutputFactory_Tab.md. Rows #1, #2,
// #3 are `use_ok` boilerplate omitted from coverage parity. Each commented
// `#[test]` block names the missing API and points at the corresponding
// entry in `future-work-vepyr.md`. When the API lands, the stub is
// uncommented, the assertion filled in (with hand-coded oracle values from
// a real-VEP-115 docker run at commit time per v2 paradigm Rule 2), and
// the coverage parity counter updated.
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #4 (~L46-47 of OutputFactory_Tab.t): `ref($of) == 'Bio::EnsEMBL::VEP::OutputFactory::Tab'`.
// Blocked: `FlatTabSink::new(cfg: &VepConfig) -> Self` constructor.
// See `future-work-vepyr.md`::`Module: tab_sink (flat-tab writer with --fields projection)`.
//
// #[test]
// fn flat_tab_sink_new_constructs() {
//     let cfg = VepConfig::default();
//     let _sink = FlatTabSink::new(&cfg);
// }

// SUBTEST #5: `$of->fields` 18-column default flat-layout list.
// Blocked: `FlatTabSink::fields(&self) -> &[CsqField]`.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[test]
// fn flat_tab_sink_default_18_columns() {
//     let sink = FlatTabSink::new(&VepConfig::default());
//     assert_eq!(sink.fields().len(), 18);
//     assert_eq!(sink.fields()[0], CsqField::UploadedVariation);
//     assert_eq!(sink.fields()[6], CsqField::Consequence);
//     assert_eq!(sink.fields()[13], CsqField::Impact);
// }

// SUBTEST #6: `headers()` 24-line preamble + `#Uploaded_variation\t...\tcustom_test`
// flat column header.
// Blocked: `FlatTabSink::headers(&self) -> Vec<String>`.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[test]
// fn flat_tab_sink_headers_preamble_shape() {
//     let sink = FlatTabSink::new(&VepConfig::default());
//     let headers = sink.headers();
//     // 24-line preamble + col-header (= 25 total) per VEP 115 default config.
//     assert!(headers.first().unwrap().starts_with("## "));
//     assert!(headers.last().unwrap().starts_with("#Uploaded_variation\t"));
// }

// SUBTEST #7: `param('merged',1)+param('custom',1)` → exactly one `SOURCE`
// field (dedup rule lives inside sink; shared with VEP-tab sink #8).
// Blocked: merged + custom dedup for `SOURCE`.
// See `future-work-vepyr.md`::`Module: tab_sink` §"Secondary gaps".
//
// #[test]
// fn flat_tab_sink_merged_custom_dedups_source() {
//     let cfg = VepConfig { merged: true, custom: true, ..Default::default() };
//     let sink = FlatTabSink::new(&cfg);
//     let count = sink.fields().iter().filter(|f| **f == CsqField::Source).count();
//     assert_eq!(count, 1);
// }

// SUBTEST #8: `headers - plugin` includes `--plugin TestPlugin` in
// command-line line. Secondary blocker: vepyr has no plugin loader.
// Blocked: `FlatTabSink::headers()` with plugin context.
// See `future-work-vepyr.md`::`Module: tab_sink` §"Secondary gaps".
//
// #[test]
// fn flat_tab_sink_headers_includes_plugin() {
//     let cfg = VepConfig { plugin: vec!["TestPlugin".into()], ..Default::default() };
//     let sink = FlatTabSink::new(&cfg);
//     let headers = sink.headers();
//     assert!(headers.iter().any(|h| h.contains("--plugin TestPlugin")));
// }

// SUBTEST #9: `output_hash_to_line({})` → 18 `-` tab-joined.
// Blocked: `FlatTabSink::format_row(&Row::default()) -> String`.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[test]
// fn flat_tab_sink_empty_row_yields_18_dashes() {
//     let sink = FlatTabSink::new(&VepConfig::default());
//     let line = sink.format_row(&CsqRow::default());
//     let cols: Vec<&str> = line.split('\t').collect();
//     assert_eq!(cols.len(), 18);
//     assert!(cols.iter().all(|c| *c == "-"));
// }

// SUBTEST #10: `output_hash_to_line({Uploaded_variation => 0})` → `0\t-\t-…`
// (preserve falsy `0` — not coerced to dash).
// Blocked: same `format_row` — falsy-id boundary handling.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[test]
// fn flat_tab_sink_falsy_zero_id_preserved() {
//     let sink = FlatTabSink::new(&VepConfig::default());
//     let row = CsqRow { uploaded_variation: Some("0".into()), ..Default::default() };
//     let line = sink.format_row(&row);
//     assert!(line.starts_with("0\t-\t"));
// }

// SUBTEST #11: `get_all_lines_by_InputBuffer` count == 744 (with
// show_ref_allele + uploaded_allele on `test_vcf`).
// Blocked: `annotate_to_flat_tab(input, cache) -> Vec<String>` end-to-end runner.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn flat_tab_sink_test_vcf_yields_744_lines() {
//     let cfg = AnnotateVcfConfig {
//         show_ref_allele: true, uploaded_allele: true, ..Default::default()
//     };
//     let lines = tab_sink::annotate_to_flat_tab(test_vcf_path(), cache_path(), &cfg, FieldSelector::default()).await.unwrap();
//     // verified via VEP 115 --tab at sink-land time.
//     assert_eq!(lines.iter().filter(|l| !l.starts_with('#')).count(), 744);
// }

// SUBTEST #12: first line content (show_ref_allele + uploaded_allele) — exact
// byte-compare. REF/UPLOADED columns are inserted between Existing_variation
// and IMPACT/DISTANCE/STRAND/FLAGS.
// Blocked: `annotate_to_flat_tab(...)` first non-header line + column insertion rule.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn flat_tab_sink_first_line_byte_compare() {
//     let cfg = AnnotateVcfConfig {
//         show_ref_allele: true, uploaded_allele: true, ..Default::default()
//     };
//     let lines = tab_sink::annotate_to_flat_tab(test_vcf_path(), cache_path(), &cfg, FieldSelector::default()).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     // verified via VEP 115:
//     assert_eq!(first, "<exact-first-line-from-real-VEP-115-tab-run>");
// }

// SUBTEST #13: last line content (FLAGS = `cds_start_NF`) byte-compare.
// Blocked: `annotate_to_flat_tab(...)` last non-header line.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn flat_tab_sink_last_line_byte_compare() {
//     let cfg = AnnotateVcfConfig {
//         show_ref_allele: true, uploaded_allele: true, ..Default::default()
//     };
//     let lines = tab_sink::annotate_to_flat_tab(test_vcf_path(), cache_path(), &cfg, FieldSelector::default()).await.unwrap();
//     let last = lines.iter().rev().find(|l| !l.starts_with('#')).unwrap();
//     // verified via VEP 115:
//     assert!(last.contains("\tcds_start_NF\t") || last.ends_with("cds_start_NF"));
// }

// SUBTEST #14 (SKIP-block): custom-annotation columns → trailing `test1`, `BAR`.
// Blocked: DOUBLE-GAP — needs (a) `tab_sink` + (b) `--custom` annotation
// source (out of current vepyr scope; A3-deferred).
// See `future-work-vepyr.md`::`Module: tab_sink` §"Secondary gaps".
//
// #[tokio::test(flavor = "multi_thread")]
// async fn flat_tab_sink_custom_annotation_trailing_columns() {
//     let cfg = AnnotateVcfConfig {
//         custom: vec![CustomAnnotation::vcf("test", "file.vcf.gz", "exact", vec!["FOO"])],
//         ..Default::default()
//     };
//     let lines = tab_sink::annotate_to_flat_tab(test_vcf_path(), cache_path(), &cfg, FieldSelector::default()).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     // verified via VEP 115 --custom:
//     assert!(first.ends_with("\ttest1\tBAR") || first.contains("\ttest1\t"));
// }

// SUBTEST #15: `--everything` mode full line — enormous line covering ALL
// CSQ columns (HGVSc, AF/gnomADe_*, SIFT/PolyPhen).
// Blocked: `annotate_to_flat_tab(everything=true)` column expansion + v115
// frequency populations + per-field stringifier for ~80 CSQ fields.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn flat_tab_sink_everything_mode_full_line() {
//     let cfg = AnnotateVcfConfig { everything: true, ..Default::default() };
//     let lines = tab_sink::annotate_to_flat_tab(test_vcf_path(), cache_path(), &cfg, FieldSelector::default()).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     // verified via VEP 115 --everything --tab:
//     assert!(first.split('\t').count() > 60); // ~80 CSQ columns in --everything.
// }

// SUBTEST #16: `--fields Location,HGVSc` projection → 2-column line.
// Blocked: `FieldSelector` API + clap `--fields` parsing.
// See `future-work-vepyr.md`::`Module: tab_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn flat_tab_sink_fields_projection_two_columns() {
//     let cfg = AnnotateVcfConfig::default();
//     let sel = FieldSelector::from_csv("Location,HGVSc").unwrap();
//     let lines = tab_sink::annotate_to_flat_tab(test_vcf_path(), cache_path(), &cfg, sel).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     assert_eq!(first.split('\t').count(), 2);
// }

// ──────────────────────────────────────────────────────────────────────────
// End of port_output_factory_tab.rs.
//
// Total substantive subtests covered: 13 (rows #1, #2, #3 are `use_ok`
// boilerplate omitted from coverage parity).
//   - 13 blocked-future-work (rows #4-#16).
//   - 0 architectural-no-analogue.
//   - 0 unit-port / integration-port / e2e-port.
//
// Coverage parity: 0 / 13 = 0% — justified in
// `detailed_plans/OutputFactory_Tab.md §"Coverage parity"`.
// Tier 5 paperwork-only port; all 13 stubs promote simultaneously when the
// `tab_sink` module lands.
// ──────────────────────────────────────────────────────────────────────────
