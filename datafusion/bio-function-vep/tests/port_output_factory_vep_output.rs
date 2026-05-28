//! Port of `ensembl-vep/t/OutputFactory_VEP_output.t` (v2 paradigm).
//!
//! See `porting-tests/detailed_plans/OutputFactory_VEP_output.md` (AUDITED
//! 2026-05-27) for the per-subtest classification table (16 behavioral
//! subtests, all `blocked-future-work` against the missing `vep_output_sink`
//! module in vepyr).
//!
//! **This file deliberately contains 0 `#[test]` functions.** Every subtest
//! is gated on a `vep_output_sink` Rust module that does not yet exist;
//! once it lands, the 16 commented-out test stubs below promote to live
//! tests simultaneously.
//!
//! ## Why this file has 0 active `#[test]` blocks
//!
//! `Bio::EnsEMBL::VEP::OutputFactory::VEP_output` emits VEP's legacy default
//! tab-separated 14-column row format: 13 fixed columns
//! (Uploaded_variation, Location, Allele, Gene, Feature, Feature_type,
//! Consequence, cDNA_position, CDS_position, Protein_position, Amino_acids,
//! Codons, Existing_variation) plus a single `Extra` column packing
//! `KEY=VALUE;...` pairs with known fields (IMPACT, DISTANCE, STRAND, FLAGS)
//! ordered before unknown.
//!
//! Vepyr's only output sink today is `vcf_sink::annotate_to_vcf` (CSQ-in-VCF
//! format) — see
//! `/Users/wojtek/Documents/vepyr/datafusion-bio-functions/datafusion/bio-function-vep/src/vcf_sink.rs`.
//! There is no `VepTabSink`, no `--output-format vep` flag, no 14-column
//! formatter, no `Extra`-column packer, and no `output_hash_to_line`
//! analogue.
//!
//! A3 RESOLVED 2026-05-25 "yes, eventually" for Stats / VEP-output / Tab
//! sinks (extended to JSON sink 2026-05-27). Per the audit, the
//! `vep_output_sink` would mirror the shape of the existing `vcf_sink` and
//! consume the same per-VFOA rows that vepyr already computes; only the
//! row formatter + header generator are missing.
//!
//! Missing API surface (see `future-work-vepyr.md::"Module: vep_output_sink"`):
//!   1. `vep_output_sink::VepTabSink::new(cfg)`.
//!   2. `VepTabSink::extra_fields() -> &[CsqField]` (default
//!      [IMPACT, DISTANCE, STRAND, FLAGS]).
//!   3. `VepTabSink::extra_field_order() -> HashMap<CsqField, usize>`.
//!   4. `VepTabSink::headers() -> Vec<String>` (26-line preamble + column
//!      header line).
//!   5. `VepTabSink::format_row(row) -> String` (14 tab-separated columns
//!      with Extra-column packer).
//!   6. `vep_output_sink::annotate_to_vep_tab(input_vcf, cache) -> Vec<String>`.
//!
//! Secondary blockers (sub-tests #10 plugin loading, #20 `--custom`
//! annotation source) are documented inline; both are A3-deferred.
//!
//! Coverage parity: 0 / 16 = 0% — by design until the sink lands.
//! Sibling Tier 5 ports (same A3-resolved pattern):
//!   - `tests/port_stats.rs` (Stats writer),
//!   - `tests/port_output_factory_tab.rs` (flat-tab sink),
//!   - `tests/port_output_factory_json.rs` (JSON sink).

// ──────────────────────────────────────────────────────────────────────────
// Blocked-future-work rows (16)
//
// Subtests numbered per detailed_plans/OutputFactory_VEP_output.md. Rows
// #1, #2, #3, #19 are `use_ok` boilerplate omitted from coverage parity.
// Each commented `#[test]` block names the missing API and points at
// the corresponding entry in `future-work-vepyr.md`. When the API lands,
// the stub is uncommented, the assertion filled in (with hand-coded oracle
// values derived from a real-VEP-115 docker run at commit time per v2
// paradigm Rule 2), and the coverage parity counter updated.
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #4 (L46-47 of OutputFactory_VEP_output.t): `ref($of) == 'Bio::EnsEMBL::VEP::OutputFactory::VEP_output'`.
// Blocked: `VepTabSink::new(cfg: &VepConfig) -> Self` constructor.
// See `future-work-vepyr.md`::`Module: vep_output_sink (legacy VEP-tab writer)`.
//
// #[test]
// fn vep_tab_sink_new_constructs() {
//     let cfg = VepConfig::default();
//     let _sink = VepTabSink::new(&cfg);
// }

// SUBTEST #5 (L49): `$of->fields == [IMPACT, DISTANCE, STRAND, FLAGS]` default.
// Blocked: `VepTabSink::extra_fields(&self) -> &[CsqField]`.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_default_extra_fields() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     assert_eq!(sink.extra_fields(),
//                &[CsqField::Impact, CsqField::Distance, CsqField::Strand, CsqField::Flags]);
// }

// SUBTEST #6 (L51): `$of->field_order == {IMPACT=>0, DISTANCE=>1, STRAND=>2, FLAGS=>3}`.
// Blocked: `VepTabSink::extra_field_order(&self) -> HashMap<CsqField, usize>`.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_extra_field_order_map() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     let order = sink.extra_field_order();
//     assert_eq!(order.get(&CsqField::Impact), Some(&0));
//     assert_eq!(order.get(&CsqField::Flags), Some(&3));
// }

// SUBTEST #7 (~L60s): `param('sift','b')` → fields gains `SIFT` at index 4.
// Blocked: dynamic Extra-field set responsive to config flag.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_with_sift_adds_sift_field() {
//     let cfg = VepConfig { sift: Some("b".into()), ..Default::default() };
//     let sink = VepTabSink::new(&cfg);
//     assert_eq!(sink.extra_field_order().get(&CsqField::Sift), Some(&4));
// }

// SUBTEST #8 (~L80s): `param('merged',1) + param('custom',1)` → exactly one
// `SOURCE` field (dedup rule lives inside sink).
// Blocked: merged + custom dedup for `SOURCE`.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_merged_custom_dedups_source_field() {
//     let cfg = VepConfig { merged: true, custom: true, ..Default::default() };
//     let sink = VepTabSink::new(&cfg);
//     let count = sink.extra_fields().iter().filter(|f| **f == CsqField::Source).count();
//     assert_eq!(count, 1);
// }

// SUBTEST #9 (~L100s-130s): `headers()` 26-line preamble + column header line
// (`## ENSEMBL VARIANT EFFECT PREDICTOR`, per-column doc lines, per-Extra-key
// doc lines, `## VEP command-line: ...`, then `#Uploaded_variation\tLocation\t…\tExtra`).
// Blocked: `VepTabSink::headers(&self) -> Vec<String>`.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_headers_preamble_shape() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     let headers = sink.headers();
//     assert!(headers.first().unwrap().starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR"));
//     assert!(headers.last().unwrap().starts_with("#Uploaded_variation\t"));
// }

// SUBTEST #10 (~L140s): `headers - plugin` includes `--plugin TestPlugin`
// in command-line line. Secondary blocker: vepyr has no plugin loader.
// Blocked: `VepTabSink::headers()` with plugin context (also requires
// `PluginLoader` infrastructure — separately blocked).
// See `future-work-vepyr.md`::`Module: vep_output_sink` §"Secondary gaps".
//
// #[test]
// fn vep_tab_sink_headers_includes_plugin_in_command_line() {
//     let cfg = VepConfig { plugin: vec!["TestPlugin".into()], ..Default::default() };
//     let sink = VepTabSink::new(&cfg);
//     let headers = sink.headers();
//     assert!(headers.iter().any(|h| h.contains("--plugin TestPlugin")));
// }

// SUBTEST #11 (~L160s): `output_hash_to_line({})` → 14 `-` tab-joined.
// Blocked: `VepTabSink::format_row(&Row::default()) -> String`.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_empty_row_yields_14_dashes() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     let line = sink.format_row(&CsqRow::default());
//     let cols: Vec<&str> = line.split('\t').collect();
//     assert_eq!(cols.len(), 14);
//     assert!(cols.iter().all(|c| *c == "-"));
// }

// SUBTEST #12 (~L180s): `output_hash_to_line({Uploaded_variation => 0})` →
// `0\t-\t-…` (preserve falsy `0` — not coerced to dash).
// Blocked: same `format_row` — falsy-id boundary handling.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_falsy_zero_id_preserved() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     let row = CsqRow { uploaded_variation: Some("0".into()), ..Default::default() };
//     let line = sink.format_row(&row);
//     assert!(line.starts_with("0\t-\t"));
// }

// SUBTEST #13 (~L200s): `output_hash_to_line({Existing_variation=>'rs123', Foo=>'bar'})`
// → 11 dashes + `rs123` + Extra=`Foo=bar`. Extra column packs unknown keys.
// Blocked: `format_row` Extra-column unknown-key passthrough.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_extra_column_passes_unknown_keys() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     let row = CsqRow {
//         existing_variation: Some("rs123".into()),
//         extra: vec![("Foo".into(), "bar".into())],
//         ..Default::default()
//     };
//     let line = sink.format_row(&row);
//     assert!(line.ends_with("\tFoo=bar"));
// }

// SUBTEST #14 (~L220s): `output_hash_to_line({Existing_variation, Foo, IMPACT})`
// → Extra=`IMPACT=HIGH;Foo=bar` (known-first order).
// Blocked: `format_row` Extra-column known-fields-first ordering invariant.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[test]
// fn vep_tab_sink_extra_orders_known_before_unknown() {
//     let sink = VepTabSink::new(&VepConfig::default());
//     let row = CsqRow {
//         impact: Some("HIGH".into()),
//         extra: vec![("Foo".into(), "bar".into())],
//         ..Default::default()
//     };
//     let line = sink.format_row(&row);
//     assert!(line.ends_with("IMPACT=HIGH;Foo=bar"));
// }

// SUBTEST #15 (~L260s): `get_all_lines_by_InputBuffer` count == 744 (with
// show_ref_allele + uploaded_allele on `test_vcf`).
// Blocked: `annotate_to_vep_tab(input_vcf, cache) -> Vec<String>` end-to-end runner.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn vep_tab_sink_test_vcf_yields_744_lines() {
//     let cfg = AnnotateVcfConfig {
//         show_ref_allele: true, uploaded_allele: true, ..Default::default()
//     };
//     let lines = vep_output_sink::annotate_to_vep_tab(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     // verified via VEP 115 at sink-land time:
//     //   docker run --rm -v $(pwd)/_cache115:/cache:ro ... vep ... test.vcf | wc -l
//     assert_eq!(lines.iter().filter(|l| !l.starts_with('#')).count(), 744);
// }

// SUBTEST #16 (~L290s): first line content (show_ref_allele + uploaded_allele)
// exact byte-compare.
// Blocked: `annotate_to_vep_tab(...)` first non-header line.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn vep_tab_sink_first_line_byte_compare() {
//     let cfg = AnnotateVcfConfig {
//         show_ref_allele: true, uploaded_allele: true, ..Default::default()
//     };
//     let lines = vep_output_sink::annotate_to_vep_tab(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     // verified via VEP 115:
//     assert_eq!(first, "<exact-first-line-from-real-VEP-115-run>");
// }

// SUBTEST #17 (~L320s): last line content (FLAGS = `cds_start_NF`) byte-compare.
// Blocked: `annotate_to_vep_tab(...)` last non-header line.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn vep_tab_sink_last_line_byte_compare() {
//     let cfg = AnnotateVcfConfig {
//         show_ref_allele: true, uploaded_allele: true, ..Default::default()
//     };
//     let lines = vep_output_sink::annotate_to_vep_tab(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let last = lines.iter().rev().find(|l| !l.starts_with('#')).unwrap();
//     // verified via VEP 115:
//     assert!(last.contains("FLAGS=cds_start_NF"));
// }

// SUBTEST #18 (~L350s-400s): `--everything` mode last column matches long
// `IMPACT=...;...;MAX_AF_POPS=gnomADe_AFR` string.
// Blocked: `annotate_to_vep_tab(everything=true)` Extra-column packer with
// v115 gnomADe populations.
// See `future-work-vepyr.md`::`Module: vep_output_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn vep_tab_sink_everything_mode_extra_column() {
//     let cfg = AnnotateVcfConfig { everything: true, ..Default::default() };
//     let lines = vep_output_sink::annotate_to_vep_tab(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     let extra = first.rsplit('\t').next().unwrap();
//     // verified via VEP 115 --everything:
//     assert!(extra.contains("IMPACT="));
//     assert!(extra.contains("MAX_AF_POPS=gnomADe_AFR"));
// }

// SUBTEST #20 (~L440s, SKIP-block): custom-annotation Extra →
// `test=test1;test_FILTER=PASS;test_FOO=BAR`.
// Blocked: DOUBLE-GAP — needs (a) `vep_output_sink` + (b) `--custom`
// annotation source (out of current vepyr scope; A3-deferred).
// See `future-work-vepyr.md`::`Module: vep_output_sink` §"Secondary gaps".
//
// #[tokio::test(flavor = "multi_thread")]
// async fn vep_tab_sink_custom_annotation_extra_keys() {
//     let cfg = AnnotateVcfConfig {
//         custom: vec![CustomAnnotation::vcf("test", "file.vcf.gz", "exact", vec!["FOO"])],
//         ..Default::default()
//     };
//     let lines = vep_output_sink::annotate_to_vep_tab(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let first = lines.iter().find(|l| !l.starts_with('#')).unwrap();
//     let extra = first.rsplit('\t').next().unwrap();
//     // verified via VEP 115 --custom:
//     assert!(extra.contains("test=test1"));
//     assert!(extra.contains("test_FILTER=PASS"));
//     assert!(extra.contains("test_FOO=BAR"));
// }

// ──────────────────────────────────────────────────────────────────────────
// End of port_output_factory_vep_output.rs.
//
// Total substantive subtests covered: 16 (rows #1, #2, #3, #19 are `use_ok`
// boilerplate omitted from coverage parity).
//   - 16 blocked-future-work (rows #4-#18, #20).
//   - 0 architectural-no-analogue.
//   - 0 unit-port / integration-port / e2e-port.
//
// Coverage parity: 0 / 16 = 0% — justified in
// `detailed_plans/OutputFactory_VEP_output.md §"Coverage parity"`.
// Tier 5 paperwork-only port; all 16 stubs promote simultaneously when the
// `vep_output_sink` module lands.
// ──────────────────────────────────────────────────────────────────────────
