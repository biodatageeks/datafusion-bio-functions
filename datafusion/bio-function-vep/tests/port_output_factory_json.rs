//! Port of `ensembl-vep/t/OutputFactory_JSON.t` (v2 paradigm).
//!
//! See `porting-tests/detailed_plans/OutputFactory_JSON.md` (AUDITED
//! 2026-05-27, promoted from EXCLUDE → PORT 2026-05-27) for the per-subtest
//! classification table (16 behavioral subtests, all `blocked-future-work`
//! against the missing `json_sink` module in vepyr).
//!
//! **This file deliberately contains 0 `#[test]` functions.** Every subtest
//! is gated on a `json_sink` module that does not yet exist; once it lands,
//! the 16 commented-out test stubs below promote to live tests simultaneously.
//!
//! A3 RESOLVED 2026-05-25 ("yes, eventually" for Stats / VEP-output / Tab
//! sinks) extended by Wojtek 2026-05-27 to cover the JSON sink. Per the
//! audit, the JSON sink shares the same architectural status as the other
//! three planned sinks: different output format, same input data model
//! (per-VFOA, per-coloc, per-regfeat rows that vepyr already computes for
//! VCF output), no architectural precludence.
//!
//! ## Why this file has 0 active `#[test]` blocks
//!
//! `Bio::EnsEMBL::VEP::OutputFactory::JSON` emits one JSON object per VCF
//! input line (vs CSQ-INFO-in-VCF, tab-separated, or VEP-tab Extra-column).
//! Per-transcript-consequence data nests as `transcript_consequences[]`,
//! co-located-variant data as `colocated_variants[]`, regulatory features
//! as `regulatory_feature_consequences[]`, and `--custom` overlays as
//! `custom_annotations{}`.
//!
//! Vepyr's only output sink today is `vcf_sink::annotate_to_vcf` — see
//! `/Users/wojtek/Documents/vepyr/datafusion-bio-functions/datafusion/bio-function-vep/src/vcf_sink.rs`.
//! There is no `JsonSink`, no `--output-format json` flag, no per-row JSON
//! serializer.
//!
//! Missing API surface (see `porting-tests/future-work-vepyr.md` entries
//! added 2026-05-27 as part of this port's audit):
//!   1. `json_sink::annotate_to_json` — main JSON sink runner.
//!   2. `json_sink::JsonOutput` — typed per-row JSON shape.
//!   3. `json_sink::format_vfoa_as_json` — per-transcript-consequence object.
//!   4. `json_sink::format_colocated_as_json` — per-colocated-variant object.
//!   5. `json_sink::format_regfeat_as_json` — per-regulatory-feature object.
//!
//! Coverage parity: 0 / 16 = 0% — by design until json_sink lands.
//! Sibling Tier 5 ports (same A3-resolved pattern):
//!   - `tests/port_stats.rs` (Stats writer),
//!   - `tests/port_output_factory_vep_output.rs` (legacy VEP-tab sink),
//!   - `tests/port_output_factory_tab.rs` (flat-tab sink).

// ──────────────────────────────────────────────────────────────────────────
// Blocked-future-work rows (16)
//
// Subtests numbered per detailed_plans/OutputFactory_JSON.md. Rows #1-#5
// and #20 are `use_ok` boilerplate omitted from coverage parity. Each
// commented `#[test]` block names the missing API and points at the
// corresponding entry in `future-work-vepyr.md`. When the API lands, the
// stub is uncommented, the assertion filled in (with hand-coded oracle
// values from a real-VEP-115 docker run at commit time per v2 paradigm
// Rule 2), and the coverage parity counter updated.
// ──────────────────────────────────────────────────────────────────────────

// SUBTEST #6 (L52 of OutputFactory_JSON.t): `ref($of) == 'Bio::EnsEMBL::VEP::OutputFactory::JSON'`.
// Blocked: `json_sink::JsonSink::new(cfg: &VepConfig) -> Self` constructor.
// See `future-work-vepyr.md`::`Module: json_sink (JSON output writer)`.
//
// #[test]
// fn json_sink_new_constructs() {
//     let cfg = VepConfig::default();
//     let _sink = JsonSink::new(&cfg);
// }

// SUBTEST #7 (L59): `is_deeply($of->headers, [])` — JSON has no preamble.
// Blocked: `JsonSink::headers(&self) -> Vec<String>` returning empty Vec.
// (Accessor must exist for symmetry with VCF/Tab/VEP-output sinks.)
// See `future-work-vepyr.md`::`Module: json_sink`.
//
// #[test]
// fn json_sink_headers_is_empty() {
//     let sink = JsonSink::new(&VepConfig::default());
//     assert!(sink.headers().is_empty());
// }

// SUBTEST #8 (L66-114): `add_VariationFeatureOverlapAllele_info($vf, {})`
// returns nested transcript_consequences + most_severe_consequence for
// rs142513484 (array of 3 transcript-consequence objects with
// gene_id/transcript_id/variant_allele/consequence_terms/impact/strand/
// cdna_*/cds_*/protein_*/codons/amino_acids).
// Blocked: `json_sink::format_vfoa_as_json(vfoa: &Vfoa, cfg) -> serde_json::Value`.
// See `future-work-vepyr.md`::`json_sink::format_vfoa_as_json`.
//
// #[test]
// fn json_sink_format_vfoa_default_shape() {
//     let vfoa: Vfoa = /* rs142513484 first VFOA */ todo!();
//     let cfg = AnnotateVcfConfig::default();
//     let json = json_sink::format_vfoa_as_json(&vfoa, &cfg);
//     // verified via VEP 115 --json at sink-land time:
//     assert_eq!(json["gene_id"], "ENSG00000154719");
//     assert_eq!(json["transcript_id"], "ENST00000307301");
//     assert_eq!(json["variant_allele"], "T");
//     assert!(json["consequence_terms"].as_array().is_some());
// }

// SUBTEST #9 (L125-166): `add_colocated_variant_info_JSON({}, [$frequency_hash], $ex)`
// returns single-allele colocated hash for rs142513484 with
// frequencies/var_synonyms/pubmed/clin_sig/minor_allele_freq/allele_string.
// Blocked: `json_sink::format_colocated_as_json(coloc, matched_alleles)`.
// See `future-work-vepyr.md`::`json_sink::format_colocated_as_json`.
//
// #[test]
// fn json_sink_format_colocated_single_allele() {
//     let coloc: ColocatedEntry = /* rs142513484 */ todo!();
//     let matched: Vec<String> = vec!["T".into()];
//     let json = json_sink::format_colocated_as_json(&coloc, &matched);
//     // verified via VEP 115 --json:
//     assert_eq!(json["id"], "rs142513484");
//     assert!(json["frequencies"]["T"]["af"].is_number());
// }

// SUBTEST #10 (L189-219): `add_colocated_variant_info_JSON` multi-allelic
// 21:25891785 G/A/T — only matched-allele A frequencies emit, not T.
// Blocked: same `format_colocated_as_json` — multi-allele frequency filtering.
// See `future-work-vepyr.md`::`json_sink::format_colocated_as_json`.
//
// #[test]
// fn json_sink_format_colocated_multi_allele_filters_by_matched() {
//     let coloc: ColocatedEntry = /* rs145564988 multi-allele */ todo!();
//     let matched: Vec<String> = vec!["A".into()];
//     let json = json_sink::format_colocated_as_json(&coloc, &matched);
//     assert!(json["frequencies"].get("A").is_some());
//     assert!(json["frequencies"].get("T").is_none()); // not matched, not emitted
// }

// SUBTEST #11 (L226): `scalar @lines == scalar @{$ib->buffer}` — one JSON
// object per VCF input row.
// Blocked: `json_sink::annotate_to_json(input_vcf, cache, cfg) -> Vec<String>`.
// See `future-work-vepyr.md`::`Module: json_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_line_count_matches_input_rows() {
//     let cfg = AnnotateVcfConfig::default();
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     // verified via VEP 115 --json: number of JSON objects = number of VCF input rows.
//     assert_eq!(lines.len(), test_vcf_row_count());
// }

// SUBTEST #12 (L231-288): `$json->decode($lines[0])` matches full expected
// hash (default mode) — input/assembly_name/seq_region_name/start/end/
// strand/id/allele_string/most_severe_consequence/transcript_consequences[3].
// Blocked: `json_sink::annotate_to_json` + `format_vfoa_as_json` integration.
// See `future-work-vepyr.md`::`Module: json_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_first_line_default_mode() {
//     let cfg = AnnotateVcfConfig::default();
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     // verified via VEP 115 --json at sink-land time:
//     assert_eq!(obj["assembly_name"], "GRCh38");
//     assert_eq!(obj["seq_region_name"], "21");
//     assert_eq!(obj["id"], "rs142513484");
//     assert_eq!(obj["most_severe_consequence"], "missense_variant");
//     assert_eq!(obj["transcript_consequences"].as_array().unwrap().len(), 3);
// }

// SUBTEST #13 (L298-310): everything+regulatory: line 28 carries
// `regulatory_feature_consequences: [{consequence_terms, variant_allele,
// regulatory_feature_id, impact, biotype}]`.
// Blocked: `json_sink::format_regfeat_as_json(regfeat) -> serde_json::Value`.
// See `future-work-vepyr.md`::`json_sink::format_regfeat_as_json`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_regulatory_feature_consequences() {
//     let cfg = AnnotateVcfConfig { everything: true, regulatory: true, ..Default::default() };
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[27]).unwrap();
//     // verified via VEP 115 --json --regulatory at sink-land time:
//     let regs = obj["regulatory_feature_consequences"].as_array().unwrap();
//     assert!(!regs.is_empty());
//     assert!(regs[0]["regulatory_feature_id"].is_string());
// }

// SUBTEST #14 (L312-445): `$json->decode($lines[0])` everything mode — first
// JSON line carries colocated_variants (with frequencies and var_synonyms),
// variant_class=SNV, transcript_consequences with HGVSc/HGVSp/SIFT/PolyPhen/
// swissprot/uniparc/biotype/canonical/ccds/exon/appris/tsl/hgnc_id/
// gene_symbol/protein_id/polyphen_score+prediction/sift_score+prediction.
// Blocked: `json_sink::annotate_to_json(everything=true)` + full per-VFOA
// field expansion (~25 keys per transcript_consequence).
// See `future-work-vepyr.md`::`Module: json_sink` + `json_sink::format_vfoa_as_json`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_first_line_everything_mode() {
//     let cfg = AnnotateVcfConfig { everything: true, ..Default::default() };
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     // verified via VEP 115 --json --everything at sink-land time:
//     assert_eq!(obj["variant_class"], "SNV");
//     assert!(obj["colocated_variants"].is_array());
//     let tc0 = &obj["transcript_consequences"][0];
//     assert!(tc0["hgvsc"].is_string());
//     assert!(tc0["sift_prediction"].is_string());
// }

// SUBTEST #15 (L456-457): minimal+pick_allele expanded count == 2, first
// allele_string == 'C/T'. `minimal=1` expands CAGAAGAAAG/TAGAAGAAAG,C into
// two CSQ groups.
// Blocked: minimal-mode expansion (shared with OutputFactory.md cluster D).
// See `future-work-vepyr.md`::`Module: json_sink` + `rejoin_variants_in_input_buffer`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_minimal_pick_allele_expanded_count() {
//     let cfg = AnnotateVcfConfig { minimal: true, pick_allele: true, ..Default::default() };
//     let lines = json_sink::annotate_to_json(minimal_input_vcf(), cache_path(), &cfg).await.unwrap();
//     // verified via VEP 115 --json --minimal --pick_allele:
//     assert_eq!(lines.len(), 2);
//     let obj0: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     assert_eq!(obj0["allele_string"], "C/T");
// }

// SUBTEST #16 (L461-462): minimal rejoined count == 1, rejoined
// allele_string == 'CAGAAGAAAG/TAGAAGAAAG/C'.
// `rejoin_variants_in_InputBuffer` collapses the two CSQ groups back into
// one VCF row with original allele_string.
// Blocked: `json_sink::rejoin_variants_in_input_buffer` (same as #15).
// See `future-work-vepyr.md`::`Module: json_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_minimal_rejoined_to_original() {
//     let cfg = AnnotateVcfConfig { minimal: true, ..Default::default() };
//     let lines = json_sink::annotate_to_json(minimal_input_vcf(), cache_path(), &cfg).await.unwrap();
//     assert_eq!(lines.len(), 1);
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     assert_eq!(obj["allele_string"], "CAGAAGAAAG/TAGAAGAAAG/C");
// }

// SUBTEST #17 (L464-508): minimal rejoined transcript_consequences mapped
// by variant_allele key (`'-'` for the deletion ALT, `'T'` for the SNV ALT);
// each carries allele_num + codons/amino_acids/consequence_terms.
// Blocked: `json_sink::annotate_to_json` minimal-mode output (same as #15).
// See `future-work-vepyr.md`::`Module: json_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_minimal_consequences_grouped_by_variant_allele() {
//     let cfg = AnnotateVcfConfig { minimal: true, ..Default::default() };
//     let lines = json_sink::annotate_to_json(minimal_input_vcf(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     let tcs = obj["transcript_consequences"].as_array().unwrap();
//     let alleles: std::collections::HashSet<&str> = tcs.iter()
//         .map(|t| t["variant_allele"].as_str().unwrap())
//         .collect();
//     assert!(alleles.contains("T"));
//     assert!(alleles.contains("-"));
// }

// SUBTEST #18 (L521-547): refseq=1 + uploaded_allele=1 + FASTA-backed —
// transcript_consequences[0] adds given_ref='-', used_ref='-' (FASTA-backfilled),
// uploaded_allele='G/GA' (original REF/ALT), refseq_match.
// Blocked: `json_sink::format_vfoa_as_json` with refseq + uploaded_allele
// branches. Cross-refs `AnnotateVcfConfig::uploaded_allele` (future-work
// line ~430) and `lookup_ref` (Parser entry).
// See `future-work-vepyr.md`::`json_sink::format_vfoa_as_json`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_refseq_uploaded_allele_fields() {
//     let cfg = AnnotateVcfConfig {
//         refseq: true, uploaded_allele: true,
//         reference_fasta_path: Some(test_fasta_path()), ..Default::default()
//     };
//     let lines = json_sink::annotate_to_json(refseq_input_vcf(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     let tc0 = &obj["transcript_consequences"][0];
//     // verified via VEP 115 --json --refseq --uploaded_allele:
//     assert_eq!(tc0["given_ref"], "-");
//     assert_eq!(tc0["used_ref"], "-");
//     assert_eq!(tc0["uploaded_allele"], "G/GA");
//     assert!(tc0["refseq_match"].is_string());
// }

// SUBTEST #19 (L556-596): total_length=1 — cdna_end == '1122/1199' format
// (cdna_start/cdna_end/cds_*/protein_* format as `<position>/<total>`).
// Blocked: `json_sink::annotate_to_json` total_length sub-feature
// (technically sink-agnostic; folded into json_sink coverage per audit).
// See `future-work-vepyr.md`::`Module: json_sink`.
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_total_length_position_with_total() {
//     let cfg = AnnotateVcfConfig { total_length: true, ..Default::default() };
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     let tc0 = &obj["transcript_consequences"][0];
//     // verified via VEP 115 --json --total_length:
//     assert_eq!(tc0["cdna_end"], "1122/1199");
// }

// SUBTEST #21 (L616-631, SKIP-block): custom exact-type — `custom_annotations.test`
// = [{name:test1, fields:{FOO:BAR,FILTER:[PASS]}, allele:T}].
// Blocked: DOUBLE-GAP — needs (a) `json_sink::format_custom_as_json` + (b)
// `--custom` annotation source (out of current vepyr scope; A3-deferred).
// See `future-work-vepyr.md`::`Module: json_sink` §"Secondary gaps".
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_custom_exact_type() {
//     let cfg = AnnotateVcfConfig {
//         custom: vec![CustomAnnotation::vcf("test", "file.vcf.gz", "exact", vec!["FOO"])],
//         ..Default::default()
//     };
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     // verified via VEP 115 --json --custom:
//     let custom_test = &obj["custom_annotations"]["test"][0];
//     assert_eq!(custom_test["name"], "test1");
//     assert_eq!(custom_test["fields"]["FOO"], "BAR");
// }

// SUBTEST #22 (L642-659, SKIP-block): custom overlap-type — `custom_annotations.test`
// = [{name:test1},{name:del1},{name:del2}] (multiple entries; no fields).
// Blocked: same DOUBLE-GAP as #21 — `format_custom_as_json` overlap branch
// + `--custom` annotation source.
// See `future-work-vepyr.md`::`Module: json_sink` §"Secondary gaps".
//
// #[tokio::test(flavor = "multi_thread")]
// async fn json_sink_custom_overlap_type() {
//     let cfg = AnnotateVcfConfig {
//         custom: vec![CustomAnnotation::vcf("test", "file.vcf.gz", "overlap", vec![])],
//         ..Default::default()
//     };
//     let lines = json_sink::annotate_to_json(test_vcf_path(), cache_path(), &cfg).await.unwrap();
//     let obj: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
//     // verified via VEP 115 --json --custom (overlap type):
//     let entries = obj["custom_annotations"]["test"].as_array().unwrap();
//     assert!(entries.len() >= 3);
//     assert_eq!(entries[0]["name"], "test1");
// }

// ──────────────────────────────────────────────────────────────────────────
// End of port_output_factory_json.rs.
//
// Total substantive subtests covered: 16 (rows #1-#5 and #20 are `use_ok`
// boilerplate omitted from coverage parity).
//   - 16 blocked-future-work (rows #6-#19, #21, #22).
//   - 0 architectural-no-analogue.
//   - 0 unit-port / integration-port / e2e-port.
//
// Coverage parity: 0 / 16 = 0% — justified in
// `detailed_plans/OutputFactory_JSON.md §"Coverage parity"`.
// Tier 5 paperwork-only port; all 16 stubs promote simultaneously when the
// `json_sink` module lands.
// ──────────────────────────────────────────────────────────────────────────
