//! v2-paradigm e2e port of `ensembl-vep/t/AnnotationSource_Cache_RegFeat.t`
//! subtests #42 + #43 (the LOAD-BEARING rs3989369 end-to-end annotation).
//!
//! Detailed plan: porting-tests/detailed_plans/AnnotationSource_Cache_RegFeat.md
//! TDD plan:      porting-tests/plans/2026-05-27-port-cache-regfeat.md
//!
//! Subtest substitution (plan §8): rs3989369 has no v115 cache regulatory
//! overlap at its v115 position (chr21:25606638), so it is replaced with
//! rs1343855353 at chr21:25039632 (T/C) — which overlaps the v115 enhancer
//! ENSR21_B6Z6N (25039631-25040274). This is a value substitution, not a
//! classification change; the subtests remain e2e-port.
//!
//! verified via VEP 115 on v115 cache at commit b97e1a2139170903c54ef9255ebc056798eff52f:
//!   - regulatory/21.parquet → ENSR21_B6Z6N at 25039631-25040274, feature_type=enhancer
//!   - variation/21.parquet  → rs1343855353 at chr21:25039632 (T/C)
//!   - reference.fa          → chr21:25039632 = T

use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::arrow::array::{Array, LargeListArray, ListArray, StringArray, StringViewArray};
use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::vcf_sink;

// ───────────────────────── shared helpers (inlined per v2 rule) ─────────────────────────

fn workspace_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

fn is_lfs_pointer(path: &Path) -> bool {
    std::fs::read_to_string(path)
        .map(|s| s.starts_with("version https://git-lfs.github.com"))
        .unwrap_or(false)
}

fn base_config(ref_fasta: &str) -> vcf_sink::AnnotateVcfConfig {
    vcf_sink::AnnotateVcfConfig {
        everything: true,
        extended_probes: true,
        reference_fasta_path: Some(ref_fasta.to_string()),
        ..Default::default()
    }
}

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
        "port_cache_regfeat_e2e: unhandled CSQ list-element type {:?}",
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
        "port_cache_regfeat_e2e: unhandled CSQ array type {:?}",
        col.data_type()
    );
}

fn most_severe_at(col: &dyn Array, row: usize) -> Option<String> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return Some(a.value(row).to_string());
    }
    if let Some(a) = col.as_any().downcast_ref::<StringViewArray>() {
        return Some(a.value(row).to_string());
    }
    None
}

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

fn parse_csq_row(csq: &str) -> Vec<Vec<String>> {
    if csq.is_empty() {
        return Vec::new();
    }
    csq.split(',')
        .map(|group| group.split('|').map(str::to_string).collect())
        .collect()
}

/// Output of one annotate_to_vcf run: per-row (CSQ string, most_severe).
async fn annotate_and_read(
    input_vcf: &Path,
    cache_path: &Path,
    ref_fasta: &Path,
) -> Vec<(String, Option<String>)> {
    let tmp = tempfile::TempDir::new().unwrap();
    let output_path = tmp.path().join("annotated.vcf");
    let config = base_config(ref_fasta.to_str().unwrap());
    let _rows = vcf_sink::annotate_to_vcf(
        input_vcf.to_str().unwrap(),
        cache_path.to_str().unwrap(),
        "parquet",
        output_path.to_str().unwrap(),
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
        .unwrap();
    let batches = ctx
        .sql("SELECT * FROM output_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    if batches.is_empty() {
        drop(tmp);
        return Vec::new();
    }
    let batch = datafusion::arrow::compute::concat_batches(&batches[0].schema(), &batches)
        .unwrap();
    let csq_idx = batch.schema().index_of("CSQ").ok();
    let mut out = Vec::new();
    for row in 0..batch.num_rows() {
        let csq = csq_idx
            .map(|i| csq_at(batch.column(i).as_ref(), row))
            .unwrap_or_default();
        // Most-severe consequence MAY come from the INFO field
        // "most_severe_consequence" — but in vcf_sink output it is part
        // of CSQ field index 1 (Consequence). We compute it from the
        // first CSQ group's Consequence column for the regulatory row
        // when present, or from any "most_severe_consequence" INFO field
        // if the VCF schema exposes it.
        let most = batch
            .schema()
            .index_of("most_severe_consequence")
            .ok()
            .and_then(|i| most_severe_at(batch.column(i).as_ref(), row));
        out.push((csq, most));
    }
    drop(tmp);
    out
}

// ───────────────────────── E2E SUBTESTS #42 + #43 ─────────────────────────

// Subtest #42 (RegFeat.t:155-163):
//   `$c->annotate_InputBuffer($ib); rs3989369 rfvs count == 1`
// vepyr analogue: annotate_to_vcf emits exactly 1 CSQ allele-group with
// `Feature_type=RegulatoryFeature` for the chr21:25039632 row (the v115
// substitute for rs3989369; see plan §8).
#[tokio::test(flavor = "multi_thread")]
async fn rs_substitute_overlaps_exactly_one_regulatory_feature() {
    let Some((cache_path, ref_fasta, input_vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_regfeat_e2e::rs_substitute_overlaps_exactly_one_regulatory_feature: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let rows = annotate_and_read(&input_vcf, &cache_path, &ref_fasta).await;
    assert_eq!(rows.len(), 1, "input.vcf has exactly one variant");
    let (csq, _) = &rows[0];
    let groups = parse_csq_row(csq);
    // CSQ Format index 5 is Feature_type (see annotate_table_function.rs:3267-3322).
    let reg_groups: Vec<&Vec<String>> = groups
        .iter()
        .filter(|g| g.len() >= 7 && g[5] == "RegulatoryFeature")
        .collect();

    // verified via VEP 115 on v115 cache at commit b97e1a2:
    //   chr21:25039632 (T>C) overlaps exactly ENSR21_B6Z6N (enhancer,
    //   25039631-25040274) — no other v115 regulatory feature spans this
    //   position. Therefore the annotate_to_vcf output must have exactly
    //   1 RegulatoryFeature CSQ group for this row.
    assert_eq!(
        reg_groups.len(),
        1,
        "chr21:25039632 must produce exactly 1 RegulatoryFeature CSQ group \
         (v115 oracle: ENSR21_B6Z6N enhancer); got {} groups: {:?}",
        reg_groups.len(),
        reg_groups,
    );
    assert_eq!(
        reg_groups[0][6], "ENSR21_B6Z6N",
        "Feature column must carry the overlapping v115 enhancer stable_id"
    );
}

// Subtest #43 (RegFeat.t:166):
//   `$vf->_finish_annotation; is($vf->display_consequence,
//    'regulatory_region_variant')`.
// vepyr analogue: the most-severe consequence (or, when missing, the CSQ
// Consequence field of the regulatory row) is `regulatory_region_variant`.
#[tokio::test(flavor = "multi_thread")]
async fn rs_substitute_display_consequence_is_regulatory_region_variant() {
    let Some((cache_path, ref_fasta, input_vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_regfeat_e2e::rs_substitute_display_consequence_is_regulatory_region_variant: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let rows = annotate_and_read(&input_vcf, &cache_path, &ref_fasta).await;
    assert_eq!(rows.len(), 1);
    let (csq, most) = &rows[0];
    let groups = parse_csq_row(csq);

    // verified via VEP 115 on v115 cache at commit b97e1a2:
    //   ENSR21_B6Z6N is an enhancer, the consequence emitted for any
    //   variant overlapping a regulatory_feature is
    //   `regulatory_region_variant` (annotate_provider.rs +
    //   transcript_consequence.rs).
    //
    // The Perl `display_consequence` is the most-severe term. vepyr's
    // CSQ output marks every regulatory-feature group with
    // Consequence=regulatory_region_variant at index 1. If the VCF sink
    // emits a separate `most_severe_consequence` INFO column, use it;
    // otherwise derive from the regulatory group's CSQ Consequence.
    if let Some(m) = most {
        assert_eq!(
            m, "regulatory_region_variant",
            "most_severe_consequence for chr21:25039632 should be \
             regulatory_region_variant; got {m:?}",
        );
    } else {
        let reg_consequence = groups
            .iter()
            .find(|g| g.len() >= 7 && g[5] == "RegulatoryFeature")
            .map(|g| g[1].as_str());
        assert_eq!(
            reg_consequence,
            Some("regulatory_region_variant"),
            "RegulatoryFeature CSQ Consequence must be \
             regulatory_region_variant; got {reg_consequence:?}",
        );
    }
}
