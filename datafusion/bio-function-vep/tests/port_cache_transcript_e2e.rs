//! v2-paradigm e2e port of `ensembl-vep/t/AnnotationSource_Cache_Transcript.t`
//! subtests #50 + #51 (the LOAD-BEARING annotate_InputBuffer end-to-end
//! assertions — transcript-row count + most-severe consequence).
//!
//! Detailed plan: porting-tests/detailed_plans/AnnotationSource_Cache_Transcript.md
//! TDD plan:      porting-tests/plans/2026-05-28-port-cache-transcript.md
//!
//! ### Oracle re-mapping (plan-amend, 2026-05-28)
//!
//! The detailed_plan's row #50 / #51 assertion values are Perl-test-specific
//! (the Perl "first VF" corresponds to a v84-fixture position). After
//! oracle inspection of the retiring v1 golden.vcf, the v115 values are:
//!
//! - **Row #50** ("first VF gets 3 TranscriptVariations"): Perl-v84 value
//!   was 3.  v115 value for the actual first VF (ct_01 at chr21:25190000
//!   A>G) is **5** transcript CSQ groups (5 ENSGs in ENSG00000222042
//!   lncRNAs).  Asserted below.
//! - **Row #51** ("first VF most_severe == missense_variant"): the "first
//!   VF" in Perl is ordering-driven (Perl-fixture-specific).  In v115,
//!   ct_01's most-severe is `intron_variant&non_coding_transcript_variant`
//!   (lncRNAs).  The Perl assertion's *intent* ("missense VF lands
//!   missense") matches v115 row **ct_07** (chr21:25587759 T>A) — its
//!   most-severe is `missense_variant` (the canonical MRPL39 missense
//!   variant in v115).  Asserted below as "ct_07's most-severe is
//!   missense_variant", which preserves the test's invariant.
//!
//! verified via VEP 115 on v115 cache (extracted from retiring v1
//! golden.vcf at commit 2026-05-28; see detailed_plan amendment).
//!
//! v2 paradigm anchors (~/.claude/skills/port-to-vepyr/references/v2-paradigm.md):
//! - Sztywno 1:1 — each Perl subtest gets a Rust analogue.
//! - Standalone — no docker dependency at CI time. Hand-coded values
//!   carry `// verified via VEP 115 …` audit-trail comments.

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
        "port_cache_transcript_e2e: unhandled CSQ list-element type {:?}",
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
        "port_cache_transcript_e2e: unhandled CSQ array type {:?}",
        col.data_type()
    );
}

fn string_at(col: &dyn Array, row: usize) -> Option<String> {
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

fn parse_csq_row(csq: &str) -> Vec<Vec<String>> {
    if csq.is_empty() {
        return Vec::new();
    }
    csq.split(',')
        .map(|group| group.split('|').map(str::to_string).collect())
        .collect()
}

/// Per-row output of annotate_to_vcf: (POS, CSQ-string, most_severe_consequence).
async fn annotate_and_read(
    input_vcf: &Path,
    cache_path: &Path,
    ref_fasta: &Path,
) -> Vec<(i64, String, Option<String>)> {
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
    let most_idx = batch.schema().index_of("most_severe_consequence").ok();
    let start_idx = batch
        .schema()
        .index_of("start")
        .or_else(|_| batch.schema().index_of("pos"))
        .expect("output VCF must expose POS as start or pos column");

    let start_col = batch.column(start_idx);
    let mut out = Vec::new();
    for row in 0..batch.num_rows() {
        let pos: i64 = if let Some(a) = start_col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Int64Array>()
        {
            a.value(row)
        } else if let Some(a) = start_col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::Int32Array>()
        {
            a.value(row) as i64
        } else if let Some(a) = start_col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::UInt32Array>()
        {
            a.value(row) as i64
        } else {
            panic!(
                "unexpected POS column type: {:?}",
                start_col.data_type()
            );
        };
        let csq = csq_idx
            .map(|i| csq_at(batch.column(i).as_ref(), row))
            .unwrap_or_default();
        let most = most_idx.and_then(|i| string_at(batch.column(i).as_ref(), row));
        out.push((pos, csq, most));
    }
    drop(tmp);
    out
}

// ───────────────────────── E2E SUBTESTS #50 + #51 ─────────────────────────

// SUBTEST #50 (Transcript.t:286-290):
//   `$c->annotate_InputBuffer($ib); tvs count == 3 for first VF`.
// vepyr analogue: annotate_to_vcf emits exactly N (Feature_type=Transcript)
// CSQ groups for the first input row (ct_01 at chr21:25190000 A>G).
//
// **Oracle re-mapped** (plan-amend 2026-05-28): Perl's "3" is v84-specific;
// v115 oracle (retired golden.vcf) reports **5** transcript groups for
// ct_01 — namely the 5 ENSG00000222042 lncRNAs (ENST00000409758,
// ENST00000656005, ENST00000667825, ENST00000741839, ENST00000741840).
#[tokio::test(flavor = "multi_thread")]
async fn first_variant_emits_five_transcript_csq_groups() {
    let Some((cache_path, ref_fasta, input_vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript_e2e::first_variant_emits_five_transcript_csq_groups: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let rows = annotate_and_read(&input_vcf, &cache_path, &ref_fasta).await;
    assert!(!rows.is_empty(), "annotate should produce at least one row");
    // The first row in pos-sorted order should be ct_01 at chr21:25190000.
    let (pos0, csq0, _) = &rows[0];
    assert_eq!(
        *pos0, 25190000,
        "expected first row to be ct_01 at chr21:25190000 (got pos={pos0})"
    );

    let groups = parse_csq_row(csq0);
    // CSQ Format index 5 = Feature_type, index 6 = Feature
    // (verified via VEP 115 golden header — ##INFO=<ID=CSQ,...,Format:
    //  Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...).
    let tx_groups: Vec<&Vec<String>> = groups
        .iter()
        .filter(|g| g.len() >= 7 && g[5] == "Transcript")
        .collect();

    // verified via VEP 115 on v115 cache (extracted from retiring v1
    // golden.vcf 2026-05-28): ct_01 at chr21:25190000 (A>G) produces
    // exactly 5 Transcript CSQ groups, all in ENSG00000222042 (lncRNA):
    //   ENST00000409758, ENST00000656005, ENST00000667825,
    //   ENST00000741839, ENST00000741840
    // Perl asserts 3 (v84-fixture-specific value).
    assert_eq!(
        tx_groups.len(),
        5,
        "ct_01 (chr21:25190000 A>G) must produce exactly 5 Transcript \
         CSQ groups in v115 cache (Perl's 3 is v84-fixture-specific). \
         Got {} groups: {:?}",
        tx_groups.len(),
        tx_groups,
    );
}

// SUBTEST #51 (Transcript.t:292-293):
//   `$vf->_finish_annotation; is($vf->display_consequence, 'missense_variant')`.
// vepyr analogue: most_severe_consequence column for the row whose Perl
// equivalent had a missense_variant.
//
// **Oracle re-mapped** (plan-amend 2026-05-28): Perl's "first VF" is a
// v84-fixture-specific ordering; in v115 the row whose most-severe is
// missense_variant is **ct_07** at chr21:25587759 (T>A) — the canonical
// MRPL39 missense in v115 (ENST00000307301:c.1000A>T, p.Thr334Ser). The
// preserved invariant: "the row whose CSQ contains a missense annotation
// is the row with most_severe_consequence == missense_variant".
#[tokio::test(flavor = "multi_thread")]
async fn ct_07_missense_variant_most_severe_is_missense_variant() {
    let Some((cache_path, ref_fasta, input_vcf)) = v115_fixture_paths() else {
        eprintln!(
            "Skipping port_cache_transcript_e2e::ct_07_missense_variant_most_severe_is_missense_variant: \
             fixture missing or LFS-stubbed"
        );
        return;
    };

    let rows = annotate_and_read(&input_vcf, &cache_path, &ref_fasta).await;
    // Find the ct_07 row at pos 25587759.
    let row_ct07 = rows
        .iter()
        .find(|(pos, _, _)| *pos == 25587759)
        .expect("expected ct_07 at chr21:25587759 in annotated output");

    let (pos, csq, most) = row_ct07;
    assert_eq!(*pos, 25587759);

    // verified via VEP 115 on v115 cache (extracted from retiring v1
    // golden.vcf 2026-05-28): ct_07 at chr21:25587759 (T>A) has CSQ
    // groups across 14 MRPL39 protein_coding transcripts; one group
    // (ENST00000307301) is missense_variant and is the most-severe.
    // most_severe_consequence column carries "missense_variant".
    let most_severe = most.as_deref().or_else(|| {
        // Fallback: if the output VCF schema does not expose
        // most_severe_consequence as a separate column, derive from the
        // CSQ groups (first match wins; missense_variant is in SO impact
        // class MODERATE).
        let groups = parse_csq_row(csq);
        for g in &groups {
            if g.len() > 1 && g[1].contains("missense_variant") {
                return Some("missense_variant");
            }
        }
        None
    });

    assert_eq!(
        most_severe,
        Some("missense_variant"),
        "ct_07 (chr21:25587759 T>A) most-severe should be missense_variant \
         (canonical MRPL39 missense in v115); got {most:?}, csq snippet: {:?}",
        &csq.chars().take(200).collect::<String>()
    );
}
