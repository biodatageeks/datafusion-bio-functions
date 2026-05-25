//! v115 port of `ensembl-vep/t/AnnotationSource_Cache_Transcript.t`.
//!
//! Asserts that vepyr's CSQ output on chr21:25–26 Mb (the canonical Perl
//! test region) matches real Ensembl VEP 115 across the transcript-cache
//! hard surface: which transcript a variant overlaps, what consequence
//! it produces, and the nearest-feature annotation when intergenic.
//!
//! Per `porting-tests/detailed_plans/AnnotationSource_Cache_Transcript.md`,
//! coverage parity is achieved against in-scope subtests; refseq /
//! merged-cache / Sereal / `--gencode_basic` / `--transcript_filter`
//! branches are `gap-justified` (out of scope or anti-goal).
//!
//! Hard CSQ fields per `porting-plan.md` Recipe 1:
//!   Consequence, Feature, Gene, BIOTYPE, IMPACT, STRAND, SIFT, PolyPhen,
//!   NEAREST, SYMBOL.
//! Everything else is soft (reported by the harness but non-fatal).

mod port_common;

const HARD_FIELDS: &[&str] = &[
    "Consequence",
    "Feature",
    "Gene",
    "BIOTYPE",
    "IMPACT",
    "STRAND",
    "SIFT",
    "PolyPhen",
    "NEAREST",
    "SYMBOL",
];

#[tokio::test(flavor = "multi_thread")]
async fn port_cache_transcript_csq_matches_golden() {
    let ran = port_common::run_and_compare_csq("cache_transcript", HARD_FIELDS)
        .await
        .unwrap();
    assert!(
        ran,
        "fixtures must be present: \
         vep-benchmark/data/port/_cache115/{{parquet/115_GRCh38_vep, reference.fa}} \
         (build via scripts/port/build_v115_parquet.sh) AND \
         datafusion/bio-function-vep/tests/data/port/cache_transcript/{{input.vcf,golden.vcf}} \
         (generate golden via scripts/port/gen_golden.sh cache_transcript)"
    );
}
