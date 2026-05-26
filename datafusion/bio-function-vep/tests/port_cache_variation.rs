//! v115 port of `ensembl-vep/t/AnnotationSource_Cache_Variation.t`.
//!
//! Asserts that vepyr's CSQ output matches real Ensembl VEP 115 across the
//! variation-cache hard surface: existing-variant identity, per-population
//! allele frequencies (1000G + gnomADe), MAX_AF/MAX_AF_POPS, clinical
//! significance, phenotype, publications, somatic flag.
//!
//! Per `porting-tests/detailed_plans/AnnotationSource_Cache_Variation.md`,
//! the matched_alleles nastiness 1–4 cases (chr21:8987004, retargeted to
//! v115 cache row rs1402623917 with C/CCGC/CGGC allele_string) are the
//! LOAD-BEARING subtests: they test vepyr's allele-trimmer against
//! multi-base indel edge cases that simpler SNV cases would not catch.
//!
//! Hard CSQ fields are Recipe 1's variation-cache hard set, adjusted for
//! the v115 cache's actual CSQ Format header (verified Task 4 §Step 4.2):
//!   - included: Existing_variation; AF + 5 1000Genomes super-pops;
//!     gnomADe_AF + 9 sub-pops including v115-only REMAINING + MID;
//!     MAX_AF + MAX_AF_POPS; CLIN_SIG, PHENO, PUBMED, SOMATIC.
//!   - dropped vs Recipe 1 (not declared in v115 golden CSQ Format):
//!     * `gnomADe_OTH_AF` — renamed to `gnomADe_REMAINING_AF` in v115.
//!     * `AA`, `EA`        — ESP-era columns removed from v115 cache header.
//! Everything else is soft (reported by the harness but non-fatal).

mod port_common;

const HARD_FIELDS: &[&str] = &[
    "Existing_variation",
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_NFE_AF",
    "gnomADe_SAS_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_MID_AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "PHENO",
    "PUBMED",
    "SOMATIC",
];

#[tokio::test(flavor = "multi_thread")]
async fn port_cache_variation_csq_matches_golden() {
    let ran = port_common::run_and_compare_csq("cache_variation", HARD_FIELDS)
        .await
        .unwrap();
    assert!(
        ran,
        "fixtures must be present: \
         vep-benchmark/data/port/_cache115/{{parquet/115_GRCh38_vep, reference.fa}} \
         (build via scripts/port/build_v115_parquet.sh) AND \
         datafusion/bio-function-vep/tests/data/port/cache_variation/{{input.vcf,golden.vcf}} \
         (generate golden via scripts/port/gen_golden.sh cache_variation)"
    );
}
