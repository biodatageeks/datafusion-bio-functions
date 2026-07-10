//! Port of the per-chrom cache-layout subtests from
//! `ensembl-vep/t/AnnotationSource_Cache_{RegFeat,Variation}.t`.
//!
//! Re-ported (Stage 0 of the port-validation pipeline, 2026-07-08) onto
//! master v0.14.0's manifest-based `PartitionedParquetCache`
//! (`parquet_cache::detect`). The pre-v0.14.0 API this originally targeted —
//! bare-`.parquet`-dir `detect` + `has_chrom` + a `CONTEXT_TYPES` const — was
//! replaced by a `chrom_manifest.json`-driven resolver. The two transformation
//! rules applied here:
//!   1. Fixtures write a `chrom_manifest.json` per entity (via
//!      `ChromManifest::new(...).write_to_entity_dir`); the old bare-file
//!      fixtures no longer satisfy `detect`, which requires a readable
//!      `variation/chrom_manifest.json`.
//!   2. `cache.has_chrom(ctx, chrom)` → `cache.context_path(ctx, chrom).is_some()`
//!      (and the negated form → `.is_none()`); `has_chrom` was removed.
//!
//! Standalone: no VEP runtime dependency, no external cache — every fixture is
//! a `tempfile::TempDir` built in-test.
//!
//! ## Flatten-order subtests (RegFeat.t #20/#21)
//!
//! Perl `AnnotationSource_Cache_RegFeat.t` subtests #20/#21 assert the order of
//! features emitted by `deserialized_obj_to_features` (`RegulatoryFeature`
//! first, `MotifFeature` last). The original quarantined port expressed this via
//! a `CONTEXT_TYPES` const that v0.14.0 removed. That ordering did not disappear
//! — it moved to the annotation layer as the public `FeatureType::rank`
//! (`transcript_consequence.rs`: `RegulatoryFeature=1`, `MotifFeature=2`), whose
//! doc ties it to the same Perl concat order
//! (`Transcript → RegulatoryFeature → MotifFeature`). The two subtests are
//! re-ported below against `FeatureType::rank` — still `unit-port` (see
//! detailed_plan rows #20/#21), and a stronger, production-governing analogue
//! than the retired const. No fake const is invented; the released engine is not
//! modified.

use datafusion_bio_function_vep::cache::manifest::{ChromDatasetEntry, ChromManifest};
use datafusion_bio_function_vep::parquet_cache::detect::PartitionedParquetCache;
use datafusion_bio_function_vep::transcript_consequence::FeatureType;
use std::path::Path;

/// Write one dummy `.parquet` shard plus a `chrom_manifest.json` for
/// `entity`/`chrom` under `base` (mirrors the engine's own fixture helper in
/// `parquet_cache::detect`'s unit tests, so `detect`/`context_path` resolve it).
fn write_entity(base: &Path, entity: &str, chrom: &str) {
    let dir = base.join(entity);
    std::fs::create_dir_all(&dir).unwrap();
    let file = format!("{chrom}.parquet");
    std::fs::write(dir.join(&file), b"parquet-shard-placeholder").unwrap();
    ChromManifest::new(vec![ChromDatasetEntry::new(chrom, file.as_str(), 1)])
        .write_to_entity_dir(&dir)
        .unwrap();
}

// -----------------------------------------------------------------------------
// Port of `ensembl-vep/t/AnnotationSource_Cache_RegFeat.t`.
// See `porting-tests/detailed_plans/AnnotationSource_Cache_RegFeat.md`.
//
// Perl `Cache::RegFeat` accessors (`serializer_type`, `file_suffix`,
// `get_dump_file_name`) probe the v84-era per-1Mb-block storable/gz layout.
// vepyr's per-chrom parquet layout collapses those concepts: the cache layout IS
// parquet (no per-block serializer dispatch), the file suffix IS `.parquet`, and
// the path-resolver is per-chrom (not per-(chrom,1Mb-region)). The unit-ports
// below assert the vepyr equivalents per the v2-paradigm "different-named
// equivalent → still unit-port" rule.
// -----------------------------------------------------------------------------

/// RegFeat fixture: `variation` (required by `detect`) + `regulatory` + `motif`,
/// chr21 only.
fn tmp_regfeat_layout() -> (tempfile::TempDir, PartitionedParquetCache) {
    let tmp = tempfile::tempdir().unwrap();
    for entity in ["variation", "regulatory", "motif"] {
        write_entity(tmp.path(), entity, "chr21");
    }
    let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap())
        .expect("detect should succeed for manifest fixture layout");
    (tmp, cache)
}

// Subtest #5 (RegFeat.t:51): `is($c->serializer_type, 'storable')`.
// vepyr's serializer-equivalent is parquet; the layout resolves through
// `.parquet` files (v2 Rule 6).
#[test]
fn regfeat_serializer_label_is_parquet_for_v115_layout() {
    let (_tmp, cache) = tmp_regfeat_layout();
    let path = cache
        .context_path("regulatory", "chr21")
        .expect("regulatory chr21 should resolve");
    assert!(
        path.extension().and_then(|e| e.to_str()) == Some("parquet"),
        "parquet IS vepyr's serializer-equivalent; resolver should produce \
         a .parquet path (got {path:?})",
    );
}

// Subtest #6 (RegFeat.t:52): `is($c->file_suffix, 'gz')`.
// vepyr's per-chrom file suffix is `.parquet`.
#[test]
fn regfeat_context_path_uses_parquet_suffix() {
    let (_tmp, cache) = tmp_regfeat_layout();
    for chrom in ["chr21"] {
        let path = cache
            .context_path("regulatory", chrom)
            .unwrap_or_else(|| panic!("regulatory/{chrom} should resolve"));
        let s = path.to_string_lossy();
        assert!(
            s.ends_with(".parquet"),
            "regulatory/{chrom} resolved path should end in .parquet (got {s})",
        );
    }
}

// Subtests #7 + #8 (RegFeat.t:54-55): get_dump_file_name → per-region path
// scheme. vepyr replaces the per-1Mb-block layout with per-chrom files; both
// Perl arities (chrom+region OR chrom+start+end) resolve to the same per-chrom
// `.parquet` path.
#[test]
fn regfeat_context_path_for_regulatory_returns_per_chrom_parquet_path() {
    let (tmp, cache) = tmp_regfeat_layout();
    let path = cache
        .context_path("regulatory", "chr21")
        .expect("regulatory chr21 should resolve");
    let expected = tmp.path().join("regulatory").join("chr21.parquet");
    assert_eq!(path, expected);
}

// Subtest #9 (RegFeat.t:57): `throws_ok { $c->get_dump_file_name() }
// qr/No chromosome/`. vepyr's resolver returns `None` when the requested chrom
// file is missing (the "no chromosome" condition).
#[test]
fn regfeat_context_path_with_empty_chrom_returns_none() {
    let (_tmp, cache) = tmp_regfeat_layout();
    // Empty chrom resolves to "<base>/regulatory/.parquet" which does not exist.
    assert!(cache.context_path("regulatory", "").is_none());
    // Missing chrom for an existing context_type → None.
    assert!(cache.context_path("regulatory", "chr99").is_none());
}

// Subtest #16 (RegFeat.t:80-84):
// `is_deeply([sort keys %{$obj->{chr}}], ['MotifFeature', 'RegulatoryFeature'])`.
// vepyr's layout exposes both `regulatory/` and `motif/` subdirs (rule 2:
// `has_chrom` → `context_path(...).is_some()`).
#[test]
fn regfeat_layout_exposes_both_regulatory_and_motif_subdirs() {
    let (_tmp, cache) = tmp_regfeat_layout();
    assert!(
        cache.context_path("regulatory", "chr21").is_some(),
        "regulatory/chr21.parquet should be detected"
    );
    assert!(
        cache.context_path("motif", "chr21").is_some(),
        "motif/chr21.parquet should be detected"
    );
}

// Subtest #34 (RegFeat.t:135-139):
//   is_deeply($c->get_all_regions_by_InputBuffer($ib), [[21, 25]])
// Perl computes (chrom, 1-Mb-block) tuples for a buffer in chr21:25-26 Mb.
// vepyr's per-chrom layout drops the 1-Mb-block dimension; the "region tuple"
// reduces to a single chrom selector (rule 2).
#[test]
fn regfeat_region_tuple_for_chr21_resolves_to_single_chrom_file() {
    let (_tmp, cache) = tmp_regfeat_layout();
    assert!(
        cache.context_path("regulatory", "chr21").is_some(),
        "buffer on chr21 → regulatory/chr21.parquet"
    );
    // No spurious resolution to other chroms.
    for other in ["chr22", "chrX", "chrMT"] {
        assert!(
            cache.context_path("regulatory", other).is_none(),
            "buffer on chr21 must not falsely resolve regulatory/{other}.parquet",
        );
    }
}

// -----------------------------------------------------------------------------
// Port of `ensembl-vep/t/AnnotationSource_Cache_Variation.t` per-chrom resolver.
// See `porting-tests/detailed_plans/AnnotationSource_Cache_Variation.md`.
//
// Perl `Cache::Variation::get_dump_file_name` returned `<dir>/<chr>/<region>_var.gz`
// — a per-(chrom, 1-Mb region) layout. vepyr's
// `context_path("variation", chrom)` produces `<base>/variation/<chrom>.parquet`
// — per-chrom only. The region-arg (1-Mb block) dimension has no vepyr analogue
// and is documented `architectural-no-analogue` in the detailed_plan.
// -----------------------------------------------------------------------------

/// Variation fixture: `variation` entity only, chr21.
fn tmp_variation_layout() -> (tempfile::TempDir, PartitionedParquetCache) {
    let tmp = tempfile::tempdir().unwrap();
    write_entity(tmp.path(), "variation", "chr21");
    let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap())
        .expect("detect should succeed for variation fixture layout");
    (tmp, cache)
}

// Subtest #18 (Variation.t:186): `get_dump_file_name(1, '1-100')` →
// `<dir>/1/1-100_var.gz`. vepyr's analogue is per-chrom-only —
// `<base>/variation/<chrom>.parquet`.
#[test]
fn variation_context_path_returns_per_chrom_parquet_path() {
    let (tmp, cache) = tmp_variation_layout();
    let path = cache
        .context_path("variation", "chr21")
        .expect("variation chr21 should resolve");
    let expected = tmp.path().join("variation").join("chr21.parquet");
    assert_eq!(path, expected);
}

// Subtest #19 (Variation.t:187): `get_dump_file_name(1, 1, 100)` — same path as
// #18 via alternate Perl arity. vepyr's resolver has a single arity
// (`context_type`, `chrom`); assert the resolved path uses the `.parquet` suffix.
#[test]
fn variation_context_path_uses_parquet_suffix() {
    let (_tmp, cache) = tmp_variation_layout();
    let path = cache
        .context_path("variation", "chr21")
        .expect("variation chr21 should resolve");
    let s = path.to_string_lossy();
    assert!(
        s.ends_with(".parquet"),
        "variation/chr21 resolved path should end in .parquet (got {s})",
    );
}

// Subtest #20 (Variation.t:189): `throws_ok { get_dump_file_name() }
// qr/No chromosome/`. vepyr's resolver returns `None` for an empty chrom (no
// file at `<base>/variation/.parquet`) instead of throwing — the equivalent
// observable behaviour.
#[test]
fn variation_context_path_with_empty_chrom_returns_none() {
    let (_tmp, cache) = tmp_variation_layout();
    // Empty chrom → "<base>/variation/.parquet" — does not exist.
    assert!(cache.context_path("variation", "").is_none());
    // Missing chrom for an existing context_type → None.
    assert!(cache.context_path("variation", "chr99").is_none());
}

// Subtest #32 (Variation.t:253-257):
//   is_deeply($c->get_all_regions_by_InputBuffer($ib), [[21, 25]])
// Perl computes (chrom, 1-Mb-block) tuples for a buffer in chr21:25-26 Mb.
// vepyr's per-chrom layout drops the 1-Mb-block dimension; the region tuple
// reduces to a single chrom selector (rule 2).
#[test]
fn variation_region_tuple_for_chr21_resolves_to_single_chrom_file() {
    let (_tmp, cache) = tmp_variation_layout();
    assert!(
        cache.context_path("variation", "chr21").is_some(),
        "buffer on chr21 → variation/chr21.parquet"
    );
    for other in ["chr22", "chrX", "chrMT"] {
        assert!(
            cache.context_path("variation", other).is_none(),
            "buffer on chr21 must not falsely resolve variation/{other}.parquet",
        );
    }
}

// -----------------------------------------------------------------------------
// Flatten-order subtests (RegFeat.t #20/#21). Perl asserts the order in which
// regulatory-family features come out of `deserialized_obj_to_features`
// (RegulatoryFeature first, MotifFeature last). In v0.14.0 that ordering is the
// public `FeatureType::rank` (Transcript=0, RegulatoryFeature=1, MotifFeature=2),
// which governs CSQ-row output order. Pure-const asserts — no fixture needed.
// -----------------------------------------------------------------------------

// Subtest #20 (RegFeat.t:100): `ref($features->[0]) == RegulatoryFeature` —
// first in the flattened list is a RegulatoryFeature. vepyr: RegulatoryFeature
// outranks (precedes) MotifFeature in `FeatureType::rank`.
#[test]
fn regfeat_regulatory_ranks_before_motif() {
    assert!(
        FeatureType::RegulatoryFeature.rank() < FeatureType::MotifFeature.rank(),
        "RegulatoryFeature (rank {}) must precede MotifFeature (rank {}) — \
         matches Perl reverse-sort-keys flatten order (RegulatoryFeature first)",
        FeatureType::RegulatoryFeature.rank(),
        FeatureType::MotifFeature.rank(),
    );
}

// Subtest #21 (RegFeat.t:101): `ref($features->[-1]) == MotifFeature` — last in
// the flattened list is a MotifFeature. Among the regulatory-family pair,
// MotifFeature carries the trailing rank.
#[test]
fn regfeat_motif_is_last_regulatory_family_feature() {
    assert!(
        FeatureType::MotifFeature.rank() > FeatureType::RegulatoryFeature.rank(),
        "MotifFeature (rank {}) must be the trailing regulatory-family feature, \
         after RegulatoryFeature (rank {}) — matches Perl $features->[-1] == MotifFeature",
        FeatureType::MotifFeature.rank(),
        FeatureType::RegulatoryFeature.rank(),
    );
}
