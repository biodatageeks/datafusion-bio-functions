//! Partitioned per-chromosome parquet cache detection and registration.
//!
//! The partitioned layout stores each context type in a subdirectory with
//! per-chromosome parquet files:
//!
//! ```text
//! 115_GRCh38_vep/
//!   variation/chr1.parquet
//!   variation/chr2.parquet
//!   transcript/chr1.parquet
//!   exon/chr1.parquet
//!   translation_core/chr1.parquet
//!   translation_sift/chr1.parquet
//!   regulatory/chr1.parquet
//!   motif/chr1.parquet
//! ```

use std::path::{Path, PathBuf};

use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{ParquetReadOptions, SessionContext};

/// Known context subdirectories in the partitioned cache layout.
pub const CONTEXT_TYPES: &[&str] = &[
    "variation",
    "transcript",
    "exon",
    "translation_core",
    "translation_sift",
    "regulatory",
    "motif",
];

/// Represents a partitioned per-chromosome parquet cache directory.
#[derive(Clone)]
pub struct PartitionedParquetCache {
    base_dir: PathBuf,
    /// Chromosomes available in the variation subdirectory (from filenames).
    available_chroms: Vec<String>,
}

impl PartitionedParquetCache {
    /// Detect a partitioned cache layout at `path`.
    ///
    /// Returns `Some` if `path` is a directory containing a `variation/`
    /// subdirectory with at least one `*.parquet` file.
    pub fn detect(path: &str) -> Option<Self> {
        let base = Path::new(path);
        if !base.is_dir() {
            return None;
        }
        let variation_dir = base.join("variation");
        if !variation_dir.is_dir() {
            return None;
        }

        let mut chroms = Vec::new();
        let entries = std::fs::read_dir(&variation_dir).ok()?;
        for entry in entries.flatten() {
            let name = entry.file_name();
            let name_str = name.to_string_lossy();
            if let Some(stem) = name_str.strip_suffix(".parquet") {
                chroms.push(stem.to_string());
            }
        }

        if chroms.is_empty() {
            return None;
        }

        // Sort chromosomes naturally: chr1, chr2, ..., chr10, ..., chrX, chrY, chrMT
        chroms.sort_by_key(|a| natural_chrom_order(a));

        Some(Self {
            base_dir: base.to_path_buf(),
            available_chroms: chroms,
        })
    }

    /// All chromosomes available in the variation cache.
    pub fn available_chroms(&self) -> &[String] {
        &self.available_chroms
    }

    /// Path to a per-chromosome parquet file for a given context type.
    /// Returns `None` if the file does not exist.
    pub fn context_path(&self, context_type: &str, chrom: &str) -> Option<PathBuf> {
        let path = self
            .base_dir
            .join(context_type)
            .join(format!("{chrom}.parquet"));
        if path.exists() { Some(path) } else { None }
    }

    /// Whether a per-chromosome parquet file exists for a given context type.
    pub fn has_chrom(&self, context_type: &str, chrom: &str) -> bool {
        self.context_path(context_type, chrom).is_some()
    }

    /// Base directory of the cache.
    pub fn base_dir(&self) -> &Path {
        &self.base_dir
    }
}

/// Register a per-chromosome parquet file as an ephemeral DataFusion table.
///
/// Returns the generated table name if registration succeeded, or `None` if
/// the parquet file does not exist for this chrom/context_type combination.
pub async fn register_chrom_parquet(
    session: &SessionContext,
    cache: &PartitionedParquetCache,
    context_type: &str,
    chrom: &str,
) -> Result<Option<String>> {
    let path = match cache.context_path(context_type, chrom) {
        Some(p) => p,
        None => return Ok(None),
    };

    let table_name = ephemeral_table_name(context_type, chrom);

    // Skip if already registered
    if session.table(&table_name).await.is_ok() {
        return Ok(Some(table_name));
    }

    session
        .register_parquet(
            &table_name,
            path.to_str().ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "non-UTF8 path for {context_type}/{chrom}.parquet"
                ))
            })?,
            ParquetReadOptions::default(),
        )
        .await?;

    Ok(Some(table_name))
}

/// Deregister an ephemeral table from the session.
pub async fn deregister_table(session: &SessionContext, name: &str) -> Result<()> {
    // deregister_table returns Option<Arc<dyn TableProvider>>; ignore it.
    let _ = session.deregister_table(name)?;
    Ok(())
}

/// Generate an ephemeral table name for a context_type + chrom combination.
fn ephemeral_table_name(context_type: &str, chrom: &str) -> String {
    format!("__vep_partitioned_{context_type}_{chrom}")
}

/// Natural chromosome ordering: numeric chroms first (sorted numerically),
/// then X, Y, MT/M, then anything else alphabetically.
fn natural_chrom_order(chrom: &str) -> (u8, u32, String) {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(n) = bare.parse::<u32>() {
        (0, n, String::new())
    } else {
        let priority = match bare {
            "X" => 1,
            "Y" => 2,
            "MT" | "M" => 3,
            _ => 4,
        };
        (priority, 0, bare.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_natural_chrom_order() {
        let mut chroms = vec![
            "chrX".to_string(),
            "chr2".to_string(),
            "chr10".to_string(),
            "chr1".to_string(),
            "chrY".to_string(),
            "chrMT".to_string(),
            "chr22".to_string(),
        ];
        chroms.sort_by_key(|a| natural_chrom_order(a));
        assert_eq!(
            chroms,
            vec!["chr1", "chr2", "chr10", "chr22", "chrX", "chrY", "chrMT"]
        );
    }

    #[test]
    fn test_detect_nonexistent_dir() {
        assert!(PartitionedParquetCache::detect("/nonexistent/path").is_none());
    }

    #[test]
    fn test_detect_file_not_dir() {
        // A regular file should not be detected as partitioned cache
        assert!(PartitionedParquetCache::detect("/dev/null").is_none());
    }

    #[test]
    fn test_ephemeral_table_name() {
        assert_eq!(
            ephemeral_table_name("variation", "chr1"),
            "__vep_partitioned_variation_chr1"
        );
    }

    #[test]
    fn test_detect_dir_without_variation() {
        let tmp = tempfile::tempdir().unwrap();
        // Directory exists but no variation/ subdirectory
        assert!(PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).is_none());
    }

    #[test]
    fn test_detect_variation_dir_empty() {
        let tmp = tempfile::tempdir().unwrap();
        std::fs::create_dir(tmp.path().join("variation")).unwrap();
        // variation/ exists but no parquet files
        assert!(PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).is_none());
    }

    #[test]
    fn test_detect_valid_layout() {
        let tmp = tempfile::tempdir().unwrap();
        let var_dir = tmp.path().join("variation");
        std::fs::create_dir(&var_dir).unwrap();
        // Create dummy parquet files (content doesn't matter for detection)
        std::fs::write(var_dir.join("chr1.parquet"), b"dummy").unwrap();
        std::fs::write(var_dir.join("chr22.parquet"), b"dummy").unwrap();
        std::fs::write(var_dir.join("chrX.parquet"), b"dummy").unwrap();

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(cache.available_chroms(), &["chr1", "chr22", "chrX"]);
    }

    #[test]
    fn test_context_path() {
        let tmp = tempfile::tempdir().unwrap();
        let var_dir = tmp.path().join("variation");
        std::fs::create_dir(&var_dir).unwrap();
        std::fs::write(var_dir.join("chr1.parquet"), b"dummy").unwrap();

        let tx_dir = tmp.path().join("transcript");
        std::fs::create_dir(&tx_dir).unwrap();
        std::fs::write(tx_dir.join("chr1.parquet"), b"dummy").unwrap();

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();

        assert!(cache.has_chrom("variation", "chr1"));
        assert!(cache.has_chrom("transcript", "chr1"));
        assert!(!cache.has_chrom("transcript", "chr99"));
        assert!(!cache.has_chrom("nonexistent", "chr1"));
    }

    // ----------------------------------------------------------------
    // Port of `ensembl-vep/t/AnnotationSource_Cache_RegFeat.t` subtests.
    // See `porting-tests/detailed_plans/AnnotationSource_Cache_RegFeat.md`
    // and `porting-tests/plans/2026-05-27-port-cache-regfeat.md`.
    //
    // Perl `Cache::RegFeat` accessors (`serializer_type`, `file_suffix`,
    // `get_dump_file_name`) probe the v84-era per-1Mb-block storable/gz
    // layout. vepyr's per-chrom parquet layout collapses those concepts:
    // the cache layout IS parquet (no per-block serializer dispatch), the
    // file suffix IS `.parquet`, and the path-resolver is per-chrom (not
    // per-(chrom,1Mb-region)). The unit-ports below assert the vepyr
    // equivalents per the v2 paradigm "different-named equivalent → still
    // unit-port" rule (see `references/v2-paradigm.md` Rule 6 / case 6).
    // ----------------------------------------------------------------

    /// Helper: build a tmp parquet layout with regulatory + motif subdirs
    /// holding one per-chrom dummy file each. Returns the resolver.
    fn tmp_regfeat_layout() -> (tempfile::TempDir, PartitionedParquetCache) {
        let tmp = tempfile::tempdir().unwrap();
        // variation/ is required by `detect` (entrypoint for chrom discovery).
        let var_dir = tmp.path().join("variation");
        std::fs::create_dir(&var_dir).unwrap();
        std::fs::write(var_dir.join("chr21.parquet"), b"dummy").unwrap();

        let reg_dir = tmp.path().join("regulatory");
        std::fs::create_dir(&reg_dir).unwrap();
        std::fs::write(reg_dir.join("chr21.parquet"), b"dummy").unwrap();

        let motif_dir = tmp.path().join("motif");
        std::fs::create_dir(&motif_dir).unwrap();
        std::fs::write(motif_dir.join("chr21.parquet"), b"dummy").unwrap();

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap())
            .expect("detect should succeed for fixture layout");
        (tmp, cache)
    }

    // Subtest #5 (RegFeat.t:51): `is($c->serializer_type, 'storable')`.
    // vepyr's serializer-equivalent is parquet; assertion takes the v2
    // form (Rule 6) — assert the layout resolves through `.parquet` files.
    #[test]
    fn regfeat_serializer_label_is_parquet_for_v115_layout() {
        let (_tmp, cache) = tmp_regfeat_layout();
        let path = cache
            .context_path("regulatory", "chr21")
            .expect("regulatory chr21 should resolve");
        assert!(
            path.extension().and_then(|e| e.to_str()) == Some("parquet"),
            "parquet IS vepyr's serializer-equivalent; resolver should produce \
             a .parquet path (got {:?})",
            path,
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

    // Subtests #7 + #8 (RegFeat.t:54-55): get_dump_file_name → per-region
    // path scheme. vepyr replaces the per-1Mb-block layout with per-chrom
    // files; both Perl arities (chrom+region OR chrom+start+end) resolve
    // to the same per-chrom .parquet path.
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
    // qr/No chromosome/`. vepyr's resolver returns `None` when the
    // requested chrom file is missing (the "no chromosome" condition).
    #[test]
    fn regfeat_context_path_with_empty_chrom_returns_none() {
        let (_tmp, cache) = tmp_regfeat_layout();
        // Empty chrom resolves to "<base>/regulatory/.parquet" which does
        // not exist → None.
        assert!(cache.context_path("regulatory", "").is_none());
        // Missing chrom for an existing context_type → None.
        assert!(cache.context_path("regulatory", "chr99").is_none());
    }

    // Subtest #16 (RegFeat.t:80-84):
    // `is_deeply([sort keys %{$obj->{chr}}], ['MotifFeature', 'RegulatoryFeature'])`.
    // vepyr's layout exposes both `regulatory/` and `motif/` subdirs.
    #[test]
    fn regfeat_layout_exposes_both_regulatory_and_motif_subdirs() {
        let (_tmp, cache) = tmp_regfeat_layout();
        assert!(
            cache.has_chrom("regulatory", "chr21"),
            "regulatory/chr21.parquet should be detected"
        );
        assert!(
            cache.has_chrom("motif", "chr21"),
            "motif/chr21.parquet should be detected"
        );
    }

    // Subtests #20 + #21 (RegFeat.t:100-101):
    //   ref($features->[0]) == RegulatoryFeature
    //   ref($features->[-1]) == MotifFeature
    // The Perl test asserts a specific flatten-order coming out of
    // `deserialized_obj_to_features`. vepyr's analogue is the
    // `CONTEXT_TYPES` const which determines the order in which context
    // types are walked (regulatory before motif).
    //
    // NB: the v115 fixture's motif/21.parquet is empty (0 rows), so a
    // record-level "last is motif" assertion is impossible. The strongest
    // hand-codable invariant is the const ordering itself — see plan §8.
    #[test]
    fn regfeat_context_types_const_orders_regulatory_before_motif() {
        let reg_pos = CONTEXT_TYPES
            .iter()
            .position(|t| *t == "regulatory")
            .expect("CONTEXT_TYPES must list 'regulatory'");
        let motif_pos = CONTEXT_TYPES
            .iter()
            .position(|t| *t == "motif")
            .expect("CONTEXT_TYPES must list 'motif'");
        assert!(
            reg_pos < motif_pos,
            "regulatory must appear before motif in CONTEXT_TYPES \
             (reg={reg_pos}, motif={motif_pos}); this preserves Perl's \
             reverse-sort-keys flatten order (RegulatoryFeature first, \
             MotifFeature last)",
        );
    }

    // Subtest #21 also asserts the LAST element of the flattened list is
    // MotifFeature; assert there is no later context_type in CONTEXT_TYPES.
    #[test]
    fn regfeat_motif_is_the_last_regulatory_context_type() {
        // The two regulatory-feature contexts in vepyr are `regulatory`
        // and `motif`. `motif` must be the trailing one (matching
        // Perl `$features->[-1]` == MotifFeature).
        let motif_pos = CONTEXT_TYPES.iter().position(|t| *t == "motif").unwrap();
        let reg_pos = CONTEXT_TYPES
            .iter()
            .position(|t| *t == "regulatory")
            .unwrap();
        // No regulatory-family context_type appears after motif.
        // (Other context_types like translation_core may appear after
        // motif in CONTEXT_TYPES — that's fine; the invariant is that
        // among {regulatory, motif}, motif is last.)
        assert!(motif_pos > reg_pos);
    }

    // Subtest #34 (RegFeat.t:135-139):
    //   is_deeply($c->get_all_regions_by_InputBuffer($ib), [[21, 25]])
    // Perl computes (chrom, 1-Mb-block) tuples for a buffer that holds
    // variants in chr21:25-26 Mb. vepyr's per-chrom layout drops the
    // 1-Mb-block dimension entirely; the "region tuple" reduces to a
    // single chrom selector.
    #[test]
    fn regfeat_region_tuple_for_chr21_resolves_to_single_chrom_file() {
        let (_tmp, cache) = tmp_regfeat_layout();
        assert!(
            cache.has_chrom("regulatory", "chr21"),
            "buffer on chr21 → regulatory/chr21.parquet"
        );
        // No spurious resolution to other chroms.
        for other in ["chr22", "chrX", "chrMT"] {
            assert!(
                !cache.has_chrom("regulatory", other),
                "buffer on chr21 must not falsely resolve regulatory/{other}.parquet",
            );
        }
    }

    // ----------------------------------------------------------------
    // Port of `ensembl-vep/t/AnnotationSource_Cache_Variation.t` subtests
    // for the per-chrom path resolver.
    // See `porting-tests/detailed_plans/AnnotationSource_Cache_Variation.md`
    // and `porting-tests/plans/2026-05-28-port-cache-variation.md`.
    //
    // The Perl `Cache::Variation::get_dump_file_name` returned a path of the
    // shape `<dir>/<chr>/<region>_var.gz` — a per-(chrom, 1-Mb region)
    // layout. vepyr's `PartitionedParquetCache::context_path("variation",
    // chrom)` produces `<base>/variation/<chrom>.parquet` — per-chrom only.
    // Following the v2 paradigm "different-named equivalent → unit-port"
    // rule, the Perl path-shape subtests map onto vepyr's per-chrom
    // resolver. The region-arg dimension (1-Mb block) has no vepyr
    // analogue and is documented `architectural-no-analogue` in the
    // detailed_plan row #21.
    // ----------------------------------------------------------------

    /// Helper: build a tmp parquet layout for a chr21-only variation
    /// fixture. Returns the resolver.
    fn tmp_variation_layout() -> (tempfile::TempDir, PartitionedParquetCache) {
        let tmp = tempfile::tempdir().unwrap();
        let var_dir = tmp.path().join("variation");
        std::fs::create_dir(&var_dir).unwrap();
        std::fs::write(var_dir.join("chr21.parquet"), b"dummy").unwrap();

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

    // Subtest #19 (Variation.t:187): `get_dump_file_name(1, 1, 100)` —
    // same path as #18 via alternate Perl arity. vepyr's resolver has a
    // single arity (`context_type`, `chrom`); the call shape is identical
    // regardless of how Perl would have spelled the region. Assert the
    // resolved path uses the `.parquet` suffix.
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
    // qr/No chromosome/`. vepyr's resolver returns `None` for an empty
    // chrom (no file exists at `<base>/variation/.parquet`) instead of
    // throwing — the equivalent observable behaviour.
    #[test]
    fn variation_context_path_with_empty_chrom_returns_none() {
        let (_tmp, cache) = tmp_variation_layout();
        // Empty chrom → "<base>/variation/.parquet" — does not exist.
        assert!(cache.context_path("variation", "").is_none());
        // Missing chrom for an existing context_type → None.
        assert!(cache.context_path("variation", "chr99").is_none());
    }

    // Subtest #32 (Variation.t:253-257): `is_deeply(
    //   $c->get_all_regions_by_InputBuffer($ib), [[21, 25]])`. Perl
    // computes (chrom, 1-Mb-block) tuples for a buffer in chr21:25-26 Mb.
    // vepyr's per-chrom layout drops the 1-Mb-block dimension; the region
    // tuple reduces to a single chrom selector. Assert `has_chrom` is
    // true for the present chrom and false for absent chroms.
    #[test]
    fn variation_region_tuple_for_chr21_resolves_to_single_chrom_file() {
        let (_tmp, cache) = tmp_variation_layout();
        assert!(
            cache.has_chrom("variation", "chr21"),
            "buffer on chr21 → variation/chr21.parquet"
        );
        for other in ["chr22", "chrX", "chrMT"] {
            assert!(
                !cache.has_chrom("variation", other),
                "buffer on chr21 must not falsely resolve variation/{other}.parquet",
            );
        }
    }
}
