//! Parquet cache backend (Phase 2, zero-Parquet migration).
//!
//! Point-lookup variation reader/writer + cache detection, replacing the Parquet
//! backend. Design + measured acceptance targets:
//! `docs/superpowers/plans/2026-07-02-parquet-cache-backend-phase2.md`.
//!
//! Encoding decisions (validated on chr1, 319,349 HG002 variants):
//! - no-dictionary zstd, ~4 KiB / 512-row data pages, footer page index;
//! - binary flags (`failed`/`somatic`/`phenotype_or_disease`) as non-nullable
//!   `Boolean` (clean schema);
//! - AF as struct-of-arrays: `<src>_alleles: List<Utf8>` + `<src>_freqs: List<f32>`
//!   (positional per population), reconstructed to exact CSQ text via `%.4g`;
//! - `variation_name` deduped (null where `== dbsnp_ids`, reconstruct via coalesce).

pub mod detect;
pub mod encode;
pub mod page_dir;
pub mod scan;
pub mod sift;
pub mod variation_lookup;
pub mod write;
