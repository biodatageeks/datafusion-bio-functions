//! Session-level configuration for VEP annotation.
//!
//! Defines the `bio.annotation` namespace. Register with a session to allow
//! SQL `SET` overrides; otherwise the compiled defaults apply.
//!
//! ```sql
//! SET bio.annotation.cache_size_mb = 2048;
//! SET bio.annotation.v5_zstd_level = 9;
//! SET bio.annotation.v5_dict_size_kb = 256;
//! ```
//!
//! Downstream (polars-bio) registers this extension on the session via
//! `SessionConfig::with_option_extension(AnnotationConfig::default())`.
//! Code that consumes the values uses [`resolve`] to read from the session
//! with fallback to defaults when the extension is not registered.

use datafusion::common::extensions_options;
use datafusion::config::ConfigExtension;
use datafusion::prelude::SessionContext;

extensions_options! {
    /// Configuration options for VEP annotation under the `bio.annotation` namespace.
    pub struct AnnotationConfig {
        /// fjall block cache size in MB for KV cache reads (default: 1024).
        ///
        /// Larger values reduce cold-start latency by caching more LSM block
        /// index pages and data blocks in memory.
        pub cache_size_mb: u64, default = 1024

        /// Zstd compression level for V5 cache writes (default: 3).
        ///
        /// Higher levels produce smaller caches at the cost of slower writes.
        /// Decompression speed is constant regardless of level (read throughput
        /// is unaffected). Recommended range: 1-19.
        /// Level 9 with dict_size_kb=256 is a good balance for write-once caches.
        pub v5_zstd_level: u64, default = 3

        /// Zstd dictionary size in KB for V5 cache writes (default: 112).
        ///
        /// The dictionary is trained from the first batch of position entries
        /// and reused for all subsequent entries. Larger dictionaries can improve
        /// compression ratio at the cost of slightly more memory during writes.
        pub v5_dict_size_kb: u64, default = 112
    }
}

impl ConfigExtension for AnnotationConfig {
    const PREFIX: &'static str = "bio.annotation";
}

/// Read annotation config from a session, falling back to defaults if the
/// extension was not registered.
pub fn resolve(session: &SessionContext) -> AnnotationConfig {
    let defaults = AnnotationConfig::default();
    let state = session.state();
    let cfg = state.config_options();
    cfg.extensions
        .get::<AnnotationConfig>()
        .cloned()
        .unwrap_or(defaults)
}
