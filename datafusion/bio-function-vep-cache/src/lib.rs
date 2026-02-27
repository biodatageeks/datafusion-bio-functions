//! fjall KV cache backend for VEP variant lookup.
//!
//! Stores one zstd-compressed entry per genomic position in a fjall LSM-tree.
//! Query variants are annotated by direct position key lookups.
//!
//! Cache creation API:
//! - [`CacheLoader`] builds a Fjall cache from a registered source table
//!   (typically Parquet-backed).
//! - Tune build behavior with [`CacheLoader::with_parallelism`],
//!   [`CacheLoader::with_zstd_level`], and [`CacheLoader::with_dict_size_kb`].
//!
//! Minimal example:
//! ```no_run
//! # use datafusion::prelude::SessionContext;
//! # use datafusion_bio_function_vep_cache::CacheLoader;
//! # async fn demo() -> datafusion::common::Result<()> {
//! let ctx = SessionContext::new();
//! ctx.register_parquet("vep_source", "/path/to/variation.parquet", Default::default()).await?;
//!
//! let _stats = CacheLoader::new("vep_source", "/path/to/variation_fjall")
//!     .with_parallelism(8)
//!     .with_zstd_level(9)
//!     .with_dict_size_kb(256)
//!     .load(&ctx)
//!     .await?;
//! # Ok(())
//! # }
//! ```

pub mod allele_index;
pub mod cache_exec;
pub mod cache_provider;
pub mod key_encoding;
pub mod kv_store;
pub mod loader;
pub mod position_entry;

pub use cache_provider::KvCacheTableProvider;
pub use kv_store::VepKvStore;
pub use loader::CacheLoader;
