//! Fjall KV cache backend for VEP variant lookup.
//!
//! Stores one zstd-compressed entry per genomic position in a fjall LSM-tree.
//! Query variants are annotated by direct position key lookups.

pub mod allele_index;
pub mod cache_exec;
pub mod cache_provider;
pub mod context_store;
pub mod key_encoding;
pub mod kv_store;
pub mod loader;
pub mod position_entry;
pub mod sift_store;

pub use cache_provider::KvCacheTableProvider;
pub use context_store::{ExonKvStore, TranslationKvStore};
pub use kv_store::VepKvStore;
pub use loader::CacheLoader;
pub use sift_store::SiftKvStore;
