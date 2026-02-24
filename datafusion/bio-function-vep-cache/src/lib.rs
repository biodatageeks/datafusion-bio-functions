//! fjall KV cache backend for VEP variant lookup.
//!
//! Provides a window-based cache strategy: group variants into ~1Mb genomic
//! windows, store each window as an Arrow IPC batch in a fjall LSM-tree,
//! and annotate query variants by loading only the windows they overlap.

pub mod allele_index;
pub mod cache_exec;
pub mod cache_provider;
pub mod key_encoding;
pub mod kv_store;
pub mod loader;
pub mod mmap_block_store;
pub mod position_entry;
pub mod position_index;
pub mod window;

pub use cache_provider::KvCacheTableProvider;
pub use kv_store::VepKvStore;
pub use loader::CacheLoader;
pub use position_index::PositionIndex;
