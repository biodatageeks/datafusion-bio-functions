//! Sequence operations for Apache DataFusion, inspired by seqkit and seqfu.
//!
//! - [`kernels`]: pure functions over byte slices (reverse complement,
//!   complement, GC content and skew, translation with CDS validation,
//!   transcription, ORF finding, IUPAC motif search, molecular weight,
//!   hashing, read quality summaries and quality trimming). They have no
//!   Arrow or DataFusion types, so other engines, such as polars-bio's
//!   Polars expression plugin, can call the same code.
//! - [`udf`]: the `seq_*` DataFusion scalar functions built on them.
//!
//! ```
//! # #[tokio::main(flavor = "current_thread")]
//! # async fn main() -> datafusion::error::Result<()> {
//! use datafusion::prelude::SessionContext;
//!
//! let ctx = SessionContext::new();
//! datafusion_bio_function_seq::register_seq_functions(&ctx);
//! let batches = ctx
//!     .sql("SELECT seq_revcomp('ACGTN') AS rc, seq_translate('ATGGCCTAA') AS aa")
//!     .await?
//!     .collect()
//!     .await?;
//! # Ok(())
//! # }
//! ```

pub mod kernels;
pub mod udf;

pub use udf::{register_seq_functions, seq_udfs};
