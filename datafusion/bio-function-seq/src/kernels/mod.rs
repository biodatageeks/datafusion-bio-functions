//! Pure sequence kernels over byte slices.
//!
//! Every function here works on `&[u8]` and writes into a caller-owned buffer,
//! with no Arrow, DataFusion or Polars types. The DataFusion functions in
//! [`crate::udf`] and the polars-bio expression plugin call the same code, so
//! the two cannot drift apart.
//!
//! The semantics follow the oracle tools (seqkit, seqfu, seqtk, Biopython);
//! polars-bio's `add-sequence-operations` change records the probes behind
//! each rule.

mod codon_tables;
mod complement;
mod gc;
mod hash;
mod metrics;
mod motif;
mod orf;
mod quality;
mod transcribe;
mod translate;
mod trim;
mod weight;
mod weights_data;

use std::fmt;

pub use complement::{Alphabet, complement_into, revcomp_into, reverse_into};
pub use gc::{
    GcCounts, GcDenominator, Windows, gc_content, gc_content_windows_into, gc_counts, gc_skew_into,
};
pub use hash::seq_hash;
pub use metrics::{ReadMetrics, read_metrics};
pub use motif::Motif;
pub use orf::{Orf, OrfFinder, StartCodons, Strands};
pub use quality::{QualityMethod, expected_errors, phred_mean, phred_min_max};
pub use transcribe::{back_transcribe_into, transcribe_into};
pub use translate::{Ambiguity, CdsError, Translator, codon_table_ids};
pub use trim::QualityTrimmer;
pub use weight::MolecularWeight;

/// An input the kernels reject.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KernelError {
    /// A character outside the IUPAC nucleotide alphabet, raised only in
    /// strict mode. `position` is a byte offset; a byte that does not start a
    /// valid UTF-8 character is reported as U+FFFD.
    InvalidBase { symbol: char, position: usize },
    /// A quality byte below the Phred offset.
    QualityBelowOffset {
        byte: u8,
        offset: u8,
        position: usize,
    },
    /// A translation frame other than ±1, ±2 or ±3.
    InvalidFrame(i64),
    /// An id that is not an NCBI translation table.
    UnknownCodonTable(i64),
    /// A sequence that is not a complete coding sequence, raised only when
    /// translating with `cds`. `detail` names the failed check and the codon.
    InvalidCds { error: CdsError, detail: String },
    /// A sequence and its quality string of different lengths.
    LengthMismatch { sequence: usize, quality: usize },
    /// An option string the kernel does not recognise.
    InvalidOption {
        name: &'static str,
        value: String,
        expected: &'static str,
    },
}

/// Formats a byte as `'Z'` when printable ASCII, otherwise as `0xNN`.
struct ShowByte(u8);

impl fmt::Display for ShowByte {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_ascii_graphic() {
            write!(f, "'{}'", self.0 as char)
        } else {
            write!(f, "0x{:02X}", self.0)
        }
    }
}

/// Formats a character as `'Z'` or `'é'` when printable, otherwise as `U+NNNN`.
struct ShowChar(char);

impl fmt::Display for ShowChar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_control() || self.0.is_whitespace() {
            write!(f, "U+{:04X}", self.0 as u32)
        } else {
            write!(f, "'{}'", self.0)
        }
    }
}

impl fmt::Display for KernelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KernelError::InvalidBase { symbol, position } => write!(
                f,
                "invalid nucleotide {} at byte offset {position}; strict mode accepts only the IUPAC codes ACGTURYSWKMBDHVN",
                ShowChar(*symbol)
            ),
            KernelError::QualityBelowOffset {
                byte,
                offset,
                position,
            } => write!(
                f,
                "quality byte {} at position {position} is below the Phred offset {offset}",
                ShowByte(*byte)
            ),
            KernelError::InvalidFrame(frame) => {
                write!(f, "invalid frame {frame}; expected 1, 2, 3, -1, -2 or -3")
            }
            KernelError::UnknownCodonTable(id) => write!(
                f,
                "unknown NCBI translation table {id}; expected one of {:?}",
                codon_table_ids().collect::<Vec<_>>()
            ),
            KernelError::InvalidCds { detail, .. } => write!(
                f,
                "invalid CDS: {detail}; cds_status reports the failing check of each row"
            ),
            KernelError::LengthMismatch { sequence, quality } => write!(
                f,
                "sequence length {sequence} differs from quality length {quality}"
            ),
            KernelError::InvalidOption {
                name,
                value,
                expected,
            } => write!(f, "invalid {name} '{value}'; expected {expected}"),
        }
    }
}

impl std::error::Error for KernelError {}
