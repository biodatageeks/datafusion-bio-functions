//! VCF file writer for processed sample QC output.
//!
//! Writes multi-sample VCF with FORMAT fields GT:GQ:DP:PL:DS.
//! Automatically applies BGZF compression when the output path ends
//! with `.vcf.gz` or `.vcf.bgz`, producing tabix-indexable output.

use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use datafusion::arrow::array::{Array, Float64Array, Int32Array, StringArray, UInt32Array};
use datafusion::common::{DataFusionError, Result};

use super::processing::ProcessedBatch;

/// VCF writer that outputs processed genotype data.
///
/// Compression is chosen automatically from the file extension:
/// - `.vcf` — plain text
/// - `.vcf.gz` / `.vcf.bgz` — BGZF (tabix-indexable blocked gzip)
pub struct VcfWriter {
    inner: WriterInner,
    /// Reusable buffer for entire variant lines. Avoids per-field write_all
    /// syscalls — a single write_all per variant instead of ~10 per sample.
    line_buf: Vec<u8>,
}

enum WriterInner {
    Plain(BufWriter<File>),
    Bgzf(Box<noodles_bgzf::Writer<File>>),
}

impl Write for WriterInner {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            WriterInner::Plain(w) => w.write(buf),
            WriterInner::Bgzf(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            WriterInner::Plain(w) => w.flush(),
            WriterInner::Bgzf(w) => w.flush(),
        }
    }
}

impl VcfWriter {
    /// Create a new VCF writer, writing the header immediately.
    pub fn new(path: &Path, sample_names: &[String]) -> Result<Self> {
        let file = File::create(path).map_err(|e| {
            DataFusionError::Execution(format!("Cannot create output VCF '{path:?}': {e}"))
        })?;

        let compressed = is_compressed_path(path);
        let mut inner = if compressed {
            WriterInner::Bgzf(Box::new(noodles_bgzf::Writer::new(file)))
        } else {
            WriterInner::Plain(BufWriter::with_capacity(1 << 20, file))
        };

        writeln!(inner, "##fileformat=VCFv4.2").map_err(io_err)?;
        writeln!(
            inner,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .map_err(io_err)?;
        writeln!(
            inner,
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Conditional genotype quality\">"
        )
        .map_err(io_err)?;
        writeln!(
            inner,
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
        )
        .map_err(io_err)?;
        writeln!(inner, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods (corrected for hom-ref 0,0,0 sites)\">").map_err(io_err)?;
        writeln!(
            inner,
            "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">"
        )
        .map_err(io_err)?;

        write!(
            inner,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        )
        .map_err(io_err)?;
        for name in sample_names {
            write!(inner, "\t{name}").map_err(io_err)?;
        }
        writeln!(inner).map_err(io_err)?;

        Ok(Self {
            inner,
            line_buf: Vec::with_capacity(128 * 1024),
        })
    }

    /// Write all variants from a processed batch.
    ///
    /// Buffers each entire variant line into `line_buf` before a single
    /// `write_all` call, avoiding ~10 syscalls per sample (fix #5).
    pub fn write_batch(&mut self, batch: &ProcessedBatch) -> Result<()> {
        let chrom = downcast_string(&batch.chrom)?;
        let start = downcast_u32_or_i32(&batch.start)?;
        let id = downcast_string(&batch.id)?;
        let ref_allele = downcast_string(&batch.ref_allele)?;
        let alt = downcast_string(&batch.alt)?;
        let qual = downcast_f64(&batch.qual)?;
        let filter = downcast_string(&batch.filter)?;

        // Scratch buffer for itoa — avoids per-integer format!() allocation
        let mut itoa_buf = itoa::Buffer::new();

        #[allow(clippy::needless_range_loop)]
        for var_idx in 0..batch.n_variants {
            self.line_buf.clear();

            // ---- Fixed columns ----
            let chrom_val = if chrom.is_null(var_idx) {
                "."
            } else {
                chrom.value(var_idx)
            };
            let pos = start[var_idx] + 1;
            let id_val = if id.is_null(var_idx) {
                "."
            } else {
                let v = id.value(var_idx);
                if v.is_empty() { "." } else { v }
            };
            let ref_val = if ref_allele.is_null(var_idx) {
                "."
            } else {
                ref_allele.value(var_idx)
            };
            let alt_val = if alt.is_null(var_idx) {
                "."
            } else {
                alt.value(var_idx)
            };
            let qual_str = if qual.is_null(var_idx) {
                ".".to_string()
            } else {
                format_qual(qual.value(var_idx))
            };
            let filter_val = if filter.is_null(var_idx) {
                "."
            } else {
                let v = filter.value(var_idx);
                if v.is_empty() { "." } else { v }
            };

            self.line_buf.extend_from_slice(chrom_val.as_bytes());
            self.line_buf.push(b'\t');
            self.line_buf
                .extend_from_slice(itoa_buf.format(pos).as_bytes());
            self.line_buf.push(b'\t');
            self.line_buf.extend_from_slice(id_val.as_bytes());
            self.line_buf.push(b'\t');
            self.line_buf.extend_from_slice(ref_val.as_bytes());
            self.line_buf.push(b'\t');
            self.line_buf.extend_from_slice(alt_val.as_bytes());
            self.line_buf.push(b'\t');
            self.line_buf.extend_from_slice(qual_str.as_bytes());
            self.line_buf.push(b'\t');
            self.line_buf.extend_from_slice(filter_val.as_bytes());
            self.line_buf.extend_from_slice(b"\t.\tGT:GQ:DP:PL:DS");

            // ---- Per-sample columns: format into line_buf ----
            let s_start = batch.sample_offsets[var_idx];
            let s_end = batch.sample_offsets[var_idx + 1];
            for s in s_start..s_end {
                self.line_buf.push(b'\t');

                // GT
                let gt = if batch.gt_final.is_null(s) {
                    "./."
                } else {
                    batch.gt_final.value(s)
                };
                self.line_buf.extend_from_slice(gt.as_bytes());

                // :GQ
                self.line_buf.push(b':');
                if batch.gq_values.is_null(s) {
                    self.line_buf.push(b'.');
                } else {
                    self.line_buf
                        .extend_from_slice(itoa_buf.format(batch.gq_values.value(s)).as_bytes());
                }

                // :DP
                self.line_buf.push(b':');
                if batch.dp_values.is_null(s) {
                    self.line_buf.push(b'.');
                } else {
                    self.line_buf
                        .extend_from_slice(itoa_buf.format(batch.dp_values.value(s)).as_bytes());
                }

                // :PL0,PL1,PL2
                self.line_buf.push(b':');
                self.line_buf
                    .extend_from_slice(itoa_buf.format(batch.pl0.value(s)).as_bytes());
                self.line_buf.push(b',');
                self.line_buf
                    .extend_from_slice(itoa_buf.format(batch.pl1.value(s)).as_bytes());
                self.line_buf.push(b',');
                self.line_buf
                    .extend_from_slice(itoa_buf.format(batch.pl2.value(s)).as_bytes());

                // :DS
                self.line_buf.push(b':');
                if batch.ds_values.is_null(s) {
                    self.line_buf.push(b'.');
                } else {
                    // Use write! for float formatting into the byte buffer via a wrapper
                    let before = self.line_buf.len();
                    write!(
                        BufToVec(&mut self.line_buf),
                        "{:.4}",
                        batch.ds_values.value(s)
                    )
                    .map_err(|e| DataFusionError::Execution(format!("fmt: {e}")))?;
                    let _ = before; // suppress unused
                }
            }

            self.line_buf.push(b'\n');

            // ---- Single write_all per variant line ----
            self.inner.write_all(&self.line_buf).map_err(io_err)?;
        }
        Ok(())
    }

    /// Flush and finalize the VCF output.
    pub fn finish(self) -> Result<()> {
        match self.inner {
            WriterInner::Plain(mut w) => w.flush().map_err(io_err),
            WriterInner::Bgzf(mut w) => w.try_finish().map_err(io_err),
        }
    }
}

/// Thin wrapper to use `write!` with a `Vec<u8>`.
struct BufToVec<'a>(&'a mut Vec<u8>);

impl std::fmt::Write for BufToVec<'_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        self.0.extend_from_slice(s.as_bytes());
        Ok(())
    }
}

// ============================================================================
// Helpers
// ============================================================================

fn is_compressed_path(path: &Path) -> bool {
    let s = path.to_string_lossy();
    s.ends_with(".gz") || s.ends_with(".bgz")
}

fn io_err(e: impl std::fmt::Display) -> DataFusionError {
    DataFusionError::Execution(format!("VCF write error: {e}"))
}

fn format_qual(q: f64) -> String {
    if q == q.round() && q.abs() < 1e15 {
        format!("{}", q as i64)
    } else {
        format!("{q:.1}")
    }
}

fn downcast_string(arr: &dyn Array) -> Result<StringArray> {
    if let Some(s) = arr.as_any().downcast_ref::<StringArray>() {
        return Ok(s.clone());
    }
    use datafusion::arrow::compute;
    use datafusion::arrow::datatypes::DataType;
    let casted = compute::cast(arr, &DataType::Utf8)?;
    casted
        .as_any()
        .downcast_ref::<StringArray>()
        .cloned()
        .ok_or_else(|| DataFusionError::Execution("Cannot cast to StringArray".into()))
}

fn downcast_u32_or_i32(arr: &dyn Array) -> Result<Vec<u32>> {
    if let Some(u) = arr.as_any().downcast_ref::<UInt32Array>() {
        return Ok((0..u.len()).map(|i| u.value(i)).collect());
    }
    if let Some(u) = arr.as_any().downcast_ref::<Int32Array>() {
        return Ok((0..u.len()).map(|i| u.value(i) as u32).collect());
    }
    use datafusion::arrow::array::UInt64Array;
    if let Some(u) = arr.as_any().downcast_ref::<UInt64Array>() {
        return Ok((0..u.len()).map(|i| u.value(i) as u32).collect());
    }
    use datafusion::arrow::array::Int64Array;
    if let Some(u) = arr.as_any().downcast_ref::<Int64Array>() {
        return Ok((0..u.len()).map(|i| u.value(i) as u32).collect());
    }
    Err(DataFusionError::Execution(
        "Cannot downcast 'start' to numeric array".into(),
    ))
}

fn downcast_f64(arr: &dyn Array) -> Result<Float64Array> {
    if let Some(f) = arr.as_any().downcast_ref::<Float64Array>() {
        return Ok(f.clone());
    }
    use datafusion::arrow::compute;
    use datafusion::arrow::datatypes::DataType;
    let casted = compute::cast(arr, &DataType::Float64)?;
    Ok(casted
        .as_any()
        .downcast_ref::<Float64Array>()
        .cloned()
        .unwrap())
}
