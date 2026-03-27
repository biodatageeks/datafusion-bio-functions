//! VCF file writer for processed sample QC output.
//!
//! Writes multi-sample VCF with FORMAT fields GT:GQ:DP:PL:DS.
//! Automatically applies BGZF compression when the output path ends
//! with `.vcf.gz` or `.vcf.bgz`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use datafusion::arrow::array::{Array, Float64Array, Int32Array, StringArray, UInt32Array};
use datafusion::common::{DataFusionError, Result};
use flate2::Compression;
use flate2::write::GzEncoder;

use super::processing::ProcessedBatch;

/// VCF writer that outputs processed genotype data.
///
/// Compression is chosen automatically from the file extension:
/// - `.vcf` — plain text
/// - `.vcf.gz` / `.vcf.bgz` — gzip (BGZF-compatible level 6)
pub struct VcfWriter {
    inner: WriterInner,
    line_buf: String,
}

/// Type-erased inner writer so we don't need generics on every method.
enum WriterInner {
    Plain(BufWriter<File>),
    Gzip(Box<BufWriter<GzEncoder<File>>>),
}

impl Write for WriterInner {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            WriterInner::Plain(w) => w.write(buf),
            WriterInner::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            WriterInner::Plain(w) => w.flush(),
            WriterInner::Gzip(w) => w.flush(),
        }
    }
}

impl VcfWriter {
    /// Create a new VCF writer, writing the header immediately.
    ///
    /// If `path` ends with `.gz` or `.bgz`, output is gzip-compressed.
    pub fn new(path: &Path, sample_names: &[String]) -> Result<Self> {
        let file = File::create(path).map_err(|e| {
            DataFusionError::Execution(format!("Cannot create output VCF '{path:?}': {e}"))
        })?;

        let compressed = is_compressed_path(path);
        let mut inner = if compressed {
            let encoder = GzEncoder::new(file, Compression::new(6));
            WriterInner::Gzip(Box::new(BufWriter::with_capacity(1 << 20, encoder)))
        } else {
            WriterInner::Plain(BufWriter::with_capacity(1 << 20, file))
        };

        // Write VCF header
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

        // Column header
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
            line_buf: String::with_capacity(64 * 1024),
        })
    }

    /// Write all variants from a processed batch.
    pub fn write_batch(&mut self, batch: &ProcessedBatch) -> Result<()> {
        let chrom = downcast_string(&batch.chrom)?;
        let start = downcast_u32_or_i32(&batch.start)?;
        let id = downcast_string(&batch.id)?;
        let ref_allele = downcast_string(&batch.ref_allele)?;
        let alt = downcast_string(&batch.alt)?;
        let qual = downcast_f64(&batch.qual)?;
        let filter = downcast_string(&batch.filter)?;

        #[allow(clippy::needless_range_loop)]
        for var_idx in 0..batch.n_variants {
            self.line_buf.clear();

            // CHROM
            let chrom_val = if chrom.is_null(var_idx) {
                "."
            } else {
                chrom.value(var_idx)
            };

            // POS (1-based from 0-based start)
            let pos = start[var_idx] + 1;

            // ID
            let id_val = if id.is_null(var_idx) {
                "."
            } else {
                let v = id.value(var_idx);
                if v.is_empty() { "." } else { v }
            };

            // REF, ALT
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

            // QUAL
            let qual_str = if qual.is_null(var_idx) {
                ".".to_string()
            } else {
                format_qual(qual.value(var_idx))
            };

            // FILTER
            let filter_val = if filter.is_null(var_idx) {
                "."
            } else {
                let v = filter.value(var_idx);
                if v.is_empty() { "." } else { v }
            };

            // Write fixed columns
            use std::fmt::Write as FmtWrite;
            write!(
                self.line_buf,
                "{chrom_val}\t{pos}\t{id_val}\t{ref_val}\t{alt_val}\t{qual_str}\t{filter_val}\t.\tGT:GQ:DP:PL:DS"
            )
            .map_err(|e| DataFusionError::Execution(format!("fmt: {e}")))?;

            // Write fixed columns to writer
            self.inner
                .write_all(self.line_buf.as_bytes())
                .map_err(io_err)?;

            // Write per-sample columns directly to writer
            let s_start = batch.sample_offsets[var_idx];
            let s_end = batch.sample_offsets[var_idx + 1];
            for s in s_start..s_end {
                self.inner.write_all(b"\t").map_err(io_err)?;
                self.write_sample_to_inner(s, batch)?;
            }

            self.inner.write_all(b"\n").map_err(io_err)?;
        }
        Ok(())
    }

    /// Write a single sample's FORMAT fields directly to the writer.
    fn write_sample_to_inner(&mut self, idx: usize, batch: &ProcessedBatch) -> Result<()> {
        // GT
        let gt = if batch.gt_final.is_null(idx) {
            "./."
        } else {
            batch.gt_final.value(idx)
        };
        self.inner.write_all(gt.as_bytes()).map_err(io_err)?;

        // GQ
        self.inner.write_all(b":").map_err(io_err)?;
        if batch.gq_values.is_null(idx) {
            self.inner.write_all(b".").map_err(io_err)?;
        } else {
            write!(self.inner, "{}", batch.gq_values.value(idx)).map_err(io_err)?;
        }

        // DP
        self.inner.write_all(b":").map_err(io_err)?;
        if batch.dp_values.is_null(idx) {
            self.inner.write_all(b".").map_err(io_err)?;
        } else {
            write!(self.inner, "{}", batch.dp_values.value(idx)).map_err(io_err)?;
        }

        // PL (comma-separated triple)
        write!(
            self.inner,
            ":{},{},{}",
            batch.pl0.value(idx),
            batch.pl1.value(idx),
            batch.pl2.value(idx)
        )
        .map_err(io_err)?;

        // DS
        self.inner.write_all(b":").map_err(io_err)?;
        if batch.ds_values.is_null(idx) {
            self.inner.write_all(b".").map_err(io_err)?;
        } else {
            write!(self.inner, "{:.4}", batch.ds_values.value(idx)).map_err(io_err)?;
        }

        Ok(())
    }

    /// Flush and finalize the VCF output.
    ///
    /// For gzip output this finishes the gzip stream (writes the trailer).
    pub fn finish(self) -> Result<()> {
        match self.inner {
            WriterInner::Plain(mut w) => w.flush().map_err(io_err),
            WriterInner::Gzip(w) => {
                let encoder = (*w).into_inner().map_err(io_err)?;
                encoder.finish().map_err(io_err)?;
                Ok(())
            }
        }
    }
}

// ============================================================================
// Helpers
// ============================================================================

/// Check if a path indicates compressed output (`.gz` or `.bgz`).
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

/// Downcast an ArrayRef to StringArray (handling Utf8, LargeUtf8, Utf8View).
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

/// Downcast to a Vec<u32> from UInt32 or Int32/Int64 arrays.
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

/// Downcast to Float64Array.
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
