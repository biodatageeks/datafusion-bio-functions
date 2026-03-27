//! Table-Valued Function implementation for `vcf_sample_qc`.
//!
//! Reads a registered VCF table, applies site filtering and per-sample genotype
//! QC transformations, writes output to VCF, and returns a summary row.

use std::any::Any;
use std::fmt::{self, Debug, Formatter};
use std::path::PathBuf;
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use async_trait::async_trait;
use datafusion::arrow::array::{Int64Array, RecordBatch, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::TableFunctionImpl;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::TableProvider;
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::logical_expr::{Expr, TableType};
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{DisplayAs, DisplayFormatType, ExecutionPlan, PlanProperties};
use datafusion::prelude::SessionContext;
use futures::Stream;

use super::processing;
use super::writer::VcfWriter;
use super::{QcConfig, extract_f64_arg, extract_string_arg};

// ============================================================================
// TableFunctionImpl
// ============================================================================

/// The `vcf_sample_qc` table function.
pub struct VcfSampleQcFunction {
    session: Arc<SessionContext>,
}

impl VcfSampleQcFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
        Self { session }
    }
}

impl Debug for VcfSampleQcFunction {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("VcfSampleQcFunction").finish()
    }
}

impl TableFunctionImpl for VcfSampleQcFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.len() < 3 {
            return Err(DataFusionError::Plan(
                "vcf_sample_qc(table_name, output_path, sample_names [, qual_min, \
                 site_gq_avg_min, site_dp_avg_min, site_dp_avg_max, \
                 sample_gq_min, sample_dp_min, sample_dp_max])"
                    .into(),
            ));
        }

        let table_name = extract_string_arg(&args[0], "table_name")?;
        let output_path = extract_string_arg(&args[1], "output_path")?;
        let sample_names_raw = extract_string_arg(&args[2], "sample_names")?;

        let sample_names: Vec<String> = if sample_names_raw == "auto" {
            Vec::new() // will be determined from data
        } else if is_sample_file_path(&sample_names_raw) {
            read_sample_names_file(&sample_names_raw)?
        } else {
            sample_names_raw
                .split(',')
                .map(|s| s.trim().to_string())
                .collect()
        };

        let mut config = QcConfig::default();
        if args.len() > 3 {
            config.qual_min = extract_f64_arg(&args[3], "qual_min")?;
        }
        if args.len() > 4 {
            config.site_gq_avg_min = extract_f64_arg(&args[4], "site_gq_avg_min")?;
        }
        if args.len() > 5 {
            config.site_dp_avg_min = extract_f64_arg(&args[5], "site_dp_avg_min")?;
        }
        if args.len() > 6 {
            config.site_dp_avg_max = extract_f64_arg(&args[6], "site_dp_avg_max")?;
        }
        if args.len() > 7 {
            config.sample_gq_min = extract_f64_arg(&args[7], "sample_gq_min")?;
        }
        if args.len() > 8 {
            config.sample_dp_min = extract_f64_arg(&args[8], "sample_dp_min")?;
        }
        if args.len() > 9 {
            config.sample_dp_max = extract_f64_arg(&args[9], "sample_dp_max")?;
        }

        Ok(Arc::new(VcfSampleQcProvider {
            session: Arc::clone(&self.session),
            table_name,
            output_path: PathBuf::from(output_path),
            sample_names,
            config,
            schema: summary_schema(),
        }))
    }
}

// ============================================================================
// TableProvider
// ============================================================================

struct VcfSampleQcProvider {
    session: Arc<SessionContext>,
    table_name: String,
    output_path: PathBuf,
    sample_names: Vec<String>,
    config: QcConfig,
    schema: SchemaRef,
}

impl Debug for VcfSampleQcProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("VcfSampleQcProvider")
            .field("table_name", &self.table_name)
            .field("output_path", &self.output_path)
            .finish()
    }
}

#[async_trait]
impl TableProvider for VcfSampleQcProvider {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn table_type(&self) -> TableType {
        TableType::Temporary
    }

    async fn scan(
        &self,
        _state: &dyn datafusion::catalog::Session,
        _projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        // Build the input execution plan from the registered table
        let df = self.session.table(&self.table_name).await?;
        let input_plan = df.create_physical_plan().await?;

        Ok(Arc::new(VcfSampleQcExec {
            input: input_plan,
            output_path: self.output_path.clone(),
            sample_names: self.sample_names.clone(),
            config: self.config.clone(),
            schema: self.schema.clone(),
            cache: PlanProperties::new(
                EquivalenceProperties::new(self.schema.clone()),
                datafusion::physical_plan::Partitioning::UnknownPartitioning(1),
                EmissionType::Final,
                Boundedness::Bounded,
            ),
        }))
    }
}

// ============================================================================
// ExecutionPlan
// ============================================================================

struct VcfSampleQcExec {
    input: Arc<dyn ExecutionPlan>,
    output_path: PathBuf,
    sample_names: Vec<String>,
    config: QcConfig,
    schema: SchemaRef,
    cache: PlanProperties,
}

impl Debug for VcfSampleQcExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("VcfSampleQcExec")
            .field("output_path", &self.output_path)
            .finish()
    }
}

impl DisplayAs for VcfSampleQcExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter) -> fmt::Result {
        write!(f, "VcfSampleQcExec: output={}", self.output_path.display())
    }
}

impl ExecutionPlan for VcfSampleQcExec {
    fn name(&self) -> &str {
        "VcfSampleQcExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn properties(&self) -> &PlanProperties {
        &self.cache
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.input]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        if children.len() != 1 {
            return Err(DataFusionError::Plan(
                "VcfSampleQcExec requires exactly 1 child".into(),
            ));
        }
        Ok(Arc::new(VcfSampleQcExec {
            input: children[0].clone(),
            output_path: self.output_path.clone(),
            sample_names: self.sample_names.clone(),
            config: self.config.clone(),
            schema: self.schema.clone(),
            cache: self.cache.clone(),
        }))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let input_stream = self.input.execute(partition, context)?;
        Ok(Box::pin(VcfSampleQcStream {
            input: input_stream,
            output_path: self.output_path.clone(),
            sample_names: self.sample_names.clone(),
            config: self.config.clone(),
            schema: self.schema.clone(),
            state: StreamState::Reading,
            writer: None,
            variants_in: 0,
            variants_out: 0,
            n_samples: 0,
        }))
    }
}

// ============================================================================
// Stream
// ============================================================================

enum StreamState {
    Reading,
    Done,
}

struct VcfSampleQcStream {
    input: SendableRecordBatchStream,
    output_path: PathBuf,
    sample_names: Vec<String>,
    config: QcConfig,
    schema: SchemaRef,
    state: StreamState,
    writer: Option<VcfWriter>,
    variants_in: i64,
    variants_out: i64,
    n_samples: i64,
}

impl VcfSampleQcStream {
    /// Lazily initialize the writer once we know the number of samples.
    fn ensure_writer(&mut self, n_samples: usize) -> Result<()> {
        if self.writer.is_some() {
            return Ok(());
        }

        // Determine sample names
        let names = if self.sample_names.is_empty() || self.sample_names == ["auto"] {
            (0..n_samples)
                .map(|i| format!("SAMPLE_{i}"))
                .collect::<Vec<_>>()
        } else if self.sample_names.len() != n_samples {
            return Err(DataFusionError::Execution(format!(
                "Provided {} sample names but data has {} samples per variant",
                self.sample_names.len(),
                n_samples
            )));
        } else {
            self.sample_names.clone()
        };

        self.n_samples = n_samples as i64;
        self.writer = Some(VcfWriter::new(&self.output_path, &names)?);
        Ok(())
    }

    fn process_input_batch(&mut self, batch: &RecordBatch) -> Result<()> {
        let n_in = batch.num_rows();
        self.variants_in += n_in as i64;

        if let Some(processed) = processing::process_batch(batch, &self.config)? {
            // Initialize writer on first successful batch
            if processed.n_variants > 0 && !processed.sample_offsets.is_empty() {
                let n_samples = if processed.sample_offsets.len() > 1 {
                    processed.sample_offsets[1] - processed.sample_offsets[0]
                } else {
                    0
                };
                self.ensure_writer(n_samples)?;
            }

            self.variants_out += processed.n_variants as i64;

            if let Some(ref mut w) = self.writer {
                w.write_batch(&processed)?;
            }
        }
        Ok(())
    }

    fn build_summary(&self) -> Result<RecordBatch> {
        let schema = self.schema.clone();
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(Int64Array::from(vec![self.variants_in])),
                Arc::new(Int64Array::from(vec![self.variants_out])),
                Arc::new(Int64Array::from(vec![self.n_samples])),
                Arc::new(StringArray::from(vec![
                    self.output_path.to_string_lossy().to_string(),
                ])),
            ],
        )?;
        Ok(batch)
    }
}

impl Stream for VcfSampleQcStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        use futures::StreamExt;

        loop {
            match self.state {
                StreamState::Reading => {
                    let batch_opt = futures::ready!(self.input.poll_next_unpin(cx));
                    match batch_opt {
                        Some(Ok(batch)) => {
                            if batch.num_rows() == 0 {
                                continue;
                            }
                            if let Err(e) = self.process_input_batch(&batch) {
                                return Poll::Ready(Some(Err(e)));
                            }
                            continue;
                        }
                        Some(Err(e)) => return Poll::Ready(Some(Err(e))),
                        None => {
                            // Input exhausted — finalize writer and emit summary
                            self.state = StreamState::Done;
                            if let Some(w) = self.writer.take()
                                && let Err(e) = w.finish()
                            {
                                return Poll::Ready(Some(Err(e)));
                            }
                            let summary = self.build_summary();
                            return Poll::Ready(Some(summary));
                        }
                    }
                }
                StreamState::Done => {
                    return Poll::Ready(None);
                }
            }
        }
    }
}

impl RecordBatchStream for VcfSampleQcStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

// ============================================================================
// Helpers
// ============================================================================

/// Summary schema returned by the TVF.
fn summary_schema() -> SchemaRef {
    Arc::new(Schema::new(vec![
        Field::new("variants_in", DataType::Int64, false),
        Field::new("variants_out", DataType::Int64, false),
        Field::new("samples", DataType::Int64, false),
        Field::new("output_path", DataType::Utf8, false),
    ]))
}

/// Check whether the sample_names argument looks like a file path rather than
/// a comma-separated list. Heuristic: starts with `/` or `.`, or ends with `.txt`.
fn is_sample_file_path(s: &str) -> bool {
    s.starts_with('/') || s.starts_with("./") || s.starts_with("~/") || s.ends_with(".txt")
}

/// Read sample names from a text file (one name per line).
/// Blank lines and lines starting with `#` are skipped.
fn read_sample_names_file(path: &str) -> Result<Vec<String>> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        DataFusionError::Execution(format!("Cannot read sample names file '{path}': {e}"))
    })?;
    let names: Vec<String> = content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.to_string())
        .collect();
    if names.is_empty() {
        return Err(DataFusionError::Execution(format!(
            "Sample names file '{path}' contains no sample names"
        )));
    }
    Ok(names)
}
