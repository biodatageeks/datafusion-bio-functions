use std::any::Any;
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::catalog::Session;
use datafusion::catalog::TableFunctionImpl;
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::{TableProvider, TableType};
use datafusion::logical_expr::Expr;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::SessionContext;

use datafusion_bio_format_bam::table_provider::BamTableProvider;

use crate::physical_exec::{PileupConfig, PileupExec};
use crate::schema::coverage_output_schema;

/// A TableProvider that wraps a child provider (e.g., BAM reader) and produces coverage output.
pub struct CoverageTableProvider {
    input: Arc<dyn TableProvider>,
    config: PileupConfig,
}

impl CoverageTableProvider {
    /// Create a new CoverageTableProvider wrapping the given input provider.
    pub fn new(input: Arc<dyn TableProvider>, config: PileupConfig) -> Self {
        Self { input, config }
    }
}

impl std::fmt::Debug for CoverageTableProvider {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CoverageTableProvider")
            .field("config", &self.config)
            .finish()
    }
}

#[async_trait]
impl TableProvider for CoverageTableProvider {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        coverage_output_schema()
    }

    fn table_type(&self) -> TableType {
        TableType::Temporary
    }

    async fn scan(
        &self,
        state: &dyn Session,
        _projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        // Only request columns needed for coverage: chrom(1), start(2), flags(4), cigar(5), mapq(6)
        let pileup_projection = vec![1, 2, 4, 5, 6];
        let child_plan = self
            .input
            .scan(state, Some(&pileup_projection), &[], None)
            .await?;
        Ok(Arc::new(PileupExec::new(child_plan, self.config.clone())))
    }
}

/// A table function that enables SQL: `SELECT * FROM coverage('path/to/file.bam')`
#[derive(Debug, Default)]
pub struct CoverageFunction;

impl TableFunctionImpl for CoverageFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.is_empty() {
            return Err(DataFusionError::Plan(
                "coverage() requires at least one argument: the BAM file path".to_string(),
            ));
        }

        // Extract file path from first argument
        let file_path = match &args[0] {
            Expr::Literal(ScalarValue::Utf8(Some(path)), _) => path.clone(),
            other => {
                return Err(DataFusionError::Plan(format!(
                    "coverage() first argument must be a string literal, got: {other}"
                )));
            }
        };

        // Create BamTableProvider — handle both tokio and non-tokio contexts
        // Always request binary CIGAR for zero-copy performance
        let bam_provider = match tokio::runtime::Handle::try_current() {
            Ok(handle) => tokio::task::block_in_place(|| {
                handle.block_on(BamTableProvider::new(file_path, None, true, None, true))
            })?,
            Err(_) => {
                let rt = tokio::runtime::Runtime::new()
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                rt.block_on(BamTableProvider::new(file_path, None, true, None, true))?
            }
        };

        let config = PileupConfig {
            binary_cigar: true,
            ..PileupConfig::default()
        };
        Ok(Arc::new(CoverageTableProvider::new(
            Arc::new(bam_provider),
            config,
        )))
    }
}

/// Register pileup/coverage table functions on a SessionContext.
pub fn register_pileup_functions(ctx: &SessionContext) {
    ctx.register_udtf("coverage", Arc::new(CoverageFunction));
}
