use std::sync::Arc;

use datafusion::catalog::TableFunctionImpl;
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;

use crate::count_overlaps::CountOverlapsProvider;
use crate::filter_op::FilterOp;

const DEFAULT_COLS: [&str; 3] = ["contig", "pos_start", "pos_end"];

type ColTriple = (String, String, String);

/// Extract a string literal from an Expr, returning an error with context on failure.
fn extract_string_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<String> {
    match arg {
        Expr::Literal(ScalarValue::Utf8(Some(val)), _) => Ok(val.clone()),
        other => Err(DataFusionError::Plan(format!(
            "{fn_name}() {name} must be a string literal, got: {other}"
        ))),
    }
}

/// Parse column and filter_op arguments from the argument list.
///
/// Supports these argument patterns (after the two required table name args):
/// - No extra args: default columns (contig, pos_start, pos_end) for both, FilterOp::Weak
/// - 3 args: shared column names for both tables
/// - 6 args: separate column names for left and right tables
/// - +1 optional 'strict'/'weak' arg at the end of any pattern above
#[allow(clippy::type_complexity)]
fn parse_col_args(args: &[Expr], fn_name: &str) -> Result<(ColTriple, ColTriple, FilterOp)> {
    let extra = &args[2..];

    // Check for trailing filter_op argument
    let (col_args, filter_op) =
        if let Some(Expr::Literal(ScalarValue::Utf8(Some(val)), _)) = extra.last() {
            match val.to_lowercase().as_str() {
                "strict" => (&extra[..extra.len() - 1], FilterOp::Strict),
                "weak" => (&extra[..extra.len() - 1], FilterOp::Weak),
                _ => (extra, FilterOp::Weak),
            }
        } else {
            (extra, FilterOp::Weak)
        };

    match col_args.len() {
        0 => {
            let cols = (
                DEFAULT_COLS[0].to_string(),
                DEFAULT_COLS[1].to_string(),
                DEFAULT_COLS[2].to_string(),
            );
            Ok((cols.clone(), cols, filter_op))
        }
        3 => {
            let c0 = extract_string_arg(&col_args[0], "column name", fn_name)?;
            let c1 = extract_string_arg(&col_args[1], "column name", fn_name)?;
            let c2 = extract_string_arg(&col_args[2], "column name", fn_name)?;
            let cols = (c0, c1, c2);
            Ok((cols.clone(), cols, filter_op))
        }
        6 => {
            let l = (
                extract_string_arg(&col_args[0], "left column name", fn_name)?,
                extract_string_arg(&col_args[1], "left column name", fn_name)?,
                extract_string_arg(&col_args[2], "left column name", fn_name)?,
            );
            let r = (
                extract_string_arg(&col_args[3], "right column name", fn_name)?,
                extract_string_arg(&col_args[4], "right column name", fn_name)?,
                extract_string_arg(&col_args[5], "right column name", fn_name)?,
            );
            Ok((l, r, filter_op))
        }
        n => Err(DataFusionError::Plan(format!(
            "{fn_name}() expects 0, 3, or 6 column name arguments (got {n}). \
             Usage: {fn_name}('left_table', 'right_table' [, col1, col2, col3 \
             [, col4, col5, col6]] [, 'strict'|'weak'])"
        ))),
    }
}

/// Internal table function that implements both coverage and count_overlaps.
struct RangeTableFunction {
    session: Arc<SessionContext>,
    coverage: bool,
    name: &'static str,
}

impl std::fmt::Debug for RangeTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for RangeTableFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.len() < 2 {
            return Err(DataFusionError::Plan(format!(
                "{}() requires at least 2 arguments: left_table and right_table names",
                self.name
            )));
        }

        let left_table = extract_string_arg(&args[0], "left_table", self.name)?;
        let right_table = extract_string_arg(&args[1], "right_table", self.name)?;
        let (cols_left, cols_right, filter_op) = parse_col_args(args, self.name)?;

        // Look up the right table schema
        let right_schema = match tokio::runtime::Handle::try_current() {
            Ok(handle) => {
                tokio::task::block_in_place(|| handle.block_on(self.session.table(&right_table)))
            }
            Err(_) => {
                let rt = tokio::runtime::Runtime::new()
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                rt.block_on(self.session.table(&right_table))
            }
        }?
        .schema()
        .as_arrow()
        .clone();

        Ok(Arc::new(CountOverlapsProvider::new(
            Arc::clone(&self.session),
            left_table,
            right_table,
            right_schema,
            vec![cols_left.0, cols_left.1, cols_left.2],
            vec![cols_right.0, cols_right.1, cols_right.2],
            filter_op,
            self.coverage,
        )))
    }
}

/// Register coverage and count_overlaps table functions on a [`SessionContext`].
///
/// After registration, these SQL table functions are available:
///
/// ```sql
/// -- Count overlapping intervals from reads for each target
/// SELECT * FROM count_overlaps('reads', 'targets')
///
/// -- Coverage (base-pair overlap) with default column names
/// SELECT * FROM coverage('reads', 'targets')
///
/// -- With custom column names (shared for both tables)
/// SELECT * FROM coverage('reads', 'targets', 'chrom', 'start', 'end')
///
/// -- With separate column names for left and right tables
/// SELECT * FROM coverage('reads', 'targets', 'chrom', 'start', 'end', 'contig', 'pos_start', 'pos_end')
///
/// -- With filter_op for 0-based half-open coordinates
/// SELECT * FROM coverage('reads', 'targets', 'contig', 'pos_start', 'pos_end', 'strict')
/// ```
pub fn register_ranges_functions(ctx: &SessionContext) {
    let session = Arc::new(ctx.clone());
    ctx.register_udtf(
        "coverage",
        Arc::new(RangeTableFunction {
            session: Arc::clone(&session),
            coverage: true,
            name: "coverage",
        }),
    );
    ctx.register_udtf(
        "count_overlaps",
        Arc::new(RangeTableFunction {
            session,
            coverage: false,
            name: "count_overlaps",
        }),
    );
}
