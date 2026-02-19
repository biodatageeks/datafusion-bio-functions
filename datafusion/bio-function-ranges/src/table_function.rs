use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::TableFunctionImpl;
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;

use crate::cluster::ClusterProvider;
use crate::complement::ComplementProvider;
use crate::count_overlaps::CountOverlapsProvider;
use crate::filter_op::FilterOp;
use crate::merge::MergeProvider;
use crate::nearest::NearestProvider;
use crate::overlap::OverlapProvider;
use crate::subtract::SubtractProvider;

/// Resolve a table's Arrow schema synchronously, handling both inside and
/// outside a tokio runtime.
fn resolve_table_schema(session: &SessionContext, table_name: &str) -> Result<Schema> {
    match tokio::runtime::Handle::try_current() {
        Ok(handle) => tokio::task::block_in_place(|| {
            let df = handle.block_on(session.table(table_name))?;
            Ok::<_, DataFusionError>(df.schema().as_arrow().clone())
        }),
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            let df = rt.block_on(session.table(table_name))?;
            Ok(df.schema().as_arrow().clone())
        }
    }
}

const DEFAULT_COLS: [&str; 3] = ["contig", "pos_start", "pos_end"];

type ColTriple = (String, String, String);

/// Extract a string literal from an Expr, returning an error with context on failure.
///
/// Rejects values containing backticks to prevent SQL injection in generated queries.
fn extract_string_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<String> {
    match arg {
        Expr::Literal(ScalarValue::Utf8(Some(val)), _) => {
            if val.contains('`') {
                return Err(DataFusionError::Plan(format!(
                    "{fn_name}() {name} must not contain backtick characters, got: {val}"
                )));
            }
            Ok(val.clone())
        }
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
    parse_col_args_extra(&args[2..], fn_name)
}

#[allow(clippy::type_complexity)]
fn parse_col_args_extra(extra: &[Expr], fn_name: &str) -> Result<(ColTriple, ColTriple, FilterOp)> {
    if extra.is_empty() {
        let cols = (
            DEFAULT_COLS[0].to_string(),
            DEFAULT_COLS[1].to_string(),
            DEFAULT_COLS[2].to_string(),
        );
        return Ok((cols.clone(), cols, FilterOp::Weak));
    }

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

fn extract_usize_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<usize> {
    let as_u64 = match arg {
        Expr::Literal(ScalarValue::Int64(Some(v)), _) => u64::try_from(*v).ok(),
        Expr::Literal(ScalarValue::Int32(Some(v)), _) => u64::try_from(*v).ok(),
        Expr::Literal(ScalarValue::UInt64(Some(v)), _) => Some(*v),
        Expr::Literal(ScalarValue::UInt32(Some(v)), _) => Some(u64::from(*v)),
        _ => None,
    };

    if let Some(v) = as_u64 {
        if v == 0 {
            return Err(DataFusionError::Plan(format!(
                "{fn_name}() {name} must be >= 1"
            )));
        }
        return usize::try_from(v).map_err(|_| {
            DataFusionError::Plan(format!("{fn_name}() {name}={v} does not fit usize"))
        });
    }

    Err(DataFusionError::Plan(format!(
        "{fn_name}() {name} must be an integer literal"
    )))
}

fn extract_bool_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<bool> {
    match arg {
        Expr::Literal(ScalarValue::Boolean(Some(v)), _) => Ok(*v),
        _ => Err(DataFusionError::Plan(format!(
            "{fn_name}() {name} must be a boolean literal"
        ))),
    }
}

/// Extract an optional leading `min_dist` integer argument from an argument slice.
///
/// Returns `(min_dist, remaining_args)`. If the first argument is not a numeric literal,
/// returns `(0, original_args)`.
fn extract_min_dist<'a>(extra: &'a [Expr], fn_name: &str) -> Result<(i64, &'a [Expr])> {
    if extra.is_empty() {
        return Ok((0i64, extra));
    }
    match &extra[0] {
        Expr::Literal(ScalarValue::Int64(Some(v)), _) => {
            if *v < 0 {
                return Err(DataFusionError::Plan(format!(
                    "{fn_name}() min_dist must be >= 0, got {v}"
                )));
            }
            Ok((*v, &extra[1..]))
        }
        Expr::Literal(ScalarValue::Int32(Some(v)), _) => {
            if *v < 0 {
                return Err(DataFusionError::Plan(format!(
                    "{fn_name}() min_dist must be >= 0, got {v}"
                )));
            }
            Ok((i64::from(*v), &extra[1..]))
        }
        Expr::Literal(ScalarValue::UInt64(Some(v)), _) => Ok((
            i64::try_from(*v).map_err(|_| {
                DataFusionError::Plan(format!("{fn_name}() min_dist value {v} does not fit i64"))
            })?,
            &extra[1..],
        )),
        Expr::Literal(ScalarValue::UInt32(Some(v)), _) => Ok((i64::from(*v), &extra[1..])),
        _ => Ok((0i64, extra)),
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

/// Internal table function for nearest interval matching.
struct NearestTableFunction {
    session: Arc<SessionContext>,
    name: &'static str,
}

impl std::fmt::Debug for NearestTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for NearestTableFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.len() < 2 {
            return Err(DataFusionError::Plan(format!(
                "{}() requires at least 2 arguments: left_table and right_table names",
                self.name
            )));
        }

        let left_table = extract_string_arg(&args[0], "left_table", self.name)?;
        let right_table = extract_string_arg(&args[1], "right_table", self.name)?;

        let mut arg_idx = 2usize;
        let mut k = 1usize;
        let mut include_overlaps = true;

        if let Some(arg) = args.get(arg_idx) {
            if matches!(
                arg,
                Expr::Literal(ScalarValue::Int64(Some(_)), _)
                    | Expr::Literal(ScalarValue::Int32(Some(_)), _)
                    | Expr::Literal(ScalarValue::UInt64(Some(_)), _)
                    | Expr::Literal(ScalarValue::UInt32(Some(_)), _)
            ) {
                k = extract_usize_arg(arg, "k", self.name)?;
                arg_idx += 1;
            }
        }

        if let Some(arg) = args.get(arg_idx) {
            if matches!(arg, Expr::Literal(ScalarValue::Boolean(Some(_)), _)) {
                include_overlaps = extract_bool_arg(arg, "overlap", self.name)?;
                arg_idx += 1;
            }
        }

        let mut compute_distance = true;
        if let Some(arg) = args.get(arg_idx) {
            if matches!(arg, Expr::Literal(ScalarValue::Boolean(Some(_)), _)) {
                compute_distance = extract_bool_arg(arg, "compute_distance", self.name)?;
                arg_idx += 1;
            }
        }

        let (cols_left, cols_right, filter_op) = parse_col_args_extra(&args[arg_idx..], self.name)?;

        let (left_schema, right_schema) = match tokio::runtime::Handle::try_current() {
            Ok(handle) => tokio::task::block_in_place(|| {
                let left = handle.block_on(self.session.table(&left_table))?;
                let right = handle.block_on(self.session.table(&right_table))?;
                Ok::<_, DataFusionError>((
                    left.schema().as_arrow().clone(),
                    right.schema().as_arrow().clone(),
                ))
            }),
            Err(_) => {
                let rt = tokio::runtime::Runtime::new()
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                let left = rt.block_on(self.session.table(&left_table))?;
                let right = rt.block_on(self.session.table(&right_table))?;
                Ok((
                    left.schema().as_arrow().clone(),
                    right.schema().as_arrow().clone(),
                ))
            }
        }?;

        Ok(Arc::new(NearestProvider::new(
            Arc::clone(&self.session),
            left_table,
            right_table,
            left_schema,
            right_schema,
            vec![cols_left.0, cols_left.1, cols_left.2],
            vec![cols_right.0, cols_right.1, cols_right.2],
            filter_op,
            include_overlaps,
            k,
            compute_distance,
        )))
    }
}

/// Internal table function for overlap (all pairs of overlapping intervals).
struct OverlapTableFunction {
    session: Arc<SessionContext>,
    name: &'static str,
}

impl std::fmt::Debug for OverlapTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for OverlapTableFunction {
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

        let (left_schema, right_schema) = match tokio::runtime::Handle::try_current() {
            Ok(handle) => tokio::task::block_in_place(|| {
                let left = handle.block_on(self.session.table(&left_table))?;
                let right = handle.block_on(self.session.table(&right_table))?;
                Ok::<_, DataFusionError>((
                    left.schema().as_arrow().clone(),
                    right.schema().as_arrow().clone(),
                ))
            }),
            Err(_) => {
                let rt = tokio::runtime::Runtime::new()
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;
                let left = rt.block_on(self.session.table(&left_table))?;
                let right = rt.block_on(self.session.table(&right_table))?;
                Ok((
                    left.schema().as_arrow().clone(),
                    right.schema().as_arrow().clone(),
                ))
            }
        }?;

        Ok(Arc::new(OverlapProvider::new(
            Arc::clone(&self.session),
            left_table,
            right_table,
            left_schema,
            right_schema,
            vec![cols_left.0, cols_left.1, cols_left.2],
            vec![cols_right.0, cols_right.1, cols_right.2],
            filter_op,
        )))
    }
}

/// Internal table function for merge (merge overlapping intervals within a single table).
struct MergeTableFunction {
    session: Arc<SessionContext>,
    name: &'static str,
}

impl std::fmt::Debug for MergeTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for MergeTableFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.is_empty() {
            return Err(DataFusionError::Plan(format!(
                "{}() requires at least 1 argument: table name",
                self.name
            )));
        }

        let table = extract_string_arg(&args[0], "table", self.name)?;

        let extra = &args[1..];
        let (min_dist, col_extra) = extract_min_dist(extra, self.name)?;
        let (cols, filter_op) = parse_merge_col_args(col_extra, self.name)?;

        Ok(Arc::new(MergeProvider::new(
            Arc::clone(&self.session),
            table,
            cols,
            min_dist,
            filter_op,
        )))
    }
}

/// Parse column + filter_op arguments for merge (single table).
///
/// Supports:
/// - No args: default columns, Weak
/// - 3 args: column names
/// - +1 optional 'strict'/'weak' at the end
fn parse_merge_col_args(
    extra: &[Expr],
    fn_name: &str,
) -> Result<((String, String, String), FilterOp)> {
    if extra.is_empty() {
        return Ok((
            (
                DEFAULT_COLS[0].to_string(),
                DEFAULT_COLS[1].to_string(),
                DEFAULT_COLS[2].to_string(),
            ),
            FilterOp::Weak,
        ));
    }

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
        0 => Ok((
            (
                DEFAULT_COLS[0].to_string(),
                DEFAULT_COLS[1].to_string(),
                DEFAULT_COLS[2].to_string(),
            ),
            filter_op,
        )),
        3 => {
            let c0 = extract_string_arg(&col_args[0], "column name", fn_name)?;
            let c1 = extract_string_arg(&col_args[1], "column name", fn_name)?;
            let c2 = extract_string_arg(&col_args[2], "column name", fn_name)?;
            Ok(((c0, c1, c2), filter_op))
        }
        n => Err(DataFusionError::Plan(format!(
            "{fn_name}() expects 0 or 3 column name arguments (got {n}). \
             Usage: {fn_name}('table' [, min_dist] [, col1, col2, col3] [, 'strict'|'weak'])"
        ))),
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

/// Internal table function for cluster (annotate intervals with cluster membership).
struct ClusterTableFunction {
    session: Arc<SessionContext>,
    name: &'static str,
}

impl std::fmt::Debug for ClusterTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for ClusterTableFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.is_empty() {
            return Err(DataFusionError::Plan(format!(
                "{}() requires at least 1 argument: table name",
                self.name
            )));
        }

        let table = extract_string_arg(&args[0], "table", self.name)?;

        let extra = &args[1..];
        let (min_dist, col_extra) = extract_min_dist(extra, self.name)?;
        let (cols, filter_op) = parse_merge_col_args(col_extra, self.name)?;

        let input_schema = resolve_table_schema(&self.session, &table)?;

        Ok(Arc::new(ClusterProvider::new(
            Arc::clone(&self.session),
            table,
            cols,
            min_dist,
            filter_op,
            input_schema,
        )))
    }
}

/// Internal table function for complement (find gaps / uncovered regions).
struct ComplementTableFunction {
    session: Arc<SessionContext>,
    name: &'static str,
}

impl std::fmt::Debug for ComplementTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for ComplementTableFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.is_empty() {
            return Err(DataFusionError::Plan(format!(
                "{}() requires at least 1 argument: table name",
                self.name
            )));
        }

        let table = extract_string_arg(&args[0], "table", self.name)?;
        let extra = &args[1..];

        // Distinguish between 1-table and 2-table patterns:
        // 1 arg:  complement('table')
        // 2 args: complement('table', 'view_table')
        // 4 args: complement('table', 'c', 's', 'e') or complement('table', 'c', 's', 'e')
        // 5 args: complement('table', 'c', 's', 'e', 'strict')
        //         OR complement('table', 'view_table', 'c', 's', 'e')
        // 8 args: complement('table', 'view_table', 'c1', 's1', 'e1', 'vc', 'vs', 've')
        // etc.
        //
        // Strategy: if extra[0] is a string, check if extra has enough args for
        // a 2-table pattern. We detect "view_table" by checking if the second
        // arg looks like a table name (string) and the total arg count >= 2.

        if extra.is_empty() {
            // complement('table') — no view, default columns
            let cols = (
                DEFAULT_COLS[0].to_string(),
                DEFAULT_COLS[1].to_string(),
                DEFAULT_COLS[2].to_string(),
            );
            return Ok(Arc::new(ComplementProvider::new(
                Arc::clone(&self.session),
                table,
                None,
                cols.clone(),
                cols,
                FilterOp::Weak,
            )));
        }

        // Try to extract the second argument as a potential view table name.
        // If the second arg is a string, it could be either a view_table name
        // or a column name. We disambiguate by count:
        //   - 1 extra + possible filter_op → view_table (no custom cols)
        //   - 3 extra + possible filter_op → custom cols, no view
        //   - 4+ extra where first is string → view_table + cols
        let second_is_string = matches!(&extra[0], Expr::Literal(ScalarValue::Utf8(Some(_)), _));

        if !second_is_string {
            return Err(DataFusionError::Plan(format!(
                "{}() unexpected argument type at position 2",
                self.name
            )));
        }

        // Check for trailing filter_op in extra
        let (work_args, filter_op) =
            if let Some(Expr::Literal(ScalarValue::Utf8(Some(val)), _)) = extra.last() {
                match val.to_lowercase().as_str() {
                    "strict" => (&extra[..extra.len() - 1], FilterOp::Strict),
                    "weak" => (&extra[..extra.len() - 1], FilterOp::Weak),
                    _ => (extra, FilterOp::Weak),
                }
            } else {
                (extra, FilterOp::Weak)
            };

        match work_args.len() {
            0 => {
                // Just filter_op was the only thing
                let cols = (
                    DEFAULT_COLS[0].to_string(),
                    DEFAULT_COLS[1].to_string(),
                    DEFAULT_COLS[2].to_string(),
                );
                Ok(Arc::new(ComplementProvider::new(
                    Arc::clone(&self.session),
                    table,
                    None,
                    cols.clone(),
                    cols,
                    filter_op,
                )))
            }
            1 => {
                // complement('table', 'view_table')
                let view_table = extract_string_arg(&work_args[0], "view_table", self.name)?;
                let cols = (
                    DEFAULT_COLS[0].to_string(),
                    DEFAULT_COLS[1].to_string(),
                    DEFAULT_COLS[2].to_string(),
                );
                Ok(Arc::new(ComplementProvider::new(
                    Arc::clone(&self.session),
                    table,
                    Some(view_table),
                    cols.clone(),
                    cols,
                    filter_op,
                )))
            }
            3 => {
                // complement('table', 'c', 's', 'e')
                let c = extract_string_arg(&work_args[0], "column name", self.name)?;
                let s = extract_string_arg(&work_args[1], "column name", self.name)?;
                let e = extract_string_arg(&work_args[2], "column name", self.name)?;
                let cols = (c, s, e);
                Ok(Arc::new(ComplementProvider::new(
                    Arc::clone(&self.session),
                    table,
                    None,
                    cols.clone(),
                    cols,
                    filter_op,
                )))
            }
            4 => {
                // complement('table', 'view_table', 'c', 's', 'e')
                let view_table = extract_string_arg(&work_args[0], "view_table", self.name)?;
                let c = extract_string_arg(&work_args[1], "column name", self.name)?;
                let s = extract_string_arg(&work_args[2], "column name", self.name)?;
                let e = extract_string_arg(&work_args[3], "column name", self.name)?;
                let cols = (c, s, e);
                Ok(Arc::new(ComplementProvider::new(
                    Arc::clone(&self.session),
                    table,
                    Some(view_table),
                    cols.clone(),
                    cols,
                    filter_op,
                )))
            }
            7 => {
                // complement('table', 'view_table', 'c1', 's1', 'e1', 'vc', 'vs', 've')
                let view_table = extract_string_arg(&work_args[0], "view_table", self.name)?;
                let c1 = extract_string_arg(&work_args[1], "column name", self.name)?;
                let s1 = extract_string_arg(&work_args[2], "column name", self.name)?;
                let e1 = extract_string_arg(&work_args[3], "column name", self.name)?;
                let vc = extract_string_arg(&work_args[4], "view column name", self.name)?;
                let vs = extract_string_arg(&work_args[5], "view column name", self.name)?;
                let ve = extract_string_arg(&work_args[6], "view column name", self.name)?;
                Ok(Arc::new(ComplementProvider::new(
                    Arc::clone(&self.session),
                    table,
                    Some(view_table),
                    (c1, s1, e1),
                    (vc, vs, ve),
                    filter_op,
                )))
            }
            n => Err(DataFusionError::Plan(format!(
                "{}() unexpected number of arguments ({n} extra args after table name). \
                 Usage: {name}('table' [, 'view_table'] [, col1, col2, col3 [, vcol1, vcol2, vcol3]] [, 'strict'|'weak'])",
                self.name,
                name = self.name
            ))),
        }
    }
}

/// Internal table function for subtract (basepair-level set difference).
struct SubtractTableFunction {
    session: Arc<SessionContext>,
    name: &'static str,
}

impl std::fmt::Debug for SubtractTableFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}Function", self.name)
    }
}

impl TableFunctionImpl for SubtractTableFunction {
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

        let left_schema = resolve_table_schema(&self.session, &left_table)?;

        Ok(Arc::new(SubtractProvider::new(
            Arc::clone(&self.session),
            left_table,
            right_table,
            cols_left,
            cols_right,
            filter_op,
            left_schema,
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
            session: Arc::clone(&session),
            coverage: false,
            name: "count_overlaps",
        }),
    );
    ctx.register_udtf(
        "nearest",
        Arc::new(NearestTableFunction {
            session: Arc::clone(&session),
            name: "nearest",
        }),
    );
    ctx.register_udtf(
        "overlap",
        Arc::new(OverlapTableFunction {
            session: Arc::clone(&session),
            name: "overlap",
        }),
    );
    ctx.register_udtf(
        "merge",
        Arc::new(MergeTableFunction {
            session: Arc::clone(&session),
            name: "merge",
        }),
    );
    ctx.register_udtf(
        "cluster",
        Arc::new(ClusterTableFunction {
            session: Arc::clone(&session),
            name: "cluster",
        }),
    );
    ctx.register_udtf(
        "complement",
        Arc::new(ComplementTableFunction {
            session: Arc::clone(&session),
            name: "complement",
        }),
    );
    ctx.register_udtf(
        "subtract",
        Arc::new(SubtractTableFunction {
            session,
            name: "subtract",
        }),
    );
}
