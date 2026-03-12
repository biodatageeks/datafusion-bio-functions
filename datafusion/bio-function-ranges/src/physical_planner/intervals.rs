use datafusion::common::JoinSide;
use datafusion::common::tree_node::{Transformed, TransformedResult, TreeNode, TreeNodeRecursion};
use datafusion::logical_expr::Operator;
use datafusion::physical_expr::PhysicalExpr;
use datafusion::physical_expr::expressions::{BinaryExpr, Column, lit};
use datafusion::physical_plan::joins::utils::{ColumnIndex, JoinFilter};
use std::collections::{BTreeSet, HashMap};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ColInterval {
    start: Arc<dyn PhysicalExpr>,
    end: Arc<dyn PhysicalExpr>,
}

impl ColInterval {
    pub fn new(start: Arc<dyn PhysicalExpr>, end: Arc<dyn PhysicalExpr>) -> Self {
        Self { start, end }
    }
    pub fn start(&self) -> Arc<dyn PhysicalExpr> {
        self.start.clone()
    }
    pub fn end(&self) -> Arc<dyn PhysicalExpr> {
        self.end.clone()
    }
}

#[derive(Debug, Clone)]
pub struct ColIntervals {
    pub left_interval: ColInterval,
    pub right_interval: ColInterval,
}

#[derive(Debug, Clone)]
pub struct ParsedIntervalJoin {
    pub intervals: ColIntervals,
    pub residual_filter: Option<JoinFilter>,
}

pub fn parse(filter: Option<&JoinFilter>) -> Option<ParsedIntervalJoin> {
    if let Some(filter) = filter {
        try_parse(filter).inspect_err(|e| log::debug!("{e}")).ok()
    } else {
        log::debug!("filter is not provided");
        None
    }
}

fn map_column_to_source_schema(
    expr: Arc<dyn PhysicalExpr>,
    indices: &[ColumnIndex],
) -> (Arc<dyn PhysicalExpr>, JoinSide) {
    let mut side: Option<JoinSide> = None;
    let message = format!("complex sub queries are not supported {expr:?}");

    expr.transform_up(|node| {
        if let Some(column) = node.as_any().downcast_ref::<Column>() {
            let new_column = Column::new(column.name(), indices[column.index()].index);
            if side.is_some() {
                panic!("{}", message);
            }
            side = Some(indices[column.index()].side);
            Ok(Transformed::yes(Arc::new(new_column)))
        } else {
            Ok(Transformed::no(node))
        }
    })
    .data()
    .map(|expr| (expr, side.expect("side not found")))
    .unwrap()
}

fn minus_one(expr: Arc<dyn PhysicalExpr>) -> Arc<dyn PhysicalExpr> {
    Arc::new(BinaryExpr::new(expr, Operator::Minus, lit(1)))
}

fn parse_condition(
    expr: &BinaryExpr,
    indices: &[ColumnIndex],
    inner: &mut IntervalBuilder,
) -> Result<(), String> {
    let is_lt = matches!(expr.op(), Operator::Lt);
    let is_lteq = matches!(expr.op(), Operator::LtEq);
    let is_gteq = matches!(expr.op(), Operator::GtEq);
    let is_gt = matches!(expr.op(), Operator::Gt);

    if !(is_lteq || is_gteq || is_lt || is_gt) {
        return Err(format!("Unsupported operator: {}", expr.op()));
    }

    match map_column_to_source_schema(expr.left().clone(), indices) {
        (rs, JoinSide::Right) if is_lteq || is_lt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (le, JoinSide::Left) => {
                    let le = if is_lt { minus_one(le) } else { le };
                    inner.with_rs(rs)?.with_le(le)?;
                    Ok(())
                }
                _ => Err("couldn't parse as rs </<= le".to_string()),
            }
        }
        (ls, JoinSide::Left) if is_lteq || is_lt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (re, JoinSide::Right) => {
                    let re = if is_lt { minus_one(re) } else { re };
                    inner.with_re(re)?.with_ls(ls)?;
                    Ok(())
                }
                _ => Err("couldn't parse as ls </<= re".to_string()),
            }
        }
        (re, JoinSide::Right) if is_gteq || is_gt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (ls, JoinSide::Left) => {
                    let re = if is_gt { minus_one(re) } else { re };
                    inner.with_re(re)?.with_ls(ls)?;
                    Ok(())
                }
                _ => Err("couldn't parse as re >/>= ls".to_string()),
            }
        }
        (le, JoinSide::Left) if is_gteq || is_gt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (rs, JoinSide::Right) => {
                    let le = if is_gt { minus_one(le) } else { le };
                    inner.with_rs(rs)?.with_le(le)?;
                    Ok(())
                }
                _ => Err("couldn't parse as le >/>= rs".to_string()),
            }
        }
        _ => Err("couldn't parse left side as neither 'rs </<=' nor 'ls </<=' nor 're >/>=' nor 'le >/>='".to_string()),
    }
}

struct IntervalBuilder {
    ls: Option<Arc<dyn PhysicalExpr>>,
    le: Option<Arc<dyn PhysicalExpr>>,
    rs: Option<Arc<dyn PhysicalExpr>>,
    re: Option<Arc<dyn PhysicalExpr>>,
}

impl IntervalBuilder {
    fn empty() -> IntervalBuilder {
        IntervalBuilder {
            ls: None,
            le: None,
            rs: None,
            re: None,
        }
    }

    fn with_ls(&mut self, ls: Arc<dyn PhysicalExpr>) -> Result<&mut Self, String> {
        if self.ls.is_some() {
            return Err("ls must not be called twice".to_string());
        }
        self.ls = Some(ls);
        Ok(self)
    }
    fn with_le(&mut self, le: Arc<dyn PhysicalExpr>) -> Result<&mut Self, String> {
        if self.le.is_some() {
            return Err("le must not be called twice".to_string());
        }
        self.le = Some(le);
        Ok(self)
    }
    fn with_rs(&mut self, rs: Arc<dyn PhysicalExpr>) -> Result<&mut Self, String> {
        if self.rs.is_some() {
            return Err("rs must not be called twice".to_string());
        }
        self.rs = Some(rs);
        Ok(self)
    }
    fn with_re(&mut self, re: Arc<dyn PhysicalExpr>) -> Result<&mut Self, String> {
        if self.re.is_some() {
            return Err("re must not be called twice".to_string());
        }
        self.re = Some(re);
        Ok(self)
    }

    fn finish(self) -> ColIntervals {
        let l_interval = ColInterval {
            start: self.ls.expect("ls must be set"),
            end: self.le.expect("le must be set"),
        };
        let r_interval = ColInterval {
            start: self.rs.expect("rs must be set"),
            end: self.re.expect("re must be set"),
        };

        ColIntervals {
            left_interval: l_interval,
            right_interval: r_interval,
        }
    }

    fn is_complete(&self) -> bool {
        self.ls.is_some() && self.le.is_some() && self.rs.is_some() && self.re.is_some()
    }
}

fn try_parse(filter: &JoinFilter) -> Result<ParsedIntervalJoin, String> {
    let mut remaining = Vec::new();
    let mut inner = IntervalBuilder::empty();
    let indices = filter.column_indices();

    for expr in split_and_conjuncts(filter.expression().clone()) {
        let Some(binary) = expr.as_any().downcast_ref::<BinaryExpr>() else {
            remaining.push(expr);
            continue;
        };

        if parse_condition(binary, indices, &mut inner).is_err() {
            remaining.push(expr);
        }
    }

    if !inner.is_complete() {
        return Err("couldn't extract interval bounds from join filter".to_string());
    }

    let residual_filter = rebuild_residual_filter(filter, remaining)?;

    Ok(ParsedIntervalJoin {
        intervals: inner.finish(),
        residual_filter,
    })
}

fn split_and_conjuncts(expr: Arc<dyn PhysicalExpr>) -> Vec<Arc<dyn PhysicalExpr>> {
    if let Some(binary) = expr.as_any().downcast_ref::<BinaryExpr>() {
        if matches!(binary.op(), Operator::And) {
            let mut out = split_and_conjuncts(binary.left().clone());
            out.extend(split_and_conjuncts(binary.right().clone()));
            return out;
        }
    }

    vec![expr]
}

fn rebuild_residual_filter(
    filter: &JoinFilter,
    conjuncts: Vec<Arc<dyn PhysicalExpr>>,
) -> Result<Option<JoinFilter>, String> {
    let mut conjuncts = conjuncts.into_iter();
    let Some(mut expr) = conjuncts.next() else {
        return Ok(None);
    };

    for rhs in conjuncts {
        expr = Arc::new(BinaryExpr::new(expr, Operator::And, rhs));
    }

    let mut referenced = BTreeSet::new();
    expr.apply(|node| {
        if let Some(column) = node.as_any().downcast_ref::<Column>() {
            referenced.insert(column.index());
        }
        Ok(TreeNodeRecursion::Continue)
    })
    .map_err(|e| e.to_string())?;

    let remap = referenced
        .iter()
        .enumerate()
        .map(|(new_idx, old_idx)| (*old_idx, new_idx))
        .collect::<HashMap<_, _>>();

    let remapped_expr = expr
        .transform_up(|node| {
            if let Some(column) = node.as_any().downcast_ref::<Column>() {
                let new_idx = *remap
                    .get(&column.index())
                    .expect("referenced columns must have a remap entry");
                Ok(Transformed::yes(
                    Arc::new(Column::new(column.name(), new_idx)) as Arc<dyn PhysicalExpr>,
                ))
            } else {
                Ok(Transformed::no(node))
            }
        })
        .data()
        .map_err(|e| e.to_string())?;

    let schema = Arc::new(datafusion::arrow::datatypes::Schema::new(
        referenced
            .iter()
            .map(|idx| filter.schema().field(*idx).clone())
            .collect::<Vec<_>>(),
    ));
    let column_indices = referenced
        .iter()
        .map(|idx| filter.column_indices()[*idx].clone())
        .collect();

    Ok(Some(JoinFilter::new(remapped_expr, column_indices, schema)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::common::tree_node::TreeNodeRecursion;
    use datafusion::error::{DataFusionError, Result};
    use datafusion::physical_plan::ExecutionPlan;
    use datafusion::physical_plan::joins::{HashJoinExec, NestedLoopJoinExec};
    use datafusion::prelude::SessionContext;

    async fn extract_filter(condition: &str) -> Result<ParsedIntervalJoin> {
        let ctx = SessionContext::new();
        ctx.sql("CREATE TABLE IF NOT EXISTS a (contig TEXT, l_start INT, l_end INT) AS VALUES ('a', 1, 2), ('b', 3, 4)").await?;
        ctx.sql("CREATE TABLE IF NOT EXISTS b (contig TEXT, name TEXT, r_end INT, r_start INT) AS VALUES ('a','x', 1, 2), ('b','x', 3, 4)").await?;

        let query = format!("SELECT * FROM a JOIN b ON a.contig = b.contig AND {condition}");

        let ds = ctx.sql(query.as_str()).await?;
        let plan = ds.create_physical_plan().await?;

        let filter: &JoinFilter = find_join_filter(&plan)?;

        try_parse(filter).map_err(DataFusionError::Internal)
    }

    #[tokio::test]
    async fn test_all_comp_combinations_for_gteq_lteq() -> Result<()> {
        let l_start = Column::new("l_start", 1);
        let l_end = Column::new("l_end", 2);
        let r_start = Column::new("r_start", 3);
        let r_end = Column::new("r_end", 2);

        let condition = "b.r_end >= a.l_start AND a.l_end >= b.r_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_end >= a.l_start AND b.r_start <= a.l_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_start <= b.r_end AND a.l_end >= b.r_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_start <= b.r_end AND b.r_start <= a.l_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_end >= b.r_start AND b.r_end >= a.l_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_end >= b.r_start AND a.l_start <= b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_start <= a.l_end AND b.r_end >= a.l_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_start <= a.l_end AND a.l_start <= b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_start <= a.l_end OR a.l_start <= b.r_end";
        let error = extract_filter(condition).await;
        assert!(error.is_err());

        Ok(())
    }

    #[tokio::test]
    async fn test_all_comp_combinations_for_gt_lt() -> Result<()> {
        let l_start = Column::new("l_start", 1);
        let l_end = BinaryExpr::new(Arc::new(Column::new("l_end", 2)), Operator::Minus, lit(1));
        let r_start = Column::new("r_start", 3);
        let r_end = BinaryExpr::new(Arc::new(Column::new("r_end", 2)), Operator::Minus, lit(1));

        let condition = "b.r_end > a.l_start AND a.l_end > b.r_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "b.r_end > a.l_start AND b.r_start < a.l_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "a.l_start < b.r_end AND a.l_end > b.r_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "a.l_start < b.r_end AND b.r_start < a.l_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "a.l_end > b.r_start AND b.r_end > a.l_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "a.l_end > b.r_start AND a.l_start < b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "b.r_start < a.l_end AND b.r_end > a.l_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "b.r_start < a.l_end AND a.l_start < b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        let condition = "b.r_start < a.l_end AND a.l_start <= b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), Column::new("r_end", 2));

        let condition = "b.r_start <= a.l_end AND a.l_start < b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?.intervals;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), Column::new("l_end", 2));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(
            format!("{:?}", to_binary(right.end)),
            format!("{:?}", r_end)
        );

        Ok(())
    }

    #[tokio::test]
    async fn test_extracts_residual_filter() -> Result<()> {
        let parsed = extract_filter(
            "b.r_end >= a.l_start AND a.l_end >= b.r_start AND a.l_start != b.r_start",
        )
        .await?;

        let residual = parsed
            .residual_filter
            .expect("expected residual filter to be preserved");
        assert_eq!(residual.column_indices().len(), 2);
        assert_eq!(format!("{residual}"), "l_start != r_start");

        Ok(())
    }

    #[tokio::test]
    #[should_panic]
    async fn test_all_comp_combinations_for_complex_query() {
        let condition = "(b.r_end - a.l_start) >= a.l_start AND a.l_end >= b.r_start";
        extract_filter(condition).await.unwrap();
    }

    fn find_join_filter(plan: &Arc<dyn ExecutionPlan>) -> Result<&JoinFilter> {
        let mut filter: Option<&JoinFilter> = None;
        plan.apply(|plan| {
            Ok(
                if let Some(hash) = plan.as_any().downcast_ref::<HashJoinExec>() {
                    filter = hash.filter();
                    TreeNodeRecursion::Stop
                } else if let Some(nested) = plan.as_any().downcast_ref::<NestedLoopJoinExec>() {
                    filter = nested.filter();
                    TreeNodeRecursion::Stop
                } else {
                    TreeNodeRecursion::Continue
                },
            )
        })
        .map(|_| filter.expect("filter not found"))
    }

    fn to_column(expr: Arc<dyn PhysicalExpr>) -> Column {
        expr.as_any().downcast_ref::<Column>().unwrap().clone()
    }
    fn to_binary(expr: Arc<dyn PhysicalExpr>) -> BinaryExpr {
        expr.as_any().downcast_ref::<BinaryExpr>().unwrap().clone()
    }
}
