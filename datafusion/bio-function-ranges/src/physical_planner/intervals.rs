use datafusion::common::JoinSide;
use datafusion::common::tree_node::{Transformed, TransformedResult, TreeNode};
use datafusion::logical_expr::Operator;
use datafusion::physical_expr::PhysicalExpr;
use datafusion::physical_expr::expressions::{BinaryExpr, Column, lit};
use datafusion::physical_plan::joins::utils::{ColumnIndex, JoinFilter};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ColInterval {
    start: Arc<dyn PhysicalExpr>,
    end: Arc<dyn PhysicalExpr>,
}

impl ColInterval {
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

pub fn parse(filter: Option<&JoinFilter>) -> Option<ColIntervals> {
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
                    inner.with_rs(rs).with_le(le);
                    Ok(())
                }
                _ => Err("couldn't parse as rs </<= le".to_string()),
            }
        }
        (ls, JoinSide::Left) if is_lteq || is_lt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (re, JoinSide::Right) => {
                    let re = if is_lt { minus_one(re) } else { re };
                    inner.with_re(re).with_ls(ls);
                    Ok(())
                }
                _ => Err("couldn't parse as ls </<= re".to_string()),
            }
        }
        (re, JoinSide::Right) if is_gteq || is_gt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (ls, JoinSide::Left) => {
                    let re = if is_gt { minus_one(re) } else { re };
                    inner.with_re(re).with_ls(ls);
                    Ok(())
                }
                _ => Err("couldn't parse as re >/>= ls".to_string()),
            }
        }
        (le, JoinSide::Left) if is_gteq || is_gt => {
            match map_column_to_source_schema(expr.right().clone(), indices) {
                (rs, JoinSide::Right) => {
                    let le = if is_gt { minus_one(le) } else { le };
                    inner.with_rs(rs).with_le(le);
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

    fn with_ls(&mut self, ls: Arc<dyn PhysicalExpr>) -> &mut Self {
        if self.ls.is_some() {
            panic!("ls must not be called twice")
        }
        self.ls = Some(ls);
        self
    }
    fn with_le(&mut self, le: Arc<dyn PhysicalExpr>) -> &mut Self {
        if self.le.is_some() {
            panic!("le must not be called twice")
        }
        self.le = Some(le);
        self
    }
    fn with_rs(&mut self, rs: Arc<dyn PhysicalExpr>) -> &mut Self {
        if self.rs.is_some() {
            panic!("rs must not be called twice")
        }
        self.rs = Some(rs);
        self
    }
    fn with_re(&mut self, re: Arc<dyn PhysicalExpr>) -> &mut Self {
        if self.re.is_some() {
            panic!("re must not be called twice")
        }
        self.re = Some(re);
        self
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
}

fn try_parse(filter: &JoinFilter) -> Result<ColIntervals, String> {
    let mut inner = IntervalBuilder::empty();

    let expr = filter
        .expression()
        .as_any()
        .downcast_ref::<BinaryExpr>()
        .ok_or("couldn't cast filter to BinaryExpr")?;

    if !matches!(expr.op(), Operator::And) {
        return Err("expr.op() is not AND".to_string());
    }

    let left = expr
        .left()
        .as_any()
        .downcast_ref::<BinaryExpr>()
        .ok_or("couldn't cast left side to BinaryExpr")?;
    let right = expr
        .right()
        .as_any()
        .downcast_ref::<BinaryExpr>()
        .ok_or("couldn't cast right side to BinaryExpr")?;
    let indices = filter.column_indices();

    parse_condition(left, indices, &mut inner)?;
    parse_condition(right, indices, &mut inner)?;

    Ok(inner.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::common::tree_node::TreeNodeRecursion;
    use datafusion::error::{DataFusionError, Result};
    use datafusion::physical_plan::ExecutionPlan;
    use datafusion::physical_plan::joins::{HashJoinExec, NestedLoopJoinExec};
    use datafusion::prelude::SessionContext;

    async fn extract_filter(condition: &str) -> Result<ColIntervals> {
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
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_end >= a.l_start AND b.r_start <= a.l_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_start <= b.r_end AND a.l_end >= b.r_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_start <= b.r_end AND b.r_start <= a.l_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_end >= b.r_start AND b.r_end >= a.l_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "a.l_end >= b.r_start AND a.l_start <= b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_start <= a.l_end AND b.r_end >= a.l_start";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(to_column(left.end), l_end);
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), r_end);

        let condition = "b.r_start <= a.l_end AND a.l_start <= b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

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
        } = extract_filter(condition).await?;

        assert_eq!(to_column(left.start), l_start);
        assert_eq!(format!("{:?}", to_binary(left.end)), format!("{:?}", l_end));
        assert_eq!(to_column(right.start), r_start);
        assert_eq!(to_column(right.end), Column::new("r_end", 2));

        let condition = "b.r_start <= a.l_end AND a.l_start < b.r_end";
        let ColIntervals {
            left_interval: left,
            right_interval: right,
        } = extract_filter(condition).await?;

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
