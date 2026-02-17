use std::any::Any;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::datatypes::{Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::Result;
use datafusion::datasource::{TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{Expr, SessionContext};

use crate::filter_op::FilterOp;

pub struct OverlapProvider {
    session: Arc<SessionContext>,
    left_table: String,
    right_table: String,
    left_schema: Schema,
    right_schema: Schema,
    columns_1: (String, String, String),
    columns_2: (String, String, String),
    filter_op: FilterOp,
    schema: SchemaRef,
}

impl OverlapProvider {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        session: Arc<SessionContext>,
        left_table: String,
        right_table: String,
        left_schema: Schema,
        right_schema: Schema,
        columns_1: Vec<String>,
        columns_2: Vec<String>,
        filter_op: FilterOp,
    ) -> Self {
        let mut fields = left_schema
            .fields()
            .iter()
            .map(|f| {
                Arc::new(Field::new(
                    format!("left_{}", f.name()),
                    f.data_type().clone(),
                    f.is_nullable(),
                ))
            })
            .collect::<Vec<_>>();
        fields.extend(right_schema.fields().iter().map(|f| {
            Arc::new(Field::new(
                format!("right_{}", f.name()),
                f.data_type().clone(),
                f.is_nullable(),
            ))
        }));
        let schema = Arc::new(Schema::new(fields));

        Self {
            session,
            left_table,
            right_table,
            left_schema,
            right_schema,
            columns_1: (
                columns_1[0].clone(),
                columns_1[1].clone(),
                columns_1[2].clone(),
            ),
            columns_2: (
                columns_2[0].clone(),
                columns_2[1].clone(),
                columns_2[2].clone(),
            ),
            filter_op,
            schema,
        }
    }
}

impl Debug for OverlapProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "OverlapProvider {{ left: {}, right: {} }}",
            self.left_table, self.right_table
        )
    }
}

#[async_trait]
impl TableProvider for OverlapProvider {
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
        _state: &dyn Session,
        _projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let select_left = self
            .left_schema
            .fields()
            .iter()
            .map(|f| format!("a.`{}` AS `left_{}`", f.name(), f.name()))
            .collect::<Vec<_>>()
            .join(", ");
        let select_right = self
            .right_schema
            .fields()
            .iter()
            .map(|f| format!("b.`{}` AS `right_{}`", f.name(), f.name()))
            .collect::<Vec<_>>()
            .join(", ");

        let sign = if self.filter_op == FilterOp::Strict {
            ""
        } else {
            "="
        };

        let (c1, s1, e1) = (&self.columns_1.0, &self.columns_1.1, &self.columns_1.2);
        let (c2, s2, e2) = (&self.columns_2.0, &self.columns_2.1, &self.columns_2.2);

        let query = format!(
            "SELECT {select_left}, {select_right} \
             FROM `{}` AS b, `{}` AS a \
             WHERE a.`{c1}` = b.`{c2}` \
             AND CAST(a.`{e1}` AS INTEGER) >{sign} CAST(b.`{s2}` AS INTEGER) \
             AND CAST(a.`{s1}` AS INTEGER) <{sign} CAST(b.`{e2}` AS INTEGER)",
            self.right_table, self.left_table,
        );

        let df = self.session.sql(&query).await?;
        df.create_physical_plan().await
    }
}
