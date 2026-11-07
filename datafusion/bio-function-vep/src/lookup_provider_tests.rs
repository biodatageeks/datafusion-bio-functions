use super::*;
use datafusion::execution::TaskContext;
use datafusion::physical_plan::empty::EmptyExec;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, PlanProperties, SendableRecordBatchStream,
};

#[derive(Debug)]
struct FailOnRead(EmptyExec);

impl DisplayAs for FailOnRead {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "FailOnRead")
    }
}

impl ExecutionPlan for FailOnRead {
    fn name(&self) -> &str {
        "FailOnRead"
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn properties(&self) -> &Arc<PlanProperties> {
        self.0.properties()
    }
    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![]
    }
    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        assert!(children.is_empty());
        Ok(self)
    }
    fn execute(
        &self,
        _partition: usize,
        _context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        Err(DataFusionError::Execution(
            "input must not be read during lookup planning".into(),
        ))
    }
}

fn cache_schema() -> Schema {
    Schema::new(
        crate::schema_contract::REQUIRED_VARIATION_COLUMNS
            .iter()
            .map(|(name, kind)| Field::new(*name, kind.clone(), false))
            .collect::<Vec<_>>(),
    )
}

fn lookup(session: Arc<SessionContext>, schema: SchemaRef) -> LookupProvider {
    let mut provider = LookupProvider::new(
        session,
        "input".into(),
        "unused_cache".into(),
        schema.as_ref().clone(),
        cache_schema(),
        vec!["variation_name".into()],
        true,
        0,
    )
    .unwrap();
    provider.set_cache_root("/nonexistent/planning-test-cache");
    provider
}

#[derive(Debug)]
struct UnreadableInput(SchemaRef);

#[async_trait]
impl TableProvider for UnreadableInput {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn schema(&self) -> SchemaRef {
        Arc::clone(&self.0)
    }
    fn table_type(&self) -> TableType {
        TableType::Temporary
    }
    async fn scan(
        &self,
        _state: &dyn Session,
        projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let plan = Arc::new(FailOnRead(EmptyExec::new(self.schema())));
        wrap_with_projection(plan, projection)
    }
}

#[tokio::test]
async fn lookup_planning_does_not_read_the_input() {
    let session = Arc::new(SessionContext::new());
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int32, false),
        Field::new("end", DataType::Int32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
    ]));
    session
        .register_table("input", Arc::new(UnreadableInput(Arc::clone(&schema))))
        .unwrap();
    lookup(Arc::clone(&session), schema)
        .scan(&session.state(), None, &[], None)
        .await
        .expect("lookup planning must not execute a prefix query");
}

/// Multiple BGZF blocks, tied starts, and enough multi-ALT rows to cross logical
/// 5000-unit buffers. Large coordinate span makes eight physical partitions real.
fn write_input(directory: &std::path::Path, chrom: &str, indexed: bool) -> PathBuf {
    use noodles_core_tabix::Position;
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    use std::io::Write;

    let path = directory.join(if indexed { "input.vcf.gz" } else { "input.vcf" });
    let header = format!(
        "##fileformat=VCFv4.2\n##contig=<ID={chrom},length=150000000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    );
    let mut rows = Vec::new();
    for i in 0..12_040 {
        let pos = (i / 2 + 1) * 16_384;
        let alt = if i % 101 == 0 { "T,G" } else { "T" };
        rows.push((pos, format!("{chrom}\t{pos}\tv{i}\tA\t{alt}\t.\tPASS\t.\n")));
    }
    let file = std::fs::File::create(&path).unwrap();
    if indexed {
        let mut writer = noodles_bgzf::io::Writer::new(file);
        writer.write_all(header.as_bytes()).unwrap();
        let mut indexer = noodles_tabix::index::Indexer::default();
        indexer.set_header(noodles_csi::binning_index::index::header::Builder::vcf().build());
        for (i, (pos, row)) in rows.iter().enumerate() {
            let start = writer.virtual_position();
            writer.write_all(row.as_bytes()).unwrap();
            if i % 100 == 99 {
                writer.flush().unwrap();
            }
            let position = Position::try_from(*pos).unwrap();
            indexer
                .add_record(
                    chrom,
                    position,
                    position,
                    Chunk::new(start, writer.virtual_position()),
                )
                .unwrap();
        }
        writer.finish().unwrap();
        let mut index = noodles_tabix::io::Writer::new(
            std::fs::File::create(path.with_extension("gz.tbi")).unwrap(),
        );
        index.write_index(&indexer.build()).unwrap();
    } else {
        let mut writer = std::io::BufWriter::new(file);
        writer.write_all(header.as_bytes()).unwrap();
        for (_, row) in rows {
            writer.write_all(row.as_bytes()).unwrap();
        }
    }
    path
}

#[tokio::test]
async fn scoped_grid_inputs_preserve_order_and_parent_parallelism() {
    use datafusion::physical_plan::ExecutionPlanProperties;
    use datafusion::prelude::{SessionConfig, col, lit};
    use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
    use futures::TryStreamExt;

    // Both chromosome spellings and the unindexed fallback retain their rows.
    for (chrom, indexed) in [("chr22", true), ("22", true), ("chr22", false)] {
        let directory = tempfile::tempdir().unwrap();
        let path = write_input(directory.path(), chrom, indexed);
        let session = Arc::new(SessionContext::new_with_config(
            SessionConfig::new()
                .with_target_partitions(8)
                .with_batch_size(640),
        ));
        let input = VcfTableProvider::new(
            path.to_str().unwrap().to_string(),
            Some(vec![]),
            Some(vec![]),
            None,
            false,
        )
        .unwrap();
        let schema = input.schema();
        session.register_table("input", Arc::new(input)).unwrap();
        let direct = session
            .table("input")
            .await
            .unwrap()
            .filter(col("chrom").eq(lit(chrom)))
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();
        let ordinary_count = direct.output_partitioning().partition_count();
        if indexed {
            assert_eq!(ordinary_count, 8, "fixture must exercise physical fan-out");
        }

        let plans = futures::future::try_join_all((0..8).map(|worker| {
            let mut provider = lookup(Arc::clone(&session), Arc::clone(&schema));
            // Include one tie group before each seam, as real warm-up windows do.
            let lo = (worker * 750_usize).saturating_sub(1) * 16_384 + 1;
            let hi = if worker == 7 {
                150_000_000
            } else {
                (worker + 1) * 750 * 16_384 + 1
            };
            provider.set_vcf_filter(Some(
                col("chrom")
                    .eq(lit(chrom))
                    .and(col("start").gt_eq(lit(lo as i64)))
                    .and(col("start").lt(lit(hi as i64))),
            ));
            provider.set_vcf_input_partitions(1);
            provider.set_target_partitions(8); // Independent variation-reader budget.
            provider.set_skip_probing(true);
            let state = session.state();
            async move {
                let plan = provider.scan(&state, None, &[], None).await?;
                Ok::<_, DataFusionError>((plan, lo, hi))
            }
        }))
        .await
        .unwrap();

        assert_eq!(session.state().config().target_partitions(), 8);
        let direct_after = session
            .table("input")
            .await
            .unwrap()
            .filter(col("chrom").eq(lit(chrom)))
            .unwrap()
            .create_physical_plan()
            .await
            .unwrap();
        assert_eq!(
            direct_after.output_partitioning().partition_count(),
            ordinary_count
        );
        // A regular LookupProvider keeps the session target too.
        let mut ordinary = lookup(Arc::clone(&session), Arc::clone(&schema));
        ordinary.set_vcf_filter(Some(col("chrom").eq(lit(chrom))));
        let ordinary = ordinary
            .scan(&session.state(), None, &[], None)
            .await
            .unwrap();
        assert_eq!(
            ordinary.output_partitioning().partition_count(),
            ordinary_count
        );

        for (plan, lo, hi) in plans {
            assert_eq!(plan.output_partitioning().partition_count(), 1);
            assert_eq!(
                plan.children()[0].output_partitioning().partition_count(),
                1
            );
            let batches: Vec<_> = plan
                .execute(0, session.task_ctx())
                .unwrap()
                .try_collect()
                .await
                .unwrap();
            let mut actual = Vec::new();
            for batch in batches {
                for row in 0..batch.num_rows() {
                    let values = ["chrom", "start", "id", "alt"].map(|name| {
                        datafusion::arrow::util::display::array_value_to_string(
                            batch.column_by_name(name).unwrap(),
                            row,
                        )
                        .unwrap()
                    });
                    actual.push(values);
                }
            }
            let expected: Vec<_> = (0..12_040)
                .filter_map(|i| {
                    let pos = (i / 2 + 1) * 16_384;
                    (lo..hi).contains(&pos).then(|| {
                        [
                            chrom.to_string(),
                            pos.to_string(),
                            format!("v{i}"),
                            // The reader's Arrow representation separates ALT
                            // alleles with '|'; their order remains significant.
                            if i % 101 == 0 { "T|G" } else { "T" }.into(),
                        ]
                    })
                })
                .collect();
            assert_eq!(actual.len(), expected.len(), "window must retain all rows");
            for (row, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
                assert_eq!(actual, expected, "input order changed at row {row}");
            }
        }
    }
}
