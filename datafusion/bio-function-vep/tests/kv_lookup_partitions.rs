use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use datafusion::arrow::array::{Int64Array, RecordBatch, StringArray};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::datasource::{MemTable, TableProvider};
use datafusion::physical_plan::{ExecutionPlan, ExecutionPlanProperties};
use datafusion::prelude::SessionContext;
use datafusion_bio_function_vep::allele::allele_matches;
use datafusion_bio_function_vep::kv_cache::cache_exec::{KvLookupExec, KvMatchMode};
use datafusion_bio_function_vep::kv_cache::position_entry::serialize_position_entry;
use datafusion_bio_function_vep::kv_cache::{KvCacheTableProvider, VepKvStore};
use datafusion_bio_function_vep::lookup_provider::LookupProvider;
use datafusion_bio_function_vep::variant_lookup_exec::{ColocatedSink, ColocatedSinkValue};
use futures::StreamExt;

fn vcf_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
    ]))
}

fn vcf_batch(pos: i64, alt: &str) -> RecordBatch {
    RecordBatch::try_new(
        vcf_schema(),
        vec![
            Arc::new(StringArray::from(vec!["1"])),
            Arc::new(Int64Array::from(vec![pos])),
            Arc::new(Int64Array::from(vec![pos])),
            Arc::new(StringArray::from(vec!["A"])),
            Arc::new(StringArray::from(vec![alt])),
        ],
    )
    .unwrap()
}

fn cache_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("variation_name", DataType::Utf8, true),
        Field::new("allele_string", DataType::Utf8, false),
        Field::new("failed", DataType::Int64, false),
    ]))
}

fn cache_batch(pos: i64, name: &str, allele: &str) -> RecordBatch {
    RecordBatch::try_new(
        cache_schema(),
        vec![
            Arc::new(StringArray::from(vec!["1"])),
            Arc::new(Int64Array::from(vec![pos])),
            Arc::new(Int64Array::from(vec![pos])),
            Arc::new(StringArray::from(vec![name])),
            Arc::new(StringArray::from(vec![allele])),
            Arc::new(Int64Array::from(vec![0])),
        ],
    )
    .unwrap()
}

fn create_store() -> (tempfile::TempDir, Arc<VepKvStore>) {
    let dir = tempfile::tempdir().unwrap();
    let store = VepKvStore::create(dir.path(), cache_schema()).unwrap();
    for (pos, name, allele) in [(100, "rs_p0", "A/G"), (200, "rs_p1", "A/T")] {
        let batch = cache_batch(pos, name, allele);
        let entry = serialize_position_entry(&[0], &batch, &[2, 3, 4, 5], 4).unwrap();
        store.put_position_entry("1", pos, &entry).unwrap();
    }
    store.persist().unwrap();
    (dir, Arc::new(store))
}

async fn collect_partition(exec: Arc<KvLookupExec>, partition: usize, ctx: &SessionContext) {
    let mut stream = exec.execute(partition, ctx.task_ctx()).unwrap();
    while let Some(batch) = stream.next().await {
        batch.unwrap();
    }
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn kv_lookup_uses_partition_specific_colocated_sinks() {
    let ctx = SessionContext::new();
    let vcf = MemTable::try_new(
        vcf_schema(),
        vec![vec![vcf_batch(100, "G")], vec![vcf_batch(200, "T")]],
    )
    .unwrap();
    ctx.register_table("vcf", Arc::new(vcf)).unwrap();
    let input_plan = ctx
        .table("vcf")
        .await
        .unwrap()
        .create_physical_plan()
        .await
        .unwrap();
    assert_eq!(input_plan.output_partitioning().partition_count(), 2);

    let (_dir, store) = create_store();
    let sinks: Vec<ColocatedSink> = (0..2)
        .map(|_| Arc::new(Mutex::new(HashMap::<_, ColocatedSinkValue>::new())))
        .collect();
    let exec = Arc::new(
        KvLookupExec::new(
            input_plan,
            store,
            vec!["variation_name".to_string()],
            KvMatchMode::Exact,
            allele_matches as fn(&str, &str, &str) -> bool,
            false,
            false,
            false,
            true,
            0,
        )
        .unwrap()
        .with_colocated_partition_sinks(sinks.clone()),
    );

    collect_partition(exec.clone(), 0, &ctx).await;
    assert_eq!(sinks[0].lock().unwrap().len(), 1);
    assert_eq!(sinks[1].lock().unwrap().len(), 0);

    collect_partition(exec, 1, &ctx).await;
    assert_eq!(sinks[0].lock().unwrap().len(), 1);
    assert_eq!(sinks[1].lock().unwrap().len(), 1);
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn lookup_provider_passes_partition_sinks_to_kv_exec() {
    let ctx = Arc::new(SessionContext::new());
    let vcf = MemTable::try_new(
        vcf_schema(),
        vec![vec![vcf_batch(100, "G")], vec![vcf_batch(200, "T")]],
    )
    .unwrap();
    ctx.register_table("vcf", Arc::new(vcf)).unwrap();

    let (_dir, store) = create_store();
    ctx.register_table("cache", Arc::new(KvCacheTableProvider::from_store(store)))
        .unwrap();

    let sinks: Vec<ColocatedSink> = (0..2)
        .map(|_| Arc::new(Mutex::new(HashMap::<_, ColocatedSinkValue>::new())))
        .collect();

    let mut provider = LookupProvider::new(
        ctx.clone(),
        "vcf".to_string(),
        "cache".to_string(),
        vcf_schema().as_ref().clone(),
        cache_schema().as_ref().clone(),
        vec!["variation_name".to_string()],
        true,
        0,
        None,
    )
    .unwrap();
    provider.set_partition_colocated_sinks(sinks.clone());

    let state = ctx.state();
    let plan = provider.scan(&state, None, &[], None).await.unwrap();
    assert_eq!(plan.output_partitioning().partition_count(), 2);

    for partition in 0..2 {
        let mut stream = plan.execute(partition, ctx.task_ctx()).unwrap();
        while let Some(batch) = stream.next().await {
            batch.unwrap();
        }
    }

    assert_eq!(sinks[0].lock().unwrap().len(), 1);
    assert_eq!(sinks[1].lock().unwrap().len(), 1);
}
