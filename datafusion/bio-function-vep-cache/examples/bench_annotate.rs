use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{Array, RecordBatch};
use datafusion::common::Result;
use datafusion::physical_plan::ExecutionPlanProperties;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use futures::StreamExt;

use datafusion_bio_function_vep_cache::allele_index::{AlleleMatcher, WindowAlleleIndex};
use datafusion_bio_function_vep_cache::key_encoding::window_id_for_position;
use datafusion_bio_function_vep_cache::kv_store::{FORMAT_V1, VepKvStore};

/// Simple exact matcher.
fn allele_matches(_vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
    let mut parts = allele_string.split('/');
    let Some(_cache_ref) = parts.next() else {
        return false;
    };
    parts.any(|a| a == vcf_alt)
}

fn get_string_col(col: &dyn Array) -> Vec<Option<String>> {
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringArray>()
    {
        (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect()
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StringViewArray>()
    {
        (0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_string())
                }
            })
            .collect()
    } else {
        (0..col.len()).map(|_| None).collect()
    }
}

fn get_u32_col(col: &dyn Array) -> Vec<u32> {
    col.as_any()
        .downcast_ref::<datafusion::arrow::array::UInt32Array>()
        .map(|a| (0..a.len()).map(|i| a.value(i)).collect())
        .unwrap_or_default()
}

struct AnnotationState {
    store: Arc<VepKvStore>,
    window_size: u64,
    format_version: u8,
    matcher: AlleleMatcher,
    current_window: Option<(String, u64)>,
    current_index: Option<WindowAlleleIndex>,
    window_loads: u64,
}

impl AnnotationState {
    fn new(store: Arc<VepKvStore>) -> Self {
        let window_size = store.window_size();
        let format_version = store.format_version();
        Self {
            store,
            window_size,
            format_version,
            matcher: allele_matches,
            current_window: None,
            current_index: None,
            window_loads: 0,
        }
    }

    fn ensure_window(&mut self, chrom: &str, pos: i64) -> Result<()> {
        let wid = window_id_for_position(pos, self.window_size);
        let need_load = match &self.current_window {
            Some((c, w)) => c != chrom || *w != wid,
            None => true,
        };

        if need_load {
            if self.format_version >= FORMAT_V1 {
                let pos_index = self.store.get_position_index(chrom, wid)?;
                self.current_window = Some((chrom.to_string(), wid));
                self.current_index = pos_index.map(WindowAlleleIndex::from_position_index);
            } else {
                let batch = self.store.get_window(chrom, wid)?;
                self.current_window = Some((chrom.to_string(), wid));
                self.current_index = match batch {
                    Some(b) => Some(WindowAlleleIndex::from_batch(b)?),
                    None => None,
                };
            }
            self.window_loads += 1;
        }
        Ok(())
    }

    fn annotate_batch(&mut self, batch: &RecordBatch) -> Result<(usize, usize)> {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom")?;
        let start_idx = schema.index_of("start")?;
        let ref_idx = schema.index_of("ref")?;
        let alt_idx = schema.index_of("alt")?;

        let chroms = get_string_col(batch.column(chrom_idx));
        let starts = get_u32_col(batch.column(start_idx));
        let refs = get_string_col(batch.column(ref_idx));
        let alts = get_string_col(batch.column(alt_idx));

        let mut matched = 0usize;
        let mut total = 0usize;

        for row in 0..batch.num_rows() {
            total += 1;
            let raw_chrom = chroms[row].as_deref().unwrap_or("");
            let chrom = raw_chrom.strip_prefix("chr").unwrap_or(raw_chrom);

            // VCF reader in 0-based half-open -> cache is 1-based closed
            let norm_start = starts[row] as i64 + 1;
            let norm_end = norm_start;

            self.ensure_window(chrom, norm_start)?;

            if let Some(index) = &self.current_index {
                let vcf_ref = refs[row].as_deref().unwrap_or("");
                let vcf_alt = alts[row].as_deref().unwrap_or("");
                let matches =
                    index.find_matches(norm_start, norm_end, vcf_ref, vcf_alt, self.matcher);
                if !matches.is_empty() {
                    matched += 1;
                }
            }
        }

        Ok((total, matched))
    }
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <fjall_cache_path> <vcf_path> [threads: 1,2,4,8] [cache_size_mb: 32]",
            args[0]
        );
        std::process::exit(1);
    }

    let cache_path = &args[1];
    let vcf_path = &args[2];
    let thread_list: Vec<usize> = if let Some(t) = args.get(3) {
        t.split(',').filter_map(|s| s.trim().parse().ok()).collect()
    } else {
        vec![1, 2, 4, 8]
    };
    // Cache size in MB (default 32)
    let cache_size_mb: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(32);
    let cache_size_bytes = cache_size_mb * 1024 * 1024;

    let store = Arc::new(VepKvStore::open_with_cache_size(
        cache_path,
        cache_size_bytes,
    )?);
    println!(
        "Cache: {} columns, window_size={}, cache_size={}MB, path={}",
        store.schema().fields().len(),
        store.window_size(),
        cache_size_mb,
        cache_path,
    );
    println!("VCF: {}", vcf_path);
    println!();
    println!(
        "{:>7} | {:>5} | {:>8} | {:>8} | {:>7} | {:>6} | {:>10} | {:>12}",
        "threads", "parts", "rows", "matched", "rate", "wins", "time (s)", "variants/s"
    );
    println!("{}", "-".repeat(85));

    for &num_threads in &thread_list {
        let config = datafusion::prelude::SessionConfig::new().with_target_partitions(num_threads);
        let ctx = SessionContext::new_with_config(config);

        let table =
            VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
        ctx.register_table("vcf", Arc::new(table))?;

        let df = ctx
            .sql("SELECT `chrom`, `start`, `ref`, `alt` FROM vcf")
            .await?;

        let physical_plan = df.create_physical_plan().await?;
        let num_partitions = physical_plan.output_partitioning().partition_count();
        let task_ctx = ctx.task_ctx();

        let start_time = Instant::now();

        // Spawn one task per partition, each with its own AnnotationState.
        let mut handles = Vec::new();
        for partition in 0..num_partitions {
            let plan = physical_plan.clone();
            let ctx = task_ctx.clone();
            let store = store.clone();

            let handle = tokio::spawn(async move {
                let mut stream = plan.execute(partition, ctx)?;
                let mut state = AnnotationState::new(store);
                let mut total_rows = 0usize;
                let mut total_matched = 0usize;

                while let Some(batch_result) = stream.next().await {
                    let batch = batch_result?;
                    let (rows, matched) = state.annotate_batch(&batch)?;
                    total_rows += rows;
                    total_matched += matched;
                }

                Ok::<_, datafusion::common::DataFusionError>((
                    total_rows,
                    total_matched,
                    state.window_loads,
                ))
            });
            handles.push(handle);
        }

        // Collect results from all partitions.
        let mut total_rows = 0usize;
        let mut total_matched = 0usize;
        let mut total_windows = 0u64;
        for handle in handles {
            let (rows, matched, wins) = handle.await.unwrap()?;
            total_rows += rows;
            total_matched += matched;
            total_windows += wins;
        }

        let elapsed = start_time.elapsed().as_secs_f64();
        let rate = if total_rows > 0 {
            total_matched as f64 / total_rows as f64 * 100.0
        } else {
            0.0
        };

        println!(
            "{:>7} | {:>5} | {:>8} | {:>8} | {:>6.1}% | {:>6} | {:>10.2} | {:>12.0}",
            num_threads,
            num_partitions,
            total_rows,
            total_matched,
            rate,
            total_windows,
            elapsed,
            total_rows as f64 / elapsed,
        );
    }

    Ok(())
}
