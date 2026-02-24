use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::Array;
use datafusion::common::Result;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use futures::StreamExt;

use datafusion_bio_function_vep_cache::allele_index::WindowAlleleIndex;
use datafusion_bio_function_vep_cache::key_encoding::window_id_for_position;
use datafusion_bio_function_vep_cache::kv_store::VepKvStore;

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

struct Stats {
    vcf_rows: u64,
    matched_rows: u64,
    window_loads: u64,
    total_window_rows_deserialized: u64,
    total_window_rows_probed: u64,       // position lookups attempted
    total_window_rows_position_hit: u64, // rows at matching positions
    total_window_rows_allele_hit: u64,   // rows passing allele match
}

impl Stats {
    fn new() -> Self {
        Self {
            vcf_rows: 0,
            matched_rows: 0,
            window_loads: 0,
            total_window_rows_deserialized: 0,
            total_window_rows_probed: 0,
            total_window_rows_position_hit: 0,
            total_window_rows_allele_hit: 0,
        }
    }
}

struct AnnotationState {
    store: Arc<VepKvStore>,
    window_size: u64,
    current_window: Option<(String, u64)>,
    current_index: Option<WindowAlleleIndex>,
    stats: Stats,
}

impl AnnotationState {
    fn new(store: Arc<VepKvStore>) -> Self {
        let window_size = store.window_size();
        Self {
            store,
            window_size,
            current_window: None,
            current_index: None,
            stats: Stats::new(),
        }
    }

    fn ensure_window(&mut self, chrom: &str, pos: i64) -> Result<()> {
        let wid = window_id_for_position(pos, self.window_size);
        let need_load = match &self.current_window {
            Some((c, w)) => c != chrom || *w != wid,
            None => true,
        };

        if need_load {
            let pos_index = self.store.get_position_index(chrom, wid)?;
            self.current_window = Some((chrom.to_string(), wid));
            self.stats.window_loads += 1;
            match pos_index {
                Some(idx) => {
                    self.stats.total_window_rows_deserialized += idx.num_rows() as u64;
                    self.current_index = Some(WindowAlleleIndex::from_position_index(idx));
                }
                None => {
                    self.current_index = None;
                }
            }
        }
        Ok(())
    }

    fn annotate_batch(&mut self, batch: &datafusion::arrow::array::RecordBatch) -> Result<()> {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom")?;
        let start_idx = schema.index_of("start")?;
        let ref_idx = schema.index_of("ref")?;
        let alt_idx = schema.index_of("alt")?;

        let chroms = get_string_col(batch.column(chrom_idx));
        let starts = get_u32_col(batch.column(start_idx));
        let refs = get_string_col(batch.column(ref_idx));
        let alts = get_string_col(batch.column(alt_idx));

        for row in 0..batch.num_rows() {
            self.stats.vcf_rows += 1;
            let raw_chrom = chroms[row].as_deref().unwrap_or("");
            let chrom = raw_chrom.strip_prefix("chr").unwrap_or(raw_chrom);
            let norm_start = starts[row] as i64 + 1;
            let norm_end = norm_start;

            self.ensure_window(chrom, norm_start)?;
            self.stats.total_window_rows_probed += 1;

            if let Some(index) = &self.current_index {
                let colocated = index.find_colocated(norm_start, norm_end);
                self.stats.total_window_rows_position_hit += colocated.len() as u64;

                let vcf_ref = refs[row].as_deref().unwrap_or("");
                let vcf_alt = alts[row].as_deref().unwrap_or("");
                let matches =
                    index.find_matches(norm_start, norm_end, vcf_ref, vcf_alt, allele_matches);
                self.stats.total_window_rows_allele_hit += matches.len() as u64;
                if !matches.is_empty() {
                    self.stats.matched_rows += 1;
                }
            }
        }

        Ok(())
    }
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <fjall_cache_path> <vcf_path>", args[0]);
        std::process::exit(1);
    }

    let cache_path = &args[1];
    let vcf_path = &args[2];

    let store = Arc::new(VepKvStore::open(cache_path)?);
    let window_size = store.window_size();
    println!(
        "Cache: {} columns, window_size={}",
        store.schema().fields().len(),
        window_size
    );
    println!("VCF: {}", vcf_path);

    let config = datafusion::prelude::SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);

    let table =
        VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
    ctx.register_table("vcf", Arc::new(table))?;

    let df = ctx
        .sql("SELECT `chrom`, `start`, `ref`, `alt` FROM vcf")
        .await?;
    let plan = df.create_physical_plan().await?;
    let task_ctx = ctx.task_ctx();

    let start_time = Instant::now();
    let mut stream = plan.execute(0, task_ctx)?;
    let mut state = AnnotationState::new(store);

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        state.annotate_batch(&batch)?;
    }

    let elapsed = start_time.elapsed().as_secs_f64();
    let s = &state.stats;

    let avg_window_rows = if s.window_loads > 0 {
        s.total_window_rows_deserialized as f64 / s.window_loads as f64
    } else {
        0.0
    };
    let vcf_per_window = if s.window_loads > 0 {
        s.vcf_rows as f64 / s.window_loads as f64
    } else {
        0.0
    };
    let utilization = if s.total_window_rows_deserialized > 0 {
        s.total_window_rows_probed as f64 / s.total_window_rows_deserialized as f64 * 100.0
    } else {
        0.0
    };

    println!("\n=== Utilization Stats ({elapsed:.1}s) ===");
    println!("VCF rows:                    {:>12}", s.vcf_rows);
    println!("Windows loaded:              {:>12}", s.window_loads);
    println!(
        "Cache rows deserialized:     {:>12}",
        s.total_window_rows_deserialized
    );
    println!(
        "VCF probes (position lookups):{:>11}",
        s.total_window_rows_probed
    );
    println!(
        "Position hits (colocated):   {:>12}",
        s.total_window_rows_position_hit
    );
    println!(
        "Allele hits (exact match):   {:>12}",
        s.total_window_rows_allele_hit
    );
    println!("VCF rows matched:            {:>12}", s.matched_rows);
    println!();
    println!("Avg cache rows/window:       {:>12.0}", avg_window_rows);
    println!("Avg VCF rows/window:         {:>12.1}", vcf_per_window);
    println!("Utilization (probed/deser):  {:>11.2}%", utilization);
    println!(
        "Amplification (deser/probed):{:>11.1}x",
        if s.total_window_rows_probed > 0 {
            s.total_window_rows_deserialized as f64 / s.total_window_rows_probed as f64
        } else {
            0.0
        }
    );
    println!(
        "Position selectivity:        {:>11.4}%",
        if s.total_window_rows_deserialized > 0 {
            s.total_window_rows_position_hit as f64 / s.total_window_rows_deserialized as f64
                * 100.0
        } else {
            0.0
        }
    );

    Ok(())
}
