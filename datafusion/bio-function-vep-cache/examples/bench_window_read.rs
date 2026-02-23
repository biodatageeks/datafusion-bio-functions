use std::sync::Arc;
use std::time::Instant;

use datafusion::common::Result;

use datafusion_bio_function_vep_cache::allele_index::WindowAlleleIndex;
use datafusion_bio_function_vep_cache::kv_store::{FORMAT_V1, VepKvStore};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <fjall_cache_path> [num_windows]", args[0]);
        std::process::exit(1);
    }

    let cache_path = &args[1];
    let num_windows: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(100);

    let store = Arc::new(VepKvStore::open(cache_path)?);
    let fmt_ver = store.format_version();
    println!(
        "Cache: {} columns, window_size={}, format=v{}",
        store.schema().fields().len(),
        store.window_size(),
        fmt_ver,
    );

    if fmt_ver >= FORMAT_V1 {
        bench_v1(&store, num_windows)
    } else {
        bench_v0(&store, num_windows)
    }
}

fn bench_v1(store: &VepKvStore, num_windows: u64) -> Result<()> {
    // Time position index reads
    let start = Instant::now();
    let mut read_count = 0u32;
    let mut total_rows = 0usize;
    for wid in 0..num_windows {
        if let Some(idx) = store.get_position_index("22", wid)? {
            total_rows += idx.num_rows();
            read_count += 1;
        }
    }
    let read_elapsed = start.elapsed();
    println!(
        "Position index read: {} windows, {} rows, {:.1}ms total, {:.2}ms/window",
        read_count,
        total_rows,
        read_elapsed.as_secs_f64() * 1000.0,
        if read_count > 0 {
            read_elapsed.as_secs_f64() * 1000.0 / read_count as f64
        } else {
            0.0
        }
    );

    // Time position index + WindowAlleleIndex build
    let start = Instant::now();
    let mut idx_count = 0u32;
    for wid in 0..num_windows {
        if let Some(idx) = store.get_position_index("22", wid)? {
            let _wai = WindowAlleleIndex::from_position_index(idx);
            idx_count += 1;
        }
    }
    let idx_elapsed = start.elapsed();
    println!(
        "Read + index: {} windows, {:.1}ms total, {:.2}ms/window",
        idx_count,
        idx_elapsed.as_secs_f64() * 1000.0,
        if idx_count > 0 {
            idx_elapsed.as_secs_f64() * 1000.0 / idx_count as f64
        } else {
            0.0
        }
    );

    // Time column read (first 3 columns, first N windows)
    let start = Instant::now();
    let mut col_count = 0u32;
    let mut col_total_rows = 0usize;
    for wid in 0..num_windows {
        if let Some(cols) = store.get_columns("22", wid, &[0, 1, 2])? {
            col_total_rows += cols[0].0.len();
            col_count += 1;
        }
    }
    let col_elapsed = start.elapsed();
    println!(
        "Column read (3 cols): {} windows, {} rows, {:.1}ms total, {:.2}ms/window",
        col_count,
        col_total_rows,
        col_elapsed.as_secs_f64() * 1000.0,
        if col_count > 0 {
            col_elapsed.as_secs_f64() * 1000.0 / col_count as f64
        } else {
            0.0
        }
    );

    Ok(())
}

fn bench_v0(store: &VepKvStore, num_windows: u64) -> Result<()> {
    // Time raw fjall reads
    let start = Instant::now();
    let mut read_count = 0u32;
    let mut total_rows = 0usize;
    let mut total_bytes = 0usize;
    for wid in 0..num_windows {
        if let Some(batch) = store.get_window("22", wid)? {
            total_rows += batch.num_rows();
            total_bytes += batch.get_array_memory_size();
            read_count += 1;
        }
    }
    let read_elapsed = start.elapsed();
    println!(
        "Raw read: {} windows, {} rows, {:.1} MB, {:.1}ms total, {:.2}ms/window",
        read_count,
        total_rows,
        total_bytes as f64 / 1_048_576.0,
        read_elapsed.as_secs_f64() * 1000.0,
        if read_count > 0 {
            read_elapsed.as_secs_f64() * 1000.0 / read_count as f64
        } else {
            0.0
        }
    );

    // Time index build
    let start = Instant::now();
    let mut idx_count = 0u32;
    for wid in 0..num_windows {
        if let Some(batch) = store.get_window("22", wid)? {
            let _idx = WindowAlleleIndex::from_batch(batch)?;
            idx_count += 1;
        }
    }
    let idx_elapsed = start.elapsed();
    println!(
        "Read + index: {} windows, {:.1}ms total, {:.2}ms/window",
        idx_count,
        idx_elapsed.as_secs_f64() * 1000.0,
        if idx_count > 0 {
            idx_elapsed.as_secs_f64() * 1000.0 / idx_count as f64
        } else {
            0.0
        }
    );

    // Index build only (read is cached by fjall)
    let overhead = idx_elapsed.as_secs_f64() * 1000.0 - read_elapsed.as_secs_f64() * 1000.0;
    println!(
        "Index build overhead: {:.1}ms total, {:.2}ms/window",
        overhead,
        if idx_count > 0 {
            overhead / idx_count as f64
        } else {
            0.0
        }
    );

    Ok(())
}
