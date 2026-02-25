//! Micro-benchmark: fjall sequential get latency for position-keyed lookups.
//!
//! Simulates the "flipped model" where each (chrom, position) is a KV entry.
//! Writes N entries with realistic value sizes, then reads them back sequentially
//! (mimicking sorted VCF position lookups) and measures per-get latency.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep-cache --example bench_fjall_gets \
//!     [num_entries] [value_size_bytes] [num_lookups] [cache_size_mb]
//!
//! Defaults: 4_000_000 entries, 4400 bytes/value, 50_000 lookups, 512 MB cache

use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let num_entries: usize = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(4_000_000);
    let value_size: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4400);
    let num_lookups: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(50_000);
    let cache_mb: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(512);

    eprintln!(
        "entries={num_entries} value_size={value_size}B lookups={num_lookups} cache={cache_mb}MB"
    );

    // Create temp directory
    let tmp = tempfile::tempdir().expect("tempdir");
    let db_path = tmp.path().join("bench_fjall");

    // --- WRITE PHASE ---
    eprintln!("\n=== WRITE PHASE ===");
    let write_start = Instant::now();
    {
        let db = fjall::Database::builder(&db_path)
            .cache_size(cache_mb * 1024 * 1024)
            .open()
            .expect("open database");
        let keyspace = db
            .keyspace("data", fjall::KeyspaceCreateOptions::default)
            .expect("open keyspace");

        // Generate a reproducible value payload
        let value_template: Vec<u8> = (0..value_size).map(|i| (i % 251) as u8).collect();

        // Write entries with keys = (chrom_code=0x0016, position as u64 BE)
        // Positions spaced ~25bp apart to mimic chr22 VEP density:
        //   chr22 has ~15M variants over ~35Mbp → ~2.3bp per variant → ~0.43 variants/bp
        //   With ~3.7 alleles per position → ~4M unique positions over 35Mbp → ~8.75bp spacing
        let chrom_prefix: [u8; 2] = [0x00, 0x16]; // chr22
        let position_step: u64 = 9; // ~9bp between positions

        for i in 0..num_entries {
            let pos = (i as u64) * position_step;
            let mut key = Vec::with_capacity(10);
            key.extend_from_slice(&chrom_prefix);
            key.extend_from_slice(&pos.to_be_bytes());

            // Vary value slightly so compression is realistic
            let mut value = value_template.clone();
            // Embed position in first 8 bytes (like real data would have varying content)
            value[0..8].copy_from_slice(&pos.to_be_bytes());

            keyspace.insert(&key, &value).expect("insert");

            if (i + 1) % 500_000 == 0 {
                eprintln!("  wrote {}/{num_entries}", i + 1);
            }
        }

        db.persist(fjall::PersistMode::SyncAll).expect("persist");
    }
    let write_elapsed = write_start.elapsed();
    eprintln!(
        "write: {:.2}s ({:.0} entries/s)",
        write_elapsed.as_secs_f64(),
        num_entries as f64 / write_elapsed.as_secs_f64()
    );

    // Measure disk usage
    let disk_bytes: u64 = walkdir(&db_path);
    eprintln!(
        "disk: {:.1}MB ({:.0} bytes/entry)",
        disk_bytes as f64 / 1_048_576.0,
        disk_bytes as f64 / num_entries as f64
    );

    // --- READ PHASE ---
    eprintln!("\n=== READ PHASE ===");

    // Re-open (simulates cold open, then cache warms up)
    let db = fjall::Database::builder(&db_path)
        .cache_size(cache_mb * 1024 * 1024)
        .open()
        .expect("open database for read");
    let keyspace = db
        .keyspace("data", fjall::KeyspaceCreateOptions::default)
        .expect("open keyspace for read");

    let chrom_prefix: [u8; 2] = [0x00, 0x16];
    let position_step: u64 = 9;
    let _max_pos = (num_entries as u64) * position_step;

    // Generate lookup positions — sequential subset (like sorted VCF)
    // Pick `num_lookups` evenly spaced positions that exist in the DB
    let lookup_step = if num_lookups > 0 {
        num_entries / num_lookups
    } else {
        1
    };
    let lookup_indices: Vec<usize> = (0..num_lookups).map(|i| i * lookup_step).collect();

    // --- Warmup: read 1000 entries to warm block cache ---
    eprintln!("warming up block cache...");
    for &idx in lookup_indices.iter().take(1000) {
        let pos = (idx as u64) * position_step;
        let mut key = Vec::with_capacity(10);
        key.extend_from_slice(&chrom_prefix);
        key.extend_from_slice(&pos.to_be_bytes());
        let _ = keyspace.get(&key).expect("warmup get");
    }

    // --- Benchmark 1: Sequential gets (sorted VCF pattern) ---
    for iter in 1..=3 {
        let start = Instant::now();
        let mut total_bytes: u64 = 0;
        let mut found: u64 = 0;
        let mut not_found: u64 = 0;

        for &idx in &lookup_indices {
            let pos = (idx as u64) * position_step;
            let mut key = Vec::with_capacity(10);
            key.extend_from_slice(&chrom_prefix);
            key.extend_from_slice(&pos.to_be_bytes());

            match keyspace.get(&key) {
                Ok(Some(val)) => {
                    total_bytes += val.len() as u64;
                    found += 1;
                }
                Ok(None) => {
                    not_found += 1;
                }
                Err(e) => panic!("get failed: {e}"),
            }
        }

        let elapsed = start.elapsed();
        let per_get_us = elapsed.as_secs_f64() * 1_000_000.0 / num_lookups as f64;
        let throughput_mbs = total_bytes as f64 / 1_048_576.0 / elapsed.as_secs_f64();

        eprintln!(
            "sequential iter={iter}: {:.3}s  {per_get_us:.1}μs/get  found={found} miss={not_found}  {throughput_mbs:.0}MB/s  ({:.0} gets/s)",
            elapsed.as_secs_f64(),
            num_lookups as f64 / elapsed.as_secs_f64()
        );
    }

    // --- Benchmark 2: Random gets (worst case) ---
    eprintln!("\n--- random access pattern ---");
    // Shuffle positions deterministically
    let mut random_indices = lookup_indices.clone();
    // Simple deterministic shuffle (LCG-based)
    let n = random_indices.len();
    let mut seed: u64 = 12345;
    for i in (1..n).rev() {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let j = (seed >> 33) as usize % (i + 1);
        random_indices.swap(i, j);
    }

    for iter in 1..=3 {
        let start = Instant::now();
        let mut total_bytes: u64 = 0;
        let mut found: u64 = 0;

        for &idx in &random_indices {
            let pos = (idx as u64) * position_step;
            let mut key = Vec::with_capacity(10);
            key.extend_from_slice(&chrom_prefix);
            key.extend_from_slice(&pos.to_be_bytes());

            if let Ok(Some(val)) = keyspace.get(&key) {
                total_bytes += val.len() as u64;
                found += 1;
            }
        }

        let elapsed = start.elapsed();
        let per_get_us = elapsed.as_secs_f64() * 1_000_000.0 / num_lookups as f64;
        let throughput_mbs = total_bytes as f64 / 1_048_576.0 / elapsed.as_secs_f64();

        eprintln!(
            "random   iter={iter}: {:.3}s  {per_get_us:.1}μs/get  found={found}  {throughput_mbs:.0}MB/s  ({:.0} gets/s)",
            elapsed.as_secs_f64(),
            num_lookups as f64 / elapsed.as_secs_f64()
        );
    }

    // --- Benchmark 3: Batch sequential with pre-built keys (eliminate key construction overhead) ---
    eprintln!("\n--- pre-built keys (no alloc in loop) ---");
    let prebuilt_keys: Vec<Vec<u8>> = lookup_indices
        .iter()
        .map(|&idx| {
            let pos = (idx as u64) * position_step;
            let mut key = Vec::with_capacity(10);
            key.extend_from_slice(&chrom_prefix);
            key.extend_from_slice(&pos.to_be_bytes());
            key
        })
        .collect();

    for iter in 1..=3 {
        let start = Instant::now();
        let mut total_bytes: u64 = 0;
        let mut found: u64 = 0;

        for key in &prebuilt_keys {
            if let Ok(Some(val)) = keyspace.get(key) {
                total_bytes += val.len() as u64;
                found += 1;
            }
        }

        let elapsed = start.elapsed();
        let per_get_us = elapsed.as_secs_f64() * 1_000_000.0 / num_lookups as f64;
        let throughput_mbs = total_bytes as f64 / 1_048_576.0 / elapsed.as_secs_f64();

        eprintln!(
            "prebuilt iter={iter}: {:.3}s  {per_get_us:.1}μs/get  found={found}  {throughput_mbs:.0}MB/s  ({:.0} gets/s)",
            elapsed.as_secs_f64(),
            num_lookups as f64 / elapsed.as_secs_f64()
        );
    }

    // --- Benchmark 4: Vary value sizes ---
    // This measures overhead per get independent of value size
    // Use a smaller DB with 100-byte values
    eprintln!("\n--- small values (100B) for overhead measurement ---");
    let small_db_path = tmp.path().join("bench_fjall_small");
    {
        let db_small = fjall::Database::builder(&small_db_path)
            .cache_size(cache_mb * 1024 * 1024)
            .open()
            .expect("open small db");
        let ks_small = db_small
            .keyspace("data", fjall::KeyspaceCreateOptions::default)
            .expect("open small keyspace");
        let small_val: Vec<u8> = vec![42u8; 100];
        for i in 0..num_lookups {
            let pos = (i as u64) * position_step;
            let mut key = Vec::with_capacity(10);
            key.extend_from_slice(&chrom_prefix);
            key.extend_from_slice(&pos.to_be_bytes());
            ks_small.insert(&key, &small_val).expect("insert small");
        }
        db_small
            .persist(fjall::PersistMode::SyncAll)
            .expect("persist small");
    }
    let db_small = fjall::Database::builder(&small_db_path)
        .cache_size(cache_mb * 1024 * 1024)
        .open()
        .expect("reopen small db");
    let ks_small = db_small
        .keyspace("data", fjall::KeyspaceCreateOptions::default)
        .expect("reopen small keyspace");

    // Warmup
    for key in prebuilt_keys.iter().take(1000) {
        let _ = ks_small.get(key);
    }

    for iter in 1..=3 {
        let start = Instant::now();
        let mut found: u64 = 0;
        for key in &prebuilt_keys {
            if let Ok(Some(_val)) = ks_small.get(key) {
                found += 1;
            }
        }
        let elapsed = start.elapsed();
        let per_get_us = elapsed.as_secs_f64() * 1_000_000.0 / num_lookups as f64;
        eprintln!(
            "small    iter={iter}: {:.3}s  {per_get_us:.1}μs/get  found={found}  ({:.0} gets/s)",
            elapsed.as_secs_f64(),
            num_lookups as f64 / elapsed.as_secs_f64()
        );
    }

    // Summary
    eprintln!("\n=== SUMMARY ===");
    eprintln!(
        "DB: {num_entries} entries × {value_size}B values = {:.1}GB uncompressed payload",
        num_entries as f64 * value_size as f64 / 1e9
    );
    eprintln!("Disk: {:.1}MB", disk_bytes as f64 / 1_048_576.0);
    eprintln!("Lookups: {num_lookups} (of {num_entries} entries)");
}

fn walkdir(path: &std::path::Path) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = std::fs::read_dir(path) {
        for entry in entries.flatten() {
            let meta = entry.metadata();
            if let Ok(m) = meta {
                if m.is_dir() {
                    total += walkdir(&entry.path());
                } else {
                    total += m.len();
                }
            }
        }
    }
    total
}
