//! Diagnose unused keyspaces in fjall databases within a VEP cache directory.
//!
//! Extra keyspaces (translations, exons, sift) in variation.fjall and
//! translation_sift.fjall degrade point-lookup performance by ~10% due
//! to fjall block cache overhead from unused keyspace metadata.
//!
//! This tool checks for extra keyspaces and reports what needs cleaning.
//! Since fjall does not support dropping keyspaces, the recommended fix
//! is to rebuild the cache from scratch using `cache_builder`.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep --example strip_fjall \
//!     --features kv-cache -- <cache_dir>

use std::path::Path;

fn check_keyspaces(db_path: &Path, expected: &[&str]) {
    let ks_dir = db_path.join("keyspaces");
    if !ks_dir.is_dir() {
        eprintln!("  {} — not a fjall database, skipping", db_path.display());
        return;
    }

    let db = match fjall::Database::builder(db_path)
        .cache_size(64 * 1024 * 1024)
        .worker_threads(1)
        .open()
    {
        Ok(db) => db,
        Err(e) => {
            eprintln!("  failed to open {}: {e}", db_path.display());
            return;
        }
    };

    let all_names = ["meta", "data", "sift", "translations", "exons"];
    let mut found_expected = Vec::new();
    let mut found_extra = Vec::new();

    for name in &all_names {
        if !db.keyspace_exists(name) {
            continue;
        }
        let ks = match db.keyspace(name, fjall::KeyspaceCreateOptions::default) {
            Ok(ks) => ks,
            Err(_) => continue,
        };
        let tables = ks.table_count();
        let empty = ks.is_empty().unwrap_or(true);
        let disk = ks.disk_space();

        if expected.contains(name) {
            found_expected.push(format!(
                "  ✓ '{}': {} tables, {:.2} GB",
                name,
                tables,
                disk as f64 / 1e9
            ));
        } else {
            found_extra.push(format!(
                "  ✗ '{}': {} tables, empty={} (UNUSED — causes ~10% perf degradation)",
                name, tables, empty
            ));
        }
    }

    // Count numeric keyspace dirs on disk
    let mut dir_count = 0;
    if let Ok(entries) = std::fs::read_dir(&ks_dir) {
        for entry in entries.flatten() {
            if entry.file_name().to_string_lossy().parse::<u32>().is_ok() {
                dir_count += 1;
            }
        }
    }

    eprintln!(
        "  {} ({} keyspace dirs on disk)",
        db_path.display(),
        dir_count
    );
    for line in &found_expected {
        eprintln!("{line}");
    }
    for line in &found_extra {
        eprintln!("{line}");
    }

    if found_extra.is_empty() {
        eprintln!("  → OK: no extra keyspaces");
    } else {
        eprintln!(
            "  → REBUILD RECOMMENDED: {} extra keyspace(s) found",
            found_extra.len()
        );
        eprintln!("    Rebuild with: cache_builder (creates only needed keyspaces)");
    }

    std::mem::forget(db);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <cache_dir>", args[0]);
        eprintln!();
        eprintln!("Checks for unused keyspaces in variation.fjall and translation_sift.fjall.");
        eprintln!("Extra keyspaces cause ~10% annotation performance degradation.");
        eprintln!("Fix: rebuild the cache from scratch using cache_builder.");
        std::process::exit(1);
    }

    let cache_dir = Path::new(&args[1]);

    eprintln!("Checking variation.fjall...");
    let var_fjall = cache_dir.join("variation.fjall");
    if var_fjall.is_dir() {
        check_keyspaces(&var_fjall, &["meta", "data"]);
    } else {
        eprintln!("  not found, skipping");
    }

    eprintln!("\nChecking translation_sift.fjall...");
    let sift_fjall = cache_dir.join("translation_sift.fjall");
    if sift_fjall.is_dir() {
        check_keyspaces(&sift_fjall, &["meta", "sift"]);
    } else {
        eprintln!("  not found, skipping");
    }
}
