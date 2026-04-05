//! Compact an existing fjall KV store in-place.
//!
//! Use this to fix uncompacted stores (all L0 SSTs, no merged levels) that
//! result from interrupted builds or missing major_compact() calls.
//!
//! Usage:
//!   cargo run --release --example compact_fjall --features kv-cache -- \
//!     <fjall_path> [keyspace_name]
//!
//! Examples:
//!   # Compact all keyspaces in a variation store
//!   compact_fjall /data/vep/variation.fjall
//!
//!   # Compact only the sift keyspace in a translation_sift store
//!   compact_fjall /data/vep/translation_sift.fjall sift

use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <fjall_path> [keyspace_name]", args[0]);
        eprintln!();
        eprintln!("Compacts all keyspaces in the fjall database, or a specific one.");
        eprintln!("This merges L0 SSTs into sorted runs for optimal read performance.");
        std::process::exit(1);
    }

    let fjall_path = &args[1];
    let target_keyspace = args.get(2).map(|s| s.as_str());

    eprintln!("Opening fjall database: {fjall_path}");
    let db = fjall::Database::builder(fjall_path)
        .cache_size(256 * 1024 * 1024)
        .open()
        .expect("failed to open fjall database");

    // Well-known keyspace names used across the VEP cache stores.
    let known_keyspaces = ["data", "meta", "sift", "translations", "exons"];

    for ks_name in &known_keyspaces {
        if target_keyspace.is_some_and(|t| *ks_name != t) {
            continue;
        }

        // Guard: only open keyspaces that already exist. db.keyspace() with
        // KeyspaceCreateOptions::default creates the keyspace if missing —
        // that would introduce the extra keyspaces this PR removes (~10% perf hit).
        if !db.keyspace_exists(ks_name) {
            continue;
        }

        let ks = match db.keyspace(ks_name, fjall::KeyspaceCreateOptions::default) {
            Ok(ks) => ks,
            Err(_) => continue,
        };

        eprintln!("Compacting keyspace '{ks_name}'...");
        let start = Instant::now();
        if let Err(e) = ks.major_compact() {
            eprintln!("  ERROR compacting '{ks_name}': {e}");
            continue;
        }
        eprintln!(
            "  '{ks_name}' compacted in {:.1}s",
            start.elapsed().as_secs_f64()
        );
    }

    eprintln!("Persisting...");
    db.persist(fjall::PersistMode::SyncAll)
        .expect("failed to persist");

    // Skip Drop to avoid deadlock (issue #86).
    std::mem::forget(db);

    eprintln!("Done.");
}
