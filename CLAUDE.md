# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

`datafusion-bio-functions` is a Rust library providing reusable Apache DataFusion UDFs (User-Defined Functions) for bioinformatics. It is the companion to `datafusion-bio-formats` (which handles bioinformatic file format readers). Both are consumed by `polars-bio` (Python genomics library) and must keep DataFusion/Arrow versions in sync.

**Ecosystem:**
- `datafusion-bio-formats` — DataFusion table providers for VCF, BAM, CRAM, BED, FASTA, FASTQ, GFF, PAIRS
- `datafusion-bio-functions` (this repo) — DataFusion scalar/aggregate/window functions for bioinformatics
- `polars-bio` — Python library (PyO3) combining both, plus interval join algorithms from `sequila-core`
- `comet-bio` — Spark extension using similar Rust core logic via JNI

## Dependency Versions (must stay in sync with polars-bio)

As of polars-bio latest:
- **DataFusion**: 50.3.0
- **Arrow**: 56.1.0
- **Rust edition**: 2024

When updating DataFusion or Arrow versions, coordinate with `datafusion-bio-formats` and `polars-bio` to ensure all three repos use the same major versions.

## Build Commands

```bash
cargo build                    # Debug build
cargo build --release          # Release build
cargo test                     # Run all tests
cargo test <test_name>         # Run a single test
cargo test -- --nocapture      # Run tests with stdout visible
cargo fmt                      # Format code
cargo fmt -- --check           # Check formatting without modifying
cargo clippy                   # Lint
cargo clippy -- -D warnings    # Lint, treat warnings as errors
```

## Architecture

This is a Cargo workspace (or single crate — check `Cargo.toml`). Functions should be registered as DataFusion UDFs following the `datafusion::logical_expr::ScalarUDF` / `AggregateUDF` / `WindowUDF` patterns.

**Key DataFusion traits to implement:**
- `ScalarUDFImpl` — for row-level functions (e.g., `reverse_complement()`, `gc_content()`, `edit_distance()`)
- `AggregateUDFImpl` — for aggregate functions
- `WindowUDFImpl` — for window functions

Functions operate on Arrow arrays (`arrow::array`) and return Arrow arrays. Use Arrow compute kernels where possible for performance.

## Conventions from Sibling Repos

- Apache 2.0 license
- Git dependencies for cross-repo references use pinned commit hashes (e.g., `git = "https://github.com/biodatageeks/...", rev = "abc123"`)
- `polars-bio` will reference this repo the same way it references `datafusion-bio-formats`
- Mutation testing is set up (cargo-mutants) — `mutants.out*/` is gitignored
- IDE: RustRover/IntelliJ (`.idea/` is gitignored)
