//! Build an exact per-chromosome position-existence index for variation Fjall lookups.
//!
//! Example:
//!   build_variation_position_index --cache /path/to/cache --chrom chr1 \
//!     --output-dir /tmp/variation.position_index

use std::path::{Path, PathBuf};

use datafusion_bio_function_vep::kv_cache::position_index::{PositionIndex, position_index_file};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Debug)]
struct Args {
    cache: PathBuf,
    chrom: String,
    cold_parquet: Option<PathBuf>,
    output_dir: Option<PathBuf>,
    batch_size: usize,
}

fn main() -> Result<()> {
    let args = parse_args()?;
    let input = args
        .cold_parquet
        .clone()
        .or_else(|| cold_variation_file_for_chrom(&args.cache, &args.chrom))
        .ok_or_else(|| format!("could not find variation parquet for {}", args.chrom))?;
    let output_dir = args
        .output_dir
        .clone()
        .unwrap_or_else(|| args.cache.join("variation.position_index"));
    let output = position_index_file(
        &output_dir,
        normalized_chrom_for_index(&args.chrom, &input)?,
    );

    eprintln!("input:  {}", input.display());
    eprintln!("output: {}", output.display());

    let index = PositionIndex::from_parquet(&input, args.batch_size)?;
    index.write_to_path(&output)?;

    eprintln!("wrote {} unique positions", index.len());
    Ok(())
}

fn parse_args() -> Result<Args> {
    let mut cache = None;
    let mut chrom = None;
    let mut cold_parquet = None;
    let mut output_dir = None;
    let mut batch_size = 65_536_usize;

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--cache" => cache = args.next().map(PathBuf::from),
            "--chrom" => chrom = args.next(),
            "--cold-parquet" => cold_parquet = args.next().map(PathBuf::from),
            "--output-dir" => output_dir = args.next().map(PathBuf::from),
            "--batch-size" => {
                batch_size = args
                    .next()
                    .ok_or("--batch-size requires a value")?
                    .parse()?
            }
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}").into()),
        }
    }

    let Some(cache) = cache else {
        print_usage();
        return Err("--cache is required".into());
    };
    let Some(chrom) = chrom else {
        print_usage();
        return Err("--chrom is required".into());
    };

    Ok(Args {
        cache,
        chrom,
        cold_parquet,
        output_dir,
        batch_size,
    })
}

fn print_usage() {
    eprintln!(
        "Usage: build_variation_position_index --cache <cache-dir> --chrom <chrN> \
         [--cold-parquet <variation/chrN_cold.parquet>] [--output-dir <dir>] [--batch-size 65536]"
    );
}

fn cold_variation_file_for_chrom(cache: &Path, chrom: &str) -> Option<PathBuf> {
    let variation_dir = cache.join("variation");
    let direct = variation_dir.join(format!("{chrom}_cold.parquet"));
    if direct.is_file() {
        return Some(direct);
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        let stripped = variation_dir.join(format!("{stripped}_cold.parquet"));
        if stripped.is_file() {
            return Some(stripped);
        }
    } else {
        let prefixed = variation_dir.join(format!("chr{chrom}_cold.parquet"));
        if prefixed.is_file() {
            return Some(prefixed);
        }
    }

    None
}

fn normalized_chrom_for_index<'a>(chrom: &'a str, path: &'a Path) -> Result<&'a str> {
    if !chrom.is_empty() {
        return Ok(chrom);
    }

    let stem = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or_else(|| format!("invalid chromosome parquet file name: {}", path.display()))?;
    Ok(stem
        .strip_suffix("_cold")
        .or_else(|| stem.strip_suffix("_warm"))
        .unwrap_or(stem))
}
