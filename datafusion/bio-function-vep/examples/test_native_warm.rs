use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;
use rayon::prelude::*;

const IO_BUFFER_SIZE: usize = 64 * 1024;

#[derive(Debug, Clone)]
struct BgzfPartition {
    start_compressed: u64,
    end_compressed: u64,
    skip_first_line: bool,
}

fn scan_bgzf_block_offsets(file_path: &Path, file_size: u64) -> Result<Vec<u64>, Box<dyn std::error::Error>> {
    use std::io::{Read, Seek, SeekFrom};
    let mut file = File::open(file_path)?;
    let mut offsets = Vec::new();
    let mut pos = 0u64;
    let mut header = [0u8; 18];

    while pos + 28 <= file_size {
        file.seek(SeekFrom::Start(pos))?;
        if file.read_exact(&mut header).is_err() {
            break;
        }
        if header[0] != 0x1f || header[1] != 0x8b || header[2] != 0x08 || header[3] != 0x04 {
            break;
        }
        let bsize = u16::from_le_bytes([header[16], header[17]]) as u64 + 1;
        let is_eof_marker = bsize == 28 && pos + bsize >= file_size;
        if !is_eof_marker {
            offsets.push(pos);
        }
        pos += bsize;
    }
    Ok(offsets)
}

fn compute_bgzf_partitions(file_path: &Path, num_partitions: usize) -> Result<Vec<BgzfPartition>, Box<dyn std::error::Error>> {
    let file_size = std::fs::metadata(file_path)?.len();
    if file_size == 0 || num_partitions <= 1 {
        return Ok(vec![BgzfPartition {
            start_compressed: 0,
            end_compressed: file_size,
            skip_first_line: false,
        }]);
    }

    let block_offsets = scan_bgzf_block_offsets(file_path, file_size)?;
    if block_offsets.len() <= 1 {
        return Ok(vec![BgzfPartition {
            start_compressed: 0,
            end_compressed: file_size,
            skip_first_line: false,
        }]);
    }

    let target_chunk = file_size / num_partitions as u64;
    let mut boundaries = vec![0u64];
    let mut current_start = 0u64;

    for i in 1..num_partitions {
        let target_offset = target_chunk * i as u64;
        let best = block_offsets
            .iter()
            .copied()
            .filter(|&off| off > current_start)
            .min_by_key(|&off| (off as i64 - target_offset as i64).unsigned_abs())
            .unwrap_or(file_size);

        if best < file_size && best > current_start {
            boundaries.push(best);
            current_start = best;
        }
    }

    // Check skip_first_line for each boundary
    let mut partitions = Vec::with_capacity(boundaries.len());
    for i in 0..boundaries.len() {
        let start = boundaries[i];
        let end = if i + 1 < boundaries.len() {
            boundaries[i + 1]
        } else {
            file_size
        };

        let skip_first_line = if start > 0 {
            // Check if preceding block ended with newline
            let prev_block_offset = block_offsets
                .iter()
                .copied()
                .filter(|&off| off < start)
                .last()
                .unwrap_or(0);
            !block_ends_with_newline(file_path, prev_block_offset)?
        } else {
            false
        };

        partitions.push(BgzfPartition {
            start_compressed: start,
            end_compressed: end,
            skip_first_line,
        });
    }

    Ok(partitions)
}

fn block_ends_with_newline(file_path: &Path, block_offset: u64) -> Result<bool, Box<dyn std::error::Error>> {
    use std::io::Read;
    let file = File::open(file_path)?;
    let mut reader = noodles_bgzf::io::Reader::new(BufReader::with_capacity(IO_BUFFER_SIZE, file));
    let vpos = noodles_bgzf::VirtualPosition::try_from((block_offset, 0))?;
    reader.seek(vpos)?;

    let mut buf = [0u8; 64 * 1024];
    let mut last_byte = None;
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        last_byte = Some(buf[n - 1]);
        if reader.virtual_position().compressed() != block_offset {
            break;
        }
    }
    Ok(last_byte == Some(b'\n'))
}

/// Fast extract max AF from a string like "A:0.05,C:0.01" or "0.05"
#[inline]
fn parse_max_af_field(s: &str) -> f64 {
    if s.is_empty() || s == "." {
        return 0.0;
    }
    let mut max_val = 0.0f64;
    for part in s.split(',') {
        let val_str = if let Some((_, v)) = part.split_once(':') {
            v
        } else {
            part
        };
        if let Ok(v) = val_str.parse::<f64>() {
            if v.is_finite() && v > max_val {
                max_val = v;
            }
        }
    }
    max_val
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new("/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/cache_bench_dev/prepared_cache/20/all_vars.gz");
    println!("Testing fast native warm starts on: {}", path.display());

    let t0 = Instant::now();
    let num_threads = 10;
    let partitions = compute_bgzf_partitions(path, num_threads)?;
    println!("Computed {} partitions in {:.2?}", partitions.len(), t0.elapsed());

    let t_scan = Instant::now();
    let results: Vec<(Vec<i64>, Vec<i64>)> = partitions
        .into_par_iter()
        .map(|part| {
            let file = File::open(path).expect("open file");
            let mut reader = noodles_bgzf::io::Reader::new(BufReader::with_capacity(IO_BUFFER_SIZE, file));
            if part.start_compressed > 0 {
                let vpos = noodles_bgzf::VirtualPosition::try_from((part.start_compressed, 0)).expect("vpos");
                reader.seek(vpos).expect("seek");
                if part.skip_first_line {
                    let mut discard = String::new();
                    reader.read_line(&mut discard).expect("skip line");
                }
            }

            let mut all_starts = Vec::with_capacity(3_000_000);
            let mut high_af_starts = Vec::with_capacity(100_000);
            let mut line = String::new();

            loop {
                if reader.virtual_position().compressed() >= part.end_compressed {
                    break;
                }
                line.clear();
                match reader.read_line(&mut line) {
                    Ok(0) => break,
                    Ok(_) => {
                        let trimmed = line.trim();
                        if trimmed.is_empty() || trimmed.starts_with('#') {
                            continue;
                        }
                        // Tab indices:
                        // 0: chr, 1: id, 2: failed, 3: somatic, 4: start, ...
                        // 10: minor_allele_freq
                        // 15: AF
                        // 21: gnomADe
                        // 31: gnomADg
                        let mut tab_idx = 0;
                        let mut start_val: i64 = 0;
                        let mut max_af = 0.0f64;

                        for field in trimmed.split('\t') {
                            match tab_idx {
                                4 => {
                                    if let Ok(v) = field.parse::<i64>() {
                                        start_val = v;
                                    }
                                }
                                10 => {
                                    if let Ok(v) = field.parse::<f64>() {
                                        if v.is_finite() && v > max_af {
                                            max_af = v;
                                        }
                                    }
                                }
                                15 | 21 | 31 => {
                                    let v = parse_max_af_field(field);
                                    if v > max_af {
                                        max_af = v;
                                    }
                                }
                                _ => {}
                            }
                            tab_idx += 1;
                            if tab_idx > 31 {
                                break;
                            }
                        }

                        if start_val > 0 {
                            all_starts.push(start_val);
                            if max_af >= 0.01 {
                                high_af_starts.push(start_val);
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("read error: {e}");
                        break;
                    }
                }
            }
            (all_starts, high_af_starts)
        })
        .collect();

    println!("Scanned all partitions in {:.2?}", t_scan.elapsed());

    let t_set = Instant::now();
    let mut total_rows = 0;
    let mut source_positions: HashSet<i64> = HashSet::new();
    let mut high_af: Vec<i64> = Vec::new();

    for (starts, highs) in results {
        total_rows += starts.len();
        source_positions.extend(starts);
        high_af.extend(highs);
    }

    let mut warm_set: HashSet<i64> = HashSet::new();
    for start in high_af {
        for offset in -1..=1 {
            let cand = start + offset;
            if source_positions.contains(&cand) {
                warm_set.insert(cand);
            }
        }
    }

    println!("Total rows scanned: {}", total_rows);
    println!("Unique source positions: {}", source_positions.len());
    println!("Warm positions count: {}", warm_set.len());
    println!("Post-processing took: {:.2?}", t_set.elapsed());
    println!("TOTAL TIME: {:.2?}", t0.elapsed());

    Ok(())
}
