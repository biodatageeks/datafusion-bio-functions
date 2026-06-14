use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use datafusion_bio_function_vep::allele::{vcf_to_vep_allele, vcf_to_vep_input_allele};
use datafusion_bio_function_vep::kv_cache::key_encoding::chrom_to_code;
use datafusion_bio_function_vep::kv_cache::position_index::{
    PositionIndex, find_position_index_file,
};
use datafusion_bio_function_vep::kv_cache::variant_bloom_index::{
    VariantBloomIndex, find_variant_bloom_index_file,
};
use datafusion_bio_function_vep::warm_cache::key::{
    position_key_from_code, variant_key_from_position,
};
use flate2::read::MultiGzDecoder;
use futures::TryStreamExt;
use lance::Dataset;
use lance::dataset::scanner::MaterializationStyle;

use crate::build::key_values_from_batch;
use crate::config::{KeyDataType, KeyMode, SandboxConfig};

#[derive(Debug)]
struct Workload {
    ordered_positions: Vec<u64>,
    variant_keys_by_position: HashMap<u64, Vec<u64>>,
}

#[derive(Debug)]
pub struct SampleColdReport {
    pub unique_cold_positions: usize,
    pub files: Vec<SampleColdFile>,
}

#[derive(Debug)]
pub struct SampleColdFile {
    pub limit: usize,
    pub path: PathBuf,
    pub rows: usize,
}

pub async fn extract_keys(config: &SandboxConfig, output_path: &Path) -> Result<()> {
    let chrom = config.dataset.chrom.as_str();
    let chrom_code = chrom_to_code(chrom.strip_prefix("chr").unwrap_or(chrom));
    let workload = load_vcf_workload(&config.benchmark.vcf_path, chrom)?;
    let warm_positions = load_warm_position_keys(config, chrom_code).await?;

    let position_index_path = find_position_index_file(
        config.dataset.cache_root.join("variation.position_index"),
        chrom,
    )
    .with_context(|| format!("missing position index for {chrom}"))?;
    let bloom_index_path = find_variant_bloom_index_file(
        config
            .dataset
            .cache_root
            .join("variation.variant_bloom_index"),
        chrom,
    )
    .with_context(|| format!("missing variant bloom index for {chrom}"))?;
    let position_index = PositionIndex::read_from_path(&position_index_path)?;
    let bloom_index = VariantBloomIndex::read_from_path(&bloom_index_path)?;

    let mut selected = Vec::new();
    for &position_key in &workload.ordered_positions {
        if warm_positions.contains(&position_key) {
            continue;
        }
        if !position_index.contains_position_key(position_key) {
            continue;
        }
        let variant_keys = workload
            .variant_keys_by_position
            .get(&position_key)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let variant_admitted = bloom_index.contains_any_variant_keys(variant_keys.iter().copied());
        let colocated_admitted = if bloom_index.supports_position_fallback_keys() {
            bloom_index.contains_position_fallback_key(position_key)
        } else {
            true
        };
        if variant_admitted || colocated_admitted {
            selected.push(match config.key.mode {
                KeyMode::Position => position_from_key(position_key),
                KeyMode::PositionKey => position_key,
            });
        }
    }
    selected.sort_unstable();
    selected.dedup();

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut file = File::create(output_path)
        .with_context(|| format!("failed to create '{}'", output_path.display()))?;
    for value in selected {
        writeln!(file, "{value}")?;
    }
    Ok(())
}

pub async fn sample_cold_positions(
    config: &SandboxConfig,
    output_dir: &Path,
    limits: &[usize],
) -> Result<SampleColdReport> {
    if limits.is_empty() {
        bail!("sample-cold requires at least one --limit value");
    }
    let positions = collect_cold_unique_positions(config).await?;
    fs::create_dir_all(output_dir)
        .with_context(|| format!("failed to create '{}'", output_dir.display()))?;

    let mut files = Vec::new();
    for &limit in limits {
        let sample = sample_evenly(&positions, limit);
        let path = output_dir.join(sample_file_name(config, limit));
        write_positions(&path, &sample)?;
        files.push(SampleColdFile {
            limit,
            path,
            rows: sample.len(),
        });
    }

    Ok(SampleColdReport {
        unique_cold_positions: positions.len(),
        files,
    })
}

async fn load_warm_position_keys(config: &SandboxConfig, chrom_code: u16) -> Result<HashSet<u64>> {
    let dataset = Dataset::open(config.dataset_path().to_string_lossy().as_ref()).await?;
    let mut scanner = dataset.scan();
    scanner
        .filter("tier = 0")?
        .project(&[config.key.column.as_str()])?
        .materialization_style(MaterializationStyle::AllLate);
    let mut stream = scanner.try_into_stream().await?;
    let mut out = HashSet::new();
    while let Some(batch) = stream.try_next().await? {
        for value in key_values_from_batch(config, &batch)? {
            let position_key = match config.key.mode {
                KeyMode::Position => position_key_from_code(chrom_code, value as i64)?,
                KeyMode::PositionKey => value,
            };
            out.insert(position_key);
        }
    }
    Ok(out)
}

fn load_vcf_workload(path: &Path, chrom: &str) -> Result<Workload> {
    let file = File::open(path).with_context(|| format!("failed to open '{}'", path.display()))?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let chrom_code = chrom_to_code(chrom_bare);
    let mut seen_target = false;
    let mut ordered_positions = Vec::new();
    let mut seen_positions = HashSet::new();
    let mut variant_keys_by_position = HashMap::<u64, Vec<u64>>::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let Some(row_chrom) = fields.next() else {
            continue;
        };
        let row_chrom_bare = row_chrom.strip_prefix("chr").unwrap_or(row_chrom);
        if row_chrom_bare != chrom_bare {
            if seen_target {
                break;
            }
            continue;
        }
        seen_target = true;
        let Some(pos_str) = fields.next() else {
            continue;
        };
        let Some(_id) = fields.next() else {
            continue;
        };
        let Some(ref_allele) = fields.next() else {
            continue;
        };
        let Some(alt_alleles) = fields.next() else {
            continue;
        };
        let pos = pos_str.parse::<i64>()?;
        let end = pos + ref_allele.len() as i64 - 1;

        for probe_start in build_probe_starts(pos, end, ref_allele, alt_alleles, true) {
            let position_key = position_key_from_code(chrom_code, probe_start)?;
            if seen_positions.insert(position_key) {
                ordered_positions.push(position_key);
            }
            let keys = variant_keys_by_position.entry(position_key).or_default();
            for key in warm_variant_key_candidates(position_key, ref_allele, alt_alleles) {
                push_unique(keys, key);
            }
        }
    }

    Ok(Workload {
        ordered_positions,
        variant_keys_by_position,
    })
}

fn warm_variant_key_candidates(position_key: u64, vcf_ref: &str, vcf_alt: &str) -> Vec<u64> {
    let mut keys = Vec::with_capacity(4);
    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (vep_ref, vep_alt) = vcf_to_vep_allele(vcf_ref, alt);
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, &vep_ref, &vep_alt),
        );
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, vcf_ref, alt),
        );
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, vcf_ref, &vep_alt),
        );
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, &vep_ref, alt),
        );
    }
    keys
}

async fn collect_cold_unique_positions(config: &SandboxConfig) -> Result<Vec<u64>> {
    let dataset_path = config.dataset_path();
    let dataset = Dataset::open(dataset_path.to_string_lossy().as_ref())
        .await
        .with_context(|| format!("failed to open Lance dataset '{}'", dataset_path.display()))?;
    let mut scanner = dataset.scan();
    scanner
        .filter("tier = 1")?
        .project(&[config.key.column.as_str()])?
        .materialization_style(MaterializationStyle::AllLate);
    let mut stream = scanner.try_into_stream().await?;

    match config.key.data_type {
        KeyDataType::Uint32 => {
            let mut bitset = PositionBitSet::default();
            while let Some(batch) = stream.try_next().await? {
                for value in key_values_from_batch(config, &batch)? {
                    let value = u32::try_from(value)
                        .with_context(|| format!("cold position {value} does not fit UInt32"))?;
                    bitset.insert(value);
                }
            }
            Ok(bitset.into_positions())
        }
        KeyDataType::Uint64 => {
            let mut positions = Vec::new();
            while let Some(batch) = stream.try_next().await? {
                positions.extend(key_values_from_batch(config, &batch)?);
            }
            positions.sort_unstable();
            positions.dedup();
            Ok(positions)
        }
    }
}

fn write_positions(path: &Path, values: &[u64]) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut file =
        File::create(path).with_context(|| format!("failed to create '{}'", path.display()))?;
    for value in values {
        writeln!(file, "{value}")?;
    }
    Ok(())
}

fn sample_file_name(config: &SandboxConfig, limit: usize) -> String {
    let key_suffix = match config.key.mode {
        KeyMode::Position => "positions_u32",
        KeyMode::PositionKey => "position_keys_u64",
    };
    format!(
        "{}_cold_sample_{}_{}.txt",
        config.dataset.chrom,
        compact_limit(limit),
        key_suffix
    )
}

fn compact_limit(limit: usize) -> String {
    if limit >= 1_000_000 && limit % 1_000_000 == 0 {
        format!("{}m", limit / 1_000_000)
    } else if limit >= 1_000 && limit % 1_000 == 0 {
        format!("{}k", limit / 1_000)
    } else {
        limit.to_string()
    }
}

fn sample_evenly(values: &[u64], limit: usize) -> Vec<u64> {
    if limit == 0 || values.is_empty() {
        return Vec::new();
    }
    if values.len() <= limit {
        return values.to_vec();
    }
    if limit == 1 {
        return vec![values[0]];
    }

    let last_input = values.len() - 1;
    let last_output = limit - 1;
    (0..limit)
        .map(|idx| values[(idx * last_input) / last_output])
        .collect()
}

fn build_probe_starts(
    norm_start_i64: i64,
    norm_end_i64: i64,
    vcf_ref: &str,
    vcf_alt: &str,
    extended_probes: bool,
) -> Vec<i64> {
    let mut probe_starts = Vec::with_capacity(6);
    push_unique(&mut probe_starts, norm_start_i64);

    if !extended_probes {
        return probe_starts;
    }

    if norm_start_i64 == norm_end_i64 {
        push_unique(&mut probe_starts, norm_start_i64.saturating_add(1));
    } else {
        push_unique(&mut probe_starts, norm_end_i64);
    }

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (_, _, input_start) = vcf_to_vep_input_allele(norm_start_i64, vcf_ref, alt);
        push_unique(&mut probe_starts, input_start);

        let shift_usize = common_prefix_len(vcf_ref, alt);
        if shift_usize == 0 {
            continue;
        }
        if let Some(shifted_start) = norm_start_i64.checked_add(shift_usize as i64) {
            push_unique(&mut probe_starts, shifted_start);
        }
    }

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (ref_event_len, alt_event_len) = canonical_event_lengths(vcf_ref, alt);
        if ref_event_len == 0 || alt_event_len != 0 {
            continue;
        }
        let del_len = ref_event_len as i64;
        let max_shift = del_len.min(32);
        for base_start in [norm_start_i64, norm_start_i64.saturating_sub(1)] {
            for shift in 0..=max_shift {
                let Some(candidate_start) = base_start.checked_add(shift) else {
                    continue;
                };
                let Some(candidate_end) = candidate_start.checked_add(del_len - 1) else {
                    continue;
                };
                if candidate_start > norm_end_i64 || candidate_end < norm_start_i64 {
                    continue;
                }
                push_unique(&mut probe_starts, candidate_start);
            }
        }
    }

    probe_starts
}

#[inline]
fn common_prefix_len(left: &str, right: &str) -> usize {
    left.as_bytes()
        .iter()
        .zip(right.as_bytes().iter())
        .take_while(|(a, b)| a == b)
        .count()
}

#[inline]
fn canonical_event_lengths(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();
    let mut ref_start = 0usize;
    let mut alt_start = 0usize;
    while ref_start < ref_bytes.len()
        && alt_start < alt_bytes.len()
        && ref_bytes[ref_start] == alt_bytes[alt_start]
    {
        ref_start += 1;
        alt_start += 1;
    }
    let mut ref_end = ref_bytes.len();
    let mut alt_end = alt_bytes.len();
    while ref_end > ref_start
        && alt_end > alt_start
        && ref_bytes[ref_end - 1] == alt_bytes[alt_end - 1]
    {
        ref_end -= 1;
        alt_end -= 1;
    }
    (ref_end - ref_start, alt_end - alt_start)
}

fn push_unique<T: Copy + PartialEq>(values: &mut Vec<T>, value: T) {
    if !values.contains(&value) {
        values.push(value);
    }
}

fn position_from_key(position_key: u64) -> u64 {
    position_key & ((1_u64 << 48) - 1)
}

#[derive(Default)]
struct PositionBitSet {
    bits: Vec<u64>,
    len: usize,
}

impl PositionBitSet {
    fn insert(&mut self, value: u32) {
        let bit = value as usize;
        let word = bit / 64;
        let mask = 1_u64 << (bit % 64);
        if word >= self.bits.len() {
            self.bits.resize(word + 1, 0);
        }
        if self.bits[word] & mask == 0 {
            self.bits[word] |= mask;
            self.len += 1;
        }
    }

    fn into_positions(self) -> Vec<u64> {
        let mut out = Vec::with_capacity(self.len);
        for (word_idx, word) in self.bits.into_iter().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let bit = remaining.trailing_zeros() as usize;
                out.push((word_idx * 64 + bit) as u64);
                remaining &= remaining - 1;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extended_probe_starts_include_shifted_insertion_positions() {
        let mut probe_starts = build_probe_starts(165_387_539, 165_387_541, "CTG", "CTCTGTG", true);
        probe_starts.sort_unstable();

        assert_eq!(probe_starts, vec![165_387_539, 165_387_540, 165_387_541]);
    }

    #[test]
    fn even_sample_spreads_across_sorted_positions() {
        let sample = sample_evenly(&[10, 20, 30, 40, 50], 3);

        assert_eq!(sample, vec![10, 30, 50]);
    }
}
