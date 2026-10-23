use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::*;
use datafusion::arrow::datatypes::*;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion_bio_format_ensembl_cache::{
    CacheSourceType, EnsemblCacheOptions, EnsemblCacheTableProvider, EnsemblEntityKind,
};
use datafusion_bio_function_vep::cache::af_bundle::af_column_order;
use datafusion_bio_function_vep::parquet_cache::write::{
    VariationParquetShardWriter, encode_variation_batch, variation_output_schema,
};
use rayon::prelude::*;

const IO_BUFFER_SIZE: usize = 64 * 1024;
const WARM_AF_THRESHOLD: f64 = 0.01;
const WARM_POSITION_RADIUS: i64 = 1;

#[derive(Debug, Clone)]
struct BgzfPartition {
    start_compressed: u64,
    end_compressed: u64,
    skip_first_line: bool,
}

fn scan_bgzf_block_offsets(file_path: &Path, file_size: u64) -> Result<Vec<u64>, Box<dyn std::error::Error>> {
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

    let mut partitions = Vec::with_capacity(boundaries.len());
    for i in 0..boundaries.len() {
        let start = boundaries[i];
        let end = if i + 1 < boundaries.len() {
            boundaries[i + 1]
        } else {
            file_size
        };

        let skip_first_line = if start > 0 {
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

fn collect_warm_starts(path: &Path) -> Result<HashSet<i64>, Box<dyn std::error::Error>> {
    let num_threads = 10;
    let partitions = compute_bgzf_partitions(path, num_threads)?;

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
                    let _ = reader.read_line(&mut discard);
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
                            if max_af >= WARM_AF_THRESHOLD {
                                high_af_starts.push(start_val);
                            }
                        }
                    }
                    Err(_) => break,
                }
            }
            (all_starts, high_af_starts)
        })
        .collect();

    let mut source_positions: HashSet<i64> = HashSet::new();
    let mut high_af: Vec<i64> = Vec::new();

    for (starts, highs) in results {
        source_positions.extend(starts);
        high_af.extend(highs);
    }

    let mut warm_set: HashSet<i64> = HashSet::new();
    for start in high_af {
        for offset in -WARM_POSITION_RADIUS..=WARM_POSITION_RADIUS {
            let cand = start + offset;
            if source_positions.contains(&cand) {
                warm_set.insert(cand);
            }
        }
    }

    Ok(warm_set)
}

#[inline]
fn parse_max_af_field(s: &str) -> f64 {
    if s.is_empty() || s == "." {
        return 0.0;
    }
    let mut max_val = 0.0f64;
    for part in s.split(',') {
        let val_str = if let Some((_, v)) = part.rsplit_once(':') {
            v
        } else {
            part
        };
        if let Ok(v) = val_str.trim().parse::<f64>() {
            if v.is_finite() && v > max_val {
                max_val = v;
            }
        }
    }
    max_val
}

struct SourceIds {
    clinvar: String,
    cosmic: String,
    dbsnp: String,
}

impl SourceIds {
    fn new() -> Self {
        Self {
            clinvar: String::with_capacity(32),
            cosmic: String::with_capacity(32),
            dbsnp: String::with_capacity(32),
        }
    }

    fn clear(&mut self) {
        self.clinvar.clear();
        self.cosmic.clear();
        self.dbsnp.clear();
    }

    fn populate(&mut self, var_name: Option<&str>, var_synonyms: Option<&str>) {
        if let Some(name) = var_name {
            self.assign_by_pattern(name);
        }
        if let Some(synonyms) = var_synonyms {
            self.parse_var_synonyms(synonyms);
        }
    }

    fn assign_by_pattern(&mut self, value: &str) {
        let value = value.trim();
        if value.is_empty() || value == "." {
            return;
        }
        if value.starts_with("rs") && value[2..].chars().all(|c| c.is_ascii_digit()) {
            Self::push_value(&mut self.dbsnp, value);
        } else if value.starts_with("COSM") || value.starts_with("COSV") {
            Self::push_value(&mut self.cosmic, value);
        } else if value.starts_with("RCV") || value.starts_with("VCV") {
            Self::push_value(&mut self.clinvar, value);
        }
    }

    fn parse_var_synonyms(&mut self, synonyms: &str) {
        let mut current_slot: Option<u8> = None;
        for token in synonyms.split(',') {
            let token = token.trim();
            if token.is_empty() {
                continue;
            }

            if let Some((source, value)) = token.split_once("::") {
                current_slot = Self::lookup_source(source);
                if let Some(slot) = current_slot {
                    self.push_slot(slot, value.trim());
                } else {
                    self.assign_by_pattern(token);
                }
                continue;
            }

            if let Some((source, value)) = token.split_once(':') {
                let slot = Self::lookup_source(source);
                if let Some(s) = slot {
                    current_slot = Some(s);
                    self.push_slot(s, value.trim());
                    continue;
                }
            }

            if let Some(slot) = current_slot {
                self.push_slot(slot, token);
            } else {
                self.assign_by_pattern(token);
            }
        }
    }

    fn lookup_source(source: &str) -> Option<u8> {
        let s = source.to_ascii_lowercase();
        if s.contains("dbsnp") {
            Some(0)
        } else if s.contains("cosmic") {
            Some(1)
        } else if s.contains("clinvar") {
            Some(2)
        } else {
            None
        }
    }

    fn push_slot(&mut self, slot: u8, val: &str) {
        match slot {
            0 => Self::push_value(&mut self.dbsnp, val),
            1 => Self::push_value(&mut self.cosmic, val),
            2 => Self::push_value(&mut self.clinvar, val),
            _ => {}
        }
    }

    fn push_value(buf: &mut String, val: &str) {
        let val = val.trim();
        if val.is_empty() || val == "." {
            return;
        }
        for existing in buf.split(',') {
            if existing == val {
                return;
            }
        }
        if !buf.is_empty() {
            buf.push(',');
        }
        buf.push_str(val);
    }
}

fn intermediate_variation_schema(has_clin_sig_ref_allele: bool) -> SchemaRef {
    let mut fields = vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new("allele_string", DataType::Utf8, false),
        Field::new("failed", DataType::Int8, true),
        Field::new("variation_name", DataType::Utf8, false),
        Field::new("clin_sig", DataType::Utf8, true),
        Field::new("clin_sig_allele", DataType::Utf8, true),
    ];
    if has_clin_sig_ref_allele {
        fields.push(Field::new("clin_sig_ref_allele", DataType::Utf8, true));
    }
    fields.extend(vec![
        Field::new("clinical_impact", DataType::Utf8, true),
        Field::new("phenotype_or_disease", DataType::Int8, true),
        Field::new("pubmed", DataType::Utf8, true),
        Field::new("somatic", DataType::Int8, true),
        Field::new("minor_allele", DataType::Utf8, true),
        Field::new("minor_allele_freq", DataType::Float64, true),
    ]);
    for name in af_column_order() {
        fields.push(Field::new(name, DataType::Utf8, true));
    }
    fields.extend(vec![
        Field::new("clinvar_ids", DataType::Utf8, true),
        Field::new("cosmic_ids", DataType::Utf8, true),
        Field::new("dbsnp_ids", DataType::Utf8, true),
        Field::new("tier", DataType::Int8, false),
    ]);
    Arc::new(Schema::new(fields))
}

struct VariationBatchBuilder {
    chrom: StringBuilder,
    start: UInt32Builder,
    end: UInt32Builder,
    allele_string: StringBuilder,
    failed: Int8Builder,
    variation_name: StringBuilder,
    clin_sig: StringBuilder,
    clin_sig_allele: StringBuilder,
    clin_sig_ref_allele: StringBuilder,
    clinical_impact: StringBuilder,
    phenotype_or_disease: Int8Builder,
    pubmed: StringBuilder,
    somatic: Int8Builder,
    minor_allele: StringBuilder,
    minor_allele_freq: Float64Builder,
    af_builders: Vec<StringBuilder>,
    clinvar_ids: StringBuilder,
    cosmic_ids: StringBuilder,
    dbsnp_ids: StringBuilder,
    tier: Int8Builder,
    pub num_rows: usize,
}

impl VariationBatchBuilder {
    fn new(capacity: usize) -> Self {
        let mut af_builders = Vec::with_capacity(27);
        for _ in 0..27 {
            af_builders.push(StringBuilder::with_capacity(capacity, capacity * 2));
        }
        Self {
            chrom: StringBuilder::with_capacity(capacity, capacity * 2),
            start: UInt32Builder::with_capacity(capacity),
            end: UInt32Builder::with_capacity(capacity),
            allele_string: StringBuilder::with_capacity(capacity, capacity * 4),
            failed: Int8Builder::with_capacity(capacity),
            variation_name: StringBuilder::with_capacity(capacity, capacity * 12),
            clin_sig: StringBuilder::with_capacity(capacity, capacity * 2),
            clin_sig_allele: StringBuilder::with_capacity(capacity, capacity * 2),
            clin_sig_ref_allele: StringBuilder::with_capacity(capacity, capacity * 2),
            clinical_impact: StringBuilder::with_capacity(capacity, capacity * 2),
            phenotype_or_disease: Int8Builder::with_capacity(capacity),
            pubmed: StringBuilder::with_capacity(capacity, capacity * 4),
            somatic: Int8Builder::with_capacity(capacity),
            minor_allele: StringBuilder::with_capacity(capacity, capacity * 2),
            minor_allele_freq: Float64Builder::with_capacity(capacity),
            af_builders,
            clinvar_ids: StringBuilder::with_capacity(capacity, capacity * 4),
            cosmic_ids: StringBuilder::with_capacity(capacity, capacity * 4),
            dbsnp_ids: StringBuilder::with_capacity(capacity, capacity * 12),
            tier: Int8Builder::with_capacity(capacity),
            num_rows: 0,
        }
    }

    #[inline]
    fn add_row(
        &mut self,
        chrom: &str,
        fields: &[&str],
        tier_val: i8,
        ids: &SourceIds,
        has_clin_sig_ref: bool,
    ) {
        self.chrom.append_value(chrom);

        let start_val: u32 = fields.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
        let end_val: u32 = fields.get(5)
            .and_then(|s| if *s == "." || s.is_empty() { None } else { s.parse().ok() })
            .unwrap_or(start_val);
        self.start.append_value(start_val);
        self.end.append_value(end_val);

        self.allele_string.append_value(fields.get(6).copied().unwrap_or("."));

        // failed
        if fields.get(2) == Some(&"1") {
            self.failed.append_value(1);
        } else {
            self.failed.append_null();
        }

        // variation_name
        self.variation_name.append_value(fields.get(1).copied().unwrap_or("."));

        // clin_sig
        append_opt_str(&mut self.clin_sig, fields.get(8).copied());
        // clin_sig_allele
        append_opt_str(&mut self.clin_sig_allele, fields.get(10).copied());
        if has_clin_sig_ref {
            append_opt_str(&mut self.clin_sig_ref_allele, fields.get(11).copied());
        }
        // clinical_impact
        append_opt_str(&mut self.clinical_impact, fields.get(12).copied());

        // phenotype_or_disease
        if fields.get(9) == Some(&"1") {
            self.phenotype_or_disease.append_value(1);
        } else {
            self.phenotype_or_disease.append_null();
        }

        // pubmed
        append_opt_str(&mut self.pubmed, fields.get(13).copied());

        // somatic
        if fields.get(3) == Some(&"1") {
            self.somatic.append_value(1);
        } else {
            self.somatic.append_null();
        }

        // minor_allele, minor_allele_freq
        self.minor_allele.append_null();
        self.minor_allele_freq.append_null();

        // 27 AF fields:
        // Canonical order:
        // af_global (6): AF(15), AFR(16), AMR(17), EAS(18), EUR(19), SAS(20)
        // af_gnomade (10): gnomADe(21), gnomADe_AFR(22), gnomADe_AMR(23), gnomADe_ASJ(24), gnomADe_EAS(25), gnomADe_FIN(26), gnomADe_NFE(27), gnomADe_SAS(29), gnomADe_MID(30), gnomADe_REMAINING(28)
        // af_gnomadg (11): gnomADg(31), gnomADg_AFR(32), gnomADg_AMI(33), gnomADg_AMR(34), gnomADg_ASJ(35), gnomADg_EAS(36), gnomADg_FIN(37), gnomADg_MID(38), gnomADg_NFE(39), gnomADg_SAS(41), gnomADg_REMAINING(40)
        const AF_INDICES: [usize; 27] = [
            15, 16, 17, 18, 19, 20,
            21, 22, 23, 24, 25, 26, 27, 29, 30, 28,
            31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 40,
        ];
        for (i, &tab_idx) in AF_INDICES.iter().enumerate() {
            append_opt_str(&mut self.af_builders[i], fields.get(tab_idx).copied());
        }

        // Source IDs
        if ids.clinvar.is_empty() {
            self.clinvar_ids.append_null();
        } else {
            self.clinvar_ids.append_value(&ids.clinvar);
        }
        if ids.cosmic.is_empty() {
            self.cosmic_ids.append_null();
        } else {
            self.cosmic_ids.append_value(&ids.cosmic);
        }
        if ids.dbsnp.is_empty() {
            self.dbsnp_ids.append_null();
        } else {
            self.dbsnp_ids.append_value(&ids.dbsnp);
        }

        self.tier.append_value(tier_val);
        self.num_rows += 1;
    }

    fn finish(mut self, schema: &SchemaRef, has_clin_sig_ref: bool) -> RecordBatch {
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(self.chrom.finish()),
            Arc::new(self.start.finish()),
            Arc::new(self.end.finish()),
            Arc::new(self.allele_string.finish()),
            Arc::new(self.failed.finish()),
            Arc::new(self.variation_name.finish()),
            Arc::new(self.clin_sig.finish()),
            Arc::new(self.clin_sig_allele.finish()),
        ];
        if has_clin_sig_ref {
            columns.push(Arc::new(self.clin_sig_ref_allele.finish()));
        }
        columns.extend(vec![
            Arc::new(self.clinical_impact.finish()) as ArrayRef,
            Arc::new(self.phenotype_or_disease.finish()),
            Arc::new(self.pubmed.finish()),
            Arc::new(self.somatic.finish()),
            Arc::new(self.minor_allele.finish()),
            Arc::new(self.minor_allele_freq.finish()),
        ]);
        for mut b in self.af_builders {
            columns.push(Arc::new(b.finish()));
        }
        columns.extend(vec![
            Arc::new(self.clinvar_ids.finish()) as ArrayRef,
            Arc::new(self.cosmic_ids.finish()),
            Arc::new(self.dbsnp_ids.finish()),
            Arc::new(self.tier.finish()),
        ]);
        RecordBatch::try_new(Arc::clone(schema), columns).expect("valid batch")
    }
}

#[inline]
fn append_opt_str(builder: &mut StringBuilder, val: Option<&str>) {
    match val {
        Some(s) if !s.is_empty() && s != "." => builder.append_value(s),
        _ => builder.append_null(),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t0 = Instant::now();
    let cache_root = "/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/cache_bench_dev/prepared_cache";
    let all_vars = Path::new(cache_root).join("20").join("all_vars.gz");

    println!("Step 1: Collecting warm starts natively...");
    let t_warm = Instant::now();
    let warm_starts = collect_warm_starts(&all_vars)?;
    println!("Warm starts collected: {} positions in {:.2?}", warm_starts.len(), t_warm.elapsed());

    println!("Step 2: Preparing schema...");
    let provider_opts = EnsemblCacheOptions::new(cache_root)
        .with_cache_source_type(CacheSourceType::Merged)
        .with_expected_cache_version("116");
    let provider = EnsemblCacheTableProvider::for_entity(EnsemblEntityKind::Variation, provider_opts)?;
    let source_schema = provider.schema();
    let out_schema = variation_output_schema(source_schema.as_ref(), "merged")?;
    let intermediate_schema = intermediate_variation_schema(true);

    println!("Step 3: Streaming Partitions to Parquet...");
    let out_file = Path::new("/tmp/chr20_native_test.parquet");
    let _ = std::fs::remove_file(out_file);
    let mut writer = VariationParquetShardWriter::create(out_file, Arc::clone(&out_schema))?;

    let partitions = compute_bgzf_partitions(&all_vars, 10)?;
    println!("Computed {} partitions", partitions.len());

    // Pass 1: Warm tier (parallel scan, small memory)
    let t_pass1 = Instant::now();
    let warm_batches: Vec<RecordBatch> = partitions
        .par_iter()
        .map(|part| {
            let file = File::open(&all_vars).expect("open file");
            let mut reader = noodles_bgzf::io::Reader::new(BufReader::with_capacity(IO_BUFFER_SIZE, file));
            if part.start_compressed > 0 {
                let vpos = noodles_bgzf::VirtualPosition::try_from((part.start_compressed, 0)).expect("vpos");
                reader.seek(vpos).expect("seek");
                if part.skip_first_line {
                    let mut discard = String::new();
                    let _ = reader.read_line(&mut discard);
                }
            }

            let mut builder = VariationBatchBuilder::new(50_000);
            let mut ids = SourceIds::new();
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
                        let mut fields = [""; 48];
                        let mut count = 0;
                        for f in trimmed.split('\t') {
                            if count < 48 {
                                fields[count] = f;
                                count += 1;
                            }
                        }
                        let start_val: i64 = fields[4].parse().unwrap_or(0);
                        if warm_starts.contains(&start_val) {
                            ids.clear();
                            ids.populate(fields.get(1).copied(), fields.get(14).copied());
                            builder.add_row("20", &fields[..count], 0, &ids, true);
                        }
                    }
                    Err(_) => break,
                }
            }
            if builder.num_rows > 0 {
                let batch = builder.finish(&intermediate_schema, true);
                Some(encode_variation_batch(&batch, &out_schema).expect("encode"))
            } else {
                None
            }
        })
        .flatten()
        .collect();

    let mut warm_rows = 0;
    for b in warm_batches {
        warm_rows += b.num_rows();
        writer.write(&b)?;
    }
    println!("Pass 1 (Warm) wrote {} rows in {:.2?}", warm_rows, t_pass1.elapsed());

    // Pass 2: Cold tier (pipelined parallel Rayon decompression/encoding + streaming Parquet write)
    let t_pass2 = Instant::now();
    let (senders, receivers): (Vec<_>, Vec<_>) = (0..partitions.len())
        .map(|_| std::sync::mpsc::channel::<RecordBatch>())
        .unzip();

    let all_vars_clone = all_vars.clone();
    let warm_starts_clone = warm_starts.clone();
    let intermediate_schema_clone = Arc::clone(&intermediate_schema);
    let out_schema_clone = Arc::clone(&out_schema);
    let partitions_clone = partitions.clone();

    let rayon_handle = std::thread::spawn(move || {
        partitions_clone
            .into_par_iter()
            .zip(senders.into_par_iter())
            .for_each(|(part, tx)| {
                let file = File::open(&all_vars_clone).expect("open file");
                let mut reader = noodles_bgzf::io::Reader::new(BufReader::with_capacity(IO_BUFFER_SIZE, file));
                if part.start_compressed > 0 {
                    let vpos = noodles_bgzf::VirtualPosition::try_from((part.start_compressed, 0)).expect("vpos");
                    reader.seek(vpos).expect("seek");
                    if part.skip_first_line {
                        let mut discard = String::new();
                        let _ = reader.read_line(&mut discard);
                    }
                }

                let mut builder = VariationBatchBuilder::new(100_000);
                let mut ids = SourceIds::new();
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
                            let mut fields = [""; 48];
                            let mut count = 0;
                            for f in trimmed.split('\t') {
                                if count < 48 {
                                    fields[count] = f;
                                    count += 1;
                                }
                            }
                            let start_val: i64 = fields[4].parse().unwrap_or(0);
                            if !warm_starts_clone.contains(&start_val) {
                                ids.clear();
                                ids.populate(fields.get(1).copied(), fields.get(14).copied());
                                builder.add_row("20", &fields[..count], 1, &ids, true);

                                if builder.num_rows >= 100_000 {
                                    let batch = builder.finish(&intermediate_schema_clone, true);
                                    let encoded = encode_variation_batch(&batch, &out_schema_clone).expect("encode");
                                    let _ = tx.send(encoded);
                                    builder = VariationBatchBuilder::new(100_000);
                                }
                            }
                        }
                        Err(_) => break,
                    }
                }
                if builder.num_rows > 0 {
                    let batch = builder.finish(&intermediate_schema_clone, true);
                    let encoded = encode_variation_batch(&batch, &out_schema_clone).expect("encode");
                    let _ = tx.send(encoded);
                }
            });
    });

    let mut cold_rows = 0;
    for rx in receivers {
        while let Ok(b) = rx.recv() {
            cold_rows += b.num_rows();
            writer.write(&b)?;
        }
    }
    rayon_handle.join().unwrap();
    println!("Pass 2 (Cold) pipelined encoded + wrote {} rows in {:.2?}", cold_rows, t_pass2.elapsed());

    let total_written = writer.finish()?;
    println!("SUCCESS! Total rows written to Parquet: {} in {:.2?}", total_written, t0.elapsed());
    Ok(())
}
