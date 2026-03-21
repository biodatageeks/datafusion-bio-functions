use std::collections::{HashMap, HashSet};

use datafusion::arrow::array::RecordBatch;
use datafusion::common::{DataFusionError, Result};

use crate::annotate_provider::string_at;

const COALESCE_GAP: i64 = 1_000_000;

/// Maximum number of interval clauses before falling back to chrom-only filtering
/// to avoid massive OR clauses that overflow DataFusion's parser/planner stack.
const MAX_INTERVAL_CLAUSES: usize = 50;

#[derive(Debug)]
pub struct MissWorklist {
    /// Bare chromosome names (without "chr" prefix), e.g. "1", "X", "MT".
    pub chroms: HashSet<String>,
    pub intervals: HashMap<String, Vec<(i64, i64)>>,
}

impl MissWorklist {
    /// Build a single-chrom worklist without scanning base_batches.
    ///
    /// Used by the partitioned path where the per-chrom parquet file IS the
    /// filter, so no variant scanning is needed. The `expanded_chroms()` set
    /// will cover both bare and "chr"-prefixed forms.
    pub fn for_chrom(chrom: &str) -> Self {
        let bare = chrom.strip_prefix("chr").unwrap_or(chrom).to_string();
        let mut chroms = HashSet::new();
        chroms.insert(bare);
        MissWorklist {
            chroms,
            intervals: HashMap::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.chroms.is_empty()
    }

    pub fn expanded_chroms(&self) -> HashSet<String> {
        let mut expanded = HashSet::new();
        for c in &self.chroms {
            let escaped = c.replace('\'', "''");
            expanded.insert(escaped.clone());
            if let Some(bare) = escaped.strip_prefix("chr") {
                expanded.insert(bare.to_string());
            } else {
                expanded.insert(format!("chr{escaped}"));
            }
        }
        expanded
    }

    pub fn chrom_filter_clause(&self) -> String {
        if self.chroms.is_empty() {
            return String::new();
        }
        let literals: Vec<String> = self
            .expanded_chroms()
            .iter()
            .map(|c| format!("'{c}'"))
            .collect();
        format!(" WHERE chrom IN ({})", literals.join(", "))
    }

    pub fn interval_filter_sql(&self) -> String {
        if self.intervals.is_empty() {
            return self.chrom_filter_clause();
        }
        // If too many intervals, fall back to chrom-only filter to avoid
        // massive OR clauses that overflow DataFusion's parser/planner stack.
        let total_intervals: usize = self.intervals.values().map(|v| v.len()).sum();
        if total_intervals > MAX_INTERVAL_CLAUSES {
            return self.chrom_filter_clause();
        }
        let mut clauses: Vec<String> = Vec::new();
        for (chrom, intervals) in &self.intervals {
            let escaped = chrom.replace('\'', "''");
            let bare = escaped.strip_prefix("chr").unwrap_or(&escaped);
            let prefixed = if escaped.starts_with("chr") {
                escaped.clone()
            } else {
                format!("chr{escaped}")
            };
            let chrom_variants = [bare.to_string(), prefixed];

            for cv in &chrom_variants {
                for (lo, hi) in intervals {
                    clauses.push(format!(
                        "(chrom = '{cv}' AND start <= {hi} AND \"end\" >= {lo})"
                    ));
                }
            }
        }
        if clauses.is_empty() {
            return self.chrom_filter_clause();
        }
        format!(" WHERE {}", clauses.join(" OR "))
    }
}

pub fn collect_miss_worklist(
    batches: &[RecordBatch],
    cache_hit_count: &mut usize,
    cache_miss_count: &mut usize,
) -> Result<MissWorklist> {
    let mut chroms = HashSet::new();
    let mut raw_intervals: HashMap<String, Vec<(i64, i64)>> = HashMap::new();

    for batch in batches {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").ok();
        let start_idx = schema.index_of("start").ok();
        let end_idx = schema.index_of("end").ok();
        let msc_idx = schema.index_of("cache_most_severe_consequence").ok();

        for row in 0..batch.num_rows() {
            let is_miss = match msc_idx {
                Some(idx) => string_at(batch.column(idx).as_ref(), row).is_none(),
                None => true,
            };

            if is_miss {
                *cache_miss_count += 1;
                if let Some(ci) = chrom_idx {
                    if let Some(c) = string_at(batch.column(ci).as_ref(), row) {
                        let norm = c.strip_prefix("chr").unwrap_or(&c).to_string();
                        chroms.insert(norm.clone());

                        if let (Some(si), Some(ei)) = (start_idx, end_idx) {
                            let s =
                                crate::annotate_provider::int64_at(batch.column(si).as_ref(), row);
                            let e =
                                crate::annotate_provider::int64_at(batch.column(ei).as_ref(), row);
                            if let (Some(s), Some(e)) = (s, e) {
                                raw_intervals.entry(norm).or_default().push((s, e));
                            }
                        }
                    }
                }
            } else {
                *cache_hit_count += 1;
            }
        }
    }

    let mut intervals = HashMap::new();
    for (chrom, mut ivs) in raw_intervals {
        ivs.sort_unstable_by_key(|&(s, _)| s);
        let mut merged: Vec<(i64, i64)> = Vec::with_capacity(ivs.len());
        for (s, e) in ivs {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 + COALESCE_GAP {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        intervals.insert(chrom, merged);
    }

    Ok(MissWorklist { chroms, intervals })
}
