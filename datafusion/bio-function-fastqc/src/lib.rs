//! Streaming FastQC quality-control modules for Apache DataFusion.
//!
//! Each module implements `QcModule` (update/merge/finalize), mirroring
//! RastQC's process_sequence/merge_from/calculate_results.

pub mod adapter_content;
pub mod basic_stats;
pub mod dup_levels;
pub mod overrepresented;
pub mod per_base_content;
pub mod per_base_n;
pub mod per_base_quality;
pub mod per_seq_gc;
pub mod per_seq_quality;
pub mod physical_exec;
pub mod seq_length;

pub use physical_exec::FastqcExec;

use std::any::Any;

/// One row of the uniform tidy output schema.
#[derive(Debug, Clone, PartialEq)]
pub struct TidyRow {
    pub module: &'static str,
    pub label: Option<String>,
    pub position: Option<i32>,
    pub metric: String,
    pub value: Option<f64>,
    pub value_str: Option<String>,
}

impl TidyRow {
    pub fn num(module: &'static str, metric: &str, value: f64) -> Self {
        TidyRow {
            module,
            label: None,
            position: None,
            metric: metric.to_string(),
            value: Some(value),
            value_str: None,
        }
    }
    pub fn status(module: &'static str, status: &str) -> Self {
        TidyRow {
            module,
            label: None,
            position: None,
            metric: "status".to_string(),
            value: None,
            value_str: Some(status.to_string()),
        }
    }
}

/// A streaming QC module accumulator.
pub trait QcModule: Send {
    /// Stable identifier, e.g. "basic_stats".
    fn name(&self) -> &'static str;
    /// Fold one read (raw sequence bytes and raw phred+33 quality bytes).
    fn update(&mut self, seq: &[u8], qual: &[u8]);
    /// Merge another accumulator of the SAME concrete type into self.
    fn merge(&mut self, other: &dyn QcModule);
    /// Emit tidy rows for this module's result.
    fn finalize(&self, out: &mut Vec<TidyRow>);
    /// For downcasting in `merge`.
    fn as_any(&self) -> &dyn Any;
}

use std::sync::Arc;

use datafusion::arrow::array::builder::{Float64Builder, Int32Builder, StringBuilder};
use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::common::{DataFusionError, Result};

use adapter_content::AdapterContent;
use basic_stats::BasicStats;
use dup_levels::DuplicationLevels;
use overrepresented::OverrepresentedSeqs;
use per_base_content::PerBaseContent;
use per_base_n::PerBaseN;
use per_base_quality::PerBaseQuality;
use per_seq_gc::PerSeqGc;
use per_seq_quality::PerSeqQuality;
use seq_length::SeqLength;

pub const ALL_MODULES: [&str; 10] = [
    "basic_stats",
    "per_base_quality",
    "per_seq_quality",
    "per_base_content",
    "per_seq_gc",
    "per_base_n",
    "seq_length",
    "overrepresented",
    "adapter_content",
    "dup_levels",
];

/// The fixed tidy output schema (independent of module selection).
pub fn tidy_schema() -> SchemaRef {
    Arc::new(Schema::new(vec![
        Field::new("module", DataType::Utf8, false),
        Field::new("label", DataType::Utf8, true),
        Field::new("position", DataType::Int32, true),
        Field::new("metric", DataType::Utf8, false),
        Field::new("value", DataType::Float64, true),
        Field::new("value_str", DataType::Utf8, true),
    ]))
}

fn make_module(name: &str) -> Result<Box<dyn QcModule>> {
    Ok(match name {
        "basic_stats" => Box::new(BasicStats::new()),
        "per_base_quality" => Box::new(PerBaseQuality::new()),
        "per_seq_quality" => Box::new(PerSeqQuality::new()),
        "per_base_content" => Box::new(PerBaseContent::new()),
        "per_seq_gc" => Box::new(PerSeqGc::new()),
        "per_base_n" => Box::new(PerBaseN::new()),
        "seq_length" => Box::new(SeqLength::new()),
        "overrepresented" => Box::new(OverrepresentedSeqs::new()),
        "adapter_content" => Box::new(AdapterContent::new()),
        "dup_levels" => Box::new(DuplicationLevels::new()),
        other => {
            return Err(DataFusionError::Plan(format!(
                "unknown fastqc module '{other}'; valid: {}",
                ALL_MODULES.join(", ")
            )))
        },
    })
}

pub struct ModuleSet {
    modules: Vec<Box<dyn QcModule>>,
}

impl ModuleSet {
    /// Build the requested modules (None => all), preserving ALL_MODULES order.
    pub fn build(selection: Option<&[String]>) -> Result<Self> {
        let names: Vec<&str> = match selection {
            None => ALL_MODULES.to_vec(),
            Some(sel) => {
                for s in sel {
                    if !ALL_MODULES.contains(&s.as_str()) {
                        return Err(DataFusionError::Plan(format!(
                            "unknown fastqc module '{s}'; valid: {}",
                            ALL_MODULES.join(", ")
                        )));
                    }
                }
                ALL_MODULES
                    .iter()
                    .copied()
                    .filter(|m| sel.iter().any(|s| s == m))
                    .collect()
            },
        };
        let modules = names
            .iter()
            .map(|n| make_module(n))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self { modules })
    }

    pub fn update_batch(&mut self, batch: &RecordBatch) -> Result<()> {
        let seq = batch
            .column_by_name("sequence")
            .ok_or_else(|| DataFusionError::Execution("fastqc input missing 'sequence'".into()))?
            .as_string::<i32>();
        let qual = batch
            .column_by_name("quality_scores")
            .ok_or_else(|| {
                DataFusionError::Execution("fastqc input missing 'quality_scores'".into())
            })?
            .as_string::<i32>();
        for i in 0..batch.num_rows() {
            if seq.is_null(i) {
                continue;
            }
            let s = seq.value(i).as_bytes();
            let q = if qual.is_null(i) {
                &[][..]
            } else {
                qual.value(i).as_bytes()
            };
            for m in self.modules.iter_mut() {
                m.update(s, q);
            }
        }
        Ok(())
    }

    /// Merge another ModuleSet (same selection/order) into self.
    pub fn merge(&mut self, other: ModuleSet) {
        for (a, b) in self.modules.iter_mut().zip(other.modules.iter()) {
            a.merge(b.as_ref());
        }
    }

    pub fn finalize(self) -> Result<RecordBatch> {
        let mut rows: Vec<TidyRow> = Vec::new();
        for m in &self.modules {
            m.finalize(&mut rows);
        }
        let mut module_b = StringBuilder::new();
        let mut label_b = StringBuilder::new();
        let mut pos_b = Int32Builder::new();
        let mut metric_b = StringBuilder::new();
        let mut value_b = Float64Builder::new();
        let mut vstr_b = StringBuilder::new();
        for r in rows {
            module_b.append_value(r.module);
            match r.label {
                Some(l) => label_b.append_value(l),
                None => label_b.append_null(),
            }
            match r.position {
                Some(p) => pos_b.append_value(p),
                None => pos_b.append_null(),
            }
            metric_b.append_value(r.metric);
            match r.value {
                Some(v) => value_b.append_value(v),
                None => value_b.append_null(),
            }
            match r.value_str {
                Some(s) => vstr_b.append_value(s),
                None => vstr_b.append_null(),
            }
        }
        RecordBatch::try_new(
            tidy_schema(),
            vec![
                Arc::new(module_b.finish()),
                Arc::new(label_b.finish()),
                Arc::new(pos_b.finish()),
                Arc::new(metric_b.finish()),
                Arc::new(value_b.finish()),
                Arc::new(vstr_b.finish()),
            ],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

#[cfg(test)]
mod set_tests {
    use datafusion::arrow::array::StringArray;

    use super::*;

    fn batch(seqs: &[&str], quals: &[&str]) -> RecordBatch {
        RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("sequence", DataType::Utf8, true),
                Field::new("quality_scores", DataType::Utf8, true),
            ])),
            vec![
                Arc::new(StringArray::from(seqs.to_vec())),
                Arc::new(StringArray::from(quals.to_vec())),
            ],
        )
        .unwrap()
    }

    #[test]
    fn selection_validates_and_orders() {
        assert!(ModuleSet::build(Some(&["bogus".to_string()])).is_err());
        let set =
            ModuleSet::build(Some(&["per_seq_gc".to_string(), "basic_stats".to_string()])).unwrap();
        assert_eq!(set.modules.len(), 2);
        assert_eq!(set.modules[0].name(), "basic_stats");
        assert_eq!(set.modules[1].name(), "per_seq_gc");
    }

    #[test]
    fn end_to_end_tidy_batch() {
        let mut set = ModuleSet::build(None).unwrap();
        set.update_batch(&batch(&["ACGT", "GGCC"], &["IIII", "IIII"]))
            .unwrap();
        let out = set.finalize().unwrap();
        assert_eq!(out.schema(), tidy_schema());
        assert!(out.num_rows() > 0);
    }
}
