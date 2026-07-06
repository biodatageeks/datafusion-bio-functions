//! Integration test for the per-transcript plugin CSQ emission chain (Task I1):
//! `PluginRegistry::open` → `take_buffer_all` → `probe_all` → `csq::field_suffix`,
//! exactly as `annotate_batch_with_transcript_engine` drives it. Verifies the
//! per-transcript gate: a matching amino-acid change emits the plugin fields; a
//! non-missense line (no aa-change) and a missing position emit empty fields of
//! the same width — reproducing the AlphaMissense golden's behavior.
#![cfg(feature = "parquet-cache")]

use std::sync::Arc;

use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion_bio_function_vep::plugin_cache::cache_manifest::{
    CacheManifest, ChromEntry, MatchColumnRecord, ValueColumnRecord,
};
use datafusion_bio_function_vep::plugin_cache::csq::{
    amino_acid_change, empty_suffix, field_suffix,
};
use datafusion_bio_function_vep::plugin_cache::registry::{EngineAttrs, PluginRegistry};
use datafusion_bio_function_vep::plugin_cache::source_manifest::{
    MatchColumn, ValueColumn, ValueType,
};
use datafusion_bio_function_vep::plugin_cache::write::{PluginShardWriter, plugin_output_schema};

fn build_alphamissense_shard(cache_root: &std::path::Path) {
    let plugin_dir = cache_root.join("plugin").join("alphamissense");
    std::fs::create_dir_all(&plugin_dir).unwrap();
    let matches = vec![MatchColumn {
        column: "protein_variant".into(),
        engine_attr: "amino_acid_change".into(),
    }];
    let vals = vec![
        ValueColumn {
            column: "am_pathogenicity".into(),
            csq_field: "am_pathogenicity".into(),
            ty: ValueType::Float32,
        },
        ValueColumn {
            column: "am_class".into(),
            csq_field: "am_class".into(),
            ty: ValueType::Utf8,
        },
    ];
    let schema = plugin_output_schema(&matches, &vals);
    // chr22:22893742 C>G, missense W (protein_variant "C17W"): am 0.4833 / ambiguous.
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["22"])),
            Arc::new(UInt32Array::from(vec![22893742u32])),
            Arc::new(UInt32Array::from(vec![22893742u32])),
            Arc::new(StringArray::from(vec!["C/G"])),
            Arc::new(StringArray::from(vec!["C17W"])),
            Arc::new(Float32Array::from(vec![0.4833f32])),
            Arc::new(StringArray::from(vec!["ambiguous"])),
            Arc::new(Int8Array::from(vec![1i8])),
        ],
    )
    .unwrap();
    let mut w = PluginShardWriter::create(&plugin_dir.join("chr22.parquet"), schema).unwrap();
    w.write(&batch).unwrap();
    w.finish().unwrap();

    let manifest = CacheManifest {
        plugin_name: "alphamissense".into(),
        source_manifest: "alphamissense.source.toml".into(),
        key_columns: vec![
            "chrom".into(),
            "start".into(),
            "end".into(),
            "allele_string".into(),
        ],
        match_columns: vec![MatchColumnRecord {
            column: "protein_variant".into(),
            engine_attr: "amino_acid_change".into(),
        }],
        value_columns: vec![
            ValueColumnRecord {
                column: "am_pathogenicity".into(),
                csq_field: "am_pathogenicity".into(),
                ty: "Float32".into(),
            },
            ValueColumnRecord {
                column: "am_class".into(),
                csq_field: "am_class".into(),
                ty: "Utf8".into(),
            },
        ],
        chroms: vec![ChromEntry {
            chrom: "chr22".into(),
            file: "chr22.parquet".into(),
            rows: 1,
            warm: 0,
            cold: 1,
        }],
        cache_source_version: None,
    };
    manifest.write(&plugin_dir).unwrap();
}

/// Reproduce the per-transcript emission for one buffer: the exact chain the
/// engine runs. Two transcript lines at the same variant — one missense (aa
/// change matches), one intron (no aa change) — plus a variant with no shard row.
#[tokio::test(flavor = "multi_thread")]
async fn plugin_csq_gates_per_transcript() {
    let dir = tempfile::tempdir().unwrap();
    build_alphamissense_shard(dir.path());

    let reg = PluginRegistry::open(dir.path(), "22").await.unwrap();
    let n = reg.csq_fields().len();
    assert_eq!(reg.csq_fields(), vec!["am_pathogenicity", "am_class"]);

    // One buffer take covering the variant position.
    let slices = reg.take_buffer_all(&[22893742]).await.unwrap();

    // Transcript line 1: missense C17W (Amino_acids "C/W", Protein_position "17").
    let attrs_missense = EngineAttrs {
        amino_acid_change: amino_acid_change("C/W", "17"),
    };
    let missense = slices.probe_all(22893742, "C/G", &attrs_missense);
    assert_eq!(field_suffix(&missense), "|0.4833|ambiguous");

    // Transcript line 2: intron (no amino-acid change) → gate → empty fields.
    let attrs_intron = EngineAttrs {
        amino_acid_change: amino_acid_change("", ""),
    };
    let intron = slices.probe_all(22893742, "C/G", &attrs_intron);
    assert_eq!(field_suffix(&intron), empty_suffix(n));
    assert_eq!(field_suffix(&intron), "||");

    // A different protein change at the same position (wrong isoform) → miss.
    let attrs_wrong = EngineAttrs {
        amino_acid_change: amino_acid_change("C/Y", "17"),
    };
    assert_eq!(
        field_suffix(&slices.probe_all(22893742, "C/G", &attrs_wrong)),
        "||"
    );

    // A variant with no shard row (different position) → empty fields.
    let none_here = slices.probe_all(99999999, "A/G", &attrs_missense);
    assert_eq!(field_suffix(&none_here), "||");
}
