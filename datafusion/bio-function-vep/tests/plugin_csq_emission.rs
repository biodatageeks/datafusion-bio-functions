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
use datafusion_bio_function_vep::plugin_cache::csq::{empty_suffix, field_suffix};
use datafusion_bio_function_vep::plugin_cache::registry::PluginRegistry;
use datafusion_bio_function_vep::plugin_cache::source_manifest::{
    MatchColumn, ValueColumn, ValueType,
};
use datafusion_bio_function_vep::plugin_cache::template::build_attr_namespace;

/// Namespace for a `C/G` variant with the given amino-acid change (rest empty).
fn ns_aa<'a>(amino_acids: &'a str, protein_pos: &'a str) -> [Option<&'a str>; 17] {
    build_attr_namespace(
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        protein_pos,
        amino_acids,
        "",
        "C",
        "G",
    )
}
use datafusion_bio_function_vep::plugin_cache::write::{PluginShardWriter, plugin_output_schema};

fn build_alphamissense_shard(cache_root: &std::path::Path) {
    let plugin_dir = cache_root.join("plugin").join("alphamissense");
    std::fs::create_dir_all(&plugin_dir).unwrap();
    let matches = vec![MatchColumn {
        column: "protein_variant".into(),
        template: "{ref_aa}{Protein_position}{alt_aa}".into(),
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
            template: "{ref_aa}{Protein_position}{alt_aa}".into(),
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
    let ns_missense = ns_aa("C/W", "17");
    let missense = slices.probe_all(22893742, "C/G", &ns_missense);
    assert_eq!(field_suffix(&missense), "|0.4833|ambiguous");

    // Transcript line 2: intron (no amino-acid change) → gate → empty fields.
    let ns_intron = ns_aa("", "");
    let intron = slices.probe_all(22893742, "C/G", &ns_intron);
    assert_eq!(field_suffix(&intron), empty_suffix(n));
    assert_eq!(field_suffix(&intron), "||");

    // A different protein change at the same position (wrong isoform) → miss.
    let ns_wrong = ns_aa("C/Y", "17");
    assert_eq!(
        field_suffix(&slices.probe_all(22893742, "C/G", &ns_wrong)),
        "||"
    );

    // A variant with no shard row (different position) → empty fields.
    let none_here = slices.probe_all(99999999, "A/G", &ns_missense);
    assert_eq!(field_suffix(&none_here), "||");
}

// A shared-anchor indel (VCF POS 100 `A>ATG`) is VEP-normalized to start 101,
// allele `-/TG`. A plugin cache keyed like the variation cache stores the row at
// the normalized start, so the lookup must use `input_start` (101), not the raw
// POS (100). Guards the PR #190 C1 fix (annotate_provider keys the take + probe
// on `input_start`); the old raw-`start_val` code would have missed here.
#[tokio::test(flavor = "multi_thread")]
async fn indel_probe_uses_normalized_start() {
    let dir = tempfile::tempdir().unwrap();
    let cache_root = dir.path();
    let plugin_dir = cache_root.join("plugin").join("demo");
    std::fs::create_dir_all(&plugin_dir).unwrap();

    // Per-variant plugin (no match discriminator).
    let matches: Vec<MatchColumn> = vec![];
    let vals = vec![ValueColumn {
        column: "score".into(),
        csq_field: "SCORE".into(),
        ty: ValueType::Float32,
    }];
    let schema = plugin_output_schema(&matches, &vals);
    // One row at the NORMALIZED coordinates: start 101, allele "-/TG".
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["1"])),
            Arc::new(UInt32Array::from(vec![101u32])),
            Arc::new(UInt32Array::from(vec![101u32])),
            Arc::new(StringArray::from(vec!["-/TG"])),
            Arc::new(Float32Array::from(vec![0.75f32])),
            Arc::new(Int8Array::from(vec![1i8])),
        ],
    )
    .unwrap();
    let mut w = PluginShardWriter::create(&plugin_dir.join("chr1.parquet"), schema).unwrap();
    w.write(&batch).unwrap();
    w.finish().unwrap();

    let manifest = CacheManifest {
        plugin_name: "demo".into(),
        source_manifest: "demo.source.toml".into(),
        key_columns: vec![
            "chrom".into(),
            "start".into(),
            "end".into(),
            "allele_string".into(),
        ],
        match_columns: vec![],
        value_columns: vec![ValueColumnRecord {
            column: "score".into(),
            csq_field: "SCORE".into(),
            ty: "Float32".into(),
        }],
        chroms: vec![ChromEntry {
            chrom: "chr1".into(),
            file: "chr1.parquet".into(),
            rows: 1,
            warm: 0,
            cold: 1,
        }],
        cache_source_version: None,
    };
    manifest.write(&plugin_dir).unwrap();

    let reg = PluginRegistry::open(cache_root, "1").await.unwrap();

    // Normalized start (101) hits.
    let hit = reg.take_buffer_all(&[101]).await.unwrap();
    assert_eq!(field_suffix(&hit.probe_all(101, "-/TG", &[])), "|0.75");
    // Raw VCF POS (100) misses — this is what the pre-fix code used.
    let miss = reg.take_buffer_all(&[100]).await.unwrap();
    assert_eq!(
        field_suffix(&miss.probe_all(100, "-/TG", &[])),
        empty_suffix(1)
    );
}
