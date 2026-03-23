//! Roundtrip correctness test with fjall backend:
//! annotate 1000 variants (fjall KV) → write VCF → read back → compare with golden VEP 115.
//!
//! Same verification as vcf_roundtrip_golden.rs but using the fjall KV backend
//! for variation lookup + SIFT instead of parquet.
//!
//! The fjall caches are generated from the trimmed parquet fixtures on first run.

use std::sync::Arc;

use datafusion::arrow::array::Array;
use datafusion::prelude::*;
use datafusion_bio_format_vcf::VcfCompressionType;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::{register_vep_functions, vcf_sink};

/// Resolve a path relative to the workspace root.
fn workspace_path(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

/// Build fjall variation cache from parquet if it doesn't exist.
async fn ensure_fjall_variation(parquet_dir: &str, fjall_path: &str) {
    if std::path::Path::new(fjall_path).exists() {
        return;
    }
    eprintln!("Building fjall variation cache from {parquet_dir}...");
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_parquet("vep_source", parquet_dir, Default::default())
        .await
        .unwrap();
    let loader = datafusion_bio_function_vep::kv_cache::CacheLoader::new("vep_source", fjall_path)
        .with_parallelism(1);
    let stats = loader.load(&ctx).await.unwrap();
    eprintln!(
        "  variation fjall: {} variants, {} positions",
        stats.total_variants, stats.total_positions
    );
}

/// Build fjall sift cache from parquet if it doesn't exist.
fn ensure_fjall_sift(parquet_path: &str, fjall_path: &str) {
    if std::path::Path::new(fjall_path).exists() {
        return;
    }
    eprintln!("Building fjall sift cache from {parquet_path}...");
    // Use the load_sift_cache example logic inline.
    use datafusion::arrow::array::{
        Float32Array, Int32Array, ListArray, StringArray, StringViewArray,
    };
    use datafusion::common::DataFusionError;
    use datafusion_bio_function_vep::kv_cache::sift_store::SiftKvStore;
    use datafusion_bio_function_vep::transcript_consequence::{
        CachedPredictions, CompactPrediction,
    };
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let db = fjall::Database::builder(fjall_path)
        .cache_size(64 * 1024 * 1024)
        .manual_journal_persist(true)
        .open()
        .unwrap();
    let sift_store = SiftKvStore::create(&db).unwrap();

    let file = std::fs::File::open(parquet_path).unwrap();
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap()
        .build()
        .unwrap();

    let mut total = 0u64;
    for batch_result in reader {
        let batch = batch_result.unwrap();
        let schema = batch.schema();
        let Some(tid_idx) = schema.index_of("transcript_id").ok() else {
            continue;
        };
        let sift_idx = schema.index_of("sift_predictions").ok();
        let poly_idx = schema.index_of("polyphen_predictions").ok();

        for row in 0..batch.num_rows() {
            let tid_col = batch.column(tid_idx);
            let transcript_id = if let Some(a) = tid_col.as_any().downcast_ref::<StringArray>() {
                if a.is_null(row) {
                    continue;
                } else {
                    a.value(row)
                }
            } else if let Some(a) = tid_col.as_any().downcast_ref::<StringViewArray>() {
                if a.is_null(row) {
                    continue;
                } else {
                    a.value(row)
                }
            } else {
                continue;
            };

            let read_preds = |col_idx: Option<usize>| -> Vec<CompactPrediction> {
                let Some(idx) = col_idx else {
                    return vec![];
                };
                let col = batch.column(idx);
                if col.is_null(row) {
                    return vec![];
                }
                let Some(list_arr) = col.as_any().downcast_ref::<ListArray>() else {
                    return vec![];
                };
                let offsets = list_arr.offsets();
                let start = offsets[row] as usize;
                let end = offsets[row + 1] as usize;
                if start == end {
                    return vec![];
                }
                let values = list_arr.values();
                let struct_arr = values
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::StructArray>()
                    .unwrap();
                let positions = struct_arr.column_by_name("position").unwrap();
                let amino_acids = struct_arr.column_by_name("amino_acid").unwrap();
                let predictions = struct_arr.column_by_name("prediction").unwrap();
                let scores = struct_arr.column_by_name("score").unwrap();
                let pos_arr = positions.as_any().downcast_ref::<Int32Array>();
                let aa_arr = amino_acids.as_any().downcast_ref::<StringArray>();
                let aa_view = amino_acids.as_any().downcast_ref::<StringViewArray>();
                let pred_arr = predictions.as_any().downcast_ref::<StringArray>();
                let pred_view = predictions.as_any().downcast_ref::<StringViewArray>();
                let score_arr = scores.as_any().downcast_ref::<Float32Array>();
                let mut out = Vec::with_capacity(end - start);
                for i in start..end {
                    let Some(pos) =
                        pos_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
                    else {
                        continue;
                    };
                    let aa = aa_arr
                        .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
                        .or_else(|| {
                            aa_view.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
                        });
                    let pred = pred_arr
                        .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
                        .or_else(|| {
                            pred_view
                                .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
                        });
                    let score =
                        score_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
                    if let (Some(aa), Some(pred), Some(sc)) = (aa, pred, score) {
                        if let Some(aa_idx) = CompactPrediction::encode_amino_acid(aa) {
                            out.push(CompactPrediction {
                                position: pos,
                                amino_acid: aa_idx,
                                prediction: CompactPrediction::encode_prediction(pred),
                                score: sc,
                            });
                        }
                    }
                }
                out
            };

            let mut preds = CachedPredictions::default();
            preds.sift = read_preds(sift_idx);
            preds.polyphen = read_preds(poly_idx);
            if preds.sift.is_empty() && preds.polyphen.is_empty() {
                continue;
            }
            preds.sort();
            total += 1;
            sift_store.put(transcript_id, &preds).unwrap();
        }
    }
    db.persist(fjall::PersistMode::SyncAll).unwrap();
    eprintln!("  sift fjall: {total} transcripts");
}

/// Create symlinks for context tables if they don't exist.
fn ensure_context_symlinks(fjall_dir: &std::path::Path, parquet_dir: &std::path::Path) {
    for dir in [
        "variation",
        "transcript",
        "exon",
        "translation_core",
        "translation_sift",
        "regulatory",
        "motif",
    ] {
        let link = fjall_dir.join(dir);
        if !link.exists() {
            let target = parquet_dir.join(dir);
            if target.exists() {
                std::os::unix::fs::symlink(&target, &link).ok();
            }
        }
    }
}

#[tokio::test(flavor = "multi_thread")]
async fn test_roundtrip_golden_fjall_all_column_values() {
    let input_vcf = workspace_path("vep-benchmark/data/golden/input_1000.vcf");
    let golden_vcf = workspace_path("vep-benchmark/data/golden/golden_1000_vep115.vcf");
    let parquet_cache = workspace_path("vep-benchmark/data/golden/cache");
    let fjall_cache = workspace_path("vep-benchmark/data/golden/fjall_cache");
    let ref_fasta = workspace_path("vep-benchmark/data/golden/reference_chr1.fa");

    if !input_vcf.exists() || !golden_vcf.exists() {
        eprintln!("Skipping: test fixtures not found");
        return;
    }

    // Ensure fjall caches exist (generated from parquet on first run).
    std::fs::create_dir_all(&fjall_cache).ok();
    let var_fjall = fjall_cache.join("variation.fjall");
    let sift_fjall = fjall_cache.join("translation_sift.fjall");
    ensure_fjall_variation(
        parquet_cache.join("variation").to_str().unwrap(),
        var_fjall.to_str().unwrap(),
    )
    .await;
    ensure_fjall_sift(
        parquet_cache
            .join("translation_sift/chr1.parquet")
            .to_str()
            .unwrap(),
        sift_fjall.to_str().unwrap(),
    );
    ensure_context_symlinks(&fjall_cache, &parquet_cache);

    let input_vcf = input_vcf.to_str().unwrap();
    let golden_vcf = golden_vcf.to_str().unwrap();
    let fjall_cache = fjall_cache.to_str().unwrap();
    let ref_fasta = ref_fasta.to_str().unwrap();

    // ── Step 1: Annotate with fjall backend and write to VCF ──
    let tmp_dir = tempfile::TempDir::new().unwrap();
    let output_path = tmp_dir.path().join("annotated_fjall.vcf");

    let ctx = SessionContext::new();
    register_vep_functions(&ctx);
    let input_str = input_vcf.to_string();
    let vcf_provider = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(input_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    ctx.register_table("vcf", Arc::new(vcf_provider)).unwrap();

    let options = format!(
        "{{\"partitioned\":true,\"everything\":true,\"extended_probes\":true,\
         \"use_fjall\":true,\"reference_fasta_path\":\"{ref_fasta}\"}}"
    );
    let rows_written = vcf_sink::annotate_to_vcf(
        &ctx,
        "vcf",
        fjall_cache,
        "parquet",
        Some(&options),
        &output_path,
        VcfCompressionType::Plain,
    )
    .await
    .unwrap();
    assert_eq!(rows_written, 1000, "Should write 1000 annotated rows");

    // ── Step 2: Read input, output, and golden ──
    let ctx2 = SessionContext::new();

    let input_str2 = input_vcf.to_string();
    let input_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(input_str2, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    let out_path_str = output_path.display().to_string();
    let output_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(out_path_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();
    let golden_str = golden_vcf.to_string();
    let golden_prov = tokio::task::spawn_blocking(move || {
        VcfTableProvider::new(golden_str, None, None, None, false).unwrap()
    })
    .await
    .unwrap();

    ctx2.register_table("input_vcf", Arc::new(input_prov))
        .unwrap();
    ctx2.register_table("output_vcf", Arc::new(output_prov))
        .unwrap();
    ctx2.register_table("golden_vcf", Arc::new(golden_prov))
        .unwrap();

    let input_batches = ctx2
        .sql("SELECT * FROM input_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let output_batches = ctx2
        .sql("SELECT * FROM output_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();
    let golden_batches = ctx2
        .sql("SELECT * FROM golden_vcf ORDER BY start")
        .await
        .unwrap()
        .collect()
        .await
        .unwrap();

    let input_rows: usize = input_batches.iter().map(|b| b.num_rows()).sum();
    let output_rows: usize = output_batches.iter().map(|b| b.num_rows()).sum();
    let golden_rows: usize = golden_batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(output_rows, 1000);
    assert_eq!(golden_rows, 1000);
    assert_eq!(input_rows, 1000);

    let input_batch =
        datafusion::arrow::compute::concat_batches(&input_batches[0].schema(), &input_batches)
            .unwrap();
    let output_batch =
        datafusion::arrow::compute::concat_batches(&output_batches[0].schema(), &output_batches)
            .unwrap();
    let golden_batch =
        datafusion::arrow::compute::concat_batches(&golden_batches[0].schema(), &golden_batches)
            .unwrap();

    // ── Step 3: Core VCF columns ──
    for col_name in [
        "chrom", "start", "end", "ref", "alt", "id", "qual", "filter",
    ] {
        if let (Ok(i), Ok(o)) = (
            input_batch.schema().index_of(col_name),
            output_batch.schema().index_of(col_name),
        ) {
            assert_eq!(
                input_batch.column(i).as_ref(),
                output_batch.column(o).as_ref(),
                "Core column '{col_name}' differs (fjall backend)"
            );
        }
    }

    // ── Step 4: INFO columns ──
    let info_columns: Vec<String> = input_batch
        .schema()
        .fields()
        .iter()
        .filter(|f| {
            f.metadata()
                .get("bio.vcf.field.field_type")
                .is_some_and(|v| v == "INFO")
        })
        .map(|f| f.name().clone())
        .collect();
    for col_name in &info_columns {
        assert!(
            output_batch.schema().index_of(col_name).is_ok(),
            "INFO '{col_name}' missing (fjall)"
        );
        let in_col = input_batch.column(input_batch.schema().index_of(col_name).unwrap());
        let out_col = output_batch.column(output_batch.schema().index_of(col_name).unwrap());
        assert_eq!(
            in_col.as_ref(),
            out_col.as_ref(),
            "INFO '{col_name}' values differ (fjall)"
        );
    }

    // ── Step 5: FORMAT columns ──
    let format_columns: Vec<String> = input_batch
        .schema()
        .fields()
        .iter()
        .filter(|f| {
            f.metadata()
                .get("bio.vcf.field.field_type")
                .is_some_and(|v| v == "FORMAT")
        })
        .map(|f| f.name().clone())
        .collect();
    for col_name in &format_columns {
        assert!(
            output_batch.schema().index_of(col_name).is_ok(),
            "FORMAT '{col_name}' missing (fjall)"
        );
        let in_col = input_batch.column(input_batch.schema().index_of(col_name).unwrap());
        let out_col = output_batch.column(output_batch.schema().index_of(col_name).unwrap());
        assert_eq!(
            in_col.as_ref(),
            out_col.as_ref(),
            "FORMAT '{col_name}' values differ (fjall)"
        );
    }

    // ── Step 6: CSQ vs golden ──
    let our_csq_name = if output_batch.schema().index_of("csq").is_ok() {
        "csq"
    } else if output_batch.schema().index_of("CSQ").is_ok() {
        "CSQ"
    } else {
        panic!("CSQ not found in fjall output")
    };
    let golden_csq_name = if golden_batch.schema().index_of("CSQ").is_ok() {
        "CSQ"
    } else if golden_batch.schema().index_of("csq").is_ok() {
        "csq"
    } else {
        panic!("CSQ not found in golden")
    };

    let our_csq = output_batch.column(output_batch.schema().index_of(our_csq_name).unwrap());
    let golden_csq = golden_batch.column(golden_batch.schema().index_of(golden_csq_name).unwrap());

    if our_csq.as_ref() == golden_csq.as_ref() {
        eprintln!("[fjall] CSQ: exact match (all 1000 rows identical)");
    } else {
        let get_str = |col: &dyn Array, row: usize| -> String {
            if col.is_null(row) {
                return String::new();
            }
            if let Some(a) = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringArray>()
            {
                return a.value(row).to_string();
            }
            if let Some(a) = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringViewArray>()
            {
                return a.value(row).to_string();
            }
            String::new()
        };
        let mut mismatched = 0;
        for row in 0..1000 {
            let ours = get_str(our_csq.as_ref(), row);
            let golden = get_str(golden_csq.as_ref(), row);
            if ours != golden {
                mismatched += 1;
                if mismatched <= 3 {
                    eprintln!("[fjall] CSQ mismatch row {row}:");
                    eprintln!("  golden: {}...", &golden[..golden.len().min(150)]);
                    eprintln!("  ours:   {}...", &ours[..ours.len().min(150)]);
                }
            }
        }
        assert_eq!(
            mismatched, 0,
            "[fjall] CSQ should have 0 mismatches vs golden VEP 115 ({mismatched} mismatched)"
        );
    }

    eprintln!(
        "[fjall] All checks passed: {} core, {} INFO, {} FORMAT columns + CSQ",
        8,
        info_columns.len(),
        format_columns.len()
    );
}
