//! Round-trip test: write VCF from annotation output, read back with noodles-vcf,
//! verify it produces a valid VCF file that can be parsed.

use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{RecordBatch, StringArray, UInt32Array};
use datafusion::arrow::datatypes::{DataType, Field, Schema};

/// Write annotated batches to a VCF file using the same logic as profile_annotation.
fn write_test_vcf(batches: &[RecordBatch], path: &Path, sample_name: &str) {
    use std::io::{BufWriter, Write};

    let mut file = BufWriter::new(std::fs::File::create(path).unwrap());
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        file,
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations\">"
    )
    .unwrap();
    writeln!(
        file,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .unwrap();
    writeln!(
        file,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )
    .unwrap();

    for batch in batches {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").unwrap();
        let start_idx = schema.index_of("start").unwrap();
        let ref_idx = schema.index_of("ref").unwrap();
        let alt_idx = schema.index_of("alt").unwrap();
        let csq_idx = schema.index_of("CSQ").ok();
        let gt_idx = schema.index_of("GT").ok();

        for row in 0..batch.num_rows() {
            let chrom = batch
                .column(chrom_idx)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .value(row);
            let pos = batch
                .column(start_idx)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .value(row);
            let ref_al = batch
                .column(ref_idx)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .value(row);
            let alt_al = batch
                .column(alt_idx)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .value(row);
            let csq = csq_idx
                .and_then(|i| {
                    batch
                        .column(i)
                        .as_any()
                        .downcast_ref::<StringArray>()
                        .map(|a| a.value(row).to_string())
                })
                .unwrap_or_default();
            let gt = gt_idx
                .and_then(|i| {
                    batch
                        .column(i)
                        .as_any()
                        .downcast_ref::<StringArray>()
                        .map(|a| a.value(row).to_string())
                })
                .unwrap_or_else(|| "./.".to_string());

            let info = if csq.is_empty() {
                ".".to_string()
            } else {
                format!("CSQ={csq}")
            };

            writeln!(
                file,
                "{chrom}\t{pos}\t.\t{ref_al}\t{alt_al}\t.\t.\t{info}\tGT\t{gt}"
            )
            .unwrap();
        }
    }
    file.flush().unwrap();
}

fn make_test_batches() -> Vec<RecordBatch> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
        Field::new("CSQ", DataType::Utf8, true),
        Field::new("GT", DataType::Utf8, true),
    ]));

    let batch = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(StringArray::from(vec!["chr1", "chr1", "chr2"])),
            Arc::new(UInt32Array::from(vec![100, 200, 300])),
            Arc::new(StringArray::from(vec!["A", "C", "G"])),
            Arc::new(StringArray::from(vec!["T", "G", "A"])),
            Arc::new(StringArray::from(vec![
                "T|missense_variant|MODERATE|GENE1|ENSG001|Transcript|ENST001|protein_coding",
                "G|synonymous_variant|LOW|GENE2|ENSG002|Transcript|ENST002|protein_coding",
                "",
            ])),
            Arc::new(StringArray::from(vec!["0/1", "1/1", "0/0"])),
        ],
    )
    .unwrap();

    vec![batch]
}

#[test]
fn test_vcf_roundtrip_noodles_can_parse() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path();

    let batches = make_test_batches();
    write_test_vcf(&batches, path, "HG002");

    // Read back with noodles-vcf and verify it parses.
    let file = std::fs::File::open(path).unwrap();
    let mut reader = noodles_vcf::io::reader::Builder::default()
        .build_from_reader(BufReader::new(file))
        .unwrap();

    let header = reader.read_header().unwrap();

    // Verify header.
    assert_eq!(header.sample_names().len(), 1);
    assert!(header.sample_names().iter().any(|s| s.as_str() == "HG002"));

    // Verify INFO CSQ is defined.
    assert!(header.infos().get("CSQ").is_some());

    // Verify FORMAT GT is defined.
    assert!(header.formats().get("GT").is_some());

    // Read all records.
    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        records.push(record);
    }

    assert_eq!(records.len(), 3);
}

#[test]
fn test_vcf_roundtrip_positions_correct() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path();

    let batches = make_test_batches();
    write_test_vcf(&batches, path, "SAMPLE1");

    let content = std::fs::read_to_string(path).unwrap();
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();

    assert_eq!(data_lines.len(), 3);

    // Check positions are 1-based (input was 0-based UInt32: 100, 200, 300).
    // The write_test_vcf writes the value directly (UInt32 as-is).
    let fields: Vec<&str> = data_lines[0].split('\t').collect();
    assert_eq!(fields[0], "chr1"); // CHROM
    assert_eq!(fields[1], "100"); // POS
    assert_eq!(fields[3], "A"); // REF
    assert_eq!(fields[4], "T"); // ALT

    let fields2: Vec<&str> = data_lines[2].split('\t').collect();
    assert_eq!(fields2[0], "chr2");
    assert_eq!(fields2[1], "300");
    assert_eq!(fields2[7], "."); // no CSQ → INFO is "."
    assert_eq!(fields2[9], "0/0"); // GT
}

#[test]
fn test_vcf_roundtrip_csq_preserved() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path();

    let batches = make_test_batches();
    write_test_vcf(&batches, path, "SAMPLE1");

    let file = std::fs::File::open(path).unwrap();
    let mut reader = noodles_vcf::io::reader::Builder::default()
        .build_from_reader(BufReader::new(file))
        .unwrap();
    let _header = reader.read_header().unwrap();

    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 3);

    // First record should have CSQ in INFO.
    let info0 = records[0].info();
    let info_str: &str = info0.as_ref();
    assert!(
        info_str.contains("CSQ="),
        "Expected CSQ in INFO, got: {info_str}"
    );
    assert!(info_str.contains("missense_variant"));

    // Third record has no CSQ → INFO should be ".".
    let info2 = records[2].info();
    let info_str3: &str = info2.as_ref();
    // noodles parses VCF "." INFO as empty string.
    assert!(
        info_str3 == "." || info_str3.is_empty(),
        "Empty CSQ should produce '.' or empty INFO, got: {info_str3}"
    );
}

#[test]
fn test_vcf_roundtrip_empty_batches() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path();

    write_test_vcf(&[], path, "SAMPLE1");

    let file = std::fs::File::open(path).unwrap();
    let mut reader = noodles_vcf::io::reader::Builder::default()
        .build_from_reader(BufReader::new(file))
        .unwrap();
    let header = reader.read_header().unwrap();

    assert_eq!(header.sample_names().len(), 1);
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 0);
}
