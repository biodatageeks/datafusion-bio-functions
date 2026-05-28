//! v2-paradigm port of `ensembl-vep/t/Parser.t`.
//!
//! Detailed plan: `porting-tests/detailed_plans/Parser.md` (audited
//! 2026-05-27, post-Phase-D with 3 Axis-B additions).
//! TDD task plan:  implemented in-context 2026-05-28 per port-to-vepyr
//! skill heuristic (12 unit-ports + 47 architectural-no-analogue + 21
//! blocked-future-work — mostly documentation, no subagent needed).
//!
//! v2 paradigm anchors (~/.claude/skills/port-to-vepyr/references/v2-paradigm.md):
//! - Sztywno 1:1 — every Perl subtest gets a Rust analogue (here:
//!   passing test, `// SUBTEST #N:` documentation comment naming the
//!   missing-by-design vepyr component, or a commented-out future-work
//!   stub).
//! - Standalone — no docker dependency, no `golden.vcf`, no
//!   `port_common::run_and_compare_csq`. Hand-coded assertion values
//!   suffice because the unit-ports in this file pin parser-layer
//!   behaviour (file-exists, gzip-detection, CRLF-tolerance, bgzip
//!   reading), not CSQ-output values.
//!
//! Where the rest of the Perl rows live:
//! - allele-minimisation rows #17, #18, #19, Axis B B1 →
//!   `src/allele.rs::tests` (pure-function tests on `vcf_to_vep_allele`).
//! - `_have_chr` rows #31, #32, #34 →
//!   `src/transcript_consequence.rs::tests` (against `normalize_chrom`).
//! - Axis B B2 (`buffer_variant_bounds` multi-chrom first-chrom-wins) →
//!   already covered by
//!   `src/annotate_provider.rs::tests::port_input_buffer_axisB1_buffer_variant_bounds_multi_chrom_binds_to_first_chrom`
//!   (the existing InputBuffer Axis B B1 test asserts the same vepyr
//!   invariant; cross-referenced here per detailed_plan).
//!
//! Coverage parity (per detailed_plan):
//! - 9 unit-port + 3 Axis B = 12 passing Rust tests across the 3 files
//!   listed above.
//! - 47 architectural-no-analogue rows + 21 blocked-future-work rows
//!   documented in the comment block at the bottom of this file plus
//!   inline at the unit-port test sites.
//! - Phase-A parity: 9 / 77 = 12% Perl-denominator unit-ports;
//!   (9 + 21) / 77 = 39% addressable once engine blockers land.

use std::io::Write;
use std::sync::Arc;

use datafusion::prelude::*;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use flate2::Compression;
use flate2::write::GzEncoder;

/// Tiny chr21 SNV row used as the body of all inline VCFs in this file.
/// Same shape as `tests/data/port/test_smoke/input.vcf` row 1.
const SMOKE_VCF: &str = "\
##fileformat=VCFv4.2
##source=vepyr-port-parser
##contig=<ID=21,length=46709983>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t25587759\ttest_smoke_snv\tT\tC\t.\tPASS\t.
";

/// Build a session context with the given file-path-backed VCF provider
/// registered as table `vcf`, then count rows via SQL. Returns the row
/// count; surfaces parse errors to the caller so `try_new`-level
/// failures propagate.
async fn count_vcf_rows(file_path: &str) -> datafusion::error::Result<usize> {
    let provider = VcfTableProvider::new(file_path.to_string(), None, None, None, false)?;
    let ctx = SessionContext::new();
    ctx.register_table("vcf", Arc::new(provider))?;
    let df = ctx.sql("SELECT COUNT(*) FROM vcf").await?;
    let batches = df.collect().await?;
    assert_eq!(batches.len(), 1);
    let batch = &batches[0];
    let col = batch.column(0);
    let arr = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .expect("COUNT(*) returns Int64");
    Ok(arr.value(0) as usize)
}

// ─────────────────────────────────────────────────────────────────────
// Unit-port rows
// ─────────────────────────────────────────────────────────────────────

/// Port of `Parser.t` subtest #7 (Perl l.52):
///   `throws_ok { $p->file('foo') } qr/File.+does not exist/`.
///
/// Vepyr analogue: `VcfTableProvider::new(<nonexistent_path>, ...)` is
/// observably an error — it panics inside `bio-format-core` object
/// storage (`object_storage.rs:183` unwrap on `kind: NotFound`) rather
/// than returning a typed `DataFusionError`. The contract Perl
/// asserts ("opening a missing path fails") holds; only the failure
/// mode differs (panic vs `Err`). Pinning the panic here documents the
/// current behaviour; a future patch can convert the unwrap to a
/// proper `DataFusionError::IoError` without breaking this test (just
/// remove the `catch_unwind` wrapper).
#[tokio::test(flavor = "multi_thread")]
async fn port_parser_subtest_07_nonexistent_path_errors() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let bogus = tmp.path().join("does-not-exist.vcf");
    // The temp dir exists but the file path does not.
    assert!(!bogus.exists());

    let bogus_str = bogus.to_str().unwrap().to_string();
    // `VcfTableProvider::new` is synchronous and panics on NotFound
    // today; use `catch_unwind` on a blocking call so the test passes
    // either way.
    let result = std::panic::catch_unwind(|| {
        VcfTableProvider::new(bogus_str.clone(), None, None, None, false)
    });
    match result {
        Ok(Ok(_)) => panic!("VcfTableProvider::new must not succeed on a nonexistent path"),
        Ok(Err(err)) => {
            // Future "fixed" world: typed error surfaces. Check it is
            // non-empty so the contract message is still informative.
            let msg = format!("{err}");
            assert!(!msg.is_empty(), "typed error must carry a message");
        }
        Err(_panic) => {
            // Current world: object-storage unwraps on NotFound. The
            // contract (observable failure on missing path) is upheld.
        }
    }
}

/// Port of `Parser.t` subtest #12 (Perl l.62):
///   `ok($p->file($test_gzvcf) =~ /GLOB/, 'gzipped VCF')`.
///
/// Vepyr analogue: `VcfTableProvider::new(<*.vcf.gz>)` succeeds and the
/// row count matches the uncompressed equivalent — i.e. gzip is
/// transparently autodetected by `noodles_vcf` (via the magic bytes
/// 0x1f 0x8b) and the reader streams decompressed records.
#[tokio::test(flavor = "multi_thread")]
async fn port_parser_subtest_12_gzipped_vcf_reads() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let gz_path = tmp.path().join("smoke.vcf.gz");
    {
        let f = std::fs::File::create(&gz_path).expect("create gz file");
        let mut encoder = GzEncoder::new(f, Compression::default());
        encoder
            .write_all(SMOKE_VCF.as_bytes())
            .expect("write gzipped VCF body");
        encoder.finish().expect("finish gz stream");
    }
    // Sanity-check: file must look gzip-magic 0x1f 0x8b on disk.
    let raw = std::fs::read(&gz_path).expect("read gz");
    assert_eq!(
        &raw[..2],
        &[0x1f, 0x8b],
        "gz file must start with gzip magic"
    );

    let rows = count_vcf_rows(gz_path.to_str().unwrap())
        .await
        .expect("gzipped VCF must read");
    assert_eq!(rows, 1, "smoke gzipped VCF has exactly one data row");
}

/// Port of `Parser.t` subtest #57 (Perl l.249-250):
///   `$p->file($windows_vcf); is($p->detect_format, 'vcf')`.
///
/// Vepyr analogue: `VcfTableProvider::new(<CRLF-VCF>)` must succeed and
/// the parser must tolerate `\r\n` line endings. Build a CRLF-line-
/// terminated tiny VCF inline and assert one row reads back. There is
/// no `detect_format` API in vepyr (single-format by design); the
/// behavioural equivalent is "the reader doesn't choke on CRLF input."
#[tokio::test(flavor = "multi_thread")]
async fn port_parser_subtest_57_crlf_line_endings_tolerated() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let crlf_path = tmp.path().join("smoke_crlf.vcf");
    let crlf_body = SMOKE_VCF.replace('\n', "\r\n");
    assert!(crlf_body.contains("\r\n"), "fixture must contain CRLF");
    std::fs::write(&crlf_path, crlf_body.as_bytes()).expect("write crlf vcf");

    let rows = count_vcf_rows(crlf_path.to_str().unwrap())
        .await
        .expect("CRLF VCF must read");
    assert_eq!(rows, 1, "CRLF-terminated VCF round-trips one row");
}

/// Port of Axis-B row B3 (post-Phase-D 2026-05-27):
///   `VcfTableProvider::try_new` with a `.vcf.bgz`-extension file
///   (bgzip not gzip) succeeds.
///
/// Detailed plan: `porting-tests/detailed_plans/Parser.md` row B3. Pin
/// extension-routing behaviour: a `.vcf.bgz` file containing a valid
/// gzip-magic stream must round-trip through the reader. We use a plain
/// gzip stream (not block-gzip) because the test only pins extension
/// handling; the underlying bgzf machinery is the responsibility of
/// `noodles-vcf` / `noodles-bgzf` and is exercised elsewhere.
#[tokio::test(flavor = "multi_thread")]
#[allow(non_snake_case)]
async fn port_parser_axisB3_vcf_bgz_extension_reads() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let bgz_path = tmp.path().join("smoke.vcf.bgz");
    {
        let f = std::fs::File::create(&bgz_path).expect("create bgz file");
        let mut encoder = GzEncoder::new(f, Compression::default());
        encoder
            .write_all(SMOKE_VCF.as_bytes())
            .expect("write gz body under .bgz extension");
        encoder.finish().expect("finish gz stream");
    }
    assert!(bgz_path.exists());
    assert_eq!(
        bgz_path.extension().and_then(|e| e.to_str()),
        Some("bgz"),
        "extension routing precondition"
    );

    // noodles_vcf accepts both `.gz` and `.bgz` as compressed inputs.
    let rows = count_vcf_rows(bgz_path.to_str().unwrap())
        .await
        .expect(".vcf.bgz must read");
    assert_eq!(rows, 1);
}

// ─────────────────────────────────────────────────────────────────────
// Blocked-future-work stubs (commented-out, in source order)
//
// Each stub names the missing vepyr API; the substantive entry lives in
// `porting-tests/future-work-vepyr.md`. These 21 rows raise vepyr's
// engine-level surface to match Parser.t once the listed APIs land.
// ─────────────────────────────────────────────────────────────────────

// SUBTESTS #20-#25 (Perl l.106-111): `get_SO_term('INS'/'DEL'/'INV'/
//   'TREP'/'BND'/'DEL_ME')` → 'insertion'/'deletion'/'inversion'/
//   'tandem_repeat'/'chromosome_breakpoint'/'mobile_element_deletion'.
// Blocked: engine blocker #2 — `sv_type_to_so` function + extended
// `SoTerm` variants. See `porting-tests/future-work-vepyr.md` SO-term
// entries and `port-status.md` engine blocker #2 cluster.
//
// #[test]
// fn port_parser_subtest_20_so_term_ins_to_insertion() {
//     assert_eq!(sv_type_to_so("INS"), Some(SoTerm::Insertion));
// }
// (similar for DEL, INV, TREP, BND, DEL_ME)

// SUBTEST #28 (Perl l.124-125): `validate_vf` synthesises a
// variation_name '1_1_G/C' when VCF ID is `.`.
// Blocked: `synthesize_variation_name(chrom, pos, ref, alt) -> String`
// missing. Future-work entry: `Uploaded_variation fallback name
// synthesis` (future-work-vepyr.md line 454+). Effort: S.
//
// #[test]
// fn port_parser_subtest_28_uploaded_variation_synthesised_from_pos_alleles() {
//     // VCF: 1 1 . G C  → CSQ Uploaded_variation should be "1_1_G/C"
//     //                  (vepyr currently emits literal ".")
//     assert_eq!(
//         synthesize_variation_name("1", 1, "G", "C"),
//         "1_1_G/C".to_string()
//     );
// }

// SUBTEST #29 (Perl l.127-129): per-port `$p->{chr} = ['1']` filter
// excludes rows on chr=2.
// Blocked: per-port `chr_filter: HashSet<String>` plumbing. SQL-level
// `WHERE chrom IN (...)` already works; this is ergonomic CLI sugar.
// Future-work entry: `Per-port chr_filter (CLI --chr analogue)`
// (future-work-vepyr.md line 463+). Effort: S.
//
// #[tokio::test]
// async fn port_parser_subtest_29_chr_filter_includes_excludes_rows() {
//     // chr_filter=["1"] → row on chr=1 stays, row on chr=2 dropped.
//     unimplemented!()
// }

// SUBTESTS #33, #35, #36, #37 — chromosome-synonyms alias table:
// #33 (l.139): `_have_chr(12)` → 1 (vepyr's `normalize_chrom` strips
//   `chr` prefix; it does NOT ADD one). Requires alias lookup.
// #35 (l.142-144): `_have_chr('M')` → 1, with rewrite chr=M → MT.
// #36 (l.146): `chromosome_synonyms($file)` loads a TSV alias table.
// #37 (l.147-149): `_have_chr('NC_000021.9')` resolves via synonyms,
//   chr is left unchanged after the lookup (only matches, never rewrites).
// All blocked on the same vepyr API gap.
// Future-work entry: `PreparedContext chromosome-synonyms alias table`
// (future-work-vepyr.md line 273+, already open from TranscriptTree
// port). Effort: M.
//
// #[test]
// fn port_parser_subtest_33_have_chr_adds_chr_prefix_via_alias() {
//     // valid set keyed "chr12"; lookup of bare "12" must succeed via alias.
//     unimplemented!()
// }
// #[test]
// fn port_parser_subtest_35_have_chr_m_aliases_to_mt() {
//     // Input "M" rewrites to "MT" via M→MT alias; vf.chrom becomes "MT".
//     unimplemented!()
// }
// (and #36/#37)

// SUBTESTS #45-#49 — `check_ref` FASTA input-validation pass:
// #45 (l.189): insertion '-/T' with check_ref=1 passes (insertions
//   bypass check_ref).
// #46 (l.191-192): 'C/T' at chr21:25585733 matches FASTA → ok.
// #47 (l.194-196): 'G/T' at chr21:25585733 mismatches FASTA → false
//   with warning 'Specified reference allele ... does not match'.
// #48 (l.198-199): 'CTT/TCC' multi-base ref matches → ok.
// #49 (l.201-203): 'TTT/T' at chr21:25585733-25585735 — regression for
//   trailing-character handling on multi-base ref.
// Blocked: vepyr `check_ref` input validation missing. Future-work
// entry: `check_ref FASTA input-validation pass` (future-work-vepyr.md
// line 472+). Effort: M (indexed FASTA reader exists at
// annotate_provider.rs:3153; just needs the UDF/table-fn surface).
//
// #[test]
// fn port_parser_subtest_46_check_ref_matches_fasta() {
//     // chr21:25585733 C/T matches reference (real FASTA value at that pos).
//     // assert!(check_ref_filter("21", 25585733, "C", &fasta).unwrap());
//     unimplemented!()
// }
// (and #45, #47, #48, #49 — same pattern, varying assertion polarity)

// SUBTESTS #50-#53 — `lookup_ref` FASTA input-REF backfill:
// #50 (l.212-214): 'N/T' at chr21:25585733 backfills N → C (FASTA);
//   allele_string becomes 'C/T'.
// #51 (l.216-218): 'N/T' on strand -1 backfills with revcomp(C) → G;
//   allele_string becomes 'G/T'.
// #52 (l.220-222): 'N/-' deletion backfills N → C; allele_string 'C/-'.
// #53 (l.224-226): insertion '-/T' lookup_ref is a no-op (no REF base
//   to fill); allele_string stays '-/T'.
// Blocked: vepyr `lookup_ref` input-REF backfill missing. Future-work
// entry: `lookup_ref FASTA input-REF backfill` (future-work-vepyr.md
// line 486+). Effort: M.
//
// #[test]
// fn port_parser_subtest_50_lookup_ref_fills_n_from_fasta() {
//     // lookup_ref_udf("21", 25585733, "N", 1) → "C"
//     unimplemented!()
// }
// (and #51, #52, #53)

// ─────────────────────────────────────────────────────────────────────
// Architectural-no-analogue rows (47 total) — documentation only.
//
// Each row below names the missing-by-design vepyr component and why
// vepyr's data model precludes it. Substantive justifications also live
// in `porting-tests/detailed_plans/Parser.md` §Architectural-no-analogue.
// ─────────────────────────────────────────────────────────────────────

// SUBTESTS #2, #3, #5, #6, #8, #9, #10, #11, #13, #14, #15, #16 —
//   base-`Parser`-class plumbing.
//
// Missing-by-design component: **multi-format `Parser` base class**.
// Perl's `Bio::EnsEMBL::VEP::Parser` is a format dispatcher that
// instantiates one of `Parser::VCF` / `Parser::VEP_input` / `Parser::ID`
// / `Parser::HGVS` / `Parser::Region` / `Parser::SPDI` / `Parser::CAID`
// based on `format` / `detect_format`. Each subclass shares the base's
// `file` / `headers` / `line_number` / `delimiter` / `valid_chromosomes`
// accessors and the gzipped-VCF GLOB-handle trick. Vepyr is VCF-only by
// project-wide design (Parser_CAID/HGVS/ID/Region/SPDI/VEP_input are
// EXCLUDE-no-analogue per port-status.md rows 20-25); `VcfTableProvider`
// IS the single parser, with no dispatch layer. STDIN input is also out
// of scope (CLI takes paths only).

// SUBTESTS #26, #27, #38-#44 — `validate_vf` per-row validation pass
//   with warning emission.
//
// Missing-by-design component: **per-row VF-validation pass with STDERR
// warning emission**. Perl's `validate_vf` runs a check chain (chr in
// valid set; start/end numeric; start ≤ end+1; allele characters; REF
// length vs coord span) and emits one warning per failed check. Vepyr
// handles each implicitly:
// - chr: silent drop via DataFusion partition filter (`discover_vcf_contigs`).
// - numeric POS/END: caught by `noodles_vcf` at parse time; surfaces as
//   a DataFusion `Err`.
// - allele characters (`Q/R` etc): NOT validated. Downstream consequence
//   prediction may produce undefined output. Could be tightened later.
// - REF length vs coord span: NOT validated.
// Adding a `validate_vf` Rust function would require a per-row warning
// sink with no current library consumer; vepyr is a DataFusion-UDF
// library where STDERR warnings don't fit the API surface.

// SUBTEST #30 — `g/c` lowercase allele uppercase coercion.
//
// Missing-by-design component: **input-allele case coercion**. Perl's
// `validate_vf` rewrites lowercase allele strings to uppercase. Vepyr
// does NOT do this — `noodles_vcf` passes alleles through verbatim and
// the consequence engine treats them case-sensitively. Adding this risks
// changing reference-comparison semantics; deliberately punt unless a
// real consumer surfaces a case.

// SUBTESTS #54, #55, #56 — `detect_format` on `<*>` and SV-shape VCF
//   rows (`SVTYPE=DUP;END=...`, `<DUP>`).
// SUBTESTS #58-#73 — `detect_format` on non-VCF formats (rsID, HGVSg,
//   HGVSc, HGVSp, VEP-input, region in 8 strand/separator variants,
//   incomplete-row sentinel).
// SUBTESTS #74, #75 — input as in-memory FH and STDIN-detect-format
//   throws.
// SUBTESTS #76, #77, #78 — `new(format=>'vcf')` / `new(format=>'foo')`
//   / `new(format=>'guess')` factory dispatch.
//
// Missing-by-design component: **format auto-detection and factory
// dispatch**. Vepyr is VCF-only with no detection layer. All non-VCF
// formats (ID, HGVS, Region, VEP-input, SPDI, CAID) are EXCLUDE-no-
// analogue at the project level (port-status.md rows 20-25). Without a
// multi-format parser, there is no "format" to detect or dispatch on.

// SUBTEST #79 — DB-mode block (11 sub-assertions, gated by
//   `SKIP: { ... unless $can_use_db }`).
//
// Missing-by-design component: **MySQL DB mode + contig→toplevel
// transform + LRG mapping**. Vepyr is offline-cache-only by design (no
// MySQL); LRG mapping is haplosaurus territory (also EXCLUDE per
// categorization.md row 19). The Perl DB block exercises:
//   - `validate_vf` on a contig name (`AP000235.3`) that transforms to
//     toplevel coords (chr=21, start=25043768).
//   - LRG round-trip: chr21:43774213 ↔ LRG_485:7166.
// Neither is reachable without a database connection.
