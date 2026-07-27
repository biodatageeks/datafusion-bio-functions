//! The port ledger gate: mechanised "sztywno 1:1" for the Ensembl-VEP → vepyr test port.
//!
//! Every ported Perl test file lives here as a pair:
//!   * `tests/port/perl/<Test>.t`        — the vendored Perl source (Apache-2.0, Ensembl)
//!   * `tests/port/<Test>.ledger.toml`   — one row per Perl assertion, classified
//!
//! This test enumerates the assertions in the vendored Perl and fails the build if any
//! assertion lacks a ledger row, if a `unit-port` row names a Rust test that does not
//! exist, or if a non-ported row lacks its justification. A port that passes proves
//! nothing about faithfulness unless coverage itself is checked — that is what this is.
//!
//! See `tests/port/README.md` for the contract, and
//! `porting-tests/docs/superpowers/specs/2026-07-27-test-port-forward-port-to-master-design.md`
//! for why it exists.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// How a Perl assertion is accounted for on the Rust side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Class {
    /// Has a direct Rust analogue; `rust` names it.
    UnitPort,
    /// Cannot exist in vepyr by design; `reason` justifies it.
    ArchitecturalNoAnalogue,
    /// Should exist eventually; `issue` tracks it.
    BlockedFutureWork,
}

impl Class {
    fn parse(s: &str) -> Result<Self, String> {
        match s {
            "unit-port" => Ok(Class::UnitPort),
            "architectural-no-analogue" => Ok(Class::ArchitecturalNoAnalogue),
            "blocked-future-work" => Ok(Class::BlockedFutureWork),
            other => Err(format!(
                "unknown class {other:?} (expected unit-port | architectural-no-analogue | blocked-future-work)"
            )),
        }
    }
}

/// One ledger row: a Perl assertion and its disposition.
#[derive(Debug, Clone)]
struct Row {
    n: usize,
    line: usize,
    desc: String,
    class: Class,
    /// `<path relative to crate root>::<fn name>`, required for `unit-port`.
    rust: Option<String>,
    /// Required for `architectural-no-analogue`.
    reason: Option<String>,
    /// Required for `blocked-future-work`.
    issue: Option<String>,
}

/// A parsed `<Test>.ledger.toml`.
#[derive(Debug)]
struct Ledger {
    /// Upstream path, e.g. `t/Utils.t` — provenance only.
    perl_file: String,
    /// Vendored file this ledger governs, relative to `tests/port/`.
    vendored: String,
    /// Upstream sha256 at vendoring time. Checked by the periodic upstream-drift
    /// script, NOT here — CI has no ensembl-vep checkout.
    #[allow(dead_code)]
    perl_sha256: String,
    /// Tripwire on the enumerator: the parse must yield exactly this many assertions.
    /// Guards against the scanner silently mis-reading a newly vendored file.
    expected_assertions: usize,
    rows: Vec<Row>,
}

fn parse_ledger(text: &str) -> Result<Ledger, String> {
    let v: toml::Value = text.parse::<toml::Value>().map_err(|e| e.to_string())?;
    let get_str = |k: &str| -> Result<String, String> {
        v.get(k)
            .and_then(toml::Value::as_str)
            .map(str::to_owned)
            .ok_or_else(|| format!("missing or non-string top-level key `{k}`"))
    };
    let expected_assertions = v
        .get("expected_assertions")
        .and_then(toml::Value::as_integer)
        .ok_or_else(|| "missing or non-integer `expected_assertions`".to_string())?;
    let raw_rows = v
        .get("row")
        .and_then(toml::Value::as_array)
        .ok_or_else(|| "missing `[[row]]` array".to_string())?;

    let mut rows = Vec::with_capacity(raw_rows.len());
    for (i, r) in raw_rows.iter().enumerate() {
        let ctx = |m: String| format!("row #{}: {m}", i + 1);
        let int = |k: &str| -> Result<usize, String> {
            r.get(k)
                .and_then(toml::Value::as_integer)
                .map(|n| n as usize)
                .ok_or_else(|| ctx(format!("missing or non-integer `{k}`")))
        };
        let opt = |k: &str| r.get(k).and_then(toml::Value::as_str).map(str::to_owned);
        rows.push(Row {
            n: int("n")?,
            line: int("line")?,
            desc: opt("desc").ok_or_else(|| ctx("missing `desc`".into()))?,
            // `&ctx`, not `ctx`: passing it by value would move a closure that `int`
            // and the `ok_or_else` calls above still borrow.
            class: Class::parse(&opt("class").ok_or_else(|| ctx("missing `class`".into()))?)
                .map_err(&ctx)?,
            rust: opt("rust"),
            reason: opt("reason"),
            issue: opt("issue"),
        });
    }

    Ok(Ledger {
        perl_file: get_str("perl_file")?,
        vendored: get_str("vendored")?,
        perl_sha256: get_str("perl_sha256")?,
        expected_assertions: expected_assertions as usize,
        rows,
    })
}

#[cfg(test)]
mod ledger_parse {
    use super::*;

    const MINIMAL: &str = r#"
perl_file           = "t/Example.t"
vendored            = "perl/Example.t"
perl_sha256         = "0000000000000000000000000000000000000000000000000000000000000000"
expected_assertions = 2

[[row]]
n     = 1
line  = 10
desc  = "does a thing"
class = "unit-port"
rust  = "tests/port_example.rs::does_a_thing"

[[row]]
n      = 2
line   = 11
desc   = "needs a database"
class  = "architectural-no-analogue"
reason = "vepyr is cache-only by design"
"#;

    #[test]
    fn parses_a_minimal_ledger() {
        let l = parse_ledger(MINIMAL).expect("should parse");
        assert_eq!(l.perl_file, "t/Example.t");
        assert_eq!(l.expected_assertions, 2);
        assert_eq!(l.rows.len(), 2);
        assert_eq!(l.rows[0].class, Class::UnitPort);
        assert_eq!(
            l.rows[0].rust.as_deref(),
            Some("tests/port_example.rs::does_a_thing")
        );
        assert_eq!(l.rows[1].class, Class::ArchitecturalNoAnalogue);
    }

    #[test]
    fn rejects_an_unknown_class() {
        let bad = MINIMAL.replace("\"unit-port\"", "\"probably-fine\"");
        let err = parse_ledger(&bad).expect_err("unknown class must be rejected");
        assert!(err.contains("probably-fine"), "unhelpful error: {err}");
    }

    #[test]
    fn rejects_a_missing_expected_assertions() {
        let bad = MINIMAL.replace("expected_assertions = 2", "");
        let err = parse_ledger(&bad).expect_err("must require expected_assertions");
        assert!(
            err.contains("expected_assertions"),
            "unhelpful error: {err}"
        );
    }
}

/// One Test::More / Test::Exception assertion found in a vendored `.t`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Assertion {
    /// 1-based ordinal in source order. This is the ledger's `n`.
    n: usize,
    /// 1-based line on which the assertion statement starts.
    line: usize,
    /// The assertion function, e.g. `is_deeply`.
    func: String,
}

/// Functions that constitute a subtest. `skip`, `plan`, `done_testing`, `diag` and
/// `note` are deliberately absent: they declare or annotate, they do not assert.
const ASSERTION_FNS: &[&str] = &[
    "ok",
    "is",
    "isnt",
    "is_deeply",
    "like",
    "unlike",
    "cmp_ok",
    "isa_ok",
    "can_ok",
    "new_ok",
    "pass",
    "fail",
    "throws_ok",
    "lives_ok",
    "dies_ok",
    "lives_and",
    "warning_like",
    "warnings_like",
    "is_passing",
];

/// Enumerate the assertions in Perl source, in source order.
///
/// Recognition rule: an assertion is a call to one of `ASSERTION_FNS` appearing at a
/// **statement boundary** — i.e. the previous significant character was `;`, `{`, `}`,
/// or nothing (start of file). This is what makes multi-line argument lists safe: the
/// continuation lines of `is_deeply(\n  foo(),\n  bar,\n)` are not at a boundary, so a
/// nested `is(` inside an argument list is never miscounted as a new subtest, while an
/// assertion inside a `SKIP: { ... }` block — which *is* at a boundary — is counted,
/// exactly as the hand numbering does.
///
/// Single-quoted and double-quoted strings, and `#` comments, are skipped when looking
/// for boundary characters, so a `;` inside a message string is not a boundary.
fn enumerate_assertions(src: &str) -> Result<Vec<Assertion>, String> {
    let mut out = Vec::new();
    let mut at_boundary = true;

    for (idx, raw) in src.lines().enumerate() {
        let line_no = idx + 1;
        let bytes = raw.as_bytes();
        let mut i = 0usize;

        while i < bytes.len() {
            let c = bytes[i] as char;

            // A `#` outside a string starts a comment: nothing else on this line matters,
            // and a comment does not change the boundary state.
            if c == '#' {
                break;
            }

            // Skip over quoted strings wholesale, honouring backslash escapes.
            if c == '\'' || c == '"' {
                let quote = c;
                i += 1;
                while i < bytes.len() {
                    let d = bytes[i] as char;
                    if d == '\\' {
                        i += 2;
                        continue;
                    }
                    i += 1;
                    if d == quote {
                        break;
                    }
                }
                at_boundary = false;
                continue;
            }

            if c.is_whitespace() {
                i += 1;
                continue;
            }

            if at_boundary && (c.is_ascii_alphabetic() || c == '_') {
                // Read the identifier and see whether it is an assertion call.
                let start = i;
                while i < bytes.len() {
                    let d = bytes[i] as char;
                    if d.is_ascii_alphanumeric() || d == '_' {
                        i += 1;
                    } else {
                        break;
                    }
                }
                let ident = &raw[start..i];
                // The call must be followed by `(` or a block `{`, allowing whitespace.
                let rest = raw[i..].trim_start();
                let is_call = rest.starts_with('(') || rest.starts_with('{');
                if is_call && ASSERTION_FNS.contains(&ident) {
                    out.push(Assertion {
                        n: out.len() + 1,
                        line: line_no,
                        func: ident.to_owned(),
                    });
                }
                at_boundary = false;
                continue;
            }

            at_boundary = matches!(c, ';' | '{' | '}');
            i += 1;
        }
    }

    if out.is_empty() {
        return Err("no assertions found — the scanner almost certainly mis-read this file".into());
    }
    Ok(out)
}

#[cfg(test)]
mod enumerate {
    use super::*;

    fn lines_of(src: &str) -> Vec<usize> {
        enumerate_assertions(src)
            .unwrap()
            .iter()
            .map(|a| a.line)
            .collect()
    }

    #[test]
    fn counts_single_line_assertions() {
        let src = "is(f(1, 1), '1', 'same');\nis(f(1, 2), '1-2', 'diff');\n";
        assert_eq!(lines_of(src), vec![1, 2]);
    }

    #[test]
    fn treats_a_multi_line_call_as_one_assertion() {
        // A nested `is(` inside the argument list must NOT be counted.
        let src = "is_deeply(\n  numberify(['0', '1']),\n  [0, 1],\n  'arrayref'\n);\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "got {found:?}");
        assert_eq!(found[0].func, "is_deeply");
        assert_eq!(found[0].line, 1);
    }

    #[test]
    fn recognises_a_block_form_assertion() {
        let src = "throws_ok {get_fh()} qr/No file/, 'no file';\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].func, "throws_ok");
    }

    #[test]
    fn is_not_confused_by_a_regex_containing_slashes_and_parens() {
        let src = r"ok(get_time =~ /\d{4}(\-\d\d){2} \d\d(\:\d\d){2}/, 'get_time');";
        assert_eq!(enumerate_assertions(src).unwrap().len(), 1);
    }

    #[test]
    fn counts_assertions_inside_a_skip_block_but_not_skip_itself() {
        let src = "SKIP: {\n  skip 'no gzip', 2 unless $HAVE_GZIP;\n  is(ref($fh), 'GLOB', 'a');\n  is(ref($gh), 'GLOB', 'b');\n}\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 2, "skip must not count: {found:?}");
        assert_eq!(lines_of(src), vec![3, 4]);
    }

    #[test]
    fn ignores_done_testing_and_comments() {
        let src = "# is(1, 1, 'commented out');\nis(1, 1, 'real');\ndone_testing();\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].line, 2);
    }

    #[test]
    fn a_semicolon_inside_a_string_is_not_a_statement_boundary() {
        let src = "my $x = 'a; is(1,1)';\nis(1, 1, 'only me');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "string contents must be inert: {found:?}");
        assert_eq!(found[0].line, 2);
    }

    #[test]
    fn errors_on_a_file_with_no_assertions() {
        assert!(enumerate_assertions("use strict;\nmy $x = 1;\n").is_err());
    }

    /// The decisive cross-check: the machine numbering must equal the numbering a human
    /// used by hand in the pre-existing `port_utils.rs` (rows 1-44, with row 30 = L218,
    /// row 31 = L219, row 32 = L220, row 33 = L221, row 44 = L310).
    #[test]
    fn matches_the_hand_numbering_of_utils_t() {
        let src = include_str!("port/perl/Utils.t");
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(
            found.len(),
            44,
            "expected 44 assertions, got {}",
            found.len()
        );
        let at = |n: usize| found[n - 1].line;
        assert_eq!(
            at(30),
            218,
            "row 30 should be the get_compressed_filehandle ok()"
        );
        assert_eq!(at(31), 219);
        assert_eq!(at(32), 220);
        assert_eq!(at(33), 221);
        assert_eq!(at(44), 310, "row 44 should be the get_version_string ok()");
    }
}

/// Index of Rust test functions declared in the crate, keyed `"<rel path>::<fn>"`.
///
/// Built by scanning source text for a `#[test]` / `#[tokio::test]` attribute followed
/// by a `fn` declaration. This is what makes the gate uncheatable: a comment claiming
/// coverage does not create an entry here.
fn index_rust_tests(roots: &[PathBuf], crate_root: &Path) -> BTreeMap<String, bool> {
    let mut idx = BTreeMap::new();
    let mut stack: Vec<PathBuf> = roots.to_vec();
    while let Some(p) = stack.pop() {
        if p.is_dir() {
            let Ok(rd) = std::fs::read_dir(&p) else {
                continue;
            };
            stack.extend(rd.filter_map(Result::ok).map(|e| e.path()));
            continue;
        }
        if p.extension().and_then(|e| e.to_str()) != Some("rs") {
            continue;
        }
        let Ok(text) = std::fs::read_to_string(&p) else {
            continue;
        };
        let rel = p
            .strip_prefix(crate_root)
            .unwrap_or(&p)
            .to_string_lossy()
            .replace('\\', "/");

        let mut pending_test = false;
        let mut pending_ignore: Option<String> = None;
        for line in text.lines() {
            let t = line.trim();
            if t == "#[test]" || t == "#[tokio::test]" || t.starts_with("#[tokio::test(") {
                pending_test = true;
                continue;
            }
            if t.starts_with("#[ignore") {
                pending_ignore = Some(t.to_owned());
                continue;
            }
            if pending_test {
                if let Some(rest) = t.strip_prefix("fn ").or_else(|| {
                    t.strip_prefix("async fn ")
                        .or_else(|| t.strip_prefix("pub fn "))
                }) {
                    let name: String = rest
                        .chars()
                        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                        .collect();
                    if !name.is_empty() {
                        // value = "is ignored with a reason"
                        let ok_ignore = pending_ignore
                            .as_deref()
                            .map(|a| a.contains('=') && a.contains('#'))
                            .unwrap_or(true);
                        idx.insert(format!("{rel}::{name}"), ok_ignore);
                    }
                    pending_test = false;
                    pending_ignore = None;
                } else if !t.is_empty() && !t.starts_with("//") && !t.starts_with("#[") {
                    pending_test = false;
                    pending_ignore = None;
                }
            }
        }
    }
    idx
}

/// All the ways a ledger can fail its contract. Returns one message per violation so a
/// contributor sees every problem at once instead of fixing them one build at a time.
fn check_ledger(
    ledger: &Ledger,
    assertions: &[Assertion],
    rust_tests: &BTreeMap<String, bool>,
) -> Vec<String> {
    let mut errs = Vec::new();
    let name = &ledger.perl_file;

    if assertions.len() != ledger.expected_assertions {
        errs.push(format!(
            "{name}: enumerator found {} assertions but the ledger declares expected_assertions = {}. \
             Either the vendored file changed (re-audit) or the scanner mis-read it (fix the scanner) \
             — do NOT simply update the number.",
            assertions.len(),
            ledger.expected_assertions
        ));
    }

    let by_n: BTreeMap<usize, &Row> = ledger.rows.iter().map(|r| (r.n, r)).collect();
    for a in assertions {
        match by_n.get(&a.n) {
            None => errs.push(format!(
                "{name}: Perl assertion #{} ({} at line {}) has NO ledger row — \
                 sztywno 1:1 violated. Add a [[row]] classifying it.",
                a.n, a.func, a.line
            )),
            Some(r) if r.line != a.line => errs.push(format!(
                "{name}: row n={} says line {} but assertion #{} is at line {} — \
                 the ledger is out of sync with the vendored Perl.",
                r.n, r.line, a.n, a.line
            )),
            Some(_) => {}
        }
    }
    for r in &ledger.rows {
        if r.n == 0 || r.n > assertions.len() {
            errs.push(format!(
                "{name}: row n={} does not correspond to any assertion (file has {})",
                r.n,
                assertions.len()
            ));
        }
        match r.class {
            Class::UnitPort => match r.rust.as_deref() {
                None => errs.push(format!(
                    "{name}: row n={} is unit-port but has no `rust` field",
                    r.n
                )),
                Some(target) if !rust_tests.contains_key(target) => errs.push(format!(
                    "{name}: row n={} names Rust test `{target}`, which does not exist. \
                     A claim of coverage must resolve to a real #[test] fn.",
                    r.n
                )),
                Some(target) if rust_tests.get(target) == Some(&false) => errs.push(format!(
                    "{name}: row n={} names `{target}`, which is #[ignore]d without a reason. \
                     Use #[ignore = \"vepyr#NN: why\"].",
                    r.n
                )),
                Some(_) => {}
            },
            Class::ArchitecturalNoAnalogue => {
                if r.reason.as_deref().unwrap_or("").trim().is_empty() {
                    errs.push(format!(
                        "{name}: row n={} is architectural-no-analogue but has no `reason`. \
                         An unjustified claim of impossibility is a hidden gap.",
                        r.n
                    ));
                }
            }
            Class::BlockedFutureWork => {
                let issue = r.issue.as_deref().unwrap_or("");
                if !issue.contains('#') {
                    errs.push(format!(
                        "{name}: row n={} is blocked-future-work but `issue` ({issue:?}) \
                         is not an issue reference — deferred work must be tracked.",
                        r.n
                    ));
                }
            }
        }
        if r.desc.trim().is_empty() {
            errs.push(format!("{name}: row n={} has an empty `desc`", r.n));
        }
    }
    errs
}

#[cfg(test)]
mod conditions {
    use super::*;

    fn assertions(n: usize) -> Vec<Assertion> {
        (1..=n)
            .map(|i| Assertion {
                n: i,
                line: 9 + i,
                func: "is".into(),
            })
            .collect()
    }

    fn ledger_with(rows: &str, expected: usize) -> Ledger {
        parse_ledger(&format!(
            "perl_file = \"t/E.t\"\nvendored = \"perl/E.t\"\n\
             perl_sha256 = \"00\"\nexpected_assertions = {expected}\n{rows}"
        ))
        .expect("fixture must parse")
    }

    fn tests_with(name: &str) -> BTreeMap<String, bool> {
        let mut m = BTreeMap::new();
        m.insert(name.to_owned(), true);
        m
    }

    #[test]
    fn flags_an_assertion_with_no_row() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"unit-port\"\nrust=\"tests/p.rs::a\"\n",
            2,
        );
        let errs = check_ledger(&l, &assertions(2), &tests_with("tests/p.rs::a"));
        assert!(
            errs.iter()
                .any(|e| e.contains("#2") && e.contains("NO ledger row")),
            "{errs:?}"
        );
    }

    #[test]
    fn flags_a_unit_port_naming_a_nonexistent_rust_test() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"unit-port\"\nrust=\"tests/p.rs::ghost\"\n",
            1,
        );
        let errs = check_ledger(&l, &assertions(1), &BTreeMap::new());
        assert!(
            errs.iter().any(|e| e.contains("does not exist")),
            "{errs:?}"
        );
    }

    #[test]
    fn flags_blocked_future_work_without_an_issue() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"blocked-future-work\"\n",
            1,
        );
        let errs = check_ledger(&l, &assertions(1), &BTreeMap::new());
        assert!(
            errs.iter().any(|e| e.contains("not an issue reference")),
            "{errs:?}"
        );
    }

    #[test]
    fn flags_architectural_no_analogue_without_a_reason() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"architectural-no-analogue\"\n",
            1,
        );
        let errs = check_ledger(&l, &assertions(1), &BTreeMap::new());
        assert!(errs.iter().any(|e| e.contains("no `reason`")), "{errs:?}");
    }

    #[test]
    fn flags_an_expected_assertions_mismatch() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"architectural-no-analogue\"\nreason=\"x\"\n",
            7,
        );
        let errs = check_ledger(&l, &assertions(1), &BTreeMap::new());
        assert!(
            errs.iter().any(|e| e.contains("expected_assertions")),
            "{errs:?}"
        );
    }

    #[test]
    fn flags_a_line_drift_between_ledger_and_perl() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=999\ndesc=\"a\"\nclass=\"architectural-no-analogue\"\nreason=\"x\"\n",
            1,
        );
        let errs = check_ledger(&l, &assertions(1), &BTreeMap::new());
        assert!(errs.iter().any(|e| e.contains("out of sync")), "{errs:?}");
    }

    #[test]
    fn accepts_a_complete_ledger() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"unit-port\"\nrust=\"tests/p.rs::a\"\n\
             [[row]]\nn=2\nline=11\ndesc=\"b\"\nclass=\"blocked-future-work\"\nissue=\"biodatageeks/vepyr#42\"\n",
            2,
        );
        let errs = check_ledger(&l, &assertions(2), &tests_with("tests/p.rs::a"));
        assert!(errs.is_empty(), "should be clean, got {errs:?}");
    }
}

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The gate. Walks every `tests/port/*.ledger.toml`, enumerates its vendored Perl, and
/// asserts the contract. Also asserts the reverse direction: every `tests/port_*.rs`
/// file must be referenced by at least one ledger row, so a port cannot land by simply
/// not writing a ledger.
#[test]
fn every_port_ledger_is_complete() {
    let root = crate_root();
    let port_dir = root.join("tests/port");
    let rust_tests = index_rust_tests(&[root.join("src"), root.join("tests")], &root);

    let mut ledgers = Vec::new();
    for entry in std::fs::read_dir(&port_dir).expect("tests/port must exist") {
        let p = entry.expect("readable dir entry").path();
        if p.to_string_lossy().ends_with(".ledger.toml") {
            ledgers.push(p);
        }
    }
    ledgers.sort();
    assert!(
        !ledgers.is_empty(),
        "no ledgers found in {} — the gate would be vacuous",
        port_dir.display()
    );

    let mut errs: Vec<String> = Vec::new();
    let mut referenced: Vec<String> = Vec::new();

    for lp in &ledgers {
        let text = std::fs::read_to_string(lp).expect("readable ledger");
        let ledger = match parse_ledger(&text) {
            Ok(l) => l,
            Err(e) => {
                errs.push(format!("{}: {e}", lp.display()));
                continue;
            }
        };
        let perl_path = port_dir.join(&ledger.vendored);
        let src = match std::fs::read_to_string(&perl_path) {
            Ok(s) => s,
            Err(e) => {
                errs.push(format!(
                    "{}: cannot read vendored Perl {}: {e}",
                    lp.display(),
                    perl_path.display()
                ));
                continue;
            }
        };
        match enumerate_assertions(&src) {
            Ok(assertions) => {
                errs.extend(check_ledger(&ledger, &assertions, &rust_tests));
                referenced.extend(
                    ledger
                        .rows
                        .iter()
                        .filter_map(|r| r.rust.clone())
                        .filter_map(|t| t.split("::").next().map(str::to_owned)),
                );
            }
            Err(e) => errs.push(format!("{}: {e}", perl_path.display())),
        }
    }

    // Reverse direction: no orphan port test files.
    for entry in std::fs::read_dir(root.join("tests")).expect("tests/ must exist") {
        let p = entry.expect("readable dir entry").path();
        let Some(fname) = p.file_name().and_then(|f| f.to_str()) else {
            continue;
        };
        if !fname.starts_with("port_") || !fname.ends_with(".rs") || fname == "port_ledger.rs" {
            continue;
        }
        let rel = format!("tests/{fname}");
        if !referenced.contains(&rel) {
            errs.push(format!(
                "{rel} is a port test file but no ledger row references it — \
                 every ported test must be accounted for by a ledger."
            ));
        }
    }

    assert!(
        errs.is_empty(),
        "port ledger gate failed with {} problem(s):\n  - {}",
        errs.len(),
        errs.join("\n  - ")
    );
}
