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
            class: Class::parse(&opt("class").ok_or_else(|| ctx("missing `class`".into()))?)
                .map_err(|e| ctx(e))?,
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
        assert!(err.contains("expected_assertions"), "unhelpful error: {err}");
    }
}
