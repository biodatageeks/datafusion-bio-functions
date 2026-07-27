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

/// Functions that constitute a subtest.
///
/// **Derivation** (do not extend this list by guessing). Every `use`d testing module in
/// `ensembl-vep/t/*.t` was enumerated; the set is exactly four:
///
/// | module          | files | assertion-emitting exports kept here                                                                          |
/// |-----------------|------:|---------------------------------------------------------------------------------------------------------------|
/// | `Test::More`    |    49 | `ok is isnt like unlike is_deeply cmp_ok isa_ok can_ok new_ok pass fail use_ok require_ok`                     |
/// | `Test::Exception`|   48 | `throws_ok dies_ok lives_ok lives_and`                                                                         |
/// | `Test::Deep`    |     3 | `cmp_deeply cmp_bag cmp_set cmp_methods`                                                                       |
/// | `Test::Warnings`|     3 | `warning warnings` (imported as `qw(warning :no_end_test)`)                                                    |
///
/// Corpus-attested (call-position occurrences over all 49 files): `is` 842,
/// `is_deeply` 588, `ok` 270, `use_ok` 148, `throws_ok` 101, `like` 11, `warning` 6,
/// `cmp_deeply` 3, `dies_ok` 2, `isa_ok` 1. The remaining names never occur, but are
/// kept because their module *is* imported corpus-wide, so the next vendored file may
/// use them — the failure mode this gate exists to prevent is a missing name.
///
/// Deliberately absent because they declare or annotate rather than assert: `plan`,
/// `done_testing`, `skip`, `todo`, `todo_skip`, `BAIL_OUT`, `diag`, `note`, `explain`,
/// and `subtest` — a `subtest` name is not a result, while assertions nested inside its
/// block are counted individually, since a block `{` opens a statement boundary.
/// Also absent: `eq_deeply`/`eq_array`/`eq_hash`/`eq_set` (bool predicates, no TAP
/// output), and `Test::Warn`'s `warning_like`/`warnings_like` plus `Test::Builder`'s
/// `is_passing` — **no** corpus file imports those, so keeping them would risk
/// mis-reading a same-named local sub as an assertion.
///
/// `warning`/`warnings` count only in *bare statement* position: `warning { ... };`
/// asserts that the block warned, whereas `my $w = warning { ... };` is a capture whose
/// assertion comes later (`like($w, qr/.../)`), and is excluded for free because `my`
/// consumes the statement boundary.
const ASSERTION_FNS: &[&str] = &[
    // Test::More
    "ok",
    "is",
    "isnt",
    "like",
    "unlike",
    "is_deeply",
    "cmp_ok",
    "isa_ok",
    "can_ok",
    "new_ok",
    "pass",
    "fail",
    "use_ok",
    "require_ok",
    // Test::Exception
    "throws_ok",
    "dies_ok",
    "lives_ok",
    "lives_and",
    // Test::Deep
    "cmp_deeply",
    "cmp_bag",
    "cmp_set",
    "cmp_methods",
    // Test::Warnings
    "warning",
    "warnings",
];

/// Declaration-only Test::More functions. **Not** assertions; recorded separately purely
/// as a parse-health landmark — see [`Scan::declarations`].
const DECLARATION_FNS: &[&str] = &[
    "done_testing",
    "plan",
    "subtest",
    "skip",
    "diag",
    "note",
    "todo",
    "todo_skip",
    "BAIL_OUT",
    "explain",
];

/// What a full scan of one `.t` file yielded.
#[derive(Debug, Clone)]
struct Scan {
    /// The assertions, in source order, numbered from 1. This is the ledger's `n`.
    assertions: Vec<Assertion>,
    /// Statement-boundary calls to [`DECLARATION_FNS`]. Not assertions. Their value is as
    /// a **parse-health landmark**: every upstream `.t` ends in `done_testing();`, so if
    /// the scanner fails to see that call at a statement boundary its boundary state
    /// desynchronised somewhere — an unterminated string or regex that ran to EOF, say.
    /// The corpus sweep turns that into a loud failure instead of a quiet under-count.
    declarations: Vec<Assertion>,
}

impl Scan {
    /// Boundary-recognised calls to a particular declaration function.
    fn declarations_named(&self, func: &str) -> usize {
        self.declarations.iter().filter(|d| d.func == func).count()
    }
}

/// Quote-like operators and how many delimited parts each takes.
///
/// Longest name first: `qq`/`qw`/`qr` must win over `q`, and `tr` over a bare `r`.
const QUOTE_OPS: &[(&str, u8)] = &[
    ("qq", 1),
    ("qw", 1),
    ("qr", 1),
    ("tr", 2),
    ("q", 1),
    ("m", 1),
    ("s", 2),
    ("y", 2),
];

/// Quote-like operators whose names are also plausible barewords or variable names (hash
/// keys, methods, regex escapes such as `\s`). These are honoured only where Perl syntax
/// actually expects a term, which keeps `{foo => 'bar', s => 5}` (9× in `TranscriptTree.t`)
/// and `->tr(...)` out of the lexer.
///
/// `qq`, `qw` and `qr` are absent because they must work where no term is expected:
/// `use FindBin qw($Bin)` follows a bareword, and `throws_ok {…} qr/…/` follows a `}`.
const AMBIGUOUS_QUOTE_OPS: &[&str] = &["q", "m", "s", "y", "tr"];

/// A Perl lexer just deep enough to locate statement boundaries reliably.
///
/// The recognition rule is unchanged and is the whole point: an assertion is a call to one
/// of [`ASSERTION_FNS`] at a **statement boundary** — the previous significant character
/// was `;`, `{`, `}`, or nothing. That is what makes multi-line argument lists safe: the
/// continuation lines of `is_deeply(\n  foo(),\n  bar,\n)` are not at a boundary, so a
/// nested `is(` inside an argument list is never miscounted, while an assertion inside a
/// `SKIP: { ... }` block — which *is* at a boundary — is counted, exactly as a hand
/// numbering does.
///
/// Everything else in here exists to stop the boundary state from being corrupted. The
/// scanner runs over the whole source rather than line-by-line, because every construct
/// that can hide a `;` — quoted strings, `q()`/`qq()`/`qw()`/`qr//`/`m//`/`s///`/`tr///`,
/// heredocs — can span lines.
struct PerlScan<'a> {
    src: &'a str,
    b: &'a [u8],
    i: usize,
    line: usize,
    /// Previous significant character was `;`, `{`, `}`, or start of file.
    at_boundary: bool,
    /// Perl syntax expects a *term* here, so `/` opens a regex rather than dividing.
    expect_term: bool,
    /// Heredoc terminators opened on the current line, with their `<<~` indent flag.
    /// Bodies begin after the line's newline, in declaration order.
    heredocs: Vec<(String, bool)>,
    assertions: Vec<Assertion>,
    declarations: Vec<Assertion>,
}

impl<'a> PerlScan<'a> {
    fn new(src: &'a str) -> Self {
        Self {
            src,
            b: src.as_bytes(),
            i: 0,
            line: 1,
            at_boundary: true,
            expect_term: true,
            heredocs: Vec::new(),
            assertions: Vec::new(),
            declarations: Vec::new(),
        }
    }

    fn at(&self, off: usize) -> Option<u8> {
        self.b.get(self.i + off).copied()
    }

    /// Advance one byte, keeping `line` in step.
    fn bump(&mut self) {
        if self.b.get(self.i) == Some(&b'\n') {
            self.line += 1;
        }
        self.i += 1;
    }

    /// Consumed a value: no longer at a statement boundary, and an operator comes next.
    fn took_term(&mut self) {
        self.at_boundary = false;
        self.expect_term = false;
    }

    fn run(mut self) -> Result<Scan, String> {
        while self.i < self.b.len() {
            let c = self.b[self.i];

            if c == b'\n' {
                self.bump();
                self.consume_heredoc_bodies()?;
                continue;
            }
            if c.is_ascii_whitespace() {
                self.bump();
                continue;
            }

            // A variable, consumed sigils-and-name together. This must happen *before* the
            // identifier branch: `$tr->get_all_Attributes` otherwise reads `tr` as a
            // substitution operator with `-` for a delimiter and eats the next 90 lines.
            if matches!(c, b'$' | b'@' | b'%' | b'&') && self.try_variable() {
                self.took_term();
                continue;
            }

            // `#` is a comment *unless* it is Perl's last-index sigil (`$#array`,
            // `$#{$ref}`, `$#$ref`). Treating it as a comment unconditionally truncates
            // the rest of the line, hiding any `;` on it — and a hidden `;` is a lost
            // statement boundary, which silently swallows the next assertion.
            if c == b'#' {
                if self.is_last_index_sigil() {
                    self.bump();
                    self.took_term();
                    continue;
                }
                while self.i < self.b.len() && self.b[self.i] != b'\n' {
                    self.i += 1; // never a newline, so `line` cannot drift
                }
                continue; // a comment does not change the boundary state
            }

            // Quoted strings and backtick commands. Multi-line by nature.
            if matches!(c, b'\'' | b'"' | b'`') {
                self.bump();
                self.consume_delimited(c, c, false)?;
                self.took_term();
                continue;
            }

            // Two-character operators worth spelling out: `=>` expects a term after it,
            // `->` does not (so `$obj->s` can never be read as a substitution).
            if c == b'=' && self.at(1) == Some(b'>') {
                self.bump();
                self.bump();
                self.at_boundary = false;
                self.expect_term = true;
                continue;
            }
            if c == b'-' && self.at(1) == Some(b'>') {
                self.bump();
                self.bump();
                self.took_term();
                continue;
            }

            if c == b'<' && self.at(1) == Some(b'<') && self.try_heredoc() {
                self.took_term();
                continue;
            }

            if is_ident_start(c) {
                let start = self.i;
                while self.i < self.b.len() && is_ident_char(self.b[self.i]) {
                    self.i += 1; // identifiers cannot contain newlines
                }
                let ident = &self.src[start..self.i];
                let line = self.line;

                if let Some(parts) = self.quote_op_here(ident) {
                    self.consume_quote_like(parts)?;
                    self.took_term();
                    continue;
                }

                if self.at_boundary && self.call_follows() {
                    if ASSERTION_FNS.contains(&ident) {
                        self.assertions.push(Assertion {
                            n: self.assertions.len() + 1,
                            line,
                            func: ident.to_owned(),
                        });
                    } else if DECLARATION_FNS.contains(&ident) {
                        self.declarations.push(Assertion {
                            n: self.declarations.len() + 1,
                            line,
                            func: ident.to_owned(),
                        });
                    }
                }
                self.at_boundary = false;
                self.expect_term = keyword_expects_term(ident);
                continue;
            }

            // A bare `/.../` regex, but only where a term is syntactically expected —
            // otherwise `/` is division. Getting this wrong in the permissive direction
            // would swallow a `;`, so the predicate is deliberately narrow.
            if c == b'/' && self.expect_term {
                self.bump();
                self.consume_delimited(b'/', b'/', false)?;
                self.consume_modifiers();
                self.took_term();
                continue;
            }

            self.at_boundary = matches!(c, b';' | b'{' | b'}');
            self.expect_term = punct_expects_term(c);
            self.bump();
        }

        if let Some((term, _)) = self.heredocs.first() {
            return Err(format!(
                "unterminated heredoc <<{term} at EOF — the scanner cannot parse this file confidently"
            ));
        }
        if self.assertions.is_empty() {
            return Err(
                "no assertions found — the scanner almost certainly mis-read this file".into(),
            );
        }
        Ok(Scan {
            assertions: self.assertions,
            declarations: self.declarations,
        })
    }

    /// Consume a Perl variable — `$x`, `@list`, `%h`, `$$ref`, `$#array`, `$Foo::Bar` —
    /// leaving `self.i` past its name. Returns false when the sigil is really an operator
    /// (`%` modulus, `&` bitwise-and, `$` before `{`), in which case `self.i` is untouched
    /// and the caller falls through to the punctuation path.
    fn try_variable(&mut self) -> bool {
        let mut j = self.i;
        while matches!(
            self.b.get(j),
            Some(b'$') | Some(b'@') | Some(b'%') | Some(b'&')
        ) {
            j += 1;
        }
        if self.b.get(j) == Some(&b'#') {
            j += 1; // `$#array`, `$#$ref`
        }
        if !self.b.get(j).copied().is_some_and(is_ident_start) {
            return false;
        }
        while self.b.get(j).copied().is_some_and(is_ident_char)
            || (self.b.get(j) == Some(&b':')
                && self.b.get(j + 1) == Some(&b':')
                && self.b.get(j + 2).copied().is_some_and(is_ident_char))
        {
            j += if self.b[j] == b':' { 2 } else { 1 };
        }
        self.i = j; // no newline can occur inside a variable name
        true
    }

    /// True when the `#` at `self.i` is Perl's `$#` last-index sigil rather than a comment.
    fn is_last_index_sigil(&self) -> bool {
        match self.i {
            0 => false,
            // `$#array`, `$#$ref`
            _ if self.b[self.i - 1] == b'$' => true,
            // `${#...}`
            n if n >= 2 => self.b[n - 1] == b'{' && self.b[n - 2] == b'$',
            _ => false,
        }
    }

    /// Does an assertion *call* follow the identifier we just read?
    ///
    /// Whitespace, newlines and comments may intervene: `FilterSet.t` writes
    /// `throws_ok\n  { ... }\n  qr/.../,\n  'desc';` seven times, and requiring the `(`
    /// or `{` on the same line lost all seven. Nothing but whitespace and comments is
    /// skipped, so the lookahead cannot run past the end of the statement.
    fn call_follows(&self) -> bool {
        let mut j = self.i;
        loop {
            match self.b.get(j) {
                Some(c) if c.is_ascii_whitespace() => j += 1,
                Some(b'#') => {
                    while j < self.b.len() && self.b[j] != b'\n' {
                        j += 1;
                    }
                }
                Some(b'(') | Some(b'{') => return true,
                _ => return false,
            }
        }
    }

    /// If `ident` is a quote-like operator *in this position*, return its part count.
    fn quote_op_here(&mut self, ident: &str) -> Option<u8> {
        let (_, parts) = QUOTE_OPS.iter().find(|(name, _)| *name == ident)?;
        if AMBIGUOUS_QUOTE_OPS.contains(&ident) && !self.expect_term {
            return None;
        }
        // Perl allows whitespace between the operator and its delimiter (except before
        // `#`); only spaces and tabs are accepted here, which is what the corpus uses.
        let mut j = self.i;
        while matches!(self.b.get(j), Some(b' ') | Some(b'\t')) {
            j += 1;
        }
        if !is_quote_delimiter(*self.b.get(j)?) {
            return None;
        }
        self.i = j;
        Some(*parts)
    }

    /// Consume a quote-like body, `self.i` on the opening delimiter.
    ///
    /// Handles both paired delimiters (`() {} [] <>`, which nest) and the arbitrary
    /// single-character form, plus `s///` / `tr///` / `y///`'s two-part shape: with a
    /// paired delimiter the second part carries its own (`s{a}{b}`), otherwise the closing
    /// delimiter of the first part doubles as the second's opener (`s/a/b/`).
    fn consume_quote_like(&mut self, parts: u8) -> Result<(), String> {
        let open = self.b[self.i];
        let close = paired_close(open);
        let nestable = close != open;
        self.bump();
        self.consume_delimited(open, close, nestable)?;

        if parts == 2 {
            if nestable {
                self.skip_gap();
                let open2 = *self
                    .b
                    .get(self.i)
                    .ok_or_else(|| format!("line {}: truncated s///-style operator", self.line))?;
                let close2 = paired_close(open2);
                self.bump();
                self.consume_delimited(open2, close2, close2 != open2)?;
            } else {
                // The first part's terminator opened the second.
                self.consume_delimited(open, close, false)?;
            }
        }
        self.consume_modifiers();
        Ok(())
    }

    /// Consume up to and including the closing delimiter, `self.i` just past the opener.
    fn consume_delimited(&mut self, open: u8, close: u8, nestable: bool) -> Result<(), String> {
        let start_line = self.line;
        let mut depth = 1usize;
        while self.i < self.b.len() {
            let c = self.b[self.i];
            if c == b'\\' {
                self.bump();
                if self.i < self.b.len() {
                    self.bump();
                }
                continue;
            }
            if nestable && c == open {
                depth += 1;
                self.bump();
                continue;
            }
            if c == close {
                self.bump();
                depth -= 1;
                if depth == 0 {
                    return Ok(());
                }
                continue;
            }
            self.bump();
        }
        Err(format!(
            "line {start_line}: unterminated {}…{} construct — the scanner cannot parse this file confidently",
            open as char, close as char
        ))
    }

    fn consume_modifiers(&mut self) {
        while self.b.get(self.i).is_some_and(|c| c.is_ascii_alphabetic()) {
            self.i += 1;
        }
    }

    /// Skip whitespace and comments between the two halves of `s{a}{b}`.
    fn skip_gap(&mut self) {
        loop {
            match self.b.get(self.i) {
                Some(c) if c.is_ascii_whitespace() => self.bump(),
                Some(b'#') => {
                    while self.i < self.b.len() && self.b[self.i] != b'\n' {
                        self.i += 1;
                    }
                }
                _ => return,
            }
        }
    }

    /// Recognise `<<TERM`, `<<"TERM"`, `<<'TERM'`, `<<~TERM` and queue the body. Returns
    /// false for a left-shift, leaving `self.i` untouched.
    fn try_heredoc(&mut self) -> bool {
        let mut j = self.i + 2;
        let indented = self.b.get(j) == Some(&b'~');
        if indented {
            j += 1;
        }
        let (term, end) = match self.b.get(j) {
            Some(&q @ (b'\'' | b'"')) => {
                let start = j + 1;
                let mut k = start;
                while k < self.b.len() && self.b[k] != q && self.b[k] != b'\n' {
                    k += 1;
                }
                if self.b.get(k) != Some(&q) {
                    return false;
                }
                (self.src[start..k].to_owned(), k + 1)
            }
            Some(&c) if is_ident_start(c) => {
                let start = j;
                let mut k = j;
                while k < self.b.len() && is_ident_char(self.b[k]) {
                    k += 1;
                }
                (self.src[start..k].to_owned(), k)
            }
            // `1 << 3`, `$x << $y`: a left-shift, not a heredoc.
            _ => return false,
        };
        self.heredocs.push((term, indented));
        self.i = end;
        true
    }

    /// Consume the bodies of every heredoc opened on the line just ended. A no-op — and
    /// in particular *not* a change of boundary state — when none were opened.
    fn consume_heredoc_bodies(&mut self) -> Result<(), String> {
        if self.heredocs.is_empty() {
            return Ok(());
        }
        let pending = std::mem::take(&mut self.heredocs);
        for (term, indented) in pending {
            loop {
                if self.i >= self.b.len() {
                    return Err(format!(
                        "unterminated heredoc <<{term} at EOF — the scanner cannot parse this file confidently"
                    ));
                }
                let start = self.i;
                while self.i < self.b.len() && self.b[self.i] != b'\n' {
                    self.i += 1;
                }
                let content = &self.src[start..self.i];
                if self.i < self.b.len() {
                    self.bump(); // the newline
                }
                let matches_term = if indented {
                    content.trim() == term
                } else {
                    content == term
                };
                if matches_term {
                    break;
                }
            }
        }
        // The heredoc *term* was already consumed at the `<<EOF` marker; the body is not a
        // token of its own, so the boundary state at the end of the opening line — which
        // may well be the `;` that terminated it — must survive intact.
        Ok(())
    }
}

fn is_ident_start(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}

fn is_ident_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_'
}

/// Delimiters accepted after a quote-like operator.
///
/// `=` and `>` are excluded so a fat comma (`s => 5`, seen 9× in `TranscriptTree.t`) is
/// never read as a substitution; `-` is excluded because `->` is pervasive; `,` `;` `)`
/// `]` `}` are excluded because they close a term rather than open one (`$h->{s}`); `#` is
/// excluded because Perl's `q#…#` form is indistinguishable here from a trailing comment.
/// Nothing in the corpus uses any of the excluded characters as a real delimiter.
fn is_quote_delimiter(c: u8) -> bool {
    c.is_ascii_punctuation()
        && !matches!(
            c,
            b'=' | b'>' | b'-' | b',' | b';' | b')' | b']' | b'}' | b'#'
        )
}

fn paired_close(open: u8) -> u8 {
    match open {
        b'(' => b')',
        b'[' => b']',
        b'{' => b'}',
        b'<' => b'>',
        other => other,
    }
}

/// After which punctuation does Perl expect a term (so `/` would open a regex)?
///
/// `)`, `]` and `}` are absent on purpose: they *close* a term, so `$h->{x} / 2` stays a
/// division. Being wrong in the permissive direction is the dangerous one — it would
/// consume a `;` into a phantom regex and lose a statement boundary.
fn punct_expects_term(c: u8) -> bool {
    matches!(
        c,
        b'(' | b'{' | b'[' | b',' | b';' | b'=' | b'~' | b'!' | b'&' | b'|' | b'?' | b':' | b'.'
    )
}

/// Named operators and keywords after which a term is expected, so that `eval q{ … }`
/// (17× in the corpus) and `$x =~ m/…/` lex correctly.
fn keyword_expects_term(ident: &str) -> bool {
    matches!(
        ident,
        "eval"
            | "return"
            | "print"
            | "printf"
            | "push"
            | "unshift"
            | "join"
            | "split"
            | "grep"
            | "map"
            | "sort"
            | "if"
            | "elsif"
            | "unless"
            | "while"
            | "until"
            | "and"
            | "or"
            | "not"
            | "xor"
            | "defined"
            | "ref"
            | "scalar"
            | "die"
            | "warn"
            | "no"
            | "use"
            | "my"
            | "our"
            | "local"
            | "x"
            | "lc"
            | "uc"
    )
}

/// Scan Perl source, returning assertions and parse-health landmarks.
fn scan_perl(src: &str) -> Result<Scan, String> {
    PerlScan::new(src).run()
}

/// Enumerate the assertions in Perl source, in source order.
fn enumerate_assertions(src: &str) -> Result<Vec<Assertion>, String> {
    scan_perl(src).map(|s| s.assertions)
}

/// A deliberately dumb second opinion: count lines whose *first* token is an assertion
/// function name followed by `(`, `{`, or nothing.
///
/// It shares only the name list with [`scan_perl`] — none of its boundary, comment,
/// string or quote-like machinery. Agreement between the two is therefore corroboration
/// rather than tautology, which is what makes the corpus sweep worth trusting. Its own
/// blind spots are the mirror image: it cannot see an assertion that is not the first
/// token on its line, and it would happily count one inside a multi-line string.
fn count_line_anchored(src: &str) -> usize {
    src.lines()
        .filter(|line| {
            let t = line.trim_start();
            let end = t
                .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
                .unwrap_or(t.len());
            let (name, rest) = t.split_at(end);
            if !ASSERTION_FNS.contains(&name) {
                return false;
            }
            let rest = rest.trim_start();
            rest.is_empty() || rest.starts_with('(') || rest.starts_with('{')
        })
        .count()
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

    // ---- defect 1: the assertion-function set ------------------------------------

    #[test]
    fn counts_use_ok_which_is_148_of_the_corpus() {
        let src = "use_ok('Bio::EnsEMBL::VEP::Runner');\nis(1, 1, 'x');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 2, "{found:?}");
        assert_eq!(found[0].func, "use_ok");
    }

    #[test]
    fn counts_cmp_deeply_from_test_deep() {
        let src = "cmp_deeply(\n  $got,\n  bag(1, 2),\n  'unordered'\n);\n";
        assert_eq!(enumerate_assertions(src).unwrap()[0].func, "cmp_deeply");
    }

    #[test]
    fn counts_a_bare_warning_block_but_not_a_captured_one() {
        // `my $w = warning {...}` is a capture whose assertion is the following `like`.
        let src = "warning { $as->annotate($ib) };\nmy $w = warning { $r->run() };\nlike($w, qr/nope/, 'warned');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(
            found.iter().map(|a| a.func.as_str()).collect::<Vec<_>>(),
            vec!["warning", "like"],
            "{found:?}"
        );
    }

    #[test]
    fn does_not_count_declarations() {
        let src = "plan tests => 2;\ndiag('hi');\nnote('there');\nsubtest 'group' => sub {\n  is(1, 1, 'inner');\n};\ndone_testing();\n";
        let scan = scan_perl(src).unwrap();
        assert_eq!(scan.assertions.len(), 1, "{:?}", scan.assertions);
        assert_eq!(scan.assertions[0].func, "is");
    }

    // ---- defect 2: quote-like operators ------------------------------------------

    #[test]
    fn an_apostrophe_inside_qr_does_not_open_a_phantom_string() {
        // Haplo_Runner.t:497 verbatim in shape. Before quote-like operators were modelled
        // the `'` in `can't` opened a string that ran off the end of the line, carrying
        // `at_boundary == false` into the next line and swallowing the assertion there.
        let src = "like($w, qr/Haplosaurus can't find transcripts/, 'warning message');\nis(1, 1, 'next line survives');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 2, "the next line was swallowed: {found:?}");
        assert_eq!(found[1].line, 2);
    }

    #[test]
    fn models_every_quote_like_operator() {
        for hidden in [
            "q(; is(1,1))",
            "qq{; is(1,1)}",
            "qw(; is(1,1))",
            "qr/; is(1,1)/",
            "my $x = 'a'; $x =~ m{; is(1,1)}",
            "my $x = 'a'; $x =~ s/; is(1,1)/x/",
            "my $x = 'a'; $x =~ s{; is(1,1)}{x}",
            "my $x = 'a'; $x =~ tr/;/x/",
            "my $x = 'a'; $x =~ y/;/x/",
        ] {
            let src = format!("my $v = {hidden};\nok(1, 'the only assertion');\n");
            let found = enumerate_assertions(&src).unwrap();
            assert_eq!(
                found.len(),
                1,
                "`{hidden}` leaked an assertion or a boundary: {found:?}"
            );
            assert_eq!(found[0].func, "ok", "`{hidden}`: {found:?}");
        }
    }

    #[test]
    fn honours_every_delimiter_shape_perl_allows() {
        for delim in [
            "( )", "[ ]", "{ }", "< >", "! !", "| |", "/ /", "' '", "~ ~",
        ] {
            let (open, close) = delim.split_once(' ').unwrap();
            let src = format!("my $v = qq{open}text; not an assertion{close};\nok(1, 'x');\n");
            let found = enumerate_assertions(&src).unwrap();
            assert_eq!(found.len(), 1, "delimiter {delim}: {found:?}");
        }
    }

    /// The one delimiter deliberately **not** modelled, recorded so the choice is visible.
    /// `q#…#` is indistinguishable here from a bareword followed by a trailing comment,
    /// and the comment reading is overwhelmingly the likelier one. Nothing in the corpus
    /// uses it; if a future vendored file does, the scan fails loudly (the truncated line
    /// loses its `;`, and either the count moves or `done_testing` stops being seen)
    /// rather than under-counting quietly.
    #[test]
    fn documents_that_a_hash_delimiter_is_not_modelled() {
        let src = "my $v = q#text#;\nok(1, 'x');\n";
        assert!(
            enumerate_assertions(src).is_err(),
            "if this starts passing, `#` can be removed from the excluded-delimiter set"
        );
    }

    #[test]
    fn nests_paired_quote_delimiters() {
        let src = "my $v = q{ outer { inner } still inside; };\nok(1, 'x');\n";
        assert_eq!(enumerate_assertions(src).unwrap().len(), 1);
    }

    #[test]
    fn does_not_read_a_variable_named_like_an_operator_as_one() {
        // `$tr->…` was consumed as `tr-…-…-`, eating 90 lines of AnnotationSource_*.t.
        let src = "$tr->stable_id('foo');\nok(!$as->apply_edits($tr), 'fail - no alignment');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].func, "ok");
    }

    #[test]
    fn does_not_read_a_fat_comma_key_as_a_substitution() {
        // TranscriptTree.t has `{foo => 'bar', s => 5, e => 10}` nine times.
        let src = "is_deeply($t->fetch('c', 7, 8), [{foo => 'bar', s => 5, e => 10}], 'obj');\nok(1, 'x');\n";
        assert_eq!(enumerate_assertions(src).unwrap().len(), 2);
    }

    #[test]
    fn treats_slash_as_division_where_no_term_is_expected() {
        let src = "my $n = $h->{count} / 2;\nis($n, 4, 'divided');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "a phantom regex swallowed a `;`: {found:?}");
        assert_eq!(found[0].line, 2);
    }

    #[test]
    fn handles_a_regex_that_begins_on_the_line_after_the_binding() {
        // Utils.t rows 44's shape: `=~` ends L311, the regex body is L312.
        let src = "ok(\n  get_version_string($f) =~\n    /ensembl\\s+\\: 86; not a boundary/,\n  'get_version_string'\n);\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].line, 1);
    }

    #[test]
    fn handles_heredocs() {
        let src = "my $t = <<'EOF';\nis(1, 1, 'not real');\nEOF\nis(2, 2, 'real');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "heredoc body was scanned: {found:?}");
        assert_eq!(found[0].line, 4);
    }

    #[test]
    fn handles_indented_and_multiple_heredocs() {
        let src = "my ($a, $b) = (<<~A, <<B);\n  ok(0, 'inside a');\n  A\nok(0, 'inside b');\nB\nis(1, 1, 'real');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].line, 6);
    }

    #[test]
    fn does_not_mistake_a_left_shift_for_a_heredoc() {
        let src = "my $mask = 1 << 3;\nis($mask, 8, 'shifted');\n";
        assert_eq!(enumerate_assertions(src).unwrap().len(), 1);
    }

    #[test]
    fn reports_an_unterminated_construct_rather_than_under_counting() {
        let err = enumerate_assertions("is(1, 1, 'x');\nmy $r = qr/never closed;\n")
            .expect_err("must not silently truncate");
        assert!(err.contains("unterminated"), "unhelpful error: {err}");
    }

    // ---- defect 3: `#` is not unconditionally a comment ---------------------------

    #[test]
    fn a_hash_inside_qw_does_not_truncate_the_line() {
        // OutputFactory.t and Haplo_Parser_VCF.t build VCF headers as `qw(#CHROM POS …)`.
        // Truncating there hid the `;` and carried a stale boundary into the next line.
        let src = "my @h = (qw(#CHROM POS ID REF ALT), 'x');\nis(scalar @h, 6, 'header');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].line, 2);
    }

    #[test]
    fn a_last_index_sigil_is_not_a_comment() {
        for src in [
            "my $last = $#list; is($last, 2, 'last index');",
            "my $last = $#{$ref}; is($last, 2, 'last index');",
            "my $last = $#$ref; is($last, 2, 'last index');",
        ] {
            let found = enumerate_assertions(src).unwrap();
            assert_eq!(found.len(), 1, "`{src}` was truncated: {found:?}");
        }
    }

    #[test]
    fn a_hash_inside_a_substitution_does_not_truncate_the_line() {
        // AnnotationSourceAdaptor.t:286.
        let src =
            "$t =~ s/MT.vcf.gz/\\#\\#\\#CHR\\#\\#\\#\\.vcf.gz/;\nis($t, 'x', 'substituted');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].line, 2);
    }

    // ---- the call may be on a later line ------------------------------------------

    #[test]
    fn finds_a_call_whose_block_is_on_the_next_line() {
        // FilterSet.t writes this seven times; the same-line rule lost all seven.
        let src = "throws_ok\n  { $fs->parse_filters(['(']) }\n  qr/incomplete parentheses/,\n  'incomplete 1';\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(
            found[0].line, 1,
            "the ordinal must point at the `throws_ok`"
        );
    }

    #[test]
    fn a_bareword_that_is_not_called_is_not_an_assertion() {
        let src = "my %h = (like => 1, ok => 2);\nis(scalar keys %h, 2, 'real');\n";
        let found = enumerate_assertions(src).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
    }

    /// The dumb second opinion must agree with the lexer on the vendored file too — this
    /// is the same corroboration the corpus sweep applies to all 49, kept hermetic.
    #[test]
    fn the_second_opinion_agrees_on_utils_t() {
        let src = include_str!("port/perl/Utils.t");
        assert_eq!(count_line_anchored(src), 44);
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

/// One `#[test]` fn found in the crate.
#[derive(Debug, Clone, Default)]
struct RustTest {
    /// False when the fn is `#[ignore]`d without an adequate reason.
    ignore_ok: bool,
    /// 1-based line of every `fn` declaration that produced this key. More than one means
    /// two identically-named `#[test]` fns in different `mod`s of the same file.
    sites: Vec<usize>,
}

/// Does an `#[ignore …]` attribute carry an adequate reason?
///
/// The contract is `#[ignore = "vepyr#NN: <reason>"]`: an issue reference *and* a human
/// reason. The old check was `a.contains('=') && a.contains('#')` over the whole attribute
/// line — but `a` is only ever set for lines starting with `#[ignore`, so `contains('#')`
/// was always true and `#[ignore = "flaky"]` sailed through. This tests the attribute
/// *value*.
fn ignore_reason_is_adequate(attr: &str) -> bool {
    let Some((_, value)) = attr.split_once('=') else {
        return false; // a bare `#[ignore]` states no reason at all
    };
    let value = value.trim().trim_end_matches(']').trim();
    let value = value.trim_matches('"').trim();
    let Some((_, after_hash)) = value.split_once('#') else {
        return false; // no issue reference
    };
    let digits = after_hash.chars().take_while(char::is_ascii_digit).count();
    if digits == 0 {
        return false; // `#NN` placeholder or a stray `#`, not a real issue number
    }
    let reason = after_hash[digits..]
        .trim_start_matches([':', '-', ' '])
        .trim();
    !reason.is_empty()
}

/// Index of Rust test functions declared in the crate, keyed `"<rel path>::<fn>"`.
///
/// Built by scanning source text for a `#[test]` / `#[tokio::test]` attribute on the
/// declaration that follows. This is what makes the gate uncheatable: a comment claiming
/// coverage does not create an entry here.
///
/// The key is *not* module-qualified, so two `#[test] fn basic()` in different `mod`s of
/// one file collide. Rather than collapse them — which would let a ledger row claim
/// coverage that resolves ambiguously — every declaration site is recorded, and
/// [`check_ledger`] refuses a row that names an ambiguous key. See
/// `rejects_an_ambiguous_rust_target` for why this is preferred over module tracking.
fn index_rust_tests(roots: &[PathBuf], crate_root: &Path) -> BTreeMap<String, RustTest> {
    let mut idx: BTreeMap<String, RustTest> = BTreeMap::new();
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
        for (name, line, ignore_ok) in scan_rust_tests(&text) {
            let e = idx.entry(format!("{rel}::{name}")).or_default();
            // A collision is reported, never merged; `false` must win so an inadequate
            // `#[ignore]` cannot be masked by a same-named sibling.
            e.ignore_ok = if e.sites.is_empty() {
                ignore_ok
            } else {
                e.ignore_ok && ignore_ok
            };
            e.sites.push(line);
        }
    }
    idx
}

/// Find `(fn name, 1-based line, ignore-is-adequate)` for every `#[test]` fn in one file.
///
/// Attributes are consumed from the front of each trimmed line, so `#[test] fn foo()` on a
/// single line is recognised as well as the attribute-on-its-own-line form, while a
/// `#[test]` appearing inside a string or trailing comment is not (it is not at the front
/// of the line). `port_ledger.rs` itself contains `#[test]` inside error-message strings,
/// which is exactly the case that must not register.
fn scan_rust_tests(text: &str) -> Vec<(String, usize, bool)> {
    let mut out = Vec::new();
    let mut pending_test = false;
    let mut pending_ignore: Option<String> = None;
    for (idx, raw) in text.lines().enumerate() {
        let mut t = raw.trim();
        while t.starts_with("#[") {
            let Some(end) = attribute_end(t) else { break };
            let attr = &t[..=end];
            if attr == "#[test]" || attr == "#[tokio::test]" || attr.starts_with("#[tokio::test(") {
                pending_test = true;
            } else if attr.starts_with("#[ignore") {
                pending_ignore = Some(attr.to_owned());
            }
            t = t[end + 1..].trim();
        }
        if t.is_empty() || t.starts_with("//") {
            continue; // attribute-only or comment line: keep the pending state
        }
        if pending_test {
            if let Some(name) = rust_fn_name(t) {
                let ignore_ok = pending_ignore
                    .as_deref()
                    .is_none_or(ignore_reason_is_adequate);
                out.push((name, idx + 1, ignore_ok));
            }
            pending_test = false;
        }
        // Never let an `#[ignore]` leak past the item it annotates.
        pending_ignore = None;
    }
    out
}

/// Index of the `]` closing the attribute that starts at byte 0, honouring nesting and
/// string literals (`#[ignore = "a]b"]`).
fn attribute_end(t: &str) -> Option<usize> {
    let b = t.as_bytes();
    let mut depth = 0usize;
    let mut in_str = false;
    let mut i = 1usize; // past the `#`
    while i < b.len() {
        match b[i] {
            b'\\' if in_str => i += 1,
            b'"' => in_str = !in_str,
            b'[' if !in_str => depth += 1,
            b']' if !in_str => {
                depth -= 1;
                if depth == 0 {
                    return Some(i);
                }
            }
            _ => {}
        }
        i += 1;
    }
    None
}

/// The declared name of a `fn` item, given the item line with attributes stripped.
fn rust_fn_name(t: &str) -> Option<String> {
    let mut tok = t.split_whitespace();
    loop {
        match tok.next()? {
            "fn" => break,
            "pub" | "async" | "const" | "unsafe" | "default" => continue,
            w if w.starts_with("pub(") => continue,
            _ => return None,
        }
    }
    let name: String = tok
        .next()?
        .chars()
        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
        .collect();
    (!name.is_empty()).then_some(name)
}

#[cfg(test)]
mod rust_index {
    use super::*;

    #[test]
    fn requires_an_issue_and_a_reason_in_an_ignore_attribute() {
        // The old check was `contains('=') && contains('#')` over the whole attribute
        // *line*, and since the line always began `#[ignore`, the `'#'` test was free.
        assert!(!ignore_reason_is_adequate("#[ignore]"));
        assert!(!ignore_reason_is_adequate("#[ignore = \"flaky\"]"));
        assert!(!ignore_reason_is_adequate("#[ignore = \"vepyr#NN: why\"]"));
        assert!(!ignore_reason_is_adequate("#[ignore = \"vepyr#42\"]"));
        assert!(!ignore_reason_is_adequate("#[ignore = \"#42:   \"]"));
        assert!(ignore_reason_is_adequate(
            "#[ignore = \"vepyr#42: needs a real parquet cache\"]"
        ));
        assert!(ignore_reason_is_adequate(
            "#[ignore = \"biodatageeks/vepyr#7 — AlphaMissense fixture is 2 GB\"]"
        ));
    }

    #[test]
    fn records_every_declaration_site_instead_of_collapsing_them() {
        // Two same-named `#[test] fn`s in different `mod`s of one file. Keying on
        // `<file>::<fn>` cannot tell them apart; silently keeping one would let a ledger
        // row claim coverage that resolves to whichever the map happened to retain.
        //
        // Written with `\n` escapes rather than as a multi-line literal on purpose: the
        // index is a *textual* scan, so an `#[test]` at the start of a physical line is
        // registered even inside a string. That is the one place this heuristic is wrong,
        // and the crate-wide collision check below would trip over this very fixture.
        let text = "mod a {\n    #[test]\n    fn basic() {}\n}\n\
                    mod b {\n    #[test]\n    fn basic() {}\n}\n";
        let found = scan_rust_tests(text);
        assert_eq!(
            found
                .iter()
                .map(|(n, l, _)| (n.as_str(), *l))
                .collect::<Vec<_>>(),
            vec![("basic", 3), ("basic", 7)]
        );

        // …and the index keeps both sites rather than one.
        let dir = std::env::temp_dir().join("port_ledger_collision_fixture");
        std::fs::create_dir_all(&dir).expect("temp dir");
        let file = dir.join("dup.rs");
        std::fs::write(&file, text).expect("write fixture");
        let idx = index_rust_tests(std::slice::from_ref(&dir), &dir);
        assert_eq!(
            idx.get("dup.rs::basic").map(|t| t.sites.clone()),
            Some(vec![3, 7]),
            "the collision was collapsed: {idx:?}"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn recognises_attributes_on_the_same_line_as_the_fn() {
        let text = "#[test] fn one() {}\n#[tokio::test(flavor = \"multi_thread\")]\nasync fn two() {}\n#[test]\n#[ignore = \"x\"]\npub fn three() {}\n";
        let found = scan_rust_tests(text);
        assert_eq!(
            found
                .iter()
                .map(|(n, _, ok)| (n.as_str(), *ok))
                .collect::<Vec<_>>(),
            vec![("one", true), ("two", true), ("three", false)]
        );
    }

    #[test]
    fn does_not_register_a_test_attribute_quoted_inside_a_string() {
        // `port_ledger.rs` itself carries `#[test]` and `#[ignore]` inside error messages.
        let text = "let msg = \"must resolve to a real #[test] fn\";\nfn helper() {}\n";
        assert!(scan_rust_tests(text).is_empty());
    }

    #[test]
    fn does_not_let_an_ignore_leak_onto_a_later_test() {
        let text = "#[ignore = \"x\"]\nstruct NotATest;\n#[test]\nfn real() {}\n";
        let found = scan_rust_tests(text);
        assert_eq!(found.len(), 1);
        assert!(found[0].2, "the stale #[ignore] leaked: {found:?}");
    }

    #[test]
    fn indexes_this_very_crate_without_a_single_collision() {
        // Cheap, honest tripwire: 870 `#[test]` fns today, 870 distinct keys. If a future
        // commit introduces a collision this fails here, before any ledger points at it.
        let root = crate_root();
        let idx = index_rust_tests(&[root.join("src"), root.join("tests")], &root);
        let collisions: Vec<_> = idx
            .iter()
            .filter(|(_, t)| t.sites.len() > 1)
            .map(|(k, t)| format!("{k} at {:?}", t.sites))
            .collect();
        assert!(
            collisions.is_empty(),
            "same-named #[test] fns in one file — rename so a ledger can name one \
             unambiguously:\n  - {}",
            collisions.join("\n  - ")
        );
        assert!(idx.len() > 500, "index looks empty: {} entries", idx.len());
    }
}

/// All the ways a ledger can fail its contract. Returns one message per violation so a
/// contributor sees every problem at once instead of fixing them one build at a time.
fn check_ledger(
    ledger: &Ledger,
    assertions: &[Assertion],
    rust_tests: &BTreeMap<String, RustTest>,
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

    // `expected_assertions` is a tripwire on the *enumerator*; this is the matching
    // tripwire on the *ledger*. Without it a ledger could declare 44 and carry 40 rows,
    // and only the four resulting "NO ledger row" messages would hint at why.
    if ledger.rows.len() != ledger.expected_assertions {
        errs.push(format!(
            "{name}: ledger declares expected_assertions = {} but carries {} [[row]] entries \
             — every assertion needs exactly one row.",
            ledger.expected_assertions,
            ledger.rows.len()
        ));
    }

    // A duplicate `n` used to be invisible: building `by_n` silently overwrote the first
    // row, and the *other* ordinal it displaced then reported a confusing "NO ledger row".
    let mut positions: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, r) in ledger.rows.iter().enumerate() {
        positions.entry(r.n).or_default().push(i + 1);
    }
    for (n, at) in &positions {
        if at.len() > 1 {
            errs.push(format!(
                "{name}: n={n} is claimed by {} rows (file positions {at:?}) — each Perl \
                 assertion must be classified exactly once, so a duplicate hides whichever \
                 ordinal it displaced.",
                at.len()
            ));
        }
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
            Class::UnitPort => match r.rust.as_deref().map(|t| (t, rust_tests.get(t))) {
                None => errs.push(format!(
                    "{name}: row n={} is unit-port but has no `rust` field",
                    r.n
                )),
                Some((target, None)) => errs.push(format!(
                    "{name}: row n={} names Rust test `{target}`, which does not exist. \
                     A claim of coverage must resolve to a real #[test] fn.",
                    r.n
                )),
                Some((target, Some(t))) if t.sites.len() > 1 => errs.push(format!(
                    "{name}: row n={} names `{target}`, but that file declares {} #[test] fns \
                     with that name (lines {:?}) — probably in different `mod`s. The claim is \
                     ambiguous: rename one so the ledger names exactly one test.",
                    r.n,
                    t.sites.len(),
                    t.sites
                )),
                Some((target, Some(t))) if !t.ignore_ok => errs.push(format!(
                    "{name}: row n={} names `{target}`, which is #[ignore]d without an adequate \
                     reason. Use #[ignore = \"vepyr#<issue number>: why\"].",
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

    fn tests_with(name: &str) -> BTreeMap<String, RustTest> {
        let mut m = BTreeMap::new();
        m.insert(
            name.to_owned(),
            RustTest {
                ignore_ok: true,
                sites: vec![1],
            },
        );
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
    fn flags_a_duplicate_n_directly() {
        // Two rows claiming n=1 used to overwrite silently in `by_n`; the *displaced*
        // ordinal (n=2) then reported "NO ledger row", which points at the wrong problem.
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"architectural-no-analogue\"\nreason=\"x\"\n\
             [[row]]\nn=1\nline=11\ndesc=\"b\"\nclass=\"architectural-no-analogue\"\nreason=\"y\"\n",
            2,
        );
        let errs = check_ledger(&l, &assertions(2), &BTreeMap::new());
        assert!(
            errs.iter().any(|e| e.contains("n=1 is claimed by 2 rows")),
            "{errs:?}"
        );
    }

    #[test]
    fn flags_a_row_count_that_does_not_match_expected_assertions() {
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"architectural-no-analogue\"\nreason=\"x\"\n",
            2,
        );
        // The enumerator agrees with `expected_assertions`; only the ledger is short.
        let errs = check_ledger(&l, &assertions(2), &BTreeMap::new());
        assert!(
            errs.iter().any(|e| e.contains("carries 1 [[row]] entries")),
            "{errs:?}"
        );
    }

    #[test]
    fn rejects_an_ambiguous_rust_target() {
        let mut tests = BTreeMap::new();
        tests.insert(
            "tests/p.rs::basic".to_owned(),
            RustTest {
                ignore_ok: true,
                sites: vec![12, 98],
            },
        );
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"unit-port\"\nrust=\"tests/p.rs::basic\"\n",
            1,
        );
        let errs = check_ledger(&l, &assertions(1), &tests);
        assert!(
            errs.iter()
                .any(|e| e.contains("ambiguous") && e.contains("[12, 98]")),
            "{errs:?}"
        );
    }

    #[test]
    fn rejects_an_ignored_target_whose_reason_is_not_a_tracked_issue() {
        let mut tests = BTreeMap::new();
        tests.insert(
            "tests/p.rs::flaky".to_owned(),
            RustTest {
                ignore_ok: false,
                sites: vec![7],
            },
        );
        let l = ledger_with(
            "[[row]]\nn=1\nline=10\ndesc=\"a\"\nclass=\"unit-port\"\nrust=\"tests/p.rs::flaky\"\n",
            1,
        );
        let errs = check_ledger(&l, &assertions(1), &tests);
        assert!(
            errs.iter().any(|e| e.contains("without an adequate")),
            "{errs:?}"
        );
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

/// Walk up from the crate towards the filesystem root looking for `rel`, so the corpus
/// sweep works from a git worktree as well as from the primary checkout.
fn find_upwards(rel: &str) -> Option<PathBuf> {
    let mut dir: &Path = &crate_root();
    loop {
        let candidate = dir.join(rel);
        if candidate.is_dir() {
            return Some(candidate);
        }
        dir = dir.parent()?;
    }
}

/// The corpus-wide sweep: proof that the enumerator is right on all 49 upstream `.t`
/// files, not merely on the one it was built against.
///
/// Not hermetic — it reads an `ensembl-vep` checkout and the hand-written port plans, so
/// it is `#[ignore]`d and CI stays offline-clean. Run it with:
///
/// ```text
/// cargo test -p datafusion-bio-function-vep --test port_ledger --offline \
///   -- --ignored --nocapture corpus::
/// ```
#[cfg(test)]
mod corpus {
    use super::*;

    /// Where an independent count came from.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Src {
        /// A numbered per-subtest mapping table in `porting-tests/detailed_plans/<T>.md`.
        /// The strongest available source: it enumerates rows one by one.
        Table,
        /// A prose claim in the same plan — `"The 52 subtests cover…"`, `"(~19 subtests)"`,
        /// `"1329 lines, 107 assertions"`. Often explicitly approximate.
        Prose,
    }

    /// One corpus file: the count the enumerator must produce, the independent count, and
    /// the adjudication of any gap. `scanner` is **frozen** — if a future change to the
    /// lexer moves any of these numbers, this test fails and the change must be justified.
    struct Expect {
        file: &'static str,
        scanner: usize,
        independent: Option<(usize, Src)>,
        /// Empty when the two agree.
        verdict: &'static str,
    }

    const COLLAPSE: &str =
        "plan table collapses assertion clusters into fewer rows; scanner counts statements";
    const PROSE: &str = "plan gives only an approximate prose total, no per-subtest table";

    #[rustfmt::skip]
    const EXPECTED: &[Expect] = &[
        Expect { file: "AnnotationSource.t", scanner: 13, independent: Some((19, Src::Table)),
            verdict: "SCANNER RIGHT: the plan expands the single `is_deeply(...) for @$inputs` \
                      statement at L97 into rows 13a-13g. 12 statements + 1 loop = 13." },
        Expect { file: "AnnotationSourceAdaptor.t", scanner: 29, independent: Some((28, Src::Table)),
            verdict: "SCANNER RIGHT: 20 table rows with the 9-assertion database SKIP block \
                      'counted as one for tally purposes' (plan L72) = 28; that tally is loose." },
        Expect { file: "AnnotationSource_Cache.t", scanner: 4, independent: Some((4, Src::Table)), verdict: "" },
        Expect { file: "AnnotationSource_Cache_RegFeat.t", scanner: 53, independent: Some((37, Src::Table)),
            verdict: "SCANNER RIGHT: the plan tabulates only the in-scope subset (its L126 \
                      reasons about '~37 subtests'), not a full enumeration." },
        Expect { file: "AnnotationSource_Cache_Transcript.t", scanner: 96, independent: Some((55, Src::Table)), verdict: COLLAPSE },
        Expect { file: "AnnotationSource_Cache_Variation.t", scanner: 54, independent: Some((40, Src::Table)), verdict: COLLAPSE },
        Expect { file: "AnnotationSource_Cache_VariationTabix.t", scanner: 52, independent: Some((48, Src::Table)), verdict: COLLAPSE },
        Expect { file: "AnnotationSource_Database_RegFeat.t", scanner: 33, independent: Some((32, Src::Prose)), verdict: PROSE },
        Expect { file: "AnnotationSource_Database_StructuralVariation.t", scanner: 20, independent: Some((19, Src::Prose)), verdict: PROSE },
        Expect { file: "AnnotationSource_Database_Transcript.t", scanner: 97, independent: Some((88, Src::Prose)), verdict: PROSE },
        Expect { file: "AnnotationSource_Database_Variation.t", scanner: 27, independent: Some((25, Src::Prose)), verdict: PROSE },
        Expect { file: "AnnotationSource_File.t", scanner: 20, independent: Some((21, Src::Prose)),
            verdict: "SCANNER RIGHT on statements: L110's `is()` sits inside a 5-format \
                      foreach, so TAP emits more results than there are statements." },
        Expect { file: "AnnotationSource_File_BED.t", scanner: 16, independent: Some((15, Src::Prose)), verdict: PROSE },
        Expect { file: "AnnotationSource_File_BigWig.t", scanner: 21, independent: Some((20, Src::Prose)), verdict: PROSE },
        Expect { file: "AnnotationSource_File_GFF.t", scanner: 62, independent: Some((52, Src::Prose)),
            verdict: "prose approximation; includes 3 bare `warning {…}` (a judgement call, \
                      see ASSERTION_FNS) and 5 use_ok" },
        Expect { file: "AnnotationSource_File_GTF.t", scanner: 29, independent: Some((25, Src::Prose)),
            verdict: "prose approximation; includes 2 bare `warning {…}` and 5 use_ok" },
        Expect { file: "AnnotationSource_File_VCF.t", scanner: 33, independent: Some((32, Src::Prose)), verdict: PROSE },
        Expect { file: "BaseVEP.t", scanner: 50, independent: Some((47, Src::Prose)),
            verdict: "plan is internally inconsistent: 39 table rows vs a prose 'all 47 \
                      Perl subtests'. Neither is an enumeration." },
        Expect { file: "CacheDir.t", scanner: 44, independent: Some((24, Src::Table)),
            verdict: "SCANNER RIGHT: the plan tabulates only its bucket-2 rows" },
        Expect { file: "Config.t", scanner: 36, independent: Some((33, Src::Table)),
            verdict: "SCANNER RIGHT: the plan claims 'all 33 Perl subtests' but collapses \
                      clusters (e.g. row 32 covers both L142 and L143)" },
        Expect { file: "FilterSet.t", scanner: 160, independent: None,
            verdict: "NO independent source (plan has no table and no total). 7 of these 160 \
                      are `throws_ok` with the block on the next line — invisible to the old \
                      same-line rule, and nothing would have caught it." },
        Expect { file: "Haplo_AnnotationSource_Cache_Transcript.t", scanner: 14, independent: Some((13, Src::Prose)), verdict: PROSE },
        Expect { file: "Haplo_AnnotationSource_Database_Transcript.t", scanner: 12, independent: Some((11, Src::Prose)), verdict: PROSE },
        Expect { file: "Haplo_AnnotationSource_File_GFF.t", scanner: 13, independent: Some((11, Src::Prose)), verdict: PROSE },
        Expect { file: "Haplo_AnnotationSource_File_GTF.t", scanner: 13, independent: Some((11, Src::Prose)), verdict: PROSE },
        Expect { file: "Haplo_InputBuffer.t", scanner: 18, independent: Some((14, Src::Prose)), verdict: PROSE },
        Expect { file: "Haplo_Parser_VCF.t", scanner: 15, independent: Some((13, Src::Prose)), verdict: PROSE },
        Expect { file: "Haplo_Runner.t", scanner: 26, independent: Some((26, Src::Prose)), verdict: "" },
        Expect { file: "InputBuffer.t", scanner: 76, independent: Some((39, Src::Table)), verdict: COLLAPSE },
        Expect { file: "OutputFactory.t", scanner: 145, independent: Some((141, Src::Prose)),
            verdict: "plan says 'approximately 141 assertions in 101 mapping rows' — it knows \
                      it is approximate" },
        Expect { file: "OutputFactory_JSON.t", scanner: 24, independent: Some((18, Src::Prose)), verdict: PROSE },
        Expect { file: "OutputFactory_Tab.t", scanner: 17, independent: Some((16, Src::Table)), verdict: COLLAPSE },
        Expect { file: "OutputFactory_VCF.t", scanner: 57, independent: Some((56, Src::Prose)),
            verdict: "plan's line count (754) is exact, so its 56 is a careful count off by one" },
        Expect { file: "OutputFactory_VEP_output.t", scanner: 20, independent: Some((19, Src::Table)), verdict: COLLAPSE },
        Expect { file: "Parser.t", scanner: 112, independent: Some((112, Src::Prose)), verdict: "" },
        Expect { file: "Parser_CAID.t", scanner: 11, independent: Some((8, Src::Prose)), verdict: PROSE },
        Expect { file: "Parser_HGVS.t", scanner: 29, independent: Some((24, Src::Prose)), verdict: PROSE },
        Expect { file: "Parser_ID.t", scanner: 14, independent: Some((11, Src::Prose)), verdict: PROSE },
        Expect { file: "Parser_Region.t", scanner: 17, independent: Some((14, Src::Prose)), verdict: PROSE },
        Expect { file: "Parser_SPDI.t", scanner: 6, independent: Some((3, Src::Prose)), verdict: PROSE },
        Expect { file: "Parser_VCF.t", scanner: 78, independent: Some((107, Src::Prose)),
            verdict: "SCANNER RIGHT: the plans transposed two numbers. Parser_VCF.md says \
                      '1329 lines, 107 assertions' and Runner.md says '1348 lines, 78 \
                      assertions'; both line counts verify exactly, but the assertion counts \
                      are each other's. Scanner: Parser_VCF 78, Runner 107." },
        Expect { file: "Parser_VEP_input.t", scanner: 18, independent: Some((16, Src::Prose)), verdict: PROSE },
        Expect { file: "Runner.t", scanner: 107, independent: Some((78, Src::Prose)),
            verdict: "SCANNER RIGHT: transposed with Parser_VCF.t — see that row" },
        Expect { file: "Stats.t", scanner: 29, independent: Some((29, Src::Table)), verdict: "" },
        Expect { file: "TranscriptTree.t", scanner: 38, independent: Some((39, Src::Table)),
            verdict: "SCANNER RIGHT: plan row 23 records the *set-up* call \
                      `$t->insert('chrobj', 5, 10, {foo => 'bar'})` at L89, which asserts \
                      nothing. Plan rows 24-39 line up with scanner 23-38." },
        Expect { file: "Utils.t", scanner: 44, independent: Some((44, Src::Table)),
            verdict: "" },
        Expect { file: "VariantRecoder.t", scanner: 29, independent: Some((29, Src::Prose)), verdict: "" },
        Expect { file: "bam_edit.t", scanner: 32, independent: Some((28, Src::Prose)), verdict: PROSE },
        Expect { file: "version.t", scanner: 3, independent: Some((3, Src::Prose)), verdict: "" },
    ];

    /// Path to the upstream `ensembl-vep/t` corpus.
    fn corpus_dir() -> PathBuf {
        std::env::var_os("VEP_PERL_CORPUS")
            .map(PathBuf::from)
            .or_else(|| find_upwards("ensembl-vep/t"))
            .expect("set VEP_PERL_CORPUS, or put an ensembl-vep checkout beside this repo")
    }

    fn t_files() -> Vec<PathBuf> {
        let mut v: Vec<PathBuf> = std::fs::read_dir(corpus_dir())
            .expect("readable corpus dir")
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| p.extension().and_then(|e| e.to_str()) == Some("t"))
            .collect();
        v.sort();
        v
    }

    #[test]
    #[ignore = "corpus sweep: reads an ensembl-vep checkout and porting-tests/detailed_plans, \
                neither of which exists in CI. Run with --ignored --nocapture."]
    fn the_enumerator_is_right_across_the_whole_corpus() {
        let files = t_files();
        assert_eq!(
            files.len(),
            EXPECTED.len(),
            "corpus has {} .t files but the adjudication table has {} rows — a new or \
             removed upstream test must be adjudicated, not skipped",
            files.len(),
            EXPECTED.len()
        );

        let mut errs: Vec<String> = Vec::new();
        let mut disagreements: Vec<String> = Vec::new();
        println!(
            "\n{:<48} {:>7} {:>6} {:>7}  source / note",
            "file", "scanner", "2nd", "indep"
        );
        println!("{:-<110}", "");

        for p in &files {
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            let src = std::fs::read_to_string(p).expect("readable .t");
            let Some(exp) = EXPECTED.iter().find(|e| e.file == name) else {
                errs.push(format!("{name}: not in the adjudication table"));
                continue;
            };

            // (1) A file the scanner cannot parse confidently fails loudly.
            let scan = match scan_perl(&src) {
                Ok(s) => s,
                Err(e) => {
                    errs.push(format!("{name}: UNPARSEABLE — {e}"));
                    continue;
                }
            };
            let got = scan.assertions.len();

            // (2) Frozen count: a silent drift here is the failure mode being fixed.
            if got != exp.scanner {
                errs.push(format!(
                    "{name}: enumerator found {got}, table says {} — a scanner change moved \
                     this count. Prove which is right before editing the table.",
                    exp.scanner
                ));
            }

            // (3) The dumb second opinion must agree. This shares only the function-name
            //     list with the lexer, so it independently rules out both a swallowed
            //     boundary and a phantom one.
            let second = count_line_anchored(&src);
            if second != got {
                errs.push(format!(
                    "{name}: lexer says {got} but the line-anchored second opinion says \
                     {second} — one of them is wrong, and neither may be trusted until \
                     the difference is explained."
                ));
            }

            // (4) Parse-health landmark: every upstream `.t` ends in `done_testing();`, so
            //     failing to see it at a statement boundary means the boundary state
            //     desynchronised — an unterminated construct that ran to EOF, typically.
            let textual = src.matches("done_testing(").count();
            let seen = scan.declarations_named("done_testing");
            if seen != textual || textual == 0 {
                errs.push(format!(
                    "{name}: {seen} of {textual} `done_testing(` calls were seen at a \
                     statement boundary — the scan desynchronised."
                ));
            }

            let (indep, src_kind) = match exp.independent {
                Some((n, s)) => (n.to_string(), format!("{s:?}")),
                None => ("-".into(), "Absent".into()),
            };
            println!(
                "{:<48} {got:>7} {second:>6} {indep:>7}  {src_kind}{}",
                name,
                if exp.verdict.is_empty() {
                    "  (agree)".to_string()
                } else {
                    String::new()
                }
            );
            if let Some((n, _)) = exp.independent
                && n != got
            {
                assert!(
                    !exp.verdict.is_empty(),
                    "{name}: scanner {got} vs independent {n} with no recorded verdict"
                );
                disagreements.push(format!("{name}: scanner {got} vs {n} — {}", exp.verdict));
            }
        }

        println!(
            "\n{} disagreement(s) with the independent source:",
            disagreements.len()
        );
        for d in &disagreements {
            println!("  - {d}");
        }
        let agreed = EXPECTED
            .iter()
            .filter(|e| e.independent.is_some_and(|(n, _)| n == e.scanner))
            .count();
        println!(
            "\n{agreed}/{} files agree exactly with the independent source; \
             {} have no independent source.",
            EXPECTED.len(),
            EXPECTED.iter().filter(|e| e.independent.is_none()).count()
        );

        assert!(
            errs.is_empty(),
            "corpus sweep failed with {} problem(s):\n  - {}",
            errs.len(),
            errs.join("\n  - ")
        );
    }

    /// Why the plans are the only independent source available: the upstream corpus
    /// declares no test plans at all, so there is no `plan tests => N` to compare against.
    #[test]
    #[ignore = "corpus sweep: reads an ensembl-vep checkout, which CI does not have"]
    fn the_corpus_declares_no_test_plans() {
        let with_plan: Vec<String> = t_files()
            .iter()
            .filter(|p| {
                std::fs::read_to_string(p)
                    .map(|s| s.contains("plan tests") || s.contains("plan skip_all"))
                    .unwrap_or(false)
            })
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        assert!(
            with_plan.is_empty(),
            "these files declare a plan and should be cross-checked against it directly, \
             which is stronger than the hand-written plans: {with_plan:?}"
        );
    }
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
