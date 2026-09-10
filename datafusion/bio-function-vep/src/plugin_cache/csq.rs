//! CSQ formatting helpers for plugin fields (spec §5.2; value escaping §5.5).
//!
//! Plugin fields are appended as trailing pipe-delimited CSQ fields. The count is
//! fixed per run (`registry.csq_fields().len()`), so every CSQ entry appends
//! exactly that many fields (empty on a miss) — keeping header/body width aligned.
//! With no plugins the suffix is empty, so output is byte-identical to before.

use crate::annotate_provider::csq_escape;
use crate::plugin_cache::lookup::PluginScalar;

/// Escape a plugin string value for the CSQ payload.
///
/// The VEP rules themselves live in one place -- `annotate_provider::csq_escape`
/// -- so the built-in and plugin paths cannot drift apart again (they had, which
/// is what vepyr#93 reported). Plugin columns get exactly the same treatment as
/// built-in ones because VEP treats them the same: `OutputFactory/VCF.pm:449`
/// pushes `get_plugin_headers` into `@fields`, and `:387` is the loop that
/// escapes every entry of that list -- including the `-`→empty rule at `:396-398`.
///
/// One deliberate deviation, applied here and nowhere else: `=` → `%3D`.
///
/// VEP 116 does **not** escape `=` in a CSQ value; its only `=`→`%3D` is HGVSp
/// at `OutputFactory.pm:1757`, and that one is `no_escape`-gated. We do it anyway
/// because ClinVar `CLNVI` carries BIC `base_change=…` values that VEP emits
/// already-`%3D`-encoded from the source VCF (VCF 4.3 requires it in an INFO
/// value); without this, 556 `ClinVar_CLNVI` entries across 7 records mismatch.
/// The real cause is an upstream INFO percent-decode, tracked separately -- when
/// that is fixed this deviation should be removed, not moved.
fn escape_csq_value(val: &str) -> String {
    let escaped = csq_escape(val);
    if escaped.contains('=') {
        // Order-safe: no VEP substitution emits an `=`, and none of `&`, `%3B`
        // or `_` contains one, so a post-pass cannot corrupt an earlier result.
        escaped.replace('=', "%3D")
    } else {
        escaped.into_owned()
    }
}

/// Format one plugin scalar for CSQ output: floats via shortest round-trip,
/// strings CSQ-escaped, `Null` → empty.
pub fn format_scalar(scalar: &PluginScalar) -> String {
    match scalar {
        PluginScalar::Str(s) => escape_csq_value(s),
        PluginScalar::F32(v) => format!("{v}"),
        PluginScalar::I32(v) => format!("{v}"),
        PluginScalar::Null => String::new(),
    }
}

/// Trailing suffix for a resolved probe: `|f1|f2|…|fN` (each field may be empty).
/// Empty string when `scalars` is empty (no plugins) → byte-identical output.
pub fn field_suffix(scalars: &[PluginScalar]) -> String {
    let mut out = String::new();
    for s in scalars {
        out.push('|');
        out.push_str(&format_scalar(s));
    }
    out
}

/// Trailing suffix of `n` empty plugin fields: `|` repeated `n` times. Used where
/// no per-transcript discriminator is available (e.g. the cached/empty paths),
/// keeping width aligned. `n == 0` → empty string.
pub fn empty_suffix(n: usize) -> String {
    "|".repeat(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_csq_delimiters_in_string_values() {
        // delimiters that would corrupt the CSQ/INFO payload are escaped
        assert_eq!(
            format_scalar(&PluginScalar::Str("a|b;c d".into())),
            "a&b%3Bc_d"
        );
        assert_eq!(format_scalar(&PluginScalar::Str("x,y".into())), "x&y");
        assert_eq!(
            format_scalar(&PluginScalar::Str("base_change=G_to_A".into())),
            "base_change%3DG_to_A"
        );
        // AlphaMissense-style values contain none → unchanged (parity-safe)
        assert_eq!(
            format_scalar(&PluginScalar::Str("likely_benign".into())),
            "likely_benign"
        );
        // floats/ints/null are unaffected
        assert_eq!(format_scalar(&PluginScalar::F32(0.2199)), "0.2199");
        assert_eq!(format_scalar(&PluginScalar::Null), "");
    }

    // ── vepyr#93 ──────────────────────────────────────────────────────────

    #[test]
    fn plugin_and_builtin_escapers_are_one_function() {
        // The two escapers must agree on everything VEP specifies. Anything
        // they disagree on is a deviation that has to be justified in the
        // spec, not an accident of having two copies.
        use crate::annotate_provider::csq_escape;
        for v in ["A  G", "x,y|z;w", "a\t\tb", "-", "  x  ", "likely_benign"] {
            assert_eq!(
                format_scalar(&PluginScalar::Str(v.to_string())),
                csq_escape(v),
                "the two escapers must be one function (value: {v:?})"
            );
        }
    }

    #[test]
    fn plugin_values_collapse_whitespace_runs() {
        // OutputFactory/VCF.pm:403 `s/\s+/\_/g` applies to plugin columns too:
        // VCF.pm:449 pushes get_plugin_headers into @fields, which is the list
        // the escaping loop at VCF.pm:387 walks.
        assert_eq!(
            format_scalar(&PluginScalar::Str("two  spaces".into())),
            "two_spaces"
        );
        assert_eq!(format_scalar(&PluginScalar::Str("a\t  b".into())), "a_b");
    }

    #[test]
    fn plugin_dash_is_blanked_like_vep() {
        // BEHAVIOUR FLIP (vepyr#93 D1). The old plugin escaper deliberately
        // preserved a bare `-`; VEP does not. VCF.pm:396-398 blanks `-` for
        // every column except Allele, and VCF.pm:449 puts plugin columns in
        // that same list. Unreachable in every published cache today (0 hits
        // across 175M rows, 5 caches), so no digest moves -- this test is the
        // record of the decision.
        assert_eq!(format_scalar(&PluginScalar::Str("-".into())), "");
        // A `-` INSIDE a value is untouched, as in VEP.
        assert_eq!(format_scalar(&PluginScalar::Str("a-b".into())), "a-b");
    }

    #[test]
    fn plugin_equals_escaping_is_a_deliberate_deviation() {
        // VEP 116 does NOT escape `=` in a CSQ value; its only `=`->%3D is
        // HGVSp (OutputFactory.pm:1757), and that one is no_escape-gated.
        // Kept anyway: ClinVar_CLNVI ships BIC `base_change=` values that VEP
        // emits already-%3D-encoded from the source VCF, and without this 556
        // entries across 7 records mismatch. Root cause is an upstream INFO
        // decode, tracked separately.
        assert_eq!(
            format_scalar(&PluginScalar::Str("base_change=G_to_A".into())),
            "base_change%3DG_to_A"
        );
        // This is the ONE value where the two escapers legitimately differ.
        use crate::annotate_provider::csq_escape;
        assert_eq!(csq_escape("a=b"), "a=b");
        assert_ne!(
            format_scalar(&PluginScalar::Str("a=b".into())),
            csq_escape("a=b")
        );
    }

    #[test]
    fn suffix_widths_and_formatting() {
        assert_eq!(field_suffix(&[]), "");
        assert_eq!(empty_suffix(0), "");
        assert_eq!(empty_suffix(2), "||");
        let scalars = vec![
            PluginScalar::F32(0.0827),
            PluginScalar::Str("likely_benign".into()),
        ];
        assert_eq!(field_suffix(&scalars), "|0.0827|likely_benign");
        // a miss row → two empty fields, same width
        assert_eq!(
            field_suffix(&[PluginScalar::Null, PluginScalar::Null]),
            "||"
        );
    }
}
