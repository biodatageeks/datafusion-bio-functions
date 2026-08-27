//! CSQ formatting helpers for plugin fields (spec §5.2).
//!
//! Plugin fields are appended as trailing pipe-delimited CSQ fields. The count is
//! fixed per run (`registry.csq_fields().len()`), so every CSQ entry appends
//! exactly that many fields (empty on a miss) — keeping header/body width aligned.
//! With no plugins the suffix is empty, so output is byte-identical to before.

use crate::plugin_cache::lookup::PluginScalar;

/// Escape a plugin string value for the CSQ payload, mirroring the engine's
/// built-in `csq_escape`: `,`/`|` → `&`, `;` → `%3B`, `=` → `%3D`,
/// whitespace → `_`. Without this a Utf8 plugin value containing a CSQ/INFO
/// delimiter would corrupt field or entry boundaries. (The built-in
/// `-`→empty convention is deliberately NOT applied to plugin values, so a
/// legitimate `-` is preserved.)
fn escape_csq_value(val: &str) -> String {
    let mut out = String::with_capacity(val.len());
    for ch in val.chars() {
        match ch {
            ',' | '|' => out.push('&'),
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            c if c.is_whitespace() => out.push('_'),
            c => out.push(c),
        }
    }
    out
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
