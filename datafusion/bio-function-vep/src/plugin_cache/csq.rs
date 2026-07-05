//! CSQ formatting helpers for plugin fields (spec §5.2).
//!
//! Plugin fields are appended as trailing pipe-delimited CSQ fields. The count is
//! fixed per run (`registry.csq_fields().len()`), so every CSQ entry appends
//! exactly that many fields (empty on a miss) — keeping header/body width aligned.
//! With no plugins the suffix is empty, so output is byte-identical to before.

use crate::plugin_cache::lookup::PluginScalar;

/// Form the AlphaMissense-style amino-acid change discriminator from the
/// transcript's `Amino_acids` (`"V/A"`) + `Protein_position` (`"550"`) →
/// `"V550A"`. `None` when there is no single-residue substitution (no
/// `Amino_acids`, not `X/Y` form, or a positional range) — a non-missense
/// consequence, which then misses the plugin lookup (the gate).
pub fn amino_acid_change(amino_acids: &str, protein_position: &str) -> Option<String> {
    let (ref_aa, alt_aa) = amino_acids.split_once('/')?;
    if ref_aa.len() != 1 || alt_aa.len() != 1 {
        return None;
    }
    if protein_position.is_empty() || protein_position.contains('-') {
        return None;
    }
    Some(format!("{ref_aa}{protein_position}{alt_aa}"))
}

/// Format one plugin scalar for CSQ output: floats via shortest round-trip,
/// strings verbatim, `Null` → empty.
pub fn format_scalar(scalar: &PluginScalar) -> String {
    match scalar {
        PluginScalar::Str(s) => s.clone(),
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
    fn amino_acid_change_forms_discriminator() {
        assert_eq!(amino_acid_change("V/A", "550").as_deref(), Some("V550A"));
        assert_eq!(amino_acid_change("W/R", "320").as_deref(), Some("W320R"));
        // non-missense / no aa-change → None (gate)
        assert_eq!(amino_acid_change("", "550"), None);
        assert_eq!(amino_acid_change("V/A", ""), None);
        assert_eq!(amino_acid_change("V", "550"), None); // no '/'
        assert_eq!(amino_acid_change("VV/A", "550"), None); // multi-residue
        assert_eq!(amino_acid_change("V/A", "550-551"), None); // range
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
