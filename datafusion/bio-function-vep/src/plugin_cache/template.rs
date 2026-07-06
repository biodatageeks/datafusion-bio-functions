//! Declarative match-discriminator templates over a fixed engine-attribute
//! namespace. A manifest `[[match_column]].template` (e.g.
//! `"{ref_aa}{Protein_position}{alt_aa}"`) is compiled once per plugin into
//! index-resolved segments and evaluated per transcript consequence. If any
//! referenced attribute is `None`, the discriminator is `None` (probe miss →
//! empty output — the per-transcript gate).

use datafusion::common::{DataFusionError, Result};

/// The engine attributes a template may reference, in fixed order. Extend only
/// when a plugin needs a value not already exposed by the transcript engine.
pub const ATTR_NAMES: &[&str] = &[
    "Consequence",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "ref_aa",
    "alt_aa",
    "ref",
    "alt",
];

/// Resolve an attribute name to its index in [`ATTR_NAMES`].
pub fn attr_index(name: &str) -> Option<usize> {
    ATTR_NAMES.iter().position(|&n| n == name)
}

#[derive(Debug, Clone, PartialEq)]
enum Segment {
    Lit(String),
    Attr(usize),
}

/// A template compiled to index-resolved segments (no per-record parsing).
#[derive(Debug, Clone)]
pub struct CompiledTemplate {
    segments: Vec<Segment>,
}

impl CompiledTemplate {
    /// Parse a `{name}` template; unknown attribute names are a hard error.
    pub fn compile(template: &str) -> Result<CompiledTemplate> {
        let mut segments = Vec::new();
        let mut rest = template;
        while let Some(open) = rest.find('{') {
            if open > 0 {
                segments.push(Segment::Lit(rest[..open].to_string()));
            }
            let close = rest[open..].find('}').ok_or_else(|| {
                DataFusionError::Execution(format!("unterminated '{{' in template '{template}'"))
            })? + open;
            let name = &rest[open + 1..close];
            let idx = attr_index(name).ok_or_else(|| {
                DataFusionError::Execution(format!("unknown template attribute '{name}'"))
            })?;
            segments.push(Segment::Attr(idx));
            rest = &rest[close + 1..];
        }
        if !rest.is_empty() {
            segments.push(Segment::Lit(rest.to_string()));
        }
        Ok(CompiledTemplate { segments })
    }

    /// Evaluate against a namespace array (same order as [`ATTR_NAMES`]). Any
    /// referenced attribute `None` → whole discriminator `None`.
    pub fn eval(&self, attrs: &[Option<&str>]) -> Option<String> {
        let mut out = String::new();
        for seg in &self.segments {
            match seg {
                Segment::Lit(s) => out.push_str(s),
                Segment::Attr(i) => out.push_str(attrs[*i]?),
            }
        }
        Some(out)
    }
}

/// Build the per-consequence namespace array from the engine's local values,
/// applying the amino-acid-change validity rules that preserve missense gating:
/// `ref_aa`/`alt_aa` are set only for a clean single-residue `X/Y` `amino_acids`;
/// `Protein_position` is passed through only when non-empty and not a range.
#[allow(clippy::too_many_arguments)]
pub fn build_attr_namespace<'a>(
    consequence: &'a str,
    gene: &'a str,
    feature_type: &'a str,
    feature: &'a str,
    biotype: &'a str,
    hgvsc: &'a str,
    hgvsp: &'a str,
    cdna_pos: &'a str,
    cds_pos: &'a str,
    protein_pos: &'a str,
    amino_acids: &'a str,
    codons: &'a str,
    ref_allele: &'a str,
    alt_allele: &'a str,
) -> [Option<&'a str>; ATTR_NAMES.len()] {
    let non_empty = |s: &'a str| if s.is_empty() { None } else { Some(s) };
    // ref_aa/alt_aa only for a clean single-residue X/Y substitution.
    let (ref_aa, alt_aa) = match amino_acids.split_once('/') {
        Some((r, a)) if r.len() == 1 && a.len() == 1 => (Some(r), Some(a)),
        _ => (None, None),
    };
    let protein = match non_empty(protein_pos) {
        Some(p) if !p.contains('-') => Some(p),
        _ => None,
    };
    [
        non_empty(consequence),
        non_empty(gene),
        non_empty(feature_type),
        non_empty(feature),
        non_empty(biotype),
        non_empty(hgvsc),
        non_empty(hgvsp),
        non_empty(cdna_pos),
        non_empty(cds_pos),
        protein,
        non_empty(amino_acids),
        non_empty(codons),
        ref_aa,
        alt_aa,
        non_empty(ref_allele),
        non_empty(alt_allele),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ns_missense() -> [Option<&'static str>; ATTR_NAMES.len()] {
        build_attr_namespace(
            "missense_variant",
            "ENSG1",
            "Transcript",
            "ENST1",
            "protein_coding",
            "",
            "",
            "",
            "",
            "320",
            "W/R",
            "TGG/CGG",
            "C",
            "T",
        )
    }

    #[test]
    fn compiles_and_evaluates_amino_acid_change() {
        let t = CompiledTemplate::compile("{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&ns_missense()).as_deref(), Some("W320R"));
    }

    #[test]
    fn feature_only_template() {
        let t = CompiledTemplate::compile("{Feature}").unwrap();
        assert_eq!(t.eval(&ns_missense()).as_deref(), Some("ENST1"));
    }

    #[test]
    fn literals_and_attrs_interleave() {
        let t = CompiledTemplate::compile("p.{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&ns_missense()).as_deref(), Some("p.W320R"));
    }

    #[test]
    fn non_missense_gates_to_none() {
        // synonymous: no aa change, no protein position
        let ns = build_attr_namespace(
            "synonymous_variant",
            "ENSG1",
            "Transcript",
            "ENST1",
            "protein_coding",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "C",
            "T",
        );
        let t = CompiledTemplate::compile("{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&ns), None);
    }

    #[test]
    fn range_position_and_multi_residue_gate() {
        let range = build_attr_namespace(
            "inframe", "", "", "", "", "", "", "", "", "550-551", "VV/A", "", "C", "T",
        );
        let t = CompiledTemplate::compile("{ref_aa}{Protein_position}{alt_aa}").unwrap();
        assert_eq!(t.eval(&range), None);
    }

    #[test]
    fn unknown_attribute_is_error() {
        assert!(CompiledTemplate::compile("{NoSuchAttr}").is_err());
    }

    #[test]
    fn unterminated_brace_is_error() {
        assert!(CompiledTemplate::compile("{ref_aa").is_err());
    }
}
