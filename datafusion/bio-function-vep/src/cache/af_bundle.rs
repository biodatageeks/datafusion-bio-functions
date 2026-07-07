//! AF-column bundling for the Parquet variation cache.
//!
//! Physically, the 27 per-population allele-frequency string columns are stored as 3
//! concatenated scalar `Utf8` columns (one per population family): each row's members are
//! joined by `|` (verified collision-free), an absent member is an empty field, and a row
//! with every member absent is null. Encoded miniblock+zstd via the shared string preset,
//! so a point-take reads ~3.6x fewer bytes than the original separate-27 layout and the AF
//! columns are ~3.8x smaller on disk, all on stock Parquet (no fullzip / no List levels).
//!
//! `bundle_af_columns` (build side, concatenate) and `unbundle_af_columns` (read side, split)
//! are exact inverses at the column-value level, so the downstream `af_values: Vec<String>`
//! interface — indexed in AF-group order — is unchanged. An absent member or a null parent
//! both materialize back to an empty string, matching the existing null→"" convention
//! (verified by the e2e parity gate: all AF fields 100%).

use std::collections::HashSet;
use std::sync::Arc;

use datafusion::arrow::array::{Array, ArrayRef, RecordBatch, StringArray, StringBuilder};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};

/// The 3 bundled list columns and, for each, the source AF column names in list order.
/// The flattened order MUST equal the AF slice of `VARIATION_REQUIRED_COLUMNS` / `AF_COL_NAMES`.
pub const AF_GROUPS: &[(&str, &[&str])] = &[
    ("af_global", &["AF", "AFR", "AMR", "EAS", "EUR", "SAS"]),
    (
        "af_gnomade",
        &[
            "gnomADe",
            "gnomADe_AFR",
            "gnomADe_AMR",
            "gnomADe_ASJ",
            "gnomADe_EAS",
            "gnomADe_FIN",
            "gnomADe_NFE",
            "gnomADe_SAS",
            "gnomADe_MID",
            "gnomADe_REMAINING",
        ],
    ),
    (
        "af_gnomadg",
        &[
            "gnomADg",
            "gnomADg_AFR",
            "gnomADg_AMI",
            "gnomADg_AMR",
            "gnomADg_ASJ",
            "gnomADg_EAS",
            "gnomADg_FIN",
            "gnomADg_MID",
            "gnomADg_NFE",
            "gnomADg_SAS",
            "gnomADg_REMAINING",
        ],
    ),
];

/// All 27 AF column names in canonical (flattened group) order.
pub fn af_column_order() -> Vec<&'static str> {
    AF_GROUPS
        .iter()
        .flat_map(|(_, cols)| cols.iter().copied())
        .collect()
}

/// The 3 bundled column names.
pub fn af_group_names() -> Vec<&'static str> {
    AF_GROUPS.iter().map(|(name, _)| *name).collect()
}

// `|` is verified collision-free in AF values (0 of 613M; `:`=allele/freq delimiter, `,`=8%).
// The build-time check in `concat_group` guards against future source changes.
const AF_CONCAT_SEP: char = '|';

/// Field for a concatenated AF group: scalar `Utf8`, encoded like the other
/// string columns.
fn concat_field(name: &str) -> Field {
    Field::new(name, DataType::Utf8, true)
}

/// Concatenate one group's members into a single Utf8 column (positional, `|`-joined,
/// null when all members absent). Errors if any member value contains the separator.
pub fn concat_group(members: &[&StringArray]) -> Result<StringArray> {
    let n = members[0].len();
    let width = members.len();
    let mut b = StringBuilder::with_capacity(n, n * 16);
    let mut s = String::new();
    for r in 0..n {
        s.clear();
        let mut all_absent = true;
        for (c, m) in members.iter().enumerate().take(width) {
            if c > 0 {
                s.push(AF_CONCAT_SEP);
            }
            if !m.is_null(r) {
                let v = m.value(r);
                if !v.is_empty() {
                    if v.contains(AF_CONCAT_SEP) {
                        return Err(DataFusionError::Execution(format!(
                            "AF value contains the '{AF_CONCAT_SEP}' bundle separator: {v:?}"
                        )));
                    }
                    all_absent = false;
                    s.push_str(v);
                }
            }
        }
        if all_absent {
            b.append_null();
        } else {
            b.append_value(&s);
        }
    }
    Ok(b.finish())
}

/// Split one bundled concat column back into its `width` member `Utf8` columns (positional;
/// absent member -> "", null parent -> "" for every member, matching downstream expectations).
pub fn split_group(col: &StringArray, width: usize) -> Vec<StringArray> {
    let n = col.len();
    let mut builders: Vec<StringBuilder> = (0..width)
        .map(|_| StringBuilder::with_capacity(n, n * 4))
        .collect();
    for r in 0..n {
        if col.is_null(r) {
            for b in builders.iter_mut() {
                b.append_value("");
            }
        } else {
            let mut parts = col.value(r).split(AF_CONCAT_SEP);
            for b in builders.iter_mut() {
                b.append_value(parts.next().unwrap_or(""));
            }
        }
    }
    builders.into_iter().map(|mut b| b.finish()).collect()
}

/// Schema with the 27 AF `Utf8` fields replaced by 3 concatenated scalar `Utf8` group fields
/// (appended after the non-AF fields, preserving non-AF order + schema metadata).
pub fn bundle_schema(schema: &Schema) -> Schema {
    let af: HashSet<&str> = af_column_order().into_iter().collect();
    let mut fields: Vec<Arc<Field>> = schema
        .fields()
        .iter()
        .filter(|f| !af.contains(f.name().as_str()))
        .cloned()
        .collect();
    for (name, _) in AF_GROUPS {
        fields.push(Arc::new(concat_field(name)));
    }
    Schema::new_with_metadata(fields, schema.metadata().clone())
}

fn string_col<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a StringArray> {
    let idx = batch
        .schema_ref()
        .index_of(name)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| DataFusionError::Execution(format!("AF column '{name}' must be Utf8")))
}

/// Build side: replace the 27 AF `Utf8` columns with 3 concatenated scalar `Utf8` columns.
pub fn bundle_af_columns(batch: &RecordBatch) -> Result<RecordBatch> {
    let af: HashSet<&str> = af_column_order().into_iter().collect();
    let schema = batch.schema();
    let mut cols: Vec<ArrayRef> = Vec::new();
    for (i, f) in schema.fields().iter().enumerate() {
        if !af.contains(f.name().as_str()) {
            cols.push(batch.column(i).clone());
        }
    }
    for (_, members) in AF_GROUPS {
        let arrays = members
            .iter()
            .map(|m| string_col(batch, m))
            .collect::<Result<Vec<_>>>()?;
        cols.push(Arc::new(concat_group(&arrays)?) as ArrayRef);
    }
    let out_schema = Arc::new(bundle_schema(schema.as_ref()));
    RecordBatch::try_new(out_schema, cols)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}

/// Read side: expand any bundled `List<Utf8>` group column back into its member
/// `Utf8` columns. Pass-through for non-group columns; robust to 0..3 groups present.
pub fn unbundle_af_columns(batch: &RecordBatch) -> Result<RecordBatch> {
    let schema = batch.schema();
    let mut fields: Vec<Arc<Field>> = Vec::new();
    let mut cols: Vec<ArrayRef> = Vec::new();
    for (i, f) in schema.fields().iter().enumerate() {
        if let Some((_, members)) = AF_GROUPS.iter().find(|(n, _)| *n == f.name()) {
            let concat = batch
                .column(i)
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "bundled column '{}' must be Utf8",
                        f.name()
                    ))
                })?;
            for (m, arr) in members.iter().zip(split_group(concat, members.len())) {
                fields.push(Arc::new(Field::new(*m, DataType::Utf8, true)));
                cols.push(Arc::new(arr) as ArrayRef);
            }
        } else {
            fields.push(f.clone());
            cols.push(batch.column(i).clone());
        }
    }
    RecordBatch::try_new(
        Arc::new(Schema::new_with_metadata(fields, schema.metadata().clone())),
        cols,
    )
    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}

/// Read-side projection rewrite: if any of the 27 AF columns is requested, drop them
/// and request the 3 bundled group columns instead (they expand back on read).
pub fn bundle_projection(projection: &[String]) -> Vec<String> {
    let af: HashSet<&str> = af_column_order().into_iter().collect();
    let mut out: Vec<String> = Vec::with_capacity(projection.len());
    let mut any_af = false;
    for c in projection {
        if af.contains(c.as_str()) {
            any_af = true;
        } else {
            out.push(c.clone());
        }
    }
    if any_af {
        for g in af_group_names() {
            if !out.iter().any(|c| c == g) {
                out.push(g.to_string());
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sa(vals: &[Option<&str>]) -> StringArray {
        StringArray::from(vals.to_vec())
    }

    #[test]
    fn af_groups_flatten_to_27_columns() {
        let order = af_column_order();
        assert_eq!(order.len(), 27);
        assert_eq!(order[0], "AF");
        assert_eq!(order[6], "gnomADe");
        assert_eq!(order[16], "gnomADg");
    }

    #[test]
    fn concat_split_round_trips_values() {
        // 6-member group, 4 rows: full, partial, all-empty, all-null
        let members_owned = vec![
            sa(&[Some("A:0.1"), Some("x"), Some(""), None]),
            sa(&[Some("A:0.2"), None, Some(""), None]),
            sa(&[Some("A:0.3"), Some(""), Some(""), None]),
            sa(&[Some("A:0.4"), Some("y"), None, None]),
            sa(&[Some("A:0.5"), Some(""), Some(""), None]),
            sa(&[Some("A:0.6"), Some("z"), Some(""), None]),
        ];
        let members: Vec<&StringArray> = members_owned.iter().collect();
        let concat = concat_group(&members).unwrap();
        assert!(!concat.is_null(0));
        assert!(!concat.is_null(1));
        assert!(concat.is_null(2)); // all-empty -> null
        assert!(concat.is_null(3)); // all-null  -> null

        let back = split_group(&concat, 6);
        assert_eq!(back.len(), 6);
        let expect: Vec<Vec<&str>> = vec![
            vec!["A:0.1", "A:0.2", "A:0.3", "A:0.4", "A:0.5", "A:0.6"],
            vec!["x", "", "", "y", "", "z"],
            vec!["", "", "", "", "", ""],
            vec!["", "", "", "", "", ""],
        ];
        for row in 0..4 {
            for col in 0..6 {
                assert_eq!(
                    back[col].value(row),
                    expect[row][col],
                    "row {row} col {col}"
                );
            }
        }
    }

    #[test]
    fn concat_group_rejects_separator_in_value() {
        let members_owned = vec![sa(&[Some("a|b")]), sa(&[Some("c")])];
        let members: Vec<&StringArray> = members_owned.iter().collect();
        assert!(concat_group(&members).is_err());
    }

    #[test]
    fn unbundle_af_columns_splits_concat() {
        let chrom = sa(&[Some("chr1"), Some("chr1")]);
        let af_global = StringArray::from(vec![Some("A:0.1|||||"), None]);
        let af_gnomade = StringArray::from(vec![Some("g|||||||||"), None]);
        let af_gnomadg = StringArray::from(vec![Some("h||||||||||"), None]);
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("af_global", DataType::Utf8, true),
            Field::new("af_gnomade", DataType::Utf8, true),
            Field::new("af_gnomadg", DataType::Utf8, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(chrom),
                Arc::new(af_global),
                Arc::new(af_gnomade),
                Arc::new(af_gnomadg),
            ],
        )
        .unwrap();
        let out = unbundle_af_columns(&batch).unwrap();
        assert_eq!(out.num_columns(), 28); // 1 non-AF + 27 AF
        let af_idx = out.schema().index_of("AF").unwrap();
        let af = out
            .column(af_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(af.value(0), "A:0.1");
        assert_eq!(af.value(1), ""); // null parent -> ""
    }
}
