//! DataFusion scalar functions over the [`crate::kernels`].
//!
//! | Function | Arguments after the sequence | Returns |
//! |---|---|---|
//! | `seq_revcomp` | `alphabet = 'dna'`, `strict = false` | same string type |
//! | `seq_reverse` | — | same string type |
//! | `seq_complement` | `alphabet = 'dna'`, `strict = false` | same string type |
//! | `seq_gc_content` | `denominator = 'length'` | `Float64` |
//! | `seq_translate` | `frame = 1`, `table = 1`, `to_stop = false`, `stop_symbol = '*'`, `cds = false`, `ambiguous = 'x'` | same string type |
//! | `seq_cds_status` | `frame = 1`, `table = 1` | `Utf8` |
//! | `seq_transcribe`, `seq_back_transcribe` | — | same string type |
//! | `seq_mean_quality` | `offset = 33`, `method = 'error_prob'` | `Float64` |
//! | `seq_expected_errors` | `offset = 33` | `Float64` |
//! | `seq_min_quality`, `seq_max_quality` | `offset = 33` | `Int32` |
//! | `seq_orfs` | `table = 1`, `min_len = 30`, `strands = 'both'`, `start_codons = 'table'`, `allow_partial = false`, `include_protein = false`, `zero_based = false` | `List<Struct>` |
//! | `seq_find_motif` | `pattern`, `both_strands = false`, `overlapping = true`, `zero_based = false` | `List<Int64>` |
//! | `seq_gc_content_windows` | `window`, `step = window`, `denominator = 'length'` | `List<Float64>` |
//! | `seq_gc_skew` | `window = 100`, `step = window` | `List<Float64>` |
//! | `seq_molecular_weight` | `seq_type = 'DNA'`, `double_stranded = false`, `circular = false`, `monoisotopic = false` | `Float64` |
//! | `seq_hash` | `canonical = false` | `UInt64` |
//! | `seq_trim_quality` | `quality` (a column), `leading`, `trailing`, `window_size`, `window_quality`, `offset = 33` | `Struct<sequence, quality_scores>` |
//!
//! The sequence (and the quality of `seq_trim_quality`) may be `Utf8`,
//! `LargeUtf8` or `Utf8View`. Later arguments are positional and must be
//! constants; all are optional except `pattern` and the `window` of
//! `seq_gc_content_windows`. A null input gives null, and so does a ratio
//! with a zero denominator. The names carry a `seq_` prefix because
//! DataFusion already has a built-in `translate`.
//!
//! `step`, `leading`, `trailing`, `window_size` and `window_quality` accept
//! `NULL` to mean "not set". `seq_trim_quality` needs at least one step.
//!
//! `seq_orfs` and `seq_find_motif` report 1-based closed coordinates unless
//! `zero_based` is true, which gives 0-based half-open ones. `seq_orfs`
//! returns a list of `{start, end, strand, frame, aa_len, has_start,
//! has_stop}` structs, with a `protein` field when `include_protein` is the
//! literal `true`.
//!
//! Lengths count bytes, so the `length` denominator of `seq_gc_content`
//! counts a non-ASCII character as several bytes. Sequence data is ASCII in
//! practice.

use std::fmt;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, AsArray, BooleanBuilder, Float64Builder, GenericStringBuilder, Int32Builder,
    Int64Builder, LargeStringArray, ListArray, ListBuilder, StringArray, StringBuilder,
    StringViewArray, StringViewBuilder, StructArray, UInt64Builder,
};
use datafusion::arrow::buffer::{NullBuffer, OffsetBuffer};
use datafusion::arrow::datatypes::{DataType, Field, FieldRef, Fields};
use datafusion::common::{DataFusionError, Result, ScalarValue, exec_err, internal_err, plan_err};
use datafusion::logical_expr::{
    ColumnarValue, ReturnFieldArgs, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature,
    Volatility,
};
use datafusion::prelude::SessionContext;

use crate::kernels::{
    Alphabet, Ambiguity, CdsError, GcDenominator, KernelError, MolecularWeight, Motif, OrfFinder,
    QualityMethod, QualityTrimmer, StartCodons, Strands, Translator, Windows, back_transcribe_into,
    complement_into, expected_errors, gc_content, gc_content_windows_into, gc_skew_into,
    phred_mean, phred_min_max, revcomp_into, reverse_into, seq_hash, transcribe_into,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Param {
    Text,
    Int,
    Bool,
    Float,
}

impl fmt::Display for Param {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Param::Text => "a string",
            Param::Int => "an integer",
            Param::Bool => "a boolean",
            Param::Float => "a number",
        })
    }
}

/// A constant argument: its name, its type, and whether `NULL` means "not
/// set" (otherwise `NULL` is an error).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Arg {
    name: &'static str,
    kind: Param,
    nullable: bool,
}

const fn arg(name: &'static str, kind: Param) -> Arg {
    Arg {
        name,
        kind,
        nullable: false,
    }
}

const fn nullable(name: &'static str, kind: Param) -> Arg {
    Arg {
        name,
        kind,
        nullable: true,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Kind {
    Revcomp,
    Reverse,
    Complement,
    GcContent,
    Translate,
    CdsStatus,
    Transcribe,
    BackTranscribe,
    MeanQuality,
    ExpectedErrors,
    MinQuality,
    MaxQuality,
    Orfs,
    FindMotif,
    GcContentWindows,
    GcSkew,
    MolecularWeight,
    Hash,
    TrimQuality,
}

/// Position of `include_protein` among the constant arguments of `seq_orfs`.
const ORFS_INCLUDE_PROTEIN: usize = 5;

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Revcomp => "seq_revcomp",
            Kind::Reverse => "seq_reverse",
            Kind::Complement => "seq_complement",
            Kind::GcContent => "seq_gc_content",
            Kind::Translate => "seq_translate",
            Kind::CdsStatus => "seq_cds_status",
            Kind::Transcribe => "seq_transcribe",
            Kind::BackTranscribe => "seq_back_transcribe",
            Kind::MeanQuality => "seq_mean_quality",
            Kind::ExpectedErrors => "seq_expected_errors",
            Kind::MinQuality => "seq_min_quality",
            Kind::MaxQuality => "seq_max_quality",
            Kind::Orfs => "seq_orfs",
            Kind::FindMotif => "seq_find_motif",
            Kind::GcContentWindows => "seq_gc_content_windows",
            Kind::GcSkew => "seq_gc_skew",
            Kind::MolecularWeight => "seq_molecular_weight",
            Kind::Hash => "seq_hash",
            Kind::TrimQuality => "seq_trim_quality",
        }
    }

    /// String columns at the start of the argument list.
    fn columns(self) -> usize {
        match self {
            Kind::TrimQuality => 2,
            _ => 1,
        }
    }

    /// Constant arguments after the columns, in order.
    fn params(self) -> &'static [Arg] {
        use Param::{Bool, Float, Int, Text};
        const NONE: &[Arg] = &[];
        const COMPLEMENT: &[Arg] = &[arg("alphabet", Text), arg("strict", Bool)];
        const GC_CONTENT: &[Arg] = &[arg("denominator", Text)];
        const TRANSLATE: &[Arg] = &[
            arg("frame", Int),
            arg("table", Int),
            arg("to_stop", Bool),
            arg("stop_symbol", Text),
            arg("cds", Bool),
            arg("ambiguous", Text),
        ];
        const CDS_STATUS: &[Arg] = &[arg("frame", Int), arg("table", Int)];
        const MEAN_QUALITY: &[Arg] = &[arg("offset", Int), arg("method", Text)];
        const OFFSET: &[Arg] = &[arg("offset", Int)];
        const ORFS: &[Arg] = &[
            arg("table", Int),
            arg("min_len", Int),
            arg("strands", Text),
            arg("start_codons", Text),
            arg("allow_partial", Bool),
            arg("include_protein", Bool),
            arg("zero_based", Bool),
        ];
        const FIND_MOTIF: &[Arg] = &[
            arg("pattern", Text),
            arg("both_strands", Bool),
            arg("overlapping", Bool),
            arg("zero_based", Bool),
        ];
        const GC_CONTENT_WINDOWS: &[Arg] = &[
            arg("window", Int),
            nullable("step", Int),
            arg("denominator", Text),
        ];
        const GC_SKEW: &[Arg] = &[arg("window", Int), nullable("step", Int)];
        const MOLECULAR_WEIGHT: &[Arg] = &[
            arg("seq_type", Text),
            arg("double_stranded", Bool),
            arg("circular", Bool),
            arg("monoisotopic", Bool),
        ];
        const HASH: &[Arg] = &[arg("canonical", Bool)];
        const TRIM_QUALITY: &[Arg] = &[
            nullable("leading", Int),
            nullable("trailing", Int),
            nullable("window_size", Int),
            nullable("window_quality", Float),
            arg("offset", Int),
        ];
        match self {
            Kind::Revcomp | Kind::Complement => COMPLEMENT,
            Kind::Reverse | Kind::Transcribe | Kind::BackTranscribe => NONE,
            Kind::GcContent => GC_CONTENT,
            Kind::Translate => TRANSLATE,
            Kind::CdsStatus => CDS_STATUS,
            Kind::MeanQuality => MEAN_QUALITY,
            Kind::ExpectedErrors | Kind::MinQuality | Kind::MaxQuality => OFFSET,
            Kind::Orfs => ORFS,
            Kind::FindMotif => FIND_MOTIF,
            Kind::GcContentWindows => GC_CONTENT_WINDOWS,
            Kind::GcSkew => GC_SKEW,
            Kind::MolecularWeight => MOLECULAR_WEIGHT,
            Kind::Hash => HASH,
            Kind::TrimQuality => TRIM_QUALITY,
        }
    }

    /// How many of [`Kind::params`] must be given.
    fn required(self) -> usize {
        match self {
            Kind::FindMotif | Kind::GcContentWindows => 1,
            _ => 0,
        }
    }
}

/// The string type a function returns for a string input of type `input`.
fn same_string_type(input: Option<&DataType>) -> DataType {
    match input {
        Some(t @ (DataType::LargeUtf8 | DataType::Utf8View)) => t.clone(),
        _ => DataType::Utf8,
    }
}

/// The fields of one ORF, as `seq_orfs` returns them.
fn orf_fields(include_protein: bool) -> Fields {
    let mut fields = vec![
        Field::new("start", DataType::Int64, true),
        Field::new("end", DataType::Int64, true),
        Field::new("strand", DataType::Utf8, true),
        Field::new("frame", DataType::Int32, true),
        Field::new("aa_len", DataType::Int64, true),
        Field::new("has_start", DataType::Boolean, true),
        Field::new("has_stop", DataType::Boolean, true),
    ];
    if include_protein {
        fields.push(Field::new("protein", DataType::Utf8, true));
    }
    Fields::from(fields)
}

/// The fields of the struct `seq_trim_quality` returns.
fn trim_fields(sequence: DataType, quality: DataType) -> Fields {
    Fields::from(vec![
        Field::new("sequence", sequence, true),
        Field::new("quality_scores", quality, true),
    ])
}

/// One of the `seq_*` scalar functions.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct SeqUdf {
    kind: Kind,
    signature: Signature,
}

impl SeqUdf {
    fn new(kind: Kind) -> Self {
        Self {
            kind,
            signature: Signature::user_defined(Volatility::Immutable),
        }
    }

    fn output_type(&self, column_types: &[DataType], include_protein: bool) -> DataType {
        match self.kind {
            Kind::Revcomp
            | Kind::Reverse
            | Kind::Complement
            | Kind::Translate
            | Kind::Transcribe
            | Kind::BackTranscribe => same_string_type(column_types.first()),
            Kind::GcContent | Kind::MeanQuality | Kind::ExpectedErrors | Kind::MolecularWeight => {
                DataType::Float64
            }
            Kind::CdsStatus => DataType::Utf8,
            Kind::MinQuality | Kind::MaxQuality => DataType::Int32,
            Kind::Hash => DataType::UInt64,
            Kind::FindMotif => DataType::new_list(DataType::Int64, true),
            Kind::GcContentWindows | Kind::GcSkew => DataType::new_list(DataType::Float64, true),
            Kind::Orfs => DataType::new_list(DataType::Struct(orf_fields(include_protein)), true),
            Kind::TrimQuality => DataType::Struct(trim_fields(
                same_string_type(column_types.first()),
                same_string_type(column_types.get(1)),
            )),
        }
    }
}

impl ScalarUDFImpl for SeqUdf {
    fn name(&self) -> &str {
        self.kind.name()
    }

    fn signature(&self) -> &Signature {
        &self.signature
    }

    fn coerce_types(&self, arg_types: &[DataType]) -> Result<Vec<DataType>> {
        let name = self.kind.name();
        let (columns, params) = (self.kind.columns(), self.kind.params());
        let (min, max) = (columns + self.kind.required(), columns + params.len());
        if arg_types.len() < min || arg_types.len() > max {
            return plan_err!(
                "{name} takes {min} to {max} arguments, got {}",
                arg_types.len()
            );
        }
        let mut coerced = Vec::with_capacity(arg_types.len());
        for (i, t) in arg_types[..columns].iter().enumerate() {
            coerced.push(match t {
                t @ (DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View) => t.clone(),
                DataType::Null => DataType::Utf8,
                other => {
                    let position = ["first", "second"][i];
                    return plan_err!(
                        "{name}: the {position} argument must be a string (Utf8, LargeUtf8 or Utf8View), got {other}"
                    );
                }
            });
        }
        for (
            t,
            Arg {
                name: param, kind, ..
            },
        ) in arg_types[columns..].iter().zip(params)
        {
            coerced.push(match (kind, t) {
                (
                    Param::Text,
                    DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View | DataType::Null,
                ) => DataType::Utf8,
                (Param::Int, t) if t.is_integer() || t.is_null() => DataType::Int64,
                (Param::Float, t) if t.is_numeric() || t.is_null() => DataType::Float64,
                (Param::Bool, DataType::Boolean | DataType::Null) => DataType::Boolean,
                _ => return plan_err!("{name}: {param} must be {kind}, got {t}"),
            });
        }
        Ok(coerced)
    }

    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        if self.kind == Kind::Orfs {
            // The type depends on the value of `include_protein`.
            return internal_err!("seq_orfs: return_field_from_args must be used");
        }
        Ok(self.output_type(arg_types, false))
    }

    fn return_field_from_args(&self, args: ReturnFieldArgs) -> Result<FieldRef> {
        let types: Vec<DataType> = args
            .arg_fields
            .iter()
            .map(|field| field.data_type().clone())
            .collect();
        let include_protein = self.kind == Kind::Orfs
            && matches!(
                args.scalar_arguments.get(1 + ORFS_INCLUDE_PROTEIN),
                Some(Some(ScalarValue::Boolean(Some(true))))
            );
        Ok(Arc::new(Field::new(
            self.kind.name(),
            self.output_type(&types, include_protein),
            true,
        )))
    }

    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let name = self.kind.name();
        let input = args.args[0].to_array(args.number_rows)?;
        let opts = Options::new(self.kind, &args.args[self.kind.columns()..])?;
        let ke = |e: KernelError| DataFusionError::Execution(format!("{name}: {e}"));
        let output = match self.kind {
            Kind::Reverse => map_strings(name, &input, OutputRatio::Same, |s, out| {
                reverse_into(s, out);
                Ok(())
            })?,
            Kind::Revcomp => {
                let alphabet = Alphabet::parse(opts.text(0, "dna")?).map_err(ke)?;
                let strict = opts.bool(1, false)?;
                map_strings(name, &input, OutputRatio::Same, |s, out| {
                    revcomp_into(s, alphabet, strict, out).map_err(ke)
                })?
            }
            Kind::Complement => {
                let alphabet = Alphabet::parse(opts.text(0, "dna")?).map_err(ke)?;
                let strict = opts.bool(1, false)?;
                map_strings(name, &input, OutputRatio::Same, |s, out| {
                    complement_into(s, alphabet, strict, out).map_err(ke)
                })?
            }
            Kind::Transcribe => map_strings(name, &input, OutputRatio::Same, |s, out| {
                transcribe_into(s, out);
                Ok(())
            })?,
            Kind::BackTranscribe => map_strings(name, &input, OutputRatio::Same, |s, out| {
                back_transcribe_into(s, out);
                Ok(())
            })?,
            Kind::Translate => {
                let ambiguity = Ambiguity::parse(opts.text(5, "x")?).map_err(ke)?;
                let mut translator =
                    Translator::new(opts.int(0, 1)?, opts.int(1, 1)?, opts.bool(2, false)?)
                        .map_err(ke)?
                        .with_stop_symbol(opts.text(3, "*")?)
                        .with_cds(opts.bool(4, false)?)
                        .with_ambiguity(ambiguity);
                map_strings(name, &input, OutputRatio::OneThird, |s, out| {
                    translator.translate_into(s, out).map_err(ke)
                })?
            }
            Kind::CdsStatus => {
                let mut translator =
                    Translator::new(opts.int(0, 1)?, opts.int(1, 1)?, false).map_err(ke)?;
                let mut builder = StringBuilder::with_capacity(input.len(), 0);
                for_each_value(name, &input, |value| {
                    let status = value.and_then(|s| translator.cds_status(s));
                    builder.append_option(status.map(CdsError::as_str));
                    Ok(())
                })?;
                Arc::new(builder.finish())
            }
            Kind::GcContent => {
                let denominator = GcDenominator::parse(opts.text(0, "length")?).map_err(ke)?;
                map_floats(name, &input, |s| Ok(gc_content(s, denominator)))?
            }
            Kind::MeanQuality => {
                let offset = opts.offset(0)?;
                let method = QualityMethod::parse(opts.text(1, "error_prob")?).map_err(ke)?;
                map_floats(name, &input, |q| phred_mean(q, offset, method).map_err(ke))?
            }
            Kind::ExpectedErrors => {
                let offset = opts.offset(0)?;
                map_floats(name, &input, |q| {
                    expected_errors(q, offset).map(Some).map_err(ke)
                })?
            }
            Kind::MinQuality | Kind::MaxQuality => {
                let offset = opts.offset(0)?;
                let max = self.kind == Kind::MaxQuality;
                let mut builder = Int32Builder::with_capacity(input.len());
                for_each_value(name, &input, |value| {
                    let range = match value {
                        Some(q) => phred_min_max(q, offset).map_err(ke)?,
                        None => None,
                    };
                    builder
                        .append_option(range.map(|(lowest, highest)| {
                            i32::from(if max { highest } else { lowest })
                        }));
                    Ok(())
                })?;
                Arc::new(builder.finish())
            }
            Kind::MolecularWeight => {
                let weigher = MolecularWeight::new(
                    opts.text(0, "DNA")?,
                    opts.bool(1, false)?,
                    opts.bool(2, false)?,
                    opts.bool(3, false)?,
                )
                .map_err(ke)?;
                map_floats(name, &input, |s| Ok(weigher.weight(s)))?
            }
            Kind::Hash => {
                let canonical = opts.bool(0, false)?;
                let mut scratch = Vec::new();
                let mut builder = UInt64Builder::with_capacity(input.len());
                for_each_value(name, &input, |value| {
                    builder.append_option(value.map(|s| seq_hash(s, canonical, &mut scratch)));
                    Ok(())
                })?;
                Arc::new(builder.finish())
            }
            Kind::FindMotif => {
                let motif =
                    Motif::new(opts.text(0, "")?, opts.bool(1, false)?, opts.bool(2, true)?)
                        .map_err(ke)?;
                let first = first_position(opts.bool(3, false)?);
                let mut hits = Vec::new();
                let mut builder = ListBuilder::new(Int64Builder::new());
                for_each_value(name, &input, |value| {
                    if let Some(s) = value {
                        motif.find_into(s, &mut hits);
                        for &hit in &hits {
                            builder.values().append_value(hit as i64 + first);
                        }
                    }
                    builder.append(value.is_some());
                    Ok(())
                })?;
                Arc::new(builder.finish())
            }
            Kind::GcContentWindows => {
                let windows = Windows::new(opts.int(0, 0)?, opts.opt_int(1)?).map_err(ke)?;
                let denominator = GcDenominator::parse(opts.text(2, "length")?).map_err(ke)?;
                map_float_lists(name, &input, |s, out| {
                    gc_content_windows_into(s, windows, denominator, out)
                })?
            }
            Kind::GcSkew => {
                let windows = Windows::new(opts.int(0, 100)?, opts.opt_int(1)?).map_err(ke)?;
                let mut skews = Vec::new();
                map_float_lists(name, &input, |s, out| {
                    skews.clear();
                    gc_skew_into(s, windows, &mut skews);
                    out.extend(skews.iter().copied().map(Some));
                })?
            }
            Kind::Orfs => {
                let include_protein = opts.bool(ORFS_INCLUDE_PROTEIN, false)?;
                let declared = matches!(
                    args.return_field.data_type(),
                    DataType::List(item) if matches!(
                        item.data_type(),
                        DataType::Struct(fields) if fields.len() == orf_fields(true).len()
                    )
                );
                if include_protein != declared {
                    return exec_err!("{name}: include_protein must be a literal true or false");
                }
                let mut finder = OrfFinder::new(
                    opts.int(0, 1)?,
                    opts.int(1, 30)?,
                    Strands::parse(opts.text(2, "both")?).map_err(ke)?,
                    StartCodons::parse(opts.text(3, "table")?).map_err(ke)?,
                    opts.bool(4, false)?,
                    include_protein,
                )
                .map_err(ke)?;
                orf_lists(
                    name,
                    &input,
                    &mut finder,
                    first_position(opts.bool(6, false)?),
                    include_protein,
                )?
            }
            Kind::TrimQuality => {
                let quality = args.args[1].to_array(args.number_rows)?;
                let window = match (opts.opt_int(2)?, opts.opt_float(3)?) {
                    (Some(size), Some(quality)) => Some((size, quality)),
                    (None, None) => None,
                    _ => {
                        return exec_err!(
                            "{name}: window_size and window_quality must be given together"
                        );
                    }
                };
                let trimmer = QualityTrimmer::new(
                    opts.opt_int(0)?,
                    opts.opt_int(1)?,
                    window,
                    opts.offset(4)?,
                )
                .map_err(ke)?;
                trim_structs(name, &input, &quality, &trimmer)?
            }
        };
        Ok(ColumnarValue::Array(output))
    }
}

/// What is added to a 0-based start to report it: 0 for 0-based half-open
/// coordinates, 1 for 1-based closed ones (whose end equals the 0-based
/// exclusive end).
fn first_position(zero_based: bool) -> i64 {
    if zero_based { 0 } else { 1 }
}

/// The constant option arguments of one call.
struct Options<'a> {
    kind: Kind,
    args: Vec<&'a ScalarValue>,
}

impl<'a> Options<'a> {
    fn new(kind: Kind, args: &'a [ColumnarValue]) -> Result<Self> {
        let params = kind.params();
        let mut scalars = Vec::with_capacity(args.len());
        for (arg, Arg { name: param, .. }) in args.iter().zip(params) {
            match arg {
                ColumnarValue::Scalar(s) => scalars.push(s),
                ColumnarValue::Array(_) => {
                    return exec_err!("{}: {param} must be a constant", kind.name());
                }
            }
        }
        Ok(Self {
            kind,
            args: scalars,
        })
    }

    fn null_err<T>(&self, i: usize) -> Result<T> {
        exec_err!(
            "{}: {} must not be null",
            self.kind.name(),
            self.kind.params()[i].name
        )
    }

    fn text(&self, i: usize, default: &'static str) -> Result<&str> {
        match self.args.get(i) {
            None => Ok(default),
            Some(ScalarValue::Utf8(Some(s))) => Ok(s),
            Some(_) => self.null_err(i),
        }
    }

    fn int(&self, i: usize, default: i64) -> Result<i64> {
        match self.args.get(i) {
            None => Ok(default),
            Some(ScalarValue::Int64(Some(v))) => Ok(*v),
            Some(_) => self.null_err(i),
        }
    }

    /// An integer argument where a missing value or `NULL` means "not set".
    fn opt_int(&self, i: usize) -> Result<Option<i64>> {
        debug_assert!(self.kind.params()[i].nullable);
        Ok(match self.args.get(i) {
            Some(ScalarValue::Int64(v)) => *v,
            _ => None,
        })
    }

    /// A numeric argument where a missing value or `NULL` means "not set".
    fn opt_float(&self, i: usize) -> Result<Option<f64>> {
        debug_assert!(self.kind.params()[i].nullable);
        Ok(match self.args.get(i) {
            Some(ScalarValue::Float64(v)) => *v,
            _ => None,
        })
    }

    fn bool(&self, i: usize, default: bool) -> Result<bool> {
        match self.args.get(i) {
            None => Ok(default),
            Some(ScalarValue::Boolean(Some(v))) => Ok(*v),
            Some(_) => self.null_err(i),
        }
    }

    /// The Phred offset argument at `i` (default 33).
    fn offset(&self, i: usize) -> Result<u8> {
        let offset = self.int(i, 33)?;
        u8::try_from(offset).or_else(|_| {
            exec_err!(
                "{}: offset must be between 0 and 255, got {offset}",
                self.kind.name()
            )
        })
    }
}

/// A string array of any of the three Arrow string types, read by index.
enum Strings<'a> {
    Utf8(&'a StringArray),
    Large(&'a LargeStringArray),
    View(&'a StringViewArray),
}

impl<'a> Strings<'a> {
    fn new(name: &str, array: &'a ArrayRef) -> Result<Self> {
        match array.data_type() {
            DataType::Utf8 => Ok(Strings::Utf8(array.as_string::<i32>())),
            DataType::LargeUtf8 => Ok(Strings::Large(array.as_string::<i64>())),
            DataType::Utf8View => Ok(Strings::View(array.as_string_view())),
            other => exec_err!("{name}: expected a string array, got {other}"),
        }
    }

    fn get(&self, i: usize) -> Option<&'a str> {
        match self {
            Strings::Utf8(a) => a.is_valid(i).then(|| a.value(i)),
            Strings::Large(a) => a.is_valid(i).then(|| a.value(i)),
            Strings::View(a) => a.is_valid(i).then(|| a.value(i)),
        }
    }
}

/// Calls `f` with every value of a string array, `None` for nulls.
fn for_each_value<F>(name: &str, input: &ArrayRef, mut f: F) -> Result<()>
where
    F: FnMut(Option<&[u8]>) -> Result<()>,
{
    match input.data_type() {
        DataType::Utf8 => input
            .as_string::<i32>()
            .iter()
            .try_for_each(|v| f(v.map(str::as_bytes))),
        DataType::LargeUtf8 => input
            .as_string::<i64>()
            .iter()
            .try_for_each(|v| f(v.map(str::as_bytes))),
        DataType::Utf8View => input
            .as_string_view()
            .iter()
            .try_for_each(|v| f(v.map(str::as_bytes))),
        other => exec_err!("{name}: expected a string array, got {other}"),
    }
}

/// Expected output size relative to the input, for pre-sizing.
#[derive(Clone, Copy)]
enum OutputRatio {
    Same,
    OneThird,
}

/// A string builder of the same type as the input.
enum StringSink {
    Utf8(GenericStringBuilder<i32>),
    LargeUtf8(GenericStringBuilder<i64>),
    Utf8View(StringViewBuilder),
}

impl StringSink {
    /// `ratio` is the expected output bytes per input byte (1 for complement,
    /// 3 input bytes per output byte for translation).
    fn for_input(input: &ArrayRef, ratio: OutputRatio) -> Self {
        let rows = input.len();
        // Size from the offsets of this array's rows, not the whole value
        // buffer, which a sliced array shares with rows outside the slice.
        let hint = |bytes: usize| match ratio {
            OutputRatio::Same => bytes,
            OutputRatio::OneThird => bytes / 3,
        };
        match input.data_type() {
            DataType::LargeUtf8 => {
                let o = input.as_string::<i64>().value_offsets();
                let bytes = (o[o.len() - 1] - o[0]) as usize;
                StringSink::LargeUtf8(GenericStringBuilder::with_capacity(rows, hint(bytes)))
            }
            DataType::Utf8View => StringSink::Utf8View(StringViewBuilder::with_capacity(rows)),
            _ => {
                let bytes = input.as_string_opt::<i32>().map_or(0, |a| {
                    let o = a.value_offsets();
                    (o[o.len() - 1] - o[0]) as usize
                });
                StringSink::Utf8(GenericStringBuilder::with_capacity(rows, hint(bytes)))
            }
        }
    }

    fn push(&mut self, value: Option<&str>) {
        match self {
            StringSink::Utf8(b) => b.append_option(value),
            StringSink::LargeUtf8(b) => b.append_option(value),
            StringSink::Utf8View(b) => b.append_option(value),
        }
    }

    fn finish(self) -> ArrayRef {
        match self {
            StringSink::Utf8(mut b) => Arc::new(b.finish()),
            StringSink::LargeUtf8(mut b) => Arc::new(b.finish()),
            StringSink::Utf8View(mut b) => Arc::new(b.finish()),
        }
    }
}

fn map_strings<F>(name: &str, input: &ArrayRef, ratio: OutputRatio, mut f: F) -> Result<ArrayRef>
where
    F: FnMut(&[u8], &mut Vec<u8>) -> Result<()>,
{
    let mut sink = StringSink::for_input(input, ratio);
    let mut buf = Vec::new();
    for_each_value(name, input, |value| {
        let Some(value) = value else {
            sink.push(None);
            return Ok(());
        };
        buf.clear();
        f(value, &mut buf)?;
        let text = std::str::from_utf8(&buf).map_err(|e| {
            DataFusionError::Internal(format!("{name} produced invalid UTF-8: {e}"))
        })?;
        sink.push(Some(text));
        Ok(())
    })?;
    Ok(sink.finish())
}

fn map_floats<F>(name: &str, input: &ArrayRef, mut f: F) -> Result<ArrayRef>
where
    F: FnMut(&[u8]) -> Result<Option<f64>>,
{
    let mut builder = Float64Builder::with_capacity(input.len());
    for_each_value(name, input, |value| {
        builder.append_option(match value {
            Some(v) => f(v)?,
            None => None,
        });
        Ok(())
    })?;
    Ok(Arc::new(builder.finish()))
}

/// One list of floats per row; `f` appends a row's values, `None` for a
/// null element.
fn map_float_lists<F>(name: &str, input: &ArrayRef, mut f: F) -> Result<ArrayRef>
where
    F: FnMut(&[u8], &mut Vec<Option<f64>>),
{
    let mut builder = ListBuilder::new(Float64Builder::new());
    let mut values = Vec::new();
    for_each_value(name, input, |value| {
        if let Some(v) = value {
            values.clear();
            f(v, &mut values);
            for x in &values {
                builder.values().append_option(*x);
            }
        }
        builder.append(value.is_some());
        Ok(())
    })?;
    Ok(Arc::new(builder.finish()))
}

/// One list of ORF structs per row.
fn orf_lists(
    name: &str,
    input: &ArrayRef,
    finder: &mut OrfFinder,
    first: i64,
    include_protein: bool,
) -> Result<ArrayRef> {
    let mut start = Int64Builder::new();
    let mut end = Int64Builder::new();
    let mut strand = StringBuilder::new();
    let mut frame = Int32Builder::new();
    let mut aa_len = Int64Builder::new();
    let mut has_start = BooleanBuilder::new();
    let mut has_stop = BooleanBuilder::new();
    let mut protein = StringBuilder::new();
    let mut offsets = vec![0i32];
    let mut validity = Vec::with_capacity(input.len());
    let mut total = 0usize;
    for_each_value(name, input, |value| {
        if let Some(seq) = value {
            finder.find(seq);
            for orf in finder.orfs() {
                start.append_value(orf.start as i64 + first);
                end.append_value(orf.end as i64);
                strand.append_value(if orf.reverse { "-" } else { "+" });
                frame.append_value(i32::from(orf.frame));
                aa_len.append_value(orf.aa_len as i64);
                has_start.append_value(orf.has_start);
                has_stop.append_value(orf.has_stop);
                if include_protein {
                    // Translation output is ASCII.
                    protein.append_value(String::from_utf8_lossy(finder.protein(orf)));
                }
            }
            total += finder.orfs().len();
        }
        offsets.push(i32::try_from(total).map_err(|_| {
            DataFusionError::Execution(format!("{name}: too many ORFs in one batch"))
        })?);
        validity.push(value.is_some());
        Ok(())
    })?;
    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(start.finish()),
        Arc::new(end.finish()),
        Arc::new(strand.finish()),
        Arc::new(frame.finish()),
        Arc::new(aa_len.finish()),
        Arc::new(has_start.finish()),
        Arc::new(has_stop.finish()),
    ];
    if include_protein {
        columns.push(Arc::new(protein.finish()));
    }
    let fields = orf_fields(include_protein);
    let item = Arc::new(Field::new_list_field(
        DataType::Struct(fields.clone()),
        true,
    ));
    let orfs = StructArray::try_new(fields, columns, None)?;
    Ok(Arc::new(ListArray::try_new(
        item,
        OffsetBuffer::new(offsets.into()),
        Arc::new(orfs),
        Some(NullBuffer::from(validity)),
    )?))
}

/// One `{sequence, quality_scores}` struct per row, both trimmed to the span
/// the trimmer keeps. A null in either input gives a null struct.
fn trim_structs(
    name: &str,
    sequence: &ArrayRef,
    quality: &ArrayRef,
    trimmer: &QualityTrimmer,
) -> Result<ArrayRef> {
    let sequences = Strings::new(name, sequence)?;
    let qualities = Strings::new(name, quality)?;
    let mut sequence_out = StringSink::for_input(sequence, OutputRatio::Same);
    let mut quality_out = StringSink::for_input(quality, OutputRatio::Same);
    let mut validity = Vec::with_capacity(sequence.len());
    for i in 0..sequence.len() {
        let (Some(seq), Some(qual)) = (sequences.get(i), qualities.get(i)) else {
            sequence_out.push(None);
            quality_out.push(None);
            validity.push(false);
            continue;
        };
        let (start, end) = trimmer
            .span(seq.as_bytes(), qual.as_bytes())
            .map_err(|e| DataFusionError::Execution(format!("{name}: {e}")))?;
        // A span is a byte range. Reads are ASCII; refuse to cut a
        // multi-byte character in two rather than emit invalid UTF-8.
        let (Some(seq), Some(qual)) = (seq.get(start..end), qual.get(start..end)) else {
            return exec_err!("{name}: cannot trim inside a multi-byte character");
        };
        sequence_out.push(Some(seq));
        quality_out.push(Some(qual));
        validity.push(true);
    }
    let fields = trim_fields(
        same_string_type(Some(sequence.data_type())),
        same_string_type(Some(quality.data_type())),
    );
    Ok(Arc::new(StructArray::try_new(
        fields,
        vec![sequence_out.finish(), quality_out.finish()],
        Some(NullBuffer::from(validity)),
    )?))
}

/// `seq_revcomp(seq [, alphabet [, strict]])`
pub fn seq_revcomp_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Revcomp))
}

/// `seq_reverse(seq)`, e.g. to reverse a quality string alongside `seq_revcomp`.
pub fn seq_reverse_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Reverse))
}

/// `seq_complement(seq [, alphabet [, strict]])`
pub fn seq_complement_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Complement))
}

/// `seq_gc_content(seq [, denominator])`
pub fn seq_gc_content_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::GcContent))
}

/// `seq_translate(seq [, frame [, table [, to_stop [, stop_symbol [, cds [, ambiguous]]]]]])`
pub fn seq_translate_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Translate))
}

/// `seq_cds_status(seq [, frame [, table]])`: null for a complete CDS,
/// otherwise the first failing check.
pub fn seq_cds_status_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::CdsStatus))
}

/// `seq_transcribe(seq)`
pub fn seq_transcribe_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Transcribe))
}

/// `seq_back_transcribe(seq)`
pub fn seq_back_transcribe_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::BackTranscribe))
}

/// `seq_mean_quality(qual [, offset [, method]])`
pub fn seq_mean_quality_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::MeanQuality))
}

/// `seq_expected_errors(qual [, offset])`
pub fn seq_expected_errors_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::ExpectedErrors))
}

/// `seq_min_quality(qual [, offset])`
pub fn seq_min_quality_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::MinQuality))
}

/// `seq_max_quality(qual [, offset])`
pub fn seq_max_quality_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::MaxQuality))
}

/// `seq_orfs(seq [, table [, min_len [, strands [, start_codons [, allow_partial [, include_protein [, zero_based]]]]]]])`
pub fn seq_orfs_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Orfs))
}

/// `seq_find_motif(seq, pattern [, both_strands [, overlapping [, zero_based]]])`
pub fn seq_find_motif_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::FindMotif))
}

/// `seq_gc_content_windows(seq, window [, step [, denominator]])`
pub fn seq_gc_content_windows_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::GcContentWindows))
}

/// `seq_gc_skew(seq [, window [, step]])`
pub fn seq_gc_skew_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::GcSkew))
}

/// `seq_molecular_weight(seq [, seq_type [, double_stranded [, circular [, monoisotopic]]]])`
pub fn seq_molecular_weight_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::MolecularWeight))
}

/// `seq_hash(seq [, canonical])`
pub fn seq_hash_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Hash))
}

/// `seq_trim_quality(seq, qual [, leading [, trailing [, window_size [, window_quality [, offset]]]]])`
pub fn seq_trim_quality_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::TrimQuality))
}

/// Every `seq_*` function.
pub fn seq_udfs() -> Vec<ScalarUDF> {
    vec![
        seq_revcomp_udf(),
        seq_reverse_udf(),
        seq_complement_udf(),
        seq_gc_content_udf(),
        seq_translate_udf(),
        seq_cds_status_udf(),
        seq_transcribe_udf(),
        seq_back_transcribe_udf(),
        seq_mean_quality_udf(),
        seq_expected_errors_udf(),
        seq_min_quality_udf(),
        seq_max_quality_udf(),
        seq_orfs_udf(),
        seq_find_motif_udf(),
        seq_gc_content_windows_udf(),
        seq_gc_skew_udf(),
        seq_molecular_weight_udf(),
        seq_hash_udf(),
        seq_trim_quality_udf(),
    ]
}

/// Registers every `seq_*` function on `ctx`.
pub fn register_seq_functions(ctx: &SessionContext) {
    for udf in seq_udfs() {
        ctx.register_udf(udf);
    }
}
