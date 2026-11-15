//! DataFusion scalar functions over the [`crate::kernels`].
//!
//! | Function | Arguments after the sequence | Returns |
//! |---|---|---|
//! | `seq_revcomp` | `alphabet = 'dna'`, `strict = false` | same string type |
//! | `seq_reverse` | — | same string type |
//! | `seq_complement` | `alphabet = 'dna'`, `strict = false` | same string type |
//! | `seq_gc_content` | `denominator = 'length'` | `Float64` |
//! | `seq_translate` | `frame = 1`, `table = 1`, `to_stop = false` | same string type |
//! | `seq_mean_quality` | `offset = 33`, `method = 'error_prob'` | `Float64` |
//!
//! The first argument may be `Utf8`, `LargeUtf8` or `Utf8View`. Later
//! arguments are optional, positional and must be constants. A null input
//! gives null, and so does a ratio with a zero denominator. The names carry a
//! `seq_` prefix because DataFusion already has a built-in `translate`.
//!
//! Lengths count bytes, so the `length` denominator of `seq_gc_content`
//! counts a non-ASCII character as several bytes. Sequence data is ASCII in
//! practice.

use std::fmt;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, AsArray, Float64Builder, GenericStringBuilder, StringViewBuilder,
};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::{DataFusionError, Result, ScalarValue, exec_err, plan_err};
use datafusion::logical_expr::{
    ColumnarValue, ScalarFunctionArgs, ScalarUDF, ScalarUDFImpl, Signature, Volatility,
};
use datafusion::prelude::SessionContext;

use crate::kernels::{
    Alphabet, GcDenominator, KernelError, QualityMethod, Translator, complement_into, gc_content,
    phred_mean, revcomp_into, reverse_into,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Param {
    Text,
    Int,
    Bool,
}

impl fmt::Display for Param {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Param::Text => "a string",
            Param::Int => "an integer",
            Param::Bool => "a boolean",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Kind {
    Revcomp,
    Reverse,
    Complement,
    GcContent,
    Translate,
    MeanQuality,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Revcomp => "seq_revcomp",
            Kind::Reverse => "seq_reverse",
            Kind::Complement => "seq_complement",
            Kind::GcContent => "seq_gc_content",
            Kind::Translate => "seq_translate",
            Kind::MeanQuality => "seq_mean_quality",
        }
    }

    /// Optional arguments after the sequence, in order.
    fn params(self) -> &'static [(&'static str, Param)] {
        match self {
            Kind::Revcomp | Kind::Complement => {
                &[("alphabet", Param::Text), ("strict", Param::Bool)]
            }
            Kind::Reverse => &[],
            Kind::GcContent => &[("denominator", Param::Text)],
            Kind::Translate => &[
                ("frame", Param::Int),
                ("table", Param::Int),
                ("to_stop", Param::Bool),
            ],
            Kind::MeanQuality => &[("offset", Param::Int), ("method", Param::Text)],
        }
    }

    fn returns_string(self) -> bool {
        matches!(
            self,
            Kind::Revcomp | Kind::Reverse | Kind::Complement | Kind::Translate
        )
    }
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
        let params = self.kind.params();
        if arg_types.is_empty() || arg_types.len() > 1 + params.len() {
            return plan_err!(
                "{name} takes 1 to {} arguments, got {}",
                1 + params.len(),
                arg_types.len()
            );
        }
        let mut coerced = Vec::with_capacity(arg_types.len());
        coerced.push(match &arg_types[0] {
            t @ (DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View) => t.clone(),
            DataType::Null => DataType::Utf8,
            other => {
                return plan_err!(
                    "{name}: the first argument must be a string (Utf8, LargeUtf8 or Utf8View), got {other}"
                );
            }
        });
        for (t, (param, kind)) in arg_types[1..].iter().zip(params) {
            coerced.push(match (kind, t) {
                (
                    Param::Text,
                    DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View | DataType::Null,
                ) => DataType::Utf8,
                (Param::Int, t) if t.is_integer() || t.is_null() => DataType::Int64,
                (Param::Bool, DataType::Boolean | DataType::Null) => DataType::Boolean,
                _ => return plan_err!("{name}: {param} must be {kind}, got {t}"),
            });
        }
        Ok(coerced)
    }

    fn return_type(&self, arg_types: &[DataType]) -> Result<DataType> {
        if !self.kind.returns_string() {
            return Ok(DataType::Float64);
        }
        Ok(match arg_types.first() {
            Some(t @ (DataType::LargeUtf8 | DataType::Utf8View)) => t.clone(),
            _ => DataType::Utf8,
        })
    }

    fn invoke_with_args(&self, args: ScalarFunctionArgs) -> Result<ColumnarValue> {
        let name = self.kind.name();
        let input = args.args[0].to_array(args.number_rows)?;
        let opts = Options::new(self.kind, &args.args[1..])?;
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
            Kind::Translate => {
                let mut translator =
                    Translator::new(opts.int(0, 1)?, opts.int(1, 1)?, opts.bool(2, false)?)
                        .map_err(ke)?;
                map_strings(name, &input, OutputRatio::OneThird, |s, out| {
                    translator.translate_into(s, out);
                    Ok(())
                })?
            }
            Kind::GcContent => {
                let denominator = GcDenominator::parse(opts.text(0, "length")?).map_err(ke)?;
                map_floats(name, &input, |s| Ok(gc_content(s, denominator)))?
            }
            Kind::MeanQuality => {
                let offset = opts.int(0, 33)?;
                let Ok(offset) = u8::try_from(offset) else {
                    return exec_err!("{name}: offset must be between 0 and 255, got {offset}");
                };
                let method = QualityMethod::parse(opts.text(1, "error_prob")?).map_err(ke)?;
                map_floats(name, &input, |q| phred_mean(q, offset, method).map_err(ke))?
            }
        };
        Ok(ColumnarValue::Array(output))
    }
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
        for (arg, (param, _)) in args.iter().zip(params) {
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
            self.kind.params()[i].0
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

    fn bool(&self, i: usize, default: bool) -> Result<bool> {
        match self.args.get(i) {
            None => Ok(default),
            Some(ScalarValue::Boolean(Some(v))) => Ok(*v),
            Some(_) => self.null_err(i),
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

/// `seq_translate(seq [, frame [, table [, to_stop]]])`
pub fn seq_translate_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::Translate))
}

/// `seq_mean_quality(qual [, offset [, method]])`
pub fn seq_mean_quality_udf() -> ScalarUDF {
    ScalarUDF::from(SeqUdf::new(Kind::MeanQuality))
}

/// Every `seq_*` function.
pub fn seq_udfs() -> Vec<ScalarUDF> {
    vec![
        seq_revcomp_udf(),
        seq_reverse_udf(),
        seq_complement_udf(),
        seq_gc_content_udf(),
        seq_translate_udf(),
        seq_mean_quality_udf(),
    ]
}

/// Registers every `seq_*` function on `ctx`.
pub fn register_seq_functions(ctx: &SessionContext) {
    for udf in seq_udfs() {
        ctx.register_udf(udf);
    }
}
