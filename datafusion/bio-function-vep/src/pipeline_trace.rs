use std::fmt::Write;
use std::sync::OnceLock;
use std::time::{Duration, Instant};

static TRACE_START: OnceLock<Instant> = OnceLock::new();

#[derive(Debug, Clone, Copy)]
pub(crate) enum PipelineTraceValue<'a> {
    Str(&'a str),
    Usize(usize),
    Duration(Duration),
}

pub(crate) fn enabled_from_env_value(value: Option<&str>) -> bool {
    value.is_some_and(|value| value != "0" && !value.eq_ignore_ascii_case("false"))
}

pub(crate) fn enabled() -> bool {
    enabled_from_env_value(std::env::var("VEP_PIPELINE_TRACE").ok().as_deref())
}

pub(crate) fn emit(stage: &str, event: &str, fields: &[(&str, PipelineTraceValue<'_>)]) {
    if !enabled() {
        return;
    }
    eprintln!(
        "{}",
        format_line_at_ms(trace_elapsed_ms(), stage, event, fields)
    );
}

pub(crate) fn format_line_at_ms(
    t_ms: f64,
    stage: &str,
    event: &str,
    fields: &[(&str, PipelineTraceValue<'_>)],
) -> String {
    let mut line = format!("[VEP_PIPELINE_TRACE] t_ms={t_ms:.3} stage={stage} event={event}");
    for (key, value) in fields {
        match value {
            PipelineTraceValue::Str(value) => {
                let _ = write!(line, " {key}={value}");
            }
            PipelineTraceValue::Usize(value) => {
                let _ = write!(line, " {key}={value}");
            }
            PipelineTraceValue::Duration(value) => {
                let _ = write!(line, " {key}_ms={:.3}", value.as_secs_f64() * 1000.0);
            }
        }
    }
    line
}

fn trace_elapsed_ms() -> f64 {
    TRACE_START
        .get_or_init(Instant::now)
        .elapsed()
        .as_secs_f64()
        * 1000.0
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use super::{PipelineTraceValue, enabled_from_env_value, format_line_at_ms};

    #[test]
    fn test_pipeline_trace_env_is_disabled_by_default_and_zero() {
        assert!(!enabled_from_env_value(None));
        assert!(!enabled_from_env_value(Some("0")));
        assert!(enabled_from_env_value(Some("1")));
        assert!(enabled_from_env_value(Some("true")));
    }

    #[test]
    fn test_pipeline_trace_line_formats_stable_fields() {
        let line = format_line_at_ms(
            12.34567,
            "annotation",
            "done",
            &[
                ("chrom", PipelineTraceValue::Str("chr1")),
                ("buffer_id", PipelineTraceValue::Usize(7)),
                ("rows", PipelineTraceValue::Usize(5000)),
                (
                    "elapsed",
                    PipelineTraceValue::Duration(Duration::from_micros(123_456)),
                ),
            ],
        );

        assert_eq!(
            line,
            "[VEP_PIPELINE_TRACE] t_ms=12.346 stage=annotation event=done chrom=chr1 buffer_id=7 rows=5000 elapsed_ms=123.456"
        );
    }
}
