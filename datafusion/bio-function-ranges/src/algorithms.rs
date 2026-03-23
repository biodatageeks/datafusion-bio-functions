use datafusion::config::ConfigField;
use std::str::FromStr;

#[derive(Debug, Eq, PartialEq, Default, Clone, Copy)]
pub enum Algorithm {
    #[default]
    Coitrees,
    IntervalTree,
    ArrayIntervalTree,
    Lapper,
    SuperIntervals,
    CoitreesNearest,
    CoitreesCountOverlaps,
}

#[derive(Debug)]
pub struct ParseAlgorithmError(String);

impl std::fmt::Display for ParseAlgorithmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ParseAlgorithmError {}

impl FromStr for Algorithm {
    type Err = ParseAlgorithmError;

    #[inline]
    fn from_str(s: &str) -> Result<Algorithm, Self::Err> {
        match s.to_lowercase().as_str() {
            "coitrees" => Ok(Algorithm::Coitrees),
            "intervaltree" => Ok(Algorithm::IntervalTree),
            "arrayintervaltree" => Ok(Algorithm::ArrayIntervalTree),
            "lapper" => Ok(Algorithm::Lapper),
            "superintervals" => Ok(Algorithm::SuperIntervals),
            "coitreesnearest" => Ok(Algorithm::CoitreesNearest),
            "coitreescountoverlaps" => Ok(Algorithm::CoitreesCountOverlaps),
            _ => Err(ParseAlgorithmError(format!(
                "Can't parse '{s}' as Algorithm"
            ))),
        }
    }
}

impl std::fmt::Display for Algorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let val = match self {
            Algorithm::Coitrees => "Coitrees",
            Algorithm::IntervalTree => "IntervalTree",
            Algorithm::ArrayIntervalTree => "ArrayIntervalTree",
            Algorithm::Lapper => "Lapper",
            Algorithm::SuperIntervals => "SuperIntervals",
            Algorithm::CoitreesNearest => "CoitreesNearest",
            Algorithm::CoitreesCountOverlaps => "CoitreesCountOverlaps",
        };
        write!(f, "{val}")
    }
}

impl From<ParseAlgorithmError> for datafusion::error::DataFusionError {
    fn from(e: ParseAlgorithmError) -> Self {
        datafusion::error::DataFusionError::External(Box::new(e))
    }
}

impl ConfigField for Algorithm {
    fn set(&mut self, _key: &str, value: &str) -> datafusion::common::Result<()> {
        *self = value.parse::<Algorithm>()?;
        Ok(())
    }

    fn visit<V: datafusion::config::Visit>(&self, visitor: &mut V, name: &str, doc: &'static str) {
        visitor.some(name, self, doc)
    }
}
