/// Filter operation for interval overlap queries.
///
/// Controls how interval boundaries are treated during overlap detection.
#[derive(Clone, PartialEq, Debug)]
pub enum FilterOp {
    /// Standard overlap: intervals overlap if they share any position.
    Weak = 0,
    /// Strict overlap: boundaries are adjusted inward (start+1, end-1).
    Strict = 1,
}
