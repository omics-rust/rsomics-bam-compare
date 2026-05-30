//! Per-bin combination operations (deeptools `getRatio` / `compute_ratio`).

use crate::CompareOpts;

/// How two scaled per-bin coverage values are combined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    Log2,
    Ratio,
    ReciprocalRatio,
    Subtract,
    Add,
    Mean,
    First,
    Second,
}

impl std::str::FromStr for Operation {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "log2" => Ok(Self::Log2),
            "ratio" => Ok(Self::Ratio),
            "reciprocal_ratio" => Ok(Self::ReciprocalRatio),
            "subtract" => Ok(Self::Subtract),
            "add" => Ok(Self::Add),
            "mean" => Ok(Self::Mean),
            "first" => Ok(Self::First),
            "second" => Ok(Self::Second),
            _ => Err(format!(
                "unknown operation '{s}'; choose log2 ratio reciprocal_ratio \
                 subtract add mean first second"
            )),
        }
    }
}

/// Mirrors deeptools `getRatio` + `compute_ratio`: ratio-family ops add the
/// pseudocount and divide (float `inf` / `nan` propagate as numpy would),
/// others operate directly.
pub(crate) fn combine(cov1: u32, cov2: u32, scale: [f64; 2], opts: &CompareOpts) -> f64 {
    let v1 = scale[0] * f64::from(cov1);
    let v2 = scale[1] * f64::from(cov2);

    match opts.operation {
        Operation::Subtract => v1 - v2,
        Operation::Add => v1 + v2,
        // Not `f64::midpoint`: deeptools/numpy compute `(v1 + v2) / 2.0`, and
        // the bedGraph must match it bit-for-bit.
        #[allow(clippy::manual_midpoint)]
        Operation::Mean => (v1 + v2) / 2.0,
        Operation::First => v1,
        Operation::Second => v2,
        Operation::Log2 | Operation::Ratio | Operation::ReciprocalRatio => {
            let num = v1 + opts.pseudocount[0];
            let den = v2 + opts.pseudocount[1];
            let ratio = num / den;
            match opts.operation {
                Operation::Log2 => ratio.log2(),
                Operation::Ratio => ratio,
                Operation::ReciprocalRatio => {
                    if ratio >= 1.0 {
                        ratio
                    } else {
                        -1.0 / ratio
                    }
                }
                _ => unreachable!(),
            }
        }
    }
}
