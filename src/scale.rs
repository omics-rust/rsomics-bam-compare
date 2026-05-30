//! Read-count scaling and geometry validation.

use rsomics_common::{Result, RsomicsError};
use rsomics_coverage_core::BinnedCoverage;

/// The two BAMs must share an identical reference set (name + length + order):
/// deeptools requires matched geometry so bins line up index-for-index.
pub(crate) fn ensure_same_geometry(a: &BinnedCoverage, b: &BinnedCoverage) -> Result<()> {
    if a.chroms.len() != b.chroms.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "BAMs differ in reference count: {} vs {}",
            a.chroms.len(),
            b.chroms.len()
        )));
    }
    for (x, y) in a.chroms.iter().zip(&b.chroms) {
        if x.name != y.name || x.chrom_len != y.chrom_len {
            return Err(RsomicsError::InvalidInput(format!(
                "BAM reference mismatch: {} (len {}) vs {} (len {})",
                x.name, x.chrom_len, y.name, y.chrom_len
            )));
        }
    }
    Ok(())
}

/// deeptools readCount scaling: `scale[i] = min(m1, m2) / mi`.
///
/// A zero read count yields 0 for the empty side and 1 for the other —
/// deeptools' `min/0` would be `inf`, which poisons every bin with NaN.
pub(crate) fn read_count_scale_factors(m1: u64, m2: u64) -> [f64; 2] {
    let min = m1.min(m2) as f64;
    let s1 = if m1 == 0 { 0.0 } else { min / m1 as f64 };
    let s2 = if m2 == 0 { 0.0 } else { min / m2 as f64 };
    [s1, s2]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_factors_equal_libs() {
        let [s1, s2] = read_count_scale_factors(100, 100);
        assert!((s1 - 1.0).abs() < f64::EPSILON);
        assert!((s2 - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn scale_factors_larger_first() {
        let [s1, s2] = read_count_scale_factors(200, 100);
        assert!((s1 - 0.5).abs() < f64::EPSILON);
        assert!((s2 - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn scale_factors_zero_guard() {
        let [s1, s2] = read_count_scale_factors(0, 100);
        assert_eq!(s1, 0.0);
        assert_eq!(s2, 0.0);
    }
}
