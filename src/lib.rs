//! Per-bin comparison of two BAMs as a bedGraph or bigWig track — deeptools
//! `bamCompare`.
//!
//! Both BAMs are binned by the shared [`rsomics_coverage_core`] primitive (same
//! tiling as `bamCoverage`), each scaled, then combined per bin by an
//! [`Operation`] (deeptools default `log2`). The genome binning lives in Layer A;
//! this crate is the bamCompare-specific layer: the readCount scale factors, the
//! per-bin two-value combination, and the bedGraph/bigWig emit.
//!
//! ## Scaling (deeptools `--scaleFactorsMethod readCount`, the default)
//!
//! For mapped read counts `m1`, `m2` (after the FLAG / MAPQ filter), each BAM's
//! scale factor is `min(m1, m2) / mi` — the larger library is scaled down to the
//! smaller (deeptools `bamCompare.get_scale_factors`). `--normalizeUsing` is
//! `None` by default and only valid with `--scaleFactorsMethod None`, so with
//! defaults the readCount factor is the only scaling applied.
//!
//! ## Combination (deeptools `getRatio` / `compute_ratio`)
//!
//! Per bin, with `v1 = scale[0] * cov1`, `v2 = scale[1] * cov2`:
//!
//! - `log2` — `log2((v1 + pc0) / (v2 + pc1))`
//! - `ratio` — `(v1 + pc0) / (v2 + pc1)`
//! - `reciprocal_ratio` — `r` if `r >= 1` else `-1 / r`, where
//!   `r = (v1 + pc0) / (v2 + pc1)`
//! - `subtract` — `v1 - v2`
//! - `add` — `v1 + v2`
//! - `mean` — `(v1 + v2) / 2`
//! - `first` — `v1`
//! - `second` — `v2`
//!
//! The pseudocount `[pc0, pc1]` (deeptools default `[1, 1]`) is added only for
//! the ratio-family operations (`log2` / `ratio` / `reciprocal_ratio`). bedGraph
//! values use Python's `{:g}` format.

#![allow(clippy::cast_precision_loss)]

mod format;
mod operation;
mod output;
mod scale;

pub use format::OutputFormat;
pub use operation::Operation;

use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use rsomics_bbi::{ChromInfo, write_bigwig};
use rsomics_common::{Result, RsomicsError};
use rsomics_coverage_core::{BinFilter, compute_coverage};

use output::{collect_chrom_intervals, write_chrom_bedgraph};
use scale::{ensure_same_geometry, read_count_scale_factors};

#[derive(Debug, Clone)]
pub struct CompareOpts {
    /// Bin size in bases (deeptools default: 50).
    pub bin_size: u32,
    /// Skip reads whose FLAG has any of these bits set (deeptools default 0).
    pub skip_flags: u16,
    /// Minimum mapping quality (deeptools default 0 = no filter).
    pub min_mapq: u8,
    pub operation: Operation,
    /// Pseudocount `[numerator, denominator]` for ratio-family operations
    /// (deeptools default `[1, 1]`).
    pub pseudocount: [f64; 2],
}

impl Default for CompareOpts {
    fn default() -> Self {
        Self {
            bin_size: 50,
            skip_flags: 0,
            min_mapq: 0,
            operation: Operation::Log2,
            pseudocount: [1.0, 1.0],
        }
    }
}

/// Bin both BAMs, combine per bin, emit bedGraph to `output`. Returns line count.
pub fn bam_compare(
    bam1: &Path,
    bam2: &Path,
    output: &mut dyn Write,
    opts: &CompareOpts,
    workers: NonZero<usize>,
) -> Result<u64> {
    let filter = BinFilter {
        skip_flags: opts.skip_flags,
        min_mapq: opts.min_mapq,
    };
    let cov1 = compute_coverage(bam1, opts.bin_size, filter, workers)?;
    let cov2 = compute_coverage(bam2, opts.bin_size, filter, workers)?;

    ensure_same_geometry(&cov1, &cov2)?;

    let scale = read_count_scale_factors(cov1.total_mapped, cov2.total_mapped);

    let bin_size = u64::from(opts.bin_size);
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut lines: u64 = 0;

    for (c1, c2) in cov1.chroms.iter().zip(&cov2.chroms) {
        if c1.bins.is_empty() {
            continue;
        }
        lines += write_chrom_bedgraph(&mut out, c1, c2, bin_size, scale, opts)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(lines)
}

/// Bin both BAMs, combine per bin, write a bigWig file to `output_path`.
pub fn bam_compare_bigwig(
    bam1: &Path,
    bam2: &Path,
    output_path: &Path,
    opts: &CompareOpts,
    workers: NonZero<usize>,
) -> Result<()> {
    let filter = BinFilter {
        skip_flags: opts.skip_flags,
        min_mapq: opts.min_mapq,
    };
    let cov1 = compute_coverage(bam1, opts.bin_size, filter, workers)?;
    let cov2 = compute_coverage(bam2, opts.bin_size, filter, workers)?;

    ensure_same_geometry(&cov1, &cov2)?;

    let scale = read_count_scale_factors(cov1.total_mapped, cov2.total_mapped);
    let bin_size = u64::from(opts.bin_size);

    let mut chroms_info: Vec<ChromInfo> = Vec::new();
    let mut intervals = Vec::new();

    for (chrom_idx, (c1, c2)) in cov1.chroms.iter().zip(&cov2.chroms).enumerate() {
        if c1.bins.is_empty() {
            continue;
        }
        let chrom_id = u32::try_from(chrom_idx)
            .map_err(|_| RsomicsError::InvalidInput("too many chromosomes".into()))?;
        chroms_info.push(ChromInfo {
            name: c1.name.clone(),
            id: chrom_id,
            length: u32::try_from(c1.chrom_len).unwrap_or(u32::MAX),
        });
        collect_chrom_intervals(c1, c2, chrom_id, bin_size, scale, opts, &mut intervals);
    }

    let mut out = std::fs::File::create(output_path).map_err(RsomicsError::Io)?;
    write_bigwig(&mut out, &chroms_info, &intervals, opts.bin_size)?;
    Ok(())
}
