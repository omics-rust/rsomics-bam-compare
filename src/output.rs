//! bedGraph / bigWig emit with run-length merging and Python-`{:g}` formatting.

use std::io::Write;

use rsomics_bbi::Interval;
use rsomics_common::{Result, RsomicsError};
use rsomics_coverage_core::ChromBins;

use crate::CompareOpts;
use crate::operation::combine;

/// bamCompare tiles each chromosome into 5-Mbp `writeBedGraph_worker` tasks; a
/// bedGraph run never merges across a tile boundary.
const CHUNK: u64 = 5_000_000;

/// Write one chromosome's combined bins as merged bedGraph lines.
///
/// Adjacent equal-value bins are collapsed (deeptools `writeBedGraph`
/// run-length encoding). All bins (including zero-valued) are written.
pub(crate) fn write_chrom_bedgraph(
    out: &mut impl Write,
    c1: &ChromBins,
    c2: &ChromBins,
    bin_size: u64,
    scale: [f64; 2],
    opts: &CompareOpts,
) -> Result<u64> {
    let mut lines: u64 = 0;
    let mut write_start: u64 = 0;
    let mut write_end: u64 = 0;
    let mut prev_val: Option<f64> = None;

    let n = c1.bins.len();
    for i in 0..n {
        let bin_start = i as u64 * bin_size;
        let bin_end = ((i as u64 + 1) * bin_size).min(c1.chrom_len);
        // bamCompare runs writeBedGraph in 5-Mbp tiles (single-process default):
        // a run never merges across a tile boundary. Unlike bamCoverage it keeps
        // tile-trailing 0 runs (a 0 here is a real log2 ratio, not absent data).
        let at_chunk_edge = bin_start.is_multiple_of(CHUNK);

        let value = combine(c1.bins[i], c2.bins[i], scale, opts);

        match prev_val {
            None => {
                write_start = bin_start;
                write_end = bin_end;
                prev_val = Some(value);
            }
            Some(pv) if values_equal(pv, value) && !at_chunk_edge => {
                write_end = bin_end;
            }
            Some(pv) => {
                write_line(out, &c1.name, write_start, write_end, pv)?;
                lines += 1;
                write_start = bin_start;
                write_end = bin_end;
                prev_val = Some(value);
            }
        }

        if i + 1 == n
            && let Some(pv) = prev_val
            && write_start != write_end
        {
            write_line(out, &c1.name, write_start, write_end, pv)?;
            lines += 1;
        }
    }

    Ok(lines)
}

/// Collect run-length-merged intervals for one chromosome into `out`.
#[allow(clippy::cast_possible_truncation)] // genomic coords fit u32; f64→f32 is the bigWig format
pub(crate) fn collect_chrom_intervals(
    c1: &ChromBins,
    c2: &ChromBins,
    chrom_id: u32,
    bin_size: u64,
    scale: [f64; 2],
    opts: &CompareOpts,
    out: &mut Vec<Interval>,
) {
    let n = c1.bins.len();
    let mut write_start: u64 = 0;
    let mut write_end: u64 = 0;
    let mut prev_val: Option<f64> = None;

    let flush = |start: u64, end: u64, val: f64, out: &mut Vec<Interval>| {
        if start == end {
            return;
        }
        out.push(Interval {
            chrom_id,
            start: start as u32,
            end: end as u32,
            value: val as f32,
        });
    };

    for i in 0..n {
        let bin_start = i as u64 * bin_size;
        let bin_end = ((i as u64 + 1) * bin_size).min(c1.chrom_len);
        let value = combine(c1.bins[i], c2.bins[i], scale, opts);

        match prev_val {
            None => {
                write_start = bin_start;
                write_end = bin_end;
                prev_val = Some(value);
            }
            Some(pv) if values_equal(pv, value) => {
                write_end = bin_end;
            }
            Some(pv) => {
                flush(write_start, write_end, pv, out);
                write_start = bin_start;
                write_end = bin_end;
                prev_val = Some(value);
            }
        }

        if i + 1 == n
            && let Some(pv) = prev_val
            && write_start != write_end
        {
            flush(write_start, write_end, pv, out);
        }
    }
}

/// Two values are "the same bin" when they format identically.
///
/// deeptools merges on the written value, so string equality is the canonical
/// test — a relative epsilon would mis-split near the format boundary.
fn values_equal(a: f64, b: f64) -> bool {
    format_g(a) == format_g(b)
}

fn write_line(out: &mut impl Write, chrom: &str, start: u64, end: u64, value: f64) -> Result<()> {
    let s = format_g(value);
    writeln!(out, "{chrom}\t{start}\t{end}\t{s}").map_err(RsomicsError::Io)
}

/// Format a float like Python's `{:g}` (6 significant digits, trailing zeros stripped).
pub(crate) fn format_g(v: f64) -> String {
    if v == 0.0 {
        // {:g} prints negative-zero as "-0"; numpy/deeptools emit "0".
        return "0".to_owned();
    }
    if v.is_nan() {
        return "nan".to_owned();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.to_owned();
    }
    python_g(v)
}

/// Python `{:g}`: 6 significant digits, exponent form outside `1e-4..1e16`,
/// trailing zeros (and bare `.`) stripped.
///
/// Exponent is taken from `{:.5e}` rather than `log10().floor()` — the
/// scientific render already rounds to 6 sig-figs, so its exponent is the
/// post-rounding one (e.g. `999999.6` rounds to `1e6`, which `log10().floor()`
/// mis-buckets as exponent 5).
fn python_g(v: f64) -> String {
    let sci = format!("{v:.5e}");
    let (_, exp_str) = sci.split_once('e').unwrap();
    let exp: i32 = exp_str.parse().unwrap();

    if !(-4..16).contains(&exp) {
        return normalise_exponential(&sci);
    }
    let decimals = usize::try_from((5 - exp).max(0)).unwrap();
    let s = format!("{v:.decimals$}");
    let s = if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        &s
    };
    s.to_owned()
}

/// Rust's `{:e}` → Python `{:g}` exponent form: sign always present, exponent ≥ 2 digits.
fn normalise_exponential(s: &str) -> String {
    let (mantissa, exp) = s.split_once('e').unwrap();
    let mantissa = if mantissa.contains('.') {
        mantissa.trim_end_matches('0').trim_end_matches('.')
    } else {
        mantissa
    };
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(rest) => ('-', rest),
        None => ('+', exp.strip_prefix('+').unwrap_or(exp)),
    };
    format!("{mantissa}e{sign}{digits:0>2}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_g_zero() {
        assert_eq!(format_g(0.0), "0");
        assert_eq!(format_g(-0.0), "0");
    }

    #[test]
    fn format_g_specials() {
        assert_eq!(format_g(f64::NAN), "nan");
        assert_eq!(format_g(f64::INFINITY), "inf");
        assert_eq!(format_g(f64::NEG_INFINITY), "-inf");
    }

    #[test]
    fn format_g_integers() {
        assert_eq!(format_g(1.0), "1");
        assert_eq!(format_g(100.0), "100");
        // Range threshold: code uses fixed-decimal for exp in [-4, 16).
        assert_eq!(format_g(1e16), "1e+16");
    }

    #[test]
    fn format_g_fractions() {
        assert_eq!(format_g(0.5), "0.5");
        assert_eq!(format_g(1.5), "1.5");
        assert_eq!(format_g(1.23456789), "1.23457");
    }

    #[test]
    fn format_g_small() {
        assert_eq!(format_g(1e-5), "1e-05");
        assert_eq!(format_g(1.5e-5), "1.5e-05");
    }
}
