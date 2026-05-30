//! Output format selection.

/// Output format for the comparison track.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    BedGraph,
    #[default]
    BigWig,
}

impl std::str::FromStr for OutputFormat {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "bedgraph" | "bedGraph" => Ok(Self::BedGraph),
            "bigwig" | "bigWig" | "BigWig" => Ok(Self::BigWig),
            _ => Err(format!(
                "unknown output format '{s}'; choose bedgraph or bigwig"
            )),
        }
    }
}
