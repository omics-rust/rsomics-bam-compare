use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bam_compare(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bam-compare");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let treat = manifest.join("tests/golden/treat.bam");
    let ctrl = manifest.join("tests/golden/ctrl.bam");
    let out = tempfile::NamedTempFile::new().unwrap();

    c.bench_function("rsomics-bam-compare golden", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([
                    "--bam1",
                    treat.to_str().unwrap(),
                    "--bam2",
                    ctrl.to_str().unwrap(),
                    "--out-file-format",
                    "bedgraph",
                    "-o",
                    out.path().to_str().unwrap(),
                    "--operation",
                    "log2",
                ])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_bam_compare);
criterion_main!(benches);
