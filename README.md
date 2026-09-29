# fpsearch-rs
Thresholded top-*k* Tanimoto search over binary fingerprints, in Rust.

[![CI](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

The index is sorted by popcount and memory-mapped. Tanimoto cannot exceed
`min(|a|,|b|) / max(|a|,|b|)`, so a threshold becomes a contiguous popcount band (the
BitBound bound) and a query never reads the records outside it.

On 2,854,800 ChEMBL 36 molecules (ECFP4, 2048 bits, Apple M4), a top-10 search at
Tanimoto ≥ 0.95 has a median of **2.6 ms on 10 threads** and 5.5 ms on one. Most of that
comes from the bound, which skips 89% of the database. One thread takes 71 ms to scan the
whole index, and 10 threads add a further 2.1× at 0.95.

```bash
pip install rdkit
curl -O https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/chembl_36_chemreps.txt.gz
python3 tools/chembl_to_fps.py chembl_36_chemreps.txt.gz chembl.fps
cargo run --release -- build chembl.fps chembl.idx
cargo run --release -- query chembl.idx <hex> --threshold 0.9 --top-k 10
cargo test          # 24 tests
```

As a library:

```rust
use fpsearch::Index;

let index = Index::open("chembl.idx")?;
let (hits, stats) = index.search_parallel(&query, 0.9, 10, 10)?;
println!("{} hits, {:.1}% pruned", hits.len(), stats.pruned_fraction() * 100.0);
```

### What is and is not claimed

- The bound is BitBound, and chemfp (Dalke 2019) is its reference implementation. No new
  algorithm is claimed. No comparison against chemfp or FPSim2 has been run, so nothing is
  claimed about how this ranks against them.
- As a sanity check, an exhaustive single-thread scan costs 25 ms per million fingerprints,
  against 177 ms for RDKit's `BulkTanimotoSimilarity` plus the Python sort in
  `tools/baseline_rdkit.py`. The RDKit figure includes that sort and a smaller subset, so
  the 7× gap is a rough ratio, not a kernel-to-kernel measurement.
- Why 10 threads give only 2.0 to 2.7× has not been measured. An earlier claim of sub-second
  search over a billion fingerprints was never measured and is withdrawn.

### More

- [Analysis](ANALYSIS.md), the reasoning, the withdrawn claim and the missing baseline
- [Benchmarks](benchmarks/RESULTS.md), every measured number and how it was produced
- [Design](docs/DESIGN.md), the index layout, the packaging plan and the traps it avoids
