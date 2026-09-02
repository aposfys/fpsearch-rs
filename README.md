# fpsearch-rs
Tanimoto similarity search over large binary fingerprint collections.

[![CI](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Memory-mapped storage, a popcount bound that skips most of the database without
comparing it, and word-parallel popcount in the inner loop.

**A top-10 search over 2,854,800 real ChEMBL molecules at Tanimoto ≥ 0.95 takes
about 2.6 ms**, skipping 89% of the database untouched. The kernel alone is
**about 7× faster than RDKit's `BulkTanimotoSimilarity`** on an identical
exhaustive scan.

```
cargo build --release
python3 tools/chembl_to_fps.py chembl_36_chemreps.txt.gz chembl.fps
cargo run --release -- build chembl.fps chembl.idx
cargo run --release -- query chembl.idx <hex> --threshold 0.9 --top-k 10
cargo test          # 21 tests
```

As a library:

```rust
use fpsearch::Index;

let index = Index::open("chembl.idx")?;
let (hits, stats) = index.search_parallel(&query, 0.9, 10, 10)?;
println!("{} hits, {:.1}% pruned", hits.len(), stats.pruned_fraction() * 100.0);
```

### How it goes fast

1. **Bit-count bound.** Tanimoto cannot exceed `min(|a|,|b|) / max(|a|,|b|)`. The
   index is sorted by popcount, so a thresholded query touches one contiguous band
   and skips the rest without a single comparison — 89% of the database at 0.95.
2. **Word-parallel popcount.** Fingerprints are `u64` words intersected with
   `count_ones`, one hardware instruction per 64 bits.
3. **Memory-mapped storage.** The index is mapped, not read, so only the band a
   query can match is ever paged in.
4. **No allocation in the hot loop**, and a scan that splits across threads with
   no shared state.

### The billion-fingerprint claim, withdrawn

An earlier version of this README claimed sub-second top-*k* over a billion
fingerprints. That was never measured. These timings are bandwidth-bound with a
730 MB index resident in page cache; a billion 2048-bit fingerprints is 256 GB,
could not be resident, and would be governed by storage rather than by anything
here. The claim is withdrawn rather than restated.

Fast Tanimoto search is not a new idea — `chemfp` has done it well for over a
decade. This repo claims no new algorithm, only a correct, well-tested
implementation with published numbers.

### More

- [Analysis](ANALYSIS.md) — what was done and why it was done that way
- [Benchmarks](benchmarks/RESULTS.md) — measured results, the baseline comparison, and what the parallel scan is worth
- [Design](docs/DESIGN.md) — packaging plan, the layout, and the traps this avoids
