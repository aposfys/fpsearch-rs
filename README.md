# fpsearch-rs
Tanimoto similarity search over large binary fingerprint collections.

[![CI](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** The kernel and its tests are real. The index, the benchmarks and the published packages are not built yet, and no timing has been measured.

Make-on-demand catalogues are past 48 billion purchasable compounds, so "find me molecules similar to this one" stopped being a loop over a list. This is a similarity search engine for binary fingerprints: memory-mapped storage, a popcount bound that skips most of the database without comparing it, and SIMD-width word operations in the inner loop.

**The claim to earn:** sub-second top-*k* Tanimoto over a billion 2048-bit fingerprints on commodity hardware, benchmarked against `chemfp` and `FPSim2` — and reported honestly if it loses.

Fast Tanimoto search is not a new idea; `chemfp` has done it well for over a decade. This repo does not claim a new algorithm. It claims a correct, well-tested, well-packaged implementation with published numbers.

### Layout
```
src/lib.rs        the kernel: bitset intersection, Tanimoto, the popcount bound
benches/          criterion benchmarks
```
Planned: `src/index.rs` (memory-mapped, popcount-sorted store), `src/main.rs` (CLI), `python/fpsearch/` (PyO3 bindings, published via maturin).

### Design notes
[How it goes fast, the packaging plan, and the traps the implementation is built to avoid](docs/DESIGN.md)
