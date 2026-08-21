# fpsearch-rs — Tanimoto search over a billion fingerprints, on one machine

[![CI](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/fpsearch-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** The kernel and its tests are real; the index, the benchmarks and
> the published packages are not built yet. No timing below is measured.

Make-on-demand catalogues are past **48 billion** purchasable compounds, so "find me
molecules similar to this one" stopped being a loop over a list. This is a similarity
search engine for binary fingerprints: memory-mapped storage, a popcount bound that skips
most of the database without comparing it, and SIMD-width word operations in the inner
loop.

**The claim to earn:** sub-second top-*k* Tanimoto over a billion 2048-bit fingerprints on
commodity hardware, benchmarked against `chemfp` and `FPSim2` — **and reported honestly if
it loses.** A repo that says "within 1.4× of the state of the art, and here is the
profiling that explains the gap" is worth more than an unverifiable speed claim.

## How it goes fast

Nothing exotic; the point is that each step is measured rather than assumed.

1. **Bit-count bound.** For fingerprints `a` and `b`, Tanimoto cannot exceed
   `min(|a|, |b|) / max(|a|, |b|)`. Sorting the database by popcount means a threshold
   query touches a contiguous band and skips the rest without a single comparison.
2. **Word-parallel popcount.** Fingerprints are stored as `u64` words and intersected with
   `count_ones`, which lowers to a single hardware instruction — 64 bits per step instead
   of one.
3. **Memory-mapped storage.** A billion 2048-bit fingerprints is 256 GB; it never fits in
   RAM, so the OS pages it and the access pattern is made sequential to keep it happy.
4. **No allocation in the hot loop.** The top-*k* heap is preallocated and reused.

## Layout

```
src/lib.rs        the kernel: bitset intersection, Tanimoto, the popcount bound
src/index.rs      memory-mapped, popcount-sorted fingerprint store
src/main.rs       CLI: build an index, query it
benches/          criterion benchmarks, plus the roofline analysis in the README
python/fpsearch/  PyO3 bindings, published to PyPI via maturin
```

## Packaging (this repo is also the release-engineering exercise)

The engine is published, not just posted:

- **crates.io** — the Rust library, semver, documented.
- **PyPI** — `pip install fpsearch`, wheels built for macOS (arm64, x86_64) and Linux via
  maturin in CI, so nobody needs a Rust toolchain to use it.
- **GitHub Releases** — tagged binaries.
- **Zenodo** — a DOI, so it is citable from a thesis.

## Honest positioning

Fast Tanimoto search is not a new idea — `chemfp` has done it well for over a decade, and
`FPSim2` is the standard Python option. This repo does not claim a new algorithm. It claims
a correct, well-tested, well-packaged implementation with published numbers, and it exists
to demonstrate systems engineering on a problem that genuinely needs it.

## Traps this implementation is built to avoid

- **Tanimoto of two empty fingerprints is undefined, not 1.0.** Returning 1.0 makes every
  malformed record a perfect match for every query. It is an explicit error here.
- **Folding collisions are silent.** Folding a sparse fingerprint into 2048 bits makes
  distinct substructures share a bit, which inflates similarity. The fold width is stored
  in the index header and a query built at a different width is rejected rather than
  silently compared.
- **The popcount bound is only valid for the standard Tanimoto.** Applied to a weighted or
  count-based variant it prunes true hits. The index records which metric it was built for.
