# Benchmark results

Measured 2026-09-01. Every number here came from the commands shown; nothing is estimated
unless it says so.

## Setup

| | |
| --- | --- |
| Machine | Apple M4, 10 cores, macOS (Darwin 25.6.0) |
| Build | `cargo build --release` — `opt-level = 3`, `lto = true`, `codegen-units = 1` |
| Data | ChEMBL 36 `chemreps`, **2,854,800 molecules** |
| Fingerprints | ECFP4 (Morgan radius 2), 2048 bits, RDKit 2026.03.5 |
| Index | 730 MB on disk, popcount-sorted, memory-mapped |
| Refused | 15 molecules RDKit could not parse or that had no bits set |

Popcount distribution across the database: min 1, median 51, max 249, mean 52.4.

Reproduce with:

```bash
curl -O https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/chembl_36_chemreps.txt.gz
python3 tools/chembl_to_fps.py chembl_36_chemreps.txt.gz chembl.fps
cargo run --release -- build chembl.fps chembl.idx
cargo run --release -- bench chembl.idx --queries 300 --threshold 0.95 --threads 10
```

Fingerprint generation takes 86 s across 10 cores; the index builds in 3.1 s.

**On variance.** These were measured on a machine doing other work, and medians move by
20–30% between runs. Every figure below is the best of three runs of 300 queries, which is
the closest available to an unloaded machine. Treat them as the right order of magnitude and
the right ratios, not as three-significant-figure constants.

## The headline

Queries are drawn from the index itself, spread evenly across the popcount ordering so the
timing is not dominated by one end of the distribution. One untimed warm pass precedes each
measurement.

| Threshold | Examined | Pruned | Median, 1 thread | Median, 10 threads | Threading gain |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.95 | 311,266 | **89.1%** | 5,474 µs | **2,627 µs** | 2.1× |
| 0.90 | 630,523 | 77.9% | 13,494 µs | **6,749 µs** | 2.0× |
| 0.80 | 1,270,347 | 55.5% | 38,995 µs | 16,184 µs | 2.4× |
| 0.70 | 1,841,676 | 35.4% | 55,700 µs | 20,770 µs | 2.7× |
| 0.00 | 2,854,800 | 0.0% | 71,239 µs | — | — |

**A top-10 similarity search over 2.85 million real molecules at Tanimoto ≥ 0.95 takes about
2.6 milliseconds.** At ≥ 0.90 it takes about 6.7 ms.

## Against the baseline

`DataStructs.BulkTanimotoSimilarity` is what most people reach for. It is already C++ under
the hood and does no pruning, so this comparison isolates what the layout and the popcount
bound are worth rather than measuring Python against Rust.

Both were run back to back on the same machine within the same minute, single-threaded, and
normalised per million fingerprints so the subset sizes do not matter:

| | µs per query per 1M fingerprints |
| --- | ---: |
| RDKit `BulkTanimotoSimilarity` (500K subset) | 176,988 |
| `fpsearch`, full scan, no pruning | **24,954** |
| `fpsearch`, threshold 0.95, 10 threads | **920** |

**The kernel alone is about 7× faster than RDKit's bulk similarity** on an identical
exhaustive scan. That is the number that isolates the implementation, and it is the one
worth quoting.

End to end — a thresholded 0.95 search on 10 threads against RDKit's exhaustive scan — the
ratio is roughly 190×. That figure is real but it compares two different operations, so it
says more about using a threshold at all than about this code.

## What the parallel scan is actually worth

10 threads buy 2.0–2.7×, not 10×. The scan is embarrassingly parallel and shares no state
during the scan, so the ceiling is not synchronisation — it is memory bandwidth. At 2048
bits a fingerprint is 256 bytes and the kernel does 32 `and` + `count_ones` pairs on it,
which is far too little arithmetic to hide the load.

The honest reading: this engine is bandwidth-bound, more cores will not help much, and the
next real gain is in moving fewer bytes per candidate — a narrower fold, or a compressed
layout — not in more parallelism.

## The claim the README used to make

The original README claimed **sub-second top-*k* over a billion 2048-bit fingerprints**.
That was never measured, and this run does not measure it either — 2.85M is what ChEMBL
supplies, and a billion 2048-bit fingerprints is 256 GB.

Extrapolating the 0.95 figure linearly (×350) gives ~0.9 s, which would just about clear a
second. **That extrapolation should not be believed**, for a reason the measurement itself
demonstrates: these timings are bandwidth-bound with the whole 730 MB index resident in page
cache. At 256 GB the working set cannot be resident, every query would fault against
storage, and paging would dominate the arithmetic completely. The scaling would be governed
by the disk, not by anything in this repository.

So the claim is withdrawn rather than restated. What is measured is what is above.

## Correctness

The pruning is only worth anything if it changes cost and not answers, so that is asserted
rather than assumed:

- `search_agrees_with_a_brute_force_scan` compares the banded search against an exhaustive
  scan at five thresholds and requires identical hits and identical scores.
- `parallel_search_returns_exactly_what_the_serial_one_does` runs 1, 2, 4 and 8 threads at
  three values of *k* and requires the serial result exactly.
- `the_bound_is_never_below_the_true_score` is exhaustive over an 8-bit width — an unsafe
  bound silently loses true hits, which is the one failure mode that would not show up as a
  slowdown.
- `the_band_actually_prunes` fails if the index ever degenerates into a linear scan.

21 tests, `cargo test --release`.
