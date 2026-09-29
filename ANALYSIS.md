# Analysis

What was built, why it was built that way, what it is not compared against, and why one
headline claim was withdrawn.

## What existed and what was added

The kernel — Tanimoto, the popcount bound, and their tests — was already correct. Missing
was everything that turns a kernel into a search engine: the store, the search, the CLI, and
any measurement at all.

Added: a memory-mapped popcount-sorted index with a versioned header, serial and parallel
top-*k* search, an `id<TAB>hex` interchange format, a four-subcommand CLI, and a benchmark
against real data rather than synthetic fingerprints.

## Design decisions, and the reasoning

**Popcount-sorted storage.** Tanimoto cannot exceed `min(|a|,|b|) / max(|a|,|b|)`, which
depends only on the two popcounts. Sorting by popcount makes the qualifying band a
*contiguous* slice that binary search finds in `O(log n)`, so a thresholded query never
touches the rest. At threshold 0.95 that skips 89% of the database.

**The header records fold width, metric and byte order.** Each exists to refuse a specific
silent error. Fingerprints folded to different widths are not comparable — different
substructures share bits — so a mismatched query is rejected rather than scored. The
popcount bound is valid only for standard binary Tanimoto, so a future count-based variant
cannot inherit a bound that would prune true hits. An index is a local cache, so rather than
pay endian conversion per word, a file from the other kind of machine is refused.

**An all-zero fingerprint is refused at build time.** Tanimoto of two empty fingerprints is
undefined, not 1.0; returning 1.0 makes every malformed record a perfect match for every
query. Refusing at build time keeps that error out of the query loop entirely, and the CLI
counts the refusals rather than silently shrinking the index.

**Parallelism with scoped threads, no thread-pool dependency.** The scan shares no state, so
each worker keeps its own heap and results merge at the end.

**The naive implementation is kept as the reference the fast one is checked against.**
`search_agrees_with_a_brute_force_scan` compares the banded search against an exhaustive
scan at five thresholds; `parallel_search_returns_exactly_what_the_serial_one_does` runs four
thread counts. A faster structure that returns different rows is not a faster structure.

## What was measured

ChEMBL 36, **2,854,800 molecules**, ECFP4 at 2048 bits, 730 MB index.

- Top-10 at Tanimoto ≥ 0.95 takes **about 2.6 ms** on 10 threads and 5.5 ms on one, with
  89.1% of the database skipped.
- Most of that is the bound. One thread scans the whole index in 71 ms and the 0.95 band in
  5.5 ms, a factor of 13.
- 10 threads buy 2.0 to 2.7×, not 10×.
- An exhaustive single-thread scan is about 7× faster than RDKit's `BulkTanimotoSimilarity`
  plus a Python sort of the scores, normalised per million fingerprints. The RDKit side
  includes the sort and used a 500K subset, so this is a sanity check, not a kernel-to-kernel
  measurement.

## The baseline this is not measured against

Fast Tanimoto search is not a new idea. The popcount bound used here is the BitBound
algorithm, and chemfp (Dalke, *Journal of Cheminformatics* 2019,
[doi:10.1186/s13321-019-0398-8](https://doi.org/10.1186/s13321-019-0398-8)) is its reference
implementation, published in part to serve as a baseline for new similarity search
implementations. This repository claims no new algorithm.

RDKit's `BulkTanimotoSimilarity` is a convenience function, not a search engine. The
comparison that would rank this implementation is against chemfp or FPSim2 on the same
fingerprints. chemfp reports a k=1000 search over 1.8M 2048-bit ChEMBL Morgan fingerprints at
27 ms/query with AVX2 popcount. That comparison has not been run here.

Dalke's paper also argues that uncompressed search on modern hardware is limited by memory
bandwidth rather than by popcount, noting that AVX2 search gains about 10% from prefetching.
That is a finding about chemfp on his hardware. It has not been measured for this code.

## Why the billion-fingerprint claim was withdrawn

The original README claimed sub-second top-*k* over a billion 2048-bit fingerprints. It was
never measured, and this run does not measure it either — ChEMBL supplies 2.85M.

Extrapolating linearly gives roughly 0.9 s, which would just clear a second. That
extrapolation should not be believed. Every timing here was taken with the whole 730 MB index
resident in page cache. A billion 2048-bit fingerprints is 256 GB. It could not be resident,
every query would fault against storage, and the scaling would be governed by the disk rather
than by anything in this repository.

Withdrawing beats restating with a caveat. A claim that survives only under an assumption
the data contradicts is not a weaker claim, it is a wrong one.

## Where the time actually goes

The measured split is clear at the level of the algorithm. At 0.95 the bound removes 89% of
the records, and that accounts for a factor of 13 on one thread. Threads add 2.0 to 2.7×.

Below that level it is not known. At 2048 bits a fingerprint is 256 bytes and the kernel does
32 `and` and `count_ones` pairs on it, which is little arithmetic per byte loaded, so memory
bandwidth is a plausible limit. So is the machine. The M4 mixes performance and efficiency
cores, the band is split into equal chunks, and threads are spawned per query. The 10-thread
run at 0.95 reads about 30 GB/s of fingerprint data. A thread sweep that reports GB/s is the
next measurement, before choosing between a narrower fold or compressed layout and better
scheduling.

`count_ones` lowers to a hardware popcount on aarch64. On x86-64 the default target has no
`popcnt`, so a build there needs `RUSTFLAGS="-C target-cpu=native"` (or `+popcnt`) to get it.
The benchmarks above were run on aarch64.
