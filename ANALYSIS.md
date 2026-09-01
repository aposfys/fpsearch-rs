# Analysis

What was built, why it was built that way, and why one headline claim was withdrawn.

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

- Top-10 at Tanimoto ≥ 0.95: **~2.6 ms**, 89.1% of the database skipped.
- The kernel alone, on an identical exhaustive scan, is **~7× faster than RDKit's
  `BulkTanimotoSimilarity`** — measured back to back, single-threaded, normalised per
  million fingerprints.
- 10 threads buy 2.0–2.7×, not 10×.

## Why the billion-fingerprint claim was withdrawn

The original README claimed sub-second top-*k* over a billion 2048-bit fingerprints. It was
never measured, and this run does not measure it either — ChEMBL supplies 2.85M.

Extrapolating linearly gives roughly 0.9 s, which would just clear a second. That
extrapolation should not be believed, and the measurements are what show why: the scan is
**bandwidth-bound with the whole 730 MB index resident in page cache**. A billion 2048-bit
fingerprints is 256 GB. It could not be resident, every query would fault against storage,
and the scaling would be governed by the disk rather than by anything in this repository.

Withdrawing beats restating with a caveat. A claim that survives only under an assumption
the data contradicts is not a weaker claim, it is a wrong one.

## Where the time actually goes

At 2048 bits a fingerprint is 256 bytes and the kernel does 32 `and` + `count_ones` pairs on
it — far too little arithmetic to hide the load. That is why more cores help so little, and
it says what to attack next: move fewer bytes per candidate (a narrower fold, or a
compressed layout), not add parallelism.
