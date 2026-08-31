# fpsearch-rs — design notes

Measured performance lives in [benchmarks/RESULTS.md](../benchmarks/RESULTS.md). This file
is the reasoning behind the implementation.

## The index format

```text
offset 0    header      64 bytes
            ids         count * u64
            popcounts   count * u32   (padded to an 8-byte boundary)
            payload     count * n_words * u64
```

Records are sorted by popcount ascending, which is the whole design: the candidate band from
`candidate_band` becomes a contiguous slice that binary search locates in `O(log n)`, and
everything outside it is provably below the threshold and never paged in.

Every section starts on an 8-byte boundary and an mmap base is page-aligned, so the `u64`
views into the payload are always correctly aligned. `Index::open` asserts it rather than
assuming it.

The header records the writer's byte order and the metric the index was built for. An index
is a local cache, not an interchange format — so rather than pay for endian conversion on
every word, a file written on the other kind of machine is refused.

## Packaging

The engine is meant to be published, not just posted. None of this is done yet:

- **crates.io** — the Rust library, semver, documented.
- **PyPI** — `pip install fpsearch`, wheels for macOS (arm64, x86_64) and Linux via maturin
  in CI, so nobody needs a Rust toolchain.
- **GitHub Releases** — tagged binaries.
- **Zenodo** — a DOI, so it is citable from a thesis.

The `.fps` text format exists partly to defer this: anything that can emit `id<TAB>hex` can
feed the index builder today, without bindings.

## Traps this implementation is built to avoid

Each of these is covered by a test, because each produces plausible output rather than an
error.

- **Tanimoto of two empty fingerprints is undefined, not 1.0.** Returning 1.0 makes every
  malformed record a perfect match for every query. It is an explicit error, and
  `IndexBuilder::push` refuses an all-zero fingerprint at build time so the error cannot
  surface later inside a query loop. The CLI counts refusals and prints them — a silent skip
  would leave the index quietly smaller than its input.
- **Folding collisions are silent.** Folding a sparse fingerprint into 2048 bits makes
  distinct substructures share a bit, which inflates similarity. The fold width is in the
  index header and a query built at a different width is rejected, not compared.
- **The popcount bound is only valid for the standard Tanimoto.** Applied to a weighted or
  count-based variant it prunes true hits. The index records which metric it was built for,
  so a future variant cannot silently inherit a bound that is wrong for it.
- **A pruned search that returns different answers is not a faster search.** The band is only
  worth having if it changes cost and nothing else, so `search_agrees_with_a_brute_force_scan`
  compares against an exhaustive scan at five thresholds, and
  `parallel_search_returns_exactly_what_the_serial_one_does` requires the threaded path to
  match the serial one exactly.
- **An unsafe bound loses hits without ever looking slow.** `the_bound_is_never_below_the_true_score`
  is exhaustive over an 8-bit width rather than sampled.
- **A truncated index should not be read as a short one.** The header carries the record
  count; `open` computes the size that implies and refuses a file shorter than that, instead
  of mapping whatever is there and returning fewer hits.

## Where the time actually goes

The scan is bandwidth-bound, not compute-bound: at 2048 bits a fingerprint is 256 bytes and
the kernel does 32 `and` + `count_ones` pairs on it, far too little arithmetic to hide the
load. That is why 10 threads buy 2.4–4.3× rather than 10×, and it is the thing to attack
next — a narrower fold, or a compressed layout that moves fewer bytes per candidate, would
help where more cores will not.
