#!/usr/bin/env python3
"""Time RDKit's own bulk Tanimoto over the same fingerprints, as the baseline to beat.

    python3 tools/baseline_rdkit.py chembl.fps --subset 500000 --queries 20

`DataStructs.BulkTanimotoSimilarity` is what most people reach for, it is already C++ under
the hood, and it does no pruning — so the comparison isolates what the popcount bound and
the memory-mapped layout are actually worth, rather than measuring Python against C.

Reports per-query wall time normalised to a million fingerprints, which is the figure that
can be compared against `fpsearch bench` regardless of subset size.
"""

from __future__ import annotations

import argparse
import random
import statistics
import sys
import time

from rdkit import DataStructs


def load(path: str, limit: int) -> list[DataStructs.ExplicitBitVect]:
    vectors = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            _, hex_fp = line.rstrip("\n").split("\t", 1)
            vectors.append(DataStructs.CreateFromBinaryText(bytes.fromhex(hex_fp)))
            if limit and len(vectors) >= limit:
                break
    return vectors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fps", help="the .fps file")
    parser.add_argument("--subset", type=int, default=500_000, help="fingerprints to load")
    parser.add_argument("--queries", type=int, default=20, help="queries to time")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    started = time.time()
    database = load(args.fps, args.subset)
    print(f"loaded {len(database):,} fingerprints in {time.time() - started:.1f}s", file=sys.stderr)
    if not database:
        sys.exit("no fingerprints loaded")

    random.seed(args.seed)
    queries = [database[random.randrange(len(database))] for _ in range(args.queries)]

    timings = []
    for query in queries:
        start = time.perf_counter()
        scores = DataStructs.BulkTanimotoSimilarity(query, database)
        # Take a top-10 so the comparison includes selection, as fpsearch's does.
        top = sorted(scores, reverse=True)[:10]
        timings.append((time.perf_counter() - start) * 1e6)
        assert top

    median = statistics.median(timings)
    per_million = median / len(database) * 1_000_000
    print(f"database      {len(database):,} fingerprints (no pruning)")
    print(f"queries       {len(timings)}")
    print(f"median        {median:,.0f} µs")
    print(f"normalised    {per_million:,.0f} µs per query per 1M fingerprints")


if __name__ == "__main__":
    main()
