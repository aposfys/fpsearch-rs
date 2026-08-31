#!/usr/bin/env python3
"""Turn a ChEMBL chemreps dump into the `.fps` text format the index builder reads.

    python3 tools/chembl_to_fps.py chembl_36_chemreps.txt.gz chembl.fps

Deliberately a separate script rather than a Rust dependency: the engine has no business
linking a chemistry toolkit, and anything that can emit `id<TAB>hex` can feed it.

Molecules RDKit refuses to parse are counted and reported, never silently dropped — a
benchmark run against a quietly smaller database is not the benchmark you think it is.
"""

from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import sys
import time

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

_GENERATOR = None
_N_BITS = 2048
_RADIUS = 2


def _init(n_bits: int, radius: int) -> None:
    global _GENERATOR, _N_BITS, _RADIUS
    _N_BITS, _RADIUS = n_bits, radius
    _GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)


def _encode(row: tuple[str, str]) -> tuple[int, str] | None:
    """One (chembl_id, smiles) pair to (numeric id, hex fingerprint)."""
    chembl_id, smiles = row
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    bitvect = _GENERATOR.GetFingerprint(mol)
    if bitvect.GetNumOnBits() == 0:
        # No bits set: Tanimoto against it is undefined, so it cannot enter an index.
        return None
    # DataStructs.BitVectToBinaryText gives bytes in the same order the Rust reader expects.
    from rdkit import DataStructs

    raw = DataStructs.BitVectToBinaryText(bitvect)
    return int(chembl_id.removeprefix("CHEMBL")), raw.hex()


def rows(path: str):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            id_col = header.index("chembl_id")
            smiles_col = header.index("canonical_smiles")
        except ValueError:
            sys.exit(f"{path}: expected chembl_id and canonical_smiles columns, got {header}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= smiles_col:
                continue
            smiles = parts[smiles_col]
            if smiles:
                yield parts[id_col], smiles


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="chembl_NN_chemreps.txt[.gz]")
    parser.add_argument("output", help="destination .fps file")
    parser.add_argument("--bits", type=int, default=2048, help="fold width (default 2048)")
    parser.add_argument("--radius", type=int, default=2, help="Morgan radius (default 2, ECFP4)")
    parser.add_argument("--limit", type=int, default=0, help="stop after N molecules")
    parser.add_argument("--processes", type=int, default=0, help="worker count (default: all cores)")
    args = parser.parse_args()

    processes = args.processes or mp.cpu_count()
    started = time.time()
    written = failed = 0

    source = rows(args.input)
    if args.limit:
        import itertools

        source = itertools.islice(source, args.limit)

    with mp.Pool(processes, initializer=_init, initargs=(args.bits, args.radius)) as pool, open(
        args.output, "w"
    ) as out:
        out.write(f"# ECFP{args.radius * 2} n_bits={args.bits} source={args.input}\n")
        for result in pool.imap(_encode, source, chunksize=2000):
            if result is None:
                failed += 1
                continue
            identifier, hex_fp = result
            out.write(f"{identifier}\t{hex_fp}\n")
            written += 1
            if written % 250_000 == 0:
                rate = written / (time.time() - started)
                print(f"  {written:,} written ({rate:,.0f}/s)", file=sys.stderr)

    elapsed = time.time() - started
    print(f"wrote {written:,} fingerprints to {args.output} in {elapsed:.1f}s", file=sys.stderr)
    if failed:
        print(f"refused {failed:,} molecules (unparseable, or no bits set)", file=sys.stderr)


if __name__ == "__main__":
    main()
