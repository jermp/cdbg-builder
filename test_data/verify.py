#!/usr/bin/env python3
"""
Correctness check for cdgb-build output.

Given the same filenames list and k used for the build, this script:
  1. Recomputes the expected canonical-k-mer -> {color ids} map from inputs.
  2. Parses <out>.fa, where each FASTA record's header is a color_set_id and
     the body is a colored unitig.
  3. Verifies that
       - every canonical input k-mer appears in exactly one unitig,
       - every k-mer inside a given unitig maps to the same input color set,
       - unitigs sharing a color_set_id all carry the same color set,
       - the number of distinct color_set_ids equals the number of distinct
         color sets in the inputs,
       - <out>.colors exists and is non-empty.

Exits 0 on success, non-zero with a message on the first failure.
"""

import os
import sys
import argparse


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def read_fasta(path: str):
    """Yield (header, seq) records. Lowercases nothing; concatenates body lines."""
    header = None
    body = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(body)
                header = line[1:]
                body = []
            else:
                body.append(line)
        if header is not None:
            yield header, "".join(body)


def kmers(seq: str, k: int):
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        sub = seq[i:i + k]
        if all(c in "ACGT" for c in sub):
            yield sub


def expected_kmer_to_colors(filenames: list, k: int) -> dict:
    """canonical kmer -> frozenset of color ids (0-indexed by file order)."""
    out = {}
    for color, path in enumerate(filenames):
        for _, seq in read_fasta(path):
            for km in kmers(seq, k):
                c = canonical(km)
                out.setdefault(c, set()).add(color)
    return {km: frozenset(s) for km, s in out.items()}


def verify(filenames_list: str, out_basename: str, k: int) -> int:
    with open(filenames_list) as fh:
        files = [ln.strip() for ln in fh if ln.strip()]

    expected = expected_kmer_to_colors(files, k)
    expected_classes = set(expected.values())

    fa_path = out_basename + ".fa"
    colors_path = out_basename + ".colors"

    if not os.path.exists(fa_path):
        print(f"FAIL: {fa_path} not found")
        return 1
    if not os.path.exists(colors_path) or os.path.getsize(colors_path) == 0:
        print(f"FAIL: {colors_path} missing or empty")
        return 1

    seen_kmers = {}            # canonical kmer -> color_set_id from output
    id_to_class = {}           # color_set_id -> frozenset of colors
    num_unitigs = 0

    for header, seq in read_fasta(fa_path):
        try:
            cid = int(header)
        except ValueError:
            print(f"FAIL: non-integer FASTA header '{header}'")
            return 1
        num_unitigs += 1

        unitig_kmers = list(kmers(seq, k))
        if not unitig_kmers:
            print(f"FAIL: unitig with header {cid} has no valid k-mers")
            return 1

        unitig_class = None
        for km in unitig_kmers:
            c = canonical(km)
            if c not in expected:
                print(f"FAIL: unitig k-mer {c} not in any input")
                return 1
            cls = expected[c]
            if unitig_class is None:
                unitig_class = cls
            elif unitig_class != cls:
                print(f"FAIL: unitig {cid} mixes color sets {unitig_class} and {cls}")
                return 1
            if c in seen_kmers:
                print(f"FAIL: k-mer {c} appears in multiple unitigs "
                      f"(ids {seen_kmers[c]} and {cid})")
                return 1
            seen_kmers[c] = cid

        if cid in id_to_class:
            if id_to_class[cid] != unitig_class:
                print(f"FAIL: color_set_id {cid} maps to "
                      f"{id_to_class[cid]} and {unitig_class}")
                return 1
        else:
            id_to_class[cid] = unitig_class

    missing = set(expected) - set(seen_kmers)
    if missing:
        print(f"FAIL: {len(missing)} input k-mer(s) missing from output, "
              f"e.g. {next(iter(missing))}")
        return 1

    distinct_ids = set(id_to_class.keys())
    if len(distinct_ids) != len(expected_classes):
        print(f"FAIL: distinct color_set_ids = {len(distinct_ids)}, "
              f"expected {len(expected_classes)}")
        return 1

    if set(id_to_class.values()) != expected_classes:
        print("FAIL: set of color sets emitted differs from expected")
        print(f"  emitted:  {set(id_to_class.values())}")
        print(f"  expected: {expected_classes}")
        return 1

    print(f"OK: {num_unitigs} unitig(s), {len(seen_kmers)} canonical k-mer(s), "
          f"{len(distinct_ids)} color set(s); {colors_path} present "
          f"({os.path.getsize(colors_path)} bytes)")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filenames", required=True, help="filenames list passed to cdgb-build (-i)")
    ap.add_argument("--out", required=True, help="output basename passed to cdgb-build (-o)")
    ap.add_argument("-k", type=int, required=True)
    args = ap.parse_args()
    sys.exit(verify(args.filenames, args.out, args.k))


if __name__ == "__main__":
    main()
