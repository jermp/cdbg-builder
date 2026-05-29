#!/usr/bin/env python3
"""
Correctness check for cdgb-build output.

Given the same filenames list and k used for the build, this script:
  1. Recomputes the expected canonical-k-mer -> color-bitmask map from inputs.
  2. Parses <out>.fa, where each FASTA record's header is a color_set_id and
     the body is a colored unitig.
  3. Verifies that
       - every canonical input k-mer appears in exactly one unitig,
       - every k-mer inside a given unitig maps to the same input color set,
       - unitigs sharing a color_set_id all carry the same color set,
       - the number of distinct color_set_ids equals the number of distinct
         color sets in the inputs,
       - <out>.color_sets exists and is non-empty.

K-mers are 2-bit packed into Python ints (A=0, C=1, G=2, T=3) and reduced to
their canonical form (min of forward / reverse complement) so the verifier
scales to tens of millions of k-mers without holding 31-byte strings.

Exits 0 on success, non-zero with a message on the first failure.
"""

import argparse
import gzip
import os
import sys


# A=0, C=1, G=2, T=3, everything else = 0xff (= "split the window")
_ACGT_TABLE = bytes(
    {0x41: 0, 0x43: 1, 0x47: 2, 0x54: 3,
     0x61: 0, 0x63: 1, 0x67: 2, 0x74: 3}.get(b, 0xff)
    for b in range(256)
)


def _open(path: str):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def stream_records(path: str):
    """Yield raw sequence bytes per record; concatenates wrapped lines."""
    body = []
    have_header = False
    with _open(path) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(b">"):
                if have_header:
                    yield b"".join(body)
                body = []
                have_header = True
            else:
                body.append(line.rstrip(b"\r\n"))
        if have_header:
            yield b"".join(body)


def canonical_kmers(seq: bytes, k: int):
    """Yield 2-bit packed canonical k-mer ints; restarts at any non-ACGT base."""
    mask = (1 << (2 * k)) - 1
    rc_high_shift = 2 * (k - 1)
    fwd = 0
    rc = 0
    filled = 0
    table = _ACGT_TABLE
    for byte in seq:
        v = table[byte]
        if v == 0xff:
            filled = 0
            fwd = 0
            rc = 0
            continue
        fwd = ((fwd << 2) | v) & mask
        rc = (rc >> 2) | ((v ^ 3) << rc_high_shift)
        if filled + 1 < k:
            filled += 1
            continue
        yield fwd if fwd <= rc else rc


def expected_colors_map(filenames: list, k: int) -> dict:
    """Returns {canonical_kmer_int: color_bitmask_int}."""
    out = {}
    for color, path in enumerate(filenames):
        bit = 1 << color
        for seq in stream_records(path):
            for km in canonical_kmers(seq, k):
                out[km] = out.get(km, 0) | bit
    return out


def verify(filenames_list: str, out_basename: str, k: int) -> int:
    with open(filenames_list) as fh:
        files = [ln.strip() for ln in fh if ln.strip()]

    print(f"computing expected k-mer set from {len(files)} input(s)...", flush=True)
    expected = expected_colors_map(files, k)
    expected_classes = set(expected.values())
    print(f"  {len(expected)} canonical k-mers, "
          f"{len(expected_classes)} distinct color set(s)", flush=True)

    fa_path = out_basename + ".fa"
    colors_path = out_basename + ".color_sets"

    if not os.path.exists(fa_path):
        print(f"FAIL: {fa_path} not found")
        return 1
    if not os.path.exists(colors_path) or os.path.getsize(colors_path) == 0:
        print(f"FAIL: {colors_path} missing or empty")
        return 1

    print(f"checking {fa_path}...", flush=True)
    seen = {}             # canonical kmer -> color_set_id from output
    id_to_class = {}      # color_set_id -> color bitmask
    num_unitigs = 0

    cid = None
    body = []

    def flush(cid, body_bytes):
        nonlocal num_unitigs
        if cid is None:
            return 0
        num_unitigs += 1
        unitig_class = None
        for km in canonical_kmers(body_bytes, k):
            cls = expected.get(km)
            if cls is None:
                print(f"FAIL: unitig k-mer (cid={cid}) not in any input")
                return 1
            if unitig_class is None:
                unitig_class = cls
            elif unitig_class != cls:
                print(f"FAIL: unitig {cid} mixes color sets "
                      f"{unitig_class:b} and {cls:b}")
                return 1
            if km in seen:
                print(f"FAIL: k-mer appears in multiple unitigs "
                      f"(ids {seen[km]} and {cid})")
                return 1
            seen[km] = cid
        if unitig_class is None:
            print(f"FAIL: unitig with header {cid} has no valid k-mers")
            return 1
        prev = id_to_class.get(cid)
        if prev is None:
            id_to_class[cid] = unitig_class
        elif prev != unitig_class:
            print(f"FAIL: color_set_id {cid} maps to "
                  f"{prev:b} and {unitig_class:b}")
            return 1
        return 0

    with open(fa_path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                rc = flush(cid, b"".join(body))
                if rc:
                    return rc
                try:
                    cid = int(line[1:].rstrip())
                except ValueError:
                    print(f"FAIL: non-integer FASTA header '{line[1:].rstrip()!r}'")
                    return 1
                body = []
            else:
                body.append(line.rstrip(b"\r\n"))
        rc = flush(cid, b"".join(body))
        if rc:
            return rc

    if len(seen) != len(expected):
        missing = len(expected) - len(seen)
        print(f"FAIL: {missing} input k-mer(s) missing from output")
        return 1

    distinct_ids = set(id_to_class.keys())
    if len(distinct_ids) != len(expected_classes):
        print(f"FAIL: distinct color_set_ids = {len(distinct_ids)}, "
              f"expected {len(expected_classes)}")
        return 1

    if set(id_to_class.values()) != expected_classes:
        print("FAIL: set of color sets emitted differs from expected")
        return 1

    print(f"OK: {num_unitigs} unitig(s), {len(seen)} canonical k-mer(s), "
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
