#!/usr/bin/env python3
"""
Strict topology check for cdgb-build output.

Beyond the existing verifier (which only confirms k-mer coverage and
monochromaticity), this script checks each unitig is a *maximal* path in
the colored dBG: every internal k-mer has out-degree 1 and in-degree 1 in
the per-color-set subgraph.

Usage:
    strict_topology_check.py --filenames <list> --out <basename> -k <k>
"""

import argparse
import gzip
import os
import sys


_ACGT_TABLE = bytes(
    {0x41: 0, 0x43: 1, 0x47: 2, 0x54: 3,
     0x61: 0, 0x63: 1, 0x67: 2, 0x74: 3}.get(b, 0xff)
    for b in range(256)
)
NUC = "ACGT"


def _open(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def stream_records(path):
    body = []
    have_header = False
    with _open(path) as fh:
        for line in fh:
            if line.startswith(b">"):
                if have_header:
                    yield b"".join(body)
                body = []
                have_header = True
            else:
                body.append(line.rstrip(b"\r\n"))
        if have_header:
            yield b"".join(body)


def revcomp_int(x, k):
    y = 0
    for _ in range(k):
        y = (y << 2) | ((x & 3) ^ 3)
        x >>= 2
    return y


def canonical_int(x, k):
    rc = revcomp_int(x, k)
    return x if x <= rc else rc


def stream_canonical_kmers(seq_bytes, k):
    mask = (1 << (2 * k)) - 1
    rc_high_shift = 2 * (k - 1)
    fwd = 0
    rc = 0
    filled = 0
    table = _ACGT_TABLE
    for byte in seq_bytes:
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


def expected_colors_map(filenames, k):
    out = {}
    for color, path in enumerate(filenames):
        bit = 1 << color
        for seq in stream_records(path):
            for km in stream_canonical_kmers(seq, k):
                out[km] = out.get(km, 0) | bit
    return out


def shift_append(fwd, nt, k):
    mask = (1 << (2 * k)) - 1
    return ((fwd << 2) | (nt & 3)) & mask


def check_unitig_topology(seq, k, kmer_colors, problems):
    """
    For each *internal* k-mer of `seq` (not the first or last), verify global
    out-degree-with-matching-color == 1 in the unitig's walking orientation,
    and in-degree-with-matching-color == 1.

    Internal k-mer at position p has:
      - predecessor k-mer at position p-1 in the unitig (must be the unique
        in-edge with matching color)
      - successor k-mer at position p+1 (must be the unique out-edge)
    """
    L = len(seq)
    if L < k:
        return
    bs = bytes(seq, "ascii")

    # Build a list of (canonical_kmer, fwd_kmer_int) along the unitig.
    chain = []
    mask = (1 << (2 * k)) - 1
    rc_high_shift = 2 * (k - 1)
    fwd = 0
    rc = 0
    filled = 0
    table = _ACGT_TABLE
    for byte in bs:
        v = table[byte]
        if v == 0xff:
            problems.append(f"non-ACGT base in unitig")
            return
        fwd = ((fwd << 2) | v) & mask
        rc = (rc >> 2) | ((v ^ 3) << rc_high_shift)
        if filled + 1 < k:
            filled += 1
            continue
        can = fwd if fwd <= rc else rc
        chain.append((can, fwd, rc))

    # K-mers at chain positions 0..N-1. Internals are 1..N-2.
    for p in range(1, len(chain) - 1):
        cur_can, cur_fwd, cur_rc = chain[p]
        cur_color = kmer_colors.get(cur_can)
        # Out-degree in walking orientation (which is fwd of cur_fwd).
        out_count = 0
        for nt in range(4):
            ext = shift_append(cur_fwd, nt, k)
            ext_can = canonical_int(ext, k)
            if ext_can in kmer_colors and kmer_colors[ext_can] == cur_color:
                out_count += 1
        # In-degree from walking orientation: predecessors of cur_fwd are
        # the k-mers ending with cur_fwd's prefix. Equivalently, RC-side
        # successors of cur.
        in_count = 0
        for nt in range(4):
            ext = shift_append(cur_rc, nt, k)
            ext_can = canonical_int(ext, k)
            if ext_can in kmer_colors and kmer_colors[ext_can] == cur_color:
                in_count += 1
        if out_count != 1 or in_count != 1:
            problems.append(
                f"  internal k-mer at pos {p}: out={out_count} in={in_count}"
                f" kmer={seq[p:p+k]} unitig_len={len(seq)}")
            if len(problems) > 20:
                return


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filenames", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("-k", type=int, required=True)
    args = ap.parse_args()

    with open(args.filenames) as fh:
        files = [ln.strip() for ln in fh if ln.strip()]
    print(f"computing global k-mer color map from {len(files)} input(s)...", flush=True)
    kmer_colors = expected_colors_map(files, args.k)
    print(f"  {len(kmer_colors)} canonical k-mers", flush=True)

    fa_path = args.out + ".fa"
    print(f"checking {fa_path}...", flush=True)
    cur_seq = None
    n_unitigs = 0
    n_bad = 0
    problems = []
    with open(fa_path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if cur_seq:
                    n_unitigs += 1
                    pre = len(problems)
                    check_unitig_topology(cur_seq.decode("ascii"),
                                          args.k, kmer_colors, problems)
                    if len(problems) > pre:
                        n_bad += 1
                cur_seq = b""
            else:
                cur_seq = (cur_seq or b"") + line.rstrip(b"\r\n")
        if cur_seq:
            n_unitigs += 1
            pre = len(problems)
            check_unitig_topology(cur_seq.decode("ascii"),
                                  args.k, kmer_colors, problems)
            if len(problems) > pre:
                n_bad += 1

    print(f"unitigs: {n_unitigs}, with internal-branch problems: {n_bad}")
    for p in problems[:30]:
        print(p)
    sys.exit(0 if n_bad == 0 else 1)


if __name__ == "__main__":
    main()
