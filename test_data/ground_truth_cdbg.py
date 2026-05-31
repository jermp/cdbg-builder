#!/usr/bin/env python3
# Naive, slow, obviously-correct colored compacted de Bruijn graph builder.
#
# Purpose: an INDEPENDENT ground truth for the unitig count / unitig set, with
# no shared code with the C++ builder. Correctness over speed.
#
# Definitions (standard colored compacted dBG, a la BCALM2/GGCAT):
#   - Nodes: every canonical k-mer occurring in an ACGT-only run of any input.
#   - Color of a k-mer: the set of input files (genomes) containing it.
#   - Edge u->v: v == u[1:]+b for some base b (k-1 overlap), in the bidirected
#     graph (canonical k-mers, both strands).
#   - Unitig: a maximal path whose internal nodes are non-branching (in-degree
#     == out-degree == 1 in the bidirected dBG) AND that is monochromatic.
#     Equivalently: extend u->v iff out-deg(u)==1, in-deg(v)==1, color(u)==color(v).
#
# Usage:
#   ground_truth_cdbg.py --filenames files.txt -k 31 [--out truth.fa] [--compare other.fa]

import argparse
import gzip
import re
import sys

A2I = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
I2A = 'ACGT'

# Reverse-complement of 4 bases packed in one byte (bits: base0 = bits 0-1 ...
# base3 = bits 6-7). Output reverses base order within the byte and complements.
RC8 = [0] * 256
for _b in range(256):
    _v = 0
    for _i in range(4):
        _base = (_b >> (2 * _i)) & 3   # base0 first -> ends at high end (reversed)
        _v = (_v << 2) | (_base ^ 3)   # complement
    RC8[_b] = _v


def revcomp(x, k):
    """Reverse-complement a 2-bit-packed k-mer (k <= 32+)."""
    nbytes = (2 * k + 7) // 8
    r = 0
    for i in range(nbytes):
        r = (r << 8) | RC8[(x >> (8 * i)) & 0xFF]
    pad_bases = nbytes * 4 - k
    return r >> (2 * pad_bases)


def decode(x, k):
    return ''.join(I2A[(x >> (2 * (k - 1 - i))) & 3] for i in range(k))


def iter_runs(path):
    """Yield each ACGT-only run (uppercased) from a FASTA(.gz)."""
    op = gzip.open(path, 'rt') if path.endswith('.gz') else open(path)
    seq = []
    with op:
        for line in op:
            if line and line[0] == '>':
                if seq:
                    for run in re.split('[^ACGT]', ''.join(seq).upper()):
                        if run:
                            yield run
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            for run in re.split('[^ACGT]', ''.join(seq).upper()):
                if run:
                    yield run


def build_colors(files, k):
    """S: dict canonical_kmer_int -> color bitmask (bit i set if file i has it)."""
    S = {}
    for gi, path in enumerate(files):
        bit = 1 << gi
        n_run = 0
        for run in iter_runs(path):
            n_run += 1
            if len(run) < k:
                continue
            # roll a 2-bit k-mer along the run
            x = 0
            mask = (1 << (2 * k)) - 1
            for i, ch in enumerate(run):
                x = ((x << 2) | A2I[ch]) & mask
                if i + 1 >= k:
                    c = x
                    rc = revcomp(x, k)
                    if rc < c:
                        c = rc
                    S[c] = S.get(c, 0) | bit
        print(f"  [{gi}] {path}: runs={n_run} cumulative_kmers={len(S)}", file=sys.stderr)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--filenames', required=True)
    ap.add_argument('-k', type=int, required=True)
    ap.add_argument('--out')
    ap.add_argument('--compare')
    args = ap.parse_args()
    k = args.k
    files = [l.strip() for l in open(args.filenames) if l.strip()]
    print(f"building colored dBG from {len(files)} files, k={k}", file=sys.stderr)

    S = build_colors(files, k)
    ncolors = len(set(S.values()))
    print(f"canonical k-mers: {len(S)}", file=sys.stderr)
    print(f"distinct colors: {ncolors}", file=sys.stderr)

    # Membership set over BOTH strands so degree counting needs no rc.
    both = set()
    for c in S:
        both.add(c)
        both.add(revcomp(c, k))

    mask = (1 << (2 * k)) - 1
    hi = 2 * (k - 1)

    def canon(x):
        # exactly one of x / rc(x) is a key in S (k odd -> no self-rc)
        return x if x in S else revcomp(x, k)

    def succ_oriented(x):
        base = (x << 2) & mask
        return [base | b for b in range(4) if (base | b) in both]

    def pred_oriented(x):
        suf = x >> 2
        return [(b << hi) | suf for b in range(4) if ((b << hi) | suf) in both]

    def out_deg(x):
        return len({canon(n) for n in succ_oriented(x)})

    def in_deg(x):
        return len({canon(n) for n in pred_oriented(x)})

    # extend x -> unique successor y iff out_deg(x)==1, in_deg(y)==1, same color
    def fwd_step(x):
        s = succ_oriented(x)
        cs = {canon(n): n for n in s}
        if len(cs) != 1:
            return None
        cy, y = next(iter(cs.items()))
        if in_deg(y) != 1:
            return None
        if S[canon(x)] != S[cy]:
            return None
        return y

    # x can be extended to the LEFT (x has a unique unary same-color predecessor)
    def back_step(x):
        p = pred_oriented(x)
        cp = {canon(n): n for n in p}
        if len(cp) != 1:
            return None
        cz, z = next(iter(cp.items()))
        if out_deg(z) != 1:
            return None
        if S[canon(x)] != S[cz]:
            return None
        return z

    visited = set()
    unitigs = []
    for c0 in S:
        if c0 in visited:
            continue
        # walk backward (oriented) to the unitig start
        x = c0
        local = {c0}
        while True:
            z = back_step(x)
            if z is None:
                break
            cz = canon(z)
            if cz in local:  # circular
                break
            x = z
            local.add(cz)
        # forward walk from start, building sequence
        seq = decode(x, k)
        visited.add(canon(x))
        while True:
            y = fwd_step(x)
            if y is None:
                break
            cy = canon(y)
            if cy in visited:  # circular closure
                break
            seq += I2A[y & 3]
            visited.add(cy)
            x = y
        unitigs.append(seq)

    print(f"UNITIGS: {len(unitigs)}", file=sys.stderr)
    print(len(unitigs))

    if args.out:
        with open(args.out, 'w') as f:
            for i, u in enumerate(unitigs):
                f.write(f">{i}\n{u}\n")

    if args.compare:
        def cseq(s):
            r = ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[c] for c in reversed(s))
            return s if s <= r else r
        truth = set(cseq(u) for u in unitigs)
        other = []
        s = []
        for line in open(args.compare):
            if line[0] == '>':
                if s:
                    other.append(''.join(s)); s = []
            else:
                s.append(line.strip())
        if s:
            other.append(''.join(s))
        oset = set(cseq(u) for u in other)
        print(f"compare {args.compare}: truth={len(truth)} other={len(oset)} "
              f"shared={len(truth & oset)} only_truth={len(truth - oset)} "
              f"only_other={len(oset - truth)}", file=sys.stderr)


if __name__ == '__main__':
    main()
