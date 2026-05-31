#!/usr/bin/env python3
"""
Authoritative correctness check for cdgb-build output.

This is the single source of truth. It builds the colored compacted de Bruijn
graph (ccdBG) independently, the naive obviously-correct way, and compares the
builder's output against it EXACTLY. Slow but correct -- correctness over speed.

Ground-truth definition (faithful to BCALM2 / GGCAT, which this tool targets):
  - Nodes: every canonical k-mer occurring in an ACGT-only run of any input.
  - Color of a k-mer: the bitmask of input files (genomes) containing it.
  - A unitig is a maximal path that is NON-BRANCHING in the underlying
    (uncolored) bidirected dBG AND monochromatic. Equivalently: extend x -> y
    iff x has a unique canonical successor y, y has a unique canonical
    predecessor, and color(x) == color(y). Splits happen at uncolored branches
    AND at color changes.

What it checks (a single exact-set comparison subsumes the first four):
  1. Coverage         -- every input k-mer is in exactly one output unitig.
  2. Monochromaticity -- every k-mer in a unitig shares one color set.
  3. Topology         -- no internal branch (no degree-2 interior node).
  4. MAXIMALITY       -- unitigs are not splittable further at simple-path
                         points (the check verify.py / strict_topology missed:
                         a pile of short fragments can pass coverage +
                         monochromaticity + colored-degree topology yet not be
                         maximal). Exact set equality catches over- AND under-
                         joining.
  5. Color-set ids    -- each .fa header cid maps to one consistent color set;
                         #distinct cids == #distinct color sets; the emitted
                         set of color sets equals the input's.
  6. .color_sets file -- present and non-empty.

Why exact set equality is authoritative: the ground-truth set IS the definition
of a correct ccdBG, so output_set == truth_set <=> the output is a correct
ccdBG. (strict_topology_check.py used COLORED degree for its branch test, which
is too weak -- it passed both an over-split output and the correct one. This
replaces it.)

Usage:
    verify.py --filenames <list> --out <basename> -k <k>
Exit 0 on success; non-zero with a message on the first failure.
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
_NUC = "ACGT"

# Reverse-complement of 4 bases packed in one byte (base0 = bits 0-1 ...).
_RC8 = [0] * 256
for _b in range(256):
    _v = 0
    for _i in range(4):
        _v = (_v << 2) | (((_b >> (2 * _i)) & 3) ^ 3)
    _RC8[_b] = _v


def _open(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def stream_records(path):
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


def revcomp_int(x, k):
    """Reverse-complement a 2-bit-packed k-mer using the byte table."""
    nbytes = (2 * k + 7) // 8
    r = 0
    for i in range(nbytes):
        r = (r << 8) | _RC8[(x >> (8 * i)) & 0xFF]
    return r >> (nbytes * 4 * 2 - 2 * k)


def decode_int(x, k):
    return "".join(_NUC[(x >> (2 * (k - 1 - i))) & 3] for i in range(k))


def canonical_kmers(seq, k):
    """Yield 2-bit packed canonical k-mer ints; restarts at any non-ACGT base."""
    mask = (1 << (2 * k)) - 1
    rc_high_shift = 2 * (k - 1)
    fwd = rc = filled = 0
    table = _ACGT_TABLE
    for byte in seq:
        v = table[byte]
        if v == 0xff:
            filled = fwd = rc = 0
            continue
        fwd = ((fwd << 2) | v) & mask
        rc = (rc >> 2) | ((v ^ 3) << rc_high_shift)
        if filled + 1 < k:
            filled += 1
            continue
        yield fwd if fwd <= rc else rc


def build_color_map(files, k):
    """{canonical_kmer_int: color_bitmask_int} over all inputs."""
    S = {}
    for color, path in enumerate(files):
        bit = 1 << color
        for seq in stream_records(path):
            for km in canonical_kmers(seq, k):
                S[km] = S.get(km, 0) | bit
    return S


def build_truth_unitigs(S, k):
    """
    Build the ground-truth ccdBG unitig set from the color map S.

    Returns {canonical_unitig_string: color_bitmask}. A unitig is a maximal
    monochromatic non-branching path (uncolored-branch + color-split).
    """
    mask = (1 << (2 * k)) - 1
    hi = 2 * (k - 1)

    # Membership over BOTH strands so neighbor existence needs no canonicalize.
    both = set()
    for c in S:
        both.add(c)
        both.add(revcomp_int(c, k))

    def canon(x):
        return x if x in S else revcomp_int(x, k)

    def succ(x):
        base = (x << 2) & mask
        return [base | b for b in range(4) if (base | b) in both]

    def pred(x):
        suf = x >> 2
        return [(b << hi) | suf for b in range(4) if ((b << hi) | suf) in both]

    def out_deg(x):
        return len({canon(n) for n in succ(x)})

    def in_deg(x):
        return len({canon(n) for n in pred(x)})

    def fwd_step(x):
        cs = {canon(n): n for n in succ(x)}
        if len(cs) != 1:
            return None
        cy, y = next(iter(cs.items()))
        if in_deg(y) != 1 or S[canon(x)] != S[cy]:
            return None
        return y

    def back_step(x):
        cp = {canon(n): n for n in pred(x)}
        if len(cp) != 1:
            return None
        cz, z = next(iter(cp.items()))
        if out_deg(z) != 1 or S[canon(x)] != S[cz]:
            return None
        return z

    def cseq(s):
        r = "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[c] for c in reversed(s))
        return s if s <= r else r

    visited = set()
    truth = {}
    for c0 in S:
        if c0 in visited:
            continue
        # Walk back (oriented) to the unitig start.
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
        # Forward walk from the start, building the sequence.
        seq = decode_int(x, k)
        col = S[canon(x)]
        visited.add(canon(x))
        while True:
            y = fwd_step(x)
            if y is None:
                break
            cy = canon(y)
            if cy in visited:  # circular closure
                break
            seq += _NUC[y & 3]
            visited.add(cy)
            x = y
        truth[cseq(seq)] = col
    return truth


def canon_str(s):
    r = "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[c] for c in reversed(s))
    return s if s <= r else r


def verify(filenames_list, out_basename, k):
    with open(filenames_list) as fh:
        files = [ln.strip() for ln in fh if ln.strip()]

    print(f"[1/3] building ground-truth color map from {len(files)} input(s)...",
          flush=True)
    S = build_color_map(files, k)
    expected_classes = set(S.values())
    print(f"      {len(S)} canonical k-mers, {len(expected_classes)} distinct "
          f"color set(s)", flush=True)

    print("[2/3] building ground-truth unitig set (naive ccdBG)...", flush=True)
    truth = build_truth_unitigs(S, k)
    print(f"      {len(truth)} ground-truth unitig(s)", flush=True)

    fa_path = out_basename + ".fa"
    colors_path = out_basename + ".color_sets"
    if not os.path.exists(fa_path):
        print(f"FAIL: {fa_path} not found")
        return 1
    if not os.path.exists(colors_path) or os.path.getsize(colors_path) == 0:
        print(f"FAIL: {colors_path} missing or empty")
        return 1

    print(f"[3/3] checking {fa_path} against ground truth...", flush=True)
    out_set = {}          # canonical unitig seq -> header cid
    id_to_class = {}      # header cid -> color bitmask (from ground truth)
    dup = 0
    cid = None
    body = []

    def flush(cid, seq_bytes):
        nonlocal dup
        if cid is None or not seq_bytes:
            return 0
        seq = seq_bytes.decode("ascii")
        cu = canon_str(seq)
        if cu in out_set:
            dup += 1
            if dup <= 5:
                print(f"FAIL: unitig emitted more than once (cid={cid})")
            return 1
        out_set[cu] = cid
        col = truth.get(cu)
        if col is None:
            print(f"FAIL: output unitig (cid={cid}, len={len(seq)}) is not a "
                  f"ground-truth unitig -- not maximal / wrong split / spurious")
            return 1
        prev = id_to_class.get(cid)
        if prev is None:
            id_to_class[cid] = col
        elif prev != col:
            print(f"FAIL: color_set_id {cid} maps to two color sets "
                  f"{prev:b} and {col:b}")
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
                    print(f"FAIL: non-integer FASTA header {line[1:].rstrip()!r}")
                    return 1
                body = []
            else:
                body.append(line.rstrip(b"\r\n"))
        rc = flush(cid, b"".join(body))
        if rc:
            return rc

    # 1-4) Exact unitig-set equality (coverage + monochromaticity + topology +
    #      maximality, all at once).
    if len(out_set) != len(truth):
        only_truth = len(set(truth) - set(out_set))
        only_out = len(set(out_set) - set(truth))
        print(f"FAIL: unitig set differs from ground truth: "
              f"output={len(out_set)} truth={len(truth)} "
              f"(missing from output={only_truth}, spurious in output={only_out})")
        return 1
    # equal sizes + every output unitig is in truth + no dups => sets are equal.

    # 5) Color-set id semantics.
    if set(id_to_class.values()) != expected_classes:
        print("FAIL: set of emitted color sets differs from the input's")
        return 1
    if len(set(id_to_class.keys())) != len(expected_classes):
        print(f"FAIL: distinct color_set_ids = {len(id_to_class)}, "
              f"expected {len(expected_classes)}")
        return 1

    print(f"OK: {len(out_set)} unitig(s) == ground truth, "
          f"{len(S)} canonical k-mer(s), {len(expected_classes)} color set(s); "
          f"{colors_path} present ({os.path.getsize(colors_path)} bytes)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--filenames", required=True,
                    help="filenames list passed to cdgb-build (-i)")
    ap.add_argument("--out", required=True,
                    help="output basename passed to cdgb-build (-o)")
    ap.add_argument("-k", type=int, required=True)
    args = ap.parse_args()
    sys.exit(verify(args.filenames, args.out, args.k))


if __name__ == "__main__":
    main()
