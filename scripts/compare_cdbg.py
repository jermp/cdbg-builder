#!/usr/bin/env python3
"""Compare two colored-cdBG builds for semantic equality.

The .fa output is emitted in cid-ascending order, header = color-class id, but
across branches/thread-counts the unitig ORDER, the unitig ORIENTATION (a unitig
vs its reverse-complement), and the cid NUMBERING can all differ legitimately.
This compares the two builds invariant to all three, in a single streaming pass:

  * unitig set     -- order-independent multiset hash of canonical(seq)
  * color partition-- multiset hash of per-color-class member sets (renumber-safe)

Two builds are semantically identical iff both accumulators match. Reads only the
.fa files (cid is in the header), so no .color_sets decode is required.

Usage:  compare_cdbg.py a.fa b.fa
Exit code 0 if identical, 1 if different.
"""
import sys
import hashlib

_COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def revcomp(s: bytes) -> bytes:
    return s.translate(_COMP)[::-1]


def canonical(s: bytes) -> bytes:
    rc = revcomp(s)
    return s if s <= rc else rc


def h128(b: bytes) -> int:
    # Deterministic 128-bit hash (stable across processes/machines, so the
    # printed hashes are reproducible and can be recorded as a golden value).
    # Summing these over a set is order-independent.
    return int.from_bytes(hashlib.blake2b(b, digest_size=16).digest(), "little")


MASK = (1 << 128) - 1


def fingerprint(path: str):
    """Return (n_unitigs, n_classes, unitig_acc, partition_acc)."""
    unitig_acc = 0          # sum of h128(canonical(seq)) over all unitigs
    partition_acc = 0       # sum of h128(group fingerprint) over all color classes
    n_unitigs = 0
    n_classes = 0
    cur_cid = None
    group_acc = 0           # sum of h128(canonical(seq)) within the current cid run
    group_n = 0

    def close_group():
        nonlocal partition_acc, n_classes
        if cur_cid is None:
            return
        # Fold the group's id-free member set + its size into one fingerprint.
        sig = h128(b"%d:%d" % (group_n, group_acc & MASK))
        partition_acc = (partition_acc + sig) & MASK
        n_classes += 1

    cid = None
    seq_parts = []

    def flush_record():
        nonlocal unitig_acc, n_unitigs, cur_cid, group_acc, group_n
        if cid is None:
            return
        c = canonical(b"".join(seq_parts))
        hc = h128(c)
        unitig_acc = (unitig_acc + hc) & MASK
        n_unitigs += 1
        if cid != cur_cid:
            close_group()
            cur_cid = cid
            group_acc = 0
            group_n = 0
        group_acc = (group_acc + hc) & MASK
        group_n += 1

    with open(path, "rb") as f:
        for line in f:
            line = line.rstrip(b"\n")
            if not line:
                continue
            if line[:1] == b">":
                flush_record()
                cid = int(line[1:].split()[0])  # header is the integer cid
                seq_parts = []
            else:
                seq_parts.append(line)
        flush_record()
        close_group()

    return n_unitigs, n_classes, unitig_acc, partition_acc


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    a = fingerprint(argv[1])
    b = fingerprint(argv[2])
    names = ("unitigs", "color_classes", "unitig_set_hash", "color_partition_hash")
    ok = True
    for name, va, vb in zip(names, a, b):
        same = va == vb
        ok &= same
        sa = va if name.startswith(("unitigs", "color_classes")) else f"{va:032x}"
        sb = vb if name.startswith(("unitigs", "color_classes")) else f"{vb:032x}"
        print(f"{'OK ' if same else 'DIFF'}  {name:22} {sa}  {sb}")
    print("IDENTICAL" if ok else "DIFFERENT")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
