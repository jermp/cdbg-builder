#!/usr/bin/env bash
# End-to-end correctness test.
# Builds a colored compacted dBG over test_data/salmonella_10/*.fasta.gz with
# k=31, then runs verify.py to validate every k-mer / unitig / color-set id.
#
# Run from anywhere; binary is expected at <repo>/build/cdgb-build.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="$REPO/build/cdgb-build"
K=31
THREADS="${THREADS:-4}"
DATA_DIR="$HERE/salmonella_10"

if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found. Build the project first:" >&2
    echo "  mkdir -p $REPO/build && cd $REPO/build && cmake .. && make -j" >&2
    exit 2
fi

shopt -s nullglob
inputs=("$DATA_DIR"/*.fasta.gz)
if (( ${#inputs[@]} == 0 )); then
    echo "error: no *.fasta.gz files in $DATA_DIR" >&2
    exit 2
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

FILELIST="$WORK/filenames.txt"
printf '%s\n' "${inputs[@]}" > "$FILELIST"

OUT="$WORK/sal"
mkdir -p "$WORK/buckets"
"$BIN" -i "$FILELIST" -k "$K" -o "$OUT" -t "$THREADS" --tmp-dir "$WORK/buckets"
python3 "$HERE/verify.py" --filenames "$FILELIST" --out "$OUT" -k "$K"
python3 "$HERE/strict_topology_check.py" --filenames "$FILELIST" --out "$OUT" -k "$K"
