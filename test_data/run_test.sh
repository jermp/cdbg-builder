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

echo "=== in-memory path ==="
OUT_MEM="$WORK/sal_mem"
"$BIN" -i "$FILELIST" -k "$K" -o "$OUT_MEM" -t "$THREADS"
python3 "$HERE/verify.py" --filenames "$FILELIST" --out "$OUT_MEM" -k "$K"

echo
echo "=== bucketed path ==="
OUT_BKT="$WORK/sal_bkt"
"$BIN" -i "$FILELIST" -k "$K" -o "$OUT_BKT" -t "$THREADS" --bucketed --tmp-dir "$WORK/buckets"
python3 "$HERE/verify.py" --filenames "$FILELIST" --out "$OUT_BKT" -k "$K"
