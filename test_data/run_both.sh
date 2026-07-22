#!/usr/bin/env bash
# Correctness + timing harness for a single dataset directory.
# Usage: run_both.sh <dataset_dir> <glob> [k] [threads] [extra cdbg-build args...]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="$REPO/build/cdbg-build"
DIR="$1"; GLOB="$2"; K="${3:-31}"; THREADS="${4:-4}"
shift 4 || shift $#
EXTRA=("$@")

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
FILELIST="$WORK/filenames.txt"
# shellcheck disable=SC2086
ls $DIR/$GLOB | grep -v '/\._' > "$FILELIST"
echo "== dataset $DIR : $(wc -l < "$FILELIST") files, k=$K, t=$THREADS, extra=${EXTRA[*]:-none} =="
OUT="$WORK/out"
"$BIN" -i "$FILELIST" -k "$K" -o "$OUT" -t "$THREADS" -d "$WORK/buckets" --verbose "${EXTRA[@]}" 2>&1 \
  | grep -E '\[stitch\]|stitch peak|unitigs after|total construction|peak resident' || true
python3 "$HERE/verify.py" --filenames "$FILELIST" --out "$OUT" -k "$K"
