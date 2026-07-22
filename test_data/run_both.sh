#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"; BIN="$REPO/build/cdbg-build"
DIR="$1"; GLOB="$2"; K="${3:-31}"; THREADS="${4:-4}"; shift 4 || shift $#; EXTRA=("$@")
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; FL="$WORK/fl.txt"
ls $DIR/$GLOB | grep -v '/\._' > "$FL"
echo "== $DIR : $(wc -l < "$FL") files k=$K t=$THREADS extra=${EXTRA[*]:-none} =="
"$BIN" -i "$FL" -k "$K" -o "$WORK/out" -t "$THREADS" -d "$WORK/b" --verbose "${EXTRA[@]}" 2>&1 \
 | grep -E 'bucket-write\]|bucket-write peak|hashmap|spill|lock_wait|bucket bytes|total constr|peak resident' || true
python3 "$HERE/verify.py" --filenames "$FL" --out "$WORK/out" -k "$K"
