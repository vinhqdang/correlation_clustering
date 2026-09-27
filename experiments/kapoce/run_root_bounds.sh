#!/usr/bin/env bash
# Root star-packing and P3-packing bounds of the KaPoCE branch-and-bound on the
# SNAP graphs with at most 2.5e4 vertices (Table "SNAP bounds", column B&B root).
# CPU only; up to about 6 GB of memory per job on the largest graph (ca-CondMat).
#
#   bash experiments/kapoce/run_root_bounds.sh            (from the repository root)
#   JOBS=1 bash experiments/kapoce/run_root_bounds.sh     (less memory)
#   GRAPHS="ca-GrQc" bash experiments/kapoce/run_root_bounds.sh   (a quick test)
#
# The results go to results/kapoce_root_snap_user.csv; the script is resumable
# (graphs already in the file are skipped).  Commit and push that file.
set -euo pipefail
cd "$(dirname "$0")/../.."
JOBS="${JOBS:-2}"
OUT=results/kapoce_root_snap_user.csv

python3 -m pip install -q -e .
BIN="$(bash experiments/kapoce/build_lbounds.sh "$PWD/kapoce_lb" | tail -n1)"
export KAPOCE_LB_BIN="$BIN"
echo "lbounds: $BIN"

# smallest graphs first; each job appends to its own file, merged at the end
if [ -n "${GRAPHS:-}" ]; then read -r -a GRAPHS <<< "$GRAPHS"; else
GRAPHS=(email-Eu-core BitcoinAlpha+ facebook ca-GrQc BitcoinOTC+ wiki-Vote ca-HepTh
        ca-HepPh ca-AstroPh ca-CondMat)
fi
pids=()
for j in $(seq 0 $((JOBS - 1))); do
    mine=()
    for i in "${!GRAPHS[@]}"; do
        [ $((i % JOBS)) -eq "$j" ] && mine+=("${GRAPHS[$i]}")
    done
    python3 experiments/kapoce_root_snap.py "$OUT.part$j" "${mine[@]}" &
    pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

python3 - "$OUT" <<'PY'
import csv, glob, sys
out = sys.argv[1]
best = {}
for f in sorted(glob.glob(out + ".part*")):
    for r in csv.DictReader(open(f)):
        best[(r["graph"], r["bound"])] = r
rows = list(best.values())
keys = ["graph", "n", "m", "bound", "lb", "time", "cpu", "limit"]
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} rows written to {out}")
PY
echo "done; now: git add $OUT && git commit -m 'KaPoCE root bounds on SNAP graphs' && git push"
