#!/usr/bin/env bash
# Re-check every archived certificate (and the cost of every archived
# clustering) against freshly downloaded instance files, from scratch.
#
#   bash experiments/recheck_all.sh [--resume] [OUT_DIR]
#
# Writes OUT_DIR/recheck.csv (results/certificates/colab) and
# OUT_DIR/recheck_local.csv (results/local/certs), default OUT_DIR
# results/mpc/fresh, and then compares every accepted value with the committed
# results/mpc/recheck*.csv.  Without --resume every certificate is checked
# again, whatever the output files contain; with --resume the rows already in
# the output files are kept, so that an interrupted run can be continued.
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=--dir
if [ "${1:-}" = "--resume" ]; then
    MODE=--resume
    shift
fi
OUT="${1:-results/mpc/fresh}"
mkdir -p "$OUT"
python experiments/fetch_data.py
python experiments/check_certificate.py "$MODE" data/raw results/certificates/colab "$OUT/recheck.csv"
python experiments/check_certificate.py "$MODE" data/raw results/local/certs "$OUT/recheck_local.csv"
python - "$OUT" <<'EOF'
import csv, sys
out = sys.argv[1]
bad = 0
for name in ("recheck.csv", "recheck_local.csv"):
    new = {r["certificate"]: r for r in csv.DictReader(open(f"{out}/{name}"))}
    old = {r["certificate"]: r for r in csv.DictReader(open(f"results/mpc/{name}"))}
    for c, r in new.items():
        o = old.get(c)
        if r["status"] != "ok" or o is None or (r["certified"], r["cost"]) != (o["certified"], o["cost"]):
            bad += 1
            print("differs:", name, c, r["status"], r["certified"], r["cost"])
    missing = set(old) - set(new)
    bad += len(missing)
    for c in sorted(missing):
        print("not checked:", name, c)
    print(f"{name}: {len(new)} checked, {sum(r['status'] == 'ok' for r in new.values())} accepted")
print("all values agree with results/mpc/" if bad == 0 else f"{bad} differences")
sys.exit(1 if bad else 0)
EOF
