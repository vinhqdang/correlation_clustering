#!/usr/bin/env bash
# Re-check every archived certificate (and the cost of every archived
# clustering) against freshly downloaded instance files.
#   bash experiments/recheck_all.sh [CERT_DIR]      (default results/certificates/colab)
set -euo pipefail
cd "$(dirname "$0")/.."
python experiments/fetch_data.py
python experiments/check_certificate.py --dir data/raw "${1:-results/certificates/colab}" \
    results/mpc/recheck.csv
