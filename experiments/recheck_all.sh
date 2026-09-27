#!/usr/bin/env bash
# Re-check every archived certificate (and the cost of every archived
# clustering) against freshly downloaded instance files.
#   bash experiments/recheck_all.sh [CERT_DIR]
# (default: results/certificates/colab and results/local/certs)
set -euo pipefail
cd "$(dirname "$0")/.."
python experiments/fetch_data.py
if [ -n "${1:-}" ]; then
    python experiments/check_certificate.py --resume data/raw "$1" results/mpc/recheck.csv
else
    python experiments/check_certificate.py --resume data/raw results/certificates/colab \
        results/mpc/recheck.csv
    python experiments/check_certificate.py --resume data/raw results/local/certs \
        results/mpc/recheck_local.csv
fi
