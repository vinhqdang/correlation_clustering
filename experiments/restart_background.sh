#!/bin/bash
# Restart the long-running background jobs after a container restart: the
# fleet supervisor, the resumable recheck and the KaPoCE root bounds.
cd "$(dirname "$0")/.."
(setsid nohup bash experiments/colab_fleet.sh > /dev/null 2>&1 &)
run() {  # name, command...
    ps -C python3 -o args= | grep -qF "$1" && return
    (nohup setsid "${@:2}" > "/tmp/$(basename "$1").out" 2>&1 < /dev/null &)
}
if [ ! -f results/mpc/recheck.done ]; then
    run recheck_work.csv python3 experiments/check_certificate.py --resume data/raw \
        results/certificates/colab results/mpc/recheck_work.csv
fi
run kapoce_root_snap_a.csv python3 experiments/kapoce_root_snap.py results/kapoce_root_snap_a.csv \
    email-Eu-core ca-GrQc BitcoinAlpha+ wiki-Vote
run kapoce_root_snap_b.csv python3 experiments/kapoce_root_snap.py results/kapoce_root_snap_b.csv \
    BitcoinOTC+ ca-HepTh facebook
run kapoce_root_snap_c.csv python3 experiments/kapoce_root_snap.py results/kapoce_root_snap_c.csv \
    ca-HepPh ca-AstroPh ca-CondMat
