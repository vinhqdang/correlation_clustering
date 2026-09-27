#!/bin/bash
# Restart the long-running background jobs after a container restart: the
# fleet supervisor and the resumable recheck.
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
# the KaPoCE root bounds run as fleet jobs (tag lkroot)
