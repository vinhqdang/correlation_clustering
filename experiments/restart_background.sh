#!/bin/bash
# Restart the long-running background jobs after a container restart: the
# fleet supervisor and the resumable recheck.
cd "$(dirname "$0")/.."
(setsid nohup bash experiments/colab_fleet.sh > /dev/null 2>&1 &)
run() {  # name, command...
    ps -C python3 -o args= | grep -qF "$1" && return
    (nohup setsid "${@:2}" > "/tmp/$(basename "$1").out" 2>&1 < /dev/null &)
}
# resumable rechecks of the Colab and the local certificates; each result is
# copied to its committed file when complete
recheck() {  # certificate directory, name
    local w="results/mpc/$2_work.csv"
    [ -f "results/mpc/$2.done" ] && return
    run "$w" bash -c "python3 experiments/check_certificate.py --resume data/raw $1 $w \
        && cp $w results/mpc/$2.csv && touch results/mpc/$2.done"
}
recheck results/certificates/colab recheck
recheck results/local/certs recheck_local
# the KaPoCE root bounds run as fleet jobs (tag lkroot)
# KaPoCE root bounds on the PACE exact instances (local, resumable)
if [ ! -f results/colab/lkroot_pace.done ]; then
    ps -eo args | grep -q "[r]un_local_kroot.sh" || \
        (nohup setsid bash experiments/run_local_kroot.sh > /tmp/run_local_kroot.out 2>&1 < /dev/null &)
fi
