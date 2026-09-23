#!/bin/bash
# Supervisor: keeps the fleet manager running.  It is restarted when it exits
# and killed and restarted when its heartbeat is older than 20 minutes.
cd "$(dirname "$0")/.."
HB=results/colab/heartbeat
while true; do
    date +%s > "$HB"
    python3 experiments/colab_fleet.py >> results/colab/fleet.out 2>&1 &
    PID=$!
    while kill -0 "$PID" 2> /dev/null; do
        sleep 60
        if [ $(( $(date +%s) - $(cut -d. -f1 "$HB") )) -gt 1200 ]; then
            echo "$(date '+%F %T') heartbeat stale, killing fleet manager" >> results/colab/fleet.log
            kill -9 "$PID"
        fi
    done
    echo "$(date '+%F %T') fleet manager exited, restarting" >> results/colab/fleet.log
    sleep 10
done
