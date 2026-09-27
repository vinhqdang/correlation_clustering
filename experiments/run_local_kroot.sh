#!/usr/bin/env bash
# KaPoCE root bounds (mode kroot, 900 s) on the PACE exact instances that have
# no result yet in results/colab; resumable, one instance at a time.
cd "$(dirname "$0")/.."
export CC_ROOT="$PWD" RESULTS_DIR="$PWD/results/colab" KAPOCE_LB_DIR="$PWD/kapoce_lb"
for i in $(seq 1 200); do
    g=$(printf "pace-exact/exact%03d.gr" "$i")
    f="results/colab/lkroot_pace-exact_exact$(printf %03d "$i").gr_0.json"
    [ -f "$f" ] && continue
    python3 experiments/colab_run_lp.py 900 0 lkroot "$g" > /dev/null 2>&1
    n=$((${n:-0} + 1))
    # commit every few results, so that they are not left untracked
    if [ $((n % 5)) -eq 0 ]; then
        git add results/colab/lkroot_pace-exact_* && \
        git -c user.name=Vinh -c user.email=dqvinh87@gmail.com commit -q \
            -m "KaPoCE root bounds on PACE exact instances" -- results/colab > /dev/null 2>&1 && \
        { git pull -q --rebase --autostash origin main && git push -q origin main; } > /dev/null 2>&1
    fi
done
git add results/colab/lkroot_pace-exact_* && \
    git -c user.name=Vinh -c user.email=dqvinh87@gmail.com commit -q \
        -m "KaPoCE root bounds on PACE exact instances" -- results/colab > /dev/null 2>&1
{ git pull -q --rebase --autostash origin main && git push -q origin main; } > /dev/null 2>&1
touch results/colab/lkroot_pace.done
