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
done
touch results/colab/lkroot_pace.done
