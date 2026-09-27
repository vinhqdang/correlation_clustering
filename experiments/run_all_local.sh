#!/usr/bin/env bash
# All jobs meant for a local workstation, at the same time:
#   1. KaPoCE root bounds on the ten small SNAP graphs   (CPU, experiments/kapoce)
#   2. RAMA multicut bounds on the 18 SNAP graphs        (GPU, experiments/rama)
#   3. PACE exact-track bound jobs, 1800 s bound and LP   (CPU, run_local_bounds.py)
#
#   conda activate cc          # Python >= 3.10, cmake >= 3.22, nvcc (CUDA toolkit)
#   bash experiments/run_all_local.sh
#
# Defaults suit 24 cores / 64 GB / one GPU: KAPOCE_JOBS=3, PACE_JOBS=16.
# Each part logs to logs_local/<part>.log and prints its lines with a prefix;
# everything is resumable: after an interruption, run the script again.
# At the end the results are committed and pushed (if git can push).
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
KAPOCE_JOBS="${KAPOCE_JOBS:-3}"
PACE_JOBS="${PACE_JOBS:-16}"
mkdir -p logs_local

# scripts checked out with Windows line endings do not run in bash
find experiments -name "*.sh" -exec sed -i 's/\r$//' {} +

say() { echo "[$(date +%H:%M:%S)] $*"; }
fail=0
command -v cmake > /dev/null || { say "cmake not found (sudo apt install cmake, or pip install cmake)"; fail=1; }
command -v g++ > /dev/null || { say "g++ not found (sudo apt install g++ make)"; fail=1; }
if [ -z "${RAMA_CPU:-}" ]; then
    command -v nvcc > /dev/null || { say "nvcc not found (conda install -y -c nvidia cuda-toolkit)"; fail=1; }
    nvidia-smi -L > /dev/null 2>&1 || { say "no GPU visible (nvidia-smi failed)"; fail=1; }
fi
python3 -c 'import sys; sys.exit(sys.version_info < (3, 10))' || { say "Python >= 3.10 needed"; fail=1; }
[ "$fail" = 0 ] || exit 1

say "installing the Python package once, before the parallel parts"
if [ -z "${VIRTUAL_ENV:-}" ] && [ -z "${CONDA_PREFIX:-}" ]; then
    [ -d .venv ] || python3 -m venv .venv
    . .venv/bin/activate
fi
python3 -m pip install -q -e . || exit 1

run_part() {  # name, command...
    local name="$1"; shift
    ( "$@" 2>&1 | tee -a "logs_local/$name.log" | sed -u "s/^/[$name] /" ) &
}

say "starting: kapoce ($KAPOCE_JOBS jobs), rama (GPU), pace ($PACE_JOBS jobs)"
say "logs: $ROOT/logs_local/{kapoce,rama,pace}.log"
export SKIP_PIP=1
run_part kapoce env JOBS="$KAPOCE_JOBS" bash experiments/kapoce/run_root_bounds.sh
P1=$!
run_part rama bash experiments/rama/run_rama.sh
P2=$!
# shellcheck disable=SC2086
run_part pace python3 experiments/run_local_bounds.py --jobs "$PACE_JOBS" ${PACE_ARGS:-}
P3=$!
st=0
for p in $P1 $P2 $P3; do wait "$p" || st=1; done
say "all parts finished (exit status $st; see logs_local/ for failures)"

[ -z "${NO_COMMIT:-}" ] || exit "$st"
git add results/kapoce_root_snap_user.csv results/rama_snap_user.csv results/local 2> /dev/null
if git commit -q -m "Local runs: KaPoCE root bounds, RAMA bounds, PACE bound jobs"; then
    git pull -q --rebase --autostash origin main && git push -q origin main \
        && say "results committed and pushed" \
        || say "committed; push failed, run: git pull --rebase && git push"
else
    say "nothing new to commit"
fi
