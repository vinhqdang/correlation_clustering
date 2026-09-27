#!/usr/bin/env bash
# Multicut dual bound of RAMA (GPU) on the distance-two support of the SNAP
# graphs, as a lower-bound baseline.  Ubuntu or WSL2 with an NVIDIA GPU.
#
#   conda activate cc                       # the same environment as for KaPoCE
#   conda install -y -c nvidia cuda-toolkit # nvcc, if `nvcc --version` fails
#   bash experiments/rama/run_rama.sh
#
#   RAMA_CPU=1 bash experiments/rama/run_rama.sh    (CPU solver, no GPU; slow)
#   GRAPHS="ca-GrQc" bash experiments/rama/run_rama.sh   (a quick test)
#
# Needs git, cmake >= 3.22.1, g++ >= 10 and, for the GPU solver, CUDA >= 11.2.
# The results go to results/rama_snap_user.csv (resumable).  Commit and push it.
set -euo pipefail
cd "$(dirname "$0")/../.."
COMMIT=ee459d6cc0c51b56f82690e91ddf769064e2a47e     # pawelswoboda/RAMA
DIR="$PWD/rama_build"
OUT=results/rama_snap_user.csv

if [ -z "${VIRTUAL_ENV:-}" ] && [ -z "${CONDA_PREFIX:-}" ]; then
    [ -d .venv ] || python3 -m venv .venv
    # shellcheck disable=SC1091
    . .venv/bin/activate
fi
python3 -m pip install -q -e .

if [ -n "${RAMA_CPU:-}" ]; then
    CUDA=OFF; TARGET=rama_text_input_cpu
else
    CUDA=ON; TARGET=rama_text_input_gpu
    command -v nvcc > /dev/null || {
        echo "nvcc not found: conda install -y -c nvidia cuda-toolkit (or RAMA_CPU=1)" >&2; exit 1; }
    nvidia-smi -L || { echo "no GPU visible (nvidia-smi failed)" >&2; exit 1; }
fi

if [ ! -d "$DIR/.git" ]; then
    git clone https://github.com/pawelswoboda/RAMA.git "$DIR"
fi
git -C "$DIR" fetch -q origin || true
git -C "$DIR" checkout -q "$COMMIT"
mkdir -p "$DIR/build_$CUDA" && cd "$DIR/build_$CUDA"
EXTRA=()
if [ "$CUDA" = ON ] && [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/nvcc" ]; then
    EXTRA=(-DCMAKE_CUDA_COMPILER="$CONDA_PREFIX/bin/nvcc")
fi
cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=$CUDA "${EXTRA[@]}" > cmake.log 2>&1 || {
    tail -n 30 cmake.log >&2; echo "cmake failed, see $PWD/cmake.log" >&2; exit 1; }
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)" "$TARGET" > make.log 2>&1 || {
    tail -n 30 make.log >&2; echo "build failed, see $PWD/make.log" >&2; exit 1; }
BIN="$(find "$PWD" -type f -name "$TARGET" -perm -u+x | head -n1)"
cd - > /dev/null
echo "RAMA: $BIN"

read -r -a G <<< "${GRAPHS:-}"
python3 experiments/rama/rama_snap.py "$BIN" "$OUT" "${G[@]}"
echo "done; now: git add $OUT && git commit -m 'RAMA bounds on SNAP graphs' && git push"
