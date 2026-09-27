#!/usr/bin/env bash
# Build `lbounds`, the root lower bounds (P3 packing and star packing) of the
# KaPoCE branch-and-bound, from the `experiments` branch of KaPoCE.
#
#   bash experiments/kapoce/build_lbounds.sh [BUILD_DIR]     (default: ./kapoce_lb)
#
# Needs git, cmake (>= 3.16), a C++17 compiler (GCC or Clang), make and Boost
# program_options (sudo apt install -y libboost-program-options-dev); no GPU.
# Prints the path of the binary at the end.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DIR="${1:-$PWD/kapoce_lb}"
COMMIT=63079a9          # branch `experiments` of kittobi1992/cluster_editing

if [ ! -d "$DIR/.git" ]; then
    git clone --recursive https://github.com/kittobi1992/cluster_editing.git "$DIR"
fi
cd "$DIR"
git fetch -q origin experiments || true
git checkout -q "$COMMIT"
git submodule update -q --init --recursive

cp "$HERE/lbounds.cc" cluster_editing/application/lbounds.cc
CML=cluster_editing/application/CMakeLists.txt
if ! grep -q "add_executable(lbounds" "$CML"; then
    cat >> "$CML" <<'EOF'

add_executable(lbounds lbounds.cc)
target_link_libraries(lbounds ${Boost_LIBRARIES})
set_property(TARGET lbounds PROPERTY CXX_STANDARD 17)
set_property(TARGET lbounds PROPERTY CXX_STANDARD_REQUIRED ON)
set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} lbounds PARENT_SCOPE)
EOF
fi

# the forced includes are needed with recent compilers only
FLAGS="-include cstdint -include limits -include string -include memory -include algorithm \
-include functional -include stdexcept -include optional -include numeric -include cassert"
EXTRA=()
# Boost (program_options) from a conda environment, if it has one
[ -n "${CONDA_PREFIX:-}" ] && EXTRA=(-DCMAKE_PREFIX_PATH="$CONDA_PREFIX")
mkdir -p build_lb && cd build_lb
cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS="$FLAGS" "${EXTRA[@]}" > cmake.log 2>&1 || {
    tail -n 30 cmake.log >&2
    if grep -qi boost cmake.log; then
        echo "Boost not found: sudo apt install -y libboost-program-options-dev" >&2
        echo "  (or, in the conda environment: conda install -y -c conda-forge boost-cpp)" >&2
    fi
    echo "cmake failed, see $PWD/cmake.log" >&2; exit 1; }
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)" lbounds > make.log 2>&1 || { tail -n 30 make.log >&2; echo "build failed, see $PWD/make.log" >&2; exit 1; }
BIN="$(find "$PWD" -type f -name lbounds -perm -u+x | head -n1)"
[ -n "$BIN" ] || { echo "build failed, see $PWD/make.log" >&2; exit 1; }
echo "$BIN"
