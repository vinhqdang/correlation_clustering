# KaPoCE baseline

KaPoCE (PACE 2021 winner, GPL-3.0) is used as an external baseline.

```bash
git clone --recursive https://github.com/kittobi1992/cluster_editing.git kapoce
cd kapoce
# heuristic solver (branch used in PACE 2021)
git checkout heuristic_submission
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_CXX_FLAGS="-include cstdint -include limits -include string -include memory -include algorithm -include functional -include stdexcept"
make ClusterEditing
export KAPOCE_BIN=$PWD/ClusterEditing
```

The forced includes are only needed with recent GCC versions.  The heuristic
reads a PACE `.gr` file from stdin and prints the best solution found when it
receives SIGTERM; `experiments/baselines.py` stops it at the time budget.

Root lower bounds of the KaPoCE branch-and-bound (P3 packing and star packing)
are computed with `lbounds.cc`, built against the `experiments` branch:

```bash
git checkout experiments
cp /path/to/this/repo/experiments/kapoce/lbounds.cc cluster_editing/application/
# append an add_executable(lbounds lbounds.cc) block to cluster_editing/application/CMakeLists.txt
mkdir build_exp && cd build_exp && cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS="... -include optional -include numeric -include cassert"
make lbounds
```

This solver stores the instance as a dense matrix, so it is only applicable to
graphs with up to about 10^4 vertices.
