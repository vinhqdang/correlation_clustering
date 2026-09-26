# SCC baseline

SCC (Hausberger, Faraj, Schulz, ALENEX 2025), commit 05d0849 of
https://github.com/ScalableCorrelationClustering/ScalableCorrelationClustering
(MIT licence).  The only change, `write_clustering.patch`, makes
`scc_evolutionary` write its final clustering (the release prints only its
cost), so that the cost can be recomputed independently.

```bash
git clone https://github.com/ScalableCorrelationClustering/ScalableCorrelationClustering scc
cd scc && git checkout 05d0849 && git apply ../experiments/scc/write_clustering.patch
./compile_withcmake.sh          # needs MPI (e.g. libopenmpi-dev)
SCC_DIR=$PWD/build python experiments/run_scc.py 60 0 evo pace-exact/exact001.gr ...
```

SCC optimises the objective of sparse signed graphs; `run_scc.py` writes the
complete instance with every non-adjacent pair as an edge of weight -1.
