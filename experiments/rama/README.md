# RAMA baseline

RAMA (Abbas and Swoboda, CVPR 2022) computes a dual lower bound for multicut.
On the multicut instance on the distance-two support P (weight +1 on edges,
-1 on pairs at distance two), its bound plus the number of such pairs is a
lower bound for correlation clustering (floating point, not certified).

```bash
bash experiments/rama/run_rama.sh              # GPU (CUDA >= 11.2, nvcc on PATH)
RAMA_CPU=1 bash experiments/rama/run_rama.sh   # CPU solver, same bound, much slower
```

Results: `results/rama_snap_user.csv`.  RAMA is pinned to commit ee459d6.
