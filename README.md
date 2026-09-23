# Certified correlation clustering at scale

Code, benchmark pipeline and manuscript for a study of min-disagreement
correlation clustering on complete signed graphs (equivalently, Cluster
Editing): the input is the positive graph G+ and every other pair is negative.

The main algorithm, **CertiFlip**, returns a clustering together with a
certified lower bound on the optimum:

* Pivot seed (3-approximation in expectation), then an iterated flipping local
  search with cluster-insertion moves (a practical version of the 2 - 2/13
  scheme of Cohen-Addad et al., STOC 2024); every later step is monotone, so
  the expected cost stays within 3 OPT.
* Lower bound from the **distance-two relaxation** (all LP variables live on
  pairs at G+-distance at most two; farther pairs are provably separated) with
  triangle, star and local subgraph inequalities, computed by **dual
  block-coordinate ascent** with exact block LPs, warm-started by a sparse star
  packing with ruin-and-recreate local search.  The bound is valid at any time.
* **Gap decomposition**: cost(C) - LB splits into non-negative per-pair and
  per-row terms; it certifies regions as locally optimal and guides a
  large-neighbourhood search with exact sub-MIPs.
* When a single block covers the graph and separation converges, CMSY rounding
  of the LP solution is added as a seed, giving a 2.06 guarantee.

## Layout

| path | contents |
|---|---|
| `ccbench/graph.py` | CSR graphs, SNAP / PACE readers, components |
| `ccbench/objective.py` | exact objective |
| `ccbench/support.py` | distance-2 support P = E+ ∪ N2, bad-triangle enumeration |
| `ccbench/pivot.py` | Pivot (KwikCluster) |
| `ccbench/localsearch.py` | vertex moves and multilevel local search (weighted objective) |
| `ccbench/insertion.py` | cluster-insertion local search |
| `ccbench/flip.py` | iterated flipping local search, 3-way pivot |
| `ccbench/lp.py` | cutting-plane LP / ILP on the distance-2 support; ACN and CMSY rounding; star and local subgraph separation |
| `ccbench/dual.py` | packing bounds (greedy, MWU, star packing with local search), MatchFlipPivot |
| `ccbench/blockdual.py` | anytime dual block-coordinate ascent |
| `ccbench/lns.py` | gap map, local certificates, exact sub-MIP neighbourhoods |
| `ccbench/certiflip.py` | the CertiFlip pipeline |
| `experiments/` | datasets, baselines (KaPoCE, Leiden-CPM), benchmark runner, analysis, ablations |
| `results/` | raw CSV results |
| `paper/` | manuscript (LaTeX) |
| `data/pace2021_exact_kapoce_bounds.csv` | published per-instance optima and bounds for the PACE 2021 exact track |

## Installation

```bash
pip install -r requirements.txt
pip install -e .
python -m pytest -q tests
```

KaPoCE (external baseline) is built as described in `experiments/kapoce/README.md`.

## Usage

```python
import ccbench as cc
from ccbench.certiflip import certiflip

g = cc.read_edgelist("graph.txt")          # or cc.read_pace("instance.gr")
res = certiflip(g, time_limit=300, rng=0)
print(res.cost, res.lower_bound, res.certified_ratio)
```

## Reproducing the experiments

```bash
# instances: SNAP graphs are downloaded on first use; PACE 2021 instances from
# https://github.com/PACE-challenge/Cluster-Editing-PACE-2021-instances into data/raw/pace/{exact,heur}
python experiments/run_bench.py --suite pace-exact --algos all --budget 60 --out results/pace_exact.csv
python experiments/run_bench.py --suite snap --algos all --budget 600 --workers 3 --out results/snap.csv
python experiments/ablation.py bounds --budget 300 --out results/ablation_lb.csv
python experiments/ablation.py scaling --out results/scaling.csv
python experiments/analyze.py --pace-exact results/pace_exact.csv --snap results/snap.csv
cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## License

BSD 3-Clause (see `LICENSE`).  KaPoCE is GPL-3.0 and is not redistributed here.
