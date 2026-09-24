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
| `ccbench/reduce.py` | exact critical-clique (twin) contraction to a weighted instance |
| `ccbench/anneal.py` | simulated annealing on (weighted) vertex moves and swaps |
| `ccbench/memetic.py` | partition crossover, PX-annealing, memetic search (`pxmem`) |
| `lean/CCProofs/` | Lean 4 proofs of the move/swap formulas, twin lemma, contraction identity and partition crossover |
| `experiments/` | datasets, baselines (KaPoCE, Leiden-CPM), benchmark runner, analysis, ablations, Colab fleet |
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
python experiments/anytime.py results/anytime_hepth.csv 600 ca-HepTh
python experiments/anytime.py results/anytime_dblp.csv 600 com-DBLP
python experiments/plot_anytime.py
python experiments/analyze.py --pace-exact results/pace_exact.csv --snap results/snap.csv \
    --ablation results/ablation_lb.csv --scaling results/scaling.csv
cd paper && pdflatex certified_correlation_clustering && bibtex certified_correlation_clustering && pdflatex certified_correlation_clustering && pdflatex certified_correlation_clustering
```

Lower-bound certificates and their independent check:

```bash
# PACE exact track: CertiFlip writes one certificate per instance, then check all of them
CERT_DIR=results/certificates/pace-exact python experiments/run_bench.py --suite pace-exact \
    --algos certiflip --budget 60 --out results/pace_exact_cert.csv
python experiments/check_certificate.py --dir results/certificates/pace-exact results/pace_exact_cert_check.csv
# a single graph (SNAP certificates, Colab tag c1): run, write and check
python experiments/colab_run_cert.py 600 0 c1 ca-GrQc
```

Head-to-head against KaPoCE (one machine, solvers run one after the other,
both pinned to CPU 0 with `taskset`; KaPoCE seeded through `KAPOCE_SEED`).
Outside Colab set `CC_ROOT` (this repository), `KAPOCE_SRC` (a KaPoCE checkout,
patched and built on first use), `RESULTS_DIR` and `SEQ_READY`:

```bash
python experiments/colab_run_seq.py 600 0 s1 ca-AstroPh      # one graph, one seed
python experiments/seq_stats.py s1 results/headtohead_s1.csv results/headtohead_s1_runs.csv
python experiments/seq_stats.py --pace s1h                     # PACE heuristic track, pooled
python experiments/make_pxmem_tables.py
```

`experiments/colab_fleet.py` distributes these jobs over Colab machines
(`experiments/colab_fleet.sh` keeps it running); results land in `results/colab/`.

## Main results (see `paper/certified_correlation_clustering.pdf`)

* PACE 2021 exact track (200 instances): the CertiFlip bound proves optimality
  on 111 of the 173 instances with known optimum (root bounds of the KaPoCE
  branch-and-bound: 79); published lower bounds improved on the open instances
  exact179 (632) and exact180 (1068); a solution of cost 2788 for exact183
  (published upper bound 2789), stored in `results/solutions/`.  All 200
  certificates are accepted by the independent checker.
* SNAP graphs with up to 10^6 edges: best known solutions certified within
  0.25%-10.1% of optimal (checked certificates), where greedy / fractional
  triangle packings certify only factors 1.17-1.81.
* PXMem against KaPoCE, 600 s each on the same core, fifteen SNAP graphs, ten
  seeds: 67 wins, 49 ties, 34 losses; better on four graphs and worse on one
  by the sign test (one each after Holm correction).
* IteratedFlip stays within 0.10%-0.44% of the PACE 2021 winner KaPoCE on
  SNAP graphs, at a fraction of its running time.

## License

BSD 3-Clause (see `LICENSE`).  KaPoCE is GPL-3.0 and is not redistributed here.
