# Certified lower bounds for correlation clustering on large sparse graphs

This repository holds the code, the raw results, the archived certificates
and the manuscript (`paper_mpc/`) of a study of lower bounds for correlation
clustering on complete graphs, also called cluster editing. The input is the
positive graph G+; every other pair is negative.

Heuristics find good clusterings of large graphs, but on such graphs nobody can
say how far a clustering is from optimal. This code computes lower bounds on
sparse graphs with up to 10^6 vertices and exports them as **certificates**.
An independent checker verifies each certificate in exact arithmetic against
the raw instance file.

## What is here

* **Distance-two support.** Optimal clusters have diameter at most two, so
  every pair at distance three or more is separated in every optimal
  clustering. The metric LP, its star and local subgraph inequalities, and
  every packing of induced stars can therefore be restricted to the pairs at
  distance at most two (`ccbench/support.py`, `ccbench/lp.py`).
* **Star local search** (`experiments/sstar/sstar.cc`). A sparse, independent
  reimplementation of the star-packing heuristic of the KaPoCE
  branch-and-bound (Bläsius et al., SEA 2022), which extends the P3-packing
  search of Gottesbüren et al. (SEA 2020).
  * Its memory is linear in the graph plus the packing.
  * It writes its stars directly as certificate rows.
  * It is not a better heuristic than the original; it runs where the dense
    original does not fit into memory.
* **Anytime block-dual bound** (`ccbench/blockdual.py`). Dual
  block-coordinate ascent with exact block LPs (HiGHS), started from any
  packing. The bound is valid after every step.
* **Independent checker** (`experiments/check_certificate.py`), about 600
  lines, sharing no code with the solver. It:
  * parses the raw file itself;
  * binds each certificate to the file by SHA-256 hashes;
  * checks the validity of every row;
  * evaluates the bound in integer arithmetic.

  The format is specified in
  [`docs/certificate_format.md`](docs/certificate_format.md). The checker
  also accepts the star packings of the KaPoCE branch-and-bound, which the
  13-line hook in `experiments/kapoce/` exports.
* **Primal methods** supply the clusterings that the bounds certify:
  * CertiFlip (`ccbench/certiflip.py`): Pivot, then iterated flipping local
    search, then block-dual bound, then a gap-guided LNS with exact sub-MIPs;
  * PXMem (`ccbench/memetic.py`): a memetic search with exact partition
    crossover.
* **Lean 4 development** (`lean/CCProofs/`). It proves the inequality that
  the checker evaluates, the separation of far pairs by optimal
  clusterings, the move, swap and contraction formulas, and partition
  crossover.

## Main results (manuscript `paper_mpc/main.pdf`)

* **PACE 2021 exact track, 173 instances with known optimum.**
  * Under the stopping rule of the KaPoCE branch-and-bound, the star local
    search proves optimality on 83 of them in a median of 0.5 s. The
    published root bound proves it on 79.
  * Run for 60 s, the same heuristic proves optimality on 122. Together with
    the CertiFlip bound, whose block LPs close 7 instances that five times the
    packing time or four more seeds do not, 130 instances are proved optimal.
  * The checked bounds improve the published root bounds of most of the 27
    open instances.
* **27 SNAP graphs, up to 1.1·10^6 vertices.** Every graph has a checked
  lower bound. On the 23 graphs with an archived clustering, the certified
  gap ranges from 0.2% to about 21%. Triangle packings leave gaps of 17% to
  90%.
* **Certificates.** Every bound in the paper comes from an archived
  certificate that the checker accepts, including the packings of the KaPoCE
  branch-and-bound.

The exact numbers, their spread over seeds and the controls are in Section 7
of the manuscript.

## Layout

| path | contents |
|---|---|
| `ccbench/` | Python package (Numba kernels): graphs and readers, objective, support, Pivot, local searches, flipping, LP and ILP on the support, packing bounds, block-dual ascent, gap map and LNS, CertiFlip, twin contraction, annealing, PXMem |
| `experiments/sstar/` | star local search (C++17, single file) |
| `experiments/check_certificate.py` | independent certificate checker |
| `experiments/colab_run_*.py` | one run of one method on one instance, with checked certificate (the scripts behind every result file) |
| `experiments/colab_fleet.py`, `colab_fleet.sh`, `colab_adopt.py`, `colab_reap.py` | job runner on a fleet of Google Colab machines |
| `experiments/mpc_tables.py`, `make_figures.py` | all tables, numbers (`paper_mpc/numbers.tex`) and figures of the manuscript, from `results/` |
| `experiments/kapoce/`, `rama/`, `scc/` | build and run scripts for the external baselines (not redistributed) |
| `experiments/run_bench.py`, `ablation.py`, `kapoce_root_snap.py`, `run_local_*.py`, `run_local_kroot.sh`, `run_all_local.sh` | benchmark runner and local jobs (triangle packings, KaPoCE root bounds, PACE bound runs) |
| `experiments/recheck_all.sh`, `restart_background.sh` | re-check all certificates; restart long-running jobs |
| `experiments/heldout.md` | held-out graphs and analysis plan of the head-to-head comparison, fixed before its runs |
| `results/colab/`, `results/local/` | one JSON file per run |
| `results/certificates/colab/`, `results/local/certs/` | the certificates and clusterings behind every reported value |
| `results/mpc/` | complete re-check of all certificates |
| `results/*.csv`, `results/scc/` | earlier benchmark results used for the triangle packings, SCC and the published PACE bounds |
| `data/` | instance manifest (SHA-256), instance table, published PACE 2021 bounds |
| `docs/` | certificate format, literature notes |
| `lean/CCProofs/` | Lean 4 development |
| `tests/` | unit tests (`python -m pytest -q tests`) |
| `paper_mpc/` | manuscript (Springer Nature template) |

## Installation

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-lock.txt && pip install -e .
python -m pytest -q tests
python experiments/fetch_data.py      # SNAP and PACE 2021 instances, checked against data/MANIFEST.sha256
```

The star local search is compiled on first use with `g++ -O2 -std=c++17`.

## Usage

```python
import ccbench as cc
from ccbench.certiflip import certiflip

g = cc.read_edgelist("graph.txt")          # or cc.read_pace("instance.gr")
res = certiflip(g, time_limit=300, rng=0, cert_path="cert.npz")
print(res.cost, res.lower_bound, res.certified_ratio)
```

Star local search and its certificate on one instance, followed by the check:

```bash
RESULTS_DIR=out CC_ROOT=$PWD python experiments/colab_run_lp.py 600 0 lsstar ca-GrQc
python experiments/check_certificate.py data/raw/ca-GrQc.txt.gz out/certs/lsstar_ca-GrQc_0.npz
```

Re-check every archived certificate (no solver needed):

```bash
bash experiments/recheck_all.sh
```

[`REPRODUCE.md`](REPRODUCE.md) maps every table and figure of the manuscript
to the command that produced it.

## License

BSD 3-Clause (see `LICENSE`). KaPoCE (GPL-3.0), RAMA and SCC are external
baselines and are not redistributed here. `experiments/kapoce/star_dump.inc`
is a 13-line hook that is compiled into the KaPoCE sources, and the resulting
program is covered by KaPoCE's licence.
