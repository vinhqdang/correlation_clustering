# Correlation clustering: certified algorithms and an evaluation pipeline

This repository contains the code, benchmark pipeline and manuscript for a
study of (min-disagreement) correlation clustering on complete signed graphs,
where the input is the positive graph G+ and every other pair is negative
(equivalently, Cluster Editing).

Main components (`ccbench/`):

| module | contents |
|---|---|
| `graph.py` | CSR graphs, SNAP / PACE readers, components |
| `objective.py` | exact objective evaluation |
| `support.py` | distance-2 support P = E+ ∪ N2, bad-triangle enumeration |
| `pivot.py` | Pivot / KwikCluster |
| `localsearch.py` | vertex-move and multilevel (Louvain-type) local search, weighted objective |
| `insertion.py` | cluster-insertion local search with neighbourhood-cleaning candidates |
| `flip.py` | iterated flipping local search and 3-way pivot |
| `lp.py` | cutting-plane metric LP / exact ILP on the distance-2 support, ACN and CMSY rounding, metric repair |
| `dual.py` | lower bounds: MWU fractional bad-triangle packing, Lagrangian triangle decomposition, greedy packing |

Install and test:

```bash
pip install -r requirements.txt
python -m pytest -q tests
```
