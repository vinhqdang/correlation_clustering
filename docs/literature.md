# Correlation clustering: literature as of September 2026

Min-disagreement correlation clustering on complete signed graphs
(equivalently, Cluster Editing).  Full bibliographic entries are in
`paper/refs.bib`.

## Approximation ratios on complete graphs

| Year | Ratio | Method | Reference |
|---|---|---|---|
| 2002/2004 | O(1) | combinatorial | Bansal, Blum, Chawla (FOCS 2002; Machine Learning 2004) |
| 2003/2005 | 4; APX-hard | metric LP, region growing | Charikar, Guruswami, Wirth (JCSS 2005) |
| 2005/2008 | 3 (Pivot), 2.5 (LP pivot) | random pivot | Ailon, Charikar, Newman (JACM 2008) |
| 2009 | deterministic 3 / 2.5 | LP-chosen pivots | van Zuylen, Williamson (MOR 2009) |
| 2015 | 2.06 | metric LP, non-linear rounding | Chawla, Makarychev, Schramm, Yaroslavtsev (STOC 2015) |
| 2022 | 1.994 + eps | Sherali-Adams | Cohen-Addad, Lee, Newman (FOCS 2022) |
| 2023 | 1.73 + eps | preclustering + SA | Cohen-Addad, Lee, Li, Newman (FOCS 2023) |
| 2024 | 2.4 in polylog MPC rounds | parallel LP/pivot | Cao, Huang, Su (SODA 2024) |
| 2024 | 1.485 + eps (STOC version claimed 1.437; corrected in arXiv v3) | cluster LP | Cao, Cohen-Addad, Lee, Li, Newman, Vogl (STOC 2024) |
| 2024 | 2 - 2/13 + eps, combinatorial, O~(n) | local search with flips | Cohen-Addad, Lolck, Pilipczuk, Thorup, Yan, Zhang (STOC 2024) |
| 2025 | 1.485 + eps in O~(2^poly(1/eps) n) | fast cluster LP | Cao et al. (STOC 2025) |
| 2025 | no pivot order beats 3 | lower bound for Pivot | Fischer, Kipouridis, Klausen, Thorup (STACS 2025) |
| 2026 | fully dynamic 1.485 | static-to-dynamic | Cao et al. (ICALP 2026) |
| 2026 | 2.9991 | pivot after removing good clusters | Lolck, Thorup, Yan (arXiv 2603.12052) |
| 2026 | deterministic 1.3865 + eps; cluster-LP gap in [4/3, 1.3865] | approximate dual separation | Garcia-Soriano, Schohn (arXiv 2607.27829, preprint) |

Hardness: NP-hard to approximate within 24/23 - eps unless P = BPP (Cao et
al. 2024).  General weighted graphs: O(log n), equivalent to multicut, no
constant factor under UGC.

## Practical and large-scale work

* Parallel / MPC / streaming / dynamic: ClusterWild!/C4 (NeurIPS 2015),
  constant-round MPC (ICML 2021), almost-3 in O(1) rounds (FOCS 2022),
  single-pass streaming (SODA 2023), dynamic and sketch-based algorithms
  (2024-2026).  Evaluations of these use small graphs or score community
  recovery rather than the objective.
* Heuristics: vote / best-one-element moves (Elsner and Schudy 2009),
  LambdaCC / Louvain / Leiden at resolution 1/2 (Veldt, Gleich, Wirth 2018),
  multilevel memetic signed-graph clustering (Hausberger, Faraj, Schulz,
  ALENEX 2025), PACE 2021 solvers (KaPoCE won both tracks).
* Exact: branch-and-bound with star, P3 and K_{a,b} packing bounds
  (Blasius et al., SEA 2022) solves 173 of the 200 PACE exact instances;
  cutting planes for clique partitioning (Grotschel and Wakabayashi 1989).
* Lower bounds at scale: matching / open-wedge packings (Veldt, ICML 2022)
  certify only about a factor of two; metric-LP projection solvers reach
  about 10^4 vertices.

## Gaps addressed in this repository

1. No implementation of the sub-3 near-linear algorithms existed; `ccbench.flip`
   implements the iterated flipping local search of the STOC 2024 paper.
2. No certificates on large graphs: the distance-two relaxation, sparse star
   packings and the anytime block dual give certified gaps of 0.2%-9.2% on SNAP
   graphs with up to 10^6 edges.
3. No benchmark of the objective at scale: `experiments/` runs ten solvers and
   five lower bounds on PACE 2021 and SNAP instances.
