# Experiments: exact CEaMP solver and P3-packing interaction on real graphs

All code is in `code/exp_*.py` and is deterministic (fixed seeds, listed below).
Raw tables: `code/results_pace.csv` (173 PACE 2021 exact instances solved in 1 h),
`code/results_snap.csv` (7 SNAP graphs), `code/results_dp.csv` (l = 0 DP on PACE),
`code/results_validation.txt`. Logs `code/run_*.log`.
Reproduce: `python3 exp_solver.py 20000; python3 exp_packings.py pace; python3 exp_packings.py snap; python3 exp_dp.py`
(about 11 + 3 + 25 minutes on 4 cores; numba and networkx 3.6.1 are installed; networkx was not needed, see below).

Status labels: [computed] = measured here; [proved] = follows from analysis_01; [heuristic] = no guarantee.

## 1. Solvers (task 1)

* `exp_solver.brute_opt`: optimum cluster editing cost by enumerating all set partitions
  (restricted growth strings, vectorised with numpy; n <= 9). `ceamp_brute(A,H,l)` decides OPT <= |H| + l.
* `exp_solver.dp_l0`: the l = 0 CSP of analysis_01 (Obs. 1, 2). Variables x_i in {0,1,2} (which pair of
  packed P3 P_i is edited). One constraint per vertex triple that contains a packed pair and that
  could be an induced P3 after edits; its scope is the set of owners of its packed pairs (size <= 3).
  Triples with no packed pair are never constrained, so the code also checks that G has no induced P3 with
  all three pairs outside the packing (`unpacked_p3`, else the l = 0 answer is NO).
  Solved by bucket elimination (feasibility tables, boolean arrays of size 3^(bag)) along an elimination
  ordering of the constraint graph. An elimination ordering is a tree decomposition with bags
  {v} + later neighbours, so this is the dynamic program over a tree decomposition with domain 3.
  I wrote the elimination heuristics (min-degree and min-fill-in) myself, in `exp_common.py`, instead of
  `networkx.approximation.treewidth_*`, because I needed an abort at a width cap and the ordering itself.
  I did not cross-check them against networkx (not done; the validation below tests only the answers).
* Validation [computed]: 20000 random instances (n = 4..9; G(n,p) with p in {.2,.35,.5,.7}, and planted
  cluster graphs with 1-4 random edge flips; random packings, 70% maximal and 30% with P3s randomly dropped),
  seed 12345. 18749 valid instances (a maximal packing needs >= 1 P3): 14132 YES, 4617 NO, 12841 maximal and
  5908 non-maximal packings. **0 mismatches** between `dp_l0` and "OPT == |H|" from brute force.
  (The l = 1, 2 comparison in the same script only checks `ceamp_brute` against OPT - |H|; it is a
  consistency test of the wrapper, not an independent test.) Planted instances make YES dominant (75%);
  all NO answers are also correct. The DP is only validated for n <= 9, which means constraint scopes of size 2
  and 3 occur, but I did not log how often bags of width > 5 occurred.

## 2. Packings (task 2): what I did and what it is not

Packings are **pure P3 packings** (induced P3 a-b-c, edges ab, bc, non-edge ac; no pair is used twice).
I did NOT use the star packings of `experiments/sstar/sstar.cc`: a star with k leaves uses its k edges and
C(k,2) leaf non-edges, but only floor(k/2) pairwise pair-disjoint P3s can be cut from it using the star's
pairs. So a star packing of value k - 1 turns into a P3 packing of value about k/2, and the P3 packings here
are weaker lower bounds than the star bounds (consistent with `low_p3 <= low_star` in the bounds csv).
* Greedy (seed 0): centres in random order; for each centre scan neighbour pairs in random order; take
  every induced P3 whose three pairs are unused. One pass gives a **maximal** packing (every induced P3 of G
  shares a pair with the packing; `is_maximal` verified on all 173, True).
* Local search (seed 1): repeat `min(300000, 100|H| + 20000)` times (SNAP: up to 1.5M): free a random packed
  P3, list all induced P3s that are then entirely free (they all use a freed pair), re-pack greedily in random
  order (best of up to 4 orders). |H| never decreases; a result with >= 2 P3s is an improvement, 1 P3 a plateau
  move. Always maximal, validated with `verify_packing` (modification-disjointness) on all 173.
  This is a (1,k)-swap plateau search, a heuristic with no optimality guarantee. Max 9 s per PACE instance.
* c-capped greedy (cap = 1, 2, 3): same greedy but a P3 is rejected if one of its vertices already lies in
  `cap` packed P3s. These are exactly 1/c-integral packings in the sense of Li-Pilipczuk-Sorge; they need not be maximal.
* I(H): nodes = packed P3s, adjacent iff they share a vertex. c = max number of packed P3s on a vertex
  (so I(H) has a clique of size c and **tw(I(H)) >= c - 1**). `tw_lb` = max(c-1, degeneracy) [proved lower bound].
  `tw_ub` = best of min-degree and min-fill-in elimination orderings, aborted at width 100 [heuristic upper bound];
  an empty entry means the heuristic exceeded 100 (or I(H) had > 4M edges), not that tw is unknown-small.

## 3. Results on the 173 PACE 2021 exact instances (n 10..620, density 0.01..1.0)

|H|/OPT is a valid lower-bound ratio (OPT known from the csv, `opt_or_empty`). [computed]

| quantity | greedy | local search |
|---|---|---|
| mean / median / min of &#124;H&#124;/OPT | 0.911 / 0.926 / 0.714 | 0.976 / 0.999 / 0.813 |
| instances with &#124;H&#124; = OPT | 20 | **77** |
| instances with ratio >= 0.99 | | 111 |
| instances with ratio >= 0.90 | | 158 |

Local-search |H| compared with the csv's `low_p3` (the KaPoCE-style P3 bound): >= on 161 of 173 (> on 80, < on 12;
max ratio 1.029, min 0.992). So these packings are as good as the ones already in the repository.

Interaction parameter of the local-search packing (same instances):

| statistic | min | 25% | median | 75% | max |
|---|---|---|---|---|---|
| c (max P3 per vertex) | 2 | 19 | 38 | 74 | 307 |
| max degree of I(H) | 2 | 44 | 90 | 193 | 892 |
| tw(I(H)) lower bound | 2 | 22 | 49 | 93 | 306 |
| tw(I(H)) upper bound (48 instances only) | 2 | 22 | 49 | 70 | 95 |
| components of I(H) | 1 | 1 | 1 | 3 | 31 |

* **The interaction parameter is not small on real graphs.** tw(I(H)) >= 15 is proved (c - 1 >= 15) for 159 of 173
  instances; only 14 have c <= 10. The heuristic upper bound exceeds 100 (the cap) on 125 instances.
  c is about half the maximum degree of G (median c / maxdeg(G) = 0.53, max 0.79), which is expected: a vertex of
  degree d is the centre or end of up to about d/2 pair-disjoint P3s.
* Components: I(H) is connected on 117 of 173; 56 have more than one component, typically sparse graphs
  (median 7 components for density <= 0.05, 1 for density >= 0.2).
* Relation to density (Spearman, 173 instances): c vs density 0.44, vs average degree of G 0.85, vs maxdeg(G) 0.91,
  vs n 0.47. So c and tw lower bound track the **degree** of G, not density as such. By density bin (medians):

| density | #inst | n | c | max deg I(H) | tw lb | components | ratio LS | ratio greedy | #(OPT=&#124;H&#124;) | tw ub found (<=100) |
|---|---|---|---|---|---|---|---|---|---|---|
| (0,.05] | 19 | 250 | 16 | 37 | 18 | 7 | 0.985 | 0.926 | 5 | 15 |
| (.05,.1] | 16 | 212 | 22 | 49 | 25 | 4.5 | 0.995 | 0.928 | 5 | 7 |
| (.1,.2] | 19 | 216 | 30 | 78 | 40 | 3 | 0.999 | 0.954 | 9 | 5 |
| (.2,.4] | 50 | 125 | 40 | 104 | 52 | 1 | 0.992 | 0.866 | 14 | 11 |
| (.4,.6] | 37 | 129 | 46 | 118 | 64 | 1 | 1.000 | 0.925 | 19 | 6 |
| (.6,.8] | 22 | 190 | 87 | 218 | 98 | 1 | 1.000 | 0.962 | 15 | 1 |
| (.8,1] | 10 | 250 | 114 | 235 | 113 | 1 | 1.000 | 1.000 | 10 | 3 |

  Sparse graphs have the smallest c, but even at density <= 0.05 the median c is 16. Dense graphs are
  near-complete, so a packing is easy to make tight (ratio 1) but c is large (about n/2).
  Small tw(I(H)) is rare: 8 instances have tw ub <= 14 (n = 10, 20, 20, 30, 159, 232, 250, 263), all with OPT = |H|.
* Tight packings are the large-c ones, not the small-c ones: among the 77 instances with OPT = |H|, median c is 37.

Trade-off against c (c-capped greedy, ratio = |H_cap| / OPT, median over 173; seed 0):

| cap | median ratio | max c | tw ub found on | max deg I(H) | median #components |
|---|---|---|---|---|---|
| 1 (vertex-disjoint) | 0.048 | 1 | 173 (tw = 0) | 0 | 36 |
| 2 | 0.097 | 2 | 173 (median 7, max 61) | 3 | 1 |
| 3 | 0.141 | 3 | 165 (median 20, max 99; 8 instances above the cap) | 6 | 1 |

Bounded c is cheap for I(H) (subcubic for c = 2) but the bound |H| collapses to 5-14% of OPT, so a
1/c-integral greedy packing is useless as a lower bound on these graphs. Even c = 2 gives tw up to 61 for I(H).
(The packings are not optimised for |H| under the cap; a better optimiser could do somewhat better, I did not try.)
The 'tw ub found on' column counts instances where the heuristic stayed within the cap 100; 'median/max' are over those.

## 4. SNAP graphs (n 986..9875; no OPT known) [computed]

| graph | n | m | maxdeg | &#124;H&#124; greedy | &#124;H&#124; LS | c (LS) | maxdeg I(H) | tw lb | comps | &#124;H&#124; cap 1 / 2 / 3 | c-cap2 tw ub |
|---|---|---|---|---|---|---|---|---|---|---|---|
| email-Eu-core | 986 | 16064 | 345 | 7941 | 8031 | 243 | 446 | 242 | 1 | 237 / 466 / 694 | 52 |
| bitcoinalpha (positive) | 3683 | 12972 | 507 | 6362 | 6474 | 267 | 412 | 266 | 2 | 477 / 832 / 1133 | 12 |
| facebook_combined | 4039 | 88234 | 1045 | 38842 | 40942 | 528 | 944 | 527 | 1 | 1145 / 2263 / 3332 | 47 |
| bitcoinotc (positive) | 5573 | 18591 | 788 | 9097 | 9261 | 433 | 731 | 432 | 3 | 676 / 1171 / 1576 | 24 |
| ca-GrQc | 5241 | 14484 | 81 | 4458 | 4916 | 39 | 89 | 38 | 101 | 855 / 1587 / 2125 | 12 |
| wiki-Vote | 7115 | 100762 | 1065 | 50280 | 50371 | 864 | 1491 | 863 | 4 | 1047 / 2010 / 2920 | > 100 |
| ca-HepTh | 9875 | 25973 | 65 | 10659 | 11448 | 39 | 101 | 38 | 93 | 1795 / 3326 / 4501 | 78 |

(n counts non-isolated vertices of the edge list. The file is read as an undirected simple graph; signed graphs
keep only positive edges.) Same picture as on PACE: c is of the order of the maximum degree / 2 on hub graphs (c = 243..864) and only
reaches tens (39) on the collaboration networks with max degree 65..81; a tw upper bound below 100 was found only for c-capped
packings. For these graphs I have no OPT, so no |H|/OPT; I only note that one pass of local search
adds 1% to 10% to the greedy packing (ca-GrQc +10%, ca-HepTh +7%).
(The sstar.cc / star lower bounds of paper 1 for these graphs are not recomputed here.)

## 5. Does the l = 0 DP decide the OPT = |H| instances quickly? (task 3) [computed]

`exp_dp.py`: local-search packing (same seeds), build the CSP (|H| <= 3000), elimination width of the **constraint
graph** (subgraph of I(H); edges only where a non-trivial triple constraint links two P3s), DP if width <= 16.

* Instances 173; CSP built for 142 (the others have |H| > 3000). Constraint-graph width upper bound (cap 150) found for 58.
  DP run on **9**, all answered **correctly** (8 with OPT = |H| answered YES, exact097 with OPT > |H| answered NO).
* Each DP run took at most 0.63 s (CSP build <= 0.06 s). Widths: 2, 3, 4, 6, 12, 5, 3, 13 (YES) and 15 (NO, exact097,
  run with the larger cap of 16).
* Of the 77 instances with OPT = |H|: 8 were decided by the DP, 31 had a constraint-graph width upper bound at all
  (median 50, max 120), the rest exceeded the cap or were too large. 66 of the 77 have tw(I(H)) >= 15 proved (c - 1 >= 15).
  For the other instances the DP with domain 3 needs 3^(w+1) states, so w around 50 is out of reach.
* On the 48 instances where both heuristic bounds exist, the constraint-graph width bound is **identical** to the
  tw(I(H)) bound (difference 0 everywhere): dropping trivial constraints did not shrink the width by this heuristic.
  (Upper bounds only; the true treewidths may differ.)
* Honest reading: l = 0 is decided instantly **whenever the width is small**, but on these instances the width is
  small only for 8 of 77 tight ones, i.e. tw(I(H)) is **not** a parameter that explains tractability of exact solving
  on PACE-sized real graphs. Exact solvers there rely on branch and bound with the same packings as lower bounds.

## 6. Caveats

* tw numbers are heuristic UPPER bounds (min-degree and min-fill-in; min-fill only for |H| <= 3000), and absent when
  they exceeded the cap. The lower bounds c - 1 and degeneracy are rigorous. The 'tw lb' column is a lower bound for tw(I(H)),
  not for the constraint graph on which the DP runs.
* c and tw depend on the packing; I maximise |H| and do not minimise c. A packing with the same |H| and smaller c may exist
  (not searched). The c-cap experiment suggests it is hard on these graphs, but it is a greedy, not an optimiser.
* LS packings reach OPT on 77 of 173 solved instances; the 'solved in 1h' set is biased to easier instances (the 27
  unsolved ones are not in the tables).
* The DP validation is on n <= 9; real-instance consistency (9 instances) is an extra check, not a proof.
* Observation 2 (scope is a clique of I(H)) was exercised implicitly: the DP uses the owners of each triple as scope, and
  the brute-force agreement on 18749 instances would fail if an interaction outside I(H) were missed.
