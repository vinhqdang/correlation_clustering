# Wave 2 (agent "d"): computer search for structure in CEaMP  -- 2026-10-19

Label: VERIFIED-BY-COMPUTER = true for exactly the search space stated, nothing more. No proofs here.

## 0. Tools and reproducibility
* `code/w2d_common.py`: exact solver. Variables e_p = edit indicator for the C(n,2) pairs, all 3*C(n,3)
  transitivity constraints, objective = number of edits (HiGHS via `scipy.optimize.milp`). Extra rows allow
  "cost <= OPT" with side objectives or side constraints, so we can ask lexicographic questions exactly
  ("among optima, minimise the number of edits of type X", "is pair p in SOME optimum?"). Also exact maximum
  P3 packing by set-packing ILP, and the graph families. n=12 solves take about 0.02 s.
  Validation: 60 random graphs n=4..7 agree with the brute-force `opt` of `fpt_tests.py` (0 mismatches).
* Scripts: `w2d_exp1.py` (l, C1, endpoints), `w2d_summ1.py` (summary), `w2d_c1b.py` (targeted C1, large R),
  `w2d_enum.py` (enumerate all optima), `w2d_profile.py` (profile experiment). All seeds are strings built
  from (family, n, index), so reruns reproduce. Results/instances: `code/instances_w2d/` (JSON; includes
  `exp1_results_9_13_40.json`, a 3-record-per-graph table, and a stray `exp1_results_9_9_3.json` from a smoke test).
* Families (n = number of vertices): gnp (p in .2,.35,.5,.65), planted clusters + flip noise (.05,.1,.2),
  twin-class blow-ups of a random graph on 3..n/2 vertices (classes cliques or independent sets),
  near-cliques (clique on n/2..n vertices with edges removed w.p. .1-.3, plus .2 attachment), random bipartite
  (p .4,.7,.9), paths of cliques (block sizes 1..4, partial joins p .15,.4,.8, sometimes closed to a ring).
  C1b adds: core + attached cliques, 2-3 overlapping cliques + 5% noise, noisy complete multipartite,
  hubs with attached cliques. "rand" packings are greedy random maximal; "max" is a maximum packing.
  Packing H is always modification-disjoint (pair-disjoint) and maximal; l = OPT - |H|, c = max P3s through a vertex.

## 1. How large is l = OPT - |H|?   VERIFIED-BY-COMPUTER
Space: 5 sizes n=9..13 x 6 families x 40 graphs = 1200 graphs; per graph 2 random maximal packings and one
maximum packing (3600 records).  Mean l (and P(l<=1)) by family:

| family | random maximal H: mean l (max) | P(l<=1) | maximum H: mean l (max) | P(l<=1) |
|---|---|---|---|---|
| pathcliques | 0.5 (3) | .91 | 0.06 (1) | 1.00 |
| planted | 0.9 (5) | .76 | 0.22 (3) | .97 |
| nearclique | 1.2 (6) | .68 | 0.28 (3) | .94 |
| blowup | 1.9 (13) | .58 | 1.23 (11) | .73 |
| gnp | 2.6 (9) | .31 | 0.98 (4) | .74 |
| bipartite | 5.9 (16) | .12 | 5.25 (15) | .14 |

* Mean l grows with n for random graphs (random maximal H: 1.3, 1.7, 2.4, 2.5, 3.1 for n=9..13).
* Maximum packings shrink l a lot in gnp/planted/nearclique/pathcliques (l=0 in 39% / 81% / 80% / 94%)
  but not in bipartite-like and twin-blow-up graphs: there OPT stays far above the packing size
  (dense bipartite: the packing number is limited because each P3 uses 3 pairs while OPT repairs many).
* Reading for FPT: l is small exactly in the "cluster-like" regimes (planted, near-clique, path of cliques)
  and with good (maximum) packings; l is large for dense bipartite-like inputs, where a parameter
  "l above packing" is not small. This is only a statement about n <= 13.

## 2. Conjecture C1 (some optimum has no unpacked edit with both endpoints outside V(H))
Exact test per instance: minimise the number of unpacked R-R edits subject to cost <= OPT (R = V \ V(H)).
* VERIFIED-BY-COMPUTER, 0 counterexamples among 469 (graph, maximal-H) pairs with l>0 and |R|>=2, drawn from
  the Section 1 space (n=9..13, 6 families, rand1/rand2/max packings; 15 (packing,n) cells of 15-45 instances).
* VERIFIED-BY-COMPUTER, 0 counterexamples among 723 further pairs (l>0, |R|>=2) in the structured large-R
  families of C1b: n=9..13, 4 families x 60 graphs x up to 3 random maximal packings; 331 have |R|>=3, 41 have |R|>=5.
  Not covered: |R|>=6 is rare (only dozens of cases), n>13, adversarial constructions, non-maximal H.
  So C1 survives, but with R mostly 2..4 the evidence is modest; it does not test long twin classes.
* The weaker "forced R-touching edits" are common: in 54 of 816 (l>0, |R|>=1) pairs (6.6%) EVERY optimum has an
  unpacked edit with exactly one endpoint in R (max 2 such edits forced). By family: gnp 18, planted 11,
  nearclique 10, bipartite 7, pathcliques 7, blowup 1; by n: 12 (n=9), 14 (10), 5 (11), 8 (12), 15 (13).
  Instances with n<=10 saved as `instances_w2d/forcedRT_*.json` (the C4 example of fpt_attempt is one pattern).
* Enumeration of ALL optima (`w2d_enum.py`; n=9,10, 5 families x 12 graphs, one random maximal H, l>0, V(H)!=V,
  complete enumeration only; 25 instances, 125 optimal solutions): 24 of the 125 optima (19%) contain an unpacked
  edit with an endpoint in R; in 7 of 25 instances some but not all optima do, in 1 instance all do, in 18 none;
  0 of the 125 contain an R-R edit. Small sample, shows only that the R-touching edit is typically one choice among
  several optima, not a rule.

## 3. Is the best location of an excess edit decided by a bounded-radius local profile?
Definition used. Incidence graph: vertices of G, G-edges, plus one node per packed P3 joined to its a,b,c by
labelled edges. Profile(p, r) of an unpacked pair p = uv: the subgraph within distance r of {u,v}, u,v marked.
Experiment: inside one instance (G,H), l>0, find unpacked pairs p, q with isomorphic Profile(., r)
(exact isomorphism test) such that p lies in SOME optimum and q in NO optimum. Such a pair shows no rule of radius r
(defined on rooted profiles) can say where excess edits go.
VERIFIED-BY-COMPUTER results:
* Space A: n=6..8, 9 families x 6 seeds: no hits.  n=9..11, 9 families x 8 seeds: 2 hits (both paths), radius 1.
* Space B (chains): n=11..15 x {path, cycle, caterpillar, pathcliques} x 14 seeds (one random maximal H or a
  maximum H per seed): 22 hits, all in path P_n; none in cycle/caterpillar/pathcliques at these sizes.
  Minimal examples (all instances in `instances_w2d/profile_hits_11_15_14_4.json`):
  - radius 1: P_11 with H = {012, 567, 8 9 10, 345} (vertex numbers 0..10 along the path), pair 24 in some optimum, pair 79 in none (OPT = 5, l = 1).
  - radius 2: n=12 path, l=1.
  - radius 3: n=13 path, H = {(1,2,3),(4,5,6),(7,8,9),(10,11,12)}, OPT = 6, l = 2:
    pairs (6,8) (in some optimum) and (3,5) (in none) have isomorphic radius-3 profiles; the graph has diameter 12.
  No radius-4 hit up to n=15 (14 seeds each); this reflects the search size, since a radius-r profile
  needs about 3r+ vertices on each side of the two pairs.
* Meaning: for H with an uncovered vertex at the end of the path and a parity defect, the position of the
  optimum edit is decided by the distance to the defect: a DP over bounded-radius profiles alone is lossy, so the
  DP has to carry global state, e.g. a counter along I(H). In these examples the global state is O(1) (chain
  parity), which is consistent with the width-w DP of Theorem A; it does NOT show that profiles plus a
  pointer cannot work, only that radius-r local rules cannot.
* Not shown: any such example in a non-path family, or one with a Lemma-2 viability filter; hits in dense random
  graphs would need larger n because their rooted profiles are almost all distinct (rigid) at n<=11.

## 4. Number of distinct optimal pointer endpoints versus (l, c)
Per (graph, H) with l>0 and |R|>=1 (816 pairs, same space as Section 1): E* = set of vertices that are an
endpoint of an unpacked edit in SOME optimum (exact: one ILP per vertex). Also the number of orbits of E* under
Aut(G,H) (orbit test by isomorphism of the incidence graph with the vertex marked).
VERIFIED-BY-COMPUTER, n=9..13:
* l=1: |E*| <= 6 for every n (max 6 at n=9,10,11,12; 4 at n=13) and every c in 1..5 (max by c: 5,6,6,5,2).
  So at l=1 the endpoint set does not grow with n over this range.
* l=2: max 9; l=3: 11; l=4..7: 12; larger l hits the n limit (13). Ratio |E*|/(2l) reaches 3.0 at l=1, so E* is
  not bounded by the 2l endpoints of a single optimum; the number of optima matters.
* Orbits equal |E*| in all but a few cells (random graphs are rigid), so symmetry does not shrink E* here;
  the twin-class reduction of fpt_attempt Lemma 4(iii) was not exercised by these families.
* Dependence on c: no monotone trend visible (l=2: max 9,8,8,8,6 for c=1..5; l=3: 8,11,10,11,8).
* Conclusion permitted: a bound of the form |E*| <= f(l) is consistent with n<=13 for l=1 (plateau 6) but
  l>=3 already saturates at n, so nothing can be said for l>=3.  A dependence on c is not detected.
  (Plateau 6 at l=1: a single P3-like defect has <=3 vertices, and the cost-neutral shifts move the edit
  to near-twins; the data shows up to 6.)

## 5. What this changes for the paper
1. C1 is a reasonable working assumption for the structure theory (0 / 1192 tested), but still unproved;
   its weaker R-touching version is false as a forcing statement (6.6% of instances force an R-touching edit).
2. A pure "local profile" DP is lossy already for a path (radius 3 shown, radius r expected for 3r+ vertices);
   state must include a global parity/counter component. Handling this with width-w state is exactly the
   Theorem A structure; the open point remains the identity matching between the two ends of a pointer.
3. Regime map for l: small for planted/near-clique/path-of-cliques inputs and for maximum packings, large for
   bipartite-like inputs; if experiments on real data are planned use maximum (not greedy) packings.
4. E*-size evidence is too weak to support a bound g(l,c); a cleaner test needs n >= 20 (needs a faster
   solver than all-optima ILP, or symmetry-reduced search) and l fixed at 2 or 3.
Suggested next steps: (a) extend radius search to long path-like I(H) families with n up to 20 (solver is fast
enough: one ILP is 0.05 s at n=20, per-pair tests are the cost); (b) run C1 on structured instances with
|R| >= 6, e.g. one long twin clique plus a small core, n up to 16.
