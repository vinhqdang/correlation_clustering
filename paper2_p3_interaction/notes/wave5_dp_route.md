# Wave 5 (prover B): the cluster-identity DP route to FPT in (l, w)

Labels: PROVED (proof here), VERIFIED-BY-COMPUTER (finite check, evidence only), CITED, CONJECTURE, FAILS.
Code: code/w5b_lib.py (library), w5b_tests.py (all set partitions, n<=11), w5b_ilp.py (HiGHS, n=11..13),
w5b_cand.py (partner locality), w5b_family.py (far-default-link family). Own code; nothing imported from other agents.

## 0. Verdict
* Route does NOT give FPT in (l, w). I have a precise statement of why and of the missing lemma (Sec. 6).
* New PROVED structure (pair-maximal H): (a) the requested cost decomposition; (b) endpoints of edited hard pairs have hard
  degree <= D+l; (c) every cluster touching an edited hard pair has <= 2D+2l+1 vertices, all of hard degree <= 2D+3l, and is
  connected in G; every other cluster is a union of whole hard components (or exactly one); (d) given ONE endpoint u of an
  edited pair, the partner lies in a set of size g(l,D) (deleted: N_h(u); inserted: a G-ball of radius 2D+2l).
  So the tiny-tiny pointer is NOT an identity problem of size n^2; only the POSITION of ~<=l bounded objects is free.
* Consequences: n^{l} * g(l,D) * 7^{w+l+1} (halves the exponent of Theorem A); FPT in (l, D, omega) where omega = treewidth
  of "soft pairs + G-edges" (Thm 8); FPT-type n*g for solutions whose edited pairs sit in one ball.
* FAILS: omega is not bounded by any f(l,w): a YES instance with l=0, w=0, c=1 has tw >= tw(cubic expander) (Sec. 5).
  Free ("default") far hard links, not budgeted pairs, break the connected-subtree/bounded-width property.
* Inserted pairs never create a far merge on their own (every cluster of an optimum is G-connected, Lemma 4), so the budget
  question "does l bound far merges" has the answer: number yes (<= l components), identities no (position unbounded, n^r).

## 1. Setting
H pair-maximal packing (no induced P3 of G has all 3 pairs hard); by MAX (audit4 Sec.1, CORRECT) this is WLOG with
w' <= w+|M|, c' <= c+|M|, D' <= D+2|M|, l' <= l; below H, w, D, l denote the primed values. soft = pair inside a member of H;
hard = other pair; N_h(v) = hard neighbours, d_h = |N_h|; D = max number of soft partners of a vertex (<= 2c <= 2(w+1)).
A solution is a partition pi of V, cost = number of pairs with (same cluster) != (edge). YES: cost <= |H|+l.
U = edited hard pairs; |U| <= cost-|H| <= l (every member of H needs its own edit, audit4 Sec.0). del = hard edge across clusters,
ins = hard non-edge inside a cluster; del_x, ins_x = number at vertex x. A hard component = component of the graph of hard edges
(= proto-cluster of LPS). It is clean if it has no hard non-edge inside (all its non-adjacent pairs are soft).

## 2. Cluster-partition cost decomposition (the requested lemma)
### Lemma 1 (CP)  PROVED, VERIFIED-BY-COMPUTER (413,258 partitions of cost <= |H|+3, n=6..11, 0 violations)
(a) cost(pi) = |E| + sum_{C in pi} f(C), f(C) = C(|C|,2) - 2 e(C), e(C) = number of G-edges inside C.
(b) cost(pi) = sum_i c_i + h(pi), c_i = number of edited pairs among the 3 pairs of P_i (c_i >= 1, depends on pi restricted to V(P_i)
    only), h(pi) = number of edited hard pairs = |U|. Excess = sum_i (c_i - 1) + |U| <= l.
(c) Let (T, beta) be any tree decomposition of the soft graph W (width <= 3w+2, audit4 Sec.5). Give each soft pair and each hard pair
    lying in a common bag to one such bag. Then cost(pi) = sum_t loc_t(pi|beta(t)) + far(pi), where loc_t counts the edited
    pairs assigned to t (depends on the partition of beta(t) and on packed-pair choices only) and far(pi) = (#far hard edges cut) +
    (#far hard non-edges inside clusters) <= |U| <= l. Every far pair is hard (soft pairs lie in bags).
Proof. (a) cost = #non-edges inside + #edges across = sum_C [C(|C|,2) - e(C)] + [|E| - sum_C e(C)]. (b) the pair sets of the P_i are
disjoint; c_i >= 1 because the triple of P_i is an induced P3 unless one of its pairs is edited. (c) partition the pairs by
location. QED
Reading: the local terms are exactly the l=0 CSP of Thm A (fpt_attempt); all non-locality is in far(pi), which is a function of
the far hard pairs only. (a) is the form used by the DP of Sec. 4.

## 3. Structure of hard pairs in a solution
### Lemma 2 (HE, H2)  PROVED, VERIFIED-BY-COMPUTER (46,900 + 46,454 configurations in the main run, 0 failures)
(i) a, b, z distinct, ab and az hard edges => zb in E or {z,b} soft.  (ii) uv a hard non-edge => N_h(u) and N_h(v) are disjoint.
Proof. Otherwise a-b-z (resp. u-z-v) is an induced P3 with all three pairs hard, contradicting pair-maximality. QED
### Lemma 3 (EX)  PROVED, VERIFIED-BY-COMPUTER (1.5M edited-pair checks on every solution of cost <= |H|+3, n<=11; ILP n=11..13)
Let pi be any solution with U the edited hard pairs, uv in U, C = cluster(u).
(a) deleted uv: d_h(u) <= D + del_u + del_v - 1.  inserted uv: d_h(u) <= D + ins_v - 1 + del_u.  In both cases d_h(u) <= D + l.
(b) |C| <= 1 + D + d_h(u) + ins_u <= 2D + 2l + 1.
(c) every y in C has d_h(y) <= |C| - 1 + del_y <= 2D + 3l =: B.
Proof. (a) Deleted: let z be a hard neighbour of u with uz kept, so z in C and v not in C. By Lemma 2(i) (centre u): zv in E or soft.
If zv is a hard edge it is cut (z in C, v not): it is a deleted edge at v other than uv; otherwise z is a soft partner of v (<= D).
So #kept <= D + del_v - 1, and d_h(u) = #kept + del_u. Inserted: u, v in C. A kept hard neighbour z lies in C, so zv is an edge of G xor pi:
a hard edge zv is impossible (z in N_h(u) cap N_h(v) = empty, Lemma 2(ii)); so zv is soft (<= D) or an inserted hard pair
(<= ins_v - 1). (b) w in C-u: uw is an edge of G xor pi, hence a soft pair (<= D), a kept hard edge (<= d_h(u)) or inserted hard (<= ins_u).
(c) hard neighbours of y inside C: <= |C|-1; outside C they are cut. QED
Tightness: slack 0 attained in the finer forms of (a) and (b) (so shaving 1 is falsified); coarse bounds are loose (ILP n<=13: size slack >= 5).
### Lemma 4 (F1)  PROVED (classical), VERIFIED-BY-COMPUTER (all 17,000+ optimal partitions in the runs)
In an optimal partition every cluster C is connected in G[C]. Proof: if C = C1 + C2 with no edge between, splitting saves |C1||C2| > 0
inserted pairs. QED
### Lemma 5 (NF, normal form)  PROVED given Lemmas 2-4, VERIFIED-BY-COMPUTER (beta part: 1.25M clusters)
Let pi be optimal, cost <= |H|+l. Call C exceptional if it contains an endpoint of a pair in U.
(alpha) exceptional C: |C| <= 2D+2l+1, C inside L_B = {v : d_h(v) <= B}, G[C] connected (edges soft or hard-kept, inside L_B).
        At most 2|U| <= 2l exceptional clusters.
(beta)  non-exceptional C: no hard edge leaves C and no hard non-edge lies in C, so C is a union of whole hard components, each clean,
        with all cross pairs soft. If C has >= 2 components then |C| <= 2D (x in K1 has every vertex of C-K1 as soft partner; same for K2).
(gamma) BIG: if d_h(v) > B then v is in a beta-cluster, and cluster(v) = hard component K(v) exactly (> 2D vertices). If K(v) is not
        clean the instance is NO. VERIFIED-BY-COMPUTER: 391 non-trivial cases.
Proof. (alpha) Lemma 3 + Lemma 4. (beta) definitions. (gamma) by Lemma 3(a),(c) v is neither an endpoint nor a member of an exceptional
cluster; K(v) has > B >= 2D vertices so it cannot share a cluster. QED
### Corollary 6 (CAND, candidate partners)  PROVED, VERIFIED-BY-COMPUTER (w5b_cand.py: 4,640 pairs of optimal partitions, 0 failures)
If uv in U (optimal pi): deleted => v in N_h(u), |N_h(u)| <= D+l; inserted => dist_{G[L_B]}(u,v) <= |C|-1 <= 2D+2l (path inside C),
so v is in Ball(u, 2D+2l) of G[L_B], which has <= (B+D+1)^{2D+2l} vertices (degrees in G[L_B] <= B+D). Observed: the distance of an
inserted pair never exceeded 2 (n<=10). CONJECTURE: always <= 2 (not proved).
This answers wave2_hardness Sec.4 / Sec.8 item 1 in the strong form for the Cover case: the locality statement holds relative to ONE endpoint.

## 4. Algorithms
### Theorem 7 (n^l)  PROVED (modulo Thm A of fpt_attempt, MAX, Lemmas 3-5)
CEaMP is solvable in time n^{l+O(1)} * g(l,D) * 7^{w+l+1}, g(l,D) = (l (B+D+1)^R)^{2l}, R = (2l-1)(2D+2l).
Proof. Take an optimal pi of cost <= |H|+l (if YES). Let X = endpoints of U and define x ~ x' if x, x' share a cluster or {x,x'} in U.
Each ~-step has G-distance <= 2D+2l inside L_B (Lemma 5 alpha; deleted pairs have distance 1), and a ~-component has <= 2l vertices, so
it lies in Ball_{G[L_B]}(x_j, R) for any of its vertices x_j. Every component contains a pair of U, so r <= |U| <= l components.
Algorithm: guess roots x_1..x_r in L_B (<= n^l), guess U' = a set of <= l pairs inside the union of the balls ((r |ball|)^{2l} choices),
run the DP of Thm A (width w' bags, domain 7) with U' as the unpacked edits, accept iff min cost + |U'| <= |H|+l. Sound for every U'
(Thm A is exact for the U' it is given); complete because the true U is enumerated. QED
Same proof: if the edited pairs of some solution fall into <= r0 balls, time n^{r0} g * 7^{w+1}. So r = 1 costs n * g: FPT-like.
Status: still XP. It shows the exponent is the number of components, not 2l, and that l bounds that number but not the cost of
placing the roots.
### Theorem 8 (FPT in (l, D, omega))  PROVED (standard partition DP over a tree decomposition; the DP itself is not coded)
Let BIG = hard components containing a vertex with d_h > B; reject if one is not clean. V' = V - BIG; Gamma' = graph on V' whose edges are
all G-edges and all soft pairs inside V'; omega = tw(Gamma'); Z = 3D+3l+1. Then OPT(G) <= |H|+l iff
|E| + sum_{K in BIG} f(K) + m <= |H|+l with m = min over partitions pi' of V' whose clusters are Gamma'-connected and have <= Z vertices of
sum_C f(C). Computable in time (Z omega)^{O(omega)} n.
Proof. (<=) every such pi' with the clusters K in BIG is a partition of V, cost by Lemma 1(a). (=>) Lemma 5: an optimal pi has the BIG
components as clusters (gamma); every other cluster is Gamma'-connected (alpha: G-connected inside L_B, hard edges lie in Gamma' because
both ends are in L_B subset V'; beta: a single clean component is hard-connected, several components are joined by soft pairs, all in
Gamma') and has <= Z vertices (exceptional: 2D+2l+1; beta with several components: 2D; one clean component K: |K| <= d_h+D+1 <= B+D+1 since its non-adjacent pairs are soft). DP: nice tree decomposition of Gamma'; state = partition of the bag +
size (<= Z) of every open cluster; introducing v into a cluster of current size s adds s - 2 (#G-neighbours of v in its block); this uses
that every G-edge lies in a bag (an edge to an already forgotten vertex cannot exist). A cluster whose last bag vertex is forgotten is
closed; connectivity guarantees no later vertex joins it. QED
Here no pointer is guessed: far merges are carried by Gamma'-connectivity (connected subtrees). omega is the only obstruction.

## 5. Where the connected-subtree property breaks
### Prop 9 (FAILS: omega is not bounded by f(l, w))  VERIFIED-BY-COMPUTER (construction) + CITED (tw of cubic expanders)
Family: F a cubic graph. Each edge e = (u->v) of F becomes a P3 (a_e, b_e, c_e) with a_e in the clique K_u, b_e, c_e in the clique K_v and
one cross edge a_e b_e. H = all these P3s: vertex-disjoint (w = 0, c = 1, D = 2), induced, pair-disjoint, pair-maximal (an induced P3
leaves a clique through a cross edge, which is soft). The partition {K_v} has cost |H| (only the cross edges), so OPT = |H|, l = 0, YES.
Every K_v is one clean hard component of <= 5 vertices (hard degree <= 4 = 2D+3l, no BIG). Contracting the connected K_v gives F, so
tw(Gamma') >= tw(F) = Omega(|F|) for cubic expanders (Grohe-Marx 2009, CITED). code/w5b_family.py checks, for K4, K33, Petersen and random
cubic graphs on 12..40 nodes (n up to 180): induced, pair-disjoint, pair-maximal, cost = |H|, hard degrees. It also checks that BFK
Rule 6.2 (S = {ab}) applies to every P3, so the family is trivial for BFK; the obstruction concerns the representation, not hardness.
Consequence. Hard edges between far vertices are free (not in the budget) and arbitrary in number; they make the footprints of clusters
disconnected in any decomposition of I(H) or W, and put an expander into Gamma'. At l = 0 they are harmless because the triple-CSP of
Thm A needs no identities; with l > 0 an exceptional cluster can be assembled from such far vertices (c = 1 example of wave2_c1_w1 3.4(iv):
two vertex-disjoint P3s joined by one hard edge that is deleted). So the DP cannot use one decomposition for both purposes.
### Where inserted pairs matter (answer to the "far merge" question)
* Pure far merges (a cluster joined only by inserted pairs): do not occur in optimal partitions (Lemma 4).
* Far inserted/deleted pairs occur (verified, w2a_far) and are bounded: <= l pairs, <= 2l exceptional clusters, each of size <= 2D+2l+1,
  all vertices of hard degree <= B, and the partner is determined up to g(l,D) choices by the other end (Cor. 6).
* The budget does bound the NUMBER r <= l of such components, but their positions are n^r (Thm 7). Guessing them needs either
  position-independent verification (Prop 10) or a bounded-width carrier (Prop 9 says none exists in general).
### Prop 10 (colour coding does not repair it)  analysis, FAILS as a method; not a theorem
* A proper colouring of the R-th power of G[L_B] with k = (B+D)^R + 1 colours is computable greedily (no randomness). It names every
  vertex of a ball by its colour, so ONE site can validate a whole exceptional component (cluster membership, all its pairs,
  C0-type triples) by G-look-ups. This is stronger than the usual "guess colours of endpoints".
* What a site cannot do is make the other vertices of the same cluster (up to Z of them, in far bags) declare the same cluster, and make
  their P3 variables x_j agree with it. Counting the cost at the minimum-colour vertex of C undercounts as soon as that vertex is
  default: the DP then accepts a configuration in which the P3s of the other members do not know that they lie in C; its cost is lower than
  that of any partition. Charging 1/|C| at every declaring site fails in the same way. Fixing it needs the identity of the other sites
  across bags, i.e. exactly the n^r guess.
* A sound variant exists only for a decomposition in which all P3s meeting one ball share a bag: tree decomposition of
  I+ = I(H) + {P ~ P' : V(P), V(P') within distance R in G[L_B]}. Then Thm A with y-variables for pairs inside bags gives
  FPT in (l, D, tw(I+)). By Prop 9, tw(I+) is not bounded by (l, w) (same family).

## 6. Exact missing lemma and best partial results
MISSING (ML, CONJECTURE). There is a polynomial-time reduction R (safe: preserves "OPT <= |H|+l", never increases l, maybe w) after
which either the instance is rejected or tw(Gamma') <= f(l, w). Equivalent form needed by Thm 8: a rule that removes "inert" hard
components (those that stay intact in every solution and whose incident P3s have a forced local solution) from G and H.
* For c = 1 (w = 0) the BFK rule 6.2 is such a rule and leaves <= 2l P3s (CITED, Thm K); it removes the family of Prop 9.
* For w >= 1 the BFK counting fails (wave2_c1_w1 3.4(i): a free external edit can be an internal pair of another member). Components of
  bounded size are covered by Thm K. The open case is the same as Q1/Q2 there (trees of unbounded size, long chains), now with the sharper
  fact that the pointer itself has a bounded candidate set and only its location is unbounded.
* Cheaper necessary subgoal (CONJECTURE, not tested at size): the ~-components of an optimal solution are within distance <= 2 in G
  (Cor. 6 observation); together with Thm 7 this would give n^{l} * 7^{w+l+1} * poly.
BEST PARTIAL RESULTS (all PROVED above): (1) Lemma 1 decomposition; (2) Lemmas 2-5, Cor. 6 structure and bounded candidate lists;
(3) Thm 7: n^{l} g(l,D) 7^{w+l+1}; (4) Thm 8: (Z omega)^{O(omega)} n, omega = tw(Gamma') (not bounded in (l,w), Prop 9);
(5) Thm K of wave2 remains the only FPT result with w small and arbitrary hard structure.

## 7. Test log
T-CP/HE/H2/EX/NF/F1/BIG  w5b_tests.py  seeds 1,2,3,11,12,13; n = 6..11; instances 5,300; all partitions of cost <= |H|+3: 450,000+;
   edited-pair checks 1.5M; optimal partitions 17,000+; non-trivial BIG cases 391; 0 failures. Families: G(n,p), planted, sparse, hub-like.
T-ILP    w5b_ilp.py    seeds 5,6; 100 instances n = 11..13; EX1 (d_h > D+l means no edit of a hard pair at v): 76 infeasibility proofs;
   EX2 (max cluster size with an edited hard pair at u): 897 maximisations, 0 above 2D+2l+1, max slack -5.
T-CAND   w5b_cand.py   seeds 7,8; 2,000 instances n = 6..10; 4,640 edited pairs of optimal partitions; max distance (inserted) 2; 0 failures.
T-FAM    w5b_family.py 7 graphs up to n = 180; all conditions hold.
Limits: n <= 13; instances random (Gamma' tw is small at this size, so Prop 9 is a construction, not an experiment);
Theorem 8's DP was not implemented; Lemma 3 slack at n <= 13 is only attained in the finer forms.
