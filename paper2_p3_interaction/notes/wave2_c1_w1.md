# Wave 2 (agent c1/w1): Conjecture C1, the cases w = 0 and w = 1

Labels: PROVED (full proof here), VERIFIED-BY-COMPUTER (small cases only, not a proof), CITED (relies on a
published statement, not re-proved here), CONJECTURE, FAILS. Notation as in analysis_01 / fpt_attempt:
H maximal modification-disjoint packing (no induced P3 of G has all three pairs unpacked), V(H), R = V \ V(H),
Q(a) = vertices sharing a packed P3 with a, q = max_a |Q(a)| <= 2(w+1), excess budget l, k = |H| + l.
Code (all in code/): w2a_lib.py (numba subset-DP exact cluster editing with an optional cluster-validity
mask; n = 11 takes ~0.1 s), w2a_c1.py, w2a_c1b.py, w2a_shrink.py, w2a_bfk.py, w2a_comp.py, w2a_far.py.
Instances: n = 9..11, "structured" (random core + attached cliques with disjoint attachment sets, then a
random maximal packing preferring core vertices) and G(n,p); about 3250 instances with |R| >= 2, half of them
with OPT > |H|.

## 0. Summary

| item | status |
|---|---|
| Twin lemma: adjacent true twins lie in one cluster in EVERY optimum | PROVED (1.1) |
| **C1**: for maximal H some optimum has no unpacked edit with both ends in R (indeed every cluster meets R in nothing or in exactly one whole twin class) | **PROVED** (1.2) + VERIFIED (3247 inst.) |
| Normal form, weights \|X\|, frozen classes (\|X\| > l), shrink rule (\|X\| >= l+2 -> l+1), kernel n <= (l+2)\|V(H)\| | PROVED (1.3) + VERIFIED |
| Pointer partners: every R-pointer is unary, or local, or reduces to a V(H)-V(H) pointer | PROVED (1.4) + VERIFIED |
| w <= 1 implies c <= 2 (forest I(H): every vertex in <= 2 P3s) | PROVED (3.1) |
| w = 0 (c = 1): consistent with van Bevern et al. (BFK); explicit pointer-enumeration FPT given BFK's reduction | PROVED modulo CITED BFK; checks VERIFIED (2) |
| **FPT in (l, s)**, s = max size of a connected component of I(H); in particular I(H) a perfect matching (w = 1, c = 2) is FPT | PROVED modulo CITED BFK Thm 6.6; VERIFIED (3.2) |
| FPT in (l, w) for trees with unbounded components | OPEN. Exact obstruction identified (3.3-3.4); no counterexample, no algorithm |
| Correction to fpt_attempt.md Sec. 3 / NEXT 1: "I(H) a perfect matching is open" | FAILS (it is FPT, 3.2) |

## 1. Conjecture C1

### 1.1 Lemma T (twins)  PROVED
If u != v and N[u] = N[v] (adjacent, same closed neighbourhood), every optimal partition P puts u, v in one cluster.
Proof. Suppose u in C_u, v in C_v, C_u != C_v. For d not in {u,v}: ud in E iff vd in E. Put
f(y, D) = #non-neighbours of y in D - #neighbours of y in D (D contains neither u nor v); then f(u, D) = f(v, D).
Moving u into C_v changes the cost by D1 = -f(u, C_u - u) + f(u, C_v - v) - 1 (the -1: the edge uv stops being deleted).
Moving v into C_u changes it by D2 = -f(v, C_v - v) + f(v, C_u - u) - 1. D1 + D2 = -2 < 0, so one move improves, contradicting
optimality. QED
This is strict, so it holds for every optimum. It is the critical-clique fact (Guo 2009) now proved, closing item (iv) of
fpt_attempt Lemma 4. VERIFIED: 2866 twin pairs in R-classes, forbidding their separation costs > OPT.

### 1.2 Theorem (C1, strong form)  PROVED
Let H be maximal, X_1..X_m the components of G[R] (twin classes, Lemma 4 of fpt_attempt), A_i = N(X_i) \ X_i. Then A_i is a subset of V(H)
and the A_i are pairwise disjoint (Lemma 4: non-adjacent r, s in R have disjoint neighbourhoods). Then some optimal partition has every cluster
containing no vertex of R, or exactly one entire class X_i. Hence some optimum has no edit with both endpoints in R.
Proof. Take an optimum P. By Lemma T each X_i lies in one cluster. An edit inside R is a deletion inside some X_i (impossible, X_i intact) or an
insertion r s, r in X_i, s in X_j, i != j (R-R edges exist only inside classes). So one exists iff some cluster C contains two classes.
Let C contain X_1, X_2, ...; x_i = |X_i|, a_i = |A_i cap C|. Moving X_1 out of C into a cluster of its own changes the cost only through the pairs
(X_1, C - X_1): the x_1 |C - X_1 - A_1| inserted non-edges disappear, the x_1 a_1 satisfied edges become deletions (N(X_1) - X_1 = A_1).
Change: Delta_1 = x_1 (x_1 + 2 a_1 - |C|). The A_i are disjoint subsets of V(H) and the X_i lie in R, so |C| >= x_1 + x_2 + a_1 + a_2.
If Delta_1 > 0 and Delta_2 > 0 then |C| < x_1 + 2 a_1 and |C| < x_2 + 2 a_2, so x_2 + a_2 < a_1 and x_1 + a_1 < a_2; adding gives x_1 + x_2 < 0, absurd.
So moving X_1 or X_2 out is free, the result is again optimal, and the number of pairs of classes sharing a cluster strictly drops. Iterate. QED
VERIFIED (w2a_c1.py, seeds 11-14): 3247 instances, n = 9..11, 1599 with OPT > |H|. The optimum under the constraint "every cluster meets R in nothing
or one whole class" equals the unconstrained optimum, and Lemma 4 + maximality hold on every instance. (T7 of fpt_attempt covered only n <= 8.)

### 1.3 Corollaries (normal form NF)  PROVED unless noted
Let S be an optimum in NF with |S| <= |H| + l. Unpacked edits number <= l (fpt_attempt Lemma 0).
 (N1) Unpacked edits have three types: HH = pair inside V(H) lying in no common P3; RH-del = deletion of X x {a}, a in A_X; RH-ins = insertion of
      X x {a}, a in V(H) \ A_X. Every pair with an endpoint in R is unpacked.
 (N2) X is inside one cluster, so X x {a} is edited for all |X| copies or for none. Excess = #HH + sum_X |X| |A_X xor B_X|, B_X = cluster(X) cap V(H).
 (N3) Frozen classes. If |X| > l then no optimum of cost <= |H| + l edits any pair at X: cluster(X) = X cup A_X exactly (any edit at X costs >= |X| unpacked edits).
      VERIFIED: 4984 tests (|X| > OPT - |H|; optimum with cluster(X) = X cup A_X forced equals OPT).
 (N4) Each a in V(H) has neighbours in at most one class. Non-edges inside A_X are packed pairs (all-unpacked P3 otherwise). For a in V(H) \ (A_X cup X) and
      z in A_X, az in E forces az packed (else {x,z,a} is an all-unpacked P3). Hence G[A_X] is a clique minus packed pairs and |N(a) cap A_X| <= q. VERIFIED.
 (N5) Shrink rule. If |X| >= l+2, deleting one vertex x of X gives an equivalent instance (G - x, H, l).
      Proof. H stays a maximal packing of G - x. (=>) a budget solution of G may be replaced by an optimal one; by N3 it does not touch X; delete x.
      (<=) take an optimal solution of G - x, cost <= |H| + l; |X - x| >= l+1 > l, so by N3 it does not touch X - x; add x as a twin, cost unchanged. QED
      A class with A_X empty is a clique component and can be deleted (OPT unchanged). VERIFIED: 10,863 shrink tests (l' = 0..4, classes up to size 5) and
      5,139 isolated-class deletions (w2a_shrink.py).
      Kernel-type consequence: classes with A_X non-empty number <= |V(H)| (disjoint non-empty A_X), so after the rules |R| <= (l+1)|V(H)| and
      n <= (l+2)|V(H)| <= 3(l+2)|H|. R is polynomially bounded by V(H) and l.

### 1.4 What C1 does to the candidate pointer partners (Lemma 4 twin classes)  PROVED, VERIFIED
 (P1) RH-del (X, a): the partner of a is cl(a), unique. It is a UNARY flag on the vertex a (cost |X|) living at a's bag node. No identity problem.
 (P2) Detachment is never useful: if B_X is non-empty and B_X cap A_X is empty, moving X out alone saves |X||B_X| > 0 (the deletions |X||A_X| are paid either way).
      So in an optimum B_X cap A_X = empty forces B_X = empty. (A computer test of this one claim is vacuous when X is already alone; the proof is the cost count.)
 (P3) RH-ins (X, a), a in B_X \ A_X: by P2 there is z in A_X cap B_X, and for every such z the pair za is a packed edge (then a in Q(z), at most q candidates,
      within I(H)-distance <= 1 of z) or an unpacked non-edge that is INSERTED, i.e. an HH pointer (z, a) is present (az edge forces packed, N4).
      So an RH-ins pointer is local to z or inherits its identity problem from an HH pointer.
 (P4) Wide classes cannot receive insertions: if |A_X| >= l + q there is an optimum with B_X subset A_X.
      Proof. alpha = |A_X cap B|, delta = |A_X \ B|, beta = |B \ A_X| >= 1. Unpacked edits: |X| beta insertions X x (B \ A_X), |X| delta deletions, and for
      pairs (z,b), z in A_X cap B, b in B \ A_X: edge implies packed (<= q per b), non-edge and unpacked implies inserted, so >= beta (alpha - q) of them.
      l >= beta (|X| + alpha - q) + delta. If alpha >= q this gives |A_X| <= l + q - 1; if alpha < q then |A_X| = alpha + delta <= q - 1 + l. QED
      VERIFIED: 900 classes with |A_X| >= l + q (w2a_c1b.py, seeds 61-64).
 Net effect: every pointer touching R is unary, local (<= q candidates), or reduces to an HH pointer. The whole identity problem sits on pairs of V(H) vertices in
 no common P3 (the "Cover case", R empty). This sharpens open problem 3 of fpt_attempt Sec. 5. HH pointers are not restricted further by C1 (Lemma 2 of fpt_attempt only).

## 2. Case w = 0 (c = 1: vertex-disjoint P3s) and consistency with BFK (arXiv 1512.04047)
CITED: BFK Reduction Rule 6.2 / Lemma 6.3 (sound for any induced subgraph G[W] with an optimal solution S of G[W] satisfying the neighbourhood conditions),
Lemmas 3.2 and 6.5 (if the rule applies to no packing graph, #packing graphs <= 2 l'), Thm 6.6 / Cor. 6.7 (FPT in (t, l') for vertex-disjoint cost-t packings
of induced subgraphs; kernel (4t+2) l' vertices).
VERIFIED-BY-COMPUTER (w2a_bfk.py E2, 5,870 random c = 1 instances, n = 8..11, seeds 21-24):
 (a) whenever Rule 6.2 applies to a single packed P3 (all three candidate S, condition checked literally), OPT(G xor S) = OPT(G) - 1;
 (b) in all 4,420 instances with H non-empty where the rule applies to no packed P3, |H| <= 2 (OPT - |H|). With n <= 11, |H| <= 3, so this only checks my
     reading of the statements; it is not a strong test.
Pointer-enumeration FPT for c = 1 (PROVED modulo CITED BFK): apply Rule 6.2 exhaustively (polynomial); each application keeps a vertex-disjoint packing and l' <= l.
Then |H| <= 2l, |V(H)| <= 6l. By C1 (1.2-1.4) the unpacked edits of some optimum are HH pairs (<= (6l)^2 candidates), RH-del flags (<= 6l), RH-ins pairs
(class, a) (<= 6l classes x 6l vertices). U ranges over <= l-subsets of O(l^2) candidates; for each U the CSP of Theorem A has independent variables (w = 0) and
is polynomial. Total l^{O(l)} n^{O(1)}: weaker than O(4^l l^3 + n + m) and consistent with it.
Why c = 1 works in our language: a pair with exactly one end in V(P) is automatically unpacked, hence an excess edit; BFK's counting then bounds |V(H)| by O(l)
and the identity problem disappears (candidate set O(l^2), independent of n). That counting needs "external pairs are unpacked", which fails once P3s share vertices.

## 3. Case w = 1 (I(H) a forest) and general components

### 3.1 w <= 1 lies inside c = 2  PROVED
Three P3s through one vertex pairwise intersect, so they form a triangle of I(H). A forest I(H) therefore has every vertex in <= 2 P3s, max degree <= 3, c <= 2.
So w = 1 is a subclass of LPS's open c = 2 case (their "grid-like structures" need cycles in I(H); the forest case is not settled in LPS).

### 3.2 Theorem K (FPT in (l, s))  PROVED modulo CITED BFK Thm 6.6
Let s = max_C |V(C)| over components C of I(H) (V(C) = union of the vertex sets of the P3s of C). CEaMP is FPT in (l, s).
Proof. (1) Components have pairwise disjoint vertex sets W_j = V(C_j): P3s of different components share no vertex. (2) Put H_j = G[W_j], tau_j = tau(G[W_j])
(brute force, exponential in s only). The P3s of C_j are pair-disjoint induced P3s of H_j, so tau_j >= |C_j|; e_j = tau_j - |C_j|. (3) "OPT(G) <= |H| + l" is the Cluster
Editing instance (G, k = |H| + l) with the vertex-disjoint cost-t packing {H_j}, t = max tau_j <= C(s,2); the excess above the lower bound h = sum tau_j is
l' = l - sum e_j (NO if negative). (4) BFK Thm 6.6(i) decides it in time f(t, l') n^{O(1)} (Cor. 6.7: 1.62^{(2t+1) l'} plus polynomial). QED
Why this escapes LPS's obstruction: a pair with exactly one end in W_j cannot be a packed pair (a packed pair lies inside one P3, so inside one W_j); BFK's
"external pair is excess" survives at component level, and dies for single P3s of a non-matching packing (LPS lines 142-144).
VERIFIED-BY-COMPUTER (w2a_comp.py, 3,179 instances, n = 8..11, components with <= 3, 5, 6, 7 vertices, Rule 6.2 applied to components with every optimal S of
G[W_j], 741 applications): the rule is sound (OPT drops by exactly tau_j); in the reduced instance with k = OPT, l' >= 0 always and #components <= 2 l' always
(ratio 2 is attained, 2,774 non-trivial reduced instances).
Special cases: s = 1 is BFK itself. I(H) a perfect matching (each P3 meets exactly one other, w = 1, c = 2, |V(C)| = 5) is FPT; likewise any packing with
components of <= s P3s. After exhaustive Rule 6.2 on components, <= 2 l' components remain and n <= (4t+2) l', so for bounded s the identity problem is solved by
kernelisation (all candidate partners among O(s^2 l) vertices). (l, s) is a weaker parameter than (l, w) since w <= s - 1, but it contains these cases.

### 3.3 Subtree modules in a forest  PROVED modulo CITED BFK Lemma 6.3
Let I(H) be a forest, F a connected sub-packing, W = V(F). A P' not in F meets W in <= 1 vertex (two vertices would lie in two members of F joined by a path in F,
closing a cycle with P' in I(H); both cannot lie in one member since pairs of one member belong to it). So no pair of P' lies in W. Hence Rule 6.2 for G[W]
(Lemma 6.3 is a statement on Cluster Editing alone) gives an equivalent instance (G xor S, H - F, l - (tau - |F|)) in which H - F is still a modification-disjoint
packing of the same induced P3s. VERIFIED-BY-COMPUTER (w2a_bfk.py E6, 5,644 forest instances, sub-packings of <= 3 P3s, 7,112 rule applications, 4,652 with |F| >= 2,
3,263 with internal excess tau > |F|): OPT drops by exactly tau; the remaining P3s are untouched. So modules can be peeled inside a tree, not only whole components.

### 3.4 Why the win-win does not give FPT for trees (exact obstruction)  analysis, not computer-checked
(i) Partition a tree component into subtrees (modules) F_1..F_r. Each vertex lies in <= 2 of the W_j (c <= 2) and |W_j cap W_k| <= 1, so BFK's counting of inner edits
goes through, EXCEPT that an edit (u,v) with u in W_j cap W_k, v in W_k \ W_j is external for j but internal for k, where it may be a mandatory non-excess edit. Such "free"
external edits occur exactly at the cut vertices between modules, and a tree forces cuts. So the certificate of Lemma 6.5 for a module can be met by a free edit, the
bound #modules <= 2 l' is lost.
(ii) At component level the win-win still leaves <= 2 l' active components, each a tree of unbounded size. Inside a component, R is empty (W = V(C)): the "Cover case".
Computing tau(G[V(C)]) is itself the problem with R empty, so the reduction is circular: deciding tau(G[W]) <= |C| + e for a tree packing covering all of W is the core open
problem, and the one case where C1 gives nothing (only HH pointers).
(iii) Independence of pointers. TRUE and PROVED (Lemma 3 of fpt_attempt): pointers y = (u,v), y' = (u',v') with disjoint vertex sets lie in no common triple (a triple
has 3 vertices), so they interact only through the x-variables of nearby P3s. The effect of y on the CSP depends on v only through (G(uv), adjacency profile of v on Q(u)) and on
u through the profile of u on Q(v): at most 2 * 2^{2w+2} values per end, plus the global condition C0 on triples {u,v,z}, z outside the environments.
It does NOT give additivity of costs: the min-sum DP transmits a finite state along I(H), and a defect in a chain of length m can be repaired by a pointer anywhere along it,
so far-apart pointers are coupled (this gives ties, not an obstruction by itself). The real obstruction is the matching of the two ends of one pointer: the relation
"(u,v) realises profiles (rho, rho') and passes C0" is defined by G and links two arbitrary nodes of the tree. For c = 1 and bounded components kernelisation puts it on O(l)-size
sets; for trees nothing bounds it. (Not tested by computer; reasoning only.)
(iv) Far HH pointers are genuinely necessary. VERIFIED (w2a_far.py): 139 of 2140 forest instances (n <= 10) have no optimum avoiding unpacked edits between P3s at I(H)-distance
>= 2; the smallest has n = 7: vertex-disjoint P3s (0,5,1), (2,3,6) joined by the edge 35, OPT = |H| + 1 only by deleting 35. So "pointers can be taken local in I(H)" FAILS, already
for c = 1; what makes c = 1 easy is |V(H)| = O(l) after reduction, not locality.

## 4. Open questions remaining (all on the Cover case, R empty)
 Q1. CONJECTURE: for a tree I(H) (subcubic forest, w = 1) CEaMP is FPT in l. Evidence: only Theorem K for bounded components; no hardness found. First target: I(H) a path
     (each P3 meets the next in one vertex): is "tau(G) <= |H| + l" FPT when G has the 2m+1 vertices of the chain? Attack: after exhaustive removal of subtree modules (3.3) show
     |V(H)| = O(l poly(w)), replacing BFK's counting by a charge to cut vertices (3.4 (i)).
 Q2. Is there a bounded-size module rule for trees (a module of <= g(l) P3s is removable once |H| > f(l))? Experiments at n <= 11 are uninformative (a chain of m P3s needs
     2m+1 vertices, so m <= 5); needs n ~ 20-30 and the Theorem A DP as exact solver (U guessed, l <= 2).
 Q3. Can HH pointers have a candidate bound f(l, w) after the R-reductions of 1.3? C1 shows R is not the source of the difficulty.
 Q4. Hardness (not attempted): any W[1]-hardness proof for (l, w) must use HH pointers on V(H), since every R-pointer is unary, local, or inherits from an HH pointer (1.4).

## 5. Test log (all pass)
 T-twin   w2a_c1.py      2866 twin pairs (Lemma T)
 T-C1     w2a_c1.py      3247 instances n = 9..11 (seeds 11-14): NF optimum = OPT; Lemma 4; N3 (4984 tests); N4
 T-D      w2a_c1b.py     seeds 61-64 x 2500: P4 on 900 classes with |A_X| >= l + q
 T-shrink w2a_shrink.py  seeds 51-54 x 1500: N5 (10,863) and isolated-clique deletion (5,139)
 T-BFK    w2a_bfk.py     E2 (c = 1): 5,870 inst.; E6 (forest modules): 5,644 inst.; seeds 21-24
 T-comp   w2a_comp.py    Theorem K reduction at component level, 3,179 inst., seeds 41-44
 T-far    w2a_far.py     far HH pointers necessary (139 / 2140, seed 1)
Limits: n <= 11; nothing here probes large w or large l; Theorem K and Sec. 2 rely on BFK's published proofs (CITED), checked only by the experiments above.
