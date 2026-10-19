# Audit 2 (hostile referee, wave 3): twin lemma, C1 strong form, normal-form consequences

Target: notes/wave2_c1_w1.md Sec. 1 (1.1-1.4), with Lemma 4 of fpt_attempt.md. Own code only, written from scratch
(`code/a3a_lib.py`, `a3a_exh.py`, `a3a_ilp.py`, `a3a_cex.py`, `a3a_n2.py`); none of the authors' scripts is imported.
Verdicts: CORRECT / GAP / WRONG. Nothing is WRONG. Several GAPs of statement, one real gap outside the three
claims (Sec. 2 of wave2_c1_w1, c = 1), and one strengthening that removes the "some optimum" and "termination" issues.

## 0. Conventions and the exact hypotheses

* Unpacked pair = pair of vertices not contained in a common member of H. Pairs with an end in R = V \ V(H) are always unpacked.
* "Maximal" is used in wave2_c1_w1 as: no induced P3 of G has all three pairs unpacked, i.e. H is maximal among
  PAIR-disjoint (modification-disjoint) induced P3 packings. This is the only hypothesis under which Lemma 4 holds.
  It is NOT maximality among vertex-disjoint packings (see 5). The text states it in the header line; the theorem
  statements 1.2-1.4 should repeat it.
* Optimum = partition of V minimising the cost (number of edited pairs); OPT >= |H| (Lemma 0), l_opt = OPT - |H|.
* Lemma 4 (re-proved by me, same argument): r,s in R adjacent, z adjacent to exactly one of them gives induced P3
  {r,s,z}, all pairs touch R; r,s non-adjacent with a common neighbour z gives induced P3 r-z-s. Both contradict maximality.
  Consequences: components of G[R] are cliques of true twins (transitivity of N[.]), N(X) \ X is inside V(H), and different
  classes have disjoint attachment sets A_X. CORRECT. Note: two classes with equal (or overlapping) non-empty attachment sets
  cannot both survive in R; with such a graph the maximal H necessarily swallows one of them (checked: a "hub" family
  produces such graphs, R then loses a class). So "classes with equal attachment sets" only occur inside V(H), where
  the twin lemma (not C1) governs them; both were tested.

## 1. Twin lemma (1.1): CORRECT

Statement: u != v, uv in E, N[u] = N[v]. Then u,v share a cluster in every optimum. No hypothesis on H.
Recomputed: with f(y,D) = |non-nbrs of y in D| - |nbrs of y in D|, moving u into C_v changes cost by
f(u,C_v - v) - 1 - f(u,C_u - u); symmetric for v; f(u,.) = f(v,.) off {u,v}; sum of the two changes is exactly -2.
So for ANY partition separating a twin pair one of the two single-vertex moves strictly improves. Strict, hence every optimum.
Edge cases (singleton C_u, C_v) fine. Computer: all true-twin pairs of G (not only those in R), exhaustive n = 7..9
(14,800 pairs in seeds 21-24; 0 separated in any optimum) and ILP "maximise separations s.t. cost <= OPT" n = 10..13 (41,900 pairs, 0).
Fix: none. (Optional: say that this is Guo's critical-clique property and that it is in fact strict.)

## 2. C1, strong form (1.2): CORRECT, and the proof can be strengthened (no ties, no iteration)

Restated: H pair-maximal, X_1..X_m the components of G[R], A_i = N(X_i) \ X_i (subsets of V(H), pairwise disjoint).
Claim: some optimum has every cluster meeting R in nothing or exactly one whole class; hence no unpacked edit with both ends in R.
Line by line:
* Edits inside R are deletions inside a class or insertions between two classes (R-R edges exist only inside classes):
  CORRECT. With the twin lemma deletions inside a class do not occur, so an R-R edit exists iff a cluster holds two classes.
* Exchange formula: moving a whole class X_1 out of cluster C into a singleton cluster changes cost by
  Delta_1 = x_1 (x_1 + 2 a_1 - |C|), a_1 = |A_1 cap C|. CORRECT: N(X_1) \ X_1 = A_1 (Lemma 4), so the inserted non-edges are
  x_1 (|C| - x_1 - a_1) and the satisfied edges that become deletions are x_1 a_1. Verified against the true cost
  difference on 1.1 million (class pair, partition) samples (a3a_exh.py, "exchange_formula_tests"), 0 mismatches;
  this includes clusters that hold two or more classes plus vertices of V(H) outside A_1 cup A_2 and vertices of A_3.
* Bound |C| >= x_1 + x_2 + a_1 + a_2: CORRECT (A_1, A_2 disjoint, inside V(H), disjoint from R).
* The authors conclude "not (Delta_1 > 0 and Delta_2 > 0), so one move is free (Delta = 0), the result is again optimal,
  iterate; the number of class pairs sharing a cluster drops". Valid, but a weaker conclusion than what the same two lines give.
  STRICT FORM (new, 3 lines): Delta_1 >= 0 and Delta_2 >= 0 give |C| <= x_1 + 2a_1 and |C| <= x_2 + 2a_2; with the bound on |C|:
  x_2 + a_2 <= a_1 and x_1 + a_1 <= a_2; adding, x_1 + x_2 <= 0, absurd. So one of the two moves has Delta < 0, i.e.
  strictly improves. Therefore NO optimum has a cluster containing two classes: C1 holds for EVERY optimum, not only for some.
  Consequences for the questions asked: (a) ties between optima are irrelevant (there are no optimal R-R merges at all);
  (b) the move is never executed on an optimum, so termination of iterating it is vacuous (the potential argument is correct
  but unnecessary); (c) "some optimum" can be replaced by "every optimum", which is what every later use (N3, P2, P4) needs.
* Computer. Exhaustive (all set partitions, n = 7..9; seeds 21-24 with the strict test: 2,825 instances with |R| >= 2, pair-maximal H
  from random and core-first greedy packings, 724 with OPT > |H|; seeds 11-14 with the weak test: 2,952 more; coverage from a3a_n2.py,
  1,436 instances: 1,104 with a class of size >= 2, 1,089 with >= 2 classes, 318 with >= 2 classes of non-empty A): (i) strict inequality
  "Delta_1 >= 0 and Delta_2 >= 0 never happens" on 4.7 million (seeds 21-24; 4.8 million more with the weak form)
  partitions in which two classes sit whole in one cluster (all partitions, not only optima): 0 failures; (ii) number of instances
  with an optimum that has an R-R edit: 0 (no tie ever); (iii) the move DAG from up to 80 optima per instance was never entered.
  ILP (HiGHS, n = 10..13, 8,150 instances, R >= 4 in 3,400, >= 3 classes in 2,000, l_opt = 0..5, families: planted packing with
  attached classes (disjoint / hub / random attachments), noisy clusters with classes on shared hubs, blow-ups, G(n,p)):
  max over optima of #R-R edits is 0 in all 8,150 (so "every optimum", not only "min over optima").
* Hypothesis necessity (a3a_cex.py, a3a_exh.py hyp=vertexmax/nonmax, a3a_ilp.py vertexmax): with H maximal only among
  VERTEX-disjoint packings C1 FAILS. Smallest: n = 5, edges 03,13,14,24,34, H = {(0,3,1)} (vertex-maximal: remaining vertices
  2,4 are 2 vertices), R = {2,4}, OPT = 2 and every optimum deletes the R-R edge 24 (also forced 03). Lemma 4 fails there:
  {3,4,2} is an induced P3 with all pairs unpacked. Rates: 86 / 2,512 vertex-maximal ILP instances (n <= 12) have no optimum
  without an R-R edit, 43 / 1,283 in the exhaustive run. With an arbitrary non-maximal H (e.g. H empty, n = 4) it fails trivially.
Verdict: CORRECT. Fix: state "maximal = no induced P3 with all three pairs unpacked"; use the strict form (every optimum),
delete the iteration/termination sentence; the proof of Lemma 4 should be in the same section (currently only cited).

## 3. Consequences (1.3, 1.4)

### N1, N2 (types of unpacked edits): CORRECT as identity, GAP in wording
N1 fine (every pair with an end in R is unpacked; R-R edits absent by Sec. 2). N2: "Excess = #HH + sum |X||A_X xor B_X|" is
the number |U| of UNPACKED edits, not "the excess"; the excess in the sense of fpt_attempt Lemma 0 is sum(|x_i|-1) + |U|.
The identity |U| = #HH + sum_X |X||A_X xor B_X|, B_X = cluster(X) cap V(H), was checked on every optimum of 1,436 exhaustive instances
(3,008 optima, 0 failures). Fix: rename to |U|; state |U| <= l - sum(|x_i|-1) <= l.

### N3 frozen classes (|X| > l): CORRECT, slightly stronger than stated
Restated: if |X| > l_opt (in particular if |X| > l >= l_opt) then EVERY optimum has cluster(X) = X cup A_X. Proof is right:
by twin lemma + Sec. 2, any edit at X is made |X| times (the pairs X x {v} are edited for all of X or none), costing >= |X| > l unpacked edits,
contradicting |U| <= l_opt. The note's phrase "no optimum of cost <= |H|+l" is fine but the quantifier is "every optimum" and the
hypothesis is "the instance is a YES instance (OPT <= |H|+l)"; for NO instances nothing is claimed or needed.
Computer: all optima of 5,100 exhaustive classes (seeds 21-24) and 9,459 ILP classes (maximise edits at x_0 s.t. cost <= OPT): 0 failures.
Fix: state the YES-instance hypothesis.

### N4 (each vertex of V(H) adjacent to at most one class; G[A_X] clique minus packed pairs; |N(a) cap A_X| <= q): CORRECT
Part 1 is just disjointness of the A_X (Lemma 4), trivial. Part 2: z,z' in A_X non-adjacent gives induced P3 z-x-z' with xz, xz' unpacked,
so zz' is packed; part 3: a not in A_X cup X, z in A_X, az in E gives induced P3 x-z-a, so az is packed. Both correct, hypothesis: pair-maximal H.
Computer: 11,600 (class, vertex) checks (seeds 21-24), 0 violations. Fix: none; say a is taken in V(H).

### N5 shrink rule (|X| >= l+2): CORRECT, threshold is tight
Restated: H pair-maximal, X a class, |X| >= l+2, l the budget parameter. Then (G,H,l) is YES iff (G-x,H,l) is YES for x in X.
(=>) restriction of a partition to G-x never costs more (the cited route through N3 is unnecessary). (<=) optimum of G-x has cost <= |H|+l,
|X-x| >= l+1 > l, N3 (valid in G-x: H still pair-maximal, X-x still a class) gives cluster(X-x) = (X-x) cup A_X, adding x as a twin costs 0.
Correct. Tightness: with |X| = l+1 it fails (deleting x can drop OPT): exhaustive runs show OPT(G-x) < OPT(G) in 5-11 instances per seed
exactly when OPT(G-x) - |H| >= |X| - 1. Tested: 2,500 exhaustive (seeds 21-24) and 4,360 ILP applicable cases, 0 failures (OPT(G-x) = OPT(G) always when
OPT(G-x) - |H| <= |X| - 2).
Gap 1 (statement): it needs H maximal, but the user-supplied H of CEaMP is not. Use Lemma 1: extend H to a pair-maximal H' (greedy,
polynomial), t = |H'| - |H|, new budget l' = l - t, answer NO if l' < 0, and use l' in the threshold l'+2. The tw(I(H')) may exceed tw(I(H)); so
the shrink rule is a preprocessing for (l, |H'|) or (l, tw(I(H'))), not for (l, tw(I(H))) with the given H.
Gap 2 (deleting clique components, A_X empty): CORRECT (X is a component of G; OPT unchanged; H untouched).

### Kernel-type bounds: CORRECT, wording
Classes with A_X non-empty are at most |V(H)| (disjoint non-empty subsets of V(H)), each of size <= l+1 after shrinking, so
|R| <= (l+1)|V(H)| and n <= (l+2)|V(H)| <= 3(l+2)|H|. Correct. It is not a kernel for the parameter l (|V(H)| is unbounded), and the
sentence "R is polynomially bounded by V(H) and l" is the accurate statement. Fix: do not call it a kernel.

### Pointer reduction (P1-P4): CORRECT, two remarks
Hypothesis for all four: instance is YES, optimum considered is any optimum (all are in normal form by Sec. 2).
* P1 (RH-del (X,a) is a unary flag at a, partner cl(a) unique): correct by N4 part 1. Remark (strengthening, not an error): deleting x-a also
  forces a to be cut from the rest of B_X cap A_X; a is adjacent to all of A_X except <= q packed partners (N4), so such an optimum pays
  >= |A_X cap B_X| - q further unpacked HH deletions. So an RH-del flag drags HH pointers with it; it is not "free unary".
* P2 (B_X cap A_X empty forces B_X empty): CORRECT by the cost count; used for every optimum. Computer: 12,000 exhaustive + 16,400 ILP
  class tests (maximise |B_X| over optima with z(x_0, A_X) = 0): 0 failures. (The authors' remark that a computer test of this is vacuous
  for isolated X is right; in the ILP it is not vacuous because B_X ranges over V(H) vertices.)
* P3 (RH-ins (X,a), a in B_X \ A_X: some z in A_X cap B_X; for every such z either a in Q(z) or the pair za is an inserted unpacked non-edge): CORRECT.
  The authors' sentence mixes "packed edge" with "packed non-edge": the dichotomy is "za packed (then a in Q(z), at most q classes
  per a) or za unpacked, then it is a non-edge by N4 part 3 and inserted since z,a share the cluster". Exhaustive 550 checks (seeds 21-24), 0 failures.
* P4 (|A_X| >= l+q implies B_X inside A_X): CORRECT, and the proof shows it for EVERY optimum, not "there is an optimum". Arithmetic rechecked:
  alpha >= q gives |A_X| <= l+q-1, alpha < q gives |A_X| <= l+q-2 (since delta <= l-1). Exhaustive 2,400 + ILP 780 classes, 0 failures.
  Needs q = max |Q(a)| over the vertices b in B \ A_X (the bound |Q(b)| <= q is what is used).
* Final summary sentence "every pointer touching R is unary, local, or reduces to an HH pointer": CORRECT as a classification of the
  pairs (X,a); it is a statement about the candidate partners, not about independence (HH pointers are also forced by P1).

## 4. Termination, ties, and clusters with two classes plus V(H) vertices (the specific worries)
* Ties: none exist. By the strict inequality of Sec. 2 no optimum has two classes in a cluster; nothing is to be chosen.
* Termination of the move: only relevant for non-optimal partitions, where the potential (number of class pairs sharing a cluster) still drops
  because a moved class becomes a singleton cluster; and every non-optimal partition of this type is strictly improved by one move.
* A cluster with X_1, X_2, A_1 cap C, A_2 cap C and other V(H) vertices: covered by the formula (the other vertices only enter |C|); tested above.
* Twin lemma applied before the exchange: needed so that the classes are whole; true for every optimum.

## 5. A real gap in the same note, outside the three claims: Sec. 2 (c = 1) uses C1 for a VERTEX-disjoint H
Sec. 2 applies 1.2-1.4 to a vertex-disjoint packing that is "reduced by BFK Rule 6.2". C1 needs pair-maximality (Sec. 2 above shows it fails
for vertex-maximal H; example n = 5). A vertex-maximal H may be extendable pair-wise: path P5, H = {(0,1,2)} is vertex-maximal and
pair-extendable by (2,3,4), so the pair-maximal H' has a vertex of degree 2 in I(H'), and the claimed "independent variables (w = 0)" CSP is
not available for H'. The conclusion survives with a fix: after BFK reductions |H| <= 2 l, extend to H' (|H'| <= |H| + l else NO),
so |V(H')| <= 9l and a Theorem-A style search over 7^{|H'|} assignments and O(l^2)-candidate U is FPT. But it is not "w = 0 + C1".
Same issue for Theorem K (components): C1 applied to the extension H' of a bounded-component packing needs I(H') components bounded.
Not tested by computer beyond the P5 and n = 5 examples (a3a_cex.py).

## 6. Table

| Statement | Verdict | Fix required before it goes in a paper |
|---|---|---|
| Twin lemma (adjacent true twins in one cluster, every optimum) | CORRECT | none; cite Guo 2009 as the known critical-clique fact, note proof is strict |
| Lemma 4 (R for pair-maximal H) | CORRECT | state hypothesis "no induced P3 with all pairs unpacked" (not vertex-maximal) and prove it in the C1 section |
| Exchange cost x_1(x_1+2a_1-\|C\|) | CORRECT | none (N(X_1) = X_1 + A_1 must be quoted from Lemma 4) |
| C1 strong form, "some optimum", iteration | CORRECT, weaker than true | replace by the strict version (x_1+x_2 <= 0): holds for EVERY optimum; delete iteration/termination; state pair-maximal H |
| C1 for H maximal only among vertex-disjoint packings | WRONG (not claimed explicitly, but used in Sec. 2) | n = 5 counterexample; keep the hypothesis visible |
| N1 types of unpacked edits | CORRECT | none |
| N2 "excess = ..." | GAP (wording) | it is \|U\|, not the excess; add the sum(\|x_i\|-1) term |
| N3 frozen classes \|X\| > l | CORRECT | state "YES instance, every optimum"; holds already for \|X\| > l_opt |
| N4 (a adjacent to one class; clique minus packed pairs; <= q) | CORRECT | none |
| N5 shrink rule \|X\| >= l+2 | CORRECT, threshold tight | add the extension step H -> H' with l' = l - t; w may grow |
| delete clique components | CORRECT | none |
| \|R\| <= (l+1)\|V(H)\|, n <= (l+2)\|V(H)\| | CORRECT | do not call it a kernel (not bounded in l alone) |
| P1 RH-del unary flag | CORRECT | remark that it forces >= \|A_X cap B_X\| - q HH deletions |
| P2 detachment never useful | CORRECT | none |
| P3 RH-ins local or HH | CORRECT | rewrite dichotomy (packed vs unpacked-non-edge) |
| P4 wide classes, \|A_X\| >= l+q | CORRECT, every optimum | say "every", define q on B \ A_X |
| "every R-pointer is unary / local / HH" | CORRECT | scope: classification of partners only |
| Sec. 2 (c = 1) uses C1 on vertex-disjoint H | GAP | use pair-maximal extension H'; |V(H')| = O(l); Theorem A search |

## 7. Reproduction
`cd code; python3 a3a_exh.py SEED 1500 7 9 [vertexmax|nonmax]` (seeds 11-14, 21-24 pair-maximal; 31 weakened hypotheses);
`python3 a3a_ilp.py SEED 100000 10 13 [vertexmax]` (600 s limit; seeds 101-103 pair-maximal n = 10..13, 104 vertex-maximal n = 10..12);
`python3 a3a_cex.py` (smallest counterexamples), `python3 a3a_n2.py 5 3000` (N2 identity, coverage statistics).
Limits: n <= 13; random-plus-planted families, with large R only through planted packings; l_opt <= 5; no instance has a class larger than 9;
the proofs, not the tests, carry the claims. The strict inequality of Sec. 2 is a pure algebra fact, verified on all partitions.
