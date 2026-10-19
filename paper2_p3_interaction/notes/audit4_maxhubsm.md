# Audit 4 (hostile referee, wave 4): MAX, HE, HUB, SM, soft-graph treewidth

Target: notes/wave2_hardness.md (F0-F2, Sec. 2 A/B/C, Sec. 3). Own code written from scratch: code/a4_lib.py (library), a4_exh.py
(exhaustive), a4_ilp.py (HiGHS), a4_adv.py (adversarial families), a4_hub2.py (strengthened HUB), a4_mut.py (mutation test).
None of the authors' or earlier auditors' scripts is imported. Verdicts: CORRECT / GAP / WRONG. Nothing is WRONG.

## 0. Conventions (what "hard", "soft", "maximal" must mean)
* packed pair (= soft) = two vertices lying in a common member of H (all 3 pairs of every member, edge or non-edge).
  hard pair = every other pair. hard edge = hard pair in E.  Two members of H share <= 1 vertex, so the packed partner sets
  of a vertex over different members are disjoint: D (max number of packed partners) equals 2c exactly (c = max members through a
  vertex), not only D <= 2c.  c <= w+1 (members through one vertex form a clique of I(H)).
* pair-maximal H: no induced P3 of G has all three pairs hard.  (Not vertex-maximality; see audit2_c1 Sec. 0.)
* YES instance with budget l: some S with |S| <= |H|+l.  For such S: U = S minus packed pairs has |U| <= |S|-|H| <= l
  (every member of H needs >= 1 edit among its own 3 pairs, members are pair-disjoint).
* Computer scope. Exhaustive (all set partitions, n = 6..9): 18,800 instances (gnp, noisy clusters, hub, sparse, Fano / AG(2,3)
  soft-clique gadgets, +noise; H maximal / partial / overlap-greedy), ~931,000 solutions of cost <= |H|+3 checked for SM, YES
  instances for l = 0..3 (about 51,000 instance-l pairs) for HUB.  ILP (n = 10..13): 1,600 instances, 6,077 solutions
  (optimal, random-objective with cost <= |H|+l, and 1,745 "forced multi-piece" ones: a, b in different components of the hard graph,
  forced in one cluster, cluster of a maximised, l = 0..3).  Mutation test (a4_mut.py): each weakened version (2l instead of 2l+1,
  Ds-1, 3(w+1)-2, ...) is falsified, so the checks are not vacuous.  0 violations anywhere.

## 1. MAX  -- CORRECT (two GAPs of wording; treewidth claim right and tight)
Restated. Let H be ANY packing of induced P3s (no maximality), M a maximal set of induced P3s of G whose pairs are pairwise
disjoint and disjoint from pairs(H).  Then (a) OPT >= |H|+|M|; (b) OPT <= |H|+l iff (|M| <= l and OPT <= |H u M| + (l-|M|));
(c) H u M is a packing and is pair-maximal; (d) tw(I(H u M)) <= tw(I(H)) + |M|.
Line by line.
* (a): H u M is a pair-disjoint packing of induced P3s; each needs its own edit.  OK.  (b): |H u M| + l - |M| = |H|+l, so (b) is an
  arithmetic identity once (a) gives |M| <= l for YES.  So the lemma's whole content is (a) + (c) + (d).  Say so.
* (c) is NOT proved in the note, it is used ("WLOG no induced P3 has all three pairs hard").  Proof: an induced P3 with all pairs hard
  w.r.t. H u M has pairs disjoint from pairs(H) and from pairs(M), so could be added to M, contradicting maximality.  Tested: 
  H u M pair-maximal in every run.  GAP 1: write (c); note "hard" afterwards means hard w.r.t. H u M, not w.r.t. H (pairs of M became soft).
* (d): put the |M| new nodes into every bag.  Every edge touching a new node is covered (the other endpoint lies in some bag, which
  also holds the new node), new nodes have connected (whole-tree) subtrees.  Width +|M|.  The worry "a new P3 is adjacent to many bags"
  is irrelevant: the bound is additive and independent of adjacency.  The bound is TIGHT in general (code/a4_adv.py star): H = {(p,v,q)}
  in a star with centre v, M = m arms (a_i,v,b_i): tw(I(H)) = 0, I(H u M) = K_{m+1}, tw = m = 0+|M| (m = 1..5 verified).  Also tw
  rises by Theta(|M|) from 0 for M-nodes arranged as a cubic graph with H-nodes subdividing its edges (K4 example: 0 -> 3).  Random runs:
  tw(I(H u M)) - tw(I(H)) <= |M| in 58,000 (H, M) pairs, equality in ~21%.  No better general bound exists; do not hunt for one.
* GAP 2 (bookkeeping the note omits): after MAX the instance is (G, H u M, l-|M|) with width <= w+|M|, c' <= c+|M|, D' <= D+2|M|.
  HUB/SM/HE must then be applied with the primed quantities (and l' = l-|M|).  Also a tree decomposition of I(H u M) is obtained in
  polynomial time from one of I(H), so MAX costs nothing algorithmically; M is found greedily in O(n^3) per pass.
Verdict CORRECT.  Fix: add proof of (c), state hypotheses (any packing H), primed parameters, tightness remark.

## 2. HE  -- CORRECT, but only for pair-maximal H (the note states it for "WLOG" after MAX; say H u M)
Statement (precise): H pair-maximal, ab and bc hard edges (in E, not packed).  Then ac in E or {a,c} packed.
Proof: otherwise a-b-c is an induced P3; its pairs ab, bc are hard by assumption; if ac is hard too, all three are hard, contradicting
pair-maximality; so ac is packed.  No YES hypothesis needed.  Also needs a,b,c distinct (trivial).  Without maximality it FAILS
(H = empty, any induced P3): the exhaustive run found 10,780 failing configurations for non-maximal H, and 0 failures over ~198,000
configurations (a,b,c) for pair-maximal H u M.  Fix: state "H pair-maximal (e.g. H u M from MAX)" in the statement itself.

## 3. HUB  -- CORRECT as stated, but it is the wrong (weak) form for the use made of it; here is the right one
Restated. YES instance (OPT <= |H|+l), H any packing, N_h(u) = {x : ux in E, {u,x} not packed}.
(a) X subset of N_h(u) with every pair a hard non-edge  =>  |X| <= 2l+1.   (b) alpha(G[N_h(u)]) <= (D+1)(2l+1).
Line by line.  (a): (x1,u,x2),(x3,u,x4),... are induced (ux_i in E, x_ix_j non-edge); all pairs hard; pair-disjoint; so a packing of
size floor(|X|/2) extending H; OPT >= |H|+floor(|X|/2) => floor(|X|/2) <= l.  OK.  Hypotheses: YES instance; H need not be maximal.
(b): an independent set I of G[N_h(u)] of size m; W[I] (packed pairs inside I) has max degree <= D, greedy gives X subset I,
|X| >= m/(D+1), pairwise non-adjacent and non-packed, i.e. all pairs hard non-edges, so |X| <= 2l+1.  OK.
GAP A (sharper, free): inside an independent set I each member P of H contributes at most ONE packed pair (a P3 with two vertices of
I has them as its ends, centre adjacent to both hence not in I), so the packed degree inside I is <= c = D/2:
alpha(G[N_h(u)]) <= (c+1)(2l+1) = (D/2+1)(2l+1).  Equality observed at l = 0 (ratio 1.0 in thousands of instances, c = 1, alpha = 2).
GAP B (important, narrative): Sec. 3 E2 says the deletion hub "needs large independent sets inside N_h(u)".  Vertex Cover on the
non-adjacency graph Gamma of N(u) is hard through large MATCHINGS of Gamma, not through large independent sets of G (= cliques of
Gamma).  Example: N_h(u) = two disjoint cliques of size m: alpha = 2, yet (a_i,u,b_i) are m pair-disjoint all-hard P3s, so m <= l.
So alpha-bound does not address E2.  The correct lemma (PROVED here, VERIFIED-BY-COMPUTER):
 HUB+ . YES instance, any H.  Let Gamma_u = graph on N_h(u) whose edges are the HARD NON-ADJACENT pairs.  Then nu(Gamma_u) <= l
 (matching x_iy_i gives pair-disjoint all-hard P3s (x_i,u,y_i)), so Gamma_u has a vertex cover T, |T| <= 2l; N_h(u) \ T has every pair
 adjacent or packed; consequently alpha(G[N_h(u)]) <= 2l + c + 1 (additive, much better than (c+1)(2l+1)).
Verified (a4_hub2.py): nu <= l and alpha <= 2l+c+1 on 64,347 (vertex, l) checks (n = 6..9) + 2,755 ILP ones (n = 10..12), 0 violations.
Interpretation: a deletion hub can have only O(l) "deletion partners" beyond a soft-clique remainder; the hub encodings E2 are dead
for this reason, and the proof is two lines.
Verdict CORRECT.  Fix before inclusion: state YES hypothesis and which packing (H u M: D', l'); replace (b) by HUB+ ; define N_h.

## 4. SM (small merge)  -- (i), (ii), (iii) CORRECT; (iii) loose, a sharp version proved and tested; narrative GAPs
Setting (precise).  S any solution (cluster partition) with |S| <= |H|+l; U = unpacked edits, u := |U| <= l; F = deleted hard edges, I =
inserted hard pairs (I = U \ E); kept hard edges = hard edges of G not in F; PIECES = connected components of the graph (V, kept hard
edges) (every vertex lies in exactly one, isolated vertices are singleton pieces); every piece lies in one cluster.  Ds = actual max
packed degree = D.  NOTE: pieces depend on S, not on (G,H); the lemma is a statement about each solution, not an instance invariant.
No maximality of H and no optimality of S is used (confirmed: tested on all solutions of cost <= |H|+3, also for non-maximal H).
(i) CORRECT.  a in A, b in B, different pieces, same cluster => ab is an edge of G xor S.  If ab hard in E and kept, a,b same piece;
 hard in E and deleted, ab is not an edge of G xor S; so hard non-edge of G, inserted (ab in I).  Otherwise ab is packed.
(ii) CORRECT.  Soft pairs between A and B <= Ds*min (count from the smaller side), all other pairs of A x B are in I:
 min(|A|,|B|)(max-Ds) <= |I cap (A x B)| <= u.  If max >= Ds+u+1 then LHS >= u+1 (min >= 1): so max <= Ds+u <= D+l.
 Statement applies to pieces that SHARE a cluster with another piece; a lone piece can be arbitrarily large (it is the whole cluster).
 Sec. 3's phrase "pieces are blobs of size <= D+l" must carry this qualifier.
(iii) CORRECT but loose.  Note's bound max(2D+l, 2l): algebra s(s-M-Ds) <= 2u with M = Ds+u is right; each piece in C has size <= M by (ii).
 SHARP version (proved here; 0 violations in the whole computer scope incl. forced multi-piece ILP clusters up to size 13):
 a cluster with >= 2 pieces has |C| <= max(D+u+1, 2D, 2u).   Proof: let A be the largest piece, a = |A|, r = |C|-a >= 1.
 Case a > D: every x in C\A has <= D soft partners so >= a-D inserted pairs to A, hence r(a-D) <= u, r <= u, a <= D+u/r,
 |C| <= D + u/r + r <= D+u+1 (convex in r on [1,u]; case impossible when u = 0).  Case a <= D: every x has >= |C|-D pairs outside its piece,
 <= D of them soft, so >= |C|-2D inserted pairs, and |I| >= |C|(|C|-2D)/2; if |C| > 2D then |C| <= 2u.
 Observed: max multi-piece cluster equals the vertex count (9 in AG(2,3) soft clique, D = 8, singleton pieces, u = 0: size D+1, valid
 because 2D = 16); hill climb (a4_adv.py pieces) came no closer than |C| = sharp bound - 3.  I do not claim sharp is tight.
GAP in Sec. 3 (not a lemma defect): E1 "adjacency by a HARD edge puts them in the same piece" ignores the <= l deleted hard edges; "the
 compatibility graph realisable by hard edges has degree <= D+l" is not what (ii) says (it bounds sizes of merged pieces, not degrees).
 The collapse discussion is a heuristic explanation, not a proof of FPT; keep it out of theorem statements.
Verdict CORRECT.  Fix: state S-dependence of pieces, the qualifier "sharing a cluster", replace (iii) by the sharp bound (or keep, with remark).

## 5. Soft-graph treewidth  -- CORRECT (edge case, and the bound is not tight beyond w = 0)
Statement: H nonempty, W = (V, packed pairs).  tw(W) <= 3(w+1)-1 with w = tw(I(H)).
Proof check.  Replace each bag by the union of the vertex sets of its members (<= 3(w+1) vertices).  Edge cover: a packed pair lies in a
member which lies in a bag.  Connectivity for vertex v: members through v are pairwise adjacent in I(H), so their subtrees pairwise
intersect, and a union of pairwise-intersecting subtrees is connected (Helly is not needed, pairwise intersection suffices).
Vertices outside V(H) are isolated in W: add as one-vertex bags.  OK.
GAP (edge case): H = empty gives W edgeless, tw 0, and the formula gives -1; write max(0, .).  Tight for w = 0 (single P3: tw(W) = 2);
for w >= 1 the largest tw(W) observed was w+1 (w = 1: 2, 2: 3, 3: 4, ..., 9: 7; n <= 12), far below 3w+2; a sharper bound may exist,
not needed for FPT purposes.  Computer: exact tw by subset DP on 18,800 + 1,600 instances, 0 violations, mutation 3(w+1)-2 falsified.
Related facts used on the way: D = 2c <= 2(w+1) CORRECT; clique number of W <= 3w+3 correct (tw+1).

## 6. Final table
| statement | verdict | fix required before inclusion | state as theorem? |
|---|---|---|---|
| MAX: OPT>=|H|+|M|; equivalence of budgets; tw(I(H u M)) <= w+|M| | CORRECT (tight) | prove pair-maximality of H u M; primed parameters (w+|M|, c+|M|, D+2|M|, l-|M|); say (b) is arithmetic | YES, as Lemma (preprocessing lemma) |
| HE: hard ab, bc => ac in E or packed | CORRECT only for pair-maximal H | put hypothesis in statement; "hard" w.r.t. H u M | Lemma / Observation (trivial) |
| HUB (a): hard-independent X has |X| <= 2l+1 | CORRECT | YES hypothesis; which packing | YES, but as a corollary of HUB+ |
| HUB (b): alpha <= (D+1)(2l+1) | CORRECT, weak and misdirected | replace by HUB+ : nu(Gamma_u) <= l, VC <= 2l, alpha <= 2l+c+1 | YES (HUB+, not (b)) |
| SM (i) | CORRECT | define pieces from S; no maximality needed | YES, Lemma |
| SM (ii) | CORRECT | qualifier "piece sharing a cluster with another piece" | YES, Lemma |
| SM (iii) |C| <= max(2D+l, 2l) | CORRECT, loose | optionally sharp max(D+u+1, 2D, 2u) (proved here) | YES (either form) |
| Sec. 3 collapse of E1-E4 | NOT a proof (heuristic; E1 ignores deleted hard edges; E2 argued with alpha instead of matching) | rewrite E2 with HUB+; mark E1 as heuristic | NO |
| tw(W) <= 3(w+1)-1 (and D = 2c <= 2(w+1)) | CORRECT | max(0,.) for H empty; pairwise intersection instead of Helly | YES, Lemma |

## 7. Recommendation
Strong enough for theorem/lemma status as they stand or after the cosmetic fixes: MAX, HE, SM (i)-(iii), soft-graph treewidth, HUB (a).
Replace HUB (b) by HUB+ (matching form, additive bound): it is stronger, equally short, and actually addresses the deletion-hub
encoding.  All five items use only: pair-disjointness, a YES budget, and counting; none needs optimality of S.  What these lemmas do
NOT give: FPT in (l, w) or any hardness.  The claim "the problem is FPT in (l, w)" (wave2 Sec. 0 line 1, Sec. 4) remains a belief
resting on an unproved locality statement for tiny-tiny pointers; do not put it in a paper beyond an open question, and do not call
Sec. 3 a proof that hub/tile/slip encodings fail (it explains why the natural gadgets die, nothing more).
