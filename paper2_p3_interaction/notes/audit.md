# Audit of the theory notes (hostile referee, wave 2)

Scope: Obs 1, Obs 2, Thm A, Thm B, B1, B2, ETH claim (hardness A1), Lemmas 0-4 (fpt_attempt), C1, C3, C4 example.
Own code (not the authors'): `code/audit_lib.py` (bitmask graphs, exhaustive set partitions), `audit_structural.py`,
`audit_dp.py`, `audit_misc.py`, `audit_conj.py`. Verdicts: CORRECT / GAP / WRONG. Nothing was WRONG.
Convention used: CEaMP asks for S with |S| <= |H| + l (LPS lines 92-97), modification-disjoint = pair sets disjoint
= any two P3s share at most one vertex (the two are equivalent: any two vertices of a P3 form one of its 3 pairs;
LPS line 88-89 uses the "share at most one vertex" form). The code (`check_packing`, `exp_common.check_packing`)
asserts induced a-b-c with ab, bc edges, ac non-edge and pairwise-disjoint pair sets, i.e. the same definition.

## 0. Computer evidence (all passes, own code)
* audit_structural.py, 4 seeds x 1500 random instances, n = 5..9, graph families gnp / hub / noisy-cliques / paths,
  packings maximal / partial / "overlap" (greedy preference for P3s sharing vertices): 5565 instances,
  all partitions enumerated (Bell(9) = 21147). Checked: Lemma 0 on 394k solutions (cost <= |H|+3), Obs 1 equivalence
  (2647 YES / 2917 NO), Obs 2 on 238k triples, Lemma 2 + Lemma 3 on 561k (solution, unpacked edit) pairs, Lemma 4 on
  4265 R-pairs (1768 instances with |R|>=2), F1 and critical-clique property on all optima. 0 failures.
* audit_dp.py, 4 seeds, 888 instances (n 5..8, |H| <= 7), l = 0,1,2: my own implementation "enumerate U, then
  bucket elimination over min-fill order of I(H), domain {1..7}": value equals brute-force OPT whenever OPT <= |H|+l
  and rejects otherwise (2664 decisions, 590 NO). Asserted inside the DP that every triple scope is a clique of I(H);
  max factor scope never exceeded (width of order on I(H)) + 1 (0 violations, widths up to 5). Domain {0..7} gives the
  same value as {1..7} (232 runs). Thm B search tree vs exact min |U| given a random complete x: 5328 checks, 0
  mismatches, node count <= (3^(b+1)-1)/2 always.
* audit_misc.py: stars, 3x3 grid, C4 example, C1 identity (300), forced multi-edit search.

## 1. Statement by statement

### Obs 1 (l = 0 is a domain-3 CSP): CORRECT
Statement: if |S| <= |H| and S is a solution then S has exactly one edit in each pairs(P_i) and no other edit; G+S is
a cluster graph iff no triple induces a P3. Proof in analysis_01 is complete (>= 1 edit per P_i because the triple
{a,b,c} survives otherwise; disjointness; counting). Tested (5564 instances). Nothing to fix except: state it as
"solution of size at most |H|" (the budget is |H|+0, the proof uses <=), and note that a CSP has a variable per P3
and one constraint per triple; the number of constraints is Theta(n^3), not O(|H|).

### Obs 2 (scopes are cliques of I(H)) and the 3^{w+1} DP: CORRECT, one GAP in the DP statement
Proof: two pairs of one triple share a vertex v; owners P, P' both contain v, so P ~ P' in I(H). Correct, also when
the owners coincide (scope shrinks). Tested on 238k triples including hub-heavy overlapping packings.
GAP (a): "tw(I(H)) = w gives 3^{w+1} n^{O(1)}" presupposes that a width-w tree decomposition is given. Computing one is
not free: exact treewidth is 2^{O(w^3)} n (Bodlaender) or Korhonen's 2-approximation 2^{O(w)} n, which changes the
base (3^{2w+2}). State "given a tree decomposition of I(H) of width w" or put the approximation factor in.
GAP (b): the notes say "the constraint graph has the same cliques as I(H)". It is a SUBGRAPH of I(H) (not every edge
of I(H) is witnessed by a triple), so a decomposition of I(H) is a decomposition of the primal graph; say that.
GAP (c): the standard CSP-on-treewidth theorem is cited as UNVERIFIED in literature_2.md; give a reference
(Cygan et al. 2015 Ch. 7 DP / Freuder 1990) or one paragraph. Helly (a clique lies in one bag) is needed to place
the constraints; add it.

### Theorem A (7^{w+1} n^{2l+O(1)}): CORRECT
Restated: given (G,H,l) and a tree decomposition of I(H) of width w, one decides CEaMP in 7^{w+1} n^{2l+O(1)}.
Line-by-line: Lemma 0 gives every x_i non-empty, |U| <= l (since |S| = sum|x_i| + |U| >= |H| + |U|). Number of U:
sum_{j<=l} C(N,j) <= n^{2l} (N <= n^2/2 unpacked pairs). Flipping U leaves each P_i an induced P3 with the same pairs
(U avoids packed pairs), so H stays modification-disjoint and Obs 2 applies in G xor U. Constraint scopes: cliques of
I(H); unpacked pairs are constants after U is fixed; triples of empty scope are tested directly (n^3 per U).
Cost |x_i| unary. Accept iff min + |U| <= |H|+l. All steps hold; own implementation agrees with brute force.
Domain size 7 vs 8: 7 is right (non-empty subsets of 3 pairs; the empty value is never feasible because P_i itself
survives). 8 (hardness B1/B2) is merely a valid but looser bound; delete it and use 7 everywhere. Do not reduce
below 7: an edited packed P3 can have 2 AND 3 edits in every optimum:
 * forced 2 edits: n=5, edges 02,04,12,23,24, H={(1,2,3)}, OPT=2 (all optima put 2 edits in P);
 * forced 3 edits (all three pairs of P flipped in every optimum): n=6, edges 01,02,03,04,12,14,15,24,25,35,45,
   H={(0,3,5)}, OPT=3.
 In random data 8-9 % of instances have every optimum with a P3 receiving >= 2 edits, ~0.2 % with 3.
Remark to delete or fix: "this is no better than the known n^{2l} for c = 2" is a non-comparison; for c=2, w is
unbounded, Theorem A is simply a different tradeoff. Also the sentence in c2_algorithm 4A ("(l,w) gives XP in l, FPT in
w for every c") is right for Thm A but its exponent is n^{2l_a+2l_b} only if one follows LPS Alg 1; Thm A's
direct enumeration of U is simpler and subsumes it (the H_b guess is unnecessary; Lemma 0 handles multi-edited P3s
inside the DP).

### Theorem B (3^l n^{O(1)} given x): CORRECT, trivial
"Given a complete x, decide whether some U with |U| <= B (B = |H|+l - sum|x_i|) completes it." Search tree: an induced
P3 of the current graph must receive a flip; packed pairs are frozen by x; branch on its <= 3 unpacked pairs, depth <= B.
If all pairs of the P3 are packed, reject. Correct (completions only flip unpacked pairs; flipping a pair twice is
useless). Tested vs exact optimum, 5328 cases. Fix: say B = excess budget left after x; say the search is over
completions with exactly the packed edits x, not over all solutions; the "consequence" (CEaMP = find x with an
FPT-checkable property) is a reformulation, not a result. Do not present Thm B as progress toward FPT: finding x is
the whole problem.

### Hardness B1 (8^{w+1}, edits on packed pairs only): CORRECT, redundant
Same proof as Thm A with U = {}, domain 8 (7 suffices). As a lemma it only restricts S to packed pairs, which is not
a restriction on solutions of cost <= |H|+l once U is guessed. Merge into Thm A.

### Hardness B2 (n^{2l} 8^{w+1} XP, any c): CORRECT, duplicate of Thm A
Proof sketch in hardness.md is valid (|U| <= l by disjointness; flipping U keeps H a packing). Two defects: base 8
instead of 7; and "recovers the c=2 XP bound of Li et al. when w is bounded" is empty: for c=2, I(H) is subcubic but
may have unbounded treewidth (LPS "grid-like structures"), so the c=2 algorithm is not recovered, only matched on
the sub-family of bounded w. Remove the sentence.

### ETH lower bound 2^{o(w)} at l = 0: CORRECT (stronger than stated), proof needs a size table
Check against lps.txt (read, lines 249-870, not only the notes):
 * Reduction is from 3-CNF (line 251); preprocessing to "exactly 3 distinct variables per clause, each variable >= 2
   occurrences" is linear (pad short clauses with fresh variables / drop pure single-occurrence variables).
 * Variable x_i: 4 m_i cliques of 5 vertices (lines 308-310), 25 P3s (p,q free in F5) per even j, 2 m_i even j: 20 m_i
   vertices and 50 m_i P3s per variable; sum_i m_i = 3m.
 * Clause: 20 vertices (Q1,Q3,Q4,Q2 = 1+4+1+14) and 114 vertices (T_a,T_c,T_b = 34+34+46) (lines 817-829), 6 + 12
   transferring P3s, padding by Lemma 1 on cliques of constant size with p <= 23: at most |V||W| <= 20*46 triangles per
   T clique, so O(1) padding P3s per clause. Totals |V| ~ 194 m, |H| = O(m) with an absolute constant (<= ~3000 m by
   crude counting; I did not compute the exact constant, none is needed). Polynomial: stated in LPS line 868-870.
 * So m = O(N) after the Sparsification Lemma (Impagliazzo-Paturi-Zane) gives |H| = O(N), and tw(I(H)) <= |H| - 1 is
   trivial (|H| nodes), hence w <= |H| - 1 = O(N). An algorithm in 2^{o(w)} n^{O(1)} would be 2^{o(N)} for 3-SAT.
   This argument works verbatim for an algorithm that computes its own decomposition (2^{o(tw)} <= 2^{o(|H|)}), so no
   decomposition needs to be given: state the theorem as "no 2^{o(tw(I(H))) n^{O(1)}}".
Caveats to write down: (1) it shows only that the exponent is Theta(w) up to constants in the exponent, not that the
base 3 (or 7) is tight; (2) it holds even though I(H) has maximum degree <= 66 and c <= 23; (3) cite Cygan et al.
Thm 14.4 style (ETH + sparsification) and LPS Thm 1 precisely (version v3), list the size table above; (4) the notes
say "not re-derived by me": I re-derived the structure (every gadget is O(1) per clause or occurrence) but not each
numerical constant of the P3 count.

### Lemma 0 (decomposition): CORRECT
Tested on 394k solutions. Needs only: S arbitrary set of pairs, G xor S cluster.

### Lemma 1 (WLOG H maximal): CORRECT but trivial; caution holds
"OPT <= |H|+l iff OPT <= |H'| + (l-t)" is the identity |H'|+l-t = |H|+l. Valid when H' extends H by t further pair-disjoint
induced P3s; if l - t < 0 the instance is a NO instance (OPT >= |H'|), consistent. Fix: say so; the warning that
tw(I(H')) may exceed tw(I(H)) is correct and important.

### Lemma 2 (viability of an unpacked edit): CORRECT
Restated: S any solution, U = unpacked part, uv in U, D(uv) = {z notin Q(u) u Q(v) u {u,v} : {u,v,z} induced P3
in G xor uv}. Then each z in D has uz or vz in U and |D| <= |U|-1. Proof checked: all three pairs unpacked, so if uz, vz
not in S the triple would stay P3; distinct z give distinct pairs (uz = vz' impossible). Tested on 561k pairs, 0 failures.
Fix: |Q(u)| = 2 (number of P3s through u) <= 2(w+1) is correct (disjoint pair sets) and the Helly step is correct; say Q(u)
is defined by "uz packed". The note "all pendant edges viable" is right, this lemma gives no candidate bound.

### Lemma 3 (locality of triples): CORRECT (a one-line observation). 
### Lemma 4 (structure of R for maximal H) and corollaries: CORRECT
Edge case rs: z with exactly one of rz, sz an edge gives an all-unpacked induced P3. Non-edge: common neighbour z may
lie in V(H) but rz, sz unpacked since r, s in R; correct. (i), (ii) correct (shortest-path argument). (iii) correct:
a permutation of a twin class is an automorphism of G fixing V(H) pointwise, so U can be moved to the first 2l class
members. (iv) = Guo's critical clique property (some optimum keeps each critical clique in one cluster; the class X
lies inside one critical clique): I tested it on all optima, it holds in every instance (all 5565); cite Guo, TCS 410
(2009) after checking the exact statement. F1: strictly stronger holds: in EVERY optimum each cluster is connected in G
(splitting saves |C1||C2| > 0), verified on all optima; the "clique component is a cluster of every optimum" corollary is
correct. Tested Lemma 4 on 4265 pairs.

### C1 of hardness (w alone / w = 0 is NP-hard): CORRECT
Take a maximal vertex-disjoint packing (greedy, polynomial), l = k - |H|; the instance asks for |S| <= k exactly.
Checked on 300 random cases. Fix: say "unbounded l"; avoid "H empty": w of the empty graph is -1 or 0 depending on
convention and H empty is a legitimate but degenerate packing; the maximal packing version is the clean statement.
### C3 (tw(G) and tw(I(H)) incomparable): CORRECT, overstated remark
(a) star K_{1,2m}: pairs disjoint, I(H) = K_m, tw(G) = 1, OPT = 2m-1 verified for m = 1..4 (l = m-1). (b) grid with
vertex-disjoint row P3s: I(H) edgeless, tw(G) = k (grid3 packing checked; tw(k x k grid) = k is standard).
Fix: the I(H) is a function of H; the statement is "for some admissible H". The consequence "hubs force large
tw(I(H))" is only true as a tradeoff: any packing on the star has all P3 through the hub, so tw(I(H)) = |H|-1 and
l = 2m-1-|H|; small tw forces large l. Write it as l + tw(I(H)) >= 2m-2.
### C4 example (fpt_attempt route (a)): CORRECT
G = C4 (01, 02, 13, 23), H = {(0,1,3)}, OPT = 2. Exactly three optima {03,12}, {01,23}, {02,13}; each has an unpacked pair
at vertex 2; H maximal; l=0 CSP unsatisfiable (all verified exhaustively). It refutes the strong locality claim only.
Reproduced statistic: among maximal-H instances with |R| >= 2 (n <= 9, my generator, 971 instances) 37 have every
optimum touching R by an unpacked edit; 0 have a forced R-R unpacked edit (supports, does not prove, Conjecture C1;
my sample is not the authors' and n <= 9 is tiny).

## 2. Table

| Statement | Verdict | Must be fixed before going in a paper |
|---|---|---|
| Obs 1 | CORRECT | say "|S| <= |H|"; constraints are Theta(n^3) |
| Obs 2 | CORRECT | primal graph is a subgraph of I(H) (not "same cliques"); Helly for bag placement |
| 3^{w+1} DP at l=0 | GAP | tree decomposition must be given (or pay 2^{O(w)} approximation, base changes); cite CSP-on-treewidth theorem |
| Thm A: 7^{w+1} n^{2l+O(1)} | CORRECT | same decomposition caveat; use 7 everywhere (explain non-empty); show the 2- and 3-edit examples; drop "no better than n^{2l}" remark |
| Thm B: 3^l given x | CORRECT (trivial) | define budget B and "completion"; do not oversell as FPT progress |
| B1: 8^{w+1} packed edits only | CORRECT, redundant | merge into Thm A with 7 |
| B2: n^{2l} 8^{w+1} | CORRECT, duplicate | delete or restate with 7; remove "recovers c=2 XP" |
| ETH 2^{o(w)}, l = 0 | CORRECT (stronger: 2^{o(tw)}) | add size table (|V| ~ 194m, |H| = O(m)), cite sparsification, state it is exponent-only (no base), cite LPS v3 Thm 1 |
| Lemma 0 | CORRECT | none |
| Lemma 1 | CORRECT, trivial | note l-t<0 means NO; keep the tw warning |
| Lemma 2 | CORRECT | none (state S arbitrary) |
| Lemma 3 | CORRECT, trivial | fold into Lemma 2 |
| Lemma 4 + (i)-(iii), F1 | CORRECT | none; cite Guo 2009 for (iv) after checking text |
| C1 (w alone) | CORRECT | use maximal vertex-disjoint packing, say l unbounded |
| C3 (incomparable) | CORRECT | "for some H"; replace "hubs force large tw" by l + tw >= 2m-2 |
| C4 counterexample | CORRECT | scope: refutes only "all unpacked edits inside V(H)" |
| "modification-disjoint" in code | CONSISTENT | pair-disjointness == share <= 1 vertex; keep both phrasings in the paper |

## 3. Not audited / residual risk
* hardness A2/A3 (constraint-language closure), c2_algorithm sections 1, 3, 4C ("components of I(H)" BFK trick) and
  fpt_attempt route (b) (21-way branching for w = 0, marked sketch by the authors): not tested by me; (b) looks right
  for w = 0 (scopes unary after U is fixed) and, if written up, would give an independent FPT proof for c = 1.
* All tests are exhaustive only for n <= 9; none probes large tw, large l or l > 2 in the DP. The proofs, not the tests,
  carry the claims; the tests catch definition mismatches (they found none).
* Conjecture C1 and the W[1]/FPT status of (l, w) remain open; nothing here supports either direction.
