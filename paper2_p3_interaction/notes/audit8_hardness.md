# Audit 8 (hostile referee): notes/wave6_hardness.md

Scope: every item labelled PROVED in `wave6_hardness.md` (Lemmas 1-8, 7b, Prop 9, the "inactive anchors cost nothing" statement, "tunnel fan-out <= D+l", the conditional FPT via d-Hitting Set / Baker-Pixley)
plus the computer claims that the paper section will quote (flip gadget G*, cost sharing, small cores, Conjecture C2).
Method: each proof re-derived line by line from the definitions and from `paper/sections/{prelim,l0,packing,dproute}.tex`; independent brute force with a fresh library
`code/a8_lib.py` (all set partitions, n <= 10, numpy only; nothing imported from `w6c_*`, `a7_*`). No git; paper untouched by this audit step.
Instances: planted cluster graphs with 1-7 random edits, `G(n,p)`, with `n = 6..10`, a random greedy pair-disjoint packing (hence always pair-maximal); a fraction with one big cluster (large hard degree).

## Verdicts

| claim | verdict | one-line reason |
|---|---|---|
| Lemma 1 (zero-excess structure) | CORRECT | proof right; (c) also gives that a hard component is a hard clique up to soft pairs |
| Lemma 2 (clean tunnel = lock) | CORRECT WITH FIX | proof right; the parenthetical "all pairs at x, x' unedited" must read "all *hard* pairs"; with all pairs unedited the stated consequence (a member through x is edited at x) is contradictory |
| "inactive anchors cost nothing" | CORRECT WITH FIX | true exactly as: `OPT(G + stub on clean x) = |H|` iff some ZE solution of G has x as a singleton (two-ended tunnel: both ends simultaneously). It is not "always free": about 60% of the clean anchors tested are not singletons in any ZE solution |
| Lemma 3 (non-respected tunnels cost >= 1/2 each) | CORRECT | (not in the requested list; used by the sharing remark) |
| Lemma 4 (release closes the cluster) | CORRECT WITH FIX | first sentence (`cluster(x) subset {x} u Q(x)`) right. Second sentence ("P is in a cost-1 state and x is together with y") is **wrong** as stated: cost 2 is possible, and the cost-1 state L has x a singleton. Needs "if P has exactly one edit and x is not a singleton" |
| Lemma 5 (forcing, EQ) | CORRECT WITH FIX | transitivity argument right; "form a triangle of I(H)" should be "form a clique of I(H) (a triangle if the three owners are distinct)" |
| Lemma 6 (occupancy) | CORRECT | |
| Lemma 7 (splice) (a), (b) | CORRECT | the "Reading" (a persistent bit that flips across the released anchor) is an interpretation, not a statement: HEURISTIC |
| Lemma 7b (majority => cores <= 2) | CORRECT | proof right; the Baker-Pixley step can be proved in 5 lines (below). But its *use* in wave6 Sec. 4.7 rests on Conjecture C2, which is FALSE as stated |
| Lemma 8(a) (released tunnel ends have hard degree <= D+l; big hard components frozen) | CORRECT WITH FIX | the numerical statements are right (Lemma EX(a), NF(gamma)); "frozen: no pair of it is edited" is false for *soft* pairs (a soft non-edge inside a big hard component can be inserted) - read "no hard pair edited". The sentence "hence the relations realised by releasable tunnels have bounded degree (partial matchings, cliques)" is a heuristic gloss |
| Lemma 8(b) (bounded-degree Grid Tiling in n d^{2k} poly) | CORRECT | |
| Prop 9 (TR model: cores <= d => d^l 3^{w+1} n^{O(1)}) | CORRECT WITH FIX | (i)-(iii) and the algorithm are right, but the model must be *defined* as: every member gets exactly one edit and the hard edits are exactly cuts of tunnels. Otherwise the minimum hitting set is only an upper bound of OPT-|H| (strictly larger in 44% of the tunnelled instances tested). Also: FPT in the model is not FPT of CEaMP |
| Contrapositive "TR-hardness needs a PFG" | HEURISTIC | Lemma 7 gives only one direction of memory (disjoint boundary sets); no equivalence is proved |
| G* is a 7-vertex tight flip gadget; mini-track | CORRECT (computer, re-verified independently) | see below |
| cost-sharing example (Sec. 4.1) | CORRECT (computer, re-verified) | OPT 2 / 4 / 3 / 2 |
| no minimal non-face of size >= 3 in tight graphs | CORRECT as evidence (re-reproduced) | atlas n <= 7, random n = 8, 9: max core size 2 |
| negative composition results (Sec. 4.5-4.6) | NOT RE-RUN | quoted as evidence of the notes only |
| "Sigma closed under M-majority for all c <= 2" / Conjecture C2 | **WRONG** (the claim); C2 as stated is FALSE | explicit counterexample with n = 6 (below); 0.4-1.4% of random c <= 2 tight instances fail; the weaker statement "Sigma is 2-decomposable for c <= 2" survives (0 failures) but is unproved |
| G* is not 2-decomposable (RRRMM) | CORRECT | re-verified |

## Details and tests

### Lemma 1 CORRECT
(a) hard pairs are unedited by ex = 0. (b) `t(ab),t(bc),t(ac)` over cluster-graph patterns {000,100,010,111} (001 and the others are not transitive), costs 2,1,1,1; the cost-1 patterns are R, L, M.
(c) a hard pair inside a cluster is an unedited edge, so a hard non-edge is never inside a cluster: each hard component lies in one cluster and pairs between different components inside a cluster are soft.
`a8_l123456.py`: cost = sum c_i + |U| and c_i >= 1 for every partition (16.4M partitions, 2800 instances n = 6..9, all pair-maximal); (a),(b),(c) on all 2832 ZE partitions. 0 violations.

### Lemma 2 and inactive anchors CORRECT WITH FIX
Proof: t in cluster(x), t != x, x'. If t not in Q(x) then xt is hard, and a non-edge since N_h(x) = {x'}; inside a cluster it would be edited. So t in Q(x), likewise t in Q(x'), contradiction with Q(x) n Q(x') = empty.
The consequence "every member through x has x at an end, edge deleted, mates together" needs ZE (uses c_i = 1); with the literal hypothesis "all pairs at x unedited" it is contradictory (a member through x has both pairs at x apart, so x's soft pairs are edited). Fix: hypothesis "all hard pairs at x and x' unedited".
Inactive anchors: (=>) split {x,x'} into singletons (no other pair changes); (<=) put the stub next to a singleton x. Tested (`a8_l123456.py`): stubs on 2743 clean anchors: `OPT = |H|` iff x is a singleton in some ZE solution of G (1104 yes, 1639 no, 0 exceptions); 1278 ZE solutions with stub all have cluster {x,x'}, no centre, mates together; paired tunnels (404 pairs: 84 yes, 320 no): same equivalence, 0 exceptions.
So "an inactive anchor costs nothing" is a statement about anchors that are singletons in some ZE solution; the paper must state it that way.

### Lemma 3 CORRECT
Proof checked: a non-respected tunnel is cut or has a third vertex t in its cluster, and then xt or x't is a hard non-edge in a cluster; each pair is charged by at most two vertex-disjoint tunnels. 18,389 partitions with cost <= |H|+2 on instances with 1-2 stubs: 0 violations.

### Lemma 4 CORRECT WITH FIX
`cluster(x) subset {x} u Q(x)` is right (7022 checks). The sentence "if P is the only member through x, P is in a cost-1 state and x is together with y" is false: in 573 tested situations c_i = 2 (cost 2: x alone and y,z apart or together at cost 2), and in 1178 situations P has c_i = 1 and x is a **singleton** (state L: x deleted from y, which is exactly what the cut of xx' makes useless). Correct: cluster(x) subset V(P); if c_i = 1 and x is not a singleton then x is together with at least one other vertex of P (cluster {x,b} or {x,b,c} when x is an end).

### Lemma 5 CORRECT WITH FIX, Lemma 6 CORRECT
Transitivity. The three owners of yz, qy, qz pairwise share a vertex, so they form a clique of I(H) (possibly fewer than three distinct members). 59,896 checks of `t(q,y) = t(q,z)` and of the forced value of a soft pair next to a hard one; 113,139 checks of Lemma 6. 0 violations. (The "Consequence" paragraph about exactly one admissible pair of neighbour states refers to the 7-vertex coupler search; not re-run.)

### Lemma 7 CORRECT (reading HEURISTIC)
(a) `Sigma` = state vectors of ZE solutions; ZE solution <-> vector (U empty, one edit per member). The ZE condition is "no vertex triple is an induced P3 of G + S", a conjunction of constraints whose scope (owners of the three pairs, at most 3 members) is a clique of I(H) (Observation obs:two), hence lies in one bag, hence in M1 u B or in M2 u B; B-agreement then gives the splice.
(b) Singleton-ness of an anchor is a function of the states of the members through it (shown below), so the splice of sigma^a on M1 u B with sigma^b on M2 is a ZE solution that makes all of S singletons, contradicting non-face-ness. Right, provided the members through S_i lie in M_i (hypothesis stated).
Tests (`a8_l7.py`, 4 x 2500 instances, n = 7..10, tree decompositions by min-fill, all tree edges): 34,262 splice checks (a) with no failure; 44 non-trivial instances of (b) (minimal non-face of size 2 with a split and B non-empty): R1 n R2 = empty each time. (b) is rarely applicable at this size; the proof is short and I re-derived it.
Whether the "persistent bit" language is justified: Lemma 7 only gives that boundary-state sets of the N sides are pairwise disjoint; it does not say that a single bit is carried. HEURISTIC.

### Lemma 7b CORRECT
Step 1 (anchors force states; face iff partial vector extends): verified, 5271 anchor subsets in 500 instances, 0 mismatches (`a8_bp.py`). Also needs the observation that an anchor with a centre role is never a face.
Step 2 (Baker-Pixley): if R' = proj_T(Sigma) is closed under a majority operation m and every binary subtuple of t is in the binary projection, then t in R'. Self-contained proof by induction on |T|: for |T| >= 3 take u_j in R' extending t restricted to T minus {j} for j = 1,2,3 (induction); m(u_1,u_2,u_3) = t coordinatewise. (Tested abstractly: 2529 random relations closed under random uniform majority rules, 612,712 non-member tuples: the criterion never accepts a non-member.)
Conclusion right: every minimal non-face has at most two anchors. Caveat: the hypothesis is closure under a *uniform* majority operation. Empirically no minimal non-face of size >= 3 occurs at all in the tested instances, so the conclusion is not tested on real instances.

### Lemma 8 (a) CORRECT WITH FIX, (b) CORRECT
(a) With H pair-maximal and S a solution with u edited hard pairs: Lemma EX(a) gives d_h(x) <= D + u for each endpoint x of an edited hard pair, and u <= |S| - |H| <= l; Lemma NF(gamma) (S optimal, |S| <= |H|+l) gives that a hard component containing a vertex of hard degree > B = 2D + 3l is a whole cluster. Tests (`a8_l8.py`, 3 x 800 instances n = 7..10 with big clusters, l = 0,1,2): 29,129 edited-hard-pair endpoints, 814 big components: no violation. In 716 of those 814 the cluster contains an *edited soft pair*: the paper's wording "no pair of it is edited" is false. Minimal example: K_9 minus the edge ac, H = {(a,b,c)}: valid, pair-maximal, OPT = 1 (insert ac), D = 2, the vertex 3 has hard degree 8 > 2D = 4, its hard component is the whole cluster and the soft pair ac is inserted. Fix: "no *hard* pair is edited".
The deduction "so relations realised by releasable tunnels have bounded degree" is not a theorem (the statements concern vertices incident to an edited hard pair); keep as HEURISTIC.
(b) Algorithm: guess a_1 (n ways), b_j among <= d partners in S_1j, then a_i among <= d partners of b_1 in S_i1: n d^{2k} checks. Tested against brute force on 9000 random instances (k = 2, 3; n = 3..5; d = 1, 2): identical answers.

### Prop 9 CORRECT WITH FIX
(i) monotone: cutting more tunnels only removes constraints - splitting each respected cluster {x,x'} into singletons changes no edit (Lemma 2: x, x' have no other cluster mates). (ii),(iii) definitions (lockable sets are closed under subsets by (i)). Algorithm: find a minimal non-lockable subset of T minus R by deleting tunnels one at a time (|T| calls of the l = 0 DP), branch on its <= d elements, depth <= l: d^l (3^{w+1} n^{O(1)} |T|) time.
Fix (model): the identity "min |R| with f(R) = 1 = optimum of the model" holds only if the model is the set of solutions with one edit per member and only tunnel cuts as hard edits. Unrestricted OPT - |H| can be smaller (411 of 931 tunnelled instances with a finite hitting-set value; e.g. a single member receives two edits instead of cutting two tunnels). Also note that cores are defined from G with the tunnel edges, w and the DP do not change when tunnel edges are cut (I(H) depends on H only) - right.
Tests (`a8_prop9.py`, 1200 tunnelled instances: 2-3 tunnels, stubs and paired tunnels, n <= 10, all partitions): monotonicity, "f(R)=1 iff R meets all cores", brute-force TR optimum = minimum hitting set, and the search-tree algorithm output = minimum: 0 exceptions.

### Independent re-verification of the computer evidence (`a8_flip.py`, `a8_core3.py`)
* G* (vertices 0,1,3,4,5,6,7; edges 01,04,05,14,15,45,03,13,36,37,06,56; members (6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1); anchor x = 7): valid, pair-maximal, OPT = 5 = |H|. Exactly five ZE solutions:
  `{0,1,3,4,5,6},{7}` (states RMMMM), `{0,1,4,5},{3,6},{7}` (RRRRL) with x a singleton; `{0,1,4,5,6},{3,7}` (LRRMM), `{0,1,4,5},{3,6,7}` (MRRRL), `{0,1,4,5},{3,7},{6}` (LRRRL) with x not a singleton.
  Locked relation of (m3,m1): {(M,M),(R,R)}; released: contains (M,R). Mini-track (t = 8 adjacent to 0,1,4,5,6): OPT 5, unique optimum `{0,1,4,5,6,8},{3,7}`; with a stub on x: OPT 6. c(G*) = 3. Sigma(G*) is not 2-decomposable (RRRMM pairwise consistent, not a solution).
* Cost sharing: two disjoint P3's (0,1,2), (3,4,5) with hard pendant edges 2-6, 5-7 (which force state R in both): OPT 2; stubs on both ends 4; one tunnel 0-3: 3; pendant edges removed (state L available) with the tunnel: 2.
* Small cores: all connected graphs of the networkx atlas with 4..7 vertices (tight = OPT equals the packing number; faces of arbitrary vertex sets via optimal partitions): maximal core size 2 in every case; n = 8, 9 random connected graphs (about 1,700 tight): maximal core size 2. Reproduces the notes' claim (evidence only).

### Conjecture C2 is false as stated (WRONG claim)
Wave 6 Sec. 4.7: "c <= 2: all instances closed under M-majority (one explicit rule: all-distinct triple -> M)". Counterexample found by `a8_c2b.py` (n = 6, c = 2, I(H) = K_4):
edges 03,05,12,14,23,24,25,34,35,45; H = (0,5,4), (1,2,5), (0,3,2), (1,4,3) (members m0..m3, each vertex in at most 2 members); OPT = 4 = |H|. The four ZE solutions have state vectors
`LLLL, LMLM, MLML, RRRR`. M-majority of (LMLM, MLML, RRRR) is MMMM, which is not a solution. (The same holds for a 7-vertex instance, `a8_c2b.py` prints it with an independent networkx cluster-graph test: Sigma = {RMMR, LRRL, MLLM, RLLR}, no uniform majority rule at all closes Sigma.)
Frequency among random tight instances with c <= 2, 3 <= |Sigma| <= 50: 77, 77, 74 of 20,000 (n = 6..9); n = 6: 82 of 6000; n = 8, 9: 1 of 6000. Wave 6's samples (n = 8, 9) were simply too rare an event to see.
What survives: in the same samples *no* c <= 2 instance failed 2-decomposability (Sigma determined by binary projections) and no minimal non-face of size >= 3 occurred. So the correct open statements are (C2') every minimal non-face has at most two anchors when c <= 2 and (C2'') Sigma is 2-decomposable when c <= 2; neither follows from Lemma 7b as it stands because the hypothesis there is a uniform majority operation (and 2-decomposability is not obviously preserved by projection to the anchor coordinates).

## Consequences for the paper section
Use: Lemma 1, Lemma 2 (fix), inactive anchors (fix), Lemma 5/6 (fix), Lemma 7 (a), (b), Lemma 7b with the self-contained proof (hypothesis: uniform majority), Lemma 8 (fix) and bounded-degree Grid Tiling, Prop 9 (fix, model defined), the verified G*, the sharing example, small-core evidence and negative composition results *as evidence only*; the persistent-flip chain and Conjecture C2 as open problems, with C2 restated as C2' (cores of size <= 2 when c <= 2) and the false uniform-M-majority form recorded as refuted.
Do not use: Lemma 4's second sentence, "frozen = no pair edited", the claim "closed under M-majority for all c <= 2", the claim that Sec. 3.1's reduction is anything but a sketch.
