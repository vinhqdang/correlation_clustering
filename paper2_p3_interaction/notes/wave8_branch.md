# Wave 8: branching on the far edits (cores of respected far pairs)

Labels: PROVED (proof here), VERIFIED-BY-COMPUTER (finite check, evidence only), HEURISTIC, OPEN.
Code (own, independent of earlier waves; numpy, scipy, networkx only): `code/w8b_lib.py` (instances, set-partition enumeration,
tree decomposition, local/non-local pairs), `w8b_branch.py` (brute-force oracle, the branching algorithm, exact OPT),
`w8b_ilp.py` (ILP oracle, FAR_m), `w8b_csp.py` (triple-dropping variant of the relaxation), `w8b_nstar.py` (N*(F)),
`w8b_stats.py`, `w8b_stats0.py`, `w8b_adv.py`, `w8b_c1.py` (core statistics), `w8b_attach.py` (Lemma A check). paper/ untouched, no git.

## 0. Verdict

1. FPT in (l,w) is NOT proved and NOT refuted. I did not find an algorithm and did not find a hardness proof.
2. What is new and PROVED: (a) a clean reformulation of zero excess ("attachment" lemma, Sec. 1); (b) a branching theorem (Sec. 3):
   CEaMP is decided by a search tree of depth <= l whose branching factor is the size d of a *minimal certificate* T, a set of
   *far hard pairs* that, when all respected, already force excess > l. Time d^l times the cost of the oracle. This is exactly the
   BFK c=1 O(4^l) scheme with "induced P3" replaced by "certificate T", and it is **complete and sound for every packing**.
3. The open question is thereby reduced to ONE precise statement (Core Bound, Sec. 5): d <= g(l,w) for every certificate encountered.
   Chain length does not enter d (FAR_m has d = 3 at l=0 and d = 5 at l=1 for all m = 3..8, Sec. 4). A family with unbounded d is
   exactly a persistent-flip chain (PFG) in the sense of wave6_hardness.md, now phrased for arbitrary far pairs, not only tunnels.
4. What I could NOT do: (i) prove the Core Bound; (ii) build a family with unbounded d (the Hall-chain target of Sec. 6 needs a centre
   with three attachment options, which a single member does not give; no gadget found); (iii) give an FPT *oracle* for the leaf test when
   local edits are allowed (the representability gap, Sec. 3.3). So the theorem of Sec. 3 is a reduction to {Core Bound, leaf oracle}.

## 1. Zero excess = attachment (PROVED, VERIFIED-BY-COMPUTER)

**Lemma A.** (G,H,0) is YES iff there is a partition pi of V with (i) every hard edge inside a block, (ii) every hard non-edge across blocks,
(iii) for every member (a,b,c) with centre b, b is in the block of a or of c.
Proof. A solution with |S| = |H| has exactly one edit per member and no hard edit (Lemma decomp). For a partition the cost of the three pairs of (a,b,c)
(ab,bc edges, ac non-edge) is [ab split]+[bc split]+[ac joined]; the patterns (ab,bc,ac) = (1,1,1), (1,0,0), (0,1,0) have cost 1, (0,0,0) costs 2,
(1,0,1)/(0,1,1) are not transitive, and (1,1,0) is not transitive. So cost 1 iff b is joined to a or to c (then the status of ac is forced by transitivity).
Hard pairs unedited gives (i),(ii). QED. Check: `w8b_attach.py`, 462 random instances (n<=8, chain and c=1 families), 170 with ZE, 0 mismatches against exact OPT = |H|.
Reading: every centre must *attach* to one of its two ends; a hard clique that is not releasable is a rigid block; the clusters are cliques of
the "allowed graph" (hard edges and soft pairs). An unreleased tunnel xx' with N_h(x)={x'} makes x unattachable (Lemma tlock).

## 2. Local and far pairs, the relaxation F(T)

Fix a tree decomposition of I(H). For v in V(H) let X_v be the set of bags containing all members through v (non-empty: they form a clique).
A hard pair {a,b} is *local* if X_a and X_b meet, *far* (non-local) otherwise; pairs with an end outside V(H) are far. Soft pairs are always counted.
Let N be the set of far hard pairs. For T subset N define F_G(T) = min over partitions pi of V of  #{counted pairs disagreeing with G}, over
partitions that respect every pair of T (a hard edge in T is joined, a hard non-edge in T is separated); far pairs outside T are free
(they cost nothing and are unconstrained by T).

**Lemma R (certificate).** (a) F is monotone: T subset T' implies F(T) <= F(T'). (b) If S is a solution with S cap T = empty, then F(T) <= |S|. (c) F(N) = min size
of a solution that uses no far edit. Hence: if F_G(T) > |H|+l then every solution of size <= |H|+l contains a pair of T.
Proof. (a) fewer constraints. (b) the partition of S respects T because the pairs of T are unedited, and its F-cost counts only a subset of S.
(c) all far pairs respected means the far part of S is empty. The last sentence is (b). QED. (The csp-variant that *drops triples* containing a far
pair outside T, `w8b_csp.py`, satisfies (a)-(c) as well and is weaker or equal.)

## 3. The branching theorem

### 3.1 Algorithm
solve(G,l): B := |H|+l. If F_G(N) <= B, accept. If l = 0, reject. Otherwise T := N and delete pairs one by one while F(T) > B stays true
(<= |N| oracle calls); the result is a *minimal certificate* T. For every e in T call solve(G xor e, l-1); accept iff one accepts.

**Theorem B1 (PROVED).** solve is correct for every packing H (no pair-maximality needed), and its search tree has depth <= l and branching factor
d := max |T| over the certificates met; the number of oracle calls is <= d^l (|N|+2).
Proof. Soundness: a solution S' of G xor e of size <= B-1 gives S' xor {e}, a solution of G of size <= B; H stays a valid packing in G xor e since
e is a hard pair. Completeness: let S be a solution of size <= B. If F(N) <= B we accept at once (this requires only that the leaf test is an exact "far-edit-free" test, Lemma R(c)).
Otherwise F(N) > B, so a minimal certificate T exists; by Lemma R every solution of size <= B contains some e in T, and S xor {e} is a solution of
G xor e of size <= B-1, found by induction. At l = 0 acceptance iff F(N) <= |H|, i.e. iff a solution of size |H| exists. QED

### 3.2 What it says
* The step "pick a P3 of the current graph and branch on its 3 pairs" (BFK) is replaced by "pick a minimal set of far pairs that jointly force excess".
* The order of far edits is irrelevant; each recursion level spends exactly one unit of l on a far hard edit. Local hard edits and the excess *inside* members
  (kappa_i >= 2) are never branched on; they are left to the oracle, which is why the long implication chains of FAR_m do not enter d (Sec. 4).

### 3.3 The oracle (what is proved, what is not)
* F_G(T) for a set T of bounded size: PROVED computable by the standard partition dynamic program on a tree decomposition of the vertex graph J
  (counted pairs, width <= 3(w+1)-1) plus one extra edge per pair of T (width <= 3w+2+|T|): time Bell(3w+3+|T|) n^{O(1)}. So the certificate tests and the shrinking are FPT in (w,d).
* The *leaf test* F_G(N) (exact search for a solution without far edits, local edits allowed) is the gap. With member variables only (no edit variables at all, i.e. N = all hard
  pairs) it is the domain-7 CSP of Theorem A, but then local pairs belong to T and d is NOT bounded: for FAR_m the l=0 certificate would contain the m-1 chain pairs {2i+1,2i+3}.
  With local pairs as bag-variables, every triple has a bag-local scope **except** one type: two local hard pairs {a,y},{b,y} and one far pair {a,b}
  (regions X_a cap X_y and X_b cap X_y can be disjoint inside X_y). Using a decomposition of I(H)^2 (width <= (w+1)(3w+1)-1, Prop. bag of the paper) all triples with a soft pair are fine.
  I do not have an FPT leaf test for the exceptional type (it is the constraint "do not make a P3 by two local edits through y against a fixed far pair"). OPEN (technical).
  The experiments below use an exact leaf test (set partitions / ILP), so they test Theorem B1 and the core sizes, not the oracle.

## 4. Computer verification (VERIFIED-BY-COMPUTER, evidence only)

Instances: "chain+extras" (path of m members, I(H) a path, c=2, w=1, random extra hard edges, valid packings), "c=1" (3m vertices, vertex-disjoint members, random extra edges),
FAR_m, N*(F) of wave6_ml (Prop 5).
* Correctness of Theorem B1 versus exact OPT, l = 0..3: n = 7: 27 instances, n = 9: 251, n = 11: 44 (set-partition oracle); n = 11, 13, 15 (m = 5,6,7): 75 instances (ILP oracle,
  l up to min(excess,2)). **0 disagreements** in all.
* FAR_m, m = 3..8 (n = 7..17, ILP): l=1 accepted after 4 nodes; the minimal certificate has size **5 for every m** ({(0,2m-1),(1,2m-2),(1,2m-1),(1,2m),(2,2m-1)} for m=5: (0,9),(1,8),(1,9),(1,10),(2,9));
  at l=0 the certificate has size **3** ({(1,2m-2),(1,2m-1),(1,2m)}). The chain length m does not enter. Same with the triple-dropping relaxation.
* N*(K4) (n = 30, w = 1, tw(Gamma') unbounded family of wave6_ml): l=1 accepted after 5 nodes, certificate size 6; l=0 rejected at the root.
* Certificate sizes (minimal T), chain+extras, ILP oracle, 60 instances each: l=0: m=4: sizes 0..4; m=6: 1..4; m=8: 1..3; adversarial search (80 instances, m=6, up to 3m extra edges): max 3.
  l=1: m=4: 3..6; m=6: 2..7; m=8: 3..7. Over all budgets (l<=3), n=9: max 9, n=11: max 11 (these are sums over several independent obstructions, one per unit of excess).
* c=1 (vertex-disjoint packings, every hard pair is far), m=4 (n=12): l=0: 3..4; l=1: 3..6; l=2: 3..9. Growth with l, none with n, consistent with the BFK counting.
No instance was found where the size of a minimal certificate at fixed l grows with n.

## 5. The precise open statement and the structure of a possible counterexample

**Core Bound (OPEN).** For pair-maximal H (or after MAX) there is g(l,w) such that every minimal certificate T with F(T) > |H|+l has |T| <= g(l,w).
Equivalent to FPT of this scheme given the oracle. It holds in all tests; it is the exact form of "locality of far edits".

**Proposition U (unbounded cores need a long carry; HEURISTIC, proof sketch only).** Suppose minimal certificates T_N with |T_N| -> infinity at l = 0. Let T_N = {q_1..q_N}. For each i pick a
partition pi_i of cost |H| respecting T_N minus q_i and violating q_i (it exists by minimality, since F(T_N minus q_i) <= |H|). Then any bag-separator of a decomposition of
I(H) splitting the far pairs q_i (their active ends) must carry, in the states of its bounded number of members, which side contains the violated pair: the splice argument of
wave6_hardness Lemma 7 holds verbatim, because the relaxed problem is again a CSP over bag-local triples. So unbounded cores need a *persistent bit* (PFG) with the
property: all N far pairs respected => infeasible, any N-1 => feasible.

Constraints on the far pairs (PROVED except for one triple type, Sec. 3.3): a triple with a far pair q has either two soft pairs (scope two adjacent members) or at most one other variable-carrying pair, so the effect of q on the member
variables is a **unary or binary-adjacent** constraint ("this soft pair of the member must be merged / split") guarded by the set of far pairs of the triple (|guard| <= 2).
This is why a tunnel acts as a lock and why guards are what T collects. Baker-Pixley (wave6_hardness Lemma 7b) bounds the core by 2 for *singleton-valued* unary guards when the base relation has a
majority operation; with multi-valued unary guards I could not prove a bound (the combination of three tuples can violate the coordinate whose own tuple was the odd one out). OPEN.

**Sharp mechanism for a counterexample (HEURISTIC).** In the attachment picture, a PFG is a Hall-type deficiency chain: centres y_1..y_N, regular partners r_1..r_{N-1} (deficiency 1),
and one releasable partner s_j (an anchor, locked = singleton) per centre; all locked: no perfect assignment; release any one: the deficiency shifts to that position and it works. A centre
with only two options (one member) cannot have the third option s_j; two members at one centre force the cluster of y_j to contain an attached vertex of each, so the two attached vertices
must be soft or hard-adjacent, which costs extra members and so far gave occupancy (wave6_hardness Lemma 6). With two options per centre the chain is bijunctive (cores <= 2, consistent with
Conjecture C2 of wave6_hardness for c <= 2). So the target object is: a three-option centre of bounded c. I did not construct it.
If it exists, the Grid-Tiling blueprint of wave6_hardness Sec. 3.1 (tracks, far tunnels as a matching, l = k^2, w = O(k)) shows that CEaMP is W[1]-hard in l+w, and the scheme of Sec. 3 has d = Theta(n).
If Core Bound holds, CEaMP is FPT in (l,w) provided the leaf oracle exists. Either way the question is now a question about *gadgets for one ZE relation*, not about identities or pointers.

## 6. Comparison with the earlier obstruction (FAR)

FAR (Prop. FAR) forces a far two-sided edit at distance m-1 with unbounded chain, which is what killed "bag-local flag" encodings. In the branching view it is harmless: the chain is
solved by the oracle (local propagation), the certificate is the 3 pairs at the end vertex 1 versus the member of 2m-1, and the branching has <= 5 children. Conjecture NF-bad was about *which members*
carry a far edit; the certificate T says *which pairs*, and its size does not depend on the distance. This is the point of the new angle.

## 7. Status table

| item | status |
|---|---|
| Lemma A (zero excess = every centre attaches to an end) | PROVED; VERIFIED (462 instances) |
| Lemma R (certificate / monotone relaxation) | PROVED |
| Theorem B1 (branching on minimal certificates, depth l, factor d, any packing) | PROVED |
| Correctness of the algorithm vs exact OPT (n <= 15, l <= 3) | VERIFIED-BY-COMPUTER (0 / ~400 disagreements) |
| Certificate test F(T) in time Bell(3w+3+|T|) n^O(1) | PROVED (standard) |
| Leaf test with local edits as bag variables | OPEN (one non-representable triple type) |
| d independent of m for FAR_m (3 at l=0, 5 at l=1), d = 6 for N*(K4), d <= 4 at l=0 on chains up to n=17 | VERIFIED-BY-COMPUTER |
| Core Bound d <= g(l,w) | OPEN (no counterexample found) |
| Unbounded cores need a persistent bit (PFG) | HEURISTIC (splice argument as in wave 6, not rewritten) |
| Existence of PFG / three-option centre gadget | OPEN |
| W[1]-hardness | OPEN; FPT in (l,w) OPEN |
