# Wave 9 (prover A): FPT in l for soft-triangle-free instances / w <= 1 -- what is proved, what is conditional, what is false

Labels: PROVED = proof in this file (finite checks that belong to a proof are named and marked "by computer"); VERIFIED-BY-COMPUTER = exhaustive or sampled
check, evidence only; HEURISTIC = argument without full proof; OPEN. Code: `code/w9_*.py`, `code/w9a_*.py`, `code/w9_eng.c` (+`w9_engpy.py`), commands in Sec. 9.
No git, `paper/` untouched. Notation as in wave7_special.md (member (a,b,c), centre b, soft/hard pairs, ZE = zero excess, soft triangle, Sigma(G,H)).
Files `w9b_*` in `code/` belong to another agent and are not used here.

## 0. Verdict

1. **FPT in (l,w) for w <= 1 is still NOT proved, and no hardness proof was found.** I did not obtain the target theorem.
2. **PROVED, new, and fully verified:** *Theorem U*: for **every soft-triangle-free instance (any width, any c)** the "far-lift zero-excess model"
   (every member edited exactly once; hard edits U avoid the cross pairs of adjacent members, |U| <= l) is decidable in time 6^l n^{O(1)} by a bounded search tree.
   It generalises the tunnel-release theorem A3 of wave 7 to *arbitrary far hard edits* (static closure included) and removes the "forked lock" gap of wave 7 for
   this model: the liftable constraints are *unit* lists, so cores have at most 2 members (Lemma H2). 1,420 instances, 0 disagreements with exhaustive search.
3. **PROVED:** a *finite exact model* of CEaMP at soft-triangle-free instances (Sec. 1): member states are the 5 set partitions of {a,b,c}, interaction at a shared
   vertex depends only on the *views* of the two members at that vertex, a hard outside vertex acts on a member through the table L(S).
   Consequences: (i) the domain of the paper's DP can be 5 instead of 7; (ii) a break (state O/X, double edit) is not a node-hit but a **cut with locks**; (iii) a hard
   vertex can never penalise a break state O without also penalising L or R (Prop. 4.4) -- this kills the natural way to realise the Grid-Tiling reduction of wave 7 with breaks.
4. **Corrections to waves 7/8:**
   * The *packing-endpoint lemma* is not needed: RM2 is FPT by a trivial deepest-LCA search tree, (2d+1)^L (Prop. 4.1; code `w9a_rm2.py`, 2,100 instances).
   * **RM2 is not a faithful abstraction.** Explicit counterexample (n = 11, |H| = 4, I(H) a star, c = 2, pair-maximal): RM2 value 2, true excess 3 (exact DP over all partitions).
     Two breaks at two members adjacent to a common member lock that member and force a third break (Prop. 4.2, `w9a_hubcex.py`). So (L1) of wave 7 is false.
   * The "duality g = max weighted packing" that I first verified on 1.4 million instances is **vacuous** (the one-region packing attains it); the informative statement is a
     bounded-certificate statement (BC, Sec. 5). The earlier "verification" must not be quoted as evidence for anything.
   * The "leaf oracle gap" of wave8_branch.md Sec. 3.3 (two local hard pairs through y and one far pair) is closed by using enlarged bags N^1[beta] (Obs. 5.3, HEURISTIC but short).
5. **Conditional theorem (Theorem F, Sec. 5):** the full CEaMP at forest I(H) (w <= 1) is FPT in l **provided the Bounded Certificate hypothesis BC holds** (this is the "Core Bound"
   of wave 8 in member form). The implication (algorithm + proof) is written out; the last level of the search tree (one remaining unit of budget) is unconditional (BC(0) is Lemma H2).
   BC is OPEN; evidence: minimal certificates had size <= t+2 in every adversarial search (Sec. 6).
6. What I could **not** do: prove BC (weighted Helly for the extended domain) -- the position of interior breaks depends on locks of neighbours (Sec. 4/6), so an exchange argument
   in the style of RM2 is false and a Helly argument for the weighted problem is missing.

## 1. The exact finite model (PROVED; checked against exact OPT)

Let H be any modification-disjoint packing, G any graph. Cross pair := hard pair {p,r} with p in V(A)\{q}, r in V(B)\{q} for two members A,B sharing q (adjacent in I(H));
Far := hard pairs that are not cross pairs. Soft-triangle-free: no vertex triple whose three pairs are packed by three distinct members.

**States.** A solution restricted to V(A) is a partition of {a,b,c}; the triple must be a cluster graph. Of the 7 non-empty edit sets only 5 give a cluster graph:
L = {a}{bc} (edit ab), R = {ab}{c} (edit bc), M = {abc} (edit ac), O = {a}{b}{c} (edit ab,bc; extra cost 1), X = {ac}{b} (edit all three; extra cost 2).
(The paper's DP of Theorem A therefore has 5^{w+1}, not 7^{w+1}, table entries: PROVED, trivial; the sets {ab,ac}, {bc,ac} leave an induced P3 inside the member.)
Final adjacency y=(y_ab,y_bc,y_ac): L=(0,1,0), R=(1,0,0), M=(1,1,1), O=(0,0,0), X=(0,0,1).
**Views.** The *view* of a state at a vertex v of the member is the pair (is v clustered with each of the other two vertices): at a: (y_ab,y_ac), at b: (y_ab,y_bc), at c: (y_bc,y_ac).
L: a:(0,0) b:(0,1) c:(1,0); R: a:(1,0) b:(1,0) c:(0,0); M: all (1,1); O: all (0,0); X: a:(0,1) b:(0,0) c:(0,1).
**Binary constraints.** Two members A,B sharing q (they share at most one vertex): the only triples carrying soft pairs of both are {p,q,r}; with u_p=[p in q's cluster], v_r likewise and eps_pr in {0,1} the status of the
(hard, cross) pair: eps_pr=1 forces u_p=v_r, eps_pr=0 forbids u_p=v_r=1. So R_AB depends on the two views only (Lemma 1 of wave 7, now for all 5 states).
**Unary constraints.** For t outside V(A) with all pairs {x,t}, x in V(A), hard (status in the current graph), N(t) cap V(A) = S: the three triples {a,b,t},{b,c,t},{a,c,t} are equivalent to x_A in L(S):

| S | {} | {a} | {b} | {c} | {a,b} | {a,c} | {b,c} | {a,b,c} |
|---|---|---|---|---|---|---|---|---|
| L(S) | all 5 | {L,O} | {O,X} | {R,O} | {R} | {X} | {L} | {M} |
| L(S) on ZE states | {L,R,M} | {L} | empty | {R} | {R} | empty | {L} | {M} |

(VERIFIED by brute force over cluster graphs on {a,b,c,t}: `w9a_listtable.py`; it is also immediate: S must be empty or a single part of the partition.)
**Lemma E (exactness) PROVED.** For every packing H: OPT(G) = |H| + min over U subset of Far of ( |U| + Ex(U) ), where Ex(U) is the minimum over F subset of Cross and
x in {L,R,M,O,X}^H such that G xor U xor F xor (soft edits of x) is a cluster graph, of |F| + sum_i extra(x_i) (infinity if none).
Proof: Lemma decomp splits a solution into soft edits (one non-empty edit set per member, i.e. a state) and hard edits (a set of cross pairs and a set of far pairs); conversely every (U,F,x) is a solution.
VERIFIED-BY-COMPUTER (`w9_lib.py` + `w7_dpy`): min over U (|U|<=2) of |U| + m_U equals the exact OPT on 60 random pair-maximal instances (n <= 11).
Lemma 2 of wave 7 survives verbatim: a far pair changes only the *vertex lists* at the members through its two ends and the static (all-hard) triples; it never changes a relation.

## 2. Unit lists have Helly number 2 (PROVED, any width)

**Fact (wave 7, Theorem A).** For soft-triangle-free (G',H) the set Sigma_perm of zero-excess state vectors satisfying all *permanent* constraints (own triples, relations, every unary triple all of whose hard pairs are cross pairs)
is closed under every majority operation. (All 40 relations of Lemma 1 are closed under all 729 uniform majorities: re-verified, `w9a_helly.py` header, 0 violations.)
**Lemma H2 (PROVED).** Let Sigma be closed under a majority operation m and let J be a family of *unit lists* "x_{v_j} = s_j" (or "empty"). If Sigma and J are jointly infeasible, some <= 2 members of J already are.
Proof. Take a minimal infeasible J, |J|=k >= 3, tuples t_j in Sigma satisfying J minus {j}. For a,b,c in J put u = m(t_a,t_b,t_c) in Sigma. For j not in {a,b,c} the three tuples all have value s_j at v_j. For j = a the tuples t_b,t_c both have s_a at v_a,
so the majority is s_a (likewise b,c). So u satisfies all of J: contradiction. QED (no condition on the domain size; this is the classical 2-SAT core argument).
On zero-excess states every **liftable** unary constraint is a unit list or "empty" or trivial (table above), because the three triples of one outside vertex t are grouped into L(S): *this is why the forked locks of wave 7 never occur for zero-excess feasibility.*
(Forked lists {L,O},{O,X},{R,O} exist only on the extended domain, i.e. for states O/X.) For arbitrary 2-element lists Helly fails (three lists {L,R},{R,M},{L,M} on one variable), so the grouping matters.
Strictness (wave 7 Lemma S) in the form needed for trees was re-verified independently: for all 40 relations and both directions, X <= Y <= {L,R,M}: img(X) != img(Y) implies |img(X)| <= 1 (2,160 checks, 910 strict, 0 violations; `w9a_strict.py`); abstract forest CSPs over the 40 relations with arbitrary lists (<= 2 per node): minimal infeasible families have size <= 2 in ~10^5 cases (`w9a_helly.py`).

## 3. Theorem U: the far-lift zero-excess model is FPT (PROVED; VERIFIED-BY-COMPUTER)

**Model FLZ.** Input (G,H,l), H a packing, soft-triangle-free. Question: is there U subset of Far, |U| <= l, such that (i) G xor U has no induced P3 all of whose pairs are hard, and (ii) Sigma(G xor U,H) != empty (every member edited exactly once, nothing else)?
A YES answer gives a solution of size |H|+|U| of CEaMP (so FLZ is an upper bound; it is exact when some optimum has no break and no cross-pair edit). Tunnel release (A3) is the special case where U consists of cuts of 2-vertex hard components.
**Theorem U.** FLZ is decidable in time 6^l * n^{O(1)} (for forests with the DP of Sec. 1; for general widths by path consistency, cited Jeavons-Cohen-Cooper, as in A1).
*Algorithm* solve(U,b): (1) if G xor U has an all-hard induced P3, branch on its (editable = non-cross, not in U) pairs, b -> b-1 (<= 3 branches); (2) else if Sigma(G xor U) != empty return YES; (3) if b = 0 return NO;
(4) the *liftable lists* are the pairs (A,t), t outside V(A), all three pairs {x,t} hard, at least one of them in Far; delete all liftable lists: if still infeasible return NO (no far edit helps, because far edits change only liftable lists);
(5) otherwise shrink greedily to a minimal infeasible family F of liftable lists; by Lemma H2 |F| <= 2; branch on the <= 3 editable pairs {a,t}, a in V(A), of each list of F (<= 6 branches), b -> b-1.
*Correctness.* Soundness is clear (YES only at a feasible U). Completeness: let U* be a solution. A far edit changes Sigma_perm not at all (it changes no relation, and no permanent unary constraint: those have only cross pairs). The family F together with Sigma_perm is infeasible, hence U*
must change at least one list of F, i.e. contain a pair {a,t} with a in V(A) for some (A,t) in F; a static all-hard P3 must be hit by a pair of U*. Branching on that pair and recursing with U* minus the pair keeps a solution. QED
*Verification.* `w9a_thmU.py` compares the algorithm with exhaustive search over all U in Far (|U| <= l): seeds 1 (120), 31,32 (2 x 400, forest I(H), n <= 10), 34 (400, I(H) with cycles allowed but no soft triangle), 35 (100, l = 3): **1,420 instances, 0 mismatches**,
largest minimal family met: 2 (the assertion "<= 2" of Lemma H2 never failed; with triple-level lists, ungrouped, the bound 3 of the majority argument also held, family size <= 2).
*Remark.* The model is what a careful reading of "cores <= 2" (Conjecture C2) buys: a lock core is a *pair of vertex lists*, so a far-lift problem is a vertex-cover-like problem and the branching factor does not depend on n, D or w.

## 4. Break semantics, RM2, and why RM2 is not faithful

**Prop 4.1 (RM2 is FPT; PROVED).** RM2 (wave 7 Sec. 5): forest, conflicts = pairs of nodes with their path, operations = break at a node (cost 1, hits every conflict whose path contains the node) and tunnel (a link of a graph of maximum degree d, cost 1, hits every
conflict with an endpoint at one of its nodes). Search tree: root every tree, let Q be an unhit conflict with deepest LCA v; some solution contains the break at v or a tunnel at one of the two ends of Q: branching <= 1 + 2d, depth L, time (2d+1)^L poly.
Proof. If Q is hit by a break at w in Q, replace it by a break at v = lca(Q): every unhit conflict R through w has lca(R) at depth <= depth(v) and contains w, hence contains the path from w to lca(R), which passes through v. If Q is hit by a tunnel, it has an endpoint of Q. QED.
VERIFIED-BY-COMPUTER: `w9a_rm2.py` (2,100 random instances, exhaustive search vs. search tree: 0 mismatches, a tunnel strictly helps in 245 of them). So the "packing-endpoint lemma" (wave 7 Sec. 5) is not needed for RM2.
**Prop 4.2 (RM2 is not the problem; PROVED/VERIFIED).** n = 11, E = {01,02,04,05,06,13,28,2-10,34,37,38,78}, H = {(2,0,1),(4,3,1),(5,0,6),(2,8,7)}: valid, pair-maximal, I(H) = star with centre (2,0,1), c = 2, no soft triangle.
Members M1=(4,3,1) and M2=(5,0,6) are each ZE-infeasible alone (their vertex lists are {O,X}), each is repaired by a single break (O), so RM2 predicts excess 2; but the exact optimum is |H|+3 = 7 (subset DP over all partitions;
the CSP optimum over U = empty is also 3): the two breaks put the views (0,0) at the shared vertices 1 and 0 of the common neighbour M0, which then cannot be ZE, so M0 breaks too. Script `w9a_hubcex.py`.
Statistics (`w9_rm2real.py`): in 19 of 300 random forest instances the hitting-set value of the *minimal ZE-infeasible regions* (breaks as node-hits) is strictly smaller than the true excess g.
**Prop 4.3 (a break is a cut with locks; PROVED from the view table).** The view of O is (0,0) at every vertex; so a neighbour at a sees O as L, at c as R, at b a new view (0,0). Hence breaking A deletes A and imposes on each neighbour the lock "compatible with the view of O at the shared vertex";
two breaks next to a common member lock it twice. Whether a break repairs a given obstruction depends on these locks: in 4,910 minimal ZE-infeasible regions (`w9_mb.py`), 3,882 are repaired by a break at *every* member (all sizes >= 3 and all but 289 regions of size 2),
and 289 weight-1 regions of size 2 are repaired by a break at only one of the two members (e.g. H = {(1,0,2),(4,0,3),(5,2,6)}, region {M0,M1}). So positions matter (HEURISTIC description: valid cut positions along a chain are {j : left-forced state compatible with the L-view and right-forced state with the R-view}, a typed hitting problem).
**Prop 4.4 (O-conflicts are dominated; PROVED by the table, `w9a_listtable.py`).** For every outside vertex t, L(S) contains O or L(S) is one of {L},{R},{X},{M}. In particular no list excludes O while allowing both L and R.
Consequently no hard vertex (hence no far edit) can make *the break state of a chain member* more expensive than both its zero-excess states; "a link that rewards two far members both being broken" cannot be built from vertex lists.
This is the exact reason why the Grid-Tiling blueprint of wave 7 (links hitting interior nodes of obstruction paths) has no realisation: interior lifts are impossible (cores have <= 2 lists, Lemma H2) and breaks carry no link-shareable conflict (Prop 4.4).
Caveat: in the weighted problem an interior kill can matter (alternating kills L,R,L,R,... along a chain need one break per adjacent pair, and lifting an interior kill removes two of them); but the conflict structure on kills is then a graph (cores of size 2) and a lift saves at most 2, so the problem is vertex-cover-like (the number of lifted kills is <= 2l): no N-way choice by a single lift is possible.

## 5. The full problem: algorithm and conditional theorem

Setting: H arbitrary packing, I(H) a forest (soft-triangle-free automatic), budget l. By Lemma E the question is whether some far set U with |U| <= l has |U| + Ex(U) <= l.
**Exact oracle (Obs 5.3, HEURISTIC, argued below).** For a fixed far set U, Ex(U) is computable by a dynamic programme over a tree decomposition of I(H) with bags enlarged to N^1[beta] (members within distance 1): variables = states of the members and the edit status of cross pairs inside the bag;
every triple that contains a soft pair, and every all-hard triple that contains a cross pair, has its three vertices inside one bag (two pairs through y: all members through y lie in one bag, their neighbours lie in its N^1; the third pair is far, hence a constant fixed by U).
(This closes the "exceptional triple type" that wave8_branch.md Sec. 3.3 lists as OPEN. I checked it only by the argument just given and by the end-to-end test below.)
**Static violations** are handled by branching (<= 3 pairs), so no pair-maximality (Lemma MAX) and no extra members are needed; the partner of an edit is bounded by Cor. CAND of the paper applied to the current graph.
**Relaxation.** For a set C of *sites* (members and edges of the forest) let Ex_C(U) be Ex(U) with all liftable unary triples at members outside C and all static triples with a far pair at edges outside C deleted.
A far edit set W with ends in the members Z(W) (members through its endpoints, and the edges at them) satisfies Ex_{relax Z(W)}(U) <= Ex(U cup W) (every constraint that survives is unchanged by W; W changes no relation).
**Hypothesis BC(t) (OPEN).** There is K(t,D) such that whenever Ex(U) > t, every inclusion-minimal set C of sites with Ex_C(U) > t has |C| <= K(t,D).
(BC(0) holds with K = 2 by Lemma H2 together with zero-excess feasibility being a CSP with unit liftable lists. BC is the **Core Bound** of wave8_branch.md Sec. 5 in member form.)
**Theorem F (conditional on BC).** If BC holds, CEaMP for forest I(H) is FPT in l: time (c * K(l,D) * Z^2)^l * poly * 5^{O(1)}, Z = 3D + 3l + 1, D = 4.
*Algorithm* A(U,r): v := Ex(U); if v <= r accept; if r = 0 reject; branch on static all-hard P3 (3 pairs) if any; otherwise C := minimal set with Ex_C(U) > r-1 (greedy deletion with the oracle);
for every far pair e with an endpoint z in V(C) and the other endpoint in the ball of radius 2 around z in G[L_B] (<= Z^2 candidates, Cor. CAND): call A(U cup {e}, r-1).
*Proof.* Soundness: every call evaluates real solutions (Lemma E). Completeness: let W be a nonempty far set with |W| + Ex(U cup W) <= r < Ex(U) (so |W| >= 1 and Ex(U cup W) <= r-1). Z(W) is a set of sites with Ex_{relax Z(W)}(U) <= r-1, so Z(W) meets C
(otherwise relax Z(W) relaxes at least the complement of C and Ex_C(U) <= r-1). Hence some e in W has an endpoint in V(C) (V(C) = vertices of the members in C and of the members at the edges in C); the partner is in the ball of radius 2 (Cor. CAND, applied to the current graph G xor U and an optimal solution that contains W); recurse with W minus {e}. Depth <= r, each node has <= 3|C| Z^2 candidates. QED
*The last level is unconditional:* for r = 1 the set C has <= 2 members by Lemma H2 and the algorithm is Theorem U's step (5). So **l <= 1 is polynomial** (already known) and the search tree for general l needs BC only for t >= 1.
**End-to-end test (VERIFIED-BY-COMPUTER, `w9a_condalg.py`, `w9a_condalg2.py`).** Brute-force oracle (breaks, cross-pair flips <= r, static ignoring in the relaxation), minimal C by deletion, branching on far pairs at ends of C (any partner), versus the exact OPT (subset DP):
310 small instances (n <= 13, |H| <= 4), including 13 whose answer requires far edits and 81 NO instances: 0 mismatches. (This tests the accounting of Lemma E and of the relaxation, not the bound K.)

## 6. Evidence on BC and negative results (all VERIFIED-BY-COMPUTER, i.e. evidence only)

* **The "duality" is vacuous.** g = max over disjoint member sets of sum of sub-instance optima is attained by the single region "all members". Tests `w9_dsearch.py` (1.4 M hill-climbed forest instances) and earlier scratch runs found 0 gaps -- meaningless. What matters is how many members a *minimal* region/certificate needs.
* **Relevant sets.** Members whose single relaxation reduces g: <= g+1 in adversarial search (`w9_relevant*.py`; g = 1..3, m <= 18, paths and trees; found (g=10, rel=11)); members in minimal reducing pairs (`w9_rel2.py`): <= 2g+1, no growth with chain length.
* **Minimal certificates** C with relaxed optimum > t (`w9_bc.py`, hill-climbing on |C|, m <= 20): sizes 1 (t=0), 3 (t=1), 4 (t=2): <= t+2; for t = 0 the Helly-2 bound is 2.
* **Single-end criticality is not enough.** Relaxing two members can reduce g by more than the sum of single reductions (`w9_subadd.py`: reduction larger than the sum of the single reductions in ~4% of the member pairs, e.g. 3 vs 2), so a branching rule must be based on sets, not on single "critical" members (this is what BC handles).
* **Break positions are not bounded by the ends of the region.** Restricting breaks to the kill members and leaves of the minimal ZE-infeasible regions of the *original* instance fails in 8 of 400 forests (cascades through a hub: the neighbour of a broken member must break too); and in 1 of 941 regions no optimal assignment breaks a leaf or kill member of that region (the optimum breaks an interior member, because breaking a leaf would force a break of its neighbour outside the region) (`w9_nk2.py`, `w9_bl3.py`). So breaks must stay inside the DP; only far lifts are branched.
* **Local flips act as lifts at distance <= 3.** A cross-pair flip {p,r} with p in A cap A'', r in B cap B'' also changes the unary lists of A'' and B'': member-disjoint region weights computed with region-local flips overshoot (29 of 200 instances, `w9_dloc.py`). So the DP must treat all cross-pair flips globally (it does, Sec. 5).
* **Finite profile algebra.** Rooted-subtree profiles (cost of the subtree per view of its root at the port vertex), normalised and saturated at CAP, form a finite closed set: CAP=1: 67 profiles, 41 induced functions; CAP=2: 115/64; CAP=3: 175/92 (`w9_closure.cpp`). So the forest problem for a fixed set of far edits is a finite-state tree automaton with costs; BC is a statement about the Carathéodory number of its unary lists. I did not manage to turn the algebra into a proof.

## 7. Why BC is plausible, and what a proof needs (HEURISTIC)

For t = 0: Lemma H2. For t >= 1 one wants: if every structure y of cost <= t (breaks, flips) leaves Sigma_y infeasible with the unit lists C, then a bounded subset of C already does for *all* y simultaneously.
Each Sigma_y is again a soft-triangle-free forest CSP (a break removes a member and adds locks), so Lemma H2 gives <= 2 lists per y; the missing step is a covering argument over the unboundedly many positions of y.
The argument sketch: a cut at a position outside the Steiner tree of the first core leaves the core intact; a cut inside it splits it into pieces, each piece is infeasible because of its *nearest forcing kill* and a lock, and a forcing kill forces the same singleton at every node of its functional chain, so the same kill serves all positions of one segment.
Tests (Sec. 6) agree (|C| <= t+2). What is not established: interaction of several cores and the locks at hubs (Prop 4.2), and the extended-domain relations are not closed under any majority that would make Helly arguments uniform.

## 8. Status table

| item | status |
|---|---|
| Lemma E: exact model (5 states, views, L(S), cross/far split), 5^{w+1} DP | PROVED; VERIFIED (60 instances vs exact OPT) |
| Lemma H2: unit lists have Helly number 2 (any width) | PROVED |
| Theorem U: far-lift ZE model FPT, 6^l n^{O(1)}, any soft-triangle-free instance | PROVED; VERIFIED (1,420 instances, 0 mismatches) |
| Prop 4.1: RM2 FPT, packing-endpoint lemma unnecessary | PROVED; VERIFIED (2,100 instances) |
| Prop 4.2: RM2 is not faithful (hub counterexample: 2 vs 3) | PROVED/VERIFIED (exact DP) |
| Prop 4.3: break = cut with locks (view table) | PROVED |
| Prop 4.4: O-conflicts dominated by L/R conflicts (no link reward for breaks) | PROVED (table, `w9a_listtable.py`) |
| leaf oracle for local edits via N^1-bags (closes wave 8 gap) | HEURISTIC (argument + end-to-end test) |
| Theorem F: FPT in l for forests given BC | PROVED implication; BC OPEN |
| BC(0) (<= 2 members) | PROVED; BC(t), t >= 1 | OPEN, evidence |
| "duality g = max packing" | VACUOUS (do not use) |
| FPT in (l,w) for w <= 1 unconditionally | **OPEN** |
| W[1]-hardness for w <= 1 | OPEN; the two natural mechanisms (interior lifts, link-rewarded breaks) are excluded (Lemma H2, Prop 4.4) |

## 9. Reproduction (all in `code/`; python3 + numpy + networkx; gcc for the engine)

`gcc -O2 -shared -fPIC -o w9_eng.so w9_eng.c` (engine: CSP tables + forest DP, used by `w9_engpy.py`); `python3 w9a_listtable.py`; `python3 w9a_strict.py`; `python3 w9a_helly.py 1 40000 7`;
`python3 w9a_thmU.py 31 400 2 10 forest` (also seeds 32, `34 400 2 10 any`, `35 100 3 8 forest`, `1 120 2 9`); `python3 w9a_rm2.py 1 600`, `w9a_rm2.py 2 1500`; `python3 w9a_hubcex.py`; `python3 w9_rm2real.py 1 300 7`;
`python3 w9a_condalg.py 1 40 2 9`, `python3 w9a_condalg2.py 1 150 1`, `w9a_condalg2.py 2 120 2`; `python3 w9_mb.py 3 3000 7 90`; `python3 w9_relevant3.py 2 60 6000 14 2 path`; `python3 w9_rel2.py 7 25 3000 18 2 tree`; `python3 w9_bc.py 11 80 4000 16 1 path`;
`python3 w9_nk2.py 1 400 7`, `python3 w9_bl3.py 2 500 8 leaf`; `python3 w9_subadd.py 1 1500 8`; `python3 w9_dloc.py 1 200 6 2`; `g++ -O2 -o w9_closure w9_closure.cpp && ./w9_closure 3`.
(`w9_dsearch.py`, `w9_bl.py`, `w9_bl2.py`, `w9_breakpos*.py`, `w9_minbreak.py`, `w9_relonly.py` are exploratory; `w9_dsearch.py` tests the vacuous duality.)

## 10. What would settle the question
* **FPT:** prove BC(t) for t >= 1 (Carathéodory number of the unit lists for the weighted, extended-domain problem) -- the induction suggested in Sec. 7 needs a statement that the set of valid cut positions of a minimal core is determined by boundedly many forcing kills; with BC, Theorem F applies; the oracle of Obs. 5.3 should be written out.
* **Hardness:** one would need an interaction that is neither an interior lift nor a link-rewarded break: by Lemma H2 and Prop 4.4 the far-acting objects are vertex lists with <= 2-element cores; any W[1]-hardness proof must therefore live in the weighted regime (heavy regions, alternating kills) where lifts are vertex-cover-like, and no construction is in sight.
