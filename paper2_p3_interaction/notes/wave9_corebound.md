# Wave 9b: the Core Bound and the leaf test, settled

Labels: PROVED (proof given here), VERIFIED-BY-COMPUTER (finite check, evidence only unless a proof says otherwise), HEURISTIC, OPEN.
Code: `code/w9b_*.py` (written from scratch; numpy, scipy, networkx, python-sat, ortools). Two independent exact oracles throughout:
`w9b_leaf.py` (bucket elimination over bag-local variables) against `w9b_lib.opt_cpsat / opt_csp_cpsat` (CP-SAT over pair variables).
paper/ untouched, no git. Reproduction commands are in Section 9.

## 0. Verdict

1. **Leaf test: PROVED FPT** (Theorem L). The "failing triple type" is an artefact of using X_a cap X_b of a decomposition of I(H) as the definition of *local*.
   With *local := the two members through the ends are equal or adjacent in I(H)* (a pair has at most 9c^2 local partners per vertex) every vertex triple of the leaf CSP
   has all its variables inside one bag N^3[beta_t] of the decomposition of I(H)^4 (radius 3 is attained, Section 2.3). So the leaf test, and every certificate test F'(T) for arbitrary |T|,
   is a bounded-treewidth min-sum CSP: time f(w) n^{O(1)}. Verified against CP-SAT on 280 instances (0 disagreements), and the whole branching algorithm B1 against exact OPT on 449 (instance,budget) pairs (0 mismatches).
2. **Core Bound, l = 0, no soft triangle: PROVED** (Theorem CB0): every minimal certificate has at most 12 far pairs (3 if it is a static violation), for every c and w. Proof: far pairs act only as
   *guarded unary lists* (Lemma G), Sigma_0 is closed under all 729 uniform majorities (Theorem A of wave 7, re-derived here at graph level), and a Helly lemma for lists of size <= 2 on a majority-closed relation over a 3-element domain gives 3 coordinates (Lemma H).
3. **Core Bound in general: FALSE** (Theorem CE). Persistent-bit chain, c = 2, tw(I(H)) = 2, pair-maximal, l = 0: for every N a *minimal* certificate T_N of **2N+4** far pairs with F'(T_N) = |H|+1
   and F'(T_N minus e) = |H| for all e. Explicit family, proof for all N by a three-line transfer argument, and VERIFIED with two independent exact oracles (full minimality check for N <= 4, F'(T) and the no-unnecessary-pair check up to N = 10), with the partition-variant F of wave 8 (N <= 4), and with
   the tree-decomposition locality of wave 8 (the minimal certificate then has 3N+3 pairs, N <= 4). The mechanism: a **two-element guard** {L,R} ("pair ac apart") switches off an escape state M of a soft-triangle block, i.e. it toggles an implication. This is exactly the
   "multi-valued unary guard" case that wave 8 could not bound; singleton guards stay bounded (kappa_1 = 2 in every chain tested).
4. **What survives**: the *smallest* certificate is not unbounded in my family (it is 3, because a far hard edge at the centre of a member always gives a 3-pair certificate at l = 0; no pendant design avoids it, Section 6).
   So B1 with a *cleverly chosen* certificate is not refuted; B1 with an arbitrary minimal certificate (as stated in wave 8) is. A bound on the minimum certificate (or a different branching object) is OPEN.
5. **Polymorphisms** (Section 4): the c = 2 language with soft triangles (163 K3-block relations, 445 diamond relations, 39 link relations, lists) has a conservative majority (unique: all-distinct -> M) and a conservative NU of arity 4, but **no** conservative NU of arity 4,5,6 with the plurality property
   needed for 2-element lists; the scope language with soft triangles has **no conservative NU of any arity** (PROVED, Lemma R5); the octahedral K4 blocks have no conservative majority/NU4/NU5. So the NU route cannot give the bound, consistent with Theorem CE.
6. **Real hardness: not obtained.** The chain is a persistent-bit chain only for the *relaxation* (partial guards drop triples). In the real problem the same instance has OPT = |H|+N (every stub must be cut, the anchors are centres), and
   singleton locks (what a clean tunnel gives) have core <= 2 on it. Searches for a chain with unbounded *singleton* core (annealing up to 9 members, octahedra, ladders, strips, K3/diamond chains) never exceed 3. So the W[1]-hardness blueprint of wave 6 is still without its gadget. HEURISTIC/OPEN.

## 1. Conventions and the certificate used here

Member (a,b,c): centre b, edges ab, bc, non-edge ac; states L (edit ab), R (edit bc), M (edit ac); Sigma = zero-excess state vectors. M(v) = members through v.
Pair {u,v} is **local** iff u,v in V(H) and some member through u equals or shares a vertex with some member through v; otherwise **far**; far hard pairs: N. (Soft pairs are local.)
For T subset N the **csp-variant** F'(T): minimum, over assignments of a status to every soft pair (a non-empty subset of its member's pairs: 7 states per member) and to every local hard pair (edited or not),
of the number of edits, such that every vertex triple *all of whose far pairs lie in T* is not an induced P3 (far pairs in T keep their G-status, triples containing a far pair outside T are dropped).
Lemma R of wave 8 holds verbatim (F' monotone; S solution with S cap T = empty gives F'(T) <= |S|; F'(N) = least size of a solution without far edit; the last because with T = N nothing is dropped).
Hence **Theorem B1 holds with F'** (branching on a minimal T with F'(T) > |H|+l; depth <= l, branching factor |T|); checked: `w9b_branch.py` agrees with exact OPT on 449 (instance,budget) pairs, budgets 0..3, n <= 12, three deletion orders.
Pair-maximality is not used in B1. The csp-variant (not the partition variant) is the right oracle because it is FPT for arbitrary T (Theorem L); the partition variant needs width growing with |T|.

## 2. The leaf test (PROVED)

### 2.1 The variables
Member i: domain the 7 non-empty subsets of {ab,bc,ac} (cost = size). Local hard pair uv: bit (cost 1). Far pairs: constants (in T) or absent. Triple constraint: not exactly two of the three pairs are edges after editing.

### 2.2 Theorem L
Given a tree decomposition of I(H) of width w (so c <= w+1), F'(T) for any T subset N, in particular the leaf test F'(N), is computable in time f(w) n^{O(1)}.

**Lemma (triple scopes).** Let tau = {x,y,z}, P(tau) its non-far pairs, W(tau) = union of M(u) over the ends u of pairs of P(tau). If P(tau) is non-empty, W(tau) lies in a ball of radius 3 of I(H) around one member.
Proof. P has one, two or three pairs. One pair uv: members A in M(u), B in M(v) with d(A,B) <= 1; every member of M(u) is within 1 of A (it contains u), every member of M(v) within 1 of B, so all within 2 of A. Two pairs uv, vw (they share v): A in M(u), B in M(v), B' in M(v), C in M(w) with d(A,B) <= 1, d(B',C) <= 1 and d(B,B') <= 1; around B: A within 1, B' within 1, C within 2, all members of M(w) within 3, all of M(u) within 2. Three pairs: the two-pair argument applied to two of them. QED.
(This is where the earlier "exceptional triple" lives: two local pairs {a,y},{y,b} and one far pair {a,b}; y is within distance 1 of both a and b, so a and b are within distance 3 and W(tau) is covered by radius 3 around any member through y. It is only *not* bag-local if local is defined through X_a cap X_y of a decomposition of I(H) itself.)

*Proof of Theorem L.* Let (T,{beta_t}) be the given decomposition and B_t = N^3[beta_t] (all members within distance 3 of beta_t). Since the sets {t : m in B_t} are subtrees (union of subtrees T_y, y within distance 3 of m, chained by the edges of I(H)), the members form a tree decomposition of I(H)^4;
|B_t| <= (w+1)(1+D+D^2+D^3) with D = 3(c-1) <= 3w. A local pair bit e(uv) lives in the bags with M(u) union M(v) subset B_t: this set is a non-empty (take t with beta_t containing a member of M(u); M(v) is within distance 3) intersection of subtrees, hence a subtree. By the Lemma every triple constraint has all its variables (the owner members of its soft pairs, the bits of its local pairs; far pairs are constants) in B_t for any t with beta_t containing the centre member.
The variables per bag are at most |B_t| members and (3|B_t|)^2 bits; min-sum dynamic programming (7^{|B_t|} 2^{(3|B_t|)^2} per bag) over the decomposition solves the CSP; the n^3 triples are assigned to bags. QED
The same CSP with some far pairs dropped is the oracle F'(T); no bound on |T| is needed.

### 2.3 Verification (`w9b_leaf.py`, `w9b_leaf_test.py`, `w9b_locality.py`)
* CSP (bucket elimination, min-fill order) against CP-SAT over pair variables with the far pairs forced: 5 random families (planted/random, n = 9..13, c <= 3), **280 instances, 0 disagreements**.
* Radius check over all triples with a non-far pair, random tree-like instances (m = 12, 16, 14; c <= 3, 2, 4; 150, 150, 100 instances): the minimal centre radius is always <= 3, and 3 occurs (46/73, 88/93, 31/46 instances), so radius 3 is needed and sufficient.
* B1 with this oracle against exact OPT: 449 (instance, budget) pairs, 0 mismatches.

## 3. The Core Bound at l = 0 without soft triangles (PROVED)

**Lemma G (guards).** Let l = 0 (so members are in L,R,M, local hard pairs are kept) and let tau be a triple containing a far pair. Then
(i) if all three pairs are far, or the non-far pairs are hard, tau is a *static* test, violated iff tau is a P3 of G, active iff its far pairs lie in T (at most 3 far pairs);
(ii) if exactly one pair xz is non-far and soft (owner M) the constraint is a *list* on M, active iff the two far pairs lie in T: x_xz(s_M) = 1 (if both far pairs are edges) or = 0 (exactly one is an edge), i.e. one of {R,M},{L,M},{M} / {L},{R},{L,R};
(iii) if two pairs uv, vw are non-far they cannot both be soft (the owners share v, so xz would be local); one soft and one local hard gives a list on the owner active iff the single far pair uw lies in T; two local hard pairs give a static test.
Triples with no far pair give Sigma_0 = the ZE set of the instance with all far pairs dropped. (Proof: read off the definitions.)
If H has no soft triangle then, by Lemma 1 of wave 7, Sigma_0 is the intersection of lists and of binary relations R_AB (one per pair of members sharing a vertex v, the conjunction of the four cross triples (p,v,r); all four cross pairs are local), nothing else.

**Theorem A (wave 7), re-derived.** Every one of the 729 uniform majorities m: {L,R,M}^3 -> {L,R,M} preserves every R_AB and every list. My own enumeration (`w9b_poly.language('sf')`: from the formula, 40 distinct binary relations, as in wave 7) and an independent graph-level computation (`w9b_thmA_graph.py`: all 6x6 orientations x 16 cross statuses of two members sharing a vertex, 32 relations) both give **0 violations over all 729 majorities** (`w9b_thmA.py`, `w9b_thmA_graph.py`).

**Lemma H (lists, domain of size 3).** Let Sigma subset D^V, |D| = 3, be closed under all uniform majorities. If lists L_i (one per coordinate i in a finite set I, |L_i| <= 2) satisfy: Sigma has no point in prod L_i, but for every i a point x^i with x^i_l in L_l for l != i, then |I| <= 3.
Proof. Suppose |I| >= 4. Then x^i_i is not in L_i. For i,j,k in I and the point m(x^i,x^j,x^k) (m a majority applied coordinatewise, it is in Sigma): at a coordinate l outside {i,j,k} the three values lie in L_l, so the majority value does (two equal: that value; all distinct: L_l = D). At i the values (x^i_i, x^j_i, x^k_i): if x^j_i = x^k_i the output is that value, in L_i. Otherwise L_i = {x^j_i,x^k_i} has exactly two elements and (x^i_i,x^j_i,x^k_i) is a permutation pi of D; the output is in L_i iff m(pi) differs from its *first* entry. Likewise j forbids the second entry and k the third entry of the permutation at that coordinate. Only the six permutations are free, so m exists unless the three coordinates all carry the same permutation, i.e. L_i,L_j,L_k are 2-lists whose missing values are pairwise distinct.
Four coordinates cannot be pairwise "distinct-missing" (3 values), so some triple has a free m, and the resulting point lies in Sigma and satisfies every list, a contradiction. QED (This is the 2-list case that the Baker-Pixley argument for singleton locks leaves open.)

**Theorem CB0.** Let H be a packing without soft triangle, l = 0. Every minimal T subset N with F'(T) > |H| has |T| <= 12 (<= 3 if T is a static violation), for all c and w.
Proof. Static violations: at most 3 far pairs. Otherwise the active guards of T are lists (Lemma G), each from at most two far pairs. Let G' be a minimal inconsistent family of guards. Merge guards on one coordinate (their intersection is a list); a minimal inconsistent subfamily of the merged lists has at most 3 coordinates (Lemma H, Sigma_0 closed under all majorities), and by minimality G' uses only those coordinates. At one coordinate the partial intersections of the guards of G' strictly shrink (3 -> <=2 -> <=1 -> 0): if some coordinate reaches the empty list the family is these <= 3 guards, else <= 2 guards per coordinate, so |G'| <= 6. The far pairs of G' form an inconsistent subset of T, equal to T by minimality: |T| <= 2*6. QED

VERIFIED-BY-COMPUTER (evidence): soft-triangle-free random tight instances (c <= 3, m = 7, 400 and c <= 2, m = 8, 386 instances, far pairs dropped): maximum minimal inconsistent list family (singleton lists / arbitrary lists) = 2 in all 786; minimal l = 0 certificates by deletion with the oracle of Section 2 (`w9b_cert0.py`): 117 tight instances (m = 5, c <= 3: 36 instances, sizes {0:2,1:4,2:12,3:16,4:2}; m = 6, c <= 2: 81 instances, sizes {0:4,1:6,2:43,3:25,4:3}), never above 4.
Random instances *with* soft triangles are almost never tight (5 tight in the run of `w9b_cert0.py`, sizes 2,3), so the random experiments say nothing about Theorem CE; the chain family is the evidence.
Status for l > 0 (soft-triangle-free): OPEN; see Section 8.

## 4. Polymorphisms (Sections 4.1 to 4.4: computer enumeration; 4.5: proved)

All searches: conservative idempotent operations (lists of size 1 and 2 are in the language because every list is a guard), SAT with counterexample-guided clause generation (`w9b_poly.py`, `w9b_nu2.py`), exhaustive check of every candidate.
4.1 *Soft-triangle-free* language (40 binary relations + 6 lists): conservative majority exists; all 729 uniform majorities preserve it.
4.2 *Free single-triple language* (every vertex triple "not exactly two pairs together" with the roles/constants as relations): no conservative NU of arity 3,4,5,6 already for the binary single-triple relations (two relations suffice, `w9b_polymin.py`).
4.3 *Triangle scope language* (relations of all graph-level configurations of 1..3 mutually adjacent members, a cross pair marked "dropped" if it is soft in an outside member: 115 binary + 8069 ternary relations, `w9b_scope.py`, `w9b_scope_run.py`): no conservative NU of arity 3,4,5,6 (`w9b_nu2.py`). 3,839 of the 8,068 ternary relations kill the majority individually and each of them needs a dropped pair (c >= 3 at the vertex).
4.4 *c = 2 language* (blocks of the triangle-connected components of a subcubic I(H): K3 (163 relations), diamond = K4 minus an edge (445 relations, all 4 members, `w9b_diamond2.py`), the 39 link relations, lists; K4 is an isolated component at c = 2 and is treated separately): conservative majority exists and is **unique** (all-distinct -> M, i.e. the "M-majority" of wave 6); a conservative NU of arity 4 exists; **no conservative NU of arity 4,5,6 whose value always has multiplicity >= 2** (the property that would make the Helly argument work for 2-element lists, `w9b_nuP.py c2`); no conservative Mal'tsev (`w9b_maltsev.py`; edge terms not tested). A relation closed under the M-majority alone can have a list core of size >= 6 (`w9b_helly.py 212100`, Baker-Pixley closure, random search), so the single available majority proves nothing about 2-element guards. Octahedral K4 blocks (521 relations): no conservative majority, NU4, NU5, Mal'tsev; 42 of them kill the majority alone; every K4 relation has kappa_1 <= 2, kappa_2 <= 3 (`w9b_k4*.py`).
4.5 **Lemma R5 (PROVED).** The language of Section 4.3 has no conservative NU polymorphism of any arity r >= 3.
Proof. It contains the triangle-scope relation Rk = {000,002,011,020,101,110,111,121,200,220} (states 0,1,2 = L,R,M; realised by the configuration printed by `w9b_r5proof.py`) and the lists. Rk intersected with {L,R}^3 is R5 = {x+y+z != 1}. A conservative f maps {L,R}^r into {L,R} and preserves R5 on {L,R}. Let N_0 = {B subset [r] : f(1_B) = 0}. For B,C in N_0 the three columns (B union C, B, C) cover every coordinate 0 or at least 2 times, so they are rows-compatible with R5 and (f(B union C),0,0) in R5 forces f(B union C) = 0. NU with r >= 3 gives f(e_i) = 0 for all i, so every singleton is in N_0 and N_0 is union-closed: f(1^r) = 0, contradicting idempotence. QED (checked for r <= 8 by SAT, `w9b_bool.py`.)
Conclusion for the question asked: **no bounded-arity NU exists for the full language with soft triangles**, and at c = 2 the NU that exists cannot bound 2-element guards; this is exactly what Theorem CE realises.

## 5. The persistent-bit chain: Core Bound is false (PROVED + VERIFIED)

### 5.1 The family G_N (`w9b_chain_cert.build`)
Vertices per block j = 1..N: p_j,q_j,r_j,c_j,b_j; a_1 extra; a_j := b_{j-1} (j >= 2). Members (centre in the middle)
A_j = (q_j, p_j, a_j), B_j = (q_j, r_j, b_j), C_j = (p_j, c_j, r_j).
Hard pairs: non-edge except the hard edges c_j q_j (j = 1..N). Pendants (vertices of no member, all their pairs far):
w_j (j = 1..N) adjacent to p_j only; w_{N+1} adjacent to a_1 only; w_{N+2} adjacent to b_N only.
|V| = 6N+3, |H| = 3N, H is a pair-maximal packing (the hard graph is a matching, so no P3 has three hard pairs), each vertex lies in at most 2 members (c = 2), I(H) is a chain of triangles joined by single edges (tw 2, N soft triangles).
Far pairs of the certificate: T_N = {w_j p_j, w_j r_j (j <= N); w_{N+1} a_1, w_{N+1} p_1; w_{N+2} b_N, w_{N+2} r_N}, |T_N| = 2N+4.

### 5.2 Structure lemma (PROVED)
With all far pairs dropped, Sigma_0(G_N) = { (s_1,..,s_N) : s_j in {LRL, RLR, RRM}, not (B_j = L and A_{j+1} = L) } where s_j lists the states of (A_j,B_j,C_j) and (B_j = L) iff s_j = RLR, (A_{j+1} = L) iff s_{j+1} = LRL.
Proof sketch: inside a block the triples give the three states (K3 scope relation {LRL,RLR,RRM}, `w9b_k3chain.py`, direct enumeration); the only triples with vertices of two blocks that carry soft pairs have both soft pairs through the glue vertex b_j = a_{j+1} and a hard non-edge cross pair, which is the single forbidden combination; the other mixed triples have one soft pair and two non-edges. Machine check: |Sigma_0| = 3, 8, 21, 55, 144, 377, 987 for N = 1..7 equals the transfer-matrix count (F_{2N+2}), `w9b_k3chain_core.py`.
The pendant guards active for T_N (a triple with its two far pairs in T_N): the triple {w_j,p_j,r_j} gives "pair p_j r_j apart", i.e. x_ac(C_j) = 0, i.e. C_j in {L,R} (not M); {w_{N+1},a_1,p_1} gives A_1 = R; {w_{N+2},b_N,r_N} gives B_N = R. No other triple is active.

### 5.3 Theorem CE (PROVED)
F'(T_N) >= |H|+1 and F'(T_N minus e) = |H| for every e, F'(empty) = |H|. So T_N is a minimal certificate of size 2N+4 at l = 0 (c = 2, w = 2, pair-maximal).
Proof. All C_j in {L,R} and A_1 = R: block 1 cannot be LRL (A = L) nor RRM (C = M), so it is RLR, B_1 = L; the link forbids A_2 = L, so block 2 is again RLR; inductively B_N = L, contradicting B_N = R. Minimality: without the A_1-guard all blocks LRL; without the B_N-guard all RLR; without the C_k-guard blocks < k are RLR, block k is RRM (B_k = R, A_{k+1} free), blocks > k are LRL. Guards are active only in pairs, so removing any of the 2N+4 far pairs deactivates exactly one guard. QED

### 5.4 Verification (all exact)
* `w9b_chain_cert.py 4`: bucket elimination and CP-SAT (`opt_csp_cpsat`) agree: F'(T_N)-|H| = 1, F'(T_N minus e)-|H| = 0 for all e, F'(empty)-|H| = 0, N = 1..4. `w9b_chain_cert2.py`: N = 6, 8, 10 (|T| = 16, 20, 24), no unnecessary element, CP-SAT confirms F'(T)-|H| = 1; elimination width stays 5.
* Partition variant of wave 8 (far pairs outside T free, transitivity kept; `w9b_chain_part.py`): same minimal certificate for N <= 4.
* Locality of wave 8 (X_u cap X_v of the path decomposition of I(H), `w9b_xlocal.py`): the link cross pairs become far, and the greedy-minimal certificate is T_N plus N-1 link pairs: 9, 12, 15 pairs for N = 2,3,4 (3N+3).
* Exhaustive list-core DFS of Sigma_0 (`w9b_k3chain.py`, `max_core`): largest minimal inconsistent list family has N+2 members (N = 3,4,5: 5,6,7, family = A_1, B_N and all C_j), kappa_1 = 2.
* Search (`w9b_k3search.py`): all K3 block configurations x choice of out/in member x 16 link statuses (1,992 feasible chains evaluated): 48 with growth, all with |Sigma_0(N=3)| = 21 and the same mechanism (4 relation types and their relabellings). Diamond chains (`w9b_diachain.py`, 448 evaluated): no growth. The grammar (`w9b_grammar3.py K3`) rediscovers the same block with growth +2 per round.

### 5.5 Why Sigma_0(G_N) escapes the majority bound
Sigma_0 contains soft triangles, so Theorem A does not apply; it is closed under the M-majority only (Section 4.4) and the guard lists have two elements. Singleton guards are bounded (kappa_1 = 2), 2-element guards are not.

## 6. Minimum certificates, real optimum, and what is still open

* **Minimum certificate = 3** for G_N (N <= 7, exact search over all guard sets, `w9b_chain_min.py`): a pendant w adjacent to p_j (the centre of A_j) with the non-edges w q_j and w a_j respected forces A_j's two pairs at its centre apart, which is not a ZE state. Lemma (cheap certificate): at l = 0 a far hard edge wx with x the centre of a member (y,x,y') and wy, wy' non-edges gives the certificate {wx,wy,wy'}. In the chain every anchor of the {L,R}-guard is such a centre (`w9b_anchor_roles.py` for all 48 variants), and no adjacency set of the pendant inside its block makes all other guards vacuous (`w9b_vacuous.py`). So Theorem CE refutes "every minimal certificate is bounded", not "some certificate is bounded"; a corrected scheme must select a small certificate (OPEN).
* **Real optimum.** OPT(G_N) = |H|+N (CP-SAT, N <= 5): in the real problem every pendant edge forces its anchor to be singleton, impossible for a centre, so all N stubs are cut. The persistent bit exists only in the relaxation (partial guards), not as a lock chain.
* **Real lock chain (wave 6 PFG): not found.** Singleton guards: kappa_1 <= 2 in K3 chains (N <= 5), octahedron chains and trees (~600 feasible instances, `w9b_octchain.py`, `w9b_octtree2.py`), diamond/triangle ladders and lattice strips (`w9b_struct_run.py`, `w9b_strip.py`), and simulated annealing over general instances with c <= 3, 4, n <= 15 and up to 9 members (`w9b_anneal.py`, 25,000 steps x 4 runs, plus short runs): kappa_1 = 3 is easy to reach (5-member gadgets, e.g. G* of wave 6), kappa_1 = 4 never; annealing over unit cells repeated along a chain (`w9b_cellanneal.py`, cells of 4 and 5 members, N = 2,3): kappa_1 stays 2. The three-option centre of wave 8 Sec. 5 is the missing object and was not found.

## 7. Does this give W[1]-hardness? (HEURISTIC / OPEN)
The Grid Tiling blueprint of wave 6 needs lock chains (singleton anchors) with core size N and equality taps. Section 6 shows what we have instead: a relaxation-chain with guards of two elements, which kills the *branching scheme* with arbitrary certificates but does not by itself change the complexity of CEaMP. Candidate next steps: (a) a gadget with kappa_1 = 4 and free ports; (b) replace pendant guards by a mechanism whose unwanted guards are vacuous (not possible for these anchors, Section 6).

## 8. Remarks on l > 0 (HEURISTIC)
Soft-triangle-free, l > 0: for a fixed pattern pi of at most l defects (non-ZE members, local edits) the remaining problem is a ZE problem on the other members with permanent extra lists, so Lemma H and Theorem CB0 give <= 12 far pairs per pattern; patterns are not bounded in number (n^l), and one needs an Erdos-Posa step: the sets of members whose defect could repair a given 12-pair core are connected in I(H), subtrees of a width-w decomposition have transversal number <= (w+1) times the packing number, which would give g(l,w) <= (l+1)(w+1) 12. Not proved (the regions are not literally subtrees; a defect member also changes lists on its neighbours).
With soft triangles l > 0 is moot: g is already unbounded at l = 0 (Theorem CE); two disjoint copies of G_N give certificates of size 2(2N+4) at l = 1.

## 9. Status table

| item | status |
|---|---|
| Leaf test with local edits: radius-3 locality, bounded-treewidth CSP (Theorem L) | PROVED; VERIFIED (280 / 449) |
| Theorem B1 with the csp-variant certificate, any packing | PROVED (as wave 8); VERIFIED (449 pairs, 0 mismatches) |
| Core Bound l = 0, no soft triangle, <= 12 (Theorem CB0, Lemmas G, H) | PROVED; VERIFIED |
| Theorem A (majority closure, no soft triangle) re-derived | VERIFIED (729 x 40 relations; 729 x 32 graph relations) |
| No conservative NU of any arity for the scope language with soft triangles (Lemma R5) | PROVED |
| c = 2 soft-triangle language: unique conservative majority, NU4 exists, no plurality-NU of arity 4,5,6, no Mal'tsev | VERIFIED-BY-COMPUTER |
| K4 octahedra: no conservative majority/NU4/NU5; kappa_1 <= 2 | VERIFIED-BY-COMPUTER |
| Core Bound false: minimal certificates of size 2N+4 (c = 2, w = 2, l = 0) | PROVED (Theorem CE); VERIFIED (N <= 10, both oracles, two locality definitions, partition variant) |
| Minimum certificate bounded? | OPEN (= 3 in the family) |
| Real lock chain (PFG) / W[1]-hardness | OPEN (kappa_1 <= 3 in all searches) |
| Core Bound soft-triangle-free, l > 0 | OPEN (Erdos-Posa sketch) |

## 10. Reproduction (all in `code/`, run from that directory)
Language and polymorphisms: `python3 w9b_scope_run.py` (writes /tmp/w9b_scopes.pkl), `w9b_nu2.py 4 all`, `w9b_nu.py 3 sf`, `w9b_thmA.py`, `w9b_thmA_graph.py`, `w9b_k4.py` then `w9b_k4nu.py`, `w9b_diamond2.py` then `w9b_diamond_an.py`, `w9b_c2lang.py`, `w9b_c2maj.py`, `w9b_nuP.py c2`, `w9b_maltsev.py`, `w9b_r5proof.py`, `w9b_bool.py`.
Leaf test and B1: `w9b_leaf_test.py`, `w9b_locality.py 1 150 12 3`, `w9b_branch.py 1 60` (also seeds 2, 3 with 150), `w9b_cert0.py 1 4000 5 3 240`, `w9b_sf_core.py 1 6000 7 3 sf`.
Chain: `w9b_k3chain.py 3`, `w9b_k3chain_core.py`, `w9b_chain_cert.py 4`, `w9b_chain_cert2.py`, `w9b_chain_part.py 4`, `w9b_xlocal.py 4`, `w9b_chain_min.py 7`, `w9b_vacuous.py 4`, `w9b_k3search.py` (about 90 s), `w9b_k3fam.py`, `w9b_anchor_roles.py`, `w9b_diachain.py`, `w9b_grammar3.py K3 4`.
Searches for lock chains: `w9b_anneal.py seed n m steps cmax`, `w9b_octchain.py`, `w9b_octtree2.py`, `w9b_struct_run.py`, `w9b_strip.py`, `w9b_cellanneal.py`.
Libraries: `w9b_lib.py` (instances, locality, `opt_cpsat`, `opt_csp_cpsat`), `w9b_sigma.py` (Sigma, `max_core`), `w9b_poly.py`, `w9b_scope.py`.
