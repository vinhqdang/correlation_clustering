# Wave 7: is CEaMP FPT in (l, w) in the smallest regime? (w <= 1, or c <= 2)

Task: settle, by a proof in either direction, the regime (a) w <= 1 (I(H) a forest, every bag <= 2 members) or (b) c <= 2 with l arbitrary,
using that Proposition FAR and the family N*(F) live at w = 1, c = 2. No git, paper/ untouched.
Code: `code/w7_*.py` (library `w7_lib.py`, written from scratch; subset-DP in C `w7_dp.c` + `w7_dpy.py`), all reproducible, commands in Section 9.

Labels. PROVED = proof given here (finite case checks that are part of a proof are marked "classification by computer", the script is named).
VERIFIED-BY-COMPUTER = exhaustive or sampled search with the named script, evidence only. HEURISTIC = argument without full proof. OPEN = not decided.
Notation: soft = packed pair, hard = unpacked; member (a,b,c) has centre b, edges ab, bc, non-edge ac; states L = edit ab, R = edit bc, M = edit ac;
Sigma(G,H) = set of zero-excess (ZE) state vectors in {L,R,M}^H (ZE solution <=> one edit per member and nothing else, Lemma decomp);
x_e(s) = 1 iff pair e of a member is together after its single edit: x_ab(L,R,M) = (0,1,1), x_bc = (1,0,1), x_ac = (0,0,1).
A **soft triangle** is a vertex triple whose three pairs are packed by three *distinct* members. For c <= 2 a soft triangle is exactly a triangle of I(H);
for every w <= 1 (I(H) a forest) there is none.

## 0. Verdict

1. **The question is not settled.** I can prove neither W[1]-hardness nor FPT of CEaMP for w <= 1 or for c <= 2. Both are still OPEN.
2. What is PROVED, and is new relative to waves 5-6:
   * **Theorem A (majority closure).** If (G,H) has no soft triangle (in particular if I(H) is a forest, i.e. w <= 1, or if c <= 2 and I(H) is triangle-free,
     *whatever the treewidth*; no maximality of H is needed), then Sigma(G,H) is closed under **every** majority operation applied coordinatewise
     (all 729 uniform ones at once; in particular under "all distinct -> M"). Consequences: (A1) l = 0 is decidable in polynomial time by path consistency
     (for any c, any w); (A2) every minimal non-face of clean anchors ("core") has at most 2 anchors, i.e. Conjecture C2 of the paper holds for triangle-free instances
     (and the false uniform-M-majority form of wave 6 is exactly true here); (A3) the tunnel-release model is FPT, time 2^l n^O(1).
     **The only possible source of non-bijunctive lock behaviour at c <= 2 is a soft triangle (a triangle of I(H)).**
   * **Strictness lemma (S).** Across any edge of the interaction forest a lock that really restricts the neighbour leaves at most one state (unit propagation);
     so lock cores are conflict *pairs*, obstructions are implication paths.
   * **Theorem B (abstract hardness, own proof).** The abstract problem "hit all obstruction paths of a forest of paths by <= L operations, each a single node or a *link*
     (a pair of nodes, link graph of maximum degree 1)" is W[1]-hard in L, by a reduction from Grid Tiling with <= (L = k^2, tree = 2k disjoint paths).
     This is exactly the 'far two-sided edit shares cost' mechanism in an abstract setting, with a unit gain per link.
   * **Theorem C (precise obstruction to realising Theorem B in CEaMP at w <= 1).** The reduction needs link operations that hit *interior* nodes of obstruction paths:
     one tunnel cut at any of the n positions of an arm must repair that arm. In CEaMP that is a core of size n; Theorem A excludes cores of size >= 3 when there is no soft triangle.
     A tunnel cut only ever repairs conflicts *at its two ends*.
3. Sharpest positive statement for the unrestricted problem: the abstract model that CEaMP at w <= 1 reduces to heuristically (endpoint-hitting links + interior-hitting single breaks, "RM2")
   is FPT *assuming a packing-endpoint lemma* that I verified on 15,000+ (instance, maximum packing) pairs but did not prove.
4. Regime (b) c <= 2: triangle-free I(H) of any width is covered by Theorem A. With triangles there is no proof; the K_4 gadget of audit 8 breaks the uniform majority but not 2-decomposability;
   exhaustive/sampled evidence (K_3 blocks with branches, all K_4 gadgets, prisms, 10^4 random instances) finds **no** non-2-decomposable Sigma and **no** core of size 3 at c <= 2.

## 1. Conventions and the exact structure of Sigma

**Lemma 1 (constraint structure) PROVED.** Let H be any packing. A state vector s is in Sigma iff no vertex triple is an induced P3 of G after editing. Classify a triple by the owners of its soft pairs:
(i) no soft pair: a static condition on G; (ii) one soft pair (owner A) and two hard pairs: a *unary* constraint on A; (iii) soft pairs from two distinct members A,B: they share a vertex q of the triple,
and the third pair {p,r} (p in A, r in B) is hard (it is not soft: soft would make it a pair of a third member C, a soft triangle, or of A or B, forcing |A cap B| >= 2); this is a *binary* constraint on (s_A,s_B);
(iv) all soft pairs from one member: the member's own state; (v) three soft pairs from three distinct members: a soft triangle, a *ternary* constraint.
For two members A,B sharing a vertex q (two members share at most one vertex, since their pair sets are disjoint) the binary constraints are the four triples (p,q,r), p in A\{q}, r in B\{q},
with static values eps_pr = [pr in E]: the conjunction R_AB(role of q in A, role of q in B, eps) in {L,R,M}^2, where the triple (p,q,r) is fine iff x_pq(s_A) + x_qr(s_B) + eps_pr != 2.
Hence **for a soft-triangle-free instance Sigma = intersection of unary lists and of the relations R_AB over pairs of members sharing a vertex.** (The cross bits eps also create unary lists on A and on B.)

**Lemma 2 (exact U-only identity; locality of effect) PROVED.** Say a solution is *U-only* if every member gets exactly one edit; its other edits form a set U of unpacked pairs. Then
OPT_U-only - |H| = min{ |U| : Sigma(G xor U, H) != empty } (apply U to G first: the members remain induced P3's because U avoids packed pairs).
Moreover Sigma(G xor U) differs from Sigma(G) only through the lists of members containing an endpoint of a pair of U and through the relations R_AB of pairs (A,B) of *adjacent*
members with p in A, r in B for some pair {p,r} in U. In particular a *far* pair of U (ends in non-adjacent members) acts only by changing a unary list at each of its two ends.
(Proof: the status of a hard pair only enters the constraints of Lemma 1 in the triples that contain it, and those involve only the members through its endpoints.)
VERIFIED-BY-COMPUTER (`w7_uonly.py`, 120 soft-triangle-free perturbed instances, n <= 12, exact subset DP): OPT-|H| <= U-only value always; equal in 119/120; one instance has OPT-|H| = 1 with U-only value 2,
so double-edited members are really needed to be exact (the U-only model is a sub-model, not the problem).

## 2. Theorem A: majority closure for soft-triangle-free instances

**Theorem A (PROVED; classification of the relations by computer, closure checked by hand per class and by computer).**
Let (G,H) have no soft triangle. Then Sigma(G,H) is closed under every majority operation m: {L,R,M}^3 -> {L,R,M} (m(a,a,b)=m(a,b,a)=m(b,a,a)=a), applied coordinatewise with the same m on all coordinates.

*Proof.* By Lemma 1, Sigma is an intersection of unary lists and binary relations R_AB; it suffices that each is closed under every majority. A unary list is a subset S of {L,R,M}; a triple of its elements
has an all-distinct coordinate only if S = {L,R,M}, which is closed. For the binary relations: R_AB depends on (role of q in A, role of q in B) in {a,b,c}^2 and eps in {0,1}^4: 144 configurations.
`w7_hub3.py`/`w7_rel6.py` enumerate them (the relation is computed from the formula of Lemma 1) and find **40 distinct relations** on {L,R,M}^2, of the following shapes only
(sizes 0,1,2,3,5: 1, 9, 18, 8, 4 relations):
 * empty, singletons and 2-element sets: closed, since among three tuples of a 2-element set only two distinct values occur in every coordinate, so m(u,v,w) is the coordinatewise majority;
 * 3-element *lines* (a full row {(p,*)} or full column {(*,q)} with p,q in {L,R}): one coordinate constant, the other arbitrary and the result is again on the line;
 * 3-element subsets of {L,R}^2 (four of them): only two values per coordinate, so m is Boolean majority, and every binary Boolean relation is closed under majority;
 * 5-element *crosses* (row p union column q, p,q in {L,R}): if the first coordinate of m(u,v,w) differs from p, either two of the tuples have that first coordinate (then they are (x,q) tuples and the second coordinate is q) or the three first coordinates are all distinct
   (then two tuples have first coordinate != p, hence are (x,q) and the second coordinate is q by majority); symmetric for columns.
Independent machine check of the closure for all 729 rules and all 96 (role pair, relation) pairs: `w7_rel5.py` (output: every rule is compatible with every relation). Independent re-derivation on real graphs, without the formula
of Lemma 1: `w7_rel7.py` builds all 576 five-vertex graphs with two members sharing a vertex and any pattern of the four cross pairs, computes Sigma by the cluster-graph test, and checks closure under **all** majorities at once
(|Sigma| in {0,1,2,3,5}: histogram 260/228/56/16/16): 0 failures. Instance-level: `w7_t4.py`, 1,568 random soft-triangle-free tight instances (c up to 4, n up to 22): 1,568 closed under all majorities, 1,568 2-decomposable. QED

*Remark (why the conjunction matters).* The 18 single triple types (e_A, e_B, eps) are **not** simultaneously closed under any uniform majority (`w7_rel2.py`: 0 of 729 rules), and
the per-coordinate rules that work for one type are incompatible for a member that uses all three of its pairs (`w7_rel4.py`: intersection of admissible rules over the pairs ab, bc, ac is empty). The four triples at one shared vertex together
always give a small, rigid relation. This is why wave 6's "all-distinct -> M" held empirically for c <= 2 and why the K_4 counterexample (audit 8) needs a soft triangle.

**Corollary A1 (l = 0 is polynomial without soft triangles) PROVED modulo the cited theorem of Jeavons-Cohen-Cooper (1998) that a language with a majority polymorphism is decided by path consistency.**
For soft-triangle-free (G,H) of arbitrary c and tw(I(H)), "is OPT = |H|?" is a binary CSP on |H| variables of domain 3 with unary lists and majority-closed binary relations; path consistency decides it in O(|H|^3) per pass.
VERIFIED-BY-COMPUTER (`w7_pc.py`, independent solver built from the vertex triples vs. brute-force enumeration of Sigma, 1,500 instances with 715 feasible and 785 infeasible): 1,500/1,500 agree.
(This is consistent with the NP-hardness of l = 0 for c <= 23 of Li-Pilipczuk-Sorge: their gadgets must contain soft triangles.)

**Corollary A2 (cores <= 2) PROVED.** For soft-triangle-free (G,H), every minimal non-face of clean anchors has at most two anchors (Lemma "majority bounds cores" of the paper, which needs a uniform majority); likewise every set of unit locks is infeasible iff two of them are
(and all projections of Sigma are 2-decomposable, Baker-Pixley). So Conjecture C2 holds when c <= 2 and I(H) is triangle-free. VERIFIED-BY-COMPUTER: no core of size 3 in 4,000 random instances (`w7_t4.py`, both classes).

**Corollary A3 (tunnel-release model is FPT in this class) PROVED.** Model TR of the paper (tunnels K_2 or stubs on clean anchors, every member singly edited, hard edits = tunnel cuts), I(H) without soft triangles. Let C be the conflict graph on the anchors
(edge xy iff {x,y} is not a face, loop iff x alone is not; computed by 2-anchor path-consistency calls). The released sets R making the instance tight are exactly the vertex covers of C (Prop. tTR (ii) + A2), so the optimum is the minimum number of tunnels whose
endpoints cover C: pick an uncovered edge of C and branch on the tunnel at either end: 2^l n^O(1) for matching tunnels (fan-out d: (2d)^l). Nothing about w or c is needed beyond Theorem A.
Illustration, VERIFIED-BY-COMPUTER (`w7_pathcores.py`): on the path P_{2N+1} with its N consecutive members, the lockable anchors are the N+1 even vertices and **every pair conflicts** (conflict graph K_{N+1}); a chain of members has to be released in all but one anchor.

## 3. The strictness lemma: locks are unit propagation

**Lemma S (PROVED by finite check; `w7_strict.py`).** For an edge between adjacent members (parent P with root vertex role q_m, child C with role q_n, cross bits eps), let LamP, LamC be the unary lists the four cross bits impose and R the relation of Lemma 1.
For all subsets X <= Y of {L,R,M} put X' = LamP cap img_R(X cap LamC), Y' likewise. If X' != Y' then |X'| <= 1. (144 configurations x 27 pairs X <= Y; 803 strict cases, 0 violations of |X'| <= 1.)
*Consequence.* A set of members whose unary lists are restricted (locks, also 'forked' locks whose list has two elements) acts on every other member only through forced single states (or a contradiction): propagation is unit propagation along tree edges.
Locks therefore interact in pairs, the obstruction behind a conflict is an *implication path* between two list-restricted members, and an unlock never creates new conflicts.
VERIFIED-BY-COMPUTER (`w7_helly.py`, fixpoint over all chains with arbitrary extra unary lists, single anchors and forked lists): the family of achievable (locked, unlocked) branch images stabilises after 3 rounds; at a hub the strict
locked images of size >= 2 never occur, hence **no hub 3-core** (a hub with three branches each carrying one lock) exists, even with forked locks. (Caution: if the unary side effects of the cross bits are ignored, the calculus predicts a 3-core with branch images
{L,M},{L,R},{M,R} (`w7_hub2.py`); the side effects of the same bits kill it.) This covers the only possible minimal 3-core shape in a forest (three locks have a median node with one lock per branch, or lie on a path with a singleton lock in the middle).
Statuses: single-lock hub case PROVED/verified; arbitrary families of forked locks with several locks per branch: OPEN as a theorem (evidence only).

## 4. The reduction attempt: from far two-sided edits to Grid Tiling, and where it breaks

### 4.1 What the far two-sided edit gives (VERIFIED-BY-COMPUTER, independent exact DP `w7_far.py`, subset DP in C)
* Cost sharing: members (0,1,2),(3,4,5), pendant hard neighbours 6~2, 7~5 force state R in both; OPT = 2 without tunnel, 4 with private stubs at 0 and 3, **3 with one tunnel 0-3** (one cut serves both members), 2 if the pendants are removed.
* FAR (Proposition FAR): path 0..2m plus chord (1,2m-1): OPT = m+1 for m = 3,4,5,6 (n up to 13), 7/9/11/13 optima, **every optimum deletes the chord**. But the chord is a *cost*, not a gain:
  with 1, 2 or 3 chords between centres (m = 6, all chord sets, `w7_chords.py`) the excess is exactly the number of chords. The 1-unit gain exists only for tunnels at *forced* members (the sharing example).
* N*(F): its weak point is only that the whole component is irreducible; for l = 1 it is solved by the trigger y. Nothing in it shares costs between positions.

### 4.2 The abstract problem and its W[1]-hardness (Theorem B)
**LPPH (linked-pair path hitting).** Input: a forest (here: disjoint paths), a family Omega of 'obstruction' paths, a set of links (pairs of nodes), L. An operation is a single node or a link; it *hits* every path containing one of its nodes.
Question: can L operations hit all paths?

**Grid Tiling with <=** (Marx 2012; Cygan et al., *Parameterized Algorithms*, Sect. 14): given S_ij subset [n]^2 (1 <= i,j <= k), choose s_ij = (a_ij,b_ij) in S_ij with a_ij <= a_{i,j+1} and b_ij <= b_{i+1,j}
(this is the book's version after transposing the grid); W[1]-hard in k.

**Theorem B (PROVED; correctness cross-checked by `w7_phlp.py`, `w7_phlp2.py`).** LPPH is W[1]-hard in L, already for a forest of 2k paths and a link set that is a perfect matching (every node in at most one link).
*Construction.* For each cell (i,j) two arms: A_ij has one node per tile of S_ij, ordered by (a,b); B_ij has one node per tile, ordered by (b,a). Row i is one long path obtained by concatenating the arms A_{i,k}, A_{i,k-1}, ..., A_{i,1} (in this order, each in increasing position);
column j likewise with B_{k,j}, ..., B_{1,j}. Obstruction paths: (1) every arm; (2) for each row i, j < k and alpha in [1,n-1] the path rho = {nodes of A_{i,j+1} with a >= alpha+1} cup {nodes of A_{i,j} with a <= alpha}
(a subpath of the row path: the suffix of A_{i,j+1} is followed directly by the prefix of A_ij); (3) the same for columns with the b-coordinate and B arms (b_ij <= b_{i+1,j}). Links: the node of tile e in A_ij with the node of tile e in B_ij. L = k^2.
*Correctness.* (=>) Given a tiling, use for every cell the link of its tile: the arms are hit; a path rho of type (2) is hit iff a_ij <= alpha or a_{i,j+1} >= alpha+1, and a_ij <= a_{i,j+1} gives this for every alpha. (<=) The 2k^2 arms are pairwise disjoint and an operation hits
at most two of them (a link only the two arms of its own cell), so exactly k^2 links are used, one per cell: a tile in every S_ij. If a_ij > a_{i,j+1}, the path rho for alpha = a_{i,j+1} is hit by no chosen node (the chosen nodes of other cells are on other arms), contradiction; columns likewise.
Size is polynomial and L = k^2, so a f(L) n^O(1) algorithm would give one for Grid Tiling with <= (ETH: no f(k) n^{o(k)}).
*Computer check* (own brute force for both problems): k=2, n=3: 300 instances (278 YES, 22 NO); k=2, n=4: 400 (209/191); k=3, n=3: 150 (8/142); **0 mismatches** in 850 instances.

### 4.3 Why it does not carry over to CEaMP at w <= 1 (Theorem C)
Dictionary of the intended embedding: arm = chain of members that is *not tight* unless some anchor on it is released; links = tunnels (K_2, fan-out 1 so HUB+ and fan-out bounds hold trivially); one cut per cell = budget l = k^2; rho = further obstructions.
**Theorem C (PROVED, for soft-triangle-free instances).** (i) In the tunnel-release model a cut of the tunnel at an anchor x repairs exactly the conflicts (cores) containing x, and all cores have <= 2 anchors; so a cut repairs only conflicts *at its two ends*.
(ii) If an arm with N anchors must be repaired by cutting a tunnel at *any* one of its anchors (as Theorem B needs: N = |S_ij| positions), then the N anchors form a minimal non-face of size N (all locked: infeasible; any N-1 locked and one released: feasible).
By A2 this forces N <= 2. Hence the reduction of 4.2 cannot be realised with tunnels at w <= 1, nor at any width as long as there is no soft triangle.
(iii) Polarity: the conflicts of a locked system are positive clauses 'release x or release y' (vertex-cover polarity); a Grid Tiling constraint forbids a *pair of choices*, a negative clause. Even with unbounded fan-out one could not forbid
pairs of cuts without an OR over a group (a core of group size).
(iv) Using breaks (double-edited members, state O = all three pairs apart, cost 2) for the interior hits: a break is an unpaired single node operation, so the tile coupling (p,q) in S_ij between two arms is lost.
HEURISTIC: the compatibility of the O state with its neighbours is controlled by the same cross bits (O is compatible with a neighbour state t iff x_qr(t) = 0 for every r that has some eps_pr = 1), so a break is not a universal decoupler;
it is a one-state phase flip whose position can be restricted by design. I did not find a way to couple two break positions across tracks at a total cost of 2.
(v) Other candidate operations: inserted pairs (T2) and cuts inside hard components of size >= 3 also act by unary changes at their two ends (Lemma 2) and a hard K_s needs s-1 cuts to isolate one vertex, so they are not better than K_2 tunnels. HEURISTIC, no proof that no combination of operations simulates interior hits.

### 4.4 Bottom line for hardness at w <= 1
Hardness would need a paired operation that repairs *interior* positions of an obstruction. Among tunnel releases (the 'far two-sided' mechanism) this is impossible (Theorem C). The remaining possibilities are
not tunnels: inserted far pairs, interactions of local breaks, soft-triangle gadgets (c = 2 needs a triangle in I(H), hence w >= 2). None is excluded, and I have no construction for any of them.

## 5. The positive direction: what an FPT proof for w <= 1 would look like (RM2)

Under Strictness (Section 3) and Lemma 2, I believe CEaMP at w <= 1 reduces to the following abstract model (HEURISTIC: it ignores that a break changes relations rather than just removing them, and hard-triple closure of U).
**RM2.** A forest; a set Omega of conflicts (pairs of nodes, with path); operations: *break* at any node (cost 1, hits every conflict whose path contains the node) and *tunnel* = link of a bounded-degree graph (cost 1, hits every conflict with an endpoint at one of its two nodes). L operations.
For tunnels F let nu(F) = maximum number of node-disjoint conflict paths among conflicts with no released endpoint = minimum number of breaks needed (tree-path duality, greedy deepest-LCA); cost = |F| + nu(F).
**Packing-endpoint lemma (OPEN; VERIFIED-BY-COMPUTER).** For every maximum packing Pi of node-disjoint conflict paths there is an optimal F all of whose tunnels have an endpoint among the (<= 2 nu) endpoints of Pi.
Test `w7_rm2.py`, `w7_rm2b.py`, `w7_rm2c.py`: 3,000 + 1,500 random forests/conflict sets/tunnel sets, and 2,500 instances against *all* maximum packings (11,109 pairs): 0 failures (tunnel gain occurred in 234 + 237 + 324 of them);
the stronger version 'both endpoints in Pi' fails (9 of 1,500), so single-endpoint tunnels do matter.
**If the lemma holds, RM2 is FPT:** nu > 2L means NO (every operation hits at most 2 of the packed paths); otherwise at most 2 nu <= 4L endpoints, <= 4Ld candidate tunnels, F ranges over the subsets of size <= L of the candidates: (4ed)^L choices, each evaluated by the greedy.
A proof of the lemma by exchange (replace a tunnel with no endpoint in Pi by a break at a node of the conflicts it alone repairs) works when only one end of the tunnel is useful or when the repaired conflicts share a node (Helly property of paths in a tree), and leaves one case open:
two disjoint conflicts repaired only by the tunnel. This is the *only* gap in the abstract model.
What is missing between RM2 and CEaMP at w <= 1 (all OPEN): (L1) that a break behaves as 'hit every conflict path through the node' (relations change; needs the exact compatibility analysis of the O state), (L2) Helly-2 for lock lists of size two in forests with several locks per branch (Section 3, evidence only),
(L3) the closure condition on U for hard triples (cascades, DEL lemma) and the cost of inserted pairs (T2) in the same calculus, (L4) pair-maximality: Lemma MAX adds <= l members, which can create soft triangles.

## 6. Regime (b): c <= 2

* Triangle-free I(H) (any width, subcubic): Theorem A and its corollaries apply (no soft triangle). In particular c <= 2 and tw(I(H)) unbounded is no obstruction to cores <= 2, to the polynomial l = 0 test, and to FPT of the tunnel-release model.
* With triangles: Sigma can fail to be closed under the uniform majority (K_4 example, audit 8). Evidence on whether lock cores >= 3 or non-2-decomposable Sigma exist at c <= 2:
  * all tight K_3 blocks (6 vertices, 3 members, all role orders and all 6 hard bits): 3,712 tight configurations, **163 distinct Sigma, all 2-decomposable**, 24 not closed under all majorities (`w7_tri.py`, `w7_tri2.py`);
  * a K_3 hub with a branch (chain of members ending in a clean anchor, arbitrary extra unary lists) at each of its three third vertices: **no 3-core** (same files);
  * all tight K_4 gadgets (every role order x 3 hard bits): 521 distinct Sigma, **0 non-2-decomposable**, 56 not closed under all majorities (`w7_k4.py`);
  * triangular prism (I(H) cubic, 6 members, 9 vertices, 3,000 random tight instances with |Sigma| >= 2): 0 non-2-decomposable (`w7_prism.py`);
  * 2,432 random c <= 4 instances with soft triangles (`w7_t4.py`): 2,377 closed under all majorities, 2,432 2-decomposable, 0 cores of size 3; 400/400 for c <= 2 (`w7_t2.py`).
  Every sampled non-2-decomposable Sigma in the literature of this project (the flip gadget G*, c = 3) has soft triangles (8 of them, checked): consistent with Theorem A.
* So the sharp statement is: *a lock core of size >= 3 (and a persistent-flip chain) needs a soft triangle; at c <= 2 that is a triangle of I(H), hence w >= 2; the evidence says triangles of bounded blocks still do not produce one.*
  Proving 2-decomposability for c <= 2 with triangles (Conjecture C2') is OPEN; it would follow from a Helly-type argument for K_3/K_4 blocks glued along a tree (finite computations above are the block case only).

## 7. Corrections and pointers to earlier notes
* wave 6 Sec. 4.7 / audit 8: 'Sigma closed under all-distinct -> M for c <= 2' is false (triangles) but **true without soft triangles** (Theorem A); Lemma 7b can be applied there.
* The abstract 'tunnel-release at w <= 1 is FPT' of Prop. tTR needs only Theorem A; cores of size >= 3, if any exist for c <= 2, live on triangles.
* FAR/N*(F): at w = 1, c = 2, l = 1 these instances are consistent with Theorems A/C: their pointer is a *single* chord or a single trigger; the width needed to share cost among many pointers comes only with unbounded cores.

## 8. Status table

| item | status |
|---|---|
| Lemma 1 (constraint structure: binary relations only without soft triangle) | PROVED |
| Lemma 2 (U-only identity, locality of effect) | PROVED; U-only is a sub-model: VERIFIED (119/120 equal) |
| Theorem A (Sigma closed under every majority, no soft triangle) | PROVED (classification by computer, per-class hand proof, independent graph-level check, 1,568 instances) |
| A1: l = 0 polynomial for soft-triangle-free, any c, any w | PROVED modulo cited Jeavons-Cohen-Cooper; VERIFIED (1,500/1,500) |
| A2: cores <= 2; Conjecture C2 for triangle-free I(H) with c <= 2 | PROVED |
| A3: tunnel-release model FPT, 2^l n^O(1) | PROVED |
| Lemma S (lock = unit propagation, strictness) | PROVED (finite check, 3,888 pairs); hub 3-core absent: VERIFIED (fixpoint) |
| Several forked locks per branch: Helly-2 | OPEN (evidence only) |
| Theorem B: LPPH W[1]-hard (perfect-matching links, 2k paths) | PROVED (own proof), cross-checked by brute force (850 instances, 0 mismatches) |
| Theorem C: tunnel releases cannot realise Theorem B (cores of size N) | PROVED for soft-triangle-free; for other operations HEURISTIC |
| cost sharing semantics (OPT 2/4/3/2), FAR m = 3..6 (7/9/11/13 optima, chord always deleted), chords carry no gain | VERIFIED-BY-COMPUTER |
| RM2 FPT (via packing-endpoint lemma) | OPEN lemma; lemma VERIFIED (11,109 pairs); RM2 is a model, relation to CEaMP HEURISTIC |
| c <= 2 with triangles: Sigma 2-decomposable, no 3-core | VERIFIED-BY-COMPUTER (K_3 blocks, K_4, prism, random); no proof |
| FPT in (l,w) for w <= 1 / for c <= 2 | **OPEN** |
| W[1]-hardness in (l,w) for w <= 1 / for c <= 2 | **OPEN** |

## 9. Reproduction (all in `code/`, python3 + numpy + networkx; gcc for the DP)
`python3 w7_rel.py` (single triple types), `w7_rel2.py` (no uniform majority for the 18 single types), `w7_rel3.py` then `w7_rel4.py` (per-coordinate rules; rel3 takes ~2 min and writes /tmp/w7_compat.pkl), `w7_rel5.py` (96 conjunction relations x 729 rules, ~2 min), `w7_rel6.py` (the 40 relations),
`w7_rel7.py` (independent graph-level check, 576 configurations), `w7_t4.py 11 4000`, `w7_pc.py 1 1500`, `w7_strict.py`, `w7_helly.py single|lists`, `w7_hub2.py`/`w7_hub3.py`/`w7_hub4.py` (abstract calculus),
`w7_tri.py`, `w7_tri2.py`, `w7_k4.py`, `w7_prism.py 1 3000`, `w7_far.py`, `w7_chords.py 6`, `w7_pathcores.py`, `w7_uonly.py 1 120`, `w7_phlp.py 1 300`, `w7_phlp2.py 2 400 2 4 3`, `w7_phlp2.py 3 150 3 3 2`, `w7_rm2.py 1 3000`, `w7_rm2b.py 5 1500`, `w7_rm2c.py 3 2500`.
Build the C DP with `gcc -O2 -shared -fPIC -o w7_dp.so w7_dp.c` (used by `w7_dpy.py`, `w7_chords.py`, `w7_uonly.py`).

## 10. What would settle the question
* **FPT** (w <= 1): prove the packing-endpoint lemma and (L1)-(L4) of Section 5, i.e. that every operation of an optimal solution is either an unpaired local deviation or a pointer whose two ends are at endpoints of conflicts in a fixed maximum packing.
* **Hardness**: exhibit a paired operation with *interior* effect (Theorem B needs it); by Theorem C it is not a tunnel release, so it must be an inserted far pair, a combination with breaks, or a soft-triangle gadget (w >= 2 at c = 2).
* **c <= 2 with triangles**: a Helly argument for K_3/K_4 blocks glued along a tree would give Conjecture C2' and the same FPT for the tunnel-release model.
