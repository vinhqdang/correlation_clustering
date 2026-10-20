# Wave 6: W[1]-hardness attempt for CEaMP in (l + tw(I(H))) with "tunnels"

Task: try to prove W[1]-hardness of Cluster-Editing-above-packing (CEaMP) parameterised by l + w, w = tw(I(H)),
by a reduction from Grid Tiling / Multicolored Clique that uses *tunnels* (hard edges between far vertices; hard
edges are transitive up to soft pairs by HE) to carry colour equalities, with l = O(k^2), w = O(k) (or poly k).
Handle the two obstacles of wave5_locality.md Sec. 6: (i) an inactive anchor must cost nothing, (ii) HUB+ must hold.
No git, paper/ untouched. Code: `code/w6c_*.py` (library `w6c_lib.py`). Outputs of long searches go to `$W6C_OUT` (default `/tmp/w6c_out`).

Labels. PROVED = proof given here (uses only the lemmas of paper/sections/packing.tex, locality.tex, dproute.tex and
analysis_01.md Obs. 1-2). VERIFIED-BY-COMPUTER = exhaustive or sampled search, the script is named, evidence only.
HEURISTIC = argument or sketch without full proof. OPEN = not decided.

Conventions: soft pair = packed pair, hard pair = unpacked pair; Q(v) = soft partners of v; N_h(v) = hard neighbours;
members of H are triples (a,b,c), centre b, edges ab, bc, non-edge ac. For a solution S, ex(S) = |S| - |H| >= 0.
A solution with ex = 0 is *zero-excess* (ZE). States of a member: L = delete ab, R = delete bc, M = insert ac.

## 0. Verdict

1. The reduction was **not completed**. I can neither prove W[1]-hardness nor an FPT algorithm.
2. Neither of the two stated obstacles is the real obstacle. (i) Inactive anchors cost nothing: PROVED (Lemma 2, and the
   instance in Sec. 4.1 verifies cost sharing). (ii) HUB+ is trivially satisfied by K_2-tunnels (hard degree 1); what HUB+/NF(gamma)
   really do is to *forbid* tunnels of large fan-out (Lemma 8), so a tile relation S_ij cannot be written into the tunnels;
   the only way left is a position encoding (one anchor per tile, lex order) and tunnels that realise a *permutation*.
   None of HUB+, SM, T1/T2-locality blocks that.
3. What a reduction of this kind needs is exactly one object, a **persistent-flip chain** (PFG): a tight instance of bounded
   I(H)-width with N clean anchors whose family of "simultaneously lockable" anchor sets has a minimal non-face of size N
   (Lemma 7 shows this is equivalent to a bit that is stored in every separator and flips exactly across the released anchor).
   Conditional statements: bounded minimal non-faces => FPT in the tunnel-release model (Prop 9, PROVED); a PFG family + comparators
   => W[1]-hardness (Sec. 3, HEURISTIC sketch of the Grid Tiling reduction).
4. A **local** flip gadget exists: found by computer, 7 non-isolated vertices, 5 members, tight (OPT = |H|), a clean anchor whose
   lock is an equality between two interface members and whose release flips it (Sec. 4.2, VERIFIED-BY-COMPUTER). Locks always act
   as unit clauses (Lemma 5), so this is non-trivial; no smaller flip gadget exists in the searched space (3 and 4 members: none).
5. **Composition** of flip gadgets into a chain failed in every way I tried (Sec. 4.5). One reason is proved (Lemma 6, occupancy):
   the two vertices a lock bridges are occupied, so they cannot be shared with the next gadget. The search for gadgets that have
   separate clean sockets found nothing (Sec. 4.6). This is the missing piece; it is a finite gadget question, not a structural theorem.
6. Evidence that minimal non-faces of size >= 3 do not exist in *small* tight instances: no violation in 723 + 35,815 + 12,351
   connected tight graphs (n <= 9) and in ILP tests up to n = 17 (Sec. 4.4). Not a proof; Li et al.'s NP-hardness at l = 0 (c <= 23)
   shows that non-bijunctive local gadgets exist, so large instances may differ.
7. Structural reason (Lemma 7b, PROVED + Sec. 4.7, VERIFIED): if the zero-excess set is closed under a majority operation, lock cores have size <= 2 and the tunnel-release model is FPT. Empirically this holds for **all** instances with c <= 2
   (one explicit rule, "M-majority"), and fails for c >= 3 (e.g. the flip gadget G* itself is not 2-decomposable). So a hardness proof needs c >= 3 and non-bijunctive gadgets; for c <= 2 the closure (Conjecture C2) is open but would settle the model.

8. Which lemma blocks (precise answer): none of HUB+, SM, T1/T2 locality blocks the tunnel/permutation design (Sec. 5). What limits it is the rigidity of zero-excess solutions (Lemmas 1, 5, 6: locks are unit clauses, release closes y, bridged vertices are occupied)
   together with the fact that long lock cores need a non-bijunctive zero-excess set (Lemma 7b), which (empirically) needs c >= 3.

## 1. Definitions

*Zero-excess structure* (used throughout). By the decomposition lemma every member P_i contains >= 1 edited pair, pair sets are
disjoint, hence |S| = sum_i c_i + |U| with c_i >= 1 edits inside P_i and U = edited hard pairs. ex(S) = 0 iff c_i = 1 for all i and U = empty.

*Tunnel.* A hard edge xx' with N_h(x) = {x'}, N_h(x') = {x} and Q(x) cap Q(x') = empty. (A *private stub* is a tunnel in which x' lies in
no member and N(x') = {x}.) x, x' are the *ends*. The tunnel is *respected* by S if cluster(x) = cluster(x') = {x,x'}. An end is *clean*
if the tunnel is its only hard edge (then all its other edges are soft).

*Lock / release.* Putting a tunnel on a vertex x is a *lock* on the members through x (Lemma 2). *Releasing* it means cutting xx'
(one hard edit). For a set A of anchors (clean ends of tunnels, each with a private stub), a set S of anchors is a *face* if some
ZE solution of G (without stubs) makes every x in S a singleton cluster, equivalently OPT(G + stubs on S) = |H|. A *minimal non-face*
is a non-face all of whose proper subsets are faces. (Equivalence with the statement about stubs: a stub on x costs 0 iff x is a singleton, Lemma 2/3.)

## 2. Proved lemmas

**Lemma 1 (zero-excess structure) PROVED.** If ex(S) = 0 then (a) every hard edge is inside a cluster and every hard non-edge across
clusters; (b) each member is in state L, R or M; (c) every cluster is a union of hard components all of whose cross pairs are soft.
Proof. ex = 0 forces U = empty, giving (a), and c_i = 1. A member (a,b,c) has final togetherness t(ab),t(bc),t(ac) in {000,100,010,001,111}
(cluster graph); the costs are 2,1,1,3,1 and the cost-1 patterns are (1,0,0)=R, (0,1,0)=L, (1,1,1)=M. (c): two vertices in one cluster
are adjacent in G+S, a hard pair inside a cluster is an unedited edge. QED.

**Lemma 2 (clean tunnel = lock) PROVED.** Let xx' be a tunnel and S a ZE solution (more generally a solution with xx' and all pairs at x,x'
unedited). Then cluster(x) = cluster(x') = {x,x'}. Consequently every member through x has x as an end, with x's edge deleted
(if x were the centre of a member both its edges would be deleted: c_i = 2, impossible at ex = 0), i.e. its state is forced; the other two vertices
y,z of a member (x,y,z) are together.
Proof. xx' is a hard edge, unedited, so x ~ x'. Let t be another vertex of the cluster. If t is not in Q(x), then xt is hard; since N_h(x) = {x'},
xt is a hard non-edge inside a cluster, hence edited. So t in Q(x); likewise t in Q(x'); but Q(x) cap Q(x') is empty. QED.

**Lemma 3 (cost of non-respected tunnels) PROVED.** Let T be a set of vertex-disjoint tunnels and S a solution; let rho(S) = number of tunnels not respected.
Then ex(S) >= rho(S)/2, and ex(S) >= number of tunnels whose edge is deleted.
Proof. A non-respected tunnel xx' either is cut (xx' in U), or x ~ x' with a third vertex t in the cluster; then by the argument of Lemma 2 one of the
hard pairs xt, x't is a hard non-edge inside a cluster, i.e. in U. Charge the tunnel to that pair of U. A pair is charged by at most two tunnels (one per end).
ex(S) >= |U|. QED. (The factor 2 is real only for inserted pairs between ends of two different tunnels; cuts are charged exactly once.)

**Lemma 4 (release closes the cluster) PROVED.** Let x be a clean end of a tunnel and S a solution in which the only edited hard pair at x is xx' (deleted).
Then cluster(x) is contained in {x} union Q(x). So if P = (x,y,z) is the only member through x, P is in a cost-1 state and x is together with y, the cluster of x is {x,y} or {x,y,z}.
Proof. Let t be in cluster(x), t != x. If xt is hard it is a non-edge (the only hard edge xx' is deleted, so x' is not in the cluster), hence edited: contradiction. QED.
So releasing frees z from y, but y becomes a closed pair {x,y} (or {x,y,z}): release never adds freedom at y.

**Lemma 5 (forcing at a lock; EQ) PROVED.** In a ZE solution with y ~ z (a locked member (x,y,z)): for every q not in {y,z},
t(q,y) = t(q,z) (t = togetherness). Hence if (q,y) is soft and (q,z) is hard, the soft pair (q,y) is *forced* to the value [qz in E]; if both
are soft, the owners of (q,y) and (q,z) are coupled by equality, and the three members through y,z,q form a triangle of I(H).
Proof: transitivity. QED.
Consequence: around a lock all members through y or z are forced unless a soft triangle {y,z,q} exists. This is why the exhaustive search of
7-vertex couplers (Sec. 4.2) shows that in every pattern with a zero-excess state with x isolated exactly one pair of neighbour states is allowed.

**Lemma 6 (occupancy) PROVED.** In a ZE solution, if (u,w) and (w,v) are both together then u ~ v, so (u,v) must be an edge or a soft pair.
In particular a vertex w cannot be an end of two locked bridges (y~z kept and z~z' kept) unless (y,z') is soft or an edge.
Proof: transitivity and Lemma 1(a). QED. This explains the failure of the naive ring of flip gadgets (Sec. 4.5).

**Lemma 7 (splice, memory) PROVED.** Let Sigma = set of ZE solutions as vectors in {L,R,M}^H, (T,X) a tree decomposition of I(H), e an edge of T,
M_1, M_2 the members that occur only on one side of e and B the members in both bags of e. Then (a) if sigma, sigma' in Sigma agree on B, the vector equal to
sigma on M_1 and B and to sigma' on M_2 is in Sigma. (b) Let A be a set of anchors, S = S_1 u S_2 a minimal non-face where all members through
anchors of S_i lie in M_i. For b in S_1 let sigma^b in Sigma satisfy all locks of S minus {b}; then
{sigma^b|_B : b in S_1} and {sigma^a|_B : a in S_2} are disjoint.
Proof. (a) By Obs. 1-2 of analysis_01.md the zero-excess condition is a conjunction of unary constraints (hard pairs) and constraints per vertex triple whose scope is a
clique of I(H), hence inside one bag, hence entirely in M_1 u B or in M_2 u B. (b) If the restrictions to B agree, splice sigma^a (a in S_2) on M_1 u B and sigma^b (b in S_1)
on M_2: all locks of S_1 hold on M_1 (a is not in S_1) and all locks of S_2 hold on M_2 (b not in S_2): S would be a face. QED.
Reading: a minimal non-face of size N forces the bag state at each of the N-1 separators between consecutive anchors to remember on which side the "unlocked" anchor lies:
a persistent bit that flips exactly across the released anchor (PFG).

**Lemma 7b (bijunctive => cores of size <= 2) PROVED.** If the set Sigma of ZE solutions (vectors in {L,R,M}^H) is closed under a majority operation m (the same m on every coordinate,
m(a,a,b) = m(a,b,a) = m(b,a,a) = a), then every minimal non-face has at most two anchors.
Proof. Anchors force the states of the members through them (Lemma 2); a set S of anchors is a face iff the resulting partial vector t on the coordinate set T extends to Sigma, i.e. t is in proj_T(Sigma). proj_T(Sigma)
is closed under m, so by the Baker-Pixley theorem (majority term => relations are determined by their binary projections) t is in proj_T(Sigma) iff every binary subtuple of t is in the corresponding binary projection of Sigma.
A binary subtuple involves at most two anchors; it extends if the pair of anchors (or the single anchor) is a face. So if all pairs are faces, S is a face. QED.
Hence a persistent-flip chain needs a Sigma that is *not* closed under any majority, i.e. non-bijunctive local structure (Sec. 4.7 shows where that starts).

**Lemma 8 (tunnel fan-out is bounded; bounded-degree Grid Tiling is FPT) PROVED.** (a) In a YES instance (OPT <= |H| + l) with pair-maximal H, by Lemma EX(a) (dproute.tex)
every endpoint of an edited hard pair has hard degree <= D + u <= D + l, and by Lemma NF(gamma) every hard component containing a vertex of hard degree > B = 2D + 3l is a whole,
non-exceptional cluster (frozen: no pair of it is edited). So a tunnel end that is ever released has hard degree <= D + l = poly(l, w). By HE the hard neighbours of a vertex are a
clique up to soft pairs, so a vertex joined by tunnels to many far anchors forces all of them into one big hard clique, which is frozen. Hence the relations realised by releasable
tunnels have bounded degree (partial matchings, or cliques of size <= D + l + 1).
(b) Grid Tiling (find a_1..a_k, b_1..b_k with (a_i,b_j) in S_ij) restricted to instances where every a and every b has at most d partners in every S_ij is solvable in time
n * d^{2k} * poly: guess a_1 (n ways), then each b_j has <= d choices, then each a_i has <= d choices given b_1. QED.
So a Grid Tiling reduction cannot write S_ij into tunnels of unbounded fan-out (this answers the wave-5 worry about "n anchors sharing hard neighbours"); it must encode a tile as one *position* on a
track and use tunnels only as a *matching* between two orderings of the same tile list.

**Prop 9 (tunnel-release model: bounded cores => FPT) PROVED.** Model TR: G has a packing H and a set T of clean tunnels (each tunnel is a hard K_2 xx' whose ends are
anchors of members, or an anchor with a private stub; tunnels pairwise vertex-disjoint), and we only look at solutions whose hard edits are cuts of tunnels. For R subset of T let f(R) = 1 iff
G minus the tunnel edges of R has a ZE solution; call S subset of T *lockable* if f(T minus S) = 1 (some ZE solution respects all tunnels of S; by Lemma 2 this makes all ends of S singletons).
Then (i) f is monotone, (ii) f(R) = 1 iff R meets every minimal non-lockable set ("core"), (iii) min |R| with f(R) = 1 is the minimum hitting set of the cores. If every core has at most d elements the
minimum can be computed in time d^l * 3^{w+1} * n^{O(1)} (bounded search tree: in the current instance find a core by deleting tunnels one at a time from T minus R, calling the DP of Theorem l0 after each deletion; branch on cutting one of its <= d tunnels).
Proof. (i) Given a ZE solution pi of G minus R' and R' subset R, split every cluster {x,x'} of a tunnel of R minus R' into two singletons; both ends were singletons apart from each other in pi (Lemma 2),
so no other pair changes. (ii),(iii) are the definitions; the algorithm uses that lockability is decided by the l = 0 DP over a width-w decomposition of I(H) (which does not change when tunnels are cut). QED.
Contrapositive: any W[1]-hardness proof inside TR must produce minimal non-lockable sets of unbounded size, i.e. (Lemma 7) a PFG. (Single-anchor faces and tunnel cores are the same notion for stubs; for paired tunnels one
tunnel is one element and releases two anchors at once: this is where the sharing of Sec. 4.1 lives.)

## 3. The reduction blueprint (what would be proved if the missing gadget existed)

**3.1 Abstract reduction (HEURISTIC sketch, all steps standard).** Source: Grid Tiling, cells (i,j), S_ij subset of [n]^2, need (a_i, b_j) in S_ij; W[1]-hard.
* Cell (i,j): list its tiles e = (a,b) in S_ij. Track A_ij: positions in lexicographic order (a,b), padded so that block a occupies positions ((a-1)n, an]; one anchor
  per tile. Track B_ij: the same tiles in order (b,a), padded likewise. For each tile a tunnel (hard K_2, a *perfect matching*, fan-out 1) between its anchor on A_ij and its anchor on B_ij.
* Track gadget = PFG chain: all N anchors locked => infeasible at ex = 0 (boundary terminals force the bit to flip); releasing any single anchor makes it feasible and the bit
  b_t = [released position <= t] is available at every position t.
* Comparators: at the block boundaries t = alpha*n, equality of the bits of A_ij and A_{i,j+1} (resp. of B_ij and B_{i+1,j}) by an equality wire between members of the two tracks
  (one shared vertex; the 2-member equality wires exist, Sec. 4.3). Then [a_ij <= alpha] = [a_{i,j+1} <= alpha] for all alpha, i.e. a_ij = a_{i,j+1}.
* Budget l = k^2. Parameters: tracks are parallel chains swept in lock step: pathwidth O(k^2 * g) where g is the width of one chain period; c bounded by the gadget.
* Correctness: a tiling gives one matched tunnel release per cell (cost 1 per cell, shared by the two tracks: the cost-sharing example of Sec. 4.1). Conversely each track needs >= 1
  release (the all-locked chain is infeasible; boundary relaxations are unshared and cost 1 per track too), a release is a tunnel cut or an unshared edit, and a cut is shared
  only between A_ij and B_ij of the same cell; budget k^2 leaves exactly one matched shared tunnel per cell, hence a tile per cell, comparators give the equalities.
  The tile identity is carried by the *matching*, never by a hard neighbourhood: fan-out 1, so HUB+ and NF(gamma) hold trivially (alpha(G[N_h(x)]) = 1).

**3.2 Where the reduction needs the gadget.** Everything above holds if a PFG family exists:
GH(N): graphs G_N, packings H_N (valid), clean anchors x_1..x_N, terminals, such that (G1) c <= c_0 and pathwidth(I(H_N)) <= w_0 independent of N; (G2) G_N without stubs is tight (OPT = |H_N|);
(G3) a stub on every x_j gives OPT = |H_N| + 1 while a stub on all but one gives OPT = |H_N| (minimal non-face = all N anchors); (G4) a tap at every position that reads the bit; (G5) equality wires between taps of two chains.
Lemma 7 shows G3 forces a persistent bit; Prop 9 shows G3 is exactly the dichotomy for the TR model. G2/G3 are what the searches below test; G4/G5 were not reached.
Status of the statement "GH => CEaMP is W[1]-hard in l + w": HEURISTIC (sketch 3.1 is routine, but I did not write out the correctness proof, since GH is not established).

## 4. Computer evidence

All scripts are in `code/`. "n" counts vertices; ILP = HiGHS via scipy (`w6c_lib.ilp_opt`), brute force = all set partitions.

### 4.1 Cost sharing and inactive anchors (VERIFIED-BY-COMPUTER, `w6c_share.py`)
P = (x,y,z), P' = (x',y',z') vertex-disjoint, z has a hard neighbour t (t-y non-edge), z' likewise, so each member's only ZE state is R (x joins y). Results (|H| = 2):
no tunnel OPT = 2; private stubs at x and x' OPT = 4; **one tunnel x-x' OPT = 3** (one cut serves both members: cost sharing); same tunnel when the members' default state L is
allowed (no forcing) OPT = 2 (**an inactive anchor costs nothing**). All instances are valid, pair-maximal packings.
(This is the c = 1 situation of Remark expscope in dproute.tex: two vertex-disjoint P3s joined by a hard edge that is deleted.)

### 4.2 A local flip gadget exists (VERIFIED-BY-COMPUTER, `w6c_flip4.py`, `w6c_flipgadget.py`, `w6c_minitrack.py`)
Search (`w6c_flip4.py n m T seed`): random configurations of m members on n vertices with an anchor member P = (x,y,z) (x clean private end), *all* 2^h patterns of the h
free hard pairs; relation of two interface members (I1, I2) over ZE solutions with x singleton (locked) and x non-singleton (released); FLIP = locked relation is a bijection between
two 2-sets D1 -> D2, released relation contains a pair (d,e), d in D1, e in D2, e != locked image.
* n = 7, m = 4: 27,130 configurations, 234,567 patterns with a ZE solution, 29,334 bijective locked relations, **0 flips**.
* n = 8, m = 4: 433 configurations, 12,055 patterns, 904 bijective, **0 flips**.
* n = 8, m = 5: 4,364 configurations, 66,874 patterns, 16,086 bijective, **128 flips** (further seeds: 1,860 flip records in total, mostly isomorphic variants).
Gadget G* (first record): vertices 0,1,3,4,5,6,7 (vertex 2 unused). K4 on {0,1,4,5}, hard edges 04 and 15 (the other four K4 edges are soft), a = 3 adjacent to 0,1,6,7,
b = 6 adjacent to 0,5,3, x = 7 the anchor. Members m0 = (6,3,7), m1 = (4,1,3), m2 = (5,0,3), m3 = (4,5,6), m4 = (6,0,1). OPT = |H| = 5, pair-maximal. The five ZE solutions:
x singleton: {0,1,3,4,5,6} (states of (m3,m1) = (M,M)) and {0,1,4,5},{3,6} (R,R); x not singleton: {0,1,4,5,6},{3,7} (M,R); {0,1,4,5},{3,6,7} (R,R); {0,1,4,5},{3,7},{6} (R,R).
So locked: (I1,I2) = EQ on {M,R}; released: I2 = R and I1 free: the flip M -> R appears only when x is released (cut), as predicted by Lemma 4 (y closed, z freed).
Mini-track (`w6c_minitrack.py`, n = 9): add t hard-adjacent to {0,1,4,5,6} (forces b into the K4 cluster). Without a stub OPT = |H| = 5 with the unique optimum {0,1,4,5,6,8},{3,7};
with a private stub on x OPT = 6. So "forced flip needs a release (cost 1)" is realised, by a tight pair-maximal instance, with c = 3 (vertex 3 lies in m0,m1,m2).

### 4.3 Two-member relations, couplers, triangles (VERIFIED-BY-COMPUTER, `w6c_rel.py`, `w6c_coupler.py`, `w6c_tri.py`, `w6c_tri3.py`, `w6c_hub3.py`)
* 5-vertex gadgets (2 members sharing one vertex): 31 distinct zero-excess relations; the 2-element ones include equality {LL,RR}, {LR,RL} and shifted bijections (wires exist).
* 7-vertex couplers (member B, lock member A = (u,u',x) and B', all 12 hard pairs, 12 shapes: 49,152 graphs): only 104 graphs have a ZE state with x isolated, and in all 104
  exactly one pair of states of (B,B') is allowed (so no flip can occur with a single shared vertex per side). This is Lemma 5 (forcing).
* Soft-triangle lock (Q,P,E pairwise sharing a vertex, 6 vertices, 144 configurations x 64 patterns): locked relations of size 2 (equality between s(q,y) and s(q,z)) occur 144 times;
  released relations are single pairs; **no flip inside the domain** (16 patterns have a new released pair, but with an input state outside the locked domain).
* Hub of three anchored members around a central member (9 vertices, 48 configurations x 512 patterns): 0 three-element non-faces.

### 4.4 Small cores (VERIFIED-BY-COMPUTER, `w6c_nonface_singletons.py`, `w6c_nonface_components.py`, `w6c_core.py`, `w6c_core2.py`, `w6c_necklace.py`)
Question: is there a tight instance with a minimal non-face of size >= 3 among lockable anchors (pairwise faces, jointly not)? For singleton anchors and also for hard components
locked as whole clusters (twin stubs):
* all 723 connected tight graphs with 5 <= n <= 7 (atlas): none; random connected tight graphs: n = 8: 35,815 (of 62,474 sampled), n = 9: 12,351 (of 27,684): none;
  hard-component locks: n = 8: 27,624 tight instances (of 44,670), n = 9: 5,059 (of 11,015): the largest minimal non-face per (instance, packing) has size 0, 1 or 2, never 3;
* ILP tests (maximal cliques of the pairwise-face graph tested for joint face-ness): planted-cluster graphs n = 11 (1,258 tight), n = 13 (833 tight), random n = 12 (251 tight), necklaces
  of N <= 8 triangles, paths/cycles with chords (n up to 17): no violation. The triangle necklace is tight for every N (OPT = max packing = floor(3N/2) - 1, N = 2..8).
* Also none for the structured searches of Sec. 4.5 (below).
Reading: lock cores of size 3 do not occur in these small instances, although non-bijunctive zero-excess sets do occur at c = 3 (Sec. 4.7, 11% of the c = 3 instances with n = 8): clean-anchor locks only force the states L/R of members through clean ends, so a
non-2-decomposable Sigma is necessary (Lemma 7b) but not sufficient for a long lock core.

### 4.5 Composition attempts (VERIFIED-BY-COMPUTER, all negative)
* Ring of N K4-joiner gadgets (`w6c_ringK4.py`): G* units, joiner w_j shared (a of unit j = b of unit j+1). Valid, pair-maximal, tight for N = 3..7 (tw(I(H)) <= 8 by min-fill, c = 6), but for N >= 4 no
  single anchor can be a singleton at zero excess (all ZE solutions attach every x_j): minimal non-faces have size 1. Reason: Lemma 6 (the shared vertex lies in two bridges).
* Other ways to chain G* (`w6c_ringsearch.py`, 144 options, N = 3 and, for the promising ones, N = 5): 108 "wire member through a new vertex" options: 85 not tight (excess +3 at N = 3: the wire vertex cannot join an occupied
  vertex), 23 tight but with no lockable anchor at all; 36 "identify one vertex of gadget j+1 with one of gadget j" options: 20 not tight, 8 tight with *all* anchors simultaneously lockable (no core; 4 of them also at N = 5),
  8 tight with no lockable anchor. `w6c_piecechain.py` (units sharing a hard K2 piece, 16 variants, N = 3, 4): never tight.
* Composition of the 1,860 flip records by identifying the right interface member with the left one (`w6c_compose2.py`) or a clean private socket vertex with another (`w6c_compose3.py`, 1,160 candidate socket
  pairs on the first 850 records; a random sample of 60 records re-checked): every composite is invalid or has ex >= 1 without stubs (no zero-excess solution).
* Alternating chain O,E,O,E,.. (`w6c_altchain.py`, `w6c_altpath.py`): O = G*, E = K4 with two slots and no bridge. Tight and compatible, but the relation between the two end slots with all anchors locked is
  the full relation {0,1}^2: E does not transmit the bit.
* Odd ring of K2 pieces with a satellite anchor per piece (`w6c_piecering.py`; idea: odd ring needs a perfect matching of pieces, releasing a satellite lets one piece leave): valid, N = 5, 7 have OPT = |H| + 2 even
  without stubs and releasing a satellite does not repair it (the clause of the departing piece still needs a ring neighbour: a piece cannot be "absorbed" by a satellite without an OR3 clause); N = 6 is tight with all anchors lockable.
* Graph-level random ring search (`w6c_graphring.py`, K4 units with slot and anchor adjacency patterns, maximum packing computed by ILP, N = 3): 1,025 tight graphs, 0 with a size-3 core.
* Rotation-symmetric ring of 4-vertex periods, 3 periods (`w6c_ring.py`): 2,087 candidates, 99 tight, 0 cores.

### 4.6 Socket gadgets for an OR-wire chain (VERIFIED-BY-COMPUTER, negative)
Design: gadgets with clean sockets l, r (ends of members, no hard edge) and flags f = "singleton"; locked => not (f_l and f_r), released allows both; joined by OR-wires (member (r_j, t_j, l_{j+1})). Such a chain gives a
Horn implication chain through all locks, i.e. a PFG (the derivation of infeasibility uses every lock as a premise). Search `w6c_sock.py` (n = 7, m = 4: 23,146 configurations, 200,325 patterns; n = 8, m = 4: 417 / 11,952;
n = 8, m = 5: 4,949 / 75,423) and the template-seeded `w6c_sock2.py` (n = 9, 8,114 configurations): **no gadget**. `w6c_flip5.py` (flip gadgets whose interface members each have a clean *private* socket vertex, followed by the N = 2 composite test): n = 8, m = 5: 71,791 configurations passing the socket condition, **0 flip gadgets** (while 128 flips
exist without the socket condition in 4,364 configurations); n = 9, m = 5: 0 (most configurations skipped for having more than 14 free hard pairs, so weak evidence); n = 9, m = 6: 323 configurations, 0.

### 4.7 Majority closure and non-bijunctive gadgets (VERIFIED-BY-COMPUTER, `w6c_majority.py`, `w6c_majority_ex.py`, `w6c_triplerel.py`)
For random connected tight instances (n = 8, 9; all optimal partitions = ZE solutions; instances with 3 <= |Sigma| <= 40) I tested closure of Sigma under each of the 729 uniform majority rules (m on three distinct states arbitrary),
grouped by c = max number of members through a vertex:
* c <= 2: **all** instances closed (n = 8: 761/761; n = 9: 225/225; a separate sample of 301 instances with n = 7, 8: 301/301). The *same single rule* works for all of them: the all-distinct triple (L,R,M in any order) maps to M
  ("M-majority"); equivalently t*(ab) and t*(bc) are the majority togetherness and t*(ac) = t*(ab) and t*(bc).
* c = 3: not closed under any rule in 265 of 2,384 instances (n = 8) and 20 of 816 (n = 9); c = 4: 75 of 443 (n = 9). Example (n = 9, c = 3, 6 members, tight, valid, pair-maximal): Sigma has 5 states and is **not 2-decomposable**
  (the vector MRLRMR is pairwise consistent but not a solution).
* The flip gadget G* of Sec. 4.2 is itself not 2-decomposable: Sigma = {LRRRL, LRRMM, RRRRL, RMMMM, MRRRL} (members m0..m4) and RRRMM (m0 = R is the lock state, m1 = R, m3 = M) is pairwise consistent but not a solution. This is the
  non-bijunctive ternary relation (lock, input, output) that Lemma 7b says a persistent flip needs.
* State-level 3-cores (`w6c_statecores.py`, `w6c_statecores2.py`): a pairwise-consistent but non-extendable assignment of L/R values to three members (the only states a clean-end lock can force) occurs in 118 of 10,658 (n = 8) and 8 of 1,332 (n = 9) tight
  instance-packings (127/14,161 and 20/1,984 in the second sample); in 79 resp. 11 of them there are >= 3 anchor-capable vertices, but **none** has a lock core of size >= 3 for singleton locks. So non-bijunctive structure with L/R-valued cores exists at n = 8, c = 3;
  what is missing is realising the forced states by isolated clean ends. First attempt (`w6c_leaflock.py`: attach to each core member a leaf member with a clean anchor sharing one vertex, 48 options per member, forced state checked by ILP): of 265 core instances (n = 8, 9),
  18 have a forcing option for all three core members; for those, 24-32 option combinations each were tested (ILP, face structure of the three anchors): **no lock 3-core** (VERIFIED-BY-COMPUTER, negative, small option space).
* Closure of Sigma is *not* a local fact: 57 of the 117 role patterns of a single triple constraint (`w6c_triplerel.py`) are individually not closed under M-majority, yet no c <= 2 instance violates it.
Reading: (a) PFG cannot exist for c <= 2 if the c <= 2 closure is a theorem (this would also give the polynomial l = 0 algorithm for c = 2 of Li et al. as a 2-SAT-type argument and, with Lemma 7b and Prop 9, FPT for c = 2 in the tunnel-release model);
(b) c >= 3 is where non-bijunctive gadgets (and so flips) begin. The closure for c <= 2 is **not proved** here (OPEN, Conjecture C2: for c <= 2 the ZE set is closed under M-majority).

## 5. Where it blocks, lemma by lemma

* **HUB+ (alpha(N_h) <= 2l + c + 1).** Not violated: K_2-tunnels have hard degree 1. It only excludes tunnels with large fan-out (Lemma 8). Consequence: tile relations cannot be put into tunnels; positions + a permutation matching
  are needed; that in turn needs comparisons of *positions* across tracks, which needs a persistent bit.
* **SM / NF(beta,gamma).** Not violated: clusters in the design are of size <= 2D. SM says that pieces that merge are small, hence frozen whole clusters impose only unary constraints: this is what makes an unreleased tunnel a lock (Lemma 2)
  and a released one a closure (Lemma 4).
* **T1/T2 locality.** Not used: the design only uses *deleted* hard edges (R2 of locality.tex), for which locality gives nothing and Cor. cand gives only the bounded candidate set for the partner,
  consistent with a matching.
* **HE / pair-maximality (MAX).** Used only to say that a clean tunnel end has all its other pairs soft; they give the rigidity of Lemma 5 (hard pairs fix the togetherness of all neighbours) which is what limits couplers.
* **What actually blocks is rigidity of zero-excess solutions (Lemmas 1, 5, 6)** combined with the polarity of tunnels: a tunnel is free exactly when its ends are singletons, so inactive anchors are *isolated*, locks are unit clauses ("y ~ z"),
  and a release only closes y and frees z. A long core therefore needs a persistent bit that crosses each lock; the bit can only pass through soft triangles {y,z,q} (Lemma 5), and gadgets built that way occupy their two bridged vertices
  (Lemma 6), so they cannot share them with the next gadget. The cheap alternatives are exactly the failed compositions of Sec. 4.5.
* **Not excluded:** a PFG with separate clean sockets and more members/vertices (n >= 10), or a PFG obtained from the l = 0 hardness gadgets of Li et al. (their clause gadgets realise non-bijunctive relations with c <= 23, and
  I(H) of their instance is up to constants the incidence graph of the formula, A2 of hardness.md). The formula Phi_N = (r_0 = 0) and (r_N = 1) and XOR(r_{j-1}, m_j, r_j) over a path of N inputs m_j has bounded pathwidth; if each input m_j could be
  clamped to 0 by a clean anchor (a member whose lock forces the clause variable), the all-locked formula would be unsatisfiable and releasing any anchor satisfiable. Two things are unchecked: pathwidth preservation of their reduction,
  and the existence of a clamp gadget with a clean end (HEURISTIC/OPEN).

## 6. Status table

| item | status |
|---|---|
| Lemma 1 zero-excess structure | PROVED |
| Lemma 2 clean tunnel = lock; inactive anchors cost nothing | PROVED (+ VERIFIED, Sec. 4.1) |
| Lemma 3 non-respected tunnels cost >= 1/2 each, cuts exactly 1 | PROVED |
| Lemma 4 release closes the cluster, frees one side | PROVED |
| Lemma 5 forcing / EQ at a lock | PROVED |
| Lemma 6 occupancy | PROVED |
| Lemma 7 splice, minimal non-face => persistent bit | PROVED |
| Lemma 7b majority-closed Sigma => lock cores <= 2 (Baker-Pixley) | PROVED |
| Sigma closed under M-majority for all tight instances with c <= 2 (n <= 9, 1,287 instances); fails for c >= 3 (G* not 2-decomposable) | VERIFIED-BY-COMPUTER (sampled) |
| Conjecture C2: c <= 2 => Sigma closed under M-majority (would give FPT for c <= 2 in the TR model) | OPEN |
| Lemma 8 tunnel fan-out bounded (NF gamma); bounded-degree Grid Tiling FPT | PROVED |
| Prop 9 TR model: bounded cores => FPT (d-Hitting Set); hardness needs unbounded cores | PROVED |
| cost sharing of a tunnel between two far members (n = 8) | VERIFIED-BY-COMPUTER |
| local flip gadget G* (7 vertices, 5 members, tight) and mini-track | VERIFIED-BY-COMPUTER |
| 3- and 4-member flip gadgets do not exist (n <= 7/8, random configs, all patterns) | VERIFIED-BY-COMPUTER (sampled configurations) |
| no minimal non-face of size >= 3 in small tight graphs (n <= 9 exhaustive/sampled, ILP up to 17) | VERIFIED-BY-COMPUTER (evidence only) |
| composition of flip gadgets into a chain (all tried compositions) | VERIFIED-BY-COMPUTER, negative |
| socket gadgets for an OR-wire chain (n <= 9) | VERIFIED-BY-COMPUTER, negative (sampled) |
| Grid Tiling reduction given GH (Sec. 3) | HEURISTIC (sketch) |
| GH (PFG chain with comparators) exists | OPEN |
| W[1]-hardness of CEaMP in l + w | OPEN |
| FPT in l + w | OPEN |

## 7. Suggestions for the next wave

1. PFG with separate sockets: extend `w6c_flip5.py` to n = 10 (anchor + two private sockets + 7 others) with a smarter generator (templates "lock P, soft triangle {y,z,q}, bit carriers", Lemma 5), or ILP/SAT synthesis of the
   hard-pair pattern instead of enumeration.
2. Use Li et al. gadgets: re-derive their variable and clause gadgets, check (a) pathwidth preservation, (b) a clamp member with a clean end; then the XOR-parity chain Phi_N above gives a PFG with explicit constants (c <= 23).
3. Prove Conjecture C2 (c <= 2: the ZE set is closed under M-majority; all 1,287 sampled instances are). With Lemma 7b and Prop 9 this gives FPT in the tunnel-release model for c <= 2, and shows that every hardness
   proof needs c >= 3. Then study c = 3 where non-bijunctive Sigma first appears (n = 9, 6 members, Sec. 4.7): search for gadgets whose non-2-decomposable part involves *lock states* (clean ends force only L/R of the members through them),
   e.g. by extending `w6c_majority_ex.py` to report minimal pairwise-consistent non-extendable partial assignments consisting of forced L/R states. Inserted pairs (T2) are a second kind of release (an inserted release extends the
   tunnel piece and does not share cost); the TR model ignores them and an FPT proof must treat them (T1 are bag-local by Cor. softfree).
4. Whatever the outcome, the sharing example (Sec. 4.1) and Prop 9 belong in the paper as the precise statement of what remains open: FPT in l + w holds in the tunnel-release model iff lock-conflicts have bounded size.
