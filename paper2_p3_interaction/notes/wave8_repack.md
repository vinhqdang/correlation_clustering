# Wave 8 (repacking): the parameter is a property of the packing, so change the packing

Task: decide whether CEaMP is FPT in (l, w = tw(I(H))) by exploiting the freedom to replace H by another packing H*
(exchange / augmentation / iterative compression), starting from the obstruction "far two-sided edits" (Prop FAR).
No git, `paper/` untouched. All code is `code/w8a_*.py` (library `w8a_lib.py`; scipy/HiGHS ILP, networkx, nauty-geng for exhaustive runs).

Labels (every claim carries one): **PROVED** (proof here), **VERIFIED-BY-COMPUTER** (finite check, script named, evidence only),
**HEURISTIC** (argument without full proof), **OPEN**. "CITED" = taken from the paper / its references (LPS23 NP-hardness at l=0).

Conventions: (G,H,l), H a modification-disjoint packing of induced P3s, pair-maximal unless said otherwise; soft pair = packed pair,
hard pair = other; Q(v) soft partners; W the soft graph; nu(G) = maximum size of a packing; OPT(G); for a solution S, U = S minus soft pairs
(hard edits), |U| <= |S|-|H|. L_B, Z = 3D+3l+1, Ball(x) = ball of radius 2 around x in G[L_B] (|Ball| <= Z^2) as in Section "optimal" of the paper.
For a tree decomposition (T,beta) of I(H): X_t = union of V(P) over P in beta_t; T_v = {t : v in X_t} (a subtree, because the members through v
form a clique of I(H), so their subtrees pairwise intersect); a,b are *co-bagged* if T_a and T_b intersect.
A hard pair {a,b} with a,b in V(H) is *A-local* if members P through a and P' through b satisfy P cap P' != empty (then P != P', and a,b are co-bagged in every
decomposition); equivalently a and b have a common neighbour in the soft graph W. A vertex r outside V(H) is a *local hub* (w.r.t. a decomposition) if the vertices of V(H) in Ball(r) are pairwise co-bagged.

## 0. Verdict

1. **No algorithm and no hardness proof [OPEN].** FPT in (l,w) stays OPEN. What the wave establishes is a precise map of where the obstruction lives once H may be changed.
2. **[PROVED, conditional on P != NP / ETH] Impossibility (sharp form of the requested target).** No polynomial procedure can output, for every YES instance, a packing H* with
   l* <= s(l) and tw(I(H*)) <= f(l) (functions of l alone) unless P = NP (Thm N); under ETH not even tw(I(H*)) = |H|^(1-eps) for l = 0.
   The only non-vacuous version is w* <= w + g(l) with a *locality* property (C1/C2 below), which Thm N does not exclude.
3. **[PROVED + VERIFIED-BY-COMPUTER] Both known obstruction families are artifacts of the choice of H.**
   * FAR_m (Prop FAR): nu = m, OPT = m+1, so excess 1 is *intrinsic* (no packing removes it), but one size-preserving exchange gives H* with
     I(H*) = cycle (tw 2) and a unique hard edit (2m-1,2m) with 2m outside V(H*) (Prop 3.1). The far chord becomes a *soft* pair.
   * Sp(F) = N*(F) (spine family, w=1, l=1, tw(Gamma') unbounded): one exchange X_0 -> Y plus one augmentation Z gives a packing of size |H|+1 = OPT,
     i.e. excess 0 and tw <= 2 (Prop 3.2). It also gives a two-line proof of OPT >= |H|+1 (replaces the two-page irreducibility argument).
     Moreover the *original* H already lies in the locality class of Thm PL (the trigger is a local hub).
4. **[PROVED; semantics VERIFIED-BY-COMPUTER] New positive result (promise-FPT), Thm PL.** If some optimal solution of size <= |H|+l has all its hard edits either co-bagged pairs of covered vertices
   (in particular all A-local / type-T1 pairs) or edits at local hubs, then an FPT algorithm (time 7^(w+1) (Z^2+1)^(3 l (w+1)) n^O(1)) decides the instance;
   it is sound on every input. FAR with the path tiling violates the promise (CLASS-OPT = OPT+1, VERIFIED), FAR with the exchanged packing satisfies it.
5. **[PROVED] Structure of the far edits.** Every hard edit of an optimal solution lies in an induced P3 of G; if some witness has two soft pairs the edit is A-local,
   otherwise every witness has exactly one soft pair and the edit can be made soft by a size-preserving exchange (Lemmas 2.3, 2.1). So "far" = "absorbable by exchange".
6. **[OPEN conjecture; VERIFIED-BY-COMPUTER evidence] Conjectures with strong evidence (not proved).** C1 (normal form NF*): there is a *maximum* packing H* and an optimum all of whose hard edits are A-local or have an endpoint outside V(H*).
   VERIFIED-BY-COMPUTER on all 11,116 connected graphs with n=8 and all smaller ones, and on about 1,350 further structured / random instances with OPT > nu up to n=20 (plus ~1,000 hill-climbing evaluations). The naive stronger statement
   "all hard edits have an endpoint outside V(H*)" is FALSE (minimal counterexample K_2 v 3K_1, n=5; 305 of 11,116 graphs at n=8).
7. **[HEURISTIC] Precise residual obstruction.** C1 is an *existence* statement and is plausible even in a hard regime; the algorithmic part (finding H*) is the problem LOC
   (locate <= l far pairs among ~n^2 candidates with memory f(l,w)). Repacking does not remove this information problem, it moves it into H*: the new members of H* are
   exactly the witnesses of the far edits (Sec. 6). Bounded-radius augmentation cannot find H* (Prop 5.2: depth m is needed with l = 1, w = 1).

## 1. Setting and what "repacking" can and cannot mean

An instance (G,H,l) asks OPT <= |H|+l. For any packing H* with |H*| >= |H|, (G,H*,l*) with l* = |H|+l-|H*| <= l is an equivalent instance (Lemma MAX of the paper is the case H* = H + hard family).
So a repacking procedure may change both parameters: l decreases when |H*| grows, w becomes tw(I(H*)).
The intrinsic excess is g(G) = OPT - nu(G); l* >= g for every packing.

## 2. Exchange calculus (all PROVED)

**Lemma 2.1 (exchange).** Let H be a packing, Q in H, and T an induced P3 of G whose pairs avoid all pairs of H - Q. Then H' = H - Q + T is a packing, |H'| = |H|, and
I(H') - T is an induced subgraph of I(H) - Q; hence tw(I(H')) <= tw(I(H)) + 1, and j simultaneous exchanges give <= w + j. If T has two soft pairs, owned by two distinct members Q1, Q2,
then H - Q1 - Q2 + T is a packing of size |H|-1.
Proof. Pair sets are disjoint by assumption. Adding a vertex to a graph raises treewidth by at most 1 (add it to every bag). QED
If H is pair-maximal and T has exactly one soft pair s in Q, the hypothesis holds with that Q. If T has one soft pair, the other two pairs are hard.

**Lemma 2.2 (closed twins are never separated).** Let S be an optimal solution and N[a] = N[b] in G. Then a and b lie in one cluster of G xor S.
Proof. Suppose a in C1, b in C2, C1 != C2. Then ab in E (twins). Let c_a (resp. c_b) be the number of disagreeing pairs {a,x} (resp. {b,x}), x not in {a,b}. These pairs plus ab account for c_a + c_b + 1
of the cost. Moving a into C2 makes the pairs at a agree exactly as those at b (same adjacency to all x), and ab is a kept edge: the cost of these pairs becomes 2 c_b; symmetric move gives 2 c_a. All other pairs are unchanged.
min(2c_a, 2c_b) <= c_a + c_b < c_a + c_b + 1, contradicting optimality. QED

**Lemma 2.3 (witness dichotomy).** Let H be pair-maximal, S optimal, e in U a hard edit. Then e lies in an induced P3 T of G ("witness"), and every witness has exactly one or two soft pairs.
(a) If some witness has two soft pairs, they are owned by distinct members P != P' sharing the middle vertex of T; so e is A-local (type T1 of the paper).
(b) Otherwise every witness has exactly one soft pair s, owned by Q, and H - Q + T is a packing of the same size in which e is soft (Lemma 2.1); we call e *absorbable*.
Proof. Inserted e = ab: Lemma half(b) of the paper gives z with az, zb in E, and T = (a,z,b) is an induced P3. Deleted e = ab: by Lemma 2.2 a,b are not closed twins, so some x is adjacent to exactly one of them, and (x,a,b) or (x,b,a) is an induced P3 containing e.
Pair-maximality gives at least one soft pair; e is hard, so at most two. Two soft pairs share a vertex z (two pairs of a triple), their owners contain z; if the owners coincided, all three pairs would be soft, contradicting that e is hard. So P != P', P cap P' contains z, a in V(P), b in V(P'). QED

**Corollary 2.4.** A hard edit that is not A-local is absorbable. After absorbing, the pair e is soft and owned by a member containing both ends; the price is one more node (tw + 1), and the pairs of Q other than s become hard (they may carry edits of S: "displacement", Lemma 2.5).

**Lemma 2.5 (displacement, single-edit case).** In Lemma 2.3(b) assume e is inserted, z lies in the cluster of a,b (kept steps), the soft step is az in Q = {a,z,r}. Then Q contains exactly one edit of S, namely a pair f in {ar, zr},
and after the exchange f is the only new hard edit; its ends are r and a vertex of T. If r lies in no other member, f has an endpoint outside V(H') (it is "OUT"). Otherwise f may again be non-local ("the defect moves").
Proof. az is an edge inside a cluster, hence kept; Q = {a,z,r} is an induced P3 in G, exactly one of ar, zr is an edge; if r is in the cluster of a,z the non-edge is inserted and the edge is kept, else the edge to r is deleted and the non-edge is untouched: in both cases exactly one edited pair among ar, zr. QED (HEURISTIC remark: iterating gives an alternating chain of exchanges; termination in general is OPEN, Sec. 6.)

**Lemma 2.6 (hard components of outside vertices; PROVED).** Let H be pair-maximal and r not in V(H). Then the hard component of r is {r} union N_h(r), N_h(r) = N(r), and every pair inside it is a hard edge or a soft pair. In particular two vertices of this component that are not co-bagged are joined by a hard edge.
Proof. r lies in no member, so all pairs at r are hard. If z in N(r) and w in N_h(z), w != r, then rz, zw are hard edges, so by Lemma HE(a) of the paper rw is an edge or soft; it is not soft, so w in N(r). Hence the component is N[r]. For z, z' in N(r): zz' is an edge or soft (HE(a) with centre r). Soft pairs lie in one member, hence co-bagged. QED
(VERIFIED-BY-COMPUTER, inline check in this session: 2,006 outside vertices of random pair-maximal packings on 6..12 vertices, 0 violations of the closure property and of the pair property.)
Consequence (HEURISTIC use): an edit at an outside vertex r with far partners always sits in a small hard component spanning far members (|N(r)| <= D+l by Lemma EX); there are at most 2l exceptional clusters (Lemma NF), so at most 2l such components are touched. Hub non-locality is therefore the same phenomenon as far hard edits: a bounded number of small hard components, spread over far members, whose cluster status the DP must keep consistent.

## 3. The two known obstruction families, resolved by repacking

**Prop 3.1 (FAR, PROVED + VERIFIED).** FAR_m: path 0..2m plus chord e* = {1,2m-1}, H = {P_i = (2i,2i+1,2i+2)}.
(a) nu(FAR_m) = m: every induced P3 contains exactly two edges, pair-disjoint members are edge-disjoint, |E| = 2m+1. With OPT = m+1 (paper): the excess is intrinsic, l* >= 1 for every packing.
(b) H* = {P_0,...,P_{m-2}, T}, T = (1, 2m-1, 2m-2) (centre 2m-1; 1 and 2m-2 are non-adjacent for m >= 3) is a packing of size m = H - P_{m-1} + T, pair-maximal, and I(H*) is a cycle of length m (P_{m-2} meets T in 2m-2, T meets P_0 in 1): tw = 2.
(c) The partition {0,1},{2,3},...,{2m-2,2m-1},{2m} has m+1 edits; relative to H* its only hard edit is {2m-1,2m} (2m is in no member of H*); the chord is the soft pair of T.
So the far pointer of the original packing becomes a hub edit at the uncovered vertex 2m; Ball(2m) cap V(H*) = {2m-1,2m-2,1} lies inside the single member T, so 2m is a local hub in every decomposition and the optimum lies in the class of Thm PL (Sec. 4) for H*.
Proof. (a) shown. (b) induced: 1~2m-1 (chord), 2m-1~2m-2 (path edge), 1 !~ 2m-2. Pair sets disjoint from P_0..P_{m-2}: T's pairs are {1,2m-1},{2m-1,2m-2},{1,2m-2}; only P_{m-1} owned {2m-2,2m-1}. (c) edits: deleted path edges {2i+1,2i+2} in P_i (i <= m-2), chord in T, {2m-1,2m}. QED
VERIFIED-BY-COMPUTER: `w8a_farproof.py` (m = 3..14: validity, pair-maximality, edits, cycle, nu and OPT by ILP); `w8a_far2.py` (m = 3,4,5: all packings of size m are enumerated: 9, 11, 13; exactly one has a far hard edit in every optimum, namely the original tiling with distance m-1; all others have none; tw(I) in {1,2});
`w8a_farilp.py` (m up to 8, joint ILP over packings and optima: min number of far hard edits = 0).

**Prop 3.2 (Sp(F), PROVED + VERIFIED).** Let Sp(F) be the spine family of Appendix "spine" (any admissible F, N >= 1), H its packing, OPT = |H|+1 (explicit partition of the paper).
Put Y = (y,u_0,b_0) and Z = (b_0,u_1,c_1), H' = (H - X_0) + Y + Z with X_0 = (u_0,b_0,u_1). Then H' is a packing, |H'| = |H|+1; hence OPT = |H|+1 = |H'| (excess 0) and tw(I(H')) <= 2.
Proof. Y: y u_0, u_0 b_0 in E, y b_0 not in E (y is adjacent only to u_0, c_0), so Y is an induced P3; its pairs {y,u_0},{y,b_0} are hard and {u_0,b_0} lies in X_0 only. Z: u_1, c_1 lie in the clique K^+_{h(1)} with h(1) in V(F) (only e_0 has head *), so u_1 c_1 is a clique edge, hard (the soft clique edges are a_ib_i, b_ic_i, u_ib_i); c_1 has no edge outside its clique K^+_{h(1)} (the vertices c_i carry neither cross nor spine edges) and b_0 is not in K^+_{h(1)} (b_0 is in K^+_*, * != h(1)), so b_0 c_1 is a hard non-edge; {b_0,u_1} lies in X_0 only. So Z is an induced P3 whose pairs avoid H - X_0 and Y. OPT >= |H'| by Lemma decomp (any packing is a lower bound), and OPT <= |H|+1 by the explicit partition. I(H) is a caterpillar; I(H) - X_0 is a forest with P_0 isolated, and Y, Z attach as P_0-Y, P_0-Z, Y-Z, Z-X_1, Z-P_1 (X_1 - P_1 adjacent): a graph whose blocks are triangles and edges, tw = 2. QED
Consequence: Sp(F) is no obstruction to any algorithm that can repack, and the proof of "OPT >= |H|+1" no longer needs Rule 6.2 / the case analysis (S6).
VERIFIED-BY-COMPUTER `w8a_spine.py`: F = triangle, K4, K33, Petersen, 8 random cubic graphs on 10..40 vertices (n up to 246 vertices): all claims (valid, pair-maximal, explicit cost |H|+1, H' valid of size |H|+1, tw(I(H')) <= 2 by min-fill) hold. `w8a_nstar.py`: nu_mod = |H|+1 = OPT (ILP) for N*(triangle) (n=17) and N*(K4) (n=29).
In addition **(PROVED)** the original H of Sp(F) is in the class of Thm PL: the optimum of the paper has the single hard edit {y,b_0}; y is outside V(H); Ball(y) cap V(H) = {u_0,c_0,b_0}, contained in V(X_0) cup V(P_0), and X_0, P_0 share b_0, so they are co-bagged: y is a local hub. VERIFIED-BY-COMPUTER `w8a_plspine.py` (CLASS-OPT = OPT = |H|+1 for triangle, K4). So the "trigger located at an unknown position" of wave 6 is not a problem for the DP of Thm PL: the position is chosen bag-locally.

## 4. Promise-FPT: Theorem PL

Define the class 𝒮(H,T,beta) of solutions S, |S| <= |H|+l, whose hard edits are (i) pairs {a,b}, a,b in V(H) in L_B, b in Ball(a), co-bagged, or (ii) pairs {r,v}, r outside V(H) a local hub, v in V(H) cap Ball(r).

**Theorem PL (PROVED; DP details routine, semantic check VERIFIED).** Let H be pair-maximal and (T,beta) a tree decomposition of I(H) of width w. There is an algorithm of running time 7^(w+1) (Z^2+1)^(3 l (w+1)) n^O(1)
(Z = 3D+3l+1 <= 6(w+1)+3l+1) that outputs a number m* >= OPT(G) with m* = OPT(G) whenever some optimal solution of size <= |H|+l lies in 𝒮(H,T,beta). In particular it never accepts a NO instance (one-sided), and it is complete if all hard edits of some optimum are A-local (type T1 included) or hub edits at local hubs.

Proof (DP and correctness). State at a bag t: (x_P) in [7]^(beta_t) (the 7 non-empty subsets of pairs(P)), and for every covered vertex v in X_t cap L_B a label E_v subset of Ball(v) of hard partners with |E_v| <= l (E_v may contain outside vertices r that are local hubs). Transitions as in Lemma dp: a vertex label is chosen when v is introduced, is carried unchanged through T_v (connected), forgotten with v; join nodes require equal labels. Cost = sum of |x_P| + (1/2) sum_v |E_v cap V(H)| + sum_v |E_v - V(H)|. Local conditions checked at a bag:
(S) symmetry: for v,v' in X_t cap V(H): v' in E_v iff v in E_{v'};
(V1) every vertex triple Delta with a soft pair: its soft pairs' owners are pairwise adjacent in I(H) (Observation two), hence (Helly) lie in one bag t; two vertices a,b of Delta lie in X_t; the hard pairs of Delta are {a,c},{b,c}, with status G xor [c in E_a], resp. [c in E_b], both determined at t; soft pairs by (x_P). Check "not an induced P3 after editing";
(V2) every triple with no soft pair and a U-pair {a,b}: if a,b in V(H) they are co-bagged, take a bag containing both; the statuses of {a,c},{b,c} are [c in E_a],[c in E_b]; if a = r is a hub, all partners of r lie in Ball(r) cap V(H), pairwise co-bagged, so (Helly) in one bag t_r; there every pair {r,c} with c in Ball(r) has status [r in E_c] visible, and {r,c} with c outside Ball(r) is not editable (unedited), {v,c} has status [c in E_v];
(V3) a triple with no soft pair and no U-pair is unedited, hence not an induced P3 because H is pair-maximal.
Soundness: a consistent labelling over all bags defines S = U union (union of x_P), U = {vp : p in E_v}; by (S) the U-pairs between covered vertices are well defined and counted once; (V1)-(V3) verify every triple, so G xor S is a cluster graph and cost(S) is the DP value. Completeness: for S in 𝒮 take E_v = {p : vp in U} (partners lie in Ball(v) by Cor cand for optimal S, are co-bagged / hub-local by definition), x_P = S cap pairs(P); all checks hold because G xor S is a cluster graph. Number of states: 7^(w+1) for members, at most 3(w+1) vertices in X_t, each with <= (Z^2+1)^l labels. QED
Note: soundness needs the hub-locality precisely because an outside vertex r has no subtree of bags: two edits at r with far partners could not be matched.

VERIFIED-BY-COMPUTER (`w8a_pl.py`, semantics of the visibility argument, not the DP itself): (1) random solutions S in the class (2,400 trials, random instances n = 7..10, min-fill decompositions): local verification agrees with global validity in 100% (0 mismatches); (2) 243 optimal solutions lying in the class were all accepted by the local verification;
(3) FAR_m, m = 3..6: with the path tiling the class excludes the chord and CLASS-OPT = OPT + 1 (promise fails, algorithm incomplete); with the exchanged packing H* of Prop 3.1, CLASS-OPT = OPT. (4) Sp(triangle), Sp(K4): CLASS-OPT = OPT with the original H. (5) `w8a_plrand.py`, 265 random instances (n = 7..10, random pair-maximal H, min-fill decomposition, l = OPT-|H| <= 2): CLASS-OPT >= OPT always (soundness; infeasible class counted as infinity); CLASS-OPT = OPT in all 253 instances where some optimum lies in the class, and CLASS-OPT > OPT in the 12 instances where none does (the promise failed for a random H in 12 of 265 = 4.5%, all with l >= 1).
Remark 4.0 (relation to the omega theorem of the paper, PROVED by the examples). Thm PL and Theorem omega (FPT in (l,D,omega), omega = tw of the vertex graph Gamma') are incomparable. Sp(F) has tw(Gamma') >= tw(F) unbounded (expander F) but lies in the PL class with w = 1, l = 1; FAR has tw(Gamma') = 2 (a cycle) but with the path tiling is outside the PL class (it enters it after one exchange).
Corollary 4.1 (static form, PROVED from PL). If a given decomposition co-bags every candidate pair (hard pairs {a,b}, a,b in L_B cap V(H), b in Ball(a), and the hub pairs), the promise holds for every instance, so CEaMP is solvable in f(l, w_Lambda) n^O(1), where w_Lambda is the width of a decomposition of I(H) that co-bags all candidate pairs (equivalently the treewidth of I(H) plus a link between the members of the two ends of every candidate pair). This is the decomposition-based form of the "Lambda" conjecture of wave 6 (stated there as open); it is only as strong as tw(I(H) + all candidate links), which is large on tunnel-rich instances, so it does not answer the main question.

## 5. Impossibility results and counterexamples

**Theorem N (PROVED, conditional on P != NP resp. ETH).** (a) If a polynomial-time algorithm maps every YES instance (G,H,l) to a packing H* with excess l* <= s(l) and tw(I(H*)) <= f(l), for any functions s, f, then P = NP.
(b) Under ETH there is no polynomial-time algorithm mapping every YES instance with l = 0 to a packing H* with |H*| = |H| and tw(I(H*)) <= |H|^(1-eps).
Proof. (a) At l = 0 the problem is NP-hard even with every vertex in <= 23 members (LPS23, cited in the paper). Given (G,H,0) run the algorithm; (G,H*,l*) is equivalent (same budget |H| = |H*|+l*); solve it by Thm A in time 7^(f(0)+1) n^(2 s(0)+O(1)), polynomial; a tree decomposition of width f(0) is found in linear time for constant width. (b) Thm l0 gives 3^(tw+1) n^O(1) = 2^(O(N^(1-eps))) for 3-SAT instances with |H| = O(N) (proof of Thm eth), contradicting ETH. QED
So the target of the task ("output H* with tw(I(H*)) <= f(l)") is impossible in this form; any usable repacking must keep w* <= w + g(l), and its output must come with a locality property (otherwise H* = H works).
The hope "if tw(I(H)) is large relative to l then a better packing exists, hence reject" is false for the same reason: every tight instance (l = 0, nu = OPT = |H|) has no better packing, and tight instances have unbounded tw under ETH.

**Prop 5.1 (the naive normal form is false; VERIFIED-BY-COMPUTER).** Let G = K_2 v 3K_1 (hubs 3,4 adjacent, leaves 0,1,2 each adjacent to both hubs; n = 5). nu = 2 (proved: the induced P3s are (x,h,y), x,y leaves; a hub supports at most one member, so nu <= 2, and a packing of size 2 uses leaf pairs {x,y},{x',y'} that are different, so covers all 5 vertices), OPT = 3, four optima. Every maximum packing covers all vertices and, with every optimum, has a hard edit with both ends covered (ILP, `w8a_ilp.py`). So "every hard edit has an endpoint outside V(H*)" fails; the surviving hard edit (e.g. {1,4}, deleted) is A-local. Counts of graphs where it fails (OUT-failures among connected graphs with at least one P3): n=5: 3 of 20; n=6: 2 of 111; n=7: 3 of 852; n=8: 305 of 11,116 (`w8a_atlas.py`, `w8a_exh.py`).

**Prop 5.2 (bounded-radius augmentation is not enough; PROVED + VERIFIED).** Path P on vertices 0..2m+2, H_0 = {(2i+1,2i+2,2i+3) : 0 <= i < m}. H_0 is a pair-maximal packing (the hard edges {0,1}, {2m+1,2m+2} are not consecutive), nu = m+1 = OPT, l = 1, w = 1. The unique packing of size m+1 is {(2i,2i+1,2i+2) : 0 <= i <= m} (it must use all 2m+2 edges; edge {0,1} forces (0,1,2), then {2,3} forces (2,3,4), ...), and it shares no member with H_0. So every improving exchange (remove r members, add r+1) has r = m. (`w8a_augchain.py`, m = 2..6: r = m.) Hence iterative compression / k-exchange local search with k bounded in l cannot find the augmentation; a long alternating path is required.
r-exchange local search on random instances (`w8a_aug.py`, seeds 1-3, n = 9..15, random pair-maximal H): of 169 instances with nu > |H|, r=1 reaches nu in 97 and r=2 in 145; 24 stay open at r=2 (VERIFIED-BY-COMPUTER).

**Remark 5.3 (neutral exchanges are essential).** In FAR_m the original tiling is already a maximum packing; no improving exchange exists. The helpful move (Prop 3.1) is *neutral* (size-preserving) and its benefit is invisible to any objective computed from (G,H) alone: it is chosen because the chord is the hard pair that the optimum cuts. This is the reason why Sec. 6 reduces repacking to locating hard edits.

**Procedure REPACK-XP (what can be stated as an algorithm; sound, completeness HEURISTIC).** Input (G,H,l).
1. Maximalise H (Lemma MAX); reject if more than l hard P3's are added; l := l - |M|, w := w + |M|.
2. Candidate far pairs: Cand = hard pairs {a,b}, a,b in L_B cap V(H), b in Ball(a), not co-bagged in the given decomposition (plus hub pairs whose hub is not local); for each such pair and each witness T with exactly one soft pair (<= n^2 Z^2 choices) the exchange H - Q + T of Lemma 2.1.
3. For every set of at most l pairwise compatible exchanges (n^O(l) sets; after each, maximalise and update l) run the DP of Thm PL on the new packing and its min-fill decomposition; accept iff some run returns a value <= |H|+l.
Soundness: PL is sound for every packing; the budget is unchanged (|H*| >= |H|). Running time n^O(l) f(l,w): XP, not FPT (it guesses the exchanges). Completeness holds whenever C1' holds with H* reachable by <= l such exchanges (single-step displacement as in Lemma 2.5 ends at an uncovered endpoint or an A-local pair); this covers FAR_m (one exchange), Sp(F) (original H), and every instance of the tests of Sec. 7 (where a good H* exists), but it is not proved. It is not asymptotically better than Thm A (both n^O(l)); its point is structural: the guessed objects are exchanges rather than hard pairs, and the DP never has to carry a far pair.
Replacing step 3 by an FPT-time choice of the exchanges is exactly the open problem C2.

**Answers to the three suggested routes.**
(i) Augmenting / rerouting along far chords: works for every far edit individually (Lemma 2.3: not A-local implies absorbable), resolves FAR (neutral exchange) and Sp(F) (exchange + augmentation); what it needs is the identity of the chord, i.e. the optimum.
(ii) Local-ratio / exchange argument "large tw implies a better packing": false (Thm N remark: tight instances); the true statement is only "a far edit can be absorbed".
(iii) Iterative compression on H: the compression step is an augmentation of unbounded depth (Prop 5.2) or a neutral exchange (Remark 5.3); neither is detectable from (G,H) with bounded search radius.

## 6. Why this does not give an algorithm: the information barrier (HEURISTIC) and the exact residual problem

Candidate for the full statement, for pair-maximal H: 
* **C1 (existence, normal form NF*), OPEN, strong evidence.** For every graph there is a maximum packing H* and an optimal solution such that every hard edit has an endpoint outside V(H*) or is A-local. Stronger variant NF** (also: all edits at one uncovered vertex have pairwise A-local-or-soft partners) is what Thm PL needs for hubs.
  Evidence (Sec. 7): exhaustive for all connected graphs n <= 8 (NF*), n <= 7 (NF**) [n = 8: half of the graphs; n = 9: random samples; see Sec. 7], ~1,350 further instances.
* **C1' (width).** In addition tw(I(H*)) <= w + g(l) and |H*| >= |H|, i.e. H* is obtained from H by exchanges along <= g(l) "far links". (OPEN; tested in Sec. 7 as RL.)
* **C2 (algorithm).** H* can be computed in time f(l,w) n^O(1). OPEN.
By Thm PL, C1' + C2 (with the decomposition of I(H*)) imply FPT in (l,w). C1 alone is not informative: it is plausible that C1' holds also for instances that are W[1]-hard, because H* can be built from an optimum (absorb its far hard edits, Lemma 2.3); the entire difficulty is C2, which has to produce the far edits without knowing the optimum.
Reduction (HEURISTIC but elementary): if S is an optimum with far edits e_1..e_k (k <= l) and H* absorbs them, then H* - H contains the witnesses T_i of the e_i, each T_i adds a node adjacent to the members through its three vertices; conversely from H*, H one reads off the e_i as the hard pairs inside the new members. So computing H* is computing the far edits.
Why the DP cannot find them (HEURISTIC): the DP of Thm PL must carry, for every far pending edit, the identity of its partner (a vertex among ~n, or among Z^2 candidates relative to a vertex that has already been forgotten) across the part of the decomposition between the two members; there are up to n potential far pairs crossing a separator, <= l of them active, so f(l,w) states cannot represent which. Hashing / colour coding of the endpoints fixes the completeness but not the soundness: two vertices of equal colour that look locally alike produce a fake configuration that is cheaper than every real solution (the DP counts one edit while the checks assumed two different pairs).
The same information appears for hubs: an outside vertex r with two edits whose partners are far apart cannot be tracked (no bag contains r); the local-hub restriction in Thm PL is exactly where soundness needs it.
What is NOT excluded: that more structure (for instance zero-excess solution sets closed under a majority operation, which empirically holds for c <= 2 in wave 6c) makes the far edits recoverable.
Hardness target (HEURISTIC, no reduction): a W[1]-hardness proof needs a bounded-width tight region with N clean sockets whose locked sets are exactly {S : S != all}, i.e. a prefix-OR bit carried across separators (the persistent-flip chain of wave 6c); then n potential far links (tunnels) crossing every separator with <= l cut would encode an index selection. None is known; the repacking viewpoint makes no statement about it.

## 7. Experiments (all VERIFIED-BY-COMPUTER; evidence, not proof)

| id | script | what | scale | result |
|---|---|---|---|---|
| E1 | w8a_farproof, w8a_far2, w8a_farilp | FAR repack | m=3..14 (explicit), m=3..5 all packings, m<=8 ILP | as Prop 3.1; only the tiling has a far edit |
| E2 | w8a_spine, w8a_nstar | Sp(F), N*(F), expanders | n up to 246 | one exchange + one augmentation, OPT=nu=|H|+1, tw<=2 |
| E3 | w8a_atlas, w8a_exh (n=8), w8a_exh3 | NF* / NF** exhaustive | all connected graphs n<=7; n=8: 11,116 graphs (4,413 with OPT>nu; OUT-failures 305; NF* (soft distance <= 2 for covered-covered hard edits): 0 violations); NF** n<=7 exhaustive, n=8 half of the graphs: 0 violations | NF*, NF**: no counterexample |
| E4 | w8a_exh3s | uniform random samples of connected graphs n=9 (2 x 350 of 261,080) | NF** | 0 violations (323 with OPT > nu) |
| E5 | w8a_ilp, w8a_rand, w8a_out, w8a_cov | random instances (gnp, tree+chords, cycle+chords, cubic, planted), n=7..18; OUT property | ~800 instances | OUT holds in all of them, fails only in tiny dense graphs (E3) |
| E6 | w8a_struct, w8a_struct3 | structured: paths/cycles/ladders/grids/caterpillars/trees with chords, planted, gnp; NF* (r=2) and NF** | 800 + 600 instances, n<=20 | 0 violations |
| E7 | w8a_hill, w8a_mut | hill-climbing / skeleton mutation (around FAR) maximising violations, r=2 and r=3 | n=9..12, ~600 evaluations | best = 0 |
| E8 | w8a_circ, w8a_theta | circulants, prisms, Moebius, wheels, Petersen, cube, dodecahedron; theta graphs | n<=20 | 0 violations |
| E9 | w8a_kappa, w8a_chain, w8a_chords | number of new members needed so that the far edits disappear; path + chord | random pair-maximal H, n=9..13; canonical tilings m=5..7 | path+chord: always 1 exchange; random H with forced far edit: kappa in {1,2,3}, kappa = 3 occurs with l = 1 (the input H is far from maximum) |
| E10 | w8a_aug, w8a_augchain | r-exchange local search | random n=9..15; chain family | chain family needs r = m (Prop 5.2) |
| E11 | w8a_pl | semantics of Thm PL | see Sec. 4 | 0 mismatches |
| E12 | w8a_rltw, w8a_plrepack | RL: min tw(I(H*)) over repackings that remove the far edits / restore the promise vs tw(I(H)) + l | random, n=8..12 | width never above tw(I(H)) + 1 (details below) |

Additional numbers (VERIFIED-BY-COMPUTER).
* NF** (A-local covered-covered hard edits and hub condition; `w8a_exh3.py`, `w8a_exh3s.py`): exhaustive for all connected graphs with n <= 7; n = 8: 5,415 graphs (four of twelve geng shards, about half of the 11,116 connected graphs), 2,710 with OPT > nu: 0 violations; n = 9: two independent uniform random samples of 350 connected graphs each (population 261,080), 323 with OPT > nu: 0 violations. NF* (without hub condition) is exhaustive for all connected graphs n <= 8 (E3).
* E12 (`w8a_rltw.py`, seed 1, 33 random instances with a forced far edit and l >= 1): min over packings H* with |H*| >= |H|, pair-maximal and no far edit of tw(I(H*)) never exceeded tw(I(H)) + 1 (so within tw(I(H)) + l since l >= 1); it was <= tw(I(H)) in 21 of the 33 and exactly tw(I(H)) + 1 in 12 (exact treewidth by elimination-order DP).
* E12' (`w8a_plrepack.py`, seeds 1-3, 71 random instances, n = 8..12, whose promise (Thm PL class, min-fill decomposition) FAILED for the random input H): in every one there is a packing H*, |H*| >= |H|, pair-maximal, for which the promise holds; min-fill upper bounds give tw(I(H*)) <= tw(I(H)) in 70 of 71 and tw(I(H)) + 1 in the remaining one (all with l >= 1). So in these tests repacking restores the promise without paying more than +1 in width (and the number of changed members was up to 3: kappa of E9).
* E7' (`w8a_hill3.py`, hill-climbing on the NF** objective with hub terms, n = 11, 12): best violation count 0.
* `w8a_hub.py`: spine variants in which the trigger y is adjacent to far vertices c_1, c_3 etc. are not pair-maximal (the P3's (y,c_i,m) are all-hard); after maximalisation they are tight (OPT = nu), i.e. far hubs adjacent to part of a far-spread hard clique are removed by pair-maximality, in line with Lemma 2.6.


## 8. Status table

| claim | status |
|---|---|
| Lemma 2.1 exchange + width; Lemma 2.2 twins; Lemma 2.3 witness dichotomy; Cor 2.4 | PROVED |
| Lemma 2.5 displacement (single-edit case) | PROVED; iteration/termination OPEN |
| Prop 3.1 FAR: nu = m, OPT = m+1, H* with cycle I(H*), hard edit at an uncovered vertex | PROVED; VERIFIED m <= 14 |
| Prop 3.2 Sp(F): one exchange + one augmentation, OPT = |H|+1 = |H'|, tw <= 2; original H in class PL | PROVED; VERIFIED up to n = 246 |
| Thm PL (promise-FPT, sound always, complete on the class) | PROVED (DP details routine); semantics VERIFIED |
| Thm N (no repacking with l-only bounds unless P = NP / ETH) | PROVED (uses LPS23) |
| Prop 5.1 OUT normal form false (K_2 v 3K_1, 305 graphs at n = 8) | PROVED/VERIFIED (nu proof; OPT, counts by computer) |
| Prop 5.2 augmentation depth m with l = 1, w = 1 | PROVED; VERIFIED m <= 6 |
| C1 normal form NF* (A-local or uncovered endpoint) | OPEN; VERIFIED exhaustive n <= 8 and ~1,350 further instances |
| C1' + C2 (repack within w + g(l), in FPT time) => FPT in (l,w) | OPEN (implication PROVED via Thm PL) |
| Information barrier / repacking = locating far edits | HEURISTIC |
| FPT in (l,w) | OPEN |
| W[1]-hardness in (l,w) | OPEN (no reduction; needs prefix-OR chain) |

## 9. Next steps (suggestions for a later wave)
1. Prove NF* for intrinsic excess 1 (OPT = nu+1): the exchange chain of Lemma 2.5 is an alternating-path walk on the edits of one optimum; a termination / exit argument (reach an uncovered endpoint or an A-local pair) would be a PROVED normal form.
2. Prove that hub edits with far partners force a far covered-covered edit (Lemma 2.6 is the first step), so that the hub clause of the promise in Thm PL can be dropped.
3. The only known route to W[1]-hardness: a bounded-width tight region with N clean sockets whose lockable sets are exactly the proper subsets (prefix-OR chain), combined with tunnels; the experiments here (every instance repackable, NF* exhaustive for n <= 8) give no hint that such a chain exists for c <= 2.
4. Combine Thm PL with the majority-closure observation of wave 6c for c <= 2.

## 10. Files
`code/w8a_lib.py` (graphs, induced P3s, packings, exact OPT by subset DP / ILP, nu by ILP, treewidth), `w8a_ilp.py` (joint packing x optimum ILP: OUT), `w8a_ilp2.py` (adds soft-distance locality), `w8a_ilp3.py` (adds hub condition),
`w8a_far*.py`, `w8a_spine.py`, `w8a_nstar.py`, `w8a_pl.py`, `w8a_plspine.py`, `w8a_exh*.py`, `w8a_atlas.py`, `w8a_struct*.py`, `w8a_hill.py`, `w8a_mut.py`, `w8a_circ.py`, `w8a_theta.py`, `w8a_hub.py`, `w8a_plrand.py`, `w8a_plrepack.py`, `w8a_exh3s.py`, `w8a_kappa.py`, `w8a_chain.py`, `w8a_chords.py`, `w8a_aug*.py`, `w8a_rltw.py`, `w8a_rand.py`, `w8a_out.py`, `w8a_cov.py`.
