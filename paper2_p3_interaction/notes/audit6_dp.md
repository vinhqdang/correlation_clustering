# Audit 6 (hostile referee, wave 5): notes/wave5_dp_route.md

Target: Lemma 1, Lemmas 2-3, Lemma 5, Cor. 6, Thm 7, Thm 8, Prop 9. Context read: audit4_maxhubsm.md, audit5_locality.md,
paper/sections/packing.tex, w5b_*.py, fpt_attempt.md (for Thm A). Own code, written from scratch (imports only a6_lib.py):
code/a6_lib.py, a6_t7.py, a6_t8.py, a6_p9.py. No git, paper/ untouched.

## Verdicts
| item | verdict |
|---|---|
| Lemma 1 (cost decomposition a, b, c) | CORRECT |
| Lemma 2 (HE, H2) | CORRECT |
| Lemma 3 (d_h bounds, |C| <= 2D+2l+1, d_h(y) <= 2D+3l) | CORRECT (sharper |C| <= 1+2D+|U| follows from the same proof) |
| Lemma 4 (clusters of an optimum are G-connected) | CORRECT |
| Lemma 5 (alpha, beta, gamma normal form) | CORRECT |
| Cor. 6 (bounded partners) | CORRECT |
| Thm 7 (n^l g 7^{w+l+1}) | CORRECT WITH FIX (algorithm is correct; three repairs below; the root/~-component machinery is unnecessary) |
| exponent halving n^{2l} -> n^l | CORRECT (genuine, and follows from Cor. 6 alone) |
| Thm 8 (FPT in l, D, omega) | CORRECT WITH FIX (one false sentence in the proof; DP description omits the join term) |
| Prop 9 (l=0, w=0, c=1, tw(Gamma') >= tw(F)) | CORRECT WITH FIX (|K_v| <= 6 not 5; hard degree <= 4 is right and is tight against B) |
Nothing is WRONG. No brute-force counterexample (n <= 10) to any stated bound.

## Lemma 1 -- CORRECT
(a) cost = sum_C [C(|C|,2) - e(C)] + [|E| - sum_C e(C)]. Right. (b) soft pairs are exactly the disjoint union of the 3-pair sets of the members, so
cost = sum c_i + |U|, c_i >= 1 because the triple is still an induced P3 if no pair of it is edited. Right. (c) pure bookkeeping; far hard pairs
are edited hard pairs, so far(pi) <= |U| <= l. Right. Only nit: width 3w+2 needs H nonempty (audit4 Sec. 5).

## Lemmas 2, 3 -- CORRECT
Lemma 2: both are the pair-maximality contradiction. Lemma 3(a) re-derived line by line.
* deleted uv: u in C, v not in C. A kept hard neighbour z of u lies in C, so zv is cut if it is an edge; by Lemma 2(i) (centre u) zv is a
  hard edge (hence a deleted edge at v, other than uv) or soft. Count: kept <= (del_v - 1) + D. d_h(u) = kept + del_u. OK.
* inserted uv: kept z in C; zv is an edge xor pi inside C, hard edge impossible by Lemma 2(ii); so soft (<= D) or inserted at v (<= ins_v - 1). OK.
* (b) all w in C-u: uw soft, kept hard edge, or inserted: |C| <= 1+D+d_h(u)+ins_u. (c) hard neighbours outside C are cut. OK.
Sharper (derived here, same proof): in the deleted case del_u+ins_u+(del_v-1) <= |U| (distinct pairs at u or v), in the inserted case
ins_u+(ins_v-1)+del_u <= |U|, so |C| <= 1+2D+|U| (<= 2D+l+1). Tested, 0 violations (slack -1 .. -4 at n <= 10).
Lemma 5: (alpha) fine; (beta) |C| <= 2D for >= 2 components: |C| = |C-K1|+|K1| <= D + |C-K2| <= 2D. Fine. (gamma): K(v) has > B+1 >= 2D vertices, v is
neither an endpoint (d_h <= D+l <= B) nor in an exceptional cluster (d_h <= B), so v is in a beta cluster that cannot hold two components. Fine.
"Not clean => NO" holds because an optimal pi would otherwise edit a hard pair inside a non-exceptional cluster. Fine.

## Cor. 6 and Theorem 7 -- algorithm is correct
How the free endpoints are handled: they are NOT identified inside a DP. Thm A (fpt_attempt, audited CORRECT in audit.md) takes U as an explicit
guess and then checks every all-hard triple directly (n^3 per U) and all triples with a packed pair as CSP constraints on cliques of I(H).
So for ANY guessed U' the DP is exact, far-apart endpoints cost nothing, and soundness is automatic (a DP value is the cost of a real partition).
The only thing to prove is completeness: some optimal pi has U inside the candidate family. That is Cor. 6 (one endpoint determines a ball
containing the other) plus Lemma 5 (alpha). I re-derived: endpoints lie in L_B (d_h <= D+l <= B); ~-steps (same cluster / U-pair) have
G[L_B]-distance <= |C|-1 <= 2D+2l (Lemma 4 + Lemma 3(c)); a ~-component has <= 2|U| vertices and every endpoint lies within (2l-1)(2D+2l) of
any root of its component; #components <= |U| <= l. Roots in L_B: <= n^l. U' inside the union of balls: (l |ball|)^{2l}. All correct.
Fixes / remarks:
1. U' must be restricted to HARD pairs of H' (Thm A enumerates unpacked pairs only; the "union of balls" contains soft pairs). One word.
2. The ~-component and radius R=(2l-1)(2D+2l) argument is unnecessary: Cor. 6 gives directly "guess ONE endpoint u of each U-pair (n^{|U|}, u in L_B),
   partner in Ball_{G[L_B]}(u, 2D+2l) (deleted: N_h(u), size <= D+l)". This is the exponent halving, with ball radius 2D+2l instead of R, so
   g(l,D) = ((B+D+1)^{2D+2l})^l-type, much smaller than (l (B+D+1)^R)^{2l}. The note's route is valid but the simple route is stronger.
3. Bookkeeping with MAX (audit4 Gap 2): D, B, l, w in the statement are the primed ones: D' <= D+2|M|, l' = l-|M|, width <= w+|M|; the bound
   n^{l'} 7^{w+|M|+1} is <= n^l 7^{w+l+1}, but g must use D' <= D+2l. Say so; also a tree decomposition of I(H') of width w+|M| is assumed/derived from one of I(H).
Computer (a6_t7.py; all set partitions, n <= 10; 1,511 instances, 8,988 optimal partitions at budgets l = OPT-|H| and +1; 6,426 edited pairs):
|U| <= l, d_h(endpoint) <= D+l, |C| <= 2D+2l+1 and <= 1+2D+|U|, endpoints in L_B, deleted pairs at distance 1, inserted pairs at G[L_B]-distance
<= 2D+2l (observed max 2, as in w5b), every ~-component inside the radius-R ball of any of its vertices, #components <= |U|: 0 failures.
(The "mutation" counter printed there is trivial; n <= 10 cannot stress the radius. This is evidence only; the proofs carry the claim.)

## Theorem 8 -- CORRECT WITH FIX
Equivalence (<=): any partition of V built from pi' and the BIG clusters has cost |E|+sum f (Lemma 1a), an upper bound on OPT. Right even without
connectivity or size restriction. (=>) needs Lemma 5. Checked pieces:
* Z = 3D+3l+1: exceptional <= 2D+2l+1, multi-component beta <= 2D, one clean non-BIG component K: each y has all other vertices of K adjacent or soft, so
  |K| <= d_h+D+1 <= B+D+1 = Z. Right (observed slack -2).
* Gamma'-connectivity: G-edges and soft pairs between V' vertices. Exceptional: G[C] connected; beta with one component: hard-connected; several
  components: all cross pairs soft. Right.
FIX 1 (a false sentence). The proof says "hard edges lie in Gamma' because both ends are in L_B subset V'". L_B is NOT a subset of V' = V - BIG:
a BIG hard component can contain vertices of small hard degree (pendant vertices of a hub component). a6_t8.py found 30 such (instance, l) pairs.
Repair, one line: by (gamma) every BIG component is a whole cluster, hence non-exceptional, so an exceptional cluster is disjoint from every BIG
component; therefore C inside L_B minus BIG inside V'. (Checked: exceptional cluster meets BIG in 0 of all cases.)
FIX 2 (DP description). "Introducing v adds s - 2k" is right (k = G-neighbours of v in its block inside the bag; a forgotten vertex cannot be a
G-neighbour of a later vertex because every G-edge is in Gamma'). But the DP is not additive at join nodes: for a block with bag part b and
exclusive parts x1, x2 in the two children, f(C) = f(C1)+f(C2)-f(bag part)+x1*x2 (no G-edges between the exclusive parts). The cross term
x1*x2 and the subtraction of f(bag part) must be stated (sizes are tracked, so both are computable). I implemented the full DP
(a6_t8.py, introduce/forget/join over a min-fill tree decomposition, per connected component): its value m satisfies
m_all <= m_DP <= m_conn (brute force over all size-<=Z partitions / over Gamma'-connected ones), and the omitted-join-term DP differs from the correct one on
~6% of the instances (134 of 2,141 DP checks), so the omission is not harmless if someone codes the note literally.
Connected-subtree claim: correct. If the last bag vertex of a cluster C is forgotten at t, then T_C (union of the vertex subtrees, connected since C is Gamma'-connected)
cannot contain the parent of t, so T_C lies below t and all of C is already introduced. The DP does not need to ENFORCE connectivity (the family of
realised partitions is a superset of the connected ones and every member is a real partition with correct cost); connectivity is used only to justify closing.
Computer: 2,778 (instance, l) pairs, n <= 10: formula  |E|+sum_{BIG} f(K)+m <= |H|+l  <=> OPT <= |H|+l in all cases (m by brute force);
every non-BIG cluster of every optimal partition is Gamma'-connected and <= Z; BIG components are whole clusters (39 non-trivial BIG clusters seen, few
because n <= 10; evidence only). 0 failures. The caveat stands: omega, Z and D are those of the primed instance.

## Proposition 9 -- CORRECT WITH FIX
Construction re-derived from the text. P3 (a_e,b_e,c_e): a_e b_e cross edge, b_e c_e inside K_v, a_e c_e non-edge: induced. Vertex-disjoint, hence pair-disjoint,
w = 0, c = 1, D = 2. Pair-maximal: every vertex has at most one cross edge, every induced P3 must use a cross edge a_e b_e (cliques contain none), and that pair is soft.
Partition {K_v} deletes exactly the cross edges: cost = |H|; OPT >= |H| by the packing; so OPT = |H|, YES with l = 0. Each K_v is a single clean hard component
(only b_e c_e is soft inside it). Contracting the K_v gives F (soft pairs a_e c_e, b_e c_e add no new adjacency), so F is a minor of Gamma' and tw(Gamma') >= tw(F) = Omega(|F|) for
cubic expanders (Grohe-Marx 2009 CITED, not re-checked).
Fix: |K_v| = (out-degree)+2(in-degree) = 3+indeg <= 6, not 5 (6 occurs when all three incident edges point into v; note the code docstring says "<= 5" too).
The claim "hard degree <= 4 = 2D+3l, no BIG" is right and tight: max d_h = 4 = B at l = 0 (a_e in a 5-clique or b_e in a 6-clique), so no BIG vertex arises, but one more
unit would change the family.
Computer (a6_p9.py): exhaustive over all set partitions for the four F with 3 edges (n = 9), both orientations: OPT = |H| = 3, induced, pair-disjoint, pair-maximal;
K4 (n=18), K33 (27), Petersen (45), random cubic 20 and 40 nodes (n = 90, 180): induced, pair-disjoint, pair-maximal, cost({K_v}) = |H|, K_v = hard components and clean,
quotient of Gamma' by the K_v is exactly F, D = 2, max |K_v| = 6, max d_h = 4. Not checked by me: the note's side claim that BFK Rule 6.2 trivialises the family (rule not in the provided files).
Scope remark (the note says it, I repeat it): this instance is l = 0, so Thm A solves it in polynomial time; Prop 9 refutes only "omega <= f(l,w)", not a route to hardness.

## Residual gaps not raised by the note (no counterexample found)
* Thm 7 and Thm 8 are stated for H pair-maximal AND with a tree decomposition of the relevant graph given; computing one is FPT (Bodlaender) but say it.
* The statement of Thm 7 should say: guess one endpoint per U-pair (simplified version, item 2); the note's version is correct but heavier.
