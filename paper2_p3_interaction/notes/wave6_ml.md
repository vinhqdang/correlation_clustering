# Wave 6 (prover): the reduction lemma ML (conjecture of dproute.tex / wave5_dp_route.md Sec. 6)

Labels: PROVED (proof given here or in the cited note, which was audited), CITED (BFK18, Grohe-Marx, taken from the literature),
VERIFIED-BY-COMPUTER (finite check, evidence only, never a proof), CONJECTURE.
Code (all own, independent of earlier waves): code/w6a_*.py, code/w6a_opt.c (exact OPT by subset DP in C; compile with `gcc -O2 -o w6a_opt w6a_opt.c`).
No git, paper/ untouched.

## 0. Verdict

ML is **neither proved nor refuted**. What I can say:

1. A counterexample to ML cannot be a single instance family. ML only asserts that SOME polynomial safe reduction exists. A family can only
   defeat a stated class of rules. I defeat the natural class (BFK Reduction Rule 6.2 applied to unions of members, Sec. 4) with a family of
   treewidth w = 1, budget l = 1 (Prop. 5). So the obvious generalisation of BFK for c = 1 cannot prove ML, already for w = 1.
2. New PROVED partial result (Thm 3): a polynomial-time safe reduction (Rule 6.2 on whole components of I(H), any size, decided by a CSP on a
   tree decomposition of I(H)) after which at most l + l' <= 2l components of I(H) remain, whatever their size. In particular l = 0 is fully reduced
   (consistent with Thm A) and ML is reduced to instances with at most 2l components of I(H).
3. New PROVED reformulation (Prop 1): tw(Gamma') is controlled, up to a factor Z*, by the treewidth of the soft graph W with its far-apart vertices
   identified along the (small) hard components. All the unbounded treewidth of Gamma' comes from identification, none from W or from G-edges as such.
4. The remaining gap is precisely the following "trigger" scenario (Sec. 6): ONE tight component whose irreducibility is caused by <= l outside
   vertices whose hard pairs must be edited, with the edited vertex at an unknown position. Branching on the trigger works in the family of Prop 5
   (bounded candidate set, l drops), but the forced-violator form of this branching is false in general (20% of tight irreducible components in the
   random tests have no forced violator, Sec. 6), so it is not a proof.

Nothing below changes the status "FPT in l + w is open". Under c <= 2 with I(H) a forest (w <= 1) the case is still open.

## 1. Conventions
(G,H,l), H a packing of modification-disjoint induced P3s, pair-maximal unless said otherwise; soft = packed pair, hard = other; W = soft graph;
HG = graph of hard edges; hard component = component of HG; D = 2c; B = 2D+3l, Z = B+D+1; BIG = hard components with a vertex of hard degree > B;
Gamma' = (V minus BIG, G-edges and soft pairs inside). A solution S is a set of pairs with G xor S a cluster graph; U = S minus packed pairs (<= l).
Components T_1..T_r of I(H), W_j = union of the vertex sets of the members of T_j, p_j = |T_j|, tau_j = OPT(G[W_j]) >= p_j.
T_j is tight if tau_j = p_j. For a set F of members, W_F = union of V(P), P in F; F is pair-closed if no member outside F has two vertices in W_F
(then no outside member has a pair inside W_F, and the instance stays a valid CEaMP instance after editing pairs inside W_F).
Rule 6.2 (BFK18, specialised): F pair-closed, S an optimal solution of G[W_F] such that for every y outside W_F the set N(y) cap W_F is empty or
exactly one cluster of G[W_F] xor S ("R-condition"; this is (R1)+(R2) of components.tex). Then OPT(G) = OPT(G xor S) + |S| (BFK Lemma 6.3, CITED,
sound for every induced subgraph). If |S| = |F| (F tight) the instance (G xor S, H minus F, l) is equivalent to (G,H,l) (members outside F have no pair inside W_F).
F is *reducible* if some such S exists, *irreducible* otherwise.

## 2. Reformulation of ML  (Prop 1)
**Prop 1 (PROVED, using Lemmas 3,5 of wave5_dp_route.md; VERIFIED-BY-COMPUTER below).** Let (G,H,l) be a YES instance, H pair-maximal.
(a) Gamma' = W[V'] union HG[V'] (every G-edge is hard or soft; soft pairs are in W).
(b) Put Z* = max{Z, (l+1)(2D+2l+1)}. Every hard component K that is not in BIG has |K| <= Z*. Consequently "reject if a non-BIG hard component is larger
    than Z*" is a safe polynomial rule.
(c) Let Q' = Gamma'/(hard components). Then tw(Q') <= tw(Gamma') <= Z*(tw(Q')+1) - 1.
Proof. (a) trivial. (b) Take an optimal S, |S| <= |H|+l, U its unpacked edits. If no edge of K is deleted, K lies in one cluster C. C exceptional:
|C| <= 2D+2l+1 <= Z (Lemma 3(b)). C non-exceptional: C is a union of whole hard components (Lemma 5 beta); several components give |C| <= 2D; a single
component K = C is clean, so every y in K has all other vertices of K as hard neighbours or soft partners and |K| <= d_h(y)+D+1 <= B+D+1 (K not BIG).
If edges of K are deleted, K minus the deleted edges has at most l+1 pieces (components of kept hard edges); K is connected so each piece contains an
endpoint of a deleted edge, hence lies in an exceptional cluster of size <= 2D+2l+1: |K| <= (l+1)(2D+2l+1). (c) hard components are connected in Gamma'
(by hard edges), so Q' is a minor of Gamma'; conversely replace each node of a tree decomposition of Q' by the vertices of its component. QED
Computer (w6a_comp.py, 1,565 instances, n = 6..10, 3,638 optimal partitions): every non-exceptional cluster is a union of whole clean components with all
cross pairs soft; every hard component split by a cluster has all its pieces in exceptional clusters and |K| <= (|U_K|+1)(1+2D+|U|); clean single-component
clusters satisfy |K| <= d_h(y)+D_y+1; exact treewidth sandwich tw(Q') <= tw(Gamma') <= max|K|(tw(Q')+1)-1 on 1,211 instances. 0 failures.
(57+58 hard components were unclean: pair-maximality does not make components clean; they all lie in exceptional clusters.)
Reading: tw(W) <= 3w+2 (packing.tex Prop soft). Gamma' has large treewidth only because far-apart vertices of W are identified (put into one small
hard component). ML says: a poly-time safe reduction can undo this identification (Prop 9 of wave5 is the model: every K_v is "inert").

## 3. Rule 6.2 on whole components: the new partial reduction

**Lemma D (violation => an edited pair at the boundary). PROVED.** Let S be any solution of G, F pair-closed, S_F = S cap C(W_F,2) and y outside W_F.
If N(y) cap W_F is not empty and not a cluster of G[W_F] xor S_F, then S contains a pair {y,u} with u in W_F.
Proof. Otherwise N_{G xor S}(y) cap W_F = N_G(y) cap W_F. In the cluster graph G xor S this set is cluster(y) cap W_F, a cluster of
(G xor S)[W_F] = G[W_F] xor S_F, or empty. QED
(This is the engine of the BFK counting, K3 in components.tex. The pair {y,u} is a packed pair of a member through u (a neighbour of F in I(H), edit "free")
or a hard pair, i.e. in U. For c = 1 and F a whole component only the second case exists, which is why BFK counts.)

**Prop 2 (polynomial test for tight components). PROVED.** For a component T of I(H), with a tree decomposition of I(H) of width w, one decides in time
3^{w+1} n^{O(1)} whether T is tight and whether it is reducible (and computes S).
Proof. Tight S = one edited pair per member, nothing else, with G[W_T] xor S a cluster graph; such an S is optimal (tau_T >= p_T). Variables: the 3 choices of
each member. Constraints: all triples of W_T containing >= 1 packed pair; the (<= 3) members owning its pairs pairwise share a vertex, so they form a clique
of I(H) and lie in a common bag; triples with no packed pair are checked statically. The R-condition for y outside W_T: A_y = N(y) cap W_T must be a cluster:
(i) every pair uv inside A_y is an edge of G xor S, (ii) every pair u in A_y, v in W_T minus A_y is a non-edge of G xor S. Each such pair is either hard
(status fixed: static test) or owned by one member (a unary constraint on its choice). So reducibility is a CSP of treewidth w. QED
(This is Thm A at l = 0 with unary constraints added. w6a_csp.py contains a DFS solver, checked against brute force on 2,861 subfamilies.)

**Thm 3 (component reduction). PROVED modulo BFK Lemma 6.3 (CITED).** Let (G,H,l) be an instance with H any packing. Let l' = l - sum_j (tau_j - p_j)
(NO if l' < 0; tau_j is only needed for non-tight components, see below). Apply Rule 6.2 to every reducible tight component. In the resulting instance,
with A = tight irreducible components and B = non-tight components, |A| <= 2l', |B| <= l - l', so the number of components r satisfies r <= l + l' <= 2l.
Proof. Rule: Lemma K2 of components.tex (S_j consists of pairs inside W_j, so R-conditions and neighbourhoods of the other components do not change).
|B| <= sum (tau_j - p_j) = l - l'. Counting for A (K3 with the non-tight components left in): let S be optimal, |S| <= |H|+l = sum_j tau_j + l'. For j in A either
|S_j| >= tau_j+1 or S_j is an optimal solution of G[W_j]; then, j being irreducible, S_j violates the R-condition, and Lemma D gives an edited pair
{y,u}, u in W_j, y outside (a pair external to j; it is hard because every packed pair with an end in W_j lies inside W_j). An external pair is external to
<= 2 components and lies in no C(W_i,2). So |S| >= sum_j tau_j + |A_a| + |A_b|/2 >= sum tau_j + |A|/2, hence |A| <= 2l'. QED
Computing tau_j for non-tight components is not needed: they are simply left alone.
Alternation with MAX (packing.tex Lemma MAX): after the rule the packing may lose pair-maximality; re-maximalising costs |M| >= 1 of l and raises w by <= |M|.
Since l decreases, the alternation stops after <= l rounds, giving a pair-maximal instance with l* <= l, w* <= w+l and <= 2 l* components of I(H*). (PROVED.)

**Corollaries.** (a) l = 0: a YES instance is reduced to the empty packing, i.e. everything is decided (agrees with Thm A). (b) If every component has vertex span
<= s, ML holds (Theorem K; kernel of size O(s^2 l)). (c) ML is equivalent to ML for instances with at most 2l components of I(H) (what remains is the
case of FEW, ARBITRARILY LARGE components, e.g. one tree of members when w = 1, Open Problem forest).
Computer (w6a_rule.py seeds 11-14, 5,921 pair-maximal instances, n = 7..11; w6a_k3.py 3 seeds x 1,200 instances x budgets l = OPT-|H|+{0,1,2}):
soundness OPT(G xor S_F) = OPT(G) - |F| for EVERY R-satisfying tight S_F of EVERY connected pair-closed F with |F| <= 3: 7,076 applications, 0 failures;
whole-component reduction of all reducible tight components in one pass: OPT(reduced) = OPT - sum p_j in all 3,552 instances,
|A| <= 2l' and r <= l + l' in all 10,656 (instance, budget) pairs; all 1,929 instances with OPT = |H| were reduced completely (every component reducible); max r after reduction = 3.
Evidence only; the proof is above.

## 4. Why Rule 6.2 on unions of members cannot prove ML (obstruction)

### 4.1 Bounded-size subfamilies are not enough even at l = 0 (VERIFIED-BY-COMPUTER, plus Prop 4 below)
Components that are reducible only as a whole exist (w6a_minred.py, seed 5: components with |T| = 2..11 of which many have s(T) = |T|, s(T) = smallest reducible
connected pair-closed subfamily; w6a_rule.py: 447 components with >= 3 members that are reducible but have no reducible subfamily of <= 2 members).
Sparse paths of members (I(H) a path, w = 1, l = 0, pair-maximal, OPT = |H|, no reducible proper connected subfamily, whole component reducible):
exhaustive search over all graphs with <= 4 (m = 3), <= 5 (m = 4) extra edges (w6a_exh.py) and extensions (w6a_ext.py, m = 5, <= 5 extra edges, 7 of 40
candidate patterns hit), re-verified by independent partition-based code (w6a_indep.py):
* m=3, n=7: edges 01 05 12 13 34 35 36 46, H=(0,1,2),(1,3,4),(5,3,6)
* m=4, n=9: edges 01 07 12 13 17 34 35 56 57 58 68, H=(0,1,2),(1,3,4),(3,5,6),(7,5,8)
* m=5, n=11: edges 01 05 08 12 13 19 29 34 35 36 46 57 58 79 7-10, H=(0,1,2),(1,3,4),(5,3,6),(7,5,8),(9,7,10)
No claim is made that such chains exist for every m on w = 1 at l = 0 (they have hubs; the family below is the scalable one). At w = 2 an infinite family exists
(w6a_fam.py: "clique necklace with expander leaves" N(F), l = 0; triangle, K4, K3,3: all connected subfamilies of size <= 3 (<= 2 for K3,3) irreducible, whole reducible).

### 4.2 The family N*(F): w = 1, l = 1, every union of members irreducible, tw(Gamma') unbounded  (Prop 5)
Construction. F a cubic graph with an orientation of in-degree <= 2 at every vertex, F' = F plus one new vertex * and one edge e_0 = (t -> *). Order the edges
e_0, e_1, ..., e_N of F' so that e_0 comes first and consecutive edges have different heads (possible: each head occurs <= 2 times). For every edge e_i take
vertices a_i,b_i,c_i,u_i, plus z and y. For v in V(F') let K_v^+ = {a_i : tail(e_i)=v} cup {b_i,c_i,u_i : head(e_i)=v} (so K_*^+ = {b_0,c_0,u_0}).
G: all pairs inside each K_v^+; cross edges a_i b_i; spine edges b_i u_{i+1} (i < N), b_N z; trigger y adjacent to u_0 and c_0 only.
H: P_i = (a_i,b_i,c_i) (centre b_i) and X_i = (u_i,b_i,u_{i+1}) (centre b_i; u_{N+1} := z).
**Prop 5 (PROVED, VERIFIED-BY-COMPUTER for F = triangle, K4, K3,3, cube, Petersen; n up to 66).**
(F1) The members are induced P3s with pairwise disjoint pair sets, H is pair-maximal, I(H) is a caterpillar (spine X_0..X_N, leaf P_i at X_i), so w = 1, c = 2, D = 4.
(F2) The partition {K_v^+} + {z} with y put into K_*^+ costs |H|+1 (deleted: the 2(N+1) cross and spine edges, all packed; inserted: yb_0). So (G,H,1) is YES.
     OPT = |H|+1 exactly (OPT >= |H|+1 follows from (F6) with Thm 3 at l = 0; also computed exactly for F = triangle: OPT = 9 = |H|+1).
(F3) Hard degrees <= 8 <= B = 11: BIG = empty. All hard components are small (|K_v^+| <= 9); the hard component of y is unclean.
(F4) Contracting every K_v^+ (connected by clique edges) and using the cross edges a_i b_i gives F as a minor of Gamma', so tw(Gamma') >= tw(F) = Omega(|V(F)|)
     for cubic expanders (Grohe-Marx, CITED). So tw(Gamma') is unbounded although l = 1, w = 1, BIG = empty.
(F5) I(H) has a single component T and it is tight: the partition of (F2), restricted to W_T, edits exactly one packed pair per member (y is outside W_T), so tau_T = |H|.
(F6) Every nonempty pair-closed set F of members (connected or not) is Rule-6.2-irreducible; in particular T is irreducible and so is every single member.
     (tau(G[W_F]) = |F| because the solution of F2 restricted to W_F has |F| edits, so only tight S_F matter.)
Proof of (F1): P3 of G need a cross edge (cliques contain none and every vertex has few cross edges), and every cross/spine edge is a soft pair; the trigger edges
yu_0, yc_0 lie in no P3 without a soft pair because y sees exactly {u_0,c_0} and (u_0,b_0), (c_0,b_0) are soft. Consecutive spine members have different heads so the
non-edge u_i u_{i+1} is real; pair-disjointness is a direct check. (F2)-(F5) direct.
Proof of (F6). Let S_F be tight (one edited pair per member of F, G[W_F] xor S_F a cluster graph).
Case 1: F contains a spine member; let X_i be the one of smallest index.
 * i >= 1. Then b_{i-1} is not in W_F (X_{i-1} not in F by minimality; P_{i-1} in F would make X_{i-1} have the pair (b_{i-1},u_i) inside W_F, against pair-closedness),
   and a_{i-1}, u_{i-1} are not in W_F either. So N(b_{i-1}) cap W_F = {u_i} cup M, M = K_{h(i-1)}^+ cap W_F. If this were a cluster C of S_F containing u_i, every m in M would
   be adjacent to u_i in G xor S_F; but (u_i,m) is a hard non-edge unless m = u_{i+1} (soft pair of X_i); then X_i inserts (u_i,u_{i+1}), so b_i is in C but not in
   {u_i} cup M: contradiction. Hence M is empty and C = {u_i}. Then X_i deletes (u_i,b_i). If P_i is not in F, c_i is outside W_F and adjacent to u_i and b_i, which are now in
   different clusters: violation. If P_i is in F, the hard edge u_i c_i is kept and C != {u_i}: contradiction. So F is irreducible.
 * i = 0. The trigger y is outside W_F and N(y) cap W_F = {u_0} (P_0 not in F) or {u_0,c_0} (P_0 in F). First case: as above C(u_0) = {u_0}, X_0 deletes (u_0,b_0),
   and c_0 (outside, adjacent to u_0 and b_0) violates. Second case: the cluster must be exactly {u_0,c_0}, so b_0 is separated from both: X_0 deletes (u_0,b_0) and P_0 deletes
   (b_0,c_0); then a_0 and u_1 (both in W_F) are adjacent to b_0 and non-adjacent to each other (hard non-edge, unedited): an induced P3 in G[W_F] xor S_F, contradiction.
Case 2: F contains only leaves. Take P_i in F. The vertices u_i, u_{i+1} are outside W_F. u_{i+1} sees b_i and the vertices of K_{h(i+1)}^+ cap W_F; a vertex m of the latter is a hard
 non-neighbour of b_i, so unless that set is empty the trace is not a cluster; if it is empty the trace is {b_i}, which is a cluster only if b_i is a singleton in S_F; but each of the three
 edits of P_i keeps b_i together with a_i or c_i, or merges all three. Violation. QED
Computer: w6a_trig3.py (F2, I(H) a forest, pair-maximal, T tight and whole-irreducible, no reducible connected subfamily with |F| <= 3-4; exact OPT = |H|+1 for the triangle),
w6a_alls.py (ALL pair-closed subfamilies, connected or not: triangle: 95 of size <= 8, K4: 1,030 of size <= 5; 0 reducible).

**Consequence (precise obstruction).** No reduction obtained by iterating Rule 6.2 on unions of members (the form used for c = 1 and in Theorem K) can prove ML:
on N*(F) it does nothing at all, l = 1, w = 1, and tw(Gamma') = Omega(|V(F)|). In particular the counting argument of BFK (Lemma K3) cannot be extended from components to members:
all 2(N+1) members of N*(F) are irreducible although l = 1. The mechanism (Lemma D): a member through u is irreducible whenever a neighbour in I(H) edits a pair at the shared
vertex and the cluster of u meets the member in a second vertex; along a path of members this gives a chain of dependencies that ends at the trigger, and the chain has unbounded length.
Alternative generalisations (Rule 6.2 for W that is not a union of members, i.e. W containing extra vertices) are NOT excluded: a W that contains the trigger and everything it
depends on works, but then it is the whole component and solving G[W] is the problem itself.

## 5. Where the structure of the obstruction points (not proved)
* In N*(F) the trigger y is the unique vertex violating every tight S_T ("forced violator"). The unpacked edits at y have a bounded candidate set (edit (y,u_0), (y,c_0), (y,b_0)
  in this example; in general Cor. 6 of wave5 bounds the partner by g(l,D)). Branching on them lowers l by one; at l = 0 Cor. (a) of Thm 3 finishes. So N*(F) is solved in
  3 * poly time by "reduce, branch on a forced violator, reduce". This is the shape a proof of ML could take (with branching replaced by a safe reduction).
* Dependency structure (PROVED by Lemma D). Let S be an optimal solution whose restriction to a member Q is tight and which has no unpacked pair at V(Q). If S_Q satisfies the R-condition,
  Q is reducible. Hence an irreducible Q has an I(H)-neighbour P' whose edited pair contains a vertex shared with Q (arc Q -> P'); at l = 0 this holds for every irreducible member
  and every optimal S. Consider a fixed point of the single-member rule at l = 0 (every member irreducible): every member has an out-arc to a member, arcs run along edges of I(H),
  in-degree <= 2(c-1) (P' has one edited pair, containing <= 2 of its vertices). So every component of I(H) contains a mutual pair P <-> P' or a cycle (functional graph), which for a forest
  means a mutual pair. For l > 0 the arcs may also end at "sinks": members with a non-tight restriction or with an unpacked edit at one of their vertices. An unbounded chain of irreducible members
  can therefore only hang off a sink or a cycle, which is what N*(F) does with the trigger (sink = X_0 with the unpacked edit at y).

## 6. The forced-violator branching is incomplete (VERIFIED-BY-COMPUTER, CONJECTURE for what remains)
For a tight irreducible component T define y (outside W_T) to be *forced* if every tight S_T violates the R-condition at y. If some y is forced, every optimal solution whose restriction to
T is tight contains a pair {y,u} (Lemma D), a bounded set of candidates (Cor. 6 of wave5), so one branch lowers l.
w6a_viol.py (6,000 random instances, n = 7..11): 592 tight irreducible components, 472 with a forced violator, **120 without** (every y individually satisfiable by some S_T, no S_T satisfies all).
w6a_core.py: in all those tests the smallest jointly unsatisfiable set of outside vertices has size 1 or 2 (distribution: size 1: 472, size 2: 120 over components with |T| <= 7; the 472 are exactly the components with a forced violator).
Also unresolved: solutions whose restriction to T is NOT tight (excess spent inside W_T); their positions are the n^l part of Thm 7 and the Thm 3 counting does not bound them.
**Conjecture TL (trigger lemma, not tested at scale).** If T is a tight irreducible component, there is a set Y of at most g(l,D) outside vertices, computable in polynomial time,
such that every solution of cost <= |H|+l whose restriction to W_T is tight edits an unpacked pair at a vertex of Y. TL + a bound for non-tight restrictions would turn Thm 3 + MAX into a
branching FPT algorithm in (l, w); TL for jointly unsatisfiable cores of size <= 2 is consistent with all random tests, but a long dependency chain could, in principle, produce larger cores.

## 7. Strongest statements (summary table)

| statement | status |
|---|---|
| Gamma' = W u HG; non-BIG hard components have <= Z* vertices; tw(Q') <= tw(Gamma') <= Z*(tw(Q')+1)-1 (Prop 1) | PROVED; VERIFIED-BY-COMPUTER (n <= 10) |
| Lemma D (violation => edited boundary pair) | PROVED |
| tight/reducible components decidable in 3^{w+1} poly (Prop 2) | PROVED |
| Component reduction leaves <= l+l' <= 2l components (Thm 3), alternation with MAX gives pair-maximal, w* <= w+l | PROVED modulo BFK Lemma 6.3 (CITED); VERIFIED-BY-COMPUTER |
| l = 0 reduces to empty (consistent with Thm A) | PROVED (Cor. a) |
| ML when all components have bounded span (Thm K) | PROVED modulo BFK (already in paper) |
| N*(F): w = 1, l = 1, BIG = empty, tw(Gamma') >= tw(F), every union of members Rule-6.2-irreducible (Prop 5) | PROVED; VERIFIED-BY-COMPUTER for 5 graphs F |
| Whole-only sparse chains, w = 1, l = 0, m = 3,4,5 | VERIFIED-BY-COMPUTER (examples only) |
| Forced-violator branching complete | FALSE as stated (120/592 components have none) |
| Trigger Lemma TL | CONJECTURE |
| ML (poly-time reduction to tw(Gamma') <= f(l,w)) | CONJECTURE; open also for w = 1 |

Precise remaining obstruction: after Thm 3 the instance has <= 2l components; a single component with an unbounded number of members can have an expander in Gamma'
(N*(F)); the only mechanism that keeps it irreducible is a bounded number (<= l) of outside "trigger" vertices whose hard pairs have to be edited; a proof of ML (or an FPT algorithm)
must locate the triggers in polynomial time (TL) or avoid locating them. Members not on a dependency chain are removed by Thm 3 / Rule 6.2; chains of dependent members hanging off a trigger
are the only thing Rule 6.2 cannot remove.

## 8. Test log
| script | what | scale | result |
|---|---|---|---|
| w6a_rule.py | soundness of Rule 6.2 on connected pair-closed F, |F| <= 3; K3-type counts; whole-only components | seeds 1, 11-14; 5,921 instances n=7..11 (+146 n=6..9) | 7,076 + 207 applications, 0 soundness failures, 0 K3 failures, 447 whole-only components |
| w6a_k3.py | whole-component reduction, OPT drop, |A| <= 2l', r <= l+l', l = 0 fully reduced | 3 x 1,200 instances, 3 budgets | 0 failures; 1,929 OPT=|H| instances fully reduced |
| w6a_comp.py | Prop 1 | 1,565 instances, 3,638 optimal partitions, n <= 10 | 0 failures |
| w6a_minred.py, w6a_chain*.py, w6a_exh.py, w6a_ext.py, w6a_indep.py | whole-only chains | exhaustive (m=3,4), partial (m=5), n <= 11 | examples above, independent recheck OK |
| w6a_fam.py | N(F) (w=2, l=0) | n up to 36 | all properties; no small reducible subfamily; whole reducible |
| w6a_trig.py (random C4 triggers) | negative result | 3,000 draws | none passed the filters (a C4 trigger does not survive pair-maximality / irreducibility on N(F)); the y-trigger of Prop 5 was found by analysis |
| w6a_trig2.py, w6a_trig3.py, w6a_alls.py | N*(F) | n up to 66; all pair-closed F for triangle (95), K4 (1,030) | Prop 5 properties hold; OPT = |H|+1 exact for triangle |
| w6a_viol.py, w6a_core.py | forced violators, core sizes | 6,000 instances n=7..11 | 120/592 without forced violator; cores <= 2 |
| w6a_csp.py | DFS solver vs brute force | 2,861 subfamilies | identical |
Limits: all brute-force data have n <= 11 (OPT by subset DP, exact), N*(F) up to n = 66 uses explicit solutions and static checks (no OPT beyond n = 18);
random instances are small, so whole-only and no-forced-violator phenomena are lower bounds on what happens at scale; w6a_opt.c needed for w6a_trig*.py.
