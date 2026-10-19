# Audit 7 (hostile referee, wave 6): notes/wave6_pointers.md and notes/wave6_ml.md

Context read: paper/sections/packing.tex, locality.tex, dproute.tex, components.tex, notes/audit5_locality.md (not re-opened beyond what the sections state), audit6_dp.md.
Own code, written from scratch, imports only numpy/networkx and `a7_lib.py`: `code/a7_lib.py`, `a7_w1.py`, `a7_w2tight.py`, `a7_lb.py`, `a7_thmP.py`, `a7_thmP_ctl.py`,
`a7_thmP_fam.py`, `a7_far.py`, `a7_ml.py`, `a7_prop1.py`, `a7_n5.py`, `a7_n5b.py`. No git, paper/ untouched.

## Verdicts

| claim | verdict | one-line reason |
|---|---|---|
| W1 conflict lemma | CORRECT | proof is two lines and right; the corollary M(z,b) <= l needs z,b merged (scope remark) |
| W2 defect bound sum c(pi) <= 2D|U| | CORRECT | proof right; "tight" only evidenced for D = 2 |
| Cor. LB (sum c_i + ceil(W/2D)) | CORRECT WITH FIX | valid lower bound; status table says "exact at l = 0": FALSE (29/300 instances) |
| F1 pruning | CORRECT | direct from W1 and |U| <= l |
| DEL rigid edges | CORRECT | |
| BAD (<= 2c eps bad members) | CORRECT | |
| Theorem P (FPT in (l,D,tw(Gamma^l))) | CORRECT | pruning is sound, DP is the dproute DP with e^+; own DP matches brute force; control with threshold l-1 fails |
| "strictly stronger than omega theorem" | CORRECT WITH FIX | holds for the K4/Petersen/cubic family only for orientations with every in-degree >= 1; for an arbitrary orientation min weight is 4 |
| Prop. FAR | CORRECT | proof checked line by line, verified m = 3..7 |
| refutation of NF-bad | CORRECT WITH FIX | genuine counterexamples, but only for "far" = I(H)-distance 2 and 3; unbounded distance is NOT refuted (FAR's far ends are the bad members) |
| (candidate lists size <= (D+l)(D+1), Sec. 0 item 3) | CORRECT WITH FIX | labelled PROVED, no proof in the note; holds for T2/deleted edits of an optimum only (T1 excluded), derivation below |
| ML Prop 1 (a)(b)(c) | CORRECT | |
| ML Lemma D | CORRECT | |
| ML Prop 2 / 3^{w+1} CSP | CORRECT | |
| ML Thm 3 (component reduction) | CORRECT (modulo cited BFK Lemma 6.3, which my tests support) | two non-algorithmic phrases (tau_j of non-tight components) are harmless |
| ML Prop 5 family N*(F), w=1, c=2, l=1, tw(Gamma') >= tw(F) | CORRECT WITH FIX | claims hold in all tests, but the proof of (F6) has two unproved sub-steps (ordering condition; Case 2) |

Nothing is WRONG. No brute-force counterexample to any proved statement.

## wave6_pointers.md

### W1 -- CORRECT
(a) merged z,b: x(z,y) = x(b,y) in G xor S; hard pairs zy, by with different G-status: exactly one is in S, hence in U. (b) unmerged, y a common hard neighbour: if neither
edge is deleted then z~y~b are one cluster. Distinct y give distinct pairs. Both parts re-derived. Brute force (`a7_w1.py`): about 342,000 partitions of cost <= |H|+3 over random pair-maximal
(and 20% non-maximal) instances n = 5..9, about 3.4e5 conflicts, 0 violations. W1 needs neither maximality nor optimality.
Scope remark: the Corollary "M(z,b) subset Delta(z,b), so |M(z,b)| <= l" is true only when z,b are in one cluster of the solution (for unmerged z,b the set M can be large: hard neighbours of z not adjacent to b
are not constrained). locality.tex Lemma near states the merged hypothesis in its first line; the note's sentence should say it too.

### W2 -- CORRECT; tightness only for D = 2
Proof: the map conflict -> edit in {zy,by} is a function; an edit e = (p,y) is the image of conflicts with pi containing p (<= |Q(p)| <= D) and conflicts (pi',p) with pi' containing y (<= D). Right.
0 violations. Tightness: `a7_w1.py`/`a7_w2tight.py` reproduce ratio exactly 1.000 (n = 6: 65 partitions in 6,000 instances, n = 7: 4), but ONLY for D = 2 (c = 1); for D = 4 the maximum I see is 0.5, for D = 6 it is 0.33 (n <= 8;
my generator is weaker than theirs, so this is weak evidence). So "the constant 2D is attained" is supported for D = 2 only; for D >= 4 attainment is open (the proof gives no reason it should be tight).
Corollary LB (cost >= sum c_i + ceil(W/2D)): follows from CP(b); 0 violations. The DP claim is plausible (counter capped at 2Dl+1). FIX: the status table says "exact at l = 0", the text of Sec. 1 says "I do not claim it is exact at l = 0".
The table is wrong: `a7_lb.py` computes the relaxation exactly (min over soft-state vectors transitive on all-soft triples of sum c_i + ceil(W/2D)): value <= OPT in all 300 instances, but in 29 of 300
instances value = |H| while OPT > |H| (the triples with two soft pairs and one hard pair are not checked), see `python3 a7_lb.py 1 300 6,7,8` (prints one instance: n = 7, H = (2,0,3),(1,5,3),(2,4,6), relaxation value 3 = |H|, OPT = 4).

### F1, DEL, BAD -- CORRECT
F1: W1 + |U| <= |S|-|H| <= l; about 5.5e6 pair-checks (l = eps..eps+2), 0 violations. DEL: |U| >= 1 + |common hard nbrs| (distinct w give distinct pairs, all different from zy): 0 violations in ~2.6e5 deleted edges.
BAD: eps = sum(c_i - 1) + |U|; bad member with c_i >= 2 counted in the first sum, bad member with c_i = 1 sits in a cost-1 state with a conflict, W1 gives an edit of U incident to V(P), each edit has two ends, each end in <= c members:
#bad <= sum(c_i-1) + 2c|U| <= 2c eps. Right. 0 violations on 3.4e5 partitions (bad count, not only the bound, was computed from the definition of "three states").
Candidate-list claim (Sec. 0 item 3, "<= (D+l)(D+1)"): no proof in the note. It is right for T2 and deleted edits of an OPTIMUM: from the hard end a, z in N_h(a) (<= D+l by Lemma EX(a)), b in Q(z) (<= D), plus deleted partners N_h(a);
from the soft end b, z in Q(b) (<= D), a in Delta(z,b) (<= l by W1). T1 pairs (z soft to both) are not in this count (they need D(D-2) more) and are bag-local anyway. State this and the optimality hypothesis (T1/T2 uses Lemma half).

### Theorem P -- CORRECT
Proof checked step by step.
(<=) f <= f~ since e^+ <= e, so any partition from the DP has true cost <= |E| + sum f(K) + m^*; no structure needed.
(=>) S optimal, |S| <= |H|+l. Every soft pair inside a cluster is merged, so |Delta| <= l by F1: it is an edge of Gamma^l; every G-edge inside a cluster of V' is hard (edge of Gamma^l) or a merged soft edge (edge of Gamma^l).
So the clusters are Gamma^l-connected (NF supplies connectivity through exactly these pairs) and f~(C) = f(C). The DP of dproute.tex only ever uses pairs lying in a bag; replacing "G-edge" by "Gamma^l-edge" everywhere (introduce: k = Gamma^l-neighbours in the block; join:
no Gamma^l-edge between exclusive parts; g(b) with e^+) gives exactly sum f~. Gamma^l is a spanning subgraph of Gamma' on V', so tw(Gamma^l) <= tw(Gamma').
Own implementation (`a7_thmP.py`): Gamma^l, BIG, Z from scratch; m_all and m_conn by brute force over all partitions of V' (cluster size <= Z), and an own tree-decomposition DP (introduce/forget/join with block sizes, min-fill decomposition, per component).
Random pair-maximal instances n = 6..10, l = 0..3: 3,340 (instance, l) pairs (about 1,100 prunes of soft pairs per 1,000 pairs; BIG non-empty in 4 instances at n = 8..9, l = 0..1):
the criterion |E| + sum f(K) + m <= |H|+l <=> OPT <= |H|+l holds for m_all, m_conn and m_dp in all of them; m_all <= m_dp <= m_conn in all of them. Sensitivity control (`a7_thmP_ctl.py`, prune with threshold l-1 at budget l): 29 of 900 pairs
give a WRONG answer, so the test does detect unsound pruning. (Z >= 3D+1 > n for most tested instances, so the size bound is not stressed; it is the same argument as in dproute, audited in audit6.)
"Strictly stronger": `a7_thmP_fam.py` builds the dproute expander family (K4, Petersen, random cubic on 20 and 40 nodes, orientation with all in-degrees in {1,2}): D = 2, BIG empty, cross soft weights {5,6,7}
(general formula |Delta| = 3 + in(p) + in(q) for the pair on the edge p -> q), for l <= 4 every cross pair is pruned, every component of Gamma^l has <= 5 vertices (tw <= 4), while the quotient of Gamma' is F
(tw(Gamma') >= tw(F)); at l = 5 the components merge (sizes 8, 20, 32, 40). FIX: the statement is true for orientations in which every vertex has in-degree >= 1; with an arbitrary orientation a source tail
gives weight 4. The generalisation to "d-regular expanders with d ~ l/2" is a sketch, as the note says.

### Proposition FAR -- CORRECT
Proof re-derived: (i) chord ends in members 0 and m-1, m >= 3 makes the cycle 1..2m-1 of length >= 5, so every member is an induced P3 and the only hard edge is the chord; every induced P3 has two edges, one of which is a soft path edge: pair-maximal.
(ii) unicyclic counting: an induced subgraph on s vertices has <= s-1 edges unless it contains the cycle (then <= s, s >= 5), so f >= -1 with equality iff s in {2,3} and connected; at most m such clusters; cost >= (2m+1) - m = m+1; the pair partition attains it.
(iii) the case analysis on the cluster C* containing both chord ends is complete (size 2: both pendants outside; size 3 with a pendant: remaining 2m-3 vertices cannot form m-1 clusters of size >= 2; f(C*) >= 0: 2m-1 vertices for m clusters).
Computer (`a7_far.py`, all set partitions): m = 3, 4, 5: OPT = m+1, 7 / 9 / 11 optima, every optimum separates the chord ends, the chord is a hard edit at I(H)-distance m-1 and two-sided, I(H) = path, locally bad members = {first, last}.
Extra: constrained subset DP for m = 6, 7: OPT = 7, 8 and best cost with the chord ends in one cluster = 8, 9, i.e. the forced deletion holds there too.

### Refutation of NF-bad -- CORRECT WITH FIX (scope)
Both examples reproduce (`a7_far.py`): n = 10, packing and pair-maximal, unique optimum, OPT = |H|+1 = 6, no locally bad member, exactly one edit (the hard edge {0,6}, resp. {0,8}) which is a deletion, two-sided (witnesses soft to both ends),
at I(H)-distance 2, resp. 3; in both I(H) has c = 3. So the literal conjecture "exists an optimum in which every far two-sided edit touches a locally bad member" is false for far = distance >= 2 and for far = distance >= 3.
FIX: neither instance bears on the relevant version. Pointers at distance <= K are bag-local in the decomposition of I(H)^K (width (w+1)(3c-2)^{K}-ish, a function of (w,c,K)); NF-bad matters only for far = distance larger than any function of (l,w).
Proposition FAR does not refute that version: its far ends are exactly the two locally bad members. So "NF-bad refuted, so bad members cannot be put in every bag" is established for K <= 3 only; for K -> infinity it is open (not proved either way).
The statement in Sec. 0 item 3 ("n = 10 instances show it also occurs at members that are not locally bad") is true as written but should add "at distance 2 or 3".

### Items not in the requested list, spot-checked
FE (kept hard edge az => x(a,t) = x(z,t)): trivial, right. AN weak form: right (flipping one pair of a cluster graph only creates P3's through that pair; all-hard triple with only e edited contradicts pair-maximality).
5.3 and Sec. 6 are evidence / ideas and are labelled so.

## wave6_ml.md

### Prop 1 -- CORRECT
(a) trivial. (b) For an optimal S with |S| <= |H|+l: no deleted edge in K: K lies in one cluster; exceptional: <= 1+2D+u <= Z; non-exceptional: union of whole hard components (NF beta), several components <= 2D, a single clean one <= d_h+D+1 <= Z.
Deleted edges in K: K minus <= l edges has <= l+1 pieces, each contains an end of a deleted edge (K connected), so lies in an exceptional cluster of size <= 2D+2l+1: |K| <= (l+1)(2D+2l+1). (Actually (l+1)(2D+l+1) suffices.)
(c) Q' is a minor of Gamma' (components are connected by hard edges in V'); the upper bound replaces each node of a decomposition of Q' by its component of size <= Z*. Right. `a7_prop1.py`: all optimal partitions of 200 instances (n = 6..8, 640 optimal partitions, 166 split hard components):
every piece of a split non-BIG component lies in an exceptional cluster, #pieces <= l+1, sizes within bounds, 0 failures (evidence only; n too small for Z* to bite).

### Lemma D -- CORRECT
Without a pair {y,u} in S, N_{G xor S}(y) cap W = N_G(y) cap W; in a cluster graph this is cluster(y) cap W, empty or a cluster of (G xor S)[W] = G[W] xor S_F. Done. `a7_ml.py`: 0 failures in about 4.0e5 (partition, y) checks, 1.3e5 of which had a violated R-condition.

### Prop 2 / the 3^{w+1} CSP -- CORRECT
Tight S = one edited pair per member, nothing else (c_i >= 1 forces equality). G[W] xor S is a cluster graph iff no induced P3 in a triple; triples with >= 1 soft pair involve <= 3 members which pairwise share a vertex (they own pairs through a common triple),
a clique of I(H), hence in one bag; triples with no soft pair are static. R-condition: every pair inside A_y, or between A_y and W minus A_y, is hard (static) or owned by one member (unary). Right. My test (`a7_ml.py`): tightness decided by enumerating 3^p states
agrees with tau_j = p_j computed by subset DP on all 1,350 instances (0 disagreements). Treewidth of the CSP primal graph is w (restriction of the decomposition to the component).

### Thm 3 -- CORRECT (modulo BFK Lemma 6.3, cited)
Counting proof checked: |S| = sum |S_j| + |S_ext| (disjoint pair sets), for tight irreducible j either |S_j| >= tau_j+1 or S_j is optimal for G[W_j], then irreducible gives a violated R-condition at some y and Lemma D an edited pair {y,u};
the pair is external to <= 2 components and hard (a packed pair with an end in W_j lies inside W_j). Non-tight components never enter the count. So |A| <= 2l', |B| <= sum(tau_j - p_j) = l - l', r <= l + l' <= 2l. Also correct: K2 (reducibility of other components
unchanged because S_F touches only pairs inside W_F), the alternation with MAX (each round with M != empty lowers l by |M| >= 1, width grows by <= |M|, final round has M = empty so the last Thm 3 pass is on a pair-maximal packing).
Remarks, not errors: (i) l' needs tau_j of non-tight components, which is not polynomial for a big component. The note says this, and the algorithm only needs the consequence "if more than 2l components remain: NO" (YES implies l' >= 0, l' <= l). Write the theorem that way.
(ii) the rule's equivalence OPT(G) = OPT(G xor S_F) + |S_F| is BFK Lemma 6.3. Own test (`a7_ml.py`, 1,350 instances n = 7..11, any packing: vertex-disjoint, pair-disjoint, non-maximal): every R-satisfying tight S_j of every component:
OPT(G xor S_j) + |S_j| = OPT(G) in 850 applications, 0 failures; one-pass reduction: OPT drops by exactly sum p_j in all instances, no remaining tight component is reducible afterwards, and for l in {OPT-|H|, +1, +2}
|A| <= 2l', |B| <= l-l', r <= l+l' in all 4,050 (instance, l) pairs (max |A| = 3 after reduction).

### Prop 5 (N*(F)) -- CORRECT WITH FIX
Own construction (`a7_n5.py`): packing, pair-maximal, I(H) a caterpillar (tw 1), c = 2, D = 4, hard degree <= 4 for the triangle (<= 7 in general) <= B = 11 so BIG empty; cost of the (F2) partition = |H|+1;
exact OPT by my own C subset DP (n = 18, triangle): OPT = 9 = |H|+1; T is tight and irreducible; all 90 pair-closed subfamilies with |F| <= 6 (triangle) irreducible; K4 (n = 30): pair-maximal, all 632 pair-closed subfamilies |F| <= 4 irreducible
for 8 random admissible orderings. tw(Gamma') >= tw(F): the quotient by the K_v^+ contains F (cross edges a_i b_i), confirmed structurally. OPT >= |H|+1 follows from (F6) for T and Lemma D counting at l' = 0, independent of BFK.
FIX (the proof of (F6), Case 1, i = 0, second sub-case, and Case 2):
1. "a_0 and u_1 ... a hard non-edge" is false when head(e_1) = tail(e_0) (a_0 in K_t^+ and u_1 in K_{head(e_1)}^+ = K_t^+ are adjacent). The note's ordering condition only asks consecutive heads to differ. Then {b_0, a_0, u_1} is a legal cluster and the written contradiction disappears.
The claim itself survives (triangle, orderings with head(e_1) = tail(e_0): T and all |F| <= 5 irreducible; K4: 4 such orderings, |F| <= 4 irreducible; a violator such as b_1 in K_t^+, adjacent to a_0 and u_1 but not to b_0, appears), but it is not proved. Cheapest repair: require head(e_1) != tail(e_0) in the ordering (always possible, t has in-degree <= 2).
2. Case 2 (only leaves, P_i in F): "m in K_{h(i+1)}^+ cap W_F is a hard non-neighbour of b_i" fails for m = a_i when tail(e_i) = head(e_{i+1}) (a_i b_i is a soft edge). Then u_{i+1} is a violator only for two of the three states of P_i;
for the state "delete b_i c_i" the violators are other vertices of K_{h(i+1)}^+ (b_{i+1}, c_{i+1}, ...): checked on the triangle with ordering e0,(1>2),(0>1),(2>0): P_1 state 'delete (5,6)' is violated by vertices 7, 9, 10, 12, not by u_{i+1} = 11.
Add this sub-case (b_{i+1} sees a_i but not b_i, so a_i must be a singleton, which forces the state "delete a_i b_i", and then u_{i+1} sees {a_i,b_i} which is not a cluster).
Other claims (F1)-(F5) verified. The consequence "Rule 6.2 on unions of members cannot prove ML" is as strong as the evidence for (F6) (|F| <= 6 triangle, <= 4 K4, whole T for the triangle) plus the proof above.

## Residual issues for the authors (ranked)
1. Status table of wave6_pointers lists Cor. LB as "exact at l = 0": false, delete.
2. NF-bad is refuted only at distance 2 and 3; present as such; no result for unbounded distance (FAR does not do it).
3. Prop 5 proof: the two sub-steps above; add condition head(e_1) != tail(e_0) and the Case 2 repair.
4. W2 tightness: state "attained for D = 2"; not shown for D >= 4.
5. "Candidate lists (D+l)(D+1)" labelled PROVED without proof; add the three-line derivation and the hypothesis (optimal S, T2/deleted).
6. Thm P strictness: all in-degrees >= 1; Thm 3: phrase the NO-rule through r > 2l instead of l' < 0.
