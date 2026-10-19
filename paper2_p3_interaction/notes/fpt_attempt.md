# FPT attempt for CEaMP in (l, w = tw(I(H)))  -- agent "fpt", 2026-10-19

Status in one line: I did NOT find an f(l,w) n^O(1) algorithm and did NOT find a hardness
proof. I proved an XP bound n^(2l) 7^(w+1), a clean "pointer" description of where the
difficulty sits, and several structural lemmas; the natural locality/wildcard routes
(a) and (d) fail as stated (concrete counterexample below). Everything proved was
tested by brute force (code/fpt_tests.py, code/fpt_tests2.py, code/fpt_examples.py).

Notation as in analysis_01. pairs(H) = union of the pairs of the P_i ("packed pairs"); any
other pair of vertices is "unpacked". Q(u) = vertices sharing some packed P3 with u.
A solution S (set of pairs to flip, G xor S cluster) with |S| <= |H|+l. R = V \ V(H).
Two vertices of one P_i form a pair of P_i, so: pair uz is packed  iff  z in Q(u).
The P3s through one vertex u pairwise intersect, so they form a clique of I(H); hence
|Q(u)| <= 2(w+1) and the clique lies in one bag (Helly).

## 1. Proved

### Lemma 0 (decomposition)  PROVED
For a solution S put x_i = S cap pairs(P_i) and U = S \ pairs(H). Then every x_i != {} and
|S| = sum |x_i| + |U|, so the excess is sum(|x_i|-1) + |U| <= l.
Proof. If S misses all three pairs of P_i, the triple {a,b,c} is still an induced P3 in G xor S.
The rest is bookkeeping because the pair sets are disjoint. QED
So a solution is (x, U): x_i in the 7 non-empty subsets of the 3 pairs of P_i, U a set of
unpacked pairs. Packed P3s edited twice or three times are therefore NOT the problem (see Thm A).

### Theorem A (XP in l, FPT in w)  PROVED, tested (T1)
CEaMP can be solved in time 7^(w+1) * n^(2l+O(1)).
Proof. Enumerate U among unpacked pairs, |U| <= l (at most n^(2l)). With U fixed, each triple
{u,v,z} of vertices is a P3 or not as a function of the status of its three pairs. Packed pairs
of the triple belong to <= 3 P_i's, any two of which share a vertex (Obs. 2), so the scope is a
clique of I(H); unpacked pairs are fixed (G-status, flipped iff in U). Thus we have a CSP on
variables x_i (domain 7, unary cost |x_i|) with constraints on cliques of I(H); bucket
elimination along a tree decomposition of width w costs 7^(w+1) n^O(1). Accept iff the min
total + |U| <= |H|+l. QED
(T1: this DP agrees with brute-force OPT <= |H|+l on 97 random instances, l = 0,1,2.)
This is no better than the known n^(2l) for c=2 (there tw(I(H)) is unbounded) but it handles
packed excess for free and shows that the whole difficulty is U.

### Theorem B (given x, the rest is FPT in l)  PROVED
Given a complete x (all packed pairs decided), deciding whether some U with |U| <= l-excess(x)
completes it takes 3^l n^O(1).
Proof. Standard search tree on G_x = G xor x. Take an induced P3 of the current graph; one of
its three pairs must be flipped; packed pairs are frozen (their status is part of x), so
branch on its (<= 3) unpacked pairs; budget drops by one. Correct since any completion must
flip an unpacked pair of every current induced P3 (no other change is allowed). QED
Consequence: CEaMP in (l,w) is "find x in a tree-structured space such that an FPT-checkable
property holds". The property is not local in x (it is where all the trouble comes from).

### Lemma 1 (WLOG H maximal, but this can raise w)  PROVED, tested (T3)
If H' = H + t further modification-disjoint induced P3s, then OPT <= |H|+l  iff  OPT <= |H'|+(l-t)
(both say OPT <= |H|+l). So for structure theory H may be assumed maximal: no induced P3 of G
has all three pairs unpacked. Warning: tw(I(H')) may exceed tw(I(H)); use only for structure,
not inside an algorithm unless the extension is chosen with bounded treewidth.

### Lemma 2 (viability of an unpacked edit)  PROVED, tested (T2: 2594 solutions)
Let S be a solution, U its unpacked part, uv in U. Let D(uv) be the set of z outside Q(u) u Q(v) u {u,v}
such that the triple {u,v,z} is an induced P3 in (G xor uv). Then for each z in D(uv) at least one of
uz, vz lies in U. Hence |D(uv)| <= |U|-1 <= l-1.
Proof. The three pairs of the triple are unpacked. If neither uz nor vz were flipped, the triple would be
an induced P3 in G xor S. Different z give different pairs. QED
Explicitly: deleting an edge uv needs < l common neighbours outside Q(u) u Q(v); inserting uv
needs |N(u) xor N(v) minus Q(u) u Q(v)| < l (near-twins). This is the familiar kernel-style condition,
and it does not bound the number of candidate pairs (all pendant edges are viable).

### Lemma 3 (an unpacked edit interacts only with the environments of its endpoints)  PROVED, tested (T4)
If uv is unpacked and a triple {u,v,z} contains a packed pair, then z in Q(u) u Q(v).
Proof. The packed pair is uz or vz. QED
Hence, for U fixed, the CSP differs from the base (U = empty) CSP only on triples containing a
pair of U; those with a packed pair live in the cliques of P3s through the endpoints of U-pairs.
A pair (u,v) in U is a "pointer" to two places (the bag nodes t(u), t(v)) of a tree decomposition.
All other triples (z outside Q) are tested on G and U alone ("C0 condition"), independent of x.

### Lemma 4 (structure of R for maximal H)  PROVED, tested (T6: 1048 pairs)
Let H be maximal and r, s in R. If rs is an edge then N[r] = N[s]; if rs is a non-edge then N(r) and
N(s) are disjoint.
Proof. Packed pairs lie inside V(H), so all pairs at r, s are unpacked. If rs is an edge and z != r,s is
adjacent to exactly one of them, {r,s,z} is an induced P3 with all pairs unpacked, contradicting
maximality; if rs is a non-edge and z is a common neighbour, {r,s,z} is an induced P3 likewise. QED
Corollaries (use also the standard fact F1: some optimal solution has all clusters connected in G,
because splitting a cluster that is disconnected in G saves |C1||C2| > 0 edits; therefore a clique
component of G is a cluster of every optimum):
 (i) R is partitioned into true-twin cliques X (component of G[R]); distinct classes have disjoint
     neighbourhoods; every neighbour of X outside X lies in V(H).
 (ii) Every vertex of R in a component meeting V(H) is adjacent to V(H) (shortest path argument:
      two adjacent R-vertices have the same neighbours, so a shortest path from r to V(H) has length 1).
 (iii) Symmetry: swapping two vertices of one class X is an automorphism of G fixing V(H) pointwise,
      hence maps solutions to solutions of equal cost. So, among the vertices of X, only the first
      min(|X|, 2l) (any fixed order) need be considered as endpoints of unpacked edits.
 (iv) (cited, not re-proved here) critical-clique property (Guo 2009): some optimum keeps a twin class in
      one cluster, so no deletion inside X is needed. [check the reference before using]
These shrink the candidate pairs to |V(H)|^2 + O(l * #classes * n) but remain polynomial of degree 2
in n for each of the l pairs, i.e. they do not give FPT.

## 2. Routes tried

### (a) "Excess edits lie in a bounded neighbourhood"  FAILS in the strong form; trivial in the weak form
Strong form "every unpacked edit has both endpoints in V(H)" is false. Concrete: G = C4 with edges
01, 02, 13, 23; H = {(0,1,3)} (induced path 0-1-3, pairs 01, 13, 03; vertex 2 is outside V(H)).
OPT = 2 = |H|+1 and the three optimal solutions are {03,12}, {01,23}, {02,13}; each contains an
unpacked pair at vertex 2 (code/fpt_examples.py). In the experiments (n <= 7) 9 of 400 random
maximal-H instances had every optimum touching R by an unpacked edit.
Weak form (edits within distance 1 of V(H)) is a theorem (Lemma 4(ii)) but says nothing: it still
leaves Theta(n^2) candidate pairs per edit. A bound on the number of candidates by g(l,w) would be
needed; I could not prove one and I did not construct a counterexample either (the expected
obstruction is a chain CSP of length m with an odd contradiction repairable at any of the m places).
 CONJECTURE C1 (supported by T5/T7: 0 forced cases in 142 instances with |R|>=2, n=7,8): some optimum
 has no unpacked edit with both endpoints in R. Not proved.

### (b) Branching on a violated P3  WORKS ONLY WHEN SCOPES ARE UNARY; fails in general
For w = 0 (I(H) edgeless; the vertex-disjoint case, FPT by van Bevern et al.) the following is a valid
FPT algorithm, I believe (sketch, not tested): with U fixed the x_i are independent, so if some P_i has no
feasible x_i, then for each of the 7 values m of x_i pick a triple violated by (m,U): a solution with
x_i = m must flip an unpacked pair of that triple, giving a 7*3 = 21-way branching of depth <= l.
For w >= 1 this breaks. Feasibility given U is a width-w CSP; its infeasibility witness (an unsatisfiable
core) can be a chain of length Theta(n): e.g. P_1..P_m a path in I(H), constraints forcing
x_1 = A, x_i = x_{i+1}, x_m = B. Any one of the m link-triples can be repaired by an edit, so a branching
set has size m, not f(l,w). So (b) gives n^l, not f(l,w).
(The unsatisfiable chain is realised in G by choosing the unpacked pairs inside each link; I did not
write the explicit graph, so this part is a sketch of the obstruction, not a verified instance.)

### (c) DP with "cluster identities"  Reformulated as a pointer CSP; no FPT proof
Lemma 3 says the true picture is: base CSP of treewidth w, plus <= l pointer variables
y = (u,v) of domain ~n^2, where y touches only the cliques of P3s at u and at v, plus the C0 condition,
plus direct interaction between two pointers only when they share a vertex (triple {u,v,v'} with
G-relation on vv') or lie in the environments of a common P_i. Adding l universal vertices gives a CSP
of treewidth <= w+l and l variables of huge domain; that is Theorem A again. A generic binary CSP of
this shape (l variables with domain n and tw <= l) is W[1]-hard, so a proof of FPT must use more
structure of Cluster Editing than the CSP shape. The structure that is available:
 * Lemma 2 (each pointer is a near-twin pair or a pair with few common neighbours),
 * closure: the cluster containing an inserted pair uv must also contain N(v), N(u) up to
   <= l compensating edits, which forces the partner sets to be "almost components".
What the DP would need: when it decides at node t(u) that u carries a half-edit with partner profile rho
(rho = adjacency of the partner to Q(u), one of 2^(2w+2) values), later it must find the other end v at t(v)
with a matching profile and a real pair uv. Profiles are finite but the matching of far-apart ends is the
identity problem: it is a pairwise relation, defined by G, between positions chosen by a tree-structured
(DP-definable) family. That is exactly the shape of Multicolored Clique, which is why I cannot rule
hardness out.

### (d) Wildcard P3s  FAILS as a literal reduction
A wildcard would be an extra induced P3 containing the unpacked edit pair, modification-disjoint from H.
By Lemma 1 we may take H maximal, and then no induced P3 of G has all pairs unpacked, so there is no
P3 of G to promote. Adding gadget vertices changes OPT in a way that is not controlled. In the C4 example
l = 1 is needed although the l = 0 CSP (3^1 values) is unsatisfiable, so the l = 0 CSP does not see the
unpacked edit at all. The correct "wildcard" is the pointer variable of (c) with n^2 values: Thm A.

## 3. Hardness: where it could come from, and why my naive reductions collapse
I looked for a W[1]-hardness reduction for (l, w) and found none. Attempts and why they fail:
 * Hub encoding (one vertex u, partners v_1..v_k pairwise adjacent, each uv_j inserted). Insertion of
   uv_j and uv_j' makes {u,v_j,v_j'} a P3 unless v_jv_j' is an edge, so the partners must be a clique,
   but Lemma 2 plus closure forces the partner set W to be almost closed under neighbours in G - N[u]
   (union of near-components), so the "selection" is not free. Deletion hubs give a Vertex Cover
   constraint (non-edges inside N(u) force one endpoint to be cut), which is FPT.
 * Chain encodings: positions along a path of I(H) carry only 2^O(w) states, so indices in [n] cannot be
   transported along I(H); comparison of two far-apart selected indices must go through G itself,
   i.e. through a shared vertex of two pointers. That is the only known source of pairwise coupling.
So: UNDECIDED. If a reduction exists it will probably take Multicolored Clique with k pointers
sharing hub vertices and a gadget forcing "partner sets are not closed" using packed P3s as filters
(each P_i can veto a pointer partner via its unary constraint). I leave this to the hardness agent.
A concrete, checkable intermediate question: for w = 1 (I(H) a path or forest), is CEaMP FPT in l?
Even the case "I(H) a perfect matching" (c=2, each P3 sharing one vertex with exactly one other) is
open, since c=2 is open in Li et al.

## 4. Test log (all pass; seeds in the files)
 T1 DP of Thm A == brute force (97 instances, l<=2)            code/fpt_tests.py
 T2 viability Lemma 2 on all 2594 solutions of cost <= |H|+l    code/fpt_tests.py
 T3 Lemma 1 equivalence; P3s of maximal H have >=2 vertices in V(H);
    non-edges inside N(r), r in R, are packed pairs (313 inst.)  code/fpt_tests.py
 T4 Lemma 3 locality of triples                                 code/fpt_tests.py
 T5 forced R-touching unpacked edits: 9 of ~400 (example C4);
    forced R-R edits: 0                                         code/fpt_tests.py
 T6 Lemma 4 on 1048 pairs of R vertices (n<=8)                  code/fpt_tests2.py
 T7 conjecture C1: 0 counterexamples in 142 instances           code/fpt_tests2.py
 C4 example printed by                                          code/fpt_examples.py
Limits: n <= 8, random sparse/medium graphs, random maximal packings; no test explores large w or
large l (instances are too small for interesting treewidth). Passing tests do not support C1 strongly.

## 5. Ranked summary
PROVED (with proofs above): Lemma 0, Thm A (n^(2l) 7^(w+1)), Thm B (3^l given x), Lemma 1, Lemma 2,
  Lemma 3, Lemma 4 and its corollaries (i)-(iii).
CITED, not re-proved: F1 (connected clusters), critical-clique property (iv).
CONJECTURE: C1 (no R-R unpacked edit needed); "CEaMP is FPT in (l,w)" (no evidence either way beyond
  the w = 0 case); "CEaMP is W[1]-hard in (l,w)" (no reduction).
FAILS: strong locality (a), wildcard reduction (d), branching on a violated P3 for w >= 1 (b): gives n^l.
NEXT (in order):
 1. Settle w = 1 / c = 2 with perfect-matching I(H) by writing the pointer CSP explicitly and trying to
    prove that pointers sharing no vertex can be handled independently (this would give FPT for forests).
 2. Hardness agent: attempt Multicolored Clique via shared-hub pointers (Section 3); first check whether
    the hub gadget can be made non-closed using filters from packed P3s.
 3. Try to prove C1 and then a bound on candidate partners for pointers whose partner is in R
    (twin classes, Lemma 4); this removes R from the identity problem.
 4. Search computationally, at n = 9..10 with structured instances (paths of P3s), for instances where
    the optimal pointer partner is not determined by its profile; that would show the DP-with-profiles
    idea is lossy.
