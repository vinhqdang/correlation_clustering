# Hardness and counterexamples (negative direction)

Labels: CLAIM (proved here), VERIFIED-BY-COMPUTER (small cases only), IDEA, OPEN.
Notation as in analysis_01.md: H packing, I(H) intersection graph, w = tw(I(H)), c = max number
of packed P3s through a vertex, l = excess budget. Scripts: code/hard_*.py.

## 0. Summary

| question | status |
|---|---|
| 2^{o(w)} n^O(1) impossible at l=0 (hence 3^{o(w)}) | CLAIM under ETH, via the size of the Li et al. reduction (A1) |
| (3-eps)^w impossible under SETH at l=0 | OPEN. No pathwidth-preserving reduction found (A2, A3) |
| XP in (l, w) | CLAIM (B2): n^{2l} * 8^{w+1} * poly, any c |
| FPT or W[1]-hard in (l, w) | OPEN. Obstruction identified (B3) |
| w alone (even w = 0) without bounding l | CLAIM: NP-hard (C1) |
| l = 0, c alone | NP-hard for c >= 23 (Li et al.); c = 2 poly; 3 <= c <= 22 OPEN (C2) |
| tw(G) and tw(I(H)) incomparable | CLAIM (C3) |
| tw(G) alone | OPEN (C4) |

## A. Question 1: l = 0

### A1. CLAIM (ETH lower bound 2^{o(w)}, hence 3^{o(w)}, for l = 0)
Statement. Assume ETH. No algorithm solves CEaMP with l = 0 in 2^{o(w)} n^O(1), where it is given
(G, H) and a tree decomposition of I(H) of width w. Since 3^{o(w)} = 2^{o(w)}, this also excludes
3^{o(w)} n^O(1). It holds even for c <= 23.

Proof. Take 3-SAT with N variables and O(N) clauses (no 2^{o(N)} algorithm under ETH, by the
Sparsification Lemma). Apply the reduction of Li, Pilipczuk, Sorge (arXiv 1910.08517, Sec. 4.1).
Sizes as read from Sec. 4.1.1 and 4.1.4 (not re-derived by me):
* variable x_i gives 4 m_i cliques of 5 vertices plus the P3s with q-p = r-q over F_5 (25 per
  triple of consecutive cliques): O(m_i) vertices and P3s;
* each clause gadget has clique sizes 1, 1, 4, 14, 34, 46, 34 (constants), is joined to three variable
  gadgets by a constant number of transferring P3s, and is padded by the triangle packing of their
  Lemma 1 / Corollary 1 with p <= 23, hence O(1) padding P3s;
* every vertex lies in <= 23 packed P3s; Phi is satisfiable iff (G, H, 0) is a YES instance.
So |V(G)| and |H| are O(N). The one-bag decomposition of I(H) has width |H|-1 = O(N); a
2^{o(w)} n^O(1) algorithm fed this decomposition would decide 3-SAT in 2^{o(N)}. QED.

Weakness. This is the cheap bound: it uses only w <= |H|. It gives no base. The risk is my reading
of the size accounting (padding counts, Sec. 4.1.2 connection P3s); it must be re-checked line by
line before it goes into the paper.

### A2. IDEA (why a base is hard to get from the Li et al. instance)
I(H) of their instance is, up to constants, the incidence graph of Phi with each variable replaced by
a cycle through its occurrences. Ordering occurrences along a path decomposition of the incidence
graph gives pw(I(H)) <= O(pw(Phi)) with a multiplicative constant (each active variable keeps the
cycle's first vertex and current end in the bag, plus clause gadget bags). A multiplicative constant
transfers only 2^{Omega(w)}, never a base; base 3 needs an additive O(1) loss. Not done.

### A3. The local constraint language (VERIFIED-BY-COMPUTER; evidence only)
Scripts: code/hard_rel.py, hard_closure.py, hard_tern.py, hard_tern2.py, hard_eq.py.
Two P3s share at most one vertex (their pairs are disjoint). I enumerated every role placement of the
shared vertex and every fixed pattern of the 4 cross pairs (5-vertex graphs, cluster test by exhaustion):
* 32 distinct binary relations on {0,1,2}^2 (0: delete ab, 1: delete bc, 2: add ac), sizes 0,1,2,3,5.
  The 2-element ones are partial bijections such as {(0,0),(1,1)}; the 5-element ones have the shape
  "x = 0 or y = 0". None has more than 5 of the 9 tuples.
* Triangle gadget: vertices y,x,z with yx in P, xz in Q, yz in R gives a ternary constraint on three
  P3s pairwise sharing a vertex (6 vertices, each vertex in at most 2 P3s). Over all 216 role placements
  and 2^6 cross patterns: 164 distinct ternary relations, sizes 0..4 out of 27. For all 164, an
  independent set-partition enumeration (n = 6) agrees that OPT = 3 iff the relation is non-empty
  (0 mismatches); this re-verifies Observation 1 on gadgets.
* Closing the binary catalogue under composition, intersection and transpose (132 relations) gives
  neither equality nor "not equal". Adding projections and slices of the ternary relations (289
  relations) gives equality and one non-trivial bijection (swap 0<->1, 2 fixed), but still no
  "not equal" and none of the 6 permuted forms y != pi(x).
Reading. This closure uses only binary pp-operations, so it proves nothing. It suggests that a direct
3-Colouring reduction ("adjacent vertices differ" as one bounded gadget between two P3s) is not
available, so a SETH-tight bound would need a cleverer encoding (for example the cycles of
proto-clusters of Li et al., where the Boolean choice is propagated by transitivity). The sparsity
(at most 5 of 9 binary tuples, at most 4 of 27 ternary tuples) also leaves open that a
(3-eps)^w algorithm exists. OPEN.
Consistency with Li et al.: all these gadgets have c = 2, and c = 2 at l = 0 is polynomial
(n^{2l+O(1)}), so the absence of NE there is what one expects. Plausibility only.

## B. Question 2: (l, w) jointly

### B1. CLAIM (edits on packed pairs only: 8^{w+1} for any l)
Restrict to solutions S contained in the 3|H| packed pairs. Then min |S| is computable in
8^{w+1} n^O(1) for every l.
Proof. Variable x_i = S ∩ pairs(P_i), 8 values. For each triple t of vertices, the unpacked pairs are
fixed and the packed pairs have owners that pairwise share a vertex (Observation 2 of analysis_01, checked
by computer there), so the constraint "t is not an induced P3 after editing" has scope a clique of I(H).
The cost is sum |x_i|, a unary term. A min-sum dynamic programme over a width-w tree decomposition of
I(H) costs 8^{w+1} times the number of constraints. QED.

### B2. CLAIM (XP in l for every c: n^{2l} * 8^{w+1} * poly)
Proof. Let S be a solution with |S| <= |H| + l. Every P in H is an induced P3 of G, so S contains an
edit among pairs(P) (otherwise P survives in G + S); these pair sets are disjoint. Hence S has >= |H|
packed edits and |S ∩ unpacked| <= l. Guess U = S ∩ unpacked: at most n^{2l} choices. Flip the pairs of
U in G. A pair in U is not a pair of any P in H, so each P is still an induced P3 and H is still a
modification-disjoint packing of the new graph. Run B1 with budget |H| + l - |U|. QED.
VERIFIED-BY-COMPUTER (code/hard_star.py, 254 random graphs, n <= 7): in every optimal partition each
packed P3 receives at least one edit and #unpacked edits <= OPT - |H|.
This recovers the c = 2 XP bound of Li et al. when w is bounded and extends to all c. What remains is
only FPT versus W[1]-hard.

### B3. OPEN: FPT or W[1]-hard in (l, w)?
Known: w = 0 (c = 1) is FPT in l (van Bevern et al.); H empty is Cluster Editing, FPT in k = l. Both
are consistent with FPT in (l, w).
Why no conflict-free FPT argument yet. An unpacked excess edit uv with u, v in no common P3 lies in the
scope of every triple {u, v, x} with x adjacent to u or v. These constraints are between the guessed
pair uv and the owners of ux, vx, so I(H) does not see them. After guessing U, B1 applies; the guess is
the only obstacle. An FPT proof needs (a) a candidate set for each unpacked edit of size depending on
(l, w) only (an exchange argument), or a different handling of unpacked edits.
Why W[1]-hardness is not obvious either. A Multicolored Clique reduction would code "colour class i
picks v_i" by a guessed unpacked pair and the adjacency tests by P3 constraints with w independent of
k. The only interaction between two guessed pairs passes through transitivity of clusters, which in Li
et al. needs long chains (unbounded w). I found no encoding. IDEA only; no gadget was built or tested.

## C. Question 3: counterexamples to weaker conditions

### C1. CLAIM (w alone is not enough, not even w = 0)
H empty is a valid packing, I(H) is empty, w = 0, l = k is the Cluster Editing budget: NP-hard
(Shamir, Sharan, Tsur 2004). Stronger: let H be a maximal vertex-disjoint P3 packing of G; then c = 1,
I(H) has no edges, l = k - |H| and the instances are again arbitrary Cluster Editing instances.
Conversely l = 0 with c <= 23 is NP-hard (Li et al.), so bounded l alone is not enough either.

### C2. Weaker condition (i): bounded degree of I(H) / bounded c
* l = 0, c <= 23: NP-hard (Li et al. Thm 1); I(H) has max degree <= 3*22 = 66.
* l = 0, c = 2: polynomial (n^{2l+O(1)}). So tw(I(H)) is sufficient but not necessary at l = 0 if
  c = 2 instances can have large tw(I(H)); I did not construct such instances, so do not assert it.
* l = 0, 3 <= c <= 22: OPEN in Li et al.; nothing new here.
Combined with A1: bounded degree of I(H) does not prevent 2^{Omega(w)}.

### C3. CLAIM (tw(G) and tw(I(H)) are incomparable)
(a) Star K_{1,2m}, hub v, P_i = (a_i, v, b_i), i = 1..m: pairs va_i, vb_i, a_ib_i are disjoint across
i, so H is modification-disjoint; all P_i share v, so I(H) = K_m with tw m-1, while tw(G) = 1. Here
OPT = 2m-1 (keep v and one leaf), so l = m-1 (VERIFIED-BY-COMPUTER for m = 2, 3: OPT = 3, 5).
(b) k x k grid with a vertex-disjoint family of induced P3s (three consecutive grid vertices on a row;
induced since a grid has no triangles): I(H) has no edges, tw 0, tw(G) = k.
Consequence: hubs force large tw(I(H)); a theory in tw(I(H)) must treat hubs separately.

### C4. OPEN: is bounded tw(G) enough at l = 0?
Not excluded. Observations only: the final graph has edges in E(G) plus chords ac of packed P3s; I
could not bound the treewidth of G plus those chords in general (for a star it stays small). The Li et
al. instances have cliques of size up to 46 but unbounded tw(G), so they do not decide it. No claim.

### C5. Bounded tw(I(H)) alone at l > 0
Not enough by C1. What holds: XP in l (B2) and the packed-edits-only variant in 8^{w+1} (B1). The
failing ingredient is exactly the unpacked excess edits.

## D. Statements usable in the paper
1. Proved: l = 0 gives 3^{w+1} n^O(1) (analysis_01, Obs. 2).
2. Proved: n^{2l} 8^{w+1} (B2) and 8^{w+1} n^O(1) for packed-pair edits (B1).
3. ETH: no 2^{o(w)} n^O(1), even for c <= 23 (A1; recheck sizes). Do not claim a SETH base.
4. Necessity of l (C1), incomparability with tw(G) (C3).
5. Open: FPT in (l, w); c = 2 FPT in l; tw(G) at l = 0; SETH base; 3 <= c <= 22 at l = 0.
