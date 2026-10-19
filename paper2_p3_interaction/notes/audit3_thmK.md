# Audit 3 (hostile referee, wave 3): Theorem K, FPT in (l, s) via BFK at component level

Sources read: README.md, notes/wave2_c1_w1.md (Sec. 3.2), notes/c2_algorithm.md (lines 148-160), and the BFK text
(arXiv 1512.04047, extracted text /tmp/claude-0/.../scratchpad/bfk/bfk.txt: Def. 1.2, Problem 1.3, Def. 3.1, Lemma 3.2,
Rule 6.2, Lemmas 6.3-6.5, Thm 6.6, Cor 6.7). Code: code/a3b_comp.py. Convention: I(H) has the P3s as nodes, two adjacent iff they
share a vertex (P3s of H are pair-disjoint, so only vertex-sharing is possible); V(C) = union of vertex sets; s = max |V(C)|.

## VERDICT (short)
Theorem K is CORRECT as a statement "FPT in (l, s)", and it can be stated with an explicit running time (below), not only
"follows from BFK". Three wording defects in the notes (none fatal): (i) "<= 2l' components with positive excess" is wrong as
worded, the true statement is "<= 2l' components that are irreducible under Rule 6.2"; (ii) t must be bounded by a function
of s explicitly (t <= C(s,2)/2); (iii) c2_algorithm.md says "<= 2l" but the sharp quantity is l' (<= l), and the exact
optimum tau_j is needed (not just a lower bound) because Rule 6.2 needs an optimal S of G[W_j].

## What BFK actually states (exact)
* Def 1.2: vertex-disjoint packing = set of INDUCED subgraphs H_i of G with mutually disjoint vertex sets. No connectivity,
  no size condition.
* Problem 1.3: input (G, H, k), H vertex-disjoint packing of induced subgraphs with 1 <= tau(H) <= t each, tau(H) = min size
  of an F-free editing set of H (for cluster editing: exact cluster-editing optimum of H). Lower bound h(H) = sum tau(H),
  parameter l := k - h(H). Parameters (t, l).
* Def 3.1: {u,v} is EXTERNAL for H iff exactly one of u,v lies in V(H). Every pair is external for <= 2 packing graphs.
* Lemma 3.2: if a size-k solution S has, for every H=(W,F) in H, (a) >= tau(H)+1 pairs of S inside W, or (b) one external pair
  of S for H, then |H| <= 2l and k <= (2t+1)l. Proof: k >= h_a+p_a + h_b + p_b/2, i.e. l >= p_a + p_b/2.
* Rule 6.2: H=(W,F) induced, optimal S of H (|S|=tau(H)) such that for u,v in W: same cluster of H△S => N_G(u)\W = N_G(v)\W;
  different clusters => N_G(u) ∩ N_G(v) ⊆ W. Then G := G△S, drop H, k -= tau(H).
* Lemma 6.3: rule is correct (proof only uses G, W, S; works for any induced subgraph, not only packed ones).
* Lemma 6.5: reduced yes-instance => k <= (2t+1)l.  Thm 6.6: (i) time O(Γ(G,(2t+1)l) + nm + sum_H Γ(H,t)); (ii) kernel
  (4t+2)l vertices. Cor 6.7: O(1.62^{(2t+1)l} + nm + n*1.62^t) and kernel (4t+2)l vertices, using Γ(G,k)=O(1.62^k+n+m) (Boecker).

## (1) Is the component packing a legal BFK packing, and is l' <= l?  CORRECT
Take C_1..C_r the components of I(H), W_j = V(C_j), H_j = G[W_j].
* Vertex-disjoint: a vertex in W_j and W_i would lie in a P3 of C_j and a P3 of C_i, which then share it and are adjacent in I(H),
  contradiction. Induced by construction. Connected in G as well (P3s are connected, chained by shared vertices), not needed.
* tau_j := tau(G[W_j]) = exact cluster-editing optimum (exponential only in s: exhaustive 3^s subset DP, or Boecker's 1.62^{tau_j}
  iterating k=0,1,..). 1 <= tau_j holds since C_j has a P3.
* tau_j >= |C_j|: the P3s of C_j are pair-disjoint induced P3s of G[W_j]; every solution edits >= 1 pair in each, all different.
* |H| = sum_j |C_j| (components partition H). Hence, with k = |H| + l:
    l' := k - sum_j tau_j = l - sum_j (tau_j - |C_j|) <= l,  equality iff every tau_j = |C_j|.
  If l' < 0 the answer is NO, justified by the lower bound OPT >= sum_j tau_j: for any solution S, G△S restricted to W_j is a
  cluster graph, so |S ∩ C(W_j,2)| >= tau_j, and the pair sets C(W_j,2) are pairwise disjoint (vertex-disjoint W_j).
  Both facts are exactly the direction needed. (t_max = max tau_j, and every P3 packing graph meets the hypothesis 1<=tau<=t.)

## (2) Do the BFK lemmas need induced / connected / bounded size?  CORRECT at component level
* Induced: needed (G△S restricted to W is H△S|_W only if H = G[W]); satisfied.  Connected: nowhere used.  Bounded size: used
  only through t (in k <= (2t+1)l and Γ(H,t)); size itself never enters.
* Lemma 6.5 core argument, re-derived: let S be optimal for G and H_j irreducible. If S has <= tau_j pairs inside W_j then exactly
  tau_j (lower bound above), so S∩C(W_j,2) is an optimal solution of H_j; by irreducibility it violates a condition of Rule 6.2.
  Violated cond. 1 (u,v same cluster, different outside neighbourhoods): if no external pair at u or v were edited, G△S would
  keep their outside neighbourhoods, so they could not share a cluster; hence an external pair is edited. Violated cond. 2
  (u,v different clusters, common outside neighbour x): without external edit at ux, vx, G△S has the induced P3 u-x-v. So
  (a) or (b) holds for every j. This uses only: pairs inside W_j are the pairs counted for tau_j, and "external" = exactly one
  end in W_j. A pair with one end in W_j and the other in R (unpacked) or in W_i, i != j, is external for j (and for i); the
  pair is counted at most twice in Lemma 3.2 (factor 1/2), which is exactly BFK's accounting. Pairs between two components
  are therefore handled. The needed non-collision "an external pair of j is not an inside pair of another i" is true because W_i, W_j disjoint
  (this is the point that FAILS for single P3s with pair-disjoint-but-vertex-sharing, a pair (u,v) external for P with u in P, v not
  in P can be an inside pair of P' and be a non-excess edit).
* Rule 6.2 on component j, soundness (Lemma 6.3) holds for any induced subgraph. After applying it, S_j changes only pairs inside
  W_j, so every other G[W_i] and its Rule-6.2 conditions (which involve only pairs incident to W_i) are unchanged; l' is unchanged
  (k and h both drop by tau_j). So one pass over all components suffices and l' stays the same.
* Wording defect: Lemma 6.5 bounds the number of IRREDUCIBLE components by 2l'. Components with "positive excess" is a different
  notion (tau_j > |C_j|) and no such bound holds for it (a reducible component can have large internal excess; it is removed
  at no cost to l').

## (3) Bound and explicit f(l, s)
Components remaining after exhaustive Rule 6.2: |H_red| <= 2l' <= 2l (BFK Lemma 3.2 + 6.5; ratio 2 attained, see (4)).
t = max_j tau_j <= min(m_j, C(|W_j|,2) - m_j) <= C(s,2)/2 (editing to the empty or the complete graph). Thus
  2t+1 <= C(s,2) + 1 = (s^2 - s + 2)/2,  kernel <= (4t+2) l' <= (s^2 - s + 2) l vertices.
Running time (Cor. 6.7 with t, l' ):
  O( 1.62^{(C(s,2)+1) * l} + n*m + n * 1.62^{C(s,2)/2} ) + (time to find components of I(H): poly)
  where the last term can be replaced by n * 3^s (subset DP) for the exact tau_j; so
  f(l, s) = 1.62^{(C(s,2)+1) l} + 1.62^{s^2/4} (times n^{O(1)}).  Doubly-in-s exponent s^2 l: FPT, not practical.
If s is instead measured as the number of P3s per component (s_P), then |V(C)| <= 2 s_P + 1 (a connected chain adds <= 2 new
vertices per P3), so the same holds with s = 2 s_P + 1. Instances with I(H) a perfect matching: |V(C)| = 5, t <= 5.
Nothing in the argument gives a lower exponent; the claim in Sec. 3.2 "t <= C(s,2)" is true but weaker than needed.

## (4) Computational test (code/a3b_comp.py, run: python3 -I a3b_comp.py 4000)
Random graphs n in 7..10, edge prob 0.2..0.6, random pair-disjoint induced-P3 packings (also thinned, non-maximal), components of I(H),
tau_j and OPT by exhaustive partitions (Bell(10) = 115,975 vectorised), 3,928 instances, 3,489 with a component of >= 2 P3s.
 * tau_j >= |C_j| : 4,163/4,163 components.
 * OPT >= sum tau_j : 3,928/3,928;  l' <= l (l = OPT-|H|, l' = OPT - sum tau_j): 3,928/3,928.
 * Decomposition cost(P) = sum_j inside_j(P) + #edited pairs not inside any W_j, and inside_j(P) >= tau_j, for EVERY partition P
   (hence every optimal one): 3,928/3,928. This is the "external pair is excess" accounting.
 * Rule 6.2 on components, with an optimal S_j chosen among all optimal solutions of G[W_j]: 3,240 applications, OPT dropped by exactly
   tau_j each time (3,240/3,240); no negative l' after reduction.
 * Reduced instances with irreducible components left: 847; #irreducible components <= 2l' in 847/847 (equality in 27);
   Lemma 3.2 property (each irreducible component has inside > tau_j or an external edited pair, in EVERY optimal partition, and
   OPT >= sum tau + p_a + p_b/2): 847/847.  Zero violations overall.
Limits: n <= 10, so s <= 10; evidence only for the lemmas, which are proved above independently.

## Verdicts per statement
1. Components of I(H) give a BFK packing of vertex-disjoint induced subgraphs G[V(C_j)], cost tau_j exact optimum: CORRECT.
2. sum tau_j >= |H| hence l' <= l (and l' < 0 => NO): CORRECT (proof above).
3. BFK lemmas need neither connectivity nor bounded size (only t); external-pair-is-excess holds at component level incl. pairs
   between components: CORRECT. (Does not transfer to single P3s of a non-matching packing: known.)
4. "# components with positive excess <= 2l'": WRONG AS WORDED. Fix: "number of components that are irreducible under Rule 6.2 is
   <= 2l'; all others are removed by Rule 6.2 without changing l'". The notes' numerical claim (#remaining components <= 2l')
   is right.
5. f(l,s): GAP (only in specification). Fix: t <= C(s,2)/2, time 1.62^{(C(s,2)+1) l'} n^{O(1)} + n*min(3^s, 1.62^{C(s,2)/2}),
   kernel <= (s^2 - s + 2) l vertices.
Statement for the paper: "Theorem K. CEaMP can be solved in time 1.62^{(C(s,2)+1)l} * n^{O(1)} + n * 3^s and admits a kernel with
<= (s^2-s+2) l vertices, where s = max_C |V(C)|. Proof: apply BFK Thm 6.6/Cor. 6.7 to the packing {G[V(C)]} with t <= C(s,2)/2 and l' <= l."
That is a citation-based proof (BFK Thm 6.6 and Lemma 3.2 are used as black boxes, hypotheses verified in (1)-(2)), which is
legitimate; the proof of the two hypotheses (vertex-disjoint, induced, tau exact, 1<=tau<=t) is what the paper must include.
Minor residual caveat: Thm 6.6 as stated runs Γ(G,(2t+1)l) on the reduced graph after kernelisation; BFK's Γ-assumption
(monotone, polynomial) holds for Boecker's algorithm. Not an issue.
