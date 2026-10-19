# c2_algorithm: what the known algorithms do, and what FPT in (l, tw I(H)) would need
Sources: LPS = arXiv 1910.08517v3, cited as `lps.txt:LINE` (scratchpad/lps/lps.txt).
BFK = van Bevern-Froese-Komusiewicz arXiv 1512.04047 (downloaded to scratchpad/bfk/, text
`bfk.txt`; line numbers refer to `pdftotext -layout` output). INFERENCE = my own, not in either paper.
Numbering clash: "LPS Obs 2" (lps.txt:1163) is a different statement from "Observation 2" of
analysis_01.md. I write "our Obs 1/2" for the analysis_01 ones.

## 1. Structure of the algorithms
### 1a. LPS c = 2 (half-integral), Theorem 2 (n^{2l+O(1)})
Two layers: Algorithm 1 (lps.txt:1089-1107) guesses everything excess; Theorem 3 (l = 0,
polynomial, lps.txt:1070-1071, proof Section 5.1) solves the rest.
Algorithm 1, for l_a = 0..l, l_b = 0..l-l_a:
* guess S_a: l_a vertex pairs {u,v} with |{u,v} ∩ V(P)| <= 1 for all P in H (line 1096), i.e. pairs not
  covered by any packed P3 (every pair inside V(P) is covered by P). Set G_a = G△S_a.
* guess H_b: l_b distinct packed P3s that get >= 2 edits "in V(P)" (line 1098).
* guess S_b: for each P in H_b a set of >= 2 of its 3 pairs (4 options per P, line 1083). G_b = G_a△S_b.
* H' = H \ H_b; ask whether G_b has a cluster-editing set of size |H'| (zero excess), by Thm 3.
* accept if |S_a|+|S_b| <= |H_b| + l (line 1101).
Correctness (lines 1110-1123) uses only: a solution S* splits as S_a* (edits on pairs in no packed
P3), S_b* (all edits in P3s having >= 2 edits), S_0* (the rest). Each P outside H_b* has exactly one
edit (>= 1 needed, <= 1 by definition of H_b*). Half-integrality is NOT used in this part (INFERENCE:
read off the proof; it is used only for the running time count and for Thm 3).
Running time (lines 1080-1084): "at most n^{2 l_a} choices for S_a"; 3|H| <= 2n (c = 2) so |H| <= n and
"O(n^{l_b}) choices for H_b"; 4^{l_b} for S_b; total O(4^{l_b} n^{2 l_a + l_b + O(1)}) <= n^{2l+O(1)}.
Where n^{2l} comes from: the S_a guess, 2 endpoints each from n vertices, l_a times. The H_b guess
contributes n^{l_b} (INFERENCE: a second XP factor, also not FPT; the paper's text only highlights S_a).
Which step blocks FPT: both guesses are blind enumerations; no rule bounds the candidate set for
S_a or H_b by a function of l. The paper says why the BFK approach does not port (lines 142-147):
"the property that packed P3s have an edit incident with them is not helpful anymore, because these
edits may be part of other packed P3s and hence not be in excess. Second, ... we need to check the
special neighborhood properties of van Bevern et al. for arbitrarily large connected components of
packed P3s". Conclusion (lines 2054-2057): "we discovered large grid-like structures that seemed
difficult to overcome in fixed-parameter time, and a corresponding W[1]-hardness result would also not
be surprising." So FPT for c = 2 is open and W[1]-hardness is considered plausible by the authors.
Theorem 3 (zero excess, c = 2), Section 5.1. Language: proto-cluster = connected component of the graph
of unpacked edges (line 1129; also Sec. 3 lines 203-207). P3 types alpha/beta/gamma/delta by how
V(P) meets proto-clusters (1152-1157). Foundation: LPS Obs 2 (1163): "A solution S ... must edit
exactly one edge or non-edge of every P3 of H, and neither non-packed edges nor non-packed non-edges
can be edited by S." (same as our Obs 1.) Rules, all "safe" = equivalence-preserving at l = 0:
* RR1 (1168): two vertices of one proto-cluster with an unpacked non-edge -> NO. (Lemma 5)
* RR2 (1177): type-beta or type-delta P3 xyz -> insert xz, drop P3. (Lemma 6)
* RR3 (1216): unpacked non-edge uv between proto-clusters A,B and packed edge xy between A,B ->
  delete xy, drop its P3. (Lemma 7)
* RR4 (1238): isolated clique proto-cluster removed. (Lemma 8)
* RR5 (1249): connected component of G of size <= 6: brute force. (Lemma 9)
* RR6 (1372): size-3 proto-cluster C with B, A singletons (structure from Lemma 13). (Lemma 14)
* RR7 (1692): five-proto-cluster configuration {z1,z2},{x},{y},{u},{v} (structure from Lemma 20),
  nine items (1700-1728). (Lemma 21)
Consequences: no proto-cluster of size >= 5 (Lemma 10, uses "v incident with at most two packed P3s"
at 1277-1278), none of size 4 (Lemmas 11, 12), none of size 3 (Lemma 13, RR6) hence all proto-clusters
have size <= 2 and all packed P3s are type gamma (Corollary 2, 1421-1423); solution clusters have size
<= 4 (Lemmas 15-22); then only deletions are needed (Lemma 23) and the instance is a 2-SAT
(Lemma 24, 1941-1946): one variable per packed edge, for each P3 (x_i or x_j) and (not x_i or not x_j),
and (x_p or x_q) for each induced P3 uvw whose two edges are in distinct packed P3s.
Why half-integrality is the engine (INFERENCE from Lemma 10): with c = 2, two big proto-clusters cannot
have all cross pairs covered, so some unpacked non-edge exists and RR3 fires. With c >= 3 or so the
padding of Lemma 1 can cover all pairs, which is exactly what the hardness uses.
### 1b. BFK, Cluster Editing above a VERTEX-DISJOINT P3 packing (FPT in l)
Two algorithms.
(i) Generic win-win (Section 3, Section 6.1): Reduction Rule 6.2 (bfk.txt:795-797): for an induced
subgraph H=(W,F) with an optimal solution S of size tau(H) such that for all u,v in W:
N(v)\W = N(u)\W if u,v are in one cluster of H△S, and N(u)∩N(v) ⊆ W otherwise, replace G by G△S,
drop H from the packing, k -= tau(H). Lemma 6.3 (816-890): there is an optimal solution containing S
(exchange argument with "good" vertices/equivalence classes of external neighbourhoods).
If no packing graph admits the rule, then any optimal solution S has for each H either (a) more than
tau(H) pairs inside W or (b) an external pair (exactly one endpoint in V(H)) (Lemma 6.5, 950-958).
Lemma 3.2 (295-322): every pair is external for at most two packing graphs (they are vertex-disjoint),
so k >= h_a + h_b + p_a + p_b/2, i.e. |H| <= 2l and k <= (2t+1) l. Then run the O(1.62^k) algorithm
(Thm 6.6, Cor 6.7: O(1.62^{(2t+1)l} + nm + n 1.62^t)).
(ii) P3-packing-specific search tree, Thm 6.13: O(4^l l^3 + n + m). Rules 6.8, 6.9 (special cases of
6.2), 6.10 (clique K with all vertices having <= 1 outside neighbour: cut it off), 6.12 (isolated
clique). Branching Rule 1: an induced P3 that contains at most one vertex of each packed P3 -> branch
3 ways, each edit is excess, l drops by 1. Branching Rules 2-4 on a packed P3 uvw with a neighbour x
(vectors (1,1,1,2), (1,1,1,1), (1,1,1,1)). End case: max degree 2, polynomial (Cases I-IV).
How excess edits touching unpacked pairs are found without guessing (BFK): in the vertex-disjoint case
an external pair is automatically in no packed P3, so it is an excess edit (INFERENCE: this is the
reason Lemma 3.2 works; the pair has one endpoint in V(P) and so lies inside no V(P'), P' != P).
In the modification-disjoint case an edit with one endpoint in V(P) can be a packed pair of another P3
(LPS lines 142-144), so the counting k >= h + |H|/2 breaks.

## 2. Excess edits and unpacked pairs
LPS: all handled by brute force. Edits on unpacked pairs are fixed in S_a before anything else;
afterwards (G_b, H') must be solved with no excess, so by LPS Obs 2 no unpacked pair is ever edited
again. P3s receiving 2 or 3 edits are removed from the packing (H_b) with their edits applied.
Remaining H' is still modification-disjoint and half-integral because applying edits to pairs of
removed P3s or to unpacked pairs changes no pair of a P3 in H' (INFERENCE, immediate).
BFK: unpacked edits are never enumerated. Either Rule 6.2 consumes a P3 at cost exactly 1 (so no excess
used), or the P3 forces an external edit (counted as excess, charged to at most 2 packing graphs), or
Branching Rule 1 hits an induced P3 with all three pairs unpacked.
INFERENCE (usable in CEaMP for any c): Branching Rule 1 carries over. An induced P3 of G whose three
pairs are all unpacked must receive an edit on one of its 3 pairs, every such edit is excess, so
3-way branch with l -> l-1. After exhaustive branching, every induced P3 of G has >= 1 packed pair.
This does not yet bound the S_a candidates: voluntary excess edits (made to avoid a worse packed
edit, not forced by an all-unpacked P3) remain and the candidate pairs are still O(n * |V(H)|).

## 3. Hardness (Theorem 1: NP-hard for l = 0, every vertex in <= 23 packed P3s)
Reduction from 3-SAT (lines 250-259). Goal (G,H,0) is YES iff one edit per packed P3 suffices.
Principle (Sec. 3, 202-246): proto-clusters (cliques joined by unpacked edges) must be merged into
larger clusters; unpacked non-edges forbid merges; a packed P3 with one vertex in T and two in R
(type alpha) can absorb either "merge" or "separate" at cost 1, which gives completeness.
* Variable gadget (282-331): cycle of 4 m_i cliques K^i_j of 5 vertices; for even j a P3 for every
  arithmetic progression q-p = r-q over F5 on (K_j, K_{j+1}, K_{j+2}), centre in K_{j+1}; this covers
  every pair between consecutive cliques exactly once, so Hvar are modification-disjoint. A solution
  must merge all even pairs or all odd pairs (Lemma 4 proof, 976-1011): truth value.
* Clause gadget (386-529): cliques Q1..Q4 in a path, unpacked non-edge Q1Q4, so one of the three
  consecutive cuts happens (three-way choice = which literal satisfies). Transferring cliques T_a,T_b,
  T_c attached by two P3s through shared centre in Q1 (resp. Q3, Q4); T attached to variable cliques by
  four P3s with centres in K_{4pi+1} (Fig. 3).
* Padding for completeness (531-830): every pair between cliques that might be merged must be covered
  by a packed P3. Merging model with levels L0..L4 (614-620), Lemma 1 (672-681): V of size p prime, W of
  size 2p, F ⊆ E a union of singletons/P3s centred in V/C8s; there is an edge-disjoint triangle packing
  covering E\F (triangles v_i w_j w'_k with j-i = k-j in F_p), W stays connected after removing used
  edges (so W stays a proto-cluster), each vertex in <= p triangles. Triangles become P3s with an edge
  inside the padded clique and one cross pair.
* Why 23 (788-830): p = smallest prime >= |V| with 2p >= |Q|. Q1,Q4: 1 vertex; Q3: p=2, 4 vertices;
  Q2: |V|=6, p=7, 14 vertices; T_a,T_c: |V|=16, p=17, 34 vertices; T_b: |V| = |Q3|+|Q4|+3*5 = 20,
  p=23, 46 vertices (824-829). Hence "each vertex is in at most 23 P3s". INFERENCE: 23 is the first prime
  >= 20, an artefact of the construction, with no claim of tightness; the paper leaves 3 <= c <= 22
  open (lines 127-128). I did not re-verify the bound for level-0 vertices (K cliques) that are
  out-neighbours of several T's; the paper asserts it at 829-830.
* Completeness (Lemma 2) and soundness (Lemmas 3, 4) as above; Lemma 3: V(H) (the cliques) are exactly the
  proto-clusters; Observation 1 (963): clusters of any solution are unions of these.

## 4. What an FPT proof in (l, tw I(H)) must establish; reusable parts
INFERENCE throughout this section.
Setting: w = tw(I(H)). Our Obs 1 and 2 give a 3^{w+1} n^{O(1)} dynamic program for (G,H,0).
A. Free XP result. Combine LPS Alg 1 (generalised: the correctness argument at lines 1114-1123 has no
use of c, only the count of H_b choices uses |H| <= n^2) with our DP in place of Theorem 3:
time 3^{w+1} * n^{2 l_a + 2 l_b + O(1)} * 4^{l_b}. The sub-packing H' has I(H') an induced subgraph of
I(H), so tw(I(H')) <= w; G_b changes no pair owned by H'. So (l, w) gives XP in l, FPT in w, for every
c. (Not stated in either paper; needs our Obs 1/2 written up as proved lemmas.) The open task is to
remove the n^{2l}.
B. To prove FPT one must show one of:
 (1) a polynomial (or f(l,w)) bound on the number of candidate unpacked pairs for S_a and candidate P3s
     for H_b, e.g. via an exchange argument; or
 (2) a dynamic programme over a tree decomposition of I(H) whose states also account for excess edits
     with budget l. Obstruction: an unpacked pair uv lies in the scope of the n-2 triples uvx, so
     the cluster structure of all vertices outside V(H) would have to be carried (Obs 1 in analysis_01
     lists this); or
 (3) a win-win with a measure like BFK Lemma 3.2 (below).
 Statements to establish in any case: normal form (S_a, H_b, S_b, S_0) [LPS 1114-1123, reusable]; the
 zero-excess DP [our Obs 1/2, LPS Obs 2 line 1163 supports]; an exchange or kernel lemma for unpacked edits.
C. Reuse candidates.
 * LPS RR1, RR2, RR3, RR4, RR5 (Lemmas 5-9) are valid for arbitrary modification-disjoint packings at
   l = 0: their proofs use only LPS Obs 2 and the safeness of RR1 (INFERENCE from reading the proofs).
   RR6, RR7 and Lemmas 10-22 use the "at most two packed P3s per vertex" count, so they are c = 2 only.
   None of the rules is valid with excess budget: every proof relies on LPS Obs 2 (exactly one edit per
   packed P3, none elsewhere). They are usable only after guessing, on (G_b, H').
 * BFK Rule 6.2 / Lemma 6.3 is stated for any induced subgraph H of G and is a statement about
   Cluster Editing alone, not about vertex-disjointness. INFERENCE: applying it to a P3 P of a
   modification-disjoint packing is safe (G△S is an equivalent instance, the other P3s keep all their
   pairs), what fails is Lemma 6.5 (bounding k).
 * Key observation (INFERENCE): connected components of I(H) are vertex-disjoint. Take the packing of
   induced subgraphs G[V(C)] for components C of I(H). BFK Lemma 3.2 and Lemma 6.5 apply verbatim to this
   packing (cost t = max_C tau(G[V(C)]), lower bound h = sum_C tau(G[V(C)]) >= |H|, so the new excess is
   <= l): if Rule 6.2 does not apply to any component, then |{components}| <= 2l. This gives FPT in
   (l, t) by BFK Thm 6.6, but t is unbounded here (a component can be large), so it does not give FPT
   in (l, w). It does suggest the plan: reduce to <= 2l components, each of treewidth <= w in I(H) but
   arbitrary size, then handle the interaction of <= 2l blobs plus the unpacked vertices X = V \ V(H).
   Checking Rule 6.2's condition for a component needs an optimal solution of G[W] of size |C| with
   prescribed neighbourhood behaviour, plausible by our DP (also INFERENCE; not proved).
 * LPS Lemma 1 (padding via F_p) is the tool for any W[1]-hardness/gadget construction that needs every
   pair between proto-clusters covered by modification-disjoint P3s; the variable gadget (291-331) is a
   bounded-degree cyclic structure. A W[1]-hardness proof in (l, tw) would have to simulate the choice
   of l unpacked pairs among n^2 on a bounded-width host: the authors' "grid-like structures" (2056).
 * BFK Branching Rule 1 (INFERENCE, see Section 2) removes forced unpacked edits.
D. Consistency checks for any claim:
 * c = 1: I(H) edgeless, w = 0; BFK gives FPT in l, so the (l, w) parameter cannot be proved hard by
   w alone.
 * c = 23 hardness instances at l = 0 must have unbounded tw(I(H)) unless P = NP (our DP). I(H) there has
   max degree <= 3*22 = 66.
 * c = 2, l = 0 is polynomial (Thm 3) although I(H) is subcubic (degree <= 3) and can have large treewidth,
   so tw(I(H)) is not necessary for tractability: the 3^{w} DP is not tight for c = 2 (INFERENCE).

## 5. Observation 2 of analysis_01 (scopes at l = 0 are cliques of I(H))
Supports:
* LPS Obs 2 (1163) is exactly the premise (one edit per packed P3, nothing else), and the constraint
  satisfaction format is implicit in Sec. 3 (proto-clusters, constraints on "merged or separated",
  lines 205-216). The 2-SAT of Lemma 24 is a triple-constraint CSP on boolean domains, its clauses
  (x_p or x_q) for P3s u-v-w with uv, vw in distinct packed P3s have scopes sharing the vertex v,
  i.e. adjacent in I(H).
* No statement in LPS mentions treewidth or any parameter of the intersection graph (grep: only "grid"
  at 2056), nor contradicts Obs 2.
Apparent tension, resolved (INFERENCE): Fig. 1 and its text (lines 229-240) describe a P3 A-B-C and a P3
C-D-E interacting "non-locally" through proto-cluster C. The two packed P3s need not share a vertex (they
may use different vertices c, c' of C), so this looks like interaction between vertex-disjoint P3s,
contradicting Obs 2's consequence. It does not: a triple {x,y,z} has all its pairs touching, so every
single constraint has a clique scope. The interaction in Fig. 1 holds only if the cross pairs between
B and C (b with the other vertices u of C) are covered by packed P3s; if bu is an unpacked non-edge the
triple (b,c,u) is a unary constraint forcing the deletion of bc. In the hardness construction all
cross pairs are covered (padding, 277-278, 641-670), and the covering P3s link the two original P3s in
I(H) (P - pad - P' share b and c'), i.e. the non-locality becomes paths and cliques in I(H). This is a
reason hardness needs huge cliques in I(H) (T_b has 46 vertices each in 23 P3s). The README's remark that
I(H) alone suffices at l = 0 is consistent with this.
Caveat: Obs 2 is a statement about single constraints only; the width bound 3^{w+1} follows from it by
the standard CSP DP (constraint graph has the same cliques as I(H) so a tree decomposition of I(H) is
one of the primal graph). Our code check (486 random instances) confirms; LPS does not check it.
