# Wave 2, negative side: is CEaMP W[1]-hard in (l + tw(I(H)))?
Agent "w2b", 2026-10-19. Scripts: code/w2b_sm.py (SM, MAX, HE checks), code/w2b_hub.py (HUB check).
Labels: PROVED (proof given), VERIFIED-BY-COMPUTER (small cases, evidence only), IDEA, FAILS.
Notation as in analysis_01 / fpt_attempt: H packing, packed pair = "soft", every other pair = "hard"
(hard edge / hard non-edge), D = max number of soft partners of a vertex, l excess budget, w = tw(I(H)).

## 0. Verdict in five lines
1. I found NO W[1]-hardness reduction, and I now have PROVED structural lemmas (SM, HUB, MAX below) that explain
   exactly why the Multicolored-Clique hub encodings collapse and that make me believe the problem is FPT in
   (l, w). I do not claim either direction as a theorem.
2. Correction of the framing in fpt_attempt.md Sec. 3 ("w independent of k"): hardness for the parameter l + w
   only needs l + w <= g(k). A reduction may use w = O(k) or O(k^2), so grids of side k in I(H) are allowed.
   The obstruction is therefore not "w is too small for a hub", it is the lemma SM below.
3. "+c" is a vacuous weakening: c <= w+1 (the P3s through one vertex pairwise share it, so form a clique of I(H)),
   hence l + w + c and l + w are the same parameter up to a function. PROVED, trivial. No weakening found
   that is hard. (tw(G) as extra parameter: see Sec. 6, only IDEA.)
4. The WLOG-maximal-H warning of fpt_attempt Lemma 1 is removable: MAX shows H can be made maximal with
   w' <= w + l. PROVED (Sec. 2, A). This is a useful new tool for the positive side too.
5. No n^{Omega(l)} ETH bound for fixed w found; SM suggests none exists if FPT holds. Only cheap bounds (Sec. 5).

## 1. Basic facts used (all PROVED)
(F0) D <= 2c <= 2(w+1): a vertex u lies in <= c packed P3s, each contributes 2 soft partners.
(F1) tw(W) <= 3(w+1) - 1 for the soft graph W = (V, packed pairs). Proof: bag of I(H) -> union of the vertex
     sets of its P3s (<= 3(w+1) vertices). Every soft pair lies in a P3, hence a bag. For a vertex v the P3s through
     v form a clique of I(H), so (Helly) lie in a common bag, and the subtrees of its P3 nodes pairwise intersect,
     so their union is connected. QED
(F2) (Lemma 0, fpt_attempt) |U| <= |S| - |H| where U = S minus packed pairs. So |U| <= l for a solution.
Pieces: for a solution S let F = U cap E (deleted hard edges), I = U \ E (inserted hard pairs). The pieces of S are the
connected components of the graph (V, E \ packed \ F) of kept hard edges. Each piece lies inside one cluster.

## 2. The three structural lemmas

### A. MAX (maximalisation costs only +l in treewidth)  PROVED, VERIFIED-BY-COMPUTER (1500 instances, n<=8)
Let M be a maximal set of induced P3s of G whose pairs are pairwise disjoint and disjoint from pairs(H).
Then OPT <= |H|+l  iff  ( |M| <= l and OPT <= |H u M| + (l-|M|) ), and tw(I(H u M)) <= w + |M| <= w + l.
Proof. H u M is a modification-disjoint packing of induced P3s, so OPT >= |H|+|M| (each packed P3 needs a
distinct edit). So |M| > l implies NO; otherwise the two budgets coincide. Adding |M| nodes to I(H) adds <= |M| to
its treewidth (put them in every bag). QED  (M is found greedily in polynomial time.)
Consequence: WLOG no induced P3 of G has all three pairs hard. This gives, for free:
(HE) hard edges ab, bc imply ac in E or ac soft. (An induced P3 a-b-c with ab, bc hard would be all-hard unless ac
     is soft.) VERIFIED-BY-COMPUTER together with MAX.
So the "all-hard P3" branching of BFK Rule 1 (3^l) is subsumed: it is exactly the enlargement H -> H u M.

### B. HUB (deletion hubs / vertex-cover gadgets have bounded size)  PROVED, VERIFIED-BY-COMPUTER (500 inst.)
If u has hard neighbours X (ux in E, ux not packed) that are pairwise non-adjacent with no soft pair among them, then
floor(|X|/2) <= OPT - |H| <= l. Proof: (x1,u,x2),(x3,u,x4),... are induced P3s with all pairs hard, pairwise
modification-disjoint and disjoint from H, so OPT >= |H| + floor(|X|/2). QED
Since each x has <= D soft partners, any independent set of N_h(u) (in G) of size m contains such an X with
|X| >= m/(D+1). So alpha(G[N_h(u)]) <= (D+1)(2l+1).
Meaning: the "hub of deletions" gadget (u adjacent to partners that must be thinned out, Vertex Cover on the
non-adjacency graph of N(u)) can have at most (D+1)(2l+1) = O(w l) leaves in a YES instance. Not just FPT-like:
it is bounded by f(l,w) before any algorithm runs. (Note it needs only D and l, not tw.)

### C. SM (small-merge lemma)  PROVED, VERIFIED-BY-COMPUTER (5683 solutions, n<=9)
Let S be a solution with |U| <= l', Ds = actual max soft degree. Let C be a cluster of G xor S and A, B two distinct
pieces inside C. Then:
 (i)  every pair {a,b}, a in A, b in B, is soft or lies in I (inserted hard pair);
 (ii) min(|A|,|B|) * (max(|A|,|B|) - Ds) <= |I cap (A x B)| <= |U|; hence max(|A|,|B|) <= Ds + |U| <= D + l;
 (iii) if C contains >= 2 pieces then |C| <= max(2D + l, 2l) (sum over piece pairs of (i)).
Proof. (i) A pair a-b with a, b in one cluster is an edge of G xor S. If it were a hard edge of G not in F it would
join a and b in one piece; hard edges in F are deleted; so it is a hard non-edge that was inserted, or a soft pair.
(ii) The soft pairs between A and B number <= Ds * min(|A|,|B|) (each vertex has <= Ds soft partners), the rest of the
|A||B| pairs are inserted by (i). If max >= D + l + 1 then min >= 1 gives LHS >= l+1 > l. (iii) Let C have pieces of
sizes s_i <= M = Ds + |U|, total s. Cross pairs = (s^2 - sum s_i^2)/2 <= Ds s/2 + |U| (soft cross pairs <= Ds s /2),
and sum s_i^2 <= M s, so s (s - M - Ds) <= 2|U|; if s > M + Ds this gives s <= 2|U|. QED
Test (code/w2b_sm.py): (i),(ii),(iii) asserted on every solution of cost <= |H|+2 on 373 random near-cluster
instances n<=9 (7360 multi-piece clusters, pieces up to size 6); no violation.

Reading of SM. Hard structure is rigid. Only <= l hard edges are cut. Two pieces share a cluster only if they are
small (<= D + l) and ALL cross pairs are soft except <= l inserted ones. So every merge is a clique in W plus l
extra pairs, and every large piece is a whole cluster by itself. The only large-scale relation the hard graph can
transmit is "same piece" (bounded size) and the only relation that connects different pieces is W-adjacency.

## 3. Why the naive Multicolored-Clique encodings collapse (precise)
Encoding E1 "hub of insertions". Hub u, classes V_1..V_k, pointer i inserts u v_{i,j_i}; compatibility = adjacency
v_{i,j}~v_{i',j'}. Selected partners lie in one cluster C. By SM(i) the pieces B_i of v_{i,j_i} and B_{i'} of
v_{i',j_{i'}} (both in C) have all cross pairs soft or inserted. Adjacency by a HARD edge puts them in the same piece,
so pieces are blobs of size <= D + l: the compatibility "graph" realisable by hard edges has degree <= D + l.
Adjacency between different pieces is only available as soft pairs (all of them), so the compatibility relation
that can be encoded is a subgraph of the piece graph of W, which has treewidth <= 3(w+1) (F1) and clique number
<= 3w+3. Clique detection in such a graph is FPT (trivial DP). The remaining case |B_i||B_{i'}| <= l has no
soft pair at all, but then the pair is "all inserted": compatible for free, no relation encoded.
In the CE language, the closure problem is: v's hard neighbour z not adjacent to u gives P3 {u,v,z}; z must then be
in C (insert uz) or the hard edge vz must be cut. Both cost excess per z, which is why N(u) and N(v) must agree up to
O(D + l) (Lemma 2 of fpt_attempt) and, by SM, why the partner piece is a bounded blob.
Encoding E2 "hub of deletions". Partners x_1..x_m hard-adjacent to u; keep a clique, delete the rest. Needs m large
independent sets inside N_h(u). By HUB the independent set size is <= (D+1)(2l+1): bounded. FAILS for the
reason the previous agent gave (Vertex Cover) and a stronger one (it cannot even be large).
Encoding E3 "tile vertex" (Grid Tiling): tile t_{ab} with hard edges to row vertex x_a and column vertex y_b.
All tiles + rows + columns become one hard component whose induced P3s x_a - t - y_b (x_a y_b non-edge, hard) are
all-hard; the family {(x_a,t_{ab},y_b)} is a modification-disjoint all-hard packing of size n^2 (pairs distinct), so
OPT >= n^2 by MAX/Lemma 1. FAILS at once. Making the pairs soft instead forces each x_a into n P3s, c = n, w >= n-1.
Encoding E4 "slip chain": a chain of cliques whose parity flips once; slip position = index j in [n]. Compatibility
across two parallel chains by sweeping gives only order/equality relations of the two positions (bounded-state
bags carry 1 bit "already slipped"). Grid Tiling with <= constraints would then be representable, but the
cell relation S_ij subset [n]^2 is arbitrary. To read it one needs a gadget touching both coordinate positions for
each of up to n^2 tiles, i.e. a subdivided K_{n,n} inside I(H) (treewidth >= n) or inside the hard graph (killed by
SM / E3). FAILS: arbitrary binary relations on [n] need a source of unbounded treewidth, and SM + (F1) say
the only candidate sources are W (tw <= 3w+2) and blobs (size <= D+l).
Escape attempts with mixed gadgets (insertions + deletions + packed edits): the hard part of such a gadget is
always either a blob (bounded) or cut by <= l deletions (SM: cut pieces are cliques up to soft pairs and <= l
inserted pairs, so big pieces are near cliques and cannot be split cheaply except at <= l vertices). I could not
construct any mixed gadget that carries information on [n] across the tree decomposition. This is evidence, not a proof of FPT.

## 4. What is the pointer after MAX + SM? (IDEA toward the positive side; labelled so)
After MAX assume H maximal (w' <= w+l). Then:
* every induced P3 of G has a soft pair; HE holds;
* for r outside V(H) the vertices of R are true-twin cliques (Lemma 4 of fpt_attempt) with only soft partners to V(H);
* by SM a pointer insertion uv puts u and v in a cluster whose pieces are all W-adjacent or fully inserted; the number of
  candidate partner pieces B for a piece A with |A||B| > l is <= |A| D <= (D+l) D, because B must contain a
  soft partner of a vertex of A (cross pairs soft or inserted, at most l inserted). Pointers between two
  pieces with |A||B| <= l need no soft pair, but by the witness argument below they are still tied to a third piece.
  Witness argument (PROVED, trivial): for an inserted pair uv in an optimal S, the cluster C has a vertex z != u,v
  (otherwise S minus uv is cheaper). With H maximal, for every z in C \ {u,v}: uz or vz is soft or inserted
  (else (u,z,v) with hard edges uz, vz and hard non-edge uv is an all-hard induced P3 of G). So every inserted pair
  has a bounded-size "support" made of soft pairs and other inserted pairs.
* A bounded candidate list for pointers (size f(l,w) per location rather than n^2) would give FPT by guessing the
  pointers in the tree decomposition (location = bag) and running Thm A inside. I did NOT prove such a list: the tiny
  case |A||B| <= l with two vertices that have a common soft neighbour z (candidates: pairs of vertices at W-distance 2,
  O(n D^2) of them, still linear in n) is where an argument is missing. So the problem remains OPEN.
Overall: SM shows that identity information in [n] cannot be moved by the hard graph, and F1 shows the soft graph
has bounded treewidth; a W[1]-hardness proof would need to move [n]-valued information through W using only bags of
size 3(w+1) and k selected positions. I now believe that is impossible (FPT), but the missing step is a locality
statement for the tiny-tiny pointers, not hardness.

## 5. Target: XP-tightness n^{Omega(l)} for fixed w, or ETH bounds
* No reduction found. If FPT in (l, w) is true then n^{Omega(l)} is false (under ETH), so this target would
  contradict my belief; SM is the reason.
* What is already available and sound: (a) 2^{o(l)} n^{O(1)} impossible under ETH even for w = 0 and c = 1 (take
  H empty... this is Cluster Editing, l = k; or a maximal vertex-disjoint packing). Cite Komusiewicz-Uhlmann 2012
  (no 2^{o(k)}, even bounded degree) -- UNVERIFIED re-read of the exact statement; literature_2.md recorded it as
  OPENED. (b) 2^{o(w)} for l = 0 (hardness.md A1, size accounting to be rechecked).
* The two combine into a "(2^{o(l)} and 2^{o(w)})" statement only; nothing multiplicative like n^{l}.

## 6. Other weakenings considered
* l + w + c: same parameter (point 3 of Sec. 0).
* l + tw(G) (instead of w): IDEA, no result. SM and MAX use only D = soft degree and the l budget; they do not
  need I(H) to have small treewidth. They give a partial statement: for FIXED D (that is, c) the hard graph is
  rigid, so any hardness for bounded c (including c = 2, the authors' conjecture) cannot come from hub or tile
  encodings either: with D = 4 the same argument kills E1-E3. The authors' "grid-like structures" for c = 2 must
  then live in W itself (max degree 4, unbounded treewidth): a hardness proof for c = 2 would have to move [n]-values
  along a planar-like W while SM forbids using the hard graph as a channel. This is the only realistic place left
  for W[1]-hardness: l + c (unbounded w), not l + w.
* Edge-disjoint (instead of modification-disjoint) packings: not examined.

## 7. Summary table
| item | status |
|---|---|
| c <= w+1, "+c" weakening is vacuous | PROVED (trivial) |
| tw(W) <= 3(w+1)-1, D <= 2(w+1) | PROVED (F0, F1) |
| MAX: H maximal at cost w' <= w+l; hard-edge lemma HE | PROVED, VERIFIED-BY-COMPUTER |
| HUB: independent hard neighbours of a vertex <= (D+1)(2l+1) | PROVED, VERIFIED-BY-COMPUTER |
| SM: merged pieces small, cross pairs soft or inserted | PROVED, VERIFIED-BY-COMPUTER |
| E1-E4 hub, deletion, tile, slip-chain encodings | FAIL (Sec. 3, reasons) |
| W[1]-hardness of (l,w) | NOT FOUND; evidence against |
| FPT in (l,w) | OPEN; reduction to a locality claim for tiny-tiny pointers (Sec. 4) |
| n^{Omega(l)} for fixed w | NOT FOUND |
| hardness for l + c (c = 2) | OPEN; only possible via W (Sec. 6) |

## 8. Suggested next steps
1. Prove or refute the locality claim of Sec. 4: with H maximal, there is an optimal solution in which every inserted
   pair has both endpoints within W-distance <= f(l,w) of V(H) or in a piece touched by a soft pair. Search for a
   counterexample with n = 9..11 using the structure of SM (many singleton pieces, a common soft neighbour).
2. If it holds, write the FPT algorithm: MAX, then guess pointers by location (bounded candidates), then Thm A.
3. For the paper's lower-bound section, use SM as the explanation of why the l + w hardness is not available and
   state c = 2 hardness (l + c) as the open problem.
