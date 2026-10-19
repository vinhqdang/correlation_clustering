# Wave 5, prover A: locality of tiny-tiny pointers

Code: code/w5a_lib.py (anatomy checker), w5a_exh.py (exhaustive), w5a_ilp.py (HiGHS, n = 11..13), w5a_loc.py (localizability),
w5a_far.py, w5a_show.py / w5a_show2.py (examples). Reuses code/a4_lib.py (partitions, pieces, packings). No git.
Labels: PROVED / VERIFIED-BY-COMPUTER (evidence only) / FAILS / CONJECTURE / IDEA.

## 0. Verdict
1. The assigned gap is CLOSED: tiny-tiny pointers are W-local. Precisely (Cor. 2 below, PROVED): if two pieces A, B of an optimal
   solution share a cluster and have NO packed pair between them, then every cross pair (a,b) has a common soft neighbour z in the
   cluster. So a and b lie in two distinct members of H that share z, i.e. adjacent nodes of I(H), i.e. in the vertex set of ONE bag.
   Candidate partners for a: <= D(D-1) <= 4c(c-1), independent of l and n. No pieces-at-distance and no profile matching needed.
   The hypothesis |A||B| <= l is not needed: "no packed pair between A and B" is the right condition (and |A||B| > l forces a packed pair, SM).
2. The proof is not an exchange of B against canonical candidates. It is the elementary fact that in an optimal solution every
   vertex has at least half of its cluster as G-neighbours, hence clusters have diameter <= 2 in G. Combined with pair-maximality
   (MAX/HE) this forces the middle vertex z to be soft on at least one side, and "no soft pair between A and B" forces both sides.
3. This does NOT give FPT in (l,w). The same argument classifies EVERY inserted pair into two types (Thm 1): T1 = common soft neighbour
   (bag-local, no identity problem), T2 = "hard step then soft step" (a --hard, same piece-- z --soft-- b). T2 pointers are NOT local:
   VERIFIED forced in about 1% of random n = 8 instances (all optima), and at infinite I(H)-distance (different components, c = 1).
   For T2 the partner is determined by a local anchor up to <= l choices (Lemma 3, PROVED), but the far end sits at hard distance 1
   from the anchor, so the identity problem moves to a hard-edge "tunnel". Whether tunnels can be handled is open (Sec. 5, 6).
4. Hardness hint (IDEA, not verified): tunnels transport equality of [n]-colours (hard cliques are near-equivalence classes by HE),
   and equality between adjacent grid cells is what Grid Tiling needs. Sec. 6.

## 1. Setting
H pair-maximal (no induced P3 of G has all three pairs unpacked; e.g. H u M from MAX). Soft = packed pair, hard = other pair.
Q(v) = soft partners of v (|Q(v)| <= D = 2c). S an OPTIMAL solution (|S| = OPT), C a cluster of G xor S, U = S minus packed pairs,
u = |U|. Inserted pair = hard non-edge inside a cluster (in U). Pieces = components of kept hard edges (hard edges of G inside a
cluster). A hard edge of G inside C is kept, hence its endpoints lie in one piece.

## 2. Lemma 1 (half-degree and diameter)  PROVED
(a) For every vertex v of an optimal solution, d_C(v) := |N_G(v) cap C| >= (|C|-1)/2.
(b) If a, b in C, ab not in E(G), then some z in C \ {a,b} has az, zb in E(G).
Proof. (a) Moving v into a new singleton cluster turns its d_C(v) kept edges into deletions and its |C|-1-d_C(v) inserted non-edges
into free pairs; other pairs are unchanged. Cost change = d_C(v) - (|C|-1-d_C(v)) >= 0 by optimality. (b) N(a) cap C and N(b) cap C
lie in C \ {a,b}, which has |C|-2 elements, while d_C(a) + d_C(b) >= |C|-1; so they meet. QED
Uses neither H nor l. It needs optimality (or just local optimality for single-vertex moves): see the control in Sec. 4.

## 3. The anatomy of an inserted pair
### Theorem 1 (T1/T2 dichotomy)  PROVED
Let (a,b) be an inserted pair of an optimal S (a, b in C, ab hard non-edge), H pair-maximal. Then there is z in C \ {a,b}, az, zb in E, and
 (T1) az and zb are both soft; then they lie in two DISTINCT members P_a, P_b of H with z in P_a cap P_b, or
 (T2) exactly one is hard, say az; then az is kept, so a and z lie in the same piece, and zb is soft (symmetric: swap a,b).
Proof. z from Lemma 1(b). If az and zb were both hard, (a,z,b) is an induced P3 with all three pairs hard (ab hard by assumption),
contradicting pair-maximality. If both soft: pairs are packed in members P_a (of az) and P_b (of zb); if P_a = P_b then a, b in one member
and ab is packed, contradiction; so P_a != P_b and both contain z. If az is hard it is a hard edge of G inside C, hence kept. QED
### Corollary 2 (tiny-tiny locality; the assigned statement)  PROVED
Let A != B be pieces inside one cluster C of an optimal S with no soft pair between A and B (in particular every tiny-tiny merge with no packed pair between
the pieces; SM gives |A||B| = number of cross pairs inserted <= l in that case). Then for every (a,b) in A x B
(i) ab is an inserted pair; (ii) there is z in C \ (A u B) with az and zb soft, in distinct members P_a != P_b, both containing z.
Proof. (i) A cross pair is hard, and a hard edge would be kept and merge the pieces. So ab is a hard non-edge in one cluster: inserted.
(ii) Theorem 1 gives z. Type T2 is impossible: say az hard, kept, so z in A; then zb in E is a pair between A and B: if hard it is kept
(z, b same piece, contradiction), if soft it is a soft pair between A and B (excluded). So T1. If z were in A, then zb soft is again a soft pair
between A and B; same for z in B. QED
Consequences (all PROVED from the corollary):
 (a) Bag locality. Let (T, X) be a tree decomposition of I(H) of width w. P_a, P_b share z, so are adjacent in I(H) and some bag t contains
     both; so a, b lie in V(X_t) = union of the vertex sets of the members in X_t, |V(X_t)| <= 3(w+1). More: clique(a), clique(b), clique(z) (members
     through a, b, z) lie in the closed I(H)-neighbourhood N[X_t]; (T, N[X]) is a tree decomposition of I(H)^2 of width <= (w+1)(3w+1)-1
     (edge of the square via a common neighbour m: both ends in N[bag of m]; a vertex v is in the bags of N[v], whose subtrees pairwise meet, so
     connected; max degree of I(H) is <= 3(c-1) <= 3w). So a tiny-tiny pointer is a pair inside one bag of a width-g(w) decomposition of I(H)^2.
 (b) Candidate list. Given a, the partners b are in Q(z) for z in Q(a): at most D(D-1) <= 4c(c-1) candidates, independent of n, l.
     Globally, candidate pairs <= n * C(D,2), each attached to a bag.
 (c) The whole blob A u B is local: for fixed b in B every a' in A is at W-distance exactly 2 from b, so A u B lies in the W-ball of radius 4 around b.
 (d) Non-identity. Because z in C is soft to both ends, the relation "a and b can be merged" is a relation inside one bag. No identity matching
     between far-apart bags occurs for these pointers; the "profile" of the partner is read off inside the bag.
### Lemma 3 (near-twin agreement; anchored candidate bound)  PROVED
Let z, b be in one cluster of any solution S (not nec. optimal) with u = |U|. Then
  |{a not in {z,b} : az, ab both hard, [az in E] != [ab in E]}| <= u.
In particular (anchor form) for a soft pair (z,b) in one cluster, M(z,b) := N_h(z) \ (N[b] u Q(b)) has |M(z,b)| <= u <= l, and M is computed from G alone.
Proof. For such a: if a in C, one of az, ab is a hard non-edge inside C, so inserted; if a not in C, the hard edge az or ab joins different clusters, so is deleted.
The pair is in U; distinct a give distinct pairs (a is an endpoint and the other endpoint is z or b, and a = z, b is excluded). QED
This is the "profile of B coincides with that of A modulo l" alternative asked for: any two vertices that share a cluster have hard neighbourhoods
(outside Q(z) u Q(b)) differing in at most l vertices. For T2 it gives: a in M(z,b), so the far endpoint is one of <= l vertices fixed by the anchor.

## 4. Computer evidence (VERIFIED-BY-COMPUTER, evidence only)
Checker anatomy() tests, for one partition: (D) every inserted pair has z in C with az, zb in E; (H) not all-hard; (T2) hard side lies in the same piece
and the other side is soft; (Cor. 2) for every pair of pieces in a cluster without soft pair between them: all cross pairs inserted and each has a
soft common neighbour z in C \ (A u B); (L3) the bound of Lemma 3 for every pair z, b in a cluster.
* Exhaustive over all set partitions, n = 6..9, ~1,900 pair-maximal instances (gnp, noisy clusters, hub, sparse, Fano/AG gadgets, +noise, random/overlap
  maximal packings): 7,082 OPTIMAL partitions, 0 violations; 1,500+ T1 pairs, 1,000+ T2 pairs, ~470 tiny-tiny cluster pairs (no soft pair between pieces).
* Control: on partitions of cost OPT+1, OPT+2 (about 110,000) (D) and (Cor. 2) fail thousands of times (n=8: 11,231 and 11,043 times); only (L3) (pure
  counting) never fails. So optimality (Lemma 1) is used and necessary; this is not a vacuous check.
* ILP (HiGHS, 7 optimal partitions per instance by random side objective with cost <= OPT), n = 11, 12, 13: path, cycle, path plus chords, caterpillar, two
  paths with cross edges, gnp, noisy, hub, sparse, Fano; plus the P13 of wave2_search (H = {123, 456, 789, 10 11 12}): ~5,700 optimal partitions,
  0 violations (cycles and caterpillars have 40-100 T1 pairs and 20-100 tiny-tiny cluster pairs per family; P13 has only a T2 pair).
* Necessity of T2 (localizability, w5a_loc.py): minimise, over optima, the number of inserted pairs with both ends in V(H) and no common soft neighbour (not
  W^2-local). n = 8, 5 families x 600 seeds: some instance forces such a pair in 25 of 2,944 (4/581 noisy, 9/597 gnp, 7/600 hub, 1/566 sparse, 4/600 Fano);
  n = 9: 2/149. Smallest example (n = 8, unique optimum): edges 02 06 14 16 17 23 25 37 45 46 47 57, H = {(0,6,4),(3,7,4),(2,5,4),(0,2,3),(6,1,7)};
  cluster {1,4,5,7} = K4 minus the inserted pair 15; common neighbours 4, 7 both give hard+soft (14 hard, 45 soft; 17 soft, 57 hard).
  So "every inserted HH pair can be taken W-local" FAILS; Theorem 1's T2 branch is real.
* Far T2 / far insertions (w5a_far.py): n = 9, H = three vertex-disjoint P3s (c = 1, I(H) edgeless) + random hard edges: 4 of 423 instances have an inserted
  pair between different P3s in EVERY optimum (e.g. 5-8 in edges 01 12 13 34 45 47 48 57 67 78). So T2 reaches infinite I(H)-distance (consistent with
  w2a_far.py for deletions). For c = 1 this is harmless only because BFK gives |V(H)| = O(l).

## 5. What is and is not established for FPT in (l,w)
Inserted pairs of an optimal solution split into
 T1 (bag-local; <= D(D-1) partners; this includes the whole tiny-tiny case): can be guessed locally in a DP over a width-g(w) decomposition of I(H)^2.
   Their interaction with the CSP of Theorem A touches only cliques in N[X_t]; their interaction with the C0 condition (triples {a,b,y}, y outside the
   environments) is a unary test on the pair (Lemma 2 of fpt_attempt) and, for two T1 pairs in a common triple, local again (pair-pair distance <= 4).
   CONJECTURE (not proved here): T1 pointers, handled as bag-local flags with a counter <= l, cost only f(l,w)-factors.
 T2 (anchored): anchor (z,b) is a soft pair (local); the partner a is in M(z,b), |M| <= l, computed from G (PROVED, Lemma 3); the far end a lies at hard
   distance one from z. Effect at a's end on the CSP: for y in Q(a) the final status of the soft pair (a,y) must equal G(b,y) (triple {a,b,y}), while the kept
   hard edge az already forces final(a,y) = G(z,y). Hence IF no other U-pair touches a, b, z or Q(a), the pointer adds NO new constraint at a's end beyond
   the base CSP, and the compatibility test is G(b,y) = G(z,y) for y in Q(a): a local-to-a test of a hard-neighbour pair. This is an observation, not a proof
   that tunnels are harmless: with several pointers the far ends can chain. CONJECTURE (T2-reduction): T2 pointers can be absorbed into the base CSP by a
   precomputed per-(z,b) unary constraint plus <= l extra pairs, again not verified.
 Deletions of hard edges between far members (n = 7 example of wave2_c1_w1 3.4(iv)): untouched by this note; they are tunnels too (kept hard edges are the same objects).
So the claim "FPT in (l,w)" is NOT proved. What is proved: the tiny-tiny pointer is not an identity problem; the whole identity problem sits in T2 and deleted-hard-edge
tunnels. The wave2 Sec. 4 worry (n D^2 candidates "still linear in n") is resolved in the sense that those candidates are bag-local; the remaining linear-in-n family
is the family of anchors, which a DP over bags handles, and the far ends, which it does not.

## 6. Suggestion about hardness (IDEA, unverified)
Any W[1]-hardness proof for (l,w) must now use T2 or deletion tunnels, since T1 is local (Cor. 2) and R-pointers reduce to HH pointers (wave2_c1_w1 1.4).
A tunnel is a kept hard edge az with a far; hard edges are transitive up to soft pairs and l exceptions (HE), so hard components are near-cliques of bounded
size (SM): tunnels realise EQUALITY of colours between far positions (blob {s_a, t_a} for colour a), not arbitrary adjacency relations. Grid Tiling needs exactly
equalities between neighbouring cells (a_ij = a_{i+1,j}, b_ij = b_{i,j+1}) plus arbitrary local membership (a_ij, b_ij) in S_ij, which the width-w CSP part can
supply with n^2 local anchors per cell. A reduction with l = O(k^2) pointers and w = O(k) is therefore not excluded by anything proved so far, and T2 is its only
candidate channel. Obstacles I could not resolve: (i) an inactive anchor must cost nothing, (ii) M(z,b) <= l forces each anchor's far end to be unique, so the
colour must be stored in WHICH anchor is active, and MAX/HUB+ (alpha(N_h) <= 2l + c + 1) must not be violated by the n anchors sharing hard neighbours. Not tested.
Suggested next experiment (not done): an ILP search at n = 12..14 for two-region instances in which every optimum activates exactly one T2 pointer whose far end
is determined by its anchor and the far-end bag rejects mismatching activations (that would show the far-end check is real).

## 7. Table
| item | status |
|---|---|
| Lemma 1: half-degree, diameter <= 2 in optimal clusters | PROVED, VERIFIED |
| Theorem 1: inserted pair is T1 (two soft steps) or T2 (hard-in-piece + soft) | PROVED, VERIFIED (7,082 + ~5,700 optimal partitions) |
| Cor. 2: no soft pair between pieces => all cross pairs have soft common neighbour; bag-local; <= D(D-1) partners; A u B in W-ball radius 4 | PROVED, VERIFIED |
| Lemma 3: near-twin agreement, |M(z,b)| <= l | PROVED, VERIFIED |
| all inserted HH pairs W-local in some optimum | FAILS (25/2,944 at n = 8; smallest example in Sec. 4) |
| T2 reaches different I(H)-components | VERIFIED (4/423, c = 1, n = 9) |
| FPT in (l,w) | OPEN; remaining difficulty = T2 and deletion tunnels |
| hardness via tunnels | IDEA only (Sec. 6) |
