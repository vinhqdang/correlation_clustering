# Wave 6 (prover B): the two remaining pointer types, T2 inserted pairs (R1) and deleted hard edges (R2)

Labels: PROVED (proof written here), VERIFIED-BY-COMPUTER (finite check, evidence only), CONJECTURE, REFUTED, IDEA.
Code (own, independent; imports only numpy and `w6b_lib.py`, nothing from earlier waves): `code/w6b_lib.py` (partitions, instances,
Delta, I(H) distances), `w6b_conflict.py`, `w6b_bad.py`, `w6b_prune.py`, `w6b_classify.py`, `w6b_gm.py`, `w6b_family.py`,
`w6b_far.py`, `w6b_nf.py`, `w6b_merge.py`. Nothing in `paper/` was touched, no git.

Notation as in packing.tex / locality.tex / dproute.tex. H pair-maximal, soft = packed, hard = otherwise, Q(v) soft partners,
N_h(v) hard neighbours, D = 2c, S a solution = a partition pi, x(p,q) = 1 iff p,q in one cluster, U = edited hard pairs,
|U| <= |S| - |H|. New notation:
* Delta(z,b) = { y not in {z,b} : zy and by hard, exactly one of them an edge of G }.   Gamma_h(z,b) = N_h(z) cap N_h(b).
* For a soft pair pi = (z,b) in a solution: weight c(pi) = |Delta(z,b)| if z,b are in one cluster ("merged"), |Gamma_h(z,b)| if not.
  A *conflict* is a pair (pi, y) with y in the relevant set. Both sets are computed from (G,H) alone.
* A member P of H is *locally bad* if each of its three states of member-cost 1 (delete xy, delete yz, insert xz) has some soft pair
  of P with weight > 0 in that state. Static.
* An edit e = (a,b) in U is *two-sided* if there are witnesses t, t' with t soft to a and t' soft to b, where a witness is a vertex t such that
  {a,b,t} is an induced P3 of G xor (S minus e). *Far* = the members of a and b are at I(H)-distance >= K (infinite allowed).

## 0. Verdict

1. FPT in (l,w) is NOT proved and not refuted. Nothing below is a proof of FPT in (l,w), and the hardness route is not completed either.
2. What is proved and new: a **conflict lemma** (W1) that makes the bounded-candidate statement sharper and gives *static weights*; a
   **defect bound** (W2: at most 2Dl conflicts); a **bad-member bound** (at most 2cl members have no conflict-free cost-1 state); a
   **rigidity lemma** for hard edges; a **pruned version of the FPT theorem of dproute.tex** (Theorem P: FPT in (l,D,tw(Gamma^l)) where
   Gamma^l is Gamma' minus the soft pairs whose merge would cost more than l; strictly stronger than the omega theorem); and an
   **exact obstruction family** (Proposition FAR): for every m>=3 an instance with I(H) a path, c=2, w=1, l=1 in which every optimum
   deletes a hard edge whose ends lie in members at I(H)-distance m-1.
3. Answer to the question "can T2 pointers and deleted hard edges be encoded as bag-local flags plus l-bounded choices?":
   * the candidate lists, YES (PROVED): from either end the other end lies in an explicit set of size <= (D+l)(D+1), computed from (G,H);
   * the flags themselves, YES locally: each end can carry "which of <= D+l hard neighbours / <= l far candidates is edited";
   * what cannot be done bag-locally is the **agreement of the two flags of one edited pair when both ends are active (two-sided) and the
     two members are far in I(H)**. This is the exact obstruction (Section 5). It is not a technicality: Proposition FAR forces it for
     l=1, w=1, and n=10 instances show it also occurs at members that are not locally bad (so "put the bad members in every bag" does not
     remove it: the conjecture NF-bad is REFUTED).
   * the **far end of a T2 insertion carries no status variable** (Lemma FE: it is transparent), so R1 is anchor-local EXCEPT for (i) the y in Delta(z,b) that are repaired by
     deletion (those are R2) and (ii) double counting when the far end is itself an anchor of the same pair (a two-sided T2 edit; Section 5.3). R1 does not disappear into R2, it
     shares R2's obstruction.
4. The identity problem is therefore concentrated in one object: a static graph of *potentially two-sided* hard pairs, which has
   unbounded treewidth in general (hard edges are free). Whether an instance-level reduction (the missing lemma ML) removes it is exactly
   the open question; I could not decide it.

## 1. Proved static tools

### W1 (conflict lemma)  PROVED, VERIFIED-BY-COMPUTER
Let S be ANY solution and pi = (z,b) a soft pair.
(a) If z,b lie in one cluster, every y in Delta(z,b) has zy in U or by in U.
(b) If z,b lie in different clusters, every y in Gamma_h(z,b) has zy in U or by in U.
Distinct y give distinct edited pairs (the pair contains y and z or b). Hence |U| >= c(pi) for every soft pair pi.
Proof. (a) In a cluster graph z and b have the same closed neighbourhood-in-clusters: x(z,y) = x(b,y) for all y. For y in Delta(z,b) the G-values
of the two hard pairs differ, so at least one of them is changed by S. (b) If y is a hard neighbour of both and neither zy nor by is edited,
then x(z,y) = x(b,y) = 1, so z,b,y are in one cluster, contradiction. Hard pairs that are edited are in U by definition. QED
Check: `w6b_conflict.py`, all partitions of cost <= |H|+3 of random pair-maximal instances n = 5..9 (132,827 partitions, about 1.1e5
conflicts), 0 violations.
Corollary (the near-twin lemma of locality.tex): M(z,b) = N_h(z) minus (N[b] cup Q(b)) is contained in Delta(z,b), so |M(z,b)| <= l.

### W2 (defect bound)  PROVED, VERIFIED-BY-COMPUTER (tight)
For every solution S: sum over soft pairs pi of c(pi) <= 2 D |U|.
Proof. Map each conflict (pi,y), pi = (z,b), to an edit of U in {zy, by} (W1). An edit e = (p,y) is chosen by conflicts (pi,y) with pi containing p
(at most |Q(p)| <= D soft pairs) and by conflicts (pi',p) with pi' containing y (at most D). QED
The constant 2D is attained (ratio exactly 1.000 observed at n = 6,7,8), so "sum <= D|U|" is false.
Consequence (Corollary LB, PROVED): cost(S) >= sum_i c_i + ceil( W(sigma) / 2D ), where sigma is the soft assignment of S and W(sigma) the total weight (c_i = number of edited
pairs of member i, a function of sigma). This is a sum of bag-local terms plus a counter bounded by 2Dl, so the minimum over assignments that are transitive on all-soft triples
is computable by a width-w dynamic program in f(l,w) n^O(1); it is a valid NO-certificate. It is only a relaxation: it ignores the triples with two soft pairs and one hard pair (those
can be added as conflicts with a single repair edit, weight <= D per edit), which I did not carry out, so I do not claim that it is exact at l = 0.

### F1 (static pruning of merges)  PROVED, VERIFIED-BY-COMPUTER
If |Delta(u,v)| > l, the soft pair (u,v) is in different clusters in every solution with |S| <= |H|+l. If |Gamma_h(u,v)| > l it is in the same
cluster in every such solution. Proof: W1 and |U| <= |S|-|H| <= l. Check: `w6b_prune.py` part 1 (49,119 partitions with l = 0..3), 0 violations.

### DEL (rigid hard edges)  PROVED, VERIFIED-BY-COMPUTER
If zy is a deleted hard edge of S, then |U| >= 1 + |N_h(z) cap N_h(y)|. Hence in a solution with |S| <= |H|+l a hard edge with at least l common
hard neighbours is never deleted ("rigid"), and the hard components of the rigid edges lie inside single clusters.
Proof. For each common hard neighbour w, z and y are in different clusters, so w cannot be kept with both: zw or yw is deleted; distinct w give distinct pairs,
different from zy. QED. Check: same run as W1, 0 violations.

### BAD (few locally bad members)  PROVED, VERIFIED-BY-COMPUTER
For every solution S with |S| = |H| + eps, the number of locally bad members is at most 2c*eps. In particular more than 2cl bad members is a NO instance
(static, polynomial test).
Proof. Cost = sum_i c_i + |U| (Lemma CP(b)), eps = sum_i (c_i - 1) + |U|. A bad member with c_i >= 2 contributes >= 1 to the first sum. A bad member with
c_i = 1 is in one of its three cost-1 states, which has a conflict (pi,y); by W1 an edit of U incident to V(P_i) exists; an edit has two ends, each in <= c members, so
at most 2c|U| such members. Total <= sum(c_i - 1) + 2c|U| <= 2c eps. QED. Check: `w6b_bad.py`, 169,254 partitions of cost <= |H|+3, 0 violations.

### FE (the far end of a T2 insertion is transparent)  PROVED
In any cluster graph, a kept hard edge az gives x(a,t) = x(z,t) for all t. So if e = (a,b) is a T2 pair (az hard and kept, zb soft and merged), the status
x(a,b) = 1 is implied by the two statuses x(a,z) = x(z,b) = 1; no variable is attached to a. The pair costs 1 and is one of the pairs of Delta(z,b); the cost of
merging the soft pair (z,b) is the number of hard non-edges it creates inside the cluster, which is c(z,b) up to the sharing discussed in Section 5.
What is NOT transparent: a y in Delta(z,b) that is repaired by DELETING zy (y leaves the cluster of z). That is a deleted hard edge, i.e. R2.

### AN (anchoring of edits in optima)  weak form PROVED, strong form VERIFIED-BY-COMPUTER
Let S be optimal and e = (a,b) in S. Then e is critical: some t makes {a,b,t} an induced P3 of G xor (S minus e). Moreover either {a,b,t} contains a packed pair,
or it contains another edited pair of S.
Proof. S minus e is not a solution (else smaller), and G xor S is a cluster graph, so every induced P3 of G xor (S minus e) contains the pair e: there is a witness t.
If all three pairs of {a,b,t} are hard and e is the only edited one, G restricted to the triple equals (G xor (S minus e)) restricted to the triple, an induced P3 with
all pairs hard, contradicting pair-maximality. QED.
Strong form (no edited pair without a soft-anchored witness) was never violated: `w6b_classify.py`, optimal partitions of random pair-maximal instances n = 5..10
(about 7,300 instances), 0 violations, witnesses always exist. Not proved.

## 2. Edit anatomy (VERIFIED-BY-COMPUTER, `w6b_classify.py`, n = 5..10)
Counts of edited hard pairs in optimal partitions (type; one-sided or two-sided; near = I(H)-distance <= 2, far = >= 3 or infinite), n = 9 / n = 10:
T1: all two-sided, all near (as proved in locality.tex). T2: one-sided far 449 / 159, one-sided near 288 / 83, two-sided far 3 / 3, two-sided near 204 / 73.
Deleted: one-sided far 850 / 292, one-sided near 417 / 99, two-sided far 80 / 8, two-sided near 1768 / 577.
So most far edits are one-sided; two-sided far edits are a small but unavoidable minority.

## 3. What is bag-local and what is not (answer to the posed question)

| edit | end carrying the anchor | other end | encoding |
|---|---|---|---|
| T1 inserted | both in V of one bag (locality.tex, Prop. bag) | | local |
| T2 inserted (a,b), anchor soft pair (z,b), az hard, one-sided | z,b in one bag; a in Delta(z,b), at most l candidates, static | a is transparent (FE) | anchor-local weight c(z,b) plus an l-bounded subset choice of Delta(z,b); no far flag |
| T2 inserted, two-sided (a is also an anchor end of the same pair) | both ends | | the pair is counted from both anchors: dedupe needs agreement of two far soft-pair states (5.3) |
| T2 inserted, y in Delta(z,b) repaired by deletion | z | y needs the flag "my hard edge to N_h(y)[i] is cut", i <= D+l | half-flags at both ends |
| deleted hard edge zy | z, flag index i <= D+l | y, flag index j <= D+l | half-flags at both ends |

* One-sided edit (witnesses only at one end): the passive end is indifferent to the status of the pair (by definition of a witness every triple through the passive end is a valid
  cluster pattern with the pair either way), so a single flag at the active end should suffice. This is an observation, not a theorem: I did not prove soundness of a
  flag-DP in the presence of cascades (triples containing two edited hard pairs; see AN).
* Two-sided edit: both ends have a conflict with the pair, both flags must be present and agree. If the two members are in a common bag the agreement is checked there.
  If they are far, it is an equality between two flags in different bags. For T1 this never happens (locality.tex). For T2 and deletions it does (Section 5).

## 4. An FPT theorem for the case where the surviving tunnels have bounded treewidth

### Theorem P (pruned omega theorem)  PROVED (modification of Theorem omega of dproute.tex; the DP is the one of that proof)
Let H be pair-maximal, l the budget (apply Lemma MAX first and read everything for (G,H',l')). Let BIG, V', Z be as in dproute.tex and let
Gamma^l be the graph on V' with edge set: all hard edges of G inside V', and all soft pairs (u,v) inside V' with |Delta(u,v)| <= l.
Let omega_l = tw(Gamma^l) and let a tree decomposition of that width be given. Then OPT(G) <= |H|+l iff |E| + sum_K f(K) + m^* <= |H|+l, where m^* is computed in
time (Z(omega_l+1))^{O(omega_l)} n + n^{O(1)}.
Proof. Theorem omega with two changes. Define f~(C) = C(|C|,2) - 2 e^+(C) where e^+(C) counts only the edges of G that are edges of Gamma^l (so all hard edges and the
unpruned soft edges); the DP of the proof computes m^* = min over Gamma^l-connected partitions of V' with clusters of <= Z vertices of sum f~(C) (the introduce and join steps only use
edges that lie in a common bag, and every edge counted by e^+ is an edge of Gamma^l).
(<=) f(C) <= f~(C) for every C (e^+ <= e), so cost(pi) = |E| + sum_K f(K) + sum f(C) <= |E| + sum_K f(K) + sum f~(C).
(=>) Let S be optimal with |S| <= |H|+l. By dproute.tex (Lemma NF) the clusters of S in V' have <= Z vertices and are connected using only hard edges kept in the cluster and soft pairs
inside the cluster. A soft pair inside a cluster is merged, so by F1 it has |Delta| <= l: it is an edge of Gamma^l. So the clusters are Gamma^l-connected, and no pruned pair lies inside
a cluster, hence f~(C) = f(C) for each cluster C of S, and cost(S) >= |E| + sum_K f(K) + m^*. QED
Why this is strictly stronger than the omega theorem: Gamma^l is a subgraph of Gamma'. In the family of Proposition expander (dproute.tex), the smallest weight of a cross soft
pair is 5 for K_4 and for the Petersen graph (`w6b_prune.py` part 2), so for l <= 4 every cross soft pair is pruned and Gamma^l has components of at most 5 vertices (largest component 5 in
both cases, checked), while tw(Gamma') >= tw(F). For larger l one would take d-regular expanders with d about l/2 (weights grow linearly in d while the hard degrees stay below 3l+4, so BIG stays empty); this is a sketch, not run.
Limit: for fixed l the same family with bounded clique size has tw(Gamma^l) unbounded once l reaches the weight, and these instances are still trivial; Theorem P does not use optimality
of the merges, only their cost, so it cannot decide the open question. It is the case "ML-type obstructions do not appear" in the precise form: all pairs of far vertices that can ever
be in a common cluster at cost <= l have bounded treewidth.

## 5. The exact obstruction

### 5.1 Pointers need not be local in any power of I(H), even for l = 1, w = 1  (Proposition FAR)  PROVED + VERIFIED-BY-COMPUTER (m = 3,4,5)
For m >= 3 let G_m be the path 0-1-...-2m with the extra edge {1, 2m-1}, and H_m = { (2i, 2i+1, 2i+2) : 0 <= i < m }.
Claims. (i) H_m is a pair-maximal packing, I(H_m) is a path on m nodes (c = 2, w = 1). (ii) OPT(G_m) = m+1 = |H_m|+1. (iii) In EVERY optimum the hard edge e* = {1,2m-1} is
deleted. Its ends lie in M_1 and M_m, at I(H)-distance m-1. (iv) The two end members are the only locally bad members.
Proof. (i) Cycle 1..2m-1 has length >= 5, so the members are induced P3's with pairwise disjoint pair sets; the only hard edge is e*, an induced P3 has two edges, so some pair is soft.
(ii) cost(pi) = |E| + sum f(C), |E| = 2m+1, f(C) = C(|C|,2) - 2e(C). G_m is connected with 2m+1 vertices and 2m+1 edges (unicyclic), so an induced subgraph on s vertices has at most s-1
edges unless it contains the cycle (then at most s, and s >= 2m-1 >= 5). Hence f(C) >= -1 for all C, with equality only if |C| in {2,3} and G[C] is connected. Such clusters have
>= 2 vertices, so there are at most m of them and sum f >= -m, cost >= m+1. The partition {0,1},{2,3},...,{2m-2,2m-1},{2m} has cost m+1.
(iii) An optimum has sum f = -m, so exactly m clusters with f = -1 (each connected, of size 2 or 3) and all others f = 0; at most one vertex lies outside the f = -1 clusters.
A cluster with f = -1 containing the pendant 0 must contain its only neighbour 1, similarly for 2m and 2m-1. If 1 and 2m-1 shared a cluster C*: if C* has f = -1 and size 2, neither pendant is in an f = -1 cluster (a pendant
there needs its only neighbour in the same cluster), contradicting 'at most one outside'; if C* has f = -1 and size 3 it contains at most one pendant, the other pendant p has its only neighbour inside C*,
so p lies in no other f = -1 cluster and the remaining 2m-3 vertices cannot form the m-1 clusters of size >= 2 that are needed; if f(C*) >= 0, at most 2m-1 vertices remain for m clusters of size >= 2. So 1 and 2m-1 are in different clusters. QED. (iv) computed.
Search: `w6b_family.py` (path P_{2m+1} plus one extra edge, all choices, m = 3,4,5) finds exactly one extra edge, {1,2m-1}, for which every optimum edits a pair at I(H)-distance >= m-1.
Check: `w6b_gm.py` (all set partitions, n = 7, 9, 11): OPT = m+1, 7 / 9 / 11 optima, all delete {1,2m-1}, distance 2 / 3 / 4, the edit is two-sided in every optimum, the subset inequality f >= -1 holds.
Consequence. No bag-local encoding in which a pointer is a pair inside one bag of a decomposition of I(H), or of I(H)^r for r < m-1, can represent this optimum. The cure is to add the one
edge M_1 M_m to the decomposed graph; that raises treewidth by 1 here, but in general the number of such potentially active tunnels is not bounded (hard edges are free).
Side-additive count: a DP that charges the conflicts of each end separately pays m+2 (`w6b_gm.py`: min over partitions of cost + number of far two-sided edits = OPT+1 in all three cases), i.e.
it is incomplete. Sharing is necessary.

### 5.2 Locally bad members do not remove the obstruction  (conjecture NF-bad REFUTED, `w6b_nf.py`)
CONJECTURE NF-bad (tested): there is an optimum in which every far two-sided edit touches a locally bad member; then adding the <= 3*2cl vertices of the bad members to every bag
would make all far two-sided edits bag-local. REFUTED: n = 10, family "path of members + random extra edges": with far meaning I(H)-distance >= 2, 7 of 4000 instances; with far meaning >= 3, 3 of 8000 instances have, in EVERY optimum, a far two-sided edit with both ends at members that are NOT bad.
Example (one of them, unique optimum, OPT = |H|+1 = 6): edges 01 06 12 23 29 34 45 56 67 78 79, H = (0,1,2) (2,3,4) (4,5,6) (6,7,8) (2,9,7). The only hard edge is {0,6}; the unique optimum is
{0,1},{2,9},{3,4},{5,6},{7,8}; it deletes {0,6}, witnesses 1 (soft to 0) and 5 (soft to 6), I(H)-distance 2, no member is locally bad. Here the frustration of the odd structure is resolved by the
hard edge, and nothing local marks the two members.
A distance-3 example (n = 10, unique optimum, OPT = |H|+1 = 6, no bad member, c = 3): edges 01 08 12 14 23 34 45 49 56 67 78, H = (0,1,2) (2,3,4) (4,5,6) (6,7,8) (1,4,9); the only hard edge is {0,8};
the unique optimum {0,1},{2,3},{4,9},{5,6},{7,8} deletes it; witnesses 1 (soft to 0) and 7 (soft to 8); the members of 0 and 8 are at I(H)-distance 3. (n = 11, K = 3, 800 instances: none found among those
sampled, evidence only.)
So "bad" (a purely local, static notion) cannot be the criterion; the members at which a far shared edit occurs are determined by the global frustration of the base CSP.

### 5.3 Far equality of soft-pair states (insertion-only merges)  VERIFIED-BY-COMPUTER (`w6b_merge.py`)
Two pieces A != B in one cluster force EVERY soft pair between them to be merged. If two such soft pairs lie in members that are far in I(H), their states must be equal although no bag contains both
members. n = 8, K = 2: 35 of 2828 instances have this in every optimum; n = 9, K = 2: 19 of 1435; none at K = 3 in the n = 9 sample (small). The base CSP enforces this equality only through
static triples, which is why it is invisible at l = 0.

### 5.4 Why a bounded-state encoding of the agreement is impossible without identity  (IDEA, not a theorem)
A far two-sided edit e = (a,b) is a static object of G (a hard pair). Both ends can name it by an index of size <= D+l from their own side, but the match of the two indices must be verified
between two bags that contain neither a nor b together. Hash/colour variants (distance-2 proper colouring of the hard graph) identify b from a only inside a ball; two far vertices with the same local
fingerprint (identical gadgets) are indistinguishable, and a false match makes the DP pay once for two edits (undercount), so it is unsound unless the hash separates all relevant vertices, which takes
log n bits per pending edit. This is the content of Prop 10 of wave5_dp_route.md, now with an explicit, minimal example (5.1) and the statement of which edits need it (two-sided, far).

## 6. Where a hardness proof would have to look  (IDEA, unverified; no reduction completed)
Mechanism suggested by 5.1 and W1: a shared edit saves exactly 1 relative to charging the two ends separately. Take k regions, each a chain of members (I(H) a path, w = 1, c = 2) that must select ONE of n
sites by the position of the site; far hard edges (tunnels, each a 2-vertex hard component, so HE and HUB are not violated) connect the site j of one region to site j of another, and the cost of the
selection is lower by 1 exactly when the two regions select tunnel-linked sites. Multicolored Clique then asks for k selections that agree on k(k-1)/2 pairs of regions; l = O(k^2), w = O(1). I did not
build it. Obstacles that are NOT resolved: (a) forcing exactly one selection per region at cost independent of the site (a chain with a token seems possible); (b) an inactive site must cost nothing
(true in 5.1: the default states are conflict-free); (c) the sharing gain must not be cancelled by the cascades (triples with two edits); (d) pair-maximality of the whole gadget; (e) BFK (c = 1) kills the
construction for vertex-disjoint H, so c = 2 is needed and the HUB-type counting limits must be rechecked. Statements that this is W[1]-hard are NOT made.

## 7. Plan / conjectures (no proofs)
* CONJECTURE Lambda: let Lambda_l be I(H) plus an edge between the members of a and b for every hard pair {a,b} that can be two-sided within budget (there are soft partners t of a and t' of b with
  b in Delta(a,t) cup Gamma_h(a,t), a in Delta(b,t') cup Gamma_h(b,t'), weights <= l). Then CEaMP is FPT in (l,D,tw(Lambda_l)). Needs: a sound flag-DP for one-sided edits (cascades), checked agreement of
  two-sided flags inside bags. Not proved; 5.2 shows that Lambda_l cannot be replaced by I(H) plus bad members.
* Status of NF-bad: REFUTED for both K = 2 and K = 3 (n = 10).
* The static graph Lambda_l is the correct "interaction graph" for pointers; Gamma^l (Theorem P) contains it up to squaring but is larger. The open question (ML) is whether instances with large tw(Lambda_l)
  can be reduced; the expander family shows tw(Lambda_l) is unbounded for l above the weight, but those instances have OPT = |H|.

## 8. Status table

| item | status |
|---|---|
| W1 conflict lemma (replaces and sharpens near-twin Lemma: a in M(z,b) subset Delta(z,b), |Delta| <= l) | PROVED, VERIFIED (1.3e5 partitions) |
| W2 defect bound sum c(pi) <= 2D|U|, tight; Corollary LB (bag-local lower bound, exact at l = 0) | PROVED, VERIFIED |
| F1 merges with |Delta| > l / unmerges with |Gamma_h| > l impossible | PROVED, VERIFIED |
| DEL rigid hard edges (>= l common hard neighbours) | PROVED, VERIFIED |
| BAD at most 2c*eps locally bad members | PROVED, VERIFIED (1.7e5 partitions) |
| FE far end of a T2 insertion is transparent | PROVED |
| AN weak anchoring (a witness has a soft pair or another edit); strong form | PROVED; VERIFIED (0 violations, n <= 10) |
| Theorem P: FPT in (l,D,tw(Gamma^l)), strictly stronger than the omega theorem | PROVED (modification of dproute.tex Thm omega) |
| Proposition FAR: l = 1, w = 1 instances whose optima all delete a hard edge at I(H)-distance m-1 | PROVED; VERIFIED m = 3,4,5 |
| side-additive (per-end) accounting is incomplete | VERIFIED m = 3,4,5 |
| NF-bad (bad members in every bag suffice) | REFUTED (n = 10: K = 2, 7 of 4000; K = 3, 3 of 8000) |
| far equality of soft-pair states in all optima | VERIFIED (n = 8, 9, K = 2) |
| FPT in (l,w) | OPEN |
| hardness via tunnels + selection (Multicolored Clique) | IDEA only |
| CONJECTURE Lambda (FPT in tw(Lambda_l)) | CONJECTURE |
