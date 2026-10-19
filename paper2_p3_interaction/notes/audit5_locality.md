# Audit 5 (hostile referee, wave 5): notes/wave5_locality.md

Target: Lemma 1, Theorem 1 (T1/T2), Corollary 2, Lemma 3, the stated n=8 W-locality counterexample.
Own code (written from scratch, imports only a5_lib.py; nothing from a4_lib or w5a_*): code/a5_lib.py, a5_check.py, a5_cex.py,
a5_small.py, a5_merge.py, a5_far.py, a5_farmerge.py, a5_l3.py, a5_tw.py. No git, paper/ untouched.

## Verdicts
| item | verdict |
|---|---|
| Lemma 1 (half-degree; diameter <= 2; middle z) | CORRECT |
| Theorem 1 (T1/T2) | CORRECT |
| Cor. 2, proof of (i),(ii) | CORRECT |
| Cor. 2, bag locality, width (w+1)(3w+1)-1 for I(H)^2, W-ball radius 4 | CORRECT |
| Cor. 2, partner bound "D(D-1) <= 4c(c-1)" | CORRECT WITH FIX (the inequality as written is FALSE; the target bound 4c(c-1) is right and has a one-line proof) |
| Cor. 2, scope / "assigned gap CLOSED" (Sec. 0 item 1, Sec. 7 table) | CORRECT WITH FIX (overclaim: the hypothesis "no soft pair between A and B" is not the general tiny-tiny merge; counterexamples below) |
| Lemma 3 (|{a}| <= u; |M(z,b)| <= l) | CORRECT |
| stated n=8 example | CORRECT as an example; the word "smallest" is WRONG (n=6 and n=7 examples exist; n<=5 none) |
| "far T2 at infinite I(H)-distance" evidence (n=9, c=1) | TRUE but IRRELEVANT to merge pointers: the inserted pair is inside one piece |

## 1. Lemma 1 -- CORRECT
(a) Moving v to a new singleton: its d_C(v) kept edges become deletions (+d_C), its |C|-1-d_C inserted non-edges become free
(-(|C|-1-d_C)); pairs to other clusters unchanged. Optimality gives 2 d_C >= |C|-1. Checked line by line, no gap (also |C|=1 trivial).
(b) a,b non-adjacent in C: N(a) and N(b) inside C are subsets of C\{a,b} (a not in N(a), a not adjacent b), sizes sum >= |C|-1 > |C|-2,
so they meet. Correct. Only local optimality under single-vertex moves is used (true; stated in the note).
Computer (a5_check.py, independent): all optimal partitions of random pair-maximal instances n=6,7,8 (291 + 149 + 40 instances, 744 + 517 + 137
optimal partitions): 0 violations of (a) or (b). Control: on cost OPT+1 partitions (a) fails in 166/291, 45/149, 27/40 instances, so the check is not vacuous.

## 2. Theorem 1 -- CORRECT
z from Lemma 1(b); all-hard (a,z,b) would be an induced P3 with all three pairs unpacked (ab hard non-edge assumed, az, zb edges): contradicts
pair-maximality. Both soft in the same member => ab packed, contradiction. A hard edge inside a cluster is kept, so its ends are in one piece.
Nothing missing (z distinct from a,b by construction). a5_check.py: 0 failures of the T1/T2 dichotomy on all optimal partitions above.

## 3. Corollary 2
### 3.1 Proof of (i),(ii): CORRECT
(i) cross pair hard; hard edge in one cluster is kept and would merge A,B; so hard non-edge, inserted. (ii) z from Thm 1. T2 case: say az hard,
then z in A (kept edge); zb in E is a pair A x B: hard => kept => z,b same piece, contradiction; soft => excluded. So T1. z in A or B gives again a
soft pair between A and B (zb or az). So z in C\(A u B), both pairs soft, distinct members both through z. Correct. Computer: every one of the
cross pairs of soft-free piece pairs (27 + 22 + 2 pairs) was inserted and had a soft common neighbour outside A u B.

### 3.2 Bag locality and the width bound: CORRECT
P_a, P_b share z, so are adjacent in I(H) (adjacency = sharing a vertex; packing.tex uses the same convention for W), hence in a common bag: a,b in V(X_t),
|V(X_t)| <= 3(w+1). Square decomposition (T, N[X_t]): vertex v lies in bags t with X_t meeting N[v]; the subtrees T_u, u in N(v), each meet T_v, so the union
is connected. Edge uv of the square with common neighbour m: any bag containing m contains N[m] and so u,v. Adjacent u,v: bag of u contains N[u]. Degree of I(H) <= 3(c-1) (3 vertices
per member, each in <= c members), c <= w+1, so |N[X_t]| <= (w+1)(1+3w). Width <= (w+1)(3w+1)-1. Verified constructively (a5_tw.py: min-degree-elimination decompositions
of I(H) for 4,000 random pair-disjoint triple systems, bags N[X_t] checked to be a valid tree decomposition of I(H)^2, 0 violations, size bound attained with equality 244 times).
Cosmetic: state that N[X_t] is a set of NODES of I(H), and that the sharper size is (w+1)(3c-2).

### 3.3 Partner bound: CORRECT WITH FIX
Written: "at most D(D-1) <= 4c(c-1)". The step D(D-1) <= 4c(c-1) is false: D = 2c exactly, D(D-1) = 2c(2c-1) = 4c^2-2c > 4c^2-4c (c=1: 2 > 0; c=2: 12 > 8). The bound 4c(c-1)
is nevertheless TRUE, by the sharper count D(D-2): given a and z in Q(a) (z shares member P_a with a), a partner b must lie in a member P_b != P_a through z (Cor. 2(ii)) and
b not in P_a (else ab packed). So b in Q(z) minus the two partners of z inside P_a: at most D-2 = 2(c-1) choices per z, total <= D(D-2) = 4c(c-1). c=1 gives 0, which is right (T1 needs
two members through z). Checked: max over a of |{b: common soft neighbour, ab hard}| never exceeds D(D-2) (a5_check.py, 0 violations; max observed 3,4,5 for n=6,7,8).
Also: the global count n*C(D,2) is fine (even n*C(D,2) minus the C(2,2)*c pairs inside members).

### 3.4 W-ball radius 4: CORRECT
Every cross pair a'b is at W-distance exactly 2 (not packed, common soft neighbour). So A is inside the radius-2 W-ball of any b in B and B inside radius 2 of any a in A;
A u B inside radius 4 of b. (Needs A,B nonempty, trivial.)

### 3.5 Scope of Cor. 2 -- overclaim, FIX required
The hypothesis is "no soft pair between A and B". The note's Sec. 0 item 1 says "The assigned gap is CLOSED: tiny-tiny pointers are W-local" and "the hypothesis |A||B| <= l is not needed".
What is proved is the soft-free case only. The general tiny-tiny merge (two pieces sharing a cluster; cross pairs split into soft pairs and <= u inserted ones, SM) is NOT covered, and it is the
generic case. Counterexamples (independent search, a5_merge.py), each with H a valid pair-maximal packing, EVERY optimum containing a MERGE-type inserted pair (a,b), a,b in V(H), a and b in
DIFFERENT pieces of one cluster, Q(a) and Q(b) disjoint (no common soft neighbour, so not T1; it is T2 with a soft pair between the pieces):
 * n=7 (smaller than the note's n=8): edges 03 04 12 13 14 15 23 35 36, H = {(4,1,5),(0,3,2),(1,3,6)}; unique optimum clusters {0,4},{1,2,3,5},{6}, OPT = 4 = |H|+1;
   pieces {1,2},{3,5} (hard edges 12, 35), cross pairs 13,15,23 soft, 25 inserted, Q(2)={0,3}, Q(5)={4,1}. Found 9 such instances in 7,815 random n=7 draws; none in 7,509 at n=6.
 * n=8: the note's own example (pieces {1,4},{5,7}, cross pairs 17,45,47 soft, 15 inserted).
In both, the offending pieces are tiny and share a cluster, so the general tiny-tiny pointer is not bag-local in the sense "common soft neighbour". (They are still I(H)-close: here members of a and b are at
I(H)-distance 2. I did NOT find, and the note does not exhibit, a MERGE-type forced pointer at unbounded I(H)-distance; see 5.) The statement actually proved (and true): soft-free piece pairs
have all cross pairs T1. Note also that in the paper's open-problem wording ("endpoints at soft distance 2 from V(H)") the candidate bound 4c(c-1) is a tautology (count of vertices at soft distance
exactly 2); the nontrivial content of Cor. 2 is only that, in the soft-free case, every cross pair is at soft distance 2. So the verdict line should read "T1 locality for soft-free merges; merges with a
soft pair between the pieces remain open", and the table row "Cor. 2 ... bag-local" must carry the qualifier. Also (d) "the profile of the partner is read off inside the bag" holds only for T1.

## 4. Lemma 3 -- CORRECT
a not in {z,b}, az, ab hard, [az in E] != [ab in E]. a in C: the non-edge among the two is a hard non-edge inside C: inserted, in U. a not in C: the edge among them joins different clusters: deleted, in U.
The map a -> that pair is injective (a is an endpoint, the other endpoint is z or b, and a is not z, b; {a,z} = {a',b} would force a in {a',b}). Hence <= u. Hypotheses: none on H or optimality; z,b in one cluster.
Anchor form: for a soft pair (z,b) in one cluster, a in M(z,b) = N_h(z) \ (N[b] u Q(b)) has az hard edge, ab non-edge and hard, a != z,b (b not in N_h(z) because zb packed). So M is contained in the set above, |M| <= u <= l,
computable from G, H alone. T2 usage: far end a of (a,b) with az hard, zb soft is in M(z,b): correct.
Computer (a5_l3.py, ALL partitions of n=7 for 60 random graphs, arbitrary pair-maximal H, 511,560 (z,b) checks, 236,292 anchor checks): 0 violations.

## 5. The stated n=8 counterexample and the "smallest"/"far" claims
* Stated example (a5_cex.py, a5_merge.py): H consists of 5 pair-disjoint induced P3s, pair-maximal (checked), OPT = 6 = |H|+1, unique optimum {0,6},{1,4,5,7},{2,3}, only U-pair is the inserted 15, Q(1)={6,7}, Q(5)={2,4}
  disjoint, "14 hard, 45 soft; 17 soft, 57 hard" confirmed; pieces {1,4} and {5,7}; so it genuinely is a merge-type T2 pair. CORRECT. It refutes only "every inserted HH pair can be taken at soft distance 2".
* "Smallest example" is WRONG. Exhaustive over all labelled graphs and ALL pair-maximal packings with n=4,5: no instance forcing an inserted pair, both ends in V(H), without a common soft neighbour. n=6 has one
  (a5_small.py): edges 02 03 04 05 12 13 15 25 35 45, H = {(3,5,4),(0,2,1)}, unique optimum {0,1,2,3,5},{4}, inserted 23 forced, Q(2)={0,1}, Q(3)={5,4} (03 hard edge kept, 02 soft: T2).
  But that pair is INSIDE one piece (hard edges 03,05,13,15,25 connect all of {0,1,2,3,5}), so it is not a merge pointer; the merge-type minimum is n=7 (above; 0 in 7,509 draws at n=6, 9 in 7,815 at n=7). Cosmetic.
* "Far T2 at infinite I(H)-distance (n=9, c=1, 5-8)": reproduced (a5_far.py: H = three disjoint P3s, 6 optima, 5-8 inserted in all). But the cluster is {4,5,7,8} with hard edges 47, 48, 57: 5 and 8 are in the SAME piece.
  The pair is determined by hard connectivity (a G-determined tunnel), not a free pointer; it does not show a non-local MERGE pointer. a5_farmerge.py (n=9, three disjoint P3s, 324 valid random instances,
  all partitions): 0 instances where every optimum has a merge-type inserted pair between different members. So verdict 3 ("T2 pointers are NOT local ... infinite distance") is supported only for intra-piece
  insertions; for merge pointers the I(H)-far question is OPEN, neither proved nor refuted (need larger instances than n=9, c=1).

## 6. Required fixes to wave5_locality.md
1. Sec. 0 item 1, Sec. 7 table: restrict the closure claim to merges with NO soft pair between the pieces; say merges with soft pairs between the pieces (generic by SM) are not covered; n=7/n=8 examples show they can be T2.
2. Cor. 2(b) and Sec. 0: replace "D(D-1) <= 4c(c-1)" by the count D(D-2) = 4c(c-1) with the argument b not in P_a.
3. "smallest example n=8": n=6 (intra-piece) / n=7 (merge-type); n<=5 none.
4. Far-T2 evidence: mark as intra-piece; do not cite as a non-local merge pointer.
5. Optional: Lemma 1, Thm 1, Lemma 3 as written can go into a paper unchanged (Lemma 3 needs "z,b in one cluster" stated; M defined for soft z,b).
