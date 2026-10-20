# Wave 11: LOC and BC (small certificates), the lock calculus, and where the obstruction really sits

Labels: **PROVED** (proof here; "by computer" marks a finite classification that is part of the proof), **VERIFIED-BY-COMPUTER** (evidence only),
**HEURISTIC**, **OPEN**.  Code: `code/w11_*.py`, written from scratch (pysat, ortools CP-SAT, numpy); nothing imported from earlier waves.
`paper/` untouched, no git.  Conventions as in `soft.tex`: member (a,b,c), centre b, ab, bc edges, ac non-edge; states L,R,M = edit ab / bc / ac;
*local* = adjacency-defined (Theorem L), *remote* = not local; F'(T), certificate T for budget t = |H|+l (F'(T) > t).

## 0. Verdict

1. **BC at l = 0 is PROVED for c <= 2, with the sharp constant 4** (Theorem MC0): if F'(N) > |H| then some certificate has at most **4** remote pairs; 4 is attained
   (explicit instance, Sec. 3.4).  No pair-maximality needed.  This answers Open Problem `op:mincert` ("can a certificate of bounded size always be chosen")
   for l = 0, c <= 2, **including soft triangles** (the family G_N of Theorem CE lives there: its minimum certificate is 3, and that is now a theorem, not an accident of the design).
   The certificate can be found in time f(w) n^O(1) (enumerate <= 2 locks, test with the leaf test of Theorem L).
2. **Reduction for every c** (Theorem R0, PROVED): at l = 0 the minimum certificate has size <= 2 kappa_1 (<= 3 if one lock is empty), where kappa_1 is the
   *unit-lock Helly number* of the zero-excess set Sigma_0.  So "bounded certificates at l = 0" **is equivalent, up to a factor 2, to "singleton lock cores of Sigma_0 are bounded"**:
   exactly the missing "real lock chain" (PFG) object of waves 6-9.  Soft-triangle-free: kappa_1 <= 2 (any c, w).  c <= 2 with soft triangles: kappa_1 <= 2 (Theorem MC0, by computer).
   c >= 3 with soft triangles: OPEN; no chain with kappa_1 >= 3 that has bounded width was found (Sec. 5).
3. **BC for l >= 1 is NOT proved.**  What is proved: a sufficient criterion (Lemma D: l+1 zero-excess cores with pairwise disjoint regions give a certificate of <= 4(l+1) pairs),
   a pattern-wise statement (every defect pattern is killed by <= 2 locks), and a **sharp obstruction to the natural proof**: the weighted gluing (sunflower + majority) step
   is *impossible* already for the simplest soft-triangle-free link relations over the five states (SAT-verified, Sec. 4.3).  Evidence that the statement is nevertheless true for c <= 2:
   hill-climbing and 12,000+ periodic-chain searches never exceed roughly 4(l+1) (Sec. 4.5).  The proof would need an "independent witnesses" bound (Sec. 4.4); that is the exact remaining gap.
4. **Kill (2) (a family where every certificate grows while l stays bounded) is impossible at l = 0, c <= 2** (it contradicts Theorem MC0); at l >= 1 and for c >= 3 no such family was found
   and Sec. 4.4 / Sec. 5 say precisely what one would have to look like (relative-kill chain; unbounded kappa_1).  So B1 with "a smallest certificate" is **not refuted anywhere**, and is proved bounded in the base case.
5. **LOC (3): no progress beyond equivalence.**  Certificate pairs *are* witnesses in all tests (157/157 pairs of minimum certificates lie in an induced P3 with exactly one soft pair, so Lemma `rpwit`(b)
   absorbs them by one exchange; not proved in general), but a bounded certificate already gives FPT by branching, with no repacking needed; BC => a list of K^l candidate packings (list-LOC).
   LOC proper (deterministic, ell-free) stays open and is, by Theorem N, not implied by anything l-only.

Status of the two statements: **LOC: OPEN.  BC: proved for l = 0 (c <= 2; sharp), open for l >= 1 and for c >= 3.**  Neither is killed.

## 1. The lock calculus at l = 0 (PROVED)

At l = 0, budget |H|: F'(T) <= |H| iff there is s in {L,R,M}^H (exactly one edit per member, no local hard edit) satisfying every vertex triple whose remote pairs all lie in T.
Let sigma(tau) = number of soft pairs of a triple.

**Lemma 1 (triples with a remote pair).** If a triple tau has a remote pair then sigma(tau) <= 1.  If sigma = 1, tau = (z; x, y) with xy owned by a member A, z not in V(A),
and zx, zy are local hard or remote.  If sigma = 0 all three pairs are hard.
*Proof.* Two soft pairs of distinct members A,B share a vertex q in V(A) cap V(B); the third pair joins a vertex of A to a vertex of B, so it is local.  Two soft pairs of one member force the third to be soft too.  QED

So every constraint of F'(T) is (i) a *static* test (sigma = 0; at l = 0 all pairs are constants), (ii) a *unary list* on a member (sigma = 1), (iii) a binary relation R_AB (sigma = 2, adjacent members) or ternary (sigma = 3, soft triangle), the last two without remote pairs.
Let Sigma_0 = the zero-excess set with every triple containing a remote pair dropped.  Then Sigma(T) = Sigma_0 cap (lists activated by T).

**Lemma 2 (lock table; also `w11_table.py`).** For an outside vertex z and A = (a,b,c) put S = N(z) cap V(A) (all pairs from z to A are constants at l = 0).  The three triples (z; x, y), xy a pair of A, allow exactly

| S | empty | a | b | c | ab | ac | bc | abc |
|---|---|---|---|---|---|---|---|---|
| allowed states | L,R,M | L | none | R | R | none | L | M |
| (five states) | all | L,4 | 4,5 | R,4 | R | 5 | L | M |

**Every singleton entry is produced by one single triple**, i.e. by at most two remote pairs: "ab apart" = {L} (needs za != zb), "bc apart" = {R}, "ac together" = {M} (za = zc = 1);
S = a: triple (z;a,b) gives {L}; S = c: (z;b,c) gives {R}; S = ab: (z;b,c) gives {R}; S = bc: (z;a,b) gives {L}; S = abc: (z;a,c) gives {M}.  The "none" entries (S = b or S = ac) need two triples, three remote pairs.
A single triple with a non-singleton list is a 2-element guard ("ac apart" = {L,R}, "ab together" = {R,M}, "bc together" = {L,M}).  *This* is where G_N lives (guard "pair p r apart" of C_j).

**Lemma 3 (strengthening).** Let T be a certificate (l = 0).  Let T^cl be the closure of T under: if a remote pair zx is in T^cl and x lies in a member A with z not in V(A), add every remote pair between z and V(A).
Then T^cl is a certificate (superset), and every guard active under T^cl is the *full* list table(S_{z,A}) of Lemma 2 of some pair (z,A): a singleton, or empty, or trivial.
*Proof.* By Lemma 1 each active guard is a triple (z; x, y) with xy in A; if one of its remote pairs zx or zy is in T^cl, closure puts all remote pairs z-A into T^cl, so all three triples of (z,A) are active, and their intersection is table(S).  QED

**Theorem R0 (reduction, any c).** Let kappa_1 be the maximum size of an inclusion-minimal family of unit locks x_A = s (including two locks on one coordinate) that is inconsistent with Sigma_0.
If F'(N) > |H| then there is a certificate with at most 3 pairs (an empty lock or a static violation) or at most 2 kappa_1 pairs.
*Proof.* Lemma 3: T^cl is a certificate whose active guards are the full locks of Lemma 2; if one is empty, its <= 3 remote pairs suffice (the triples of that lock only need the pairs of (z,A)).
Otherwise the full locks form a family of unit locks inconsistent with Sigma_0 (or Sigma_0 = empty: certificate of size 0, or a static violation: <= 3 pairs); take a minimal inconsistent subfamily, size <= kappa_1; replace each lock by the single triple that produces it
(Lemma 2), <= 2 remote pairs; the active guards under the resulting set contain these locks, hence it is inconsistent with Sigma_0.  QED
Consequence: for l = 0 the question "is a bounded certificate always available" is *exactly* the Helly number of unit locks of Sigma_0 (lower bound too: a family of kappa_1 independent locks cannot be shortcut).

**Soft-triangle-free (any c, w).** Sigma_0 is closed under every majority (Theorem `sfA`); Lemma `sfH2` gives kappa_1 <= 2; hence a certificate of <= 4 pairs exists (the paper's Theorem CB0 bounds *every* minimal certificate by 12; this is the complementary "some certificate has <= 4" statement).

## 2. Theorem MC0 (c <= 2, with soft triangles): kappa_1 <= 2, minimum certificate <= 4 (PROVED; finite classification by computer)

**Structure of Sigma_0 at c <= 2.**  Members share at most one vertex, every vertex is in at most two members, so I(H) has maximum degree 3 and its triangle-connected components
are K3 (one triangle), the diamond K4 minus an edge, and K4 (octahedron; it is a whole connected component of I(H)).  By Lemma 1 every constraint of Sigma_0 is a unary list, or a relation on one triangle-connected component
(its members' vertices, all pairs between them constants or soft), or a link relation R_AB for an edge of I(H) in no triangle.  (An "outside" pair that is soft for a member E not in the component would make E part of the component, except
for the pair of the two private vertices of a diamond; the code allows that case.)  So Sigma_0 is an intersection of lifted relations of the types **link, K3, diamond, K4**, plus lists.

**Classification (`w11_maj2.py`, exhaustive over centre choices and all hard bits).**

| type | configurations | distinct relations over {L,R,M}^k | closed under the M-majority m(x,y,z) = x if two equal, else M | majorities (of 729) preserving all |
|---|---|---|---|---|
| link (2 members) | 9 x 16 | 32 | all | 729 |
| K3 | 27 x 64 | 92 | all | 1 (the M-majority) |
| diamond | 81 x 2^8 x 3 | 332 | all | 1 |
| K4 | 81 x 8 | 254 | 248 (6 fail) | 0 |

(The numbers 163/445/39/521 of wave 9b count ordered role orders; the present counts are up to the L<->R relabelling of a member, which the M-majority commutes with.)
For K4 the M-majority fails on 6 relations, but a **K4 component has only four coordinates and its unit-lock Helly number is 2 anyway**: `w11_kappa.py` computes kappa_1 of all 254 relations and of all
840 distinct sub-relations R cap (product of arbitrary non-empty lists): kappa_1 in {1,2} for K4 (and for the other types).

**Theorem MC0.** Let c <= 2.  Then kappa_1(Sigma_0) <= 2, hence if F'(N) > |H| there is a certificate with at most 4 remote pairs (at most 3 if it is an empty lock).
*Proof.* Sigma_0 is a product over the connected components of I(H), and lists are applied inside components, so an inconsistent family of unit locks lives in one component.  A K4 component: kappa_1 <= 2 by the table above.
Any other component: Sigma_0 restricted to it is an intersection of link, K3 and diamond relations and lists, each closed under the M-majority (table; lists are preserved by every conservative operation, and the M-majority is conservative).
By Lemma `sfH2` (Helly number two for unit lists under a majority) kappa_1 <= 2.  Apply Theorem R0.  QED

### 2.1 Evidence / checks (all VERIFIED-BY-COMPUTER, independent of the proof)
* `w11_t1.py`: random c <= 2 instances (general slot matchings, triangles, diamonds, K4, random hard clusters, pendants; about 25,000 instances in total across the logs): Sigma_0 was closed under the M-majority in
  every instance without K4 component and in every K4 instance tried except the expected rare ones (`w11_t2.py`: 4 of 5,333 octahedron instances are not closed); exact minimum certificate sizes (CEGAR with CP-SAT master and SAT checker) were in {0,1,2,3,4},
  never above 4.
* `w11_t3.py`: same without pair-maximality (independent random hard pairs): sizes <= 3, no violation.
* Hill-climbing (`w11_hc.py`, l = 0, 5-8 members, up to 6 pendants): maximum 4.
* **Tightness of 4** (`w11_tight.py`): c = 2, pair-maximal, no soft triangle, n = 13, H = {(1,0,3),(5,0,4),(2,6,1),(2,7,8)}, clusters [-1,1,-1,0,1,0,-1,-1,-1,1,-1,-1,0]: minimum certificate has exactly 4 pairs {(0,9),(0,12),(4,9),(5,12)}; all 24,804 three-subsets of its 54 remote pairs are feasible.

### 2.2 What this says about G_N
`w11_gn.py` reproduces Theorem CE with new code (valid, pair-maximal, c = 2, N soft triangles, OPT - |H| = N for N <= 2 by an independent subset DP, T_N minimal for N <= 4),
and the exact minimum certificate is 3 for N = 1..4 (pendant at a centre).  Theorem MC0 shows that this is forced: **no choice of pendants can make every certificate large at l = 0 when c <= 2.**
A systematic search over *all* pair-maximal pendant attachments of the K3 chain (`w11_pat2.py`, 11 valid neighbourhood patterns for middle/start/end pendants, N = 4 and N = 6, 336+38+69 combinations):
minimum certificate in {3,4} independent of N, never growing.  Mechanism: a pendant whose neighbourhood makes the 2-guard {L,R} on C_j also acts on A_j (centre p_j), so the *full* list of the same vertex is a singleton or empty one;
any benign choice turns the guards into locks, which have Helly number 2.

## 3. Why the hard direction has to look different (what an l = 0 counterexample needs)
By Theorem R0 an unbounded minimum certificate at l = 0 needs **unbounded unit-lock cores of Sigma_0** at fixed (c, w): a chain of blocks whose zero-excess set carries a bit through the chain such that a *singleton* lock
kills only the escape state of its block.  At c <= 2 this is impossible (MC0); soft-triangle-free impossible (sfA).  At c >= 3 Lemma R5 shows the scope language has no conservative near-unanimity polymorphism, so by Baker-Pixley
*some* pp-definition has unbounded singleton cores; whether one can be realised with bounded treewidth and as a packing is open (Sec. 5).  This is the same object that a W[1]-hardness proof through tunnels would need
(`op:chain`), now seen from the certificate side: **BC at l = 0 and the chain-object question are the same question for c >= 3.**

## 4. l >= 1

### 4.1 Setting: patterns and regions
For budget |H|+l the objects that cost are *defect objects*: a member edited at least twice (extra cost >= 1) and an edited local hard pair (cost 1); cost = |H| + sum of extras + number of flips >= |H| + #defect objects.
A *pattern* pi is a set of at most l defect objects with values.  For a pattern, Sigma_0^pi is the zero-excess set of the instance with the pattern applied (flips applied; defect members replaced by their pair statuses, which become hard constants).

**Lemma P (pattern-wise Helly, PROVED by computer).**  For c <= 2 and every pattern pi (every F-flip set, every frozen defect member), kappa_1(Sigma_0^pi) <= 2.
*Proof.* Flips change hard bits of enumerated configurations; freezing a member turns its soft pairs into hard constants, i.e. replaces a K3/diamond/K4 relation by a slice that is a relation of the same enumerated families
(slices are intersections with lists); lists are arbitrary; the classification of Sec. 2 covers all of them.  QED
So **every single pattern is killed by <= 2 locks (<= 4 remote pairs)**; the difficulty is to kill *all* patterns with one bounded set.

### 4.2 Sufficient criterion
For a zero-excess core C (two locks, minimal inconsistent sub-CSP on a member set M(C)) put D(C) = M(C) cup {local hard pairs occurring in a triple with sigma >= 1 whose soft owners lie in M(C)}.

**Lemma D (disjoint regions, PROVED).**  If T contains l+1 zero-excess cores with pairwise disjoint D(C), then F'(T) > |H|+l.  Hence if there are l+1 such cores among all locks, a certificate with at most 4(l+1) pairs exists.
*Proof.* An assignment x feasible for T with no defect object in D(C_i) restricts to M(C_i) as a state vector satisfying all constraints used by the core and its locks: contradiction.  So x has a defect object in each D(C_i); they are distinct objects, cost >= |H| + l + 1.  QED
So a minimum certificate can exceed 4(l+1) pairs only if the maximum number nu of disjoint-region cores is <= l.  (HEURISTIC, not proved: if core regions are connected sets of members, they are subtrees of a decomposition, whose transversal number equals the packing number, and a transversal by members has size <= (w+1) times that; then **all remaining conflicts pass through a bounded part of the decomposition**.)

### 4.3 The natural proof fails: no weighted gluing (SAT-verified)
Weighted Helly (the "increment lemma": if F'(N) > F'(Phi) then a bounded C has F'(Phi cup C) > F'(Phi); with it B1 would run in (4(l+1))^l f(w) n^O(1) by iterating increments) would follow from a *gluing* of three
optimal assignments x^a, x^b, x^c violating one lock each, with defect sets forming a sunflower (Erdos-Rado: m > l!(k-1)^l sets contain k petals), if coordinatewise majority produced a feasible assignment whose defects lie in the kernel.
`w11_maj5.py`, `w11_maj5d.py`: over the five states {L,R,M,4,5} (4 = edits ab,bc; 5 = all three), all coordinates of the three vectors are either equal or have at most one defective value (the sunflower pattern).
**No operation on those patterns closes the relations**, even if its output may be defective: UNSAT for the K3 relations (53 of 785 distinct five-state K3 relations are individually unclosable), and for the *soft-triangle-free link relations*:
with non-defective outputs UNSAT (a single link relation already, e.g. ZE part {LL,LR,LM,RR,MR} with defect states), with defective outputs SAT (cost then grows by up to 2 per petal coordinate: no budget-preserving gluing).
This is the Proposition `hub` phenomenon (breaks lock their neighbours) in its sharpest form and explains why exchange arguments fail.  With the six hard pairs of a K3 as flip variables the relation is not closed under the pair (M-majority on states, Boolean majority on flips)
(`w11_maj3.py`: all 64 hard-bit configurations of the first role order fail; the cause is cluster editing on the three private vertices: three flip bits under 'not exactly two edges').

### 4.4 What remains: independent witnesses (HEURISTIC reduction, OPEN statement)
Let T be a minimum certificate with m locks, x^j (cost <= l) a witness violating only lock j, D_j its defect set.  Sunflower lemma on {D_j}: if the kernel is non-empty, freezing one kernel object o puts all petal locks in a certificate of the frozen instance
(budget l - c(o), same class), which is where an induction on l would act (HEURISTIC: the induction needs the frozen instance's *minimum* certificate to be comparable with sub-certificates of T; I could not close it).
If the kernel is empty the witnesses carry pairwise disjoint defects: "k locks with independent single-lock witnesses".  At l = 0 this is Helly number 2 (witnesses are zero-excess states: MC0).  For l >= 1 it is exactly the relative-kill chain:
the defect at position i is killed only by a lock whose conflict is with the *list induced by the defect*; a long chain needs (a) a rigid forced path from the first lock (so that zero-excess regions are long) and (b) defect-induced information that does *not* propagate (so that each position needs a private lock).
Argument that (a)+(b) are incompatible for c <= 2 (HEURISTIC): a lock consistent with the forced path, placed in the interior, forces the same rigid continuation and therefore creates a *shorter* core with the far conflicting lock; a 2-list lock does not, but its *full* vertex pattern is a singleton one of the two values, giving either a prefix core [1,j] or a suffix core [j,k];
among the two families the smallest prefix and the largest suffix core are either disjoint (nu >= 2: Lemma D) or leave a gap of consecutive positions without locks, where nothing kills a defect (then F'(N) <= |H|+l and no certificate is needed).
I could not turn this into a proof because defects change the interface in more ways than "cut" (they induce lists on neighbours, which are information-preserving along rigid chains).

### 4.5 Evidence for l >= 1, c <= 2 (VERIFIED-BY-COMPUTER, not a proof)
* Exact minimum certificates (CEGAR) found by hill-climbing over (tree-like member structure, hard clusters, pendants), n <= 20:  l = 0: **4**;  l = 1: **8** (several seeds; members 5, 6, 7, 8; stays 7-8 for m = 5..10);  l = 2: **10**; l = 3: **10** (still climbing at the time of writing); l = 4: **14**
  (see the final numbers in Sec. 8).  Compatible with K(l) = 4(l+1) = bound of Lemma D; none exceeds it.  Certificates found at l = 1 are unions of two to four 2-pair guards.
* Instance with intrinsic baseline excess F'(empty) = |H|+3 (`w11_inc3.py`): minimum certificates 1, 2, 11 pairs for l = 3, 4, 5, i.e. effective l' = 0, 1, 2 give 1, 2, 11 (11 <= 4(l'+1)+... compatible with 12).
* Increment sizes along greedy chains (`w11_inc2.py`): hill-climbed maximum increment 9 (chain [1,1,9]) and 7 (chain [2,2,3,7,3]); the increment lemma for arbitrary Phi is therefore *not* a-priori <= 4, but the minimum certificate of the same instance is 11 pairs.  The increment lemma is not needed for BC itself.
* Periodic chains (`w11_cell.py`: random cell with k = 3,4 members, in-cell sharing of multiplicity <= c, ports, cluster labels with offset to the next cell, pendants), N = 3 vs 5: l = 0, c = 2 and c = 3: 24,000 cells, **no growth** (sizes (1,1),(2,2),(3,3),(4,4) only);
  l = 1, c = 2: 3,000 cells; one instance went 5 -> 6 between N = 3 and 5 and then stays at 6 for N = 6, 7, 9, 12 (`w11_growth.py`).
* K3-chain with *all* pair-maximal pendant patterns, l = 1 (`w11_pat3.py`, `w11_pat4.py`): minimum certificates 3,5,6,7,8; for 120 random valid pattern combinations N = 3 and N = 5 give identical sizes.

## 5. c >= 3
* Theorem R0 holds for every c: BC(l=0) <=> bounded unit-lock cores of Sigma_0 (up to 2).  No bound is known with soft triangles.
* Lemma R5: no conservative NU polymorphism for the c >= 3 scope language.  By Baker-Pixley, in the pp-closure there are relations with unit-lock cores of every size; the open point is bounded treewidth + realisability as a packing.
* Search: 12,000 random c <= 3 periodic cells at l = 0 (`w11_cell.py 3 3`, `4 3`): maximum minimum certificate 4, no growth between N = 3 and 5.  The cells do not specifically target the R5 relation, so this is weak evidence only.
* Wave 9b's annealing for kappa_1 (<= 3 everywhere) remains the best evidence that BC(l=0) is also true for c = 3.  A proof or a counterexample needs a direct analysis of the Rk-type relations with bounded treewidth.

## 6. LOC (item 3)
* BC => FPT (Theorem B1 with the leaf test, no repacking).  With the *list* version: B1 produces <= K^l edit sets; each remote pair e of a certificate is a witness of an absorbable far edit (Lemma `rpwit`(b)) whenever it lies in an
  induced P3 with exactly one soft pair; so BC gives a list of <= K^l repackings H*, one of which has the optimum in class S.  `w11_wit.py`: in 81 random c = 2 instances all 157 pairs of minimum certificates lie in such a P3 (VERIFIED-BY-COMPUTER, not proved;
  a lock triple that is a triangle in G is not a witness of the edit it forbids, so the general statement needs the extra argument that a minimum certificate contains a P3-type triple for every pair).
* Hence LOC is *not* advanced independently: (a) the information problem "which of ~n^2 pairs" is solved by BC only for l = 0 (where nothing needs to be located), (b) Theorem N still forbids l-only bounds for any deterministic H*.
* One structural by-product relevant for LOC: at l = 0 locks are *vertex-local* objects (z,A); a remote edit repairs the lock (z,A), so the first end of a far edit ranges over the guard vertices z of a bounded number of locks, not over n positions.

## 7. Status table

| item | status |
|---|---|
| Lemma 1-3, Theorem R0 (l = 0, any c: min certificate <= 2 kappa_1 / 3) | PROVED |
| Classification of the c <= 2 constraint relations (32 link, 92 K3, 332 diamond, 254 K4 relations; M-majority closes link, K3, diamond; kappa_1 <= 2 on 840 K4 sub-relations) | VERIFIED-BY-COMPUTER (exhaustive; part of the proof of MC0) |
| **Theorem MC0** (c <= 2, l = 0: a certificate with <= 4 pairs exists; 4 attained) | **PROVED** (finite check by computer) |
| Soft-triangle-free, any c: some certificate <= 4 | PROVED (Theorem sfA + Lemma sfH2 + R0) |
| Lemma P (every pattern killed by <= 2 locks, c <= 2) | PROVED (by computer) |
| Lemma D (l+1 disjoint-region cores => certificate <= 4(l+1)) | PROVED |
| No budget-preserving majority gluing over five states (links; K3) | VERIFIED-BY-COMPUTER (SAT, exhaustive over the 42 patterns) |
| BC for l >= 1, c <= 2 (K(l) = O(l)) | OPEN; evidence: max found 4, 8, 10, 10, 14 for l = 0..4 |
| Independent-witness bound (empty-kernel sunflower) | OPEN (heuristic argument for c <= 2, Sec. 4.4) |
| BC / kappa_1 for c >= 3 | OPEN; no growth in 12,000 cells; Lemma R5 predicts unbounded cores in some pp-definition |
| Existence of a family killing B1 (every certificate large, l bounded) | none found; impossible for l = 0, c <= 2 (MC0) |
| LOC | OPEN (BC => list-LOC; certificate pairs are P3 witnesses in all tests) |

## 8. Final numbers of the long runs
(see the end of this file; filled in after the runs completed)

## 9. Reproduction (from `code/`; python3 + pysat + ortools + numpy + networkx)
* Library: `w11_lib.py` (instance, F' oracle with incremental SAT, `min_certificate` = CEGAR with CP-SAT master, `min_increment`, `fprime_value`, exact OPT by subset DP, n <= 16), `w11_gen.py` (pair-maximal generators: members + hard clusters; c <= 2 slot matchings; c <= 3).
* G_N and chains: `w11_gn.py 4`, `w11_pat.py`/`w11_pat2.py 4` and `6`, `w11_pat3.py 3 1`, `w11_pat4.py 1 3,5 120`.
* Theorem MC0: `w11_table.py`, `w11_maj.py` (single-scope view, shows why scopes must be taken per triangle-connected component), `w11_maj2.py` (the classification), `w11_kappa.py`,
  stress: `w11_t1.py 6 4 4000 101`, `w11_t2.py 7 6000 5`, `w11_t3.py 1 1500 5 3`, `w11_tight.py`, hill-climb `w11_hc.py 6 6 0 1 1500 3`.
* l >= 1: `w11_hc.py m npend l seed iters restarts`, `w11_inc.py`, `w11_inc2.py`, `w11_inc3.py`, `w11_cell.py k cmax samples seed Ns l`, `w11_growth.py`, `w11_maj3.py`, `w11_maj5.py 3|5`, `w11_maj5b.py`, `w11_maj5d.py`, `w11_abs.py` (abstract chain model; **discarded**: states-only chains have spurious infinite costs because the hard constants forbid recovery, flips are needed).
* c = 3: `w11_cell.py 3 3 6000 11 3,5`, `w11_hc3.py`.  LOC: `w11_wit.py 1 400 5 4`.
