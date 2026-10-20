# Audit 10: wave 7 (notes/wave7_special.md, code/w7_*.py)

Hostile referee pass. Everything below was re-derived or re-run with **new code written from scratch** (`code/a10_*.py`, library `a10_lib.py`, numpy only, no import of any `w7_*`/`a*_` file).
Ground truth for every zero-excess (ZE) set Sigma is **full enumeration of all set partitions of V** (Bell(10) = 115,975), ZE = partitions of cost |H|; the decomposition lemma is only asserted, never used.
No git, `paper/` untouched. (The task text says "audit notes/wave7_special.md"; the file is `notes/wave7_special.md`.)

## Verdicts

| claim | verdict | one line |
|---|---|---|
| Lemma 1 (Sigma = unary lists + binary relations R_AB when no soft triangle) | CORRECT | case analysis re-done; end-to-end confirmed by own CSP-from-triples solver (below) |
| Theorem A (Sigma closed under every majority, no soft triangle) | CORRECT | 40-relation list complete and exact; hand proofs per class valid; 20,920 random instances, 0 failures |
| A1 (l = 0 polynomial, soft-triangle-free, any c, any w) | CORRECT (modulo Jeavons-Cohen-Cooper) | own path-consistency solver agrees with brute force on 20,920 + 3,400 instances; no conflict with LPS23 (their gadgets are full of soft triangles) |
| A2 (cores <= 2; C2 for triangle-free, c <= 2) | CORRECT | follows from Thm A + paper's Lemma `tmaj`; 0 three-cores in 1,602 (instance, anchor-set) brute-force checks |
| A3 (tunnel-release model FPT, 2^l n^O(1)) | CORRECT | the vertex-cover characterisation verified by direct brute force on G + stubs |
| Strictness Lemma S | CORRECT | re-checked on real graphs, 0 violations; but it is a one-edge statement, see remark |
| Lemma 2 (exact value for U-only model; locality) | CORRECT WITH FIX | identity exact (1,166/1,166); locality statement omits the static (all-hard) triples |
| Theorem B (LPPH W[1]-hard from Grid Tiling with <=) | CORRECT | hand proof checked line by line; 1,400+ instances against two independent solvers, 0 mismatches |
| Theorem C (tunnel route closed w/o soft triangles) | CORRECT WITH FIX | (i),(ii) hold once (ii) is stated through the vertex-cover characterisation; (iii) is not an independent obstruction; (iv),(v) are heuristics, not part of the theorem |

No claim is WRONG.

## 1. Is the 40-relation enumeration complete? (a10_rel.py)

Independent construction: every pair of members (a,b,c), (a',b',c') sharing exactly vertex q on 5 vertices, all 6 x 6 role orders, all 16 patterns of the four hard cross pairs: **576 real graphs**, all valid packings.
For each graph and each of the 9 state pairs I edited the graph and tested the four triples (p,q,r) for being an induced P3 on the edited graph (no x/eps formula).
Result: **40 distinct binary relations, sizes 0,1,2,3,5 with multiplicities 1, 9, 18, 8, 4**, exactly as claimed. Shapes: 3-sets = 2 rows + 2 columns (rows/columns p,q in {L,R}) + 4 subsets of {L,R}^2; 5-sets = the four crosses row p + column q with p,q in {L,R}. Nothing else.
Completeness of the *configuration space*: R_AB depends only on (role of q in A, role of q in B, four eps bits), because any two soft pairs of a triple share a vertex, and for two distinct members that vertex is q; the third pair (p,r) cannot be soft (owner A or B would force |A cap B| >= 2; a third owner is a soft triangle). Unary constraints from triples with one soft pair are subsets of {L,R,M}, closed under every majority. So the 576 graphs are the complete list.
Closure: all **729 majority rules enumerated explicitly** (3^6 choices on the six all-distinct patterns), applied to all triples of tuples of each of the 40 relations (and of the 32 distinct full 5-vertex Sigma): **0 failures**.
Hand arguments per class (note, Sect. 2): 0/1/2-element sets (at most two values per coordinate, majority is forced), lines (one coordinate constant), 3-subsets of {L,R}^2 (Boolean majority, every binary Boolean relation is closed under it), crosses (the case analysis on the first coordinate, using that two of three all-distinct first coordinates differ from p): I checked each; all valid. Note that the cross proof is the only one that is not a one-liner and it is right.
Robustness remark worth adding to the note: the state labelling L,R,M of a member depends on the orientation (a,b,c) vs (c,b,a) (swaps L,R). Because Sigma is closed under *all* majorities, which is a set invariant under relabelling by any permutation of {L,R,M}, closure does not depend on the orientation conventions. (Closure under one fixed m would.)

## 2. Closure on random instances, independent brute force (a10_closure.py, a10_tree.py)

* `a10_closure.py`: planted random instances, n = 5..10, 2..6 members, optional noise on hard pairs, c up to 4. Sigma from partition enumeration. Eight runs, 45,900 instances, **20,920 soft-triangle-free** (13,300 with Sigma non-empty, 1,300+ with |Sigma| >= 3): **all 20,920 closed under all 729 rules, all 2-decomposable, own path-consistency solver agreed with brute force on all 20,920** (feasible and infeasible). cmax histogram includes c = 3 and 4 (hubs), which are soft-triangle-free but not triangle-free in I(H).
* `a10_tree.py`: tree-like packings (each member shares one vertex with the earlier ones, hubs of 4 members), n <= 9, random hard pairs: 3,400 instances, 700+ non-empty, 186 with |Sigma| >= 3: all closed, 2-decomposable, path consistency exact.
* Negative control (harness can fail): the same code on instances **with** soft triangles finds 451 not closed under all rules and 75 not 2-decomposable (in the four runs that recorded c: 34 not 2-decomposable, **all at c = 3, none at c = 2**; 233 not closed, 33 of them at c = 2: consistent with wave 7 Sect. 6 and with audit 8), and the paper's Example `ex:maj` (n = 6, c = 2, four soft triangles) is reported not closed, also not under M-majority. My PC solver (which ignores ternary constraints) is wrong in a handful of those, as it should be.
So Theorem A's hypothesis is exactly the right dividing line in the tests: soft-triangle-free never fails, soft triangles do.

## 3. Hypotheses versus the paper (csp.tex, prelim.tex, app_tunnels.tex)

* ZE: paper = solution with |S| = |H| (app_tunnels, start); by Lemma `decomp` exactly one edit per member and no hard edit; states 1,2,3 = pair ab, bc, ac = the note's L,R,M. **Match.**
* Packing: modification-disjoint (pair sets disjoint, members share <= 1 vertex), each member an induced P3 (ab, bc edges, ac non-edge). Lemma 1 uses exactly this. **Match.**
* Pair-maximality: the paper says (app_tunnels intro) that none of the ZE lemmas except `tfan` needs it. Theorem A does not use it either: an all-hard induced P3 is a static violation, Sigma = empty, which is closed. **Match; the note's "H not necessarily pair-maximal" is right.**
* Soft pair = pair inside some member (paper: `soft` = union of pairs(P_i), owner unique), hard = other. The note's "soft triangle" (three pairs of a triple owned by three distinct members) is **new terminology, not in the paper**; the paper's closest notion is the "triangle if distinct" remark in Lemma `trigid`(a). It must be defined if Theorem A is imported.
* soft triangle vs I(H): a soft triangle gives a triangle of I(H) (owners pairwise share a vertex). Conversely for c <= 2 every triangle of I(H) is a soft triangle (three distinct intersection vertices p,q,r, the pairs pq, qr, pr lie in the three members). For c >= 3 a triangle of I(H) need not be soft (three members through one vertex): so "I(H) triangle-free (any c) => no soft triangle" and "forest => none" are right; the note's "c <= 2 and I(H) triangle-free" is exactly "c <= 2 and no soft triangle". **Correct as stated.**
* The paper's `tmaj` needs only *some* majority; Theorem A supplies all of them, so A2 follows without any change to the paper's lemma. Locks are unit values (paper: "a lock fixes the states ... to unit values"), which is what `tmaj` uses. The Helly-type obstruction for *forked* (2-element) locks that I noticed (three lists {L,R},{L,M},{R,M} defeat any Baker-Pixley Helly argument when all three coordinates carry the same value triple) is correctly left OPEN in wave 7 Sect. 3; nothing in A1-A3 needs it.

## 4. Does polynomial l = 0 contradict NP-hardness at l = 0 (LPS23, Thm 1)? (a10_lps.py)

No. I rebuilt the variable gadget of LPS23 Sect. 4.1.1 from `lps/lps.txt` (three cliques K_0,K_1,K_2 of 5 vertices over F_5, all edges between consecutive cliques, one member v_{0,p} v_{1,q} v_{2,r} for each arithmetic progression q-p = r-q): 25 members, valid modification-disjoint packing, **every cross pair is soft** (this is the design goal of the gadget: "all vertex pairs between consecutive proto-clusters are covered by P3's"), c = 5, and it contains **100 soft triangles** (e.g. pairs of u=v_{0,0}, v=v_{1,0}, w=v_{2,1} are owned by three different members (0,5,10), (4,5,11), (0,8,11)). Because every pair between the cliques is packed, any triple with one vertex in each clique that is not an arithmetic progression is a soft triangle, so the construction cannot avoid them. A1 therefore only says that the hardness must come from ternary (soft-triangle) constraints, which is also what 3-SAT clause structure needs; it is consistent with Theorem 1 and with the paper's `eth` theorem (which lives on those instances). Also consistent with LPS23's open question: at c = 2 hardness for l = 0 would need triangles in I(H).

## 5. Lemma 2 (a10_lemma2.py, a10_local.py)

Identity OPT_{U-only} - |H| = min{|U| : Sigma(G xor U, H) != empty}: OPT_{U-only} computed as min cost over partitions that put **exactly one edit in each member**; right side by enumerating U over unpacked pairs and recomputing Sigma by partition enumeration. 1,166 random instances (n = 5..7, 1..3 members): **identity holds in 1,166/1,166**. OPT < U-only in 18 of them (e.g. n=5, E={02,03,12,23,24}, H={(1,2,4)}: OPT 2, U-only 3), so "NOT exact" is correct and the same phenomenon as the paper's Remark `domain`.
Locality: for random soft-triangle-free instances, for every hard pair e and every triple, I compared the constraint tables (over the owners' states) before and after flipping e: 5,529 changed constraints. Zero violations of: changes only in triples containing both ends of e; scope (owners of the soft pairs) inside the members through the two ends; a binary change on (A,B) needs one end in A, the other in B, and A,B adjacent; **a far pair (no member through p adjacent to a member through r) changes only unary constraints**.
**Fix needed:** 2,187 of the 5,529 changes are *static* triples (all three pairs hard, scope empty). Flipping a hard pair can create or destroy an all-hard induced P3, which empties or un-empties Sigma globally. The sentence "Sigma(G xor U) differs from Sigma(G) only through the lists of members ... and the relations R_AB" must add "and through the static triples (scope empty), which gate non-emptiness". The identity itself already accounts for them (they are part of "Sigma non-empty"). Ternary constraints are untouched (all pairs soft), so nothing to add for soft triangles.

## 6. Strictness Lemma S (a10_strict.py)

Own reconstruction on the 576 graphs (109 distinct triples (R, Lam_A, Lam_B) after de-duplication), all X <= Y subsets: 2,943 pairs, 767 strict, **0 violations of |X'| <= 1**. CORRECT. Remark: the lemma is per tree edge; "propagation is unit propagation along a path" needs that a restricted child list is a singleton at every step, which is true for unit locks (the only ones used in A2/A3) and is what the note says; the many-forked-locks case stays OPEN there. Fine.

## 7. Theorem B, by hand and by brute force (a10_thmB.py)

Hand check of the construction: arms A_ij (ordered by (a,b)), B_ij (ordered by (b,a)); row i = A_{i,k}, ..., A_{i,1}; column j = B_{k,j}, ..., B_{1,j}; the 2k paths are disjoint and every obstruction is a subpath: an arm, or (suffix of A_{i,j+1} with a >= alpha+1) followed directly by (prefix of A_ij with a <= alpha), which is contiguous because A_{i,j+1} immediately precedes A_ij in the row. Links = (tile e in A_ij, tile e in B_ij), a perfect matching. L = k^2.
(=>) a tiling gives k^2 links hitting all arms; rho_alpha is hit iff a_ij <= alpha or a_{i,j+1} >= alpha+1, true for every alpha when a_ij <= a_{i,j+1}. (<=) the 2k^2 arms are pairwise disjoint and an operation meets at most two of them (a link only the two arms of its own cell), so with budget k^2 every operation is a link, one per cell: a tile (same tile in A and B) in every S_ij. If a_ij > a_{i,j+1}, rho for alpha = a_{i,j+1} in [1,n-1] contains neither chosen node, other chosen nodes lie on other arms. Columns identical with b. An empty S_ij gives an empty arm (unhittable): NO, correct. An empty rho (both prefix and suffix empty) forces a_ij > a_{i,j+1} in every choice: NO, correct. Polynomial size; parameter k^2. All steps valid. Nit: the book's Grid Tiling with <= has the inequalities on the other coordinates (a along i, b along j); the note says "after transposing", which is legitimate (rename i<->j, a<->b).
Brute force, own construction and **two solvers that know nothing about the construction** (generic hitting-set branching on a smallest unhit path; and naive enumeration of all L-subsets of operations on tiny instances): k=2,n=3: 300 instances (262 yes/38 no); k=2,n=4: 200 (188/12); k=3,n=3: 150 (69/81); k=2,n=2 with naive enumeration: 150 (116/34); k=2,n=3 with naive enumeration: 60 (26/34). **860 instances, 0 mismatches** (additional to the notes' 850).

## 8. Theorem C

(i) f monotone, a cut at x repairs exactly the cores containing x, cores have <= 2 anchors (A2). A tunnel core is a minimal set of tunnels whose anchors cannot all be locked; its anchors contain a non-face of <= 2 anchors, so it has <= 2 tunnels. OK.
Direct test of the underlying characterisation (a10_tunnel.py): random soft-triangle-free instances (n <= 7), stubs on up to 4 anchors (anchor = in a member, no hard edge), brute force over all partitions of G + stubs, for every subset S of stubbed anchors: "G + stubs(S) has a ZE solution" equals "Sigma non-empty and every single anchor and every pair of S is a face" in **1,602/1,602 checks (0 mismatches)**, and equals the direct face definition (Lemma `tinactive`) in all; 216 non-faces of size >= 3 occurred, **none** was pairwise consistent (no 3-core).
(ii) Needs a rewording. As written ("an arm that must be repaired by cutting a tunnel at any one of its N anchors forces a minimal non-face of size N") the minimality tacitly assumes the rest of the instance is satisfiable. Clean version: by A3 the feasible release sets are exactly the vertex covers of the conflict graph C. If an arm of N >= 3 anchors is repaired by releasing any single one of them (rest of the release set unchanged), then locking all N is infeasible, so some edge xy of C lies in the arm, and releasing a third arm anchor leaves xy uncovered, infeasible: contradiction. This gives the claimed obstruction for N >= 3 rigorously (N <= 2 only).
(iii) "Polarity" is not an independent obstruction: Theorem B's obstructions are themselves **positive** clauses (hit one of N nodes); the incompatibility is (ii) (size N >= 3 positive clauses). Monotonicity of f only excludes negative clauses, which Theorem B does not use. Harmless but mis-sold; suggest deleting or rephrasing.
(iv),(v) remain HEURISTIC as labelled; they are not part of what "Theorem C" proves. The status table is right ("PROVED for soft-triangle-free; other operations HEURISTIC").

## 9. Other statements checked in passing
* A2 "cores <= 2" gives C2 for c <= 2 triangle-free; wave-6 "all-distinct -> M" is then one of the 729 rules. Consistent with audit 8's counterexample (it has 4 soft triangles, rechecked: 4 soft triangles, Sigma = {0000,0202,1111,2020} in 0-indexed states, not closed under M-majority).
* Section 6 of the note (c <= 2 with triangles: 2-decomposable): the 12,568 soft-triangle instances of the second batch found non-2-decomposable Sigma only at c = 3 (34 cases), never at c = 2, and non-closed ones at c = 2 (33 cases), in line with the note.
* The notes' "mod JCC" for A1 is accurate: majority closure of all constraints => strong 3-consistency (path consistency) decides.

## 10. Requested fixes (all minor, none touches a theorem)
1. Lemma 2: add the static triples to the locality sentence.
2. Theorem C(ii): use the vertex-cover form above; drop or rephrase (iii).
3. Theorem A: define "soft triangle" explicitly if imported into the paper; mention orientation invariance.
4. Optional: state in A1 that the hardness gadgets of LPS23 contain soft triangles (100 in one variable gadget), which closes the "no contradiction" point.

## 11. Reproduction (code/, python3 + numpy)
`a10_rel.py` (40 relations, 729 rules), `a10_closure.py SEED COUNT` (e.g. `11 8000`), `a10_tree.py 5 3000`, `a10_lemma2.py 5 1200`, `a10_local.py 1 150`, `a10_strict.py`, `a10_thmB.py 1 300 2 3 0.4` (add any 6th argument for naive enumeration), `a10_tunnel.py 7 1500`, `a10_lps.py`. Library `a10_lib.py`.
