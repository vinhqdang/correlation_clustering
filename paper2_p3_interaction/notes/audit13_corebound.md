# Audit 13: wave 9b (Core Bound, leaf test, G_N), hostile referee

Scope: `notes/wave9_corebound.md`, compared with `notes/wave8_branch.md`. The brief said "audit notes/", the folder is `notes/`.
Independent code (written from the note's definitions, imports nothing from w8*/w9*): `code/a13_lib.py`, `a13_dp.c` (exact subset DP, n <= 22, compiled to `a13_dp`), `a13_t1_valid.py` ... `a13_t11_G.py`.
No git, `paper/` untouched. Labels: CORRECT / CORRECT WITH FIX / WRONG.

## Verdicts

| claim | verdict |
|---|---|
| Theorem L (leaf test FPT, radius 3, min-sum CSP in f(w) n^O(1)) | CORRECT WITH FIX |
| Theorem L covers the wave-8 failing triple type | CORRECT, but by changing the definition of "local"; the wave-8 gap is bypassed, not closed |
| Theorem B1 with F' (csp-variant) | CORRECT; F' is NOT the function of wave 8 |
| Lemma G, Theorem A closure, Lemma H, Theorem CB0 (<= 12) | CORRECT WITH FIX (l = 0 only; B1 never branches at l = 0) |
| G_N is a valid input (induced, modification-disjoint, pair-maximal, c = 2, tw 2) | CORRECT |
| Theorem CE (T_N minimal, size 2N+4, for all N) | CORRECT WITH FIX |
| Lemma R5 (no conservative NU of any arity) | CORRECT WITH FIX (needs c >= 3; remark on c = 2 over-reaches) |
| Minimum-certificate remark (min = 3, B1 with a clever certificate not refuted) | CORRECT WITH FIX |

## 1. G_N is a valid input (`a13_t1_valid.py`)
Built from Sec. 5.1 only. For N = 1, 2, 3 (n = 9, 15, 21): every member is an induced P3 with the stated centre; the 9N member pairs are pairwise distinct (modification-disjoint); there is **no** induced P3 of G whose three pairs all lie outside the member pairs (pair-maximal, all triples tested); c = 2; |E(I(H))| = 4N-1 (N triangles joined by N-1 edges, tw = 2). Exact OPT by subset DP: OPT(G_N) = |H| + N for N = 1, 2, 3 (n = 21 took 28 s), as Sec. 6 says. The pendants w_{N+1} (on a_1) and w_{N+2} (on b_N) cost nothing because A_1 = R, B_N = R detach them; that is why OPT is |H|+N and not |H|+N+2.

## 2. Is F' the function used in wave 8's B1? No, and that is fine
Wave 8's F: partition variant, local := X_a cap X_b nonempty in a tree decomposition of I(H), far pairs outside T are *free but still transitive*.
Wave 9's F': csp-variant, local := members through the ends equal or adjacent in I(H), triples containing a far pair outside T are *dropped*.
Hence F' <= F (every partition is an F'-assignment) and the two N differ. Lemma R (a)-(c) and Theorem B1 hold for any split of the hard pairs into local and far and for any such relaxation, so B1 holds verbatim for F'. I re-derived Lemma R(b),(c) for F' (the non-empty-subset condition on a member is implied by the member's own triple, which has no far pair). `a13_t8_b1.py`: my own B1 (F' by CP-SAT over pair variables, certificates by randomly ordered deletion, recursion on G xor e) against exact OPT, random graphs n <= 10, random non-maximal packings, l = 0..2: 237 (instance, budget) pairs over two seeds, 0 mismatches.
Consequence for the quoted claims: a minimal certificate of F' is a certificate of F, not conversely, so Theorem CE must be (and was) checked for both. See Sec. 5.

## 3. Theorem L
Proof re-read line by line. Triple-scope lemma: correct (one pair: radius 2; two pairs share a vertex: radius 3 around B; three pairs reduce to two). Bag sets: {t : m in N^3[beta_t]} is a union of subtrees T_y chained along paths of I(H), hence a subtree; the bit of a local pair lives on an intersection of subtrees, nonempty by the lemma; |B_t| <= (w+1)(1+D+D^2+D^3), D <= 3(c-1) <= 3w. CSP DP is f(w) n^O(1) (f doubly exponential). Checks of mine (`a13_t7_L.py`):
* explicit instance: 6 members in a path A'-A-B-B'-C-C', u in A cap A', v in B cap B', w in C cap C'; uv, vw local, uw far; W(tau) is all six members and the minimal covering radius is exactly 3 (eccentricities 5,4,3,3,4,5). Radius 3 is attained for the lemma. Random tree-like packings: maximum radius 3, never 4.
* primal graph of the leaf-test CSP on G_N (members + local-pair bits, one clique per triple with a non-far pair): min-fill width 11 for N = 3, 6, 10, 14, constant.
FIX needed. (i) The claim "the failing triple type is an artefact" is right, but the cure is a new definition of local. With wave 8's X-locality, local pairs only need M(u) cup M(v) inside one bag; the members can be arbitrarily far in I(H), so the radius argument does not apply and the wave-8 leaf test stays unproved. The note should say "bypassed by changing N", and say that B1 is then run with N defined by member adjacency. (ii) Radius 3 being attained shows tightness of the ball construction, not that no other decomposition could do better. (iii) Pairs with an end outside V(H) are far, so |N| is Theta(n^2); harmless.

## 4. Theorem CB0 (l = 0, soft-triangle-free)
* Lemma G: classification re-checked by brute force on random packings (`a13_t11_G.py`, 300 instances): every triple with a far pair has (#far, #soft non-far, #local-hard non-far) in {(3,0,0), (2,0,1), (2,1,0), (1,0,2), (1,1,1)}; never two soft non-far pairs, never one soft pair with fewer than two far pairs unless a local hard pair is also present. Matches (i)-(iii).
* Theorem A at graph level (`a13_t6_maj.py`): two members sharing exactly one vertex, all 36 orientation pairs x 16 cross statuses, relation = conjunction of all triples on the 5 vertices: 32 distinct relations (one is empty), **all 729 uniform majorities preserve all of them**, and all lists. In the soft-triangle-free case every cross pair is hard (a soft cross pair p r with p v, v r soft is a soft triangle), so these are all the binary scopes; every other triple gives a list or a static test. I did not re-prove wave 7 Lemma 1 itself; the case analysis above replaces it.
* Lemma H proof checked by hand (the failing case is exactly three 2-lists with pairwise distinct missing values; a fourth coordinate always breaks it). Random test: 294 closed relations of arity 4, 5, 10,647 minimal inconsistent list systems, max 2 coordinates, 0 above 3.
* CB0 counting (<= 3 coordinates, <= 2 guards each unless one coordinate empties in <= 3 guards, <= 2 far pairs per guard => 12) is correct. 12 is not claimed tight.
FIX: (a) at l = 0 B1 does not branch (accept iff F'(N) <= |H|), so CB0, like CE, is a statement about the *structure* of l = 0 certificates; its algorithmic value is only through an l >= 1 version, which is Sec. 8 (open, Erdos-Posa sketch). The note's Verdict 2 reads as if B1 were rescued; it is not. (b) The static-violation bound "<= 3" needs the three pairs of one triple only; fine.

## 5. Theorem CE
Everything below is my own code (`a13_t2_cert.py`, `a13_t3_sigma.py`, `a13_t4_union.py`, `a13_t5_xloc.py`, `a13_t10_conn.py`).
* F'(empty) = |H|, F'(T_N) = |H|+1, F'(T_N - e) = |H| for every e: N = 1..6, 8, 10 (CP-SAT, pair variables, definition as in Sec. 1). Monotonicity makes the single-deletion test sufficient for minimality. T_N consists of far pairs (all checked).
* Partition variant with member-adjacency locality, exact subset DP: N = 1, 2: F(T_N) = |H|+1, F(T_N - e) = |H| for all e. (My DP validated against CP-SAT on 30 random graphs; F >= F' on 30 random splits.)
* Sigma_0 structure lemma: brute force over all 3^{3N} state vectors, N = 1..4: |Sigma_0| = 3, 8, 21, 55 and equal to the claimed set {s_j in {LRL,RLR,RRM}, not (B_j = L and A_{j+1} = L)}. Scope analysis: the only constraint scopes meeting two blocks are {B_j, A_{j+1}}; active-triple count is 28N-8; interior blocks are isomorphic, so the lemma holds for all N (the note's "proof sketch" is adequate with this isomorphism remark; I add it). The inductive transfer argument (A_1 = R, C_j != M => RLR => B_1 = L => A_2 != L => ... => B_N = L, against B_N = R) and the three minimality constructions (all LRL; all RLR; RLR..RLR RRM LRL..LRL) are correct, checked by hand.
* Wave-8 locality (X-sets of the path decomposition with bags {A_j,B_j,C_j}, {B_j,A_{j+1}}), partition F: T_N alone is not a certificate (F = |H| for N >= 2), but T_N plus the "link" pairs is, and greedy minimalisation gives minimal certificates of size 9, 12, 15, 18 for N = 2..5 = 3N+3, as claimed. So the unboundedness does not depend on which locality is used.
FIXES.
1. The theorem proves F'(T_N) >= |H|+1; the "= |H|+1" in the abstract is verified (N <= 10) but not proved. It is not needed.
2. At l = 0, T_N is never branched on. The algorithmic consequence is the l >= 1 statement of Sec. 8, which the note leaves unverified. I verified it (`a13_t4_union.py`): k = l+1 disjoint copies of G_N (valid, c = 2, tw 2): T = union of the T_N is a minimal certificate for budget |H|+l of size (l+1)(2N+4) (N = 1, 2, 3 with k = 2; N = 1 with k = 3). Disjoint copies can be removed by solving connected components separately, so I built a **connected** family (`a13_t10_conn.py`): chain G_{2N} plus one extra pendant w_m on the glue vertex b_N, far pairs w_m b_N, w_m r_N, w_m p_{N+1}; T = T_{2N} + these three, |T| = 4N+7, F'(T)-|H| = 2 > 1, every F'(T - e) - |H| = 1 (N = 1..4; valid, pair-maximal, connected, c = 2). So "B1 with an arbitrary minimal certificate has unbounded branching" is true already at l = 1 for connected instances. Add this to the note.
3. Nothing is claimed beyond B1: CE does not touch the complexity of CEaMP (G_N has OPT = |H|+N).

## 6. Lemma R5
The relation Rk = {000,002,011,020,101,110,111,121,200,220} (L,R,M = 0,1,2) is realised, in my reconstruction (`a13_t9_R5.py`), by members A=(2,1,0), B=(4,3,1), C=(5,0,3), hard non-edges 25, 45, and four outside members, each owning one "dropped" pair (24, 23, 40, 51); instance valid, c = 3. The conjunction of the triples not touching an outside member's pair is exactly Rk. The proof (union-closure of N_0 = {B : f(1_B) = 0} on R5 = {x+y+z != 1}, NU forces all singletons into N_0, contradiction with idempotence) is correct; brute force confirms no Boolean conservative NU preserving R5 for r = 3, 4.
FIX: Rk needs c >= 3 (vertex 2 lies in three members). The "c = 2 language has an NU of arity 4 but no plurality NU" part is computer-only and I did not re-run it. The sentence "consistent with Theorem CE" is overstated: CE is a c = 2 phenomenon, where R5 is unavailable; what explains CE is that an NU of arity m only bounds list cores when the other m-1 values at a coordinate coincide (singleton lists), which fails for 2-element lists. I checked that reasoning by hand.

## 7. Minimum-certificate remark
True in the family: a pendant w on the centre p_j of A_j with wq_j, wa_j far non-edges gives {wp_j, wq_j, wa_j}, forcing both pairs at the centre to be edited (excess >= 1); l+1 pendants at distinct centres give a certificate of size 3(l+1), and G_N has N such pendants. FIX: the lemma needs wy and wy' to be **far** (if one is local it is a counted pair and the argument changes); it is a property of the instance design, not a bound on the minimum certificate in general, and computing a minimum certificate is itself not shown FPT. So "B1 with a cleverly chosen certificate is not refuted" is right and the question stays OPEN.

## 8. Not re-verified (taken from the note)
Counts of K3/diamond/link relations (163/445/39), uniqueness of the conservative majority, NU4 at c = 2, K4 results, kappa_1 searches, the 280-instance leaf-test and 449-pair B1 runs (replaced by my own 237-pair B1 run), Sec. 8 Erdos-Posa sketch (labelled open by the note, still unproved).

## 9. Reproduction (from `code/`)
`python3 a13_t1_valid.py` (about 30 s), `a13_t2_cert.py 6`, `a13_t3_sigma.py`, `a13_t4_union.py`, `a13_t5_xloc.py 5` (several minutes), `a13_t6_maj.py`, `a13_t7_L.py`, `a13_t8_b1.py 1 40` and `7 40`, `a13_t9_R5.py`, `a13_t10_conn.py`, `a13_t11_G.py` (noisy: prints validate warnings). Build the DP first: `gcc -O2 -o a13_dp a13_dp.c`.
