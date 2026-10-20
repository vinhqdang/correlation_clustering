# Audit 12: hostile referee report on notes/wave9_fpt_forest.md

Scope: wave9_fpt_forest.md (code w9_*.py, w9a_*.py, w9_eng.c). `paper/` not touched, no git. Independent code: `code/a12_lib.py` (exact OPT by subset DP, cluster test, packings),
`a12_t1.py` (claim 1), `a12_t2.py`/`a12_t2b.py` (relations/majorities), `a12_t3.py`/`a12_t3b.py` (Helly-2, Theorem U), `a12_t4.py` (RM2), `a12_t5.py` (hub example),
`a12_t6.py` (core step of Theorem F), `a12_t7.py`/`a12_t7b.py` (claim 7). None imports any w9*/a10/a11 file. All instances n <= 11.

## Verdicts

| # | claim | verdict |
|---|---|---|
| 1 | 5 member states instead of 7; 5^{w+1} DP; compatible with the paper's Theorem A | CORRECT (no contradiction; same objects, 5 is a sub-domain) |
| 2 | Helly number 2 for unit lists over a majority-closed base (Lemma H2) | CORRECT |
| 3 | Theorem U: far-lift zero-excess model in 6^l n^O(1) on soft-triangle-free instances | CORRECT (the model is a restricted one; see caveats) |
| 4 | RM2 is FPT, (2d+1)^L, trivial deepest-LCA search tree | CORRECT |
| 5 | n=11 example: exact OPT = |H|+3, RM2 predicts 2, so RM2 unfaithful | CORRECT |
| 6 | Theorem F (conditional on BC); l <= 1 polynomial | CORRECT WITH FIX (wording/definitions; BC is not used for correctness at all; l<=1 is trivially polynomial) |
| 7 | a hard vertex never excludes O unless it also excludes L or R | CORRECT (for unary lists of an outside vertex, and also for binary slices); the "consequences" drawn from it are interpretation |

## 1. Five states vs seven (CORRECT)

Meaning. The paper's variable is the edit set x_i, a non-empty subset of pairs(P_i) (7 values). Edit set <-> resulting graph on {a,b,c} is a bijection, so the notes' "set partition of
{a,b,c}" is the same object seen through the result. Exhaustive check on 3 vertices (`a12_t1.py` header run): of the 7 non-empty edit sets, {ab,ac} and {bc,ac} leave an induced P3 on the
member's own triple (masks 5, 6 not cluster), the other five are cluster graphs: ab=L, bc=R, ac=M, {ab,bc}=O, {ab,bc,ac}=X. The paper's constraint C_t for t = V(P_i) (scope {i}) rejects
exactly the two bad values, so deleting them from the domain changes no solution and no cost. Hence the DP of Theorem A (and the bases 7^{w+1}, 7^{w+l+1} of Theorems A and nl) can be read
with base 5. It is a valid improvement, not a different model. The paper's Remark (domain {1,2,3} does not suffice) is consistent: its two examples are O (two edits) and X (three edits).
Check against exact OPT (`a12_t1.py`, seeds 1,2,3; 74 random instances n=5..8, |H| <= 4, maximal and non-maximal packings, 0 mismatches): OPT = min over hard sets W with |W| <= OPT-|H|
of |W| + min over x in {5 states}^H, identical to the 7-state value and to the subset-DP OPT. In all 195 optimal partitions enumerated, every member is edited, the member edit sets that occur
are in {ab, bc, ac, {ab,bc}, {ab,bc,ac}} (seed 3 contains X), and |W| + sum|x_i| = OPT. Note: the table of lists L(S) (Section 1 of the notes) was re-derived by brute force over a 4-vertex graph: identical
(`a12_t7.py`). Minor: the notes' sentence "5, NOT 7" reads like a conflict with the paper; it is not one, only the paper's phrase "domain size 7" can be tightened to 5.

## 2. Helly-2 for unit lists (CORRECT)

Proof (minimal infeasible family J, |J| >= 3, tuples t_j satisfying J minus {j}, u = m(t_a,t_b,t_c)) is the standard argument and is right: for j outside {a,b,c} all three have value s_j, for j = a two of the
three have s_a, and m(x,y,y)=y. It needs only one majority operation m applied coordinatewise and closure of the base set Sigma under m; the domain size is irrelevant.
Base closure re-derived independently. A single cross-pair triple gives a relation that is NOT majority-closed (14 of 18 such relations violate many of the 729 majorities: `a12_t2.py`); what is
closed is the conjunction over the four cross pairs at a shared vertex: I recomputed all 3x3x16 configurations from the formula (`a12_t2b.py`): 40 distinct relations on {L,R,M}^2, sizes {0,1,2,3,5},
0 violations under all 729 majorities. This matches wave 7 Thm A; the sentence in wave9 Sec. 2 must be read as "the 40 conjunction relations", not the single triples. On the 5-state domain: 99 conjunction relations,
no median (for any of the 120 linear orders) closes them all, consistent with the notes' remark that the extended domain is not majority closed (I did not test non-median majorities).
Instance-level (`a12_t3.py helly`, 4 seeds, ~15,000 random soft-triangle-free instances, n 6..10, |H| <= 5, cyclic I(H) allowed; I also asserted that every triple carrying a soft and a far pair
is a single-soft-pair triple of a liftable group (A,t)): among instances where Sigma_perm is feasible but the liftable lists kill it, the minimal infeasible list family had size 1 (4,621 times) or 2 (1,265 times), never 3.

## 3. Theorem U (CORRECT)

Proof read. Soundness: trivial. Completeness: a far edit changes no relation and no permanent unary list (cross pair statuses are untouched), only liftable lists and static triples; so any solution must hit a
pair {a,t} of one of the <= 2 lists of a minimal infeasible family (Lemma H2); the branching is <= 3 (static P3) or <= 6 (two lists x 3 pairs); depth <= l; 6^l. One check the notes skip: a far pair {a,t} with t in V(B)
changes both (A,t) and (B,a); this is covered because branching is on the pair itself. For general width the CSP is binary + unary and majority-closed, so path consistency (Jeavons-Cohen-Cooper) decides it:
polynomial, no dependence on w. Correct.
Independent test (`a12_t3.py U`, `a12_t3b.py`): my own implementation of the algorithm (lists computed from the triples of the current graph, greedy minimal family, branching on far pairs of the family) versus brute force
over all U subset of Far with the ZE test `G xor U xor (L/R/M edits)` is a cluster graph: l=2: 900 instances (n 6..11, |H| <= 4/5, I(H) with cycles in 84 of 400), 122 of the last 400 answer YES only with |U| >= 1, 185 NO, 0 mismatches;
l=3: 150 instances, 0 mismatches (a further l=3 run of 400 instances, `a12_t3b.py 12 400 3`, see the addendum at the end). Largest minimal family seen 2, largest branching 6.
Caveats (not errors): (a) the question answered is FLZ, a restricted model (no breaks, no cross-pair edits); Theorem U is FPT for FLZ only, it does not decide CEaMP; the notes say so. (b) Static all-hard P3 that contain a cross pair are
branched on (including the cross pair) rather than left to F: harmless (soundness unaffected, completeness unaffected), but contradicts the "U subset of Far" formulation; state it.

## 4. RM2 FPT (CORRECT)

Exchange argument verified: Q unhit conflict with deepest LCA v; if some solution hits Q by a break at w in Q, every other unhit conflict R through w has top node an ancestor of w with depth <= depth(v), hence an ancestor-or-equal
of v on the root path of w, so R passes through v (R contains the whole path from w to its top). Tunnel case: the tunnel has an end that is an endpoint of Q: <= 2d branches. Total 1+2d, depth L. Right.
Independent test (`a12_t4.py`): random rooted forests (m <= 10), <= 7 conflicts, tunnel graph with max degree d in {1,2,3}: minimum number of operations by exhaustive search equals the least L at which the (2d+1)-way search tree succeeds
in 5,000 of 5,000 instances. The packing-endpoint lemma of wave 7 is indeed not needed for RM2 (but RM2 is only a model, see 5).

## 5. The n=11 example (CORRECT)

`a12_t5.py`: edges 01,02,04,05,06,13,28,2-10,34,37,38,78 (12 edges). H = {(2,0,1),(4,3,1),(5,0,6),(2,8,7)}: all four are induced P3, 12 pairs pairwise disjoint, 22 induced P3 in G, every one has a soft pair (pair-maximal),
I(H) = star with centre (2,0,1), no soft triangle. **Exact OPT by subset DP over all partitions = 7 = |H|+3** (confirmed). Search over all hard-edit sets of size <= 2 with all 5-state vectors: minimum 7 (no hard edit helps).
Cheapest state vectors with no hard edits: (L,X,O,L) and (O,O,O,L), both cost 7. Each of M1=(4,3,1), M2=(5,0,6) cannot be zero-excess even when all other members take arbitrary states (outside vertices 7,8 with N cap V(M1) = {3} (the centre), resp. 4 with N cap V(M2) = {0}: list {O,X}),
while M0 and M3 can. Individually each is repaired by one break; the pair of breaks forces M0 to break (the only cheapest vectors with M1=M2=O have M0=O). So RM2 (breaks as node hits, independent regions {M1},{M2}) gives excess 2, the true excess is 3.
This falsifies exactly the claim (L1) of wave 7 Sec. 5 (a break behaves as "hit every conflict path through the node", with no side effects on neighbours), and hence "CEaMP at w <= 1 is RM2". Remarks: H is not a maximum packing of this G (greedy search finds 6 pairwise modification-disjoint induced P3), which is allowed (CEaMP takes any packing) but means l = 3 is measured against |H|, not against a maximum packing.

## 6. Theorem F and the l <= 1 claim (CORRECT WITH FIX)

Proof checked step by step.
* Soundness: every call evaluates real solutions by Lemma E. Fine.
* Completeness. (i) Take an overall optimum; its remainder over U is optimal for G xor U (otherwise improve the whole), so Cor. CAND applies to it, as the notes need; choose W minimising |W|+Ex(U cup W); if Ex(U) > r then W is non-empty. (ii) Z(W) must meet C: constraints not touched by W are identical in U and U cup W, so Ex_{complement of Z(W)}(U) <= Ex(U cup W) <= r-1; if C were inside the complement, monotonicity (fewer constraints, smaller value) would give Ex_C <= r-1, contradicting the choice of C. Right. (iii) Z(W) meets C in a member or an edge at a member containing an endpoint of a pair of W, so that pair has an endpoint in V(C). (iv) The partner is in the radius-2 ball by CAND on the current graph (pair-maximality of H in G xor U holds since static all-hard P3 are branched away first). (v) Recurse with W minus the pair: |W'|+Ex <= r-1. All consistent.
* The role of BC. Correctness does NOT use BC; only the running time does (|C| <= K). The notes' phrase "BC implies correctness" is therefore a misreading; what is true: the algorithm is correct unconditionally (given the oracle), FPT iff the minimal sets C have bounded size. This should be stated so, since BC(t) is an assumption about every inclusion-minimal C reachable, for every U. The bound also needs D <= 4 (true: D(H) = 2c <= 4 for w <= 1) and the bag/bounded-degree argument for the oracle (degree of I(H) <= 3(c-1) = 3).
* Fix 1: define "sites", V(C) and the deletion rule precisely: which static triples belong to which edge when a triple contains cross pairs of different edges (the notes give no rule). I tested the following rule: a triple with a far pair is kept iff one of its associated sites (member owning the soft pair, or edge of a cross pair) is in C. `a12_t6.py` checks the heart of (ii)-(iii): for U = empty, instances with Ex > r (r = 1, 2), minimal C by deletion, every nonempty far W with |W|+Ex(W) <= r has an endpoint in V(C): 288 eligible instances in total (r=1: 127 with 99 useful W; r=2: 161 with 621 useful W), 0 failures; |C| seen 0..2. This tests the combinatorial core, not BC, not the DP oracle, and not the full recursion (I did not re-run the end-to-end loop of w9a_condalg; its 310 instances remain the notes' own evidence).
* Fix 2: Obs. 5.3 (N^1-bag oracle) is only HEURISTIC, so "PROVED implication" holds modulo the oracle; the notes label this correctly in the status table but Sec. 0 item 5 calls the implication written out; keep the label.
* Fix 3: the cross-pair branching in static P3 (see 3(b)).
* l <= 1 polynomial: TRUE and trivial. For fixed w the paper's Theorem A already gives 5^{w+1} n^{2l+O(1)}, polynomial for l <= 1 (and n^{O(1)} for any constant l); so "already known" is right and the "unconditional last level" adds only a speed-up (Theorem U step 5 at r=1: a minimal C has <= 2 members by H2, which I confirmed above). Do not present it as progress on FPT.
Open point re-confirmed: BC(t), t >= 1, is not proved, and no hardness is proved. Verdict on the claim as a conditional statement: CORRECT WITH FIX (wording and definitions).

## 7. "A hard vertex can never exclude O unless it also excludes L or R" (CORRECT)

Brute force on the 4-vertex graph (`a12_t7.py`): the 8 lists L(S) are {all 5}, {L,O}, {O,X}, {O,R}, {R}, {X}, {L}, {M}. Every list excluding O is {R}, {X}, {L} or {M}, which excludes L or R (or both). Reason: O has all final statuses 0 and a single-triple constraint excludes it only by forcing a status to 1, and y_ab = 1 kills L, y_bc = 1 kills R, y_ac = 1 kills both. Extended: all 14 distinct lists obtainable from one outside vertex using any subset of its three triples, and all intersections of up to 3 of them: 0 cases exclude O but keep both L and R. Also for the binary relations (99 conjunction relations, each slice for a fixed neighbour state): 0 cases (`a12_t7b.py`).
Limits: this is about constraints from vertices outside the member and (as checked) from neighbouring members given their state; it does not say that breaks cannot be forced (the intersection {L,O} cap {O,R} = {O} forces O), and the consequences drawn in Prop. 4.4 / Sec. 0 item 3(iii) ("kills the Grid-Tiling blueprint", "no link can reward breaks") are interpretation and should be labelled HEURISTIC; the n=11 example shows that breaks do interact through locks at a hub.

## Addendum: runs pending when this file was written
`a12_t3b.py 12 400 3` (Theorem U with l=3): see the last line below if present.
