# Audit 14: hostile referee report on notes/wave10_hardness_proof.md (W[1]-hardness of CEaMP in l + tw(I(H)))

Scope: wave10_hardness_proof.md, wave8_hardness.md, wave6_hardness.md (skimmed), code w10_*, w8c_*; paper sections prelim/nl/components/open read for definitions and consistency.
No git, paper/ untouched. New code: `code/a14_chain.py`, `a14_thin.py`, `a14_params.py`, `a14_exact2.py`.

## 0. Overall verdict: SOUND WITH FIXES (no gap in the logic found; several statements need correcting or proving in writing)

I tried hardest on item (1), the restricted-model gap. There is none: the notes do not assume tunnel-release (TR) semantics in the soundness direction. Theorem S derives it for arbitrary edit sets (Lemma 5 + decomposition lemma + Lemma 6 + region counting). I re-derived the counting (Sec. 1) and it holds. Independent computation agrees (Sec. 6).

| component | verdict |
|---|---|
| Chain lemma (Lemma 1), all m >= 7 | SOUND (computer-assisted; independently re-verified to m = 100) |
| Profile table (Lemma 2) | SOUND (re-verified for m = 14, 26, 42, 62 with own code) |
| Blow-up (Lemma 3) | SOUND |
| Prop. 4 tunnel-release semantics, OR-switch | SOUND (independent model, 0 mismatches; exact comparator test passes) |
| No-split Lemma 5 | SOUND WITH FIX (statement over-claims, see F1) |
| Region bound Lemma 6, Theorem S, Theorem C | SOUND (re-derived) |
| Hypotheses (valid, modification-disjoint, pair-maximal) | SOUND; pair-maximality proof should be written (F3) |
| Parameters / polynomiality | SOUND WITH FIX (pw and u_max are measured, not proved; F2, F4) |
| Grid Tiling with <= source and encoding | SOUND (theorem number unchecked, F5) |
| Consistency with the paper | SOUND (Sec. 5) |

## 1. Item (1): soundness for arbitrary edit sets

Re-derivation. S any solution, ex(S) <= k^2 = l.
* Decomposition lemma (paper prelim): each member has >= 1 edit, |U| <= ex <= l.
* Lemma 5: if K_v meets >= 2 clusters, the smaller part X (|X| <= p/2) and rest Y give |X||Y| cross pairs. At most u_max per vertex are soft, so >= |X|(|Y|-u_max) hard edges are edited. The bound is concave in |X|, minimised at |X| = 1 (p-1-u_max) or p/2 (p/2 (p/2-u_max)). Both exceed l for p >= k^2+5, u_max <= 3 (checked for k = 1..4: 3 > 1, 28 > 9, 88 > 16). So no K_v is split.
* Regions: R_c, T_c, Pi_c as in the notes. A pair lies in two Pi's only if it is a cross-anchor pair (checked: u in two of R_A+T_A, R_B+T_B forces u in R_A cap T_B or T_A cap R_B). Pad pairs of comparator members lie in some Pi_c but are soft edges inside a whole K_v, hence unedited. All other pairs of H_cmp members (y-pairs, z_A z_B) lie in no Pi, so each H_cmp member contributes one edit to S_rest. Members of H_c are pair-disjoint, so |S cap Pi_c| >= |H_c| + h_c.
* Lemma 6: if |S cap Pi_c| = |H_c| then the partition restricted to R_c + T_c is ZE for G_c. Stub-stub pairs are outside Pi_c but cannot be together (the anchors would then share a cluster across an unedited hard non-edge). All tunnels at c are uncut, so every anchor is locked, and every position is forbidden by a lock or a permanent lock: contradiction. So |S cap Pi_c| >= |H_c|+1.
* Counting: ex >= sum_c max(1,h_c) - q + (|S_rest|-|H_cmp|), with sum h_c = h_in + 2q. This gives ex >= |Z| + h_in + q and |Z| >= 2k^2 - h_in - 2q. Hence ex >= max(q, 2k^2-q) >= k^2. At ex = k^2: q = k^2, Z empty, h_in = 0, h_c = 1, S_rest = exactly one edit per H_cmp member. The single non-member edit per chain is a cross-anchor pair, and it is a tunnel deletion (a cross-anchor insertion with g' != sigma(g) would put two B anchors together via the unedited tunnel at g, across the unedited hard non-edge in R_B). So S = k^2 tunnel deletions (one per cell, since each tunnel joins A_ij to B_ij) + one edit per member, i.e. S is a ZE partition of I - C, and Prop. 4 applies. Thus the "restricted model" is a consequence, not an assumption. I found no hole.

Fixes:
* F1 (wording, harmless). Lemma 5 / Theorem S proof say "no edit lies on a pair inside a super-vertex". False: host pairs (h,i)(h,k) are soft non-edges inside K_h and every ZE partition (z joined to K_h) edits them (state M, the single edit of the host member). The proof only needs "pad pairs of comparator members (soft edges inside a whole K_v) are unedited". Reword.
* F6 (clarify). Lemma 6's "ZE for G_c" must say that comparator pad pairs are hard edges of G_c (they are, since the comparator guests are absent); the notes say "plain edges", fine.

## 2. Item (2): size and parameters (my own measurements, thick instance, `a14_params.py`)

| k | n | p | \|V\| | \|H\| | valid packing | all-hard induced P3 | c | u_max | components of I(H) | pw(I(H)) upper bound (vertex separation, explicit order) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 6 | 550 | 3298 | yes | 0 | 26 = 4p+2 | 2 | 2 | 497 |
| 1 | 3 | 6 | 1006 | 6202 | yes | 0 | 26 | 2 | 2 | 497 |
| 1 | 4 | 6 | 1614 | 10074 | yes | 0 | 26 | 2 | 2 | 491 |
| 1 | 5 | 6 | 2378 | 14958 | yes | 0 | 26 | 2 | 2 | 497 |
| 1 | 2 | 10 | 902 | 8854 | yes | 0 | 42 | 2 | 2 | 1349 |
| 2 | 2 | 10 | 3616 | 35612 | yes | 0 | 42 | 2 | 4 | 5049 (claimed bound 8800) |

* l = k^2 exactly, c = 4p+2, pw independent of n: confirmed. Size polynomial: |V| = O(k^2 n^2 p), |H| = O(k^2 n^2 p^2), p = O(k^2). Parameters are functions of k only: O(k^6) for pw.
* Observation (not in the notes): I(H) has only 2k components (k row groups of A chains linked by comparators, k column groups of B chains). Tunnels are hard pairs and give no I(H) edges. So pw really is O(k p^2) per component. The stated 22 k^2 p^2 is valid but loose.
* F2. "pw <= 22k^2p^2" is a measurement (interval cover) plus a sketch. A proof: sweep cells left to right. Every member lies within 3 consecutive cells of its own chain, guests and comparator members sit at cell 2u+1 and a comparator joins two chains at the same cell. Bag size <= (chains per component) x (members through a cell) = O(k p^2). Add this to the write-up.
* F4. "u_max <= 3 by construction" holds by direct count (host: host-matching + lock guest, or host-matching + two comparator guests = at most 3, equal to 2 for k = 2 because end chains have one neighbour; non-host: at most 2 pads from neighbouring z's / lock / comparator). It is measured 2 (k = 1, 2). Write the count. A bound of p-1 matchings is also needed (p >= 6 fine).

## 3. Item (3): hypotheses

* Valid induced packing, pair-disjoint: verified on all instances above (own code, checks edges ab, bc, non-edge ac, and no pair reused).
* Pair-maximal: 0 all-hard induced P3 for k = 2, n = 2, p = 10 (notes had only k = 1, n = 2, p = 8). Easy proof exists: the hard graph is a disjoint union of (a) K_v glue cliques minus pad/host pairs, (b) joins K_u + K_v along the 4 thin hard edges per chain (a matching), (c) tunnel edges x-x' (anchors have no other hard pair). In every component every non-adjacent pair is a pad/host pair, hence soft, so no all-hard induced P3. Write this (F3).
* The paper's `pair-maximal` is "no induced P3 with all three pairs hard": same notion as tested.

## 4. Item (4): chain lemma

Independent reimplementation (`a14_chain.py`: own validity check, own SAT with possible-together variables only, own member patterns L/R/M, all-solution enumeration). ZE partitions of OH_m are exactly the m partitions with spine word R^a M L^(m-1-a) for every m in {7..60 step 1 for 7..14, then step 3 up to 60, 80, 100}: 0 failures; partitions distinct. The transfer-matrix argument is logically fine: interior transfer matrix is cell-independent by translation invariance, boundary vectors depend on m only via the exponent, M has Jordan block of size 2 at 1 (char poly x^4 (x-1)^2), so N_m is affine for m >= m1 and N_m = m on 7..40 pins it. Caveat for the write-up: state precisely that left cap windows are m-independent and right cap windows depend only on distance to the right end; my direct count to m = 100 is extra evidence, not a proof for all m. The reduction uses m = 2(1+n(n+1)) >= 14 only. Verdict SOUND (computer-assisted certificate).

## 5. Items (5)-(7)

(5) Profile/blow-up generality. Own profile code (from my own enumeration of P_a, no use of w10_prof): for m = 14, 26 (n = 2, 3, in `a14_thin.py` Model assertions) and m = 42, 62 (n = 4, 5 runs): partner sets are constant on classes A0..A3 (A0 empty at u = 0, A3 empty at u = U-1) and the effective forbidden sets of pad / X / Y / all are exactly {2u,2u+1}, {a >= 2u}, {a <= 2u+1}, {} for every group. Lemma 3 is a clean proof; thick-vs-thin verified by the notes plus my exact test below.

(6) Grid Tiling with <=. The standard definition (Marx; used in Cygan et al., confirmed against papers and MIT 6.890/6.892 notes via web search) takes the first coordinate non-decreasing between vertically adjacent cells and the second between horizontally adjacent ones. The notes compare first coordinates horizontally and second coordinates vertically. That is the transposed convention, equivalent by renaming i <-> j. W[1]-hard in k and no f(k) n^{o(k)} under ETH. Encoding: A_ij ordered by (a,b) so block = a, B_ij ordered by (b,a) so block = b; clauses (beta(c1) <= alpha) or (beta(c2) > alpha) for alpha = 1..n-1 are equivalent to beta(c1) <= beta(c2) (if beta(c1) > beta(c2) then alpha = beta(c2) <= n-1 violates). Correct. F5: the notes cite "Thm 14.28" for the <= version; I could not verify the numbering (I believe 14.28 is plain Grid Tiling). Cite as "Grid Tiling with <= [Marx 2007; Cygan et al. Sec. 14.4]" and check the number in the book.

(7) Consistency with the paper.
* XP: 5^{w+1} n^{2l} with l = k^2 and n^l bounds are consistent (GT is solvable in n^{O(k^2)}; lower bound is n^{o(k)}).
* Theorem A / sfA: the instances contain soft triangles (already thin: e_j, e_{j+1}, e_{j+2} via G_j, F_j, F_{j+1}). l = 0 is the 3^w CSP, fine.
* Theorem omega (FPT in (l, D, omega)): here D = 2c = O(p) and l = k^2 are bounded in k, so omega must be unbounded, otherwise FPT = W[1]. It is: Gamma' contains the tunnel edges (hard edges of degree <= B), which join A_ij and B_ij at transposed positions, so tw(Gamma') grows with n. Consistent (the paper itself notes omega is not bounded by (l,w)).
* Theorem K / three (components): I(H) has 2k components, r <= 2l = 2k^2 holds (k = 1: r = 2, l = 1, equality). Each row/column component is tight in isolation (anchors have no hard edge without tunnels) and irreducible (all anchors cannot be singletons), consistent with the counting lemma. Component span is unbounded in n, so Theorem K gives no FPT. Soft-triangle-free closure and Corollary on c = 2: not contradicted (c = 4p+2, l >= 1).
* LOC/repacking: hardness here implies LOC cannot be FPT-solvable unless FPT = W[1]; the paper only lists LOC as open, not as solved.
* Hardness at bounded c remains open; the glue argument needs p > l. Note also "c <= w+1" makes "l + c + pw" redundant with l + pw.

## 6. Independent computation

1. `a14_chain.py`: Lemma 1 for 30+ values of m up to 100 (above).
2. `a14_thin.py`: independent builder of the thin-guest model from the notes' Sec. 2 text (profiles derived from own enumeration, guests, locks, OR-switch, tunnels as assumption literals). ZE feasibility of I - C vs "C is a tiling":
   * k=2, n=2: 6 random instances (dens 0.6) + 30 instances (dens 0.85, up to 256 choices): every choice of one tile per cell, 0 mismatches. Also ALL cut sets of size <= k^2 = 4 (up to 2517 per instance): feasible exactly for the tilings (one cut per cell), none of size < 4, no non-tiling set: 0 mismatches.
   * k=2, n=3 (8 instances, up to 420 choices), k=3, n=2 (6 instances, 18 chains), k=1, n=4, 5: 0 mismatches. (k=3 samples had 0 tilings; the YES side was checked at k = 1, 2.)
3. `a14_exact2.py`: EXACT full-problem test (all 80,200 pairs editable, lazy-triple SAT + cardinality bound, UNSAT is a proof) of the comparator on the thick instance with p = 6: two chains with m = 10 (groups: pad, slot block 1, separator, slot block 2, pad), private stub tunnels, one comparator, |V| = 401, |H| = 2349, c = 26, u_max = 2, l = 2 (p-1-u_max = 3 > 2):
   * tiles (block1 | block2): ([1],[3]) -> SAT, cost |H|+2 (expected YES); ([3],[3]) -> SAT; ([1,3],[1]) -> SAT; ([3],[1]) (block(P)=2 > block(Q)=1) -> UNSAT at l = 2 (expected NO). All four match. This closes the gap that the notes' own row-gadget run (w10_ex2, p = 8) never finished. l = 1 runs were stopped (not claimed).
   * Together with the notes' exact k = 1 runs (6.5), cheaper solutions using hard edits, splits or insertions do not occur on the instances solved. No full k = 2 exact solve was possible (about 3,600 vertices at p = 10); the proof of Theorem S carries that case.

## 7. Remaining caveats for the write-up (none invalidates the reduction)
1. Lemma 1(b) and Lemma 2/profile rest on finite computations; present them as certificates with the translation-invariance argument spelled out.
2. Fix wording F1; add proofs F2 (pw sweep), F3 (pair-maximality), F4 (u_max count); fix citation F5.
3. Lemma 5 needs p >= l + u_max + 2 with u_max = 3 for k >= 3; state p = k^2+5 rounded up to even.
4. The ETH remark (no f(l+w) n^{o((l+w)^{1/6})}) is right only if p = poly(k) and pw = poly(k); fine as stated.
