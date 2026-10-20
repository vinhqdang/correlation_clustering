# `P3Proofs.Hard`: combinatorial core of the W[1]-hardness reduction

Source: `paper/sections/hardness.tex`, `paper/sections/app_hardness.tex`,
`notes/wave10_hardness_proof.md`, `notes/audit14_hardness.md`, `notes/audit15_hardness_second.md`,
`notes/wave8_hardness.md`, `code/w8c_oh.py`, `code/w10_chain.py`, `code/a14_chain.py`.

**Honest scope.** The graph-level proof (the thick instance with `10^4`-`10^5` vertices, Lemmas
"Blow-up", "Profile of a host", "Tunnel release", "Region bound", Theorem "Soundness" as a statement
about edit sets of `I_tile`, Theorem "Completeness", the pathwidth bound, pair-maximality) is
**NOT formalised**. In particular **Theorem S at graph level is not formalised**: only the arithmetic
of its counting argument is, as an abstract inequality system whose hypotheses (decomposition lemma,
`r_c >= max(1,h_c)`, `sum h_c = h_in + 2q`, the formula for `ex`) are taken for granted. Likewise the
chain lemma is proved only for `m = 7, 8, 9`, not for all `m >= 7`.

Nothing here uses `sorry`, `native_decide` or an extra axiom. Every theorem depends on at most
`propext, Classical.choice, Quot.sound` (checked in `Axioms.lean`, which prints `#print axioms` for
all main theorems during `lake build`). The finite computations use `decide +kernel`, which is
checked by the Lean *kernel* and does not add `Lean.ofReduceBool` (unlike `native_decide`).

## Files

| file | content |
|---|---|
| `GridTiling.lean` | Grid Tiling with `<=`, block encoding, comparator equivalence |
| `Counting.lean` | abstract counting system of Theorem S, tightness |
| `Chain.lean` | thin chain `OH_m`, ZE partitions, finite chain lemma `m = 7,8,9` |
| `Transfer.lean` | transfer-matrix algebra (`x^4 (x-1)^2`, affine growth) |
| `Params.lean` | parameter arithmetic (`p`, `c = 4p+2`, no-split inequality) |
| `All.lean`, `Axioms.lean` | umbrella import, axiom audit |

## Theorem -> paper statement -> status

Status legend: **proved** = Lean theorem, kernel-checked, no sorry, no extra axiom;
**finite check only** = proved for the listed finite cases by kernel computation, not for the general
statement; **not formalised**.

### (1) Grid Tiling and the comparator encoding (`GridTiling.lean`)

| Lean | paper | status |
|---|---|---|
| `GTInstance`, `GTInstance.IsSolution` | Grid Tiling with `<=` (hardness.tex, "Grid Tiling with `<=` is the following problem"): `a_ij <= a_{i,j+1}`, `b_ij <= b_{i+1,j}` (the paper's transposed convention, copied exactly) | definition |
| `solvable_iff_swap_alt` | "We use the version of the cited chapter after exchanging the two coordinates" | proved (exchanging coordinates turns the paper's convention into the other one; I did not re-check the book's convention) |
| `le_iff_forall_clause` | comparator clause family: `beta(c1) <= beta(c2)` iff for every `alpha in [n-1]`, `beta(c1) <= alpha or beta(c2) > alpha` (proof of Prop. "Tunnel release", comparators) | proved |
| `slotGroup`, `sepGroup`, `numGroups`, `slot_lt_sep`, `sep_lt_slot`, `pos_le_sep_iff`, `sep_le_pos_iff` | block encoding: slot group `u = (alpha-1)(n+1)+s`, separator `u_alpha = alpha(n+1)`, `U = 1+n(n+1)`; `a_c <= 2u_alpha-1 iff beta(c) <= alpha`, `a_c >= 2u_alpha+2 iff beta(c) >= alpha+1` (Prop. "Tunnel release") | proved |
| `slotGroup_inj`, `slot_ne_sep`, `group_cover`, `slot_lt_numGroups`, `sep_lt_numGroups` | the groups `1..n(n+1)` are partitioned into slots and separators, encoding injective, all `< U` | proved |
| `posClauses_iff` | the `n-1` position clauses of the OR-switch guests hold iff the blocks compare (`alpha <= beta`) | proved |
| `comparatorConditions_iff`, `solvable_iff_comparators` | last step of Prop. "Tunnel release" (=>): the comparator conditions on the A-chains (rows, order `(a,b)`) and on the B-chains (columns, order `(b,a)`) are equivalent to the Grid Tiling constraints | proved |

Not formalised in this part: that a ZE partition of the graph realises these clauses (Lemmas
"Profile of a host", "Blow-up", Prop. "Tunnel release" at graph level).

### (2) Counting core of soundness (`Counting.lean`)

Abstract data: `h_c, E_c` (`r_c = h_c + E_c`), `q`, `h_in`, `sl = |S_rest| - |H_cmp|`, `ex`, with
hypotheses `sum h_c = h_in + 2q`, `r_c >= 1` (Lemma "Region bound", with `r_c >= max(1,h_c)`),
`ex + q = sum r_c + sl` (Theorem S, display with the corrections F2/F3 of audit 15).

| Lean | paper | status |
|---|---|---|
| `CountSystem.ex_eq` (`ex = h_in + q + sum E + sl`), `q_le_ex` | identity behind Thm S; `ex >= q` (it is implied by the system, no separate "hard edit" hypothesis is needed) | proved |
| `CountSystem.two_k_sq_le_ex_add_q`, `ex_ge_max`, `k_sq_le_ex` | Theorem S: `ex >= max(q, 2k^2-q) >= k^2` | proved (arithmetic only) |
| `CountSystem.tight` | Theorem S, equality case: `ex = k^2` forces `q = k^2`, `h_in = 0`, `|S_rest| = |H_cmp|`, all `E_c = 0` (one edit per member of every `H_c`), all `h_c = 1` | proved (arithmetic only) |
| `tightSystem`, `tightSystem_ex` | the bound is attained by a consistent system | proved |
| `CountSystem.ex_ge_Z`, `Z_ge` | the paper's route: `ex >= |Z| + h_in + q`, `|Z| >= 2k^2 - sum h_c` | proved |
| `printed_inequality_fails`, `CountSystem.ex_ge_Z_E` | **remark on the paper**: the displayed `ex >= |Z| + h_in + q + sum_c E_c` is not valid for all systems (it counts `E_c` twice for `c in Z`, for which `r_c = E_c`; counterexample `k=1`, `h=(0,0)`, `E=(1,1)`, `ex=2` vs right side `4`). The correct form is `ex >= |Z| + h_in + q + sum_{c notin Z} E_c` (`ex_ge_Z_E`). The conclusions of Thm S (equality case) are unaffected: they only use `ex >= |Z|+h_in+q` and `ex = h_in+q+sum E+sl`, and at equality `Z = {}`. The paper text of app_hardness.tex (proof of Thm "Soundness") should be corrected | proved |
| Theorem S at graph level (regions `Pi_c`, `Lemma "Region bound"`, decomposition lemma, cross-anchor pairs are tunnels, application of Prop. "Tunnel release") | Theorem "Soundness" | **not formalised** |

### (3) Chain lemma, finite version, and transfer-matrix algebra (`Chain.lean`, `Transfer.lean`)

| Lean | paper | status |
|---|---|---|
| `PInst`, `OH m` | thin chain `OH_m` of Sec. "The thin chain" / `oh_instance(m)` of `w8c_oh.py` (`3m+2` vertices, `3m` members, 4 hard edges; vertices numbered cell by cell, a renaming of the python numbering) | definition |
| `OH_valid` (`validB`, `degB`) | members are pair-disjoint induced `P_3`'s (`ab, bc` edges, `ac` non-edge), hard edges are not soft, every vertex lies in `<= 4` members | finite check only (`m = 7,8,9`) |
| `PInst.IsZE`, `SamePart` | ZE partition: hard pairs respected, every member has exactly one edit (states `L,R,M`) | definition |
| `isZE_iff` | `IsZE` is equivalent to the constraint-list form used by the search | proved (general `PInst`) |
| `dfs_complete`, `ze_in_dfs` | the restricted-growth DFS with pruning finds every ZE partition (up to renaming of labels) | proved (general `PInst`) |
| `Certificate m`, `certificate_7/8/9`, `chain_classification`, `chain_7/8/9` | Lemma "Chain lemma": the ZE partitions of `OH_m` are exactly `P_0..P_{m-1}`, pairwise distinct, with spine words `R^a M L^(m-1-a)`; here: every `P_a` is ZE; every ZE labelling `f : N -> N` induces on the vertices the same partition as exactly one `P_a`, `a < m`; its spine word is `R^a M L^(m-1-a)` | **finite check only: `m = 7, 8, 9`** (kernel computation `decide +kernel`, no `native_decide`). All `m >= 7`: not formalised (see next lines) |
| `Mt`, `Mt_charpoly` | transfer matrix `M` (6x6), `chi_M = x^4 (x-1)^2` | proved (Mathlib `Matrix.charpoly`) |
| `Mt_ch`, `Mt_pow_affine` | `M^4 (M-1)^2 = 0` and `M^t = M^4 + (t-4)(M^5-M^4)` for `t >= 4` ("`M^t` is an affine function of `t`") | proved |
| `Ncount`, `Ncount_affine`, `Ncount_eq_id` | for any vectors `u, v` and shift `c0`: `N_m = u^T M^(m-c0) v` is affine for `m >= c0+4`; if `N = m` at two consecutive arguments `>= c0+4` then `N_m = m` for all `m >= c0+4` | proved |
| (link graph -> matrix) | "the cell-by-cell DP of `OH_m` has this `M` as transition matrix at every interior cell of every `m`, with `u` independent of `m` and `v_m` depending only on the distance to the right end", and `N_m = m` for `m = 7..40` | **not formalised** (computer output of `w8c_transfer2.py`, `w10_chain.py transfer`). `Transfer.lean` proves only the algebra; this link is the only place where Lemma 1(b) rests on the matrix argument, and `N_m = m` for `m = 10..40` is not checked in Lean (only `m = 7,8,9`, by `Chain.lean`) |
| `P_a` ZE for all `m` (part (a) of the Chain lemma, "window argument", `m <= 60`) | | **not formalised** (finite check only for `m = 7,8,9`) |

Comment on `decide +kernel`: the search `dfs` is evaluated by the kernel (about 100 s for `Chain.lean`).
`dfs_complete` is proved in general, so the finite run is a genuine proof for these `m`; soundness of
the search is not needed because every output is checked against `P_a` explicitly.

### (4) Parameter arithmetic (`Params.lean`)

| Lean | paper | status |
|---|---|---|
| `pOf`, `pOf_even`, `pOf_ge`, `pOf_le`, `pOf_least` | `p` = least even integer `>= k^2+5`, hence `p <= k^2+6` (hardness.tex, Thm statement) | proved |
| `cOf`, `cOf_le` | `c = 4p+2 <= 4k^2+26` | proved |
| `p_ge_l_add`, `no_split_x1`, `no_split_x1_pOf` | Lemma "No super-vertex is split" / Lemma 5: `p - 1 - u > l = k^2` for `u <= 3` (i.e. `p >= l + u_max + 2`) | proved |
| `no_split_all` | full arithmetic of the proof: for `p` even, `p >= k^2+5`, `u <= 3` (`u <= 2` if `k <= 2`), `x + y + u = p`, `1 <= x <= p/2`: `x*y > k^2` (concavity of `x(p-x-u)`) | proved |
| `host_count_le` | `u(v) <= 3` at a host: `1 + 1`, or `1 + 2` | proved (trivial arithmetic; the combinatorial count of matchings per super-vertex is **not formalised**) |
| `c(H) <= 4p+2` as a statement about the graph, the pathwidth bound `pw(I(H)) < 24 k p^2`, pair-maximality | Lemmas "Parameters and validity", "Pathwidth", "Pair-maximality" | **not formalised** |

## Not formalised (overview)

Blow-up lemma, host profile lemma, Prop. "Tunnel release" at graph level, region bound, Theorem S
(graph level), Theorem C, pathwidth bound, pair-maximality, the reduction itself and Theorem H.
Also not formalised: the finiteness/transfer-matrix link for all `m >= 7`.
