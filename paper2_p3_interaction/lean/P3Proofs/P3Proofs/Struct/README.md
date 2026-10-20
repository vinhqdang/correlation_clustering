# Struct: structural lemmas of the P3 interaction paper

Self-contained (imports Mathlib only; definitions in `Defs2.lean`, `Packing.lean`).
Root file: `P3Proofs/Struct.lean`.  Namespace `P3S`.

## Conventions
* Graph `G : SimpleGraph V`, `V` a `Fintype`; a clustering is a labelling `c : V → α`.
* `cost G c` is the sum over **ordered** pairs of disagreements, so `cost = 2 * (paper cost)`.
  `Sol G α k` means "OPT <= k" (some labelling with `cost <= 2k`).  `IsOpt G c`: `c` is optimal among
  all labellings `V → α`.
* Single-vertex moves to a fresh label need an unused label `z` (`hz : ∀ x, c x ≠ z`);
  automatic for `[Infinite α]` (`exists_fresh`; the primed lemma versions).  Twin lemma needs none.
* A packing member is a triple `(a,b,c)` = induced P3 with centre `b` (`IsP3`); `pr P` its six
  ordered pairs; `IsPacking` = induced P3s with pairwise disjoint `pr`; `Soft H x y` = pair in some
  `pr P`; hard = not soft.  Pair-maximal = `PairMaximal`.

## Status
| Lean theorem (file) | Paper statement | Status |
|---|---|---|
| `twin_moves_sum`, `twin_split_improvable`, `optimal_keeps_twins` (Twin) | `lem:twin` (app_structure), strict form: the two moves sum to exactly -2 (here -4 in ordered units) so one is a strict improvement; twins in one cluster in every optimum | proved; no hypothesis on H. `SameClosed u v` = `u ≠ v ∧ N[u]=N[v]` |
| `half_degree`, `half_degree'` (Half) | `lem:half` (a): `d_C(v) >= (|C|-1)/2` in an optimum | proved (as `|C \ v| <= 2 d_C(v)`) |
| `cluster_common_neighbour` (Half) | `lem:half` (b): non-adjacent `a,b` in a cluster have a common neighbour in the cluster | proved |
| `cluster_dist_le_two` (Half) | diameter of a cluster `<= 2` in `G` (`G.dist a b <= 2`) | proved |
| `optimal_separates_far`, `optimal_separates_far'` (Half) | pairs at distance `>= 3` (distinct, non-adjacent, no common neighbour) are separated in every optimum | proved (corollary of the diameter lemma) |
| `IsPacking.insert`, `IsPacking.union` (Packing) | `lem:MAX` (c)/(first claim of `lem:fam`): adding induced P3(s) with all pairs unpacked/pairwise disjoint gives a packing | proved |
| `two_mul_card_le_cost`, `sol_ge_card` (Packing, Max) | `OPT >= |H|` (the lower-bound half of `lem:decomp`) | proved |
| `max_opt_ge`, `hardFam_card_le`, `insert_opt_ge` (Max) | `lem:fam`, `lem:MAX` (a): `OPT >= |H u M|`; `|M| > l` implies NO | proved |
| `max_budget_equiv` (Max) | `lem:MAX` (b): `OPT <= |H|+l  iff  |M| <= l and OPT <= |H u M| + (l-|M|)` | proved (OPT as `Sol`, `l-|M|` in Z) |
| `max_pairMaximal` (Max) | `lem:MAX` (c): `H u M` is pair-maximal for inclusion-maximal `M` | proved |
| `hard_edges_close`, `hard_nonedge_no_common` (Packing) | `lem:HE`/`lem:H2` (a),(b) | proved |
| `hub_matching` (Hub) | `lem:HUB` (a): `nu(Lambda_u) <= l` | proved (matching expressed by recursive predicate `HasMatch`) |
| `hub_cover` (Hub) | `lem:HUB` (b): vertex cover `T`, `|T| <= 2l`, remaining pairs edges/soft | proved |
| `hub_alpha` (Hub) | `lem:HUB` (c): `alpha(G[N_h(u)]) <= 2l+c+1` | proved; `c` is a hypothesis bounding the members through a vertex |
| `hub_clique` (Hub) | `cor:HUB` `|Y| <= 2l+1` | proved |
| `sm_cross_count`, `sm_size_bound` (SmallMerge) | `lem:SM` (ii) | proved, **abstracted**: cross-pair condition (i), `|U+| <= u`, soft degree `<= D` are hypotheses |
| `sm_union` (SmallMerge) | `lem:SM` (iii): `|C| <= max{D+u+1, 2D, 2u}` for a cluster with two pieces | proved, **abstracted** as above |

## Not formalised
* `lem:CP` (cost formula), `lem:decomp` (the equality `|S| = sum |x_i| + |U|`; only the
  inequality `OPT >= |H|` is proved), `lem:deg`.
* `lem:EX`, `lem:conn`, `lem:NF`, `cor:cand` (optimal.tex), `lem:cross` (that pieces share only soft/inserted
  pairs, i.e. SM (i) from the piece definition), `lem:viab`.
* `lem:R`, `lem:exch`, `thm:C1`, `cor:nf`, `prop:pointer` (app_structure.tex).
* `prop:soft`, `prop:B`, and everything in csp.tex (`lem:five` , DPs, theorems A/l0) and locality.tex.
