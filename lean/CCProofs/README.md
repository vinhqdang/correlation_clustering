# Machine-checked formulas (Lean 4 + Mathlib)

Formal proofs of the mathematical statements used by the solvers in `ccbench/`.
Only the formulas are verified; the Python/numba code is checked against them by
the tests in `tests/`.

| Lean theorem | Statement | Used in |
|---|---|---|
| `CC.cost_update` | relabelling one vertex changes only the pairs through it | all local moves |
| `CC.cost_move` | moving `v` from `A` to `B` changes the disagreements by `(|B| - 2k_B) - (|A| - 1 - 2k_A)` | `ccbench/anneal.py` (`_anneal`) |
| `CC.costW_move` | weighted move: `s (S_B - S_A + s) + 2 (w(v, A - v) - w(v, B))` | `_anneal_w` in `ccbench/anneal.py` (`_local_move` in `ccbench/localsearch.py` uses a different model with three parameters `w`, `cnt`, `kappa`; its formula reduces to this one for `w = cnt`, `kappa = 1`, and the general case is not covered) |
| `CC.costW_swap` | exchange of two nodes = two moves, with the updated cluster sizes and weights | `_anneal_w` (swap proposals) |
| `CC.twins_split_improvable`, `CC.optimal_keeps_twins` | twins (equal closed neighbourhoods) are together in every optimal clustering | `ccbench/reduce.py` |
| `CC.pivot_round_twins` | a Pivot round never separates twins | seeds in `ccbench/memetic.py` |
| `CC.cost_contract` | the contracted weighted instance has the same cost as the expanded clustering | `ccbench/reduce.py` |
| `CC.d_child`, `CC.cost_child_a`, `CC.cost_child_b`, `CC.crossover_le` | partition crossover: block-wise choice of the better parent is never worse than either parent | `partition_crossover` in `ccbench/memetic.py` |
| `CC.packing_bound` | weak duality for packings of vertex sets with per-pair load at most 1 | triangle / star / subgraph packings in `ccbench/dual.py` |
| `CC.bad_triangle_costIn`, `CC.star_bound` | every clustering pays at least 1 on a bad triangle and at least k - 1 on a star with k independent leaves | the same |
| `CC.lagrangian_weak_duality`, `CC.cost_eq_affine`, `CC.cc_dual_bound`, `CC.sep_triangle` | the Lagrangian bound `K + bᵀy + Σ min(0, w - Aᵀy)` is at most the cost of every clustering, for all rows valid for clusterings (triangle rows are) | `ccbench/blockdual.py` |
| `CC.optimal_separates_far` | vertices that are non-adjacent and have no common neighbour are separated by every optimal clustering | distance-2 support in `ccbench/support.py` |
| `CC.far_dual_bound`, `CC.cost_eq_supp`, `CC.sepFar_of_optimal`, `CC.far_dual_bound_optimal`, `CC.far_dual_le_opt` | on the distance-two support `P` (far pairs fixed apart), `bᵀy + Σ_{p∈P} [min(0, w_p - (Aᵀy)_p) - min(0, w_p)]` is at most the cost of every far-separating clustering, for rows supported on `P` and valid for far-separating clusterings; hence at most OPT | `experiments/check_certificate.py`, `ccbench/blockdual.py` |
| `CC.sep_triangle`, `CC.sep_far_triangle`, `CC.star_row_valid` | row families valid for far-separating clusterings: triangles `x_uw ≤ x_uv + x_vw`, `1 ≤ x_uv + x_vw` for a far pair `(u, w)`, and star rows `Σ_T x_vt - Σ_R x_tt' ≥ |T| - |R| - 1` (`R` = non-far leaf pairs `t < t'`) | `ccbench/blockdual.py`, row checks in `experiments/check_certificate.py` |
| `CC.crossoverP_le`, `CC.costW_eq_pcost`, `CC.crossoverW_le` | partition crossover for every pairwise-additive objective, in particular for the weighted cost `costW` of the contracted instance | `partition_crossover` in `ccbench/memetic.py` |
| `CC.cost_congr`, `CC.cost_comp_injective`, `CC.costW_comp_injective`, `CC.cost_canon` | the cost depends only on the partition; relabelling by an injective map (e.g. child labels `α ⊕ β` into the population's label type) does not change it | `partition_crossover` in `ccbench/memetic.py` (child labels `a` / `b + ka`, then `_compact`) |
| `CC.cost_eq_two_mul_unordered`, `CC.cost_eq_two_mul_card` | `cost` is twice the number of unordered disagreements (pairs `u < v`) | `ccbench/objective.py` |
| `CC.reach_le`, `CC.output_le`, `CC.telescoping_le`, `CC.search_guarantee` | no operator makes the best member worse, so the output costs at most the Pivot seed of every outcome, and E[cost] <= 3 OPT follows from the Pivot guarantee (used as a hypothesis) | `memetic_w`, `local_ils`, `pxmem` |

`cost` counts every unordered pair twice (sum over ordered pairs), so the formulas
appear multiplied by 2.  The 3-approximation of Pivot in expectation (Ailon,
Charikar and Newman, JACM 2008) is a hypothesis of `CC.search_guarantee`; it is not
re-proved here.

## Building

```bash
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
cd lean/CCProofs
lake exe cache get      # prebuilt Mathlib
lake build
lake env lean Axioms.lean   # every theorem depends only on propext, Classical.choice, Quot.sound
```
