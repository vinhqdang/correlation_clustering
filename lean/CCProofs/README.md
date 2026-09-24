# Machine-checked formulas (Lean 4 + Mathlib)

Formal proofs of the mathematical statements used by the solvers in `ccbench/`.
Only the formulas are verified; the Python/numba code is checked against them by
the tests in `tests/`.

| Lean theorem | Statement | Used in |
|---|---|---|
| `CC.cost_update` | relabelling one vertex changes only the pairs through it | all local moves |
| `CC.cost_move` | moving `v` from `A` to `B` changes the disagreements by `(|B| - 2k_B) - (|A| - 1 - 2k_A)` | `ccbench/anneal.py` (`_anneal`) |
| `CC.costW_move` | weighted move: `s (S_B - S_A + s) + 2 (w(v, A - v) - w(v, B))` | `_anneal_w`, `_local_move` |
| `CC.costW_swap` | exchange of two nodes = two moves, with the updated cluster sizes and weights | `_anneal_w` (swap proposals) |
| `CC.twins_split_improvable`, `CC.optimal_keeps_twins` | twins (equal closed neighbourhoods) are together in every optimal clustering | `ccbench/reduce.py` |
| `CC.pivot_round_twins` | a Pivot round never separates twins | seeds in `ccbench/memetic.py` |
| `CC.cost_contract` | the contracted weighted instance has the same cost as the expanded clustering | `ccbench/reduce.py` |
| `CC.d_child`, `CC.cost_child_a`, `CC.cost_child_b`, `CC.crossover_le` | partition crossover: block-wise choice of the better parent is never worse than either parent | `partition_crossover` in `ccbench/memetic.py` |

`cost` counts every unordered pair twice (sum over ordered pairs), so the formulas
appear multiplied by 2.

## Building

```bash
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
cd lean/CCProofs
lake exe cache get      # prebuilt Mathlib
lake build
lake env lean Axioms.lean   # every theorem depends only on propext, Classical.choice, Quot.sound
```
