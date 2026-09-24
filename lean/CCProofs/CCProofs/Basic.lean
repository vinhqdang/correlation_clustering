/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import Mathlib

/-!
# Correlation clustering: objective and local changes

`G` is the positive graph on a finite vertex set `V`; every other pair is negative.
A clustering is a labelling `c : V → α`.  `d G c u v` is `1` when `(u, v)` with
`u ≠ v` is a disagreement (a positive pair that is split or a negative pair that
is joined) and `0` otherwise.  `cost G c` sums `d` over ordered pairs, so it is
twice the number of disagreements.

`cost_update` is the basic locality fact: relabelling one vertex `v` changes only
the terms of the pairs through `v`.  It is proved for an arbitrary symmetric pair
function (`sum_update_pairs`), which also covers the weighted objective.
-/

open Finset

namespace CC

section Generic

variable {U : Type*} [Fintype U] [DecidableEq U] {α : Type*}

/-- Pair functions that are symmetric and vanish on the diagonal change, under a
relabelling of one node `v`, only on the pairs through `v`. -/
theorem sum_update_pairs (F : U → U → α → α → ℤ)
    (hsymm : ∀ i j x y, F i j x y = F j i y x) (hdiag : ∀ i x y, F i i x y = 0)
    (c : U → α) (v : U) (y : α) :
    (∑ i, ∑ j, F i j (Function.update c v y i) (Function.update c v y j)) -
        ∑ i, ∑ j, F i j (c i) (c j) =
      2 * ∑ j, (F v j y (Function.update c v y j) - F v j (c v) (c j)) := by
  set c' := Function.update c v y with hc'
  set δ : U → U → ℤ := fun i j => F i j (c' i) (c' j) - F i j (c i) (c j) with hδ
  have hsplit : ∀ i j, δ i j = (if i = v then δ v j else 0) + (if j = v then δ i v else 0) := by
    intro i j
    by_cases hi : i = v <;> by_cases hj : j = v
    · subst hi; subst hj; simp [hδ, hdiag]
    · subst hi; simp [hj]
    · subst hj; simp [hi]
    · have h1 : c' i = c i := Function.update_of_ne hi y c
      have h2 : c' j = c j := Function.update_of_ne hj y c
      simp [hδ, hi, hj, h1, h2]
  have hsym : ∀ i, δ i v = δ v i := by
    intro i; simp only [hδ]; rw [hsymm i v, hsymm i v]
  have hL : (∑ i, ∑ j, F i j (c' i) (c' j)) - ∑ i, ∑ j, F i j (c i) (c j) =
      ∑ i, ∑ j, δ i j := by
    simp only [hδ, Finset.sum_sub_distrib]
  rw [hL]
  have hcv : c' v = y := Function.update_self v y c
  calc ∑ i, ∑ j, δ i j
      = ∑ i, ∑ j, ((if i = v then δ v j else 0) + (if j = v then δ i v else 0)) := by
        refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => hsplit i j))
    _ = (∑ j, δ v j) + ∑ i, δ i v := by
        simp [Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.sum_ite_irrel]
    _ = 2 * ∑ j, δ v j := by
        rw [show (∑ i, δ i v) = ∑ j, δ v j from Finset.sum_congr rfl (fun i _ => hsym i)]
        ring
    _ = 2 * ∑ j, (F v j y (c' j) - F v j (c v) (c j)) := by
        simp [hδ, hcv]

end Generic

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- Pair function of the objective: `pd G u w x y` is the disagreement of the pair
`(u, w)` when `u` has label `x` and `w` has label `y`. -/
def pd (u w : V) (x y : α) : ℤ :=
  if u = w then 0 else if (G.Adj u w ↔ x = y) then 0 else 1

/-- Disagreement indicator of the pair `(u, v)` under the clustering `c`. -/
def d (c : V → α) (u v : V) : ℤ := pd G u v (c u) (c v)

/-- Twice the number of disagreements of `c`. -/
def cost (c : V → α) : ℤ := ∑ u, ∑ v, d G c u v

omit [Fintype V] in
lemma pd_symm (u w : V) (x y : α) : pd G u w x y = pd G w u y x := by
  unfold pd
  by_cases h : u = w
  · subst h; simp
  · have h' : w ≠ u := Ne.symm h
    simp only [h, h', ite_false, G.adj_comm u w, @eq_comm _ x y]

omit [Fintype V] in
lemma pd_diag (u : V) (x y : α) : pd G u u x y = 0 := by simp [pd]

omit [Fintype V] in
lemma d_comm (c : V → α) (u v : V) : d G c u v = d G c v u := pd_symm G u v (c u) (c v)

omit [Fintype V] in
lemma d_self (c : V → α) (u : V) : d G c u u = 0 := pd_diag G u _ _

/-- Relabelling `v` changes the cost only through the pairs `(v, w)`. -/
theorem cost_update (c : V → α) (v : V) (y : α) :
    cost G (Function.update c v y) - cost G c =
      2 * ∑ w, (d G (Function.update c v y) v w - d G c v w) := by
  have h := sum_update_pairs (pd G) (pd_symm G) (pd_diag G) c v y
  simp only [cost, d]
  rw [h, Function.update_self]

end CC
