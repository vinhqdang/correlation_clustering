/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Basic

/-!
# Ordered and unordered counts

`cost` sums the disagreement indicator over ordered pairs.  With a linear order
on the vertices (the vertex numbering of the code), every unordered pair is
written once as `(u, v)` with `u < v`, and `cost` is exactly twice the number of
unordered disagreements (`cost_eq_two_mul_unordered`, `cost_eq_two_mul_card`).
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [LinearOrder V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- `cost` is twice the sum of the disagreements over the pairs `u < v`. -/
theorem cost_eq_two_mul_unordered (c : V → α) :
    cost G c = 2 * ∑ u, ∑ v, if u < v then d G c u v else 0 := by
  have hpt : ∀ u v, d G c u v =
      (if u < v then d G c u v else 0) + (if v < u then d G c v u else 0) := by
    intro u v
    rcases lt_trichotomy u v with h | h | h
    · simp [h, not_lt.mpr h.le]
    · subst h; simp [d_self]
    · simp [h, not_lt.mpr h.le, d_comm G c u v]
  have hswap : (∑ u, ∑ v, if v < u then d G c v u else 0) =
      ∑ u, ∑ v, if u < v then d G c u v else 0 := Finset.sum_comm
  simp only [cost]
  rw [Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => hpt u v))]
  simp only [Finset.sum_add_distrib]
  rw [hswap]
  ring

/-- `cost` is twice the number of unordered disagreements `{u, v}` (written with
`u < v`): positive pairs that are split and negative pairs that are joined. -/
theorem cost_eq_two_mul_card (c : V → α) :
    cost G c = 2 * ((univ.filter fun p : V × V =>
      p.1 < p.2 ∧ ¬ (G.Adj p.1 p.2 ↔ c p.1 = c p.2)).card : ℤ) := by
  rw [cost_eq_two_mul_unordered, Finset.natCast_card_filter, ← Finset.univ_product_univ,
    Finset.sum_product]
  congr 1
  refine Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))
  by_cases h : u < v
  · have hne : u ≠ v := h.ne
    by_cases e : (G.Adj u v ↔ c u = c v) <;> simp [h, hne, e, d, pd]
  · simp [h]

end CC
