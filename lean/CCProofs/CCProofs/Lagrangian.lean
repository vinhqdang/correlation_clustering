/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Basic

/-!
# Lagrangian lower bounds

The block dual ascent certifies `LB(y) = K + bᵀy + ∑ₚ min(0, wₚ - (Aᵀy)ₚ)`.

* `lagrangian_weak_duality`: for any `x ∈ [0, 1]^P` with `A x ≥ b` and any
  `y ≥ 0`, `bᵀy + ∑ₚ min(0, wₚ - (Aᵀy)ₚ) ≤ wᵀx`.
* A clustering `c` is encoded by the separation indicators `sep c u v ∈ {0, 1}`;
  its cost is affine in them (`cost_eq_affine`), and they satisfy every
  triangle inequality (`sep_triangle`).
* `cc_dual_bound`: for rows that every clustering satisfies (triangle rows,
  star and local-subgraph rows), `LB(y) ≤ cost c` for every clustering `c` and
  every `y ≥ 0`; hence `LB(y) ≤ OPT`.
-/

open Finset

namespace CC

theorem lagrangian_weak_duality {P R : Type*} [Fintype P] [Fintype R]
    (A : R → P → ℝ) (b : R → ℝ) (w x : P → ℝ) (y : R → ℝ)
    (hx0 : ∀ p, 0 ≤ x p) (hx1 : ∀ p, x p ≤ 1) (hAx : ∀ r, b r ≤ ∑ p, A r p * x p)
    (hy : ∀ r, 0 ≤ y r) :
    ∑ r, y r * b r + ∑ p, min 0 (w p - ∑ r, y r * A r p) ≤ ∑ p, w p * x p := by
  have h1 : ∀ p, min 0 (w p - ∑ r, y r * A r p) ≤ (w p - ∑ r, y r * A r p) * x p := by
    intro p
    rcases le_or_gt 0 (w p - ∑ r, y r * A r p) with h | h
    · rw [min_eq_left h]; exact mul_nonneg h (hx0 p)
    · rw [min_eq_right h.le]; nlinarith [hx1 p]
  have h2 : ∑ r, y r * b r ≤ ∑ p, (∑ r, y r * A r p) * x p := by
    calc ∑ r, y r * b r ≤ ∑ r, y r * ∑ p, A r p * x p :=
          Finset.sum_le_sum (fun r _ => mul_le_mul_of_nonneg_left (hAx r) (hy r))
      _ = ∑ p, (∑ r, y r * A r p) * x p := by
          simp only [Finset.mul_sum, Finset.sum_mul]
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl (fun p _ => Finset.sum_congr rfl (fun r _ => ?_))
          ring
  have h3 : ∑ p, min 0 (w p - ∑ r, y r * A r p) ≤
      ∑ p, w p * x p - ∑ p, (∑ r, y r * A r p) * x p := by
    rw [← Finset.sum_sub_distrib]
    refine Finset.sum_le_sum (fun p _ => ?_)
    have := h1 p
    linarith [show (w p - ∑ r, y r * A r p) * x p =
      w p * x p - (∑ r, y r * A r p) * x p by ring]
  linarith

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- Separation indicator of the pair `p` under the clustering `c`. -/
def sep (c : V → α) (p : V × V) : ℝ := if c p.1 = c p.2 then 0 else 1

/-- Objective coefficients of the separation variables (ordered pairs). -/
def wcoef (p : V × V) : ℝ := if p.1 = p.2 then 0 else if G.Adj p.1 p.2 then 1 else -1

/-- Constant term: the number of ordered negative pairs. -/
def Kconst : ℝ := ∑ p : V × V, if p.1 ≠ p.2 ∧ ¬ G.Adj p.1 p.2 then 1 else 0

omit [Fintype V] [DecidableEq V] in
lemma sep_triangle (c : V → α) (u v w : V) :
    sep c (u, w) ≤ sep c (u, v) + sep c (v, w) := by
  unfold sep
  by_cases a : c u = c v <;> by_cases b : c v = c w <;> by_cases e : c u = c w <;>
    simp_all

omit [Fintype V] [DecidableEq V] in
lemma sep_mem (c : V → α) (p : V × V) : 0 ≤ sep c p ∧ sep c p ≤ 1 := by
  unfold sep; split_ifs <;> norm_num

/-- The cost is affine in the separation indicators. -/
theorem cost_eq_affine (c : V → α) :
    (cost G c : ℝ) = Kconst G + ∑ p : V × V, wcoef G p * sep c p := by
  have hpt : ∀ u v, (d G c u v : ℝ) =
      (if u ≠ v ∧ ¬ G.Adj u v then 1 else 0) + wcoef G (u, v) * sep c (u, v) := by
    intro u v
    by_cases huv : u = v
    · subst huv; simp [d_self, wcoef]
    · simp only [d, pd, wcoef, sep, huv, ite_false, ne_eq, not_false_eq_true, true_and]
      by_cases ha : G.Adj u v <;> by_cases hc : c u = c v <;> simp [ha, hc]
  simp only [cost, Int.cast_sum, hpt, Finset.sum_add_distrib, Kconst]
  rw [← Finset.univ_product_univ, Finset.sum_product, Finset.sum_product]

/-- Lagrangian lower bound: valid for every clustering, hence for the optimum. -/
theorem cc_dual_bound {R : Type*} [Fintype R] (A : R → V × V → ℝ) (b : R → ℝ)
    (y : R → ℝ) (hy : ∀ r, 0 ≤ y r) (c : V → α)
    (hvalid : ∀ r, b r ≤ ∑ p, A r p * sep c p) :
    Kconst G + ∑ r, y r * b r + ∑ p, min 0 (wcoef G p - ∑ r, y r * A r p) ≤
      (cost G c : ℝ) := by
  have h := lagrangian_weak_duality A b (wcoef G) (sep c) y
    (fun p => (sep_mem c p).1) (fun p => (sep_mem c p).2) hvalid hy
  rw [cost_eq_affine G c]
  linarith

end CC
