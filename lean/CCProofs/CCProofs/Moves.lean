/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Basic

/-!
# Move and swap formulas

* `cost_move`: moving a vertex `v` from its cluster `A` to a cluster `B` changes
  the number of disagreements by `(|B| - 2 k_B) - (|A| - 1 - 2 k_A)`, where `k_X` is
  the number of positive neighbours of `v` in `X`.
* Weighted instances (`WInst`): node `i` stands for `s i` vertices and `w i j`
  counts the positive pairs between the vertex sets of `i` and `j`.  Moving node
  `v` (size `s`) from `A` to `B` changes the cost by
  `s (S_B - S_A + s) + 2 (w(v, A - v) - w(v, B))` (`costW_move`), and exchanging
  `v ∈ A` with `u ∈ B` changes it by the sum of two such terms evaluated before and
  after the first move (`costW_swap`).  These are the formulas used by the
  annealing and local search code.
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- Change of `2 * (number of disagreements)` when `v` moves to the cluster with label `y`. -/
theorem cost_move (c : V → α) (v : V) (y : α) (hy : y ≠ c v) :
    cost G (Function.update c v y) - cost G c =
      2 * ((((univ.filter fun w => c w = y).card : ℤ) -
              2 * (univ.filter fun w => c w = y ∧ G.Adj v w).card) -
           (((univ.filter fun w => w ≠ v ∧ c w = c v).card : ℤ) -
              2 * (univ.filter fun w => c w = c v ∧ G.Adj v w).card)) := by
  have key : ∀ w, d G (Function.update c v y) v w - d G c v w =
      ((if c w = y then 1 else 0) - 2 * (if c w = y ∧ G.Adj v w then 1 else 0)) -
        ((if w ≠ v ∧ c w = c v then 1 else 0) -
          2 * (if c w = c v ∧ G.Adj v w then 1 else 0)) := by
    intro w
    by_cases hw : w = v
    · subst hw
      have : c w ≠ y := fun h => hy h.symm
      simp [d_self, this]
    · have hu : Function.update c v y w = c w := Function.update_of_ne hw y c
      have hne : v ≠ w := Ne.symm hw
      simp only [d, pd, hne, ite_false, Function.update_self, hu, hw, ne_eq, not_false_eq_true,
        true_and]
      by_cases ha : G.Adj v w <;> by_cases h1 : c w = y <;> by_cases h2 : c w = c v <;>
        simp_all [eq_comm]
  rw [cost_update, Finset.sum_congr rfl (fun w _ => key w)]
  simp only [Finset.card_filter, Nat.cast_sum, Nat.cast_ite, Nat.cast_one, Nat.cast_zero,
    Finset.sum_sub_distrib, Finset.mul_sum]

section Weighted

variable {U : Type*} [Fintype U] [DecidableEq U]

/-- A weighted instance on the node set `U`. -/
structure WInst (U : Type*) where
  s : U → ℤ
  w : U → U → ℤ
  w_symm : ∀ i j, w i j = w j i

/-- Pair function of the weighted objective (per ordered pair of distinct nodes):
joined pairs pay the negative pairs `s i s j - w i j`, separated pairs pay the
positive pairs `w i j`. -/
def wpd (I : WInst U) (i j : U) (x y : α) : ℤ :=
  if i = j then 0 else if x = y then I.s i * I.s j - I.w i j else I.w i j

/-- Weighted cost (each unordered pair of nodes counted twice; the pairs inside a
node contribute a constant and are omitted). -/
def costW (I : WInst U) (c : U → α) : ℤ := ∑ i, ∑ j, wpd I i j (c i) (c j)

/-- Total size of the nodes with label `x`. -/
def S (I : WInst U) (c : U → α) (x : α) : ℤ := ∑ j, if c j = x then I.s j else 0

/-- Weight between node `i` and the other nodes with label `x`. -/
def W (I : WInst U) (c : U → α) (i : U) (x : α) : ℤ :=
  ∑ j, if j ≠ i ∧ c j = x then I.w i j else 0

omit [Fintype U] in
lemma wpd_symm (I : WInst U) (i j : U) (x y : α) : wpd I i j x y = wpd I j i y x := by
  unfold wpd
  by_cases h : i = j
  · subst h; simp
  · have h' : j ≠ i := Ne.symm h
    simp only [h, h', ite_false, @eq_comm _ x y, I.w_symm i j, mul_comm (I.s i)]

omit [Fintype U] in
lemma wpd_diag (I : WInst U) (i : U) (x y : α) : wpd I i i x y = 0 := by simp [wpd]

/-- Moving node `v` (size `s v`) to the cluster with label `y`. -/
theorem costW_move (I : WInst U) (c : U → α) (v : U) (y : α) (hy : y ≠ c v) :
    costW I (Function.update c v y) - costW I c =
      2 * (I.s v * (S I c y - S I c (c v) + I.s v) + 2 * (W I c v (c v) - W I c v y)) := by
  have h := sum_update_pairs (wpd I) (wpd_symm I) (wpd_diag I) c v y
  simp only [costW]
  rw [h]
  congr 1
  have key : ∀ j, wpd I v j y (Function.update c v y j) - wpd I v j (c v) (c j) =
      I.s v * ((if c j = y then I.s j else 0) - (if c j = c v then I.s j else 0) +
        (if j = v then I.s v else 0)) +
      2 * ((if j ≠ v ∧ c j = c v then I.w v j else 0) -
        (if j ≠ v ∧ c j = y then I.w v j else 0)) := by
    intro j
    by_cases hj : j = v
    · subst hj
      have : c j ≠ y := fun h => hy h.symm
      simp [wpd, this]
    · have hu : Function.update c v y j = c j := Function.update_of_ne hj y c
      have hne : v ≠ j := Ne.symm hj
      simp only [wpd, hne, ite_false, hu, hj, ne_eq, not_false_eq_true, true_and]
      by_cases h1 : y = c j <;> by_cases h2 : c v = c j <;>
        simp_all [eq_comm] <;> ring
  rw [Finset.sum_congr rfl (fun j _ => key j)]
  simp only [S, W, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum,
    Finset.sum_ite_eq', Finset.mem_univ, ite_true]

lemma S_update (I : WInst U) (c : U → α) (v : U) (y x : α) :
    S I (Function.update c v y) x =
      S I c x - (if c v = x then I.s v else 0) + (if y = x then I.s v else 0) := by
  have key : ∀ j, (if Function.update c v y j = x then I.s j else 0) =
      (if c j = x then I.s j else 0) - (if j = v then (if c v = x then I.s v else 0) else 0) +
        (if j = v then (if y = x then I.s v else 0) else 0) := by
    intro j
    by_cases hj : j = v
    · subst hj; simp
    · simp [hj]
  simp only [S]
  rw [Finset.sum_congr rfl (fun j _ => key j)]
  simp [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_ite_eq']

lemma W_update (I : WInst U) (c : U → α) (v i : U) (hi : i ≠ v) (y x : α) :
    W I (Function.update c v y) i x =
      W I c i x - (if c v = x then I.w i v else 0) + (if y = x then I.w i v else 0) := by
  have key : ∀ j, (if j ≠ i ∧ Function.update c v y j = x then I.w i j else 0) =
      (if j ≠ i ∧ c j = x then I.w i j else 0) -
        (if j = v then (if c v = x then I.w i v else 0) else 0) +
        (if j = v then (if y = x then I.w i v else 0) else 0) := by
    intro j
    by_cases hj : j = v
    · subst hj; simp [Ne.symm hi]
    · simp [hj]
  simp only [W]
  rw [Finset.sum_congr rfl (fun j _ => key j)]
  simp [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_ite_eq']

/-- Exchanging node `v` (label `a`) with node `u` (label `b`): the change is the
move of `v` to `b` plus the move of `u` to `a` evaluated after the first move,
with `S'_a = S_a - s_v`, `S'_b = S_b + s_v`, `w'(u, b - u) = w(u, b - u) + w(u, v)` and
`w'(u, a) = w(u, a) - w(u, v)`. -/
theorem costW_swap (I : WInst U) (c : U → α) (u v : U) (a b : α)
    (hv : c v = a) (hu : c u = b) (hab : a ≠ b) (huv : u ≠ v) :
    costW I (Function.update (Function.update c v b) u a) - costW I c =
      2 * (I.s v * (S I c b - S I c a + I.s v) + 2 * (W I c v a - W I c v b)) +
      2 * (I.s u * ((S I c a - I.s v) - (S I c b + I.s v) + I.s u) +
        2 * ((W I c u b + I.w u v) - (W I c u a - I.w u v))) := by
  set c1 := Function.update c v b with hc1
  have hc1u : c1 u = b := by rw [hc1, Function.update_of_ne huv, hu]
  have h1 := costW_move I c v b (by rw [hv]; exact Ne.symm hab)
  have h2 := costW_move I c1 u a (by rw [hc1u]; exact hab)
  rw [hc1u] at h2
  rw [hv] at h1
  have e1 : S I c1 a = S I c a - I.s v := by
    rw [hc1, S_update]; simp [hv, Ne.symm hab]
  have e2 : S I c1 b = S I c b + I.s v := by
    rw [hc1, S_update]; simp [hv, hab]
  have e3 : W I c1 u b = W I c u b + I.w u v := by
    rw [hc1, W_update I c v u huv]; simp [hv, hab]
  have e4 : W I c1 u a = W I c u a - I.w u v := by
    rw [hc1, W_update I c v u huv]; simp [hv, Ne.symm hab]
  rw [e1, e2, e3, e4] at h2
  linarith

end Weighted

end CC
