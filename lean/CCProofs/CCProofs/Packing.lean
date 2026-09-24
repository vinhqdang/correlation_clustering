/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Basic

/-!
# Packing lower bounds

Let `H k` (`k : ι`) be vertex sets with weights `y k ≥ 0` such that every pair of
distinct vertices is covered with total weight at most `1`, and let `m k` be a
lower bound on the cost that a clustering pays inside `H k`.  Then every
clustering costs at least `∑ k, y k * m k` (`packing_bound`).  This is the lower
bound certificate behind the triangle, star and local-subgraph packings.

* `bad_triangle`: a clustering makes at least one mistake on every bad triangle
  (two positive pairs and one negative pair), i.e. `m = 1` (`2` in the ordered
  count used by `cost`).
* `star_bound`: on a star with centre `v` and `k ≥ 1` pairwise non-adjacent
  leaves a clustering makes at least `k - 1` mistakes.
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- Cost of `c` on the (ordered) pairs inside `H`. -/
def costIn (c : V → α) (H : Finset V) : ℤ := ∑ u ∈ H, ∑ v ∈ H, d G c u v

omit [Fintype V] in
lemma d_nonneg (c : V → α) (u v : V) : 0 ≤ d G c u v := by
  unfold d pd; split_ifs <;> norm_num

omit [Fintype V] in
lemma d_le_one (c : V → α) (u v : V) : d G c u v ≤ 1 := by
  unfold d pd; split_ifs <;> norm_num

/-- Weak duality for packings of vertex sets. -/
theorem packing_bound {ι : Type*} [Fintype ι] (H : ι → Finset V) (y m : ι → ℝ)
    (hy : ∀ k, 0 ≤ y k)
    (hload : ∀ u v, u ≠ v → ∑ k, (if u ∈ H k ∧ v ∈ H k then y k else 0) ≤ 1)
    (c : V → α) (hm : ∀ k, m k ≤ (costIn G c (H k) : ℝ)) :
    ∑ k, y k * m k ≤ (cost G c : ℝ) := by
  -- the load-weighted cost is at most the cost, pair by pair
  have hpair : ∀ u v, (d G c u v : ℝ) *
      (∑ k, (if u ∈ H k ∧ v ∈ H k then y k else 0)) ≤ (d G c u v : ℝ) := by
    intro u v
    by_cases huv : u = v
    · subst huv; simp [d_self]
    · have h0 : (0 : ℝ) ≤ d G c u v := by exact_mod_cast d_nonneg G c u v
      have := hload u v huv
      nlinarith
  have hin : ∀ k, (costIn G c (H k) : ℝ) =
      ∑ u, ∑ v, (if u ∈ H k ∧ v ∈ H k then (d G c u v : ℝ) else 0) := by
    intro k
    have hrow : ∀ u, (∑ v, (if u ∈ H k ∧ v ∈ H k then (d G c u v : ℝ) else 0)) =
        if u ∈ H k then ∑ v ∈ H k, (d G c u v : ℝ) else 0 := by
      intro u
      by_cases hu : u ∈ H k
      · simp only [hu, true_and, ite_true]
        rw [← Finset.sum_filter, Finset.filter_mem_eq_inter, Finset.univ_inter]
      · simp [hu]
    simp only [hrow]
    rw [← Finset.sum_filter, Finset.filter_mem_eq_inter, Finset.univ_inter]
    simp [costIn]
  calc ∑ k, y k * m k
      ≤ ∑ k, y k * (costIn G c (H k) : ℝ) :=
        Finset.sum_le_sum (fun k _ => mul_le_mul_of_nonneg_left (hm k) (hy k))
    _ = ∑ u, ∑ v, (d G c u v : ℝ) * (∑ k, (if u ∈ H k ∧ v ∈ H k then y k else 0)) := by
        simp only [hin, Finset.mul_sum]
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl (fun u _ => ?_)
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl (fun v _ => ?_)
        refine Finset.sum_congr rfl (fun k _ => ?_)
        split_ifs <;> ring
    _ ≤ ∑ u, ∑ v, (d G c u v : ℝ) :=
        Finset.sum_le_sum (fun u _ => Finset.sum_le_sum (fun v _ => hpair u v))
    _ = (cost G c : ℝ) := by simp [cost]

omit [Fintype V] in
/-- Every clustering makes a mistake on a bad triangle. -/
theorem bad_triangle (c : V → α) (u v w : V) (huv : u ≠ v) (hvw : v ≠ w) (huw : u ≠ w)
    (h1 : G.Adj u v) (h2 : G.Adj v w) (h3 : ¬ G.Adj u w) :
    1 ≤ d G c u v + d G c v w + d G c u w := by
  unfold d pd
  simp only [huv, hvw, huw, ite_false, h1, h2, h3, true_iff, false_iff]
  by_cases a : c u = c v <;> by_cases b : c v = c w <;> by_cases e : c u = c w <;> simp_all

omit [Fintype V] in
/-- The same in the ordered count used by `cost`: at least `2` inside the triangle. -/
theorem bad_triangle_costIn (c : V → α) (u v w : V) (huv : u ≠ v) (hvw : v ≠ w)
    (huw : u ≠ w) (h1 : G.Adj u v) (h2 : G.Adj v w) (h3 : ¬ G.Adj u w) :
    2 ≤ costIn G c {u, v, w} := by
  have hb := bad_triangle G c u v w huv hvw huw h1 h2 h3
  have hu : u ∉ ({v, w} : Finset V) := by simp [huv, huw]
  have hv : v ∉ ({w} : Finset V) := by simp [hvw]
  simp only [costIn, Finset.sum_insert hu, Finset.sum_insert hv, Finset.sum_singleton,
    d_self]
  have e1 := d_comm G c u v
  have e2 := d_comm G c v w
  have e3 := d_comm G c u w
  linarith

omit [Fintype V] in
/-- Star bound: centre `v`, leaves `L` (`k = |L|`), every leaf adjacent to `v`
and no two leaves adjacent.  Every clustering makes at least `k - 1` mistakes on
the pairs of the star (`2 (k - 1)` in the ordered count); for `k = 0` the bound is trivial. -/
theorem star_bound (c : V → α) (v : V) (L : Finset V) (hv : v ∉ L)
    (hadj : ∀ l ∈ L, G.Adj v l) (hind : ∀ l ∈ L, ∀ l' ∈ L, ¬ G.Adj l l') :
    2 * ((L.card : ℤ) - 1) ≤ costIn G c (insert v L) := by
  classical
  -- J: the leaves in the cluster of the centre
  set J := L.filter (fun l => c l = c v) with hJ
  -- mistakes between the centre and the leaves outside J
  have hcentre : ∑ l ∈ L, d G c v l = ((L.card : ℤ) - J.card) := by
    have : ∀ l ∈ L, d G c v l = if c l = c v then 0 else 1 := by
      intro l hl
      have hne : v ≠ l := fun e => hv (e ▸ hl)
      unfold d pd
      simp only [hne, ite_false, hadj l hl, true_iff]
      by_cases h : c l = c v
      · simp [h]
      · simp [h, Ne.symm h]
    rw [Finset.sum_congr rfl this, Finset.sum_ite, Finset.sum_const_zero, zero_add,
      Finset.sum_const, nsmul_eq_mul, mul_one]
    have hc := Finset.card_filter_add_card_filter_not (s := L) (fun l => c l = c v)
    rw [← hJ] at hc
    omega
  -- mistakes among the leaves of J (negative pairs joined)
  have hleaves : ∑ l ∈ L, ∑ l' ∈ L, d G c l l' ≥ (J.card : ℤ) * J.card - J.card := by
    have hsub : ∀ l ∈ L, ∑ l' ∈ L, d G c l l' ≥
        ∑ l' ∈ J, (if l ∈ J ∧ l' ≠ l then 1 else 0) := by
      intro l hl
      calc ∑ l' ∈ L, d G c l l' ≥ ∑ l' ∈ J, d G c l l' := by
            apply Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
            intro i _ _; exact d_nonneg G c l i
        _ ≥ ∑ l' ∈ J, (if l ∈ J ∧ l' ≠ l then 1 else 0) := by
            apply Finset.sum_le_sum
            intro l' hl'
            by_cases h : l ∈ J ∧ l' ≠ l
            · have hlJ := (Finset.mem_filter.mp h.1)
              have hl'J := (Finset.mem_filter.mp hl')
              have hne : l ≠ l' := fun e => h.2 e.symm
              have heq : c l = c l' := by rw [hlJ.2, hl'J.2]
              simp only [h, d, pd, hne, ite_false,
                hind l hlJ.1 l' hl'J.1, heq, false_iff, not_true_eq_false]
              simp [Ne.symm hne]
            · simp only [h, ite_false]; exact d_nonneg G c l l'
    calc ∑ l ∈ L, ∑ l' ∈ L, d G c l l'
        ≥ ∑ l ∈ L, ∑ l' ∈ J, (if l ∈ J ∧ l' ≠ l then (1 : ℤ) else 0) :=
          Finset.sum_le_sum (fun l hl => hsub l hl)
      _ ≥ ∑ l ∈ J, ∑ l' ∈ J, (if l ∈ J ∧ l' ≠ l then (1 : ℤ) else 0) := by
          apply Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
          intro i _ _; apply Finset.sum_nonneg; intro j _; split_ifs <;> norm_num
      _ = ∑ l ∈ J, ((J.card : ℤ) - 1) := by
          refine Finset.sum_congr rfl (fun l hl => ?_)
          have : ∑ l' ∈ J, (if l ∈ J ∧ l' ≠ l then (1 : ℤ) else 0) =
              ∑ l' ∈ J.erase l, (1 : ℤ) := by
            rw [← Finset.sum_filter]
            congr 1
            ext x; simp [hl, Finset.mem_erase, and_comm]
          rw [this, Finset.sum_const, Finset.card_erase_of_mem hl]
          simp [Nat.cast_sub (Finset.card_pos.mpr ⟨l, hl⟩)]
      _ = (J.card : ℤ) * J.card - J.card := by
          rw [Finset.sum_const, nsmul_eq_mul]; ring
  -- total: 2 (k - j) + j (j - 1) ≥ 2 (k - 1) for every integer j ≥ 0
  have hsplit : costIn G c (insert v L) =
      d G c v v + ∑ l ∈ L, d G c v l + ∑ l ∈ L, d G c l v + ∑ l ∈ L, ∑ l' ∈ L, d G c l l' := by
    simp only [costIn, Finset.sum_insert hv, Finset.sum_add_distrib]
    ring
  have hsym : ∑ l ∈ L, d G c l v = ∑ l ∈ L, d G c v l :=
    Finset.sum_congr rfl (fun l _ => d_comm G c l v)
  rw [hsplit, hsym, d_self, hcentre]
  have hj : (0 : ℤ) ≤ J.card := by positivity
  nlinarith [sq_nonneg ((J.card : ℤ) - 1), sq_nonneg ((J.card : ℤ) - 2)]

end CC
