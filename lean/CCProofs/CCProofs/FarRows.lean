/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.FarDual

/-!
# Rows valid for far-separating clusterings

The row families of `ccbench/blockdual.py`, written in the separation variables
`x = sep c` and checked for every clustering `c` that separates the far pairs
(`SepFar`); with `far_dual_bound` they give certified lower bounds.

* (t1) `sep_triangle` (in `Lagrangian.lean`): `x_uw ≤ x_uv + x_vw`, valid for every
  clustering.
* (t2) `sep_far_triangle`: if `(u, w)` is far, `1 ≤ x_uv + x_vw`.
* (star) `star_row_valid`: for a centre `v`, a finset `T` of leaves and the set `R`
  of the leaf pairs that are not far,
  `∑_{t ∈ T} x_vt - ∑_{{t, t'} ∈ R} x_tt' ≥ |T| - |R| - 1`.
  Unordered leaf pairs are written `(t, t')` with `t < t'` for a linear order on
  the vertices (`starPairs`), as in the certificate files.

Proof of the star row: let `S` be the leaves in the cluster of `v`.  Then
`∑_T x_vt = |T| - |S|`; the pairs of `S` are joined, hence not far, hence in `R`,
and `x = 0` on them.  With `m` the least element of `S`, the `|S| - 1` pairs
`(m, t)`, `t ∈ S`, `t ≠ m`, lie in `R` and cost nothing, so
`∑_R x ≤ |R| - |S| + 1` (also when `S` is empty).
-/

open Finset

namespace CC

section Triangle

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- (t2) A path `u - v - w` between a far pair is cut at least once. -/
theorem sep_far_triangle (c : V → α) (hc : SepFar G c) (u v w : V) (h : far G u w) :
    1 ≤ sep c (u, v) + sep c (v, w) := by
  have huw := hc u w h
  unfold sep
  by_cases a : c u = c v
  · have b : c v ≠ c w := fun e => huw (a.trans e)
    simp [a, b]
  · have := (sep_mem c (v, w)).1
    simp only [a, ite_false]
    unfold sep at this
    linarith

end Triangle

section Star

variable {V : Type*} [Fintype V] [LinearOrder V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- The leaf pairs `(t, t')`, `t < t'`, of `T` that are not far. -/
def starPairs (T : Finset V) : Finset (V × V) :=
  (T ×ˢ T).filter fun p => p.1 < p.2 ∧ ¬ far G p.1 p.2

omit [Fintype V] [LinearOrder V] in
lemma sum_sep_centre (c : V → α) (v : V) (T : Finset V) :
    ∑ t ∈ T, sep c (v, t) = (T.card : ℝ) - (T.filter fun t => c t = c v).card := by
  have h : ∀ t ∈ T, sep c (v, t) = 1 - (if c t = c v then (1 : ℝ) else 0) := by
    intro t _
    by_cases h : c t = c v
    · simp [sep, h]
    · simp [sep, h, Ne.symm h]
  rw [Finset.sum_congr rfl h, Finset.sum_sub_distrib, Finset.natCast_card_filter]
  simp

/-- The joined leaf pairs of the star cost nothing, and there are at least `|S| - 1`
of them in `R`. -/
lemma sum_sep_leaves (c : V → α) (hc : SepFar G c) (v : V) (T : Finset V) :
    ∑ p ∈ starPairs G T, sep c p ≤
      ((starPairs G T).card : ℝ) - (T.filter fun t => c t = c v).card + 1 := by
  set S := T.filter fun t => c t = c v with hS
  set R := starPairs G T with hR
  have hle : ∀ Q : Finset (V × V), Q ⊆ R → ∑ p ∈ R \ Q, sep c p ≤ ((R \ Q).card : ℝ) := by
    intro Q _
    calc ∑ p ∈ R \ Q, sep c p ≤ ∑ p ∈ R \ Q, (1 : ℝ) :=
          Finset.sum_le_sum (fun p _ => (sep_mem c p).2)
      _ = ((R \ Q).card : ℝ) := by simp
  rcases S.eq_empty_or_nonempty with hS0 | hS1
  · rw [hS0]
    have := hle ∅ (Finset.empty_subset _)
    simp only [Finset.sdiff_empty] at this
    simp only [Finset.card_empty, Nat.cast_zero]
    linarith
  · set m := S.min' hS1 with hm
    have hmS : m ∈ S := Finset.min'_mem S hS1
    set Q := (S.erase m).image (fun t => (m, t)) with hQ
    have hQR : Q ⊆ R := by
      intro p hp
      simp only [hQ, Finset.mem_image, Finset.mem_erase] at hp
      obtain ⟨t, ⟨htm, htS⟩, rfl⟩ := hp
      have hmT := (Finset.mem_filter.mp hmS)
      have htT := (Finset.mem_filter.mp htS)
      have hlt : m < t := lt_of_le_of_ne (Finset.min'_le S t htS) (Ne.symm htm)
      have hnf : ¬ far G m t := fun hf => hc m t hf (by rw [hmT.2, htT.2])
      simp [hR, starPairs, hmT.1, htT.1, hlt, hnf]
    have hQ0 : ∑ p ∈ Q, sep c p = 0 := by
      refine Finset.sum_eq_zero (fun p hp => ?_)
      simp only [hQ, Finset.mem_image, Finset.mem_erase] at hp
      obtain ⟨t, ⟨_, htS⟩, rfl⟩ := hp
      have hmT := (Finset.mem_filter.mp hmS)
      have htT := (Finset.mem_filter.mp htS)
      simp [sep, hmT.2, htT.2]
    have hQc : Q.card + 1 = S.card := by
      rw [hQ, Finset.card_image_of_injective _ (fun a b h => (Prod.mk.inj h).2),
        Finset.card_erase_add_one hmS]
    have hsplit := Finset.sum_sdiff (f := fun p => sep c p) hQR
    have hcard := Finset.card_sdiff_add_card_eq_card hQR
    have h1 := hle Q hQR
    have hc1 : ((R \ Q).card : ℝ) + Q.card = R.card := by exact_mod_cast hcard
    have hc2 : (Q.card : ℝ) + 1 = S.card := by exact_mod_cast hQc
    linarith

/-- (star) The star row is valid for every far-separating clustering. -/
theorem star_row_valid (c : V → α) (hc : SepFar G c) (v : V) (T : Finset V) :
    (T.card : ℝ) - (starPairs G T).card - 1 ≤
      ∑ t ∈ T, sep c (v, t) - ∑ p ∈ starPairs G T, sep c p := by
  rw [sum_sep_centre]
  have := sum_sep_leaves G c hc v T
  linarith

end Star

end CC
