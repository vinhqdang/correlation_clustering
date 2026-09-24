/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Moves

/-!
# Far pairs are separated by every optimal clustering

If `u` and `w` are not adjacent and have no common neighbour (distance at least
three in `G`), no optimal clustering puts them in the same cluster
(`optimal_separates_far`).  Proof: if they share a cluster `C`, their
neighbourhoods inside `C` are disjoint subsets of `C \ {u, w}`, so the two moves
"`u` to a new cluster" and "`w` to a new cluster" have changes summing to at
most `-4` (ordered count); one of them is a strict improvement.  This justifies
restricting the LP to the pairs at distance at most two (far pairs fixed apart).
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- Moving `v` to an unused label `z` changes the cost by
`-2 (|A| - 1 - 2 k_A)`. -/
lemma cost_isolate (c : V → α) (v : V) (z : α) (hz : ∀ x, c x ≠ z) :
    cost G (Function.update c v z) - cost G c =
      -2 * (((univ.filter fun w => w ≠ v ∧ c w = c v).card : ℤ) -
        2 * (univ.filter fun w => c w = c v ∧ G.Adj v w).card) := by
  rw [cost_move G c v z (fun h => hz v h.symm)]
  have e1 : (univ.filter fun w => c w = z) = ∅ := by
    ext x; simp [hz x]
  have e2 : (univ.filter fun w => c w = z ∧ G.Adj v w) = ∅ := by
    ext x; simp [hz x]
  rw [e1, e2]
  simp
  ring

theorem optimal_separates_far (c : V → α) (hopt : ∀ c' : V → α, cost G c ≤ cost G c')
    (z : α) (hz : ∀ x, c x ≠ z) (u w : V) (huw : u ≠ w) (hnadj : ¬ G.Adj u w)
    (hcommon : ∀ x, ¬ (G.Adj u x ∧ G.Adj w x)) : c u ≠ c w := by
  intro hc
  set C := univ.filter fun x => c x = c u with hC
  have hu : u ∈ C := by simp [hC]
  have hw : w ∈ C := by simp [hC, hc]
  -- sizes of the cluster minus one
  have nu : (univ.filter fun x => x ≠ u ∧ c x = c u) = C.erase u := by
    ext x; simp [hC, and_comm]
  have nw : (univ.filter fun x => x ≠ w ∧ c x = c w) = C.erase w := by
    ext x; simp [hC, hc, and_comm]
  -- neighbourhoods inside the cluster
  set Ku := univ.filter fun x => c x = c u ∧ G.Adj u x with hKu
  set Kw := univ.filter fun x => c x = c w ∧ G.Adj w x with hKw
  have hdisj : Disjoint Ku Kw := by
    rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [hKu, hKw, Finset.mem_filter, Finset.mem_univ, true_and] at hx hx'
    exact hcommon x ⟨hx.2, hx'.2⟩
  have hsub : Ku ∪ Kw ⊆ (C.erase u).erase w := by
    intro x hx
    simp only [Finset.mem_union, hKu, hKw, Finset.mem_filter, Finset.mem_univ,
      true_and] at hx
    simp only [Finset.mem_erase, hC, Finset.mem_filter, Finset.mem_univ, true_and]
    rcases hx with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · refine ⟨?_, h2.ne.symm, h1⟩
      rintro rfl; exact hnadj h2
    · refine ⟨h2.ne.symm, ?_, by rw [h1, hc]⟩
      rintro rfl; exact hnadj h2.symm
  have hcard : (Ku.card : ℤ) + Kw.card ≤ (C.card : ℤ) - 2 := by
    have h1 := Finset.card_le_card hsub
    rw [Finset.card_union_of_disjoint hdisj] at h1
    have h2 : ((C.erase u).erase w).card = C.card - 2 := by
      rw [Finset.card_erase_of_mem (Finset.mem_erase.mpr ⟨Ne.symm huw, hw⟩),
        Finset.card_erase_of_mem hu]
      omega
    have h3 : 2 ≤ C.card := by
      have : ({u, w} : Finset V) ⊆ C := by
        intro x hx; simp at hx; rcases hx with rfl | rfl; exact hu; exact hw
      have := Finset.card_le_card this
      rw [Finset.card_pair huw] at this
      exact this
    omega
  have du := cost_isolate G c u z hz
  have dw := cost_isolate G c w z hz
  rw [nu] at du
  rw [nw] at dw
  rw [← hc] at dw
  have ou := hopt (Function.update c u z)
  have ow := hopt (Function.update c w z)
  have cu : ((C.erase u).card : ℤ) = C.card - 1 := by
    rw [Finset.card_erase_of_mem hu]; have : 1 ≤ C.card := Finset.card_pos.mpr ⟨u, hu⟩
    omega
  have cw : ((C.erase w).card : ℤ) = C.card - 1 := by
    rw [Finset.card_erase_of_mem hw]; have : 1 ≤ C.card := Finset.card_pos.mpr ⟨w, hw⟩
    omega
  rw [cu] at du
  rw [cw] at dw
  have hKw' : (univ.filter fun x => c x = c u ∧ G.Adj w x) = Kw := by
    rw [hKw, hc]
  rw [hKw'] at dw
  linarith

end CC
