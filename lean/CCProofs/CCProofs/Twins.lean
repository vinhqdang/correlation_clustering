/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Basic

/-!
# Twins (critical cliques) are never split by an optimal clustering

Two vertices are twins when they are adjacent and have the same neighbours
(equal closed neighbourhoods).  If a clustering separates twins `u` and `v`,
moving `v` to the cluster of `u` and moving `u` to the cluster of `v` together
decrease `2 * (number of disagreements)` by exactly `4`, so one of the two moves
is a strict improvement (`twins_split_improvable`).  Hence every optimal
clustering keeps every pair of twins, and so every critical clique, together
(`optimal_keeps_twins`).  Pivot never separates twins either
(`pivot_round_twins`): in every round both or neither join the pivot's cluster.
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- `u` and `v` are twins: adjacent, and with the same other neighbours. -/
def Twins (u v : V) : Prop := G.Adj u v ∧ ∀ w, w ≠ u → w ≠ v → (G.Adj u w ↔ G.Adj v w)

theorem twins_split_improvable (c : V → α) (u v : V) (h : Twins G u v) (hc : c u ≠ c v) :
    (cost G (Function.update c v (c u)) - cost G c) +
      (cost G (Function.update c u (c v)) - cost G c) = -4 := by
  have huv : u ≠ v := h.1.ne
  rw [cost_update, cost_update]
  have key : ∀ w,
      (d G (Function.update c v (c u)) v w - d G c v w) +
        (d G (Function.update c u (c v)) u w - d G c u w) =
      (if w = u then -1 else 0) + (if w = v then -1 else 0) := by
    intro w
    by_cases hwu : w = u
    · subst hwu
      simp only [d, pd, Function.update_self, Function.update_of_ne huv, Ne.symm huv, huv,
        ite_false, ite_true, h.1.symm, Ne.symm hc]
      simp
    · by_cases hwv : w = v
      · subst hwv
        simp only [d, pd, Function.update_self, Function.update_of_ne (Ne.symm huv), huv,
          Ne.symm huv, ite_false, ite_true, h.1, hc]
        simp
      · have hadj : G.Adj u w ↔ G.Adj v w := h.2 w hwu hwv
        simp only [d, pd, Function.update_self, Function.update_of_ne hwu,
          Function.update_of_ne hwv, Ne.symm hwu, Ne.symm hwv, hwu, hwv, ite_false, hadj]
        by_cases ha : G.Adj v w <;> by_cases h1 : c u = c w <;> by_cases h2 : c v = c w <;>
          simp [ha, h1, h2]
  have hsum : ∑ w, ((d G (Function.update c v (c u)) v w - d G c v w) +
        (d G (Function.update c u (c v)) u w - d G c u w)) = -2 := by
    rw [Finset.sum_congr rfl (fun w _ => key w)]
    simp [Finset.sum_add_distrib, Finset.sum_ite_eq']
  rw [Finset.sum_add_distrib] at hsum
  linarith

/-- Every optimal clustering puts twins in the same cluster. -/
theorem optimal_keeps_twins (c : V → α) (hopt : ∀ c' : V → α, cost G c ≤ cost G c')
    (u v : V) (h : Twins G u v) : c u = c v := by
  by_contra hc
  have e := twins_split_improvable G c u v h hc
  have h1 := hopt (Function.update c v (c u))
  have h2 := hopt (Function.update c u (c v))
  linarith

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- A vertex `p` is `u` or a neighbour of `u` iff it is `v` or a neighbour of `v`. -/
theorem twins_closed_nbhd (u v : V) (h : Twins G u v) (p : V) :
    (p = u ∨ G.Adj p u) ↔ (p = v ∨ G.Adj p v) := by
  by_cases hpu : p = u
  · subst hpu; simp [h.1]
  · by_cases hpv : p = v
    · subst hpv; simp [h.1.symm]
    · have e := h.2 p hpu hpv
      simp only [hpu, hpv, false_or, G.adj_comm p u, G.adj_comm p v]
      exact e

omit [Fintype V] in
/-- One Pivot round on the remaining set `S` with pivot `p` forms the cluster of
`p` and its neighbours in `S`; twins in `S` are both in it or both outside. -/
theorem pivot_round_twins (S : Finset V) (p u v : V) (h : Twins G u v)
    (hu : u ∈ S) (hv : v ∈ S) :
    u ∈ S.filter (fun w => w = p ∨ G.Adj p w) ↔ v ∈ S.filter (fun w => w = p ∨ G.Adj p w) := by
  simp only [Finset.mem_filter, hu, hv, true_and]
  have e := twins_closed_nbhd G u v h p
  constructor
  · intro hx
    exact (e.mp (hx.imp Eq.symm id)).imp Eq.symm id
  · intro hx
    exact (e.mpr (hx.imp Eq.symm id)).imp Eq.symm id

end CC
