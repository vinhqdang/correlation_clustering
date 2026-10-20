import P3Proofs.Struct.Defs2

/-!
# Twin lemma (paper: `lem:twin`, app_structure.tex)

Vertices `u ≠ v` with `N[u] = N[v]` (they are then adjacent).  For every labelling
`c` with `c u ≠ c v` the two cost changes (move `v` into the cluster of `u`, move `u` into
the cluster of `v`) sum to exactly `-4` in ordered units, i.e. `-2` in the paper's units;
so one of them is a strict improvement, and `u`, `v` share a cluster in every optimum.
No hypothesis on a packing is used.
-/

open Finset

namespace P3S

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- `u ≠ v` with equal closed neighbourhoods `N[u] = N[v]`. -/
def SameClosed (u v : V) : Prop :=
  u ≠ v ∧ ∀ w, (w = u ∨ G.Adj u w) ↔ (w = v ∨ G.Adj v w)

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
lemma SameClosed.adj {u v : V} (h : SameClosed G u v) : G.Adj u v := by
  have := (h.2 v).mpr (Or.inl rfl)
  rcases this with e | e
  · exact absurd e.symm h.1
  · exact e

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
lemma SameClosed.other {u v : V} (h : SameClosed G u v) (w : V) (hu : w ≠ u) (hv : w ≠ v) :
    G.Adj u w ↔ G.Adj v w := by
  have := h.2 w
  simp only [hu, hv, false_or] at this
  exact this

/-- Strict form: the two moves change the cost by `-4` in total (ordered units). -/
theorem twin_moves_sum (c : V → α) (u v : V) (h : SameClosed G u v) (hc : c u ≠ c v) :
    (cost G (Function.update c v (c u)) - cost G c) +
      (cost G (Function.update c u (c v)) - cost G c) = -4 := by
  have huv : u ≠ v := h.1
  have hadj : G.Adj u v := h.adj G
  rw [cost_update, cost_update]
  have key : ∀ w,
      (d G (Function.update c v (c u)) v w - d G c v w) +
        (d G (Function.update c u (c v)) u w - d G c u w) =
      (if w = u then -1 else 0) + (if w = v then -1 else 0) := by
    intro w
    by_cases hwu : w = u
    · subst hwu
      simp only [d, pd, Function.update_self, Function.update_of_ne huv, Ne.symm huv, huv,
        ite_false, ite_true, hadj.symm, Ne.symm hc]
      simp
    · by_cases hwv : w = v
      · subst hwv
        simp only [d, pd, Function.update_self, Function.update_of_ne (Ne.symm huv), huv,
          Ne.symm huv, ite_false, ite_true, hadj, hc]
        simp
      · have hadj' : G.Adj u w ↔ G.Adj v w := h.other G w hwu hwv
        simp only [d, pd, Function.update_self, Function.update_of_ne hwu,
          Function.update_of_ne hwv, Ne.symm hwu, Ne.symm hwv, hwu, hwv, ite_false, hadj']
        by_cases ha : G.Adj v w <;> by_cases h1 : c u = c w <;> by_cases h2 : c v = c w <;>
          simp [ha, h1, h2]
  have hsum : ∑ w, ((d G (Function.update c v (c u)) v w - d G c v w) +
        (d G (Function.update c u (c v)) u w - d G c u w)) = -2 := by
    rw [Finset.sum_congr rfl (fun w _ => key w)]
    simp [Finset.sum_add_distrib, Finset.sum_ite_eq']
  rw [Finset.sum_add_distrib] at hsum
  linarith

/-- One of the two moves is a strict improvement (by at least `2` in ordered units). -/
theorem twin_split_improvable (c : V → α) (u v : V) (h : SameClosed G u v) (hc : c u ≠ c v) :
    cost G (Function.update c v (c u)) < cost G c ∨
      cost G (Function.update c u (c v)) < cost G c := by
  have e := twin_moves_sum G c u v h hc
  by_contra hcon
  push_neg at hcon
  linarith [hcon.1, hcon.2]

/-- **Twin lemma.** Every optimal clustering puts `u`, `v` with `N[u] = N[v]` in one cluster. -/
theorem optimal_keeps_twins (c : V → α) (hopt : IsOpt G c) (u v : V) (h : SameClosed G u v) :
    c u = c v := by
  by_contra hc
  have e := twin_moves_sum G c u v h hc
  have h1 := hopt (Function.update c v (c u))
  have h2 := hopt (Function.update c u (c v))
  linarith

end P3S
