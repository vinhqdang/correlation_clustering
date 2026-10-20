import Mathlib

/-!
# Self-contained definitions for the structural lemmas

`G` is the graph on a finite vertex set `V`; a clustering is a labelling `c : V → α`.
`d G c u v` is `1` if the ordered pair `(u, v)`, `u ≠ v`, is a disagreement (an edge that
is split or a non-edge that is joined).  `cost G c` sums `d` over ordered pairs, hence
`cost = 2 * (number of disagreements)`; "OPT ≤ k" is therefore `cost ≤ 2 * k`
(`Sol`).  The single-vertex move formulas `cost_update`, `cost_move`, `cost_isolate`
are the tools for all optimality arguments.

Packings: a triple `(a, b, c)` is an induced `P₃` with centre `b` (`IsP3`).  `pr P` is the
set of the six ordered pairs inside `P`.  A packing is a finite set of induced `P₃`s
with pairwise disjoint `pr` (modification-disjointness).  A pair is soft if it is in
`pr P` for a member, hard otherwise.
-/

open Finset

namespace P3S

section Generic

variable {U : Type*} [Fintype U] [DecidableEq U] {α : Type*}

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

/-- Disagreement of the pair `(u, w)` when `u` has label `x` and `w` has label `y`. -/
def pd (u w : V) (x y : α) : ℤ :=
  if u = w then 0 else if (G.Adj u w ↔ x = y) then 0 else 1

/-- Disagreement indicator of the ordered pair `(u, v)` under `c`. -/
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

omit [Fintype V] in
lemma d_nonneg (c : V → α) (u v : V) : 0 ≤ d G c u v := by
  unfold d pd; split_ifs <;> simp

/-- Relabelling `v` changes the cost only through the pairs `(v, w)`. -/
theorem cost_update (c : V → α) (v : V) (y : α) :
    cost G (Function.update c v y) - cost G c =
      2 * ∑ w, (d G (Function.update c v y) v w - d G c v w) := by
  have h := sum_update_pairs (pd G) (pd_symm G) (pd_diag G) c v y
  simp only [cost, d]
  rw [h, Function.update_self]

/-- Change of `cost` when `v` moves to the cluster with label `y ≠ c v`. -/
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

/-- Moving `v` to an unused label `z` changes the cost by `-2 (|A| - 1 - 2 k_A)`. -/
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

/-- `c` is an optimal clustering among all labellings by `α`. -/
def IsOpt (c : V → α) : Prop := ∀ c' : V → α, cost G c ≤ cost G c'

/-- `OPT ≤ k` in the paper's units (number of edits): some labelling has at most `k`
disagreeing pairs, i.e. `cost ≤ 2 k`. -/
def Sol (α : Type*) [DecidableEq α] (k : ℤ) : Prop := ∃ c : V → α, cost G c ≤ 2 * k

/-- A label of `α` not used by `c` (automatic when `α` is infinite). -/
lemma exists_fresh [Infinite α] (c : V → α) : ∃ z, ∀ x, c x ≠ z := by
  by_contra h
  have hs : Function.Surjective c := by
    intro z
    by_contra hz
    exact h ⟨z, fun x hx => hz ⟨x, hx⟩⟩
  have := Finite.of_surjective c hs
  exact not_finite α

end P3S
