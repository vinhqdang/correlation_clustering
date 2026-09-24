/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Moves

/-!
# The objective depends only on the partition

The cost of a clustering depends only on which pairs share a label, not on the
labels themselves.

* `cost_congr`, `costW_congr`: two labellings (possibly with different label
  types) that put the same pairs together have the same cost, unweighted and
  weighted.
* `cost_comp_injective`, `costW_comp_injective`: relabelling by an injective map
  (for example the embedding of the labels `α ⊕ β` of a crossover child into the
  label type of the population) does not change the cost.
* `canon`: every labelling has a labelling by vertices with the same clusters
  (each vertex is labelled by a representative of its cluster), so minimising over
  `V → V` is minimising over all labellings (`cost_canon`).
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α β : Type*} [DecidableEq α] [DecidableEq β]

omit [Fintype V] in
lemma d_congr (c : V → α) (c' : V → β) (h : ∀ u v, c u = c v ↔ c' u = c' v) (u v : V) :
    d G c u v = d G c' u v := by
  simp only [d, pd, h u v]

/-- Labellings with the same clusters have the same cost. -/
theorem cost_congr (c : V → α) (c' : V → β) (h : ∀ u v, c u = c v ↔ c' u = c' v) :
    cost G c = cost G c' := by
  simp only [cost, d_congr G c c' h]

/-- Relabelling invariance: an injective map of the labels does not change the cost. -/
theorem cost_comp_injective (e : α → β) (he : Function.Injective e) (c : V → α) :
    cost G (e ∘ c) = cost G c :=
  cost_congr G _ _ (fun _ _ => he.eq_iff)

/-- A labelling by vertices with the same clusters as `c`: every vertex is
labelled by a representative of its cluster. -/
noncomputable def canon {δ : Type*} (c : V → δ) : V → V :=
  fun u => @Classical.epsilon V ⟨u⟩ (fun v => c v = c u)

omit [Fintype V] [DecidableEq V] in
lemma canon_spec {δ : Type*} (c : V → δ) (u : V) : c (canon c u) = c u :=
  Classical.epsilon_spec (p := fun v => c v = c u) ⟨u, rfl⟩

omit [Fintype V] [DecidableEq V] in
lemma canon_iff {δ : Type*} (c : V → δ) (u v : V) : canon c u = canon c v ↔ c u = c v := by
  constructor
  · intro h
    rw [← canon_spec c u, h, canon_spec c v]
  · intro h
    unfold canon
    rw [h]

theorem cost_canon (c : V → α) : cost G (canon c) = cost G c :=
  cost_congr G _ _ (fun u v => canon_iff c u v)

section Weighted

variable {U : Type*} [Fintype U] [DecidableEq U]

/-- Weighted labellings with the same clusters have the same weighted cost. -/
theorem costW_congr (I : WInst U) (c : U → α) (c' : U → β)
    (h : ∀ i j, c i = c j ↔ c' i = c' j) : costW I c = costW I c' := by
  simp only [costW, wpd, h]

/-- Relabelling invariance of the weighted cost. -/
theorem costW_comp_injective (I : WInst U) (e : α → β) (he : Function.Injective e)
    (c : U → α) : costW I (e ∘ c) = costW I c :=
  costW_congr I _ _ (fun _ _ => he.eq_iff)

end Weighted

end CC
