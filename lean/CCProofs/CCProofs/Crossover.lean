/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Basic

/-!
# Partition crossover

Let `a` and `b` be two clusterings and `blk` a block labelling that is constant on
every cluster of `a` and on every cluster of `b` (every block is a union of
clusters of both parents; the connected components of the graph that links each
cluster of `a` to the clusters of `b` it meets are the finest such blocks).  The
child takes, independently in every block `X`, the clustering of `b` if `T X` and
the clustering of `a` otherwise.

* `d_child`: on every pair the child has the disagreement of the chosen parent;
  pairs in different blocks are separated by `a`, `b` and the child alike.
* `cost_child_a`, `cost_child_b`: the cost of the child is the cost of a parent
  plus the sum, over the blocks taken from the other parent, of the difference of
  the two parents' block costs.
* `crossover_le`: choosing `T X := bcost b X < bcost a X`, the child costs no more
  than either parent.
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α β γ : Type*} [DecidableEq α] [DecidableEq β] [DecidableEq γ]

/-- Every block of `blk` is a union of clusters of `a` and a union of clusters of `b`. -/
structure Compatible (a : V → α) (b : V → β) (blk : V → γ) : Prop where
  ha : ∀ u v, a u = a v → blk u = blk v
  hb : ∀ u v, b u = b v → blk u = blk v

/-- The child: clusters of `b` on the blocks `X` with `T X`, clusters of `a` elsewhere. -/
def child (a : V → α) (b : V → β) (blk : V → γ) (T : γ → Prop) [DecidablePred T] :
    V → α ⊕ β :=
  fun u => if T (blk u) then Sum.inr (b u) else Sum.inl (a u)

/-- Cost of the clustering `p` on the pairs inside block `X`. -/
def bcost {δ : Type*} [DecidableEq δ] (p : V → δ) (blk : V → γ) (X : γ) : ℤ :=
  ∑ u, ∑ v, if blk u = X ∧ blk v = X then d G p u v else 0

variable {a : V → α} {b : V → β} {blk : V → γ}

omit [Fintype V] [DecidableEq γ] in
/-- Pairs in different blocks are separated by both parents, so the parents agree there. -/
lemma d_cross (h : Compatible a b blk) (u v : V) (huv : blk u ≠ blk v) :
    d G a u v = d G b u v := by
  have ha : a u ≠ a v := fun e => huv (h.ha u v e)
  have hb : b u ≠ b v := fun e => huv (h.hb u v e)
  simp [d, pd, ha, hb]

omit [Fintype V] in
lemma d_child (h : Compatible a b blk) (T : γ → Prop) [DecidablePred T] (u v : V) :
    d G (child a b blk T) u v = if T (blk u) then d G b u v else d G a u v := by
  by_cases hbl : blk u = blk v
  · by_cases hT : T (blk u)
    · have hT' : T (blk v) := hbl ▸ hT
      simp [d, pd, child, hT, hT']
    · have hT' : ¬ T (blk v) := hbl ▸ hT
      simp [d, pd, child, hT, hT']
  · have ha : a u ≠ a v := fun e => hbl (h.ha u v e)
    have hb : b u ≠ b v := fun e => hbl (h.hb u v e)
    have hc : child a b blk T u ≠ child a b blk T v := by
      unfold child
      by_cases h1 : T (blk u) <;> by_cases h2 : T (blk v) <;> simp [h1, h2, ha, hb]
    have e := d_cross G h u v hbl
    by_cases hT : T (blk u)
    · simp only [hT, ite_true]
      simp [d, pd, hc, hb]
    · simp only [hT, ite_false]
      simp [d, pd, hc, ha]

omit [DecidableEq V] in
/-- Regrouping a double sum of pair terms supported inside blocks by blocks. -/
lemma sum_blocks (e : V → V → ℤ) (T : γ → Prop) [DecidablePred T] :
    ∑ X ∈ univ.image blk, (if T X then ∑ u, ∑ v, (if blk u = X ∧ blk v = X then e u v else 0)
      else 0) =
      ∑ u, ∑ v, if blk u = blk v ∧ T (blk u) then e u v else 0 := by
  have step : ∀ X, (if T X then ∑ u, ∑ v, (if blk u = X ∧ blk v = X then e u v else 0) else 0) =
      ∑ u, ∑ v, (if blk u = X then (if blk u = blk v ∧ T (blk u) then e u v else 0) else 0) := by
    intro X
    by_cases hT : T X
    · simp only [hT, ite_true]
      refine Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))
      by_cases h1 : blk u = X
      · subst h1; by_cases h2 : blk v = blk u
        · simp [h2, hT]
        · simp [h2, Ne.symm h2]
      · simp [h1]
    · simp only [hT, ite_false]
      symm
      refine Finset.sum_eq_zero (fun u _ => Finset.sum_eq_zero (fun v _ => ?_))
      by_cases h1 : blk u = X
      · subst h1; simp [hT]
      · simp [h1]
  rw [Finset.sum_congr rfl (fun X _ => step X), Finset.sum_comm]
  refine Finset.sum_congr rfl (fun u _ => ?_)
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun v _ => ?_)
  rw [Finset.sum_ite_eq]
  simp [Finset.mem_image_of_mem blk (Finset.mem_univ u)]

lemma bcost_diff {δ ε : Type*} [DecidableEq δ] [DecidableEq ε] (p : V → δ) (q : V → ε)
    (X : γ) :
    ∑ u, ∑ v, (if blk u = X ∧ blk v = X then d G p u v - d G q u v else 0) =
      bcost G p blk X - bcost G q blk X := by
  simp only [bcost, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))
  split_ifs <;> simp

/-- Cost of the child relative to parent `a`. -/
theorem cost_child_a (h : Compatible a b blk) (T : γ → Prop) [DecidablePred T] :
    cost G (child a b blk T) =
      cost G a + ∑ X ∈ univ.image blk, if T X then bcost G b blk X - bcost G a blk X else 0 := by
  have hpt : ∀ u v, d G (child a b blk T) u v =
      d G a u v + if blk u = blk v ∧ T (blk u) then d G b u v - d G a u v else 0 := by
    intro u v
    rw [d_child G h T u v]
    by_cases hT : T (blk u)
    · by_cases hbl : blk u = blk v
      · rw [ite_eq_left hT, ite_eq_left ⟨hbl, hT⟩]; ring
      · simp [hT, hbl, d_cross G h u v hbl]
    · simp [hT]
  have hs : ∀ X, (if T X then bcost G b blk X - bcost G a blk X else 0) =
      (if T X then ∑ u, ∑ v, (if blk u = X ∧ blk v = X then d G b u v - d G a u v else 0)
        else 0) := by
    intro X; rw [bcost_diff]
  simp only [cost]
  rw [Finset.sum_congr rfl (fun X _ => hs X), sum_blocks]
  simp only [hpt, Finset.sum_add_distrib]

/-- Cost of the child relative to parent `b`. -/
theorem cost_child_b (h : Compatible a b blk) (T : γ → Prop) [DecidablePred T] :
    cost G (child a b blk T) =
      cost G b + ∑ X ∈ univ.image blk,
        if ¬ T X then bcost G a blk X - bcost G b blk X else 0 := by
  have hpt : ∀ u v, d G (child a b blk T) u v =
      d G b u v + if blk u = blk v ∧ ¬ T (blk u) then d G a u v - d G b u v else 0 := by
    intro u v
    rw [d_child G h T u v]
    by_cases hT : T (blk u)
    · simp [hT]
    · by_cases hbl : blk u = blk v
      · rw [ite_eq_right hT, ite_eq_left ⟨hbl, hT⟩]; ring
      · simp [hT, hbl, d_cross G h u v hbl]
  have hs : ∀ X, (if ¬ T X then bcost G a blk X - bcost G b blk X else 0) =
      (if ¬ T X then ∑ u, ∑ v, (if blk u = X ∧ blk v = X then d G a u v - d G b u v else 0)
        else 0) := by
    intro X; rw [bcost_diff]
  simp only [cost]
  rw [Finset.sum_congr rfl (fun X _ => hs X), sum_blocks (T := fun X => ¬ T X)]
  simp only [hpt, Finset.sum_add_distrib]

/-- Partition crossover: taking in every block the parent with the smaller block
cost gives a child that is at least as good as both parents. -/
theorem crossover_le (h : Compatible a b blk) :
    cost G (child a b blk (fun X => bcost G b blk X < bcost G a blk X)) ≤ cost G a ∧
    cost G (child a b blk (fun X => bcost G b blk X < bcost G a blk X)) ≤ cost G b := by
  constructor
  · rw [cost_child_a G h]
    have : ∑ X ∈ univ.image blk, (if bcost G b blk X < bcost G a blk X then
        bcost G b blk X - bcost G a blk X else 0) ≤ 0 := by
      apply Finset.sum_nonpos
      intro X _
      by_cases hX : bcost G b blk X < bcost G a blk X
      · simp only [hX, ite_true]; linarith
      · simp [hX]
    linarith
  · rw [cost_child_b G h]
    have : ∑ X ∈ univ.image blk, (if ¬ bcost G b blk X < bcost G a blk X then
        bcost G a blk X - bcost G b blk X else 0) ≤ 0 := by
      apply Finset.sum_nonpos
      intro X _
      by_cases hX : bcost G b blk X < bcost G a blk X
      · simp [hX]
      · simp only [hX, not_false_eq_true, ite_true]; linarith [not_lt.mp hX]
    linarith

end CC
