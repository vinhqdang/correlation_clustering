/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Crossover
import CCProofs.Moves

/-!
# Partition crossover for pairwise-additive objectives

The memetic search applies partition crossover to the contracted weighted
instance.  The argument of `Crossover.lean` uses only that the objective is a sum
over pairs of a term that depends on whether the pair is joined; here it is
proved for every such objective and instantiated for `costW`.

* `pcost J S c = ∑ u, ∑ v, (if c u = c v then J u v else S u v)`: the pair
  `(u, v)` pays `J u v` when joined and `S u v` when separated.  Both `cost G`
  (`cost_eq_pcost`) and `costW I` (`costW_eq_pcost`) have this form.
* `pterm_child`, `cost_childP_a`, `cost_childP_b`, `crossoverP_le`: the analogues of
  `d_child`, `cost_child_a`, `cost_child_b` and `crossover_le`.
* `crossoverW_le`: the weighted partition crossover (block-wise choice of the
  parent with the smaller weighted block cost) is never worse than either parent.
-/

open Finset

namespace CC

section Generic

variable {V : Type*} [Fintype V] [DecidableEq V]
variable {α β γ δ : Type*} [DecidableEq α] [DecidableEq β] [DecidableEq γ] [DecidableEq δ]

omit [Fintype V] [DecidableEq V] in
/-- Term of the pair `(u, v)`: `J u v` when joined, `S u v` when separated. -/
def pterm (J S : V → V → ℤ) (c : V → δ) (u v : V) : ℤ := if c u = c v then J u v else S u v

omit [DecidableEq V] in
/-- A pairwise-additive objective. -/
def pcost (J S : V → V → ℤ) (c : V → δ) : ℤ := ∑ u, ∑ v, pterm J S c u v

omit [DecidableEq V] in
/-- Its part on the pairs inside block `X`. -/
def bpcost (J S : V → V → ℤ) (p : V → δ) (blk : V → γ) (X : γ) : ℤ :=
  ∑ u, ∑ v, if blk u = X ∧ blk v = X then pterm J S p u v else 0

variable {J S : V → V → ℤ} {a : V → α} {b : V → β} {blk : V → γ}

omit [Fintype V] [DecidableEq V] [DecidableEq γ] in
lemma pterm_cross (h : Compatible a b blk) (u v : V) (huv : blk u ≠ blk v) :
    pterm J S a u v = pterm J S b u v := by
  have ha : a u ≠ a v := fun e => huv (h.ha u v e)
  have hb : b u ≠ b v := fun e => huv (h.hb u v e)
  simp [pterm, ha, hb]

omit [Fintype V] [DecidableEq V] in
lemma pterm_child (h : Compatible a b blk) (T : γ → Prop) [DecidablePred T] (u v : V) :
    pterm J S (child a b blk T) u v = if T (blk u) then pterm J S b u v else pterm J S a u v := by
  by_cases hbl : blk u = blk v
  · by_cases hT : T (blk u)
    · have hT' : T (blk v) := hbl ▸ hT
      simp [pterm, child, hT, hT']
    · have hT' : ¬ T (blk v) := hbl ▸ hT
      simp [pterm, child, hT, hT']
  · have ha : a u ≠ a v := fun e => hbl (h.ha u v e)
    have hb : b u ≠ b v := fun e => hbl (h.hb u v e)
    have hc : child a b blk T u ≠ child a b blk T v := by
      unfold child
      by_cases h1 : T (blk u) <;> by_cases h2 : T (blk v) <;> simp [h1, h2, ha, hb]
    by_cases hT : T (blk u)
    · simp [pterm, hc, hb, hT]
    · simp [pterm, hc, ha, hT]

omit [DecidableEq V] [DecidableEq α] [DecidableEq β] in
lemma bpcost_diff {ε ζ : Type*} [DecidableEq ε] [DecidableEq ζ] (p : V → ε) (q : V → ζ)
    (X : γ) :
    ∑ u, ∑ v, (if blk u = X ∧ blk v = X then pterm J S p u v - pterm J S q u v else 0) =
      bpcost J S p blk X - bpcost J S q blk X := by
  simp only [bpcost, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))
  split_ifs <;> simp

omit [DecidableEq V] in
/-- Cost of the child relative to parent `a`. -/
theorem cost_childP_a (h : Compatible a b blk) (T : γ → Prop) [DecidablePred T] :
    pcost J S (child a b blk T) =
      pcost J S a + ∑ X ∈ univ.image blk,
        if T X then bpcost J S b blk X - bpcost J S a blk X else 0 := by
  have hpt : ∀ u v, pterm J S (child a b blk T) u v =
      pterm J S a u v +
        if blk u = blk v ∧ T (blk u) then pterm J S b u v - pterm J S a u v else 0 := by
    intro u v
    rw [pterm_child h T u v]
    by_cases hT : T (blk u)
    · by_cases hbl : blk u = blk v
      · rw [ite_eq_left hT, ite_eq_left ⟨hbl, hT⟩]; ring
      · simp [hT, hbl, pterm_cross (J := J) (S := S) h u v hbl]
    · simp [hT]
  have hs : ∀ X, (if T X then bpcost J S b blk X - bpcost J S a blk X else 0) =
      (if T X then ∑ u, ∑ v,
        (if blk u = X ∧ blk v = X then pterm J S b u v - pterm J S a u v else 0) else 0) := by
    intro X; rw [bpcost_diff]
  simp only [pcost]
  rw [Finset.sum_congr rfl (fun X _ => hs X), sum_blocks]
  simp only [hpt, Finset.sum_add_distrib]

omit [DecidableEq V] in
/-- Cost of the child relative to parent `b`. -/
theorem cost_childP_b (h : Compatible a b blk) (T : γ → Prop) [DecidablePred T] :
    pcost J S (child a b blk T) =
      pcost J S b + ∑ X ∈ univ.image blk,
        if ¬ T X then bpcost J S a blk X - bpcost J S b blk X else 0 := by
  have hpt : ∀ u v, pterm J S (child a b blk T) u v =
      pterm J S b u v +
        if blk u = blk v ∧ ¬ T (blk u) then pterm J S a u v - pterm J S b u v else 0 := by
    intro u v
    rw [pterm_child h T u v]
    by_cases hT : T (blk u)
    · simp [hT]
    · by_cases hbl : blk u = blk v
      · rw [ite_eq_right hT, ite_eq_left ⟨hbl, hT⟩]; ring
      · simp [hT, hbl, pterm_cross (J := J) (S := S) h u v hbl]
  have hs : ∀ X, (if ¬ T X then bpcost J S a blk X - bpcost J S b blk X else 0) =
      (if ¬ T X then ∑ u, ∑ v,
        (if blk u = X ∧ blk v = X then pterm J S a u v - pterm J S b u v else 0) else 0) := by
    intro X; rw [bpcost_diff]
  simp only [pcost]
  rw [Finset.sum_congr rfl (fun X _ => hs X), sum_blocks (T := fun X => ¬ T X)]
  simp only [hpt, Finset.sum_add_distrib]

omit [DecidableEq V] in
/-- Partition crossover for a pairwise-additive objective. -/
theorem crossoverP_le (h : Compatible a b blk) :
    pcost J S (child a b blk (fun X => bpcost J S b blk X < bpcost J S a blk X)) ≤
        pcost J S a ∧
    pcost J S (child a b blk (fun X => bpcost J S b blk X < bpcost J S a blk X)) ≤
        pcost J S b := by
  constructor
  · rw [cost_childP_a h]
    have : ∑ X ∈ univ.image blk, (if bpcost J S b blk X < bpcost J S a blk X then
        bpcost J S b blk X - bpcost J S a blk X else 0) ≤ 0 := by
      apply Finset.sum_nonpos
      intro X _
      by_cases hX : bpcost J S b blk X < bpcost J S a blk X
      · simp only [hX, ite_true]; linarith
      · simp [hX]
    linarith
  · rw [cost_childP_b h]
    have : ∑ X ∈ univ.image blk, (if ¬ bpcost J S b blk X < bpcost J S a blk X then
        bpcost J S a blk X - bpcost J S b blk X else 0) ≤ 0 := by
      apply Finset.sum_nonpos
      intro X _
      by_cases hX : bpcost J S b blk X < bpcost J S a blk X
      · simp [hX]
      · simp only [hX, not_false_eq_true, ite_true]; linarith [not_lt.mp hX]
    linarith

end Generic

section Instances

variable {V : Type*} [Fintype V] [DecidableEq V]
variable {α β γ δ : Type*} [DecidableEq α] [DecidableEq β] [DecidableEq γ] [DecidableEq δ]

/-- The unweighted objective is pairwise additive. -/
theorem cost_eq_pcost (G : SimpleGraph V) [DecidableRel G.Adj] (c : V → δ) :
    cost G c = pcost (fun u v => if u = v then 0 else if G.Adj u v then 0 else 1)
      (fun u v => if u = v then 0 else if G.Adj u v then 1 else 0) c := by
  simp only [cost, pcost]
  refine Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))
  by_cases huv : u = v
  · simp [d, pd, pterm, huv]
  · by_cases ha : G.Adj u v <;> by_cases hc : c u = c v <;> simp [d, pd, pterm, huv, ha, hc]

variable {U : Type*} [Fintype U] [DecidableEq U]

/-- Joined and separated pair costs of the weighted objective. -/
def Jw (I : WInst U) (i j : U) : ℤ := if i = j then 0 else I.s i * I.s j - I.w i j

/-- Separated pair cost of the weighted objective. -/
def Sw (I : WInst U) (i j : U) : ℤ := if i = j then 0 else I.w i j

/-- The weighted objective is pairwise additive. -/
theorem costW_eq_pcost (I : WInst U) (c : U → δ) : costW I c = pcost (Jw I) (Sw I) c := by
  simp only [costW, pcost]
  refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
  by_cases hij : i = j
  · simp [wpd, pterm, Jw, hij]
  · by_cases hc : c i = c j <;> simp [wpd, pterm, Jw, Sw, hij, hc]

/-- Weighted cost of the clustering `p` on the pairs of nodes inside block `X`. -/
def bcostW (I : WInst U) (p : U → δ) (blk : U → γ) (X : γ) : ℤ :=
  bpcost (Jw I) (Sw I) p blk X

/-- Weighted partition crossover: on the contracted instance, taking in every block
the parent with the smaller weighted block cost gives a child that is at least
as good as both parents. -/
theorem crossoverW_le (I : WInst U) {a : U → α} {b : U → β} {blk : U → γ}
    (h : Compatible a b blk) :
    costW I (child a b blk (fun X => bcostW I b blk X < bcostW I a blk X)) ≤ costW I a ∧
    costW I (child a b blk (fun X => bcostW I b blk X < bcostW I a blk X)) ≤ costW I b := by
  simp only [costW_eq_pcost, bcostW]
  exact crossoverP_le h

end Instances

end CC
