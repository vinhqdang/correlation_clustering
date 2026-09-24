/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Lagrangian
import CCProofs.Distance
import CCProofs.Relabel

/-!
# The Lagrangian bound on the distance-two support

The relaxation of `ccbench/blockdual.py` has one variable per pair at distance at
most two; the far pairs (non-adjacent, no common neighbour) are fixed apart.
`Psupp G` is the set of ordered pairs of distinct vertices that are not far.

* `cost_eq_supp`: for a clustering that separates every far pair (`SepFar`),
  `cost = ∑ₚ∈P (wₚ xₚ - min(0, wₚ))`: a far pair is negative and separated, so it
  costs nothing, and the constant `K` of `cost_eq_affine` cancels against
  `-∑ₚ∈P min(0, wₚ)` up to the far pairs.
* `far_dual_bound`: for rows supported on `P` that the clustering satisfies and
  `y ≥ 0`,
  `bᵀy + ∑ₚ∈P (min(0, wₚ - (Aᵀy)ₚ) - min(0, wₚ)) ≤ cost`.
  This is the quantity evaluated by `experiments/check_certificate.py` (there over
  unordered pairs, here over ordered pairs, so both sides are doubled).
* `sepFar_of_optimal`: an optimal clustering with an unused label separates the
  far pairs (`optimal_separates_far`).
* `far_dual_bound_optimal`, `far_dual_le_opt`: for rows valid for every
  far-separating clustering the bound is at most the optimum; the second form
  needs no unused label and bounds the cost of every clustering of any label type.
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- `u` and `v` are far: distinct, non-adjacent and without a common neighbour. -/
def far (u v : V) : Prop := u ≠ v ∧ ¬ G.Adj u v ∧ ∀ w, ¬ (G.Adj u w ∧ G.Adj w v)

instance far.decidable (u v : V) : Decidable (far G u v) := by
  unfold far; infer_instance

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] [DecidableEq α] in
/-- The clustering `c` separates every far pair. -/
def SepFar (c : V → α) : Prop := ∀ u v, far G u v → c u ≠ c v

/-- The support `P`: ordered pairs of distinct vertices at distance at most two. -/
def Psupp : Finset (V × V) := univ.filter fun p => p.1 ≠ p.2 ∧ ¬ far G p.1 p.2

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
lemma far_symm {u v : V} (h : far G u v) : far G v u :=
  ⟨Ne.symm h.1, fun e => h.2.1 e.symm, fun w hw => h.2.2 w ⟨hw.2.symm, hw.1.symm⟩⟩

lemma d_eq_supp (c : V → α) (hc : SepFar G c) (p : V × V) :
    (d G c p.1 p.2 : ℝ) =
      if p ∈ Psupp G then wcoef G p * sep c p - min 0 (wcoef G p) else 0 := by
  obtain ⟨u, v⟩ := p
  by_cases huv : u = v
  · subst huv; simp [d_self, Psupp]
  by_cases hf : far G u v
  · have hcu := hc u v hf
    have hna := hf.2.1
    simp [d, pd, Psupp, huv, hf, hcu, hna]
  · have hmem : (u, v) ∈ Psupp G := by simp [Psupp, huv, hf]
    simp only [hmem, ite_true]
    by_cases ha : G.Adj u v <;> by_cases hcu : c u = c v <;>
      simp [d, pd, wcoef, sep, huv, ha, hcu]

/-- The cost of a far-separating clustering, written on the support `P`. -/
theorem cost_eq_supp (c : V → α) (hc : SepFar G c) :
    (cost G c : ℝ) = ∑ p ∈ Psupp G, (wcoef G p * sep c p - min 0 (wcoef G p)) := by
  have h : (cost G c : ℝ) = ∑ p : V × V, (d G c p.1 p.2 : ℝ) := by
    simp only [cost, Int.cast_sum]
    rw [← Finset.univ_product_univ, Finset.sum_product]
  rw [h, Finset.sum_congr rfl (fun p _ => d_eq_supp G c hc p), Finset.sum_ite_mem,
    Finset.univ_inter]

/-- Weak duality on the distance-two support, for one far-separating clustering. -/
theorem far_dual_bound {R : Type*} [Fintype R] (A : R → V × V → ℝ) (b : R → ℝ)
    (y : R → ℝ) (hy : ∀ r, 0 ≤ y r) (hA : ∀ r p, p ∉ Psupp G → A r p = 0)
    (c : V → α) (hc : SepFar G c) (hvalid : ∀ r, b r ≤ ∑ p, A r p * sep c p) :
    ∑ r, y r * b r +
        ∑ p ∈ Psupp G, (min 0 (wcoef G p - ∑ r, y r * A r p) - min 0 (wcoef G p)) ≤
      (cost G c : ℝ) := by
  set a : V × V → ℝ := fun p => ∑ r, y r * A r p with ha
  have h2 : ∑ r, y r * b r ≤ ∑ p, a p * sep c p := by
    calc ∑ r, y r * b r ≤ ∑ r, y r * ∑ p, A r p * sep c p :=
          Finset.sum_le_sum (fun r _ => mul_le_mul_of_nonneg_left (hvalid r) (hy r))
      _ = ∑ p, a p * sep c p := by
          simp only [ha, Finset.mul_sum, Finset.sum_mul]
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl (fun p _ => Finset.sum_congr rfl (fun r _ => ?_))
          ring
  have h2' : ∑ p, a p * sep c p = ∑ p ∈ Psupp G, a p * sep c p := by
    symm
    apply Finset.sum_subset (Finset.subset_univ _)
    intro p _ hp
    simp [ha, hA _ p hp]
  have h1 : ∀ p, min 0 (wcoef G p - a p) ≤ (wcoef G p - a p) * sep c p := by
    intro p
    have hx := sep_mem c p
    rcases le_or_gt 0 (wcoef G p - a p) with h | h
    · rw [min_eq_left h]; exact mul_nonneg h hx.1
    · rw [min_eq_right h.le]; nlinarith [hx.2]
  have h3 : ∑ p ∈ Psupp G, min 0 (wcoef G p - a p) ≤
      ∑ p ∈ Psupp G, (wcoef G p - a p) * sep c p :=
    Finset.sum_le_sum (fun p _ => h1 p)
  have e : ∑ p ∈ Psupp G, (wcoef G p * sep c p - min 0 (wcoef G p)) =
      ∑ p ∈ Psupp G, (wcoef G p - a p) * sep c p + ∑ p ∈ Psupp G, a p * sep c p -
        ∑ p ∈ Psupp G, min 0 (wcoef G p) := by
    rw [← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun p _ => ?_)
    ring
  rw [cost_eq_supp G c hc, e, Finset.sum_sub_distrib]
  linarith

/-- An optimal clustering with an unused label separates every far pair. -/
theorem sepFar_of_optimal (c : V → α) (hopt : ∀ c' : V → α, cost G c ≤ cost G c')
    (z : α) (hz : ∀ x, c x ≠ z) : SepFar G c := by
  intro u v ⟨huv, hna, hcom⟩
  exact optimal_separates_far G c hopt z hz u v huv hna
    (fun x hx => hcom x ⟨hx.1, hx.2.symm⟩)

/-- The bound is at most the cost of an optimal clustering (with an unused label),
hence at most the cost of every clustering. -/
theorem far_dual_bound_optimal {R : Type*} [Fintype R] (A : R → V × V → ℝ) (b : R → ℝ)
    (y : R → ℝ) (hy : ∀ r, 0 ≤ y r) (hA : ∀ r p, p ∉ Psupp G → A r p = 0)
    (hvalid : ∀ c : V → α, SepFar G c → ∀ r, b r ≤ ∑ p, A r p * sep c p)
    (c : V → α) (hopt : ∀ c' : V → α, cost G c ≤ cost G c') (z : α) (hz : ∀ x, c x ≠ z)
    (c' : V → α) :
    ∑ r, y r * b r +
        ∑ p ∈ Psupp G, (min 0 (wcoef G p - ∑ r, y r * A r p) - min 0 (wcoef G p)) ≤
      (cost G c' : ℝ) := by
  have hs := sepFar_of_optimal G c hopt z hz
  have h := far_dual_bound G A b y hy hA c hs (hvalid c hs)
  have h' : (cost G c : ℝ) ≤ cost G c' := by exact_mod_cast hopt c'
  linarith

/-- `LB ≤ OPT` without an unused label: if the rows are valid for every
far-separating clustering (clusterings are labelled by vertices, which loses
nothing by `cost_canon`), the bound is at most the cost of every clustering, of
any label type. -/
theorem far_dual_le_opt {R : Type*} [Fintype R] (A : R → V × V → ℝ) (b : R → ℝ)
    (y : R → ℝ) (hy : ∀ r, 0 ≤ y r) (hA : ∀ r p, p ∉ Psupp G → A r p = 0)
    (hvalid : ∀ c : V → V, SepFar G c → ∀ r, b r ≤ ∑ p, A r p * sep c p)
    (c' : V → α) :
    ∑ r, y r * b r +
        ∑ p ∈ Psupp G, (min 0 (wcoef G p - ∑ r, y r * A r p) - min 0 (wcoef G p)) ≤
      (cost G c' : ℝ) := by
  -- a cheapest labelling by vertices
  obtain ⟨c0, -, hc0⟩ := Finset.exists_min_image (univ : Finset (V → V)) (cost G)
    ⟨id, Finset.mem_univ _⟩
  have hmin : ∀ e : V → Option V, cost G c0 ≤ cost G e := by
    intro e
    rw [← cost_canon G e]
    exact hc0 _ (Finset.mem_univ _)
  have hmin' : cost G c0 ≤ cost G c' := by
    rw [← cost_canon G c']
    exact hc0 _ (Finset.mem_univ _)
  -- embedded into `Option V` it is optimal and leaves the label `none` unused
  have hopt : ∀ e : V → Option V, cost G (some ∘ c0) ≤ cost G e := by
    intro e
    rw [cost_comp_injective G some (Option.some_injective V)]
    exact hmin e
  have hs1 := sepFar_of_optimal G (some ∘ c0) hopt none (fun x => Option.some_ne_none _)
  have hs : SepFar G c0 := fun u v h e => hs1 u v h (congrArg some e)
  have h := far_dual_bound G A b y hy hA c0 hs (hvalid c0 hs)
  have h' : (cost G c0 : ℝ) ≤ cost G c' := by exact_mod_cast hmin'
  linarith

end CC
