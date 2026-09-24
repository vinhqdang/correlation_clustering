/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Moves

/-!
# Critical-clique contraction preserves the objective

Let `π : V → U` map every vertex to its class, where every class is a clique of
`G` and adjacency between different classes depends only on the classes
(`G.Adj u v ↔ R (π u) (π v)`); critical cliques have both properties.  The
contracted weighted instance has node sizes `s i = |π⁻¹(i)|` and weights
`w i j = s i * s j` if `R i j` and `0` otherwise.  For every clustering `c` of
the classes, the cost of the expanded clustering `c ∘ π` on `G` equals the
weighted cost of `c` (`cost_contract`).  Together with `optimal_keeps_twins`
this shows that optimising the contracted instance optimises `G`.
-/

open Finset

namespace CC

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]
variable {U : Type*} [Fintype U] [DecidableEq U]

/-- Size of the class `i`. -/
def sz (π : V → U) (i : U) : ℤ := ((univ.filter fun v => π v = i).card : ℤ)

/-- The contracted weighted instance. -/
def contracted (π : V → U) (R : U → U → Prop) [DecidableRel R] (hR : ∀ i j, R i j ↔ R j i) :
    WInst U where
  s := sz π
  w := fun i j => if R i j then sz π i * sz π j else 0
  w_symm := by
    intro i j
    by_cases h : R i j
    · have h' : R j i := (hR i j).mp h
      simp [h, h', mul_comm]
    · have h' : ¬ R j i := fun q => h ((hR i j).mpr q)
      simp [h, h']

omit [DecidableEq V] in
lemma sum_comp_eq (π : V → U) (F : U → ℤ) : ∑ u, F (π u) = ∑ i, sz π i * F i := by
  rw [← Finset.sum_fiberwise (s := univ) (g := π) (f := fun u => F (π u))]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Finset.sum_congr rfl (fun u hu => by rw [(Finset.mem_filter.mp hu).2])]
  simp [sz, Finset.sum_const]

theorem cost_contract (π : V → U) (R : U → U → Prop) [DecidableRel R]
    (hR : ∀ i j, R i j ↔ R j i)
    (hclique : ∀ u v, u ≠ v → π u = π v → G.Adj u v)
    (hcross : ∀ u v, π u ≠ π v → (G.Adj u v ↔ R (π u) (π v))) (c : U → α) :
    cost G (c ∘ π) = costW (contracted π R hR) c := by
  -- the disagreement of a vertex pair depends only on the classes of its ends
  set q : U → U → ℤ := fun i j =>
    if i = j then 0 else if (R i j ↔ c i = c j) then 0 else 1 with hq
  have hpt : ∀ u v, d G (c ∘ π) u v = q (π u) (π v) := by
    intro u v
    by_cases huv : u = v
    · subst huv; simp [d_self, hq]
    · by_cases hc : π u = π v
      · have ha := hclique u v huv hc
        simp [d, pd, huv, hq, hc, ha]
      · have e := hcross u v hc
        simp only [d, pd, huv, ite_false, hq, hc, Function.comp_apply, e]
  have hw : ∀ i j, sz π i * sz π j * q i j = wpd (contracted π R hR) i j (c i) (c j) := by
    intro i j
    simp only [hq, wpd, contracted]
    by_cases hij : i = j
    · simp [hij]
    · by_cases hr : R i j <;> by_cases hc : c i = c j <;> simp [hij, hr, hc]
  simp only [cost, hpt]
  rw [sum_comp_eq π (fun i => ∑ v, q i (π v))]
  simp only [costW]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [sum_comp_eq π (fun j => q i j), Finset.mul_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [← hw i j]
  ring

end CC
