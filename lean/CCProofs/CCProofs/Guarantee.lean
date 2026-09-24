/-
Copyright (c) 2026 Quang-Vinh Dang. All rights reserved.
Released under the BSD 3-Clause license as described in the file LICENSE.
Authors: Quang-Vinh Dang
-/
import CCProofs.Crossover

/-!
# The approximation guarantee of the memetic search

The search keeps a population of clusterings.  Every operator replaces a
member by a clustering that is not more expensive than the member it replaces:
annealing returns the best state it visited (the start included), localized
annealing undoes a step whose total change is positive, partition crossover is
never worse than either parent (`crossover_le`), and offspring enter the
population only by replacing a member that is not cheaper.

* `Reach`: the populations obtainable by such replacements.
* `reach_le`: after any number of steps, some member costs at most as much as
  any given member of the initial population.
* `expected_guarantee`: if, for every outcome of the random seeds, the output
  costs at most the Pivot clustering of that outcome, and Pivot is a
  `ρ`-approximation in expectation (`ρ = 3`, Ailon, Charikar and Newman, JACM
  2008; used as a hypothesis), the output is a `ρ`-approximation in expectation.
-/

open Finset

namespace CC

section Population

variable {β : Type*} {P : Type*} [DecidableEq P] (f : β → ℤ)

/-- Populations reachable from `pop` by replacing a member with a clustering
that costs at most as much. -/
inductive Reach (pop : P → β) : (P → β) → Prop
  | refl : Reach pop pop
  | step {q : P → β} (j : P) (x : β) :
      Reach pop q → f x ≤ f (q j) → Reach pop (Function.update q j x)

/-- The best member never gets worse. -/
theorem reach_le {pop q : P → β} (h : Reach f pop q) (i : P) : ∃ j, f (q j) ≤ f (pop i) := by
  induction h with
  | refl => exact ⟨i, le_refl _⟩
  | step j x _ hx ih =>
    obtain ⟨k, hk⟩ := ih
    by_cases hkj : k = j
    · subst hkj
      exact ⟨k, by rw [Function.update_self]; exact le_trans hx hk⟩
    · exact ⟨k, by rw [Function.update_of_ne hkj]; exact hk⟩

/-- The output (a cheapest member of the final population) costs at most any
member of the initial population. -/
theorem output_le {pop q : P → β} (h : Reach f pop q) (out : β)
    (hout : ∀ j, f out ≤ f (q j)) (i : P) : f out ≤ f (pop i) := by
  obtain ⟨j, hj⟩ := reach_le f h i
  exact le_trans (hout j) hj

end Population

section Operators

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α β γ : Type*} [DecidableEq α] [DecidableEq β] [DecidableEq γ]

/-- Annealing returns the cheapest state of its trajectory, which starts at the
input; the result is never worse than the input. -/
theorem best_of_trajectory_le (traj : List (V → α)) (start : V → α) (best : V → α)
    (hmem : start ∈ traj) (hbest : ∀ s ∈ traj, cost G best ≤ cost G s) :
    cost G best ≤ cost G start := hbest start hmem

/-- A sequence of moves whose changes sum to a non-positive number (the
acceptance test of the localized annealing) does not increase the cost. -/
theorem telescoping_le (traj : ℕ → V → α) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, (cost G (traj (i + 1)) - cost G (traj i)) ≤ 0) :
    cost G (traj n) ≤ cost G (traj 0) := by
  rw [Finset.sum_range_sub (fun i => cost G (traj i))] at hsum
  linarith

end Operators

section Expectation

/-- Pointwise domination transfers an approximation guarantee in expectation.
`p` is the distribution of the random seeds (a finite sample space). -/
theorem expected_guarantee {Ω : Type*} [Fintype Ω] (p : Ω → ℝ) (hp : ∀ ω, 0 ≤ p ω)
    (out seed : Ω → ℝ) (hdom : ∀ ω, out ω ≤ seed ω) (ρ opt : ℝ)
    (hseed : ∑ ω, p ω * seed ω ≤ ρ * opt) :
    ∑ ω, p ω * out ω ≤ ρ * opt := by
  have : ∑ ω, p ω * out ω ≤ ∑ ω, p ω * seed ω :=
    Finset.sum_le_sum (fun ω _ => mul_le_mul_of_nonneg_left (hdom ω) (hp ω))
  linarith

/-- The guarantee of the whole search: for every outcome `ω`, the initial
population contains a clustering no more expensive than the Pivot clustering
`piv ω` (the Pivot seed itself, or its local-search improvement), every step is a
non-worsening replacement, and the output is a cheapest final member.  If Pivot
is a `ρ`-approximation in expectation, so is the output. -/
theorem search_guarantee {Ω β P : Type*} [Fintype Ω] [DecidableEq P]
    (p : Ω → ℝ) (hp : ∀ ω, 0 ≤ p ω) (f : β → ℤ)
    (pop0 popT : Ω → P → β) (out piv : Ω → β) (i0 : Ω → P)
    (hreach : ∀ ω, Reach f (pop0 ω) (popT ω))
    (hseed : ∀ ω, f (pop0 ω (i0 ω)) ≤ f (piv ω))
    (hout : ∀ ω j, f (out ω) ≤ f (popT ω j))
    (ρ opt : ℝ) (hpivot : ∑ ω, p ω * (f (piv ω) : ℝ) ≤ ρ * opt) :
    ∑ ω, p ω * (f (out ω) : ℝ) ≤ ρ * opt := by
  apply expected_guarantee p hp (fun ω => (f (out ω) : ℝ)) (fun ω => (f (piv ω) : ℝ)) _ ρ opt
    hpivot
  intro ω
  have h1 := output_le f (hreach ω) (out ω) (hout ω) (i0 ω)
  exact_mod_cast le_trans h1 (hseed ω)

end Expectation

end CC
