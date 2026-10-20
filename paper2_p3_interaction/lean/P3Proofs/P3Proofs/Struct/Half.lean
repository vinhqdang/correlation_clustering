import P3Proofs.Struct.Defs2

/-!
# Half-degree, diameter and far pairs (paper: `lem:half`, optimal.tex)

Let `c` be optimal (among labellings by `α`, with an unused label available; this is
automatic when `α` is infinite) and `C(v) = {w | c w = c v}`.

* `half_degree`: `|C(v)| - 1 ≤ 2 · |N(v) ∩ C(v)|` (Lemma half (a)).
* `cluster_common_neighbour`: two non-adjacent distinct vertices of one cluster have a
  common neighbour in that cluster (Lemma half (b)).
* `cluster_dist_le_two`: hence `G.dist a b ≤ 2` for `a, b` in a common cluster.
* `optimal_separates_far`: two distinct non-adjacent vertices without common neighbour
  (distance `≥ 3`) lie in different clusters of every optimum.
-/

open Finset

namespace P3S

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- **Half-degree lemma (a).** `d_C(v) ≥ (|C| - 1)/2`. -/
theorem half_degree (c : V → α) (hopt : IsOpt G c) (z : α) (hz : ∀ x, c x ≠ z) (v : V) :
    ((univ.filter fun w => w ≠ v ∧ c w = c v).card : ℤ) ≤
      2 * (univ.filter fun w => c w = c v ∧ G.Adj v w).card := by
  have h := cost_isolate G c v z hz
  have h2 := hopt (Function.update c v z)
  linarith

/-- Half-degree lemma for an infinite label type. -/
theorem half_degree' [Infinite α] (c : V → α) (hopt : IsOpt G c) (v : V) :
    ((univ.filter fun w => w ≠ v ∧ c w = c v).card : ℤ) ≤
      2 * (univ.filter fun w => c w = c v ∧ G.Adj v w).card := by
  obtain ⟨z, hz⟩ := exists_fresh c
  exact half_degree G c hopt z hz v

/-- **Half-degree lemma (b)** (diameter two inside a cluster). -/
theorem cluster_common_neighbour (c : V → α) (hopt : IsOpt G c) (z : α) (hz : ∀ x, c x ≠ z)
    (a b : V) (hab : a ≠ b) (hc : c a = c b) (hnadj : ¬ G.Adj a b) :
    ∃ x, x ≠ a ∧ x ≠ b ∧ c x = c a ∧ G.Adj a x ∧ G.Adj x b := by
  by_contra hno
  push Not at hno
  set C := univ.filter fun x => c x = c a with hC
  have ha : a ∈ C := by simp [hC]
  have hb : b ∈ C := by simp [hC, hc]
  set Ka := univ.filter fun x => c x = c a ∧ G.Adj a x with hKa
  set Kb := univ.filter fun x => c x = c b ∧ G.Adj b x with hKb
  have hdisj : Disjoint Ka Kb := by
    rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [hKa, hKb, Finset.mem_filter, Finset.mem_univ, true_and] at hx hx'
    have hxa : x ≠ a := hx.2.ne.symm
    have hxb : x ≠ b := hx'.2.ne.symm
    exact absurd hx'.2.symm (by simpa using hno x hxa hxb hx.1 hx.2)
  have hsub : Ka ∪ Kb ⊆ (C.erase a).erase b := by
    intro x hx
    simp only [Finset.mem_union, hKa, hKb, Finset.mem_filter, Finset.mem_univ, true_and] at hx
    simp only [Finset.mem_erase, hC, Finset.mem_filter, Finset.mem_univ, true_and]
    rcases hx with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · refine ⟨?_, h2.ne.symm, h1⟩
      rintro rfl; exact hnadj h2
    · refine ⟨h2.ne.symm, ?_, by rw [h1, hc]⟩
      rintro rfl; exact hnadj h2.symm
  have hcard : (Ka.card : ℤ) + Kb.card ≤ (C.card : ℤ) - 2 := by
    have h1 := Finset.card_le_card hsub
    rw [Finset.card_union_of_disjoint hdisj] at h1
    have h2 : ((C.erase a).erase b).card = C.card - 2 := by
      rw [Finset.card_erase_of_mem (Finset.mem_erase.mpr ⟨Ne.symm hab, hb⟩),
        Finset.card_erase_of_mem ha]
      omega
    have h3 : 2 ≤ C.card := by
      have : ({a, b} : Finset V) ⊆ C := by
        intro x hx; simp at hx; rcases hx with rfl | rfl; exact ha; exact hb
      have := Finset.card_le_card this
      rw [Finset.card_pair hab] at this
      exact this
    omega
  have da := half_degree G c hopt z hz a
  have db := half_degree G c hopt z hz b
  have nA : (univ.filter fun x => x ≠ a ∧ c x = c a) = C.erase a := by
    ext x; simp [hC, and_comm]
  have nB : (univ.filter fun x => x ≠ b ∧ c x = c b) = C.erase b := by
    ext x; simp [hC, hc, and_comm]
  rw [nA] at da
  rw [nB] at db
  have cu : ((C.erase a).card : ℤ) = C.card - 1 := by
    rw [Finset.card_erase_of_mem ha]; have : 1 ≤ C.card := Finset.card_pos.mpr ⟨a, ha⟩
    omega
  have cw : ((C.erase b).card : ℤ) = C.card - 1 := by
    rw [Finset.card_erase_of_mem hb]; have : 1 ≤ C.card := Finset.card_pos.mpr ⟨b, hb⟩
    omega
  rw [cu] at da
  rw [cw] at db
  have hKb' : Kb = (univ.filter fun x => c x = c a ∧ G.Adj b x) := by
    rw [hKb, hc]
  rw [← hKb'] at db
  have hKa' : (univ.filter fun w => c w = c a ∧ G.Adj a w) = Ka := rfl
  rw [hKa'] at da
  have hb' : (univ.filter fun w => c w = c b ∧ G.Adj b w) = Kb := rfl
  rw [← hc] at hb'
  linarith

/-- Two vertices of one cluster of an optimum are equal, adjacent, or at distance two
through a vertex of the cluster; in particular `G.dist a b ≤ 2`. -/
theorem cluster_dist_le_two (c : V → α) (hopt : IsOpt G c) (z : α) (hz : ∀ x, c x ≠ z)
    (a b : V) (hc : c a = c b) : G.dist a b ≤ 2 := by
  by_cases hab : a = b
  · subst hab; simp
  by_cases hadj : G.Adj a b
  · have := SimpleGraph.dist_le hadj.toWalk
    simp at this
    omega
  · obtain ⟨x, -, -, -, h1, h2⟩ := cluster_common_neighbour G c hopt z hz a b hab hc hadj
    have := SimpleGraph.dist_le (h1.toWalk.append h2.toWalk)
    simpa using this

/-- Far pairs (distinct, non-adjacent, no common neighbour) are separated in every optimum. -/
theorem optimal_separates_far (c : V → α) (hopt : IsOpt G c) (z : α) (hz : ∀ x, c x ≠ z)
    (u w : V) (huw : u ≠ w) (hnadj : ¬ G.Adj u w)
    (hcommon : ∀ x, ¬ (G.Adj u x ∧ G.Adj x w)) : c u ≠ c w := by
  intro hc
  obtain ⟨x, -, -, -, h1, h2⟩ := cluster_common_neighbour G c hopt z hz u w huw hc hnadj
  exact hcommon x ⟨h1, h2⟩

/-- Versions for an infinite label type. -/
theorem optimal_separates_far' [Infinite α] (c : V → α) (hopt : IsOpt G c)
    (u w : V) (huw : u ≠ w) (hnadj : ¬ G.Adj u w)
    (hcommon : ∀ x, ¬ (G.Adj u x ∧ G.Adj x w)) : c u ≠ c w := by
  obtain ⟨z, hz⟩ := exists_fresh c
  exact optimal_separates_far G c hopt z hz u w huw hnadj hcommon

end P3S
