import Mathlib

/-!
# Small merge (paper: `lem:SM`, app_counting.tex), abstracted

The paper's pieces are the components of the graph of kept hard edges of a solution and
`Up = U⁺` is the set of inserted hard non-edges.  What the counting argument uses is
only the following, and we take it as hypotheses:

* `piece : V → β` assigns a piece to every vertex;
* `Sf` is a symmetric "soft partner" relation with at most `D` partners inside `C`
  (`D ≥ max |Q(v)|`);
* (SM (i), the cross-pair lemma) two vertices of `C` in different pieces are soft partners
  or form an inserted pair (`s(x, y) ∈ Up`);
* `|Up| ≤ u`.

Then `sm_cross_count` is SM (ii) and `sm_union` is SM (iii):
`|C| ≤ max {D + u + 1, 2 D, 2 u}` when `C` meets at least two pieces.
The link "pieces = components of kept hard edges" (SM (i), `lem:cross`) is *not*
formalised here.
-/

open Finset

namespace P3S

variable {V β : Type*} [DecidableEq V] [DecidableEq β]
variable (Sf : V → V → Prop) [DecidableRel Sf] (Up : Finset (Sym2 V)) (D : ℕ)

/-- **SM (ii).** For disjoint `A`, `B` all of whose cross pairs are soft or inserted:
`|B| (|A| - D) ≤ |Up|`. -/
theorem sm_cross_count (A B : Finset V) (hAB : Disjoint A B)
    (hsoft : ∀ b ∈ B, (A.filter fun a => Sf b a).card ≤ D)
    (hcross : ∀ a ∈ A, ∀ b ∈ B, Sf b a ∨ s(a, b) ∈ Up) :
    (B.card : ℤ) * ((A.card : ℤ) - D) ≤ Up.card := by
  set P := (A ×ˢ B).filter (fun p => ¬ Sf p.2 p.1) with hP
  have hPcard : P.card = ∑ b ∈ B, (A.filter fun a => ¬ Sf b a).card := by
    rw [hP, Finset.card_filter, Finset.sum_product_right]
    refine Finset.sum_congr rfl (fun b _ => ?_)
    rw [Finset.card_filter]
  have hterm : ∀ b ∈ B, (A.card : ℤ) - D ≤ (A.filter fun a => ¬ Sf b a).card := by
    intro b hb
    have h1 := Finset.card_filter_add_card_filter_not (s := A) (fun a => Sf b a)
    have h2 := hsoft b hb
    omega
  have hge : (B.card : ℤ) * ((A.card : ℤ) - D) ≤ P.card := by
    rw [hPcard, Nat.cast_sum]
    calc (B.card : ℤ) * ((A.card : ℤ) - D) = ∑ b ∈ B, ((A.card : ℤ) - D) := by rw [Finset.sum_const, nsmul_eq_mul]
      _ ≤ _ := Finset.sum_le_sum hterm
  have hinj : Set.InjOn (fun p : V × V => s(p.1, p.2)) (P : Set (V × V)) := by
    intro p hp q hq hpq
    simp only [hP, Finset.coe_filter, Finset.mem_product, Set.mem_setOf_eq] at hp hq
    simp only [Sym2.eq_iff] at hpq
    rcases hpq with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact Prod.ext h1 h2
    · exfalso
      have : p.1 ∈ A ∧ p.1 ∈ B := ⟨hp.1.1, h1 ▸ hq.1.2⟩
      exact Finset.disjoint_left.mp hAB this.1 this.2
  have hmaps : ∀ p ∈ P, (fun p : V × V => s(p.1, p.2)) p ∈ Up := by
    intro p hp
    simp only [hP, Finset.mem_filter, Finset.mem_product] at hp
    rcases hcross p.1 hp.1.1 p.2 hp.1.2 with h | h
    · exact absurd h hp.2
    · exact h
  have := Finset.card_le_card_of_injOn _ hmaps hinj
  calc (B.card : ℤ) * ((A.card : ℤ) - D) ≤ P.card := hge
    _ ≤ Up.card := by exact_mod_cast this

/-- Corollary of SM (ii): if `B` is non-empty then `|A| ≤ D + |Up|`. -/
theorem sm_size_bound (A B : Finset V) (hAB : Disjoint A B) (hB : B.Nonempty)
    (hsoft : ∀ b ∈ B, (A.filter fun a => Sf b a).card ≤ D)
    (hcross : ∀ a ∈ A, ∀ b ∈ B, Sf b a ∨ s(a, b) ∈ Up) :
    A.card ≤ D + Up.card := by
  have h := sm_cross_count Sf Up D A B hAB hsoft hcross
  have hb : (1 : ℤ) ≤ B.card := by exact_mod_cast hB.card_pos
  by_contra hno
  push Not at hno
  have : (1 : ℤ) * ((A.card : ℤ) - D) ≤ (B.card : ℤ) * ((A.card : ℤ) - D) :=
    mul_le_mul_of_nonneg_right hb (by omega)
  omega

/-- **SM (iii).** If `C` meets two pieces then `|C| ≤ max {D + u + 1, 2 D, 2 u}`. -/
theorem sm_union (piece : V → β) (C : Finset V) (u : ℕ) (hUp : Up.card ≤ u)
    (hSf : ∀ x y, Sf x y → Sf y x)
    (hsoft : ∀ x ∈ C, (C.filter fun y => Sf x y).card ≤ D)
    (hcross : ∀ x ∈ C, ∀ y ∈ C, piece x ≠ piece y → Sf x y ∨ s(x, y) ∈ Up)
    (x₀ y₀ : V) (hx₀ : x₀ ∈ C) (hy₀ : y₀ ∈ C) (hne : piece x₀ ≠ piece y₀) :
    (C.card : ℤ) ≤ max ((D : ℤ) + u + 1) (max (2 * D) (2 * u)) := by
  by_cases hbig : ∃ z ∈ C, D < (C.filter fun y => piece y = piece z).card
  · -- Case `a > D`
    obtain ⟨z, hz, hzD⟩ := hbig
    set A := C.filter (fun y => piece y = piece z) with hA
    set B := C.filter (fun y => piece y ≠ piece z) with hB
    have hAB : Disjoint A B := by
      rw [Finset.disjoint_left]
      intro a ha hb
      simp only [hA, hB, Finset.mem_filter] at ha hb
      exact hb.2 ha.2
    have hBne : B.Nonempty := by
      by_cases h : piece x₀ = piece z
      · refine ⟨y₀, ?_⟩
        simp only [hB, Finset.mem_filter]
        exact ⟨hy₀, fun h' => hne (h.trans h'.symm)⟩
      · exact ⟨x₀, by simp only [hB, Finset.mem_filter]; exact ⟨hx₀, h⟩⟩
    have hcnt := sm_cross_count Sf Up D A B hAB
      (fun b hb => le_trans (Finset.card_le_card (Finset.filter_subset_filter _
        (Finset.filter_subset _ _))) (hsoft b (Finset.mem_filter.mp hb).1))
      (fun a ha b hb => by
        simp only [hA, hB, Finset.mem_filter] at ha hb
        rcases hcross a ha.1 b hb.1 (by rw [ha.2]; exact fun h => hb.2 h.symm) with h | h
        · exact Or.inl (hSf _ _ h)
        · exact Or.inr h)
    have hsum : A.card + B.card = C.card := Finset.card_filter_add_card_filter_not _
    have hr : (1 : ℤ) ≤ B.card := by exact_mod_cast hBne.card_pos
    have hs : (1 : ℤ) ≤ (A.card : ℤ) - D := by omega
    have hu : (Up.card : ℤ) ≤ u := by exact_mod_cast hUp
    have hsum' : (A.card : ℤ) + B.card = C.card := by exact_mod_cast hsum
    have : (B.card : ℤ) + ((A.card : ℤ) - D) ≤ u + 1 := by
      nlinarith [mul_nonneg (sub_nonneg.mpr hr) (sub_nonneg.mpr hs)]
    exact le_max_of_le_left (by linarith)
  · -- Case: every piece has at most `D` vertices in `C`
    push Not at hbig
    set W := (C ×ˢ C).filter (fun p => piece p.1 ≠ piece p.2 ∧ ¬ Sf p.1 p.2) with hW
    have hWcard : W.card = ∑ x ∈ C, (C.filter fun y => piece x ≠ piece y ∧ ¬ Sf x y).card := by
      rw [hW, Finset.card_filter, Finset.sum_product]
      refine Finset.sum_congr rfl (fun x _ => ?_)
      rw [Finset.card_filter]
    have hterm : ∀ x ∈ C, (C.card : ℤ) - 2 * D ≤
        (C.filter fun y => piece x ≠ piece y ∧ ¬ Sf x y).card := by
      intro x hx
      have hsub : C ⊆ ((C.filter fun y => piece x ≠ piece y ∧ ¬ Sf x y) ∪
          (C.filter fun y => piece y = piece x)) ∪ (C.filter fun y => Sf x y) := by
        intro y hy
        simp only [Finset.mem_union, Finset.mem_filter]
        by_cases h1 : piece y = piece x
        · exact Or.inl (Or.inr ⟨hy, h1⟩)
        · by_cases h2 : Sf x y
          · exact Or.inr ⟨hy, h2⟩
          · exact Or.inl (Or.inl ⟨hy, fun h => h1 h.symm, h2⟩)
      have h1 := Finset.card_le_card hsub
      have h2 := Finset.card_union_le (C.filter fun y => piece x ≠ piece y ∧ ¬ Sf x y)
        (C.filter fun y => piece y = piece x)
      have h3 := Finset.card_union_le ((C.filter fun y => piece x ≠ piece y ∧ ¬ Sf x y) ∪
          (C.filter fun y => piece y = piece x)) (C.filter fun y => Sf x y)
      have h4 := hbig x hx
      have h5 := hsoft x hx
      omega
    have hWge : (C.card : ℤ) * ((C.card : ℤ) - 2 * D) ≤ W.card := by
      rw [hWcard, Nat.cast_sum]
      calc (C.card : ℤ) * ((C.card : ℤ) - 2 * D) = ∑ x ∈ C, ((C.card : ℤ) - 2 * D) := by rw [Finset.sum_const, nsmul_eq_mul]
        _ ≤ _ := Finset.sum_le_sum hterm
    have hmaps : ∀ p ∈ W, (fun p : V × V => s(p.1, p.2)) p ∈ Up := by
      intro p hp
      simp only [hW, Finset.mem_filter, Finset.mem_product] at hp
      rcases hcross p.1 hp.1.1 p.2 hp.1.2 hp.2.1 with h | h
      · exact absurd h hp.2.2
      · exact h
    have hfib : ∀ e ∈ W.image (fun p : V × V => s(p.1, p.2)),
        (W.filter fun p => s(p.1, p.2) = e).card ≤ 2 := by
      intro e _
      induction e using Sym2.ind with
      | _ a b =>
        have hsub : (W.filter fun p => s(p.1, p.2) = s(a, b)) ⊆ {(a, b), (b, a)} := by
          intro p hp
          simp only [Finset.mem_filter, Sym2.eq_iff] at hp
          simp only [Finset.mem_insert, Finset.mem_singleton]
          rcases hp.2 with ⟨h1, h2⟩ | ⟨h1, h2⟩
          · exact Or.inl (Prod.ext h1 h2)
          · exact Or.inr (Prod.ext h1 h2)
        calc _ ≤ ({(a, b), (b, a)} : Finset (V × V)).card := Finset.card_le_card hsub
          _ ≤ 2 := Finset.card_le_two
    have h2 := Finset.card_le_mul_card_image (f := fun p : V × V => s(p.1, p.2)) W 2 hfib
    have h3 : (W.image (fun p : V × V => s(p.1, p.2))).card ≤ Up.card :=
      Finset.card_le_card (by
        intro e he
        obtain ⟨p, hp, rfl⟩ := Finset.mem_image.mp he
        exact hmaps p hp)
    have hu : (Up.card : ℤ) ≤ u := by exact_mod_cast hUp
    have hW2 : (W.card : ℤ) ≤ 2 * u := by
      have : W.card ≤ 2 * u := by omega
      exact_mod_cast this
    have key : (C.card : ℤ) * ((C.card : ℤ) - 2 * D) ≤ 2 * u := le_trans hWge hW2
    by_cases hcase : (C.card : ℤ) ≤ 2 * D
    · exact le_max_of_le_right (le_max_of_le_left hcase)
    · push Not at hcase
      have hC : (C.card : ℤ) * 1 ≤ (C.card : ℤ) * ((C.card : ℤ) - 2 * D) :=
        mul_le_mul_of_nonneg_left (by omega) (by positivity)
      exact le_max_of_le_right (le_max_of_le_right (by linarith))

end P3S
