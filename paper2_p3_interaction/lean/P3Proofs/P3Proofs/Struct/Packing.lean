import P3Proofs.Struct.Defs2

/-!
# Packings, soft/hard pairs, pair-maximality, and `OPT ≥ |H|`

A triple `P = (a, b, c)` is an induced `P₃` with centre `b` (`IsP3`).  `pr P` is the set of
the six ordered pairs inside `V(P)`; a packing is a finite set of induced `P₃`s with
pairwise disjoint `pr` (modification-disjoint; paper `sec:prelim`).  `Soft H x y` says
`xy ∈ ⋃ pairs(P)`; otherwise the pair is hard.

* `two_mul_card_le_cost`: every labelling of `G` has at least `|H|` disagreeing pairs
  (the lower-bound half of Lemma `decomp` / `fam`).
* HE / H2 (`hard_edges_close`, `hard_nonedge_no_common`): paper `lem:HE`.
-/

open Finset

namespace P3S

variable {V : Type*} [DecidableEq V]

/-- `(a, b, c)` is an induced `P₃` with centre `b`. -/
def IsP3 (G : SimpleGraph V) (P : V × V × V) : Prop :=
  G.Adj P.1 P.2.1 ∧ G.Adj P.2.1 P.2.2 ∧ ¬ G.Adj P.1 P.2.2 ∧ P.1 ≠ P.2.2

/-- The six ordered pairs inside `V(P)`. -/
def pr (P : V × V × V) : Finset (V × V) :=
  {(P.1, P.2.1), (P.2.1, P.1), (P.2.1, P.2.2), (P.2.2, P.2.1), (P.1, P.2.2), (P.2.2, P.1)}

/-- A modification-disjoint packing of induced `P₃`s. -/
structure IsPacking (G : SimpleGraph V) (H : Finset (V × V × V)) : Prop where
  p3 : ∀ P ∈ H, IsP3 G P
  disj : (H : Set (V × V × V)).PairwiseDisjoint pr

/-- The pair `(x, y)` is soft: it lies inside a member of `H`. -/
def Soft (H : Finset (V × V × V)) (x y : V) : Prop := ∃ P ∈ H, (x, y) ∈ pr P

/-- No induced `P₃` has all its pairs hard. -/
def PairMaximal (G : SimpleGraph V) (H : Finset (V × V × V)) : Prop :=
  ∀ Q, IsP3 G Q → ∃ x ∈ pr Q, Soft H x.1 x.2

/-- A family of induced `P₃`s hard over `H`: modification-disjoint and avoiding soft pairs. -/
structure HardFam (G : SimpleGraph V) (H M : Finset (V × V × V)) : Prop where
  p3 : ∀ Q ∈ M, IsP3 G Q
  disj : (M : Set (V × V × V)).PairwiseDisjoint pr
  hard : ∀ Q ∈ M, ∀ x ∈ pr Q, ¬ Soft H x.1 x.2

/-- `M` is an inclusion-maximal family hard over `H`. -/
def MaxHardFam (G : SimpleGraph V) (H M : Finset (V × V × V)) : Prop :=
  HardFam G H M ∧ ∀ Q, Q ∉ M → ¬ HardFam G H (insert Q M)

lemma mem_pr {P : V × V × V} {x : V × V} :
    x ∈ pr P ↔ x = (P.1, P.2.1) ∨ x = (P.2.1, P.1) ∨ x = (P.2.1, P.2.2) ∨
      x = (P.2.2, P.2.1) ∨ x = (P.1, P.2.2) ∨ x = (P.2.2, P.1) := by
  simp [pr]

lemma pr_symm {P : V × V × V} {x y : V} (h : (x, y) ∈ pr P) : (y, x) ∈ pr P := by
  simp only [mem_pr, Prod.mk.injEq] at h ⊢
  tauto

lemma Soft.symm {H : Finset (V × V × V)} {x y : V} (h : Soft H x y) : Soft H y x := by
  obtain ⟨P, hP, hx⟩ := h
  exact ⟨P, hP, pr_symm hx⟩

lemma Soft.mono {H H' : Finset (V × V × V)} (hsub : H ⊆ H') {x y : V} (h : Soft H x y) :
    Soft H' x y := by
  obtain ⟨P, hP, hx⟩ := h
  exact ⟨P, hsub hP, hx⟩

lemma fst_mem_pr (P : V × V × V) : (P.1, P.2.1) ∈ pr P := by simp [pr]

/-- Adding one induced `P₃` all of whose pairs are unpacked yields a packing. -/
theorem IsPacking.insert {G : SimpleGraph V} {H : Finset (V × V × V)} (hH : IsPacking G H)
    (Q : V × V × V) (hQ : IsP3 G Q) (hhard : ∀ x ∈ pr Q, ¬ Soft H x.1 x.2) :
    IsPacking G (Insert.insert Q H) := by
  refine ⟨?_, ?_⟩
  · intro P hP
    rcases Finset.mem_insert.mp hP with rfl | hP
    · exact hQ
    · exact hH.p3 P hP
  · rw [Finset.coe_insert]
    refine hH.disj.insert ?_
    intro j hj _
    rw [Finset.disjoint_left]
    intro x hxQ hxj
    exact hhard x hxQ ⟨j, hj, hxj⟩

/-- Union of a packing with a family hard over it is a packing. -/
theorem IsPacking.union {G : SimpleGraph V} {H M : Finset (V × V × V)} (hH : IsPacking G H)
    (hM : HardFam G H M) : IsPacking G (H ∪ M) := by
  refine ⟨?_, ?_⟩
  · intro P hP
    rcases Finset.mem_union.mp hP with h | h
    · exact hH.p3 P h
    · exact hM.p3 P h
  · rw [Finset.coe_union]
    refine Set.PairwiseDisjoint.union hH.disj hM.disj ?_
    intro i hi j hj _
    rw [Finset.disjoint_left]
    intro x hxi hxj
    exact hM.hard j hj x hxj ⟨i, hi, hxi⟩

lemma HardFam.disjoint {G : SimpleGraph V} {H M : Finset (V × V × V)} (hM : HardFam G H M) :
    Disjoint H M := by
  rw [Finset.disjoint_left]
  intro P hH hMem
  exact hM.hard P hMem _ (fst_mem_pr P) ⟨P, hH, fst_mem_pr P⟩

/-! ### HE and H2 -/

/-- Lemma HE (a). If `ab`, `bc` are hard edges (`a ≠ c`) then `ac ∈ E` or `ac` is soft. -/
theorem hard_edges_close {G : SimpleGraph V} {H : Finset (V × V × V)} (hmax : PairMaximal G H)
    (a b c : V) (hac : a ≠ c) (hab : G.Adj a b) (hbc : G.Adj b c)
    (hs1 : ¬ Soft H a b) (hs2 : ¬ Soft H b c) : G.Adj a c ∨ Soft H a c := by
  by_contra hno
  push Not at hno
  obtain ⟨x, hx, hxs⟩ := hmax (a, b, c) ⟨hab, hbc, hno.1, hac⟩
  have h1 := hs1; have h2 := hs2
  have h3 : ¬ Soft H b a := fun h => h1 h.symm
  have h4 : ¬ Soft H c b := fun h => h2 h.symm
  have h5 : ¬ Soft H c a := fun h => hno.2 h.symm
  rcases x with ⟨p, q⟩
  simp only [mem_pr, Prod.mk.injEq] at hx
  rcases hx with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
    simp_all

/-- Lemma H2 (b). If `uv` is a hard non-edge, no `z` has both `uz`, `zv` hard edges. -/
theorem hard_nonedge_no_common {G : SimpleGraph V} {H : Finset (V × V × V)}
    (hmax : PairMaximal G H) (u v : V) (huv : u ≠ v) (hn : ¬ G.Adj u v) (hs : ¬ Soft H u v)
    (z : V) (h1 : G.Adj u z) (h2 : G.Adj z v) (hs1 : ¬ Soft H u z) (hs2 : ¬ Soft H z v) :
    False := by
  obtain ⟨x, hx, hxs⟩ := hmax (u, z, v) ⟨h1, h2, hn, huv⟩
  have h3 : ¬ Soft H z u := fun h => hs1 h.symm
  have h4 : ¬ Soft H v z := fun h => hs2 h.symm
  have h5 : ¬ Soft H v u := fun h => hs h.symm
  rcases x with ⟨p, q⟩
  simp only [mem_pr, Prod.mk.injEq] at hx
  rcases hx with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
    simp_all

/-! ### Each member of a packing forces a disagreement -/

section Cost

variable [Fintype V] (G : SimpleGraph V) [DecidableRel G.Adj] {α : Type*} [DecidableEq α]

lemma sum_pr_ge_two (c : V → α) (P : V × V × V) (hP : IsP3 G P) :
    2 ≤ ∑ x ∈ pr P, d G c x.1 x.2 := by
  rcases P with ⟨a, b, e⟩
  obtain ⟨hab, hbc, hae, hne⟩ := hP
  simp only at hab hbc hae hne
  have n1 : a ≠ b := hab.ne
  have n2 : b ≠ e := hbc.ne
  have s1 : G.Adj b a := hab.symm
  have s2 : G.Adj e b := hbc.symm
  have s3 : ¬ G.Adj e a := fun h => hae h.symm
  have e1 : (a, b) ∉ ({(b, a), (b, e), (e, b), (a, e), (e, a)} : Finset (V × V)) := by
    simp [Prod.ext_iff, n1, n2, hne, n1.symm, n2.symm, hne.symm]
  have e2 : (b, a) ∉ ({(b, e), (e, b), (a, e), (e, a)} : Finset (V × V)) := by
    simp [Prod.ext_iff, n1, n2, hne, n1.symm, n2.symm, hne.symm]
  have e3 : (b, e) ∉ ({(e, b), (a, e), (e, a)} : Finset (V × V)) := by
    simp [Prod.ext_iff, n1, n2, hne, n1.symm, n2.symm, hne.symm]
  have e4 : (e, b) ∉ ({(a, e), (e, a)} : Finset (V × V)) := by
    simp [Prod.ext_iff, n1, n2, hne, n1.symm, n2.symm, hne.symm]
  have e5 : (a, e) ∉ ({(e, a)} : Finset (V × V)) := by
    simp [Prod.ext_iff, hne]
  unfold pr
  simp only
  rw [Finset.sum_insert e1, Finset.sum_insert e2, Finset.sum_insert e3, Finset.sum_insert e4,
    Finset.sum_insert e5, Finset.sum_singleton]
  simp only [d, pd, n1, n2, hne, n1.symm, n2.symm, hne.symm, if_false, hab, hbc, hae, s1, s2, s3]
  by_cases h1 : c a = c b <;> by_cases h2 : c b = c e <;> simp [h1, h2, eq_comm] <;>
    first | omega | (have := h1.trans h2; simp_all)

lemma d_sum_nonneg (c : V → α) (s : Finset (V × V)) : 0 ≤ ∑ x ∈ s, d G c x.1 x.2 :=
  Finset.sum_nonneg (fun x _ => d_nonneg G c x.1 x.2)

/-- Every labelling has at least `|H|` disagreeing pairs: `2 |H| ≤ cost`. -/
theorem two_mul_card_le_cost {H : Finset (V × V × V)} (hH : IsPacking G H) (c : V → α) :
    2 * (H.card : ℤ) ≤ cost G c := by
  have h1 : 2 * (H.card : ℤ) ≤ ∑ P ∈ H, ∑ x ∈ pr P, d G c x.1 x.2 := by
    have : ∑ P ∈ H, (2 : ℤ) ≤ ∑ P ∈ H, ∑ x ∈ pr P, d G c x.1 x.2 :=
      Finset.sum_le_sum (fun P hP => sum_pr_ge_two G c P (hH.p3 P hP))
    simpa [mul_comm] using this
  have h2 : ∑ P ∈ H, ∑ x ∈ pr P, d G c x.1 x.2 = ∑ x ∈ H.biUnion pr, d G c x.1 x.2 :=
    (Finset.sum_biUnion hH.disj).symm
  have h3 : ∑ x ∈ H.biUnion pr, d G c x.1 x.2 ≤ ∑ x ∈ (univ : Finset V) ×ˢ (univ : Finset V),
      d G c x.1 x.2 :=
    Finset.sum_le_sum_of_subset_of_nonneg (by intro x _; simp)
      (fun x _ _ => d_nonneg G c x.1 x.2)
  have h4 : cost G c = ∑ x ∈ (univ : Finset V) ×ˢ (univ : Finset V), d G c x.1 x.2 := by
    rw [Finset.sum_product]; rfl
  linarith

end Cost

end P3S
