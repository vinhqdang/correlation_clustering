import P3Proofs.Struct.Packing

/-!
# Maximalisation (paper: `lem:fam`, `lem:MAX`, prelim.tex)

`H` a packing, `M` a family hard over `H` (induced `P₃`s, pairwise modification-disjoint,
with all pairs unpacked).  `H ∪ M` is a packing (`IsPacking.union`, one-member version
`IsPacking.insert`), `OPT ≥ |H| + |M|` (`max_opt_ge`), the budget equivalence
`OPT ≤ |H| + ℓ ↔ |M| ≤ ℓ ∧ OPT ≤ |H ∪ M| + (ℓ - |M|)` holds (`max_budget_equiv`),
and an inclusion-maximal such `M` makes `H ∪ M` pair-maximal (`max_pairMaximal`).
`OPT ≤ k` is `Sol G α k` (some labelling with `cost ≤ 2 k`, i.e. at most `k` edits).
-/

open Finset

namespace P3S

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- `OPT ≥ |H|` for every packing `H` (Lemma `fam`, second sentence, for `Φ = ∅`). -/
theorem sol_ge_card {H : Finset (V × V × V)} (hH : IsPacking G H) {k : ℤ} (h : Sol G α k) :
    (H.card : ℤ) ≤ k := by
  obtain ⟨c, hc⟩ := h
  have := two_mul_card_le_cost G hH c
  linarith

/-- Lemma `fam`, general form / `MAX` (a): `OPT ≥ |H| + |M| = |H ∪ M|`. -/
theorem max_opt_ge {H M : Finset (V × V × V)} (hH : IsPacking G H) (hM : HardFam G H M)
    {k : ℤ} (h : Sol G α k) : (H.card : ℤ) + M.card ≤ k := by
  have := sol_ge_card G (hH.union hM) h
  rwa [Finset.card_union_of_disjoint hM.disjoint, Nat.cast_add] at this

/-- Lemma `fam`: in a YES instance a hard family has at most `ℓ` members; so
`|M| > ℓ` implies NO. -/
theorem hardFam_card_le {H M : Finset (V × V × V)} (hH : IsPacking G H) (hM : HardFam G H M)
    (ℓ : ℕ) (h : Sol G α ((H.card : ℤ) + ℓ)) : M.card ≤ ℓ := by
  have := max_opt_ge G hH hM h
  omega

/-- `MAX` (b): budget equivalence `ℓ ↦ ℓ - |M|`. -/
theorem max_budget_equiv {H M : Finset (V × V × V)} (hH : IsPacking G H) (hM : HardFam G H M)
    (ℓ : ℕ) :
    Sol G α ((H.card : ℤ) + ℓ) ↔
      (M.card ≤ ℓ ∧ Sol G α (((H ∪ M).card : ℤ) + ((ℓ : ℤ) - M.card))) := by
  have e : (((H ∪ M).card : ℤ) + ((ℓ : ℤ) - M.card)) = (H.card : ℤ) + ℓ := by
    rw [Finset.card_union_of_disjoint hM.disjoint, Nat.cast_add]; ring
  rw [e]
  constructor
  · intro h; exact ⟨hardFam_card_le G hH hM ℓ h, h⟩
  · intro h; exact h.2

/-- `MAX` (c): `H ∪ M` is a pair-maximal packing when `M` is inclusion-maximal hard over `H`. -/
theorem max_pairMaximal {H M : Finset (V × V × V)} (hmax : MaxHardFam G H M) :
    PairMaximal G (H ∪ M) := by
  intro Q hQ
  by_contra hno
  push Not at hno
  have hQM : Q ∉ M := by
    intro hm
    exact hno _ (fst_mem_pr Q) ⟨Q, Finset.mem_union_right _ hm, fst_mem_pr Q⟩
  apply hmax.2 Q hQM
  refine ⟨?_, ?_, ?_⟩
  · intro P hP
    rcases Finset.mem_insert.mp hP with rfl | hP
    · exact hQ
    · exact hmax.1.p3 P hP
  · rw [Finset.coe_insert]
    refine hmax.1.disj.insert ?_
    intro j hj _
    rw [Finset.disjoint_left]
    intro x hxQ hxj
    exact hno x hxQ ⟨j, Finset.mem_union_right _ hj, hxj⟩
  · intro P hP x hx hs
    rcases Finset.mem_insert.mp hP with rfl | hP
    · exact hno x hx (hs.mono Finset.subset_union_left)
    · exact hmax.1.hard P hP x hx hs

/-- `MAX` (c) for one added member: an induced `P₃` with all pairs unpacked extends the
packing and the lower bound `OPT ≥ |H| + 1` follows. -/
theorem insert_opt_ge {H : Finset (V × V × V)} (hH : IsPacking G H) (Q : V × V × V)
    (hQ : IsP3 G Q) (hhard : ∀ x ∈ pr Q, ¬ Soft H x.1 x.2) {k : ℤ} (h : Sol G α k) :
    (H.card : ℤ) + 1 ≤ k := by
  have hp := hH.insert Q hQ hhard
  have hQH : Q ∉ H := fun hm => hhard _ (fst_mem_pr Q) ⟨Q, hm, fst_mem_pr Q⟩
  have := sol_ge_card G hp h
  rwa [Finset.card_insert_of_notMem hQH, Nat.cast_add, Nat.cast_one] at this

end P3S
