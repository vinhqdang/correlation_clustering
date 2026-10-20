import Mathlib

/-!
# Grid Tiling with `≤` and the block encoding used by the comparator gadget

Paper: `sections/hardness.tex` (statement of Grid Tiling with `≤`), `app_hardness.tex`
(Sec. "The reduction", Proposition "Tunnel release", last step of the proof of (⇒)),
`notes/wave10_hardness_proof.md` Sec. 2.1, 2.4, 2.5.

Convention (exactly the paper's, i.e. the "transposed" one): an instance is `k`, `n` and non-empty sets
`S i j ⊆ [n]²` (`1 ≤ i,j ≤ k`).  A solution chooses `s i j = (a_ij, b_ij) ∈ S i j` with
`a_ij ≤ a_{i,j+1}` (first coordinate, along a ROW, `j ↦ j+1`) and
`b_ij ≤ b_{i+1,j}` (second coordinate, along a COLUMN, `i ↦ i+1`).

Everything in this file is proved without `sorry` and without axioms beyond the Lean core ones.
-/

namespace P3.Hard

open Finset

/-! ## 1. The problem -/

/-- An instance of Grid Tiling with `≤`: `S i j ⊆ [n] × [n]`, non-empty. -/
structure GTInstance where
  k : ℕ
  n : ℕ
  S : Fin k → Fin k → Finset (ℕ × ℕ)
  sub : ∀ i j, S i j ⊆ Icc 1 n ×ˢ Icc 1 n
  nonempty : ∀ i j, (S i j).Nonempty

/-- A solution (paper convention): `a_ij ≤ a_{i,j+1}` and `b_ij ≤ b_{i+1,j}`. -/
def GTInstance.IsSolution (I : GTInstance) (s : Fin I.k → Fin I.k → ℕ × ℕ) : Prop :=
  (∀ i j, s i j ∈ I.S i j) ∧
  (∀ i j j', j'.val = j.val + 1 → (s i j).1 ≤ (s i j').1) ∧
  (∀ i i' j, i'.val = i.val + 1 → (s i j).2 ≤ (s i' j).2)

def GTInstance.Solvable (I : GTInstance) : Prop := ∃ s, I.IsSolution s

/-- The same problem in the other coordinate convention (first coordinate along columns,
second along rows). -/
def GTInstance.IsSolutionAlt (I : GTInstance) (s : Fin I.k → Fin I.k → ℕ × ℕ) : Prop :=
  (∀ i j, s i j ∈ I.S i j) ∧
  (∀ i i' j, i'.val = i.val + 1 → (s i j).1 ≤ (s i' j).1) ∧
  (∀ i j j', j'.val = j.val + 1 → (s i j).2 ≤ (s i j').2)

/-- Exchange the two coordinates of every set. -/
def GTInstance.swap (I : GTInstance) : GTInstance where
  k := I.k
  n := I.n
  S i j := (I.S i j).image Prod.swap
  sub i j := by
    intro x hx
    obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hx
    have := I.sub i j hy
    simp only [Finset.mem_product, Finset.mem_Icc] at this ⊢
    exact ⟨this.2, this.1⟩
  nonempty i j := ((I.nonempty i j).image _)

/-- Exchanging the coordinates turns the paper's convention into the other one (and back). -/
theorem solvable_iff_swap_alt (I : GTInstance) :
    I.Solvable ↔ ∃ s, I.swap.IsSolutionAlt s := by
  constructor
  · rintro ⟨s, hs1, hs2, hs3⟩
    refine ⟨fun i j => (s i j).swap, ?_, ?_, ?_⟩
    · intro i j
      exact Finset.mem_image.2 ⟨s i j, hs1 i j, rfl⟩
    · intro i i' j h; exact hs3 i i' j h
    · intro i j j' h; exact hs2 i j j' h
  · rintro ⟨s, hs1, hs2, hs3⟩
    refine ⟨fun i j => (s i j).swap, ?_, ?_, ?_⟩
    · intro i j
      obtain ⟨y, hy, hy2⟩ := Finset.mem_image.1 (hs1 i j)
      have : (s i j).swap = y := by rw [← hy2]; exact Prod.swap_swap y
      show (s i j).swap ∈ I.S i j
      rw [this]; exact hy
    · intro i j j' h; exact hs3 i j j' h
    · intro i i' j h; exact hs2 i i' j h

/-! ## 2. Pure OR-clause family -/

/-- the clause of threshold `t`: `a ≤ t` or `b > t`. -/
def Clause (a b t : ℕ) : Prop := a ≤ t ∨ t < b

/-- **Comparator equivalence (block level).**  For `a, b ∈ [n]`:
`a ≤ b` iff for every threshold `t ∈ [1, n-1]` the clause `(a ≤ t) ∨ (b > t)` holds. -/
theorem le_iff_forall_clause (n a b : ℕ) (ha : a ∈ Icc 1 n) (hb : b ∈ Icc 1 n) :
    a ≤ b ↔ ∀ t ∈ Icc 1 (n - 1), Clause a b t := by
  simp only [Finset.mem_Icc, Clause] at *
  constructor
  · intro h t _
    by_cases hta : a ≤ t
    · exact Or.inl hta
    · exact Or.inr (by omega)
  · intro h
    by_contra hab
    have := h b ⟨hb.1, by omega⟩
    omega

/-! ## 3. Block encoding by separator thresholds -/

/-- Group of the slot `(α, s)` (`α` = block, `s` = slot inside the block), `α, s ∈ [n]`:
`u = (α-1)(n+1) + s`. -/
def slotGroup (n α s : ℕ) : ℕ := (α - 1) * (n + 1) + s

/-- Separator group `u_t = t (n+1)` after block `t`. -/
def sepGroup (n t : ℕ) : ℕ := t * (n + 1)

/-- Number of groups `U = 1 + n(n+1)` (group `0` is an end pad). -/
def numGroups (n : ℕ) : ℕ := 1 + n * (n + 1)

/-- chain positions: group `u` consists of the positions `2u, 2u+1`. A chain whose cut tunnel is in
group `u` sits in position `2u` (the other choice `2u+1` behaves identically). -/
def groupPos (u : ℕ) : ℕ := 2 * u

section arith
variable {n α s β s' t : ℕ}

theorem slot_lt_sep (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (_ht : t ∈ Icc 1 n) :
    slotGroup n α s < sepGroup n t ↔ α ≤ t := by
  simp only [Finset.mem_Icc, slotGroup, sepGroup] at *
  constructor
  · intro h
    by_contra hlt
    have : t + 1 ≤ α := by omega
    have h2 : t * (n + 1) ≤ (α - 1) * (n + 1) := by
      apply Nat.mul_le_mul_right; omega
    omega
  · intro h
    have h2 : (α - 1) * (n + 1) + (n + 1) ≤ t * (n + 1) := by
      have : α * (n + 1) ≤ t * (n + 1) := Nat.mul_le_mul_right _ h
      have e : α * (n + 1) = (α - 1) * (n + 1) + (n + 1) := by
        obtain ⟨a', rfl⟩ : ∃ a', α = a' + 1 := ⟨α - 1, by omega⟩
        simp [Nat.add_mul]
      omega
    omega

theorem sep_lt_slot (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (ht : t ∈ Icc 1 n) :
    sepGroup n t < slotGroup n α s ↔ t < α := by
  have := slot_lt_sep hα hs ht
  have hne : slotGroup n α s ≠ sepGroup n t := by
    intro h
    simp only [Finset.mem_Icc, slotGroup, sepGroup] at *
    -- slot ≡ s (mod n+1) with 1 ≤ s ≤ n, separator ≡ 0
    have h1 : ((α - 1) * (n + 1) + s) % (n + 1) = s := by
      rw [Nat.mul_comm, Nat.mul_add_mod]; exact Nat.mod_eq_of_lt (by omega)
    have h2 : (t * (n + 1)) % (n + 1) = 0 := Nat.mul_mod_left _ _
    rw [h] at h1; omega
  constructor
  · intro h
    by_contra hle
    have : α ≤ t := by omega
    have := (slot_lt_sep hα hs ht).2 this
    omega
  · intro h
    by_contra hle
    have : slotGroup n α s ≤ sepGroup n t := by omega
    have hlt : slotGroup n α s < sepGroup n t := lt_of_le_of_ne this hne
    have := (slot_lt_sep hα hs ht).1 hlt
    omega

/-- slot groups are never separator groups, never group `0`. -/
theorem slot_ne_sep (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (ht : t ∈ Icc 1 n) :
    slotGroup n α s ≠ sepGroup n t := by
  intro h
  have h1 := (slot_lt_sep hα hs ht)
  have h2 := (sep_lt_slot hα hs ht)
  rw [h] at h1 h2
  simp at h1 h2
  omega

theorem slot_pos (hs : s ∈ Icc 1 n) : 1 ≤ slotGroup n α s := by
  simp only [Finset.mem_Icc, slotGroup] at *; omega

theorem slot_lt_numGroups (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) :
    slotGroup n α s < numGroups n := by
  simp only [Finset.mem_Icc, slotGroup, numGroups] at *
  have : (α - 1) * (n + 1) + s ≤ (n - 1) * (n + 1) + n := by
    have : (α - 1) * (n + 1) ≤ (n - 1) * (n + 1) := Nat.mul_le_mul_right _ (by omega)
    omega
  have e : (n - 1) * (n + 1) + n + 1 = n * (n + 1) := by
    obtain ⟨a', rfl⟩ : ∃ a', n = a' + 1 := ⟨n - 1, by omega⟩
    simp [Nat.add_mul, Nat.mul_add]; ring
  omega

theorem sep_lt_numGroups (ht : t ∈ Icc 1 n) : sepGroup n t < numGroups n := by
  simp only [Finset.mem_Icc, sepGroup, numGroups] at *
  have : t * (n + 1) ≤ n * (n + 1) := Nat.mul_le_mul_right _ ht.2
  omega

/-- injectivity of the slot encoding. -/
theorem slotGroup_inj (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (hβ : β ∈ Icc 1 n)
    (hs' : s' ∈ Icc 1 n) (h : slotGroup n α s = slotGroup n β s') : α = β ∧ s = s' := by
  simp only [Finset.mem_Icc, slotGroup] at *
  have h1 : ((α - 1) * (n + 1) + s) % (n + 1) = s := by
    rw [Nat.mul_comm, Nat.mul_add_mod]; exact Nat.mod_eq_of_lt (by omega)
  have h2 : ((β - 1) * (n + 1) + s') % (n + 1) = s' := by
    rw [Nat.mul_comm, Nat.mul_add_mod]; exact Nat.mod_eq_of_lt (by omega)
  have hss : s = s' := by rw [← h1, ← h2, h]
  subst hss
  have : (α - 1) * (n + 1) = (β - 1) * (n + 1) := by omega
  have := Nat.eq_of_mul_eq_mul_right (by omega) this
  omega

/-- Coverage: every group `1 ≤ u ≤ n(n+1)` is a separator or a slot (and not both). -/
theorem group_cover (u : ℕ) (h1 : 1 ≤ u) (h2 : u ≤ n * (n + 1)) :
    (∃ t ∈ Icc 1 n, u = sepGroup n t) ∨ (∃ α ∈ Icc 1 n, ∃ s ∈ Icc 1 n, u = slotGroup n α s) := by
  have hq0 := Nat.zero_le (u / (n + 1))
  by_cases hm : u % (n + 1) = 0
  · left
    have hd := Nat.div_add_mod u (n + 1)
    rw [hm] at hd
    refine ⟨u / (n + 1), ?_, ?_⟩
    · simp only [Finset.mem_Icc]
      constructor
      · by_contra h0
        have h0' : u / (n + 1) = 0 := by omega
        rw [h0'] at hd
        omega
      · have : (n + 1) * (u / (n + 1)) ≤ (n + 1) * n := by
          rw [Nat.mul_comm (n + 1) n]; omega
        exact Nat.le_of_mul_le_mul_left this (Nat.succ_pos n)
    · simp only [sepGroup]
      rw [Nat.mul_comm]; omega
  · right
    refine ⟨u / (n + 1) + 1, ?_, u % (n + 1), ?_, ?_⟩
    · simp only [Finset.mem_Icc]
      have hlt : u < n * (n + 1) := by
        have := Nat.div_add_mod u (n + 1)
        have := Nat.mod_lt u (show n + 1 > 0 by omega)
        by_contra hge
        have : u = n * (n + 1) := by omega
        subst this
        simp at hm
      constructor
      · omega
      · have : u / (n + 1) < n := (Nat.div_lt_iff_lt_mul (by omega)).2 (by linarith [hlt, Nat.mul_comm n (n + 1)])
        omega
    · simp only [Finset.mem_Icc]
      have := Nat.mod_lt u (show n + 1 > 0 by omega)
      omega
    · simp only [slotGroup, Nat.add_sub_cancel]
      have := Nat.div_add_mod u (n + 1)
      rw [Nat.mul_comm]; omega

end arith

/-! ## 4. The comparator clause in terms of chain positions -/

/-- The comparator clause of threshold `t` between two chains at positions `pa`, `pb`:
the guest `y` can sit with `z_A` iff `pa ≤ 2u_t - 1`, with `z_B` iff `pb ≥ 2u_t + 2`
(paper: table of effective forbidden sets, X / Y). -/
def PosClause (n pa pb t : ℕ) : Prop :=
  pa + 1 ≤ 2 * sepGroup n t ∨ 2 * sepGroup n t + 2 ≤ pb

section
variable {n α s β s' : ℕ}

/-- position of a slot group versus a separator (only the position ≤ / ≥ matters). -/
theorem pos_le_sep_iff (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (ht : t ∈ Icc 1 n) :
    groupPos (slotGroup n α s) + 1 ≤ 2 * sepGroup n t ↔ α ≤ t := by
  rw [← slot_lt_sep hα hs ht]; unfold groupPos; omega

theorem sep_le_pos_iff (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (ht : t ∈ Icc 1 n) :
    2 * sepGroup n t + 2 ≤ groupPos (slotGroup n α s) ↔ t < α := by
  rw [← sep_lt_slot hα hs ht]; unfold groupPos; omega

/-- **Comparator equivalence (position level).**  If two chains sit in the slot groups of blocks
`α` and `β`, then all `n-1` clauses `PosClause` hold iff `α ≤ β`. -/
theorem posClauses_iff (hα : α ∈ Icc 1 n) (hs : s ∈ Icc 1 n) (hβ : β ∈ Icc 1 n)
    (hs' : s' ∈ Icc 1 n) :
    (∀ t ∈ Icc 1 (n - 1), PosClause n (groupPos (slotGroup n α s)) (groupPos (slotGroup n β s')) t)
      ↔ α ≤ β := by
  rw [le_iff_forall_clause n α β hα hβ]
  have key : ∀ t ∈ Icc 1 (n - 1),
      (PosClause n (groupPos (slotGroup n α s)) (groupPos (slotGroup n β s')) t ↔ Clause α β t) := by
    intro t ht
    have ht' : t ∈ Icc 1 n := by simp only [Finset.mem_Icc] at *; omega
    unfold PosClause Clause
    rw [pos_le_sep_iff hα hs ht', sep_le_pos_iff hβ hs' ht']
  exact ⟨fun h t ht => (key t ht).1 (h t ht), fun h t ht => (key t ht).2 (h t ht)⟩

end

/-! ## 5. Comparator conditions = Grid Tiling constraints -/

/-- position of chain `A_ij` for the tile `(a,b)`: slot group `(α = a, s = b)`. -/
def posA (t : ℕ × ℕ) (n : ℕ) : ℕ := groupPos (slotGroup n t.1 t.2)

/-- position of chain `B_ij` for the tile `(a,b)`: slot group `(α = b, s = a)`
(B-chains are ordered by `(b,a)`). -/
def posB (t : ℕ × ℕ) (n : ℕ) : ℕ := groupPos (slotGroup n t.2 t.1)

/-- The comparator conditions of the reduction for a choice `s` of tiles: for every pair of
horizontally adjacent A-chains and every threshold, the clause; for every pair of vertically
adjacent B-chains and every threshold, the clause. -/
def ComparatorConditions (I : GTInstance) (s : Fin I.k → Fin I.k → ℕ × ℕ) : Prop :=
  (∀ i j j', j'.val = j.val + 1 → ∀ t ∈ Icc 1 (I.n - 1),
      PosClause I.n (posA (s i j) I.n) (posA (s i j') I.n) t) ∧
  (∀ i i' j, i'.val = i.val + 1 → ∀ t ∈ Icc 1 (I.n - 1),
      PosClause I.n (posB (s i j) I.n) (posB (s i' j) I.n) t)

/-- **Combinatorial core of Prop. "Tunnel release" (⇒), last step.**  For every choice of tiles
(`s i j ∈ S i j`), the comparator conditions hold iff the choice satisfies the Grid Tiling
constraints `a_ij ≤ a_{i,j+1}`, `b_ij ≤ b_{i+1,j}`. -/
theorem comparatorConditions_iff (I : GTInstance) (s : Fin I.k → Fin I.k → ℕ × ℕ)
    (hs : ∀ i j, s i j ∈ I.S i j) :
    ComparatorConditions I s ↔ I.IsSolution s := by
  have hmem : ∀ i j, (s i j).1 ∈ Icc 1 I.n ∧ (s i j).2 ∈ Icc 1 I.n := by
    intro i j
    have := I.sub i j (hs i j)
    simp only [Finset.mem_product] at this
    exact this
  unfold ComparatorConditions GTInstance.IsSolution
  constructor
  · rintro ⟨hrow, hcol⟩
    refine ⟨hs, ?_, ?_⟩
    · intro i j j' h
      have := hrow i j j' h
      exact (posClauses_iff (hmem i j).1 (hmem i j).2 (hmem i j').1 (hmem i j').2).1 this
    · intro i i' j h
      have := hcol i i' j h
      exact (posClauses_iff (hmem i j).2 (hmem i j).1 (hmem i' j).2 (hmem i' j).1).1 this
  · rintro ⟨-, hrow, hcol⟩
    refine ⟨?_, ?_⟩
    · intro i j j' h
      exact (posClauses_iff (hmem i j).1 (hmem i j).2 (hmem i j').1 (hmem i j').2).2 (hrow i j j' h)
    · intro i i' j h
      exact (posClauses_iff (hmem i j).2 (hmem i j).1 (hmem i' j).2 (hmem i' j).1).2 (hcol i i' j h)

/-- Corollary: Grid Tiling is solvable iff some choice of tiles satisfies the comparator
conditions. -/
theorem solvable_iff_comparators (I : GTInstance) :
    I.Solvable ↔ ∃ s, (∀ i j, s i j ∈ I.S i j) ∧ ComparatorConditions I s := by
  constructor
  · rintro ⟨s, h⟩
    exact ⟨s, h.1, (comparatorConditions_iff I s h.1).2 h⟩
  · rintro ⟨s, h1, h2⟩
    exact ⟨s, (comparatorConditions_iff I s h1).1 h2⟩

end P3.Hard
