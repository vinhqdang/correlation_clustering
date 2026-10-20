import Mathlib

/-!
# The thin one-hot chain `OH_m`, zero-excess (ZE) partitions, finite version of the chain lemma

Paper: `app_hardness.tex`, Sec. "The thin chain", Lemma "Chain lemma" (a),(b);
`notes/wave8_hardness.md` Sec. 3, `code/w8c_oh.py` (the construction), `code/w10_chain.py`
(`explicit_partition`), `code/a14_chain.py` (independent SAT enumeration).

## What is defined
* `PInst`: a "packing instance" on the vertices `0..n-1`: a list of members `(a,b,c)` (induced `P₃`
  with centre `b`: edges `ab, bc`, non-edge `ac`; the three pairs of a member are *soft*), and a list of
  extra hard edges; every other pair is a hard non-edge.
* `OH m`: the thin chain, vertex `s_k ↦ k` (`0 ≤ k ≤ 2m`), `e_j ↦ 2m+1+j` (`0 ≤ j ≤ m`),
  members exactly as in `oh_instance(m)` of `w8c_oh.py`.
* `IsZE I f`: the partition given by the labelling `f : ℕ → ℕ` is *zero-excess*:
  every hard pair is respected (hard edge together, hard non-edge apart) and every member has
  exactly one edit (cost `[a≁b] + [b≁c] + [a∼c] = 1`, i.e. states `L`, `R`, `M`).
* `Pa m a`: the explicit partition `P_a` (formula `explicit_partition` of `w10_chain.py`).

## What is proved (all with `decide +kernel`: the kernel evaluates the finite check; NO
`native_decide`, hence no `Lean.ofReduceBool` axiom; see `#print axioms` at the end of the file)
* `OH_valid`: the member list is a pair-disjoint packing of induced `P₃`'s, hard pairs are not soft,
  every vertex lies in at most 4 members (`m = 7..12`);
* for `m = 7,…,12` (the cases `m = 7,8,9` requested, `10,11,12` as a bonus):
  1. every `P_a` (`a < m`) is ZE;
  2. every ZE partition (arbitrary labelling `f : ℕ → ℕ`) induces on the `3m+2` vertices the same
     partition as exactly one `P_a`;
  3. the `P_a` are pairwise distinct as partitions;
  4. the spine word of every ZE partition is `R^a M L^(m-1-a)`.
  The proof of 2 is an exhaustive search (restricted-growth-string DFS with pruning) whose
  *completeness* is proved in general (`dfs_complete`) and whose finite run is checked by the kernel.

NOT proved here: the statement for all `m ≥ 7` (that is the transfer-matrix argument of
`Transfer.lean`, whose link to the graph is not formalised).
-/

namespace P3.Hard

/-! ## 1. Packing instances and ZE partitions -/

structure PInst where
  n : ℕ
  mems : List (ℕ × ℕ × ℕ)
  hardE : List (ℕ × ℕ)

/-- the unordered pair `{u,v}` equals `{x,y}`. -/
def pairEq (u v x y : ℕ) : Bool := (u == x && v == y) || (u == y && v == x)

namespace PInst

/-- soft pair: lies in a member. -/
def softB (I : PInst) (u v : ℕ) : Bool :=
  I.mems.any fun t => pairEq u v t.1 t.2.1 || pairEq u v t.2.1 t.2.2 || pairEq u v t.1 t.2.2

/-- edge of the graph: a member edge `ab`, `bc`, or a hard edge. -/
def edgeB (I : PInst) (u v : ℕ) : Bool :=
  I.mems.any (fun t => pairEq u v t.1 t.2.1 || pairEq u v t.2.1 t.2.2) ||
    I.hardE.any (fun e => pairEq u v e.1 e.2)

/-- the number of edits of the member `t = (a,b,c)` under the partition of `f`:
`[a≁b] + [b≁c] + [a∼c]`. -/
def edits (f : ℕ → ℕ) (t : ℕ × ℕ × ℕ) : ℕ :=
  (if f t.1 = f t.2.1 then 0 else 1) + (if f t.2.1 = f t.2.2 then 0 else 1) +
    (if f t.1 = f t.2.2 then 1 else 0)

/-- **Zero-excess partition** (labelling `f`): hard pairs respected, one edit per member. -/
def IsZE (I : PInst) (f : ℕ → ℕ) : Prop :=
  (∀ u v, u < v → v < I.n → I.softB u v = false → (I.edgeB u v = true ↔ f u = f v)) ∧
  (∀ t ∈ I.mems, edits f t = 1)

end PInst

open PInst

/-- Two labellings induce the same partition of `0..n-1`. -/
def SamePart (n : ℕ) (f g : ℕ → ℕ) : Prop :=
  ∀ u, u < n → ∀ v, v < n → (f u = f v ↔ g u = g v)

instance (n : ℕ) (f g : ℕ → ℕ) : Decidable (SamePart n f g) := by
  unfold SamePart; infer_instance

/-! ## 2. Constraint lists and the general DFS-completeness theorem -/

inductive Con
  | eq (u v : ℕ)
  | ne (u v : ℕ)
  | mem (a b c : ℕ)

namespace Con

def last : Con → ℕ
  | eq u v => max u v
  | ne u v => max u v
  | mem a b c => max a (max b c)

def holdsB (g : ℕ → ℕ) : Con → Bool
  | eq u v => g u == g v
  | ne u v => g u != g v
  | mem a b c => edits g (a, b, c) == 1

end Con

/-- the hard-pair constraints and the member constraints of an instance. -/
def PInst.cons (I : PInst) : List Con :=
  (List.range I.n).flatMap (fun v => (List.range v).filterMap (fun u =>
      if I.softB u v then none else if I.edgeB u v then some (Con.eq u v) else some (Con.ne u v))) ++
    I.mems.map (fun t => Con.mem t.1 t.2.1 t.2.2)

theorem isZE_iff (I : PInst) (f : ℕ → ℕ) :
    I.IsZE f ↔ ∀ c ∈ I.cons, c.holdsB f = true := by
  constructor
  · rintro ⟨h1, h2⟩ c hc
    simp only [PInst.cons, List.mem_append, List.mem_flatMap, List.mem_filterMap, List.mem_range,
      List.mem_map] at hc
    rcases hc with ⟨v, hv, u, hu, hcu⟩ | ⟨t, ht, rfl⟩
    · by_cases hs : I.softB u v = true
      · simp [hs] at hcu
      · have hs' : I.softB u v = false := by simpa using hs
        by_cases he : I.edgeB u v = true
        · simp [hs, he] at hcu
          subst hcu
          have := (h1 u v hu hv hs').1 he
          simp [Con.holdsB, this]
        · simp [hs, he] at hcu
          subst hcu
          have : ¬ f u = f v := fun h => he ((h1 u v hu hv hs').2 h)
          simp [Con.holdsB, this]
    · simpa [Con.holdsB] using h2 t ht
  · intro h
    refine ⟨?_, ?_⟩
    · intro u v huv hv hs
      have hmem : (if I.edgeB u v then Con.eq u v else Con.ne u v) ∈ I.cons := by
        simp only [PInst.cons, List.mem_append, List.mem_flatMap, List.mem_filterMap, List.mem_range,
          List.mem_map]
        left
        refine ⟨v, hv, u, huv, ?_⟩
        by_cases he : I.edgeB u v = true <;> simp [hs, he]
      have := h _ hmem
      by_cases he : I.edgeB u v = true
      · simp [he, Con.holdsB] at this ⊢
        exact this
      · have he' : I.edgeB u v = false := by simpa using he
        simp [he', Con.holdsB] at this ⊢
        exact this
    · intro t ht
      have hmem : Con.mem t.1 t.2.1 t.2.2 ∈ I.cons := by
        simp only [PInst.cons, List.mem_append, List.mem_map]
        right
        exact ⟨t, ht, rfl⟩
      have := h _ hmem
      simpa [Con.holdsB] using this

/-! ### canonical form of a labelling (first occurrence) -/

/-- label of `i` := the least index with the same label. -/
noncomputable def canon (f : ℕ → ℕ) (i : ℕ) : ℕ := Nat.find (⟨i, rfl⟩ : ∃ u, f u = f i)

theorem canon_spec (f : ℕ → ℕ) (i : ℕ) : f (canon f i) = f i :=
  Nat.find_spec (⟨i, rfl⟩ : ∃ u, f u = f i)

theorem canon_le (f : ℕ → ℕ) (i : ℕ) : canon f i ≤ i :=
  Nat.find_min' (⟨i, rfl⟩ : ∃ u, f u = f i) rfl

theorem canon_eq_iff (f : ℕ → ℕ) (i j : ℕ) : canon f i = canon f j ↔ f i = f j := by
  constructor
  · intro h
    have h1 := canon_spec f i
    have h2 := canon_spec f j
    rw [h] at h1
    exact h1.symm.trans h2
  · intro h
    apply le_antisymm
    · exact Nat.find_min' (⟨i, rfl⟩ : ∃ u, f u = f i) (by rw [canon_spec, h])
    · exact Nat.find_min' (⟨j, rfl⟩ : ∃ u, f u = f j) (by rw [canon_spec, h])

theorem canon_idem (f : ℕ → ℕ) (i : ℕ) : canon f (canon f i) = canon f i :=
  (canon_eq_iff f _ _).2 (canon_spec f i)

theorem edits_canon (f : ℕ → ℕ) (t : ℕ × ℕ × ℕ) : edits (canon f) t = edits f t := by
  unfold edits
  simp only [canon_eq_iff]

theorem holdsB_canon (f : ℕ → ℕ) (c : Con) : c.holdsB (canon f) = c.holdsB f := by
  cases c with
  | eq u v => simp [Con.holdsB, canon_eq_iff]
  | ne u v =>
    by_cases h : f u = f v
    · have h' := (canon_eq_iff f u v).2 h
      simp [Con.holdsB, h, h']
    · have h' : ¬ canon f u = canon f v := fun x => h ((canon_eq_iff f u v).1 x)
      simp [Con.holdsB, h, h']
  | mem a b c => simp [Con.holdsB, edits_canon]

theorem holdsB_congr (g g' : ℕ → ℕ) (c : Con) (h : ∀ i, i ≤ c.last → g i = g' i) :
    c.holdsB g = c.holdsB g' := by
  cases c with
  | eq u v =>
    simp only [Con.last] at h
    simp [Con.holdsB, h u (le_max_left _ _), h v (le_max_right _ _)]
  | ne u v =>
    simp only [Con.last] at h
    simp [Con.holdsB, h u (le_max_left _ _), h v (le_max_right _ _)]
  | mem a b c =>
    simp only [Con.last] at h
    have ha := h a (le_max_left _ _)
    have hb := h b (le_trans (le_max_left _ _) (le_max_right _ _))
    have hc := h c (le_trans (le_max_right _ _) (le_max_right _ _))
    simp [Con.holdsB, edits, ha, hb, hc]

/-! ### the search -/

/-- check of the constraints whose last index is `k` against a list `w` (labels of `0..`). -/
def okAt (cs : List Con) (k : ℕ) (w : List ℕ) : Bool :=
  (cs.filter (fun c => c.last == k)).all (fun c => c.holdsB (fun i => w.getD i 0))

/-- candidate labels of vertex `k`: `k` (a new class) or an existing representative `r < k`
with `w[r] = r`. -/
def cand (k : ℕ) (w : List ℕ) : List ℕ :=
  (List.range (k + 1)).filter (fun r => r == k || w.getD r 0 == r)

def step (cs : List Con) (k : ℕ) (ws : List (List ℕ)) : List (List ℕ) :=
  ws.flatMap fun w => ((cand k w).map (fun r => w ++ [r])).filter (okAt cs k)

/-- all restricted-growth labellings of `0..k-1` passing all constraints with last index `< k`. -/
def dfs (cs : List Con) : ℕ → List (List ℕ)
  | 0 => [[]]
  | k + 1 => step cs k (dfs cs k)

/-- the canonical labelling of `0..k-1` as a list. -/
noncomputable def canonList (f : ℕ → ℕ) (k : ℕ) : List ℕ := (List.range k).map (canon f)

theorem canonList_getD (f : ℕ → ℕ) (k i : ℕ) (hi : i < k) : (canonList f k).getD i 0 = canon f i := by
  simp [canonList, List.getD_eq_getElem?_getD, hi]

theorem canonList_succ (f : ℕ → ℕ) (k : ℕ) : canonList f (k + 1) = canonList f k ++ [canon f k] := by
  simp [canonList, List.range_succ]

/-- **Completeness of the search**: if `f` satisfies all constraints `cs`, then the canonical
labelling of the first `k` vertices is found at level `k`, for every `k`. -/
theorem dfs_complete (cs : List Con) (f : ℕ → ℕ) (hf : ∀ c ∈ cs, c.holdsB f = true) :
    ∀ k, canonList f k ∈ dfs cs k := by
  intro k
  induction k with
  | zero => simp [canonList, dfs]
  | succ k ih =>
    simp only [dfs, step, List.mem_flatMap]
    refine ⟨canonList f k, ih, ?_⟩
    rw [List.mem_filter]
    constructor
    · rw [List.mem_map]
      refine ⟨canon f k, ?_, canonList_succ f k |>.symm⟩
      simp only [cand, List.mem_filter, List.mem_range]
      refine ⟨Nat.lt_succ_of_le (canon_le f k), ?_⟩
      by_cases h : canon f k = k
      · simp [h]
      · have hlt : canon f k < k := lt_of_le_of_ne (canon_le f k) h
        rw [canonList_getD f k _ hlt, canon_idem]
        simp
    · unfold okAt
      rw [List.all_eq_true]
      intro c hc
      rw [List.mem_filter] at hc
      have hl : c.last = k := by simpa using hc.2
      have : c.holdsB (fun i => (canonList f (k + 1)).getD i 0) = c.holdsB (canon f) := by
        apply holdsB_congr
        intro i hi
        rw [canonList_getD f (k + 1) i (by omega)]
      rw [this, holdsB_canon]
      exact hf c hc.1

/-- Every ZE labelling induces (on `0..n-1`) the same partition as some labelling found by the
search. -/
theorem ze_in_dfs (I : PInst) (f : ℕ → ℕ) (hf : I.IsZE f) :
    ∃ w ∈ dfs I.cons I.n, ∀ u, u < I.n → ∀ v, v < I.n →
      (f u = f v ↔ w.getD u 0 = w.getD v 0) := by
  refine ⟨canonList f I.n, dfs_complete I.cons f ((isZE_iff I f).1 hf) I.n, ?_⟩
  intro u hu v hv
  rw [canonList_getD f _ _ hu, canonList_getD f _ _ hv, canon_eq_iff]

/-! ## 3. The chain `OH_m` -/

/-- vertex `s_k`. -/
def sV (k : ℕ) : ℕ := k
/-- vertex `e_j` (of `OH_m`). -/
def eV (m j : ℕ) : ℕ := 2 * m + 1 + j

/-- The thin one-hot chain `OH_m` (`m ≥ 7`), exactly `oh_instance(m)` of `code/w8c_oh.py`:
spine members `P_j = (s_2j, s_2j+1, s_2j+2)`, left cap, interior `G_j`, `F_j`, right cap, four hard
edges. -/
def OH (m : ℕ) : PInst where
  n := 3 * m + 2
  mems :=
    (List.range m).map (fun j => (sV (2 * j), sV (2 * j + 1), sV (2 * j + 2))) ++
    [(eV m 0, sV 2, eV m 1), (sV 1, eV m 1, sV 3), (eV m 0, sV 4, eV m 2)] ++
    (List.range (m - 5)).map (fun i => (eV m (i + 2), sV (2 * (i + 2) + 1), eV m (i + 2 + 2))) ++
    (List.range (m - 4)).map (fun i => (eV m (i + 2), sV (2 * (i + 2) + 2), eV m (i + 2 + 1))) ++
    [(sV (2 * m - 5), eV m (m - 3), sV (2 * m - 1)), (eV m (m - 2), sV (2 * m - 2), eV m m),
     (sV (2 * m - 4), eV m (m - 1), eV m m), (sV (2 * m - 3), eV m (m - 1), sV (2 * m)),
     (eV m (m - 2), eV m (m - 1), sV (2 * m - 1)), (sV (2 * m - 3), eV m m, sV (2 * m - 1))]
  hardE := [(sV 0, eV m 1), (eV m 0, sV 3), (sV (2 * m - 3), eV m (m - 2)), (sV (2 * m), eV m m)]

/-- the normalised list of all `3|H|` soft pairs. -/
def softPairs (I : PInst) : List (ℕ × ℕ) :=
  I.mems.flatMap fun t =>
    [(min t.1 t.2.1, max t.1 t.2.1), (min t.2.1 t.2.2, max t.2.1 t.2.2), (min t.1 t.2.2, max t.1 t.2.2)]

/-- Validity of a packing instance: all vertices in range, members pair-disjoint (the `3|H|` soft
pairs are distinct), `ab, bc` edges, `ac` non-edge, hard edges are not soft and in range. -/
def validB (I : PInst) : Bool :=
  I.mems.all (fun t => decide (t.1 < I.n ∧ t.2.1 < I.n ∧ t.2.2 < I.n) &&
      I.edgeB t.1 t.2.1 && I.edgeB t.2.1 t.2.2 && !I.edgeB t.1 t.2.2) &&
  decide (softPairs I).Nodup &&
  I.hardE.all (fun e => decide (e.1 < I.n ∧ e.2 < I.n) && !I.softB e.1 e.2)

/-- every vertex lies in at most 4 members (`c = 4` for the thin chain). -/
def degB (I : PInst) : Bool :=
  (List.range I.n).all fun v =>
    decide ((I.mems.filter (fun t => t.1 == v || t.2.1 == v || t.2.2 == v)).length ≤ 4)

theorem OH_valid : ∀ m ∈ [7, 8, 9, 10, 11, 12], validB (OH m) = true ∧ degB (OH m) = true := by
  decide +kernel

/-! ### The explicit partitions `P_a` -/

/-- label of `s_k` in `P_a`: cell index. -/
def sLab (a k : ℕ) : ℕ := if k % 2 = 1 then k / 2 else if k / 2 ≤ a then k / 2 else k / 2 - 1

/-- label of `e_j` in `P_a` (`explicit_partition` of `w10_chain.py`, with its cap exceptions). -/
def eLab (m a j : ℕ) : ℕ :=
  if j = 0 then 1 else if j = 1 then 0 else if j ≤ m - 2 then j
  else if j = m - 1 then (if a = m - 1 then m - 2 else m - 1) else m - 1

/-- the explicit ZE partition `P_a` as a labelling of the vertices `0..3m+1`. -/
def Pa (m a i : ℕ) : ℕ := if i ≤ 2 * m then sLab a i else eLab m a (i - (2 * m + 1))

/-! ### spine words -/

inductive St | L | R | M | X deriving DecidableEq, Repr

def stOf (f : ℕ → ℕ) (x y z : ℕ) : St :=
  if f x = f y ∧ f y = f z then St.M else if f x = f y then St.R else if f y = f z then St.L else St.X

def spineWord (m : ℕ) (f : ℕ → ℕ) : List St :=
  (List.range m).map fun j => stOf f (2 * j) (2 * j + 1) (2 * j + 2)

theorem spineWord_congr (m : ℕ) (f g : ℕ → ℕ) (h : SamePart (3 * m + 2) f g) :
    spineWord m f = spineWord m g := by
  unfold spineWord
  apply List.map_congr_left
  intro j hj
  rw [List.mem_range] at hj
  unfold stOf
  have h1 := h (2 * j) (by omega) (2 * j + 1) (by omega)
  have h2 := h (2 * j + 1) (by omega) (2 * j + 2) (by omega)
  simp only [h1, h2]

/-- The word `R^a M L^(m-1-a)`. -/
def oneHotWord (m a : ℕ) : List St := List.replicate a St.R ++ [St.M] ++ List.replicate (m - 1 - a) St.L

/-! ### Main finite statement -/

/-- The finite certificate for a given `m`: 1. every `P_a` is ZE; 2. every labelling found by the
search is a `P_a`; 3. the `P_a` are pairwise distinct; 4. spine word of `P_a` is `R^a M L^(m-1-a)`. -/
def Certificate (m : ℕ) : Prop :=
  (∀ a, a < m → ∀ c ∈ (OH m).cons, c.holdsB (Pa m a) = true) ∧
  (∀ w ∈ dfs (OH m).cons (OH m).n, ∃ a, a < m ∧ ∀ u, u < (OH m).n → ∀ v, v < (OH m).n →
      (w.getD u 0 = w.getD v 0 ↔ Pa m a u = Pa m a v)) ∧
  (∀ a, a < m → ∀ b, b < m → SamePart (OH m).n (Pa m a) (Pa m b) → a = b) ∧
  (∀ a, a < m → spineWord m (Pa m a) = oneHotWord m a)

theorem certificate_7 : Certificate 7 := by unfold Certificate; decide +kernel
theorem certificate_8 : Certificate 8 := by unfold Certificate; decide +kernel
theorem certificate_9 : Certificate 9 := by unfold Certificate; decide +kernel

/-- **Chain lemma, finite version.**  If `Certificate m` holds then
(i) every `P_a` is a ZE partition of `OH_m`,
(ii) every ZE partition of `OH_m` (arbitrary labelling `f`) induces on the vertices the same
partition as `P_a` for exactly one `a < m`,
(iii) its spine word is `R^a M L^(m-1-a)` for that `a`. -/
def ChainStatement (m : ℕ) : Prop :=
    (∀ a, a < m → (OH m).IsZE (Pa m a)) ∧
    (∀ f, (OH m).IsZE f → ∃! a, a < m ∧ SamePart (OH m).n f (Pa m a)) ∧
    (∀ f, (OH m).IsZE f → ∃ a, a < m ∧ spineWord m f = oneHotWord m a)

theorem chain_classification (m : ℕ) (hc : Certificate m) : ChainStatement m := by
  unfold ChainStatement
  obtain ⟨h1, h2, h3, h4⟩ := hc
  refine ⟨fun a ha => (isZE_iff _ _).2 (h1 a ha), ?_, ?_⟩
  · intro f hf
    obtain ⟨w, hw, hfw⟩ := ze_in_dfs _ f hf
    obtain ⟨a, ha, hwa⟩ := h2 w hw
    have hfa : SamePart (OH m).n f (Pa m a) := fun u hu v hv => by
      rw [hfw u hu v hv]; exact hwa u hu v hv
    refine ⟨a, ⟨ha, hfa⟩, ?_⟩
    rintro b ⟨hb, hfb⟩
    apply h3 b hb a ha
    intro u hu v hv
    rw [← hfb u hu v hv, hfa u hu v hv]
  · intro f hf
    obtain ⟨w, hw, hfw⟩ := ze_in_dfs _ f hf
    obtain ⟨a, ha, hwa⟩ := h2 w hw
    have hfa : SamePart (3 * m + 2) f (Pa m a) := fun u hu v hv => by
      rw [hfw u hu v hv]; exact hwa u hu v hv
    exact ⟨a, ha, by rw [spineWord_congr m f _ hfa]; exact h4 a ha⟩

theorem chain_7 : ChainStatement 7 := chain_classification 7 certificate_7
theorem chain_8 : ChainStatement 8 := chain_classification 8 certificate_8
theorem chain_9 : ChainStatement 9 := chain_classification 9 certificate_9

end P3.Hard
