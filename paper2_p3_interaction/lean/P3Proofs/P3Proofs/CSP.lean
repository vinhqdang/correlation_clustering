import P3Proofs.Defs

/-!
# Lower bound, the `ℓ = 0` CSP and the five states

* `needs_edit`, `card_le_size`, `card_le_OPT`: every solution edits a pair of every member;
  members have disjoint pairs, so `OPT ≥ |H|`.
* `obs_one` (Observation 3.1, first part): a solution of size `≤ |H|` has exactly one pair in
  each `pairsOf (P i)` and no other pair; it equals `editOf P x` for some `x : ι → Fin 3`.
* `csp_iff` (Observation 3.1): `OPT ≤ |H|` iff some `x : ι → Fin 3` makes `G △ editOf P x`
  free of induced `P₃`'s (the constraints of the CSP are one per triple, `indP3_agree`
  shows each depends only on the owners' variables); `opt_eq_card_iff` is the `OPT = |H|` form.
* `owners_clique` (Observation 3.2): the owners of the pairs of a triple pairwise share a
  vertex (or coincide), i.e. form a clique of the intersection graph.
* `five_states` (Lemma 3.5, "five states"): reduction to a finite check
  (`five_states_fin`, by `decide`).
-/

namespace P3

universe u
variable {V : Type u} {ι : Type*}

/-! ### Pairs of a member -/

/-- The three pairs of a member: `0 ↦ ab`, `1 ↦ bc`, `2 ↦ ac`. -/
def pairOf (t : V × V × V) : Fin 3 → Sym2 V :=
  ![s(t.1, t.2.1), s(t.2.1, t.2.2), s(t.1, t.2.2)]

theorem pairOf_mem (t : V × V × V) (k : Fin 3) : pairOf t k ∈ pairsOf t := by
  fin_cases k <;> simp [pairOf, pairsOf]

theorem exists_pairOf {t : V × V × V} {e : Sym2 V} (h : e ∈ pairsOf t) :
    ∃ k, pairOf t k = e := by
  simp only [pairsOf, Set.mem_insert_iff, Set.mem_singleton_iff] at h
  rcases h with h | h | h
  · exact ⟨0, by simp [pairOf, h]⟩
  · exact ⟨1, by simp [pairOf, h]⟩
  · exact ⟨2, by simp [pairOf, h]⟩

/-- The edit graph that flips, for each member `i`, the single pair `pairOf (P i) (x i)`. -/
def editOf (P : ι → V × V × V) (x : ι → Fin 3) : SimpleGraph V :=
  SimpleGraph.fromEdgeSet {e | ∃ i, e = pairOf (P i) (x i)}

theorem editOf_adj (P : ι → V × V × V) (x : ι → Fin 3) (u v : V) :
    (editOf P x).Adj u v ↔ (∃ i, s(u, v) = pairOf (P i) (x i)) ∧ u ≠ v := by
  simp [editOf, SimpleGraph.fromEdgeSet_adj]

/-! ### Lower bound `OPT ≥ |H|` -/

section LowerBound
variable {G : SimpleGraph V} {P : ι → V × V × V}

/-- Every solution flips at least one pair of every member. -/
theorem needs_edit (hP : IsPacking G P) {E : SimpleGraph V} (hE : IsSolution G E) (i : ι) :
    ∃ e ∈ pairsOf (P i), e ∈ E.edgeSet := by
  by_contra hcon
  simp only [not_exists, not_and] at hcon
  obtain ⟨hab, hbc, hac, hn⟩ := hP.p3 i
  apply hE (P i).1 (P i).2.1 (P i).2.2
  have h1 := hcon s((P i).1, (P i).2.1) (by simp [pairsOf])
  have h2 := hcon s((P i).2.1, (P i).2.2) (by simp [pairsOf])
  have h3 := hcon s((P i).1, (P i).2.2) (by simp [pairsOf])
  simp only [SimpleGraph.mem_edgeSet] at h1 h2 h3
  refine ⟨?_, ?_, hac, ?_⟩ <;> simp [edit_adj, hab, hbc, hn, h1, h2, h3]

theorem choose_edits (hP : IsPacking G P) {E : SimpleGraph V} (hE : IsSolution G E) :
    ∃ f : ι → Sym2 V, (∀ i, f i ∈ pairsOf (P i) ∧ f i ∈ E.edgeSet) ∧ Function.Injective f := by
  choose f hf using needs_edit hP hE
  refine ⟨f, fun i => ⟨(hf i).1, (hf i).2⟩, ?_⟩
  intro i j h
  by_contra hne
  refine Set.disjoint_left.1 (hP.disj hne) (hf i).1 ?_
  rw [h]; exact (hf j).1

/-- Members have disjoint pairs and each needs an edit: `|H| ≤ |S|`. -/
theorem card_le_size [Finite V] [Finite ι] (hP : IsPacking G P) {E : SimpleGraph V}
    (hE : IsSolution G E) : Nat.card ι ≤ size E := by
  obtain ⟨f, hf, hinj⟩ := choose_edits hP hE
  let g : ι → E.edgeSet := fun i => ⟨f i, (hf i).2⟩
  have hg : Function.Injective g := fun i j h => hinj (congrArg Subtype.val h)
  have := Nat.card_le_card_of_injective g hg
  rwa [Nat.card_coe_set_eq] at this

/-- `OPT ≥ |H|` for every modification-disjoint packing. -/
theorem card_le_OPT [Finite V] [Finite ι] (hP : IsPacking G P) : Nat.card ι ≤ OPT G := by
  obtain ⟨E, hE, hk⟩ := exists_opt G
  rw [← hk]
  exact card_le_size hP hE

end LowerBound

/-! ### Observation 3.1: the `ℓ = 0` CSP -/

section Obs
variable {G : SimpleGraph V} {P : ι → V × V × V}

/-- A solution of size at most `|H|` has exactly one pair in each member, no other pair, and
is the edit graph of a choice `x : ι → Fin 3`. -/
theorem obs_one [Finite V] [Finite ι] (hP : IsPacking G P) {E : SimpleGraph V}
    (hE : IsSolution G E) (h : size E ≤ Nat.card ι) :
    (∀ i, ∃ e, E.edgeSet ∩ pairsOf (P i) = {e}) ∧ (∀ e ∈ E.edgeSet, IsSoft P e) ∧
      ∃ x : ι → Fin 3, E = editOf P x := by
  obtain ⟨f, hf, hinj⟩ := choose_edits hP hE
  let g : ι → E.edgeSet := fun i => ⟨f i, (hf i).2⟩
  have hg : Function.Injective g := fun i j h => hinj (congrArg Subtype.val h)
  have hbij : Function.Bijective g := by
    refine hg.bijective_of_nat_card_le ?_
    rw [Nat.card_coe_set_eq]
    exact h
  have hsurj : ∀ e ∈ E.edgeSet, ∃ j, f j = e := by
    intro e he
    obtain ⟨j, hj⟩ := hbij.2 ⟨e, he⟩
    exact ⟨j, congrArg Subtype.val hj⟩
  refine ⟨?_, ?_, ?_⟩
  · intro i
    refine ⟨f i, ?_⟩
    ext e
    constructor
    · rintro ⟨heE, heP⟩
      obtain ⟨j, hj⟩ := hsurj e heE
      have hji : j = i := by
        by_contra hne
        exact Set.disjoint_left.1 (hP.disj hne) (hj ▸ (hf j).1) heP
      rw [← hj, hji]; rfl
    · intro he
      rw [Set.mem_singleton_iff] at he
      rw [he]; exact ⟨(hf i).2, (hf i).1⟩
  · intro e he
    obtain ⟨j, hj⟩ := hsurj e he
    exact ⟨j, hj ▸ (hf j).1⟩
  · choose x hx using fun i => exists_pairOf (hf i).1
    refine ⟨x, ?_⟩
    ext u v
    rw [editOf_adj]
    constructor
    · intro huv
      obtain ⟨j, hj⟩ := hsurj s(u, v) huv
      exact ⟨⟨j, by rw [hx j, hj]⟩, huv.ne⟩
    · rintro ⟨⟨i, hi⟩, -⟩
      rw [← SimpleGraph.mem_edgeSet, hi, hx i]
      exact (hf i).2

theorem size_editOf_le [Finite V] [Finite ι] (P : ι → V × V × V) (x : ι → Fin 3) :
    size (editOf P x) ≤ Nat.card ι := by
  have hsub : (editOf P x).edgeSet ⊆ (fun i => pairOf (P i) (x i)) '' Set.univ := by
    intro e he
    induction e using Sym2.ind with
    | h u v =>
      rw [SimpleGraph.mem_edgeSet, editOf_adj] at he
      obtain ⟨⟨i, hi⟩, -⟩ := he
      exact ⟨i, Set.mem_univ _, hi.symm⟩
  calc size (editOf P x) ≤ ((fun i => pairOf (P i) (x i)) '' Set.univ).ncard :=
        Set.ncard_le_ncard hsub (Set.toFinite _)
    _ ≤ (Set.univ : Set ι).ncard := Set.ncard_image_le (Set.toFinite _)
    _ = Nat.card ι := Set.ncard_univ ι

/-- Observation 3.1.  `(G, H, 0)` is a YES instance (`OPT ≤ |H|`) iff the CSP with variables
`x i ∈ {1,2,3}` (here `Fin 3`) and one constraint per vertex triple, "not an induced `P₃`
after flipping", is satisfiable. -/
theorem csp_iff [Finite V] [Finite ι] (hP : IsPacking G P) :
    OPT G ≤ Nat.card ι ↔
      ∃ x : ι → Fin 3, ∀ a b c, ¬ IndP3 (edit G (editOf P x)) a b c := by
  constructor
  · intro h
    obtain ⟨E, hE, hk⟩ := exists_opt G
    obtain ⟨-, -, x, hx⟩ := obs_one hP hE (hk ▸ h)
    exact ⟨x, by subst hx; exact hE⟩
  · rintro ⟨x, hx⟩
    exact (OPT_le (E := editOf P x) hx).trans (size_editOf_le P x)

/-- `OPT = |H|` iff some choice of one edit per member leaves no induced `P₃`. -/
theorem opt_eq_card_iff [Finite V] [Finite ι] (hP : IsPacking G P) :
    OPT G = Nat.card ι ↔
      ∃ x : ι → Fin 3, ∀ a b c, ¬ IndP3 (edit G (editOf P x)) a b c := by
  rw [← csp_iff hP]
  exact ⟨le_of_eq, fun h => le_antisymm h (card_le_OPT hP)⟩

/-! ### Scopes of the constraints (Observation 3.2) -/

/-- The owners of the soft pairs of a triple `t`. -/
def owners (P : ι → V × V × V) (t : V × V × V) : Set ι :=
  {i | ¬ Disjoint (pairsOf (P i)) (pairsOf t)}

theorem editOf_adj_of_mem (hP : IsPacking G P) (x : ι → Fin 3) {u v : V} (huv : u ≠ v) {j : ι}
    (hj : s(u, v) ∈ pairsOf (P j)) :
    (editOf P x).Adj u v ↔ s(u, v) = pairOf (P j) (x j) := by
  rw [editOf_adj]
  constructor
  · rintro ⟨⟨i, hi⟩, -⟩
    have hij : i = j := by
      by_contra hne
      refine Set.disjoint_left.1 (hP.disj hne) ?_ hj
      rw [hi]; exact pairOf_mem _ _
    subst hij; exact hi
  · intro h; exact ⟨⟨j, h⟩, huv⟩

theorem not_editOf_adj_of_hard (x : ι → Fin 3) {u v : V} (h : ¬ IsSoft P s(u, v)) :
    ¬ (editOf P x).Adj u v := by
  rw [editOf_adj]
  rintro ⟨⟨i, hi⟩, -⟩
  exact h ⟨i, by rw [hi]; exact pairOf_mem _ _⟩

/-- The status of a pair `uv` of the triple `t` in `G △ editOf P x` depends only on the
variables of the owners of `t`. -/
theorem edit_adj_agree (hP : IsPacking G P) {x y : ι → Fin 3} {t : V × V × V} {u v : V}
    (hm : s(u, v) ∈ pairsOf t) (hxy : ∀ i ∈ owners P t, x i = y i) :
    (edit G (editOf P x)).Adj u v ↔ (edit G (editOf P y)).Adj u v := by
  by_cases huv : u = v
  · subst huv; simp
  have : (editOf P x).Adj u v ↔ (editOf P y).Adj u v := by
    by_cases hs : IsSoft P s(u, v)
    · obtain ⟨j, hj⟩ := hs
      have hjo : j ∈ owners P t := Set.not_disjoint_iff.2 ⟨_, hj, hm⟩
      rw [editOf_adj_of_mem hP x huv hj, editOf_adj_of_mem hP y huv hj, hxy j hjo]
    · exact iff_of_false (not_editOf_adj_of_hard x hs) (not_editOf_adj_of_hard y hs)
  simp only [edit_adj, this]

/-- Constraint scope: whether a triple `(a, b, c)` whose pairs lie in `pairsOf t` is an
induced `P₃` of `G △ editOf P x` depends only on `x` restricted to `owners P t`. -/
theorem indP3_agree (hP : IsPacking G P) {x y : ι → Fin 3} {t : V × V × V} {a b c : V}
    (hab : s(a, b) ∈ pairsOf t) (hbc : s(b, c) ∈ pairsOf t) (hac : s(a, c) ∈ pairsOf t)
    (hxy : ∀ i ∈ owners P t, x i = y i) :
    IndP3 (edit G (editOf P x)) a b c ↔ IndP3 (edit G (editOf P y)) a b c := by
  unfold IndP3
  rw [edit_adj_agree hP hab hxy, edit_adj_agree hP hbc hxy, edit_adj_agree hP hac hxy]

/-- Two pairs inside the same triple share a vertex. -/
theorem pairs_share {t : V × V × V} {e e' : Sym2 V} (he : e ∈ pairsOf t) (he' : e' ∈ pairsOf t) :
    ∃ v, v ∈ e ∧ v ∈ e' := by
  simp only [pairsOf, Set.mem_insert_iff, Set.mem_singleton_iff] at he he'
  rcases he with rfl | rfl | rfl <;> rcases he' with rfl | rfl | rfl <;>
    first
    | (refine ⟨t.1, ?_, ?_⟩ <;> simp; done)
    | (refine ⟨t.2.1, ?_, ?_⟩ <;> simp; done)
    | (refine ⟨t.2.2, ?_, ?_⟩ <;> simp)

theorem mem_verts_of_mem_pairsOf {t : V × V × V} {e : Sym2 V} {v : V} (he : e ∈ pairsOf t)
    (hv : v ∈ e) : v ∈ verts t := by
  simp only [pairsOf, Set.mem_insert_iff, Set.mem_singleton_iff] at he
  simp only [verts, Set.mem_insert_iff, Set.mem_singleton_iff]
  rcases he with rfl | rfl | rfl <;> simp [Sym2.mem_iff] at hv <;> tauto

/-- Observation 3.2: any two distinct owners of pairs of one triple share a vertex, i.e. the
owner set is a clique of the intersection graph `I(H)`. -/
theorem owners_clique {t : V × V × V} {i j : ι} (hi : i ∈ owners P t) (hj : j ∈ owners P t) :
    i = j ∨ ∃ v, v ∈ verts (P i) ∧ v ∈ verts (P j) := by
  by_cases hij : i = j
  · exact Or.inl hij
  right
  obtain ⟨e, he1, he2⟩ := Set.not_disjoint_iff.1 hi
  obtain ⟨e', he1', he2'⟩ := Set.not_disjoint_iff.1 hj
  obtain ⟨v, hv, hv'⟩ := pairs_share he2 he2'
  exact ⟨v, mem_verts_of_mem_pairsOf he1 hv, mem_verts_of_mem_pairsOf he1' hv'⟩

end Obs

/-! ### Five states -/

/-- Number of edges on `{a, b, c}` after editing the pairs `x ⊆ {ab, bc, ac}` of an induced
`P₃` `(a, b, c)` (encoded `0 ↦ ab`, `1 ↦ bc`, `2 ↦ ac`; `ab, bc` edges, `ac` non-edge). -/
def edgeCount (x : Finset (Fin 3)) : ℕ :=
  (if (0 : Fin 3) ∈ x then 0 else 1) + (if (1 : Fin 3) ∈ x then 0 else 1) +
    (if (2 : Fin 3) ∈ x then 1 else 0)

/-- The five states: `{ab}`, `{bc}`, `{ac}`, `{ab, bc}`, `{ab, bc, ac}`. -/
def fiveStates : Finset (Finset (Fin 3)) := {{0}, {1}, {2}, {0, 1}, {0, 1, 2}}

/-- The finite case check: a non-empty edit set leaves no triple with exactly two edges iff it
is one of the five states. -/
theorem five_states_fin :
    ∀ x : Finset (Fin 3), (x.Nonempty ∧ edgeCount x ≠ 2) ↔ x ∈ fiveStates := by decide

/-- Exactly the two sets `{ab, ac}` and `{bc, ac}` are non-empty and leave an induced `P₃`. -/
theorem two_bad_states :
    ∀ x : Finset (Fin 3), (x.Nonempty ∧ x ∉ fiveStates) ↔ (x = {0, 2} ∨ x = {1, 2}) := by
  decide

/-- Costs of the five states are `1, 1, 1, 2, 3`. -/
theorem five_states_cost : fiveStates.val.map Finset.card = {1, 1, 1, 2, 3} := by decide

theorem edgeCount_ne_two_nonempty : ∀ x : Finset (Fin 3), edgeCount x ≠ 2 → x.Nonempty := by
  decide

/-- `(a, b, c)` has exactly two edges, i.e. one of the three centres gives an induced `P₃`. -/
def TripleBad (G : SimpleGraph V) (a b c : V) : Prop :=
  IndP3 G a b c ∨ IndP3 G b a c ∨ IndP3 G a c b

open scoped Classical in
/-- A triple of distinct vertices is an induced `P₃` (for some centre) iff it has exactly two
edges. -/
theorem tripleBad_iff (G : SimpleGraph V) {a b c : V} (hab : a ≠ b) (hbc : b ≠ c)
    (hac : a ≠ c) :
    TripleBad G a b c ↔
      (if G.Adj a b then 1 else 0) + (if G.Adj b c then 1 else 0) +
        (if G.Adj a c then 1 else 0) = 2 := by
  unfold TripleBad IndP3
  simp only [G.adj_comm b a, G.adj_comm c b]
  by_cases h1 : G.Adj a b <;> by_cases h2 : G.Adj b c <;> by_cases h3 : G.Adj a c <;>
    simp [h1, h2, h3, hab, hbc, hac]

open scoped Classical in
/-- The set of pairs of the triple `t` that `E` flips. -/
noncomputable def pattern (E : SimpleGraph V) (t : V × V × V) : Finset (Fin 3) :=
  Finset.univ.filter fun k => pairOf t k ∈ E.edgeSet

theorem mem_pattern (E : SimpleGraph V) (t : V × V × V) (k : Fin 3) :
    k ∈ pattern E t ↔ pairOf t k ∈ E.edgeSet := by
  simp [pattern]

theorem pattern_zero (E : SimpleGraph V) (a b c : V) : (0 : Fin 3) ∈ pattern E (a, b, c) ↔ E.Adj a b := by
  rw [mem_pattern]; simp [pairOf]

theorem pattern_one (E : SimpleGraph V) (a b c : V) : (1 : Fin 3) ∈ pattern E (a, b, c) ↔ E.Adj b c := by
  rw [mem_pattern]; simp [pairOf]

theorem pattern_two (E : SimpleGraph V) (a b c : V) : (2 : Fin 3) ∈ pattern E (a, b, c) ↔ E.Adj a c := by
  rw [mem_pattern]; simp [pairOf]

/-- Lemma (Five states), graph version: for an induced `P₃` `(a, b, c)` of `G` and any edit
graph `E`, the triple has no induced `P₃` in `G △ E` iff the set of flipped pairs of the
triple is one of the five states. -/
theorem five_states {G E : SimpleGraph V} {a b c : V} (h : IndP3 G a b c) :
    ¬ TripleBad (edit G E) a b c ↔ pattern E (a, b, c) ∈ fiveStates := by
  rw [← five_states_fin]
  obtain ⟨hab, hbc, hac, hn⟩ := h
  have key : TripleBad (edit G E) a b c ↔ edgeCount (pattern E (a, b, c)) = 2 := by
    rw [tripleBad_iff _ hab.ne hbc.ne hac]
    by_cases h1 : E.Adj a b <;> by_cases h2 : E.Adj b c <;> by_cases h3 : E.Adj a c <;>
      simp [edgeCount, pattern_zero, pattern_one, pattern_two, hab, hbc, hn, h1, h2, h3]
  rw [key]
  constructor
  · intro h2; exact ⟨edgeCount_ne_two_nonempty _ h2, h2⟩
  · intro h2; exact h2.2

/-- In every solution, each member's edit set is one of the five states. -/
theorem member_state_five {G : SimpleGraph V} {P : ι → V × V × V} (hP : IsPacking G P)
    {E : SimpleGraph V} (hE : IsSolution G E) (i : ι) : pattern E (P i) ∈ fiveStates := by
  have := (five_states (E := E) (hP.p3 i)).1
  exact this (fun hb => by
    rcases hb with hb | hb | hb <;> exact hE _ _ _ hb)

end P3
