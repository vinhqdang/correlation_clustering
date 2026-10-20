import P3Proofs.Struct.Max

/-!
# HUB⁺ (paper: `lem:HUB`, `cor:HUB`, app_counting.tex)

`u` a vertex, `N_h(u)` the hard neighbours, `Λ_u` the graph on `N_h(u)` whose edges are the
hard non-edges `{x, y}` (`¬ G.Adj x y`, `x ≠ y`, `¬ Soft H x y`).  A matching of `Λ` inside
`J` of size `t` is expressed by the recursive predicate `HasMatch J t`.

* `hub_matching` (a): in a YES instance `Sol G α (|H| + ℓ)`, every matching of `Λ_u` has
  size at most `ℓ` (`ν(Λ_u) ≤ ℓ`).
* `hub_cover` (b): there is `T ⊆ J` with `|T| ≤ 2ℓ` such that every pair in `J \ T` is an edge
  of `G` or soft (`T` is a vertex cover of `Λ`).
* `hub_alpha` (c): `J ⊆ N_h(u)` independent in `G` has `|J| ≤ 2ℓ + c + 1`, where `c` bounds
  the number of members of `H` through a vertex.
* `hub_clique` (Cor. HUB): if every pair of `Y ⊆ N_h(u)` is a hard non-edge, `|Y| ≤ 2ℓ + 1`.
-/

open Finset

namespace P3S

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]
variable {α : Type*} [DecidableEq α]

/-- `HasMatch H J t`: `Λ[J]` (hard non-edges inside `J`) has a matching with `t` edges. -/
def HasMatch (H : Finset (V × V × V)) : Finset V → ℕ → Prop
  | _, 0 => True
  | J, t + 1 => ∃ x ∈ J, ∃ y ∈ J, x ≠ y ∧ ¬ G.Adj x y ∧ ¬ Soft H x y ∧
      HasMatch H ((J.erase x).erase y) t

lemma hasMatch_mono {H : Finset (V × V × V)} :
    ∀ (t : ℕ) {J J' : Finset V}, J ⊆ J' → HasMatch G H J t → HasMatch G H J' t
  | 0, _, _, _, _ => trivial
  | t + 1, J, J', hJ, ⟨x, hx, y, hy, hxy, hn, hs, h⟩ =>
      ⟨x, hJ hx, y, hJ hy, hxy, hn, hs,
        hasMatch_mono G t (fun z hz => by
          simp only [Finset.mem_erase] at hz ⊢
          exact ⟨hz.1, hz.2.1, hJ hz.2.2⟩) h⟩

/-- A matching of `Λ_u` yields a family of `P₃`s `(x, u, y)` hard over `H`. -/
lemma hasMatch_hardFam {H : Finset (V × V × V)} (u : V) :
    ∀ (t : ℕ) (J : Finset V), (∀ x ∈ J, G.Adj u x ∧ ¬ Soft H u x) → HasMatch G H J t →
      ∃ M : Finset (V × V × V), HardFam G H M ∧ M.card = t ∧
        ∀ Q ∈ M, Q.1 ∈ J ∧ Q.2.1 = u ∧ Q.2.2 ∈ J
  | 0, J, _, _ => ⟨∅, ⟨by simp, by simp, by simp⟩, by simp, by simp⟩
  | t + 1, J, hJ, ⟨x, hx, y, hy, hxy, hn, hs, hrest⟩ => by
    set J' := (J.erase x).erase y with hJ'
    have hJ'sub : ∀ z ∈ J', z ∈ J ∧ z ≠ x ∧ z ≠ y := by
      intro z hz
      simp only [hJ', Finset.mem_erase] at hz
      exact ⟨hz.2.2, hz.2.1, hz.1⟩
    obtain ⟨M', hM', hcard, hmem⟩ :=
      hasMatch_hardFam u t J' (fun z hz => hJ z (hJ'sub z hz).1) hrest
    have hux : G.Adj u x := (hJ x hx).1
    have huy : G.Adj u y := (hJ y hy).1
    have hsx : ¬ Soft H u x := (hJ x hx).2
    have hsy : ¬ Soft H u y := (hJ y hy).2
    have hsx' : ¬ Soft H x u := fun h => hsx h.symm
    have hsy' : ¬ Soft H y u := fun h => hsy h.symm
    have hs' : ¬ Soft H y x := fun h => hs h.symm
    have hQnot : (x, u, y) ∉ M' := by
      intro hm
      have := (hmem _ hm).1
      have h2 := hJ'sub _ this
      exact h2.2.1 rfl
    refine ⟨insert (x, u, y) M', ⟨?_, ?_, ?_⟩, ?_, ?_⟩
    · intro Q hQ
      rcases Finset.mem_insert.mp hQ with rfl | hQ
      · exact ⟨hux.symm, huy, hn, hxy⟩
      · exact hM'.p3 Q hQ
    · rw [Finset.coe_insert]
      refine hM'.disj.insert ?_
      intro j hj _
      obtain ⟨hj1, hj2, hj3⟩ := hmem j hj
      obtain ⟨h1, h2, h3⟩ := hJ'sub _ hj1
      obtain ⟨h4, h5, h6⟩ := hJ'sub _ hj3
      have hjx := (hJ _ h1).1
      have hjy := (hJ _ h4).1
      rw [hj2] at hjx hjy
      rw [Finset.disjoint_left]
      rintro ⟨p, q⟩ hp hq
      simp only [mem_pr, Prod.mk.injEq] at hp hq
      have n1 : u ≠ x := hux.ne
      have n2 : u ≠ y := huy.ne
      have n3 : u ≠ j.1 := hjx.ne
      have n4 : u ≠ j.2.2 := hjy.ne
      rw [hj2] at hq
      rcases hp with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ |
        ⟨rfl, rfl⟩ <;>
      rcases hq with ⟨h7, h8⟩ | ⟨h7, h8⟩ | ⟨h7, h8⟩ | ⟨h7, h8⟩ | ⟨h7, h8⟩ | ⟨h7, h8⟩ <;>
      simp_all [eq_comm]
    · intro Q hQ x' hx'
      rcases Finset.mem_insert.mp hQ with rfl | hQ
      · simp only [mem_pr, Prod.mk.injEq] at hx'
        rcases x' with ⟨p, q⟩
        simp only at hx' ⊢
        rcases hx' with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ |
          ⟨rfl, rfl⟩ <;> assumption
      · exact hM'.hard Q hQ x' hx'
    · rw [Finset.card_insert_of_notMem hQnot, hcard]
    · intro Q hQ
      rcases Finset.mem_insert.mp hQ with rfl | hQ
      · exact ⟨hx, rfl, hy⟩
      · obtain ⟨a, b, c⟩ := hmem Q hQ
        exact ⟨(hJ'sub _ a).1, b, (hJ'sub _ c).1⟩

/-- **HUB (a).** In a YES instance every matching of `Λ_u` has at most `ℓ` edges. -/
theorem hub_matching {H : Finset (V × V × V)} (hH : IsPacking G H) (ℓ : ℕ)
    (hsol : Sol G α ((H.card : ℤ) + ℓ)) (u : V) (J : Finset V)
    (hJ : ∀ x ∈ J, G.Adj u x ∧ ¬ Soft H u x) (t : ℕ) (hm : HasMatch G H J t) : t ≤ ℓ := by
  obtain ⟨M, hM, hcard, -⟩ := hasMatch_hardFam G u t J hJ hm
  have := hardFam_card_le G hH hM ℓ hsol
  omega

/-- **HUB (b).** A matching bound `k` gives a vertex cover `T` of `Λ[J]` with `|T| ≤ 2k`. -/
theorem hub_cover {H : Finset (V × V × V)} :
    ∀ (k : ℕ) (J : Finset V), (∀ t, HasMatch G H J t → t ≤ k) →
      ∃ T ⊆ J, T.card ≤ 2 * k ∧
        ∀ x ∈ J \ T, ∀ y ∈ J \ T, x ≠ y → G.Adj x y ∨ Soft H x y
  | k, J, hk => by
    by_cases hex : ∃ x ∈ J, ∃ y ∈ J, x ≠ y ∧ ¬ G.Adj x y ∧ ¬ Soft H x y
    · obtain ⟨x, hx, y, hy, hxy, hn, hs⟩ := hex
      cases k with
      | zero =>
        have := hk 1 ⟨x, hx, y, hy, hxy, hn, hs, trivial⟩
        omega
      | succ k =>
        have hrec := hub_cover k ((J.erase x).erase y) (by
          intro t ht
          have := hk (t + 1) ⟨x, hx, y, hy, hxy, hn, hs, ht⟩
          omega)
        obtain ⟨T', hT'sub, hT'card, hT'⟩ := hrec
        refine ⟨insert x (insert y T'), ?_, ?_, ?_⟩
        · intro z hz
          simp only [Finset.mem_insert] at hz
          rcases hz with rfl | rfl | hz
          · exact hx
          · exact hy
          · have := hT'sub hz
            simp only [Finset.mem_erase] at this
            exact this.2.2
        · calc (insert x (insert y T')).card ≤ (insert y T').card + 1 := Finset.card_insert_le _ _
            _ ≤ T'.card + 1 + 1 := by
              have := Finset.card_insert_le y T'
              omega
            _ ≤ 2 * (k + 1) := by omega
        · intro a ha b hb hab
          simp only [Finset.mem_sdiff, Finset.mem_insert, not_or] at ha hb
          apply hT' a _ b _ hab
          · simp only [Finset.mem_sdiff, Finset.mem_erase]
            exact ⟨⟨ha.2.1, ha.2.2.1, ha.1⟩, ha.2.2.2⟩
          · simp only [Finset.mem_sdiff, Finset.mem_erase]
            exact ⟨⟨hb.2.1, hb.2.2.1, hb.1⟩, hb.2.2.2⟩
    · refine ⟨∅, by simp, by simp, ?_⟩
      intro x hx y hy hxy
      simp only [Finset.sdiff_empty] at hx hy
      by_contra hno
      push Not at hno
      exact hex ⟨x, hx, y, hy, hxy, hno.1, hno.2⟩

/-- Soft non-edges: the two ends of a soft non-adjacent pair are the two ends of a member. -/
lemma soft_nonadj {H : Finset (V × V × V)} (hH : IsPacking G H) {x y : V} (hs : Soft H x y)
    (hn : ¬ G.Adj x y) : ∃ P ∈ H, (x = P.1 ∧ y = P.2.2) ∨ (x = P.2.2 ∧ y = P.1) := by
  obtain ⟨P, hP, hx⟩ := hs
  refine ⟨P, hP, ?_⟩
  obtain ⟨h1, h2, h3, h4⟩ := hH.p3 P hP
  simp only [mem_pr, Prod.mk.injEq] at hx
  rcases hx with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
  · exact absurd h1 hn
  · exact absurd h1.symm hn
  · exact absurd h2 hn
  · exact absurd h2.symm hn
  · exact Or.inl ⟨rfl, rfl⟩
  · exact Or.inr ⟨rfl, rfl⟩

/-- **HUB (c).** Let `J ⊆ N_h(u)` be independent in `G` and let at most `c` members of `H`
pass through any vertex.  Then `|J| ≤ 2ℓ + c + 1`. -/
theorem hub_alpha {H : Finset (V × V × V)} (hH : IsPacking G H) (ℓ : ℕ)
    (hsol : Sol G α ((H.card : ℤ) + ℓ)) (c : ℕ)
    (hc : ∀ v, (H.filter fun P => v = P.1 ∨ v = P.2.1 ∨ v = P.2.2).card ≤ c)
    (u : V) (J : Finset V) (hJ : ∀ x ∈ J, G.Adj u x ∧ ¬ Soft H u x)
    (hind : ∀ x ∈ J, ∀ y ∈ J, ¬ G.Adj x y) : J.card ≤ 2 * ℓ + c + 1 := by
  obtain ⟨T, hTsub, hTcard, hT⟩ := hub_cover G ℓ J
    (fun t ht => hub_matching G hH ℓ hsol u J hJ t ht)
  set J' := J \ T with hJ'
  have hsoft : ∀ x ∈ J', ∀ y ∈ J', x ≠ y → Soft H x y := fun x hx y hy hxy =>
    (hT x hx y hy hxy).resolve_left (hind x (Finset.mem_sdiff.mp hx).1 y (Finset.mem_sdiff.mp hy).1)
  have hJ'card : J'.card ≤ c + 1 := by
    rcases J'.eq_empty_or_nonempty with he | ⟨x, hx⟩
    · rw [he]; simp
    · have hmap : ∀ y ∈ J'.erase x, ∃ P ∈ H.filter (fun P => x = P.1 ∨ x = P.2.1 ∨ x = P.2.2),
          (x = P.1 ∧ y = P.2.2) ∨ (x = P.2.2 ∧ y = P.1) := by
        intro y hy
        rw [Finset.mem_erase] at hy
        obtain ⟨P, hP, hP'⟩ := soft_nonadj G hH (hsoft x hx y hy.2 (Ne.symm hy.1))
          (hind x (Finset.mem_sdiff.mp hx).1 y (Finset.mem_sdiff.mp hy.2).1)
        refine ⟨P, ?_, hP'⟩
        rw [Finset.mem_filter]
        refine ⟨hP, ?_⟩
        rcases hP' with ⟨h, _⟩ | ⟨h, _⟩
        · exact Or.inl h
        · exact Or.inr (Or.inr h)
      choose! f hf1 hf2 using hmap
      have hinj : Set.InjOn f (J'.erase x : Set V) := by
        intro y hy y' hy' hyy
        have a := hf2 y hy
        have b := hf2 y' hy'
        rw [hyy] at a
        obtain ⟨h1, h2, h3, h4⟩ := hH.p3 (f y') (Finset.mem_filter.mp (hf1 y' hy')).1
        rcases a with ⟨a1, a2⟩ | ⟨a1, a2⟩ <;> rcases b with ⟨b1, b2⟩ | ⟨b1, b2⟩
        · rw [a2, b2]
        · exact absurd (a1.symm.trans b1) h4
        · exact absurd (a1.symm.trans b1) h4.symm
        · rw [a2, b2]
      have := Finset.card_le_card_of_injOn f (s := J'.erase x)
        (t := H.filter fun P => x = P.1 ∨ x = P.2.1 ∨ x = P.2.2) (fun y hy => hf1 y hy) hinj
      have h2 := hc x
      have h3 := Finset.card_erase_add_one hx
      omega
  have := Finset.card_le_card_sdiff_add_card (s := J) (t := T)
  omega

/-- **Corollary HUB.** If every pair of `Y ⊆ N_h(u)` is a hard non-edge then `|Y| ≤ 2ℓ + 1`. -/
theorem hub_clique {H : Finset (V × V × V)} (hH : IsPacking G H) (ℓ : ℕ)
    (hsol : Sol G α ((H.card : ℤ) + ℓ)) (u : V) (Y : Finset V)
    (hY : ∀ x ∈ Y, G.Adj u x ∧ ¬ Soft H u x)
    (hpair : ∀ x ∈ Y, ∀ y ∈ Y, x ≠ y → ¬ G.Adj x y ∧ ¬ Soft H x y) : Y.card ≤ 2 * ℓ + 1 := by
  -- a clique of `Λ` of size `2m` contains a matching of size `m`
  have key : ∀ m : ℕ, ∀ Y' : Finset V, Y' ⊆ Y → 2 * m ≤ Y'.card → HasMatch G H Y' m := by
    intro m
    induction m with
    | zero => intro Y' _ _; trivial
    | succ m ih =>
      intro Y' hY' hcard
      have h1 : 0 < Y'.card := by omega
      obtain ⟨x, hx⟩ := Finset.card_pos.mp h1
      have h2 : 0 < (Y'.erase x).card := by rw [Finset.card_erase_of_mem hx]; omega
      obtain ⟨y, hy⟩ := Finset.card_pos.mp h2
      have hyx : y ≠ x := (Finset.mem_erase.mp hy).1
      have hy' : y ∈ Y' := (Finset.mem_erase.mp hy).2
      refine ⟨x, hx, y, hy', Ne.symm hyx, (hpair x (hY' hx) y (hY' hy') (Ne.symm hyx)).1,
        (hpair x (hY' hx) y (hY' hy') (Ne.symm hyx)).2, ?_⟩
      apply ih
      · intro z hz; simp only [Finset.mem_erase] at hz; exact hY' hz.2.2
      · rw [Finset.card_erase_of_mem hy, Finset.card_erase_of_mem hx]; omega
  by_contra hno
  push Not at hno
  have := hub_matching G hH ℓ hsol u Y hY (ℓ + 1) (key (ℓ + 1) Y (Finset.Subset.refl _) (by omega))
  omega

end P3S
