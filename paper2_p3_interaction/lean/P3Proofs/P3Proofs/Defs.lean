import Mathlib

/-!
# Definitions: graphs, edits, cluster graphs, cost, packings

* A graph is a `SimpleGraph V`.  An *edit set* is itself a `SimpleGraph E` on the same
  vertex type (its edges are the flipped pairs); `edit G E` is `G △ E`, and the size of an
  edit set is `size E = E.edgeSet.ncard`.
* `IndP3 G a b c`: `(a, b, c)` is an induced `P₃` with centre `b`.
* `IsCluster G`: no induced `P₃` (paper: "no vertex triple with exactly two edges", see
  `tripleBad_iff` in `CSP.lean`).
* `clusterGraph c`: the cluster graph of a clustering `c : V → α` (cliques = label classes).
* `cost G c`: number of pairs on which `G` and `clusterGraph c` differ
  (= size of the edit set `G △ clusterGraph c`).
* `OPT G`: minimum size of a solution (edit set `E` with `G △ E` a cluster graph).
* Packings: `P : ι → V × V × V`, every member an induced `P₃`, pair sets pairwise disjoint
  (`IsPacking`); `IsPairMaximal`: no induced `P₃` of `G` has all three pairs hard.
-/

open Classical

namespace P3

universe u
variable {V : Type u}

/-- `G △ E`: flip the adjacency of every pair that is an edge of `E`. -/
def edit (G E : SimpleGraph V) : SimpleGraph V where
  Adj u v := ¬ (G.Adj u v ↔ E.Adj u v)
  symm := ⟨fun u v h => by simpa only [G.adj_comm u v, E.adj_comm u v] using h⟩
  loopless := ⟨fun u h => h (by simp)⟩

@[simp] theorem edit_adj (G E : SimpleGraph V) (u v : V) :
    (edit G E).Adj u v ↔ ¬ (G.Adj u v ↔ E.Adj u v) := Iff.rfl

theorem edit_edit (G E : SimpleGraph V) : edit G (edit G E) = E := by
  ext u v
  simp only [edit_adj]
  tauto

/-- Size of an edit set: the number of flipped pairs. -/
noncomputable def size (E : SimpleGraph V) : ℕ := E.edgeSet.ncard

/-- `(a, b, c)` is an induced `P₃` with centre `b`: `ab, bc ∈ E`, `ac ∉ E`. -/
def IndP3 (G : SimpleGraph V) (a b c : V) : Prop :=
  G.Adj a b ∧ G.Adj b c ∧ a ≠ c ∧ ¬ G.Adj a c

/-- A cluster graph is a graph without induced `P₃`. -/
def IsCluster (G : SimpleGraph V) : Prop := ∀ a b c, ¬ IndP3 G a b c

/-- The cluster graph of a labelling: distinct vertices with equal labels are adjacent. -/
def clusterGraph {α : Type*} (c : V → α) : SimpleGraph V where
  Adj u v := u ≠ v ∧ c u = c v
  symm := ⟨fun _ _ h => ⟨h.1.symm, h.2.symm⟩⟩
  loopless := ⟨fun _ h => h.1 rfl⟩

/-- Cost of a clustering: number of pairs `ab` with `[ab ∈ E] ≠ [a, b in one cluster]`. -/
noncomputable def cost {α : Type*} (G : SimpleGraph V) (c : V → α) : ℕ :=
  size (edit G (clusterGraph c))

/-- `E` is a solution: `G △ E` is a cluster graph. -/
def IsSolution (G E : SimpleGraph V) : Prop := IsCluster (edit G E)

/-- `OPT(G)`: minimum size of a solution. -/
noncomputable def OPT (G : SimpleGraph V) : ℕ :=
  sInf {k | ∃ E, IsSolution G E ∧ size E = k}

/-- `OPT` computed over clusterings with labels in `V` (every partition of `V` can be
labelled by vertices). -/
noncomputable def OPTc (G : SimpleGraph V) : ℕ :=
  sInf {k | ∃ c : V → V, cost G c = k}

/-! ### Cluster graphs: transitivity and clusterings -/

theorem isCluster_iff_trans (G : SimpleGraph V) :
    IsCluster G ↔ ∀ a b c, G.Adj a b → G.Adj b c → a ≠ c → G.Adj a c := by
  constructor
  · intro h a b c hab hbc hac
    by_contra hn
    exact h a b c ⟨hab, hbc, hac, hn⟩
  · intro h a b c ⟨hab, hbc, hac, hn⟩
    exact hn (h a b c hab hbc hac)

theorem isCluster_clusterGraph {α : Type*} (c : V → α) : IsCluster (clusterGraph c) := by
  rw [isCluster_iff_trans]
  intro a b d hab hbd had
  exact ⟨had, hab.2.trans hbd.2⟩

/-- A graph is a cluster graph iff it is the cluster graph of some clustering. -/
theorem isCluster_iff_exists_clustering (G : SimpleGraph V) :
    IsCluster G ↔ ∃ c : V → V, G = clusterGraph c := by
  constructor
  · intro h
    rw [isCluster_iff_trans] at h
    let s : Setoid V :=
      { r := fun u v => u = v ∨ G.Adj u v
        iseqv :=
          { refl := fun u => Or.inl rfl
            symm := by
              rintro u v (rfl | h)
              · exact Or.inl rfl
              · exact Or.inr h.symm
            trans := by
              rintro u v w (rfl | h1) (rfl | h2)
              · exact Or.inl rfl
              · exact Or.inr h2
              · exact Or.inr h1
              · by_cases huw : u = w
                · exact Or.inl huw
                · exact Or.inr (h u v w h1 h2 huw) } }
    refine ⟨fun u => (Quotient.mk s u).out, ?_⟩
    ext u v
    have key : (Quotient.mk s u).out = (Quotient.mk s v).out ↔ u = v ∨ G.Adj u v := by
      rw [Quotient.out_inj]
      exact Quotient.eq (r := s)
    simp only [clusterGraph, key]
    constructor
    · intro huv
      exact ⟨huv.ne, Or.inr huv⟩
    · rintro ⟨hne, rfl | hadj⟩
      · exact absurd rfl hne
      · exact hadj
  · rintro ⟨c, rfl⟩
    exact isCluster_clusterGraph c

/-! ### Cost of a clustering = size of an edit set -/

theorem cost_eq_size {α : Type*} (G : SimpleGraph V) (c : V → α) :
    cost G c = size (edit G (clusterGraph c)) := rfl

theorem isSolution_edit_clustering {α : Type*} (G : SimpleGraph V) (c : V → α) :
    IsSolution G (edit G (clusterGraph c)) := by
  unfold IsSolution
  rw [edit_edit]
  exact isCluster_clusterGraph c

theorem exists_clustering_of_solution {G E : SimpleGraph V} (h : IsSolution G E) :
    ∃ c : V → V, cost G c = size E ∧ edit G E = clusterGraph c := by
  obtain ⟨c, hc⟩ := (isCluster_iff_exists_clustering _).1 h
  refine ⟨c, ?_, hc⟩
  rw [cost_eq_size, ← hc, edit_edit]

theorem solution_nonempty (G : SimpleGraph V) :
    {k | ∃ E, IsSolution G E ∧ size E = k}.Nonempty := by
  refine ⟨size G, G, ?_, rfl⟩
  unfold IsSolution
  have : edit G G = ⊥ := by ext u v; simp
  rw [this, isCluster_iff_trans]
  intro a b c h; exact absurd h (by simp)

theorem OPT_le {G E : SimpleGraph V} (h : IsSolution G E) : OPT G ≤ size E :=
  Nat.sInf_le ⟨E, h, rfl⟩

theorem exists_opt (G : SimpleGraph V) : ∃ E, IsSolution G E ∧ size E = OPT G :=
  Nat.sInf_mem (solution_nonempty G)

/-- Optimum over edit sets equals optimum over clusterings. -/
theorem OPT_eq_OPTc (G : SimpleGraph V) : OPT G = OPTc G := by
  apply le_antisymm
  · have hne : {k | ∃ c : V → V, cost G c = k}.Nonempty := ⟨_, fun u => u, rfl⟩
    obtain ⟨c, hc⟩ := Nat.sInf_mem hne
    calc OPT G ≤ size (edit G (clusterGraph c)) := OPT_le (isSolution_edit_clustering G c)
      _ = cost G c := rfl
      _ = OPTc G := hc
  · obtain ⟨E, hE, hk⟩ := exists_opt G
    obtain ⟨c, hc, -⟩ := exists_clustering_of_solution hE
    rw [← hk, ← hc]
    exact Nat.sInf_le ⟨c, rfl⟩

/-! ### Packings -/

/-- The three pairs inside the vertex set of a triple `(a, b, c)`. -/
def pairsOf (t : V × V × V) : Set (Sym2 V) := {s(t.1, t.2.1), s(t.2.1, t.2.2), s(t.1, t.2.2)}

/-- The vertex set of a triple. -/
def verts (t : V × V × V) : Set V := {t.1, t.2.1, t.2.2}

/-- A modification-disjoint packing indexed by `ι`: members are induced `P₃`'s of `G`
with pairwise disjoint pair sets. -/
structure IsPacking {ι : Type*} (G : SimpleGraph V) (P : ι → V × V × V) : Prop where
  p3 : ∀ i, IndP3 G (P i).1 (P i).2.1 (P i).2.2
  disj : Pairwise fun i j => Disjoint (pairsOf (P i)) (pairsOf (P j))

/-- A pair is soft if it lies in a member. -/
def IsSoft {ι : Type*} (P : ι → V × V × V) (e : Sym2 V) : Prop := ∃ i, e ∈ pairsOf (P i)

/-- Pair-maximal: no induced `P₃` of `G` has all three pairs hard. -/
def IsPairMaximal {ι : Type*} (G : SimpleGraph V) (P : ι → V × V × V) : Prop :=
  IsPacking G P ∧ ∀ a b c, IndP3 G a b c → ∃ e ∈ pairsOf (a, b, c), IsSoft P e

/-- A packing is injective on members (follows from disjointness). -/
theorem IsPacking.injective {ι : Type*} {G : SimpleGraph V} {P : ι → V × V × V}
    (hP : IsPacking G P) : Function.Injective P := by
  intro i j h
  by_contra hne
  have hd := hP.disj hne
  have h1 : s((P i).1, (P i).2.1) ∈ pairsOf (P i) := by simp [pairsOf]
  have h2 : s((P i).1, (P i).2.1) ∈ pairsOf (P j) := h ▸ h1
  exact Set.disjoint_left.1 hd h1 h2

end P3
