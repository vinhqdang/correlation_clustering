import Mathlib

/-!
# The counting core of soundness (Theorem S) as an abstract inequality system

Paper: `app_hardness.tex`, Lemma "Region bound" and Theorem "Soundness" (proof, the displays
with `r_c = h_c + E_c ≥ max{1,h_c}`), including the corrections F2/F3 of
`notes/audit15_hardness_second.md`.

What is formalised: ONLY the arithmetic.  The graph-level facts (decomposition lemma, Lemma 6
`r_c ≥ 1`, `Σ h_c = h_in + 2q`, the formula for `ex(S)`) are the *hypotheses* of the abstract system
below; they are NOT derived here.

Data (all natural numbers), `C = 2k²` chains indexed by `Fin (2*k^2)`:
* `h c`   : number of edits in `S ∩ Π_c` on pairs of no member of `H_c`,
* `E c`   : member-excess of `H_c` (edits beyond one per member),  `r_c = h_c + E_c`,
* `q`     : number of edited cross-anchor pairs,
* `hin`   : number of hard edits lying in exactly one `Π_c`,
* `sl`    : `|S_rest| - |H_cmp|` (slack of the comparator members),
* `ex`    : the excess `|S| - |H|`.
Hypotheses:
* (hsum)  `Σ_c h_c = h_in + 2q`        (a cross-anchor edit lies in two regions),
* (hreg)  `1 ≤ r_c = h_c + E_c`        (Lemma "Region bound": `r_c ≥ max{1,h_c}`),
* (hex)   `ex + q = Σ_c r_c + sl`      (i.e. `ex = Σ r_c - q + (|S_rest|-|H_cmp|)`).
-/

namespace P3.Hard

open Finset

/-- The abstract inequality system of the proof of Theorem S (arithmetic part only). -/
structure CountSystem (k : ℕ) where
  h : Fin (2 * k ^ 2) → ℕ
  E : Fin (2 * k ^ 2) → ℕ
  q : ℕ
  hin : ℕ
  sl : ℕ
  ex : ℕ
  hsum : ∑ c, h c = hin + 2 * q
  hreg : ∀ c, 1 ≤ h c + E c
  hex : ex + q = ∑ c, (h c + E c) + sl

namespace CountSystem

variable {k : ℕ} (D : CountSystem k)

/-- Closed form of the excess: `ex = h_in + q + ΣE + sl`. -/
theorem ex_eq : D.ex = D.hin + D.q + ∑ c, D.E c + D.sl := by
  have h1 := D.hex
  rw [Finset.sum_add_distrib, D.hsum] at h1
  omega

/-- every edited cross-anchor pair is a hard edit: `q ≤ ex` (follows from the system). -/
theorem q_le_ex : D.q ≤ D.ex := by
  have := D.ex_eq; omega

/-- `Σ_c r_c ≥ 2k²` since every region needs `r_c ≥ 1`. -/
theorem sum_r_ge : 2 * k ^ 2 ≤ ∑ c, (D.h c + D.E c) := by
  calc 2 * k ^ 2 = ∑ _c : Fin (2 * k ^ 2), 1 := by simp
    _ ≤ ∑ c, (D.h c + D.E c) := Finset.sum_le_sum (fun c _ => D.hreg c)

/-- **`ex ≥ 2k² - q`** (stated additively). -/
theorem two_k_sq_le_ex_add_q : 2 * k ^ 2 ≤ D.ex + D.q := by
  have := D.hex; have := D.sum_r_ge; omega

/-- **Theorem S, counting part.** `ex ≥ max(q, 2k² - q) ≥ k²`. -/
theorem ex_ge_max : max D.q (2 * k ^ 2 - D.q) ≤ D.ex := by
  have h1 := D.q_le_ex
  have h2 := D.two_k_sq_le_ex_add_q
  exact max_le h1 (by omega)

theorem k_sq_le_max (q : ℕ) : k ^ 2 ≤ max q (2 * k ^ 2 - q) := by
  by_cases h : k ^ 2 ≤ q
  · exact le_trans h (le_max_left _ _)
  · exact le_trans (by omega) (le_max_right _ _)

theorem k_sq_le_ex : k ^ 2 ≤ D.ex := le_trans (k_sq_le_max D.q) D.ex_ge_max

/-- **Tightness analysis (the exact form of Theorem S).**  If `ex = k²` then
`q = k²`, `h_in = 0`, `sl = 0` (so `|S_rest| = |H_cmp|`), every `E_c = 0` (one edit per member of
every `H_c`), and every chain has exactly one non-member edit: `h_c = 1`. -/
theorem tight (hx : D.ex = k ^ 2) :
    D.q = k ^ 2 ∧ D.hin = 0 ∧ D.sl = 0 ∧ (∀ c, D.E c = 0) ∧ (∀ c, D.h c = 1) := by
  have h1 := D.ex_eq
  have h2 := D.two_k_sq_le_ex_add_q
  have hq : D.q = k ^ 2 := by
    have := D.q_le_ex; omega
  have hsumE : ∑ c, D.E c = 0 := by omega
  have hhin : D.hin = 0 := by omega
  have hsl : D.sl = 0 := by omega
  have hE : ∀ c, D.E c = 0 := by
    intro c
    have := (Finset.sum_eq_zero_iff.1 hsumE) c (Finset.mem_univ c)
    exact this
  have hge : ∀ c, 1 ≤ D.h c := by
    intro c; have := D.hreg c; rw [hE c] at this; omega
  have hh : ∀ c, D.h c = 1 := by
    by_contra hcon
    push_neg at hcon
    obtain ⟨c0, hc0⟩ := hcon
    have hgt : 1 < D.h c0 := lt_of_le_of_ne (hge c0) (Ne.symm hc0)
    have : ∑ _c : Fin (2 * k ^ 2), 1 < ∑ c, D.h c :=
      Finset.sum_lt_sum (fun c _ => hge c) ⟨c0, Finset.mem_univ _, hgt⟩
    have h3 := D.hsum
    simp at this
    omega
  exact ⟨hq, hhin, hsl, hE, hh⟩

end CountSystem

/-- Converse (the bound is attained by a consistent system): `h_c = 1`, `E_c = 0`, `q = k²`,
`h_in = 0`, `sl = 0` gives `ex = k²`. -/
def tightSystem (k : ℕ) : CountSystem k where
  h := fun _ => 1
  E := fun _ => 0
  q := k ^ 2
  hin := 0
  sl := 0
  ex := k ^ 2
  hsum := by simp; ring
  hreg := by intro c; simp
  hex := by simp; ring

theorem tightSystem_ex (k : ℕ) : (tightSystem k).ex = k ^ 2 := rfl

/-! ### The paper's route through `Z = {c : h_c = 0}` -/

namespace CountSystem

variable {k : ℕ} (D : CountSystem k)

/-- `Z = {c : h_c = 0}`. -/
def Z : Finset (Fin (2 * k ^ 2)) := Finset.univ.filter (fun c => D.h c = 0)

/-- `r_c ≥ max{1, h_c}`, written as `r_c ≥ h_c + [h_c = 0]`. -/
theorem r_ge_max (c : Fin (2 * k ^ 2)) : D.h c + (if D.h c = 0 then 1 else 0) ≤ D.h c + D.E c := by
  have := D.hreg c
  split_ifs with h <;> omega

/-- `ex ≥ |Z| + h_in + q` (the displayed inequality of the paper's proof, the correct version
without the term `ΣE_c`). -/
theorem ex_ge_Z : D.Z.card + D.hin + D.q ≤ D.ex := by
  have h1 : ∑ c, (D.h c + (if D.h c = 0 then 1 else 0)) ≤ ∑ c, (D.h c + D.E c) :=
    Finset.sum_le_sum (fun c _ => D.r_ge_max c)
  rw [Finset.sum_add_distrib, D.hsum] at h1
  have h2 : ∑ c, (if D.h c = 0 then 1 else 0) = D.Z.card := by
    simp [Z, Finset.card_filter]
  rw [h2] at h1
  have := D.hex
  omega

/-- `|Z| ≥ 2k² - Σ h_c` (the chains outside `Z` number at most `Σ h_c`). -/
theorem Z_ge : 2 * k ^ 2 ≤ D.Z.card + ∑ c, D.h c := by
  have h1 : ∑ _c : Fin (2 * k ^ 2), 1 ≤ ∑ c, (D.h c + (if D.h c = 0 then 1 else 0)) :=
    Finset.sum_le_sum (fun c _ => by split_ifs with h <;> omega)
  rw [Finset.sum_add_distrib] at h1
  have h2 : ∑ c, (if D.h c = 0 then 1 else 0) = D.Z.card := by
    simp [Z, Finset.card_filter]
  rw [h2] at h1
  simp at h1
  omega

end CountSystem

/-- **The printed inequality `ex ≥ |Z| + h_in + q + Σ_c E_c` is not valid in general**
(it counts `E_c` twice for the chains of `Z`, for which `r_c = E_c`).  The abstract system
admits `k = 1`, `h = (0,0)`, `E = (1,1)`, `q = h_in = sl = 0`, `ex = 2`, where the right-hand
side is `4`.  The conclusions of Theorem S are unaffected (see `CountSystem.tight`, which uses
only valid steps; at `ex = k²` one has `Z = ∅`). -/
def counterexampleSystem : CountSystem 1 where
  h := fun _ => 0
  E := fun _ => 1
  q := 0
  hin := 0
  sl := 0
  ex := 2
  hsum := by simp
  hreg := by intro c; simp
  hex := by simp

theorem printed_inequality_fails :
    ¬ ((counterexampleSystem.Z.card + counterexampleSystem.hin + counterexampleSystem.q
        + ∑ c, counterexampleSystem.E c) ≤ counterexampleSystem.ex) := by
  simp [counterexampleSystem, CountSystem.Z]

/-- The correct version: only chains outside `Z` contribute `E_c`
(for `c ∈ Z` the slack is `E_c - 1`). -/
theorem CountSystem.ex_ge_Z_E {k : ℕ} (D : CountSystem k) :
    D.Z.card + D.hin + D.q + ∑ c ∈ (D.Z)ᶜ, D.E c ≤ D.ex := by
  have e := D.ex_eq
  have hZ : D.Z.card ≤ ∑ c ∈ D.Z, D.E c := by
    calc D.Z.card = ∑ _c ∈ D.Z, 1 := by simp
      _ ≤ ∑ c ∈ D.Z, D.E c := by
        apply Finset.sum_le_sum
        intro c hc
        have hc0 : D.h c = 0 := by simpa [CountSystem.Z] using hc
        have := D.hreg c; omega
  have hs : ∑ c, D.E c = ∑ c ∈ D.Z, D.E c + ∑ c ∈ (D.Z)ᶜ, D.E c := by
    rw [Finset.sum_add_sum_compl]
  omega

end P3.Hard
