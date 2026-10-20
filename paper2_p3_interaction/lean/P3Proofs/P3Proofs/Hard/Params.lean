import Mathlib

/-!
# Parameter arithmetic of the reduction

Paper: `app_hardness.tex`, Lemma "Parameters and validity" (b),(d) and Lemma "No super-vertex is
split" (the arithmetic of `φ(x) = x(p - x - u(v))`); `hardness.tex` (`p ≤ k²+6`, `c ≤ 4k²+26`);
`notes/wave10_hardness_proof.md` Sec. 4, 5 (Lemma 5).

Only pure arithmetic is formalised: the combinatorial meaning of `u(v)` (number of matchings used at a
super-vertex, `≤ 3`, `≤ 2` for `k ≤ 2`) and of `c(H) ≤ 4p+2` is NOT formalised; they are hypotheses.
-/

namespace P3.Hard

/-- the least even `p` with `p ≥ k² + 5`. -/
def pOf (k : ℕ) : ℕ := k ^ 2 + 5 + (k ^ 2 + 5) % 2

theorem pOf_even (k : ℕ) : Even (pOf k) := by
  unfold pOf
  rw [Nat.even_iff]; omega

theorem pOf_ge (k : ℕ) : k ^ 2 + 5 ≤ pOf k := by unfold pOf; omega

theorem pOf_le (k : ℕ) : pOf k ≤ k ^ 2 + 6 := by unfold pOf; omega

/-- `pOf k` is the least even number `≥ k²+5`. -/
theorem pOf_least (k p : ℕ) (hp : Even p) (h : k ^ 2 + 5 ≤ p) : pOf k ≤ p := by
  unfold pOf
  rw [Nat.even_iff] at hp
  omega

/-- connection degree: `c = 4p + 2`. -/
def cOf (p : ℕ) : ℕ := 4 * p + 2

/-- `c ≤ 4k² + 26` for `p = pOf k` (paper: "`c ≤ 4k²+26`"). -/
theorem cOf_le (k : ℕ) : cOf (pOf k) ≤ 4 * k ^ 2 + 26 := by
  have := pOf_le k; unfold cOf; omega

/-- `l = k²`, `p ≥ k²+5` gives `p ≥ l + u_max + 2` with `u_max = 3`. -/
theorem p_ge_l_add (k p : ℕ) (h : k ^ 2 + 5 ≤ p) : k ^ 2 + 3 + 2 ≤ p := by omega

/-- **No-split inequality (the case `x = 1` of Lemma "no super-vertex is split").**
`p - 1 - u > l = k²` for `u ≤ 3` and `p ≥ k²+5` (stated without subtraction). -/
theorem no_split_x1 (k p u : ℕ) (hp : k ^ 2 + 5 ≤ p) (hu : u ≤ 3) : k ^ 2 + 1 + u < p := by
  omega

/-- The same for `p = pOf k`. -/
theorem no_split_x1_pOf (k u : ℕ) (hu : u ≤ 3) : k ^ 2 < pOf k - 1 - u := by
  have := pOf_ge k; omega

/-- **Full arithmetic of Lemma "no super-vertex is split".**  Let `p` be even, `p ≥ k²+5`,
`u ≤ 3` (and `u ≤ 2` if `k ≤ 2`).  A part `X` of size `x` (`1 ≤ x ≤ p/2`) of a super-vertex has
at least `x * y` hard edges to the rest, `y = p - x - u`; then `x * y > k²`,
so more than `l = k²` hard edits would be needed. -/
theorem no_split_all (k p u x y : ℕ) (hp : Even p) (hpk : k ^ 2 + 5 ≤ p) (hu : u ≤ 3)
    (hu2 : k ≤ 2 → u ≤ 2) (hxy : x + y + u = p) (hx1 : 1 ≤ x) (hx2 : 2 * x ≤ p) :
    k ^ 2 < x * y := by
  rw [Nat.even_iff] at hp
  -- p/2 ≥ 1 + u
  have hhalf : 1 + u ≤ p / 2 := by
    by_cases hk : k ≤ 2
    · have := hu2 hk
      have h0 := Nat.zero_le (k ^ 2)
      omega
    · have hk3 : 3 ≤ k := by omega
      have : 9 ≤ k ^ 2 := by nlinarith
      omega
  have hy1 : 1 ≤ y := by omega
  -- (x-1)(y-1) ≥ 0 gives x*y ≥ x + y - 1 = p - u - 1 (concavity: the minimum is at the endpoints)
  have h1 : x + y ≤ x * y + 1 := by
    obtain ⟨a, rfl⟩ : ∃ a, x = a + 1 := ⟨x - 1, by omega⟩
    obtain ⟨b, rfl⟩ : ∃ b, y = b + 1 := ⟨y - 1, by omega⟩
    nlinarith [Nat.zero_le (a * b)]
  have : k ^ 2 + 1 + u < p := no_split_x1 k p u hpk hu
  omega

/-! ### Small facts used for `u(v) ≤ 3` -/

/-- the count at a host: permanent guest (1) + lock guest (1) in a tile group, or
permanent guest (1) + at most two comparator guests (2) in a separator group: at most `3`;
at `k ≤ 2` a chain has at most one horizontal neighbour so at most `1 + 1 = 2`. -/
theorem host_count_le (z y cmp : ℕ) (hz : z = 1) (hcmp : cmp ≤ 2) (hy : y ≤ 1) (hex : y = 0 ∨ cmp = 0) :
    z + y + cmp ≤ 3 := by omega

end P3.Hard
