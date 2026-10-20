import Mathlib

/-!
# The transfer-matrix algebra of the chain lemma, Lemma "Chain lemma" (b)

Paper: `app_hardness.tex`, proof of Lemma "Chain lemma", part (b);
`notes/wave10_hardness_proof.md` Lemma 1(b).

The paper's argument: the number `N_m` of ZE partitions of `OH_m` is `u^T M^{m-c0} v_m` for the
fixed 6x6 integer matrix `M` below; `x^4 | χ_M`, `χ_M = x^4 (x-1)^2`, so `M^t` is an affine function
of `t` for `t ≥ 4`, hence `N_m` is affine in `m` for large `m`; since `N_m = m` for `m = 7..40`,
`N_m = m` for all `m ≥ 7`.

What is formalised here (no `sorry`, no axiom except the core ones):
1. `charpoly Mt = X^4 (X-1)^2` for the explicit matrix (Mathlib `Matrix.charpoly`);
2. `Mt^4 (Mt-1)^2 = 0` (Cayley-Hamilton);
3. `Mt^t = Mt^4 + (t-4) (Mt^5 - Mt^4)` for `t ≥ 4` (affine growth);
4. for ANY vectors `u, v` and shift `c0`, the sequence `N m = u^T Mt^(m-c0) v` is affine for
   `m ≥ c0+4`, and if `N` equals `m` at two consecutive arguments `≥ c0+4`, then `N m = m` for all
   `m ≥ c0 + 4`.

NOT formalised (this is the part of the chain lemma that rests on the computer-generated transfer
matrix): that the cell-by-cell dynamic programme of `OH_m` really has this 6x6 matrix `M` as its
transition matrix at every interior cell of every `m`, with left vector `u` independent of `m` and
right vector `v_m` depending only on the distance to the right end.  The step "`N_m` is affine in `m`
hence equals `m`" is the ONLY place of Lemma 1(b) where the matrix argument is used.
-/

namespace P3.Hard

open Polynomial Matrix

/-- The interior transfer matrix `M` of the note (rows/columns: the six relative states). -/
def Mt : Matrix (Fin 6) (Fin 6) ℤ :=
  !![0,0,0,0,1,1;
     0,0,0,0,0,0;
     1,1,1,1,0,0;
     0,0,0,0,0,0;
     0,0,0,0,1,1;
     0,0,0,0,0,0]

/-- **Characteristic polynomial** of the transfer matrix: `x^4 (x-1)^2`. -/
theorem Mt_charpoly : Mt.charpoly = X ^ 4 * (X - 1) ^ 2 := by
  rw [Matrix.charpoly, Matrix.det_succ_row_zero]
  simp [Fin.sum_univ_succ, Matrix.charmatrix, Mt, Matrix.det_succ_row_zero, Fin.succAbove]
  ring

/-- Cayley-Hamilton: `M^4 (M-1)^2 = 0`. -/
theorem Mt_ch : Mt ^ 4 * (Mt - 1) ^ 2 = 0 := by
  have h := Matrix.aeval_self_charpoly Mt
  rw [Mt_charpoly] at h
  simpa [map_mul, map_pow, map_sub] using h

/-- `M^(t+2) - 2 M^(t+1) + M^t = 0` for `t ≥ 4`. -/
theorem Mt_second_diff (t : ℕ) (ht : 4 ≤ t) : Mt ^ (t + 2) - 2 • Mt ^ (t + 1) + Mt ^ t = 0 := by
  obtain ⟨s, rfl⟩ : ∃ s, t = s + 4 := ⟨t - 4, by omega⟩
  have e : Mt ^ (s + 4 + 2) - 2 • Mt ^ (s + 4 + 1) + Mt ^ (s + 4)
      = Mt ^ s * (Mt ^ 4 * (Mt - 1) ^ 2) := by
    have e1 : Mt ^ (s + 4 + 2) = Mt ^ s * (Mt ^ 4 * Mt ^ 2) := by
      rw [← pow_add, ← pow_add]
    have e2 : Mt ^ (s + 4 + 1) = Mt ^ s * (Mt ^ 4 * Mt) := by
      rw [← pow_succ, ← pow_add]
    have e3 : Mt ^ (s + 4) = Mt ^ s * Mt ^ 4 := by
      rw [← pow_add]
    rw [e1, e2, e3, two_smul]
    noncomm_ring
  rw [e, Mt_ch, mul_zero]

/-- **Affine growth** of the powers: `M^t = M^4 + (t-4)(M^5 - M^4)` for `t ≥ 4`. -/
theorem Mt_pow_affine (t : ℕ) (ht : 4 ≤ t) :
    Mt ^ t = Mt ^ 4 + (t - 4) • (Mt ^ 5 - Mt ^ 4) := by
  induction t, ht using Nat.le_induction with
  | base => simp
  | succ t ht ih =>
    -- M^(t+1) - M^t = M^5 - M^4 (by induction on the second differences)
    have hdiff : ∀ s, 4 ≤ s → Mt ^ (s + 1) - Mt ^ s = Mt ^ 5 - Mt ^ 4 := by
      intro s hs
      induction s, hs using Nat.le_induction with
      | base => rfl
      | succ s hs ih2 =>
        have := Mt_second_diff s hs
        have e : Mt ^ (s + 1 + 1) - Mt ^ (s + 1) = Mt ^ (s + 1) - Mt ^ s := by
          have : Mt ^ (s + 2) - 2 • Mt ^ (s + 1) + Mt ^ s = 0 := this
          rw [two_smul] at this
          have e2 : s + 1 + 1 = s + 2 := rfl
          rw [e2]
          have := sub_eq_zero.mp (show (Mt ^ (s + 2) - Mt ^ (s + 1)) - (Mt ^ (s + 1) - Mt ^ s) = 0 by
            rw [← this]; abel)
          exact this
        rw [e, ih2]
    have h1 := hdiff t ht
    have e3 : Mt ^ (t + 1) = Mt ^ t + (Mt ^ 5 - Mt ^ 4) := by
      rw [← h1]; abel
    rw [e3, ih]
    have : (t + 1 - 4) = (t - 4) + 1 := by omega
    rw [this, add_smul, one_smul]
    abel

/-! ### The counting sequence `N_m = u^T M^(m-c0) v` -/

/-- the functional `A ↦ u ⬝ A v`. -/
def phi (u v : Fin 6 → ℤ) (A : Matrix (Fin 6) (Fin 6) ℤ) : ℤ := u ⬝ᵥ (A *ᵥ v)

theorem phi_add (u v : Fin 6 → ℤ) (A B : Matrix (Fin 6) (Fin 6) ℤ) :
    phi u v (A + B) = phi u v A + phi u v B := by
  simp [phi, Matrix.add_mulVec, dotProduct_add]

theorem phi_sub (u v : Fin 6 → ℤ) (A B : Matrix (Fin 6) (Fin 6) ℤ) :
    phi u v (A - B) = phi u v A - phi u v B := by
  simp [phi, Matrix.sub_mulVec, dotProduct_sub]

theorem phi_nsmul (u v : Fin 6 → ℤ) (n : ℕ) (A : Matrix (Fin 6) (Fin 6) ℤ) :
    phi u v (n • A) = n * phi u v A := by
  induction n with
  | zero => simp [phi]
  | succ n ih =>
    rw [succ_nsmul, phi_add, ih]; push_cast; ring

/-- the transfer-matrix count `N_m`. -/
def Ncount (u v : Fin 6 → ℤ) (c0 m : ℕ) : ℤ := phi u v (Mt ^ (m - c0))

/-- `N_m` is affine in `m` for `m ≥ c0 + 4`. -/
theorem Ncount_affine (u v : Fin 6 → ℤ) (c0 m : ℕ) (hm : c0 + 4 ≤ m) :
    Ncount u v c0 m = Ncount u v c0 (c0 + 4) + (m - (c0 + 4) : ℕ) *
      (Ncount u v c0 (c0 + 5) - Ncount u v c0 (c0 + 4)) := by
  unfold Ncount
  rw [Mt_pow_affine (m - c0) (by omega), phi_add, phi_nsmul, phi_sub]
  have h4 : c0 + 4 - c0 = 4 := by omega
  have h5 : c0 + 5 - c0 = 5 := by omega
  have h6 : m - c0 - 4 = m - (c0 + 4) := by omega
  rw [h4, h5, h6]

/-- **Affine ⇒ identity.**  If `N_{c0+4} = c0+4` and `N_{c0+5} = c0+5`, then `N_m = m` for all
`m ≥ c0 + 4`.  (In the paper `c0 + 4 ≥ m1` and the two values are among `m = 7..40`.) -/
theorem Ncount_eq_id (u v : Fin 6 → ℤ) (c0 : ℕ)
    (h4 : Ncount u v c0 (c0 + 4) = (c0 + 4 : ℕ)) (h5 : Ncount u v c0 (c0 + 5) = (c0 + 5 : ℕ))
    (m : ℕ) (hm : c0 + 4 ≤ m) : Ncount u v c0 m = m := by
  rw [Ncount_affine u v c0 m hm, h4, h5]
  push_cast
  have : ((m - (c0 + 4) : ℕ) : ℤ) = (m : ℤ) - (c0 + 4) := by
    rw [Nat.cast_sub hm]; push_cast; ring
  rw [this]; ring

end P3.Hard
