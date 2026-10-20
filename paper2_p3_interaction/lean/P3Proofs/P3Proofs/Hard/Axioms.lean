import P3Proofs.Hard.All

/-! Axiom audit of the main theorems of `P3Proofs.Hard.*`.  Every line must print only
`propext`, `Classical.choice`, `Quot.sound` (or fewer); in particular no `sorryAx` and no
`Lean.ofReduceBool` (the latter would indicate `native_decide`). -/

open P3.Hard

#print axioms le_iff_forall_clause
#print axioms posClauses_iff
#print axioms comparatorConditions_iff
#print axioms solvable_iff_comparators
#print axioms solvable_iff_swap_alt
#print axioms slotGroup_inj
#print axioms group_cover
#print axioms CountSystem.ex_ge_max
#print axioms CountSystem.k_sq_le_ex
#print axioms CountSystem.tight
#print axioms tightSystem_ex
#print axioms CountSystem.ex_ge_Z
#print axioms CountSystem.Z_ge
#print axioms printed_inequality_fails
#print axioms CountSystem.ex_ge_Z_E
#print axioms OH_valid
#print axioms chain_7
#print axioms chain_8
#print axioms chain_9
#print axioms dfs_complete
#print axioms Mt_charpoly
#print axioms Mt_ch
#print axioms Mt_pow_affine
#print axioms Ncount_affine
#print axioms Ncount_eq_id
#print axioms pOf_least
#print axioms cOf_le
#print axioms no_split_x1
#print axioms no_split_all
