# P3Proofs: Lean 4 + Mathlib formalisation of the core of paper2_p3_interaction

Files: `P3Proofs/Defs.lean` (definitions, cluster graphs, cost/OPT equivalence), `P3Proofs/CSP.lean`
(lower bound, `l = 0` CSP, scopes, five states). No `sorry`, no extra axioms (only `propext`,
`Classical.choice`, `Quot.sound`; checked with `#print axioms`).

## Modelling choices
* Graph: `SimpleGraph V` on a finite type `V` (`[Finite V]`). An edit set is a `SimpleGraph E`
  (its edges are the flipped pairs); `edit G E` is `G sd E`; `size E = E.edgeSet.ncard`.
* `IndP3 G a b c`: induced P3 with centre `b` (`ab, bc` edges, `ac` non-edge, `a != c`).
  `IsCluster G := forall a b c, not IndP3 G a b c`.
* Clustering: labelling `c : V -> alpha`; `cost G c = size (edit G (clusterGraph c))`, i.e. the number of
  pairs `ab` with `[ab in E] != [a,b in one cluster]`. `OPT` = min size of a solution; `OPTc` = min cost over
  clusterings with labels in `V`.
* Packing: `P : iota -> V x V x V` (`[Finite iota]`, `|H| = Nat.card iota`), each `P i` an induced P3,
  `pairsOf (P i)` pairwise disjoint (`IsPacking`). Injectivity of `P` is derived (`IsPacking.injective`).
  `IsPairMaximal` as in the paper. `IsSoft` = lies in some `pairsOf (P i)`.
* CSP variables `x : iota -> Fin 3` (0 = ab, 1 = bc, 2 = ac); `editOf P x` flips exactly the chosen pairs.

## Status table

| Lean name | Paper statement | Status |
|---|---|---|
| `isCluster_iff_trans` | cluster graph iff no induced P3 iff transitive (prelim) | proved |
| `isCluster_iff_exists_clustering` | solutions <-> partitions into cliques | proved (labels in `V`) |
| `tripleBad_iff` | "no vertex triple with exactly two edges" = no induced P3 (any centre) | proved (distinct triples) |
| `cost_eq_size`, `isSolution_edit_clustering`, `exists_clustering_of_solution` | cost of partition = size of edit set; solutions <-> partitions | proved |
| `OPT_eq_OPTc` | OPT over edit sets = OPT over partitions | proved |
| `needs_edit`, `choose_edits` | each member needs an edit (first half of Lemma 2.3, "Solutions decompose": every `x_i` nonempty) | proved |
| `card_le_size`, `card_le_OPT` | OPT >= abs(H) for every modification-disjoint packing | proved |
| `obs_one` | Obs. 3.1, first part: a solution with abs(S) <= abs(H) has exactly one pair in each `pairs(P_i)` and no other pair | proved (as: `E` meets each `pairsOf (P i)` in a singleton, all of `E` soft, `E = editOf P x`) |
| `csp_iff`, `opt_eq_card_iff` | Obs. 3.1: `(G,H,0)` YES iff the CSP (`x_i` in {1,2,3}, one "no induced P3 after flipping" constraint per triple) is satisfiable; also `OPT = abs(H)` form | proved |
| `edit_adj_agree`, `indP3_agree` | Obs. 3.1 "truth depends on owners' variables of the soft pairs of `t` only (hard pairs keep their status)" | proved (constraint of a triple depends only on `x` restricted to `owners P t`) |
| `owners_clique` | Obs. 3.2: owner set of a triple is a clique of `I(H)` | proved, as: two distinct owners share a vertex. NOT formalised: "at most three nodes", "primal graph subset of I(H)" as a graph statement, tree decompositions/treewidth |
| `five_states_fin`, `two_bad_states`, `five_states_cost` | Lemma 3.5 (five states) | proved, finite case check by `decide` |
| `five_states`, `member_state_five` | same, on graphs: for an induced P3 `(a,b,c)`, `G sd E` has no P3 on the triple iff the flipped-pair set is one of the five; hence in every solution each `x_i` is one of the five | proved |
| `IsPairMaximal` | pair-maximal packing | definition only |
| Lemma 2.3 (decomposition), second half | `abs(S) = sum abs(x_i) + abs(U)`, `abs(U) <= abs(S) - abs(H)` | not formalised (only `abs(S) >= abs(H)`) |
| Lemma 2.1 (degrees), Lemma 2.4 (hard families), Lemma 2.5 (maximalisation) | `abs(Q(v)) = 2 * #members through v`; `OPT >= abs(H)+abs(Phi)`; maximalisation (a)-(e) | not formalised |
| Lemma 3.4 (DP over `I(H)`), Thm 3.6 (`3^(w+1) n^O(1)`), Thm 3.7 (`5^(w+1) n^(2l+O(1))`) | algorithms / running times / treewidth | not formalised (no algorithmic or treewidth content) |
| Remark 3.8 (domain {1,2,3} insufficient for `l>0`) | two concrete instances | not formalised |
| Thm 4.1 (ETH), Prop 4.2 (NP-hard `w=0`), Prop 4.3 (stars/grids) | complexity-theoretic / treewidth | not formalised |
| Prop 4.4 (vertex-maximality), Prop 4.5 (C4) | small explicit instances (OPT, optima by enumeration) | not formalised |

## Weakenings / caveats
* Obs. 3.1 is stated for `OPT <= abs(H)` and `OPT = abs(H)`; the paper's "S a solution with `abs(S) <= abs(H)`" is `obs_one`.
* The five-state lemma: `five_states_fin` is the finite combinatorial fragment; `five_states` lifts it to graphs.
  Cost (size) of each state is recorded in `five_states_cost`.
* "Constraints only involve triples whose pairs lie in members sharing a vertex" is formalised as
  `indP3_agree` (dependence only on owners) plus `owners_clique` (owners pairwise share a vertex), not as a statement about
  treewidth.
