Dear Editors,

I submit the manuscript "Cluster Editing above a P3 packing: W[1]-hardness in the excess and the treewidth of the intersection graph" for consideration in Algorithmica.

**Summary.** Li, Pilipczuk and Sorge parameterize Cluster Editing by the excess l of the solution size over a modification-disjoint packing H of induced P3's, and ask for which packings the problem is fixed-parameter tractable. The manuscript uses the treewidth w of the intersection graph of the packing as the structural measure. The main result is negative: the problem is W[1]-hard parameterized by l + w, by l + pathwidth and by l + c (c is the maximum number of packed P3's through a vertex), even for pair-maximal packings. Under ETH there is no algorithm in time f(l,w) n^{o(sqrt l)}. The reduction is from Grid Tiling with <=. On the positive side, a constraint-satisfaction formulation gives algorithms in time 3^{w+1} n^{O(1)} for l = 0 and 5^{w+1} n^{2l+O(1)} in general, and the manuscript identifies a configuration (the soft triangle) whose absence makes l = 0 polynomial for every c and w. Open questions are listed.

**Relation to earlier work.** The problem and the parameter l are due to Li, Pilipczuk and Sorge; their question is the starting point. For vertex-disjoint packings, fixed-parameter tractability in l follows from van Bevern, Froese and Komusiewicz (BFK); I use their reduction rule on components of the intersection graph and cite its soundness proof. The Cluster Editing results of Blaha et al. (arXiv:2602.15683) concern the treewidth of the input graph in the fair setting and a different parameter; the treewidth of the intersection graph used here is unrelated.

**Verification.** The reduction is computer-assisted at four finite steps, which the paper lists, and its counting core and several structural lemmas are machine-checked in Lean 4; the Lean files do not cover the graph-level statements of the reduction. Code and Lean files are at https://github.com/vinhqdang/correlation_clustering (paper2_p3_interaction/code and lean/P3Proofs).

**Length.** The full version is long (about 100 pages including appendices). The main text states every result and describes the hardness reduction and the algorithms; the long proofs of the reduction lemmas, the soft-triangle analysis and the further structural results (components, optimal solutions, locality, repacking) are in the appendices after the references. I am glad to shorten or reorganize if the editors prefer.

The manuscript is not under consideration elsewhere. I have no competing interests.

Sincerely,
Quang-Vinh Dang
British University Vietnam, Hung Yen, Vietnam
vinh.dq4@buv.edu.vn, ORCID 0000-0002-3877-8024

*Suggested reviewers:* omitted; I could not verify names and affiliations from official pages.
