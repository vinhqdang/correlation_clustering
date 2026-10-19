# Literature (verified 2026-10-19 against arXiv 1910.08517v3 and its abstract)

Li, Pilipczuk, Sorge. *Cluster Editing parameterized above modification-disjoint
P3-packings.* STACS 2021 (LIPIcs 187, 49); ACM Trans. Algorithms 20(1), Art. 3,
2023 (doi 10.1145/3626526). arXiv:1910.08517.

Problem CEaMP: input G, a packing H of induced P3s that share no edge and no
non-edge, and l >= 0; question: is there a cluster-editing set S with
|S| - |H| <= l?

Results as stated in the paper:
* Thm 1. NP-hard even for l = 0 and when every vertex lies in at most 23 packed
  P3s. (Answers a question of van Bevern, Froese, Komusiewicz, CSR 2016 / ToCS 2018.)
* Define a packing to be 1/c-integral if each vertex is in at most c packed P3s.
  - c = 1 (vertex-disjoint packing): FPT in l (van Bevern et al.).
  - c = 2 (half-integral): XP, time n^(2l + O(1)) (their Thm 2). Not known to be FPT.
  - 3 <= c <= 22: open in that paper.
* They remark that the reductions known for constant maximum degree give
  instances without tight P3-packing bounds (LP argument), so their hardness
  needed new gadgets.

Searches (WebSearch, 2026-10-19) found no later paper that extends the c = 2
result or introduces a treewidth-type parameter for this problem. A search is
not proof of absence: check citations of the paper on Semantic Scholar /
Google Scholar before claiming novelty.

Related (cited in paper 1, not yet re-read for this question):
* van Bevern, Froese, Komusiewicz. Parameterizing edge modification problems
  above lower bounds. ToCS 62(3), 2018 (arXiv 1512.04047).
* Komusiewicz, Uhlmann. Cluster editing with locally bounded modifications. DAM 2012
  (NP-hardness for small maximum degree; relevant to what bounded-degree
  structure cannot give).
* Blaesius et al., SEA 2022 (packings as lower bounds inside branch and bound):
  the practical motivation, since their packings are the instances of CEaMP.
