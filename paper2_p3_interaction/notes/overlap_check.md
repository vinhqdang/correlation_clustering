# Overlap check: two flagged items and a fresh literature search

Date of check: 2026-10-20. Method: web fetch of arXiv abs and pdf (PDF text extracted locally with pdftotext), the DROPS volume page, the Semantic Scholar citations API, the arXiv query API, and extended web searches. No git, paper/ not touched.

Tags: V = read on the primary page or in the paper text. S = search snippet or aggregator only. N = negative result of a search (absence of evidence, not proof).

## Verdict summary

| Item | Verdict |
|---|---|
| (1) arXiv 2602.15683 | Merely relates. Different problem (Fair Correlation Clustering), different parameter (treewidth of the input graph). No contradiction with our results. The Cluster Editing corollary is about tw(G), not tw(I(H)) and not above a P3-packing. |
| (2) WG 2026 LIPIcs vol. 376 paper 17 | Merely relates. Different problem variant (weighted Correlation Clustering with weight-0 "fuzzy" pairs), parameter k plus degeneracy or closure of the fuzzy graph, results are polynomial kernels and hardness. No packing, no treewidth, no overlap with our hardness or XP results. |
| Fresh search | No paper found that studies Cluster Editing above a P3-packing other than Li-Pilipczuk-Sorge (LPS), none that settles c=2 beyond LPS's n^{2l+O(1)}, none with "CEaMP", none with treewidth of the conflict graph or intersection graph I(H). Negative results only (N). |

## 1. arXiv 2602.15683

- Title (V): "Fair Correlation Clustering Meets Graph Parameters".
- Authors (V): Johannes Blaha, Robert Ganian, Katharina Gillig, Jonathan S. Hojlev, Simon Wietheger (TU Wien; Univ. Bonn; DTU). Written "Højlev" in the paper.
- arXiv:2602.15683v1 [cs.DS], 17 Feb 2026. DOI 10.48550/arXiv.2602.15683. No venue stated in the PDF (V). No later version checked.
- Abstract (V, near verbatim): they study the generalization of Correlation Clustering with fairness constraints via fairlets; Fair Correlation Clustering "has so far lacked a detailed analysis from the parameterized complexity paradigm"; they give tractability results under structural graph parameterizations including treewidth, treedepth and vertex cover number.
- Problem (V): input graph G, budget B, kappa-colouring chi; every cluster must be partitionable into fairlets (a fairlet is a minimum multiset of colours with the same colour ratio as the whole graph; "fairlet size" c~ is its number of vertices). Correlation Clustering is defined there as modifying at most k edges to get disjoint cliques, i.e. unweighted Cluster Editing.
- Main theorems (V):
  - Thm 1: FCC is FPT in the vertex cover number of G.
  - Thm 2: FCC is in XP in tw(G) plus fairlet size.
  - Thm 3: FCC is FPT in tw(G) when the fairlet size is at most 2.
  - Thm 4: FCC is FPT in treedepth of G plus fairlet size.
  - Cor 1: "Fairness-oblivious Correlation Clustering (that is, Cluster Editing) is in FPT w.r.t. the treewidth of G." Fairlet size 1 is the unfair case, reduced to size 2 with tw increased by at most a factor 2.
- Open-question claim (V): the intro says FPT of Cluster Editing in tw "was not known and previously posed as an open question [10]" and Section 6 calls it "to the best of our knowledge" open. Reference [10] is not identified in my extraction (not checked further).
- Relation to LPS and P3-packings (V for absence in bibliography index): no Li/Pilipczuk/Sorge entry and no mention of P3-packings or modification-disjoint packings in the extracted text (grep on "P3", "packing", "Sorge" found nothing).
- Relation to our work: our parameters are l (excess above a modification-disjoint induced P3-packing H) and tw(I(H)), the intersection graph of the packing, not tw(G). tw(G) is a different, global parameter; Cluster Editing is FPT in k anyway, and this corollary is a statement about the unparameterized-solution-size version. Their Cor 1 does not imply or contradict any of: W[1]-hardness in l + tw(I(H)), the 5^{w+1} n^{2l+O(1)} algorithm, or the soft-triangle polynomial case. Their result is also not a statement about tw(I(H)) of a packing of a solution-size-bounded instance.
- Caution: the claim that tw-FPT for Cluster Editing was open is theirs; we have not independently verified the history of that question and our text should attribute it to them if used at all. Also note that I could not verify the correctness of their proof.
- Verdict: merely relates. Cite only if we want a sentence on tw(G) as a different structural parameter.

## 2. WG 2026, LIPIcs vol. 376, paper 17

- Title (V): "Clustering with Locally Bounded Ignorance".
- Authors (V): Jaroslav Garvardt and Christian Komusiewicz (University of Jena).
- Venue (V, DROPS volume page): LIPIcs Vol. 376, "52nd International Workshop on Graph-Theoretic Concepts in Computer Science (WG 2026)", paper 17, DOI 10.4230/LIPIcs.WG.2026.17. Page/article numbers not recorded; the DROPS fetch covered only the first 100000 of about 155000 characters of the page, but paper 17 was found there.
- arXiv version (V): arXiv:2605.13917v1 [cs.DS, cs.CC], 13 May 2026. The arXiv page lists no comment or journal reference.
- Abstract (V, PDF text): Correlation Clustering with weight function omega on pairs, omega >= 0 integer, cost = weight of edges between clusters plus missing edges inside clusters. Deciding cost at most k is FPT in k (via equivalence with Edge Multicut) but a polynomial kernel is presumably impossible. They study the "fuzzy edge graph" (weight-0 pairs): polynomial kernel for k + d (d = degeneracy of fuzzy graph) and for k + c (c = closure of fuzzy graph); hardness when the graph induced by edges and nonedges is very restricted.
- Main theorems (V):
  - Thm 1.1: kernel with 30 k^3 d vertices for k + d.
  - Thm 1.2: kernel with 30 k^2 c^2 vertices for k + c.
  - Thm 1.3: NP-hard and no polynomial kernel in k unless NP in coNP/poly, even if the nonedges form a matching and the real-edge graph has degeneracy 2 and closure 2 (unit weights).
  - Thm 1.4: NP-hard even if the real-edge graph is a binary tree or a star.
- Relation to LPS and packings (V): the reference list has no Li/Sorge entry and the text contains no "packing" or "above lower bound" discussion (grep on the extracted text). It cites van Bevern-Froese-Komusiewicz, ToCS 62(3):739-770, 2018 in the bibliography; I did not check where it is discussed in the text.
- Relation to our work: the parameter is the structure of the weight-0 pairs, our parameter is the excess over a P3-packing. Cluster Editing is the unweighted special case with no fuzzy pairs (fuzzy graph empty, d = c = 0) and is already FPT with a polynomial kernel in k, so their positive results say nothing about parameters smaller than k. Their hardness applies to weighted/fuzzy inputs, ours to plain Cluster Editing.
- Verdict: merely relates (kernelization of the weighted/fuzzy generalization). No overlap. Optional one-sentence citation as recent structural-parameter work on Correlation Clustering.
- venues.md says the WG paper is about "kernels by structural parameters of a fuzzy edge graph": correct (V).

## 3. Fresh search results

Queries run (web, extended mode): "cluster editing above packing P3 parameterized 2026"; "Cluster Editing Parameterized above Modification-disjoint P3-packings cited by 2025 2026"; "cluster editing parameterized treewidth of conflict graph OR modification-disjoint P3 packing 2025"; "CEaMP cluster editing above modification-disjoint packing c=2 ..."; the fuzzy-graph query for item 2. Also the arXiv query API sorted by date for "cluster editing" (40 most recent entries) and the Semantic Scholar citations endpoint for 1910.08517.

- LPS itself (V, arXiv abstract page): arXiv 1910.08517, Shaohua Li, Marcin Pilipczuk, Manuel Sorge, submitted 2019-10-18, online v3 2023-10-26, 41 pages 13 figures (arXiv comments). Abstract: NP-hard even for l = 0 and each vertex in at most 23 packed P3s; answers a question of van Bevern, Froese, Komusiewicz; asks for the largest c with tractability; c = 1 is FPT in l (van Bevern et al.); c = 2 solvable in |V|^{2l+O(1)} time. Journal: ACM TALG 20(1) Art. 3, DOI 10.1145/3626526 (from venues.md and search snippets, S). Note page counts differ between sources (41 pages on arXiv, 43 pages in venues.md for TALG); do not quote one without checking.
- Semantic Scholar citations of 1910.08517 (V, API output, 7 entries, offset 7 empty): Garvardt et al. "Graph Clustering Problems Under the Lens of Parameterized Local Search" (IPEC 2023 / JCSS 2026); Cluster Vertex Splitting and Company (arXiv 2309.00504); Bocci et al., "A New Temporal Interpretation of Cluster Editing" (arXiv 2202.01103, JCSS 2024); Crespelle et al. survey of edge modification (Comput. Sci. Rev. 2023); PACE 2021 Cluster Editing report (IPEC 2021); a thesis-like record; one duplicate record of LPS. None is by title about c = 2 or packings beyond LPS. I did not open them, so I cannot rule out a sentence about c = 2 inside any of them. Semantic Scholar coverage of 2026 papers may lag (N).
- arXiv recent listing (V, API): the newest "cluster editing" items are "Cluster Deletion is as Hard to Approximate as Vertex Cover" (arXiv 2608.04883, 2026-08-05), and "Fantastic Flips and Where to Find Them: A General Framework for Parameterized Local Search on Partitioning Problems" (2025-06-30). Neither is about packings or treewidth of a conflict graph (judged by title only).
- Other items surfaced by search, judged by title and snippet only (S): "On the Parameterized Complexity of s-Club Cluster Edge Deletion" (arXiv 2510.07065; W[1]-hard for pathwidth, a different problem); "Generalized Graph Packing Problems Parameterized by Treewidth" (arXiv 2509.06091, ESA 2025; H-Packing on bounded-treewidth graphs, not Cluster Editing); "Parameterized Algorithms for Editing to Uniform Cluster Graph" (arXiv 2404.10023); "Cluster Editing on Cographs and Related Classes" (STACS 2025).
- "CEaMP", "cluster editing above packing" with 2025-2026 date, "treewidth of conflict graph": no hits (N). I found no paper that settles c = 2 (FPT, polynomial or hard), and none that studies treewidth of I(H) (N).
- Not done: Google Scholar and the Semantic Scholar web UI were not reachable through these tools, only the API. A citing-works check via Google Scholar "cited by" before submission is still recommended.

## 4. Suggested related-work sentences (only verified facts)

1. Li, Pilipczuk and Sorge [LPS] showed that Cluster Editing parameterized by the excess l of the budget over the size of a modification-disjoint induced P3-packing H is NP-hard already for l = 0 with every vertex in at most 23 packed P3s, answering a question of van Bevern, Froese and Komusiewicz, and gave an algorithm running in time |V|^{2l+O(1)} when every vertex lies in at most two packed P3s (c = 2).
2. Van Bevern, Froese and Komusiewicz showed that the case c = 1 is fixed-parameter tractable in l [LPS abstract].
3. Blaha et al. [arXiv:2602.15683] study Fair Correlation Clustering and, as a corollary, state that unweighted Cluster Editing is FPT in the treewidth of the input graph, which they describe as an open question to their knowledge; this parameter is unrelated to the treewidth of the intersection graph of a packing used here.
4. Garvardt and Komusiewicz [WG 2026, LIPIcs 376, 17] give polynomial kernels for weighted Correlation Clustering parameterized by k plus the degeneracy or the closure of the graph of weight-0 pairs, and hardness for restricted edge and nonedge structure; their setting is kernelization in k, not parameterization above a packing.

Suggested bib entries (check page numbers and the exact DROPS metadata before use):
- Blaha, Ganian, Gillig, Hojlev, Wietheger. Fair Correlation Clustering Meets Graph Parameters. arXiv:2602.15683, 2026.
- Garvardt, Komusiewicz. Clustering with Locally Bounded Ignorance. WG 2026, LIPIcs vol. 376, art. 17, doi:10.4230/LIPIcs.WG.2026.17. Preprint arXiv:2605.13917.
- Li, Pilipczuk, Sorge. Cluster Editing Parameterized above Modification-disjoint P3-packings. STACS 2021 (LIPIcs 187, 49:1-49:16); ACM Trans. Algorithms 20(1), Art. 3, doi:10.1145/3626526; arXiv:1910.08517.

## 5. Caveats and suggested follow-ups

- Do not state that we settle LPS's "largest tractable c" question: our paper says c = 2 is open and our hardness uses c = O(l).
- Our claim that tw(I(H)) for c = 2 and general l is XP only matches LPS's |V|^{2l+O(1)} for c = 2; relation to LPS Theorem on c = 2 should be stated carefully (we did not re-read LPS's proof here).
- The notes in venues.md section 0 item 5 can be updated: the WG 2026 paper does not overlap. A one-line addition suggesting a "no overlap" conclusion would be appropriate but I have not edited venues.md.
- Before submission: repeat the check on Google Scholar "cited by" for 1910.08517 and 2602.15683, and re-check the arXiv listing for new Cluster Editing papers (the check above is dated 2026-10-20).
