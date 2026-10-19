# Literature and novelty check (2026-10-19)

Convention: [OPENED] = I fetched and read the URL (or a PDF downloaded from it) in this session.
[SEARCH-ONLY] = seen only as a web-search snippet. UNVERIFIED = not opened.
Nothing below is from memory.

## 1. Papers citing Li-Pilipczuk-Sorge (arXiv 1910.08517)

Source of the list: Semantic Scholar API, citations endpoint, queried 2026-10-19 (first calls returned
HTTP 429, a retry succeeded) [OPENED]:
https://api.semanticscholar.org/graph/v1/paper/arXiv:1910.08517/citations?fields=title,year,venue,externalIds&limit=100
It returned 7 records (one is a duplicate record of Li et al. itself, CorpusId 270434514). Semantic Scholar
coverage is incomplete; OpenAlex gives cited_by_count = 2 for the arXiv record
(https://api.openalex.org/works/doi:10.48550/arXiv.1910.08517 [OPENED]); its cited-by list could not be
fetched (daily budget exhausted; I did not create an API key, which would email the user). Google Scholar
was not accessible.

| # | Citing paper | What I checked | Relevant to c=2 / 3<=c<=22 / treewidth? |
|---|---|---|---|
| 1 | Crespelle, Drange, Fomin, Golovach, "A survey of parameterized algorithms and the complexity of edge modification", arXiv 2001.06867 | [OPENED] https://arxiv.org/pdf/2001.06867 , text lines 505-520 | Only restates: Li et al. show NP-completeness for l=0. No improvement on c=2, no treewidth. |
| 2 | "The Complexity of Cluster Vertex Splitting and Company", arXiv 2309.00504 | [OPENED] https://arxiv.org/pdf/2309.00504 | Uses modification-disjoint P3 packings only as a lower-bound tool (Prop. 6.3). Nothing on c or treewidth. |
| 3 | Bocci et al., "A New Temporal Interpretation of Cluster Editing", arXiv 2202.01103 | [OPENED] https://arxiv.org/pdf/2202.01103 | One sentence in the introduction citing the parameter "above the packing lower bound". Nothing new. |
| 4 | Garvardt, Morawietz, Nichterlein, Weller, "Graph Clustering Problems Under the Lens of Parameterized Local Search", IPEC 2023 (JCSS 2026) | [OPENED] Dagstuhl PDF https://drops.dagstuhl.de/storage/00lipics/lipics-vol285-ipec2023/LIPIcs.IPEC.2023.20/LIPIcs.IPEC.2023.20.pdf ; grep shows citation in "above lower-bound parameterizations [6, 37]" and no treewidth result | Cites only. Local-search parameterization, unrelated. |
| 5 | Kellerhals, Koana, Nichterlein, Zschoche, "The PACE 2021 Challenge: Cluster Editing", IPEC 2021 | UNVERIFIED (not opened; Semantic Scholar record only) | Presumably practical; cannot say. |
| 6 | Thesis "Parameterized algorithms for connectivity, separation, and modification problems in graphs" (2021, Semantic Scholar CorpusId 259380955) | UNVERIFIED (not opened; no author or text seen) | Cannot say. |
| 7 | Duplicate record of Li et al. (CorpusId 270434514) | n/a | n/a |

Result of (1): among the citing papers I could open (1-4), none improves c=2, resolves 3<=c<=22, or uses
treewidth-like parameters. Items 5 and 6 are UNVERIFIED, and the citation lists are incomplete, so this is a
"not found" and not a proof of absence. Earlier WebSearch queries on the title also returned only the
paper itself (https://arxiv.org/abs/1910.08517 [SEARCH-ONLY]).

What Li et al. say about their own open end (from the arXiv PDF, https://arxiv.org/pdf/1910.08517 [OPENED]):
* Thm 1: NP-hard for l=0 with every vertex incident to at most 23 packed P3s.
* Thm 2: half-integral packing (c=2) is in XP, n^(2l+O(1)).
* Conclusions, verbatim in substance: they leave open whether c=2 is FPT; "we discovered large grid-like
  structures that seemed difficult to overcome in fixed-parameter time, and a corresponding W[1]-hardness
  result would also not be surprising." They suggest examining which substructures of the hardness
  reduction are rare in practice. No treewidth statement anywhere (grep for "treewidth" finds only
  reference-list noise). This is useful for us: the authors themselves suspect W[1]-hardness for c=2 in l,
  so an FPT result in (l, tw(I(H))) must be consistent with grids appearing in c=2 instances
  (I(H) is subcubic when c=2 and a grid-like I(H) has large treewidth, so there is no contradiction yet).

## 2. Cluster Editing / Deletion and treewidth: what is known

Honest statement: I did NOT find any paper that states the complexity of plain Cluster Editing (or Cluster
Deletion) parameterized by treewidth alone, or by tw + k, tw + number of clusters.

Searched: several WebSearch queries combining "cluster editing", "treewidth", "W[1]-hard", "FPT", "bounded
treewidth". Findings that are adjacent but not the answer:
* s-Club Cluster Edge Deletion: FPT in s + tw (Montecchiani et al.); W[1]-hard in pathwidth, hence
  treewidth, even with pathwidth + number of clubs (arXiv 2510.07065 [OPENED, text read]
  https://arxiv.org/pdf/2510.07065). It says the s=1 case is "extensively investigated" and lists twin cover,
  neighborhood diversity, cluster vertex deletion number as FPT parameters for Cluster Deletion, and that
  Cluster Deletion is NP-hard on planar graphs and on C4-free graphs of max degree four. It does not
  state a treewidth result for s=1. Caution: that is a deletion problem with a diameter bound; it does not
  transfer to Cluster Editing, where added edges can join far-apart vertices.
* Survey 2001.06867 [OPENED]: mentions treewidth for graph classes and for other problems. The sentence
  "Dvorak and Knop showed W[1]-hard in pathwidth, FPT in treedepth" (survey lines 1126-1127) is about
  length-bounded cuts, not Cluster Editing. I checked this; do not misattribute it.
* Cluster Editing on cographs / clique-width + number of clusters (arXiv 2412.12454 [OPENED, grep]):
  NP-complete on cographs; no treewidth statement found in the text.
* Komusiewicz-Uhlmann, DAM 2012 [OPENED] https://www.uni-marburg.de/de/fb12/arbeitsgruppen/algorith/paper/ce-bounded-degree.pdf :
  Cluster Editing NP-hard on max degree 6; Cluster Deletion NP-hard on 4-regular graphs, polynomial
  (O(n^1.5 log^2 n)) at max degree 3; no 2^o(k) algorithm under ETH even at bounded degree; kernels of
  4dt (editing) and 2dt (deletion) vertices for d clusters and t modifications per vertex. They state no
  treewidth result. Consequence for us: bounded degree (hence bounded-degree I(H)) does not give
  tractability of the full problem, consistent with Li et al. Thm 1.
* Fomin, Kratsch, Pilipczuk, Pilipczuk, Villanger, JCSS 2014 (p-Cluster Editing, 2^O(sqrt(pk)) and ETH
  tight): [SEARCH-ONLY] https://arxiv.org/pdf/1112.4419 (preprint "Subexponential fixed-parameter
  tractability of cluster editing"), https://folk.uib.no/nmiff/articles/2014/2014d.pdf . Not opened. Number of
  clusters, not treewidth.
* Gramm et al. and Boecker-Damaschke: only seen as survey citations (survey text line ~506: 2.27^k / 1.76^k
  Gramm, 1.76^k / 1.41^k Boecker-Damaschke). Not opened as primary sources.
* Crespelle et al. survey (full survey): no Cluster Editing + treewidth row found by grep.

My own reasoning, NOT a citation (to be proved or checked before use): Cluster Deletion is
MSO2-expressible, so FPT in treewidth by Courcelle is plausible; Cluster Editing allows insertions between
vertices at arbitrary distance in G, so tw(G) is not a natural parameter and may not give FPT. Treat as
a question to settle, not a fact.

## 3. van Bevern, Froese, Komusiewicz (arXiv 1512.04047, ToCS 62(3) 2018)

[OPENED] https://arxiv.org/abs/1512.04047 and PDF https://arxiv.org/pdf/1512.04047 (pdftotext, lines cited).
* Setting: F-free Editing with a given VERTEX-DISJOINT packing H of induced F's; parameter l = k - h(H)
  where h(H) is the lower bound from the packing (arXiv abstract page).
* Cluster Editing results (intro, lines ~201-203): O(1.62^((2t+1)l) + nm + n*1.62^t)-time algorithm and an
  O(t*l)-vertex kernel for cost-t packings (t in O(log n) needed for kernel time), and a 4^l * n^O(1)
  algorithm for P3 packings.
* Precisely, for P3 packings (Thm 6.13, line 1046): O(4^l * l^3 + m + n) time. The generic route via
  Reduction Rule 6.2 with t=1 gives O(4.26^l + n + m) (line 995). This is the c=1 FPT result. The proof is a
  branching algorithm with Reduction Rules 6.8, 6.9, 6.10, 6.12 and Branching Rules 1-4. Branching Rule 1:
  an induced P3 containing at most one vertex of each packed P3 is branched into its 3 repairs (parameter
  drops by one in each).
* Hardness (Sect. 7): K3-free Editing NP-hard at l=0 for edge-disjoint K3 packings; K_q-free editing
  (q >= 6) NP-hard at l=0 even for vertex-disjoint packings (abstract page).
* Open questions (Sect. 8, lines ~1360-1375): edit problems for 4-vertex forbidden graphs (Cograph, K4-free,
  Claw-free, Diamond-free); Triangle Vertex Deletion; "it is open to determine the complexity of Cluster
  Editing and Feedback Arc Set in Tournaments parameterized above the size of edge-disjoint packings of
  forbidden induced subgraphs"; better running times. The question about modification-disjoint P3
  packings is framed this way by Li et al. (their text: van Bevern et al. "asked whether Cluster Editing is
  fixed-parameter tractable when parameterized above a stronger lower bound, the size of a
  modification-disjoint packing of P3s"), answered NO for l=0 (NP-hard), with the c-integral refinement left
  partly open.
* No treewidth or intersection-graph parameter appears in this paper.

## 4. Closest prior work on "above guarantee / packing" and CSP-on-treewidth

* Mahajan, Raman, Sikdar, "Parameterizing above or below guaranteed values", JCSS 75(2):137-153, 2009
  [SEARCH-ONLY] https://tcs.rwth-aachen.de/~sikdar/index_files/above_guar_jour.pdf (not opened). Origin of the
  "tight lower bound" parameterization idea. Mahajan-Raman 1999 (MaxSat, MaxCut) likewise search-only.
* van Bevern et al. (Sect. 3 above) is the direct packing-bound predecessor; Li et al. is the successor.
* Survey 2001.06867 [OPENED] Open problem 2.3 (Cograph Editing above P4 packing) shows the above-packing
  line is considered open in the survey; no treewidth variant mentioned.
* CSP of bounded domain on treewidth: my Observation 1-2 (analysis_01.md) reduce l=0 to a domain-3 CSP with
  clique scopes in I(H). The generic facts used (a CSP with domain d and primal graph of treewidth w is
  solvable in d^(w+1) poly time) are standard but I did NOT open a source for them in this session:
  UNVERIFIED citation, please add a textbook reference (e.g. Cygan et al., Parameterized Algorithms,
  Springer 2015, which is cited in 2309.00504, ref 16; the specific theorem was not opened).
* I found no paper that defines the intersection graph of a P3 packing (or of a packing of forbidden
  subgraphs) as a parameter for Cluster Editing/Deletion or related edit problems.

## Verdict

As far as this search can establish, the question (l, tw(I(H))) for CEaMP is novel: none of the four
citing papers I could open, nor the van Bevern et al. paper, nor Komusiewicz-Uhlmann, introduces a treewidth
or intersection-graph parameter for packings, and no paper resolves c=2 or 3<=c<=22. Caveats: citation
coverage is incomplete (Semantic Scholar 7 records, OpenAlex count 2; two citing items, the PACE 2021 report
and a 2021 thesis, were not opened; Google Scholar unavailable), and I found no source at all on the
complexity of plain Cluster Editing parameterized by treewidth, so claims about tw(G) must be proved rather
than cited. Note the l=0 case with tw(I(H)) is easy by Observation 2, so the real novelty and difficulty
lie in l > 0, where Li et al. themselves suspect W[1]-hardness for c=2. The three closest papers are:
(1) Li, Pilipczuk, Sorge, arXiv 1910.08517 (the problem, Thm 1 and 2, the open c-range);
(2) van Bevern, Froese, Komusiewicz, arXiv 1512.04047 (c=1 FPT in 4^l, framework, origin of the question);
(3) Komusiewicz-Uhlmann, DAM 2012 (hardness at bounded degree, which limits what a degree or locality
condition on I(H) can buy). The survey arXiv 2001.06867 is the best map of the area but not a close
result.
