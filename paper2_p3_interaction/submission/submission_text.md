# Submission details (plain text, copy and paste)

## Journal
Algorithmica (Springer)

## Article type
Research article / Original research

## Title
Cluster Editing above a P3 packing: W[1]-hardness in the excess and the treewidth of the intersection graph

## Abstract (about 243 words)
Cluster Editing asks for the fewest edge insertions and deletions that turn a graph into a disjoint union of cliques. A known parameterization measures the solution size above the size of a modification-disjoint packing H of induced paths on three vertices, called the excess ℓ, and asks for which packings the problem is fixed-parameter tractable. We use the treewidth w of the intersection graph of the packing as the structural measure. Our main result is negative: the problem is W[1]-hard parameterized by ℓ+w, also by ℓ plus the pathwidth of the intersection graph, and by ℓ+c, where c is the maximum number of packed paths through a vertex, even for pair-maximal packings. Hence there is no algorithm with running time f(ℓ,w)n^O(1) unless fixed-parameter tractability collapses to W[1], and under the Exponential Time Hypothesis none with running time f(ℓ,w)n^o(sqrt ℓ). The reduction, from Grid Tiling with ≤, is computer-assisted at four finite steps, which we list; its counting core and several structural lemmas are also machine-checked in the proof assistant Lean 4, which does not cover the graph-level statements of the reduction. Matching upper bounds come from a constraint-satisfaction view: time 3^(w+1) n^O(1) for ℓ=0 and 5^(w+1) n^(2ℓ+O(1)) in general; under the Exponential Time Hypothesis there is no 2^o(w) n^O(1) algorithm even for ℓ=0. We isolate a configuration, the soft triangle, as the source of difficulty: without soft triangles the case ℓ=0 is polynomial for every c and w. Further structural results and open questions are stated, with proofs in the appendices.

## Keywords (6)
Cluster Editing; Parameterized complexity; W[1]-hardness; Treewidth; Packing of induced paths; Constraint satisfaction

## MSC 2020 classification
68Q27 (Parameterized complexity), 68W40 (Analysis of algorithms), 05C85 (Graph algorithms)

## Author
Quang-Vinh Dang (corresponding author)
Affiliation: British University Vietnam, Hung Yen, Vietnam
Email: vinh.dq4@buv.edu.vn
ORCID: 0000-0002-3877-8024

## Declarations (as in the manuscript)
- Funding: No funding was received for conducting this study.
- Competing interests: The author has no relevant financial or non-financial interests to disclose.
- Ethics approval and consent to participate: Not applicable.
- Consent for publication: Not applicable.
- Materials availability: Not applicable.
- Data availability and code availability: all code, data and Lean files are at https://github.com/vinhqdang/correlation_clustering (directories paper2_p3_interaction/code and paper2_p3_interaction/lean/P3Proofs).
- Author contributions: single author.

## Files to upload
1. manuscript_latex.zip (editable LaTeX source: main.tex, sections/, refs.bib, main.bbl, sn-jnl.cls, sn-mathphys-num.bst, main.pdf; compiles to 113 pages)
2. manuscript.pdf (compiled PDF of the manuscript)
3. cover_letter.pdf

## Things to check on the portal before pressing submit
- The portal asks for Author contribution and Competing interest in its own form; only that information is used in the published version.
- Funding must also be entered in the portal (the wording above is a placeholder, valid only if there truly is no funding).
- Suggested reviewers (optional): the corresponding author must give an institutional email or a verifiable link for each; none are suggested in the cover letter.
- Policy on AI and LLM use (Springer Nature): use of a large language model beyond copy editing should be documented in the Methods section or a suitable alternative part of the manuscript. The manuscript currently contains no such statement. This is your decision to make before submitting; see the message in the chat.
- No simultaneous submission elsewhere (do not submit to a conference or another journal in parallel).
- There is no submission or publication fee. Open access is optional after acceptance.
