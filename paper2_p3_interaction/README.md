# Paper 2: when do local repair choices combine into a global optimum?

Working folder for the second paper. The first paper (certified lower bounds)
is in `../paper_mpc/`. This one is a theory paper; experiments are secondary.

## Question

Given a graph G, a packing H of modification-disjoint induced P3s and a
budget l, decide whether G can be turned into a cluster graph with at most
|H| + l edits (problem CEaMP of Li, Pilipczuk and Sorge). Find a measurable
interaction parameter w with an algorithm in f(l, w) * n^O(1) time, plus a
matching hardness or counterexample for weaker conditions.

## Status (see notes/)

* `notes/literature.md`: what is known, verified against the arXiv full text.
* `notes/analysis_01.md`: first structural analysis. It corrects one premise
  of the proposal: at l = 0 there is no transitivity-mediated interaction
  between disjoint P3s; the global difficulty sits in the excess edits (l > 0).
* `code/check_csp_view.py`: brute-force check of the l = 0 reduction.

Nothing here is a claimed result yet. Every statement in `notes/` is marked
as proved (with the proof), checked by computer on small cases, or open.
